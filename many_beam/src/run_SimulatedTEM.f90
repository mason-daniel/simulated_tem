
    program run_SimulatedTEM
!---^^^^^^^^^^^^^^^^^^^^^^^^
!*      read in a file, compute the simulated TEM image
!*      and output as png
!*
!*      Daniel Mason
!*      (c) UKAEA October 2024
!*      
!*      version history
!*          v.0.0.1         Oct 24      First working version

       
        use iso_fortran_env
        use Lib_Lattices
        use Lib_Elements
        use Lib_XYZFiles
        use Lib_Filenames
        use Lib_Png
        use Lib_IntegrateManyBeams
        use Lib_Quaternions
        use Lib_Gvectors
        use Lib_Callipers
        use Lib_CommandLineArguments
#ifdef MPI
        use mpi_f08
#endif
        implicit none


        character(len=8),parameter          ::      VERSION = "0.0.1"


        real(kind=real64),parameter         ::      PI = 3.14159265390d0
        integer                             ::      rank = 0,nProcs = 1


    !---    command line parameter input        
        type(CommandLineArguments)          ::      cla
        
        character(len=256)                  ::      filename = ""           !   input filename
        character(len=256)                  ::      outfile = ""            !   output filename prefix
        character(len=8)                    ::      latticename = UNKNOWN_LATTICE
        real(kind=real64)                   ::      T = 300.0d0             !   temperature (K)
        real(kind=real64)                   ::      V = 200.0d0             !   accelerator voltage (kV)
        real(kind=real64),dimension(3)      ::      a0_in = 0.0             !   indicative lattice parameter (A)
        integer                             ::      nPrec = 1               !   precession angles
        real(kind=real64)                   ::      precAngle = 0.003d0         !   precession angle (rad)
        logical                             ::      opPng = .true.
        logical                             ::      columnar = .true.       !   use columnar approximation
        logical                             ::      lossy = .false.         !   use imaginary parts of crystal structure factor
        logical                             ::      twobeam = .false.       !   two-beam mode
        integer,dimension(:),allocatable    ::      hkl
        logical                             ::      yflip = .true.          !   Right-hand rule for png output - y=0 at bottom.
        logical                             ::      pngBlur = .true.        !   blurring radius sigma/a
        logical                             ::      pngNorm = .true.        !   normalise output intensity
        real(kind=real64)                   ::      n_g = 1                 !   which reflection to make bright

        real(kind=real64)                   ::      tilt_max = 5.0          !   maximum stage tilt (deg)
        logical                             ::      opDiffPatt = .true.     !   output diffraction pattern
        logical                             ::      df = .true.             !   find dark field
        

    !---    physical variables
        type(IntegrateManyBeams)            ::      imb
        integer                             ::      Nx,Ny           !   image size
        real(kind=real64),dimension(:,:),allocatable        ::      img
        !integer,dimension(:),allocatable    ::      hkl_bright      !   = n_g*hkl - which spot should be bright
        !real(kind=real64),dimension(3,3)    ::      foil_tilt,foil_tilt2
        integer,dimension(4)            ::      dat
        character(len=256)              ::      dummy
        integer             ::      ierror
        integer             ::      ii,ix,iy,iyp
        real(kind=real64)   ::      dd!,I_g

    !---    timing
        integer,parameter               ::      TIMER = kind(1_int32)
        integer(kind=TIMER),parameter   ::      T_INIT      = 1_TIMER
        integer(kind=TIMER),parameter   ::      T_CTOR      = 2_TIMER
        integer(kind=TIMER),parameter   ::      T_TILT      = 3_TIMER
        integer(kind=TIMER),parameter   ::      T_OUTPUT    = 4_TIMER
        integer(kind=TIMER),parameter   ::      T_PHASE     = 5_TIMER
        integer(kind=TIMER),parameter   ::      T_INTEGRATE = 6_TIMER
        integer(kind=TIMER),parameter   ::      T_TOTAL     = 7_TIMER
        type(Callipers),dimension(T_TOTAL)      ::      tt




        tt(T_TOTAL) = Callipers_ctor() 
        tt(T_INIT) = Callipers_ctor()
#ifdef MPI        
        call MPI_INIT(ierror)
        call MPI_COMM_SIZE(MPI_COMM_WORLD, nProcs, ierror)
        call MPI_COMM_RANK(MPI_COMM_WORLD, rank, ierror)
#endif              
        


!******************************************************************************    
!*    
!*      COMMAND LINE ARGUMENTS
!*    
!*    
!******************************************************************************    
    
        
        
    !---    read command line arguments
        cla = CommandLineArguments_ctor(30)  
         
        call setProgramDescription( cla, "run_SimulatedTEM" )
        call setProgramVersion( cla, VERSION )   
          
        
        call get( cla,"f",filename ,LIB_CLA_REQUIRED,"            input filename" )           
        
        ii=0; dat = 0; call get( cla,"g",dat,ii ,LIB_CLA_REQUIRED,"            g-vector reflection to output" )
        if ( (ii==3).or.(ii==4) ) then
            allocate(hkl(ii))
            hkl(1:ii) = dat(1:ii)       
        end if


        outfile = filename                                       
        ii=0 ; call get( cla,"a0",a0_in,ii ,LIB_CLA_OPTIONAL,"         lattice parameter(s)" )
        if (ii==1) a0_in(2:3) = a0_in(1)
        call get( cla,"o",outfile ,LIB_CLA_OPTIONAL,"          output filename prefix" )           
        if ( (trim(getSuffix(outfile))=="xyz").or.(trim(getSuffix(outfile))=="lammps").or.(trim(getSuffix(outfile))=="lmp").or.(trim(getSuffix(outfile))=="dat").or.(trim(getSuffix(outfile))=="png") ) then
            outfile = trim(removeSuffix(outfile))
        end if
        call get( cla,"lattice",latticename ,LIB_CLA_OPTIONAL,"    lattice type" )                                                        
        call get( cla,"T",T ,LIB_CLA_OPTIONAL,"          temperature (K)" )                                                        
        call get( cla,"V",V ,LIB_CLA_OPTIONAL,"          accelerating voltage (kV)" )                                                        
        call get( cla,"png",opPng ,LIB_CLA_OPTIONAL,"        output result as .png file" )                           
        call get( cla,"blur",pngBlur ,LIB_CLA_OPTIONAL,"       blur .png file by smoothing width sigma" )                                                        
        call get( cla,"norm",pngNorm ,LIB_CLA_OPTIONAL,"       normalise .png file intensity" )                                                        
        call get( cla,"nPrecAngle",nPrec ,LIB_CLA_OPTIONAL," number of precession images to take" )                                                        
        call get( cla,"precAngle",precAngle ,LIB_CLA_OPTIONAL,"  precession angle (rad)" )             
        call get( cla,"col",columnar ,LIB_CLA_OPTIONAL,"        use columnar approximation" )                 
        call get( cla,"loss",lossy ,LIB_CLA_OPTIONAL,"       use complex crystal structure factor for inelastic electron scattering" )    
        call get( cla,"2",twobeam ,LIB_CLA_OPTIONAL,"          compute in two beam approximation" )                 
        call get( cla,"yflip",yflip ,LIB_CLA_OPTIONAL,"      output png with y=0 at bottom, to match ovito orientation z out of screen" )                 
        call get( cla,"ng",n_g ,LIB_CLA_OPTIONAL,"         tilt foil until this reflection bright" )                 
        df = (n_g /= 1.0d0)
        call get( cla,"dark",df ,LIB_CLA_OPTIONAL,"       tweak foil tilt to find close dark field" )                 
        call get( cla,"theta",tilt_max,LIB_CLA_OPTIONAL,"      maximum foil tilt angle" )           
        
        
        if (rank==0) call report(cla)
        if (hasHelpArgument(cla)) call errorExit("done")
        if (.not. allRequiredArgumentsSet(cla)) call errorExit("error")
        call delete(cla)
        
        if (.not. allocated(hkl)) call errorExit("error - expected -g h,k,l or -g h,i,k,l")


!******************************************************************************    
!*    
!*      WELCOME
!*    
!*    
!******************************************************************************    
        
        if (rank==0) then
            print *,"run_SimulatedTEM"
            print *,"^^^^^^^^^^^^^^^^"
            call get_environment_variable("HOSTNAME",dummy)
            if (len_trim(dummy)==0) dummy="localhost"           !   ??
            print *,"   running on ",nProcs," processors on "//trim(dummy)
            print *,"   input atom positions    """//trim(filename)//""""
            print *,"   output file prefix      """//trim(outfile)//""""
            if (trim(latticename)==UNKNOWN_LATTICE) then
                print *,"   lattice type TBD"
            else
                print *,"   lattice type            """//trim(latticename)//""""
            end if
            print *,"   output g-vector         ",hkl
            print *,"   accelerator voltage     ",V," (kV)"
            print *,"   temperature             ",T," (K)"
            if (any(a0_in<=0)) then
                print *,"   indicative latt param TBD"
            else
                print *,"   indicative latt param   ",a0_in
            end if
            print *,"   use columnar approx     ",columnar
            print *,"   inelastic scattering    ",lossy
            print *,"   two beam mode           ",twobeam
            print *,"   output .png             ",oppng," with blur? ",pngBlur," normalised? ",pngNorm
            print *,"   bright reflection n_g = ",n_g
            print *,"   dark field tweak?       ",df
            print *,"   max foil tilt angle     ",tilt_max," (deg)"
            print *,"   number precession imgs  ",nPrec
            print *,"   precession angle        ",precAngle," (rad)"
            print *,""  
        end if


!******************************************************************************    
!*    
!*      CONSTRUCTOR
!*    
!*    
!******************************************************************************    

        ! if (rank==0) then
        !     print *,""
        !     print *,"initialising"
        !     print *,"^^^^^^^^^^^^"
        ! end if
        call Lib_IntegrateManyBeams_init_MPI()
        call pause(tt(T_INIT))


    !---    constructor - read input file and determine parameters of calculation
        tt(T_CTOR) = Callipers_ctor()
        imb = IntegrateManyBeams_ctor( latticename,a0_in,filename,T,V, nPrec = nPrec, precAngle = precAngle, columnar = columnar, lossy = lossy )
        call pause(tt(T_CTOR))

        
    !---    tilt the foil to "good" diffration conditions
        tt(T_TILT) = Callipers_ctor()
        call tiltFoil(imb,tilt_max*PI/180,hkl,n_g,df)

        if (opDiffPatt) then
            call outputDiffractionPattern( imb,trim(outfile)//"_before",1024 )
        end if

        if (twobeam) then
            call changeGvectors(imb,rule = LIB_IMB_BLOCK_HIGH_SG,hkl_in = hkl) 
        else
            call changeGvectors(imb,rule = LIB_IMB_BLOCK_HIGH_SG,rho_in = 2.0d0)
        end if
        call setImagingSpace(imb)                        
        call pause(tt(T_TILT))


    !---    compute the phase field x = exp[ - ig.u(r) ]        
        tt(T_PHASE) = Callipers_ctor()
        call computePhaseFields(imb)     
        call pause(tt(T_PHASE))



        if (rank==0) then
            print *,""
            print *,"running"
            print *,"^^^^^^^"
            print *,""
            call report(imb)
            print *,""
        end if


!******************************************************************************    
!*    
!*      INTEGRATE
!*    
!*    
!******************************************************************************    
        if (rank==0) then
            print *,""
            print *,"propagating beams"
            print *,"^^^^^^^^^^^^^^^^^"
        end if
        tt(T_INTEGRATE) = Callipers_ctor()
        call integrate( imb )
        call pause(tt(T_INTEGRATE))
        
!******************************************************************************    
!*    
!*      OUTPUT
!*    
!*    
!******************************************************************************    


        tt(T_OUTPUT) = Callipers_ctor()
        if (opDiffPatt) then
        !    call outputDiffractionPattern( imb,trim(outfile),1024 )
        end if

        if (opPng) then
            if (rank==0) then
                print *,""
                print *,"output .png"
                print *,"^^^^^^^^^^^"
            end if
    
            call getIntensity( imb,hkl,img )      !   gathers result on rank 0
 


            if (rank==0) then
                Nx = size(img,dim=1)
                Ny = size(img,dim=2)
                print *,"run_SimulatedTEM info - image size ",Nx,"x",Ny 
                write(*,fmt='(a,3f10.6)') " run_SimulatedTEM info - minmaxavg img ",minval(img),maxval(img),sum(img)/(Nx*Ny)


                if ((rank==0).and. pngBlur) call blurImg( img, getSigma( imb )/getA( imb ) )
 


                if (yflip) then
                    !   flip image in y-axis
                    do iy = 0,int( (Ny-1)/2 )      !   round down
                        iyp = Ny-1-iy
                        do ix = 0,Nx-1
                            dd = img(ix,iy)
                            img(ix,iy) = img(ix,iyp)
                            img(ix,iyp) = dd
                        end do
                    end do
                end if
                call writePng( trim(outfile)//".png",img,normalise=pngNorm )
            end if
        end if
        call pause(tt(T_OUTPUT))


!******************************************************************************    
!*    
!*      BYE BYE
!*    
!*    
!******************************************************************************    

        if (rank==0) then
            print *,""
            print *,"timing data (s)"
            print *,"^^^^^^^^^^^^^^^"
            write(*,fmt='(a,f16.6)') "   initialisation ",elapsed(tt(T_INIT))
            write(*,fmt='(a,f16.6)') "   constructor    ",elapsed(tt(T_CTOR))
            write(*,fmt='(a,f16.6)') "   tilt foil      ",elapsed(tt(T_TILT))
            write(*,fmt='(a,f16.6)') "   phase factor   ",elapsed(tt(T_PHASE))
            write(*,fmt='(a,f16.6)') "   integration    ",elapsed(tt(T_INTEGRATE))
            write(*,fmt='(a,f16.6)') "   output         ",elapsed(tt(T_OUTPUT))
            write(*,fmt='(a,f16.6)') "   total          ",elapsed(tt(T_TOTAL))
        end if

        if (rank==0) then
            print *,""
            print *,"bye bye"
            print *,"^^^^^^^"
            print *,""
        end if
        call errorExit("done")
        if (rank==0) print *,""

               
    contains
!---^^^^^^^^
!        
        
        
        subroutine blurImg(img,sigmaona)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      perform a simple 2d blur operation on an image.
    !*      This can be useful in low contrast images, where tiny fluctuations in 
    !*      atom density or grad x lead to (apparent) differences in the output image.
    !*      blur radius is sigma/a pixels.
            real(kind=real64),dimension(0:,0:),intent(inout)        ::      img
            real(kind=real64),intent(in)                            ::      sigmaona      !   blurring
            
            integer             ::      Nx,Ny
            integer             ::      ix,iy,ik,kk
            real(kind=real64),dimension(:),allocatable      ::      img_tmp
            real(kind=real64),dimension(:),allocatable      ::      kernel
            real(kind=real64)   ::      ff,ww
            Nx = size(img,dim=1)
            Ny = size(img,dim=2)
            allocate(img_tmp(0:max(Nx,Ny)-1))
            write(*,fmt='(a,f8.3,a)') " run_SimulatedTEM::blurImg info - blur ",sigmaona," (px)"
            kk = ceiling(3*sigmaona)
            allocate(kernel(-kk:kk))
            do ik = -kk,kk
                kernel(ik) = exp( -ik*ik*0.5d0 )
            end do

        !---    blur in x direction for each row iy
            do iy = 0,Ny-1
                img_tmp(0:Nx-1) = img(0:Nx-1,iy)
                do ix = 0,Nx-1
                    ff = 0 ; ww = 0
                    do ik = max(0,ix-kk),min(Nx-1,ix+kk)
                        ww = ww + kernel(ik-ix)
                        ff = ff + kernel(ik-ix)*img_tmp(ik)
                    end do
                    img(ix,iy) = ff/ww
                end do                
            end do

        !---    blur in y direction for each column ix
            do ix = 0,Nx-1
                img_tmp(0:Ny-1) = img(ix,0:Ny-1)
                do iy = 0,Ny-1
                    ff = 0 ; ww = 0
                    do ik = max(0,iy-kk),min(Ny-1,iy+kk)
                        ww = ww + kernel(ik-iy)
                        ff = ff + kernel(ik-iy)*img_tmp(ik)
                    end do
                    img(ix,iy) = ff/ww
                end do                
            end do
            return
        end subroutine blurImg

        subroutine tiltFoil(imb,tilt_max,hkl,n_g,dark)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ 
    !*      find a good tilt angle 0:tilt_max (radians) so that the reflection
    !*          ng [hkl] is bright. Optionally then tweak until the reflection [hkl] is dark

            type(IntegrateManyBeams),intent(inout)              ::      imb
            real(kind=real64),intent(in)                        ::      tilt_max
            integer,dimension(:),intent(in)                     ::      hkl
            real(kind=real64),intent(in)                        ::      n_g
            logical,intent(in)                                  ::      dark
            real(kind=real64),dimension(3,3)    ::      foil_tilt

            integer,parameter                   ::      NTILT = 100
            real(kind=real64),parameter         ::      TWEAK = 0.01d0

            logical                             ::      MillerBravais 
            integer,dimension(size(hkl))        ::      hkl_bright
            real(kind=real64)                   ::      I_g
            real(kind=real64),dimension(3,3)    ::      foil_tilt2
            character(len=16)                   ::      form
            logical                             ::      maxIntensity = .false.      !   select tilt based on foil exit intensity
             


            MillerBravais = (size(hkl)==4)
            form = "(a,3i4,a,f16.8)"
            if (MillerBravais) form = "(a,4i4,a,f16.8)"

            if (rank==0) then
                call perfectLatticeIntensity(imb,hkl,maxIntensity,I_g)
                write(*,fmt=form) "tiltFoil info - perfect lattice intensity ",hkl," before tilt ",I_g
            end if

            if (tilt_max<=0) return
 
    
            if (.not. ( (n_g==1.0).and.df) ) then
                if (abs(n_g - nint(n_g)) < 1.0d-6) then
                    if (rank==0) then
                        print *,""
                        write(*,fmt='(a,i2)') " selecting bright reflection n_g = ",nint(n_g)
                        print *,"^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^"
                        print *,""
                    end if
                    hkl_bright = nint( hkl * n_g )
                    call selectFoilTilt( imb,hkl_bright,tilt_max,NTILT,.true., foil_tilt )
                else
                    if (rank==0) then
                        print *,""
                        write(*,fmt='(a,f10.6)') " selecting bright reflections for n_g = ",n_g
                        print *,"^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^"
                        print *,""
                    end if
                    hkl_bright = floor( hkl * n_g )
                    call selectFoilTilt( imb,hkl_bright,tilt_max,NTILT,.true., foil_tilt )
                    hkl_bright = ceiling( hkl * n_g )
                    call selectFoilTilt( imb,hkl_bright,tilt_max,NTILT,.true., foil_tilt2 )
                    foil_tilt = quaternionToRotMat( slerp( Quaternion_ctor(foil_tilt),Quaternion_ctor(foil_tilt2), n_g - int(n_g) ) )
                end if            
                if (rank==0) print *,""
                call applyFoilTilt(imb,foil_tilt)
                if (rank==0) then
                    call perfectLatticeIntensity(imb,hkl_bright,maxIntensity,I_g)
                    write(*,fmt=form) "tiltFoil info - perfect lattice intensity ",hkl_bright," after tilt ",I_g
                end if
            end if
            
            if (dark) then
            !   add an additional tweak to make hkl extra dark
                call selectFoilTilt( imb,hkl,tilt_max*TWEAK,NTILT,.false., foil_tilt )
                call applyFoilTilt(imb,foil_tilt)
                if (rank==0) then
                    call perfectLatticeIntensity(imb,hkl,maxIntensity,I_g)
                    write(*,fmt=form)  "tiltFoil info - perfect lattice intensity ",hkl," after tilt to dark field ",I_g
                end if
            end if

            return
        end subroutine tiltFoil


        subroutine outputDiffractionPattern( imb,prefix,Nx )
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      output the diffraction pattern as a .png
            type(IntegrateManyBeams),intent(inout)              ::      imb
            character(len=*),intent(in)                         ::      prefix
            integer,intent(in)                                  ::      Nx
            real(kind=real64),dimension(:,:),allocatable        ::      img
            type(Gvectors)                  ::      gv
            integer                         ::      ix,iy,ii
            integer                         ::      ng
            real(kind=real64),dimension(:,:),allocatable        ::      gg
            real(kind=real64)               ::      modg
            real(kind=real64)               ::      spotRadius,ss
            real(kind=real64),dimension(2)  ::      dx
            real(kind=real64),dimension(:),allocatable          ::      I_g

        !---    get the g vectors
            gv = getGvectors(imb)
            ng = getn(gv)
            ii = whichg(gv,(/0.0d0,0.0d0,0.0d0/))       !   is [000] one of the set?
            if (ii == LIB_GVECTORS_UNSET) then
                allocate(gg(3,0:nG))
                gg(:,0) = 0
                do ii = 1,ng
                    gg(:,ii) = getG(gv,ii)
                end do
            else
                allocate(gg(3,1:nG))                
                do ii = 1,ng
                    gg(:,ii) = getG(gv,ii)
                end do
            end if


        !---    find the largest magnitude g-vector
            modg = -huge(1.0)
            do ii = lbound(gg,dim=2),nG
                modg = max(modg,norm2(gg(:,ii)))
            end do

        !---    find the scaling of the image
            spotRadius = Nx * sqrt(0.1d0/(nG*PI))
            ss = ( Nx/2 - 2*spotRadius )/modg
            spotRadius = min( spotRadius , modg*ss/4 )

        !---    construct the image            
            allocate(img(0:Nx-1,0:Nx-1))
            allocate(I_g(0:nG))
            call getIntensity(imb,I_g)
            
            do ii = lbound(gg,dim=2),nG
                print *,"spot ",ii,gethkl(gv,ii),gg(:,ii),I_g(ii)
            end do

            img = 0
            do iy = 0,Nx-1
                do ix = 0,Nx-1
                    do ii = lbound(gg,dim=2),nG
                        dx(1:2) = (/ix,iy/) - ( gg(1:2,ii)*ss + Nx/2 )
                        if (norm2(dx)<=spotRadius) then
                            img(ix,iy) = I_g(ii)
                        end if
                    end do
                end do
            end do

        !---    write image
            if (rank==0) print *,"run_SimulatedTEM::outputDiffractionPattern info - drawing ",(nG+1)," spots to """//trim(prefix)//".diff.png"//""""
            call writePng( trim(prefix)//".diff.png",img )
            return
        end subroutine outputDiffractionPattern



        subroutine errorExit(message)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
            character(len=*),intent(in),optional             ::      message
#ifdef MPI            
            integer             ::      ierror
#endif
            if (rank==0) then
                if (present(message)) print *,"run_SimulatedTEM::"//trim(message)
            end if
#ifdef MPI            
            call MPI_FINALIZE(ierror)
#endif          
            stop
        end subroutine errorExit
        
    end program run_SimulatedTEM
