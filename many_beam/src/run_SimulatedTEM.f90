
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

!*      good break points to set (GDB)
!*      b run_SimulatedTEM.f90:274            !   after command line params
!*      b run_SimulatedTEM.f90:296            !   after initialisation
!*      b run_SimulatedTEM.f90:352            !   after foil tilt



       
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
        use Lib_ComputePhaseFactor
        use Lib_ManyBeam
#ifdef MPI
        use mpi_f08
#endif
        implicit none


        character(len=8),parameter          ::      VERSION = "0.0.1"


        real(kind=real64),parameter         ::      PI = 3.14159265390d0
        integer                             ::      rank = 0,nProcs = 1

        integer,parameter                   ::      CLA_FILE   = 1
        integer,parameter                   ::      CLA_MICRO  = 2
        integer,parameter                   ::      CLA_DIFF   = 3
        integer,parameter                   ::      CLA_CALC   = 4
        integer,parameter                   ::      CLA_OUTPUT = 5

    !---    command line parameter input        
        type(CommandLineArguments)          ::      cla
        
        character(len=256)                  ::      filename = ""           !   input filename
        character(len=256)                  ::      outfile = ""            !   output filename prefix
        character(len=8)                    ::      latticename = UNKNOWN_LATTICE
        real(kind=real64)                   ::      T = 300.0d0             !   temperature (K)
        real(kind=real64)                   ::      V = 200.0d0             !   accelerator voltage (kV)
        real(kind=real64),dimension(3)      ::      a0_in = 0.0             !   indicative lattice parameter (A)
        integer                             ::      nPrec = 1               !   precession angles
        real(kind=real64)                   ::      precAngle = 0.003d0     !   precession angle (rad)
        logical                             ::      opPng = .true.          !   should we automatically output a png file
        logical                             ::      columnar = .true.       !   use columnar approximation
        logical                             ::      lossy = .false.         !   use imaginary parts of crystal structure factor
        logical                             ::      twobeam = .false.       !   two-beam mode
        integer,dimension(:),allocatable    ::      hkl                     !   the g-vector for the aperture. allocatable because its either hkl or hjkl 
        logical                             ::      yflip = .true.          !   Right-hand rule for png output - y=0 at bottom.
        logical                             ::      pngBlur = .true.        !   blurring radius sigma/a
        logical                             ::      pngNorm = .true.        !   normalise output intensity
        real(kind=real64)                   ::      n_g = 1                 !   which reflection to make bright
        logical                             ::      df = .true.             !   find dark field        
        real(kind=real64)                   ::      tilt_max = 5.0          !   maximum stage tilt (deg)
        logical                             ::      opDiffPatt = .true.     !   output diffraction pattern
        integer                             ::      ntilt = 500             !   number of tilt angles searched for good diffraction pattern
        real(kind=real64)                   ::      tweak = 0.001d0         !   tweak in tilt angle permitted to optimise diffraction pattern
        real(kind=real64)                   ::      imin = 0.001d0          !   intensity threshold for "important" beam
        logical                             ::      dbg = .false.   

    !---    physical variables
        type(IntegrateManyBeams)            ::      imb
        integer                             ::      Nx,Ny           !   image size
        real(kind=real64),dimension(:,:),allocatable        ::      img
        integer,dimension(4)                ::      dat
        character(len=256)                  ::      dummy
        integer,dimension(:,:),allocatable  ::      hkl_important

    !---    dummy variables        
        integer             ::      ierror
        integer             ::      ii,ix,iy,iyp
        real(kind=real64)   ::      dd

    !---    timing
        integer,parameter               ::      TIMER = kind(1_int32)
        integer(kind=TIMER),parameter   ::      T_INIT      = 1_TIMER
        integer(kind=TIMER),parameter   ::      T_TILT      = 2_TIMER
        integer(kind=TIMER),parameter   ::      T_OUTPUT    = 3_TIMER
        integer(kind=TIMER),parameter   ::      T_INTEGRATE = 4_TIMER
        integer(kind=TIMER),parameter   ::      T_TOTAL     = 5_TIMER
        type(Callipers),dimension(T_TOTAL)      ::      tt




!******************************************************************************    
!*    
!*      START
!*    
!*    
!******************************************************************************    
    

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
        call setCategories(cla,(/   "file handling         ",       &
                                    "microscope settings   ",       &
                                    "diffraction conditions",       &
                                    "calculation           ",       &
                                    "output                " /))
        call setProgramDescription( cla, "run_SimulatedTEM" )
        call setProgramVersion( cla, VERSION )   
          
    
        call get( cla,"f",filename ,LIB_CLA_REQUIRED,"            input filename",CLA_FILE )           
        ii=0 ; call get( cla,"a0",a0_in,ii ,LIB_CLA_OPTIONAL,"         lattice parameter(s)",CLA_FILE )
        if (ii==1) a0_in(2:3) = a0_in(1)
        call get( cla,"lattice",latticename ,LIB_CLA_OPTIONAL,"    lattice type",CLA_FILE )            
        outfile = filename                                       
        call get( cla,"o",outfile ,LIB_CLA_OPTIONAL,"          output filename prefix",CLA_FILE )           
        if ( (trim(getSuffix(outfile))=="xyz").or.(trim(getSuffix(outfile))=="lammps").or.(trim(getSuffix(outfile))=="lmp").or.(trim(getSuffix(outfile))=="dat").or.(trim(getSuffix(outfile))=="png") ) then
            outfile = trim(removeSuffix(outfile))
        end if
        
        
        call get( cla,"T",T ,LIB_CLA_OPTIONAL,"          temperature (K)",CLA_MICRO )                                                        
        call get( cla,"V",V ,LIB_CLA_OPTIONAL,"          accelerating voltage (kV)",CLA_MICRO )                                                        
        call get( cla,"theta",tilt_max,LIB_CLA_OPTIONAL,"      maximum foil tilt angle (deg)",CLA_MICRO )           
        call get( cla,"ntilt",ntilt,LIB_CLA_OPTIONAL,"      number of foil tilt angles considered",CLA_MICRO )           
        call get( cla,"tweak",tweak,LIB_CLA_OPTIONAL,"      fine-tune tweak angle as fraction of theta",CLA_MICRO )           


        ii=0; dat = LIB_CLA_NODEFAULT_I; call get( cla,"g",dat,ii ,LIB_CLA_REQUIRED,"            g-vector reflection to output",CLA_DIFF )
        if ( (ii==3).or.(ii==4) ) then
            allocate(hkl(ii))
            hkl(1:ii) = dat(1:ii)       
        end if
        call get( cla,"ng",n_g ,LIB_CLA_OPTIONAL,"         tilt foil until this reflection bright",CLA_DIFF )                 
        call get( cla,"nPrecAngle",nPrec ,LIB_CLA_OPTIONAL," number of precession images to take",CLA_DIFF )                                                        
        call get( cla,"precAngle",precAngle ,LIB_CLA_OPTIONAL,"  precession angle (rad)",CLA_DIFF )             
        df = (n_g /= 1.0d0)
        call get( cla,"dark",df ,LIB_CLA_OPTIONAL,"       tweak foil tilt to find close dark field",CLA_DIFF )                 

        call get( cla,"col",columnar ,LIB_CLA_OPTIONAL,"        use columnar approximation",CLA_CALC )                 
        call get( cla,"loss",lossy ,LIB_CLA_OPTIONAL,"       use complex crystal structure factor for inelastic electron scattering",CLA_CALC )    
        call get( cla,"2",twobeam ,LIB_CLA_OPTIONAL,"          compute in two beam approximation",CLA_CALC )                 
        call get( cla,"imin",imin,LIB_CLA_OPTIONAL,"       intensity threshold for beam selection",CLA_CALC )           



        call get( cla,"png",opPng ,LIB_CLA_OPTIONAL,"        output result as .png file",CLA_OUTPUT)                           
        call get( cla,"blur",pngBlur ,LIB_CLA_OPTIONAL,"       blur .png file by smoothing width sigma",CLA_OUTPUT )                                                        
        call get( cla,"norm",pngNorm ,LIB_CLA_OPTIONAL,"       normalise .png file intensity",CLA_OUTPUT )                                                        
        call get( cla,"yflip",yflip ,LIB_CLA_OPTIONAL,"      output png with y=0 at bottom, to match ovito orientation z out of screen",CLA_OUTPUT )                 
        call get( cla,"opDP",opDiffPatt ,LIB_CLA_OPTIONAL,"       output diffraction pattern",CLA_OUTPUT )                 
        call get( cla,"dbg",dbg ,LIB_CLA_OPTIONAL,"        additional debug output",CLA_OUTPUT )                 

    !---    check for input sanity        
        if (hasArgument(cla,"imin") .and. twobeam) then
            if (rank==0) then
                print *,"run_SimulatedTEM warning - argument conflict - attempt to set -2 and -imin parameters"
                print *,"run_SimulatedTEM warning - overriding two beam mode"
            end if
            twobeam = .false.
        end if
        if ((rank==0) .and. .not. hasArgument(cla,"f")) print *,"error - expected argument '-f <filename>'"
        if ((rank==0) .and. .not. allocated(hkl)) print *,"error - expected argument '-g h,k,l' or '-g h,i,k,l'"


        if (rank==0) call report(cla)
        if (hasHelpArgument(cla)) call errorExit("done")
        if (.not. allRequiredArgumentsSet(cla)) call errorExit("error - required arguments unset")
        call delete(cla)
        

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
            if (len_trim(dummy)==0) then
                call get_environment_variable("HOST",dummy)
                if (len_trim(dummy)==0) dummy="localhost"           !   ??
            end if
            print *,"   running on ",nProcs," processors on "//trim(dummy)
            print *,""
            print *,"   file handling"
            print *,"   ^^^^^^^^^^^^^"
            print *,"       input atom positions    """//trim(filename)//""""
            print *,"       output file prefix      """//trim(outfile)//""""
            if (trim(latticename)==UNKNOWN_LATTICE) then
                print *,"       lattice type            TBD"
            else
                print *,"       lattice type            """//trim(latticename)//""""
            end if
            if (any(a0_in<=0)) then
                print *,"       indicative latt param   TBD"
            else
                print *,"       indicative latt param   ",a0_in
            end if
            print *,""

            print *,"   microscope settings"
            print *,"   ^^^^^^^^^^^^^^^^^^^"
            print *,"       accelerator voltage     ",V," (kV)"
            print *,"       temperature             ",T," (K)"
            print *,"       max foil tilt angle     ",tilt_max," (deg)"
            print *,"       number of foil tilts    ",ntilt
            print *,"       tweak foil tilt         ",tweak*tilt_max," (deg)"
            print *,""

            print *,"   diffraction condition"
            print *,"   ^^^^^^^^^^^^^^^^^^^^^"
            print *,"       output g-vector         ",hkl
            print *,"       number precession imgs  ",nPrec
            print *,"       precession angle        ",precAngle," (rad)"
            print *,"       bright reflection n_g = ",n_g
            print *,"       dark field tweak?       ",df
            print *,""

            print *,"   calculation settings"
            print *,"   ^^^^^^^^^^^^^^^^^^^^"
            print *,"       use columnar approx     ",columnar
            print *,"       inelastic scattering    ",lossy
            print *,"       two beam mode           ",twobeam
            if (.not. twobeam) &
            print *,"       beam intensity thresh   ",imin            
            print *,""

            print *,"   output settings"
            print *,"   ^^^^^^^^^^^^^^^"
            print *,"       output .png             ",oppng," with blur? ",pngBlur," normalised? ",pngNorm
            print *,"       output diff patt        ",opDiffPatt 
            print *,"       debug mode              ",dbg
            print *,""  
        end if





!******************************************************************************    
!*    
!*      CONSTRUCTOR
!*    
!*    
!******************************************************************************    

        if (rank==0) then
            print *,""
            print *,"initialising"
            print *,"^^^^^^^^^^^^"
        end if

        if (dbg) then
            LIB_IMB_DBG = dbg
            !Lib_ComputePhaseFactor_DBG = dbg
        end if

        call Lib_IntegrateManyBeams_init_MPI()
        !call pause(tt(T_INIT))


    !---    constructor - read input file and determine parameters of calculation
        !tt(T_CTOR) = Callipers_ctor()
        imb = IntegrateManyBeams_ctor( latticename,a0_in,filename,T,V, nPrec = nPrec, precAngle = precAngle, columnar = columnar, lossy = lossy, rho_in = max(4.0d0,norm2(n_g*hkl)+4)  )
        call pause(tt(T_INIT))

        
        if ((rank==0) .and. opDiffPatt) then
            call outputDiffractionPattern( imb,trim(outfile)//"_before",1024 )
        end if

    !---    tilt the foil to "good" diffraction conditions
        
        if (rank==0) then
            print *,""
            print *,"finding diffraction condition"
            print *,"^^^^^^^^^^^^^^^^^^^^^^^^^^^^^"
        end if
        tt(T_TILT) = Callipers_ctor()        
        call tiltFoil(imb,tilt_max*PI/180,hkl,n_g,df)
        call pause(tt(T_TILT))


        tt(T_OUTPUT) = Callipers_ctor()
        call pause(tt(T_OUTPUT))
        if ((rank==0) .and. opDiffPatt) then
            call start(tt(T_OUTPUT))
            call outputDiffractionPattern( imb,trim(outfile)//"_ready",1024 )
            call pause(tt(T_OUTPUT))
        end if
        
        call start(tt(T_TILT))
        if (twobeam) then
            call changeGvectors( imb , hkl_in = hkl ) 
        else
            call selectImportantBeams( imb, imin = imin, hkl = hkl_important )
            call changeGvectors( imb , rule = LIB_IMB_BLOCK_HIGH_SG , hkl_in = hkl_important )
        end if

        call setImagingSpace(imb)                        
        call pause(tt(T_TILT))


    ! !---    compute the phase field x = exp[ - ig.u(r) ]        
    !     tt(T_PHASE) = Callipers_ctor()
        
    !     call pause(tt(T_PHASE))



        if (rank==0) then
            print *,""
            print *,"run parameters"
            print *,"^^^^^^^^^^^^^^"
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


        
        
        if (opDiffPatt) then
            call start(tt(T_OUTPUT))
            call outputDiffractionPattern( imb,trim(outfile)//"_after",1024 )
            call pause(tt(T_OUTPUT))
        end if

        ! if (opDiffPatt) then
        ! !    call outputDiffractionPattern( imb,trim(outfile),1024 )
        ! end if

        if (opPng) then
            call start(tt(T_OUTPUT))
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


                if ((rank==0).and. pngBlur) call blurImg( img, 2/getA( imb ) )
 


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
            call pause(tt(T_OUTPUT))
        end if
        


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
            write(*,fmt='(a,f16.6)') "   tilt foil      ",elapsed(tt(T_TILT))
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
            integer             ::      ix,iy,ik,kernel_width
            real(kind=real64),dimension(:),allocatable      ::      img_tmp
            real(kind=real64),dimension(:),allocatable      ::      kernel
            real(kind=real64)   ::      ff,ww,i2s2

            Nx = size(img,dim=1)
            Ny = size(img,dim=2)
            allocate(img_tmp(0:max(Nx,Ny)-1))
            
            write(*,fmt='(a,f8.3,a)') " run_SimulatedTEM::blurImg info - blur ",sigmaona," (px)"
            kernel_width = ceiling(3*sigmaona)            
            allocate(kernel(-kernel_width:kernel_width))
            i2s2 = 1/(2*sigmaona*sigmaona)
            do ik = -kernel_width,kernel_width
                kernel(ik) = exp( -ik*ik*i2s2 )
            end do

        !---    blur in x direction for each row iy
            do iy = 0,Ny-1
                img_tmp(0:Nx-1) = img(0:Nx-1,iy)
                do ix = 0,Nx-1
                    ff = 0 ; ww = 0
                    do ik = max(0,ix-kernel_width),min(Nx-1,ix+kernel_width)
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
                    do ik = max(0,iy-kernel_width),min(Ny-1,iy+kernel_width)
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
    !*      ng [hkl] is bright. Optionally then tweak until the reflection [hkl] is dark

            type(IntegrateManyBeams),intent(inout)              ::      imb
            real(kind=real64),intent(in)                        ::      tilt_max
            integer,dimension(:),intent(in)                     ::      hkl
            real(kind=real64),intent(in)                        ::      n_g
            logical,intent(in)                                  ::      dark
            real(kind=real64),dimension(3,3)    ::      foil_tilt

            logical                             ::      MillerBravais 
            integer,dimension(size(hkl))        ::      hkl_bright
            real(kind=real64)                   ::      I_g
            real(kind=real64),dimension(3,3)    ::      foil_tilt2
            character(len=16)                   ::      form
            logical                             ::      maxIntensity = .false.      !   select tilt based on foil exit intensity
             


            MillerBravais = (size(hkl)==4)
            form = "(a,3i4,a,f16.8)"
            if (MillerBravais) form = "(a,4i4,a,f10.6)"

            if (rank==0) then
                call perfectLatticeIntensity(imb,hkl,maxIntensity,I_g)
                write(*,fmt=form) " tiltFoil info - perfect lattice intensity ",hkl," before tilt ",I_g
            end if

            if (tilt_max<=0) return
 
    
            if (.not. ( (n_g==1.0).and.df) ) then
                if (abs(n_g - nint(n_g)) < 1.0d-6) then
                    if (rank==0) then
                        print *,""
                        write(*,fmt=form) " selecting bright reflection g,n_g = ",hkl,",",n_g
                        print *,""
                    end if
                    hkl_bright = nint( hkl * n_g )
                    call selectFoilTilt( imb,hkl_bright,tilt_max,NTILT,.true., foil_tilt )
                else
                    if (rank==0) then
                        print *,""
                        write(*,fmt=form) " selecting bright reflections for g,n_g = ",hkl,",",n_g
                        print *,""
                    end if
                    hkl_bright = floor( hkl * n_g )
                    call selectFoilTilt( imb,hkl_bright,tilt_max,NTILT,.true., foil_tilt )
                    hkl_bright = ceiling( hkl * n_g )
                    call selectFoilTilt( imb,hkl_bright,tilt_max,NTILT,.true., foil_tilt2 )
                    foil_tilt = quaternionToRotMat( slerp( Quaternion_ctor(foil_tilt),Quaternion_ctor(foil_tilt2), ceiling(n_g) - n_g ) )       !   ceiling(x)-x = 1 for x = 3.00 and = 0 for 3.99
                end if            
                if (rank==0) print *,""
                call applyFoilTilt(imb,foil_tilt)
                if (rank==0) then
                    call perfectLatticeIntensity(imb,hkl_bright,maxIntensity,I_g)
                    write(*,fmt=form) " tiltFoil info - perfect lattice intensity ",hkl_bright," after tilt ",I_g
                end if
            end if
            
            if (dark) then
            !   add an additional tweak to make hkl extra dark
                call selectFoilTilt( imb,hkl,tilt_max*TWEAK,NTILT,.false., foil_tilt )
                call applyFoilTilt(imb,foil_tilt)
                if (rank==0) then
                    call perfectLatticeIntensity(imb,hkl,maxIntensity,I_g)
                    write(*,fmt=form)  " tiltFoil info - perfect lattice intensity ",hkl," after tilt to dark field ",I_g
                end if
            end if

            return
        end subroutine tiltFoil

        subroutine selectImportantBeams( imb, imin, hkl )
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      select a beam as "important" if its intensity ever gets above imin
    !*      allocate and return the array of reflections hkl
            type(IntegrateManyBeams),intent(in)                 ::      imb
            real(kind=real64),intent(in)                        ::      imin
            integer,dimension(:,:),allocatable,intent(out)      ::      hkl

            type(Gvectors)                  ::      gv
            real(kind=real64)               ::      i_important
            real(kind=real64),dimension(:),allocatable          ::      I_g
            logical                         ::      MillerBravais 
            integer                         ::      nG
            integer                         ::      kk,nn
            character(len=32)               ::      form
             


        !---    get the set of g-vectors            
            gv = getGvectors(imb)
            nG = getn(gv)
            MillerBravais = isMillerBravais(gv)

        !---    find max intensity of beam inside the (assumed perfect) foil            
            allocate(I_g(0:nG))
            call getIntensity(imb,I_g,maxIntensity = .true.)

        !---    count the important beams
            nn = 0
            do kk = 0,nG
                if (I_g(kk)>=imin) nn = nn + 1
            end do

        !---    is this a sensible number?             
            if (nn<2) then
                kk = maxloc( I_g(1:),dim=1 )            !   given that I_g(0) = I_g(hkl=[000]) is important, what is number 2?
                i_important = I_g(kk) * 0.50d0          !   I'll say that anything within 50% of the number 2 beam is important. 
                if (rank == 0) then
                    print *,"selectImportantBeams info - did not find any beams (except hkl=[000]) over minimum intensity threshold ",imin
                    if (MillerBravais) then
                        print *,"    second most important beam is [hjkl] = ",gethjkl(gv,kk)," with intensity ",I_g(kk)
                    else
                        print *,"    second most important beam is [hkl] = ",gethkl(gv,kk)," with intensity ",I_g(kk)
                    end if
                    print *,"    using intensity threshold ",i_important
                end if
            end if
                    

        !---    allocate memory for important reflections            
            if (MillerBravais) then
                allocate(hkl(4,nn))
            else
                allocate(hkl(3,nn))
            end if
            hkl = 0

        !---    store the important reflections. The first is hkl=[000], which is not in gv.
            form = "(a,i4,a,3i4,a,f16.8)"
            if (MillerBravais) form = "(a,i4,a,4i4,a,f10.6)"

            nn = 1
            if (rank == 0) write (*,fmt=form) " selectImportantBeams info - selected beam ",nn," reflection ",hkl(:,nn)," intensity ",I_g(0)
            do kk = 1,nG
                if (I_g(kk)>=imin) then
                    nn = nn + 1
                    if (MillerBravais) then
                        hkl(:,nn) = gethjkl(gv,kk)
                    else
                        hkl(:,nn) = gethkl(gv,kk)
                    end if
                    if (rank == 0) write (*,fmt=form) " selectImportantBeams info - selected beam ",nn," reflection ",hkl(:,nn)," intensity ",I_g(kk)
                end if
            end do
            

            return
        end subroutine selectImportantBeams

        subroutine outputDiffractionPattern( imb,prefix,Nx )
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      output the diffraction pattern as a .png
            type(IntegrateManyBeams),intent(inout)              ::      imb
            character(len=*),intent(in)                         ::      prefix
            integer,intent(in)                                  ::      Nx
            real(kind=real64),dimension(:,:),allocatable        ::      img
            type(Gvectors)                  ::      gv
            integer                         ::      ix,iy,ii
            integer                         ::      nG
            real(kind=real64),dimension(:,:),allocatable        ::      gg
            real(kind=real64)               ::      modg
            real(kind=real64)               ::      spotRadius,ss
            real(kind=real64),dimension(2)  ::      dx
            real(kind=real64),dimension(:),allocatable          ::      I_g
            real(kind=real64),parameter     ::      MIN_GREYSCALE = 1.0d0/256

        !---    get the g vectors
            gv = getGvectors(imb)
            nG = getn(gv)
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


        !---    find the intensity of the beams at the foil exit            
            allocate(I_g(0:nG))
            call getIntensity(imb,I_g,maxIntensity = .false.)
            if ((rank==0).and.dbg) then
                do ii = lbound(gg,dim=2),nG
                    print *,"spot ",ii,gethkl(gv,ii),gg(:,ii),I_g(ii)
                end do
            end if

        !---    find a decent scaling of the image
            spotRadius = Nx * sqrt(0.1d0/(nG*PI))
            ss = ( Nx/2 - 2*spotRadius )/modg
            spotRadius = min( spotRadius , modg*ss/4 )
            if (rank==0) print *,"run_SimulatedTEM::outputDiffractionPattern info -  minmax I_g ",minval(I_g),maxval(I_g)


        !---    construct the image            
            allocate(img(0:Nx-1,0:Nx-1))
            img = 0
            do iy = 0,Nx-1
                do ix = 0,Nx-1
                    do ii = lbound(gg,dim=2),nG
                        dx(1:2) = (/ix,iy/) - ( gg(1:2,ii)*ss + Nx/2 )
                        if (I_g(ii) >= MIN_GREYSCALE) then
                            if (norm2(dx)<=spotRadius)   img(ix,iy) = I_g(ii) 
                        else
                            if (norm2(dx)<=spotRadius/2) img(ix,iy) = MIN_GREYSCALE
                        end if
                    end do
                end do
            end do

        !---    write image
            if (rank==0) then
                print *,"run_SimulatedTEM::outputDiffractionPattern info - drawing ",(nG+1)," spots to """//trim(prefix)//".diff.png"//""""
                print *,"run_SimulatedTEM::outputDiffractionPattern info -  minmax img ",minval(img),maxval(img)
                call writePng( trim(prefix)//".diff.png",img )
            end if
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
            if (rank==0) print *,""
            stop
        end subroutine errorExit
        
    end program run_SimulatedTEM
