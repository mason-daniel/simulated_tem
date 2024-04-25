    program radiallyAveragedPowerFunction
!---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
!*      compute the radially averaged power function for an input .png file
!*      return the first and second moments of the power function - which will give a characteristic size and "error bar"
!*      note: this program works only in pixel units - you'll need to convert to real space units yourself
!*      Daniel Mason, UKAEA 2021

        use Lib_FFTW3f                      !   fast fourier transform
        use Lib_Png                         !   read a png
        use Lib_CommandLineArguments        !   generate -help and read in parameters
        use Lib_Exteriors                   !   handle periodic boundary conditions.
        use Lib_RAPFanalysis
        use Lib_RidlerCalvard
        use iso_fortran_env
        implicit none
        

    !---    command line options
        type(CommandLineArguments)      ::      cla
        character(len=256)              ::      file = ""                           !       .png filename to read
        character(len=256)              ::      outfile = ""                        !       .png filename to read
        integer                         ::      Nq = LIB_CLA_NODEFAULT_I            !       number of q radius divisions 
        real(kind=real64)               ::      r_min = LIB_CLA_NODEFAULT_R         !       minimum radius  
        real(kind=real64)               ::      r_max = LIB_CLA_NODEFAULT_R         !       maximum radius  
        logical                         ::      pbc      = .true.                   !       assume periodic boundary conditions
        logical                         ::      optable  = .true.                   !       output C(q) 
        logical                         ::      opfft    = .false.                  !       output fft
        logical                         ::      dbg      = .false.
        logical                         ::      negative = .false.
        logical                         ::      cprime   = .false.                  !       use C'(q), where 2 pi C'(q) = C(q)
         
        
    !---    details about the file        
        real(kind=real64),dimension(:,:),allocatable        ::      img_in        
        real(kind=real64),dimension(:,:,:),allocatable      ::      img3d_in        
        integer                                             ::      Nx,Ny,Nz
        logical                                             ::      file2d
        
    !---    physics-based variables     
        integer                                             ::      border          !   defines how the exterior of the image is handled
        real(kind=real64),dimension(:,:),allocatable        ::      img_ext         !   a larger image, including exterior region
        real(kind=real64),dimension(:),allocatable          ::      rpf,rps         !   radial power function
        real(kind=real64),parameter                         ::      PI = 3.141592654d0 
        real(kind=real64)                                   ::      q_min,q_max     !   min/max q vector length, corresponding to max/min radius
        real(kind=real64)                                   ::      rpfsum          !   int rpf(q) d2q                     
        real(kind=real64)                                   ::      qrpfbar         !   int |q| rpf(q) d2q/int rpf(q) d2q  
        real(kind=real64)                                   ::      q2rpfbar        !   int |q^2| rpf(q) d2q/int rpf(q) d2q
        real(kind=real64)                                   ::      rbar,rstd       !   mean and std dev radius
        real(kind=real64)                                   ::      r2bar           !   mean square radius 
        real(kind=real64)                                   ::      b,t,f           !   Ridler Calvard background,threshold,foreground
        real(kind=real64)                                   ::      bstd,fbar,snr,fot         !   signal to noise ratio, fraction over threshold
        
        
        real(kind=real64),dimension(:,:,:),allocatable      ::        ft3d
    !---    dummy variables
        integer                                     ::      ii , ix,iy,iz , Mx,My,Mz , jx,jy
        real(kind=real64)                           ::      aa 
        real(kind=real64)                           ::      rr , qq  
        character(len=256)                          ::      dummy
        logical                                     ::      ok
        
    !---    read in command line arguments
        cla = CommandLineArguments_ctor(20)         !   allocates space for up to 10 arguments
        call setProgramDescription( cla, "radiallyAveragedPowerFunction.exe \n finds radially averaged power spectrum from FFT" )
        call setProgramVersion( cla, "1.0" )

        call get( cla,"f",file ,LIB_CLA_REQUIRED,"         input filename (.png or .dat format)" )  
        outfile=file
        call get( cla,"o",outfile ,LIB_CLA_OPTIONAL,"       output filename prefix ( defaults to input file )" )  
        if (.not. hasArgument(cla,"o")) outfile = file
        
        call get( cla,"N",Nq ,LIB_CLA_OPTIONAL,"       number of q radius divisions" )  
        call get( cla,"r_min",r_min ,LIB_CLA_OPTIONAL,"   minimum radius" )  
        call get( cla,"r_max",r_max ,LIB_CLA_OPTIONAL,"   maximum radius" )  
        call get( cla,"pbc",pbc ,LIB_CLA_OPTIONAL,"     assume periodic boundary conditions in input image" )  
        call get( cla,"optable",optable ,LIB_CLA_OPTIONAL,"     output c(q) table" )  
        call get( cla,"opfft",opfft ,LIB_CLA_OPTIONAL,"     output fft image" )  
        call get( cla,"dbg",dbg ,LIB_CLA_OPTIONAL,"     debug mode" )  
        call get( cla,"negative",negative ,LIB_CLA_OPTIONAL,"     read file in, scale 0:1 then take negative [1:0] " )  
        call get( cla,"cprime",cprime ,LIB_CLA_OPTIONAL,"     use C'(q), where 2 pi q C'(q) = C(q) " )  
        
        call report(cla)              
        if (.not. allRequiredArgumentsSet(cla)) stop
        if (hasHelpArgument(cla)) stop
        call delete(cla)
                
        
        if (dbg) then
            FFTW_DBG = .true.
            LIB_RANDC_DBG = .true.
            RAPFanalysis_dbg = .true.
        end if
        
    !---    input the image
        ii = index( file,".png",back = .true. )
        if (ii>0) then
            call readPng( file,img_in )
            Nx = size(img_in,dim=1)
            Ny = size(img_in,dim=2)
            Nz = 1
            file2d = .true.
        else
            ok = .false.
            open(unit=600,file=trim(file),action="read")
            do
                read(unit=600,fmt='(a)',iostat=ii) dummy
                if (ii/=0)   exit
             
                                        
                dummy = adjustl(dummy) ; if (dummy(1:1) == "#") cycle                    
                    
                read(dummy,fmt=*,iostat=ii) Nx,Ny,Nz
                if (ii==0) then
                    print *,"3d .dat file ",Nx,Ny,Nz
                    file2d = .false. ; ok = .true.
                else
                    read(dummy,fmt=*,iostat=ii) Nx,Ny
                    if (ii==0) then
                        print *,"2d .dat file ",Nx,Ny
                        Nz = 1 ; file2d = .true. ; ok = .true.
                    end if
                end if              
                exit  
            end do
            if (.not. ok) then
                print *,"radiallyAveragedPowerFunction error - could not read file """//trim(file)//""""
                print *,"expecting generic 2d/3d dat format"
                print *," ""# one or more comment lines "" "
                print *,"Nx Ny [Nz]"
                print *,"data[:,1] or data[:,1,1]"
                print *,"data[:,2] or data[:,2,1]"
                print *,"etc"
                stop
            end if
            if (file2d) then
                allocate(img_in(Nx,Ny))
                do iy = 1,Ny
                    read(600,fmt=*) img_in(:,iy)
                end do
            else
                allocate(img3d_in(Nx,Ny,Nz))
                do iz = 1,Nz
                    do iy = 1,Ny
                        read(600,fmt=*) img3d_in(:,iy,iz)
                    end do
                end do
                pbc = .false. ; opfft = .false.
                print *,"note: -pbc, -opfft options not coded for 3d .dat file format"
            end if
            close(unit=600)
        end if   
        
        if (negative) then
            if (file2d) then
                !   use q_max and q_min as handy dummies. 
                !   note all I;m doing is rescaling max value to 0 and min value to 1.
                q_max = maxval(img_in)
                q_min = minval(img_in)
                q_min = 1/max(1.0d-16,q_max-q_min)                
                img_in = (q_max-img_in)*q_min
            else
                q_max = maxval(img3d_in)
                q_min = minval(img3d_in)
                q_min = 1/max(1.0d-16,q_max-q_min)                
                img3d_in = (q_max-img3d_in)*q_min
            end if
        end if
        
        
    !---    what to do with the image exterior?
    !   note that a fast fourier transform assumes that the input is periodic. 
    !   If it is not, then it will generate artefacts at the boundaries
    !   So give the option to add a smoothly varying blur to the edges of the image
        if (file2d) then
            if (.not. pbc) then
                print *,""
                print *,"radiallyAveragedPowerFunction info - adding exterior"
             !   border = max(Nx,max(Ny,Nz))/2 /4          
             !   allocate(img_ext(-2*border:Nx-1+2*border,-2*border:Ny-1+2*border))
             !   aa = sum(img_in)/(Nx*Ny)
             !   Exterior_dbg = .true.
             !   call setExterior( img_in,EXTERIOR_BLURAVG,border,img_ext(-border:Nx-1+border,-border:Ny-1+border),aa )
             !   deallocate(img_in)
             !   allocate(img_in(-border:Nx-1+border,-border:Ny-1+border))
             !   img_in = img_ext(-border:Nx-1+border,-border:Ny-1+border)
             !   call setExterior( img_in,EXTERIOR_AVERAGE,border,img_ext,aa  )
             !    
             !    
             !   Nx = Nx + 4*border
             !   Ny = Ny + 4*border
                
             !    border = max(Nx,max(Ny,Nz))/2 /4     
             !    allocate(img_ext(-border:Nx-1+border,-border:Ny-1+border))
             !    call setExterior( img_in,EXTERIOR_BLURPBC,border,img_ext )
             !    
             !    
             !    Nx = Nx + 2*border
             !    Ny = Ny + 2*border
                
                allocate(img_ext(0:2*Nx-1,0:2*Ny-1))
                img_ext = 0.0d0            
                
                Mx = int((Nx+1)/2) 
                My = int((Ny+1)/2) 
                
                do iy = 0,Ny-1
                    do ix = 0,Nx-1                         
                        img_ext(ix+Mx,iy+My) = img_in(ix+1,iy+1)
                    end do
                    
                    do ix = 0,Mx-1                                                       
                        img_ext(ix      ,iy+My) = img_in( Mx-ix,iy+1 )
                        img_ext(ix+Nx+Mx,iy+My) = img_in( Nx-ix,iy+1 )
                    end do
                end do            
                 
                do iy = 0,My-1
                
                    do ix = 0,Nx-1
                        img_ext(ix+Mx,iy) = img_in( ix+1,My-iy)
                        img_ext(ix+Mx,iy+Ny+My) = img_in( ix+1,Ny-iy)
                    end do
                    
                    do ix = 0,Mx-1
                        img_ext(ix      ,iy)       = img_in( Mx-ix,My-iy)
                        img_ext(ix+Nx+Mx,iy)       = img_in( Nx-ix,My-iy)
                        img_ext(ix+Nx+Mx,iy+Ny+My) = img_in( Nx-ix,Ny-iy)
                        img_ext(ix      ,iy+Ny+My) = img_in( Mx-ix,Ny-iy)
                    end do
                    
                    
                end do                                    
                
                
                
                
                
                Nx = Nx * 2
                Ny = Ny * 2
                
                !if (len_trim(outfile)/=0) &                         
                !call writePng( trim(outfile)//".ext.png" , img_ext )
                
            else            
                allocate(img_ext(0:Nx-1,0:Ny-1))     
                img_ext = img_in
            end if
        end if
        
    !---    determine the range of the radial power function
        if (Nq == LIB_CLA_NODEFAULT_I) Nq = max(Nx,max(Ny,Nz))/2         
        if (r_min == LIB_CLA_NODEFAULT_R) r_min = 2
        if (r_max == LIB_CLA_NODEFAULT_R) r_max = max(Nx,max(Ny,Nz))/2
        q_max = 2*PI/r_min
        q_min = 2*PI/r_max
        
        
    !---    friendly message    
        print *,"radiallyAveragedPowerFunction"
        print *,"   input file  "//trim(file)
        if (file2d) then
            print *,"   Nx,Ny       ",Nx,Ny
        else
            print *,"   Nx,Ny,Nz    ",Nx,Ny,Nz
        end if
        print *,"   r_min,q_max ",r_min,q_max
        print *,"   r_max,q_min ",r_max,q_min        
        print *,"   N           ",Nq
        print *,"   pbc         ",pbc
        if (len_trim(outfile)/=0) &
        print *,"   output file "//trim(outfile)
        print *,""

        
                   
    !---    compute the radial power function         
        print *,"radial power function"
        allocate(rpf(0:Nq))
        
        if (file2d) then
            !aa = sum(img_ext)/(size(img_ext,dim=1)*size(img_ext,dim=2))         !   average value of the image 
            !img_ext = img_ext - aa                                               !   remove the average value
                                                                                  
            call findImageIntensityFeatures( img_in, b,t,f,bstd,fbar,snr,fot )
        
            call radialPowerFunction( img_ext ,q_min,q_max, rpf,rpfsum,qrpfbar,q2rpfbar,snr )
        else
            aa = sum(img3d_in)/(Nx*Ny*Nz)         !   average value of the image 
            img3d_in = img3d_in - aa                                              !   remove the average value
                                                      
            call radialPowerFunction( img3d_in ,q_min,q_max, rpf,rpfsum,qrpfbar,q2rpfbar  )
            
            
       ! !---    DEBUG - TO DELETE
       !     Mx = int(Nx/2)                      !   note real array gives complex FT with N/2 output values! 
       !     My = int(Ny/2)
       !     Mz = int(Nz/2)            
       !     allocate(ft3d(-Mx:Mx,-My:My,-Mz:Mz))
       !     call FFT3d( img3d_in,ft3d )
       !     
       !     call writePng("test_x3d.png",log(ft3d(0,:,:)),normalise=.true. )
       !     call writePng("test_y3d.png",log(ft3d(:,0,:)),normalise=.true. )
       !     call writePng("test_z3d.png",log(ft3d(:,:,0)),normalise=.true. )
            
            
        end if    
        
    !---    output the radial power function         
        if (optable) then
            write(*,fmt='(a8,100a24)') "# ","|q|","|r|","C(q)","q C(q)" ,"q^2 C(q)" ,"q^3 C(q)"
            do ii = 0,Nq
                 
                qq = ii*q_max / Nq
                rr = 2*3.141592654d0 * qq * rpf(ii)
                 
                 
                if (ii == 0) then
                    write(*,fmt='(i8,100g24.12)') ii,qq,0.0d0,rpf(ii),rr , qq*rr,qq*qq*rr                   
                else
                    write(*,fmt='(i8,100g24.12)') ii,qq,(2*PI)/qq,rpf(ii),rr ,qq*rr,qq*qq*rr
                end if
                             
            end do
        end if        
    !---    output the integrated quantities
        print *,"<in>                     = ",aa
        print *,"integral 2 pi |q| rpf(q) = ",rpfsum
        print *,"< rpf(q) >               = ",rpfsum/(Nx*Ny*Nz)
        print *,"<q>                      = ",qrpfbar 
        print *,"<q2>                     = ",q2rpfbar
        print *,"single Gaussian sigma    = ",sqrt(PI)/(2*qrpfbar)
        
        call findMeanAndStdDev(qrpfbar,q2rpfbar,rbar,rstd,r2bar)
        print *,"lognormal sigma          = ",rbar," +/- ",rstd
                
        
    !---    find some more image properties
        if (file2d) then
            print *,"Ridler and Calvard b,t,f = ",b,t,f
            print *,"fraction over threshold  = ",fot," = ",fot*Nx*Ny," px"
            print *,"signal to noise ratio    = ",snr
            print *,"background noise,<f>     = ",bstd,fbar
            print *,"Area over thresh/pi<r^2> = ",fot*Nx*Ny/(PI*r2bar)
        end if         
    !---    bye bye
        
        
        if (opfft) then
            if (.not. pbc) then   
                deallocate(img_in)                                   
                allocate(img_in(size(img_ext,dim=1),size(img_ext,dim=2)))
                img_in = img_ext
            end if             
            
            
            
            
            if (len_trim(outfile)/=0) then
                outfile = trim(outfile)//".ft.png"
            else
                outfile = trim(file)//".ft.png"
            end if
            
            
            Nx = int(Nx/2)        !   note real array gives complex FT with N/2 output values!           
            Ny = int(Ny/2)        !               
            deallocate(img_ext)
                         
            allocate(img_ext(-Nx:Nx,-Ny:Ny))
            
            
             
            call FFT2d( img_in,img_ext )
            img_ext = log(img_ext)
            call writePng(outfile,img_ext,normalise=.true. )
            
        end if        
        
        
        if (cprime) then
            allocate(rps(0:Nq))
            call radialPowerSpectrum( img_ext ,q_min,q_max, rps,rpf )
            print *,"radial power spectrum ",sum(rps),sum(rpf),sqrt(PI)*sum(rps)/sum(rpf)
            deallocate(rps)
        end if                   
        deallocate(rpf)
        
        
        print *,""
        print *,"done"
        print *,""
        
    end program radiallyAveragedPowerFunction
        
      