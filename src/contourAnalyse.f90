!   gfortran -ffree-line-length-256 -O2 ${MYF90LIB}/NBAX_StringTokenizers.f90 ${MYF90LIB}/Lib_Callipers.f90 ${MYF90LIB}/Lib_LocalFilter2d.f90 ${MYF90LIB}/Lib_CommandLineArguments.f90 ${MYF90LIB}/Lib_Png.f90 ${MYF90LIB}/Lib_GaussianBlurs.f90 ${MYF90LIB}/Lib_Exteriors.f90  ${MYF90LIB}/Lib_MarchingSquares.f90 contourAnalyse.f90 -o contourAnalyse.exe -llapack

    program contourAnalyse
!---^^^^^^^^^^^^^^^^^^^^^^
        use iso_fortran_env
        use Lib_LowPassFilter3d
        use Lib_MaxlikelihoodFilter
        use Lib_CommandLineArguments      
        use NBAX_StringTokenizers
        use Lib_Exteriors
        use Lib_Filenames
        use Lib_MarchingSquares
        use Lib_ColourScale
        use Lib_GaussianBlurs
        use Lib_ImageCorrelationFunction
        use Lib_RidlerCalvard
        use Lib_LocalFilter2d
        use Lib_Quicksort
        use Lib_Callipers
        use Lib_Png     
        !use Lib_Sinc
        implicit none


    !---    required input variables
        character(len=256)                  ::      filename = ""
        character(len=256)                  ::      outfile_prefix = "[input file prefix]"
        logical                             ::      negative = .false.
        real(kind=real64)                   ::      sigma = LIB_CLA_NODEFAULT_R
        real(kind=real64)                   ::      lambda = LIB_CLA_NODEFAULT_R
        integer                             ::      nLevels = 0          
        real(kind=real64)                   ::      iso_min = 0.0d0    
        real(kind=real64)                   ::      iso_max = 1.0d0     
        integer                             ::      nIso_spec = 0
        real(kind=real64),dimension(1000)   ::      iso_spec = 0 
        real(kind=real64)                   ::      flattenLinearThresh = LIB_CLA_NODEFAULT_R
        real(kind=real64)                   ::      lengthperpx = 1.0d0
        logical                             ::      flattenLinear = .false.           
        logical                             ::      ignoreLinear = .false.           
        logical                             ::      iso_auto = .false.
        logical                             ::      op_tmp = .true.
        logical                             ::      op_contour = .false.
        logical                             ::      op_contourPositive = .false.
        logical                             ::      pbc = .false.
        integer                             ::      nCorrelationFunctionBins = LIB_CLA_NODEFAULT_I
        real(kind=real64)                   ::      correlationFunctionrMax = LIB_CLA_NODEFAULT_R
        real(kind=real64)                   ::      localFilterSigma = LIB_CLA_NODEFAULT_R
        logical                             ::      FTfilter = .false.
        integer,dimension(2)                ::      nFourierCoeffs = LIB_CLA_NODEFAULT_I
        logical                             ::      removeAvg = .false.
        logical                             ::      opLPF = .false.
        logical                             ::      opHPF = .false.
        logical                             ::      normOp = .false.
        logical                             ::      rgbOp = .false.
        logical                             ::      absOp = .false.
        logical                             ::      sqrtOp = .false.
        logical                             ::      sandp = .false.
        logical                             ::      quiet = .false.
        logical                             ::      mf = .false.
        logical                             ::      denoisefilt = .false.
        integer                             ::      nhstripe = 0
        integer                             ::      nvstripe = 0
        integer                             ::      pspot_nBins = LIB_CLA_NODEFAULT_I
        real(kind=real64)                   ::      pspot_max = LIB_CLA_NODEFAULT_R
        
        integer,dimension(128)              ::      hstripe = LIB_CLA_NODEFAULT_I
        integer,dimension(128)              ::      vstripe = LIB_CLA_NODEFAULT_I
!        logical                             ::      strainToTEM = .false.  
!        logical                             ::      TEMtoStrain = .false.        
        integer                             ::      nRepeat = 1
        integer,dimension(2)                ::      dbg_px = -1
        character(len=8)                    ::      VERSION = "1.1.0"
        integer                             ::      colourScale
        type(CommandLineArguments)          ::      cla
        
    !---    constants
        real(kind=real64),parameter         ::      PI = 3.141592654d0        
        
        
    
    !---    image
        real(kind=real64),dimension(:,:),allocatable        ::      img,img_blur,img_iso,img_flatline,img_cfn
    
    
    !---    contour analysis
        real(kind=real64)                                   ::      isolevel, area,length,curv,curv2 , areap,lengthp,curvp,curv2p , length2onareabar
        integer                                             ::      nLines,nPaths
        real(kind=real64),dimension(:,:),allocatable        ::      line
        integer,dimension(:),allocatable                    ::      pathLen
        real(kind=real64),dimension(:,:,:),allocatable      ::      path
        integer,dimension(:,:),allocatable                  ::      indx

    !---    correlation function
        real(kind=real64),dimension(:),allocatable          ::      correlationFunction
        real(kind=real64)                                   ::      delta
        
    !---    timing
        integer,parameter   ::      T_TOTAL   = 0
        integer,parameter   ::      T_READ    = 1
        integer,parameter   ::      T_BLUR    = 2
        integer,parameter   ::      T_FT      = 3
        integer,parameter   ::      T_CORR    = 4
        integer,parameter   ::      T_ANALYSE = 5
        type(Callipers),dimension(0:5)       ::  tt
    
    !---    dummy
        integer             ::      nx,ny ,ix,iy , jx,jy
        type(GaussianBlur)  ::      gb
        integer             ::      ii,jj,kk,mm , oldnp,np, newnp
        logical             ::      ok
        real(kind=real64)   ::      amax,lmax,cmax,c2max , apbar,dpbar,lpbar,a2pbar 
        real(kind=real64)   ::      aa,ll,bb,cc,cc2 , th,ff,bstd,fbar,snr,fot , iso , sqrtap,a2p
        character(len=256)  ::      dummy
        real(kind=real64),dimension(:,:),allocatable        ::      x,dat
        real(kind=real64),dimension(:,:,:),allocatable      ::      img1
        real(kind=real64),dimension(:,:,:,:),allocatable    ::      qdat
        real(kind=real64),dimension(2)                      ::      xmax
        real(kind=real64)               ::      img0avg          
!        real(kind=real64)               ::      xig = LIB_CLA_NODEFAULT_R,sg = LIB_CLA_NODEFAULT_R,gg = LIB_CLA_NODEFAULT_R,zmax = LIB_CLA_NODEFAULT_R
!        real(kind=real64),dimension(4)  ::      dummy_r
        real(kind=real64),dimension(9)  ::      medianFilterLocal3x3
        
        integer,dimension(:),allocatable                    ::      nspots
        real(kind=real64),dimension(:,:),pointer            ::      pspotarea,psa_tmp
        real(kind=real64),dimension(:),allocatable          ::      pspot_hist 
        real(kind=real64),dimension(:,:),allocatable        ::      nspot_hist 
        real(kind=real64),dimension(:),allocatable          ::      peakHeight
        integer,dimension(:),allocatable                    ::      nEnclosure,nAreap
        
        
        
        
        
        
        tt(T_TOTAL) = Callipers_ctor()
        
        
        
    !---    input command line parameters
        cla = CommandLineArguments_ctor(50)
        call setProgramDescription( cla, "contourAnalyse.exe "    &
            //"\n    performs filtering and contour analysis of a .png file \n    or .dat file with format \n    # comment lines \n    Nx Ny \n    dat(1:Nx,1) \n    dat(1:Nx,2) \n    ..." )
        call setProgramVersion( cla, VERSION )
        call setCategories( cla,(/  "input/output    ",     &
                                    "filtering       ",     &
                                    "correlation fn  ",     &
                                    "contour analysis" /) )
        call get( cla,"f",filename ,LIB_CLA_REQUIRED,"                    filename",1 )  
        call get( cla,"o",outfile_prefix ,LIB_CLA_OPTIONAL,"                  output filename",1 )          
        if (.not. hasArgument(cla,"o")) outfile_prefix = removeSuffix(filename)
        call get( cla,"negative",negative ,LIB_CLA_OPTIONAL,"           take negative (.png only)",1 )          
        
        call get( cla,"sigma",sigma ,LIB_CLA_OPTIONAL           ,"              Gaussian blur radius",2 )  
        call get( cla,"lambda",lambda ,LIB_CLA_OPTIONAL           ,"             max likelihood curvature lengthscale ",2 )  
        call get( cla,"sandp",sandp ,LIB_CLA_OPTIONAL,           "              remove salt and pepper noise",2 )  
        call get( cla,"mf",mf ,LIB_CLA_OPTIONAL,           "                 apply 3x3 median filter",2 )  
        call get( cla,"hstripe",hstripe,nhstripe ,LIB_CLA_OPTIONAL,           "            remove horizontal stripes",2 )  
        call get( cla,"vstripe",vstripe,nvstripe ,LIB_CLA_OPTIONAL,           "            remove vertical stripes",2 )  
        call get( cla,"localFilter",localFilterSigma ,LIB_CLA_OPTIONAL,"        set range for max likelihood filter ( eg 1.0 )",2 )  
        call get( cla,"denoise",denoisefilt ,LIB_CLA_OPTIONAL,"            iterative denoise using max likelihood filter ( eg 1.0 )",2 )  
        
        call get( cla,"n",nFourierCoeffs,ii,LIB_CLA_OPTIONAL,"                  nFourierCoeffs to remove/retain",2 ) 
        FTfilter = hasArgument(cla,"n")  
        if (FTfilter .and. (ii == 1))  nFourierCoeffs(2) = nFourierCoeffs(1)
        
        call get( cla,"removeAvg",removeAvg ,LIB_CLA_OPTIONAL,"          remove average",2 ) 
        call get( cla,"opLPF",opLPF ,LIB_CLA_OPTIONAL,"              output low pass filter",2 ) 
        call get( cla,"opHPF",opHPF ,LIB_CLA_OPTIONAL,"              output high pass filter",2 ) 
        call get( cla,"norm",normop ,LIB_CLA_OPTIONAL,"               normalise output pngs to 0:1",2 ) 
        call get( cla,"rgb",dummy ,LIB_CLA_OPTIONAL,"               output rgb png",2 ) 
        
        
        if (hasArgument(cla,"rgb")) then
            colourScale = getColourScale(dummy)
            if (colourScale == COLOURSCALE_UNSET) then
                print *,"contourAnalyse.exe error - did not recognise colour scale name """//trim(dummy)//""" choose from"
                call listAvailableColourScales()
                stop
            end if
            rgbOp = .true.
        end if
        
        
        
        
        
        
        
        
        call get( cla,"q",quiet ,LIB_CLA_OPTIONAL,"               suppress most normal output",2 ) 
        
        
!        call get( cla,"abs",absop ,LIB_CLA_OPTIONAL,"                prefilter f=abs(f-<f>)",2 ) 
!        call get( cla,"sqrt",sqrtop ,LIB_CLA_OPTIONAL,"               prefilter f=sqrt(abs(f-<f>))",2 ) 
!        call get( cla,"strainToTEM",strainToTEM ,LIB_CLA_OPTIONAL,"               strain to TEM ",2 ) 
!        ii= 4 ; call get( cla,"TEM",dummy_r ,ii,LIB_CLA_OPTIONAL,"               strain<->TEM parameters (xig,sg,|g|,zmax)",2 ) 
!        if (hasArgument( cla,"TEM" )) then
!            if (ii/=4) then
!                stop "-TEM parsing error - expect parameters (xig,sg,|g|,zmax) in A units"
!            else    
!                xig  = dummy_r(1)
!                sg   = dummy_r(2)
!                gg   = dummy_r(3)
!                zmax = dummy_r(4)
!            end if
!        end if
!        call get( cla,"TEMtostrain",TEMtoStrain ,LIB_CLA_OPTIONAL,"                TEM to strain",2 ) 
        
        call get( cla,"flattenLinearThresh",flattenLinearThresh ,LIB_CLA_OPTIONAL,"flatten image where features circularity over this level (eg 4.0)",2 )  
        flattenLinear = hasArgument(cla,"flattenLinearThresh")
        call get( cla,"ignoreLinearThresh",flattenLinearThresh ,LIB_CLA_OPTIONAL,"flatten image where features circularity over this level (eg 4.0)",2 )  
        ignoreLinear = hasArgument(cla,"flattenLinearThresh")
        
        
        ii = 2 ; call get( cla,"dbg_px",dbg_px ,ii,LIB_CLA_OPTIONAL,"             debug pixel",2 )
        if (hasArgument( cla,"dbg_px" )) then
            if (ii/=2) then
                stop "-dbg_px parsing error - expect parameters (x,y) "
            else    
                LIB_LF2D_DBGX = dbg_px(1)
                LIB_LF2D_DBGY = dbg_px(2)
            end if
        end if
        
        
        call get( cla,"cf_nbins",nCorrelationFunctionBins,LIB_CLA_OPTIONAL         ,"           number of correlation function bins",3 )  
        call get( cla,"cf_rmax",correlationFunctionRmax,LIB_CLA_OPTIONAL         ,"            correlation function max radius (default N/2 px)",3 )  
        
        call get( cla,"nLevels",nLevels ,LIB_CLA_OPTIONAL         ,"            number of isolevels (set 0 for off)",4 )  
        call get( cla,"iso_min",iso_min ,LIB_CLA_OPTIONAL         ,"            minimum isolevel",4 )  
        call get( cla,"iso_max",iso_max ,LIB_CLA_OPTIONAL         ,"            maximum isolevel",4 )  
        call get( cla,"iso",iso_spec,nIso_spec ,LIB_CLA_OPTIONAL         ,"            specify isolevels",4 )  
        call get( cla,"iso_auto",iso_auto ,LIB_CLA_OPTIONAL         ,"           auto set isolevels",4 )  
        call get( cla,"pbc",pbc ,LIB_CLA_OPTIONAL,"                periodic boundary conditions ( warning - currently buggy if shapes thread cell )",4 )  
        call get( cla,"nRepeat",nRepeat ,LIB_CLA_OPTIONAL,"            pseudo periodic boundary image repeats ( workaround for buggy pbc )",4 )  
        call get( cla,"op_tmp",op_tmp ,LIB_CLA_OPTIONAL          ,"             output temporary images",4 )  
        call get( cla,"op_contour",op_contour ,LIB_CLA_OPTIONAL      ,"         output contour plot",4 )  
        call get( cla,"op_contour+",op_contourPositive ,LIB_CLA_OPTIONAL      ,"         output contour plot, only for +ve area contours",4 )  
        call get( cla,"s",lengthPerPx ,LIB_CLA_OPTIONAL,"               scale (eg nm/px)",4 ) 
        call get( cla,"nBin",pSpot_nBins ,LIB_CLA_OPTIONAL,"               number of bins for histogram",4 ) 
        call get( cla,"maxD",pSpot_max   ,LIB_CLA_OPTIONAL,"               maximum diameter in histogram",4 ) 

        
        if (rgbop .and. .not. (hasArgument(cla,"iso_min") .and. hasArgument(cla,"iso_max")) ) iso_auto = .true.
        
        
        if (.not. quiet) call report(cla)
        if (.not. allRequiredArgumentsSet(cla)) stop
        if (hasHelpArgument(cla)) stop
        call delete(cla)
        
        if ( (nLevels>0).and.(nIso_spec>0) ) then
            print *,"contourAnalyse.exe warning - have specified evenly spaced iso levels _and_ specific iso levels."
            print *,"                             revert to using specified only!"
            nLevels = 0
        end if
        
        if ( (nIso_spec==0).and.(nLevels>0) ) then
            nIso_spec = nLevels
            if (nLevels == 1) then
                iso_spec(1) = iso_min
            else
                do ii = 1,nLevels
                !---    set isolevel
                    iso_spec(ii) = iso_min + (iso_max-iso_min)*(ii-1)/(nLevels-1)
                end do
            end if
        end if
        
        
        if ( (nIso_spec==0).and.(flattenLinear) ) then
            print *,"contourAnalyse.exe error - have specified flatten image, but need defined isolevels"
            stop
        end if
        
        if ( (op_contourPositive).and.(op_contour) ) then
            print *,"contourAnalyse.exe warning - -op_contour implies -op_contour+"
            op_contourPositive = .false.
        end if
        
        
        
        absop = absop .or. sqrtop
        
    !---    output input 
        if (.not. quiet) then
            print *,""
            print *,"contourAnalyse.exe v"//trim(VERSION)
            print *,"^^^^^^^^^^^^^^^^^^^^"//repeat("^",len_trim(VERSION))
            print *,""
        
            print *,"filename               : "//trim(filename)
            print *,"output filename prefix : "//trim(outfile_prefix)
            if (sigma==LIB_CLA_NODEFAULT_R) then
                print *,"gaussian blur          : no"   
            else
                print *,"gaussian blur sigma    : ",sigma
            end if
            if (lambda==LIB_CLA_NODEFAULT_R) then
                print *,"log likelihood filt    : no"   
            else
                print *,"log likelihood lambda  : ",lambda
            end if
            
            print *,"nLevels                : ",nLevels
            print *,"iso_spec               : ",iso_spec(1:nIso_spec)
            print *,"pbc? nRepeat           : ",pbc,nRepeat
            if (iso_auto) then
                print *,"iso_auto               : ",iso_auto
            else
                print *,"iso_min                : ",iso_min
                print *,"iso_max                : ",iso_max
            end if
            print *,"length per px          : ",lengthPerPx
            print *,"prefilter f=abs(f-<f>) : ",absop
            print *,"prefilter f=sqrt(abs(f-<f>)) : ",sqrtop
            print *,"normalise .png [0:1]   : ",normop     
            print *,"op_tmp                 : ",op_tmp    
            print *,"op_contour             : ",op_contour," +ve only? ",op_contourPositive     
            if (nCorrelationFunctionBins==LIB_CLA_NODEFAULT_I) then
                print *,"correlation fn         : no"   
            else
                print *,"correlation fn bins    : ",nCorrelationFunctionBins    
                if (correlationFunctionrMax == LIB_CLA_NODEFAULT_R) then
                    print *,"correlation fn rMax    : N/2 px (tbd)"                 
                else
                    print *,"correlation fn rMax    : ",correlationFunctionrMax    
                end if
            end if
            print *,"sandp                  : ",sandp    
            print *,"3x3 median filter      : ",mf
            if (nhstripe > 0) &
            print *,"hstripe                : ",hstripe(1:nhstripe)    
            if (nvstripe > 0) &
            print *,"vstripe                : ",vstripe(1:nvstripe)        
            if (localFilterSigma==LIB_CLA_NODEFAULT_R) then
                print *,"localFilter            : no"   
            else
                print *,"localFilter            : ",localFilterSigma    
            end if
            if (flattenLinear) then
                print *,"flattenLinearThresh    : ",flattenLinearThresh
            else
                print *,"flattenLinearThresh    : no"
            end if
            if (FTfilter) then
                print *,"FT filter, n coeffs    : ",nFourierCoeffs
                print *,"    remove avg         : ",removeAvg
                print *,"    opLPF              : ",opLPF
                print *,"    opHPF              : ",opHPF
            end if
    !        if (strainToTEM) then
    !            print *,"attempting strainToTEM"
    !        end if
    !        if (TEMtoStrain) then
    !            print *,"attempting TEMtoStrain"
    !        end if
                
            print *,""
        end if
        
    !---    check for nonsense input
        !filename = removePngExtension(filename)
        ! call setClockpath()
        
        inquire(file=trim(filename),exist=ok)
        if (.not. ok) then
            print *,"contourAnalyse error - file """//trim(filename)//""" not found"
            stop
        end if
        
            
        if (trim( getSuffix(filename) )/="png") then
            tt(T_READ) = Callipers_ctor()
            if (.not. quiet) then
                print *,""
                print *,"reading input file"
                print *,""
                !print *,"contourAnalyse warning - file """//trim(filename)//".png"" not found"
            end if
    !         inquire(file=trim(filename),exist=ok)
    !         if (.not. ok) then
    !             print *,"contourAnalyse error - file """//trim(filename)//""" not found"
    !             stop
    !         end if            
            open(unit=700,file=trim(filename),action="read")
            
!                if (strainToTEM .or. TEMtoStrain) then
!                    do ii = 1,1000
!                        read(unit=700,fmt='(a)') dummy            
!                        dummy = adjustl(dummy)
!                        
!                        jj = index(dummy,"s_g (1/A)")
!                        if (jj/=0) then
!                            call parse( dummy(jj+9:),sg,ok )
!                            if (ok) then
!                                print *,"strainToTEM info - s_g = ",sg,"(1/A)"
!                            end if
!                        end if
!                        
!                        jj = index(dummy,"d2bi_xi (A)")
!                        if (jj/=0) then
!                            call parse( dummy(jj+11:),dummy_r,ix )
!                            if (ix == 2) then
!                                xig = dummy_r(2)                                
!                                print *,"strainToTEM info - xig = ",xig,"(A)"
!                            end if
!                        end if
!                        
!                        jj = index(dummy,"foil (A)")
!                        if (jj/=0) then
!                            call parse( dummy(jj+8:),dummy_r,ix )
!                            if (ix == 3) then   
!                                zmax = dummy_r(3)
!                                print *,"strainToTEM info - zmax = ",zmax,"(A)"
!                            end if
!                        end if
!                        
!                        jj = index(dummy,"d2bi_g")
!                        if (jj/=0) then
!                            call parse( dummy(jj+6:),dummy_r,ix )
!                            if (ix == 4) then
!                                gg = dummy_r(4)
!                                print *,"strainToTEM info - |g| = ",gg,"(1/A)"
!                            end if
!                        end if
!                        if (dummy(1:1) == "#") cycle
!                        exit
!                    end do
!                    if (.not. (strainToTEM .or. TEMtoStrain)) stop "contourAnalyse error - failed to read required parameters s_g,xi_g,zmax from file"
!                    if (any( (/xig,sg,zmax,gg/) == LIB_CLA_NODEFAULT_R)) stop "contourAnalyse error - failed to read required parameters s_g,xi_g,zmax"
!            
!                else
                    do ii = 1,1000
                        read(unit=700,fmt='(a)') dummy            
                        dummy = adjustl(dummy)
                        if (dummy(1:1) == "#") cycle
                        exit
                    end do
!                end if
                read(dummy,fmt=*,iostat=jj) nx,ny
                if (jj/=0) then
                    print *,"contourAnalyse error - file header expects comment lines ""#"" then nx,ny"
                    stop 
                end if
                allocate(img(1:nx,1:ny))
                read(700,fmt=*,iostat=jj) img
                if (jj/=0) then
                    print *,"contourAnalyse error - error reading image data."
                    stop 
                end if
            close(unit=700)
            call pause(tt(T_READ))
        else
    !---    read in .png file
            tt(T_READ) = Callipers_ctor()
            if (.not. quiet) then
                print *,""
                print *,"reading input .png file"
                print *,""                
            end if
            ! call setClockPath()            
            call readPng(trim(filename),img,negative = negative)            
            img = max(0.0d0,min(1.0d0,img))
            filename = removePngExtension(filename)
            nx = size(img,dim=1)
            ny = size(img,dim=2)
            call pause(tt(T_READ))
        end if
        
        if (.not. quiet) then
            print *,""
            print *,"Ridler-Calvard analysis before filtering"
            print *,""
        end if
        call findImageIntensityFeatures( img, bb,th,ff,bstd,fbar,snr,fot )
        if (.not. quiet) then
            write (*,fmt='(a,3f12.5)') "   back,thresh,fore ",bb,th,ff
            write (*,fmt='(a,3f12.5)') "   b std.dev, <f>   ",bstd,fbar
            write (*,fmt='(a,3f12.5)') "   snr,fot          ",snr,fot
            print *,"total area (px)   ",Nx*Ny
            print *,"total area (L^2)  ",Nx*Ny*lengthPerPx*lengthPerPx
            print *,""
        end if
        
        
         
        if (absop) then
            aa = sum(img)/(nx*ny)
            img = abs(img-aa)
            if (sqrtop) img = sqrt(img)            
        end if
        
        !if (sigma < 0.5d0) then
        !    print *,"contourAnalyse warning - sigma < 0.5?"                        
        !end if

        
        if (FTfilter) then        
            tt(T_FT) = Callipers_ctor()
            if (.not. quiet) then
                print *,""
                print *,"Fourier transform filter"
                print *,""
            end if
            allocate(x(2,nx*ny))
            allocate(dat(1,nx*ny))
            do iy = 1,ny
                do ix = 1,nx
                    x(1:2,ix+nx*(iy-1)) = (/ix,iy/)
                    dat(1,ix+nx*(iy-1)) = img(ix,iy)
                end do
            end do
            xmax(1:2) = (/nx,ny/)
            allocate(qdat(1,4,0:nFourierCoeffs(1),0:nFourierCoeffs(2)))        
            call generateFourierCoefficients(x,dat,xmax,nFourierCoeffs,qdat)
            deallocate(x)
            deallocate(dat)
            allocate(img1(1,nx,ny))
            
            img0avg = sum(img)/(nx*ny)   
            
            img1 = interpolateFromFourierCoefficients(nx,ny,nFourierCoeffs,qdat)
            deallocate(qdat)
            
            if (opLPF) then
                call writePng( trim(outfile_prefix)//".lpf.png",img1(1,:,:),normalise = normop)
            end if
            img(:,:) = img(:,:) - img1(1,:,:)
            if (.not. removeAvg) img(:,:) = img(:,:) + img0avg
            if (opHPF) then        
                if (.not. normop) then
                    do iy = 1,ny    
                        do ix = 1,nx
                            aa = img(ix,iy)
                            aa = (aa-iso_min)/(iso_max-iso_min)
                            aa = max( 0.0d0,min( 1.0d0,aa ) )
                            img1(1,ix,iy) = aa
                        end do
                    end do
                    call writePng( trim(outfile_prefix)//".hpf.png",img1(1,:,:),normalise = .false. )
                else
                    call writePng( trim(outfile_prefix)//".hpf.png",img(:,:),normalise = .true. )
                end if
            end if                
            deallocate(img1)
            call pause(tt(T_FT))
        end if        
        
        if (iso_auto) then
            iso_min = minval(img)
            iso_max = maxval(img)
            if (.not. quiet) print *,"contourAnalyse info - img data minmax ",iso_min,iso_max
        else
            if (.not. quiet) print *,"contourAnalyse info - img data minmax ",minval(img),maxval(img)
        end if
        
        
        if (nRepeat>1) then
            allocate(dat(0:nx-1,0:ny-1))
            dat = img
            deallocate(img)
            allocate( img(1:nRepeat*nx,1:nRepeat*ny) )
            do iy = 0,nRepeat-1
                do ix = 0,nRepeat-1
                    img( ix*nx+1:ix*nx+nx,iy*ny+1:iy*ny+ny ) = dat(0:nx-1,0:ny-1)
                end do
            end do
            nx = nRepeat*nx
            ny = nRepeat*ny
            deallocate(dat)
        end if
        
        
        if (flattenLinear) then
            allocate(img_flatLine(0:nx-1,0:ny-1))     
            img_flatLine = img
        end if
        
        
        
    !---    blur file
        tt(T_BLUR) = Callipers_ctor()
        if (pbc) then
            allocate(img_blur(0:nx-1,0:ny-1))                
        else
            allocate(img_blur(-1:nx,-1:ny))        
        end if
        if (sandp) then
            if (.not. quiet) then
                print *,""
                print *,"salt and pepper filter"
                print *,"" 
            end if
            if ((localFilterSigma/=LIB_CLA_NODEFAULT_R).and.(localFilterSigma>0.5)) then 
                call setLocalFilter2d_sigma(localFilterSigma)
            else
                call setLocalFilter2d_sigma(0.707d0)
            end if
            call denoise( img,3.0d0,img_blur(0:nx-1,0:ny-1),aa )  
            img = min(1.0d0,max(0.0d0,img_blur(0:nx-1,0:ny-1)))
        end if
                
        if (denoisefilt) then
            if (.not. quiet) then
                print *,""
                print *,"denoise local filter"
                print *,"" 
            end if
            if ((localFilterSigma/=LIB_CLA_NODEFAULT_R).and.(localFilterSigma>0.5)) then 
                call setLocalFilter2d_sigma(localFilterSigma)
            else
                call setLocalFilter2d_sigma(0.707d0)
            end if
            call denoise( img,0.0d0,img_blur(0:nx-1,0:ny-1),aa )  
            img = min(1.0d0,max(0.0d0,img_blur(0:nx-1,0:ny-1)))
        end if        
        
        
        if (mf) then
            if (.not. quiet) then
                print *,""
                print *,"3x3 median filter"
                print *,""
            end if
            do iy = 0,ny-1
                do ix = 0,nx-1
                    ii = 0
                    do jy = max(0,iy-1),min(ny-1,iy+1)
                        do jx = max(0,ix-1),min(nx-1,ix+1)
                            ii = ii + 1
                            medianFilterLocal3x3(ii) = img(jx+1,jy+1)
                        end do
                    end do              
                                           
                    call quicksort( medianFilterLocal3x3(1:ii) )
                    if (mod(ii,2)==1) then
                        img_blur(ix,iy) = medianFilterLocal3x3( (ii+1)/2 )
                    else
                        img_blur(ix,iy) = ( medianFilterLocal3x3( ii/2 ) + medianFilterLocal3x3( ii/2+1 ) )/2
                    end if                        
                    !write (*,fmt='(3i6,a,9f10.4)',advance="no") ix,iy,ii," img ",img(max(1,ix):min(nx,ix+2),max(1,iy):min(ny,iy+2))
                    !write (*,fmt='(a,9f10.4,a,f10.4)') " -> ",medianFilterLocal3x3," = ",img_blur(ix,iy)
                end do
            end do
            img(1:nx,1:ny) = img_blur(0:nx-1,0:ny-1)
        end if 
        
    
        do ix = 1,nhstripe
            if (.not. quiet) then
                print *,""
                print *,"horizontal filter at pixel stripe ",ix," y = ",hstripe(ix) 
                print *,"" 
            end if
            if ((localFilterSigma/=LIB_CLA_NODEFAULT_R).and.(localFilterSigma>0.5)) then 
                call setLocalFilter2d_sigma(localFilterSigma)
            else
                call setLocalFilter2d_sigma(0.707d0)
            end if
            call denoise( img(0:nx-1,max(0,hstripe(ix)-3):min(ny-1,hstripe(ix)+3)),1.5d0,img_blur(0:nx-1,max(0,hstripe(ix)-3):min(ny-1,hstripe(ix)+3)),"V" )  
            img(0:nx-1,max(0,hstripe(ix)-3):min(ny-1,hstripe(ix)+3)) = min(1.0d0,max(0.0d0,img_blur(0:nx-1,max(0,hstripe(ix)-3):min(ny-1,hstripe(ix)+3))))
         end do
     
        
    
        do iy = 1,nvstripe
            if (.not. quiet) then
                print *,""
                print *,"vertical filter at pixel stripe x = ",vstripe(iy) 
                print *,"" 
            end if
            if ((localFilterSigma/=LIB_CLA_NODEFAULT_R).and.(localFilterSigma>0.5)) then 
                call setLocalFilter2d_sigma(localFilterSigma)
            else
                call setLocalFilter2d_sigma(0.707d0)
            end if
            call denoise( img(max(0,vstripe(iy)-3):min(nx-1,vstripe(iy)+3),0:ny-1),1.5d0,img_blur(max(0,vstripe(iy)-3):min(nx-1,vstripe(iy)+3),0:ny-1),"H" )  
            img(max(0,vstripe(iy)-3):min(nx-1,vstripe(iy)+3),0:ny-1) = min(1.0d0,max(0.0d0,img_blur(max(0,vstripe(iy)-3):min(nx-1,vstripe(iy)+3),0:ny-1)))
        end do      
    
        
!        if ((localFilterSigma/=LIB_CLA_NODEFAULT_R).and.(localFilterSigma>0.5)) then        
!            print *,""
!            print *,"local denoise filter file"
!            print *,""
!            call setLocalFilter2d_sigma(localFilterSigma)
!            call localFilter( img,img_blur(0:nx-1,0:ny-1) )
!            img = min(1.0d0,max(0.0d0,img_blur(0:nx-1,0:ny-1)))
!        end if
!        
        
        
        if ( (sigma /= LIB_CLA_NODEFAULT_R).and.(sigma>0) ) then
            if (.not. quiet) then
                print *,""
                print *,"Gaussian blurring file"
                print *,""
            end if
            gb = GaussianBlur_ctor(sigma)
            call report(gb)
            img_blur = 0.0d0
            do iy = 0,ny-1
                do ix = 0,nx-1
                    img_blur(ix,iy) = blur( gb,img,ix,iy,pbc )
                end do
            end do
            call delete(gb)
            img = img_blur(0:nx-1,0:ny-1)
        else
            img_blur = 0.0d0
            img_blur(0:nx-1,0:ny-1) = img(1:nx,1:ny) 
        end if   
        
        if ( (lambda /= LIB_CLA_NODEFAULT_R).and.(lambda>0) ) then
            if (.not. quiet) then
                print *,""
                print *,"log likelihood filtering file"
                print *,""
            end if
            call findImageIntensityFeatures( img, bb,th,ff,bstd,fbar,snr,fot )
            th = bstd/lambda
            
            call maxLikelihoodFilter( img,img_blur(0:,0:), bstd,th )
            
            img = img_blur(0:nx-1,0:ny-1)
        else
            img_blur = 0.0d0
            img_blur(0:nx-1,0:ny-1) = img(1:nx,1:ny) 
        end if   
        
        
                                           
        
        

        if (.not. pbc) then
            Exterior_dbg = .true.
            call setExterior(img,EXTERIOR_LOHI,1,img_blur)
        end if
        
        
        
!        if (strainToTEM) then   
!            print *,""  
!            print *,"constructing TEM image from strain"       
!            print *,""
!             
!            cc = zmax/xig
!            do iy = 0,ny-1
!                do ix = 0,nx-1
!                    aa = sg + img_blur(ix,iy)*gg/(2*PI)
!                    aa = xig*aa
!                    aa = sqrt( 1 + aa*aa )
!                    aa = cc*sinc(cc*aa)
!                    img_blur(ix,iy) = aa*aa
!                end do
!            end do    
!        end if
        
        
        
!        if (TEMtoStrain) then   
!            print *,""  
!            print *,"constructing strain image from TEM"       
!            print *,""
!             
!            cc = xig/zmax
!            do iy = 0,ny-1
!                do ix = 0,nx-1
!                    aa = sqrt( img_blur(ix,iy) )*cc
!                    aa = cc*isinc( aa )
!                    aa = aa*aa - 1
!                    aa = sqrt(aa)/xig - sg
!                    aa = aa*(2*PI/gg)
!                    img_blur(ix,iy) = aa
!                end do
!            end do    
!        end if 
        
        
        
        
        
        if (.not. quiet) then
            print *,""
            print *,"Ridler-Calvard analysis after filtering"
            print *,""
        end if
        call findImageIntensityFeatures( img, bb,th,ff,bstd,fbar,snr,fot )
        if (.not. quiet) then
            write (*,fmt='(a,3f12.5)') "   back,thresh,fore ",bb,th,ff
            write (*,fmt='(a,3f12.5)') "   b std.dev, <f>   ",bstd,fbar
            write (*,fmt='(a,3f12.5)') "   snr,fot          ",snr,fot
            print *,""
        end if        
        
        
        
        
        
        if (op_tmp) then 
            if (rgbop) then
                allocate(img1(3,0:nx-1,0:ny-1))
                
                do iy = 0,ny-1    
                    do ix = 0,nx-1  
                        aa = img_blur(ix,iy)
                        aa = (aa-iso_min)/(iso_max-iso_min)
                        aa = max( 0.0d0,min( 1.0d0,aa ) )
                        img1(:,ix,iy) = getRGB_double( colourScale,aa )
                        
                        
                     !   if ( max(ix,iy)<5 ) then
                     !       print *,ix,iy,aa,int(256*img1(:,ix,iy))
                     !   end if
                        
                    end do
                end do                    
                call write_rgb_png( trim(outfile_prefix)//".png",img1(:,:,:) )
                deallocate(img1)
            else
                if (.not. normop) then
                    allocate(img1(1,0:nx-1,0:ny-1))
                    do iy = 0,ny-1    
                        do ix = 0,nx-1  
                            aa = img_blur(ix,iy)
                            aa = (aa-iso_min)/(iso_max-iso_min)
                            aa = max( 0.0d0,min( 1.0d0,aa ) )
                            img1(1,ix,iy) = aa
                        end do
                    end do
                    call writePng( trim(outfile_prefix)//".png",img1(1,:,:),normalise = .false. )
                    deallocate(img1)
                else 
                    call writePng(trim(outfile_prefix)//".png",img_blur(0:nx-1,0:ny-1),normalise=.true.)
                end if
            end if
        end if
        call pause(tt(T_BLUR))
        
!     !---    construct tiny file from blurred file.
!         print *,""
!         print *,"constructing tiny file"
!         print *,""    
!         sigma = max(1.0d0,sigma)                !   don't want to go bigger than full size- haven't coded interpolation
!         mx = ceiling(nx/sigma)
!         my = ceiling(ny/sigma)
        
        
        if (nCorrelationFunctionBins /= LIB_CLA_NODEFAULT_I) then
            tt(T_CORR) = Callipers_ctor()
            allocate( correlationFunction(0:nCorrelationFunctionBins) )
            if (correlationFunctionrMax == LIB_CLA_NODEFAULT_R) correlationFunctionrMax = max( 1, max( Nx,Ny )/2 )
            
            aa = correlationFunctionrMax / nCorrelationFunctionBins         !   pixel width per bin.
            print *,"contourAnalyse info - radial correlation function pixel width per bin ",aa
            cc = min( 1.0d0, 2/aa )                                         !   no point in having higher resolution in image than 2x pixel width
            ix = nint( nx*cc ) ; iy = nint( ny*cc )
            allocate(img_cfn(0:ix+1,0:iy+1))
            print *,"contourAnalyse info - shrink image to ",ix,"x",iy," (scale factor ",cc,") for radial correlation function"
            call shrinkImage( img_blur(0:nx-1,0:ny-1),ix,iy, 1.0d0,(/0.0d0,0.0d0/), img_cfn )
            call findCorrelationFunction( img_cfn(1:ix,1:iy), pbc, correlationFunctionrMax*cc,nCorrelationFunctionBins, correlationFunction )
            deallocate(img_cfn)
                       
!            call findCorrelationFunction( img_blur(0:nx-1,0:ny-1), pbc, correlationFunctionrMax,nCorrelationFunctionBins, correlationFunction )
            open(unit=902,file=trim(outfile_prefix)//".cfn.out",action="write")
            write(unit=902,fmt='(a)') "# contourAnalyse v."//trim(VERSION)
            delta = correlationFunctionrMax/nCorrelationFunctionBins
            write(unit=902,fmt='(a,i8,2f16.6)') "# nBins,rMax,delta ",nCorrelationFunctionBins,correlationFunctionrMax,delta

        !---    quick post process analysis: min max avg stdev
            aa = huge(1.0) ; amax = -huge(1.0) ; cc = 0.0d0 ; cc2 = 0.0d0
            do iy = 0,ny-1
                do ix = 0,nx-1
                    ll = img_blur(ix,iy)
                    cc = cc + ll
                    cc2 = cc2 + ll*ll
                    aa = min(aa,ll)
                    amax = max(aa,ll)
                end do
            end do
            cc = cc/(nx*ny)
            cc2 = cc2/(nx*ny)
            cc2 = sqrt( cc2 - cc*cc )
            write(unit=902,fmt='(a,4f16.6)') "# min,max,avg,stdev ",aa,amax,cc,cc2
           
                    
                        
        !---    quick post process analysis: find 50% point 
            aa = 0.0d0
            do ii = 0,nCorrelationFunctionBins-1
                if ( (correlationFunction(ii)-0.5d0)*(correlationFunction(ii+1)-0.5d0) < 0) then
                    !   passes through 50% point in this interval
                    aa = ii + (0.5d0-correlationFunction(ii))/(correlationFunction(ii+1)-correlationFunction(ii))
                    exit
                end if
            end do
            if (aa>0)  write(unit=902,fmt='(a,2f16.6)') "# 50% point ",aa*delta            
            do ii = 1,nCorrelationFunctionBins-1
                if ( correlationFunction(ii)>max( correlationFunction(ii-1),correlationFunction(ii+1) ) ) then
                    aa = ( correlationFunction(ii-1) - 2*correlationFunction(ii) + correlationFunction(ii+1) )/2
                    bb = ( correlationFunction(ii+1) - correlationFunction(ii-1) )/2
                    cc = - bb/(2*aa)
                    ll = aa*cc*cc + bb*cc + correlationFunction(ii)
                    write(unit=902,fmt='(a,2f16.6)') "# maximum ",(ii+cc)*delta,ll
                else if ( correlationFunction(ii)<min( correlationFunction(ii-1),correlationFunction(ii+1) ) ) then
                    aa = ( correlationFunction(ii-1) - 2*correlationFunction(ii) + correlationFunction(ii+1) )/2
                    bb = ( correlationFunction(ii+1) - correlationFunction(ii-1) )/2
                    cc = - bb/(2*aa)
                    ll = aa*cc*cc + bb*cc + correlationFunction(ii)
                    write(unit=902,fmt='(a,2f16.6)') "# minimum ",(ii+cc)*delta,ll
                end if
            end do
            ii = 0
            do 
                if ( correlationFunction(ii) == correlationFunction(ii+1) ) then
                    do jj = ii+1,nCorrelationFunctionBins
                        if ( correlationFunction(ii) /= correlationFunction(jj) ) then
                            write(unit=902,fmt='(a,2f16.6)') "# plateau ",(ii+jj-1)*delta/2,correlationFunction(ii)
                            ii = jj-1
                            exit
                        end if
                    end do
                end if
                ii = ii + 1
                if (ii>=nCorrelationFunctionBins-1) exit
            end do
            
            write(unit=902,fmt='(i8)') nCorrelationFunctionBins
            do ii = 0,nCorrelationFunctionBins
                write(unit=902,fmt='(i6,2f16.6)') ii,ii*delta,correlationFunction(ii)
            end do
            close(unit=902)
            call pause(tt(T_CORR))
        end if

        
!         
        
    !---    run through isolevels   
        
        if ( (nLevels>0) .or. (nIso_spec>0) ) then
            tt(T_ANALYSE) = Callipers_ctor()
            print *,""
            print *,"starting marching squares analysis"
            print *,""
            open(unit=901,file=trim(outfile_prefix)//".out",action="write",position="append")
            write(unit=901,fmt='(a)') "# contourAnalyse v."//trim(VERSION)
            if (FTfilter) write(unit=901,fmt='(a,i4)') "# Fourier Series high pass filter n = ",nFourierCoeffs
            write(unit=901,fmt='(a,f12.6)') "# Gaussian blur sigma = ",sigma
            write(unit=901,fmt='(a,f12.6)') "# max likelihood local filter sigma = ",localFilterSigma
            write(unit=901,fmt='(a,f12.6)') "# scale per px ",lengthPerPx            
            if (negative) then
                write(unit=901,fmt='(a,f12.6)') "# NOTE: using negative image"
                write(unit=*,fmt='(a,f12.6)') "NOTE: using negative image"
            end if
            write(unit=901,fmt='(a,l4,a,i4)') "# PBC ",pbc," n repeat = ",nRepeat
            write(unit=901,fmt='(a6,a16,a6,12a16,a6,a16,a12,a6)') "#level","iso"," nPath","max a","max l","max c","max c2","total a","total l","total c","total c2","a+","l(a+)","c(a+)","c2(a+)","n(a+)","<l^2/4 pi a>","npx (f>iso)","new"
            write(*,fmt='(a6,a16,a6,12a16,a6,a16,a12,a6)')        " level","iso"," nPath","max a","max l","max c","max c2","total a","total l","total c","total c2","a+","l(a+)","c(a+)","c2(a+)","n(a+)","<l^2/4 pi a>","npx (f>iso)","new"
            if (op_contour .or. op_contourPositive) then
                if (pbc) then
                    allocate(img_iso(0:nx-1,0:ny-1))                 
                else
                    allocate(img_iso(-1:nx,-1:ny))     
                end if
                img_iso = 0.0d0
                
               ! do ii = 1,max(nIso_spec,nLevels)
               ! 
               ! !---    set isolevel
               !     if (nLevels == 1) then
               !         isolevel = iso_min
               !         iso = isolevel
               !     else if (nLevels > 1) then
               !         isolevel = iso_min + (iso_max-iso_min)*(ii-1)/(nLevels-1)
               !         iso = isolevel
               !     else    
               !         isolevel = iso_spec(ii)
               !         iso = real(ii)/nIso_spec
               !     end if
               !     
               !     !if (negative) iso = 1 - iso 
               !     
               !     where (img_blur >= isolevel)
               !         img_iso = iso/2
               !     end where
               !     
               ! end do    
                    
            end if
    
            oldnp = 0 ; np = 0 ; apbar = 0.0d0 ; dpbar = 0.0d0 ; lpbar = 0.0d0  
            allocate(pspotarea(1000, max(nIso_spec,nLevels))) ; pspotarea = 0.0d0
            allocate(nSpots(max(nIso_spec,nLevels))) ; nSpots = 0
            
            
        !---    test peaks
            allocate(indx(size(img_blur,dim=1),size(img_blur,dim=2)))
            allocate(nAreap(nIso_spec))
            nAreap = 0
            call contiguousPeaks( img_blur,iso_spec(1:nIso_spec),indx,peakHeight,nEnclosure )
!            call writePng(trim(outfile_prefix)//"_indx.png",real(indx,kind=real64),normalise=.true.)
        
            
            
            do ii = 1,nIso_spec
            
            !---    set isolevel
!                if (nLevels == 1) then
!                    isolevel = iso_min
!                    iso = isolevel
!                else if (nLevels > 1) then
!                    isolevel = iso_min + (iso_max-iso_min)*(ii-1)/(nLevels-1)
!                    iso = isolevel
!                else    
                    isolevel = iso_spec(ii)                 !   this is the isolevel we are working on
                    iso = real(ii)/nIso_spec                !   ... but this is the greyscale we will write to a .png ( if one is wanted )
!                end if
            
                
                
               ! if (negative) iso = 1 - iso 
                
                sqrtap = 0.0d0              !   average sqrt area for all +ve enclosures
                a2p = 0.0d0                 !   avg area squared for all +ve enclosures
                
            !---    compute marching squares analysis
                nLines = nLinesNeeded( img_blur,isolevel,.false. )                
                nPaths = 0
                if (nLines>0) then
                    allocate(line(4,nLines))
                    call marchingSquares( img_blur,isolevel,nLines,line )                       
                    
                    
                    
                    
                    if (pbc) then
                        call contiguousLines( nLines,line , nPaths,pathLen,path,Nx,Ny )
                    else
                        call contiguousLines( nLines,line , nPaths,pathLen,path,Nx+2,Ny+2 )
                        line = line - 1 !    because img_blur(-1:Nx,-1:) but marching squares outputs (0:,0:)
                        path = path - 1
                    end if
                    
                    if (op_contour) then
                        
                        do kk = 1,nLines
                            do jj = 0,ceiling(2*sigma)-1
                                ix = int( line(1,kk) + ( line(3,kk)-line(1,kk) )*jj/ceiling(2*sigma) )
                                iy = int( line(2,kk) + ( line(4,kk)-line(2,kk) )*jj/ceiling(2*sigma) )
                                !print *,"line ",ii,kk,jj,ix,iy
                                if ( (ix>=0).and.(ix<nx).and.(iy>=0).and.(iy<ny) ) &
                                img_iso(ix,iy) = iso 
                            end do
                        end do
                        
                    end if
                    
                    deallocate(line)
                    
                end if
                
                
            !---    analyse paths
                area = 0.0d0 ; length = 0.0d0 ; curv = 0.0d0 ; curv2 = 0.0d0
                areap = 0.0d0 ; lengthp = 0.0d0 ; curvp = 0.0d0 ; curv2p = 0.0d0
                length2onareabar = 0.0d0
                amax = -huge(1.0) ; lmax = -huge(1.0) ; jj = 0
                do kk = 1,nPaths
                    call areaAndLengthAndCurvContiguousLine(path(:,:,kk),pathlen(kk),aa,ll,cc,cc2)            
                    
                
                        
                    if (aa>amax) then
                        amax = aa
                        lmax = ll
                        cmax = cc
                        c2max = cc2
                        !jj = kk
                    end if
                    
                    if (aa>0) then
                        !   perimeter l = pi d. area A = pi d^2/4 = l^2 / (4 pi )
                        if ( (.not. ignoreLinear).or.(ll*ll<4*PI*flattenLinearThresh) ) then                        
                            
                            nAreap(ii) = nAreap(ii) + 1
                            areap = areap + aa
                            lengthp = lengthp + ll
                            length2onareabar = length2onareabar + ll*ll/aa
                            curvp = curvp + cc
                            curv2p = curv2p + cc2
                            
                            sqrtap = sqrtap + sqrt(aa)
                            a2p = a2p + aa*aa
                            
                            if ( nSpots(ii) + 1> size(pSpotArea,dim=1) ) then
                                allocate(psa_tmp(nSpots(ii)*2,size(pSpotArea,dim=2) ))
                                psa_tmp(1:nSpots(ii),:) = pSpotArea(1:nSpots(ii),:)
                                psa_tmp(nSpots(ii)+1:,:) = 0.0d0
                                
                                deallocate(pSpotArea)
                                pSpotArea => psa_tmp
                            end if
                            
                            nSpots(ii) = nSpots(ii) + 1
                            pSpotArea( nSpots(ii),ii ) = aa*lengthPerPx
                            
                        end if
                        
                    end if
                    
                    area = area + aa
                    length = length + ll
                    curv = curv + cc
                    curv2 = curv2 + cc2
                    
                    if (flattenLinear) then
                        if ( (aa>0).and.(ll*ll > 4*3.141592654d0*aa*flattenLinearThresh) ) then
                            !print *,"flattenInsidePath ",isolevel,kk,ll*ll/(4*3.141592654d0*aa)," at ", &
                             !   sum(path(1,0:pathLen(kk)-1,kk))/pathLen(kk),sum(path(2,0:pathLen(kk)-1,kk))/pathLen(kk)                              
                             call flattenInsidePath( img_flatline ,pathLen(kk),path(:,0:,kk), isolevel,imgThresh=(iso_max-iso_min)/(nLevels-1) )
                        end if
                    end if
                    
                    if (op_contourPositive .and. (aa>0)) then                        
                     
                        do mm = 0,pathLen(kk)-1
                            do jj = 0,ceiling(2*sigma)-1
                                ix = int( path(1,mm,kk) + ( path(1,mod(mm+1,pathLen(kk)),kk)-path(1,mm,kk) )*jj/ceiling(2*sigma) )
                                iy = int( path(2,mm,kk) + ( path(2,mod(mm+1,pathLen(kk)),kk)-path(2,mm,kk) )*jj/ceiling(2*sigma) )
                                 
                                if ( (ix>=0).and.(ix<nx).and.(iy>=0).and.(iy<ny) ) &
                                img_iso(ix,iy) = iso 
                            end do
                        end do
                        
                    end if
                    
                    
                    
                    
                    
                end do
                
                
                 
                     
                     
                   
                
                
                if (nPaths==0) then
                    amax = 0.0d0
                    lmax = 0.0d0
                    cmax = 0.0d0
                    c2max = 0.0d0
                else
                    deallocate(pathLen)
                    deallocate(path)
                end if
                if (nAreap(ii)>0) length2onareabar = length2onareabar/(4*3.141592654d0*nAreap(ii))
                
            !---    output result
            
                amax = amax*lengthPerPx*lengthPerPx ; area = area*lengthPerPx*lengthPerPx ; areap = areap*lengthPerPx*lengthPerPx
                lmax = lmax*lengthPerPx             ; length = length*lengthPerPx         ; lengthp = lengthp*lengthPerPx           
                cmax = cmax                         ; curv = curv                         ; curvp = curvp                        
                c2max = c2max/lengthPerPx           ; curv2 = curv2/lengthPerPx           ; curv2p = curv2p/lengthPerPx          
                                          
            
                sqrtap = sqrtap*lengthPerPx  
                a2p = a2p*lengthPerPx*lengthPerPx*lengthPerPx*lengthPerPx  
                
                
                !newnp = max(0,nAreap - oldnp)        !   number of new +ve enclosures.
                newnp = max(0,nAreap(ii) - nEnclosure(ii-1))
                apbar = apbar + ( areap/max(1,nAreap(ii)) )*newnp
                a2pbar = a2pbar + ( a2p/max(1,nAreap(ii)) )*newnp
                dpbar = dpbar + ( sqrtap/ max(1,nAreap(ii)) )*newnp     !   will become "<diameter>", but for now just sum sqrt(area)
                lpbar = lpbar + ( lengthp/max(1,nAreap(ii)) )*newnp
                np = np + newnp
                write(*,fmt='(i6,f16.8,i6,12f16.4,i6,f16.4,i12,i6)')        ii,isolevel,nPaths,amax,lmax,cmax,c2max,area,length,curv,curv2,areap,lengthp,curvp,curv2p,nAreap(ii),length2onareabar,count(img_blur(0:Nx-1,0:Ny-1)>=isolevel),newnp                                                                              
                write(unit=901,fmt='(i6,f16.8,i6,12f16.4,i6,f16.4,i12,i6)') ii,isolevel,nPaths,amax,lmax,cmax,c2max,area,length,curv,curv2,areap,lengthp,curvp,curv2p,nAreap(ii),length2onareabar,count(img_blur(0:Nx-1,0:Ny-1)>=isolevel),newnp
                oldnp = nAreap(ii)
                
                
            end do
            
            a2pbar = a2pbar / max(np,1)
            apbar = apbar / max(np,1)
            dpbar = ( 2 *dpbar / sqrt(PI) ) / max(np,1)         !   note - need to convert sum sqrt(area) into diameter : D = 2 sqrt(A/PI) 
            lpbar = lpbar / max(np,1)
            
            write(*,fmt='(a,i6)')    "  +ve area enclosures ",np
            write(*,fmt='(a,f16.4,10f16.8)') "  <a^2>,<a>,<d>,<l>         ",a2pbar,apbar,dpbar,lpbar
            write(unit=901,fmt='(a,i6)')    "#  +ve area enclosures ",np
            write(unit=901,fmt='(a,f16.4,10f16.8)') "# <a^2>,<a>,<d>,<l>         ",a2pbar,apbar,dpbar,lpbar
            
            
        !---    construct a histogram of the spot sizes ( diameters ) 
            if (nLevels>1) then
               
                if (pSpot_nBins == LIB_CLA_NODEFAULT_I) pSpot_nBins = ceiling( snapTo(sqrt( np*1.0d0 ) ) )
                aa = maxval( pSpotArea )
                if (pSpot_max == LIB_CLA_NODEFAULT_R) pSpot_max = snapTo( sqrt(4*aa/PI) )
                print *,"contourAnalyse.exe info - constructing frequency-diameter histogram with ",pSpot_nBins," to ",pSpot_max
                allocate(pSpot_hist(0:pSpot_nBins-1)) ; pSpot_hist = 0.0d0
                allocate(nSpot_hist(2,0:pSpot_nBins-1)) ; nSpot_hist = 0      !   nSpots_hist(1,k) = number of enclosures with diameter in kth bin. nSpots_hist(2,k) = number of _new spots_ with diameter in kth bin. 
                oldnp = 0
                do ii = nLevels,1,-1
                    !newnp = max(0,nSpots(ii)-oldnp)         !   number of spots which appear at this intensity 
                    newnp = max(0,nAreap(ii) - nEnclosure(ii-1))

                    do jj = 1,nSpots(ii)
                        aa = pSpotArea(jj,ii)
                        cc = sqrt( 4*aa/PI )
                        kk = int(  pSpot_nBins * cc / pSpot_max )
                        if (kk>=pSpot_nBins) cycle
                        nSpot_hist(1,kk) = nSpot_hist(1,kk) + 1
                        nSpot_hist(2,kk) = nSpot_hist(2,kk) + real(newnp)/nSpots(ii)
                        pSpot_hist(kk) = pSpot_hist(kk) + cc
                    end do
                    oldnp = nSpots(ii)
                    print *,"new ",ii,nSpots(ii),nAreap(ii),nEnclosure(ii),newnp
                end do
            !---    normalise histogram, find the new spots and the averages
                print *,""
                write (unit=901,fmt='(a)') ""
                write (*,fmt='(a6,6a16)') "bin","min","max","avg","n contour","spots"
                write (unit=901,fmt='(a6,6a16)') "# bin","min","max","avg","n contour","spots"
                do kk = 0,pSpot_nBins-1
                    pSpot_hist(kk) = pSpot_hist(kk) / max(1.0d0,nSpot_hist(1,kk))  
                    aa = kk*pSpot_max/pSpot_nBins
                    cc = (kk+1)*pSpot_max/pSpot_nBins
                    write (*,fmt='(i6,3f16.6,i16,f16.6)') kk,aa,cc,pSpot_hist(kk),nint(nSpot_hist(1,kk)),nSpot_hist(2,kk)
                    write (unit=901,fmt='(i6,3f16.6,i16,f16.6)') kk,aa,cc,pSpot_hist(kk),nint(nSpot_hist(1,kk)),nSpot_hist(2,kk)
                end do
            !---
            
            
            
            end if                
            
            
            
            
            write(unit=901,fmt='(a)') ""
            close(unit=901)
            if (op_contour .or. op_contourPositive) call writePng(trim(outfile_prefix)//"_contour.png",img_iso(0:nx-1,0:ny-1),normalise = normop, negative = negative)
            if (op_contour .or. op_contourPositive) deallocate(img_iso)        
            if (flattenLinear) call writePng(trim(outfile_prefix)//"_flatline.png",img_flatline,normalise = normop)
            if (flattenLinear) deallocate(img_flatline)   
            call pause( tt(T_ANALYSE) )
        end if
        
       
        
        
        
    !---    bye bye
    
        deallocate(img)
        deallocate(img_blur)
        
        if (.not. quiet) then
            print *,""
            print *,"time taken (read)      : ",elapsed(tt(T_READ))
            print *,"time taken (blur)      : ",elapsed(tt(T_BLUR))
            if (nLevels+nIso_spec>0) &
            print *,"time taken (analyse)   : ",elapsed(tt(T_ANALYSE))
            if (FTfilter) &
            print *,"time taken (FT)        : ",elapsed(tt(T_FT))
            if (nCorrelationFunctionBins /= LIB_CLA_NODEFAULT_I) &
            print *,"time taken (corr fn)   : ",elapsed(tt(T_CORR))
            print *,"time taken (total)     : ",elapsed(tt(T_TOTAL))
            print *,""
            print *,"done"
            print *,""
        end if
        
     
     contains
 !---^^^^^^^^
 
        pure function snapTo( x ) result( y )
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      snap value to 1,2,5,10,20,50... so that y>=x
            real(kind=real64),intent(in)            ::      x
            real(kind=real64)                       ::      y
            integer         ::      ii
            y = 1d-6
            do ii = 1,100        !   don't loop forever
                if (y>=x) return
                if (2*y>=x) then
                    y = 2*y ; return
                end if
                if (5*y>=x) then
                    y = 5*y ; return
                end if
                y = 10*y                
            end do
            
            return
        end function snapTo
    
    
        subroutine flattenInsidePath( img ,pathLen,path, isolevel,imgThresh )
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      I have a path defined in the image
    !*      reduce the img level by (up to) imgThresh ( usually will be difference between isolevels )
    !*      within the path
            real(kind=real64),dimension(0:,0:),intent(inout)    ::      img
            integer,intent(in)                                  ::      pathlen
            real(kind=real64),dimension(:,0:),intent(in)        ::      path
            real(kind=real64),intent(in)                        ::      isolevel,imgThresh
            
            logical     ::      done
            integer     ::      ixmin,ixmax,iymin,iymax,ix,iy
            integer     ::      jx,jy
            integer     ::      kk
            real(kind=real64)   ::      imgsum4,ff
            logical,dimension(:,:),allocatable  ::      mask
            
        !---    Now mask every point inside the path
                
            ixmin = int( minval(path(1,0:pathLen-1)) )
            ixmax = int( maxval(path(1,0:pathLen-1)) )
            iymin = int( minval(path(2,0:pathLen-1)) )
            iymax = int( maxval(path(2,0:pathLen-1)) )
            allocate(mask(ixmin:ixmax,iymin:iymax))
            mask = .false.          !   mask set to true inside contour
            
        !   set mask true for each point on path. This should draw a ring, importantly it is a good start for the masked region.    
            do kk = 0,pathLen-1     
                ix = int( path(1,kk) )
                iy = int( path(2,kk) )
                mask(ix,iy) = .true.
            end do
        
        !   flood fill mask inside path                
            do 
                done = .true.                               !   iterate until complete
                do iy = iymin,iymax                         !   only search range covered by path
                    do ix = ixmin,ixmax
                        if ( .not. mask(ix,iy)) cycle       !   look for a lit pixel
                        
                        do jy = max(iymin,iy-1),min(iymax,iy+1)         !   check neighbours of lit pixel
                            do jx = max(ixmin,ix-1),min(ixmax,ix+1)   
                                if (mask(jx,jy)) cycle                  !   already done.
                                
                                if (img(jx,jy)>=isolevel) then
                                    !   possibly need to set this point inside
                                    if ((jx-ix)*(jx-ix) + (jy-iy)*(jy-iy) == 1) then
                                        !   yes, this is a nearest neighbour site, and so must be masked.
                                        mask(jx,jy) = .true.
                                        done = .false.
                                    else
                                        !   this is a next nearest neighbour site. Check for topology.
                                        imgsum4 = img(ix,iy) + img(jx,iy) + img(ix,jy) + img(jx,jy) 
                                        if (imgsum4 >= 4*isolevel) then
                                            !   the interpolation at the corner between sites is lit. Mask
                                            mask(jx,jy) = .true.
                                            done = .false.
                                        end if
                                    end if
                                end if
                                
                            end do
                        end do
                    end do
                end do
                if (done) exit
            end do
            
        !---    now for each masked point, drop the intensity level by up to  imgThresh
            do iy = iymin,iymax                         !   only search range covered by path
                do ix = ixmin,ixmax
                    if ( .not. mask(ix,iy)) cycle       !   look for a lit pixel
                    ff = img(ix,iy)
                    img(ix,iy) = ff - min(imgThresh,ff-isolevel)
                end do
            end do        
            
            
            
            return
        end subroutine flattenInsidePath
        
        
        
        subroutine findCorrelationFunction( img, pbc, rMax,nBins, c )
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      compute the correlation function 
    !*          c(r) = int_x int_y | img(x)-img(y) | delta( |x-y|-r ) d2x d2y / ( 2 pi r )
    !*      
            real(kind=real64),dimension(0:,0:),intent(in)       ::      img
            logical,intent(in)                                  ::      pbc
            real(kind=real64),intent(in)                        ::      rMax
            integer,intent(in)                                  ::      nBins
            real(kind=real64),dimension(0:),intent(out)         ::      c           !   (0:nBins) with c(0) = 0:delta , c(1) = delta:2 delta etc, delta = rMax/nBins
            
            integer             ::      Nx,Ny
            integer             ::      ix,iy,jx,jy,kx,ky 
            real(kind=real64),dimension(-1:nBins+2)          ::      c_tmp
            integer             ::      nrmax,rr,ii,nrmax2 
            integer,dimension(:),allocatable                 ::      bin  !, ninbin         
            real(kind=real64),dimension(:,:),allocatable     ::      weight 
            real(kind=real64),dimension(:),allocatable       ::      winbin   
            real(kind=real64)   ::      ww,delta,idelta,xx,yy,pm,p0,pp
            real(kind=real64)   ::      fi,fj,absfij,absfijbar,fbar
            
            Nx = size(img,dim=1)
            Ny = size(img,dim=2)
            nrmax = ceiling(rMax)
            delta = rMax/nBins
            idelta = 1/delta
            
            print *,"contourAnalyse::findCorrelationFunction() info - nBins ",nBins," rMax ",rMax," delta ",delta
            
            
            
            nrmax2 =  (nrmax+1)*(nrmax+1)                           !   max range searched for (ix^2 + iy^2)
            allocate(weight(-1:1,0:nrmax2)) ; weight = 0            !   a pixel at (ix^2 + iy^2) adds to three bins, with Gaussian weighting
            allocate(bin(0:2*nrmax*nrmax)) ; bin = -1               !   which bin is at centre for pixel at (ix^2 + iy^2)
!            allocate(ninbin(0:nrmax2)) ; ninbin = 0
            
            
        !---    first find weight on bins  
            !ninbin = 0             
            do ix = 0,nrmax
                do iy = 0,nrmax
                
                    rr = ix*ix + iy*iy  
                    if (rr > nrmax2) cycle               !   out of range
                    ii = bin(rr)
                    if (ii /= -1) then                      !   have done this bin. record we've seen it again                
!                         if (ix+iy == 0) then
!                             ninbin(rr) = ninbin(rr) + 1
!                         else if (ix*iy == 0) then
!                             ninbin(rr) = ninbin(rr) + 2   
!                         else
!                             ninbin(rr) = ninbin(rr) + 4
!                         end if
                        cycle
                    end if
                          
                    
                    yy = sqrt( real(rr,kind=real64) )       !   distance as real number
                    ii = nint( yy*idelta )                   !   bin containing this distance
                                        
                    xx = ii*delta - yy                                      
                    p0 =      exp( -2*xx*xx )
                    pm = p0 * exp( -2*delta*(-2*xx + delta) )
                    pp = p0 * exp( -2*delta*(+2*xx + delta) )
                    
                    if (ii==0) pm = 0.0d0
                   
                  
                    ww = 1/(pm + p0 + pp)
                    pm = pm * ww
                    p0 = p0 * ww
                    pp = pp * ww
                    
                    bin(rr) = ii    
                    weight(-1:1,rr) = (/pm,p0,pp/)
                    
!                     if (ix+iy == 0) then
!                         ninbin(rr) = ninbin(rr) + 1
!                     else if (ix*iy == 0) then
!                         ninbin(rr) = ninbin(rr) + 2
!                     else
!                         ninbin(rr) = ninbin(rr) + 4
!                     end if
                 
                end do
            end do


        !---    now have which bin each radius squared hits, and its weight.
            allocate(winbin(-1:nBins+2))
!             winbin = 0
!             do rr = 0,nrmax2
!                 ii = bin(rr)
!                 if (ii /= -1) then                
!                     winbin(ii-1) = winbin(ii-1) + weight(-1,rr)*ninbin(rr)
!                     winbin(ii  ) = winbin(ii  ) + weight( 0,rr)*ninbin(rr)
!                     winbin(ii+1) = winbin(ii+1) + weight(+1,rr)*ninbin(rr)                
!                 end if
!             end do                
!             
!             winbin(-1) = huge(1.0)
!             do rr = 0,nrmax2
!                 ii = bin(rr)
!                 if (ii /= -1) then                
!                     weight(-1,rr) = weight(-1,rr) / winbin(ii-1)
!                     weight( 0,rr) = weight( 0,rr) / winbin(ii  )
!                     weight(+1,rr) = weight(+1,rr) / winbin(ii+1)
!                     
!                     
!                     if (rr<100) print *,rr,ii,winbin(ii),weight(:,rr)
!                     
!                 end if              
!                 
!                
!                   
!             end do                
            
!             deallocate( winbin )
!             deallocate( ninbin )    
            
        !---    now can compute correlation function
            absfijbar = 0.0d0 
            fbar = 0.0d0
            do iy = 0,Ny-1
                do ix = 0,Nx-1
                    fi = img( ix,iy )
                    fbar = fbar + fi
                    do jy = 0,Ny-1
                        do jx = 0,Nx-1
                            fj = img( jx,jy )        
                            absfij = abs( fi - fj ) **2
                            absfijbar = absfijbar + absfij
                        end do
                    end do
                end do
            end do
            ww = 1.0d0/(Nx*Ny)
            fbar = fbar * ww
            absfijbar = absfijbar *ww*ww
            print *,"<fi>      ", fbar     ," in ",Nx*Ny," px"  
            print *,"<|fi-fj|> ", absfijbar   
            
            
                    
            winbin = 0
            c_tmp = 0
!             absfijbar = 0
!             nabsfij = 0
!             fbar = sum(img)/(Nx*Ny)
            do iy = 0,Ny-1
                do ix = 0,Nx-1
                    fi = img( ix,iy )
!                     if (pbc) then
!                         do ky = 0,nrmax
!                             jy = mod( iy + Ny + ky , Ny )
!                             do kx = -nrmax,nrmax
!                             
!                                 rr = kx*kx + ky*ky
!                                 ii = bin(rr)
!                                 if (ii==-1) cycle
!                                                         
!                                 jx = mod( ix + Nx + kx , Nx )                               
!                                 fj = img( jx,jy )
!                                 absfij = ( fi - fj )*( fi - fj )
!                                 
!                                 if (ky > 0) then 
!                                     c_tmp(ii-1:ii+1) = c_tmp(ii-1:ii+1) + 2*weight(-1:1,rr)*absfij
!                                     winbin( ii-1:ii+1 ) = winbin( ii-1:ii+1 ) + 2*weight(-1:1,rr)
!                                     absfijbar = absfijbar + 2*absfij
!                                     nabsfij = nabsfij + 2                                    
!                                 else
!                                     c_tmp(ii-1:ii+1) = c_tmp(ii-1:ii+1) + weight(-1:1,rr)*absfij
!                                     winbin( ii-1:ii+1 ) = winbin( ii-1:ii+1 ) + weight(-1:1,rr)
!                                     absfijbar = absfijbar + absfij
!                                     nabsfij = nabsfij + 1                                   
!                                 end if
!                             end do
!                         end do
!                    else
                        do ky = -nrmax,nrmax
                            jy = iy+ky
                            do kx = -nrmax,nrmax
                            
                                rr = kx*kx + ky*ky
                                ii = bin(rr)
                                if (ii==-1) cycle
                                                        
                                jx = ix+kx
                                
                                
                                fj = fbar
                                if (pbc) then
                                    fj = img( mod(jx+Nx,Nx),mod(jy+Ny,Ny) )
                                else if ( (jx*(Nx-1-jx)>=0).and.(jy*(Ny-1-jy)>=0) ) then
                                    fj = img( jx,jy )
                                end if
                                
                                absfij = abs( fi - fj ) **2                        
                                c_tmp(ii-1:ii+1) = c_tmp(ii-1:ii+1) + weight(-1:1,rr)*absfij
                                winbin( ii-1:ii+1 ) = winbin( ii-1:ii+1 ) + weight(-1:1,rr)

                            end do
                        end do
!                     else
!                         do jy = iy,min(iy+nrmax,Ny-1)
!                             ky = jy - iy
!                             do jx = max(0,ix-nrmax),min(ix+nrmax,Nx-1)
!                                 kx = jx - ix
!                                 
!                                 rr = kx*kx + ky*ky
!                                 ii = bin(rr)
!                                 if (ii==-1) cycle
!                                                         
!                                 fj = img( jx,jy )
!                                 absfij = abs( fi - fj )
!                                 
!                                 
!                                 if (ky > 0) then 
!                                     c_tmp(ii-1:ii+1) = c_tmp(ii-1:ii+1) + 2*weight(-1:1,rr)*absfij
!                                     winbin( ii-1:ii+1 ) = winbin( ii-1:ii+1 ) + 2*weight(-1:1,rr)
!                                     absfijbar = absfijbar + 2*absfij
!                                     nabsfij = nabsfij + 2                                    
!                                 else
!                                     c_tmp(ii-1:ii+1) = c_tmp(ii-1:ii+1) + weight(-1:1,rr)*absfij
!                                     winbin( ii-1:ii+1 ) = winbin( ii-1:ii+1 ) + weight(-1:1,rr)
!                                     absfijbar = absfijbar + absfij
!                                     nabsfij = nabsfij + 1                                   
!                                 end if
!                             end do                                                                             
!                         end do
!                    end if
                end do
            end do   
                           
        !---    now normalise the bins
            !absfijbar = absfijbar/nabsfij

            do ii = 0,nBins
                print *,ii,ii*delta,winbin(ii),c_tmp(ii),c_tmp(ii)/max(1.0d-8,absfijbar * winbin(ii)) 
            end do
            
            
            
            do ii = 0,nBins
                ww =  max(1.0d-8,absfijbar * winbin(ii)) 
                c(ii) = c_tmp(ii)/ww
            end do                    
          !  c(0:nBins) = c_tmp(0:nBins)
          !  ww = 1/c(nBins)
          !  c = c * ww
            
            
            return
        end subroutine findCorrelationFunction
            
!             
!             
!         subroutine findCorrelationFunction( img, pbc, rMax,nBins, c )
!     !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
!     !*      compute the correlation function 
!     !*          c(r) = int_x int_y | img(x)-img(y) | delta( |x-y|-r ) d2x d2y / ( 2 pi r )
!     !*      
!             real(kind=real64),dimension(0:,0:),intent(in)       ::      img
!             logical,intent(in)                                  ::      pbc
!             real(kind=real64),intent(in)                        ::      rMax
!             integer,intent(in)                                  ::      nBins
!             real(kind=real64),dimension(0:),intent(out)         ::      c           !   (0:nBins) with c(0) = 0:delta , c(1) = delta:2 delta etc, delta = rMax/nBins
!             
!             integer             ::      Nx,Ny
!             integer             ::      ix,iy,jx,jy,kx,ky 
!             real(kind=real64),dimension(-1:nBins+2)          ::      c_tmp
!             integer             ::      nrmax,rr,ii,nrmax2
!             integer,dimension(:),allocatable                 ::      bin  , ninbin         
!             real(kind=real64),dimension(:,:),allocatable     ::      weight 
!             real(kind=real64),dimension(:),allocatable       ::      winbin   
!             real(kind=real64)   ::      ww,delta,idelta,xx,yy,pm,p0,pp
!             real(kind=real64)   ::      fi,fj,absfij
!             
!             Nx = size(img,dim=1)
!             Ny = size(img,dim=2)
!             nrmax = ceiling(rMax)
!             delta = rMax/nBins
!             idelta = 1/delta
!             
!             print *,"contourAnalyse::findCorrelationFunction() info - nBins ",nBins," rMax ",rMax," delta ",delta
!             
!             
!             
!             nrmax2 =  (nrmax+1)*(nrmax+1)
!             allocate(weight(-1:1,0:nrmax2)) ; weight = 0
!             allocate(bin(0:2*nrmax*nrmax)) ; bin = -1            
!             allocate(ninbin(0:nrmax2)) ; ninbin = 0
!             
!             
!         !---    first find weight on bins  
!             ninbin = 0             
!             do ix = 0,nrmax
!                 do iy = 0,nrmax
!                 
!                     rr = ix*ix + iy*iy  
!                     if (rr > nrmax2) cycle               !   out of range
!                     ii = bin(rr)
!                     if (ii /= -1) then                      !   have done this bin. record we've seen it again                
!                         if (ix+iy == 0) then
!                             ninbin(rr) = ninbin(rr) + 1
!                         else if (ix*iy == 0) then
!                             ninbin(rr) = ninbin(rr) + 2   
!                         else
!                             ninbin(rr) = ninbin(rr) + 4
!                         end if
!                         cycle
!                     end if
!                           
!                     
!                     yy = sqrt( real(rr,kind=real64) )       !   distance as real number
!                     ii = nint( yy*idelta )                   !   bin containing this distance
!                                         
!                     xx = ii*delta - yy                                      
!                     p0 =      exp( -2*xx*xx )
!                     pm = p0 * exp( -2*delta*(-2*xx + delta) )
!                     pp = p0 * exp( -2*delta*(+2*xx + delta) )
!                     
!                     if (ii==0) pm = 0.0d0
!                    
!                   
!                     ww = 1/(pm + p0 + pp)
!                     pm = pm * ww
!                     p0 = p0 * ww
!                     pp = pp * ww
!                     
!                     bin(rr) = ii    
!                     weight(-1:1,rr) = (/pm,p0,pp/)
!                     
!                     if (ix+iy == 0) then
!                         ninbin(rr) = ninbin(rr) + 1
!                     else if (ix*iy == 0) then
!                         ninbin(rr) = ninbin(rr) + 2
!                     else
!                         ninbin(rr) = ninbin(rr) + 4
!                     end if
!                  
!                 end do
!             end do
! 
! 
!         !---    now have a count in each bin, and the weight for each bin and +/- 1.
!             allocate(winbin(-1:nBins+2))
!             winbin = 0
!             do rr = 0,nrmax2
!                 ii = bin(rr)
!                 if (ii /= -1) then                
!                     winbin(ii-1) = winbin(ii-1) + weight(-1,rr)*ninbin(rr)
!                     winbin(ii  ) = winbin(ii  ) + weight( 0,rr)*ninbin(rr)
!                     winbin(ii+1) = winbin(ii+1) + weight(+1,rr)*ninbin(rr)                
!                 end if
!             end do                
!             
!             winbin(-1) = huge(1.0)
!             do rr = 0,nrmax2
!                 ii = bin(rr)
!                 if (ii /= -1) then                
!                     weight(-1,rr) = weight(-1,rr) / winbin(ii-1)
!                     weight( 0,rr) = weight( 0,rr) / winbin(ii  )
!                     weight(+1,rr) = weight(+1,rr) / winbin(ii+1)
!                     
!                     
!                     if (rr<100) print *,rr,ii,winbin(ii),weight(:,rr)
!                     
!                 end if              
!                 
!                
!                   
!             end do                
!             
!             deallocate( winbin )
!             deallocate( ninbin )    
!             
!         !---    now can compute correlation function
!         
!             
!             c_tmp = 0
!             do iy = 0,Ny-1
!                 do ix = 0,Nx-1
!                     fi = img( ix,iy )
!                     if (pbc) then
!                         do ky = 0,nrmax
!                             jy = mod( iy + Ny + ky , Ny )
!                             do kx = -nrmax,nrmax
!                             
!                                 rr = kx*kx + ky*ky
!                                 ii = bin(rr)
!                                 if (ii==-1) cycle
!                                                         
!                                 jx = mod( ix + Nx + kx , Nx )                               
!                                 fj = img( jx,jy )
!                                 absfij = ( fi - fj )*( fi - fj )
!                                 
!                                 if (ky > 0) absfij = 2*absfij
!                                 c_tmp(ii-1:ii+1) = c_tmp(ii-1:ii+1) + weight(-1:1,rr)*absfij
!                             end do
!                         end do
!                     else
!                         do jy = iy,min(iy+nrmax,Ny-1)
!                             ky = jy - iy
!                             do jx = max(0,ix-nrmax),min(ix+nrmax,Nx-1)
!                                 kx = jx - ix
!                                 
!                                 rr = kx*kx + ky*ky
!                                 ii = bin(rr)
!                                 if (ii==-1) cycle
!                                                         
!                                 fj = img( jx,jy )
!                                 absfij = abs( fi - fj )
!                                 if (ky > 0) absfij = 2*absfij
!                                 c_tmp(ii-1:ii+1) = c_tmp(ii-1:ii+1) + weight(-1:1,rr)*absfij
!                             end do                                                                             
!                         end do
!                     end if
!                 end do
!             end do   
!                                 
!             c(0:nBins) = c_tmp(0:nBins)
!             ww = 1/c(nBins)
!             c = c * ww
!           !  do ii = 0,nBins
!           !      print *,ii,ii*delta,c(ii)
!           !  end do
!             
!             
!             return
!         end subroutine findCorrelationFunction
!             
            
        
        subroutine contiguousPeaks( img_in,iso,indx,peakHeight, nEnclosure )
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      starting with isolevel 1, define contiguous peaks 
    !*      then at isolevel 2, add pixels to the lowest number contiguous peak,
    !*      and so on, until each pixel has a peak number.
    !*      make peaks for regions _higher_ than isolevel
    !*
    !*      return the number of peaks, the height of each peak
    !*      and for each isolevel the number of spots which appear
    
    !*      iso 1
    !*               __  
    !*              /  \
    !*              | 1|
    !*              \_ /
    !*
    !*
    !*      ----------------------------------------------
    !*      iso 2    __  
    !*              /  \ 
    !*              | 1| 
    !*              \_ /          ___  
    !*                           /   \ 
    !*                           | 2 | 
    !*                           \__ / 
    !*
    !*      ----------------------------------------------
    !*      iso 3___________________
    !*          /    __             \                      
    !*          |   /  \             \        ___               
    !*          |   | 1|              \      /   \              
    !*          |   \_ /          ___  \     | 3 |              
    !*          |                /   \  \    \__ /              
    !*          |          1     | 2 |  |                       
    !*           \               \__ /  |                 
    !*            \                    /                   
    !*             \__________________/
    !*
    !*
            real(kind=real64),dimension(:,:),intent(in)             ::      img_in
            real(kind=real64),dimension(:),intent(in)               ::      iso 
            integer,dimension(:,:),intent(out)                      ::      indx
            real(kind=real64),dimension(:),allocatable,intent(out)  ::      peakHeight
            
            integer,dimension(:),allocatable,intent(out)            ::      nEnclosure
            
            
            integer,dimension(:),allocatable                        ::      nRegions,nNegative,nSpots
            
            
            integer             ::      Nx,Ny,Niso
            
            integer             ::      ix,iy,kk,ki,kj
            integer             ::      nPeaks
            integer,dimension(:,:),allocatable                      ::      img 
            logical             ::      ok
            
            
            integer,dimension(:),allocatable                        ::      reindx,rereindx
            
            real(kind=real64),dimension(:,:),allocatable                      ::      img_tmp
            integer,dimension(:,:),allocatable                      ::      indx_tmp
            
            Nx = size(img_in,dim=1)
            Ny = size(img_in,dim=2)
            Niso = size(iso)
        
        !---    convert real map into integer map, with each pixel being the highest isolevel reached.    
            allocate(img(0:Nx+1,0:Ny+1))
            allocate(img_tmp(0:Nx+1,0:Ny+1))
            img_tmp(0,:)        = -huge(1.0)
            img_tmp(Nx+1,:)     = -huge(1.0)
            img_tmp(:,0)        = -huge(1.0)
            img_tmp(0,Ny+1)     = -huge(1.0)
            img_tmp(1:Nx,1:Ny)  = img_in(1:Nx,1:Ny)
            allocate(indx_tmp(0:Nx+1,0:Ny+1))
            
            
            allocate(nRegions(nIso))
            allocate(nNegative(nIso))
            allocate(nEnclosure(0:nIso))
            allocate(nSpots(nIso))
             
            
            
            
            !print *,"min max ",minval(img_in),maxval(img_in)

            img = 0            
            do iy = 1,Ny                    
                do ix = 1,Nx
                    do kk = 1,Niso
                        if (img_in(ix,iy)>=iso(kk)) img(ix,iy) = kk
                    end do
                end do
            end do
            
            !do kk = 1,Niso
            !    print *,"level ",kk," iso ",iso(kk)," count ",count(img==kk),count(img_in>=iso(kk))
            !end do
            
            
            
            
           ! print *,"px 293,576",img_in(293,576),img(293,576),indx_tmp(293,576)
           ! print *,img(292:294,575)
           ! print *,img(292:294,576)
           ! print *,img(292:294,577)
           ! print *,indx_tmp(292:294,575)
           ! print *,indx_tmp(292:294,576)
           ! print *,indx_tmp(292:294,577)
            
                     
        !---    now add the peaks
            nPeaks = 0
            indx_tmp = 0
            do kk = Niso,1,-1
            
                !   at isolevel k, number all pixels at this level.
                do iy = 1,Ny                    
                    do ix = 1,Nx                     
                        if (img(ix,iy)==kk) then
                            nPeaks = nPeaks + 1
                            indx_tmp(ix,iy) = nPeaks
                        end if
                    end do
                end do
           !     print *,"level ",kk," iso ",iso(kk)," pixels currently ",nPeaks
                     
                !   now make all pixels at this level contiguous regions
                !   do this in two passes - LRUD then corners
                do
                    ok = .true.    
                    do iy = 1,Ny                    
                        do ix = 1,Nx          
                             
                            if (img(ix,iy) == kk) then
                                !   this pixel is at the right level.
                                ki = indx_tmp(ix,iy)
                                
                                
                            !---    contiguity LRUD
                                if (img(ix-1,iy)==kk) then      !   left pixel at same isolevel
                                    kj = indx_tmp(ix-1,iy)      !   current indexing of left pixel ...
                                    if (ki /= kj) then          !   ... is different to indexing of pixel (ix,iy)
                                        ki = min( ki,kj )       !   which has the lower index? I'll use that.
                                        indx_tmp(ix-1,iy) = ki  !   set left pixel ...
                                        indx_tmp(ix,iy) = ki    !   ... and central pixel to same index
                                        ok = .false.            !   and we'll have to make another pass.
                                    end if
                                end if
                                   
                                if (img(ix+1,iy)==kk) then      !   right pixel etc.
                                    kj = indx_tmp(ix+1,iy)
                                    if (ki /= kj) then
                                        ki = min( ki,kj )
                                        indx_tmp(ix+1,iy) = ki
                                        indx_tmp(ix,iy) = ki
                                        ok = .false.
                                    end if
                                end if
                                
                                if (img(ix,iy-1)==kk) then
                                    kj = indx_tmp(ix,iy-1)
                                    if (ki /= kj) then
                                        ki = min( ki,kj )
                                        indx_tmp(ix,iy-1) = ki
                                        indx_tmp(ix,iy) = ki
                                        ok = .false.
                                    end if
                                end if
                                
                                if (img(ix,iy+1)==kk) then
                                    kj = indx_tmp(ix,iy+1)
                                    if (ki /= kj) then
                                        ki = min( ki,kj )
                                        indx_tmp(ix,iy+1) = ki
                                        indx_tmp(ix,iy) = ki
                                        ok = .false.
                                    end if
                                end if
                                
                            end if  !   pixel(ix,iy) is at right level
                                
      !                end do
      !            end do           !   scan over ix,iy
      !           
      !            if (ok) exit     !   no more changes to contiguity this scan
      !        end do         
      !                        
      !        do
      !            ok = .true.    
      !            do iy = 1,Ny                    
      !                do ix = 1,Nx          
      !                     
                            if (img(ix,iy) == kk) then
                                !   this pixel is at the right level.
                                ki = indx_tmp(ix,iy)
                                
                                
                            !---    contiguity corners: note we need to check the interpolated level at the corner of the pixel 
                            !       in order to ensure topological correctness.
                                if (img(ix-1,iy-1)==kk) then                        !   top left pixel same isolevel
                                    kj = indx_tmp(ix-1,iy-1)
                                    if (ki /= kj) then                              !   ... but different index
                                        if ( sum( img_tmp(ix-1:ix,iy-1:iy) ) >= 4*iso(kk) ) then    !   linear interpolation of midpoint between ix,iy and ix-1,iy-1 > isolevel.
                                            ki = min( ki,kj )
                                            indx_tmp(ix-1,iy-1) = ki
                                            indx_tmp(ix,iy) = ki
                                            ok = .false.
                                        end if
                                    end if
                                end if
                                
                                if (img(ix+1,iy-1)==kk) then
                                    kj = indx_tmp(ix+1,iy-1)
                                    if (ki /= kj) then
                                        if ( sum( img_tmp(ix:ix+1,iy-1:iy) ) >= 4*iso(kk) ) then
                                            ki = min( ki,kj )
                                            indx_tmp(ix+1,iy-1) = ki
                                            indx_tmp(ix,iy) = ki
                                            ok = .false.
                                        end if
                                    end if
                                end if
                                 
                                if (img(ix-1,iy+1)==kk) then
                                    kj = indx_tmp(ix-1,iy+1)
                                    if (ki /= kj) then
                                        if ( sum( img_tmp(ix-1:ix,iy:iy+1) ) >= 4*iso(kk) ) then
                                            ki = min( ki,kj )
                                            indx_tmp(ix-1,iy+1) = ki
                                            indx_tmp(ix,iy) = ki
                                            ok = .false.
                                        end if
                                    end if
                                end if
                                 
                                if (img(ix+1,iy+1)==kk) then
                                    kj = indx_tmp(ix+1,iy+1)
                                    if (ki /= kj) then
                                        if ( sum( img_tmp(ix:ix+1,iy:iy+1) ) >= 4*iso(kk) ) then
                                            ki = min( ki,kj )
                                            indx_tmp(ix+1,iy+1) = ki
                                            indx_tmp(ix,iy) = ki
                                            ok = .false.
                                        end if
                                    end if
                                end if
                                     
                            end if  !   pixel(ix,iy) is at right level
                                
                        end do
                    end do           !   scan over ix,iy
                   
                    if (ok) exit     !   no more changes to contiguity this scan
                end do         
            end do
            
           ! print *,"px 293,576",img_in(293,576),img(293,576),indx_tmp(293,576)
           ! print *,img(292:294,575)
           ! print *,img(292:294,576)
           ! print *,img(292:294,577)
           ! print *,indx_tmp(292:294,575)
           ! print *,indx_tmp(292:294,576)
           ! print *,indx_tmp(292:294,577)
           ! 
           ! 
           ! call writePng("test1_img.png",real(img,kind=real64),normalise=.true.)
            
            
                
        !---    at this point I have topologically correctly indexed regions of similar intensity level.
        !       But not peaks.
        !       I do, however, know that a region at iso(k) has a lower region number than a region at iso(k-1)
        !       ( because previous scan was from nIso to 1 )
        !       so I scan from bottom to top 1 to nIso, looking for the lowest region number touched at a higher iso level
        
           ! call writePng("test2_indx_tmp.png",real(indx_tmp,kind=real64),normalise=.true.)
            
            
        !---    count the number of regions at each isolevel            
        
!            print *,"minmax 1:Nx,1:Ny img_in   ",minval(img_in(1:Nx,1:Ny))  ,maxval(img_in(1:Nx,1:Ny))  
!            print *,"minmax 1:Nx,1:Ny img      ",minval(img(1:Nx,1:Ny))     ,maxval(img(1:Nx,1:Ny))     
!            print *,"minmax 1:Nx,1:Ny indx_tmp ",minval(indx_tmp(1:Nx,1:Ny)),maxval(indx_tmp(1:Nx,1:Ny))
!        
        
            allocate(reindx(0:nPeaks))     
            reindx = 0
            do iy = 1,Ny                    
                do ix = 1,Nx         
                    ki = indx_tmp(ix,iy)
                    reindx(ki) = img(ix,iy)     !   record the isolevels                 
                end do
            end do
            
            nRegions = 0             
            do ix = 1,nPeaks
                kk = reindx(ix)
                if (kk>0) nRegions(kk) = nRegions(kk) + 1
            end do
            do kk = 1,nIso
                if (nRegions(kk) == 0) then
                    nRegions(kk) = 1
                else
                    exit
                end if
            end do
            
            
            
        !---    now combine isolevel regions into peaks                           
            do ix = 1,nPeaks
                reindx(ix) = ix             !   it may be that a region does not touch any higher iso level aka its a maximum. So reindex to self.
            end do
            
            
           ! print *,"px 524,573 ",img_in(524,573),img(524,573),indx_tmp(524,573)
           ! print *,img(523:525,572)
           ! print *,img(523:525,573)
           ! print *,img(523:525,574)
           ! print *,indx_tmp(523:525,572)
           ! print *,indx_tmp(523:525,573)
           ! print *,indx_tmp(523:525,574)
            
        !   lowest index touching pixel LRUD 
            do iy = 1,Ny                    
                do ix = 1,Nx          
                         
                    ki = indx_tmp(ix,iy)                    !   index of pixel (ix,iy)
                    
                    kj = indx_tmp(ix-1,iy)                  !   index of pixel to left
                    if ((ki*kj>0).and.(kj<ki)) reindx(ki) = min(reindx(ki),reindx(kj))              !   note: can't touch a lower index unless pixel to left _also_ has a higher isolevel
                    
                         
                    kj = indx_tmp(ix+1,iy)
                    if ((ki*kj>0).and.(kj<ki)) reindx(ki) = min(reindx(ki),reindx(kj))
                    
                    kj = indx_tmp(ix,iy-1)
                    if ((ki*kj>0).and.(kj<ki)) reindx(ki) = min(reindx(ki),reindx(kj))
                    
                    kj = indx_tmp(ix,iy+1)
                    if ((ki*kj>0).and.(kj<ki)) reindx(ki) = min(reindx(ki),reindx(kj))
                        
!                    end do
!                end do                          
!                 
!            
!        !   lowest index touching pixel corner 
!            do iy = 1,Ny                    
!                do ix = 1,Nx          
!                     
!                    ki = indx_tmp(ix,iy)
                    
                    kj = indx_tmp(ix-1,iy-1)
                    if ((ki*kj>0).and.(kj<ki)) then           !   To be topologically correct, have to check the value of the corner of the pixel
                        if ( sum( img_tmp(ix-1:ix,iy-1:iy) ) >= 4*img_in(ix,iy) ) reindx(ki) = min(reindx(ki),reindx(kj))          !  value increases between ix,iy and ix-1,iy-1.              
                    end if
                    
                    kj = indx_tmp(ix+1,iy-1)
                    if ((ki*kj>0).and.(kj<ki)) then
                        kk = img(ix,iy)
                        if ( sum( img_tmp(ix:ix+1,iy-1:iy) ) >= 4*img_in(ix,iy) ) reindx(ki) = min(reindx(ki),reindx(kj))                         
                    end if
                    
                    kj = indx_tmp(ix-1,iy+1)
                    if ((ki*kj>0).and.(kj<ki)) then
                        kk = img(ix,iy)
                        if ( sum( img_tmp(ix-1:ix,iy:iy+1) ) >= 4*img_in(ix,iy) ) reindx(ki) = min(reindx(ki),reindx(kj))                            
                    end if
                    
                    kj = indx_tmp(ix+1,iy+1)
                    if ((ki*kj>0).and.(kj<ki)) then
                        kk = img(ix,iy)
                        if ( sum( img_tmp(ix:ix+1,iy:iy+1) ) >= 4*img_in(ix,iy) ) reindx(ki) = min(reindx(ki),reindx(kj))                           
                    end if
                end do
            end do                          
                        
           ! print *,"px 524,573 ",img_in(524,573),img(524,573),indx_tmp(524,573)
           ! print *,reindx(indx_tmp(523,572)),reindx(indx_tmp(524,572)),reindx(indx_tmp(525,572))
           ! print *,reindx(indx_tmp(523,573)),reindx(indx_tmp(524,573)),reindx(indx_tmp(525,573))
           ! print *,reindx(indx_tmp(523,574)),reindx(indx_tmp(524,574)),reindx(indx_tmp(525,574))
           
            
            
                        
        !   I now know that a pixel with indx_tmp(ix,iy) = ki touches a pixel with index reindx(ki) <= ki
            do iy = 1,Ny    
                do ix = 1,Nx
                    ki = indx_tmp(ix,iy)
                    indx_tmp(ix,iy) = reindx(ki)
                end do
            end do
           ! call writePng("test4_indx_tmp.png",real(indx_tmp,kind=real64),normalise=.true.)
           ! print *,"px 524,573 ",img_in(524,573),img(524,573),indx_tmp(524,573)
           ! print *,reindx(indx_tmp(523,572)),reindx(indx_tmp(524,572)),reindx(indx_tmp(525,572))
           ! print *,reindx(indx_tmp(523,573)),reindx(indx_tmp(524,573)),reindx(indx_tmp(525,573))
           ! print *,reindx(indx_tmp(523,574)),reindx(indx_tmp(524,574)),reindx(indx_tmp(525,574))
           
            
        !   now each peak is correctly marked, but the peaks are not yet ordered.
            allocate(rereindx(0:nPeaks))
            nPeaks = 0
            rereindx = 0
            do iy = 1,Ny    
                do ix = 1,Nx
                    ki = indx_tmp(ix,iy)
                    if (ki==0) cycle
                    if (.not. any(reindx(1:nPeaks)==ki)) then
                        nPeaks = nPeaks + 1
                        reindx(nPeaks) = ki
                        rereindx(ki) = nPeaks
                    end if
                end do
            end do
            do iy = 1,Ny    
                do ix = 1,Nx
                    ki = indx_tmp(ix,iy)
                    indx_tmp(ix,iy) = rereindx(ki)
                end do
            end do
            
            indx(1:Nx,1:Ny) = indx_tmp(1:Nx,1:Ny)
            !print *,"contourAnalyse::contiguousPeaks info - number of peaks ",nPeaks
            
            allocate(peakHeight(nPeaks))
            peakHeight = -huge(1.0)
            do iy = 1,Ny    
                do ix = 1,Nx
                    ki = indx(ix,iy)
                    if (ki==0) cycle
                    peakHeight(ki) = max( peakHeight(ki),img_in(ix,iy) )
                end do
            end do
            
!           print *,""
!           write(*,fmt='(a8,a12,a16)') "peak","area (px)","height"
!           do ki = 1,nPeaks
!               write(*,fmt='(i8,i12,f16.8)') ki,count(indx(:,:)==ki),peakHeight(ki)
!           end do
!           print *,""
            
            
           ! print *,""
           ! write(*,fmt='(a8,a32,a12)') "isolevel","range","peak count"
           ! write(*,fmt='(i8,a16,f16.8,i12)') 0,"low",iso(1),count( peakHeight(:)<iso(1) )
           ! do kk = 1,nIso-1
           !     write(*,fmt='(i8,2f16.8i12)') kk,iso(kk),iso(kk+1),count( (peakHeight(:)>=iso(kk)).and.(peakHeight(:)<iso(kk+1)) )
           ! end do
           ! write(*,fmt='(i8,f16.8,a16,i12)') nIso,iso(nIso),"high",count( (peakHeight(:)>=iso(nIso)) )
           ! print *,""
            
            
        !---    count the number of peaks which vanish between isolevels
            nNegative = 0
            do ki = 1,nPeaks
                if ( peakHeight(ki) < iso(1) ) then
                    !nNegative(0) = nNegative(0) + 1 
                else if ( peakHeight(ki) > iso(nIso) ) then
                    nNegative(nIso) = nNegative(nIso) + 1                
                else
                    do kk = 1,nIso-1
                        if (iso(kk+1) > peakHeight(ki)) then
                            nNegative(kk) = nNegative(kk) + 1
                            exit
                        end if
                    end do
                end if
            end do

        !---    combine regions and disappearing regions to find number of enclosures
            nEnclosure(0) = 1
            do kk = 1,nIso           
                nEnclosure(kk) = nRegions(kk) - nNegative(kk)
            end do 
            
        !---    combine regions and enclosures to find number of spots
            do kk = 1,nIso           
                nSpots(kk) = nRegions(kk) - nEnclosure(kk-1)
            end do
            
        !---    spot report
        !    write(*,fmt='(a8,a12,100a12)') "level","iso","nRegion","nNegative","nEnclosure","nSpots"
        !    do kk = 1,nIso
        !        write(*,fmt='(i8,f12.4,100i12)') kk,iso(kk),nRegions(kk),nNegative(kk),nEnclosure(kk),nSpots(kk)
        !    end do
        !    print *,""
            
            
            return
        end subroutine contiguousPeaks
    
    
    
    
    


    end program contourAnalyse