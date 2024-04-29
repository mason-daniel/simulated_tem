
program run_d2bi
    !---^^^^^^^^^^^^^^^^
    !*      a code for simulated TEM in the dynamical two beam imaging approximation
    !*      Daniel Mason, UKAEA, Feb 2021
    !*
     
    
            use iso_fortran_env
            use mpi_f08 
            
    
            use Lib_CommandLineArguments
            use Lib_XYZFiles
            use Lib_Filenames
             
            use Lib_SimpleProgressBar
            use Lib_Callipers
            use Lib_SimpleSupercells
            use Lib_ComplexSupercells
            use VoidIsosurfaces
            use NBAX_StringTokenizers
            use Lib_LinkCell3d
            use Lib_DeformationGradients
            use Lib_Lattices
            use Lib_RelativisticElectrons
            use Lib_DiffractionConditions
            use Lib_DynamicalTwoBeamImaging
            use Lib_ReadExtinctionDistances
            use Lib_Png
    
            implicit none
    
             
        !       version
        !   v.1.1.5     final version for paper
        !   v.1.1.6     cosmetic tweaks to screen output format, some unused functions removed
    
            character(len=8)                ::      VERSION = "1.1.6"
            real(kind=real64),parameter     ::      PI = 3.141592653589790d0 
            
            
        !---    command line run-time parameters
            
            type(CommandLineArguments)      ::      cla
    
            character(len=256)              ::      filename = ""
            character(len=256)              ::      outfile0 = ""
            real(kind=real64)               ::      a0 = 3.0d0
            character(len=8)                ::      latticeName = "bcc"
            real(kind=real64),dimension(3)  ::      a0in = LIB_CLA_NODEFAULT_R
            real(kind=real64),dimension(3)  ::      addOffset = LIB_CLA_NODEFAULT_R
            !logical                         ::      findXYZOffset = .false.
            !logical                         ::      removeOffset = .true.
            logical                         ::      xpad = .false.
            logical                         ::      ypad = .false.
            logical                         ::      zpad = .false.
            real(kind=real64),dimension(3)  ::      d0 = LIB_CLA_NODEFAULT_R
            logical                         ::      surf = .false.
            
            
            character(len=256)              ::      orientationFile = ""
            !character(len=256)              ::      defGradFile = ""
    
            real(kind=real64)               ::      V = 200.0d0                     !   electron voltage in keV
            real(kind=real64)               ::      xi0 = 103.889565d0              !   extinction distance in perfect crystal ( set to bcc W )
            real(kind=real64)               ::      xig = 207.533203d0              !   extinction distance in perfect crystal ( set to bcc W g = [110] )
            real(kind=real64),dimension(3)  ::      g3 = (/1,1,0/)
            real(kind=real64),dimension(3)  ::      k3 = (/0,0,1/)
            character(len=256)              ::      extinction_distances_filename = ""
            !real(kind=real64),dimension(3)  ::      z = LIB_CLA_NODEFAULT_R
            
            real(kind=real64),dimension(3,3)::      U = LIB_CLA_NODEFAULT_R
            real(kind=real64),dimension(3,3)::      R = LIB_CLA_NODEFAULT_R
            real(kind=real64)               ::      Sg = LIB_CLA_NODEFAULT_R        
            real(kind=real64)               ::      ng = LIB_CLA_NODEFAULT_R
            logical                         ::      dark = .true.
            
            real(kind=real64)               ::      tomoAngle = 20.0d0        !   degrees
            integer                         ::      ntomoAngle = 1
            
            integer                         ::      nPrecAngle = 1
            real(kind=real64)               ::      precAngle = 0.005d0       !   5 mrad
            
            
            integer                         ::      nAperture = 1
            real(kind=real64)               ::      apertureAngle = 0.0001d0   !   0.1 mrad
            
            integer                         ::      nDiffConditions = 1
            integer,parameter               ::      NDIFFMAX = 100
            real(kind=real64),dimension(4,NDIFFMAX) ::  k4array = 0 ,g4array = 0
            real(kind=real64),dimension(NDIFFMAX)             ::  ngarray = LIB_CLA_NODEFAULT_R
            real(kind=real64),dimension(NDIFFMAX)   ::  Sgarray = LIB_CLA_NODEFAULT_R
            logical,dimension(NDIFFMAX)             ::  darkarray = .true.
             
            logical                         ::      useDensity = .false.
            integer,dimension(2)            ::      M = (/4,4/)
            logical                         ::      opxyz = .false.
            logical                         ::      alterg = .false.
            
            
            logical                         ::      opPng = .true.
            real(kind=real64)               ::      png_min = 0.0d0,png_max = 0.0d0,png_blur = 2.50d0
            
            
            
    
        !---    timing data
    
            integer,parameter       ::      T_READ = 1
            integer,parameter       ::      T_CTOR = 2
            integer,parameter       ::      T_DEFGRAD = 3
            integer,parameter       ::      T_D2BI_INIT = 6
            integer,parameter       ::      T_OUTPUT = 7
            integer,parameter       ::      T_D2BI = 8
            integer,parameter       ::      T_TOTAL = 10
            type(Callipers),dimension(T_TOTAL)    ::      tt
    
    
        !---   important class instances
            type(XYZFile)                       ::      xyz
            type(SimpleSupercell)               ::      super       !   voxel spacing periodicity etc
            type(Lattice)                       ::      latt        !   expected lattice
            type(DynamicalTwoBeamImaging)       ::      d2bi
            
    
        !---    information about the super cell
            real(kind=real64),dimension(3,3)    ::      a_super,a_cell,a_vox , dfg  !   ,ia_cell 
            real(kind=real64),dimension(:,:),pointer            ::  x               !   (3,nAtoms) , atom positions
            real(kind=real64),dimension(3)      ::      offset 
            integer                             ::      Mx,My,Mz                     !   size of phase field
            integer                             ::      nAtoms
            logical                             ::      unitDefGrad = .false.
            logical                             ::      miller4 = .false.
    
        !---    dynamical 2 beam imaging
            real(kind=real64),dimension(:,:),allocatable    ::  d2bimg                  !   intensity map (0:w-1,0:h-1) , w,h is size of ouput image, usually order Mx,My
    
    
        !---    dummies
            logical                             ::      ok
            real(kind=real64),dimension(12)     ::      dat
            real(kind=real64),dimension(4*NDIFFMAX)     ::      darray
            character(len=256)                  ::      outfile
            integer                             ::      nDiff
            integer                             ::      ix,ii,kk,ll , nn
            character(len=5)                    ::      aaaa
            character(len=40)                   ::      bbbb
            character(len=8)                    ::      cccc
            real(kind=real64),dimension(3)      ::      xmin,xmax,yy,d
            real(kind=real64)                   ::      aperpx,modg,weight,pad , lenaxis1,lenaxis2
            character(len=256)                  ::      dummy
            integer,dimension(3,2)              ::      idat
            real(kind=real64),dimension(3,3)    ::      eps_bar,dfg_bar,rot_bar , U0,R0
            real(kind=real64)                   ::      theta
            
            integer,dimension(:),allocatable    ::      nGrain
            !logical                             ::      xyzHasDistance = .false.
            real(kind=real64),dimension(4)      ::      g4 
            real(kind=real64),dimension(4)      ::      k4 
            
            
            
            
            
            integer                             ::      nProcs,rank,ierror
            
            
            
            
            call MPI_INIT(ierror)
            call MPI_COMM_SIZE(MPI_COMM_WORLD, nProcs, ierror)
            call MPI_COMM_RANK(MPI_COMM_WORLD, rank, ierror)
          
            
            
            
        !---    start of program
            tt(T_TOTAL) = Callipers_ctor()
    
    
    
        !---    read command line arguments
            cla = CommandLineArguments_ctor(100)
            call setProgramDescription( cla, "d2bi  \n    simulated TEM in the dynamical two-beam imaging approximation " )
            call setProgramVersion( cla, VERSION )
            call setCategories(cla,(/   "input/output         ",       &
                                        "imaging conditions   ",       &
                                        "cell manipulation    ",       &
                                        ".png image output    ",       &
                                        "calculation          ",       &
                                        "microscope properties"        &
                                         /))
            call get( cla,"f",filename ,LIB_CLA_REQUIRED,"             input filename",1 )
            
    
            ii = 0
            call get( cla,"a0",a0in,ii ,LIB_CLA_REQUIRED,"            unit cell dimension (A)",1 )
            call get( cla,"lattice",latticeName ,LIB_CLA_OPTIONAL,"     lattice type",1 )
            if (getLatticeType(latticeName) == LATTICE_HCP) then
                a0 = a0in(1)
                if (ii==1) then
                    a0in(2) = a0in(1)
                    a0in(3) = a0in(1)*sqrt(8.0d0/3)
                end if
            else
                if (ii==1) then
                    a0in(2:3) = a0in(1)
                    a0 = a0in(1)
                else
                    a0 = ( a0in(1)+a0in(2)*a0in(3) )/3
                end if
            end if
    
            outfile0 = filename
            call get( cla,"o",outfile0 ,LIB_CLA_OPTIONAL,"           output filename",1)
    
           ! call get( cla,"DG",defGradFile,LIB_CLA_OPTIONAL,"          location of def grad field file ",1)
             
            call get( cla,"dp",orientationFile,LIB_CLA_OPTIONAL,"          location of diffraction pattern file ",1)
            call get( cla,"png",opPng,LIB_CLA_OPTIONAL,"         output .png file ",4)
            call get( cla,"png_min",png_min,LIB_CLA_OPTIONAL,"     .png intensity 0 level ",4)
            call get( cla,"png_max",png_max,LIB_CLA_OPTIONAL,"     .png intensity 1 level (set 0 for automatic)",4)
            call get( cla,"png_blur",png_blur,LIB_CLA_OPTIONAL,"    .png gaussian blur radius (A) ",4)
            
            
            
                            
    
           

           call get( cla,"xifile",extinction_distances_filename,LIB_CLA_OPTIONAL,"       location of extinction distances file",5  )
           call get( cla,"xi0",xi0,LIB_CLA_OPTIONAL,"          xi_0 (A)",5 )
           call get( cla,"xig",xig,LIB_CLA_OPTIONAL,"          xi_g (A)",5 )
    
           
            ii=0 ; darray=LIB_CLA_NODEFAULT_R ; call get( cla,"g",darray,ii,LIB_CLA_REQUIRED,"              g-vector (reduced recip space direction [hkl] or [hkil])" ,2)
            if (getLatticeType(latticeName) == LATTICE_HCP) then
                !   expect 4 miller-bravais indices
                if (mod(ii,4) == 0) then
                    nDiffConditions = int( ii / 4 )
                    do ii = 1,nDiffConditions
                        g4array(1:4,ii) = darray( ii*4-3:ii*4 )
                    end do
                else
                    call errorExit("d2bi error - can't parse array of miller-bravais index g-vectors.")
                end if
            else 
                !   expect 3 miller indices
                if (mod(ii,3) == 0) then
                    nDiffConditions = int( ii / 3 )
                    do ii = 1,nDiffConditions
                        g4array(1:3,ii) = darray( ii*3-2:ii*3 )
                    end do
                else
                    call errorExit("d2bi error - can't parse array of miller index g-vectors.")
                end if
            end if   
    !         end if
    
            ii=0 ; darray=LIB_CLA_NODEFAULT_R ; call get( cla,"k",darray,ii,LIB_CLA_REQUIRED,"              zone axis (cell space direction, [uvw] or [uvtw])",2 )
            if (getLatticeType(latticeName) == LATTICE_HCP) then
                !   expect 4 miller-bravais indices
                miller4 = .true.
                if (mod(ii,4) == 0) then
                    if ( int( ii / 4 )==nDiffConditions ) then
                        do ii = 1,nDiffConditions
                            k4array(1:4,ii) = darray( ii*4-3:ii*4 )
                        end do
                    else if ( int( ii / 4 )==1 ) then
                        do ii = 1,nDiffConditions
                            k4array(1:4,ii) = darray( 1:4 )
                        end do
                    else if (nDiffConditions == 1) then
                        nDiffConditions = int( ii / 4 )
                        do ii = 1,nDiffConditions
                            k4array(1:4,ii) = darray( ii*4-3:ii*4 )
                            g4array(1:4,ii) = g4array( 1:4,ii )
                        end do
                    else
                        call errorExit("d2bi error - array of miller-bravais index k-vectors must match g-vectors. ")
                    end if
                else
                    call errorExit("d2bi error - can't parse array of miller-bravais index k-vectors.")
                end if
            else 
                !   expect 3 miller indices
                miller4 = .false.
                if (mod(ii,3) == 0) then
                    if ( int( ii / 3 )==nDiffConditions ) then
                        do ii = 1,nDiffConditions
                            k4array(1:3,ii) = darray( ii*3-2:ii*3 )
                        end do
                    else if ( int( ii / 3 )==1 ) then
                        do ii = 1,nDiffConditions
                            k4array(1:3,ii) = darray( 1:3 )
                        end do
                    else if (nDiffConditions == 1) then
                        nDiffConditions = int( ii / 3 )
                        do ii = 1,nDiffConditions
                            k4array(1:3,ii) = darray( ii*3-2:ii*3 )
                            g4array(1:3,ii) = g4array( 1:3,ii )
                        end do                    
                    else 
                        call errorExit("d2bi error - array of miller-bravais index k-vectors must match g-vectors. ")
                    end if
                else
                    call errorExit("d2bi error - can't parse array of miller-bravais index k-vectors.")
                end if
            end if    
    
    
    
    
    
           !ii=3 ; call get( cla,"z",z,ii,LIB_CLA_OPTIONAL,"     electron direction (real space direction to be tweaked by deviation parameter)",2 )
           
    
            !call get( cla,"sg",Sg,LIB_CLA_OPTIONAL,  "    value for deviation parameter s_g (leave unset for fixed k-vector)",2  )
            !call get( cla,"ng",ng,LIB_CLA_OPTIONAL,"    value for deviation parameter (g,ng)",2  )
            
            ii=0 ; call get( cla,"sg",sgarray,ii,LIB_CLA_OPTIONAL,  "           value for deviation parameter s_g (leave unset for fixed k-vector)",2  )              
            if ( ii==1 ) then
                sgarray(1:nDiffConditions) = sgarray(1)        
            else if ( (nDiffConditions==1).and.(ii>1) ) then
                nDiffConditions = ii
                do ii = 1,nDiffConditions
                    k4array(1:4,ii) = k4array( 1:4,1 )
                    g4array(1:4,ii) = g4array( 1:4,1 )
                end do          
            else if ( ii*(ii-nDiffConditions)/=0 ) then
                call errorExit("d2bi error - array of -sg deviation parameters must match k,g-vectors ")
            end if    
        
            
            ii=0 ; call get( cla,"ng",ngarray,ii,LIB_CLA_OPTIONAL,"           value for deviation parameter defined by (g,ng)",2  )
            if ( ii==1 )  then
                ngarray(1:nDiffConditions) = ngarray(1)   
            else if ( (nDiffConditions==1).and.(ii>1) ) then
                nDiffConditions = ii
                do ii = 1,nDiffConditions
                    k4array(1:4,ii) = k4array( 1:4,1 )
                    g4array(1:4,ii) = g4array( 1:4,1 )
                end do
                sgarray(1:nDiffConditions) = sgarray(1)       
            else if ( ii*(ii-nDiffConditions)/=0 ) then
                call errorExit("d2bi error - array of -ng parameters must match k,g-vectors ")
            end if    
            
            ii=0 ; darkarray=.true. ; call get( cla,"dark",darkarray,ii,LIB_CLA_OPTIONAL,"         dark field imaging mode ",2  )
            if ( ii==1 )  then
                darkarray(1:nDiffConditions) = darkarray(1)   
            else if ( (nDiffConditions==1).and.(ii>1) ) then
                nDiffConditions = ii
                do ii = 1,nDiffConditions
                    k4array(1:4,ii) = k4array( 1:4,1 )
                    g4array(1:4,ii) = g4array( 1:4,1 )
                end do
                sgarray(1:nDiffConditions) = sgarray(1)       
                ngarray(1:nDiffConditions) = ngarray(1)              
            else if ( ii*(ii-nDiffConditions)/=0 ) then
                call errorExit("d2bi error - array of -dark settings must match k,g-vectors ")
            end if 
            call get( cla,"theta",THETA_MAX,LIB_CLA_OPTIONAL,"        maximum tilt search angle (deg) ",2 )
            !if (hasArgument(cla,"theta")) THETA_MAX = THETA_MAX*PI/180.0d0
 
            
            call get( cla,"V",V,LIB_CLA_OPTIONAL,"            accelerating voltage (keV)",6)
            call get( cla,"nPrecAngle",nPrecAngle ,LIB_CLA_OPTIONAL,"   number of points take for convergent beam precession",6  )
            call get( cla,"precAngle",precAngle,LIB_CLA_OPTIONAL,"    angle (radians) for convergent beam precession",6  )
           
            ! call get( cla,"nAperture ",nAperture ,LIB_CLA_OPTIONAL,"    number of points for aperture (try 1,8,17)",6  )
            ! call get( cla,"apertureAngle",apertureAngle,LIB_CLA_OPTIONAL,"angle (radians) for aperture",6  )
           
            call get( cla,"nTomoAngle",ntomoAngle,LIB_CLA_OPTIONAL,"   tomography tilt steps",6 )
            call get( cla,"tomoAngle",tomoAngle,LIB_CLA_OPTIONAL,"    tomography tilt half angle (deg)",6 )
             
    
            
            
            
            call get( cla,"density",useDensity,LIB_CLA_OPTIONAL,"      use atomic density field in TEM calc",5  )
            call get( cla,"alterg",alterg,LIB_CLA_OPTIONAL,"       use deformation gradient/diffraction pattern to adjust reciprocal lattice vectors ",5  )
            
    
            ii = 0; call get( cla,"M",M,ii ,LIB_CLA_OPTIONAL,"            voxel grid per unit cell ",5 )
            if (hasArgument(cla,"M")) then
                if (ii == 1) then
                    M(1:2) = M(1)
                else if (ii /= 2) then
                    call errorExit("d2bi error - expected -M xdiv or -M xdiv,zdiv")
                end if
            end if     

            call get( cla,"xpad",xpad,LIB_CLA_OPTIONAL,"         while reading input file, add padding in x-direction to create a surface ",3  )
            call get( cla,"ypad",ypad,LIB_CLA_OPTIONAL,"         while reading input file, add padding in y-direction to create a surface ",3  )
            call get( cla,"zpad",zpad,LIB_CLA_OPTIONAL,"         while reading input file, add padding in z-direction to create a surface ",3  )
            
            call get( cla,"surf",surf,LIB_CLA_OPTIONAL,"         check atom extents for foil thickness instead of box size ",3  )
            ii=0 ; dat = LIB_CLA_NODEFAULT_R ; call get( cla,"d",dat,ii,LIB_CLA_OPTIONAL,"            imaging space extents. Set 0 to use supercell, unset to use atom extents." ,3)
            if (ii == 1) then
                d0(1:2) = 0.0d0
                d0(3) = dat(1)
            else if (ii==2) then
                d0(1:2) = dat(1:2)
                d0(3) = 0.0d0
            else if (ii==3) then
                d0(1:3) = dat(1:3)
            end if


            !call get( cla,"findXYZOffset",findXYZOffset ,LIB_CLA_OPTIONAL,"find offset of atoms in .xyz file",3 )
            !call get( cla,"offset",removeOffset ,LIB_CLA_OPTIONAL,"       remove any global supercell offsets generated internally",3 )
            ii=3 ; call get( cla,"addOffset",addOffset ,ii,LIB_CLA_OPTIONAL,"    add a global atom position offset (in order to stitch together images)",3 )
    
            ii=9 ; dat = LIB_CLA_NODEFAULT_R ; call get( cla,"U",dat,ii,LIB_CLA_OPTIONAL,"            force crystal rotation matrix in column-major order" ,3)
            if (hasArgument(cla,"U")) then
                U = reshape( dat , (/3,3/) )
                if (abs( determinant3Mat(U)-1 )>0.001d0) then
                    call errorExit("d2bi error - input U matrix has a determinant /= 1")           
                end if
            end if
    
            ii=9 ; dat = LIB_CLA_NODEFAULT_R ; call get( cla,"R",dat,ii,LIB_CLA_OPTIONAL,"            force foil tilt rotation matrix in column-major order" ,3)
            if (hasArgument(cla,"R")) then
                R = reshape( dat , (/3,3/) )
                if (abs( determinant3Mat(R)-1 )>0.001d0) then
                    call errorExit("d2bi error - input R matrix has a determinant  /= 1")                         
                end if
            end if
            call get( cla,"opxyz",opxyz,LIB_CLA_OPTIONAL,"        output .xyz file of atom positions after rotation to two-beam condition",3  )

            
    


            !call get( cla,"bright",ok,LIB_CLA_OPTIONAL,"  bright field imaging mode ",2  )
            !if (hasArgument(cla,"bright")) dark = .not. ok
    
    
    
            if (rank==0) &
            call report(cla)
            if (.not. allRequiredArgumentsSet(cla)) call errorExit("")
            if (hasHelpArgument(cla)) call errorExit("")
    
    
    
    
    
    
     
    
    
    !    !---    check for dependencies in input options
    !
            inquire(file=trim(filename),exist=ok)
            if (.not. ok) then
                call errorExit("d2bi error - could not find input file """//trim(filename)//"""")            
            end if
            
           if (hasArgument( cla,"png_min" ).or. hasArgument( cla,"png_max" ).or. hasArgument( cla,"png_blur" )) then
               opPng = .true.
           end if
            
            
            if (hasArgument( cla,"xifile" )) then
                 if (hasArgument( cla,"xi0" ) .or. hasArgument( cla,"xig" )) call errorExit("d2bi error - attempt to define xi_0 or xi_g and read from file")
            end if
    
    !
     
           if (hasArgument( cla,"ng" )) then
               if (hasArgument( cla,"sg" )) stop "d2bi error - attempt to define both n_g and s_g"
           end if
    
    
            if ( nPrecAngle  < 1 ) then
                if (rank==0) write(*,fmt='(a)') "d2bi warning - attempt to define < 1 precession angle?"
                if (rank==0) write(*,fmt='(a)') "    setting -nPrecAngle  1"            
                nPrecAngle  = 1
            end if
            
           
            if ( (len_trim(outfile0)==0).and. (opxyz .or. oppng) ) then
                if (rank==0) write(*,fmt='(a)') "d2bi warning - conflicting options -o """" and -opxyz / -png"
                if (rank==0) write(*,fmt='(a)') "    setting -noopxyz -nopng"
                opxyz = .false.
                oppng = .false.
            end if
             
            
            oppng = oppng .and. (rank==0)
            opxyz = opxyz .and. (rank==0)
            
    
            if (.not. dark) then
                if ( (ng/=LIB_CLA_NODEFAULT_R).and.(ng>1) ) then
                    if (rank==0) write(*,fmt='(a)') "d2bi warning - conflicting options '-nodark/-bright' and '-ng > 1'"
                    if (rank==0) write(*,fmt='(a)') "   bright field computed by computing [hkl] point at max intensity, and inverting image."
                    if (rank==0) write(*,fmt='(a)') "   setting ng = 1."
                    ng = 1 
                end if            
            end if
            
    
            if ((xpad .or. ypad .or. zpad) .and. .not. surf) then
                if (rank==0) write(*,fmt='(a)') "d2bi warning - conflicting options '-(xyz)pad' and '-nosurf'" 
                if (rank==0) write(*,fmt='(a)') "   generating an unrelaxed surface during file input, requires foil thickness to be calculated from atom positions. Applying '-surf'."
                surf = .true.
            end if        
            
            if (xpad .and. (d0(1) == LIB_CLA_NODEFAULT_R) ) then
                if (rank==0) write(*,fmt='(a)') "d2bi warning - conflicting options '-xpad' and '-d' not fully set" 
                if (rank==0) write(*,fmt='(a)') "    using x supercell extent after padding" 
                d0(1) = 0            
            end if        
            
            if (ypad .and. (d0(2) == LIB_CLA_NODEFAULT_R) ) then
                if (rank==0) write(*,fmt='(a)') "d2bi warning - conflicting options '-ypad' and '-d' not fully set" 
                if (rank==0) write(*,fmt='(a)') "    using y supercell extent after padding" 
                d0(2) = 0            
            end if      
            
            if (zpad .and. (d0(3) == LIB_CLA_NODEFAULT_R) ) then
                if (rank==0) write(*,fmt='(a)') "d2bi warning - conflicting options '-zpad' and '-d' not fully set" 
                if (rank==0) write(*,fmt='(a)') "    using z supercell extent after padding" 
                d0(3) = 0            
            end if      
            
            
            if ((ntomoAngle > 1).and. .not. (xpad .or. ypad .or. zpad) ) then
                if (rank==0) write(*,fmt='(a)') "d2bi warning - options '-ntomoAngle' and '-noxpad -noypad -nozpad'" 
                if (rank==0) write(*,fmt='(a)') "   Need an explicit surface calculation to do tomography sensibly. Assuming that the input file has buffer space in direction parallel to k-vector"
            end if        
            
            
            if ((ntomoAngle > 1).and. .not. surf) then
                if (rank==0) write(*,fmt='(a)') "d2bi warning - conflicting options '-ntomoAngle' and '-nosurf'" 
                if (rank==0) write(*,fmt='(a)') "   Need an explicit surface calculation to do tomography sensibly. Setting -surf"
                surf = .true.
            end if        
     
            
            ! if ( len_trim(defGradFile) /= 0 ) then
            !     inquire(file=trim(defGradFile),exist=ok)
            !     if (.not. ok) then
            !         call errorExit("d2bi error - could not find def grad file -DG """//trim(defGradFile)//"""")            
            !     end if
            !     if (.not. alterg) then
            !         if (rank==0) write(*,fmt='(a)') "d2bi warning - conflicting options '-DG' and '-noalterg'" 
            !         if (rank==0) write(*,fmt='(a)') "    no point reading the deformation gradient and not using it" 
            !         if (rank==0) write(*,fmt='(a)') "    setting -alterg" 
            !         alterg = .true.
            !     end if
            ! end if
            
            
            
            if ( len_trim(orientationFile) /= 0 ) then
                inquire(file=trim(orientationFile),exist=ok)
                if (.not. ok) then
                    call errorExit("d2bi error - could not find orientation file -dp """//trim(orientationFile)//"""")            
                end if
                if (.not. alterg) then
                    if (rank==0) write(*,fmt='(a)') "d2bi warning - conflicting options '-dp' and '-noalterg'" 
                    if (rank==0) write(*,fmt='(a)') "    no point reading the diffraction pattern and not using it" 
                    if (rank==0) write(*,fmt='(a)') "    setting -alterg" 
                    alterg = .true.
                end if
            end if
            
            
            call delete(cla)
    
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            tt(T_OUTPUT) = Callipers_ctor()
            call pause( tt(T_OUTPUT) )
    
        !---    friendly welcome give the runtime parameters
            if (rank==0) then
                print *,""
                print *,"d2bi v"//trim(VERSION)
                print *,"^^^^^^"//repeat("^",len_trim(VERSION))
                print *,""
                call get_environment_variable("HOSTNAME",dummy)
                print *,"running on ",nProcs," processors on """//trim(dummy)//""""
                print *,""
                print *,"options"
                print *,"^^^^^^^"
                print *,""
                print *,"input filename     : ",trim(filename)
                if (len_trim(outfile0)==0) then
                    print *,"no output files"
                else
                    print *,"output file prefix : ",trim(outfile0)
                end if
                write(*,fmt='(a,3f10.5,a)') " a0                 : ",a0in," (A)"
                print *,"lattice            : """//trim(latticeName)//""""
                print *,"M (voxels per cell): ",M
        
                !print *,"findXYZOffset      : ",findXYZOffset
     
                if (len_trim(orientationFile)/=0) then
                    print *,"using orientation file """//trim(orientationFile)//""""
                else 
                    print *,"no orientation file set"
                end if
                ! if (len_trim(defGradFile)/=0) then
                !     print *,"using def grad file """//trim(defGradFile)//""""
                ! else 
                !     print *,"no deformation gradient file set"
                ! end if
                print *,"output .png ?      : ",opPng
                if (opPng) &
                write(*,fmt='(a,3f12.5)') " .png min/max/blur  : ",png_min,png_max,png_blur
                
                print *,"adjust recip latt  : ",alterg
        
        
                write (*,fmt='(a,f10.3,a)') " electron voltage  : ",V," (keV)"
                if (extinction_distances_filename=="") then
                    write(*,fmt='(a,2f16.6,a)') "    xi_g(A),xi_0(A) ",xig,xi0," (A)"
                else
                    write(*,fmt='(a)') " extinction distance file """//trim(extinction_distances_filename)//""""
                end if
                
                print *,"number of diffraction conditions ",nDiffConditions
                do ii = 1,nDiffConditions   
                    print *,"    condition ",ii
                    if (miller4) then
                        write(*,fmt='(a,4f12.3)') "                     zone axis k = [uvtw] ",k4array(1:4,ii)
                        write(*,fmt='(a,4f12.3)') "                     g-vector [hkil]      ",g4array(1:4,ii)
                    else
                        write(*,fmt='(a,3f12.3)') "                     zone axis k = [uvw] ",k4array(1:3,ii)
                        write(*,fmt='(a,3f12.3)') "                     g-vector [hkl]      ",g4array(1:3,ii)
                    end if
                    if (Sgarray(ii) /= LIB_CLA_NODEFAULT_R) write(*,fmt='(a,f12.3)') "                   s_g                   ",Sgarray(ii)
                    if (ngarray(ii) /= LIB_CLA_NODEFAULT_R) write(*,fmt='(a,f12.3)') "                    ng                   ",ngarray(ii)
                    if (darkarray(ii)) then
                        print *,"                    dark field imaging  "               
                    else
                        print *,"                    bright field imaging   "
                    end if
        
                end do
                
                if (nPrecAngle>1) then
                    print *,"    number of angles for convergent beam precession ",nPrecAngle 
                    write (*,fmt='(a,f10.5,a)') "    angle for convergent beam precession  ",precAngle," (rad)"
                else
                    print *,"    no convergent beam precession"
                end if
               
                ! if (nAperture>1) then
                !     print *,"    number of points sampling aperture ",nAperture
                !     write (*,fmt='(a,f10.5,a)') "     angle for aperture ",apertureAngle," (rad)"
                ! else
                !     print *,"    using pinprick aperture "
                ! end if
               
               if (xpad) print *,"    creating surface by padding normal to x- direction in input .xyz file"
               if (ypad) print *,"    creating surface by padding normal to y- direction in input .xyz file"
               if (zpad) print *,"    creating surface by padding normal to z- direction in input .xyz file"
        
    
                if (d0(1) == LIB_CLA_NODEFAULT_R) then
                    print *,"    voxel extent x- direction set by atomic extent"
                else if (d0(1) == 0.0d0) then
                    if (xpad) then
                        print *,"    voxel extent x- direction set by supercell after padding"           
                    else             
                        print *,"    voxel extent x- direction set by input supercell"           
                    end if
                else
                    print *,"    voxel extent x- direction set to ",d0(1)                    
                end if
                
    
                if (d0(2) == LIB_CLA_NODEFAULT_R) then
                    print *,"    voxel extent y- direction set by atomic extent"
                else if (d0(2) == 0.0d0) then
                    if (ypad) then
                        print *,"    voxel extent y- direction set by supercell after padding"                
                    else             
                        print *,"    voxel extent y- direction set by input supercell"           
                    end if                    
                else
                    print *,"    voxel extent y- direction set to ",d0(2)                    
                end if
    
                if (d0(3) == LIB_CLA_NODEFAULT_R) then
                    print *,"    voxel extent z- direction set by atomic extent"
                else if (d0(3) == 0.0d0) then
                    if (zpad) then
                        print *,"    voxel extent z- direction set by supercell after padding"                
                    else             
                        print *,"    voxel extent z- direction set by input supercell"           
                    end if                    
                else
                    print *,"    voxel extent z- direction set to ",d0(3)                    
                end if                        
                
                
    !  
    
                 if (surf) print *,"                    atom extents used to determine true foil thickness for dark field tilting"
        
                 
        
                if (U(1,1) /= LIB_CLA_NODEFAULT_R) then
                    print *,"    crystal rotation"
                    write(*,fmt='(a,3f16.12)') "                         ",U(1,:)
                    write(*,fmt='(a,3f16.12)') "                         ",U(2,:)
                    write(*,fmt='(a,3f16.12)') "                         ",U(3,:)
                end if
        
                if (R(1,1) /= LIB_CLA_NODEFAULT_R) then
                    print *,"    foil tilt"
                    write(*,fmt='(a,3f16.12)') "                         ",R(1,:)
                    write(*,fmt='(a,3f16.12)') "                         ",R(2,:)
                    write(*,fmt='(a,3f16.12)') "                         ",R(3,:)
                else
                    write (*,fmt='(a,f10.5,a)') "     tilt stage max angle ",THETA_MAX," (deg)"
                end if
                if (ntomoAngle>1) then
                    write (*,fmt='(a,f10.5,a)') "     tomography max angle ",tomoAngle," (deg)"
                    print *,"    tomography steps     ",ntomoAngle
                else
                    print *,"    single orientation (no tomography)     "
                end if
        
         
               print *,"    output .xyz file after rotation + tilt ",opxyz
               
               print *,"    use atomic density field   ",useDensity
                
        
               
                if (addOffset(1) /= LIB_CLA_NODEFAULT_R) then
                    print *,"add global offset : ",addOffset
                end if
                print *,""
            end if
    
    
    
            if (getLatticeType(latticeName) == LATTICE_HCP) then
                latt = Lattice_ctor( latticeName,(sqrt(8.0d0/3)+sqrt(3.0d0))/2 )
                a0 = a0in(1)
                call setCoverA(latt,a0in(3)/a0in(1))       
            else
                latt = Lattice_ctor( latticeName,(1.0d0+sqrt(2.0d0))/2 )
                a0 = ( a0in(1) + a0in(2) + a0in(3) )/3
            end if
    
            
    
    
    
        !---    read in the .xyz file - rank 0 does the reading then broadcasts
            tt(T_READ) = Callipers_ctor()
            if (rank==0) then
                
                print *,""
                print *,"reading input file"
                print *,"^^^^^^^^^^^^^^^^^^"
                print *,""
                
                
                xyz = XYZFile_ctor(filename)
                call readHeader(xyz,ok)
                call input(xyz,verbose=(rank==0))
                offset = 0
                ! if (findXYZOffset) then
                !     offset = findOffset(xyz,a0,offset=0.25d0)
                !     call adjustOffset(xyz,a0,offset=0.25d0,verbose=(rank==0))
                !     print *,"d2bi info - offset adjusted from file ",offset
                ! else
                !     print *,"d2bi info - not adjusting offset in .xyz file"
                !     offset = 0
                ! end if
                call report(xyz)
                call getSupercell(xyz,a_super,ok)
                if (ok) then            
                    print *,"d2bi info - supercell read from file"
                    write(*,fmt='(3f16.6)') a_super(1,:)
                    write(*,fmt='(3f16.6)') a_super(2,:)
                    write(*,fmt='(3f16.6)') a_super(3,:)
                    write(*,fmt='(a,f16.4,a)') "d2bi info - supercell volume ",superCellVolume(a_super)," (A^3)"
                else
                    call errorExit( "d2bi error - supercell not read from file" )
                end if
                nAtoms = getNatoms(xyz)
                
                        
                if (getnColumns(xyz)<14) then
                    print *,"d2bi warning - expecting .xyz format"
                    print *,"   atom pos_x pos_y pos_z grain dg11,dg21,dg31,dg12,dg22,dg32,dg13,dg23,dg33 sublattice [distance]"
                    print *,"assuming all grain = 1"
                    print *,"assuming all def grad = (/ 1,0,0,0,1,0,0,0,1 /)"
                    print *,"assuming all sublattice = 1"
                    call setNcolumns(xyz,14)
                    call getColumnsp( xyz,x )
                    do ii = 1,nAtoms
                        x(4,ii) = 1
                        x(5:13,ii) = (/ 1,0,0,0,1,0,0,0,1 /)
                        x(14,ii) = 1
                    end do
                    unitDefGrad = .true.
                end if
                call getColumnsp( xyz,x )
                if (addOffset(1) /= LIB_CLA_NODEFAULT_R) then
                    do ii = 1,nAtoms
                        x(1:3,ii) = x(1:3,ii) + addOffset(1:3)
                    end do
                end if
                
                
            !---    check for periodic boundary conditions - useful for PARCAS format!          
                
                Mx = max(3,floor( M(1)*norm2(a_super(:,1))/a0 ))
                My = max(3,floor( M(1)*norm2(a_super(:,2))/a0 ))
                Mz = max(3,floor( M(2)*norm2(a_super(:,3))/a0 ))
                a_vox(1:3,1) = a_super(1:3,1) / Mx
                a_vox(1:3,2) = a_super(1:3,2) / My
                a_vox(1:3,3) = a_super(1:3,3) / Mz
                super = SimpleSupercell_ctor(a_vox,Mx,My,Mz)
                  
                do ii = 1,nAtoms
                    x(:,ii) = wrapPBC( super, x(:,ii) )
                end do
    
                
                
            !---    optional request for extra padding of input .xyz file to create a surface            
                 
                if (xpad .or. ypad .or. zpad) then
    !                 if (nTomoAngle > 1) then
    !                     if (nPrecAngle  > 1) then
    !                         theta = tomoAngle*PI/180.0 + THETA_MAX + precAngle  
    !                     else
    !                         theta = tomoAngle*PI/180.0 + THETA_MAX  
    !                     end if                                      
    !                 else
    !                     if (nPrecAngle  > 1) then
    !                         theta = THETA_MAX + precAngle  
    !                     else
    !                         theta = THETA_MAX  
    !                     end if
    !                 end if
                    
                    if (nTomoAngle > 1) then
                        theta = ( tomoAngle + THETA_MAX )*PI/180.0 + precAngle  
                    else
                        theta = THETA_MAX*PI/180.0 + precAngle  
                    end if
                    
                    
                    print *,""
                    print *,"d2bi info - creating a surface"
                    print *,"    assuming max angle ",theta*180/PI," (deg)"
                end if
            
            
                if (xpad) then
                    !   increase x- supercell repeat length by 50%, and centre atoms
    !                 print *,"d2bi info - creating a surface normal to x"
    !                 pad = norm2(a_super(1,:)) / 4
    !                 x(1,1:nAtoms) = x(1,1:nAtoms) + pad
    !                 a_super(1,:) = a_super(1,:)*1.5d0   
    
                    print *,"d2bi info - creating a surface normal to x"
                    lenaxis1 = norm2(a_super(1,:))
                    lenaxis2 = max( norm2(a_super(2,:)),norm2(a_super(3,:)) )
                    if (d0(1)/=LIB_CLA_NODEFAULT_R) lenaxis2 = max(lenaxis2,d0(1)/cos(theta))
                    if (d0(2)/=LIB_CLA_NODEFAULT_R) lenaxis2 = max(lenaxis2,d0(2)/cos(theta))
                    write (*,fmt='(4(a,f12.3))') "    pad dimension ",lenaxis1," (A) perpendicular dimension ",lenaxis2
                    pad = ( lenaxis1*cos(theta) + lenaxis2*sin(theta) )                  
                    x(1,1:nAtoms) = x(1,1:nAtoms) + ( pad - lenaxis1 )*1.05/2
                    a_super(1,:) = a_super(1,:)*(pad*1.05d0/lenaxis1)
                end if
                
                if (ypad) then
                    !   increase y- supercell repeat length by 50%, and centre atoms
                    print *,"d2bi info - creating a surface normal to y"
                    lenaxis1 = norm2(a_super(2,:))
                    lenaxis2 = max( norm2(a_super(1,:)),norm2(a_super(3,:)) )
                    if (d0(1)/=LIB_CLA_NODEFAULT_R) lenaxis2 = max(lenaxis2,d0(1)/cos(theta))
                    if (d0(2)/=LIB_CLA_NODEFAULT_R) lenaxis2 = max(lenaxis2,d0(2)/cos(theta))
                    write (*,fmt='(4(a,f12.3))') "    pad dimension ",lenaxis1," (A) perpendicular dimension ",lenaxis2
                    pad = ( lenaxis1*cos(theta) + lenaxis2*sin(theta) )                  
                    x(2,1:nAtoms) = x(2,1:nAtoms) + ( pad - lenaxis1 )*1.05/2
                    a_super(2,:) = a_super(2,:)*(pad*1.05d0/lenaxis1)
                    
                end if
    
                if (zpad) then
    !                 !   increase z- supercell repeat length by 50%, and centre atoms
    !                 print *,"d2bi info - creating a surface normal to z"
    !                 pad = norm2(a_super(3,:)) / 4
    !                 x(3,1:nAtoms) = x(3,1:nAtoms) + pad
    !                 a_super(3,:) = a_super(3,:)*1.5d0   
    
                    print *,"d2bi info - creating a surface normal to z"
                    lenaxis1 = norm2(a_super(3,:))
                    lenaxis2 = max( norm2(a_super(1,:)),norm2(a_super(2,:)) )
                    if (d0(1)/=LIB_CLA_NODEFAULT_R) lenaxis2 = max(lenaxis2,d0(1)/cos(theta))
                    if (d0(2)/=LIB_CLA_NODEFAULT_R) lenaxis2 = max(lenaxis2,d0(2)/cos(theta))
                    write (*,fmt='(4(a,f12.3))') "    pad dimension ",lenaxis1," (A) perpendicular dimension ",lenaxis2
                    pad = ( lenaxis1*cos(theta) + lenaxis2*sin(theta) )                  
                    x(3,1:nAtoms) = x(3,1:nAtoms) + ( pad - lenaxis1 )*1.05/2
                    a_super(3,:) = a_super(3,:)*(pad*1.05d0/lenaxis1)
    
                end if
                
                if (xpad .or. ypad .or. zpad) then
                    write (*,fmt='(4(a,f12.3))')"    padding needed ", ( pad - lenaxis1 )," (+5% = ",( pad - lenaxis1 )*1.05," )"
                    print *,"d2bi info - supercell after padding"
                    print *,a_super(1,:)
                    print *,a_super(2,:)
                    print *,a_super(3,:)
                    print *,"d2bi info - supercell volume ",superCellVolume(a_super)      
                    if ( (d0(3) == 0).or.(d0(3) == LIB_CLA_NODEFAULT_R) ) then
                        print *,"d2bi info - setting voxel supercell dimension 3 to ",pad
                        d0(3) = pad 
                    end if
                end if
                
                
                print *,""
            end if
            
        !---    broadcast  
            call MPI_BCAST( nAtoms, 1, MPI_INTEGER , 0, MPI_COMM_WORLD, ierror )        
            call MPI_BCAST( a_super, 9, MPI_DOUBLE_PRECISION , 0, MPI_COMM_WORLD, ierror )      
            if (rank/=0) allocate( x(3,nAtoms) )
            call MPI_BCAST( x(1:3,1:nAtoms), 3*nAtoms, MPI_DOUBLE_PRECISION , 0, MPI_COMM_WORLD, ierror )       
            if (xpad .or. ypad .or. zpad)  call MPI_BCAST( d0, 3, MPI_DOUBLE_PRECISION , 0, MPI_COMM_WORLD, ierror )       
           
            call pause(tt(T_READ))
            
     
    
         !---    set up d2bi calculation
             
            tt(T_CTOR) = Callipers_ctor()
    
            Mx = max(3,floor( M(1)*norm2(a_super(:,1))/a0 ))
            My = max(3,floor( M(1)*norm2(a_super(:,2))/a0 ))
            Mz = max(3,floor( M(2)*norm2(a_super(:,3))/a0 ))
            a_vox(1:3,1) = a_super(1:3,1) / Mx
            a_vox(1:3,2) = a_super(1:3,2) / My
            a_vox(1:3,3) = a_super(1:3,3) / Mz
            super = SimpleSupercell_ctor(a_vox,Mx,My,Mz)        
            
            a_cell = a0*getConventionalCell(latt)
            if (rank==0) then
                print *,"d2bi info - conventional unit cell expected from lattice type"
                write (*,fmt='(3f12.6)') a_cell(1,:)
                write (*,fmt='(3f12.6)') a_cell(2,:)
                write (*,fmt='(3f12.6)') a_cell(3,:)
                print *,""
            end if
    
            call pause(tt(T_CTOR))
    
      
    
        !---    may need to compute the average def grad, if alter_g is set.
            dfg_bar = reshape( (/1,0,0,0,1,0,0,0,1/),(/3,3/) )    
            eps_bar = reshape( (/0,0,0,0,0,0,0,0,0/),(/3,3/) )
            rot_bar = reshape( (/1,0,0,0,1,0,0,0,1/),(/3,3/) )
        
            tt(T_DEFGRAD) = Callipers_ctor()
            if (alterg) then
                    
    
                ! if ( len_trim(defGradFile) /= 0 ) then
                
                !     !   read spatially varying deformation gradient field from disk. Rank 0 reads then broadcasts
                !     if (rank==0) then
                !         print *,""
                !         print *,"reading def grad field"
                !         print *,"^^^^^^^^^^^^^^^^^^^^^^"
                !         print *,""
                !         print *,"filename """//trim(defGradFile)//""" (binary file)"
                   
                         
                !         call computeAvgDefGrad(eps_bar,dfg_bar,weight)
             
                !         open(unit=808,file=trim(defGradFile),action="read",form="unformatted")
                !         !---    old line 22/05/23            
                !         !   read(unit=808) ii
                !         !   if (ii /= Mx*My*Mz) then
                !         !       new line
                !             read(unit=808) dummy(1:10)
                !             xyzHasDistance = (dummy(9:10) == "+d")
                !             read(unit=808) ix,iy,iz
                !             if ( .not. ( (ix==Mx).and.(iy==My).and.(iz==Mz) ) ) then                    
                !                 if (rank==0) print *,"d2bi warning - phase field allocated has ",Mx,",",My,",",Mz," nodes, phase field in file has ",ix,",",iy,",",iz," nodes"                        
                !             end if    
                            
                !             if (xyzHasDistance) then
                !                 if (rank==0) print *,"reading column for distance"                        
                !                 do iz = 0,Mz-1
                !                     do iy = 0,My-1
                !                         do ix = 0,Mx-1
                !                             call progressBar( ix + 1+ Mx*(iy + My*iz)+1,Mx*My*Mz )
                !                             read (unit=808) dat(1:9),dp 
                !                             call computeAvgDefGrad( dat(1:9),1.0d0, eps_bar,dfg_bar,weight)
                !                         end do
                !                     end do
                !                 end do
                !             else
                !                 do iz = 0,Mz-1
                !                     do iy = 0,My-1
                !                         do ix = 0,Mx-1
                !                             call progressBar( ix + 1+ Mx*(iy + My*iz)+1,Mx*My*Mz )
                !                             read (unit=808) dat(1:9)
                !                             call computeAvgDefGrad( dat(1:9),1.0d0, eps_bar,dfg_bar,weight)
                !                         end do
                !                     end do
                !                 end do
                !             end if
                !         close(unit=808)
                !         call computeAvgDefGrad(eps_bar,dfg_bar,weight , dfg) ; dfg_bar = dfg
                !         print *,"d2bi info - average deformation gradient, strain, rot read from binary file "
                        
                     
                !      end if
                      
                     
                !else if (len_trim(orientationFile)/=0) then
                if (len_trim(orientationFile)/=0) then
                
                    !   read average orientation from disk
                    if (rank==0) then
                        print *,""
                        print *,"reading orientation file"
                        print *,"^^^^^^^^^^^^^^^^^^^^^^^^"
                        print *,""
                        print *,"filename """//trim(orientationFile)//""""
                        open(unit=809,file=trim(orientationFile),action="read")
                            read(unit=809,fmt='(i10)') kk           !   number of grains
                            print *,"d2bi info - read number of grains = ",kk
                            allocate(nGrain(0:kk))
                            nGrain = 0
                            do ii = 1,nAtoms
                                ll = nint( x(4,ii) ) ; ll = max(0,min(kk,ll))
                                nGrain( ll ) = nGrain( ll ) + 1
                            end do
                            print *,"d2bi info - atom count in each grain"
                            do ii = 1,kk
                                print *,"    grain ",ii," count ",nGrain(ii)
                            end do
                            ll = maxloc(nGrain,dim=1) - 1
                            if (rank==0) print *,"using orientation for most populous grain ",ll," frac ", nGrain(ll)*100.0/nAtoms," %"
                            deallocate(nGrain)
                            do ii = 1,kk
                                read(unit=809,fmt=*) ix,dfg 
                                if (ii==ll) then
                                    dfg_bar = dfg
                                    print *,"d2bi info - average deformation gradient, strain, rot read from orientation file "                            
                                end if
                                
                            end do
                            
                        close(unit=809)
                    end if                 
                     
                else if (.not. unitDefGrad) then
                
                    !   reading deformation gradient per atom, using this to construct field
                    if (rank==0) then
                        print *,""
                        print *,"deformation gradient per atom"
                        print *,"^^^^^^^^^^^^^^^^^^^^^^^^^^^^^"
                        print *,""   
                    
                        kk = maxval(nint( x(4,:) ))
                        print *,"d2bi info - read number of grains from .xyz file = ",kk
                        allocate(nGrain(0:kk))
                        nGrain = 0
                        do ii = 1,nAtoms
                            ll = nint( x(4,ii) ) ; ll = max(0,min(kk,ll))
                            nGrain( ll ) = nGrain( ll ) + 1
                        end do
                        
                        print *,"findVoids.exe info - atom count in each grain"
                        do ii = 1,kk
                            print *,"    grain ",ii," count ",nGrain(ii)
                        end do
        
                        ll = maxloc(nGrain,dim=1) - 1
                        print *,"using average deformation gradient for most populous grain ",ll," frac ", nGrain(ll)*100.0/nAtoms," %"
                        deallocate(nGrain)
            
                        
                        call computeAvgDefGrad(eps_bar,dfg_bar,weight)
                        do ii = 1,nAtoms
                            if ( nint(x(4,ii)) == ll) call computeAvgDefGrad( x(5:13,ii),1.0d0, eps_bar,dfg_bar,weight)
                        end do
                        call computeAvgDefGrad(eps_bar,dfg_bar,weight , dfg) ; dfg_bar = dfg
                        print *,"d2bi info - average deformation gradient, strain, rot read from extended .xyz file "
                    
                    end if
                      
                         
                    
                else
                        
                    !   assuming unit deformation gradient - note that columns were set above.
                     
                    dfg_bar = reshape( (/1,0,0,0,1,0,0,0,1/) , (/3,3/) )
                     
                 
                end if    
                
                call MPI_BCAST( dfg_bar, 9, MPI_DOUBLE_PRECISION , 0, MPI_COMM_WORLD, ierror )  
                call DefGradToStrainAndRotMat(dfg_bar,eps_bar,rot_bar)
                
                
                if (rank==0) then
                    write (*,fmt='(a36,a,a36,a,a36)') "average def grad","    ","average strain","    ","average rotation"
                    write (*,fmt='(3f12.6,a,3f12.6,a,3f12.6)') dfg_bar(1,:),"    ",eps_bar(1,:),"    ",rot_bar(1,:)
                    write (*,fmt='(3f12.6,a,3f12.6,a,3f12.6)') dfg_bar(2,:),"    ",eps_bar(2,:),"    ",rot_bar(2,:)
                    write (*,fmt='(3f12.6,a,3f12.6,a,3f12.6)') dfg_bar(3,:),"    ",eps_bar(3,:),"    ",rot_bar(3,:)
                end if            
                
                
    
                if (alterg) a_cell = matmul( dfg_bar,a_cell )
                if (rank==0) then
                    print *,"d2bi info - conventional unit cell expected including rotation/strain"
                    write (*,fmt='(3f12.6)') a_cell(1,:)
                    write (*,fmt='(3f12.6)') a_cell(2,:)
                    write (*,fmt='(3f12.6)') a_cell(3,:)
                    print *,""
                end if
                
            end if      !   computation of average def grad for alter_g calculation
            call pause(tt(T_DEFGRAD))
            
            
         
     
            
            
            
            
     
            
        !---    D2BI calculation        
    
            
            tt(T_D2BI_INIT) = Callipers_ctor()
            
            if (rank==0) then
                print *,""
                print *,"find dynamical 2 beam image"
                print *,"^^^^^^^^^^^^^^^^^^^^^^^^^^^"
                print *,""
            end if
            
            
            do nDiff = 1,nDiffConditions
                if (rank==0) then
                    print *,""
                    print *,"diffraction condition ",nDiff,"/",nDiffConditions
                    print *,""
                end if
                
                
                !print *,"ready rank ",rank," a"
            
                if (miller4) then
                    g4(1:4) = g4array(1:4,nDiff)
                    k4(1:4) = k4array(1:4,nDiff)
                    g3(1:3) = MillerBravaisToMiller_plane(g4)    
                    k3(1:3) = MillerBravaisToMiller_direction(k4)
                    if (rank==0) then 
                        write(*,fmt='(a,4f8.2)') "Miller-Bravais g = ",g4
                        write(*,fmt='(a,4f8.2)') "Miller-Bravais k = ",k4
                    end if
                else            
                    g3(1:3) = g4array(1:3,nDiff)
                    k3(1:3) = k4array(1:3,nDiff)
                end if            
                if (rank==0) then 
                    write(*,fmt='(a,4f8.2)') "Miller index   g = ",g3
                    write(*,fmt='(a,4f8.2)') "Miller index   k = ",k3
                end if
                sg = sgarray(nDiff)
                ng = ngarray(nDiff)
                dark = darkarray(nDiff)
                
                
                if (sg /= LIB_CLA_NODEFAULT_R) then
                    modg = (2*3.141592654d0/a0)*norm2(g3)
                    if (abs(2*3.141592654d0*Sg)>modg) call errorExit("d2bi error - input parameter s_g/|g| must be < 1/(2pi) ( 1 - |g|/(2|k|) )")               
                end if
                
        
               if (len_trim(extinction_distances_filename)/=0) then
                    idat(1:3,1) = (/ 0,0,0 /)
                    idat(1:3,2) = nint( g3(1:3) )
                    call readExtinctionDistance( extinction_distances_filename,idat(1:3,1:2),getSymmetry(latt),dat(1:2), ok )
                    if (ok) then
                        xi0 = dat(1)
                        xig = dat(2)
                    else
                        stop "d2bi error - failed to find xi_g for [000] and g=[hkl] in extinction distances file"
                    end if
                end if
                
            
            
                
                d(1:3) = d0(1:3)
                if (rank==0) print *,""
                if (rank==0) print *,"d2bi info - finding extent of atoms in cell"
                call rotateToStandardFrame( k3,a_cell, g3 ,rot_bar )
                if (rank==0) then
                    write (*,fmt='(a)') "d2bi info - rotation to bring g vector into 1- direction and zone axis into 3- direction"
                    write (*,fmt='(3f12.6)') rot_bar(1,:)
                    write (*,fmt='(3f12.6)') rot_bar(2,:)
                    write (*,fmt='(3f12.6)') rot_bar(3,:)
                end if
                xmin = huge(1.0) ; xmax = -huge(1.0)
                do ii = 1,nAtoms
                    yy(1:3) = x(1:3,ii)
                    yy(1:3) = rot_bar(1:3,1)*yy(1) + rot_bar(1:3,2)*yy(2) + rot_bar(1:3,3)*yy(3)
                    xmin(1:3) = min( xmin(1:3),yy(1:3) )
                    xmax(1:3) = max( xmax(1:3),yy(1:3) )
                end do
                if (rank==0) then
                    write (*,fmt='(a)') "d2bi info - min/max/extent ( bounding rectangle ) of atoms after rotation to zone axis"                
                    write (*,fmt='(a,3f12.4)') "x: ",xmin(1),xmax(1),xmax(1)-xmin(1)
                    write (*,fmt='(a,3f12.4)') "y: ",xmin(2),xmax(2),xmax(2)-xmin(2)
                    write (*,fmt='(a,3f12.4)') "z: ",xmin(3),xmax(3),xmax(3)-xmin(3)
                end if
                do ii = 1,3
                    if ( d0(ii)==LIB_CLA_NODEFAULT_R ) d(ii) = xmax(ii) - xmin(ii)
                end do
                 
                d2bi = DynamicalTwoBeamImaging_ctor( super,latt,a0 , k3,g3,V, m(1),m(2))
                call setConventionalUnitCell( d2bi,a_cell )            
                call setVoxelSupercellSize( d2bi,d )
                 
          
                if (rank==0) then
                    print *,""
                    print *,"d2bi info - establishing diffraction conditions for imaging"
                 end if
        
                if (sg /= LIB_CLA_NODEFAULT_R) then
                    if (rank==0) print *,"d2bi info - setting deviation parameter and finding crystal tilt"
                    call setDeviationParameter(d2bi,sg)
                    call setSupercellTilt(d2bi)
                
                else if (ng /= LIB_CLA_NODEFAULT_R) then
                    if (rank==0) print *,"d2bi info - setting (g,ng) condition and finding crystal rotation and tilt"
                    if (surf) then
                        call setSupercellTilt( d2bi,ng, extinction_distances_filename,x )          !   finds rotated frame voxels, using actual extent of atoms to generate true thickness
                    else
                        call setSupercellTilt( d2bi,ng, extinction_distances_filename )               !   finds rotated frame voxels
                    end if 
                else if (R(1,1) /= LIB_CLA_NODEFAULT_R) then
                    if (U(1,1) /= LIB_CLA_NODEFAULT_R) then
                        if (rank==0) print *,"d2bi info - setting diffraction condition using previously determined crystal rotation and tilt"
                        call setSupercellTilt( d2bi, U,R )                              !   finds rotated frame voxels
                    else
                        if (rank==0) print *,"d2bi info - setting diffraction condition using no crystal rotation, but previously determined tilt"
                        U = reshape( (/ 1,0,0,0,1,0,0,0,1 /),(/3,3/) )
                        call setSupercellTilt( d2bi, U,R )                              !   finds rotated frame voxels
                    end if                
                else
                    if (rank==0) print *,"d2bi info - setting deviation parameter sg=0 and finding crystal tilt"
                    call setDeviationParameter(d2bi,sg=0.0d0)
                    call setSupercellTilt(d2bi)       
                end if
     
        
            !---    the diffraction conditions are now set. 
                call getRotation(d2bi,U0,R0) 
                sg = getDeviationParameter( d2bi )              !   could have been set by rotations
                call setExtinctionDistances( d2bi,xi0,xig )                                             !   sets values for xi_0 and xi_g
                call findImageParameters( d2bi,aperpx,xmax,modg,Mx,My,Mz )                              !   returns angstroms per pixel and number of rotated frame voxel subdivisions
                call setuseDensityField(d2bi,useDensity)
                
         
            !---    output to screen information about the deviation parameter and the crystal and foil tilts used
                if (rank==0) then
                    print *,"d2bi info - deviation parameter ",sg," (A)"
                    
                    print *,"d2bi info - crystal rotation matrix , foil tilt matrix"
                    write (*,fmt='(3f12.5,a,3f12.5)') U0(1,:),"    ",R0(1,:) 
                    write (*,fmt='(3f12.5,a,3f12.5)') U0(2,:),"    ",R0(2,:) 
                    write (*,fmt='(3f12.5,a,3f12.5)') U0(3,:),"    ",R0(3,:) 
                                                      
                    write(*,fmt='(a)',advance="no") "d2bi info - crystal rot -U "
                    write(bbbb,fmt='(f16.12)') U0(1,1) ; write(*,fmt='(a)',advance="no")  trim(adjustl(bbbb))//","     
                    write(bbbb,fmt='(f16.12)') U0(2,1) ; write(*,fmt='(a)',advance="no")  trim(adjustl(bbbb))//","
                    write(bbbb,fmt='(f16.12)') U0(3,1) ; write(*,fmt='(a)',advance="no")  trim(adjustl(bbbb))//","
                    write(bbbb,fmt='(f16.12)') U0(1,2) ; write(*,fmt='(a)',advance="no")  trim(adjustl(bbbb))//","
                    write(bbbb,fmt='(f16.12)') U0(2,2) ; write(*,fmt='(a)',advance="no")  trim(adjustl(bbbb))//","
                    write(bbbb,fmt='(f16.12)') U0(3,2) ; write(*,fmt='(a)',advance="no")  trim(adjustl(bbbb))//","
                    write(bbbb,fmt='(f16.12)') U0(1,3) ; write(*,fmt='(a)',advance="no")  trim(adjustl(bbbb))//","
                    write(bbbb,fmt='(f16.12)') U0(2,3) ; write(*,fmt='(a)',advance="no")  trim(adjustl(bbbb))//","
                    write(bbbb,fmt='(f16.12)') U0(3,3) ; write(*,fmt='(a)',advance="yes") trim(adjustl(bbbb)) 
                                                                                
                    write(*,fmt='(a)',advance="no") "d2bi info - foil tilt -R "
                    write(bbbb,fmt='(f16.12)') R0(1,1) ; write(*,fmt='(a)',advance="no")  trim(adjustl(bbbb))//","
                    write(bbbb,fmt='(f16.12)') R0(2,1) ; write(*,fmt='(a)',advance="no")  trim(adjustl(bbbb))//","
                    write(bbbb,fmt='(f16.12)') R0(3,1) ; write(*,fmt='(a)',advance="no")  trim(adjustl(bbbb))//","
                    write(bbbb,fmt='(f16.12)') R0(1,2) ; write(*,fmt='(a)',advance="no")  trim(adjustl(bbbb))//","
                    write(bbbb,fmt='(f16.12)') R0(2,2) ; write(*,fmt='(a)',advance="no")  trim(adjustl(bbbb))//","
                    write(bbbb,fmt='(f16.12)') R0(3,2) ; write(*,fmt='(a)',advance="no")  trim(adjustl(bbbb))//","
                    write(bbbb,fmt='(f16.12)') R0(1,3) ; write(*,fmt='(a)',advance="no")  trim(adjustl(bbbb))//","
                    write(bbbb,fmt='(f16.12)') R0(2,3) ; write(*,fmt='(a)',advance="no")  trim(adjustl(bbbb))//","
                    write(bbbb,fmt='(f16.12)') R0(3,3) ; write(*,fmt='(a)',advance="yes") trim(adjustl(bbbb)) 
                end if
                                    
                
                
            !---    construct an output file name which reflects the diffraction conditions
                if (len_trim(outfile0)>0) then
                    if (miller4) then
                        write(bbbb,fmt='(i2)') nint(g4(1)) ; bbbb = adjustl(bbbb)
                        write(aaaa,fmt='(i2)') nint(g4(2)) ; bbbb = trim(bbbb)//"_"//adjustl(aaaa)
                        write(aaaa,fmt='(i2)') nint(g4(3)) ; bbbb = trim(bbbb)//"_"//adjustl(aaaa)            
                        write(aaaa,fmt='(i2)') nint(g4(4)) ; bbbb = trim(bbbb)//"_"//adjustl(aaaa)            
                    else
                        write(bbbb,fmt='(i2)') nint(g3(1)) ; bbbb = adjustl(bbbb)
                        write(aaaa,fmt='(i2)') nint(g3(2)) ; bbbb = trim(bbbb)//"_"//adjustl(aaaa)
                        write(aaaa,fmt='(i2)') nint(g3(3)) ; bbbb = trim(bbbb)//"_"//adjustl(aaaa)
                    end if                
                    !print *,"d2bi info - writing simulated TEM image"
                    dummy = trim(outfile0)//".g"//trim(bbbb)
            
                    if (miller4) then
                        write(bbbb,fmt='(f5.2)') k4(1) ; bbbb = adjustl(bbbb)
                        write(aaaa,fmt='(f5.2)') k4(2) ; bbbb = trim(bbbb)//"_"//adjustl(aaaa)
                        write(aaaa,fmt='(f5.2)') k4(3) ; bbbb = trim(bbbb)//"_"//adjustl(aaaa)                
                        write(aaaa,fmt='(f5.2)') k4(4) ; bbbb = trim(bbbb)//"_"//adjustl(aaaa)                
                    else
                        write(bbbb,fmt='(f5.2)') k3(1) ; bbbb = adjustl(bbbb)
                        write(aaaa,fmt='(f5.2)') k3(2) ; bbbb = trim(bbbb)//"_"//adjustl(aaaa)
                        write(aaaa,fmt='(f5.2)') k3(3) ; bbbb = trim(bbbb)//"_"//adjustl(aaaa)
                    end if
                    dummy = trim(dummy)//".k"//trim(bbbb)
            
                    if ( ng /= LIB_CLA_NODEFAULT_R ) then
                        write(cccc,fmt='(f8.2)') ng   ; bbbb = "_ng"//adjustl(cccc)
                    else if (dark) then                                     
                        write(cccc,fmt='(f8.5)') sg   ; bbbb = "_df_sg"//adjustl(cccc)
                    else
                        write(cccc,fmt='(f8.5)') sg   ; bbbb = "_bf_sg"//adjustl(cccc)                
                    end if
                    outfile = trim(dummy)//trim(bbbb) 
                    if (nPrecAngle > 1) outfile = trim(outfile)//"_p"
                    !if (nAperture > 1)  outfile = trim(outfile)//"_a"
                end if        
                        
                
                
                allocate(d2bimg(0:Mx-1,0:My-1)) 
                call pause(tt(T_D2BI_INIT))
                tt(T_D2BI) = Callipers_ctor()
                call pause(tt(T_D2BI))
                
                
                do nn = 0,ntomoAngle-1
                
                 
                    
                    if (ntomoAngle == 1) then           
                                               
                        if (opxyz) call opXyzFile( d2bi, x, opxyzfilename = trim(outfile)//".xyz")
                    
                    else
                    
                    !---    construct the rotation matrix for this tomography step.
                        !theta = -tomoAngle + 2*tomoAngle*nn/(ntomoAngle-1)
                        
                        theta = tomoAngle - 2*tomoAngle*nn/(ntomoAngle-1)
                        if (rank==0) then
                            print *,""
                            print *,"d2bi info - tomography step ",(nn+1),"/",ntomoAngle," theta = ",theta," deg"
                        end if
                        theta = theta * 3.141592654d0 / 180.0d0
                        R = reshape( (/ 1.0d0,0.0d0,0.0d0 , 0.0d0,cos(theta),sin(theta) , 0.0d0,-sin(theta),cos(theta) /),(/ 3,3 /) )
                        R = matmul( R,R0 )
                        
                        if (rank==0) then
                            print *,"d2bi info - rotation to standard frame, rotation matrix tilt"    
                            write (*,fmt='(3f12.5,a,3f12.5)') U0(1,:),"    ",R(1,:)
                            write (*,fmt='(3f12.5,a,3f12.5)') U0(2,:),"    ",R(2,:)
                            write (*,fmt='(3f12.5,a,3f12.5)') U0(3,:),"    ",R(3,:)
                            print *,""
                        end if                
                        call setSupercellTilt( d2bi, U0,R ) 
                                        
                        if (opxyz) call opXyzFile( d2bi, x, opxyzfilename = numberFile( trim(outfile),nn,"xyz" ))                
                
                    end if                
                    
                    !print *,"ready rank ",rank," b"
                    
                    
                !---    compute the simulated TEM image            
                    call MPI_BARRIER(MPI_COMM_WORLD,ierror)
                    call start(tt(T_D2BI))
                                
                    call computeD2BI_img( d2bi,x, d2bimg , precAngle,nPrecAngle,apertureAngle,nAperture  )           
        
                    if (.not. dark)  d2bimg = 1-d2bimg
                      
                    if (rank==0) then     
                        print *,""            
                        call report(d2bi)
                        print *,""     
                        if (.not. dark) then
                            write(*,fmt='(a,3f12.6,a)') "d2bi info - minmaxavg d2bimg     ",minval(d2bimg),maxval(d2bimg),sum(d2bimg)/(size(d2bimg,dim=1)*size(d2bimg,dim=2))," (transmitted beam) "
                        else                    
                            write(*,fmt='(a,3f12.6,a)') "d2bi info - minmaxavg d2bimg     ",minval(d2bimg),maxval(d2bimg),sum(d2bimg)/(size(d2bimg,dim=1)*size(d2bimg,dim=2))," (diffracted beam) "
                        end if
                        print *,""     
                    end if                        
                    
                    call pause(tt(T_D2BI) )
                     
                    
        
                 
                    if (rank==0) then
                        if (len_trim(outfile0)>0) then
                            call start(tt(T_OUTPUT))
                            
                            print *,""
                            if (ntomoAngle == 1) then
                                dummy = trim(outfile) 
                            else
                                dummy = numberFile( trim(outfile),nn,"" )
                            end if
                            write(*,fmt='(a)') "writing output to """//trim(dummy)//".dat"""
                            write(*,fmt='(a)') repeat("^",len_trim(dummy)+24)
                           ! print *,""
                            
                            open( unit=903,file=trim(dummy)//".tmp",action="write" )
                                write(unit=903,fmt='(a)') "# d2bi v."//trim(VERSION)
                                if (miller4) then                        
                                    write(unit=903,fmt='(a,4f10.3)') "# k       ",k4
                                    write(unit=903,fmt='(a,4f10.3,f12.5)') "# g       ",g4,modg
                                else
                                    write(unit=903,fmt='(a,3f10.3)') "# k       ",k3
                                    write(unit=903,fmt='(a,3f10.3,f12.5)') "# g       ",g3,modg
                                end if
                                write(unit=903,fmt='(a,3i8)')    "# M       ",M
                                write(unit=903,fmt='(a,3f10.0)') "# V (keV) ",V
                                write(unit=903,fmt='(a,3f12.5)') "# xi (A)  ",xi0,xig
                                if (ng /= LIB_CLA_NODEFAULT_R) then
                                    write(unit=903,fmt='(a,f8.3)') "# n_g     ",ng 
                                end if
                                write(unit=903,fmt='(a,3f12.5)') "# s_g (1/A)    ",getDeviationParameter( d2bi )
                              
                                write(unit=903,fmt='(a,3f12.3)') "# img space (A)",xmax
                                write(unit=903,fmt='(a,3f12.5)') "# scale (A/px) ",aperpx
                                
                                
                                call getRotation(d2bi,U,R) 
                                write(unit=903,fmt='(a,9f12.6)') "# rotation ",matmul( R,U )
                                
                                
                                write(unit=903,fmt='(2i8)') size(d2bimg,dim=1),size(d2bimg,dim=2)
                                write(unit=903,fmt='(1000f12.7)') d2bimg
                            close(unit=903)
                            call system( "mv "//trim(dummy)//".tmp "//trim(dummy)//".dat " )
                            print *,""
                            
                            if (opPng) then
                                call outputImgAsPng( d2bimg, png_min,png_max,png_blur/aperpx, dummy )
                                print *,""
                            end if
                            
                            call pause(tt(T_OUTPUT))
                        end if
                    end if            
                    
                    if (rank==0) print *,"" 
                    call MPI_BARRIER(MPI_COMM_WORLD,ierror)
                    
                end do      !   n tomography steps            
                deallocate(d2bimg)
               
                
            end do      !   nDiff diffraction conditions     
            
            
            
            if (rank==0) print *,""            
    
            
            
    
            if (rank==0) then
                print *,""
                print *,"timing data (s)"
                print *,"^^^^^^^^^^^^^^^"
                print *,""
                write(*,fmt='(a,f16.6)') "   read file ",elapsed(tt(T_READ))
                write(*,fmt='(a,f16.6)') "   ctor      ",elapsed(tt(T_CTOR))
                write(*,fmt='(a,f16.6)') "   def grad  ",elapsed(tt(T_DEFGRAD))
                write(*,fmt='(a,f16.6)') "   D2BI init ",elapsed(tt(T_D2BI_INIT))
                write(*,fmt='(a,f16.6)') "   D2BI calc ",elapsed(tt(T_D2BI))
                write(*,fmt='(a,f16.6)') "   output    ",elapsed(tt(T_OUTPUT))
                write(*,fmt='(a,f16.6)') "   total     ",elapsed(tt(T_TOTAL))
                print *,""
        
        
                print *,""
                
            end if
            call errorExit("done")
            
            
        contains
    !---^^^^^^^^
    
            subroutine outputImgAsPng( d2bimg, iso_min,iso_max, blur, filename )
        !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
                real(kind=real64),dimension(:,:),intent(in)         ::      d2bimg
                real(kind=real64),intent(in)                        ::      iso_max,iso_min
                real(kind=real64),intent(in)                        ::      blur
                character(len=*),intent(in)                         ::      filename
                real(kind=real64),dimension(size(d2bimg,dim=1),size(d2bimg,dim=2))      ::      img
                
                
                write(*,fmt='(a)') "writing output to """//trim(filename)//".png"""
                write(*,fmt='(a)') repeat("^",len_trim(filename)+24)
                if (blur>0) then
                    write(*,fmt='(a,f8.2,a)') "d2bi::outputImgAsPng info - applying gaussian blur = ",blur," px"
                    call gaussianBlurImage( d2bimg,blur, img )
                else
                    img = d2bimg
                end if
                if ( abs(iso_max-iso_min) == 0 ) then
                    !   normalise output
                    print *,"d2bi::outputImgAsPng info - normalise intensity"  
                    call writePng( trim(filename)//".png",img,normalise=.true. )
                else
                    !   contrast scale output
                    img = (img - iso_min)/(iso_max-iso_min)
                    img = max(0.0d0,min(1.0d0,img))
                    write(*,fmt='(3(a,f16.8))') "d2bi::outputImgAsPng info - intensity range ",iso_min,":",iso_max
                    call writePng( trim(filename)//".png",img,normalise=.false. )                
                end if
                return
            end subroutine outputImgAsPng
                
    
    
            
             subroutine gaussianBlurImage( f_in,t , f_out )
        !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
        !*      produce a Gaussian blur of input image, width t
                real(kind=real64),dimension(0:,0:),intent(in)       ::      f_in
                real(kind=real64),dimension(0:,0:),intent(inout)    ::      f_out
                real(kind=real64),intent(in)                        ::      t
                
                
                
                real(kind=real64),dimension(:),allocatable      ::      kernel ,f_stripe,w_stripe
                  
                
                integer             ::      Nx,Ny
                integer             ::      ix,iy,jj,kk
                integer             ::      Nk
                real(kind=real64)   ::      i2s2,ww,wf,ws!,wx
                
                Nx = size(f_in,dim=1)
                Ny = size(f_in,dim=2)
                 
            !---    compute the unnormalised kernel
                Nk = max(5,ceiling( t*5 ))             !   range of pixels to search is +/- Nk
                allocate(kernel( 0:Nk ))
                i2s2 = 1/(2*t*t)             
                kernel(0) = 1.0d0   
                do ix = 1,Nk                
                    kernel( ix )  = exp( -ix*ix*i2s2 )                                
                end do
                
                     
            !---    compute the output blurred image. First do x strips.
                allocate(f_stripe(0:Ny-1))
                allocate(w_stripe(0:Nx-1))
             
                do ix = 0,Nx-1
                    wf = 0.0d0 ; ws = 0.0d0
                    do jj = max(0,ix-Nk),min(Nx-1,ix+Nk)
                        kk = abs(jj-ix)           
                        ww = kernel( kk )                                         
                        wf = wf + ww*f_in(jj,0)
                        ws = ws + ww
                    end do
                    f_out(ix,0) = wf
                    w_stripe(ix) = 1/ws
                end do
                do iy = 1,Ny-1
                    do ix = 0,Nx-1
                        wf = 0.0d0 ; ws = 0.0d0
                        do jj = max(0,ix-Nk),min(Nx-1,ix+Nk)
                            kk = abs(jj-ix)           
                            ww = kernel( kk )                                         
                            wf = wf + ww*f_in(jj,iy)
                            ws = ws + ww
                        end do
                        f_out(ix,iy) = wf                     
                    end do
                end do
     
            !   at this point, f_out has blurring in the x-direction only. weight stores the kernel weighting from this op.            
                
                                 
            !---    now do y strips
                do ix = 0,Nx-1
                !   make a copy of this stripe
                    f_stripe(0:Ny-1) = f_out(ix,0:Ny-1)
                    do iy = 0,Ny-1
                        wf = 0.0d0 ; ws = 0.0d0
                        do jj = max(0,iy-Nk),min(Ny-1,iy+Nk)
                            kk = abs(jj-iy)                                  
                            ww = kernel( kk )
                            wf = wf + ww*f_stripe(jj)
                            ws = ws + ww
                        end do
                        f_out(ix,iy) = wf *w_stripe(ix) /(ws)                !   note: ws /= 0 because there is always at least one pixel contributing. Several really.
                                            
                    end do
                end do
            !---    now have a normalised Gaussian blur function
                          
                
                return
            end subroutine gaussianBlurImage        
            
    
            subroutine errorExit(message)
        !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
                character(len=*),intent(in)             ::      message
                integer                                 ::      ierror
                if (rank==0) print *,trim(message)
                call MPI_FINALIZE(ierror)
                stop
            end subroutine errorExit
                
            
        end program run_d2bi