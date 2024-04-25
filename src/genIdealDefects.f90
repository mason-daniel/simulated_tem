! gfortran ${MYF90LIB}/NBAX_StringTokenizers.f90 ${MYF90LIB}/Lib_XYZFiles.f90 src/Lib_Perlin.f90 src/genIdealDefects.f90  -o bin/genIdealDefects.exe -ffree-line-length-256 -O0 -fbacktrace -fbounds-check -g -DDEBUG
! ./bin/genIdealDefects.exe


    program generateIdealDefects
!---^^^^^^^^^^^^^^^^^^^^^^^^^^^^
!*      rewrite of a standalone program to generate idealised defects in a lattice
        use NBAX_StringTokenizers
        use Lib_XYZFiles
        use Lib_RandomSeed
        use Lib_Perlin2dPlane
        use iso_fortran_env
        implicit none

        integer,parameter       ::      NLATTICE        = 4
        integer,parameter       ::      LATTICE_SC      = 0
        integer,parameter       ::      LATTICE_BCC     = 1
        integer,parameter       ::      LATTICE_BCC111  = 2
        integer,parameter       ::      LATTICE_FCC     = 3
        character(len=6),dimension(0:NLATTICE-1),parameter       ::      LATTICE_NAME = (/ "sc    ","bcc   ","bcc111","fcc   " /)


        integer,parameter       ::      DEFECT_NDATA    = 22     !   maximum number of real data entries to describe defect = habit plane, normal, compression etc.

        integer,parameter       ::      DEFECT_TYPE_SPHERE      = 0
        integer,parameter       ::      DEFECT_TYPE_LOOP        = 1
        integer,parameter       ::      DEFECT_TYPE_HALF        = 2
        integer,parameter       ::      DEFECT_TYPE_RSS         = 3
        integer,parameter       ::      DEFECT_TYPE_PERLIN2D    = 4
        integer,parameter       ::      DEFECT_TYPE_ANNULUS     = 5
        character(len=8),dimension(0:5),parameter               ::      DEFECT_TYPE_NAME = (/ "sphere  ","loop    ","half    ","rss     ","Perlin2D","annulus " /)

        integer,parameter       ::      DEFECT_CHARACTER_NORMAL = 0
        integer,parameter       ::      DEFECT_CHARACTER_VAC    = 1
        integer,parameter       ::      DEFECT_CHARACTER_INT    = 2
        integer,parameter       ::      DEFECT_CHARACTER_2INT   = 3
        integer,parameter       ::      DEFECT_CHARACTER_QHULL  = 4
        character(len=6),dimension(0:4),parameter               ::      DEFECT_CHAR_NAME = (/ "normal","vac   ","int   ","2 int ","qhull " /)
        
        real(kind=real64),parameter                         ::      DIPOLE_OFFSET = 1/4.0d0


    !---    the cell
        integer                                             ::      Nx,Ny,Nz            !   cell repeats
        logical                                             ::      px,py,pz            !   periodicity in x,y,z
        real(kind=real64)                                   ::      a0                  !   lattice parameter
        integer                                             ::      nMotif              !   number of motif points
        real(kind=real64),dimension(3,8)                    ::      motif               !   motif points
        real(kind=real64),dimension(3,3)                    ::      cell                !   cell lattice parameters
        real(kind=real64),dimension(3,3)                    ::      icell               !   inverse cell lattice parameters
        real(kind=real64),dimension(3,3)                    ::      strain              !   cell strain parameters
        integer                                             ::      latticeType
        
        
        integer                                             ::      wsnMotif            !   number of motif points for W-S cell
        real(kind=real64),dimension(3,24)                   ::      wsmotif             !   W-S cell motif points
        logical                                             ::      perfectLattice
        character(len=256)                                  ::      inputFile

    !---    the atoms
        integer                                             ::      nAtoms,nAtoms0,nAtoms1     !   number of atoms, in perfect lattice
        real(kind=real64),dimension(:,:),allocatable        ::      x0,x1                   !   positions, in perfect lattice
        integer,dimension(:),allocatable                    ::      t0,t1                   !   atom types
        integer                                             ::      atomCharacter

    !---    the defects
        integer                                             ::      m                   !   number of defects
        integer,dimension(:),allocatable                    ::      defect_type         !   type of defect
        integer,dimension(:),allocatable                    ::      defect_char         !   character of defect
        real(kind=real64),dimension(:,:),allocatable        ::      defect_centre
        real(kind=real64),dimension(:,:),allocatable        ::      defect_data

        integer                                             ::      nV,nI
        integer,dimension(:),allocatable                    ::      defect_nPointDefect
        integer,dimension(:),allocatable                    ::      x0InDefect


    !---    the output file
        integer                                             ::      nAtomNames          !   number of different types of atoms registered
        character(len=8),dimension(1000)                    ::      atomNames           !   = "Fe","C" etc
        character(len=256)                                  ::      prefix
        type(XYZFile)                                       ::      xyz
        character(len=XYZFILE_COMMENTLINELENGTH),dimension(:),allocatable      ::     headerLine 
        real(kind=real64),dimension(:,:),pointer            ::      colp
        
    !---
        integer                 ::      ii,jj,kk,ll
        character(len=256)      ::      dummy
        integer                 ::      ix,iy,iz,ik
        real(kind=real64),dimension(3)      ::      xx,yy,zz,dz
        real(kind=real64),dimension(4)      ::      qq
        real(kind=real64),dimension(100)    ::      distanceToMotifPoint
        real(kind=real64)       ::      dd,ee
        logical                 ::      ok

        print *,"generateIdealDefects"
        print *,"^^^^^^^^^^^^^^^^^^^^"
        print *,""
        

    !---    establish output file
        print *,"enter filename prefix"
        read(*,fmt='(a)') prefix
        print *,"generateIdealDefects info - prefix = """//trim(prefix)//""""
        print *,""

        print *,"enter atom names"
        read(*,fmt='(a)') dummy
        call parse(dummy,atomNames,nAtomNames)
        write(*,fmt='(a)',advance="no") " generateIdealDefects info - registered atom names "
        do ii = 1,nAtomNames-1
            write(*,fmt='(a7)',advance="no") atomNames(ii)//","
        end do
        write(*,fmt='(a)',advance="yes")  atomNames(nAtomNames)
        print *,""




    !---    establish lattice
        write(*,fmt='(a)',advance="no") " enter lattice name ("
        do ii = 0,NLATTICE-2
            write(*,fmt='(a7)',advance="no") LATTICE_NAME(ii)//","
        end do
        write(*,fmt='(a6,a)',advance="yes")  LATTICE_NAME(NLATTICE-1),") or input file name"
        read(*,fmt='(a)') dummy
        call setMotif( dummy,perfectLattice )
        if (.not. perfectLattice) then
            inputFile = dummy
            write(*,fmt='(a)',advance="no") " now enter lattice name"
            read(*,fmt='(a)') dummy
            call setMotif( dummy,ok )
        end if
            
        
    
    !---    find cell repeats
        print *,"enter cell repeats Nx,Ny,Nz"
        if (latticeType == LATTICE_BCC111) then
            print *,"note (5x3x8) is approximately cuboidal for bcc111 lattice"
        end if
        read(*,fmt=*) Nx,Ny,Nz
        print *,"generateIdealDefects info - Nx,Ny,Nz = ",Nx,Ny,Nz

        print *,"enter cell periodicity px,py,pz"
        read(*,fmt=*) px,py,pz
        print *,"generateIdealDefects info - px,py,pz = ",px,py,pz
        
        
    !---    construct perfect lattice in cell coords
        nAtoms0 = Nx*Ny*Nz*nMotif
        allocate(x0(3,nAtoms0))
        allocate(t0(nAtoms0))
        print *,"generateIdealDefects info - nAtoms (perfect lattice)   = ",nAtoms0
        do iz = 0,Nz-1
            do iy = 0,Ny-1
                do ix = 0,Nx-1
                    do ik = 1,nMotif
                        ii = ik + nMotif*(ix + Nx*(iy + Ny*iz))
                        x0(1:3,ii) = (/ ix,iy,iz /) + motif(1:3,ik)
                        t0(ii) = 1
                    end do
                end do
            end do
        end do
        print *,""            
            
        call init_random_seed()

    !---    read in lattice parameter
        print *,"enter lattice parameter a0"
        read(*,fmt=*) a0
        print *,"enter strain matrix"
        read(*,fmt=*) strain(1:3,1:3)
         
        
        
        cell(1:3,1) = a0*(/1,0,0/)
        cell(1:3,2) = a0*(/0,1,0/)
        cell(1:3,3) = a0*(/0,0,1/)
        
        call inverse3Mat( cell,icell )
        
        xx(1:3) = (/1,0,0/)
        yy(1:3) = (/0,1,0/)
        zz(1:3) = (/0,0,1/)

        if (latticeType == LATTICE_BCC111) then
            xx(1:3) = xx(1:3)*sqrt(2.0d0)
            yy(1:3) = yy(1:3)*sqrt(6.0d0)
            zz(1:3) = zz(1:3)*sqrt(0.75d0)
        end if
        print *,"generateIdealDefects info - lattice parameters"
        write(*,fmt='(a4,3f12.5)') "a",cellToReal( xx  )
        write(*,fmt='(a4,3f12.5)') "b",cellToReal( yy  )
        write(*,fmt='(a4,3f12.5)') "c",cellToReal( zz  )
        print *,""
        
        cell(1:3,1) = a0*xx(1:3)
        cell(1:3,2) = a0*yy(1:3)
        cell(1:3,3) = a0*zz(1:3)
        
        
        
        
        call inverse3Mat( cell,icell )
        
        
        if (perfectLattice) then                    
            xyz = XYZFile_ctor( prefix )
            call setAtomNames( xyz,atomNames(1:nAtomNames) )
            call setNColumns( xyz,4 )
            call setColumn_Description( xyz,"atom  position_x position_y position_z  charge" )
            nAtoms1 = nAtoms0
            allocate(x1(3,nAtoms1))
            allocate(t1(nAtoms1))
            x1(:,:) = x0(:,:)
            t1(:) = t0(:)
            call setColumns( xyz,x1 )
            call setAtomTypes( xyz,t1 )
            call report(xyz)
        else
            xyz = XYZFile_ctor( inputFile )     
            call readHeader(xyz,ok)  
            if (.not. ok) then
                print *,"generateIdealDefects error - could not read .xyz file """//trim(inputFile)//""""
                stop
            end if 
            call input(xyz,ok)
            call report(xyz)
            nAtoms1 = getNatoms(xyz)
            allocate(x1(3,nAtoms1))
            allocate(t1(nAtoms1))
            call getColumns( xyz,x1 )
            do ii = 1,nAtoms1
                xx(1:3) = icell(1:3,1)*x1(1,ii) + icell(1:3,2)*x1(2,ii) + icell(1:3,3)*x1(3,ii)
                x1(1:3,ii) = xx(1:3)
            end do
            call setColumns( xyz,x1 )
            
            call getAtomTypes( xyz,t1 )
            ii = getnAtomNames(xyz)
            if (ii /= nAtomNames) then
                print *,"generateIdealDefects error - incompatible number of atom types: input file says ",nAtomNames," .xyz file says ",ii
                stop
            end if      
            call setFilename(xyz,prefix)
        end if
        

        
    !---    read in number of defects
        print *,"enter number of defects M"
        read(*,fmt=*) m
        print *,"generateIdealDefects info - M = ",m
        allocate(defect_type(0:m))
        allocate(defect_char(0:m))
        allocate(defect_centre(3,0:m))
        allocate(defect_data(DEFECT_NDATA,0:m))
        allocate(defect_nPointDefect(0:m))
        allocate(x0InDefect(m))
        defect_nPointDefect = 0
        print *,""
        do ii = 1,m
            print *,"defect ",ii,"/",m
            print *,"enter character of defect ( "//DEFECT_CHAR_NAME(1)//","//DEFECT_CHAR_NAME(2)//","//DEFECT_CHAR_NAME(4)//")"
            read(*,fmt='(a)') dummy
            defect_char(ii) = getCharacter( dummy )
            write (*,fmt='(a)',advance="no") " enter type of defect      ( "
            do jj = 0,4
                write (*,fmt='(a)',advance="no") DEFECT_TYPE_NAME(jj)//","
            end do
            write (*,fmt='(a)',advance="yes") DEFECT_TYPE_NAME(5)//")"
!            print *,"enter type of defect      ( "//DEFECT_TYPE_NAME(0)//",",DEFECT_TYPE_NAME(1)//",",DEFECT_TYPE_NAME(2)//",",DEFECT_TYPE_NAME(3)//",",DEFECT_TYPE_NAME(4),")"
            read(*,fmt='(a)') dummy
            defect_type(ii) = getType( dummy )

            if ( (defect_type(ii) == DEFECT_TYPE_LOOP).or.(defect_type(ii) == DEFECT_TYPE_HALF).or.(defect_type(ii) == DEFECT_TYPE_PERLIN2D).or.(defect_type(ii) == DEFECT_TYPE_ANNULUS) ) then
                print *,"enter burgers vector of defect ( cell coords )"
                read(*,fmt=*) xx
                call setBurgers( ii,xx )
            end if
            if ( (defect_type(ii) == DEFECT_TYPE_LOOP).or.(defect_type(ii) == DEFECT_TYPE_ANNULUS) ) then
                print *,"enter habit plane normal of defect "
                read(*,fmt=*) xx
                call setNormal( ii,xx )                        
            else if ( (defect_type(ii) == DEFECT_TYPE_HALF).or.(defect_type(ii) == DEFECT_TYPE_PERLIN2D) ) then
                print *,"enter second direction "
                read(*,fmt=*) xx
                call setNormal( ii,xx )      
            end if  
            if ( (defect_type(ii) == DEFECT_TYPE_LOOP).or.(defect_type(ii) == DEFECT_TYPE_HALF).or.(defect_type(ii) == DEFECT_TYPE_ANNULUS) ) then                                   
                print *,"enter displacement scaling of defect "
                read(*,fmt=*) dd
                call setScaling( ii,dd )
            else if ( defect_type(ii) == DEFECT_TYPE_PERLIN2D ) then                                   
                print *,"enter threshold level -1:1 "
                read(*,fmt=*) dd
                call setScaling( ii,dd )
            end if
            if ( defect_type(ii) == DEFECT_TYPE_RSS ) then
                print *,"enter rss count"
                read(*,fmt=*) dd
                call setRssCount( ii,nint(dd) )
            else 
                print *,"enter centre of defect ( cell coords - note motif )"
                read(*,fmt=*) xx
                defect_centre(1:3,ii) = xx
                if ( defect_type(ii) == DEFECT_TYPE_PERLIN2D ) then                   
                    print *,"enter lengthscale ( cell coords )"
                else
                    print *,"enter radius of defect ( cell coords )"
                end if
                read(*,fmt=*) dd
                call setRadius( ii,dd )
            end if
            if ( defect_type(ii) == DEFECT_TYPE_ANNULUS ) then
                print *,"enter second radius"
                read(*,fmt=*) dd
                call setSecondRadius( ii,dd )
            end if
            
            call setCaptureScaling( ii,1.0d0 ) 
            call setLoopBasis( ii )

            
            
            
        !---    try to get the correct number of point defects by considering a trial in a perfect lattice            
            if ( defect_type(ii) == DEFECT_TYPE_RSS ) then
            
            else if ( (defect_type(ii) == DEFECT_TYPE_LOOP).or.(defect_type(ii) == DEFECT_TYPE_ANNULUS) ) then
            !---    construct a trial loop defect
                print *,""
                print *,"considering trial defect with burgers vector and normal aligned"
                defect_type(0) = defect_type(ii)
                defect_char(0) = defect_char(ii)
                defect_data(:,0) = defect_data(:,ii)
                call setNormal( 0,getBurgers(ii) )     
                call setLoopBasis( 0 )            
                call setCaptureScaling( 0,1.0d0 ) 
                call planeOffset(0,zz,toggle=.false.)        
                defect_centre(1:3,0) = (/Nx/2,Ny/2,Nz/2/) + zz
                
                    
            !---    count number of point defects in trial
                kk = 0
                do jj = 1,nAtoms0                
                    if (inDefect( x0(1:3,jj),0 )) kk = kk + 1                    
                end do
                defect_nPointDefect(0) = kk
                call reportDefect(6,0)            
                print *,""
                
            !---    now try again but with correct normal
                call setNormal( 0,getNormal(ii) )             
                call setLoopBasis( 0 )        
                
                dd = 1.0d0
                do ll = 1,12                
                !---    count number of point defects in trial
                    call setCaptureScaling( 0,dd )                
                    kk = 0
                    do jj = 1,nAtoms1                
                        if (inDefect( x1(1:3,jj),0 )) kk = kk + 1                    
                    end do
                    dd = (dd + defect_nPointDefect(0)*1.0d0/max(0.5,real(kk)))/2
                    if (defect_nPointDefect(0) == kk) then
                        print *,"point defect count with desired normal ",kk
                        call reportDefect(6,0)    
                        exit
                    end if
                    print *,"point defect count with desired normal ",kk," scaling ",dd
                    if (ll == 6) then
                        print *,"genIdealDefects warning - could not fix correct number of point defects with this normal"
                        
                        ee = huge(1.0)
                        do jj = 1,nAtoms1       
                        
                            dz(1:3) = x1(1:3,jj) - defect_centre(1:3,ii)                                      !   separation in cell coordinates
                            dz = applypbc(dz)
                            dz = cellToReal(dz)         
                            if (norm2(dz)<ee) then
                                ee = norm2(dz)
                                kk = jj
                            end if
                        end do 
                        print *,"nearest atom to centre ",kk, x1(1:3,kk) , ee
                        
                        
                        
                        defect_centre(1:3,0) = defect_centre(1:3,0) - zz
                        call planeOffset(0,zz,toggle=.true.)        
                        defect_centre(1:3,0) = defect_centre(1:3,0) + zz
                    end if                
                end do
                print *,""
                defect_data(:,ii) = defect_data(:,0)
                
                
                
            !---    set desired defect
                do jj = 1,nMotif
                    xx = defect_centre(1:3,ii) - ( int( defect_centre(1:3,ii) ) + motif(1:3,jj) )
                    distanceToMotifPoint(jj) = norm2(xx)
                    print *,"distanceToMotifPoint ",jj,xx,distanceToMotifPoint(jj)
                end do
                jj = minloc( distanceToMotifPoint(1:nMotif) , dim = 1 )            
                print *,"best motif ",jj
                print *,int( defect_centre(1:3,ii) )+ motif(1:3,jj)
                print *,zz
                defect_centre(1:3,ii) = ( int( defect_centre(1:3,ii) ) + motif(1:3,jj) ) + zz
            end if    
            
        end do
        
        
        
        

    !---    count actual number of atoms in defects
        nV = 0 ; nI = 0
        defect_nPointDefect = 0
        do jj = 1,nAtoms1
            kk = 1      !   number of atoms on this site
            do ii = 1,m
                if (inDefect( x1(1:3,jj),ii )) then
                    if (defect_char(ii) == DEFECT_CHARACTER_VAC) then
                        kk = kk - 1
                    else if (defect_char(ii) == DEFECT_CHARACTER_INT) then
                        kk = kk + 1
                    else if (defect_char(ii) == DEFECT_CHARACTER_QHULL) then
                        kk = kk + 1
                    end if
                    defect_nPointDefect(ii) = defect_nPointDefect(ii) + 1
                else if (inRss(jj,ii)) then
                    if (defect_char(ii) == DEFECT_CHARACTER_VAC) then
                        kk = kk - 1
                    end if
                    defect_nPointDefect(ii) = defect_nPointDefect(ii) + 1
                else if (inPerlin2d(x1(1:3,jj),ii)) then
                    if (defect_char(ii) == DEFECT_CHARACTER_VAC) then
                        kk = kk - 1
                    else if (defect_char(ii) == DEFECT_CHARACTER_INT) then
                        kk = kk + 1
                    end if
                    defect_nPointDefect(ii) = defect_nPointDefect(ii) + 1
                end if                    
            end do
            
            if (any(defect_char(1:m) == DEFECT_CHARACTER_QHULL)) then
                if (kk > 1) then
                    nI = nI + 1
                end if
            else                 
                if (kk <= 0) then
                    nV = nV + 1
                else if (kk == 2) then
                    nI = nI + 1
                else if (kk == 3) then
                    nI = nI + 2
                else if (kk > 3) then
                    print *,"generateIdealDefects warning - attempt to place 3+ interstitials into same site - reduce to 2."
                    nI = nI + 3
                end if
            end if
            
        end do
        do ii = 1,m
            print *,"defect ",ii,"/",m
            call reportDefect(6,ii)
            print *,""
        end do
        if (any(defect_char(1:m) == DEFECT_CHARACTER_QHULL)) then
            print *,"total point defect count ",nI
            nAtoms = (1+wsnMotif)*nI
        else
            print *,"total point defect count nV = ",nV," nI = ",nI
            nAtoms = nAtoms1 + nI - nV
        end if


    !---    reconstruct lattice in cell coords

        
        call setNAtoms( xyz,nAtoms )
        print *,"generateIdealDefects info - nAtoms (including defects)   = ",nAtoms
!        kk = 0      !   number of header lines
!        allocate(headerLine(m*8))
!        do ii = 1,m
!            call reportDefectHeaderLines( ii,headerLine(kk+1:) )
!            do ll = 1,8
!                dummy = headerLine(kk+1)
!                if (dummy(1:1) .eq. "#") then
!                    kk = kk + 1
!                else
!                    exit
!                end if
!            end do
!        end do
!        call setnHeaderLines(xyz,kk)
!        call setHeaderLines(xyz,headerLine(1:kk))
!        call setHeaderLines(xyz,Nx,Ny,Nz,cell)        
        
        nAtoms = 0
        do jj = 1,nAtoms1
            atomCharacter = DEFECT_CHARACTER_NORMAL
            kk = 1      !   number of atoms on this site
            x0inDefect = 0
            ll = 0      !   number of defects atom j belongs to
            dz = 0.0d0
            do ii = 1,m
                if (inDefect( x1(1:3,jj),ii )) then
                    if (defect_char(ii) == DEFECT_CHARACTER_VAC) then
                        kk = kk - 1
                    else if (defect_char(ii) == DEFECT_CHARACTER_INT) then
                        kk = kk + 1
                    else if (defect_char(ii) == DEFECT_CHARACTER_QHULL) then
                        kk = kk + 1
                    end if
                    ll = ll + 1
                    x0inDefect(ll) = ii
                else if (inRss( jj,ii )) then
                    if (defect_char(ii) == DEFECT_CHARACTER_VAC) then
                        kk = kk - 1
                    end if
                    ll = ll + 1
                    x0inDefect(ll) = ii                    
                else if (inPerlin2d( x1(1:3,jj),ii )) then
                    if (defect_char(ii) == DEFECT_CHARACTER_VAC) then
                        kk = kk - 1
                    else if (defect_char(ii) == DEFECT_CHARACTER_INT) then
                        kk = kk + 1
                    end if
                    ll = ll + 1
                    x0inDefect(ll) = ii                    
                else
                    dz = dz + deltaz( x1(1:3,jj),ii )
                end if

            end do
            if (any(defect_char(1:m) == DEFECT_CHARACTER_QHULL)) then
                if (kk <= 1) then
                    atomCharacter = DEFECT_CHARACTER_VAC
                else                    
                    atomCharacter = DEFECT_CHARACTER_QHULL
                end if
            else                 
                if (kk <= 0) then
                    atomCharacter = DEFECT_CHARACTER_VAC
                else if (kk == 2) then
                    atomCharacter = DEFECT_CHARACTER_INT
                else if (kk >= 3) then
                    atomCharacter = DEFECT_CHARACTER_2INT
                end if
            end if

            select case(atomCharacter)
                case(DEFECT_CHARACTER_VAC)
                    !   no atom to place
                case(DEFECT_CHARACTER_QHULL)
                    nAtoms = nAtoms + 1
                    qq(1:3) = cellToReal( x1(1:3,jj)  )  
                    qq(4)   = 1.0d0
                    call setColumns( xyz,nAtoms,qq )
                    call setAtomType( xyz,nAtoms,t1(jj) )
                    do ii = 1,wsnMotif
                        nAtoms = nAtoms + 1
                        qq(1:3) = cellToReal( x1(1:3,jj)+wsMotif(1:3,ii)  )  
                        qq(4)   = 0.0d0
                        call setColumns( xyz,nAtoms,qq )
                        call setAtomType( xyz,nAtoms,t1(jj) )
                    end do
                case(DEFECT_CHARACTER_NORMAL)
                    nAtoms = nAtoms + 1
                    qq(1:3) = cellToReal( x1(1:3,jj)  ) + dz(1:3)
                    qq(4)   = norm2(dz)/a0
                    call setColumns( xyz,nAtoms,qq )
                    call setAtomType( xyz,nAtoms,t1(jj) )
                case(DEFECT_CHARACTER_INT)
                    zz = getBurgers( x0inDefect(1) )
                    qq(4)   = 1.0d0
                    qq(1:3) = cellToReal( x1(1:3,jj)+zz(1:3)*DIPOLE_OFFSET  ) + dz(1:3)
                    call setColumns( xyz,nAtoms+1,qq )
                    qq(1:3) = cellToReal( x1(1:3,jj)-zz(1:3)*DIPOLE_OFFSET  ) + dz(1:3)
                    call setColumns( xyz,nAtoms+2,qq )
                    call setAtomType( xyz,nAtoms+1,t1(jj) )
                    call setAtomType( xyz,nAtoms+2,t1(jj) )
                    nAtoms = nAtoms + 2
                case(DEFECT_CHARACTER_2INT)
                    zz = getBurgers( x0inDefect(1) )
                    xx = getBurgers( x0inDefect(2) )
                    call completeBasis( zz,xx,yy )
                    qq(4)   = 2.0d0
                    qq(1:3) = cellToReal( x1(1:3,jj)+zz(1:3)*DIPOLE_OFFSET                         ) + dz(1:3)
                    call setColumns( xyz,nAtoms+1,qq )
                    qq(1:3) = cellToReal( x1(1:3,jj)-zz(1:3)*DIPOLE_OFFSET/2-xx(1:3)*DIPOLE_OFFSET*sqrt(0.75d0)  ) + dz(1:3)
                    call setColumns( xyz,nAtoms+2,qq )
                    qq(1:3) = cellToReal( x1(1:3,jj)-zz(1:3)*DIPOLE_OFFSET/2+xx(1:3)*DIPOLE_OFFSET*sqrt(0.75d0)  ) + dz(1:3)
                    call setColumns( xyz,nAtoms+3,qq )
                    call setAtomType( xyz,nAtoms+1,t1(jj) )
                    call setAtomType( xyz,nAtoms+2,t1(jj) )
                    call setAtomType( xyz,nAtoms+3,t1(jj) )
                    nAtoms = nAtoms + 3
            end select


        end do
        
        
    !---    adjust pbc        
              
        
        call getColumnsp( xyz,colp )
        if (.not. px) then
            do ii = 1,nAtoms
                colp(1:3,ii) = colp(1:3,ii) + cell(1:3,1)*Nx/4
            end do
        end if
        if (.not. py) then
            do ii = 1,nAtoms
                colp(1:3,ii) = colp(1:3,ii) + cell(1:3,2)*Ny/4
            end do
        end if
        if (.not. pz) then
            do ii = 1,nAtoms
                colp(1:3,ii) = colp(1:3,ii) + cell(1:3,3)*Nz/4
            end do
        end if
        colp(1:3,1:nAtoms) = colp(1:3,1:nAtoms) + a0/4     

        cell = cell + matmul( strain,cell )
        
        
        if (px) then
            if (py) then
                if (pz) then                     
                    call setColumn_Description(xyz,Nx,Ny,Nz,cell) 
                else
                    call setColumn_Description(xyz,Nx,Ny,Nz+Nz/2,cell) 
                end if
            else
                if (pz) then                     
                    call setColumn_Description(xyz,Nx,Ny+Ny/2,Nz,cell) 
                else
                    call setColumn_Description(xyz,Nx,Ny+Ny/2,Nz+Nz/2,cell) 
                end if            
            end if
        else
            if (py) then
                if (pz) then                     
                    call setColumn_Description(xyz,Nx+Nx/2,Ny,Nz,cell) 
                else
                    call setColumn_Description(xyz,Nx+Nx/2,Ny,Nz+Nz/2,cell) 
                end if
            else
                if (pz) then                     
                    call setColumn_Description(xyz,Nx+Nx/2,Ny+Ny/2,Nz,cell) 
                else
                    call setColumn_Description(xyz,Nx+Nx/2,Ny+Ny/2,Nz+Nz/2,cell) 
                end if            
            end if
        end if   
         
        do ii = 1,nAtoms
            xx = colp(1:3,ii)
            colp(1:3,ii) = xx(1:3) + strain(1:3,1)*xx(1) + strain(1:3,2)*xx(2) + strain(1:3,3)*xx(3) 
        end do
        
        
    !---    output file
        print *,""
        print *,"output file"
        call setnHeaderLines(xyz,0) 
        
        
        
        
        call report(xyz)
        if (any(defect_char(1:m) == DEFECT_CHARACTER_QHULL)) then
            call outputQhull(xyz)            
        else
            print *,"writing file ..."
            call output(xyz)
            print *,"... writing file done"
        end if
        
        

    !---    goodbye
        
!         call delete(xyz)
!         deallocate(defect_type)
!         deallocate(defect_char)
!         deallocate(defect_centre)
!         deallocate(defect_data)
!         deallocate(x0)
!         deallocate(t0)
!         deallocate(x1)
!         deallocate(t1)
!         deallocate(headerLine)
! 
        print *,""
        print *,"done"
        print *,""




    contains
!---^^^^^^^^

    !---    geometry
    
    
        subroutine planeOffset(i,dz,toggle)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ 
            integer,intent(in)                  ::      i
            real(kind=real64),dimension(3),intent(out)  ::      dz
            logical,intent(in)                  ::      toggle
            real(kind=real64),dimension(8,-1:1,-1:1,-1:1)      ::      planeOffsetb,planeOffsetn
            real(kind=real64),dimension(3)      ::      nn,bb 
            integer                             ::      jj
            integer         ::      ix,iy,iz
            real(kind=real64)   ::      dd,b0
        !---    debugging - check number of planes 
            nn = getNormal(i)
            bb = getBurgers(i)
            b0 = norm2(bb) ; bb = bb/b0
            planeOffsetn = huge(1.0)
            planeOffsetb = huge(1.0)
            do iz = -1,1
              do iy = -1,1
                do ix = -1,1
                    do jj = 1,nMotif
                        dd = dot_product( nn,motif(1:3,jj) + (/ix,iy,iz/) ) 
                        if (abs(dd)>1.0d-8) planeOffsetn( jj,ix,iy,iz ) = dd
                        dd = dot_product( bb,motif(1:3,jj) + (/ix,iy,iz/) ) 
                        if (abs(dd)>1.0d-8) planeOffsetb( jj,ix,iy,iz ) = dd
!                        print *,"plane offset ",jj,ix,iy,iz,planeOffsetn( jj,ix,iy,iz ),planeOffsetb( jj,ix,iy,iz ),ceiling( planeOffsetb( jj,ix,iy,iz )/(planeOffsetn( jj,ix,iy,iz )+1.0d-8) )
                    end do               
                end do
              end do
            end do          

            dd = minval( abs(planeOffsetn) ) 
            print *,"smallest plane spacing n ",dd
            b0 = b0
            print *,"burgers vector |b|       ",b0
            jj = floor( b0/(dd-1.0d-8) )
            print *,"normal planes per burgers vector ",jj
            
            if (mod(jj,2)==1) then
                dz = 0.0d0          !   odd number of planes per burgers vector means should align on a lattice point
            else
                dz = dd*bb(1:3)/2   !   offset by half a plane spacing in direction of burgers vector
            end if
            
            if (toggle) then
                dz(3) = dz(3) + 0.5d0
            end if
         
            print *,"plane offset ",dz                      
            return
        end subroutine planeOffset      
                
                
    

        function inDefect( x,i ) result(is)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      returns true if atom j is part of the defect region i

            real(kind=real64),dimension(3),intent(in)   ::      x
            integer,intent(in)      ::      i
            logical                 ::      is



            real(kind=real64),dimension(3)  ::      dx,xx,yy,nn,bb
            real(kind=real64)       ::      rr,zz,xxxx,yyyy,dxxx,dxyy ,zeta ,rho,rho2

            
            is = .false.

            dx(1:3) = x(1:3) - defect_centre(1:3,i)                                       !   separation in cell coordinates
            dx = applypbc(dx)
            dx = cellToReal(dx)
            rho = getRadius(i)*a0
            select case( defect_type(i) )
                case(DEFECT_TYPE_SPHERE)
                    rr = norm2(dx)
                    is = (rr <= rho+1.0d-8)
                case(DEFECT_TYPE_LOOP)

                
                    call getLoopBasis( i , xx,yy,nn,bb )  !   now xx and yy are in the plane of the loop, with length >= rho. nn,bb is normal,burgers in real space

                    rr = norm2( dx - bb*dot_product(dx,bb)/dot_product(bb,bb) )
                    is = (rr <= rho+1.0d-8)
                    
                    
                    xxxx = dot_product( xx,xx )     !   length2 x
                    yyyy = dot_product( yy,yy )     !   length2 y
                    dxxx = dot_product( dx,xx )
                    dxyy = dot_product( dx,yy )

                    !is = (dxxx*dxxx*yyyy*yyyy + dxyy*dxyy*xxxx*xxxx <= xxxx*xxxx*yyyy*yyyy*(1+1.0d-8))       !   equation of an ellipse  Px^2/a^2 + Py^2/b^2 <= 1
                    if (is) then
                         zz = dot_product( nn,dx )      !   projection on normal vector as real space length
                         zz = zz*dot_product(bb,nn)     !   projection of distance from point in plane to point along burgers vector

                     !---   need to take atoms in a layer whose thickness depends on habit plane
                         zeta = dot_product(bb,nn)*0.49d0*getCaptureScaling(i)
                         is = abs(zz)<norm2(bb)*zeta
                    end if
                    
                case(DEFECT_TYPE_ANNULUS)

                    rho2 = getSecondRadius(i)*a0
                    call getLoopBasis( i , xx,yy,nn,bb )  !   now xx and yy are in the plane of the loop, with length >= rho. nn,bb is normal,burgers in real space

                    rr = norm2( dx - bb*dot_product(dx,bb)/dot_product(bb,bb) )
                    is = ( (rr <= rho+1.0d-8).and.(rr >= rho2-1.0d-8) )
                    
                    
                    xxxx = dot_product( xx,xx )     !   length2 x
                    yyyy = dot_product( yy,yy )     !   length2 y
                    dxxx = dot_product( dx,xx )
                    dxyy = dot_product( dx,yy )

                    !is = (dxxx*dxxx*yyyy*yyyy + dxyy*dxyy*xxxx*xxxx <= xxxx*xxxx*yyyy*yyyy*(1+1.0d-8))       !   equation of an ellipse  Px^2/a^2 + Py^2/b^2 <= 1
                    if (is) then
                         zz = dot_product( nn,dx )      !   projection on normal vector as real space length
                         zz = zz*dot_product(bb,nn)     !   projection of distance from point in plane to point along burgers vector

                     !---   need to take atoms in a layer whose thickness depends on habit plane
                         zeta = dot_product(bb,nn)*0.49d0*getCaptureScaling(i)
                         is = abs(zz)<norm2(bb)*zeta
                    end if

                case(DEFECT_TYPE_HALF)
                
                    call getLoopBasis( i , xx,yy,nn,bb )  !   now xx and nn are in the plane. bb is burgers vector
                
                        
                    xxxx = dot_product( xx,xx )     !   length2 x
                    dxxx = dot_product( dx,xx )
                    is = abs(dxxx/xxxx) < rho
                    
                    if (norm2(dx) < a0*0.5) then
                        print *,"inDefect info - getLoopBasis xx = ",xx
                        print *,"inDefect info - getLoopBasis yy = ",yy
                        print *,"inDefect info - getLoopBasis nn = ",nn
                        print *,"inDefect info - getLoopBasis bb = ",bb
                        print *,"inDefect info - getLoopBasis x = ",x
                        print *,"inDefect info - getLoopBasis dx = ",dx
                        print *,"inDefect info - getLoopBasis xxxx,dxxx,rho = ",xxxx,dxxx,rho
                    end if
                    
                    
                    
                    if (is) then
                         zz = dot_product( nn,dx )      !   projection on normal vector as real space length
                         zz = zz*dot_product(bb,nn)     !   projection of distance from point in plane to point along burgers vector

                     !---   need to take atoms in a layer whose thickness depends on habit plane
                         zeta = dot_product(bb,nn)*0.49d0*getCaptureScaling(i)
                         is = abs(zz)<norm2(bb)*zeta
                    end if
                    
                
                    
            end select


            return
        end function inDefect


        function inRss( j,i ) result(is)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ 
    !*      returns true if atom j is part of the defect region i

            integer,intent(in)      ::      j
            integer,intent(in)      ::      i
            logical                 ::      is




            real(kind=real64)       ::     dd

            integer,save                                ::      rss_n
            integer,dimension(:),allocatable,save       ::      rss_loc
            logical,save                                ::      rss_loc_firstCall = .true.
            integer                 ::      jj,kk
            
            is = .false.

            select case( defect_type(i) )
                case(DEFECT_TYPE_RSS)
                    
                    if (rss_loc_firstCall) then
                        rss_n = getRssCount( i )
                        allocate(rss_loc(rss_n))
                        rss_loc = 0
                        do jj = 1,rss_n
                            do
                                call random_number(dd)
                                kk = int( dd*nAtoms1 ) + 1
                                if (any(rss_loc(1:jj)==kk)) cycle
                                rss_loc(jj) = kk
                                exit
                            end do                            
                        end do                        
                        rss_loc_firstCall = .false.                        
                    end if                  
                    
                    is = any(rss_loc(1:rss_n)==j)
                    
            end select


            return
        end function inRss


        function inPerlin2D( x,i ) result(is)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ 
    !*      returns true if atom j is part of the defect region i
            real(kind=real64),dimension(3),intent(in)   ::      x
            
            integer,intent(in)      ::      i
            logical                 ::      is



            real(kind=real64),dimension(3)      ::     dx,xx,yy,nn,bb

            real(kind=real64)       ::     dd,xxxx,dxxx,yyyy,dyyy,zz,zeta

            integer,save                                ::      p2d_n
            real(kind=real64),dimension(:,:),allocatable,save     ::      p2d_f
            logical,save                                ::      p2d_loc_firstCall = .true.
            integer,save                                ::      minx,maxx,miny,maxy
            
            
            real(kind=real64),dimension(2,2)            ::      aa 
            integer                 ::      jj,kk,nFineGrid,ix,iy
            
            real(kind=real64),dimension(3,8)            ::      corner 
            
            is = .false.

            select case( defect_type(i) )
                case(DEFECT_TYPE_PERLIN2D)
                    
                    nFineGrid = max(Nx,max(Ny,Nz))
                
                    if (p2d_loc_firstCall) then
                    
                        corner(1:3,1) = (/0 ,0 ,0 /)
                        corner(1:3,2) = (/Nx,0 ,0 /)
                        corner(1:3,3) = (/0 ,Ny,0 /)
                        corner(1:3,4) = (/Nx,Ny,0 /)
                        corner(1:3,5) = (/0 ,0 ,Nz/) 
                        corner(1:3,6) = (/Nx,0 ,Nz/)
                        corner(1:3,7) = (/0 ,Ny,Nz/) 
                        corner(1:3,8) = (/Nx,Ny,Nz/)
                    
                        minx = huge(1)
                        maxx = -huge(1)
                        miny = huge(1)
                        maxy = -huge(1)
                        call getLoopBasis( i , xx,yy,nn,bb )  !   now xx and nn are in the plane. bb is burgers vector
                        
                        print *,"xx = ",xx
                        print *,"yy = ",yy
                        do jj = 1,8    
                            dx(1:3) = corner(1:3,jj) - defect_centre(1:3,i)                                            
                            !dx = applypbc(dx)
                            dx = cellToReal(dx)
               
                            dxxx = dot_product( dx,xx )/(a0/2)
                            dyyy = dot_product( dx,yy )/(a0/2)
                            
                            ix = floor( dxxx )
                            iy = floor( dyyy )
                            
                            minx = min(minx,ix)
                            maxx = max(maxx,ix)
                            miny = min(miny,iy)
                            maxy = max(maxy,iy)
                            print *,dx,ix,iy
                            
                        end do
                        
                        print *,"minmax ",minx,maxx,miny,maxy
                        
                        allocate(p2d_f(minx:maxx,miny:maxy))
                        aa(1:2,1) = (/ 0.5d0,0.0d0 /)
                        aa(1:2,2) = (/ 0.0d0,0.5d0 /)
                        call perlinNoise( aa,maxx-minx,maxy-miny, getRadius(i)  , p2d_f )
                        p2d_loc_firstCall = .false.
                        

                        print *,""
                        do iy = -miny,maxy-1
                            do ix = -minx,maxx-1
                                if (p2d_f(ix,iy)>getScaling(i) ) then
                                    write(*,fmt='(a2)',advance="no") "XX"
                                else
                                    write(*,fmt='(a2)',advance="no") " "
                                end if
                            end do
                            write(*,fmt='(a2)',advance="yes") ""            
                        end do
                                                
                        
                        
                        
                        
                        
                        
                        
                        
                    end if                        
                        

                    dx(1:3) = x(1:3) - defect_centre(1:3,i)                                       !   separation in cell coordinates
                    !dx = applypbc(dx)
                    dx = cellToReal(dx)
               
                    call getLoopBasis( i , xx,yy,nn,bb )  !   now xx and nn are in the plane. bb is burgers vector
                                               
                    zz = dot_product( nn,dx )      !   projection on normal vector as real space length
                    zz = zz*dot_product(bb,nn)     !   projection of distance from point in plane to point along burgers vector

                !---   need to take atoms in a layer whose thickness depends on habit plane
                    zeta = dot_product(bb,nn)*0.49d0*getCaptureScaling(i)
                    is = abs(zz)<norm2(bb)*zeta
                    
                    

                    if (is) then                            
                    
                        !print *,"  inPerlin2D ",dx,zz,is
                        
                        dxxx = dot_product( dx,xx )/(a0/2)
                        dyyy = dot_product( dx,yy )/(a0/2)
                        
                        ix = floor( dxxx )
                        iy = floor( dyyy )
                        xxxx = dxxx - ix
                        yyyy = dyyy - iy
                        
                        !print *,"ix,xxxx,iy,yyyy ",ix,xxxx,iy,yyyy 
                                     
                                    
                        !print  *,ix,iy,minx,maxx,miny,maxy
                        if ( (ix >= minx).and.( ix+1 <= maxx).and.(iy >= miny).and.( iy+1 <= maxy) ) then
                            
                            zeta = p2d_f( ix  ,iy   )*(1-xxxx)*(1-yyyy)              &
                                 + p2d_f( ix+1,iy   )*(  xxxx)*(1-yyyy)              &
                                 + p2d_f( ix  ,iy+1 )*(1-xxxx)*(  yyyy)              &
                                 + p2d_f( ix+1,iy+1 )*(  xxxx)*(  yyyy) 
                            
                                 
                            is = zeta > getScaling(i)                             
                        else
                            is = .false.
                        end if    
                                                     
                    end if
                    
                        
                        
                    
            end select


            return
        end function inPerlin2D



        function deltaz( x,i ) result(dz)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      returns a displacement in the direction of the burgers vector

            real(kind=real64),dimension(3),intent(in)   ::      x
            integer,intent(in)      ::      i
            real(kind=real64),dimension(3)              ::      dz


            real(kind=real64),dimension(3)  ::      dx,xx,yy,nn,bb
            real(kind=real64)       ::      rr,zz,xxxx,yyyy,dxxx,dxyy ,z0,rho,rho2

            dz = 0.0d0

            dx(1:3) = x(1:3) - defect_centre(1:3,i)                                       !   separation in cell coordinates
            dx = applypbc(dx)
            dx = cellToReal(dx)
            rho = getRadius(i)*a0
            
            select case( defect_type(i) )
                case(DEFECT_TYPE_SPHERE)

                case(DEFECT_TYPE_LOOP)

                    call getLoopBasis( i , xx,yy,nn,bb )  !   now xx and yy are in the plane of the loop, with length >= rho. nn,bb is normal,burgers in real space

                    rr = norm2( dx - bb*dot_product(dx,bb)/dot_product(bb,bb) )/rho
!                     rr = dxxx*dxxx/(xxxx*xxxx) + dxyy*dxyy/(yyyy*yyyy)       !   equation of an ellipse  Px^2/a^2 + Py^2/b^2 <= 1
!                     rr = sqrt(rr)
                    if (rr < 2) then
                        zz = dot_product( nn,dx )           !   projection on normal vector as real space length
                        zz = zz*dot_product(bb,nn)          !   |b| x projection of distance from point in plane to point along burgers vector

                        yy = ( Nx*cell(1:3,1) + Ny*cell(1:3,2) + Nz*cell(1:3,3) )/4      !    vector to half furthest corner

                        z0 = ( abs(nn(1)*yy(1)) + abs(nn(2)*yy(2)) + abs(nn(3)*yy(3)) ) !   projection to half furthest corner
                        z0 = z0*abs(dot_product(bb,nn))          !   |b| x projection of distance from point in plane to point along burgers vector


!                        print *,zz,yy,nn,z0

                        if (zz*zz < z0*z0) then
                            dz(1:3) = 0.5d0*(1 - abs(zz)/z0)*bb(1:3)
                        end if

                        dz(1:3) = dz(1:3)* min(1.0d0,2.0d0-rr)

                        if (zz<0) dz = -dz
                    end if

                case(DEFECT_TYPE_ANNULUS)
                    rho2 = getSecondRadius(i)*a0
                    call getLoopBasis( i , xx,yy,nn,bb )  !   now xx and yy are in the plane of the loop, with length >= rho. nn,bb is normal,burgers in real space

                    rr = norm2( dx - bb*dot_product(dx,bb)/dot_product(bb,bb) )
!                     rr = dxxx*dxxx/(xxxx*xxxx) + dxyy*dxyy/(yyyy*yyyy)       !   equation of an ellipse  Px^2/a^2 + Py^2/b^2 <= 1
!                     rr = sqrt(rr)
                    if ( (rr>=rho2).and.(rr < 2*rho) ) then
                        zz = dot_product( nn,dx )           !   projection on normal vector as real space length
                        zz = zz*dot_product(bb,nn)          !   |b| x projection of distance from point in plane to point along burgers vector

                        yy = ( Nx*cell(1:3,1) + Ny*cell(1:3,2) + Nz*cell(1:3,3) )/4      !    vector to half furthest corner

                        z0 = ( abs(nn(1)*yy(1)) + abs(nn(2)*yy(2)) + abs(nn(3)*yy(3)) ) !   projection to half furthest corner
                        z0 = z0*abs(dot_product(bb,nn))          !   |b| x projection of distance from point in plane to point along burgers vector


!                        print *,zz,yy,nn,z0

                        if (zz*zz < z0*z0) then
                            dz(1:3) = 0.5d0*(1 - abs(zz)/z0)*bb(1:3)
                        end if

                        dz(1:3) = dz(1:3)* min(1.0d0,2.0d0-rr/rho)

                        if (zz<0) dz = -dz
                    end if


                case(DEFECT_TYPE_HALF)                               
                
                
            end select

            !dz(1:3) = dz(1:3) + strain(1:3,1)*dz(1) + strain(1:3,2)*dz(2) + strain(1:3,3)*dz(3)


            select case(defect_char(i))
                case(DEFECT_CHARACTER_VAC)
                    dz = -dz
            end select

            dz = dz * getScaling(i)

            return
        end function deltaz





        pure function cellToReal(x) result(y)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      convert cell space coordinates to real space coordinates
    !*      optinally apply final strain
            real(kind=real64),dimension(3),intent(in)       ::      x
            !logical,intent(in),optional                     ::      strained
            real(kind=real64),dimension(3)                  ::      y

            y(1:3) = cell(1:3,1)*x(1) + cell(1:3,2)*x(2) + cell(1:3,3)*x(3)
!            if (present(strained)) then
!                if (strained) then
!                    y(1:3) = y(1:3) + strain(1:3,1)*y(1) + strain(1:3,2)*y(2) + strain(1:3,3)*y(3)
!                end if
!            end if
            return
        end function cellToReal

        pure function applyPBC(dx) result(dy)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      apply periodic boundaries to cell space coordinates separation dx
            real(kind=real64),dimension(3),intent(in)       ::      dx
            real(kind=real64),dimension(3)                  ::      dy

            dy = dx
            if (px) then
                if (dy(1)*2>Nx) then
                    dy(1) = dy(1) - Nx
                else if (dy(1)*2<-Nx) then
                    dy(1) = dy(1) + Nx
                end if
            end if

            if (py) then
                if (dy(2)*2>Ny) then
                    dy(2) = dy(2) - Ny
                else if (dy(2)*2<-Ny) then
                    dy(2) = dy(2) + Ny
                end if
            end if

            if (pz) then
                if (dy(3)*2>Nz) then
                    dy(3) = dy(3) - Nz
                else if (dy(3)*2<-Nz) then
                    dy(3) = dy(3) + Nz
                end if
            end if
            return
        end function applyPBC


    !---

        subroutine completeBasis( z,x,y )
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      given the z-axis vector, complete the basis to provide a triplet x,y,z
    !*      if on input x is set, then use this as a hint to attempt to place x along this direction
            real(kind=real64),dimension(3),intent(in)       ::      z
            real(kind=real64),dimension(3),intent(inout)    ::      x
            real(kind=real64),dimension(3),intent(out)      ::      y

            real(kind=real64),dimension(3)      ::      zz        !   normalised
            real(kind=real64)                   ::      xxxx,zzzz  !   vector lengths
            real(kind=real64)                   ::      zdotx

        !---    check we have a normalised z as input
            zzzz = norm2(z)
            if (zzzz == 0) then
                !   can't do anything with zero input vector
                x = (/ 1,0,0 /)
                y = (/ 0,1,0 /)
                return
            end if
            zz = z/zzzz

        !---    check for a sensible hint for the x-direction
            xxxx = norm2(x)
            zdotx = dot_product( zz,x )
            if ( (xxxx < 0.001d0 ).or.(abs(zdotx) >= xxxx*0.999d0) ) then
                print *,"completeBasis info - haven't got a good hint ",zz,x
                !   haven't got a good hint. Make random hint                                
                if (abs(zz(3))>max(abs(zz(1)),abs(zz(2)))) then
                    !   z points along 3-axis
                    x(1) = 1.0d0
                    x(3) = 0.0d0
                    if (abs(zz(1))>0) then
                        x(2) = - zz(2)/zz(1)
                    else
                        x(2) = 0.0d0
                    end if
                    zdotx = dot_product( zz,x )                    
                else if (abs(z(2))>max(abs(z(1)),abs(z(3)))) then
                    !   z points along 2-axis
                    x(3) = 1.0d0
                    x(2) = 0.0d0
                    if (abs(zz(3))>0) then
                        x(1) = - zz(1)/zz(3)
                    else
                        x(1) = 0.0d0
                    end if
                    zdotx = dot_product( zz,x )                    
                else                    
                    !   z points along 1-axis
                    x(2) = 1.0d0
                    x(1) = 0.0d0
                    if (abs(zz(2))>0) then
                        x(3) = - zz(3)/zz(2)
                    else
                        x(3) = 0.0d0
                    end if
                    zdotx = dot_product( zz,x )                    
                end if
                print *,"suggest ",x,norm2(x),zdotx
            end if
            

        !---    remove projection of z on x-direction
            x = x - zz*zdotx
            xxxx = norm2(x)
            x = x/xxxx

        !---    construct y-direction
            y = cross_product( z,x )

!                     !---    debug
!                     print *,"genIdealDefects::completeBasis debug"
!                     print *,x(1),y(1),zz(1)
!                     print *,x(2),y(2),zz(2)
!                     print *,x(3),y(3),zz(3)
!                     print *,dot_product(x,x),dot_product(x,y),dot_product(x,zz)
!                     print *,dot_product(y,x),dot_product(y,y),dot_product(y,zz)
!                     print *,dot_product(zz,x),dot_product(zz,y),dot_product(zz,zz)


            return
        end subroutine completeBasis


        subroutine completeLoopBasis( n,b,rho , x,y )
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      given the loop normal n, and burgers vector b
    !*      and loop radius rho, compute the set of basis vectors
    !*          x,y,n
    !*      so that a position on the loop is given by
    !*          p(theta) = x cos(theta) + y sin(theta)
    !*      and | p - (p.b)/(b.b) b | = rho
    !*      ie the projection on the plane b.r = 0 is a circle radius rho.
    !*      on input n should be normalised, b should have non-zero length
    !*      and n.b should not be zero.

            real(kind=real64),dimension(3),intent(in)       ::      n,b
            real(kind=real64),intent(in)                    ::      rho
            real(kind=real64),dimension(3),intent(out)      ::      x,y

            real(kind=real64)       ::      dd,ib2

        !---    attempt to find a vector x normal to both n and b
            x = cross_product(n,b)
            call completeBasis( n,x,y )


        !---    now find the correct length for x and y
        !   require | x - (x.b)/(b.b) b | = rho
            ib2 = 1/( b(1)*b(1) + b(2)*b(2) + b(3)*b(3) )

            dd = x(1)*b(1) + x(2)*b(2) + x(3)*b(3)
            dd = rho/sqrt( 1 - dd*dd*ib2 )
            x(1:3) = x(1:3) * dd

            dd = y(1)*b(1) + y(2)*b(2) + y(3)*b(3)
            dd = rho/sqrt( 1 - dd*dd*ib2 )
            y(1:3) = y(1:3) * dd
!             
!             
!             print *,"completeLoopBasis info - "
!             print *,x
!             print *,y
!             print *,z
!             print *,dot_product(x,x),dot_product(x,y),dot_product(x,n)
!             print *,dot_product(y,x),dot_product(y,y),dot_product(y,n)
!             print *,dot_product(n,x),dot_product(n,y),dot_product(n,n)
!             

            return
        end subroutine completeLoopBasis

        pure function cross_product(x,y) result(z)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      Ronseal. z = x ^ y
            real(kind=real64),dimension(3),intent(in)       ::      x,y
            real(kind=real64),dimension(3)                  ::      z
            z(1) = x(2)*y(3) - x(3)*y(2)
            z(2) = x(3)*y(1) - x(1)*y(3)
            z(3) = x(1)*y(2) - x(2)*y(1)
            return
        end function cross_product



    !---    help interpretting input

        subroutine setMotif( lattice,ok )
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
            character(len=*),intent(in)     ::      lattice
            logical,intent(out)             ::      ok


            real(kind=real64),dimension(3,1),parameter  ::  SC_MOTIF     = reshape( (/0,0,0/),(/3,1/) )
            real(kind=real64),dimension(3,2),parameter  ::  BCC_MOTIF    = reshape( (/0,0,0 , 1,1,1 /)*0.5d0,(/3,2/) )
            real(kind=real64),dimension(3,6),parameter  ::  BCC111_MOTIF = reshape( (/0,0,0 , 3,5,2 , 0,2,2 , 3,3,0 , 0,4,4 , 3,1,4 /)/6.0d0,(/3,6/) )
            real(kind=real64),dimension(3,4),parameter  ::  FCC_MOTIF    = reshape( (/0,0,0 , 1,1,0 , 1,0,1 , 0,1,1/)*0.5d0,(/3,4/) )

            real(kind=real64),dimension(3,24),parameter ::  BCC_WS_MOTIF = reshape( (/ 2,1,0, 2,-1,0, 2,0,1, 2,0,-1,  -2,1,0, -2,-1,0, -2,0,1, -2,0,-1,     &
                                                                                       1,2,0, -1,2,0, 0,2,1, 0,2,-1,  1,-2,0, -1,-2,0, 0,-2,1, 0,-2,-1,     &
                                                                                       1,0,2, -1,0,2, 0,1,2, 0,-1,2,  1,0,-2, -1,0,-2, 0,1,-2, 0,-1,-2      &                                                                                  
                                                                                       /)*0.125d0,(/3,24/) )            
            
            real(kind=real64),dimension(3,14),parameter ::  FCC_WS_MOTIF = reshape( (/1,1,1, 1,1,-1, 1,-1,1, 1,-1,-1, -1,1,1, -1,1,-1, -1,-1,1, -1,-1,-1,   &
                                                                                      2,0,0, -2,0,0, 0,2,0, 0,-2,0, 0,0,2, 0,0,-2   /)*0.125d0,(/3,14/) )

                                                                                      
            real(kind=real64),dimension(3,8),parameter ::  SC_WS_MOTIF = reshape( (/ 1,1,1, 1,1,-1, 1,-1,1, 1,-1,-1 , -1,1,1, -1,1,-1, -1,-1,1, -1,-1,-1    &                                                                                  
                                                                                       /)*0.25d0,(/3,8/) )            
                                                                                      
            integer     ::      ii

            ok = .true.
            if (trim(lattice)==trim(LATTICE_NAME(LATTICE_BCC))) then
                nMotif = size(BCC_MOTIF,dim=2)
                motif(1:3,1:nMotif) = BCC_MOTIF(1:3,1:nMotif)
                latticeType = LATTICE_BCC
                wsnMotif = size(BCC_WS_MOTIF,dim=2)
                wsmotif(1:3,1:wsnMotif) = BCC_WS_MOTIF(1:3,1:wsnMotif)
!                wsnMotif = 0
            else if (trim(lattice)==trim(LATTICE_NAME(LATTICE_BCC111))) then
                nMotif = size(BCC111_MOTIF,dim=2)
                motif(1:3,1:nMotif) = BCC111_MOTIF(1:3,1:nMotif)
                latticeType = LATTICE_BCC111
                print *,"generateIdealDefects::setMotif warning - WS motif not set"
                wsnMotif = 0
                wsmotif = 0
            else if (trim(lattice)==trim(LATTICE_NAME(LATTICE_FCC))) then
                nMotif = size(FCC_MOTIF,dim=2)
                motif(1:3,1:nMotif) = FCC_MOTIF(1:3,1:nMotif)
                latticeType = LATTICE_FCC
                wsnMotif = size(FCC_WS_MOTIF,dim=2)
                wsmotif(1:3,1:wsnMotif) = FCC_WS_MOTIF(1:3,1:wsnMotif)
            else if (trim(lattice)==trim(LATTICE_NAME(LATTICE_SC))) then
                nMotif = size(SC_MOTIF,dim=2)
                motif(1:3,1:nMotif) = SC_MOTIF(1:3,1:nMotif)
                latticeType = LATTICE_SC
                wsnMotif = size(SC_WS_MOTIF,dim=2)
                wsmotif(1:3,1:wsnMotif) = SC_WS_MOTIF(1:3,1:wsnMotif)
            else
                nMotif = size(SC_MOTIF,dim=2)
                motif(1:3,1:nMotif) = SC_MOTIF(1:3,1:nMotif)
                latticeType = LATTICE_SC
                wsnMotif = size(SC_WS_MOTIF,dim=2)
                wsmotif(1:3,1:wsnMotif) = SC_WS_MOTIF(1:3,1:wsnMotif)
                ok = .false.
                !print *,"generateIdealDefects::setMotif error - did not recognise lattice name """//trim(lattice)//""""
                !stop
            end if


            print *,"generateIdealDefects::setMotif info - lattice = "//trim(LATTICE_NAME(latticeType))
            if (ok) then
                print *,"generateIdealDefects::setMotif info - nMotif = ",nMotif
                do ii = 1,nMotif
                    print *,"motif point ",ii, motif(1:3,ii)
                end do
            end if


            return
        end subroutine setMotif

        function getCharacter( defect ) result (c )
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
            character(len=*),intent(in)     ::      defect
            integer            ::      c

            if (trim(defect)==trim(DEFECT_CHAR_NAME(DEFECT_CHARACTER_VAC))) then
                c = DEFECT_CHARACTER_VAC
            else if (trim(defect)==trim(DEFECT_CHAR_NAME(DEFECT_CHARACTER_INT))) then
                c = DEFECT_CHARACTER_INT
            else if (trim(defect)==trim(DEFECT_CHAR_NAME(DEFECT_CHARACTER_QHULL))) then
                c = DEFECT_CHARACTER_QHULL
            else
                print *,"generateIdealDefects::getCharacter error - did not recognise character """//trim(defect)//""""
                stop
            end if

            return
        end function getCharacter


        function getType( defect ) result (t )
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
            character(len=*),intent(in)     ::      defect
            integer             ::      t

            if (trim(defect)==trim(DEFECT_TYPE_NAME(DEFECT_TYPE_SPHERE))) then
                t = DEFECT_TYPE_SPHERE
            else if (trim(defect)==trim(DEFECT_TYPE_NAME(DEFECT_TYPE_LOOP))) then
                t = DEFECT_TYPE_LOOP
            else if (trim(defect)==trim(DEFECT_TYPE_NAME(DEFECT_TYPE_ANNULUS))) then
                t = DEFECT_TYPE_ANNULUS
            else if (trim(defect)==trim(DEFECT_TYPE_NAME(DEFECT_TYPE_HALF))) then
                t = DEFECT_TYPE_HALF
            else if (trim(defect)==trim(DEFECT_TYPE_NAME(DEFECT_TYPE_PERLIN2D))) then
                t = DEFECT_TYPE_PERLIN2D
            else if (trim(defect)==trim(DEFECT_TYPE_NAME(DEFECT_TYPE_RSS))) then
                t = DEFECT_TYPE_RSS
            else
                print *,"generateIdealDefects::getType error - did not recognise type """//trim(defect)//""""
                stop
            end if

            return
        end function getType

    !----   get/set methods

        subroutine setNormal( i,n )
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^
            integer,intent(in)                                          ::      i       !   number of defect
            real(kind=real64),dimension(3),intent(in)                   ::      n
            real(kind=real64)       ::      nn
            select case( defect_type(i) )
                case(DEFECT_TYPE_SPHERE)
                case(DEFECT_TYPE_LOOP)
                    nn = norm2(n(1:3))
                    if (nn==0) then
                        print *,"generateIdealDefects::setNormal warning - normal has zero magnitude?"
                        defect_data(1:3,i) = 0.0d0
                    else
                        defect_data(1:3,i) = n(1:3) / nn
                    end if
                case(DEFECT_TYPE_ANNULUS)
                    nn = norm2(n(1:3))
                    if (nn==0) then
                        print *,"generateIdealDefects::setNormal warning - normal has zero magnitude?"
                        defect_data(1:3,i) = 0.0d0
                    else
                        defect_data(1:3,i) = n(1:3) / nn
                    end if
                case(DEFECT_TYPE_HALF)
                    nn = norm2(n(1:3))
                    if (nn==0) then
                        print *,"generateIdealDefects::setNormal warning - normal has zero magnitude?"
                        defect_data(1:3,i) = 0.0d0
                    else
                        defect_data(1:3,i) = n(1:3) / nn
                    end if
                case(DEFECT_TYPE_PERLIN2D)
                    nn = norm2(n(1:3))
                    if (nn==0) then
                        print *,"generateIdealDefects::setNormal warning - normal has zero magnitude?"
                        defect_data(1:3,i) = 0.0d0
                    else
                        defect_data(1:3,i) = n(1:3) / nn
                    end if
            end select
            return
        end subroutine setNormal

        function getNormal( i ) result(n)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
            integer,intent(in)                                          ::      i       !   number of defect
            real(kind=real64),dimension(3)                              ::      n
            n(1:3) = defect_data(1:3,i)
            return
        end function getNormal


        subroutine setBurgers( i,b )
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^
            integer,intent(in)                                          ::      i       !   number of defect
            real(kind=real64),dimension(3),intent(in)                   ::      b
            select case( defect_type(i) )
                case(DEFECT_TYPE_SPHERE)
                case(DEFECT_TYPE_LOOP)
                    defect_data(4:6,i) = b(1:3)
                case(DEFECT_TYPE_ANNULUS)
                    defect_data(4:6,i) = b(1:3)
                case(DEFECT_TYPE_HALF)
                    defect_data(4:6,i) = b(1:3)
                case(DEFECT_TYPE_PERLIN2D)
                    defect_data(4:6,i) = b(1:3)
            end select
            return
        end subroutine setBurgers

        function getBurgers( i ) result(b)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
            integer,intent(in)                                          ::      i       !   number of defect
            real(kind=real64),dimension(3)                              ::      b
            real(kind=real64)       ::      dd
            
            b(1:3) = defect_data(4:6,i)
            if ( defect_type(i) == DEFECT_TYPE_RSS ) then
                !   random unit vector
                do
                    call random_number(b) ; b = 2*b - 1
                    dd = dot_product(b,b)
                    if (dd*(1-dd)>0) then
                        b = b / sqrt(dd)
                        exit
                    end if
                end do
            end if
            return
        end function getBurgers

        subroutine setScaling( i,s )
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^
            integer,intent(in)                                          ::      i       !   number of defect
            real(kind=real64),intent(in)                                ::      s

            select case( defect_type(i) )
                case(DEFECT_TYPE_SPHERE)
                case(DEFECT_TYPE_LOOP)
                    defect_data(7,i) = s
                case(DEFECT_TYPE_ANNULUS)
                    defect_data(7,i) = s
                case(DEFECT_TYPE_HALF)
                    defect_data(7,i) = s
                case(DEFECT_TYPE_PERLIN2D)
                    defect_data(7,i) = s
            end select
            return
        end subroutine setScaling

        function getScaling( i ) result(s)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
            integer,intent(in)                                          ::      i       !   number of defect
            real(kind=real64)                                           ::      s

            s = defect_data(7,i)
            return
        end function getScaling
        
        subroutine setRadius( i,r )
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^
            integer,intent(in)                                          ::      i       !   number of defect
            real(kind=real64),intent(in)                                ::      r

            select case( defect_type(i) )
                case(DEFECT_TYPE_SPHERE)
                    defect_data(8,i) = r
                case(DEFECT_TYPE_LOOP)
                    defect_data(8,i) = r
                case(DEFECT_TYPE_ANNULUS)
                    defect_data(8,i) = r
                case(DEFECT_TYPE_HALF)
                    defect_data(8,i) = r
                case(DEFECT_TYPE_PERLIN2D)
                    defect_data(8,i) = r
            end select
            return
        end subroutine setRadius

        function getRadius( i ) result(r)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
            integer,intent(in)                                          ::      i       !   number of defect
            real(kind=real64)                                           ::      r

            r = defect_data(8,i)
            return
        end function getRadius

        subroutine setCaptureScaling( i,c )
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
            integer,intent(in)                                          ::      i       !   number of defect
            real(kind=real64),intent(in)                                ::      c

            select case( defect_type(i) )
                case(DEFECT_TYPE_SPHERE)
                    defect_data(21,i) = 1.0d0
                case(DEFECT_TYPE_LOOP)
                    defect_data(21,i) = c
                case(DEFECT_TYPE_ANNULUS)
                    defect_data(21,i) = c
                case(DEFECT_TYPE_HALF)
                    defect_data(21,i) = c
                case(DEFECT_TYPE_PERLIN2D)
                    defect_data(21,i) = c
            end select
            return
        end subroutine setCaptureScaling

        function getCaptureScaling( i ) result(c)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
            integer,intent(in)                                          ::      i       !   number of defect
            real(kind=real64)                                           ::      c
            c = defect_data(21,i)
            return
        end function getCaptureScaling

        subroutine setRssCount( i,c )
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
            integer,intent(in)                                          ::      i       !   number of defect
            integer,intent(in)                                ::      c

            select case( defect_type(i) )
                case(DEFECT_TYPE_RSS)
                    defect_data(21,i) = c
            end select
            return
        end subroutine setRssCount

        function getRssCount( i ) result(c)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
            integer,intent(in)                                          ::      i       !   number of defect
            integer                                                     ::      c
            c = nint( defect_data(21,i) )
            return
        end function getRssCount

        subroutine setSecondRadius( i,r )
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
            integer,intent(in)                                          ::      i       !   number of defect
            real(kind=real64),intent(in)                                ::      r

            select case( defect_type(i) )
                case(DEFECT_TYPE_ANNULUS)
                    defect_data(22,i) = r
            end select
            return
        end subroutine setSecondRadius

        function getSecondRadius( i ) result(r)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
            integer,intent(in)                                          ::      i       !   number of defect
            real(kind=real64)                                           ::      r
            r = defect_data(22,i)
            return
        end function getSecondRadius


        subroutine setLoopBasis( i )
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^
            integer,intent(in)                                          ::      i       !   number of defect
            real(kind=real64),dimension(3)      ::      nn,bb
            select case( defect_type(i) )
                case(DEFECT_TYPE_SPHERE)
                case(DEFECT_TYPE_LOOP)
                    nn = getNormal(i)
                    nn = cellToReal(nn)
                    nn = nn / norm2(nn)

                    bb = getBurgers(i)
                    bb = cellToReal(bb)
                    call completeLoopBasis( nn,bb,getRadius(i)*a0,defect_data(9:11,i),defect_data(12:14,i) )
                    defect_data(15:17,i) = nn(1:3)
                    defect_data(18:20,i) = bb(1:3)
                case(DEFECT_TYPE_ANNULUS)
                    nn = getNormal(i)
                    nn = cellToReal(nn)
                    nn = nn / norm2(nn)

                    bb = getBurgers(i)
                    bb = cellToReal(bb)
                    call completeLoopBasis( nn,bb,getRadius(i)*a0,defect_data(9:11,i),defect_data(12:14,i) )
                    defect_data(15:17,i) = nn(1:3)
                    defect_data(18:20,i) = bb(1:3)

                case(DEFECT_TYPE_HALF)
                    nn = getNormal(i)
                    nn = cellToReal(nn)
                    nn = nn / norm2(nn)
                    bb = getBurgers(i)
                    bb = cellToReal(bb)
                    call completeLoopBasis( nn,bb,1.0d0,defect_data(9:11,i),defect_data(12:14,i) )
                    defect_data(15:17,i) = nn(1:3)
                    defect_data(18:20,i) = bb(1:3)
                    
                case(DEFECT_TYPE_PERLIN2D)
                    nn = getNormal(i)
                    nn = cellToReal(nn)
                    nn = nn / norm2(nn)
                    bb = getBurgers(i)
                    bb = cellToReal(bb)
                    call completeLoopBasis( nn,bb,1.0d0,defect_data(9:11,i),defect_data(12:14,i) )
                    defect_data(9:11,i) = defect_data(9:11,i)/norm2(defect_data(9:11,i))
                    defect_data(12:14,i) = defect_data(12:14,i)/norm2(defect_data(12:14,i))
                    defect_data(15:17,i) = nn(1:3)
                    defect_data(18:20,i) = bb(1:3)
            end select
            return
        end subroutine setLoopBasis



        subroutine getLoopBasis( i,x,y,n,b )
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
            integer,intent(in)                                          ::      i       !   number of defect
            real(kind=real64),dimension(3),intent(out)                  ::      x,y,n,b
            select case( defect_type(i) )
                case(DEFECT_TYPE_SPHERE)
                    x(1:3) = (/1,0,0/)
                    y(1:3) = (/0,1,0/)
                    n(1:3) = (/0,0,1/)
                    b(1:3) = (/0,0,1/)
                case(DEFECT_TYPE_LOOP)
                    x(1:3) = defect_data(9:11,i)
                    y(1:3) = defect_data(12:14,i)
                    n(1:3) = defect_data(15:17,i)
                    b(1:3) = defect_data(18:20,i)
                case(DEFECT_TYPE_ANNULUS)
                    x(1:3) = defect_data(9:11,i)
                    y(1:3) = defect_data(12:14,i)
                    n(1:3) = defect_data(15:17,i)
                    b(1:3) = defect_data(18:20,i)
                case(DEFECT_TYPE_HALF)
                    x(1:3) = defect_data(9:11,i)
                    y(1:3) = defect_data(12:14,i)
                    n(1:3) = defect_data(15:17,i)
                    b(1:3) = defect_data(18:20,i)
                case(DEFECT_TYPE_PERLIN2D)
                    x(1:3) = defect_data(9:11,i)
                    y(1:3) = defect_data(12:14,i)
                    n(1:3) = defect_data(15:17,i)
                    b(1:3) = defect_data(18:20,i)
            end select
            return
        end subroutine getLoopBasis
    !---    utility

        subroutine reportDefect( u,i )
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
            integer,intent(in)                          ::      u
            integer,intent(in)                          ::      i       !   number of defect

            write (unit=u,fmt='(a)') DEFECT_TYPE_NAME(defect_type(i))//"["//DEFECT_CHAR_NAME(defect_char(i))//"]"
            select case( defect_type(i) )
                case(DEFECT_TYPE_SPHERE)
                    write (unit=u,fmt='(a,3f12.5)') "centre (uc)           ",defect_centre(1:3,i)
                    write (unit=u,fmt='(a,i12)')    "nPointDefect          ",defect_nPointDefect(i)
                    write (unit=u,fmt='(a,f12.5)') "radius (uc)           ",getRadius(i)
                case(DEFECT_TYPE_LOOP)
                    write (unit=u,fmt='(a,3f12.5)') "centre (uc)           ",defect_centre(1:3,i)
                    write (unit=u,fmt='(a,i12)')    "nPointDefect          ",defect_nPointDefect(i)
                    write (unit=u,fmt='(a,3f12.5,a,f12.5)') "Burgers vector (uc)   ",getBurgers(i)," |b| ",norm2( getBurgers(i) )
                    write (unit=u,fmt='(a,3f12.5)') "habit plane normal    ",getNormal(i)
                    write (unit=u,fmt='(a,f12.5)')  "displacement scaling  ",getScaling(i)
                    write (unit=u,fmt='(a,f12.5)')  "radius (uc)           ",getRadius(i)
                case(DEFECT_TYPE_ANNULUS)
                    write (unit=u,fmt='(a,3f12.5)') "centre (uc)           ",defect_centre(1:3,i)
                    write (unit=u,fmt='(a,i12)')    "nPointDefect          ",defect_nPointDefect(i)
                    write (unit=u,fmt='(a,3f12.5,a,f12.5)') "Burgers vector (uc)   ",getBurgers(i)," |b| ",norm2( getBurgers(i) )
                    write (unit=u,fmt='(a,3f12.5)') "habit plane normal    ",getNormal(i)
                    write (unit=u,fmt='(a,f12.5)')  "displacement scaling  ",getScaling(i)
                    write (unit=u,fmt='(a,f12.5)')  "outer radius (uc)     ",getRadius(i)
                    write (unit=u,fmt='(a,f12.5)')  "inner radius (uc)     ",getSecondRadius(i)
                case(DEFECT_TYPE_HALF)
                    write (unit=u,fmt='(a,3f12.5)') "centre (uc)           ",defect_centre(1:3,i)
                    write (unit=u,fmt='(a,i12)')    "nPointDefect          ",defect_nPointDefect(i)
                    write (unit=u,fmt='(a,3f12.5)') "Burgers vector (uc)   ",getBurgers(i)
                    write (unit=u,fmt='(a,3f12.5)') "habit plane normal    ",getNormal(i)
                    write (unit=u,fmt='(a,f12.5)')  "displacement scaling  ",getScaling(i)
                case(DEFECT_TYPE_PERLIN2D)
                    write (unit=u,fmt='(a,3f12.5)') "centre (uc)           ",defect_centre(1:3,i)
                    write (unit=u,fmt='(a,i12)')    "nPointDefect          ",defect_nPointDefect(i)
                    write (unit=u,fmt='(a,3f12.5)') "Burgers vector (uc)   ",getBurgers(i)
                    write (unit=u,fmt='(a,3f12.5)') "habit plane normal    ",getNormal(i)
                    write (unit=u,fmt='(a,f12.5)')  "threshold             ",getScaling(i)
                    write (unit=u,fmt='(a,f12.5)')  "lengthscale (uc)      ",getRadius(i)
                case(DEFECT_TYPE_RSS)
                    write (unit=u,fmt='(a,i12)')    "count                 ",getRssCount(i)
            end select
            return
        end subroutine reportDefect


        subroutine reportDefectHeaderLines( i,headerLine )
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
            integer,intent(in)                          ::      i       !   number of defect
            character(len=XYZFILE_COMMENTLINELENGTH),dimension(:),intent(out)       ::     headerLine 

            headerLine(1) = "# "//DEFECT_TYPE_NAME(defect_type(i))//"["//DEFECT_CHAR_NAME(defect_char(i))//"]"
            select case( defect_type(i) )
                case(DEFECT_TYPE_SPHERE)
                    write (headerLine(2),fmt='(a,3f12.5)') "# "//"centre (uc)           ",defect_centre(1:3,i)
                    write (headerLine(3),fmt='(a,i12)')    "# "//"nPointDefect          ",defect_nPointDefect(i)
                    write (headerLine(4),fmt='(a,f12.5)') "# "//"radius (uc)           ",getRadius(i)
                case(DEFECT_TYPE_LOOP)
                    write (headerLine(2),fmt='(a,3f12.5)') "# "//"centre (uc)           ",defect_centre(1:3,i)
                    write (headerLine(3),fmt='(a,i12)')    "# "//"nPointDefect          ",defect_nPointDefect(i)
                    write (headerLine(4),fmt='(a,3f12.5,a,f12.5)') "# "//"Burgers vector (uc)   ",getBurgers(i)," |b| ",norm2( getBurgers(i) )
                    write (headerLine(5),fmt='(a,3f12.5)') "# "//"habit plane normal    ",getNormal(i)
                    write (headerLine(6),fmt='(a,f12.5)')  "# "//"displacement scaling  ",getScaling(i)
                    write (headerLine(7),fmt='(a,f12.5)')  "# "//"radius (uc)           ",getRadius(i)
                case(DEFECT_TYPE_ANNULUS)
                    write (headerLine(2),fmt='(a,3f12.5)') "# "//"centre (uc)           ",defect_centre(1:3,i)
                    write (headerLine(3),fmt='(a,i12)')    "# "//"nPointDefect          ",defect_nPointDefect(i)
                    write (headerLine(4),fmt='(a,3f12.5,a,f12.5)') "# "//"Burgers vector (uc)   ",getBurgers(i)," |b| ",norm2( getBurgers(i) )
                    write (headerLine(5),fmt='(a,3f12.5)') "# "//"habit plane normal    ",getNormal(i)
                    write (headerLine(6),fmt='(a,f12.5)')  "# "//"displacement scaling  ",getScaling(i)
                    write (headerLine(7),fmt='(a,2f12.5)') "# "//"radii (uc)           ",getSecondRadius(i),getRadius(i)
                case(DEFECT_TYPE_HALF)
                    write (headerLine(2),fmt='(a,3f12.5)') "# "//"centre (uc)           ",defect_centre(1:3,i)
                    write (headerLine(3),fmt='(a,i12)')    "# "//"nPointDefect          ",defect_nPointDefect(i)
                    write (headerLine(4),fmt='(a,3f12.5)') "# "//"Burgers vector (uc)   ",getBurgers(i)
                    write (headerLine(5),fmt='(a,3f12.5)') "# "//"habit plane normal    ",getNormal(i)
                    write (headerLine(6),fmt='(a,f12.5)')  "# "//"displacement scaling  ",getScaling(i)
                case(DEFECT_TYPE_PERLIN2D)
                    write (headerLine(2),fmt='(a,3f12.5)') "# "//"centre (uc)           ",defect_centre(1:3,i)
                    write (headerLine(3),fmt='(a,i12)')    "# "//"nPointDefect          ",defect_nPointDefect(i)
                    write (headerLine(4),fmt='(a,3f12.5)') "# "//"Burgers vector (uc)   ",getBurgers(i)
                    write (headerLine(5),fmt='(a,3f12.5)') "# "//"habit plane normal    ",getNormal(i)
                    write (headerLine(6),fmt='(a,f12.5)')  "# "//"threshold             ",getScaling(i)
                    write (headerLine(7),fmt='(a,f12.5)')  "# "//"lengthscale (uc)      ",getRadius(i)
                case(DEFECT_TYPE_RSS)
                    write (headerLine(2),fmt='(a,i12)')    "# "//"count                 ",getRssCount(i)                    
            end select
            return
        end subroutine reportDefectHeaderLines


    !---

        pure subroutine inverse3Mat(M,N)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      monadic operator
    !*      finds the inverse of a general three matrix
            real(kind=real64),dimension(3,3),intent(in)       ::  M
            real(kind=real64),dimension(3,3),intent(out)      ::  N
            real(kind=real64)            ::      idd
            real(kind=real64),dimension(3,3),parameter        ::      &
            IDENTITY3MAT = reshape( (/ 1.0d0,0.0d0,0.0d0 , 0.0d0,1.0d0,0.0d0 , 0.0d0,0.0d0,1.0d0 /) &
                                   ,(/ 3,3 /) )

            idd = determinant3Mat(M)
            if (abs(idd) < tiny(1.0d0)) then
                N = IDENTITY3MAT
                return
            end if
            idd = 1.0/idd

            N(1,1)   = ( M(2,2)*M(3,3) - M(2,3)*M(3,2) ) * idd
            N(2,1)   = ( M(2,3)*M(3,1) - M(2,1)*M(3,3) ) * idd
            N(3,1)   = ( M(2,1)*M(3,2) - M(2,2)*M(3,1) ) * idd

            N(1,2)   = ( M(1,3)*M(3,2) - M(1,2)*M(3,3) ) * idd
            N(2,2)   = ( M(1,1)*M(3,3) - M(1,3)*M(3,1) ) * idd
            N(3,2)   = ( M(1,2)*M(3,1) - M(1,1)*M(3,2) ) * idd

            N(1,3)   = ( M(1,2)*M(2,3) - M(1,3)*M(2,2) ) * idd
            N(2,3)   = ( M(1,3)*M(2,1) - M(1,1)*M(2,3) ) * idd
            N(3,3)   = ( M(1,1)*M(2,2) - M(1,2)*M(2,1) ) * idd

            return
        end subroutine inverse3Mat

        pure function determinant3Mat(M) result(d)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      returns the determinant of M
            real(kind=real64),dimension(3,3),intent(in)      ::      M
            real(kind=real64)                                ::      d
            real(kind=real64),dimension(9)       ::      dd
            dd(1) = M(1,1)*( M(2,2)*M(3,3) - M(2,3)*M(3,2) )
            dd(2) = M(1,2)*( M(2,3)*M(3,1) - M(2,1)*M(3,3) )
            dd(3) = M(1,3)*( M(2,1)*M(3,2) - M(2,2)*M(3,1) )
            dd(4) = M(2,1)*( M(3,2)*M(1,3) - M(3,3)*M(1,2) )
            dd(5) = M(2,2)*( M(3,3)*M(1,1) - M(3,1)*M(1,3) )
            dd(6) = M(2,3)*( M(3,1)*M(1,2) - M(3,2)*M(1,1) )
            dd(7) = M(3,1)*( M(1,2)*M(2,3) - M(1,3)*M(2,2) )
            dd(8) = M(3,2)*( M(1,3)*M(2,1) - M(1,1)*M(2,3) )
            dd(9) = M(3,3)*( M(1,1)*M(2,2) - M(1,2)*M(2,1) )
            d = (1.0d0/3.0d0) * sum(dd)
            return
        end function determinant3Mat
        
    end program generateIdealDefects