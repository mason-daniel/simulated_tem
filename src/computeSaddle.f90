
    program computeSaddle
!---^^^^^^^^^^^^^^^^^^^^^
!
!       version 0.0.2 - includes option to output electron density and embedding energy
!
!
!

        use Lib_CommandLineArguments
        use Lib_XYZFiles
        use simpleMD
        use EAM_PotentialTables
        use Lib_SimpleSupercells
        use Lib_NBAX
        use NBAX_StringTokenizers
        use Lib_Quicksort
        use Lib_Filenames
        use Lib_Callipers
        use iso_fortran_env
        implicit none
        
        
        external                ::      DSYSV 
        

        real(kind=real64),parameter     ::      KB = 0.00008617d0       !   Boltzmann const eV/K
        real(kind=real64),parameter     ::      PI = 3.1415926535898d0
        character(len=8),parameter      ::      VERSION = "0.0.2"       
        type(CommandLineArguments)      ::      cla

    !---    input parameters
        character(len=256)              ::      filename1 = "" ,filename2 = ""         !   input filenames
        character(len=256)              ::      outfile = ""
        character(len=256)              ::      potential = ""
        real(kind=real64)               ::      R0 = 10.0d0             !   radius to compute moving atoms
        integer                         ::      centralAtom = 1
        real(Kind=real64),dimension(3)  ::      deltax = LIB_CLA_NODEFAULT_R
        logical                         ::      relax = .true.
        logical                         ::      elasticConst = .false.
        logical                         ::      dipoleTens = .false.
        real(kind=real64),dimension(3,3)    ::      strain = 0.0d0
        integer                         ::      nReplicas = 7
        real(kind=real64)               ::      disp_scale = LIB_CLA_NODEFAULT_R             !   permitted displacement scale
        logical                         ::      opElectronDensity = .false.
        
    !---    physically meaningful parameters deduced from input file
        real(kind=real64),dimension(3,3)    ::      a_super
        type(XYZFile)                       ::      xyz
        type(EAM_Alloy)                     ::      eam
        real(kind=real64)                   ::      rc
        logical                             ::      opxyz = .true.
        integer,dimension(:),pointer        ::      at      !   atom types
        integer                             ::      nAtomNames
        logical                             ::      secondFile

    !---    analysis tools
        type(Callipers),dimension(10)       ::      tt        
        integer,parameter                   ::      TIMER_ALL = 1
        integer,parameter                   ::      TIMER_CG = 3
        integer,parameter                   ::      TIMER_IO = 5
        integer,parameter                   ::      TIMER_SADDLE = 6
        
    !---    dummy parameter

        logical                             ::      ok
        type(NBAX)                          ::      xml
        character(len=256)                  ::      dummy
        type(MD)                            ::      smd
        real(kind=real64),dimension(:,:),pointer    ::      xx,vv 
        integer                 ::      ii,jj ,nAtoms !,ni,nv , ki,kj
 
        real(kind=real64)                   ::      ee , c11,c12,c44, vol
        real(kind=real64),dimension(9)      ::      dat
        real(kind=real64),dimension(3)      ::      yy
        character(len=16),dimension(10)     ::      atomNames
        real(kind=real64),dimension(6,6)    ::      CC
        real(kind=real64),dimension(3,3)    ::      PP,FF
        real(kind=real64),dimension(6)      ::      ss
        
            real(kind=real64),dimension(66*6)       ::      work
            integer,dimension(6)                    ::      ipiv

!        integer                             ::      nColumns



    !---    important class instances and physical properties

        real(kind=real64)                           ::       e_old,e_new       !   total (excess) potential energy recorded
        real(kind=real64),dimension(:),allocatable  ::      ee_out


        real(kind=real64),dimension(:,:),allocatable    ::      x_old,x_new

      ! real(kind=real64),dimension(3)          ::      dx

        
        
        
    !---    set up timers        
        do ii = 1,size(tt)
            tt(ii) = Callipers_ctor()
            call pause(tt(ii))
        end do
        call start(tt(TIMER_ALL))    
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
    !---    read command line arguments
        cla = CommandLineArguments_ctor(30)
 

        call setProgramDescription( cla, "computeSaddle.exe" )
        call setProgramVersion( cla, VERSION )

        call get( cla,"f1",filename1 ,LIB_CLA_REQUIRED,"       input filename 1" )
        call get( cla,"f2",filename2 ,LIB_CLA_OPTIONAL,"     input filename 2" )
        outfile = removeSuffix(filename1)
        call get( cla,"o",outfile ,LIB_CLA_OPTIONAL,"      output filename" )
         
        call get( cla,"p",potential,LIB_CLA_REQUIRED,"        potential filename" )

        call get( cla,"e",opElectronDensity ,LIB_CLA_OPTIONAL,"      output electron density and embedding energy data" )
        opxyz = opElectronDensity
        call get( cla,"opxyz",opxyz ,LIB_CLA_OPTIONAL,"  output relaxed configs " )
        call get( cla,"R0",R0,LIB_CLA_OPTIONAL,"     relaxation radius" )
        
        call get( cla,"d",disp_scale,LIB_CLA_OPTIONAL,"      maximum displacement ( default to rc/4 )" )
        
        ii = 3 ; call get( cla,"deltax",deltax,ii,LIB_CLA_OPTIONAL," move central atom by this displacement to generate after state" )
        
        secondFile = hasArgument(cla,"f2")
        if (secondFile .or. hasArgument(cla,"deltax")) then
            call get( cla,"atom",centralAtom,LIB_CLA_REQUIRED,"   central moving atom" )
        else
            call get( cla,"atom",centralAtom,LIB_CLA_OPTIONAL,"   central moving atom" )
        end if
        call get( cla,"CG",relax,LIB_CLA_OPTIONAL,"     CG relax before/after states first" )
        elasticConst = relax
        dipoleTens = relax
        call get( cla,"elastic",elasticConst,LIB_CLA_OPTIONAL,"compute elastic constants" )
        call get( cla,"dipole",dipoleTens,LIB_CLA_OPTIONAL," compute dipole tensor" )
        
        call get( cla,"n",nReplicas,LIB_CLA_OPTIONAL,"      number of NEB replicas" )
        call get( cla,"dbg",SIMPLEMD_DBG,LIB_CLA_OPTIONAL,"    debug" )
        ii=0; dat = LIB_CLA_NODEFAULT_R
        call get( cla,"strain",dat,ii,LIB_CLA_OPTIONAL," add strain ( 1 component hydrostatic, 3 = iso, 9 = aniso )" )
        if (hasArgument(cla,"strain")) then 
            if (ii==1) then
                strain = reshape( (/ dat(1),0.0d0,0.0d0,0.0d0,dat(1),0.0d0,0.0d0,0.0d0,dat(1) /),(/3,3/) )
            else if (ii==3) then
                strain = reshape( (/ dat(1),0.0d0,0.0d0,0.0d0,dat(2),0.0d0,0.0d0,0.0d0,dat(3) /),(/3,3/) )
            else if (ii==9) then
                strain = reshape( dat(1:9),(/3,3/) )
            else
                stop "computeSaddle.exe error - need 1,3,or 9 components for strain"
            end if
        end if
        call report(cla)
        if (.not. allRequiredArgumentsSet(cla)) stop
        if (hasHelpArgument(cla)) stop
        call delete(cla)
   
        if (getSuffix(outfile)=="xyz") outfile = removeSuffix(outfile)


    !---    welcome message
        print *,""
        print *,"ComputeSaddle.exe v"//trim(VERSION)
        print *,"^^^^^^^^^^^^^^^^^^^"//repeat("^",len_trim(VERSION))
        print *,""
        print *,"input filename1     : """//trim(filename1)//""""
        if (secondFile) then
            print *,"input filename2     : """//trim(filename2)//""""
            print *,"central atom        : ",centralAtom
            print *,"sphere radius       : ",R0
        end if
        if (deltax(1) /= LIB_CLA_NODEFAULT_R) then
            print *,"central atom        : ",centralAtom
            print *,"displace atom       : ",deltax
        end if
        print *,"output filename     : """//trim(outfile)//""""
        print *,"potential filename  : """//trim(potential)//""""
        print *,"start with CG relax : ",relax
        print *,"debug               : ",SIMPLEMD_DBG
        print *,"op electron density : ",opElectronDensity
        
        if (any(abs(strain)>0)) then
            print *,"strain          : ",strain(1,:)
            print *,"                : ",strain(2,:)
            print *,"                : ",strain(3,:)            
        end if    
         
        print *,""


 


    !---
        call start(tt(TIMER_IO))
        print *,""
        print *,"reading potential file"
        print *,"^^^^^^^^^^^^^^^^^^^^^^"
        print *,""
        call input(xml,potential,ok)

        if (ok) then
            call inputFromXML(xml,eam,ok)
            print *,"computeSaddle.exe info - potential read from file"
            !call report(eam)
        else
            print *,"computeSaddle.exe warning - potential not read from file"
            stop
        end if
        rc = getCutoff(eam)
        if (disp_scale == LIB_CLA_NODEFAULT_R) disp_scale = rc/4
        call delete(xml)
        print *,""
        call pause(tt(TIMER_IO))
    !---





    !---
        call start(tt(TIMER_IO))
        print *,""
        print *,"reading input file 1"
        print *,"^^^^^^^^^^^^^^^^^^^^"
        print *,""
        xyz = XYZFile_ctor(filename1)
        call readHeader(xyz,ok)
        call input(xyz,verbose=.true.)
        call getSupercell(xyz,a_super,ok)
        if (ok) then
            print *,"computeSaddle.exe info - supercell read from file"
            print *,a_super(1,:)
            print *,a_super(2,:)
            print *,a_super(3,:)
        else
            print *,"computeSaddle.exe warning - supercell not read from file"
        end if
        call pause(tt(TIMER_IO))
        call report(xyz)
        
        if (opxyz) then
            call setLammpsFormat(xyz,.false.)
            if (opElectronDensity) then
                call setNColumns(xyz,6)
            else
                call setNColumns(xyz,4)
            end if
        end if
        call getColumnsp(xyz,xx)
        call getTypesp(xyz,at)
        nAtoms = getNatoms(xyz)
        allocate(x_old(3,nAtoms))
       
        if (any(abs(strain)>0)) then
            a_super = a_super + matmul(strain,a_super)
            print *,"computeSaddle.exe info - strained supercell"
            print *,a_super(1,:)
            print *,a_super(2,:)
            print *,a_super(3,:)
            x_old(1:3,1:nAtoms) = xx(1:3,1:nAtoms)
            do ii = 1,nAtoms
                yy(1:3) = strain(1:3,1)*xx(1,ii) + strain(1:3,2)*xx(2,ii) + strain(1:3,3)*xx(3,ii)
                xx(1:3,ii) = xx(1:3,ii) + yy(1:3)
            end do
            
            if (.not. secondfile) then
                print *,"atoms "
                do ii = 1,min(nAtoms,10)
                    print *,ii,at(ii),x_old(1:3,ii),xx(1:3,ii) 
                end do
                print *,"..."
                print *,""
            end if            
            
        end if
        x_old(1:3,1:nAtoms) = xx(1:3,1:nAtoms)
                
    !   check that input element ordering is the same as the potential file ordering
        nAtomNames = getnAtomNames(xyz)
        do ii = 1,getNTypes( eam )
            atomNames(ii) = getName(eam,ii)
        end do
        at = -at        !   set all atom types negative in order to renumber them
        do jj = 1,nAtomNames
            ok = .false.
            do ii = 1,getNTypes( eam )
                if ( trim(getAtomName(xyz,jj))==trim(atomNames(ii)) ) then
                    print *,"computeSaddle.exe info - atom name """//trim(getAtomName(xyz,jj))//""" in .xyz file is type ",ii," in potential file"
                    where (at == -jj)
                        at = ii
                    end where
                    ok = .true.
                end if
            end do
            if (.not. ok) then
                print *,"computeSaddle.exe warning - atom name """//trim(getAtomName(xyz,jj))//""" in .xyz file not matched in potential file"
                if (any(at==-jj)) then
                    print *,"computeSaddle.exe ERROR - atoms of type """//trim(getAtomName(xyz,jj))//""" detected in .xyz file"
                    stop
                end if
            end if
        end do
        call setAtomNames( xyz,atomNames(1:getNTypes( eam )) )

        ii = getNcolumns(xyz)
        nAtoms = getNatoms(xyz)

    !---
        print *,""


    
    
        
        
        
        if (secondFile) then
            call start(tt(TIMER_IO))
            print *,""
            print *,"reading input file 2"
            print *,"^^^^^^^^^^^^^^^^^^^^"
            print *,""
            xyz = XYZFile_ctor(filename2)
            call readHeader(xyz,ok)
            call input(xyz,verbose=.true.)
            !call getSupercell(xyz,a_super,ok)
            !if (ok) then
            !    print *,"computeSaddle.exe info - supercell read from file"
            !    print *,a_super(1,:)
            !    print *,a_super(2,:)
            !    print *,a_super(3,:)
            !else
            !    print *,"computeSaddle.exe warning - supercell not read from file"
            !end if
            call pause(tt(TIMER_IO))        
            if (getNatoms(xyz)/=nAtoms) stop "computeSaddle.exe error - number of atoms doesn't agree in files" 
            print *,"computeSaddle.exe info - nColumns file 1 ",size(xx,dim=1)," nColumns file 2 ",getNcolumns(xyz)
            
            call getColumnsp(xyz,xx)
            !call getTypesp(xyz,at)
            allocate(x_new(3,nAtoms))
            
            if (any(abs(strain)>0)) then
                do ii = 1,nAtoms
                    yy(1:3) = strain(1:3,1)*xx(1,ii) + strain(1:3,2)*xx(2,ii) + strain(1:3,3)*xx(3,ii)
                    xx(1:3,ii) = xx(1:3,ii) + yy(1:3)
                end do
            end if
            
            x_new(1:3,1:nAtoms) = xx(1:3,1:nAtoms)
            
    
            print *,"atoms "
            do ii = 1,min(nAtoms,10)
                write (*,fmt='(a,i8,a,i4,a,3f12.4,a,3f12.4)') "atom ",ii," type ",at(ii)," from ",x_old(1:3,ii)," to ",x_new(1:3,ii) 
            end do
            print *,"..."
            if (centralAtom>10) then
                write (*,fmt='(a,i8,a,i4,a,3f12.4,a,3f12.4)') "atom ",centralAtom," type ",at(centralAtom)," from ",x_old(1:3,centralAtom)," to ",x_new(1:3,centralAtom) 

                print *,"..."
            end if
            print *,""
            
        
        
        
        end if
            
            
    !---    construct MD object
        print *,""
        print *,"construct simpleMD"
        print *,"^^^^^^^^^^^^^^^^^^"
        print *,""
        allocate(vv(3,nAtoms))
        vv = 0.0d0
        smd = MD_ctor(xx,at,vv,a_super,eam,dt=1.0d0)
        call report(smd,verbose=.true.) 
        print *,""
    !---
        
 
         
     
        if (relax) then
            print *,""
            print *,"relax file 1"
            print *,"^^^^^^^^^^^^"
            print *,""
            call setPositions(smd,x_old)
            call CGrelax(smd,e_old,disp_scale = disp_scale)
            x_old(1:3,1:nAtoms) = xx(1:3,1:nAtoms)
            print *,"computeSaddle.exe info - CG relax returns before ",e_old
            
            if (elasticConst) then
                CC = elasticConstants(smd)
                write(unit=*,fmt='(a)')   "computeSaddle.exe info -  elastic constants (eV/A^3)"
                write(unit=*,fmt='(6f16.5)') CC(1,1:6)
                write(unit=*,fmt='(6f16.5)') CC(2,1:6)
                write(unit=*,fmt='(6f16.5)') CC(3,1:6)
                write(unit=*,fmt='(6f16.5)') CC(4,1:6)
                write(unit=*,fmt='(6f16.5)') CC(5,1:6)
                write(unit=*,fmt='(6f16.5)') CC(6,1:6)

                c11 = (cc(1,1)+cc(2,2)+cc(3,3))/3
                c12 = (cc(1,2)+cc(1,3)+cc(2,3))/3
                c44 = (cc(4,4)+cc(5,5)+cc(6,6))/3
    
                print *,"c11          ",c11," = ",c11*160.2176530d0,"GPa"
                print *,"c12          ",c12," = ",c12*160.2176530d0,"GPa"
                print *,"c44          ",c44," = ",c44*160.2176530d0,"GPa"
                print *,"mu           ",(c11-c12)/2," = ",(c11-c12)/2*160.2176530d0,"GPa"
                print *,"nu           ",c12/(c11+c12)
                print *,""
            end if
            
            if (dipoleTens) then
                PP = dipoleTensor(smd)
                write(unit=*,fmt='(a)')   "computeSaddle.exe info -  dipole tensor (eV)"
                write(unit=*,fmt='(6f16.5)') PP(1,1:3)
                write(unit=*,fmt='(6f16.5)') PP(2,1:3)
                write(unit=*,fmt='(6f16.5)') PP(3,1:3)
                
                if (elasticConst) then
                
                    ss(1) = PP(1,1)
                    ss(2) = PP(2,2)
                    ss(3) = PP(3,3)
                    ss(4) = ( PP(2,3) + PP(3,2) )/ 2
                    ss(5) = ( PP(3,1) + PP(1,3) )/ 2
                    ss(6) = ( PP(1,2) + PP(2,1) )/ 2 
                   
                    vol =  determinant3Mat(a_super)
                                       
                    ss = ss / vol
                    
                    call DSYSV( "L",6,1,CC,6,ipiv,ss,6,work,size(work),ii )
                    
                !---    stretch tensor = 1 + strain  
                    FF(1:3,1) = (/ 1+ss(1),ss(6),ss(5) /)
                    FF(1:3,2) = (/ ss(6),1+ss(2),ss(4) /)
                    FF(1:3,3) = (/ ss(5),ss(4),1+ss(3) /)
                    
                    
                    ee = determinant3Mat(FF)
                    
                    write(unit=*,fmt='(a,f16.5,a)') "volume supercell ",vol," A^3"
                    write(unit=*,fmt='(a,f16.5,a)') "volume change    ",(ee-1)*vol," A^3"
                    
                    
                !---    strain tensor
                    FF(1,1) = FF(1,1) - 1.0d0
                    FF(2,2) = FF(2,2) - 1.0d0
                    FF(3,3) = FF(3,3) - 1.0d0
                    
                    write(unit=*,fmt='(a)')   "computeSaddle.exe info -  strain tensor (eV)"
                    write(unit=*,fmt='(6f16.5)') FF(1,1:3)
                    write(unit=*,fmt='(6f16.5)') FF(2,1:3)
                    write(unit=*,fmt='(6f16.5)') FF(3,1:3)
                    
                    
                    FF = matmul( FF,PP )
                    ee = -( FF(1,1) + FF(2,2) + FF(3,3) )/(2*vol)
                    write(unit=*,fmt='(a,f16.8,a)') "energy change    ",ee," eV"
                    write(unit=*,fmt='(a,f16.5,a)') "energy inc elas  ",e_old + ee," eV"
                                                              
                    
                end if                            
                    
            end if
            
            
            
        end if
        
        if (opxyz) then
            
            call setLammpsFormat(xyz,.false.)
            if (opElectronDensity) then
                call setNColumns(xyz,6)
            else
                call setNColumns(xyz,4)
            end if
!            call setNColumns(xyz,4)
            call getColumnsp(xyz,xx)     
            xx(1:3,1:nAtoms) = x_old(1:3,1:nAtoms)    
            call setPositions(smd,xx)
            write(dummy,fmt='(9f14.5)') a_super
            if (opElectronDensity) then
                dummy="Lattice="""//trim(dummy)//""" Properties=Species:S:1:Pos:R:3:Energy:R:1:rho:R:1:F:R:1"
            else
                dummy="Lattice="""//trim(dummy)//""" Properties=Species:S:1:Pos:R:3:Energy:R:1"
            end if
            call setColumn_description(xyz,dummy)
            call setnHeaderLines(xyz,0)
            call potentialEnergyArray(smd,xx(4,:))
            if (opElectronDensity) then
                call electronDensityArray(smd,xx(5,:))
                call embeddingEnergyArray(smd,xx(6,:))
            end if
            call start(tt(TIMER_IO))
            call setFilename(xyz,trim(outfile)//".relax1.xyz" )
            call output(xyz)
            call pause(tt(TIMER_IO))
        end if
        if ((.not. secondFile) .and. (deltax(1)==LIB_CLA_NODEFAULT_R)) then
            print *,""
            print *,"done"
            print *,""
            stop
        end if
    
        if (.not. secondFile) then
            allocate(x_new(3,nAtoms))
            x_new(1:3,1:nAtoms) = x_old(1:3,1:nAtoms)
            x_new(1:3,centralAtom) = x_new(1:3,centralAtom) + deltax(1:3)
            !print *,"displace atom ",centralAtom," by ",deltax
            write (*,fmt='(a,i8,a,i4,a,3f12.4,a,3f12.4)') "atom ",centralAtom," type ",at(centralAtom)," from ",x_old(1:3,centralAtom)," to ",x_new(1:3,centralAtom) 
        end if
        
        if (relax .or. .not. secondFile) then
        
            print *,""
            print *,"relax file 2"
            print *,"^^^^^^^^^^^^"
            print *,""
            
            
            call setPositions(smd,x_new(1:3,1:nAtoms))
            call CGrelax(smd,e_new,disp_scale = disp_scale)
            x_new(1:3,1:nAtoms) = xx(1:3,1:nAtoms)
            
            print *,"computeSaddle.exe info - CG relax returns after ",e_new," (eV)"
            write (*,fmt='(a,i8,a,i4,a,3f12.4,a,3f12.4)') "atom ",centralAtom," type ",at(centralAtom)," from ",x_old(1:3,centralAtom)," to ",x_new(1:3,centralAtom) 
        end if
        
        
        if (opxyz) then
            call setNColumns(xyz,4)
            call getColumnsp(xyz,xx)     
            xx(1:3,1:nAtoms) = x_new(1:3,1:nAtoms)    
            call setPositions(smd,xx)
            write(dummy,fmt='(9f14.5)') a_super
            !dummy="Lattice="""//trim(dummy)//""" Properties=Species:S:1:Pos:R:3:Energy:R:1"
            if (opElectronDensity) then
                dummy="Lattice="""//trim(dummy)//""" Properties=Species:S:1:Pos:R:3:Energy:R:1:rho:R:1:F:R:1"
            else
                dummy="Lattice="""//trim(dummy)//""" Properties=Species:S:1:Pos:R:3:Energy:R:1"
            end if
            call setColumn_description(xyz,dummy)
            call setnHeaderLines(xyz,0)
            call potentialEnergyArray(smd,xx(4,:))
            if (opElectronDensity) then
                call electronDensityArray(smd,xx(5,:))
                call embeddingEnergyArray(smd,xx(6,:))
            end if
            call start(tt(TIMER_IO))
            call setFilename(xyz,trim(outfile)//".relax2.xyz" )
            call output(xyz)
            call pause(tt(TIMER_IO))
        end if
                    
        
        ee = 1000.0d0       !   "insane" barrier energy level
        allocate(ee_out(nReplicas))
        
            print *,""
            print *,"compute saddle"
            print *,"^^^^^^^^^^^^^^"
            print *,""
        if (opxyz) then
            call computeSaddlePointEnergy( smd,x_old,x_new,centralAtom,abs(R0),e_old,nReplicas, ee,leaveInSaddle = .true.,ee_out = ee_out,outfile = outfile,disp_scale = disp_scale )
        else
            call computeSaddlePointEnergy( smd,x_old,x_new,centralAtom,abs(R0),e_old,nReplicas, ee,leaveInSaddle = .true.,ee_out = ee_out,disp_scale = disp_scale )
        end if        
        x_new(1:3,1:nAtoms) = xx(1:3,1:nAtoms)
        print *,""
        write(*,fmt='(a8,2a16)') "replica","reaction coord","energy"
        do ii = 1,nReplicas
            write(*,fmt='(i8,2f16.5)') ii,(ii-1)*1.0d0/(nReplicas-1),ee_out(ii)
        end do
        print *,""
        
        
        if (opxyz) then
!            call setNColumns(xyz,4)
!            call getColumnsp(xyz,xx)     
            !xx(1:3,1:nAtoms) = x_new(1:3,1:nAtoms)    
            call setPositions(smd,x_new(1:3,1:nAtoms))
            write(dummy,fmt='(9f14.5)') a_super
            if (opElectronDensity) then
                dummy="Lattice="""//trim(dummy)//""" Properties=Species:S:1:Pos:R:3:Energy:R:1:rho:R:1:F:R:1"
            else
                dummy="Lattice="""//trim(dummy)//""" Properties=Species:S:1:Pos:R:3:Energy:R:1"
            end if
            !dummy="Lattice="""//trim(dummy)//""" Properties=Species:S:1:Pos:R:3:Energy:R:1"
            call setColumn_description(xyz,dummy)
            call setnHeaderLines(xyz,0)
            call potentialEnergyArray(smd,xx(4,:))
            if (opElectronDensity) then
                call electronDensityArray(smd,xx(5,:))
                call embeddingEnergyArray(smd,xx(6,:))
            end if
            call start(tt(TIMER_IO))
            call setFilename(xyz,trim(outfile)//".neb.xyz" )
            call output(xyz)
            call pause(tt(TIMER_IO))
        end if
 

        call delete(smd)
        call delete(xyz)

    !---
    
    
        call pause(tt(TIMER_ALL)) 
        print *,""        
        print *,"timing data:"
        write(*,fmt='(a,f12.3)') "   total     ",elapsed(tt(TIMER_ALL))
        if (relax) then
            write(*,fmt='(a,f12.3)') "   CG        ",elapsed(tt(TIMER_CG))
        end if
        write(*,fmt='(a,f12.3)') "   saddle    ",elapsed(tt(TIMER_SADDLE))
        write(*,fmt='(a,f12.3)') "   IO        ",elapsed(tt(TIMER_IO))
        print *,""        
        
    
    
    
    
        print *,""
        print *,"done"
        print *,""
        
        
    contains
!---^^^^^^^^


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
                    
 
    end program computeSaddle