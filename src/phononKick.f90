
    program phononKick
!---^^^^^^^^^^^^^^^^^^
        use Lib_CommandLineArguments
        use Lib_XYZFiles
        use Lib_SimpleSupercells
        use Lib_LinkCell3d
        use Lib_NBAX
        use Lib_Callipers
        use EAM_PotentialTables
        use iso_fortran_env
        implicit none 
        
        external                        ::      SSYSV                   !   lapack
        
        
        type(CommandLineArguments)      ::      cla

    !---    input parameters        
        character(len=256)              ::      filename = ""           !   input filename
        character(len=256)              ::      outfile = ""            !   output filename        
        character(len=256)              ::      potential = ""          !   potential filename
        integer                         ::      atom = LIB_CLA_NODEFAULT_I              !   number of atom to kick
        real(kind=real64),dimension(3)  ::      v = LIB_CLA_NODEFAULT_R                 !   number of random thermal configs to take
        real(kind=real64)               ::      R0 = 10.0d0             !   radius to compute moving atoms
        integer                         ::      vcol = 5                !   column where atom velocities are found/stored (default 5 = atom x y z vx vy vz)
        real(kind=real64)               ::      deltax = 1.0d0          !   distance to shift atom by velocity v before recalc
        
    !---    physically meaningful parameters deduced from input file       
        real(kind=real64),dimension(3,3)    ::      a_super,a_cell      !   size of periodic supercell, size of link cell
        integer                             ::      nAtoms              !   number of atoms
        integer                             ::      nColumns            !   number of .xyz file columns
        integer                             ::      Nx,Ny,Nz            !   number of link cells in supercell
        type(XYZFile)                       ::      xyz
        type(LinkCell3d)                    ::      lc3d
        type(NBAX)                          ::      xml
        type(EAM_PotentialTable)            ::      eam
        real(kind=real64),dimension(:,:),pointer    ::      colp,colp0
        integer,dimension(:),pointer        ::      tp
        real(kind=real64)                   ::      rc                  !   potential cutoff range used for link cell list  
        real(kind=real64)                   ::      tt                  !   time at which to compute velocity nudge.
        real(kind=real64)                   ::      ke_before,ke_after
        real(kind=real64)                   ::      mass
                
    !---    output values computed from input file
        integer                             ::      nSphere,nCore             !   number of atoms in cutoff sphere 
        real(Kind=real64),dimension(:,:),allocatable    ::      u       !   velocity of atoms in cutoff sphere
        
        type(Callipers)                     ::      calc_hess,calc_ssysv,calc
        
        
    !---    dummy parameter     
        logical                             ::      ok
        integer                             ::      ii,jj,kk,nn,k2,step,coreatom
        
        integer,dimension(:),allocatable    ::      id,id2
        real(kind=real64),dimension(:,:),allocatable        ::      dx,dx2
        real(kind=real64),dimension(:),allocatable          ::      dr,dr2
        real(kind=real32),dimension(:,:),allocatable        ::      hess 
        real(kind=real64),dimension(:,:),allocatable        ::      force,f1
        real(kind=real64)                                   ::      e0,e1,e2
        logical,dimension(:),allocatable                    ::      core
        integer,dimension(:),allocatable                    ::      indx
        real(kind=real32),dimension(:),allocatable          ::      bb,work
        integer,dimension(:),allocatable                    ::      ipiv
        real(kind=real64),dimension(3)                      ::      ubar
        real(kind=real64)                                   ::      d1,d2,d3
        
    !---    read command line arguments
        cla = CommandLineArguments_ctor(10)  
         
        call setProgramDescription( cla, "phononKick" )
        call setProgramVersion( cla, "0.1" )               
        
        call setCategories(cla,(/ "files handling ",  &
                                  "calculation    ",  &
                                  "debug          " /))
        
        call get( cla,"f",filename ,LIB_CLA_REQUIRED,"       input filename",1 )        
        outfile = trim(filename)//".u.xyz"
        call get( cla,"o",outfile  ,LIB_CLA_OPTIONAL,"       output filename",1 )     
        call get( cla,"p",potential,LIB_CLA_REQUIRED,"       potential filename",1 )     
        call get( cla,"vcol",vcol  ,LIB_CLA_OPTIONAL,"       column where atom velocities are found/stored (default 5 = atom x y z vx vy vz)",1 )     

        call get( cla,"atom",atom  ,LIB_CLA_REQUIRED,"       atom to kick",2 )     
        ii = 3 ; call get( cla,"v",v  ,ii,LIB_CLA_REQUIRED,"       velocity of atom to kick",2 )     
        
        call get( cla,"R0",R0  ,LIB_CLA_OPTIONAL,"       sphere radius for velocity tweaks",2 )     
        call get( cla,"deltax",deltax  ,LIB_CLA_OPTIONAL,"       displacement for velocity calc",2 )     
        call report(cla)
        if (.not. allRequiredArgumentsSet(cla)) stop
        if (hasHelpArgument(cla)) stop
        call delete(cla)
    !---
    
      
    !---    screen dump input parameters
        print *,""
        print *,"file handling"
        print *,"    input filename         "//trim(filename)
        print *,"    output filename        "//trim(outfile)
        print *,"    potential filename     "//trim(potential)
        print *,"    velocity column        ",vcol
        print *,""
        print *,"calculation"
        print *,"    atom to kick           ",atom
        print *,"    velocity (A/fs)        ",v
        print *,"    sphere radius (A)      ",R0
        print *,"    deltax (A)             ",deltax
        print *,""
        
         
         
    !---     
        print *,""
        print *,"reading input file"
        print *,"^^^^^^^^^^^^^^^^^^"
        print *,""
        xyz = XYZFile_ctor(filename)
        call readHeader(xyz,ok)
        call input(xyz,verbose=.true.)
         
        call getSupercell(xyz,a_super,ok) 
        if (ok) then
            print *,"phononKick info - supercell read from file"
            print *,a_super(1,:)
            print *,a_super(2,:)
            print *,a_super(3,:)
        else
            print *,"phononKick warning - supercell not read from file"
        end if                
        print *,""
    !---      
    
    !---     
        print *,""
        print *,"reading potential file"
        print *,"^^^^^^^^^^^^^^^^^^^^^^"
        print *,""
        call input(xml,potential,ok)
         
        if (ok) then
            call inputFromXML(xml,eam,ok)
            print *,"phononKick info - potential read from file"
            call report(eam)
        else
            print *,"phononKick warning - potential not read from file"
            stop
        end if   
        call delete(xml)     
        print *,""
    !---       
    
    
    !---    
        print *,"checking input file"
        print *,"^^^^^^^^^^^^^^^^^^"
        print *,""
        nColumns = getNColumns(xyz)
        if (nColumns < vcol+1) then
            !   not enough columns to read velocity
            print *,"phononKick error - read ",nColumns," data columns after atom type. Expected velocity in cols ",vcol-1,":",vcol+1
            stop
        end if
        nAtoms = getNAtoms(xyz)     
        call getColumnsp(xyz,colp)    
        call getTypesp(xyz,tp)
        call report(xyz)   
        allocate(colp0(nColumns,nAtoms))
        colp0 = colp
        print *,""
        print *,"phononKick info - atom to kick"
        print *,"    number   ",atom,tp(atom)
        print *,"    type     ",tp(atom)," """//trim(getAtomName(xyz,tp(atom)))//""""
        print *,"    position ",colp(1:3,atom)," A"
        print *,"    v_before ",colp(vcol-1:vcol+1,atom)," A/fs"
        print *,"    v_after  ",v," A/fs"
        mass = getMass(eam)
        print *,"    mass     ",mass," eV (fs/A)^2"
        ke_before = dot_product(colp(vcol-1:vcol+1,atom),colp(vcol-1:vcol+1,atom))*mass/2
        print *,"    ke_before",ke_before," eV"
        
        colp(vcol-1:vcol+1,atom) = v(1:3)
        ke_after = dot_product(v,v)*mass/2
        print *,"    ke_after ",ke_after," eV"
        print *,""
    !---    
         
    
    
    !---   
        print *,""
        print *,"constructing link cell list"
        print *,"^^^^^^^^^^^^^^^^^^^^^^^^^^^"
        print *,"" 
        rc = getCutoff(eam)
        Nx = max(3,floor(norm2(a_super(1:3,1)/rc)))
        Ny = max(3,floor(norm2(a_super(1:3,2)/rc)))
        Nz = max(3,floor(norm2(a_super(1:3,3)/rc)))
        
        nn = estimateMaxCountPerCell(Nx*Ny*Nz,nAtoms) 
        a_cell(1:3,1) = a_super(1:3,1)/Nx
        a_cell(1:3,2) = a_super(1:3,2)/Ny
        a_cell(1:3,3) = a_super(1:3,3)/Nz         
        lc3d = LinkCell3D_ctor(a_cell,nn,Nx,Ny,Nz)
        do ii = 1,nAtoms
            call add( lc3d,ii,colp(1:3,ii))            
        end do
        call report(lc3d)
        print *,""
    !---
    
    
    !---  
        print *,""
        print *,"computing force, energy at initial position"
        print *,"^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^"
        print *,""    
        
        nn = getnNeighMax(lc3d)
        allocate(dx(3,nn))
        allocate(dr(nn))
        allocate(id(nn))
        
        
!       e0 = 0.0d0  
!       do ii = 1,nAtoms
!           
!           call neighbourList( lc3d,colp(1:3,ii),rc, nn,id,dx,dr )
!       
!       !   find self as a neighbour, and remove from neighbour list
!           jj = minloc(dr(1:nn),dim=1)       
!           do kk = jj+1,nn                    
!               id(kk-1) = id(kk) ; dx(:,kk-1) = dx(:,kk) ; dr(kk-1) = dr(kk)
!           end do
!           nn = nn - 1
!           
!       !   check which core atom is referenced, which neighbour     
!           do kk = 1,nn            
!           !   convert square distance and vector to distance and normal vector
!               dr(kk) = sqrt(dr(kk))
!               dx(:,kk) = dx(:,kk) / dr(kk)
!           end do
!          
!           e0 = e0 + getEnergy(eam, nn,dr(1:nn) )        
!       end do
!       print *,"total energy ",e0
!       
!       !LinkCell3D_dbg = .true.
!       do step = 0,20
!           tt = step*0.05*2.74/norm2(v)
!        
!           !print *,""
!           !print *,"step ",step
!           call cut(lc3d,atom,colp(1:3,atom),ok)        
!           colp(1:3,atom) = colp(1:3,atom)+v(1:3)*tt  
!           call add(lc3d,atom,colp(1:3,atom))
!           !print *,"atom moved: ",tt,colp(1:3,atom),ok
!           e1 = 0
!           
!           do ii = 1,nAtoms
!               call neighbourList( lc3d,colp(1:3,ii), rc, nn,id,dr )
!               
!           !   find self as a neighbour, and remove from neighbour list
!               jj = minloc(dr(1:nn),dim=1)       
!               do kk = jj+1,nn                    
!                   dr(kk-1) = dr(kk)
!               end do
!               nn = nn - 1
!                   
!           !   convert square distance and vector to distance and normal vector
!               do kk = 1,nn                            
!                   dr(kk) = sqrt(dr(kk))
!               end do
!               !if (ii==atom) print *,"central atom ",ii,nn,dr(1:nn)
!                               
!               e1 = e1 + getEnergy(eam, nn,dr(1:nn) )        
!               !if (ii==atom) print *,"central atom ",ii,nn,getEnergy(eam, nn,dr(1:nn) ) 
!           end do
!                
!           !print *,""  
!           call cut(lc3d,atom,colp(1:3,atom))         
!           colp(1:3,atom) = colp0(1:3,atom)
!           call add(lc3d,atom,colp(1:3,atom))
!           print *,"potential energy drag  ",tt,norm2(v)*tt,e1-e0
!           !print *,""
!       end do
!       print *,""
!           
!               
        
        
       ! stop
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        !print *,"long neighbour list calc..."
        calc = Callipers_ctor()
        
        call neighbourList_long( lc3d,colp(1:3,atom),R0, nSphere )
        allocate(dr2(nSphere))
        allocate(id2(nSphere))
        call neighbourList_long( lc3d,colp(1:3,atom),R0, nSphere,id2,dr2 )
        !print *,"...done"
        print *,"number of atoms in sphere ",nSphere      
        
        allocate(indx(nAtoms)) ; indx(:) = -1
        
    !---    identify the atoms in the core region
        nCore = 0        
        do k2 = 1,nSphere     
            if ( dr2(k2) < (R0-rc)*(R0-rc) ) then
                ii = id2(k2)
                nCore = nCore + 1
                indx(ii) = nCore
                if (ii==atom) coreatom = nCore
            end if
        end do
        where (indx==-1)
            indx = nCore+1
        end where
        print *,"core atom count ",nCore," central atom ",coreatom
        
        
    !---    compute energy in sphere region
        allocate(force(3,nCore+1))
        force = 0.0d0  
        e0 = 0.0d0  
        do k2 = 1,nSphere
            
            call neighbourList( lc3d,colp(1:3,id2(k2)),rc, nn,id,dx,dr )
        
        !   find self as a neighbour, and remove from neighbour list
            jj = minloc(dr(1:nn),dim=1)       
            do kk = jj+1,nn                    
                id(kk-1) = id(kk) ; dx(:,kk-1) = dx(:,kk) ; dr(kk-1) = dr(kk)
            end do
            nn = nn - 1
            
        !   check which core atom is referenced, which neighbour     
            do kk = 1,nn            
                id(kk) = indx( id(kk) )             !   either pointing to a core atom or to nCore+1            
                
            !   convert square distance and vector to distance and normal vector
                dr(kk) = sqrt(dr(kk))
                dx(:,kk) = dx(:,kk) / dr(kk)
            end do
            
            !write(*,fmt='(a,2i6,a,1000i6)') "sphere ",k2,ii," neigh ",id(1:nn)
                         
            call addForce(eam, indx(k2) ,nn,dr(1:nn),dx(1:3,1:nn),id(1:nn), force )
            e0 = e0 + getEnergy(eam, nn,dr(1:nn) )        
        end do
        call pause(calc)
        !call pause(calc_hess)
        
      !  do k2 = 1,nSphere
      !      print *,k2,id2(k2),core(k2),sqrt(dr2(k2))
      !  end do
        
        
         
       
       
    !---    cut all from sphere
       ! do k2 = 1,nSphere
       !     ii = id2(k2)
       !     call cut( lc3d,ii,colp(1:3,ii) , ok) 
       ! end do
        
        do step = 0,20
            tt = step*0.05*2.74/norm2(v)
        
        !---    add all atoms to sphere, in slightly displaced positions    
           ! do k2 = 1,nSphere
           !     ii = id2(k2)                                
           !     call add( lc3d,ii,colp(1:3,ii) + colp(vcol-1:vcol+1,ii)*tt)
           ! end do
            
           call cut(lc3d,atom,colp(1:3,atom)) 
           colp(1:3,atom) = colp(1:3,atom) +  v(1:3)*tt        
           call add(lc3d,atom,colp(1:3,atom))
           
            e1 = 0
            
            do k2 = 1,nSphere
                call neighbourList( lc3d,colp(1:3,id2(k2)),rc, nn,id,dx,dr )
                
            !   find self as a neighbour, and remove from neighbour list
                jj = minloc(dr(1:nn),dim=1)       
                do kk = jj+1,nn                    
                    id(kk-1) = id(kk) ; dx(:,kk-1) = dx(:,kk) ; dr(kk-1) = dr(kk)
                end do
                nn = nn - 1
                    
            !   check which core atom is referenced, which neighbour     
                do kk = 1,nn            
                !   convert square distance and vector to distance and normal vector
                    dr(kk) = sqrt(dr(kk))
                end do
            
                                
                e1 = e1 + getEnergy(eam, nn,dr(1:nn) )        
            end do
                   
        !---    cut all atoms from sphere again
          !  do k2 = 1,nSphere
          !      ii = id2(k2)                                
          !      call cut( lc3d,ii,colp(1:3,ii) + colp(vcol-1:vcol+1,ii)*tt)
          !  end do
            call cut(lc3d,atom,colp(1:3,atom))         
            colp(1:3,atom) = colp0(1:3,atom) 
            call add(lc3d,atom,colp(1:3,atom))
            print *,"potential energy drag  ",tt,norm2(v)*tt,e1-e0
        end do
        print *,""
            
                
                
    !---    add all to sphere
       ! do k2 = 1,nSphere
       !     ii = id2(k2)
       !     call add( lc3d,ii,colp(1:3,ii)) 
       ! end do
        
        
       
       
       
        call start(calc)
        calc_hess = Callipers_ctor()
        tt = deltax/norm2(v)
        print *,"computing atom positions at time ",tt," fs"
        print *,"central atom delta = ",tt*v(1:3)," A"
        call cut(lc3d,atom,colp(1:3,atom))        
        colp(1:3,atom) = colp(1:3,atom)+v(1:3)*tt
        call add(lc3d,atom,colp(1:3,atom))

        allocate(hess(3*(nCore+1),3*(nCore+1)))
        allocate(f1(3,nCore+1))
        !allocate(u(3,nSphere+1))
        !u = 0.0d0        
        f1 = 0
        e1 = 0
        hess = 0
        do k2 = 1,nSphere
            
            call neighbourList( lc3d,colp(1:3,id2(k2)),rc, nn,id,dx,dr )
            
        !   find self as a neighbour, and remove from neighbour list
            jj = minloc(dr(1:nn),dim=1)       
            do kk = jj+1,nn                    
                id(kk-1) = id(kk) ; dx(:,kk-1) = dx(:,kk) ; dr(kk-1) = dr(kk)
            end do
            nn = nn - 1
            
        !   check which core atom is referenced, which neighbour     
            do kk = 1,nn            
                id(kk) = indx( id(kk) )             !   either pointing to a core atom or to nCore+1            
                
            !   convert square distance and vector to distance and normal vector
                dr(kk) = sqrt(dr(kk))
                dx(:,kk) = dx(:,kk) / dr(kk)
            end do
            
            call addForce(eam, indx(k2),nn,dr(1:nn),dx(1:3,1:nn),id(1:nn), f1 )
            call addHessian(eam, indx(k2),nn,dr(1:nn),dx(1:3,1:nn),id(1:nn), hess )
            e1 = e1 + getEnergy(eam, nn,dr(1:nn) )        
        end do
        call pause(calc)       
        call pause(calc_hess)
        print *,""
        print *,"potential energy before ",e0
        print *,"potential energy after  ",e1
        print *,"number of atoms in sphere ",nSphere
        print *,"number of atoms in core   ",nCore
        
         
        
    !---    construct Lagrange multipler problem hess uu = bb
        call start(calc)
        calc_ssysv = Callipers_ctor()
        allocate(bb(3*nCore+2))
        allocate(work(66*(3*nCore+2)))
        allocate(ipiv(3*nCore+2))
        bb(:) = 0.0d0
        hess(:,3*nCore+1:3*nCore+3) = 0.0d0
        do k2 = 1,nCore
            bb(k2*3-2:k2*3) =  real( f1(:,k2)/tt,kind=real32 ) 
            hess(3*k2-2:3*k2,3*nCore+2) = real( force(1:3,k2),kind=real32 )
        end do    
        hess(3*coreAtom-2:3*coreAtom,3*nCore+1) = real( v(1:3),kind=real32 )                        
        
        call SSYSV( "U",3*nCore+2,1,hess,3*nCore+3,ipiv,bb,3*nCore+2,work,size(work),ii )      !   note: LDA = 3*nCore+3
        print *,"SSYSV returns info = ",ii
        call pause(calc)
        call pause(calc_ssysv)
        
    !---    check constraints
        bb(3*nCore+1:) = 0.0d0
        ke_after = sum( bb*bb )*mass/2
        print *,"ke in u ",ke_after," eV"
        ubar = 0.0d0 ; d1 = 0.0d0 ; d2 = 0.0d0 ; d3 = 0.0d0
        do k2 = 1,nSphere
            ii = id2(k2)        
            jj = indx(ii)
            if (jj <= nCore) d1 = d1 + dot_product( bb(jj*3-2:jj*3),colp(vcol-1:vcol+1,ii) )
        end do
        do k2 = 1,nCore             
            ubar = ubar + bb(k2*3-2:k2*3)
            d2 = d2 + dot_product( bb(k2*3-2:k2*3),f1(:,k2) )
            d3 = d3 + dot_product( bb(k2*3-2:k2*3),force(:,k2) )
        end do
        ubar = ubar / nCore
        print *,"<u>     ",ubar
        print *,"u.v     ",d1
        print *,"u.f     ",d2
        print *,"u.f0    ",d3
        
        
        call start(calc)
        bb(coreatom*3-2:coreatom*3) = bb(coreatom*3-2:coreatom*3) + v(1:3)
        do k2 = 1,nSphere
            ii = id2(k2)        
            jj = indx(ii)
            if (jj <= nCore) colp(vcol-1:vcol+1,ii) = colp(vcol-1:vcol+1,ii) + bb(jj*3-2:jj*3)
        end do   
        call pause(calc) 
        !colp(vcol-1:vcol+1,atom) = v(1:3)
        
        
        
        
        
        
    !---    full output of forces and velocities
    !   write(*,fmt='(3a6,a16,100a48)') "#n ","id","core","r","f_before","f_after","u","velocity"
    !   do k2 = 1,nSphere
    !       ii = id2(k2)
    !        if (core(k2)) &
    !        write(*,fmt='(2i6,l6,100f16.8)') k2,ii,core(k2),sqrt(dr2(k2)),force(:,k2),f1(:,k2),bb(k2*3-2:k2*3),colp(vcol-1:vcol+1,ii)
    !   end do
    !   print *,""
        
                    
    !---    now move all the atoms in the core using their new velocities.
    !       note: not MD, I'm just making a position shift
        call cut(lc3d,atom,colp(1:3,atom))
        colp(1:3,atom) = colp0(1:3,atom)
        call add(lc3d,atom,colp(1:3,atom))
        print *,"potential energy before ",e0
        print *,"potential energy after  ",e1
        
        
     !---    cut all from sphere
         do k2 = 1,nSphere
            ii = id2(k2)        
            jj = indx(ii)
            if (jj <= nCore) call cut( lc3d,ii,colp(1:3,ii) , ok) 
         end do
        
        do step = 0,20
            tt = step*0.05*2.74/norm2(v)
         !  call cut(lc3d,atom,colp(1:3,atom))       
         !  call add(lc3d,atom,colp(1:3,atom)+v(1:3)*tt)  
         
         
     !---    add all from sphere
         do k2 = 1,nSphere
            ii = id2(k2)        
            jj = indx(ii)
            
            if (jj <= nCore) then
                colp(1:3,ii) = colp0(1:3,ii) + bb(jj*3-2:jj*3)*tt
                call add( lc3d,ii,colp(1:3,ii))
            end if
         end do
                       
          ! do k2 = 1,nSphere
          !     ii = id2(k2)                                
          !     call add( lc3d,ii,colp(1:3,ii) + colp(vcol-1:vcol+1,ii)*tt)
          ! end do
          !
            e2 = 0
            
            do k2 = 1,nSphere
                ii = id2(k2)
                call neighbourList( lc3d,colp(1:3,ii),rc, nn,id,dx,dr )
                
            !   find self as a neighbour, and remove from neighbour list
                jj = minloc(dr(1:nn),dim=1)       
                do kk = jj+1,nn                    
                    id(kk-1) = id(kk) ; dx(:,kk-1) = dx(:,kk) ; dr(kk-1) = dr(kk)
                end do
                nn = nn - 1
                    
            !   check which core atom is referenced, which neighbour     
                do kk = 1,nn            
                !   convert square distance and vector to distance and normal vector
                    dr(kk) = sqrt(dr(kk))
                end do
            
                
                
                e2 = e2 + getEnergy(eam, nn,dr(1:nn) )        
            end do
                
     !---    cut all from sphere
         do k2 = 1,nSphere
            ii = id2(k2)        
            jj = indx(ii)
            if (jj <= nCore) call cut( lc3d,ii,colp(1:3,ii) , ok) 
         end do
        
            print *,"potential energy shift  ",tt,norm2(v)*tt,e2-e0
        end do
        print *,""
            
        
        
        
                
    !---    add all to sphere
         do k2 = 1,nSphere
             ii = id2(k2)
             colp(1:3,ii) = colp0(1:3,ii)
             call add( lc3d,ii,colp(1:3,ii)) 
         end do
        
        
        
    !---    output
        
        call setFilename(xyz,outfile)
        call output(xyz)
        call report(xyz)
    
    
        print *,"core calculation time ( forces + hessian + SSYSV ) "
        call report(calc)
        print *,"Hessian calculation time   "
        call report(calc_hess)
        print *,"SSYSV calculation time     " 
        call report(calc_ssysv)
        print *,""
        write(unit=0,fmt='(a)') "done"
        print *,""
        
        
        

    end program phononKick
         
         
         
         
         
         
         