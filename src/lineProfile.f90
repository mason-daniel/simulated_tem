
!
!    note: debye waller factors for elemental crystals can be found in Peng et al Acta Cryst. (1996). A52, 456-470
!       B = 8 pi^2 <u_x^2>
!    values at T=280K:
!           Li      Be       C      Na      Mg      A1      Si      K       Ca      Ca      Sc      Ti      V       Cr                                                                                   
!           4.5958  0.4042  0.1422  6.3648  1.7319  0.7460  0.4634  10.2398 1.9086  2.4847  0.7125  0.4911  0.5465  0.2398
!
!           Fe(bcc) Fe(fcc) Ni      Cu      Zn      Ge      Rb      Sr      Y       Zr      Nb      Mo      Pd
!           0.3106  0.5330  0.3401  0.5261  1.0915  0.5734  12.6041 3.6221  0.8176  0.5435  0.4297  0.2059  0.4259
!
!           Ag      Sn      Cs      Ba      La      Tb      Ho      Ta      W       Pt      Au      Pb      Th
!           0.7017  1.0825  16.7698 2.9279  1.7739  0.9660  0.8023  0.3078  0.1526  0.3557  0.5925  2.0341  0.6977
!
!   so reasonable thermal factors to include at room temp would be
!   
!       <u_x^2>_W  = 0.0019
!       <u_x^2>_Zr = 0.0069
!       <u_x^2>_Fe = 0.0039
!

    program LineProfile
!---^^^^^^^^^^^^^^^^^^^
!*      compute an x-ray power diffraction line profile from an input xyz file
!*  
        use Lib_CommandLineArguments
        use Lib_XYZFiles
        use Lib_Callipers
!        use Lib_TriangulatedSurface
        use Lib_SimpleProgressBar
        use Lib_RandomSeed
        use Lib_RotationMatrices
        use Lib_Voxelise
        use Lib_SimpleSupercells
        use Lib_FFTW3f
!        use OMP_LIB
        use iso_fortran_env
        implicit none
         
        type(CommandLineArguments)      ::      cla
        real(kind=real64),parameter     ::      PI = 3.141592653590d0
        complex(kind=real64),parameter  ::      I = cmplx( 0.0d0,1.0d0,kind=real64 )
        complex(kind=real64),parameter  ::      ZERO = cmplx( 0.0d0,0.0d0,kind=real64 )
        complex(kind=real64),parameter  ::      ONE = cmplx( 1.0d0,0.0d0,kind=real64 )
        character(len=8)                ::      VERSION = "1.0.0"
    !---    input parameters        
        character(len=256)              ::      filename = ""           !   input filename
        real(kind=real64)               ::      sigma = 0.0d0           !   periodic supercell repeat characteristic length ~ grain boundary size
        !real(kind=real64)               ::      u2x = 0.0d0             !   Debye-Waller vibration length
        real(kind=real64)               ::      rmin = 1.0d0 , rmax = 10.0d0
        real(kind=real64)               ::      qmin , qmax             !   minimum q-vector length (1/A) = 2 pi/a
        integer                         ::      nq = 100                 !   number of q subdivisions to take
!        integer                         ::      ng = 100                !   number of geodesic subdivisions to take
!        logical                         ::      opGeodesic = .false.
        logical                         ::      vox = .true.
        logical                         ::      soas = .false.
        logical                         ::      soas_linint = .false.
        
        real(kind=real64)               ::      a0 = 1.0d0              !   voxel cell size
        
        real(kind=real64),dimension(3)  ::      lineDir = (/0,0,1/)
        
        
    !---    physically meaningful parameters deduced from input file       
        real(kind=real64),dimension(3,3)    ::      a_super 
        type(XYZFile)                       ::      xyz 
        real(kind=real64),dimension(:,:),pointer    ::      colp
        integer                         ::      nAtoms
        
!        type(TriangulatedSurface)       ::      geodesic
         real(Kind=real64)               ::      qmod            !   modulus of q vector
        
!        complex(kind=real64),dimension(:,:),allocatable ::      pbc_q
!        complex(kind=real64),dimension(:,:),allocatable ::      dw_q
!        complex(kind=real64),dimension(:,:),allocatable ::      atom_q
!        real(kind=real64),dimension(:,:),allocatable    ::      unitq           !   unit q vector components
!        real(kind=real64),dimension(:,:),allocatable    ::      pos_private
!       complex(kind=real64),dimension(:,:),allocatable ::      q_private
!        real(kind=real64),dimension(:),allocatable      ::      Iq
        
        
    !---    dummy parameter     
        
        real(kind=real64)               ::      dd,i2s2,ff , aa,aa_sum !, zz_r,zz_i,sqrtu2x
        integer                         ::      ii,jj ,kk,ll , mm,nn   ! , i_progress  , hh
        integer                         ::      ix,iy,iz
        real(kind=real64),dimension(3)  ::      rr
        complex(kind=real64)            ::      zz , ss!, s1,s2,s3
        logical                         ::      ok
        
!        integer,dimension(3)            ::      vertex
!        real(kind=real64),dimension(3,3)    ::      tri
        
        
        integer                         ::      mx,my,mz
        real(kind=real64),dimension(:,:,:),allocatable  ::      vox_density
        real(kind=real64),dimension(:),allocatable      ::      rpf
        real(kind=real64)               ::      rpfsum,qrpfbar,q2rpfbar 
        real(kind=real64),dimension(3,3)    ::      a_cell,r_cell,ir_cell
        
        
        type(SimpleSupercell)               ::      super,rot_super
        real(kind=real64),dimension(3,27)   ::      Ruvw
        real(kind=real64),dimension(:),allocatable      ::      slice , sliceDensity 
        real(kind=real64),dimension(:),allocatable      ::      sliceFT
        
        
        
        
        integer,parameter               ::      T_INPUT = 1
        integer,parameter               ::      T_VOX = 2
        integer,parameter               ::      T_FFT = 3
        integer,parameter               ::      T_TOTAL = 4
        
        type(Callipers),dimension(T_TOTAL)    ::    timer
        
        
        
        
        timer(T_TOTAL) = Callipers_ctor()
        
    !---    read command line arguments
        cla = CommandLineArguments_ctor(20)  
         
        call setProgramDescription( cla, "LineProfile.exe" )
        call setProgramVersion( cla, VERSION )               
        
         
        
        call get( cla,"f",filename ,LIB_CLA_REQUIRED,"       input filename" )                                                        
        !call get( cla,"sigma",sigma ,LIB_CLA_OPTIONAL," grain size" )                                                        
        !call get( cla,"u2x",u2x     ,LIB_CLA_OPTIONAL,"   debye waller vibration variance (A^2)" )                                                        
        call get( cla,"rmin",rmin ,LIB_CLA_OPTIONAL,"  minimum real space separation considered" )                                                        
        call get( cla,"rmax",rmax ,LIB_CLA_OPTIONAL,"  maximum real space separation considered" )         
        qmin = 2*PI/rmax ; qmax = 2*PI/rmin                                               
        call get( cla,"qmin",qmin ,LIB_CLA_OPTIONAL,"  minimum q point considered" )                                                        
        call get( cla,"qmax",qmax ,LIB_CLA_OPTIONAL,"  maximum q point considered" )                                                        
        call get( cla,"nq",nq ,LIB_CLA_OPTIONAL,"    number of q points" )              
        call get( cla,"vox",vox,LIB_CLA_OPTIONAL,"   use voxelisation method" )
        if (.not. vox) then
            ii = 3 ; call get( cla,"k",lineDir,ii,LIB_CLA_REQUIRED,"     use 1d method with specified z direction" )
        else
            ii = 3 ; call get( cla,"k",lineDir,ii,LIB_CLA_OPTIONAL,"     use 1d method with specified z direction" )
        end if
         
        call get( cla,"soas",soas,LIB_CLA_OPTIONAL,"   use soas method" )
        call get( cla,"soas_linint",soas_linint,LIB_CLA_OPTIONAL,"   use soas method with linear interpolation" )
        if (soas_linint) soas = .true.
        call get( cla,"a0",a0,LIB_CLA_OPTIONAL,"    voxelisation method cell side" )
        call get( cla,"sr",SIGMA_RANGE ,LIB_CLA_OPTIONAL,"    voxelisation method sigma multiples" )   
!        call get( cla,"ng",ng ,LIB_CLA_OPTIONAL,"    number of angles per q point" )       
!        call get( cla,"opGeo",opGeodesic ,LIB_CLA_OPTIONAL," output geodesic of k points used" )                 
                                                   
           
!        if ( hasArgument(cla,"ng").or. hasArgument(cla,"opGeo") ) vox = .false.
!        if ( hasArgument(cla,"a0") ) vox = .true.
        
        if ( vox .and. hasArgument(cla,"k") ) then
            print *,"LineProfile.exe warning - conficting arguments. Choose -vox for 3d or -novox -k <x,y,z> for 1d"
            vox = .false.
        end if
        
        
        call report(cla)
        if (.not. allRequiredArgumentsSet(cla)) stop
        if (hasHelpArgument(cla)) stop
        call delete(cla)
        
        
    !---
        print *,""
        print *,"lineProfile.exe v"//trim(VERSION)
        print *,"^^^^^^^^^^^^^^^^^"//repeat("^",len_trim(VERSION))
        print *,""
    
        
        timer(T_INPUT) = Callipers_ctor()
        print *,""
        print *,"reading input file"
        print *,"^^^^^^^^^^^^^^^^^^"
        print *,""
        xyz = XYZFile_ctor(filename)
        call readHeader(xyz,ok)
        call input(xyz,verbose=.true.)
        call getSupercell(xyz,a_super,ok) 
        if (ok) then
            print *,"LineProfile.exe info - supercell read from file"
            print *,a_super(1,:)
            print *,a_super(2,:)
            print *,a_super(3,:)                                                                   
        else
            print *,"LineProfile.exe warning - supercell not read from file"                                        
        end if
        call report(xyz)        
        call getColumnsp(xyz,colp)
        nAtoms = getNAtoms(xyz)
        call pause( timer(T_INPUT) )
        
    !---    
    
        mx = nint( norm2( a_super(:,1) / a0 ) )               
        my = nint( norm2( a_super(:,2) / a0 ) )
        mz = nint( norm2( a_super(:,3) / a0 ) )
        a_cell(1:3,1) = a_super(1:3,1)/mx
        a_cell(1:3,2) = a_super(1:3,2)/my
        a_cell(1:3,3) = a_super(1:3,3)/mz
        sigma = a0/2
        
    
    
        if (vox) then
                         
            timer(T_VOX) = Callipers_ctor()
            print *,"LineProfile.exe info - using voxel method "
            print *,"LineProfile.exe info - voxel count ",mx,my,mz," broadening",sigma," range ",SIGMA_RANGE
            allocate( vox_density(0:mx-1,0:my-1,0:mz-1) )                      
            
            if (soas) then
                call SOAS_vox( nAtoms,colp, A_super,mx,my,mz, vox_density , soas_linint)
            else                             
                call voxeliseAtomicDensity( colp, sigma, a_super, vox_density )
            end if
            call delete(xyz)
            call pause(timer(T_VOX))
            
            
            !deallocate( colp )
            
            timer(T_FFT) = Callipers_ctor()
            allocate(rpf(0:nq))
            
            
            
            FFTW_DBG = .true.
            call radialPowerFunction( vox_density,qmin,qmax, rpf,rpfsum,qrpfbar,q2rpfbar,a_cell )
            
            print *,""
            write(*,fmt='(a16,a24)') "q (1/A)","I(q) (arb)"           
            do ii = 0,nq
                write(*,fmt='(f16.8,g24.8)') ( qmin + ii*(qmax-qmin)/nq ),rpf(ii)
            end do
            call pause(timer(T_FFT))
          
            
        else
        
            super = SimpleSupercell_ctor( a_cell,Mx,My,Mz )        
            call suggestSupercellOrientedWithN( super, lineDir , rot_super )
            
            r_cell = getSuperA(rot_super)            
            mz = ceiling( norm2(r_cell(:,3))/a0 )
            r_cell(:,3) = r_cell(:,3)/mz
            rot_super = SimpleSupercell_ctor( r_cell,1,1,Mz )   
            
            print *,"rotated supercell aligned with desired line direction"
            call report( rot_super )
            !r_cell = getSuperA(rot_super)
            
           ! mz = ceiling( norm2(r_cell(:,3))/a0 )
            print *,"LineProfile.exe info - using line method "
            print *,"LineProfile.exe info - slice count ",mz," broadening",sigma," range ",SIGMA_RANGE
            
            
        !---    find vectors to periodic replicas
            mm = 0
            do kk = -1,1
                do jj = -1,1
                    do ii = -1,1
                        mm = mm + 1
                        Ruvw(:,mm) = a_super(:,1)*ii + a_super(:,2)*jj + a_super(:,3)*kk
                        !print *,"Ruvw ",mm,ii,jj,kk,Ruvw(:,mm)
                    end do
                end do
            end do
                                                             
            
            
            
            call inverse3Mat( r_cell,ir_cell )
            dd = norm2(r_cell(:,3))                         !   slice thickness
            mm = ceiling( SIGMA_RANGE*sigma / dd )          !   number of slices to add for each atom
            
            allocate(sliceDensity(-mm:mm))
            allocate(slice(0:mz-1))
            !nq = int(mz/2) + 1
            allocate(sliceFT(0:nq))
            
            print *,"LineProfile.exe info - slice width ",dd," slices per atom ",(2*mm+1)
            
            
            do jj = 1,27
                rr(1:3) = Ruvw(1:3,jj)
                Ruvw(1:3,jj) = ir_cell(1:3,1)*rr(1) + ir_cell(1:3,2)*rr(2) + ir_cell(1:3,3)*rr(3) 
            end do
            
            
            do ii = 1,nAtoms
                rr(1:3) = colp(1:3,ii)
                colp(1:3,ii) = ir_cell(1:3,1)*rr(1) + ir_cell(1:3,2)*rr(2) + ir_cell(1:3,3)*rr(3) 
            end do
            i2s2 = 0.5d0
            
            
            slice = 0
            nn = 0
            do ii = 1,nAtoms
                
                do jj = 1,27
                    rr(1:3) = colp(1:3,ii) + Ruvw(1:3,jj)
                    
                !   find position of atom in slice space
                    ix = floor( rr(1) )
                    iy = floor( rr(2) )
                    iz = floor( rr(3) )
                    
                    
                    
                    if ( (ix*ix + iy*iy == 0) .and. (iz*(mz-1-iz)>=0) ) then
                        !   in the first replica of slice space
                        
                    !---    compute the gaussian projection on each slice
                        !ee = rr(3) - iz         !   offset of atom from 0th slice in z-direction                        
                        aa_sum = 0
                        do kk = -mm,mm                        
                            
                            ff = rr(3) - (iz+kk)                      !   offset of atom from kth slice in z-direction
                            
                            aa = exp( - ff*ff*i2s2 )
                            aa_sum = aa_sum + aa
                            sliceDensity(kk) = aa                            
                        end do
                        aa = 1/aa_sum
                        sliceDensity(-mm:mm) = sliceDensity(-mm:mm) * aa
                        
                        
                        !if (ii < 100) print *,"ii,jj,rr,ix,iy,iz ",ii,jj,rr,ix,iy,iz,aa_sum 
                                                
                    !---    add normalised gaussian projection to each slice                        
                        do kk = -mm,mm
                            ll = mod( iz + mz + kk , mz )
                            slice(ll) = slice(ll) + sliceDensity(kk)
                        end do     
                        nn = nn + 1                  
                          
                    end if
                end do
                    
                
            end do
            print *,"rotated box uses ",nn," atoms vs ",nAtoms," in file"
            slice = slice / nn
             
            
            sliceFT = 0.0d0
             
            
            
            
            do ii = 0,nq
                qmod = ( qmin + ii*(qmax-qmin)/nq )
                zz = cmplx( 0.0d0,-qmod*dd,kind=real64 ) 
                ss = 0
                do jj = 0,mz-1
                    ss = ss + slice(jj)*exp( zz*jj )
                end do
                sliceFT(ii) = real(ss)**2 + aimag(ss)**2
            end do
            
            
                                                              
            
            
            
            
            
          !  print *,""
          !  write(*,fmt='(a16,a24)') "z (Z)","f(z)"           
          !  do ii = 0,mz-1
          !      write(*,fmt='(f16.8,g24.8)') dd*ii,slice(ii)
          !  end do
            
            print *,""
            write(*,fmt='(a16,a24)') "q (1/A)","I(q) (arb)"           
            do ii = 0,nq
                write(*,fmt='(f16.8,g24.8)') ( qmin + ii*(qmax-qmin)/nq ),sliceFT(ii)
            end do
                
           
!            stop
        end if
        
        
        
        
        print *,""        
        print *,"timing data:"
        print *,"   read file ",elapsed(timer(T_INPUT))
        print *,"   vox       ",elapsed(timer(T_VOX))
        print *,"   fft       ",elapsed(timer(T_FFT))
        print *,"   total     ",elapsed(timer(T_TOTAL))
        print *,""        
        
        print *,""
        print *,"done"
        print *,""         
        
    !---
    
    
    
    
    
    
    
    
!    
!    
!        !   compute a geodesic of q-points
!        
!        geodesic = GeodesicGrid(ng*2)
!        call report(geodesic)
!        ng = getN_node(geodesic)
!        call init_random_seed( 12345 )
!        dd = 0.01 * 3.14159 / sqrt(real(ng))            !   1% * PI / points
!        allocate(unitq(3,ng))
!        if (opGeodesic) then
!            print *,"distorted geodesic vectors dd = ",dd
!            print *,""
!            print *,""
!            print *,ng
!            print *,"Lattice=""2 0 0 0 2 0 0 0 2"" Properties=Species:S:1:Pos:R:3"            
!        end if
!        do ii = 1,ng
!            rr = getNode(geodesic,ii)
!            unitq(1:3,ii) = rotateVector( rndRotMat(dd),rr )
!            if (opGeodesic) print *,"Du ",unitq(1:3,ii)
!        end do
!        !call delete(geodesic)
!        if (opGeodesic) print *,""
!        
!        
!        allocate(q_private(0:nq,ng))
!        
!        
!        !   compute effect of periodic supercells on broadening q-vector
!        allocate(pbc_q(0:nq,ng))
!        
!        
!        
!        
!        if (sigma == 0) then
!            pbc_q = ONE
!        else
!            
!            
!            pbc_q = ZERO
!!$OMP PARALLEL DEFAULT(PRIVATE) , SHARED(a_super,pbc_q,sigma,nq,ng,qmin,qmax,unitq,i_progress)            
!            
!            q_private = ZERO
!            dd = maxval( (/norm2(a_super(:,1)),norm2(a_super(:,2)),norm2(a_super(:,3))/) )          !   longest supercell repeat vector
!            mm = ceiling( 2*sigma/dd )                                                                 !   number of periodic cell repeats to consider
!            if (OMP_GET_THREAD_NUM()==0) print *,"LineProfile.exe info - max supercell repeat length ",dd," sigma = ",sigma," [hkl] repeats ",mm
!            i2s2 = 1/(2*sigma*sigma)
!            i_progress = 0
!        !$OMP DO    
!            do ii = 1,ng  
!            
!                i_progress = i_progress + 1                          
!                if ((OMP_GET_THREAD_NUM()==0).or.(i_progress > 1)) call progressBar( i_progress,ng )
!                do ll = -mm,mm
!                    do kk = -mm,mm                                                                           
!                        do hh = -mm,mm
!                            rr(1:3) = hh*a_super(1:3,1) + kk*a_super(1:3,1) + ll*a_super(1:3,3) 
!                            dd = rr(1)*rr(1) + rr(2)*rr(2) + rr(3)*rr(3) 
!                            dd = dd*i2s2
!                            if (dd>2.0d0) cycle                 !   compute in spherical region out to 2 sigma only.
!                            
!                            ff = unitq(1,ii)*rr(1) + unitq(2,ii)*rr(2) + unitq(3,ii)*rr(3)       !   dot product q.R_hkl/|q|
!                            
!                            
!                            do jj = 0,nq
!                                qmod = qmin + jj*(qmax-qmin)/nq                            
!                                ee = exp( - qmod*dd )        !       !   gaussian weight of this replica
!                                                                
!                                zz = Exp( -I * qmod*ff )
!                                
!                                q_private(jj,ii) = q_private(jj,ii) + ee*zz
!                            end do
!                            
!                        end do
!                    end do
!                end do
!            end do   
!        !$OMP END DO    
!            
!        !$OMP CRITICAL
!            pbc_q(0:nq,1:ng) = pbc_q(0:nq,1:ng) + q_private(0:nq,1:ng)
!        !$OMP END CRITICAL  
!            
!!$OMP END PARALLEL   
!            print *,""                  
!        end if
!       
!        
!    !---    compute Debye-Waller contribution to broadening
!        allocate(dw_q(0:nq,ng))
!        if (u2x == 0) then
!            dw_q = ONE
!        else
!            
!            dw_q = ZERO                                              
!!$OMP PARALLEL DEFAULT(PRIVATE) , SHARED(a_super,dw_q,u2x,nq,ng,qmin,qmax,unitq,i_progress)    
!
!            q_private = ZERO   
!            mm = 10                                                                !   number of position subdivisions to consider
!            if (OMP_GET_THREAD_NUM()==0) print *,"LineProfile.exe info - Debye-Waller variance u2x = ",u2x,"(A^2) repeats ",mm
!            i2s2 = 9.0d0/(2*mm*mm)  !   so that at [hkl]=[m00] have factor ( -9/2 )
!            sqrtu2x = sqrt(u2x)
!            i_progress = 0
!            
!        !$OMP DO              
!            do ii = 1,ng
!                
!                i_progress = i_progress + 1                          
!                if ((OMP_GET_THREAD_NUM()==0).or.(i_progress > 1)) call progressBar( i_progress,ng )
!
!                do ll = -mm,mm
!                    do kk = -mm,mm
!                        do hh = -mm,mm
!                            
!                            dd = hh*hh + kk*kk + ll*ll
!                            dd = dd*i2s2
!                            if (dd>4.5d0) cycle                 !   compute in spherical region out to 3 sigma only.
!                            
!                            ff = unitq(1,ii)*hh + unitq(2,ii)*kk + unitq(3,ii)*ll       !   dot product q.u/|q|
!                            ff = ff*sqrtu2x
!                            
!                            do jj = 0,nq
!                                qmod = qmin + jj*(qmax-qmin)/nq                            
!                                ee = exp( - qmod*dd )        !       !   gaussian weight of this position
!                                                                
!                                zz = Exp( -I * qmod*ff )
!                                
!                                q_private(jj,ii) = q_private(jj,ii) + ee*zz
!                            end do
!                        end do
!                    end do
!                end do
!            end do  
!            
!        !$OMP END DO    
!            
!        !$OMP CRITICAL
!            dw_q(0:nq,1:ng) = dw_q(0:nq,1:ng) + q_private(0:nq,1:ng)
!        !$OMP END CRITICAL  
!          
!!$OMP END PARALLEL      
!            print *,""                        
!        end if
!        
!        
!        
!        
!        
!        
!        
!    !---    compute atom contribution to broadening        
!        print *,"LineProfile.exe info - compute atom contribution to broadening"
!        allocate(atom_q(0:nq,ng))
!        !allocate(pos_private(3,nAtoms))         
!        atom_q = ZERO
!       ! pos_private(1:3,1:nAtoms) = colp(1:3,1:nAtoms)  
!!$OMP PARALLEL DEFAULT(PRIVATE) , SHARED(atom_q,nq,ng,qmin,qmax,unitq,colp,nAtoms,i_progress)    
!            
!         
!        q_private = ZERO
!        i_progress = 0
!    !$OMP DO    
!        do ii = 1,ng
!                        
!            i_progress = i_progress + 1                          
!            if ((OMP_GET_THREAD_NUM()==0).or.(i_progress > 1)) call progressBar( i_progress,ng )
!              
!            do kk = 1,nAtoms
!                                                    
!                ff = unitq(1,ii)*colp(1,kk) + unitq(2,ii)*colp(2,kk) + unitq(3,ii)*colp(3,kk)       !   dot product q.x/|q|
!                    do jj = 0,nq
!                        qmod = qmin + jj*(qmax-qmin)/nq                            
!                                                         
!                        !zz = Exp( -I * qmod*ff )
!                        
!                        zz_r = cos( qmod*ff )
!                        zz_i = sin( qmod*ff )
!                        
!                        q_private(jj,ii) = q_private(jj,ii) + cmplx( zz_r,zz_i,kind=real64 )
!                    end do
!            end do                    
!                
! 
!        end do             
!    !$OMP END DO    
!    
!        
!            
!    !$OMP CRITICAL
!        atom_q(0:nq,1:ng) = atom_q(0:nq,1:ng) + q_private(0:nq,1:ng)
!    !$OMP END CRITICAL  
!               
!     
!!$OMP END PARALLEL                              
!        
!           
!        print *,""
!        print *,"LineProfile.exe info - combining result"
!    !---    add everything together and find the result
!    !       S(q) = sum_hkljm p(R_hkl) g(x_m) Exp[ -i q.(r_j+x_m+R_hkl) ]
!    !            = (sum_hkl p(R_hkl) Exp [ - iq.R_hkl ] ) (sum_m g(x_m) Exp[ -i q.x_m ]) (sum_j Exp[ -i q.r_j ])
!    !            = ^------- pbc_q ----------------------^ ^---------- dw_q ------------^ ^-----atom_q----------^
!        
!        
!        dd = 1.0d0/(nAtoms*ng*ng)
!        
!
!        do jj = 0,nq
!                               
!            zz = 0.0d0
!            do ii = 1,ng
!                zz =  pbc_q(jj,ii)*dw_q(jj,ii)*atom_q(jj,ii)
!                atom_q(jj,ii) = zz*dd
!            end do
!                              
!         end do         
!                    
!        
!        allocate(Iq(0:nq))
!        Iq = 0
!        aa_sum = 0.0d0
!        do ii = 1,getn_triangle( geodesic )
!            call getTriangleInfo(geodesic,ii,tri,vertex )
!            aa = findArea(tri)
!            aa_sum = aa_sum + aa
!            
!            do jj = 0,nq
!                s1 = atom_q(jj,vertex(1))
!                s2 = atom_q(jj,vertex(2))
!                s3 = atom_q(jj,vertex(3))
!                zz = s1 + s2 + s3
!                Iq(jj) = Iq(jj) + real( (s1*conjg(s1) + s2*conjg(s2) + s3*conjg(s3) + zz*conjg(zz))*aa , kind=real64 )
!            end do    
!                
!                
!        end do
!        Iq = Iq / (12 * aa_sum)
!        print *,"triangulated area sum /4 Pi = ",aa_sum/(4*3.141592654d0)
!        
!        print *,""
!        write (*,fmt='(3a16)') "q","r=2 pi/q","I(q)"
!        do jj = 0,nq
!                    
!            qmod = qmin + jj*(qmax-qmin)/nq                            
!            
!            write (*,fmt='(3g16.4)') qmod,(2*PI/qmod),Iq(jj)
!         end do         
!            
!            
!    !---    
!       print *,""
!       print *,"done"
!       print *,""
!       
         
    contains
!---^^^^^^^^



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
        
        
        subroutine SOAS_vox( nAtoms,x, A_super,Nx,Ny,Nz, rho , linint)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      given the atom positions x(1:3,1:nAtoms)
    !*      and the periodic supercell A_super(1:3,1:3)
    !*      construct a density voxelisation rho(0:Nx-1,0:Ny-1,0:Nz-1)
    !*      with a Gaussian smoothing kernel width sigma range 3 sigma
    !*      sigma is chosen to be half the smallest voxel side length
    !*      
    !*      if linint then performs a trilinear interpolation
    !*
            integer,intent(in)                                  ::      nAtoms
            real(kind=real64),dimension(:,:),intent(in)         ::      x
            real(kind=real64),dimension(3,3),intent(in)         ::      A_super
            integer,intent(in)                                  ::      Nx,Ny,Nz
            real(kind=real64),dimension(0:,0:,0:),intent(out)   ::      rho
            logical,intent(in)                                  ::      linint
            real(kind=real64)                   ::      sigma
            real(kind=real64),dimension(3,3)    ::      A_cell,iA_cell 
            real(kind=real64),dimension(3,3)    ::      A_fine,iA_fine
            
            
            
        !---    the acceleration here is given by the way the kernel is constructed
        !       an atom is found to be between voxels ix:ix+1,iy:iy+1,iz:iz+1
        !       and so it contributes to nodes ix-1:ix+2 etc ( range = 3 sigma and sigma = cell/2 )
        !       we find kernel(-1:2,-1:2,-1:2) which gives the contribution on each node
        !       now here's the clever bit:
        !       I know that I am doing a Gaussian smoothing, so I don't really need to be that accurate with the offset of the atom in the kernel
        !       In fact, it probably suffices to find the position of the atom within 1/10th of the cell side.
        !       So I can _precompute_ kernel(-1:2,-1:2,-1:2, :,:,: )
        !   
        !       I have the option of computing the kernel assuming the atom is in the centre of the fine-spaced voxel grid ( linint=.false. )
        !       or at the nodes of a fine-spaced grid ( linint=.true. ), then taking a trilinear interpolation of the weighting of the atom at each node.
        !
        
            integer,parameter           ::      NFINE = 10
            real(kind=real64)           ::      iNFINE
            
            real(kind=real64),allocatable    ::      kernel(:,:,: , :,:,:)
            real(kind=real64)           ::      kernel_linint(-1:2,-1:2,-1:2)
            real(kind=real64)           ::      xx,yy,zz , dx,dy,dz
            
            real(kind=real64)           ::      i2s2 , dd , ddsum
            
            integer                     ::      ix,iy,iz
            integer                     ::      jx,jy,jz
            integer                     ::      kx,ky,kz
            integer                     ::      lx,ly,lz
            integer                     ::      ii
            
            
             
            
        !---    find the unit cell for the voxels and for the fine-mesh subdivision of the voxels                
            A_cell(1:3,1) = A_super(1:3,1)/Nx
            A_cell(1:3,2) = A_super(1:3,2)/Ny
            A_cell(1:3,3) = A_super(1:3,3)/Nz
            
            call inverse3Mat(A_cell,iA_cell)
                  
            iNFINE = 1.0d0/NFINE
            A_fine = A_cell * iNFINE
            iA_fine = iA_cell * NFINE                   !   fine spacing cell
            
            sigma = minval( (/ norm2( A_cell(1:3,1) ),norm2( A_cell(1:3,2) ),norm2( A_cell(1:3,3) ) /) )/2      
            i2s2 = 1/(2*sigma*sigma)    
            
            
            
        !---    construct the kernel        
            if (linint) then
                print *,"lineProfile::SOAS_vox info - constructing kernels for linear interpolation"
            !---   assume kernel fine space points are on fine space nodes
                allocate( kernel(-1:2,-1:2,-1:2,0:NFINE,0:NFINE,0:NFINE) )
                
    
                do iz = 0,NFINE                
                    do iy = 0,NFINE
                        do ix = 0,NFINE
                        
                            call progressBar( 1+ix + (NFINE+1)*(iy + (NFINE+1)*iz) , (NFINE+1)**3 )
                        
                        !---    construct position of point at node in fine cell
                            xx = A_fine(1,1)*ix + A_fine(1,2)*iy + A_fine(1,3)*iz
                            yy = A_fine(2,1)*ix + A_fine(2,2)*iy + A_fine(2,3)*iz
                            zz = A_fine(3,1)*ix + A_fine(3,2)*iy + A_fine(3,3)*iz
                                
                            ddsum = 0.0d0
                            
                            do jz = -1,2
                                do jy = -1,2
                                    do jx = -1,2
                                    
                                    !---    construct separation of node from point in fine cell
                                        dx = A_cell(1,1)*jx + A_cell(1,2)*jy + A_cell(1,3)*jz - xx
                                        dy = A_cell(2,1)*jx + A_cell(2,2)*jy + A_cell(2,3)*jz - yy
                                        dz = A_cell(3,1)*jx + A_cell(3,2)*jy + A_cell(3,3)*jz - zz
                            
                                       
                                        dd = dx*dx + dy*dy + dz*dz     
                                        dd = dd*i2s2
                                        if (dd>4.5d0) then  !   3 sigma condition:  r^/( 2 sigma^2 ) = 4.5 when r = 3 sigma 
                                            dd = 0
                                        else
                                            dd = Exp( -dd )                                   
                                        end if
                                        kernel( jx,jy,jz , ix,iy,iz ) = dd
                                        ddsum = ddsum + dd        
                                        
                                    end do
                                end do
                            end do      !   node offset
    
                        !---    normalise                        
                            ddsum = 1/ddsum
                            kernel( :,:,: , ix,iy,iz ) = kernel( :,:,: , ix,iy,iz ) * ddsum
                            
                            
                        end do
                    end do
                end do                 
                
            else        !   not linint
            
                print *,"lineProfile::SOAS_vox info - constructing kernel"    
            !---   assume kernel fine space points are in centre of fine space voxels
                allocate( kernel(-1:2,-1:2,-1:2,0:NFINE-1,0:NFINE-1,0:NFINE-1) )
                    
                do iz = 0,NFINE-1                
                    do iy = 0,NFINE-1
                        do ix = 0,NFINE-1
                        
                            call progressBar( 1+ix + NFINE*(iy + NFINE*iz) , NFINE**3 ) 
                        
                        !---    construct position of point at centre of fine cell
                            xx = A_fine(1,1)*(ix+0.5d0) + A_fine(1,2)*(iy+0.5d0) + A_fine(1,3)*(iz+0.5d0)
                            yy = A_fine(2,1)*(ix+0.5d0) + A_fine(2,2)*(iy+0.5d0) + A_fine(2,3)*(iz+0.5d0)
                            zz = A_fine(3,1)*(ix+0.5d0) + A_fine(3,2)*(iy+0.5d0) + A_fine(3,3)*(iz+0.5d0)
                                
                            ddsum = 0.0d0
                            
                            do jz = -1,2
                                do jy = -1,2
                                    do jx = -1,2
                                    
                                    !---    construct separation of node from point in fine cell
                                        dx = A_cell(1,1)*jx + A_cell(1,2)*jy + A_cell(1,3)*jz - xx
                                        dy = A_cell(2,1)*jx + A_cell(2,2)*jy + A_cell(2,3)*jz - yy
                                        dz = A_cell(3,1)*jx + A_cell(3,2)*jy + A_cell(3,3)*jz - zz
                            
                                       
                                        dd = dx*dx + dy*dy + dz*dz     
                                        dd = dd*i2s2
                                        if (dd>4.5d0) then  !   3 sigma condition:  r^/( 2 sigma^2 ) = 4.5 when r = 3 sigma 
                                            dd = 0
                                        else
                                            dd = Exp( -dd )                                   
                                        end if
                                        kernel( jx,jy,jz , ix,iy,iz ) = dd
                                        ddsum = ddsum + dd        
                                        
                                    end do
                                end do
                            end do      !   node offset
    
                        !---    normalise                        
                            ddsum = 1/ddsum
                            kernel( :,:,: , ix,iy,iz ) = kernel( :,:,: , ix,iy,iz ) * ddsum
                            
                            
                        end do
                    end do
                end do       
            end if
          
        !---    kernel done
        
        
        
        
        !---    now loop through each atom. 
            print *,"lineProfile::SOAS_vox info - voxelising"    
            rho = 0.0
            do ii = 1,nAtoms
                 
                call progressBar( ii,nAtoms )       
            
            !---    find the position of the atom in voxel space ( should be 0:Nx-1, but the atom may be outside the periodic supercell bounds. Will fix below ) 
                xx = iA_cell(1,1)*x(1,ii) + iA_cell(1,2)*x(2,ii) + iA_cell(1,3)*x(3,ii)
                yy = iA_cell(2,1)*x(1,ii) + iA_cell(2,2)*x(2,ii) + iA_cell(2,3)*x(3,ii)
                zz = iA_cell(3,1)*x(1,ii) + iA_cell(3,2)*x(2,ii) + iA_cell(3,3)*x(3,ii)
                
            
            !---    find which voxel the atom sits in. Note: haven't yet made sure this cell is within range...
                jx = floor( xx )                        ! eg xx = 12.3456 -> jx = 12 )
                jy = floor( yy )  
                jz = floor( zz )  
            
            
            !---    now find which fine voxel the atom sits in. 
                ix = floor( (xx - jx)*NFINE )           !  xx = 12.3456 -> ix = 3   ( if NFINE = 10 )
                iy = floor( (yy - jy)*NFINE )
                iz = floor( (zz - jz)*NFINE )
                
                
                if (linint) then
                
                !   I know the atom is somewhere between fine-spaced nodes (ix:ix+1 , iy+iy+1 , iz+iz+1)
                !   find the weighting on each and construct a tri-linear interpolation of the kernel 
                    xx = (xx - jx)*NFINE - ix           ! xx = 12.3456 -> xx = 0.456
                    yy = (yy - jy)*NFINE - iy
                    zz = (zz - jz)*NFINE - iz
                               
                    kernel_linint( :,:,: ) = kernel( :,:,: , ix  ,iy  ,iz   )*(1-xx)*(1-yy)*(1-zz)     &
                                           + kernel( :,:,: , ix+1,iy  ,iz   )*(  xx)*(1-yy)*(1-zz)     &
                                           + kernel( :,:,: , ix  ,iy+1,iz   )*(1-xx)*(  yy)*(1-zz)     &
                                           + kernel( :,:,: , ix+1,iy+1,iz   )*(  xx)*(  yy)*(1-zz)     &
                                           + kernel( :,:,: , ix  ,iy  ,iz+1 )*(1-xx)*(1-yy)*(  zz)     &
                                           + kernel( :,:,: , ix+1,iy  ,iz+1 )*(  xx)*(1-yy)*(  zz)     &
                                           + kernel( :,:,: , ix  ,iy+1,iz+1 )*(1-xx)*(  yy)*(  zz)     &
                                           + kernel( :,:,: , ix+1,iy+1,iz+1 )*(  xx)*(  yy)*(  zz)      
                    
                              
                !---    now can add kernel            
                    do kz = -1,2
                        lz = mod( jz + kz + Nz , Nz )       !   now this is within 1st periodic replica
                        do ky = -1,2
                            ly = mod( jy + ky + Ny , Ny )        
                            do kx = -1,2
                                lx = mod( jx + kx + Nx , Nx )  
                                
                                rho(lx,ly,lz) = rho(lx,ly,lz) + kernel_linint( kx,ky,kz )
                                
                            end do
                        end do
                    end do    
                    
                
                else
                    
                !---    now can add kernel assuming the atom is in the centre of fine spaced voxel (ix,iy,iz)    
                    do kz = -1,2
                        lz = mod( jz + kz + Nz , Nz )       !   now this is within 1st periodic replica
                        do ky = -1,2
                            ly = mod( jy + ky + Ny , Ny )        
                            do kx = -1,2
                                lx = mod( jx + kx + Nx , Nx )  
                                
                                rho(lx,ly,lz) = rho(lx,ly,lz) + kernel( kx,ky,kz , ix,iy,iz )
                                
                            end do
                        end do
                    end do    
                    
                end if                
                
            end do     
            
            return
        end subroutine SOAS_vox                   
                                 
        
    end program LineProfile    