
    module Lib_DiffractionConditions
!---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
!*      This module is intended to find "good" two-beam diffraction conditions, 
!*      based on a request for a particular g-vector and zone axis from the user.
!*

        use iso_fortran_env
#ifdef MPI           
        use mpi_f08 
#endif

        use Lib_RelativisticElectrons
        use Lib_ReadExtinctionDistances
        use Lib_RotationMatrices 
        use Lib_Lattices
        use Lib_GoldenSection
        implicit none
        private
        
        
        real(kind=real64),private,parameter         ::      PI = 3.141592653589790d0  

        integer,public,parameter                    ::      NTHETA  = 2000
        real(kind=real64),public                    ::      RHO_MAX = 10.0d0
        real(kind=real64),public                    ::      THETA_MAX = 5.0d0  
        integer(kind=int64),private,parameter       ::      BADF00D = int( z'BADF00D',kind=int64 )        
        real(kind=real64),private,parameter         ::      SG_UNSET = transfer( (BADF00D+ishft(BADF00D,32_int64)),1.0d0 )
        real(kind=real64),private,parameter         ::      DEGTORAD = 3.141592654d0/180

        public              ::      rotateToStandardFrame
        public              ::      findBestRotationMatrix
        public              ::      tiltToDarkField
        public              ::      tiltToBrightField
        public              ::      intensityAtDepthZ
        public              ::      energy_deviation_parameter
        public              ::      energyToLengthDeviationParameter
        public              ::      findCrystalliteSideLengths
        public              ::      findTransformedSupercell
        public              ::      getReciprocalLatticeVectors
        public              ::      reflection
        
        interface   scoreForIncidentVector
            module procedure    scoreForIncidentVector1
            module procedure    scoreForIncidentVector2
        end interface
        
        interface   reflection
            module procedure    reflection0
            module procedure    reflection1            
        end interface
        
        interface   rotateToStandardFrame
            module procedure    rotateToStandardFrame0
            module procedure    rotateToStandardFrame1            
        end interface
        
        
        
        
    contains
!---^^^^^^^^


        subroutine rotateToStandardFrame0( z0,A_cell, hkl ,U )
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ 
    !*      given a zone axis to aim for z0 , and the lattice vectors describing the conventional unit cell A_cell
    !*      and a target reflection to hit bragg condition hkl
    !*      return the best rotation matrix U which puts the zone axis along (0,0,1) and the g-vector along (1,0,0)
            
            real(kind=real64),dimension(3),intent(in)       ::      z0
            real(kind=real64),dimension(3,3),intent(in)     ::      A_cell              !   lattice vectors of conventional unit cell
            real(kind=real64),dimension(3),intent(in)       ::      hkl
            real(kind=real64),dimension(3,3),intent(out)    ::      U
            
            
            real(kind=real64),dimension(3,3)        ::      BB      !   reciprocal lattice vectors
            
            real(kind=real64),dimension(3)          ::      zz,gg,xx,yy
              
            
        !---    find the reciprocal lattice vectors in the lab frame _before_ rotation
        !       I want to use the definition that the columns of B are the reciprocal lattice vectors.
            BB = getReciprocalLatticeVectors( A_cell )
            !call inverse3Mat(transpose(A_cell),BB) ; BB = 2*PI*BB  
            
        !       ... and from this find the g-vector corresponding to [hkl]             
        !       this will be consistent with my choice for reciprocal lattice vectors.
            gg = reflection( BB, hkl )
             
        !   find the normalised zone axis in lab frame
            zz = matmul( A_cell,z0 )  
            zz(1:3) = zz(1:3)/norm2( zz )
            
            
        !---    find the rotation matrix which puts the desired zone axis along (0,0,1) and the g-vector along (1,0,0)
            xx = gg(1:3) - dot_product(gg,zz) * zz(1:3)                 
            xx = xx / norm2(xx)         
            yy = cross_product(zz,xx)
        
        
           U(1,1:3) = xx(1:3)
           U(2,1:3) = yy(1:3)
           U(3,1:3) = zz(1:3)
            
        !    
            return
        end subroutine rotateToStandardFrame0    
            

        subroutine rotateToStandardFrame1( z0,A_cell, hkl ,U )
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ 
    !*      given a zone axis to aim for z0 , and the lattice vectors describing the conventional unit cell A_cell
    !*      and a target reflection to hit bragg condition hkl
    !*      return the best rotation matrix U which puts the zone axis along (0,0,1) and the g-vector along (1,0,0)
            
            real(kind=real64),dimension(3),intent(in)       ::      z0
            real(kind=real64),dimension(3,3),intent(in)     ::      A_cell              !   lattice vectors of conventional unit cell
            integer,dimension(3),intent(in)                 ::      hkl
            real(kind=real64),dimension(3,3),intent(out)    ::      U
            
            call rotateToStandardFrame0( z0,A_cell, real(hkl,kind=real64) ,U )
            
            return
        end subroutine rotateToStandardFrame1    
            
        
        
        subroutine findBestRotationMatrix( v,z0,A_cell,latt, hkl ,U, R , extinction_distances_filename,foil_thickness, nNearBragg,hkl_nearBragg,eps_nearBragg )
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ 
    !*      given the electron velocity v ( in A/fs )
    !*      a zone axis to aim for z0 , and the lattice vectors describing the conventional unit cell A_cell
    !*      and a target reflection to hit bragg condition hkl
    !*      and a maximum tilt angle permitted
    !*      and a lattice type ( so that permitted reflections are found correctly )
    !*      and reciprocal lattice vectors B
    !*      and a small energy parameter gamma to prevent scores blowing up
    !*      and a maximum distance in reflection space to search rho
    !*      return the best rotation matrix U which puts the zone axis along (0,0,1) and the g-vector along (1,0,0)
    !*      and the tilt-stage tweak R which rotates the crystal into a good two-beam condition.
    !*      which has the greatest energy deviation for all other reflections
    
    !*      optionally return a list of those reflections which are near the Bragg condition
            
            real(kind=real64),intent(in)                    ::      v
            real(kind=real64),dimension(3),intent(in)       ::      z0
            real(kind=real64),dimension(3,3),intent(in)     ::      A_cell              !   lattice vectors of conventional unit cell
            type(Lattice),intent(in)                        ::      latt 
            integer,dimension(3),intent(in)                 ::      hkl
            real(kind=real64),dimension(3,3),intent(out)    ::      U,R
            character(len=*),intent(in)                     ::      extinction_distances_filename
            real(kind=real64),intent(in)                    ::      foil_thickness
            integer,intent(out),optional                    ::      nNearBragg
            integer,dimension(:,:),allocatable,intent(out),optional         ::      hkl_nearBragg
            real(kind=real64),dimension(:),allocatable,intent(out),optional ::      eps_nearBragg    
            
            
            
            
            real(kind=real64),dimension(3,3)        ::      BB      !   reciprocal lattice vectors
             
            real(kind=real64),dimension(3)          ::      zz,ghkl0,gg 
            real(kind=real64),dimension(:,:),allocatable    ::      khat
             
            integer                                 ::      nn      !   number of permitted reflections
            integer                                 ::      nk      !   number of permitted incident vectors
            real(kind=real64),dimension(:,:),allocatable    ::      permitted_ghkl
            real(kind=real64),dimension(:),allocatable      ::      xi_g
            integer,dimension(:,:),allocatable      ::      permitted_hkl
            
            real(kind=real64)                       ::      gamma
            real(kind=real64)                       ::      theta0, phi0
            
            real(kind=real64)                       ::      theta, phi
            real(kind=real64)                       ::      ss,ss_best,ee 
            real(kind=real64),dimension(:),allocatable      ::      ss_list
            integer                                 ::      ii,ii_best
!             character(len=16)                       ::      aaaa
            
            integer,dimension(3)                    ::      uvw_h,uvw_m,uvw_l
            logical                                 ::      ok
            
            
            
            
            integer                             ::      nProcs,rank,ierror
            
            call MPI_COMM_SIZE(MPI_COMM_WORLD, nProcs, ierror)
            call MPI_COMM_RANK(MPI_COMM_WORLD, rank, ierror)
            
            
            
            
            
            
            
        !---    find the reciprocal lattice vectors in the lab frame _before_ rotation
        !       I want to use the definition that the columns of B are the reciprocal lattice vectors.
            BB = getReciprocalLatticeVectors( A_cell )
            !call inverse3Mat(transpose(A_cell),BB) ; BB = 2*PI*BB  
            
        !       ... and from this find the g-vector corresponding to [hkl]             
        !       this will be consistent with my choice for reciprocal lattice vectors.
            gg = reflection( BB, hkl )
            
        !   find the normalised zone axis in lab frame
            zz = matmul( A_cell,z0 )  
            ss = 1/norm2( zz )
            zz(1:3) = zz(1:3) * ss
            
            call rotateToStandardFrame( z0,A_cell, hkl ,U ) 
            
        !---    find the g-vector in this new rotated frame
            ghkl0(1:3) = U(1:3,1)*gg(1) + U(1:3,2)*gg(2) + U(1:3,3)*gg(3) 
            if (rank==0) then
                write (*,fmt='(a)') "Lib_DiffractionConditions::findBestRotationMatrix info - zone axis, unit vec in lab frame, [hkl], g-vec in lab frame before, g-vec after"
                write (*,fmt='(f12.5,a,f12.5,a,i4,a,f12.5,a,f12.5)') z0(1),"    ",zz(1),"    ",hkl(1),"    ",gg(1),"    ",ghkl0(1)
                write (*,fmt='(f12.5,a,f12.5,a,i4,a,f12.5,a,f12.5)') z0(2),"    ",zz(2),"    ",hkl(2),"    ",gg(2),"    ",ghkl0(2)
                write (*,fmt='(f12.5,a,f12.5,a,i4,a,f12.5,a,f12.5)') z0(3),"    ",zz(3),"    ",hkl(3),"    ",gg(3),"    ",ghkl0(3)
            end if            
            
            
        !---    find all permitted g-vectors close to [000]            
            
            call find_gvectors_permittedReflections( latt,BB,hkl,RHO_MAX, nn,permitted_ghkl,permitted_hkl )
            if (rank==0) write (*,fmt='(a,i8,a,f12.3)') "Lib_DiffractionConditions::findBestRotationMatrix info - number of permitted reflections ",nn," max radius ", RHO_MAX
            
        !---    read extinction distances for all permitted g-vectors
            allocate(xi_g(nn))            
            call readExtinctionDistance( extinction_distances_filename,permitted_hkl,getSymmetry(latt),xi_g, ok )
            if (rank==0) write (*,fmt='(a,l2)') "Lib_DiffractionConditions::readExtinctionDistance info - ",ok
            
            
        !---    find all rotations which have ghkl0 at the Bragg condition...
            allocate( khat(3,4*NTHETA) )
                            
            theta0 = 0 ; phi0 = 0
            call findAllEulerAngles( v,ghkl0,theta0,phi0,THETA_MAX*DEGTORAD,khat,nk )
            if (nk == 0) then
                if (rank==0) write (*,fmt='(a)') "Lib_DiffractionConditions::findBestRotationMatrix warning - no tilt angles meet diffraction condition - trying double tilt..."
                call findAllEulerAngles( v,ghkl0,theta0,phi0,2*THETA_MAX*DEGTORAD,khat,nk )
                if (nk == 0) then
                    if (rank==0) write (*,fmt='(a)') "Lib_DiffractionConditions::findBestRotationMatrix WARNING - no tilt angles meet diffraction condition - trying no phi restriction..."
                    call findAllEulerAngles( v,ghkl0,theta0,phi0,THETA_MAX*DEGTORAD,khat,nk , nolimitPhi = .true.)
                end if
            end if
            if (rank==0) write (*,fmt='(a,i8,a,f12.3,a)') "Lib_DiffractionConditions::findBestRotationMatrix info - number of tilt angles to test ",nk," max theta ",THETA_MAX," (deg)"
            
        !---    note that testing a rotated g vector against a constant k vector is equivalent to a constant g-vector and rotating k.            
             
            
             
        !---    give each incident vector a score. Find the best
            allocate( ss_list(nk) ) 
            ss_best = huge(1.0) ; ii_best = 0
            do ii = 1,nk                                                !   for each possible incident vector ...
    !            ss = scoreForIncidentVector( v,khat(:,ii),gamma,permitted_ghkl )  !   ... find the score ... 
                ss = scoreForIncidentVector( v,khat(:,ii),permitted_ghkl , xi_g,foil_thickness )  !   ... find the score ... 
                ss_list(ii) = ss
                if (ss<ss_best) then
                    ss_best = ss ; ii_best = ii                         !   ... and find the best.
                end if
            end do
    
            deallocate( ss_list )
                        
        !---    v.0.1       - assume rotation theta about x then phi about z                    
            !theta = acos( khat(3,ii_best) ) ; phi = acos( khat(1,ii_best)/sin(theta) )
        !---    v.0.2       - rotation theta about x then phi about y        
        !    theta = -asin( khat(2,ii_best) ) ; phi = acos( khat(3,ii_best)/cos(theta) )
        !---    v.0.3           
        !   khat(:,nkhat+jj) = (/ -sinp, sint*cosp, cost*cosp /)
            phi = -asin( khat(1,ii_best) ) ; theta = acos( khat(3,ii_best)/cos(phi) )
             
              
            
            if (rank==0) then
                print *,"Lib_DiffractionConditions::findBestRotationMatrix info - best score ",ss_best!," zero deviation score ",score( 0.0d0,gamma )*gamma/PI
                print *,"Lib_DiffractionConditions::findBestRotationMatrix info - tilt angles"
                write(*,fmt='(a,f12.5,a,f12.5,a)') "   tilt about [100] theta     ",theta," (rad) = ",theta*180.0d0/PI," (deg)"
                write(*,fmt='(a,f12.5,a,f12.5,a)') "   tilt about [010] phi       ",phi," (rad) = ",phi*180.0d0/PI," (deg)"            
            end if             
         
         
            
            gamma = 2*PI*HBAR*v*0.01d0
            if (present(nNearBragg)) then
                if (rank==0) print *,"Lib_DiffractionConditions::findBestRotationMatrix info - Lorentzian broadening gamma (sg=0.01) = ",gamma," (eV)"            
            ss_best = 1/score( 0.0d0,gamma )
            
                nNearBragg = 0
                do ii = 1,nn
                    ee = energy_deviation_parameter( v,khat(:,ii_best),permitted_ghkl(:,ii) )
                    ss = score( ee,gamma )*ss_best
                    if ( (ii<=2).or.(ss>0.25d0).or.(scoreForIncidentVector( v,khat(:,ii_best),permitted_ghkl(:,ii), xi_g(1),xi_g(ii),foil_thickness )>0.1)  ) then
                        nNearBragg = nNearBragg+1
                    end if
                end do
                allocate(hkl_nearBragg(3,nNearBragg))
                allocate(eps_nearBragg(nNearBragg))
                call inverse3Mat(BB,R)      !   find inverse recip latt vecs to find reflection indices
                nNearBragg = 0
                do ii = 1,nn
                    ee = energy_deviation_parameter( v,khat(:,ii_best),permitted_ghkl(:,ii) )
                    ss = score( ee,gamma )*ss_best
                    if ( (ii<=2).or.(ss>0.25d0).or.(scoreForIncidentVector( v,khat(:,ii_best),permitted_ghkl(:,ii), xi_g(1),xi_g(ii),foil_thickness )>0.1) ) then
                        nNearBragg = nNearBragg+1
                        hkl_nearBragg(1:3,nNearBragg) = nint(matmul(R,permitted_ghkl(:,ii)))
                        eps_nearBragg(nNearBragg) = ee
                    end if
                end do
            end if
            
            
            
            
        !!---    recover the best incident vector direction as the z-axis for the rotation matrix....
        
            R(1,1:3) = (/  cos(phi) , sin(theta)*sin(phi) , cos(theta)*sin(phi) /)
            R(2,1:3) = (/  0.0d0    , cos(theta)          , -sin(theta)         /)
            R(3,1:3) = (/ -sin(phi) , sin(theta)*cos(phi) , cos(theta)*cos(phi) /)
           
             
            ee = energy_deviation_parameter( v,khat(:,ii_best),ghkl0 )
            if (rank==0) then
                if (abs(ee)<gamma * 1.0d-8) then
                    write(*,fmt='(a,f16.10,a)') " Lib_DiffractionConditions::findBestRotationMatrix info - deviation param best ",ee," (eV) "
                else if (abs(ee)<gamma) then
                    write(*,fmt='(a,f16.10,a)') " Lib_DiffractionConditions::findBestRotationMatrix warning - deviation param best ",ee," (eV)"
                else  
                    write(*,fmt='(a,f16.10,a)') " Lib_DiffractionConditions::findBestRotationMatrix ERROR - deviation param best ",ee," (eV)"
                end if
            end if            
            if (R(3,3)<0) then
                R(:,3) = -R(:,3)
                R(:,2) = -R(:,2)
            end if
            
            
            
            
!         !       ... and complete the basis 
!             
            if (rank==0) print *,"Lib_DiffractionConditions::findBestRotationMatrix info - best rotation matrix, [uvw] (low, med, high precision)"
            
            zz = R(3,1:3)                                           !   lab frame direction of k-vector
            zz = reflection( transpose(BB),zz )                     !   expressed in cell space
            zz = zz/norm2(zz)                                       !   normalised to unit vector                                    
            
            uvw_h(1:3) = nint( zz(1:3) * 256 )
            uvw_m(1:3) = nint( zz(1:3) * 16 )
            uvw_l(1:3) = nint( zz(1:3) * 1 )
            
            
            
            return
        end subroutine findBestRotationMatrix
!
        subroutine tiltToBrightField( v,A_cell, hkl , xi0,xig,foil_thickness,U, R )
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      given the electron velocity v ( in A/fs )
    !*      and the lattice vectors describing the conventional unit cell A_cell
    !*      and a g-vector reflection
    !*      and a lattice type ( so that permitted reflections are found correctly )
    !*      and the tilt-stage orientation R which gives good 2 beam conditions
    !*      and the tilt-stage tweak R which rotates the crystal into a good dark field condition
            
            real(kind=real64),intent(in)                    ::      v
            real(kind=real64),dimension(3,3),intent(in)     ::      A_cell              !   lattice vectors of conventional unit cell
            integer,dimension(3),intent(in)                 ::      hkl
            real(kind=real64),intent(in)                    ::      xi0,xig,foil_thickness
            real(kind=real64),dimension(3,3),intent(in)     ::      U
            real(kind=real64),dimension(3,3),intent(inout)  ::      R
            
            
                        
            
            real(kind=real64)               ::      bg,b0,sg , sgp 
            real(kind=real64)               ::      uu 
            real(kind=real64)               ::      cost,sint


            real(kind=real64),dimension(3,3)    ::      BB,RR
            real(kind=real64),dimension(3)      ::      kk,gg
            
            real(kind=real64)               ::      theta,lp,const,phiL2
 
            
            type(GoldenSection)             ::      gold
            real(kind=real64)               ::      x1,x2,x3,y1,y2,y3
            logical                         ::      isConverged,isWithinTimeLimit
        
!            character(len=16)                       ::      aaaa
            
            
            integer                             ::      nProcs,rank,ierror
            
            call MPI_COMM_SIZE(MPI_COMM_WORLD, nProcs, ierror)
            call MPI_COMM_RANK(MPI_COMM_WORLD, rank, ierror)
            
            
        !----   compute g-vector in rotated frame. This will be close to [100], but not identical to it.
            BB = getReciprocalLatticeVectors( A_cell )
            !call inverse3Mat(transpose(A_cell),BB) ; BB = 2*PI*BB                                              
            kk = reflection( BB, hkl )                                                                  
            RR = matmul( R,U )
            gg(1:3) = RR(1:3,1)*kk(1) + RR(1:3,2)*kk(2) + RR(1:3,3)*kk(3) 
            !print *,"gg = ",gg
            
            if (rank==0) then
                print *,"Lib_DiffractionConditions::tiltToBrightField() info - [hkl], g-vec in lab frame, g-vec after rotation"
                write (*,fmt='(i4,a,f12.5,a,f12.5)') hkl(1),"    ",kk(1),"    ",gg(1)
                write (*,fmt='(i4,a,f12.5,a,f12.5)') hkl(2),"    ",kk(2),"    ",gg(2)
                write (*,fmt='(i4,a,f12.5,a,f12.5)') hkl(3),"    ",kk(3),"    ",gg(3)
            end if
                                
        !---    propagation of electrons in perfect xtal
        !               
        !       d/dz (  phi_0 ) = PI i  ( b0    bg        )      !   with b0 = 1/xi0, bg = 1/xig
        !            (  phi_g )         ( bg    b0 + 2 sg )      !         
        !            
            bg = 1/xig
            b0 = 1/xi0
            uu = energy_deviation_parameter( v,(/0.0d0,0.0d0,1.0d0/),gg )
            sg = energyToLengthDeviationParameter(v,uu)
            phiL2 = (bg*sin( PI* sqrt( sg*sg + bg*bg )* foil_thickness ))**2 / ( sg*sg + bg*bg )          
            if (rank==0) then
                print *,"Lib_DiffractionConditions::tiltToBrightField info - deviation parameter (before) ",uu," (eV) ",sg," (1/A)"
                print *,"Lib_DiffractionConditions::tiltToBrightField info - foil thickness L  = ",foil_thickness  
                print *,"Lib_DiffractionConditions::tiltToBrightField info - | phi_g(L) |^2 (before)  ", (bg*sin( PI* sqrt( sg*sg + bg*bg )* foil_thickness ))**2 / ( sg*sg + bg*bg )
                print *,"Lib_DiffractionConditions::tiltToBrightField info - b0,bg,sg = ",b0,bg,sg
            end if            
            
            const = - HBAR*(gg(1)*gg(1)+gg(2)*gg(2)+gg(3)*gg(3))/(4*PI*ME*v)
            
            x1 = -THETA_MAX*0.01d0*DEGTORAD
            cost = cos(x1) ; sint = sin(x1)            
            Lp  = foil_thickness/cost
            sgp = ( sint*gg(1) - cost*gg(3) )/(2*PI)            
            y1 = -(bg*sin( PI* sqrt( sgp*sgp + bg*bg )*Lp ))**2 / ( sgp*sgp + bg*bg )
            
            x3 = +THETA_MAX*0.01d0*DEGTORAD
            cost = cos(x3) ; sint = sin(x3)            
            Lp  = foil_thickness/cost
            sgp = ( sint*gg(1) - cost*gg(3) )/(2*PI)            
            y3 = -(bg*sin( PI* sqrt( sgp*sgp + bg*bg )*Lp ))**2 / ( sgp*sgp + bg*bg )
            
                                   
            gold = GoldenSection_ctor( x1,x3, y1,y3)
            do
                call nextPoint(gold,x2)
                cost = cos(x2) ; sint = sin(x2)            
                Lp  = foil_thickness/cost
                sgp = ( sint*gg(1) - cost*gg(3) )/(2*PI)            
                y2 = -(bg*sin( PI* sqrt( sgp*sgp + bg*bg )*Lp ))**2 / ( sgp*sgp + bg*bg )
            
                call minimise(gold,x2,y2,isConverged,isWithinTimeLimit)
                if (isConverged) exit
                if (.not. isWithinTimeLimit) exit
            end do
            
            y2 = -y2
            if (y2 > phiL2) then
                theta = x2 ; cost = cos(x2) ; sint = sin(x2)  
            else
                theta = 0.0d0; cost = 1.0d0 ; sint = 0.0d0
            end if
           ! print *,"best rotation angle ",theta
            if (rank==0) print *,"Lib_DiffractionConditions::tiltToBrightField info - | phi_g(L) |^2 (after)  ",-y2
            
            
            RR(1,1:3) = (/ cost , 0.0d0 , sint /)
            RR(2,1:3) = (/ 0.0d0 , 1.0d0 , 0.0d0 /)
            RR(3,1:3) = (/-sint , 0.0d0 , cost /)
            
            
            
            
            
!             
!             
!             
!             
!         !---    solution at opposite edge of foil: 
!         !
!         !           | phi_g(L) |^2 = bg^2 Sin^2[  pi Sqrt[ sg^2 + bg^2 ] L ] / ( sg^2 + bg^2 )
!         !    
!         !
!         !       so dark field points are where 
!         !           Sqrt[ sg^2 + bg^2 ] L == n 
!         !   
!         !       bright field points where
!         !           Sqrt[ sg^2 + bg^2 ] L == n + 1/2 == n'
!         !    
!         
!             print *,"dbg    sg*sg + bg*bg,sqrt( sg*sg + bg*bg )*foil_thickness ",sg*sg + bg*bg,sqrt( sg*sg + bg*bg )*foil_thickness
!         
!             nn = max( 0, nint( sqrt( sg*sg + bg*bg )*foil_thickness ) )
!            
!         !
!         !       to reach this point therefore we want to choose a new sg'
!         !       
!         !           sg' = +/- Sqrt[ (n'/L)^2 - bg^2 ]
!         !
!         !       with sign chosen so sg sg' > 0
!         !
!         !       so what is the best n' to choose?
!         !       I want to minimise | sg' - sg | with abs(n') > L bg 
!         !               
!             
!         !   try n' = n + 1/2
!             if ( abs(nn+0.5) >= foil_thickness*bg ) then
!                 alpha = sign( sqrt( (nn*nn+nn+0.25d0)/(foil_thickness*foil_thickness) - bg*bg ) , sg )
!             else
!                 alpha = huge(1.0)
!             end if
!             
!         !   try n' = n - 1/2
!             if ( abs(nn-0.5) >= foil_thickness*bg ) then
!                 beta = sign( sqrt( (nn*nn-nn+0.25d0)/(foil_thickness*foil_thickness) - bg*bg ) , sg )
!             else
!                 beta = huge(1.0)
!             end if
!         
!             if ( min(alpha,beta)==huge(1.0) ) then
!                 print *,"Lib_DiffractionConditions::tiltToBrightField ERROR - can't find a good n' = n +/- 1/2"
!                 print *,"       bg,sg,L,n = ",bg,sg,foil_thickness,nn 
!                 stop
!             else if ( abs(alpha-sg) < abs(beta-sg) ) then
!                 sgp = alpha
!                 print *,"Lib_DiffractionConditions::tiltToBrightField() info - n (DF) = ",nn," -> n (BF) = ",nn+0.5d0," sg' = ",sgp
!             else 
!                 sgp = beta
!                 print *,"Lib_DiffractionConditions::tiltToBrightField() info - n (DF) = ",nn," -> n (BF) = ",nn-0.5d0," sg' = ",sgp
!             end if
!               
!             
!             
!             
!         !   now find a small rotation about the y-axis which transforms sg to sg'
!         !           (  cos(theta)   0       sin(theta) )
!         !       R = (     0         1          0       )
!         !           ( -sin(theta)   0       cos(theta) )
!         !
!         !   under this rotation 
!         !
!         !       g' = R g
!         !
!         !       sg' = -1/( 2 pi |k| ) ( 2 k.g' + g'.g' )
!         !   
!         !           = -1/( 2 pi ) ( -2 sin(theta) g1 + 2 cos(theta) g3 + |g|^2/|k| )
!         !
!         !   which is solved by 
!         !
!         !       cos(theta) = (a +/- sqrt b)/c
!         !
!         !       with a = g3 |k| ( |g|^2  + 2 |k| pi sg' )
!         !            b = g1^2 |k|^2 ( 4 |k|^2( g1^2 + g3^2) - ( |g|^2 + 2|k| pi sg' )^2 )
!         !            c = - 2 |k|^2( g1^2 + g3^2 )
!         !
!         !       but note g2 = 0, |g|^2 = g1^2 + g3^2, and so
!         !
!             g2      = gg(1)*gg(1) + gg(3)*gg(3)
!             kmod    = ME*v/HBAR
!             
!             alpha = kmod * gg(3) * ( g2 + 2*kmod*PI*sgp )
!             beta  = gg(1)*gg(1) * kmod*kmod * ( 4*kmod*kmod*g2 - ( g2 + 2*kmod * PI * sgp )**2 )
!             gamma = -2 * kmod*kmod * g2             
!             
!             print *,"dbg sgp,g2,kmod,alpha,beta,gamma,HBAR*HBAR/(2*ME),v ",sgp,g2,kmod,alpha,beta,gamma,HBAR*HBAR/(2*ME),v
!             
!             
!             print *,"eps = ", (HBAR*HBAR/(2*ME))*( 2*kmod*gg(3) + g2 )
!             
!             if (beta<0) then
!                 print *,"Lib_DiffractionConditions::tiltToBrightField ERROR - can't find a good rotation angle to take sg->sg'"
!                 print *,"    cos(theta) = (a +/- sqrt b)/c"
!                 print *,"             a = ",alpha
!                 print *,"             b = ",beta
!                 print *,"             c = ",gamma
!                 stop
!             else
!                 beta = sqrt(beta)
!             end if
!             
!             
!         !   want to take the smallest angle, which minimises |g3'|
!         !   compute all four options...    
!                  
!             cost1 = (alpha + beta)/gamma
!             sint1 = sqrt( max(0.0d0,1 - cost1*cost1) )     !   max is to avoid sqrt( -0 ), actually |cos theta|<=1 
!         
!             g31 = abs( -sint1*gg(1) + cost1*gg(3) )
!             g32 = abs( sint1*gg(1) + cost1*gg(3) ) 
!             
!             cost2 = (alpha - beta)/gamma
!             sint2 = sqrt( max(0.0d0,1 - cost2*cost2) )      
!                     
!             g33 = abs( -sint2*gg(1) + cost2*gg(3) )
!             g34 = abs( sint2*gg(1) + cost2*gg(3) ) 
!             
!             
!             
!             if ( g31 <= minval( (/g32,g33,g34/) ) ) then
!                 RR(1,1:3) = (/ cost1 , 0.0d0 , sint1 /)
!                 RR(2,1:3) = (/ 0.0d0 , 1.0d0 , 0.0d0 /)
!                 RR(3,1:3) = (/-sint1 , 0.0d0 , cost1 /)
!             else if ( g32 <= minval( (/g31,g33,g34/) ) ) then
!                 RR(1,1:3) = (/ cost1 , 0.0d0 ,-sint1 /)
!                 RR(2,1:3) = (/ 0.0d0 , 1.0d0 , 0.0d0 /)
!                 RR(3,1:3) = (/ sint1 , 0.0d0 , cost1 /)
!             else if ( g33 <= minval( (/g31,g32,g34/) ) ) then
!                 RR(1,1:3) = (/ cost2 , 0.0d0 , sint2 /)
!                 RR(2,1:3) = (/ 0.0d0 , 1.0d0 , 0.0d0  /)
!                 RR(3,1:3) = (/-sint2 , 0.0d0 , cost2 /)
!             else  
!                 RR(1,1:3) = (/ cost2 , 0.0d0 ,-sint2 /)
!                 RR(2,1:3) = (/ 0.0d0 , 1.0d0 , 0.0d0 /)
!                 RR(3,1:3) = (/ sint2 , 0.0d0 , cost2 /)
!             end if
!             
!             
! !             
!         !---    not sure about this condition 20/01/23 
!         !       can give a -ve sqrt             
!         !    if ( (sqrt( aa*aa + bg*bg )*foil_thickness/PI )<nn .and. nn>0 ) then
!       !        print *,"Lib_DiffractionConditions::tiltToBrightField() info - n (BF)           = ",sqrt( aa*aa + bg*bg )*foil_thickness/PI,"->",nn," -1/2"
!       !        ap = sign( sqrt( ((nn-0.5d0)*PI/foil_thickness)**2 - bg*bg ), aa )             
!         !    else
!       !        print *,"Lib_DiffractionConditions::tiltToBrightField() info - n (BF)           = ",sqrt( aa*aa + bg*bg )*foil_thickness/PI,"->",nn," +1/2"
!       !        ap = sign( sqrt( ((nn+0.5d0)*PI/foil_thickness)**2 - bg*bg ), aa )
!       !    end if           
!           
!              
!              
!              !print *,"Lib_DiffractionConditions::tiltToBrightField() info - a,a'           = ",aa,ap         
!              uu = (bg*bg/(aa*aa + bg*bg)) * Sin( sqrt(aa*aa+bg*bg)*foil_thickness )**2
!              print *,"Lib_DiffractionConditions::tiltToBrightField() info - |phi_g(z)^2|     =  ",uu," in good (n,ng) diffraction conditions"
!  
!              
!       !---    k-vector in rotated frame is exactly parallel to [001], which may not be exactly down the zone axis.    
!           kmod    = (ME*v/HBAR)       
!               kk(1:3) = kmod*(/ 0,0,1 /)
!            !  print *,"|k| ",kmod    
!                
!       !---    find desired deviation parameter            
!           uu = -4*kmod*aa
!           up = -4*kmod*ap
!           k2 = kmod*kmod
!           kg = kmod*gg(3)
!           g2 = gg(1)*gg(1) + gg(2)*gg(2) + gg(3)*gg(3)
!             
!               
!            !  print *,"2k.g + |g|^2 ", 2*kg + g2, uu ," -> ",up
!                       
!            
!           dd = (kg*kg - g2*k2)
!            
!           print *,"dbg??  dd*(g2*g2 + up*up - 2*g2*(up + 2*k2)) = ", dd*(g2*g2 + up*up - 2*g2*(up + 2*k2))
!           
!           
!           alpha = sqrt( dd*(g2*g2 + up*up - 2*g2*(up + 2*k2)) )
!           dd = 1/(2*dd)
!           alpha = abs( alpha * dd )
!           if (dd*kg>0) then                               
!               beta  = ( (up - g2)/2 + abs(kg)*alpha )/g2      !!
!           else
!               beta  = ( (up - g2)/2 - abs(kg)*alpha )/g2
!           end if
!           
!                           
!           !print *," k -> alpha k + beta g ",alpha,beta
!           
!           kk = alpha*kk + beta*gg
!           k2 = dd*dd
!           kg = dot_product( kk,gg )
!           !print *," |k| ",kmod,norm2(kk)
!           
!           uu = (HBAR*HBAR/(2*ME))*(2*kg + g2)
!           zz = kk/norm2(kk)
!           
!           !
!           aa = PI*energyToLengthDeviationParameter( v,uu )
!           
!           print *,"dbg??  aa*aa + bg*bg = ", aa*aa + bg*bg
!           !print *,"half periods           ",sqrt( aa*aa + bg*bg )*foil_thickness/PI 
!           uu = ( (bg*bg)/(aa*aa + bg*bg) ) * Sin( sqrt(aa*aa+bg*bg)*foil_thickness )**2
!             print *,"Lib_DiffractionConditions::tiltToBrightField() info - |phi_g(z)^2|     =  ",uu," tilted to bright field"
!                   
!             xx = gg - dot_product(gg,zz)*zz
!             xx = xx/norm2(xx)
!             yy = cross_product( zz,xx ) 
!             RR(1:3,1) = xx
!             RR(1:3,2) = yy
!             RR(1:3,3) = zz
            

            R = matmul(RR,R)
        
          ! if (rank == 0) then
          !     print *,"Lib_DiffractionConditions::tiltToBrightField info - tilt rotation matrix to bring to bright field conditions"
          !     write (*,fmt='(3f12.5,a,3f12.5)') RR(1,:) 
          !     write (*,fmt='(3f12.5,a,3f12.5)') RR(2,:) 
          !     write (*,fmt='(3f12.5,a,3f12.5)') RR(3,:) 
          !     
          !     print *,"Lib_DiffractionConditions::tiltToBrightField info - total rotation matrix to bring to bright field conditions"
          !     write (*,fmt='(3f12.5,a,3f12.5)') R(1,:) 
          !     write (*,fmt='(3f12.5,a,3f12.5)') R(2,:) 
          !     write (*,fmt='(3f12.5,a,3f12.5)') R(3,:) 
          !      
          !         
!         !      write(*,fmt='(a)',advance="no") "Lib_DiffractionConditions::tiltToBrightField info - foil tilt -D2BI_R "
!         !      write(aaaa,fmt='(f16.12)') R(1,1) ; write(*,fmt='(a)',advance="no")  trim(adjustl(aaaa))//","
!         !      write(aaaa,fmt='(f16.12)') R(2,1) ; write(*,fmt='(a)',advance="no")  trim(adjustl(aaaa))//","
!         !      write(aaaa,fmt='(f16.12)') R(3,1) ; write(*,fmt='(a)',advance="no")  trim(adjustl(aaaa))//","
!         !      write(aaaa,fmt='(f16.12)') R(1,2) ; write(*,fmt='(a)',advance="no")  trim(adjustl(aaaa))//","
!         !      write(aaaa,fmt='(f16.12)') R(2,2) ; write(*,fmt='(a)',advance="no")  trim(adjustl(aaaa))//","
!         !      write(aaaa,fmt='(f16.12)') R(3,2) ; write(*,fmt='(a)',advance="no")  trim(adjustl(aaaa))//","
!         !      write(aaaa,fmt='(f16.12)') R(1,3) ; write(*,fmt='(a)',advance="no")  trim(adjustl(aaaa))//","
!         !      write(aaaa,fmt='(f16.12)') R(2,3) ; write(*,fmt='(a)',advance="no")  trim(adjustl(aaaa))//","
!         !      write(aaaa,fmt='(f16.12)') R(3,3) ; write(*,fmt='(a)',advance="yes") trim(adjustl(aaaa)) 
          ! end if            
            
            RR = matmul( R,U )                                                   
            kk = reflection( BB, hkl )                                          !   using kk as a dummy briefly...                                                  
            gg(1:3) = RR(1:3,1)*kk(1) + RR(1:3,2)*kk(2) + RR(1:3,3)*kk(3)       !   ...because I don't need it otherwise         
             
            
        !---      compute g-vector in rotated frame. This will be close to [100], but not identical to it.
            
            if (rank == 0) then
                print *,"Lib_DiffractionConditions::tiltToBrightField() info - [hkl], g-vec in lab frame, g-vec after rotation"
                write (*,fmt='(i4,a,f12.5,a,f12.5)') hkl(1),"    ",kk(1),"    ",gg(1)
                write (*,fmt='(i4,a,f12.5,a,f12.5)') hkl(2),"    ",kk(2),"    ",gg(2)
                write (*,fmt='(i4,a,f12.5,a,f12.5)') hkl(3),"    ",kk(3),"    ",gg(3)
            end if    
            
            
            
            uu = energy_deviation_parameter( v,(/0.0d0,0.0d0,1.0d0/),gg )
            sg = energyToLengthDeviationParameter( v,uu )
            if (rank == 0) print *,"Lib_DiffractionConditions::tiltToBrightField info - deviation parameter (after)  ",uu," (eV) ",sg," (1/A)"
            if (rank == 0) print *,"Lib_DiffractionConditions::tiltToBrightField info - | phi_g(L) |^2 (after)   ", (bg*sin( PI* sqrt( sg*sg + bg*bg )* foil_thickness ))**2 / ( sg*sg + bg*bg )
            
            
            !stop
            
            return
        end subroutine tiltToBrightField
        
        subroutine tiltToDarkField( v,A_cell, hkl , xi0,xig,foil_thickness,U, R )
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      given the electron velocity v ( in A/fs )
    !*      and the lattice vectors describing the conventional unit cell A_cell
    !*      and a g-vector reflection
    !*      and a lattice type ( so that permitted reflections are found correctly )
    !*      and the tilt-stage orientation R which gives good 2 beam conditions
    !*      and the tilt-stage tweak R which rotates the crystal into a good dark field condition
            
            real(kind=real64),intent(in)                    ::      v
            real(kind=real64),dimension(3,3),intent(in)     ::      A_cell              !   lattice vectors of conventional unit cell
            integer,dimension(3),intent(in)                 ::      hkl
            real(kind=real64),intent(in)                    ::      xi0,xig,foil_thickness
            real(kind=real64),dimension(3,3),intent(in)     ::      U
            real(kind=real64),dimension(3,3),intent(inout)  ::      R
            
            
                        
            
            real(kind=real64)               ::      bg,b0,sg , sgp 
            real(kind=real64)               ::      uu 
            real(kind=real64)               ::      cost,sint


            real(kind=real64),dimension(3,3)    ::      BB,RR,RU
            real(kind=real64),dimension(3)      ::      kk,gg
            
            real(kind=real64)               ::      theta,lp,const,phiL2
 
            
            type(GoldenSection)             ::      gold
            real(kind=real64)               ::      x1,x2,x3,y1,y2,y3
            logical                         ::      isConverged,isWithinTimeLimit
        
!            character(len=16)                       ::      aaaa
            
            
            integer                             ::      nProcs,rank,ierror
            
            call MPI_COMM_SIZE(MPI_COMM_WORLD, nProcs, ierror)
            call MPI_COMM_RANK(MPI_COMM_WORLD, rank, ierror)
            
            
            
        !----   compute g-vector in rotated frame. This will be close to [100], but not identical to it.
            BB = getReciprocalLatticeVectors( A_cell )
            !call inverse3Mat(transpose(A_cell),BB) ; BB = 2*PI*BB                                              
            kk = reflection( BB, hkl )                                                                  
            RU = matmul( R,U )
            gg(1:3) = RU(1:3,1)*kk(1) + RU(1:3,2)*kk(2) + RU(1:3,3)*kk(3) 
            !print *,"gg = ",gg
            
            if (rank==0) then
                write (*,fmt='(a)') "Lib_DiffractionConditions::tiltToDarkField() info - [hkl], g-vec in lab frame, g-vec after rotation"
                write (*,fmt='(i4,a,f12.5,a,f12.5)') hkl(1),"    ",kk(1),"    ",gg(1)
                write (*,fmt='(i4,a,f12.5,a,f12.5)') hkl(2),"    ",kk(2),"    ",gg(2)
                write (*,fmt='(i4,a,f12.5,a,f12.5)') hkl(3),"    ",kk(3),"    ",gg(3)
            end if
                            
        !---    propagation of electrons in perfect xtal
        !               
        !       d/dz (  phi_0 ) = PI i  ( b0    bg        )      !   with b0 = 1/xi0, bg = 1/xig
        !            (  phi_g )         ( bg    b0 + 2 sg )      !         
        !            
            bg = 1/xig
            b0 = 1/xi0
            uu = energy_deviation_parameter( v,(/0.0d0,0.0d0,1.0d0/),gg )
            sg = energyToLengthDeviationParameter(v,uu)
            phiL2 = (bg*sin( PI* sqrt( sg*sg + bg*bg )* foil_thickness ))**2 / ( sg*sg + bg*bg )    
            if (rank==0) then      
                write (*,fmt='(3(a,f16.6))') "Lib_DiffractionConditions::tiltToDarkField info - deviation parameter (before) ",uu," (eV) ",sg," (1/A)"
                write (*,fmt='(3(a,f16.6))') "Lib_DiffractionConditions::tiltToDarkField info - foil thickness L  = ",foil_thickness  
                write (*,fmt='(3(a,f16.6))') "Lib_DiffractionConditions::tiltToDarkField info - | phi_g(L) |^2 (before)  ", (bg*sin( PI* sqrt( sg*sg + bg*bg )* foil_thickness ))**2 / ( sg*sg + bg*bg )
                !write (*,fmt='(3(a,f16.6))') "Lib_DiffractionConditions::tiltToDarkField info - b0,bg,sg                ",b0,bg,sg
            end if            
            
            const = - HBAR*(gg(1)*gg(1)+gg(2)*gg(2)+gg(3)*gg(3))/(4*PI*ME*v)
            
            x1 = -THETA_MAX*0.01d0*DEGTORAD
            cost = cos(x1) ; sint = sin(x1)            
            Lp  = foil_thickness/cost
            sgp = ( sint*gg(1) - cost*gg(3) )/(2*PI)            
            y1 = (bg*sin( PI* sqrt( sgp*sgp + bg*bg )*Lp ))**2 / ( sgp*sgp + bg*bg )
            
            x3 = +THETA_MAX*0.01d0*DEGTORAD
            cost = cos(x3) ; sint = sin(x3)            
            Lp  = foil_thickness/cost
            sgp = ( sint*gg(1) - cost*gg(3) )/(2*PI)            
            y3 = (bg*sin( PI* sqrt( sgp*sgp + bg*bg )*Lp ))**2 / ( sgp*sgp + bg*bg )
            
                                   
            gold = GoldenSection_ctor( x1,x3, y1,y3)
            do
                call nextPoint(gold,x2)
                cost = cos(x2) ; sint = sin(x2)            
                Lp  = foil_thickness/cost
                sgp = ( sint*gg(1) - cost*gg(3) )/(2*PI)            
                y2 = (bg*sin( PI* sqrt( sgp*sgp + bg*bg )*Lp ))**2 / ( sgp*sgp + bg*bg )
            
                call minimise(gold,x2,y2,isConverged,isWithinTimeLimit)
                if (isConverged) exit
                if (.not. isWithinTimeLimit) exit
            end do
            
      
            if (y2 < phiL2) then
                theta = x2 ; cost = cos(x2) ; sint = sin(x2)  
            else
                theta = 0.0d0; cost = 1.0d0 ; sint = 0.0d0
            end if
            sgp = ( sint*gg(1) - cost*gg(3) )/(2*PI)      
            if (rank==0) then                
                write (*,fmt='(3(a,f16.6))') "Lib_DiffractionConditions::tiltToDarkField info - best rotation angle     ",theta," (rad)"
                !print *,"Lib_DiffractionConditions::tiltToDarkField info - b0,bg,sg'               ",b0,bg,sgp
                !print *,"Lib_DiffractionConditions::tiltToDarkField info - intensity variation     ",bg*bg/( sgp*sgp+bg*bg )
                write (*,fmt='(3(a,f16.6))') "Lib_DiffractionConditions::tiltToDarkField info - | phi_g(L) |^2 (after)  ",y2
            end if            
            
            RR(1,1:3) = (/ cost , 0.0d0 , sint /)
            RR(2,1:3) = (/ 0.0d0 , 1.0d0 , 0.0d0 /)
            RR(3,1:3) = (/-sint , 0.0d0 , cost /)
            
            
            R = matmul(RR,R)
            
            
            if (rank==0) then   
                write(*,fmt='(a)') "Lib_DiffractionConditions::tiltToDarkField info - tilt rotation matrix to bring to Dark field conditions"
                write (*,fmt='(3f12.5,a,3f12.5)') RR(1,:) 
                write (*,fmt='(3f12.5,a,3f12.5)') RR(2,:) 
                write (*,fmt='(3f12.5,a,3f12.5)') RR(3,:) 
               
               ! write(*,fmt='(a)') "Lib_DiffractionConditions::tiltToDarkField info - total rotation matrix used"
               ! write (*,fmt='(3f12.5,a,3f12.5)') R(1,:) 
               ! write (*,fmt='(3f12.5,a,3f12.5)') R(2,:) 
               ! write (*,fmt='(3f12.5,a,3f12.5)') R(3,:) 
                 
                    
                ! write(*,fmt='(a)',advance="no") "Lib_DiffractionConditions::tiltToDarkField info - foil rotation -R "
                ! write(aaaa,fmt='(f16.12)') R(1,1) ; write(*,fmt='(a)',advance="no")  trim(adjustl(aaaa))//","
                ! write(aaaa,fmt='(f16.12)') R(2,1) ; write(*,fmt='(a)',advance="no")  trim(adjustl(aaaa))//","
                ! write(aaaa,fmt='(f16.12)') R(3,1) ; write(*,fmt='(a)',advance="no")  trim(adjustl(aaaa))//","
                ! write(aaaa,fmt='(f16.12)') R(1,2) ; write(*,fmt='(a)',advance="no")  trim(adjustl(aaaa))//","
                ! write(aaaa,fmt='(f16.12)') R(2,2) ; write(*,fmt='(a)',advance="no")  trim(adjustl(aaaa))//","
                ! write(aaaa,fmt='(f16.12)') R(3,2) ; write(*,fmt='(a)',advance="no")  trim(adjustl(aaaa))//","
                ! write(aaaa,fmt='(f16.12)') R(1,3) ; write(*,fmt='(a)',advance="no")  trim(adjustl(aaaa))//","
                ! write(aaaa,fmt='(f16.12)') R(2,3) ; write(*,fmt='(a)',advance="no")  trim(adjustl(aaaa))//","
                ! write(aaaa,fmt='(f16.12)') R(3,3) ; write(*,fmt='(a)',advance="yes") trim(adjustl(aaaa)) 
            end if            
            
            RU = matmul( R,U )                                                   
            kk = reflection( BB, hkl )                                          !   using kk as a dummy briefly...                                                  
            gg(1:3) = RU(1:3,1)*kk(1) + RU(1:3,2)*kk(2) + RU(1:3,3)*kk(3)       !   ...because I don't need it otherwise         
             
            
        !---      compute g-vector in rotated frame. This will be close to [100], but not identical to it.
            
            if (rank==0) then   
                write(*,fmt='(a)') "Lib_DiffractionConditions::tiltToDarkField() info - [hkl], g-vec in lab frame, g-vec after rotation"
                write (*,fmt='(i4,a,f12.5,a,f12.5)') hkl(1),"    ",kk(1),"    ",gg(1)
                write (*,fmt='(i4,a,f12.5,a,f12.5)') hkl(2),"    ",kk(2),"    ",gg(2)
                write (*,fmt='(i4,a,f12.5,a,f12.5)') hkl(3),"    ",kk(3),"    ",gg(3)
            end if   
            
            
            
            uu = energy_deviation_parameter( v,(/0.0d0,0.0d0,1.0d0/),gg )
            sg = energyToLengthDeviationParameter( v,uu )
            if (rank==0) then
                write (*,fmt='(3(a,f16.6))') "Lib_DiffractionConditions::tiltToDarkField info - deviation parameter (after)  ",uu," (eV) ",sg," (1/A)"
                !print *,"Lib_DiffractionConditions::tiltToDarkField info - | phi_g(L) |^2 (after)   ", (bg*sin( PI* sqrt( sg*sg + bg*bg )* foil_thickness ))**2 / ( sg*sg + bg*bg )
            end if            
            
            !stop
            
            return
        end subroutine tiltToDarkField
        
    
        subroutine findCrystalliteSideLengths( A_super,R, d )
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      given the original periodic supercell A_super
    !*      and the suggestion that we rotate it by R
    !*      to get the best diffraction condition
    !*      find a suspended cuboidal crystallite 
    !*      whose sides are length d1,d2,d3
    !*      ( described here by the matrix D )
    !*      and a translation delta
    !*      and a set of periodic replicas required of the original supercell to fill D
    !*      given the transformation
    !*          x' = R (x + A[uvw]) + delta
    
            real(kind=real64),dimension(3,3),intent(in)     ::          A_super
            real(kind=real64),dimension(3,3),intent(in)     ::          R
            real(kind=real64),dimension(3),intent(inout)    ::          d
            
            real(kind=real64),dimension(3)          ::      xx 
            real(kind=real64),dimension(3,50)       ::      corner = reshape(           &
                              (/  0,0,0   ,   4,0,0   ,   0,4,0   ,   4,4,0   ,         &       !   corners
                                  0,0,4   ,   4,0,4   ,   0,4,4   ,   4,4,4   ,         &
                                  0,0,2 ,   0,2,0   ,   0,0,2   ,   4,0,2   ,           &       !   edges
                                  4,2,0 ,   0,4,2   ,   2,4,0   ,   0,2,4   ,           &
                                  2,0,4 ,   4,4,2   ,   4,2,4   ,   2,4,4   ,           &  
                                  2,2,0 ,   2,2,4   ,   2,0,2   ,   2,4,0   ,   0,2,2   ,   4,2,2   ,   &       !   face centres
                                  0,1,1 ,   0,1,3   ,   0,3,1   ,   0,3,3   ,           &       !   face quarter points
                                  4,1,1 ,   4,1,3   ,   4,3,1   ,   4,3,3   ,           &        
                                  1,0,1 ,   1,0,3   ,   3,0,1   ,   3,0,3   ,           &        
                                  1,4,1 ,   1,4,3   ,   3,4,1   ,   3,4,3   ,           &        
                                  1,1,0 ,   1,3,0   ,   3,1,0   ,   3,3,0   ,           &        
                                  1,1,4 ,   1,3,4   ,   3,1,4   ,   3,3,4               &                                   
                                   /),(/3,50/) ) * 0.25d0
!            real(kind=real64),dimension(3,8)       ::      corner = reshape(           &
!                            (/  0,0,0   ,   2,0,0   ,   0,2,0   ,   2,2,0   ,           &
!                                0,0,2   ,   2,0,2   ,   0,2,2   ,   2,2,2      /),(/3,8/) ) * 0.5d0
!                     
            real(kind=real64)                   ::      xj,tmin,tmax
            integer     ::      ii,jj         
                                
        !---    step 1: find acceptable side lengths for the suspended crystallite
                      
            
        
            do jj = 1,3
                if ( (d(jj) == 0).or.(d(jj) == SG_UNSET) ) then                             
                    print *,"Lib_DiffractionConditions::findCrystalliteSideLengths info - setting dimension ",jj," from supercell corners"  
                    tmin = huge(1.0)
                    tmax = -huge(1.0)                    
                    do ii = 1,8
                        xx(1:3) = corner(1:3,ii) - 0.5d0
                        xx(1:3) = A_super(1:3,1)*xx(1) + A_super(1:3,2)*xx(2) + A_super(1:3,3)*xx(3)
                        xj = R(jj,1)*xx(1) + R(jj,2)*xx(2) + R(jj,3)*xx(3)
                        tmin = min( tmin,xj )
                        tmax = max( tmax,xj )
                    end do
                    d(jj) = tmax - tmin                         
                end if
            end do
            
            !print *,"Lib_DiffractionConditions::findCrystalliteSideLengths info - crystallite side lengths ",dd
            
            return
        end subroutine findCrystalliteSideLengths
        
         
            
        subroutine findTransformedSupercell( A_super,R,D_super,delta,uvw,nuvw , foilthickness )
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      given the original periodic supercell A_super
    !*      and the suggestion that we rotate it by R
    !*      to get the best diffraction condition
    !*      find a suspended cuboidal crystallite 
    !*      whose sides are length d1,d2,d3
    !*      ( described here by the matrix D )
    !*      and a translation delta
    !*      and a set of periodic replicas required of the original supercell to fill D
    !*      given the transformation
    !*          x' = R (x + A[uvw]) + delta
    
            real(kind=real64),dimension(3,3),intent(in)     ::          A_super
            real(kind=real64),dimension(3,3),intent(in)     ::          R
            real(kind=real64),dimension(3,3),intent(out)    ::          D_super
            real(kind=real64),dimension(3),intent(out)      ::          delta
            integer,dimension(3,125),intent(out)            ::          uvw
            integer,intent(out)                             ::          nuvw
            real(kind=real64),dimension(3),intent(in)       ::          foilthickness
            
            real(kind=real64),dimension(3)          ::      dd,xx 
!              real(kind=real64),dimension(3,50)     ::      corner = reshape(           &
!                              (/  0,0,0   ,   4,0,0   ,   0,4,0   ,   4,4,0   ,         &       !   corners
!                                  0,0,4   ,   4,0,4   ,   0,4,4   ,   4,4,4   ,         &
!                                  0,0,2 ,   0,2,0   ,   0,0,2   ,   4,0,2   ,           &       !   edges
!                                  4,2,0 ,   0,4,2   ,   2,4,0   ,   0,2,4   ,           &
!                                  2,0,4 ,   4,4,2   ,   4,2,4   ,   2,4,4   ,           &  
!                                  2,2,0 ,   2,2,4   ,   2,0,2   ,   2,4,0   ,   0,2,2   ,   4,2,2   ,   &       !   face centres
!                                  0,1,1 ,   0,1,3   ,   0,3,1   ,   0,3,3   ,           &       !   face quarter points
!                                  4,1,1 ,   4,1,3   ,   4,3,1   ,   4,3,3   ,           &        
!                                  1,0,1 ,   1,0,3   ,   3,0,1   ,   3,0,3   ,           &        
!                                  1,4,1 ,   1,4,3   ,   3,4,1   ,   3,4,3   ,           &        
!                                  1,1,0 ,   1,3,0   ,   3,1,0   ,   3,3,0   ,           &        
!                                  1,1,4 ,   1,3,4   ,   3,1,4   ,   3,3,4               &                                   
!                                   /),(/3,50/) ) * 0.25d0
!            real(kind=real64),dimension(3,8)       ::      corner = reshape(           &
!                            (/  0,0,0   ,   2,0,0   ,   0,2,0   ,   2,2,0   ,           &
!                                0,0,2   ,   2,0,2   ,   0,2,2   ,   2,2,2      /),(/3,8/) ) * 0.5d0
!                     
           
            integer                             ::      ix,iy,iz
           
                                 
                                
        !---    step 1: find acceptable side lengths for the suspended crystallite
                      
            dd = foilthickness
            !call findCrystalliteSideLengths( A_super,R, dd )
        
            !do jj = 1,3
            !   if ( (foilthickness(jj) == 0).or.(foilthickness(jj) == SG_UNSET) ) then                             
            !        print *,"Lib_DiffractionConditions::findTransformedSupercell info - setting dimension ",jj," from supercell corners"   
            !        tmin = huge(1.0)
            !        tmax = -huge(1.0)                    
            !        do ii = 1,8
            !            xx(1:3) = corner(1:3,ii) - 0.5d0
            !            xx(1:3) = A_super(1:3,1)*xx(1) + A_super(1:3,2)*xx(2) + A_super(1:3,3)*xx(3)
            !            xj = R(jj,1)*xx(1) + R(jj,2)*xx(2) + R(jj,3)*xx(3)
            !            tmin = min( tmin,xj )
            !            tmax = max( tmax,xj )
            !        end do
            !       dd(jj) = tmax - tmin            
            !   else
            !       !print *,"Lib_DiffractionConditions::findTransformedSupercell info - setting dimension ",jj," to ",foilthickness(jj)                            
            !       dd(jj) = foilthickness(jj) 
            !   end if
            !end do
            !
            !print *,"Lib_DiffractionConditions::findTransformedSupercell info - crystallite side lengths ",dd
!             write (*,fmt='(3f12.4)') tmin(1),tmax(1),dd(1)
!             write (*,fmt='(3f12.4)') tmin(2),tmax(2),dd(2)
!             write (*,fmt='(3f12.4)') tmin(3),tmax(3),dd(3)
            
            D_super = 0.0d0
            D_super(1,1) = dd(1)
            D_super(2,2) = dd(2)
            D_super(3,3) = dd(3)
            
        !   print *,"Lib_DiffractionConditions::findTransformedSupercell info - supercell, rotation, box side lengths"
        !   write (*,fmt='(3f16.8,a,3f16.8,a,f16.8)') A_super(1,:),"    ",R(1,:),"    ",dd(1)                                               
        !   write (*,fmt='(3f16.8,a,3f16.8,a,f16.8)') A_super(2,:),"    ",R(2,:),"    ",dd(2)                                               
        !   write (*,fmt='(3f16.8,a,3f16.8,a,f16.8)') A_super(3,:),"    ",R(3,:),"    ",dd(3)                                               
        !---    step 2: determine the translation vector delta
            xx(1:3) = 0.5d0
            xx(1:3) = A_super(1:3,1)*xx(1) + A_super(1:3,2)*xx(2) + A_super(1:3,3)*xx(3)
            xx(1:3) = R(1:3,1)*xx(1) + R(1:3,2)*xx(2) + R(1:3,3)*xx(3)
            delta(1:3) = dd(1:3)/2 - xx(1:3)
        !    write(*,fmt='(a,3f12.4)') " Lib_DiffractionConditions::findTransformedSupercell info - delta ",delta(1:3)
            
        !---    step 3: find which periodic replicas are required
        
        
        !   21/09/22 - I can't get this to work reliably every time. 
        !   surely the algorithm below is correct?? But sometimes it misses ??
        
!-----------   GARBAGE SLOW ALGORITHM
            nuvw = 0
            do iz = -1,1
                do iy = -1,1
                    do ix = -1,1
                        nuvw = nuvw + 1
                        uvw(1:3,nuvw) = (/ix,iy,iz/)
                    end do
                end do
            end do
            return
!-----------           
        
        
        
        
        
        
        
        
!        
!            call inverse3Mat(A_super,iA_super)
!            nuvw = 1
!            uvw(1:3,nuvw) = 0
!            !write(*,fmt='(a,i4,a,3i4)') " Lib_DiffractionConditions::findTransformedSupercell info - [uvw] replica ",1," required ",0,0,0
!            !   solving 
!            !           R(x + u) + delta = D y              where D are the corners of the desired box, and u are periodic replica vectors
!           
!            
!            do iz = -2,2
!                do iy = -2,2
!                    do ix = -2,2
!                        do ii = 1,size(corner,dim=2)
!                            yy(1:3) = corner(1:3,ii) + (/ix,iy,iz/)                                             !   position of corner in periodic replica in cell space
!                            xx(1:3) = A_super(1:3,1)*yy(1) + A_super(1:3,2)*yy(2) + A_super(1:3,3)*yy(3)        !   position in real space
!                            yy(1:3) = R(1:3,1)*xx(1) + R(1:3,2)*xx(2) + R(1:3,3)*xx(3)                          !   position rotated
!                            xx(1:3) = yy(1:3) + delta(1:3)  
!                            yy(1) = xx(1) / dd(1) ; yy(2) = xx(2) / dd(2) ; yy(3) = xx(3) / dd(3)               !   position in image space
!                            uvw_trial(1:3) = floor( yy(1:3) )
!                            if (.not. all(uvw_trial == 0)) cycle
!                            uvw_trial(1:3) = (/ ix,iy,iz /)
!                            
!                            ok = .false.        !   do we have this [uvw] replica already stored?
!                            do kk = 1,nuvw                                                          
!                                ok = ok .or. ( (uvw(1,kk) == uvw_trial(1)).and.(uvw(2,kk) == uvw_trial(2)).and.(uvw(3,kk) == uvw_trial(3)) )
!                            end do
!                            if (.not. ok) then  
!                                nuvw = nuvw + 1
!                                uvw(1:3,nuvw) = uvw_trial(1:3)
!                            end if
!                        end do
!                    end do
!                end do
!            end do                  
!            write(*,fmt='(a,i4,a)') " Lib_DiffractionConditions::findTransformedSupercell info - [uvw] replicas required ",nuvw
!            return
!            
!            
!            
!            
!            
!            
!            
!            
!            
!            
!            
!            
!            do ii = 1,size(corner,dim=2)
!
!            !---    find position of the corner in the coordinate basis of the periodic supercell held  
!                !xx(1:3) = A_super(1:3,1)*corner(1,ii) + A_super(1:3,2)*corner(2,ii) + A_super(1:3,3)*corner(3,ii)           
!                xx(1:3) = dd(1:3)*corner(1:3,ii) - delta(1:3)
!                xx(1:3) = R(1,1:3)*xx(1) + R(2,1:3)*xx(2) + R(3,1:3)*xx(3) 
!                xx(1:3) = iA_super(1:3,1)*xx(1) + iA_super(1:3,2)*xx(2) + iA_super(1:3,3)*xx(3)
!            
!                do jj = 1,8
!                    
!                    yy(1:3) = xx(1:3) + ( 2.0d0 * corner(1:3,jj) - 1.0d0 )*0.001d0      !   add +/- 0.1% to the corner
!                !---    which periodic repeat is this corner in?        
!                    !write(*,fmt='(a,3f12.4)') " Lib_DiffractionConditions::findTransformedSupercell info - test corner ",xx
!                    uvw_trial(1:3) = floor( yy(1:3) )
!                    
!                    ok = .false.        !   do we have this [uvw] replica already stored?
!                    do kk = 1,nuvw                                                          
!                        ok = ok .or. ( (uvw(1,kk) == uvw_trial(1)).and.(uvw(2,kk) == uvw_trial(2)).and.(uvw(3,kk) == uvw_trial(3)) )
!                    end do
!                    if (.not. ok) then  
!                        nuvw = nuvw + 1
!                        uvw(1:3,nuvw) = uvw_trial(1:3)
!                        !write(*,fmt='(a,i4,a,3i4)') " Lib_DiffractionConditions::findTransformedSupercell info - [uvw] replica ",nuvw," required ",uvw_trial(1:3)
!                    end if
!                    
!                end do
!            end do
!            write(*,fmt='(a,i4,a)') " Lib_DiffractionConditions::findTransformedSupercell info - [uvw] replicas required ",nuvw
            return
        end subroutine findTransformedSupercell
                   
        
        
        
        pure function energyToLengthDeviationParameter(v,eps) result(sg)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      convert electron velocity (A/fs) and energy (eV)
    !*      to deviation parameter (1/A)
            real(kind=real64),intent(in)            ::      v,eps
            real(kind=real64)                       ::      sg
            
            sg = -eps/(2*PI*HBAR*v)
            
            return
        end function energyToLengthDeviationParameter
           
        
            
            
        pure function intensityAtDepthZ( xi0,xig,sg, z ) result( phig2 )
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      according to the Howie Whelan equations for a perfect xtal
    !*  
    !*      d/dz ( phi_0 ) = i ( b_0   b_g     ) ( phi_0 )
    !*           ( phi_g )     ( b_g   b_0 + a ) ( phi_g )
    !*
    !*      We can find the value of phi_g(z) analytically by solving the PDE
    !*      with boundary conditions phi_0(0) = 1, phi_g(0) = 0
    !*
    !*          phi_g(z) = ( b_g/sqrt( a^2 + 4 b_g^2 ) ) x
    !*                                  ( Exp[ i lambda+ z ] - Exp[ i lambda- z ] )
    !*
    !*      where lambda+/- are the eigenvalues of the matrix above
    !*
    !*          lambda+/- = (b_0 + a/2) +/- Sqrt[ b_g^2 + (a/2)^2 ]
    !*
    !*      The coefficients are bg = (PI/xig) 
    !*                           b0 = (PI/xi0) 
    !*                           a  = (2*PI*sg)
    !*
            real(kind=real64),intent(in)            ::      xi0,xig             !   extinction distances ( A )
            real(kind=real64),intent(in)            ::      sg                  !   deviation parameter ( 1/A )
            real(kind=real64),intent(in)            ::      z                   !   depth ( crystal thickness ) (A)
            real(kind=real64)                       ::      phig2
            
            complex(kind=real64)        ::      phig            
            real(kind=real64)           ::      bg,b0,aa,lambda1,lambda2,dd
            
            b0 = PI/xi0
            bg = PI/xig
            aa = 2*PI*sg
            
        !---    discriminant 
            dd = sqrt( aa*aa/4 + bg*bg )
            
        !       ... and hence eigenvalues
            lambda1 = (b0 + aa/2) + dd
            lambda2 = (b0 + aa/2) - dd
            
        !---    compute the complex-valued diffracted beam 
            phig = exp( cmplx( 0.0d0,lambda1*z,kind=real64 ) ) - exp( cmplx( 0.0d0,lambda2*z,kind=real64 ) ) 
            phig = phig*bg/(2*dd)
            
       
        !       ... and find intensity 
            phig2 = real( phig*conjg(phig),kind=real64 )        !   note the cast to real is unnecessary mathematically.
            
            return
        end function intensityAtDepthZ
         
!         
!                         
!             
        

        pure function energy_deviation_parameter( v,khat,g ) result( eps )
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      given the electron velocity v ( in A/fs )
    !*      and the incident beam vector direction khat ( unit vector )
    !*      and the g-vector g ( in A^-1 )
    !*      find the (energy) deviation from the Bragg condition eps
    !*
    !*          eps = (HBAR^2/2ME) |k+g|^2 - (HBAR^2/2ME) |k|^2  
    !*          
            real(kind=real64),intent(in)                    ::      v
            real(kind=real64),dimension(3),intent(in)       ::      khat,g
            real(kind=real64)                               ::      eps
            real(kind=real64),parameter         ::      HBAR2ON2M = HBAR*HBAR/(2*ME)        !   note: constants defined in Lib_RelativisticElectrons
          
            real(kind=real64)                   ::      const
            
            const = 2*ME*v/HBAR           !   ME v / HBAR scales the incident beam direction to wavevector in 1/A. constant 2 is used in definition of eps.
            
 
            eps = khat(1)*g(1) + khat(2)*g(2) + khat(3)*g(3)  
                            
            eps = const*eps 
            eps = eps + g(1)*g(1) + g(2)*g(2) + g(3)*g(3)  
            
            eps = eps * HBAR2ON2M
            
            
!           eps = (const*khat(1)+g(1))*g(1)     &
!               + (const*khat(2)+g(2))*g(2)     &
!               + (const*khat(3)+g(3))*g(3)  
!                   
!           eps = eps * HBAR2ON2M
             
            return
        end function energy_deviation_parameter
           
    !---
    
!         subroutine findRotationsForBraggCondition( v,g,cosq,sinq , phi,nphi )
!     !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
!     !*      given a tilt angle 0 < theta < pi/2
!     !*      and a reflection g              ( 1/A )
!     !*      and the electron velocity       ( A/fs )
!     !*      return (four of the eight) possible azimuthal rotations phi which will bring the reflection into Bragg condition (eps=0)
!     !*      note that the solutions "-phi" are also valid.
!     !*
!     !*      This is done by finding the condition where khat.ghkl + (HBAR/2Mv) |ghkl|^2 = 0
!     !*      The unit vector khat is given by
! !********** OLD              khat = ( sin theta cos phi , sin theta sin phi , cos theta )
! 
!     !*      first rotate by theta about x axis then by phi about y axis
!     !*          khat = ( cos theta sin phi, - sin theta, cos theta cos phi )
! 
!     !*      Here theta is fixed, find the values phi
!     !*      
!     !*      note that if theta = 0 or ghkl(1) = ghkl(2) = 0, there is no solution here.
!     !*
!             real(kind=real64),intent(in)                ::      v
!             real(kind=real64),dimension(3),intent(in)   ::      g
!             real(kind=real64),intent(in)                ::      cosq,sinq       !   cos(theta) sin(theta)
!             real(kind=real64),dimension(4),intent(out)  ::      phi
!             integer,intent(out)                         ::      nphi
!             
!             real(kind=real64)           ::      aa , uu, vv , pp
!             real(kind=real64)           ::      xx  !, secq, cscq 
!                          
!             nphi = 0 
!             
!             xx = g(1)*g(1) + g(3)*g(3)             
!             
!             aa = xx + g(2)*g(2)
!             aa = -aa*HBAR/(2*ME*v)          !   constants defined in Lib_RelativisticElectrons
!                                             !   aa = - hbar |g|^2 /( 2 m v )
!             
!             
!                         
!             uu = aa + g(2)*sinq
!             
!             vv = xx*cosq*cosq - uu*uu       !   (g1^2 + g3^2)cos^2 theta - (aa + g2*sin theta)^2
!             
!             
!             if (vv<0) return
!                                             !   sec(theta)/(g1^2 + g3^2)
!             xx = 1/(cosq*xx)
!             
!             vv = sqrt(vv)*g(1)*xx
!             
!             uu = uu*g(3)*xx
!             
!             
!             
!           !  if (max(abs(uu+vv),abs(uu-vv))>1) return
!             
!             
!             if ( abs(uu+vv)<=1 ) then
!                 pp = acos( uu + vv )            
!                 if (abs(pp) <= THETA_MAX) then
!                     phi(1) = pp
!                     phi(2) = -pp
!                     nphi = 2
!                 end if
!             end if 
!             
!             if ( abs(uu-vv)<=1 ) then
!                 pp = acos( uu - vv )
!                 if (abs(pp) <= THETA_MAX) then
!                     phi(nphi+1) = pp
!                     phi(nphi+2) = -pp
!                     nphi = nphi+2
!                 end if
!             end if    
!              
!             
!             
!             return
!         end subroutine findRotationsForBraggCondition
!         
        subroutine findRotationsForBraggCondition2( v,g,cosq,sinq,phi0 , phi,nphi ,nolimitPhi)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      given a tilt angle 0 < theta < pi/2
    !*      and a reflection g              ( 1/A )
    !*      and the electron velocity       ( A/fs )
    !*      return (four of the eight) possible azimuthal rotations phi which will bring the reflection into Bragg condition (eps=0)
    !*      note that the solutions "-phi" are also valid.
    !*
    !*      This is done by finding the condition where khat.ghkl + (HBAR/2Mv) |ghkl|^2 = 0
    !*      The unit vector khat is given by
!********** OLD              khat = ( sin theta cos phi , sin theta sin phi , cos theta )

    !*      first rotate by theta about x axis then by phi about y axis
    !*          khat = ( cos theta sin phi, - sin theta, cos theta cos phi )

    !*      Here theta is fixed, find the values phi
    !*      
    !*      note that if theta = 0 or ghkl(1) = ghkl(2) = 0, there is no solution here.
    !*
            real(kind=real64),intent(in)                ::      v
            real(kind=real64),dimension(3),intent(in)   ::      g
            real(kind=real64),intent(in)                ::      cosq,sinq       !   cos(theta) sin(theta)
            real(kind=real64),intent(in)                ::      phi0
            real(kind=real64),dimension(4),intent(out)  ::      phi
            integer,intent(out)                         ::      nphi
            logical,intent(in),optional                 ::      nolimitPhi
            
            real(kind=real64)           ::      aa , uu, vv , pp
            real(kind=real64)           ::      xx  !, secq, cscq 
            
            integer     ::      jj
             
            nphi = 0 
            
            xx = g(2)*sinq + g(3)*cosq  
            
            aa = g(1)*g(1) + g(2)*g(2) + g(3)*g(3) 
            aa = -aa*HBAR/(2*ME*v)          !   constants defined in Lib_RelativisticElectrons
                                            !   aa = - hbar |g|^2 /( 2 m v )
            
            
                        
            uu = g(1)*g(1) + xx*xx          
            
            vv = uu - aa*aa
            
            
            if (vv<0) then  
               ! print *,"findRotationsForBraggCondition2 vv = ",vv
                return
            end if
                                            !   sec(theta)/(g1^2 + g3^2)
            uu = 1/uu
            
            vv = sqrt(vv)*g(1)*uu
            
            uu = aa*xx*uu
            
            
            
          !  if (max(abs(uu+vv),abs(uu-vv))>1) return
            
            if (present(nolimitPhi)) then
                if (nolimitPhi) then
                    if ( abs(uu+vv)<=1 ) then
                        pp = acos( uu + vv )            
                        jj = nint( (pp - phi0)/(2*PI) )
                        pp = pp - 2*jj*PI
                        phi(1) = pp
                        phi(2) = -pp
                        nphi = 2   
                    end if 
                    
                    if ( abs(uu-vv)<=1 ) then
                        pp = acos( uu - vv )
                        jj = nint( (pp - phi0)/(2*PI) )
                        pp = pp - 2*jj*PI
                        phi(nphi+1) = pp
                        phi(nphi+2) = -pp
                        nphi = nphi+2
                    end if    
                else
                        
                    if ( abs(uu+vv)<=1 ) then
                        pp = acos( uu + vv )            
                        jj = nint( (pp - phi0)/(2*PI) )
                        pp = pp - 2*jj*PI
                        
                        if (abs(pp-phi0) <= 2*THETA_MAX*DEGTORAD) then
                            phi(1) = pp
                            phi(2) = -pp
                            nphi = 2   
                        end if
                    end if 
                    
                    if ( abs(uu-vv)<=1 ) then
                        pp = acos( uu - vv )
                        jj = nint( (pp - phi0)/(2*PI) )
                        pp = pp - 2*jj*PI
                        if (abs(pp-phi0) <= 2*THETA_MAX*DEGTORAD) then
                            phi(nphi+1) = pp
                            phi(nphi+2) = -pp
                            nphi = nphi+2
                        end if
                    end if    
                        
                end if
                return
            end if
            
            if ( abs(uu+vv)<=1 ) then
                pp = acos( uu + vv )            
                jj = nint( (pp - phi0)/(2*PI) )
                pp = pp - 2*jj*PI
                
                if (abs(pp-phi0) <= THETA_MAX*DEGTORAD) then
                    phi(1) = pp
                    phi(2) = -pp
                    nphi = 2   
                end if
            end if 
            
            if ( abs(uu-vv)<=1 ) then
                pp = acos( uu - vv )
                jj = nint( (pp - phi0)/(2*PI) )
                pp = pp - 2*jj*PI
                if (abs(pp-phi0) <= THETA_MAX*DEGTORAD) then
                    phi(nphi+1) = pp
                    phi(nphi+2) = -pp
                    nphi = nphi+2
                end if
            end if    
             
            
            
            return
        end subroutine findRotationsForBraggCondition2
        
!         subroutine findAllIncidentVectors( v,g,theta_max,khat,nkhat )
!     !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
!     !*      first rotate by theta about x axis then by phi about y axis
!     !*          khat = ( cos theta sin phi, - sin theta, cos theta cos phi )
! 
!     
!             real(kind=real64),intent(in)                    ::      v
!             real(kind=real64),dimension(3),intent(in)       ::      g   
!             real(kind=real64),intent(in)                    ::      theta_max         
!             real(kind=real64),dimension(:,:),intent(out)    ::      khat
!             integer,intent(out)                             ::      nkhat
!             
!             integer                         ::      ntheta
!             real(kind=real64)               ::      theta,pp,cost,sint,cosp,sinp 
!             real(kind=real64),dimension(4)  ::      phi
!             integer                         ::      nphi
!             integer         ::      ii,jj
!              
!             
!             ntheta = size(khat,dim=2)/4
!             nkhat = 0
!             do ii = 1,ntheta
!                 theta = theta_max*ii/ntheta
!                 cost = cos(theta) ; sint = sin(theta)
!                 call findRotationsForBraggCondition( v,g,cost,sint,phi,nphi )
!                 do jj = 1,nphi
!                     pp = phi(jj) ; cosp = cos(pp) ; sinp = sin(pp)
!                     khat(:,nkhat+jj) = (/ cost*sinp, -sint, cost*cosp /)
!                 end do
!                 nkhat = nkhat + nphi
!                         
!             end do    
!                 
!                 
!             
!           ! print *,"e dev parm all khat"
!           ! do ii = 1,nkhat
!           !     print *,ii,khat(:,ii),energy_deviation_parameter( v,khat(:,ii),g )
!           ! end do
!             
!             
!             return
!         end subroutine findAllIncidentVectors
!         
        
        subroutine findAllEulerAngles( v,g,theta0,phi0,theta_max,khat,nkhat , nolimitPhi )
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      first rotate by theta about x axis then by phi about y axis
    !*          khat = ( cos theta sin phi, - sin theta, cos theta cos phi )

    
            real(kind=real64),intent(in)                    ::      v
            real(kind=real64),dimension(3),intent(in)       ::      g   
            real(kind=real64),intent(in)                    ::      theta0,phi0
            real(kind=real64),intent(in)                    ::      theta_max         
            real(kind=real64),dimension(:,:),intent(out)    ::      khat
            integer,intent(out)                             ::      nkhat
            logical,intent(in),optional                     ::      nolimitPhi
            integer                         ::      ntheta
            real(kind=real64)               ::      theta_min,theta,pp,cost,sint,cosp,sinp 
            real(kind=real64),dimension(4)  ::      phi
            integer                         ::      nphi
            integer         ::      ii,jj
            
            
            
            integer                             ::      nProcs,rank,ierror
            
            call MPI_COMM_SIZE(MPI_COMM_WORLD, nProcs, ierror)
            call MPI_COMM_RANK(MPI_COMM_WORLD, rank, ierror)
            
             
            theta_min = max(0.0d0,theta0 - theta_max/2)
            
            if (rank == 0) write(*,fmt='(3(a,f12.3))') "Lib_DiffractionConditions::findAllEulerAngles info - checking theta range ",theta_min,":",theta_min+theta_max," (rad)"
            
            ntheta = size(khat,dim=2)/4
            nkhat = 0
            do ii = 1,ntheta
                theta = theta_min + theta_max*ii/ntheta
                cost = cos(theta) ; sint = sin(theta)
                if (present(nolimitPhi)) then
                    call findRotationsForBraggCondition2( v,g,cost,sint,phi0,phi,nphi,nolimitPhi)                
                else
                    call findRotationsForBraggCondition2( v,g,cost,sint,phi0,phi,nphi )
                end if
                do jj = 1,nphi
                    pp = phi(jj) ; cosp = cos(pp) ; sinp = sin(pp)
                    khat(:,nkhat+jj) = (/ -sinp, sint*cosp, cost*cosp /)
                end do
                nkhat = nkhat + nphi
                        
            end do    
                
                
            
          ! print *,"e dev parm all khat"
          ! do ii = 1,nkhat
          !     print *,ii,khat(:,ii),energy_deviation_parameter( v,khat(:,ii),g )
          ! end do
            
            
            return
        end subroutine findAllEulerAngles
        
        
        
    !---
    
    
        pure function getReciprocalLatticeVectors( A_cell ) result( B_cell )    
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      given the matrix A_cell whose columns are the unit cell vectors
    !*      return the matrix B_cell whose columns are the reciprocal lattice vectors
            real(kind=real64),dimension(3,3),intent(in)     ::      A_cell              !   lattice vectors of conventional unit cell
            real(kind=real64),dimension(3,3)                ::      B_cell              !   conventional reciprocal lattice vectors
            
            call inverse3Mat(transpose(A_cell),B_cell) 
            B_cell = 2*PI*B_cell  
            return
        end function getReciprocalLatticeVectors

            
        pure function reflection0( B, hkl ) result(g)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ 
            real(kind=real64),dimension(3,3),intent(in)     ::      B       !   reciprocal lattice vectors 2 pi A^-1            ( 1/A ) 
            integer,dimension(3),intent(in)                 ::      hkl
            real(kind=real64),dimension(3)                  ::      g
            
            g(1) = B(1,1)*hkl(1) + B(1,2)*hkl(2) + B(1,3)*hkl(3)
            g(2) = B(2,1)*hkl(1) + B(2,2)*hkl(2) + B(2,3)*hkl(3)
            g(3) = B(3,1)*hkl(1) + B(3,2)*hkl(2) + B(3,3)*hkl(3)
                                        
            return
        end function reflection0
        
        pure function reflection1( B, hkl ) result(g)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ 
            real(kind=real64),dimension(3,3),intent(in)     ::      B       !   reciprocal lattice vectors 2 pi A^-1            ( 1/A ) 
            real(kind=real64),dimension(3),intent(in)       ::      hkl
            real(kind=real64),dimension(3)                  ::      g
            
            g(1) = B(1,1)*hkl(1) + B(1,2)*hkl(2) + B(1,3)*hkl(3)
            g(2) = B(2,1)*hkl(1) + B(2,2)*hkl(2) + B(2,3)*hkl(3)
            g(3) = B(3,1)*hkl(1) + B(3,2)*hkl(2) + B(3,3)*hkl(3)
                                        
            return
        end function reflection1
    !---    
    
!             
!         subroutine permittedReflections( lattice,n,hkl,rho_in )
!     !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ 
!     !*      return the permitted reflections [hkl]
!     !*      such that (h)^2 + (k)^2 + (l)^2 <= rho^2
!     !*      given the lattice type
!             integer,intent(in)                  ::      lattice            
!             integer,intent(out)                 ::      n
!             integer,dimension(:,:),allocatable,intent(out)  ::      hkl
!             real(kind=real64),intent(in),optional        ::      rho_in
!              
!             integer             ::      mm,rr,pass
!             integer             ::      hp,kp,lp,hx,kx,lx
!             logical             ::      ok
!             real(kind=real64)   ::      rho
!             
!             rho = RHO_MAX ; if (present(rho_in)) rho = rho_in
!             
!             mm = ceiling( rho + 1.0d-12 )        !   how far to look in each direction
!             
!             do pass = 1,2                       !   first pass - count. second pass - store.
!             
!                 n = 0
!                 do lp = -mm,+mm
!                     do kp = -mm,+mm
!                         do hp = -mm,+mm
!                             rr = lp*lp + kp*kp + hp*hp
!                             if (rr > rho*rho + 1.0d-6) cycle
!                                                     
!                             select case(lattice)
!                                 case( LATTICE_FCC )
!                                     hx = mod( hp, 2 )
!                                     kx = mod( kp, 2 )
!                                     lx = mod( lp, 2 )
!                                     ok = ( hx*hx + kx*kx + lx*lx ) - ( hx*kx + kx*lx + lx*hx ) == 0     !   true if all odd or all even.
!                                 case( LATTICE_BCC )
!                                     ok = mod( hp+kp+lp,2 ) == 0     !   true if sum is even
!                                 case default
!                                     ok = .true.
!                             end select
!                                                         
!                             if (ok) then
!                                 n = n + 1 
!                                 if (pass==2) hkl(:,n) = (/ hp,kp,lp /)
!                             end if
!                             
!                         end do
!                     end do
!                 end do
!                 if (pass==1) allocate(hkl(3,n))
!                     
!             end do        
!             return
!         end subroutine permittedReflections
!                                   
    !---    
    
        subroutine find_gvectors_permittedReflections( latt,B,hkl,rho, n,ghkl,hklp  )
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ 
    !*      find the g-vectors corresponding to permitted reflections centred on [hkl]
            type(Lattice),intent(in)            ::      latt            
            real(kind=real64),intent(in)        ::      rho
            real(kind=real64),dimension(3,3),intent(in)     ::      B
            integer,dimension(3),intent(in)     ::      hkl
            integer,intent(out)                 ::      n
            real(kind=real64),dimension(:,:),allocatable,intent(out)  ::      ghkl
            integer,dimension(:,:),allocatable,intent(out)            ::      hklp
            
            !integer,dimension(:,:),allocatable  ::      hklp
            integer                             ::      ii,i000,ihkl
            integer,dimension(3)                ::      hhkkll
            real(kind=real64),dimension(3)      ::      gswp
            
            call permittedReflections( latt,n,hklp,rho ) !   note: this allocates hklp(3,n)
             
            i000 = 0
            ihkl = 0
            allocate( ghkl(3,n) )
            !   generate set of g-vectors centred on [hkl]
!             do ii = 1,n
!                 if (hklp(1,ii)*hklp(1,ii) + hklp(2,ii)*hklp(2,ii) + hklp(3,ii)*hklp(3,ii)==0) ihkl = ii
!                 hhkkll(1:3) = hklp(1:3,ii)+hkl(1:3)
!                 if (hhkkll(1)*hhkkll(1) + hhkkll(2)*hhkkll(2) + hhkkll(3)*hhkkll(3)==0) i000 = ii
!                 ghkl(:,ii) = reflection( B,hklp(1:3,ii)+hkl(1:3) )
!             end do
            !   generate set of g-vectors centred on [000]
            do ii = 1,n
                if (hklp(1,ii)*hklp(1,ii) + hklp(2,ii)*hklp(2,ii) + hklp(3,ii)*hklp(3,ii)==0) i000 = ii
                hhkkll(1:3) = hklp(1:3,ii)-hkl(1:3)
                if (hhkkll(1)*hhkkll(1) + hhkkll(2)*hhkkll(2) + hhkkll(3)*hhkkll(3)==0) ihkl = ii
                ghkl(:,ii) = reflection( B,hklp(1:3,ii)  )
            end do
             
            
            
        
        !---    debugging checks only         
            if (i000 /= 0) then
                gswp(1:3) = ghkl(1:3,1)
                ghkl(1:3,1) = ghkl(1:3,i000)  
                ghkl(1:3,i000) = gswp(1:3)
                
                hhkkll(1:3) = hklp(1:3,1)
                hklp(1:3,1) = hklp(1:3,i000)
                hklp(1:3,i000) = hhkkll(1:3)
                
            end if
            
            if (ihkl /= 0) then            
                gswp(1:3) = ghkl(1:3,2)
                ghkl(1:3,2) = ghkl(1:3,ihkl)  
                ghkl(1:3,ihkl) = gswp(1:3)
                
                hhkkll(1:3) = hklp(1:3,1)
                hklp(1:3,1) = hklp(1:3,ihkl)
                hklp(1:3,ihkl) = hhkkll(1:3)
            end if
                        
            return
        end subroutine find_gvectors_permittedReflections
        
        
    !---    
    
        pure function score( eps,gamma ) result( s )
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      A Lorentzian score function based on energy deviation 
    !*          s = 1/(gamma^2 + eps^2)
            real(kind=real64),intent(in)            ::      eps         !   energy deviation parameter
            real(kind=real64),intent(in)            ::      gamma       !   small energy deviation parameter
                        
            real(kind=real64)                       ::      s   
            
            s = gamma*gamma + eps*eps
            s = 1/s
            
            return
        end function score
        
!         pure function scoreForIncidentVector0(v,khat,gamma,ghkl) result(s)
!     !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
!     !*      find the total score for all considered reflections
!          
!             real(kind=real64),intent(in)                    ::      v
!             real(kind=real64),dimension(3),intent(in)       ::      khat 
!             real(kind=real64),intent(in)                    ::      gamma
!             real(kind=real64),dimension(:,:),intent(in)     ::      ghkl
!             
!             real(kind=real64)                               ::      s
!             real(kind=real64)       ::      eps,score0
!             integer                 ::      ii,nn
!             nn = size(ghkl,dim=2)
!             
!             
!             score0 = score( 0.0d0,gamma ) 
!             s = 0
!             do ii = 1,nn
!                 eps = energy_deviation_parameter( v,khat,ghkl(:,ii) )
!                 s = s + score( eps,gamma )
!             end do
!             s = s / score0 - 2
!             
!             return
!         end function scoreForIncidentVector0
        
        pure function scoreForIncidentVector1(v,khat,ghkl , xi_g,foil_thickness) result(s)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ 
    !*      find the total score for all considered reflections
         
            real(kind=real64),intent(in)                    ::      v
            real(kind=real64),dimension(3),intent(in)       ::      khat 

            real(kind=real64),dimension(:,:),intent(in)     ::      ghkl
            
            real(kind=real64),dimension(:),intent(in)       ::      xi_g
            real(kind=real64),intent(in)                    ::      foil_thickness
           
            real(kind=real64)                               ::      s
            !real(kind=real64)       ::      eps,sg
            integer                 ::      ii!,nn
            
            
            
!             s = 0
!             nn = size(ghkl,dim=2)            
!             do ii = 3,nn            !   because ignore score value of [000] and [hkl] which are first two permitted vectors
!                 eps = energy_deviation_parameter( v,khat,ghkl(:,ii) )
!                 sg = energyToLengthDeviationParameter(v,eps)
!                 s = s + intensityAtDepthZ( xi_g(1),xi_g(ii),sg, foil_thickness )
!             end do

            s = 0
            do ii = 3,size(ghkl,dim=2)            !   because ignore score value of [000] and [hkl] which are first two permitted vectors
                s = s + scoreForIncidentVector2( v,khat,ghkl(:,ii) , xi_g(1),xi_g(ii), foil_thickness )
            end do
            
            return
        end function scoreForIncidentVector1
        
        pure function scoreForIncidentVector2( v,khat,ghkl , xi_0,xi_g,foil_thickness) result(s)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      find the total score for all considered reflections
         
            real(kind=real64),intent(in)                    ::      v
            real(kind=real64),dimension(3),intent(in)       ::      khat 
 
            real(kind=real64),dimension(:),intent(in)       ::      ghkl
            
            real(kind=real64),intent(in)                    ::      xi_0,xi_g
            real(kind=real64),intent(in)                    ::      foil_thickness
           
            real(kind=real64)                               ::      s
            real(kind=real64)       ::      eps,sg

            eps = energy_deviation_parameter( v,khat,ghkl )
            sg  = energyToLengthDeviationParameter(v,eps)
            s   = intensityAtDepthZ( xi_0,xi_g,sg, foil_thickness )
            
            return
        end function scoreForIncidentVector2
        
    !---
    
        

        pure subroutine inverse3Mat(M,N)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
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
    !*      returns the determinant of real matrix M
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
        
        pure function cross_product(x,y) result(z)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
            real(kind=real64),dimension(3),intent(in)       ::      x,y
            real(kind=real64),dimension(3)                  ::      z
            z(1) = x(2)*y(3) - x(3)*y(2)
            z(2) = x(3)*y(1) - x(1)*y(3)
            z(3) = x(1)*y(2) - x(2)*y(1)
            return
        end function cross_product
        
         
    end module Lib_DiffractionConditions
        
    
! !   gfortran -ffree-line-length-256 ${MYF90LIB}/NBAX_StringTokenizers.f90 ${MYF90LIB}/Lib_CommandLineArguments.f90 ${MYF90LIB}/Lib_RotationMatrices.f90 ${MYF90LIB}/Lib_RelativisticElectrons.f90 Lib_ReadExtinctionDistances.f90 Lib_DiffractionConditions.f90 -o testLib_DiffractionConditions.exe        
!     
!     program testLib_DiffractionConditions
! !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
!         use Lib_CommandLineArguments
!         use Lib_ReadExtinctionDistances
!         use Lib_DiffractionConditions
!         use Lib_RelativisticElectrons
!         use iso_fortran_env
!         implicit none
!         
!         type(CommandLineArguments)      ::      cla
!         character(len=256)                          ::      extinction_distances_filename = ""
!         real(kind=real64)                           ::      EV = 150000d0
!         real(kind=real64),parameter                 ::      a0 = 3.1652d0
!         real(kind=real64)                           ::      Z = 1000.0d0
!         real(kind=real64),dimension(3,3),parameter  ::      AA = reshape( (/1,0,0 , 0,1,0 , 0,0,1/),(/3,3/) )*a0
!         integer,parameter                           ::      lattice = LATTICE_BCC
!         integer,dimension(3)                        ::      HKL = (/ 4,4,0 /)
!         real(kind=real64),dimension(3)              ::      zz = (/ 0.0,0.0,1.0 /)
!         
!         
!     !---    read extinction distances from file
!         integer,dimension(:,:),allocatable          ::      hkl_permitted
!         real(kind=real64),dimension(:),allocatable  ::      xi_g
!         integer                                     ::      nReflections
!         
!     !---    find intensity of peaks near Bragg condition
!         integer                                     ::      nNearBragg
!         integer,dimension(:,:),allocatable          ::      hkl_nearBragg
!         real(kind=real64),dimension(:),allocatable  ::      eps_nearBragg    
!         
!         
!         real(kind=real64)                   ::      vv,ss,pp
!         real(kind=real64),dimension(3,3)    ::      RR
!         character(len=256)          ::      dummy
!         logical                     ::      ok
!         integer                     ::      ii
!         
!         
!     !---    read command line arguments
!         cla = CommandLineArguments_ctor(10)  
!          
!         call setProgramDescription( cla, "testLib_DiffractionConditions.exe" )
!         call setProgramVersion( cla, "0.1.0" )               
!         
!          
!         
!         call get( cla,"f",extinction_distances_filename ,LIB_CLA_REQUIRED," extinction distances file" )                                                        
!         call get( cla,"V",EV ,LIB_CLA_REQUIRED," electron acceleration voltage" )                                                        
!         ii=3 ; call get( cla,"hkl",hkl ,ii,LIB_CLA_REQUIRED," reflection [hkl]" )                                                        
!         ii=3 ; call get( cla,"zz",zz , ii,LIB_CLA_OPTIONAL," zone axis" )                                                        
!         call get( cla,"Z",Z ,LIB_CLA_OPTIONAL," foil thickness" )       
!         pp = THETA_MAX*180.0d0/3.141592654d0
!         call get( cla,"theta",pp ,LIB_CLA_OPTIONAL," max tilt angle " )                 
!         THETA_MAX = pp*3.141592654d0/180.0d0
!         
!         
!         
!         call report(cla)
!         if (.not. allRequiredArgumentsSet(cla)) stop
!         if (hasHelpArgument(cla)) stop
!         call delete(cla)
!         
!         
!         
!        ! print *,"usage:   testLib_DiffractionConditions.exe filename EV h,k,l z1,z2,z3 theta_max(deg) foil_thickness"
!        ! call get_command_argument( 1, extinction_distances_filename )
!         print *,"extinction distances file """//trim(extinction_distances_filename)//""""
!                                
!       !  call get_command_argument( 2, dummy )
!       !  read(dummy,*) EV
!         print *,"electron acceleration voltage ",EV," (V)"
!                                
!       !  call get_command_argument( 3, dummy )
!       !  read(dummy,*) hkl
!         print *,"reflection [hkl] ",hkl
!                                
!       !  call get_command_argument( 4, dummy )
!       !  read(dummy,*) zz
!         print *,"zone axis ",zz
!         print *,"lattice symmetry ",lattice
!       !  print *,"foil thickness ",Z," (A)"
!         
!       !  call get_command_argument( 5, dummy )
!       !  read(dummy,*) pp
!       !  THETA_MAX = pp*3.141592654d0/180.0d0
!         print *,"max tilt angle ",pp," (deg) = ",THETA_MAX," (rad)"
!         
!         
!       !  call get_command_argument( 6, dummy )
!      !   read(dummy,*) Z        
!         print *,"foil thickness ",Z," (A)"
!         
!         print *,""
!         print *,""
!         
!         
!         !call listPermittedReflections( lattice,20.0d0 )
!         
!         
!         !call permittedReflections( lattice,nReflections,hkl_permitted )
!         
!         vv = velocity( EV )
!         
!          
!         call findBestRotationMatrix( vv,zz,AA,lattice, hkl , RR , extinction_distances_filename,Z, nNearBragg,hkl_nearBragg,eps_nearBragg )
! 
!             allocate(xi_g(nNearBragg))
!             call readExtinctionDistance( extinction_distances_filename,hkl_nearBragg,xi_g, ok )
!             print *,"readExtinctionDistance() info - ",ok
! 
!             print *,"findBestRotationMatrix info - reflections near Bragg condition "
!             write (*,fmt='(a12,100a12)') "reflection","eps_g (eV)","s_g (1/A)","xi_g (A)","|phig(Z)|^2"
!             do ii = 1,nNearBragg                
!                 ss = energyToLengthDeviationParameter( vv,eps_nearBragg(ii) )
!                 pp = intensityAtDepthZ( xi_g(1),xi_g(ii),ss, Z )        !   note: hkl(1) = [000]
!                 write (*,fmt='(3i4,100f12.5)') hkl_nearBragg(:,ii),eps_nearBragg(ii),ss,xi_g(ii),pp
!             end do
! 
!            
!            
!             
!             
!             
!         print *,""
!         print *,"done"
!         print *,""
!     end  program testLib_DiffractionConditions