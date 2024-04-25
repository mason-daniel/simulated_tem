
    module Lib_DynamicalTwoBeam
!---^^^^^^^^^^^^^^^^^^^^^^^^^^^
        use Lib_RelativisticElectrons
        use iso_fortran_env 
        implicit none
        private


        real(kind=real64),private,parameter     ::      PI = 3.1415926539d0
        real(kind=real64),private,parameter     ::      HBAR = 1.054571817d-34     !    J.s
        real(kind=real64),private,parameter     ::      ME = 9.10938356d-31        !    Kg
        real(kind=real64),private,parameter     ::      KEV = 1.602176634d-16      !    J

        public      ::      dynamicalTwoBeamIntegration
        
        interface       dynamicalTwoBeamIntegration
            module procedure        dynamicalTwoBeamIntegration_disp
            module procedure        dynamicalTwoBeamIntegration_phase
            module procedure        dynamicalTwoBeamIntegration_strain
            module procedure        dynamicalTwoBeamIntegration_pos
            module procedure        dynamicalTwoBeamIntegration_phasecol
        end interface
        
    contains
!---^^^^^^^^


    
        subroutine dynamicalTwoBeamIntegration_phase( E,xig,xi0,g,x,xmax , phig2 )
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      integrate dynamical two beam equations through depth zmax
    !*      given the phase factor x
    !*      and the scattering vector g
    !*      and two parameters electron beam energy E and extinction distance xig
    
            real(kind=real64),intent(in)                        ::      E           !   electron beam energy (keV)
            real(kind=real64),intent(in)                        ::      xig         !   extinction distance (A)   
            real(kind=real64),intent(in)                        ::      xi0         !   chaaracteristic length of average crystal potential (A)   
            real(kind=real64),dimension(3),intent(in)           ::      g           !   (1/A)
            complex(kind=real64),dimension(:,:,0:),intent(in)   ::      x
            real(kind=real64),dimension(3),intent(in)           ::      xmax        !   foil extent (A)
            
            real(kind=real64),dimension(size(x,dim=1),size(x,dim=2)),intent(inout)    ::      phig2
                        
            
            integer                 ::      nx,ny,nz
            integer                 ::      ix,iy
            real(kind=real64)       ::      aa,bg,b0
            real(kind=real64)       ::      vv,eps,kk , deltaz,dphig2 !,xi
            
!            complex(kind=real64)    ::      xbar
            
        !---    find the number of x,y,z divisions
            nx = size(x,dim=1)
            ny = size(x,dim=2)
            nz = size(x,dim=3)-1
            if (mod(nz,2) /= 0) then
                print *,"Lib_DynamicalTwoBeam::dynamicalTwoBeamIntegration_phase() error - RK4 algorithm wants to find a derivative at z=0 and z=zmax"
                print *,"and so requires an odd number of position points in z, found ",nz
               ! phig2 = 0.0d0
                return
            end if
            
            
        !---    compute the coefficients needed for the column integration
            vv = velocity( E*1000 )             !   velocity of electron A/fs
            kk = 2*PI/wavelength( E*1000 )      !   wave vector of electron (A^-1)
            eps = ( 2*kk*g(3) + g(1)*g(1) + g(2)*g(2) + g(3)*g(3) )     !   wavevector deviation^2 (A^-2)
            deltaz = xmax(3)/nz                 !   length of z division
            eps = eps/(2*kk)            
!             xi = xig*nz/xmax(3)
!             b0 = (PI/xi0)
!             bb = -eps*nz/xmax(3)
                
            bg = (PI/xig) * deltaz
            b0 = (PI/xi0) * deltaz
            aa = -eps     * deltaz
            
            print *,"Lib_DynamicalTwoBeam::dynamicalTwoBeamIntegration_phase() info - "
            write(*,fmt='(a,g16.8)') "    electron velocity       (A/fs) ",vv
            write(*,fmt='(a,g16.8)') "    electron wavevector      (1/A) ",kk
            write(*,fmt='(a,3g16.8)')"    scattering wavevector    (1/A) ",g
            write(*,fmt='(a,3g16.8)')"    z-division                 (A) ",deltaz
            write(*,fmt='(a,g16.8)') "    PI/xi0                (1/cell) ",b0
            write(*,fmt='(a,g16.8)') "    PI/xig                (1/cell) ",bg
            write(*,fmt='(a,g16.8)') "    - eps/(hbar v)        (1/cell) ",aa
            
             
           ! do ix = 0,nz
           !     xbar = sum(x(:,:,ix)) / (nx*ny)
           !     print *,"mean x ",ix,xbar
           ! end do
           ! print *,"xig/xi0",xig/xi0
             
            
            
            do iy = 1,ny
                do ix = 1,nx
                
                    call dynamicalTwoBeamIntegrationColumn_phase( aa,b0,bg,x(ix,iy,:),dphig2  )     
                    phig2(ix,iy) = phig2(ix,iy) + dphig2
                    ! stop
               end do
            end do
            
            return
            
            
            
        end subroutine dynamicalTwoBeamIntegration_phase
            

        subroutine dynamicalTwoBeamIntegration_phasecol( E,xig,xi0,g,x,xmax , phig2 )
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      integrate dynamical two beam equations through depth zmax
    !*      given the phase factor x
    !*      and the scattering vector g
    !*      and two parameters electron beam energy E and extinction distance xig
    
            real(kind=real64),intent(in)                        ::      E           !   electron beam energy (keV)
            real(kind=real64),intent(in)                        ::      xig         !   extinction distance (A)   
            real(kind=real64),intent(in)                        ::      xi0         !   chaaracteristic length of average crystal potential (A)   
            real(kind=real64),dimension(3),intent(in)           ::      g           !   (1/A)
            complex(kind=real64),dimension(0:),intent(in)       ::      x
            real(kind=real64),dimension(3),intent(in)           ::      xmax        !   foil extent (A)
            
            real(kind=real64),intent(inout)                     ::      phig2
                        
            
            integer                 ::      nz
            real(kind=real64)       ::      aa,bg,b0
            real(kind=real64)       ::      vv,eps,kk , deltaz,dphig2 !,xi
            
!            complex(kind=real64)    ::      xbar
            
        !---    find the number of x,y,z divisions
            nz = size(x)-1
            if (mod(nz,2) /= 0) then
                print *,"Lib_DynamicalTwoBeam::dynamicalTwoBeamIntegration_phasecol() error - RK4 algorithm wants to find a derivative at z=0 and z=zmax"
                print *,"and so requires an odd number of position points in z, found ",nz
               ! phig2 = 0.0d0
                return
            end if
            
            
        !---    compute the coefficients needed for the column integration
            vv = velocity( E*1000 )             !   velocity of electron A/fs
            kk = 2*PI/wavelength( E*1000 )      !   wave vector of electron (A^-1)
            eps = ( 2*kk*g(3) + g(1)*g(1) + g(2)*g(2) + g(3)*g(3) )     !   wavevector deviation^2 (A^-2)
            deltaz = xmax(3)/nz                 !   length of z division
            eps = eps/(2*kk)            
                
            bg = (PI/xig) * deltaz
            b0 = (PI/xi0) * deltaz
            aa = -eps     * deltaz
            
!             print *,"Lib_DynamicalTwoBeam::dynamicalTwoBeamIntegration_phase() info - "
!             write(*,fmt='(a,g16.8)') "    electron velocity       (A/fs) ",vv
!             write(*,fmt='(a,g16.8)') "    electron wavevector      (1/A) ",kk
!             write(*,fmt='(a,3g16.8)')"    scattering wavevector    (1/A) ",g
!             write(*,fmt='(a,3g16.8)')"    z-division                 (A) ",deltaz
!             write(*,fmt='(a,g16.8)') "    PI/xi0                (1/cell) ",b0
!             write(*,fmt='(a,g16.8)') "    PI/xig                (1/cell) ",bg
!             write(*,fmt='(a,g16.8)') "    - eps/(hbar v)        (1/cell) ",aa
                         
            
                
            call dynamicalTwoBeamIntegrationColumn_phase( aa,b0,bg,x(:),dphig2  )     
            phig2 = phig2 + dphig2
    
            return
                                   
        end subroutine dynamicalTwoBeamIntegration_phasecol
            
    
        subroutine dynamicalTwoBeamIntegration_pos( E,xig,xi0,g,r,xmax , phi2 )
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      integrate dynamical two beam equations through depth zmax
    !*      given the displacement field u
    !*      and the scattering vector g
    !*      and two parameters electron beam energy E and extinction distance xig
    
            real(kind=real64),intent(in)                        ::      E           !   electron beam energy (keV)
            real(kind=real64),intent(in)                        ::      xig         !   extinction distance (A)   
            real(kind=real64),intent(in)                        ::      xi0         !   chaaracteristic length of average crystal potential (A)   
            real(kind=real64),dimension(3),intent(in)           ::      g           !   (1/A)
            real(kind=real64),dimension(:,:,:,0:),intent(in)    ::      r
            real(kind=real64),dimension(3),intent(in)           ::      xmax        !   foil extent (A)
            
            real(kind=real64),dimension(size(r,dim=2),size(r,dim=3)),intent(out)                  ::      phi2
                        
            !real(kind=real64),dimension(3,0:size(r,dim=4)-1)    ::      rcol
            real(kind=real64),dimension(:,:),allocatable    ::      rcol
            
            integer                 ::      nx,ny,nz
            integer                 ::      ix,iy
            real(kind=real64)       ::      aa,bb,b0
            real(kind=real64)       ::      vv,eps,kk,xi
        
            
        !---    find the number of x,y,z divisions
            nx = size(r,dim=2)
            ny = size(r,dim=3)
            nz = size(r,dim=4)-1
            if (mod(nz,2) /= 0) then
                print *,"Lib_DynamicalTwoBeam::dynamicalTwoBeamIntegration_pos() error - RK4 algorithm wants to find a derivative at z=0 and z=zmax"
                print *,"and so requires an odd number of position points in z, found ",nz
                phi2 = 0.0d0
                return
            end if
            
            
        !---    compute the coefficients needed for the column integration
            vv = velocity( E*1000 )             !   velocity of electron A/fs
            kk = 2*PI/wavelength( E*1000 )      !   wave vector of electron (A^-1)
            eps = ( 2*kk*g(3) + g(1)*g(1) + g(2)*g(2) + g(3)*g(3) )     !   wavevector deviation^2 (A^-2)
            eps = eps/(2*kk)            
            xi = xig*nz/xmax(3)
            b0 = PI/xi0
            bb = -eps*nz/xmax(3)
                
            aa = PI/xi
            
            print *,"Lib_DynamicalTwoBeam::dynamicalTwoBeamIntegration_pos() info - "
            write(*,fmt='(a,g16.8)') "    electron velocity        (m/s) ",vv
            write(*,fmt='(a,g16.8)') "    electron wavevector      (1/A) ",kk
            write(*,fmt='(a,3g16.8)')"    scattering wavevector    (1/A) ",g
            write(*,fmt='(a,g16.8)') "    PI/xi0                (1/cell) ",b0
            write(*,fmt='(a,g16.8)') "    - eps/(hbar v)        (1/cell) ",bb
            write(*,fmt='(a,g16.8)') "    PI/xig                (1/cell) ",aa
            
            allocate( rcol(3,-4:nz) )
            do iy = 1,ny
                do ix = 1,nx
                    !ix = nx/2 ; iy = ny/2
                    !rcol(1:3,0:nz) = r(1:3,ix,iy,0:nz)
                    call extractColumnAndFix( r,ix,iy, rcol )
                    call dynamicalTwoBeamIntegrationColumn_disp( aa,b0,bb,g,rcol, phi2(ix,iy) )    !       note- same equation position or displacement
               end do
           end do
            
            return
            
        contains
    !---^^^^^^^^
    
            subroutine extractColumnAndFix( r,ix,iy, rcol )
        !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
        !*      it is much better for RK4 if displacements are smoothly varying 
        
        !       write f(x) = a + bx + cx^2 + dx^3
        !        with f(0) = a
        !            f'(0) = b
        !            f(-3) = 0
        !           f'(-3) = 0
                real(kind=real64),dimension(:,:,:,0:),intent(in)    ::      r
                integer,intent(in)                                  ::      ix,iy
                real(kind=real64),dimension(:,-4:),intent(out)      ::      rcol
               
                real(kind=real64),dimension(3)           ::      f0,f1,f2
                
                real(kind=real64),dimension(6,6),parameter  ::      IM = reshape( (/    & 
                                    512,0,0,80,30,3,                                    &
                                    0,512,0,-192,-64,-6,                                &
                                    0,0,256,192,48,4,                                   &
                                    0,0,0,-80,-30,-3,                                   &
                                    0,0,0,-128,-56,-6,                                  &
                                    0,0,0,-64,-32,-4 /),(/6,6/) )/512.0d0
                
                real(kind=real64),dimension(6)      ::      aa
                                                    
                integer             ::      ii,nz
                
                rcol = 0
                f0 = r(1:3,ix,iy,0)
                f1 = ( -3*r(1:3,ix,iy,0) + 4*r(1:3,ix,iy,1) - r(1:3,ix,iy,2) )/2
                f2 = ( 2*r(1:3,ix,iy,0) - 5*r(1:3,ix,iy,1) + 4*r(1:3,ix,iy,2) - r(1:3,ix,iy,3) )
                
                do ii = 1,3
                    aa(1:3) = (/f0(ii),f1(ii),f2(ii)/)
                    aa(4:6) = iM(4:6,1)*f0(ii) + iM(4:6,2)*f1(ii) + iM(4:6,3)*f2(ii)  
                    
                    rcol(ii,-1) = aa(1) -   aa(2) +   aa(3) -    aa(4) +    aa(5) -     aa(6)
                    rcol(ii,-2) = aa(1) - 2*aa(2) + 4*aa(3) -  8*aa(4) + 16*aa(5) -  64*aa(6)
                    rcol(ii,-3) = aa(1) - 3*aa(2) + 9*aa(3) - 27*aa(4) + 81*aa(5) - 243*aa(6)
                end do            
                
                
               ! rcol(1:3,-1) = 9*(3*f0 - 2*f1)/32
               ! rcol(1:3,-2) =   (  f0 -   f1)/2
               ! rcol(1:3,-3) =   (5*f0 - 6*f1)/32
               ! rcol(1:3,-4) = 0
                
                nz = size(r,dim=4)-1
                    
                do ii = 0,nz
                    rcol(1:3,ii) = r(1:3,ix,iy,ii)    
                end do
                
               ! f0 = r(1:3,ix,iy,nz)
               ! f1 = ( 3*r(1:3,ix,iy,nz) - 4*r(1:3,ix,iy,nz-1) + r(1:3,ix,iy,nz-2) )/2
               ! 
               ! rcol(1:3,nz+1) = 9*(3*f0 + 2*f1)/32
               ! rcol(1:3,nz+2) =   (  f0 +   f1)/2
               ! rcol(1:3,nz+3) =   (5*f0 + 6*f1)/32
               ! rcol(1:3,nz+4) = 0
               ! 
                return
            end subroutine extractColumnAndFix
            
            
            
        end subroutine dynamicalTwoBeamIntegration_pos
            

    
        subroutine dynamicalTwoBeamIntegration_disp( E,xig,g,u,xmax , phi2 )
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      integrate dynamical two beam equations through depth zmax
    !*      given the displacement field u
    !*      and the scattering vector g
    !*      and two parameters electron beam energy E and extinction distance xig
    
            real(kind=real64),intent(in)                        ::      E           !   electron beam energy (keV)
            real(kind=real64),intent(in)                        ::      xig         !   extinction distance (A)   
            real(kind=real64),dimension(3),intent(in)           ::      g           !   (1/A)
            real(kind=real64),dimension(:,:,:,0:),intent(in)    ::      u
            real(kind=real64),dimension(3),intent(in)           ::      xmax        !   foil extent (A)
            
            real(kind=real64),dimension(size(u,dim=2),size(u,dim=3)),intent(out)                  ::      phi2
                        
            real(kind=real64),dimension(3,0:size(u,dim=4)-1)    ::      ucol
            
            integer                 ::      nx,ny,nz
            integer                 ::      ix,iy
            real(kind=real64)       ::      aa,bb
            real(kind=real64)       ::      vv,eps,kk,xi
        
            
        !---    find the number of x,y,z divisions
            nx = size(u,dim=2)
            ny = size(u,dim=3)
            nz = size(u,dim=4)-1
            if (mod(nz,2) /= 0) then
                print *,"Lib_DynamicalTwoBeam::dynamicalTwoBeamIntegration_disp() error - RK4 algorithm wants to find a derivative at z=0 and z=zmax"
                print *,"and so requires an odd number of displacement points in z, found ",nz
                phi2 = 0.0d0
                return
            end if
            
            
        !---    compute the coefficients needed for the column integration
            vv = sqrt( 2*E*KEV/ME )             !   velocity of electron ( non-relativistic approximation ) m/s
            kk = (ME*vv/HBAR)*1d-10             !   wave vector of electron (A^-1)
            eps = ( 2*kk*g(3) + g(1)*g(1) + g(2)*g(2) + g(3)*g(3) )     !   wavevector deviation^2 (A^-2)
            eps = eps/(2*kk)            
                           print *,"eps",eps
            xi = xig*nz/xmax(3)
                           print *,"HBAR*v ",HBAR*vv
            bb = -eps*nz/xmax(3)
                
            aa = PI/xi
            
            print *,"Lib_DynamicalTwoBeam::dynamicalTwoBeamIntegration_disp() info - "
            write(*,fmt='(a,g16.8)') "    electron velocity     (m/s) ",vv
            write(*,fmt='(a,g16.8)') "    electron wavevector   (1/A) ",kk
            write(*,fmt='(a,3g16.8)')"    scattering wavevector (1/A) ",g
            write(*,fmt='(a,g16.8)') "    eps/(hbar v)       (1/cell) ",-bb
            write(*,fmt='(a,g16.8)') "    PI/xi              (1/cell) ",aa
            
            
            do iy = 1,ny
                do ix = 1,nx
                    ucol(1:3,0:nz) = u(1:3,ix,iy,0:nz)
                    call dynamicalTwoBeamIntegrationColumn_disp( aa,bb,0.0d0,g,ucol, phi2(ix,iy) )
                end do
            end do
            
            return
        end subroutine dynamicalTwoBeamIntegration_disp
            
            
        subroutine dynamicalTwoBeamIntegration_strain( E,xig,g,dudz,xmax , phi2 ,strain)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      integrate dynamical two beam equations through depth zmax
    !*      given the displacement field u
    !*      and the scattering vector g
    !*      and two parameters electron beam energy E and extinction distance xig
    
            real(kind=real64),intent(in)                        ::      E           !   electron beam energy (keV)
            real(kind=real64),intent(in)                        ::      xig         !   extinction distance (A)   
            real(kind=real64),dimension(3),intent(in)           ::      g           !   (1/A)
            real(kind=real64),dimension(:,:,:,0:),intent(in)    ::      dudz     
            real(kind=real64),dimension(3),intent(in)           ::      xmax        !   foil extent (A)
            
            real(kind=real64),dimension(size(dudz,dim=2),size(dudz,dim=3)),intent(out)                  ::      phi2
            
            logical,intent(in)      ::      strain
                        
            real(kind=real64),dimension(3,0:size(dudz,dim=4)-1)    ::      dudzcol
            
            integer                 ::      nx,ny,nz
            integer                 ::      ix,iy
            real(kind=real64)       ::      aa,bb
            real(kind=real64)       ::      vv,eps,kk,xi
        
            
            if (.not. strain) then
                call dynamicalTwoBeamIntegration_disp( E,xig,g,dudz,xmax , phi2 )           !   dudz is really just u
                return
            end if
            
        !---    find the number of x,y,z divisions
            nx = size(dudz,dim=2)
            ny = size(dudz,dim=3)
            nz = size(dudz,dim=4)-1
            if (mod(nz,2) /= 0) then
                print *,"Lib_DynamicalTwoBeam::dynamicalTwoBeamIntegration_strain() error - RK4 algorithm wants to find a derivative at z=0 and z=zmax"
                print *,"and so requires an odd number of displacement points in z, found ",nz
                phi2 = 0.0d0
                return
            end if
            
            
        !---    compute the coefficients needed for the column integration
            vv = sqrt( 2*E*KEV/ME )             !   velocity of electron ( non-relativistic approximation ) m/s
            kk = (ME*vv/HBAR)*1d-10             !   wave vector of electron (A^-1)
            eps = ( 2*kk*g(3) + g(1)*g(1) + g(2)*g(2) + g(3)*g(3) )     !   wavevector deviation^2 (A^-2)
            eps = eps/(2*kk)            
                           print *,"eps",eps
            xi = xig*nz/xmax(3)
                           print *,"HBAR*v ",HBAR*vv
            bb = -eps*nz/xmax(3)
                
            aa = PI/xi
            
            print *,"Lib_DynamicalTwoBeam::dynamicalTwoBeamIntegration_strain() info - "
            write(*,fmt='(a,g16.8)') "    electron velocity     (m/s) ",vv
            write(*,fmt='(a,g16.8)') "    electron wavevector   (1/A) ",kk
            write(*,fmt='(a,3g16.8)')"    scattering wavevector (1/A) ",g
            write(*,fmt='(a,g16.8)') "    eps/(hbar v)       (1/cell) ",-bb
            write(*,fmt='(a,g16.8)') "    PI/xi              (1/cell) ",aa
            
            
            do iy = 1,ny
                do ix = 1,nx
                    dudzcol(1:3,0:nz) = dudz(1:3,ix,iy,0:nz)
                    call dynamicalTwoBeamIntegrationColumn_strain( aa,bb,g,dudzcol, phi2(ix,iy) )
                end do
            end do
            
            return
        end subroutine dynamicalTwoBeamIntegration_strain
            
    
        subroutine dynamicalTwoBeamIntegrationColumn_phase( a,b0,bg,x , phig2 )
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      integrate dynamical two beam equations through depth zmax
    !*      given the displacement field u(z) 
    !*      and the scattering vector g ( in reciprocal cell units )
    !*      and three parameters  b0 = pi / xi_0 , bg = pi / xi_g
    !*                            a  = - eps / (h v)
    !*      with xi_g = pi h v / | U_g |
    !*      eps = hbar^2 ( k+g )^2 / (2m) - hbar^2 k^2/(2m)
    !*      v = hbar k / m
    !*
    !*      return the intensity squared phi2
    !*
    !*      phi_0' = i b0 phi_0 + i bg exp[ i g.u ] phi_g
    !*      phi_g' = i (a + b0) phi_g + i bg exp[ -i g.u ] phi_0
     
    !*      on input x = exp[ i g.u ]
     
            real(kind=real64),intent(in)                    ::      a,bg,b0
            complex(kind=real64),dimension(0:),intent(in)   ::      x
            
            real(kind=real64),intent(out)                   ::      phig2
            
            
            
            integer             ::      nz
            real(kind=real64)   ::      dz!,gdotu
            !real(kind=real64)   ::      lambda = 0.10d0
            integer             ::      ii
            
            
            real(kind=real64)   ::      phi0r ,phi0i ,phigr ,phigi          !   these are the field evaluated at z ( real and imaginary parts )
            real(kind=real64)   ::      dphi0r,dphi0i,dphigr,dphigi         !   derivative of field evalulated at z,z+dz/2,z+dz
            real(kind=real64)   ::      phi0r1,phi0i1,phigr1,phigi1         !   temporary 1 needed for RK4 integration
            real(kind=real64)   ::      phi0r2,phi0i2,phigr2,phigi2         !   temporary 2 needed for RK4 integration
            
            real(kind=real64)   ::      phimag , bgs,bgc
            
        !---    find the number of z divisions : note this should be an odd number ( checked in area() call )
            nz = size(x)-1
            dz = 2.0d0 / nz      !   note we use double steps because the derivative is computed at z,z+dz/2,z+dz
            
                  
        !    next initialise the propagating fields - wlog set real component of phi0 to 1, everything else zero
            phi0r = 1.0d0 ; phi0i = 0.0d0  ; phigr = 0.0d0 ; phigi = 0.0d0
            
            
        !   now integrate in z        
            !print *,"column disp ",0,phi0r,phi0i,phigr,phigi,(phi0r+phigr)**2 + (phi0i+phigi)**2
            do ii = 0,nz-2,2      !   note that we use double steps
            
            !---    RK4 algorithm step 1:
            !   compute derivative at z
            !   phi_0' = ( - acosgu phigi - asingu phigr ) + i ( acosgu phigr - asingu phigi )
            !   phi_g' = ( - acosgu phi0i + asingu phi0r - b phigi ) + i ( acosgu phi0r + asingu phi0i + b phigr)   
                bgc = bg*real(x(ii)) ; bgs = bg*aimag(x(ii))
                dphi0r =       b0 *phi0i - bgs*phigr + bgc*phigi
                dphi0i =      -b0 *phi0r + bgc*phigr - bgs*phigi 
                dphigr =  (a + b0)*phi0i + bgs*phigr + bgc*phigi 
                dphigi = -(a + b0)*phi0r + bgc*phigr + bgs*phigi 
                
            !   store temps and update to z+dz/2
                phi0r2 = phi0r
                phi0i2 = phi0i
                phigr2 = phigr
                phigi2 = phigi
                
                phi0r1 = phi0r + dphi0r * dz/6
                phi0i1 = phi0i + dphi0i * dz/6
                phigr1 = phigr + dphigr * dz/6
                phigi1 = phigi + dphigi * dz/6
               
                phi0r  = phi0r + dphi0r * dz/2
                phi0i  = phi0i + dphi0i * dz/2
                phigr  = phigr + dphigr * dz/2
                phigi  = phigi + dphigi * dz/2
               
                
            !---    RK4 algorithm step 2:
            !   compute derivative at z+dz/2
                !dphi0r = -acosgu(ii+1)*phigi - asingu(ii+1)*phigr           - b0*phi0i !- lambda*phi0r 
                !dphi0i = acosgu(ii+1)*phigr - asingu(ii+1)*phigi            + b0*phi0r !- lambda*phi0i
                !dphigr = -bg*phigi - acosgu(ii+1)*phi0i + asingu(ii+1)*phi0r - b0*phigi !- lambda*phigr
                !dphigi = bg*phigr + acosgu(ii+1)*phi0r + asingu(ii+1)*phi0i  + b0*phi0r !- lambda*phigi
                bgc = bg*real(x(ii+1)) ; bgs = bg*aimag(x(ii+1))
                dphi0r =       b0 *phi0i - bgs*phigr + bgc*phigi
                dphi0i =      -b0 *phi0r + bgc*phigr - bgs*phigi 
                dphigr =  (a + b0)*phi0i + bgs*phigr + bgc*phigi 
                dphigi = -(a + b0)*phi0r + bgc*phigr + bgs*phigi 
                
            !   store temps and recompute z+dz/2
                phi0r  = phi0r2 + dphi0r * dz/2
                phi0i  = phi0i2 + dphi0i * dz/2
                phigr  = phigr2 + dphigr * dz/2
                phigi  = phigi2 + dphigi * dz/2
                                   
                phi0r1 = phi0r1 + dphi0r * dz/3
                phi0i1 = phi0i1 + dphi0i * dz/3
                phigr1 = phigr1 + dphigr * dz/3
                phigi1 = phigi1 + dphigi * dz/3
               
                
            !---    RK4 algorithm step 3:
            !   recompute derivative at z+dz/2
!                 dphi0r = -acosgu(ii+1)*phigi - asingu(ii+1)*phigr           - b0*phi0i !- lambda*phi0r
!                 dphi0i = acosgu(ii+1)*phigr - asingu(ii+1)*phigi            + b0*phi0r !- lambda*phi0i
!                 dphigr = -bg*phigi - acosgu(ii+1)*phi0i + asingu(ii+1)*phi0r - b0*phigi !- lambda*phigr
!                 dphigi = bg*phigr + acosgu(ii+1)*phi0r + asingu(ii+1)*phi0i  + b0*phi0r !- lambda*phigi
                bgc = bg*real(x(ii+1)) ; bgs = bg*aimag(x(ii+1))
                dphi0r =       b0 *phi0i - bgs*phigr + bgc*phigi
                dphi0i =      -b0 *phi0r + bgc*phigr - bgs*phigi 
                dphigr =  (a + b0)*phi0i + bgs*phigr + bgc*phigi 
                dphigi = -(a + b0)*phi0r + bgc*phigr + bgs*phigi        
                         
            !   store temps and compute z+dz
                phi0r  = phi0r2 + dphi0r * dz
                phi0i  = phi0i2 + dphi0i * dz
                phigr  = phigr2 + dphigr * dz
                phigi  = phigi2 + dphigi * dz
                                   
                phi0r1 = phi0r1 + dphi0r * dz/3
                phi0i1 = phi0i1 + dphi0i * dz/3
                phigr1 = phigr1 + dphigr * dz/3
                phigi1 = phigi1 + dphigi * dz/3
               
            
            !---    RK4 algorithm step 4:
            !   recompute derivative at z+dz
!             dphi0r = -acosgu(ii+2)*phigi - asingu(ii+2)*phigr           - b0*phi0i !- lambda*phi0r
!             dphi0i = acosgu(ii+2)*phigr - asingu(ii+2)*phigi            + b0*phi0r !- lambda*phi0i
!             dphigr = -bg*phigi - acosgu(ii+2)*phi0i + asingu(ii+2)*phi0r - b0*phigi !- lambda*phigr
!             dphigi = bg*phigr + acosgu(ii+2)*phi0r + asingu(ii+2)*phi0i  + b0*phi0r !- lambda*phigi
                bgc = bg*real(x(ii+2)) ; bgs = bg*aimag(x(ii+2))
                dphi0r =       b0 *phi0i - bgs*phigr + bgc*phigi
                dphi0i =      -b0 *phi0r + bgc*phigr - bgs*phigi 
                dphigr =  (a + b0)*phi0i + bgs*phigr + bgc*phigi 
                dphigi = -(a + b0)*phi0r + bgc*phigr + bgs*phigi    
                             
            !   recompute z+dz
                phi0r  = phi0r1 + dphi0r * dz/6
                phi0i  = phi0i1 + dphi0i * dz/6
                phigr  = phigr1 + dphigr * dz/6
                phigi  = phigi1 + dphigi * dz/6
                
                
                !  print *,"column ",ii,phi0r,phi0i,phigr,phigi,(phi0r+phigr)**2 + (phi0i+phigi)**2
                             
                phi0r1 = (phi0r+phigr)
                phi0i1 = (phi0i+phigi)
                      
                phimag =  phi0r1*phi0r1 +  phi0i1*phi0i1
                phimag = 1/sqrt(phimag)
                
                phi0r = phi0r * phimag
                phi0i = phi0i * phimag
                phigr = phigr * phimag
                phigi = phigi * phimag
                
            end do           
            
             
            
        !---    finally compute the intensity by combining the two beams
           ! phi0r = phi0r + phigr
           ! phi0i = phi0i + phigi        
            phig2 = phigr*phigr + phigi*phigi
            
            return
        end subroutine dynamicalTwoBeamIntegrationColumn_phase
            
    
    
        subroutine dynamicalTwoBeamIntegrationColumn_disp( a,b,b0,g,u , phi2 )
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      integrate dynamical two beam equations through depth zmax
    !*      given the displacement field u(z) 
    !*      and the scattering vector g ( in reciprocal cell units )
    !*      and two parameters  a = pi / xi_g
    !*                          b = - eps / (h v)
    !*      with xi_g = pi h v / | U_g |
    !*      eps = hbar^2 ( k+g )^2 / (2m) - hbar^2 k^2/(2m)
    !*      v = hbar k / m
    !*
    !*      return the intensity squared phi2
    !*
    !*      phi_0' = ia exp[ i g.u ] phi_g
    !*      phi_g' = ib phi_g + ia exp[ -i g.u ] phi_0
    
    !*      phi_0' = i (acosgu + i asingu) (phigr + i phigi)
    !*             = ( - acosgu phigi - asingu phigr ) + i ( acosgu phigr - asingu phigi )
    !*      phi_g' = i (b phigr + i b phigi) + i (acosgu - i asingu) (phi0r + i phi0i)
    !*             = ( - acosgu phi0i + asingu phi0r - b phigi ) + i ( acosgu phi0r + asingu phi0i + b phigr)
    !*          
    
            real(kind=real64),intent(in)                    ::      a,b,b0
            real(kind=real64),dimension(3),intent(in)       ::      g
            real(kind=real64),dimension(:,0:),intent(in)    ::      u
            
            real(kind=real64),intent(out)                   ::      phi2
            
            real(kind=real64),dimension(0:size(u)-1)        ::  asingu,acosgu,bb,bb0
            
            
            integer             ::      nz
            real(kind=real64)   ::      dz,gdotu
            real(kind=real64)   ::      lambda = 0.01d0
            integer             ::      ii
            
            
            real(kind=real64)   ::      phi0r ,phi0i ,phigr ,phigi          !   these are the field evaluated at z ( real and imaginary parts )
            real(kind=real64)   ::      dphi0r,dphi0i,dphigr,dphigi         !   derivative of field evalulated at z,z+dz/2,z+dz
            real(kind=real64)   ::      phi0r1,phi0i1,phigr1,phigi1         !   temporary 1 needed for RK4 integration
            real(kind=real64)   ::      phi0r2,phi0i2,phigr2,phigi2         !   temporary 2 needed for RK4 integration
            
        !---    find the number of z divisions : note this should be an odd number ( checked in area() call )
            nz = size(u,dim=2)-1
            dz = 2.0d0 / nz      !   note we use double steps because the derivative is computed at z,z+dz/2,z+dz
            
                    
            
        !---    first precompute the exponential term 
        !*      this can be done quickly using a modern compiler's sincos function
            do ii = 0,nz
                gdotu = g(1)*u(1,ii) + g(2)*u(2,ii) + g(3)*u(3,ii)
                asingu(ii) = sin( gdotu )
                acosgu(ii) = cos( gdotu )                
            !   print *,ii,u(:,ii),a*asingu(ii),a*acosgu(ii)
            end do
            asingu = asingu*a
            acosgu = acosgu*a
            bb = b
            bb0 = b0
                                                                                                       
            asingu(0) = 0                    ; acosgu(0) = 0                    ; bb(0) = 0.0d0        ; bb0(0) = 0.0d0       
            asingu(1) = asingu(1) * 0.1035d0 ; acosgu(1) = acosgu(1) * 0.1035d0 ; bb(1) = b * 0.1035d0 ; bb0(1) = b0 * 0.1035d0
            asingu(2) = asingu(2) * 0.5000d0 ; acosgu(2) = acosgu(2) * 0.5000d0 ; bb(2) = b * 0.5000d0 ; bb0(2) = b0 * 0.5000d0
            asingu(3) = asingu(3) * 0.8965d0 ; acosgu(3) = acosgu(3) * 0.8965d0 ; bb(3) = b * 0.8965d0 ; bb0(3) = b0 * 0.8965d0
            
        !    next initialise the propagating fields - wlog set real component of phi0 to 1, everything else zero
            phi0r = 1.0d0 ; phi0i = 0.0d0  ; phigr = 0.0d0 ; phigi = 0.0d0
            
            
        !   now integrate in z        
            !print *,"column disp ",0,phi0r,phi0i,phigr,phigi,(phi0r+phigr)**2 + (phi0i+phigi)**2
            do ii = 0,nz-2,2      !   note that we use double steps
            
            !---    RK4 algorithm step 1:
            !   compute derivative at z
            !   phi_0' = ( - acosgu phigi - asingu phigr ) + i ( acosgu phigr - asingu phigi )
            !   phi_g' = ( - acosgu phi0i + asingu phi0r - b phigi ) + i ( acosgu phi0r + asingu phi0i + b phigr)
                dphi0r = -acosgu(ii)*phigi - asingu(ii)*phigr                - bb0(ii)*phi0i - lambda*phi0r
                dphi0i = acosgu(ii)*phigr - asingu(ii)*phigi                 + bb0(ii)*phi0r - lambda*phi0i
                dphigr = -bb(ii)*phigi - acosgu(ii)*phi0i + asingu(ii)*phi0r - bb0(ii)*phigi - lambda*phigr
                dphigi = bb(ii)*phigr + acosgu(ii)*phi0r + asingu(ii)*phi0i  + bb0(ii)*phi0r - lambda*phigi
                
            !   store temps and update to z+dz/2
                phi0r2 = phi0r
                phi0i2 = phi0i
                phigr2 = phigr
                phigi2 = phigi
                
                phi0r1 = phi0r + dphi0r * dz/6
                phi0i1 = phi0i + dphi0i * dz/6
                phigr1 = phigr + dphigr * dz/6
                phigi1 = phigi + dphigi * dz/6
               
                phi0r  = phi0r + dphi0r * dz/2
                phi0i  = phi0i + dphi0i * dz/2
                phigr  = phigr + dphigr * dz/2
                phigi  = phigi + dphigi * dz/2
               
                
            !---    RK4 algorithm step 2:
            !   compute derivative at z+dz/2
                dphi0r = -acosgu(ii+1)*phigi - asingu(ii+1)*phigr                  - bb0(ii+1)*phi0i - lambda*phi0r 
                dphi0i = acosgu(ii+1)*phigr - asingu(ii+1)*phigi                   + bb0(ii+1)*phi0r - lambda*phi0i
                dphigr = -bb(ii+1)*phigi - acosgu(ii+1)*phi0i + asingu(ii+1)*phi0r - bb0(ii+1)*phigi - lambda*phigr
                dphigi = bb(ii+1)*phigr + acosgu(ii+1)*phi0r + asingu(ii+1)*phi0i  + bb0(ii+1)*phi0r - lambda*phigi
                
            !   store temps and recompute z+dz/2
                phi0r  = phi0r2 + dphi0r * dz/2
                phi0i  = phi0i2 + dphi0i * dz/2
                phigr  = phigr2 + dphigr * dz/2
                phigi  = phigi2 + dphigi * dz/2
                                   
                phi0r1 = phi0r1 + dphi0r * dz/3
                phi0i1 = phi0i1 + dphi0i * dz/3
                phigr1 = phigr1 + dphigr * dz/3
                phigi1 = phigi1 + dphigi * dz/3
               
                
            !---    RK4 algorithm step 3:
            !   recompute derivative at z+dz/2
                dphi0r = -acosgu(ii+1)*phigi - asingu(ii+1)*phigr                  - bb0(ii+1)*phi0i - lambda*phi0r
                dphi0i = acosgu(ii+1)*phigr - asingu(ii+1)*phigi                   + bb0(ii+1)*phi0r - lambda*phi0i
                dphigr = -bb(ii+1)*phigi - acosgu(ii+1)*phi0i + asingu(ii+1)*phi0r - bb0(ii+1)*phigi - lambda*phigr
                dphigi = bb(ii+1)*phigr + acosgu(ii+1)*phi0r + asingu(ii+1)*phi0i  + bb0(ii+1)*phi0r - lambda*phigi
                
            !   store temps and compute z+dz
                phi0r  = phi0r2 + dphi0r * dz
                phi0i  = phi0i2 + dphi0i * dz
                phigr  = phigr2 + dphigr * dz
                phigi  = phigi2 + dphigi * dz
                                   
                phi0r1 = phi0r1 + dphi0r * dz/3
                phi0i1 = phi0i1 + dphi0i * dz/3
                phigr1 = phigr1 + dphigr * dz/3
                phigi1 = phigi1 + dphigi * dz/3
               
            
            !---    RK4 algorithm step 4:
            !   recompute derivative at z+dz
                dphi0r = -acosgu(ii+2)*phigi - asingu(ii+2)*phigr                - bb0(ii+2)*phi0i - lambda*phi0r
                dphi0i = acosgu(ii+2)*phigr - asingu(ii+2)*phigi                 + bb0(ii+2)*phi0r - lambda*phi0i
                dphigr = -bb(ii+2)*phigi - acosgu(ii+2)*phi0i + asingu(ii+2)*phi0r - bb0(ii+2)*phigi - lambda*phigr
                dphigi = bb(ii+2)*phigr + acosgu(ii+2)*phi0r + asingu(ii+2)*phi0i  + bb0(ii+2)*phi0r - lambda*phigi
                
            !   recompute z+dz
                phi0r  = phi0r1 + dphi0r * dz/6
                phi0i  = phi0i1 + dphi0i * dz/6
                phigr  = phigr1 + dphigr * dz/6
                phigi  = phigi1 + dphigi * dz/6
                
                
                !print *,"column disp ",ii,phi0r,phi0i,phigr,phigi,(phi0r+phigr)**2 + (phi0i+phigi)**2
                                   
            end do           
            
            
            !stop
            
        !---    finally compute the intensity by combining the two beams
            phi0r = phi0r + phigr
            phi0i = phi0i + phigi        
            phi2 = phi0r*phi0r + phi0i*phi0i
            
            return
        end subroutine dynamicalTwoBeamIntegrationColumn_disp
            
    
    
    
        subroutine dynamicalTwoBeamIntegrationColumn_strain( a,b,g,dudz , phi2 )
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      integrate dynamical two beam equations through depth zmax
    !*      given the displacement field u(z) 
    !*      and the scattering vector g ( in reciprocal cell units )
    !*      and two parameters  a = pi / xi_g
    !*                          b = - eps / (h v)
    !*      with xi_g = pi h v / | U_g |
    !*      eps = hbar^2 ( k+g )^2 / (2m) - hbar^2 k^2/(2m)
    !*      v = hbar k / m
    !*
    !*      return the intensity squared phi2
    !*
    !*      phi_0' = ia phi_g
    !*      phi_g' = i[ g. dudz + b ] phi_g + ia phi_0
    
    
            real(kind=real64),intent(in)                    ::      a,b
            real(kind=real64),dimension(3),intent(in)       ::      g
            real(kind=real64),dimension(:,0:),intent(in)    ::      dudz
            
            real(kind=real64),intent(out)                   ::      phi2
            
            real(kind=real64),dimension(0:size(dudz)-1)        ::      gdotdudzb
            
            
            integer             ::      nz
            real(kind=real64)   ::      dz
            integer             ::      ii
            
            
            real(kind=real64)   ::      phi0r ,phi0i ,phigr ,phigi          !   these are the field evaluated at z ( real and imaginary parts )
            real(kind=real64)   ::      dphi0r,dphi0i,dphigr,dphigi         !   derivative of field evalulated at z,z+dz/2,z+dz
            real(kind=real64)   ::      phi0r1,phi0i1,phigr1,phigi1         !   temporary 1 needed for RK4 integration
            real(kind=real64)   ::      phi0r2,phi0i2,phigr2,phigi2         !   temporary 2 needed for RK4 integration
            
        !---    find the number of z divisions : note this should be an odd number ( checked in area() call )
            nz = size(dudz,dim=2)-1
            dz = 2.0d0 / nz      !   note we use double steps because the derivative is computed at z,z+dz/2,z+dz
            
                    
            
            do ii = 0,nz
                gdotdudzb = g(1)*dudz(1,ii) + g(2)*dudz(2,ii) + g(3)*dudz(3,ii)
            end do
            gdotdudzb = gdotdudzb + b
            
            
        !    next initialise the propagating fields - wlog set real component of phi0 to 1, everything else zero
            phi0r = 1.0d0 ; phi0i = 0.0d0  ; phigr = 0.0d0 ; phigi = 0.0d0
            
            
        !   now integrate in z        
            do ii = 0,nz-2,2      !   note that we use double steps
            
            !---    RK4 algorithm step 1:
            !   compute derivative at z
                dphi0r = -a*phigi
                dphi0i = a*phigr
                dphigr = -gdotdudzb(ii)*phigi - a*phi0i 
                dphigi = gdotdudzb(ii)*phigr + a*phi0r 
                
            !   store temps and update to z+dz/2
                phi0r2 = phi0r
                phi0i2 = phi0i
                phigr2 = phigr
                phigi2 = phigi
                
                phi0r1 = phi0r + dphi0r * dz/6
                phi0i1 = phi0i + dphi0i * dz/6
                phigr1 = phigr + dphigr * dz/6
                phigi1 = phigi + dphigi * dz/6
               
                phi0r  = phi0r + dphi0r * dz/2
                phi0i  = phi0i + dphi0i * dz/2
                phigr  = phigr + dphigr * dz/2
                phigi  = phigi + dphigi * dz/2
               
                
            !---    RK4 algorithm step 2:
            !   compute derivative at z+dz/2
                dphi0r = -a*phigi
                dphi0i = a*phigr
                dphigr = -gdotdudzb(ii+1)*phigi - a*phi0i 
                dphigi = gdotdudzb(ii+1)*phigr + a*phi0r 
                
            !   store temps and recompute z+dz/2
                phi0r  = phi0r2 + dphi0r * dz/2
                phi0i  = phi0i2 + dphi0i * dz/2
                phigr  = phigr2 + dphigr * dz/2
                phigi  = phigi2 + dphigi * dz/2
                                   
                phi0r1 = phi0r1 + dphi0r * dz/3
                phi0i1 = phi0i1 + dphi0i * dz/3
                phigr1 = phigr1 + dphigr * dz/3
                phigi1 = phigi1 + dphigi * dz/3
               
                
            !---    RK4 algorithm step 3:
            !   recompute derivative at z+dz/2
                dphi0r = -a*phigi
                dphi0i = a*phigr
                dphigr = -gdotdudzb(ii+1)*phigi - a*phi0i 
                dphigi = gdotdudzb(ii+1)*phigr + a*phi0r 
                
            !   store temps and compute z+dz
                phi0r  = phi0r2 + dphi0r * dz
                phi0i  = phi0i2 + dphi0i * dz
                phigr  = phigr2 + dphigr * dz
                phigi  = phigi2 + dphigi * dz
                                   
                phi0r1 = phi0r1 + dphi0r * dz/3
                phi0i1 = phi0i1 + dphi0i * dz/3
                phigr1 = phigr1 + dphigr * dz/3
                phigi1 = phigi1 + dphigi * dz/3
               
            
            !---    RK4 algorithm step 4:
            !   recompute derivative at z+dz
                dphi0r = -a*phigi
                dphi0i = a*phigr
                dphigr = -gdotdudzb(ii+2)*phigi - a*phi0i 
                dphigi = gdotdudzb(ii+2)*phigr + a*phi0r 
                
            !   recompute z+dz
                phi0r  = phi0r1 + dphi0r * dz/6
                phi0i  = phi0i1 + dphi0i * dz/6
                phigr  = phigr1 + dphigr * dz/6
                phigi  = phigi1 + dphigi * dz/6
                                   
            end do           
            
            
            
            
        !---    finally compute the intensity by combining the two beams
            phi0r = phi0r + phigr
            phi0i = phi0i + phigi        
            phi2 = phi0r*phi0r + phi0i*phi0i
            
            return
        end subroutine dynamicalTwoBeamIntegrationColumn_strain
            
    
    end module Lib_DynamicalTwoBeam         