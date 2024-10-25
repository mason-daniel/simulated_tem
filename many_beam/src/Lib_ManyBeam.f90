
    module Lib_ManyBeam
!---^^^^^^^^^^^^^^^^^^^
!*      code to propagate many electron beams through a defected crystal
!*
!*          (k+g)/|k| . grad Phi = i M Phi 
!*      where
!*          M = 2 pi s_g - pi rho iXi - Y
!*      s_g is diagonal matrix holding deviation parameters
!*      rho is atom density scaled to 1 in crystal and 0 in vacuum
!*      iXi is matrix of inverse (complex) extinction distances
!*      Y is a diagonal matrix containing the strain contribution
!*          y_g = Im( x_g* ( (k+g)/|k| ) . grad x_g )
!*      
!*
!*      Daniel Mason
!*      (c) UKAEA Sept 2024
!*

        use iso_fortran_env
        use Lib_Gvectors
        use Lib_RelativisticElectrons

        implicit none
        private

        external        ::      DSYEV



        real(kind=real64),parameter                 ::      PI = 3.14159265390d0
        complex(kind=real64),parameter              ::      EYE = complex( 0.0d0,1.0d0 )

        logical,public                              ::      LIB_MANYBEAM_DBG = .false.

        public          ::      ManyBeam_ctor
        public          ::      report
        public          ::      delete

        public          ::      setPhaseFactors
        public          ::      sanityCheck
        public          ::      deviationParameter
        public          ::      finddPhidz
        public          ::      getAngle
        public          ::      setOrientationDependence
        public          ::      perfectLatticeIntensity



        type,public     ::      ManyBeam
            private
            type(Gvectors),pointer                              ::      gv
            integer                                             ::      nG                  !   number of g-vectors            
            real(kind=real64),dimension(3)                      ::      kvec                !   propagating beam is usually z-direction, ie kvec = (/ 0,0,k /), unless specifically tweaked for precession method
            real(kind=real64)                                   ::      k                   !   propagating beam is along z-direction with magnitude k
            real(kind=real64),dimension(:,:),pointer            ::      kplusgonk           !   (3,0:nG)      (k + g)/|k|
            real(kind=real64),dimension(:),pointer              ::      Sg                  !   (0:nG)      deviation parameter 
            complex(kind=real64),dimension(:,:),pointer         ::      iXi                 !   (0:nG,0:nG) matrix of inverse (complex) extinction distances
            real(kind=real64)                                   ::      a                   !   cell spacing, needed to find finite difference gradients with correct length scaling
            real(kind=real64),dimension(:,:,:,:,:),pointer      ::      grad_arg_x          !   (3,n+veG,lbx:ubx,lby:uby,lbz:ubz) pointer to chunk of phase factor gradients            
            real(kind=real64),dimension(:,:,:),pointer          ::      rho                 !   (lbx:ubx,lby:uby,lbz:ubz) pointer to chunk of atom densities
        end type 


        interface           ManyBeam_ctor
            module procedure        ManyBeam_null
            module procedure        ManyBeam_ctor0
        end interface

        interface           report
            module procedure        report0
        end interface

        interface           delete
            module procedure        delete0
        end interface


        interface           sanityCheck
            module procedure        sanityCheck0
        end interface


        interface           perfectLatticeIntensity
            module procedure        perfectLatticeIntensity0
            module procedure        perfectLatticeIntensity1
        end interface
        


    contains
!---^^^^^^^^

        function ManyBeam_null() result(this)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      empty constructor
            type(ManyBeam)          ::      this
            this%nG = 0
            this%k = 1.0d0
            this%a = 1.0d0
            nullify(this%kplusgonk)
            nullify(this%Sg)
            nullify(this%ixi)
            nullify(this%gv)
            nullify(this%grad_arg_x)
            nullify(this%rho)
            return
        end function ManyBeam_null


        function ManyBeam_ctor0(k,gv,gv2,xi) result(this)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      default constructor - sets up the extinction distances and deviation parameter
    !*      but does not link pointers to any atomistic information. See also setPhaseFactors() to complete construction.
    !*
            real(kind=real64),dimension(3),intent(in)       ::      k                   !   incident electron beam vector. Typically (0,0,|k|), unless specifically chosen otherwise for precession method
            type(Gvectors),intent(in),target                ::      gv                  !   set of g-vectors
            type(Gvectors),intent(in)                       ::      gv2                 !   set of g-vectors including all g" = g - g'
            complex(kind=real64),dimension(:),intent(in)    ::      xi                  !   (1:ng2) (complex) extinction distances. Must have same ordering as g-vectors in g2
            type(ManyBeam)                                  ::      this          
            
           ! real(kind=real64),dimension(3)  ::      gg
            integer         ::      ii,jj,kk

            this = ManyBeam_null()
            this%gv => gv
            this%nG = getn(this%gv)                     
            this%kvec = k
            this%k = norm2(k)                               !   only need to store modulus

        !---    set orientation dependence
            allocate(this%Sg(0:this%nG))
            allocate(this%kplusgonk(3,0:this%nG))
            call setOrientationDependence(this)

        !---    set extinction distances
            allocate(this%ixi(0:this%nG,0:this%nG))
            do jj = 0,this%nG                
                do ii = 0,this%nG

                !---    matrix element iXi(i,j) = 1/( xi_{g_i - g_j} )
                    kk = whichg( gv2, getG(this%gv,jj) - getG(this%gv,ii) )
                    if (kk == LIB_GVECTORS_UNSET) then
                        print *,"ManyBeam_ctor0 info - ",ii,getG(this%gv,ii),jj,getG(this%gv,jj)
                        stop "Lib_ManyBeam::ManyBeam_ctor0 Error - incomplete g-vector set. Can't find required xi_g"
                    end if
  
                    this%iXi(ii,jj) = 1/xi(kk)
                end do
            end do
 
            return
        end function ManyBeam_ctor0


        

        subroutine delete0(this)
    !---^^^^^^^^^^^^^^^^^^^^^^^^
            type(ManyBeam),intent(inout)       ::      this
            if (this%nG==0) return
            deallocate( this%kplusgonk )
            deallocate( this%Sg )
            deallocate( this%iXi )
            this = ManyBeam_null()
            return
        end subroutine delete0

        
        subroutine report0(this,u,o)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^
            type(ManyBeam),intent(in)           ::      this
            integer,intent(in),optional         ::      u,o
            integer     ::      uu,oo
            integer     ::      ii
            
            uu = 6 ; if (present(u)) uu = u
            oo = 0 ; if (present(o)) oo = o
            write(unit=uu,fmt='(a,i6,a,3f12.6,a)') repeat(" ",oo)//"ManyBeam [nG=",this%nG,", k = ",this%kvec," (1/A)]"
            if (isMillerBravais(this%gv)) then
                write(unit=uu,fmt='(a,a16,3a12,a36)') repeat(" ",oo+4),"Miller-Bravais","S_g (1/A)"," Re[xi_g]"," Im[xi_g]"," (k+g)/|k|"
                write(unit=uu,fmt='(a,4i4,6f12.6)') repeat(" ",oo+4),(/0,0,0,0/),this%Sg(0),real(1/this%iXi(0,0)),aimag(1/this%iXi(0,0)),this%kplusgonk(:,0)
                do ii = 1,this%nG
                    write(unit=uu,fmt='(a,4i4,6f12.6)') repeat(" ",oo+4),gethjkl(this%gv,ii),this%Sg(ii),real(1/this%iXi(ii,0)),aimag(1/this%iXi(ii,0)),this%kplusgonk(:,ii)
                end do
            else
                write(unit=uu,fmt='(a,a12,3a12,a36)') repeat(" ",oo+4),"Miller index","S_g (1/A)"," Re[xi_g]"," Im[xi_g]"," (k+g)/|k|"
                write(unit=uu,fmt='(a,3i4,6f12.6)') repeat(" ",oo+4),(/0,0,0/),this%Sg(0),real(1/this%iXi(0,0)),aimag(1/this%iXi(0,0)),this%kplusgonk(:,0)
                do ii = 1,this%nG
                    write(unit=uu,fmt='(a,3i4,6f12.6)') repeat(" ",oo+4),gethkl(this%gv,ii),this%Sg(ii),real(1/this%iXi(ii,0)),aimag(1/this%iXi(ii,0)),this%kplusgonk(:,ii)
                end do
            end if

            !call report(this%gv,uu,oo+4)
             
             
            return
        end subroutine report0

    !---

    
        subroutine setPhaseFactors( this,a, grad_arg_x,rho )
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      link pointers to the atomistic information, phase field and atomic density
            type(ManyBeam),intent(inout)                                    ::      this
            real(kind=real64),intent(in)                                    ::      a                   !   cell side length, needed for finite difference gradients
            !complex(kind=real64),dimension(:,:,:,:),pointer,intent(in)      ::      x                   !   (n+veG,lbx:ubx,lby:uby,lbz:ubz) pointer to chunk of phase factors
            real(kind=real64),dimension(:,:,:,:,:),pointer,intent(in)       ::      grad_arg_x              !   (3,n+veG,lbx:ubx,lby:uby,lbz:ubz) pointer to chunk of phase factor gradientss            
            real(kind=real64),dimension(:,:,:),pointer,intent(in)           ::      rho                 !   (lbx:ubx,lby:uby,lbz:ubz) pointer to chunk of atom densities

            this%a = a
            !this%x => x
            this%grad_arg_x => grad_arg_x
            this%rho => rho

            return
        end subroutine setPhaseFactors


        pure logical function sanityCheck0(this)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      returns true if all the pointers are set
            type(ManyBeam),intent(in)           ::      this            
            sanityCheck0 = (this%nG>0)

            !sanityCheck0 = sanityCheck0 .and. associated(this%x)
            sanityCheck0 = sanityCheck0 .and. associated(this%grad_arg_x)
            sanityCheck0 = sanityCheck0 .and. associated(this%rho)

            !sanityCheck0 = sanityCheck0 .and. (ubound(this%x,dim=1)==nPositiveg(this%gv)) .and. (ubound(this%grad_arg_x,dim=1)==ubound(this%x,dim=1))  
            sanityCheck0 = sanityCheck0 .and. (ubound(this%kplusgonk,dim=2)==this%nG) .and. (ubound(this%Sg,dim=1)==this%nG) .and. (ubound(this%iXi,dim=1)==this%nG)   
            

            return
        end function sanityCheck0


        pure real(kind=real64) function getAngle(this,i)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      return the absolute angle from the z-axis tan(theta) = max( |(k+g).(xy)| ) / (k+g).z for beam i
            type(ManyBeam),intent(in)           ::      this 
            integer,intent(in)                  ::      i
            real(kind=real64),dimension(3)      ::      xy

            xy(1:3) = (/ this%kplusgonk(1,i),this%kplusgonk(2,i),1.0d-32 /)
            xy = xy/norm2(xy)

            getAngle = atan2( this%kplusgonk(1,i)*xy(1) + this%kplusgonk(2,i)*xy(2) , this%kplusgonk(3,i) )
            
            return
        end function getAngle

!-------

         
        subroutine setOrientationDependence(this)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      A foil tilt has been applied to the AtomSpace and Gvectors
    !*      Change the orientation dependence in the beams        
            type(ManyBeam),intent(inout)                            ::      this
            integer                             ::      ii
            real(kind=real64),dimension(3)      ::      gg
            this%Sg(0) = 0.0d0
            this%kplusgonk(:,0) = (/ 0,0,1 /)
            do ii = 1,this%nG
                gg = getG(this%gv,ii)
                this%kplusgonk(:,ii) = (this%kvec(:) + gg(:)) / this%k
                this%Sg(ii) = deviationParameter(this%kvec(:),gg)                 
            end do
            return
        end subroutine setOrientationDependence



        subroutine perfectLatticeIntensity1(this, L_R, R, maxIntensity, Ig )
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      use the many-beam equations in the columnar/loss free approximation
    !*      to find the intensity at the back of the foil
    !*      assuming a tilt R.
            type(ManyBeam),intent(in)                               ::      this
            real(kind=real64),intent(in)                            ::      L_R     !   foil thickness after rotation
            real(kind=real64),dimension(3,3),intent(in)             ::      R
            logical,intent(in)                                      ::      maxIntensity
            real(kind=real64),dimension(0:this%nG),intent(out)      ::      Ig
 
            integer                             ::      ii
            real(kind=real64),dimension(3)      ::      gg
            real(kind=real64),dimension(0:this%nG)  ::      Sg

            
            Sg = 0
        !---    find deviation parameters when rotated to new position            
            do ii = 1,this%nG
                gg = getG(this%gv,ii)
                gg(:) = R(:,1)*gg(1) + R(:,2)*gg(2) + R(:,3)*gg(3)                 
                Sg(ii) = deviationParameter(this%kvec(:),gg)                 
            end do
    
            call perfectLatticeIntensity0(this, L_R, maxIntensity, Ig, Sg)


            return
        end subroutine perfectLatticeIntensity1

        subroutine perfectLatticeIntensity0(this, L, maxIntensity, Ig  , Sg)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      use the many-beam equations in the columnar/loss free approximation
    !*      to find the maximum intensity in the foil, 
    !*      or the intensity at the back of the foil
            type(ManyBeam),intent(in)                               ::      this
            real(kind=real64),intent(in)                            ::      L
            logical,intent(in)                                      ::      maxIntensity
            real(kind=real64),dimension(0:this%nG),intent(out)      ::      Ig
            real(kind=real64),dimension(0:this%nG),intent(in),optional      ::      Sg
            
            complex(kind=real64),dimension(0:this%nG,0:this%nG)     ::      Psi
        !    real(kind=real64),dimension(0:this%nG,0:this%nG)        ::      M_r
        !    real(kind=real64)                                       ::      rho
            !real(kind=real64),dimension(3,0:this%nG)                ::      grad_arg_x              !   (3,0:nG)

            real(kind=real64),dimension(0:this%nG)                  ::      lambda_g
            !real(kind=real64),dimension( 66*(1+this%nG) )           ::      work
            
            integer                                                 ::      ii,jj,kk            
            real(kind=real64)                   ::      zz
            complex(kind=real64)                ::      phi
 
            if (present(Sg)) then
                call perfectLatticeEigendecomp( this, lambda_g, Psi , Sg )
            else
                call perfectLatticeEigendecomp( this, lambda_g, Psi , Sg )
            end if

        ! !---    add the extinction distance term
        !     M_r(0:this%nG,0:this%nG) = - PI * real( this%iXi(0:this%nG,0:this%nG) )
            
        ! !---    add the Sg term to the diagonal. Note that Sg(0) = 0
        !     if (present(Sg)) then
        !         do ii = 1,this%nG               
        !             M_r(ii,ii) = M_r(ii,ii) + 2 * PI * Sg(ii)            
        !         end do
        !     else
        !         do ii = 1,this%nG               
        !             M_r(ii,ii) = M_r(ii,ii) + 2 * PI * this%Sg(ii)            
        !         end do
        !     end if
 

        !     if (LIB_MANYBEAM_DBG) then
        !         print *,"perfectLatticeIntensity0 info - propagationMatrix Re(M)"                 
        !         do ii = 0,this%nG
        !             write (*,fmt='(100f16.8)') M_r(ii,:)
        !         end do
                 
        ! !     end if

        ! !---    find eigenvalues and eigenvectors
        !     call DSYEV( "V","U",(1+this%nG),M_r,(1+this%nG),lambda,work,size(work),ii )

        !---    what is the thickness of the foil?            

            if (LIB_MANYBEAM_DBG) print *,"perfectLatticeIntensity0 info - lambda ",lambda_g
            
        !---    compute solution



            if (maxIntensity) then

                !   redefine Psi_gg' = Psi_gg' Psi_g'0^H
                do jj = 0,this%nG
                    do ii = 1,this%nG
                        Psi(ii,jj) =  Psi(ii,jj) * conjg( Psi(0,jj) )
                    end do
                    Psi(0,jj) =  Psi(0,jj) * conjg( Psi(0,jj) )
                end do
                
                !   compute max intensity at points z = 0,a,2a ... L 
                Ig = -huge(1.0)
                do kk = 0,floor( L/this%a )
                    zz = kk*this%a
                    do ii = 0,this%nG
                        phi = 0.0d0
                        do jj = 0,this%nG
                            phi = phi + Psi(ii,jj) * exp( EYE * lambda_g(jj) * zz ) 
                        end do
                        Ig(ii) = max( Ig(ii) , real(phi)*real(phi) + aimag(phi)*aimag(phi) )
                    end do
                end do

            else

                !   compute intensity at point z = L only
                do ii = 0,this%nG
                    phi = 0.0d0
                    do jj = 0,this%nG
                        phi = phi + Psi(ii,jj) * exp( EYE * lambda_g(jj) * L ) * conjg( Psi(0,jj) )
                    end do
                    Ig(ii) = real(phi)*real(phi) + aimag(phi)*aimag(phi)
                end do
    
            end if

            return
        end subroutine perfectLatticeIntensity0
  


        subroutine perfectLatticeEigendecomp( this, lambda_g, Psi , Sg )
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      construct the perfect lattice propagation matrix  M_c
    !*      and then perform the eigendecomposition  M_c = Psi Lambda Psi^T
    !*      return the eigenvalues lambda_g and the eigenvectors 

            type(ManyBeam),intent(in)                               ::      this
            real(kind=real64),dimension(0:this%nG),intent(out)      ::      lambda_g
            complex(kind=real64),dimension(0:this%nG,0:this%nG),intent(out)     ::      Psi
            real(kind=real64),dimension(0:this%nG),intent(in),optional          ::      Sg
            
            complex(kind=real64),dimension( 66*(1+this%nG) )        ::      work
            real(kind=real64),dimension( 3*(1+this%nG)-2 )          ::      rwork
            
            integer                                                 ::      ii           
 
        !---    add the extinction distance term
            Psi(0:this%nG,0:this%nG) = - PI * this%iXi(0:this%nG,0:this%nG)
            
        !---    add the Sg term to the diagonal. Note that Sg(0) = 0
            if (present(Sg)) then
                do ii = 1,this%nG               
                    Psi(ii,ii) = Psi(ii,ii) + 2 * PI * Sg(ii)            
                end do
            else
                do ii = 1,this%nG               
                    Psi(ii,ii) = Psi(ii,ii) + 2 * PI * this%Sg(ii)            
                end do
            end if
  

        !---    find eigenvalues and eigenvectors
            call ZHEEV( "V","U",(1+this%nG),Psi,(1+this%nG),lambda_g,work,size(work),rwork,ii )

            return
        end subroutine perfectLatticeEigendecomp



        pure function propagationMatrix( this, rho, grad_arg_x ) result( M )
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      (k+g)/|k| . grad Phi = i M Phi 
    !*      compute the matrix M at position r
    !*      where
    !*          M = 2 pi s_g - pi rho iXi - Y
    !*      s_g is diagonal matrix holding deviation parameters
    !*      rho is atom density scaled to 1 in crystal and 0 in vacuum
    !*      iXi is matrix of inverse (complex) extinction distances
    !*      Y is a diagonal matrix containing the strain contribution
    !*          y_g = Im( x_g* ( (k+g)/|k| ) . grad x_g )
            type(ManyBeam),intent(in)                           ::      this
            real(kind=real64),intent(in)                        ::      rho
            !complex(kind=real64),dimension(0:),intent(in)       ::      x                   !   (0:nG)
            real(kind=real64),dimension(:,:),intent(in)         ::      grad_arg_x              !   (3,1:nG)
            complex(kind=real64),dimension(0:this%nG,0:this%nG) ::      M                  

            integer                 ::      ii  
            complex(kind=real64)    ::      y_g         !   y_g = x_g* ( (k+g)/|k| ) . grad x_g

        !---    add the extinction distance term
            M(0:this%nG,0:this%nG) = - PI * rho * this%iXi(0:this%nG,0:this%nG)
            
        !---    add the Sg term to the diagonal. Note that Sg(0) = 0
            do ii = 1,this%nG               
                M(ii,ii) = M(ii,ii) + 2 * PI * this%Sg(ii)            
            end do

        !---    add the elastic strain term
            do ii = 1,this%nG
                y_g = this%kplusgonk(1,ii)*grad_arg_x(1,ii) + this%kplusgonk(2,ii)*grad_arg_x(2,ii) + this%kplusgonk(3,ii)*grad_arg_x(3,ii)
                M(ii,ii) =  M(ii,ii) - y_g
            end do
 

            return
        end function propagationMatrix

        pure function getGradPhi( phi, ix,iy , a ) result( gradphi )
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      compute the 2d gradient in the x-y- plane of the complex phase field x at node point (ix,iy,iz) with node spacing a
    !*      note that I will have phi on bounds lbx:ubx, but requesting gradient in interior (lbx+1:ubx-1)
            complex(kind=real64),dimension(:,:,:),pointer,intent(in)            ::      phi                 !   (0:nG,lbx:ubx,lby:uby)      a 2d slice of phi
            integer,intent(in)                                                  ::      ix,iy
            real(kind=real64),intent(in)                                        ::      a
            complex(kind=real64),dimension(2,0:size(phi,dim=1)-1)               ::      gradphi
            integer             ::      nG
            integer             ::      ii

        !---    find the size of the problem                 
            nG = ubound(phi,dim=1)
 
        !---    6 point kernel             
            do ii = 0,nG
                gradphi(1,ii) = - phi(ii,ix-1,iy-1) + phi(ii,ix+1,iy-1) - 2*phi(ii,ix-1,iy) + 2*phi(ii,ix+1,iy) - phi(ii,ix-1,iy+1) + phi(ii,ix+1,iy+1)
                gradphi(2,ii) = - phi(ii,ix-1,iy-1) - 2*phi(ii,ix,iy-1) - phi(ii,ix+1,iy-1) + phi(ii,ix-1,iy+1) + 2*phi(ii,ix,iy+1) + phi(ii,ix+1,iy+1)
            end do
            
            gradphi = gradphi / (8*a)

            return
        end function getGradPhi


        subroutine finddPhidz( this , phi, iz, dphidz, columnar )
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ 
    !*      compute derivative of phi in the z direction.
    !*      note that I will have phi on bounds lbx:ubx, but am requesting gradient in interior (lbx+1:ubx-1)
            type(ManyBeam),intent(in)                                       ::      this
            complex(kind=real64),dimension(:,:,:),pointer,intent(in)        ::      phi                     !   (0:nG,lbx:ubx,lby:uby)
            integer,intent(in)                                              ::      iz                      !   slice- needed to find corresponding x_g,rho etc
            complex(kind=real64),dimension(:,:,:),pointer,intent(out)       ::      dphidz                  !   (0:nG,lbx:ubx,lby:uby)
            logical,intent(in)                                              ::      columnar

            integer                                             ::      ix,iy
            complex(kind=real64),dimension(2,0:this%nG)         ::      gradphi_r
            !complex(kind=real64),dimension(0:this%nG)           ::      x_r               
            real(kind=real64),dimension(3,this%nG)              ::      grad_arg_x_r          
            real(kind=real64)                                   ::      rho_r
            complex(kind=real64),dimension(0:this%nG,0:this%nG) ::      M_r                
            complex(kind=real64)                                ::      dphi_g_r
            integer                                             ::      lbx,ubx,lby,uby         !   bounds of phi

            integer                                             ::      ii,jj,kk
            integer,dimension(this%nG)          ::          indx

            !print *,"integrate0 finddPhidz "!,size(phi,dim=1),size(dphidz,dim=1)

        !---    find the bounds of the set of nodes I am computing
            lbx = lbound(phi,dim=2) ; ubx = ubound(phi,dim=2)
            lby = lbound(phi,dim=3) ; uby = ubound(phi,dim=3)

            !  print *,"integrate0 finddPhidz phi    ",lbx,ubx,lby,uby," at ",iz
            !  print *,"integrate0 finddPhidz dphidz ",lbound(dphidz,dim=2),ubound(dphidz,dim=2),lbound(dphidz,dim=3),ubound(dphidz,dim=3)
            ! print *,"finddPhidz x      ",lbound(this%x,dim=2),ubound(this%x,dim=2),lbound(this%x,dim=3),ubound(this%x,dim=3),lbound(this%x,dim=4),ubound(this%x,dim=4)
            ! print *,"finddPhidz grad_arg_x ",lbound(this%grad_arg_x,dim=3),ubound(this%grad_arg_x,dim=3),lbound(this%grad_arg_x,dim=4),ubound(this%grad_arg_x,dim=4),lbound(this%grad_arg_x,dim=5),ubound(this%grad_arg_x,dim=5)
            ! print *,"finddPhidz rho    ",lbound(this%rho,dim=1),ubound(this%rho,dim=1),lbound(this%rho,dim=2),ubound(this%rho,dim=2),lbound(this%rho,dim=3),ubound(this%rho,dim=3)

        !---    only the "+ve" g-vectors are stored in this%x
        !       so find the index of each g-vector to the positive set...
            jj = 0
            indx = 0
            do ii = 1,this%nG
                if (isPositiveg(this%gv,ii)) then
                    jj = jj + 1
                    indx(ii) = jj
                end if
            end do
        !   ... and then find the index of g-vectors in the negative set. 
            do ii = 1,this%nG
                if (.not. isPositiveg(this%gv,ii)) then
                    do kk = 1,this%nG
                        if (indx(getMinusg(this%gv,ii)) == kk) then
                            indx(ii) = -kk
                            exit
                        end if
                    end do
                end if
            end do
           ! print *,"finddPhidz ",iz,columnar

            !x_r(0) = 1.0d0
            !grad_arg_x_r(:,0) = 0.0d0
            gradphi_r = 0.0d0
            do iy = lby+1,uby-1
                do ix = lbx+1,ubx-1


                    !LIB_MANYBEAM_DBG = ((abs(ix-7.5)<1).and.(abs(iy-7.5)<1)) 


                !---    extract the crystal properties at node ix,iy,iz
                    do ii = 1,this%nG
                        jj = indx(ii)
                        if (jj>0) then
                            grad_arg_x_r(1:3,ii) = this%grad_arg_x(1:3,jj , ix,iy,iz )                            
                        else    !   I need to extract x_g = x_{-g}*
                            grad_arg_x_r(1:3,ii) = - this%grad_arg_x(1:3,-jj , ix,iy,iz )
                        end if
                    end do
                    rho_r = this%rho(ix,iy,iz)                                     
                    !print *,"ix,iy,x_r,rho_r ",ix,iy,x_r,rho_r


                    if (.not. columnar) gradphi_r = getGradPhi( phi, ix,iy , this%a ) 
                    

                !---    compute the propagation matrix
                    
                    M_r = propagationMatrix( this, rho_r, grad_arg_x_r ) 
                    

                    if (LIB_MANYBEAM_DBG) then
                        print *,"phi ",ix,iy, phi(:,ix,iy)
                        !write (*,fmt='(100f16.8)') real(phi(:,ix,iy)),aimag(phi(:,ix,iy))
                        print *,"propagationMatrix Re(M)",ix,iy                        
                        do kk = 0,this%nG
                            write (*,fmt='(100f16.8)') real(M_r(kk,:))
                        end do

                        ! print *,"propagationMatrix Im(M)",ix,iy                        
                        ! do kk = 0,this%nG
                        !     write (*,fmt='(100f16.8)') aimag(M_r(kk,:))
                        ! end do
                    end if

                !---    compute the derivative in the z-direction
                !       (k+g)/|k| . grad Phi = i M Phi                                         
                    do kk = 0,this%nG
                        dphi_g_r = EYE * dot_product( M_r(kk,:),phi(:,ix,iy) )
                        dphi_g_r = dphi_g_r - this%kplusgonk(1,kk)*gradphi_r(1,kk) - this%kplusgonk(2,kk)*gradphi_r(2,kk)
                        dphidz(kk,ix,iy) = dphi_g_r / this%kplusgonk(3,kk) 
                    end do

                    if (LIB_MANYBEAM_DBG) then
                        print *,"dphidz ",ix,iy , dphidz(:,ix,iy)
                    end if


                end do
            end do

            
            return
        end subroutine finddPhidz

 
!-------

        pure real(kind=real64) function deviationParameter(k,g)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      given a g-vector and the magnitude of k,
    !*      compute the deviation paramter s_g
            real(kind=real64),dimension(3),intent(in)   ::      k
            real(kind=real64),dimension(3),intent(in)   ::      g
            real(kind=real64)           ::      eps_g,hbarv

            eps_g = ( dot_product( g,2*k+g ) ) * HBAR*HBAR/(2*ME)
            hbarv = HBAR*HBAR*norm2(k)/ME
            deviationParameter = - eps_g/(2*PI*hbarv)

            return
        end function deviationParameter

!-------         
                       

    end module Lib_ManyBeam