
    module Lib_ManyBeam
!---^^^^^^^^^^^^^^^^^^^
!*      code to propagate many electron beams through a defected crystal
!*
!*          (k+g)/|k| . grad Phi = M Phi 
!*      where
!*          M = i 2 pi s_g - i pi rho iXi - Y
!*      s_g is diagonal matrix holding deviation parameters
!*      rho is atom density scaled to 1 in crystal and 0 in vacuum
!*      iXi is matrix of inverse (complex) extinction distances
!*      Y is a diagonal matrix containing the strain contribution
!*          y_g = x_g* ( (k+g)/|k| ) . grad x_g
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



        real(kind=real64),parameter                 ::      PI = 3.14159265390d0
        complex(kind=real64),parameter              ::      EYE = complex( 0.0d0,1.0d0 )

 

        public          ::      ManyBeam_ctor
        public          ::      report
        public          ::      delete

        public          ::      setPhaseFactors
        public          ::      sanityCheck

        public          ::      finddPhidz
        public          ::      getAngle




        type,public     ::      ManyBeam
            private
            integer                                             ::      nG                  !   number of g-vectors
            real(kind=real64)                                   ::      k                   !   propagating beam is along z-direction with magnitude k
            real(kind=real64),dimension(:,:),pointer            ::      kplusgonk           !   (3,0:nG)      (k + g)/|k|
            real(kind=real64),dimension(:),pointer              ::      Sg                  !   (0:nG)      deviation parameter 
            complex(kind=real64),dimension(:,:),pointer         ::      iXi                 !   (0:nG,0:nG) matrix of inverse (complex) extinction distances

            real(kind=real64)                                   ::      a                   !   cell spacing, needed to find finite difference gradients with correct length scaling
            complex(kind=real64),dimension(:,:,:,:),pointer     ::      x                   !   (nG,lbx:ubx,lby:uby,lbz:ubz) pointer to chunk of phase factors
            complex(kind=real64),dimension(:,:,:,:,:),pointer   ::      grad_x              !   (3,nG,lbx:ubx,lby:uby,lbz:ubz) pointer to chunk of phase factor gradientss            
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
        


    contains
!---^^^^^^^^

        function ManyBeam_null() result(this)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      empty constructor
            type(ManyBeam)          ::      this
            this%nG = 0
            this%k = 1
            nullify(this%kplusgonk)
            nullify(this%Sg)
            nullify(this%iXi)
            this%a = 1.0d0
            nullify(this%x)
            nullify(this%grad_x)
            nullify(this%rho)
            return
        end function ManyBeam_null


        function ManyBeam_ctor0(k,g,xi) result(this)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      default constructor - sets up the extinction distances and deviation parameter
    !*      but does not link pointers to any atomistic information. See also setPhaseFactors() to complete construction.
    !*
            real(kind=real64),dimension(3),intent(in)       ::      k                   !   incident electron beam vector. Typically (0,0,|k|), unless specifically chosen otherwise for precession method
            type(Gvectors),intent(in)                       ::      g                   !   set of g-vectors
            complex(kind=real64),dimension(0:),intent(in)   ::      xi                  !   (complex) extinction distances. Must have same ordering as g-vectors 
            type(ManyBeam)                                  ::      this          
            
            real(kind=real64),dimension(3)  ::      gg
            integer         ::      ii,jj,kk

            this = ManyBeam_null()
            this%nG = getn(g)
            this%k = norm2(k)                               !   only need to store modulus
            allocate(this%kplusgonk(3,0:this%nG))
            allocate(this%Sg(0:this%nG))
            allocate(this%iXi(0:this%nG,0:this%nG))
            this%Sg(0) = 0.0d0
            do ii = 1,this%nG
                gg = getG(g,ii)
                this%kplusgonk(:,ii) = (k(:) + gg(:)) / this%k
                this%Sg(ii) = deviationParameter(k,gg)
            end do
            do jj = 0,this%nG                
                do ii = 0,this%nG
                    kk = whichg( g, getG(g,ii) - getG(g,jj) )
                    this%iXi(ii,jj) = xi(kk)
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
            type(ManyBeam),intent(in)          ::      this
            integer,intent(in),optional         ::      u,o
            integer     ::      uu,oo
            uu = 6 ; if (present(u)) uu = u
            oo = 0 ; if (present(o)) oo = o
            write(unit=uu,fmt='(a,i6,a)') repeat(" ",oo)//"ManyBeam [nG=",this%nG,"]"
            return
        end subroutine report0

    !---

        subroutine setPhaseFactors( this,a, x,grad_x,rho )
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      link pointers to the atomistic information, phase field and atomic density
            type(ManyBeam),intent(inout)                                    ::      this
            real(kind=real64),intent(in)                                    ::      a                   !   cell side length, needed for finite difference gradients
            complex(kind=real64),dimension(:,:,:,:),pointer,intent(in)      ::      x                   !   (nG,lbx:ubx,lby:uby,lbz:ubz) pointer to chunk of phase factors
            complex(kind=real64),dimension(:,:,:,:,:),pointer,intent(in)    ::      grad_x              !   (3,nG,lbx:ubx,lby:uby,lbz:ubz) pointer to chunk of phase factor gradientss            
            real(kind=real64),dimension(:,:,:),pointer,intent(in)           ::      rho                 !   (lbx:ubx,lby:uby,lbz:ubz) pointer to chunk of atom densities

            this%a = a
            this%x => x
            this%grad_x => grad_x
            this%rho => rho

            return
        end subroutine setPhaseFactors


        pure logical function sanityCheck0(this)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      returns true if all the pointers are set
            type(ManyBeam),intent(in)           ::      this            
            sanityCheck0 = (this%nG>0)

            sanityCheck0 = sanityCheck0 .and. associated(this%x)
            sanityCheck0 = sanityCheck0 .and. associated(this%grad_x)
            sanityCheck0 = sanityCheck0 .and. associated(this%rho)

            sanityCheck0 = sanityCheck0 .and. (ubound(this%x,dim=1)==this%nG) .and. (ubound(this%grad_x,dim=1)==this%nG)  
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

            getAngle = atan2( this%kplusgonk(3,i) , this%kplusgonk(1,i)*xy(1) + this%kplusgonk(2,i)*xy(2) )

            return
        end function getAngle

!-------

        pure function propagationMatrix( this, rho, x, grad_x ) result( M )
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      (k+g)/|k| . grad Phi = M Phi 
    !*      compute the matrix M at position r
    !*      where
    !*          M = i 2 pi s_g - i pi rho iXi - Y
    !*      s_g is diagonal matrix holding deviation parameters
    !*      rho is atom density scaled to 1 in crystal and 0 in vacuum
    !*      iXi is matrix of inverse (complex) extinction distances
    !*      Y is a diagonal matrix containing the strain contribution
    !*          y_g = x_g* ( (k+g)/|k| ) . grad x_g            
            type(ManyBeam),intent(in)                           ::      this
            real(kind=real64),intent(in)                        ::      rho
            complex(kind=real64),dimension(:),intent(in)        ::      x                   !   (1:nG)
            complex(kind=real64),dimension(:,:),intent(in)      ::      grad_x              !   (3,1:nG)
            complex(kind=real64),dimension(0:this%nG,0:this%nG) ::      M                  

            integer                 ::      ii
            complex(kind=real64)    ::      y_g         !   y_g = x_g* ( (k+g)/|k| ) . grad x_g

            M(0:this%nG,0:this%nG) = - EYE * PI * rho * this%iXi(0:this%nG,0:this%nG)

            do ii = 1,this%nG
                M(ii,ii) =  M(ii,ii) + 2 * EYE * PI * this%Sg(ii)            
            end do

            do ii = 1,this%nG
                y_g = this%kplusgonk(1,ii)*grad_x(1,ii) + this%kplusgonk(2,ii)*grad_x(2,ii) + this%kplusgonk(3,ii)*grad_x(3,ii)
                y_g = conjg( x(ii) ) * y_g 
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


        subroutine finddPhidz( this , phi, iz, dphidz )
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      compute derivative of phi in the z direction.
    !*      note that I will have phi on bounds lbx:ubx, but am requesting gradient in interior (lbx+1:ubx-1)
            type(ManyBeam),intent(in)                                       ::      this
            complex(kind=real64),dimension(:,:,:),pointer,intent(in)        ::      phi                     !   (0:nG,lbx:ubx,lby:uby)
            integer,intent(in)                                              ::      iz                      !   slice- needed to find corresponding x_g,rho etc
            complex(kind=real64),dimension(:,:,:),pointer,intent(out)       ::      dphidz                  !   (0:nG,lbx:ubx,lby:uby,lbz:ubz)


            integer                                             ::      ix,iy
            complex(kind=real64),dimension(2,0:this%nG)         ::      gradphi_r
            complex(kind=real64),dimension(this%nG)             ::      x_r               
            complex(kind=real64),dimension(3,this%nG)           ::      grad_x_r          
            real(kind=real64)                                   ::      rho_r
            complex(kind=real64),dimension(0:this%nG,0:this%nG) ::      M_r                
            complex(kind=real64)                                ::      dphi_g_r
            integer                                             ::      lbx,ubx,lby,uby         !   bounds of phi

            integer                                             ::      kk

        !---    find the bounds of the set of nodes I am computing
            lbx = lbound(phi,dim=2) ; ubx = ubound(phi,dim=2)
            lby = lbound(phi,dim=3) ; uby = ubound(phi,dim=3)

            do iy = lby+1,uby-1
                do ix = lbx+1,ubx-1

                !---    extract the crystal properties at node ix,iy,iz
                    x_r(1:this%nG) = this%x(1:this%nG , ix,iy,iz )
                    grad_x_r(1:3,1:this%nG) = this%grad_x(1:3,1:this%nG , ix,iy,iz )
                    rho_r = this%rho(ix,iy,iz)                    
                    gradphi_r = getGradPhi( phi, ix,iy , this%a ) 
                    
                !---    compute the propagation matrix
                    M_r = propagationMatrix( this, rho_r, x_r, grad_x_r ) 

                !---    compute the derivative in the z-direction
                !       (k+g)/|k| . grad Phi = M Phi                                         
                    do kk = 0,this%nG
                        dphi_g_r = dot_product( M_r(kk,:),phi(:,ix,iy) )
                        dphi_g_r = dphi_g_r - this%kplusgonk(1,kk)*gradphi_r(1,kk) - this%kplusgonk(2,kk)*gradphi_r(2,kk)
                        dphidz(kk,ix,iy) = dphi_g_r / this%kplusgonk(3,kk) 
                    end do


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