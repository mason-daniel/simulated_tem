
    module Lib_RK4
!---^^^^^^^^^^^^^^
!*      implementation of RK4 method.
!*      currently 2d complex fields phi(w,z) integrated in z direction
!*      note that this module is serial. If you need parallel integration, it is up to you to share phi across processes.
!*      
!*      Daniel Mason
!*      (c) UKAEA Oct 2024
!*
!*      if
!*          phi' = f(phi)
!*      then
!*          phi(z+dz) = phi(z) + dz/6( k1 + 2 k2 + 2 k3 + k4 )
!*      with
!*          k1 = f( phi(z) )
!*          k2 = f( phi(z) + (dz/2) k1 )
!*          k3 = f( phi(z) + (dz/2) k2 )
!*          k4 = f( phi(z) + dz k3 )
!*
!*      usage
!*          integrator = rk4_ctor(n,dz)
!*          call getphip_dphip(integrator , phip,dphip)
!*          [ set phip = phi(z=0) ]
!*          do ii = 0,nz-1,2
!*              call update(integrator)             !   initialises values 
!*              dz = 0
!*              do step = 1,4
!*                  [ compute dphip = f(phip,ii+dz) ]       !   dz = 0,1,2
!*                  call update(integrator,step,dz)
!*              end do
!*          end do
!*          [ extract phi(z) = phip ]
!*          call delete(integrator)


        use iso_fortran_env
        implicit none
        private

        public      ::      RK4_ctor         
        public      ::      report           
        public      ::      delete           

        public      ::      getphip_dphip
        public      ::      update

        type,public     ::      RK4
            private
            integer                 ::      n      !   number of nodes to integrate
            real(kind=real64)                               ::      deltaz
            complex(kind=real64),dimension(:),pointer       ::      phi0,phip,phiz      !   (1:n) phi0 is a copy of the original input state. phip is a pointer to original input state ie overwrites it.
            complex(kind=real64),dimension(:),pointer       ::      dphip               !   f(phi')
        end type


        interface   RK4_ctor
            module procedure        RK4_null
            module procedure        RK4_ctor0
        end interface

        interface   report
            module procedure        report0
        end interface
        
        interface   delete
            module procedure        delete0
        end interface
        
        interface   update
            module procedure        update0
            module procedure        update1
        end interface
    contains
!---^^^^^^^^^

        function RK4_null() result(this)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      empty constructor
            type(RK4)           ::      this
            this%n = 0
            this%deltaz = 0
            nullify(this%phi0)
            nullify(this%phip)
            nullify(this%phiz)
            nullify(this%dphip)
            return
        end function RK4_null

        function RK4_ctor0(n,deltaz) result(this)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      default constructor
            integer,intent(in)                  ::      n
            real(kind=real64),intent(in)        ::      deltaz
            type(RK4)           ::      this
            this%n = n
            this%deltaz = deltaz
                        
            allocate(this%phi0(this%n))
            allocate(this%phip(this%n))
            allocate(this%phiz(this%n)) 
            allocate(this%dphip(this%n))
            return
        end function RK4_ctor0
        
        subroutine report0(this,u,o)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^
            type(RK4),intent(in)          ::      this
            integer,intent(in),optional         ::      u,o
            integer     ::      uu,oo
            uu = 6 ; if (present(u)) uu = u
            oo = 0 ; if (present(o)) oo = o
            write(unit=uu,fmt='(a,i8,a)') repeat(" ",oo)//"RK4 [n=",this%n,"]"
            return
        end subroutine report0

        subroutine delete0(this)
    !---^^^^^^^^^^^^^^^^^^^^^^^^
            type(RK4),intent(inout)       ::      this
            if (this%n==0) return
            deallocate(this%phi0)
            deallocate(this%phip)
            deallocate(this%phiz)
            deallocate(this%dphip)
            this = RK4_null()
            return
        end subroutine delete0
        
    !---

        subroutine getphip_dphip( this,phip,dphip )
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      return pointers to phi' and dphi' which can be used to generate state and derivative at z
            type(RK4),intent(inout)                                 ::      this
            complex(kind=real64),dimension(:),pointer,intent(out)   ::      phip
            complex(kind=real64),dimension(:),pointer,intent(out)   ::      dphip
            phip => this%phip
            dphip => this%dphip
            return
        end subroutine getphip_dphip

        subroutine update0( this )
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^ 
    !*      Have computed the value of f(phi') externally, so add it to the RK4 accumulator
    !*
    !*          phi' = f(phi)
    !*      then
    !*          phi(z+dz) = phi(z) + dz/6( k1 + 2 k2 + 2 k3 + k4 )
    !*      with
    !*          k1 = f( phi(z) )
    !*          k2 = f( phi(z) + (dz/2) k1 )
    !*          k3 = f( phi(z) + (dz/2) k2 )
    !*          k4 = f( phi(z) + dz k3 )

            type(RK4),intent(inout)                             ::      this
           
            this%phi0(:) = this%phip(:)
             
            return
        end subroutine update0


        subroutine update1( this,step,dz )
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      Have computed the value of f(phi') externally, so add it to the RK4 accumulator
    !*
    !*          phi' = f(phi)
    !*      then
    !*          phi(z+dz) = phi(z) + dz/6( k1 + 2 k2 + 2 k3 + k4 )
    !*      with
    !*          k1 = f( phi(z) )
    !*          k2 = f( phi(z) + (dz/2) k1 )
    !*          k3 = f( phi(z) + (dz/2) k2 )
    !*          k4 = f( phi(z) + dz k3 )

            type(RK4),intent(inout)                             ::      this
            integer,intent(in)                                  ::      step
            integer,intent(out)                                 ::      dz

            dz = 0

            select case(step)
                case(1)
                    !   have initialised, and have dphip = k1 
                    this%phiz(:) = this%phi0(:) + (this%deltaz/6)*this%dphip(:)
                    this%phip(:) = this%phi0(:) + (this%deltaz/2)*this%dphip(:)
                    dz = 1
                case(2)
                    !   have dphip = k2  
                    this%phiz(:) = this%phiz(:) + (this%deltaz/3)*this%dphip(:)
                    this%phip(:) = this%phi0(:) + (this%deltaz/2)*this%dphip(:)
                    dz = 1
                case(3)
                    !   have dphip = k3
                    this%phiz(:) = this%phiz(:) + (this%deltaz/3)*this%dphip(:)
                    this%phip(:) = this%phi0(:) + this%deltaz*this%dphip(:)
                    dz = 2
                case(4)
                    !   have dphip = k4. Have completed the step
                    this%phip(:) = this%phiz(:) + (this%deltaz/6)*this%dphip(:)                    
                    dz = 0
            end select 
            return
        end subroutine update1

    end module Lib_RK4