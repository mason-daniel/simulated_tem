
    module Lib_FibonacciSphere
!---^^^^^^^^^^^^^^^^^^^^^^^^^^^
!*      A simple implementation of the Fibonacci sphere algorithm
!*      to distribute points evenly on a sphere

        use iso_fortran_env
        implicit none
        private


        real(kind=real64),private,parameter             ::      PI = 3.141592653589790d0
        real(kind=real64),private,parameter             ::      PHI = PI * ( sqrt(5.0d0) - 1.0d0 )      !   golden angle in radians
        

        public      ::          fibonacci_sphere
        public      ::          fibonacci_cap


    contains
!--^^^^^^^^^


        subroutine fibonacci_sphere(x)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      place points on the unit sphere
            real(kind=real64),dimension(:,:),intent(out)        ::      x           !   (3,n)

            integer             ::      nn
            integer             ::      ii
            real(kind=real64)   ::      xx,yy,zz
            real(kind=real64)   ::      rho,theta

            nn = size(x,dim=2)

            do ii = 1,nn

                zz = 1 - 2*((ii-1)/real(nn-1))          !   z from 1 to -1

                rho = sqrt( max(0.0d0,1-zz*zz) )        !   radius of disc at z

                theta = PHI * (ii-1)

                xx = rho*cos(theta) 
                yy = rho*sin(theta)

                x(:,ii) = (/ xx,yy,zz /)
            end do

            return
        end subroutine fibonacci_sphere


        subroutine fibonacci_cap(x,psi)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      place points on the unit sphere with z > cos(psi) 
            real(kind=real64),dimension(:,:),intent(out)        ::      x           !   (3,n)
            real(kind=real64),intent(in)                        ::      psi         !   0:Pi
            integer             ::      nn
            integer             ::      ii
            real(kind=real64)   ::      xx,yy,zz
            real(kind=real64)   ::      rho,theta
            real(kind=real64)   ::      zmin,dz

            nn = size(x,dim=2)
            zmin = cos(psi)
            dz = (1.0d0 - zmin)/(nn-1)
            
            do ii = 1,nn

                zz = zmin + (ii-1)*dz

                rho = sqrt( max(0.0d0,1.0d0-zz*zz) )        !   radius of disc at z

                theta = PHI * (ii-1)

                xx = rho*cos(theta) 
                yy = rho*sin(theta)

                x(:,ii) = (/ xx,yy,zz /)
            end do
            


            return
        end subroutine fibonacci_cap

    end module Lib_FibonacciSphere





 