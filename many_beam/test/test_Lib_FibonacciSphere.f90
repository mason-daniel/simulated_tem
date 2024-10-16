
    program test_Lib_FibonacciSphere
!---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
!*      simple test program
        use iso_fortran_env
        use Lib_FibonacciSphere
        implicit none

        real(kind=real64),parameter             ::      PI = 3.141592653589790d0
        integer             ::          n = 500                     !   number of points
        real(kind=real64)   ::          theta = PI                   !   max angle
        character(len=256)  ::      dummy

        real(kind=real64),dimension(:,:),allocatable        ::      x
        real(kind=real64),dimension(3)              ::      dx
        integer             ::      ii,jj
        real(kind=real64)   ::      dmin
        logical             ::      ok

        print *,"usage: test_Lib_FibonacciSphere [n [theta]]"
        call get_command_argument(1,dummy)
        if (len_trim(dummy)>0) read(dummy,fmt=*) n
        print *,"using n = ",n

        call get_command_argument(2,dummy)
        if (len_trim(dummy)>0) read(dummy,fmt=*) theta        
        print *,"using theta = ",theta


        allocate(x(3,n))

        if (theta==PI) then
            call fibonacci_sphere(x)
        else
            call fibonacci_cap(x,theta)
        end if
         
        open(unit=500,file="test_Lib_FibonacciSphere.xyz",action="write")
            write(unit=500,fmt='(i8)') n
            write(unit=500,fmt='(a)') "Lattice=""1 0 0 0 1 0 0 0 1"" Properties=species:S:1:pos:R:3"
            do ii = 1,n
                write(unit=500,fmt='(a,3f12.6)') "Du ",x(:,ii)
            end do            
        close(unit=500)

        dmin = 5.0d0
        if ((theta==PI).and.(n<=10000))  then
            dmin = huge(1.0)
            do ii = 2,n
                do jj = 1,ii-1
                    dx = x(:,jj) - x(:,ii)
                    dmin = min(dmin, dot_product(dx,dx))
                end do
            end do
            print *,"minimum separation ",sqrt(dmin),dmin*n
        end if

        ok = (dmin*n > 4.0d0 + 4.0d0/n ) 
        print *,dmin*n , 4.0d0 + 4.0d0/n , ok

        if (ok) then
            print *,"PASS"
        else
            print *,"FAIL"
        end if



        print *,""
        print *,"done"
        print *,""

    end program test_Lib_FibonacciSphere
