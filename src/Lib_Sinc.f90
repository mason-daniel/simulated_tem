
    module Lib_Sinc
!---^^^^^^^^^^^^^^^
!*      simple implementation of a regular sinc function
!*          sinc(x) = sin(pi x)/(pi x)

        use iso_fortran_env
        implicit none
        
        real(kind=real64),private,parameter     ::      EPS = 1.0d-12
        real(kind=real64),private,parameter     ::      PI = 3.14159265358979d0
        real(kind=real64),private,parameter     ::      IPI = 1/PI
        
        public      ::      sinc
        public      ::      isinc
        
        
    contains
!---^^^^^^^^


        elemental function sinc( x ) result( y )
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      returns sin(pi x)/(pi x)
    
            real(kind=real64),intent(in)        ::      x
            real(kind=real64)                   ::      y
            
            y = PI*x
            if (abs(y)<EPS) then
                y = 1.0d0 - y*y/6
            else
                y = sin(y)/y
            end if
            
            return
        end function sinc
        
        
                   
        elemental function isinc( y ) result( x )
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      returns 0<=x<=1, where y = sin(pi x)/(pi x)
            real(kind=real64),intent(in)        ::      y
            real(kind=real64)                   ::      x
             
            real(kind=real64)       ::      dx,dy,ss,cc,ix
            integer                 ::      ii
            
        !---    simple cases are easy to return
            if (y>=1) then
                x = 0.0d0
            else if (y<=0) then
                x = 1.0d0
            else if (y>1-EPS) then
                x = sqrt( 6 - 6*y )*IPI
            else if (y<EPS) then
                x = 1 - y
            else    
                !   use a Pade approximation to find sinc, 
                !   and an inverse to find an approximate value
                x = 6*(sqrt(35.0d0)*sqrt( -3985*y*y+130862*y+25583)-455*y-1855)/(75*y-551)
                x = sqrt(x)     !   note this is y = sin(pi x)/(pi x)
                
                
                !   now use Newton-Raphson to find an improved solution.
                !   note that I know the derivative is non zero.                
                do ii = 1,20
                    ss = sin(x) ; cc = cos(x) ; ix = 1/x
                    ss = ss*ix
                    if (abs(ss - y)<1.0d-12) exit                
                    cc = cc*ix
                    dy = cc - ss*ix
                    dx = ( y - ss )/dy
                    x  = x + dx
                    !print *,(sin(x+1.0d-8)/(x+1.0d-8) - sin(x-1.0d-8)/(x-1.0d-8))/(2d-8),dy
                end do
                
                x  = IPI*(x+dx)
                
            end if
            
            
            return            
        end function isinc
        
        
    end module Lib_Sinc            

!   gfortran -ffree-line-length-256 ${MYF90LIB}/Lib_Sinc.f90      
!                
!    program testLib_Sinc
!!---^^^^^^^^^^^^^^^^^^^^
!        use Lib_Sinc
!        use iso_fortran_env
!        implicit none
!        
!        character(len=256)      ::      dummy
!        real(kind=real64)       ::      xx,yy,zz
!        
!        call get_command_argument( 1,dummy )
!        read(dummy,fmt=*)  xx
!        
!        yy = sinc(xx)
!        zz = isinc(yy)
!        print *,"xx = ",xx," sinc( pi x ) ",yy," isinc( sinc( pi x ) ) ",zz," err ",abs(zz-xx)
!        
!    end program testLib_Sinc
        