
    module Lib_SimpleProgressBar
!---^^^^^^^^^^^^^^^^^^^^^^^^^^^^
!*
!*      Very simple module which outputs to screen progress through a long operation
!*      usage:
!*          do i = 1,n
!*              call progressBar( i,n )
!*          end do
!*      don't try to use i<1 by mistake.
!*
!*      Daniel Mason, UKAEA
!*      April 2022
!*

        use iso_fortran_env
        private
        
        
        public      ::      progressBar
        
        
        interface       progressBar       
            module procedure        progressBar_32
            module procedure        progressBar_64
        end interface
        
    contains
!---^^^^^^^^

        subroutine progressBar_32( i,n )
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      output 10%,20% etc as iteration counter i gets towards n
            integer(kind=int32),intent(in)          ::      i       !   1<=i<=n
            integer(kind=int32),intent(in)          ::      n       !   end point 100%
            
            integer             ::  non10,pct
            
            non10 = max(1,int(n*0.1d0))
            
            pct = 0
            
            if (i == 1) write(*,fmt='(i3,a)',advance="no") pct,"% " 
            
            if ( mod(i,non10) == 0) then
                pct = int( i/real(non10) )*10  
                if (pct == 100) then
                    write(*,fmt='(i3,a)',advance="yes") pct,"% "
                    return
                else
                    write(*,fmt='(i3,a)',advance="no") pct,"% "
                end if
       
            end if 
            
            return
        end subroutine progressBar_32
        
        
        subroutine progressBar_64( i,n )
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      output 10%,20% etc as iteration counter i gets towards n
            integer(kind=int64),intent(in)          ::      i       !   1<=i<=n
            integer(kind=int64),intent(in)          ::      n       !   end point 100%
            
            integer(kind=int64)             ::  non10
            integer                         ::  pct
            
            non10 = max(1,int(n*0.1d0))
            
            pct = 0
            
            if (i == 1) write(*,fmt='(i3,a)',advance="no") pct,"% " 
            
            if ( mod(i,non10) == 0) then
                pct = int( i/real(non10,kind=real64),kind=int32 )*10  
                if (pct == 100) then
                    write(*,fmt='(i3,a)',advance="yes") pct,"% "
                    return
                else
                    write(*,fmt='(i3,a)',advance="no") pct,"% "
                end if
       
            end if
             
            
            return
        end subroutine progressBar_64
                
    end module Lib_SimpleProgressBar        