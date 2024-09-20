
    module Lib_Callipers
!---^^^^^^^^^^^^^^^^^^^
!*      simple timer function
!*      It is not expected that this code will have any serious impact on scientific results,
!*      it is really just to be able to time how long sections of code take to run.
!*
!*      Daniel Mason
!*      UKAEA
!*      April 2022
!*

        use iso_fortran_env
        implicit none
        private
        
    !---
    
        public      ::      Callipers_ctor  !   create new stopwatch
        public      ::      report          !   report status of stopwatch
        public      ::      delete          !   does nothing - no dynamic memory
        
        public      ::      start           !   start the stopwatch
        public      ::      pause           !   pause the stopwatch
        public      ::      reset           !   reset the stopwatch
        public      ::      elapsed         !   return time taken in s
        
    !---
    
        type,public     ::      Callipers
            integer(kind=int64)     ::      t,t_elapsed
        end type        
        
    !--- 
    
        interface   start
            module procedure        start0
        end interface
        
        interface   pause
            module procedure        pause0
        end interface
        
        interface   reset
            module procedure        reset0
        end interface
        
    
        interface   report
            module procedure        report0
        end interface
        
        interface   delete
            module procedure        delete0
        end interface
        
        interface   elapsed
            module procedure        elapsed0
        end interface
        
    contains
!---^^^^^^^^

        function Callipers_ctor() result(this)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^    
            type(Callipers)                     ::      this
            call system_clock( count = this%t )
            this%t_elapsed = 0
            return
        end function Callipers_ctor
        
        subroutine delete0(this)
    !---^^^^^^^^^^^^^^^^^^^^^^^^
            type(Callipers),intent(inout)       ::      this
            this%t = 0
            this%t_elapsed = 0
            return
        end subroutine delete0
        
        subroutine start0(this)
    !---^^^^^^^^^^^^^^^^^^^^^^^
            type(Callipers),intent(inout)       ::      this            
            call system_clock( count = this%t )
            return
        end subroutine start0
        
        subroutine pause0(this)
    !---^^^^^^^^^^^^^^^^^^^^^^^
            type(Callipers),intent(inout)       ::      this
            integer(kind=int64)                 ::      tt
            call system_clock( count = tt )
            this%t_elapsed = this%t_elapsed + (tt - this%t)
            return
        end subroutine pause0
        
        subroutine reset0(this)
    !---^^^^^^^^^^^^^^^^^^^^^^^
            type(Callipers),intent(inout)       ::      this
            call system_clock( count = this%t )
            this%t_elapsed = 0
            return
        end subroutine reset0
        
        subroutine report0(this,u,o)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^
            type(Callipers),intent(in)          ::      this
            integer,intent(in),optional         ::      u,o
            integer     ::      uu,oo
            uu = 6 ; if (present(u)) uu = u
            oo = 0 ; if (present(o)) oo = o
            write(unit=uu,fmt='(a,f24.9,a)') repeat(" ",oo)//"Callipers [t=",elapsed0(this)," s]"
            return
        end subroutine report0
            
        function elapsed0(this) result(t)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      returns the elapsed time in seconds. 
    !*      note that it asks the system clock what the count rate is.
            type(Callipers),intent(in)          ::      this
            real(kind=real64)                   ::      t
            integer(kind=int64) ::      tt,rr
            call system_clock( count = tt, count_rate = rr )
            
            if (this%t_elapsed == 0) then
                tt = tt - this%t
            else
                tt = this%t_elapsed
            end if
            
            if (rr == 1000) then
                t = tt*1.0d-3
            else if (rr == 1000000) then
                t = tt*1.0d-6
            else if (rr == 1000000000) then
                t = tt*1.0d-9
            else
                t = tt
            end if
                
            return
        end function elapsed0
        
        
    end module Lib_Callipers        