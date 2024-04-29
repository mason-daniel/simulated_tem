
    module Lib_GoldenSection
!---^^^^^^^^^^^^^^^^^^^^^^^^
!*      generic golden section search
!*      operation should be
!*

!*          y1 = func(x1) ; y3 = func(x3)                                       !   ideally x1,x3 bound the minimum, but algorithm checks.
!*          gold = GoldenSection_ctor( x1,x3, y1,y3)
!*          do
!*              call nextPoint(gold,x2)
!*              y2 = func(x2)
!*              call minimise(gold,x2,y2,isConverged,isWithinTimeLimit)
!*              if (isConverged) exit
!*              if (.not. isWithinTimeLimit) stop "failed"
!*          end do



        use iso_fortran_env
        implicit none
        private


    !---
        public      ::      GoldenSection_ctor
        public      ::      report
        public      ::      delete

        public      ::      nextPoint
        public      ::      minimise

    !---

        real(kind=real64),private,parameter     ::      GOLD = (1 + sqrt(5.0d0))/2              !    = 1.618034
        real(kind=real64),public                ::      GOLD_TOL = 1.0d-6
        real(kind=real64),public                ::      GOLD_MINSTEP = 50
        real(kind=real64),public                ::      GOLD_BOUNDSTEP = 10

        logical,public          ::      LIB_GOLD_DBG = .false.

    !---

        type,public     ::      GoldenSection
            private
            real(kind=real64)       ::      x1,x2,x3                !   points at which function is evaluated
            real(kind=real64)       ::      y1,y2,y3                !   f(x1)...
            logical                 ::      bounded
            integer                 ::      step                    !   number of function calcs
            real(kind=real64)       ::      xtol,ytol
        end type

    !---

        interface   GoldenSection_ctor
            module procedure    GoldenSection_null
            module procedure    GoldenSection_ctor0

        end interface

        interface   report
            module procedure    report0
        end interface

        interface   delete
            module procedure    delete0
        end interface

        interface   checkBounds
            module procedure    checkBounds0
        end interface

        interface   nextPoint
            module procedure    nextPoint0
        end interface

        interface   minimise
            module procedure    minimise0
        end interface


    contains
!---^^^^^^^^


        subroutine report0(this,u,o)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^
            type(GoldenSection),intent(inout)           ::      this
            integer,intent(in),optional                 ::      u,o
            integer                                     ::      uu,oo

            uu = 6 ; if (present(u)) uu = u
            oo = 0 ; if (present(o)) oo = o
            if (.not. this%bounded) then
                write(unit=uu,fmt='(a,i6,4(a,g16.6),a,2l6)') repeat(" ",oo),this%step," ",this%x1,repeat(" ",16),this%x3," ",this%y1,repeat(" ",16),this%y3," ",this%bounded
            else
                write(unit=uu,fmt='(a,i6,a,3g16.6,a,3g16.6,a,2l6)') repeat(" ",oo),this%step," ",this%x1,this%x2,this%x3," ",this%y1,this%y2,this%y3," ",this%bounded
            end if
            return
        end subroutine report0


        subroutine delete0(this)
    !---^^^^^^^^^^^^^^^^^^^^^^^^
            type(GoldenSection),intent(inout)         ::      this
            this = GoldenSection_ctor()
            return
        end subroutine delete0


        function GoldenSection_null() result(this)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
            type(GoldenSection)                       ::      this
            this%x1 = 0.0d0
            this%x3 = 1.0d0
            this%x2 = this%x1 + (this%x3 - this%x1)*(1-1/GOLD)
            this%y1 = 0.0d0
            this%y2 = 0.0d0
            this%y3 = 0.0d0
            this%step = 0
            this%bounded = .false.
            return
        end function GoldenSection_null


        function GoldenSection_ctor0(x1,x3 , y1,y3) result(this)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
            real(kind=real64),intent(in)                ::      x1,x3,y1,y3
            type(GoldenSection)                         ::      this
            if (x1<x3) then
                this%x1 = x1
                this%x3 = x3
                this%y1 = y1
                this%y3 = y3
            else
                this%x3 = x3
                this%x1 = x1
                this%y3 = y3
                this%y1 = y1
            end if
            this%x2 = 0.0d0
            this%y2 = huge(1.0)
            this%step = 0
            this%bounded = .false.
            this%xtol = (this%x3-this%x1)*GOLD_TOL
            this%ytol = max(GOLD_TOL,abs(this%y3-this%y1))*GOLD_TOL
            return
        end function GoldenSection_ctor0

    !---

        subroutine checkBounds0( this, xt,yt )
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      this routine is called if y1 or y3 is a minimum, or on first step where it isnt known
    !*      add the point (xt,yt). Is the function bounded now?
            type(GoldenSection),intent(inout)           ::      this
            real(kind=real64),intent(in)                ::      xt,yt

            if (this%step == 1) then
                if (LIB_GOLD_DBG) print *,"Lib_GoldenSection::checkBounds0 info - x1-xt-x3 = ",this%x1,xt,this%x3,":",this%y1,yt,this%y3
                if (yt < min(this%y1,this%y3)) then
                        !   y1    y3
                        !     \  /              BOUNDED
                        !      yt
                        !   x1-xt-x3
                    this%bounded = .true.
                    this%x2 = xt
                    this%y2 = yt
                else if ( (yt > this%y1).and. (this%y1 > this%y3) ) then
                        !      yt
                        !     /  \
                        !   y1    \             NOT BOUNDED- TRY xt-x3-
                        !          y3
                        !   x1-xt--x3
                    this%bounded = .false.
                    this%x1 = xt
                    this%y1 = yt
                else if ( (yt > this%y3).and. (this%y3 > this%y1) ) then
                        !      yt
                        !     /  \
                        !    /    y3             NOT BOUNDED- TRY -x1-xt
                        !   y1
                        !   x1-xt--x3
                    this%bounded = .false.
                    this%x3 = xt
                    this%y3 = yt
                else if ( (this%y1 > yt).and. (yt > this%y3) ) then
                        !   y1
                        !     \
                        !      yt        NOT BOUNDED- TRY xt-x3-
                        !       \
                        !        y3
                        !   x1-xt-x3
                    this%bounded = .false.
                    this%x1 = xt
                    this%y1 = yt
                else if ( (this%y3 > yt).and. (yt > this%y1) ) then
                        !         y3
                        !        /
                        !      yt        NOT BOUNDED- TRY -x1-xt
                        !     /
                        !   y1
                        !   x1-xt-x3
                    this%bounded = .false.
                    this%x3 = xt
                    this%y3 = yt
                end if
                return
            end if



            if (this%y1 < this%y3) then
                !   if y1 < y3, then I will have been offered a point xt-x1-x3
                if (LIB_GOLD_DBG) print *,"Lib_GoldenSection::checkBounds0 info - xt-x1-x3 = ",xt,this%x1,this%x3,":",yt,this%y1,this%y3
                if (yt < this%y1) then
                        !         y3
                        !        /
                        !      y1               NOT BOUNDED
                        !     /
                        !   yt
                        !   xt-x1-x3
                    this%bounded = .false.
                else
                        !   yt    y3
                        !     \  /              BOUNDED
                        !      y1
                        !   xt-x1-x3
                    this%bounded = .true.
                end if
                this%x2 = this%x1
                this%x1 = xt
                this%y2 = this%y1
                this%y1 = yt

            else
                !   if y3 < y1, then I will have been offered a point x1-x3-xt
                if (LIB_GOLD_DBG) print *,"Lib_GoldenSection::checkBounds0 info - x1-x3-xt = ",this%x1,this%x3,xt,":",this%y1,this%y3,yt
                if (yt < this%y3) then
                        !   y1
                        !     \
                        !      y3         NOT BOUNDED
                        !        \
                        !         yt
                        !   x1-x3-xt
                    this%bounded = .false.
                else
                        !   y1    yt
                        !     \  /              BOUNDED
                        !      y3
                        !   x1-x3-xt
                    this%bounded = .true.
                end if
                this%x2 = this%x3
                this%x3 = xt
                this%y2 = this%y3
                this%y3 = yt
            end if




            return
        end subroutine checkBounds0

        subroutine nextPoint0(this,xt)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      request position of next trial point
            type(GoldenSection),intent(inout)           ::      this
            real(kind=real64),intent(out)               ::      xt

!            print *,"nextPoint0 this%bounded,this%expanding ", this%bounded,this%expanding,this%step
            if (.not. this%bounded) then

                if (this%step == 0) then
                    !   this is the first call. We don't actually know that [x1,x3] does not bound the minimum ...
                    xt = this%x1 + (this%x3 - this%x1)* (GOLD-1)        !   note GOLD-1 = GOLD/(GOLD+1)
                    if (LIB_GOLD_DBG) print *,"Lib_GoldenSection::nextPoint0 info - not bounded, first call, try x1-xt-x3 = ",this%x1,xt,this%x3

                else
                    !   this is not the first call, and we know [x1,x2,x3] do not bound the minimum.

                    if (this%y1 < this%y3) then
                        !   look for a point xt-x1-x3
                        xt = this%x1 - (this%x3 - this%x1)/GOLD
                        if (LIB_GOLD_DBG) print *,"Lib_GoldenSection::nextPoint0 info - not bounded, y1<y3, try xt-x1-x3 = ",xt,this%x1,this%x3
                    else
                        !   look for a point x1-x3-xt
                        xt = this%x3 + (this%x3 - this%x1)/GOLD
                        if (LIB_GOLD_DBG) print *,"Lib_GoldenSection::nextPoint0 info - not bounded, y3<y1, try x1-x3-xt = ",this%x1,this%x3,xt
                    end if

                end if
            else
                !  y1                       y3
                !    \    y3         y1    /
                !     \  /      or     \  /
                !      y2               y2
                !  x1--x2-x3         x1-x2--x3

                xt = this%x1 + (this%x3 - this%x2)
            end if


            return
        end subroutine nextPoint0

        subroutine minimise0(this,xt,yt,is,ok)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      add the point (xt,yt). Is the function minimized now?
    !*      return minimum point if it is
    !*      returns ok = false if timed out
            type(GoldenSection),intent(inout)           ::      this
            real(kind=real64),intent(inout)             ::      xt,yt
            logical,intent(out)                         ::      is,ok

            this%step = this%step + 1
            !print *,"minimise step ",this%step," bounded ",this%bounded," trial ",xt,yt

            if (.not. this%bounded) then

                call checkBounds( this, xt,yt )
                is = .false.
                ok = this%bounded .or. ( this%step <= GOLD_BOUNDSTEP )
                return
            end if

            is = .false.
            ok = ( this%step <= GOLD_MINSTEP )
            if (.not. ok) return


            if (yt<this%y2) then
                !   have found a new minimum

                    !            y3                     y3          y1                    y1
                    !  y1       /            y1        /              \         y3          \     y3
                    !    \    y2        or     \      /       or       \       /             y2   /
                    !     \  /                  y2   /                  \    y2     or        \  /
                    !      yt                    \  /                    \  /                  yt
                    !                             yt                      yt
                    !  x1--xt-x2-x3          x1-x2-xt---x3          x1---xt-x2-x3        x1-x2-xt--x3
                if (xt<this%x2) then
                    this%x3 = this%x2
                    this%y3 = this%y2
                    this%x2 = xt
                    this%y2 = yt
                else
                    this%x1 = this%x2
                    this%y1 = this%y2
                    this%x2 = xt
                    this%y2 = yt
                end if

           else
               !   have found a new bound

                !            y3                     y3       y1                    y1
                !  y1       /            y1        /           \         y3          \     y3
                !    \    yt        or     \      /             \       /             y2   /
                !     \  /                  yt   /       or      \    y2     or        \  /
                !      y2                    \  /                 \  /                  yt
                !                             y2                   yt
                !  x1--x2-xt-x3          x1-xt-x2---x3       x1---xt-x2-x3        x1-x2-xt--x3
                if (xt>this%x2) then
                    this%x3 = xt
                    this%y3 = yt
                else
                    this%x1 = xt
                    this%y1 = yt
                end if
            end if

            is = ( abs(this%x3 - this%x1) < this%xtol ).or.( abs(this%y2 - max(this%y1,this%y3)) < this%ytol )

            if (is) then
                xt = this%x2 ; yt = this%y2
            end if
            !this%step = this%step + 1

            return
        end subroutine minimise0





    end module Lib_GoldenSection


!     program testLib_GoldenSection
! !---^^^^^^^^^^^^^^^^^^^^
!         use Lib_Gold
!         use iso_fortran_env
!         implicit none

!         integer                 ::      ii
!         real(kind=real64)       ::      xx

!         real(kind=real64)       ::      xt,yt ,x1,x3
!         type(GoldenSection)     ::      gold
!         logical     ::          is,ok

!     !---    demonstrate presence of actual minimum

!         do ii = 1,10000
!             xx = ii*0.001d0
!             if ( func(xx)<min(func(xx-1.0d-3),func(xx+1.0d-3)) ) then
!                 print *,"minimum of func found at ",xx,func(xx)
!             end if
!         end do


!    !----    guess between 1 and 1.5
!         x1 = 1.0d0 ; x3 = 1.5d0
!         print *,"guess between ",x1,x3
!         gold = GoldenSection_ctor( x1,x3 , func(x1),func(x3) )
!         call report(gold)
!         do
!             call nextPoint(gold,xt)
!             yt = func(xt)
!             call minimise(gold,xt,yt ,is,ok)
!             call report(gold)
!             if (is) exit
!             if (.not. ok) then
!                 print *,"testLib_Gold error - minimise failed"
!                 exit
!             end if
!         end do
!         print *,"golden section result ",xt,yt
!         print *,""


!    !----    guess between 1.5 and 2
!         x1 = 1.5d0 ; x3 = 2.0d0
!         print *,"guess between ",x1,x3
!         gold = GoldenSection_ctor( x1,x3 , func(x1),func(x3) )
!         call report(gold)
!         do
!             call nextPoint(gold,xt)
!             yt = func(xt)
!             call minimise(gold,xt,yt ,is,ok)
!             call report(gold)
!             if (is) exit
!             if (.not. ok) then
!                 print *,"testLib_Gold error - minimise failed"
!                 exit
!             end if
!         end do
!         print *,"golden section result ",xt,yt
!         print *,""

!    !----    guess between 2 and 2.5
!         x1 = 2.0d0 ; x3 = 2.5d0
!         print *,"guess between ",x1,x3
!         gold = GoldenSection_ctor( x1,x3 , func(x1),func(x3) )
!         call report(gold)
!         do
!             call nextPoint(gold,xt)
!             yt = func(xt)
!             call minimise(gold,xt,yt ,is,ok)
!             call report(gold)
!             if (is) exit
!             if (.not. ok) then
!                 print *,"testLib_Gold error - minimise failed"
!                 exit
!             end if
!         end do
!         print *,"golden section result ",xt,yt
!         print *,""

!    !----    guess between 3 and 4
!         x1 = 3.0d0 ; x3 = 4.0d0
!         print *,"guess between ",x1,x3
!         gold = GoldenSection_ctor( x1,x3 , func(x1),func(x3) )
!         call report(gold)
!         do
!             call nextPoint(gold,xt)
!             yt = func(xt)
!             call minimise(gold,xt,yt ,is,ok)
!             call report(gold)
!             if (is) exit
!             if (.not. ok) then
!                 print *,"testLib_Gold error - minimise failed"
!                 exit
!             end if
!         end do
!         print *,"golden section result ",xt,yt
!         print *,""

!    !----    guess between 4 and 5
!         x1 = 4.0d0 ; x3 = 5.0d0
!         print *,"guess between ",x1,x3
!         gold = GoldenSection_ctor( x1,x3 , func(x1),func(x3) )
!         call report(gold)
!         do
!             call nextPoint(gold,xt)
!             yt = func(xt)
!             call minimise(gold,xt,yt ,is,ok)
!             call report(gold)
!             if (is) exit
!             if (.not. ok) then
!                 print *,"testLib_Gold error - minimise failed"
!                 exit
!             end if
!         end do
!         print *,"golden section result ",xt,yt
!         print *,""

!    !----    guess between 5 and 6
!         x1 = 5.0d0 ; x3 = 6.0d0
!         print *,"guess between ",x1,x3
!         gold = GoldenSection_ctor( x1,x3 , func(x1),func(x3) )
!         call report(gold)
!         do
!             call nextPoint(gold,xt)
!             yt = func(xt)
!             call minimise(gold,xt,yt ,is,ok)
!             call report(gold)
!             if (is) exit
!             if (.not. ok) then
!                 print *,"testLib_Gold error - minimise failed"
!                 exit
!             end if
!         end do
!         print *,"golden section result ",xt,yt
!         print *,""



!    !----    guess between 1 and 3
!         x1 = 1.0d0 ; x3 = 3.0d0
!         print *,"guess between ",x1,x3
!         gold = GoldenSection_ctor( x1,x3 , func(x1),func(x3) )
!         call report(gold)
!         do
!             call nextPoint(gold,xt)
!             yt = func(xt)
!             call minimise(gold,xt,yt ,is,ok)
!             call report(gold)
!             if (is) exit
!             if (.not. ok) then
!                 print *,"testLib_Gold error - minimise failed"
!                 exit
!             end if
!         end do
!         print *,"golden section result ",xt,yt
!         print *,""




!         print *,""
!         print *,"done"
!         print *,""


!     contains
! !---^^^^^^^^

!         pure real(kind=real64) function func(x)
!     !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

!             real(kind=real64),intent(in)        ::      x

!             if (x<0) then
!                 func = huge(1.0)
!             else
!                 func = 1.0d0 + x*(1.0d0 + x*(-0.1d0 + x*0.2d0)) * sin(3*x) + cos(x)*sqrt(x)
!             end if
!             return
!         end function func
!     end program testLib_GoldenSection