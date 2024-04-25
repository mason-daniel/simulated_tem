
    module Lib_QuinticSplinesEven
!---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^  
!*      A quintic spline is an extension of the more familiar cubic spline,
!*      but with continuous derivatives up to fourth order ( cf cubic to 2nd ).
!*      The mathematics of the construction is very similar to that of the cubic spline
!*      ( see numerical recipes section 3.3 )
!*      Write the fourth derivative as continuous linear variation
!*          y"" = A v_i + B v_i+1
!*      with A = (x_i+1 - x)/(x_i+1-x_i) and B = 1-A
!*      then integrate.
!*      We expect to be given {x_i} and {y_i} as input, so set constants with
!*      condition y(x_i) = y_i and y"(x_i) = u_i
!*      giving
!*
!*          y = A y_i + B y_i+1
!*            + C u_i + D u_i+1
!*            + E v_i + F v_i+1
!*      with
!*          C = A/6 ( A^2 - 1 ) ( x_i+1 - x_i )^2
!*          D = B/6 ( B^2 - 1 ) ( x_i+1 - x_i )^2
!*          E = A/360 ( 3A^4 - 10A^2 + 7 ) ( x_i+1 - x_i )^4
!*          F = B/360 ( 3B^4 - 10B^2 + 7 ) ( x_i+1 - x_i )^4
!*
!*
!*      Complete by solving for {u_i} and {v_i} by insisting y',y"' continuous.
!*      ( note: y,y",y"" are continuous by construction )
!*
!*      It is more expensive to compute than cubic, as the matrix equation that you end
!*      up with is band diagonal, rather than tridiagonal.
!*      However, once the coefficients u_i and v_i are found, it is quick to interpolate
!*
!*      This version assumes evenly spaced x_i. This avoids a division, and so speeds up interpolation somewhat.
!*
!*      Author      :   Daniel Mason
!*      Version     :   1.0
!*      Revision    :   Feb 2012
!*
!*      note        :   must be compiled against LAPACK
!*
        use NBAX_StringTokenizers
        use Lib_NBAX
        use iso_Fortran_env
        implicit none
        private
!         include "mpif.h"

        external    ::      DGESV,SGESV             !**********     NOTE - ISSUES WITH DGESV - WORKAROUND WITH SGESV ************!

    !---

        public      ::      QuinticSplineEven_ctor
        public      ::      report
        public      ::      delete
        public      ::      clear
        public      ::      clone           !   deep copy with allocation clone(this,that) => this = that

        public      ::      refitQSpline
        public      ::      splint
        public      ::      qsplint

        public      ::      qsplint_inrange
        public      ::      qsplint_array_inrange
        public      ::      qsplint_array_inrange_sum


        public      ::      getXmin,setXmin
        public      ::      getXmax,setXmax

        public      ::      getnKnots
        public      ::      gety,sety
        public      ::      isNull

        public      ::      outputAsXML
        public      ::      inputFromXML

!         public      ::      MPI_exchange

    !---

        type,public     ::      QuinticSplineEven
            private
            integer                                 ::          N       !   number of knots
            real(kind=real64)                       ::          xmin,xmax,dx,idx        !   dx = (xmax-xmin)/(N-1) ; idx = (N-1)/(xmax-xmin)
            integer,dimension(2,4)                  ::          key     !   what to do with derivatives at end
            real(kind=real64),dimension(2,4)        ::          value   !   fixed values of derivatives at end ( where appropriate )
            real(kind=real64),dimension(:,:),pointer  ::          y       !   (1:3,1:N) value at knots (1:3 = (/y,(dx^2/6)*y",(dx^4/360)*y""/))
!             real(kind=real64),dimension(:),pointer  ::          u       !   y"(1:N) value at knots
!             real(kind=real64),dimension(:),pointer  ::          v       !   y""(1:N) value at knots
        end type QuinticSplineEven

    !---

        interface   QuinticSplineEven_ctor
            module procedure    QuinticSplineEven_null
            module procedure    QuinticSplineEven_ctor0
            module procedure    QuinticSplineEven_ctor1
        end interface

        interface   report
            module procedure    report0
        end interface

        interface   delete
            module procedure    delete0
        end interface

        interface   clear
            module procedure    clear0
        end interface

        interface   clone
            module procedure    clone0
        end interface
        
        interface   gety
            module procedure    gety0
        end interface

        interface   sety
            module procedure    sety0
        end interface

        interface  isNull
            module procedure    isNull0
        end interface

        interface   qsplint
            module procedure    qsplint0
            module procedure    qsplint0a
            module procedure    qsplint1
            module procedure    qsplint2
        end interface

        interface   qsplint_inrange
            module procedure    qsplint_inrange0
            module procedure    qsplint_inrange0a
            module procedure    qsplint_inrange1
            module procedure    qsplint_inrange2
        end interface


        interface   qsplint_array_inrange
            module procedure    qsplint_array_inrange0
            module procedure    qsplint_array_inrange0a
            module procedure    qsplint_array_inrange0b
        end interface

        interface   qsplint_array_inrange_sum
            module procedure    qsplint_array_inrange_sum0
            module procedure    qsplint_array_inrange_sum1
        end interface

        interface   inputFromXML
            module procedure    inputFromXML0
        end interface

        interface   outputAsXML
            module procedure    outputAsXML0
        end interface

!         interface   MPI_exchange
!             module procedure    MPI_exchange0
!         end interface



    !---

        integer,public,parameter        ::      QS_DERIV_FREE = 0
        integer,public,parameter        ::      QS_DERIV_FIXED = -1

    contains
!---^^^^^^^^

        function QuinticSplineEven_null( ) result( this )
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      quintic spline memory allocation and fitting
            type(QuinticSplineEven)                     ::      this
            this%N = 0
            nullify(this%y)
!             nullify(this%u)
!             nullify(this%v)
            return
        end function QuinticSplineEven_null

        function QuinticSplineEven_ctor0( N,xmin,xmax,y  ,key,value) result( this )
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      quintic spline memory allocation and fitting

    !*      on input key(1:2,1:4) defines what to do with the derivatives at either end
    !*              and can be FREE,FIXED ( to a value given ) or interpolated to order key = 1,2,3,4
    !*              ( suggest always used accuracy order = 4 unless there is a good reason not to )
    !*      value ( if present ) defines up to four FIXED derivatives
    !*      note: only four derivs max may be set, this is the limit of the degrees of freedom.

            integer,intent(in)                                  ::      N
            real(kind=real64),intent(in)                        ::      xmin,xmax
            real(kind=real64),dimension(:),intent(in)           ::      y
            integer,dimension(2,4),intent(in)                   ::      key     !   what to do with derivatives at end
            real(kind=real64),dimension(:),intent(in),optional  ::      value


            type(QuinticSplineEven)                     ::      this

            real(kind=real64),dimension(0:7,4,4),parameter        ::      FINITE_DIFF_COEFF = reshape( (/                         &
                            -1.0d0,     1.0d0,      0.0d0,      0.0d0,      0.0d0,      0.0d0,      0.0d0,      0.0d0,          &
                            -1.5d0,     2.0d0,      -.5d0,      0.0d0,      0.0d0,      0.0d0,      0.0d0,      0.0d0,          &
                       -11.d0/6.d0,     3.0d0,     -1.5d0,  1.d0/3.d0,      0.0d0,      0.0d0,      0.0d0,      0.0d0,          &
                      -25.d0/12.d0,     4.0d0,     -3.d0,   4.d0/3.d0,     -.25d0,      0.0d0,      0.0d0,      0.0d0,          &
                             1.0d0,     -2.d0,      1.0d0,      0.0d0,      0.0d0,      0.0d0,      0.0d0,      0.0d0,          &
                             2.0d0,     -5.d0,      4.0d0,      -1.d0,      0.0d0,      0.0d0,      0.0d0,      0.0d0,          &
                         35./12.d0, -26./3.d0,   19./2.d0,  -14./3.d0,   11./12.d0,     0.0d0,      0.0d0,      0.0d0,          &
                          15./4.d0, -77./6.d0,  107./6.d0,     -13.d0,   61./12.d0,  -5./6.d0,      0.0d0,      0.0d0,          &
                            -1.0d0,     3.0d0,     -3.0d0,      1.0d0,      0.0d0,      0.0d0,      0.0d0,      0.0d0,          &
                            -2.5d0,     9.0d0,    -12.0d0,      7.0d0,     -1.5d0,      0.0d0,      0.0d0,      0.0d0,          &
                         -17./4.d0,  71./4.d0,  -59./2.d0,   49./2.d0,   41./4.d0,    7./4.d0,      0.0d0,      0.0d0,          &
                         -49./8.d0,    29.0d0, -461./8.d0,      62.d0, -307./8.d0,     13.0d0,  -15./8.d0,      0.0d0,          &
                             1.0d0,    -4.0d0,      6.0d0,     -4.0d0,      1.0d0,      0.0d0,      0.0d0,      0.0d0,          &
                             3.0d0,   -14.0d0,     26.0d0,    -24.0d0,     11.0d0,     -2.0d0,      0.0d0,      0.0d0,          &
                          35./6.d0,   -31.0d0,  137./2.d0, -242./3.d0,  107./2.d0,    -19.0d0,   17./6.d0,      0.0d0,          &
                          28./3.d0,-111./2.d0,    142.0d0,-1219./6.d0,    176.0d0, -185./2.d0,   82./3.d0,   -7./2.d0 /),(/8,4,4/) )

            integer     ::      ii,jj,kk,ll
            real(kind=real64)   ::      xx

        !---    sanity check: are the derivatives sensible?
            ii = count( key(:,:)==QS_DERIV_FREE )
            if (ii < 4) then
                write (unit=0,fmt='(a,8i6)') "key=",key
                stop "Lib_QuinticSplinesEven::QuinticSplineEven_ctor0 error - derivatives at ends overspecified"
            end if
            if (any( key(:,:)>4 )) then
                write (unit=0,fmt='(a,8i6)') "key=",key
                stop "Lib_QuinticSplinesEven::QuinticSplineEven_ctor0 error - finite difference interpolated derivatives over maximum accuracy of 4"
            end if
            ii = count( key(:,:)==QS_DERIV_FIXED )
            if (ii>0) then
                if (.not. present(value)) then
                    write (unit=0,fmt='(a,8i6)') "key=",key
                    stop "Lib_QuinticSplinesEven::QuinticSplineEven_ctor0 error - fixed derivative values not given as input"
                else if (size(value)<ii) then
                    write (unit=0,fmt='(a,8i6)') "key=",key
                    stop "Lib_QuinticSplinesEven::QuinticSplineEven_ctor0 error - insufficient fixed derivative values given as input"
                end if
            end if

            this%N = N
            if (isNull(this)) then
                this = QuinticSplineEven_ctor()
                return
            end if
            if (size(y) < this%N) then
                write (unit=0,fmt='(a,2i12)') "size(y) , N",size(y),this%N
                stop "Lib_QuinticSplinesEven::QuinticSplineEven_ctor0 error - size(y) < N"
            end if


            allocate( this%y(3,this%N) )
!             allocate( this%y(2,this%N) )
!             allocate( this%y(3,this%N) )

            this%xmin = xmin
            this%xmax = xmax
            this%dx = (this%xmax - this%xmin)/(this%N-1)
            this%idx = (this%N-1)/(this%xmax - this%xmin)
            this%y = 0.0d0
            do ii = 1,this%N
                this%y(1,ii) = y(ii)
            end do

!             print *,"dx = ",this%dx
!
!             do jj = 1,2
!                 do ii = 1,4         !   derivative
!                     do kk = 1,4     !   accuracy
!                         xx = 0.0d0
!                         if (jj==1) then
!                             do ll = 0,7
!                                 xx = xx + FINITE_DIFF_COEFF(ll, kk, ii)*this%y(1,ll+1)
!                             end do
!                         else
!                             do ll = 0,7
!                                 xx = xx + FINITE_DIFF_COEFF(ll, kk, ii)*this%y(1,this%N-ll)
!                             end do
!                         end if
!                         print *,jj,ii,kk, xx*(this%idx**ii)
!                     end do
!                 end do
!             end do
!
!
!             print *,key,count( this%key(:,:)==QS_DERIV_FREE )

            this%key = key
            if (count( this%key(:,:)==QS_DERIV_FREE )>4) then                   !   note: could be 5,6,7,8
                ! underfitted
                if (count( this%key(1,:)==QS_DERIV_FREE )==4) then
                    !   no fitting on left hand end:
                    this%key(1,1:4) = (/ QS_DERIV_FREE,4,QS_DERIV_FREE,4 /)
                else if (count( this%key(1,:)==QS_DERIV_FREE )==3) then
                    !   single fitting on left hand end: improve this
                    if (this%key(1,1)==QS_DERIV_FREE) then
                        this%key(1,1) = 4
                    else if (this%key(1,2)==QS_DERIV_FREE) then
                        this%key(1,2) = 4
                    end if
                else
                    !   at least double fitting on left hand end. It is the right end that is the problem.
                end if
            end if
            if (count( this%key(:,:)==QS_DERIV_FREE )>4) then
                ! underfitted
                if (count( this%key(2,:)==QS_DERIV_FREE )==4) then
                    !   no fitting on right hand end:
                    this%key(2,1:4) = (/ QS_DERIV_FREE,4,QS_DERIV_FREE,4 /)
                else if (count( this%key(2,:)==QS_DERIV_FREE )==3) then
                    !   single fitting on right hand end: improve this
                    if (this%key(2,1)==QS_DERIV_FREE) then
                        this%key(2,1) = 4
                    else if (this%key(2,2)==QS_DERIV_FREE) then
                        this%key(2,2) = 4
                    end if
                end if
            end if

            kk = 0
            do ii = 1,4
                do jj = 1,2
                    if (this%key(jj,ii) == QS_DERIV_FREE) then
                        cycle
                    else if (this%key(jj,ii) == QS_DERIV_FIXED) then
                        kk = kk + 1
                        this%value(jj,ii) = value(kk)
                        ! print *,"fixed derivative ",ii," at end ",jj," is ",this%value(jj,ii)
                    else
                        xx = 0.0d0
                        if (jj==1) then
                            do ll = 0,7
                                xx = xx + FINITE_DIFF_COEFF(ll, this%key(jj,ii), ii)*this%y(1,ll+1)
                            end do
                        else
                            do ll = 0,7
                                xx = xx + FINITE_DIFF_COEFF(ll, this%key(jj,ii), ii)*this%y(1,this%N-ll)
                            end do
                        end if
                        this%value(jj,ii) = xx*(this%idx**ii)
                        ! print *,"interpolated derivative ",ii," at end ",jj," is ",this%value(jj,ii)
                    end if
                end do
            end do

            call refitQSpline(this)

            return
        end function QuinticSplineEven_ctor0



        function QuinticSplineEven_ctor1( N,xmin,xmax,y) result( this )
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      quintic spline memory allocation and fitting

    !*      on input key(1:2,1:4) defines what to do with the derivatives at either end
    !*              and can be FREE,FIXED ( to a value given ) or interpolated to order key = 1,2,3,4
    !*      value ( if present ) defines up to four FIXED derivatives
    !*      note: only four derivs max may be set, this is the limit of the degrees of freedom.

            integer,intent(in)                                  ::      N
            real(kind=real64),intent(in)                        ::      xmin,xmax
            real(kind=real64),dimension(:),intent(in)           ::      y
            type(QuinticSplineEven)                             ::      this
            integer,dimension(2,4)                      ::      key     !   what to do with derivatives at end
            key = QS_DERIV_FREE
            key(1,2) = 4
            key(1,4) = 4
            key(2,2) = 4
            key(2,4) = 4
            this = QuinticSplineEven_ctor0( N,xmin,xmax,y , key)
            return
        end function QuinticSplineEven_ctor1

        subroutine refitQSpline( this )
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      perform quintic spline fitting.
    !*      currently do it the slow way with full LAPACK linear equation solver
    !*              ie  solve A.B = C
    !*      this could be improved dramatically with some thought.
    !*
    !*      B order (u1,v1,u2,v2...uN,vN)




            type(QuinticSplineEven),intent(inout)       ::      this

            real(kind=real64),dimension(2*this%N,2*this%N)       ::      A
            real(kind=real64),dimension(2*this%N)                ::      B
            integer,dimension(2*this%N)             ::      IPIV
            integer                                 ::      NRHS,LDA,LDB,INFO
             

            integer             ::      ii,jj,kk
            real(kind=real64)                ::      dx

        !---    sanity check: arrays allocated?
            if (isNull(this)) then
                write (unit=0,fmt='(a)') "Lib_QuinticSplinesEven::refit - warning refitting unallocated spline?"
                return
            end if

        !---    sanity check: x in ascending order
            if ( this%xmin >= this%xmax) then
                write (unit=0,fmt='(a,i6,2g16.8)') "N,xmin,xmax=",this%N,this%xmin,this%xmax
                stop "Lib_QuinticSplinesEven::refit - error xmin >= xmax"
            end if


            !print *,"Lib_QuinticSplinesEven::refitQSpline info - this%N,this%dx = ",this%N,this%dx,this%idx

            NRHS = 1
            LDA = 2*this%N
            LDB = 2*this%N

            A = 0.0d0
            B = 0.0d0


            kk = 0              !   position of constraints
            do jj = 1,2     !   end
                do ii = 1,4         !   derivative
                    if (this%key(jj,ii) /= QS_DERIV_FREE) then
                        kk = kk + 1  ; if (kk==3) kk = 2*this%N-1
                        select case(ii)
                            case(1)
                                if (jj==1) then
                                    A( kk,1:4 )                 = (this%dx/360.0d0)*(/ -2.0d0 ,  8*this%dx , -1.0d0 ,  7*this%dx /)
                                    B(kk) = this%value(jj,ii) - this%idx*(this%y(1,2)-this%y(1,1))
                                else
                                    A( kk,2*this%N-3:2*this%N ) = (this%dx/360.0d0)*(/  1.0d0 , -7*this%dx ,  2.0d0 , -8*this%dx /)
                                    B(kk) = this%value(jj,ii) - this%idx*(this%y(1,this%N)-this%y(1,this%N-1))
                                end if
                            case(2)
                                if (jj==1) then
                                    A( kk,1 )                   = 1.0d0
                                else
                                    A( kk,2*this%N-1 )          = 1.0d0
                                end if
                                B(kk) = this%value(jj,ii)
                            case(3)
                                if (jj==1) then
                                    A( kk,1:4 )                 = (/ -6*this%idx ,   -2*this%dx , 6*this%idx ,  -1*this%dx /)/6.0d0
                                else
                                    A( kk,2*this%N-3:2*this%N ) = (/ -6*this%idx ,    1*this%dx , 6*this%idx ,   2*this%dx /)/6.0d0
                                end if
                                B(kk) = this%value(jj,ii)
                            case(4)
                                if (jj==1) then
                                    A( kk,2 )                   = 1.0d0
                                else
                                    A( kk,2*this%N )            = 1.0d0
                                end if
                                B(kk) = this%value(jj,ii)
                        end select
                    end if
                end do
            end do



            do ii = 2,this%N-1

            !---    fix first derivatives at x_i
                A( 2*ii-1, 2*ii-3:2*ii+2 ) = (/ 60*this%dx  ,-7*this%dx*this%dx*this%dx,            &
                                                240*this%dx ,-16*this%dx*this%dx*this%dx,           &
                                                60*this%dx  ,-7*this%dx*this%dx*this%dx /)
                B( 2*ii-1 ) = 360*this%idx*( (this%y(1,ii+1) - 2*this%y(1,ii) + this%y(1,ii-1)) )

            !---    fix third derivatives at x_i
                A( 2*ii  , 2*ii-3:2*ii+2 ) = (/ -6*this%idx ,this%dx,                               &
                                                12*this%idx ,4*this%dx,                             &
                                                -6*this%idx ,this%dx /)
                B( 2*ii ) = 0.0d0

            end do
!             A(2*this%N-1,2*this%N-3:2*this%N)   &
!                                     = (/ -6.0d0,this%dx*this%dx,6.0d0,0.0d0 /)
!             A(2*this%N,2*this%N)    = 1


        !---    for debugging only, here is the matrix and vector. I think they're OK.
        !                do ii = 1,min(20,2*this%N)
        !                    write (*,fmt='(1000f10.2)',advance="no") A(ii,1:min(20,2*this%N))
        !                    write (*,fmt='(a,g10.2)') "            ",B(ii)
        !                end do
        !---

            call DGESV( 2*this%N, NRHS, A, LDA, IPIV, B, LDB, INFO )
            
            
            if (INFO /= 0) then
                write (unit=0,fmt='(a,i6)') "Lib_QuinticSplinesEven::refitQSpline - warning DGESV returned INFO=",INFO
            end if

            dx = this%dx*this%dx/6

            do ii = 1,this%N
                this%y(2,ii) = dx*B(2*ii-1)
                this%y(3,ii) = 0.1d0*dx*dx*B(2*ii)
            end do


            return
        end subroutine refitQSpline

!         subroutine refitQSpline( this ,d2y0,d4y0 )
!     !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
!     !*      perform quintic spline fitting.
!     !*      currently do it the slow way with full LAPACK linear equation solver
!     !*              ie  solve A.B = C
!     !*      this could be improved dramatically with some thought.
!     !*
!     !*      B order (u1,v1,u2,v2...uN,vN)
!             type(QuinticSplineEven),intent(inout)       ::      this
!             real(kind=real64),intent(in),optional     ::      d2y0,d4y0       !   second,fourth derivs at x = xmin
!
!             real(kind=real64),dimension(2*this%N,2*this%N)       ::      A
!             real(kind=real64),dimension(2*this%N)                ::      B
!             integer,dimension(2*this%N)             ::      IPIV
!             integer                                 ::      NRHS,LDA,LDB,INFO
!
!             integer             ::      ii
!             real(kind=real64)                ::      dy1,dy2
!
!         !---    sanity check: arrays allocated?
!             if (this%N == 0) then
!                 write (unit=0,fmt='(a)') "Lib_QuinticSplinesEven::refit - warning refitting unallocated spline?"
!                 return
!             end if
!
!         !---    sanity check: x in ascending order
!             if ( this%xmin >= this%xmax) then
!                 write (unit=0,fmt='(a,i6,2g16.8)') "N,xmin,xmax=",this%N,this%xmin,this%xmax
!                 stop "Lib_QuinticSplinesEven::refit - error xmin >= xmax"
!             end if
!
!             NRHS = 1
!             LDA = 2*this%N
!             LDB = 2*this%N
!
!             A = 0.0d0
!             B = 0.0d0
!
!             A(1,1:2) = (/ 0.0d0,1.0d0 /)                                  !   sets y"" = 0.0d0
!             A(2,1:4) = (/ -6.0d0,0.0d0,6.0d0,-this%dx*this%dx /)  !   sets y"' = 0.0d0
!
!
!             if (present(d2y0)) then
!               A(1,1:2) = (/ 1.0d0,0.0d0 /)
!               B(1)   = d2y0
!             end if
!
!             if (present(d4y0)) then
!               A(2,1:4) = (/ 0.0d0,1.0d0,0.0d0,0.0d0 /)
!               B(2)   = d4y0
!             end if
!
!
!
!             do ii = 2,this%N-1
!
!             !---  fix first derivatives at x_i
!                 A( 2*ii-1, 2*ii-3:2*ii+2 ) = (/ 60*this%dx    ,-7*this%dx*this%dx*this%dx,            &
!                                                 240*this%dx   ,-16*this%dx*this%dx*this%dx,           &
!                                                 60*this%dx    ,-7*this%dx*this%dx*this%dx /)
!                 B( 2*ii-1 ) = 360*this%idx*( (this%y(1,ii+1) - 2*this%y(1,ii) + this%y(1,ii-1)) )
!
!           !---    fix third derivatives at x_i
!                 A( 2*ii  , 2*ii-3:2*ii+2 ) = (/ -6*this%idx   ,this%dx,                               &
!                                               12*this%idx ,4*this%dx,                             &
!                                               -6*this%idx ,this%dx /)
!               B( 2*ii ) = 0.0d0
!
!             end do
!             A(2*this%N-1,2*this%N-3:2*this%N)   &
!                                     = (/ -6.0d0,this%dx*this%dx,6.0d0,0.0d0 /)
!             A(2*this%N,2*this%N)    = 1
!
!
!         !---    for debugging only, here is the matrix and vector. I think they're OK.
!                       do ii = 1,2*this%N
!                           write (*,fmt='(1000f12.5)',advance="no") A(ii,:)
!                           write (*,fmt='(a,f12.5)') "            ",B(ii)
!                       end do
!         !---
!
!             call DGESV( 2*this%N, NRHS, A, LDA, IPIV, B, LDB, INFO )
!             if (INFO /= 0) then
!                 write (unit=0,fmt='(a,i6)') "Lib_QuinticSplinesEven::refit - warning DGESV returned INFO=",INFO
!             end if
!
!             do ii = 1,this%N
!                 this%y(2,ii) = B(2*ii-1)
!                 this%y(3,ii) = B(2*ii)
!             end do
!
!
!             return
!         end subroutine refit

        pure function splint(this,x) result(y)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      interpolate function to find value at x
            type(QuinticSplineEven),intent(in)                  ::      this
            real(kind=real64),intent(in)                        ::      x
            real(kind=real64)                                   ::      y

            integer                         ::      ii
            real(kind=real64)               ::      aa,bb,cc,dd,dx,idx2
            real(kind=real64),dimension(6)  ::      vv


            if (isNull(this)) then
                y = 0.0d0
                return
            end if
            idx2 = 6*this%idx*this%idx
            if (x < this%xmin) then

                vv(1:2) =  (/ this%y(1,1),this%y(1,2) /)
                vv(3:4) =  (/ this%y(2,1),this%y(2,2) /)
                vv(5:6) =  (/ this%y(3,1),this%y(3,2) /)

                aa = this%idx*( (-vv(1) + vv(2)) - (2*vv(3) + vv(4)) + (8*vv(5)+7*vv(6)) )
                bb = idx2*vv(3)
                cc = idx2*this%idx*( (-vv(3)+vv(4)) - 10*(2*vv(5)+vv(6)) )
                dd = 10*idx2*idx2*vv(5)

                dx = x - this%xmin
                y = this%y(1,1) + dx*(aa + 0.5d0*dx*(bb + 0.33333333333333d0*dx*(cc + dd*dx/4)))

                return
            else if (x >= this%xmax) then

                vv(1:2) =  (/ this%y(1,this%N-1),this%y(1,this%N) /)
                vv(3:4) =  (/ this%y(2,this%N-1),this%y(2,this%N) /)
                vv(5:6) =  (/ this%y(3,this%N-1),this%y(3,this%N) /)

                aa = this%idx*( (-vv(1)+vv(2)) + (vv(3)+2*vv(4)) - (7*vv(5)+8*vv(6)) )
                bb = idx2*vv(4)
                cc = idx2*this%idx*( (-vv(3)+vv(4)) + 10*(vv(5)+2*vv(6)) )
                dd = 10*idx2*idx2*vv(6)

                dx = x - this%xmax
                y = this%y(1,this%N) + dx*(aa + 0.5d0*dx*(bb + 0.33333333333333d0*dx*(cc + dd*dx/4)))

                return
            else
                aa = (x-this%xmin)*this%idx
                ii = int( aa ) + 1
                aa = ii - aa
                bb = 1.0d0-aa

            !---    explicit prefetch
            !   3/11/22 - gprof is complaining this is one of the most expensive lines in my MD. Try to rewrite?
!                vv(1:2) =  (/ this%y(1,ii),this%y(1,ii+1) /)
!                vv(3:4) =  (/ this%y(2,ii),this%y(2,ii+1) /)
!                vv(5:6) =  (/ this%y(3,ii),this%y(3,ii+1) /)
!                cc = aa*aa - 1.0d0 ; dd = bb*bb - 1.0d0
!
!
!                y  =        aa*( vv(1) + cc*(vv(3) + (3*cc-4.0d0)*vv(5)) )                          &
!                          + bb*( vv(2) + dd*(vv(4) + (3*dd-4.0d0)*vv(6)) )
!

!                vv(1:2) =  (/ this%y(1,ii),this%y(1,ii+1) /)
!                vv(3:4) =  (/ this%y(2,ii),this%y(2,ii+1) /)
!                vv(5:6) =  (/ this%y(3,ii),this%y(3,ii+1) /)
                cc = aa*aa - 1.0d0 ; dd = bb*bb - 1.0d0


                y =     aa*( this%y(1,ii)   + cc*(this%y(2,ii)   + (3*cc-4.0d0)*this%y(3,ii))   )                          
                y = y + bb*( this%y(1,ii+1) + dd*(this%y(2,ii+1) + (3*dd-4.0d0)*this%y(3,ii+1)) )


            end if

            return
        end function splint


        pure subroutine qsplint0a(this,x ,dy)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      interpolate function to find first derivative at x
            type(QuinticSplineEven),intent(in)              ::      this
            real(kind=real64),intent(in)                ::      x
            real(kind=real64),intent(out)               ::      dy

            integer                         ::      ii
            real(kind=real64)               ::      aa,bb,cc,dd,dx,idx2
            real(kind=real64),dimension(6)  ::      vv

            if (isNull(this)) then
                dy = 0.0d0
                return
            end if
            idx2 = 6*this%idx*this%idx
            if (x < this%xmin) then

                vv(1:2) =  (/ this%y(1,1),this%y(1,2) /)
                vv(3:4) =  (/ this%y(2,1),this%y(2,2) /)
                vv(5:6) =  (/ this%y(3,1),this%y(3,2) /)

                aa = this%idx*( (-vv(1) + vv(2)) - (2*vv(3) + vv(4)) + (8*vv(5)+7*vv(6)) )
                bb = idx2*vv(3)
                cc = idx2*this%idx*( (-vv(3)+vv(4)) - 10*(2*vv(5)+vv(6)) )
                dd = 10*idx2*idx2*vv(5)

                dx = x - this%xmin

                dy = aa + dx*(bb + 0.5d0*dx*( cc + 0.33333333333333d0*dx*dd ))
                return
            else if (x >= this%xmax) then

                vv(1:2) =  (/ this%y(1,this%N-1),this%y(1,this%N) /)
                vv(3:4) =  (/ this%y(2,this%N-1),this%y(2,this%N) /)
                vv(5:6) =  (/ this%y(3,this%N-1),this%y(3,this%N) /)

                aa = this%idx*( (-vv(1)+vv(2)) + (vv(3)+2*vv(4)) - (7*vv(5)+8*vv(6)) )
                bb = idx2*vv(4)
                cc = idx2*this%idx*( (-vv(3)+vv(4)) + 10*(vv(5)+2*vv(6)) )
                dd = 10*idx2*idx2*vv(6)

                dx = x - this%xmax
                dy = aa + dx*(bb + 0.5d0*dx*( cc + 0.33333333333333d0*dx*dd ))
                return
            else
                aa = (x-this%xmin)*this%idx
                ii = int( aa ) + 1
                aa = ii - aa
                bb = 1.0d0-aa

            !---    explicit prefetch
                vv(1:2) =  (/ this%y(1,ii),this%y(1,ii+1) /)
                vv(3:4) =  (/ this%y(2,ii),this%y(2,ii+1) /)
                vv(5:6) =  (/ this%y(3,ii),this%y(3,ii+1) /)
                cc = aa*aa - 1.0d0 ; dd = bb*bb - 1.0d0



                dy  = this%idx*( vv(2)-vv(1)                                                                    &
                               + (2.0d0+3*dd)*vv(4) - (2.0d0+3*cc)*vv(3)                                        &
                               + (15*dd*dd - 8.0d0)*vv(6) - (15*cc*cc - 8.0d0)*vv(5)  )

            end if


            return
        end subroutine qsplint0a


        pure subroutine qsplint0(this,x ,y,dy)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      interpolate function to find value and first derivative at x
            type(QuinticSplineEven),intent(in)              ::      this
            real(kind=real64),intent(in)                ::      x
            real(kind=real64),intent(out)               ::      y
            real(kind=real64),intent(out)               ::      dy

            integer                         ::      ii
            real(kind=real64)               ::      aa,bb,cc,dd,dx,idx2
            real(kind=real64),dimension(6)  ::      vv


            if (isNull(this)) then
                y = 0.0d0
                dy = 0.0d0
                return
            end if
            idx2 = 6*this%idx*this%idx
            if (x < this%xmin) then

                vv(1:2) =  (/ this%y(1,1),this%y(1,2) /)
                vv(3:4) =  (/ this%y(2,1),this%y(2,2) /)
                vv(5:6) =  (/ this%y(3,1),this%y(3,2) /)

                aa = this%idx*( (-vv(1) + vv(2)) - (2*vv(3) + vv(4)) + (8*vv(5)+7*vv(6)) )
                bb = idx2*vv(3)
                cc = idx2*this%idx*( (-vv(3)+vv(4)) - 10*(2*vv(5)+vv(6)) )
                dd = 10*idx2*idx2*vv(5)

                dx = x - this%xmin
                y = this%y(1,1) + dx*(aa + 0.5d0*dx*(bb + 0.33333333333333d0*dx*(cc + 0.25d0*dd*dx)))

                dy = aa + dx*(bb + 0.5d0*dx*( cc + 0.33333333333333d0*dx*dd ))

                return
            else if (x >= this%xmax) then

                vv(1:2) =  (/ this%y(1,this%N-1),this%y(1,this%N) /)
                vv(3:4) =  (/ this%y(2,this%N-1),this%y(2,this%N) /)
                vv(5:6) =  (/ this%y(3,this%N-1),this%y(3,this%N) /)

                aa = this%idx*( (-vv(1)+vv(2)) + (vv(3)+2*vv(4)) - (7*vv(5)+8*vv(6)) )
                bb = idx2*vv(4)
                cc = idx2*this%idx*( (-vv(3)+vv(4)) + 10*(vv(5)+2*vv(6)) )
                dd = 10*idx2*idx2*vv(6)

                dx = x - this%xmax
                y = this%y(1,this%N) + dx*(aa + 0.5d0*dx*(bb + 0.33333333333333d0*dx*(cc + 0.25d0*dd*dx)))
                dy = aa + dx*(bb + 0.5d0*dx*( cc + 0.33333333333333d0*dx*dd ))

                return
            else
                aa = (x-this%xmin)*this%idx
                ii = int( aa ) + 1
                aa = ii - aa
                bb = 1.0d0-aa

            !---    explicit prefetch
                vv(1:2) =  (/ this%y(1,ii),this%y(1,ii+1) /)
                vv(3:4) =  (/ this%y(2,ii),this%y(2,ii+1) /)
                vv(5:6) =  (/ this%y(3,ii),this%y(3,ii+1) /)
                cc = aa*aa - 1.0d0 ; dd = bb*bb - 1.0d0


                y  =        aa*( vv(1) + cc*(vv(3) + (3*cc-4.0d0)*vv(5)) )                          &
                          + bb*( vv(2) + dd*(vv(4) + (3*dd-4.0d0)*vv(6)) )



                dy  = this%idx*( vv(2)-vv(1)                                                                    &
                               + (2.0d0+3*dd)*vv(4) - (2.0d0+3*cc)*vv(3)                                        &
                               + (15*dd*dd - 8.0d0)*vv(6) - (15*cc*cc - 8.0d0)*vv(5)  )

            end if


            return
        end subroutine qsplint0


        pure subroutine qsplint1(this,x ,y,dy,d2y)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      interpolate function to find value and higher derivatives at x
            type(QuinticSplineEven),intent(in)          ::      this
            real(kind=real64),intent(in)            ::      x
            real(kind=real64),intent(out)           ::      y
            real(kind=real64),intent(out)           ::      dy
            real(kind=real64),intent(out)           ::      d2y

            integer                         ::      ii
            real(kind=real64)               ::      aa,bb,cc,dd,dx,idx2
            real(kind=real64),dimension(6)  ::      vv


            if (isNull(this))then
                y = 0.0d0
                dy = 0.0d0
                d2y = 0.0d0
                return
            end if
            idx2 = 6*this%idx*this%idx
            if (x < this%xmin) then

                vv(1:2) =  (/ this%y(1,1),this%y(1,2) /)
                vv(3:4) =  (/ this%y(2,1),this%y(2,2) /)
                vv(5:6) =  (/ this%y(3,1),this%y(3,2) /)

                aa = this%idx*( (-vv(1) + vv(2)) - (2*vv(3) + vv(4)) + (8*vv(5)+7*vv(6)) )
                bb = idx2*vv(3)
                cc = idx2*this%idx*( (-vv(3)+vv(4)) - 10*(2*vv(5)+vv(6)) )
                dd = 10*idx2*idx2*vv(5)

                dx = x - this%xmin
                y = this%y(1,1) + dx*(aa + 0.5d0*dx*(bb + 0.33333333333333d0*dx*(cc + 0.25d0*dd*dx)))

                dy = aa + dx*(bb + 0.5d0*dx*( cc + 0.33333333333333d0*dx*dd ))
                d2y = bb + dx*(cc + 0.5d0*dx*dd)

                return
            else if (x >= this%xmax) then

                vv(1:2) =  (/ this%y(1,this%N-1),this%y(1,this%N) /)
                vv(3:4) =  (/ this%y(2,this%N-1),this%y(2,this%N) /)
                vv(5:6) =  (/ this%y(3,this%N-1),this%y(3,this%N) /)

                aa = this%idx*( (-vv(1)+vv(2)) + (vv(3)+2*vv(4)) - (7*vv(5)+8*vv(6)) )
                bb = idx2*vv(4)
                cc = idx2*this%idx*( (-vv(3)+vv(4)) + 10*(vv(5)+2*vv(6)) )
                dd = 10*idx2*idx2*vv(6)

                dx = x - this%xmax
                y = this%y(1,this%N) + dx*(aa + 0.5d0*dx*(bb + 0.33333333333333d0*dx*(cc + 0.25d0*dd*dx)))
                dy = aa + dx*(bb + 0.5d0*dx*( cc + 0.33333333333333d0*dx*dd ))
                d2y = bb + dx*(cc + 0.5d0*dx*dd)

                return
            else
                aa = (x-this%xmin)*this%idx
                ii = int( aa ) + 1
                aa = ii - aa
                bb = 1.0d0-aa

            !---    explicit prefetch
                vv(1:2) =  (/ this%y(1,ii),this%y(1,ii+1) /)
                vv(3:4) =  (/ this%y(2,ii),this%y(2,ii+1) /)
                vv(5:6) =  (/ this%y(3,ii),this%y(3,ii+1) /)
                cc = aa*aa - 1.0d0 ; dd = bb*bb - 1.0d0


                y  =        aa*( vv(1) + cc*(vv(3) + (3*cc-4.0d0)*vv(5)) )                          &
                          + bb*( vv(2) + dd*(vv(4) + (3*dd-4.0d0)*vv(6)) )



                dy  = this%idx*( vv(2)-vv(1)                                                                    &
                               + (2.0d0+3*dd)*vv(4) - (2.0d0+3*cc)*vv(3)                                        &
                               + (15*dd*dd - 8.0d0)*vv(6) - (15*cc*cc - 8.0d0)*vv(5)  )


                d2y =  idx2*( aa*vv(3) + bb*vv(4) + 10*( aa*cc*vv(5) + bb*dd*vv(6) ) )
!

               !print *,"qsplint1 ",x,y,dy,d2y  ,aa,bb,cc,dd,idx2 , "vv ",vv
            end if

            return
        end subroutine qsplint1


        pure subroutine qsplint2(this,x ,y,dy,d2y,d3y,d4y)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      interpolate function to find value and all derivatives at x
            type(QuinticSplineEven),intent(in)          ::      this
            real(kind=real64),intent(in)            ::      x
            real(kind=real64),intent(out)           ::      y
            real(kind=real64),intent(out)           ::      dy
            real(kind=real64),intent(out)           ::      d2y
            real(kind=real64),intent(out)           ::      d3y
            real(kind=real64),intent(out)           ::      d4y

            integer                         ::      ii
            real(kind=real64)               ::      aa,bb,cc,dd,dx,idx2
            real(kind=real64),dimension(6)  ::      vv


            if (isNull(this)) then
                y = 0.0d0
                dy = 0.0d0
                d2y = 0.0d0
                d3y = 0.0d0
                d4y = 0.0d0
                return
            end if
            idx2 = 6*this%idx*this%idx
            if (x < this%xmin) then

                vv(1:2) =  (/ this%y(1,1),this%y(1,2) /)
                vv(3:4) =  (/ this%y(2,1),this%y(2,2) /)
                vv(5:6) =  (/ this%y(3,1),this%y(3,2) /)

                aa = this%idx*( (-vv(1) + vv(2)) - (2*vv(3) + vv(4)) + (8*vv(5)+7*vv(6)) )
                bb = idx2*vv(3)
                cc = idx2*this%idx*( (-vv(3)+vv(4)) - 10*(2*vv(5)+vv(6)) )
                dd = 10*idx2*idx2*vv(5)

                dx = x - this%xmin
                y = this%y(1,1) + dx*(aa + 0.5d0*dx*(bb + 0.33333333333333d0*dx*(cc + 0.25d0*dd*dx)))

                dy = aa + dx*(bb + 0.5d0*dx*( cc + 0.33333333333333d0*dx*dd ))
                d2y = bb + dx*(cc + 0.5d0*dx*dd)
                d3y = cc + dx*dd
                d4y = dd
                return
            else if (x >= this%xmax) then

                vv(1:2) =  (/ this%y(1,this%N-1),this%y(1,this%N) /)
                vv(3:4) =  (/ this%y(2,this%N-1),this%y(2,this%N) /)
                vv(5:6) =  (/ this%y(3,this%N-1),this%y(3,this%N) /)

                aa = this%idx*( (-vv(1)+vv(2)) + (vv(3)+2*vv(4)) - (7*vv(5)+8*vv(6)) )
                bb = idx2*vv(4)
                cc = idx2*this%idx*( (-vv(3)+vv(4)) + 10*(vv(5)+2*vv(6)) )
                dd = 10*idx2*idx2*vv(6)

                dx = x - this%xmax
                y = this%y(1,this%N) + dx*(aa + 0.5d0*dx*(bb + 0.33333333333333d0*dx*(cc + 0.25d0*dd*dx)))
                dy = aa + dx*(bb + 0.5d0*dx*( cc + 0.33333333333333d0*dx*dd ))
                d2y = bb + dx*(cc + 0.5d0*dx*dd)
                d3y = cc + dx*dd
                d4y = dd
                return
            else
                aa = (x-this%xmin)*this%idx
                ii = int( aa ) + 1
                aa = ii - aa
                bb = 1.0d0-aa

            !---    explicit prefetch
                vv(1:2) =  (/ this%y(1,ii),this%y(1,ii+1) /)
                vv(3:4) =  (/ this%y(2,ii),this%y(2,ii+1) /)
                vv(5:6) =  (/ this%y(3,ii),this%y(3,ii+1) /)
                cc = aa*aa - 1.0d0 ; dd = bb*bb - 1.0d0


                y  =        aa*( vv(1) + cc*(vv(3) + (3*cc-4.0d0)*vv(5)) )                          &
                          + bb*( vv(2) + dd*(vv(4) + (3*dd-4.0d0)*vv(6)) )



                dy  = this%idx*( vv(2)-vv(1)                                                                    &
                               + (2.0d0+3*dd)*vv(4) - (2.0d0+3*cc)*vv(3)                                        &
                               + (15*dd*dd - 8.0d0)*vv(6) - (15*cc*cc - 8.0d0)*vv(5)  )


                d2y =  idx2*( aa*vv(3) + bb*vv(4)                             &
                             + 10*( aa*cc*vv(5) + bb*dd*vv(6) ) )
!


                d3y = idx2*this%idx*( vv(4) - vv(3)                             &
                                + 10*( (3*dd+2)*vv(6) - (3*cc+2)*vv(5)  ) )


                d4y = 10*idx2*idx2*( aa*vv(5) + bb*vv(6) )

            end if

            return
        end subroutine qsplint2

    !---


        pure subroutine qsplint_inrange0a(this,x ,dy)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      interpolate function to find first derivative at x
            type(QuinticSplineEven),intent(in)              ::      this
            real(kind=real64),intent(in)                ::      x
            real(kind=real64),intent(out)               ::      dy

            integer                         ::      ii
            real(kind=real64)               ::      aa,cc,dd 
           !real(kind=real64),dimension(6)  ::      vv

            if (isNull(this)) then
                dy = 0.0d0
                return
            end if
             
            if ( (this%xmin-x)*(this%xmax-x)>=0 ) then

                dy = 0.0d0
                return
            else
                aa = (x-this%xmin)*this%idx
                ii = ceiling(aa)  
                
                              
                dy = this%y(1,ii+1) - this%y(1,ii)
                cc = aa*aa - 1.0d0
                dy = dy - (2.0d0+3*cc)*this%y(2,ii) - (15*cc*cc - 8.0d0)*this%y(3,ii)
                
                dd = aa*(aa-2.0d0)
                dy = dy + (2.0d0+3*dd)*this%y(2,ii+1) + (15*dd*dd - 8.0d0)*this%y(3,ii+1)
                
                dy = dy * this%idx
                
            
               

!                 ii = int(aa) + 1
!                 aa = ii - aa
!                 bb = 1.0d0-aa
!                                               
!             !---    explicit prefetch
!                 vv(1:2) =  (/ this%y(1,ii),this%y(1,ii+1) /)
!                 vv(3:4) =  (/ this%y(2,ii),this%y(2,ii+1) /)
!                 vv(5:6) =  (/ this%y(3,ii),this%y(3,ii+1) /)
!                 cc = aa*aa - 1.0d0 ; dd = bb*bb - 1.0d0
! 
! 
!                 dy  = this%idx*( vv(2)-vv(1)                                                                    &
!                                + (2.0d0+3*dd)*vv(4) - (2.0d0+3*cc)*vv(3)                                        &
!                                + (15*dd*dd - 8.0d0)*vv(6) - (15*cc*cc - 8.0d0)*vv(5)  )
! 
                                               
                               
            end if


            return
        end subroutine qsplint_inrange0a


        pure subroutine qsplint_inrange0(this,x ,y,dy)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      interpolate function to find value and first derivative at x
            type(QuinticSplineEven),intent(in)              ::      this
            real(kind=real64),intent(in)                ::      x
            real(kind=real64),intent(out)               ::      y
            real(kind=real64),intent(out)               ::      dy

            integer                         ::      ii
            real(kind=real64)               ::      aa,bb,cc,dd,idx2
            real(kind=real64),dimension(6)  ::      vv


            if (isNull(this)) then
                y = 0.0d0
                dy = 0.0d0
                return
            end if
            idx2 = 6*this%idx*this%idx
            if ( (this%xmin-x)*(this%xmax-x)>=0 ) then
                y = 0.0d0
                dy = 0.0d0

                return
            else
                aa = (x-this%xmin)*this%idx
                ii = int( aa ) + 1
                aa = ii - aa
                bb = 1.0d0-aa

            !---    explicit prefetch
                vv(1:2) =  (/ this%y(1,ii),this%y(1,ii+1) /)
                vv(3:4) =  (/ this%y(2,ii),this%y(2,ii+1) /)
                vv(5:6) =  (/ this%y(3,ii),this%y(3,ii+1) /)
                cc = aa*aa - 1.0d0 ; dd = bb*bb - 1.0d0


                y  =        aa*( vv(1) + cc*(vv(3) + (3*cc-4.0d0)*vv(5)) )                          &
                          + bb*( vv(2) + dd*(vv(4) + (3*dd-4.0d0)*vv(6)) )



                dy  = this%idx*( vv(2)-vv(1)                                                                    &
                               + (2.0d0+3*dd)*vv(4) - (2.0d0+3*cc)*vv(3)                                        &
                               + (15*dd*dd - 8.0d0)*vv(6) - (15*cc*cc - 8.0d0)*vv(5)  )

            end if


            return
        end subroutine qsplint_inrange0


        pure subroutine qsplint_inrange1(this,x ,y,dy,d2y)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      interpolate function to find value and higher derivatives at x
            type(QuinticSplineEven),intent(in)          ::      this
            real(kind=real64),intent(in)            ::      x
            real(kind=real64),intent(out)           ::      y
            real(kind=real64),intent(out)           ::      dy
            real(kind=real64),intent(out)           ::      d2y

            integer                         ::      ii
            real(kind=real64)               ::      aa,bb,cc,dd,idx2
            real(kind=real64),dimension(6)  ::      vv


            if (isNull(this)) then
                y = 0.0d0
                dy = 0.0d0
                d2y = 0.0d0
                return
            end if
            idx2 = 6*this%idx*this%idx
            if ( (this%xmin-x)*(this%xmax-x)>=0 ) then

                y = 0.0d0
                dy = 0.0d0
                d2y = 0.0d0
                return

            else
                aa = (x-this%xmin)*this%idx
                ii = int( aa ) + 1
                aa = ii - aa
                bb = 1.0d0-aa

            !---    explicit prefetch
                vv(1:2) =  (/ this%y(1,ii),this%y(1,ii+1) /)
                vv(3:4) =  (/ this%y(2,ii),this%y(2,ii+1) /)
                vv(5:6) =  (/ this%y(3,ii),this%y(3,ii+1) /)
                cc = aa*aa - 1.0d0 ; dd = bb*bb - 1.0d0


!                y  =        aa*( vv(1) + cc*(vv(3) + (3*cc-4.0d0)*vv(5)) )                          &
!                          + bb*( vv(2) + dd*(vv(4) + (3*dd-4.0d0)*vv(6)) )
!
!
!
!                dy  = this%idx*( vv(2)-vv(1)                                                                    &
!                               + (2.0d0+3*dd)*vv(4) - (2.0d0+3*cc)*vv(3)                                        &
!                               + (15*dd*dd - 8.0d0)*vv(6) - (15*cc*cc - 8.0d0)*vv(5)  )
!
!
!                d2y =  idx2*( aa*vv(3) + bb*vv(4)                             &
!                             + 10*( aa*cc*vv(5) + bb*dd*vv(6) ) )
!!
!


                y = aa*(vv(1) + (aa*aa-1)*(vv(3) + (3*aa*aa-7)*vv(5)))          &
                  + bb*(vv(2) + (bb*bb-1)*(vv(4) + (3*bb*bb-7)*vv(6)))

                dy = this%idx*( (vv(2)-vv(1)) - (vv(4)-vv(3)) + 7*(vv(6)-vv(5)) &
                                - 3*aa*aa*(vv(3) + 5*(aa*aa-2)*vv(5))           &
                                + 3*bb*bb*(vv(4) + 5*(bb*bb-2)*vv(6)) )            
                                
                d2y = 6*this%idx*this%idx*( aa*(vv(3) + 10*(aa*aa-1)*vv(5))     &
                                          + bb*(vv(4) + 10*(bb*bb-1)*vv(6)) )
                  
                  
!
            end if

            return
        end subroutine qsplint_inrange1


        pure subroutine qsplint_inrange2(this,x ,y,dy,d2y,d3y,d4y)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      interpolate function to find value and all derivatives at x
            type(QuinticSplineEven),intent(in)          ::      this
            real(kind=real64),intent(in)            ::      x
            real(kind=real64),intent(out)           ::      y
            real(kind=real64),intent(out)           ::      dy
            real(kind=real64),intent(out)           ::      d2y
            real(kind=real64),intent(out)           ::      d3y
            real(kind=real64),intent(out),optional           ::      d4y

            integer                         ::      ii
            real(kind=real64)               ::      aa,bb,cc,dd,idx2
            real(kind=real64),dimension(6)  ::      vv


            if (isNull(this)) then
                y = 0.0d0
                dy = 0.0d0
                d2y = 0.0d0
                d3y = 0.0d0
                if (present(d4y)) d4y = 0.0d0
                return
            end if
            idx2 = 6*this%idx*this%idx
            if ( (this%xmin-x)*(this%xmax-x)>=0 ) then

                y = 0.0d0
                dy = 0.0d0
                d2y = 0.0d0
                d3y = 0.0d0
                if (present(d4y)) d4y = 0.0d0
                return
            else
                aa = (x-this%xmin)*this%idx
                ii = int( aa ) + 1
                aa = ii - aa
                bb = 1.0d0-aa

            !---    explicit prefetch
                vv(1:2) =  (/ this%y(1,ii),this%y(1,ii+1) /)
                vv(3:4) =  (/ this%y(2,ii),this%y(2,ii+1) /)
                vv(5:6) =  (/ this%y(3,ii),this%y(3,ii+1) /)
                cc = aa*aa - 1.0d0 ; dd = bb*bb - 1.0d0


                y  =        aa*( vv(1) + cc*(vv(3) + (3*cc-4.0d0)*vv(5)) )                          &
                          + bb*( vv(2) + dd*(vv(4) + (3*dd-4.0d0)*vv(6)) )



                dy  = this%idx*( vv(2)-vv(1)                                                                    &
                               + (2.0d0+3*dd)*vv(4) - (2.0d0+3*cc)*vv(3)                                        &
                               + (15*dd*dd - 8.0d0)*vv(6) - (15*cc*cc - 8.0d0)*vv(5)  )


                d2y =  idx2*( aa*vv(3) + bb*vv(4)                             &
                             + 10*( aa*cc*vv(5) + bb*dd*vv(6) ) )
!


                d3y = idx2*this%idx*( vv(4) - vv(3)                             &
                                + 10*( (3*dd+2)*vv(6) - (3*cc+2)*vv(5)  ) )


                if (present(d4y)) d4y = 10*idx2*idx2*( aa*vv(5) + bb*vv(6) )

            end if

            return
        end subroutine qsplint_inrange2



    !---

        pure subroutine qsplint_array_inrange0(this,n,x ,y,dy)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      interpolate function to find value and first derivative at array of positions x
    !*      where all x are guaranteed in range
            type(QuinticSplineEven),intent(in)              ::      this
            integer,intent(in)                              ::      n
            real(kind=real64),dimension(:),intent(in)       ::      x
            real(kind=real64),dimension(:),intent(out)      ::      y
            real(kind=real64),dimension(:),intent(out)      ::      dy

            integer             ::      ii,jj
            real(kind=real64)                ::      aa,bb,cc,dd

            real(kind=real64),dimension(6)  ::  vv

            y = 0.0d0
            dy = 0.0d0
            do jj = 1,n

                if ( (this%xmin-x(jj))*(this%xmax-x(jj))>0 ) cycle

                aa = (x(jj)-this%xmin)*this%idx
                ii = int( aa ) + 1
                aa = ii - aa
                bb = 1.0d0-aa


            !---    explicit prefetch
                vv(1:2) = (/ this%y(1,ii),this%y(1,ii+1) /)
                vv(3:4) = (/ this%y(2,ii),this%y(2,ii+1) /)
                vv(5:6) = (/ this%y(3,ii),this%y(3,ii+1) /)
                cc = aa*aa - 1.0d0 ; dd = bb*bb - 1.0d0



                y(jj)   =  aa*( vv(1) + cc*(vv(3) + (3*cc-4.0d0)*vv(5)) )                           &
                         + bb*( vv(2) + dd*(vv(4) + (3*dd-4.0d0)*vv(6)) )


                dy(jj)  = this%idx*( vv(2)-vv(1)                                                                    &
                               + (2.0d0+3*dd)*vv(4) - (2.0d0+3*cc)*vv(3)                                        &
                               + (15*dd*dd - 8.0d0)*vv(6) - (15*cc*cc - 8.0d0)*vv(5)  )



            end do

            return
        end subroutine qsplint_array_inrange0

        pure subroutine qsplint_array_inrange0a(this,n,x ,dy,derivative)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      interpolate function to find first derivative at array of positions x
    !*      where all x are guaranteed in range
            type(QuinticSplineEven),intent(in)              ::      this
            integer,intent(in)                              ::      n
            real(kind=real64),dimension(:),intent(in)       ::      x
            real(kind=real64),dimension(:),intent(out)      ::      dy
            logical,intent(in)                              ::      derivative

            integer                             ::      ii,jj
            real(kind=real64)                   ::      aa,bb 

            real(kind=real64),dimension(6)  ::  vv
            
            ! real(kind=real64)                   ::      v1,v2,v3,v4,v5,v6

            if (.not. derivative) call qsplint_array_inrange0b(this,n,x ,dy)

            !idx2 = 6*this%idx*this%idx
            dy = 0.0d0
            do jj = 1,n

                if ( (this%xmin-x(jj))*(this%xmax-x(jj))>0 ) cycle


                aa = (x(jj)-this%xmin)*this%idx
               ! ii = int( aa ) + 1
                ii = ceiling(aa)
                aa = ii - aa
                bb = 1.0d0-aa

 
            !---    explicit prefetch
                vv(1:2) = (/ this%y(1,ii),this%y(1,ii+1) /)
                vv(3:4) = (/ this%y(2,ii),this%y(2,ii+1) /)
                vv(5:6) = (/ this%y(3,ii),this%y(3,ii+1) /)
                !cc = aa*aa - 1.0d0 ; dd = bb*bb - 1.0d0
                dy(jj)  = this%idx*( (vv(2)-vv(1)) - (vv(4)-vv(3)) + 7*(vv(6)-vv(5))        & 
                                     - 3*aa*aa*( vv(3) + 5*(aa*aa-2)*vv(5) )                &
                                     + 3*bb*bb*( vv(4) + 5*(bb*bb-2)*vv(6) )  )        
 
!                dy(jj)  = this%idx*( vv(2)-vv(1)                                                                    &
!                               + (2.0d0+3*dd)*vv(4) - (2.0d0+3*cc)*vv(3)                                        &
!                               + (15*dd*dd - 8.0d0)*vv(6) - (15*cc*cc - 8.0d0)*vv(5)  )
 
                

                
            end do

            return
        end subroutine qsplint_array_inrange0a


        pure subroutine qsplint_array_inrange0b(this,n,x ,y)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      interpolate function to find value at array of positions x
    !*      where all x are guaranteed in range
            type(QuinticSplineEven),intent(in)              ::      this
            integer,intent(in)                              ::      n
            real(kind=real64),dimension(:),intent(in)       ::      x
            real(kind=real64),dimension(:),intent(out)      ::      y

            integer             ::      ii,jj
            real(kind=real64)                ::      aa,bb,cc,dd

            real(kind=real64),dimension(6)  ::  vv

            y = 0.0d0
            do jj = 1,n

                if ( (this%xmin-x(jj))*(this%xmax-x(jj))>0 ) cycle


                aa = (x(jj)-this%xmin)*this%idx
                ii = int( aa ) + 1
                aa = ii - aa
                bb = 1.0d0-aa


            !---    explicit prefetch
                vv(1:2) = (/ this%y(1,ii),this%y(1,ii+1) /)
                vv(3:4) = (/ this%y(2,ii),this%y(2,ii+1) /)
                vv(5:6) = (/ this%y(3,ii),this%y(3,ii+1) /)
                cc = aa*aa - 1.0d0 ; dd = bb*bb - 1.0d0

                y(jj)   =   aa*( vv(1) + cc*(vv(3) + (3*cc-4.0d0)*vv(5)) )                          &
                          + bb*( vv(2) + dd*(vv(4) + (3*dd-4.0d0)*vv(6)) )



            end do

            return
        end subroutine qsplint_array_inrange0b



    !---



        pure subroutine qsplint_array_inrange_sum0(this,n,x ,y)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      interpolate function to find value at array of positions x
    !*      where all x are guaranteed in range
    !*      return the sum of the answer
            type(QuinticSplineEven),intent(in)              ::      this
            integer,intent(in)                              ::      n
            real(kind=real64),dimension(:),intent(in)       ::      x
            real(kind=real64),intent(out)                   ::      y

            integer                         ::      ii,jj
            real(kind=real64)               ::      aa,bb,cc,dd


!             real(kind=real64),dimension(6)  ::      vv


            y = 0.0d0
            do jj = 1,n

                if ( (this%xmin-x(jj))*(this%xmax-x(jj))>0 ) cycle


                aa = (x(jj)-this%xmin)*this%idx
                ii = int( aa ) + 1


                aa = ii - aa ; cc = aa*aa - 1.0d0





                y   =  y  + aa*( this%y(1,ii) + cc*(this%y(2,ii) + (3*cc-4.0d0)*this%y(3,ii)) )

                bb = 1.0d0-aa ; dd = aa*(aa-2.0d0)      !   dd = bb*bb-1

                y   =  y  + bb*( this%y(1,ii+1) + dd*(this%y(2,ii+1) + (3*dd-4.0d0)*this%y(3,ii+1)) )


!             !---    explicit prefetch
!                 vv(1:2) = (/ this%y(1,ii),this%y(1,ii+1) /)
!                 vv(3:4) = (/ this%y(2,ii),this%y(2,ii+1) /)
!                 vv(5:6) = (/ this%y(3,ii),this%y(3,ii+1) /)
!
!                 aa = ii - aa
!                 bb = 1.0d0-aa
!
!
!                 cc = aa*aa - 1.0d0 ; dd = bb*bb - 1.0d0
!
!
!                 y   =  y  + aa*( vv(1) + cc*(vv(3) + (3*cc-4.0d0)*vv(5)) )                          &
!                           + bb*( vv(2) + dd*(vv(4) + (3*dd-4.0d0)*vv(6)) )

            end do

            return
        end subroutine qsplint_array_inrange_sum0

        pure subroutine qsplint_array_inrange_sum1(this,n,x ,y,dy)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      interpolate function to find value at array of positions x
    !*      where all x are guaranteed in range
    !*      return the sum of the answer and derivative at each point
            type(QuinticSplineEven),intent(in)              ::      this
            integer,intent(in)                              ::      n
            real(kind=real64),dimension(:),intent(in)       ::      x
            real(kind=real64),intent(out)                   ::      y
            real(kind=real64),dimension(:),intent(out)      ::      dy


            integer                         ::      ii,jj
            real(kind=real64)               ::      aa,bb,cc,dd
            real(kind=real64),dimension(6)  ::      vv
            

            y = 0.0d0
            dy = 0.0d0
            do jj = 1,n

                if ( (this%xmin-x(jj))*(this%xmax-x(jj))>0 ) cycle


                aa = (x(jj)-this%xmin)*this%idx
                !ii = int( aa ) + 1
                ii = ceiling(aa)
                aa = ii - aa
                bb = 1.0d0-aa

                     
            !---    explicit prefetch         
                vv(1:2) = (/ this%y(1,ii),this%y(1,ii+1) /)
                vv(3:4) = (/ this%y(2,ii),this%y(2,ii+1) /)
                vv(5:6) = (/ this%y(3,ii),this%y(3,ii+1) /)
                
                
            !    y = y + aa*(vv(1) + (aa*aa-1)*(vv(3) + (3*aa*aa-7)*vv(5)))   &
            !          + bb*(vv(2) + (bb*bb-1)*(vv(4) + (3*bb*bb-7)*vv(6)))    
            !    
            !    
            !    dy(jj)  = this%idx*( (vv(2)-vv(1)) - (vv(4)-vv(3)) + 7*(vv(6)-vv(5))        & 
            !         - 3*aa*aa*( vv(3) + 5*(aa*aa-2)*vv(5) )                &
            !         + 3*bb*bb*( vv(4) + 5*(bb*bb-2)*vv(6) )  )        

               
                
                 cc = aa*aa - 1.0d0 ; dd = bb*bb - 1.0d0
              
              
              
                 y      =  y  + aa*( vv(1) + cc*(vv(3) + (3*cc-4.0d0)*vv(5)) )                           &
                              + bb*( vv(2) + dd*(vv(4) + (3*dd-4.0d0)*vv(6)) )
              
                 dy(jj) = this%idx*( vv(2)-vv(1)                                                                     &
                                + (2.0d0+3*dd)*vv(4) - (2.0d0+3*cc)*vv(3)                                        &
                                + (15*dd*dd - 8.0d0)*vv(6) - (15*cc*cc - 8.0d0)*vv(5)  )

 
            end do

            return
        end subroutine qsplint_array_inrange_sum1

!-------

        subroutine delete0(this)
    !---^^^^^^^^^^^^^^^^^^^^^^^^ 
    !*      deallocate dynamically allocated memory.
            type(QuinticSplineEven),intent(inout)   ::      this
!            logical,intent(in),optional             ::      scrub
            if (isNull(this)) return
            deallocate(this%y)
!             deallocate(this%u)
!             deallocate(this%v)
            nullify(this%y)
!             nullify(this%u)
!             nullify(this%v)
            this%N = 0
            return
        end subroutine delete0
        
        
        
        subroutine clone0(this,that)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^    
    !*      deep copy with allocation this = that
            type(QuinticSplineEven),intent(inout)   ::      this
            type(QuinticSplineEven),intent(in)      ::      that
            
            if (that%N == 0) then
                call delete(this)
                return
            end if
            
            if (this%N /= 0) then
                if (this%N < that%N) deallocate(this%y)
            end if
            
            this%N = that%N
            this%xmin = that%xmin
            this%xmax = that%xmax
            this%dx = that%dx
            this%idx = that%idx
            this%key = that%key
            this%value = that%value
            
            allocate(this%y(3,this%N))
            this%y(1:3,1:this%N) = that%y(1:3,1:this%N)
            
            return
        end subroutine clone0
        
        
         

        
        subroutine clear0(this)
    !---^^^^^^^^^^^^^^^^^^^^^^^
    !*      zeroes memory without deallocation
            type(QuinticSplineEven),intent(inout)   ::      this
            if (isNull(this)) return
            this%y = 0
            return
        end subroutine clear0
        
        
        subroutine report0( this,u,o )
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
            type(QuinticSplineEven),intent(in)      ::      this
            integer,intent(in),optional             ::      u,o
            integer             ::      uu,oo
            integer         ::      ii,jj
            uu = 6 ; if (present(u)) uu = u
            oo = 0 ; if (present(o)) oo = o
            if (isNull(this)) then
                write (unit=uu,fmt='(a,i6)') repeat(" ",oo)//"QuinticSplineEven(null)"
                return
            else
                write (unit=uu,fmt='(a,i6)') repeat(" ",oo)//"QuinticSplineEven N = ",this%N
            end if
            write (unit=uu,fmt='(a,a8,2a16)') repeat(" ",oo+2),"deriv ","left","right"
            do ii = 1,4
                write (unit=uu,fmt='(a,i8)',advance="no") repeat(" ",oo+2),ii
                do jj = 1,2
                    if (this%key(jj,ii)==QS_DERIV_FREE) then
                        write (unit=uu,fmt='(a16)',advance="no") "        free"
                    else if (this%key(jj,ii)==QS_DERIV_FIXED) then
                        write (unit=uu,fmt='(g16.8)',advance="no") this%value(jj,ii)
                    else
                        write (unit=uu,fmt='(a14,i2)',advance="no") "        interp",this%key(jj,ii)
                    end if
                end do
                write (unit=uu,fmt='(a)',advance="yes") ""
            end do


            write (unit=uu,fmt='(a,a6,4a16)') repeat(" ",oo),"knot","x","y","y""","y"""""
            do ii = 1,this%N
                write (unit=uu,fmt='(a,i6,4f16.8)') repeat(" ",oo),ii,(this%xmin+(ii-1)*this%dx),this%y(1,ii),this%y(2,ii),this%y(3,ii)
            end do
            return
        end subroutine report0

    !---

        pure function getXmin(this) result(xmin)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
            type(QuinticSplineEven),intent(in)      ::      this
            real(kind=real64)                       ::      xmin
            xmin = this%xmin
            return
        end function getXmin

        pure function getXmax(this) result(xmax)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
            type(QuinticSplineEven),intent(in)      ::      this
            real(kind=real64)                       ::      xmax
            xmax = this%xmax
            return
        end function getXmax


        subroutine setXmin(this,xmin,refit)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
            type(QuinticSplineEven),intent(inout)   ::      this
            real(kind=real64),intent(in)            ::      xmin
            logical,intent(in)                      ::      refit
            this%xmin = xmin
            this%dx = (this%xmax - this%xmin)/(this%N-1)
            this%idx = (this%N-1)/(this%xmax - this%xmin)

            if (refit) call refitQSpline(this)
            return
        end subroutine setXmin

        subroutine setXmax(this,xmax,refit)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
            type(QuinticSplineEven),intent(inout)   ::      this
            real(kind=real64),intent(in)            ::      xmax
            logical,intent(in)                      ::      refit
            this%xmax = xmax
            this%dx = (this%xmax - this%xmin)/(this%N-1)
            this%idx = (this%N-1)/(this%xmax - this%xmin)

            if (refit) call refitQSpline(this)
            return
        end subroutine setXmax

    !---

        pure function getnKnots(this) result(nKnots)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
            type(QuinticSplineEven),intent(in)      ::      this
            integer                                 ::      nKnots
            nKnots = this%N
            return
        end function getnKnots

        pure function isNull0(this) result(is)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
            type(QuinticSplineEven),intent(in)      ::      this
            logical                                 ::      is
            is = (this%N == 0)
            return
        end function isNull0

    !---


        subroutine gety0(this,y)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^
            type(QuinticSplineEven),intent(in)          ::      this
            real(kind=real64),dimension(this%N),intent(out) ::      y
            y(:) = this%y(1,:)
            return
        end subroutine gety0

        subroutine sety0(this,y,refit)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^
            type(QuinticSplineEven),intent(inout)       ::      this
            real(kind=real64),dimension(this%N),intent(in)  ::      y
            logical,intent(in)                      ::      refit
            this%y(1,:) = y(:)
            if (refit) call refitQSpline(this)
            return
        end subroutine sety0

    !---


        subroutine inputFromXML0(xml,this,ok)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
            type(QuinticSplineEven),intent(out)     ::      this
            type(NBAX),intent(in)                   ::      xml
            logical,intent(out)                     ::      ok


            logical                                         ::      allok
            real(kind=real64),dimension(:),allocatable      ::      y_in
            integer                                         ::      ii,nn,jj,kk,ll,mm
            real(kind=real64)                               ::      xmin,xmax,xx
            character(len=2048)                             ::      line
            real(kind=real64),dimension(256)                ::      ra
            type(NBAX),pointer                              ::      xmlp
            logical                                         ::      ok2


            integer,dimension(2,4)          ::      key
            real(kind=real64),dimension(8)  ::      value

            call getAttributeValue(xml,"N",nn,ok) ; allok = ok
            if (.not. ok) write(unit=0,fmt='(a)') "QuinticSplinesEven::inputFromXML0 error - could not read attribute ""N"""
            xmin = 0.0d0 ; call getAttributeValue(xml,"xmin",xmin,ok)
            if ((nn>0) .and. .not. ok) stop "QuinticSplinesEven::inputFromXML0 error - could not read attribute ""xmin"""
            xmax = 1.0d0 ; call getAttributeValue(xml,"xmax",xmax,ok)
            if ((nn>0) .and. .not. ok) stop "QuinticSplinesEven::inputFromXML0 error - could not read attribute ""xmax"""


        !---    read in nn reals from text lines
            allocate(y_in(nn))
            kk = 0
            do ii = 1,getNText(xml)
                call getText(xml,ii,line,ok)
                call parse( line,ra,jj )
                if (jj>0) then
                    y_in(kk+1:min(nn,kk+jj)) = ra(1:min(nn-kk,jj))
                    kk = kk + min(nn-kk,jj)
                    if (kk >= nn) exit
                end if
            end do
            if (kk < nn) stop "QuinticSplinesEven::inputFromXML0 error - unable to parse N reals from xml"

!             print *,"allok ",allok
            key = QS_DERIV_FREE
            kk = 0
            do ii = 1,4
                call getChild( xml,"Derivative",ii,xmlp,ok)
                if (.not. ok) then
!                     write(unit=0,fmt='(a,i6)') "QuinticSplinesEven::inputFromXML0 warning - could not read <Derivative/> ",ii
                else
                    call getAttributeValue(xmlp,"N",ll,ok) ; allok = allok .and. ok
                    if (.not. ok) write(unit=0,fmt='(a,i6)') "QuinticSplinesEven::inputFromXML0 error - could not read attribute ""N"" for <Derivative/> ",ii

!             print *,"allok ",allok

                    line = "" ; xx = 0.0d0
                    call getAttributeValue(xmlp,"left",xx,ok2)
                    call getAttributeValue(xmlp,"left",line,ok) ; allok = allok  .and. (ok .or. ok2)
                    if (.not. (ok .or. ok2)) write(unit=0,fmt='(a,i6)') "QuinticSplinesEven::inputFromXML0 warning - could not read attribute ""left"" for <Derivative/> ",ii
                    if (ok2) then
                        key(1,ll) = QS_DERIV_FIXED
                        kk = kk + 1
                        value(kk) = xx
                    else if ( (line(1:3)=="INT").or.(line(1:3)=="int").or.(line(1:3)=="Int") ) then
                        mm = index(line," ")
                        if (mm>0) then
                            call parse(line(mm:),mm,ok)
                        else
                            ok = .false.
                        end if
                        if (.not. ok) then
                            write(unit=0,fmt='(a,i6)') "QuinticSplinesEven::inputFromXML0 error - could not read attribute value ""left"" for <Derivative/> ",ii
                        else
                            key(1,ll) = mm
                        end if
                    end if
!             print *,"allok ",allok

                    line = "" ; xx = 0.0d0
                    call getAttributeValue(xmlp,"right",xx,ok2)
                    call getAttributeValue(xmlp,"right",line,ok) ; allok = allok .and. (ok .or. ok2)
                    if (.not. (ok .or. ok2)) write(unit=0,fmt='(a,i6)') "QuinticSplinesEven::inputFromXML0 warning - could not read attribute ""right"" for <Derivative/> ",ii
                    if (ok2) then
                        key(2,ll) = QS_DERIV_FIXED
                        kk = kk + 1
                        value(kk) = xx
                    else if ( (line(1:3)=="INT").or.(line(1:3)=="int").or.(line(1:3)=="Int") ) then
                        mm = index(line," ")
                        if (mm>0) then
                            call parse(line(mm:),mm,ok)
                        else
                            ok = .false.
                        end if
                        if (.not. ok) then
                            write(unit=0,fmt='(a,i6)') "QuinticSplinesEven::inputFromXML0 error - could not read attribute value ""right"" for <Derivative/> ",ii
                        else
                            key(2,ll) = mm
                        end if
                    end if

                end if
            end do


            this = QuinticSplineEven_ctor( nn,xmin,xmax , y_in , key,value)
            !call report(this)
            deallocate(y_in)

            ok = allok
            if (.not. ok) then
                stop "QuinticSplinesEven::inputFromXML0 error - unable to read QuinticSplineEven from xml"
            end if
            return
        end subroutine inputFromXML0


        subroutine outputAsXML0(this,xml)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
            type(QuinticSplineEven),intent(in)      ::      this
            type(NBAX),intent(out)                  ::      xml

            integer                                         ::      ii
            character(len=2048)                             ::      line
            type(NBAX),dimension(:),pointer                 ::      xmlp
            character(len=8)                                ::      key

            xml = NBAX_ctor("QuinticSplineEven")
            call addAttribute(xml,"N",this%N)
            if (isNull(this)) return
            call addAttribute(xml,"xmin",this%xmin)
            call addAttribute(xml,"xmax",this%xmax)


            allocate(xmlp(4))
            do ii = 1,4
                xmlp(ii) = NBAX_ctor("Derivative")
                call addAttribute(xmlp(ii),"N",ii)
                select case(this%key(1,ii))
                    case (QS_DERIV_FREE)
                        call addAttribute(xmlp(ii),"left","free")
                    case (QS_DERIV_FIXED)
                        call addAttribute(xmlp(ii),"left",this%value(1,ii))
                    case default
                        write(key,fmt='(a6,i2)') "interp",this%key(1,ii)
                        call addAttribute(xmlp(ii),"left",key)
                end select
                select case(this%key(2,ii))
                    case (QS_DERIV_FREE)
                        call addAttribute(xmlp(ii),"right","free")
                    case (QS_DERIV_FIXED)
                        call addAttribute(xmlp(ii),"right",this%value(2,ii))
                    case default
                        write(key,fmt='(a6,i2)') "interp",this%key(2,ii)
                        call addAttribute(xmlp(ii),"right",key)
                end select
                call addChild(xml,xmlp(ii))
            end do



            do ii = 1,this%N,8
                write(line,fmt='(8g24.16)') this%y(1,ii:min(this%N,ii+7))
                call addText(xml,line)
            end do


            return
        end subroutine outputAsXML0

!         subroutine MPI_exchange0(this)
!     !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
!             type(QuinticSplineEven),intent(inout)       ::      this
!             integer,dimension(9)                ::  ia
!             real(kind=real64),dimension(10)     ::  ra
!             real(kind=real64),dimension(:),allocatable  ::  y
!
!
!             integer     ::      ii
!             integer     ::      nProcs,rank
!
!             call MPI_COMM_SIZE(MPI_COMM_WORLD,nProcs,ii)
!             if (nProcs == 1) return
!             call MPI_COMM_RANK(MPI_COMM_WORLD,rank,ii)
!
!             if (rank == 0) then
!                 ia(1) = this%N
!                 ia(2:9) = pack(this%key,.true.)
!                 call MPI_BCAST( ia,9,MPI_INTEGER,0,MPI_COMM_WORLD,ii)
!                 allocate(y(this%N))
!                 do ii = 1,this%N
!                     y(ii) = this%y(1,ii)
!                 end do
!             else
!                 call MPI_BCAST( ia,9,MPI_INTEGER,0,MPI_COMM_WORLD,ii)
!                 this%N = ia(1)
!                 this%key(:,:) = reshape( ia(2:9),(/2,4/) )
!                 allocate(y(this%N))
!             end if
!
!             ra(1:2) = (/ this%xmin,this%xmax /)
!             ra(3:10) = pack(this%value,.true.)
!             call MPI_BCAST( ra,10,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ii)
!
!             call MPI_BCAST( y,this%N,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ii)
!             this = QuinticSplineEven_ctor0( this%N,ra(1),ra(2),y  ,this%key,ra(3:10))
!             deallocate(y)
!             return
!         end subroutine MPI_exchange0

    end module Lib_QuinticSplinesEven

!
!     program testLib_QuinticSplinesEven
! !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
!         use Lib_QuinticSplinesEven
!         use Lib_Timers
!         use Lib_NBAX
!
!         use iso_Fortran_env
!         implicit none
!
!         integer,parameter           ::      N = 20
!
!         integer                 ::      ii
!         real(kind=real64),dimension(N)       ::      xx,yy
!         real(kind=real64)                    ::      xxxx,yyyy,dy,d2y,d3y,d4y
!         real(kind=real64),dimension(20)        ::     y_array,dy_array
!         real(kind=real64)                    ::      xmin,xmax
!
!         type(NBAX)                  ::      xml
!
!
!         integer                                     ::      Nx
!         real(kind=real64),dimension(:),allocatable  ::      xxx
!         type(Timer)                                 ::      tt
!         logical                                     ::      ok
!         integer,dimension(2,4)      ::      key
!         real(kind=real64),dimension(4)  ::  value
!         type(QuinticSplineEven)     ::      q
!
!         print *,""
!         print *,"QuinticSplineEven test program"
!         print *,""
!         do ii = 1,N
!            xxxx = (ii-1)*2.0d0*3.141592654d0/(N-1)
!            xx(ii) = xxxx
!            yy(ii) = cos( xxxx )
!         end do
!         xmin = 0.0d0
!         xmax = 2.0d0*3.141592654d0
!
!         q = QuinticSplineEven_ctor( N,xmin,xmax,yy  )
!         call report(q)
!         print *,""
!
!         key = QS_DERIV_FREE
!         key(:,2) = QS_DERIV_FIXED
!         key(:,4) = QS_DERIV_FIXED
!         value(:) = (/ -1.0d0,-1.0d0,1.0d0,1.0d0 /)
!         q = QuinticSplineEven_ctor( N,xmin,xmax,yy , key,value )
!         call report(q)
!         print *,""
!
!         call outputAsXML(q,xml)
!         open(unit=400,file="test.xml",action="write")
!             call output(xml,400,.true.)
!         close(unit=400)
!
!
!         write (*,fmt='(8a16)') "theta","cos(theta)","-sin(theta)","y","dy","d2y","d3y","d4y"
!
!         do ii = -10,N*10+11
!             xxxx = (ii-1)*2.0*3.141592654/(N*10)
!             yyyy = splint( q,xxxx )
!             call qsplint( q,xxxx, yyyy,dy,d2y,d3y,d4y )
!             write (*,fmt='(8f16.8)') xxxx,cos( xxxx ),-sin( xxxx ),yyyy,dy,d2y,d3y,d4y
!         end do
!
!
!         call delete(q)
!         xml = NBAX_ctor()
!         open(unit=500,file="test.xml",action="read")
!             call input(xml,500,ok)
!             print *,"input from disk ",ok
!             call inputFromXML( xml,q,ok )
!             print *,"inputFromXML ",ok
!         close(unit=500)
!         call report(q)
!         print *,""
!
!
!
! !         print *,"enter number of trials"
!         Nx = 100000000
!         tt = Timer_ctor()
!         allocate(xxx(Nx))
!         call random_number(xxx)
!         print *,"elapsed time (number gen) ",elapsed(tt)
!
!         call reset(tt)
!         d2y = 0.0d0     !   sum answers to avoid optimizing out...
! !         tt = Timer_ctor()
!         do ii = 1,Nx
!             yyyy = splint( q,xxx(ii) )
!             d2y = d2y + yyyy
!         end do
!         d3y = elapsed(tt)
!         print *,"elapsed time (val) ",d3y/Nx
!         print *,"(avoid opt val sum)",d2y/Nx
!
!         call reset(tt)
!         d2y = 0.0d0     !   sum answers to avoid optimizing out...
! !         tt = Timer_ctor()
!         do ii = 1,Nx
!             call qsplint( q,xxx(ii), dy )
!             d2y = d2y + dy
!         end do
!         d3y = elapsed(tt)
!         print *,"elapsed time (deriv) ",d3y/Nx
!         print *,"(avoid opt deriv sum)",d2y/Nx
!
!
!         call reset(tt)
!         d2y = 0.0d0     !   sum answers to avoid optimizing out...
! !         tt = Timer_ctor()
!         do ii = 1,Nx,20
!             call qsplint_array_inrange( q,xxx(ii:ii+19), y_array,dy_array )
!             d2y = d2y + sum(y_array)
!         end do
!         d3y = elapsed(tt)
!         print *,"elapsed time (val + deriv array) ",d3y/Nx
!         print *,"(avoid opt val sum)",d2y/Nx
!
!         call reset(tt)
!         d2y = 0.0d0     !   sum answers to avoid optimizing out...
! !         tt = Timer_ctor()
!         do ii = 1,Nx,20
!             call qsplint_array_inrange( q,xxx(ii:ii+19), dy_array , derivative = .true. )
!             d2y = d2y + sum(dy_array)
!         end do
!         d3y = elapsed(tt)
!         print *,"elapsed time (deriv array) ",d3y/Nx
!         print *,"(avoid opt deriv sum)",d2y/Nx
!
!
!
!
!
!         call reset(tt)
!         d2y = 0.0d0     !   sum answers to avoid optimizing out...
! !         tt = Timer_ctor()
!         do ii = 1,Nx
!             call qsplint( q,xxx(ii), yyyy,dy )
!             d2y = d2y + yyyy
!         end do
!         d3y = elapsed(tt)
!         print *,"elapsed time (val + deriv) ",d3y/Nx
!         print *,"(avoid opt val sum)",d2y/Nx
!
!
!         call delete(q)
!         print *,""
!         print *,"done"
!         print *,""
!
!     end program testLib_QuinticSplinesEven
