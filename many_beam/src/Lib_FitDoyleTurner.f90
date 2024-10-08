
    module Lib_FitDoyleTurner
!---^^^^^^^^^^^^^^^^^^^^^^^^^
!*
!*      This module fits a scattering function f(s)
!*      to a Doyle-Turner approximation
!*          f_DT(s) = sum_j a_j Exp( -b)j s^2 )
!*      
!*      This is acheived by minimising in a least squares sense 
!*          fit = sum_i ( f(s_i) - f_DT(s_i) )^2
!*      where f(s_i) is a tabulated set of data points for the scattering function
!*
!*      The fit is non-linear, so is solved using steepest descent (Golden Section search)
!*

        use Lib_GoldenSection
        use iso_fortran_env
        implicit none
        private

        real(kind=real64),public        ::      Lib_FitDoyleTurner_DERIV_TOL = 1.0d-16 
        real(kind=real64),public        ::      Lib_FitDoyleTurner_FIT_TOL = 1.0d-16
        integer,public                  ::      Lib_FitDoyleTurner_N_GOLD_DIRECTIONS = 100
        real(kind=real64),public        ::      Lib_FitDoyleTurner_MINA = -5.0d0               !   minimum value of a DT coeff starting point in search
        real(kind=real64),public        ::      Lib_FitDoyleTurner_MAXA = +5.0d0               !   minimum value of a DT coeff starting point in search

        public          ::      getFit
        public          ::      fitDoyleTurner
        public          ::      DoyleTurnerSum

        interface       fitDoyleTurner
            module procedure        fitDoyleTurner0             !   fit with known starting point
            module procedure        fitDoyleTurner1             !   random search starting points
        end interface
            

        
    contains
!---^^^^^^^^

        subroutine fitDoyleTurner1( s,f,n,a,b, fit )       
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      given the tabulated scattering function f_i at points s_i
    !*      try different starting points for Doyle Turner coefficients, and find best solution
            real(kind=real64),dimension(:),intent(in)           ::      s       !   (1:nData) scattering points
            real(kind=real64),dimension(:),intent(in)           ::      f       !   (1:nData) scattering amplitudes
            real(kind=real64),dimension(:),intent(out)          ::      a       !   Doyle Turner a coefficients  
            real(kind=real64),dimension(:),intent(out)          ::      b       !   Doyle Turner b coefficients
            integer,intent(in)                                  ::      n
            real(kind=real64),intent(out)                       ::      fit     !   mean square error

            integer             ::      ii
            real(kind=real64)   ::      bestFit
            real(kind=real64),dimension(size(a))        ::      besta,bestb

            bestFit = huge(1.0)
            do ii = 1,n
                call random_number(a) 
                call random_number(b) 
                a = Lib_FitDoyleTurner_MINA + a*(Lib_FitDoyleTurner_MAXA-Lib_FitDoyleTurner_MINA)
                b = b*(Lib_FitDoyleTurner_MAXA)
                call fitDoyleTurner0( s,f,a,b, fit )
                if (fit<bestFit) then
                    bestFit = fit
                    besta = a
                    bestb = b
                    write (*,fmt='(a,i6,a,g12.3,a,100f12.6)') "Lib_FitDoyleTurner::fitDoyleTurner1 info - step ",ii," fit ",fit," a,b ",a,b
                end if
            end do
            return
        end subroutine fitDoyleTurner1


        subroutine fitDoyleTurner0( s,f,a,b, fit )
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      given the tabulated scattering function f_i at points s_i
    !*      refine the Doyle Turner coefficients
            real(kind=real64),dimension(:),intent(in)           ::      s       !   (1:nData) scattering points
            real(kind=real64),dimension(:),intent(in)           ::      f       !   (1:nData) scattering amplitudes
            real(kind=real64),dimension(:),intent(inout)        ::      a       !   Doyle Turner a coefficients
            real(kind=real64),dimension(:),intent(inout)        ::      b       !   Doyle Turner b coefficients
            real(kind=real64),intent(out)                       ::      fit     !   mean square error

            real(kind=real64),dimension(size(a))                ::      dfitda,dfitdb
            real(kind=real64),dimension(2*size(a))              ::      dfitdab
            real(kind=real64)                                   ::      mod_dfitdab
            integer                                             ::      nCoeff,nData
            real(kind=real64)                                   ::      lambda
            real(kind=real64)                                   ::      x1,x2,x3,f1,f2,f3
            type(GoldenSection)                     ::      gold
            logical                                 ::      isConverged,isWithinTimeLimit
            integer                                 ::      loop,step,ii


            nData = size(s)
            nCoeff = size(a)

            do loop = 1,Lib_FitDoyleTurner_N_GOLD_DIRECTIONS        !   try downhill search, but don't search forever

            !---    compute the fit and the derivative at this point
                call getFitAndDeriv( s,f,a,b, fit,dfitda,dfitdb )
                dfitdab(1:nCoeff) = dfitda
                dfitdab(nCoeff+1:2*nCoeff) = dfitdb
                mod_dfitdab = norm2(dfitdab)

            !---    have we converged?
                if (fit < Lib_FitDoyleTurner_FIT_TOL*nData ) return                        !   yes, absolute fit is less than target leve
                if (mod_dfitdab < Lib_FitDoyleTurner_DERIV_TOL*nData) return               !   yes, derivative is less than target
                if (fit > 1.0d8*nData ) return                                             !    no, going horribly wrong!
                
            !---    we need to go in the direction -dfitdab to improve the fit. Estimate how far
            !       fit(0) = fit(dx) + fit'(dx).dx
            !       say dx = lambda fit'(dx)
                lambda = -fit / (mod_dfitdab*mod_dfitdab)
                !print *,"fitDoyleTurner0 loop,fit,mod_dfitdab,lambda ",loop,fit,mod_dfitdab,lambda
                x1 = 0.0d0 ; f1 = fit
                x3 = lambda ; f3 = getFit( s,f, a + x3*dfitda, max(0.0d0,b + x3*dfitdb ))
                GOLD_TOL = 1.0d-16 
                gold = GoldenSection_ctor( x3,x1, f3,f1 )
                do step = 1,Lib_FitDoyleTurner_N_GOLD_DIRECTIONS
                    call nextPoint(gold,x2)
                    f2 = getFit( s,f, a + x2*dfitda, max(0.0d0,b + x2*dfitdb ))
                    call minimise(gold,x2,f2,isConverged,isWithinTimeLimit)                    
                    !call report(gold)
                    !print *,x2,f2,isConverged,isWithinTimeLimit
                    if (isConverged) exit
                    if (.not. isWithinTimeLimit) exit
                end do
                call delete(gold)

                a = a + x2*dfitda
                b = b + x2*dfitdb 

            end do
            fit = fit / nData

        !---    as a courtesy, bubble sort coefficients into absolute size order
            do 
                isConverged = .true.
                do ii = 1,nCoeff-1
                    if (abs(a(ii))<abs(a(ii+1))) then
                        x1 = a(ii) 
                        a(ii) = a(ii+1) 
                        a(ii+1) = x1
                        x2 = b(ii) 
                        b(ii) = b(ii+1) 
                        b(ii+1) = x2
                        isConverged = .false.
                    end if
                end do
                if (isConverged) exit
            end do

            return
        end subroutine fitDoyleTurner0


        pure subroutine getFitAndDeriv(s,f,a,b, fit,dfitda,dfitdb )
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !       fit = sum_i ( f(s_i) - f_DT(s_i) )^2
    !       comptue the fit and its derivatives with respect to a and b
            real(kind=real64),dimension(:),intent(in)           ::      s       !   (1:nData) scattering points
            real(kind=real64),dimension(:),intent(in)           ::      f       !   (1:nData) scattering amplitudes
            real(kind=real64),dimension(:),intent(in)           ::      a       !   Doyle Turner a coefficients
            real(kind=real64),dimension(:),intent(in)           ::      b       !   Doyle Turner b coefficients
            real(kind=real64),intent(out)                       ::      fit     !   mean square error
            real(kind=real64),dimension(:),intent(out)          ::      dfitda
            real(kind=real64),dimension(:),intent(out)          ::      dfitdb

            integer                 ::      nData,nCoeff
            integer                 ::      ii,jj
            real(kind=real64)       ::      xx,ff,df

            nData = size(s)
            nCoeff = size(a)

            fit = 0
            dfitda = 0
            dfitdb = 0
            do ii = 1,nData
                ff = DoyleTurnerSum(s(ii),a,b)
                do jj = 1,nCoeff
                    xx = exp( -b(jj)*s(ii)*s(ii) )
                    df = ff - f(ii)
                    fit = fit + df*df
                    dfitda(jj) = dfitda(jj) + 2*df*xx
                    dfitdb(jj) = dfitdb(jj) - 2*df*a(jj)*xx*s(ii)*s(ii)
                end do
            end do

            return
        end subroutine getFitAndDeriv

        pure real(kind=real64) function getFit(s,f,a,b )
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !       fit = sum_i ( f(s_i) - f_DT(s_i) )^2
    !       comptue the fit and its derivatives with respect to a and b
            real(kind=real64),dimension(:),intent(in)           ::      s       !   (1:nData) scattering points
            real(kind=real64),dimension(:),intent(in)           ::      f       !   (1:nData) scattering amplitudes
            real(kind=real64),dimension(:),intent(in)           ::      a       !   Doyle Turner a coefficients
            real(kind=real64),dimension(:),intent(in)           ::      b       !   Doyle Turner b coefficients

            integer                 ::      nData,nCoeff
            integer                 ::      ii
            real(kind=real64)       ::      df,ff

            nData = size(s)
            nCoeff = size(a)

            getFit = 0
            do ii = 1,nData
                ff = DoyleTurnerSum(s(ii),a,b)
                df = ff - f(ii)
                getFit = getFit + df*df                
            end do

            return
        end function getFit

        pure real(kind=real64) function DoyleTurnerSum(s,a,b)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      compute Doyle-Turner sum at scattering point s
            real(kind=real64),intent(in)                        ::      s       !   scattering angle s = q/(4 pi) = sin(theta)/lambda
            real(kind=real64),dimension(:),intent(in)           ::      a       !   Doyle Turner a coefficients
            real(kind=real64),dimension(:),intent(in)           ::      b       !   Doyle Turner b coefficients

            integer                 ::      jj

            DoyleTurnerSum = 0
            do jj = 1,size(a)
                DoyleTurnerSum = DoyleTurnerSum + a(jj) * exp( -b(jj)*s*s ) 
            end do

            return
        end function DoyleTurnerSum

    end module Lib_FitDoyleTurner

        