
    program test_Lib_REALfit
!---^^^^^^^^^^^^^^^^^^^^^^^^
!*      simple test program to demonstrate functioning of Lib_REALfit
        use iso_fortran_env
        use Lib_REALfit
        use Lib_inputReal
        use Lib_Elements
        use Lib_FitDoyleTurner
        use Lib_ColouredTerminal
        implicit none


        character(len=8)        ::      element = "Fe"
        character(len=256)      ::      dummy
        integer                 ::      nS
        real(kind=real64)       ::      smin,smax,ds        
        real(kind=real64),dimension(5)              ::      a0,b0           !   Doyle-Turner coefficients (before)
        real(kind=real64),dimension(5)              ::      a,b             !   Doyle-Turner coefficients
        real(kind=real64),dimension(:),allocatable  ::      f,s               !   Scattering amplitudes
        !real(kind=real64)       ::    chisq_bar  
        real(kind=real64)       ::      mse,mse0,mse1
        character(len=256)      ::  outline
        character(len=256)      ::  expectedline = ""
        logical         ::      ok
        

        print *,"usage: test_Lib_REALfit [element = ""Fe""]"

        call get_command_argument(1,dummy)
        if (len_trim(dummy)>0) element = trim(dummy)

        smin = LIB_INPUTREAL_SMIN
        smax = LIB_INPUTREAL_SMAX
        nS = LIB_INPUTREAL_NS
        ds = (smax-smin)/nS
        allocate(s(0:ns))
        allocate(f(0:nS))
        call getScatteringAmplitude(element,s,f)
        a0 = getrealDoyleTurner_a(element)
        b0 = getrealDoyleTurner_b(element)
        a = a0
        b = b0
!        call fitDoyleTurnerCoefficients( a,b, f, chisq_bar )

!        print *,"test_Lib_REALfit info - "//trim(element)," chi^2 = ",chisq_bar
        write (*,fmt='(a,10f12.6)') "before ",a0,b0


        mse0 = getFit(s,f,a0,b0 )
        mse1 = getFit(s,f,a,b )
        a = a0
        b = b0
        call fitDoyleTurner( s,f,a,b, mse )
        write (*,fmt='(a,10f12.6)') "after  ",a,b
        print *,"test_Lib_REALfit info - "//trim(element)," mean square error before = ",mse0 ," after ",mse1,mse



        write (outline,fmt='(a6,a4,5f12.6,a4,5f12.6)') element,"a",a,"b",b
        print *,trim(outline)
        expectedLine = "Fe       a    2.377012    2.337033    1.748293    0.559985    0.133168   b   67.523278   16.782146    3.380222    0.693041    0.077519"

        ok = (len_trim(dummy)>0) .or. (mse < 1.0d-4)

        if (ok) then
            print *,colour(GREEN,"PASS")
        else
            print *,colour(RED,"FAIL")
        end if


    end program test_Lib_REALfit