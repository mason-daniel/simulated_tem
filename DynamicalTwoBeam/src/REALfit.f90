
    program REALfit
!---^^^^^^^^^^^^^^^
!*      program to fit Doyle-Turner coefficients using Lib_REALfit
        use iso_fortran_env
        use Lib_REALfit
        use Lib_inputReal
        use Lib_Elements
        !use Lib_ColouredTerminal
        implicit none


        character(len=8)        ::      element  
        integer                 ::      ee
        integer                 ::      nS
        real(kind=real64)       ::      smin,smax,ds        
        real(kind=real64),dimension(5)              ::      a0,b0           !   Doyle-Turner coefficients (before)
        real(kind=real64),dimension(5)              ::      a,b             !   Doyle-Turner coefficients
        real(kind=real64),dimension(:),allocatable  ::      f               !   Scattering amplitudes
        real(kind=real64)       ::    chisq_bar  
         

        print *,"usage: REALfit"

        smin = LIB_INPUTREAL_SMIN
        smax = LIB_INPUTREAL_SMAX
        nS = LIB_INPUTREAL_NS
        ds = (smax-smin)/nS
        allocate(f(0:nS))

        do ee = 1,numberOfElements()

            element = getElementName(ee)

            call getScatteringAmplitude(element,f)
            a0 = getDoyleTurner_a(element)
            b0 = getDoyleTurner_b(element)
            a = a0
            b = b0
            call fitDoyleTurnerCoefficients( a,b, f, chisq_bar )

            !print *,"test_Lib_REALfit info - "//trim(element)," chi^2 = ",chisq_bar
            !write (*,fmt='(a,10f12.6)') "before ",a0,b0
            write (*,fmt='(a4,11f12.6)') trim(element),a,b,chisq_bar

        end do
              
        print *,""
        print *,"done"
        print *,""


    end program REALfit