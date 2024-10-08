
    program test_Lib_TDSfit
!---^^^^^^^^^^^^^^^^^^^^^^^^
!*      simple test program to demonstrate functioning of Lib_TDSfit
        use iso_fortran_env
        use Lib_TDSfit
        use Lib_inputTDS
        use Lib_Elements
        use Lib_ColouredTerminal
        use Lib_FitDoyleTurner
        implicit none


        character(len=8)        ::      element = "Fe"
        character(len=256)      ::      dummy
        integer                 ::      nS
        real(kind=real64)       ::      smin,smax,ds        
        real(kind=real64)       ::      V = 200000
        real(kind=real64)       ::      T = 300
        real(kind=real64),dimension(5)              ::      a0,b0           !   Doyle-Turner coefficients (before)
        real(kind=real64),dimension(5)              ::      a,b             !   Doyle-Turner coefficients
        real(kind=real64),dimension(:),allocatable  ::      s,fimag           !   Scattering amplitudes
        real(kind=real64)       ::      mse,mse0
       ! character(len=256)      ::      outline
       ! character(len=256)      ::      expectedline = ""
        logical                 ::      ok
        real(kind=real64)       ::      bdw
        integer                 ::      ii

        print *,"usage: test_Lib_TDSfit [element = ""Fe"" [voltage = 200000 [temperature = 300]]]"

        call get_command_argument(1,dummy)
        if (len_trim(dummy)>0) element = trim(dummy)

        smin = LIB_INPUTTDS_SMIN
        smax = LIB_INPUTTDS_SMAX
        nS = LIB_INPUTTDS_NS
        ds = (smax-smin)/nS
        allocate(s(0:ns))
    

        allocate(fimag(0:nS))
        bdw = DebyeWallerFactor(element,T)

        print *,"test_Lib_inputTDS info - V = ",V," (V) T = ",T," (K)"
        print *,"    element "//trim(element)
        print *,"    D-W factor B ",bdw," U2 = ",bdw/(8*3.14159265390d0**2)
        call getImagScatteringAmplitude(element,V,bdw,s,fimag)


        a0 = getscatDoyleTurner_a(element,T)
        b0 = getscatDoyleTurner_b(element,T)
        

        mse0 = getFit(s,fimag,a0,b0 )

        a = a0
        b = b0
        !print *,"test_Lib_inputTDS info - a ",a0
        !print *,"test_Lib_inputTDS info - b ",b0

        call fitDoyleTurner( s,fimag,a,b, mse )

        !call TDSfit( a,b, fimag, chisq_bar )
       
         

        print *,"test_Lib_TDSfit info - "//trim(element)," mean square error before = ",mse0 ," after ",mse

        
        write (*,fmt='(a,10f12.6)') "before ",a0,b0
        write (*,fmt='(a,10f12.6)') "after  ",a,b
        



       ! write (outline,fmt='(a6,a4,5f12.6,a4,5f12.6)') element,"a",a,"b",b
       ! print *,trim(outline)
       ! expectedLine = "Fe       a    2.377012    2.337033    1.748293    0.559985    0.133168   b   67.523278   16.782146    3.380222    0.693041    0.077519"

        

        write(*,fmt='(4a16)') "s","f(s)","before","after"
        do ii = 0,nS
            write (*,fmt='(4f16.8)') s(ii),fimag(ii),DoyleTurnerSum(s(ii),a0,b0),DoyleTurnerSum(s(ii),a,b)
        end do

        ok = (len_trim(dummy)>0) .or. (mse<1.0d-4)

        if (ok) then
            print *,colour(GREEN,"PASS")
        else
            print *,colour(RED,"FAIL")
        end if


    end program test_Lib_TDSfit