
    program test_Lib_inputTDS
!---^^^^^^^^^^^^^^^^^^^^^^^^^^
        use iso_fortran_env
        use Lib_inputReal
        use Lib_inputTDS
        use Lib_Elements
        use Lib_linint
        use Lib_ColouredTerminal
        implicit none

        include "FE_imag.h"

        real(kind=real64),parameter     ::      EPS = 0.0005d0

        integer                 ::      nS
        real(kind=real64)       ::      smin,smax,ds
        character(len=8)        ::      element = "Fe"
        real(kind=real64)       ::      V = 200000
        real(kind=real64)       ::      T = 300
        character(len=256)      ::      dummy
        real(kind=real64),dimension(:),allocatable      ::      s,fimag
        real(kind=real64)       ::      bdw
        logical         ::      ok
        integer         ::      ii,jj
        real(kind=real64)       ::      ss,aa,ff
        

        print *,"usage: test_Lib_inputTDS [element = ""Fe"" [voltage = 200000 [temperature = 300]]]"

        call get_command_argument(1,dummy)
        if (len_trim(dummy)>0) element = trim(dummy)

        call get_command_argument(2,dummy)
        if (len_trim(dummy)>0) read(dummy,fmt=*) V

        call get_command_argument(3,dummy)
        if (len_trim(dummy)>0) read(dummy,fmt=*) T

        smin = LIB_INPUTTDS_SMIN
        smax = LIB_INPUTTDS_SMAX
        nS = LIB_INPUTTDS_NS
        ds = (smax-smin)/nS 
        !print *,"hello?",smin,smax,ds,ns
        allocate(s(0:nS))
        allocate(fimag(0:nS))
        bdw = DebyeWallerFactor(element,T)

        print *,"test_Lib_inputTDS info - V = ",V," (V) T = ",T," (K)"
        print *,"    element "//trim(element)
        print *,"    D-W factor B ",bdw," U2 = ",bdw/(8*3.14159265390d0**2)

        call getImagScatteringAmplitude(element,V,bdw,s,fimag)
        ok = .true.
        write(*,fmt='(a6,a6,3a16,a4)') "elem","indx","s=q/4pi","f(s)","table(Fe)","ok?"
        do ii = 0,nS
            ss = smin + ii*ds                       !   s value of my scattering amplitude
            jj = floor( ss/DSFE )                   !   where to find s in the table "Fe.h". 
            aa = ( ss - jj*DSFE ) / DSFE            !   Fraction of interval to interpolate
            ff = linint( ffe_ref,jj,aa )
            ok = ok .and. abs(ff-fimag(ii))<EPS
            write(*,fmt='(a6,i6,3f16.8,l4)') element,ii,ss,fimag(ii) , ff,(len_trim(dummy)>0).or.(abs(ff-fimag(ii))<EPS)
        end do


        ok = ok .or. (len_trim(dummy)>0) 

        if (ok) then
            print *,colour(GREEN,"PASS")
        else
            print *,colour(RED,"FAIL")
        end if

    end program test_Lib_inputTDS
