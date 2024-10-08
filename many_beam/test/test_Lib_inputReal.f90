
    program test_Lib_inputReal
!---^^^^^^^^^^^^^^^^^^^^^^^^^^
        use iso_fortran_env
        use Lib_inputReal
        use Lib_linint
        use Lib_ColouredTerminal
        implicit none

        include "FE.h"

        real(kind=real64),parameter     ::      EPS = 0.0005d0

        integer                 ::      nS
        real(kind=real64)       ::      smin,smax,ds
        character(len=8)        ::      element = "Fe"
        character(len=256)      ::      dummy
        real(kind=real64),dimension(:),allocatable      ::      f,s

        logical         ::      ok
        integer         ::      ii,jj
        real(kind=real64)       ::      ss,aa,ff


        print *,"usage: test_Lib_inputReal [element = ""Fe""]"

        call get_command_argument(1,dummy)
        if (len_trim(dummy)>0) element = trim(dummy)

        smin = LIB_INPUTREAL_SMIN
        smax = LIB_INPUTREAL_SMAX
        nS = LIB_INPUTREAL_NS
        ds = (smax-smin)/nS 

        allocate(f(0:nS))
        allocate(s(0:nS))
        call getScatteringAmplitude(element,s,f)
        ok = .true.
        write(*,fmt='(a6,a6,3a12,a4)') "elem","indx","s=q/4pi","f(s)","table(Fe)","ok?"
        do ii = 0,nS
            ss = smin + ii*ds                       !   s value of my scattering amplitude
            jj = floor( ss/DSFE )                   !   where to find s in the table "Fe.h". 
            aa = ( ss - jj*DSFE ) / DSFE            !   Fraction of interval to interpolate
            ff = linint( ffe_ref,jj,aa )
            ok = ok .and. abs(ff-f(ii))<EPS
            write(*,fmt='(a6,i6,3f12.6,l4)') element,ii,ss,f(ii) , ff,(len_trim(dummy)>0).or.(abs(ff-f(ii))<EPS)
        end do


        ok = ok .or. (len_trim(dummy)>0) 

        if (ok) then
            print *,colour(GREEN,"PASS")
        else
            print *,colour(RED,"FAIL")
        end if

    end program test_Lib_inputReal
