
    module Lib_DoyleTurnerCoefficients
!---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
!*      return the real and imaginary Doyle-Turner coefficients 
!*      for the elements tabulated by  
!*      S.L. Dudarev et al./Surface Science 330 (1995) 86-100 table 1

        use iso_fortran_env
        use Lib_Elements
        implicit none
        private

        include "DoyleTurnerCoefficients.h"


        public          ::          getDoyleTurner_a
        public          ::          getDoyleTurner_b


    contains
!---^^^^^^^^


        function getDoyleTurner_a(el) result(a)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      return the real Doyle-Turner coefficients a(1:5) for element el
            character(len=*),intent(in)         ::      el
            real(kind=real64),dimension(5)      ::      a

            character(len=ELNAME_LEN)           ::      el_tidy
            integer         ::      ee

            a = 0

            el_tidy = getElementName(whichElement(el))          !   convert to standard format.
            if (el_tidy == "") then
                print *,"Lib_DoyleTurnerCoefficients::getDoyleTurner_a error - did not recognise input element """//trim(el)//""""                
                return
            end if


            do ee = 1,N_DOYLETURNER_ELEMENTS
                if (el_tidy == getElementName(whichElement(DoyleTurner_elements(ee)))) then
                    !   the input element is indexed ee in the tables
                    a(1:5) = DoyleTurner_a(1:5,ee)
                    return
                end if        
            end do

            return
        end function getDoyleTurner_a


        function getDoyleTurner_b(el) result(b)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      return the imaginary Doyle-Turner coefficients b(1:5) for element el
            character(len=*),intent(in)         ::      el
            real(kind=real64),dimension(5)      ::      b

            character(len=ELNAME_LEN)           ::      el_tidy
            integer         ::      ee

            b = 0

            el_tidy = getElementName(whichElement(el))          !   convert to standard format.
            if (el_tidy == "") then
                print *,"Lib_DoyleTurnerCoefficients::getDoyleTurner_b error - did not recognise input element """//trim(el)//""""                
                return
            end if


            do ee = 1,N_DOYLETURNER_ELEMENTS
                if (el_tidy == getElementName(whichElement(DoyleTurner_elements(ee)))) then
                    !   the input element is indexed ee in the tables
                    b(1:5) = DoyleTurner_b(1:5,ee)
                    return
                end if        
            end do

            return
        end function getDoyleTurner_b        


    end module Lib_DoyleTurnerCoefficients

