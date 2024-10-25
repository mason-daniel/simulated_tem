
    module Lib_RelativisticElectrons
!---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
        use iso_fortran_env
        implicit none
        private


        real(kind=real64),parameter,public      ::      ME =  0.056856301036d0      !   electron mass in weird units eV A^-2 fs^2
        real(kind=real64),parameter,private     ::      E  =  1.0d0                 !   coulomb in eV 
        real(kind=real64),parameter,private     ::      C  =  2997.92458d0          !   speed of light in A fs^-1
        real(kind=real64),parameter,public      ::      H  =  4.135667696d0         !   Planck constant in eV fs
        real(kind=real64),parameter,public      ::      HBAR =  H/(2*3.14159265390d0)         !   Dirac constant in eV fs
        
        
        
        public          ::      wavelength
        public          ::      velocity
        
    contains
!---^^^^^^^^


        pure function wavelength( V ) result( lambda )
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      returns electron wavelength given accelerating voltage V
            real(kind=real64),intent(in)            ::      V
            real(kind=real64)                       ::      lambda
            
            if (V<=0) then
                lambda = 0
            else
                lambda = 2*ME*E*V*( 1 + E*V/(2*ME*C*C) )
                lambda = H/sqrt( lambda )
            end if
            
            return
        end function wavelength
        
        pure function velocity( V ) result( u )
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      returns electron velocity given accelerating voltage V
            real(kind=real64),intent(in)            ::      V
            real(kind=real64)                       ::      u
            
            u = 1 + E*V/(ME*C*C)
            u = C*sqrt( 1 - 1/(u*u) )
            
            return
        end function velocity
        
            
    end module Lib_RelativisticElectrons      
         