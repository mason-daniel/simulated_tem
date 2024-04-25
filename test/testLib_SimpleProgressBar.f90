!   gfortran -ffree-line-length-256 src/Lib_SimpleProgressBar.f90 src/testLib_SimpleProgressBar.f90 -o Test/testLib_SimpleProgressBar.exe

    program testLib_SimpleProgressBar
!---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
!*      simple program to test functioning of Lib_SimpleProgressBar
!*
!*      correct functioning
!*
!*          $ ./Test/testLib_SimpleProgressBar.exe                
!*            0%  10%  20%  30%  40%  50%  60%  70%  80%  90% 100%
!*            0%  10%  20%  30%  40%  50%  60%  70%  80%  90% 100%
!*            0%  10%  20%  30% 100%                              
!*                                                                
!*           done                                                                                                  
                                                               
        use Lib_SimpleProgressBar
        implicit none
        
        
        integer         ::       nn
         
        integer                     ::      ii          
        
        
        nn = 1000
        do ii = 1,nn
            call progressBar(ii,nn)
        end do
         
        
        nn = 3917
        do ii = 1,nn
            call progressBar(ii,nn)
        end do
        
        
        nn = 3
        do ii = 1,nn
            call progressBar(ii,nn)
        end do
        
        print *,""
        print *,"done"
        print *,""
        
    end program testLib_SimpleProgressBar