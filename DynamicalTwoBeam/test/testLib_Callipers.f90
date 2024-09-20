 
    program testLib_Callipers
!---^^^^^^^^^^^^^^^^^^^^^^^^^
!*      A simple test program illustrating the use of Lib_Callipers for timing purposes
!*      A successful run should look (something) like
!*
!*      $ time ./Test/testLib_Callipers.exe     
!*      Callipers [t=             1.007790200 s]
!*                                              
!*       done                                   
!*                                              
!*                                              
!*      real    0m1.083s                        
!*      user    0m1.015s                        
!*      sys     0m1.015s                        
!*                                              
!*                                              
!*
!*
        use Lib_ColouredTerminal
        use Lib_Callipers
        implicit none
        
        type(Callipers)     ::      c
        real(kind=8)        ::      xx
        logical             ::      ok
        
        c = Callipers_ctor()
        
        call system("sleep 1") 
        
        call report(c)
        
        xx = elapsed(c)
        ok = ( ( xx > 0.9d0 ).and.( xx < 9.0d0 ) )      !   check sleep 1 is about 1 second.
        
    !---    output the result "PASS" or "FAIL"    
        if (ok) then
            print *,colour(LIGHT_GREEN,"PASS")
        else
            print *,colour(RED,"FAIL")
        end if
         
    end program testLib_Callipers