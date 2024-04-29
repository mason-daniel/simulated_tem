
!   gfortran -ffree-line-length-256 src/Lib_Quicksort.f90 src/testLib_Quicksort.f90 -o Test/testLib_Quicksort.exe

    program testLib_Quicksort
!---^^^^^^^^^^^^^^^^^^^^^^^^^^
!*      test correct functioning of Lib_Quicksort
!*
!*      correct functioning 
!*
!*          $ ./Test/testLib_Quicksort.exe                                                                                                                                     
!*           original array                                                                                                                                                 
!*                0.36128190      0.60206431      0.94111145      0.95778430      0.91291261      0.31338051      0.24320246      0.32609701      0.48893872      0.09445132
!*           element ordering (small to large)                                                                                                                              
!*                        10               7               6               8               1               9               2               5               3               4
!*           reordering list (large to small)                                                                                                                               
!*                0.95778430      0.94111145      0.91291261      0.60206431      0.48893872      0.36128190      0.32609701      0.31338051      0.24320246      0.09445132
!*                                                                                                                                                                          
!*           done                                                                                                                                                           
!*

        use iso_fortran_env
        use Lib_ColouredTerminal
        use Lib_Quicksort
        implicit none
        
        
        real(kind=real64),dimension(10)     ::      array = (/0.361281889,0.602064309,0.941111459,0.957784303,0.912912632,0.313380521,0.243202457,0.326097002,0.488938725,0.09445132 /)
        integer,dimension(10)               ::      indx
         
        

        character(len=256),dimension(3)            ::      output
        character(len=*),dimension(3),parameter   ::      output0 = (/  "0.36128190      0.60206431      0.94111145      0.95778430      0.91291261      0.31338051      0.24320246      0.32609701      0.48893872      0.09445132  ", &
                                                                        "10               7               6               8               1               9               2               5               3               4          ", &
                                                                        "0.95778430      0.94111145      0.91291261      0.60206431      0.48893872      0.36128190      0.32609701      0.31338051      0.24320246      0.09445132  "  /)
                                                                        
                                                                       
                                                                        
        logical                     ::      ok 
        integer                     ::      ii                                                             
                                                                   
                               
        
        
        print *,"original array"
        write (output(1),fmt='(10f16.8)') array(:)
        
        
        print *,"element ordering (small to large)"
        do ii = 1,10 ; indx(ii) = ii ; end do
        call quicksort(array,indx)
        write (output(2),fmt='(10i16)') indx(:)
        
        print *,"reordering list (large to small)"
        call quicksort(array,back=.true.)
        write (output(3),fmt='(10f16.8)') array(:)
        
        
        
        
    !---    here is the simple test: does the output look like the stored output?
        ok = .true.
        do ii = 1,size(output)            
            if ( trim(cutSpaces(output(ii))) == trim(cutSpaces(output0(ii))) ) then
                write (*,fmt='(a)') trim(cutSpaces(output(ii)))
            else
                ok = .false.
                write (*,fmt='(a)') trim(cutSpaces(output0(ii)))//"    "//colour(RED,trim(cutSpaces(output(ii))))
            end if
        end do   
        
       
    !---    output the result "PASS" or "FAIL"    
        if (ok) then
            print *,colour(LIGHT_GREEN,"PASS")
        else
            print *,colour(RED,"FAIL")
        end if
        
        print *,""
        print *,"done"
        print *,""
        
    end program testLib_Quicksort
        

        