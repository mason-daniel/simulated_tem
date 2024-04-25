     
    
    program testLib_RelativisticElectrons
!---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
!*  
!*      simple program to test functioning of Lib_RelativisticElectrons
!*      
!*          successful result        
!*          
!*              accelerating voltage V =      150000.0000
!*              wavelength       0.02957041
!*              velocity      1901.64295838
!*              
!*               done
!*              

        use iso_fortran_env
        use Lib_ColouredTerminal
        use Lib_RelativisticElectrons
        implicit none
        
        real(kind=real64)       ::      VV
        
        

        character(len=256),dimension(3)           ::      output
        character(len=*),dimension(3),parameter   ::      output0 = (/  "accelerating voltage V =      150000.0000", &
                                                                        "wavelength       0.02957041              ", &
                                                                        "velocity      1901.64295838              "  /)
                                                                       
                                                                        
        logical                     ::      ok 
        integer                     ::      ii                                                             
                                     
        
        
        VV = 150000.0d0
        write (output(1),fmt='(a,f16.4)') "accelerating voltage V = ",VV
        write (output(2),fmt='(a,f16.8)') "wavelength ",wavelength(VV)
        write (output(3),fmt='(a,f16.8)') "velocity   ",velocity(VV)
        
        
        
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
        
    end program testLib_RelativisticElectrons
