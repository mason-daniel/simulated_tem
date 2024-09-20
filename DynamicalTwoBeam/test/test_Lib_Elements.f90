
    program test_Lib_Elements
!---^^^^^^^^^^^^^^^^^^^^^^^^^
!*      simple test code illustrating functioning of Lib_Elements
!*
!*      Daniel Mason
!*      (c) UKAEA Spet 2024
!*

        use iso_fortran_env
        use Lib_ColouredTerminal
        use Lib_Elements
        implicit none
        


        character(len=8)        ::  element = "Fe"
        character(len=256)      ::  outline,dummy 
        character(len=256)      ::  expectedline = "   15    Fe   bcc      2.8665      2.8665      2.8665     55.8400" 
        logical                 ::  ok
    
        print *,"usage: test_Lib_Elements [element = ""Fe""]"

        call get_command_argument(1,dummy)
        if (len_trim(dummy)>0) element = trim(dummy)
          
        
        
        write(outline,fmt='(i6,2a6,4f12.4)') whichElement(element)                  &
                                ,getElementName(whichElement(element))              &
                                ,getLatticeName(element)                            &
                                ,getLatticeConstant(element)                        &
                                ,getMass(element)

        write(unit=*,fmt='(3a6,4a12)') "indx","name","latt","latt const","","","mass"
        write(unit=*,fmt='(a)') trim(outline)

        ok = (len_trim(dummy)>0) .or. (cutspaces(outline) == cutspaces(expectedline))

        if (ok) then
            print *,colour(GREEN,"PASS")
        else
            print *,colour(RED,"FAIL")
        end if

    end program test_Lib_Elements


