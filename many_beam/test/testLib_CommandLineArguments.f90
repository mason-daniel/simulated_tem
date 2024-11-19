
!   gfortran -ffree-line-length-256 src/NBAX_StringTokenizers.f90 src/Lib_CommandLineArguments.f90 src/testLib_CommandLineArguments.f90 -o Test/testLib_CommandLineArguments.exe

    program testLib_CommandLineArguments
!---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
!*      A simple test program to show correct working of Lib_CommandLineArguments
!*          
!*      correct functioning ( example 1 )
!*  
!*          $ ./Test/testLib_CommandLineArguments.exe                                                          
!*          testLib_CommandLineArguments.exe                                                                   
!*          tests functioning of Lib_CommandLineArguments                                                      
!*                                                                                                             
!*          usage                                                                                              
!*          ^^^^^                                                                                              
!*              scalars                                                                                        
!*                  -a <int>            integer value                                                          
!*                  [-c <float> ]       real value [ default 1.000 ]                                           
!*                  [-e ]               logical value [ default F ]                                            
!*                  [-f <char> ]        string value [ default "sandwiches, prawn" ]                           
!*                                                                                                             
!*              arrays                                                                                         
!*                  [-b <int_array> ]   integer array [ default 0,0,0,... ]                                    
!*                  [-d <float_array> ] real array value with min 3 values [ default 0.0000,0.0000,0.0000,... ]
!*                                                                                                             
!*      correct functioning ( example 2 )        
!*          
!*          $ ./Test/testLib_CommandLineArguments.exe -a 12345
!*          testLib_CommandLineArguments.exe
!*          tests functioning of Lib_CommandLineArguments
!*          
!*          a =       12345
!*          b =
!*          c =        1.00000000
!*          d =         0.0000000       0.0000000       0.0000000
!*          e =       F
!*          f = "sandwiches, prawn"
!*          
!*           done
!*          
!*      correct functioning ( example 3 )     
!*  
!*          dmason@L1233 ~/Programs/Grains                                                                                                                
!*          $ ./Test/testLib_CommandLineArguments.exe -a 12345 -b 1,2,3,4,5 -c 1.2345E+06 -d 1.0e-5,-0.0,0.0,12345.6789,1.0e+5 -e -f "sandwiches, herring"
!*          testLib_CommandLineArguments.exe                                                                                                              
!*          tests functioning of Lib_CommandLineArguments                                                                                                 
!*                                                                                                                                                        
!*          a =       12345                                                                                                                               
!*          b =           1         2         3         4         5                                                                                       
!*          c =  1234500.00000000                                                                                                                         
!*          d =         0.0000100      -0.0000000       0.0000000   12345.6789000  100000.0000000                                                         
!*          e =       T                                                                                                                                   
!*          f = "sandwiches, herring"                                                                                                                     
!*                                                                                                                                                         
!*          done                                                                                                                                         
!*                                                                                                                                                        
       
        use Lib_CommandLineArguments
        use Lib_ColouredTerminal
        use iso_Fortran_env
        implicit none
        
        character(len=8)                ::      VERSION =  "1.0.0"
        type(CommandLineArguments)      ::      cla
        
        integer                         ::      a = LIB_CLA_NODEFAULT_I
        integer,dimension(10)           ::      b = 0
        real(kind=real64)               ::      c = 1.0d0
        real(kind=real64),dimension(10) ::      d = 0.0d0
        logical                         ::      e = .false.
        character(len=256)              ::      f = "sandwiches, prawn"
        
        integer                         ::      nb = 0
        integer                         ::      nd = 3
        
        character(len=256),dimension(6)             ::      output
        character(len=*),dimension(6),parameter     ::      output0 = (/ "a =       12345                                                                      ",           &
                                                                         "b =           1         2         3         4         5                              ",           &
                                                                         "c =  1234500.00000000                                                                ",           &
                                                                         "d =         0.0000100       0.0000000       0.0000000   12345.6789000  100000.0000000",           &
                                                                         "e =       T                                                                          ",           &
                                                                         "f = ""sandwiches, herring""                                                            "        /) 
        logical             ::      ok
        integer             ::      ii 
        
        
        
        cla = CommandLineArguments_ctor(10)         !   allocates space for up to 10 arguments
        
        !   optionally add a number of subcategories to organise options
        call setCategories(cla,(/ "scalars","arrays " /))
        
        !   optionally request that the first few arguments are handled elsewhere, and not of form "-key value" 
        !   (eg might want ./a.exe <filename> [-option1 ... ] - in which case ignore the first argument with)
        
        !   optionally might want to set the program name and some description          
        call setProgramDescription( cla, "testLib_CommandLineArguments.exe \ntests functioning of Lib_CommandLineArguments" )
        call setProgramVersion( cla, VERSION )
        
        
        !   add an option, return its value
        call get( cla,"a",a,LIB_CLA_REQUIRED,"   integer value",1 )
        call get( cla,"b",b,nb,LIB_CLA_OPTIONAL," integer array ",2 )
        call get( cla,"c",c,LIB_CLA_OPTIONAL," real value",1 )
        call get( cla,"d",d,nd,LIB_CLA_OPTIONAL," real array value with min 3 values",2 )
        call get( cla,"e",e,LIB_CLA_OPTIONAL," logical value",1 )
        call get( cla,"f",f,LIB_CLA_OPTIONAL," string value",1 )
        
        !   report the options available with 
        call report(cla)               
        
        !   optionally eg exit
        if (.not. allRequiredArgumentsSet(cla)) stop
        if (hasHelpArgument(cla)) stop
        
        call delete(cla)    
            
        
        write(output(1),fmt='(a,i10)')      "a =  ",a
        write(output(2),fmt='(a,10i10)')    "b =  ",b(1:nb)
        write(output(3),fmt='(a,f16.8)')    "c =  ",c
        write(output(4),fmt='(a,10f16.7)')  "d =  ",d(1:nd)
        write(output(5),fmt='(a,l6)')       "e =  ",e
        write(output(6),fmt='(a,a)')        "f = """//trim(f)//""""    
        
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
         
        
    end program testLib_CommandLineArguments