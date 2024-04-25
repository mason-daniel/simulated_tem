!   gfortran -ffree-line-length-256 src/NBAX_StringTokenizers.f90 src/testNBAX_StringTokenizers.f90 -o Test/testNBAX_StringTokenizers.exe

    program testNBAX_StringTokenizers
!---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
!*      test correct working of NBAX_StringTokenizers 
!*
!*          correct functioning 
!*                  
!*              $ ./Test/testNBAX_StringTokenizers.exe                                                                                                      
!*               test string 'hello,   "Freddy". How are you?' parsed with delimiters ', .\"'                                                               
!*               hello                                                                                                                                      
!*               Freddy                                                                                                                                     
!*               How                                                                                                                                        
!*               are                                                                                                                                        
!*               you?                                                                                                                                       
!*               test string '1 -21 0 123456' converted to integers                                                                                         
!*                         1         -21           0      123456                                                                                            
!*               converted back to string '1 -21 0 123456 '                                                                                                 
!*               test string '1.0d0 -2.0e-2 123456 0 -0.0 123456789.987654321' converted to reals                                                           
!*                     1.00000000000         -0.200000000000E-01       123456.000000           0.00000000000           0.00000000000           123456789.988
!*               converted back to string '  1.000000 -0.020000  123456.000000  0.000000  0.000000  123456789.987654'                                       
!*                                                                                                                                                          
!*               done                                                                                                                                       
!*                                                                                                                                                          

       use NBAX_StringTokenizers
       use Lib_ColouredTerminal   
       use iso_fortran_env
       implicit none

       type(StringTokenizer)       ::      st
       character(len=256)          ::      dummy  
       character(len=64)           ::      token
       integer                          ::      nn
       integer,dimension(100)           ::      iarray
       real(kind=real64),dimension(100) ::      rarray
       
       
        character(len=256),dimension(12)           ::      output
        character(len=*),dimension(12),parameter   ::      output0 = (/ "test string 'hello,   ""Freddy"". How are you?' parsed with delimiters ', .\""'                                                               ", &
                                                                        "hello                                                                                                                                      ", &
                                                                        "Freddy                                                                                                                                     ", &
                                                                        "How                                                                                                                                        ", &
                                                                        "are                                                                                                                                        ", &
                                                                        "you?                                                                                                                                       ", &   
                                                                        "test string '1 -21 0 123456' converted to integers                                                                                         ", &
                                                                        "          1         -21           0      123456                                                                                            ", &
                                                                        "converted back to string '1 -21 0 123456 '                                                                                                 ",  &
                                                                        "test string '1.0d0 -2.0e-2 123456 0 -0.0 123456789.987654321' converted to reals                                                           ",  &
                                                                        "      1.00000000000         -0.200000000000E-01       123456.000000           0.00000000000           0.00000000000           123456789.988",  &
                                                                        "converted back to string '  1.000000 -0.020000  123456.000000  0.000000  0.000000  123456789.987654'                                       " /)
                                                                        
                                                                        
                                                                        
                                                                        
                                                                        
                                                                         
        logical                     ::      ok 
        integer                     ::      ii          
        
       
       
       
       
       
       
       
       
       
       
       
       dummy = 'hello,   "Freddy". How are you?'
       output(1) = "test string '"//trim(dummy)//"' parsed with delimiters ', .\""' "
       st = StringTokenizer_ctor(dummy,', .\"')
       ii = 1
       do
           call nextToken( st,token )
           ii = ii + 1
           output(ii) = trim(token)
           if (.not. hasMoreTokens(st)) exit
       end do
       
       dummy = '1 -21 0 123456'
       output(7) = "test string '"//trim(dummy)//"' converted to integers"       
       call parse( dummy,iarray,nn )
       write (output(8),fmt='(100i12)') iarray(1:nn)
       call itoa(iarray(1:nn),12,dummy,nn) 
       output(9) = "converted back to string '"//dummy(1:nn)//"'"       
       
       
       dummy = '1.0d0 -2.0e-2 123456 0 -0.0 123456789.987654321'
       output(10) = "test string '"//trim(dummy)//"' converted to reals"       
       call parse( dummy,rarray,nn )
       write (output(11),fmt='(100g24.12)') rarray(1:nn)
       call ftoa(rarray(1:nn),12,6,dummy,nn) 
       output(12) = "converted back to string '"//dummy(1:nn)//"'"       
       
       
        
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
       
   end program testNBAX_StringTokenizers
