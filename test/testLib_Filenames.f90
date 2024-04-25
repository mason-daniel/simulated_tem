
!   gfortran -ffree-line-length-256 src/Lib_Filenames.f90 src/testLib_Filenames.f90 -o Test/testLib_Filenames.exe

   
    program testLib_Filenames
!---^^^^^^^^^^^^^^^^^^^^^^^^^
!*      test code to prove correct working of Lib_Filenames
!*
!*      successful operation
!*
!*          $ ./Test/testLib_Filenames.exe                                              
!*                                                                                      
!*           getSuffix()                                                                
!*           getSuffix("directory/foo")                  = ""              
!*           getSuffix("directory/foo.")                 = ""             
!*           getSuffix("directory/foo.bar")              = "bar"          
!*           getSuffix("foo.bar")                        = "bar"                    
!*           getSuffix("directory.dir/foo.bar")          = "bar"      
!*           getSuffix("directory/foo.bar.txt")          = "txt"      
!*                                                                                      
!*           removeSuffix()                                                             
!*           removeSuffix("directory/foo")               = "directory/foo"              
!*           removeSuffix("directory/foo.")              = "directory/foo"             
!*           removeSuffix("directory/foo.bar")           = "directory/foo"          
!*           removeSuffix("foo.bar")                     = "foo"                    
!*           removeSuffix("directory.dir/foo.bar")       = "directory.dir/foo"      
!*           removeSuffix("directory/foo.bar.txt")       = "directory/foo.bar"      
!*                                                                                      
!*           numberFile()                                                               
!*           numberFile("directory/foo",123)             = "directory/foo.00123"       
!*           numberFile("directory/foo",123,"")          = "directory/foo.00123"       
!*           numberFile("directory/foo",123,"png")       = "directory/foo.00123.png"    
!*           numberFile("directory/foo.bar",123)         = "directory/foo.00123.bar"   
!*           numberFile("directory/foo.bar",123,"")      = "directory/foo.bar.00123"   
!*           numberFile("directory/foo.bar",123,"png")   = "directory/foo.bar.00123.png"
!*                                                                                      
!*           done                                                                       
!*                                                                                      

       use Lib_Filenames
       use Lib_ColouredTerminal
       implicit none
       
       
       
       

        character(len=256),dimension(21) ::      output
        character(len=*),dimension(21),parameter   ::      output0 = (/  "getSuffix()                                                                       "      ,         &
                                                                        "getSuffix(""directory/foo"")                  = """"                                  "  ,         &
                                                                        "getSuffix(""directory/foo."")                 = """"                                  "  ,         &
                                                                        "getSuffix(""directory/foo.bar"")              = ""bar""                               "  ,         &
                                                                        "getSuffix(""foo.bar"")                        = ""bar""                               "  ,         &
                                                                        "getSuffix(""directory.dir/foo.bar"")          = ""bar""                               "  ,         &
                                                                        "getSuffix(""directory/foo.bar.txt"")          = ""txt""                               "  ,         &                                                            
                                                                        "removeSuffix()                                                                    "      ,         &
                                                                        "removeSuffix(""directory/foo"")               = ""directory/foo""                     "  ,         &
                                                                        "removeSuffix(""directory/foo."")              = ""directory/foo""                     "  ,         &
                                                                        "removeSuffix(""directory/foo.bar"")           = ""directory/foo""                     "  ,         &
                                                                        "removeSuffix(""foo.bar"")                     = ""foo""                               "  ,         &
                                                                        "removeSuffix(""directory.dir/foo.bar"")       = ""directory.dir/foo""                 "  ,         &
                                                                        "removeSuffix(""directory/foo.bar.txt"")       = ""directory/foo.bar""                 "  ,         &                                                                           
                                                                        "numberFile()                                                                      "      ,         &
                                                                        "numberFile(""directory/foo"",123)             = ""directory/foo.00123""               "  ,         &
                                                                        "numberFile(""directory/foo"",123,"""")          = ""directory/foo.00123""               ",         &
                                                                        "numberFile(""directory/foo"",123,""png"")       = ""directory/foo.00123.png""           ",         &
                                                                        "numberFile(""directory/foo.bar"",123)         = ""directory/foo.00123.bar""           "  ,         &
                                                                        "numberFile(""directory/foo.bar"",123,"""")      = ""directory/foo.bar.00123""           ",         &
                                                                        "numberFile(""directory/foo.bar"",123,""png"")   = ""directory/foo.bar.00123.png""       " /)
        
         
        integer                     ::      ii                                                             
        logical                     ::      ok                                                             
                                                      
               
       
       
       
       
       
       output(1)  = "getSuffix()"
       output(2)  = "getSuffix(""directory/foo"")                  = """//trim(getSuffix("directory/foo"))//""""
       output(3)  = "getSuffix(""directory/foo."")                 = """//trim(getSuffix("directory/foo."))//""""
       output(4)  = "getSuffix(""directory/foo.bar"")              = """//trim(getSuffix("directory/foo.bar"))//""""
       output(5)  = "getSuffix(""foo.bar"")                        = """//trim(getSuffix("foo.bar"))//""""
       output(6)  = "getSuffix(""directory.dir/foo.bar"")          = """//trim(getSuffix("directory.dir/foo.bar"))//""""
       output(7)  = "getSuffix(""directory/foo.bar.txt"")          = """//trim(getSuffix("directory/foo.bar.txt"))//""""
       
       output(8)  = "removeSuffix()"
       output(9)  = "removeSuffix(""directory/foo"")               = """//trim(removeSuffix("directory/foo"))//""""
       output(10) = "removeSuffix(""directory/foo."")              = """//trim(removeSuffix("directory/foo."))//""""
       output(11) = "removeSuffix(""directory/foo.bar"")           = """//trim(removeSuffix("directory/foo.bar"))//""""
       output(12) = "removeSuffix(""foo.bar"")                     = """//trim(removeSuffix("foo.bar"))//""""
       output(13) = "removeSuffix(""directory.dir/foo.bar"")       = """//trim(removeSuffix("directory.dir/foo.bar"))//""""
       output(14) = "removeSuffix(""directory/foo.bar.txt"")       = """//trim(removeSuffix("directory/foo.bar.txt"))//""""

       output(15) = "numberFile()"
       output(16) = "numberFile(""directory/foo"",123)             = """//trim(numberFile("directory/foo",123))//""""
       output(17) = "numberFile(""directory/foo"",123,"""")          = """//trim(numberFile("directory/foo",123,""))//""""
       output(18) = "numberFile(""directory/foo"",123,""png"")       = """//trim(numberFile("directory/foo",123,"png"))//""""
       output(19) = "numberFile(""directory/foo.bar"",123)         = """//trim(numberFile("directory/foo.bar",123))//""""
       output(20) = "numberFile(""directory/foo.bar"",123,"""")      = """//trim(numberFile("directory/foo.bar",123,""))//""""
       output(21) = "numberFile(""directory/foo.bar"",123,""png"")   = """//trim(numberFile("directory/foo.bar",123,"png"))//""""

       
       
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
       
   end program testLib_Filenames
   
       
    
        
        
    