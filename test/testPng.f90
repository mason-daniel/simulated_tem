
    program testPng
!---^^^^^^^^^^^^^^^^
!*
!*      Daniel Mason
!*      (c) UKAEA          
!*      June 2023
!*

!*      unit test for Lib_Png

        use iso_fortran_env
        use Lib_ColouredTerminal
        use Lib_CommandLineArguments
        use Lib_Png
        implicit none
         
         
        type(CommandLineArguments)      ::      cla
        character(len=256)              ::      filename = "data/img1.png"             
        logical                         ::      quiet = .false.
        
        real(kind=real64),dimension(:,:),allocatable        ::      img
        integer                                             ::      Nx,Ny
        logical                                             ::      ok
                                       
    !---    check command line args                                                                                                                                                                                                                      
        cla = CommandLineArguments_ctor(10)
        call setProgramDescription( cla, "testPng.exe" )
        call setProgramVersion( cla, "1.0.0" )
        call get( cla,"f",filename,LIB_CLA_REQUIRED," filename" )
        call get( cla,"q",quiet,LIB_CLA_OPTIONAL," quiet mode" )
        call report(cla)              !   gives full output iff there is an error, or "-h" is a cla
        if (.not. allRequiredArgumentsSet(cla)) stop
        if (hasHelpArgument(cla)) stop
        call delete(cla)
        
        
        
        
    !---    read input .png        
        if (.not. quiet) print *,"reading input """//trim(filename)//""""
        inquire(file=trim(filename),exist=ok)
        if (.not. ok) then  
            print *,"file not found error """//trim(filename)//""" "//colour(RED,"FAIL")
            stop
        end if
        call readPng(filename,img)
        
        
        
    !---    here is the test - did I return an allocated array?        
        ok = allocated(img)
        if (ok .and. .not. quiet) then
            Nx = size(img,dim=1)
            Ny = size(img,dim=2)
            print *,"image size ",Nx,"x",Ny
            print *,""
        end if        
        
        
        
    !---    output the result "PASS" or "FAIL"    
        if (ok) then
            print *,colour(LIGHT_GREEN,"PASS")
        else
            print *,colour(RED,"FAIL")
        end if
        
         
    !---    bye bye   
        if (.not. quiet) then
            print *,""
            print *,"done"
            print *,""
        end if
        
    end program testPng