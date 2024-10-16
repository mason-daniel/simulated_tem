
    program test_Lib_IntegrateManyBeams
!---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
!*      
        use iso_fortran_env
        use Lib_Lattices
        use Lib_Elements
        use Lib_XYZFiles
        use Lib_IntegrateManyBeams
        use Lib_CommandLineArguments
#ifdef MPI
        use mpi_f08
#endif
        implicit none




        integer             ::      rank = 0,nProcs = 1


    !---    command line parameter input        
        type(CommandLineArguments)      ::      cla
        
        character(len=256)              ::      filename = ""           !   input filename
!        character(len=256)              ::      xifilename = ""             !   input extinction distance filename
        character(len=8)                ::      latticename = UNKNOWN_LATTICE
        real(kind=real64)               ::      T = 300.0d0             !   temperature (K)
        real(kind=real64)               ::      V = 200.0d0             !   accelerator voltage (kV)
        real(kind=real64),dimension(3)  ::      a0_in = 3.0             !   indicative lattice parameter (A)
        integer                         ::      nPrec = 1               !   precession angles
        real(kind=real64)               ::      theta = 0.003d0         !   precession angle


    !---    physical variables
        type(IntegrateManyBeams)        ::      imb

        integer         ::      ierror
        integer         ::      ii

#ifdef MPI        
        call MPI_INIT(ierror)
        call MPI_COMM_SIZE(MPI_COMM_WORLD, nProcs, ierror)
        call MPI_COMM_RANK(MPI_COMM_WORLD, rank, ierror)
#endif              
        call Lib_IntegrateManyBeams_init_MPI()



!******************************************************************************    
!*    
!*      COMMAND LINE ARGUMENTS
!*    
!*    
!******************************************************************************    
    
        
        
    !---    read command line arguments
        cla = CommandLineArguments_ctor(30)  
         
        call setProgramDescription( cla, "test_Lib_IntegrateManyBeams" )
        call setProgramVersion( cla, "0.0.1" )   
          
        
        call get( cla,"f",filename ,LIB_CLA_REQUIRED,"          input filename" )                                                                  
!        call get( cla,"xi",xifilename ,LIB_CLA_REQUIRED,"         extinction distance input filename" )                                                        
        call get( cla,"lattice",latticename ,LIB_CLA_OPTIONAL,"  lattice type" )                                                        
        call get( cla,"T",T ,LIB_CLA_OPTIONAL,"   temperature (K)" )                                                        
        call get( cla,"V",V ,LIB_CLA_OPTIONAL,"   accelerating voltage (kV)" )                                                        
        ii=0
        call get( cla,"a0",a0_in,ii ,LIB_CLA_REQUIRED,"         lattice parameter(s)" )
        if (ii==1) a0_in(2:3) = a0_in(1)
        
        if (rank==0) call report(cla)
        if (.not. allRequiredArgumentsSet(cla)) call errorExit("")
        if (hasHelpArgument(cla)) call errorExit("")
        call delete(cla)
        



!******************************************************************************    
!*    
!*      WELCOME
!*    
!*    
!******************************************************************************    
        
        if (rank==0) then
            print *,"test_Lib_IntegrateManyBeams"
            print *,"^^^^^^^^^^^^^^^^^^^^^^^^^^^"
            print *,"   input atom positions    """//trim(filename)//""""
 !           print *,"   extinction distances    """//trim(xifilename)//""""
            print *,"   lattice type            """//trim(latticename)//""""
            print *,"   accelerator voltage     ",V," (kV)"
            print *,"   temperature             ",T," (K)"
            print *,"   indicative latt param   ",a0_in
            print *,""  
        end if


!******************************************************************************    
!*    
!*      CONSTRUCTOR
!*    
!*    
!******************************************************************************    
            
        imb = IntegrateManyBeams_ctor( latticename,a0_in,filename,T,V,nPrec = nPrec,theta = theta )
        if (rank==0) call report(imb)



        if (rank==0) print *,""
        call errorExit("PASS")
        if (rank==0) print *,""

               
    contains
!---^^^^^^^^
!        

        subroutine errorExit(message)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
            character(len=*),intent(in),optional             ::      message
#ifdef MPI            
            integer             ::      ierror
#endif
            if (rank==0) then
                if (present(message)) print *,"test_Lib_IntegrateManyBeams::"//trim(message)
            end if
#ifdef MPI            
            call MPI_FINALIZE(ierror)
#endif          
            stop
        end subroutine errorExit
        
        
 
        
    end program test_Lib_IntegrateManyBeams
