
add_executable(test_Lib_AtomSpace test_Lib_AtomSpace.f90)
target_link_libraries(test_Lib_AtomSpace ManyBeam_Utils CommandLine )
add_test(NAME Lib_AtomSpace_test COMMAND test_Lib_AtomSpace )
set_tests_properties(Lib_AtomSpace_test 
    PROPERTIES
    LABELS "EXAMPLE;SHORT"
    #ctest will look for "PASS" in the stdout and consider it a pass if it finds it
    PASS_REGULAR_EXPRESSION "PASS") 
     
    
add_executable(test_Lib_ComputePhaseFactor test_Lib_ComputePhaseFactor.f90)
target_link_libraries(test_Lib_ComputePhaseFactor ManyBeam_Utils CommandLine )
add_test(NAME Lib_ComputePhaseFactor_test COMMAND test_Lib_ComputePhaseFactor )
set_tests_properties(Lib_ComputePhaseFactor_test 
    PROPERTIES
    LABELS "EXAMPLE;SHORT"
    #ctest will look for "PASS" in the stdout and consider it a pass if it finds it
    PASS_REGULAR_EXPRESSION "PASS") 
     
    
add_executable(test_Lib_FibonacciSphere test_Lib_FibonacciSphere.f90)
target_link_libraries(test_Lib_FibonacciSphere ManyBeam_Utils CommandLine )
add_test(NAME Lib_FibonacciSphere_test COMMAND test_Lib_FibonacciSphere 1000 )
set_tests_properties(Lib_FibonacciSphere_test
    PROPERTIES
    LABELS "EXAMPLE;SHORT"
    #ctest will look for "PASS" in the stdout and consider it a pass if it finds it
    PASS_REGULAR_EXPRESSION "PASS") 
     
    
 
add_executable(test_Lib_ImagingSpace test_Lib_ImagingSpace.F90)
target_link_libraries(test_Lib_ImagingSpace ManyBeam_Utils CommandLine )
if ( CMAKE_SERIAL )
    add_test(NAME Lib_ImagingSpace_test 
        COMMAND test_Lib_ImagingSpace
    )
else()
    add_test(NAME Lib_ImagingSpace_test
        COMMAND mpirun -n ${NPROCS_TEST} test_Lib_ImagingSpace
    )
endif()
set_tests_properties(Lib_ImagingSpace_test 
    PROPERTIES
    LABELS "EXAMPLE;SHORT"
    #ctest will look for "PASS" in the stdout and consider it a pass if it finds it
    PASS_REGULAR_EXPRESSION "PASS") 
     

add_executable(test_Lib_Gvectors test_Lib_Gvectors.f90)
target_link_libraries(test_Lib_Gvectors ManyBeam_Utils CommandLine )
add_test(NAME Lib_Gvectors_test COMMAND test_Lib_Gvectors )
set_tests_properties(Lib_Gvectors_test 
    PROPERTIES
    LABELS "EXAMPLE;SHORT"
        #ctest will look for "PASS" in the stdout and consider it a pass if it finds it
    PASS_REGULAR_EXPRESSION "PASS") 

add_executable(test_Lib_IntegrateManyBeams test_Lib_IntegrateManyBeams.f90)
target_link_libraries(test_Lib_IntegrateManyBeams ManyBeam_Utils CommandLine )
add_test(NAME Lib_IntegrateManyBeams_test COMMAND test_Lib_IntegrateManyBeams -f ${PROJECT_SOURCE_DIR}/data/test.i111_loop_W.lammps -a0 3.1652 )
set_tests_properties(Lib_IntegrateManyBeams_test 
    PROPERTIES
    LABELS "EXAMPLE;SHORT"
        #ctest will look for "PASS" in the stdout and consider it a pass if it finds it
    PASS_REGULAR_EXPRESSION "PASS") 
         

    add_executable(test_Lib_Elements test_Lib_Elements.f90)
    target_link_libraries(test_Lib_Elements CSF_Utils)
    add_test(NAME Lib_Elements_test COMMAND test_Lib_Elements)
    set_tests_properties(Lib_Elements_test
        PROPERTIES
        LABELS "SHORT"
        #ctest will look for "PASS" in the stdout and consider it a pass if it finds it
        PASS_REGULAR_EXPRESSION "PASS") 
         
    add_executable(test_Lib_inputReal test_Lib_inputReal.f90)
    target_link_libraries(test_Lib_inputReal CSF_Utils)
    add_test(NAME Lib_inputReal_test COMMAND test_Lib_inputReal)
    set_tests_properties(Lib_inputReal_test
        PROPERTIES
        LABELS "SHORT"
        #ctest will look for "PASS" in the stdout and consider it a pass if it finds it
        PASS_REGULAR_EXPRESSION "PASS") 
         
    add_executable(test_Lib_REALfit test_Lib_REALfit.f90)
    target_link_libraries(test_Lib_REALfit CSF_Utils)
    add_test(NAME Lib_REALfit_test COMMAND test_Lib_REALfit)
    set_tests_properties(Lib_REALfit_test
        PROPERTIES
        LABELS "SHORT"
        #ctest will look for "PASS" in the stdout and consider it a pass if it finds it
        PASS_REGULAR_EXPRESSION "PASS") 
        
        
    add_executable(test_Lib_InputTDS test_Lib_InputTDS.f90)
    target_link_libraries(test_Lib_InputTDS CSF_Utils)
    add_test(NAME Lib_InputTDS_test COMMAND test_Lib_InputTDS)
    set_tests_properties(Lib_InputTDS_test
        PROPERTIES
        LABELS "SHORT"
        #ctest will look for "PASS" in the stdout and consider it a pass if it finds it
        PASS_REGULAR_EXPRESSION "PASS") 
    
        
    add_executable(test_Lib_TDSfit test_Lib_TDSfit.f90)
    target_link_libraries(test_Lib_TDSfit CSF_Utils)
    add_test(NAME Lib_TDSfit_test COMMAND test_Lib_TDSfit)
    set_tests_properties(Lib_TDSfit_test
        PROPERTIES
        LABELS "SHORT"
        #ctest will look for "PASS" in the stdout and consider it a pass if it finds it
        PASS_REGULAR_EXPRESSION "PASS") 
    
        
    add_executable(test_Lib_CrystalStructureFactor test_Lib_CrystalStructureFactor.f90)
    target_link_libraries(test_Lib_CrystalStructureFactor CSF_Utils)
    add_test(NAME Lib_CrystalStructureFactor_test COMMAND test_Lib_CrystalStructureFactor -hkl 2,0,0 )
    set_tests_properties(Lib_CrystalStructureFactor_test
        PROPERTIES
        LABELS "SHORT"
        #ctest will look for "PASS" in the stdout and consider it a pass if it finds it
        PASS_REGULAR_EXPRESSION "PASS") 
        
     