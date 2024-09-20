
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
add_test(NAME Lib_IntegrateManyBeams_test COMMAND test_Lib_IntegrateManyBeams )
set_tests_properties(Lib_IntegrateManyBeams_test 
    PROPERTIES
    LABELS "EXAMPLE;SHORT"
        #ctest will look for "PASS" in the stdout and consider it a pass if it finds it
    PASS_REGULAR_EXPRESSION "PASS") 
         
