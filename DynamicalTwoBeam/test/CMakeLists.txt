 
# test utils

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
 
add_executable(testcallipers testLib_Callipers.f90)
target_link_libraries(testcallipers CSF_Utils )
add_test(NAME calliperstest COMMAND testcallipers)
set_tests_properties(calliperstest
    PROPERTIES
    LABELS "SHORT"
    #ctest will look for "PASS" in the stdout and consider it a pass if it finds it
    PASS_REGULAR_EXPRESSION "PASS") 
     


add_executable(testFilenames testLib_Filenames.f90)
target_link_libraries(testFilenames CSF_Utils)
add_test(NAME Filenamestest COMMAND testFilenames )
set_tests_properties(Filenamestest
    PROPERTIES
    LABELS "SHORT"
    #ctest will look for "PASS" in the stdout and consider it a pass if it finds it
    PASS_REGULAR_EXPRESSION "PASS") 
 
add_executable(testCLA testLib_CommandLineArguments.f90)
target_link_libraries(testCLA CSF_Utils )
add_test(NAME CLAtest COMMAND testCLA -a 12345 -b 1,2,3,4,5 -c 1.2345E+06 -d 1.0e-5,-0.0,0.0,12345.6789,1.0e+5 -e -f "sandwiches, herring")
set_tests_properties(CLAtest
    PROPERTIES
    LABELS "SHORT"
    #ctest will look for "PASS" in the stdout and consider it a pass if it finds it
    PASS_REGULAR_EXPRESSION "PASS") 
     

add_executable(testTokenizer testNBAX_StringTokenizers.f90)
target_link_libraries(testTokenizer CSF_Utils )
add_test(NAME Tokenizertest COMMAND testTokenizer)
set_tests_properties(Tokenizertest
    PROPERTIES
    LABELS "SHORT"
    #ctest will look for "PASS" in the stdout and consider it a pass if it finds it
    PASS_REGULAR_EXPRESSION "PASS")     
 
add_executable(testLattice testLib_Lattices.f90)
target_link_libraries(testLattice CSF_Utils )
add_test(NAME Latticetest COMMAND testLattice )
set_tests_properties(Latticetest
    PROPERTIES
    LABELS "SHORT"
    #ctest will look for "PASS" in the stdout and consider it a pass if it finds it
    PASS_REGULAR_EXPRESSION "PASS") 
     

add_executable(testRelativistic testLib_RelativisticElectrons.f90)
target_link_libraries(testRelativistic CSF_Utils )
add_test(NAME Relativistictest COMMAND testRelativistic )
set_tests_properties(Relativistictest
    PROPERTIES
    LABELS "SHORT"
    #ctest will look for "PASS" in the stdout and consider it a pass if it finds it
    PASS_REGULAR_EXPRESSION "PASS") 
        
    
        