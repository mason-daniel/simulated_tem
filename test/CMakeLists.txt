
find_package(PNG REQUIRED)
include_directories(${PNG_INCLUDE_DIR})
 
 
# test utils

add_executable(testcallipers testLib_Callipers.f90)
target_link_libraries(testcallipers PCA_Utils )
add_test(NAME calliperstest COMMAND testcallipers)
set_tests_properties(calliperstest
    PROPERTIES
    LABELS "SHORT"
    #ctest will look for "PASS" in the stdout and consider it a pass if it finds it
    PASS_REGULAR_EXPRESSION "PASS") 
     


add_executable(testFilenames testLib_Filenames.f90)
target_link_libraries(testFilenames PCA_Utils)
add_test(NAME Filenamestest COMMAND testFilenames )
set_tests_properties(Filenamestest
    PROPERTIES
    LABELS "SHORT"
    #ctest will look for "PASS" in the stdout and consider it a pass if it finds it
    PASS_REGULAR_EXPRESSION "PASS") 
 
add_executable(testCLA testLib_CommandLineArguments.f90)
target_link_libraries(testCLA PCA_Utils )
add_test(NAME CLAtest COMMAND testCLA -a 12345 -b 1,2,3,4,5 -c 1.2345E+06 -d 1.0e-5,-0.0,0.0,12345.6789,1.0e+5 -e -f "sandwiches, herring")
set_tests_properties(CLAtest
    PROPERTIES
    LABELS "SHORT"
    #ctest will look for "PASS" in the stdout and consider it a pass if it finds it
    PASS_REGULAR_EXPRESSION "PASS") 
     
    


add_executable(testQuicksort testLib_Quicksort.f90)
target_link_libraries(testQuicksort PCA_Utils )
add_test(NAME Quicksorttest COMMAND testQuicksort)
set_tests_properties(Quicksorttest
    PROPERTIES
    LABELS "SHORT"
    #ctest will look for "PASS" in the stdout and consider it a pass if it finds it
    PASS_REGULAR_EXPRESSION "PASS") 

    

add_executable(testRandom testLib_RandomSeed.f90)
target_link_libraries(testRandom PCA_Utils )
add_test(NAME Randomtest COMMAND testRandom)
set_tests_properties(Randomtest
    PROPERTIES
    LABELS "SHORT"
    #ctest will look for "PASS" in the stdout and consider it a pass if it finds it
    PASS_REGULAR_EXPRESSION "PASS") 


add_executable(testTokenizer testNBAX_StringTokenizers.f90)
target_link_libraries(testTokenizer PCA_Utils )
add_test(NAME Tokenizertest COMMAND testTokenizer)
set_tests_properties(Tokenizertest
    PROPERTIES
    LABELS "SHORT"
    #ctest will look for "PASS" in the stdout and consider it a pass if it finds it
    PASS_REGULAR_EXPRESSION "PASS") 

add_executable(testGolden testLib_GoldenSection.f90)
target_link_libraries(testGolden PCA_Utils )
add_test(NAME Goldentest COMMAND testGolden)
set_tests_properties(Goldentest
    PROPERTIES
    LABELS "SHORT"
    #ctest will look for "PASS" in the stdout and consider it a pass if it finds it
    PASS_REGULAR_EXPRESSION "PASS") 

    
    

# test atoms
 
add_executable(testcomplexSS testLib_ComplexSupercells.f90)
target_link_libraries(testcomplexSS PCA_Utils PCA_Atoms )
add_test(NAME complexSStest COMMAND testcomplexSS )
set_tests_properties(complexSStest
    PROPERTIES
    LABELS "SHORT"
    #ctest will look for "PASS" in the stdout and consider it a pass if it finds it
    PASS_REGULAR_EXPRESSION "PASS") 
     
     
add_executable(testLattice testLib_Lattices.f90)
target_link_libraries(testLattice PCA_Utils PCA_Atoms )
add_test(NAME Latticetest COMMAND testLattice )
set_tests_properties(Latticetest
    PROPERTIES
    LABELS "SHORT"
    #ctest will look for "PASS" in the stdout and consider it a pass if it finds it
    PASS_REGULAR_EXPRESSION "PASS") 
     
    
     
add_executable(testLc3d testLib_LinkCell3D.f90)
target_link_libraries(testLc3d PCA_Utils PCA_Atoms )
add_test(NAME Lc3dtest COMMAND testLc3d )
set_tests_properties(Lc3dtest
    PROPERTIES
    LABELS "SHORT"
    #ctest will look for "PASS" in the stdout and consider it a pass if it finds it
    PASS_REGULAR_EXPRESSION "PASS") 
     
    
    
add_executable(testQuaternion testLib_Quaternions.f90)
target_link_libraries(testQuaternion PCA_Utils PCA_Atoms )
add_test(NAME Quaterniontest COMMAND testQuaternion )
set_tests_properties(Quaterniontest
    PROPERTIES
    LABELS "SHORT"
    #ctest will look for "PASS" in the stdout and consider it a pass if it finds it
    PASS_REGULAR_EXPRESSION "PASS") 
     
    
    
add_executable(testRotationMat testLib_RotationMatrices.f90)
target_link_libraries(testRotationMat PCA_Utils PCA_Atoms )
add_test(NAME RotationMattest COMMAND testRotationMat )
set_tests_properties(RotationMattest
    PROPERTIES
    LABELS "SHORT"
    #ctest will look for "PASS" in the stdout and consider it a pass if it finds it
    PASS_REGULAR_EXPRESSION "PASS") 
     
    
add_executable(testSupercell testLib_SimpleSupercells.f90)
target_link_libraries(testSupercell PCA_Utils PCA_Atoms )
add_test(NAME Supercelltest COMMAND testSupercell )
set_tests_properties(Supercelltest
    PROPERTIES
    LABELS "SHORT"
    #ctest will look for "PASS" in the stdout and consider it a pass if it finds it
    PASS_REGULAR_EXPRESSION "PASS") 
      
    
    
    
# test png    
    

#Here is a trivial example test that just checks libpng still works
add_executable(testpng testPng.f90)
target_link_libraries(testpng PCA_Utils PCA_Png ${PNG_LIBRARY} )
add_test(NAME pngtest COMMAND testpng -f ${PROJECT_SOURCE_DIR}/Data/img1.png )
set_tests_properties(pngtest
    PROPERTIES
    LABELS "SHORT"
    #ctest will look for "PASS" in the stdout and consider it a pass if it finds it
    PASS_REGULAR_EXPRESSION "PASS") 
    
    
    
    
# test defgrad 


add_executable(testDefGrad testLib_DeformationGradients.f90)
target_link_libraries(testDefGrad PCA_Utils PCA_DefGrad )
add_test(NAME DefGradtest COMMAND testDefGrad )
set_tests_properties(DefGradtest
    PROPERTIES
    LABELS "SHORT"
    #ctest will look for "PASS" in the stdout and consider it a pass if it finds it
    PASS_REGULAR_EXPRESSION "PASS")     
    

    
    
# add_executable(testDistance testLib_DistancePointToTriangle.f90)
# target_link_libraries(testDistance PCA_Utils PCA_DefGrad )
# add_test(NAME Distancetest COMMAND testDistance )
# set_tests_properties(Distancetest
#     PROPERTIES
#     LABELS "SHORT"
#     #ctest will look for "PASS" in the stdout and consider it a pass if it finds it
#     PASS_REGULAR_EXPRESSION "PASS") 
        
    
    
    
# add_executable(testFitTanh testLib_FitTanhFunction.f90)
# target_link_libraries(testFitTanh PCA_Utils PCA_DefGrad )
# add_test(NAME FitTanhtest COMMAND testFitTanh )
# set_tests_properties(FitTanhtest
#     PROPERTIES
#     LABELS "SHORT"
#     #ctest will look for "PASS" in the stdout and consider it a pass if it finds it
#     PASS_REGULAR_EXPRESSION "PASS") 
        
    
    
    
# add_executable(testGraph testLib_Graphs.f90)
# target_link_libraries(testGraph PCA_Utils PCA_DefGrad )
# add_test(NAME Graphtest COMMAND testGraph )
# set_tests_properties(Graphtest
#     PROPERTIES
#     LABELS "SHORT"
#     #ctest will look for "PASS" in the stdout and consider it a pass if it finds it
#     PASS_REGULAR_EXPRESSION "PASS") 
        
add_executable(testClustering testLib_TopologicallyCorrectClustering.f90)
target_link_libraries(testClustering PCA_Utils PCA_DefGrad )
add_test(NAME Clusteringtest COMMAND testClustering )
set_tests_properties(Clusteringtest
    PROPERTIES
    LABELS "SHORT"
    #ctest will look for "PASS" in the stdout and consider it a pass if it finds it
    PASS_REGULAR_EXPRESSION "PASS") 
        
    
    
# test d2bi
    
add_executable(testRelativistic testLib_RelativisticElectrons.f90)
target_link_libraries(testRelativistic PCA_Utils PCA_D2BI )
add_test(NAME Relativistictest COMMAND testRelativistic )
set_tests_properties(Relativistictest
    PROPERTIES
    LABELS "SHORT"
    #ctest will look for "PASS" in the stdout and consider it a pass if it finds it
    PASS_REGULAR_EXPRESSION "PASS") 
        
    
    
    
    
    
# integration test
     
    
# find_program(BASH_PROGRAM bash)
# if (BASH_PROGRAM)
#   add_test(voidtest ${BASH_PROGRAM} ${CMAKE_CURRENT_SOURCE_DIR}/testVoids.sh ${PROJECT_SOURCE_DIR}/Data/)
#   set_tests_properties(voidtest
#     PROPERTIES
#     LABELS "LONG"
#     #ctest will look for "PASS" in the stdout and consider it a pass if it finds it
#     PASS_REGULAR_EXPRESSION "PASS") 
# endif()


if (BASH_PROGRAM)
  add_test(d2bitest ${BASH_PROGRAM} ${CMAKE_CURRENT_SOURCE_DIR}/testD2BI.sh ${PROJECT_SOURCE_DIR}/Data/)
  set_tests_properties(d2bitest
    PROPERTIES
    LABELS "LONG"
    #ctest will look for "PASS" in the stdout and consider it a pass if it finds it
    PASS_REGULAR_EXPRESSION "PASS") 
endif()


    