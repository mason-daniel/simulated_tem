

#Configure the install
 
find_package(LAPACK REQUIRED)
 

#Split the command line parsing code into its own lib
#And do the same thing
add_library(CommandLine STATIC
    Lib_CommandLineArguments.f90
    Lib_Filenames.f90
    Lib_Callipers.f90
    Lib_ColouredTerminal.f90
    Lib_SimpleProgressBar.f90
    NBAX_StringTokenizers.f90
    NBAX3.f90 )
set_target_properties(CommandLine PROPERTIES 
    Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/include)
target_include_directories(CommandLine
    PUBLIC
    $<BUILD_INTERFACE:${CMAKE_BINARY_DIR}/include>
    $<INSTALL_INTERFACE:include>)

    
#Split utility code into its own lib
add_library(ManyBeam_Utils STATIC
    Lib_RandomSeed.f90
    Lib_XYZFiles.f90
    Lib_FactoriseParallel.f90
    Lib_RotationMatrices.f90
    Lib_Quaternions.f90
    Lib_ImagingSpace.F90
    Lib_ComputePhaseFactor.f90
    Lib_AtomSpace.f90
    Lib_Lattices.f90
    Lib_Gvectors.f90
    Lib_DeformationGradients.f90
    Lib_ReadExtinctionDistances.f90
    Lib_RelativisticElectrons.f90
    Lib_ManyBeam.f90 
    Lib_IntegrateManyBeams.f90
    )
set_target_properties(ManyBeam_Utils PROPERTIES 
    Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/include)
target_include_directories(ManyBeam_Utils
    PUBLIC
    $<BUILD_INTERFACE:${CMAKE_BINARY_DIR}/include>
    $<INSTALL_INTERFACE:include>)
target_link_libraries(ManyBeam_Utils ${LAPACK_LIBRARIES})
    
    
     





#Copy libraries to the lib directory of
#the install directory

install(TARGETS ManyBeam_Utils
    EXPORT ${PROJECT_NAME}_Targets
    ARCHIVE
    DESTINATION lib
    INCLUDES 
    DESTINATION
    include)
    
    
install(TARGETS CommandLine
    EXPORT ${PROJECT_NAME}_Targets
    ARCHIVE
    DESTINATION lib
    INCLUDES 
    DESTINATION
    include)
    
#...and the fortran module to the include directory
install(DIRECTORY
    ${CMAKE_BINARY_DIR}/include
    DESTINATION
    ${CMAKE_INSTALL_PREFIX})


 