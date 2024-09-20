

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
add_library(CSF_Utils STATIC
    Lib_RelativisticElectrons.f90
    Lib_Elements.f90
    Lib_Lattices.f90
    Lib_Linint.f90
    DebyeWallerFactor.f90
    Lib_inputREAL.f90
    Lib_REALfit.f90
    Lib_inputTDS.f90
    Lib_CrystalStructureFactor.f90
    )
set_target_properties(CSF_Utils PROPERTIES 
    Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/include)
target_include_directories(CSF_Utils
    PUBLIC
    $<BUILD_INTERFACE:${CMAKE_BINARY_DIR}/include>
    $<INSTALL_INTERFACE:include>)
target_link_libraries(CSF_Utils CommandLine ${LAPACK_LIBRARIES})
    
    
     





#Copy libraries to the lib directory of
#the install directory

install(TARGETS CSF_Utils
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


 