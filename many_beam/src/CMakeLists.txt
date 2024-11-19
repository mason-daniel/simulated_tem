

#Configure the install
 
find_package(LAPACK REQUIRED)
 

find_package(PNG REQUIRED)
include_directories(${PNG_INCLUDE_DIR})
 

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
    Lib_RandomSeed.f90
    Lib_GoldenSection.f90
    Lib_RelativisticElectrons.f90
    Lib_Elements.f90
    Lib_Lattices.f90
    Lib_Linint.f90    
    Lib_FitDoyleTurner.f90
    Lib_inputREAL.f90
    Lib_REALfit.f90
    Lib_inputTDS.f90
    Lib_TDSfit.f90
    Lib_CrystalStructureFactor.f90
    )
set_target_properties(CSF_Utils PROPERTIES 
    Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/include)
target_include_directories(CSF_Utils
    PUBLIC
    $<BUILD_INTERFACE:${CMAKE_BINARY_DIR}/include>
    $<INSTALL_INTERFACE:include>)
target_link_libraries(CSF_Utils CommandLine ${LAPACK_LIBRARIES})
    

#Split png handling code into its own lib
add_library(PCA_Png STATIC
    Lib_Png.f90
    Lib_ColourScale.f90 
    Lib_Greyscale.c
)
set_target_properties(PCA_Png PROPERTIES 
    Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/include)
target_include_directories(PCA_Png
    PUBLIC
    $<BUILD_INTERFACE:${CMAKE_BINARY_DIR}/include>
    $<INSTALL_INTERFACE:include>)
target_link_libraries(PCA_Png ${PNG_LIBRARY})
    
    
    
#Split utility code into its own lib
add_library(ManyBeam_Utils STATIC
    Lib_XYZFiles32.f90
    Lib_FactoriseParallel.f90
    Lib_RotationMatrices.f90
    Lib_Quaternions.f90
    Lib_FibonacciSphere.f90
    Lib_ImagingSpace.F90
    Lib_ComputePhaseFactor.f90
    Lib_AtomSpace.f90
    Lib_Lattices.f90
    Lib_Gvectors.f90
    Lib_DeformationGradients.f90
    Lib_ReadExtinctionDistances.f90
    Lib_RelativisticElectrons.f90
    Lib_ManyBeam.f90 
    Lib_RK4.f90
    Lib_IntegrateManyBeams.f90
    )
set_target_properties(ManyBeam_Utils PROPERTIES 
    Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/include)
target_include_directories(ManyBeam_Utils
    PUBLIC
    $<BUILD_INTERFACE:${CMAKE_BINARY_DIR}/include>
    $<INSTALL_INTERFACE:include>)
target_link_libraries(ManyBeam_Utils CSF_Utils CommandLine ${LAPACK_LIBRARIES})
    
    
     





#Copy libraries to the lib directory of
#the install directory

install(TARGETS CSF_Utils
    EXPORT ${PROJECT_NAME}_Targets
    ARCHIVE
    DESTINATION lib
    INCLUDES 
    DESTINATION
    include)
    

install(TARGETS ManyBeam_Utils
    EXPORT ${PROJECT_NAME}_Targets
    ARCHIVE
    DESTINATION lib
    INCLUDES 
    DESTINATION
    include)
    

install(TARGETS PCA_Png
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


 