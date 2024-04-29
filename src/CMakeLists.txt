

#Configure the install
 
find_package(LAPACK REQUIRED)

find_package(PNG REQUIRED)
include_directories(${PNG_INCLUDE_DIR})
 

 
# include(CMakePrintHelpers)
# cmake_print_variables(CMAKE_INSTALL_PREFIX)
    
#Split utility code into its own lib
add_library(PCA_Utils STATIC
    NBAX_StringTokenizers.f90
    NBAX3.f90
    Lib_CommandLineArguments.f90
    Lib_SimpleProgressBar.f90
    Lib_Callipers.f90
    Lib_RandomSeed.f90
    Lib_GoldenSection.f90
    Lib_FactoriseParallel.f90    
    Lib_Filenames.f90
    Lib_ColouredTerminal.f90
#    Lib_Quicksort.f90
    )
set_target_properties(PCA_Utils PROPERTIES 
    Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/include)
target_include_directories(PCA_Utils
    PUBLIC
    $<BUILD_INTERFACE:${CMAKE_BINARY_DIR}/include>
    $<INSTALL_INTERFACE:include>)

    
    
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
    
    
#Split image processing code into its own lib
add_library(PCA_Atoms STATIC
    Lib_RotationMatrices.f90
    Lib_Quaternions.f90
    Lib_Lattices.f90
    Lib_SimpleSupercells.f90
    Lib_ComplexSupercells.f90
    Lib_LinkCell3D.f90
    Lib_XYZFiles.f90
    # Lib_Voxelise.f90
    # Lib_QuinticSplinesEven2.f90
    # Lib_ReadVTK.f90
    # EAM_PotentialTables.f90
    # Lib_ConjugateGradient.f90
    # simpleMD.f90    
)
set_target_properties(PCA_Atoms PROPERTIES 
    Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/include)
target_include_directories(PCA_Atoms
    PUBLIC
    $<BUILD_INTERFACE:${CMAKE_BINARY_DIR}/include>
    $<INSTALL_INTERFACE:include>)
target_link_libraries(PCA_Atoms PCA_Utils ${LAPACK_LIBRARIES})
    
    
    
# #Split PCA processing code into its own lib
# add_library(PCA_DefGrad STATIC

#     # GrainNodes.f90
#     # Lib_ClosestVector.f90
#     # SolveTripleJunction.f90
#     # Lib_DistancePointToTriangle.f90
#     # Lib_DiffractionSpot.omp.f90
#     # fitDeformationGradient.f90
#     # Lib_FitTanhFunction.f90
#     # Lib_Graphs.f90
#     # Lib_FFTW3f.f90
# )
# set_target_properties(PCA_DefGrad PROPERTIES 
#     Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/include)
# target_include_directories(PCA_DefGrad
#     PUBLIC
#     $<BUILD_INTERFACE:${CMAKE_BINARY_DIR}/include>
#     $<INSTALL_INTERFACE:include>)
# target_link_libraries(PCA_DefGrad PCA_Utils PCA_Atoms fftw3)



    
#Split d2bi into its own lib
add_library(PCA_D2BI STATIC
    Lib_RelativisticElectrons.f90    
    Lib_DiffractionConditions.f90
    Lib_DeformationGradients.f90
    Lib_TopologicallyCorrectClustering.f90
    Lib_MarchingCubes.f90
    VoidIsosurfaces.f90
    Lib_DynamicalTwoBeamImaging.f90
    Lib_ReadExtinctionDistances.f90
    Lib_EvenlySpacedPointsInCircle.f90    
#    Lib_Splines1.f90    
    #Lib_DataDoubler.f90
)
set_target_properties(PCA_D2BI PROPERTIES 
    Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/include)
target_include_directories(PCA_D2BI
    PUBLIC
    $<BUILD_INTERFACE:${CMAKE_BINARY_DIR}/include>
    $<INSTALL_INTERFACE:include>)
target_link_libraries(PCA_D2BI PCA_Atoms PCA_Utils)
    


#Copy libraries to the lib directory of
#the install directory

install(TARGETS PCA_Utils PCA_Png PCA_Atoms PCA_D2BI
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


 