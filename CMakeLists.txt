#Change this to whatever cmake version you have
cmake_minimum_required(VERSION 3.14) 

#Change the name and the version as you like
project(SIMULATED_TEM
    VERSION 0.1.0
    #Add C, CXX etc if interface moves to the same project
    LANGUAGES C Fortran)

     
find_package(MPI)
include_directories(SYSTEM ${MPI_INCLUDE_PATH})    

#find_package(OpenMP)

include(CMakePrintHelpers)
cmake_print_variables(PROJECT_NAME)
cmake_print_variables(CMAKE_INSTALL_PREFIX)
cmake_print_variables(PROJECT_SOURCE_DIR)
cmake_print_variables(CMAKE_Fortran_COMPILER)
 


# check that parallel version should be built
# I'm not sure this is working exactly as expected: 
#   I'd think that "MPI_Fortran_HAVE_F08_MODULE" means "I can find the mpi_f08 module" 
#   
if ( CMAKE_SERIAL )
    message ("forcing serial compilation")
else()
    message ("parallel compilation")
    #   Grab MPI?
    find_package(MPI)
    include_directories(SYSTEM ${MPI_INCLUDE_PATH})      
    if (MPI_Fortran_HAVE_F08_MODULE)
        message(STATUS "MPI F08-compliant compiler found, will add -DMPI flag to compilation") 
        add_definitions(-DMPI)
    else()
        message(WARNING "MPI F08-compliant compiler not found")    
    endif()
endif()


        
# if (OPENMP_FOUND)
#     set (CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${OpenMP_Fortran_FLAGS}")
#     set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
# endif()



#Add the option to turn tests off i.e. cmake -DENABLE_TESTS=OFF ....
option(ENABLE_TESTS "Build unit tests" ON)


#Don't know if the flags work for ifort
if(CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
    #Set some compiler flags
    message(STATUS "setting compile/link options for GNU fortran compiler")    
    set (CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -ffree-line-length-256")
    add_compile_options(
        #Added regardless of build options
        "-I/usr/include;-cpp;-Wall"
        #For CMAKE_BUILD_TYPE=Release
        "$<$<CONFIG:Release>:-O2;-march=native>"
        #For general profiling, prevents inlining, and adds debug info but still optimised
        "$<$<CONFIG:Perf>:-O2;-fno-inline;-fno-omit-frame-pointer;-g>"
        #For gprof, prevents inlining, and adds debug info but still optimised
        "$<$<CONFIG:Gprof>:-O0;-fno-inline;-g;>"
        )
    add_link_options( 
         "$<$<CONFIG:Debug>:-fbacktrace;-fbounds-check;-ffpe-trap=invalid,zero,overflow;-g;-O0;-Wsurprising>"
         "$<$<CONFIG:Gprof>:-p;-pg>"
         )
elseif(CMAKE_Fortran_COMPILER_ID MATCHES "Intel")    
    #Add flags for intel (am assuming the old and LLVM compilers have the same flags 
    #I think can use STREQUAL "Intel" and STREQUAL "IntelLLVM" to distinguish
    #if needed
    message(STATUS "setting compile/link options for Intel fortran compiler")    
    add_compile_options(
        #Added regardless of build options
        "-cpp"
        )
endif()

#Normally you want to add -DCMAKE_INSTALL_PREFIX=/path/to/install/location
#but this will default to the build directory for convenience
if (CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
    set(CMAKE_INSTALL_PREFIX ${CMAKE_BINARY_DIR} CACHE PATH "..." FORCE)
endif()

add_subdirectory(src)
add_subdirectory(bin)

if (ENABLE_TESTS)
    enable_testing()
    add_subdirectory(test)
endif()


#This stuff is boiler plate to export the details of
#the project to the install directory
#once done another cmake project can
#easily link to the library
include(CMakePackageConfigHelpers)
write_basic_package_version_file(
    ${PROJECT_NAME}ConfigVersion.cmake
    VERSION ${PACKAGE_VERSION}
    COMPATIBILITY AnyNewerVersion)

install(FILES 
    ${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}ConfigVersion.cmake
    DESTINATION lib/cmake)

install(EXPORT ${PROJECT_NAME}_Targets
    FILE ${PROJECT_NAME}Config.cmake
    NAMESPACE
    ${PROJECT_NAME}::
    DESTINATION
    lib/cmake)
export(TARGETS PCA_Utils PCA_Png PCA_Atoms FILE ${PROJECT_NAME}Config.cmake)
