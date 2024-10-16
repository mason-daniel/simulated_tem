 



add_executable(run_SimulatedTEM ${CMAKE_SOURCE_DIR}/src/run_SimulatedTEM.f90)
target_link_libraries(run_SimulatedTEM ManyBeam_Utils CSF_Utils PCA_Png CommandLine)

#Install the exec to bin folder of install path
install(TARGETS run_SimulatedTEM
    RUNTIME
    DESTINATION bin)


