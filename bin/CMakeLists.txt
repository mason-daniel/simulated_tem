add_executable(d2bi ${PROJECT_SOURCE_DIR}/src/d2bi.f90)
target_link_libraries(d2bi PCA_Utils PCA_Atoms PCA_DefGrad PCA_D2BI PCA_Png ${LAPACK_LIBRARIES})
 
#Install the exec to bin folder of install path
install(TARGETS d2bi
    RUNTIME
    DESTINATION bin)
    
