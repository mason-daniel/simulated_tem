     
    program testLib_Voxelise
!---^^^^^^^^^^^^^^^^^^^^^^^^ 
!*      simple program to test functioning of Lib_Voxelise
!*      
!*          successful result        
!*
!*           writing data for validation to "Output/testTopologicallyCorrectClustering.xyz"
!*           Lib_topologicallyCorrectClustering::topologicallyCorrectCluster() info - quick forwards cluster pass
!*           Lib_topologicallyCorrectClustering::topologicallyCorrectCluster() info - topologically correct cluster pass
!*           number of clusters          191
!*           clusters over size 10
!*           cluster            1  size         2328
!*           cluster            7  size           17
!*           cluster           45  size           12
!*           cluster           62  size           14
!*           total in clusters         2668  cf f>=0         2668
!*          
!*           done
!*          
!*                      
!*


        use iso_fortran_env
        use Lib_Voxelise
        implicit none

        
        
        print *,""
        print *,"done"
        print *,""  
        
    end program testLib_Voxelise