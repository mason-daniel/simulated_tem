
   
!   gfortran -ffree-line-length-256 src/Lib_Quaternions.f90 src/testLib_Quaternions.f90 -I./src/ -o Test/testLib_Quaternions.exe -llapack

    
    program testLib_Quaternions
!---^^^^^^^^^^^^^^^^^^^^^^^^^^^
!*
!*      simple test program to show the functioning of quaternion code
!*      
!*      successful output
!*   
!*      $ ./Test/testLib_Quaternions.exe                                                        
!*      Quaternion [w,x,y,z=      0.78575224,      0.23021127,      0.09061727,     -0.56690802]
!*                0.34080762      0.93262072     -0.11861178                                    
!*               -0.84917625      0.25123615     -0.46452136                                    
!*               -0.40342268      0.25903472      0.87758256                                    
!*      rotation by 45 degrees - returns      45.00000000                                       
!*                                                                                              
!*       quaternion rotated back again                                                          
!*      Quaternion [w,x,y,z=      0.78575224,      0.23021127,      0.09061727,     -0.56690802]
!*                                                                                              
!*       quaternion in cubic symmetry group closest to identity                                 
!*      Quaternion [w,x,y,z=     -0.95647524,     -0.22686004,      0.09870786,     -0.15474624]
!*                0.93262072     -0.34080762     -0.11861178                                    
!*                0.25123615      0.84917625     -0.46452136                                    
!*                0.25903472      0.40342268      0.87758256                                    
!*                                                                                              
!*       done                                                                                   
!*                                                                                                                                                                                                      
!*                                                                                              
                                                                                                                
        use iso_fortran_env
        use Lib_ColouredTerminal
        use Lib_Quaternions
        implicit none
        
        type(Quaternion)                    ::      q1,q2 
        real(kind=real64),dimension(3,3)    ::      R1,R2 
        real(kind=real64)       ::      theta = 45.0d0  * ( 3.141592653590d0/180.0d0 )  !   45 degrees

       
        
        

        character(len=256),dimension(6)            ::      output
        character(len=*),dimension(6),parameter   ::      output0 = (/  "0.78575224 0.23021127 0.09061727 -0.56690802   ", &
                                                                        "45.00000000                                    ", &
                                                                        " 0.78575224 0.23021127 0.09061727 -0.56690802  ", &
                                                                        "0.93262072 -0.34080762 -0.11861178             ", &
                                                                        "0.25123615 0.84917625 -0.46452136              ", &
                                                                        "0.25903472 0.40342268 0.87758256               "  /)
                                                                        
                                                                       
                                                                        
        logical                     ::      ok 
        integer                     ::      ii                                                             
                                                                   
                               
        
        
        
        
        
        
        
        
    !---    construct an arbitrary rotation matrix. 
    !       Note fortran uses column-major order, so the actual matrix "looks like" the transpose of how its written here.
        R1(1:3,1:3) = reshape(  (/ 0.340807622914d0,-0.849176251043d0,-0.403422680111d0,     &                               
                                   0.932620721761d0,0.251236146093d0,0.259034724000d0,      &    
                                  -0.118611776418d0,-0.464521359639d0,0.877582561890d0 /) , (/3,3/) )    
          
                                                            
        q1 = Quaternion_ctor( R1 )
        
    !---    test conversion back-forth between rot mat representation        
        call report( q1,asRotMat = .true. )
        write(output(1) ,fmt='(4f16.8)') getComponents(q1)
        
    !---    rotate quaternion by a recognisable angle
                             
        R2(1:3,1:3) = reshape(  (/ cos(theta),-sin(theta),0.0d0, sin(theta),cos(theta),0.0d0, 0.0d0,0.0d0,1.0d0 /) , (/3,3/) )                   
        q2 = Quaternion_ctor( R2 )  * q1
        write (*,fmt='(a,f16.8)') "rotation by 45 degrees - returns ",quaternionDistanceInDegrees( distance(q1,q2) )
    
        write(output(2) ,fmt='(4f16.8)') quaternionDistanceInDegrees( distance(q1,q2) )
        
    !---    check we can rotate back again    
        print *,""
        print *,"quaternion rotated back again"   
        q2 = inverse( Quaternion_ctor( R2 ) ) * q2
        call report( q2 )
        write(output(3) ,fmt='(4f16.8)') getComponents(q2) 
       
    !---    find a symmetry related rotation matrix
        print *,""
        print *,"quaternion in cubic symmetry group closest to identity"
        call imposeSymmetry( q1,QUATERNION_SYMMETRY_CUBIC )
        call report( q1,asRotMat = .true. )
        R1 = quaternionToRotMat(q1)
        write(output(4) ,fmt='(4f16.8)') R1(1,1:3)
        write(output(5) ,fmt='(4f16.8)') R1(2,1:3)
        write(output(6) ,fmt='(4f16.8)') R1(3,1:3)
        
        
        
        
        
    !---    here is the simple test: does the output look like the stored output?
        ok = .true.
        do ii = 1,size(output)            
            if ( trim(cutSpaces(output(ii))) == trim(cutSpaces(output0(ii))) ) then
                write (*,fmt='(a)') trim(cutSpaces(output(ii)))
            else
                ok = .false.
                write (*,fmt='(a)') trim(cutSpaces(output0(ii)))//"    "//colour(RED,trim(cutSpaces(output(ii))))
            end if
        end do   
        
       
    !---    output the result "PASS" or "FAIL"    
        if (ok) then
            print *,colour(LIGHT_GREEN,"PASS")
        else
            print *,colour(RED,"FAIL")
        end if
        
        
        print *,""
        print *,"done"
        print *,""
        
        
    end program testLib_Quaternions
       
       
       
       
       
       
       
       
       
       
       