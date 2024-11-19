
!   gfortran -ffree-line-length-256 src/Lib_RotationMatrices.f90 src/testLib_RotationMatrices.f90 -o Test/testLib_RotationMatrices.exe

    program testLib_RotationMatrices
!---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
!*      test code to check the correct working of Lib_RotationMatrices
!*      successful operation
!*  
!*      $ ./Test/testLib_RotationMatrices.f90                                                           
!*      construct rotation matrix with Euler angles phi,theta,chi =      1.00000    -0.50000     0.25000
!*      RotationMatrix                                                                                  
!*            0.340808    0.932621   -0.118612                                                          
!*           -0.849176    0.251236   -0.464521                                                          
!*           -0.403423    0.259035    0.877583                                                          
!*                                                                                                      
!*       find z-axis z = R [ 1,2,3 ]                                                                    
!*       complete basis using new z-vector                                                              
!*       det(R)   1.0000000000000000                                                                    
!*       det(R)   1.0000000000000000                                                                    
!*       new rotation matrix                                                                            
!*      RotationMatrix                                                                                  
!*            0.708490   -0.503509    0.494490                                                          
!*            0.705060    0.535319   -0.465106                                                          
!*           -0.030525    0.678169    0.734272                                                          
!*                                                                                                      
!*       done                                                                                           
!*                                                                                                      
!*                                                                                                      
        use iso_fortran_env
        use Lib_ColouredTerminal
        use Lib_RotationMatrices
        implicit none
        
        real(kind=real64),dimension(3,3)            ::      R
        
        real(kind=real64),parameter         ::      PHI = 1.0d0
        real(kind=real64),parameter         ::      THETA = -0.5d0
        real(kind=real64),parameter         ::      CHI = 0.25d0
        
        real(kind=real64),dimension(3)      ::      x,y,z
        
        

        character(len=256),dimension(6)           ::      output
        character(len=*),dimension(6),parameter   ::      output0 = (/  " 0.340808    0.932621   -0.118612  ", &
                                                                        "-0.849176    0.251236   -0.464521  ", &
                                                                        "-0.403423    0.259035    0.877583  ", &
                                                                        " 0.708490   -0.503509    0.494490  ", &
                                                                        " 0.705060    0.535319   -0.465106  ", &
                                                                        "-0.030525    0.678169    0.734272  "  /)
                                                                         
        logical                     ::      ok 
        integer                     ::      ii          
        
        
        write(*,fmt='(a,3f12.5)') "construct rotation matrix with Euler angles phi,theta,chi = ",PHI,THETA,CHI
        R = RotationMatrix_ctor( PHI,THETA,CHI )
        call report(R)
        write(output(1),fmt='(3f16.6)') R(1,1:3)
        write(output(2),fmt='(3f16.6)') R(2,1:3)
        write(output(3),fmt='(3f16.6)') R(3,1:3)
        
        
        print *,""
        print *,"find z-axis z = R [ 1,2,3 ]"
        z = rotateVector( R, (/1.0d0,2.0d0,3.0d0/) )
        print *,"complete basis using new z-vector"        
        z = z/norm2(z)
        x = 0 ; call completeBasis( z,x,y )
        R = RotationMatrix_ctor(x,y,z)
        print *,"new rotation matrix"
        call report(R)
        write(output(4),fmt='(3f16.6)') R(1,1:3)
        write(output(5),fmt='(3f16.6)') R(2,1:3)
        write(output(6),fmt='(3f16.6)') R(3,1:3)
        
        
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
    end program testLib_RotationMatrices