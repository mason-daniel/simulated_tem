

!   gfortran -ffree-line-length-256 -Og -g ${MYF90LIB}/Lib_Quaternions.f90 ${MYF90LIB}/Lib_RotationMatrices.f90 -I${MYF90LIB} src/testCubicDistance.f90 -o testCubicDistance.exe -llapack


    program testCubicDistance
!---^^^^^^^^^^^^^^^^^^^^^^^^^   
        use Lib_Quaternions
        use Lib_RotationMatrices
        use iso_fortran_env
        implicit none
        
        
        include "cubic_rot_q_data.h"
         real(kind=real64),dimension(3,39),parameter       ::      v0 = reshape(            &
                                                                (/  0.0d0,0.0d0,0.0d0,      &
                                                                    0.5d0,0.5d0,0.5d0,      &
                                                                    -.5d0,0.5d0,0.5d0,      &
                                                                    0.5d0,-.5d0,0.5d0,      &
                                                                    0.5d0,0.5d0,-.5d0,      &
                                                                    -.5d0,-.5d0,0.5d0,      &
                                                                    -.5d0,0.5d0,-.5d0,      &
                                                                    0.5d0,-.5d0,-.5d0,      &                                                                                                                                                                        
                                                                    -.5d0,-.5d0,-.5d0,      &
                                                                    1.0d0,0.0d0,0.0d0,      &
                                                                    -1.d0,0.0d0,0.0d0,      &
                                                                    0.0d0,1.0d0,0.0d0,      &
                                                                    0.0d0,-1.d0,0.0d0,      &
                                                                    0.0d0,0.0d0,1.0d0,      &
                                                                    0.0d0,0.0d0,-1.d0,      &
                                                                    1.5d0,0.5d0,0.5d0,      &
                                                                    1.5d0,-.5d0,0.5d0,      &                                                                    
                                                                    1.5d0,0.5d0,-.5d0,      &                                                                    
                                                                    1.5d0,-.5d0,-.5d0,      &                                                                    
                                                                   -1.5d0,0.5d0,0.5d0,      &
                                                                   -1.5d0,-.5d0,0.5d0,      &                                                                    
                                                                   -1.5d0,0.5d0,-.5d0,      &                                                                    
                                                                   -1.5d0,-.5d0,-.5d0,      &                                                                    
                                                                   0.5d0, 1.5d0,0.5d0,      &
                                                                   -.5d0, 1.5d0,0.5d0,      &                                                                    
                                                                   0.5d0, 1.5d0,-.5d0,      &                                                                    
                                                                   -.5d0, 1.5d0,-.5d0,      &                                                                    
                                                                   0.5d0,-1.5d0,0.5d0,      &
                                                                   -.5d0,-1.5d0,0.5d0,      &                                                                    
                                                                   0.5d0,-1.5d0,-.5d0,      &                                                                    
                                                                   -.5d0,-1.5d0,-.5d0,      &                                                                    
                                                                   0.5d0,0.5d0, 1.5d0,      &
                                                                   -.5d0,0.5d0, 1.5d0,      &                                                                    
                                                                   0.5d0,-.5d0, 1.5d0,      &                                                                    
                                                                   -.5d0,-.5d0, 1.5d0,      &                                                                    
                                                                   0.5d0,0.5d0,-1.5d0,      &
                                                                   -.5d0,0.5d0,-1.5d0,      &                                                                    
                                                                   0.5d0,-.5d0,-1.5d0,      &                                                                    
                                                                   -.5d0,-.5d0,-1.5d0   /),(/3,39/) )                                                                    
                                                                    
                                                                    
        
        integer     ::      ii,jj
        real(kind=real64),dimension(3,3)            ::      R1,R2,U
        real(kind=real64),dimension(3,39)           ::      vv
        
        real(kind=real64),dimension(3,3,24)         ::      RU,URU
        real(kind=real64),dimension(24)         ::      fn3
        R1 = rndRotMat()
        call report(R1)
        print *,""
        
        open(unit=600,file="test.xyz",action="write")
        write(600,fmt='(i8)') (24*39)
        write(600,fmt='(a)') "atom position_x position_y position_z charge"
        do ii = 1,24
            U = quaternionToRotMat(cubic_rot_q(ii))
            vv = rotateVector( U,v0 )
            do jj = 1,39    
                write(600,fmt='(a4,100f12.6)') "Du ",vv(1:3,jj),ii*1.0
            end do
!             
            RU(:,:,ii) = matmul( R1,U )
            URU(:,:,ii) = matmul( transpose(U),RU(:,:,ii) )
            print *,"matrix ",ii
            call report(U)
            call report(RU(:,:,ii))
            call report(URU(:,:,ii))
        end do
        close(unit=600)     
        
        
        do ii = 1,24
            do jj = 1,24    
                fn3(jj) = frobeniusNorm3x3(RU(:,:,ii),URU(:,:,jj))                 
            end do
            write(*,fmt='(i6,f8.3,i6,24f8.3)') ii,minval(fn3),minloc(fn3),fn3
        end do
        
        print *,""
        print *,"done"
        print *,""
        
    contains
!---^^^^^^^^
    
        pure function frobeniusNorm3x3(A,B) result(f)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
            real(kind=real64),dimension(3,3),intent(in)     ::      A,B
            real(kind=real64)                               ::      f
            
            f = (A(1,1) - B(1,1))*(A(1,1) - B(1,1))     &
              + (A(1,2) - B(1,2))*(A(1,2) - B(1,2))     &
              + (A(1,3) - B(1,3))*(A(1,3) - B(1,3))     &
              + (A(2,1) - B(2,1))*(A(2,1) - B(2,1))     &
              + (A(2,2) - B(2,2))*(A(2,2) - B(2,2))     &
              + (A(2,3) - B(2,3))*(A(2,3) - B(2,3))     &
              + (A(3,1) - B(3,1))*(A(3,1) - B(3,1))     &
              + (A(3,2) - B(3,2))*(A(3,2) - B(3,2))     &
              + (A(3,3) - B(3,3))*(A(3,3) - B(3,3))
            
            f = sqrt(max(0.0d0,f))
                
            return              
        end function frobeniusNorm3x3       
                                          
            
            
    end program testCubicDistance 
        
        
        
        
        