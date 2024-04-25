
    module Lib_DeformationMatrixToStrainRotation
!---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
!*      given the general deformation matrix M = (1 + e)R
!*      recover the strain e and rotation matrix R
!*
!*      This is done simply by constructing the Left Cauchy-Green Deformation Tensor
!*         M M^T = (1+e)(1+e)^T 
!*      You can then find the eigendecomposition M M^T = U D U^T
!*      and then (1+e) = U sqrt(D) U^T
!*      
!*      So we have a solution ( if it exists )  
!*          e = U sqrt(D) U^T - 1
!*          R = (1 + e)^-1 M


        use iso_fortran_env
        implicit none
        private
        
        external    ::      DSYEV
        
        public      ::      DeformationMatrixToStrainRotation

    contains
!---^^^^^^^^

        subroutine DeformationMatrixToStrainRotation(M,e,R)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
            real(kind=real64),dimension(3,3),intent(in)     ::      M
            real(kind=real64),dimension(3,3),intent(out)    ::      e,R
            
            real(kind=real64),dimension(3,3)        ::      MMT,UU,i1pluse
            real(kind=real64),dimension(3)          ::      dd,idd
            real(kind=real64),dimension(3*66)       ::      work
            integer                                 ::      ii
            
        !---    compute symmetrised matrix from M            
            MMT(1:3,1) = M(1:3,1)*M(1,1) + M(1:3,2)*M(1,2) + M(1:3,3)*M(1,3)
            MMT(1:3,2) = M(1:3,1)*M(2,1) + M(1:3,2)*M(2,2) + M(1:3,3)*M(2,3)
            MMT(1:3,3) = M(1:3,1)*M(3,1) + M(1:3,2)*M(3,2) + M(1:3,3)*M(3,3)
            
        !---    find eivals and eivecs using Lapack
            call DSYEV( "V","U",3,MMT,3,dd,work,size(work),ii)         
            if (ii/=0) then
                print *,"Lib_DeformationMatrixToStrainRotation::DeformationMatrixToStrainRotation() error - DSYEV failed with error code ",ii
                e = 0
                R = reshape( (/1,0,0 , 0,1,0 , 0,0,1/),(/3,3/) )
                return
            end if
            UU = MMT        !   this will be optimised out. The point is that MMT now stores normalised eivecs.            
            
            
        !---    find the square roots of the eigenvalues. If any are negative, throw a warning and return zero.
            do ii = 1,3
                if (dd(ii)>0) then
                    dd(ii) = sqrt(dd(ii))
                    idd(ii) = 1/dd(ii)
                else if (dd(ii)>-1.0d-8) then
                    !   actually this is probably just a machine precision rounding error. Fix quietly.                    
                    dd(ii) = 0.0d0
                    idd(ii) = 0.0d0
                else
                    print *,"Lib_DeformationMatrixToStrainRotation::DeformationMatrixToStrainRotation() warning - negative eigenvalue ",dd(ii)," setting zero"
                    dd(ii) = 0.0d0
                    idd(ii) = 0.0d0                    
                end if
            end do
                    
        !---    can now construct the matrix (1+e)^-1
            i1pluse(1:3,1) = UU(1:3,1)*idd(1)*UU(1,1) + UU(1:3,2)*idd(2)*UU(1,2) + UU(1:3,3)*idd(3)*UU(1,3)
            i1pluse(1:3,2) = UU(1:3,1)*idd(1)*UU(2,1) + UU(1:3,2)*idd(2)*UU(2,2) + UU(1:3,3)*idd(3)*UU(2,3)
            i1pluse(1:3,3) = UU(1:3,1)*idd(1)*UU(3,1) + UU(1:3,2)*idd(2)*UU(3,2) + UU(1:3,3)*idd(3)*UU(3,3)
            
        !---    ... and use this to construct R
            R(1:3,1) = i1pluse(1:3,1)*M(1,1) + i1pluse(1:3,2)*M(2,1) + i1pluse(1:3,3)*M(3,1)             
            R(1:3,2) = i1pluse(1:3,1)*M(1,2) + i1pluse(1:3,2)*M(2,2) + i1pluse(1:3,3)*M(3,2)             
            R(1:3,3) = i1pluse(1:3,1)*M(1,3) + i1pluse(1:3,2)*M(2,3) + i1pluse(1:3,3)*M(3,3)    
            
        !---    finally construct the matrix e
            e(1:3,1) = UU(1:3,1)*dd(1)*UU(1,1) + UU(1:3,2)*dd(2)*UU(1,2) + UU(1:3,3)*dd(3)*UU(1,3)
            e(1:3,2) = UU(1:3,1)*dd(1)*UU(2,1) + UU(1:3,2)*dd(2)*UU(2,2) + UU(1:3,3)*dd(3)*UU(2,3)
            e(1:3,3) = UU(1:3,1)*dd(1)*UU(3,1) + UU(1:3,2)*dd(2)*UU(3,2) + UU(1:3,3)*dd(3)*UU(3,3)
            
            e(1,1)   = e(1,1) - 1.0d0
            e(2,2)   = e(2,2) - 1.0d0
            e(3,3)   = e(3,3) - 1.0d0
            
            return
        end subroutine DeformationMatrixToStrainRotation
        
    end module Lib_DeformationMatrixToStrainRotation
    
!!   gfortran -ffree-line-length-256 ${MYF90LIB}/Lib_RotationMatrices.f90 ${MYF90LIB}/Lib_DeformationMatrixToStrainRotation.f90 -llapack    
!    
!    program testLib_DeformationMatrixToStrainRotation
!!---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
!        use Lib_DeformationMatrixToStrainRotation
!        use Lib_RotationMatrices
!        use iso_fortran_env
!        implicit none
!        
!        real(kind=real64),dimension(3,3)        ::      R,e,I
!        real(kind=real64),dimension(3,3)        ::      M
!        
!        
!        I = reshape( (/1,0,0 , 0,1,0 , 0,0,1/) , (/3,3/) )
!        
!        e = reshape( (/0.1,0.0,-0.01 , 0.0,-0.1,0.02 , -0.01,0.02,0.2/) , (/3,3/) )
!        R = RotationMatrix_ctor( 0.1d0,0.2d0,-0.3d0 ) 
!        
!        
!        
!        M = I + e
!        M = matmul( M,R )
!        print *,"before"
!        write (*,fmt='(a,3f16.8)') "strain ",e(1,:)
!        write (*,fmt='(a,3f16.8)') "       ",e(2,:)
!        write (*,fmt='(a,3f16.8)') "       ",e(3,:)
!        write (*,fmt='(a,3f16.8)') "rot    ",R(1,:)
!        write (*,fmt='(a,3f16.8)') "       ",R(2,:)
!        write (*,fmt='(a,3f16.8)') "       ",R(3,:)
!        
!        e = 0 ; R = I
!
!        print *,"after"        
!        call DeformationMatrixToStrainRotation(M,e,R)
!        write (*,fmt='(a,3f16.8)') "strain ",e(1,:)
!        write (*,fmt='(a,3f16.8)') "       ",e(2,:)
!        write (*,fmt='(a,3f16.8)') "       ",e(3,:)
!        write (*,fmt='(a,3f16.8)') "rot    ",R(1,:)
!        write (*,fmt='(a,3f16.8)') "       ",R(2,:)
!        write (*,fmt='(a,3f16.8)') "       ",R(3,:)
!                    
!        print *,""
!        print *,"done"
!        print *,""
!        
!     end program testLib_DeformationMatrixToStrainRotation
!        