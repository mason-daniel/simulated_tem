
    module distributeQuaternions
!---^^^^^^^^^^^^^^^^^^^^^^^^^^^^
!*      distribute quaternions according to Mitchell's best candidate algorithm
!*      imposing cubic symmetry
        use Lib_MitchellsBestCandidate
        use Lib_RotationMatrices
        use Lib_Quaternions       
        use iso_fortran_env
        implicit none
        private
        
        include "cubic_rot_q_data.h"
        
       ! integer,public,parameter        ::      SYMMETRY_NONE = 0
       ! integer,public,parameter        ::      SYMMETRY_CUBIC = 1
        
        integer,public                  ::      distributeQuaternions_sym = SYMMETRY_CUBIC
        
        public          ::         distribute
        
        
    contains
!---^^^^^^^^

        subroutine distribute( n,a, q )
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
            integer,intent(in)          ::      n       !   make n points
            integer,intent(in)          ::      a       !   candidates at each step m = a n + 1
            type(Quaternion),dimension(:),intent(out)   ::      q       
            
            integer                                     ::      ii
            real(kind=real64),dimension(4,n)            ::     x_fixed 
            real(kind=real64),dimension(4)              ::     x_best
             
            q(1) = quaternion_identity
            x_fixed(:,1) = getComponents(q(1))
            do ii = 1,n-1
                call bestCandidate( x_fixed(:,1:ii), a, quatDistance, quatTrial, x_best )
                x_fixed(:,ii+1) = x_best(:)
                q(ii+1) =  Quaternion_ctor(x_best(:))
            end do
            
            do ii = 1,n-1
                print *,"min dist ",ii,minval( distance( q(ii),q(ii+1:) ) )
            end do
            
            
            return 
        end subroutine distribute
            
            
        function quatDistance(x1,x2) result(d)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      find distance between quaternions expressed as vectors
            real(kind=real64),dimension(:),intent(in)       ::      x1,x2
            real(kind=real64)                               ::      d
            
             
            type(Quaternion)                    ::      q1,q2 
            real(kind=real64),dimension(3,3)    ::      R1,R2
            integer             ::      ii,jj 
            real(kind=real64)   ::      dd
            q1 = Quaternion_ctor(x1)
            q2 = Quaternion_ctor(x2)
            
            R1 = quaternionToRotMat(q1) 
            R2 = quaternionToRotMat(q2) 
            
            d = -huge(1.0)
            do ii = 1,3
                do jj = 1,3
                    
                    dd = R1(1,jj)*R2(1,ii) + R1(2,jj)*R2(2,ii) + R1(3,jj)*R2(3,ii)  
                    d  = max(abs(dd),d)          !   find the maximum dot product between all arms of the jack
                    
                end do
            end do                 
            
            d = (1 - d)/2            
            
            return
        end function quatDistance
        

        
        function quatTrial(d) result(x)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      generate a trial quaternion with symmetry imposed.
            integer,intent(in)                  ::      d
            real(kind=real64),dimension(d)      ::      x
            real(kind=real64),dimension(3,3)    ::      RR
            type(Quaternion)                    ::      qq
            RR = rndRotMat()
            qq = Quaternion_ctor(RR)
            call imposeSymmetry( qq , distributeQuaternions_sym )
            x = getComponents(qq)
            return
        end function quatTrial
        
! 
!         subroutine imposeSymmetry( q )
!     !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
!             type(Quaternion),intent(inout)      ::      q
!             
!             integer             ::      ss,sbest
!             real(kind=real64)   ::      dd,dbest
!             type(Quaternion)    ::      qs
!             type(Quaternion),dimension(:),pointer   ::      qsym
!             
!             
!             select case(distributeQuaternions_sym)
!                 case(SYMMETRY_NONE)
!                     return
!                 case(SYMMETRY_CUBIC)
!                     qsym => cubic_rot_q                    
!             end select
!             
!             
!             
!             dbest = huge(1.0) ; sbest = 0             
!             do ss = 1,size(qsym)
!                 qs = q * qsym(ss) 
!                 dd = distance( qs,quaternion_identity )
!                 if (dd < dbest) then
!                     dbest = dd
!                     sbest = ss
!                 end if
!             end do                
!             q = qsym(sbest) * q
!             
!             
!             
!             return
!         end subroutine imposeSymmetry
! 
!         


    end module distributeQuaternions    
            
    
!   gfortran -ffree-line-length-256 -Og -g ${MYF90LIB}/Lib_Quaternions.f90 ${MYF90LIB}/Lib_RotationMatrices.f90 ${MYF90LIB}/Lib_MitchellsBestCandidate.f90 src/distributeQuaternions.f90 -o testDistributeQuaternions.exe -llapack
    
    program testDistributeQuaternions
!---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
        use distributeQuaternions
        use Lib_RotationMatrices
        use Lib_Quaternions
        use iso_fortran_env
        implicit none
        
        include "cubic_rot_q_data.h"
        
        
        integer                                         ::      N,a
        type(Quaternion),dimension(:),allocatable       ::      q
        
        integer                                 ::      ii
        real(kind=real64),dimension(3,3)        ::      RR
        N = 40
        a = 10
        
        print *,"enter N,a"
        read(*,*) N,a
        
        allocate(q(N))
        
      !q(1) = Quaternion_ctor( rndRotMat() )
      !do ii = 1,24
      !    print *,ii
      !    call report( quaternionToRotMat( q(1) * cubic_rot_q(ii) ) )
      !end do
      !stop
            
        
        
        
        
        
        
        
        
        
        
        
        call distribute( n,a, q )
        
        open(unit=600,file="testDistribute.xyz",action="write")
            write(unit=600,fmt=*) N*6
            !write(unit=600,fmt='(a)') "Lattice=""2.0 0.0 0.0 0.0 2.0 0.0 0.0 0.0 2.0"" Properties=species:S:1:pos:R:3:conc:R:1"
            write(unit=600,fmt='(a)') "atom position_x position_y position_z charge"
            
            do ii = 1,N
                RR = quaternionToRotMat(q(ii))
                write(unit=600,fmt='(a3,100f12.5)')  "J ", RR(:,1),ii*1.0
                write(unit=600,fmt='(a3,100f12.5)')  "J ", RR(:,2),ii*1.0
                write(unit=600,fmt='(a3,100f12.5)')  "J ", RR(:,3),ii*1.0
                write(unit=600,fmt='(a3,100f12.5)')  "J ",-RR(:,1),ii*1.0
                write(unit=600,fmt='(a3,100f12.5)')  "J ",-RR(:,2),ii*1.0
                write(unit=600,fmt='(a3,100f12.5)')  "J ",-RR(:,3),ii*1.0                
            end do
        close(unit=600)
        
        return
    end program testDistributeQuaternions
        