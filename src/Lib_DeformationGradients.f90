
    module Lib_DeformationGradients
!---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
!*      A few simple routines to work with deformation gradients 
!*      constructed as T = (I + e)R 
        use iso_fortran_env
        use Lib_Quaternions
        use Lib_RotationMatrices
        implicit none
        private
        
        external        ::      DSYEV
        
    !---

        public          ::      StrainAndRotMatToDefGrad
        public          ::      StrainAndQuaternionToDefGrad
        public          ::      DefGradToStrainAndRotMat
        public          ::      DefGradToRotMat
        public          ::      DefGradToStrain
        public          ::      DefGradToQuaternion
        public          ::      distance 
        public          ::      imposeSymmetry
        public          ::      addDefGrad,addDefGrad_array
        public          ::      scaleDefGrad,scaleDefGrad_array
        
        public          ::      computeAvgDefGrad
        public          ::      strainInvariants
        
    !---
    
        logical,public          ::      DeformationGradient_dbg = .false.
        
        integer,public,parameter        ::      DEFGRAD_SYMMETRY_NONE = 0
        integer,public,parameter        ::      DEFGRAD_SYMMETRY_CUBIC = 1
        real(kind=real64),dimension(3,3),public,parameter        ::      DEFGRAD_IDENTITY = reshape( (/1.d0,0.d0,0.d0,0.d0,1.d0,0.d0,0.d0,0.d0,1.d0/),(/3,3/) )
                              
    !---
    
        interface   distance 
            module procedure    distanceBetweenDefGrad
        end interface
        
        interface imposeSymmetry
            module procedure    imposeSymmetry0
        end interface
        
        interface computeAvgDefGrad
            module procedure    computeAvgDefGradx
            module procedure    computeAvgDefGradxa
            module procedure    computeAvgDefGrad0
            module procedure    computeAvgDefGrad0a
            module procedure    computeAvgDefGrad0b
            module procedure    computeAvgDefGrad0c
            module procedure    computeAvgDefGrad0d
            module procedure    computeAvgDefGrad0e
            module procedure    computeAvgDefGrad1
            module procedure    computeAvgDefGrad1a
        end interface
        
         
    contains
!---^^^^^^^^

        subroutine strainInvariants(eps,I1,I2,I3)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      compute the three strain invariants
            real(kind=real64),dimension(3,3),intent(in)     ::      eps
            real(kind=real64),intent(out)                   ::      I1,I2,I3
            I1 = eps(1,1)+eps(2,2)+eps(3,3)
            I2 = eps(1,1)*eps(2,2)+eps(2,2)*eps(3,3)+eps(3,3)*eps(1,1)-eps(1,2)*eps(2,1)-eps(2,3)*eps(3,2)-eps(3,1)*eps(1,3) 
            I3 = eps(1,1)*eps(2,2)*eps(3,3)+eps(1,2)*eps(2,3)*eps(3,1)+eps(3,2)*eps(2,1)*eps(1,3)-eps(3,2)*eps(1,1)*eps(2,3)-eps(1,3)*eps(2,2)*eps(3,1)-eps(1,2)*eps(3,3)*eps(2,1)
            return
        end subroutine strainInvariants
    

        subroutine StrainAndRotMatToDefGrad( eps,R , T )
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      T = (I + e)R
            real(kind=real64),dimension(3,3),intent(in)     ::      eps,R
            real(kind=real64),dimension(3,3),intent(out)    ::      T
            
            T(1,1) = (        1 + eps(1,1) )*R(1,1)     + ( eps(2,1) + eps(1,2) )*R(2,1)/2  + ( eps(3,1) + eps(1,3) )*R(3,1)/2
            T(2,1) = ( eps(2,1) + eps(1,2) )*R(1,1)/2   + (        1 + eps(2,2) )*R(2,1)    + ( eps(2,3) + eps(3,2) )*R(3,1)/2  
            T(3,1) = ( eps(3,1) + eps(1,3) )*R(1,1)/2   + ( eps(3,2) + eps(2,3) )*R(2,1)/2  + (        1 + eps(3,3) )*R(3,1)   
            T(1,2) = (        1 + eps(1,1) )*R(1,2)     + ( eps(2,1) + eps(1,2) )*R(2,2)/2  + ( eps(3,1) + eps(1,3) )*R(3,2)/2  
            T(2,2) = ( eps(2,1) + eps(1,2) )*R(1,2)/2   + (        1 + eps(2,2) )*R(2,2)    + ( eps(2,3) + eps(3,2) )*R(3,2)/2  
            T(3,2) = ( eps(3,1) + eps(1,3) )*R(1,2)/2   + ( eps(3,2) + eps(2,3) )*R(2,2)/2  + (        1 + eps(3,3) )*R(3,2)  
            T(1,3) = (        1 + eps(1,1) )*R(1,3)     + ( eps(2,1) + eps(1,2) )*R(2,3)/2  + ( eps(3,1) + eps(1,3) )*R(3,3)/2  
            T(2,3) = ( eps(2,1) + eps(1,2) )*R(1,3)/2   + (        1 + eps(2,2) )*R(2,3)    + ( eps(2,3) + eps(3,2) )*R(3,3)/2  
            T(3,3) = ( eps(3,1) + eps(1,3) )*R(1,3)/2   + ( eps(3,2) + eps(2,3) )*R(2,3)/2  + (        1 + eps(3,3) )*R(3,3)  

            return
        end subroutine StrainAndRotMatToDefGrad
        

 
        subroutine StrainAndQuaternionToDefGrad( eps,q , T )
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      T = (I + e)R
            real(kind=real64),dimension(3,3),intent(in)     ::      eps
            type(Quaternion),intent(in)                     ::      q
            real(kind=real64),dimension(3,3),intent(out)    ::      T
            real(kind=real64),dimension(3,3)        ::      RR
            RR = quaternionToRotMat(q)
            call StrainAndRotMatToDefGrad( eps,RR , T )
            return
        end subroutine StrainAndQuaternionToDefGrad
        

        subroutine DefGradToStrainAndRotMat( T , eps,R )
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      T = (I + e)R
    !*      so T T* = (1+e)R R* (1+e)* = (1+e)(1+e)* = (1+e)^2
            real(kind=real64),dimension(3,3),intent(in)     ::      T
            real(kind=real64),dimension(3,3),intent(out)    ::      eps,R
             
            real(kind=real64),dimension(3,3)        ::      XX,iOnePlusEps
            real(kind=real64),dimension(3)          ::      BB  
            real(kind=real64),dimension(3*66)       ::      work
            integer                                 ::      ii
            
        !---    compute XX = T T* and find its eigendecomposition            
            XX(1:3,1) = T(1:3,1)*T(1,1) + T(1:3,2)*T(1,2) +  T(1:3,3)*T(1,3)
            XX(1:3,2) = T(1:3,1)*T(2,1) + T(1:3,2)*T(2,2) +  T(1:3,3)*T(2,3)
            XX(1:3,3) = T(1:3,1)*T(3,1) + T(1:3,2)*T(3,2) +  T(1:3,3)*T(3,3)
            
         !    print *,"square  matrix "
         !    print *,XX(1,1:3)
         !    print *,XX(2,1:3)
         !    print *,XX(3,1:3)
            
            
            call DSYEV("V","U",3,XX,3,BB,work,size(work),ii)
            
            
            
         !   print *,"DefGradToStrainAndRotMat DSYEV ",ii
         !   print *,"lambda1 ",BB(1)," ",XX(1:3,1)
         !   print *,"lambda2 ",BB(2)," ",XX(1:3,2)
         !   print *,"lambda3 ",BB(3)," ",XX(1:3,3)
            
        !---    now find the square root of T T*, and from this the strain matrix 
            BB(1) = sqrt(BB(1))
            BB(2) = sqrt(BB(2))
            BB(3) = sqrt(BB(3))        
                           
            eps(1:3,1) = BB(1)*XX(1:3,1)*XX(1,1) + BB(2)*XX(1:3,2)*XX(1,2) + BB(3)*XX(1:3,3)*XX(1,3) 
            eps(1:3,2) = BB(1)*XX(1:3,1)*XX(2,1) + BB(2)*XX(1:3,2)*XX(2,2) + BB(3)*XX(1:3,3)*XX(2,3) 
            eps(1:3,3) = BB(1)*XX(1:3,1)*XX(3,1) + BB(2)*XX(1:3,2)*XX(3,2) + BB(3)*XX(1:3,3)*XX(3,3) 
            
        !   print *,"square root matrix "
        !   print *,eps(1,1:3)
        !   print *,eps(2,1:3)
        !   print *,eps(3,1:3)
        !   
        !   print *,"square-of-square-root"
        !   R = matmul( eps,transpose(eps) )
        !   print *,R(1,1:3)
        !   print *,R(2,1:3)
        !   print *,R(3,1:3)
            
            
            eps(1,1) = eps(1,1) - 1.0d0
            eps(2,2) = eps(2,2) - 1.0d0
            eps(3,3) = eps(3,3) - 1.0d0
            
        !---    now find the inverse square root of T T*             
             BB(1) = 1/BB(1)
             BB(2) = 1/BB(2)
             BB(3) = 1/BB(3)        
                                        
             iOnePlusEps(1:3,1) = BB(1)*XX(1:3,1)*XX(1,1) + BB(2)*XX(1:3,2)*XX(1,2) + BB(3)*XX(1:3,3)*XX(1,3)
             iOnePlusEps(1:3,2) = BB(1)*XX(1:3,1)*XX(2,1) + BB(2)*XX(1:3,2)*XX(2,2) + BB(3)*XX(1:3,3)*XX(2,3)
             iOnePlusEps(1:3,3) = BB(1)*XX(1:3,1)*XX(3,1) + BB(2)*XX(1:3,2)*XX(3,2) + BB(3)*XX(1:3,3)*XX(3,3)
            
        !---    using this we can find the rotation matrix 
            R(1:3,1) = iOnePlusEps(1:3,1)*T(1,1) + iOnePlusEps(1:3,2)*T(2,1) + iOnePlusEps(1:3,3)*T(3,1)    
            R(1:3,2) = iOnePlusEps(1:3,1)*T(1,2) + iOnePlusEps(1:3,2)*T(2,2) + iOnePlusEps(1:3,3)*T(3,2)   
            R(1:3,3) = iOnePlusEps(1:3,1)*T(1,3) + iOnePlusEps(1:3,2)*T(2,3) + iOnePlusEps(1:3,3)*T(3,3)   
        
           ! print *,"DefGradToStrainAndRotMat ",dot_product(R(1:3,1),R(1:3,2)) ,dot_product(R(1:3,2),R(1:3,3)),dot_product(R(1:3,3),R(1:3,1))
        !   
        !  print *,"DefGradToStrainAndRotMat |T| = ",determinant3Mat(T) 
        !   print *,"DefGradToStrainAndRotMat |eps| = ",determinant3Mat(eps) 
        !   print *,"DefGradToStrainAndRotMat |R| = ",determinant3Mat(R) 
        !   
            return
        end subroutine DefGradToStrainAndRotMat
        
        

!         pure function determinant3Mat(M) result(d)
!     !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
!     !*      returns the determinant of M
!             real(kind=real64),dimension(3,3),intent(in)      ::      M
!             real(kind=real64)                                ::      d
!             real(kind=real64),dimension(9)       ::      dd
!             dd(1) = M(1,1)*( M(2,2)*M(3,3) - M(2,3)*M(3,2) )
!             dd(2) = M(1,2)*( M(2,3)*M(3,1) - M(2,1)*M(3,3) )
!             dd(3) = M(1,3)*( M(2,1)*M(3,2) - M(2,2)*M(3,1) )
!             dd(4) = M(2,1)*( M(3,2)*M(1,3) - M(3,3)*M(1,2) )
!             dd(5) = M(2,2)*( M(3,3)*M(1,1) - M(3,1)*M(1,3) )
!             dd(6) = M(2,3)*( M(3,1)*M(1,2) - M(3,2)*M(1,1) )
!             dd(7) = M(3,1)*( M(1,2)*M(2,3) - M(1,3)*M(2,2) )
!             dd(8) = M(3,2)*( M(1,3)*M(2,1) - M(1,1)*M(2,3) )
!             dd(9) = M(3,3)*( M(1,1)*M(2,2) - M(1,2)*M(2,1) )
!             d = (1.0d0/3.0d0) * sum(dd)
!             return
!         end function determinant3Mat
!         
        
        function addDefGrad( T1,T2,x1,x2 ) result(T)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      find T = T1 + T2
    !*      optionally T = x1 T1 + x2 T2
            real(kind=real64),dimension(3,3),intent(in)     ::      T1,T2
            real(kind=real64),intent(in),optional           ::      x1,x2
            real(kind=real64),dimension(3,3)                ::      T 
            
            real(kind=real64),dimension(3,3)                ::      eps1,eps2,RR
            type(Quaternion)                                ::      q1,q2,qx 
            real(kind=real64)                               ::      xx
            
            if (present(x1)) then
                if (abs(x1)<1.0d-16) then
                    if (abs(x2)<1.0d-16) then
                        T(1:3,1) = (/ 1,0,0 /)
                        T(1:3,2) = (/ 0,1,0 /)
                        T(1:3,3) = (/ 0,0,1 /)      !   x1 = x2 = 0
                    else
                        call DefGradToStrainAndRotMat( T2 , eps2,RR )
                        q2 = Quaternion_ctor(RR)
                        call StrainAndQuaternionToDefGrad( x2*eps2,q2 , T )                   !   x1 = 0, x2 /= 0
                    end if
                else
                    if (abs(x2)<1.0d-16) then
                        call DefGradToStrainAndRotMat( T1 , eps1,RR )
                        q1 = Quaternion_ctor(RR)
                        call StrainAndQuaternionToDefGrad( x1*eps1,q1 , T )                   !   x1 /= 0, x2 = 0
                    else
                        call DefGradToStrainAndRotMat( T1 , eps1,RR )
                        q1 = Quaternion_ctor(RR)                
                        call DefGradToStrainAndRotMat( T2 , eps2,RR )
                        q2 = Quaternion_ctor(RR)
                        xx = x2/(x1+x2)
                        qx = slerp( q1,q2, xx )
                        eps1 = eps1*x1 + eps2*x2
                        call StrainAndQuaternionToDefGrad( eps1,qx , T )
                    end if
                end if               
            else
                call DefGradToStrainAndRotMat( T1 , eps1,RR )
                q1 = Quaternion_ctor(RR)                
                call DefGradToStrainAndRotMat( T2 , eps2,RR )
                q2 = Quaternion_ctor(RR)                
                qx = midpoint( q1,q2 )
                eps1 = (eps1 + eps2)/2
                call StrainAndQuaternionToDefGrad( eps1,qx , T )
            end if
            
            return
        end function addDefGrad
            
            
        function addDefGrad_array( T1,T2,x1,x2 ) result(T)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      find T = T1 + T2
    !*      optionally T = x1 T1 + x2 T2
            real(kind=real64),dimension(:),intent(in)       ::      T1,T2
            real(kind=real64),intent(in)                    ::      x1,x2
            real(kind=real64),dimension(size(T1))           ::      T 
            
            real(kind=real64),dimension(3,3)                ::      eps1,eps2,RR
            type(Quaternion)                                ::      q1,q2,qx 
            real(kind=real64)                               ::      xx
            
            if (abs(x1)<1.0d-16) then
                if (abs(x2)<1.0d-16) then
                    T(1:9) = (/ 1,0,0 , 0,1,0 , 0,0,1 /)      !   x1 = x2 = 0
                else
                    call DefGradToStrainAndRotMat( reshape( T2,(/3,3/)) , eps2,RR )
                    q2 = Quaternion_ctor(RR)
                    call StrainAndQuaternionToDefGrad( x2*eps2,q2 , RR )                   !   x1 = 0, x2 /= 0
                    T(1:3) = RR(1:3,1)
                    T(4:6) = RR(1:3,2)
                    T(7:9) = RR(1:3,3)
                end if
            else
                if (abs(x2)<1.0d-16) then
                    call DefGradToStrainAndRotMat( reshape( T1,(/3,3/)) , eps1,RR )
                    q1 = Quaternion_ctor(RR)
                    call StrainAndQuaternionToDefGrad( x1*eps1,q1 , RR )                   !   x1 /= 0, x2 = 0
                    T(1:3) = RR(1:3,1)
                    T(4:6) = RR(1:3,2)
                    T(7:9) = RR(1:3,3)
                else
                    call DefGradToStrainAndRotMat( reshape( T1,(/3,3/)) , eps1,RR )
                    q1 = Quaternion_ctor(RR)                
                    call DefGradToStrainAndRotMat( reshape( T2,(/3,3/)) , eps2,RR )
                    q2 = Quaternion_ctor(RR)
                    xx = x2/(x1+x2)
                    qx = slerp( q1,q2, xx )
                    eps1 = eps1*x1 + eps2*x2
                    call StrainAndQuaternionToDefGrad( eps1,qx , RR )
                    T(1:3) = RR(1:3,1)
                    T(4:6) = RR(1:3,2)
                    T(7:9) = RR(1:3,3)
                end if
            end if               
            
            return
        end function addDefGrad_array
            
        function scaleDefGrad( T1,x ) result(T)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ 
    !*      find T = x T1 
  
            real(kind=real64),dimension(3,3),intent(in)     ::      T1
            real(kind=real64),intent(in)                    ::      x
            real(kind=real64),dimension(3,3)                ::      T 
            
            real(kind=real64),dimension(3,3)                ::      eps1,RR 
    
            
            if (abs(x)<1.0d-16) then
                T(1:3,1:3) = reshape( (/ 1,0,0 , 0,1,0 , 0,0,1 /),(/3,3/) )
            else
                call DefGradToStrainAndRotMat( T1 , eps1,RR )
                call StrainAndRotMatToDefGrad( x*eps1,RR , T )
               
            end if               
            
            return
        end function scaleDefGrad
            
        function scaleDefGrad_array( T1,x ) result(T)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      find T = x T1 
  
            real(kind=real64),dimension(:),intent(in)       ::      T1
            real(kind=real64),intent(in)                    ::      x
            real(kind=real64),dimension(size(T1))           ::      T 
            
            real(kind=real64),dimension(3,3)                ::      eps1,RR,TT
    
            
            if (abs(x)<1.0d-16) then
                T(1:9) = (/ 1,0,0 , 0,1,0 , 0,0,1 /)       
            else
                call DefGradToStrainAndRotMat( reshape( T1,(/3,3/)) , eps1,RR )
                call StrainAndRotMatToDefGrad( x*eps1,RR , TT )
                T(1:3) = TT(1:3,1)
                T(4:6) = TT(1:3,2)
                T(7:9) = TT(1:3,3)
            end if               
            
            return
        end function scaleDefGrad_array
            
        
        subroutine DefGradToRotMat( T , R )
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      T = (I + e)R
    !*      so T T* = (1+e)R R* (1+e)* = (1+e)(1+e)* = (1+e)^2
            real(kind=real64),dimension(3,3),intent(in)     ::      T
            real(kind=real64),dimension(3,3),intent(out)    ::      R
             
            real(kind=real64),dimension(3,3)        ::      XX,iOnePlusEps
            real(kind=real64),dimension(3)          ::      BB  
            real(kind=real64),dimension(3*66)       ::      work
            integer                                 ::      ii
            
        !---    compute XX = T T* and find its eigendecomposition            
            XX(1:3,1) = T(1,1:3)*T(1,1) + T(2,1:3)*T(2,1) + T(3,1:3)*T(3,1) 
            XX(1:3,2) = T(1,1:3)*T(1,2) + T(2,1:3)*T(2,2) + T(3,1:3)*T(3,2) 
            XX(1:3,3) = T(1,1:3)*T(1,3) + T(2,1:3)*T(2,3) + T(3,1:3)*T(3,3)             
            call DSYEV("V","U",3,XX,3,BB,work,size(work),ii)
            
            
            
        !---    now find the inverse square root of T T*             
             BB(1) = 1/sqrt(BB(1))    
             BB(2) = 1/sqrt(BB(2))    
             BB(3) = 1/sqrt(BB(3))         
                                        
             iOnePlusEps(1:3,1) = BB(1)*XX(1:3,1)*XX(1,1) + BB(2)*XX(1:3,2)*XX(1,2) + BB(3)*XX(1:3,3)*XX(1,3)
             iOnePlusEps(1:3,2) = BB(1)*XX(1:3,1)*XX(2,1) + BB(2)*XX(1:3,2)*XX(2,2) + BB(3)*XX(1:3,3)*XX(2,3)
             iOnePlusEps(1:3,3) = BB(1)*XX(1:3,1)*XX(3,1) + BB(2)*XX(1:3,2)*XX(3,2) + BB(3)*XX(1:3,3)*XX(3,3)
            
        !---    using this we can find the rotation matrix 
            R(1:3,1) = iOnePlusEps(1:3,1)*T(1,1) + iOnePlusEps(1:3,2)*T(2,1) + iOnePlusEps(1:3,3)*T(3,1)    
            R(1:3,2) = iOnePlusEps(1:3,1)*T(1,2) + iOnePlusEps(1:3,2)*T(2,2) + iOnePlusEps(1:3,3)*T(3,2)   
            R(1:3,3) = iOnePlusEps(1:3,1)*T(1,3) + iOnePlusEps(1:3,2)*T(2,3) + iOnePlusEps(1:3,3)*T(3,3)   
        
            
            return
        end subroutine DefGradToRotMat
        
        
        subroutine DefGradToQuaternion( T , q )
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      T = (I + e)R
    !*      so T T* = (1+e)R R* (1+e)* = (1+e)(1+e)* = (1+e)^2
            real(kind=real64),dimension(3,3),intent(in)     ::      T
            type(Quaternion),intent(out)                    ::      q
            real(kind=real64),dimension(3,3)        ::      RR
            call DefGradToRotMat( T , RR )  
            q = Quaternion_ctor(RR)
            return
        end subroutine DefGradToQuaternion
        
        
        subroutine DefGradToStrain( T , eps )
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      T = (I + e)R
    !*      so T T* = (1+e)R R* (1+e)* = (1+e)(1+e)* = (1+e)^2
            real(kind=real64),dimension(3,3),intent(in)     ::      T
            real(kind=real64),dimension(3,3),intent(out)    ::      eps
             
            real(kind=real64),dimension(3,3)        ::      XX 
            real(kind=real64),dimension(3)          ::      BB  
            real(kind=real64),dimension(3*66)       ::      work
            integer                                 ::      ii
            
        !---    compute XX = T T* and find its eigendecomposition            
            XX(1:3,1) = T(1:3,1)*T(1,1) + T(1:3,2)*T(1,2) +  T(1:3,3)*T(1,3)
            XX(1:3,2) = T(1:3,1)*T(2,1) + T(1:3,2)*T(2,2) +  T(1:3,3)*T(2,3)
            XX(1:3,3) = T(1:3,1)*T(3,1) + T(1:3,2)*T(3,2) +  T(1:3,3)*T(3,3)
             
            call DSYEV("V","U",3,XX,3,BB,work,size(work),ii)
            
             
            
        !---    now find the square root of T T*, and from this the strain matrix 
            BB(1) = sqrt(BB(1))
            BB(2) = sqrt(BB(2))
            BB(3) = sqrt(BB(3))        
                           
            eps(1:3,1) = BB(1)*XX(1:3,1)*XX(1,1) + BB(2)*XX(1:3,2)*XX(1,2) + BB(3)*XX(1:3,3)*XX(1,3) 
            eps(1:3,2) = BB(1)*XX(1:3,1)*XX(2,1) + BB(2)*XX(1:3,2)*XX(2,2) + BB(3)*XX(1:3,3)*XX(2,3) 
            eps(1:3,3) = BB(1)*XX(1:3,1)*XX(3,1) + BB(2)*XX(1:3,2)*XX(3,2) + BB(3)*XX(1:3,3)*XX(3,3) 
            
            
            eps(1,1) = eps(1,1) - 1.0d0
            eps(2,2) = eps(2,2) - 1.0d0
            eps(3,3) = eps(3,3) - 1.0d0
            
        !   
            return
        end subroutine DefGradToStrain 
        
    !---
        
        
        
        function distanceBetweenDefGrad( T1,T2,symmetry ) result(d)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      find a distance between two deformation gradients, based on the distance between the rotation matrices
            real(kind=real64),dimension(3,3),intent(in)     ::      T1,T2
            integer,intent(in)                              ::      symmetry
            real(kind=real64)                               ::      d
        
            integer                                 ::      ss
            type(Quaternion)                        ::      q1,q2,qs
    
            
            call DefGradToQuaternion( T1 , q1 )            
            call DefGradToQuaternion( T2 , q2 )
            
            d = huge(1.0d0)
            do ss = 1,getNSymmetries( symmetry )
                qs = symmetryRelatedQuaternion(q1,symmetry,ss)    
                d = min( d, distance(qs,q2) )
            end do
            
            return
        end function distanceBetweenDefGrad
            
            
         
            
        subroutine imposeSymmetry0( T,symmetry ) 
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      
            real(kind=real64),dimension(3,3),intent(inout)  ::      T
            integer,intent(in)                              ::      symmetry

            real(kind=real64),dimension(3,3)        ::      RR,eps
            type(Quaternion)                        ::      qq!,q2
!           integer             ::      ii,ibest
!            real(kind=real64)   ::      dd,dbest
          !  
!             real(kind=real64),dimension(3,3,48),parameter   ::      CUBSYMMAT = reshape( (/                                 &
!                     1, 0, 0 ,  0, 1, 0 ,  0, 0, 1,    0, 1, 0 ,  1, 0, 0 ,  0, 0, 1,     0, 1, 0 ,  0, 0, 1 , 1, 0, 0,      &
!                     1, 0, 0 ,  0, 0, 1 ,  0, 1, 0,    0, 0, 1 ,  1, 0, 0 ,  0, 1, 0,     0, 0, 1 ,  0, 1, 0 , 1, 0, 0,      &
!                     1, 0, 0 ,  0,-1, 0 ,  0, 0, 1,    0,-1, 0 ,  1, 0, 0 ,  0, 0, 1,     0,-1, 0 ,  0, 0, 1 , 1, 0, 0,      &
!                     1, 0, 0 ,  0, 0,-1 ,  0, 0, 1,    0, 0,-1 ,  1, 0, 0 ,  0, 0, 1,     0, 0,-1 ,  0, 0, 1 , 1, 0, 0,      &
!                     1, 0, 0 ,  0, 1, 0 ,  0, 0,-1,    0, 1, 0 ,  1, 0, 0 ,  0, 0,-1,     0, 1, 0 ,  0, 0,-1 , 1, 0, 0,      &
!                     1, 0, 0 ,  0, 0, 1 ,  0,-1, 0,    0, 0, 1 ,  1, 0, 0 ,  0,-1, 0,     0, 0, 1 ,  0,-1, 0 , 1, 0, 0,      &
!                     1, 0, 0 ,  0,-1, 0 ,  0, 0,-1,    0,-1, 0 ,  1, 0, 0 ,  0, 0,-1,     0,-1, 0 ,  0, 0,-1 , 1, 0, 0,      &
!                     1, 0, 0 ,  0, 0,-1 ,  0, 0,-1,    0, 0,-1 ,  1, 0, 0 ,  0, 0,-1,     0, 0,-1 ,  0, 0,-1 , 1, 0, 0,      &
!                    -1, 0, 0 ,  0, 1, 0 ,  0, 0, 1,    0, 1, 0 , -1, 0, 0 ,  0, 0, 1,     0, 1, 0 ,  0, 0, 1 ,-1, 0, 0,      &
!                    -1, 0, 0 ,  0, 0, 1 ,  0, 1, 0,    0, 0, 1 , -1, 0, 0 ,  0, 1, 0,     0, 0, 1 ,  0, 1, 0 ,-1, 0, 0,      &
!                    -1, 0, 0 ,  0,-1, 0 ,  0, 0, 1,    0,-1, 0 , -1, 0, 0 ,  0, 0, 1,     0,-1, 0 ,  0, 0, 1 ,-1, 0, 0,      &
!                    -1, 0, 0 ,  0, 0,-1 ,  0, 0, 1,    0, 0,-1 , -1, 0, 0 ,  0, 0, 1,     0, 0,-1 ,  0, 0, 1 ,-1, 0, 0,      &
!                    -1, 0, 0 ,  0, 1, 0 ,  0, 0,-1,    0, 1, 0 , -1, 0, 0 ,  0, 0,-1,     0, 1, 0 ,  0, 0,-1 ,-1, 0, 0,      &
!                    -1, 0, 0 ,  0, 0, 1 ,  0,-1, 0,    0, 0, 1 , -1, 0, 0 ,  0,-1, 0,     0, 0, 1 ,  0,-1, 0 ,-1, 0, 0,      &
!                    -1, 0, 0 ,  0,-1, 0 ,  0, 0,-1,    0,-1, 0 , -1, 0, 0 ,  0, 0,-1,     0,-1, 0 ,  0, 0,-1 ,-1, 0, 0,      &
!                    -1, 0, 0 ,  0, 0,-1 ,  0, 0,-1,    0, 0,-1 , -1, 0, 0 ,  0, 0,-1,     0, 0,-1 ,  0, 0,-1 ,-1, 0, 0   /),(/3,3,48/) )                                 
!             
!           !  
!             ! if (symmetry == QUATERNION_SYMMETRY_CUBIC) then
!                  do ii = 1,size(CUBSYMMAT,dim=3)
!                      qq = Quaternion_ctor( CUBSYMMAT(:,:,ii) )
!                      call report(qq,asRotMat=.true.)
!                      
!                  end do
!                 stop                    
!                  return
!              !end if
            
            
            
            
            
            
            
            
            !return
            call DefGradToStrainAndRotMat( T , eps, RR )           
             
              
             
            qq = Quaternion_ctor(RR)
            call imposeSymmetry( qq,symmetry )   
            
            !!   
            !dd = distance(qq)
            !RR(1:3,1) = -RR(1:3,1)      !   mirror flip
            !q2 = Quaternion_ctor(RR)
            !call imposeSymmetry( q2,symmetry )   
            !if (distance(q2)<dd) then
            !    eps(1:3,1) = -eps(1:3,1)
            !    call StrainAndQuaternionToDefGrad( eps, q2,T )      
            !    return
            !end if
            !
           
            call StrainAndQuaternionToDefGrad( eps, qq,T )  
            
            !print *,"imposesym ",determinant3Mat(RR),determinant3Mat(T)
                        
            return
        end subroutine imposeSymmetry0
            
            
        subroutine computeAvgDefGradx(eps_sum,dg_sum,weight_sum)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ 
    !*      initialise constructing an average deformation gradient
            real(kind=real64),dimension(3,3),intent(inout)  ::      eps_sum,dg_sum
            real(kind=real64),intent(inout)                 ::      weight_sum
            
            eps_sum = 0.0d0
            dg_sum = 0.0d0
            weight_sum = 0.0d0
            
            return
        end subroutine computeAvgDefGradx
            
        subroutine computeAvgDefGradxa(eps_sum,dg_sum,weight_sum)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      initialise constructing an average deformation gradient
            real(kind=real64),dimension(3,3),intent(inout)  ::      dg_sum
            real(kind=real64),dimension(6),intent(inout)    ::      eps_sum
            real(kind=real64),intent(inout)                 ::      weight_sum
            
            eps_sum = 0.0d0
            dg_sum = 0.0d0
            weight_sum = 0.0d0
            
            return
        end subroutine computeAvgDefGradxa
            
        subroutine computeAvgDefGrad0(dg_in,weight,eps_sum,dg_sum,weight_sum)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      start/continue constructing an average deformation gradient
            real(kind=real64),dimension(3,3),intent(in)     ::      dg_in
            real(kind=real64),intent(in)                    ::      weight
            
            real(kind=real64),dimension(3,3),intent(inout)  ::      eps_sum,dg_sum
            real(kind=real64),intent(inout)                 ::      weight_sum
            
            real(kind=real64),dimension(3,3)            ::      eps
            
            call DefGradToStrain( dg_in , eps )
            eps_sum = eps_sum + weight*eps
            dg_sum = dg_sum + weight*dg_in
            weight_sum = weight_sum + weight
            
            return
        end subroutine computeAvgDefGrad0
            
        subroutine computeAvgDefGrad0a(dg_in,weight,eps_sum,dg_sum,weight_sum)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      start/continue constructing an average deformation gradient
            real(kind=real64),dimension(3,3),intent(in)     ::      dg_in
            real(kind=real64),intent(in)                    ::      weight
            
            real(kind=real64),dimension(3,3),intent(inout)  ::      dg_sum
            real(kind=real64),dimension(6),intent(inout)    ::      eps_sum
            real(kind=real64),intent(inout)                 ::      weight_sum
            
            real(kind=real64),dimension(3,3)            ::      eps
            
            call DefGradToStrain( dg_in , eps )
            eps_sum(1) = eps_sum(1) + weight*eps(1,1)
            eps_sum(2) = eps_sum(2) + weight*eps(2,2)
            eps_sum(3) = eps_sum(3) + weight*eps(3,3)
            eps_sum(4) = eps_sum(4) + weight*eps(1,2)
            eps_sum(5) = eps_sum(5) + weight*eps(2,3)
            eps_sum(6) = eps_sum(6) + weight*eps(3,1)
            dg_sum = dg_sum + weight*dg_in
            weight_sum = weight_sum + weight
            
            return
        end subroutine computeAvgDefGrad0a
            
        
        
        
        subroutine computeAvgDefGrad0b(dg_in,weight,eps_sum,dg_sum,weight_sum)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      start/continue constructing an average deformation gradient
            real(kind=real64),dimension(9),intent(in)       ::      dg_in
            real(kind=real64),intent(in)                    ::      weight
            
            real(kind=real64),dimension(3,3),intent(inout)  ::      eps_sum,dg_sum
            real(kind=real64),intent(inout)                 ::      weight_sum
            
            real(kind=real64),dimension(3,3)            ::      dg,eps
            dg = reshape(dg_in,(/3,3/))
            call DefGradToStrain( dg , eps )
            eps_sum = eps_sum + weight*eps
            dg_sum = dg_sum + weight*dg
            weight_sum = weight_sum + weight
            
            return
        end subroutine computeAvgDefGrad0b
            
        subroutine computeAvgDefGrad0c(dg_in,weight,eps_sum,dg_sum,weight_sum)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      start/continue constructing an average deformation gradient
            real(kind=real64),dimension(9),intent(in)       ::      dg_in
            real(kind=real64),intent(in)                    ::      weight
            
            real(kind=real64),dimension(3,3),intent(inout)  ::      dg_sum
            real(kind=real64),dimension(6),intent(inout)    ::      eps_sum
            real(kind=real64),intent(inout)                 ::      weight_sum
            
            real(kind=real64),dimension(3,3)            ::      dg,eps
            dg = reshape(dg_in,(/3,3/))
            call DefGradToStrain( dg , eps )
            eps_sum(1) = eps_sum(1) + weight*eps(1,1)
            eps_sum(2) = eps_sum(2) + weight*eps(2,2)
            eps_sum(3) = eps_sum(3) + weight*eps(3,3)
            eps_sum(4) = eps_sum(4) + weight*eps(1,2)
            eps_sum(5) = eps_sum(5) + weight*eps(2,3)
            eps_sum(6) = eps_sum(6) + weight*eps(3,1)
            dg_sum = dg_sum + weight*dg 
            weight_sum = weight_sum + weight
            
            return
        end subroutine computeAvgDefGrad0c
            
        
        
        
        subroutine computeAvgDefGrad0d(dg_in,weight,eps_sum,dg_sum,weight_sum)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      start/continue constructing an average deformation gradient
            real(kind=real32),dimension(9),intent(in)       ::      dg_in
            real(kind=real64),intent(in)                    ::      weight
            
            real(kind=real64),dimension(3,3),intent(inout)  ::      eps_sum,dg_sum
            real(kind=real64),intent(inout)                 ::      weight_sum
            
            real(kind=real64),dimension(3,3)            ::      dg,eps
            dg = reshape(dg_in,(/3,3/))
            call DefGradToStrain( dg , eps )
            eps_sum = eps_sum + weight*eps
            dg_sum = dg_sum + weight*dg
            weight_sum = weight_sum + weight
            
            return
        end subroutine computeAvgDefGrad0d
        
        subroutine computeAvgDefGrad0e(dg_in,weight,eps_sum,dg_sum,weight_sum)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      start/continue constructing an average deformation gradient
            real(kind=real32),dimension(9),intent(in)       ::      dg_in
            real(kind=real64),intent(in)                    ::      weight
            
            real(kind=real64),dimension(3,3),intent(inout)  ::      dg_sum
            real(kind=real64),dimension(6),intent(inout)    ::      eps_sum
            real(kind=real64),intent(inout)                 ::      weight_sum
            
            real(kind=real64),dimension(3,3)            ::      dg,eps
            dg = reshape(dg_in,(/3,3/))
            call DefGradToStrain( dg , eps )
            eps_sum(1) = eps_sum(1) + weight*eps(1,1)
            eps_sum(2) = eps_sum(2) + weight*eps(2,2)
            eps_sum(3) = eps_sum(3) + weight*eps(3,3)
            eps_sum(4) = eps_sum(4) + weight*eps(1,2)
            eps_sum(5) = eps_sum(5) + weight*eps(2,3)
            eps_sum(6) = eps_sum(6) + weight*eps(3,1)
            dg_sum = dg_sum + weight*dg 
            weight_sum = weight_sum + weight
            
            return
        end subroutine computeAvgDefGrad0e
            
        
        
        
        subroutine computeAvgDefGrad1(eps_sum,dg_sum,weight_sum,dg_bar)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      finish constructing an average deformation gradient
            real(kind=real64),dimension(3,3),intent(in)     ::      eps_sum,dg_sum
            real(kind=real64),intent(in)                    ::      weight_sum
            real(kind=real64),dimension(3,3),intent(out)    ::      dg_bar
            
            real(kind=real64),dimension(3,3)            ::      eps,rot
            real(kind=real64)                           ::      iw
            if (weight_sum == 0) then            
                dg_bar = reshape( (/1,0,0 , 0,1,0 , 0,0,1/),(/3,3/) )
                return
            end if            
            
            iw = 1/weight_sum
            dg_bar(1:3,1:3) = dg_sum(1:3,1:3)*iw            
            call DefGradToRotMat( dg_bar , rot )  
            eps(1:3,1:3) = eps_sum(1:3,1:3)*iw
            
            call StrainAndRotMatToDefGrad( eps,rot , dg_bar )
            return
        end subroutine computeAvgDefGrad1
        
        
        subroutine computeAvgDefGrad1a(eps_sum,dg_sum,weight_sum,dg_bar)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ 
    !*      finish constructing an average deformation gradient
            real(kind=real64),dimension(6),intent(in)       ::      eps_sum
            real(kind=real64),dimension(3,3),intent(in)     ::      dg_sum
            real(kind=real64),intent(in)                    ::      weight_sum
            real(kind=real64),dimension(3,3),intent(out)    ::      dg_bar
            
            real(kind=real64),dimension(3,3)            ::      eps,rot
            real(kind=real64)                           ::      iw
            if (weight_sum == 0) then            
                dg_bar = reshape( (/1,0,0 , 0,1,0 , 0,0,1/),(/3,3/) )
                return
            end if            
            
            iw = 1/weight_sum
            
            dg_bar(1:3,1:3) = dg_sum(1:3,1:3)*iw            
            call DefGradToRotMat( dg_bar , rot )  
            
            eps(1,1) = eps_sum(1)*iw
            eps(2,2) = eps_sum(2)*iw
            eps(3,3) = eps_sum(3)*iw
            eps(1,2) = eps_sum(4)*iw
            eps(2,3) = eps_sum(5)*iw
            eps(3,1) = eps_sum(6)*iw
            eps(2,1) = eps(1,2)
            eps(3,2) = eps(2,3)
            eps(1,3) = eps(3,1)
            
            call StrainAndRotMatToDefGrad( eps,rot , dg_bar )
            return
        end subroutine computeAvgDefGrad1a
        
    end module Lib_DeformationGradients
    
    
!   gfortran -ffree-line-length-256 -Og -g ${MYF90LIB}/Lib_Quaternions.f90 ${MYF90LIB}/Lib_RotationMatrices.f90 ${MYF90LIB}/Lib_DeformationGradients.f90 -o testDeformationGradients.exe

    
!    program testDeformationGradients
!!---^^^^^^^^^^^^^^^^^^^^^^^^
!        use iso_fortran_env
!        use Lib_RotationMatrices
!        use Lib_Quaternions
!        use Lib_DeformationGradients
!        implicit none
!        
!        real(kind=real64),dimension(3,3)    ::      T1,T2,T3
!        
!        real(kind=real64),dimension(3,3)    ::      eps1,eps2,eps3, R1,R2,R3
!        
!        call random_number(eps1) ; eps1 = ( eps1 + transpose(eps1) - 2 )*0.1       
!        R1 = rndRotMat()
!        
!        
!        print *,"eps1 " ; call report(eps1)        
!        print *,"R1   " ; call report(R1  )
!        print *,""
!        print *,"distance R1-R2 ",quaternionDistanceInDegrees( distance( Quaternion_ctor(R1),Quaternion_ctor(R2) ) )
!        
!        call StrainAndRotMatToDefGrad(eps1,R1,T1)
!        
!                          
!        print *,"T1   " ; call report(T1  )
!        print *,""
!        
!        print *,"test extract"
!        call DefGradToStrainAndRotMat( T1 , eps1,R1 )
!        print *,"eps1' " ; call report(eps1)    
!        print *,"R1'   " ; call report(R1  )
!        print *,""
!        print *,"...and reconstruct"
!        call StrainAndRotMatToDefGrad(eps1,R1,T1)
!        print *,"T1'   " ; call report(T1  )
!        print *,""
!        
!        print *,"test scale"
!        T3 = scaleDefGrad( T1,0.1d0 )
!        call DefGradToStrainAndRotMat( T3 , eps3,R3 )
!        print *,"eps3 " ; call report(eps3)    
!        print *,"R3   " ; call report(R3  )
!        print *,""
!        
!        
!        
!        call random_number(eps2) ; eps2 = ( eps2 + transpose(eps2) - 2 )*0.1
!        R2 = rndRotMat()
!        call StrainAndRotMatToDefGrad(eps2,R2,T2)
!        print *,"eps2 " ; call report(eps2)
!        print *,"R2   " ; call report(R2  )
!        print *,""
!        print *,"distance R1-R2 ",quaternionDistanceInDegrees( distance( Quaternion_ctor(R1),Quaternion_ctor(R2) ) )
!        
!        
!        print *,"test add"
!        T3 = addDefGrad( T1,T2,0.5d0,0.5d0 )
!        call DefGradToStrainAndRotMat( T3 , eps3,R3 )
!        print *,"eps3 " ; call report(eps3)    
!        print *,"R3   " ; call report(R3  )
!        print *,""
!        print *,"distance R1-R3 ",quaternionDistanceInDegrees( distance( Quaternion_ctor(R1),Quaternion_ctor(R3) ) )
!        print *,"distance R2-R3 ",quaternionDistanceInDegrees( distance( Quaternion_ctor(R2),Quaternion_ctor(R3) ) )
!        
!        
!        
!        print *,""
!        print *,"done"
!        print *,""
!        
!    end program testDeformationGradients
!    
!        
!    
!        
!        
!    