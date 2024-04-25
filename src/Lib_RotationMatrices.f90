
    module Lib_RotationMatrices
!---^^^^^^^^^^^^^^^^^^^^^^^^^^^
!*
!*      A simple module to generate a rotation matrix
!*      from Euler angles
!*      using the x-convention:
!*          First rotate by phi about the z-axis
!*              D = ( cos(phi)  sin(phi)    0   )
!*                  (-sin(phi)  cos(phi)    0   )
!*                  (   0       0           1   )
!*          Then rotate by theta about the x' axis
!*              C = (   1       0           0   )
!*                  (   0   cos(the)    sin(the))
!*                  (   0   -sin(the)   cos(the))
!*          Finally rotate by chi about the z' axis
!*              B = ( cos(chi)  sin(chi)    0   )
!*                  (-sin(chi)  cos(chi)    0   )
!*                  (   0       0           1   )
!*      The rotation matrix A = B C D 
!*      See   https://mathworld.wolfram.com/EulerAngles.html
!*
!*      also provides a few utility functions
!*
!*      Daniel Mason, UKAEA
!*      April 2022
        use Lib_RandomSeed
        use iso_fortran_env
        implicit none
        private
        
    !---
    
        real(kind=real64),dimension(3,3),public,parameter   ::  RotationMatrix_identity &
                                = reshape( (/1.0d0,0.0d0,0.0d0,0.0d0,1.0d0,0.0d0,0.0d0,0.0d0,1.0d0/),(/3,3/) )
            
    !---   
        
        public      ::      RotationMatrix_ctor             !   create a new 3x3 matrix to represent a rotation matrix
        public      ::      report                          !   simple reporting of the rotation matrix    
        
        public      ::      inverseRotationMatrix           !   = transpose 
        
        public      ::      rotateVector                    !   find v = R u       
        public      ::      completeBasis                   !   given the z- direction, and maybe the x-, complete a basis set
        public      ::      tidyRotationMatrix              !   mirror flips and 90 degree rotations to make sure column 1 is closest to x, column 2 is closest to y, column 3 is z. ( cubic symmetry )
        public      ::      rndUnitVec                      !   generate a random unit vector in the sphere
        public      ::      rndRotMat                       !   generate a random rotation matrix
        
    !---
        
        interface   RotationMatrix_ctor
            module procedure        RotationMatrix_null
            module procedure        RotationMatrix_ctor0
            module procedure        RotationMatrix_ctor1
            module procedure        RotationMatrix_ctor2
            module procedure        RotationMatrix_ctor3
        end interface
        
        interface   completeBasis
            module procedure        completeBasis0
        end interface
        
        interface report
            module procedure    report0
        end interface
        
    !---        
        
        interface rotateVector
            module procedure    rotateVector0
            module procedure    rotateVector1
        end interface
        
        
        interface rndRotMat
            module procedure    rndRotMat0
            module procedure    rndRotMat1
        end interface
        
        
        interface rndUnitVec
            module procedure    rndUnitVec0
            module procedure    rndUnitVec1
        end interface
        
    !---        
        
    contains
!---^^^^^^^^    
        pure function RotationMatrix_null( ) result(R)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
            real(kind=real64),dimension(3,3)    ::      R
            R(1:3,1) = (/ 1.0d0,0.0d0,0.0d0 /)
            R(1:3,2) = (/ 0.0d0,1.0d0,0.0d0 /)
            R(1:3,3) = (/ 0.0d0,0.0d0,1.0d0 /)
            return
        end function RotationMatrix_null
        
        pure function RotationMatrix_ctor0( phi,theta,chi ) result(R)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
            real(kind=real64),intent(in)        ::      phi             !   angle to x-axis
            real(kind=real64),intent(in)        ::      theta           !   angle to z-axis
            real(kind=real64),intent(in)        ::      chi             !   rotation azimuthal angle
            real(kind=real64),dimension(3,3)    ::      R
            
            real(kind=real64)       ::      ct,st,cp,sp,cc,sc
            
            ct = cos(theta)
            st = sin(theta)
            cp = cos(phi)
            sp = sin(phi)
            cc = cos(chi)
            sc = sin(chi)
                        
            R(1,1) = cc*cp - ct*sp*sc
            R(1,2) = cc*sp + ct*cp*sc
            R(1,3) = sc*st
            R(2,1) = -sc*cp - ct*sp*cc 
            R(2,2) = -sc*sp + ct*cp*cc 
            R(2,3) = cc*st 
            R(3,1) = st*sp 
            R(3,2) = -st*cp 
            R(3,3) = ct            
            
            return
        end function RotationMatrix_ctor0
            

        pure function RotationMatrix_ctor1( n,theta ) result(R)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      construct a rotation matrix from axis and angle. 
            real(kind=real64),dimension(3),intent(in)        ::      n  !   axis
            real(kind=real64),intent(in)        ::      theta           !   angle
            real(kind=real64),dimension(3,3)    ::      R
            
            real(kind=real64)       ::      ct,st 
            real(kind=real64)       ::      ux,uy,uz
           
            ct = 1/norm2(n)         !   use as a dummy for now...
            ux = n(1)*ct 
            uy = n(2)*ct
            uz = n(3)*ct
            
            ct = cos(theta)
            st = sin(theta)
            
            
            R(1,1) = ct + ux*ux*(1-ct)
            R(1,2) = ux*uy*(1-ct) - uz*st
            R(1,3) = ux*uz*(1-ct) + uy*st
            R(2,1) = uy*ux*(1-ct) + uz*st 
            R(2,2) = ct + uy*uy*(1-ct) 
            R(2,3) = uy*uz*(1-ct) - ux*st
            R(3,1) = uz*ux*(1-ct) - uy*st
            R(3,2) = uz*uy*(1-ct) + ux*st
            R(3,3) = ct + uz*uz*(1-ct)
            
            return
        end function RotationMatrix_ctor1
                    
        function RotationMatrix_ctor2( u9 ) result(R)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      construct a rotation matrix from nine numbers
            real(kind=real64),dimension(9),intent(in)        ::      u9
            real(kind=real64),dimension(3,3)    ::      R
            R(1:3,1) = u9(1:3)
            R(1:3,2) = u9(4:6)
            R(1:3,3) = u9(7:9)
            call tidyRotationMatrix(R)            
            return
        end function RotationMatrix_ctor2
        
        function RotationMatrix_ctor3( x,y,z ) result(R)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      construct a rotation matrix from nine numbers
            real(kind=real64),dimension(3),intent(in)        ::      x,y,z
            real(kind=real64),dimension(3,3)    ::      R
            R(1:3,1) = x(1:3)
            R(1:3,2) = y(1:3)
            R(1:3,3) = z(1:3)
            print *,"det(R)",determinant3Mat(R)
            call restoreUnitaryness(R)         
            print *,"det(R)",determinant3Mat(R)               
            return
        end function RotationMatrix_ctor3
        
    !---
                    
    
         subroutine tidyRotationMatrix(M,norescale)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
            real(kind=real64),dimension(3,3),intent(inout)    ::      M
            logical,intent(in),optional         ::      norescale
!            real(kind=real64)       ::      dd      !,mx,my,mz
            real(kind=real64),dimension(3,3)    ::      mtmp
            
            mtmp = M
            
            if (abs(mtmp(3,3))>=max(abs(mtmp(3,1)),abs(mtmp(3,2)))) then
                !   largest z component is in the z- direction
                !   *,*,z or *,*,-z                
                if (mtmp(3,3) > 0) then
                    !   *,*,z
                    !   ... and z component is +ve direction.
                    if (abs(mtmp(1,1))>=abs(mtmp(1,2))) then
                        !   largest x component in x-direction
                        if (mtmp(1,1) > 0) then
                            !   x,y,z
                            M(1:3,1) =  (/ mtmp(1,1),mtmp(2,1),mtmp(3,1) /)    
                            M(1:3,2) =  (/ mtmp(1,2),mtmp(2,2),mtmp(3,2) /)
                            M(1:3,3) =  (/ mtmp(1,3),mtmp(2,3),mtmp(3,3) /)
                        else
                            !  -x,-y,z
                            M(1:3,1) = -(/ mtmp(1,1),mtmp(2,1),mtmp(3,1) /)    
                            M(1:3,2) = -(/ mtmp(1,2),mtmp(2,2),mtmp(3,2) /)
                            M(1:3,3) =  (/ mtmp(1,3),mtmp(2,3),mtmp(3,3) /)
                        end if
                    else
                        !   largest x component in y-direction
                        if (mtmp(1,2) > 0) then
                            !   y,-x,z
                            M(1:3,1) =  (/ mtmp(1,2),mtmp(2,2),mtmp(3,2) /)    
                            M(1:3,2) = -(/ mtmp(1,1),mtmp(2,1),mtmp(3,1) /)
                            M(1:3,3) =  (/ mtmp(1,3),mtmp(2,3),mtmp(3,3) /)                            
                        else
                            !  -y, x,z
                            M(1:3,1) = -(/ mtmp(1,2),mtmp(2,2),mtmp(3,2) /)    
                            M(1:3,2) =  (/ mtmp(1,1),mtmp(2,1),mtmp(3,1) /)
                            M(1:3,3) =  (/ mtmp(1,3),mtmp(2,3),mtmp(3,3) /)
                        end if
                    end if
                else
                    !   *,*,-z
                    if (abs(mtmp(1,1))>=abs(mtmp(1,2))) then
                        !   largest x component in x-direction
                        if (mtmp(1,1) > 0) then
                            !   x,-y,-z
                            M(1:3,1) =  (/ mtmp(1,1),mtmp(2,1),mtmp(3,1) /)    
                            M(1:3,2) = -(/ mtmp(1,2),mtmp(2,2),mtmp(3,2) /)
                            M(1:3,3) = -(/ mtmp(1,3),mtmp(2,3),mtmp(3,3) /)
                        else
                            !  -x, y,-z
                            M(1:3,1) = -(/ mtmp(1,1),mtmp(2,1),mtmp(3,1) /)    
                            M(1:3,2) =  (/ mtmp(1,2),mtmp(2,2),mtmp(3,2) /)
                            M(1:3,3) = -(/ mtmp(1,3),mtmp(2,3),mtmp(3,3) /)
                        end if
                    else
                        !   largest x component in y-direction
                        if (mtmp(1,2) > 0) then
                            !   y,x,-z
                            M(1:3,1) =  (/ mtmp(1,2),mtmp(2,2),mtmp(3,2) /)    
                            M(1:3,2) =  (/ mtmp(1,1),mtmp(2,1),mtmp(3,1) /)
                            M(1:3,3) = -(/ mtmp(1,3),mtmp(2,3),mtmp(3,3) /)                            
                        else
                            !  -y,-x,-z
                            M(1:3,1) = -(/ mtmp(1,2),mtmp(2,2),mtmp(3,2) /)    
                            M(1:3,2) = -(/ mtmp(1,1),mtmp(2,1),mtmp(3,1) /)
                            M(1:3,3) = -(/ mtmp(1,3),mtmp(2,3),mtmp(3,3) /)
                        end if
                    end if
                end if
                
            else if (abs(mtmp(3,2))>=abs(mtmp(3,1))) then    
                !   largest z component in y direction
                !   *,*,y or *,*,-y
                 if (mtmp(3,2) > 0) then
                    !   ... and z component is +ve direction.
                    !   *,*,y                     
                    if (abs(mtmp(1,1))>=abs(mtmp(1,3))) then
                        !   largest x component in x-direction
                        if (mtmp(1,1) > 0) then
                            !   x,-z,y
                            M(1:3,1) =  (/ mtmp(1,1),mtmp(2,1),mtmp(3,1) /)    
                            M(1:3,2) = -(/ mtmp(1,3),mtmp(2,3),mtmp(3,3) /)
                            M(1:3,3) =  (/ mtmp(1,2),mtmp(2,2),mtmp(3,2) /)
                        else
                            !  -x,z,y
                            M(1:3,1) = -(/ mtmp(1,1),mtmp(2,1),mtmp(3,1) /)    
                            M(1:3,2) =  (/ mtmp(1,3),mtmp(2,3),mtmp(3,3) /)
                            M(1:3,3) =  (/ mtmp(1,2),mtmp(2,2),mtmp(3,2) /)
                        end if
                    else
                        !   largest x component in z-direction
                        if (mtmp(1,3) > 0) then
                            !   z, x,y
                            M(1:3,1) =  (/ mtmp(1,3),mtmp(2,3),mtmp(3,3) /)    
                            M(1:3,2) =  (/ mtmp(1,1),mtmp(2,1),mtmp(3,1) /)
                            M(1:3,3) =  (/ mtmp(1,2),mtmp(2,2),mtmp(3,2) /)                            
                        else
                            !  -z,-x,z
                            M(1:3,1) = -(/ mtmp(1,3),mtmp(2,3),mtmp(3,3) /)    
                            M(1:3,2) = -(/ mtmp(1,1),mtmp(2,1),mtmp(3,1) /)
                            M(1:3,3) =  (/ mtmp(1,2),mtmp(2,2),mtmp(3,2) /)
                        end if
                    end if
                else
                    !   *,*,-y                                     
                    if (abs(mtmp(1,1))>=abs(mtmp(1,3))) then
                        !   largest x component in x-direction
                        if (mtmp(1,1) > 0) then
                            !   x,z,-y
                            M(1:3,1) =  (/ mtmp(1,1),mtmp(2,1),mtmp(3,1) /)    
                            M(1:3,2) =  (/ mtmp(1,3),mtmp(2,3),mtmp(3,3) /)
                            M(1:3,3) = -(/ mtmp(1,2),mtmp(2,2),mtmp(3,2) /)
                        else
                            !  -x,-z,-y
                            M(1:3,1) = -(/ mtmp(1,1),mtmp(2,1),mtmp(3,1) /)    
                            M(1:3,2) = -(/ mtmp(1,3),mtmp(2,3),mtmp(3,3) /)
                            M(1:3,3) = -(/ mtmp(1,2),mtmp(2,2),mtmp(3,2) /)
                        end if
                    else
                        !   largest x component in z-direction
                        if (mtmp(1,3) > 0) then
                            !   z,-x,-y
                            M(1:3,1) =  (/ mtmp(1,3),mtmp(2,3),mtmp(3,3) /)    
                            M(1:3,2) = -(/ mtmp(1,1),mtmp(2,1),mtmp(3,1) /)
                            M(1:3,3) = -(/ mtmp(1,2),mtmp(2,2),mtmp(3,2) /)                            
                        else
                            !  -z, x,-y
                            M(1:3,1) = -(/ mtmp(1,3),mtmp(2,3),mtmp(3,3) /)    
                            M(1:3,2) =  (/ mtmp(1,1),mtmp(2,1),mtmp(3,1) /)
                            M(1:3,3) = -(/ mtmp(1,2),mtmp(2,2),mtmp(3,2) /)
                        end if
                    end if
                end if
                 
            else          
                !   largest z component in x direction
                !   *,*,x or *,*,-x
                 if (mtmp(3,1) > 0) then
                    !   ... and z component is +ve direction.
                    !   *,*,x                     
                    if (abs(mtmp(1,2))>=abs(mtmp(1,3))) then
                        !   largest x component in y-direction
                        if (mtmp(1,2) > 0) then
                            !   y,z,x
                            M(1:3,1) =  (/ mtmp(1,2),mtmp(2,2),mtmp(3,2) /)    
                            M(1:3,2) =  (/ mtmp(1,3),mtmp(2,3),mtmp(3,3) /)
                            M(1:3,3) =  (/ mtmp(1,1),mtmp(2,1),mtmp(3,1) /)
                        else
                            !  -y,-z,x
                            M(1:3,1) = -(/ mtmp(1,2),mtmp(2,2),mtmp(3,2) /)    
                            M(1:3,2) = -(/ mtmp(1,3),mtmp(2,3),mtmp(3,3) /)
                            M(1:3,3) =  (/ mtmp(1,1),mtmp(2,1),mtmp(3,1) /)
                        end if
                    else
                        !   largest x component in z-direction
                        if (mtmp(1,3) > 0) then
                            !   z,-y,x
                            M(1:3,1) =  (/ mtmp(1,3),mtmp(2,3),mtmp(3,3) /)    
                            M(1:3,2) = -(/ mtmp(1,2),mtmp(2,2),mtmp(3,2) /)
                            M(1:3,3) =  (/ mtmp(1,1),mtmp(2,1),mtmp(3,1) /)                      
                        else
                            !  -z,y,x
                            M(1:3,1) = -(/ mtmp(1,3),mtmp(2,3),mtmp(3,3) /)    
                            M(1:3,2) =  (/ mtmp(1,2),mtmp(2,2),mtmp(3,2) /)
                            M(1:3,3) =  (/ mtmp(1,1),mtmp(2,1),mtmp(3,1) /)
                        end if
                    end if
                else
                    !   *,*,-x                                     
                    if (abs(mtmp(1,2))>=abs(mtmp(1,3))) then
                        !   largest x component in y-direction
                        if (mtmp(1,2) > 0) then
                            !   y,-z,-x
                            M(1:3,1) =  (/ mtmp(1,2),mtmp(2,2),mtmp(3,2) /)    
                            M(1:3,2) = -(/ mtmp(1,3),mtmp(2,3),mtmp(3,3) /)
                            M(1:3,3) = -(/ mtmp(1,1),mtmp(2,1),mtmp(3,1) /)    
                        else
                            !  -y, z,-x
                            M(1:3,1) = -(/ mtmp(1,2),mtmp(2,2),mtmp(3,2) /)    
                            M(1:3,2) =  (/ mtmp(1,3),mtmp(2,3),mtmp(3,3) /)
                            M(1:3,3) = -(/ mtmp(1,1),mtmp(2,1),mtmp(3,1) /)
                        end if
                    else
                        !   largest x component in z-direction
                        if (mtmp(1,3) > 0) then
                            !   z,y,-x
                            M(1:3,1) =  (/ mtmp(1,3),mtmp(2,3),mtmp(3,3) /)    
                            M(1:3,2) =  (/ mtmp(1,2),mtmp(2,2),mtmp(3,2) /)
                            M(1:3,3) = -(/ mtmp(1,1),mtmp(2,1),mtmp(3,1) /)                            
                        else
                            !  -z,-y,-x
                            M(1:3,1) = -(/ mtmp(1,3),mtmp(2,3),mtmp(3,3) /)    
                            M(1:3,2) = -(/ mtmp(1,2),mtmp(2,2),mtmp(3,2) /)
                            M(1:3,3) = -(/ mtmp(1,1),mtmp(2,1),mtmp(3,1) /)
                        end if
                    end if
                end if
                 
            end if
                
      
            

            if (present(norescale)) then
                if (norescale) return
            end if
             
            call completeBasis( M(1:3,3),M(1:3,1),M(1:3,2) )
            call restoreUnitaryness( M )
            
            return
        end subroutine tidyRotationMatrix
        
        pure subroutine restoreUnitaryness( M )
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      make sure that M is a unitary matrix
            real(kind=real64),dimension(3,3),intent(inout)    ::  M
            real(kind=real64)       ::      dd
            
            dd = determinant3Mat(M)
            dd = dd**(-1/3.0d0)
            M = M*dd                
            
            return
        end subroutine restoreUnitaryness
        
        pure function inverseRotationMatrix(M) result(iM)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      note that the inverse of a rotation matrix is its transpose.
            real(kind=real64),dimension(3,3),intent(in)       ::  M
            real(kind=real64),dimension(3,3)                  ::  iM
            iM(1:3,1:3) = transpose( M(1:3,1:3) )
            return
        end function inverseRotationMatrix

       
        pure function rotateVector0(R,x) result(y)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      return y = Rx
            real(kind=real64),dimension(3,3),intent(in)     ::  R
            real(kind=real64),dimension(3),intent(in)       ::  x
            real(kind=real64),dimension(3)                  ::  y
            y(1:3) = R(1:3,1)*x(1) + R(1:3,2)*x(2) + R(1:3,3)*x(3) 
            return
        end function rotateVector0    
        
        pure function rotateVector1(R,x) result(y)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      return y = Rx
            real(kind=real64),dimension(3,3),intent(in)     ::  R
            real(kind=real64),dimension(:,:),intent(in)     ::  x
            real(kind=real64),dimension(3,size(x,dim=2))    ::  y
            integer         ::      ii
            do ii = 1,size(x,dim=2)
                y(1:3,ii) = R(1:3,1)*x(1,ii) + R(1:3,2)*x(2,ii) + R(1:3,3)*x(3,ii) 
            end do
            return
        end function rotateVector1    
        
     
        subroutine completeBasis0( z,x,y )
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      given the z-axis vector, complete the basis to provide a triplet x,y,n
    !*      if on input x is set, then use this as a hint to attempt to place x along this direction
            real(kind=real64),dimension(3),intent(in)       ::      z
            real(kind=real64),dimension(3),intent(inout)    ::      x
            real(kind=real64),dimension(3),intent(out)      ::      y

            real(kind=real64)                   ::      xxxx    !   vector lengths
            real(kind=real64)                   ::      zdotx
            
        
        !---    check for a sensible hint for the x-direction
            xxxx = norm2(x)
            zdotx = dot_product( z,x )
            if ( (xxxx < 0.001d0 ).or.(abs(zdotx) >= xxxx*0.999d0) ) then
                !   haven't got a good hint. Make random hint
                if (abs(z(3))>max(abs(z(1)),abs(z(2)))) then
                    !   z points along 3-axis
                    x(1) = 1.0d0
                    x(3) = 0.0d0
                    if (abs(z(1))>0) then
                        x(2) = - z(2)/z(1)
                    else
                        x(2) = 0.0d0
                    end if
                    zdotx = dot_product( z,x )
                else if (abs(z(2))>max(abs(z(1)),abs(z(3)))) then
                    !   z points along 2-axis
                    x(3) = 1.0d0
                    x(2) = 0.0d0
                    if (abs(z(3))>0) then
                        x(1) = - z(1)/z(3)
                    else
                        x(1) = 0.0d0
                    end if
                    zdotx = dot_product( z,x )
                else
                    !   z points along 1-axis
                    x(2) = 1.0d0
                    x(1) = 0.0d0
                    if (abs(z(2))>0) then
                        x(3) = - z(3)/z(2)
                    else
                        x(3) = 0.0d0
                    end if
                    zdotx = dot_product( z,x )
                end if
            end if


        !---    remove projection of z on x-direction
            x = x - z*zdotx
            xxxx = norm2(x)
            x = x/xxxx

        !---    construct y-direction
            y = cross_product( z,x )

            return
        end subroutine completeBasis0
                          
        
        pure function cross_product(x,y) result(z)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
            real(kind=real64),dimension(3),intent(in)       ::      x,y
            real(kind=real64),dimension(3)                  ::      z
            z(1) = x(2)*y(3) - x(3)*y(2)
            z(2) = x(3)*y(1) - x(1)*y(3)
            z(3) = x(1)*y(2) - x(2)*y(1)
            return
        end function cross_product
        

        pure function determinant3Mat(M) result(d)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      returns the determinant of M
            real(kind=real64),dimension(3,3),intent(in)      ::      M
            real(kind=real64)                                ::      d
            real(kind=real64),dimension(9)       ::      dd
            dd(1) = M(1,1)*( M(2,2)*M(3,3) - M(2,3)*M(3,2) )
            dd(2) = M(1,2)*( M(2,3)*M(3,1) - M(2,1)*M(3,3) )
            dd(3) = M(1,3)*( M(2,1)*M(3,2) - M(2,2)*M(3,1) )
            dd(4) = M(2,1)*( M(3,2)*M(1,3) - M(3,3)*M(1,2) )
            dd(5) = M(2,2)*( M(3,3)*M(1,1) - M(3,1)*M(1,3) )
            dd(6) = M(2,3)*( M(3,1)*M(1,2) - M(3,2)*M(1,1) )
            dd(7) = M(3,1)*( M(1,2)*M(2,3) - M(1,3)*M(2,2) )
            dd(8) = M(3,2)*( M(1,3)*M(2,1) - M(1,1)*M(2,3) )
            dd(9) = M(3,3)*( M(1,1)*M(2,2) - M(1,2)*M(2,1) )
            d = (1.0d0/3.0d0) * sum(dd)
            return
        end function determinant3Mat
        
    !---
        
        subroutine report0(this,u,o)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      simple output. Defaults to unit 6 = screen.
    !*      optional argument o determines left hand margin
            real(kind=real64),dimension(3,3),intent(in)    ::      this
            integer,intent(in),optional     ::      u,o
            integer     ::      uu,oo
            uu = 6 ; if (present(u)) uu = u
            oo = 0 ; if (present(o)) oo = o
            write(unit=uu,fmt='(a)') repeat(" ",oo)//"RotationMatrix"            
            write(unit=uu,fmt='(a,3f12.6)') repeat(" ",oo+2),this(1,1:3)
            write(unit=uu,fmt='(a,3f12.6)') repeat(" ",oo+2),this(2,1:3)
            write(unit=uu,fmt='(a,3f12.6)') repeat(" ",oo+2),this(3,1:3)
            return
        end subroutine report0
    
    
        function rndRotMat0(x) result(R)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      set up randomly oriented rotation matrix by setting a random z and x axis direction
    !*      if x present, then make z = (0,0,1+/-x) and x = (1+/-x,0,0)
            real(kind=real64),intent(in),optional   ::      x
            
            real(kind=real64),dimension(3,3)    ::      R
            real(kind=real64),dimension(3)      ::      vv,vabs
            real(kind=real64)                   ::      dd
            if (present(x)) then
            
                !do 
                !    vv = rndUnitVec() 
                !    vabs = abs(vv)
                !    dd = maxval( vabs )
                !    if (dd > 1-x) exit
                !end do
                !if ( dd >= max(vabs(2),vabs(3)) ) then
                !    R(1,1) = dd ; R(2,1) = vv(2) ; R(3,1) = vv(3)
                !else if (dd >= max(vabs(1),vabs(3)) ) then
                !    R(1,1) = dd ; R(2,1) = vv(1) ; R(3,1) = vv(3)
                !else 
                !    R(1,1) = dd ; R(2,1) = vv(1) ; R(3,1) = vv(2)
                !end if
                
                call random_number(dd)
                dd = (2*dd - 1)*x
                R(1:3,1) = (/ 1 + dd,0.0d0,0.0d0 /)
                
                do 
                    vv = rndUnitVec()
                    vabs = abs(vv) 
                    dd =  maxval( vabs )
                    if (dd > 1-x) exit
                end do
                
                if (vabs(1) >= max(vabs(2),vabs(3)) ) then
                    R(1,3) = vv(2) ; R(2,3) = vv(3) ; R(3,3) = vabs(1)
                else if (vabs(2) >= max(vabs(1),vabs(3)) ) then
                    R(1,3) = vv(1) ; R(2,3) = vv(3) ; R(3,3) = vabs(2) 
                else 
                    R(1,3) = vv(1) ; R(2,3) = vv(2) ; R(3,3) = vabs(3) 
                end if
               ! print *," rnd ",R(1:3,1)," ,d ",dd," z ",R(1:3,3),(dd >= max(vabs(2),vabs(3)) ),(dd >= max(vabs(1),vabs(3)) )," v ",vv," vabs ",vabs
                call completeBasis( R(1:3,3),R(1:3,1),R(1:3,2) )
                return
            end if
            
        
            R(1:3,1) = rndUnitVec() 
            do 
                R(1:3,3) = rndUnitVec() 
                if (abs( dot_product( R(1:3,3),R(1:3,1) ) )<0.9 ) exit
            end do
            call completeBasis( R(1:3,3),R(1:3,1),R(1:3,2) )
            
            return                       
        end function rndRotMat0
        
        function rndUnitVec0() result(x)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
            real(kind=real64),dimension(3)  ::      x
            real(kind=real64)               ::      dd
            do 
                call random_number(x)
                x = 2*x - 1
                dd = x(1)*x(1) + x(2)*x(2) + x(3)*x(3)
                if (dd*(1-dd) > 0.0d0) exit

            end do 
            x = x/sqrt(dd)
            
            return
        end function rndUnitVec0 
        
        

        function rndRotMat1(seed,x) result(R)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      set up randomly oriented rotation matrix by setting a random z and x axis direction
    !*      if x present, then make z = (0,0,1+/-x) and x = (1+/-x,0,0)
            integer,intent(inout)                   ::      seed
            real(kind=real64),intent(in),optional   ::      x
            
            real(kind=real64),dimension(3,3)    ::      R
            real(kind=real64),dimension(3)      ::      vv,vabs
            real(kind=real64)                   ::      dd
            if (present(x)) then
            
                !do 
                !    vv = rndUnitVec() 
                !    vabs = abs(vv)
                !    dd = maxval( vabs )
                !    if (dd > 1-x) exit
                !end do
                !if ( dd >= max(vabs(2),vabs(3)) ) then
                !    R(1,1) = dd ; R(2,1) = vv(2) ; R(3,1) = vv(3)
                !else if (dd >= max(vabs(1),vabs(3)) ) then
                !    R(1,1) = dd ; R(2,1) = vv(1) ; R(3,1) = vv(3)
                !else 
                !    R(1,1) = dd ; R(2,1) = vv(1) ; R(3,1) = vv(2)
                !end if
                
                dd = ran0(seed)
                dd = (2*dd - 1)*x
                R(1:3,1) = (/ 1 + dd,0.0d0,0.0d0 /)
                
                do 
                    vv = rndUnitVec(seed)
                    vabs = abs(vv) 
                    dd =  maxval( vabs )
                    if (dd > 1-x) exit
                end do
                
                if (vabs(1) >= max(vabs(2),vabs(3)) ) then
                    R(1,3) = vv(2) ; R(2,3) = vv(3) ; R(3,3) = vabs(1)
                else if (vabs(2) >= max(vabs(1),vabs(3)) ) then
                    R(1,3) = vv(1) ; R(2,3) = vv(3) ; R(3,3) = vabs(2) 
                else 
                    R(1,3) = vv(1) ; R(2,3) = vv(2) ; R(3,3) = vabs(3) 
                end if
               ! print *," rnd ",R(1:3,1)," ,d ",dd," z ",R(1:3,3),(dd >= max(vabs(2),vabs(3)) ),(dd >= max(vabs(1),vabs(3)) )," v ",vv," vabs ",vabs
                call completeBasis( R(1:3,3),R(1:3,1),R(1:3,2) )
                return
            end if
            
        
            R(1:3,1) = rndUnitVec() 
            do 
                R(1:3,3) = rndUnitVec(seed) 
                if (abs( dot_product( R(1:3,3),R(1:3,1) ) )<0.9 ) exit
            end do
            call completeBasis( R(1:3,3),R(1:3,1),R(1:3,2) )
            
            return                       
        end function rndRotMat1
        
        function rndUnitVec1(seed) result(x)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
            integer,intent(inout)           ::      seed
            real(kind=real64),dimension(3)  ::      x
            
            real(kind=real64)               ::      dd
            do 
                x = ran0(seed)
                x = 2*x - 1
                dd = x(1)*x(1) + x(2)*x(2) + x(3)*x(3)
                if (dd*(1-dd) > 0.0d0) exit

            end do 
            x = x/sqrt(dd)
            
            return
        end function rndUnitVec1     
            
    end module Lib_RotationMatrices