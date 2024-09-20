
    module Lib_Quaternions
!---^^^^^^^^^^^^^^^^^^^^^^
!*      utility module to privide functionality for quaternions
!*      see eg https://en.wikipedia.org/wiki/Quaternion
!*  
!*      Daniel Mason, UKAEA
!*      April 2022
!*
!*      My quaternions are defined as (/ w, x,y,z /)


        use iso_fortran_env
        implicit none
        private
        
        external    ::      DSYEV                       !   note: to convert rotation matrix to quaternion robustly requires an eigendecomposition.
                                                        
    !---                                                
                                                        
        public      ::      Quaternion_ctor             !   construct a quaternion object
        public      ::      delete                      !   no dynamic memory
        public      ::      report                      !   simple dump 
     
        public      ::      normalise                   !   convert to unit quaternion
        public      ::      norm                        !   convert to unit quaternion
        public      ::      operator(*)                 !   multiply two quaternions                                
        public      ::      operator(+)                 !   add two quaternions                                
        public      ::      operator(-)                 !   subtract two quaternions                                
        public      ::      operator(.dot.)             !   dot product
        public      ::      distance                    !   quaternion distance measure - note not an angular measure
        public      ::      quaternionDistanceInDegrees !   convert distance to great circle rotation
        public      ::      degreesAsQuaternionDistance !   convert great circle rotation to distance
        public      ::      inverse                     !   inverse quaternion
        public      ::      quaternionToRotMat          !   convert to 3x3 rotation matrix
        public      ::      rotateVector                !   rotate vector
        public      ::      getComponents               !   simple get operation - returns (/ w, x,y,z /)
        public      ::      imposeSymmetry              !   reduce rotation to nearest in symmetric group
        public      ::      getNSymmetries              !   return number of symmetries in group
        public      ::      symmetryRelatedQuaternion   !   return one of the symmetry related quaternions 
        public      ::      midpoint                    !   midpoint between two quaternions
        public      ::      slerp                       !   interpolation between two quaternions    
        
    !---
    
        logical,public          ::      Quaternion_dbg = .false.
        
    !---     
    
        real(kind=real64),private,parameter     ::      PI = 3.141592653590d0
        
    !---   
        
    
        
        integer,public,parameter        ::      QUATERNION_SYMMETRY_NONE = 0        
        integer,public,parameter        ::      QUATERNION_SYMMETRY_CUBIC = 1                  
        integer,public,parameter        ::      QUATERNION_SYMMETRY_HEX = 2                  
        
    !---
    
        type,public     ::      Quaternion     
            real(kind=real64)           ::      w , x,y,z
        end type Quaternion
        
    !---  
         
        include     "cubic_rot_q_data.h"
        include     "hex_rot_q_data.h"
        type(Quaternion),public,parameter       ::      Quaternion_identity = Quaternion(1.0d0,0.0d0,0.0d0,0.0d0)
        
    !---
    
        
        interface Quaternion_ctor
            module procedure    Quaternion_null
            module procedure    Quaternion_ctor0
            module procedure    Quaternion_ctor1
            module procedure    RotMatToQuaternion
        end interface
                
        interface delete
            module procedure    delete0
        end interface
        
        interface report
            module procedure    report0
            module procedure    report1
        end interface
        
        interface normalise
            module procedure    normalise0
        end interface
        
        interface norm
            module procedure    norm0
        end interface
        
        interface operator(*)
            module procedure    times
        end interface
        
        interface operator(+)
            module procedure    add
        end interface

        interface operator(-)
            module procedure    minus
        end interface
        
        interface operator(.dot.)
            module procedure    dotProduct
        end interface
        
        interface distance
            module procedure    quaternionDistance0
            module procedure    quaternionDistance1
            module procedure    quaternionDistance2
        end interface
        
        interface inverse
            module procedure    inverseQuaternion
        end interface
                                    
        interface rotateVector
            module procedure    rotateVector0
        end interface          
        
        interface imposeSymmetry
            module procedure    imposeSymmetry0
        end interface
        
    contains
!---^^^^^^^^

        pure function Quaternion_null() result(this)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      returns a unit quaternion equivalent to identity rotation
            type(Quaternion)           ::      this
            this%w = 1
            this%x = 0
            this%y = 0
            this%z = 0
            return
        end function Quaternion_null
                         
        pure function Quaternion_ctor0(w,x,y,z) result(this)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      return a quaternion object from the four components in order w , x,y,z
            type(Quaternion)                ::      this
            real(kind=real64),intent(in)    ::      w, x,y,z 
            this%w = w
            this%x = x
            this%y = y
            this%z = z                        
            return
        end function Quaternion_ctor0
                       
        pure function Quaternion_ctor1(x) result(this)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      return a quaternion object from the four components in order (/w , x,y,z/)
            type(Quaternion)           ::      this
            real(kind=real64),dimension(4),intent(in)     ::      x       
            this%w = x(1)
            this%x = x(2)
            this%y = x(3)
            this%z = x(4)                        
            return
        end function Quaternion_ctor1
    !---
    
        subroutine delete0(this)
    !---^^^^^^^^^^^^^^^^^^^^^^^^
    !*      note: no dynamic memory
            type(Quaternion),intent(inout)    ::      this            
            this = Quaternion_null()
            return
        end subroutine delete0
        
    !---
    
        subroutine report0(this,u,o)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      simple output. Defaults to unit 6 = screen.
    !*      optional argument o determines left hand margin
            type(Quaternion),intent(in)    ::      this
            integer,intent(in),optional     ::      u,o
            integer     ::      uu,oo
            uu = 6 ; if (present(u)) uu = u
            oo = 0 ; if (present(o)) oo = o
            write(unit=uu,fmt='(5(a,f16.12))') repeat(" ",oo)//"Quaternion [w,x,y,z=",this%w,",",this%x,",",this%y,",",this%z,"]"    
            return
        end subroutine report0
    
        subroutine report1(this,asRotMat,u,o)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      alternate output as a rotation matrix. Defaults to unit 6 = screen.
    !*      optional argument o determines left hand margin    
            type(Quaternion),intent(in)    ::      this
            logical,intent(in)              ::      asRotMat
            integer,intent(in),optional     ::      u,o
            integer     ::      uu,oo
            real(kind=real64),dimension(3,3)    ::      RR
            uu = 6 ; if (present(u)) uu = u
            oo = 0 ; if (present(o)) oo = o
            call report0(this,uu,oo)
            if (.not. asRotMat) return
            RR = quaternionToRotMat(this)
            write(unit=uu,fmt='(a,3f16.8)') repeat(" ",oo+4),RR(1,:)
            write(unit=uu,fmt='(a,3f16.8)') repeat(" ",oo+4),RR(2,:)
            write(unit=uu,fmt='(a,3f16.8)') repeat(" ",oo+4),RR(3,:)
            
            return
        end subroutine report1
    
    !---
        
        elemental subroutine normalise0(this)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      ensure this is a unit quaternion
            type(Quaternion),intent(inout)    ::      this
            real(kind=real64)       ::      tt
            tt = this%w*this%w + this%x*this%x + this%y*this%y + this%z*this%z
            if (tt == 0) then
                this = Quaternion_ctor()
            else
                tt = 1/sqrt(tt)
                this%w = this%w * tt
                this%x = this%x * tt
                this%y = this%y * tt
                this%z = this%z * tt
            end if
            return
        end subroutine normalise0

        elemental function norm0(that) result(this)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      ensure this is a unit quaternion
            type(Quaternion),intent(in)         ::      that
            type(Quaternion)                    ::      this
            real(kind=real64)       ::      tt
            tt = that%w*that%w + that%x*that%x + that%y*that%y + that%z*that%z
            if (tt == 0) then
                this = Quaternion_ctor()
            else
                tt = 1/sqrt(tt)
                this%w = that%w * tt
                this%x = that%x * tt
                this%y = that%y * tt
                this%z = that%z * tt
            end if
            return
        end function norm0
    
        pure function times(a,b) result(this)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      this = a * b
            type(Quaternion),intent(in)     ::      a,b
            type(Quaternion)                ::      this
            this%w = a%w*b%w - a%x*b%x - a%y*b%y - a%z*b%z
            this%x = a%w*b%x + a%x*b%w + a%y*b%z - a%z*b%y
            this%y = a%w*b%y - a%x*b%z + a%y*b%w + a%z*b%x
            this%z = a%w*b%z + a%x*b%y - a%y*b%x + a%z*b%w 
            return
        end function times
        
        pure function add(a,b) result(this)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      this = a + b
            type(Quaternion),intent(in)     ::      a,b
            type(Quaternion)                ::      this
            this%w = a%w + b%w 
            this%x = a%x + b%x 
            this%y = a%y + b%y 
            this%z = a%z + b%z 
            return
        end function add

        pure function minus(a,b) result(this)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      this = a - b
            type(Quaternion),intent(in)     ::      a,b
            type(Quaternion)                ::      this
            this%w = a%w - b%w 
            this%x = a%x - b%x 
            this%y = a%y - b%y 
            this%z = a%z - b%z 
            return
        end function minus

        elemental function dotProduct(a,b) result(d)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      d = a . b*
            type(Quaternion),intent(in)     ::      a,b
            real(kind=real64)               ::      d
            d = a%w*b%w + a%x*b%x + a%y*b%y + a%z*b%z
            return
        end function dotProduct
        
        
    
    !---
    
        function RotMatToQuaternion(R) result(q)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      return a quaternion q = (/ w, x,y,z  /)
    !*      https://d3cw3dd2w32x2b.cloudfront.net/wp-content/uploads/2015/01/matrix-to-quat.pdf
    !*      https://en.wikipedia.org/wiki/Quaternions_and_spatial_rotation
            real(kind=real64),dimension(3,3),intent(in)         ::      R
            type(Quaternion)                                    ::      q
        
            real(kind=real64),dimension(4,4)        ::      KK
            real(kind=real64),dimension(4)          ::      ww
            real(kind=real64),dimension(4*66)       ::      work
            integer                                 ::      ii
            
            KK = 0
            KK(1,1) = R(1,1) - R(2,2) - R(3,3)
            KK(2,1) = R(2,1) + R(1,2)
            KK(3,1) = R(3,1) + R(1,3)
            KK(4,1) = R(2,3) - R(3,2)
            
            KK(2,2) = R(2,2) - R(1,1) - R(3,3)
            KK(3,2) = R(3,2) + R(2,3)
            KK(4,2) = R(3,1) - R(1,3)
            
            KK(3,3) = R(3,3) - R(1,1) - R(2,2)
            KK(4,3) = R(1,2) - R(2,1)
            
            KK(4,4) = R(1,1) + R(2,2) + R(3,3)                      
            
            call DSYEV( "V","L",4,KK,4,ww,work,size(work),ii )
                       
            if (KK(4,4)>0) KK(1:4,4) = -KK(1:4,4)       !   because this is an eigenvalue problem, I am free to take the -ve value. Helps with slerp. 
            
            
            q%x = KK(1,4)
            q%y = KK(2,4)
            q%z = KK(3,4)
            q%w = -KK(4,4)
            
            
                
            
            return
        end function RotMatToQuaternion
                    
        pure function quaternionToRotMat(q) result(R)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      https://en.wikipedia.org/wiki/Quaternions_and_spatial_rotation
            type(Quaternion),intent(in)                 ::      q
            real(kind=real64),dimension(3,3)            ::      R
        
            real(kind=real64)       ::      ss 
              
            ss = ( q%x*q%x + q%y*q%y + q%z*q%z + q%w*q%w )/2
            
            R(1,1) = ss - ( q%y*q%y + q%z*q%z )
            R(2,1) =      ( q%x*q%y + q%z*q%w )
            R(3,1) =      ( q%x*q%z - q%y*q%w )
            R(1,2) =      ( q%x*q%y - q%z*q%w )
            R(2,2) = ss - ( q%x*q%x + q%z*q%z )
            R(3,2) =      ( q%y*q%z + q%x*q%w )
            R(1,3) =      ( q%x*q%z + q%y*q%w )
            R(2,3) =      ( q%y*q%z - q%x*q%w )
            R(3,3) = ss - ( q%x*q%x + q%y*q%y )
                        
            ss = 1/ss
            R(1:3,1:3) = R(1:3,1:3) * ss

            !if (determinant3Mat(R)<0) R = -R
                        
            return
        end function quaternionToRotMat
        
    !    
        
        
        
        
    !
    
        elemental function quaternionDistanceInDegrees( d ) result(theta)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      returns distance as an angle in degrees
            real(kind=real64),intent(in)            ::      d
            real(kind=real64)                       ::      theta
            theta = acos(1-2*d)
            theta = theta*180.0d0/PI
            return
        end function quaternionDistanceInDegrees
    
        elemental function degreesAsQuaternionDistance( theta ) result(d )
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      returns angle in degrees as a distance
            real(kind=real64),intent(in)            ::      theta
            real(kind=real64)                       ::      d
            d = theta*PI/180.0d0
            d = (1 - cos(d))/2
            return
        end function degreesAsQuaternionDistance
    
            
            
        elemental function quaternionDistance0(q1,q2,normalise) result(d)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      returns a distance metric for quaternions
    !*          d = 1 - |q1.q2|^2 = (1-cos theta)/2
    !*      if normalise, then ensure q1 and q2 are unit quaternions first.
    !*      returns d=0 for q1 = q2 and d=1 for orientations 180 apart
            type(Quaternion),intent(in)                 ::      q1,q2
            logical,intent(in),optional                 ::      normalise
            real(kind=real64)                           ::      d 
        
            real(kind=real64)       ::      q1dotq1,q2dotq2,q1dotq2 
                          
            q1dotq2 = q1 .dot. q2
            d = q1dotq2*q1dotq2
            
            if (present(normalise)) then
                if (normalise) then
                    q1dotq1 = q1 .dot. q1
                    q2dotq2 = q2 .dot. q2
                    if (q1dotq1*q2dotq2 == 0) then
                        d = 1 ; return
                    end if
                    d = d/(q1dotq1*q2dotq2)
                end if
            end if
            
            d = 1 - d
         
            return
        end function quaternionDistance0
        

        function quaternionDistance1(q1,q2,symmetry) result(d)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      returns a distance metric for quaternions
    !*          d = 1 - |q1.q2|^2 = (1-cos theta)/2  
    !*      taking into account all rotations in the symmetry group  
            type(Quaternion),intent(in)         ::      q1,q2
            integer,intent(in)                  ::      symmetry
            real(kind=real64)                   ::      d 
            integer                             ::      ss
            type(Quaternion)                    ::      qs 
            type(Quaternion),dimension(:),pointer   ::      qsym
            real(kind=real64)                   ::      q1dotq2,q1dotq2best
            
            select case(symmetry)
                case(QUATERNION_SYMMETRY_CUBIC)
                    qsym => cubic_rot_q                    
                case(QUATERNION_SYMMETRY_HEX)
                    qsym => hex_rot_q                    
                case default
                    d = quaternionDistance0(q1,q2)
                    return
            end select

           ! print *,"q1"            
           ! call report(q1)
           ! write (*,fmt='(9f12.4)') QuaternionToRotMat(q1)
           ! print *,"q2"
           ! call report(q2)
           ! write (*,fmt='(9f12.4)') QuaternionToRotMat(q2)
            q1dotq2best = 0.0d0
            do ss = 1,size(qsym)
                qs = q2 * qsym(ss)
                !call report(qs)
                !write (*,fmt='(9f12.4)') QuaternionToRotMat(qs)
                q1dotq2 = q1 .dot. qs
                q1dotq2best = max(q1dotq2best,q1dotq2*q1dotq2)
                !print *,"qdist ",ss,1-q1dotq2*q1dotq2
            end do
             
            d = 1 - q1dotq2best
            
            
            return
        end function quaternionDistance1
                
        elemental function quaternionDistance2(q1) result(d)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      returns a distance metric for a unit quaternion from the identity
    !*          d = 1 - |q1.q2|^2 = (1-cos theta)/2
    !*      note  a.b = a%w*b%w + a%x*b%x + a%y*b%y + a%z*b%z
    !*      and Quaternion_identity = Quaternion(1.0d0,0.0d0,0.0d0,0.0d0)
            type(Quaternion),intent(in)                 ::      q1            
            real(kind=real64)                           ::      d 
        
            d = 1 - q1%w*q1%w
         
            return
        end function quaternionDistance2
        
            
    !---                        
        
        elemental function inverseQuaternion(that) result(this)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      this * that = that * this = unit
            type(Quaternion),intent(in)                 ::      that
            type(Quaternion)                            ::      this
            real(kind=real64)       ::      tt
            tt = ( that%x*that%x + that%y*that%y + that%z*that%z + that%w*that%w )
            if (tt==0) then
                this = Quaternion_ctor()
            else
                tt = 1/sqrt(tt)
                this%x = - that%x * tt 
                this%y = - that%y * tt
                this%z = - that%z * tt
                this%w =   that%w * tt
            end if
            return
        end function inverseQuaternion
        
        pure function rotateVector0(this,a) result(b)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      https://en.wikipedia.org/wiki/Quaternions_and_spatial_rotation
    !*      b = q a q^-1
            type(Quaternion),intent(in)                 ::      this
            real(kind=real64),dimension(3),intent(in)   ::      a
            real(kind=real64),dimension(3)              ::      b
        
            real(kind=real64)       ::      ss 
              
            
            ss = ( this%x*this%x + this%y*this%y + this%z*this%z + this%w*this%w )/2
            
            b(1:3) = ss*a(1:3) 
            
            b(1) = b(1) - ( this%y*this%y + this%z*this%z )*a(1)    &
                        + ( this%x*this%y - this%z*this%w )*a(2)    &
                        + ( this%x*this%z + this%y*this%w )*a(3)
            
            b(2) = b(2) + ( this%x*this%y + this%z*this%w )*a(1)    &
                        - ( this%x*this%x + this%z*this%z )*a(2)    &
                        + ( this%y*this%z - this%x*this%w )*a(3)
                        
            b(3) = b(3) + ( this%x*this%z - this%y*this%w )*a(1)    &
                        + ( this%y*this%z + this%x*this%w )*a(2)    &
                        - ( this%x*this%x + this%y*this%y )*a(3)
            
            
            ss = 1/ss
            b(1:3) = b(1:3) * ss
            
            return
        end function rotateVector0
            
        
        pure function getComponents(this) result(x)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      return components of the quaternion as a vector (/ w , x,y,z /)
            type(Quaternion),intent(in)     ::      this
            real(kind=real64),dimension(4)  ::      x
            x = (/ this%w,this%x,this%y,this%z /)
            return
        end function getComponents
        
        pure function midpoint(q1,q2) result(qm)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^  
    !*      find the midpoint between two quaternions          
            type(Quaternion),intent(in)     ::      q1,q2
            type(Quaternion)                ::      qm
            real(kind=real64)       ::      dd
            
            qm%w = q1%w + q2%w
            qm%x = q1%x + q2%x
            qm%y = q1%y + q2%y
            qm%z = q1%z + q2%z
            dd = qm%w*qm%w + qm%x*qm%x + qm%y*qm%y + qm%z*qm%z
            if (dd>0) then
                dd = 1/sqrt(dd)
                qm%w = qm%w * dd
                qm%x = qm%x * dd
                qm%y = qm%y * dd
                qm%z = qm%z * dd
            end if
            return
        end function midpoint
            
        pure function slerp(q1,q2,x) result(qx)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^     
    !*      return q1(1-x) + q2 x
    !*      spherical linear interpolation       
            type(Quaternion),intent(in)     ::      q1,q2
            real(kind=real64),intent(in)    ::      x
            type(Quaternion)                ::      qx
            real(kind=real64)       ::      omega
            real(kind=real64)       ::      x1,x2
            
            x1 = q1 .dot. q2
            if (x1 >= 1.0d0) then
                !   the angle between the quaternions is zero, which means they are the same, which means the interpolation is either end
                qx = q1
                return
            else if (x1 <= -1.0d0) then
                !   the angle between the quaternions is pi, which means they are antiparallel, which means the interpolation is undefined. 
                !   pick the midpoint.                
                qx = midpoint(q1,q2)
                return
            else
                x2 = 2*x1*x1 - 1
                omega = acos( x2 )     !   angle subtended by arc in 4D
            end if
             
            
            x2 = 1/sin(omega)            
            x1 = sin( (1-x)*omega )*x2
            x2 = sin( x*omega )*x2
            
            qx%w = x1*q1%w + x2*q2%w
            qx%x = x1*q1%x + x2*q2%x
            qx%y = x1*q1%y + x2*q2%y
            qx%z = x1*q1%z + x2*q2%z
            
            call normalise(qx)
            
            return 
        end function slerp
            

        subroutine imposeSymmetry0( this,symmetry, q_best_sym )
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      try rotating the quaternion through all options in the symmetry group
    !*      and find the one closest to the unit quaternion
    !*      optionally return the symmetrising quaternion such that this_output = this_input * q_best_sym
            type(Quaternion),intent(inout)      ::      this
            integer,intent(in)                  ::      symmetry
            type(Quaternion),intent(out),optional   ::      q_best_sym
            integer                                 ::      ss
            real(kind=real64)                       ::      dd,dbest
            type(Quaternion)                        ::      qs,qbest
            type(Quaternion),dimension(:),pointer   ::      qsym
             
            select case(symmetry)
                
                case(QUATERNION_SYMMETRY_CUBIC)
                    qsym => cubic_rot_q       
                case(QUATERNION_SYMMETRY_HEX)
                    qsym => hex_rot_q       
                case default
                    return
                                 
            end select
             
            
            dbest = huge(1.0)   
            do ss = 1,size(qsym)
                qs = this * qsym(ss) 
                dd = distance( qs )
                 
                if (dd < dbest) then
                    dbest = dd
                    qbest = qs
                    if (present(q_best_sym)) q_best_sym = qsym(ss) 
                end if
            end do                

            this = qbest
             
            
            return
        end subroutine imposeSymmetry0
        
         
        pure function getNSymmetries(symmetry) result(n)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ 
    !*      return number of rotation symmetries in the group
            integer,intent(in)                  ::      symmetry
            integer                             ::      n
            select case(symmetry)                
                case(QUATERNION_SYMMETRY_CUBIC)
                    n = 24
                case(QUATERNION_SYMMETRY_HEX)
                    n = 6
                case default
                    n = 1
            end select
            return
        end function getNSymmetries
        
        function symmetryRelatedQuaternion(this,symmetry,s) result(that)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ 
    !*      return a quaternion rotated by one of the rotation symmetries in the group
            type(Quaternion),intent(in)         ::      this 
            integer,intent(in)                  ::      symmetry
            integer,intent(in)                  ::      s
            type(Quaternion)                    ::      that 
            
            type(Quaternion),dimension(:),pointer   ::      qsym
            select case(symmetry)                
                case(QUATERNION_SYMMETRY_CUBIC)
                    qsym => cubic_rot_q           
                case(QUATERNION_SYMMETRY_HEX)
                    qsym => hex_rot_q           
                case default
                    that = this
                    return
            end select
            that = this * qsym(s)
            return
        end function symmetryRelatedQuaternion
        
        
        
    end module Lib_Quaternions 
!     
!     
! !   gfortran -ffree-line-length-256  src/Lib_Quaternions.f90 src/Lib_RotationMatrices.f90 -llapack
!    
!     program testLib_Quaternions
! !---^^^^^^^^^^^^^^^^^^^^^^^^^^^
!         use Lib_Quaternions
!         use Lib_RotationMatrices
!         use iso_fortran_env
!         implicit none
!         
!         real(kind=real64),parameter     ::      PI = 3.14159265358979d0
!         real(kind=real64)       ::      theta
!         integer         ::      ii
!         real(kind=real64),dimension(3,3)    ::      RR
!         type(Quaternion)                    ::      qq
!         
!         qq = Quaternion_ctor( 1.0d0,0.0d0,0.0d0,0.0d0 )
!             call report(qq)
!             call report(inverse(qq))
!             call report(quaternionToRotMat(qq))
!         
!         
!         do ii = -2,3
!             theta = 2*PI*ii/6           !   -120,-60,0,60,120,180 
!             
!             print *,"theta ",theta
!             
!             RR(1,1:3) = (/ cos(theta), sin(theta),0.0d0 /)
!             RR(2,1:3) = (/-sin(theta), cos(theta),0.0d0 /)
!             RR(3,1:3) = (/ 0.0d0,      0.0d0,    1.0d0 /)
!             
!             qq = Quaternion_ctor(RR)
!             call report(qq)
!             call report(quaternionToRotMat(qq))
! !                          
!         end do
!         
!         do ii = -2,3
!             theta = 2*PI*ii/6           !   -120,-60,0,60,120,180 
!             
!             print *,"theta ",theta
!             
!             RR(1,1:3) = (/-cos(theta), sin(theta),0.0d0 /)
!             RR(2,1:3) = (/ sin(theta), cos(theta),0.0d0 /)
!             RR(3,1:3) = (/ 0.0d0,      0.0d0,   -1.0d0 /)
!             
!             qq = Quaternion_ctor(RR)
!             call report(qq)
!             call report(quaternionToRotMat(qq))
! !                          
!         end do
!         
!         do ii = -2,3
!             theta = 2*PI*ii/6           !   -120,-60,0,60,120,180 
!             
!             print *,"theta ",theta
!             
!             RR(1,1:3) = (/ cos(theta),-sin(theta),0.0d0 /)
!             RR(2,1:3) = (/-sin(theta),-cos(theta),0.0d0 /)
!             RR(3,1:3) = (/ 0.0d0,      0.0d0,   -1.0d0 /)
!             
!             qq = Quaternion_ctor(RR)
!             call report(qq)
!             call report(quaternionToRotMat(qq))
! !                          
!         end do
!         print *,""
!         print *,"done"
!         print *,""
!      
!             
!                         
!     
!     end program testLib_Quaternions
!     