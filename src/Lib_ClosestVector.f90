
    module Lib_ClosestVector
!---^^^^^^^^^^^^^^^^^^^^^^^^
!*      given a set of M target vectors v_k ( M < 100 )
!*      approximately evenly spaced with |v_k| < R
!*      find the closest target vector ( in terms of Euclidiean distance )
!*      to a test vector w with |w| < R

        
        
        use iso_fortran_env
        implicit none
        private
    
    !---    
            
        public      ::      ClosestVector_ctor
        public      ::      report
        public      ::      delete    
        public      ::      clone           !   deep copy with allocation
    
        public      ::      getNsets,getNtargetVectors
        public      ::      targetVector
        public      ::      zeroCentreOfMass
        public      ::      closestTargetVector
        public      ::      closestTargetVectorSet
        public      ::      distanceTargetVectorSet
        
        logical,public      ::      LIB_CV_DBG = .false.
        
    !---
    
        type,public     ::      ClosestVector
            private
            integer                                         ::      m           !   number of target vectors per set
            integer                                         ::      n           !   number of target vector sets
            real(kind=real64),dimension(:,:,:),pointer      ::      v           !   target vector (3,m,n)
            
            real(kind=real64)                     ::      rMax            !   an indicative maximum lengthscale for target and test vectors
            real(kind=real64)                     ::      idelta          !   fine cubic grid subdivision inverse length
            integer                               ::      nDiv            !   half-number of cubic grid mesh points per side ( ix = -nDiv:nDiv )
            integer,dimension(:,:,:,:,:),pointer  ::      indx            !   (0:12,n,-nDiv:nDiv,..)  0 = number of neighbours to test for this mesh point
        end type
        
        
    !---
    
        interface       ClosestVector_ctor
            module procedure        ClosestVector_null
            module procedure        ClosestVector_ctor0  
            module procedure        ClosestVector_ctor0a  
            module procedure        ClosestVector_ctor1 
            module procedure        ClosestVector_ctor1a        
        end interface
                
        interface       report
            module procedure        report0           
        end interface
                
        interface       delete
            module procedure        delete0           
        end interface
        
        
        interface       clone
            module procedure        clone0
        end interface
        
        interface       closestTargetVector
            !module procedure        closestTargetVector_naive     
            module procedure        closestTargetVector_mesh
            module procedure        closestTargetVector_meshall
        end interface
         
    contains
!---^^^^^^^^



        function ClosestVector_null() result(this)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !       null constructor
            type(ClosestVector)     ::      this
            this%m = 0
            this%n = 0
            nullify(this%v)
            this%nDiv = 0
            this%rMax = 0.0d0
            this%idelta = 0.0d0
            nullify(this%indx)
            return
        end function ClosestVector_null

        function ClosestVector_ctor0(v) result(this)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      standard constructor
            real(kind=real64),dimension(:,:),intent(in)     ::      v
            type(ClosestVector)     ::      this
            this = ClosestVector_null()
            this%m = size(v,dim=2)
            this%n = 1
            allocate(this%v(3,this%m,1))
            this%v(1:3,1:this%m,1) = v(1:3,1:this%m)
            return
        end function ClosestVector_ctor0
        
        function ClosestVector_ctor0a(v) result(this)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      standard constructor
            real(kind=real64),dimension(:,:,:),intent(in)     ::      v
            type(ClosestVector)     ::      this
            this = ClosestVector_null()
            this%m = size(v,dim=2)
            this%n = size(v,dim=3)
            allocate(this%v(3,this%m,this%n))
            this%v(1:3,1:this%m,1:this%n) = v(1:3,1:this%m,1:this%n)
            return
        end function ClosestVector_ctor0a

        function ClosestVector_ctor1(v,nDiv) result(this)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      complex constructor - allocate a fine cubic mesh, and determine which target vector is closest to each point
            real(kind=real64),dimension(:,:),intent(in)     ::      v
            integer,intent(in)                              ::      nDiv        !   number of fine grid points per side.
            type(ClosestVector)     ::      this
            integer                 ::      ii,nn,jj,kk                      
            real(kind=real64)       ::      dd,delta,dmin
            integer                 ::      ix,iy,iz
            real(kind=real64),dimension(3)  ::      xx
            
            this = ClosestVector_ctor0(v)
            
        !---    find longest target vector
            this%rMax = 0.0d0
            do ii = 1,this%m                
                this%rMax = max( this%rMax,norm2(this%v(:,ii,1)) )
            end do
            print *,"Lib_ClosestVector::ClosestVector_ctor1 info - longest target vector ",this%rMax
            delta = this%rMax/nDiv
            
            
        !---    now find typical spacing of target vectors    
            dd = (4*3.141592654d0/3)*this%rMax**3           !   volume 
            dd = this%m/dd                                  !   number density
            dd = dd**(1/3.0d0)                              !   length scale
            this%rMax = this%rMax + dd                      !   add buffer space.
            print *,"Lib_ClosestVector::ClosestVector_ctor1 info - long range inc buffer ",this%rMax
                        
        !---    find the lengthscale of the fine mesh and allocate
            nn = ceiling(this%rMax/delta)                   !   new larger number of divisions includes buffer
            delta = this%rMax/nn                            !   final spacing of grid                                
            this%nDiv = nn
            allocate(this%indx(0:12,1,-this%nDiv:this%nDiv,-this%nDiv:this%nDiv,-this%nDiv:this%nDiv))
            this%idelta = 1/delta
            print *,"Lib_ClosestVector::ClosestVector_ctor1 info - fine grid nDiv,space  ",this%nDiv,delta
            
        !---    run through each point on the fine mesh, find distance to each target. 
        !       If it is within delta of the optimum distance, need to store as a possible.
            this%indx = 0
            do iz = -this%nDiv,this%nDiv
                do iy = -this%nDiv,this%nDiv
                    do ix = -this%nDiv,this%nDiv
                        xx = (/ix,iy,iz/)*delta                             !   position of the mesh point
                        call closestTargetVector_naive(this,xx,1,jj,dmin)   !   find the nearest target
                        dmin = sqrt(dmin)
                        
                        do kk = 1,this%m
                            dd = norm2( this%v(:,kk,1)-xx )               !   distance to this target
                            if (dd<dmin + sqrt(3.0d0)*delta) then       !   this target is close to the optimum distance
                                jj = this%indx(0,1,ix,iy,iz) + 1
                                this%indx(0,1,ix,iy,iz) = jj
                                this%indx(jj,1,ix,iy,iz) = kk
                            end if                         
                        end do
                    end do
                end do
            end do                                        
            
            !do ii = 0,12
            !    print *,"number of fine mesh with neigh count ",ii," = ",count( this%indx(0,1,:,:,:) == ii )
            !end do
            
            
            
            return
        end function ClosestVector_ctor1

        function ClosestVector_ctor1a(v,nDiv) result(this)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      complex constructor - allocate a fine cubic mesh, and determine which target vector is closest to each point
            real(kind=real64),dimension(:,:,:),intent(in)   ::      v
            integer,intent(in)                              ::      nDiv        !   number of fine grid points per side.
            type(ClosestVector)     ::      this
            integer                 ::      ii,nn,jj,kk                      
            real(kind=real64)       ::      dd,delta,dmin
            integer                 ::      ix,iy,iz
            real(kind=real64),dimension(3)  ::      xx
            
            this = ClosestVector_ctor0a(v)
            
        !---    find longest target vector
            this%rMax = 0.0d0
            do jj = 1,this%n
                do ii = 1,this%m                
                    this%rMax = max( this%rMax,norm2(this%v(:,ii,jj)) )
                end do
            end do
            !print *,"Lib_ClosestVector::ClosestVector_ctor1a info - longest target vector ",this%rMax
            delta = this%rMax/nDiv
            
            
        !---    now find typical spacing of target vectors    
            dd = (4*3.141592654d0/3)*this%rMax**3           !   volume 
            dd = this%m/dd                                  !   number density
            dd = dd**(1/3.0d0)                              !   length scale
            this%rMax = this%rMax + dd                      !   add buffer space.
            !print *,"Lib_ClosestVector::ClosestVector_ctor1a info - long range inc buffer ",this%rMax
                        
        !---    find the lengthscale of the fine mesh and allocate
            nn = ceiling(this%rMax/delta)                   !   new larger number of divisions includes buffer
            delta = this%rMax/nn                            !   final spacing of grid                                
            this%nDiv = nn
            allocate(this%indx(0:12,this%n,-this%nDiv:this%nDiv,-this%nDiv:this%nDiv,-this%nDiv:this%nDiv))
            this%idelta = 1/delta
            !print *,"Lib_ClosestVector::ClosestVector_ctor1a info - fine grid nDiv,space  ",this%nDiv,delta
            
        !---    run through each point on the fine mesh, find distance to each target. 
        !       If it is within delta of the optimum distance, need to store as a possible.
            this%indx = 0
            do iz = -this%nDiv,this%nDiv
                do iy = -this%nDiv,this%nDiv
                    do ix = -this%nDiv,this%nDiv
                        xx = (/ix,iy,iz/)*delta                             !   position of the mesh point
                        
                        
                        do nn = 1,this%n
                        
                            call closestTargetVector_naive(this,xx,nn,jj,dmin)   !   find the nearest target
                            dmin = sqrt(dmin)
                            
                            do kk = 1,this%m
                                dd = norm2( this%v(:,kk,nn)-xx )               !   distance to this target
                                if (dd<dmin + sqrt(3.0d0)*delta) then          !   this target is close to the optimum distance
                                    jj = this%indx(0,nn,ix,iy,iz) + 1
                                    this%indx(0,nn,ix,iy,iz) = jj
                                    this%indx(jj,nn,ix,iy,iz) = kk
                                end if                         
                            end do
                            
                        end do
                    end do
                end do
            end do                                        
            
           ! do ii = 0,12
           !     print *,"number of fine mesh with neigh count ",ii," = ",count( this%indx(0,:,:,:,:) == ii )/real(this%n)
           ! end do
            
            
            
            return
        end function ClosestVector_ctor1a


        
        subroutine report0(this,u,o)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      simple screen dump to unit u[=screen], offset o[=0]
            type(ClosestVector),intent(in)      ::      this
            integer,intent(in),optional         ::      u,o
            integer         ::      uu,oo
            uu = 6 ; if (present(u)) uu = u
            oo = 0 ; if (present(o)) oo = o    
            if (this%ndiv==0) then        
                write (unit=uu,fmt='(3(a,i6))') repeat(" ",oo)//"ClosestVector[m=",this%m,",n=",this%n,"]"
            else
                write (unit=uu,fmt='(a,i6,a,i6,a,f12.5,a)') repeat(" ",oo)//"ClosestVector[m=",this%m,",n=",this%n,",rMax=",this%rMax,"]"
            end if
            return
        end subroutine report0
        

        
        subroutine delete0(this)
    !---^^^^^^^^^^^^^^^^^^^^^^^^
    !*      simple destructor
            type(ClosestVector),intent(inout)   ::      this
            if (this%n==0) return
            deallocate(this%v)
            deallocate(this%indx)
            this = ClosestVector_ctor()
            return
        end subroutine delete0
        
        subroutine clone0(this,that)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      deep copy with allocation this = that
            type(ClosestVector),intent(out)   ::      this
            type(ClosestVector),intent(in)   ::      that
            this = ClosestVector_ctor()
            if (that%n==0) return
             
            this%m = that%m
            this%n = that%n
            this%rMax = that%rMax
            this%iDelta = that%iDelta
            this%nDiv = that%nDiv
            
            allocate(this%v(3,this%m,this%n))
            this%v(:,:,:) = that%v(:,:,:)
            
            allocate(this%indx(0:12,this%n,-this%nDiv:this%nDiv,-this%nDiv:this%nDiv,-this%nDiv:this%nDiv))
            this%indx(:,:,:,:,:) = that%indx(:,:,:,:,:)
            
            return
        end subroutine clone0
        
    !---
         
        
    
    
        pure function getNsets(this) result(n)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      return number of target sets
            type(ClosestVector),intent(in)      ::      this
            integer                             ::      n
            n = this%n
            return
        end function getNsets
        
    
        pure function getNtargetVectors(this) result(m)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      return the number of target vectors
            type(ClosestVector),intent(in)      ::      this
            integer                             ::      m
            m = this%m
            return
        end function getNtargetVectors
        
    !---
        
    
        pure function targetVector(this,k,n) result(v)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      return the kth target vector in the nth set
            type(ClosestVector),intent(in)      ::      this
            integer,intent(in)                  ::      k,n
            real(kind=real64),dimension(3)      ::      v
            v = this%v(:,k,n)
            return
        end function targetVector
        
        
        pure subroutine closestTargetVector_naive(this,w,n,k,d2)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      find the closest target vector in set n to test vector w
    !*      return its distance squared d2
            type(ClosestVector),intent(in)              ::      this
            integer,intent(in)                          ::      n
            real(kind=real64),dimension(3),intent(in)   ::      w
            integer,intent(out)                         ::      k
            real(kind=real64),intent(out)               ::      d2
            
            integer             ::      ii
            real(kind=real64)   ::      d2min,dx 
          
            d2min = huge(1.0)
            do ii = 1,this%m
                dx = this%v(1,ii,n) - w(1) ; d2 = dx*dx
                dx = this%v(2,ii,n) - w(2) ; d2 = d2 + dx*dx
                dx = this%v(3,ii,n) - w(3) ; d2 = d2 + dx*dx
                
                if (d2<d2min) then
                    k = ii
                    d2min = d2
                end if
            end do
              
            
!            if (LIB_CV_DBG) then
!                do ii = 1,this%m
!                    dx = this%v(1,ii) - w(1) ; d2 = dx*dx
!                    dx = this%v(2,ii) - w(2) ; d2 = d2 + dx*dx
!                    dx = this%v(3,ii) - w(3) ; d2 = d2 + dx*dx
!                
!                    if (ii==k) then
!                        print *,"closestTargetVector_naive dbg * ",ii,d2," *"
!                    else                   
!                        print *,"closestTargetVector_naive dbg   ",ii,d2
!                    end if
!                end do
!            end if
            d2 = d2min    
            return
        end subroutine closestTargetVector_naive
        
        pure subroutine closestTargetVector_mesh(this,w,n,k,d2 )!,mesh)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      find the closest target vector in set n to test vector w
    !*      return its distance squared d2
            type(ClosestVector),intent(in)              ::      this
            real(kind=real64),dimension(3),intent(in)   ::      w
            integer,intent(in)                          ::      n
            integer,intent(out)                         ::      k
            real(kind=real64),intent(out)               ::      d2
            !logical,intent(in)              ::      mesh                
            integer             ::      ix,iy,iz,ii,jj 
            real(kind=real64)   ::      d2min,dx  
             
             
        !---    check if this test vector is within the acceptable range
            d2 = w(1)*w(1) + w(2)*w(2) + w(3)*w(3)
            if (d2 > this%rMax*this%rMax) then
                call closestTargetVector_naive(this,w,n,k,d2)
                return
            end if
            
        !---    find the close mesh point            
            ix = nint( this%idelta*w(1) ) 
            iy = nint( this%idelta*w(2) ) 
            iz = nint( this%idelta*w(3) ) 
                        
        !---    test possible neighbours             
            k = this%indx(1,n,ix,iy,iz)
            dx = this%v(1,k,n) - w(1) ; d2min = dx*dx
            dx = this%v(2,k,n) - w(2) ; d2min = d2min + dx*dx
            dx = this%v(3,k,n) - w(3) ; d2min = d2min + dx*dx    
            do jj = 2,this%indx(0,n,ix,iy,iz)
                ii = this%indx(jj,n,ix,iy,iz)
                dx = this%v(1,ii,n) - w(1) ; d2 = dx*dx
                dx = this%v(2,ii,n) - w(2) ; d2 = d2 + dx*dx
                dx = this%v(3,ii,n) - w(3) ; d2 = d2 + dx*dx
                if (d2<d2min) then
                    k = ii
                    d2min = d2
                end if
            end do
              
            
            
!            if (LIB_CV_DBG) then
!                print *,"closestTargetVector_mesh dbg    ",ix,iy,iz
!                print *,"closestTargetVector_mesh dbg    ",(/ix,iy,iz/)/this%idelta
!                print *,"closestTargetVector_mesh dbg    ",w
!                print *,"closestTargetVector_mesh dbg    ",w-(/ix,iy,iz/)/this%idelta
!                
!                do jj = 1,this%indx(0,ix,iy,iz)
!                    ii = this%indx(jj,ix,iy,iz)
!                    dx = this%v(1,ii) - w(1) ; d2 = dx*dx
!                    dx = this%v(2,ii) - w(2) ; d2 = d2 + dx*dx
!                    dx = this%v(3,ii) - w(3) ; d2 = d2 + dx*dx
!                    if (ii==k) then
!                        print *,"closestTargetVector_mesh  dbg * ",ii,d2," *"
!                    else                   
!                        print *,"closestTargetVector_mesh  dbg   ",ii,d2
!                    end if
!                end do
!            end if
            
            
            d2 = d2min    
            return
        end subroutine closestTargetVector_mesh
        
        pure subroutine closestTargetVector_meshall(this,w,k,d2 )
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      find the closest target vector in all sets to test vector w
    !*      return its distance squared d2
            type(ClosestVector),intent(in)              ::      this
            real(kind=real64),dimension(3),intent(in)   ::      w
            integer,dimension(:),intent(out)            ::      k
            real(kind=real64),dimension(:),intent(out)  ::      d2
            !logical,intent(in)              ::      mesh                
            integer             ::      ix,iy,iz,ii,jj,nn
            real(kind=real64)   ::      dx,dd,d2min
            integer,dimension(0:12) ::  indx_tmp
             
        !---    check if this test vector is within the acceptable range
            dd = w(1)*w(1) + w(2)*w(2) + w(3)*w(3)
            if (dd > this%rMax*this%rMax) then
                do nn = 1,this%n
                    call closestTargetVector_naive(this,w,nn,k(nn),d2(nn))
                end do
                return
            end if
            
        !---    find the close mesh point            
            ix = nint( this%idelta*w(1) ) 
            iy = nint( this%idelta*w(2) ) 
            iz = nint( this%idelta*w(3) ) 
                        
        !---    test possible neighbours             
            do nn = 1,this%n
                indx_tmp(0:12) = this%indx(0:12,nn,ix,iy,iz)
                ii = indx_tmp(1)
                k(nn) = ii
                dx = this%v(1,ii,nn) - w(1) ; d2min = dx*dx
                dx = this%v(2,ii,nn) - w(2) ; d2min = d2min + dx*dx
                dx = this%v(3,ii,nn) - w(3) ; d2min = d2min + dx*dx    
                do jj = 2,indx_tmp(0)
                    ii = indx_tmp(jj)
                    dx = this%v(1,ii,nn) - w(1) ; dd = dx*dx
                    dx = this%v(2,ii,nn) - w(2) ; dd = dd + dx*dx
                    dx = this%v(3,ii,nn) - w(3) ; dd = dd + dx*dx
                    if (dd<d2min) then
                        k(nn) = ii
                        d2min = dd
                    end if
                end do                  
                d2(nn) = d2min    
            end do
            return
        end subroutine closestTargetVector_meshall
        
        
        subroutine closestTargetVectorSet(this,w,rMin,k,dd,np )
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      find the closest target vector set k to input test vectors w
    !*      ignore any vectors in w further than rMin from a neighbour
    !*      also return distance dd and number of pairs np
    
            type(ClosestVector),intent(in)                   ::      this
            real(kind=real64),dimension(:,:),intent(in)     ::      w
            real(kind=real64),intent(in)                    ::      rMin
            integer,intent(out)                             ::      k
            real(kind=real64),intent(out)                   ::      dd
            integer,intent(out)                             ::      np
                    
            real(kind=real64),dimension(this%n) ::      d2          !   distance of set w to target set.
            integer                 ::      nn 
            
            integer,dimension(1:this%n) ::  nPair  
            
             

            
            do nn = 1,this%n
                call distanceTargetVectorSet(this,w,nn,rMin,d2(nn),nPair(nn) ) 
            end do
 
            k = minloc(d2,dim=1)        !   which set of targets has the smallest squared projection?
            
            dd = d2(k)
            np = nPair(k)
           
            return
        end subroutine closestTargetVectorSet
        
        
        
        subroutine distanceTargetVectorSet(this,w,n,rMin,d2,nPair )
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      find the distance to a target vector set n to test vectors w
    !*      ignore any vectors in w further than rMin from a neighbour
            type(ClosestVector),intent(in)                  ::      this
            real(kind=real64),dimension(:,:),intent(in)     ::      w
            real(kind=real64),intent(in)                    ::      rMin
            integer,intent(in)                              ::      n             
            real(kind=real64),intent(out)                   ::      d2          !   distance of set w to target set.
            integer,intent(out)                             ::      nPair
            
            integer                 ::      ii,jj,kk!,nn
            
            real(kind=real64)       ::      dd,d2min 
            integer                 ::      ix,iy,iz        ! , kkkk
            integer,dimension(0:12) ::      indx_tmp
            
            real(kind=real32)       ::       w1,w2,w3  , dx
            real(kind=real32),dimension(3,1:this%m) ::      vv
             
            d2 = 0.0d0 ; nPair = 0 
            vv(1:3,1:this%m) = real(this%v(1:3,1:this%m,n),kind=real32)           
            
            do ii = 1,size(w,dim=2)
            
                w1 = real(w(1,ii),kind=real32)
                w2 = real(w(2,ii),kind=real32)
                w3 = real(w(3,ii),kind=real32)
            
            !---    check if this test vector is in bounds.                
                dd = w1*w1 + w2*w2 + w3*w3
                if (dd > this%rMax*this%rMax) cycle     !   out of bounds- all must take same penalty ie can ignore test vector.

            !---    find location of test vector                            
                ix = nint( this%idelta*w1 ) 
                iy = nint( this%idelta*w2 ) 
                iz = nint( this%idelta*w3 )                 
                
            !---    compare to possible targets - find minimum separation
                indx_tmp(0:12) = this%indx(0:12,n,ix,iy,iz)
                d2min = rMin*rMin                   !   don't count atoms more than rMin away
                do jj = 1,indx_tmp(0)
                    kk = indx_tmp(jj)
                    dx = vv(1,kk) - w1 ; dd = dx*dx
                    dx = vv(2,kk) - w2 ; dd = dd + dx*dx
                    dx = vv(3,kk) - w3 ; dd = dd + dx*dx                    
                    d2min = min(d2min,dd)                    
                end do    
                                 
                d2 = d2 + d2min
                if (d2min<rMin*rMin-1.0d-8) nPair = nPair + 1                       
            end do
                        
            return
        end subroutine distanceTargetVectorSet
        
        
        
        
    !---    helper methods
    
        pure subroutine zeroCentreOfMass( v )
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      remove any linear offset
            real(kind=real64),dimension(:,:),intent(inout)      ::      v
            
            real(kind=real64)       ::      vx 
            integer                 ::      mm
            
            mm = size(v,dim=2)
            if (mm == 0) then
                return
            else if (mm == 1) then
                v = 0
                return
            else
                vx = sum(v(1,:))/mm ; v(1,:) = v(1,:) - vx 
                vx = sum(v(2,:))/mm ; v(2,:) = v(2,:) - vx 
                vx = sum(v(3,:))/mm ; v(3,:) = v(3,:) - vx 
            end if
            return
        end subroutine zeroCentreOfMass                            
             
        
        
        

    end module Lib_ClosestVector      
        