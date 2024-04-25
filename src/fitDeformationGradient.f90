
     
    module fitDeformationGradient 
!---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
!*
!*      given an expected set of atom positions r_k
!*      and an observed set of atom positions r'_k
!*      find the optimal transformation
!*          r'_k = T r_k + d
!*      
        use Lib_DeformationGradients
        use Lib_ClosestVector
        use Lib_RotationMatrices
        use Lib_Quaternions
        use Lib_Lattices
        use iso_fortran_env
        use Lib_RandomSeed
        
 
        implicit none
        private
        
        include "cubic0072.h"
        include "hex0064.h"
        include "sym_breaker.h"
        
        external    ::      DSYSV
       
    !---
        
        public      ::      FitDG_ctor
        public      ::      report
        public      ::      delete
        public      ::      clone           !   deep copy with allocation
        
        
        public      ::      fit
        public      ::      getRotationMatrix
        public      ::      getnRotationMatrices
        
        
       ! public      ::      defGradDistance
        
       logical,public       ::      FitDG_DBG = .false.
       
       
    !---
        
        type,public          ::    FitDG
            private
            type(Lattice)                               ::      latt            
            integer                                     ::      N       !   number of orientation trials in target set
            real(kind=real64),dimension(:,:,:),pointer  ::      U       !   orientation trials (3,3,N)
            type(ClosestVector),dimension(:),pointer    ::      cv      !   storage of target vectors
        end type
        
    !---

        interface       FitDG_ctor
            module procedure        FitDG_null
            module procedure        FitDG_ctor0
            module procedure        FitDG_ctor1
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
        
        
        interface       fit
            module procedure        fit0        
        end interface
        
        
        
            
                
        
                
        
    contains
!---^^^^^^^^

        function FitDG_null() result(this)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
            type(FitDG)     ::      this
            this%N = 0
            this%latt = Lattice_ctor()
            nullify(this%U)
            nullify(this%cv)
            return
        end function FitDG_null

        function FitDG_ctor0(latt) result(this)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ 
    !*      read the default rotation matrices 
            type(Lattice),intent(in)                        ::      latt
           

            type(FitDG)     ::      this
             
            !real(kind=real64),dimension(3,3)        ::      UU
             
            integer             ::      nn,ii,kk,mm 
            real(kind=real64),dimension(:,:,:),allocatable  ::      rr
            real(kind=real64),dimension(:,:,:),pointer      ::      uu_tmp
            type(Quaternion)    ::      qq
            integer,parameter   ::      NDIV = 10
            
            this%latt = latt
            
            
            select case(getSymmetry(this%latt))
                case(QUATERNION_SYMMETRY_HEX)
                    uu_tmp => u_hex0064                
                    !print *,"fitDeformationGradient::FitDG_ctor0 info - using angle dataset ""u_hex0064"""  
                case default
                    uu_tmp => u_cubic0072           
                    !print *,"fitDeformationGradient::FitDG_ctor0 info - using angle dataset ""u_cubic0072"""
            end select
            
            
            
            this%N = size(uu_tmp,dim=3)
            allocate(this%U(3,3,this%N))
            do ii = 1,this%N                    
                qq = Quaternion_ctor( uu_tmp(1:3,1:3,ii) )      
                call imposeSymmetry( qq,getSymmetry(this%latt) )
                this%U(:,:,ii) = QuaternionToRotMat(qq)
            end do
            
            
            nn = getNMotif(this%latt)
            allocate(this%cv(nn))
            
            
            allocate(rr(3,getNneighbours(this%latt),this%N))
            do nn = 1,getNMotif(this%latt)
                mm = getNneighbours(this%latt,nn)
                do ii = 1,this%N
                    do kk = 1,mm
                        rr(:,kk,ii) = rotateVector( this%U(:,:,ii),getNeighbour(this%latt,kk,nn) )
                    end do
                end do
                this%cv(nn) = ClosestVector_ctor(rr(:,1:mm,1:this%N),NDIV)
            end do
            deallocate(rr)
            
             
             
            return
        end function FitDG_ctor0


        function FitDG_ctor1(latt,filename) result(this)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      read the rotation matrices in from a given file
            type(Lattice),intent(in)                        ::      latt
            
            character(len=*),intent(in)                     ::      filename

            type(FitDG)     ::      this
             
            real(kind=real64),dimension(3,3)        ::      UU
            type(Quaternion)                        ::      qq
            integer             ::      nn,ii,kk,mm 
            real(kind=real64),dimension(:,:,:),allocatable  ::      rr
           
            integer,parameter   ::      NDIV = 10
            
            this%latt = latt
            open(unit=991,file=trim(filename),action="read")            
                read(unit=991,fmt=*)        this%N
                allocate(this%U(3,3,this%N))
                do ii = 1,this%N
                    read(unit=991,fmt=*)  UU
                    qq = Quaternion_ctor(UU)
                    call imposeSymmetry( qq,getSymmetry(this%latt) )
                    this%U(:,:,ii) = quaternionToRotMat(qq)
                end do
            close(unit=991)
            
            
            nn = getNMotif(this%latt)
            allocate(this%cv(nn))
            
            
            allocate(rr(3,getNneighbours(this%latt),this%N))
            do nn = 1,getNMotif(this%latt)
                mm = getNneighbours(this%latt,nn)
                do ii = 1,this%N
                    do kk = 1,mm
                        rr(:,kk,ii) = rotateVector( this%U(:,:,ii),getNeighbour(this%latt,kk,nn) )
                    end do
                end do
                this%cv(nn) = ClosestVector_ctor(rr(:,1:mm,1:this%N),NDIV)
            end do
            deallocate(rr)
            
             
             
            return
        end function FitDG_ctor1

    !---
        
        subroutine report0(this,u,o)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      
            type(FitDG),intent(in)          ::      this
            integer,intent(in),optional     ::      u,o
            
            integer             ::      uu,oo,ii
            uu = 6 ; if (present(u)) uu = u
            oo = 0 ; if (present(o)) oo = o
            
            write(unit=uu,fmt='(3(a,i6))') repeat(" ",oo)//"FitDG[n=",this%N,"]"
            call report(this%latt,uu,oo+4)
            do ii = 1,getNMotif(this%latt)
                write(unit=uu,fmt='(3(a,i6))') repeat(" ",oo)//"lattice motif ",ii
                call report(this%cv(ii),uu,oo+4)
            end do
            return
        end subroutine report0

        subroutine delete0(this)
    !---^^^^^^^^^^^^^^^^^^^^^^^^
            type(FitDG),intent(inout)       ::      this
            integer         ::      ii
            if (this%N == 0) return
            
            deallocate(this%U)
            do ii = 1,getNMotif(this%latt)
                call delete(this%cv(ii))
            end do            
            call delete(this%latt)
            this = FitDG_ctor()
            
            return
        end subroutine delete0
        
        subroutine clone0(this,that)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      deep copy with allocation this = that
            type(FitDG),intent(out)       ::      this
            type(FitDG),intent(in)        ::      that
            integer         ::      ii
            this = FitDG_ctor()
            if (that%N == 0) return
            this%N = that%N
            allocate(this%U(3,3,this%N))
            this%U(:,:,:) = that%U(:,:,:)
            this%latt = Lattice_ctor()
            call clone(this%latt,that%latt)
            allocate(this%cv(getNMotif(this%latt)))
            do ii = 1,getNMotif(this%latt)
                this%cv(ii) = ClosestVector_ctor()
                call clone(this%cv(ii),that%cv(ii))
            end do                        
            return
        end subroutine clone0
        
        
        
        
    !---
    
        subroutine fit0( this, rmin, v , T ,m, dd )
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      given an input data set v
    !*      find the best transformation  v = T v0 + delta
    !*      and return motif m
            type(FitDG),intent(in)                          ::      this
            real(kind=real64),intent(in)                    ::      rmin
            real(Kind=real64),dimension(:,:),intent(in)     ::      v            
            real(kind=real64),dimension(3,3),intent(out)    ::      T
            integer,intent(out)                             ::      m
            real(kind=real64),intent(out)                   ::      dd 
            
            real(kind=real64)                   ::      dm
            real(kind=real64),dimension(3,3)    ::      Tm
            integer                             ::      mm
            
            dd = huge(1.0)
  
           
            do mm = 1,getNmotif(this%latt)
                call fit2( this,mm, v, Tm , rmin , dm )
                if (dm<dd) then
                    dd = dm
                    T = Tm
                    m = mm
                end if
            end do
 
            return
!              
        end subroutine fit0
        
    
        subroutine fit2( this,m, v, T , rmin , dd)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      given an input data set v
    !*      find the best transformation  v = T v0 + delta
    !*      assuming starting from motif m
            type(FitDG),intent(in)                          ::      this
            integer,intent(in)                              ::      m
            real(Kind=real64),dimension(:,:),intent(in)     ::      v
            real(kind=real64),dimension(3,3),intent(out)    ::      T
            real(kind=real64),intent(in)                    ::      rmin
            real(kind=real64),intent(out)                   ::      dd 
            
            real(kind=real64),dimension(:,:),allocatable    ::      vv           !   offset removed
 
            
            integer                             ::      mm,mPaired,bestMat,mmmm , step ,jj
            integer,parameter                   ::      UNPAIRED = -1
            integer,parameter                   ::      NSTEPS = 10
            real(kind=real64)                   ::      ff,bestf

            real(kind=real64),dimension(3,3)    ::      XX,RR    
!            integer                             ::      seed
            
            
            
             
            
            
        !---    allocate memory for the working set, and zero the centre of position.
            mm = size(v,dim=2)
            allocate(vv(3,mm))
             
            vv = v
            call zeroCentreOfMass( vv )
            
            if (FitDG_DBG) then
                print *,"fitDeformationGradient::fit2 info - fitting to ",mm," vectors in target set ",m," rmin ",rmin
                do jj = 1,mm
                    write (*,fmt='(i6,f16.8,a,100f16.8)') jj,norm2(vv(:,jj))," x ",vv(:,jj)
                end do
            end if
            
            
            
        !   We want to find Um v0 ~ R v
        !   The def grad is found by minimising 
        !       S = | T v0 + delta - (Um^T R v) |^2
            
        !---    step 1: coarse fitting to the stored set of vectors              
        !   this will usually return a good set, but can fail
            bestf = huge(1.0) ; bestMat = 0 ; XX = RotationMatrix_identity 
            do step = 1,NSTEPS
            
            !---    break symmetry. Use rMin/10 as an indicative "small" scale
                if (step>1) then
                
                    !call random_number( vv )
                    
!                    seed = SYM_BREAKER( mod(mm*step,size(SYM_BREAKER)) )
!                    call ran0array(seed,vv)
!                                                     
!                    vv = (2*vv - 1)*rMin*0.1d0
!                    XX = rndRotMat(seed,0.1d0)
!                    
                    call init_random_seed( SYM_BREAKER( mod(mm*step,size(SYM_BREAKER)) ) )
                    call random_number( vv )
                    vv = (2*vv - 1)*rMin*0.1d0
                    XX = rndRotMat(0.1d0)
                    
                    vv = vv + rotateVector( XX,v )
                    call zeroCentreOfMass( vv )
                end if
                            
            !---    find closest target vector set            
                call closestTargetVectorSet(this%cv(m),vv,rMin,mmmm,dd,mPaired )
                
                
                
                
            !---    give this target vector set a score - fairly arbitrarily chosen metric                                               
                ff = (1 - mPaired*1.0/mm)**2 + (dd/mm)*(2/rMin)**2
                
                if (FitDG_DBG) print *,"fitDeformationGradient::fit2 info - step ",step," motif ",m," paired ",mPaired," set ",mmmm," score ",ff," good? ",ff<0.5
                
                 
                if (ff < bestf) then
                    bestMat = mmmm ; bestf = ff ; RR = XX
                end if
                
                if (ff < 0.25d0)  exit        !   good fit. Escape
                
            end do    
                        
             
        !---    now we know the best rotation is Um v0 ~ R v, can back rotate through Um
            if (FitDG_DBG) then
                print *,"fitDeformationGradient::fit2 returns ",ff," for set ",bestMat," in motif ",m
                write(*,fmt='(3f12.4,a,3f12.4)') RR(1:3,1) ," , ", this%U(:,1,bestMat)
                write(*,fmt='(3f12.4,a,3f12.4)') RR(1:3,2) ," , ", this%U(:,2,bestMat)
                write(*,fmt='(3f12.4,a,3f12.4)') RR(1:3,3) ," , ", this%U(:,3,bestMat)
            end if

        
        
            RR = matmul( transpose(this%U(:,:,bestMat)),RR )
            vv = rotateVector( RR,v )   
            call zeroCentreOfMass( vv )
          
            !print *,"bestMat, bestf ",bestMat,bestf
            !call report(RR)
            
        !---    find the deformation gradient for this rotation                                
            call fit1( this,m, vv, T , rmin, dd  )
            if (FitDG_DBG) then
                print *,"fitDeformationGradient::fit1 returns ",dd," for set ",bestMat," in motif ",m
                write(*,fmt='(3f12.4,a,3f12.4)') RR(1:3,1) ," , ", T(1,1:3)
                write(*,fmt='(3f12.4,a,3f12.4)') RR(1:3,2) ," , ", T(2,1:3)
                write(*,fmt='(3f12.4,a,3f12.4)') RR(1:3,3) ," , ", T(3,1:3)
                call DefGradToRotMat( T , XX )
                print *,"rotation distance ",distance( quaternion_ctor(XX) ) 
            end if
                
            
            
            
            
        !---    finally need to rotate the deformation gradient back  
            T = matmul( transpose(RR),T )
!             
!             
!             !call tidyDeformationGradient( this, T  )
          !call tidyRotationMatrix(T,norescale=.true.) 
            
          
          
          return
!              
        end subroutine fit2
        
        
        
        
        subroutine fit1( this,m, v, T , rmin, dd  )
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      given an input data set v which is known to be a good match to 
    !*      default data set v0 , ie has already been rotated by one of the stored matrices ,
    !*      find the best transformation  v = T v0 + delta
            type(FitDG),intent(in)                          ::      this
            integer,intent(in)                              ::      m
            real(Kind=real64),dimension(:,:),intent(in)     ::      v
            real(kind=real64),dimension(3,3),intent(out)    ::      T
            real(kind=real64),intent(in)                    ::      rmin
            real(kind=real64),intent(out)                   ::      dd                         
            
            real(kind=real64),dimension(3)      ::      delta,dx
            
            integer             ::      ii,mm,jj,kk,mPaired!    ,mmmm , step     ! ,bestMat
            integer,dimension(:),allocatable    ::      indx
            integer,parameter                   ::      UNPAIRED = -1
            
            real(kind=real64),dimension(12,12)      ::      AA
            real(kind=real64),dimension(12)         ::      BB
            real(kind=real64),dimension(12*64)      ::      work
            integer,dimension(12)                   ::      ipiv
            real(Kind=real64),dimension(:,:),allocatable        ::      rr
             
            
            mm = getNNeighbours(this%latt,m)            
            allocate(rr(3,mm))
            rr(1:3,1:mm) = getNeighbours(this%latt,m)
            
            mm = size(v,dim=2)          
            
        !---    step 1: match each vector in the set to the nearest in target set 1 ( identity matrix )            
            allocate(indx(mm)) ; indx = UNPAIRED
            mPaired = 0
            
            
            
            do kk = 1,mm
                call closestTargetVector(this%cv(m),v(:,kk),1,jj,dd )
                if (dd < rmin*rmin ) then
                    mPaired = mPaired + 1
                    indx(kk) = jj
                end if
            end do
            
                        
        !---    step 2: find least squares fit.
        !   S = sum_k (T r_k + d - v_k)^2
           
            AA = 0.0d0
            BB = 0.0d0      !   ordering T11,T21,T31,T12,T22,T32,T13,T23,T33,d1,d2,d3
            
                    
            AA(10,10) = mPaired
            AA(11,11) = mPaired
            AA(12,12) = mPaired                
            
            do kk = 1,mm           !   loop over matched pairs
                jj = indx(kk)
                if (jj == UNPAIRED) cycle
                 
                                
            !---    dS/dT  ( lower corner ) 
                dd = rr(1,jj)*rr(1,jj) ; AA(1,1) = AA(1,1) + dd ; AA(2,2) = AA(2,2) + dd ; AA(3,3) = AA(3,3) + dd                 
                dd = rr(2,jj)*rr(1,jj) ; AA(4,1) = AA(4,1) + dd ; AA(5,2) = AA(5,2) + dd ; AA(6,3) = AA(6,3) + dd      
                dd = rr(3,jj)*rr(1,jj) ; AA(7,1) = AA(7,1) + dd ; AA(8,2) = AA(8,2) + dd ; AA(9,3) = AA(9,3) + dd  
                dd = rr(2,jj)*rr(2,jj) ; AA(4,4) = AA(4,4) + dd ; AA(5,5) = AA(5,5) + dd ; AA(6,6) = AA(6,6) + dd 
                dd = rr(3,jj)*rr(2,jj) ; AA(7,4) = AA(7,4) + dd ; AA(8,5) = AA(8,5) + dd ; AA(9,6) = AA(9,6) + dd 
                dd = rr(3,jj)*rr(3,jj) ; AA(7,7) = AA(7,7) + dd ; AA(8,8) = AA(8,8) + dd ; AA(9,9) = AA(9,9) + dd 
                         
            !---    dS/ddelta                                       
                dd = rr(1,jj) ; BB(1:3) = BB(1:3) + v(1:3,kk)*dd ; AA(10,1) = AA(10,1) + dd ; AA(11,2) = AA(11,2) + dd ; AA(12,3) = AA(12,3) + dd                
                dd = rr(2,jj) ; BB(4:6) = BB(4:6) + v(1:3,kk)*dd ; AA(10,4) = AA(10,4) + dd ; AA(11,5) = AA(11,5) + dd ; AA(12,6) = AA(12,6) + dd  
                dd = rr(3,jj) ; BB(7:9) = BB(7:9) + v(1:3,kk)*dd ; AA(10,7) = AA(10,7) + dd ; AA(11,6) = AA(11,6) + dd ; AA(12,9) = AA(12,9) + dd  
                                
                BB(10:12) = BB(10:12) + v(1:3,kk)
                                            
            end do 
            if (FitDG_DBG) then
                print *,"LAPACK D/S SYSV:"
                do ii = 1,12
                    write (*,fmt='(12f12.6,a,f12.6)') AA(ii,:)," , ",BB(ii)
                end do
            end if
            
            call DSYSV( "L",12,1,AA,12,ipiv,BB,12,work,size(work),ii )                    
            
            !   check for plausible output
            if ( (ii == 0) .and. (maxval( (/BB(1),BB(5),BB(9)/))<1.5d0) .and. (minval( (/BB(1),BB(5),BB(9)/))>0.5d0) ) then                             
                T(1:3,1)   = BB(1:3)
                T(1:3,2)   = BB(4:6)
                T(1:3,3)   = BB(7:9)
                delta(1:3) = BB(10:12)
            else
                if (FitDG_DBG) write (*,fmt='(a,i6,13f12.6)') " D/S SYSV returns ",ii,  BB ,work(1) 
                
                T(1:3,1)   = (/ 1.0d0,0.0d0,0.0d0 /)
                T(1:3,2)   = (/ 0.0d0,1.0d0,0.0d0 /)
                T(1:3,3)   = (/ 0.0d0,0.0d0,1.0d0 /)
                delta(1:3) = (/ 0.0d0,0.0d0,0.0d0 /)
            end if                        
            
            
        !---    find distance from deformed canonical set to input.                 
            dd = 0 
            do kk = 1,mm           !   loop over matched pairs
                jj = indx(kk)
                if (jj == UNPAIRED) cycle       
                dx(1:3) = T(1:3,1)*rr(1,jj) + T(1:3,2)*rr(2,jj) + T(1:3,3)*rr(3,jj) + delta(1:3) - v(1:3,kk)
                dd = dd + dx(1)*dx(1) + dx(2)*dx(2) + dx(3)*dx(3)  
            end do 
            
            
            return           
        end subroutine fit1
        
 
        
        pure function getRotationMatrix(this,i) result(U)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      return the ith stored rotation matrix
            type(FitDG),intent(in)                          ::      this
            integer,intent(in)                              ::      i
            real(kind=real64),dimension(3,3)                ::      U
            U = this%U(:,:,i)
            return
        end function getRotationMatrix
        
        pure function getnRotationMatrices(this) result(n)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      return the number of stored rotation matrices
            type(FitDG),intent(in)                          ::      this
            integer                                         ::      n
            n = this%n
            return
        end function getnRotationMatrices
        
        
    end module fitDeformationGradient
                    

    
    
    