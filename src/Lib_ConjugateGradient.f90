
    module Lib_ConjugateGradient  
!---^^^^^^^^^^^^^^^^^^^^^^^^^^^^
!*      perform a conjugate gradient relaxation to solve the symmetric equations
!*          A x = b
!*      where A is sparse and symmetric and real

!        use OMP_LIB
        use iso_fortran_env
        implicit none
        private
        
        
        public      ::      conjgrad

        
        interface           conjgrad
            module procedure            conjgrad0
            module procedure            conjgrad1
        end interface
                
        interface           sparseMatVec
            module procedure            sparseMatVec0
            module procedure            sparseMatVec1
        end interface
        
    contains
!---^^^^^^^^


        subroutine conjgrad0( A,indx, x,b , eps )
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      On input A_ik is the kth non-zero column of row i of the sparse symmetric matrix A
    !*      indx_0i is the number of non-zero columns of row i
    !*      indx_ki is the column index of the kth non-zero column
    !*      x is the initial guess, b is the target
    !*      eps is the tolerance for the solution - want |Ax-b|^2 < N eps where N is the number of rows
    !*      on output
    !*      x is the solution, eps is the error |Ax-b|^2/N
    
            real(kind=real64),dimension(:,:),intent(in)         ::      A
            integer,dimension(0:,:),intent(in)                  ::      indx
            real(kind=real64),dimension(:),intent(in)           ::      b
            real(kind=real64),dimension(:),intent(inout)        ::      x
            real(kind=real64),intent(inout)                     ::      eps
            
            real(kind=real64),dimension(:),allocatable      ::      rr,pp,qq
            real(kind=real64)                               ::      r2,aa,bb
                   
            integer                     ::      ii,kk,NN,maxSteps
                      
            NN = size(x)
            allocate(rr(NN))
            allocate(pp(NN))
            allocate(qq(NN))            
            
            call sparseMatVec( A,indx,x, rr )
            rr(1:NN) = b(1:NN) - rr(1:NN)
            pp(1:NN) = rr(1:NN)
            maxSteps =  ceiling(sqrt(1.0d0*NN))
            
            do kk = 1,maxSteps
            
            !---    r2 = r.r
                r2 = 0.0d0 ; do ii = 1,NN ; r2 = r2 + rr(ii)*rr(ii) ; end do                
                !print *,"Lib_ConjugateGradient::conjgrad() info - step ",kk,"/",maxSteps," residual ",r2/NN,"/",eps
                if (r2 <= NN*eps) exit
                       
            !---    q = A p
                call sparseMatVec( A,indx,pp, qq )
               !do ii = 1,NN
               !    write (*,fmt='(2f16.8)') pp(ii),qq(ii)
               !end do 
                
                
                
            !---    a = p.q
                aa = 0.0d0 ; do ii = 1,NN ; aa = aa + pp(ii)*qq(ii) ; end do
                 
                
                !print *," p.q = ",aa," r.r/p.q = ",r2/aa
            !---    a = r2/a
                if (abs(aa)<1.0d-12) exit
                    !   weird result: p.Ap = 0 can only be true if p is a zero eigenmode of A, or if p=0
                    !   if this is true, the best solution is probably to try again with a different guess vector
                aa = r2/aa
                
                
            !---    x = x + a p
                x(1:NN) = x(1:NN) + aa*pp(1:NN)
               
            !---    r = r - a q
                rr(1:NN) = rr(1:NN) - aa*qq(1:NN) 
                
                
            !---    aa = r.r/r2
                bb = 0.0d0 ; do ii = 1,NN ; bb = bb + rr(ii)*rr(ii) ; end do
                !print *," r'.r' = ",bb," b = ",bb/r2
                bb = bb/r2          !   note r2>0
                 
                            
             !---   p = r + ap
                 pp(1:NN) = rr(1:NN) + bb*pp(1:NN)     
                 
            end do
            
        !---    at this point we should have the solution
            call sparseMatVec( A,indx,x, rr )
            rr(1:NN) = b(1:NN) - rr(1:NN)
            r2 = 0.0d0 ; do ii = 1,NN ; r2 = r2 + rr(ii)*rr(ii) ; end do                
            eps = r2/NN
            
            deallocate(rr)
            deallocate(pp)
            deallocate(qq)           
          
            
            return
        end subroutine conjgrad0           
            
            
        subroutine sparseMatVec0( A,indx,x, y )
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      compute sparse matrix multiply Ax = y
    !*      On input A_ik is the kth non-zero column of row i of the sparse symmetric matrix A
    !*      indx_0i is the number of non-zero columns of row i
    !*      indx_ki is the column index of the kth non-zero column
    
            real(kind=real64),dimension(:,:),intent(in)         ::      A
            integer,dimension(0:,:),intent(in)                  ::      indx
            real(kind=real64),dimension(:),intent(in)           ::      x
            real(kind=real64),dimension(:),intent(out)          ::      y
            real(kind=real64)           ::      yi
            integer                     ::      ii,jj,kk,NN       
            NN = size(x)
            
!$OMP PARALLEL PRIVATE(ii,jj,kk,yi) SHARED(x,y,A,indx,NN)
!$OMP DO
            do ii = 1,NN                
                yi = 0                              !   compute y_i = A_ij x_j
                do kk = 1,indx(0,ii)                !   for each non-zero column at row i
                    jj = indx(kk,ii)                !   kth non-zero column is indexed j
                    yi = yi + A(kk,ii)*x(jj)        !   sum A_ij x_j
                end do
                y(ii) = yi
            end do
!$OMP END DO            
!$OMP END PARALLEL
            
            return
        end subroutine sparseMatVec0


        subroutine conjgrad1( A,indx, x,b , eps , restartable)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      solve Ax = b
    !*      On input A_ik is the kth non-zero column of row i of the sparse symmetric matrix A
    !*      indx_0i is the number of non-zero columns of row i
    !*      indx_ki is the column index of the kth non-zero column
    !*      x is the initial guess, b is the target
    !*      eps is the tolerance for the solution - want |Ax-b|^2 < N eps where N is the number of rows
    !*      on output
    !*      x is the solution, eps is the error |Ax-b|^2/N
    !*      if restartable, then will exit if the residual gets larger.
    
            real(kind=real32),dimension(:,:,:,:),intent(in)     ::      A
            integer,dimension(0:,:),intent(in)                  ::      indx
            real(kind=real64),dimension(:,:),intent(in)         ::      b
            real(kind=real64),dimension(:,:),intent(inout)      ::      x
            real(kind=real64),intent(inout)                     ::      eps
            logical,intent(in),optional                         ::      restartable
            
            real(kind=real64),dimension(:,:),allocatable        ::      rr,pp,qq
            real(kind=real64)                                   ::      r2,aa,bb,oldr2
                   
            integer                     ::      kk,NN,maxSteps 
                      
            NN = size(x,dim=2)
            allocate(rr(3,NN))
            allocate(pp(3,NN))
            allocate(qq(3,NN))            
            
            call sparseMatVec( A,indx,x, rr )
            rr(1:3,1:NN) = b(1:3,1:NN) - rr(1:3,1:NN)
            pp(1:3,1:NN) = rr(1:3,1:NN)
            maxSteps =  ceiling(sqrt(1.0d0*NN))
            oldr2 = huge(1.0)
            do kk = 1,maxSteps
            
            !---    r2 = r.r
!                 r2 = 0.0d0
!                 do ii = 1,NN
!                     r2 = r2 + rr(1,ii)*rr(1,ii) + rr(2,ii)*rr(2,ii) + rr(3,ii)*rr(3,ii)
!                 end do                
                r2 = dot_product3(rr) 
                
                !print *,"Lib_ConjugateGradient::conjgrad() info - step ",kk,"/",maxSteps," residual ",r2/NN,"/",eps
                
                if (r2 <= NN*eps) exit
                
                if (present(restartable)) then
                    if (restartable) then
                        if ( (kk>2).and.(r2 > oldr2) ) exit
                            !   residual getting bigger. Escape
                        oldr2 = r2
                    end if
                end if                            
                       
            !---    q = A p
                call sparseMatVec( A,indx,pp, qq )
                !do ii = 1,NN
                !    write (*,fmt='(3f16.8,a,3f16.8)') pp(:,ii),",",qq(:,ii)
                !end do 
                
                
                
            !---    a = p.q
!                 aa = 0.0d0
!                 do ii = 1,NN
!                     aa = aa + pp(1,ii)*qq(1,ii) + pp(2,ii)*qq(2,ii) + pp(3,ii)*qq(3,ii)
!                 end do
                aa = dot_product3(pp,qq)
                
                !print *," p.q = ",aa," r.r/p.q = ",r2/aa
            !---    a = r2/a
                if (abs(aa)<1.0d-12) exit
                    !   weird result: p.Ap = 0 can only be true if p is a zero eigenmode of A, or if p=0
                    !   if this is true, the best solution is probably to try again with a different guess vector
                aa = r2/aa
                
                
            !---    x = x + a p
                x(1:3,1:NN) = x(1:3,1:NN) + aa*pp(1:3,1:NN)
               
            !---    r = r - a q
                rr(1:3,1:NN) = rr(1:3,1:NN) - aa*qq(1:3,1:NN) 
                
                
            !---    aa = r.r/r2
                bb = dot_product3(rr)
!                 bb = 0.0d0
!                 do ii = 1,NN
!                     bb = bb + rr(1,ii)*rr(1,ii) + rr(2,ii)*rr(2,ii) + rr(3,ii)*rr(3,ii)
!                 end do
                !print *," r'.r' = ",bb," b = ",bb/r2
                bb = bb/r2          !   note r2>0
                 
                            
             !---   p = r + ap
                 pp(1:3,1:NN) = rr(1:3,1:NN) + bb*pp(1:3,1:NN)     
                 
            end do
            
        !---    at this point we should have the solution
            call sparseMatVec( A,indx,x, rr )
            rr(1:3,1:NN) = b(1:3,1:NN) - rr(1:3,1:NN)
            r2 = dot_product3(rr)
!             r2 = 0.0d0
!             do ii = 1,NN
!                 r2 = r2 + rr(1,ii)*rr(1,ii) + rr(2,ii)*rr(2,ii) + rr(3,ii)*rr(3,ii)
!             end do                
            eps = r2/NN
            
            deallocate(rr)
            deallocate(pp)
            deallocate(qq)           
          
            
            return
        end subroutine conjgrad1           
            
        function dot_product3( x,y ) result(d)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^    
    !*      return d = x.y , or d = x.x
            real(kind=real64),dimension(:,:),intent(in)             ::      x
            real(kind=real64),dimension(:,:),intent(in),optional    ::      y
            real(kind=real64)                                   ::      d  
            
            integer             ::      ii,nn
            real(kind=real64)   ::      privated
            
            d = 0.0d0  
            
            if (present(y)) then
!$OMP PARALLEL PRIVATE(ii,nn,privated) SHARED(x,y,d)
                privated = 0.0d0            
                nn = size(x,dim=2)          
    !$OMP DO
                do ii = 1,nn    
                    privated = privated + x(1,ii)*y(1,ii) + x(2,ii)*y(2,ii) + x(3,ii)*y(3,ii)            
                end do
    !$OMP END DO
    !$OMP CRITICAL
                d = d + privated
    !$OMP END CRITICAL    
!$OMP END PARALLEL    
            else
!$OMP PARALLEL PRIVATE(ii,nn,privated) SHARED(x,d)
                privated = 0.0d0            
                nn = size(x,dim=2)          
    !$OMP DO
                do ii = 1,nn    
                    privated = privated + x(1,ii)*x(1,ii) + x(2,ii)*x(2,ii) + x(3,ii)*x(3,ii)            
                end do
    !$OMP END DO
    !$OMP CRITICAL
                d = d + privated
    !$OMP END CRITICAL    
!$OMP END PARALLEL    
            end if            


            return
        end function dot_product3
            
                  
        
        subroutine sparseMatVec1( A,indx,x, y )
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      compute sparse matrix multiply Ax = y
    !*      On input A_ik is the kth non-zero column of row i of the sparse symmetric matrix A
    !*      indx_0i is the number of non-zero columns of row i
    !*      indx_ki is the column index of the kth non-zero column
    
            real(kind=real32),dimension(:,:,:,:),intent(in)     ::      A
            integer,dimension(0:,:),intent(in)                  ::      indx
            real(kind=real64),dimension(:,:),intent(in)         ::      x
            real(kind=real64),dimension(:,:),intent(out)        ::      y
            real(kind=real64)           ::      y1,y2,y3
            integer                     ::      ii,jj,kk,NN       
            
!$OMP PARALLEL PRIVATE(ii,jj,kk,NN,y1,y2,y3) SHARED(A,x,y,indx)
            NN = size(x,dim=2)
!$OMP DO
            do ii = 1,NN                
                y1 = 0.0d0                          !   compute y_i = A_ij x_j
                y2 = 0.0d0
                y3 = 0.0d0
                do kk = 1,indx(0,ii)                !   for each non-zero column at row i
                    jj = indx(kk,ii)                !   kth non-zero column is indexed j
                    y1 = y1 + A(1,1,kk,ii)*x(1,jj) + A(1,2,kk,ii)*x(2,jj) + A(1,3,kk,ii)*x(3,jj)        !   sum A_ij x_j
                    y2 = y2 + A(2,1,kk,ii)*x(1,jj) + A(2,2,kk,ii)*x(2,jj) + A(2,3,kk,ii)*x(3,jj)        !   sum A_ij x_j
                    y3 = y3 + A(3,1,kk,ii)*x(1,jj) + A(3,2,kk,ii)*x(2,jj) + A(3,3,kk,ii)*x(3,jj)        !   sum A_ij x_j
                end do
                y(1,ii) = y1
                y(2,ii) = y2
                y(3,ii) = y3
            end do
!$OMP END DO            
!$OMP END PARALLEL
            
            return
        end subroutine sparseMatVec1

    end module Lib_ConjugateGradient
    
    
!!   gfortran Lib_ConjugateGradient.F90 ${MYF90LIB}/Lib_RandomSeed.f90 ${MYF90LIB}/Lib_Callipers.f90 -fopenmp -o testLib_ConjugateGradient.exe -llapack   
!!   export OMP_NUM_THREADS=2                                      
!!   ./testLib_ConjugateGradient.exe                                          
!    
!    program testLib_ConjugateGradient 
!!---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
!        use Lib_ConjugateGradient 
!        use Lib_RandomSeed
!        use Lib_Callipers
!        use iso_fortran_env
!        implicit none          
!                
!        integer,parameter           ::      N = 200
!        integer,parameter           ::      L = 10
!        
!        
!        real(kind=real32),dimension(N,N)    ::      A,A0
!        real(kind=real32),dimension(N)      ::      b,b0,r
!        real(kind=real32),dimension(N*66)   ::      work
!        integer,dimension(N)                ::      ipiv
!        
!        integer     ::      ii,kk,jj,mi,mj
!        real(kind=real64)       ::      dd
!       
!        integer,dimension(0:L,N)            ::      indx
!        real(kind=real64),dimension(L,N)    ::      AS
!        logical     ::      ok
!        real(kind=real64),dimension(N)      ::      bS,xs
!        real(kind=real64)                   ::      eps
!        
!        type(Callipers)     ::      t1,t2
!        
!        
!        call init_random_seed(12345)
!        A = 0.0d0
!         
!        do ii = 1,N
!            indx(0,ii) = 1
!            indx(1,ii) = ii
!        end do
!        
!        do ii = 1,N
!            
!            do kk = 2,L 
!                do 
!                    call random_number(dd)
!                    jj = floor(1+dd*N) 
!                    !print *,ii,kk,jj,indx(0,ii),indx(0,jj)       
!                    if (jj==ii) cycle    
!                    mi = indx(0,ii)  
!                    if (any(indx(1:mi,ii)==jj)) cycle       !   already got entry                      
!                    mj = indx(0,jj)           
!                    if ( (L-mi)*(L-mj) > 0) then
!                        !   free slot for Aij
!                        call random_number(dd)
!                        A(ii,jj) = dd
!                        A(jj,ii) = dd
!                        
!                        mi = mi + 1
!                        indx(0,ii) = mi
!                        indx(mi,ii) = jj
!                        AS(mi,ii) = dd
!                        
!                        mj = mj + 1
!                        indx(0,jj) = mj
!                        indx(mj,jj) = ii
!                        AS(mj,jj) = dd
!                        
!                        exit
!                    end if
!                    if (L==mi) exit
!                end do
!             end do
!        end do
!        
!        do ii = 1,N
!            dd = sum(A(:,ii))
!            dd = 1 - dd
!            A(ii,ii) = dd
!            AS(1,ii) = dd            
!        end do
!        
!        print *,"dense matrix "
!        do ii = 1,N
!            write (*,fmt='(1000f10.5)') A(ii,:)
!        end do
!        
!        print *,""
!        print *,"sparse matrix"
!        do ii = 1,N
!            do jj = 1,N
!                ok = .false.
!                do kk = 1,indx(0,ii)
!                    if (indx(kk,ii)==jj) then
!                        write (*,fmt='(f10.5)',advance="no") AS(kk,ii)
!                        ok = .true.
!                        exit
!                    end if
!                end do
!                if (.not. ok) write (*,fmt='(f10.5)',advance="no") 0.0d0
!            end do
!            write (*,fmt='(a)',advance="yes") ""
!        end do
!        
!        call random_number(b)
!        b0 = b
!        A0 = A
!        
!        t1 = Callipers_ctor()
!        call SSYSV( "U",N,1,A,N,IPIV,b,N,work,size(work),ii)
!        call pause(t1)
!        print *,"SSYSV returns ",ii,elapsed(t1)
!        
!        
!        do ii = 1,N            
!            r(ii) = dot_product( A0(ii,:),b(:) )
!            write(*,fmt='(3f12.6)') b0(ii),b(ii),r(ii)
!        end do
!                                     
!        r = b0 - r
!        print *,"SSYSV error ",norm2(r)/N
!        print *,""
!        
!        eps = 1.0d-8
!        bS = real(b0,kind=real64)
!        xS = 0.0d0
!        t2 = Callipers_ctor()
!        call conjgrad( AS,indx, xS,bS , eps )
!        
!        print *,"conjgrad returns ",eps,elapsed(t2)
!        
!        do ii = 1,N                        
!            r(ii) = dot_product( A0(ii,:),xs(:) )
!            write(*,fmt='(3f12.6)') xs(ii),r(ii)
!        end do
!        r = b0 - r                       
!        print *,"CG error ",norm2(r)/N
!        print *,""
!        print *,"done"
!        print *,""
!        
!        
!    end  program testLib_ConjugateGradient        
!        
         