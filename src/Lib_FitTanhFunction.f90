
    module Lib_FitTanhFunction
!---^^^^^^^^^^^^^^^^^^^^^^^^^^  
!*      given data ( x_i,f_i ) where x_i is a vector position in d dimensions
!*      and f_i is a function value between -1 and +1
!*      fit to f = tanh( x.b + c )
!*      where b is a vector direction in d dimensions

!*      This is done in two parts - first an approximate solution is found to the linear fit f = x.b + c
!*      then this is iterated to convergence.
!*      Both parts are least squares fitting
!*      with the first part being linear and the second non-linear. 

        use iso_fortran_env
        implicit none
        private
        
        external    ::      DSYSV 
                
        
        public          ::      fitTanhFunction
        
        interface       fitTanhFunction 
            module procedure    fitTanhFunction1
            module procedure    fitTanhFunction2
        end interface
        
        
    contains
!---^^^^^^^^

        subroutine linearFit( x,f , b,c )
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      perform a linear fit to the data f = x.b + c
    !*      This can be done with linear least squares and becomes a simple set of linear equations
    !*      and so is a good way to find a preliminary guess
            real(kind=real64),dimension(:,:),intent(in)             ::      x
            real(kind=real64),dimension(:),intent(in)               ::      f
            real(kind=real64),dimension(size(x,dim=1)),intent(out)  ::      b
            real(kind=real64),intent(out)                           ::      c
            
            integer             ::      nn,dd
            real(kind=real64),dimension(:,:),allocatable        ::      aa
            real(kind=real64),dimension(:),allocatable          ::      bb
            real(kind=real64),dimension(:),allocatable          ::      work
            integer,dimension(:),allocatable                    ::      ipiv
             
            
            
            
            
            
            
            integer             ::      ii,jj,kk
            
            dd = size(x,dim=1)      !   dimensionality
            nn = size(x,dim=2)      !   number of data points
            
        !---    solve as linear least squares
        !   S = 1/2 sum_k ( ( x_k.b + c ) - f_k )^2
        !   find dS/db, dS/dc and set to zero - gives linear equations
            allocate(aa(dd+1,dd+1)) ; aa = 0 
            allocate(bb(dd+1))      ; bb = 0    
            do kk = 1,nn
                do jj = 1,dd
                    do ii = 1,jj
                        aa(ii,jj) = aa(ii,jj) + x(ii,kk)*x(jj,kk)
                    end do
                    aa(jj,dd+1) = aa(jj,dd+1) + x(jj,kk)
                    bb(jj) = bb(jj) + f(kk)*x(jj,kk)
                end do
            end do
            bb(dd+1) = sum( f(1:nn) )
            aa(dd+1,dd+1) = nn
            
        !---    solve
            allocate(ipiv(dd+1))
            allocate(work(66*(dd+1)))   
        !   do ii = 1,dd+1
        !       write (*,fmt='(f16.8,a,100f16.8)') bb(ii)," , ",aa(ii,:)
        !   end do

        
 
            call DSYSV( "U",dd+1,1,aa,dd+1,ipiv,bb,dd+1,work,size(work),ii )
 


            b = bb(1:dd) 
            c = bb(dd+1)
                
            
            return
        end subroutine linearFit
    
        
        
        subroutine fitTanhFunction1( x,f, b,c )
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      perform a fit to the data f = x.b + c
            real(kind=real64),dimension(:,:),intent(in)             ::      x
            real(kind=real64),dimension(:),intent(in)               ::      f
            real(kind=real64),dimension(size(x,dim=1)),intent(out)  ::      b
            real(kind=real64),intent(out)                           ::      c
            
            integer             ::      nn,dd
            real(kind=real64),dimension(:),allocatable          ::      tt
            real(kind=real64),dimension(:,:),allocatable        ::      aa
            real(kind=real64),dimension(:),allocatable          ::      bb
            real(kind=real64),dimension(:),allocatable          ::      work
            integer,dimension(:),allocatable                    ::      ipiv
            integer             ::      ii,jj,kk,step
            integer,parameter   ::      NSTEPS = 20
            real(kind=real64),parameter     ::      EPS = 1.0d-6
            real(kind=real64)   ::      qq,ss,s_old
            
            real(kind=real64),dimension(size(x,dim=1))  ::      b_old
            real(kind=real64)                           ::      c_old
             
            
            
            dd = size(x,dim=1)      !   dimensionality
            nn = size(x,dim=2)      !   number of data points
            
            
        !---    make a preliminary linear fit
            call linearFit( x,f , b,c )
           ! return
                        
            allocate(tt(nn))
            allocate(aa(dd+1,dd+1))  
            allocate(bb(dd+1))       
            allocate(work(66*(dd+1)))   
            allocate(ipiv(dd+1))

        !---    compute tanh( x.b + c ) and the error
            call fit( b,c,tt,ss )
            
        !---    store the initial guess
            s_old = ss
            b_old = b
            c_old = c       
            
            !print *,"solution step ",0,b_old,c_old,s_old
            
            do step = 1,NSTEPS*dd
                
                
                    
            !---    fit by minimising S = 1/2 sum_k ( t_k - f_k )^2
            !---    compute derivatives wrt b,c
                aa = 0 ; bb = 0
                do kk = 1,nn
                
                !---    first derivative. Note d/dx tanh x = (1-tanh^2 x)
                    qq = (tt(kk) - f(kk))*(1 - tt(kk)*tt(kk))
                    do jj = 1,dd
                        bb(jj) = bb(jj) + qq*x(jj,kk)
                    end do
                    bb(dd+1) = bb(dd+1) + qq
                    
                    
                    
                !---    second derivative. Note that x.b + c is linear, so this is really easy.
                !   note: only need upper triangle.
                    qq = 1 + tt(kk)*(2*f(kk) - 3*tt(kk))
                    do jj = 1,dd
                        do ii = 1,jj
                            aa(ii,jj) = aa(ii,jj) + qq*x(ii,kk)*x(jj,kk)
                        end do
                        aa(jj,dd+1) = aa(jj,dd+1) + qq*x(jj,kk)
                        aa(dd+1,dd+1) = aa(dd+1,dd+1) + qq
                    end do
                    
                end do
                
                !do ii = 1,dd+1
                !   write (*,fmt='(f16.8,a,100f16.8)') bb(ii)," , ",aa(ii,:)
                !end do
                
                
            !---    find an improved guess
                bb = - bb       !   want to solve S = S0 + S' dx + 1/2 S" dx dx
                                !   gives S" dx = - S'
                                
             call DSYSV( "U",dd+1,1,aa,dd+1,ipiv,bb,dd+1,work,size(work),ii )
    
                
                b = b_old + bb(1:dd)
                c = c_old + bb(dd+1)
                
            !---    compute tanh( x.b + c ) and the error
                call fit( b,c,tt,ss )
            
                
            !---    is this improved guess actually any better??
                if (ss>s_old) then
                    !   nope. Try a shorter step?                   
                    do kk = 1,4
                        bb = bb / 2
                        b = b_old + bb(1:dd) 
                        c = c_old + bb(dd+1) 
                        call fit( b,c,tt,ss )
                        if (ss<s_old) exit
                    end do
                end if
                
                !print *,"solution step ",step,b_old,c_old,s_old
                
                if (ss<s_old-EPS*EPS*nn) then
                    !   the solution is improving.
                    s_old = ss
                    b_old = b
                    c_old = c
                else if (ss<s_old) then
                    !   the solution is converged
                    exit
                else    
                    !   the solution is not improving
                    b = b_old
                    c = c_old
                    ss = s_old
                    exit
                end if
                
                
                
            end do
                
                
                            
            return
            
        contains
    !---^^^^^^^^
    
            subroutine fit( b,c,t,s )
        !---^^^^^^^^^^^^^^^^^^^^^^^^^
                real(kind=real64),dimension(:),intent(in)           ::      b
                real(kind=real64),intent(in)                        ::      c
                real(kind=real64),dimension(:),intent(out)          ::      t
                real(kind=real64),intent(out)                       ::      s
            
                real(kind=real64)           ::      qq
                integer                     ::      kk,jj
                s = 0.0
                do kk = 1,nn                             
                    qq = c
                    do jj = 1,dd
                        qq = qq + x(jj,kk)*b(jj)
                    end do
                    t(kk) = tanh( qq )
                    qq = t(kk) - f(kk)
                    s = s + qq*qq
                end do  
                s = s/2
                return
            end subroutine fit
            
            
        end subroutine fitTanhFunction1

        
        subroutine fitTanhFunction2( x,f, b,c , rc,weight)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      perform a fit to the data f = x.b + c
    !*      with a hint that the width should be ~ rc
            real(kind=real64),dimension(:,:),intent(in)             ::      x
            real(kind=real64),dimension(:),intent(in)               ::      f
            real(kind=real64),dimension(size(x,dim=1)),intent(out)  ::      b
            real(kind=real64),intent(out)                           ::      c
            real(kind=real64),intent(in)                            ::      rc
            real(kind=real64),intent(in),optional                   ::      weight
            
            integer             ::      nn,dd 
            real(kind=real64),dimension(:),allocatable          ::      tt
            real(kind=real64),dimension(:,:),allocatable        ::      aa
            real(kind=real64),dimension(:),allocatable          ::      bb  
            real(kind=real64),dimension(:),allocatable          ::      work
            integer,dimension(:),allocatable                    ::      ipiv
            
            real(kind=real64),dimension(size(x,dim=1))      ::      avec
        
            real(kind=real64)                               ::      a0
            integer             ::      ii,jj,kk,step
            integer,parameter   ::      NSTEPS = 20
            real(kind=real64),parameter     ::      EPS = 1.0d-6
            real(kind=real64)   ::      qq,ss,s_old,iq2,ww 
            
            real(kind=real64),dimension(size(x,dim=1))  ::      b_old
            real(kind=real64)                           ::      c_old
            
             
            dd = size(x,dim=1)      !   dimensionality
            nn = size(x,dim=2)      !   number of data points
            ww = 1.0d0
            if (present(weight)) ww = weight
                 
        !---    make a preliminary linear fit
            call linearFit( x,f , b,c )
            
                        
            allocate(tt(nn))
            allocate(aa(dd+1,dd+1))  
            allocate(bb(dd+1))
            allocate(work(66*(dd+1)))   
            allocate(ipiv(dd+1))
            
           

        !---    compute tanh( x.b + c ) and the error
            call fit2( b,c,tt,ss,rc,ww )
            
        !---    store the initial guess
            s_old = ss
            b_old = b
            c_old = c       
            
            !print *,"solution step ",0,b_old,c_old,s_old
            
            do step = 1,NSTEPS*dd
                
                
                    
            !---    fit by minimising S = 1/2 sum_k ( t_k - f_k )^2
            !---    compute derivatives wrt b,c
                aa = 0 ; bb = 0
                a0 = 0.0d0 ; avec = 0.0d0 
                do kk = 1,nn
                
                !---    first derivative. Note d/dx tanh x = (1-tanh^2 x)
                    qq = (tt(kk) - f(kk))*(1 - tt(kk)*tt(kk))
                    bb(1:dd) = bb(1:dd) + qq*x(1:dd,kk)
                    bb(dd+1) = bb(dd+1) + qq
                    
                    
                    
                !---    second derivative. Note that x.b + c is linear, so this is really easy.
                !   note: only need upper triangle.
                    qq = 1 + tt(kk)*(2*f(kk) - 3*tt(kk))
                    
                     do jj = 1,dd
 !                         do ii = 1,jj
 !                             aa(ii,jj) = aa(ii,jj) + qq*x(ii,kk)*x(jj,kk)
 !                         end do
                         aa(1:jj,jj)   = aa(1:jj,jj) + qq*x(1:jj,kk)*x(jj,kk)
                        ! aa(jj,dd+1)   = aa(jj,dd+1)             + qq*x(jj,kk)
                        ! aa(dd+1,dd+1) = aa(dd+1,dd+1) + qq
                         
                     end do
 
                    a0 = a0 + dd*qq
                    avec(1:dd) = avec(1:dd) + qq*x(1:dd,kk)
                    
                    
                end do
                aa(1:dd,dd+1) = avec(1:dd)
                aa(dd+1,dd+1) = a0
                
          
             !--    add contribution S += weight/2 (1/|b| - rc)^2  
                qq = norm2(bb)
                if (qq>0) then      
                    iq2 = 1/(qq*qq)
                    
                    do jj = 1,dd
                        bb(jj) = bb(jj) + ww*(qq*rc - 1)*b(jj)*iq2
                        do ii = 1,jj-1
                            aa(ii,jj) = aa(ii,jj) + ww*b(ii)*b(jj)*(4-3*qq*rc)*iq2*iq2*iq2
                        end do
                        aa(jj,jj) = aa(jj,jj) + ww*(qq*qq*(qq*rc-1)+b(jj)*b(jj)*(4-3*qq*rc))*iq2*iq2*iq2
                        
                    end do
                end if
                
                
                !do ii = 1,dd+1
                !   write (*,fmt='(f16.8,a,100f16.8)') bb(ii)," , ",aa(ii,:)
                !end do
                
                
            !---    find an improved guess
                bb = - bb       !   want to solve S = S0 + S' dx + 1/2 S" dx dx
                                !   gives S" dx = - S'
 
                call DSYSV( "U",dd+1,1,aa,dd+1,ipiv,bb,dd+1,work,size(work),ii )
         
                
                b = b_old + bb(1:dd)
                c = c_old + bb(dd+1)
                
            !---    compute tanh( x.b + c ) and the error
                call fit2( b,c,tt,ss,rc,ww )
            
                
            !---    is this improved guess actually any better??
                if (ss>s_old) then
                    !   nope. Try a shorter step?                   
                    do kk = 1,4
                        bb = bb / 2
                        b = b_old + bb(1:dd) 
                        c = c_old + bb(dd+1) 
                        call fit2( b,c,tt,ss,rc,ww )
                        if (ss<s_old) exit
                    end do
                end if
                
                !print *,"solution step ",step,b_old,c_old,s_old
                
                if (ss<s_old-EPS*EPS*nn) then
                    !   the solution is improving.
                    s_old = ss
                    b_old = b
                    c_old = c
                else if (ss<s_old) then
                    !   the solution is converged
                    exit
                else    
                    !   the solution is not improving
                    b = b_old
                    c = c_old
                    ss = s_old
                    exit
                end if
                
                
                
            end do
                 
                
                            
            return
            
        contains
    !---^^^^^^^^
    
            subroutine fit2( b,c,t,s,rc,w )
        !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
                real(kind=real64),dimension(:),intent(in)           ::      b
                real(kind=real64),intent(in)                        ::      c
                real(kind=real64),dimension(:),intent(out)          ::      t
                real(kind=real64),intent(out)                       ::      s
                real(kind=real64),intent(in)                        ::      rc,w
                real(kind=real64)           ::      qq
                integer                     ::      kk
                s = 0.0
                do kk = 1,nn                             
                    qq = c + dot_product(x(1:dd,kk),b(1:dd))                   
                    t(kk) = tanh( qq )
                    qq = t(kk) - f(kk)
                    s = s + qq*qq
                end do 
                qq = norm2(b)
                if (qq>0) then
                    qq = (1/qq - rc)
                    s  = s + w*qq*qq/2
                end if
                s = s/2
                return
            end subroutine fit2
            
            
        end subroutine fitTanhFunction2

    
    end module Lib_FitTanhFunction  
    
    
!   gfortran -ffree-line-length-256 Lib_FitTanhFunction.f90 -llapack
