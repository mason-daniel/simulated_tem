
    module Lib_FactoriseParallel
!---^^^^^^^^^^^^^^^^^^^^^^^^^^^^
!*      break up a rectangle size A x B 
!*      into N roughly equal-sized parts               

!*      Daniel Mason
!*      (c) UKAEA 2023

        use iso_fortran_env
        implicit none
        
        public          ::      factoriseParallel
        
        interface   factoriseParallel 
            module procedure        FactoriseParallel0
            module procedure        FactoriseParallel1
            module procedure        FactoriseParallel2
            module procedure        FactoriseParallel3 
        end interface
        
        interface   score
            module procedure        score0
            module procedure        score1
            module procedure        score2
            module procedure        score3
        end interface
        
        
    contains
!---^^^^^^^^
    
        subroutine FactoriseParallel0( A,B,N,C,D )
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      divide rectangle size AxB into (C,D) blocks with CxD = N
            real(kind=real64),intent(in)        ::      A,B
            integer,intent(in)                  ::      N
            integer,intent(out)                 ::      C,D
            
            integer,dimension(:),allocatable    ::      factor
            
            integer                             ::      ii,besti
            real(kind=real64)                   ::      ss,bestScore
            
            select case(N)
                case(1)
                    C = 1 ; D = 1 
                case(2)
                    if (A>=B) then
                        C = 2 ; D = 1
                    else
                        C = 1 ; D = 2 
                    end if
                case(3) 
                    if (A>=B) then
                        C = 3 ; D = 1 
                    else
                        C = 1 ; D = 3 
                    end if
                case(4)
                    if (A>=2*B) then
                        C = 4 ; D = 1 
                    else if (B>=2*A) then
                        C = 1 ; D = 4 
                    else 
                        C = 2 ; D = 2 
                    end if
                case default                                
                    call factorise( N,factor )
                    
                    bestScore = huge(1.0)
                    besti = 0
                    do ii = 1,size(factor)
                        C = factor(ii)
                        D = int(N/C)
                        ss = score(A,B,C,D)
                        if (ss < bestScore) then
                            bestScore = ss
                            besti = ii
                        end if
                    end do
                    C = factor(besti)
                    D = int(N/C)
                    deallocate(factor)
            end select
            
            return
        end subroutine FactoriseParallel0
                
    
        subroutine FactoriseParallel1( A,B,N,C,D )
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ 
    !*      divide rectangle size AxB into (C,D) blocks with CxD = N
            integer,intent(in)                  ::      A,B
            integer,intent(in)                  ::      N
            integer,intent(out)                 ::      C,D
            
            integer,dimension(:),allocatable    ::      factor
            
            integer                             ::      ii,besti
            real(kind=real64)                   ::      ss,bestScore
            
            select case(N)
                case(1)
                    C = 1 ; D = 1 
                case(2)
                    if (A>=B) then
                        C = 2 ; D = 1
                    else
                        C = 1 ; D = 2 
                    end if
                case(3) 
                    if (A>=B) then
                        C = 3 ; D = 1 
                    else
                        C = 1 ; D = 3 
                    end if
                case(4)
                    if (A>=2*B) then
                        C = 4 ; D = 1 
                    else if (B>=2*A) then
                        C = 1 ; D = 4 
                    else 
                        C = 2 ; D = 2 
                    end if
                case default                                
                    call factorise( N,factor )
                    
                    bestScore = huge(1.0)
                    besti = 0
                    do ii = 1,size(factor)
                        C = factor(ii)
                        D = int(N/C)
                        ss = score(A,B,C,D)
                        if (ss < bestScore) then
                            bestScore = ss
                            besti = ii
                        end if
                    end do
                    C = factor(besti)
                    D = int(N/C)
                    deallocate(factor)
            end select
            
            return
        end subroutine FactoriseParallel1
                
            
            
        subroutine FactoriseParallel2( A,B,C , N , D,E,F )
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      divide rectangle size AxBxC into (D,E,F) blocks with DxExF = N
            real(kind=real64),intent(in)        ::      A,B,C
            integer,intent(in)                  ::      N
            integer,intent(out)                 ::      D,E,F
            
            integer,dimension(:),allocatable    ::      factor1,factor2
            
            integer                             ::      ii,jj,besti,bestj
            real(kind=real64)                   ::      ss,bestScore
            
            select case(N)
                case(1)
                    D = 1 ; E = 1 ; F = 1
                case default                                
                    call factorise( N,factor1 )
                    
                    bestScore = huge(1.0)
                    besti = 0 ; bestj = 0
                    do ii = 1,size(factor1)
                        D = factor1(ii)
                        call factorise( int(N/D),factor2 )
                        do jj = 1,size(factor2)
                            E = factor2(jj)
                            F = int( int(N/D)/E )
                            ss = score(A,B,C , D,E,F)
                            if (ss < bestScore) then
                                bestScore = ss
                                besti = ii
                                bestj = jj
                            end if
                        end do
                        deallocate(factor2)
                    end do
                    D = factor1(besti)
                    call factorise( int(N/D),factor2 )
                    E = factor2(bestj)
                    F = int( int(N/D)/E )
                    deallocate(factor1)
            end select
            
            return
        end subroutine FactoriseParallel2
                
    
        subroutine FactoriseParallel3( A,B,C , N , D,E,F )
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ 
    !*      divide rectangle size AxBxC into (D,E,F) blocks with DxExF = N
            integer,intent(in)                  ::      A,B,C
            integer,intent(in)                  ::      N
            integer,intent(out)                 ::      D,E,F
            
            integer,dimension(:),allocatable    ::      factor1,factor2
            
            integer                             ::      ii,besti,jj,bestj
            real(kind=real64)                   ::      ss,bestScore
            
            select case(N)
                case(1)
                    D = 1 ; E = 1 ; F = 1
                case default                                
                    call factorise( N,factor1 )
                    
                    bestScore = huge(1.0)
                    besti = 0 ; bestj = 0
                    
                    !print *,"factors ",factor1
                    
                    do ii = 1,size(factor1)
                        D = factor1(ii)
                        call factorise( int(N/D),factor2 )
                        
                        !print *,"factor ",ii,factor1(ii)," factors ",factor2
                        do jj = 1,size(factor2)
                            E = factor2(jj)
                            F = int( int(N/D)/E )
                            !print *,"trial D,E,F = ",D,E,F
                            ss = score(A,B,C , D,E,F)
                            if (ss < bestScore) then
                                bestScore = ss
                                besti = ii
                                bestj = jj
                            end if
                        end do
                        deallocate(factor2)
                    end do
                    D = factor1(besti)
                    call factorise( int(N/D),factor2 )
                    E = factor2(bestj)
                    F = int( int(N/D)/E )
                    deallocate(factor1)
            end select
            
            return
        end subroutine FactoriseParallel3
                
            
              

        subroutine factorise( N,factor )
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      find the possible factors of N
    !*      not expecting N very large
        
            integer,intent(in)                  ::      N
            integer,dimension(:),allocatable    ::      factor
            
            logical,dimension(N)        ::      ok
            integer                     ::      ii,jj
            
            ok = .false.
            
            do ii = 1,ceiling( sqrt( real(N) ) )
                 
                if (mod(N,ii)==0) then
                    ok(ii) = .true.
                    jj = int( N/ii )
                    ok(jj) = .true.
                end if
            end do
            
            allocate(factor(count(ok)))
            jj = 0
            do ii = 1,N
                if (ok(ii)) then
                    jj = jj + 1
                    factor(jj) = ii
                end if
            end do
            
            return
        end subroutine factorise
        
        
        pure real(kind=real64) function score0( A,B,C,D )        
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      divide rectangle size AxB into (C,D) blocks
    !*      return score = max dimension / min dimension
            real(kind=real64),intent(in)            ::      A,B
            integer,intent(in)                      ::      C,D
            
            real(kind=real64)       ::      len1,len2
            
            len1 = A/C
            len2 = B/D
            
            if (len1>=len2) then
                score0 = len1/len2
            else
                score0 = len2/len1
            end if
            
            return
        end function score0
            
        pure real(kind=real64) function score1( A,B,C,D )        
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      divide rectangle size AxB into (C,D) blocks
    !*      return score = max dimension / min dimension
            integer,intent(in)                  ::      A,B
            integer,intent(in)                  ::      C,D
            
          ! integer                 ::      len1,len2 
            
            integer         ::      len1,len2,len1m,len2m
            len1 = max(1,int(A/C))
            len2 = max(1,int(B/D)) 
            
            len1m = len1 ; if (mod(A,C)/=0) len1m = len1m + 1
            len2m = len2 ; if (mod(B,D)/=0) len2m = len2m + 1
             
            score1 = real( max(len1m,len2m)) / min(len1,len2) 
                     
            
            return
        end function score1

        
        pure real(kind=real64) function score2( A,B,C , D,E,F )        
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      divide rectangle size AxBxC into (C,D,E) blocks
    !*      return score = max dimension / min dimension
            real(kind=real64),intent(in)            ::      A,B,C
            integer,intent(in)                      ::      D,E,F
            
            real(kind=real64)       ::      len1,len2,len3
            
            len1 = A/D
            len2 = B/E
            len3 = C/F
            
            
             
            
            score2 = real(maxval( (/len1,len2,len3/) ))/minval( (/len1,len2,len3/) )
            
            return
        end function score2
            
        pure real(kind=real64) function score3( A,B,C , D,E,F )        
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      divide rectangle size AxBxC into (C,D,E) blocks
    !*      return score = max dimension / min dimension
            integer,intent(in)                  ::      A,B,C
            integer,intent(in)                  ::      D,E,F
            
            integer         ::      len1,len2,len3,len1m,len2m,len3m
            len1 = max(1,int(A/D))
            len2 = max(1,int(B/E)) 
            len3 = max(1,int(C/F)) 
            
            len1m = len1 ; if (mod(A,D)/=0) len1m = len1m + 1
            len2m = len2 ; if (mod(B,E)/=0) len2m = len2m + 1
            len3m = len3 ; if (mod(C,F)/=0) len3m = len3m + 1
            !print *,"score3 ABC",A,B,C,"DEF",D,E,F," len ",len1,len2,len3,len1m,len2m,len3m
            
            score3 = real(maxval( (/len1m,len2m,len3m/) ))/minval( (/len1,len2,len3/) )
                    
!              
!                         
            
            return
        end function score3

        
        
        
        
        
    end module Lib_FactoriseParallel                            