
    module Lib_LocalFilter2d 
!---^^^^^^^^^^^^^^^^^^^^^^^^
!*      given the values of points on a 2d grid 
!*      find the interpolated value at the central point
!*      
!*  
!*      |---|---|---|---|         fit the values on the nodes surrounding 0             8---5---4---5---8
!*      |   |   |   |   |         by minimising the least squares fit                   |   |   |   |   |
!*      |---|---|---|---|           S = sum_i w_i (f(x_i) - f_i)^2                      5---2---1---2---5
!*      |   |   |   |   |                                                               |   |   |   |   |
!*      |---|---0---|---|         where f(x) = f0 + f'.x + 1/2 x. f" x                  4---1---0---1---4
!*      |   |   |   |   |                                                               |   |   |   |   |
!*      |---|---|---|---|         and w_i = exp[ - x^2/(2 s^2) ]                        5---2---1---2---5
!*      |   |   |   |   |                                                               |   |   |   |   |
!*      |---|---|---|---|         ignore any pixel for which the value                  8---5---4---5---8
!*                                 = LOCALFILTER2D_IGNORE                              values of x^2 on nodes
!*
!*      [ optional first step ]
!*          there is optional weighing w_i in the least squares fit.
!*          set the weighting to be gaussian with
!*      call setLocalFilter2d_sigma( [s=1] )
!*          where if no argument is given sigma is set to 1,
!*          or if argument is zero then weighting is set to w_i = 1
!*          or sigma is set to max( s,1/3 )
!*
!*      f0 = interpolate55( f )     
!*          makes an interpolation using the 5x5 grid f
!*
!*      call denoise( f_in, sig, f_out, stdev )
!*          finds the standard deviation difference between input and smoothed pixels
!*          and flattens any pixels greater than sig std devs.
!*          note call localFilter( f_in, f_out ) is a denoise with sig = 0 - ie flatten all pixels.
!*      

        use iso_fortran_env
        implicit none
        private
        
        
    !---
    
        external        ::      DGESV,SGESV           !   lapack solver for linear-least-squares problem. ******* WARNING - DGESV PROBLEMS WORKAROUND WITH SGESV ********
        
        
    !---    function calls provided
    
        public      ::      setLocalFilter2d_sigma
        public      ::      interpolate55
        public      ::      denoise
        public      ::      localFilter
        public      ::      doubleSizeImage,halfSizeImage
        
    !---    define a 64bit floating constant "LOCALFILTER2D_IGNORE" with a really unlikely bit pattern        
        integer(kind=int64),private,parameter               ::      BADF00D = int( z'BADF00D',kind=int64 )
        real(kind=real64),public,parameter                  ::      LOCALFILTER2D_IGNORE = transfer( (BADF00D+ishft(BADF00D,32_int64)),1.0d0 )
        real(kind=real64),public,parameter                  ::      LOCALFILTER2D_IGNORE32 = transfer( (BADF00D),1.0_real32 )
    
                                                                 
    !   weighing of grid nodes                                                                                     
        real(kind=real64),private                           ::      sigma = 1.0d0
        real(kind=real64),dimension(-2:2,-2:2),private      ::      w = reshape( (/   0.018315639d0, 0.082084999d0, 0.135335283d0, 0.082084999d0, 0.018315639d0,  &                 !   weights if sigma = 1
                                                                                      0.082084999d0, 0.367879441d0, 0.606530660d0, 0.367879441d0, 0.082084999d0,  &
                                                                                      0.135335283d0, 0.606530660d0, 0.000000000d0, 0.606530660d0, 0.135335283d0,  &
                                                                                      0.082084999d0, 0.367879441d0, 0.606530660d0, 0.367879441d0, 0.082084999d0,  &
                                                                                      0.018315639d0, 0.082084999d0, 0.135335283d0, 0.082084999d0, 0.018315639d0 /),(/5,5/) )
      
    !---
    
        interface   localFilter
            module procedure    localFilter0
            module procedure    localFilter1
        end interface                                                                                      
                              
        interface   interpolate55
            module procedure    interpolate550
            module procedure    interpolate551
            module procedure    interpolate55_stripe
        end interface                                                                                      
                                                                                      
        
        interface   denoise
            module procedure    denoise0
            module procedure    denoise1
            module procedure    denoise_stripe
        end interface
    !---
    
        logical,private     ::      dbg = .false.
        integer,public     ::      LIB_LF2D_DBGX = -1 , LIB_LF2D_DBGY = -1
        
    contains
!---^^^^^^^^

        subroutine setLocalFilter2d_sigma(s)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      change the Gaussian sampling window.
    !*      note there is no need to normalise
            real(kind=real64),intent(in),optional        ::  s
            real(kind=real64)       ::      i2s2
            integer     ::      ii,jj
            if (.not. present(s)) then
                sigma = 1.0d0
                do jj = -2,2
                    do ii = -2,2
                        w(ii,jj) = exp( - (ii*ii+jj*jj)*0.5d0 )
                    end do
                end do          
            else
                if (s == 0) then
                    sigma = 0.0d0
                    w = 1.0d0
                else 
                    sigma = s
                    sigma = max( s,1.0d0/3 )
                    i2s2 = 1/(2*sigma*sigma)
                    do jj = -2,2
                        do ii = -2,2
                            w(ii,jj) = exp( - (ii*ii+jj*jj)*i2s2 )
                        end do
                    end do          
                end if
            end if
            w(0,0) = 0.0d0
            return
        end subroutine setLocalFilter2d_sigma             
        
    !---    
    
        function interpolate550( f ) result( f0 )
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^                       
            real(kind=real64),dimension(-2:2,-2:2),intent(in)       ::      f
            real(kind=real64)                                       ::      f0
            
            real(kind=real64),dimension(-2:2,-2:2)      ::      ww
            
            real(kind=real64),dimension(6,6)        ::      AA
            real(kind=real64),dimension(6)          ::      bb
            integer,dimension(6)                    ::      ipiv
            integer                                 ::      ii,jj,ierror
            integer                                 ::      ignoreCount
            
            real(kind=real64)       ::      sumw,sumwx,sumwy
            real(kind=real64)       ::      sumwxx,sumwxy,sumwyy
            real(kind=real64)       ::      sumwxxx,sumwxxy,sumwxyy,sumwyyy
            real(kind=real64)       ::      sumwxxxx,sumwxxxy,sumwxxyy,sumwxyyy,sumwyyyy
            real(kind=real64)       ::      sumwf,sumwfx,sumwfy
            real(kind=real64)       ::      sumwfxx,sumwfxy,sumwfyy
            
            real(kind=real64)       ::      wwww
            
                    real(kind=real64)                           ::      ss
                    real(kind=real64),dimension(-2:2,-2:2)      ::      gg
                    
                    
        !********   PROBLEMS WITH DGESV WORKAROUND
            real(kind=real32),dimension(6,6)        ::      AA32
            real(kind=real32),dimension(6)          ::      bb32                   
                    
                    
                    
                    
                    
                    
                    
                    
                    
        !---    set up the weighting on the nodes to use.
            ignoreCount = 0
            do jj = -2,2
                do ii = -2,2
                    if (f(ii,jj)==LOCALFILTER2D_IGNORE) then
                        ww(ii,jj) = 0.0d0
                        ignoreCount = ignoreCount + 1
                    else
                        ww(ii,jj) = w(ii,jj)
                    end if
                end do
            end do
            
            if (dbg) print *,"ignoreCount ",ignoreCount
            if (f(0,0)/= LOCALFILTER2D_IGNORE) then
                if (ignoreCount >= 10) ww(0,0) = 0.5d0
                if (ignoreCount >= 12) ww(0,0) = 1.0d0
            end if
            if (ignoreCount >= 14) then
                sumw     = 0; sumwf    = 0 
                do jj = -2,2
                    do ii = -2,2
                        sumw     = sumw     + ww(ii,jj)
                        sumwf    = sumwf    + ww(ii,jj)*f(ii,jj) 
                    end do
                end do
                if (sumw>0) then
                    f0 = sumwf/sumw
                else
                    f0 = LOCALFILTER2D_IGNORE
                end if
                return
            end if
                
        !---    now need to solve the set of linear equations A x = b for x
        !   find matrix A and rhs b
            sumw     = 0; sumwx    = 0; sumwy    = 0 
            sumwxx   = 0; sumwxy   = 0; sumwyy   = 0 
            sumwxxx  = 0; sumwxxy  = 0; sumwxyy  = 0 ; sumwyyy  = 0
            sumwxxxx = 0; sumwxxxy = 0; sumwxxyy = 0 ; sumwxyyy = 0 ; sumwyyyy = 0
            sumwf    = 0; sumwfx   = 0; sumwfy   = 0 
            sumwfxx  = 0; sumwfxy  = 0; sumwfyy  = 0 
            do jj = -2,2
                do ii = -2,2
                    wwww     = ww(ii,jj)
                    sumw     = sumw     + wwww
                    sumwx    = sumwx    + wwww*ii
                    sumwxx   = sumwxx   + wwww*ii*ii
                    sumwxxx  = sumwxxx  + wwww*ii*ii*ii
                    sumwxxxx = sumwxxxx + wwww*ii*ii*ii*ii
                    sumwy    = sumwy    + wwww*jj    
                    sumwyy   = sumwyy   + wwww*jj*jj
                    sumwyyy  = sumwyyy  + wwww*jj*jj*jj
                    sumwyyyy = sumwyyyy + wwww*jj*jj*jj*jj
                    sumwxy   = sumwxy   + wwww*ii*jj
                    sumwxyy  = sumwxyy  + wwww*ii*jj*jj
                    sumwxyyy = sumwxyyy + wwww*ii*jj*jj*jj
                    sumwxxy  = sumwxxy  + wwww*ii*ii*jj
                    sumwxxyy = sumwxxyy + wwww*ii*ii*jj*jj
                    sumwxxxy = sumwxxxy + wwww*ii*ii*ii*jj
                    wwww     = wwww*f(ii,jj)
                    sumwf    = sumwf    + wwww 
                    sumwfx   = sumwfx   + wwww*ii
                    sumwfxx  = sumwfxx  + wwww*ii*ii
                    sumwfxy  = sumwfxy  + wwww*ii*jj
                    sumwfy   = sumwfy   + wwww*jj 
                    sumwfyy  = sumwfyy  + wwww*jj*jj
                end do
            end do
            AA(1:6,1) = (/ sumw    ,sumwx    ,sumwy    ,sumwxx    ,sumwxy    ,sumwyy     /)
            AA(1:6,2) = (/ sumwx   ,sumwxx   ,sumwxy   ,sumwxxx   ,sumwxxy   ,sumwxyy    /)
            AA(1:6,3) = (/ sumwy   ,sumwxy   ,sumwyy   ,sumwxxy   ,sumwxyy   ,sumwyyy    /)
            AA(1:6,4) = (/ sumwxx/2,sumwxxx/2,sumwxxy/2,sumwxxxx/2,sumwxxxy/2,sumwxxyy/2 /)
            AA(1:6,5) = (/ sumwxy/2,sumwxxy/2,sumwxyy/2,sumwxxxy/2,sumwxxyy/2,sumwxyyy/2 /)
            AA(1:6,6) = (/ sumwyy/2,sumwxyy/2,sumwyyy/2,sumwxxyy/2,sumwxyyy/2,sumwyyyy/2 /)
            bb(1:6)   = (/ sumwf   ,sumwfx   ,sumwfy   ,sumwfxx   ,sumwfxy   ,sumwfyy    /)        
                    
            if (dbg) then
                print *,"Lib_LocalFilter2d::interpolate55 info " 
                do ii = -2,2
                    write (*,fmt='(5f16.8,a,5f16.8)') f(ii,:)," , ",ww(ii,:)
                end do                
                do ii = 1,6
                    write (*,fmt='(6f16.8,a,f16.8)') AA(ii,:)," , ",bb(ii)
                end do
            end if             
                                
            
!***********        PROBLEMS WITH DGESV 15/03/22
!*          WORKAROUND WITH SGESV
#ifndef LAPACK_DSV_WORKAROUND
            call DGESV( 6,1,AA,6,ipiv,bb,6, ierror ) 
#else            
            AA32 = real(AA,kind=real32)
            bb32 = real(bb,kind=real32)
            call SGESV( 6,1,AA32,6,ipiv,bb32,6, ierror )             
            bb = real(bb32,kind=real64)
#endif            
!           END OF WORKAROUND
!***********           
            
            
            if (dbg) then
                print *,"Lib_LocalFilter2d::interpolate55 info - DGESV ",ierror,bb                               
                ss = 0.0d0
                do jj = -2,2
                    do ii = -2,2
                        gg(ii,jj) = (bb(1)+bb(2)*ii+bb(3)*jj+(bb(4)*ii*ii+2*bb(5)*ii*jj+bb(6)*jj*jj)/2)
                        ss = ss + ww(ii,jj)*(f(ii,jj)-gg(ii,jj))**2 
                    end do
                end do
                print *,"Lib_LocalFilter2d::interpolate55 info - error ",ss
                 do ii = -2,2               
                     write (*,fmt='(3(5f8.3,a))') f(ii,:)," , ",ww(ii,:)," , ",gg(ii,:)
                 end do
            end if           
            
            
            if (ierror /= 0) then
            !   could not find a good solution. This could be because there aren't enough non-zero coefficients to fit quadratic...
            !   try again, this time just with simple linear gradient solution.
                AA(1:4,1) = (/ sumw    ,sumwx    ,sumwy    ,sumwxx     /)
                AA(1:4,2) = (/ sumwx   ,sumwxx   ,sumwxy   ,sumwxxx    /)
                AA(1:4,3) = (/ sumwy   ,sumwxy   ,sumwyy   ,sumwxxy    /)
                AA(1:4,4) = (/ sumwxx/2,sumwxxx/2,sumwxxy/2,sumwxxxx/2 /)
                bb(1:4)   = (/ sumwf ,sumwfx ,sumwfy ,sumwfxx    /)     
                
             
            
!***********        PROBLEMS WITH DGESV
!           WORKAROUND WITH SGESV      

#ifndef LAPACK_DSV_WORKAROUND
            call DGESV( 4,1,AA(1:4,1:4),4,ipiv,bb(1:4),4, ierror ) 
#else            
            AA32 = real(AA,kind=real32)
            BB32 = real(BB,kind=real32)
            call SGESV( 4,1,AA32(1:4,1:4),4,ipiv,bb32(1:4),4, ierror )  
            bb = real(BB32,kind=real64)
#endif        
            
!           END OF WORKAROUND
!***********
                            
                
                
                
                if (dbg)  print *,"Lib_LocalFilter2d::interpolate55 info - DGESV ",ierror,bb(1:4)
                if (ierror /= 0) then
                !   still can't find a good solution. Try to find a weighted average?
                    bb(1)   = LOCALFILTER2D_IGNORE
                end if
            end if                       
                    
            f0 = bb(1)        
            !stop
            return
        end function interpolate550
            
    
    
        function interpolate551( f ) result( f0 )
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^                       
            real(kind=real32),dimension(-2:2,-2:2),intent(in)       ::      f
            real(kind=real32)                                       ::      f0
            
            real(kind=real32),dimension(-2:2,-2:2)      ::      ww
            
            real(kind=real32),dimension(6,6)        ::      AA
            real(kind=real32),dimension(6)          ::      bb
            integer,dimension(6)                    ::      ipiv
            integer                                 ::      ii,jj,ierror
            integer                                 ::      ignoreCount
            
            real(kind=real32)       ::      sumw,sumwx,sumwy
            real(kind=real32)       ::      sumwxx,sumwxy,sumwyy
            real(kind=real32)       ::      sumwxxx,sumwxxy,sumwxyy,sumwyyy
            real(kind=real32)       ::      sumwxxxx,sumwxxxy,sumwxxyy,sumwxyyy,sumwyyyy
            real(kind=real32)       ::      sumwf,sumwfx,sumwfy
            real(kind=real32)       ::      sumwfxx,sumwfxy,sumwfyy
            
            real(kind=real32)       ::      wwww
            
                    real(kind=real32)                           ::      ss
                    real(kind=real32),dimension(-2:2,-2:2)      ::      gg
                    
        !---    set up the weighting on the nodes to use.
            ignoreCount = 0
            do jj = -2,2
                do ii = -2,2
                    if (f(ii,jj)==LOCALFILTER2D_IGNORE32) then
                        ww(ii,jj) = 0.0d0
                        ignoreCount = ignoreCount + 1
                    else
                        ww(ii,jj) = real( w(ii,jj),kind=real32 )
                    end if
                end do
            end do
            
            if (dbg) print *,"ignoreCount ",ignoreCount
            if (f(0,0)/= LOCALFILTER2D_IGNORE32) then
                if (ignoreCount >= 10) ww(0,0) = 0.5 
                if (ignoreCount >= 12) ww(0,0) = 1.0 
            end if
            if (ignoreCount >= 14) then
                sumw     = 0; sumwf    = 0 
                do jj = -2,2
                    do ii = -2,2
                        sumw     = sumw     + ww(ii,jj)
                        sumwf    = sumwf    + ww(ii,jj)*f(ii,jj) 
                    end do
                end do
                if (sumw>0) then
                    f0 = sumwf/sumw
                else
                    f0 = LOCALFILTER2D_IGNORE32
                end if
                return
            end if
                
        !---    now need to solve the set of linear equations A x = b for x
        !   find matrix A and rhs b
            sumw     = 0; sumwx    = 0; sumwy    = 0 
            sumwxx   = 0; sumwxy   = 0; sumwyy   = 0 
            sumwxxx  = 0; sumwxxy  = 0; sumwxyy  = 0 ; sumwyyy  = 0
            sumwxxxx = 0; sumwxxxy = 0; sumwxxyy = 0 ; sumwxyyy = 0 ; sumwyyyy = 0
            sumwf    = 0; sumwfx   = 0; sumwfy   = 0 
            sumwfxx  = 0; sumwfxy  = 0; sumwfyy  = 0 
            do jj = -2,2
                do ii = -2,2
                    wwww     = ww(ii,jj)
                    sumw     = sumw     + wwww
                    sumwx    = sumwx    + wwww*ii
                    sumwxx   = sumwxx   + wwww*ii*ii
                    sumwxxx  = sumwxxx  + wwww*ii*ii*ii
                    sumwxxxx = sumwxxxx + wwww*ii*ii*ii*ii
                    sumwy    = sumwy    + wwww*jj    
                    sumwyy   = sumwyy   + wwww*jj*jj
                    sumwyyy  = sumwyyy  + wwww*jj*jj*jj
                    sumwyyyy = sumwyyyy + wwww*jj*jj*jj*jj
                    sumwxy   = sumwxy   + wwww*ii*jj
                    sumwxyy  = sumwxyy  + wwww*ii*jj*jj
                    sumwxyyy = sumwxyyy + wwww*ii*jj*jj*jj
                    sumwxxy  = sumwxxy  + wwww*ii*ii*jj
                    sumwxxyy = sumwxxyy + wwww*ii*ii*jj*jj
                    sumwxxxy = sumwxxxy + wwww*ii*ii*ii*jj
                    wwww     = wwww*f(ii,jj)
                    sumwf    = sumwf    + wwww 
                    sumwfx   = sumwfx   + wwww*ii
                    sumwfxx  = sumwfxx  + wwww*ii*ii
                    sumwfxy  = sumwfxy  + wwww*ii*jj
                    sumwfy   = sumwfy   + wwww*jj 
                    sumwfyy  = sumwfyy  + wwww*jj*jj
                end do
            end do
            AA(1:6,1) = (/ sumw    ,sumwx    ,sumwy    ,sumwxx    ,sumwxy    ,sumwyy     /)
            AA(1:6,2) = (/ sumwx   ,sumwxx   ,sumwxy   ,sumwxxx   ,sumwxxy   ,sumwxyy    /)
            AA(1:6,3) = (/ sumwy   ,sumwxy   ,sumwyy   ,sumwxxy   ,sumwxyy   ,sumwyyy    /)
            AA(1:6,4) = (/ sumwxx/2,sumwxxx/2,sumwxxy/2,sumwxxxx/2,sumwxxxy/2,sumwxxyy/2 /)
            AA(1:6,5) = (/ sumwxy/2,sumwxxy/2,sumwxyy/2,sumwxxxy/2,sumwxxyy/2,sumwxyyy/2 /)
            AA(1:6,6) = (/ sumwyy/2,sumwxyy/2,sumwyyy/2,sumwxxyy/2,sumwxyyy/2,sumwyyyy/2 /)
            bb(1:6)   = (/ sumwf   ,sumwfx   ,sumwfy   ,sumwfxx   ,sumwfxy   ,sumwfyy    /)        
                    
            if (dbg) then
                print *,"Lib_LocalFilter2d::interpolate55 info " 
                do ii = -2,2
                    write (*,fmt='(5f16.8,a,5f16.8)') f(ii,:)," , ",ww(ii,:)
                end do                
                do ii = 1,6
                    write (*,fmt='(6f16.8,a,f16.8)') AA(ii,:)," , ",bb(ii)
                end do
            end if             
                                
            call SGESV( 6,1,AA,6,ipiv,bb,6, ierror ) 
            
            if (dbg) then
                print *,"Lib_LocalFilter2d::interpolate55 info - SGESV ",ierror,bb                               
                ss = 0.0 
                do jj = -2,2
                    do ii = -2,2
                        gg(ii,jj) = (bb(1)+bb(2)*ii+bb(3)*jj+(bb(4)*ii*ii+2*bb(5)*ii*jj+bb(6)*jj*jj)/2)
                        ss = ss + ww(ii,jj)*(f(ii,jj)-gg(ii,jj))**2 
                    end do
                end do
                print *,"Lib_LocalFilter2d::interpolate55 info - error ",ss
                 do ii = -2,2               
                     write (*,fmt='(3(5f8.3,a))') f(ii,:)," , ",ww(ii,:)," , ",gg(ii,:)
                 end do
            end if           
            
            
            if (ierror /= 0) then
            !   could not find a good solution. This could be because there aren't enough non-zero coefficients to fit quadratic...
            !   try again, this time just with simple linear gradient solution.
                AA(1:4,1) = (/ sumw    ,sumwx    ,sumwy    ,sumwxx     /)
                AA(1:4,2) = (/ sumwx   ,sumwxx   ,sumwxy   ,sumwxxx    /)
                AA(1:4,3) = (/ sumwy   ,sumwxy   ,sumwyy   ,sumwxxy    /)
                AA(1:4,4) = (/ sumwxx/2,sumwxxx/2,sumwxxy/2,sumwxxxx/2 /)
                bb(1:4)   = (/ sumwf ,sumwfx ,sumwfy ,sumwfxx    /)     
                call SGESV( 4,1,AA(1:4,1:4),4,ipiv,bb(1:4),4, ierror ) 
                if (dbg)  print *,"Lib_LocalFilter2d::interpolate55 info - DGESV ",ierror,bb(1:4)
                if (ierror /= 0) then
                !   still can't find a good solution. Try to find a weighted average?
                    bb(1)   = LOCALFILTER2D_IGNORE32
                end if
            end if                       
                    
            f0 = bb(1)        
            !stop
            return
        end function interpolate551
            
    
        function interpolate55_stripe( f,horv ) result( f0 )
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^                       
            real(kind=real64),dimension(-2:2,-2:2),intent(in)       ::      f
            character(len=1),intent(in)                             ::      horv
            real(kind=real64)                                       ::      f0
            
            real(kind=real64),dimension(-2:2,-2:2)      ::      ww
            
            real(kind=real64),dimension(3,3)        ::      AA
            real(kind=real64),dimension(3)          ::      bb
            integer,dimension(3)                    ::      ipiv
            integer                                 ::      ii,jj,ierror
            integer                                 ::      ignoreCount
            
            real(kind=real64)       ::      sumw,sumwx      !,sumwy
            real(kind=real64)       ::      sumwxx          !,sumwxy,sumwyy
            real(kind=real64)       ::      sumwxxx         !   ,sumwxxy,sumwxyy,sumwyyy
            real(kind=real64)       ::      sumwxxxx        !   ,sumwxxxy,sumwxxyy,sumwxyyy,sumwyyyy
            real(kind=real64)       ::      sumwf,sumwfx    !    ,sumwfy
            real(kind=real64)       ::      sumwfxx         !,sumwfxy,sumwfyy
            
            real(kind=real64)       ::      wwww,ffff
            
                    real(kind=real64)                           ::      ss
                    real(kind=real64),dimension(-2:2,-2:2)      ::      gg
                    
       
!***********    PROBLEMS WITH DGESV WORKAROUND
            real(kind=real32),dimension(3,3)        ::      AA32
            real(kind=real32),dimension(3)          ::      bb32                                
                    
                    
                    
                    
                    
                    
                    
                    
                    
        !---    set up the weighting on the nodes to use.
            ignoreCount = 0
            do jj = -2,2
                do ii = -2,2
                    if (f(ii,jj)==LOCALFILTER2D_IGNORE) then
                        ww(ii,jj) = 0.0d0
                        ignoreCount = ignoreCount + 1
                    else
                        ww(ii,jj) = w(ii,jj)
                    end if
                end do
            end do
            
            if (dbg) print *,"ignoreCount ",ignoreCount
            if (f(0,0)/= LOCALFILTER2D_IGNORE) then
                if (ignoreCount >= 10) ww(0,0) = 0.5d0
                if (ignoreCount >= 12) ww(0,0) = 1.0d0
            end if
            if (ignoreCount >= 14) then
                sumw     = 0; sumwf    = 0 
                do jj = -2,2
                    do ii = -2,2
                        sumw     = sumw     + ww(ii,jj)
                        sumwf    = sumwf    + ww(ii,jj)*f(ii,jj) 
                    end do
                end do
                if (sumw>0) then
                    f0 = sumwf/sumw
                else
                    f0 = LOCALFILTER2D_IGNORE
                end if
                return
            end if
                 
            
        !---    now need to solve the set of linear equations A x = b for x
        !   find matrix A and rhs b
            sumw     = 0; sumwx    = 0 
            sumwxx   = 0   
            sumwxxx  = 0 
            sumwxxxx = 0 
            sumwf    = 0; sumwfx   = 0 
            sumwfxx  = 0 
            do jj = -2,2
                do ii = -2,2
                    if (horv=="V") then             !   removing vertical stripes means smooth over x values
                        wwww     = ww(ii,jj)
                        ffff     = f(ii,jj)
                    else                            !   removing horizontal stripes means smooth over y values
                        wwww     = ww(jj,ii)
                        ffff     = f(jj,ii)
                    end if
                        
                    sumw     = sumw     + wwww
                    sumwx    = sumwx    + wwww*ii
                    sumwxx   = sumwxx   + wwww*ii*ii
                    sumwxxx  = sumwxxx  + wwww*ii*ii*ii
                    sumwxxxx = sumwxxxx + wwww*ii*ii*ii*ii
                    wwww     = wwww*ffff
                    sumwf    = sumwf    + wwww 
                    sumwfx   = sumwfx   + wwww*ii
                    sumwfxx  = sumwfxx  + wwww*ii*ii
                end do
            end do
            AA(1:3,1) = (/ sumw    ,sumwx    ,sumwxx     /)
            AA(1:3,2) = (/ sumwx   ,sumwxx   ,sumwxxx    /)           
            AA(1:3,3) = (/ sumwxx/2,sumwxxx/2,sumwxxxx/2 /)
            bb(1:3)   = (/ sumwf   ,sumwfx   ,sumwfxx    /)        
                    
            
            
            
            
            if (dbg) then
                print *,"Lib_LocalFilter2d::interpolate55_stripe info " 
                do ii = -2,2
                    write (*,fmt='(5f16.8,a,5f16.8)') f(ii,:)," , ",ww(ii,:)
                end do                
                do ii = 1,3
                    write (*,fmt='(3f16.8,a,f16.8)') AA(ii,:)," , ",bb(ii)
                end do
            end if             
                                
            
!***********        PROBLEMS WITH DGESV
!           WORKAROUND WITH SGESV     
        
#ifndef LAPACK_DSV_WORKAROUND
            call DGESV( 3,1,AA,3,ipiv,bb,3, ierror )
#else
            !
            AA32 = real(AA,kind=real32)
            BB32 = real(BB,kind=real32)
            call SGESV( 3,1,AA32,3,ipiv,bb32,3, ierror ) 
            bb = real(BB32,kind=real64)
#endif            
!           END OF WORKAROUND
!***********
            
            
            if (dbg) then
                print *,"Lib_LocalFilter2d::interpolate55_stripe info - DGESV ",ierror,bb                               
                ss = 0.0d0
                do jj = -2,2
                    do ii = -2,2
                        gg(ii,jj) = (bb(1)+bb(2)*ii+bb(3)*ii*ii/2)
                        ss = ss + ww(ii,jj)*(f(ii,jj)-gg(ii,jj))**2 
                    end do
                end do
                print *,"Lib_LocalFilter2d::interpolate55_stripe info - error ",ss
                 do ii = -2,2               
                     write (*,fmt='(3(5f8.3,a))') f(ii,:)," , ",ww(ii,:)," , ",gg(ii,:)
                 end do
            end if           
            
            
            if (ierror /= 0) then
            !   could not find a good solution. This could be because there aren't enough non-zero coefficients to fit quadratic...
            !   try again, this time just with simple linear gradient solution.
                bb(1)   = LOCALFILTER2D_IGNORE
            end if                       
                    
            f0 = bb(1)        
            !stop
            return
        end function interpolate55_stripe
            
    
    
    !---
    
        subroutine localFilter0( f_in,f_out )
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      perform a local filter on f_in
            real(kind=real64),dimension(:,:),intent(in)         ::      f_in
            real(kind=real64),dimension(:,:),intent(out)        ::      f_out
            
            real(kind=real64)       ::      stdev
            call denoise( f_in,0.0d0,f_out,stdev )
            
            return
        end subroutine localFilter0
        
        subroutine localFilter1( f_in,f_out )
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      perform a local filter on f_in
            real(kind=real32),dimension(:,:),intent(in)         ::      f_in
            real(kind=real32),dimension(:,:),intent(out)        ::      f_out
            
            real(kind=real64)       ::      stdev
            call denoise( f_in,0.0d0,f_out,stdev )
            
            return
        end subroutine localFilter1
        
    !---

    
        recursive subroutine denoise0( f_in,sig,f_out,stdev )
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      find standard deviation of noise, defined as rms difference berween input value and filtered value
    !*      and smooth any dodgy pixels with > sig standard deviations.
            real(kind=real64),dimension(:,:),intent(in)         ::      f_in
            real(kind=real64),intent(in)                        ::      sig
            real(kind=real64),dimension(:,:),intent(out)        ::      f_out
            real(kind=real64),intent(out)                       ::      stdev
            real(kind=real64),dimension(:,:),allocatable        ::      f_tmp
            
            integer         ::      nx,ny
            integer         ::      ix,iy
            real(kind=real64),dimension(-2:2,-2:2)      ::      f_local
            real(kind=real64)           ::      f0
            real(kind=real64)           ::      df2sum,fsum
            integer                     ::      nGoodPix
            
            
            nx = size(f_in,dim=1)
            ny = size(f_in,dim=2)
            
            df2sum = 0.0d0
            fsum = 0.0d0
            nGoodPix = 0
            do iy = 1,ny
                do ix = 1,nx
                    f0 = f_in(ix,iy)
                    if (f0 == LOCALFILTER2D_IGNORE) cycle
                    fsum = fsum + f0
                    df2sum = df2sum + f0*f0
                    nGoodPix = nGoodPix + 1
                end do
            end do
             
            
            if (nGoodPix>1) then
                df2sum = df2sum/nGoodPix
                fsum = fsum/nGoodPix
                stdev = sqrt( df2sum - fsum*fsum )
            else
                stdev = 0                
            end if
            if (dbg) print *,"Lib_LocalFilter2d::denoise0 - info global stdev,ngoodpix ",stdev,ngoodpix 
            
                    
            
            df2sum = 0.0d0
            nGoodPix = 0
            
            do iy = 1,ny
            
                f_local = LOCALFILTER2D_IGNORE
                if (iy == 1) then   
                    f_local(0:2,0:2) = f_in(1:3,1:3)
                else if (iy == 2) then
                    f_local(0:2,-1:2) = f_in(1:3,1:4)
                else if (iy == ny-1) then
                    f_local(0:2,-2:1) = f_in(1:3,ny-3:ny)
                else if (iy == ny) then
                    f_local(0:2,-2:0) = f_in(1:3,ny-2:ny)
                else
                    f_local(0:2,-2:2) = f_in(1:3,iy-2:iy+2)                
                end if
            
                do ix = 1,nx
                
                    dbg = ( ((ix-LIB_LF2D_DBGX)*(ix-LIB_LF2D_DBGX) + (iy-LIB_LF2D_DBGY)*(iy-LIB_LF2D_DBGY) )==0 )
                    if (dbg) then
                        print *,""
                        
                        write (*,fmt='(2i4,1000f8.3)') ix,iy,f_local(-2,:)
                        write (*,fmt='(a8,1000f8.3)') "",f_local(-1,:)
                        write (*,fmt='(a8,1000f8.3)') "",f_local(0,:)
                        write (*,fmt='(a8,1000f8.3)') "",f_local(1,:)
                        write (*,fmt='(a8,1000f8.3)') "",f_local(2,:)
                    end if
                    f0 = interpolate55( f_local )
                    f_out(ix,iy) = f0
                    if ( (f0 /= LOCALFILTER2D_IGNORE).and.(f_in(ix,iy) /= LOCALFILTER2D_IGNORE) ) then
                        f0 = f0 - f_in(ix,iy)
                        df2sum = df2sum + f0*f0
                        nGoodPix = nGoodPix + 1
                    end if
                    if (dbg) then
                        print *,"f_out = ",f_out(ix,iy)," df = ",f0
                        print *,""
                    end if         
                !---    slow lines                               
                     f_local(-2,-2:2) = f_local(-1,-2:2)
                     f_local(-1,-2:2) = f_local( 0,-2:2)
                     f_local( 0,-2:2) = f_local( 1,-2:2)
                     f_local( 1,-2:2) = f_local( 2,-2:2)
                !---    let compiler handle it                
!                    f_local = eoshift( f_local,shift=1,dim=1 )

                    if (ix >= nx-2) then
                        f_local(2,-2:2) = LOCALFILTER2D_IGNORE
                    else 
                        if (iy == 1) then   
                            f_local(2,0:2) = f_in(ix+3,1:3)
                        else if (iy == 2) then
                            f_local(2,-1:2) = f_in(ix+3,1:4)
                        else if (iy == ny-1) then
                            f_local(2,-2:1) = f_in(ix+3,ny-3:ny)
                        else if (iy == ny) then
                            f_local(2,-2:0) = f_in(ix+3,ny-2:ny)
                        else
                            f_local(2,-2:2) = f_in(ix+3,iy-2:iy+2)                
                        end if
                    end if                        
            
                end do
            end do    

            if (nGoodPix>1) then
                stdev = sqrt( df2sum / nGoodPix )
            else
                stdev = 0
                return
            end if
            if (dbg) print *,"Lib_LocalFilter2d::denoise0 - info noise  stdev,ngoodpix,sig ",stdev,ngoodpix,sig
            
            if (sig == 0) return    !   note that if called recursively this is where we return.
            
            
            
        !   set f_tmp to be f_in, but now with any dodgy pixel > sig standard deviations set to IGNORE.
            allocate( f_tmp(nx,ny) )    
            f_tmp = f_in        
            do iy = 1,ny
                do ix = 1,nx             
                    dbg = ( ((ix-LIB_LF2D_DBGX)*(ix-LIB_LF2D_DBGX) + (iy-LIB_LF2D_DBGY)*(iy-LIB_LF2D_DBGY) )==0 )
                    if ( (f_out(ix,iy) /= LOCALFILTER2D_IGNORE).and.(f_in(ix,iy) /= LOCALFILTER2D_IGNORE) ) then
                        f0 = f_out(ix,iy) - f_in(ix,iy)
                        if (f0*f0 > stdev*stdev*sig*sig) then
                            f_tmp(ix,iy) = LOCALFILTER2D_IGNORE
                            if (dbg) print *,"pixel ",ix,iy," input ",f_in(ix,iy)," output ",f_out(ix,iy)," set to ignore, std dev ",f0/stdev
                        end if
                    end if                        
                end do
            end do
            
        !   now recompute the filtered image, this time with the dodgy pixels not contributing to the standard deviation
            call denoise( f_tmp,0.0d0,f_out,stdev )
            deallocate(f_tmp)
            
        !   now I have a good guess of what the filtered image would look like if there were no dodgy pixels,
        !   and hence a good idea of what the true background noise is.
        !   set the pixels < sig standard deviations back to the original values - this is a denoise not a smoothing filter.
        
            nGoodPix = 0
            df2sum = 0.0d0
            do iy = 1,ny
                do ix = 1,nx                    
                    dbg = ( ((ix-LIB_LF2D_DBGX)*(ix-LIB_LF2D_DBGX) + (iy-LIB_LF2D_DBGY)*(iy-LIB_LF2D_DBGY) )==0 )
                    if ( (f_out(ix,iy) /= LOCALFILTER2D_IGNORE).and.(f_in(ix,iy) /= LOCALFILTER2D_IGNORE) ) then
                        f0 = f_out(ix,iy) - f_in(ix,iy)
                        if (f0*f0 < stdev*stdev*sig*sig) then
                            if (dbg) print *,"pixel ",ix,iy," input ",f_in(ix,iy)," output ",f_out(ix,iy)," set to revert, std dev ",f0/stdev
                            f_out(ix,iy) = f_in(ix,iy)
                            f0 = 0
                        else
                           if (dbg)  print *,"pixel ",ix,iy," input ",f_in(ix,iy)," output ",f_out(ix,iy)," set to flatten, std dev ",f0/stdev
                        end if
                        nGoodPix = nGoodPix + 1
                        df2sum = df2sum + f0*f0
                    else  
                        f_out(ix,iy) = f_in(ix,iy)
                    end if                        
                end do
            end do
            
            if (nGoodPix>1) then
                stdev = sqrt( df2sum / nGoodPix )
            else
                stdev = 0
            end if
            if (dbg) print *,"Lib_LocalFilter2d::denoise0 - info noise  stdev,ngoodpix,sig ",stdev,ngoodpix,sig
            return
        end subroutine denoise0        
        
    
        recursive subroutine denoise1( f_in,sig,f_out,stdev )
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      find standard deviation of noise, defined as rms difference berween input value and filtered value
    !*      and smooth any dodgy pixels with > sig standard deviations.
            real(kind=real32),dimension(:,:),intent(in)         ::      f_in
            real(kind=real64),intent(in)                        ::      sig
            real(kind=real32),dimension(:,:),intent(out)        ::      f_out
            real(kind=real64),intent(out)                       ::      stdev
            real(kind=real32),dimension(:,:),allocatable        ::      f_tmp
            
            integer         ::      nx,ny
            integer         ::      ix,iy
            real(kind=real32),dimension(-2:2,-2:2)      ::      f_local
            real(kind=real32)           ::      f0
            real(kind=real32)           ::      df2sum,fsum
            integer                     ::      nGoodPix
            
            
            nx = size(f_in,dim=1)
            ny = size(f_in,dim=2)
            
            df2sum = 0.0d0
            fsum = 0.0d0
            nGoodPix = 0
            do iy = 1,ny
                do ix = 1,nx
                    f0 = f_in(ix,iy)
                    if (f0 == LOCALFILTER2D_IGNORE32) cycle
                    fsum = fsum + f0
                    df2sum = df2sum + f0*f0
                    nGoodPix = nGoodPix + 1
                end do
            end do
             
            
            if (nGoodPix>1) then
                df2sum = df2sum/nGoodPix
                fsum = fsum/nGoodPix
                stdev = sqrt( df2sum - fsum*fsum )
            else
                stdev = 0                
            end if
            if (dbg) print *,"Lib_LocalFilter2d::denoise1 - info global stdev,ngoodpix ",stdev,ngoodpix 
            
                    
            
            df2sum = 0.0d0
            nGoodPix = 0
            
            do iy = 1,ny
            
                f_local = LOCALFILTER2D_IGNORE32
                if (iy == 1) then   
                    f_local(0:2,0:2) = f_in(1:3,1:3)
                else if (iy == 2) then
                    f_local(0:2,-1:2) = f_in(1:3,1:4)
                else if (iy == ny-1) then
                    f_local(0:2,-2:1) = f_in(1:3,ny-3:ny)
                else if (iy == ny) then
                    f_local(0:2,-2:0) = f_in(1:3,ny-2:ny)
                else
                    f_local(0:2,-2:2) = f_in(1:3,iy-2:iy+2)                
                end if
            
                do ix = 1,nx
                
                    dbg = ( ((ix-LIB_LF2D_DBGX)*(ix-LIB_LF2D_DBGX) + (iy-LIB_LF2D_DBGY)*(iy-LIB_LF2D_DBGY) )==0 )
                    if (dbg) then
                        print *,""
                        
                        write (*,fmt='(2i4,1000f8.3)') ix,iy,f_local(-2,:)
                        write (*,fmt='(a8,1000f8.3)') "",f_local(-1,:)
                        write (*,fmt='(a8,1000f8.3)') "",f_local(0,:)
                        write (*,fmt='(a8,1000f8.3)') "",f_local(1,:)
                        write (*,fmt='(a8,1000f8.3)') "",f_local(2,:)
                    end if
                    f0 = interpolate55( f_local )
                    f_out(ix,iy) = f0
                    if ( (f0 /= LOCALFILTER2D_IGNORE32).and.(f_in(ix,iy) /= LOCALFILTER2D_IGNORE32) ) then
                        f0 = f0 - f_in(ix,iy)
                        df2sum = df2sum + f0*f0
                        nGoodPix = nGoodPix + 1
                    end if
                    if (dbg) then
                        print *,"f_out = ",f_out(ix,iy)," df = ",f0
                        print *,""
                    end if         
                !---    slow lines                               
                     f_local(-2,-2:2) = f_local(-1,-2:2)
                     f_local(-1,-2:2) = f_local( 0,-2:2)
                     f_local( 0,-2:2) = f_local( 1,-2:2)
                     f_local( 1,-2:2) = f_local( 2,-2:2)
                !---    let compiler handle it                
!                    f_local = eoshift( f_local,shift=1,dim=1 )

                    if (ix >= nx-2) then
                        f_local(2,-2:2) = LOCALFILTER2D_IGNORE32
                    else 
                        if (iy == 1) then   
                            f_local(2,0:2) = f_in(ix+3,1:3)
                        else if (iy == 2) then
                            f_local(2,-1:2) = f_in(ix+3,1:4)
                        else if (iy == ny-1) then
                            f_local(2,-2:1) = f_in(ix+3,ny-3:ny)
                        else if (iy == ny) then
                            f_local(2,-2:0) = f_in(ix+3,ny-2:ny)
                        else
                            f_local(2,-2:2) = f_in(ix+3,iy-2:iy+2)                
                        end if
                    end if                        
            
                end do
            end do    

            if (nGoodPix>1) then
                stdev = sqrt( df2sum / nGoodPix )
            else
                stdev = 0
                return
            end if
            if (dbg) print *,"Lib_LocalFilter2d::denoise1 - info noise  stdev,ngoodpix,sig ",stdev,ngoodpix,sig
            
            if (sig == 0) return    !   note that if called recursively this is where we return.
            
            
            
        !   set f_tmp to be f_in, but now with any dodgy pixel > sig standard deviations set to IGNORE.
            allocate( f_tmp(nx,ny) )    
            f_tmp = f_in        
            do iy = 1,ny
                do ix = 1,nx             
                    dbg = ( ((ix-LIB_LF2D_DBGX)*(ix-LIB_LF2D_DBGX) + (iy-LIB_LF2D_DBGY)*(iy-LIB_LF2D_DBGY) )==0 )
                    if ( (f_out(ix,iy) /= LOCALFILTER2D_IGNORE32).and.(f_in(ix,iy) /= LOCALFILTER2D_IGNORE32) ) then
                        f0 = f_out(ix,iy) - f_in(ix,iy)
                        if (f0*f0 > stdev*stdev*sig*sig) then
                            f_tmp(ix,iy) = LOCALFILTER2D_IGNORE32
                            if (dbg) print *,"pixel ",ix,iy," input ",f_in(ix,iy)," output ",f_out(ix,iy)," set to ignore, std dev ",f0/stdev
                        end if
                    end if                        
                end do
            end do
            
        !   now recompute the filtered image, this time with the dodgy pixels not contributing to the standard deviation
            call denoise( f_tmp,0.0d0,f_out,stdev )
            deallocate(f_tmp)
            
        !   now I have a good guess of what the filtered image would look like if there were no dodgy pixels,
        !   and hence a good idea of what the true background noise is.
        !   set the pixels < sig standard deviations back to the original values - this is a denoise not a smoothing filter.
        
            nGoodPix = 0
            df2sum = 0.0d0
            do iy = 1,ny
                do ix = 1,nx                    
                    dbg = ( ((ix-LIB_LF2D_DBGX)*(ix-LIB_LF2D_DBGX) + (iy-LIB_LF2D_DBGY)*(iy-LIB_LF2D_DBGY) )==0 )
                    if ( (f_out(ix,iy) /= LOCALFILTER2D_IGNORE32).and.(f_in(ix,iy) /= LOCALFILTER2D_IGNORE32) ) then
                        f0 = f_out(ix,iy) - f_in(ix,iy)
                        if (f0*f0 < stdev*stdev*sig*sig) then
                            if (dbg) print *,"pixel ",ix,iy," input ",f_in(ix,iy)," output ",f_out(ix,iy)," set to revert, std dev ",f0/stdev
                            f_out(ix,iy) = f_in(ix,iy)
                            f0 = 0
                        else
                           if (dbg)  print *,"pixel ",ix,iy," input ",f_in(ix,iy)," output ",f_out(ix,iy)," set to flatten, std dev ",f0/stdev
                        end if
                        nGoodPix = nGoodPix + 1
                        df2sum = df2sum + f0*f0
                    else  
                        f_out(ix,iy) = f_in(ix,iy)
                    end if                        
                end do
            end do
            
            if (nGoodPix>1) then
                stdev = sqrt( df2sum / nGoodPix )
            else
                stdev = 0
            end if
            if (dbg) print *,"Lib_LocalFilter2d::denoise1 - info noise  stdev,ngoodpix,sig ",stdev,ngoodpix,sig
            return
        end subroutine denoise1        
        
    
        recursive subroutine denoise_stripe( f_in,sig,f_out,horv )
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      find standard deviation of noise, defined as rms difference berween input value and filtered value
    !*      and smooth any dodgy pixels with > sig standard deviations.
            real(kind=real64),dimension(:,:),intent(in)         ::      f_in
            real(kind=real64),intent(in)                        ::      sig
            real(kind=real64),dimension(:,:),intent(out)        ::      f_out
            character(len=1),intent(in)                         ::      horv
            real(kind=real64),dimension(:,:),allocatable        ::      f_tmp
            
            integer         ::      nx,ny
            integer         ::      ix,iy
            real(kind=real64),dimension(-2:2,-2:2)      ::      f_local
            real(kind=real64)           ::      f0
            real(kind=real64)           ::      df2sum,fsum,stdev
            integer                     ::      nGoodPix
            
            dbg = LIB_LF2D_DBGX /= -1
            
            
            nx = size(f_in,dim=1)
            ny = size(f_in,dim=2)
            
            df2sum = 0.0d0
            fsum = 0.0d0
            nGoodPix = 0
            do iy = 1,ny
                do ix = 1,nx
                    f0 = f_in(ix,iy)
                    if (f0 == LOCALFILTER2D_IGNORE) cycle
                    fsum = fsum + f0
                    df2sum = df2sum + f0*f0
                    nGoodPix = nGoodPix + 1
                end do
            end do
             
            
            if (nGoodPix>1) then
                df2sum = df2sum/nGoodPix
                fsum = fsum/nGoodPix
                stdev = sqrt( df2sum - fsum*fsum )
            else
                stdev = 0                
            end if
            if (dbg) print *,"Lib_LocalFilter2d::denoise_stripe - info global stdev,ngoodpix ",stdev,ngoodpix 
            
                    
            
            df2sum = 0.0d0
            nGoodPix = 0
            
            do iy = 1,ny
            
                f_local = LOCALFILTER2D_IGNORE
                if (iy == 1) then   
                    f_local(0:2,0:2) = f_in(1:3,1:3)
                else if (iy == 2) then
                    f_local(0:2,-1:2) = f_in(1:3,1:4)
                else if (iy == ny-1) then
                    f_local(0:2,-2:1) = f_in(1:3,ny-3:ny)
                else if (iy == ny) then
                    f_local(0:2,-2:0) = f_in(1:3,ny-2:ny)
                else
                    f_local(0:2,-2:2) = f_in(1:3,iy-2:iy+2)                
                end if
            
                do ix = 1,nx
                
                    dbg = ( ((ix-LIB_LF2D_DBGX)*(ix-LIB_LF2D_DBGX) + (iy-LIB_LF2D_DBGY)*(iy-LIB_LF2D_DBGY) )==0 )
                    if (dbg) then
                        print *,""
                        
                        write (*,fmt='(2i4,1000f8.3)') ix,iy,f_local(-2,:)
                        write (*,fmt='(a8,1000f8.3)') "",f_local(-1,:)
                        write (*,fmt='(a8,1000f8.3)') "",f_local(0,:)
                        write (*,fmt='(a8,1000f8.3)') "",f_local(1,:)
                        write (*,fmt='(a8,1000f8.3)') "",f_local(2,:)
                    end if
                    f0 = interpolate55( f_local,horv )
                    f_out(ix,iy) = f0
                    if ( (f0 /= LOCALFILTER2D_IGNORE).and.(f_in(ix,iy) /= LOCALFILTER2D_IGNORE) ) then
                        f0 = f0 - f_in(ix,iy)
                        df2sum = df2sum + f0*f0
                        nGoodPix = nGoodPix + 1
                    end if
                    if (dbg) then
                        print *,"f_out = ",f_out(ix,iy)," df = ",f0
                        print *,""
                    end if         
                !---    slow lines                               
                     f_local(-2,-2:2) = f_local(-1,-2:2)
                     f_local(-1,-2:2) = f_local( 0,-2:2)
                     f_local( 0,-2:2) = f_local( 1,-2:2)
                     f_local( 1,-2:2) = f_local( 2,-2:2)
                !---    let compiler handle it                
!                    f_local = eoshift( f_local,shift=1,dim=1 )

                    if (ix >= nx-2) then
                        f_local(2,-2:2) = LOCALFILTER2D_IGNORE
                    else 
                        if (iy == 1) then   
                            f_local(2,0:2) = f_in(ix+3,1:3)
                        else if (iy == 2) then
                            f_local(2,-1:2) = f_in(ix+3,1:4)
                        else if (iy == ny-1) then
                            f_local(2,-2:1) = f_in(ix+3,ny-3:ny)
                        else if (iy == ny) then
                            f_local(2,-2:0) = f_in(ix+3,ny-2:ny)
                        else
                            f_local(2,-2:2) = f_in(ix+3,iy-2:iy+2)                
                        end if
                    end if                        
            
                end do
            end do    

            if (nGoodPix>1) then
                stdev = sqrt( df2sum / nGoodPix )
            else
                stdev = 0
                return
            end if
            if (dbg) print *,"Lib_LocalFilter2d::denoise0 - info noise  stdev,ngoodpix,sig ",stdev,ngoodpix,sig
            
            if (sig == 0) return    !   note that if called recursively this is where we return.
            
            
            
        !   set f_tmp to be f_in, but now with any dodgy pixel > sig standard deviations set to IGNORE.
            allocate( f_tmp(nx,ny) )    
            f_tmp = f_in        
            do iy = 1,ny
                do ix = 1,nx             
                    dbg = ( ((ix-LIB_LF2D_DBGX)*(ix-LIB_LF2D_DBGX) + (iy-LIB_LF2D_DBGY)*(iy-LIB_LF2D_DBGY) )==0 )
                    if ( (f_out(ix,iy) /= LOCALFILTER2D_IGNORE).and.(f_in(ix,iy) /= LOCALFILTER2D_IGNORE) ) then
                        f0 = f_out(ix,iy) - f_in(ix,iy)
                        if (f0*f0 > stdev*stdev*sig*sig) then
                            f_tmp(ix,iy) = LOCALFILTER2D_IGNORE
                            if (dbg) print *,"pixel ",ix,iy," input ",f_in(ix,iy)," output ",f_out(ix,iy)," set to ignore, std dev ",f0/stdev
                        end if
                    end if                        
                end do
            end do
            
        !   now recompute the filtered image, this time with the dodgy pixels not contributing to the standard deviation
            call denoise( f_tmp,0.0d0,f_out,horv )
            deallocate(f_tmp)
            
        !   now I have a good guess of what the filtered image would look like if there were no dodgy pixels,
        !   and hence a good idea of what the true background noise is.
        !   set the pixels < sig standard deviations back to the original values - this is a denoise not a smoothing filter.
        
            nGoodPix = 0
            df2sum = 0.0d0
            do iy = 1,ny
                do ix = 1,nx                    
                    dbg = ( ((ix-LIB_LF2D_DBGX)*(ix-LIB_LF2D_DBGX) + (iy-LIB_LF2D_DBGY)*(iy-LIB_LF2D_DBGY) )==0 )
                    if ( (f_out(ix,iy) /= LOCALFILTER2D_IGNORE).and.(f_in(ix,iy) /= LOCALFILTER2D_IGNORE) ) then
                        f0 = f_out(ix,iy) - f_in(ix,iy)
                        if (f0*f0 < stdev*stdev*sig*sig) then
                            if (dbg) print *,"pixel ",ix,iy," input ",f_in(ix,iy)," output ",f_out(ix,iy)," set to revert, std dev ",f0/stdev
                            f_out(ix,iy) = f_in(ix,iy)
                            f0 = 0
                        else
                           if (dbg)  print *,"pixel ",ix,iy," input ",f_in(ix,iy)," output ",f_out(ix,iy)," set to flatten, std dev ",f0/stdev
                        end if
                        nGoodPix = nGoodPix + 1
                        df2sum = df2sum + f0*f0
                    else  
                        f_out(ix,iy) = f_in(ix,iy)
                    end if                        
                end do
            end do
            
            if (nGoodPix>1) then
                stdev = sqrt( df2sum / nGoodPix )
            else
                stdev = 0
            end if
            if (dbg) print *,"Lib_LocalFilter2d::denoise0 - info noise  stdev,ngoodpix,sig ",stdev,ngoodpix,sig
            return
        end subroutine denoise_stripe     
        
    !---
    
    
        
        subroutine halfSizeImage( img_in,img_out )
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      shrink image by 50% . note img_out must have dimensions ceiling(mx*0.5)                
            real(kind=real64),dimension(:,:),intent(in)       ::      img_in
            real(kind=real64),dimension(:,:),intent(inout)    ::      img_out
            integer             ::      ix,iy , mx,my
            real(kind=real64)   ::      gbar
            integer             ::      nn
            
            
            mx = size(img_in,dim=1)
            my = size(img_in,dim=2)
            img_out = 0.0d0
            
            do iy = 1,my-1,2
                do ix = 1,mx-1,2
                
                    nn = 0 ; gbar = 0.0d0
                    if (img_in(ix,iy) /= LOCALFILTER2D_IGNORE ) then
                        gbar = gbar + img_in(ix,iy)
                        nn = nn + 1
                    end if
                    if (img_in(ix,iy+1) /= LOCALFILTER2D_IGNORE ) then
                        gbar = gbar + img_in(ix,iy+1)
                        nn = nn + 1
                    end if
                    if (img_in(ix+1,iy) /= LOCALFILTER2D_IGNORE ) then
                        gbar = gbar + img_in(ix+1,iy)
                        nn = nn + 1
                    end if
                    if (img_in(ix+1,iy+1) /= LOCALFILTER2D_IGNORE ) then
                        gbar = gbar + img_in(ix+1,iy+1)
                        nn = nn + 1
                    end if
                    if (nn == 0) then   
                        gbar = LOCALFILTER2D_IGNORE
                    else
                        gbar = gbar/nn
                    end if
                    
                    img_out( ix/2+1,iy/2+1 ) = gbar                     !   note: gfortran issues a Warning: Integer division truncated to constant ‘0’ at (1) [-Winteger-division]
                                                                        !   while perfectly true I'm using integer division, when ix>1, it doesn't equal zero.
                    
                end do
                
                if (mod(mx,2)==1) then
                    
                    nn = 0 ; gbar = 0.0d0
                    if (img_in(mx,iy) /= LOCALFILTER2D_IGNORE ) then
                        gbar = gbar + img_in(mx,iy)
                        nn = nn + 1
                    end if
                    if (img_in(mx,iy+1) /= LOCALFILTER2D_IGNORE ) then
                        gbar = gbar + img_in(mx,iy+1)
                        nn = nn + 1
                    end if
                    if (nn == 0) then   
                        gbar = LOCALFILTER2D_IGNORE
                    else
                        gbar = gbar/nn
                    end if
                    
                    img_out( mx/2+1,iy/2+1 ) = gbar             
                end if
                
            end do          
            if (mod(my,2)==1) then
                do ix = 1,mx-1,2
                
                    nn = 0 ; gbar = 0.0d0
                    if (img_in(ix,my) /= LOCALFILTER2D_IGNORE ) then
                        gbar = gbar + img_in(ix,my)
                        nn = nn + 1
                    end if                                
                    if (img_in(ix+1,my) /= LOCALFILTER2D_IGNORE ) then
                        gbar = gbar + img_in(ix+1,my)
                        nn = nn + 1                                       
                    end if
                    if (nn == 0) then   
                        gbar = LOCALFILTER2D_IGNORE
                    else
                        gbar = gbar/nn
                    end if
                    
                    img_out( ix/2+1,my/2+1 ) = gbar                   
                    
                end do
                
                if (mod(mx,2)==1) then                    
                    gbar = LOCALFILTER2D_IGNORE
                    if (img_in(mx,my) /= LOCALFILTER2D_IGNORE ) gbar = img_in(mx,my)
                    img_out( mx/2+1,my/2+1 ) = gbar             
                end if
            
            end if
            
                  

            return
        end subroutine halfSizeImage
            
            
    !---    
        
        subroutine doubleSizeImage( img_in,img_out )
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      enlarge image by 100% . note img_in must have dimensions ceiling(mx*0.5)
            real(kind=real64),dimension(0:,0:),intent(in)           ::      img_in
            real(kind=real64),dimension(0:,0:),intent(inout)        ::      img_out
            
            
            integer         ::      Mx,My
            integer         ::      ix,iy
            
            Mx = size(img_out,dim=1)
            My = size(img_out,dim=2)
          
             
        !---    copy across the even number points, fill in the odd
            do iy = 0,My-1,2
                do ix = 0,Mx-1,2        !   mx=4; ix=0,2   mx=5; ix=0,2,4
                    img_out( ix,iy ) = img_in( ix/2,iy/2 )
                end do
                if (mod(mx,2)==0) then
                    call fillInOdd( img_out(0:mx-2,iy) )
                    img_out( mx-1,iy ) = img_in( mx/2+1,iy/2+1 )
                else
                    call fillInOdd( img_out(:,iy) )
                end if
            end do 
            do ix = 0,2*Mx-2
                call fillInOdd( img_out(ix,:) )
            end do            
            
            return
            
        contains
    !---^^^^^^^^
        
    
            subroutine fillInOdd( line )
        !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^ 
        !*      use a quadratic interpolation to fill in blanks along a line
                real(kind=real64),dimension(0:),intent(inout)       ::   line                
                integer         ::      mm,ii
                
                integer                         ::   indx
                real(kind=real64),dimension(4)  ::   yy
                
                mm = size(line)
                if (mm == 3) then
                    indx = 0 
                    if (line(0) /= LOCALFILTER2D_IGNORE) indx = indx + 1
                    if (line(2) /= LOCALFILTER2D_IGNORE) indx = indx + 2
                    yy(1:2) = (/ line(0),line(2) /)
                    line(1) = interpolateUsingIndx(yy,indx) 
                else
                
                    indx = 0 
                    if (line(0) /= LOCALFILTER2D_IGNORE) indx = indx + 2
                    if (line(2) /= LOCALFILTER2D_IGNORE) indx = indx + 4
                    if (line(4) /= LOCALFILTER2D_IGNORE) indx = indx + 8
                    yy(2:4) = (/ line(0),line(2),line(4) /)
                    line(1) = interpolateUsingIndx(yy,indx) 
                    
                    do ii = 3,mm-4,2
                        indx = 0                    
                        if (line(ii-3) /= LOCALFILTER2D_IGNORE) indx = indx + 1
                        if (line(ii-1) /= LOCALFILTER2D_IGNORE) indx = indx + 2
                        if (line(ii+1) /= LOCALFILTER2D_IGNORE) indx = indx + 4
                        if (line(ii+3) /= LOCALFILTER2D_IGNORE) indx = indx + 8
                        yy(1:4) = (/ line(ii-3),line(ii-1),line(ii+1),line(ii+3) /)
                        line(ii) = interpolateUsingIndx(yy,indx) 
                    end do
                    
                    indx = 0 
                    if (line(mm-5) /= LOCALFILTER2D_IGNORE) indx = indx + 1
                    if (line(mm-3) /= LOCALFILTER2D_IGNORE) indx = indx + 2
                    if (line(mm-1) /= LOCALFILTER2D_IGNORE) indx = indx + 4
                    yy(1:3) = (/ line(mm-5),line(mm-3),line(mm-1) /)
                    line(mm-2) = interpolateUsingIndx(yy,indx) 
                    
                end if
                
                return
            end subroutine fillInOdd
                       
                                
        end subroutine doubleSizeImage
            

        pure function interpolateUsingIndx( y,indx ) result(x)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !       find the interpolated value between y(2) and y(3)
    !       when some of the entries of y are unknown.
    !       note that this should only be used for very rate IGNOREd pixels
    !       and could cause problems if there are larger patches.
    
            real(kind=real64),dimension(4),intent(in)       ::      y
            integer,intent(in)                              ::      indx
            real(kind=real64)                               ::      x
            
            real(kind=real64),dimension(4,0:15),parameter   ::      kernel = reshape(           &
                                                           (/       0, 0, 0, 0,                 &            !      0
                                                                   16, 0, 0, 0,                 &            !      1
                                                                    0,16, 0, 0,                 &            !      2
                                                                   -8,24, 0, 0,                 &            !      3
                                                                    0, 0,16, 0,                 &            !      4
                                                                    4, 0,12, 0,                 &            !      5
                                                                    0, 8, 8, 0,                 &            !      6
                                                                   -2,12, 6, 0,                 &            !      7
                                                                    0, 0, 0,16,                 &            !      8
                                                                    8, 0, 0, 8,                 &            !      9
                                                                    0,12, 0, 4,                 &            !      10
                                                                   -4,18, 0, 2,                 &            !      11
                                                                    0, 0,24,-8,                 &            !      12
                                                                    2, 0,18,-4,                 &            !      13
                                                                    0, 6,12,-2,                 &            !      14
                                                                   -1, 9, 9,-1     /) , (/4,16/)  )/16.0d0             
            x = kernel(1,indx)*y(1) + kernel(2,indx)*y(2) + kernel(3,indx)*y(3) + kernel(4,indx)*y(4)                                              
            return
        end function interpolateUsingIndx                           
                
        
        
        
    end module Lib_LocalFilter2d       

        