
    module SolveTripleJunction
!---^^^^^^^^^^^^^^^^^^^^^^^^^^
        use iso_fortran_env
        implicit none
        private
        
        public      ::      SolveTripleJunctionLine
        public      ::      SolveTripleJunctionFace
        
        
        logical,public      ::      dbg_SolveTripleJunction = .false.
        
    contains
!---^^^^^^^^ 
        
    
    
        subroutine SolveTripleJunctionLine( f,g,h , n,x )
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      Given three sets of values on the corners of a cube x,y,z = (-1:1)
    !*      Find faces where the triple junction line passes through
    !*      return the line(s) between the faces.
    !*      return n = 0 if no lines are found
    
        !*          8-------7
        !*         /|      /|
        !*        / |     / |       z
        !*       5-------6  |       | y
        !*       |  4----|--3       |/
        !*       | /     | /        *---x
        !*       |/      |/
        !*       1-------2
        !*   
    
            real(kind=real64),dimension(8),intent(in)       ::      f,g,h
            real(kind=real64),dimension(:,:,:),intent(out)  ::      x       !   (1:3,1:2,1:n)
            integer,intent(out)                             ::      n
            
            real(kind=real64),dimension(4)      ::      fface,gface,hface
            real(kind=real64)                   ::      xface,yface
            
            logical                             ::      ok
            real(kind=real64),dimension(3,6)    ::      xx                                                                                
            integer,dimension(4,6)      ::      face = reshape( (/ 1,2,3,4,             &       !   note: direction of indices points into cube.
                                                                   1,5,6,2,             &
                                                                   1,4,8,5,             &
                                                                   2,6,7,3,             &
                                                                   4,3,7,8,             &
                                                                   5,8,7,6  /),(/4,6/) )
            integer             ::      ii,jj,nn                                                       
             
           
            xx = 0 ; nn = 0
            do ii = 1,6 
                do jj = 1,4
                    fface(jj) = f( face(jj,ii) )
                    gface(jj) = g( face(jj,ii) )
                    hface(jj) = h( face(jj,ii) )
                end do
                if (dbg_SolveTripleJunction) print *,"face ",ii," corners ",face(:,ii)
                
                call SolveTripleJunctionFace( fface,gface,hface, xface,yface, ok )
                if (ok) then
                    nn = nn + 1
                    select case(ii)
                        case(1)
                            xx(1:3,nn) = (/ xface,yface,0.0d0 /)                                                              
                        case(2)                                                               !*          8-------7            
                            xx(1:3,nn) = (/ yface,0.0d0,xface /)                              !*         /|      /|            
                        case(3)                                                               !*        / |     / |       z    
                            xx(1:3,nn) = (/ 0.0d0,xface,yface /)                              !*       5-------6  |       | y  
                        case(4)                                                               !*       |  4----|--3       |/   
                            xx(1:3,nn) = (/ 1.0d0,yface,xface /)                              !*       | /     | /        *---x
                        case(5)                                                               !*       |/      |/              
                            xx(1:3,nn) = (/ xface,1.0d0,yface /)                              !*       1-------2               
                        case(6)                             
                            xx(1:3,nn) = (/ yface,xface,1.0d0 /)
                    end select 
                end if
                
            end do
            
            if (nn == 2) then
                x(1:3,1:2,1) = xx(1:3,1:2)
                n = 1
            else if (nn == 3) then
                x(1:3,1,1:3) = xx(1:3,1:3)
                x(1:3,2,1:3) = 0.50d0
                n = 3
       !     else if (nn == 4) then
       !         x(1:3,1:2,1) = xx(1:3,1:2)
       !         x(1:3,1:2,2) = xx(1:3,2:3)
       !         x(1:3,1:2,3) = xx(1:3,3:4)
       !         x(1:3,1,4)   = xx(1:3,4)
       !         x(1:3,2,4)   = xx(1:3,1)
       !         n = 4
       !     else if (nn >= 5) then
       !         x(1:3,1,1:nn) = xx(1:3,1:nn)
       !         x(1:3,2,1:nn) = 0.5d0
       !         n = nn 
            else if (nn==1) then
                x(1:3,1,1) = xx(1:3,1)
                x(1:3,2,1) = 0.5d0
                n = 1         
            else if (nn<1) then
                n = 0
            else 
                if (dbg_SolveTripleJunction) print *,"SolveTripleJunctionLine error - found ",nn," faces with the line ??"
                n = 0
                !stop
            end if
                                
            
            
            
            return
        end subroutine SolveTripleJunctionLine
    
    
        subroutine SolveTripleJunctionFace( f,g,h, x,y, ok )
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      Given three sets of values on the corners of a square x,y = [0:1)
    !*      Find my triple junction point.
            real(kind=real64),dimension(4),intent(in)   ::      f,g,h
            real(kind=real64),intent(out)               ::      x,y
            logical,intent(out)                         ::      ok
             
            integer             ::      code
            
            real(kind=real64),dimension(-1:1,-1:1)  ::      aa
            real(kind=real64),dimension(2)      ::      da,xx
            real(kind=real64),dimension(2,2)    ::      d2a
            real(kind=real64)                   ::      a0,amin , xms,yms
            real(kind=real64),parameter     ::      delta = 1/6.0d0
            
            
            real(kind=real64),dimension(2)      ::      u1,v1,u2,v2
            real(kind=real64),dimension(2,6)    ::      uu,vv
            integer                             ::      nn,mm
            
            real(kind=real64),dimension(4)      ::      phif,phig,phih
            
         !   call marchingSquaresLine( f-0.5d0, nn,u1,v1,u2,v2 )
         !   if (nn>=1) then
         !       uu(1:2,1) = u1(1:2) ; vv(1:2,1) = v1(1:2) - u1(1:2)
         !   end if
         !   if (nn==2) then                                                                                                           
         !       uu(1:2,2) = u2(1:2) ; vv(1:2,2) = v2(1:2) - u2(1:2)
         !   end if            
         !   
         !   call marchingSquaresLine( g-0.5d0, mm,u1,v1,u2,v2 )
         !   if (mm>=1) then
         !       uu(1:2,nn+1) = u1(1:2) ; vv(1:2,nn+1) = v1(1:2) - u1(1:2) ; nn = nn + 1
         !   end if
         !   if (mm==2) then                                                                                                           
         !       uu(1:2,nn+1) = u2(1:2) ; vv(1:2,nn+1) = v2(1:2) - u2(1:2) ; nn = nn + 1
         !   end if
         !   
         !   
         !   call marchingSquaresLine( h-0.5d0, mm,u1,v1,u2,v2 )
         !   if (mm>=1) then
         !       uu(1:2,nn+1) = u1(1:2) ; vv(1:2,nn+1) = v1(1:2) - u1(1:2) ; nn = nn + 1
         !   end if                                                                     
         !   if (mm==2) then                                                                                                                         
         !       uu(1:2,nn+1) = u2(1:2) ; vv(1:2,nn+1) = v2(1:2) - u2(1:2) ; nn = nn + 1
         !   end if
         !   
         !   call closestPointLines( uu(1:2,1:nn),vv(1:2,1:nn), x,y )
         !   
         !   
         !   print *,"closestPointLines ",x,y
         !   
         !   ok = (x*(1-x)>=0).and.(y*(1-y)>=0)
         !   return
            
          !  phif = min(f,-h)
          !  phig = min(-f,g)
          !  phih = min(-g,h)
          
            phif = f
            phig = g
            phih = h
            
            ok = .false.
           
            if (minval( (/maxval(phif),maxval(phig),maxval(phih)/) )<0.0) then
                !if (dbg_SolveTripleJunction) print *,"nodes not in grains ",maxval(phif),maxval(phig),maxval(phih)
                return      !   one of the nodes is not in a grain proper
            end if
            
            code =   whichVertexGreatest( phif(1),phig(1),phih(1) )    &
                 + 4*whichVertexGreatest( phif(2),phig(2),phih(2) )    &
                 +16*whichVertexGreatest( phif(3),phig(3),phih(3) )    &
                 +64*whichVertexGreatest( phif(4),phig(4),phih(4) )
                        
            
            if (dbg_SolveTripleJunction) then
                
                    write(*,fmt='(a,4f10.5)') "f ",f
                    write(*,fmt='(a,4f10.5)') "g ",g
                    write(*,fmt='(a,4f10.5)') "h ",h
                    write(*,fmt='(a,4f10.5)') "phif ",phif
                    write(*,fmt='(a,4f10.5)') "phig ",phig
                    write(*,fmt='(a,4f10.5)') "phih ",phih
                    print *,"code ",code
            end if     
                 
                 
            ok = .true. ; x = 0.5d0 ; y = 0.5d0
            select case (code)
            
            !   ffgh
                case ( 0 + 4*0 + 16*1 + 64*2 ) ; y = 1/3.
                case ( 0 + 4*0 + 16*2 + 64*1 ) ; y = 1/3.         
                case ( 0 + 4*1 + 16*0 + 64*2 )  
                case ( 0 + 4*2 + 16*0 + 64*1 )           
                case ( 0 + 4*1 + 16*2 + 64*0 ) ; x = 1/3.
                case ( 0 + 4*2 + 16*1 + 64*0 ) ; x = 1/3.          
                case ( 1 + 4*0 + 16*0 + 64*2 ) ; x = 2/3.
                case ( 2 + 4*0 + 16*0 + 64*1 ) ; x = 2/3.          
                case ( 1 + 4*0 + 16*2 + 64*0 )  
                case ( 2 + 4*0 + 16*1 + 64*0 )           
                case ( 1 + 4*2 + 16*0 + 64*0 ) ; y = 2/3.
                case ( 2 + 4*1 + 16*0 + 64*0 ) ; y = 2/3.
             
            !   ggfh
                case ( 1 + 4*1 + 16*2 + 64*0 ) ; y = 1/3. 
                case ( 1 + 4*1 + 16*0 + 64*2 ) ; y = 1/3. 
                case ( 1 + 4*2 + 16*1 + 64*0 )            
                case ( 1 + 4*0 + 16*1 + 64*2 )            
                case ( 1 + 4*2 + 16*0 + 64*1 ) ; x = 1/3. 
                case ( 1 + 4*0 + 16*2 + 64*1 ) ; x = 1/3. 
                case ( 2 + 4*1 + 16*1 + 64*0 ) ; x = 2/3. 
                case ( 0 + 4*1 + 16*1 + 64*2 ) ; x = 2/3. 
                case ( 2 + 4*1 + 16*0 + 64*1 )            
                case ( 0 + 4*1 + 16*2 + 64*1 )            
                case ( 2 + 4*0 + 16*1 + 64*1 ) ; y = 2/3. 
                case ( 0 + 4*2 + 16*1 + 64*1 ) ; y = 2/3. 
                
            !   hhfg
                case ( 2 + 4*2 + 16*0 + 64*1 ) ; y = 1/3. 
                case ( 2 + 4*2 + 16*1 + 64*0 ) ; y = 1/3. 
                case ( 2 + 4*0 + 16*2 + 64*1 )            
                case ( 2 + 4*1 + 16*2 + 64*0 )            
                case ( 2 + 4*0 + 16*1 + 64*2 ) ; x = 1/3. 
                case ( 2 + 4*1 + 16*0 + 64*2 ) ; x = 1/3. 
                case ( 0 + 4*2 + 16*2 + 64*1 ) ; x = 2/3. 
                case ( 1 + 4*2 + 16*2 + 64*0 ) ; x = 2/3. 
                case ( 0 + 4*2 + 16*1 + 64*2 )            
                case ( 1 + 4*2 + 16*0 + 64*2 )            
                case ( 0 + 4*1 + 16*2 + 64*2 ) ; y = 2/3. 
                case ( 1 + 4*0 + 16*2 + 64*2 ) ; y = 2/3. 
                   
                
                
                case default
                    ok = .false.
                    x = - 1; y = -1    
                    if (dbg_SolveTripleJunction) print *,"bad code"
                    return      
            end select
            
            !return
            uu = 0 ; vv = 0
            call marchingSquaresLine( f, nn,u1,v1,u2,v2 )
            if (nn>=1) then
                uu(1:2,1) = u1(1:2) ; vv(1:2,1) = v1(1:2) - u1(1:2)
            end if
            if (nn==2) then                                                                                                           
                uu(1:2,2) = u2(1:2) ; vv(1:2,2) = v2(1:2) - u2(1:2)
            end if            
            
            call marchingSquaresLine( g, mm,u1,v1,u2,v2 )
            if (mm>=1) then
                uu(1:2,nn+1) = u1(1:2) ; vv(1:2,nn+1) = v1(1:2) - u1(1:2) ; nn = nn + 1
            end if
            if (mm==2) then                                                                                                           
                uu(1:2,nn+1) = u2(1:2) ; vv(1:2,nn+1) = v2(1:2) - u2(1:2) ; nn = nn + 1
            end if
            
            
            call marchingSquaresLine( h, mm,u1,v1,u2,v2 )
            if (mm>=1) then
                uu(1:2,nn+1) = u1(1:2) ; vv(1:2,nn+1) = v1(1:2) - u1(1:2) ; nn = nn + 1
            end if                                                                     
            if (mm==2) then                                                                                                                         
                uu(1:2,nn+1) = u2(1:2) ; vv(1:2,nn+1) = v2(1:2) - u2(1:2) ; nn = nn + 1
            end if
            
            call closestPointLines( uu(1:2,1:nn),vv(1:2,1:nn), xms,yms )
            
            
            !print *,"closestPointLines ",x,y
            
            
             ok = .true.
            ! ok = (xms*(1-xms)>=0).and.(yms*(1-yms)>=0)
            if ((xms*(1-xms)>=0).and.(yms*(1-yms)>=0)) then
                x = xms
                y = yms
            else
            !    print *,"Number of lines ",nn
            !    do mm = 1,nn
            !        print *,uu(1:2,mm),":",uu(1:2,mm) + vv(1:2,mm)
            !    end do
            !    print *,xms,yms
            !
            !    print *,"f ",f
            !    print *,"g ",g
            !    print *,"h ",h
            !    
            !    
            !    print *,"phif ",phif
            !    print *,"phig ",phig
            !    print *,"phih ",phih
            !    
            !    
            !    
            !
            !    call marchingSquaresLine( f, mm,u1,v1,u2,v2 ) ; print *,"f ",mm,u1,v1,u2,v2               
            !    call marchingSquaresLine( g, mm,u1,v1,u2,v2 ) ; print *,"g ",mm,u1,v1,u2,v2
            !    call marchingSquaresLine( h, mm,u1,v1,u2,v2 ) ; print *,"h ",mm,u1,v1,u2,v2
            !    !stop
            end if
            
            if (dbg_SolveTripleJunction) then
                print *,"Number of lines ",nn
                do mm = 1,nn
                   print *,uu(1:2,mm),":",uu(1:2,mm) + vv(1:2,mm)
                end do
                print *,xms,yms
            end if      
            
            return
            
            
            
           ! print *,"fixed geom ",x,y
            return
            
            
        !---    try to improve the solution
            aa(-1,-1) = distance( f,g,h , -delta,-delta )
            aa( 0,-1) = distance( f,g,h ,  0.0d0,-delta )
            aa( 1,-1) = distance( f,g,h ,  delta,-delta )
            aa(-1, 0) = distance( f,g,h , -delta, 0.0d0 )
            aa( 0, 0) = distance( f,g,h ,  0.0d0, 0.0d0 )
            aa( 1, 0) = distance( f,g,h ,  delta, 0.0d0 )
            aa(-1, 1) = distance( f,g,h , -delta, delta )
            aa( 0, 1) = distance( f,g,h ,  0.0d0, delta )
            aa( 1, 1) = distance( f,g,h ,  delta, delta )
            call stencil( aa, a0,da,d2a , xx, amin )
            xx = xx * delta
            if ( ( (x+xx(1))*(1 - x-xx(1)) > 0 ).and.( (y+xx(2))*(1 - y-xx(2)) > 0 ) ) then
                x = x + xx(1)
                y = y + xx(2)
            end if
            
           ! print *,"iterated geom ",x,y
            return
            
        contains
    !---^^^^^^^^
    
            pure function whichVertexGreatest( f,g,h ) result( v )
        !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
        !*      return 0 if f > max(g,h), 1 for g, 2 for h
                real(kind=real64),intent(in)        ::      f,g,h
                integer                             ::      v
                v = 0                   !   assume f**
                if (g>f) then           !   g** or *gf
                    v = 1
                    if (h>g) v = 2           !   hgf
                else if (h>f) then
                    v = 2
                    if (g>h) v = 1
                end if
                return
            end function whichVertexGreatest
            
            pure function distance( f,g,h , x,y ) result( d )
        !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
        !*      return squared distance from triple j point
                real(kind=real64),dimension(4),intent(in)       ::      f,g,h
                real(kind=real64),intent(in)                    ::      x,y
                real(kind=real64)                               ::      d
                real(kind=real64)       ::      dd
                
                dd = f(1)*(1-x)*(1-y) + f(2)*(x)*(1-y) + f(3)*(x)*(y) + f(4)*(1-x)*(y)  
                d = dd*dd
                dd = g(1)*(1-x)*(1-y) + g(2)*(x)*(1-y) + g(3)*(x)*(y) + g(4)*(1-x)*(y)  
                d = d + dd*dd                                            
                dd = h(1)*(1-x)*(1-y) + h(2)*(x)*(1-y) + h(3)*(x)*(y) + h(4)*(1-x)*(y)  
                d = d + dd*dd
                
                return
            end function distance
            
        end subroutine SolveTripleJunctionFace
    
        
        
        subroutine stencil( f, f0,df,d2f , x, fmin )
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      compute solution of 9 point stencil, find derivatives and minimum point
            real(kind=real64),dimension(-1:1,-1:1),intent(In)       ::      f
            real(kind=real64),intent(out)                           ::      f0
            real(kind=real64),dimension(2),intent(out)              ::      df
            real(kind=real64),dimension(2,2),intent(out)            ::      d2f
            real(kind=real64),dimension(2),intent(out)              ::      x
            real(kind=real64),intent(out)                           ::      fmin
            
            integer             ::      ix,iy
            real(kind=real64)   ::      dd
            f0 = ( 5*f(0,0) + 2*(f(0,1)+f(0,-1)+f(-1,0)+f(1,0)) - (f(1,1)+f(-1,1)+f(1,-1)+f(-1,-1)) )/9

            df(1) = ( (f(1,0)+f(1,1)+f(1,-1)) - (f(-1,0)+f(-1,1)+f(-1,-1)) )/6
            df(2) = ( (f(0,1)+f(1,1)+f(-1,1)) - (f(0,-1)+f(1,-1)+f(-1,-1)) )/6
            
            d2f(1,1) = ( (f(-1,0)+f(-1,-1)+f(-1,1)+f(1,0)+f(1,-1)+f(1,1)) - 2*(f(0,0)+f(0,-1)+f(0,1)) )/3
            d2f(2,2) = ( (f(0,-1)+f(0,1)+f(-1,-1)+f(-1,1)+f(1,-1)+f(1,1)) - 2*(f(0,0)+f(-1,0)+f(1,0)) )/3
            d2f(2,1) = ( (f(-1,-1)+f(1,1)) - (f(-1,1)+f(1,-1)) )/4
            d2f(1,2) = d2f(2,1)
            
            dd = d2f(1,1)*d2f(2,2) - d2f(1,2)*d2f(1,2)
            
            fmin = huge(1.0)
            do iy = -1,1
                do ix = -1,1
                    if (f(ix,iy) < fmin) then
                        x = (/ ix,iy /)
                        fmin = f(ix,iy)
                    end if
                end do
            end do
            
            if (abs(dd)<1.0d-16) then
                fmin = huge(1.0)
                do iy = -1,1
                    do ix = -1,1
                        if (f(ix,iy) < fmin) then
                            x = (/ ix,iy /)
                            fmin = f(ix,iy)
                        end if
                    end do
                end do                
                return
            end if
            
            dd = 1/dd
            x(1) = ( -d2f(2,2)*df(1) + d2f(1,2)*df(2) )*dd
            x(2) = (  d2f(2,1)*df(1) - d2f(1,1)*df(2) )*dd
            
            fmin = f0 + df(1)*x(1) + df(2)*x(2) + ( d2f(1,1)*x(1)*x(1) + 2*d2f(1,2)*x(1)*x(2) + d2f(2,2)*x(2)*x(2) )/2
            do iy = -1,1
                do ix = -1,1
                    if (f(ix,iy) < fmin) then
                        x = (/ ix,iy /)
                        fmin = f(ix,iy)
                    end if
                end do
            end do                
            
            return
        end subroutine stencil               
                
    
        
        recursive subroutine marchingSquaresLine( f, n,u1,v1,u2,v2 )
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      given the four values on the corners of a square
    !*      find zero, one or two lines where the value is zero by bilinear interpolation
    !*      
            real(kind=real64),dimension(4),intent(in)       ::      f
            integer,intent(out)                             ::      n           !    = 0,1,2 = number of lines found
            real(kind=real64),dimension(2),intent(out)      ::      u1,v1,u2,v2     !   line from u to v
            
            integer             ::      code
            real(kind=real64)   ::      fc
            
            code = 0
            if (f(1)>=0) code = code + 1
            if (f(2)>=0) code = code + 2
            if (f(3)>=0) code = code + 4
            if (f(4)>=0) code = code + 8
    
           ! print *,"MSL ",code,f
            if (code>=8) then
!                 code = 15 - code
!                 fc = -fc
                call marchingSquaresLine( -f, n,u1,v1,u2,v2 )
                return
            end if
            
            fc = f(1)+f(2)+f(3)+f(4)    
            n = 1 ; u1 = -1 ; v1 = -1 ; u2 = -1 ; v2 = -1
            select case(code)
                case(0)
                    n = 0
                case(1)
                    u1(1:2) = (/ 0.0d0,linint(f(1),f(4)) /) 
                    v1(1:2) = (/ linint(f(1),f(2)),0.0d0 /) 
                case(2) 
                    u1(1:2) = (/ linint(f(1),f(2)),0.0d0 /) 
                    v1(1:2) = (/ 1.0d0,linint(f(2),f(3)) /) 
                case(3)                      
                    u1(1:2) = (/ 0.0d0,linint(f(1),f(4)) /) 
                    v1(1:2) = (/ 1.0d0,linint(f(2),f(3)) /) 
                case(4)                      
                    u1(1:2) = (/ linint(f(4),f(3)),1.0d0 /) 
                    v1(1:2) = (/ 1.0d0,linint(f(2),f(3)) /) 
                case(5)                      
                    n = 2                    
                    if (fc>=0) then
                        u1(1:2) = (/ 0.0d0,linint(f(1),f(4)) /) 
                        v1(1:2) = (/ linint(f(4),f(3)),1.0d0 /)                    
                        u2(1:2) = (/ linint(f(1),f(2)),0.0d0 /)  
                        v2(1:2) = (/ 1.0d0,linint(f(2),f(3)) /)                     
                    else
                        u1(1:2) = (/ 0.0d0,linint(f(1),f(4)) /) 
                        v1(1:2) = (/ linint(f(1),f(2)),0.0d0 /)                     
                        u2(1:2) = (/ linint(f(4),f(3)),1.0d0 /) 
                        v2(1:2) = (/ 1.0d0,linint(f(2),f(3)) /)                                         
                    end if
                case(6)
                    u1(1:2) = (/ linint(f(1),f(2)),0.0d0 /) 
                    v1(1:2) = (/ linint(f(4),f(3)),1.0d0 /) 
                case(7)
                    u1(1:2) = (/ 0.0d0,linint(f(1),f(4)) /) 
                    v1(1:2) = (/ linint(f(4),f(3)),1.0d0 /) 
                  
            end select
            
            return
            
        contains
    !---^^^^^^^^
    
            pure function linint( f1,f2 ) result(x)
        !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
        !*      form linear interpolation f = f1(1-x) + f2x
        !*      find f = 0 at x = f1/(f1-f2)
                real(kind=real64),intent(in)        ::      f1,f2
                real(kind=real64)                   ::      x
                x = 0.5d0
                if (abs(f1-f2)>1.0d-16) x = f1/(f1-f2)
                
                return
            end function linint
          
        end subroutine marchingSquaresLine
        
        
        subroutine closestPointLines( u,v, p1,p2 )
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      given the lines x_i = u_i + lambda_i v_i
    !*      find the point p which is closest in the least squares sense
    !*      ie minimise    S = ( x_i - p )^2
    !*      subject to  p in [0,1]
    !*      Minimum distance point to line is ( p - u - v [(p-u).v/v^2] )^2
    !*      which gives a simple linear least squares fit for closest point.
            real(kind=real64),dimension(:,:),intent(in)             ::      u,v
            real(kind=real64),intent(out)              ::      p1,p2
            
            integer             ::      nn
            integer             ::      ii
            real(kind=real64)   ::      uv,vv
            real(kind=real64)   ::      b1,b2
            real(kind=real64)   ::      a11,a12,a22
            
            p1 = -1.0d0
            p2 = -1.0d0
            nn = size(u,dim=2)
            if (nn<=1) return
                
            
            a11 = 0.0d0
            a12 = 0.0d0
            a22 = 0.0d0
            b1  = 0.0d0
            b2  = 0.0d0
            
            
          !  print *,"closestPointLines ",nn,"u: ",u,"v: ",v
            
            
            do ii = 1,nn
            
                uv = u(1,ii)*v(1,ii) + u(2,ii)*v(2,ii)  
                vv = v(1,ii)*v(1,ii) + v(2,ii)*v(2,ii)  
                
                a11 = a11 + vv - v(1,ii)*v(1,ii)  
                a12 = a12      - v(1,ii)*v(2,ii)  
                a22 = a22 + vv - v(2,ii)*v(2,ii)  
                b1  = b1 + u(1,ii)*vv - v(1,ii)*uv  
                b2  = b2 + u(2,ii)*vv - v(2,ii)*uv  
                 

            end do
            
        !---    find solution of 2x2 matrix equation
            vv = a11*a22 - a12*a12              !   determinant
            if (abs(vv)<1.0d-16)  return
          
            
            vv = 1/vv
            p1 = ( a22*b1 - a12*b2)*vv
            p2 = (-a12*b1 + a11*b2)*vv
            
           ! print *,"det ",1/vv," a ",a11,a12,a22," b ",b1,b2," p ",p1,p2
              
   !     !---    strictly speaking should do a constrained optimisation if out of bounds.
   !     !   but could just bring the point back to the boundary as a quick fix.
   !         p1 = max( 0.0d0, min( 1.0d0,p1 ) )
   !         p2 = max( 0.0d0, min( 1.0d0,p2 ) )
            
            
            return
        end subroutine closestPointLines
        
        
        
        
        
    end module SolveTripleJunction
    