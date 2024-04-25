
    module Lib_DistancePointToTriangle
!---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
!*      compute the signed distance from a point in 3d to a triangle in 3d.
        use iso_fortran_env
        implicit none
        private 
        
        
        public          ::      distancePointToTriangle
        
        interface   distancePointToTriangle
            module procedure    distancePointToTriangle0
            module procedure    distancePointToTriangle1
        end interface
            
        
        
    contains
!---^^^^^^^^

 
        pure function distancePointToTriangle0( p0,p1,p2,p3 ) result( d )
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*               /                        .
    !*       6     /                          By placing z' along p13 
    !*           /                            and p12 in the y-z plane
    !*     -----*   1                         we have a set of regions where p10 might be 
    !*          |\      /
    !*          | \   /                       We need to test where we are (0-6) 
    !*      5   |0 \/   2                     and then find the distance accordingly.
    !*          |  /\                         .
    !*          | /   \                       The starred points are the triangle vertices, which in the x',y',z' frame are
    !*          |/      \                         ( 0, 0, 0 )
    !*     -----*   3                             ( 0, |p13|² |p12|² - (p12.p13)²  , (p12.p13) )
    !*           \                                ( 0, 0, |p13|² )
    !*       4     \                           .
    !*               \                         .
                
            real(kind=real64),dimension(3),intent(in)       ::      p0,p1,p2,p3
            real(kind=real64)                               ::      d          
    
            real(kind=real64),dimension(3)      ::      p12,p13,p10 
            real(kind=real64),dimension(3)      ::      xp,yp
            
            real(kind=real64)                   ::      qq,q2,ypp2,q3,xpp1,xpp2,xpp3
            real(kind=real64)                   ::      p12dotp12,p12dotp13,p13dotp13
            
            
        !---    subtract p1 offset so that p1 is at the origin.
            p10(1:3) = p0(1:3) - p1(1:3)
            p12(1:3) = p2(1:3) - p1(1:3)
            p13(1:3) = p3(1:3) - p1(1:3)
            
        !---    find new x- axis direction 
            xp(1) = p12(2)*p13(3) - p12(3)*p13(2)
            xp(2) = p12(3)*p13(1) - p12(1)*p13(3)
            xp(3) = p12(1)*p13(2) - p12(2)*p13(1)
        
           ! print *,"p10 ",p10
           ! print *,"p12 ",p12
           ! print *,"p13 ",p13
           ! print *,"x'  ",xp
            
            
            
        !---    check - does this triangle in fact have zero area? If so its a problem...
            qq = xp(1)*xp(1) + xp(2)*xp(2) + xp(3)*xp(3)        ! = 4*area^2           
            if (qq < 1.0d-16) then
                !   yes. Zero area - big problem. 
                !   Can't return a signed area properly. Try our best anyway...
                !   is the "triangle" a line or a point?
                qq = p12(1)*p12(1) + p12(2)*p12(2) + p12(3)*p12(3)      !   distance squared p1 to p2.
                if (qq < 1.0d-16) then
                    !   p1 = p2
                    q2 = p13(1)*p13(1) + p13(2)*p13(2) + p13(3)*p13(3)      !   distance squared p1 to p3.
                    if (q2 < 1.0d-16) then
                        !   p1 = p2 = p3. Point
                        d = norm2( p10 ) ; return
                    else
                        !   p1 = p2 /= p3. Line.
                        !   line is x = lambda p13. Point is p10
                        !   min distance when lambda = p10.p13/p13.p13
                        d = sqrt( distanceSquaredPointToLine( p13,p10 ) )
                    end if
                else
                    !   p1 /= p2. Is p1=p3 or p2=p3?
                    q2 = p13(1)*p13(1) + p13(2)*p13(2) + p13(3)*p13(3)      !   distance squared p1 to p3.   
                    if (q2 < 1.0d-16) then
                        !   p1 = p3 /= p2. Line.
                        !   line is x = lambda p12. Point is p10
                        !   min distance when lambda = p10.p12/p12.p12
                        d = sqrt( distanceSquaredPointToLine( p12,p10 ) )
                    else
                        !   p1 /= p2 = p3. Line.
                        !   line is x = lambda p13. Point is p10
                        !   min distance when lambda = p10.p13/p13.p13
                        d = sqrt( distanceSquaredPointToLine( p13,p10 ) )
                    end if    
                end if  !   p1=p2?
                
                
            end if      !   zero area triangle?
        !---    end of pathological point/line question

    
    
    !---    so now I know p1,p2,p3 form a triangle.
    !       which sector am I in?
            p12dotp12 = p12(1)*p12(1) + p12(2)*p12(2) + p12(3)*p12(3) 
            p12dotp13 = p12(1)*p13(1) + p12(2)*p13(2) + p12(3)*p13(3) 
            p13dotp13 = p13(1)*p13(1) + p13(2)*p13(2) + p13(3)*p13(3)  
            
             
            yp(1:3) = p12(1:3)*p13dotp13 - p12dotp13*p13(1:3)         
                            
            xpp1 = xp(1)*p10(1) + xp(2)*p10(2) + xp(3)*p10(3) 
            xpp2 = yp(1)*p10(1) + yp(2)*p10(2) + yp(3)*p10(3)          
            xpp3 = p13(1)*p10(1) + p13(2)*p10(2) + p13(3)*p10(3)
            
               
            ypp2 = p12dotp12*p13dotp13  - p12dotp13*p12dotp13        
             
            if (xpp2 >= 0) then    !   right of z-axis                 
                if ( xpp2*p12dotp13 - xpp3*ypp2 <= 0 ) then
                    !   left of p12                     
                    if (xpp2*(p12dotp13-p13dotp13) - (xpp3-p13dotp13)*ypp2 > 0) then
                        d = xpp1/norm2(xp)
                        return
                    end if
                end if
            end if
                     
            qq = distanceSquaredPointToLine( p2-p3,p0-p3 )            
            q2 = distanceSquaredPointToLine( p12,p10 )
            q3 = distanceSquaredPointToLine( p13,p10 )
           
            
            qq = min( qq,min(q2,q3) )
            d  = sign( sqrt(qq),xpp1 )
            
            
            
            
            return
           
             
        end function distancePointToTriangle0
    
        pure function distancePointToTriangle1( p0,p1,p2,p3 ) result( d )
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*               /                       .
    !*       6     /                         By placing z' along p13 
    !*           /                           and p12 in the y-z plane
    !*     -----*   1                        we have a set of regions where p10 might be 
    !*          |\      /                    .
    !*          | \   /                      We need to test where we are (0-6) 
    !*      5   |0 \/   2                    and then find the distance accordingly.
    !*          |  /\                        .
    !*          | /   \                      The starred points are the triangle vertices, which in the x',y',z' frame are
    !*          |/      \                         ( 0, 0, 0 )
    !*     -----*   3                             ( 0, |p13|² |p12|² - (p12.p13)²  , (p12.p13) )
    !*           \                                ( 0, 0, |p13|² )
    !*       4     \                         .
    !*               \                       .
                
            real(kind=real64),dimension(:,:),intent(in)       ::      p0
            real(kind=real64),dimension(3),intent(in)         ::      p1,p2,p3
            real(kind=real64),dimension(size(p0,dim=2))       ::      d          
    
            real(kind=real64),dimension(3)      ::      p12,p13,p10 
            real(kind=real64),dimension(3)      ::      xp,yp
            
            real(kind=real64)                   ::      qq,q2,ypp2,q3,xpp1,xpp2,xpp3
            real(kind=real64)                   ::      p12dotp12,p12dotp13,p13dotp13 
            
            integer         ::      ii,nn
            
            nn = size(p0,dim=2)
            
            
        !---    subtract p1 offset so that p1 is at the origin.
           ! p10(1:3) = p0(1:3) - p1(1:3)
            p12(1:3) = p2(1:3) - p1(1:3)
            p13(1:3) = p3(1:3) - p1(1:3)
            
        !---    find new x- axis direction 
            xp(1) = p12(2)*p13(3) - p12(3)*p13(2)
            xp(2) = p12(3)*p13(1) - p12(1)*p13(3)
            xp(3) = p12(1)*p13(2) - p12(2)*p13(1)
        
            
        !---    check - does this triangle in fact have zero area? If so its a problem...
            qq = xp(1)*xp(1) + xp(2)*xp(2) + xp(3)*xp(3)        ! = 4*area^2           
            if (qq < 1.0d-16) then
                !   yes. Zero area - big problem. 
                !   Can't return a signed area properly. Try our best anyway...
                !   is the "triangle" a line or a point?
                qq = p12(1)*p12(1) + p12(2)*p12(2) + p12(3)*p12(3)      !   distance squared p1 to p2.
                if (qq < 1.0d-16) then
                    !   p1 = p2
                    q2 = p13(1)*p13(1) + p13(2)*p13(2) + p13(3)*p13(3)      !   distance squared p1 to p3.
                    if (q2 < 1.0d-16) then
                        !   p1 = p2 = p3. Point
                        do ii = 1,nn
                            p10(1:3) = p0(1:3,ii) - p1(1:3)
                            d(ii) = norm2( p10 )
                        end do    
                        
                    else
                        !   p1 = p2 /= p3. Line.
                        !   line is x = lambda p13. Point is p10
                        !   min distance when lambda = p10.p13/p13.p13
                        do ii = 1,nn
                            p10(1:3) = p0(1:3,ii) - p1(1:3)
                            d(ii) = sqrt( distanceSquaredPointToLine( p13,p10 ) )
                        end do    
                        
                    end if
                else
                    !   p1 /= p2. Is p1=p3 or p2=p3?
                    q2 = p13(1)*p13(1) + p13(2)*p13(2) + p13(3)*p13(3)      !   distance squared p1 to p3.   
                    if (q2 < 1.0d-16) then
                        !   p1 = p3 /= p2. Line.
                        !   line is x = lambda p12. Point is p10
                        !   min distance when lambda = p10.p12/p12.p12
                        do ii = 1,nn
                            p10(1:3) = p0(1:3,ii) - p1(1:3)
                            d(ii) = sqrt( distanceSquaredPointToLine( p12,p10 ) )
                        end do    
                        
                    else
                        !   p1 /= p2 = p3. Line.
                        !   line is x = lambda p13. Point is p10
                        !   min distance when lambda = p10.p13/p13.p13
                        do ii = 1,nn
                            p10(1:3) = p0(1:3,ii) - p1(1:3)
                            d(ii) = sqrt( distanceSquaredPointToLine( p13,p10 ) )
                        end do    
                        
                    end if    
                end if  !   p1=p2?
                
                return
            end if      !   zero area triangle?
        !---    end of pathological point/line question

    
    
    !---    so now I know p1,p2,p3 form a triangle.
    !       which sector am I in?
            p12dotp12 = p12(1)*p12(1) + p12(2)*p12(2) + p12(3)*p12(3) 
            p12dotp13 = p12(1)*p13(1) + p12(2)*p13(2) + p12(3)*p13(3) 
            p13dotp13 = p13(1)*p13(1) + p13(2)*p13(2) + p13(3)*p13(3)  
            
             
            yp(1:3) = p12(1:3)*p13dotp13 - p12dotp13*p13(1:3)         
                            
            ypp2 = p12dotp12*p13dotp13  - p12dotp13*p12dotp13        
            
            do ii = 1,nn
                p10(1:3) = p0(1:3,ii) - p1(1:3)
                xpp1 = xp(1)*p10(1) + xp(2)*p10(2) + xp(3)*p10(3) 
                xpp2 = yp(1)*p10(1) + yp(2)*p10(2) + yp(3)*p10(3)          
                xpp3 = p13(1)*p10(1) + p13(2)*p10(2) + p13(3)*p10(3)
                 
                 
                if (xpp2 >= 0) then    !   right of z-axis                 
                    if ( xpp2*p12dotp13 - xpp3*ypp2 <= 0 ) then
                        !   left of p12                     
                        if (xpp2*(p12dotp13-p13dotp13) - (xpp3-p13dotp13)*ypp2 > 0) then
                            d(ii) = xpp1/norm2(xp)
                            cycle                            
                        end if
                    end if
                end if
                         
                qq = distanceSquaredPointToLine( p2-p3,p0(:,ii)-p3 )            
                q2 = distanceSquaredPointToLine( p12,p10 )
                q3 = distanceSquaredPointToLine( p13,p10 )
               
                
                qq    = min( qq,min(q2,q3) )
                d(ii) = sign( sqrt(qq),xpp1 )
            end do            
                
            
            return
           
       
        end function distancePointToTriangle1
 
        pure function distanceSquaredPointToLine( v,x ) result( d2 )
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      find the closest distance from a point x to the line 0 + lambda v.
    !*      if lambda < 0 or lambda > 1 then return distance to vertex.
            real(kind=real64),dimension(3),intent(in)       ::      v,x
            real(kind=real64)                               ::      d2
            
            real(kind=real64)           ::      xdotx,vdotx,vdotv
            xdotx =   x(1)*x(1) + x(2)*x(2) + x(3)*x(3)     
            vdotx =   x(1)*v(1) + x(2)*v(2) + x(3)*v(3)     
            
            d2 = xdotx
            vdotv = v(1)*v(1) + v(2)*v(2) + v(3)*v(3) 
            if (vdotx>=vdotv) then 
                d2 = d2 - 2*vdotx + vdotv
            else if (vdotx>0) then
                d2 = d2 - vdotx*vdotx/vdotv
            end if 
            
            
            d2 = max(0.0d0,d2)      !   to avoid -0 rounding errors later on.
            return
        end function distanceSquaredPointToLine
             
        
        

    end module Lib_DistancePointToTriangle
    
!    