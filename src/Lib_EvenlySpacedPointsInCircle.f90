
    module Lib_EvenlySpacedPointsInCircle
!---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
!*      Place n points in a disk radius 1, so that they are maximally spaced.
!*      If the number of points is large, use a modification of the sunflower seed placement
!*      But what if the number of points is small?
!*      Then I have found some probably suboptimal selections here
!*  
!*      Daniel Mason
!*      (c) UKAEA August 2023
!*

        use iso_fortran_env
        implicit none
        private
         
        include "evenlySpaced.h"
        
        real(kind=real64),parameter     ::      PI = 3.141592654d0
        
        public          ::      evenlySpacedPointsInCircle
        public          ::      findEvenlySpacedPointsInCircle

    contains
!---^^^^^^^^

    

        subroutine evenlySpacedPointsInCircle(x,file,force)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
            real(kind=real64),dimension(:,:),intent(out)        ::      x
            character(len=*),intent(in),optional                ::      file
            logical,intent(in),optional                         ::      force
            
            integer                 ::      nn,nKnown,ii,mm 
            logical                 ::      ok
            character(len=4096)     ::      dummy
            real(kind=real64)       ::      e0,ee
            
            
            
            
            nn = size(x,dim=2)
            x = 0
            
            if (present(file)) then
                open(unit=500,file=trim(file),action="read")
                    read(500,fmt=*) nKnown
                    ok = .false.                     !   assume we haven't found solution
                    x(1:2,1:nn) = sunflower(nn, alpha=0.0d0)      !   default to sunflower solution
                    e0 = energy( x ) 
                    do 
                        read(500,fmt='(i6,a)',iostat=ii) mm,dummy
                        if (ii /= 0) exit
                        if (mm==nn) then
                            read(dummy,fmt=*) e0,x(1:2,1:nn)
                            print *,"Lib_EvenlySpacedPointsInCircle::findEvenlySpacedPointsInCircle info - solution ",e0
                            ok = .true.            !   have got a solution
                        end if
                    end do
                close(unit=500)  
                if (ok) then
                    !   we have a solution
                    if (present(force)) then
                        if (force) then
                            !   try to find a better one anyway.
                            call findEvenlySpacedPointsInCircle(x,ee)
                            if (ee < e0) then   
                                !   have improved the solution.
                                print *,"Lib_EvenlySpacedPointsInCircle::findEvenlySpacedPointsInCircle info - improved ",ee    
                                open(unit=600,file=trim(file),action="write",position="append")
                                   write(600,fmt='(i6,f16.8,10000f8.4)') nn,ee,x(1:2,1:nn)
                                close(unit=600)  
                            end if
                        end if
                    end if
                else
                    !   we don't have a solution. Stick with sunflower?
                    if (nn<=MAXKNOWN) then
                        call findEvenlySpacedPointsInCircle(x,ee)
                        open(unit=600,file=trim(file),action="write",position="append")
                           write(600,fmt='(i6,f16.8,10000f8.4)') nn,ee,x(1:2,1:nn)
                        close(unit=600)  
                    end if
                end if                
            else
                if (nn<=MAXKNOWN) then  
                    x(1:2,1:nn) = xKnown(1:2,1:nn,nn) 
                else
                    x(1:2,1:nn) = sunflower(nn, alpha=0.0d0)      !   default to sunflower solution
                end if
            end if
            
            
            
            
            
            return
        end subroutine evenlySpacedPointsInCircle
            
    

        subroutine findEvenlySpacedPointsInCircle(x,e)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
            real(kind=real64),dimension(:,:),intent(inout)      ::      x
            real(kind=real64),intent(out)                       ::      e
            integer                 ::      nn
            
            
            real(kind=real64)       ::      iT,dAngle
            real(kind=real64),dimension(2,size(x,dim=2))    ::      yy,xbest
            integer                 ::      ii,loop,trial,accept
            real(kind=real64)       ::      ebest
            logical                 ::      ok
            real(kind=real64)       ::      ff,gg
            
            
            nn = size(x,dim=2)
            
           ! x = sunflower(nn, alpha=1.0d0) 
            xbest = x
            ff = 10000.0**(0.1)
            gg = 0.025**(0.1)
            e = areal_energy( x,NGRID=50 )  
            ebest = e
            print *,"Lib_EvenlySpacedPointsInCircle::findEvenlySpacedPointsInCircle info - initial ",e
            
            
            
            do trial = 1,10
            
        !---    convert from real space to (/phi,theta/) so that radius stays in range 0:1          
                x = xbest
                do ii = 1,nn    
                    call backChangeOfCoords( x(:,ii),yy(:,ii) )
                end do
                
                iT = 4/ebest
                e = areal_energy( x,NGRID=50 )  
 
                dAngle = 0.4d0
                
                do loop = 1,10
                    
                    accept = 0
                    do ii = 1,100*nn*nn
                        call MetropolisMove( x,yy,e,dAngle,iT,ok )
                        if (ok) accept = accept + 1
                        !print *,ii,ee,iT
                    end do
                    
                    iT = iT*ff
                    dAngle = dAngle*gg 
                    
                    !print *,"loop ",e,iT,dAngle,real(accept)/(100*nn*nn)
                    
                end do            
                print *,trial,e 
                
                if (e < ebest) then
                    xbest = x
                    ebest = e                    
                end if
                 
                
                
            end do
            
            
            
            do trial = 1,3
            
        !---    convert from real space to (/phi,theta/) so that radius stays in range 0:1          
                x = xbest
                do ii = 1,nn    
                    call backChangeOfCoords( x(:,ii),yy(:,ii) )
                end do
                
                iT = 20/ebest
                e = areal_energy( x,NGRID=50 )  
 
                dAngle = 0.1d0
                
                do loop = 1,10
                    
                    accept = 0
                    do ii = 1,100*nn*nn
                        call MetropolisMove( x,yy,e,dAngle,iT,ok )
                        if (ok) accept = accept + 1
                        !print *,ii,ee,iT
                    end do
                    
                    iT = iT*ff
                    dAngle = dAngle*gg 
                    
                    !print *,"loop ",e,iT,dAngle,real(accept)/(100*nn*nn)
                    
                end do            
                print *,trial,e 
                
                if (e < ebest) then
                    xbest = x
                    ebest = e                    
                end if
                 
                
                
            end do
            
            
            x = xbest
            e = ebest
            print *,"Lib_EvenlySpacedPointsInCircle::findEvenlySpacedPointsInCircle info - best ",ebest    
            return
        end subroutine findEvenlySpacedPointsInCircle
            
        
        subroutine MetropolisMove( x,y,e,dAngle,iT,ok )
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ 
    !*      make a random angle selection +/- dAngle
    !*      and accept with probability Exp[ -deltaE/T ]
            real(kind=real64),dimension(:,:),intent(inout)      ::      x,y
            real(kind=real64),intent(inout)                     ::      e
            real(kind=real64),intent(in)                        ::      dAngle,iT
            logical,intent(out)                                 ::      ok
            
            real(kind=real64)               ::      trialEnergy,deltaEonT,pp
            
            integer                         ::      ii,nn
            real(kind=real64),dimension(2)  ::      xi0     !,yi0
            real(kind=real64),dimension(4)  ::      zeta
            
            nn = size(x,dim=2)
            
            call random_number(zeta)
            zeta(1:2) = dAngle*( 2*zeta(1:2) - 1.0d0 )          !   range +/- dangle
            ii = floor( zeta(3)*nn )+ 1                         !   range 1:nn
            
            xi0 = x(:,ii)                                       !   store position in case of rejection
            call changeOfCoords( y(1:2,ii)+zeta(1:2) , x(1:2,ii) )       
            trialEnergy = areal_energy( x )     
            deltaEonT = ( trialEnergy - e )*iT
            !print *,"mm ",ii,zeta(1:2),e,trialEnergy,deltaEonT,exp( - deltaEonT ),zeta(4)
            if (deltaEonT>0) then
                pp = exp( - deltaEonT )
                if (pp < zeta(4)) then
                    x(1:2,ii) = xi0
                    ok = .false.
                    return
                end if
            end if
            e = trialEnergy
            y(1:2,ii) = y(1:2,ii)+zeta(1:2)
            ok = .true.
            return
        end subroutine MetropolisMove
            
        
         
        
        
            
        pure real(kind=real64) function energy( x ) 
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      compute energy of set of points x
    !*          E = 1/2 sum_ij  |r_ij|^{-4}
            real(kind=real64),dimension(:,:),intent(in)     ::  x
            integer             ::      nn
            integer             ::      ii,jj
            real(kind=real64)   ::      dx,dy,dd
            
            nn = size(x,dim=2)
            energy = 0.0d0
            do jj = 1,nn-1
                do ii = jj+1,nn
                    dx = x(1,ii) - x(1,jj)
                    dy = x(2,ii) - x(2,jj)
                    dd = dx*dx + dy*dy
                    dd = 1/max(1.0d-8,dd)
                    
                    energy = energy + dd*dd
                end do
            end do
            return
        end function energy
        
        
        
            
        real(kind=real64) function areal_energy( x,NGRID ) 
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      compute energy of set of points x 
    !*          E = <d^2> = mean square distance of point x_i to point in circle
    !*      
            real(kind=real64),dimension(:,:),intent(in)             ::      x            
            integer,intent(in),optional                             ::      NGRID
            integer,save                                            ::      nArealPoints = 0
            real(kind=real64),dimension(:,:),allocatable,save       ::      arealPoint
            
            
            
            integer             ::      nn,ng
            integer             ::      ii,jj,j_best
            real(kind=real64)   ::      dx,dy,dd,d_best
            real(kind=real64),dimension(size(x,dim=2))      ::      d2sum
            integer,dimension(size(x,dim=2))                ::      nnear

            ng = 20 
            if (present(NGRID)) then
                ng = NGRID
                if (nArealPoints/=0) then
                    deallocate(arealPoint)
                    nArealPoints = 0
                end if
            end if
            if (nArealPoints==0) then
            !---    save the area points
                do jj = -ng,ng
                    dy = real(jj)/ng
                    do ii = -ng,ng
                        dx = real(ii)/ng
                        dd = dx*dx + dy*dy
                        if (dd<1.0d0) nArealPoints = nArealPoints + 1
                    end do
                end do
                allocate(arealPoint(2,nArealPoints))
                nArealPoints=0
                do jj = -ng,ng
                    dy = real(jj)/ng
                    do ii = -ng,ng
                        dx = real(ii)/ng
                        dd = dx*dx + dy*dy
                        if (dd<1.0d0) then
                            nArealPoints = nArealPoints + 1
                            arealPoint(:,nArealPoints) = (/ dx,dy /)
                        end if                            
                    end do
                end do
            end if                
            
            
            
        !---    For each point on the grid, find the nearest x point
        !       and record its distance squared   
            nn = size(x,dim=2)
            d2sum = 0.0d0
            nnear = 0
            do ii = 1,nArealPoints
                d_best = huge(1.0)
                do jj = 1,nn
                    dx = arealPoint(1,ii) - x(1,jj)
                    dy = arealPoint(2,ii) - x(2,jj)
                    dd = dx*dx + dy*dy
                    if (dd<d_best) then
                        j_best = jj
                        d_best = dd
                    end if
                end do
                
                d2sum(j_best) = d2sum(j_best) + d_best
                nnear(j_best) = nnear(j_best) + 1
                
            end do
            
            areal_energy = 0.0d0
            do jj = 1,nn
                dd = d2sum(jj) / max( nnear(jj),1 )
                areal_energy = areal_energy + dd
            end do
            
            return
        end function areal_energy
        
        
        
        
        
        
        
        pure subroutine changeOfCoords( y, x,dx )
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      given coordinate system y = (/phi,theta/)
    !*      return x = cos^2 phi (/ cos(theta),sin(theta) /)
    !*      and dx(:,1) = (/ d x1/dphi,d x1/dtheta /)
    !*          dx(:,2) = (/  d x2/dphi, d x2/dtheta /)
            real(kind=real64),dimension(2),intent(in)       ::  y
            real(kind=real64),dimension(2),intent(out)      ::  x
            real(kind=real64),dimension(2,2),intent(out),optional    ::  dx
            
            real(kind=real64)       ::      cosp,sinp,cost,sint
            real(kind=real64)       ::      phi,theta
            
            phi = y(1)
            theta = y(2)
            cosp = cos(phi)
            sinp = sin(phi)
            cost = cos(theta)
            sint = sin(theta)
            
            x = (/ cosp*cosp*cost,cosp*cosp*sint /)
            if (present(dx)) then
                dx(:,1) = (/ 2*cosp*sinp*cost , -cosp*cosp*sint /)
                dx(:,2) = (/ 2*cosp*sinp*sint , cosp*cosp*cost /)
            end if            
            
            return
        end subroutine changeOfCoords
        
        pure subroutine backChangeOfCoords( x,y )
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      given x = cos^2 phi (/ cos(theta),sin(theta) /)
    !*      return y = (/phi,theta)
            real(kind=real64),dimension(2),intent(in)       ::  x
            real(kind=real64),dimension(2),intent(out)      ::  y
            
            real(kind=real64)       ::      phi4,phi2
            
            
            phi4 = x(1)*x(1) + x(2)*x(2)
            
            if (phi4>=1) then
                y(1) = 0.0d0
                y(2) = acos( x(1) / sqrt(phi4) )
            else if (phi4>1.0d-16) then
                phi2 = sqrt(phi4)
                y(1) = acos( sqrt(phi2) )
                y(2) = acos( x(1)/phi2 )
            else
                y(1) = 0.0d0
                y(2) = 0.0d0
            end if
            
            return
        end subroutine backChangeOfCoords
        
                
        
        pure function sunflower(n, alpha) result(x) 
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      modified sunflower algorithm
    !*      Peter Mortensen
    !*      https://stackoverflow.com/questions/28567166/uniformly-distribute-x-points-inside-a-circle
            integer,intent(in)                  ::      n
            real(kind=real64),intent(in)        ::      alpha
            real(kind=real64),dimension(2,n)    ::      x
            
            integer                         ::      bb                              !   number of boundary points
            real(kind=real64),parameter     ::      PHI = (SQRT(5.0d0)+1.0d0)/2     !   Golden Ratio
            
            
            integer             ::      kk
            real(kind=real64)   ::      rr,theta
            
             bb = nint( sqrt( alpha*alpha*n ) )         !   number of boundary points
            
            do kk = 1,n
                rr = radius(kk,n,bb)
                theta = 2*PI*kk/(PHI*PHI)
                x(1:2,kk) = rr*(/ cos(theta),sin(theta) /)
            end do
            
            return
        end function sunflower
        
        pure real(kind=real64) function radius( k,n,b )
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
            integer,intent(in)          ::      k       !   point number
            integer,intent(in)          ::      n       !   max point number
            integer,intent(in)          ::      b       !   number of points on boundary
              
            if (k > n - b) then 
                radius = 1                       !   boudary point
            else
                radius = real(2*k - 1,kind=real64)/(2*n-b-1)
                radius = sqrt(radius)
            end if
            return
        end function radius
                
 



    end module Lib_EvenlySpacedPointsInCircle
    

    

        