
    module Lib_MarchingSquares
!---^^^^^^^^^^^^^^^^^^^^^^^^^^
        use Lib_LinkCell3D      !   to find contiguous lines
        use iso_fortran_env
        implicit none
        private
        
    !---
    
        public      ::      marchingSquares  
        public      ::      marchingSquaresArea
        public      ::      nLinesNeeded
        public      ::      contiguousLines
        public      ::      areaAndLengthAndCurvContiguousLine
        
        
    !---
    
        logical,public          ::      MarchingSquare_dbg = .false.
                    
        
    contains
!---^^^^^^^^ 

        pure function getCode( f1,f2,f3,f4 ) result(c)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !
    !       f1---f2    1----2
    !       |    |     |    |
    !       |    |     |    |
    !       f3---f4    4----8
    !
    !       note:   o----2  o----2               1----o  1----o
    !               | / /|  |   \|               |\ \ |  |/   |
    !               |/ / |  |\   |               | \ \|  |   /|
    !               3----o  3----o               o----4  o----4
    !                 6       16                    9      17   
    !
    
            real(kind=real64),intent(in)        ::      f1,f2,f3,f4
            integer                             ::      c
            c = 0
            if (f1 >= 0) c = 1
            if (f2 >= 0) c = c + 2
            if (f3 >= 0) c = c + 4
            if (f4 >= 0) c = c + 8
            
            if (c == 6) then
                if (f1+f2+f3+f4<0) c = 16
            else if (c == 9) then
                if (f1+f2+f3+f4<0) c = 17
            end if
            
            return
        end function getCode
             

        function nLinesNeeded( img,isolevel,exact ) result(nLines)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      (optionally estimate) the number of lines needed at this isolevel
            real(kind=real64),dimension(0:,0:),intent(in)       ::      img
            real(kind=real64),intent(in)                        ::      isolevel
            logical,intent(in)                                  ::      exact
            integer             ::      nLines
            integer             ::      nx,ny,ix,iy,jx,jy , code
            real(kind=real64)   ::      f1,f2,f3,f4
            !integer             ::      cc
            
            integer,dimension(0:17),parameter       ::      NL = (/ 0,1,1,1, 1,1,2,1, 1,2,1,1, 1,1,1,0 , 2,2 /)
            
            nx = size(img,dim=1)
            ny = size(img,dim=2)
            nLines = 0
            do iy = 0,ny-1  
                jy = iy + 1 ; if (iy==ny-1) jy = 0
                do ix = 0,nx-1
                    jx = ix + 1 ; if (ix==nx-1) jx = 0

                    f1 = img(ix,iy) - isolevel ; f2 = img(jx,iy) - isolevel 
                    f3 = img(ix,jy) - isolevel ; f4 = img(jx,jy) - isolevel       
                    
                    code = getCode( f1,f2,f3,f4 )
                    
                    nLines = nLines + NL(code)
                 
                    
                end do
            end do
            
            
            
            
            
!            
!            
!            if (.not. exact) then
!                nLines = 4*count( img>isolevel )        !   worst case scenario is all points separated surrounded by quadrilaterals
!            else
!                nx = size(img,dim=1)
!                ny = size(img,dim=2)
!                nLines = 0
!                do iy = 0,ny-2
!                    f1 = img(nx-1,iy) - isolevel ; f3 = img(nx-1,iy+1) - isolevel 
!                    f2 = img(0,iy) - isolevel ; f4 = img(0,iy+1) - isolevel 
!                    cc = getCode( f1,f2,f3,f4 )
!                    nLines = nLines + NL(cc)
!                    do ix = 0,nx-2
!                        f1 = f2 ; f3 = f4
!                        f2 = img(ix+1,iy) - isolevel ; f4 = img(ix+1,iy+1) - isolevel
!                        cc = getCode( f1,f2,f3,f4 )
!                        nLines = nLines + NL(cc)
!                    end do
!                end do
!                f1 = img(nx-1,ny-1) - isolevel ; f3 = img(nx-1,0) - isolevel 
!                f2 = img(0,ny-1) - isolevel ; f4 = img(0,0) - isolevel 
!                cc = getCode( f1,f2,f3,f4 )
!                nLines = nLines + NL(cc)
!                do ix = 0,nx-2
!                    f1 = f2 ; f3 = f4
!                    f2 = img(ix+1,ny-1) - isolevel ; f4 = img(ix+1,0) - isolevel
!                    cc = getCode( f1,f2,f3,f4 )
!                    nLines = nLines + NL(cc)
!                end do              
!            end if      
            return
        end function nLinesNeeded   

        subroutine marchingSquares( img,isolevel,nLines,line )
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      Note: assumes PBC
    !*      If PBC are not needed, then make sure img(0,:) = img(:,0) < isolevel, then will have correctly closed contours
    !*
    
            real(kind=real64),dimension(0:,0:),intent(in)       ::      img
            real(kind=real64),intent(in)                        ::      isolevel
            integer,intent(out)                                 ::      nLines
            real(kind=real64),dimension(:,:),intent(out)        ::      line            !   (4,nLines) (fromx,fromy,tox,toy)
            integer             ::      nx,ny,ix,iy , jx,jy , code
            real(kind=real64)   ::      f1,f2,f3,f4 
            
            nx = size(img,dim=1)
            ny = size(img,dim=2)
            nLines = 0
            line = 0.
            !print *,"Lib_MarchingSquares::marchingSquares info - nx,ny = ",nx,ny
            
            
!            do iy = 0,ny-2
!                f1 = img(nx-1,iy) - isolevel ; f3 = img(nx-1,iy+1) - isolevel 
!                f2 = img(0,iy) - isolevel ; f4 = img(0,iy+1) - isolevel 
!                call addLines( getCode( f1,f2,f3,f4 ),real(nx-1,kind=real64),real(iy,kind=real64) )             
!                do ix = 0,nx-2
!                    f1 = f2 ; f3 = f4
!                    f2 = img(ix+1,iy) - isolevel ; f4 = img(ix+1,iy+1) - isolevel                   
!                    call addLines( getCode( f1,f2,f3,f4 ),real(ix,kind=real64),real(iy,kind=real64) )
!                end do
!            end do
!            f1 = img(nx-1,ny-1) - isolevel ; f3 = img(nx-1,0) - isolevel 
!            f2 = img(0,ny-1) - isolevel ; f4 = img(0,0) - isolevel 
!            call addLines( getCode( f1,f2,f3,f4 ),real(nx-1,kind=real64),real(ny-1,kind=real64) )
!            do ix = 0,nx-2
!                f1 = f2 ; f3 = f4
!                f2 = img(ix+1,ny-1) - isolevel ; f4 = img(ix+1,0) - isolevel
!                call addLines( getCode( f1,f2,f3,f4 ),real(ix,kind=real64),real(ny-1,kind=real64) )
!            end do

            do iy = 0,ny-1  
                jy = iy + 1 ; if (iy==ny-1) jy = 0
                do ix = 0,nx-1
                    jx = ix + 1 ; if (ix==nx-1) jx = 0

                    f1 = img(ix,iy) - isolevel ; f2 = img(jx,iy) - isolevel 
                    f3 = img(ix,jy) - isolevel ; f4 = img(jx,jy) - isolevel       
                    
                    code = getCode( f1,f2,f3,f4 )
                    
                    call addLines( code,real(ix,kind=real64),real(iy,kind=real64) )
                    
                 !  if (nLines <= 10) then
                 !      print *,"marching squares ",nLines,ix,iy,jx,jy,f1,f2,f3,f4,code ,line(1:4,nLines)
                 !  end if
                    
                    
                end do
            end do
            
            return

        contains
    !---^^^^^^^^
    
            subroutine addLines(c,x,y)
        !---^^^^^^^^^^^^^^^^^^^^^^^^^^
        !       codes:
        !
        !       f1---f2    1----2
        !       |    |     |    |
        !       |    |     |    |
        !       f3---f4    4----8
        !
        !       note:   o----2  o----2               1----o  1----o
        !               | / /|  |   \|               |\ \ |  |/   |
        !               |/ / |  |\   |               | \ \|  |   /|
        !               3----o  3----o               o----4  o----4
        !                 6       16                    9      17   
        !
        !       attempt to close curves _anti clockwise_
        
                integer,intent(in)                  ::      c
                real(kind=real64),intent(in)        ::      x,y
                select case(c)
                    case(1)
                        nLines = nLines + 1
                        line(1:4,nLines) = (/ x                      ,y+f1/(f1-f3),x+f1/(f1-f2)         ,y            /)                                    
                    case(2)                                                                                                                       
                        nLines = nLines + 1                                                                                                       
                        line(1:4,nLines) = (/ x+f1/(f1-f2),y                      ,x+1.0d0               ,y+f2/(f2-f4) /)                       
                    case(3)
                        nLines = nLines + 1
                        line(1:4,nLines) = (/ x                      ,y+f1/(f1-f3),x+1.0d0               ,y+f2/(f2-f4) /)                       
                    case(4)                                                                    
                        nLines = nLines + 1                                                    
                        line(1:4,nLines) = (/ x+f3/(f3-f4),y+1.0d0                ,x                     ,y+f1/(f1-f3) /)
                    case(5)
                        nLines = nLines + 1
                        line(1:4,nLines) = (/ x+f3/(f3-f4),y+1.0d0                ,x+f1/(f1-f2)         ,y      /)
                    case(6)
                        nLines = nLines + 1
                        line(1:4,nLines) = (/ x+f3/(f3-f4),y+1.0d0                ,x+1.0d0               ,y+f2/(f2-f4) /)
                        nLines = nLines + 1
                        line(1:4,nLines) = (/ x+f1/(f1-f2),y                      ,x                     ,y+f1/(f1-f3) /)
                    case(7)
                        nLines = nLines + 1
                        line(1:4,nLines) = (/ x+f3/(f3-f4),y+1.0d0                ,x+1.0d0               ,y+f2/(f2-f4) /)                       
                    case(8)
                        nLines = nLines + 1
                        line(1:4,nLines) = (/ x+1.0d0                ,y+f2/(f2-f4),x+f3/(f3-f4),y+1.0d0                 /)
                    case(9)
                        nLines = nLines + 1
                        line(1:4,nLines) = (/ x                      ,y+f1/(f1-f3),x+f3/(f3-f4),y+1.0d0                 /)
                        nLines = nLines + 1
                        line(1:4,nLines) = (/ x+1.0d0                ,y+f2/(f2-f4),x+f1/(f1-f2),y                       /)
                    case(10)
                        nLines = nLines + 1
                        line(1:4,nLines) = (/ x+f1/(f1-f2),y                      ,x+f3/(f3-f4),y+1.0d0                 /)
                    case(11)
                        nLines = nLines + 1
                        line(1:4,nLines) = (/ x                      ,y+f1/(f1-f3),x+f3/(f3-f4),y+1.0d0                 /)
                    case(12)
                        nLines = nLines + 1
                        line(1:4,nLines) = (/ x+1.0d0                ,y+f2/(f2-f4),x                      ,y+f1/(f1-f3) /)
                    case(13)
                        nLines = nLines + 1
                        line(1:4,nLines) = (/ x+1.0d0                ,y+f2/(f2-f4),x+f1/(f1-f2),y                       /)
                    case(14)
                        nLines = nLines + 1
                        line(1:4,nLines) = (/ x+f1/(f1-f2),y                      ,x                      ,y+f1/(f1-f3) /)
                    case(16)
                        nLines = nLines + 1
                        line(1:4,nLines) = (/ x+f3/(f3-f4),y+1.0d0                ,x                      ,y+f1/(f1-f3) /)
                        nLines = nLines + 1
                        line(1:4,nLines) = (/ x+f1/(f1-f2),y                      ,x+1.0d0                ,y+f2/(f2-f4) /)
                    case(17)
                        nLines = nLines + 1
                        line(1:4,nLines) = (/ x                      ,y+f1/(f1-f3),x+f1/(f1-f2),y                       /)
                        nLines = nLines + 1
                        line(1:4,nLines) = (/ x+1.0d0                ,y+f2/(f2-f4),x+f3/(f3-f4),y+1.0d0                 /)
                        
                end select                  
                return
            end subroutine addLines         
            
                        
        end subroutine marchingSquares

    !---
    
    
        subroutine marchingSquaresArea( img,isolevel,area )
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
            real(kind=real64),dimension(0:,0:),intent(in)       ::      img
            real(kind=real64),intent(in)                        ::      isolevel
            real(kind=real64),intent(out)                       ::      area
            integer             ::      nx,ny,ix,iy
            real(kind=real64)   ::      f1,f2,f3,f4 
            
            nx = size(img,dim=1)
            ny = size(img,dim=2)
            area = 0
            do iy = 0,ny-2
                f1 = img(nx-1,iy) - isolevel ; f3 = img(nx-1,iy+1) - isolevel 
                f2 = img(0,iy) - isolevel ; f4 = img(0,iy+1) - isolevel 
                area = area + addArea( getCode( f1,f2,f3,f4 ) )             
                do ix = 0,nx-2
                    f1 = f2 ; f3 = f4
                    f2 = img(ix+1,iy) - isolevel ; f4 = img(ix+1,iy+1) - isolevel                   
                    area = area + addArea( getCode( f1,f2,f3,f4 ) )
                end do
            end do
            f1 = img(nx-1,ny-1) - isolevel ; f3 = img(nx-1,0) - isolevel 
            f2 = img(0,ny-1) - isolevel ; f4 = img(0,0) - isolevel 
            area = area + addArea( getCode( f1,f2,f3,f4 ) )
            do ix = 0,nx-2
                f1 = f2 ; f3 = f4
                f2 = img(ix+1,ny-1) - isolevel ; f4 = img(ix+1,0) - isolevel
                area = area + addArea( getCode( f1,f2,f3,f4 ) )
            end do
            
            return

        contains
    !---^^^^^^^^
    
            function addArea(c) result(a)
        !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
                integer,intent(in)                  ::      c
                real(kind=real64)                   ::      a
                a = 0.0d0
                select case(c)
                    case(1)
                        a = f1*f1/( 2*(f1-f3)*(f1-f2) )
                    case(2)                                                                                                                       
                        a = f2*f2/( 2*(f2-f4)*(f2-f1) )
                    case(3)
                        a = ( f1*(f2-f4) + f2*(f1-f3) )/(2*(f2-f4)*(f1-f3))
                    case(4)                                                                    
                        a = f3*f3/(2*(f3-f1)*(f3-f4)) 
                    case(5)
                        a = ( f1*(f3-f4) + f3*(f1-f2) )/( 2*(f1-f2)*(f3-f4) )
                    case(6)
                        a = 1 - f1*f1/( 2*(f1-f3)*(f1-f2) ) - f4*f4/(2*(f4-f2)*(f4-f3)) 
                    case(7)
                        a = 1 - f4*f4/(2*(f4-f2)*(f4-f3)) 
                    case(8)
                        a = f4*f4/(2*(f4-f2)*(f4-f3)) 
                    case(9)
                        a = 1 - f1*f4/( 2*(f4-f2)*(f1-f2) ) - f1*f4/(2*(f1-f3)*(f4-f3)) 
                    case(10)
                        a = ( f4*(f2-f1) + f2*(f4-f3) )/(2*(f2-f1)*(f4-f3))
                    case(11)
                        a = 1 - f3*f3/(2*(f3-f1)*(f3-f4)) 
                    case(12)
                        a = (f3*(f4-f2) + f4*(f3-f1))/(2*(f3-f1)*(f4-f2))
                    case(13)
                        a = 1 - f2*f2/( 2*(f2-f4)*(f2-f1) )
                    case(14)
                        a = 1 - f1*f1/( 2*(f1-f3)*(f1-f2) )
                    case(15)
                        a = 1
                    case(16)
                        a = f2*f2/( 2*(f2-f4)*(f2-f1) ) + f3*f3/(2*(f3-f1)*(f3-f4)) 
                    case(17)
                        a = f1*f1/( 2*(f1-f3)*(f1-f2) ) + f4*f4/(2*(f4-f2)*(f4-f3)) 
                end select  
                 
                                
                return
            end function addArea            
            
                        
        end subroutine marchingSquaresArea

        subroutine contiguousLines( nLines,line , nPaths,pathLen,path ,Nx,Ny)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      break up the lines into contiguous segments
            integer,intent(in)                  ::      nLines
            real(kind=real64),dimension(:,:),intent(in)   ::        line            !    (1:4,1:nLines) = xFrom,yFrom, xTo,yTo
            integer,intent(out)                 ::      nPaths
            integer,dimension(:),allocatable,intent(out)    ::      pathLen                     !   nPaths
            real(kind=real64),dimension(:,:,:),allocatable,intent(out)  ::        path          !    (1:2,0:pathLen,nPaths) = (x,y). Note (:,0) = (:,pathLen)
            
            integer,intent(in)          ::      Nx,Ny       !   link cell subdivisions
            
            
            integer                     ::      ii,jj,kk,ll,i1,i2,i3,mm,nn
            integer,dimension(nLines)   ::      indx
            real(kind=real64)           ::      dij
            logical     ::      done
            real(kind=real64),dimension(3,3)    ::      aa = reshape( (/1.0,0.0,0.0 , 0.0,1.0,0.0 , 0.0,0.0,4.0/),(/3,3/) )
            type(LinkCell3d)        ::      lc3d
            real(kind=real64),dimension(3)      ::      xx
            real(kind=real64),dimension(:),allocatable      ::      dr
            integer,dimension(:),allocatable      ::      id
            logical         ::      ok
!            logical         ::      dbg
            
!             print *,"Lib_MarchingSquares::contiguousLines info - nLines = ",nLines
           ! print *,"Lib_MarchingSquares::contiguousLines info - creating contiguous paths from ",nLines," lines."
        !---    construct a linkcell list and add the "from" points        
            mm = estimateMaxCountPerCell(Nx*Ny,nLines)
            lc3d = LinkCell3D_ctor(aa,mm,Nx,Ny,1)
            xx(3) = 0.5d0
            
        
            do ii = 1,nLines
                xx(1:2) = line(1:2,ii)
                call add( lc3d,ii,xx)                        
            end do
            
            nn = getnNeighMax(lc3d)
            allocate(dr(nn))
            allocate(id(nn))
            
            
           !     call report( lc3d )
           !     print *,"minmax ",minval(line(1,1:nLines)),maxval(line(1,1:nLines))
           !     print *,"minmax ",minval(line(2,1:nLines)),maxval(line(2,1:nLines))
           !     print *,"minmax ",minval(line(3,1:nLines)),maxval(line(3,1:nLines))
           !     print *,"minmax ",minval(line(4,1:nLines)),maxval(line(4,1:nLines))
           !     print *,"number of lines ",nLines
           !  
           !     
           !     print *,"line(1:4,1) ",line(1:4,1)
           !     print *,"line(1:4,2) ",line(1:4,2)
           !     print *,"line(1:4,3) ",line(1:4,3)
           !     print *,"line(1:4,4) ",line(1:4,4)
                
                
            indx = 0                    !   line index is set to zero if it isn't connected yet, otherwise it's set to the path number
            nPaths = 0
            do mm = 1,nLines            !   loop until done, but not forever
                done = .true.
                
                
                do ii = 1,nLines        
                    i1 = indx(ii)
                    xx(1:2) = line(3:4,ii)                                  !   3:4 is a "to" point   
                    call neighbourList( lc3d,xx,1.1d0, nn,id,dr )         
                                         
                    do kk = 1,nn
                        
                        if (dr(kk)>1.0d-12) cycle       !   this point is too far away to be connected.
                                                
                        jj = id(kk)                 !   line number
                        if (ii==jj) cycle           !   can this actually happen???
                        
                        !   line i "to" point connects with line j "from"
                        i2 =  indx(jj)
                        if ( (i1+i2>0).and.(i1==i2) ) cycle     !   done already
                        
                      !  print *,"    connects ",i1,ii,line(:,ii)," to ",i2,jj,line(:,jj)," d ",dr(kk)
                                        
                        if (i1+i2 == 0) then
                            !   neither i nor j have been put on a path
                            nPaths = nPaths + 1
                            indx(ii) = nPaths
                            indx(jj) = nPaths
                            
                        else if (i1 == 0) then
                            !   j is on a path i2 , i is not
                            indx(ii) = i2 
                        else if (indx(jj) == 0) then
                            !   i is on a path i1, j is not
                            indx(jj) = i1 
                        else 
                            !   i and j on different paths atm. combine
                            i3 = min(i1,i2)
                            do ll = 1,nLines
                                if ( (indx(ll)-i1)*(indx(ll)-i2)==0 ) indx(ll) = i3     !   if line is i1 or i2
                            end do
                        end if
                        
                        done = .false.
                        exit
                        
                        
                        
                    end do
                end do
                if (done) exit                
            end do
            deallocate(dr)
            deallocate(id)
            call delete(lc3d)
 
            
        !---    some paths will have been cut completely by this process. Reorder sequentially
            nPaths = 0
            do i1 = 1,nLines  
                if (indx(i1) > nPaths) then
                    nPaths = nPaths + 1
                    where (indx == indx(i1))
                        indx = nPaths
                    end where
                end if
            end do            
          !  print *,"Lib_MarchingSquares::contiguousLines info - nPaths = ",nPaths  
            
        !---    can now put the lines into the paths. Find the longest path
            allocate(pathLen(nPaths))
            kk = 0 
            do i1 = 1,nPaths  
                kk = max(kk,count(indx==i1))
            end do
            allocate(path(2,0:kk,nPaths))            
           ! print *,"Lib_MarchingSquares::contiguousLines info - longest pathLen = ",kk
           ! 
          !  print *,"paths for lines 1:10 ",indx(1:10)
          !  print *,"paths for lines nLines-10: ",indx(nLines-10:)
            
            
            pathLen = 0
            path = 0
            do ii = 1,nLines
                i1 = indx(ii)       
                if (i1 == 0) then
                    print *,"did not ascribe ",i1,ii,line(1:4,ii)," to a path?"
                    cycle
                end if
                
                if (pathLen(i1)==0) then
                    path(1:2,0,i1) = line(1:2,ii)               !   path 0 is "from" of first line
                    path(1:2,1,i1) = line(3:4,ii)               !   path 1 is "to" of first line
                    indx(ii) = 0                                !   remove this line, it's now in a path
                end if
                pathLen(i1) = pathLen(i1) + 1
                
                
            end do
            
           ! print *,"pathlen (1:10) ", pathLen(1:10)
            
            
            
!             print *,"Lib_MarchingSquares::contiguousLines info - pathLen = ",pathLen
        !---    now order the paths so they make a circuit
            do i1 = 1,nPaths
                
                do kk = 2,pathLen(i1) 
                    ok = .false.
                    do ii = 1,nLines
                        if (indx(ii) /= i1) cycle
                        dij = ( path(1,kk-1,i1) - line(1,ii) )**2 + ( path(2,kk-1,i1) - line(2,ii) )**2     !   compare "from" point of next line to end of path
                        if (dij < 1.0d-12) then
                            path(1:2,kk,i1) = line(3:4,ii)
                            indx(ii) = 0                                !   remove this line, it's now in a path
                            ok = .true.
                            exit
                        end if
                    end do
                    if (.not. ok) then
                        dij = ( path(1,kk-1,i1) - path(1,0,i1) )**2 + ( path(2,kk-1,i1) - path(2,0,i1) )**2
                        if (dij < 1) then
                            pathlen(i1) = kk-1
                            exit
                        else                    
                            print *,"Lib_MarchingSquares::contiguousLines error - failed to close curve after ",kk," expected ",pathlen(i1)," start ",path(1:2,0,i1)," end ",path(1:2,kk-1,i1)," |d| ",sqrt(dij)
                            stop
                        end if
                    end if
                    !print *," path ",i1," length ",kk," expected path length ",pathlen(i1)," ok ",ok," last point ",path(1:2,kk,i1)," line ",line(1:4,ii)
                end do
            end do
             
             
                    
          !  do i1 = 1,min(10,nPaths)
          !      do kk = 0,min(10,pathLen(i1))
          !          print *,"path ",i1," point ",kk,path(1:2,kk,i1)
          !      end do
          !      do kk = max(10,pathLen(i1)-10),pathLen(i1)
          !          print *,"path ",i1," point ",kk,path(1:2,kk,i1)
          !      end do
          !  end do
!            stop                        
            return
        end subroutine contiguousLines
             
        
        subroutine areaAndLengthAndCurvContiguousLine(path,pathlen,area,length,curv,curv2)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      given a path 0:pathlen
    !*      compute its area and length and integrated curvature and integrated curvature squared
            real(kind=real64),dimension(:,0:),intent(in)    ::      path
            integer,intent(in)                              ::      pathlen
            real(kind=real64),intent(out)                   ::      area,length,curv,curv2
            
            integer                             ::      kk 
            real(kind=real64),dimension(2)      ::      p0,p1,p2
            real(kind=real64)                   ::      iR,ll,ll0,aa
            
            area = 0
            length = 0
            curv = 0
            curv2 = 0
            
            
            
            p1 = path(1:2,pathlen)        !   close circuit by using "before" point of first segment to be the last point in the path
            p2 = path(1:2,0)
            ll = norm2( p2-p1 ) 
            
            
            do kk = 1,pathlen
                ll0 = ll
                p0 = p1
                p1 = p2
                p2 = path(1:2,kk)
                ll = norm2( p2-p1 )    
                length = length + ll   
                aa =  p1(2)*p2(1) - p1(1)*p2(2)        
                area = area + aa
!                print *,"area and len ",kk,p0,p1,p2,aa,ll
                iR = MengerCurvature(p0,p1,p2) 
                 
                
                !if (abs(iR)*(ll+ll0)>1) print *,"areaAndLengthAndCurvContiguousLine warning ",kk,pathlen," l ",ll,ll0,iR," p ",p0,p1,p2
                ! if (iR*(ll+ll0)>sqrt(4/3.0d0)) iR = sqrt(4/3.0d0)/(ll+ll0)
                curv = curv + iR*(ll+ll0)
                curv2 = curv2 + iR*iR*(ll+ll0)
!                 print *,"area and len ",kk,p0,p1,p2,aa,ll,iR,length,area,curv
            end do
            area = area / 2
            curv = curv / 2
            curv2 = curv2 / 2
            
!             stop
            
            return
        end subroutine areaAndLengthAndCurvContiguousLine
        
        
        pure function MengerCurvature(u,v,w) result(c)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      given three points in 2D, u,v,w
    !*      compute the circle whose perimeter intersects all three points
    !*      return the curvature c = 1/R
    !*      using the Menger formula
    !*          c = 1/R = 4 A / ( |u-v||v-w||w-u| )
            real(kind=real64),dimension(2),intent(in)           ::      u,v,w
            real(kind=real64)                                   ::      c
            
            real(kind=real64)       ::      twoaa,uv2,vw2,wu2
            
            twoaa = ( u(1)-v(1) )*( w(2)-v(2) ) - ( u(2)-v(2) )*( w(1)-v(1) )       !   2 x signed area
            
            uv2 = ( u(1)-v(1) )*( u(1)-v(1) ) + ( u(2)-v(2) )*( u(2)-v(2) )         !   squared distance between points
            vw2 = ( w(1)-v(1) )*( w(1)-v(1) ) + ( w(2)-v(2) )*( w(2)-v(2) )
            wu2 = ( w(1)-u(1) )*( w(1)-u(1) ) + ( w(2)-u(2) )*( w(2)-u(2) )
            
            c = uv2*vw2*wu2
            if (c > 0) then
                c = 2*twoaa / sqrt( c )
            !   else coincident points. Area must be zero too. Curvature undefined, but for our purposes probably ok to return zero.                
            end if
            
            return
        end function MengerCurvature

    end module Lib_MarchingSquares
     