
    module lib_linint
!---^^^^^^^^^^^^^^^^^
!*      simple linear interpolation functions
!*      Daniel Mason
!*      (c) UKAEA Dec 2023

        use iso_fortran_env
        implicit none
        private
        
        
        public      ::      linint
        public      ::      bilinint
        public      ::      trilinint
        public      ::      quadlinint


        public      ::      dlinint                 !   derivative of function
        public      ::      dbilinint
        public      ::      dtrilinint
        public      ::      dquadlinint

        
        interface   linint
            module procedure    linint64
            module procedure    linintv32
            module procedure    linintv64
            module procedure    lininta32
            module procedure    lininta64
        end interface
                
        interface   bilinint
            module procedure    bilinint64
            module procedure    bilinintv32
            module procedure    bilinintv64
            module procedure    bilininta32
            module procedure    bilininta64
        end interface
        
        interface   trilinint
            module procedure    trilinintv32
            module procedure    trilinintv64
            module procedure    trilininta32
            module procedure    trilininta64
        end interface
        
        interface   quadlinint
            module procedure    quadlinintv32
            module procedure    quadlinintv64
            module procedure    quadlininta32
            module procedure    quadlininta64
        end interface
        
        
        
        interface   dlinint
            module procedure    dlinintv32
            module procedure    dlinintv64
        end interface
                
        interface   dbilinint
            module procedure    dbilinintv32
            module procedure    dbilinintv64
        end interface
        
        interface   dtrilinint
            module procedure    dtrilinintv32
            module procedure    dtrilinintv64
        end interface
        
        interface   dquadlinint
            module procedure    dquadlinintv32
            module procedure    dquadlinintv64
        end interface
        
        
    contains
!---^^^^^^^^

 
        pure function linintv32( y,i,a ) result(f)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      linearly interpolate y a fraction a of the distance between knots i and i+1
            real(kind=real32),dimension(:,0:),intent(in)        ::      y           !   (0:n-1)
            integer,intent(in)                                  ::      i
            real(Kind=real64),intent(in)                        ::      a
            real(kind=real32),dimension(size(y,dim=1))          ::      f
            integer             ::      nn
            nn = size(y,dim=2)
            
            if (i < 0) then
                f = y(:,0)
            else if (i>=nn-1) then
                f = y(:,nn-1)
            else
                f = real( y(:,i)*(1-a) + y(:,i+1)*a , kind=real32 )
            end if
            !write(*,fmt='(a,100g12.4)') "linintv32    ", y(:,i),y(:,i+1),a,f
            
            return
        end function linintv32
        
        
        pure function bilinintv32( y,i,j,a,b ) result(f)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ 
    !*      bilinearly interpolate y a fraction a of the distance between knots i and i+1
    !*      and a fraction b between knots j and j+1
            real(kind=real32),dimension(:,0:,0:),intent(in)     ::      y           !   (0:nx-1,0:ny-1)
            integer,intent(in)                                  ::      i,j
            real(kind=real64),intent(in)                        ::      a,b
            real(kind=real32),dimension(size(y,dim=1))          ::      f
            
            real(kind=real32),dimension(size(y,dim=1))       ::      linj,linjp1
            integer                 ::      ny
            ny = size(y,dim=3)
            
            if (j < 0) then                 
                linj   = linint( y(:,:,0),i,a )
                linjp1 = linint( y(:,:,0),i,a )
            else if (j >= ny-1) then
                linj   = linint( y(:,:,ny-1),i,a )
                linjp1 = linint( y(:,:,ny-1),i,a )
            else
                linj   = linint( y(:,:,j)  ,i,a )
                linjp1 = linint( y(:,:,j+1),i,a )
            end if
                                   
            
            f = real( linj*(1-b) + linjp1*b , kind=real32 )
            !write(*,fmt='(a,100g12.4)') "bilinintv32  ",linj,linjp1,b,f
            
            return
        end function bilinintv32
        
        pure function trilinintv32( y,i,j,k,a,b,c ) result(f) 
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ 
    !*      trilinearly interpolate y a fraction a of the distance between knots i and i+1
    !*      and a fraction b between knots j and j+1
    !*      and a fraction c between k and k+1
            real(kind=real32),dimension(:,0:,0:,0:),intent(in)    ::      y           !   (0:nx-1,0:ny-1,0:nz-1)
            integer,intent(in)                                  ::      i,j,k
            real(kind=real64),intent(in)                        ::      a,b,c
            real(kind=real32),dimension(size(y,dim=1))          ::      f
            real(kind=real32),dimension(size(y,dim=1))       ::      link,linkp1
            integer                 ::      nz
            nz = size(y,dim=4)
            
            if (k < 0) then                 
                link   = bilinint( y(:,:,:,0),i,j,a,b )
                linkp1 = bilinint( y(:,:,:,0),i,j,a,b )
            else if (k >= nz-1) then
                link   = bilinint( y(:,:,:,nz-1),i,j,a,b )
                linkp1 = bilinint( y(:,:,:,nz-1),i,j,a,b )
            else
                link   = bilinint( y(:,:,:,k)  ,i,j,a,b )
                linkp1 = bilinint( y(:,:,:,k+1),i,j,a,b )
            end if
                      
                         
            f = real( link*(1-c) + linkp1*c , kind=real32 )
            !write(*,fmt='(a,100g12.4)') "trilinintv32 ",link,linkp1,c,f
            
            return
        end function trilinintv32
        
        
        pure function quadlinintv32( y,i,j,k,l,a,b,c,d ) result(f) 
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ 
    !*      quadrilinearly interpolate y a fraction a of the distance between knots i and i+1
    !*      and a fraction b between knots j and j+1
    !*      and a fraction c between k and k+1
    !*      and a fraction d between l and l+1
            real(kind=real32),dimension(:,0:,0:,0:,0:),intent(in) ::      y           !   (0:nx-1,0:ny-1,0:nz-1,0:nw-1)
            integer,intent(in)                                  ::      i,j,k,l
            real(kind=real64),intent(in)                        ::      a,b,c,d
            real(kind=real32),dimension(size(y,dim=1))          ::      f
            
            real(kind=real32),dimension(size(y,dim=1))       ::      linl,linlp1
            integer                 ::      nw
            nw = size(y,dim=4)
            
            if (l < 0) then                 
                linl   = trilinint( y(:,:,:,:,0),i,j,k,a,b,c )
                linlp1 = trilinint( y(:,:,:,:,0),i,j,k,a,b,c )
            else if (l >= nw-1) then
                linl   = trilinint( y(:,:,:,:,nw-1),i,j,k,a,b,c )
                linlp1 = trilinint( y(:,:,:,:,nw-1),i,j,k,a,b,c )
            else
                linl   = trilinint( y(:,:,:,:,l)  ,i,j,k,a,b,c )
                linlp1 = trilinint( y(:,:,:,:,l+1),i,j,k,a,b,c )
            end if
                                   
            f = real( linl*(1-d) + linlp1*d , kind=real32 )
            
            return
        end function quadlinintv32
        
        
 
        pure function dlinintv32( y,h,i ) result(df)    
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      derivative of y at a fraction a of the distance between knots i and i+1
    !*      assuming even knot spacing h.
            real(kind=real32),dimension(:,0:),intent(in)        ::      y           !   (0:n-1)
            real(kind=real64),intent(in)                        ::      h
            integer,intent(in)                                  ::      i
            real(kind=real32),dimension(size(y,dim=1))          ::      df
            integer             ::      nn
            nn = size(y,dim=2)
            
            if (i < 0) then
                df = 0
            else if (i>=nn-1) then
                df = 0
            else
                df = real( ( y(:,i+1) - y(:,i) )/h , kind=real32 )
            end if
            
            return
        end function dlinintv32
        
        
        pure function dbilinintv32( y,h1,h2,i,j,a,b ) result(df)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      derivative of y at a fraction a of the distance between knots i and i+1
    !*      assuming even knot spacing h1 in x and h2 in y
            real(kind=real32),dimension(:,0:,0:),intent(in)     ::      y           !   (0:n-1)
            real(kind=real64),intent(in)                        ::      h1,h2
            integer,intent(in)                                  ::      i,j
            real(kind=real64),intent(in)                        ::      a,b
            real(kind=real32),dimension(size(y,dim=1),2)        ::      df
            
            real(kind=real64),dimension(size(y,dim=1),4)  ::      yint
            
            yint(:,1) = bilinint(y,i,j  ,0.0d0,b)
            yint(:,2) = bilinint(y,i+1,j,0.0d0,b)
            yint(:,3) = bilinint(y,i,j  ,a,0.0d0) 
            yint(:,4) = bilinint(y,i,j+1,a,0.0d0) 
             
            df(:,1) = real( ( yint(:,2) - yint(:,1) )/h1 , kind=real32 )
            df(:,2) = real( ( yint(:,4) - yint(:,3) )/h2 , kind=real32 )
            
            return
        end function dbilinintv32
        
        
        pure function dtrilinintv32( y,h1,h2,h3,i,j,k,a,b,c ) result(df)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      derivative of y at a fraction a of the distance between knots i and i+1
    !*      assuming even knot spacing h1 in x and h2 in y
            real(kind=real32),dimension(:,0:,0:,0:),intent(in)  ::      y           !   (0:n-1)
            real(kind=real64),intent(in)                        ::      h1,h2,h3
            integer,intent(in)                                  ::      i,j,k
            real(kind=real64),intent(in)                        ::      a,b,c
            real(kind=real32),dimension(size(y,dim=1),3)        ::      df
            
            real(kind=real64),dimension(size(y,dim=1),6)        ::      yint            
            
            yint(:,1) = trilinint(y,i,j  ,k,0.0d0,b,c)
            yint(:,2) = trilinint(y,i+1,j,k,0.0d0,b,c)
            yint(:,3) = trilinint(y,i,j  ,k,a,0.0d0,c) 
            yint(:,4) = trilinint(y,i,j+1,k,a,0.0d0,c) 
            yint(:,5) = trilinint(y,i,j,k  ,a,b,0.0d0) 
            yint(:,6) = trilinint(y,i,j,k+1,a,b,0.0d0) 
            
            df(:,1) = real( ( yint(:,2) - yint(:,1) )/h1 , kind=real32 )
            df(:,2) = real( ( yint(:,4) - yint(:,3) )/h2 , kind=real32 )
            df(:,3) = real( ( yint(:,6) - yint(:,5) )/h3 , kind=real32 )
            
            return
        end function dtrilinintv32
        
        
        pure function dquadlinintv32( y,h1,h2,h3,h4,i,j,k,l,a,b,c,d ) result(df)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      derivative of y at a fraction a of the distance between knots i and i+1
    !*      assuming even knot spacing h1 in x and h2 in y
            real(kind=real32),dimension(:,0:,0:,0:,0:),intent(in)     ::      y           !   (0:n-1)
            real(kind=real64),intent(in)                        ::      h1,h2,h3,h4
            integer,intent(in)                                  ::      i,j,k,l
            real(kind=real64),intent(in)                        ::      a,b,c,d
            real(kind=real32),dimension(size(y,dim=1),4)        ::      df
            
            real(kind=real64),dimension(size(y,dim=1),8)  ::      yint            
            
            yint(:,1) = quadlinint(y,i,j  ,k,l,0.0d0,b,c,d)
            yint(:,2) = quadlinint(y,i+1,j,k,l,0.0d0,b,c,d)
            yint(:,3) = quadlinint(y,i,j  ,k,l,a,0.0d0,c,d) 
            yint(:,4) = quadlinint(y,i,j+1,k,l,a,0.0d0,c,d) 
            yint(:,5) = quadlinint(y,i,j,k  ,l,a,b,0.0d0,d) 
            yint(:,6) = quadlinint(y,i,j,k+1,l,a,b,0.0d0,d) 
            yint(:,7) = quadlinint(y,i,j,k,l  ,a,b,c,0.0d0) 
            yint(:,8) = quadlinint(y,i,j,k,l+1,a,b,c,0.0d0) 
             
            df(:,1) = real( ( yint(:,2) - yint(:,1) )/h1 , kind=real32 )
            df(:,2) = real( ( yint(:,4) - yint(:,3) )/h2 , kind=real32 )
            df(:,3) = real( ( yint(:,6) - yint(:,5) )/h3 , kind=real32 )
            df(:,4) = real( ( yint(:,8) - yint(:,7) )/h4 , kind=real32 )
            
            return
        end function dquadlinintv32
        
        
        pure function dlinintv64( y,h,i ) result(df)       
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      derivative of y at a fraction a of the distance between knots i and i+1
    !*      assuming even knot spacing h.
            real(kind=real64),dimension(:,0:),intent(in)        ::      y           !   (0:n-1)
            real(kind=real64),intent(in)                        ::      h
            integer,intent(in)                                  ::      i
            real(kind=real64),dimension(size(y,dim=1))          ::      df
            integer             ::      nx
            nx = size(y,dim=2)
            
            if (i < 0) then
                df = 0
            else if (i>=nx-1) then
                df = 0
            else
                df =( y(:,i+1) - y(:,i) )/h 
            end if
            
            return
        end function dlinintv64
        
        
        pure function dbilinintv64( y,h1,h2,i,j,a,b ) result(df)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      derivative of y at a fraction a of the distance between knots i and i+1
    !*      assuming even knot spacing h1 in x and h2 in y
            real(kind=real64),dimension(:,0:,0:),intent(in)     ::      y           !   (0:n-1)
            real(kind=real64),intent(in)                        ::      h1,h2
            integer,intent(in)                                  ::      i,j
            real(kind=real64),intent(in)                        ::      a,b
            real(kind=real64),dimension(size(y,dim=1),2)        ::      df
            
            real(kind=real64),dimension(size(y,dim=1),4)  ::      yint
            
            yint(:,1) = bilinint(y,i,j  ,0.0d0,b)
            yint(:,2) = bilinint(y,i+1,j,0.0d0,b)
            yint(:,3) = bilinint(y,i,j  ,a,0.0d0) 
            yint(:,4) = bilinint(y,i,j+1,a,0.0d0) 
             
            df(:,1) = ( yint(:,2) - yint(:,1) )/h1 
            df(:,2) = ( yint(:,4) - yint(:,3) )/h2 
            
            return
        end function dbilinintv64
        
        
        pure function dtrilinintv64( y,h1,h2,h3,i,j,k,a,b,c ) result(df)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      derivative of y at a fraction a of the distance between knots i and i+1
    !*      assuming even knot spacing h1 in x and h2 in y
            real(kind=real64),dimension(:,0:,0:,0:),intent(in)  ::      y           !   (0:n-1)
            real(kind=real64),intent(in)                        ::      h1,h2,h3
            integer,intent(in)                                  ::      i,j,k
            real(kind=real64),intent(in)                        ::      a,b,c
            real(kind=real64),dimension(size(y,dim=1),3)        ::      df
            
            real(kind=real64),dimension(size(y,dim=1),6)        ::      yint            
            
            yint(:,1) = trilinint(y,i,j  ,k,0.0d0,b,c)
            yint(:,2) = trilinint(y,i+1,j,k,0.0d0,b,c)
            yint(:,3) = trilinint(y,i,j  ,k,a,0.0d0,c) 
            yint(:,4) = trilinint(y,i,j+1,k,a,0.0d0,c) 
            yint(:,5) = trilinint(y,i,j,k  ,a,b,0.0d0) 
            yint(:,6) = trilinint(y,i,j,k+1,a,b,0.0d0) 
            
            df(:,1) = ( yint(:,2) - yint(:,1) )/h1 
            df(:,2) = ( yint(:,4) - yint(:,3) )/h2 
            df(:,3) = ( yint(:,6) - yint(:,5) )/h3 
            
            return
        end function dtrilinintv64
        
        
        pure function dquadlinintv64( y,h1,h2,h3,h4,i,j,k,l,a,b,c,d ) result(df)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      derivative of y at a fraction a of the distance between knots i and i+1
    !*      assuming even knot spacing h1 in x and h2 in y
            real(kind=real64),dimension(:,0:,0:,0:,0:),intent(in)     ::      y           !   (0:n-1)
            real(kind=real64),intent(in)                        ::      h1,h2,h3,h4
            integer,intent(in)                                  ::      i,j,k,l
            real(kind=real64),intent(in)                        ::      a,b,c,d
            real(kind=real64),dimension(size(y,dim=1),4)        ::      df
            
            real(kind=real64),dimension(size(y,dim=1),8)  ::      yint            
            
            yint(:,1) = quadlinint(y,i,j  ,k,l,0.0d0,b,c,d)
            yint(:,2) = quadlinint(y,i+1,j,k,l,0.0d0,b,c,d)
            yint(:,3) = quadlinint(y,i,j  ,k,l,a,0.0d0,c,d) 
            yint(:,4) = quadlinint(y,i,j+1,k,l,a,0.0d0,c,d) 
            yint(:,5) = quadlinint(y,i,j,k  ,l,a,b,0.0d0,d) 
            yint(:,6) = quadlinint(y,i,j,k+1,l,a,b,0.0d0,d) 
            yint(:,7) = quadlinint(y,i,j,k,l  ,a,b,c,0.0d0) 
            yint(:,8) = quadlinint(y,i,j,k,l+1,a,b,c,0.0d0) 
             
            df(:,1) = ( yint(:,2) - yint(:,1) )/h1 
            df(:,2) = ( yint(:,4) - yint(:,3) )/h2 
            df(:,3) = ( yint(:,6) - yint(:,5) )/h3 
            df(:,4) = ( yint(:,8) - yint(:,7) )/h4 
            
            return
        end function dquadlinintv64
        
        
        
        
!-------



 
        pure function linint64( y,i,a ) result(f)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      linearly interpolate y a fraction a of the distance between knots i and i+1
            real(kind=real64),dimension(0:),intent(in)          ::      y           !   (0:nx-1)
            integer,intent(in)                                  ::      i
            real(Kind=real64),intent(in)                        ::      a
            real(kind=real64)                                   ::      f
            integer             ::      nx
            nx = size(y,dim=1)
            
            if (i < 0) then
                f = y(0)
            else if (i>=nx-1) then
                f = y(nx-1)
            else
                f = y(i)*(1-a) + y(i+1)*a
            end if
            
            return
        end function linint64
        
        
        pure function bilinint64( y,i,j,a,b ) result(f)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ 
    !*      bilinearly interpolate y a fraction a of the distance between knots i and i+1
    !*      and a fraction b between knots j and j+1
            real(kind=real64),dimension(0:,0:),intent(in)       ::      y           !   (0:nx-1,0:ny-1)
            integer,intent(in)                                  ::      i,j
            real(kind=real64),intent(in)                        ::      a,b
            real(kind=real64)                                   ::      f
            
            real(kind=real64)       ::      linj,linjp1
            integer                 ::      ny
            ny = size(y,dim=2)
            
            if (j < 0) then                 
                linj   = linint( y(:,0),i,a )
                linjp1 = linint( y(:,0),i,a )
            else if (j >= ny-1) then
                linj   = linint( y(:,ny-1),i,a )
                linjp1 = linint( y(:,ny-1),i,a )
            else
                linj   = linint( y(:,j)  ,i,a )
                linjp1 = linint( y(:,j+1),i,a )
            end if
                                   
            
            f = linj*(1-b) + linjp1*b
            
            return
        end function bilinint64
                
        
        
        
        
        
        
        pure function lininta32( y,i,a ) result(f)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      linearly interpolate y a fraction a of the distance between knots i and i+1
            real(kind=real32),dimension(:,:,0:),intent(in)      ::      y           !   (0:nx-1)
            integer,intent(in)                                  ::      i
            real(Kind=real64),intent(in)                        ::      a
            real(kind=real32),dimension(size(y,dim=1),size(y,dim=2))          ::      f
            integer             ::      nx
            nx = size(y,dim=3)
            
            if (i < 0) then
                f = y(:,:,0)
            else if (i>=nx-1) then
                f = y(:,:,nx-1)
            else
                f = real( y(:,:,i)*(1-a) + y(:,:,i+1)*a , kind=real32 )
            end if
             
            return
        end function lininta32
        
        
        pure function bilininta32( y,i,j,a,b ) result(f)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ 
    !*      bilinearly interpolate y a fraction a of the distance between knots i and i+1
    !*      and a fraction b between knots j and j+1
            real(kind=real32),dimension(:,:,0:,0:),intent(in)   ::      y           !   (0:nx-1,0:ny-1)
            integer,intent(in)                                  ::      i,j
            real(kind=real64),intent(in)                        ::      a,b
            real(kind=real32),dimension(size(y,dim=1),size(y,dim=2))          ::      f
            
            real(kind=real32),dimension(size(y,dim=1),size(y,dim=2))       ::      linj,linjp1
            integer                 ::      nx,ny
            real(kind=real32)       ::      a32,b32

            ny = size(y,dim=4)
            if (j < 0) then                 
                f = lininta32(y(:,:,:,0),i,a)
            else if (j >= ny-1) then
                f = lininta32(y(:,:,:,ny-1),i,a)
            else
                linj = lininta32(y(:,:,:,j),i,a)
                linjp1 = lininta32(y(:,:,:,j+1),i,a)
                f(:,:) = real( linj(:,:)*(1.0-b) + linjp1(:,:)*b , kind=real32 )
            end if
            return

            nx = size(y,dim=3)
            ny = size(y,dim=4)
            
            a32 = real(a,kind=real32)
            b32 = real(b,kind=real32)
            if (j < 0) then                 
                if (i < 0) then
                    f = y(:,:,0,0)
                else if (i>=nx-1) then
                    f = y(:,:,nx-1,0)
                else
                    f = y(:,:,i,0)*(1.0-a32) + y(:,:,i+1,0)*a32 
                end if
            else if (j >= ny-1) then
                if (i < 0) then
                    f = y(:,:,0,ny-1)
                else if (i>=nx-1) then
                    f = y(:,:,nx-1,ny-1)
                else
                    f = y(:,:,i,ny-1)*(1.0-a32) + y(:,:,i+1,ny-1)*a32
                end if
            else
                if (i < 0) then
                    linj(:,:)   = y(:,:,0,j)
                    linjp1(:,:) = y(:,:,0,j+1)
                else if (i>=nx-1) then
                    linj(:,:)   = y(:,:,nx-1,j)
                    linjp1(:,:) = y(:,:,nx-1,j+1)
                else
                    linj(:,:)   = y(:,:,i,j)*(1.0-a32) + y(:,:,i+1,j)*a32 
                    linjp1(:,:) = y(:,:,i,j+1)*(1.0-a32) + y(:,:,i+1,j+1)*a32 
                end if
                f(:,:) = linj(:,:)*(1.0-b32) + linjp1(:,:)*b32
            end if
                                   
            
            
            !write(*,fmt='(a,100g12.4)') "bilininta32  ",linj,linjp1,b,f
            
            return
        end function bilininta32
        
        pure function trilininta32( y,i,j,k,a,b,c ) result(f) 
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ 
    !*      trilinearly interpolate y a fraction a of the distance between knots i and i+1
    !*      and a fraction b between knots j and j+1
    !*      and a fraction c between k and k+1
            real(kind=real32),dimension(:,:,0:,0:,0:),intent(in)    ::      y           !   (0:nx-1,0:ny-1,0:nz-1)
            integer,intent(in)                                  ::      i,j,k
            real(kind=real64),intent(in)                        ::      a,b,c
            real(kind=real32),dimension(size(y,dim=1),size(y,dim=2))          ::      f
            real(kind=real32),dimension(size(y,dim=1),size(y,dim=2))       ::      link,linkp1
            integer                 ::      nz
            nz = size(y,dim=5)
            
            if (k < 0) then                 
                link   = bilinint( y(:,:,:,:,0),i,j,a,b )
                linkp1 = bilinint( y(:,:,:,:,0),i,j,a,b )
            else if (k >= nz-1) then
                link   = bilinint( y(:,:,:,:,nz-1),i,j,a,b )
                linkp1 = bilinint( y(:,:,:,:,nz-1),i,j,a,b )
            else
                link   = bilinint( y(:,:,:,:,k)  ,i,j,a,b )
                linkp1 = bilinint( y(:,:,:,:,k+1),i,j,a,b )
            end if
                      
                         
            f = real( link*(1-c) + linkp1*c , kind=real32 )
            !write(*,fmt='(a,100g12.4)') "trilininta32 ",link,linkp1,c,f
            
            return
        end function trilininta32
        
        
        pure function quadlininta32( y,i,j,k,l,a,b,c,d ) result(f) 
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ 
    !*      quadrilinearly interpolate y a fraction a of the distance between knots i and i+1
    !*      and a fraction b between knots j and j+1
    !*      and a fraction c between k and k+1
    !*      and a fraction d between l and l+1
            real(kind=real32),dimension(:,:,0:,0:,0:,0:),intent(in) ::      y           !   (0:nx-1,0:ny-1,0:nz-1,0:nw-1)
            integer,intent(in)                                  ::      i,j,k,l
            real(kind=real64),intent(in)                        ::      a,b,c,d
            real(kind=real32),dimension(size(y,dim=1),size(y,dim=2))          ::      f
            
            real(kind=real32),dimension(size(y,dim=1),size(y,dim=2))       ::      linl,linlp1
            integer                 ::      nw
            nw = size(y,dim=6)
            
            if (l < 0) then                 
                linl   = trilinint( y(:,:,:,:,:,0),i,j,k,a,b,c )
                linlp1 = trilinint( y(:,:,:,:,:,0),i,j,k,a,b,c )
            else if (l >= nw-1) then
                linl   = trilinint( y(:,:,:,:,:,nw-1),i,j,k,a,b,c )
                linlp1 = trilinint( y(:,:,:,:,:,nw-1),i,j,k,a,b,c )
            else
                linl   = trilinint( y(:,:,:,:,:,l)  ,i,j,k,a,b,c )
                linlp1 = trilinint( y(:,:,:,:,:,l+1),i,j,k,a,b,c )
            end if
                                   
            f = real( linl*(1-d) + linlp1*d , kind=real32 )
            
            return
        end function quadlininta32
        
        
        
        pure function lininta64( y,i,a ) result(f)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      linearly interpolate y a fraction a of the distance between knots i and i+1
            real(kind=real64),dimension(:,:,0:),intent(in)      ::      y           !   (0:n-1)
            integer,intent(in)                                  ::      i
            real(Kind=real64),intent(in)                        ::      a
            real(kind=real64),dimension(size(y,dim=1),size(y,dim=2))          ::      f
            integer             ::      nn
            nn = size(y,dim=3)
            
            if (i < 0) then
                f = y(:,:,0)
            else if (i>=nn-1) then
                f = y(:,:,nn-1)
            else
                f =  y(:,:,i)*(1-a) + y(:,:,i+1)*a
            end if
            !write(*,fmt='(a,100g12.4)') "lininta64    ", y(:,i),y(:,i+1),a,f
            
            return
        end function lininta64
        
        
        pure function bilininta64( y,i,j,a,b ) result(f)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ 
    !*      bilinearly interpolate y a fraction a of the distance between knots i and i+1
    !*      and a fraction b between knots j and j+1
            real(kind=real64),dimension(:,:,0:,0:),intent(in)     ::      y           !   (0:nx-1,0:ny-1)
            integer,intent(in)                                  ::      i,j
            real(kind=real64),intent(in)                        ::      a,b
            real(kind=real64),dimension(size(y,dim=1),size(y,dim=2))          ::      f
            
            real(kind=real64),dimension(size(y,dim=1),size(y,dim=2))       ::      linj,linjp1
            integer                 ::      ny
            ny = size(y,dim=4)
            
            if (j < 0) then                 
                linj   = linint( y(:,:,:,0),i,a )
                linjp1 = linint( y(:,:,:,0),i,a )
            else if (j >= ny-1) then
                linj   = linint( y(:,:,:,ny-1),i,a )
                linjp1 = linint( y(:,:,:,ny-1),i,a )
            else
                linj   = linint( y(:,:,:,j)  ,i,a )
                linjp1 = linint( y(:,:,:,j+1),i,a )
            end if
                                   
            
            f = linj*(1-b) + linjp1*b
            !write(*,fmt='(a,100g12.4)') "bilininta64  ",linj,linjp1,b,f
            
            return
        end function bilininta64
        
        pure function trilininta64( y,i,j,k,a,b,c ) result(f) 
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ 
    !*      trilinearly interpolate y a fraction a of the distance between knots i and i+1
    !*      and a fraction b between knots j and j+1
    !*      and a fraction c between k and k+1
            real(kind=real64),dimension(:,:,0:,0:,0:),intent(in)    ::      y           !   (0:nx-1,0:ny-1,0:nz-1)
            integer,intent(in)                                  ::      i,j,k
            real(kind=real64),intent(in)                        ::      a,b,c
            real(kind=real64),dimension(size(y,dim=1),size(y,dim=2))          ::      f
            real(kind=real64),dimension(size(y,dim=1),size(y,dim=2))       ::      link,linkp1
            integer                 ::      nz
            nz = size(y,dim=5)
            
            if (k < 0) then                 
                link   = bilinint( y(:,:,:,:,0),i,j,a,b )
                linkp1 = bilinint( y(:,:,:,:,0),i,j,a,b )
            else if (k >= nz-1) then
                link   = bilinint( y(:,:,:,:,nz-1),i,j,a,b )
                linkp1 = bilinint( y(:,:,:,:,nz-1),i,j,a,b )
            else
                link   = bilinint( y(:,:,:,:,k)  ,i,j,a,b )
                linkp1 = bilinint( y(:,:,:,:,k+1),i,j,a,b )
            end if
                      
                         
            f = link*(1-c) + linkp1*c 
            !write(*,fmt='(a,100g12.4)') "trilininta64 ",link,linkp1,c,f
            
            return
        end function trilininta64
        
        
        pure function quadlininta64( y,i,j,k,l,a,b,c,d ) result(f) 
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ 
    !*      quadrilinearly interpolate y a fraction a of the distance between knots i and i+1
    !*      and a fraction b between knots j and j+1
    !*      and a fraction c between k and k+1
    !*      and a fraction d between l and l+1
            real(kind=real64),dimension(:,:,0:,0:,0:,0:),intent(in) ::      y           !   (0:nx-1,0:ny-1,0:nz-1,0:nw-1)
            integer,intent(in)                                  ::      i,j,k,l
            real(kind=real64),intent(in)                        ::      a,b,c,d
            real(kind=real64),dimension(size(y,dim=1),size(y,dim=2))          ::      f
            
            real(kind=real64),dimension(size(y,dim=1),size(y,dim=2))       ::      linl,linlp1
            integer                 ::      nw
            nw = size(y,dim=6)
            
            if (l < 0) then                 
                linl   = trilinint( y(:,:,:,:,:,0),i,j,k,a,b,c )
                linlp1 = trilinint( y(:,:,:,:,:,0),i,j,k,a,b,c )
            else if (l >= nw-1) then
                linl   = trilinint( y(:,:,:,:,:,nw-1),i,j,k,a,b,c )
                linlp1 = trilinint( y(:,:,:,:,:,nw-1),i,j,k,a,b,c )
            else
                linl   = trilinint( y(:,:,:,:,:,l)  ,i,j,k,a,b,c )
                linlp1 = trilinint( y(:,:,:,:,:,l+1),i,j,k,a,b,c )
            end if
                                   
            f = linl*(1-d) + linlp1*d
            
            return
        end function quadlininta64
        
        
        
        
        
        
        
        pure function linintv64( y,i,a ) result(f)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      linearly interpolate y a fraction a of the distance between knots i and i+1
            real(kind=real64),dimension(:,0:),intent(in)        ::      y           !   (0:n-1)
            integer,intent(in)                                  ::      i
            real(Kind=real64),intent(in)                        ::      a
            real(kind=real64),dimension(size(y,dim=1))          ::      f
            integer             ::      nn
            nn = size(y,dim=2)
            
            if (i < 0) then
                f = y(:,0)
            else if (i>=nn-1) then
                f = y(:,nn-1)
            else
                f =  y(:,i)*(1-a) + y(:,i+1)*a 
            end if
            !write(*,fmt='(a,100g12.4)') "linintv64    ", y(:,i),y(:,i+1),a,f
            
            return
        end function linintv64
        
        
        pure function bilinintv64( y,i,j,a,b ) result(f)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ 
    !*      bilinearly interpolate y a fraction a of the distance between knots i and i+1
    !*      and a fraction b between knots j and j+1
            real(kind=real64),dimension(:,0:,0:),intent(in)     ::      y           !   (0:nx-1,0:ny-1)
            integer,intent(in)                                  ::      i,j
            real(kind=real64),intent(in)                        ::      a,b
            real(kind=real64),dimension(size(y,dim=1))          ::      f
            
            real(kind=real64),dimension(size(y,dim=1))       ::      linj,linjp1
            integer                 ::      ny
            ny = size(y,dim=3)
            
            if (j < 0) then                 
                linj   = linint( y(:,:,0),i,a )
                linjp1 = linint( y(:,:,0),i,a )
            else if (j >= ny-1) then
                linj   = linint( y(:,:,ny-1),i,a )
                linjp1 = linint( y(:,:,ny-1),i,a )
            else
                linj   = linint( y(:,:,j)  ,i,a )
                linjp1 = linint( y(:,:,j+1),i,a )
            end if
                                   
            
            f = linj*(1-b) + linjp1*b 
            !write(*,fmt='(a,100g12.4)') "bilinintv64  ",linj,linjp1,b,f
            
            return
        end function bilinintv64
        
        pure function trilinintv64( y,i,j,k,a,b,c ) result(f) 
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ 
    !*      trilinearly interpolate y a fraction a of the distance between knots i and i+1
    !*      and a fraction b between knots j and j+1
    !*      and a fraction c between k and k+1
            real(kind=real64),dimension(:,0:,0:,0:),intent(in)    ::      y           !   (0:nx-1,0:ny-1,0:nz-1)
            integer,intent(in)                                  ::      i,j,k
            real(kind=real64),intent(in)                        ::      a,b,c
            real(kind=real64),dimension(size(y,dim=1))          ::      f
            real(kind=real64),dimension(size(y,dim=1))       ::      link,linkp1
            integer                 ::      nz
            nz = size(y,dim=4)
            
            if (k < 0) then                 
                link   = bilinint( y(:,:,:,0),i,j,a,b )
                linkp1 = bilinint( y(:,:,:,0),i,j,a,b )
            else if (k >= nz-1) then
                link   = bilinint( y(:,:,:,nz-1),i,j,a,b )
                linkp1 = bilinint( y(:,:,:,nz-1),i,j,a,b )
            else
                link   = bilinint( y(:,:,:,k)  ,i,j,a,b )
                linkp1 = bilinint( y(:,:,:,k+1),i,j,a,b )
            end if
                      
                         
            f =  link*(1-c) + linkp1*c 
            !write(*,fmt='(a,100g12.4)') "trilinintv64 ",link,linkp1,c,f
            
            return
        end function trilinintv64
        
        
        pure function quadlinintv64( y,i,j,k,l,a,b,c,d ) result(f) 
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ 
    !*      quadrilinearly interpolate y a fraction a of the distance between knots i and i+1
    !*      and a fraction b between knots j and j+1
    !*      and a fraction c between k and k+1
    !*      and a fraction d between l and l+1
            real(kind=real64),dimension(:,0:,0:,0:,0:),intent(in) ::      y           !   (0:nx-1,0:ny-1,0:nz-1,0:nw-1)
            integer,intent(in)                                  ::      i,j,k,l
            real(kind=real64),intent(in)                        ::      a,b,c,d
            real(kind=real64),dimension(size(y,dim=1))          ::      f
            
            real(kind=real64),dimension(size(y,dim=1))       ::      linl,linlp1
            integer                 ::      nw
            nw = size(y,dim=4)
            
            if (l < 0) then                 
                linl   = trilinint( y(:,:,:,:,0),i,j,k,a,b,c )
                linlp1 = trilinint( y(:,:,:,:,0),i,j,k,a,b,c )
            else if (l >= nw-1) then
                linl   = trilinint( y(:,:,:,:,nw-1),i,j,k,a,b,c )
                linlp1 = trilinint( y(:,:,:,:,nw-1),i,j,k,a,b,c )
            else
                linl   = trilinint( y(:,:,:,:,l)  ,i,j,k,a,b,c )
                linlp1 = trilinint( y(:,:,:,:,l+1),i,j,k,a,b,c )
            end if
                                   
            f = linl*(1-d) + linlp1*d
            
            return
        end function quadlinintv64
        
        
        
        
        
        
    end module Lib_linint        