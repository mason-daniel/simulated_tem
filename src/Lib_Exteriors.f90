
    module Lib_Exteriors
!---^^^^^^^^^^^^^^^^^^^^
!*      Take a 2D dataset [0:mx-1,0:my-1]
!*      and extrapolate the exterior, based on one of several models
        use Lib_GaussianBlurs
        use iso_fortran_env
        implicit none
        private
        
    !---
    
        public          ::      setExterior
    
    !---
    
    
        logical,public          ::      Exterior_dbg = .false.
        
    !---    define a 64bit floating constant "EXTERIOR_IGNORE" with a really unlikely bit pattern        
        integer(kind=int64),private,parameter               ::      BADF00D = int(z'BADF00D',kind=int64)
        real(kind=real64),public,parameter                  ::      EXTERIOR_IGNORE = transfer( (BADF00D+ishft(BADF00D,32_int64)),1.0d0 )

        integer,public,parameter        ::      EXTERIOR_ZERO = 0             !   
        integer,public,parameter        ::      EXTERIOR_PBC = 1              !   
        integer,public,parameter        ::      EXTERIOR_AVERAGE = 2          !   
        integer,public,parameter        ::      EXTERIOR_EXTRAPOLATE = 3      !     quadratic extrapolation
        integer,public,parameter        ::      EXTERIOR_BLUR = 4             !     gaussian blur with increasing radius 
        integer,public,parameter        ::      EXTERIOR_BLUR0 = 5            !     gaussian blur tend to zero
        integer,public,parameter        ::      EXTERIOR_BLURAVG = 6            !     gaussian blur tend to avg
        integer,public,parameter        ::      EXTERIOR_BLURPBC = 7            !     gaussian blur tend to pbc
        integer,public,parameter        ::      EXTERIOR_HILO = 8             !     -huge(1) if boundary point <0 , +huge(1) if boundary > 0
        integer,public,parameter        ::      EXTERIOR_LOHI = 9             !     +huge(1) if boundary point <0 , -huge(1) if boundary > 0
        
                           
    !---
        
    
        interface setExterior
            module procedure    setExterior2D
        end interface
        
    
    contains
!---^^^^^^^^

        subroutine setExterior2D( f_in,model,border,f_out, f_level )
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      apply the exterior model to fill in pixels
            real(kind=real64),dimension(0:,0:),intent(in)       ::      f_in
            integer,intent(in)                                  ::      model
            integer,intent(in)                                  ::      border
            real(kind=real64),dimension(-border:,-border:),intent(out)  ::      f_out
            real(kind=real64),intent(in),optional               ::      f_level
            
            integer             ::      nx,ny,npix
            integer             ::      ix,iy,jx,jy,kk
            real(kind=real64)   ::      fbar
            real(kind=real64),dimension(:,:,:),allocatable      ::      df      !   (1:5,0:nx-1,0:ny-1) (/ df/dx,df/dy,d2f/dx2,d2ydy2,d2f/dxdy /)
!            real(kind=real64),dimension(:,:),allocatable      ::      df2,df3      !   (1:5,0:nx-1,0:ny-1) (/ df/dx,df/dy,d2f/dx2,d2ydy2,d2f/dxdy /)
            
            type(GaussianBlur)  ::      gg
            
            nx = size(f_in,dim=1)
            ny = size(f_in,dim=2)
            
            fbar = 0.0d0
            npix = 0
            do iy = 0,ny-1
                do ix = 0,nx-1
                    if (f_in(ix,iy)==EXTERIOR_IGNORE) cycle
                    fbar = fbar + f_in(ix,iy)
                    npix = npix + 1
                end do
            end do  
            if (npix>0) fbar = fbar/npix
                
            if (Exterior_dbg) print *,"Lib_Exteriors::setExterior2D info - average ",fbar," in ",npix," pixels"                 
            if (Exterior_dbg) print *,"Lib_Exteriors::setExterior2D info - input image ",nx,ny," border ",border," output image ",size(f_out,dim=1),size(f_out,dim=2)         
            if ((model == EXTERIOR_ZERO).or.(model == EXTERIOR_BLUR0)) fbar = 0
            if (present(f_level)) fbar = f_level
            if (model == EXTERIOR_HILO) fbar = huge(1.0)
            if (model == EXTERIOR_LOHI) fbar = -huge(1.0)
                    
        !---    set default exterior            
            f_out = fbar
            
        !---    set central region.  
            f_out(0:nx-1,0:ny-1) = f_in(0:nx-1,0:ny-1)

            
            
            
                        
            select case(model)
                case (EXTERIOR_HILO)
                    do ix = 0,nx-1
                        if (f_in(ix,0)<0) f_out(ix,-1) = -huge(1.0)                                    
                        if (f_in(ix,ny-1)<0) f_out(ix,ny) = -huge(1.0)            
                    end do
                    do iy = 0,ny-1
                        if (f_in(0,iy)<0) f_out(-1,iy) = -huge(1.0)            
                        if (f_in(nx-1,iy)<0) f_out(nx,iy) = -huge(1.0)            
                    end do
                    if (f_in(0,0)<0)       f_out(-1,-1) = -huge(1.0)
                    if (f_in(nx-1,0)<0)    f_out(nx,-1) = -huge(1.0)
                    if (f_in(0,ny-1)<0)    f_out(-1,ny) = -huge(1.0)
                    if (f_in(nx-1,ny-1)<0) f_out(nx,ny) = -huge(1.0)
                case (EXTERIOR_LOHI)
                    do ix = 0,nx-1
                        if (f_in(ix,0)<0) f_out(ix,-1) = +huge(1.0)                                    
                        if (f_in(ix,ny-1)<0) f_out(ix,ny) = +huge(1.0)            
                    end do
                    do iy = 0,ny-1
                        if (f_in(0,iy)<0) f_out(-1,iy) = +huge(1.0)            
                        if (f_in(nx-1,iy)<0) f_out(nx,iy) = +huge(1.0)            
                    end do
                    if (f_in(0,0)<0)       f_out(-1,-1) = +huge(1.0)
                    if (f_in(nx-1,0)<0)    f_out(nx,-1) = +huge(1.0)
                    if (f_in(0,ny-1)<0)    f_out(-1,ny) = +huge(1.0)
                    if (f_in(nx-1,ny-1)<0) f_out(nx,ny) = +huge(1.0)
                case (EXTERIOR_PBC)
                    if (border>=ny) print *,"Lib_Exteriors::setExterior2D warning - border >= ny"
                    if (border>=nx) print *,"Lib_Exteriors::setExterior2D warning - border >= nx"
                    do iy = -border,ny-1+border
                        jy = mod( iy + 8*ny,ny )
                        do ix = -border,nx-1+border 
                            if ( ( ix*(nx-1-ix)>=0 ).and.( iy*(ny-1-iy)>=0 ) ) cycle
                            jx = mod( ix + 8*nx,nx )
                            f_out(ix,iy) = f_in(jx,jy)
                        end do
                    end do                                      
                case (EXTERIOR_EXTRAPOLATE)
                    if (border>=ny) print *,"Lib_Exteriors::setExterior2D warning - border >= ny"
                    if (border>=nx) print *,"Lib_Exteriors::setExterior2D warning - border >= nx"
                    allocate(df(5,0:nx-1,0:ny-1))
                    df = 0
                    
                !---    find first and second derviatives at each point
                                
                    select case(nx)
                        case (1)
                            !   there can be no gradient in the x-direction
                        case (2)
                            !   there is a gradient in the x-drection, but it is fixed.
                            do iy = 0,ny-1
                                df(1,0,iy) = f_in(1,iy)-f_in(0,iy)
                            end do
                            df(1,1,:) = df(1,0,:)
                        case (3)
                            !   there is a gradient in the x-drection, but the second deriv is fixed.
                            do iy = 0,ny-1
                                df(1,0,iy) = (-3*f_in(0,iy)+4*f_in(1,iy)-f_in(2,iy))/2
                                df(1,2,iy) = (3*f_in(0,iy)-4*f_in(1,iy)+f_in(2,iy))/2
                                df(3,0,iy) = (f_in(0,iy)-2*f_in(1,iy)+f_in(2,iy))
                                df(3,2,iy) = df(3,0,iy)
                            end do                          
                        case default        
                            !   there is a varying gradient in the x-drection
                            do iy = 0,ny-1
                                df(1,0,iy) = (-3*f_in(0,iy)+4*f_in(1,iy)-f_in(2,iy))/2
                                df(3,0,iy) = (2*f_in(0,iy)-5*f_in(1,iy)+4*f_in(2,iy)-f_in(3,iy))
                                df(1,nx-1,iy) = (3*f_in(nx-3,iy)-4*f_in(nx-2,iy)+f_in(nx-1,iy))/2
                                df(3,nx-1,iy) = (-f_in(nx-4,iy)+4*f_in(nx-3,iy)-5*f_in(nx-2,iy)+2*f_in(nx-1,iy))
                            end do
                    end select      
                    
                    select case(ny)
                        case (1)
                            !   there can be no gradient in the y-direction
                        case (2)
                            !   there is a gradient in the y-drection, but it is fixed.
                            do ix = 0,nx-1
                                df(2,ix,0) = f_in(ix,1)-f_in(ix,1)
                            end do
                            df(1,:,1) = df(1,:,0)
                            select case(nx)
                                case (1)
                                    !   no gradient in x-direction
                                case (2)
                                    !   fixed gradient in x-direction
                                    df(5,0:1,0) = df(2,1,0) - df(2,0,0)
                                case default                                                                        
                                    df(5,0,0) = (-3*df(2,0,0) + 4*df(2,1,0) - df(2,2,0))/2
                                    do ix = 1,nx-2
                                        df(5,ix,0) = (df(2,ix+1,0) - df(2,ix-1,0))/2
                                    end do
                                    df(5,nx-1,0) = (3*df(2,nx-3,0) - 4*df(2,nx-2,0) + df(2,nx-1,0))/2
                            end select
                            df(5,0:nx-1,1) = df(5,0:nx-1,0)
                        case (3)
                            !   there is a gradient in the y-drection, but the second deriv is fixed.
                            do ix = 0,nx-1
                                df(2,ix,0) = (-3*f_in(ix,0)+4*f_in(ix,1)-f_in(ix,2))/2
                                df(2,ix,2) = (3*f_in(ix,0)-4*f_in(ix,12)+f_in(ix,2))/2
                                df(4,ix,0) = (f_in(ix,0)-2*f_in(ix,1)+f_in(ix,2))
                                df(4,ix,2) = df(3,ix,0)
                            end do                          
                            select case(nx)
                                case (1)
                                    !   no gradient in x-direction
                                case (2)
                                    !   fixed gradient in x-direction
                                    df(5,0:1,0) = df(2,1,0) - df(2,0,0)
                                    df(5,0:1,ny-1) = df(2,1,ny-1) - df(2,0,ny-1)
                                case default                                                                        
                                    df(5,0,0) = (-3*df(2,0,0) + 4*df(2,1,0) - df(2,2,0))/2
                                    df(5,0,ny-1) = (-3*df(2,0,ny-1) + 4*df(2,1,ny-1) - df(2,2,ny-1))/2
                                    do ix = 1,nx-2
                                        df(5,ix,0) = (df(2,ix+1,0) - df(2,ix-1,0))/2
                                        df(5,ix,ny-1) = (df(2,ix+1,ny-1) - df(2,ix-1,ny-1))/2
                                    end do
                                    df(5,nx-1,0) = (3*df(2,nx-3,0) - 4*df(2,nx-2,0) + df(2,nx-1,0))/2
                                    df(5,nx-1,ny-1) = (3*df(2,nx-3,ny-1) - 4*df(2,nx-2,ny-1) + df(2,nx-1,ny-1))/2
                            end select
                        case default        
                            !   there is a varying gradient in the y-drection
                            do ix = 0,nx-1
                                df(2,ix,0) = (-3*f_in(ix,0)+4*f_in(ix,1)-f_in(ix,2))/2
                                df(4,ix,0) = (2*f_in(ix,0)-5*f_in(ix,1)+4*f_in(ix,2)-f_in(ix,3))
                                df(2,ix,ny-1) = (3*f_in(ix,ny-3)-4*f_in(ix,ny-2)+f_in(ix,ny-1))/2
                                df(4,ix,ny-1) = (-f_in(ix,ny-4)+4*f_in(ix,ny-3)-5*f_in(ix,ny-2)+2*f_in(ix,ny-1))
                            end do
                            select case(nx)
                                case (1)
                                    !   no gradient in x-direction
                                case (2)
                                    !   fixed gradient in x-direction
                                    df(5,0:1,0) = df(2,1,0) - df(2,0,0)
                                    df(5,0:1,ny-1) = df(2,1,ny-1) - df(2,0,ny-1)
                                case default                                                                        
                                    df(5,0,0) = (-3*df(2,0,0) + 4*df(2,1,0) - df(2,2,0))/2
                                    df(5,0,ny-1) = (-3*df(2,0,ny-1) + 4*df(2,1,ny-1) - df(2,2,ny-1))/2
                                    do ix = 1,nx-2
                                        df(5,ix,0) = (df(2,ix+1,0) - df(2,ix-1,0))/2
                                        df(5,ix,ny-1) = (df(2,ix+1,ny-1) - df(2,ix-1,ny-1))/2
                                    end do
                                    df(5,nx-1,0) = (3*df(2,nx-3,0) - 4*df(2,nx-2,0) + df(2,nx-1,0))/2
                                    df(5,nx-1,ny-1) = (3*df(2,nx-3,ny-1) - 4*df(2,nx-2,ny-1) + df(2,nx-1,ny-1))/2
                            end select
                    end select      
                    
                    !df(3:,:,:) = 0.0d0
                    
                !---    fill in borders
                    
                    do iy = 0,ny-1
                        do ix = -border,-1
                            f_out(ix,iy) = f_in(0,iy) + ix*(df(1,0,iy) + ix*df(3,0,iy)/2)
                        end do 
                        do ix = 1,border
                            f_out(nx-1+ix,iy) = f_in(nx-1,iy) + ix*(df(1,nx-1,iy) + ix*df(3,nx-1,iy)/2)
                        end do 
                    end do
                    
                    do ix = 0,nx-1
                        do iy = -border,-1
                            f_out(ix,iy) = f_in(ix,0) + iy*(df(2,ix,0) + iy*df(4,ix,0)/2)
                        end do 
                        do iy = 1,border
                            f_out(ix,ny-1+iy) = f_in(ix,ny-1) + iy*(df(2,ix,ny-1) + iy*df(4,ix,ny-1)/2)
                        end do 
                    end do
                    
                !---    fill in corners
                    
                    do iy = -border,-1
                        do ix = -border,-1
                            f_out(ix,iy) = f_in(0,0) + ix*(df(1,0,0) + ix*df(3,0,0)/2 + iy*df(5,0,0)) &
                                                     + iy*(df(2,0,0) + iy*(df(4,0,0)/2))
                        end do 
                        do ix = 1,border
                            f_out(nx-1+ix,iy) = f_in(nx-1,0) + ix*(df(1,nx-1,0) + ix*df(3,nx-1,0)/2 + iy*df(5,nx-1,0)) &
                                                     + iy*(df(2,nx-1,0) + iy*(df(4,nx-1,0)/2))
                        end do 
                    end do
                    
                    do iy = 1,border
                        do ix = -border,-1
                            f_out(ix,ny-1+iy) = f_in(0,ny-1) + ix*(df(1,0,ny-1) + ix*df(3,0,ny-1)/2 + iy*df(5,0,ny-1)) &
                                                     + iy*(df(2,0,ny-1) + iy*(df(4,0,ny-1)/2))
                        end do 
                        do ix = 1,border
                            f_out(nx-1+ix,ny-1+iy) = f_in(nx-1,ny-1) + ix*(df(1,nx-1,ny-1) + ix*df(3,nx-1,ny-1)/2 + iy*df(5,nx-1,ny-1)) &
                                                     + iy*(df(2,nx-1,ny-1) + iy*(df(4,nx-1,ny-1)/2))
                        end do 
                    end do
                    
                    
                    
                case (EXTERIOR_BLUR)
                
                    allocate(df(-border:nx-1+border,-border:ny-1+border,1))
                    df = GAUSSIANBLUR_IGNORE
                   ! df(-border,:,1) = fbar
                   ! df(nx-1+border,:,1) = fbar
                   ! df(:,-border,1) = fbar
                   ! df(:,ny-1+border,1) = fbar
                     
                    df(0:nx-1,0:ny-1,1) = f_in(0:nx-1,0:ny-1)
                    do kk = 1,border
                        
                        gg = GaussianBlur_ctor( kk * 5.0d0/12 )     !   compromise between 1/3, which is minimum for 3 sig boundary, and 1/2, which is 2 sig
                        iy = -kk
                        do ix = -kk,nx-1+kk
                            f_out(ix,iy) = blur( gg,df(:,:,1),ix+border,iy+border )
                        end do                                      
                        ix = -kk
                        do iy = 1-kk,ny-2+kk
                            f_out(ix,iy) = blur( gg,df(:,:,1),ix+border,iy+border )
                        end do                                      
                        ix = nx-1+kk
                        do iy = 1-kk,ny-2+kk
                            f_out(ix,iy) = blur( gg,df(:,:,1),ix+border,iy+border )
                        end do                                      
                        iy = ny-1+kk
                        do ix = -kk,nx-1+kk
                            f_out(ix,iy) = blur( gg,df(:,:,1),ix+border,iy+border )
                        end do                                      
                        call delete(gg)
                        
                    end do
                            
                    
                    
                case (EXTERIOR_BLURPBC)
                
                    
                
                
               !    allocate(df2(-border:nx-1+border,-border:ny-1+border))
               !    df2 = GAUSSIANBLUR_IGNORE
               !    df2(0:nx-1,0:ny-1) = f_in(0:nx-1,0:ny-1)
               !    
               !    allocate(df3(-border:nx-1+border,-border:ny-1+border))
               !    gg = GaussianBlur_ctor( border * 0.5d0 )
               !    call GaussianBlurImage( gg,df2,df3 )
               !    where (df2 == GAUSSIANBLUR_IGNORE)
               !        f_out = df3
               !    end where
                    
                    
                     allocate(df(-border:nx-1+border,-border:ny-1+border,1))
                     df = GAUSSIANBLUR_IGNORE
                    ! df(-border,:,1) = fbar
                    ! df(nx-1+border,:,1) = fbar
                    ! df(:,-border,1) = fbar
                    ! df(:,ny-1+border,1) = fbar
                      
                     df(0:nx-1,0:ny-1,1) = f_in(0:nx-1,0:ny-1)
                     do kk = 1,border
                         
                         gg = GaussianBlur_ctor( kk * 0.5d0 )     !   compromise between 1/3, which is minimum for 3 sig boundary, and 1/2, which is 2 sig
                         iy = -kk
                         do ix = -kk,nx-1+kk
                             f_out(ix,iy) = blur( gg,df(:,:,1),ix+border,iy+border,pbc=.true. )
                         end do                                      
                         ix = -kk
                         do iy = 1-kk,ny-2+kk
                             f_out(ix,iy) = blur( gg,df(:,:,1),ix+border,iy+border,pbc=.true. )
                         end do                                      
                         ix = nx-1+kk
                         do iy = 1-kk,ny-2+kk
                             f_out(ix,iy) = blur( gg,df(:,:,1),ix+border,iy+border,pbc=.true.)
                         end do                                      
                         iy = ny-1+kk
                         do ix = -kk,nx-1+kk
                             f_out(ix,iy) = blur( gg,df(:,:,1),ix+border,iy+border,pbc=.true. )
                         end do                                      
                         call delete(gg)
                         
                     end do
                   !        
                case (EXTERIOR_BLUR0)
                
                    allocate(df(-border:nx-1+border,-border:ny-1+border,1))
                    df = 0.0d0
                    df(1-border:nx-2+border,1-border:ny-2+border,1) = GAUSSIANBLUR_IGNORE
                    df(0:nx-1,0:ny-1,1) = f_in(0:nx-1,0:ny-1)
                    
                    do kk = 1,border-1
                        
                        gg = GaussianBlur_ctor( kk * 5.0d0/12 )     !   compromise between 1/3, which is minimum for 3 sig boundary, and 1/2, which is 2 sig
                        iy = -kk
                        do ix = -kk,nx-1+kk
                            f_out(ix,iy) = blur( gg,df(:,:,1),ix+border,iy+border )
                        end do                                      
                        ix = -kk
                        do iy = 1-kk,ny-2+kk
                            f_out(ix,iy) = blur( gg,df(:,:,1),ix+border,iy+border )
                        end do                                      
                        ix = nx-1+kk
                        do iy = 1-kk,ny-2+kk
                            f_out(ix,iy) = blur( gg,df(:,:,1),ix+border,iy+border )
                        end do                                      
                        iy = ny-1+kk
                        do ix = -kk,nx-1+kk
                            f_out(ix,iy) = blur( gg,df(:,:,1),ix+border,iy+border )
                        end do                                      
                        call delete(gg)
                        
                    end do
                    
                case (EXTERIOR_BLURAVG)
                
                    allocate(df(-border:nx-1+border,-border:ny-1+border,1))
                    !df = fbar
                    !df(1-border:nx-2+border,1-border:ny-2+border,1) = GAUSSIANBLUR_IGNORE
                    
                    df = GAUSSIANBLUR_IGNORE
                    df(-border,0:ny-1,1) = fbar
                    df(nx-1+border,0:ny-1,1) = fbar
                    df(0:nx-1,-border,1) = fbar
                    df(0:nx-1,ny-1+border,1) = fbar
                    
                    df(0:nx-1,0:ny-1,1) = f_in(0:nx-1,0:ny-1)
                    
                    do kk = 1,border-1
                        
                        gg = GaussianBlur_ctor( kk * 5.0d0/12 )     !   compromise between 1/3, which is minimum for 3 sig boundary, and 1/2, which is 2 sig
                        iy = -kk
                        do ix = -kk,nx-1+kk
                            f_out(ix,iy) = blur( gg,df(:,:,1),ix+border,iy+border )
                        end do                                      
                        ix = -kk
                        do iy = 1-kk,ny-2+kk
                            f_out(ix,iy) = blur( gg,df(:,:,1),ix+border,iy+border )
                        end do                                      
                        ix = nx-1+kk
                        do iy = 1-kk,ny-2+kk
                            f_out(ix,iy) = blur( gg,df(:,:,1),ix+border,iy+border )
                        end do                                      
                        iy = ny-1+kk
                        do ix = -kk,nx-1+kk
                            f_out(ix,iy) = blur( gg,df(:,:,1),ix+border,iy+border )
                        end do                                      
                        call delete(gg)
                        
                    end do
                            
            end select
            
            return
        end subroutine setExterior2D

        
    end module Lib_Exteriors
        