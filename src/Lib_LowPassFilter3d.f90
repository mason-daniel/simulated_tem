
    module Lib_LowPassFilter3d
!---^^^^^^^^^^^^^^^^^^^^^^^^^^^
        use iso_fortran_env
        implicit none
        private
        
    !---
    
        public      ::      generateFourierCoefficients
        public      ::      interpolateFromFourierCoefficients        
     
    !---
    
        logical,public          ::      LowPassFilter3d_dbg = .false.
        real(kind=real64),private,parameter                     ::      PI = 3.14159265358979d0
                   
    !---
    
        interface         interpolateFromFourierCoefficients
            module procedure    interpolateFromFourierCoefficients21
            module procedure    interpolateFromFourierCoefficients22
            module procedure    interpolateFromFourierCoefficients31
            module procedure    interpolateFromFourierCoefficients32
        end interface
    
        interface         generateFourierCoefficients
            module procedure    generateFourierCoefficients2
            module procedure    generateFourierCoefficients3
        end interface
            
        
        
    contains
!---^^^^^^^^


        subroutine generateFourierCoefficients2(x,dat,xmax,nFourierCoeffs,qdat)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      given the data dat(1:d,1:n) at 2d points x(1:2:n) ranged from 0:xmax(1:2)
    !*      compute the first few fourier coefficients
            real(kind=real64),dimension(:,:),intent(in)     ::      x
            real(kind=real64),dimension(:,:),intent(in)     ::      dat
            real(kind=real64),dimension(2),intent(in)       ::      xmax
            integer,dimension(2),intent(in)                 ::      nFourierCoeffs
            real(kind=real64),dimension(:,:,0:,0:),intent(out)      ::      qdat          !   (1:d,1:4,0:nFourierCoeffs,0:nFourierCoeffs)
            
            
            integer                     ::          ix,iy
            real(kind=real64)           ::          qx,qy
            
            integer                     ::          ii,nn,dd
            real(kind=real64)           ::          cx,sx,cy,sy , xx,yy
            real(kind=real64)           ::          d2,d1
            
    
            dd = size(dat,dim=1)
            nn = size(dat,dim=2)
            
            qdat = 0.0d0
            d1 = 2*PI/xmax(1)          
            d2 = 2*PI/xmax(2)   
            
            do ii = 1,nn
                xx = x(1,ii)
                yy = x(2,ii)
           
                do iy = 0,nFourierCoeffs(2)
                    qy = yy*d2*iy
                    cy = cos( qy )
                    sy = sin( qy )
                    
                    do ix = 0,nFourierCoeffs(1)
                        qx = xx*d1*ix
                        cx = cos( qx )
                        sx = sin( qx )
                                                      
                        qdat( 1:dd,1,ix,iy ) = qdat( 1:dd,1,ix,iy ) + cx*cy*dat(1:dd,ii)    
                        qdat( 1:dd,2,ix,iy ) = qdat( 1:dd,2,ix,iy ) + sx*cy*dat(1:dd,ii) 
                        qdat( 1:dd,3,ix,iy ) = qdat( 1:dd,3,ix,iy ) + cx*sy*dat(1:dd,ii) 
                        qdat( 1:dd,4,ix,iy ) = qdat( 1:dd,4,ix,iy ) + sx*sy*dat(1:dd,ii) 
                                      
                    end do
                end do     
            end do        
                                    
            xx = 1.0d0/nn
            
            do iy = 0,nFourierCoeffs(2)
                d2 = 2 ; if (iy==0) d2 = 1
                do ix = 0,nFourierCoeffs(1)
                    d1 = 2 ; if (ix==0) d1 = 1
                    qdat(:,:,ix,iy) = qdat(:,:,ix,iy) * d1*d2*xx
                end do
            end do
            
             
                        
            return
        end subroutine generateFourierCoefficients2                            
                        
        subroutine generateFourierCoefficients3(x,dat,xmax,nFourierCoeffs,qdat)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      given the data dat(1:d,1:n) at 3d points x(1:3,1:n) ranged from 0:xmax(1:3)
    !*      compute the first few fourier coefficients
            real(kind=real64),dimension(:,:),intent(in)     ::      x
            real(kind=real64),dimension(:,:),intent(in)     ::      dat
            real(kind=real64),dimension(3),intent(in)       ::      xmax
            integer,dimension(3),intent(in)                 ::      nFourierCoeffs
            real(kind=real64),dimension(:,:,0:,0:,0:),intent(out)      ::      qdat          !   (1:d,1:8,0:nFourierCoeffs,0:nFourierCoeffs,0:nFourierCoeffs)
            
            
            integer                     ::          ix,iy,iz
            real(kind=real64)           ::          qx,qy,qz
            
            integer                     ::          ii,nn,dd
            real(kind=real64)           ::          cx,sx,cy,sy,cz,sz , xx,yy,zz
            real(kind=real64)           ::          d3,d2,d1
            real(kind=real64),dimension(1:size(dat,dim=1))  ::      dcz,dsz
    
            dd = size(dat,dim=1)
            nn = size(dat,dim=2)
            
            qdat = 0.0d0
            d1 = 2*PI/xmax(1)          
            d2 = 2*PI/xmax(2)   
            d3 = 2*PI/xmax(3)   
            
            do ii = 1,nn
                xx = x(1,ii)
                yy = x(2,ii)
                zz = x(3,ii)
               
            
                do iz = 0,nFourierCoeffs(3)
                    qz = zz*d3*iz
                    cz = cos( qz )
                    sz = sin( qz )
                    dcz(1:dd) = dat(1:dd,ii)*cz
                    dsz(1:dd) = dat(1:dd,ii)*sz
                    
                    do iy = 0,nFourierCoeffs(2)
                        qy = yy*d2*iy
                        cy = cos( qy )
                        sy = sin( qy )
                        
                        do ix = 0,nFourierCoeffs(1)
                            qx = xx*d1*ix
                            cx = cos( qx )
                            sx = sin( qx )
                                                          
                            qdat( 1:dd,1,ix,iy,iz ) = qdat( 1:dd,1,ix,iy,iz ) + cx*cy*dcz(1:dd)     
                            qdat( 1:dd,2,ix,iy,iz ) = qdat( 1:dd,2,ix,iy,iz ) + sx*cy*dcz(1:dd) 
                            qdat( 1:dd,3,ix,iy,iz ) = qdat( 1:dd,3,ix,iy,iz ) + cx*sy*dcz(1:dd) 
                            qdat( 1:dd,4,ix,iy,iz ) = qdat( 1:dd,4,ix,iy,iz ) + sx*sy*dcz(1:dd) 
                            qdat( 1:dd,5,ix,iy,iz ) = qdat( 1:dd,5,ix,iy,iz ) + cx*cy*dsz(1:dd) 
                            qdat( 1:dd,6,ix,iy,iz ) = qdat( 1:dd,6,ix,iy,iz ) + sx*cy*dsz(1:dd) 
                            qdat( 1:dd,7,ix,iy,iz ) = qdat( 1:dd,7,ix,iy,iz ) + cx*sy*dsz(1:dd) 
                            qdat( 1:dd,8,ix,iy,iz ) = qdat( 1:dd,8,ix,iy,iz ) + sx*sy*dsz(1:dd) 
                                          
                        end do
                    end do     
                end do
            end do        
                                    
            xx = 1.0d0/nn
            
            do iz = 0,nFourierCoeffs(3)
                d3 = 2 ; if (iz==0) d3 = 1
                do iy = 0,nFourierCoeffs(2)
                    d2 = 2 ; if (iy==0) d2 = 1
                    do ix = 0,nFourierCoeffs(1)
                        d1 = 2 ; if (ix==0) d1 = 1
                        qdat(:,:,ix,iy,iz) = qdat(:,:,ix,iy,iz) * d1*d2*d3*xx
                    end do
                end do
            end do
            
             
                        
            return
        end subroutine generateFourierCoefficients3                            
        
    !---        
        
                        
        function interpolateFromFourierCoefficients21(x,xmax,nFourierCoeffs,qdat) result(dat)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
            real(kind=real64),dimension(2),intent(in)       ::      x
            real(kind=real64),dimension(2),intent(in)       ::      xmax
            integer,dimension(2),intent(in)                 ::      nFourierCoeffs
            real(kind=real64),dimension(:,:,0:,0:),intent(in)      ::      qdat                 !   (1:d,1:4,0:nFourierCoeffs,0:nFourierCoeffs)
            real(kind=real64),dimension(size(qdat,dim=1))   ::      dat
            
            
            integer                     ::          ix,iy 
            real(kind=real64)           ::          qx,qy 
            
            integer                     ::          dd
            real(kind=real64)           ::          cx,sx,cy,sy 
            real(kind=real64)           ::          d2,d1
    
            dd = size(dat)
                        
            d1 = x(1)*2*PI/xmax(1)          
            d2 = x(2)*2*PI/xmax(2)   
            
            dat = 0
        
            do iy = 0,nFourierCoeffs(2)
                qy = d2*iy
                sy = sin( qy )
                cy = cos( qy )
                                                       
                do ix = 0,nFourierCoeffs(1)
                    qx = d1*ix
                    sx = sin( qx )
                    cx = cos( qx )                        
                    
                   dat(1:dd) = dat(1:dd) + ( ( qdat( 1:dd,1,ix,iy ) * cx             & 
                                         +     qdat( 1:dd,2,ix,iy ) * sx )*cy        & 
                                         +   ( qdat( 1:dd,3,ix,iy ) * cx             & 
                                         +     qdat( 1:dd,4,ix,iy ) * sx )*sy )       
                end do
            end do
                            
            return
        end function interpolateFromFourierCoefficients21        
                            
                        
        function interpolateFromFourierCoefficients22(nx,ny,nFourierCoeffs,qdat) result(dat)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      interpolate at each point of a regularly spaced lattice
            integer,intent(in)                              ::      nx,ny
            integer,dimension(2),intent(in)                 ::      nFourierCoeffs
            real(kind=real64),dimension(:,:,0:,0:),intent(in)      ::      qdat                 !   (1:d,1:4,0:nFourierCoeffs,0:nFourierCoeffs)
            real(kind=real64),dimension(size(qdat,dim=1),0:nx-1,0:ny-1)   ::      dat
            
            
            integer                     ::          ix,iy 
            integer                     ::          jx,jy 
            real(kind=real64)           ::          qx,qy 
            
            integer                     ::          dd
            real(kind=real64)           ::          cx,sx,cy,sy 
            real(kind=real64)           ::          d2,d1  
    
            real(kind=real64),dimension(size(qdat,dim=1))       ::      ddat
            
            dd = size(qdat,dim=1)
                       
            d1 = 2*PI/nx          
            d2 = 2*PI/ny   
    
            
            dat = 0
            
        
            do jy = 0,ny-1
            do iy = 0,nFourierCoeffs(2)
                qy = d2*iy*jy
                sy = sin( qy )
                cy = cos( qy )
                             
                do jx = 0,nx-1                          
                do ix = 0,nFourierCoeffs(1)
                    qx = d1*ix*jx
                    sx = sin( qx )
                    cx = cos( qx )                        
                    
                    ddat(1:dd) = ( ( qdat( 1:dd,1,ix,iy ) * cx             & 
                               +     qdat( 1:dd,2,ix,iy ) * sx )*cy        & 
                               +   ( qdat( 1:dd,3,ix,iy ) * cx             & 
                               +     qdat( 1:dd,4,ix,iy ) * sx )*sy ) 
                                
                    dat(1:dd,jx,jy) = dat(1:dd,jx,jy) + ddat(1:dd)                                   
                               
                end do
                end do
            end do
            end do
                            
            return
        end function interpolateFromFourierCoefficients22        
                            
                          
            
              
        function interpolateFromFourierCoefficients31(x,xmax,nFourierCoeffs,qdat) result(dat)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
            real(kind=real64),dimension(3),intent(in)       ::      x
            real(kind=real64),dimension(3),intent(in)       ::      xmax
            integer,dimension(3),intent(in)                 ::      nFourierCoeffs
            real(kind=real64),dimension(:,:,0:,0:,0:),intent(in)      ::      qdat                 !   (1:d,1:8,0:nFourierCoeffs,0:nFourierCoeffs,0:nFourierCoeffs)
            real(kind=real64),dimension(size(qdat,dim=1))   ::      dat
            
            
            integer                     ::          ix,iy,iz
            real(kind=real64)           ::          qx,qy,qz
            
            integer                     ::          dd
            real(kind=real64)           ::          cx,sx,cy,sy,cz,sz
            real(kind=real64)           ::          d3,d2,d1
    
            dd = size(dat)
                        
            d1 = x(1)*2*PI/xmax(1)          
            d2 = x(2)*2*PI/xmax(2)   
            d3 = x(3)*2*PI/xmax(3)   
            
            dat = 0
            
            do iz = 0,nFourierCoeffs(3)
                qz = d3*iz
                sz = sin( qz )
                cz = cos( qz )                
                
                do iy = 0,nFourierCoeffs(2)
                    qy = d2*iy
                    sy = sin( qy )
                    cy = cos( qy )
                                                           
                    do ix = 0,nFourierCoeffs(1)
                        qx = d1*ix
                        sx = sin( qx )
                        cx = cos( qx )                        
                        
                       dat(1:dd) = dat(1:dd) + ( ( qdat( 1:dd,1,ix,iy,iz ) * cx             & 
                                             +     qdat( 1:dd,2,ix,iy,iz ) * sx )*cy        & 
                                             +   ( qdat( 1:dd,3,ix,iy,iz ) * cx             & 
                                             +     qdat( 1:dd,4,ix,iy,iz ) * sx )*sy )*cz   & 
                                             + ( ( qdat( 1:dd,5,ix,iy,iz ) * cx             & 
                                             +     qdat( 1:dd,6,ix,iy,iz ) * sx )*cy        & 
                                             +   ( qdat( 1:dd,7,ix,iy,iz ) * cx             & 
                                             +     qdat( 1:dd,8,ix,iy,iz ) * sx )*sy )*sz     
                        
                    end do
                end do
            end do      
                            
            return
        end function interpolateFromFourierCoefficients31        
                            
                        
        function interpolateFromFourierCoefficients32(nx,ny,nz,nFourierCoeffs,qdat) result(dat)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      interpolate at each point of a regularly spaced lattice
            integer,intent(in)                              ::      nx,ny,nz
            integer,dimension(3),intent(in)                 ::      nFourierCoeffs
            real(kind=real64),dimension(:,:,0:,0:,0:),intent(in)      ::      qdat                 !   (1:d,1:8,0:nFourierCoeffs,0:nFourierCoeffs,0:nFourierCoeffs)
            real(kind=real64),dimension(size(qdat,dim=1),0:nx-1,0:ny-1,0:nz-1)   ::      dat
            
            
            integer                     ::          ix,iy,iz
            integer                     ::          jx,jy,jz
            real(kind=real64)           ::          qx,qy,qz
            
            integer                     ::          dd
            real(kind=real64)           ::          cx,sx,cy,sy,cz,sz
            real(kind=real64)           ::          d3,d2,d1  
    
            real(kind=real64),dimension(size(qdat,dim=1))       ::      ddat
            
            dd = size(qdat,dim=1)
                        
            
            
            
            d1 = 2*PI/nx          
            d2 = 2*PI/ny   
            d3 = 2*PI/nz   
            
            dat = 0
            
            do jz = 0,nz-1                
            do iz = 0,nFourierCoeffs(3)
                qz = d3*iz*jz
                sz = sin( qz )
                cz = cos( qz )                
                
                do jy = 0,ny-1
                do iy = 0,nFourierCoeffs(2)
                    qy = d2*iy*jy
                    sy = sin( qy )
                    cy = cos( qy )
                                 
                    do jx = 0,nx-1                          
                    do ix = 0,nFourierCoeffs(1)
                        qx = d1*ix*jx
                        sx = sin( qx )
                        cx = cos( qx )                        
                        
                        ddat(1:dd) = ( ( qdat( 1:dd,1,ix,iy,iz ) * cx             & 
                                   +     qdat( 1:dd,2,ix,iy,iz ) * sx )*cy        & 
                                   +   ( qdat( 1:dd,3,ix,iy,iz ) * cx             & 
                                   +     qdat( 1:dd,4,ix,iy,iz ) * sx )*sy )*cz   & 
                                   + ( ( qdat( 1:dd,5,ix,iy,iz ) * cx             & 
                                   +     qdat( 1:dd,6,ix,iy,iz ) * sx )*cy        & 
                                   +   ( qdat( 1:dd,7,ix,iy,iz ) * cx             & 
                                   +     qdat( 1:dd,8,ix,iy,iz ) * sx )*sy )*sz     
                        
                        dat(1:dd,jx,jy,jz) = dat(1:dd,jx,jy,jz) + ddat(1:dd)                                   
                                   
                    end do
                    end do
                end do
                end do
            end do      
            end do      
                            
            return
        end function interpolateFromFourierCoefficients32        
                            
                          
            
            

        
    end module Lib_LowPassFilter3d
    
!       gfortran -ffree-line-length-256  Lib_LowPassFilter3d.f90   
!    
!    program testLib_LowPassFilter3d
!!---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
!        use Lib_LowPassFilter3d
!        use iso_fortran_env
!        implicit none
!        
!        integer             ::      Nx = 6
!        real(kind=real64)   ::      L = 1.0d0
!        real(kind=real64)   ::      PI = 3.14159265358979d0
!                                                         
!        integer,dimension(3)     ::      nFourierCoeffs = 1
!        integer     ::      d = 1
!        
!        
!        integer     ::      ix,iy,iz,nn
!        real(kind=real64),dimension(3)   ::      xmax 
!        real(kind=real64),dimension(:,:),allocatable            ::      x,dat
!        real(kind=real64),dimension(:,:,:,:,:),allocatable      ::      qdat
!        real(kind=real64),dimension(1)  ::      yy
!        xmax(1:3) = L
!        
!        allocate(dat(d,Nx*Nx*Nx))
!        allocate(x(3,Nx*Nx*Nx))
!        
!        nn = 0
!        do iz = 0,Nx-1
!            do iy = 0,Nx-1
!                do ix = 0,Nx-1
!                    nn = nn + 1
!                    x(1:3,nn) = (/ ix,iy,iz /)*L/Nx
!                    
!                    dat(1:d,nn) = 1.0d0 + 0.3d0*cos( 2*PI*x(1,nn)/L ) + 0.2d0*sin( 2*PI*x(2,nn)/L ) + 0.1d0*sin( 2*PI*x(3,nn)/L )*cos( 2*PI*x(1,nn)/L )
!                    
!               end do
!           end do
!        end do
!        
!        allocate(qdat(d,8,0:nFourierCoeffs(1),0:nFourierCoeffs(2),0:nFourierCoeffs(3)))
!        call generateFourierCoefficients(x,dat,xmax,nFourierCoeffs,qdat)
!        
!        do iz = 0,nFourierCoeffs(3)
!            do iy = 0,nFourierCoeffs(2)
!                do ix = 0,nFourierCoeffs(1)
!                    write(*,fmt='(a,3i4,8f16.8)') "FC ",ix,iy,iz,qdat(1,:,ix,iy,iz)
!                end do
!            end do
!        end do
!        
!        do ix = 1,nn
!            print *,x(1:3,ix),dat(1,ix),interpolateFromFourierCoefficients(x(:,ix),xmax,nFourierCoeffs,qdat) 
!        end do
!        
!        
!        
!        
!        
!        print *,""
!        print *,"done"
!        print *,""
!        
!        
!    end program testLib_LowPassFilter3d
!        
!        