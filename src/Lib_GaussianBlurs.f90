
    module Lib_GaussianBlurs
!---^^^^^^^^^^^^^^^^^^^^^^^^
!*      perform a pixel-by-pixel blur on a 2D dataset
!*      ignore pixels labelled by GAUSSIANBLUR_IGNORE
        use iso_fortran_env
        implicit none
        private
        
    !---
    
        public      ::      GaussianBlur_ctor
        public      ::      delete
        public      ::      report
        
        public      ::      blur
        public      ::      gaussianBlurImage
        
    !---
    
        logical,public          ::      GaussianBlur_dbg = .false.
        integer(kind=int64),private,parameter               ::      BADF00D = int(z'BADF00D',kind=int64)
        real(kind=real64),public,parameter                  ::      GAUSSIANBLUR_IGNORE = transfer( (BADF00D+ishft(BADF00D,32_int64)),1.0d0 )
                   
    !---
    
        type,public     ::      GaussianBlur
            private
            real(kind=real64)               ::      sigma           !   width of blur in pixel units
            integer                         ::      m               !   size of blur region 
            real(kind=real64),dimension(:,:),pointer    ::      g   !   [-m:m,-m:m]
        end type GaussianBlur
        
    !---
    
        interface GaussianBlur_ctor
            module procedure    GaussianBlur_null
            module procedure    GaussianBlur_ctor0
        end interface
                
        interface delete
            module procedure    delete0
        end interface
        
        interface report
            module procedure    report0
        end interface
        
        interface blur
            module procedure    blur0
            module procedure    blur1
        end interface
        
    contains
!---^^^^^^^^

        function GaussianBlur_null() result(this)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
            type(GaussianBlur)           ::      this
            this%sigma = 0.0d0
            this%m = 0 
            nullify(this%g)            
            return
        end function GaussianBlur_null
                         
        function GaussianBlur_ctor0(sigma) result(this)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
            type(GaussianBlur)           ::      this
            real(kind=real64),intent(in) ::      sigma
            integer             ::      ix,iy
            real(kind=real64)   ::      i2s2,dd
            this%sigma = abs(sigma)
            this%m = ceiling(3*sigma)
            allocate(this%g(-this%m:this%m,-this%m:this%m))
            if (this%m == 0) then
                this%g = 1.0d0      !   no blur
            else
                i2s2 = 1/(2*this%sigma*this%sigma)
                do iy = -this%m,this%m
                    do ix = -this%m,this%m
                        dd = real(ix*ix + iy*iy,kind=real64)
                        dd = dd*i2s2
                        if (dd>9.0d0) then
                            this%g(ix,iy) = 0.0d0
                        else
                            this%g(ix,iy) = exp(-dd)
                        end if
                    end do
                end do
            end if
            return
        end function GaussianBlur_ctor0
                       
    !---
    
        subroutine delete0(this)
    !---^^^^^^^^^^^^^^^^^^^^^^^^
            type(GaussianBlur),intent(inout)    ::      this
            if (this%m == 0) return
            deallocate(this%g)
            this = GaussianBlur_null()
            return
        end subroutine delete0
        
    !---
    
        subroutine report0(this,u,o)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^
            type(GaussianBlur),intent(in)    ::      this
            integer,intent(in),optional     ::      u,o
            integer     ::      uu,oo
            uu = 6 ; if (present(u)) uu = u
            oo = 0 ; if (present(o)) oo = o
            write(unit=uu,fmt='(a,f16.4,a)') repeat(" ",oo)//"GaussianBlur [sigma=",this%sigma,"]"            
            return
        end subroutine report0
    
    !---
    
        function blur0( this,img,ix,iy ) result( gimg) 
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      return gaussian blur of 2d img centred on ix,iy
            type(GaussianBlur),intent(in)    ::      this
            real(kind=real64),dimension(0:,0:),intent(in)   ::  img
            integer,intent(in)              ::      ix,iy
            real(kind=real64)               ::  gimg
            integer                 ::      nx,ny,kx,ky 
            real(kind=real64)       ::      ff,gg,gsum
            nx = size(img,dim=1)
            ny = size(img,dim=2)
            gimg = 0.0d0 ; gsum = 0.0d0 
            do ky = max(0,iy-this%m),min(ny-1,iy+this%m)
                do kx = max(0,ix-this%m),min(nx-1,ix+this%m)
                    ff = img(kx,ky)
                    if (ff == GAUSSIANBLUR_IGNORE) cycle
                    gg = this%g(kx-ix,ky-iy)
                    gsum = gsum + gg
                    gimg = gimg + ff*gg
                end do
            end do
            if (gsum > 0) then
                gimg = gimg / gsum
            else
                gimg = GAUSSIANBLUR_IGNORE
            end if
            return
        end function blur0
                    
        function blur1( this,img,ix,iy,pbc ) result( gimg) 
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      return gaussian blur of 2d img centred on ix,iy
    !*      periodic boundary conditions
            type(GaussianBlur),intent(in)    ::      this
            real(kind=real64),dimension(0:,0:),intent(in)   ::  img
            integer,intent(in)              ::      ix,iy
            logical,intent(in)              ::      pbc
            real(kind=real64)               ::  gimg
            integer                 ::      nx,ny,kx,ky,jx,jy
            real(kind=real64)       ::      ff,gg,gsum
            
            if (.not. pbc) then
                gimg = blur0( this,img,ix,iy )
                return
             
            end if
            
            nx = size(img,dim=1)
            ny = size(img,dim=2)
            gimg = 0.0d0 ; gsum = 0.0d0 
            do jy = -this%m,this%m
                ky = mod ( iy + jy + ny, ny )
                do jx = -this%m,this%m
                    kx = mod ( ix + jx + nx, nx )
                    ff = img(kx,ky)
                    if (ff == GAUSSIANBLUR_IGNORE) cycle
                    gg = this%g(jx,jy)
                    gsum = gsum + gg
                    gimg = gimg + ff*gg
                end do
            end do
            if (gsum > 0) then
                gimg = gimg / gsum
            else
                gimg = GAUSSIANBLUR_IGNORE
            end if
            return
        end function blur1
        
        
         subroutine gaussianBlurImage( this, f_in, f_out )
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ 
    !*      produce a Gaussian blur of input image, width t
            type(GaussianBlur),intent(in)    ::      this
            real(kind=real64),dimension(0:,0:),intent(in)       ::      f_in
            real(kind=real64),dimension(0:,0:),intent(inout)    ::      f_out
            
            
            
            
            real(kind=real64),dimension(:),allocatable      ::      kernel ,f_stripe
            real(kind=real64),dimension(:,:),allocatable      ::      w_stripe
              
            
            integer             ::      Nx,Ny
            integer             ::      ix,iy,jj,kk
            integer             ::      Nk
            real(kind=real64)   ::      i2s2,ww,wf,ws,ff
            
            Nx = size(f_in,dim=1)
            Ny = size(f_in,dim=2)
             
        !---    compute the unnormalised kernel
            Nk = max(3,ceiling( this%sigma*3 ))             !   range of pixels to search is +/- Nk
            allocate(kernel( 0:Nk ))
            i2s2 = 1/(2*this%sigma*this%sigma)             
            kernel(0) = 1.0d0   
            do ix = 1,Nk                
                kernel( ix )  = exp( -ix*ix*i2s2 )                                
            end do
            
                 
        !---    compute the output blurred image. First do x strips.
            allocate(f_stripe(0:max(Nx,Ny)-1))
            allocate(w_stripe(0:Nx-1,0:Ny-1))
          
          ! do ix = 0,Nx-1
          !     wf = 0.0d0 ; ws = 0.0d0
          !     do jj = max(0,ix-Nk),min(Nx-1,ix+Nk)
          !         ff = f_in(jj,0)
          !         if ( ff == GAUSSIANBLUR_IGNORE ) cycle
          !         kk = abs(jj-ix)           
          !         ww = kernel( kk )           
          !                                      
          !         wf = wf + ww*ff
          !         ws = ws + ww
          !     end do
          !     f_out(ix,0) = wf
          !     w_stripe(ix,0) = ws
          ! end do
                         
            !f_out = GAUSSIANBLUR_IGNORE
            !w_stripe = 0.0d0
            do iy = 0,Ny-1
                f_stripe(0:Nx-1) = f_in(0:Nx-1,iy)
                do ix = 0,Nx-1
                    wf = 0.0d0 ; ws = 0.0d0
                    do jj = max(0,ix-Nk),min(Nx-1,ix+Nk)
                        ff = f_stripe(jj)
                        if ( ff == GAUSSIANBLUR_IGNORE ) cycle
                        kk = abs(jj-ix)           
                        ww = kernel( kk )                                         
                        wf = wf + ww*ff
                        ws = ws + ww
                    end do
                            
                    if (ws>0) then
                        f_out(ix,iy) = wf   
                        w_stripe(ix,iy) = ws      
                    else
                        f_out(ix,iy) = GAUSSIANBLUR_IGNORE
                        w_stripe(ix,iy) = 1.0d0     
                    end if      
                end do
            end do
 
 
        !   at this point, f_out has blurring in the x-direction only. weight stores the kernel weighting from this op.            
            
                             
        !---    now do y strips
            
            do ix = 0,Nx-1
            !   make a copy of this stripe
                f_stripe(0:Ny-1) = f_out(ix,0:Ny-1)
                do iy = 0,Ny-1
                    wf = 0.0d0 ; ws = 0.0d0
                    do jj = max(0,iy-Nk),min(Ny-1,iy+Nk)
                        ff = f_stripe(jj)
                        if ( ff == GAUSSIANBLUR_IGNORE ) cycle
                        kk = abs(jj-iy)                                  
                        ww = w_stripe(ix,iy)*kernel( kk )
                        wf = wf + ww*ff
                        ws = ws + ww
                    end do
                    !ws = w_stripe(ix,iy)*ws
                    if (ws>0) then
                        f_out(ix,iy) = wf/ws
                    else
                        f_out(ix,iy) = GAUSSIANBLUR_IGNORE
                    end if
                     
                                     
                                        
                end do
            end do
    
        !---    now have a normalised Gaussian blur function
                      
            
            return
        end subroutine gaussianBlurImage
                          
         
        
    end module Lib_GaussianBlurs
    
    
!!   gfortran -ffree-line-length-256 -Og -g Lib_GaussianBlurs.f90 -o testGaussianBlurs.exe
!
!    
!    program testGaussianBlurs
!!---^^^^^^^^^^^^^^^^^^^^^^^^
!        use iso_fortran_env
!        use Lib_GaussianBlurs
!        implicit none
!        
!        type(GaussianBlur)           ::      this
!        
!        this = GaussianBlur_ctor()
!        
!        call report(this)
!        
!        call delete(this)
!        
!        print *,""
!        print *,"done"
!        print *,""
!        
!    end program testGaussianBlurs
!    
!        
    
        
        
    