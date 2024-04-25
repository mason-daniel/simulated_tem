
    module Lib_FFTW3f
!---^^^^^^^^^^^^^^^^^
!*      wrapper module to call fftw3
!*      constructs the plans etc


        use iso_c_binding
        use iso_fortran_env
        use Lib_SimpleProgressBar
        implicit none
        include "fftw3.f03"

        private

        real(kind=real64),private,parameter     ::      PI = 3.14159265358979d0
        logical,public                          ::      FFTW_DBG = .false.

        public      ::      FFT1d
        public      ::      FFT2d
        public      ::      FFT3d

        public      ::      radialPowerFunction
        public      ::      radialPowerSpectrum

        interface   FFT1d
            module procedure        dft_r2c_1d
            module procedure        dft_r2I_1d
        end interface
        
        interface   FFT2d
            module procedure        dft_r2c_2d
            module procedure        dft_r2I_2d
        end interface

        interface   FFT3d
            module procedure        dft_r2c_3d
            module procedure        dft_r2I_3d
            module procedure        dft_r2I_3d_inplace_old
        end interface

        interface   radialPowerFunction
            module procedure        radialPowerFunction2d
            module procedure        radialPowerFunction3d
        end interface

        interface   radialPowerSpectrum
            module procedure        radialPowerSpectrum2d
        end interface


    contains
!---^^^^^^^^




        subroutine dft_r2c_1d( in,out )
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      performs a 1d FFT
    !*      returns a 1d complex valued solution - note out(i,j) is complex conjugate of out(Nx-i,j)
            real(kind=real64),dimension(:),intent(inout)      ::        in
            complex(kind=real64),dimension(:),intent(out)     ::        out

            type(c_ptr)                             ::      plan

            integer(kind=c_int)                     ::      Nx,Mx



            Nx = size(in,dim=1)
           
            Mx = int(Nx/2+1)        !   note real array gives complex FT with N/2 output values!

        !---    forward FT
            plan = fftw_plan_dft_r2c_1d( Nx, in, out, FFTW_ESTIMATE )                                    !   note: C has row-major order
            call fftw_execute_dft_r2c(plan,in, out)

        !---    tidy up
            call fftw_destroy_plan(plan)


            return
        end subroutine dft_r2c_1d

        subroutine dft_r2I_1d( in,out )
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      performs the 1d FFT on in
    !*      returns a 1d real valued power spectrum with dimension( 2 (Nx/2+1) - 1 , Ny )
            real(kind=real64),dimension(:),intent(inout)      ::        in
           ! real(kind=real64),dimension(-(size(in,dim=1)/2):,-(size(in,dim=2)/2):),intent(inout)        ::        out
            real(kind=real64),dimension(0:),intent(inout)        ::        out

            integer(kind=c_int)                     ::      Nx,Mx
            complex(kind=real64),dimension(:),allocatable     ::        in_ft
            integer                     ::      ix                                            
            real(kind=real64)           ::      aa,bb,rr
            complex(kind=real64)        ::      ft

            Nx = size(in,dim=1)
   

            Mx = int(Nx/2)                      !   note real array gives complex FT with N/2 output values!                             
            allocate(in_ft(0:Mx))
            call dft_r2c_1d( in,in_ft )

            if (FFTW_DBG) print *,"Lib_FFTW3f::dft_r2I_1d info - mean ",sum(in)/(Nx),in_ft(0)

            do ix = 0,Mx
            
                ft = in_ft( ix ) 

                aa = real(ft)
                bb = aimag(ft)
                rr = aa*aa + bb*bb


                out( ix ) = rr
            end do

             
            aa = 1.0d0/(Nx)
            out = out * aa

        !---    tidy up
            deallocate(in_ft)
            return
        end subroutine dft_r2I_1d





        subroutine dft_r2c_2d( in,out )
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      performs a 2d FFT
    !*      returns a 2d complex valued solution - note out(i,j) is complex conjugate of out(Nx-i,j)
            real(kind=real64),dimension(:,:),intent(inout)      ::        in
            complex(kind=real64),dimension(:,:),intent(out)     ::        out

            type(c_ptr)                             ::      plan

            integer(kind=c_int)                     ::      Nx,Ny,Mx



            Nx = size(in,dim=1)
            Ny = size(in,dim=2)

            Mx = int(Nx/2+1)        !   note real array gives complex FT with N/2 output values!

        !---    forward FT
            plan = fftw_plan_dft_r2c_2d( Ny,Nx, in, out, FFTW_ESTIMATE )                                    !   note: C has row-major order
            call fftw_execute_dft_r2c(plan,in, out)

        !---    tidy up
            call fftw_destroy_plan(plan)


            return
        end subroutine dft_r2c_2d

        subroutine dft_r2I_2d( in,out )
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      performs the 2d FFT on in
    !*      returns a 2d real valued power spectrum with dimension( 2 (Nx/2+1) - 1 , Ny )
            real(kind=real64),dimension(:,:),intent(inout)      ::        in
           ! real(kind=real64),dimension(-(size(in,dim=1)/2):,-(size(in,dim=2)/2):),intent(inout)        ::        out
            real(kind=real64),dimension(0:,0:),intent(inout)        ::        out

            integer(kind=c_int)                     ::      Nx,Ny,Mx
            complex(kind=real64),dimension(:,:),allocatable     ::        in_ft
            integer                     ::      ix,iy
            real(kind=real64)           ::      aa,bb,rr
            complex(kind=real64)        ::      ft

            Nx = size(in,dim=1)
            Ny = size(in,dim=2)

            Mx = int(Nx/2)                      !   note real array gives complex FT with N/2 output values!
            allocate(in_ft(0:Mx,0:Ny-1))
            call dft_r2c_2d( in,in_ft )

            if (FFTW_DBG) print *,"Lib_FFTW3f::dft_r2I_2d info - mean ",sum(in)/(Nx*Ny),in_ft(0,0)

            do iy = -Ny/2,Ny/2
                do ix = -Mx,Mx

                    if (ix<0) then
                        if (iy<=0) then
                            ft = conjg( in_ft( -ix,-iy ) )
                        else
                            ft = conjg( in_ft( -ix,Ny-iy ) )
                        end if
                    else
                        if (iy<0) then
                            ft = ( in_ft( ix,Ny+iy ) )
                        else if (iy == 0) then
                            ft = ( in_ft( ix,0 ) )
                        else
                            ft = ( in_ft( ix,iy ) )
                        end if
                    end if


                    aa = real(ft)
                    bb = aimag(ft)
                    rr = aa*aa + bb*bb


                    !out( ix,iy ) = rr
                    out( ix+Mx,iy+Ny/2 ) = rr
                end do
            end do

            aa = 1.0d0/(Nx*Ny)
            out = out * aa

        !---    tidy up
            deallocate(in_ft)
            return
        end subroutine dft_r2I_2d




        subroutine dft_r2c_3d( in,out )
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      performs a 3d FFT
    !*      returns a 3d complex valued solution - note out(i,j) is complex conjugate of out(Nx-i,j)
            real(kind=real64),dimension(:,:,:)                    ::        in
            complex(kind=real64),dimension(:,:,:),intent(out)     ::        out

            type(c_ptr)                             ::      plan

            integer(kind=c_int)                     ::      Nx,Ny,Nz    !,Mx



            Nx = size(in,dim=1)
            Ny = size(in,dim=2)
            Nz = size(in,dim=3)


        !---    forward FT
            plan = fftw_plan_dft_r2c_3d( Nz,Ny,Nx, in, out, FFTW_ESTIMATE )                                    !   note: C has row-major order
            call fftw_execute_dft_r2c(plan,in, out)

        !---    tidy up
            call fftw_destroy_plan(plan)


            return
        end subroutine dft_r2c_3d




        subroutine dft_r2I_3d( in,out )
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      performs the 2d FFT on in
    !*      returns a 2d real valued power spectrum with dimension( 2 (Nx/2+1) - 1 , Ny )
            real(kind=real64),dimension(:,:,:),intent(inout)      ::        in

            real(kind=real64),dimension(0:,0:,0:),intent(inout)        ::        out

            integer(kind=c_int)                     ::      Nx,Ny,Nz,Mx
            complex(kind=real64),dimension(:,:,:),allocatable     ::        in_ft
            integer                     ::      ix,iy,iz
            real(kind=real64)           ::      aa,bb,rr , dx,dy,dz

            complex(kind=real64)        ::      fft,fft100,fftJ00,fft010,fft0J0,fft001,fft00J


            Nx = size(in,dim=1)
            Ny = size(in,dim=2)
            Nz = size(in,dim=3)

            Mx = int(Nx/2 )                      !   note real array gives complex FT with N/2 output values!


            if (FFTW_DBG) print *,"Lib_FFTW3f::dft_r2I_3d info - minmaxavg ",minval(in(:,:,:)),",",maxval(in(:,:,:)),",",sum(in(:,:,:))/(Nx*Ny*Nz)


            allocate(in_ft(0:Mx,0:Ny-1,0:Nz-1))
            call dft_r2c_3d( in,in_ft )

            if (FFTW_DBG) print *,"Lib_FFTW3f::dft_r2I_3d info - mean ",sum(in)/(Nx*Ny*Nz),in_ft(0,0,0)

            do iz = -Nz/2,Nz/2
                do iy = -Ny/2,Ny/2
                    do ix = -Mx,Mx

                        fft = ft( ix,iy,iz )

                        if ( (abs(iz)==Nz/2).or.(abs(iy)==Ny/2).or.(abs(ix)==Nx/2) ) then


                            aa = real(fft)
                            bb = aimag(fft)
                            rr = aa*aa + bb*bb

                        else

                            fft100 = ft( ix+1,iy,iz )
                            fftJ00 = ft( ix-1,iy,iz )
                            fft010 = ft( ix,iy+1,iz )
                            fft0J0 = ft( ix,iy-1,iz )
                            fft001 = ft( ix,iy,iz+1 )
                            fft00J = ft( ix,iy,iz-1 )

                         ! aa = real( fft + fft100 + fftJ00 + fft010 + fft0J0 + fft001 + fft00J )/7
                         ! dx = ( real( fft100 ) - real( fftJ00 ) )/2
                         ! dy = ( real( fft010 ) - real( fft0J0 ) )/2
                         ! dz = ( real( fft001 ) - real( fft00J ) )/2
                         !
                         ! rr = aa*aa + (dx*dx+dy*dy+dz*dz)/12
                         !
                         ! bb = aimag( fft + fft100 + fftJ00 + fft010 + fft0J0 + fft001 + fft00J )/7
                         !  dx = ( aimag( fft100 ) - aimag( fftJ00 ) )/2
                         !  dy = ( aimag( fft010 ) - aimag( fft0J0 ) )/2
                         !  dz = ( aimag( fft001 ) - aimag( fft00J ) )/2
                         !
                         !  rr = rr + bb*bb + (dx*dx+dy*dy+dz*dz)/12

                            fft = ( fft + fft100 + fftJ00 + fft010 + fft0J0 + fft001 + fft00J )/7
                            fft100 = (fft100 - fftJ00)/2
                            fft010 = (fft010 - fft0J0)/2
                            fft001 = (fft001 - fft00J)/2

                            aa = real( fft )
                            dx = real( fft100 )
                            dy = real( fft010 )
                            dz = real( fft001 )

                            rr = aa*aa + (dx*dx+dy*dy+dz*dz)/12

                            bb = aimag( fft )
                            dx = aimag( fft100 )
                            dy = aimag( fft010 )
                            dz = aimag( fft001 )

                            rr = rr + bb*bb + (dx*dx+dy*dy+dz*dz)/12

                        end if


                        !out( ix,iy,iz ) = rr
                        out( ix+Mx,iy+Ny/2,iz+Nz/2 ) = rr
                    end do
                end do
            end do

            aa = (1.0d0/Nx) * (1.0d0/Ny) * (1.0d0/Nz)       !   because Nx Ny Nz can be > 3 G
            out = out * aa

        !---    tidy up
            deallocate(in_ft)
            return


        contains
    !---^^^^^^^^

            pure complex(kind=real64) function ft( ix,iy,iz )
        !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
                integer,intent(in)          ::      ix,iy,iz

                if (ix<0) then
                    if (iy<0) then
                        if (iz<0) then
                            ft = conjg( in_ft( -ix,-iy,-iz ) )
                        !else if (iz==0) then
                        !    ft = conjg( in_ft( -ix,-iy,0 ) )
                        else
                            ft = conjg( in_ft( -ix,-iy,Nz-iz ) )
                        end if
                    else if (iy == 0) then
                        !ft = conjg( in_ft( -ix,0 ) )
                        if (iz<0) then
                            ft = conjg( in_ft( -ix,0,-iz ) )
                        !else if (iz==0) then
                        !    ft = conjg( in_ft( -ix,0,0 ) )
                        else
                            ft = conjg( in_ft( -ix,0,Nz-iz ) )
                        end if
                    else
                        !ft = conjg( in_ft( -ix,Ny-iy ) )
                        if (iz<0) then
                            ft = conjg( in_ft( -ix,Ny-iy,-iz ) )
                        !else if (iz==0) then
                        !    ft = conjg( in_ft( -ix,Ny-iy,0 ) )
                        else
                            ft = conjg( in_ft( -ix,Ny-iy,Nz-iz ) )
                        end if
                    end if
                else
                    if (iy<0) then
                        !ft = ( in_ft( ix,Ny+iy ) )
                        if (iz<0) then
                            ft = ( in_ft( ix,Ny+iy,-iz ) )
                        !else if (iz==0) then
                        !    ft = ( in_ft( ix,Ny+iy,0 ) )
                        else
                            ft = ( in_ft( ix,Ny+iy,Nz-iz ) )
                        end if
                    else if (iy == 0) then
                        !ft = ( in_ft( ix,0 ) )
                        if (iz<0) then
                            ft = ( in_ft( ix,0,-iz ) )
                        !else if (iz==0) then
                        !    ft = ( in_ft( ix,0,0 ) )
                        else
                            ft = ( in_ft( ix,0,Nz-iz ) )
                        end if
                    else
                        !ft = ( in_ft( ix,iy ) )
                        if (iz<0) then
                            ft = ( in_ft( ix,iy,-iz ) )
                        !else if (iz==0) then
                        !    ft = ( in_ft( ix,iy,0 ) )
                        else
                            ft = ( in_ft( ix,iy,Nz-iz ) )
                        end if
                    end if
                end if
                return
            end function ft




        end subroutine dft_r2I_3d





        subroutine dft_r2I_3d_inplace_27( inout )
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      performs the 2d FFT on in
    !*      returns a 2d real valued power spectrum with dimension( 2 (Nx/2+1) - 1 , Ny )
            real(kind=real64),dimension(:,:,:),allocatable,intent(inout)            ::        inout


            integer(kind=c_int)                     ::      Nx,Ny,Nz,Mx,My,Mz
            complex(kind=real64),dimension(:,:,:),allocatable     ::        in_ft
            integer                     ::      ix,iy,iz
            real(kind=real64)           ::      aa,bb,rr , invox , dx,dy,dz,dxx,dxy,dyy,dyz,dzz,dzx

            complex(kind=real64),dimension(-1:1,-1:1,-1:1)      ::  fft
            complex(kind=real64)        ::      u0,ux,uy,uz,uxx,uxy,uyy,uyz,uzz,uzx
            integer                     ::      jy,jz

            real(kind=real64),dimension(-1:1,-1:1,-1:1),parameter     ::      u0_kernel = reshape( (/         &
                            -2, 1,-2    ,    1, 4, 1    ,   -2, 1,-2,                               &
                             1, 4, 1    ,    4, 7, 4    ,    1, 4, 1,                               &
                            -2, 1,-2    ,    1, 4, 1    ,   -2, 1,-2                     /),(/3,3,3/) )/27.0d0

            real(kind=real64),dimension(-1:1,-1:1,-1:1),parameter     ::      ux_kernel = reshape( (/         &
                             1, 0,-1    ,     1, 0,-1   ,    1, 0,-1,                               &
                             1, 0,-1    ,     1, 0,-1   ,    1, 0,-1,                               &
                             1, 0,-1    ,     1, 0,-1   ,    1, 0,-1                     /),(/3,3,3/) )/18.0d0

            real(kind=real64),dimension(-1:1,-1:1,-1:1),parameter     ::      uy_kernel = reshape( (/         &
                            -1,-1,-1    ,     0, 0, 0   ,    1, 1, 1,                               &
                            -1,-1,-1    ,     0, 0, 0   ,    1, 1, 1,                               &
                            -1,-1,-1    ,     0, 0, 0   ,    1, 1, 1                     /),(/3,3,3/) )/18.0d0

            real(kind=real64),dimension(-1:1,-1:1,-1:1),parameter     ::      uz_kernel = reshape( (/         &
                            -1,-1,-1    ,    -1,-1,-1   ,   -1,-1,-1,                               &
                             0, 0, 0    ,     0, 0, 0   ,    0, 0, 0,                               &
                             1, 1, 1    ,     1, 1, 1   ,    1, 1, 1                     /),(/3,3,3/) )/18.0d0


            real(kind=real64),dimension(-1:1,-1:1,-1:1),parameter     ::      uxx_kernel = reshape( (/            &
                             1,-2, 1    ,     1,-2, 1   ,    1,-2, 1,                               &
                             1,-2, 1    ,     1,-2, 1   ,    1,-2, 1,                               &
                             1,-2, 1    ,     1,-2, 1   ,    1,-2, 1                     /),(/3,3,3/) )/9.0d0

            real(kind=real64),dimension(-1:1,-1:1,-1:1),parameter     ::      uyy_kernel = reshape( (/            &
                             1, 1, 1    ,    -2,-2,-2   ,    1, 1, 1,                               &
                             1, 1, 1    ,    -2,-2,-2   ,    1, 1, 1,                               &
                             1, 1, 1    ,    -2,-2,-2   ,    1, 1, 1                     /),(/3,3,3/) )/9.0d0

            real(kind=real64),dimension(-1:1,-1:1,-1:1),parameter     ::      uzz_kernel = reshape( (/            &
                             1, 1, 1    ,     1, 1, 1   ,    1, 1, 1,                               &
                            -2,-2,-2    ,    -2,-2,-2   ,   -2,-2,-2,                               &
                             1, 1, 1    ,     1, 1, 1   ,    1, 1, 1                     /),(/3,3,3/) )/9.0d0


            real(kind=real64),dimension(-1:1,-1:1,-1:1),parameter     ::      uxy_kernel = reshape( (/            &
                            -1, 0, 1    ,     0, 0, 0   ,    1, 0,-1,                               &
                            -1, 0, 1    ,     0, 0, 0   ,    1, 0,-1,                               &
                            -1, 0, 1    ,     0, 0, 0   ,    1, 0,-1                     /),(/3,3,3/) )/12.0d0

            real(kind=real64),dimension(-1:1,-1:1,-1:1),parameter     ::      uyz_kernel = reshape( (/            &
                             1, 1, 1    ,     0, 0, 0   ,   -1,-1,-1,                               &
                             0, 0, 0    ,     0, 0, 0   ,    0, 0, 0,                               &
                            -1,-1,-1    ,     0, 0, 0   ,    1, 1, 1                     /),(/3,3,3/) )/12.0d0

            real(kind=real64),dimension(-1:1,-1:1,-1:1),parameter     ::      uzx_kernel = reshape( (/            &
                            -1, 0, 1    ,    -1, 0, 1   ,   -1, 0, 1,                               &
                             0, 0, 0    ,     0, 0, 0   ,    0, 0, 0,                               &
                             1, 0,-1    ,     1, 0,-1   ,    1, 0,-1                     /),(/3,3,3/) )/12.0d0



          !  complex(kind=real64)        ::      ft

            Nx = size(inout,dim=1)
            Ny = size(inout,dim=2)
            Nz = size(inout,dim=3)

            Mx = int(Nx/2)        !   note real array gives complex FT with N/2 output values!
            My = int(Ny/2)        !
            Mz = int(Nz/2)        !

            invox = 1.0d0/( real(Nx,kind=real64)*real(Ny,kind=real64)*real(Nz,kind=real64) )


        !---    compute the power spectrum
            if (FFTW_DBG) then
                print *,"Lib_FFTW3f::radialPowerFunction3d_inplace info - input data size Nx,Ny,Nz   = ",Nx,Ny,Nz
                print *,"Lib_FFTW3f::radialPowerFunction3d_inplace info - FT data size    Mx,My,Mz   = ",Mx,My,Mz
            end if

            if (FFTW_DBG) print *,"Lib_FFTW3f::dft_r2I_3d_inplace info - minmaxavg ",minval(inout(:,:,:)),",",maxval(inout(:,:,:)),",",sum(inout(:,:,:))*invox


            allocate(in_ft(0:Mx,0:Ny-1,0:Nz-1))
            call dft_r2c_3d( inout,in_ft )
            deallocate(inout)



            allocate(inout(-Mx:Mx,-My:My,-Mz:Mz))


            if (FFTW_DBG) print *,"Lib_FFTW3f::dft_r2I_3d_inplace info - mean ",in_ft(0,0,0)*invox

        !---    edge cases
            iz = -Mz
            do iy = -My,My
                do ix = -Mx,Mx
                    u0 = ft( ix,iy,iz )
                    aa = real(u0) ; bb = aimag(u0)
                    inout( ix,iy,iz ) = aa*aa + bb*bb
                end do
            end do
            iz = Mz
            do iy = -My,My
                do ix = -Mx,Mx
                    u0 = ft( ix,iy,iz )
                    aa = real(u0) ; bb = aimag(u0)
                    inout( ix,iy,iz ) = aa*aa + bb*bb
                end do
            end do
            iy = -My
            do iz = -Mz,Mz
                do ix = -Mx,Mx
                    u0 = ft( ix,iy,iz )
                    aa = real(u0) ; bb = aimag(u0)
                    inout( ix,iy,iz ) = aa*aa + bb*bb
                end do
            end do
            iy = My
            do iz = -Mz,Mz
                do ix = -Mx,Mx
                    u0 = ft( ix,iy,iz )
                    aa = real(u0) ; bb = aimag(u0)
                    inout( ix,iy,iz ) = aa*aa + bb*bb
                end do
            end do
            ix = -Mx
            do iz = -Mz,Mz
                do iy = -My,My
                    u0 = ft( ix,iy,iz )
                    aa = real(u0) ; bb = aimag(u0)
                    inout( ix,iy,iz ) = aa*aa + bb*bb
                end do
            end do
            ix = Mx
            do iz = -Mz,Mz
                do iy = -My,My
                    u0 = ft( ix,iy,iz )
                    aa = real(u0) ; bb = aimag(u0)
                    inout( ix,iy,iz ) = aa*aa + bb*bb
                end do
            end do



            do iz = 1-Mz,Mz-1
                do iy = 1-My,My-1

                !---    construct 3x3x3 sliding box of values from FT
                    do jz = -1,1
                        do jy = -1,1
                            fft(0,jy,jz) = ft( -Mx,iy+jy,iz+jz)
                            fft(1,jy,jz) = ft(1-Mx,iy+jy,iz+jz)
                        end do
                    end do
                    
                    call progressBar( (iy+My)+(iz+Mz-1)*(2*My-1) , (2*My-1)+(2*Mz-2)*(2*My-1) )

                    do ix = 1-Mx,Mx-1


                    !---    shift entries in 3x3x3 box one step to left, and fill in new values
                        fft(-1,:,:) = fft( 0,:,:)
                        fft( 0,:,:) = fft( 1,:,:)
                        do jz = -1,1
                            do jy = -1,1
                                fft(1,jy,jz) = ft(ix+1,iy+jy,iz+jz)
                            end do
                        end do


                        u0  = sum( u0_kernel*fft  )
                        ux  = sum( ux_kernel*fft  )
                        uy  = sum( uy_kernel*fft  )
                        uz  = sum( uz_kernel*fft  )
                        uxx = sum( uxx_kernel*fft )
                        uxy = sum( uxy_kernel*fft )
                        uyy = sum( uyy_kernel*fft )
                        uyz = sum( uyz_kernel*fft )
                        uzz = sum( uzz_kernel*fft )
                        uzx = sum( uzx_kernel*fft )                                                                                  

                        aa  = real(u0 )
                        dx  = real(ux )
                        dy  = real(uy )
                        dz  = real(uz )
                        dxx = real(uxx)
                        dxy = real(uxy)
                        dyy = real(uyy)
                        dyz = real(uyz)
                        dzz = real(uzz)
                        dzx = real(uzx)

                        rr = aa*aa                              &
                            + (dx*dx + dy*dy + dz*dz)/12        &
                            + aa*(dxx+dyy+dzz)/12               &
                            + (dxy*dxy+dyz*dyz+dzx*dzx)/144     &
                            + (dxx*dyy+dyy*dzz+dzz*dxx)/288     &
                            + (dxx*dxx+dyy*dyy+dzz*dzz)/320

                        aa  = aimag(u0 )
                        dx  = aimag(ux )
                        dy  = aimag(uy )
                        dz  = aimag(uz )
                        dxx = aimag(uxx)
                        dxy = aimag(uxy)
                        dyy = aimag(uyy)
                        dyz = aimag(uyz)
                        dzz = aimag(uzz)
                        dzx = aimag(uzx)

                        rr = rr + aa*aa                         &
                            + (dx*dx + dy*dy + dz*dz)/12        &
                            + aa*(dxx+dyy+dzz)/12               &
                            + (dxy*dxy+dyz*dyz+dzx*dzx)/144     &
                            + (dxx*dyy+dyy*dzz+dzz*dxx)/288     &
                            + (dxx*dxx+dyy*dyy+dzz*dzz)/320



                        inout( ix,iy,iz ) = rr



                    end do
                end do
            end do


            inout = inout * invox

        !---    tidy up
            deallocate(in_ft)
            return



        contains
    !---^^^^^^^^

            pure complex(kind=real64) function ft( ix,iy,iz )
        !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
                integer,intent(in)          ::      ix,iy,iz

                if (ix<0) then
                    if (iy<=0) then
                        if (iz<=0) then
                            ft = conjg( in_ft( -ix,-iy,-iz ) )
                        else
                            ft = conjg( in_ft( -ix,-iy,Nz-iz ) )
                        end if
                    else
                        if (iz<=0) then
                            ft = conjg( in_ft( -ix,Ny-iy,-iz ) )
                        else
                            ft = conjg( in_ft( -ix,Ny-iy,Nz-iz ) )
                        end if
                    end if
                else
                    if (iy<0) then
                        if (iz<=0) then
                            ft = ( in_ft( ix,Ny+iy,-iz ) )
                        else
                            ft = ( in_ft( ix,Ny+iy,Nz-iz ) )
                        end if
                    else if (iy == 0) then
                        if (iz<=0) then
                            ft = ( in_ft( ix,0,-iz ) )
                        else
                            ft = ( in_ft( ix,0,Nz-iz ) )
                        end if
                    else
                        if (iz<=0) then
                            ft = ( in_ft( ix,iy,-iz ) )
                        else
                            ft = ( in_ft( ix,iy,Nz-iz ) )
                        end if
                    end if
                end if
                return
            end function ft



        end subroutine dft_r2I_3d_inplace_27





        subroutine dft_r2I_3d_inplace_old( inout )
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      performs the 2d FFT on in
    !*      returns a 2d real valued power spectrum with dimension( 2 (Nx/2+1) - 1 , Ny )
            real(kind=real64),dimension(:,:,:),allocatable,intent(inout)            ::        inout


            integer(kind=c_int)                     ::      Nx,Ny,Nz,Mx,My,Mz
            complex(kind=real64),dimension(:,:,:),allocatable     ::        in_ft
            integer                     ::      ix,iy,iz
            real(kind=real64)           ::      aa,bb,rr , invox
            complex(kind=real64)        ::      fft


            Nx = size(inout,dim=1)
            Ny = size(inout,dim=2)
            Nz = size(inout,dim=3)

            Mx = int(Nx/2)        !   note real array gives complex FT with N/2 output values!
            My = int(Ny/2)        !
            Mz = int(Nz/2)        !

            invox = 1.0d0/( real(Nx,kind=real64)*real(Ny,kind=real64)*real(Nz,kind=real64) )


        !---    compute the power spectrum
            if (FFTW_DBG) then
                print *,"Lib_FFTW3f::radialPowerFunction3d_inplace info - input data size Nx,Ny,Nz   = ",Nx,Ny,Nz
                print *,"Lib_FFTW3f::radialPowerFunction3d_inplace info - FT data size    Mx,My,Mz   = ",Mx,My,Mz
            end if

            if (FFTW_DBG) print *,"Lib_FFTW3f::dft_r2I_3d_inplace info - minmaxavg ",minval(inout(:,:,:)),",",maxval(inout(:,:,:)),",",sum(inout(:,:,:))*invox


            allocate(in_ft(0:Mx,0:Ny-1,0:Nz-1))
            call dft_r2c_3d( inout,in_ft )
            deallocate(inout)



            allocate(inout(-Mx:Mx,-My:My,-Mz:Mz))


            if (FFTW_DBG) print *,"Lib_FFTW3f::dft_r2I_3d_inplace info - mean ",in_ft(0,0,0)*invox

            do iz = -Mz,Mz
                do iy = -My,My
                    do ix = -Mx,Mx
! 
!                         if (ix<0) then
!                             if (iy<=0) then
!                                 if (iz<=0) then
!                                     ft = conjg( in_ft( -ix,-iy,-iz ) )
!                                 else
!                                     ft = conjg( in_ft( -ix,-iy,Nz-iz ) )
!                                 end if
!                             else
!                                 !ft = conjg( in_ft( -ix,Ny-iy ) )
!                                 if (iz<=0) then
!                                     ft = conjg( in_ft( -ix,Ny-iy,-iz ) )
!                                 else
!                                     ft = conjg( in_ft( -ix,Ny-iy,Nz-iz ) )
!                                 end if
!                             end if
!                         else
!                             if (iy<0) then
!                                 !ft = ( in_ft( ix,Ny+iy ) )
!                                 if (iz<=0) then
!                                     ft = ( in_ft( ix,Ny+iy,-iz ) )
!                                 else
!                                     ft = ( in_ft( ix,Ny+iy,Nz-iz ) )
!                                 end if
!                             else if (iy == 0) then
!                                 !ft = ( in_ft( ix,0 ) )
!                                 if (iz<=0) then
!                                     ft = ( in_ft( ix,0,-iz ) )
!                                 else
!                                     ft = ( in_ft( ix,0,Nz-iz ) )
!                                 end if
!                             else
!                                 !ft = ( in_ft( ix,iy ) )
!                                 if (iz<=0) then
!                                     ft = ( in_ft( ix,iy,-iz ) )
!                                 else
!                                     ft = ( in_ft( ix,iy,Nz-iz ) )
!                                 end if
!                             end if
!                         end if
! 
                        fft = ft( ix,iy,iz )
                        aa = real(fft)
                        bb = aimag(fft)
                        rr = aa*aa + bb*bb

                       !




                        inout( ix,iy,iz ) = rr
                    end do
                end do
            end do


            inout = inout * invox

        !---    tidy up
            deallocate(in_ft)
            return

        contains
    !---^^^^^^^^^            


            pure complex(kind=real64) function ft( ix,iy,iz )
        !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
                integer,intent(in)          ::      ix,iy,iz

                if (ix<0) then
                    if (iy<=0) then
                        if (iz<=0) then
                            ft = conjg( in_ft( -ix,-iy,-iz ) )
                        else
                            ft = conjg( in_ft( -ix,-iy,Nz-iz ) )
                        end if
                    else
                        if (iz<=0) then
                            ft = conjg( in_ft( -ix,Ny-iy,-iz ) )
                        else
                            ft = conjg( in_ft( -ix,Ny-iy,Nz-iz ) )
                        end if
                    end if
                else
                    if (iy<0) then
                        if (iz<=0) then
                            ft = ( in_ft( ix,Ny+iy,-iz ) )
                        else
                            ft = ( in_ft( ix,Ny+iy,Nz-iz ) )
                        end if
                    else if (iy == 0) then
                        if (iz<=0) then
                            ft = ( in_ft( ix,0,-iz ) )
                        else
                            ft = ( in_ft( ix,0,Nz-iz ) )
                        end if
                    else
                        if (iz<=0) then
                            ft = ( in_ft( ix,iy,-iz ) )
                        else
                            ft = ( in_ft( ix,iy,Nz-iz ) )
                        end if
                    end if
                end if
                return
            end function ft

        end subroutine dft_r2I_3d_inplace_old





        subroutine dft_r2I_3d_inplace_7( inout )
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      performs the 2d FFT on in
    !*      returns a 2d real valued power spectrum with dimension( 2 (Nx/2+1) - 1 , Ny )
            real(kind=real64),dimension(:,:,:),allocatable,intent(inout)            ::        inout


            integer(kind=c_int)                     ::      Nx,Ny,Nz,Mx,My,Mz
            complex(kind=real64),dimension(:,:,:),allocatable     ::        in_ft
            integer                     ::      ix,iy,iz
            real(kind=real64)           ::      aa,bb,rr , invox ,dx,dy,dz
            complex(kind=real64)        ::      fft0,fft100,fftJ00,fft010,fft0J0,fft001,fft00J


            Nx = size(inout,dim=1)
            Ny = size(inout,dim=2)
            Nz = size(inout,dim=3)

            Mx = int(Nx/2)        !   note real array gives complex FT with N/2 output values!
            My = int(Ny/2)        !
            Mz = int(Nz/2)        !

            invox = 1.0d0/( real(Nx,kind=real64)*real(Ny,kind=real64)*real(Nz,kind=real64) )


        !---    compute the power spectrum
            if (FFTW_DBG) then
                print *,"Lib_FFTW3f::radialPowerFunction3d_inplace_7 info - input data size Nx,Ny,Nz   = ",Nx,Ny,Nz
                print *,"Lib_FFTW3f::radialPowerFunction3d_inplace_7 info - FT data size    Mx,My,Mz   = ",Mx,My,Mz
            end if

            if (FFTW_DBG) print *,"Lib_FFTW3f::dft_r2I_3d_inplace_7 info - minmaxavg ",minval(inout(:,:,:)),",",maxval(inout(:,:,:)),",",sum(inout(:,:,:))*invox


            allocate(in_ft(0:Mx,0:Ny-1,0:Nz-1))
            call dft_r2c_3d( inout,in_ft )
            deallocate(inout)



            allocate(inout(-Mx:Mx,-My:My,-Mz:Mz))


            if (FFTW_DBG) print *,"Lib_FFTW3f::dft_r2I_3d_inplace_7 info - mean ",in_ft(0,0,0)*invox

            do iz = -Mz,Mz
                do iy = -My,My
                
                    call progressBar( (iy+My+1)+(iz+Mz)*(2*My+1) , (2*My+1)+(2*Mz)*(2*My+1) )
                
                    do ix = -Mx,Mx
!           
! 
                        fft0 = ft( ix,iy,iz )
                        
                        if ( (abs(ix)==Mx).or.(abs(iy)==My).or.(abs(iz)==Mz) ) then
                        
                            aa = real(fft0)
                            bb = aimag(fft0)
                            rr = aa*aa + bb*bb

                        else
                       !
                            fft100 = ft( ix+1,iy,iz )
                            fftJ00 = ft( ix-1,iy,iz )
                            fft010 = ft( ix,iy+1,iz )
                            fft0J0 = ft( ix,iy-1,iz )
                            fft001 = ft( ix,iy,iz+1 )
                            fft00J = ft( ix,iy,iz-1 )
                            
                            fft0 = (fft0 + fft100 + fftJ00 + fft010 + fft0J0 + fft001 + fft00J)/7
                            fft100 = (fft100 - fftJ00)/2
                            fft010 = (fft010 - fft0J0)/2
                            fft001 = (fft001 - fft00J)/2
                            
                            aa = real( fft0 )
                            dx = real( fft100 ) 
                            dy = real( fft010 ) 
                            dz = real( fft001 ) 
                            
                            rr = aa*aa + (dx*dx+dy*dy+dz*dz)/12
                            
                            
                            aa = aimag( fft0 )
                            dx = aimag( fft100 ) 
                            dy = aimag( fft010 ) 
                            dz = aimag( fft001 ) 
                            
                            rr = rr + aa*aa + (dx*dx+dy*dy+dz*dz)/12
                        end if

                        inout( ix,iy,iz ) = rr
                    end do
                end do
            end do


            inout = inout * invox

        !---    tidy up
            deallocate(in_ft)
            return

        contains
    !---^^^^^^^^^            


            pure complex(kind=real64) function ft( ix,iy,iz )
        !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
                integer,intent(in)          ::      ix,iy,iz

                if (ix<0) then
                    if (iy<=0) then
                        if (iz<=0) then
                            ft = conjg( in_ft( -ix,-iy,-iz ) )
                        else
                            ft = conjg( in_ft( -ix,-iy,Nz-iz ) )
                        end if
                    else
                        if (iz<=0) then
                            ft = conjg( in_ft( -ix,Ny-iy,-iz ) )
                        else
                            ft = conjg( in_ft( -ix,Ny-iy,Nz-iz ) )
                        end if
                    end if
                else
                    if (iy<0) then
                        if (iz<=0) then
                            ft = ( in_ft( ix,Ny+iy,-iz ) )
                        else
                            ft = ( in_ft( ix,Ny+iy,Nz-iz ) )
                        end if
                    else if (iy == 0) then
                        if (iz<=0) then
                            ft = ( in_ft( ix,0,-iz ) )
                        else
                            ft = ( in_ft( ix,0,Nz-iz ) )
                        end if
                    else
                        if (iz<=0) then
                            ft = ( in_ft( ix,iy,-iz ) )
                        else
                            ft = ( in_ft( ix,iy,Nz-iz ) )
                        end if
                    end if
                end if
                return
            end function ft

        end subroutine dft_r2I_3d_inplace_7









        subroutine radialPowerSpectrum2d( in,q_min,q_max, cprime,rpf )
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      compute the radial power spectrum between q_min and q_max
    !*      also return some simple integrals associated with the power function
    !*      2 pi q C'(q) = int delta( q - |k| ) |F(k)|^2 d2k
            real(kind=real64),dimension(:,:),intent(inout)      ::      in
            real(kind=real64),intent(in)                        ::      q_min,q_max
            real(kind=real64),dimension(0:),intent(out)         ::      cprime          !   int delta( q - |k| ) |F(k)|^2 d2k  / (2 pi q)
            real(kind=real64),dimension(0:),intent(out)         ::      rpf             !   int delta( q - |k| ) |F(k)|^2 d2k
            !real(kind=real64),intent(out)                       ::      cprimesum       !   int cprime dq
            !real(kind=real64),intent(out)                       ::      rpfsum          !   int 2 pi q cprime dq = int  |F(k)|^2 d2k

            real(kind=real64),dimension(:,:),allocatable        ::      out

            integer             ::      Nx,Ny,Mx,My,nBins
            integer             ::      ix,iy,ii,nBinMiss !, imax
            integer,dimension(:),allocatable                    ::     winbin
!            real(kind=real64),dimension(:),allocatable          ::     qrpf,q2rpf, medprod
            real(kind=real64)   ::      qq,deltaqx,deltaqy,ideltaq,inxny  !,noise,dd

!            logical             ::      ok

            Nx = size(in,dim=1)
            Ny = size(in,dim=2)

            Mx = int(Nx/2)        !   note real array gives complex FT with N/2 output values!
            My = int(Ny/2)        !

        !---    compute the power spectrum
            if (FFTW_DBG) then

                print *,"Nx,Ny   = ",Nx,Ny
                print *,"Mx,My   = ",Mx,My
            end if
            allocate(out(-Mx:Mx,-My:My))
            call dft_r2I_2d( in,out )


            if (FFTW_DBG) then
                print *,"Lib_FFTW3f::radialPowerSpectrum2d() dbg - proof of Parseval's theorem"
                print *,"Nx Ny   = ",Nx*Ny
                print *,"<in^2>  = ",sum(in**2)/(Nx*Ny)
                print *,"<out^2> = ",sum(out)  /(Nx*Ny)
            end if


            nBins = size(rpf)-1
            allocate( winbin(0:nBins) )

           ! allocate( qrpf(0:nBins) )                !
           ! allocate( q2rpf(0:nBins) )
            deltaqx = PI/Mx     !    2*PI/Nx
            deltaqy = PI/My     !    2*PI/Ny
            ideltaq = nBins/q_max
            inxny = 1.0d0/(Nx*Ny)
            print *,"Lib_FFTW3f::radialPowerSpectrum2d() info - deltaqx,deltaqy,deltaq ",deltaqx,deltaqy,1/ideltaq

            winbin = 0
            rpf = 0              !   compute as sum delta( |q| - q ) out(q) d2q
            !rpfsum = 0           !   compute as sum out(q) d2q
            !qrpfbar = 0         !   compute as sum |q| out(q) d2q, then normalises
            !q2rpfbar = 0        !   compute as sum |q^2| out(q) d2q, then normalises
            nBinMiss = 0
            !qrpf = 0
            !q2rpf = 0
            do iy = -My,My
                do ix = -Mx,Mx
                    qq = sqrt( (ix*deltaqx)**2 + (iy*deltaqy)**2 )
                    if (qq < q_min) cycle
                    if (qq > q_max) cycle
                    !rpfsum = rpfsum + out(ix,iy)
                    !qrpfbar = qrpfbar + qq*out(ix,iy)
                    !q2rpfbar = q2rpfbar + qq*qq*out(ix,iy)
                    ii = nint( qq*ideltaq )
                    if (ii*(nBins-ii)<0) then
                        nBinMiss = nBinMiss + 1
                        cycle
                    end if
                    winbin(ii) = winbin(ii) + 1

                    rpf(ii) = rpf(ii) + out(ix,iy)
                    !qrpf(ii) = qrpf(ii) + qq*out(ix,iy)
                    !q2rpf(ii) = q2rpf(ii) + qq*qq*out(ix,iy)
                end do
            end do
            !qrpfbar = qrpfbar/max(1.0d-16,rpfsum)
            !q2rpfbar = q2rpfbar/max(1.0d-16,rpfsum)
            !if (FFTW_DBG) print *,"Lib_FFTW3f::radialPowerSpectrum2d() dbg - sum over all q ",rpfsum!,qrpfbar,q2rpfbar
            if (FFTW_DBG) print *,"Lib_FFTW3f::radialPowerSpectrum2d() dbg - number bin misses ",nBinMiss,"/",(2*Mx+1)*(2*My+1)


        !---    construct radial average

            cprime = 0
            do ii = 0,nBins
                if (winbin(ii)>0) then
                    cprime(ii) = rpf(ii)/winbin(ii)
                end if
            end do
            !cprimesum = sum(cprime)




       ! !---    find largest value of rpf
       !     qq = -huge(1.0)
       !     do ii = 0,nBins
       !         !print *,ii,winbin(ii),rpf(ii),qrpf(ii),q2rpf(ii)
       !         if (winbin(ii)>0) qq = max( qq, rpf(ii)/winbin(ii) )
       !     end do
       !     qq = qq/(Nx*Ny)
       !     if (FFTW_DBG) print *,"Lib_FFTW3f::radialPowerFunction() dbg - largest value of rpf() = ",qq
       !     noise = qq*1.0d-7
       !     if (present(snr)) then
       !         if (snr>0) noise = max(noise,qq/(snr*snr) )
       !     end if
       !    !print *,"Lib_FFTW3f::radialPowerFunction() info - using noise level ",noise

       !    do ii = 0,nBins
       !        dd = inxny/max(1,winbin(ii))
       !        rpf(ii) = rpf(ii)*dd
       !   !     qrpf(ii) = qrpf(ii)*dd
       !   !     q2rpf(ii) = q2rpf(ii)*dd
       !    end do
       !


      ! !----    find a good estimate for the level at which pixel noise takes over
      ! !    compute a 5 point median filter for qrpf and q2rpf
      ! !    then find the product/q^3.
      !      allocate( medprod(0:nBins) ) ; medprod = huge(1.0)
      !      do ii = 1,nBins-2
      !          qq = ii*q_max/nBins
      !          dd = 1/(144*qq*qq*qq)
      !          medprod(ii) = ( qrpf(ii-2) +  qrpf(ii+2) + 2*qrpf(ii-1) + 2*qrpf(ii+1) + 6*qrpf(ii) )*( q2rpf(ii-2) + q2rpf(ii+2) + 2*q2rpf(ii-1) + 2*q2rpf(ii+1) + 6*q2rpf(ii) )*dd
      !          ! print *,ii,medprod(ii)
      !      end do
      !
      !  !---    find first maximum in this distribution
      !      do ii = 2,nBins-2
      !          if (  medprod(ii) > max(medprod(ii-1),medprod(ii+1)) ) then
      !              imax = ii
      !              exit
      !          end if
      !      end do
      !      medprod(1:imax) = huge(1.0)
      !
      !  !---    now find minimum
      !     imax = minloc( medprod(:),dim=1 ) - 1
      !
      !     !do ii = imax,nBins-2
      !     !     if (  medprod(ii) < min(medprod(ii-1),medprod(ii+1)) ) then
      !     !         imax = ii
      !     !         exit
      !     !     end if
      !     ! end do
      !
      !
      !
      !    noise = rpf(imax)*inxny/max(1,winbin(imax))
      !     if (FFTW_DBG) print *,"Lib_FFTW3f::radialPowerSpectrum2d() info - using noise level ",imax," q= ",imax*q_max/nBins," rpf(q)= ",noise ," medprod ",medprod(imax)
      !
      !  !---    normalise and return
      !      rpfsum = 0          !   compute as sum out(q) d2q
      !      qrpfbar = 0         !   compute as sum |q| out(q) d2q, then normalises
      !      q2rpfbar = 0        !   compute as sum |q^2| out(q) d2q, then normalises
      !
      !      do ii = 0,imax
      !          rpfsum = rpfsum + rpf(ii)
      !          qrpfbar = qrpfbar + qrpf(ii)
      !          q2rpfbar = q2rpfbar + q2rpf(ii)
      !      end do
      !
      !       do ii = 0,nBins
      !          dd = inxny/max(1,winbin(ii))
      !          rpf(ii) = rpf(ii)*dd
      !      end do
      !
      !
      !      qrpfbar = qrpfbar/max(1.0d-16,rpfsum)
      !      q2rpfbar = q2rpfbar/max(1.0d-16,rpfsum)
      !      if (FFTW_DBG) print *,"Lib_FFTW3f::radialPowerSpectrum2d() dbg - sum over large ",rpfsum,qrpfbar,q2rpfbar

        !---    tidy up
            deallocate(winbin)


            return
        end subroutine radialPowerSpectrum2d





















        subroutine radialPowerFunction2d( in,q_min,q_max, rpf,rpfsum,qrpfbar,q2rpfbar , snr )
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      compute the radial power function between q_min and q_max
    !*      also return some simple integrals associated with the power function
            real(kind=real64),dimension(:,:),intent(inout)      ::      in
            real(kind=real64),intent(in)                        ::      q_min,q_max
            real(kind=real64),dimension(0:),intent(out)         ::      rpf             !   | F(q) |^2
            real(kind=real64),intent(out)                       ::      rpfsum          !   int rpf(q) d2q
            real(kind=real64),intent(out)                       ::      qrpfbar         !   int |q| rpf(q) d2q/int rpf(q) d2q
            real(kind=real64),intent(out)                       ::      q2rpfbar        !   int |q^2| rpf(q) d2q/int rpf(q) d2q
            real(kind=real64),intent(in),optional               ::      snr
            real(kind=real64),dimension(:,:),allocatable        ::      out

            integer             ::      Nx,Ny,Mx,My,nBins
            integer             ::      ix,iy,ii,nBinMiss , imax
            integer,dimension(:),allocatable                    ::     winbin
            real(kind=real64),dimension(:),allocatable          ::     qrpf,q2rpf , medprod
            real(kind=real64)   ::      qq,deltaqx,deltaqy,deltaq,ideltaq,noise,inxny,dd

!            logical             ::      ok

            Nx = size(in,dim=1)
            Ny = size(in,dim=2)

            Mx = int(Nx/2)        !   note real array gives complex FT with N/2 output values!
            My = int(Ny/2)        !

        !---    compute the power spectrum
            if (FFTW_DBG) then

                print *,"Nx,Ny   = ",Nx,Ny
                print *,"Mx,My   = ",Mx,My
            end if
            allocate(out(-Mx:Mx,-My:My))
            call dft_r2I_2d( in,out )


            if (FFTW_DBG) then
                print *,"Lib_FFTW3f::radialPowerFunction() dbg - proof of Parseval's theorem"
                print *,"Nx Ny   = ",Nx*Ny
                print *,"<in^2>  = ",sum(in**2)/(Nx*Ny)
                print *,"<out^2> = ",sum(out)  /(Nx*Ny)
            end if


            nBins = size(rpf)-1
            allocate( winbin(0:nBins) )

            allocate( qrpf(0:nBins) )
            allocate( q2rpf(0:nBins) )
            deltaqx = PI/Mx     !    2*PI/Nx
            deltaqy = PI/My     !    2*PI/Ny
            deltaq = (q_max-q_min)/nBins
            ideltaq = 1/deltaq
            inxny = 1.0d0/(Nx*Ny)
            !print *,"deltaqx,deltaqy,deltaq ",deltaqx,deltaqy,1/ideltaq

            winbin = 0
            rpf = 0
            rpfsum = 0          !   compute as sum out(q) d2q
            qrpfbar = 0         !   compute as sum |q| out(q) d2q, then normalises
            q2rpfbar = 0        !   compute as sum |q^2| out(q) d2q, then normalises
            nBinMiss = 0
            qrpf = 0
            q2rpf = 0
            do iy = -My,My
                do ix = -Mx,Mx
                    qq = sqrt( (ix*deltaqx)**2 + (iy*deltaqy)**2 )
                    !if (qq < q_min) cycle
                    !if (qq > q_max) cycle
                    rpfsum = rpfsum + out(ix,iy)
                    qrpfbar = qrpfbar + qq*out(ix,iy)
                    q2rpfbar = q2rpfbar + qq*qq*out(ix,iy)
                    ii = nint( (qq-q_min)*ideltaq )
                    if (ii*(nBins-ii)<0) then
                        nBinMiss = nBinMiss + 1
                        cycle
                    end if
                    winbin(ii) = winbin(ii) + 1
                    rpf(ii) = rpf(ii) + out(ix,iy)
                    qrpf(ii) = qrpf(ii) + qq*out(ix,iy)
                    q2rpf(ii) = q2rpf(ii) + qq*qq*out(ix,iy)
                end do
            end do
            qrpfbar = qrpfbar/max(1.0d-16,rpfsum)
            q2rpfbar = q2rpfbar/max(1.0d-16,rpfsum)
            if (FFTW_DBG) print *,"Lib_FFTW3f::radialPowerFunction() dbg - sum over all q ",rpfsum,qrpfbar,q2rpfbar
            if (FFTW_DBG) print *,"Lib_FFTW3f::radialPowerFunction() dbg - number bin misses ",nBinMiss,"/",(2*Mx+1)*(2*My+1)

       ! !---    find largest value of rpf
       !     qq = -huge(1.0)
       !     do ii = 0,nBins
       !         !print *,ii,winbin(ii),rpf(ii),qrpf(ii),q2rpf(ii)
       !         if (winbin(ii)>0) qq = max( qq, rpf(ii)/winbin(ii) )
       !     end do
       !     qq = qq/(Nx*Ny)
       !     if (FFTW_DBG) print *,"Lib_FFTW3f::radialPowerFunction() dbg - largest value of rpf() = ",qq
       !     noise = qq*1.0d-7
       !     if (present(snr)) then
       !         if (snr>0) noise = max(noise,qq/(snr*snr) )
       !     end if
       !    !print *,"Lib_FFTW3f::radialPowerFunction() info - using noise level ",noise

       !    do ii = 0,nBins
       !        dd = inxny/max(1,winbin(ii))
       !        rpf(ii) = rpf(ii)*dd
       !   !     qrpf(ii) = qrpf(ii)*dd
       !   !     q2rpf(ii) = q2rpf(ii)*dd
       !    end do
       !


       !----    find a good estimate for the level at which pixel noise takes over
       !    compute a 5 point median filter for qrpf and q2rpf
       !    then find the product/q^3.
            allocate( medprod(0:nBins) ) ; medprod = huge(1.0)
            do ii = 2,nBins-2
                qq = ii*q_max/nBins
                dd = 1/(144*qq*qq*qq)
                medprod(ii) = ( qrpf(ii-2) +  qrpf(ii+2) + 2*qrpf(ii-1) + 2*qrpf(ii+1) + 6*qrpf(ii) )*( q2rpf(ii-2) + q2rpf(ii+2) + 2*q2rpf(ii-1) + 2*q2rpf(ii+1) + 6*q2rpf(ii) )*dd
                ! print *,ii,medprod(ii)
            end do

        !---    find first maximum in this distribution
            imax = 0
            do ii = 2,nBins-2
                if (  medprod(ii) > max(medprod(ii-1),medprod(ii+1)) ) then
                    imax = ii
                    exit
                end if
            end do
            medprod(1:imax) = huge(1.0)

        !---    now find minimum
           imax = minloc( medprod(:),dim=1 ) - 1

           !do ii = imax,nBins-2
           !     if (  medprod(ii) < min(medprod(ii-1),medprod(ii+1)) ) then
           !         imax = ii
           !         exit
           !     end if
           ! end do



           noise = rpf(imax)*inxny/max(1,winbin(imax))
           if (FFTW_DBG) print *,"Lib_FFTW3f::radialPowerFunction() info - using noise level ",imax," q= ",imax*q_max/nBins," rpf(q)= ",noise ," medprod ",medprod(imax)

        !---    normalise and return
            rpfsum = 0          !   compute as sum out(q) d2q
            qrpfbar = 0         !   compute as sum |q| out(q) d2q, then normalises
            q2rpfbar = 0        !   compute as sum |q^2| out(q) d2q, then normalises

            do ii = 0,imax
                rpfsum = rpfsum + rpf(ii)
                qrpfbar = qrpfbar + qrpf(ii)
                q2rpfbar = q2rpfbar + q2rpf(ii)
            end do

             do ii = 0,nBins
                dd = inxny/max(1,winbin(ii))
                rpf(ii) = rpf(ii)*dd
            end do


            qrpfbar = qrpfbar/max(1.0d-16,rpfsum)
            q2rpfbar = q2rpfbar/max(1.0d-16,rpfsum)
            if (FFTW_DBG) print *,"Lib_FFTW3f::radialPowerFunction() dbg - sum over large ",rpfsum,qrpfbar,q2rpfbar

        !---    tidy up
            deallocate(winbin)


            return
        end subroutine radialPowerFunction2d





        subroutine radialPowerFunction3d( inout,q_min,q_max, rpf,rpfsum,qrpfbar,q2rpfbar,a_cell_in )
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      compute the radial power function between q_min and q_max
    !*      also return some simple integrals associated with the power function
            real(kind=real64),dimension(:,:,:),allocatable,intent(inout)    ::      inout
            real(kind=real64),intent(in)                        ::      q_min,q_max
            real(kind=real64),dimension(0:),intent(out)         ::      rpf             !   | F(q) |^2
            real(kind=real64),intent(out)                       ::      rpfsum          !   int rpf(q) d2q
            real(kind=real64),intent(out)                       ::      qrpfbar         !   int |q| rpf(q) d2q/int rpf(q) d2q
            real(kind=real64),intent(out)                       ::      q2rpfbar        !   int |q^2| rpf(q) d2q/int rpf(q) d2q

            real(kind=real64),dimension(3,3),intent(in),optional    ::      a_cell_in      !   optionally add size and shape of input cells, otherwise assume identity


            !real(kind=real64),dimension(:,:,:),allocatable      ::        out

            integer             ::      Nx,Ny,Nz,Mx,My,Mz,nBins
            integer             ::      ix,iy,iz,ii,nBinMiss , imax  
            real(kind=real64),dimension(:),allocatable          ::     winbin
            real(kind=real64),dimension(:),allocatable          ::     qrpf,q2rpf , medprod
            real(kind=real64)                           ::      qq,deltaq,ideltaq,noise,dd , invox , aa,i2s2 
            real(kind=real64),dimension(3,3)            ::      ia_cell
            real(kind=real64),dimension(3)              ::      qvec
             
            Nx = size(inout,dim=1)
            Ny = size(inout,dim=2)
            Nz = size(inout,dim=3)

            Mx = int(Nx/2)        !   note real array gives complex FT with N/2 output values!
            My = int(Ny/2)        !
            Mz = int(Nz/2)        !

        !---    compute the power spectrum
            if (FFTW_DBG) then
                print *,"Lib_FFTW3f::radialPowerFunction3d info - input data size Nx,Ny,Nz   = ",Nx,Ny,Nz
            end if


            call dft_r2I_3d_inplace_old( inout )

            nBins = size(rpf)-1
            allocate( winbin(0:nBins) )

            allocate( qrpf(0:nBins) )
            allocate( q2rpf(0:nBins) )

        !---    convert real space cell into vector q separation.
            ia_cell = reshape( (/1,0,0,0,1,0,0,0,1/),(/3,3/) )
            if (present(a_cell_in)) call inverse3Mat(a_cell_in,ia_cell)


            ia_cell(1:3,1) = ia_cell(1:3,1) * PI / Mx
            ia_cell(1:3,2) = ia_cell(1:3,2) * PI / My
            ia_cell(1:3,3) = ia_cell(1:3,3) * PI / Mz


            deltaq = (q_max-q_min)/nBins
            ideltaq = 1/deltaq
            invox = 1.0d0/( real(Nx,kind=real64)*real(Ny,kind=real64)*real(Nz,kind=real64) )

            print *,"Lib_FFTW3f::radialPowerFunction3d() info - deltaq          = ",deltaq
            aa = determinant3Mat(ia_cell)**(1/3.0d0)            !   typical voxel lengthscale
            print *,"Lib_FFTW3f::radialPowerFunction3d() info - FFTW pixel size = ",aa
            i2s2 = 1/(2*aa*aa)          !
            !mm = max(1,ceiling( 2*aa*ideltaq ))
            !print *,"Lib_FFTW3f::radialPowerFunction3d() info - q intervals contributed per pixel ",mm
            !allocate(weight(-mm:mm))

          !  print *,""
          !  do ix = -Mx,Mx
          !      write (*,fmt='(i6,10g16.6)') ix,inout(ix,0,0),inout(ix,My/4,0),inout(ix,My/2,0),inout(ix,My/4,Mz/4),inout(ix,My/2,Mz/4)
          !  end do
          !  print *,""


            winbin = 0
            rpf = 0
            rpfsum = 0          !   compute as sum out(q) d2q
            qrpfbar = 0         !   compute as sum |q| out(q) d2q, then normalises
            q2rpfbar = 0        !   compute as sum |q^2| out(q) d2q, then normalises
            nBinMiss = 0
            qrpf = 0
            q2rpf = 0
            do iz = -Mz,Mz
                do iy = -My,My
                
                    call progressBar( (iy+My+1)+(iz+Mz)*(2*My+1) , (2*My+1)+(2*Mz)*(2*My+1) )
                    do ix = -Mx,Mx

                    
                    
                    
                        qvec(1:3) = ia_cell(1:3,1)*ix + ia_cell(1:3,2)*iy + ia_cell(1:3,3)*iz
                        qq = norm2( qvec )

                        dd = inout(ix,iy,iz)

                        rpfsum = rpfsum + dd
                        qrpfbar = qrpfbar + qq*dd
                        q2rpfbar = q2rpfbar + qq*qq*dd
                        ii = nint( (qq-q_min)*ideltaq )

                        if (ii*(nBins-ii)<0) then
                            nBinMiss = nBinMiss + 1
                            cycle
                        end if


                    !!---   use voxel weightings to antialias
                    !    weight = 0.0
                    !    do jj = max(0,ii-mm),min(nBins,ii+mm)
                    !       ww = ( qq - q_min - jj*deltaq )
                    !       weight( jj-ii ) =  exp( - ww*ww*i2s2 )
                    !   end do
                    !   ww = 1/max(1.0d-16,sum(weight))
                    !   weight = weight*ww
                    !
                    !
                    !    do jj = max(0,ii-mm),min(nBins,ii+mm)
                    !       ww = weight( jj-ii )
                    !        winbin(ii) = winbin(ii) + ww
                    !        rpf(ii) = rpf(ii) + ww*dd
                    !        qrpf(ii) = qrpf(ii) + ww*qq*dd
                    !        q2rpf(ii) = q2rpf(ii) + ww*qq*qq*dd
                    !
                    !    end do

                    !---    no antialiasing
                         winbin(ii) = winbin(ii) + 1
                         rpf(ii) = rpf(ii) + dd
                         qrpf(ii) = qrpf(ii) + qq*dd
                         q2rpf(ii) = q2rpf(ii) + qq*qq*dd

                    end do
                end do
            end do
            qrpfbar = qrpfbar/max(1.0d-16,rpfsum)
            q2rpfbar = q2rpfbar/max(1.0d-16,rpfsum)
            if (FFTW_DBG) print *,"Lib_FFTW3f::radialPowerFunction3d() dbg - sum over all q ",rpfsum,qrpfbar,q2rpfbar
            if (FFTW_DBG) print *,"Lib_FFTW3f::radialPowerFunction3d() dbg - bin misses     ",nBinMiss*invox*100.0,"%"


       !----    find a good estimate for the level at which pixel noise takes over
       !    compute a 5 point median filter for qrpf and q2rpf
       !    then find the product/q^3.
            allocate( medprod(0:nBins) ) ; medprod = huge(1.0)
            do ii = 2,nBins-2
                qq = q_min + ii*deltaq
                dd = 1/(144*qq*qq*qq)
                medprod(ii) = ( qrpf(ii-2) + 2*qrpf(ii-1) + 6*qrpf(ii) + 2*qrpf(ii+1) + qrpf(ii+2) )                &
                            *( q2rpf(ii-2) + 2*q2rpf(ii-1) + 6*q2rpf(ii) + 2*q2rpf(ii+1) + q2rpf(ii+2) )*dd
            end do

        !---    find first maximum in this distribution
            imax = 0
            do ii = 2,nBins-2
                if (  medprod(ii) > max(medprod(ii-1),medprod(ii+1)) ) then
                    imax = ii
                    exit
                end if
            end do
            medprod(1:imax) = huge(1.0)

        !---    now find minimum
           imax = minloc( medprod(:),dim=1 ) - 1

           noise = rpf(imax)*invox/max(1.0d0,winbin(imax))
           if (FFTW_DBG) print *,"Lib_FFTW3f::radialPowerFunction3d() info - using noise level ",imax," q= ",imax*q_max/nBins," rpf(q)= ",noise ," medprod ",medprod(imax)," winbin ",winbin(imax)
           deallocate( medprod )


        !---    normalise and return
            rpfsum = 0          !   compute as sum out(q) d2q
            qrpfbar = 0         !   compute as sum |q| out(q) d2q, then normalises
            q2rpfbar = 0        !   compute as sum |q^2| out(q) d2q, then normalises

            do ii = 0,imax
                rpfsum = rpfsum + rpf(ii)
                qrpfbar = qrpfbar + qrpf(ii)
                q2rpfbar = q2rpfbar + q2rpf(ii)
            end do

            invox = sqrt(invox)
             do ii = 0,nBins
                dd = invox/max(1.0d0,winbin(ii))
                rpf(ii) = rpf(ii)*dd
            end do


            qrpfbar = qrpfbar/max(1.0d-16,rpfsum)
            q2rpfbar = q2rpfbar/max(1.0d-16,rpfsum)
            if (FFTW_DBG) print *,"Lib_FFTW3f::radialPowerFunction3d() dbg - sum over large ",rpfsum,qrpfbar,q2rpfbar

        !---    tidy up
            deallocate(winbin)


            return
        end subroutine radialPowerFunction3d


        subroutine radialPowerFunction3d_old( in,q_min,q_max, rpf,rpfsum,qrpfbar,q2rpfbar,a_cell_in )
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      compute the radial power function between q_min and q_max
    !*      also return some simple integrals associated with the power function
            real(kind=real64),dimension(:,:,:),intent(inout)    ::      in
            real(kind=real64),intent(in)                        ::      q_min,q_max
            real(kind=real64),dimension(0:),intent(out)         ::      rpf             !   | F(q) |^2
            real(kind=real64),intent(out)                       ::      rpfsum          !   int rpf(q) d2q
            real(kind=real64),intent(out)                       ::      qrpfbar         !   int |q| rpf(q) d2q/int rpf(q) d2q
            real(kind=real64),intent(out)                       ::      q2rpfbar        !   int |q^2| rpf(q) d2q/int rpf(q) d2q

            real(kind=real64),dimension(3,3),intent(in),optional    ::      a_cell_in      !   optionally add size and shape of input cells, otherwise assume identity


            real(kind=real64),dimension(:,:,:),allocatable      ::        out

            integer             ::      Nx,Ny,Nz,Mx,My,Mz,nBins
            integer             ::      ix,iy,iz,ii,nBinMiss , imax
            integer,dimension(:),allocatable                    ::     winbin
            real(kind=real64),dimension(:),allocatable          ::     qrpf,q2rpf , medprod
            real(kind=real64)   ::      qq,deltaq,ideltaq,noise,inxnynz,dd
            real(kind=real64),dimension(3,3)            ::      ia_cell
            real(kind=real64),dimension(3)              ::      pionm,qvec

            Nx = size(in,dim=1)
            Ny = size(in,dim=2)
            Nz = size(in,dim=3)

            Mx = int(Nx/2)        !   note real array gives complex FT with N/2 output values!
            My = int(Ny/2)        !
            Mz = int(Nz/2)        !

        !---    compute the power spectrum
            if (FFTW_DBG) then
                print *,"Lib_FFTW3f::radialPowerFunction3d info - input data size Nx,Ny,Nz   = ",Nx,Ny,Nz
                !print *,"Lib_FFTW3f::radialPowerFunction3d info - FT data size    Mx,My,Mz   = ",Mx,My,Mz
            end if

            allocate(out(-Mx:Mx,-My:My,-Mz:Mz))
            call dft_r2I_3d( in,out )

!            if (FFTW_DBG) then
!                print *,"x = 0 plane"
!                do iz = -Mz,Mz
!                    write (*,fmt='(1000f14.4)') out(0,:,iz)
!                end do
!                print *,"y = 0 plane"
!                do ix = -Mx,Mx
!                    write (*,fmt='(1000f14.4)') out(ix,0,:)
!                end do
!                print *,"z = 0 plane"
!                do iy = -My,My
!                    write (*,fmt='(1000f14.4)') out(:,iy,0)
!                end do
!
!            end if
!
!




           if (FFTW_DBG) then
               print *,"Lib_FFTW3f::radialPowerFunction3d() dbg - proof of Parseval's theorem"
               !print *,"Nx Ny Nz= ",Nx*Ny*Nz
               print *,"<in^2>  = ",sum(in**2)/(Nx*Ny*Nz)
               print *,"<out^2> = ",sum(out)  /(Nx*Ny*Nz)
           end if


            nBins = size(rpf)-1
            allocate( winbin(0:nBins) )

            allocate( qrpf(0:nBins) )
            allocate( q2rpf(0:nBins) )

        !---    convert real space cell into vector q separation.
            ia_cell = reshape( (/1,0,0,0,1,0,0,0,1/),(/3,3/) )
            if (present(a_cell_in)) call inverse3Mat(a_cell_in,ia_cell)

            pionm(1:3) = PI/(/Mx,My,Mz/)










            !   deltaqx = PI/Mx     !    2*PI/Nx
            !   deltaqy = PI/My     !    2*PI/Ny
            !   deltaqz = PI/Mz     !    2*PI/Ny


            deltaq = (q_max-q_min)/nBins
            ideltaq = 1/deltaq


            winbin = 0
            rpf = 0
            rpfsum = 0          !   compute as sum out(q) d2q
            qrpfbar = 0         !   compute as sum |q| out(q) d2q, then normalises
            q2rpfbar = 0        !   compute as sum |q^2| out(q) d2q, then normalises
            nBinMiss = 0
            qrpf = 0
            q2rpf = 0
            do iz = -Mz,Mz
                do iy = -My,My
                    do ix = -Mx,Mx

                        qvec(1:3) = ia_cell(1:3,1)*ix*pionm(1) + ia_cell(1:3,2)*iy*pionm(2) + ia_cell(1:3,3)*iz*pionm(3)
                        qq = norm2( qvec )

!                        qq = sqrt( (ix*deltaqx)**2 + (iy*deltaqy)**2 + (iz*deltaqz)**2 )

                        rpfsum = rpfsum + out(ix,iy,iz)
                        qrpfbar = qrpfbar + qq*out(ix,iy,iz)
                        q2rpfbar = q2rpfbar + qq*qq*out(ix,iy,iz)
                        ii = nint( (qq-q_min)*ideltaq )

                        if (ii*(nBins-ii)<0) then
                            nBinMiss = nBinMiss + 1
                            cycle
                        end if

                        winbin(ii) = winbin(ii) + 1
                        rpf(ii) = rpf(ii) + out(ix,iy,iz)
                        qrpf(ii) = qrpf(ii) + qq*out(ix,iy,iz)
                        q2rpf(ii) = q2rpf(ii) + qq*qq*out(ix,iy,iz)
                    end do
                end do
            end do
            qrpfbar = qrpfbar/max(1.0d-16,rpfsum)
            q2rpfbar = q2rpfbar/max(1.0d-16,rpfsum)
            if (FFTW_DBG) print *,"Lib_FFTW3f::radialPowerFunction3d() dbg - sum over all q ",rpfsum,qrpfbar,q2rpfbar
            if (FFTW_DBG) print *,"Lib_FFTW3f::radialPowerFunction3d() dbg - number bin misses ",nBinMiss,"/",(2*Mx+1)*(2*My+1)*(2*Mz+1)


       !----    find a good estimate for the level at which pixel noise takes over
       !    compute a 5 point median filter for qrpf and q2rpf
       !    then find the product/q^3.
            allocate( medprod(0:nBins) ) ; medprod = huge(1.0)
            do ii = 2,nBins-2
                qq = q_min + ii*deltaq
                dd = 1/(144*qq*qq*qq)
                medprod(ii) = ( qrpf(ii-2) + 2*qrpf(ii-1) + 6*qrpf(ii) + 2*qrpf(ii+1) + qrpf(ii+2) )                &
                            *( q2rpf(ii-2) + 2*q2rpf(ii-1) + 6*q2rpf(ii) + 2*q2rpf(ii+1) + q2rpf(ii+2) )*dd
            end do

        !---    find first maximum in this distribution
            imax = 0
            do ii = 2,nBins-2
                if (  medprod(ii) > max(medprod(ii-1),medprod(ii+1)) ) then
                    imax = ii
                    exit
                end if
            end do
            medprod(1:imax) = huge(1.0)

        !---    now find minimum
           imax = minloc( medprod(:),dim=1 ) - 1

           !do ii = imax,nBins-2
           !     if (  medprod(ii) < min(medprod(ii-1),medprod(ii+1)) ) then
           !         imax = ii
           !         exit
           !     end if
           ! end do


           inxnynz = (1.0d0/Nx)*(1.0d0/Ny)*(1.0d0/Nz)           !   why? because Nx Ny Nz can be > 3.2B.
           noise = rpf(imax)*inxnynz/max(1,winbin(imax))
           if (FFTW_DBG) print *,"Lib_FFTW3f::radialPowerFunction3d() info - using noise level ",imax," q= ",imax*q_max/nBins," rpf(q)= ",noise ," medprod ",medprod(imax)

        !---    normalise and return
            rpfsum = 0          !   compute as sum out(q) d2q
            qrpfbar = 0         !   compute as sum |q| out(q) d2q, then normalises
            q2rpfbar = 0        !   compute as sum |q^2| out(q) d2q, then normalises

            do ii = 0,imax
                rpfsum = rpfsum + rpf(ii)
                qrpfbar = qrpfbar + qrpf(ii)
                q2rpfbar = q2rpfbar + q2rpf(ii)
            end do

             do ii = 0,nBins
                dd = inxnynz/max(1,winbin(ii))
                rpf(ii) = rpf(ii)*dd
            end do


            qrpfbar = qrpfbar/max(1.0d-16,rpfsum)
            q2rpfbar = q2rpfbar/max(1.0d-16,rpfsum)
            if (FFTW_DBG) print *,"Lib_FFTW3f::radialPowerFunction3d() dbg - sum over large ",rpfsum,qrpfbar,q2rpfbar

        !---    tidy up
            deallocate(winbin)


            return
        end subroutine radialPowerFunction3d_old


        subroutine radialPowerFunction_( in,qmax , rpf,rpfsum,qrpfbar,q2rpfbar )
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      INCORRECT IMPLEMENTATION. WHY DON'T I DELETE IT???!!!
            real(kind=real64),dimension(:,:),intent(inout)      ::      in
            integer,intent(in)                                  ::      qmax
            real(kind=real64),dimension(0:),intent(out)         ::      rpf
            real(kind=real64),intent(out)                       ::      rpfsum,qrpfbar,q2rpfbar         !   int rpf(q) d2q, int |q| rpf(q) d2q/int rpf(q) d2q , int |q^2| rpf(q) d2q/int rpf(q) d2q

            real(kind=real64),dimension(:,:),allocatable        ::        out

            integer             ::      Nx,Ny,Mx,My,nBins
            integer             ::      ix,iy,ii,qmax2,jj
            real(kind=real64),dimension(:,:),allocatable        ::      weight
            real(kind=real64),dimension(:),allocatable          ::      rpf_tmp,winbin
            integer,dimension(:),allocatable                    ::      bin
            real(kind=real64)   ::      qq,ww,deltaq,ideltaq

            real(kind=real64)   ::      i2s2 = 1/(2*0.333d0**2)
            Nx = size(in,dim=1)
            Ny = size(in,dim=2)

            Mx = int(Nx/2)        !   note real array gives complex FT with N/2 output values!
            My = int(Ny/2)        !   note real array gives complex FT with N/2 output values!

        !---    compute the power spectrum
            allocate(out(-Mx:Mx,-My:My))
            call dft_r2I_2d( in,out )


        !---    find indexing of radius squared to bin
            !nBins = max(Mx,My)      !   not going further than half the image size.
            !qmax =  max(Mx,My) ! sqrt( real(Mx*Mx + My*My) )
            nBins = size(rpf)-1
            deltaq = qmax / nBins
            ideltaq = nBins / qmax
            qmax2 = qmax*qmax


        !---    proof of Parseval's theorem.
        !       sum_i | x_i^2 | = (1/2pi) int_q | X(q) |^2 dq
        !   note for my Fourier transform qx = 2 pi kx / Nx
        !   so dq = 2 pi/Nx dk
        !
            print *,"Nx Ny   = ",Nx*Ny
            print *,"<in^2>  = ",sum(in**2)/(Nx*Ny)
            print *,"<out^2> = ",sum(out)  /(Nx*Ny)


            allocate(bin(0:qmax2))
            allocate(weight(-2:2,0:qmax2))
            bin = -1
            weight = 0
            do ix = 0,nBins
                do iy = 0,ix
                    ii = ix*ix + iy*iy

                    if (ii>qmax2) cycle
                    if (bin(ii)>=0) cycle
                    qq = sqrt( real(ii,kind=real64) )*ideltaq
                    jj = nint( qq )
                    bin(ii) = jj

                    qq = qq-jj                                      !   distance of separation to bin
                    if (jj>1) weight(-2,ii) = exp( -i2s2*(qq-2)*(qq-2) )
                    if (jj>0) weight(-1,ii) = exp( -i2s2*(qq-1)*(qq-1) )
                    weight( 0,ii) = exp( -i2s2* qq   * qq )
                    weight( 1,ii) = exp( -i2s2*(qq+1)*(qq+1) )
                    weight( 2,ii) = exp( -i2s2*(qq+2)*(qq+2) )
                    ww = sum( weight(:,ii) )
                    ww = 1/ww
                    weight(:,ii) = weight(:,ii)*ww



                end do
            end do

            !print *,"nbins ",nbins, maxval(bin),qmax


        !---    construct radial power function
            allocate( rpf_tmp(-2:nBins+2) )
            allocate( winbin(-2:nBins+2) )
            rpf_tmp = 0 ; winbin = 0.0d0
            rpfsum = 0          !   compute as sum out(q) d2q
            qrpfbar = 0         !   compute as sum |q| out(q) d2q, then normalise
            q2rpfbar = 0

            do iy = -My,My
                do ix = -Mx,Mx
                    ii = ix*ix + iy*iy
                    if (ii > qmax2) cycle
                    jj = bin(ii)
                    if (jj == -1) cycle
                    rpf_tmp(jj-2:jj+2) = rpf_tmp(jj-2:jj+2) + weight(-2:2,ii)*out(ix,iy)
                    winbin(jj-2:jj+2)  = winbin(jj-2:jj+2) + weight(-2:2,ii)
                    qq = sqrt( real(ii,kind=real64) )       !   |q|
                    rpfsum = rpfsum + out(ix,iy)
                    qrpfbar = qrpfbar + qq*out(ix,iy)
                    q2rpfbar = q2rpfbar + qq*qq*out(ix,iy)
                end do
            end do
            q2rpfbar = q2rpfbar/max(1.0d-16,rpfsum)
            qrpfbar = qrpfbar/max(1.0d-16,rpfsum)


        !---    normalise and return
            do jj = 0,nBins
                ww = 0.0d0
                if (winbin(jj)>1.0d-12) ww = 1/( Nx*Ny*winbin(jj) )
                rpf(jj) = rpf_tmp(jj)*ww
            end do

        !---    tify up
            deallocate(winbin)
            deallocate(rpf_tmp)
            deallocate(weight)
            deallocate(bin)
            deallocate(out)

            return
        end subroutine radialPowerFunction_




        pure subroutine inverse3Mat(M,N)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      monadic operator
    !*      finds the inverse of a general three matrix
            real(kind=real64),dimension(3,3),intent(in)       ::  M
            real(kind=real64),dimension(3,3),intent(out)      ::  N
            real(kind=real64)            ::      idd
            real(kind=real64),dimension(3,3),parameter        ::      &
            IDENTITY3MAT = reshape( (/ 1.0d0,0.0d0,0.0d0 , 0.0d0,1.0d0,0.0d0 , 0.0d0,0.0d0,1.0d0 /) &
                                   ,(/ 3,3 /) )

            idd = determinant3Mat(M)
            if (abs(idd) < tiny(1.0d0)) then
                N = IDENTITY3MAT
                return
            end if
            idd = 1.0/idd

            N(1,1)   = ( M(2,2)*M(3,3) - M(2,3)*M(3,2) ) * idd
            N(2,1)   = ( M(2,3)*M(3,1) - M(2,1)*M(3,3) ) * idd
            N(3,1)   = ( M(2,1)*M(3,2) - M(2,2)*M(3,1) ) * idd

            N(1,2)   = ( M(1,3)*M(3,2) - M(1,2)*M(3,3) ) * idd
            N(2,2)   = ( M(1,1)*M(3,3) - M(1,3)*M(3,1) ) * idd
            N(3,2)   = ( M(1,2)*M(3,1) - M(1,1)*M(3,2) ) * idd

            N(1,3)   = ( M(1,2)*M(2,3) - M(1,3)*M(2,2) ) * idd
            N(2,3)   = ( M(1,3)*M(2,1) - M(1,1)*M(2,3) ) * idd
            N(3,3)   = ( M(1,1)*M(2,2) - M(1,2)*M(2,1) ) * idd

            return
        end subroutine inverse3Mat

        pure function determinant3Mat(M) result(d)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      returns the determinant of M
            real(kind=real64),dimension(3,3),intent(in)      ::      M
            real(kind=real64)                                ::      d
            real(kind=real64),dimension(9)       ::      dd
            dd(1) = M(1,1)*( M(2,2)*M(3,3) - M(2,3)*M(3,2) )
            dd(2) = M(1,2)*( M(2,3)*M(3,1) - M(2,1)*M(3,3) )
            dd(3) = M(1,3)*( M(2,1)*M(3,2) - M(2,2)*M(3,1) )
            dd(4) = M(2,1)*( M(3,2)*M(1,3) - M(3,3)*M(1,2) )
            dd(5) = M(2,2)*( M(3,3)*M(1,1) - M(3,1)*M(1,3) )
            dd(6) = M(2,3)*( M(3,1)*M(1,2) - M(3,2)*M(1,1) )
            dd(7) = M(3,1)*( M(1,2)*M(2,3) - M(1,3)*M(2,2) )
            dd(8) = M(3,2)*( M(1,3)*M(2,1) - M(1,1)*M(2,3) )
            dd(9) = M(3,3)*( M(1,1)*M(2,2) - M(1,2)*M(2,1) )
            d = (1.0d0/3.0d0) * sum(dd)
            return
        end function determinant3Mat




    end module Lib_FFTW3f

!   gcc -c ${MYF90LIB}/Lib_Greyscale.c -lpng ; gfortran -ffree-line-length-256 -c ${MYF90LIB}/Lib_Png.f90 testFFTWf.f90 -lfftw3 -I/usr/include ; gfortran Lib_Greyscale.o  Lib_Png.o testFFTWf.o -o testLib_FFTW3f.exe -lpng -lfftw3

!     program testLib_FFTW3f
! !---^^^^^^^^^^^^^^^^^^^^^^
!         use Lib_FFTW3f
!         use Lib_Png
!         use iso_fortran_env
!         implicit none
!
!         character(len=256)                                  ::      dummy
!         real(kind=real64),dimension(:,:),allocatable        ::      img_in,img_out
!
!
!         integer             ::      Nx,Ny,Mx,My
!         integer             ::      ix  !,iy,ii
!         real(kind=real64)   ::      aa , qq,rpfqsum,qmax   ! , rr
!         real(kind=real64),dimension(:),allocatable  ::  rpf
!         !integer,dimension(:),allocatable  ::  nrpf
!         integer             ::      nBins
!
!          call get_command_argument(1,dummy)
!          call readPng( dummy,img_in )
!
!        ! dummy = "test.png"
!        ! Nx = 512 ; Ny = 512
!        ! allocate(img_in(Nx,Ny))
!        ! do iy = 1,Ny
!        !     aa = 2*3.141592654d0 * iy / 64
!        !     img_in(:,iy) = cos(aa)
!        ! end do
!        ! call writePng( dummy, img_in,normalise = .true. )
!
!
!         Nx = size(img_in,dim=1)
!         Ny = size(img_in,dim=2)
!
!         Mx = int( Nx/2 )
!         My = int( Ny/2 )
!
!         allocate(img_out(-Mx:Mx,-My:My))
!
!         call FFT2d( img_in,img_out )
!
!         img_out = log( max(1.0d-12,img_out) )
!
!         call writePng(trim(dummy)//".ft.png",img_out,normalise=.true.)
!
!
!         print *,"radial power function"
!         qmax =  max(Mx,My) !sqrt( real(Mx*Mx + My*My) )
!         nBins = 2*ceiling( qmax )
!         allocate(rpf(0:nBins ))
!         aa = sum(img_in)/(Nx*Ny)
!         img_in = img_in - aa
!         aa = sum(img_in**2)/(Nx*Ny)
!         print *,"<in^2> = ",aa
!         call radialPowerFunction( img_in, rpf )
!         rpfqsum = 0.0d0
!         write(*,fmt='(a8,100a24)') "# ","|q|","|r|","C(q)","2 pi q C(q)"
!         do ix = 0,nBins
!             !aa = ix*real(nBins)/(nBins-1)
!             qq = ix*qmax/nBins
!             aa = 2*3.141592654d0 * qq * rpf(ix)
!
!
!             if (rpf(ix)>1e-12) then
!                 if (ix == 0) then
!                     write(*,fmt='(i8,100f24.12)') ix,qq,0.0d0,rpf(ix),aa
!                 else
!                     write(*,fmt='(i8,100f24.12)') ix,qq,max(Nx,Ny)/qq,rpf(ix),aa
!                 end if
!             end if
!             rpfqsum = rpfqsum + aa
!             if (mod(ix,2)==1) rpfqsum = rpfqsum + aa            !   1-4-1 rule is Simpson's rule.
!         end do
!         rpfqsum = rpfqsum * 2*qmax/(3*nBins)
!         print *,"integral 2 pi rpf(q) = ",rpfqsum
!
!
!
! !
! !         deallocate(rpf)
! !         nBins = min( Nx,Ny )/4
! !         nBins = nBins*nBins
! !         allocate(rpf(0:nBins))
! !         allocate(nrpf(0:nBins))
! !         rpf = 0  ; nrpf = 0
! !         do iy = 1,Ny
! !             do ix = 1,Mx*2-1
! !                 ii = ( ix - Mx )**2 + ( iy - Ny/2 )**2
! !                 if (ii>nBins) cycle
! !                 nrpf( ii ) = nrpf( ii ) + 1
! !                 aa = exp( img_out(ix,iy) )
! !                 rpf( ii ) = rpf( ii ) + aa
! !             end do
! !         end do
! !
!        !print *,""
!        !do ii = 0,nBins
!        !    if (nrpf(ii)>0) &
!        !    print *,ii,sqrt( real(ii) ),rpf(ii)/nrpf(ii)
!        !end do
!
!
!
!         print *,""
!         print *,"done"
!         print *,""
!     end program testLib_FFTW3f
!
!