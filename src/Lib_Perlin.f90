

    module Lib_Perlin2dPlane
!---^^^^^^^^^^^^^^^^^^^^^^^^
!*      produce Perlin noise on a 2d plane

        use iso_fortran_env
        implicit none            
        private
                      
        public      ::      perlinNoise
        
        logical,private     ::      LIB_PERLIN_DBG = .true.       
           
        integer,parameter,private   ::      LIB_PERLIN_LININT = 0
        integer,parameter,private   ::      LIB_PERLIN_GAUSSINT = 1
        integer,parameter,private   ::      LIB_PERLIN_SMOOTHINT = 2
        integer,public              ::      LIB_PERLIN_INTERPOLATION = LIB_PERLIN_GAUSSINT
        
        
    contains
!---^^^^^^^^

        subroutine perlinNoise( a,Nx,Ny, lengthscale , f )
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      produce Perlin noise on a 2d PBC plane grid 0:Nx-1, 0:Ny-1
    !*      with 2d lattice vectors a
    !*      using the given lengthscale 
    
            real(kind=real64),dimension(2,2),intent(in)         ::      a
            integer,intent(in)                                  ::      Nx,Ny
            real(kind=real64),intent(in)                        ::      lengthscale
            real(kind=real64),dimension(0:,0:),intent(out)      ::      f
            
            
            integer         ::      Mx,My           !   coarse grain cell size
            
            integer         ::      ix,iy
            integer         ::      jx,jy,jjx,jjy
            integer         ::      kx,ky
        
            real(kind=real64),dimension(:,:,:),allocatable      ::      coarseGrad
            real(kind=real64),dimension(:,:),allocatable        ::      coarseDP
            
            real(kind=real64),dimension(2)      ::      gg,dx,alpha,xx
            real(kind=real64)                   ::      dd,x1,x2,aa,ww,wsum
            real(kind=real64),dimension(2,2)    ::      bb
            
                
        !---    first establish the number of coarse grain points
            Mx = ceiling( Nx*norm2(a(:,1)) / lengthscale )
            My = ceiling( Ny*norm2(a(:,2)) / lengthscale )
            bb(1:2,1) = a(1:2,1)*Nx/Mx
            bb(1:2,2) = a(1:2,2)*Ny/My
            aa = 1.0d0/sqrt( abs( bb(1,1)*bb(2,2) - bb(2,1)*bb(1,2) ) )              !    inverse lengthscale of coarse lattice
            bb = bb*aa            
            if (LIB_PERLIN_DBG) then
                print *,"Lib_Perlin2dPlane::perlinNoise info - fine grain points   ",Nx,Ny
                print *,"Lib_Perlin2dPlane::perlinNoise info - lengthscale         ",lengthscale
                write(*,fmt='(a,2f12.5,a,2f12.5)') " Lib_Perlin2dPlane::perlinNoise info - fine lattice        ",a(1:2,1),",",a(1:2,2)
                print *,"Lib_Perlin2dPlane::perlinNoise info - coarse grain points ",Mx,My
                write(*,fmt='(a,2f12.5,a,2f12.5)') " Lib_Perlin2dPlane::perlinNoise info - coarse lattice      ",bb(1:2,1),",",bb(1:2,2)            
            end if
            allocate(coarseGrad(2,0:Mx,0:My ))
            allocate(coarseDP(0:Mx,0:My ))
            
        
    !---    select random unit vectors on coarse grid     
    !
    !        /
    !       o            o
    !                   /         the gradient is evaluated at o
    !                             find the dot product between gradient and 
    !                             vector to the point x
    !                     
    !                   \  
    !       o--          o
         
            do jy = 0,My-1
                do jx = 0,Mx-1
                    do
                        call random_number(gg)
                        gg = gg*2 - 1
                        dd = norm2(gg)
                        if (dd*(1-dd)>0) exit
                    end do
                    coarseGrad(1:2,jx,jy) = gg(1:2)*sqrt(0.125d0)/dd
                end do
            end do
            coarseGrad(1:2,:,My)  = coarseGrad(1:2,:,0)
            coarseGrad(1:2,Mx,:)  = coarseGrad(1:2,0,:)
            coarseGrad(1:2,Mx,My) = coarseGrad(1:2,0,0)
                
                
                
        !---    find dot product between gradient and distance vector to points on offset mesh
        !
        !
        !       x------------x
        !       |            |        the gradient is evaluated at o
        !       |      /     |        find the dot product between gradient and 
        !       |     o      |        vector to the point x
        !       |            |
        !       |            |
        !       x------------x
        !
                
            coarseDP(:,:) = 0.0d0
            do iy = 0,My-1                                                 !   (ix,iy) are the x- points
                do ix = 0,Mx-1                     
                    do ky = 0,1
                        jy = iy + ky          !   (jx,jy) are the o- points
                        do kx = 0,1
                            jx = ix + kx  
            
                            gg = coarseGrad(1:2,jx,jy) 
                            dx(1:2) = bb(1:2,1)*(kx-0.5d0) + bb(1:2,2)*(ky-0.5d0)
                            
                            coarseDP(ix,iy) = coarseDP(ix,iy) + dx(1)*gg(1) + dx(2)*gg(2) 
                        end do
                    end do
                end do
            end do
            coarseDP(:,My) = coarseDP(:,0)
            coarseDP(Mx,:) = coarseDP(0,:)
            coarseDP(Mx,My) = coarseDP(0,0)
!            print *,"maxmin ",maxval(coarseDP),minval(coarseDP)
                    
            
        !---    interpolate the dot products onto the fine grid      
            alpha(1) = real(Mx)/Nx
            alpha(2) = real(My)/Ny
            
            do ky = 0,Ny - 1
                do kx = 0,Nx - 1
                
                    xx(1:2) = (/kx,ky/)*alpha(1:2)
                    ix = int(xx(1)) ; iy = int(xx(2))
    
    
                        !print *,kx,ky,xx,ix,iy,gg    
                    dd = coarseDP(ix  ,iy  )        
                    if (LIB_PERLIN_INTERPOLATION == LIB_PERLIN_LININT) then
                                       
                        gg(1) = xx(1) - ix ; gg(2) = xx(2) - iy
                        dd = coarseDP(ix  ,iy  )*(1-gg(1))*(1-gg(2))                 &
                           + coarseDP(ix+1,iy  )*(  gg(1))*(1-gg(2))                 &
                           + coarseDP(ix  ,iy+1)*(1-gg(1))*(  gg(2))                 &
                           + coarseDP(ix+1,iy+1)*(  gg(1))*(  gg(2)) 
                           
                    else if (LIB_PERLIN_INTERPOLATION == LIB_PERLIN_GAUSSINT) then               
                        
                        dd = 0.0d0 ; wsum = 0.0d0
                        do jy = iy-2,iy+3
                            jjy = mod(jy + My,My)
                            do jx = ix-2,ix+3
                                jjx = mod(jx + Mx,Mx)
                                ww = (jx-xx(1))*(jx-xx(1)) + (jy-xx(2))*(jy-xx(2))      !   distance squared
                                ww = exp( -  ww  )                          !   gaussian with width 1/2 
                                dd = dd + ww*coarseDP(jjx,jjy)
                                wsum = wsum + ww
!                                  if ( (kx==100).and.(ky==100) ) then
!                                      print *,ix,iy,jx,jy,(jx-ix)*(jx-ix) + (jy-iy)*(jy-iy),ww,dd,wsum
!                                  end if
                            end do
                        end do
                        dd = 1.4 *dd/wsum
                    
                    else if (LIB_PERLIN_INTERPOLATION == LIB_PERLIN_SMOOTHINT) then               
                                 
                        gg(1) = xx(1) - ix ; gg(2) = xx(2) - iy
                        
                        !gg = gg*gg*gg*(10 - 15*gg + 6*gg*gg)
                        gg = gg*gg*(3-2*gg)
                        dd = coarseDP(ix  ,iy  )*(1-gg(1))*(1-gg(2))                 &
                           + coarseDP(ix+1,iy  )*(  gg(1))*(1-gg(2))                 &
                           + coarseDP(ix  ,iy+1)*(1-gg(1))*(  gg(2))                 &
                           + coarseDP(ix+1,iy+1)*(  gg(1))*(  gg(2)) 
                    
                    end if
                    f(kx,ky) = dd 
                
                end do
            end do
        
            return
        end subroutine perlinNoise
            
    end module  Lib_Perlin2dPlane                   
        
!        
!!    gcc -c ${MYF90LIB}/Lib_Greyscale.c -lpng16 ; gfortran -ffree-line-length-256 -c ${MYF90LIB}/Lib_Png.f90 ${MYF90LIB}/Lib_Perlin.f90 ; gfortran Lib_Greyscale.o  Lib_Png.o Lib_Perlin.o -o testLib_Perlin.exe -lpng16    
!        
!    program testLib_Perlin
!!---^^^^^^^^^^^^^^^^^^^^^^
!        use Lib_Perlin2dPlane
!        use Lib_Png
!        use iso_fortran_env
!        implicit none
!            
!        
!        integer         ::      Nx,Ny,Mx,My
!        logical         ::      pbc
!        real(kind=real64),dimension(:,:),allocatable        ::      p
!        real(kind=real64),dimension(2,2)    ::      a = reshape( (/1,0,0,1/),(/2,2/) )
!        real(kind=real64)                   ::      lambda = 64.0d0
!        character(len=5)                    ::      aaaa
!        integer                             ::      ii,iy
!        Mx = 2048
!        My = 2048
!        Nx = nint(Mx/lambda)
!        Ny = nint(My/lambda)
!        pbc = .false.
!        print *,"test old lib"
!        
!        allocate(p(0:Mx-1,0:My-1))
!        
!        !call PerlinNoise( Nx,Ny,p , pbc)
!        call perlinNoise( a,Mx,My, lambda , p )
!        
!            where (p>=0)
!                p = 1.0d0
!            else where 
!                p = 0.0d0
!            end where       
!        call writePng( "test.png",p ,normalise=.true. )     
!        do ii = 1,10  
!            call perlinNoise( a,Mx,My, lambda , p )
!            print *,"minmaxavg ",minval(p),maxval(p),sum(p)/(Mx*My)," > 0",count(p>=0)/(Mx*My*0.01),"%"
!                    
!            where (p>=0)
!                p = 1.0d0
!            else where 
!                p = -1.0d0
!            end where            
!            
!            open(unit=400,file="test.dat",action="write")
!            write(unit=400,fmt='(a,i8,f16.8)') "# Perlin noise ",Mx,lambda
!            write(unit=400,fmt='(2i8)') Mx,My            
!            do iy = 0,My-1
!                write(400,fmt='(4096f12.4)') p(:,iy)
!            end do 
!            close(unit=400)
!            !call writePng( "test.png",p ,normalise=.true. )
!            call system( "../2dFiltering/bin/contourAnalyse.exe -f test.dat -nLevels 1 -iso_min 0.0 -iso_max 0.0 | grep ""level"" -A 1 " ) 
!        end do                
!        
!        print *,""
!        print *,"done"
!        print *,""
!        
!    end program testLib_Perlin