
    module Lib_MaxLikelihoodFilter
!---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
        use Lib_ConjugateGradient
        use iso_fortran_env
        implicit none
        private
        
        public      ::      maxLikelihoodFilter
        
        
        
    contains
!---^^^^^^^^

        subroutine maxLikelihoodFilter( f,g, s,t )
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      compute a maximum likelihood filtering with range given by t
    !*      The probability of finding scene f given that the true solution is g
    !*      is  P = p1 p2
    !*      where   p1 = Exp[ - (f-g)^2/(2 s^2 ) ]
    !*              p2 = Exp[ - ( del^2(g) )^2/(2 t^2) ]
    !*      so the log likelihood is proportional to
    !*          lambda = -(t^2/2) (f-g)^2 - (s^2/2) ( del^2(g) )^2     
    !*      which we can maximise wrt g with
    !*          d lambda / dg = 0
    !*      which is a set of linear equations in f,g,s,t
    
            real(kind=real64),dimension(:,:),intent(in)     ::      f
            real(kind=real64),dimension(:,:),intent(out)    ::      g
            real(kind=real64),intent(in)                    ::      s,t
            
            real(kind=real64),dimension(:,:),allocatable    ::      AA              !   (1:13,Nx Ny)    -   filter kernel 
            integer,dimension(:,:),allocatable              ::      indx            !   (0:13,Nx Ny)    -   index of non-zero rows/cols    
            
            integer         ::      Nx,Ny
            real(kind=real64),dimension(:),allocatable      ::      ff,gg
            real(kind=real64)       ::       ww,xx,eps,s2ont2
            integer         ::      ii,jj,ix,iy,jx,jy,kx,ky , dd,mi
            logical         ::      ok 
            
            real(kind=real64),dimension(0:4),parameter      ::      kernel4 = (/ 20,-8,2,0,1 /)
            real(kind=real64),dimension(0:8),parameter      ::      kernel8 = (/ 36,0,-14,0,9,-4,0,0,4 /)/9.0d0
            integer,parameter                               ::      bandwidth_kernel4 = 13
            integer,parameter                               ::      bandwidth_kernel8 = 25
            
            logical             ::      useKernel8 = .false.
            integer                                         ::      bandwidth            
            real(kind=real64),dimension(:),allocatable      ::      kernel
            
            
            if (useKernel8) then
                allocate(kernel(0:size(kernel8)-1))
                kernel = kernel8
                bandwidth = bandwidth_kernel8
            else 
                allocate(kernel(0:size(kernel4)-1))
                kernel = kernel4
                bandwidth = bandwidth_kernel4
            end if
            
           
            
            Nx = size(f,dim=1)
            Ny = size(f,dim=2)
            allocate(AA(bandwidth,Nx*Ny))
            allocate(indx(0:bandwidth,Nx*Ny))
            allocate(ff(Nx*Ny))
            allocate(gg(Nx*Ny))
        
            
      
        !---    pack f into vector    
            do jj = 1,Ny
                ff( (jj-1)*Nx+1:jj*Nx ) = f(1:Nx,jj)
            end do
            gg = ff
            
        !---    compute laplacian squared kernel    
            
            !kernel(0) = 20*(s*s)/(t*t)      !   note I'm setting central element at end        !   0   0   1   0   0
            !kernel(1) = -8*(s*s)/(t*t)                                                         !   0   2  -8   2   0
            !kernel(2) = 2*(s*s)/(t*t)                                                          !   1  -8  20  -8   1
            !kernel(3) = 0                                                                      !   0   2  -8   2   0
            !kernel(4) = (s*s)/(t*t)                                                            !   0   0   1   0   0 
            !
            
            !kernel(0:4) = (/ -20,-8,2,0,1 /) * (s*s)/(t*t)
            
            !kernel(0:8) = ((/ 36,0,-14,0,9,-4,0,0,4 /)/9.0d0) * (s*s)/(t*t)
            
            kernel = kernel * (s*s)/(t*t)
            
            
        !---    compute the non-local filter kernel as a sparse matrix
            do iy = 1,Ny
                do ix = 1,Nx
                
                    ii = ix + (iy-1)*Nx     
                    mi = 0
                    
                    do ky = -2,2
                        jy = ky + iy
                        ok = ( jy*(Ny+1-jy)>0 )             !   inside region
                        do kx = -2,2
                            dd = kx*kx + ky*ky
                            if (dd >= size(kernel)) cycle               !   outside kernel
                            
                            jx = kx + ix
                            if (ok .and. ( jx*(Nx+1-jx)>0 )) then   !   inside region
                                jj = jx + (jy-1)*Nx   
                                if (ii==jj) cycle           !   note: am going to set central element at end...  
                                mi = mi + 1
                                xx = kernel(dd)
                                ww = ww + xx
                                AA(mi,ii) = xx
                                indx(mi,ii) = jj
                            end if
                        end do
                    end do
                    
                !---    normalise the kernel, add diagonal
                    ww = sum(AA(1:mi,ii))
                    mi = mi + 1
                    indx(mi,ii) = ii
                    AA(mi,ii) = 1.0d0 - ww
                    indx(0,ii) = mi
                    
                end do
            end do
            
        !---    now use conjugate gradients to find maximum 
            do ii = 1,3
                eps = 1.0d-8
                call conjgrad( AA,indx, gg,ff , eps )     
                if ( eps < 1.0d-8 ) exit       
            end do
                  
       
            g(1:Nx,1:Ny) = reshape(gg,(/Nx,Ny/))
            
            
            return
        end subroutine maxLikelihoodFilter     
            
            
        
    end module Lib_MaxlikelihoodFilter    
    
    
! !   gcc -c ${MYF90LIB}/Lib_Greyscale.c; gfortran -c -ffree-line-length-256 ${MYF90LIB}/Lib_Png.f90 ${MYF90LIB}/Lib_RandomSeed.f90 -lpng ; gfortran -c -ffree-line-length-256 -fopenmp Lib_ConjugateGradient.F90 -O2 -march=tigerlake
! !   gfortran -c -ffree-line-length-256 ${MYF90LIB}/Lib_RidlerCalvard.F90 Lib_MaxlikelihoodFilter.f90 -O2 -march=tigerlake; gfortran Lib_Greyscale.o Lib_Png.o Lib_RandomSeed.o Lib_ConjugateGradient.o Lib_RidlerCalvard.o Lib_MaxlikelihoodFilter.o -o testLib_MaxlikelihoodFilter.exe -fopenmp -lpng -llapack
! 
!     
!     program testLib_MaxlikelihoodFilter
! !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
!         use Lib_MaxlikelihoodFilter
!         use Lib_Png
!         use Lib_RandomSeed
!         use Lib_RidlerCalvard
!         use iso_fortran_env
!         implicit none
!         
!         integer             ::      Nx,Ny
!         
!         integer             ::      ii,jj,ix,iy
!         integer             ::      i0,j0
!         real(kind=real64)   ::      sigma,i2s2,dd,fbar,f2bar
!         
!         real(kind=real64),dimension(:,:),allocatable      ::      f,g
!         
!         real(kind=real64)   ::      ss, lambda ,bb,tt,ff,bstd,snr,fot
!         
!         character(len=256)  ::      dummy,filename
!         
!         print *,"testLib_MaxlikelihoodFilter.exe"
!         print *,"^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^"
!         print *,""
!         print *,"usage:"
!         print *,"   ./testLib_MaxlikelihoodFilter.exe filename lengthscale"
!          print *,""
!          call get_command_argument( 1,filename )
!          call get_command_argument( 2,dummy )
!          read(dummy,fmt=*) lambda
!          print *,"file = """//trim(filename)//""""
!          print *,"length = ",lambda
!          print *,""
!          
!          call readPng( filename,f )
!          Nx = size(f,dim=1)
!          Ny = size(f,dim=2)
!          allocate(g(Nx,Ny))
!          
!         
!         
! !        i0 = Nx/2
! !        j0 = Ny/2
! !         
! !        sigma = min(Nx,Ny)/8
! !        i2s2 = 1/(2*sigma*sigma)
! !                                        
! !        do jj = 1,Ny
! !            f(:,jj) = gaussianVariate( Nx )            
! !        end do 
! !        f = f * 0.05d0 + 0.1d0
! !           
! !        do iy = 1,Ny
! !            do ix = 1,Ny
! !            
! !                dd = (ix-i0)*(ix-i0) + (iy-j0)*(iy-j0)
! !                dd = dd*i2s2
! !                f(ix,iy) = f(ix,iy) + 0.8*exp( -dd )
! !                
! !            end do
! !        end do
! !                
! !        
! !        f = max( 0.0d0,min(1.0d0,f) )        
! !        call writePng("noisy.png",f)
!         
! !        fbar = 0.0d0
! !        f2bar = 0.0d0
! !        do iy = 1,Ny
! !            do ix = 1,Ny
! !                fbar = fbar + f(ix,iy)
! !                f2bar = f2bar + f(ix,iy)*f(ix,iy)
! !            end do               
! !        end do
! !        fbar = fbar/(Nx*Ny)
! !        f2bar = f2bar/(Nx*Ny)
! !       ss = sqrt( max(0.0d0,f2bar-fbar*fbar) ) 
! 
!         call findImageIntensityFeatures( f, bb,tt,ff,bstd,fbar,snr,fot )
! 
!         ss = bstd
!         
!         tt = fbar/lambda
!         
!         print *,"ss   = ",ss
!         print *,"fbar = ",fbar
!         print *,"tt   = ",tt
!         
!         
!         call maxLikelihoodFilter( f,g, ss,tt )
!         g = max( 0.0d0,min(1.0d0,g) )        
!         call writePng(trim(removePngExtension(filename))//".filt.png",g)
!         
!         print *,""
!         print *,"done"
!         print *,""
!     end program testLib_MaxlikelihoodFilter    
!         
        
        
        
        
        
        
        
        
        
        
        
        