
    module Lib_ComputePhaseFactor
!---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
!*      Given a set of atom positions r_j 
!*      and a set of g-vectors g_k
!*      compute the complex phase factors x_kj = exp(- i g_k.r_j )
!*
!*      To call computePhaseFactor
!*      note that we are computing x on an evenly spaced cubic grid, with grid spacing a.
!*      We therefore need to affine translate the position vectors of the atoms to 
!*          rt = R ( r - delta )/a
!*      where R is a rotation matrix so that the atom positions are covered by the grid 
!*      We will also ask for a section of grid to be computed, say ( lbound:ubound,lbound:ubound, 0:Nz-1 )
!*      so it is not necessary to pass all the atoms, only those with  lbound<=rt<=ubound
!*


!*      Daniel Mason
!*      (c) UKAEA September 2024
!*
        use Lib_SimpleProgressBar
        use Lib_ImagingSpace
        use iso_fortran_env
#ifdef MPI
        use mpi_f08
#endif
        implicit none
        private

        real(kind=real64),parameter                 ::      PI = 3.14159265390d0

        integer,private                 ::      rank = 0, nProcs = 1

        logical,public          ::      Lib_ComputePhaseFactor_DBG = .false.      
 
        public                  ::      Lib_ComputePhaseFactor_init_MPI

        public                  ::      computePhaseFactor
      





    contains
!---^^^^^^^^


        subroutine computePhaseFactor( mynAtoms, rt,g, img, x, gradx,rho )
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      compute a smoothed approximation x ~ exp( - i g.u )  
    !*      and at the same time a density estimation rho(r) and grad x
    !*      I compute a correctly weighted but unnormalised x~ for each atom j
    !*                     exp( -i g.Rj ) exp( -s (Rj-rk)² ) 
    !*          x~(rk) =   ---------------------------------
    !*                          sum_k exp( -s (Rj-rk)² ) 
    !*      where s = a²/(2 sig²) and the denominator shares atom j with unit weight over nearby nodes k
    !*      and the density function
    !*                         exp( -s (Rj-rk)² ) 
    !*          rho(rk) =    -----------------------
    !*                      sum_k exp( -s (Rj-rk)² )   
    !*      and the gradient function 
    !*                          2 s (Rj-rk) exp( -i g.Rj ) exp( -s (Rj-rk)² ) 
    !*          grad~(rk)   =   ----------------------------------------------
    !*                                     sum_k exp( -s (Rj-rk)² ) 
    !*      From these I then construct the phase factor
    !*          x       = x~/|x~|
    !*      and the phase factor gradient    
    !*          ∇x      = ∇( x~/|x~| ) 
    !*                  = grad~ /|x~|  - 1/|x~|² ∇|x~|
    !*                  = grad~ /|x~|  - x~/|x~|² (  x~ grad~* + x~* grad~ )/(2 |x~|)
    !*                  = ( grad~ - (x~ x~ grad~* /|x~|² )/(2 |x~|)
    !*                  = ( grad~ - x² grad~* )/(2 |x~|)
    !*

            integer,intent(in)                                          ::      mynAtoms
            real(kind=real64),dimension(1:,1:),intent(in)               ::      rt      !   (1:3,1:myNatoms) atom positions, scaled to rt = R (r-delta)/a. Typically in range  lbound(x):ubound(x)
            real(kind=real64),dimension(1:,1:),intent(in)               ::      g       !   (1:3,1:nGvec), reciprocal lattice vectors, with lengths in 1/A
            type(ImagingSpace),intent(in)                               ::      img 
            complex(kind=real64),dimension(:,:,:,:),pointer,intent(inout)       ::      x       !   (1:nGvec,lbx:ubx,lby:uby,lbz:ubz)
            complex(kind=real64),dimension(:,:,:,:,:),pointer,intent(inout)     ::      gradx   !   (3,1:nGvec,lbx:ubx,lby:uby,lbz:ubz)
            real(kind=real64),dimension(:,:,:),pointer,intent(inout)            ::      rho     !   (lbx:ubx,lby:uby,lbz:ubz)

            real(kind=real64),dimension(size(g,dim=2))          ::      gdotoff
            integer                                             ::      nGvec                  
            integer                                             ::      nn              !   number of points to search in kernel
            real(kind=real64),dimension(:,:),allocatable        ::      dist_kernel     !   (0:3,1:nn)
            complex(kind=real64),dimension(:,:),allocatable     ::      x_weight        !   (1:nGvec,1:nn)          weight to add to each nearby node.  
            complex(kind=real64),dimension(:,:,:),allocatable   ::      grad_x_weight   !   (3,1:nGvec,1:nn)        weight to add to each nearby node.  
            real(kind=real64),dimension(:),allocatable          ::      rho_weight      !   (1:nn)                  weight to add to each nearby node.  
            integer,dimension(:,:),allocatable                  ::      dx_kernel       !   (1:3,1:nn)      vector to neighbour nodes      
            real(kind=real64)                                   ::      aa,dd,a2on2s2,rmy1,rmy2,rmy3 ,d0 , weightsum , gdotr  
            complex(kind=real64),dimension(size(g,dim=2))       ::      eigdotr
            complex(kind=real64)                                ::      xx
            complex(kind=real64),dimension(3)                   ::      zz
            integer             ::      lbx,ubx,lby,uby,lbz,ubz         !   bounds of the region covered by the nodes x
            integer             ::      ix,iy,iz , jx,jy,jz             !   node numbers
            integer             ::      ii,jj,kk 

        !---    find the size of the problem                 
            nGvec = size(g,dim=2)

        !---    find the bounds of the set of nodes I am adding to
            lbx = lbound(x,dim=2) ; ubx = ubound(x,dim=2)
            lby = lbound(x,dim=3) ; uby = ubound(x,dim=3)
            lbz = lbound(x,dim=4) ; ubz = ubound(x,dim=4)



            x = complex( 1.0d-32,0.0d0 )        !   this trick of setting x to nearly zero means I don't have to test for exact zero to normalise, and means "empty" cells are given value x = 1+i0.
            gradx = 0.0d0
            rho = 0.0d0
            print *,"ComputePhaseFactor lbx,ubx,lby,uby,lbz,ubz ",lbx,ubx,lby,uby,lbz,ubz," nGvec ",nGvec


        !---    find the Gaussian kernel
            aa = geta(img)
            dd = aa/getsigma(img)
            a2on2s2 = dd*dd/2
            call computeKernels( geta(img),getsigma(img), nn, dx_kernel,dist_kernel )
            allocate(x_weight(nGvec,nn))                            !   to collect exp[ -i g.r - |r|^2/(2sig^2) ] on each node
            allocate(grad_x_weight(3,nGvec,nn))                     !   to collect grad exp[ -i g.r - |r|^2/(2sig^2) ] on each node
            
            allocate(rho_weight(nn))                                !   to collect exp[ - |r|^2/(2sig^2) ] on each node

        !---    compute g.offset, a constant phase factor for all atoms in the supercell. 
        !       strictly speaking shouldn't change much, unless something weird is done to break up the input files.
            do jj = 1,nGvec
                gdotoff(jj) = aa*dot_product( g(:,jj),getdelta(img) )
            end do


            do ii = 1,myNatoms

                if (rank==0) call progressBar( ii,mynAtoms )

            !---    compute exp( - i g.r ) for this atom
                do jj = 1,nGvec
                    dd = g(1,jj)*rt(1,ii) + g(2,jj)*rt(2,ii) + g(3,jj)*rt(3,ii)             !   g.(r - delta)/a)
                    gdotr = dd*aa + gdotoff(jj)        
                    eigdotr(jj) = complex( cos(gdotr),-sin(gdotr) )                    
                end do

            !---    find nearest node to atom i , remembering that nodes are at (1/2,1/2,1/2) positions
                ix = floor( rt(1,ii) ) ; rmy1 = rt(1,ii) - (ix+0.5d0) 
                iy = floor( rt(2,ii) ) ; rmy2 = rt(2,ii) - (iy+0.5d0)
                iz = floor( rt(3,ii) ) ; rmy3 = rt(3,ii) - (iz+0.5d0)
                
            !---    compute Gaussian weighting of exp( -i g.r ) on neighbour nodes
                
            !       we want to compute
            !           d^2 = ( r - (y+dy) )^2 (a/s)^2 /2 
            !           = ((a/s)^2/2)  *  ( (r-y)^2 - 2 r.dy + (dy)^2 )
            !       
            !       which we could do as q = (r - y) - dy, d = a2on2s2 * |q|^2 
            !       ... this is 9 flops.
            !       or we can do it as
            !           = (rmy1*rmy1 + rmy2*rmy2 + rmy3*rmy3)*a2on2s2        
            !               - rmy.dist_kernel
            !               + dist_kernel(0)
            !       ... which is 7.
            !       Yes. This is mission critical. So yes, 20% improvement is worth it.

 
                x_weight = 0.0d0
                rho_weight = 0.0d0
!                grad_x_weight = 0.0d0
                
                d0 = (rmy1*rmy1 + rmy2*rmy2 + rmy3*rmy3)*a2on2s2        
                 
                do kk = 1,nn
                    dd = d0 + dist_kernel(0,kk) + rmy1*dist_kernel(1,kk) + rmy2*dist_kernel(2,kk) + rmy3*dist_kernel(3,kk)   !   7 flops.
                    if (dd <= Lib_ImagingSpace_SIGMA_MULT*Lib_ImagingSpace_SIGMA_MULT) then
                        dd = exp( - dd )
                        rho_weight(kk) = dd                                                             !   exp( -s (Rj-rk)² ) 
                        do jj = 1,nGvec

                            x_weight(jj,kk) = dd * eigdotr(jj)                                          !   exp( -i g.Rj ) exp( -s (Rj-rk)² ) 

!                            grad_x_weight(1,jj,kk) = dd * eigdotr(jj) * rmy1*dist_kernel(1,kk)          !   2 s (Rj-rk) exp( -i g.Rj ) exp( -s (Rj-rk)² )  
!                            grad_x_weight(2,jj,kk) = dd * eigdotr(jj) * rmy2*dist_kernel(2,kk)          
!                            grad_x_weight(3,jj,kk) = dd * eigdotr(jj) * rmy3*dist_kernel(3,kk)  

                        end do
                    end if
                end do

            !---    normalise weighting so each atom actually has unit weight spread across each node, regardless of position wrt nodes.             
                weightsum = sum(rho_weight)   
                d0 = 1/max(1.0d-32,weightsum)                                                           !   sum_k exp( -s (Rj-rk)² ) 
                do kk = 1,nn
                    jx = ix + dx_kernel(1,kk)
                    jy = iy + dx_kernel(2,kk)
                    jz = iz + dx_kernel(3,kk)
                    if (inBounds(jx,jy,jz)) then
                        rho(jx,jy,jz) = rho(jx,jy,jz) + rho_weight(kk) * d0                              
                        x(1:nGvec,jx,jy,jz) = x(1:nGvec,jx,jy,jz) + x_weight(1:nGvec,kk) * d0
!                        gradx(1:3,1:nGvec,jx,jy,jz) = gradx(1:3,1:nGvec,jx,jy,jz) +  grad_x_weight(1:3,1:nGvec,kk) * d0
                    end if
                end do

            end do

        !---    normalise x = x~/|x~| so it is a phase factor everywhere.
            do jz = lbz,ubz
                do jy = lby,uby
                    do jx = lbx,ubx
                        do jj = 1,nGvec
                            xx = x(jj,jx,jy,jz)
                            dd = 1/abs(xx)    
                            xx = xx*dd                   
                            x(jj,jx,jy,jz) = xx
!                            zz = gradx(1:3,jj,jx,jy,jz) 
!                            gradx(1:3,jj,jx,jy,jz) = ( zz - xx*xx*conjg( zz ) )*dd/(2*aa)
                        end do
                    end do
                end do
            end do

        !---    compute gradx
            zz = complex(1.0d0,0.0d0)
            do jz = lbz,ubz
                do jy = lby,uby
                    do jx = lbx,ubx

                        if (jx==lbx) then
                            gradx(1,:,jx,jy,jz) = ( -3*x(:,jx,jy,jz) + 4*x(:,jx+1,jy,jz) - x(:,jx+2,jy,jz) )
                        else if (jx==ubx) then
                            gradx(1,:,jx,jy,jz) = ( x(:,jx-2,jy,jz) - 4*x(:,jx-1,jy,jz) + 3*x(:,jx,jy,jz) )
                        else
                            gradx(1,:,jx,jy,jz) = ( -x(:,jx-1,jy,jz) + x(:,jx+1,jy,jz) )
                        end if

                        if (jy==lby) then
                            gradx(2,:,jx,jy,jz) = ( -3*x(:,jx,jy,jz) + 4*x(:,jx,jy+1,jz) - x(:,jx,jy+2,jz) )
                        else if (jy==uby) then
                            gradx(2,:,jx,jy,jz) = ( x(:,jx,jy-2,jz) - 4*x(:,jx,jy-1,jz) + 3*x(:,jx,jy,jz) )
                        else
                            gradx(2,:,jx,jy,jz) = ( -x(:,jx,jy-1,jz) + x(:,jx,jy+1,jz) )
                        end if
    
                        if (jz==lbz) then
                            gradx(3,:,jx,jy,jz) = ( -zz + x(:,jx,jy,jz+1) )
                        else if (jz==ubz) then
                            gradx(3,:,jx,jy,jz) = ( -x(:,jx,jy,jz-1) +zz )
                        else
                            gradx(3,:,jx,jy,jz) = ( -x(:,jx,jy,jz-1) + x(:,jx,jy,jz+1) )
                        end if

                    end do
                end do
            end do
            gradx = gradx/(2*aa)



            return

        contains
    !---^^^^^^^^

            pure logical function inBounds(jx,jy,jz)
        !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
        !*      is the node (jx,jy,jz) in bounds of the array x(:,:,:,:)
                integer,intent(in)          ::      jx,jy,jz
                inBounds = (jx>=lbx).and.(jx<=ubx)
                inBounds = inBounds.and.(jy>=lby).and.(jy<=uby)
                inBounds = inBounds.and.(jz>=lbz).and.(jz<=ubz)
                return
            end function inBounds


        end subroutine computePhaseFactor





        subroutine computeKernels( a,sigma, n, dx_kernel,dist_kernel )
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      given the cell side and broadening a,sigma
    !*      compute the kernel of neighbour cells which could be reached by an atom in cell (0,0,0)
    !*
    !*      +-------+-------+-------+-------+-------+
    !*      |       |       |       |       |       |
    !*      |       |   x   |   x   |   x   |       |
    !*      |       |       |       |       |       |
    !*      +-------+-------+-------+-------+-------+           central cell is (0,0,0), marked 0
    !*      |       |       |       |       |       |           an atom might be at any of the points * in the central cell
    !*      |   x   |   x   |   x   |   x   |   x   |           find the neighbour cells x within range 3 sigma of any possible atom in (0,0,0)
    !*      |       |       |       |       |       |
    !*      +-------+-------+-------+-------+-------+
    !*      |       |       | * * * |       |       |
    !*      |   x   |   x   | * 0 * |   x   |   x   |
    !*      |       |       | * * * |       |       |
    !*      +-------+-------+-------+-------+-------+
    !*      |       |       |       |       |       |
    !*      |   x   |   x   |   x   |   x   |   x   |
    !*      |       |       |       |       |       |
    !*      +-------+-------+-------+-------+-------+
    !*      |       |       |       |       |       |
    !*      |       |   x   |   x   |   x   |       |
    !*      |       |       |       |       |       |
    !*      +-------+-------+-------+-------+-------+
    !*
   
            real(kind=real64),intent(in)            ::      a,sigma
            integer,intent(out)                     ::      n
            integer,dimension(:,:),allocatable,intent(out)              ::      dx_kernel
            real(kind=real64),dimension(:,:),allocatable,intent(out)    ::      dist_kernel



            real(kind=real64),dimension(3,27),parameter         ::      corner3d = reshape (  (/ -1,-1,-1 ,  0,-1,-1 ,  1,-1,-1,    &
                                                                                                 -1, 0,-1 ,  0, 0,-1 ,  1, 0,-1,    &    
                                                                                                 -1, 1,-1 ,  0, 1,-1 ,  1, 1,-1,    &    
                                                                                                 -1,-1, 0 ,  0,-1, 0 ,  1,-1, 0,    &
                                                                                                 -1, 0, 0 ,  0, 0, 0 ,  1, 0, 0,    &    
                                                                                                 -1, 1, 0 ,  0, 1, 0 ,  1, 1, 0,    &    
                                                                                                 -1,-1, 1 ,  0,-1, 1 ,  1,-1, 1,    &
                                                                                                 -1, 0, 1 ,  0, 0, 1 ,  1, 0, 1,    &    
                                                                                                 -1, 1, 1 ,  0, 1, 1 ,  1, 1, 1 /) , (/ 3,27 /)  ) * 0.5d0
        
            ! real(kind=real64),dimension(2,9),parameter          ::      corner2d = reshape (  (/ -1,-1 ,  0,-1 ,  1,-1,     &
            !                                                                                      -1, 0 ,  0, 0 ,  1, 0,     &    
            !                                                                                      -1, 1 ,  0, 1 ,  1, 1 /) , (/ 2,9 /)  ) * 0.5d0

            integer             ::      ix,iy,iz , jj
            integer             ::      nbuf
            real(kind=real64)   ::      xmy1,xmy2,xmy3,dd,a2on2s2
            logical             ::      ok


            nbuf = buffer(a,sigma)
            dd = a/sigma
            a2on2s2 = dd*dd/2

            n = 0                
            do iz = -nbuf,nbuf
                do iy = -nbuf,nbuf
                    do ix = -nbuf,nbuf
                        ok = .false.
                        do jj = 1,ubound(corner3d,dim=2)
                            xmy1 = ix + corner3d(1,jj)
                            xmy2 = iy + corner3d(2,jj)
                            xmy3 = iz + corner3d(3,jj)
                            dd = xmy1*xmy1 + xmy2*xmy2 + xmy3*xmy3
                            if (dd*a2on2s2 <= Lib_ImagingSpace_SIGMA_MULT*Lib_ImagingSpace_SIGMA_MULT*0.5000001d0) ok = .true.
                        end do
                        if (ok) n = n + 1
                    end do
                end do
            end do
             
            
            allocate(dist_kernel(0:3,n))                 !   term 0 is used for (d_i^2)/(2s^2)  
            allocate(dx_kernel(3,n))                     !   used to store positions of neighbouring nodes in order


            n = 0
            do iz = -nbuf,nbuf
                do iy = -nbuf,nbuf
                    do ix = -nbuf,nbuf
                        ok = .false.
                        do jj = 1,ubound(corner3d,dim=2)
                            xmy1 = ix + corner3d(1,jj)
                            xmy2 = iy + corner3d(2,jj)
                            xmy3 = iz + corner3d(3,jj)
                            dd = xmy1*xmy1 + xmy2*xmy2 + xmy3*xmy3
                            if (dd*a2on2s2 <= Lib_ImagingSpace_SIGMA_MULT*Lib_ImagingSpace_SIGMA_MULT*0.5000001d0) ok = .true.
                        end do
                        if (ok) then
                            n = n + 1
                            dd = ix*ix+iy*iy+iz*iz
                            dist_kernel(0:3,n) = (/ dd,-2.0d0*ix,-2.0d0*iy,-2.0d0*iz /)*a2on2s2
                            dx_kernel(1:3,n) = (/ix,iy,iz/)
                        end if                            
                    end do
                end do
            end do  
            

            if (Lib_ComputePhaseFactor_DBG) then
                print *,"Lib_ComputePhaseFactor::computeKernels info"
                print *,"a      ",a
                print *,"sigma  ",sigma
                print *,"nbuf   ",nbuf
                write(*,fmt='(a4,3a8,100a16)') "n","dx","",""," |dx|^2 (a/s)^2/2 "," (-2.dx) (a/s)^2/2 " 
                do jj = 1,n
                    write(*,fmt='(i4,3i8,100f16.8)') jj,dx_kernel(:,jj),dist_kernel(0:3,jj)
                end do
            end if
            
            return
        end subroutine computeKernels

            
        subroutine Lib_ComputePhaseFactor_init_MPI()
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
#ifdef MPI
            integer             ::      ierror
            logical             ::      ok
            call MPI_Initialized(ok,ierror)
            if (.not. ok) call MPI_INIT(ierror)
            call MPI_COMM_SIZE(MPI_COMM_WORLD, nProcs, ierror)
            call MPI_COMM_RANK(MPI_COMM_WORLD, rank, ierror)
            if (rank == 0) print *,"Lib_ComputePhaseFactor::Lib_ComputePhaseFactor_init_MPI info - initialised with ",nProcs," processes"  
#else
            print *,"Lib_ComputePhaseFactor::Lib_ComputePhaseFactor_init_MPI info - serial mode"  
#endif            
            return
        end subroutine Lib_ComputePhaseFactor_init_MPI

 

    end module Lib_ComputePhaseFactor     
