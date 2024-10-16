
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
        public                  ::      densityScale


        interface   computePhaseFactor
            module procedure    computePhaseFactor0
            module procedure    computePhaseFactor1
        end interface

        interface   twentySevenPointStencil
            module procedure    real_twentySevenPointStencil
            module procedure    complex_twentySevenPointStencil
        end interface

        
    contains
!---^^^^^^^^

        subroutine computePhaseFactor0( mynAtoms, rt,g, img, x,grad_x,rho )
            
            integer,intent(in)                                          ::      mynAtoms
            real(kind=real64),dimension(1:,1:),intent(in)               ::      rt      !   (1:3,1:myNatoms) atom positions, scaled to rt = R (r-delta)/a. Typically in range  lbound(x):ubound(x)
            real(kind=real64),dimension(1:,1:),intent(in)               ::      g       !   (1:3,1:nGvec), reciprocal lattice vectors, with lengths in 1/A
            type(ImagingSpace),intent(in)                               ::      img 
            
            complex(kind=real64),dimension(:,:,:,:),pointer,intent(out)         ::      x           !   (1:nGvec,lbx:ubx,lby:uby,lbz:ubz)
            complex(kind=real64),dimension(:,:,:,:,:),pointer,intent(out)       ::      grad_x      !   (3,1:nGvec,lbx:ubx,lby:uby,lbz:ubz)
            real(kind=real64),dimension(:,:,:),pointer,intent(inout)            ::      rho         !   (lbx:ubx,lby:uby,lbz:ubz)

            x = 1.0d0
            grad_x = 0.0d0
            rho = 1.0d0
            return
        end subroutine computePhaseFactor0


        subroutine computePhaseFactor1( mynAtoms, rt,g, img, grad_arg_x,rho,x )
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
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
    !*
    !*      The value required by the calculation is x* (k+g)/|k| ∇x
    !*      Now it could be that we are doing a precession or tomography calc, so there may be multiple (k+g)/|k|
    !*      but I still can store and reuse x* ∇x
    !*      Since x = exp( - ig.p )                 where p ~ u(r) is a multivalued displacement field
    !*          ∇x  = ∇( - ig.p ) x = i x ∇ arg(x) 
    !*        x* ∇x = i ∇ arg(x) 
    !*      So compute the complex field x, but return the _real_ field ∇ arg(x) = Im( x* ∇ x )

            integer,intent(in)                                          ::      mynAtoms
            real(kind=real64),dimension(1:,1:),intent(in)               ::      rt      !   (1:3,1:myNatoms) atom positions, scaled to rt = R (r-delta)/a. Typically in range  lbound(x):ubound(x)
            real(kind=real64),dimension(1:,1:),intent(in)               ::      g       !   (1:3,1:nGvec), reciprocal lattice vectors, with lengths in 1/A
            type(ImagingSpace),intent(in)                               ::      img 
            
            real(kind=real64),dimension(:,:,:,:,:),pointer,intent(inout)        ::      grad_arg_x   !   (3,1:nGvec,lbx:ubx,lby:uby,lbz:ubz)
            real(kind=real64),dimension(:,:,:),pointer,intent(inout)            ::      rho     !   (lbx:ubx,lby:uby,lbz:ubz)
            complex(kind=real64),dimension(:,:,:,:),pointer,intent(inout),optional   ::      x


            complex(kind=real64),dimension(:,:,:,:),pointer     ::      xx       !   (1:nGvec,lbx:ubx,lby:uby,lbz:ubz)            
            real(kind=real64),dimension(size(g,dim=2))          ::      gdotoff
            integer                                             ::      nGvec                  
            integer                                             ::      nn              !   number of points to search in kernel
            real(kind=real64),dimension(:,:),allocatable        ::      dist_kernel     !   (0:3,1:nn)
            complex(kind=real64),dimension(:,:),allocatable     ::      x_weight        !   (1:nGvec,1:nn)          weight to add to each nearby node.  
            real(kind=real64),dimension(:),allocatable          ::      rho_weight      !   (1:nn)                  weight to add to each nearby node.  
            integer,dimension(:,:),allocatable                  ::      dx_kernel       !   (1:3,1:nn)      vector to neighbour nodes      
            real(kind=real64)                                   ::      aa,dd,a2on2s2,rmy1,rmy2,rmy3 ,d0 , weightsum , gdotr  
            complex(kind=real64),dimension(size(g,dim=2))       ::      eigdotr
            complex(kind=real64),dimension(3)                   ::      grad_x
            complex(kind=real64),dimension(size(g,dim=2),-1:1,-1:1,-1:1)    ::      local_x
            integer             ::      lbx,ubx,lby,uby,lbz,ubz         !   bounds of the region covered by the nodes x
            integer             ::      ix,iy,iz , jx,jy,jz             !   node numbers
            integer             ::      ii,jj,kk 
 

        !---    find the size of the problem                 
            nGvec = size(g,dim=2)

        !---    find the bounds of the set of nodes I am adding to
            lbx = lbound(rho,dim=1) ; ubx = ubound(rho,dim=1)
            lby = lbound(rho,dim=2) ; uby = ubound(rho,dim=2)
            lbz = lbound(rho,dim=3) ; ubz = ubound(rho,dim=3)

            allocate(xx(1:nGvec,lbx-1:ubx+1,lby-1:uby+1,lbz-1:ubz+1))
            xx = 0.0d0
            grad_arg_x = 0.0d0
            rho = 0.0d0
            print *,"ComputePhaseFactor lbx,ubx,lby,uby,lbz,ubz ",lbx,ubx,lby,uby,lbz,ubz," nGvec ",nGvec


        !---    find the Gaussian kernel
            aa = geta(img)
            dd = aa/getsigma(img)
            a2on2s2 = dd*dd/2
            call computeKernels( geta(img),getsigma(img), nn, dx_kernel,dist_kernel )
            allocate(x_weight(nGvec,nn))                            !   to collect exp[ -i g.r - |r|^2/(2sig^2) ] on each node            
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
                
                d0 = (rmy1*rmy1 + rmy2*rmy2 + rmy3*rmy3)*a2on2s2        
                 
                do kk = 1,nn
                    dd = d0 + dist_kernel(0,kk) + rmy1*dist_kernel(1,kk) + rmy2*dist_kernel(2,kk) + rmy3*dist_kernel(3,kk)   !   7 flops.
                    if (dd <= Lib_ImagingSpace_SIGMA_MULT*Lib_ImagingSpace_SIGMA_MULT) then
                        dd = exp( - dd )
                        rho_weight(kk) = dd                                                             !   exp( -s (Rj-rk)² ) 
                        do jj = 1,nGvec

                            x_weight(jj,kk) = dd * eigdotr(jj)                                          !   exp( -i g.Rj ) exp( -s (Rj-rk)² ) 
 
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
                        xx(1:nGvec,jx,jy,jz) = xx(1:nGvec,jx,jy,jz) + x_weight(1:nGvec,kk) * d0
                    end if
                end do

            end do
 

         !---    add the padding around the x stored region by extrapolating the derivative           
             do jy = lby,uby
                 do jx = lbx,ubx
                     xx(:,jx,jy,lbz-1) =  3*xx(:,jx,jy,lbz) - 3*xx(:,jx,jy,lbz+1) + xx(:,jx,jy,lbz+2)
                 end do
             end do
             do jy = lby,uby
                 do jx = lbx,ubx
                     xx(:,jx,jy,ubz+1) =  xx(:,jx,jy,ubz-2) - 3*xx(:,jx,jy,ubz-1) + 3*xx(:,jx,jy,ubz) 
                 end do
             end do

             do jz = lbz-1,ubz+1
                 do jx = lbx,ubx
                     xx(:,jx,lby-1,jz) =  3*xx(:,jx,lby,jz) - 3*xx(:,jx,lby+1,jz) + xx(:,jx,lby+2,jz)
                 end do
             end do
             do jz = lbz-1,ubz+1
                 do jx = lbx,ubx
                     xx(:,jx,uby+1,jz) =  xx(:,jx,uby-2,jz) - 3*xx(:,jx,uby-1,jz) + 3*xx(:,jx,uby,jz) 
                 end do
             end do

             do jz = lbz-1,ubz+1
                 do jy = lby-1,uby+1
                     xx(:,lbx-1,jy,jz) =  3*xx(:,lbx,jy,jz) - 3*xx(:,lbx+1,jy,jz) + xx(:,lbx+2,jy,jz)
                 end do
             end do
             do jz = lbz-1,ubz+1
                 do jy = lby-1,uby+1
                     xx(:,ubx+1,jy,jz) =  xx(:,ubx-2,jy,jz) - 3*xx(:,ubx-1,jy,jz) + 3*xx(:,ubx,jy,jz) 
                 end do
             end do



        !---    compute x* grad x everywhere      
            do jz = lbz,ubz
                if (rank==0) call progressBar( jz+1-lbz,ubz+1-lbz )
                do jy = lby,uby

                    local_x(:,0:1,-1:1,-1:1) = xx(:,lbx-1:lbx,jy-1:jy+1,jz-1:jz+1) 
                    do jx = lbx,ubx
                        local_x(:,-1,:,:) = local_x(:,0,:,:)
                        local_x(:,0,:,:) = local_x(:,1,:,:)
                        local_x(:,1,:,:) = xx(:,jx+1,jy-1:jy+1,jz-1:jz+1)

                        do jj = 1,nGvec
                            call twentySevenPointStencil( f = local_x(jj,:,:,:) , df = grad_x )
                            grad_arg_x(1:3,jj,jx,jy,jz) = aimag( conjg( xx(jj,jx,jy,jz) )*grad_x(1:3) )
                        end do
                        
                    end do
                end do
            end do

            grad_arg_x = grad_arg_x / aa        !   normalise by length

            if (present(x)) then
                x => xx
            else
                deallocate(xx)
            end if
 
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


        end subroutine computePhaseFactor1


        subroutine real_twentySevenPointStencil( f, f0,df )
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ 
    !*      approximate f = f0 + dfx x + dfy y + dfz z
    !*      using a 27 point stencil

            real(kind=real64),dimension(-1:1,-1:1,-1:1),intent(in)              ::          f
            real(kind=real64),intent(out),optional                              ::          f0
            real(kind=real64),dimension(3),intent(out),optional                 ::          df
            

        !---    weights of points in the stencil. w(r) = exp( - 2 r^2  )  
            real(kind=real64),parameter             ::      W0 = 1.0d0          
            real(kind=real64),parameter             ::      W1 = exp( - 2.0d0 )
            real(kind=real64),parameter             ::      W2 = exp( - 4.0d0 )
            real(kind=real64),parameter             ::      W3 = exp( - 6.0d0 )

        !---    stencils            
            real(kind=real64),dimension(-1:1,-1:1,-1:1),parameter       ::      STENCIL_DX = reshape( (/        &
                                                                                             -W3,.0d0,  W3,     &               !   -1,-1,-1 : 1,-1,-1
                                                                                             -W2,.0d0,  W2,     &               !   -1, 0,-1 : 1, 0,-1
                                                                                             -W3,.0d0,  W3,     &               !   -1, 1,-1 : 1, 1,-1
                                                                                             -W2,.0d0,  W2,     &               !   -1,-1, 0 : 1,-1, 0
                                                                                             -W1,.0d0,  W1,     &               !   -1, 0, 0 : 1, 0, 0
                                                                                             -W2,.0d0,  W2,     &               !   -1, 1, 0 : 1, 1, 0
                                                                                             -W3,.0d0,  W3,     &               !   -1,-1, 1 : 1,-1, 1    
                                                                                             -W2,.0d0,  W2,     &               !   -1, 0, 1 : 1, 0, 1
                                                                                             -W3,.0d0,  W3      &               !   -1, 1, 1 : 1, 1, 1
                                                                                             /)  , (/3,3,3/) ) / ( 2*W1 + 8*W2 + 8*W3 )

            real(kind=real64),dimension(-1:1,-1:1,-1:1),parameter       ::      STENCIL_DY = reshape( (/        &
                                                                                             -W3, -W2, -W3,     &               !   -1,-1,-1 : 1,-1,-1
                                                                                            .0d0,.0d0,.0d0,     &               !   -1, 0,-1 : 1, 0,-1
                                                                                              W3,  W2,  W3,     &               !   -1, 1,-1 : 1, 1,-1
                                                                                             -W2, -W1 ,-W2,     &               !   -1,-1, 0 : 1,-1, 0
                                                                                            .0d0,.0d0,.0d0,     &               !   -1, 0, 0 : 1, 0, 0
                                                                                              W2,  W1,  W2,     &               !   -1, 1, 0 : 1, 1, 0
                                                                                             -W3, -W2, -W3,     &               !   -1,-1, 1 : 1,-1, 1    
                                                                                            .0d0,.0d0,.0d0,     &               !   -1, 0, 1 : 1, 0, 1
                                                                                              W3,  W2,  W3      &               !   -1, 1, 1 : 1, 1, 1
                                                                                             /)  , (/3,3,3/) ) / ( 2*W1 + 8*W2 + 8*W3 )

            real(kind=real64),dimension(-1:1,-1:1,-1:1),parameter       ::      STENCIL_DZ = reshape( (/        &
                                                                                             -W3, -W2, -W3,     &               !   -1,-1,-1 : 1,-1,-1
                                                                                             -W2, -W1, -W2,     &               !   -1, 0,-1 : 1, 0,-1
                                                                                             -W3, -W2, -W3,     &               !   -1, 1,-1 : 1, 1,-1
                                                                                            .0d0,.0d0,.0d0,     &               !   -1,-1, 0 : 1,-1, 0
                                                                                            .0d0,.0d0,.0d0,     &               !   -1, 0, 0 : 1, 0, 0
                                                                                            .0d0,.0d0,.0d0,     &               !   -1, 1, 0 : 1, 1, 0
                                                                                              W3,  W2,  W3,     &               !   -1,-1, 1 : 1,-1, 1    
                                                                                              W2,  W1,  W2,     &               !   -1, 0, 1 : 1, 0, 1
                                                                                              W3,  W2,  W3      &               !   -1, 1, 1 : 1, 1, 1
                                                                                             /)  , (/3,3,3/) ) / ( 2*W1 + 8*W2 + 8*W3 )

                                               
            real(kind=real64),dimension(-1:1,-1:1,-1:1),parameter       ::      STENCIL_0  = reshape( (/        &
                                                                                              W3,  W2,  W3,     &               !   -1,-1,-1 : 1,-1,-1
                                                                                              W2,  W1,  W2,     &               !   -1, 0,-1 : 1, 0,-1
                                                                                              W3,  W2,  W3,     &               !   -1, 1,-1 : 1, 1,-1
                                                                                              W2,  W1,  W2,     &               !   -1,-1, 0 : 1,-1, 0
                                                                                              W1,  W0,  W1,     &               !   -1, 0, 0 : 1, 0, 0
                                                                                              W2,  W1,  W2,     &               !   -1, 1, 0 : 1, 1, 0
                                                                                              W3,  W2,  W3,     &               !   -1,-1, 1 : 1,-1, 1    
                                                                                              W2,  W1,  W2,     &               !   -1, 0, 1 : 1, 0, 1
                                                                                              W3,  W2,  W3      &               !   -1, 1, 1 : 1, 1, 1
                                                                                             /)  , (/3,3,3/) ) / ( W0 + 6*W1 + 12*W2 + 8*W3 )

            if (present(f0)) then
                f0 = sum( STENCIL_0*f )
            end if

            if (present(df)) then
                df(1) = sum( STENCIL_DX*f )
                df(2) = sum( STENCIL_DY*f )
                df(3) = sum( STENCIL_DZ*f )
            end if

            return
        end subroutine real_twentySevenPointStencil
 



        subroutine complex_twentySevenPointStencil( f, f0,df )
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ 
    !*      approximate f = f0 + dfx x + dfy y + dfz z
    !*      using a 27 point stencil

            complex(kind=real64),dimension(-1:1,-1:1,-1:1),intent(in)              ::          f
            complex(kind=real64),intent(out),optional                              ::          f0
            complex(kind=real64),dimension(3),intent(out),optional                 ::          df
            

        !---    weights of points in the stencil. w(r) = exp( - 2 r^2  )  
            real(kind=real64),parameter             ::      W0 = 1.0d0          
            real(kind=real64),parameter             ::      W1 = exp( - 2.0d0 )
            real(kind=real64),parameter             ::      W2 = exp( - 4.0d0 )
            real(kind=real64),parameter             ::      W3 = exp( - 6.0d0 )

        !---    stencils            
            real(kind=real64),dimension(-1:1,-1:1,-1:1),parameter       ::      STENCIL_DX = reshape( (/        &
                                                                                             -W3,.0d0,  W3,     &               !   -1,-1,-1 : 1,-1,-1
                                                                                             -W2,.0d0,  W2,     &               !   -1, 0,-1 : 1, 0,-1
                                                                                             -W3,.0d0,  W3,     &               !   -1, 1,-1 : 1, 1,-1
                                                                                             -W2,.0d0,  W2,     &               !   -1,-1, 0 : 1,-1, 0
                                                                                             -W1,.0d0,  W1,     &               !   -1, 0, 0 : 1, 0, 0
                                                                                             -W2,.0d0,  W2,     &               !   -1, 1, 0 : 1, 1, 0
                                                                                             -W3,.0d0,  W3,     &               !   -1,-1, 1 : 1,-1, 1    
                                                                                             -W2,.0d0,  W2,     &               !   -1, 0, 1 : 1, 0, 1
                                                                                             -W3,.0d0,  W3      &               !   -1, 1, 1 : 1, 1, 1
                                                                                             /)  , (/3,3,3/) ) / ( 2*W1 + 8*W2 + 8*W3 )

            real(kind=real64),dimension(-1:1,-1:1,-1:1),parameter       ::      STENCIL_DY = reshape( (/        &
                                                                                             -W3, -W2, -W3,     &               !   -1,-1,-1 : 1,-1,-1
                                                                                            .0d0,.0d0,.0d0,     &               !   -1, 0,-1 : 1, 0,-1
                                                                                              W3,  W2,  W3,     &               !   -1, 1,-1 : 1, 1,-1
                                                                                             -W2, -W1 ,-W2,     &               !   -1,-1, 0 : 1,-1, 0
                                                                                            .0d0,.0d0,.0d0,     &               !   -1, 0, 0 : 1, 0, 0
                                                                                              W2,  W1,  W2,     &               !   -1, 1, 0 : 1, 1, 0
                                                                                             -W3, -W2, -W3,     &               !   -1,-1, 1 : 1,-1, 1    
                                                                                            .0d0,.0d0,.0d0,     &               !   -1, 0, 1 : 1, 0, 1
                                                                                              W3,  W2,  W3      &               !   -1, 1, 1 : 1, 1, 1
                                                                                             /)  , (/3,3,3/) ) / ( 2*W1 + 8*W2 + 8*W3 )

            real(kind=real64),dimension(-1:1,-1:1,-1:1),parameter       ::      STENCIL_DZ = reshape( (/        &
                                                                                             -W3, -W2, -W3,     &               !   -1,-1,-1 : 1,-1,-1
                                                                                             -W2, -W1, -W2,     &               !   -1, 0,-1 : 1, 0,-1
                                                                                             -W3, -W2, -W3,     &               !   -1, 1,-1 : 1, 1,-1
                                                                                            .0d0,.0d0,.0d0,     &               !   -1,-1, 0 : 1,-1, 0
                                                                                            .0d0,.0d0,.0d0,     &               !   -1, 0, 0 : 1, 0, 0
                                                                                            .0d0,.0d0,.0d0,     &               !   -1, 1, 0 : 1, 1, 0
                                                                                              W3,  W2,  W3,     &               !   -1,-1, 1 : 1,-1, 1    
                                                                                              W2,  W1,  W2,     &               !   -1, 0, 1 : 1, 0, 1
                                                                                              W3,  W2,  W3      &               !   -1, 1, 1 : 1, 1, 1
                                                                                             /)  , (/3,3,3/) ) / ( 2*W1 + 8*W2 + 8*W3 )

                                               
            real(kind=real64),dimension(-1:1,-1:1,-1:1),parameter       ::      STENCIL_0  = reshape( (/        &
                                                                                              W3,  W2,  W3,     &               !   -1,-1,-1 : 1,-1,-1
                                                                                              W2,  W1,  W2,     &               !   -1, 0,-1 : 1, 0,-1
                                                                                              W3,  W2,  W3,     &               !   -1, 1,-1 : 1, 1,-1
                                                                                              W2,  W1,  W2,     &               !   -1,-1, 0 : 1,-1, 0
                                                                                              W1,  W0,  W1,     &               !   -1, 0, 0 : 1, 0, 0
                                                                                              W2,  W1,  W2,     &               !   -1, 1, 0 : 1, 1, 0
                                                                                              W3,  W2,  W3,     &               !   -1,-1, 1 : 1,-1, 1    
                                                                                              W2,  W1,  W2,     &               !   -1, 0, 1 : 1, 0, 1
                                                                                              W3,  W2,  W3      &               !   -1, 1, 1 : 1, 1, 1
                                                                                             /)  , (/3,3,3/) ) / ( W0 + 6*W1 + 12*W2 + 8*W3 )

            if (present(f0)) then
                f0 = sum( STENCIL_0*f )
            end if

            if (present(df)) then
                df(1) = sum( STENCIL_DX*f )
                df(2) = sum( STENCIL_DY*f )
                df(3) = sum( STENCIL_DZ*f )
            end if

            return
        end subroutine complex_twentySevenPointStencil
 

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

            
        elemental real(kind=real64) function smoothstep( x )
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      smoothstep function returns 0 for x<=0, 1 for x>=1, and smooth polynomial between
            real(kind=real64),intent(in)        ::      x
            if (x<=0) then
                smoothstep = 0.0d0
            else if (x>=1) then
                smoothstep = 1.0d0
            else
                smoothstep = x*x*x*(10-15*x+6*x*x)
            end if
        end function smoothstep

        elemental real(kind=real64) function densityScale( rho,omega )
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      smooth density function, returns 0 for rho*omega <= 0.5, 1 for rho*omega >= 1
            real(kind=real64),intent(in)            ::      rho
            real(kind=real64),intent(in)            ::      omega
            real(kind=real64)           ::          xx
            !xx = 2*rho*omega - 1
            xx = rho*omega
            densityScale =  smoothstep( xx )
            return
        end function densityScale 

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
