

    module Lib_ImagingSpace
!---^^^^^^^^^^^^^^^^^^^^^^^
!*      A short module defining a derived type containing a definition of the imaging space
!*
!*      The imaging space is a set of nodes 0:Nx-1,0:Ny-1,0:Nz-1
!*      with cubic spacing, cell side a
!*      and smearing sigma
!*  
!*      atoms will be given an affine transform
!*          rt = (r - delta)/a
!*      and then "placed" at a grid node int(rt)
!*      In order to get correct imaging at the grid node 0, we need atoms to be at positions 0.5-nBuf 
!*  
!*      For send/recv operations, the ordering of the neighbours is
!*
!*        8   1   2             y
!*          +---+               ^
!*        7 |   | 3             |
!*          +---+               +---->x
!*        6   5   4      

!*      Daniel Mason
!*      (c) UKAEA September 2024
!*
#ifdef MPI
        use mpi_f08
#endif
        use Lib_FactoriseParallel
        use iso_fortran_env
        implicit none
        private

        public      ::      Lib_ImagingSpace_init_MPI

        public      ::      ImagingSpace_ctor
        public      ::      report
        public      ::      delete
        
        public      ::      buffer
        public      ::      inMyCell
        
        

        public      ::      geta
        public      ::      getsigma
        public      ::      getdelta
        public      ::      getnBuf
        public      ::      getNx,getNy,getNz
        public      ::      getMx,getMy        
        public      ::      getBounds
        public      ::      sendrecv
        public      ::      whoseBlock
        public      ::      blockExtent

        integer,private                     ::      rank = 0
        integer,private                     ::      nProcs = 1         
        integer,public                      ::      Lib_ImagingSpace_SIGMA_MULT = 3

        integer,private,parameter           ::      N  = 1
        integer,private,parameter           ::      NE = 2
        integer,private,parameter           ::      E  = 3
        integer,private,parameter           ::      SE = 4
        integer,private,parameter           ::      S  = 5
        integer,private,parameter           ::      SW = 6
        integer,private,parameter           ::      W  = 7
        integer,private,parameter           ::      NW = 8
        


        type,public     ::      ImagingSpace
            private
            real(kind=real64)               ::      a                   !   grid side length
            real(kind=real64)               ::      sigma               !   smoothing length scale
            real(kind=real64),dimension(3)  ::      delta
            integer                         ::      nBuf                !   number of buffing cells required to ensure smoothness at boundaries
            integer                         ::      Nx,Ny,Nz            !   total number of image space cells, divided across all processors
            integer                         ::      Mx,My               !   number of image space blocks, one each per processor. Note Mz = 1.
            integer                         ::      ix,iy               !   image space block I am responsible for (0:Mx-1)
            integer                         ::      lbx,ubx,lby,uby     !   the region of image space cells I am responsible for 
            integer,dimension(:,:),pointer  ::      p                   !   (0:Mx-1,0:My-1) the processor who is responsible for each block
            integer,dimension(8)            ::      neigh               !   the processor that hold my neighbour

        end type 


        interface   ImagingSpace_ctor
            module procedure        ImagingSpace_null
            module procedure        ImagingSpace_ctor0
            module procedure        ImagingSpace_ctor1
        end interface
        
        interface   report
            module procedure        report0
        end interface
        
        interface   delete
            module procedure        delete0
        end interface

        interface   geta 
            module procedure        geta0
        end interface
 
        interface   getsigma 
            module procedure        getsigma0
        end interface
 
        interface   getdelta
            module procedure        getdelta0
        end interface
 
        interface   getnBuf
            module procedure        getnBuf0
        end interface

        interface   getNx
            module procedure        getNx0
        end interface

        interface   getNy
            module procedure        getNy0
        end interface

        interface   getNz
            module procedure        getNz0
        end interface

        interface   getMx
            module procedure        getMx0
        end interface

        interface   getMy
            module procedure        getMy0
        end interface

 
        interface   inMyCell
            module procedure        inMyCell0
            module procedure        inMyCell1
        end interface
 

        interface   getBounds
            module procedure        getBounds0
            module procedure        getBounds1
        end interface
 

    contains
!---^^^^^^^^


!******************************************************************************
!
!           standard functions
!
!******************************************************************************

    

        function ImagingSpace_null() result(this)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^    
            type(ImagingSpace)                     ::      this
            this%a = 1.0d0
            this%sigma = 0.5d0
            this%delta = 0.0d0
            this%nBuf = buffer(this%a,this%sigma)
            this%Nx = 1
            this%Ny = 1
            this%Nz = 1            
            this%Mx = 1
            this%My = 1
            this%ix = 0
            this%iy = 0
            nullify(this%p)
            this%neigh = rank
            this%Mx = 1
            this%My = 1
            !call setMyBounds(this)
            return
        end function ImagingSpace_null
        

        function ImagingSpace_ctor0(a,sigma,delta,Nx,Ny,Nz) result(this)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^    
            type(ImagingSpace)                          ::      this
            real(kind=real64),intent(in)                ::      a
            real(kind=real64),intent(in)                ::      sigma
            real(kind=real64),dimension(3),intent(in)   ::      delta
            integer,intent(in)                          ::      Nx,Ny,Nz
            this = ImagingSpace_null()
            this%a = a
            this%sigma = sigma
            this%delta = delta
            this%nBuf = buffer(this%a,this%sigma)
            this%Nx = Nx
            this%Ny = Ny
            this%Nz = Nz
            call FactoriseParallel( this%Nx,this%Ny,nProcs,this%Mx,this%My )            
            allocate(this%p(0:this%Mx-1,0:this%My-1))

            call setMyBounds(this)
            return
        end function ImagingSpace_ctor0

        function ImagingSpace_ctor1(a,delta,Nx,Ny,Nz) result(this)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ 
            type(ImagingSpace)                          ::      this
            real(kind=real64),intent(in)                ::      a
            real(kind=real64),dimension(3),intent(in)   ::      delta
            integer,intent(in)                          ::      Nx,Ny,Nz
            this = ImagingSpace_ctor0(a,a/2,delta,Nx,Ny,Nz)
            return
        end function ImagingSpace_ctor1


        subroutine delete0(this)
    !---^^^^^^^^^^^^^^^^^^^^^^^^
            type(ImagingSpace),intent(inout)       ::      this
            if (associated(this%p)) deallocate(this%p)
            this = ImagingSpace_null()
            return
        end subroutine delete0

        
        subroutine report0(this,u,o)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^
            type(ImagingSpace),intent(in)          ::      this
            integer,intent(in),optional         ::      u,o
            integer     ::      uu,oo
            integer     ::      ix,iy,Cx,Cy,pp
            uu = 6 ; if (present(u)) uu = u
            oo = 0 ; if (present(o)) oo = o
            write(unit=uu,fmt='(2(a,f16.8),5(a,i8))') repeat(" ",oo)//"ImagingSpace [a,sigma= ",this%a,",",this%sigma," , Nx,Ny,Nz= ",this%Nx,",",this%Ny,",",this%Nz,", nBuf=",this%nBuf," ]"
            !write(unit=uu,fmt='(6(a,i8))') repeat(" ",oo+4)//"my block ",this%lbx,":",this%ubx," , ",this%lby,":",this%uby," , ",0,":",this%Nz-1

            Cx = ceiling(real(this%Nx)/this%Mx)
            Cy = ceiling(real(this%Ny)/this%My)

            do pp = 0,nProcs-1
                iy = pp/this%Mx
                ix = pp - this%Mx*iy
                write (*,fmt='(10(a,i4))') repeat(" ",oo+4)//"proc ",pp,                        &
                                " x ",ix*Cx,":",min( this%Nx-1,ix*Cx + Cx - 1 ),                &
                                " y ",iy*Cy,":",min( this%Ny-1,iy*Cy + Cy - 1 )                  
            end do            
            return
        end subroutine report0




!******************************************************************************
!
!           member functions
!
!******************************************************************************

        
        pure logical function inMyCell0(this,xt,buffered,border)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      given the reduced space position
    !*          xt = R (x-delta)/a
    !*      return true if the atom is in a cell I am responsible for.
    !*      optionally return true if within the buffered region
            type(ImagingSpace),intent(in)               ::      this
            real(kind=real64),dimension(:),intent(in)   ::      xt
            logical,intent(in)                          ::      buffered
            integer,intent(in)                          ::      border      !   required for non-columnar apprxo
            integer         ::          ix,iy,iz
            integer         ::          pad

            pad = this%nBuf          !   x-y directions need buffer + 1 extra cells + optionally x-y needs extra border for non-columnar approx
            if (size(xt)==2) then
                !   2d test - check xt(1:2) is in the x-y bounds, but do not test for z
                
                ix = floor( xt(1) )            
                iy = floor( xt(2) )

                if (buffered) then
                    inMyCell0 = (ix>=this%lbx-pad-border) .and. (ix<=this%ubx+pad+border) 
                    inMyCell0 = inMyCell0 .and. (iy>=this%lby-pad-border) .and. (iy<=this%uby+pad+border) 
                else
                    inMyCell0 = (ix>=this%lbx) .and. (ix<=this%ubx) 
                    inMyCell0 = inMyCell0 .and. (iy>=this%lby) .and. (iy<=this%uby) 
                end if

            else
                !   3d test - check that in z bounds too

                ix = floor( xt(1) )            
                iy = floor( xt(2) )
                iz = floor( xt(3) )

                if (buffered) then
                    inMyCell0 = (ix>=this%lbx-pad-border) .and. (ix<=this%ubx+pad+border) 
                    inMyCell0 = inMyCell0 .and. (iy>=this%lby-pad-border) .and. (iy<=this%uby+pad+border) 
                    inMyCell0 = inMyCell0 .and. (iz>=-this%nBuf) .and. (iz<this%Nz+this%nBuf) 
                else
                    inMyCell0 = (ix>=this%lbx) .and. (ix<=this%ubx) 
                    inMyCell0 = inMyCell0 .and. (iy>=this%lby) .and. (iy<=this%uby) 
                    inMyCell0 = inMyCell0 .and. (iz>=0) .and. (iz<this%Nz) 
                end if
            end if

            return
        end function inMyCell0

        pure logical function inMyCell1(this,xt,buffered,border)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      given the reduced space position
    !*          xt = R (x-delta)/a
    !*      return true if the atom is in a cell I am responsible for.
    !*      optionally return true if within the buffered region
    
            type(ImagingSpace),intent(in)               ::      this
            real(kind=real32),dimension(:),intent(in)   ::      xt
            logical,intent(in)                          ::      buffered
            integer,intent(in)                          ::      border      !   required for non-columnar apprxo
            inMyCell1 = inMyCell0(this,real(xt,kind=real64),buffered,border)

            return
        end function inMyCell1        

        subroutine setMyBounds(this)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      given my rank and number of processes, determine the block of cells I am responsible for 
    !*      note that this will _not_ include any buffer
            type(ImagingSpace),intent(inout)            ::      this

            integer         ::          Cx,Cy                       !   size of each block
            integer         ::          pp,ix,iy                    !   processor number, position in grid
            
            Cx = ceiling(real(this%Nx)/this%Mx)
            Cy = ceiling(real(this%Ny)/this%My)

            do pp = 0,nProcs-1
                iy = pp/this%Mx
                ix = pp - this%Mx*iy
                this%p(ix,iy) = pp
                if (pp == rank) then
                    this%lbx = ix*Cx
                    this%ubx = min( this%Nx-1,this%lbx + Cx - 1 )
                    this%lby = iy*Cy
                    this%uby = min( this%Ny-1,this%lby + Cy - 1 )  
                    this%ix = ix
                    this%iy = iy                  
                end if
                if (rank == 0) then
                    write (*,fmt='(10(a,i4))') " Lib_ImagingSpace::setMyBounds info - proc ",pp,    &
                                " x ",ix*Cx,":",min( this%Nx-1,ix*Cx + Cx - 1 ),                 &
                                " y ",iy*Cy,":",min( this%Ny-1,iy*Cy + Cy - 1 )                 
                end if
            end do

            call setMyNeighbours( this )


            return
        end subroutine setMyBounds


        subroutine setMyNeighbours( this )
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      For send/recv operations, the ordering of the neighbours is
    !*
    !*        8   1   2     y
    !*          +---+       ^
    !*        7 |   | 3     |
    !*          +---+       +---->x
    !*        6   5   4      
    !*      return processor number for neighbours 1:8
            type(ImagingSpace),intent(inout)        ::      this
            !integer,dimension(8),intent(out)        ::      neigh

            integer             ::      jx,jy

        !---    no periodic boundaries. Assume that I hold the edge, unless proved otherwise.
            this%neigh(:) = rank             
            ! jx = this%ix     ; jy = this%iy + 1
            ! if ( (jx>=0).and.(jx<this%Mx).and.(jy>=0).and.(jy<this%My) ) this%neigh(N)  = this%p( jx,jy )
            ! jx = this%ix + 1 ; jy = this%iy + 1 
            ! if ( (jx>=0).and.(jx<this%Mx).and.(jy>=0).and.(jy<this%My) ) this%neigh(NE) = this%p( jx,jy )
            ! jx = this%ix + 1 ; jy = this%iy 
            ! if ( (jx>=0).and.(jx<this%Mx).and.(jy>=0).and.(jy<this%My) ) this%neigh(E)  = this%p( jx,jy )
            ! jx = this%ix + 1 ; jy = this%iy - 1
            ! if ( (jx>=0).and.(jx<this%Mx).and.(jy>=0).and.(jy<this%My) ) this%neigh(SE) = this%p( jx,jy )
            ! jx = this%ix     ; jy = this%iy - 1
            ! if ( (jx>=0).and.(jx<this%Mx).and.(jy>=0).and.(jy<this%My) ) this%neigh(S)  = this%p( jx,jy )
            ! jx = this%ix - 1 ; jy = this%iy - 1
            ! if ( (jx>=0).and.(jx<this%Mx).and.(jy>=0).and.(jy<this%My) ) this%neigh(SW) = this%p( jx,jy )
            ! jx = this%ix - 1 ; jy = this%iy  
            ! if ( (jx>=0).and.(jx<this%Mx).and.(jy>=0).and.(jy<this%My) ) this%neigh(W)  = this%p( jx,jy )
            ! jx = this%ix - 1 ; jy = this%iy + 1
            ! if ( (jx>=0).and.(jx<this%Mx).and.(jy>=0).and.(jy<this%My) ) this%neigh(NW) = this%p( jx,jy )

            jx = mod( this%ix     + this%Mx , this%Mx ) ; jy = mod( this%iy + 1 + this%My , this%My )
            this%neigh(N)  = this%p( jx,jy )
            jx = mod( this%ix + 1 + this%Mx , this%Mx ) ; jy = mod( this%iy + 1 + this%My , this%My )
            this%neigh(NE) = this%p( jx,jy )
            jx = mod( this%ix + 1 + this%Mx , this%Mx ) ; jy = mod( this%iy     + this%My , this%My )
            this%neigh(E)  = this%p( jx,jy )
            jx = mod( this%ix + 1 + this%Mx , this%Mx ) ; jy = mod( this%iy - 1 + this%My , this%My )
            this%neigh(SE) = this%p( jx,jy )
            jx = mod( this%ix     + this%Mx , this%Mx ) ; jy = mod( this%iy - 1 + this%My , this%My )
            this%neigh(S)  = this%p( jx,jy )
            jx = mod( this%ix - 1 + this%Mx , this%Mx ) ; jy = mod( this%iy - 1 + this%My , this%My )
            this%neigh(SW) = this%p( jx,jy )
            jx = mod( this%ix - 1 + this%Mx , this%Mx ) ; jy = mod( this%iy     + this%My , this%My )
            this%neigh(W)  = this%p( jx,jy )
            jx = mod( this%ix - 1 + this%Mx , this%Mx ) ; jy = mod( this%iy + 1 + this%My , this%My )
            this%neigh(NW) = this%p( jx,jy )

        !---    if the corner neighbours are the same as one of the edge neighbours, then I don't need to sendrecv. So set neigh to self, and no excnahge is sone
            if ( (this%neigh(NE) == this%neigh(N)).or.(this%neigh(NE) == this%neigh(E)) ) this%neigh(NE) = rank
            if ( (this%neigh(NW) == this%neigh(N)).or.(this%neigh(NW) == this%neigh(W)) ) this%neigh(NW) = rank
            if ( (this%neigh(SE) == this%neigh(S)).or.(this%neigh(SE) == this%neigh(E)) ) this%neigh(SE) = rank
            if ( (this%neigh(SW) == this%neigh(S)).or.(this%neigh(SW) == this%neigh(W)) ) this%neigh(SW) = rank


            write(*,fmt='(a,i6,a,8i6)') "Lib_ImagingSpace::setMyNeighbours rank ",rank," neighbours ",this%neigh

            ! print *,"setMyNeighbours ",rank,this%neigh
            ! neigh(1) = this%p( mod(this%ix+this%Mx+0,this%Mx),mod(this%iy+this%My+1,this%My) )
            ! neigh(2) = this%p( mod(this%ix+this%Mx+1,this%Mx),mod(this%iy+this%My+1,this%My) )
            ! neigh(3) = this%p( mod(this%ix+this%Mx+1,this%Mx),mod(this%iy+this%My+0,this%My) )
            ! neigh(4) = this%p( mod(this%ix+this%Mx+1,this%Mx),mod(this%iy+this%My-1,this%My) )
            ! neigh(5) = this%p( mod(this%ix+this%Mx+0,this%Mx),mod(this%iy+this%My-1,this%My) )
            ! neigh(6) = this%p( mod(this%ix+this%Mx-1,this%Mx),mod(this%iy+this%My-1,this%My) )
            ! neigh(7) = this%p( mod(this%ix+this%Mx-1,this%Mx),mod(this%iy+this%My+0,this%My) )
            ! neigh(8) = this%p( mod(this%ix+this%Mx-1,this%Mx),mod(this%iy+this%My+1,this%My) )

            return

        end subroutine setMyNeighbours

        pure integer function whoseBlock( this,ix,iy )
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      which processor is responsible for block (ix,iy) = (0:Mx-1,0:My-1)
            type(ImagingSpace),intent(in)                                       ::      this
            integer,intent(in)                                                  ::      ix,iy
            whoseBlock = this%p(ix,iy)
            return
        end function whoseBlock


        pure subroutine blockExtent( this,p,lbx,ubx,lby,uby )
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      return the extent of the block held by processor p
            type(ImagingSpace),intent(in)                                       ::      this
            integer,intent(in)                                                  ::      p
            integer,intent(out)                                                 ::      lbx,ubx,lby,uby

            integer                 ::  ix,iy,Cx,Cy
            integer                 ::  pp

            Cx = ceiling(real(this%Nx)/this%Mx)
            Cy = ceiling(real(this%Ny)/this%My)

            do pp = 0,nProcs-1
                iy = pp/this%Mx
                ix = pp - this%Mx*iy
                if (pp == p) then
                    lbx = ix*Cx
                    ubx = min( this%Nx-1,this%lbx + Cx - 1 )
                    lby = iy*Cy
                    uby = min( this%Ny-1,this%lby + Cy - 1 )  
                end if
            end do
            return
        end subroutine blockExtent



        subroutine sendrecv( this,phi,border )
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      complete the update of the diffracted electron beams at slice z
    !*      by exchanging information across processors
    !*      For send/recv operations, the ordering of the neighbours is
    !*
    !*        8   1   2     y
    !*          +---+       ^
    !*        7 |   | 3     |
    !*          +---+       +---->x
    !*        6   5   4                  
            type(ImagingSpace),intent(in)                                       ::      this
            complex(kind=real64),dimension(:,:,:),pointer,intent(inout)         ::      phi                 !   (0:nG,lbx-border:ubx+border,lby-border:uby+border)      a 2d slice of phi, of which (lbx+1:ubx-1) is interior
            integer,intent(in)                                                  ::      border
            integer                         ::      nG ,dx,dy! ,lbx,ubx,lby,uby

            logical,save                                                        ::      firstcall = .true.
            complex(kind=real64),dimension(:,:,:),allocatable,save              ::      phi_send_N , phi_recv_N 
            complex(kind=real64),dimension(:,:,:),allocatable,save              ::      phi_send_E , phi_recv_S 
            complex(kind=real64),dimension(:,:,:),allocatable,save              ::      phi_send_S , phi_recv_E 
            complex(kind=real64),dimension(:,:,:),allocatable,save              ::      phi_send_W , phi_recv_W 
            complex(kind=real64),dimension(:,:,:),allocatable,save              ::      phi_send_NE, phi_recv_NE
            complex(kind=real64),dimension(:,:,:),allocatable,save              ::      phi_send_SE, phi_recv_SE
            complex(kind=real64),dimension(:,:,:),allocatable,save              ::      phi_send_SW, phi_recv_SW
            complex(kind=real64),dimension(:,:,:),allocatable,save              ::      phi_send_NW, phi_recv_NW


#ifdef MPI            
            integer                             ::      ierror
            integer,dimension(8),parameter      ::      tag     = 10000 + (/ N ,NE,E ,SE,S ,SW,W ,NW /)         
            type(MPI_Request),dimension(16)     ::      request
#endif
            if (nProcs==1) return



          
        !---    determine size of problem
            nG  = size( phi,dim=1 )
            dx = this%ubx+1-this%lbx
            dy = this%uby+1-this%lby

            !print *,"sendrecv rank ",rank," start ",size(phi,dim=2),size(phi,dim=3),lbound(phi,dim=2),ubound(phi,dim=2),lbound(phi,dim=3),ubound(phi,dim=3),dx,dy


        !---    allocate explicit buffers if necessary
            if (firstcall) then
                if (this%neigh(N) /= rank) then
                    allocate(phi_send_N(nG,dx,border))
                    allocate(phi_recv_N(nG,dx,border))
                end if
                if (this%neigh(E) /= rank) then
                    allocate(phi_send_E(nG,border,dy))
                    allocate(phi_recv_E(nG,border,dy))
                end if
                if (this%neigh(S) /= rank) then
                    allocate(phi_send_S(nG,dx,border))
                    allocate(phi_recv_S(nG,dx,border))
                end if
                if (this%neigh(W) /= rank) then
                    allocate(phi_send_W(nG,border,dy))
                    allocate(phi_recv_W(nG,border,dy))
                end if
                if (this%neigh(NE)/= rank) then
                    allocate(phi_send_NE(nG,border,border))
                    allocate(phi_recv_NE(nG,border,border))
                end if
                if (this%neigh(SE)/= rank) then
                    allocate(phi_send_SE(nG,border,border))
                    allocate(phi_recv_SE(nG,border,border))
                end if
                if (this%neigh(SW)/= rank) then
                    allocate(phi_send_SW(nG,border,border))
                    allocate(phi_recv_SW(nG,border,border))
                end if
                if (this%neigh(NW)/= rank) then
                    allocate(phi_send_NW(nG,border,border))
                    allocate(phi_recv_NW(nG,border,border))
                end if
                firstcall = .false.
            end if

 

#ifdef MPI
             
   
            request = MPI_REQUEST_NULL


        !---    post the receives 
            if (this%neigh(N) /= rank) call MPI_IRECV( phi_recv_N  , size(phi_recv_N ) , MPI_DOUBLE_COMPLEX,this%neigh(N) ,tag(N) ,MPI_COMM_WORLD,request(1), ierror )
            if (this%neigh(E) /= rank) call MPI_IRECV( phi_recv_E  , size(phi_recv_E ) , MPI_DOUBLE_COMPLEX,this%neigh(E) ,tag(E) ,MPI_COMM_WORLD,request(2), ierror )
            if (this%neigh(S) /= rank) call MPI_IRECV( phi_recv_S  , size(phi_recv_S ) , MPI_DOUBLE_COMPLEX,this%neigh(S) ,tag(S) ,MPI_COMM_WORLD,request(3), ierror )
            if (this%neigh(W) /= rank) call MPI_IRECV( phi_recv_W  , size(phi_recv_W ) , MPI_DOUBLE_COMPLEX,this%neigh(W) ,tag(W) ,MPI_COMM_WORLD,request(4), ierror )
            if (this%neigh(NE)/= rank) call MPI_IRECV( phi_recv_NE , size(phi_recv_NE) , MPI_DOUBLE_COMPLEX,this%neigh(NE),tag(NE),MPI_COMM_WORLD,request(5), ierror )
            if (this%neigh(SE)/= rank) call MPI_IRECV( phi_recv_SE , size(phi_recv_SE) , MPI_DOUBLE_COMPLEX,this%neigh(SE),tag(SE),MPI_COMM_WORLD,request(6), ierror )
            if (this%neigh(SW)/= rank) call MPI_IRECV( phi_recv_SW , size(phi_recv_SW) , MPI_DOUBLE_COMPLEX,this%neigh(SW),tag(SW),MPI_COMM_WORLD,request(7), ierror )
            if (this%neigh(NW)/= rank) call MPI_IRECV( phi_recv_NW , size(phi_recv_NW) , MPI_DOUBLE_COMPLEX,this%neigh(NW),tag(NW),MPI_COMM_WORLD,request(8), ierror )

        !---    place information to send into send buffers
            if (this%neigh(N) /= rank) phi_send_N (:,:,:) = phi(:,this%lbx:this%ubx,this%uby-border+1:this%uby)        
            if (this%neigh(E) /= rank) phi_send_E (:,:,:) = phi(:,this%ubx-border+1:this%ubx,this%lby:this%uby)        
            if (this%neigh(S) /= rank) phi_send_S (:,:,:) = phi(:,this%lbx:this%ubx,this%lby:this%lby+border-1)        
            if (this%neigh(W) /= rank) phi_send_W (:,:,:) = phi(:,this%lbx:this%lbx+border-1,this%lby:this%uby)        
            if (this%neigh(NE)/= rank) phi_send_NE(:,:,:) = phi(:,this%ubx-border+1:this%ubx,this%uby-border+1:this%uby)
            if (this%neigh(SE)/= rank) phi_send_SE(:,:,:) = phi(:,this%ubx-border+1:this%ubx,this%lby:this%lby+border-1)
            if (this%neigh(SW)/= rank) phi_send_SW(:,:,:) = phi(:,this%lbx:this%lbx+border-1,this%lby:this%lby+border-1)
            if (this%neigh(NW)/= rank) phi_send_NW(:,:,:) = phi(:,this%lbx:this%lbx+border-1,this%uby-border+1:this%uby)

        !---    post the sends
            if (this%neigh(N) /= rank) call MPI_ISEND( phi_send_N  , size(phi_send_N ) , MPI_DOUBLE_COMPLEX,this%neigh(N) ,tag(S) ,MPI_COMM_WORLD,request(9) , ierror )
            if (this%neigh(E) /= rank) call MPI_ISEND( phi_send_E  , size(phi_send_E ) , MPI_DOUBLE_COMPLEX,this%neigh(E) ,tag(W) ,MPI_COMM_WORLD,request(10), ierror )
            if (this%neigh(S) /= rank) call MPI_ISEND( phi_send_S  , size(phi_send_S ) , MPI_DOUBLE_COMPLEX,this%neigh(S) ,tag(N) ,MPI_COMM_WORLD,request(11), ierror )
            if (this%neigh(W) /= rank) call MPI_ISEND( phi_send_W  , size(phi_send_W ) , MPI_DOUBLE_COMPLEX,this%neigh(W) ,tag(E) ,MPI_COMM_WORLD,request(12), ierror )
            if (this%neigh(NE)/= rank) call MPI_ISEND( phi_send_NE , size(phi_send_NE) , MPI_DOUBLE_COMPLEX,this%neigh(NE),tag(SW),MPI_COMM_WORLD,request(13), ierror )
            if (this%neigh(SE)/= rank) call MPI_ISEND( phi_send_SE , size(phi_send_SE) , MPI_DOUBLE_COMPLEX,this%neigh(SE),tag(NW),MPI_COMM_WORLD,request(14), ierror )
            if (this%neigh(SW)/= rank) call MPI_ISEND( phi_send_SW , size(phi_send_SW) , MPI_DOUBLE_COMPLEX,this%neigh(SW),tag(NE),MPI_COMM_WORLD,request(15), ierror )
            if (this%neigh(NW)/= rank) call MPI_ISEND( phi_send_NW , size(phi_send_NW) , MPI_DOUBLE_COMPLEX,this%neigh(NW),tag(SE),MPI_COMM_WORLD,request(16), ierror )

        !---    wait for messages to be received.
            call MPI_WAITALL( size(request),request,MPI_STATUSES_IGNORE,ierror )

        !---    unpack the received information 
            if (this%neigh(N) /= rank) phi(:,this%lbx:this%ubx,this%uby+1:this%uby+border)          = phi_recv_N (:,:,:)
            if (this%neigh(E) /= rank) phi(:,this%ubx+1:this%ubx+border,this%lby:this%uby)          = phi_recv_E (:,:,:)
            if (this%neigh(S) /= rank) phi(:,this%lbx:this%ubx,this%lby-border:this%lby-1)          = phi_recv_S (:,:,:)
            if (this%neigh(W) /= rank) phi(:,this%lbx-border:this%lbx-1,this%lby:this%uby)          = phi_recv_W (:,:,:)
            if (this%neigh(NE)/= rank) phi(:,this%ubx+1:this%ubx+border,this%uby+1:this%uby+border) = phi_recv_NE(:,:,:)
            if (this%neigh(SE)/= rank) phi(:,this%ubx+1:this%ubx+border,this%lby-border:this%lby-1) = phi_recv_SE(:,:,:)
            if (this%neigh(SW)/= rank) phi(:,this%lbx-border:this%lbx-1,this%lby-border:this%lby-1) = phi_recv_SW(:,:,:)
            if (this%neigh(NW)/= rank) phi(:,this%lbx-border:this%lbx-1,this%uby+1:this%uby+border) = phi_recv_NW(:,:,:)


 



            !call MPI_BARRIER( MPI_COMM_WORLD,ierror )
            !print *,"sendrecv rank ",rank," done "
#endif
  
            return
        end subroutine sendrecv

!******************************************************************************
!
!           Accessors
!
!******************************************************************************



        
        pure real(kind=real64) function geta0(this) 
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
            type(ImagingSpace),intent(in)          ::      this
            geta0 = this%a
            return
        end function geta0

        pure real(kind=real64) function getsigma0(this) 
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
            type(ImagingSpace),intent(in)          ::      this
            getsigma0 = this%sigma
            return
        end function getsigma0        

        pure function getdelta0(this) result (delta)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
            type(ImagingSpace),intent(in)           ::      this
            real(kind=real64),dimension(3)          ::      delta
            delta = this%delta
            return
        end function getdelta0
        
        pure integer function getnBuf0(this) 
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
            type(ImagingSpace),intent(in)          ::      this
            getnBuf0 = this%nBuf
            return
        end function getnBuf0

        pure integer function getNx0(this) 
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
            type(ImagingSpace),intent(in)          ::      this
            getNx0 = this%Nx
            return
        end function getNx0

        pure integer function getNy0(this) 
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
            type(ImagingSpace),intent(in)          ::      this
            getNy0 = this%Ny
            return
        end function getNy0

        pure integer function getNz0(this) 
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
            type(ImagingSpace),intent(in)          ::      this
            getNz0 = this%Nz
            return
        end function getNz0

        pure integer function getMx0(this) 
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
            type(ImagingSpace),intent(in)          ::      this
            getMx0 = this%Mx
            return
        end function getMx0

        pure integer function getMy0(this) 
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
            type(ImagingSpace),intent(in)          ::      this
            getMy0 = this%My
            return
        end function getMy0


        pure integer function getlbx(this) 
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
            type(ImagingSpace),intent(in)          ::      this
            getlbx = this%lbx
            return
        end function getlbx
 
        pure integer function getubx(this) 
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
            type(ImagingSpace),intent(in)          ::      this
            getubx = this%ubx
            return
        end function getubx
 
        pure integer function getlby(this) 
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
            type(ImagingSpace),intent(in)          ::      this
            getlby = this%lby
            return
        end function getlby
 
        pure integer function getuby(this) 
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
            type(ImagingSpace),intent(in)          ::      this
            getuby = this%uby
            return
        end function getuby
 
        pure subroutine getBounds0(this,b) 
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
            type(ImagingSpace),intent(in)           ::      this
            integer,dimension(2,3),intent(out)      ::      b
            b(:,1) = (/ this%lbx,this%ubx /)
            b(:,2) = (/ this%lby,this%uby /)
            b(:,3) = (/ 0,this%Nz-1 /)
            return
        end subroutine getBounds0
 
        pure subroutine getBounds1(this,lbx,ubx,lby,uby) 
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
            type(ImagingSpace),intent(in)           ::      this
            integer,intent(out)                     ::      lbx,ubx,lby,uby
            lbx = this%lbx
            ubx = this%ubx
            lby = this%lby
            uby = this%uby
            return
        end subroutine getBounds1



!******************************************************************************
!
!           static functions
!
!******************************************************************************




        pure integer function buffer(a,sigma)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      given the range of the broadening, return the buffer range of cells required to 
    !*      catch all the broadening
            real(kind=real64),intent(in)                ::      a
            real(kind=real64),intent(in)                ::      sigma
            buffer = ceiling( Lib_ImagingSpace_SIGMA_MULT*sigma / a )
            return
        end function buffer



        
            
        subroutine Lib_ImagingSpace_init_MPI()
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
#ifdef MPI
            integer             ::      ierror
            logical             ::      ok
            call MPI_Initialized(ok,ierror)
            if (.not. ok) call MPI_INIT(ierror)
            call MPI_COMM_SIZE(MPI_COMM_WORLD, nProcs, ierror)
            call MPI_COMM_RANK(MPI_COMM_WORLD, rank, ierror)
            if (rank == 0) print *,"Lib_ImagingSpace::Lib_ImagingSpace_init_MPI info - initialised with ",nProcs," processes"  
#else
            print *,"Lib_ImagingSpace::Lib_ImagingSpace_init_MPI info - serial mode"  
#endif            
            return
        end subroutine Lib_ImagingSpace_init_MPI

                      
    end module Lib_ImagingSpace