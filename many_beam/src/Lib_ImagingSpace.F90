

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
            uu = 6 ; if (present(u)) uu = u
            oo = 0 ; if (present(o)) oo = o
            write(unit=uu,fmt='(2(a,f16.8),5(a,i8))') repeat(" ",oo)//"ImagingSpace [a,sigma= ",this%a,",",this%sigma," , Nx,Ny,Nz= ",this%Nx,",",this%Ny,",",this%Nz,", nBuf=",this%nBuf," ]"
            write(unit=uu,fmt='(6(a,i8))') repeat(" ",oo+4)//"my block ",this%lbx,":",this%ubx," , ",this%lby,":",this%uby," , ",0,":",this%Nz-1
            return
        end subroutine report0




!******************************************************************************
!
!           member functions
!
!******************************************************************************

        
        pure logical function inMyCell0(this,xt,buffered)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      given the reduced space position
    !*          xt = R (x-delta)/a
    !*      return true if the atom is in a cell I am responsible for.
    !*      optionally return true if within the buffered region
            type(ImagingSpace),intent(in)               ::      this
            real(kind=real64),dimension(:),intent(in)   ::      xt
            logical,intent(in)                          ::      buffered
            integer         ::          ix,iy,iz

            if (size(xt)==2) then
                !   2d test - check xt(1:2) is in the x-y bounds, but do not test for z
                
                ix = floor( xt(1) )            
                iy = floor( xt(2) )

                if (buffered) then
                    inMyCell0 = (ix>=this%lbx-this%nBuf) .and. (ix<=this%ubx+this%nBuf) 
                    inMyCell0 = inMyCell0 .and. (iy>=this%lby-this%nBuf) .and. (iy<=this%uby+this%nBuf) 
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
                    inMyCell0 = (ix>=this%lbx-this%nBuf) .and. (ix<=this%ubx+this%nBuf) 
                    inMyCell0 = inMyCell0 .and. (iy>=this%lby-this%nBuf) .and. (iy<=this%uby+this%nBuf) 
                    inMyCell0 = inMyCell0 .and. (iz>=-this%nBuf) .and. (iz<this%Nz+this%nBuf) 
                else
                    inMyCell0 = (ix>=this%lbx) .and. (ix<=this%ubx) 
                    inMyCell0 = inMyCell0 .and. (iy>=this%lby) .and. (iy<=this%uby) 
                    inMyCell0 = inMyCell0 .and. (iz>=0) .and. (iz<this%Nz) 
                end if
            end if

            return
        end function inMyCell0

        pure logical function inMyCell1(this,xt,buffered)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      given the reduced space position
    !*          xt = R (x-delta)/a
    !*      return true if the atom is in a cell I am responsible for.
    !*      optionally return true if within the buffered region
            type(ImagingSpace),intent(in)               ::      this
            real(kind=real32),dimension(:),intent(in)   ::      xt
            logical,intent(in)                          ::      buffered
            inMyCell1 = inMyCell0(this,real(xt,kind=real64),buffered)

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
            jx = this%ix     ; jy = this%iy + 1
            if ( (jx>=0).and.(jx<this%Mx).and.(jy>=0).and.(jy<this%My) ) this%neigh(N)  = this%p( jx,jy )
            jx = this%ix + 1 ; jy = this%iy + 1 
            if ( (jx>=0).and.(jx<this%Mx).and.(jy>=0).and.(jy<this%My) ) this%neigh(NE) = this%p( jx,jy )
            jx = this%ix + 1 ; jy = this%iy 
            if ( (jx>=0).and.(jx<this%Mx).and.(jy>=0).and.(jy<this%My) ) this%neigh(E)  = this%p( jx,jy )
            jx = this%ix + 1 ; jy = this%iy - 1
            if ( (jx>=0).and.(jx<this%Mx).and.(jy>=0).and.(jy<this%My) ) this%neigh(SE) = this%p( jx,jy )
            jx = this%ix     ; jy = this%iy - 1
            if ( (jx>=0).and.(jx<this%Mx).and.(jy>=0).and.(jy<this%My) ) this%neigh(S)  = this%p( jx,jy )
            jx = this%ix - 1 ; jy = this%iy - 1
            if ( (jx>=0).and.(jx<this%Mx).and.(jy>=0).and.(jy<this%My) ) this%neigh(SW) = this%p( jx,jy )
            jx = this%ix - 1 ; jy = this%iy  
            if ( (jx>=0).and.(jx<this%Mx).and.(jy>=0).and.(jy<this%My) ) this%neigh(W)  = this%p( jx,jy )
            jx = this%ix - 1 ; jy = this%iy + 1
            if ( (jx>=0).and.(jx<this%Mx).and.(jy>=0).and.(jy<this%My) ) this%neigh(NW) = this%p( jx,jy )


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



        subroutine sendrecv( this,phi )
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
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
            complex(kind=real64),dimension(:,:,:),pointer,intent(inout)         ::      phi                 !   (0:nG,lbx:ubx,lby:uby)      a 2d slice of phi, of which (lbx+1:ubx-1) is interior
            integer                         ::      nG,lbx,ubx,lby,uby

#ifdef MPI            
            integer                         ::      ierror
            integer,dimension(8)            ::      tag           
            type(MPI_Request),dimension(8)  ::      request
#endif

        !---    determine size of problem
            nG  = ubound( phi,dim=1 )
            lbx = lbound( phi,dim=2 )
            ubx = ubound( phi,dim=2 )
            lby = lbound( phi,dim=3 )
            uby = ubound( phi,dim=3 )

            print *,"sendrecv rank ",rank,(ubx-lbx-1),(uby-lby-1)

        !---    complete the boundary update where I hold edge, by interpolating out
            if (this%neigh(S) == rank) phi(:,lbx+1:ubx-1,lby) = 3*phi(:,lbx+1:ubx-1,lby+1) - 3*phi(:,lbx+1:ubx-1,lby+2) + phi(:,lbx+1:ubx-1,lby+3) 
            if (this%neigh(N) == rank) phi(:,lbx+1:ubx-1,uby) = phi(:,lbx+1:ubx-1,uby-3) - 3*phi(:,lbx+1:ubx-1,uby-2) + 3*phi(:,lbx+1:ubx-1,uby-1)
            if (this%neigh(W) == rank) phi(:,lbx,lby+1:uby-1) = 3*phi(:,lbx+1,lby+1:uby-1) - 3*phi(:,lbx+1,lby+1:uby-1) + phi(:,lbx+2,lby+1:uby-1)
            if (this%neigh(E) == rank) phi(:,ubx,lby+1:uby-1) = phi(:,ubx-3,lby+1:uby-1) - 3*phi(:,ubx-2,lby+1:uby-1) + 3*phi(:,ubx-1,lby+1:uby-1)



#ifdef MPI
            
            request = MPI_REQUEST_NULL
        !---    find my neighbours
            !call myNeighbours( this, neigh )            
 
            tag(1:8) = 10000 + (/ N ,NE,E ,SE,S ,SW,W ,NW /)

            if (this%neigh(N) /= rank) print *,"sendrecv rank ",rank," post recv from ",this%neigh(N),tag(N)
            if (this%neigh(E) /= rank) print *,"sendrecv rank ",rank," post recv from ",this%neigh(E),tag(E)
            if (this%neigh(S) /= rank) print *,"sendrecv rank ",rank," post recv from ",this%neigh(S),tag(S)
            if (this%neigh(W) /= rank) print *,"sendrecv rank ",rank," post recv from ",this%neigh(W),tag(W)


        !---    post the receives
            if (this%neigh(N) /= rank) call MPI_IRECV( phi(:,lbx+1:ubx-1,uby),(nG+1)*(ubx-lbx-1),MPI_DOUBLE_COMPLEX,this%neigh(N),tag(N),MPI_COMM_WORLD,request(1), ierror )
            if (this%neigh(E) /= rank) call MPI_IRECV( phi(:,ubx,lby+1:uby-1),(nG+1)*(uby-lby-1),MPI_DOUBLE_COMPLEX,this%neigh(E),tag(E),MPI_COMM_WORLD,request(2), ierror )
            if (this%neigh(S) /= rank) call MPI_IRECV( phi(:,lbx+1:ubx-1,lby),(nG+1)*(ubx-lbx-1),MPI_DOUBLE_COMPLEX,this%neigh(S),tag(S),MPI_COMM_WORLD,request(3), ierror )
            if (this%neigh(W) /= rank) call MPI_IRECV( phi(:,lbx,lby+1:uby-1),(nG+1)*(uby-lby-1),MPI_DOUBLE_COMPLEX,this%neigh(W),tag(W),MPI_COMM_WORLD,request(4), ierror )


        !---    post the sends
    !*        8   1   2     y
    !*          +---+       ^
    !*        7 |   | 3     |
    !*          +---+       +---->x
    !*        6   5   4  
            if (this%neigh(N) /= rank) print *,"sendrecv rank ",rank," post send to ",this%neigh(N),tag(S)
            if (this%neigh(E) /= rank) print *,"sendrecv rank ",rank," post send to ",this%neigh(E),tag(W)
            if (this%neigh(S) /= rank) print *,"sendrecv rank ",rank," post send to ",this%neigh(S),tag(N)
            if (this%neigh(W) /= rank) print *,"sendrecv rank ",rank," post send to ",this%neigh(W),tag(E)


            if (this%neigh(N) /= rank) call MPI_ISEND( phi(:,lbx+1:ubx-1,lby+1),(nG+1)*(ubx-lbx-1),MPI_DOUBLE_COMPLEX,this%neigh(N),tag(S),MPI_COMM_WORLD,request(5) ,ierror )
            if (this%neigh(E) /= rank) call MPI_ISEND( phi(:,lbx+1,lby+1:uby-1),(nG+1)*(uby-lby-1),MPI_DOUBLE_COMPLEX,this%neigh(E),tag(W),MPI_COMM_WORLD,request(6), ierror )
            if (this%neigh(S) /= rank) call MPI_ISEND( phi(:,lbx+1:ubx-1,uby-1),(nG+1)*(ubx-lbx-1),MPI_DOUBLE_COMPLEX,this%neigh(S),tag(N),MPI_COMM_WORLD,request(7), ierror )
            if (this%neigh(W) /= rank) call MPI_ISEND( phi(:,ubx-1,lby+1:uby-1),(nG+1)*(uby-lby-1),MPI_DOUBLE_COMPLEX,this%neigh(W),tag(E),MPI_COMM_WORLD,request(8), ierror )



        !---    wait for the send=recv to complete
            call MPI_WAITALL( size(request),request,MPI_STATUSES_IGNORE,ierror )
#endif
 

        !---    complete the corners
            print *,"sendrecv rank ",rank," corners"
            phi(:,ubx,uby) = ( phi(:,ubx,uby-3) - 3*phi(:,ubx,uby-2) + 3*phi(:,ubx,uby-1) + phi(:,ubx-3,uby) - 3*phi(:,ubx-2,uby) + 3*phi(:,ubx-1,uby) )/2
            phi(:,ubx,lby) = ( 3*phi(:,ubx,lby+1) - 3*phi(:,ubx,lby+2) + phi(:,ubx,lby+3) + phi(:,ubx-3,lby) - 3*phi(:,ubx-2,lby) + 3*phi(:,ubx-1,lby) )/2
            phi(:,lbx,lby) = ( 3*phi(:,lbx,lby+1) - 3*phi(:,lbx,lby+2) + phi(:,lbx,lby+3) + 3*phi(:,lbx+1,lby) - 3*phi(:,lbx+2,lby) + phi(:,lbx+3,lby) )/2
            phi(:,lbx,uby) = ( phi(:,lbx,uby-3) - 3*phi(:,lbx,uby-2) + 3*phi(:,lbx,uby-1) + 3*phi(:,lbx+1,uby) - 3*phi(:,lbx+2,uby) + phi(:,lbx+3,uby) )/2
            print *,"sendrecv rank ",rank," done"


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