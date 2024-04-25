
    module Lib_LinkCell3D
!---^^^^^^^^^^^^^^^^^^^^^^
!*      A very simple, lightweight link cell implementation.
!*      assumes pbc, so add buffer cells if necessary.
!*      assumes cell size > max neighbour cutoff range, so searches 27 neighbour cells only.
!*      note: use neighbourList_long if search range > cell size
        use Lib_SimpleSupercells
        use iso_fortran_env
        implicit none
        private
        
    !---
    
        public      ::      LinkCell3D_ctor
        public      ::      delete
        public      ::      report
     
        public      ::      clone       !   deep copy with allocate
        public      ::      clear
        public      ::      add
        public      ::      cut
        public      ::      move
        
        public      ::      neighbourList
        public      ::      neighbourList_long
        public      ::      getnNeighMax
        public      ::      getnPoints
        public      ::      getnMax
        public      ::      getMaxPerCell
        public      ::      getPosition,getId
        public      ::      getSuper,getA
        public      ::      getNx,getNy,getNz
        public      ::      translate
        public      ::      getMaxId
        
        public      ::      estimateMaxCountPerCell
        public      ::      estimateMaxCountPerVolume
        public      ::      nearestNeighbour

    !---
    
        logical,public          ::      LinkCell3D_dbg = .false.
        integer(kind=int32),public,parameter        ::      NOTANID = -987654321    
        
    !---
    
        type        ::      Cell
            integer                                     ::      nMax    !   max count of points in cell
            integer                                     ::      nindx   !   count of points in each cell
            integer,dimension(:),pointer                ::      indx    !   (1:nMax)   
            real(kind=real64),dimension(:,:),pointer    ::      x       !   (3,nMax) reduced cell coords position wrt lfb of points in cell
        end type Cell
    
        type,public     ::      LinkCell3D
            private
            type(SimpleSupercell)                   ::      super
            integer                                 ::      nMax        !   max count of points in cell
            integer                                 ::      nMax0       !   sanity testing - I'll reallocate indx, but only up to 16x this value. After that I'm assuming there's a bug.
            
            type(Cell),dimension(:,:,:),pointer     ::      c
            
        end type LinkCell3D
        
    !---
    
        interface Cell_ctor
            module procedure    Cell_null
            module procedure    Cell_ctor0
        end interface
        
        interface LinkCell3D_ctor
            module procedure    LinkCell3D_null
            module procedure    LinkCell3D_ctor0
            module procedure    LinkCell3D_ctor0a
            module procedure    LinkCell3D_ctor1            
        end interface
                
        interface delete
            module procedure    delete0
            module procedure    delete1
        end interface
        
        interface report
            module procedure    report0
            module procedure    report1
        end interface
        
        interface clone
            module procedure    clone0
            module procedure    clone1
        end interface
        
        interface clear
            module procedure    clear0
            module procedure    clear1
        end interface
        
        interface add
            module procedure    add0
            module procedure    add1
        end interface
        
        interface cut
            module procedure    cut0
        end interface
         
        interface move
            module procedure    move0
        end interface
         
        interface getPosition
            module procedure    getPosition0
        end interface
         
         
        interface getSuper
            module procedure    getSuper0
        end interface
        interface getA
            module procedure    getA0
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
         
        interface translate
            module procedure    translate0
        end interface
         
         
        interface neighbourList
            module procedure    neighbourList0
            module procedure    neighbourList0a
!            module procedure    neighbourList0c
            module procedure    neighbourList1
            module procedure    neighbourList1a
            module procedure    neighbourList1b
            module procedure    neighbourList1c
            module procedure    neighbourList2c
            module procedure    neighbourList1d
            module procedure    neighbourList1e
            module procedure    neighbourList2             
            module procedure    neighbourList2b            
            module procedure    neighbourList3c
        end interface
         
         
         
        interface nearestNeighbour
            module  procedure    nearestNeighbour0
            module  procedure    nearestNeighbour1
        end interface
        
        interface neighbourList_long
            module  procedure    neighbourList_long1b
            module  procedure    neighbourList_long1c
            module  procedure    neighbourList_long1d
        end interface
         
        interface getnPoints
            module  procedure    getnPoints0
            module  procedure    getnPoints1
        end interface
         
        interface estimateMaxCountPerCell
            module  procedure    estimateMaxCountPerCell0
            module  procedure    estimateMaxCountPerCell1
        end interface
         
        
    contains
!---^^^^^^^^

        function LinkCell3D_null() result(this)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
            type(LinkCell3D)           ::      this
            this%super = SimpleSupercell_ctor()
            this%nMax0 = 0
            this%nMax = 0
            nullify(this%c)
            return
        end function LinkCell3D_null
                         
        function LinkCell3D_ctor0(a,nMax,Nx,Ny,Nz) result(this)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
            real(kind=real64),dimension(3,3),intent(in)     ::  a
            integer,intent(in)         ::      nMax,Nx,Ny,Nz
            type(LinkCell3D)           ::      this
            integer     ::      ix,iy,iz
            !real        ::      mm
            this%super = SimpleSupercell_ctor(a,Nx,Ny,Nz)
            this%nMax0 = nMax
            allocate(this%c(0:this%super%Nx-1,0:this%super%Ny-1,0:this%super%Nz-1))
            !print *,"Lib_LinkCell3D::LinkCell3D_ctor0() info - allocating ",this%super%Nx*this%super%Ny*this%super%Nz," cells with ",nMax," atoms per cell"
           
            do iz = 0,this%super%Nz-1
                do iy = 0,this%super%Ny-1
                    do ix = 0,this%super%Nx-1
                        this%c(ix,iy,iz) = Cell_ctor(nMax) 
                    end do
                end do
            end do
            this%nMax = nMax 

            return
        end function LinkCell3D_ctor0
        
        function LinkCell3D_ctor0a(nMax,Nx,Ny,Nz) result(this)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
            integer,intent(in)         ::      nMax,Nx,Ny,Nz
            type(LinkCell3D)           ::      this
            integer     ::      ix,iy,iz
            real(kind=real64),dimension(3,3)    ::      identity3x3 = reshape( (/1,0,0,0,1,0,0,0,1/),(/3,3/) )
            this%super = SimpleSupercell_ctor(identity3x3,Nx,Ny,Nz)
            this%nMax0 = nMax
            allocate(this%c(0:this%super%Nx-1,0:this%super%Ny-1,0:this%super%Nz-1))
            do iz = 0,this%super%Nz-1
                do iy = 0,this%super%Ny-1
                    do ix = 0,this%super%Nx-1
                        this%c(ix,iy,iz) = Cell_ctor(nMax) 
                    end do
                end do
            end do
            this%nMax = nMax 

            return
        end function LinkCell3D_ctor0a
        
        function LinkCell3D_ctor1(nMax,super) result(this)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
            integer,intent(in)         ::      nMax
            type(SimpleSupercell),intent(in)    ::      super
            type(LinkCell3D)           ::      this
            integer     ::      ix,iy,iz
            
            
            this%nMax0 = nMax
            this%super = super
            
            allocate(this%c(0:this%super%Nx-1,0:this%super%Ny-1,0:this%super%Nz-1))
            do iz = 0,this%super%Nz-1
                do iy = 0,this%super%Ny-1
                    do ix = 0,this%super%Nx-1
                        this%c(ix,iy,iz) = Cell_ctor(nMax) 
                    end do
                end do
            end do
            this%nMax = nMax

            return
        end function LinkCell3D_ctor1
        
    !---    
    
        
        function Cell_null() result(this)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ 
            type(Cell)           ::      this
            this%nMax = 0
            this%nindx = 0
            nullify(this%indx)
            nullify(this%x)
            return
        end function Cell_null
                         
        function Cell_ctor0(nMax) result(this)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
            integer,intent(in)         ::      nMax 
            type(Cell)           ::      this
            this%nMax = nMax
            this%nindx = 0
            allocate(this%indx(this%nMax))
            allocate(this%x(3,this%nMax))
            this%indx = NOTANID
            this%x = 0            
            return
        end function Cell_ctor0
        
    !---
    
        function estimateMaxCountPerCell0(nCells,nTotal) result(nMax)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      if nTotal atoms are distributed homogeneously among the cells
    !*      what is the maximum density per cell?
    !*      return level where we are 99% sure a cell will be under this level
    !*      nCells = this%super%Nx*this%super%Ny*this%super%Nz
    !*      but entered this way as we expect to need this function before the lc3d is constructed.
    
            integer,intent(in)                  ::      nCells,nTotal
            integer                             ::      nMax
            
            real(kind=real64)           ::      rho,pp,psum
            integer                     ::      kk
            
            rho = real(nTotal,kind=real64)/(nCells)          !   average count per cell
            
            pp = exp( - rho )
            psum = pp
            do kk = 1,nTotal
                if (psum >= 0.99d0) then               
                    nMax = kk
                    return
                end if
                pp = pp * rho / kk
                psum = psum + pp
            end do
            nMax = nTotal
            return
        end function estimateMaxCountPerCell0
            
    
        function estimateMaxCountPerCell1(nCells,nTotal) result(nMax)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      if nTotal atoms are distributed homogeneously among the cells
    !*      what is the maximum density per cell?
    !*      return level where we are 99% sure a cell will be under this level
    !*      nCells = this%super%Nx*this%super%Ny*this%super%Nz
    !*      but entered this way as we expect to need this function before the lc3d is constructed.
    
            real(kind=real64),intent(in)        ::      nCells
            integer,intent(in)                  ::      nTotal
            integer                             ::      nMax
            
            real(kind=real64)           ::      rho,pp,psum
            integer                     ::      kk
            
            rho = nTotal/nCells         !   average count per cell
            
            pp = exp( - rho )
            psum = pp
            do kk = 1,nTotal
                if (psum >= 0.99d0) then               
                    nMax = kk
                    return
                end if
                pp = pp * rho / kk
                psum = psum + pp
            end do
            nMax = nTotal
            return
        end function estimateMaxCountPerCell1
        
        function estimateMaxCountPerVolume(this,V) result(nMax)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      return level where we are 99.9% sure a volume V will have less than or equal to nMax atoms
            type(LinkCell3d),intent(in)         ::      this
            real(kind=real64),intent(in)        ::      V
            integer                             ::      nMax
            
            real(kind=real64)           ::      rho,pp,psum
            integer                     ::      kk
            
            nMax = getNpoints(this)
            rho = nMax/superCellVolume(this%super)           !   atom density
            rho = rho * V                                           !   expected count in volume V
            
        !---    
            pp = exp( - rho )
            psum = pp
            do kk = 1,nMax
                if (psum >= 0.999d0) then               
                    nMax = kk
                    return
                end if
                pp = pp * rho / kk
                psum = psum + pp
            end do
            return
        end function estimateMaxCountPerVolume
        
        
        
    !---    
        
        
        
        
        
        
        subroutine doubleCell(this)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      double storage capacity of cell
            type(Cell),intent(inout)      ::      this
            integer,dimension(:),pointer        ::      indx_tmp
            real(kind=real64),dimension(:,:),pointer        ::      x_tmp
            integer     ::      newNmax
            if (this%nMax == 0) then
                this = Cell_ctor(2)
            else                
                newNmax = max(2,this%nMax*2)
                allocate(indx_tmp(newNmax))
                allocate(x_tmp(3,newNmax))
                indx_tmp(1:this%nMax) = this%indx(1:this%nMax)
                x_tmp(1:3,1:this%nMax) = this%x(1:3,1:this%nMax)
                deallocate(this%indx)
                deallocate(this%x)
                this%indx => indx_tmp
                this%x => x_tmp
                this%nMax = newNmax
            end if
            return
        end subroutine doubleCell
        
    !     subroutine halfCell(this)   !   for future use.
    ! !---^^^^^^^^^^^^^^^^^^^^^^^^^^^
    ! !*      halve storage capacity of cell
    !         type(Cell),intent(inout)      ::      this
    !         integer,dimension(:),pointer        ::      indx_tmp
    !         real(kind=real64),dimension(:,:),pointer        ::      x_tmp
    !         integer     ::      newNmax
    !         if (this%nMax == 0) then
                 
    !         else if (this%nMax == 1) then
    !             deallocate(this%indx)
    !             deallocate(this%x)
    !             nullify(this%indx)
    !             nullify(this%x)
    !             this%nMax = 0
    !             this%nIndx = 0                 
    !         else                
    !             newNmax = (this%nMax+1)/2       !   round up
    !             allocate(indx_tmp(newNmax))
    !             allocate(x_tmp(3,newNmax))
    !             indx_tmp(1:newNmax) = this%indx(1:newNmax)
    !             x_tmp(1:3,1:newNmax) = this%x(1:3,1:newNmax)
    !             deallocate(this%indx)
    !             deallocate(this%x)
    !             this%indx => indx_tmp
    !             this%x => x_tmp
    !             this%nMax = newNmax
    !         end if
    !         return
    !     end subroutine halfCell
            
            
    !---
    
        subroutine delete0(this)
    !---^^^^^^^^^^^^^^^^^^^^^^^^
    !*      deallocate dynamic memory
            type(LinkCell3D),intent(inout)    ::      this
            integer     ::      ix,iy,iz
            do iz = 0,this%super%Nz-1
                do iy = 0,this%super%Ny-1
                    do ix = 0,this%super%Nx-1
                        call delete(this%c(ix,iy,iz))
                    end do
                end do
            end do
            deallocate(this%c)
            this = LinkCell3D_null()
            return
        end subroutine delete0
        
        subroutine delete1(this)
    !---^^^^^^^^^^^^^^^^^^^^^^^^
    !*      deallocate dynamic memory    
            type(Cell),intent(inout)    ::      this
            if (this%nMax == 0) return
            deallocate(this%indx)
            deallocate(this%x)            
            this = Cell_null()
            return
        end subroutine delete1
    !---
    
        subroutine report0(this,u,o)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      simple dump unit u ( default to screen )
    !*      with optional left hand margin o spaces
            type(LinkCell3D),intent(in)    ::      this
            integer,intent(in),optional     ::      u,o
            integer     ::      uu,oo
            real        ::      mm
            
            uu = 6 ; if (present(u)) uu = u
            oo = 0 ; if (present(o)) oo = o
            write(unit=uu,fmt='(a,i4,a,i8,4(a,i6))') &
                repeat(" ",oo)//"LinkCell3D [Nmax0/cell = ",this%Nmax0,",nPoints = ",getNpoints(this),",nNodes = ",this%super%Nx,",",this%super%Ny,",",this%super%Nz,"]"            
            mm = (this%nMax*7 + 3)*4
            mm = (mm*this%super%Nx)*this%super%Ny*this%super%Nz/(1024.0*1024.0)
            write(unit=uu,fmt='(a,f12.3,a)')  repeat(" ",oo+4)//"memory alloc ",mm," Mb"
            call report(this%super,uu,oo+4)
            return
        end subroutine report0
    
        subroutine report1(this,u,o)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      simple dump unit u ( default to screen )
    !*      with optional left hand margin o spaces
            type(Cell),intent(in)            ::      this
            integer,intent(in),optional     ::      u,o
            integer     ::      uu,oo,ii
            uu = 6 ; if (present(u)) uu = u
            oo = 0 ; if (present(o)) oo = o
            write(unit=uu,fmt='(3(a,i4))')  repeat(" ",oo)//"Cell [nIndx/nMax = ",this%nIndx,"/",this%nMax,"]"
            do ii = 1,this%nIndx
                write(unit=uu,fmt='(a,i4,i8,3f10.6)')  repeat(" ",oo+2),ii,this%indx(ii),this%x(1:3,ii)
            end do
            return
        end subroutine report1
    
    
        subroutine clone0( this,that )
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      this = that deep copy with allocate
            type(Cell),intent(inout)            ::      this
            type(Cell),intent(in)               ::      that
            if (this%nMax < that%nIndx) then
                call delete(this)
                this%nMax = that%nMax
                allocate(this%indx(this%nMax))
                allocate(this%x(3,this%nMax))
            end if
            this%nIndx = that%nIndx
            this%indx(1:this%nIndx) = that%indx(1:this%nIndx)
            this%x(1:3,1:this%nIndx) = that%x(1:3,1:this%nIndx)
            return
        end subroutine clone0
         
        subroutine clone1( this,that )
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      this = that deep copy with allocate
            type(LinkCell3D),intent(inout)            ::      this
            type(LinkCell3D),intent(in)               ::      that
            integer         ::      Nx,Ny,Nz , ix,iy,iz
            
            this%super = that%super 
            if (.not. associated(this%c)) then
                allocate(this%c(0:getNx(this%super)-1,0:getNy(this%super)-1,0:getNz(this%super)-1))
                do iz = 0,size(this%c,dim=3)-1
                    do iy = 0,size(this%c,dim=2)-1
                        do ix = 0,size(this%c,dim=1)-1
                            this%c(ix,iy,iz) = Cell_ctor()
                        end do
                    end do
                end do
            else
                Nx = size(this%c,dim=1)
                Ny = size(this%c,dim=2)
                Nz = size(this%c,dim=3)
                if (any( (/Nx,Ny,Nz/) < (/ getNx(that%super),getNy(that%super),getNz(that%super) /) )) then
                    deallocate(this%c)
                    allocate(this%c(0:getNx(this%super)-1,0:getNy(this%super)-1,0:getNz(this%super)-1))
                    do iz = 0,size(this%c,dim=3)-1
                        do iy = 0,size(this%c,dim=2)-1
                            do ix = 0,size(this%c,dim=1)-1
                                this%c(ix,iy,iz) = Cell_ctor()
                            end do
                        end do
                    end do
                end if
            end if
            this%nMax = that%nMax
            this%nMax0 = that%nMax0
            do iz = 0,size(this%c,dim=3)-1
                do iy = 0,size(this%c,dim=2)-1
                    do ix = 0,size(this%c,dim=1)-1
                        call clone( this%c(ix,iy,iz),that%c(ix,iy,iz) )
                    end do
                end do
            end do
            return
        end subroutine clone1
    !---

    
        subroutine clear0(this)
    !---^^^^^^^^^^^^^^^^^^^^^^^
    !*      remove all points from all cells. Do not deallocate memory
            type(LinkCell3D),intent(inout)   ::      this
            integer     ::      ix,iy,iz
            do iz = 0,this%super%Nz-1
                do iy = 0,this%super%Ny-1
                    do ix = 0,this%super%Nx-1
                        call clear(this%c(ix,iy,iz))
                    end do
                end do
            end do
            this%nMax = 0
            return
        end subroutine clear0
    
        subroutine clear1(this)
    !---^^^^^^^^^^^^^^^^^^^^^^^
    !*      remove all points from all cells. Do not deallocate memory
            type(Cell),intent(inout)   ::      this
            this%nindx = 0
            if (this%nMax > 0) this%indx = NOTANID            
            return
        end subroutine clear1
    
    !---
    
        pure function getSuper0(this) result(super)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      return simple supercell associtaed with link cell list
            type(LinkCell3D),intent(in)      ::      this
            type(SimpleSupercell)            ::      super
            super = this%super
            return
        end function getSuper0
    
        pure function getA0(this) result(a)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      return the 3x3 matrix describing extent of a link cell 
            type(LinkCell3D),intent(in)      ::      this
            real(kind=real64),dimension(3,3) ::      a
            a = getA(this%super)
            return
        end function geta0
        

        
        pure function getNx0(this) result(Nx)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      number of link cell repeats in x direction
            type(LinkCell3D),intent(in)     ::      this
            integer                         ::      Nx
            Nx = getNx(this%super)
            return
        end function getNx0


        pure function getNy0(this) result(Ny)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      number of link cell repeats in y direction
            type(LinkCell3D),intent(in)     ::      this
            integer                         ::      Ny
            Ny = getNy(this%super)
            return
        end function getNy0

        pure function getNz0(this) result(Nz)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      number of link cell repeats in z direction    
            type(LinkCell3D),intent(in)     ::      this
            integer                         ::      Nz
            Nz = getNz(this%super)
            return
        end function getNz0
        
    
        pure subroutine translate0(this,dx)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      translate all atoms by a small rigid vector dx
            type(LinkCell3D),intent(inout)                  ::      this
            real(kind=real64),dimension(3),intent(in)       ::      dx
            
            real(kind=real64),dimension(3)      ::      dy
            integer     ::      ix,iy,iz,ik
            dy = realSpaceToCell( this%super,dx )
            do iz = 0,this%super%Nz-1
                do iy = 0,this%super%Ny-1
                    do ix = 0,this%super%Nx-1
                        do ik = 1,this%c(ix,iy,iz)%nindx
                            this%c(ix,iy,iz)%x(1:3,ik) = this%c(ix,iy,iz)%x(1:3,ik) + dy(1:3)
                        end do
                    end do
                end do
            end do
            
            return
        end subroutine translate0
     
     
     
    
    !---        
          
        subroutine add1(this,x)
    !---^^^^^^^^^^^^^^^^^^^^^^^ 
    !*      add many points to real space positions x with default sequential id
            type(LinkCell3D),intent(inout)              ::      this
            real(kind=real64),dimension(:,:),intent(in) ::      x
            integer         ::      id
            do id = 1,size(x,dim=2)
                call add0(this,id,x(1:3,id))
            end do
            return
        end subroutine add1

        subroutine add0(this,id,x)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      add a point labelled id to real space position x
            type(LinkCell3D),intent(inout)      ::      this
            integer,intent(in)                  ::      id
            real(kind=real64),dimension(:),intent(in)   ::      x
            
            integer         ::      ix,iy,iz,nn
            real(kind=real64),dimension(3)      ::      yy,zz

        !---    find position in reduced cell coords
            yy = realSpaceToCell( this%super,x(1:3) )

        !---    find offset from cell origin in cell coords           
            zz(1:3) = yy(1:3) - floor(yy(1:3))
                                    
        !---    find cell
            call whichCell( this%super,yy,ix,iy,iz,cellSpace = .true. )
            
        !---    sanity check: too many points per cell?
            nn = this%c(ix,iy,iz)%nindx + 1                        
            
            if (nn > this%c(ix,iy,iz)%nMax) then
                !   problem. Need to reallocate.
                
                if (this%c(ix,iy,iz)%nMax == 16*this%nMax0) then
                    !   real problem. I've reallocated and reallocated and am still having trouble. There's a bug somewhere.
                    print *,"Lib_LinkCell3d::add0 error - have reallocated several times and am still finding too many points per cell"
                    print *,"    cell ",ix,iy,iz
                    print *,"    original guess points per cell ",this%nMax0," want ",2*this%c(ix,iy,iz)%nMax," in cell ",ix,iy,iz
                    print *,"    assuming a bug and crashing here. Sorry."
                    stop "Lib_LinkCell3d::add0 fatal error"
                end if
                
                if (LinkCell3D_dbg) then
                    print *,"Lib_LinkCell3d::add0 info - reallocating cell ",ix,iy,iz," to store ",this%c(ix,iy,iz)%nMax*2
                end if
                    
                call doubleCell(this%c(ix,iy,iz))                                 
                
            end if
            
            
        !---    add point to cell
        
        
            this%c(ix,iy,iz)%nindx = nn           
            this%c(ix,iy,iz)%indx(nn) = id
            this%c(ix,iy,iz)%x(1:3,nn) = zz(1:3)
            this%nMax = max(this%nMax,nn)
  
            if (LinkCell3D_dbg) then
                print *,"add0 id,x ",id,x," to cell ",ix,iy,iz," number ",nn
                call report(this%c(ix,iy,iz))
            end if
            
            return
        end subroutine add0
        

        subroutine cut0(this,id,x , ok)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      cut a point labelled id from real space position approximately x (should be correct cell)
            type(LinkCell3D),intent(inout)      ::      this
            integer,intent(in)                  ::      id
            real(kind=real64),dimension(:),intent(in)   ::      x
            logical,intent(out),optional        ::      ok
            
            
            integer         ::      ix,iy,iz,ik,nn
            real(kind=real64),dimension(3)      ::      yy!,zz

            if (present(ok)) ok = .true.
            
        !---    find position in reduced cell coords
            yy = realSpaceToCell( this%super,x(1:3) )

        !---    find cell
            call whichCell( this%super,yy,ix,iy,iz,cellSpace = .true. )
            
        !---    how many points in this cell?
            nn = this%c(ix,iy,iz)%nindx
             
        !---    which point is this?
            do ik = 1,nn
                if (this%c(ix,iy,iz)%indx(ik) == id) then
                !   swap last point in cell with this one.       
                    if (LinkCell3D_dbg) then
                        print *,"cut0 id,x ",id,x," from cell ",ix,iy,iz," number ",ik
                        call report(this%c(ix,iy,iz))
                    end if      
                    this%c(ix,iy,iz)%indx(ik) = this%c(ix,iy,iz)%indx(nn)
                    this%c(ix,iy,iz)%x(1:3,ik) = this%c(ix,iy,iz)%x(1:3,nn)                    
                    this%c(ix,iy,iz)%nindx = nn-1
                                   
                    return
                end if
            end do
            
            
            
            
        !---    problem - shouldn't have got this far.
            if (present(ok)) then
                !   it is possible, I suppose, that we meant a fail is possible???
                ok = .false.
                return
            else
                
                print *,"Lib_LinkCell3d::cut0 error - have searched cell and can't find atom with id ",id
                call report(this)
                print *,"    cell ",ix,iy,iz
                print *,"    x    ",x
                print *,"    y    ",yy
                call report( this%c(ix,iy,iz) )
                
                print *,"    assuming a bug and crashing here. Sorry."
                stop "Lib_LinkCell3d::cut0 fatal error"
            end if    
            return
        end subroutine cut0        
        
        
        subroutine move0(this,id,x , dx , ok)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      cut a point labelled id from real space position approximately x (should be correct cell)
    !*      replace at position x+dx
    !*      on exit ok = .true. if the atom moved within one cell and .false. if changed cell
            type(LinkCell3D),intent(inout)      ::      this
            integer,intent(in)                  ::      id
            real(kind=real64),dimension(:),intent(in)   ::      x
            real(kind=real64),dimension(3),intent(in)   ::      dx
            logical,intent(out)                         ::      ok
            
            
            integer                             ::      ix,iy,iz,ik ,jx,jy,jz  ,ni,nj , Nx,Ny,Nz
            real(kind=real64),dimension(3)      ::      yy 
            integer                             ::      f1,f2,f3

            
            real(kind=real64),dimension(3,3)    ::      ia_cell

        !---    unroll routines from supercells here for efficiency            
            ia_cell = getiA(this%super) 
            Nx = getNx(this%super)
            Ny = getNy(this%super)
            Nz = getNz(this%super)
            
            
        !---    find cell before move
            yy(1:3) = ia_cell(1:3,1)*x(1) + ia_cell(1:3,2)*x(2) + ia_cell(1:3,3)*x(3)
            ix = int( mod( floor( yy(1) ) + Nx*8,Nx ) ) 
            iy = int( mod( floor( yy(2) ) + Ny*8,Ny ) ) 
            iz = int( mod( floor( yy(3) ) + Nz*8,Nz ) )  
            
            
        !---    find cell after move
            yy(1:3) = yy(1:3) + ia_cell(1:3,1)*dx(1) + ia_cell(1:3,2)*dx(2) + ia_cell(1:3,3)*dx(3)
            call whichCell( this%super,yy,jx,jy,jz,cellSpace = .true. )
            f1 = floor( yy(1) )
            f2 = floor( yy(2) )
            f3 = floor( yy(3) )
            jx = int( mod( f1 + Nx*8,Nx ) )  !  note jx might not be equal to f1 due to periodic boundaries
            jy = int( mod( f2 + Ny*8,Ny ) )  
            jz = int( mod( f3 + Nz*8,Nz ) )  
            yy(1) = yy(1)-f1
            yy(2) = yy(2)-f2
            yy(3) = yy(3)-f3
            
            
            
            if ( (ix==jx).and.(iy==jy).and.(iz==jz) ) then             
                    
            !---    atom shifts position but does not change cell 
                do ik = 1,this%c(ix,iy,iz)%nindx
                    if (this%c(ix,iy,iz)%indx(ik) == id) then                         
                        this%c(ix,iy,iz)%x(1:3,ik) = yy(1:3)
                        exit
                    end if
                end do
                ok = .true.
                
            else
            
            !---    atom changes cell from (ix,iy,iz) to (jx,jy,jz)
                ni = this%c(ix,iy,iz)%nindx
                do ik = 1,ni             
                    if (this%c(ix,iy,iz)%indx(ik) == id) then      
                        !   this is the atom in cell (ix,iy,iz) to cut
                        this%c(ix,iy,iz)%indx(ik)  = this%c(ix,iy,iz)%indx(ni)
                        this%c(ix,iy,iz)%x(1:3,ik) = this%c(ix,iy,iz)%x(1:3,ni)                    
                        this%c(ix,iy,iz)%nindx     = ni-1
                        exit
                    end if
                end do
                
                nj = this%c(jx,jy,jz)%nindx + 1
                
                if (nj > this%c(jx,jy,jz)%nMax) then
                    !   problem. Need to reallocate.
                    
                    if (this%c(jx,jy,jz)%nMax == 16*this%nMax0) then
                        !   real problem. I've reallocated and reallocated and am still having trouble. There's a bug somewhere.
                        print *,"Lib_LinkCell3d::move0 error - have reallocated several times and am still finding too many points per cell"
                        print *,"    cell ",jx,jy,jz
                        print *,"    original guess points per cell ",this%nMax0," want ",2*this%c(jx,jy,jz)%nMax," in cell ",jx,jy,jz
                        print *,"    assuming a bug and crashing here. Sorry."
                        stop "Lib_LinkCell3d::move0 fatal error"
                    end if
                        
                    call doubleCell(this%c(jx,jy,jz)) 
                    
                               
                end if                
                
                this%c(jx,jy,jz)%indx(nj)  = id  
                this%c(jx,jy,jz)%x(1:3,nj) = yy(1:3)
                this%c(jx,jy,jz)%nindx     = nj
                this%nMax = max( this%nMax,nj )
                ok = .false.
                
            end if
            
             
            return
        end subroutine move0  
        
        
        pure function getPosition0( this,ix,iy,iz,ik,includeCellOffset ) result(x)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      returns the position of the ik th atom in the (ix,iy,iz) link cell.
    !*      optionally return the offset of the link cell origin to give a real space location.
            type(LinkCell3D),intent(in)                     ::      this
            integer,intent(in)                              ::      ix,iy,iz,ik
            real(kind=real64),dimension(3)                  ::      x
            logical,intent(in)                              ::      includeCellOffset
            x = this%c(ix,iy,iz)%x(1:3,ik)
            if (includeCellOffset) x = x + (/ix,iy,iz/)
            x = cellToRealSpace( this%super,x )
            return
        end function getPosition0
    

        pure function getId( this,ix,iy,iz,ik ) result(id)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      returns id of ik th atom in the (ix,iy,iz) link cell.
            type(LinkCell3D),intent(in)                     ::      this
            integer,intent(in)                              ::      ix,iy,iz,ik
            integer                                         ::      id
            id = this%c(ix,iy,iz)%indx(ik)
            return
        end function getId        
        
    !---
    
    
        subroutine nearestNeighbour0( this,x, id,dr2min,self,dx )
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      return the nearest neighbour about a point in real space x  
            type(LinkCell3D),intent(in)                     ::      this
            real(kind=real64),dimension(3),intent(in)       ::      x
            logical,intent(in)                              ::      self    !   should we count a point exactly at x?
            
            integer,intent(out)                             ::      id
            real(kind=real64),intent(out)                   ::      dr2min       
            real(kind=real64),dimension(3),intent(out),optional     ::      dx      
            integer                             ::      ix,iy,iz
            integer                             ::      jx,jy,jz , jk
            integer                             ::      kx,ky,kz 
            real(kind=real64),dimension(3)      ::      yy,zi,dzi,dzij
            real(kind=real64)                   ::      dd,ddcutoff
            
            
            integer                             ::      Nx,Ny,Nz
            real(kind=real64),dimension(3,3)    ::      aa 
            
            
            
        !---    find position in reduced cell coords
            aa = getA(this%super)
            Nx = getNx(this%super)
            Ny = getNy(this%super)
            Nz = getNz(this%super)
              
            if (self) then
                ddcutoff = -1
            else    
                ddcutoff = 1.0d-8
            end if
            id = 0
            dr2min = huge(1.0)
            
            
        !---    find position in reduced cell coords
            yy = realSpaceToCell( this%super,x )          

        !---    find offset from cell origin in cell coords
            zi(1:3) = yy(1:3) - floor(yy(1:3))
                                    
        !---    find cell
            call whichCell( this%super,yy,ix,iy,iz,cellSpace = .true. )
            
            if (LinkCell3D_dbg) then
                print *,"nneigh ",x," i ",ix,iy,iz
                call report(this)
                call report(this%c(ix,iy,iz))
            end if
            
            
         !---    find neighbour cells
            do kz = -1,1
                jz = mod( iz + Nz + kz , Nz )
                dzi(3) = kz - zi(3)
                do ky = -1,1
                    jy = mod( iy + Ny + ky , Ny )
                    dzi(2) = ky - zi(2)
                    do kx = -1,1
                        jx = mod( ix + Nx + kx , Nx )
                        dzi(1) = kx - zi(1)
                                                                           
                        do jk = 1,this%c(jx,jy,jz)%nindx
                        
                        
                        !---    find separation between points in cell space
                            dzij(1:3) = this%c(jx,jy,jz)%x(1:3,jk) + dzi(1:3)
            
                        !---    find separation between points in real space
                            !dzij = cellToRealSpace( this%super,dzij )
                            dzij(1:3) = aa(1:3,1)*dzij(1) + aa(1:3,2)*dzij(2) + aa(1:3,3)*dzij(3)
    
                        !---    find square distance
                            dd = dzij(1)*dzij(1) + dzij(2)*dzij(2) + dzij(3)*dzij(3) 
    
!                             if (LinkCell3D_dbg) then
!                                 print *," j ",jx,jy,jz,jk," at ",this%c(jx,jy,jz)%x(1:3,jk)
!                                 print *," dx ",dzij,dd
!                             end if
                            
                            if ((dd>ddcutoff).and.(dd<dr2min)) then
                                dr2min = dd
                                id = this%c(jx,jy,jz)%indx(jk)
                                if (present(dx)) dx = dzij
                            end if
                            
                        end do
                    end do 
                end do
            end do
            return
        end subroutine nearestNeighbour0     
        
        
        
        subroutine nearestNeighbour1( this,x, id,dr2min )
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ 
    !*      return the nearest neighbour about a point in real space x  
            type(LinkCell3D),intent(in)                     ::      this
            real(kind=real64),dimension(3),intent(in)       ::      x
            
            integer,intent(out)                             ::      id
            real(kind=real64),intent(out)                   ::      dr2min       
            integer                             ::      ix,iy,iz
            integer                             ::      jx,jy,jz , jk
            integer                             ::      kx,ky,kz 
            real(kind=real64),dimension(3)      ::      yy,zi,dzi,dzij
            real(kind=real64)                   ::      dd 
            
      
            integer                             ::      Nx,Ny,Nz
            real(kind=real64),dimension(3,3)    ::      aa 
            
            id = 0
            dr2min = huge(1.0)
            aa = getA(this%super)
            Nx = getNx(this%super)
            Ny = getNy(this%super)
            Nz = getNz(this%super)
            
        !---    find position in reduced cell coords
            yy = realSpaceToCell( this%super,x )          

        !---    find offset from cell origin in cell coords
            zi(1:3) = yy(1:3) - floor(yy(1:3))
                                    
        !---    find cell
            call whichCell( this%super,yy,ix,iy,iz,cellSpace = .true. )
            
            if (LinkCell3D_dbg) then
                print *,"nneigh ",x," i ",ix,iy,iz
                call report(this)
                call report(this%c(ix,iy,iz))
            end if
            
            
         !---    find neighbour cells
            do kz = -1,1
                jz = mod( iz + Nz + kz , Nz )
                dzi(3) = kz - zi(3)
                do ky = -1,1
                    jy = mod( iy + Ny + ky , Ny )
                    dzi(2) = ky - zi(2)
                    do kx = -1,1
                        jx = mod( ix + Nx + kx , Nx )
                        dzi(1) = kx - zi(1)
                                                                           
                        do jk = 1,this%c(jx,jy,jz)%nindx
                        
                        
                        !---    find separation between points in cell space
                            dzij(1:3) = this%c(jx,jy,jz)%x(1:3,jk) + dzi(1:3)
            
                        !---    find separation between points in real space
                            !dzij = cellToRealSpace( this%super,dzij )
                            dzij(1:3) = aa(1:3,1)*dzij(1) + aa(1:3,2)*dzij(2) + aa(1:3,3)*dzij(3)
    
                        !---    find square distance
                            dd = dzij(1)*dzij(1) + dzij(2)*dzij(2) + dzij(3)*dzij(3) 
    
                            if (LinkCell3D_dbg) then
                                print *," j ",jx,jy,jz,jk," at ",this%c(jx,jy,jz)%x(1:3,jk)
                                print *," dx ",dzij,dd
                            end if
                            
                            if (dd<dr2min) then
                                dr2min = dd
                                id = this%c(jx,jy,jz)%indx(jk)
                                 
                            end if
                            
                        end do
                    end do 
                end do
            end do
            return
        end subroutine nearestNeighbour1     
        
        
    !---
    
    
        subroutine neighbourList0( this,x,sigma, n,id,dx )
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      return the neighbours about a point in real space x within sphere radius sigma
            type(LinkCell3D),intent(in)                     ::      this
            real(kind=real64),dimension(3),intent(in)       ::      x
            real(kind=real64),intent(in)                    ::      sigma
            integer,intent(out)                             ::      n
            integer,dimension(:),intent(out)                ::      id
            real(kind=real64),dimension(:,:),intent(out)    ::      dx
            
            integer                             ::      ix,iy,iz
            integer                             ::      jx,jy,jz , jk
            integer                             ::      kx,ky,kz 
            real(kind=real64),dimension(3)      ::      yy,zi  
            real(kind=real64)                   ::      dd,dzi1,dzi2,dzi3 ,dzij1,dzij2,dzij3 ,dzj1,dzj2,dzj3
            real(kind=real64)                   ::      a11,a21,a31,a12,a22,a32,a13,a23,a33
            n = 0
            
            
        !---    find position in reduced cell coords
            yy = realSpaceToCell( this%super,x )
            a11 = this%super%a(1,1)
            a21 = this%super%a(2,1)
            a31 = this%super%a(3,1)
            a12 = this%super%a(1,2)
            a22 = this%super%a(2,2)
            a32 = this%super%a(3,2)
            a13 = this%super%a(1,3)
            a23 = this%super%a(2,3)
            a33 = this%super%a(3,3) 

        !---    find offset from cell origin in cell coords
            zi(1:3) = yy(1:3) - floor(yy(1:3))
                                    
        !---    find cell
            call whichCell( this%super,yy,ix,iy,iz,cellSpace = .true. )

        !---    find neighbour cells
            do kz = -1,1
                jz = mod( iz + this%super%Nz + kz , this%super%Nz )
                dzi3 = kz - zi(3)
                
                do ky = -1,1
                    jy = mod( iy + this%super%Ny + ky , this%super%Ny )
                    dzi2 = ky - zi(2)
                    
                    do kx = -1,1
                        jx = mod( ix + this%super%Nx + kx , this%super%Nx )
                        dzi1 = kx - zi(1)
                                                                                                   
                        do jk = 1,this%c(jx,jy,jz)%nindx
                        
                        
                        !---    find separation between points in cell space
                            dzj1 = this%c(jx,jy,jz)%x(1,jk) + dzi1
                            dzj2 = this%c(jx,jy,jz)%x(2,jk) + dzi2
                            dzj3 = this%c(jx,jy,jz)%x(3,jk) + dzi3
            
                        !---    find separation between points in real space
                            dzij1 = a11*dzj1 + a12*dzj2 + a13*dzj3
                            dzij2 = a21*dzj1 + a22*dzj2 + a23*dzj3
                            dzij3 = a31*dzj1 + a32*dzj2 + a33*dzj3
                    
                        !---    find square distance
                            dd = dzij1*dzij1 + dzij2*dzij2 + dzij3*dzij3 
                                                      
                        !---    compare to cutoff
                            if (dd <= sigma*sigma) then
                                n = n + 1
                                id(n) = this%c(jx,jy,jz)%indx(jk)   
                                dx(1,n) = dzij1                           
                                dx(2,n) = dzij2                            
                                dx(3,n) = dzij3                           
                            end if
                        
                        end do
                    end do
                end do
            end do
            return
        end subroutine neighbourList0
     
    
        subroutine neighbourList0a( this,x,sigma, n,id,dx,dr2 )
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      return the neighbours about a point in real space x within sphere radius sigma
            type(LinkCell3D),intent(in)                     ::      this
            real(kind=real64),dimension(3),intent(in)       ::      x
            real(kind=real64),intent(in)                    ::      sigma
            integer,intent(out)                             ::      n
            integer,dimension(:),intent(out)                ::      id
            real(kind=real64),dimension(:,:),intent(out)    ::      dx
            real(kind=real64),dimension(:),intent(out)      ::      dr2
            
            integer                             ::      Nx,Ny,Nz
            integer                             ::      ix,iy,iz
            integer                             ::      jx,jy,jz , jk
            integer                             ::      kx,ky,kz 
            real(kind=real64),dimension(3)      ::      yy,zi,dzi,dzij
            real(kind=real64)                   ::      dd
            real(kind=real64),dimension(3,3)    ::      aa 
            
            n = 0
            
            
        !---    find position in reduced cell coords
            aa = getA(this%super)
            Nx = getNx(this%super)
            Ny = getNy(this%super)
            Nz = getNz(this%super)
            yy = realSpaceToCell( this%super,x )
            !yy(1:3) = this%ia(1:3,1)*x(1) + this%ia(1:3,2)*x(2) + this%ia(1:3,3)*x(3)   

        !---    find offset from cell origin in cell coords
            zi(1:3) = yy(1:3) - floor(yy(1:3))
                                    
        !---    find cell
            !call whichCell( this%super,yy,ix,iy,iz,cellSpace = .true. ) 
            ix = mod( int( yy(1) + 256*Nx ), Nx )     !    pbc
            iy = mod( int( yy(2) + 256*Ny ), Ny )
            iz = mod( int( yy(3) + 256*Nz ), Nz )

        !---    find neighbour cells
            do kz = -1,1
                jz = mod( iz + Nz + kz , Nz )
                do ky = -1,1
                    jy = mod( iy + Ny + ky , Ny )
                    do kx = -1,1
                        jx = mod( ix + Nx + kx , Nx )
                        dzi(1:3) = (/ kx,ky,kz /) - zi(1:3)
                                                                           
                        do jk = 1,this%c(jx,jy,jz)%nindx
                        
                        
                        !---    find separation between points in cell space
                            dzij(1:3) = this%c(jx,jy,jz)%x(1:3,jk) + dzi(1:3)
            
                        !---    find separation between points in real space
                            ! dzij = cellToRealSpace( this%super,dzij )
                              dzij(1:3) = aa(1:3,1)*dzij(1) + aa(1:3,2)*dzij(2) + aa(1:3,3)*dzij(3)
    
                        !---    find square distance
                            dd = dzij(1)*dzij(1) + dzij(2)*dzij(2) + dzij(3)*dzij(3) 
    
                        !---    compare to cutoff
                            if (dd <= sigma*sigma) then
                                n = n + 1
                                id(n) = this%c(jx,jy,jz)%indx(jk)
                                dx(1:3,n) = dzij(1:3)
                                dr2(n) = dd
                            end if
                        end do
                    end do
                end do
            end do
            return
        end subroutine neighbourList0a 
    

        subroutine neighbourList0c_27( this,x, n,id,dx )
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      return all neighbours about a point in real space x, imposing minimum image convention  
            type(LinkCell3D),intent(in)                     ::      this
            real(kind=real64),dimension(3),intent(in)       ::      x
            integer,intent(out)                             ::      n
            integer,dimension(:),intent(out)                ::      id
            real(kind=real64),dimension(:,:),intent(out)    ::      dx
            
           ! integer                             ::      ix,iy,iz
            integer                             ::      jx,jy,jz , jk
           ! integer                             ::      kx,ky,kz 
            real(kind=real64),dimension(3)      ::      yy!,zi!,dzi
            real(kind=real64)                   ::      dzi1,dzi2,dzi3 ,dzij1,dzij2,dzij3 ,dzj1,dzj2,dzj3      !,dd
            real(kind=real64)                   ::      a11,a21,a31,a12,a22,a32,a13,a23,a33
            
            n = 0
            
             
        !---    find position in reduced cell coords
            yy = realSpaceToCell( this%super,x )
            a11 = this%super%a(1,1)
            a21 = this%super%a(2,1)
            a31 = this%super%a(3,1)
            a12 = this%super%a(1,2)
            a22 = this%super%a(2,2)
            a32 = this%super%a(3,2)
            a13 = this%super%a(1,3)
            a23 = this%super%a(2,3)
            a33 = this%super%a(3,3) 
 
        !---    find neighbour cells
            do jz = 0,this%super%Nz-1               
                dzi3 = jz - yy(3) 
                do jy = 0,this%super%Ny-1 
                    dzi2 = jy - yy(2)
                    do jx = 0,this%super%Nx-1     
                        dzi1 = jx - yy(1)                                                                        
                        
                        do jk = 1,this%c(jx,jy,jz)%nindx
                        
                        
                        !---    find separation between points in cell space
                            dzj1 = this%c(jx,jy,jz)%x(1,jk) + dzi1
                            dzj2 = this%c(jx,jy,jz)%x(2,jk) + dzi2
                            dzj3 = this%c(jx,jy,jz)%x(3,jk) + dzi3
                            
                        !---    minimum image convention
                            if ( 2*dzj1 > this%super%Nx ) then
                                dzj1 = dzj1 - this%super%Nx
                            else if ( 2*dzj1 <= -this%super%Nx ) then
                                dzj1 = dzj1 + this%super%Nx
                            end if
                            if ( 2*dzj2 > this%super%Ny ) then
                                dzj2 = dzj2 - this%super%Ny
                            else if ( 2*dzj2 <= -this%super%Ny ) then
                                dzj2 = dzj2 + this%super%Ny
                            end if
                            if ( 2*dzj3 > this%super%Nz ) then
                                dzj3 = dzj3 - this%super%Nz
                            else if ( 2*dzj3 <= -this%super%Nz ) then
                                dzj3 = dzj3 + this%super%Nz
                            end if
                             
                            
                            
                        !---    find separation between points in real space
                            dzij1 = a11*dzj1 + a12*dzj2 + a13*dzj3
                            dzij2 = a21*dzj1 + a22*dzj2 + a23*dzj3
                            dzij3 = a31*dzj1 + a32*dzj2 + a33*dzj3
     
                            n = n + 1
                            id(n) = this%c(jx,jy,jz)%indx(jk)
                            dx(1:3,n) = (/dzij1,dzij2,dzij3/)
                        
                        end do
                    end do
                end do
            end do
            return
        end subroutine neighbourList0c_27        
    
        subroutine neighbourList1( this,x,sigma, n,dx )
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      return the neighbours about a point in real space x within sphere radius sigma
            type(LinkCell3D),intent(in)                     ::      this
            real(kind=real64),dimension(3),intent(in)       ::      x
            real(kind=real64),intent(in)                    ::      sigma
            integer,intent(out)                             ::      n
            real(kind=real64),dimension(:,:),intent(out)    ::      dx
            
            integer                             ::      ix,iy,iz
            integer                             ::      jx,jy,jz , jk
            integer                             ::      kx,ky,kz 
            real(kind=real64),dimension(3)      ::      yy,zi!,dzi
            real(kind=real64)                   ::      dd,dzi1,dzi2,dzi3 ,dzij1,dzij2,dzij3 ,dzj1,dzj2,dzj3
            real(kind=real64)                   ::      a11,a21,a31,a12,a22,a32,a13,a23,a33
            n = 0
            
            
        !---    find position in reduced cell coords
            yy = realSpaceToCell( this%super,x )
            a11 = this%super%a(1,1)
            a21 = this%super%a(2,1)
            a31 = this%super%a(3,1)
            a12 = this%super%a(1,2)
            a22 = this%super%a(2,2)
            a32 = this%super%a(3,2)
            a13 = this%super%a(1,3)
            a23 = this%super%a(2,3)
            a33 = this%super%a(3,3) 

        !---    find offset from cell origin in cell coords
            zi(1:3) = yy(1:3) - floor(yy(1:3))
                                    
        !---    find cell
            call whichCell( this%super,yy,ix,iy,iz,cellSpace = .true. )

        !---    find neighbour cells
            do kz = -1,1
                jz = mod( iz + this%super%Nz + kz , this%super%Nz )
                dzi3 = kz - zi(3)
                
                do ky = -1,1
                    jy = mod( iy + this%super%Ny + ky , this%super%Ny )
                    dzi2 = ky - zi(2)
                    
                    do kx = -1,1
                        jx = mod( ix + this%super%Nx + kx , this%super%Nx )
                        dzi1 = kx - zi(1)
                                                                                                   
                        do jk = 1,this%c(jx,jy,jz)%nindx
                        
                        
                        !---    find separation between points in cell space
                            dzj1 = this%c(jx,jy,jz)%x(1,jk) + dzi1
                            dzj2 = this%c(jx,jy,jz)%x(2,jk) + dzi2
                            dzj3 = this%c(jx,jy,jz)%x(3,jk) + dzi3
            
                        !---    find separation between points in real space
                            dzij1 = a11*dzj1 + a12*dzj2 + a13*dzj3
                            dzij2 = a21*dzj1 + a22*dzj2 + a23*dzj3
                            dzij3 = a31*dzj1 + a32*dzj2 + a33*dzj3
    
                        !---    find square distance
                            dd = dzij1*dzij1 + dzij2*dzij2 + dzij3*dzij3 
    
                        !---    compare to cutoff
                            if (dd <= sigma*sigma) then
                                n = n + 1
                                dx(1,n) = dzij1
                                dx(2,n) = dzij2
                                dx(3,n) = dzij3
                            end if
                        
                        end do
                    end do
                end do
            end do
            return
        end subroutine neighbourList1 
        
        subroutine neighbourList1a( this,x,sigma, n,r2 )
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      return distance squared to the neighbours about a point in real space x within sphere radius sigma
            type(LinkCell3D),intent(in)                     ::      this
            real(kind=real64),dimension(3),intent(in)       ::      x
            real(kind=real64),intent(in)                    ::      sigma
            integer,intent(out)                             ::      n
            real(kind=real64),dimension(:),intent(out)      ::      r2
            
            integer                             ::      ix,iy,iz
            integer                             ::      jx,jy,jz , jk
            integer                             ::      kx,ky,kz 
            real(kind=real64),dimension(3)      ::      yy,zi,dzi,dzij
            real(kind=real64)                   ::      dd
            
            n = 0
            
            
        !---    find position in reduced cell coords
            yy = realSpaceToCell( this%super,x )

        !---    find offset from cell origin in cell coords
            zi(1:3) = yy(1:3) - floor(yy(1:3))
                                    
        !---    find cell
            call whichCell( this%super,yy,ix,iy,iz,cellSpace = .true. )

        !---    find neighbour cells
            do kz = -1,1
                jz = mod( iz + this%super%Nz + kz , this%super%Nz )
                do ky = -1,1
                    jy = mod( iy + this%super%Ny + ky , this%super%Ny )
                    do kx = -1,1
                        jx = mod( ix + this%super%Nx + kx , this%super%Nx )
                        dzi(1:3) = (/ kx,ky,kz /) - zi(1:3)
                                                                           
                        do jk = 1,this%c(jx,jy,jz)%nindx
                        
                        
                        !---    find separation between points in cell space
                            dzij(1:3) = this%c(jx,jy,jz)%x(1:3,jk) + dzi(1:3)
            
                        !---    find separation between points in real space
                            dzij = cellToRealSpace( this%super,dzij )
    
                        !---    find square distance
                            dd = dzij(1)*dzij(1) + dzij(2)*dzij(2) + dzij(3)*dzij(3) 
    
                        !---    compare to cutoff
                            if (dd <= sigma*sigma) then
                                n = n + 1
                                r2(n) = dd
                            end if
                        end do
                    end do
                end do
            end do
            return
        end subroutine neighbourList1a 
        

        
        
        subroutine neighbourList1b( this,x,sigma, n,id,r2 )
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      return the neighbours about a point in real space x within sphere radius sigma
            type(LinkCell3D),intent(in)                     ::      this
            real(kind=real64),dimension(3),intent(in)       ::      x
            real(kind=real64),intent(in)                    ::      sigma
            integer,intent(out)                             ::      n
            integer,dimension(:),intent(out)                ::      id
            real(kind=real64),dimension(:),intent(out)      ::      r2
            
            integer                             ::      Nx,Ny,Nz
            integer                             ::      ix,iy,iz
            integer                             ::      jx,jy,jz , jk
            integer                             ::      kx,ky,kz 
            real(kind=real64),dimension(3)      ::      zi,dzij   
            real(kind=real64)                   ::      dd,dzi1,dzi2,dzi3  ,dzj1,dzj2,dzj3
            
            real(kind=real64),dimension(3,3)    ::      a_cell 
            
            n = 0
            
            Nx = this%super%Nx
            Ny = this%super%Ny
            Nz = this%super%Nz
            
            a_cell = getA(this%super)
            
            
        !---    find position in reduced cell coords
        
             zi = realSpaceToCell( this%super,x )
        

        !---    find offset from cell origin in cell coords
            
             
                                    
        !---    find cell
        !    call whichCell( this%super,yy,ix,iy,iz,cellSpace = .true. )
            ix = mod( floor(zi(1)) + 16*Nx,Nx )  
            iy = mod( floor(zi(2)) + 16*Ny,Ny )  
            iz = mod( floor(zi(3)) + 16*Nz,Nz )  
        
        
            zi(1) = zi(1) - ix
            zi(2) = zi(2) - iy
            zi(3) = zi(3) - iz
        

        !---    find neighbour cells
            do kz = -1,1
                jz = mod( iz + Nz + kz , Nz )
                dzi3 = kz - zi(3)
                
                do ky = -1,1
                    jy = mod( iy + Ny + ky , Ny )
                    dzi2 = ky - zi(2)
                    
                    do kx = -1,1
                        jx = mod( ix + Nx + kx , Nx )
                        dzi1 = kx - zi(1)
                                                                                                   
                        do jk = 1,this%c(jx,jy,jz)%nindx
                        
                        
                        !---    find separation between points in cell space
                            dzj1 = this%c(jx,jy,jz)%x(1,jk) + dzi1
                            dzj2 = this%c(jx,jy,jz)%x(2,jk) + dzi2
                            dzj3 = this%c(jx,jy,jz)%x(3,jk) + dzi3
            
                        !---    find separation between points in real space
                            dzij(1:3) = a_cell(1:3,1)*dzj1 + a_cell(1:3,2)*dzj2 + a_cell(1:3,3)*dzj3
                            !dzij2 = a21*dzj1 + a22*dzj2 + a23*dzj3
                            !dzij3 = a31*dzj1 + a32*dzj2 + a33*dzj3
    
                        !---    find square distance
                            dd = dzij(1)*dzij(1) + dzij(2)*dzij(2) + dzij(3)*dzij(3) 
    
                        !---    compare to cutoff
                            if (dd <= sigma*sigma) then
                                n = n + 1
                                id(n) = this%c(jx,jy,jz)%indx(jk)
                                r2(n) = dd
                            end if
                        
                        end do
                    end do
                end do
            end do
            return
        end subroutine neighbourList1b 
        
        subroutine neighbourList1d( this,x,sigma, n,id )
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      return the neighbours about a point in real space x within sphere radius sigma
            type(LinkCell3D),intent(in)                     ::      this
            real(kind=real64),dimension(3),intent(in)       ::      x
            real(kind=real64),intent(in)                    ::      sigma
            integer,intent(out)                             ::      n
            integer,dimension(:),intent(out)                ::      id
            
            integer                             ::      ix,iy,iz
            integer                             ::      jx,jy,jz , jk
            integer                             ::      kx,ky,kz 
            real(kind=real64),dimension(3)      ::      yy,zi  
            real(kind=real64)                   ::      dd,dzi1,dzi2,dzi3 ,dzij1,dzij2,dzij3 ,dzj1,dzj2,dzj3
            real(kind=real64)                   ::      a11,a21,a31,a12,a22,a32,a13,a23,a33
            n = 0
            
            
        !---    find position in reduced cell coords
            yy = realSpaceToCell( this%super,x )
            a11 = this%super%a(1,1)
            a21 = this%super%a(2,1)
            a31 = this%super%a(3,1)
            a12 = this%super%a(1,2)
            a22 = this%super%a(2,2)
            a32 = this%super%a(3,2)
            a13 = this%super%a(1,3)
            a23 = this%super%a(2,3)
            a33 = this%super%a(3,3) 

        !---    find offset from cell origin in cell coords
            zi(1:3) = yy(1:3) - floor(yy(1:3))
                                    
        !---    find cell
            call whichCell( this%super,yy,ix,iy,iz,cellSpace = .true. )

        !---    find neighbour cells
            do kz = -1,1
                jz = mod( iz + this%super%Nz + kz , this%super%Nz )
                dzi3 = kz - zi(3)
                
                do ky = -1,1
                    jy = mod( iy + this%super%Ny + ky , this%super%Ny )
                    dzi2 = ky - zi(2)
                    
                    do kx = -1,1
                        jx = mod( ix + this%super%Nx + kx , this%super%Nx )
                        dzi1 = kx - zi(1)
                                                                                                   
                        do jk = 1,this%c(jx,jy,jz)%nindx
                        
                        
                        !---    find separation between points in cell space
                            dzj1 = this%c(jx,jy,jz)%x(1,jk) + dzi1
                            dzj2 = this%c(jx,jy,jz)%x(2,jk) + dzi2
                            dzj3 = this%c(jx,jy,jz)%x(3,jk) + dzi3
            
                        !---    find separation between points in real space
                            dzij1 = a11*dzj1 + a12*dzj2 + a13*dzj3
                            dzij2 = a21*dzj1 + a22*dzj2 + a23*dzj3
                            dzij3 = a31*dzj1 + a32*dzj2 + a33*dzj3
                    
                        !---    find square distance
                            dd = dzij1*dzij1 + dzij2*dzij2 + dzij3*dzij3 
                            
                           
    
    
                        !---    compare to cutoff
                            if (dd <= sigma*sigma) then
                                n = n + 1
                                id(n) = this%c(jx,jy,jz)%indx(jk)                                
                            end if
                        
                        end do
                    end do
                end do
            end do
            return
        end subroutine neighbourList1d 
        

        subroutine neighbourList1e( this,x,sigma, n )
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ 
    !*      return the number of neighbours about a point in real space x within sphere radius sigma
            type(LinkCell3D),intent(in)                     ::      this
            real(kind=real64),dimension(3),intent(in)       ::      x
            real(kind=real64),intent(in)                    ::      sigma
            integer,intent(out)                             ::      n
            
            integer                             ::      ix,iy,iz
            integer                             ::      jx,jy,jz , jk
            integer                             ::      kx,ky,kz 
            real(kind=real64),dimension(3)      ::      yy,zi!,dzi
            real(kind=real64)                   ::      dd,dzi1,dzi2,dzi3 ,dzij1,dzij2,dzij3 ,dzj1,dzj2,dzj3
            real(kind=real64)                   ::      a11,a21,a31,a12,a22,a32,a13,a23,a33
            n = 0
            
            
        !---    find position in reduced cell coords
            yy = realSpaceToCell( this%super,x )
            a11 = this%super%a(1,1)
            a21 = this%super%a(2,1)
            a31 = this%super%a(3,1)
            a12 = this%super%a(1,2)
            a22 = this%super%a(2,2)
            a32 = this%super%a(3,2)
            a13 = this%super%a(1,3)
            a23 = this%super%a(2,3)
            a33 = this%super%a(3,3) 

        !---    find offset from cell origin in cell coords
            zi(1:3) = yy(1:3) - floor(yy(1:3))
                                    
        !---    find cell
            call whichCell( this%super,yy,ix,iy,iz,cellSpace = .true. )

        !---    find neighbour cells
            do kz = -1,1
                jz = mod( iz + this%super%Nz + kz , this%super%Nz )
                dzi3 = kz - zi(3)
                
                do ky = -1,1
                    jy = mod( iy + this%super%Ny + ky , this%super%Ny )
                    dzi2 = ky - zi(2)
                    
                    do kx = -1,1
                        jx = mod( ix + this%super%Nx + kx , this%super%Nx )
                        dzi1 = kx - zi(1)
                                                                                                   
                        do jk = 1,this%c(jx,jy,jz)%nindx
                        
                        
                        !---    find separation between points in cell space
                            dzj1 = this%c(jx,jy,jz)%x(1,jk) + dzi1
                            dzj2 = this%c(jx,jy,jz)%x(2,jk) + dzi2
                            dzj3 = this%c(jx,jy,jz)%x(3,jk) + dzi3
            
                        !---    find separation between points in real space
                            dzij1 = a11*dzj1 + a12*dzj2 + a13*dzj3
                            dzij2 = a21*dzj1 + a22*dzj2 + a23*dzj3
                            dzij3 = a31*dzj1 + a32*dzj2 + a33*dzj3
    
                        !---    find square distance
                            dd = dzij1*dzij1 + dzij2*dzij2 + dzij3*dzij3 
    
                        !---    compare to cutoff
                            if (dd <= sigma*sigma)  n = n + 1
                        
                        end do
                    end do
                end do
            end do
            return
        end subroutine neighbourList1e 
                
        subroutine neighbourList3c( this,ix,iy,iz,ik,sigma, n,id )
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ 
    !*      return the neighbours about a point 
            type(LinkCell3D),intent(in)                     ::      this
            integer,intent(in)                              ::      ik,ix,iy,iz  
            real(kind=real64),intent(in)                    ::      sigma
            integer,intent(out)                             ::      n
            integer,dimension(:),intent(out)                ::      id
 
            
            
            integer                             ::      jx,jy,jz , jk
            integer                             ::      kx,ky,kz 
            real(kind=real64),dimension(3)      ::      zi!,yy,dzi
            real(kind=real64)                   ::      dzi1,dzi2,dzi3 ,dzij1,dzij2,dzij3 ,dzj1,dzj2,dzj3
            real(kind=real64)                   ::      a11,a21,a31,a12,a22,a32,a13,a23,a33 , dd
             
            
            n = 0
            
        !---    find position in reduced cell coords
             
            a11 = this%super%a(1,1)
            a21 = this%super%a(2,1)
            a31 = this%super%a(3,1)
            a12 = this%super%a(1,2)
            a22 = this%super%a(2,2)
            a32 = this%super%a(3,2)
            a13 = this%super%a(1,3)
            a23 = this%super%a(2,3)
            a33 = this%super%a(3,3) 

        !---    find offset from cell origin in cell coords
            zi(1:3) = this%c(ix,iy,iz)%x(1:3,ik)
            
            if (LinkCell3D_dbg) print *,"lc3d ",ix,iy,iz,ik,zi
                         
        !---    find neighbour cells
            do kz = -1,1
                jz = mod( iz + this%super%Nz + kz , this%super%Nz )
                dzi3 = kz - zi(3)
                
                do ky = -1,1
                    jy = mod( iy + this%super%Ny + ky , this%super%Ny )
                    dzi2 = ky - zi(2)
                    
                    do kx = -1,1
                        jx = mod( ix + this%super%Nx + kx , this%super%Nx )
                        dzi1 = kx - zi(1)
                                                                                                   
                        do jk = 1,this%c(jx,jy,jz)%nindx
                            
                            if ( (kx*kx + ky*ky + kz*kz == 0).and.(ik==jk) ) cycle      !   do not add self to neighbour list
                        
                        !---    find separation between points in cell space
                            dzj1 = this%c(jx,jy,jz)%x(1,jk) + dzi1
                            dzj2 = this%c(jx,jy,jz)%x(2,jk) + dzi2
                            dzj3 = this%c(jx,jy,jz)%x(3,jk) + dzi3
                            
                            
            
                        !---    find separation between points in real space
                            dzij1 = a11*dzj1 + a12*dzj2 + a13*dzj3
                            dzij2 = a21*dzj1 + a22*dzj2 + a23*dzj3
                            dzij3 = a31*dzj1 + a32*dzj2 + a33*dzj3
    
                        !---    find square distance
                            
                            dd = dzij1*dzij1 + dzij2*dzij2 + dzij3*dzij3 
                            
                            
                            
                            if (dd<=sigma*sigma) then
                                n = n + 1
                                id(n) = this%c(jx,jy,jz)%indx(jk)
                                if (LinkCell3D_dbg) print *,"lc3d ",jx,jy,jz,"neigh ",n,id(n),dzij1,dzij2,dzij3,dd
                            end if
                        end do
                    end do
                end do
            end do
            return
        end subroutine neighbourList3c 
 
        
        subroutine neighbourList1c( this,x, n,id,r2 )
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ 
    !*      return the neighbours about a point in real space x  in all 27 cells, no cutoff
            type(LinkCell3D),intent(in)                     ::      this
            real(kind=real64),dimension(3),intent(in)       ::      x            
            integer,intent(out)                             ::      n
            integer,dimension(:),intent(out)                ::      id
            real(kind=real64),dimension(:),intent(out)      ::      r2
            
            integer                             ::      ix,iy,iz
            integer                             ::      jx,jy,jz , jk
            integer                             ::      kx,ky,kz 
            real(kind=real64),dimension(3)      ::      yy,zi!,dzi
            real(kind=real64)                   ::      dzi1,dzi2,dzi3 ,dzij1,dzij2,dzij3 ,dzj1,dzj2,dzj3
            real(kind=real64)                   ::      a11,a21,a31,a12,a22,a32,a13,a23,a33
            
            if (getNnodes(this%super)==27) then
                call neighbourList1c_27( this,x, n,id,r2 )
                return
            end if
            
            n = 0
            
            
            
            
        !---    find position in reduced cell coords
            yy = realSpaceToCell( this%super,x )
            a11 = this%super%a(1,1)
            a21 = this%super%a(2,1)
            a31 = this%super%a(3,1)
            a12 = this%super%a(1,2)
            a22 = this%super%a(2,2)
            a32 = this%super%a(3,2)
            a13 = this%super%a(1,3)
            a23 = this%super%a(2,3)
            a33 = this%super%a(3,3) 

        !---    find offset from cell origin in cell coords
            zi(1:3) = yy(1:3) - floor(yy(1:3))
                                    
        !---    find cell
            call whichCell( this%super,yy,ix,iy,iz,cellSpace = .true. )

        !---    find neighbour cells
            do kz = -1,1
                jz = mod( iz + this%super%Nz + kz , this%super%Nz )
                dzi3 = kz - zi(3)
                
                do ky = -1,1
                    jy = mod( iy + this%super%Ny + ky , this%super%Ny )
                    dzi2 = ky - zi(2)
                    
                    do kx = -1,1
                        jx = mod( ix + this%super%Nx + kx , this%super%Nx )
                        dzi1 = kx - zi(1)
                                                                                                   
                        do jk = 1,this%c(jx,jy,jz)%nindx
                        
                        
                        !---    find separation between points in cell space
                            dzj1 = this%c(jx,jy,jz)%x(1,jk) + dzi1
                            dzj2 = this%c(jx,jy,jz)%x(2,jk) + dzi2
                            dzj3 = this%c(jx,jy,jz)%x(3,jk) + dzi3
            
                        !---    find separation between points in real space
                            dzij1 = a11*dzj1 + a12*dzj2 + a13*dzj3
                            dzij2 = a21*dzj1 + a22*dzj2 + a23*dzj3
                            dzij3 = a31*dzj1 + a32*dzj2 + a33*dzj3
    
                        !---    find square distance
                            
                            n = n + 1
                            id(n) = this%c(jx,jy,jz)%indx(jk)
                            r2(n) = dzij1*dzij1 + dzij2*dzij2 + dzij3*dzij3 
                        
                        end do
                    end do
                end do
            end do
            return
        end subroutine neighbourList1c 
 

        
        subroutine neighbourList2c( this,ix,iy,iz, n,id,dx )
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ 
    !*      return the neighbours about a point in real space x  in all 27 cells, no cutoff
            type(LinkCell3D),intent(in)                     ::      this
            integer,intent(in)                              ::      ix,iy,iz            
            integer,intent(out)                             ::      n
            integer,dimension(:),intent(out)                ::      id
            real(kind=real64),dimension(:,:),intent(out)    ::      dx
            
            
            integer                             ::      jx,jy,jz , jk
            integer                             ::      kx,ky,kz 
            real(kind=real64)                   ::      dzj1,dzj2,dzj3
            real(kind=real64)                   ::      a11,a21,a31,a12,a22,a32,a13,a23,a33
             
            
            n = 0
            
            
            
            
        !---    find position in reduced cell coords
           
            a11 = this%super%a(1,1)
            a21 = this%super%a(2,1)
            a31 = this%super%a(3,1)
            a12 = this%super%a(1,2)
            a22 = this%super%a(2,2)
            a32 = this%super%a(3,2)
            a13 = this%super%a(1,3)
            a23 = this%super%a(2,3)
            a33 = this%super%a(3,3) 
 
        !---    find neighbour cells
            do kz = -1,1
                jz = mod( iz + this%super%Nz + kz , this%super%Nz )
                
                
                do ky = -1,1
                    jy = mod( iy + this%super%Ny + ky , this%super%Ny )
  
                    
                    do kx = -1,1
                        jx = mod( ix + this%super%Nx + kx , this%super%Nx )
             
                                                                                                   
                        do jk = 1,this%c(jx,jy,jz)%nindx
                         
                            n = n + 1
                            id(n) = this%c(jx,jy,jz)%indx(jk)
                        
                        !---    find separation between points in cell space
                            dzj1 = this%c(jx,jy,jz)%x(1,jk) + kx
                            dzj2 = this%c(jx,jy,jz)%x(2,jk) + ky
                            dzj3 = this%c(jx,jy,jz)%x(3,jk) + kz
            
                        !---    find separation between points in real space
                            dx(1,n) = a11*dzj1 + a12*dzj2 + a13*dzj3
                            dx(2,n) = a21*dzj1 + a22*dzj2 + a23*dzj3
                            dx(3,n) = a31*dzj1 + a32*dzj2 + a33*dzj3
    
                        end do
                    end do
                end do
            end do
            return
        end subroutine neighbourList2c 
 

        subroutine neighbourList1c_27( this,x, n,id,r2 )
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      return the neighbours about a point in real space x in all 27 cells, no cutoff, when there is only 27 cells 
            type(LinkCell3D),intent(in)                     ::      this
            real(kind=real64),dimension(3),intent(in)       ::      x            
            integer,intent(out)                             ::      n
            integer,dimension(:),intent(out)                ::      id
            real(kind=real64),dimension(:),intent(out)      ::      r2
            
            !integer                             ::      ix,iy,iz
            integer                             ::      jx,jy,jz , jk
            !integer                             ::      kx,ky,kz 
            real(kind=real64),dimension(3)      ::      yy!,zi!,dzi
            real(kind=real64)                   ::      dzi1,dzi2,dzi3 ,dzij1,dzij2,dzij3 ,dzj1,dzj2,dzj3      !,dd
            real(kind=real64)                   ::      a11,a21,a31,a12,a22,a32,a13,a23,a33
            
            n = 0
            
             
        !---    find position in reduced cell coords
            yy = realSpaceToCell( this%super,x )
            a11 = this%super%a(1,1)
            a21 = this%super%a(2,1)
            a31 = this%super%a(3,1)
            a12 = this%super%a(1,2)
            a22 = this%super%a(2,2)
            a32 = this%super%a(3,2)
            a13 = this%super%a(1,3)
            a23 = this%super%a(2,3)
            a33 = this%super%a(3,3) 
 
        !---    find neighbour cells
            do jz = 0,this%super%Nz-1               
                dzi3 = jz - yy(3) 
                do jy = 0,this%super%Ny-1 
                    dzi2 = jy - yy(2)
                    do jx = 0,this%super%Nx-1     
                        dzi1 = jx - yy(1)                                                                        
                        
                        do jk = 1,this%c(jx,jy,jz)%nindx
                        
                        
                        !---    find separation between points in cell space
                            dzj1 = this%c(jx,jy,jz)%x(1,jk) + dzi1
                            dzj2 = this%c(jx,jy,jz)%x(2,jk) + dzi2
                            dzj3 = this%c(jx,jy,jz)%x(3,jk) + dzi3
                            
                        !---    minimum image convention
                            if ( 2*dzj1 > this%super%Nx ) then
                                dzj1 = dzj1 - this%super%Nx
                            else if ( 2*dzj1 <= -this%super%Nx ) then
                                dzj1 = dzj1 + this%super%Nx
                            end if
                            if ( 2*dzj2 > this%super%Ny ) then
                                dzj2 = dzj2 - this%super%Ny
                            else if ( 2*dzj2 <= -this%super%Ny ) then
                                dzj2 = dzj2 + this%super%Ny
                            end if
                            if ( 2*dzj3 > this%super%Nz ) then
                                dzj3 = dzj3 - this%super%Nz
                            else if ( 2*dzj3 <= -this%super%Nz ) then
                                dzj3 = dzj3 + this%super%Nz
                            end if
                             
                            
                            
                        !---    find separation between points in real space
                            dzij1 = a11*dzj1 + a12*dzj2 + a13*dzj3
                            dzij2 = a21*dzj1 + a22*dzj2 + a23*dzj3
                            dzij3 = a31*dzj1 + a32*dzj2 + a33*dzj3
    
                        !---    find square distance
                            
                            n = n + 1
                            id(n) = this%c(jx,jy,jz)%indx(jk)
                            r2(n) = dzij1*dzij1 + dzij2*dzij2 + dzij3*dzij3 
                        
                        end do
                    end do
                end do
            end do
            return
        end subroutine neighbourList1c_27


        subroutine neighbourList2( this,ix,iy,iz, n,dx )
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      return the neighbours about cell ix,iy,iz
            type(LinkCell3D),intent(in)                     ::      this
            integer,intent(in)                              ::      ix,iy,iz
            integer,intent(out)                             ::      n
            real(kind=real64),dimension(:,:),intent(out)    ::      dx
            
            integer                             ::      jx,jy,jz , jk
            integer                             ::      kx,ky,kz 
            real(kind=real64),dimension(3)      ::      dzij
            
            n = 0            
           
        !---    find neighbour cells
            do kz = -1,1
                jz = mod( iz + this%super%Nz + kz , this%super%Nz )
                do ky = -1,1
                    jy = mod( iy + this%super%Ny + ky , this%super%Ny )
                    do kx = -1,1
                        jx = mod( ix + this%super%Nx + kx , this%super%Nx )
                                                                           
                        do jk = 1,this%c(jx,jy,jz)%nindx
                                               
                        !---    find separation between points in cell space
                            dzij(1:3) = this%c(jx,jy,jz)%x(1:3,jk) + (/ kx,ky,kz /) 
            
                        !---    find separation between points in real space
                            n = n + 1
                            dx(1:3,n) = cellToRealSpace( this%super,dzij )
                            
                        end do
                    end do
                end do
            end do
            return
        end subroutine neighbourList2         
        
        
        subroutine neighbourList2b( this,x, n,id,dx )
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ 
    !*      return the neighbours about a point in real space x within all neighbour cells, no cutoff
            type(LinkCell3D),intent(in)                     ::      this
            real(kind=real64),dimension(3),intent(in)       ::      x
            integer,intent(out)                             ::      n
            integer,dimension(:),intent(out)                ::      id
            real(kind=real64),dimension(:,:),intent(out)    ::      dx
            
            integer                             ::      ix,iy,iz
            integer                             ::      jx,jy,jz , jk
            integer                             ::      kx,ky,kz 
            real(kind=real64),dimension(3)      ::      yy,zi!,dzi
            real(kind=real64)                   ::      dzi1,dzi2,dzi3   ,dzj1,dzj2,dzj3 !,dd
            real(kind=real64)                   ::      a11,a21,a31,a12,a22,a32,a13,a23,a33
            n = 0
            
            
            if (getNnodes(this%super)==27) then
                call neighbourList0c_27( this,x, n,id,dx )
                return
            end if
        !---    find position in reduced cell coords
            yy = realSpaceToCell( this%super,x )
            a11 = this%super%a(1,1)
            a21 = this%super%a(2,1)
            a31 = this%super%a(3,1)
            a12 = this%super%a(1,2)
            a22 = this%super%a(2,2)
            a32 = this%super%a(3,2)
            a13 = this%super%a(1,3)
            a23 = this%super%a(2,3)
            a33 = this%super%a(3,3) 

        !---    find offset from cell origin in cell coords
            zi(1:3) = yy(1:3) - floor(yy(1:3))
                                    
        !---    find cell
            call whichCell( this%super,yy,ix,iy,iz,cellSpace = .true. )

        !---    find neighbour cells
            do kz = -1,1
                jz = mod( iz + this%super%Nz + kz , this%super%Nz )
                dzi3 = kz - zi(3)
                
                do ky = -1,1
                    jy = mod( iy + this%super%Ny + ky , this%super%Ny )
                    dzi2 = ky - zi(2)
                    
                    do kx = -1,1
                        jx = mod( ix + this%super%Nx + kx , this%super%Nx )
                        dzi1 = kx - zi(1)
                                                                                                   
                        do jk = 1,this%c(jx,jy,jz)%nindx
                        
                        
                        !---    find separation between points in cell space
                            dzj1 = this%c(jx,jy,jz)%x(1,jk) + dzi1
                            dzj2 = this%c(jx,jy,jz)%x(2,jk) + dzi2
                            dzj3 = this%c(jx,jy,jz)%x(3,jk) + dzi3
            
    
                        !---    don't compare to cutoff- take all
                            n = n + 1
                            id(n) = this%c(jx,jy,jz)%indx(jk)
                            dx(1,n) = a11*dzj1 + a12*dzj2 + a13*dzj3
                            dx(2,n) = a21*dzj1 + a22*dzj2 + a23*dzj3
                            dx(3,n) = a31*dzj1 + a32*dzj2 + a33*dzj3
                            
                        
                        end do
                    end do
                end do
            end do
            return
        end subroutine neighbourList2b 
        
    !---
    
    

        subroutine neighbourList_long1b( this,x,sigma, n,id,r2 )
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      return distance squared to the neighbours about a point in real space x within sphere radius sigma
            type(LinkCell3D),intent(in)                     ::      this
            real(kind=real64),dimension(3),intent(in)       ::      x
            real(kind=real64),intent(in)                    ::      sigma
            integer,intent(out)                             ::      n
            integer,dimension(:),intent(out)                ::      id
            real(kind=real64),dimension(:),intent(out)      ::      r2
            
            integer                             ::      ix,iy,iz
            integer                             ::      jx,jy,jz , jk
            integer                             ::      kx,ky,kz 
            integer                             ::      mx,my,mz
            real(kind=real64),dimension(3)      ::      yy,zi,dzi,dzij
            real(kind=real64)                   ::      dd
            
            n = 0
            
        !---    assuming a fairly cubic cell, how many cell repeats does sigma cover?
             
            mx = ceiling( sigma / norm2( this%super%a(1:3,1) ) ) 
            my = ceiling( sigma / norm2( this%super%a(1:3,2) ) ) 
            mz = ceiling( sigma / norm2( this%super%a(1:3,3) ) ) 
            
            
        !---    find position in reduced cell coords
            yy = realSpaceToCell( this%super,x )

        !---    find offset from cell origin in cell coords
            zi(1:3) = yy(1:3) - floor(yy(1:3))
                                    
        !---    find cell
            call whichCell( this%super,yy,ix,iy,iz,cellSpace = .true. )

            
            if (LinkCell3D_dbg) then
                print *,"Lib_LinkCell3d::neighbourList_long1b() info - range ",sigma," cells ",mx,my,mz
                print *,"Lib_LinkCell3d::neighbourList_long1b() info - central point ",x," in cell ",ix,iy,iz
            end if
            
            
        !---    find neighbour cells
            do kz = -mz,mz
                jz = mod( iz + this%super%Nz + kz , this%super%Nz )
                do ky = -my,my
                    jy = mod( iy + this%super%Ny + ky , this%super%Ny )
                    do kx = -mx,mx
                        jx = mod( ix + this%super%Nx + kx , this%super%Nx )
                        dzi(1:3) = (/ kx,ky,kz /) - zi(1:3)
                                                                           
                        do jk = 1,this%c(jx,jy,jz)%nindx
                        
                        
                        !---    find separation between points in cell space
                            dzij(1:3) = this%c(jx,jy,jz)%x(1:3,jk) + dzi(1:3)
            
                        !---    find separation between points in real space
                            dzij = cellToRealSpace( this%super,dzij )
    
                        !---    find square distance
                            dd = dzij(1)*dzij(1) + dzij(2)*dzij(2) + dzij(3)*dzij(3) 
    
                        !---    compare to cutoff
                            if (dd <= sigma*sigma) then
                            
                                if (n == size(id)) then
                                    print *,"Lib_LinkCell3d::neighbourList_long1b() error - inadequate size for id array, need > ",n
                                    return
                                end if
                            
                                n = n + 1
                                id(n) = this%c(jx,jy,jz)%indx(jk)
                                r2(n) = dd
                                
                                !print *,"match ",yy," to ",jx,jy,jz,jk," ",dzij," ",id(n),dd
                                
                                
                            end if
                        end do
                    end do
                end do
            end do
            return
        end subroutine neighbourList_long1b 
        

        subroutine neighbourList_long1c( this,x,sigma, n,id,dx,r2 )
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      return distance squared to the neighbours about a point in real space x within sphere radius sigma
            type(LinkCell3D),intent(in)                     ::      this
            real(kind=real64),dimension(3),intent(in)       ::      x
            real(kind=real64),intent(in)                    ::      sigma
            integer,intent(out)                             ::      n
            integer,dimension(:),intent(out)                ::      id
            real(kind=real64),dimension(:,:),intent(out)    ::      dx
            real(kind=real64),dimension(:),intent(out)      ::      r2
            
            integer                             ::      ix,iy,iz
            integer                             ::      jx,jy,jz , jk
            integer                             ::      kx,ky,kz 
            integer                             ::      mx,my,mz
            real(kind=real64),dimension(3)      ::      yy,zi,dzi,dzij
            real(kind=real64)                   ::      dd
            
            n = 0
            
        !---    assuming a fairly cubic cell, how many cell repeats does sigma cover?
             
            mx = ceiling( sigma / norm2( this%super%a(1:3,1) ) ) 
            my = ceiling( sigma / norm2( this%super%a(1:3,2) ) ) 
            mz = ceiling( sigma / norm2( this%super%a(1:3,3) ) ) 
            
            
        !---    find position in reduced cell coords
            yy = realSpaceToCell( this%super,x )

        !---    find offset from cell origin in cell coords
            zi(1:3) = yy(1:3) - floor(yy(1:3))
                                    
        !---    find cell
            call whichCell( this%super,yy,ix,iy,iz,cellSpace = .true. )

            
            if (LinkCell3D_dbg) then
                print *,"Lib_LinkCell3d::neighbourList_long1c() info - range ",sigma," cells ",mx,my,mz
                print *,"Lib_LinkCell3d::neighbourList_long1c() info - central point ",x," in cell ",ix,iy,iz
            end if
            
            
        !---    find neighbour cells
            do kz = -mz,mz
                jz = mod( iz + this%super%Nz + kz , this%super%Nz )
                do ky = -my,my
                    jy = mod( iy + this%super%Ny + ky , this%super%Ny )
                    do kx = -mx,mx
                        jx = mod( ix + this%super%Nx + kx , this%super%Nx )
                        dzi(1:3) = (/ kx,ky,kz /) - zi(1:3)
                                                                           
                        do jk = 1,this%c(jx,jy,jz)%nindx
                        
                        
                        !---    find separation between points in cell space
                            dzij(1:3) = this%c(jx,jy,jz)%x(1:3,jk) + dzi(1:3)
            
                        !---    find separation between points in real space
                            dzij = cellToRealSpace( this%super,dzij )
    
                        !---    find square distance
                            dd = dzij(1)*dzij(1) + dzij(2)*dzij(2) + dzij(3)*dzij(3) 
    
                        !---    compare to cutoff
                            if (dd <= sigma*sigma) then
                             
                                n = n + 1
                                id(n) = this%c(jx,jy,jz)%indx(jk)
                                r2(n) = dd
                                dx(1:3,n) = dzij(1:3)
                                
                                
                            end if
                        end do
                    end do
                end do
            end do
            return
        end subroutine neighbourList_long1c 
        

        
        
        subroutine neighbourList_long1d( this,x,sigma, n  )
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ 
    !*      return count of neighbours in range sigma
            type(LinkCell3D),intent(in)                     ::      this
            real(kind=real64),dimension(3),intent(in)       ::      x
            real(kind=real64),intent(in)                    ::      sigma
            integer,intent(out)                             ::      n
            
            integer                             ::      ix,iy,iz
            integer                             ::      jx,jy,jz , jk
            integer                             ::      kx,ky,kz 
            integer                             ::      mx,my,mz
            real(kind=real64),dimension(3)      ::      yy,zi,dzi,dzij
            real(kind=real64)                   ::      dd
            
            n = 0
            
        !---    assuming a fairly cubic cell, how many cell repeats does sigma cover?
             
            mx = ceiling( sigma / norm2( this%super%a(1:3,1) ) ) 
            my = ceiling( sigma / norm2( this%super%a(1:3,2) ) ) 
            mz = ceiling( sigma / norm2( this%super%a(1:3,3) ) ) 
            
            
        !---    find position in reduced cell coords
            yy = realSpaceToCell( this%super,x )

        !---    find offset from cell origin in cell coords
            zi(1:3) = yy(1:3) - floor(yy(1:3))
                                    
        !---    find cell
            call whichCell( this%super,yy,ix,iy,iz,cellSpace = .true. )

            
            if (LinkCell3D_dbg) then
                print *,"Lib_LinkCell3d::neighbourList_long1d() info - range ",sigma," cells ",mx,my,mz
                print *,"Lib_LinkCell3d::neighbourList_long1d() info - central point ",x," in cell ",ix,iy,iz
            end if
            
            
        !---    find neighbour cells
            do kz = -mz,mz
                jz = mod( iz + this%super%Nz + kz , this%super%Nz )
                do ky = -my,my
                    jy = mod( iy + this%super%Ny + ky , this%super%Ny )
                    do kx = -mx,mx
                        jx = mod( ix + this%super%Nx + kx , this%super%Nx )
                        dzi(1:3) = (/ kx,ky,kz /) - zi(1:3)
                                                                           
                        do jk = 1,this%c(jx,jy,jz)%nindx
                        
                        
                        !---    find separation between points in cell space
                            dzij(1:3) = this%c(jx,jy,jz)%x(1:3,jk) + dzi(1:3)
            
                        !---    find separation between points in real space
                            dzij = cellToRealSpace( this%super,dzij )
    
                        !---    find square distance
                            dd = dzij(1)*dzij(1) + dzij(2)*dzij(2) + dzij(3)*dzij(3) 
    
                        !---    compare to cutoff
                            if (dd <= sigma*sigma) n = n + 1
                        end do
                    end do
                end do
            end do
            return
        end subroutine neighbourList_long1d 
        

         
        
        
    !---
        
        
        pure function getnNeighMax(this) result(n)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      return a bounding upper limit for the number of possible neighbours returned by neighbourList
    !*      assuming all cells are maximally densely populated.
            type(LinkCell3D),intent(in)                     ::      this
            integer                                         ::      n
            n = min( sum( this%c(:,:,:)%nindx ),this%nMax*27 )
            return
        end function getnNeighMax
         
        
        pure function getnPoints0(this) result(n)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ 
    !*      return total count of all points in all cells
            type(LinkCell3D),intent(in)                     ::      this
            integer                                         ::      n
            n = sum( this%c(:,:,:)%nindx )
            return
        end function getnPoints0
        
        
        pure function getnPoints1(this,ix,iy,iz) result(n)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      return count of points in cell (ix,iy,iz)
            type(LinkCell3D),intent(in)                     ::      this
            integer,intent(in)                              ::      ix,iy,iz
            integer                                         ::      n
            n = this%c(ix,iy,iz)%nindx 
            return
        end function getnPoints1
        
        
        pure function getMaxId( this ) result(id)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      returns max id in any link cell.
            type(LinkCell3D),intent(in)                     ::      this
            integer                                         ::      ix,iy,iz,ik
            integer                                         ::      id
            id = -huge(1)
            do iz = 0,this%super%Nz-1
                do iy = 0,this%super%Ny-1
                    do ix = 0,this%super%Nz-1
                        do ik = 1,this%c(ix,iy,iz)%nindx
                            id = max( id,this%c(ix,iy,iz)%indx(ik) )
                        end do
                    end do
                end do
            end do
            return
        end function getMaxId        
         
        
        
        pure function getnMax(this) result(n)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      return highest point count in any one cell
            type(LinkCell3D),intent(in)                     ::      this
            integer                                         ::      n
            n = this%nMax
            return
        end function getnMax
                            
        pure function getMaxPerCell(this) result(n)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      return highest storage in any one cell
            type(LinkCell3D),intent(in)                     ::      this
            integer                                         ::      n
            n = maxval( this%c(:,:,:)%nMax )
            return
        end function getMaxPerCell
        
        
    end module Lib_LinkCell3D
    
        
        
    