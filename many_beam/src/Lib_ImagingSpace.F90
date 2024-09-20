

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
        public      ::      getBounds


        integer,private                     ::      rank = 0
        integer,private                     ::      nProcs = 1         
        integer,public                      ::      Lib_ImagingSpace_SIGMA_MULT = 2


        type,public     ::      ImagingSpace
            private
            real(kind=real64)               ::      a                   !   grid side length
            real(kind=real64)               ::      sigma               !   smoothing length scale
            real(kind=real64),dimension(3)  ::      delta
            integer                         ::      nBuf                !   number of buffing cells required to ensure smoothness at boundaries
            integer                         ::      Nx,Ny,Nz
            integer                         ::      lbx,ubx,lby,uby     !   the region of image space cells I am responsible for
            integer,dimension(:,:),pointer  ::      p                   !   (0:Mx-1,0:My-1) the processor who is responsible for each block
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
 
        interface   inMyCell
            module procedure        inMyCell0
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
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^    
            type(ImagingSpace)                     ::      this
            this%a = 1.0d0
            this%sigma = 0.5d0
            this%delta = 0.0d0
            this%nBuf = buffer(this%a,this%sigma)
            this%Nx = 1
            this%Ny = 1
            this%Nz = 1            
            nullify(this%p)
            call setMyBounds(this)
            return
        end function ImagingSpace_null
        

        function ImagingSpace_ctor0(a,sigma,delta,Nx,Ny,Nz) result(this)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^    
            type(ImagingSpace)                          ::      this
            real(kind=real64),intent(in)                ::      a
            real(kind=real64),intent(in)                ::      sigma
            real(kind=real64),dimension(3),intent(in)   ::      delta
            integer,intent(in)                          ::      Nx,Ny,Nz
            this%a = a
            this%sigma = sigma
            this%delta = delta
            this%nBuf = buffer(this%a,this%sigma)
            this%Nx = Nx
            this%Ny = Ny
            this%Nz = Nz
            nullify(this%p)
            call setMyBounds(this)
            return
        end function ImagingSpace_ctor0

        function ImagingSpace_ctor1(a,delta,Nx,Ny,Nz) result(this)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^    
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

        subroutine setMyBounds(this)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      given my rank and number of processes, determine the block of cells I am responsible for 
    !*      note that this will _not_ include any buffer
            type(ImagingSpace),intent(inout)            ::      this

            integer         ::          Mx,My                       !   number of blocks in x and y dimension 
            integer         ::          Cx,Cy                       !   size of each block
            integer         ::          pp,ix,iy                    !   processor number, position in grid
            
            call FactoriseParallel( this%Nx,this%Ny,nProcs,Mx,My )
            
            allocate(this%p(0:Mx-1,0:My-1))
            Cx = ceiling(real(this%Nx)/Mx)
            Cy = ceiling(real(this%Ny)/My)

            do pp = 0,nProcs-1
                iy = pp/Mx
                ix = pp - Mx*iy
                this%p(ix,iy) = pp
                if (pp == rank) then
                    this%lbx = ix*Cx
                    this%ubx = min( this%Nx-1,this%lbx + Cx - 1 )
                    this%lby = iy*Cy
                    this%uby = min( this%Ny-1,this%lby + Cy - 1 )                    
                end if
            end do

            return
        end subroutine setMyBounds



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