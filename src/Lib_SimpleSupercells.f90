
    module Lib_SimpleSupercells
!---^^^^^^^^^^^^^^^^^^^^^^^^^^^
!*      provide bare functionality for a simple periodic 3d supercell
!*      constructed from a set of unit cells
!*      offers periodic boundaries, cell space to real space functions
!*
!*      Daniel Mason, UKAEA
!*      April 2022
!*

        use Lib_RotationMatrices        !   for complete basis call
        use iso_fortran_env
        implicit none
        private
        
    !---
    
        public      ::      SimpleSupercell_ctor
        public      ::      delete
        public      ::      report
     
    !---
    
        public      ::      cellToRealSpace
        public      ::      realSpaceToCell
        public      ::      getNx,getNy,getNz 
        public      ::      getNnodes 
        public      ::      setNx,setNy,setNz
        public      ::      getA,getSuperA,getiA
        public      ::      getLatticeVector,getSuperVector
        public      ::      getSuperCellSideLength
        public      ::      getCellSideLength
        public      ::      setA    
        public      ::      wrapPBC
        public      ::      whichCell
        public      ::      unitCellVolume
        public      ::      superCellVolume
        public      ::      minimumImage
        public      ::      pointInSupercell
        public      ::      suggestSupercellOrientedWithN
        
        public      ::      linearInterpolation
        
        
        
    !---
    
        logical,public          ::      SimpleSupercell_dbg = .false.
                   
    !---
    
        type,public     ::      SimpleSupercell
            !   note that the fields are accessible publically
            integer                                 ::      Nx,Ny,Nz        !   number of repeats
            real(kind=real64),dimension(3,3)        ::      a               !   unit cell vectors - vector 1 is a(1:3,1) etc
            real(kind=real64),dimension(3,3)        ::      ia              !   inverse unit cell vectors             
        end type SimpleSupercell
        
    !---
    
        interface SimpleSupercell_ctor
            module procedure    SimpleSupercell_null
            module procedure    SimpleSupercell_ctor0
        end interface
                
        interface delete
            module procedure    delete0
        end interface
        
        interface report
            module procedure    report0
        end interface
        

        interface   getA
            module procedure        getA0
        end interface
        
        interface   getiA
            module procedure        getiA0
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
        
        interface   getNnodes
            module procedure        getNnodes0
        end interface   
                 

        interface   cellToRealSpace
            module procedure        cellToRealSpace0
            module procedure        cellToRealSpace0a
            module procedure        cellToRealSpace1
            module procedure        cellToRealSpace2
        end interface

        interface   realSpaceToCell
            module procedure        realSpaceToCell0
        end interface

        interface   unitCellVolume
            module procedure        unitCellVolume0
        end interface
        interface   superCellVolume
            module procedure        superCellVolume0
            module procedure        superCellVolume1
        end interface
        
        
        interface   setA
            module procedure        setA0
        end interface
!        
        
        interface   wrapPBC
            module procedure        wrapPBC0
            module procedure        wrapPBC1
            module procedure        wrapPBC2
        end interface
        
         
        interface   pointInSupercell
            module procedure        pointInSupercell0
        end interface
!        
        interface   linearInterpolation
            module procedure        linearInterpolation0
            module procedure        linearInterpolation2d
        end interface
!        
        
        interface   suggestSupercellOrientedWithN
            module procedure        suggestSupercellOrientedWithN0
            module procedure        suggestSupercellOrientedWithN1
            module procedure        suggestSupercellOrientedWithNandG
        end interface
        
        interface   suggestPlane
            module procedure        suggestPlane0
            module procedure        suggestPlane1
        end interface
        
        

        interface   getCellSideLength
            module procedure    getCellSideLength0
        end interface
        
        interface   getSupercellSideLength
            module procedure    getSupercellSideLength0
        end interface
        
        interface   getLatticeVector
            module procedure    getLatticeVector0
        end interface
        
        interface   getSuperVector
            module procedure    getSuperVector0
        end interface
        
        
    contains
!---^^^^^^^^

        function SimpleSupercell_null() result(this)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      default supercell of side 1
            type(SimpleSupercell)           ::      this
            this = SimpleSupercell_ctor0( reshape( (/1.0d0,0.0d0,0.0d0 , 0.0d0,1.0d0,0.0d0 , 0.0d0,0.0d0,1.0d0/),(/3,3/) ),1,1,1 )
            return
        end function SimpleSupercell_null
                         
        pure function SimpleSupercell_ctor0(a,Nx,Ny,Nz) result(this)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      simple supercell with Nx,Ny,Nz unit cells
    !*      whose lattice vectors a_ij is the ith cartesian component of the jth vector
            type(SimpleSupercell)                           ::      this
            real(kind=real64),dimension(3,3),intent(in)     ::      a
            integer,intent(in)          ::      Nx,Ny,Nz
            this%a = a
            this%Nx = Nx
            this%Ny = Ny
            this%Nz = Nz
            call inverse3Mat( this%a,this%ia )             
            return
        end function SimpleSupercell_ctor0
                       
    !---
    
        subroutine delete0(this)
    !---^^^^^^^^^^^^^^^^^^^^^^^^
    !*      no dynamic memory
            type(SimpleSupercell),intent(inout)    ::      this            
            this = SimpleSupercell_null()
            return
        end subroutine delete0
        
    !---
    
        subroutine report0(this,u,o,verbose)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      simple output. Defaults to unit 6 = screen.
    !*      optional argument o determines left hand margin
    !*      optionally gives verbose output    
            type(SimpleSupercell),intent(in)    ::      this
            integer,intent(in),optional         ::      u,o
            logical,intent(in),optional         ::      verbose

            integer         ::      uu,oo
            
            real(kind=real64),parameter     ::      DEGINRAD = 180.0d0/3.141592654d0
            uu = 6 ; if (present(u)) uu = u
            oo = 0 ; if (present(o)) oo = o


            write(unit=uu,fmt='(3(a,i4),a)')  repeat(" ",oo)//"Supercell[Nx,Ny,Nz = ",this%Nx,",",this%Ny,",",this%Nz," ]"
            write(unit=uu,fmt='(a)')        repeat(" ",oo+2)//"lattice vectors"
            write(unit=uu,fmt='(a,3f14.8)') repeat(" ",oo+2),this%a(1,1),this%a(1,2),this%a(1,3)
            write(unit=uu,fmt='(a,3f14.8)') repeat(" ",oo+2),this%a(2,1),this%a(2,2),this%a(2,3)
            write(unit=uu,fmt='(a,3f14.8)') repeat(" ",oo+2),this%a(3,1),this%a(3,2),this%a(3,3)                
            if (present(verbose)) then
                if (verbose) then
                    write(unit=uu,fmt='(a)')        repeat(" ",oo+2)//"lattice lengths"
                    write(unit=uu,fmt='(a,3f14.8)') repeat(" ",oo+2),norm2(this%a(:,1)),norm2(this%a(:,2)),norm2(this%a(:,3))
                    write(unit=uu,fmt='(a)')        repeat(" ",oo+2)//"lattice angles (deg)"
                    write(unit=uu,fmt='(a,3f14.8)') repeat(" ",oo+2),DEGINRAD*acos( dot_product(this%a(:,2),this%a(:,3))/(norm2(this%a(:,2))*norm2(this%a(:,3))) )          &
                                                                    ,DEGINRAD*acos( dot_product(this%a(:,3),this%a(:,1))/(norm2(this%a(:,3))*norm2(this%a(:,1))) )          &
                                                                    ,DEGINRAD*acos( dot_product(this%a(:,1),this%a(:,2))/(norm2(this%a(:,1))*norm2(this%a(:,2))) )
                        
                end if
            end if
            return
        end subroutine report0

    
    !---
    

    !---
    
        pure function unitCellVolume0( this ) result(V)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      returns the volume of a single unit cell 
            type(SimpleSupercell),intent(in)                ::      this
            real(kind=real64)                               ::      V
            V = determinant3Mat(this%a)
            return
        end function unitCellVolume0

    
        pure function superCellVolume0( this ) result(V)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      returns the volume of the whole supercell
            type(SimpleSupercell),intent(in)                ::      this
            real(kind=real64)                               ::      V
            V = unitCellVolume0( this )
            V = V * this%Nx*this%Ny*this%Nz
            return
        end function superCellVolume0

        pure function superCellVolume1( a ) result(V)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      returns the volume of a cell defined by a 3x3 matrix 
            real(kind=real64),dimension(3,3),intent(in)     ::      a
            real(kind=real64)                               ::      V
            V = determinant3Mat( a ) 
            return
        end function superCellVolume1

    !---
                
    
        pure function cellToRealSpace0( this,x ) result( y )
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      returns the position in real space given a position in cell space
            type(SimpleSupercell),intent(in)                ::      this
            real(kind=real64),dimension(3),intent(in)       ::      x
            real(kind=real64),dimension(3)                  ::      y
            y(1:3) = this%a(1:3,1)*x(1) + this%a(1:3,2)*x(2) + this%a(1:3,3)*x(3)
            return
        end function cellToRealSpace0

        pure function cellToRealSpace0a( this,x,wrapPBC ) result( y )
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      returns the position in real space given a position in cell space    
    !*      wraps using periodic boundary conditions
            type(SimpleSupercell),intent(in)                ::      this
            real(kind=real64),dimension(3),intent(in)       ::      x
            logical,intent(in)                              ::      wrapPBC
            real(kind=real64),dimension(3)                  ::      y
            
            real(kind=real64),dimension(3)      ::      xx
            integer                             ::      ii
             
            if (.not. wrapPBC) then
                y(1:3) = this%a(1:3,1)*x(1) + this%a(1:3,2)*x(2) + this%a(1:3,3)*x(3)
                return
            end if
            
            xx = x     
            ii = floor(xx(1))
            xx(1) = xx(1) - ii
            ii = mod( ii+ishft(this%Nx,8),this%Nx )
            xx(1) = (ii+xx(1))
            ii = floor(xx(2))
            xx(2) = xx(2) - ii
            ii = mod( ii+ishft(this%Ny,8),this%Ny )
            xx(2) = (ii+xx(2))
            ii = floor(xx(3))
            xx(3) = xx(3) - ii
            ii = mod( ii+ishft(this%Nz,8),this%Nz )
            xx(3) = (ii+xx(3))
            
            y(1:3) = this%a(1:3,1)*xx(1) + this%a(1:3,2)*xx(2) + this%a(1:3,3)*xx(3)
            return
        end function cellToRealSpace0a

        pure function cellToRealSpace1( this,c ) result( y )
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      returns the position in real space given a position in cell space    
            type(SimpleSupercell),intent(in)                        ::      this
            integer,dimension(3),intent(in)                         ::      c
            real(kind=real64),dimension(3)                  ::      y
            y(1:3) = this%a(1:3,1)*c(1) + this%a(1:3,2)*c(2) + this%a(1:3,3)*c(3)
            return
        end function cellToRealSpace1

        pure function cellToRealSpace2( this,ix,iy,iz ) result( y )
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      returns the position in real space given a position in cell space    
            type(SimpleSupercell),intent(in)                        ::      this
            integer,intent(in)                                      ::      ix,iy,iz
            real(kind=real64),dimension(3)                  ::      y
            y(1:3) = this%a(1:3,1)*ix + this%a(1:3,2)*iy + this%a(1:3,3)*iz
            return
        end function cellToRealSpace2
        
        pure function realSpaceToCell0( this,y ) result( x )
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      returns the position in cell space given a position in real space    
            type(SimpleSupercell),intent(in)                ::      this
            real(kind=real64),dimension(:),intent(in)       ::      y
            real(kind=real64),dimension(3)                  ::      x
            x(1:3) = this%ia(1:3,1)*y(1) + this%ia(1:3,2)*y(2) + this%ia(1:3,3)*y(3)
            return
        end function realSpaceToCell0
        
        pure subroutine whichCell( this,y, ix,iy,iz , cellspace)  
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      which unit cell (ix,iy,iz) does position y correspond to? 
    !*      if cellspace then this is basically just checking for periodic boundary conditions
    !*      otherwise find the cell space position
            type(SimpleSupercell),intent(in)                ::      this
            real(kind=real64),dimension(:),intent(in)       ::      y
            logical,intent(in),optional                     ::      cellspace
            integer,intent(out)                             ::      ix,iy,iz
            real(kind=real64)                  ::      xx
            
            if (present(cellspace)) then
                if (cellspace) then
                    ix = mod( floor( y(1) ) + ishft(this%Nx,4),this%Nx )  
                    iy = mod( floor( y(2) ) + ishft(this%Ny,4),this%Ny )  
                    iz = mod( floor( y(3) ) + ishft(this%Nz,4),this%Nz )  
                    return
                end if
            end if
            xx = this%ia(1,1)*y(1) + this%ia(1,2)*y(2) + this%ia(1,3)*y(3)
            ix = mod( floor( xx ) + ishft(this%Nx,4),this%Nx ) 
            xx = this%ia(2,1)*y(1) + this%ia(2,2)*y(2) + this%ia(2,3)*y(3)
            iy = mod( floor( xx ) + ishft(this%Ny,4),this%Ny ) 
            xx = this%ia(3,1)*y(1) + this%ia(3,2)*y(2) + this%ia(3,3)*y(3)
            iz = mod( floor( xx ) + ishft(this%Nz,4),this%Nz ) 
            return
        end subroutine whichCell



        pure function minimumImage( this,dy ) result( dyp )
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      use the minimum image convention to transform the input vector dy 
    !*      to one which is inside the supercell 
            type(SimpleSupercell),intent(in)                ::      this
            real(kind=real64),dimension(3),intent(in)       ::      dy
            real(kind=real64),dimension(3)                  ::      dyp
            integer             ::      ii
            real(kind=real64),dimension(3)                  ::      dx
            dx(1:3) = realSpaceToCell( this,dy )
            dyp(1:3) = dy(1:3)
            
            ii = floor(dx(1))
            if (2*ii>this%nx)  dyp(1:3) = dyp(1:3) - this%nx*this%a(1:3,1)
            if (2*ii<-this%nx) dyp(1:3) = dyp(1:3) + this%nx*this%a(1:3,1)
            
            ii = floor(dx(2))
            if (2*ii>this%ny)  dyp(1:3) = dyp(1:3) - this%ny*this%a(1:3,2)
            if (2*ii<-this%ny) dyp(1:3) = dyp(1:3) + this%ny*this%a(1:3,2)
            
            ii = floor(dx(3))
            if (2*ii>this%nz)  dyp(1:3) = dyp(1:3) - this%nz*this%a(1:3,3)
            if (2*ii<-this%nz) dyp(1:3) = dyp(1:3) + this%nz*this%a(1:3,3)
            
            return
        end function minimumImage       
        
    !---        
        
         subroutine suggestSupercellOrientedWithN0( this, n , super , mx,my,mz )
    !----^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      suggest a new supercell with (as close as possible) the same volume
    !*      but oriented with z-axis along n
    !*      Aim for unit cell side lengths of a0, but optionally allow (mx) subdivisions per side
            type(SimpleSupercell),intent(in)            ::      this 
            real(kind=real64),dimension(3),intent(in)   ::      n
            type(SimpleSupercell),intent(out)           ::      super
            integer,intent(in),optional                 ::      mx,my,mz
        
            real(kind=real64),dimension(3)      ::      dx1,dx2,dx3
            real(kind=real64)                   ::      a0 
            real(kind=real64),dimension(3,3)    ::      aa
            integer             ::      nx,ny,nz

        !---    suggest a plane            
            a0 = unitCellVolume( this )**(1/3.0d0)
            call suggestPlane( getSuperA(this) ,n, dx1,dx2,dx3 )
            
            nx = max(1,nint( norm2(dx1)/a0 ))           !   unit cell equivalent count in direction-1
            ny = max(1,nint( norm2(dx2)/a0 ))           !    
            nz = max(1,nint( norm2(dx3)/a0 ))           !    
            
            if (present(mz)) then
                nx = nx * mx
                ny = ny * my
                nz = nz * mz
            end if
            
            aa(1:3,1) = dx1(1:3)/nx
            aa(1:3,2) = dx2(1:3)/ny
            aa(1:3,3) = dx3(1:3)/nz
            
            super = SimpleSupercell_ctor( aa,nx,ny,nz )
            
            return
        end subroutine suggestSupercellOrientedWithN0
                
         subroutine suggestSupercellOrientedWithN1( supera, n ,a0, super , mx,my,mz )
    !----^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      suggest a new supercell with (as close as possible) the same volume
    !*      but oriented with z-axis along n
    !*      Aim for unit cell side lengths of a0, but optionally allow (mx) subdivisions per side
            real(kind=real64),dimension(3,3),intent(in) ::      supera
            real(kind=real64),dimension(3),intent(in)   ::      n
            real(kind=real64),intent(in)                ::      a0 
            type(SimpleSupercell),intent(out)           ::      super
            integer,intent(in),optional                 ::      mx,my,mz
        
            real(kind=real64),dimension(3)      ::      dx1,dx2,dx3
            
            real(kind=real64),dimension(3,3)    ::      aa
            integer             ::      nx,ny,nz 

        !---    suggest a plane            
             
            call suggestPlane( supera ,n, dx1,dx2,dx3 )
            ! print *,"suggestSupercellOrientedWithN1 dbg - n ",n
            ! print *,"suggestSupercellOrientedWithN1 dbg - supera"
            ! print *,supera(1,:)
            ! print *,supera(2,:)
            ! print *,supera(3,:)
            ! print *,"suggestSupercellOrientedWithN1 dbg - plane"
            ! print *,"dx1 ",dx1," norm ",norm2(dx1)
            ! print *,"dx2 ",dx2," norm ",norm2(dx2)
            ! print *,"dx3 ",dx3," norm ",norm2(dx3)
            
            
            nx = max(1,nint( norm2(dx1)/a0 ))           !   unit cell equivalent count in direction-1
            ny = max(1,nint( norm2(dx2)/a0 ))           !    
            nz = max(1,nint( norm2(dx3)/a0 ))           !    
            
            !print *,"nx,ny,nz ",nx,ny,nz," mxx,myy,mzz ",mx,my,mz
            
            
            if (present(mz)) then
                nx = nx * mx
                ny = ny * my
                nz = nz * mz
            end if
            
            aa(1:3,1) = dx1(1:3)/nx
            aa(1:3,2) = dx2(1:3)/ny
            aa(1:3,3) = dx3(1:3)/nz
            
            super = SimpleSupercell_ctor( aa,nx,ny,nz )
            
            return
        end subroutine suggestSupercellOrientedWithN1
                
         subroutine suggestSupercellOrientedWithNandG( supera, n,g ,a0, super , mx,my,mz )
    !----^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      suggest a new supercell with (as close as possible) the same volume
    !*      but oriented with z-axis along n and x-axis along g
    !*      Aim for unit cell side lengths of a0, but optionally allow (mx) subdivisions per side
            real(kind=real64),dimension(3,3),intent(in) ::      supera
            real(kind=real64),dimension(3),intent(in)   ::      n,g
            real(kind=real64),intent(in)                ::      a0 
            type(SimpleSupercell),intent(out)           ::      super
            integer,intent(in),optional                 ::      mx,my,mz
        
            real(kind=real64),dimension(3)      ::      dx1,dx2,dx3
            
            real(kind=real64),dimension(3,3)    ::      aa
            integer             ::      nx,ny,nz 

        !---    suggest a plane            
             
            call suggestPlane( supera ,n,g, dx1,dx2,dx3 )
            
            nx = max(1,nint( norm2(dx1)/a0 ))           !   unit cell equivalent count in direction-1
            ny = max(1,nint( norm2(dx2)/a0 ))           !    
            nz = max(1,nint( norm2(dx3)/a0 ))           !    
            
            
          !print *,"suggestSupercellOrientedWithNandG dbg - n ",n 
          !print *,"suggestSupercellOrientedWithNandG dbg - g ",g  
          !print *,"suggestSupercellOrientedWithNandG dbg - supera"
          !print *,supera(1,:)                                  
          !print *,supera(2,:)                                  
          !print *,supera(3,:)                                  
          !print *,"suggestSupercellOrientedWithNandG dbg - plane" 
          !print *,"dx1 ",dx1," norm ",norm2(dx1)               
          !print *,"dx2 ",dx2," norm ",norm2(dx2)               
          !print *,"dx3 ",dx3," norm ",norm2(dx3)               
            
            if (present(mz)) then
                nx = nx * mx
                ny = ny * my
                nz = nz * mz
            end if
            
            aa(1:3,1) = dx1(1:3)/nx
            aa(1:3,2) = dx2(1:3)/ny
            aa(1:3,3) = dx3(1:3)/nz
            
            super = SimpleSupercell_ctor( aa,nx,ny,nz )
            
            return
        end subroutine suggestSupercellOrientedWithNandG
                
    !---
            

        subroutine suggestPlane0( a_super,n, dx1,dx2,dx3 )
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      given that we want to sample points on a plane r.n = 0
    !*      suggest three vectors dx1,dx2,dx3
    !*      which cover the supercell
            real(kind=real64),dimension(3,3),intent(in)     ::      a_super
            real(kind=real64),dimension(3),intent(in)       ::      n
            real(kind=real64),dimension(3),intent(out)      ::      dx1,dx2,dx3

            real(kind=real64),dimension(3)      ::      nn,pp
            real(kind=real64)                   ::      d1,d2,d3,vol
            
        !---    which face of the supercell is closest to n?
        !   use this to decide dx1 and dx2.
        
        !   find projection of given normal on each face normal
            nn = cross_product( a_super(1:3,1),a_super(1:3,2) )
            pp(1) = dot_product( nn(1:3),n(1:3) ) / norm2(nn)
            nn = cross_product( a_super(1:3,2),a_super(1:3,3) )
            pp(2) = dot_product( nn(1:3),n(1:3) ) / norm2(nn)
            nn = cross_product( a_super(1:3,3),a_super(1:3,1) )
            pp(3) = dot_product( nn(1:3),n(1:3) ) / norm2(nn)
            
            
        !   which has gratest projection?
            pp = abs(pp)            
            nn(1:3) = n(1:3) / norm2(n(1:3))
            
            !print *,"Lib_SimpleSupercells::suggestPlane() dbg - n ",nn," p ",pp
            
            
            if (pp(1)>=max(pp(2),pp(3))) then
                !   n has greatest projection on a1,a2 face
                dx1 = a_super(1:3,1)                !   hint to push vector 1 along a1
                call completeBasis( nn,dx1,dx2 )                                
                !   find projection along a1,a2 vectors
                d1 = dot_product( dx1,a_super(1:3,1) )
                d2 = dot_product( dx2,a_super(1:3,2) )                

            else if (pp(2)>max(pp(3),pp(1))) then
                !   n has greatest projection on a2,a3 face
                dx1 = a_super(1:3,2)                !   hint to push vector 1 along a2
                call completeBasis( nn,dx1,dx2 )                
                !   find projection along a2,a3 vectors
                d1 = dot_product( dx1,a_super(1:3,2) )
                d2 = dot_product( dx2,a_super(1:3,3) )                
            
            else !if (pp(3)>max(pp(1),pp(2))) then
                !   n has greatest projection on a3,a1 face
                dx1 = a_super(1:3,3)                !   hint to push vector 1 along a3
                call completeBasis( nn,dx1,dx2 )                
                !   find projection along a3,a1 vectors
                d1 = dot_product( dx1,a_super(1:3,3) )
                d2 = dot_product( dx2,a_super(1:3,1) )                
            
            
            end if

            !print *,"Lib_SimpleSupercells::suggestPlane() dbg - dx1,dx2,d1,d2 ",dx1,",",dx2,",",d1,d2
            
        !---    can now complete the job by ensuring resultant box has same volume as initial box            
            vol = determinant3Mat(a_super)
            d3 = vol / (d1*d2)
            
            dx1 = dx1*d1
            dx2 = dx2*d2
            dx3 = nn*d3
            
            

            return
        end subroutine suggestPlane0



        subroutine suggestPlane1( a_super,n,g, dx1,dx2,dx3 )
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      given that we want to sample points on a plane r.n = 0
    !*      suggest three vectors dx1,dx2,dx3
    !*      which cover the supercell, with dx1 pointing along the g direction.
    !*      
            real(kind=real64),dimension(3,3),intent(in)     ::      a_super
            real(kind=real64),dimension(3),intent(in)       ::      n,g
            real(kind=real64),dimension(3),intent(out)      ::      dx1,dx2,dx3

            real(kind=real64),dimension(3)      ::      nn,pp,qq
            real(kind=real64)                   ::      d1,d2,d3,vol
             
            
            
        !---    which face of the supercell is closest to g?
        !   use this to decide length of new x- direction
        
            
        !   normalise z- direction
            dx3(1:3) = n(1:3) / norm2(n(1:3))
        !   x- direction hint will be along g
            dx1(1:3) = g(1:3) / norm2(g(1:3))
        !   complete the basis to find a normalised x- and y- direction.
            call completeBasis( dx3,dx1,dx2 )       
            
          !  print *,"Lib_SimpleSupercells::suggestPlane1 info - completed basis"
          !  write (*,fmt='(3f12.6)') dx1(1),dx2(1),dx3(1)
          !  write (*,fmt='(3f12.6)') dx1(2),dx2(2),dx3(2)
          !  write (*,fmt='(3f12.6)') dx1(3),dx2(3),dx3(3)
            

        !   find projection of x- and y- on each face normal
            nn = cross_product( a_super(1:3,1),a_super(1:3,2) )  
            pp(1) = dot_product( nn(1:3),dx1(1:3) ) / norm2(nn)
            qq(1) = dot_product( nn(1:3),dx2(1:3) ) / norm2(nn)
            nn = cross_product( a_super(1:3,2),a_super(1:3,3) )
            pp(2) = dot_product( nn(1:3),dx1(1:3) ) / norm2(nn)
            qq(2) = dot_product( nn(1:3),dx2(1:3) ) / norm2(nn)
            nn = cross_product( a_super(1:3,3),a_super(1:3,1) )
            pp(3) = dot_product( nn(1:3),dx1(1:3) ) / norm2(nn)
            qq(3) = dot_product( nn(1:3),dx2(1:3) ) / norm2(nn)
            pp = abs(pp)    
            qq = abs(qq)    
            
          !  print *,"Lib_SimpleSupercells::suggestPlane1 info - projection of x' on plane xy,yz,zx = ",pp
          !  print *,"Lib_SimpleSupercells::suggestPlane1 info - projection of y' on plane xy,yz,zx = ",qq
                        
            
            
            
            if (pp(1)>=max(pp(2),pp(3))) then
                !   dx1 has greatest projection on a1,a2 face
                d1 = dot_product( dx1,a_super(1:3,3) )            
                if (qq(2)>=max(qq(3),qq(1))) then
                    !   and dx2 has greatest projection on a2,a3 face       x',y',z' ~ z,x,y
                    d2 = dot_product( dx2,a_super(1:3,1) )            
                else
                    !   assume dx2 projection is along a3,a1 face           x',y',z' ~ z,-y,x                     
                    d2 = -dot_product( dx2,a_super(1:3,2) )   
                end if    
            else if (pp(2)>max(pp(3),pp(1))) then
                !   dx1 has greatest projection on a2,a3 face               
                d1 = dot_product( dx1,a_super(1:3,1) )          
                if (qq(3)>=max(qq(1),qq(2))) then
                    !   and dx2 has greatest projection on a3,a1 face       x',y',z' ~ x,y,z
                    d2 = dot_product( dx2,a_super(1:3,2) )            
                else
                    !   assume dx2 projection is along a1,a2 face           x',y',z' ~ x,-z,y
                    d2 = -dot_product( dx2,a_super(1:3,3) )   
                end if                
            else  
                !   dx1 has greatest projection on a3,a1 face
                d1 = dot_product( dx1,a_super(1:3,2) )    
                if (qq(1)>=max(qq(2),qq(3))) then
                    !   and dx2 has greatest projection on a1,a2 face       x',y',z' ~ y,z,x
                    d2 = dot_product( dx2,a_super(1:3,3) )            
                else
                    !   assume dx2 projection is along a2,a3 face           x',y',z' ~ y,-x,z
                    d2 = -dot_product( dx2,a_super(1:3,1) )   
                end if                     
            end if

            
            
        !---    can now complete the job by ensuring resultant box has same volume as initial box            
            vol = determinant3Mat(a_super)
            d3 = vol / (d1*d2)
            
            dx1 = dx1*d1
            dx2 = dx2*d2
            dx3 = dx3*d3
            
            

            return
        end subroutine suggestPlane1



            
!         subroutine completeBasis( z,x,y )
!     !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
!     !*      given the z-axis vector, complete the basis to provide a triplet x,y,z
!     !*      if on input x is set, then use this as a hint to attempt to place x along this direction
!             real(kind=real64),dimension(3),intent(in)       ::      z
!             real(kind=real64),dimension(3),intent(inout)    ::      x
!             real(kind=real64),dimension(3),intent(out)      ::      y
! 
!             real(kind=real64),dimension(3)      ::      zz        !   normalised
!             real(kind=real64)                   ::      xxxx,zzzz  !   vector lengths
!             real(kind=real64)                   ::      zdotx
! 
!         !---    check we have a normalised z as input
!             zzzz = norm2(z)
!             if (zzzz == 0) then
!                 !   can't do anything with zero input vector
!                 x = (/ 1,0,0 /)
!                 y = (/ 0,1,0 /)
!                 return
!             end if
!             zz = z/zzzz
! 
!         !---    check for a sensible hint for the x-direction
!             xxxx = norm2(x)
!             zdotx = dot_product( zz,x )
!             if ( (xxxx < 0.001d0 ).or.(abs(zdotx) >= xxxx*0.999d0) ) then
!                 if (SimpleSupercell_dbg) print *,"Lib_SimpleSupercells::completeBasis info - haven't got a good hint for x-vector"
!                 !   haven't got a good hint. Make random hint
!                 if (abs(zz(3))>max(abs(zz(1)),abs(zz(2)))) then
!                     !   z points along 3-axis
!                     x(1) = 1.0d0
!                     x(3) = 0.0d0
!                     if (abs(zz(1))>0) then
!                         x(2) = - zz(2)/zz(1)
!                     else
!                         x(2) = 0.0d0
!                     end if
!                     zdotx = dot_product( zz,x )
!                 else if (abs(zz(2))>max(abs(zz(1)),abs(zz(3)))) then
!                     !   z points along 2-axis
!                     x(3) = 1.0d0
!                     x(2) = 0.0d0
!                     if (abs(zz(3))>0) then
!                         x(1) = - zz(1)/zz(3)
!                     else
!                         x(1) = 0.0d0
!                     end if
!                     zdotx = dot_product( zz,x )
!                 else
!                     !   z points along 1-axis
!                     x(2) = 1.0d0
!                     x(1) = 0.0d0
!                     if (abs(zz(2))>0) then
!                         x(3) = - zz(3)/zz(2)
!                     else
!                         x(3) = 0.0d0
!                     end if
!                     zdotx = dot_product( zz,x )
!                 end if
!                 !print *,"suggest ",x,norm2(x),zdotx
!             end if
! 
! 
!         !---    remove projection of z on x-direction
!             x = x - zz*zdotx
!             xxxx = norm2(x)
!             x = x/xxxx
! 
!         !---    construct y-direction
!             y = cross_product( zz,x )
! 
!             if (SimpleSupercell_dbg) then
!                 print *,"Lib_SimpleSupercells::completeBasis info -"
!                 print *,"    x :",x
!                 print *,"    y :",y
!                 print *,"    z :",z
!             end if
!             return
!         end subroutine completeBasis
!         

            
            
    !---
        
        pure function getA0(this) result(a)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      return unit cell
            type(SimpleSupercell),intent(in)          ::      this
            real(kind=real64),dimension(3,3)    ::      a
            a = this%a
            return
        end function getA0
        
        pure function getiA0(this) result(ia)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      return inverse unit cell
            type(SimpleSupercell),intent(in)          ::      this
            real(kind=real64),dimension(3,3)    ::      ia
            ia = this%ia
            return
        end function getiA0
        

        pure subroutine setA0(this,a)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      set unit cell
            type(SimpleSupercell),intent(inout)          ::      this
            real(kind=real64),dimension(3,3),intent(in)    ::      a
            this%a = a
            call inverse3Mat(this%a,this%ia)
            return
        end subroutine setA0

        pure function getSuperA(this) result(a)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      return super cell
            type(SimpleSupercell),intent(in)          ::      this
            real(kind=real64),dimension(3,3)    ::      a
            a(1:3,1) = this%a(1:3,1)*this%Nx
            a(1:3,2) = this%a(1:3,2)*this%Ny
            a(1:3,3) = this%a(1:3,3)*this%Nz
            return
        end function getSuperA
        
        pure function getSuperCellSideLength0(this,i) result(a)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      return super cell side i length
            type(SimpleSupercell),intent(in)        ::      this
            integer,intent(in)                      ::      i
            real(kind=real64)                       ::      a
            if (i==1) then
                a = norm2(this%a(1:3,1))*this%Nx
            else if (i==2) then
                a = norm2(this%a(1:3,2))*this%Ny
            else
                a = norm2(this%a(1:3,3))*this%Nz                                
            end if
            return
        end function getSuperCellSideLength0
        
        elemental function getCellSideLength0(this,i) result(a)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      return unit cell side i length
            type(SimpleSupercell),intent(in)        ::      this
            integer,intent(in)                      ::      i
            real(kind=real64)                       ::      a
            a = norm2(this%a(1:3,i))
            return
        end function getCellSideLength0
        
        pure function getLatticeVector0(this,i) result(a)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      return lattice vector i
            type(SimpleSupercell),intent(in)        ::      this
            integer,intent(in)                      ::      i
            real(kind=real64),dimension(3)          ::      a
            if (i==1) then
                a = this%a(1:3,1)
            else if (i==2) then
                a = this%a(1:3,2)
            else
                a = this%a(1:3,3)                              
            end if
            return
        end function getLatticeVector0
        
        pure function getSuperVector0(this,i) result(a)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      return supercell vector i
            type(SimpleSupercell),intent(in)        ::      this
            integer,intent(in)                      ::      i
            real(kind=real64),dimension(3)          ::      a
            if (i==1) then
                a = this%a(1:3,1)*this%Nx
            else if (i==2) then
                a = this%a(1:3,2)*this%Ny
            else
                a = this%a(1:3,3)*this%Nz                                
            end if
            return
        end function getSuperVector0
        
        
    !---

        pure function getNnodes0(this) result(n)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      return the total number of cells in the supercell
            type(SimpleSupercell),intent(in)      ::      this
            integer                         ::      n
            n = this%Nx*this%Ny*this%Nz
            return
        end function getNnodes0

        
        
        
    !---        
        
        pure function getNx0(this) result(Nx)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
            type(SimpleSupercell),intent(in)      ::      this
            integer                         ::      Nx
            Nx = this%Nx
            return
        end function getNx0


        pure function getNy0(this) result(Ny)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
            type(SimpleSupercell),intent(in)      ::      this
            integer                         ::      Ny
            Ny = this%Ny
            return
        end function getNy0

        pure function getNz0(this) result(Nz)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
            type(SimpleSupercell),intent(in)      ::      this
            integer                         ::      Nz
            Nz = this%Nz
            return
        end function getNz0
 


        pure subroutine setNx(this,Nx)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
            type(SimpleSupercell),intent(inout)      ::      this
            integer,intent(in)                         ::      Nx
            this%Nx = Nx
            return
        end subroutine setNx


        pure subroutine setNy(this,Ny)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
            type(SimpleSupercell),intent(inout)      ::      this
            integer,intent(in)                         ::      Ny
            this%Ny = Ny
            return
        end subroutine setNy

        pure subroutine setNz(this,Nz)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
            type(SimpleSupercell),intent(inout)      ::      this
            integer,intent(in)                         ::      Nz
            this%Nz = Nz
            return
        end subroutine setNz
        
    !---
    
    
        pure function wrapPBC0( this,y ) result( yp )
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      wrap real space position y to one which is within the super cell
            type(SimpleSupercell),intent(in)                      ::      this
            real(kind=real64),dimension(3),intent(in)       ::      y
            real(kind=real64),dimension(3)                  ::      yp
            integer             ::      ii
            real(kind=real64),dimension(3)                  ::      xx
            
            xx(1:3) = this%ia(1:3,1)*y(1) + this%ia(1:3,2)*y(2) + this%ia(1:3,3)*y(3)
            !xx(1:3) = realSpaceToCell( this,y )
            
            ii = floor(xx(1))
            xx(1) = xx(1) - ii
            ii = mod( ii+ishft(this%Nx,8),this%Nx )
            xx(1) = (ii+xx(1))
            
            ii = floor(xx(2))
            xx(2) = xx(2) - ii
            ii = mod( ii+ishft(this%Ny,8),this%Ny )
            xx(2) = (ii+xx(2))
            
            ii = floor(xx(3))
            xx(3) = xx(3) - ii
            ii = mod( ii+ishft(this%Nz,8),this%Nz )
            xx(3) = (ii+xx(3))
            
            yp(1:3) = cellToRealSpace( this,xx )
            return
        end function wrapPBC0
        
        pure function wrapPBC1( this,x,cellspace ) result( xp )
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      wrap position y to one which is within the super cell 
            type(SimpleSupercell),intent(in)                      ::      this
            real(kind=real64),dimension(3),intent(in)       ::      x
            logical,intent(in)                              ::      cellspace
            real(kind=real64),dimension(3)                  ::      xp
            integer             ::      ii
            
            if (.not. cellspace) then
                xp = wrapPBC0( this,x ) 
                return
            end if
            ii = floor(xp(1))
            xp(1) = xp(1) - ii
            ii = mod( ii+ishft(this%Nx,8),this%Nx )
            xp(1) = (ii+xp(1))
            ii = floor(xp(2))
            xp(2) = xp(2) - ii
            ii = mod( ii+ishft(this%Ny,8),this%Ny )
            xp(2) = (ii+xp(2))
            ii = floor(xp(3))
            xp(3) = xp(3) - ii
            ii = mod( ii+ishft(this%Nz,8),this%Nz )
            xp(3) = (ii+xp(3))
            return
        end function wrapPBC1

        pure function wrapPBC2( this,c ) result( cp )
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      wrap position c to one which is within the super cell 
            type(SimpleSupercell),intent(in)                ::      this
            integer,dimension(3),intent(in)                 ::      c
            integer,dimension(3)                            ::      cp          
            cp(1) = mod( c(1)+ishft(this%Nx,8),this%Nx )
            cp(2) = mod( c(2)+ishft(this%Ny,8),this%Ny )
            cp(3) = mod( c(3)+ishft(this%Nz,8),this%Nz )
            return
        end function wrapPBC2
                
        pure function pointInSupercell0( this,x ) result(is)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      !returns true if wrapping x back into the supercell does not move it
    !*      returns true if point is in the (0,0,0) replica
            type(SimpleSupercell),intent(in)                ::      this
            real(kind=real64),dimension(3),intent(in)       ::      x
            logical                                         ::      is
    
            integer             ::      ix,iy,iz
!           real(kind=real64)   ::      y1,y2,y3
!            y1 = this%iA(1,1)*x(1) + this%iA(1,2)*x(2) + this%iA(1,3)*x(3) 
!            y2 = this%iA(2,1)*x(1) + this%iA(2,2)*x(2) + this%iA(2,3)*x(3) 
!            y3 = this%iA(3,1)*x(1) + this%iA(3,2)*x(2) + this%iA(3,3)*x(3)  
!            ix = floor(y1)
!            iy = floor(y2)
!            iz = floor(y3)
!            ix = ix - mod(ix+ishft(this%Nx,8),this%Nx )
!            iy = iy - mod(iy+ishft(this%Ny,8),this%Ny )
!            iz = iz - mod(iz+ishft(this%Nz,8),this%Nz )            
!            is = (ix*ix + iy*iy + iz*iz == 0)

            ix = floor( this%iA(1,1)*x(1) + this%iA(1,2)*x(2) + this%iA(1,3)*x(3) )     
            is = ( ix == mod(ix+ishft(this%Nx,8),this%Nx ) )
            
            iy = floor( this%iA(2,1)*x(1) + this%iA(2,2)*x(2) + this%iA(2,3)*x(3) )     
            is = is .and. ( iy == mod(iy+ishft(this%Ny,8),this%Ny ) )
            
            iz = floor( this%iA(3,1)*x(1) + this%iA(3,2)*x(2) + this%iA(3,3)*x(3) )     
            is = is .and. ( iz == mod(iz+ishft(this%Nz,8),this%Nz ) )
            
            
            



!            is =          (floor( (this%iA(1,1)*x(1) + this%iA(1,2)*x(2) + this%iA(1,3)*x(3))/this%Nx ) == 0)
!            is = is .and. (floor( (this%iA(2,1)*x(1) + this%iA(2,2)*x(2) + this%iA(2,3)*x(3))/this%Ny ) == 0)
!            is = is .and. (floor( (this%iA(3,1)*x(1) + this%iA(3,2)*x(2) + this%iA(3,3)*x(3))/this%Nz ) == 0)  
            return
        end function pointInSupercell0
            
        
        pure function linearInterpolation0( this,x,dat ) result( d )
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      given a vector field dat defined at the unit cell nodes
    !*      compute a linear interpolation of the value at x
            type(SimpleSupercell),intent(in)                    ::      this
            real(kind=real64),dimension(3),intent(in)           ::      x
            real(kind=real64),dimension(:,0:,0:,0:),intent(in)  ::      dat
            real(kind=real64),dimension(size(dat,dim=1))        ::      d
            
            integer                         ::      ix,iy,iz
            integer                         ::      jx,jy,jz
            real(kind=real64)               ::      ux,uy,uz
            real(kind=real64),dimension(3)  ::      yy
                                              
        !---    find x in cell space
            yy(1:3) = this%ia(1:3,1)*x(1) + this%ia(1:3,2)*x(2) + this%ia(1:3,3)*x(3)
            ix = floor( yy(1) )     !   note: could be outside cell range 0:Nx-1
            iy = floor( yy(2) )
            iz = floor( yy(3) )
            
        !       find distance across cell
            ux = yy(1) - ix
            uy = yy(2) - iy
            uz = yy(3) - iz
                        
        !       find left-front-bottom cell accounting for pbc
            ix = mod( ix + ishft(this%Nx,8),this%Nx ) 
            iy = mod( iy + ishft(this%Ny,8),this%Ny ) 
            iz = mod( iz + ishft(this%Nz,8),this%Nz ) 
            
        !       find right-back-top cell accounting for pbc
            jx = mod( ix+1,this%Nx )
            jy = mod( iy+1,this%Ny )
            jz = mod( iz+1,this%Nz )
            
            d(:) = dat(:,ix,iy,iz)*(1-ux)*(1-uy)*(1-uz)             &
                 + dat(:,jx,iy,iz)*(  ux)*(1-uy)*(1-uz)             &
                 + dat(:,ix,jy,iz)*(1-ux)*(  uy)*(1-uz)             &
                 + dat(:,jx,jy,iz)*(  ux)*(  uy)*(1-uz)             &
                 + dat(:,ix,iy,jz)*(1-ux)*(1-uy)*(  uz)             &
                 + dat(:,jx,iy,jz)*(  ux)*(1-uy)*(  uz)             &
                 + dat(:,ix,jy,jz)*(1-ux)*(  uy)*(  uz)             &
                 + dat(:,jx,jy,jz)*(  ux)*(  uy)*(  uz)    
                 
           return
       end function linearInterpolation0 
                
        pure function linearInterpolation2d( this,x,dat ) result( d )
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      given a vector field dat defined at the unit cell nodes
    !*      compute a linear interpolation of the value at x
    !*      special case: Nz=1
            type(SimpleSupercell),intent(in)                    ::      this
            real(kind=real64),dimension(3),intent(in)           ::      x
            real(kind=real64),dimension(:,0:,0:),intent(in)     ::      dat
            real(kind=real64),dimension(size(dat,dim=1))        ::      d
            
            integer                         ::      ix,iy 
            integer                         ::      jx,jy 
            real(kind=real64)               ::      ux,uy 
            real(kind=real64),dimension(2)  ::      yy
                                              
        !---    find x in cell space
            yy(1:2) = this%ia(1:2,1)*x(1) + this%ia(1:2,2)*x(2) + this%ia(1:2,3)*x(3)  
            ix = floor( yy(1) )     !   note: could be outside cell range 0:Nx-1
            iy = floor( yy(2) )
            
        !       find distance across cell
            ux = yy(1) - ix
            uy = yy(2) - iy
                        
        !       find left-front-bottom cell accounting for pbc
            ix = mod( ix + ishft(this%Nx,8),this%Nx ) 
            iy = mod( iy + ishft(this%Ny,8),this%Ny ) 
            
        !       find right-back-top cell accounting for pbc
            jx = mod( ix+1,this%Nx )
            jy = mod( iy+1,this%Ny )
            
            d(:) = dat(:,ix,iy)*(1-ux)*(1-uy)              &
                 + dat(:,jx,iy)*(  ux)*(1-uy)              &
                 + dat(:,ix,jy)*(1-ux)*(  uy)              &
                 + dat(:,jx,jy)*(  ux)*(  uy)               
                                                            
           return
       end function linearInterpolation2d 
                
!------

        pure function cross_product(x,y) result(z)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
            real(kind=real64),dimension(3),intent(in)       ::      x,y
            real(kind=real64),dimension(3)                  ::      z
            z(1) = x(2)*y(3) - x(3)*y(2)
            z(2) = x(3)*y(1) - x(1)*y(3)
            z(3) = x(1)*y(2) - x(2)*y(1)
            return
        end function cross_product


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
            
        
        
        
    end module Lib_SimpleSupercells
    