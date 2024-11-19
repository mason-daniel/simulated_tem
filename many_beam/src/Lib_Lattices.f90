
    module Lib_Lattices
!---^^^^^^^^^^^^^^^^^^^
!*      Provides some useful functions for working with crystal lattice types
!*      At present, does not support all Bravais lattices. Not much demand for some of them... but they could be implemented here.
!*      Note that a constructor exists to produce them if necessary.
!*
!*      note that hcp cell has default side lengths |a| = a0 / 2^1/6, |c| = a0 sqrt(8/3) / 2^1/6
!*      so that the primitive unit cell has volume 2 a0^3 
!*      a different c/a ratio can be handled by assuming an elastic strain in the z-direction or using the setCoverA() call
!*
!*      Daniel Mason, UKAEA
!*      April 2022



        use iso_fortran_env
        implicit none
        private
        
    !---
    
        public      ::      Lattice_ctor                    !   construct a Lattice option
        public      ::      delete                          !   deallocates dynamic memory       
        public      ::      report                          !   simple output of lattice properties
        public      ::      clone                           !   deep copy with allocation
     
        public      ::      getNmotif
        public      ::      getMotif        
        public      ::      getNneighbours
        public      ::      getNeighbour
        public      ::      getNeighbours
        public      ::      getSymmetry
        public      ::      getPrimitiveCell
        public      ::      getConventionalCell
        public      ::      getConventionalnMotif
        public      ::      getConventionalMotif
        public      ::      getReciprocalLatticeVectors
        public      ::      getOmega0
        
        public      ::      getLatticeType
        public      ::      getLatticeName
        public      ::      listAvailableLattices
        public      ::      getSymmetryName
        public      ::      listAvailablesymmetries
        
        public      ::      listPermittedReflections
        public      ::      permittedReflections    
        
        public      ::      getNDiffractionSpots
        public      ::      getDiffractionSpot
        public      ::      findConventionalCellFromDiffractionSpots
        
        public      ::      setCoverA
        public      ::      MillerToMillerBravais_direction
        public      ::      MillerBravaisToMiller_direction
        public      ::      MillerToMillerBravais_plane
        public      ::      MillerBravaisToMiller_plane
 
        
    !---
    
        logical,public          ::      Lattice_dbg = .false.
        
        integer,private,parameter      ::      NLATTICE        = 4
        integer,public,parameter       ::      LATTICE_CUSTOM  = 0
        integer,public,parameter       ::      LATTICE_SC      = 1
        integer,public,parameter       ::      LATTICE_BCC     = 2 
        integer,public,parameter       ::      LATTICE_FCC     = 3
        integer,public,parameter       ::      LATTICE_HCP     = 4
        character(len=6),dimension(0:NLATTICE),public,parameter       ::      LATTICE_NAME = (/ "custom","sc    ","bcc   ","fcc   ","hcp   " /)

        real(kind=real64),private,parameter                 ::  PI = 3.141592653589790d0


        real(kind=real64),dimension(3,1),private,parameter  ::  SC_MOTIF     = reshape( (/1,1,1/)*0.50d0,(/3,1/) )
        real(kind=real64),dimension(3,1),private,parameter  ::  BCC_MOTIF    = reshape( (/1,1,1/)*0.25d0,(/3,1/) )
        real(kind=real64),dimension(3,1),private,parameter  ::  FCC_MOTIF    = reshape( (/1,1,1/)*0.25d0,(/3,1/) )
        real(kind=real64),dimension(3,2),private,parameter  ::  HCP_MOTIF    = reshape( (/1,1,1 , 5,9,7 /)/12.0d0,(/3,2/) )  

        real(kind=real64),dimension(3,3),parameter,private  ::  PRIMITIVE_SC  = reshape( (/  1,0,0 , 0,1,0  , 0,0,1  /) , (/3,3/) )
        real(kind=real64),dimension(3,3),parameter,private  ::  PRIMITIVE_BCC = reshape( (/ -1,1,1 , 1,-1,1 , 1,1,-1 /) , (/3,3/) )*0.5d0
        real(kind=real64),dimension(3,3),parameter,private  ::  PRIMITIVE_FCC = reshape( (/  0,1,1 , 1,0,1  , 1,1,0  /) , (/3,3/) )*0.5d0
        real(kind=real64),dimension(3,3),parameter,private  ::  PRIMITIVE_HCP = reshape( (/  0.5d0,-sqrt(0.75d0),0.0d0 , 0.5d0,sqrt(0.75d0),0.0d0 , 0.0d0,0.0d0,sqrt(8.0d0/3)  /) , (/3,3/) )  
        
        real(kind=real64),dimension(3,3),parameter,private  ::  CONVENTIONAL_SC  = reshape( (/  1,0,0 , 0,1,0  , 0,0,1  /) , (/3,3/) )
        real(kind=real64),dimension(3,3),parameter,private  ::  CONVENTIONAL_BCC = reshape( (/  1,0,0 , 0,1,0  , 0,0,1  /) , (/3,3/) )
        real(kind=real64),dimension(3,3),parameter,private  ::  CONVENTIONAL_FCC = reshape( (/  1,0,0 , 0,1,0  , 0,0,1  /) , (/3,3/) )
!        real(kind=real64),dimension(3,3),parameter,private  ::  CONVENTIONAL_HCP = reshape( (/  1.0d0,0.0d0,0.0d0, 0.0d0,sqrt(3.0d0),0.0d0, 0.0d0,0.0d0,sqrt(8.0d0/3) /) , (/3,3/) ) 
        real(kind=real64),dimension(3,3),parameter,private  ::  CONVENTIONAL_HCP = PRIMITIVE_HCP
        
        real(kind=real64),dimension(3,3),parameter,private  ::  DIFFPATT_SC_1  = reshape( (/ 1,0,0 , 0,1,0 , 0,0,1 /) , (/3,3/) ) 
        real(kind=real64),dimension(3,3),parameter,private  ::  DIFFPATT_BCC_1 = reshape( (/ 2,0,0 , 0,2,0 , 0,0,2 /) , (/3,3/) ) 
        real(kind=real64),dimension(3,12),parameter,private ::  DIFFPATT_BCC_2 = reshape( (/ 1,1,0 , 1,0,1 , 0,1,1 , -1,1,0,-1,0,1,0,-1,1 , 1,-1,0,1,0,-1,0,1,-1 , -1,-1,0,-1,0,-1,0,-1,-1 /) , (/3,12/) ) 
        real(kind=real64),dimension(3,3),parameter,private  ::  DIFFPATT_FCC_1 = reshape( (/ 2,0,0 , 0,2,0 , 0,0,2 /) , (/3,3/) ) 
        real(kind=real64),dimension(3,8),parameter,private  ::  DIFFPATT_FCC_2 = reshape( (/ 1,1,1 , -1,1,1 , 1,-1,1 , 1,1,-1 , -1,-1,1 , -1,1,-1 , 1,-1,-1 , -1,-1,-1 /) , (/3,8/) ) 
        real(kind=real64),dimension(3,3),parameter,private  ::  DIFFPATT_HCP_1 = reshape( (/ 1,1,0 , -1,1,0 , 0,0,1 /) , (/3,3/) ) 
        real(kind=real64),dimension(3,10),parameter,private ::  DIFFPATT_HCP_2 = reshape( (/ 1,0,1 , 0,1,1 , -1,0,1,0,-1,1 , 1,-1,0,1,0,-1,0,1,-1 , -1,-1,0,-1,0,-1,0,-1,-1 /) , (/3,10/) ) 
        
        integer,private,parameter      ::      NSYMMETRY       = 2
        integer,public,parameter       ::      LATTICE_SYM_CUSTOM      = 0
        integer,public,parameter       ::      LATTICE_SYM_CUBIC       = 1
        integer,public,parameter       ::      LATTICE_SYM_HEXAGONAL   = 2
        character(len=12),dimension(0:NSYMMETRY),public,parameter       ::      SYMMETRY_NAME = (/ "custom      ","cubic       ","hexagonal   " /)


    !---
    
        type,public     ::      Lattice
            private
            integer                                         ::      latt
            integer                                         ::      sym
            integer                                         ::      nMotif                  
            real(kind=real64),dimension(:,:),pointer        ::      motif                   !   (3,this%nMotif)
            real(kind=real64),dimension(3,3)                ::      b
            real(kind=real64),dimension(3,3)                ::      u                       !   conventional cell c = u b
            
            integer,dimension(:),pointer                    ::      nNeighbours             !   (this%nMotif)
            real(kind=real64),dimension(:,:,:),pointer      ::      neighbour               !   (3,maxval(nNeighbours),this%nMotif)
        end type Lattice
        
    !---
    
        interface Lattice_ctor
            module procedure        Lattice_null
            module procedure        Lattice_ctor0
            module procedure        Lattice_ctor1

        end interface
                
        interface delete
            module procedure        delete0
        end interface
        
        interface report
            module procedure        report0
            module procedure        report1
        end interface
        
        interface       clone
            module procedure        clone0
        end interface
        
        interface   getNneighbours
            module procedure        getNneighbours0
            module procedure        getNneighbours1
        end interface
            

        interface getLatticeName
            module procedure        getLatticeName0
        end interface
        
        

        interface getOmega0
            module procedure        getOmega00
        end interface
        
        interface getMotif
            module procedure        getMotif0
        end interface
        
        interface   getConventionalCell
            module procedure        getConventionalCell0
        end interface

        interface   getReciprocalLatticeVectors
            module procedure        getReciprocalLatticeVectors0
        end interface
        
        interface   getLatticeType
            module procedure        getLatticeType0
            module procedure        getLatticeType1
        end interface

        
        interface MillerToMillerBravais_direction
            module procedure        MillerToMillerBravais_direction0
            module procedure        MillerToMillerBravais_direction1
        end interface
        
        interface MillerBravaisToMiller_direction
            module procedure        MillerBravaisToMiller_direction0
            module procedure        MillerBravaisToMiller_direction1
        end interface
        
        interface MillerToMillerBravais_plane
            module procedure        MillerToMillerBravais_plane0
            module procedure        MillerToMillerBravais_plane1
        end interface
        
        interface MillerBravaisToMiller_plane
            module procedure        MillerBravaisToMiller_plane0
            module procedure        MillerBravaisToMiller_plane1
        end interface

        
         
    contains
!---^^^^^^^^

        function Lattice_null() result(this)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ 
    !*      default object has no memory allocation
            type(Lattice)           ::      this
            this%nMotif = 0
            nullify(this%motif)     
            nullify(this%nNeighbours)
            nullify(this%neighbour)               
            this%b = PRIMITIVE_SC
            this%u = reshape( (/1,0,0,0,1,0,0,0,1/),(/3,3/) )
            return
        end function Lattice_null
                         
        function Lattice_ctor0(motif,b,c,rc) result(this)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      compute a lattice of arbitrary type given an input motif, primitive cell, and cutoff range for the near neighbours
            real(kind=real64),dimension(:,:),intent(in)     ::      motif
            real(kind=real64),dimension(3,3),intent(in)     ::      b,c
            real(kind=real64),intent(in)                    ::      rc
            type(Lattice)           ::      this
            
            integer                             ::      ix,iy,iz,ik,ii,mm,jj
            real(kind=real64),dimension(3)      ::      dx,dy 
            real(kind=real64),dimension(:,:,:),allocatable  ::      neigh_tmp
!            logical                             ::      ok
            
            this%latt = LATTICE_CUSTOM
            this%sym = LATTICE_SYM_CUSTOM
            this%b = b
            
        !---    conventional cell c = u b
            call inverse3Mat( this%b,this%u )
            this%u = matmul( c,this%u )
            
            this%nMotif = size(motif,dim=2)
            allocate(this%motif(3,this%nMotif))
            this%motif = motif
            
            allocate(this%nNeighbours(this%nMotif))
            this%nNeighbours = 0
            mm = ceiling(rc+2)            
            
            allocate( neigh_tmp(3,8*mm*mm*mm*this%nMotif,this%nMotif) )
             
            
            
            neigh_tmp = 0
            do ii = 1,this%nMotif
                jj = 0          !   will be number of neighbours for this sublattice
                do iz = -mm,mm
                    do iy = -mm,mm
                        do ix = -mm,mm
                            do ik = 1,this%nMotif
                                
                                dy = (/ix,iy,iz/) + this%motif(:,ik) - this%motif(:,ii)
                                
                                dx(1:3) = this%b(1:3,1)*dy(1) + this%b(1:3,2)*dy(2) + this%b(1:3,3)*dy(3) 
                                
                                if (norm2(dx)<rc+1.0d-8) then                       
                                    jj = jj + 1
                                    neigh_tmp(1:3,jj,ii) = dx(1:3)
!                                    print *,"motif ",ii,jj," neigh ",ik,ix,iy,iz," d ",norm2(dx)," x ",dx        
                                    !this%nNeighbours(ii) = this%nNeighbours(ii) + 1
                                    !neigh_tmp(1:3,this%nNeighbours(ii),ii) = dx(1:3)
                                    
                                    !print *,"motif ",ii,this%nNeighbours(ii)," neigh ",ik,ix,iy,iz," d ",norm2(dx)," x ",dx
                                end if
                            end do
                        end do
                    end do
                end do
                this%nNeighbours(ii) = jj
                
                
            end do
                
            mm = maxval(this%nNeighbours)
            allocate(this%neighbour(3,mm,this%nMotif))
            this%neighbour = 0
            do ii = 1,this%nMotif
                jj = this%nNeighbours(ii)
                this%neighbour(1:3,1:jj,ii) = neigh_tmp(1:3,1:jj,ii)
            end do
            deallocate(neigh_tmp)
            
            
            !this%nNeighbours = 0
            !do iz = -mm,mm
            !    do iy = -mm,mm
            !        do ix = -mm,mm
            !            do ik = 1,this%nMotif
            !                do ii = 1,this%nMotif
            !                    dx = (/ix,iy,iz/) + this%motif(:,ik) - this%motif(:,ii)
            !                    dx = matmul( this%b,dx )
            !                    if (norm2(dx)<rc+1.0d-8) then
            !                        this%nNeighbours(ii) = this%nNeighbours(ii) + 1
            !                        this%neighbour(:,this%nNeighbours(ii),ii) = dx(:)
            !                    end if
            !                end do
            !            end do
            !        end do
            !    end do
            !end do
            
!            print *,"unsorted"
!            call report(this)
            
            call sortNeighbourList(this)
            
!            print *,"sorted"
!            call report(this)
            
            return
        end function Lattice_ctor0
                       
    !---
    
        function Lattice_ctor1(name,rc_in) result(this)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      produce a named lattice type, with a cutoff for the near neighbours
            character(len=*),intent(in)                     ::      name
            real(kind=real64),intent(in),optional           ::      rc_in
            type(Lattice)           ::      this
            
            real(Kind=real64)       ::      rc
            
            rc = 1.001d0 ; if (present(rc_in)) rc = rc_in
            
            select case(getLatticeType( name ))
                case(LATTICE_BCC)
                    this = Lattice_ctor0(BCC_MOTIF,PRIMITIVE_BCC,CONVENTIONAL_BCC,rc)                                        
                    this%latt = LATTICE_BCC
                    this%sym = LATTICE_SYM_CUBIC
                case(LATTICE_FCC)
                    this = Lattice_ctor0(FCC_MOTIF,PRIMITIVE_FCC,CONVENTIONAL_FCC,rc)                    
                    this%latt = LATTICE_FCC
                    this%sym = LATTICE_SYM_CUBIC
                case(LATTICE_HCP)
                    this = Lattice_ctor0(HCP_MOTIF,PRIMITIVE_HCP,CONVENTIONAL_HCP,rc)                    
                    this%latt = LATTICE_HCP
                    this%sym = LATTICE_SYM_HEXAGONAL                    
                case default
                    this = Lattice_ctor0(SC_MOTIF,PRIMITIVE_SC,CONVENTIONAL_SC,rc)                    
                    this%latt = LATTICE_SC
                    this%sym = LATTICE_SYM_CUBIC
            end select
            
            return
        end function Lattice_ctor1
    
        
        
    !---
        
        subroutine delete0(this)
    !---^^^^^^^^^^^^^^^^^^^^^^^^
    !*      deallocate dynamic memory
            type(Lattice),intent(inout)    ::      this
            if (this%nMotif == 0) return
            deallocate(this%motif)
            deallocate(this%nNeighbours)
            deallocate(this%neighbour)
            this = Lattice_null()
            return
        end subroutine delete0
        
        
        subroutine clone0(this,that)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      deep copy with allocation this = that
            type(Lattice),intent(out)   ::      this
            type(Lattice),intent(in)   ::      that
            this = Lattice_ctor()
            this%latt = that%latt
            this%sym  = that%sym
            this%nMotif = that%nMotif
            this%b(:,:) = that%b(:,:)
            this%u(:,:) = that%u(:,:)
            allocate(this%motif(3,this%nMotif))
            this%motif(:,:) = that%motif(:,:)
            allocate(this%nNeighbours(this%nMotif))
            this%nNeighbours(:) = that%nNeighbours(:)
            allocate(this%neighbour(3,maxval(this%nNeighbours),this%nMotif))
            this%neighbour(:,:,:) = that%neighbour(:,:,:)            
            return
        end subroutine clone0
         
               
        
    !---
    
        subroutine report0(this,u,o)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      simple output. Defaults to unit 6 = screen.
    !*      optional argument o determines left hand margin    
            type(Lattice),intent(in)        ::      this
            integer,intent(in),optional     ::      u,o
            integer     ::      uu,oo!,ii,kk
            uu = 6 ; if (present(u)) uu = u
            oo = 0 ; if (present(o)) oo = o
            write(unit=uu,fmt='(3(a,i4))') repeat(" ",oo)//"Lattice ["//LATTICE_NAME(this%latt)//","//SYMMETRY_NAME(this%sym)//",nMotif=",this%nMotif,",nNeigh(max)=",maxval(this%nNeighbours),"]"
            
            ! write (unit=uu,fmt='(a)',advance="no") repeat(" ",oo+4)
            ! do ii = 1,this%nMotif
            !     write (unit=uu,fmt='(a4,a18,i2,a3)',advance="no") "","neighbours, motif ",ii," "
            ! end do
            ! write (unit=uu,fmt='(a)',advance="yes") ""
            ! do kk = 1,maxval(this%nNeighbours)
            !     write (unit=uu,fmt='(a)',advance="no") repeat(" ",oo+4)
            !     do ii = 1,this%nMotif
            !         if (kk<=this%nNeighbours(ii)) then      !   if one motif point has fewer neighbours, write a blank space
            !             write (unit=uu,fmt='(3f8.3)',advance="no") this%neighbour(:,kk,ii)
            !         else
            !             write (unit=uu,fmt='(a24)',advance="no")  ""
            !         end if                    
            !         if (ii<this%nMotif)  write (unit=uu,fmt='(a)',advance="no") " , "                
            !     end do
            !     write (unit=uu,fmt='(a)',advance="yes") ""
            ! end do
            return
        end subroutine report0
    
    
        subroutine report1(this,verbose,u,o)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      simple output. Defaults to unit 6 = screen.
    !*      optional argument o determines left hand margin    
            type(Lattice),intent(in)        ::      this
            logical,intent(in)              ::      verbose
            integer,intent(in),optional     ::      u,o
            integer     ::      uu,oo,ii,kk
            uu = 6 ; if (present(u)) uu = u
            oo = 0 ; if (present(o)) oo = o
            call report0(this,uu,oo)
            if (.not. verbose) return
            write (unit=uu,fmt='(a)',advance="no") repeat(" ",oo+4)
            do ii = 1,this%nMotif
                write (unit=uu,fmt='(a4,a18,i2,a3)',advance="no") "","neighbours, motif ",ii," "
            end do
            write (unit=uu,fmt='(a)',advance="yes") ""
            do kk = 1,maxval(this%nNeighbours)
                write (unit=uu,fmt='(a)',advance="no") repeat(" ",oo+4)
                do ii = 1,this%nMotif
                    if (kk<=this%nNeighbours(ii)) then      !   if one motif point has fewer neighbours, write a blank space
                        write (unit=uu,fmt='(3f8.3)',advance="no") this%neighbour(:,kk,ii)
                    else
                        write (unit=uu,fmt='(a24)',advance="no")  ""
                    end if                    
                    if (ii<this%nMotif)  write (unit=uu,fmt='(a)',advance="no") " , "                
                end do
                write (unit=uu,fmt='(a)',advance="yes") ""
            end do
            return
        end subroutine report1
    
    !---
        
        pure function getNmotif(this) result(n)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      returns the number of motif points in the primitive lattice.
    !*      note that the number in the conventional lattice may be different.
            type(Lattice),intent(in)    ::      this
            integer                     ::      n
            n = this%nMotif
            return
        end function getNmotif
        
        
        pure function getConventionalNmotif(this) result(n)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      returns the number of motif points in the conventional lattice.
    !*      note that the number in the primitive lattice may be different.
            type(Lattice),intent(in)    ::      this
            integer                     ::      n
            
            select case(this%latt)
                case(LATTICE_BCC)
                    n = 2                                
                case(LATTICE_FCC)
                    n = 4 
                case(LATTICE_HCP)
                    n = 2
                    !n = 4
                case default
                    n = 1
            end select
             
            return
        end function getConventionalNmotif
        
        pure function getMotif0(this,i) result(x)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      returns the motif points in the primitive unit cell
    !*      note that the motif points in the conventional cell might be different
            type(Lattice),intent(in)    ::      this
            integer,intent(in)          ::      i
            real(kind=real64),dimension(3)  ::      x
            x = this%motif(:,i)
            return
        end function getMotif0
        

        pure function getConventionalMotif(this,i,removeOffset) result(c)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      returns the motif points in the conventional cell
    !*      by default adds a [1/4,1/4,1/4] offset to atomic positions
    !*      but can optionally switch this off

            type(Lattice),intent(in)                            ::      this
            integer,intent(in)                                  ::      i
            logical,intent(in),optional                         ::      removeOffset
            real(kind=real64),dimension(3)                      ::      c
            
            if (present(removeOffset)) then
                if (removeOffset) then

                    select case(this%latt)
                        case(LATTICE_BCC)
                            select case(i)
                                case(1)
                                    c = (/ 0.00d0,0.00d0,0.00d0 /)
                                case(2)
                                    c = (/ 0.50d0,0.50d0,0.50d0 /)
                            end select   
                        case(LATTICE_FCC)
                            select case(i)
                                case(1)
                                    c = (/ 0.00d0,0.00d0,0.00d0 /)
                                case(2)
                                    c = (/ 0.00d0,0.50d0,0.50d0 /)
                                case(3)
                                    c = (/ 0.50d0,0.00d0,0.50d0 /)
                                case(4)
                                    c = (/ 0.50d0,0.50d0,0.00d0 /)
                            end select   
                        case(LATTICE_HCP)
                            select case(i)
                                case(1)
                                    c =  (/ 0.00d0,0.00d0,0.00d0 /)
                                case(2)
                                    c = (/ 2,4,3 /)/6.0d0
                            end select   
                        case default
                            c = 0.0d0
                    end select
                    return
                end if
            end if

            select case(this%latt)
                case(LATTICE_BCC)
                    select case(i)
                        case(1)
                            c = (/ 0.25d0,0.25d0,0.25d0 /)
                        case(2)
                            c = (/ 0.75d0,0.75d0,0.75d0 /)
                     end select   
                case(LATTICE_FCC)
                    select case(i)
                        case(1)
                            c = (/ 0.25d0,0.25d0,0.25d0 /)
                        case(2)
                            c = (/ 0.25d0,0.75d0,0.75d0 /)
                        case(3)
                            c = (/ 0.75d0,0.25d0,0.75d0 /)
                        case(4)
                            c = (/ 0.75d0,0.75d0,0.25d0 /)
                     end select   
                case(LATTICE_HCP)
                    select case(i)
                        case(1)
                            c = (/ 1,1,1 /)/12.0d0
                        case(2)
                            c = (/ 5,9,7 /)/12.0d0
                     end select   
                case default
                    c = 0.25d0
            end select
             
            return
        end function getConventionalMotif
        
!-------                
        
        pure function getNneighbours0(this) result(n)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      returns maximum number of neighbours within rc for any sublattice
            type(Lattice),intent(in)    ::      this
            integer                     ::      n
            n = maxval(this%Nneighbours)
            return
        end function getNneighbours0
        
        
        pure function getNneighbours1(this,i) result(n)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      returns number of neighbours within rc for sublattice i   
            type(Lattice),intent(in)    ::      this
            integer,intent(in)          ::      i
            integer                     ::      n
            n = this%Nneighbours(i)
            return
        end function getNneighbours1
        
        
        pure function getNeighbour(this,k,i) result(x)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      returns vector to neighbour k in sublattice i
            type(Lattice),intent(in)    ::      this
            integer,intent(in)          ::      i,k
            real(kind=real64),dimension(3)  ::      x
            x = this%neighbour(:,k,i)
            return
        end function getNeighbour
        
    
        pure function getNeighbours(this,i) result(x)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      returns vectors to all neighbours in sublattice i
            type(Lattice),intent(in)    ::      this
            integer,intent(in)          ::      i
            real(kind=real64),dimension(3,this%nNeighbours(i))  ::      x
            x = this%neighbour(:,1:this%nNeighbours(i),i)
            return
        end function getNeighbours
        
        
    !---    help interpretting input
    
        subroutine listAvailableLattices()
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      screen dump all known lattices
            integer     ::      ii
            print *,"Lib_Lattices::listAvailableLattices() info"
            do ii = 1,NLATTICE
                write(*,fmt='(a)') """"//trim(LATTICE_NAME(ii))//""""
            end do
            return
        end subroutine listAvailableLattices


        subroutine listAvailableSymmetries()
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      screen dump all known symmetries
            integer     ::      ii
            print *,"Lib_Lattices::listAvailableSymmetries() info"
            do ii = 1,NSYMMETRY
                write(*,fmt='(a)') """"//trim(SYMMETRY_NAME(ii))//""""
            end do
            return
        end subroutine listAvailableSymmetries

        function getSymmetryName( sym ) result(symmetryType)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      convert character string to integer code
            character(len=*),intent(in)     ::      sym
            integer                         ::      symmetryType
            integer         ::      ii
            do ii = 1,NSYMMETRY
                if (trim(sym)==trim(SYMMETRY_NAME(ii))) then
                    symmetryType = ii
                    return
                end if
            end do
            symmetryType = LATTICE_SYM_CUSTOM

            return
        end function getSymmetryName

        function getLatticeName0( this ) result(name)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      find name of lattice from type
            type(Lattice),intent(in)        ::      this
            character(len=6)                ::      name
            name = LATTICE_NAME(this%latt)
            return
        end function getLatticeName0
        
        function getLatticeType0( lattice ) result(latticeType)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      convert character string to integer code
            character(len=*),intent(in)     ::      lattice
            integer                         ::      latticeType
            integer         ::      ii
            do ii = 1,NLATTICE
                if (trim(lattice)==trim(LATTICE_NAME(ii))) then
                    latticeType = ii
                    return
                end if
            end do
            latticeType = LATTICE_CUSTOM
            return
        end function getLatticeType0

        pure function getLatticeType1( this ) result(latticeType)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      return integer code
            type(Lattice),intent(in)        ::      this
            integer                         ::      latticeType            
            latticeType = this%latt
            return
        end function getLatticeType1

!-------        
        
    !---
    
    
        subroutine listPermittedReflections( this,rho_in )
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      A helper routine - just gives the indices of permitted reflections 
    !*      ( in the positive octant )
            type(Lattice),intent(in)                ::      this
            real(kind=real64),intent(in),optional   ::      rho_in
            
            real(kind=real64)   ::      rho 
            integer             ::      mm,rr
            integer             ::      hp,kp,lp,hx,kx,lx
            logical             ::      ok
            
            rho = 10.0d0 ; if (present(rho_in)) rho = rho_in
            
            mm = ceiling( rho + 1.0d-12 )        !   how far to look in each direction
            
            print *,"Lib_Lattices::listPermittedReflections info - "
            select case(this%latt)
                case( LATTICE_SC:LATTICE_FCC )
                    do lp = 0,+mm
                        do kp = 0,+mm
                            do hp = 0,+mm
                                rr = lp*lp + kp*kp + hp*hp
                                if (rr > rho*rho + 1.0d-6) cycle
                                
                                if (.not. ((hp>=kp).and.(kp>=lp)) ) cycle                                                                         
                                
                                select case(this%latt)
                                    case( LATTICE_FCC )
                                        hx = mod( hp, 2 ) 
                                        kx = mod( kp, 2 )
                                        lx = mod( lp, 2 )
                                        ok = ( hx*hx + kx*kx + lx*lx ) - ( hx*kx + kx*lx + lx*hx ) == 0     !   true if all odd or all even.
                                    case( LATTICE_BCC )
                                        ok = mod( hp+kp+lp,2 ) == 0     !   true if sum is even
                                    case default
                                        ok = .true.
                                end select
                                                            
                                if (ok) then
                                    write (*,fmt='(3i4)') hp,kp,lp 
                                end if
                                
                            end do
                        end do
                    end do
                case( LATTICE_HCP )
                    do lp = 0,+mm
                        do kp = -mm,+mm
                            do hp = 0,+mm
                                rr = floor( (hp+kp)*(hp+kp) + (hp-kp)*(hp-kp)/3 + lp*lp*0.375d0 )
                                if (rr > rho*rho + 1.0d-6) cycle
                                
                                if (.not. (abs(hp)>=abs(kp)) ) cycle                                                                         
                                ok = .not. ( (mod( hp+2*kp,3 )==0) .and. (mod( abs(lp),2 )==1) )         !   false if h+2k is modulo 3 and l is odd   
                                                            
                                if (ok) then
                                    write (*,fmt='(3i4)') hp,kp,lp 
                                end if
                                
                            end do
                        end do
                    end do
            end select            
            return
        end subroutine listPermittedReflections


        subroutine permittedReflections( this,n,hkl,rho_in )
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      return the permitted reflections [hkl]
    !*      such that (h)^2 + (k)^2 + (l)^2 <= rho^2
    !*      given the lattice type
            type(Lattice),intent(in)                        ::      this        
            integer,intent(out)                             ::      n
            integer,dimension(:,:),allocatable,intent(out)  ::      hkl
            real(kind=real64),intent(in),optional           ::      rho_in
             
            integer             ::      mm,rr,pass
            integer             ::      hp,kp,lp,hx,kx,lx
            logical             ::      ok
            real(kind=real64)   ::      rho
            
            rho = 10.0d0 ; if (present(rho_in)) rho = rho_in
            
            mm = ceiling( rho + 1.0d-12 )        !   how far to look in each direction
            
            do pass = 1,2                       !   first pass - count. second pass - store.
            
                n = 0
                do lp = -mm,+mm
                    do kp = -mm,+mm
                        do hp = -mm,+mm
                            rr = lp*lp + kp*kp + hp*hp                        
                            select case(this%latt)
                                case( LATTICE_FCC )
                                    if (rr > rho*rho + 1.0d-6) cycle
                                    hx = mod( hp, 2 )
                                    kx = mod( kp, 2 )
                                    lx = mod( lp, 2 )
                                    ok = ( hx*hx + kx*kx + lx*lx ) - ( hx*kx + kx*lx + lx*hx ) == 0     !   true if all odd or all even.
                                case( LATTICE_BCC )                                 
                                    if (rr > rho*rho + 1.0d-6) cycle
                                    ok = mod( hp+kp+lp,2 ) == 0     !   true if sum is even
                                case( LATTICE_HCP )
                                    rr = floor( (hp+kp)*(hp+kp) + (hp-kp)*(hp-kp)/3 + lp*lp*0.375d0 )
                                    if (rr > rho*rho + 1.0d-6) cycle                                                        
                                    ok = .not. ( (mod( hp+2*kp,3 )==0) .and. (mod( abs(lp),2 )==1) )         !   false if h+2k is modulo 3 and l is odd            
                                case default
                                    if (rr > rho*rho + 1.0d-6) cycle
                                    ok = .true.
                            end select
                                                        
                            if (ok) then
                                n = n + 1 
                                if (pass==2) hkl(:,n) = (/ hp,kp,lp /)
                            end if
                            
                        end do
                    end do
                end do
                if (pass==1) allocate(hkl(3,n))
                    
            end do        
            return
        end subroutine permittedReflections
                                  
    !---
        
        subroutine setCoverA(this,covera)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      set the c/a ratio for the lattice. Intended for hcp ( ie set with covera = sqrt(8/3) 
            type(Lattice),intent(inout)     ::      this
            real(kind=real64),intent(in)    ::      covera
            real(kind=real64)       ::      oldcovera
            real(kind=real64),dimension(3,3)        ::      cc,iu
            
!            print *,"Lib_Lattices::setCoverA info - old primitive cell"
!            print *,this%b(1,:)
!            print *,this%b(2,:)
!            print *,this%b(3,:)
            
        !---    expect conventional cell c = u b
            cc = matmul( this%u,this%b )            !   actual conventional cell
            
!            print *,"Lib_Lattices::setCoverA info - old conventional cell"
!            print *,cc(1,:)
!            print *,cc(2,:)
!            print *,cc(3,:)
            
            
        !---    find previous c/a
            oldcovera = norm2(cc(:,3))/norm2(cc(:,1))
!            print *,"old c/a ",oldcovera," new c/a ",covera            
            
                        
                        
                        
        !---    modify conventional cell by adjusting c
            cc(:,3) = cc(:,3) * covera/oldcovera
            
!            print *,"Lib_Lattices::setCoverA info - new conventional cell"
!            print *,cc(1,:)
!            print *,cc(2,:)
!            print *,cc(3,:)
            
            
        !---    adjust primitive cell accordingly   
            call inverse3Mat( this%u,iu )
            this%b = matmul( iu,cc )            !   new primitive cell
            
            
            !print *,"Lib_Lattices::setCoverA info - old primitive cell"
            !print *,this%b(1,:)
            !print *,this%b(2,:)
            !print *,this%b(3,:)
            
            
            
            this%neighbour(3,:,:) = this%neighbour(3,:,:) * covera/oldcovera
            
            return
        end subroutine setCoverA
            
        
        pure function MillerBravaisToMiller_direction0( uvtw ) result( uvw )
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      HCP specific routine
    !       The HCP primitivelattice vectors 
    !       PRIMITIVE_HCP = reshape( (/  0.5d0,-sqrt(0.75d0),0.0d0 , 0.5d0,sqrt(0.75d0),0.0d0 , 0.0d0,0.0d0,sqrt(8.0d0/3)  /) , (/3,3/) ) 
    !           a1 =  (/ 0.5d0,-sqrt(0.75d0),0.0d0 /)
    !           a2 =  (/ 0.5d0,+sqrt(0.75d0),0.0d0 /)
    !           a3 =  - a1 - a2  =  (/ 1.0d0,0.0d0,0.0d0 /)
    !           c  =  (/ 0.0d0,0.0d0,1.0d0 /)
    !
    !       so to get a 3-vector out of 4 miller indices 
    !           x = M [uvtw] = A [uvw]
    !       with
    !           M = (   a1  a2  a3      0   )     A = (   a1  a2  a3  )
    !               (   |   |   |       0   )         (   |   |   |   )
    !               (   |   |   |       1   )         (   |   |   |   )
    !    
    !       so 
    !           [uvw] = A^-1 M [ uvtw ]
    !
    !                 = ( 1  0 -1  0 )
    !                   ( 0  1 -1  0 )
    !                   ( 0  0  0  1 )
     
    
            integer,dimension(4),intent(in)   ::      uvtw
            integer,dimension(3)              ::      uvw
 
            uvw(1) = uvtw(1) - uvtw(3)     
            uvw(2) = uvtw(2) - uvtw(3)     
            uvw(3) = uvtw(4)     
            
            
            return
        end function MillerBravaisToMiller_direction0                     
                    
        
        
        pure function MillerToMillerBravais_direction0( uvw ) result( uvtw )
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      HCP specific routine
    !       The HCP primitivelattice vectors 
    !       PRIMITIVE_HCP = reshape( (/  0.5d0,-sqrt(0.75d0),0.0d0 , 0.5d0,sqrt(0.75d0),0.0d0 , 0.0d0,0.0d0,sqrt(8.0d0/3)  /) , (/3,3/) ) 
    !           a1 =  (/ 0.5d0,-sqrt(0.75d0),0.0d0 /)
    !           a2 =  (/ 0.5d0,+sqrt(0.75d0),0.0d0 /)
    !           a3 =  - a1 - a2  =  (/ 1.0d0,0.0d0,0.0d0 /)
    !           c  =  (/ 0.0d0,0.0d0,1.0d0 /)
    !
    !       so to get a 3-vector out of 4 miller indices 
    !           x = M [uvtw] = A [uvw]
    !       with
    !           M = (   a1  a2  a3      0   )     A = (   a1  a2  a3  )
    !               (   |   |   |       0   )         (   |   |   |   )
    !               (   |   |   |       1   )         (   |   |   |   )
    !    
    !       so 
    !           [uvw] = A^-1 M [ uvtw ]
    !
    !                 = ( 1  0 -1  0 )
    !                   ( 0  1 -1  0 )
    !                   ( 0  0  0  1 )
     
    
            integer,dimension(3),intent(in)   ::      uvw
            integer,dimension(4)              ::      uvtw
 
            uvtw(1) = ( 2*uvw(1) - uvw(2))/3
            uvtw(2) = (-uvw(1) + 2*uvw(2))/3
            uvtw(3) = (-uvw(1)   - uvw(2))/3
            uvtw(4) = uvw(3)     
            
            
            return
        end function MillerToMillerBravais_direction0                     
                    
        pure function MillerBravaisToMiller_plane0( hkil ) result( hkl )
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      HCP specific routine
    !
    !       if a plane is defined by  d = [hkl].[uvw]
    !       then  
    !           d = [hkl]. ( 1  0 -1  0 ) [uvtw]
    !                      ( 0  1 -1  0 )
    !                      ( 0  0  0  1 )
    !        so
    !           [hkil] = (  1  0  0 ) [hkl]
    !                    (  0  1  0 )
    !                    ( -1 -1  0 )
    !                    (  0  0  1 )
    !       left-inverse to give
    !                                    
    !           [hkl]  = (  1  0  0  0 ) [hkil]
    !                    (  0  1  0  0 )
    !                    (  0  0  0  1 )
    
            integer,dimension(4),intent(in)   ::      hkil
            integer,dimension(3)              ::      hkl
 
            hkl(1) = hkil(1) 
            hkl(2) = hkil(2)  
            hkl(3) = hkil(4)     
                        
            return
        end function MillerBravaisToMiller_plane0                     
                    
        
        
        pure function MillerToMillerBravais_plane0( hkl ) result( hkil )
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      HCP specific routine
    
            integer,dimension(3),intent(in)   ::      hkl
            integer,dimension(4)              ::      hkil
 
            hkil(1) = hkl(1)
            hkil(2) = hkl(2)
            hkil(3) = (-hkl(1) - hkl(2))
            hkil(4) = hkl(3)     
            
            
            return
        end function MillerToMillerBravais_plane0                     
                    
        
        
        pure function MillerBravaisToMiller_direction1( uvtw ) result( uvw )
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      HCP specific routine
    !       The HCP primitivelattice vectors 
    !       PRIMITIVE_HCP = reshape( (/  0.5d0,-sqrt(0.75d0),0.0d0 , 0.5d0,sqrt(0.75d0),0.0d0 , 0.0d0,0.0d0,sqrt(8.0d0/3)  /) , (/3,3/) ) 
    !           a1 =  (/ 0.5d0,-sqrt(0.75d0),0.0d0 /)
    !           a2 =  (/ 0.5d0,+sqrt(0.75d0),0.0d0 /)
    !           a3 =  - a1 - a2  =  (/ 1.0d0,0.0d0,0.0d0 /)
    !           c  =  (/ 0.0d0,0.0d0,1.0d0 /)
    !
    !       so to get a 3-vector out of 4 miller indices 
    !           x = M [uvtw] = A [uvw]
    !       with
    !           M = (   a1  a2  a3      0   )     A = (   a1  a2  a3  )
    !               (   |   |   |       0   )         (   |   |   |   )
    !               (   |   |   |       1   )         (   |   |   |   )
    !    
    !       so 
    !           [uvw] = A^-1 M [ uvtw ]
    !
    !                 = ( 1  0 -1  0 )
    !                   ( 0  1 -1  0 )
    !                   ( 0  0  0  1 )
     
    
            real(kind=real64),dimension(4),intent(in)   ::      uvtw
            real(kind=real64),dimension(3)              ::      uvw
 
            uvw(1) = uvtw(1) - uvtw(3)     
            uvw(2) = uvtw(2) - uvtw(3)     
            uvw(3) = uvtw(4)     
            
            
            return
        end function MillerBravaisToMiller_direction1                     
                    
        
        
        pure function MillerToMillerBravais_direction1( uvw ) result( uvtw )
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      HCP specific routine
    !       The HCP primitivelattice vectors 
    !       PRIMITIVE_HCP = reshape( (/  0.5d0,-sqrt(0.75d0),0.0d0 , 0.5d0,sqrt(0.75d0),0.0d0 , 0.0d0,0.0d0,sqrt(8.0d0/3)  /) , (/3,3/) ) 
    !           a1 =  (/ 0.5d0,-sqrt(0.75d0),0.0d0 /)
    !           a2 =  (/ 0.5d0,+sqrt(0.75d0),0.0d0 /)
    !           a3 =  - a1 - a2  =  (/ 1.0d0,0.0d0,0.0d0 /)
    !           c  =  (/ 0.0d0,0.0d0,1.0d0 /)
    !
    !       so to get a 3-vector out of 4 miller indices 
    !           x = M [uvtw] = A [uvw]
    !       with
    !           M = (   a1  a2  a3      0   )     A = (   a1  a2  a3  )
    !               (   |   |   |       0   )         (   |   |   |   )
    !               (   |   |   |       1   )         (   |   |   |   )
    !    
    !       so 
    !           [uvw] = A^-1 M [ uvtw ]
    !
    !                 = ( 1  0 -1  0 )
    !                   ( 0  1 -1  0 )
    !                   ( 0  0  0  1 )
     
    
            real(kind=real64),dimension(3),intent(in)   ::      uvw
            real(kind=real64),dimension(4)              ::      uvtw
 
            uvtw(1) = ( 2*uvw(1) - uvw(2))/3
            uvtw(2) = (-uvw(1) + 2*uvw(2))/3
            uvtw(3) = (-uvw(1)   - uvw(2))/3
            uvtw(4) = uvw(3)     
            
            
            return
        end function MillerToMillerBravais_direction1                     
                    
        pure function MillerBravaisToMiller_plane1( hkil ) result( hkl )
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      HCP specific routine
    !
    !       if a plane is defined by  d = [hkl].[uvw]
    !       then  
    !           d = [hkl]. ( 1  0 -1  0 ) [uvtw]
    !                      ( 0  1 -1  0 )
    !                      ( 0  0  0  1 )
    !        so
    !           [hkil] = (  1  0  0 ) [hkl]
    !                    (  0  1  0 )
    !                    ( -1 -1  0 )
    !                    (  0  0  1 )
    !       left-inverse to give
    !                                    
    !           [hkl]  = (  1  0  0  0 ) [hkil]
    !                    (  0  1  0  0 )
    !                    (  0  0  0  1 )
    
            real(kind=real64),dimension(4),intent(in)   ::      hkil
            real(kind=real64),dimension(3)              ::      hkl
 
            hkl(1) = hkil(1) 
            hkl(2) = hkil(2)  
            hkl(3) = hkil(4)     
                        
            return
        end function MillerBravaisToMiller_plane1                     
                    
        
        
        pure function MillerToMillerBravais_plane1( hkl ) result( hkil )
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      HCP specific routine
    
            real(kind=real64),dimension(3),intent(in)   ::      hkl
            real(kind=real64),dimension(4)              ::      hkil
 
            hkil(1) = hkl(1)
            hkil(2) = hkl(2)
            hkil(3) = (-hkl(1) - hkl(2))
            hkil(4) = hkl(3)     
            
            
            return
        end function MillerToMillerBravais_plane1                     
                    
        
        
        
        pure function getOmega00(this) result(Omega0)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      returns volume _per motif point_ 
    !*      this is generally a volume per atom
            type(Lattice),intent(in)    ::      this
            real(kind=real64)           ::      Omega0
            
            Omega0 = this%b(1,1)*( this%b(2,2)*this%b(3,3) - this%b(2,3)*this%b(3,2) ) &
                   + this%b(1,2)*( this%b(2,3)*this%b(3,1) - this%b(2,1)*this%b(3,3) ) &
                   + this%b(1,3)*( this%b(2,1)*this%b(3,2) - this%b(2,2)*this%b(3,1) ) &
                   + this%b(2,1)*( this%b(3,2)*this%b(1,3) - this%b(3,3)*this%b(1,2) ) &
                   + this%b(2,2)*( this%b(3,3)*this%b(1,1) - this%b(3,1)*this%b(1,3) ) &
                   + this%b(2,3)*( this%b(3,1)*this%b(1,2) - this%b(3,2)*this%b(1,1) ) &
                   + this%b(3,1)*( this%b(1,2)*this%b(2,3) - this%b(1,3)*this%b(2,2) ) &
                   + this%b(3,2)*( this%b(1,3)*this%b(2,1) - this%b(1,1)*this%b(2,3) ) &
                   + this%b(3,3)*( this%b(1,1)*this%b(2,2) - this%b(1,2)*this%b(2,1) )
            Omega0 = Omega0 / (3*this%nMotif)
       
            return
        end function getOmega00
        
        pure function getConventionalCell0(this) result(c)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ 
    !*      return a conventional unit cell vectors as a 3x3 matrix c
    !*      where c_ij is the ith Cartesian component of the jth vector
            type(Lattice),intent(in)                    ::      this
             
            real(kind=real64),dimension(3,3)            ::      c
            
            c = matmul( this%u,this%b )
             
            return
        end function getConventionalCell0

        pure function getReciprocalLatticeVectors0(this) result(B)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      return reciprocal lattice vectors as a 3x3 matrix B
    !*      where B_ij is the ith Cartesian component of the jth vector
    !*      note: not scaled by lattice parameter
            type(Lattice),intent(in)                    ::      this             
            real(kind=real64),dimension(3,3)            ::      B
            real(kind=real64),dimension(3,3)            ::      C
            C = getConventionalCell0(this)
            call inverse3Mat(transpose(C),B) 
            B = 2*PI*B  
            return
        end function getReciprocalLatticeVectors0

        pure function getPrimitiveCell( this ) result(b)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      returns the primitive unit cell vectors as a 3x3 matrix b
    !*      where b_ij is the ith Cartesian component of the jth vector    
            type(Lattice),intent(in)            ::      this
            real(kind=real64),dimension(3,3)    ::      b
            b = this%b
            return
        end function getPrimitiveCell


        pure function getSymmetry( this ) result(symmetryType)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      returns the integer code for the symmetry 
            type(Lattice),intent(in)        ::      this
            integer                         ::      symmetryType
            symmetryType = this%sym

            return
        end function getSymmetry

    !---
        
        pure function getNDiffractionSpots(this,longlist) result(n)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      returns the number of diffraction spots to consider.
    !*      optionally adds more as a longer list
            type(Lattice),intent(in)    ::      this
            logical,intent(in),optional ::      longlist            
            integer                     ::      n
            
            select case(this%latt)
                case(LATTICE_BCC)
                    n = size( DIFFPATT_BCC_1,dim=2 )                        
                case(LATTICE_FCC)
                    n = size( DIFFPATT_FCC_1,dim=2 )  
                case(LATTICE_HCP)
                    n = size( DIFFPATT_HCP_1,dim=2 )  
                case default
                    n = size( DIFFPATT_SC_1,dim=2 )      
            end select
            
            if (present(longlist)) then
                if (longlist) then                
                    select case(this%latt)
                        case(LATTICE_BCC)
                            n = n + size( DIFFPATT_BCC_2,dim=2 )                        
                        case(LATTICE_FCC)
                            n = n + size( DIFFPATT_FCC_2,dim=2 )  
                        case(LATTICE_HCP)
                            n = n + size( DIFFPATT_HCP_2,dim=2 )  
                        case default
                            n = n + size( DIFFPATT_BCC_2,dim=2 )        
                    end select
                end if
            end if
             
            return
        end function getNDiffractionSpots
        
        
        pure function getDiffractionSpot(this,i) result(x)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      returns a g-vector in the conventional unit cell basis 
    !*      eg g=[200] rather than g=2 pi/a0 [200] 
            type(Lattice),intent(in)        ::      this
            integer,intent(in)              ::      i           
            real(kind=real64),dimension(3)  ::      x
            integer             ::      nn
            nn = getNDiffractionSpots(this)
            select case(this%latt)
                case(LATTICE_BCC)
                    if (i <= nn ) then
                        x(1:3) = DIFFPATT_BCC_1(1:3,i)
                    else
                        x(1:3) = DIFFPATT_BCC_2(1:3,i-nn)
                    end if
                case(LATTICE_FCC)
                    if (i <= nn ) then
                        x(1:3) = DIFFPATT_FCC_1(1:3,i)
                    else
                        x(1:3) = DIFFPATT_FCC_2(1:3,i-nn)
                    end if
                case(LATTICE_HCP)
                    if (i <= nn ) then
                        x(1:3) = DIFFPATT_HCP_1(1:3,i)
                    else
                        x(1:3) = DIFFPATT_HCP_2(1:3,i-nn)
                    end if
                case default
                    if (i <= nn ) then
                        x(1:3) = DIFFPATT_SC_1(1:3,i)
                    else
                        x(1:3) = DIFFPATT_BCC_2(1:3,i-nn)
                    end if                        
            end select
            
            return
        end function getDiffractionSpot
        
        
        subroutine findConventionalCellFromDiffractionSpots( this,q , a)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      given diffraction spot positions eg q = 2 pi/a0 [ 200 ]
    !*      find the best conventional cell eg a0 [ 100 010 001 ]
    !*      note that this conventional cell will probably be rotated and strained
            type(Lattice),intent(in)                        ::      this
            real(kind=real64),dimension(:,:),intent(in)     ::      q
            real(kind=real64),dimension(3,3),intent(out)    ::      a
            
            real(kind=real64),dimension(3,3)        ::  q0,iq
            integer         ::      nn
            
        !---    find the expected positions of the diffraction spots
            do nn = 1,3
                q0(1:3,nn) = getDiffractionSpot(this,nn)
            end do
            
            call inverse3Mat(q,iq) 
            
            
            a = matmul( q0,iq )
            
            a = (2*PI)*transpose(a)
            
            nn = size(q,dim=2)              !   number of spots returned in diffraction pattern
            
            if (nn > 3) then                
                !   have overfitted data required to recover unit cell
                !   start by finding an approximate ( but probably very good ) estimate
            end if            
            
            return
        end subroutine findConventionalCellFromDiffractionSpots
        
        
        pure subroutine inverse3Mat(M,N)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      utility monadic operator
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
    !*      utility returns the determinant of M
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

!-------

        pure subroutine sortNeighbourList(this)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      helper for ease of debugging output
    !*      sorts the neighbour list into length order.
            type(Lattice),intent(inout)         ::      this
            
            integer             ::      ii,mm
            logical             ::      ok
            real(kind=real64)   ::      dd
            real(kind=real64),dimension(3)  ::      neigh_tmp 
            
            
            do mm = 1,this%nMotif
            
                do  !   bubble sort
                    ok = .true.
                    do ii = 1,this%nNeighbours(mm)-1
                        neigh_tmp(1:3) = this%neighbour(1:3,ii,mm)
                        dd = norm2(neigh_tmp)
                        if ( (dd > norm2(this%neighbour(1:3,ii+1,mm))+1.0d-8) ) then
                            this%neighbour(1:3,ii,mm) = this%neighbour(1:3,ii+1,mm)
                            this%neighbour(1:3,ii+1,mm) = neigh_tmp(1:3)
                            ok = .false.
                        end if
                    end do      
                    if (ok) exit         
                end do
            
            
            end do
                    
            return
        end subroutine sortNeighbourList 
            
        
        
    end module Lib_Lattices
    
    