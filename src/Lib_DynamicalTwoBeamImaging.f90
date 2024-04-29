
    module Lib_DynamicalTwoBeamImaging
!---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
!*      produce a simulated TEM image from atom positions
!*      using the two beam columnar approximation
!       MPI version 06/07/23   
!*
        use iso_fortran_env
        use mpi_f08 

        use Lib_RelativisticElectrons
        use Lib_SimpleProgressBar        
        use Lib_LinkCell3D
        use Lib_Lattices
        use Lib_SimpleSupercells
        use Lib_ComplexSupercells
        use Lib_DiffractionConditions
        use Lib_Quaternions
        !use Lib_DataDoubler
        !use Lib_Voxelise   
        use Lib_XYZFiles                    !   because findGdotRphaseField has the option to dump atom positions tilted to diffraction conditions
        use Lib_FactoriseParallel
        use Lib_EvenlySpacedPointsInCircle
        
        implicit none
        private
        
!          
    !---
    
        public      ::      DynamicalTwoBeamImaging_ctor
        public      ::      delete
        public      ::      report
     
        public      ::      setDeviationParameter
        public      ::      getDeviationParameter
        public      ::      setIntegrationPoints
        public      ::      setMicroscopeParameters
        public      ::      setExtinctionDistances
        public      ::      setConventionalUnitCell
        public      ::      setVoxelSupercellSize
        public      ::      setUseDensityField
        public      ::      findImageParameters
        public      ::      computeD2BI_img
        public      ::      setSupercellTilt
        public      ::      setTestMode
        public      ::      getRotation
        public      ::      opXyzFile
        
        
        
        
        
        
    !---
    
        logical,public          ::      DynamicalTwoBeamImaging_dbg = .false.
        real(kind=real64),private,parameter                 ::      PI = 3.141592653589790d0  
        integer(kind=int64),private,parameter               ::      BADF00D = int( z'BADF00D',kind=int64 )        
        real(kind=real64),private,parameter                 ::      SG_UNSET = transfer( (BADF00D+ishft(BADF00D,32_int64)),1.0d0 )
                                                                              
    !---
    
        type,public     ::      DynamicalTwoBeamImaging
            private
            real(kind=real64),dimension(3)              ::      k               !   direction of electron beam ( unnormalised ) 
            real(kind=real64),dimension(3)              ::      g               !   g vector expressed in reciprocal vectors ( eg [110] )
            !real(kind=real64),dimension(3)              ::      kk              !   electron beam correct units
            !real(kind=real64),dimension(3)              ::      gg              !   g vector correct units
            real(kind=real64)                           ::      a0              !   characteristic lengthscale
            type(SimpleSupercell)                       ::      super           !   oriented supercell
            type(Lattice)                               ::      latt            !   lattice 
            real(kind=real64)                           ::      sg              !   deviation parameter
!             real(kind=real64)                           ::      sg_sig          !   deviation parameter broadening
!             integer                                     ::      sg_n            !   number of deviation parameter samples to take
            integer                                     ::      mxy,mzz         !   number of subdivisions to take per unit cell axis
            real(kind=real64)                           ::      xi0,xig         !   extinction parameters 
            real(kind=real64)                           ::      E               !   electron accelerating voltage (keV)
            logical                                     ::      useDensityField   !   should we use void calc?
            type(SimpleSupercell)                       ::      gksuper         !   supercell oriented with electron beam parameters
            real(kind=real64),dimension(3)              ::      delta           !   offset required by tilting of k vector to accommodate sg       
            real(kind=real64)                           ::      foilThickness   !   desired foil thickness
            real(kind=real64),dimension(3)              ::      voxelSupercellDimensions    !   size of the voxel supercell 
            integer                                     ::      nuvw            !   number of periodic replicas used to generate crystallite for imaging
            integer,dimension(3,125)                    ::      uvw             !   which periodic replicas are required to generate crystallite for imaging
            real(kind=real64),dimension(3,3)            ::      a_cell          !   conventional unit cell vectors
            real(kind=real64),dimension(3,3)            ::      U               !   rotation matrix used to align zone axis with (0,0,1)                             
            real(kind=real64),dimension(3,3)            ::      R               !   rotation matrix used to tilt into good diffraction conditions                    
            logical                                     ::      test            !   put in test mode - compute rotations, output xyz, but do not do any of the heavy voxelisation work or produce an image.
        end type DynamicalTwoBeamImaging        
        
    !---
    
        interface DynamicalTwoBeamImaging_ctor
            module procedure    DynamicalTwoBeamImaging_null
            module procedure    DynamicalTwoBeamImaging_ctor0
        end interface
                
        interface delete
            module procedure    delete0
        end interface
        
        interface report
            module procedure    report0
        end interface
        
        interface   setSupercellTilt
            module procedure    setSupercellTilt0
            module procedure    setSupercellTilt1
            module procedure    setSupercellTilt1a
            module procedure    setSupercellTilt2
        end interface
        
        
        interface   setVoxelSupercellSize
            module procedure    setVoxelSupercellSize2
        end interface
        
        interface   setUseDensityField
!            module procedure    setUseDensityField0
!             module procedure    setUseDensityField1
             module procedure    setUseDensityField2
        end interface        
        
    contains
!---^^^^^^^^

        function DynamicalTwoBeamImaging_null() result(this)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
            type(DynamicalTwoBeamImaging)           ::      this
            this%super = SimpleSupercell_ctor()
            this%latt = Lattice_ctor()
            this%gksuper = SimpleSupercell_ctor()
            this%a0 = 1.0d0
!            nullify(this%rho)
            call setIntegrationPoints( this,0,0 )
            call setMicroscopeParameters( this,(/0.0d0,0.0d0,1.0d0/),(/1.0d0,1.0d0,0.0d0/),150.0d0 )
            call setExtinctionDistances( this, 103.9d0,207.5d0 )
            call setUseDensityField(this , .true.)
            !this%gg = (/ sqrt(0.5d0),sqrt(0.5d0),0.0d0 /)
            !this%kk = (/ 0,0,1 /)
            this%voxelSupercellDimensions(1:3) = 0.0d0     !    this is code for set it from input supercell
            call setConventionalUnitCell( this, this%a0 * reshape( (/1,0,0 , 0,1,0 , 0,0,1/),(/3,3/) ) )
            this%sg = SG_UNSET
!             this%sg_n = 0
!             this%sg_sig = 0.0d0
            this%delta = 0.0d0
            this%nuvw = 1
            this%uvw(1:3,1) = (/ 0,0,0 /)
            this%R(1:3,1:3) = reshape( (/1,0,0 , 0,1,0 , 0,0,1/),(/3,3/) )
            this%U(1:3,1:3) = reshape( (/1,0,0 , 0,1,0 , 0,0,1/),(/3,3/) )
            call setTestMode(this,.true.)                
            return
        end function DynamicalTwoBeamImaging_null
                         
        function DynamicalTwoBeamImaging_ctor0(super,latt,a0 , k,g,E, mxy,mzz) result(this)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
            type(DynamicalTwoBeamImaging)           ::      this
            type(SimpleSupercell),intent(in)        ::      super
            type(Lattice),intent(in)                ::      latt
            real(kind=real64),intent(in)            ::      a0
            real(kind=real64),dimension(3),intent(in)       ::      k,g
            real(kind=real64),intent(in)                    ::      E
            integer,intent(in)                              ::      mxy,mzz
             
            this = DynamicalTwoBeamImaging_null()
            this%super = super
            this%latt = latt
            this%a0 = a0    
            
            call setTestMode(this,.false.)                
            call setMicroscopeParameters( this,k,g,E )
            call setIntegrationPoints( this,mxy,mzz )
            
            return
        end function DynamicalTwoBeamImaging_ctor0
                       
        subroutine setTestMode(this,test)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
            type(DynamicalTwoBeamImaging),intent(inout)     ::      this
            logical,intent(in)                              ::      test                
            this%test = test
            return
        end subroutine setTestMode                                                    
            
    !---
    
        subroutine delete0(this)
    !---^^^^^^^^^^^^^^^^^^^^^^^^
    !       note: do not deallocate  or latt, the memory for this is stored elsewhere, I just have a copy of the pointers.
            type(DynamicalTwoBeamImaging),intent(inout)    ::      this
            if (this%mxy == 0) return
!            if (this%useDensityField) deallocate(this%rho)
            this = DynamicalTwoBeamImaging_null()
            return
        end subroutine delete0
        
    !---
    
        subroutine report0(this,u,o)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^
            type(DynamicalTwoBeamImaging),intent(in)    ::      this
            integer,intent(in),optional     ::      u,o
            integer     ::      uu,oo
            
            real(kind=real64),dimension(3)          ::      gg
            real(kind=real64),dimension(3,3)        ::      a_cell,b_cell,RU
            real(kind=real64)                       ::      modk
            
            uu = 6 ; if (present(u)) uu = u
            oo = 0 ; if (present(o)) oo = o
            write(unit=uu,fmt='(a)') repeat(" ",oo)//"DynamicalTwoBeamImaging"        
            write(unit=uu,fmt='(a)') repeat(" ",oo+2)//"Supercell associated with input data"
            call report(this%super,uu,oo+4)        
            
            write(unit=uu,fmt='(a)') repeat(" ",oo+2)//"Voxels oriented with electron beam"
            write(unit=uu,fmt='(a)') repeat(" ",oo+2)//"    ( gives the x- y- spacing of the output image pixels and the z- depth resolution of the integration )"
            call report(this%gksuper,uu,oo+4)
            write(unit=uu,fmt='(4(a,f14.6))') repeat(" ",oo+2)//"extinction distances xi0,xig = ",this%xi0,",",this%xig," (1/A)"
            write(unit=uu,fmt='(4(a,3f8.2))') repeat(" ",oo+2)//"k-vector [uvw]  = ",this%k 
            write(unit=uu,fmt='(4(a,3f8.2))') repeat(" ",oo+2)//"g-vector [hkl]  = ",this%g
            write(unit=uu,fmt='(4(a,f14.6))') repeat(" ",oo+2)//"dev param sg    = ",this%sg," (1/A)"
            
            
            write(unit=uu,fmt='(a,f16.6,a)')  repeat(" ",oo+2)//" v              = ",velocity( this%E*1000 )," (m/s)"
            
            modk = 2*PI/wavelength( this%E*1000 )       !   wave vector of electron (A^-1)            
            
            
            write(unit=uu,fmt='(a,3f14.6,a)') repeat(" ",oo+2)//" aperture (lab) = ",getAperture(this)," (1/A)"
            write(unit=uu,fmt='(a,3f14.6,a)') repeat(" ",oo+2)//" beam dir (lab) = ",modk*(/0,0,1/)," (1/A)"
            
            a_cell = matmul( this%U,this%a_cell )
            b_cell = getReciprocalLatticeVectors( a_cell )      
            gg(1:3) = reflection( b_cell,this%g )
            write(unit=uu,fmt='(a,3f14.6,a)') repeat(" ",oo+2)//" g-vec (xtal)   = ",gg," (1/A)"
            
            RU = matmul( this%R,this%U )
            write(unit=uu,fmt='(a,3f14.6,a)') repeat(" ",oo+2)//" zone ax (xtal) = ",matmul(transpose(RU),(/0.0d0,0.0d0,modk/))," (1/A)"
            
            
           ! write(unit=uu,fmt='(a,f16.6,a)')  repeat(" ",oo+2)//" s_g           = ",this%sg," (A)"
            write(unit=uu,fmt='(a,f16.6,a)')  repeat(" ",oo+2)//" w = xi_g s_g   = ",this%sg*this%xig
                
            
            
            
            
            
            return
        end subroutine report0
    
    !---
        
        subroutine getRotation( this,U,R )
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
            type(DynamicalTwoBeamImaging),intent(in)        ::      this
            real(kind=real64),dimension(3,3),intent(out)    ::      U,R
            
            U = this%U
            R = this%R
            
            return
        end subroutine getRotation
    
    
        subroutine setDeviationParameter( this,sg )
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ 
            type(DynamicalTwoBeamImaging),intent(inout)     ::      this
            real(kind=real64),intent(in)                    ::      sg
            this%sg = sg

            return
        end subroutine setDeviationParameter
!     
!         subroutine setDeviationParameter( this,sg,sg_n,sg_sig )
!     !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ 
!             type(DynamicalTwoBeamImaging),intent(inout)     ::      this
!             real(kind=real64),intent(in)                    ::      sg
!             integer,intent(in),optional                     ::      sg_n
!             real(kind=real64),intent(in),optional           ::      sg_sig
!             this%sg = sg
! !             if (present(sg_n)) then
! !                 this%sg_n = sg_n
! !                 this%sg_sig = sg_sig
! !             end if
!             !call setSupercellTilt( this )
!             return
!         end subroutine setDeviationParameter
    
        
        pure real(kind=real64) function getDeviationParameter( this,g )
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
            type(DynamicalTwoBeamImaging),intent(in)                ::      this
            real(Kind=real64),dimension(3),intent(in),optional      ::      g
            real(kind=real64)                       ::      modk,vv,eps
            
            if (present(g)) then
                modk = 2*PI/wavelength( this%E*1000 )
                vv = HBAR*modk/ME
                eps = energy_deviation_parameter( vv,(/0.0d0,0.0d0,1.0d0/),g )      
                getDeviationParameter = energyToLengthDeviationParameter( vv,eps )
            else
                getDeviationParameter = this%sg
            end if
            return
        end function getDeviationParameter
    
        subroutine setIntegrationPoints( this,mxy,mzz )
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ 
            type(DynamicalTwoBeamImaging),intent(inout)     ::      this
            integer,intent(in)                              ::      mxy,mzz
            this%mxy = mxy
            this%mzz = mzz
            return
        end subroutine setIntegrationPoints
    
        subroutine setMicroscopeParameters( this,k,g,E )
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ 
            type(DynamicalTwoBeamImaging),intent(inout)     ::      this
            real(kind=real64),dimension(3),intent(in)       ::      k,g
            real(kind=real64),intent(in)                    ::      E
            !real(kind=real64)                       ::      modk
            
            this%k = k
            this%g = g
            this%E = E
            
        !---    find the electron beam vector - will always be along the z- direction in lab space, regardless of the orientation of the crystal and the tilt
!            modk = 2*PI/wavelength( this%E*1000 )       !   wave vector of electron (A^-1)            
!            this%kk(1:3) = modk * (/ 0,0,1 /)           !   beam vector direction is always along [001]
            
            
            
            return
        end subroutine setMicroscopeParameters
    
        subroutine setExtinctionDistances( this,xi0,xig )
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ 
            type(DynamicalTwoBeamImaging),intent(inout)     ::      this
            real(kind=real64),intent(in)                    ::      xi0,xig
            this%xi0 = xi0
            this%xig = xig
            return
        end subroutine setExtinctionDistances
     
        
        
        subroutine setVoxelSupercellSize2( this,cuboidSides )
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ 
    !*      set the voxel supercell dimensions
    !*      The voxel supercell (this%gksuper) is aligned with the g- and k- vectors.
    !*      Its size is either determined by the original input supercell if cuboidSides=0
    !*      or set to be equal to cuboidSides if specified
            type(DynamicalTwoBeamImaging),intent(inout)     ::      this
            real(kind=real64),dimension(3),intent(in)       ::      cuboidSides
            
            real(kind=real64),dimension(3)              ::      zz,gg
            real(kind=real64),dimension(3,3)            ::      a_super , b_cell
            real(kind=real64)                           ::      xlen,ylen,zlen
            type(SimpleSupercell)                       ::      gksuper
            integer,dimension(3)                        ::      hkl
            
            
            integer                             ::      nProcs,rank,ierror
            
            call MPI_COMM_SIZE(MPI_COMM_WORLD, nProcs, ierror)
            call MPI_COMM_RANK(MPI_COMM_WORLD, rank, ierror)
            
            
            
        !--- this%k is given as lattice vector [uvw] - convert to real space vector
            hkl(1:3) = nint(this%g(1:3))
            call rotateToStandardFrame( this%k,this%a_cell,hkl , this%U  )
                     
            

            
        !---    find matrix of reciprocal lattice vectors   
            a_super = matmul( this%U,getSuperA(this%super) )                    !   crystal axes in rotated frame
            b_cell = getReciprocalLatticeVectors( this%a_cell ) 
            
        !---    find vector direction of desired g-vector (eg [200]) in crystal frame   g = B [ hkl ]
            gg(1:3) = reflection( b_cell,this%g )
            gg(1:3) = this%U(1:3,1)*gg(1) + this%U(1:3,2)*gg(2) + this%U(1:3,3)*gg(3)                        
            gg = gg/norm2(gg)
            zz(1:3) = (/ 0.0d0,0.0d0,1.0d0 /)                                   !   zone axis in rotated frame is always [001]
            
            
        !---    generate a new supercell, whose axes are oriented with the g- and k- vectors.            
            call suggestSupercellOrientedWithN( a_super,zz,gg,this%a0,gksuper )   
                            
            
        !---    find the new dimensions of the voxel supercell                               
            xlen = getSupercellSideLength( gksuper,1 )   
            ylen = getSupercellSideLength( gksuper,2 )   
            zlen = getSupercellSideLength( gksuper,3 )   
            
            if (cuboidSides(3) == 0) then
                !   set the voxel depth from the supercell            
                this%voxelSupercellDimensions(3) = zlen                        
                
            else            
                !   set the voxel depth by hand
                this%voxelSupercellDimensions(3) = cuboidSides(3)                
            end if
             
            if (cuboidSides(1) == 0) then
                if (cuboidSides(2) == 0) then
                    !   set x- and y- extents of TEM image from the supercell
                    this%voxelSupercellDimensions(1) = xlen * sqrt( zlen/this%voxelSupercellDimensions(3) )
                    this%voxelSupercellDimensions(2) = ylen * sqrt( zlen/this%voxelSupercellDimensions(3) )
                else
                    this%voxelSupercellDimensions(1) = xlen * (ylen*zlen)/( cuboidSides(2)*this%voxelSupercellDimensions(3) )
                    this%voxelSupercellDimensions(2) = cuboidSides(2)   
                end if
            else    
                if (cuboidSides(2) == 0) then
                    !   set y- extents of TEM image from the supercell
                    this%voxelSupercellDimensions(1) = cuboidSides(1)   
                    this%voxelSupercellDimensions(2) = ylen * (xlen*zlen)/( cuboidSides(1)*this%voxelSupercellDimensions(3) )
                else
                    this%voxelSupercellDimensions(1) = cuboidSides(1)
                    this%voxelSupercellDimensions(2) = cuboidSides(2)
                end if
            end if
            
            if (rank==0) then
                print *,"Lib_DynamicalTwoBeamImaging::setVoxelSupercellSize2 info -"
                write(*,fmt='(5(a,f12.3))') "    input hints                              ",cuboidSides(1),",",cuboidSides(2),",",cuboidSides(3)," A"
                write(*,fmt='(5(a,f12.3))') "    side lengths of oriented supercell       ",xlen,",",ylen,",",zlen," A"
                write(*,fmt='(5(a,f12.3))') "    TEM image size                           ",this%voxelSupercellDimensions(1),",",this%voxelSupercellDimensions(2)," A"
                write(*,fmt='(5(a,f12.3))') "    supercell thickness                      ",this%voxelSupercellDimensions(3)," A"
            end if            
                     
            return
        end subroutine setVoxelSupercellSize2
    
        subroutine setConventionalUnitCell( this,a_cell )
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
            type(DynamicalTwoBeamImaging),intent(inout)     ::      this
            real(kind=real64),dimension(3,3),intent(in)     ::      a_cell
            
            this%a_cell(1:3,1:3) = a_cell(1:3,1:3)
            
            return
        end subroutine setConventionalUnitCell
    
    
        subroutine findImageParameters( this,aperpx,xmax,modg,Mx,My,Mz )
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      given the input settings, compute the scale bar (Angstrom per pixel)
    !*      the dimension of the image in A 
    !*      the magnitude of the g-vector 
    !*      and the number of pixels
            type(DynamicalTwoBeamImaging),intent(in)    ::      this
            real(kind=real64),intent(out)               ::      aperpx
            real(kind=real64),dimension(3),intent(out)  ::      xmax
            real(kind=real64),intent(out)               ::      modg
            integer,intent(out)                         ::      Mx,My,Mz
            
            real(kind=real64),dimension(3)          ::      gg
            real(kind=real64),dimension(3,3)        ::      b_cell
            
            
        !---    find vector direction of desired g-vector aperture (eg [200]) in crystal frame   g = B [ hkl ]
            b_cell = getReciprocalLatticeVectors( this%a_cell ) 
            gg(1:3) = reflection( b_cell,this%g )
            modg = norm2(gg)

            xmax(1) = getSuperCellSideLength( this%gksuper,1 )
            xmax(2) = getSuperCellSideLength( this%gksuper,2 )
            xmax(3) = getSuperCellSideLength( this%gksuper,3 )
            
            
            Mx = getNx( this%gksuper ) 
            My = getNy( this%gksuper ) 
            Mz = getNz( this%gksuper ) 
            aperpx = sqrt( xmax(1)*xmax(2) / (Mx*My) )
            
            return
        end subroutine findImageParameters
             
            
            
        subroutine findFoilThickness( this,x,RU,d,top,thick )
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      compute atom extents under assumption of rotation RU
    !*  
            type(DynamicalTwoBeamImaging),intent(inout)     ::      this
            real(kind=real64),dimension(:,:),intent(in)     ::      x           !   (3,nAtoms)  positions of atoms
            real(kind=real64),dimension(3,3),intent(in)     ::      RU
            real(kind=real64),dimension(3),intent(out)      ::      d
            real(kind=real64),dimension(:,:),allocatable,intent(out)    ::      top,thick
            
            
             
            real(kind=real64)                       ::      id3
             
            real(kind=real64),dimension(3,3)        ::      D_super     
            real(kind=real64),dimension(3)          ::      dd,xx,Ruvw
           
            integer                                 ::      Nx,Ny,ii , nAtoms,nn , ix,iy,iz
              
            integer                             ::      nProcs,rank,ierror
            
            call MPI_COMM_SIZE(MPI_COMM_WORLD, nProcs, ierror)
            call MPI_COMM_RANK(MPI_COMM_WORLD, rank, ierror)
             
        !---    check the true thickness of the foil - if there is a surface this may not be the box dimension.
             
            nAtoms = size(x,dim=2)
             
        !---    find the side lengths of the supercell in the rotated frame
            dd(1:3) = this%voxelSupercellDimensions
            call findTransformedSupercell( getSuperA(this%super),RU,D_super,this%delta,this%uvw,this%nuvw , dd )
            
        !---    allocate thickness array in x-y
            dd(1) = D_super(1,1)
            dd(2) = D_super(2,2)
            dd(3) = D_super(3,3)
            d(1:3) = dd(1:3)
            
            if (rank==0) write(*,fmt='(5(a,f12.3))') "Lib_DynamicalTwoBeamImaging::findFoilThickness info - standard frame cell sides ",dd(1),"x",dd(2),"x",dd(3)," (A)"
            Nx = max(1,nint(dd(1)/(3*this%a0)))
            Ny = max(1,nint(dd(2)/(3*this%a0)))            
            dd(1) = Nx/dd(1)
            dd(2) = Ny/dd(2)
            id3 = 1/dd(3)
            allocate(top(0:Nx-1,0:Ny-1))        ;   top = huge(1.0)
            allocate(thick(0:Nx-1,0:Ny-1))      ;   thick =-huge(1.0)
            !if (rank==0) print *,"Lib_DynamicalTwoBeamImaging::findFoilThickness info - finding thickness on mesh ",Nx,Ny 
            
            do nn = 1,this%nuvw     !   periodic replicas
            
                Ruvw(1:3) = this%uvw(1,nn)*getSuperVector(this%super,1)    &
                          + this%uvw(2,nn)*getSuperVector(this%super,2)    &
                          + this%uvw(3,nn)*getSuperVector(this%super,3)
                      
                do ii = 1,nAtoms
                    xx(1:3) = x(1:3,ii) + Ruvw(1:3)         !   position of the atom including periodic repeat.
                    xx(1:3) = RU(1:3,1)*xx(1) + RU(1:3,2)*xx(2) + RU(1:3,3)*xx(3) 
                    xx(1:3) = xx(1:3) + this%delta(1:3)
                    ix = floor( dd(1)*xx(1) ) ; if ( ix*(Nx-1-ix)<0 ) cycle
                    iy = floor( dd(2)*xx(2) ) ; if ( iy*(Ny-1-iy)<0 ) cycle
                    iz = floor( id3*xx(3) ) ; if (iz /= 0) cycle
                    top(ix,iy) = min( top(ix,iy),xx(3) )
                    thick(ix,iy) = max( thick(ix,iy),xx(3) )
                end do           
                
            end do
            
            !if (rank==0) print *,"Lib_DynamicalTwoBeamImaging::findFoilThickness info - top,bottom atoms after zone axis rotation found at ",minval(top),maxval(thick)
            
            where (thick > top)
                thick = thick - top
            elsewhere
                thick = 0
            end where
             
            
            if (rank==0) write(*,fmt='(5(a,f12.3))') "Lib_DynamicalTwoBeamImaging::findFoilThickness info - min/max thickness found ",minval(thick),",",maxval(thick)," (A)"
            return
        end subroutine findFoilThickness
        
        
        subroutine setSupercellTilt0( this )
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      set the tilt angle of the foil according to desired diffraction conditions.
    !*      - this function sets zero tilt angle if this%sg = SG_UNSET 
    !*        or it rotates about an axis orthogonal to g- & k- until sg is the desired value.
    
    
            type(DynamicalTwoBeamImaging),intent(inout)     ::      this
            
            real(kind=real64)                       ::      vv,modg,bb,aa,cc,eps! , theta1,theta2 !,modk
             
            real(kind=real64),dimension(3,3)        ::      a_super,a_cell,a_cell0 ,RU,d_super,b_cell
            real(kind=real64),dimension(3)          ::      gg  !,kk
            integer                                 ::      Nx,Ny,Nz , ii
            integer,dimension(3)                    ::      hkl  

            real(kind=real64)                       ::      modk  , theta1,theta2,theta3,theta4,sg1,sg2,sg3,sg4
            
            integer                             ::      nProcs,rank,ierror
           ! character(len=16)                   ::      aaaa
            call MPI_COMM_SIZE(MPI_COMM_WORLD, nProcs, ierror)
            call MPI_COMM_RANK(MPI_COMM_WORLD, rank, ierror)
                                  
            
        !---    find rotation to bring box to standard frame ( k -> [001], g -> [100] )            
            hkl(1:3) = nint(this%g(1:3))               
            call rotateToStandardFrame( this%k,this%a_cell,hkl , this%U  )
            
        !---    k vector after rotation is always [001] in lab frame.
            !modk = 2*PI/wavelength( this%E*1000 )       !   wave vector of electron (A^-1)            
            !this%kk(1:3) = modk * (/ 0,0,1 /)           !   beam vector direction is always along [001]
        
            
        !---    g vector will vary after rotation, but its magnitude doesn't. Start with the standard frame g-vector.
        !       note: g = 2 pi B ~g, with B^-1 = ( U A )^T
            a_cell0 = matmul( this%U,this%a_cell )
            b_cell = getReciprocalLatticeVectors( a_cell0 )      
            gg(1:3) = reflection( b_cell,this%g )
            modg = norm2(gg)            
            
                                                        
                                                         
                                                         
            
          !  print *,"gg before rotation ",gg
            
        !---    this will now give us energy deviation parameter _before_ rotation
            !vv = velocity( this%E*1000 )                 !   velocity of electron A/fs
            modk = 2*PI/wavelength( this%E*1000 )
            vv = HBAR*modk/ME
            eps = energy_deviation_parameter( vv,(/0.0d0,0.0d0,1.0d0/),gg )            
            if (rank==0) write(*,fmt='(a,f12.5,a)') "Lib_DynamicalTwoBeamImaging::setSupercellTilt0 info - deviation parameter before rotation ",energyToLengthDeviationParameter( vv,eps )," (1/A)"
            
        !---    sg = - g_3/(2 pi) - hbar |g^2| /(4 pi m v )       
        
        !       if UA -> RUA then B -> R B
        !       so 2 pi sg -> sin(theta) gg(1) - cos(theta) gg(3) - HBAR*modg*modg/( 2 ME * vv )
            if (this%sg == SG_UNSET) then
                if (rank==0) write(*,fmt='(a)') "Lib_DynamicalTwoBeamImaging::setSupercellTilt0 info - desired deviation parameter unset"
                this%sg = energyToLengthDeviationParameter( vv,eps )
                this%R = reshape( (/1,0,0,0,1,0,0,0,1/),(/3,3/) )
            else
            
                !   solve for desired sg
                !   - sin(theta) gg(1) + cos(theta) gg(3) = aa
                !   with aa = - 2 pi |k| sg - |g|^2
!                if (rank==0) write(*,fmt='(a,f12.5,a)') "Lib_DynamicalTwoBeamImaging::setSupercellTilt0 info - desired deviation parameter     ",this%sg," (1/A)"
!                !modk = 2*PI/wavelength( this%E*1000 ) 
!                aa = - 2*PI*this%sg - modg*modg/(2*modk)
!                
!
!                cc = gg(1)*gg(1) + gg(3)*gg(3)
!                bb = cc - aa*aa
!                bb = sqrt( max(0.0d0,bb) )
!                    
!                costheta  = (gg(3)*aa + gg(1)*bb)/cc
!                costheta2 = (gg(3)*aa - gg(1)*bb)/cc
!                
!                !if (rank==0) print *,"v,k,a,b,c",vv,modk,aa,bb,cc,"g",gg,"costheta",costheta,costheta2
!                
!                if (costheta2>costheta) costheta = costheta2
!                sintheta = sqrt( max(0.0d0, 1 - costheta*costheta ) )
!                
!                this%R = reshape( (/costheta,0.0d0,-sintheta , 0.0d0,1.0d0,0.0d0 , +sintheta,0.0d0,costheta/),(/3,3/) )
                
            
!                !   solve for desired sg
!                !   sin(theta) gg(1) - cos(theta) gg(3) = aa
                aa = 2 * PI * this%sg + HBAR*modg*modg/( 2 * ME * vv )
                cc = gg(1)*gg(1) + gg(3)*gg(3)
                bb = cc - aa*aa
                if (rank==0) write(*,fmt='(a,f12.3,a)') "Lib_DynamicalTwoBeamImaging::setSupercellTilt0 info - desired deviation parameter ",this%sg
                if (bb<0) then
                    if (rank==0) write(*,fmt='(a)') "Lib_DynamicalTwoBeamImaging::setSupercellTilt0() warning - s_g/|g| too large."
                    this%sg = energyToLengthDeviationParameter( vv,eps )
                    this%R = reshape( (/1,0,0,0,1,0,0,0,1/),(/3,3/) )
                else
                    bb = sqrt(bb)
                    cc = 1/cc
                    
                !---    four possible solutions for angle theta.                    
                    theta1 = acos( - (gg(3)*aa + gg(1)*bb )*cc )
                    theta2 = acos( - (gg(3)*aa - gg(1)*bb )*cc )
 
!                     if (abs(theta2)<abs(theta1)) theta1 = theta2
!                     !theta1 = min(abs(theta1),abs(theta2))
!                     !if (this%sg<0) theta1 = -theta1
!                     
                    if (rank==0) print *,"theta1,theta2,theta ",acos( - (gg(3)*aa + gg(1)*bb )*cc ),acos( - (gg(3)*aa - gg(1)*bb )*cc ),theta1
                    
                    theta3 = -theta1
                    theta4 = -theta2
                    
                !---    find deviation parameter for each option                    
                    this%R = reshape( (/cos(theta1),0.0d0,-sin(theta1)    ,0.0d0,1.0d0,0.0d0,    sin(theta1),0.0d0,cos(theta1)/),(/3,3/) )
                    a_cell = matmul( this%R,a_cell0 )
                    b_cell = getReciprocalLatticeVectors( a_cell )      
                    gg(1:3) = reflection( b_cell,this%g )                    
                    eps = energy_deviation_parameter( vv,(/0.0d0,0.0d0,1.0d0/),gg )
                    sg1 = energyToLengthDeviationParameter( vv,eps )
                    

                    this%R = reshape( (/cos(theta2),0.0d0,-sin(theta2)    ,0.0d0,1.0d0,0.0d0,    sin(theta2),0.0d0,cos(theta2)/),(/3,3/) )
                    a_cell = matmul( this%R,a_cell0 )
                    b_cell = getReciprocalLatticeVectors( a_cell )      
                    gg(1:3) = reflection( b_cell,this%g )                    
                    eps = energy_deviation_parameter( vv,(/0.0d0,0.0d0,1.0d0/),gg )
                    sg2 = energyToLengthDeviationParameter( vv,eps )
                    

                    this%R = reshape( (/cos(theta3),0.0d0,-sin(theta3)    ,0.0d0,1.0d0,0.0d0,    sin(theta3),0.0d0,cos(theta3)/),(/3,3/) )
                    a_cell = matmul( this%R,a_cell0 )
                    b_cell = getReciprocalLatticeVectors( a_cell )      
                    gg(1:3) = reflection( b_cell,this%g )                    
                    eps = energy_deviation_parameter( vv,(/0.0d0,0.0d0,1.0d0/),gg )
                    sg3 = energyToLengthDeviationParameter( vv,eps )
                    

                    this%R = reshape( (/cos(theta4),0.0d0,-sin(theta4)    ,0.0d0,1.0d0,0.0d0,    sin(theta4),0.0d0,cos(theta4)/),(/3,3/) )
                    a_cell = matmul( this%R,a_cell0 )
                    b_cell = getReciprocalLatticeVectors( a_cell )      
                    gg(1:3) = reflection( b_cell,this%g )                    
                    eps = energy_deviation_parameter( vv,(/0.0d0,0.0d0,1.0d0/),gg )
                    sg4 = energyToLengthDeviationParameter( vv,eps )
                                                                                
                !---    which is the best?
                    if (rank==0) write(*,fmt='(3(a,f12.5))') "Lib_DynamicalTwoBeamImaging::setSupercellTilt0() deviation parameter using theta = ",theta1," sg = ",sg1," (1/A)"
                    if (rank==0) write(*,fmt='(3(a,f12.5))') "Lib_DynamicalTwoBeamImaging::setSupercellTilt0() deviation parameter using theta = ",theta2," sg = ",sg2," (1/A)"
                    if (rank==0) write(*,fmt='(3(a,f12.5))') "Lib_DynamicalTwoBeamImaging::setSupercellTilt0() deviation parameter using theta = ",theta3," sg = ",sg3," (1/A)"
                    if (rank==0) write(*,fmt='(3(a,f12.5))') "Lib_DynamicalTwoBeamImaging::setSupercellTilt0() deviation parameter using theta = ",theta4," sg = ",sg4," (1/A)"
                                            
                !---    remove possibility of choosing a 180 degree rotation by fudging the answer so its never selected                   
                    if (abs(theta1)*2 > PI) sg1 = huge(1.0)
                    if (abs(theta2)*2 > PI) sg2 = huge(1.0)
                    if (abs(theta3)*2 > PI) sg3 = huge(1.0)
                    if (abs(theta4)*2 > PI) sg4 = huge(1.0)                    
                    
                !---    select best option                    
                    ii = minloc( abs( (/ this%sg-sg1,this%sg-sg2,this%sg-sg3,this%sg-sg4 /) ) , dim=1 )
                    select case(ii)
                        case(1)
                            theta1 = theta1
                        case(2)
                            theta1 = theta2                            
                        case(3)
                            theta1 = theta3
                        case(4)
                            theta1 = theta4
                    end select 
                    
                     
                    this%R = reshape( (/cos(theta1),0.0d0,-sin(theta1)    ,0.0d0,1.0d0,0.0d0,    sin(theta1),0.0d0,cos(theta1)/),(/3,3/) )
                   
                end if        
                


            !---    CHECK
                a_cell = matmul( this%R,a_cell0 )
                b_cell = getReciprocalLatticeVectors( a_cell )      
                gg(1:3) = reflection( b_cell,this%g )
                
                eps = energy_deviation_parameter( vv,(/0.0d0,0.0d0,1.0d0/),gg )
                this%sg = energyToLengthDeviationParameter( vv,eps )
                
                 
                
                if (rank==0) write(*,fmt='(3(a,f12.5))') "Lib_DynamicalTwoBeamImaging::setSupercellTilt0 info - rotation angle ",theta1," (rad) = ",theta1*180.0d0/PI," (deg)"
!                if (rank==0) write(*,fmt='(a,f12.3,a)') "Lib_DynamicalTwoBeamImaging::setSupercellTilt0 info - desired deviation parameter     ",this%sg," (1/A)"
                if (rank==0) write(*,fmt='(a,f12.5,a)') "Lib_DynamicalTwoBeamImaging::setSupercellTilt0 info - deviation parameter set         ",this%sg," (1/A)"       
             
            end if 
            if (rank==0) then
                print *,"Lib_DynamicalTwoBeamImaging::setSupercellTilt0 info - rotation matrix to tilt to desired sg"    
                write (*,fmt='(3f12.5,3(a,i4))') this%R(1,:)
                write (*,fmt='(3f12.5,3(a,i4))') this%R(2,:)
                write (*,fmt='(3f12.5,3(a,i4))') this%R(3,:)
            end if                  
            
            
            RU = matmul( this%R,this%U )
            a_super = getsupera(this%super)
            
            
            call findtransformedsupercell( a_super,RU,d_super,this%delta,this%uvw,this%nuvw , this%voxelSupercellDimensions )
            
        !---    note that d is a cuboidal crystallite, so the side lengths are given by the diagonal entries.
            nx = max(this%mxy,nint(this%mxy*d_super(1,1)/this%a0))
            ny = max(this%mxy,nint(this%mxy*d_super(2,2)/this%a0))
            nz = max(this%mzz,nint(this%mzz*d_super(3,3)/this%a0))
            
            d_super(1,1) = d_super(1,1) / nx
            d_super(2,2) = d_super(2,2) / ny
            d_super(3,3) = d_super(3,3) / nz            
            
            this%gksuper = simplesupercell_ctor( d_super,nx,ny,nz )
            
             
            gg = getAperture( this )
            this%sg = getDeviationParameter( this,gg )
            if (rank==0) print *,"Lib_DynamicalTwoBeamImaging::setSupercellTilt0 info - deviation parameter ",this%sg," (1/A)"
            
             
                        
            return
        end subroutine setSupercellTilt0
        
            
        subroutine setSupercellTilt1( this,ng, extinction_distances_filename,x )
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      use the extinction distances to compute the best tilt for the desired (g,ng) two-beam condition 
    !*      if x is present, then the foil thickness is determined from the positions of the atoms rather than the raw supercell extents.
    
            type(DynamicalTwoBeamImaging),intent(inout)     ::      this
            integer,intent(in)                              ::      ng
            
            character(len=*),intent(in)                     ::      extinction_distances_filename
            real(kind=real64),dimension(:,:),intent(in),optional     ::      x           !   (3,nAtoms)  positions of atoms
            
            real(kind=real64)                       ::      vv   , z_thickness, zbar,z2bar!, modk     , eps
             
            real(kind=real64),dimension(3,3)        ::      a_super,D_super , RU !, ia_super     ,b_cell    
            real(kind=real64),dimension(3)          ::      gg,dd
            integer,dimension(3)                    ::      hkl
            real(kind=real64),dimension(:,:),allocatable    ::      top,thick
            integer                                 ::      Nx,Ny,Nz ,nn , ix,iy
            
            
            
            integer                             ::      nProcs,rank,ierror
            
            call MPI_COMM_SIZE(MPI_COMM_WORLD, nProcs, ierror)
            call MPI_COMM_RANK(MPI_COMM_WORLD, rank, ierror)
            
            
            
            hkl(1:3) = nint(ng*this%g(1:3))              !  (g,ng) spot ie g=[200] ng=3 -> [600]
            
            !a_super = getSuperA(this%super)
!            call inverse3Mat(a_super,ia_super)
            
            
            
            if (present(x)) then
            !---    check the true thickness of the foil - if there is a surface this may not be the box dimension.
                if (rank==0) write(*,fmt='(a,f12.3,a)') "Lib_DynamicalTwoBeamImaging::setSupercellTilt1 info - comparing box dimensions with atom extents "
                call rotateToStandardFrame( this%k,this%a_cell,hkl , this%U  )
                call findFoilThickness( this,x,this%U,dd,top,thick )
                zbar = 0.0d0 ; z2bar = 0.0d0 ; nn = 0 ; Nx = size(top,dim=1) ; Ny = size(top,dim=2)
                do iy = 0,Ny-1
                    do ix = 0,Nx-1
                        if (thick(ix,iy) == 0) cycle
                        zbar = zbar + thick(ix,iy)
                        z2bar = z2bar + thick(ix,iy)*thick(ix,iy)
                        nn = nn + 1
                    end do
                end do
                zbar = zbar / nn                !   mean thickness
                z2bar = z2bar / nn              !   mean square thickness
                z2bar = sqrt( max(0.0d0,z2bar - zbar*zbar) )    !   std dev thickness
                 
                if (rank==0) write(*,fmt='(2(a,f12.3),3(a,i6))') "Lib_DynamicalTwoBeamImaging::setSupercellTilt1 info - avg thickness found ",zbar," +/- ",z2bar," in ",nn,"/",Nx*Ny," cells"
                
                if (zbar + 2*this%a0 > dd(3)) then
                    if (rank==0) write(*,fmt='(a)') "Lib_DynamicalTwoBeamImaging::setSupercellTilt1 info - atoms span depth of supercell"
                    zbar = dd(3)                    
                else 
                    if (rank==0) then
                        if (z2bar > 10*this%a0) then
                            write(*,fmt='(a,f12.3,a)') "Lib_DynamicalTwoBeamImaging::setSupercellTilt1 WARNING - very rough surface found = ",z2bar/this%a0," u.c."
                        else if (z2bar > 2*this%a0) then
                            write(*,fmt='(a,f12.3,a)') "Lib_DynamicalTwoBeamImaging::setSupercellTilt1 warning - rough surface found = ",z2bar/this%a0," u.c."
                        else
                            write(*,fmt='(a,f12.3,a)') "Lib_DynamicalTwoBeamImaging::setSupercellTilt1 info - smooth surface found = ",z2bar/this%a0," u.c."
                        end if
                    end if
                end if
                deallocate(top)
                deallocate(thick)
                
            else            
                !   have not been provided with atom positions, so will take foil_thickness(3) as a bona-fide thickness of the foil ( before tilt )
                if (rank==0) write(*,fmt='(a)') "Lib_DynamicalTwoBeamImaging::setSupercellTilt1 info - using input supercell dimensions for foil thickness, not checking actual atom extents"              
                zbar = this%voxelSupercellDimensions(3)
            end if
            
            this%foilThickness = zbar
            
            
            
            vv = velocity( this%E*1000 )                 !   velocity of electron A/fs
            z_thickness = abs(this%foilThickness)
            if (rank==0) write(*,fmt='(a,f12.3,a)') "Lib_DynamicalTwoBeamImaging::setSupercellTilt1 info - foil thickness before rotation to 2 beam condition ",z_thickness
            call findBestRotationMatrix( vv,this%k,this%a_cell,this%latt, hkl , this%U,this%R , extinction_distances_filename,z_thickness ) 
            z_thickness = abs( this%foilThickness / this%R(3,3) )
            if (rank==0) write(*,fmt='(a,f12.3,a)') "Lib_DynamicalTwoBeamImaging::setSupercellTilt1 info - foil thickness after rotation to 2 beam condition  ",z_thickness      
            call tiltToDarkField( velocity( this%E*1000 ) ,this%a_cell, nint(this%g(1:3)) , this%xi0,this%xig,z_thickness,this%U,this%R )   
            z_thickness = abs( this%foilThickness / this%R(3,3) ) 
            if (rank==0) write(*,fmt='(a,f12.3,a)') "Lib_DynamicalTwoBeamImaging::setSupercellTilt1 info - foil thickness after tilt to dark field            ",z_thickness      
            this%voxelSupercellDimensions(3) = abs( this%voxelSupercellDimensions(3) / this%R(3,3) )
            if (rank==0) write(*,fmt='(a,f12.3,a)') "Lib_DynamicalTwoBeamImaging::setSupercellTilt1 info - setting box depth after rotation to                ",this%voxelSupercellDimensions(3)
            
            
            RU = matmul( this%r,this%u )
            a_super = getsupera(this%super)
            
            
            call findTransformedSupercell( a_super,RU,d_super,this%delta,this%uvw,this%nuvw , this%voxelSupercellDimensions )
            
        !---    note that d is a cuboidal crystallite, so the side lengths are given by the diagonal entries.
            nx = max(this%mxy,nint(this%mxy*d_super(1,1)/this%a0))
            ny = max(this%mxy,nint(this%mxy*d_super(2,2)/this%a0))
            nz = max(this%mzz,nint(this%mzz*d_super(3,3)/this%a0))
            
            d_super(1,1) = d_super(1,1) / nx
            d_super(2,2) = d_super(2,2) / ny
            d_super(3,3) = d_super(3,3) / nz            
            
            this%gksuper = SimpleSupercell_ctor( d_super,nx,ny,nz )
            
            !if (rank==0) then
            !    print *,"Lib_DynamicalTwoBeamImaging::setSupercellTilt1 info - voxels to use"
            !    call report(this%gksuper)
            !end if
            
            
        !---    find the electron beam vector
            !modk = 2*PI/wavelength( this%E*1000 )       !   wave vector of electron (A^-1)            
            !this%kk(1:3) = modk * (/ 0,0,1 /)           !   beam vector direction is always along [001]
            

        !---    find the g-vector in the rotated crystal       
            !call inverse3Mat( transpose(this%a_cell),b_cell )                                       
            !gg(1:3) = 2*PI*( b_cell(1:3,1)*this%g(1) + b_cell(1:3,2)*this%g(2) + b_cell(1:3,3)*this%g(3) )            
!             b_cell = getReciprocalLatticeVectors( this%a_cell )      
!             gg(1:3) = reflection( b_cell,this%g )
!             gg(1:3) = RU(1:3,1)*gg(1) + RU(1:3,2)*gg(2) + RU(1:3,3)*gg(3) 
            
            
        !---    find the deviation parameter
!             eps = energy_deviation_parameter( vv,(/0.0d0,0.0d0,1.0d0/),gg )
!             this%sg = energyToLengthDeviationParameter( vv,eps )
!             if (rank==0) write(*,fmt='(a,f12.3,a)') "Lib_DynamicalTwoBeamImaging::setSupercellTilt1 info - deviation parameter ",this%sg," (A)"
!             


            gg = getAperture( this )
            this%sg = getDeviationParameter( this,gg )
            if (rank==0) write(*,fmt='(a,f12.6,a)') "Lib_DynamicalTwoBeamImaging::setSupercellTilt1 info - deviation parameter ",this%sg," (1/A)"
            
            
            return
        end subroutine setSupercellTilt1
        
        
            
        subroutine setSupercellTilt1a( this,ng, extinction_distances_filename,x )
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      use the extinction distances to compute the best tilt for the desired (g,ng) two-beam condition 
    !*      if x is present, then the foil thickness is determined from the positions of the atoms rather than the raw supercell extents.
    
            type(DynamicalTwoBeamImaging),intent(inout)     ::      this
            real(kind=real64),intent(in)                    ::      ng
            
            character(len=*),intent(in)                     ::      extinction_distances_filename
            real(kind=real64),dimension(:,:),intent(in),optional     ::      x           !   (3,nAtoms)  positions of atoms
             
            real(kind=real64),dimension(3,3)        ::      a_super,D_super , RU !, ia_super     ,b_cell    
            real(kind=real64),dimension(3)          ::      gg, dd
 
            integer                                 ::      Nx,Ny,Nz
            
            type(Quaternion)                        ::      q1,q2
            
            integer                             ::      nProcs,rank,ierror
            
            call MPI_COMM_SIZE(MPI_COMM_WORLD, nProcs, ierror)
            call MPI_COMM_RANK(MPI_COMM_WORLD, rank, ierror)
            
            dd = this%voxelSupercellDimensions
            
            if (present(x)) then
                if (rank==0) print *,"Lib_DynamicalTwoBeamImaging::setSupercellTilt1a info - finding tilt for ng = ",floor(ng)," with atom positions"
                this%voxelSupercellDimensions = dd
                call setSupercellTilt1( this,floor(ng), extinction_distances_filename,x )
                q1 = Quaternion_ctor(this%R)
                if (rank==0) print *,""
                if (ng-floor(ng)>1.0d-4) then
                    if (rank==0) print *,""
                    if (rank==0) print *,"Lib_DynamicalTwoBeamImaging::setSupercellTilt1a info - finding tilt for ng = ",ceiling(ng)," with atom positions"
                    this%voxelSupercellDimensions = dd
                    call setSupercellTilt1( this,ceiling(ng), extinction_distances_filename,x )
                    q2 = Quaternion_ctor(this%R)
                    if (rank==0) print *,""
                else
                    q2 = q1
                end if
            else
                if (rank==0) print *,"Lib_DynamicalTwoBeamImaging::setSupercellTilt1a info - finding tilt for ng = ",floor(ng)," without atom positions"
                this%voxelSupercellDimensions = dd
                call setSupercellTilt1( this,floor(ng), extinction_distances_filename )
                q1 = Quaternion_ctor(this%R)
                if (rank==0) print *,""
                if (ng-floor(ng)>1.0d-4) then
                    if (rank==0) print *,"Lib_DynamicalTwoBeamImaging::setSupercellTilt1a info - finding tilt for ng = ",ceiling(ng)," without atom positions"
                    this%voxelSupercellDimensions = dd
                    call setSupercellTilt1( this,ceiling(ng), extinction_distances_filename )
                    q2 = Quaternion_ctor(this%R)
                    if (rank==0) print *,""
                else
                    q2 = q1
                end if
                
            end if    
            
            ! if (rank==0) then
            !     print *,""
            !     print *,"quaternion ",floor(ng)
            !     call report(q1,asRotMat=.true.)
            !     print *,"quaternion ",ceiling(ng)
            !     call report(q2,asRotMat=.true.)
            !     write(*,fmt='(a,f8.3)')"quaternion ng = ",ng
            !     call report(slerp(q1,q2,ng - floor(ng)),asRotMat=.true.)
            ! end if
            
            q1 = slerp(q1,q2,ng - floor(ng))        !   interpolate rotation between points for floor(ng) and ceiling(ng)
             
            this%R = quaternionToRotMat( q1 )
                
                
            this%voxelSupercellDimensions(3) = abs( dd(3) / this%R(3,3) )
            if (rank==0) write(*,fmt='(a,f12.3,a)') "Lib_DynamicalTwoBeamImaging::setSupercellTilt1a info - setting box depth after rotation to               ",this%voxelSupercellDimensions(3)
            
            
            RU = matmul( this%R,this%U )
            a_super = getsupera(this%super)
            
            
            call findTransformedSupercell( a_super,RU,d_super,this%delta,this%uvw,this%nuvw , this%voxelSupercellDimensions )
            
        !---    note that d is a cuboidal crystallite, so the side lengths are given by the diagonal entries.
            nx = max(this%mxy,nint(this%mxy*d_super(1,1)/this%a0))
            ny = max(this%mxy,nint(this%mxy*d_super(2,2)/this%a0))
            nz = max(this%mzz,nint(this%mzz*d_super(3,3)/this%a0))
            
            d_super(1,1) = d_super(1,1) / nx
            d_super(2,2) = d_super(2,2) / ny
            d_super(3,3) = d_super(3,3) / nz            
            
            this%gksuper = SimpleSupercell_ctor( d_super,nx,ny,nz )
             

            gg = getAperture( this )
            this%sg = getDeviationParameter( this,gg )
            if (rank==0) write(*,fmt='(a,f12.6,a)') "Lib_DynamicalTwoBeamImaging::setSupercellTilt1a info - deviation parameter ",this%sg," (1/A)"
            if (rank==0) print *,""
            
            return
        end subroutine setSupercellTilt1a
                
            
        subroutine setSupercellTilt2( this, U,R )
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ 
    !*      set a previously determined tilt 
    
            type(DynamicalTwoBeamImaging),intent(inout)     ::      this
            
            real(kind=real64),dimension(3,3),intent(in)     ::      U,R
            
            
!            real(kind=real64)                       ::      vv  , eps   !, modk
             
            real(kind=real64),dimension(3,3)        ::      a_super,D_super , RU           !,b_cell
            real(kind=real64),dimension(3)          ::      gg
            
            integer                                 ::      Nx,Ny,Nz 
            
            
            
            integer                             ::      nProcs,rank,ierror
            
            call MPI_COMM_SIZE(MPI_COMM_WORLD, nProcs, ierror)
            call MPI_COMM_RANK(MPI_COMM_WORLD, rank, ierror)
             
            this%U = U
            this%R = R
            
            a_super = getSuperA(this%super)
            RU = matmul( this%R,this%U )
            
            
!           if (this%voxelSupercellDimensions(3)>0) then
!               this%voxelSupercellDimensions(3) = this%voxelSupercellDimensions(3) / this%R(3,3)      
!               if (rank==0) print *,"Lib_DynamicalTwoBeamImaging::setSupercellTilt2 info - new z-thickness after rotation is ",this%voxelSupercellDimensions(3)
!           end if  
            
            call findTransformedSupercell( a_super,RU,D_super,this%delta,this%uvw,this%nuvw , this%voxelSupercellDimensions )
            
        !---    note that D is a cuboidal crystallite, so the side lengths are given by the diagonal entries.
            Nx = max(this%mxy,nint(this%mxy*D_super(1,1)/this%a0))
            Ny = max(this%mxy,nint(this%mxy*D_super(2,2)/this%a0))
            Nz = max(this%mzz,nint(this%mzz*D_super(3,3)/this%a0))
            
            D_super(1,1) = D_super(1,1) / Nx
            D_super(2,2) = D_super(2,2) / Ny
            D_super(3,3) = D_super(3,3) / Nz            
            
            this%gksuper = SimpleSupercell_ctor( D_super,Nx,Ny,Nz )
            
        !---    find the electron beam vector
           !modk = 2*PI/wavelength( this%E*1000 )       !   wave vector of electron (A^-1)            
           !this%kk(1:3) = modk * (/ 0,0,1 /)           !   beam vector direction is always along [001]
            

        !---    find the g-vector in the rotated crystal       
            !call inverse3Mat( transpose(this%a_cell),b_cell )                                       
            !gg(1:3) = 2*PI*( b_cell(1:3,1)*this%g(1) + b_cell(1:3,2)*this%g(2) + b_cell(1:3,3)*this%g(3) )            
!             b_cell = getReciprocalLatticeVectors( this%a_cell )      
!             gg(1:3) = reflection( b_cell,this%g )
!             gg(1:3) = RU(1:3,1)*gg(1) + RU(1:3,2)*gg(2) + RU(1:3,3)*gg(3) 
            
            
        !---    find the deviation parameter
        !   vv = velocity( this%E*1000 )                 !   velocity of electron A/fs
        !   eps = energy_deviation_parameter( vv,(/0.0d0,0.0d0,1.0d0/),gg )
        !   this%sg = energyToLengthDeviationParameter( vv,eps )
        !   if (rank==0) write(*,fmt='(a,f12.6,a)') "Lib_DynamicalTwoBeamImaging::setSupercellTilt2 info - deviation parameter ",this%sg," (1/A)"
            
             
            gg = getAperture( this )
            this%sg = getDeviationParameter( this,gg )
            if (rank==0) print *,"Lib_DynamicalTwoBeamImaging::setSupercellTilt2 info - deviation parameter ",this%sg," (1/A)"
            
             
            return
        end subroutine setSupercellTilt2
        
        
    !---
    
        subroutine setUseDensityField2( this,is )
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
            type(DynamicalTwoBeamImaging),intent(inout)     ::      this            
            logical,intent(in)                              ::      is
            this%useDensityField = is
            return
        end subroutine setUseDensityField2
        
!     

        pure function getAperture( this, BeamTilt ) result( g ) 
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
            type(DynamicalTwoBeamImaging),intent(in)                ::      this          
            real(kind=real64),dimension(3,3),intent(in),optional    ::      BeamTilt     !   additional tilt to perform convergent weak-beam dark-field imaging, or beam precession
            real(kind=real64),dimension(3)                          ::      g
            
            real(kind=real64),dimension(3,3)    ::      b_cell,RU
            
    
            b_cell = getReciprocalLatticeVectors( this%a_cell )             
            g(1:3) = reflection( b_cell,this%g )
            
            RU = matmul( this%R,this%U ) 
            if (present(BeamTilt)) RU = matmul(BeamTilt,RU)
                        
            g(1:3) = RU(1:3,1)*g(1) + RU(1:3,2)*g(2) + RU(1:3,3)*g(3)             
            
            return
        end function getAperture
        
        
            
!      
        subroutine findGdotR( this, x, lc3d,lc3d_delta,gdotr , BeamTilt )
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ 
    !*      compute g.r at each atomic position in rotated and replicated basis
    !*      return a link-cell list with atom positions and exp( i g.r )
            
            type(DynamicalTwoBeamImaging),intent(inout)             ::      this
            real(kind=real64),dimension(:,:),intent(in)             ::      x           !   (3,nAtoms)  positions of atoms
            type(LinkCell3d),intent(inout)                          ::      lc3d
            real(kind=real64),dimension(3),intent(in)               ::      lc3d_delta       !   offset of link cell list 
            real(kind=real64),dimension(:),pointer                  ::      gdotr
            
            real(kind=real64),dimension(3,3),intent(in),optional    ::      BeamTilt     !   additional tilt to perform convergent weak-beam dark-field imaging, or beam precession
            
            
            
            
        !---    spatial extent of lc3d            
            type(SimpleSupercell)               ::      super
            
            real(kind=real64),dimension(:),pointer          ::      gdotr_tmp       !   (1:nAtoms) phase factor computed at atom sites
            
            integer             ::      nAtoms,nPoints,jj,ii,nn
            real(kind=real64),dimension(3)      ::      Ruvw,yy,gg,zz
            real(kind=real64),dimension(3,3)    ::      RU      !,b_cell
            logical             ::      ok
            
            
            
            !print *,"Lib_DynamicalTwoBeamImaging::findGdotR"
            
            super = getSuper(lc3d)
            call clear(lc3d)             
            nAtoms = size(x,dim=2)
            
            !call report(lc3d)
            !print *,"lc3d_delta ",lc3d_delta
            if (.not. associated(gdotr)) then
                nPoints = int(nAtoms*1.5)
                allocate(gdotr(nPoints))
            else
                nPoints = size(gdotr)
            end if
                          
            
            RU = matmul( this%R,this%U ) 
            
        !---    finally add a beam tilt for convergent weak beam or precession             
            if (present(BeamTilt)) RU = matmul(BeamTilt,RU)
            
             
        !---    find vector direction of desired g-vector aperture (eg [200]) in crystal frame   g = B [ hkl ]    
            if (present(BeamTilt)) then                
                gg = getAperture( this,BeamTilt )
            else
                gg = getAperture( this )
            end if
            
        !--- 
              
            ii = 0                  !   number of points recorded in lc3d
            do nn = 1,this%nuvw     !   periodic replicas
                
                Ruvw(1:3) = this%uvw(1,nn)*getSuperVector(this%super,1)    &
                          + this%uvw(2,nn)*getSuperVector(this%super,2)    &
                          + this%uvw(3,nn)*getSuperVector(this%super,3)
                  
                do jj = 1,nAtoms
                !   compute y = RU ( x + A [uvw] ) + delta
                    yy(1:3) = x(1:3,jj) + Ruvw(1:3)
                    yy(1:3) = RU(1:3,1)*yy(1) + RU(1:3,2)*yy(2) + RU(1:3,3)*yy(3) 
                    yy(1:3) = yy(1:3) + this%delta(1:3)
                    zz(1:3) = yy(1:3) - lc3d_delta(1:3)
                    
                !   check if this atom is within lc3d                    
                    
                    ok = pointInSupercell( super,zz )
                    
                    

                    if (ok) then
                        !   yy is within central periodic replica of lc3d  
                        
                        ii = ii + 1
                        
                        if (ii > nPoints) then
                            allocate(gdotr_tmp(nPoints*2))
                            gdotr_tmp(1:nPoints) = gdotr(1:nPoints)
                            deallocate(gdotr)
                            gdotr => gdotr_tmp
                            nPoints = nPoints*2                                
                        end if
                        
                        call add(lc3d,ii,zz)
                        gdotr(ii) = gg(1)*yy(1) + gg(2)*yy(2) + gg(3)*yy(3)
                        
                    end if    
                end do        ! atoms                
            end do            ! replicas
            
            nPoints = ii
          !  print *,"Lib_DynamicalTwoBeamImaging::findGdotR() info - atom count ",nAtoms," -> ",nPoints," after rotation and replicas"
!             
            
            
            return
        end subroutine findGdotR
            
            
    
        subroutine findGdotR_aperture( this, x, lc3d,lc3d_delta,gdotr,apertureShift , BeamTilt )
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      compute g.r at each atomic position in rotated and replicated basis
    !*      return a link-cell list with atom positions and exp( i g.r )
    !*      this version finds gdotr for multiple g-vectors in aperture            
            type(DynamicalTwoBeamImaging),intent(inout)             ::      this
            real(kind=real64),dimension(:,:),intent(in)             ::      x           !   (3,nAtoms)  positions of atoms
            type(LinkCell3d),intent(inout)                          ::      lc3d
            real(kind=real64),dimension(3),intent(in)               ::      lc3d_delta       !   offset of link cell list 
            real(kind=real64),dimension(:,:),pointer                ::      gdotr
            
            real(kind=real64),dimension(:,:),intent(in)             ::      apertureShift     !   additional shift on g-vector to sample the aperture
            
            
            real(kind=real64),dimension(3,3),intent(in),optional    ::      BeamTilt     !   additional tilt to perform convergent weak-beam dark-field imaging, or beam precession
            
            
            
        !---    spatial extent of lc3d            
            type(SimpleSupercell)               ::      super
            
            real(kind=real64),dimension(:,:),pointer          ::      gdotr_tmp       !   (1:nAtoms) phase factor computed at atom sites
            
            integer             ::      nAtoms,nPoints,jj,ii,nn,kk
            real(kind=real64),dimension(3)      ::      Ruvw,yy,gg,zz
            real(kind=real64),dimension(3,3)    ::      RU      !,b_cell
            logical             ::      ok
            integer             ::      nAperture
            real(kind=real64),dimension(3,size(apertureShift,dim=2))        ::      gg_shift
            
            
            !print *,"Lib_DynamicalTwoBeamImaging::findGdotR"
            
            super = getSuper(lc3d)
            call clear(lc3d)             
            nAtoms = size(x,dim=2)
            nAperture = size(apertureShift,dim=2)
            
            !call report(lc3d)
            !print *,"lc3d_delta ",lc3d_delta
            if (.not. associated(gdotr)) then
                nPoints = int(nAtoms*1.5)
                allocate(gdotr(nPoints,nAperture))
            else
                nPoints = size(gdotr)
            end if
                          
            
!             
!         !---    find matrix of reciprocal lattice vectors   
!             b_cell = getReciprocalLatticeVectors( this%a_cell ) 
!             
!         !---    find vector direction of desired g-vector aperture (eg [200]) in crystal frame   g = B [ hkl ]
!             gg(1:3) = reflection( b_cell,this%g )
            
            RU = matmul( this%R,this%U ) 
            
        !---    finally add a beam tilt for convergent weak beam or precession             
            if (present(BeamTilt)) RU = matmul(BeamTilt,RU)
            
            !this%gg(1:3) = this%U(1:3,1)*gg(1) + this%U(1:3,2)*gg(2) + this%U(1:3,3)*gg(3)           
            !this%gg(1:3) = RU(1:3,1)*gg(1) + RU(1:3,2)*gg(2) + RU(1:3,3)*gg(3)           
            
            
            if (present(BeamTilt)) then                
                gg = getAperture( this,BeamTilt )
            else
                gg = getAperture( this )
            end if
            
            
            
            do kk = 1,nAperture 
                gg_shift(1:2,kk) = gg(1:2) + apertureShift(1:2,kk)
                gg_shift(3,kk) = gg(3)
            end do
             
        !--- 
              
            ii = 0                  !   number of points recorded in lc3d
            do nn = 1,this%nuvw     !   periodic replicas
                
                Ruvw(1:3) = this%uvw(1,nn)*getSuperVector(this%super,1)    &
                          + this%uvw(2,nn)*getSuperVector(this%super,2)    &
                          + this%uvw(3,nn)*getSuperVector(this%super,3)
                  
                do jj = 1,nAtoms
                !   compute y = RU ( x + A [uvw] ) + delta
                    yy(1:3) = x(1:3,jj) + Ruvw(1:3)
                    yy(1:3) = RU(1:3,1)*yy(1) + RU(1:3,2)*yy(2) + RU(1:3,3)*yy(3) 
                    yy(1:3) = yy(1:3) + this%delta(1:3)
                    zz(1:3) = yy(1:3) - lc3d_delta(1:3)
                    
                !   check if this atom is within lc3d                    
                    
                    ok = pointInSupercell( super,zz )
                    
                    

                    if (ok) then
                        !   yy is within central periodic replica of lc3d  
                        
                        ii = ii + 1
                        
                        if (ii > nPoints) then
                            allocate(gdotr_tmp(nPoints*2,nAperture))
                            gdotr_tmp(1:nPoints,1:nAperture) = gdotr(1:nPoints,1:nAperture)
                            deallocate(gdotr)
                            gdotr => gdotr_tmp
                            nPoints = nPoints*2                                
                        end if
                        
                        call add(lc3d,ii,zz)
                        do kk = 1,nAperture 
                            gdotr(ii,kk) = gg_shift(1,kk)*yy(1) + gg_shift(2,kk)*yy(2) + gg_shift(3,kk)*yy(3)
                        end do
                        
                    end if    
                end do        ! atoms                
            end do            ! replicas
            
            nPoints = ii
          !  print *,"Lib_DynamicalTwoBeamImaging::findGdotR() info - atom count ",nAtoms," -> ",nPoints," after rotation and replicas"
!             
            
            
            return
        end subroutine findGdotR_aperture
            
            
            
    
        subroutine opXyzFile( this, x, opxyzfilename)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      output .xyz format file with atoms in correctly oriented positions
            
            type(DynamicalTwoBeamImaging),intent(inout)             ::      this
            real(kind=real64),dimension(:,:),intent(in)             ::      x           !   (3,nAtoms)  positions of atoms
            character(len=*),intent(in)                             ::      opxyzfilename 
            
            integer             ::      nAtoms,nPoints,jj,nn
            integer             ::      ix,iy,iz
            real(kind=real64)   ::      id1,id2,id3,gdoty,gTz
            real(kind=real64),dimension(3)      ::      Ruvw,yy,gg,zz,Tz,gnorm
            
            real(kind=real64),dimension(3,3)    ::      d_super,RU , TT     !,b_cell
            logical             ::      ok,opgTz
            
            character(len=256)  ::      dummy
            type(XYZFile)       ::      xyz
            
            print *,"Lib_DynamicalTwoBeamImaging::opXyzFile"
            
            nAtoms = size(x,dim=2)
            nPoints = int(nAtoms*1.5)
              
            
            d_super = getSuperA(this%gksuper)
            id1 = 1/d_super(1,1)       !   inverse side length of cuboidal crystallite
            id2 = 1/d_super(2,2)
            id3 = 1/d_super(3,3)
            
        !---    find matrix of reciprocal lattice vectors   
!           b_cell = getReciprocalLatticeVectors( this%a_cell ) 
            
        !---    find vector direction of desired g-vector (eg [200]) in crystal frame   g = B [ hkl ]
!         gg(1:3) = reflection( b_cell,this%g )
!         this%gg(1:3) = this%U(1:3,1)*gg(1) + this%U(1:3,2)*gg(2) + this%U(1:3,3)*gg(3)            
!         write(*,fmt='(a,3f16.3)')  "  [hkl]      ",this%g            
!         write(*,fmt='(a,3f16.6)')  "  g-vec crys ",gg
!         write(*,fmt='(a,3f16.6)')  "  aperture   ",this%gg
            gg = getAperture( this )
            gnorm = gg/norm2(gg)
            
        !---    find vector direction of z-vector [001] in crystal frame
            RU = matmul( this%R,this%U )
            zz(1:3) = RU(3,1:3)
                         
            write(*,fmt='(a,3f16.3)')  "  beam dir   ",(/ 0.0d0,0.0d0,1.0d0 /)
            write(*,fmt='(a,3f16.6)')  "  k-vec crys ",zz            
            write(*,fmt='(a,3f16.6)')  "  beam dir   ",reflection(RU,zz)
            
            
            opgTz = (size(x,dim=1)>=13)
            
            
        !---    first count number of atoms in rotated frame + replicas 
        !       have to do this first, because top of .xyz file needs it.            
            nPoints = 0             !   number of ig.r points
            do nn = 1,this%nuvw     !   periodic replicas
                
                Ruvw(1:3) = this%uvw(1,nn)*getSuperVector(this%super,1)    &
                          + this%uvw(2,nn)*getSuperVector(this%super,2)    &
                          + this%uvw(3,nn)*getSuperVector(this%super,3)
                  
                do jj = 1,nAtoms
                !   compute y = RU ( x + A [uvw] ) + delta
                    yy(1:3) = x(1:3,jj) + Ruvw(1:3)
                    yy(1:3) = RU(1:3,1)*yy(1) + RU(1:3,2)*yy(2) + RU(1:3,3)*yy(3) 
                    yy(1:3) = yy(1:3) + this%delta(1:3)
                    
                    ok = (floor( id1*yy(1) )==0).and.(floor( id2*yy(2) )==0).and.(floor( id3*yy(3) )==0)

                    if (ok) nPoints = nPoints + 1
                end do        ! atoms                
            end do            ! replicas            
            write(*,fmt='(3(a,i10))') "Lib_DynamicalTwoBeamImaging::opXyzFile() info - atom count ",nAtoms," -> ",nPoints," after rotation and replicas"
            
            
        !---    set up .xyz file            
            if (opgTz) then
                write(dummy,fmt='(a,9f10.3,a)') "Lattice=""",d_super,""" Properties=species:S:1:pos:R:3:real:R:1:imag:R:1:gTz:R:1"
                xyz = XYZFile_ctor(nPoints+8,nAtomNames=2,nColumns=6,nHeaderLines=0)
            else
                write(dummy,fmt='(a,9f10.3,a)') "Lattice=""",d_super,""" Properties=species:S:1:pos:R:3:real:R:1:imag:R:1"
                xyz = XYZFile_ctor(nPoints+8,nAtomNames=2,nColumns=5,nHeaderLines=0)
            end if            
            call setAtomNames( xyz,(/ "Eigr","C   " /) )
            call setAtomTypes( xyz,1 )
            call setFilename(xyz,opxyzfilename)
            call setColumn_Description(xyz,trim(dummy))                
            
        !---    now output atom positions             
            nPoints = 0             !   number of ig.r points
            do nn = 1,this%nuvw     !   periodic replicas
                
                Ruvw(1:3) = this%uvw(1,nn)*getSuperVector(this%super,1)    &
                          + this%uvw(2,nn)*getSuperVector(this%super,2)    &
                          + this%uvw(3,nn)*getSuperVector(this%super,3)
                  
                do jj = 1,nAtoms
                !   compute y = RU ( x + A [uvw] ) + delta
                    yy(1:3) = x(1:3,jj) + Ruvw(1:3)
                    yy(1:3) = RU(1:3,1)*yy(1) + RU(1:3,2)*yy(2) + RU(1:3,3)*yy(3) 
                    yy(1:3) = yy(1:3) + this%delta(1:3)
                    
                    ok = (floor( id1*yy(1) )==0).and.(floor( id2*yy(2) )==0).and.(floor( id3*yy(3) )==0)

                    if (ok) then
                        nPoints = nPoints + 1
                        
                        gdoty = gg(1)*yy(1) + gg(2)*yy(2) + gg(3)*yy(3)
                        
                        if (opgTz) then
                            TT(1:3,1:3) = reshape( x(5:13,jj),(/3,3/) )
                            Tz(1:3) = TT(1:3,1)*zz(1) + TT(1:3,2)*zz(1) + TT(1:3,3)*zz(3) 
                            gTz     = gnorm(1)*Tz(1) + gnorm(2)*Tz(2) + gnorm(3)*Tz(3) 
                            call setColumns( xyz,nPoints, (/ yy(1),yy(2),yy(3),cos(gdoty),sin(gdoty),gTz /) )
                        else
                            call setColumns( xyz,nPoints, (/ yy(1),yy(2),yy(3),cos(gdoty),sin(gdoty) /) )                            
                        end if                      
                    end if
                    
                end do        ! atoms                
            end do            ! replicas
            
        !---    ... and output corners of the box            
            
            do iz = 0,1
                do iy = 0,1
                    do ix = 0,1                             
                        yy(1:3) = ix*getSuperVector(this%super,1) + iy*getSuperVector(this%super,2) + iz*getSuperVector(this%super,3) 
                        yy(1:3) = RU(1:3,1)*yy(1) + RU(1:3,2)*yy(2) + RU(1:3,3)*yy(3) + this%delta(1:3)
                        nPoints = nPoints + 1 
                        if (opgTz) then
                            call setColumns( xyz,nPoints, (/yy(1),yy(2),yy(3),1.0d0,0.0d0,0.0d0 /) )
                        else
                            call setColumns( xyz,nPoints, (/yy(1),yy(2),yy(3),1.0d0,0.0d0 /) )
                        end if
                        call setAtomType( xyz,nPoints, 2 )
                    end do
                end do
            end do                                  
                        
        !---    write file
            call report(xyz)
            call output(xyz)
            call delete(xyz)            
            
            
            
            
            return
        end subroutine opXyzFile
            
            
            
            
            
            
            
            
    !---
    
        subroutine computeD2BI_img( this,x,d2big2 , precessionAngle,nPrecessionAngle, apertureAngle,nAperturePoints )
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
            type(DynamicalTwoBeamImaging),intent(inout)     ::      this
            real(kind=real64),dimension(:,:),intent(in)     ::      x           !   (3,nAtoms)  positions of atoms
            real(kind=real64),dimension(0:,0:),intent(out)  ::      d2big2
           
            real(kind=real64),intent(in),optional           ::      precessionAngle         !   beam tilt angle in radians    
            integer,intent(in),optional                     ::      nPrecessionAngle        !   number of beam tilts to average over
            
            real(kind=real64),intent(in),optional           ::      apertureAngle           !   aperture angle in radians    
            integer,intent(in),optional                     ::      nAperturePoints         !   number of aperture points to average over
            
            
            !real(kind=real64)                   ::      aa
            real(kind=real64)                   ::      deltaz,bg,b0
                     
            complex(kind=real64),parameter      ::      I = cmplx( 0.0d0,1.0d0,kind=real64 )
            
            real(kind=real64)                   ::      aperpx
            real(kind=real64),dimension(3)      ::      xmax  
            real(kind=real64)                   ::      modg,modk
            !integer                             ::      Mx,My,Mz
                        
            integer                             ::      ix,iy,iz


            complex(kind=real64),dimension(:,:),allocatable     ::      expigdoty       !   (1:nAtoms) phase factor computed at atom sites
            real(kind=real64),dimension(:),allocatable          ::      rho
            
            type(LinkCell3d)                    ::      lc3d
            type(SimpleSupercell)               ::      super
            real(kind=real64),dimension(3)      ::      lc3d_delta
            integer                             ::      nMax,nAtoms
            real(kind=real64),dimension(:),pointer       ::      gdotr
            real(kind=real64),dimension(:,:),pointer     ::      gdotr_aperture
            integer                             ::      deltaM
            real(kind=real64)                   ::      sigma
            real(kind=real64),dimension(3)      ::      xx,gg
            real(kind=real64),dimension(:),allocatable          ::      r2!,aa_aperture
            integer,dimension(:),allocatable                    ::      id
            integer                             ::      nn,kk          !,ii
            complex(kind=real64)                ::      sumigr
            complex(kind=real64),dimension(:),allocatable                ::      sumigr_aperture
            real(kind=real64)                   ::      i2s2,ww,r2_min 
            real(kind=real64),dimension(3,3)    ::      a_vox, a_column
            
            integer                             ::      Mx,My,Mz            !   number of pixels
            integer                             ::      Cx,Cy               !   number of columns to divide into 
            integer                             ::      Dx,Dy               !   number of voxels in each column, excluding DeltaM buffer
            integer                             ::      Ex,Ey,Ez            
            
            integer                             ::      jx,jy , jj  
            integer                             ::      prec,nPrec,app,nApp
            real(kind=real64)                   ::      phi,costheta,sintheta,ux,uy,d2big2_aperture ,sg ,sg_bar,sg2_bar
            real(kind=real64),dimension(3,3)    ::      BeamTilt
!            real(kind=real64),dimension(2)      ::      apertureShift
            real(kind=real64),dimension(:,:),allocatable      ::      dgAperturePoint,aa
            integer                             ::      nProcs,rank,ierror
            real(kind=real64),dimension(:,:),allocatable      ::  d2big2_tmp
            
            call MPI_COMM_SIZE(MPI_COMM_WORLD, nProcs, ierror)
            call MPI_COMM_RANK(MPI_COMM_WORLD, rank, ierror)
            
            
            
            if (this%test) then
                if (rank == 0) print *,"Lib_DynamicalTwoBeamImaging::computeD2BI_img() info - test mode set, not computing image"
                d2big2 = 0.0d0      !   that'll do
                return
            end if
            
            
            call findImageParameters( this,aperpx,xmax,modg,Mx,My,Mz )
            
            
            allocate( d2big2_tmp(0:Mx-1,0:My-1) )
            
            
        !---    find precession parameters            
            nPrec = 0
            costheta = 1.0d0
            sintheta = 0.0d0
            if (present(nPrecessionAngle)) then
                nPrec=nPrecessionAngle
                if (nPrec>1) then                    
                    costheta = cos(precessionAngle) 
                    sintheta = sin(precessionAngle)
                    if (rank==0) write(*,fmt='(a,f12.6,3(a,f12.4))') "Lib_DynamicalTwoBeamImaging::computeD2BI_img() info - precession angle ",precessionAngle," (rad)"
                end if
            end if 
            
            
        !---    find aperture parameters            
            if (present(nAperturePoints)) then
                nApp = nAperturePoints
                modk = 2*PI/wavelength( this%E*1000 )
                if ( (rank==0).and.(nApp>1) ) write(*,fmt='(a,f12.6,3(a,f12.4))') "Lib_DynamicalTwoBeamImaging::computeD2BI_img() info - aperture angle ",apertureAngle," (rad) gives |g| = ",modg,"+/-",apertureAngle * modk," (1/A)"
                modg = apertureAngle * modk                
            else
                nApp = 1
                modg = 0.0d0
            end if            
            allocate(dgAperturePoint(3,nApp))
            call evenlySpacedPointsInCircle(dgAperturePoint)
            dgAperturePoint = modg * dgAperturePoint
            
            

        !---    version July 2023.
            
        
            nAtoms = size(x,dim=2)
            sigma = this%a0 / 2                                                                                     !   spreading width for atomic phase factor
            i2s2 = 1/(2*sigma*sigma)
            deltaM = ceiling( 3*sigma/min( getCellSideLength(this%gksuper,1),getCellSideLength(this%gksuper,2) ) )  !   number of gksuper cells needed for 3 sigma
            !nMax = ceiling( 2*deltaM*deltaM*nAtoms/superCellVolume(this%gksuper) )
            nMax = estimateMaxCountPerCell( real(Mx,kind=real64)*My*Mz,nAtoms ) *deltaM*deltaM 
            if (rank==0) then
                write(*,fmt='(5(a,i6))')"Lib_DynamicalTwoBeamImaging::computeD2BI_img() info - voxelisation of phase field for image construction"
                call report(this%gksuper)
                write(*,fmt='(a,i10,i6,a,f12.3,i6)')"Lib_DynamicalTwoBeamImaging::computeD2BI_img() info - nAtoms,nMax/cell ",nAtoms,nMax," sigma, deltaM ",sigma,deltaM
            end if
            
            allocate(rho(0:Mz))     
            allocate(expigdoty(0:Mz,nApp))
                        
            nn = 27*nMax
            allocate(r2(nn))                            
            allocate(id(nn))
            
            a_vox = getA( this%gksuper )
            
            
            
            
        !---    start by setting link cell list to cover all gksuper.                
           ! super = SimpleSupercell_ctor( getA( this%gksuper ), Mx,My,Mz )
           ! lc3d = LinkCell3D_ctor(nMax,super)         
           ! lc3d_delta = 0.0d0
           
!             
!             if (Mx>My) then
!                 Cx = max(1,floor( sqrt( real(nProcs*Mx)/My ) ))
!                 Cy = floor( real(nProcs)/Cx )
!             else
!                 Cy = max(1,floor( sqrt( real(nProcs*My)/Mx ) ))
!                 Cx = floor( real(nProcs)/Cy )
!             end if                

            call factoriseParallel( Mx,My , nProcs , Cx,Cy )

            
            
            Dx = ceiling(real(Mx)/Cx)
            Dy = ceiling(real(My)/Cy)
            
            !   example: Mx=12,My=7
            !   Cx = 3,Cy = 2
            !   Dx = 4,Dy = 4
            !   DeltaM = 1
            
            !              0 1 2 3 4 5 6 7 8 9   11         0 1 2 3 4 5 6 7 8 9   11     
            !         _______________                           _______________           
            !        |_|_|_|_|_|_|_|_|                         |   |   |   |   |          
            !        |_|_|_|_|_|_|_|_|                         |___|___|___|___|          
            !    0   |_|_|X|X|X|X|_|_|� � � � � �          � � |   |X X X X|   |� �       
            !    1   |_|_|X|X|X|X|_|_|� � � � � �          � � |___|X X X X|___|� �       
            !    2   |_|_|X|X|X|X|_|_|� � � � � �          � � |   |X X X X|   |� �       
            !    3   |_|_|X|X|X|X|_|_|� � � � � �          � � |___|X_X_X_X|___|� �       
            !    4   |_|_|_|_|_|_|_|_|� � � � � �          � � |   |   |   |   |� �       
            !    5   |_|_|_|_|_|_|_|_|� � � � � �          � � |___|___|___|___|� �       
            !    6        � � � � � � � � � � � �          � � � � � � � � � � � �                 
            !                                                                            
                                                                                       
            !             0 1 2 3 4 5 6 7 8 9   11      
            !                                           
            !                                           
            !    0        � � � � � � � � � � � �       
            !    1    ____�_�_�_�_�_� � � � � � �       
            !    2   |_|_|_|_|_|_|_|_|� � � � � �       
            !    3   |_|_|_|_|_|_|_|_|� � � � � �       
            !    4   |_|_|X|X|X|X|_|_|� � � � � �       
            !    5   |_|_|X|X|X|X|_|_|� � � � � �       
            !    6   |_|_|X|X|X|X|_|_|� � � � � �              
            !        |_|_|_|_|_|_|_|_|                  
            !        |_|_|_|_|_|_|_|_|
                    
            
            Ex = max(3, floor( real(Dx)/deltaM + 2 ))
            Ey = max(3, floor( real(Dy)/deltaM + 2 ))
            Ez = max(3, floor( real(Mz)/deltaM) )
            
            a_column(1:3,1) = (Dx + 2*deltaM)*a_vox(1:3,1)/Ex
            a_column(1:3,2) = (Dy + 2*deltaM)*a_vox(1:3,2)/Ey
           ! a_column(1:3,3) = a_vox(1:3,3)       
            a_column(1:3,3) = Mz*a_vox(1:3,3)/Ez
            if (rank==0) then
                write(*,fmt='(5(a,i6))') "Lib_DynamicalTwoBeamImaging::computeD2BI_img() info - x vox ",Mx," patches ",Cx," width ",Dx," columns/patch ",Ex
                write(*,fmt='(5(a,i6))') "Lib_DynamicalTwoBeamImaging::computeD2BI_img() info - y vox ",My," patches ",Cy," width ",Dy," columns/patch ",Ey
                
                if (Cx * Cy < nProcs) write(*,fmt='(5(a,i6))') "Lib_DynamicalTwoBeamImaging::computeD2BI_img() WARNING - number of columns < number of processes"
            end if            
           
            !super = SimpleSupercell_ctor( a_column, Ex,Ey,Mz )           
            super = SimpleSupercell_ctor( a_column, Ex,Ey,Ez )           
            
            
            lc3d = LinkCell3D_ctor(nMax,super)
            lc3d_delta = 0.0d0
            
            if (nApp>1) then
                nullify(gdotr) 
                nullify(gdotr_aperture)            
                allocate(gdotr_aperture( ceiling(nAtoms*1.5/nProcs),nApp ))            
                               
            else
                nullify(gdotr)            
                allocate(gdotr( ceiling(nAtoms*1.5/nProcs) ))                
            end if     
            allocate(sumigr_aperture(nApp))
            
            ! if (rank==0) then
            !     print *,"Lib_DynamicalTwoBeamImaging::computeD2BI_img() info - link cell list for each processor"
            !     call report(lc3d,6,4)
            ! end if
            
            
        !---                
            !call findGdotR( this, x, lc3d,lc3d_delta,gdotr )    
            
            rho = 1
            d2big2 = 0                                      
            deltaz = xmax(3)/Mz                 !   length of z division
            bg = (PI/this%xig)  
            b0 = (PI/this%xi0)  
            
            
            !aa = (2*PI*this%sg)   
            
            
            allocate(aa(nApp,nPrec))            
            sg_bar = 0.0d0
            sg2_bar = 0.0d0
            do prec = 1,nPrec
            
                if (prec==1) then
                !   prec = 1 is always zero beam tilt.
                    gg = getAperture(this)
                else
                    phi = (prec-1)*2*PI / (nPrec-1)
                    ux = sin(phi)
                    uy = cos(phi)
                    BeamTilt(1,1:3) = (/ costheta + ux*ux*(1-costheta)  ,  ux*uy*(1-costheta)               ,   uy*sintheta     /) 
                    BeamTilt(2,1:3) = (/ ux*uy*(1-costheta)             ,  costheta + uy*uy*(1-costheta)    ,  -ux*sintheta     /) 
                    BeamTilt(3,1:3) = (/ -uy*sintheta                   ,  ux*sintheta                      ,  costheta         /) 
                    gg = getAperture(this,BeamTilt)
                end if
                            
                do app = 1,nApp
                    sg = getDeviationParameter( this,gg+dgAperturePoint(:,app) )
                    aa(app,prec) = (2*PI*sg)   
                    sg_bar = sg_bar + sg
                    sg2_bar = sg2_bar + sg*sg
                    !if (rank==0) write(*,fmt='(a,i6,a,3f10.5,a,f10.5,a)') "aperture ",app+nApp*prec," g = ",gg," s_g = ",sg," (1/A)"
                     
                end do
            end do
            sg_bar = sg_bar / (nApp*nPrec)
            sg2_bar = sg2_bar / (nApp*nPrec)
            if (nApp*nPrec>2) then
                if (rank==0) write(*,fmt='(a,i6,4(a,f10.5))') "number of apertures ",nApp*nPrec," sg = ",sg_bar," +/- ",sqrt( max(0.0d0,sg2_bar - sg_bar*sg_bar) )*sqrt( real(nApp*nPrec)/(nApp*nPrec+1) )," (1/A)"
            end if
            
            
            !print *,"rank ",rank," computeD2BI_img ",nApp,nPrec,sg_bar,aa(1,1)
            
            
            d2big2 = 0.0d0   
            do jj = 0,Cx*Cy-1
                jy = int( jj/Cx )
                jx = jj - jy*Cx 
                
                
                if (mod( jj,nProcs ) == rank) then
                     
                    
                    !   lc3d_delta = cellToRealSpace( this%gksuper,jx*Dx-DeltaM,jy*Dy-DeltaM,0 )
                    lc3d_delta = cellToRealSpace( this%gksuper,jx*Dx,jy*Dy,0 ) - cellToRealSpace( super,1,1,0 )

                           
                    if (nProcs<=16) then    
                        write(*,fmt='(8(a,i6))') "Lib_DynamicalTwoBeamImaging::computeD2BI_img() info - patch ",jx,",",jy," offset vox",jx*Dx,",",jy*Dy," running on process ",rank
                    else if (rank==0) then
                        write(*,fmt='(8(a,i6))') "Lib_DynamicalTwoBeamImaging::computeD2BI_img() info - finding g.r"
                    end if   
                    d2big2 = 0.d0 
                    
                    if (nApp==1) then                                                        
                        call findGdotR( this, x, lc3d,lc3d_delta,gdotr )    
                    else
                        call findGdotR_aperture( this, x, lc3d,lc3d_delta,gdotr_aperture,dgAperturePoint )    
                    end if             
                         
                         
            
                    do iy = jy*Dy,min( (jy+1)*Dy,My)-1                                       !   (ix,iy) is the pixel position...
                        do ix = jx*Dx,min( (jx+1)*Dx,Mx )-1
                        
                            if (rank==0) call progressBar( prec + nPrec*((ix-jx*Dx)+Dx*(iy-jy*Dy)) , nPrec*Dx*Dy )

                            if (nApp == 1) then
                                app = 1
                                do iz = 0,Mz-1                                                    !   ... (ix,iy,iz) is the voxel position in this%gksuper
                                    xx = cellToRealSpace( this%gksuper,ix,iy,iz ) - lc3d_delta
                                    call neighbourList( lc3d,xx,3*sigma, nn,id,r2 )
                                    
                                    r2_min = huge(1.0)
                                    sumigr = 0.0d0
                                    do kk = 1,nn
                                        r2_min = min( r2_min,r2(kk) )               !   shortest scaled range to an atom                  
                                        sumigr = sumigr + exp( cmplx( -r2(kk)*i2s2 , gdotr(id(kk)) , kind=real64 ) )         
                                    end do
                                    if (this%useDensityField) rho( iz ) = convertScaledRangeToRho( r2_min,sigma )  
                                    ww = abs( sumigr )
                                    ww = 1/max( 1.0d-16,ww )
                                    expigdoty( iz,app ) = sumigr*ww     
                                    
                                end do               
                                
                                expigdoty( Mz,app ) = expigdoty( 0,app ) 
                                rho( Mz ) = rho( 0 ) 
                                 
                                do prec = 1,nPrec
                                    call dynamicalTwoBeamIntegrationColumn_phaseRho( aa(app,prec),b0,bg,deltaz,expigdoty(:,app),rho,d2big2_aperture )                                       
                                    d2big2(ix,iy) = d2big2(ix,iy) + d2big2_aperture                             
                                end do                         
                                
                            else
                            
                                
                                do iz = 0,Mz-1                                                        !   ... (ix,iy,iz) is the voxel position in this%gksuper
                                    xx = cellToRealSpace( this%gksuper,ix,iy,iz ) - lc3d_delta
                                    call neighbourList( lc3d,xx,3*sigma, nn,id,r2 )
                                    
                                    r2_min = huge(1.0)
                                    sumigr_aperture = 0.0d0
                                    do kk = 1,nn
                                        r2_min = min( r2_min,r2(kk) )               !   shortest scaled range to an atom                  
                                        do app = 1,nApp
                                            sumigr_aperture(app) = sumigr_aperture(app) + exp( cmplx( -r2(kk)*i2s2 , gdotr_aperture(id(kk),app) , kind=real64 ) )         
                                        end do
                                    end do
                                    if (this%useDensityField) rho( iz ) = convertScaledRangeToRho( r2_min,sigma )  
                                    do app = 1,nApp
                                        ww = abs( sumigr_aperture(app) )
                                        ww = 1/max( 1.0d-16,ww )
                                        expigdoty( iz,app ) = sumigr_aperture(app)*ww   
                                    end do  
                                end do      

                                expigdoty( Mz,: ) = expigdoty( 0,: ) 
                                rho( Mz ) = rho( 0 ) 
                                
                                    
                                do prec = 1,nPrec
                                    do app = 1,nApp
                                        call dynamicalTwoBeamIntegrationColumn_phaseRho( aa(app,prec),b0,bg,deltaz,expigdoty(:,app),rho,d2big2_aperture )                                       
                                        d2big2(ix,iy) = d2big2(ix,iy) + d2big2_aperture                             
                                    end do
                                end do   
                                                                  
                            end if      !   aperture
                            
                            
                        end do      !   pixel ix,iy
                    end do
                             
                                    
                end if          !   MPI rank
                    
                 
            end do              !   jj
                         
             
            
            call MPI_Reduce(d2big2, d2big2_tmp, Mx*My, MPI_DOUBLE_PRECISION, MPI_SUM, 0,MPI_COMM_WORLD,ierror )
            if (rank == 0) d2big2 = d2big2_tmp / (nApp*nPrec)
            
            !   MPI_Op op, int root, MPI_Comm comm)
            
                 
            if (nApp==1) then
                deallocate(gdotr)
            else
                deallocate(gdotr_aperture)
            end if
            deallocate(rho)     
            deallocate(expigdoty)
            deallocate(r2)                            
            deallocate(id)
            call delete(lc3d)
            
         
            return
            
        contains
    !---^^^^^^^^
    
            elemental real(kind=real64) function convertScaledRangeToRho( r2_min,sig )             
        !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
        !       a simple method for finding density from minimum distance to atom
        !       - void tapers at sig and is certain at 2 sig
                real(kind=real64),intent(in)                ::      r2_min  
                real(kind=real64),intent(in)                ::      sig
                
                real(kind=real64)           ::      xx
                
                
                
                if (r2_min < sig*sig) then
                    convertScaledRangeToRho = 1.0d0
                else if (r2_min > 4*sig*sig) then
                    convertScaledRangeToRho = 0.0d0
                else
                    xx = (sqrt(r2_min)/sig - 1.0d0)                                     !   0 at sig and 1 at 2 sig
                    convertScaledRangeToRho = 1 - xx*xx*xx*( 10 - 15*xx + 6*xx*xx )  
                end if
                
                return
            end function convertScaledRangeToRho
                    
                
            
        end subroutine computeD2BI_img
                                                                     
    
      
               
!      
    
!         subroutine dynamicalTwoBeamIntegrationColumn_phase( a,b0,bg,dz, x , phig2 )
!     !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
!     !*      integrate dynamical two beam equations through depth zmax
!     !*      given the phase factor x = exp[ i g.u ]
!     !*      and the scattering vector g ( in reciprocal cell units )
!     !*      and three parameters  b0 = pi / xi_0 , bg = pi / xi_g
!     !*                            a  = - eps / (h v)
!     !*      with xi_g = pi h v / | U_g |
!     !*      eps = hbar^2 ( k+g )^2 / (2m) - hbar^2 k^2/(2m)
!     !*      v = hbar k / m
!     !*
!     !*      return the intensity squared phi2
!     !*
!     !*      phi_0' = i b0 phi_0               + i bg exp[ i g.u ] phi_g
!     !*      phi_g' = i bg exp[ -i g.u ] phi_0 + i (b0+a)          phi_g 
!      
!     !*      on input x = exp[ i g.u ]
!      
!             real(kind=real64),intent(in)                    ::      a,bg,b0,dz
!             complex(kind=real64),dimension(0:),intent(in)   ::      x
!             
!             real(kind=real64),intent(out)                   ::      phig2
!             
!             
!             
!             integer             ::      nz
!              
!             integer             ::      ii
!             
!             
!             real(kind=real64)   ::      phi0r ,phi0i ,phigr ,phigi          !   these are the field evaluated at z ( real and imaginary parts )
!             real(kind=real64)   ::      dphi0r,dphi0i,dphigr,dphigi         !   derivative of field evalulated at z,z+dz/2,z+dz
!             real(kind=real64)   ::      phi0r1,phi0i1,phigr1,phigi1         !   temporary 1 needed for RK4 integration
!             real(kind=real64)   ::      phi0r2,phi0i2,phigr2,phigi2         !   temporary 2 needed for RK4 integration
!             
!             real(kind=real64)   ::      phimag , bgs,bgc
!             
!             
!         !---    find the number of z divisions : note this should be an odd number 
!             nz = size(x)-1
!            
!             
!                   
!         !    next initialise the propagating fields - wlog set real component of phi0 to 1, everything else zero
!             phi0r = 1.0d0 ; phi0i = 0.0d0  ; phigr = 0.0d0 ; phigi = 0.0d0
!             
!             !print *,0.0d0,phi0r,phi0i,phigr,phigi 
!         !   now integrate in z    
!             
!             if (DynamicalTwoBeamImaging_dbg) write(*,fmt='(100f12.5)') 0*dz,real(x(0)),aimag(x(0)),phi0r,phi0i,phigr,phigi,phigr*phigr + phigi*phigi 
!                 
!             
!             do ii = 0,nz-2,2      !   note that we use double steps
!             
!             !---    RK4 algorithm step 1:
!             !   compute derivative at z
!                 bgc = bg*real(x(ii)) ; bgs = bg*aimag(x(ii))
!                 dphi0r =           -  b0*phi0i -   bgs*phigr -   bgc*phigi
!                 dphi0i =  b0*phi0r             +   bgc*phigr -   bgs*phigi 
!                 dphigr = bgs*phi0r - bgc*phi0i               -(b0+a)*phigi
!                 dphigi = bgc*phi0r + bgs*phi0i +(b0+a)*phigr
!                 
!                !  print *,ii*dz,dphi0r,dphi0i,dphigr,dphigi , " , ", phi0r,phi0i,phigr,phigi
!                 
!             !   store temps and update to z+2 dz 
!                 phi0r2 = phi0r
!                 phi0i2 = phi0i
!                 phigr2 = phigr
!                 phigi2 = phigi
!                 
!                 phi0r1 = phi0r + dphi0r * dz/3
!                 phi0i1 = phi0i + dphi0i * dz/3
!                 phigr1 = phigr + dphigr * dz/3
!                 phigi1 = phigi + dphigi * dz/3
!                
!                 phi0r  = phi0r + dphi0r * dz
!                 phi0i  = phi0i + dphi0i * dz
!                 phigr  = phigr + dphigr * dz
!                 phigi  = phigi + dphigi * dz
!                
!                 
!             !---    RK4 algorithm step 2:
!             !   compute derivative at z+dz
!                 bgc = bg*real(x(ii+1)) ; bgs = bg*aimag(x(ii+1))
!                 dphi0r =           -  b0*phi0i -   bgs*phigr -   bgc*phigi 
!                 dphi0i =  b0*phi0r             +   bgc*phigr -   bgs*phigi 
!                 dphigr = bgs*phi0r - bgc*phi0i               -(b0+a)*phigi 
!                 dphigi = bgc*phi0r + bgs*phi0i +(b0+a)*phigr      
!                 
!             !   store temps and recompute z+2dz
!                 phi0r  = phi0r2 + dphi0r * dz
!                 phi0i  = phi0i2 + dphi0i * dz
!                 phigr  = phigr2 + dphigr * dz
!                 phigi  = phigi2 + dphigi * dz
!                                    
!                 phi0r1 = phi0r1 + dphi0r * 2*dz/3
!                 phi0i1 = phi0i1 + dphi0i * 2*dz/3
!                 phigr1 = phigr1 + dphigr * 2*dz/3
!                 phigi1 = phigi1 + dphigi * 2*dz/3
!                
!                 
!             !---    RK4 algorithm step 3:
!             !   recompute derivative at z+2dz 
!                 bgc = bg*real(x(ii+1)) ; bgs = bg*aimag(x(ii+1))
!                 dphi0r =           -  b0*phi0i -   bgs*phigr -   bgc*phigi 
!                 dphi0i =  b0*phi0r             +   bgc*phigr -   bgs*phigi 
!                 dphigr = bgs*phi0r - bgc*phi0i               -(b0+a)*phigi 
!                 dphigi = bgc*phi0r + bgs*phi0i +(b0+a)*phigr
!                          
!             !   store temps and compute z+dz
!                 phi0r  = phi0r2 + dphi0r * 2*dz
!                 phi0i  = phi0i2 + dphi0i * 2*dz
!                 phigr  = phigr2 + dphigr * 2*dz
!                 phigi  = phigi2 + dphigi * 2*dz
!                                    
!                 phi0r1 = phi0r1 + dphi0r * 2*dz/3
!                 phi0i1 = phi0i1 + dphi0i * 2*dz/3
!                 phigr1 = phigr1 + dphigr * 2*dz/3
!                 phigi1 = phigi1 + dphigi * 2*dz/3
!                
!             
!             !---    RK4 algorithm step 4:
!             !   recompute derivative at z+2 dz
!                 bgc = bg*real(x(ii+2)) ; bgs = bg*aimag(x(ii+2))
!                 dphi0r =           -  b0*phi0i -   bgs*phigr -   bgc*phigi 
!                 dphi0i =  b0*phi0r             +   bgc*phigr -   bgs*phigi 
!                 dphigr = bgs*phi0r - bgc*phi0i               -(b0+a)*phigi 
!                 dphigi = bgc*phi0r + bgs*phi0i +(b0+a)*phigr
!                              
!             !   recompute z+dz
!                 phi0r  = phi0r1 + dphi0r * dz/3
!                 phi0i  = phi0i1 + dphi0i * dz/3
!                 phigr  = phigr1 + dphigr * dz/3
!                 phigi  = phigi1 + dphigi * dz/3
!                 
!                                              
!                 !phi0r1 = (phi0r+phigr)
!                 !phi0i1 = (phi0i+phigi)
!                       
!                  if(DynamicalTwoBeamImaging_dbg) write(*,fmt='(100f12.5)') (ii+2)*dz,real(x(ii+2)),aimag(x(ii+2)),phi0r,phi0i,phigr,phigi,phigr*phigr + phigi*phigi 
!                 
!                 
!             !---    this normalisation step is probably unnecessary - in testing I can't see any difference to 
!             !       several d.p. But it can't hurt...                 
!                   phimag =  phi0r*phi0r +  phi0i*phi0i + phigr*phigr +  phigi*phigi
!                   phimag = 1/sqrt(phimag)
!                   
!                   phi0r = phi0r * phimag
!                   phi0i = phi0i * phimag
!                   phigr = phigr * phimag
!                   phigi = phigi * phimag
!                
!             end do           
!                          
!         !---    finally compute the intensity of the diffracted beam
!                phig2 = phigr*phigr + phigi*phigi
!             
!             
!             return
!         end subroutine dynamicalTwoBeamIntegrationColumn_phase
!                 
    
!         subroutine dynamicalTwoBeamIntegrationColumn_phaserange( a,w,b0,bg,dz, x , phig2 )
!     !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
!     !*      integrate dynamical two beam equations through depth zmax
!     !*      given the phase factor x = exp[ i g.u ]
!     !*      and the scattering vector g ( in reciprocal cell units )
!     !*      and three parameters  b0 = pi / xi_0 , bg = pi / xi_g
!     !*                            a  = - eps / (h v)
!     !*      with xi_g = pi h v / | U_g |
!     !*      eps = hbar^2 ( k+g )^2 / (2m) - hbar^2 k^2/(2m)
!     !*      v = hbar k / m
!     !*
!     !*      return the intensity squared phi2
!     !*
!     !*      phi_0' = i b0 phi_0               + i bg exp[ i g.u ] phi_g
!     !*      phi_g' = i bg exp[ -i g.u ] phi_0 + i (b0+a)          phi_g 
!      
!     !*      on input x = exp[ i g.u ]
!                                                                      
!             real(kind=real64),dimension(:),intent(in)       ::      a,w
!             real(kind=real64),intent(in)                    ::      bg,b0,dz
!             complex(kind=real64),dimension(0:),intent(in)   ::      x
!             
!             real(kind=real64),intent(out)                   ::      phig2
!             
!             
!             
!             integer             ::      nz
!              
!             integer             ::      ii 
!             
!             
!             real(kind=real64),dimension(size(a))   ::      phi0r ,phi0i ,phigr ,phigi          !   these are the field evaluated at z ( real and imaginary parts )
!             real(kind=real64),dimension(size(a))   ::      dphi0r,dphi0i,dphigr,dphigi         !   derivative of field evalulated at z,z+dz/2,z+dz
!             real(kind=real64),dimension(size(a))   ::      phi0r1,phi0i1,phigr1,phigi1         !   temporary 1 needed for RK4 integration
!             real(kind=real64),dimension(size(a))   ::      phi0r2,phi0i2,phigr2,phigi2         !   temporary 2 needed for RK4 integration
!             
!             real(kind=real64)                           ::      bgs,bgc,ww
!             
!             
!         !---    find the number of z divisions : note this should be an odd number 
!             nz = size(x)-1
!             
!             
!                   
!         !    next initialise the propagating fields - wlog set real component of phi0 to 1, everything else zero
!             phi0r = 1.0d0 ; phi0i = 0.0d0  ; phigr = 0.0d0 ; phigi = 0.0d0
!             
!             !print *,0.0d0,phi0r,phi0i,phigr,phigi 
!         !   now integrate in z    
!             
!               
!             do ii = 0,nz-2,2      !   note that we use double steps
!             
!             !---    RK4 algorithm step 1:
!             !   compute derivative at z
!                 bgc = bg*real(x(ii)) ; bgs = bg*aimag(x(ii))
!                 dphi0r(:) =              -  b0*phi0i(:) -   bgs*phigr(:) -   bgc*phigi(:)
!                 dphi0i(:) =  b0*phi0r(:)                +   bgc*phigr(:) -   bgs*phigi(:) 
!                 dphigr(:) = bgs*phi0r(:) - bgc*phi0i(:)               -(b0+a(:))*phigi(:)
!                 dphigi(:) = bgc*phi0r(:) + bgs*phi0i(:) +(b0+a(:))*phigr(:)
!                 
!                !  print *,ii*dz,dphi0r,dphi0i,dphigr,dphigi , " , ", phi0r,phi0i,phigr,phigi
!                 
!             !   store temps and update to z+2 dz 
!                 phi0r2(:) = phi0r(:)
!                 phi0i2(:) = phi0i(:)
!                 phigr2(:) = phigr(:)
!                 phigi2(:) = phigi(:)
!                 
!                 phi0r1(:) = phi0r(:) + dphi0r(:) * dz/3
!                 phi0i1(:) = phi0i(:) + dphi0i(:) * dz/3
!                 phigr1(:) = phigr(:) + dphigr(:) * dz/3
!                 phigi1(:) = phigi(:) + dphigi(:) * dz/3
!                
!                 phi0r(:)  = phi0r(:) + dphi0r(:) * dz
!                 phi0i(:)  = phi0i(:) + dphi0i(:) * dz
!                 phigr(:)  = phigr(:) + dphigr(:) * dz
!                 phigi(:)  = phigi(:) + dphigi(:) * dz
!                
!                 
!             !---    RK4 algorithm step 2:
!             !   compute derivative at z+dz
!                 bgc = bg*real(x(ii+1)) ; bgs = bg*aimag(x(ii+1))
!                 dphi0r(:) =              -  b0*phi0i(:) -   bgs*phigr(:) -   bgc*phigi(:)
!                 dphi0i(:) =  b0*phi0r(:)                +   bgc*phigr(:) -   bgs*phigi(:) 
!                 dphigr(:) = bgs*phi0r(:) - bgc*phi0i(:)               -(b0+a(:))*phigi(:)
!                 dphigi(:) = bgc*phi0r(:) + bgs*phi0i(:) +(b0+a(:))*phigr(:)
!                 
!             !   store temps and recompute z+2dz
!                 phi0r(:)  = phi0r2(:) + dphi0r(:) * dz
!                 phi0i(:)  = phi0i2(:) + dphi0i(:) * dz
!                 phigr(:)  = phigr2(:) + dphigr(:) * dz
!                 phigi(:)  = phigi2(:) + dphigi(:) * dz
!                                    
!                 phi0r1(:) = phi0r1(:) + dphi0r(:) * 2*dz/3
!                 phi0i1(:) = phi0i1(:) + dphi0i(:) * 2*dz/3
!                 phigr1(:) = phigr1(:) + dphigr(:) * 2*dz/3
!                 phigi1(:) = phigi1(:) + dphigi(:) * 2*dz/3
!                
!                 
!             !---    RK4 algorithm step 3:
!             !   recompute derivative at z+2dz 
!                 bgc = bg*real(x(ii+1)) ; bgs = bg*aimag(x(ii+1))
!                 dphi0r(:) =              -  b0*phi0i(:) -   bgs*phigr(:) -   bgc*phigi(:)
!                 dphi0i(:) =  b0*phi0r(:)                +   bgc*phigr(:) -   bgs*phigi(:) 
!                 dphigr(:) = bgs*phi0r(:) - bgc*phi0i(:)               -(b0+a(:))*phigi(:)
!                 dphigi(:) = bgc*phi0r(:) + bgs*phi0i(:) +(b0+a(:))*phigr(:)
!                          
!             !   store temps and compute z+dz
!                 phi0r(:)  = phi0r2(:) + dphi0r(:) * 2*dz
!                 phi0i(:)  = phi0i2(:) + dphi0i(:) * 2*dz
!                 phigr(:)  = phigr2(:) + dphigr(:) * 2*dz
!                 phigi(:)  = phigi2(:) + dphigi(:) * 2*dz
!                                    
!                 phi0r1(:) = phi0r1(:) + dphi0r(:) * 2*dz/3
!                 phi0i1(:) = phi0i1(:) + dphi0i(:) * 2*dz/3
!                 phigr1(:) = phigr1(:) + dphigr(:) * 2*dz/3
!                 phigi1(:) = phigi1(:) + dphigi(:) * 2*dz/3
!                
!             
!             !---    RK4 algorithm step 4:
!             !   recompute derivative at z+2 dz
!                 bgc = bg*real(x(ii+2)) ; bgs = bg*aimag(x(ii+2))
!                 dphi0r(:) =              -  b0*phi0i(:) -   bgs*phigr(:) -   bgc*phigi(:)
!                 dphi0i(:) =  b0*phi0r(:)                +   bgc*phigr(:) -   bgs*phigi(:) 
!                 dphigr(:) = bgs*phi0r(:) - bgc*phi0i(:)               -(b0+a(:))*phigi(:)
!                 dphigi(:) = bgc*phi0r(:) + bgs*phi0i(:) +(b0+a(:))*phigr(:)
!                              
!             !   recompute z+dz
!                 phi0r(:)  = phi0r1(:) + dphi0r(:) * dz/3
!                 phi0i(:)  = phi0i1(:) + dphi0i(:) * dz/3
!                 phigr(:)  = phigr1(:) + dphigr(:) * dz/3
!                 phigi(:)  = phigi1(:) + dphigi(:) * dz/3
!                 
!                
!             end do           
!                          
!         !---    finally compute the intensity of the diffracted beam
!             ww = 0.0d0
!             phig2 = 0.0d0
!             do ii = 1,size(a)
!                 ww = ww + w(ii)
!                 phig2 = phig2 + w(ii)*( phigr(ii)*phigr(ii) + phigi(ii)*phigi(ii) )
!             end do
!             phig2 = phig2/ww
!             
!             return
!         end subroutine dynamicalTwoBeamIntegrationColumn_phaserange
!                 
        subroutine dynamicalTwoBeamIntegrationColumn_phaserho( a,b0,bg,dz,x,rho , phig2 )
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      integrate dynamical two beam equations through depth zmax
    !*      given the phase factor x = exp[ i g.u ]
    !*      and the scattering vector g ( in reciprocal cell units )
    !*      and three parameters  b0 = pi / xi_0 , bg = pi / xi_g
    !*                            a  = - eps / (h v)
    !*      with xi_g = pi h v / | U_g |
    !*      eps = hbar^2 ( k+g )^2 / (2m) - hbar^2 k^2/(2m)
    !*      v = hbar k / m
    !*
    !*      return the intensity squared phi2
    !*
    !*      phi_0' = i b0 phi_0               + i bg exp[ i g.u ] phi_g
    !*      phi_g' = i bg exp[ -i g.u ] phi_0 + i (b0+a)          phi_g 
     
    !*      on input x = exp[ i g.u ] 
    !*      uses the density field rho
    !*
     
            real(kind=real64),intent(in)                    ::      a,bg,b0,dz
            complex(kind=real64),dimension(0:),intent(in)   ::      x
            real(kind=real64),dimension(0:),intent(in)      ::      rho
            
            real(kind=real64),intent(out)                   ::      phig2
            
            
            
            integer             ::      nz
            !real(kind=real64)   ::      dz
            integer             ::      ii
            
            
            real(kind=real64)   ::      phi0r ,phi0i ,phigr ,phigi          !   these are the field evaluated at z ( real and imaginary parts )
            real(kind=real64)   ::      dphi0r,dphi0i,dphigr,dphigi         !   derivative of field evalulated at z,z+dz/2,z+dz
            real(kind=real64)   ::      phi0r1,phi0i1,phigr1,phigi1         !   temporary 1 needed for RK4 integration
            real(kind=real64)   ::      phi0r2,phi0i2,phigr2,phigi2         !   temporary 2 needed for RK4 integration
            
            real(kind=real64),dimension(0:size(x)-1)    ::      c0,cgr,cgi,ca
            
            real(kind=real64)   ::      phimag !, bgs,bgc
            
        !---    find the number of z divisions : note this should be an odd number 
            nz = size(x)-1
              
!       LINES 29/11/22            
!        !       find the input b0,bg scaled by density
!            do ii = 0,nz
!                c0(ii)  = b0*rho(ii)
!                !  OLD LINE 140922                
!!                ca(ii)  = c0(ii) + a
!               !   NEW LINE 
!               ca(ii) = (b0 + a)*rho(ii)
!                cgr(ii) = bg*rho(ii)*real(x(ii))
!                cgi(ii) = bg*rho(ii)*aimag(x(ii)) 
!            end do                    
!            
        !       find the input b0,bg scaled by density
!       LINES 30/11/22        
            do ii = 0,nz
                c0(ii)  = b0*rho(ii)  
                ca(ii)  = b0*rho(ii) + a
                cgr(ii) = bg*rho(ii)*real(x(ii))
                cgi(ii) = bg*rho(ii)*aimag(x(ii)) 
            end do                    
            
            
            
        !    next initialise the propagating fields - wlog set real component of phi0 to 1, everything else zero
            phi0r = 1.0d0 ; phi0i = 0.0d0  ; phigr = 0.0d0 ; phigi = 0.0d0
            
            
        !   now integrate in z        
            do ii = 0,nz-2,2      !   note that we use double steps
            
            !---    RK4 algorithm step 1:
            !   compute derivative at z
                dphi0r =               -  c0(ii)*phi0i -   cgi(ii)*phigr -   cgr(ii)*phigi
                dphi0i =  c0(ii)*phi0r                 +   cgr(ii)*phigr -   cgi(ii)*phigi 
                dphigr = cgi(ii)*phi0r - cgr(ii)*phi0i                   -    ca(ii)*phigi
                dphigi = cgr(ii)*phi0r + cgi(ii)*phi0i +   ca(ii)*phigr
                
            !   store temps and update to z+dz 
                phi0r2 = phi0r                 
                phi0i2 = phi0i
                phigr2 = phigr
                phigi2 = phigi
                
                phi0r1 = phi0r + dphi0r * dz/3
                phi0i1 = phi0i + dphi0i * dz/3
                phigr1 = phigr + dphigr * dz/3
                phigi1 = phigi + dphigi * dz/3
                                             
                phi0r  = phi0r + dphi0r * dz 
                phi0i  = phi0i + dphi0i * dz 
                phigr  = phigr + dphigr * dz 
                phigi  = phigi + dphigi * dz 
               
                
            !---    RK4 algorithm step 2:
            !   compute derivative at z+dz 
                dphi0r =                 -  c0(ii+1)*phi0i -   cgi(ii+1)*phigr -   cgr(ii+1)*phigi
                dphi0i =  c0(ii+1)*phi0r                   +   cgr(ii+1)*phigr -   cgi(ii+1)*phigi 
                dphigr = cgi(ii+1)*phi0r - cgr(ii+1)*phi0i                     -    ca(ii+1)*phigi
                dphigi = cgr(ii+1)*phi0r + cgi(ii+1)*phi0i +   ca(ii+1)*phigr
                
            !   store temps and recompute z+dz 
                phi0r  = phi0r2 + dphi0r * dz 
                phi0i  = phi0i2 + dphi0i * dz 
                phigr  = phigr2 + dphigr * dz 
                phigi  = phigi2 + dphigi * dz 
                                   
                phi0r1 = phi0r1 + dphi0r * 2*dz/3
                phi0i1 = phi0i1 + dphi0i * 2*dz/3
                phigr1 = phigr1 + dphigr * 2*dz/3
                phigi1 = phigi1 + dphigi * 2*dz/3
               
                
            !---    RK4 algorithm step 3:
            !   recompute derivative at z+dz 
                dphi0r =                 -  c0(ii+1)*phi0i -   cgi(ii+1)*phigr -   cgr(ii+1)*phigi
                dphi0i =  c0(ii+1)*phi0r                   +   cgr(ii+1)*phigr -   cgi(ii+1)*phigi 
                dphigr = cgi(ii+1)*phi0r - cgr(ii+1)*phi0i                     -    ca(ii+1)*phigi
                dphigi = cgr(ii+1)*phi0r + cgi(ii+1)*phi0i +   ca(ii+1)*phigr
                         
            !   store temps and compute z+2dz
                phi0r  = phi0r2 + dphi0r * 2*dz
                phi0i  = phi0i2 + dphi0i * 2*dz
                phigr  = phigr2 + dphigr * 2*dz
                phigi  = phigi2 + dphigi * 2*dz
                                   
                phi0r1 = phi0r1 + dphi0r * 2*dz/3
                phi0i1 = phi0i1 + dphi0i * 2*dz/3
                phigr1 = phigr1 + dphigr * 2*dz/3
                phigi1 = phigi1 + dphigi * 2*dz/3
               
            
            !---    RK4 algorithm step 4:
            !   recompute derivative at z+2dz
                dphi0r =                 -  c0(ii+2)*phi0i -   cgi(ii+2)*phigr -   cgr(ii+2)*phigi
                dphi0i =  c0(ii+2)*phi0r                   +   cgr(ii+2)*phigr -   cgi(ii+2)*phigi 
                dphigr = cgi(ii+2)*phi0r - cgr(ii+2)*phi0i                     -    ca(ii+2)*phigi
                dphigi = cgr(ii+2)*phi0r + cgi(ii+2)*phi0i +   ca(ii+2)*phigr
                             
            !   recompute z+dz
                phi0r  = phi0r1 + dphi0r * dz/3
                phi0i  = phi0i1 + dphi0i * dz/3
                phigr  = phigr1 + dphigr * dz/3
                phigi  = phigi1 + dphigi * dz/3
                                                        
                phimag =  phi0r*phi0r +  phi0i*phi0i + phigr*phigr +  phigi*phigi
                phimag = 1/sqrt(phimag)                                          
              
                phi0r = phi0r * phimag
                phi0i = phi0i * phimag
                phigr = phigr * phimag
                phigi = phigi * phimag
                
            end do           
                         
        !---    finally compute the intensity of the diffracted beam
            phig2 = phigr*phigr + phigi*phigi
            
            return
        end subroutine dynamicalTwoBeamIntegrationColumn_phaserho
                
    
    
    end module Lib_DynamicalTwoBeamImaging
    
    
!!   gfortran -ffree-line-length-256 -Og -g Lib_DynamicalTwoBeamImaging.f90 -o testDynamicalTwoBeamImagings.exe
!
!    
!    program testDynamicalTwoBeamImagings
!!---^^^^^^^^^^^^^^^^^^^^^^^^
!        use iso_fortran_env
!        use Lib_DynamicalTwoBeamImaging
!        implicit none
!        
!        type(DynamicalTwoBeamImaging)           ::      this
!        
!        this = DynamicalTwoBeamImaging_ctor()
!        
!        call report(this)
!        
!        call delete(this)
!        
!        print *,""
!        print *,"done"
!        print *,""
!        
!    end program testDynamicalTwoBeamImagings
!    
!        
    
        
        
    