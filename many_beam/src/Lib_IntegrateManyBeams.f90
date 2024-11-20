
    module Lib_IntegrateManyBeams
!---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
!*      works with one or more ManyBeam objects and integrates electron beams in parallel
!*      works with one or more ManyBeam objects and integrates electron beams in parallel
!*
!*      order to construct
!*          
!*          call setLattice( this,a0_in,latticename,rho_in )            !   determines the set of g-vectors and the lattice parameter
!*          call readInputXyzFile( this,filename,dfg_bar )              !   reads in the atoms, determines the atom space and finds the average deformation gradient
!*              call setR(this%gv,dfg_bar)                                    
!*              mb = ManyBeams_ctor()
!*          call setImagingSpace(this)                                  !   sets the imaging space according to the atom space and dispersion angle
!*          call computePhaseFields(this)                               !   computes the phase factor and atom density


!*      Daniel Mason
!*      (c) UKAEA Sept 2024
!*

        use iso_fortran_env
        use Lib_AtomSpace
        use Lib_ImagingSpace
        use Lib_Lattices
        use Lib_Gvectors
        use Lib_FactoriseParallel
        use Lib_ComputePhaseFactor
        use Lib_ManyBeam
        use Lib_XYZFiles 
        use Lib_SimpleProgressBar
        use NBAX_StringTokenizers
        use Lib_RelativisticElectrons
        use Lib_RotationMatrices
        use Lib_Elements
        use Lib_ReadExtinctionDistances
        use Lib_DeformationGradients
        use Lib_FibonacciSphere
        use Lib_CrystalStructureFactor        
        use Lib_RK4
        use Lib_Filenames
#ifdef MPI
        use mpi_f08
#endif
        implicit none
        private

        integer,private,parameter                       ::      SELECTION_RULE = kind( 1_int32 )
        integer(kind=SELECTION_RULE),public,parameter   ::      LIB_IMB_KEEP_ALL_GVEC = 0_SELECTION_RULE
        integer(kind=SELECTION_RULE),public,parameter   ::      LIB_IMB_BLOCK_HIGH_SG = 1_SELECTION_RULE
        integer(kind=SELECTION_RULE),public,parameter   ::      LIB_IMB_ZOLZ          = 2_SELECTION_RULE            !   keep only zeroth order laue diffraction spots
        integer(kind=SELECTION_RULE),public,parameter   ::      LIB_IMB_FOLZ          = 3_SELECTION_RULE            !   keep 0th and 1st order laue diffraction spots
        integer(kind=SELECTION_RULE),public,parameter   ::      LIB_IMB_SOLZ          = 4_SELECTION_RULE            !   keep 0th and 1st and 2nd order laue diffraction spots


        integer,private,parameter                       ::      STATUS = kind( 1_int32 )
        integer(kind=STATUS),private,parameter          ::      LIB_IMB_STATUS_UNSET    = 0_STATUS
        integer(kind=STATUS),private,parameter          ::      LIB_IMB_FILE_READ       = 1_STATUS
        integer(kind=STATUS),private,parameter          ::      LIB_IMB_ATOMSPACE_SET   = 2_STATUS
        integer(kind=STATUS),private,parameter          ::      LIB_IMB_GVECTORS_SET    = 4_STATUS
        integer(kind=STATUS),private,parameter          ::      LIB_IMB_MANYBEAM_SET    = 8_STATUS
        integer(kind=STATUS),private,parameter          ::      LIB_IMB_CSTRUCTFACT_SET = 16_STATUS
        integer(kind=STATUS),private,parameter          ::      LIB_IMB_IMAGESPACE_SET  = 32_STATUS
        integer(kind=STATUS),private,parameter          ::      LIB_IMB_SIM_RESULT      = 128_STATUS
        

        logical,public                      ::      LIB_IMB_DBG = .false.
        real(kind=real64),parameter         ::      PI = 3.14159265390d0

        integer,private                     ::      rank = 0, nProcs = 1

        public          ::      Lib_IntegrateManyBeams_init_MPI

        public          ::      IntegrateManyBeams_ctor
        public          ::      report
        public          ::      delete
        public          ::      integrate
        public          ::      getImagingSpace
        public          ::      getGvectors
        public          ::      getIntensity
        public          ::      getA
        public          ::      getSigma
        public          ::      selectFoilTilt
        public          ::      applyFoilTilt
        public          ::      perfectLatticeIntensity
        public          ::      changeGvectors
        public          ::      setImagingSpace        
        public          ::      computePhaseFields
        public          ::      sanityCheck



        type,public     ::      phi_slice
            private
            complex(kind=real64),dimension(:,:,:),pointer       ::      phi                 !   (0:nG,lbx:ubx,lby:uby ) the solution 
        end type

        type,public     ::      IntegrateManyBeams
            private
            
            type(AtomSpace)                                     ::      as
            type(ImagingSpace)                                  ::      is
            type(Gvectors),pointer                              ::      gv
            
            real(kind=real64)                                   ::      T                   !   calculation temperature
            real(kind=real64)                                   ::      V                   !   electron acceleration voltage (V)
            complex(kind=real64),dimension(:),pointer           ::      xi                  !   (1:nG2) complex extinction distances for all vectors g" = g - g'
            integer                                             ::      nPrec               !   how many sets used for precession averaging?
            real(kind=real64)                                   ::      precAngle           !   precession angle (rad)
            real(kind=real64)                                   ::      a , sigma           !   cell side , imaging blur
            real(kind=real64)                                   ::      a0                  !   characteristic lattice size used to scale g-vectors
            integer                                             ::      M = 4               !   a = a0/M
            integer                                             ::      zlayers = 20        !   allocate rho,grad_arg_x in blocks of Zlayers

            integer                                             ::      nAtoms
            character(len=XYZFILE_ATOMNAMELENGTH)               ::      element
            type(Lattice)                                       ::      latt
            real(kind=real64),dimension(3,3)                    ::      dfg_bar             !   mean deformation gradient, used to set g-vectors 
            real(kind=real32),dimension(:,:),pointer            ::      r                   !   (3,nAtoms) positions of atoms, read by readInputXyzFile(), converted to imaging space units by setImagingSpace()
            real(kind=real64)                                   ::      omega0              !   volume per atom
            
            logical                                             ::      columnar            !   use columnar approximation
            logical                                             ::      lossy               !   use imaginary parts of crystal srtucture factor
            integer(kind=STATUS)                                ::      status              !   status of calculation

            type(ManyBeam),dimension(:),pointer                 ::      mb
            real(kind=real32),dimension(:,:,:,:,:),pointer      ::      grad_arg_x          !   (3,nG,lbx:ubx,lby:uby,lbz:ubz) pointer to chunk of phase factor gradients        
            real(kind=real32),dimension(:,:,:),pointer          ::      rho                 !   (lbx:ubx,lby:uby,lbz:ubz) pointer to chunk of atom densities            
            type(phi_slice),dimension(:),pointer                ::      slice               !   (1:nPrec) the solution

        end type


        interface       IntegrateManyBeams_ctor
            module procedure            IntegrateManyBeams_null
            module procedure            IntegrateManyBeams_ctor1
        end interface

        interface       report
            module procedure            report0
        end interface

        interface       delete
            module procedure            delete0
        end interface

        interface       integrate
            module procedure            integrate0
        end interface


        interface       getIntensity
            module procedure            getIntensity0
            module procedure            getIntensity1
        end interface

        interface       perfectLatticeIntensity
            module procedure            perfectLatticeIntensity0
        end interface
        

        interface       getSigma
            module procedure            getSigma0
        end interface

        interface       getA
            module procedure            getA0
        end interface

        interface       changeGvectors
            module procedure            changeGvectors0
            module procedure            changeGvectors1
        end interface
        

        interface       sanityCheck
            module procedure            sanityCheck0
        end interface
        

    contains
!---^^^^^^^^

        function IntegrateManyBeams_null() result(this)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      empty constructor
            type(IntegrateManyBeams)          ::      this
            this%latt = Lattice_ctor()
            this%as = AtomSpace_ctor()
            this%is = ImagingSpace_ctor()
            this%nPrec = 0
            this%element = ""
            nullify(this%mb)
            this%T = 300.0d0
            this%V = 200000d0
            this%dfg_bar = reshape( (/1,0,0,0,1,0,0,0,1/),(/3,3/) )
            nullify(this%xi)
            nullify(this%gv)
            nullify(this%grad_arg_x)
            nullify(this%rho)
            this%nAtoms = 0
            this%a0 = 0
            this%omega0 = 0
            nullify(this%r)
            this%columnar = .true.
            this%lossy = .true.
            nullify(this%slice)
            this%status = LIB_IMB_STATUS_UNSET
            return
        end function IntegrateManyBeams_null



        function IntegrateManyBeams_ctor1( latticename,a0_in,filename,T,V,nPrec,precAngle,columnar,lossy,rho_in  ) result(this)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^                  
    !*      default constructor. Intended to read the .xyz input file and determine basic parameters of the calculation.
    !*      Sets a default large number of g-vectors and no foil tilt.
    !*      Sets the default atom space assuming no foil tilt.
    !*      Does not compute the optimal foil tilt
    !*      Does not set the final short list of g-vectors
    !*      Does not set the imaging space 
    !*      Does not compute the phase factors
    !        
            character(len=*),intent(inout)              ::      latticename             !   eg "bcc" the name of the lattice type 
            real(kind=real64),dimension(3),intent(inout)::      a0_in                   !   lattice parameter(s) - may be set to defaults if a0_in = 0
            character(len=*),intent(in)                 ::      filename                !   atom position file
            real(kind=real64),intent(in)                ::      T                       !   temperature (K)
            real(kind=real64),intent(in)                ::      V                       !   accelerating voltage (kV)
            integer,intent(in)                          ::      nPrec                   !   number of many beam copies to use for precession
            real(kind=real64),intent(in)                ::      precAngle               !   precession angle (rad)
            logical,intent(in),optional                 ::      columnar                !   use columnar approximation
            logical,intent(in),optional                 ::      lossy                   !   use imaginary parts of crystal structure factor
            real(kind=real64),intent(in),optional       ::      rho_in                  !   radius of reflections to select
            type(IntegrateManyBeams)                    ::      this

            real(kind=real64),dimension(3)              ::      xyz_offset              !   offset of atom positions read from xyz file
            real(kind=real64),dimension(3,3)            ::      a_super                 !   periodic supercell lattice vectors
            integer                                     ::      ii  
            !real(kind=real64)                           ::      a0                      !   characteristic lattice lengthscale = a0_in(1)

            this = IntegrateManyBeams_null()

        !---    store some basic facts about the calculation            
            if (present(columnar)) this%columnar = columnar
            if (present(lossy)) this%lossy = lossy            
            this%T = T
            this%V = V * 1000
            this%precAngle = precAngle
            this%nPrec = max(1,nPrec)
            


        !---    read input file            
            if (rank==0) then
                print *,""
                print *,"reading input atomic positions file"
                print *,"^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^"
            end if
            call readInputXyzFile( this,filename,xyz_offset,a_super )                
            ii = whichElement(this%element)
            if (ii == 0) then
                if (rank==0) print *,"Lib_IntegrateManyBeams::IntegrateManyBeams_ctor1 warning - do not recognise element """//trim(this%element)//""" - expect problems"
            end if


        !---    set the default lattice parameter, if not already set            
            if (any(a0_in<=0)) then            
                a0_in = getLatticeConstant(this%element)
                if (whichElement(this%element) == 0) then
                    call errorExit("error - can't determine lattice parameter from element name """//trim(this%element)//"""")
                else
                    if (rank==0) write(*,fmt='(a,3f12.5,a)') " Lib_IntegrateManyBeams::IntegrateManyBeams_ctor1 info - setting lattice parameter ",a0_in," from element name"
                end if
            end if
            this%a0 = a0_in(1)                      

        !---    set the default lattice type, if not already set
            if (latticename == UNKNOWN_LATTICE) then
                latticename = getLatticeName(this%element)
                if (latticename == UNKNOWN_LATTICE) then
                    call errorExit("error - can't determine lattice type from element name """//trim(this%element)//"""")
                else
                    if (rank==0) print *,"Lib_IntegrateManyBeams::IntegrateManyBeams_ctor1 info - setting lattice type """//trim(latticename)//""" from element name"
                end if
            end if
            this%latt = Lattice_ctor(latticename)            
            if (getLatticeType(this%latt) == LATTICE_HCP) call setCoverA(this%latt,a0_in(3)/a0_in(1))                                         
            this%omega0 = getOmega0(this%latt) * this%a0**3
            call setFileRead(this,.true.)
 

        !---    set the default atom space assuming no foil tilt
            this%as = AtomSpace_ctor( geta(this) , xyz_offset,a_super)
            call setR(this%as,RotationMatrix_identity)          
            call setThickness(this%as)   
            call setAtomSpaceSet(this,.true.)
!            if (rank==0) call report(this%as)            
            !if (rank==0) print *,""


        !---    set the default g-vectors, assuming no foil tilt
            if (rank==0) then
                print *,""
                print *,"setting g-vectors"
                print *,"^^^^^^^^^^^^^^^^^"
            end if
            allocate(this%gv)
            call changeGvectors( this,rule = LIB_IMB_FOLZ, rho_in = rho_in )  
             
            return
        end function IntegrateManyBeams_ctor1



        subroutine delete0(this)
    !---^^^^^^^^^^^^^^^^^^^^^^^^
            type(IntegrateManyBeams),intent(inout)       ::      this
            integer             ::      mm
            if (this%nPrec==0) return
            call delete(this%latt)
            call delete(this%as)
            call delete(this%is)
            call delete(this%gv)
            deallocate(this%gv)
            deallocate(this%xi)
            do mm = 1,this%nPrec
                call delete(this%mb(mm))
                deallocate(this%slice(mm)%phi)
            end do
            deallocate(this%slice)
            deallocate(this%mb)
      
            if (this%nAtoms>0)  deallocate(this%r)
          
            
            if (associated(this%grad_arg_x)) deallocate(this%grad_arg_x)            
            if (associated(this%rho)) deallocate(this%rho)
            this = IntegrateManyBeams_null()
            return
        end subroutine delete0

 

        
        subroutine report0(this,u,o)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^
            type(IntegrateManyBeams),intent(in)          ::      this
            integer,intent(in),optional         ::      u,o
            integer     ::      uu,oo
            integer     ::      mm
            uu = 6 ; if (present(u)) uu = u
            oo = 0 ; if (present(o)) oo = o
            write(unit=uu,fmt='(a,i6,a)') repeat(" ",oo)//"IntegrateManyBeams [nPrec=",this%nPrec,"]"
            write(unit=uu,fmt='(a,f10.2,a)') repeat(" ",oo+4)//"Temperature      : ",this%T," (K)"
            write(unit=uu,fmt='(a,f10.2,a)') repeat(" ",oo+4)//"Electron voltage : ",this%V," (keV)"
            call report(this%latt,uu,oo+4)
            call report(this%gv,uu,oo+4)
            call report(this%as,uu,oo+4)
            call report(this%is,uu,oo+4)
            do mm = 1,this%nPrec
                call report(this%mb(mm),uu,oo+4)
            end do
            return
        end subroutine report0

    !---


        subroutine changeGvectors0(this,hkl_in)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ 
    !*      given the gvectors may or may not be allocated, 
    !*      update the list of g-vectors to a set containing out to |hkl| <= rho_in
    !*      without changing the input atom positions or lattice.
    !*      This subroutine _will_ construct the crystal structure factors and the manybeams 
    !*      This subroutine _will not_ construct the imaging space, phase field/density, or allocate memory for solution
    !*      special case- set two beam
            type(IntegrateManyBeams),intent(inout)      ::      this
            integer,dimension(:),intent(in)             ::      hkl_in
            integer,dimension(size(hkl_in,dim=1),1)     ::      hkl_out
            hkl_out(:,1) = hkl_in(:)
            call changeGvectors1(this,rule = LIB_IMB_KEEP_ALL_GVEC,hkl_in = hkl_out)
            return
        end subroutine changeGvectors0

        subroutine changeGvectors1(this,rule,rho_in,hkl_in)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      given the gvectors may or may not be allocated, 
    !*      update the list of g-vectors to a set containing out to |hkl| <= rho_in
    !*      without changing the input atom positions or lattice.
    !*      This subroutine _will_ construct the crystal structure factors and the manybeams 
    !*      This subroutine _will not_ construct the imaging space, phase field/density, or allocate memory for solution
            type(IntegrateManyBeams),intent(inout)      ::      this
            integer(kind=SELECTION_RULE),intent(in)     ::      rule
            real(kind=real64),intent(in),optional       ::      rho_in
            integer,dimension(:,:),intent(in),optional  ::      hkl_in

            type(Gvectors)                      ::      gv2             !   a set of g-vectors with double the radius, used to compute all the crystal structure factors.
            integer                             ::      nG2
            integer                             ::      ii,mm
            real(kind=real64),dimension(3)      ::      kk
            real(kind=real64)                               ::      dphi
            integer,dimension(:,:),allocatable              ::      hkl                 !   reflections
            complex(kind=real64),dimension(:),allocatable   ::      Vg                  !   crystal structure factor
            real(kind=real64),dimension(3,3)    ::      RR!,R0

            !R0 = RotationMatrix_Identity

            if (associated(this%gv)) then
                call delete(this%gv)
            end if

            if (present(rho_in)) then
                call setGvectors( this,rule,rho_in=rho_in )  
            else if (present(hkl_in)) then
                call setGvectors( this,rule,hkl_in=hkl_in )  
            else
                call errorExit("changeGvectors1 error - need to call with either radius |hkl|<rho_in or specified vector [hkl_in]")
            end if
            
            !call report(this%gv)
            
 
            gv2 = Gvectors_ctor(this%gv,doubleSet=.true.)
            nG2 = getn(gv2)
            if (rank==0) print *,"Lib_IntegrateManyBeams::changeGvectors1 info - computing ",nG2," crystal structure factors"


            if (associated(this%xi)) deallocate(this%xi)
            allocate(this%xi(nG2))
            if (isMillerBravais(gv2)) then
                allocate(hkl(4,nG2))
                do ii = 1,nG2
                    hkl(:,ii) = gethjkl(gv2,ii)
                end do
            else
                allocate(hkl(3,nG2))
                do ii = 1,nG2
                    hkl(:,ii) = gethkl(gv2,ii)
                end do
            end if
            allocate(Vg(nG2))
            if (rank == 0) then
                call computeCrystalStructureFactors( this%element,hkl, this%T,this%V, Vg )       !    note computeCrystalStructureFactors wants accelerator voltage in V not kV.
            end if
#ifdef MPI
            call MPI_BCAST( Vg,nG2,MPI_DOUBLE_COMPLEX,0,MPI_COMM_WORLD,ii)
#endif        
            call setCStructFactSet(this,.true.)
            
            do ii = 1,nG2                
                this%xi(ii) = PI * HBAR * velocity( this%V ) / Vg(ii)
                if (.not. this%lossy) this%xi(ii) = real(this%xi(ii))
            end do
            if (rank==0) then
                if (.not. this%lossy) then
                    if (isMillerBravais(this%gv)) then
                        write(*,fmt='(a16,a38,a16)') "reflection","V_g","xi_g"
                    else
                        write(*,fmt='(a12,a38,a16)') "reflection","V_g","xi_g"
                    end if
                else
                    if (isMillerBravais(this%gv)) then
                        write(*,fmt='(a16,a39,a39)') "reflection","V_g","xi_g"
                    else
                        write(*,fmt='(a12,a39,a39)') "reflection","V_g","xi_g"
                    end if
                end if
                if ( LIB_IMB_DBG .or. (nG2<=10) ) then
                    do ii = 1,nG2  
                        if (.not. this%lossy) then
                            if (isMillerBravais(this%gv)) then
                                write(*,fmt='(4i4,5(a,f16.8))') hkl(:,ii)," (",real(Vg(ii))," +i",aimag(Vg(ii))," ) ",real(this%xi(ii))
                            else
                                write(*,fmt='(3i4,5(a,f16.8))') hkl(:,ii)," (",real(Vg(ii))," +i",aimag(Vg(ii))," ) ",real(this%xi(ii))
                            end if
                        else
                            if (isMillerBravais(this%gv)) then
                                write(*,fmt='(4i4,5(a,f16.8))') hkl(:,ii)," (",real(Vg(ii))," +i",aimag(Vg(ii))," )  (",real(this%xi(ii))," +i",aimag(this%xi(ii)),")"
                            else
                                write(*,fmt='(3i4,5(a,f16.8))') hkl(:,ii)," (",real(Vg(ii))," +i",aimag(Vg(ii))," )  (",real(this%xi(ii))," +i",aimag(this%xi(ii)),")"
                            end if
                        end if
                    end do
                else
                    do ii = 1,nG2  
                        if (.not. this%lossy) then
                            if (isMillerBravais(this%gv)) then
                                if (norm2( real(hkl(:,ii)) ) <= 2.0)    &
                                write(*,fmt='(4i4,5(a,f16.8))') hkl(:,ii)," (",real(Vg(ii))," +i",aimag(Vg(ii))," ) ",real(this%xi(ii))
                            else
                                if (norm2( real(hkl(:,ii)) ) <= 2.0)    &
                                write(*,fmt='(3i4,5(a,f16.8))') hkl(:,ii)," (",real(Vg(ii))," +i",aimag(Vg(ii))," ) ",real(this%xi(ii))
                            end if
                        else
                            if (isMillerBravais(this%gv)) then
                                if (norm2( real(hkl(:,ii)) ) <= 2.0)    &
                                write(*,fmt='(4i4,5(a,f16.8))') hkl(:,ii)," (",real(Vg(ii))," +i",aimag(Vg(ii))," )  (",real(this%xi(ii))," +i",aimag(this%xi(ii)),")"
                            else
                                if (norm2( real(hkl(:,ii)) ) <= 2.0)    &
                                write(*,fmt='(3i4,5(a,f16.8))') hkl(:,ii)," (",real(Vg(ii))," +i",aimag(Vg(ii))," )  (",real(this%xi(ii))," +i",aimag(this%xi(ii)),")"
                            end if
                        end if
                    end do
                end if
                print *,""
            end if


            if (rank==0) print *,"Lib_IntegrateManyBeams::changeGvectors1 info - initialising many-beam calculations"

           
            if (associated(this%mb)) then
                do mm = 1,this%nPrec
                    call delete(this%mb(mm))
                end do
            else
                allocate(this%mb(this%nPrec))
            end if
            kk = (/0,0,1/) * 2 * PI / wavelength( this%V ) 
            this%mb(1) = ManyBeam_ctor( geta(this),kk,this%gv,gv2,this%xi)
            if (this%nPrec > 2) then
                !   precession
                dphi = 2*PI / (this%nPrec-1)
                do mm = 2,this%nPrec
                    RR = RotationMatrix_ctor( (/0.0d0,1.0d0,0.0d0/),this%precAngle )             !   rotation about y axis by precession angle 
                    RR = matmul( RotationMatrix_ctor( (/0.0d0,0.0d0,1.0d0/),(mm-2)*dphi ),RR )   !   rotation about z axis by angle phi 
                    this%mb(mm) = ManyBeam_ctor( geta(this),matmul(RR,kk),this%gv,gv2,this%xi )
                end do
            end if
            call setManyBeamSet(this,.true.)
            
            return
        end subroutine changeGvectors1

        



        subroutine perfectLatticeIntensity0(this,hkl,maxIntensity, I_g)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      compute the perfect lattice intensity estimate at reflection hkl
    !*      find either the max intensity in the foil, or the intensity at the foil exit
            type(IntegrateManyBeams),intent(in)         ::      this
            integer,dimension(:),intent(in)             ::      hkl    
            logical,intent(in)                          ::      maxIntensity
            real(kind=real64),intent(out)               ::      I_g

            real(kind=real64)           ::      LL
            integer                     ::      mm,nG,kk
            real(kind=real64),dimension(:),allocatable      ::      Ig,Igm
 
        !---    find the perfect lattice intensities
            LL = getThickness(this%as)
            if(rank==0) print *,"Lib_IntegrateManyBeams::perfectLatticeIntensity0 info - foil thickness ",LL
            nG = getn(this%gv)
            allocate(Ig(0:nG))
            allocate(Igm(0:nG))
            Ig = 0
            do mm = 1,this%nPrec
                call perfectLatticeIntensity(this%mb(mm), LL, maxIntensity ,Igm)
                Ig = Ig + Igm
            end do
            Ig = Ig / this%nPrec
             

            kk = whichg(this%gv,hkl)
            I_g = Ig(kk)

            return
        end subroutine perfectLatticeIntensity0
            


        subroutine applyFoilTilt(this,R)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      apply a foil tilt R
    !*      this will change the G-vectors and the AtomSpace
    !*      and needs to be propagated to the ManyBeams
    !*      note that x = exp[ -i g.u ]
    !*      so under a rotation R x' = exp[ -i (R g)T (R u) ] = x       ie x is invariant
            type(IntegrateManyBeams),intent(inout)          ::      this
            real(kind=real64),dimension(3,3),intent(in)     ::      R
            !real(kind=real64),dimension(3,3)     ::      RR
            integer             ::      mm
            if (rank==0) then
                print *,"Lib_IntegrateManyBeams::applyFoilTilt info - rotation matrix"
                write(*,fmt='(3f16.8)') R(1,:)
                write(*,fmt='(3f16.8)') R(2,:)
                write(*,fmt='(3f16.8)') R(3,:)
            end if
            !call report(this%gv)
            call rotate(this%gv,R)                                  !   update g-vectors - note that A -> R A, so reciprocal lattice B -> R B. 
            !call report(this%gv)
            call setR(this%as,R)                                    !   update atom space

            this%dfg_bar = matmul( R,this%dfg_bar )                         !   rotate crystal orientation

            do mm = 1,this%nPrec
                call setOrientationDependence(this%mb(mm))          !   mb has a pointer to gv
            end do
            call setImageSpaceSet(this,.false.)
            call setSimResultSet(this,.false.)
            return
        end subroutine applyFoilTilt
            

        subroutine getIntensity0(this,hkl,I)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      get the result at a given g-vector.
    !*      gather on rank 0
            type(IntegrateManyBeams),intent(in)     ::      this
            integer,dimension(:),intent(in)         ::      hkl
            real(kind=real64),dimension(:,:),allocatable    ::      I

            integer                     ::      ii
            integer                     ::      ix,iy,nx,ny
            integer                     ::      jx,jy,mx,my
            integer                     ::      pp,lbx,ubx,lby,uby
            
            integer                     ::      mm
            real(kind=real64)           ::      Ig
            complex(kind=real64)        ::      phig

#ifdef MPI
            real(kind=real64),dimension(:,:),allocatable    ::      I_tmp
            integer                     ::      ierror
#endif

        !---    which g-vector is this?
            ii = whichg(this%gv,hkl)
            if (ii==0) then
                if (rank==0) print *,"Lib_IntegrateManyBeams::getIntensity0 error - looking for reflection ",hkl
                call errorExit("getIntensity0 error - did not recognise g-vector as one of set")
            end if

        !---    allocate space for output image
            nx = getNx(this%is)
            ny = getNy(this%is)
            allocate(I(0:nx-1,0:ny-1))
            I = 0
            
            
        
            mx = getMx(this%is)
            my = getMy(this%is)
            !print *,"Lib_IntegrateManyBeams::getIntensity0 info - rank ",rank," nx,ny,mx,my,nPrec ",nx,ny,mx,my,this%nPrec
            do jy = 0,my-1
                do jx = 0,mx-1

                !---    find who has the blocks of this image
                    pp = whoseBlock( this%is,jx,jy )
                    call blockExtent( this%is,pp,lbx,ubx,lby,uby )

                   ! print *,"Lib_IntegrateManyBeams::getIntensity0 info - rank ",pp," computing ",lbx,":",ubx,",",lby,":",uby
                    if (pp == rank) then
                !---    compute intensity in block
                        !print *,"Lib_IntegrateManyBeams::getIntensity0 info - rank ",pp," computing ",lbx,":",ubx,",",lby,":",uby
                        do iy = lby,uby
                            do ix = lbx,ubx
                                Ig = 0.0d0
                                do mm = 1,this%nPrec
                                    phig = this%slice(mm)%phi( ii,ix,iy )
                                    Ig = Ig + real(phig)*real(phig) + aimag(phig)*aimag(phig)
                                end do
                                I(ix,iy) = Ig
                            end do
                        end do
                        
                    end if
 
                end do
            end do

#ifdef MPI
            if (nprocs>1) then
                allocate(I_tmp(0:nx-1,0:ny-1))
                I_tmp = I
                call MPI_REDUCE( I_tmp,I,nx*ny,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,ierror )
            end if
#endif        

            I = I / this%nPrec
           
            return
        end subroutine getIntensity0

        subroutine getIntensity1(this,I_g,maxIntensity)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      get the result at all g-vectors at the foil exit
    !*      gather on rank 0. Equivalent to perfectLatticeIntensity, but this will also produce solution if it is found
            type(IntegrateManyBeams),intent(in)             ::      this
            real(kind=real64),dimension(0:),intent(out)     ::      I_g
            logical,intent(in)                              ::      maxIntensity                !   if true and result not set, find maximum intensity inside foil, otherwise compute intensity at foil exit

            integer                     ::      ii,nG
            integer                     ::      ix,iy,nx,ny
            integer                     ::      jx,jy,mx,my
            integer                     ::      pp,lbx,ubx,lby,uby
          
            integer                     ::      mm
            real(kind=real64)           ::      phig2
            complex(kind=real64)        ::      phig

#ifdef MPI
            real(kind=real64),dimension(:),allocatable    ::      I_tmp
            integer                     ::      ierror
#endif

            if (.not. isGvectorsSet(this).and.isManyBeamSet(this)) then
                I_g = 0
                return
            end if

            nG = getN(this%gV)
            nx = getNx(this%is)
            ny = getNy(this%is)
            mx = getMx(this%is)
            my = getMy(this%is)

            if (.not. isSimResultSet(this)) then
                call perfectLatticeIntensity( this%mb(1), getL( this%as ) , maxIntensity , I_g )            !   note getL() returns true foil thickness at current rotation angle
                return
            end if


            !print *,"Lib_IntegrateManyBeams::getIntensity1 info - rank ",rank," nx,ny,mx,my,nPrec ",nx,ny,mx,my,this%nPrec
            I_g = 0
            do jy = 0,my-1
                do jx = 0,mx-1

                !---    find who has the blocks of this image
                    pp = whoseBlock( this%is,jx,jy )
                    call blockExtent( this%is,pp,lbx,ubx,lby,uby )

                   ! print *,"Lib_IntegrateManyBeams::getIntensity0 info - rank ",pp," computing ",lbx,":",ubx,",",lby,":",uby
                    if (pp == rank) then
                !---    compute intensity in block
                        !print *,"Lib_IntegrateManyBeams::getIntensity0 info - rank ",pp," computing ",lbx,":",ubx,",",lby,":",uby
                        do iy = lby,uby
                            do ix = lbx,ubx
                                do ii = 0,nG
                                    phig2 = 0.0d0
                                    do mm = 1,this%nPrec
                                        phig = this%slice(mm)%phi( ii,ix,iy )
                                        phig2 = phig2 + real(phig)*real(phig) + aimag(phig)*aimag(phig)
                                    end do
                                    I_g(ii) = I_g(ii) + phig2
                                end do
                            end do
                        end do
                        
                    end if
 
                end do
            end do

#ifdef MPI
            if (nprocs>1) then
                allocate(I_tmp(0:nG))
                I_tmp = I_g
                call MPI_REDUCE( I_tmp,I_g,size(I_tmp),MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,ierror )
            end if
#endif        

            I_g = I_g / (nx*ny*this%nPrec)

            return
        end subroutine getIntensity1

        subroutine selectFoilTilt( this,hkl,theta_max,nTrial,bright, R )
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      select a good foil orientation so that reflection hkl is bright/dark
    !*      
            type(IntegrateManyBeams),intent(inout)  ::      this                !   needs to be inout because this%mb(1) is changed ( and changed back )
            integer,dimension(:),intent(in)         ::      hkl                 !   reflection
            real(kind=real64),intent(in)            ::      theta_max           !   max tilt angle (radians)
            integer,intent(in)                      ::      nTrial              !   number of tilts to try
            logical,intent(in)                      ::      bright              !   if true, make [hkl] intensity high. Otherwise make low.
            real(kind=real64),dimension(3,3),intent(out)    ::      R


            real(kind=real64),dimension(3,3,nTrial)     ::      RR
            real(kind=real64)                           ::      L_R             !   foil thickness after rotation
            integer                                     ::      ii,nG,kk,bestRank
            real(kind=real64),dimension(:),allocatable  ::      Ig_R            !   intensity
            real(kind=real64)                           ::      qq,Ik_R 

            integer,dimension(:),allocatable            ::      besti,besti_tmp
            real(kind=real64),dimension(:),allocatable  ::      bestq,bestq_tmp
            logical             ::      maxIntensity = .false.                  !   compute foil tilt based on exit intensity
        


        !---    which reflection are we looking at? Do we have it in our set?
            nG = getn(this%gv)
            allocate(Ig_R(0:nG))
            kk = whichg(this%gv,hkl)
            if (kk == LIB_GVECTORS_UNSET) then
                if (rank==0) print *,"Lib_IntegrateManyBeams::selectFoilTilt error - looking for reflection ",hkl
                call errorExit("selectFoilTilt error - could not find desired reflection in set")
            end if

            !print *,"selectFoilTilt ",hkl,bright," g-vector ",kk
        !---    select the orientation of the foil tilts
            call fibonacci_cap(RR(:,3,:),theta_max)            
            

            do ii = 1,nTrial
                if (mod(ii,nProcs)/=rank) cycle
                RR(:,1,ii) = (/1,0,0/)
                call completeBasis( RR(:,3,ii),RR(:,1,ii),RR(:,2,ii) )
            end do


        !---    find the intensities of the reflections
            allocate(besti(0:nProcs-1))
            besti = 0
            allocate(bestq(0:nProcs-1))
            besti = 0
            bestq = 0 ; bestq(rank) = -huge(1.0)
            do ii = 1,nTrial
                if (mod(ii,nProcs)/=rank) cycle
                L_R = getThickness( this%as,RR(:,:,ii) )
                call perfectLatticeIntensity( this%mb(1), L_R, RR(:,:,ii), maxIntensity , Ig_R )         !   compute fo
                Ik_R = Ig_R(kk)
                if (bright) then
                    Ig_R(kk) = 0.0d0                
                    qq = Ik_R - maxval(Ig_R(1:nG))              !   high if selected reflection is bright                    
                else
                    qq = 1 - Ik_R                               !   high if selected reflection is dark
                end if
                if (qq > bestq(rank)) then
                    bestq(rank) = qq
                    besti(rank) = ii                    
                end if
                !write(*,fmt='(a,i6,100f12.6)') "selectFoilTilt ",ii,qq,RR(:,:,ii),Ik_R
            end do

#ifdef MPI
            if (nProcs>1) then
                allocate(besti_tmp(0:nProcs-1))
                besti_tmp = besti
                call MPI_ALLREDUCE(besti_tmp,besti,nProcs,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,ii )
                deallocate(besti_tmp)
                allocate(bestq_tmp(0:nProcs-1))
                bestq_tmp = bestq
                call MPI_ALLREDUCE(bestq_tmp,bestq,nProcs,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ii )
                deallocate(bestq_tmp)
            end if
#endif

            !if (rank==0) print *,"bestq ",bestq(rank)
            bestRank = maxloc(bestq,dim=1)-1      !   which rank found the best quality?
            qq = bestq(bestRank)
            !if (rank == 0) print *,"besti ",qq,bestRank,besti(bestRank)
            ii = besti(bestRank)                  !   which trial was it?
            if (bestRank==rank) then
                R(:,:) = RR(:,:,ii)
            else
                R(:,3) = RR(:,3,ii)
                R(:,1) = (/1,0,0/)
                call completeBasis( R(:,3),R(:,1),R(:,2) )
            end if
            L_R = getThickness( this%as,R )                 !   returns new thickness, assuming that a virtual tilt is applied
            if (rank==0) print *,R,L_R
            call perfectLatticeIntensity( this%mb(1), L_R, R, maxIntensity, Ig_R )
            Ik_R = Ig_R(kk)
            if (bright) then                
                Ig_R(kk) = 0.0d0                
                qq = Ik_R - maxval(Ig_R(1:nG))              !   high if selected reflection is bright                    
            else
                qq = 1 - Ik_R                               !   high if selected reflection is dark
            end if

        !---    return best solution
            if (rank==0) then
                write(*,fmt='(a,9f10.6,2(a,f10.6))') " foil tilt "
                write(*,fmt='(3f20.12)')  R(1,:)
                write(*,fmt='(3f20.12)')  R(2,:)
                write(*,fmt='(3f20.12)')  R(3,:)
                write(*,fmt='(2(a,f10.6))') " quality ",qq," I_g ",Ik_R
            end if
            !R = RR(:,:,ii)


            return
        end subroutine selectFoilTilt













        function getImagingSpace(this) result(is)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
            type(IntegrateManyBeams),intent(in)     ::      this
            type(imagingSpace)                      ::      is
            is = this%is
            return
        end function getImagingSpace



        function getGvectors(this) result(gv)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
            type(IntegrateManyBeams),intent(in)     ::      this
            type(Gvectors)                          ::      gv
            gv = this%gv
            return
        end function getGvectors


        subroutine integrate0(this)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      propagate many beam equations in z through the system
            type(IntegrateManyBeams),intent(inout)     ::      this

            type(RK4)           ::      integrator
            complex(kind=real64),dimension(:),pointer       ::      rk_phip
            complex(kind=real64),dimension(:),pointer       ::      rk_dphip
            complex(kind=real64),dimension(:,:,:),pointer   ::      dphidz
            
            integer,dimension(2,3)                          ::      bb
            integer             ::      iz,rk_step,iz_layer
            integer             ::      mm
            integer             ::      nx,ny,nz,ng
            integer             ::      my_nx,my_ny
            integer             ::      dz
            real(kind=real64)   ::      deltaz

        !---    sanity test
            if (.not. isFileRead(this)) call errorExit("integrate0 error - file read not complete")
            if (.not. isAtomSpaceSet(this)) call errorExit("integrate0 error - AtomSpace not set")
            if (.not. isGvectorsSet(this)) call errorExit("integrate0 error - g-vectors not set")
            if (.not. isManyBeamSet(this)) call errorExit("integrate0 error - ManyBeam conditions not set")
            if (.not. isCStructFactSet(this)) call errorExit("integrate0 error - structure factors not set")
            if (.not. isImageSpaceSet(this)) call errorExit("integrate0 error - image space not set")


            
            !call computePhaseFields(this,0) 


        !---    allocate memory and temporaries needed for RK4 integration
            nx = getNx(this%is)
            ny = getNy(this%is)
            nz = getNz(this%is)
            nG = getn(this%gv)
            call getBounds(this%is,bb)
            bb(1,:) = bb(1,:) - 1           !   buffer pixel
            bb(2,:) = bb(2,:) + 1           !   buffer pixel
            my_nx = bb(2,1)+1-bb(1,1)
            my_ny = bb(2,2)+1-bb(1,2)
            
            deltaz = 2*geta(this%is)                            !   dz = 2a because the RK4 integrator needs to evaluate at half steps
            integrator = RK4_ctor( n=(nG+1) * my_nx * my_ny, deltaz=deltaz  )     
            if (rank==0) call report(integrator)
            if (rank==0) print *,"Lib_IntegrateManyBeams::integrate0 info - bounds: nx,ny,nz = ",nx,ny,nz," lbx,ubx,lby,uby = ",bb(1:2,1),bb(1:2,2)," my_nx,my_ny = ",my_nx,my_ny

        !---    get pointers to current state and derivative stored in the integrator
            call getphip_dphip(integrator , rk_phip,rk_dphip)           
            allocate(dphidz(0:nG,bb(1,1):bb(2,1),bb(1,2):bb(2,2)))


        !---    set boundary conditions at z=0
            dphidz = 0      
            rk_phip = 0
            rk_dphip = 0
            do mm = 1,this%nPrec
                this%slice(mm)%phi(0   ,:,:) = 1.0d0
                this%slice(mm)%phi(1:nG,:,:) = 0.0d0
            end do

            
        !---    integrate!   
            do iz = 0,getNz(this%is)-3,2
                !if (rank == 0) call progressBar(iz/2+1,(getNz(this%is)-1)/2)

                if (mod(iz,this%zlayers)==0) then
                    iz_layer = iz
                    call computePhaseFields(this,iz_layer) 
                end if


                do mm = 1,this%nPrec

                    rk_phip = pack(this%slice(mm)%phi,.true.)
                    call update(integrator)
                    dz = 0
                    do rk_step = 1,4
                    !   ensure that phi is distributed correctly 
                        this%slice(mm)%phi(0:nG,bb(1,1):bb(2,1),bb(1,2):bb(2,2)) = reshape( rk_phip, (/nG+1,my_nx,my_ny/) )
                    !   compute dphidz
                        call finddPhidz( this%mb(mm) , this%slice(mm)%phi , iz-iz_layer+dz, dphidz, this%columnar )
                        rk_dphip = pack(dphidz,.true.)
                    !   update accumulators in RK4
                        call update(integrator,rk_step,dz)
                       
                    end do
                    this%slice(mm)%phi(0:nG,bb(1,1):bb(2,1),bb(1,2):bb(2,2)) = reshape( rk_phip, (/nG+1,my_nx,my_ny/) )
                    if (.not. this%columnar) call sendrecv( this%is,this%slice(mm)%phi )

                end do

            end do
            call setSimResultSet(this,.true.)

        

            call delete(integrator)
            return
        end subroutine integrate0

!-------

        subroutine setGvectors( this,rule,rho_in,hkl_in )
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      given the name of the lattice, and a cutoff radius rho for 
    !*      the reflections, construct a set of g-vectors.

            type(IntegrateManyBeams),intent(inout)          ::      this
            integer(kind=SELECTION_RULE),intent(in)         ::      rule
            real(kind=real64),intent(in),optional           ::      rho_in          !   compute all g-vectors within given radius
            integer,dimension(:,:),intent(in),optional      ::      hkl_in          !   compute sepcified reflections

            real(kind=real64),dimension(3,3)        ::      a_conventional_cell
            integer                                 ::      nn
            integer,dimension(:,:),allocatable      ::      hkl,hjkl
            logical,dimension(:),allocatable        ::      keep
            integer,dimension(3)                    ::      hkl_swp
            real(kind=real64),dimension(3)          ::      kk,gg
            integer                                 ::      ii,jj
            logical                 ::      done



            if (rank==0) then
                print *,"Lib_IntegrateManyBeams::setGvectors info - previous g-vectors ",isGvectorsSet(this)
                print *,"Lib_IntegrateManyBeams::setGvectors info - deformation gradient"
                write (*,fmt='(3f16.8)') this%dfg_bar(1,:)
                write (*,fmt='(3f16.8)') this%dfg_bar(2,:)
                write (*,fmt='(3f16.8)') this%dfg_bar(3,:)
                if (present(rho_in)) then
                    print *,"Lib_IntegrateManyBeams::setGvectors info - many beam mode, radius |[hkl]|<=",rho_in
                else if (present(hkl_in)) then
                    print *,"Lib_IntegrateManyBeams::setGvectors info - selecting [000] + additional ",size(hkl_in,dim=2)," beams"
                end if
            end if

            



            if (present(rho_in)) then
            !---    find all reflections within certain radius
                call permittedReflections( this%latt,nn,hkl,rho_in )            
            else if (present(hkl_in)) then
                nn = size(hkl_in,dim=2) + 1
                allocate(hkl(3,nn))
                hkl(:,1) = 0

                if (size(hkl_in,dim=1)==3) then
                    hkl(:,2:) = hkl_in(:,1:)
                else
                    do ii = 1,size(hkl_in,dim=2)
                        hkl(:,ii+1) = MillerBravaisToMiller_plane( hkl_in(:,ii) )
                    end do
                end if
            end if

        !---    find the transmitted beam vector
            kk = (/0,0,1/) * ME * velocity( this%V )/ HBAR                 


        !---    compute the correct conventional cell according to hints from the input file: lattice type, deformation gradient, and lattice constant.
        !       note that this respects the strain in the input file, but does not apply a foil tilt
            a_conventional_cell = getConventionalCell(this%latt)
            a_conventional_cell = matmul( this%dfg_bar,a_conventional_cell )
            a_conventional_cell = this%a0 * a_conventional_cell            


        !---    check which reflections we will keep 

            allocate(keep(nn))
            this%gv = Gvectors_ctor( a_conventional_cell,hkl(:,2:nn) )        !   make a temporary set of g-vectors, so that we can test for Laue zones
            do ii = 1,nn                
                gg = getg( a_conventional_cell,hkl(:,ii) )
                keep(ii) = keepGvector( this,gg,rule )
            end do
            jj = 0
            do ii = 1,nn
                if (keep(ii)) then
                    jj = jj + 1
                    hkl(:,jj) = hkl(:,ii)
                !else
                !    print *,"rejecting ",hkl(:,ii),deviationParameter( (/0,0,1/) * ME * velocity( this%V ) / HBAR , getg( a_conventional_cell,hkl(:,ii) ) )
                end if
                
            end do
            call delete(this%gv)                                            !   delete the temporary set of g-vectors.
            nn = count(keep)
            if (rank==0) print *,"Lib_IntegrateManyBeams::setGvectors info -keeping ",nn,"/",size(keep)
            
 
        !---    bubble sort permitted reflections
            do
                done = .true.
                do ii = 1,nn-1
                    if ( dot_product(hkl(:,ii),hkl(:,ii)) > dot_product(hkl(:,ii+1),hkl(:,ii+1)) ) then
                        hkl_swp(:)  = hkl(:,ii)
                        hkl(:,ii)   = hkl(:,ii+1)
                        hkl(:,ii+1) = hkl_swp(:)
                        done = .false.
                    end if
                end do
                if (done) exit
            end do
        !---    now permitted reflection 1 should be [000]

         
    !---    construct a g-vector object            
            if (getLatticeType(this%latt) == LATTICE_HCP) then
                allocate(hjkl(4,nn))
                do ii = 1,nn
                    hjkl(:,ii) = nint( MillerToMillerBravais_plane(real(hkl(:,ii),kind=real64)) )
                end do
                this%gv = Gvectors_ctor(a_conventional_cell,hjkl(:,2:nn))
            else
                this%gv = Gvectors_ctor(a_conventional_cell,hkl(:,2:nn))
            end if
 
            if (rank==0) call report(this%gv)
            if (rank==0) print *,""
            call setGvectorsSet(this,.true.)
            return
        end subroutine setGvectors


 
        logical function keepGvector( this,gg,rule )
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      should this g-vector be excluded by my self-imposed rules?
            type(IntegrateManyBeams),intent(in)         ::      this
            real(kind=real64),dimension(3),intent(in)   ::      gg
            integer(kind=SELECTION_RULE),intent(in)     ::      rule
            real(kind=real64),dimension(3)          ::      kk
            real(kind=real64)                       ::      sg
            integer                                 ::      zz
            
            logical,save                            ::      alpha_found = .false.
            real(kind=real64),save                  ::      alpha = 0.0d0
            integer                                 ::      ii
            integer                                 ::      closest_g_to_001
            real(kind=real64)                       ::      dd,dmax,modk
            
            

            select case(rule)

                case (LIB_IMB_BLOCK_HIGH_SG)

                    kk = (/0,0,1/) * ME * velocity( this%V ) / HBAR     
                    sg = deviationParameter( kk,gg )                
                    keepGvector = (abs(2*PI*sg*geta(this)) < 0.5d0) 

                case (LIB_IMB_ZOLZ:LIB_IMB_SOLZ)
 
                    if (.not. alpha_found) then                          
                        
                    !---    find g-vector closest to [001]
                        dmax = 0.1d0                   !   don't accept any old rubbish, has to be within 10%
                        closest_g_to_001 = LIB_GVECTORS_UNSET
                        do ii = 1,getn(this%gv)
                            kk = getg(this%gv,ii)       !   find a g-vector ...
                            modk = norm2(kk)
                            if (modk>0) then
                                dd = kk(3) / modk           !   ... find its normalisd projection on [001] direction
                                if (dd > dmax) then         !   this is closest to z-direction
                                    dmax = dd
                                    closest_g_to_001 = ii
                                end if
                            end if
                        end do
                        if (closest_g_to_001 /= LIB_GVECTORS_UNSET) then
                            alpha = 1/(dmax)              !   normalising factor so that alpha kk.[001] = 1
                        end if
                        alpha_found = .true.
                        
                    end if                 

                    zz = nint( alpha * gg(3) )              !   zone order
                    !print *,"vec ",gg," Laue zone order ",zz

                    if (rule == LIB_IMB_ZOLZ) then
                        keepGvector = (zz==0)
                    else if (rule == LIB_IMB_FOLZ) then
                        keepGvector = (abs(zz)<=1)
                    else 
                        keepGvector = (abs(zz)<=2)
                    end if

                case default

                    keepGvector = .true.

            end select
            return
        end function keepGvector
 


        subroutine readInputXyzFile( this,filename,xyz_offset,a_super )
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      reads in the atom positions ( shared by all processes )
    !*      computes the average deformation gradient (this%dfg_bar) if possible
    !*      sets the most populous element name (this%element)
    !*      returns the offset of atoms in file, and periodic supercell lattice vectors
     
            type(IntegrateManyBeams),intent(inout)              ::      this
            character(len=*),intent(in)                         ::      filename
            real(kind=real64),dimension(3),intent(out)          ::      xyz_offset  !   atom position offset read from file
            real(kind=real64),dimension(3,3),intent(out)        ::      a_super

            type(XYZFile)       ::      xyz
            logical             ::      ok,ok_tmp
            
            integer             ::      ii,kk,nGrains
            integer,dimension(:),allocatable        ::  nGrain
            real(kind=real64),dimension(3,3)    ::      eps_bar,rot_bar,dfg
            real(kind=real64)                   ::      weight
            integer,dimension(:),pointer        ::      tt
            integer                             ::      mostPopType,nTypeMax
             

            if (rank==0) print *,"Lib_IntegrateManyBeams::readInputXyzFile info"

 

            if (this%nAtoms > 0) deallocate(this%r)
            inquire(file=trim(filename),exist=ok)
            if (.not. ok) call errorExit("Lib_IntegrateManyBeams::readInputXyzFile ERROR - file not found """//trim(filename)//"""")


        !---    read the file into rank 0            
            if (rank == 0) then    
                xyz = XYZFile_ctor(filename)
                call readHeader(xyz,ok)
 

                if (ok) then
                    call input(xyz,verbose=.true.)
                    call report(xyz)        
                    this%nAtoms = getNAtoms(xyz) 
                    call getColumnsp(xyz,this%r)
                    this%r = this%r + 1.0e-8
                    call getSupercell(xyz,a_super,ok) 
                    if (.not. ok) print *,"Lib_IntegrateManyBeams::readInputXyzFile error - could not read supercell information"
                    if (getNColumns(xyz)>=13) then
                        print *,"Lib_IntegrateManyBeams::readInputXyzFile info - assuming column 4 is grain number, 5:13 are deformation gradient"
                        nGrains = maxval(nint( this%r(4,:) ))
                        print *,"Lib_IntegrateManyBeams::readInputXyzFile info - read number of grains from .xyz file = ",nGrains
                        allocate(nGrain(0:nGrains))
                        nGrain = 0
                        do ii = 1,this%nAtoms
                            kk = nint( this%r(4,ii) )  
                            nGrain( kk ) = nGrain( kk ) + 1
                        end do
                        
                        print *,"Lib_IntegrateManyBeams::readInputXyzFile info - atom count in each grain"
                        do kk = 0,nGrains
                            print *,"    grain ",kk," count ",nGrain(kk)
                        end do
        
                        kk = maxloc(nGrain,dim=1) - 1
                        print *,"Lib_IntegrateManyBeams::readInputXyzFile info - using average deformation gradient for most populous grain ",kk," frac ", nGrain(kk)*100.0/this%nAtoms," %"
                        deallocate(nGrain)
                                    
                        call computeAvgDefGrad(eps_bar,this%dfg_bar,weight)
                        do ii = 1,this%nAtoms
                            if (  nint( this%r(4,ii) ) == kk) call computeAvgDefGrad( this%r(5:13,ii),1.0d0, eps_bar,this%dfg_bar,weight)
                        end do
                        call computeAvgDefGrad(eps_bar,this%dfg_bar,weight , dfg) ; this%dfg_bar = dfg
                        print *,"d2bi info - average deformation gradient, strain, rot read from extended .xyz file "

                        call DefGradToStrainAndRotMat(this%dfg_bar,eps_bar,rot_bar)
        
                    else        !   don't have columns data
                        if (rank==0) then
                            print *,"Lib_IntegrateManyBeams::readInputXyzFile info - don't have deformation gradient data in input file"
                            print *,"   assuming lattice orientation aligned with xyz coordinates"
                        end if
                        this%dfg_bar = reshape( (/1,0,0,0,1,0,0,0,1/),(/3,3/) )
                        eps_bar = 0
                        rot_bar = reshape( (/1,0,0,0,1,0,0,0,1/),(/3,3/) )
                    
                    end if
                    if (rank==0) then
                        write (*,fmt='(a36,a,a36,a,a36)') "average def grad","    ","average strain","    ","average rotation"
                        write (*,fmt='(3f12.6,a,3f12.6,a,3f12.6)') this%dfg_bar(1,:),"    ",eps_bar(1,:),"    ",rot_bar(1,:)
                        write (*,fmt='(3f12.6,a,3f12.6,a,3f12.6)') this%dfg_bar(2,:),"    ",eps_bar(2,:),"    ",rot_bar(2,:)
                        write (*,fmt='(3f12.6,a,3f12.6,a,3f12.6)') this%dfg_bar(3,:),"    ",eps_bar(3,:),"    ",rot_bar(3,:)
                    end if   

                !---    attempt to id the most populous element
                    call getTypesp(xyz,tt)
                    nTypeMax = -huge(1) 
                    mostPopType = 1
                    do ii = 1,getnAtomNames(xyz)
                        kk = count(tt == ii)
                        if (kk>nTypeMax) then
                            mostPopType = ii
                            nTypeMax = kk
                        end if
                    end do
                    this%element = getAtomName(xyz,mostPopType)
                    print *,"Lib_IntegrateManyBeams::readInputXyzFile info - most common atom type """//trim(this%element)//""""

                !---    attempt to read atom offset
                    call getOrigin(xyz,xyz_offset,ok_tmp)
                    if (ok_tmp) then
                        print *,"Lib_IntegrateManyBeams::readInputXyzFile info - read atom offset ",xyz_offset
                    else
                        print *,"Lib_IntegrateManyBeams::readInputXyzFile info - zero atom offset"
                    end if

                else
                    print *,"Lib_IntegrateManyBeams::readInputXyzFile error - could not read file header"
                end if  !   header read ok
                print *,""
                
            end if  !   rank 0
#ifdef MPI            
            call MPI_BCAST(ok,1,MPI_LOGICAL,0,MPI_COMM_WORLD,ii)
            if (.not. ok) call errorExit("readInputXyzFile ERROR - could not read file """//trim(filename)//"""")
            call MPI_BCAST(this%nAtoms,1,MPI_INTEGER,0,MPI_COMM_WORLD,ii)
            call MPI_BCAST(this%element,XYZFILE_ATOMNAMELENGTH,MPI_CHARACTER,0,MPI_COMM_WORLD,ii)
            if (rank/=0) allocate(this%r(3,this%nAtoms))
            call MPI_BCAST(this%r(1:3,1:this%nAtoms),3*this%nAtoms,MPI_REAL,0,MPI_COMM_WORLD,ii)                
            call MPI_BCAST(a_super,9,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ii)
            call MPI_BCAST(xyz_offset,3,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ii)
            call MPI_BCAST(this%dfg_bar,9,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ii)
#else
            if (.not. ok) call errorExit("readInputXyzFile ERROR - could not read file """//trim(filename)//"""")
            if (.not. ok) call errorExit("readInputXyzFile ERROR - could not read file """//trim(filename)//"""")
#endif  
            
            

            return
        end subroutine readInputXyzFile
            



        subroutine setImagingSpace(this)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ 
    !*      once the g-vectors and the atom space have been determined, can set the imaging space
    !*      We can also make transformation on the atom positions, and store all the periodic copies I need.
    !*      will also allocate memory for the solution.

            type(IntegrateManyBeams),intent(inout)          ::      this
        !    real(kind=real64),intent(in)                    ::      sigma
            integer             ::      Nx,Ny,Nz
            real(kind=real64)   ::      theta           !   max angle required for beam dispersion
            integer             ::      ii,jj,mm
            integer             ::      nG
            real(kind=real32),dimension(:,:),pointer        ::      rt_tmp
            real(kind=real64),dimension(3,9)                ::      xtp
            integer                                         ::      np
            integer                                         ::      mynAtoms
            integer,dimension(2,3)                          ::      bounds
           ! if (rank==0) print *,"Lib_IntegrateManyBeams::setImagingSpace"

            if ( isImageSpaceSet(this) ) call errorExit( "setImagingSpace error - duplicate call image space set already.")

        !---    find the maximum divergence angle for the beams used            
            theta = 0
            do mm = 1,this%nPrec
                do ii = 1,getn(this%gv)
                    theta = max(theta,getAngle(this%mb(mm),ii))
                   ! print *,"Lib_IntegrateManyBeams::setImagingSpace info - angle prec,g ",mm,ii," = ",getAngle(this%mb(mm),ii)
                end do
            end do
            
            call suggestImagingSpace(this%as,Nx,Ny,Nz,theta)
           ! if (rank==0) print *,"Lib_IntegrateManyBeams::setImagingSpace info - theta = ",theta," Nx,Ny,Nz = ",Nx,Ny,Nz
            call setDelta(this%as,Nx,Ny,Nz)             !   can't know the offset until the imaging space is known.

            this%is = ImagingSpace_ctor( geta(this),getSigma(this), getdelta(this%as),Nx,Ny,Nz )
            if (rank==0) call report(this%is) 
             

        !---    can now make my transformation on the atom positions, and store all the periodic copies I need.
        !       first count the number of atoms after periodic copies...
            mynAtoms = 0
            do ii = 1,this%nAtoms
                call periodicCopies(this%as,this%is,this%r(:,ii),np,xtp)
                do jj = 1,np
                    if (inMyCell(this%is,xtp(:,jj),buffered=.true.)) mynAtoms = mynAtoms + 1
                end do
            end do

        !       then store the atoms after periodic copies
            allocate(rt_tmp(3,mynAtoms))
            mynAtoms = 0
            do ii = 1,this%nAtoms
                call periodicCopies(this%as,this%is,this%r(:,ii),np,xtp)
                do jj = 1,np
                    if (inMyCell(this%is,xtp(:,jj),buffered=.true.)) then
                        mynAtoms = mynAtoms + 1
                        rt_tmp(1:3,mynAtoms) = real( xtp(1:3,jj),kind=real32 )
                    end if
                end do
            end do
            print *,"Lib_IntegrateManyBeams::setImagingSpace info - rank ",rank," mynAtoms = ",mynAtoms,"/",this%nAtoms
            print *,"Lib_IntegrateManyBeams::setImagingSpace info - rank ",rank," minmax x ",minval(rt_tmp(1,:)),maxval(rt_tmp(1,:))
            print *,"Lib_IntegrateManyBeams::setImagingSpace info - rank ",rank," minmax y ",minval(rt_tmp(2,:)),maxval(rt_tmp(2,:))
            print *,"Lib_IntegrateManyBeams::setImagingSpace info - rank ",rank," minmax z ",minval(rt_tmp(3,:)),maxval(rt_tmp(3,:))

            deallocate(this%r)                                  !   note that this deallocates the memory read in from the .xyz file
            this%nAtoms = mynAtoms
            this%r => rt_tmp
            if (rank==0) print *,""


        !---    allocate memory for solution
            nG = getn(this%gv)
            call getBounds(this%is,bounds)
            allocate(this%slice(this%nPrec))
            do mm = 1,this%nPrec                
                allocate(this%slice(mm)%phi(0:nG,bounds(1,1)-1:bounds(2,1)+1,bounds(1,2)-1:bounds(2,2)+1))
                this%slice(mm)%phi = 0
            end do            
            if (rank==0) print *,"Lib_IntegrateManyBeams::setImagingSpace info - memory for Phi_g(r) ", size(this%slice(1)%phi)*16.0*this%nPrec/(1024*1024)," (Mb)"

            call setImageSpaceSet(this,.true.)
            return
        end subroutine setImagingSpace


        subroutine computePhaseFields(this , iz)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      once the imaging space is set, I can compute the phase factors in my region
     
            type(IntegrateManyBeams),intent(inout)          ::      this
            integer,intent(in)                              ::      iz              !   compute and store rho, grad_arg_x for iz:iz+zlayers-1

            integer,dimension(2,3)                          ::      bb              !   bounds of imaging space
            integer                                         ::      nG              !   number of g-vectors in many beams, not equal to number of phase fields computed
            integer                                         ::      ii,mm ,jx,jy,jz
            integer                                         ::      iz_high
            real(kind=real64)                               ::      mingrad_arg_x,maxgrad_arg_x,sumgrad_arg_x,modgrad_arg_x
            real                                            ::      nvox

        !---    static variables which I don't need to compute every time
            logical,save                                        ::      firstCall = .true.
            integer,save                                        ::      nGcalc 
            real(kind=real64),dimension(:,:),allocatable,save   ::      gg          !   (3,nG)      g-vectors



            if (rank==0) print *,"Lib_IntegrateManyBeams::computePhaseFields info - zlayer         ",iz,":",iz+this%zlayers-1," (+2)"

        !---    size of problem                        
            call getBounds(this%is,bb)             
            nG = getn(this%gv)           
            nvox = real(bb(2,1)+1-bb(1,1))*real(bb(2,2)+1-bb(1,2))*(this%zlayers + 2)               !   note: need to compute zlayers + 2 for rk4 to work


        !---    allocate memory for a z-slice of rho            
            if (firstCall) then
                if (rank==0) then
                    print *,"Lib_IntegrateManyBeams::computePhaseFields info - z bounds       ",bb(1,3),":",bb(2,3)
                    !print *,"Lib_IntegrateManyBeams::computePhaseFields info - computing      ",iz,":",iz+this%zlayers-1
                    print *,"Lib_IntegrateManyBeams::computePhaseFields info - memory for rho ", nvox*4.0/(1024.0*1024.0)," Mb"
                    !print *,"Lib_IntegrateManyBeams::computePhaseFields info - bb = ",bb," nG = ",nG
                end if
                allocate(this%rho( bb(1,1):bb(2,1),bb(1,2):bb(2,2),0:this%zlayers+1 ))                    
            end if

        !---    we don't need to compute phase factors for all the g-vectors.
        !       g=0 gives a constant, x(g=0) = exp[ -i g.u ] = 1
        !       and x(-g) =  exp[ -i (-g).u ] = x(g)*
            if (firstCall) then
                nGcalc = nPositiveg(this%gv)
                allocate(gg(3,nGcalc))

                nGcalc = 0
                nG = getn(this%gv)
                do ii = 1,nG
                    if (isPositiveg(this%gv,ii)) then
                        nGcalc = nGcalc + 1
                        gg(:,nGcalc) = getG(this%gv,ii)
                        if ((rank==0).and.LIB_IMB_DBG) print *,"Lib_IntegrateManyBeams::computePhaseFields info - compute grad_arg_x ",nGcalc," g-vec ",ii,gg(:,nGcalc)
                    else
                        !   already have computed -g, so don't need to compute +g
                    end if                
                end do
            end if
            
        !---    allocate memory for a z-slice of grad_arg_x
            if (firstCall) then
                if (rank==0) then
                    print *,"Lib_IntegrateManyBeams::computePhaseFields info - memory for grad_arg_x ", (3*nGcalc)*nvox*8.0/(1024.0*1024.0)," Mb"
                end if
                allocate(this%grad_arg_x(3,nGcalc,bb(1,1):bb(2,1),bb(1,2):bb(2,2),0:this%zlayers+1))        !   note: 
            end if


        !---    compute the phase factor for my z-slice            
            iz_high = min( getnz(this%is) , iz + this%zlayers + 1 )
            call computePhaseFactor( this%nAtoms,this%r,gg, this%is, getDelta(this%as),iz, iz_high, this%grad_arg_x,this%rho )
            this%rho = real( this%rho/geta(this)**3 , kind=real32 )
            !if ((rank==0).and.LIB_IMB_DBG) print *,"Lib_IntegrateManyBeams::computePhaseFields info - minmaxavg rho      ",minval(abs(this%rho)),maxval(abs(this%rho)),sum(abs(this%rho))/(size(this%rho))," (1/A^3) "
            this%rho = densityScale( this%rho, real(this%omega0,kind=real32 ) )
            if (rank==0) print *,"Lib_IntegrateManyBeams::computePhaseFields info - minmaxavg rho'     ",minval(abs(this%rho)),maxval(abs(this%rho)),sum(abs(this%rho))/(size(this%rho))


        !---    dump of atom positions and z-slice
            if (LIB_IMB_DBG .and. (iz==0)) then                
                open(unit=500,file=trim(numberFile("test",rank,"xyz")),action="write")
                    write(unit=500,fmt='(i8)') this%nAtoms + (bb(2,1)+1-bb(1,1))*(bb(2,2)+1-bb(1,2))*this%zlayers
                    write(unit=500,fmt='(a,9f12.6,a)') "Lattice=""",getA_super(this%as),""" Properties=Species:S:1:Pos:R:3:rho:R:1:GradArgX:R:3"
                    do ii = 1,this%nAtoms
                        write(unit=500,fmt='(a,10f16.8)') this%element,fromImagingSpace( this%as,this%r(:,ii) ),0.0d0,0.0d0,0.0d0,0.0d0
                    end do
                    do jz = 0,this%zlayers-1
                        do jy = bb(1,2),bb(2,2)
                            do jx = bb(1,1),bb(2,1)
                                write(unit=500,fmt='(a,10f16.6)') "g",fromImagingSpace( this%as,(/jx,jy,iz+jz/) ),this%rho(jx,jy,jz),this%grad_arg_x(:,1,jx,jy,jz)
                            end do
                        end do
                    end do                    
                close(unit=500)
            end if


        !---    unpack the computed g-vectors

           if ((rank==0).and.LIB_IMB_DBG) then
                do ii = 1,nGcalc
                    mingrad_arg_x = huge(1.0)
                    maxgrad_arg_x = -huge(1.0)
                    sumgrad_arg_x = 0
                    do jz = lbound(this%grad_arg_x,dim=5),ubound(this%grad_arg_x,dim=5)
                        do jy = lbound(this%grad_arg_x,dim=4),ubound(this%grad_arg_x,dim=4)
                            do jx = lbound(this%grad_arg_x,dim=3),ubound(this%grad_arg_x,dim=3)
                                modgrad_arg_x = norm2( this%grad_arg_x(:,ii,jx,jy,jz ))
                                mingrad_arg_x = min( mingrad_arg_x,modgrad_arg_x )
                                maxgrad_arg_x = max( maxgrad_arg_x,modgrad_arg_x )
                                sumgrad_arg_x = sumgrad_arg_x + modgrad_arg_x
                                !if ( (ii==4).and.(modgrad_arg_x==0) ) print *,"zero? ",ii,jx,jy,jz,this%grad_arg_x(:,ii,jx,jy,jz ),this%rho(jx,jy,jz)
                            end do
                        end do
                    end do
                    
                    print *,"Lib_IntegrateManyBeams::computePhaseFields info - minmaxavg |grad_arg_x(",ii,")| ",mingrad_arg_x,maxgrad_arg_x,sumgrad_arg_x/nvox
                end do
            end if
            

        !---    set the phase factors for my z-slice           
            do mm = 1,this%nPrec
                call setPhaseFactors( this%mb(mm),iz,iz + this%zlayers - 1, this%grad_arg_x,this%rho )
            end do
            

            firstCall = .false.
            !call setPhaseFieldSet(this,.true.)
            return
        end subroutine computePhaseFields


        pure real(kind=real64) function getSigma0(this)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      returns width of blurring used to construct phase fields
            type(IntegrateManyBeams),intent(in)           ::      this            
            getSigma0 = this%a0/this%M
            return
        end function getSigma0

        pure real(kind=real64) function getA0(this)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      returns pixel size
            type(IntegrateManyBeams),intent(in)           ::      this            
            getA0 = this%a0/this%M
            return
        end function getA0
 

    !---    calculation status

        pure logical function sanityCheck0(this)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      returns true if all the pointers are set
            type(IntegrateManyBeams),intent(in)           ::      this            
            integer         ::      mm

            sanityCheck0 = (getn(this%gv)>0)

            sanityCheck0 = sanityCheck0 .and. associated(this%grad_arg_x)
            sanityCheck0 = sanityCheck0 .and. associated(this%rho)

            do mm = 1,this%nPrec
                sanityCheck0 = sanityCheck0 .and. sanityCheck(this%mb(mm))
            end do
            

            return
        end function sanityCheck0


        pure subroutine setFileRead(this,is)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
            type(IntegrateManyBeams),intent(inout)      ::      this            
            logical,intent(in)                          ::      is
            this%status = this%status - iand(this%status,LIB_IMB_FILE_READ)
            if (is) this%status = this%status + LIB_IMB_FILE_READ 
            return
        end subroutine setFileRead

        pure logical function isFileRead(this)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
            type(IntegrateManyBeams),intent(in)         ::      this            
            isFileRead = iand(this%status,LIB_IMB_FILE_READ) == LIB_IMB_FILE_READ
            return
        end function isFileRead

        pure subroutine setAtomSpaceSet(this,is)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
            type(IntegrateManyBeams),intent(inout)      ::      this          
            logical,intent(in)                          ::      is
            this%status = this%status - iand(this%status,LIB_IMB_ATOMSPACE_SET)
            if (is) this%status = this%status + LIB_IMB_ATOMSPACE_SET 
            return
        end subroutine setAtomSpaceSet

        pure logical function isAtomSpaceSet(this)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
            type(IntegrateManyBeams),intent(in)         ::      this            
            isAtomSpaceSet = iand(this%status,LIB_IMB_ATOMSPACE_SET) == LIB_IMB_ATOMSPACE_SET
            return
        end function isAtomSpaceSet

        pure subroutine setGvectorsSet(this,is)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
            type(IntegrateManyBeams),intent(inout)      ::      this          
            logical,intent(in)                          ::      is
            this%status = this%status - iand(this%status,LIB_IMB_GVECTORS_SET)
            if (is) this%status = this%status + LIB_IMB_GVECTORS_SET 
            return
        end subroutine setGvectorsSet

        pure logical function isGvectorsSet(this)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
            type(IntegrateManyBeams),intent(in)         ::      this            
            isGvectorsSet = iand(this%status,LIB_IMB_GVECTORS_SET) == LIB_IMB_GVECTORS_SET
            return
        end function isGvectorsSet

        pure subroutine setManyBeamSet(this,is)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
            type(IntegrateManyBeams),intent(inout)      ::      this          
            logical,intent(in)                          ::      is
            this%status = this%status - iand(this%status,LIB_IMB_MANYBEAM_SET)
            if (is) this%status = this%status + LIB_IMB_MANYBEAM_SET 
            return
        end subroutine setManyBeamSet

        pure logical function isManyBeamSet(this)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
            type(IntegrateManyBeams),intent(in)         ::      this            
            isManyBeamSet = iand(this%status,LIB_IMB_MANYBEAM_SET) == LIB_IMB_MANYBEAM_SET
            return
        end function isManyBeamSet

        pure subroutine setCStructFactSet(this,is)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
            type(IntegrateManyBeams),intent(inout)      ::      this          
            logical,intent(in)                          ::      is
            this%status = this%status - iand(this%status,LIB_IMB_CSTRUCTFACT_SET)
            if (is) this%status = this%status + LIB_IMB_CSTRUCTFACT_SET 
            return
        end subroutine setCStructFactSet

        pure logical function isCStructFactSet(this)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
            type(IntegrateManyBeams),intent(in)         ::      this            
            isCStructFactSet = iand(this%status,LIB_IMB_CSTRUCTFACT_SET) == LIB_IMB_CSTRUCTFACT_SET
            return
        end function isCStructFactSet

        pure subroutine setImageSpaceSet(this,is)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
            type(IntegrateManyBeams),intent(inout)      ::      this          
            logical,intent(in)                          ::      is
            this%status = this%status - iand(this%status,LIB_IMB_IMAGESPACE_SET)
            if (is) this%status = this%status + LIB_IMB_IMAGESPACE_SET 
            return
        end subroutine setImageSpaceSet

        pure logical function isImageSpaceSet(this)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
            type(IntegrateManyBeams),intent(in)         ::      this            
            isImageSpaceSet = iand(this%status,LIB_IMB_IMAGESPACE_SET) == LIB_IMB_IMAGESPACE_SET
            return
        end function isImageSpaceSet

    !     pure subroutine setPhaseFieldSet(this,is)
    ! !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !         type(IntegrateManyBeams),intent(inout)      ::      this          
    !         logical,intent(in)                          ::      is
    !         this%status = this%status - iand(this%status,LIB_IMB_PHASEFIELD_SET)
    !         if (is) this%status = this%status + LIB_IMB_PHASEFIELD_SET 
    !         return
    !     end subroutine setPhaseFieldSet

    !     pure logical function isPhaseFieldSet(this)
    ! !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !         type(IntegrateManyBeams),intent(in)         ::      this            
    !         isPhaseFieldSet = iand(this%status,LIB_IMB_PHASEFIELD_SET) == LIB_IMB_PHASEFIELD_SET
    !         return
    !     end function isPhaseFieldSet

        pure subroutine setSimResultSet(this,is)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
            type(IntegrateManyBeams),intent(inout)      ::      this           
            logical,intent(in)                          ::      is
            this%status = this%status - iand(this%status,LIB_IMB_SIM_RESULT)
            if (is) this%status = this%status + LIB_IMB_SIM_RESULT 
            return
        end subroutine setSimResultSet

        pure logical function isSimResultSet(this)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
            type(IntegrateManyBeams),intent(in)         ::      this            
            isSimResultSet = iand(this%status,LIB_IMB_SIM_RESULT) == LIB_IMB_SIM_RESULT
            return
        end function isSimResultSet

!-------
            
        subroutine Lib_IntegrateManyBeams_init_MPI()
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
#ifdef MPI
            integer             ::      ierror
            logical             ::      ok
            call MPI_Initialized(ok,ierror)
            if (.not. ok) call MPI_INIT(ierror)
            call MPI_COMM_SIZE(MPI_COMM_WORLD, nProcs, ierror)
            call MPI_COMM_RANK(MPI_COMM_WORLD, rank, ierror)
            if (rank == 0) print *,"Lib_IntegrateManyBeams::Lib_IntegrateManyBeams_init_MPI info - initialised with ",nProcs," processes"  
#else
            print *,"Lib_IntegrateManyBeams::Lib_IntegrateManyBeams_init_MPI info - serial mode"  
#endif            
            call Lib_AtomSpace_init_MPI()
            call Lib_ImagingSpace_init_MPI()
            call Lib_ComputePhaseFactor_init_MPI()
            return
        end subroutine Lib_IntegrateManyBeams_init_MPI


        subroutine errorExit(message)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
            character(len=*),intent(in),optional             ::      message
#ifdef MPI            
            integer         ::      ierror
#endif
            if (rank==0) then
                if (present(message)) print *,"Lib_IntegrateManyBeams::"//trim(message)
            end if
#ifdef MPI            
            call MPI_FINALIZE(ierror)
#endif
            stop
        end subroutine errorExit


    end module Lib_IntegrateManyBeams
