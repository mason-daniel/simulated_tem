
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
        use NBAX_StringTokenizers
        use Lib_RelativisticElectrons
        use Lib_RotationMatrices
        use Lib_Elements
        use Lib_ReadExtinctionDistances
        use Lib_DeformationGradients
        use Lib_FibonacciSphere
        use Lib_CrystalStructureFactor
        use Lib_RK4
#ifdef MPI
        use mpi_f08
#endif
        implicit none
        private

        real(kind=real64),parameter         ::      PI = 3.14159265390d0

        integer,private                     ::      rank = 0, nProcs = 1

        public          ::      Lib_IntegrateManyBeams_init_MPI

        public          ::      IntegrateManyBeams_ctor
        public          ::      report
        public          ::      delete
        public          ::      integrate
        public          ::      getImagingSpace
        public          ::      getIntensity
        public          ::      getA
        public          ::      getSigma
        public          ::      selectFoilTilt
        public          ::      applyFoilTilt
        public          ::      perfectLatticeIntensity
        public          ::      sanityCheck

        type,public     ::      phi_slice
            private
            complex(kind=real64),dimension(:,:,:),pointer     ::      phi                 !   (0:nG,lbx:ubx,lby:uby,lbz:ubz ) the solution 
        end type

        type,public     ::      IntegrateManyBeams
            private
            type(Lattice)                                       ::      latt
            type(AtomSpace)                                     ::      as
            type(ImagingSpace)                                  ::      is
            type(Gvectors),pointer                              ::      gv
            integer                                             ::      nPrec            !   how many sets used for precession averaging?
            type(ManyBeam),dimension(:),pointer                 ::      mb
            real(kind=real64)                                   ::      a , sigma           !   cell side , imaging blur
            real(kind=real64)                                   ::      omega0              !   volume per atom
            complex(kind=real64),dimension(:),pointer           ::      xi                  !   (1:nG2) complex extinction distances for all vectors g" = g - g'
            

            !complex(kind=real64),dimension(:,:,:,:),pointer     ::      x                   !   (0:3,nG,lbx:ubx,lby:uby,lbz:ubz) pointer to chunk of phase factors
            real(kind=real64),dimension(:,:,:,:,:),pointer      ::      grad_arg_x          !   (3,nG,lbx:ubx,lby:uby,lbz:ubz) pointer to chunk of phase factor gradients        
            real(kind=real64),dimension(:,:,:),pointer          ::      rho                 !   (lbx:ubx,lby:uby,lbz:ubz) pointer to chunk of atom densities

            integer                                             ::      nAtoms
            real(kind=real64),dimension(:,:),pointer            ::      r                   !   (3,nAtoms) positions of atoms

            type(phi_slice),dimension(:),pointer                ::      slice               !   (1:nPrec) the solution
            

            logical                                             ::      columnar            !   use columnar approximation
            logical                                             ::      lossy               !   use imaginary parts of crystal srtucture factor
        end type


        interface       IntegrateManyBeams_ctor
            module procedure            IntegrateManyBeams_null
            !module procedure            IntegrateManyBeams_ctor0
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

        interface           sanityCheck
            module procedure        sanityCheck0
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
            !this%gv = Gvectors_ctor()
            this%nPrec = 0
            nullify(this%mb)
            this%a = 1.0d0
            this%sigma = 1.0d0
            nullify(this%xi)
            !nullify(this%x)
            nullify(this%gv)
            nullify(this%grad_arg_x)
            nullify(this%rho)
            this%nAtoms = 0
            this%omega0 = 0
            nullify(this%r)
            this%columnar = .true.
            this%lossy = .true.
            nullify(this%slice)
            return
        end function IntegrateManyBeams_null


    !     function IntegrateManyBeams_ctor0( xifilename, latticename,a0_in,filename,V,nPrec,theta ) result(this)
    ! !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^           
    ! !*      default constructor
    !         character(len=*),intent(in)                 ::      xifilename              !   where to find crystal potential data
    !         character(len=*),intent(in)                 ::      latticename             !   eg "bcc" the name of the lattice type 
    !         real(kind=real64),dimension(3),intent(in)   ::      a0_in                   !   lattice parameter(s)
    !         character(len=*),intent(in)                 ::      filename                !   atom position file
    !         real(kind=real64),intent(in)                ::      V                       !   accelerating voltage (keV)
    !         integer,intent(in)                          ::      nPrec                   !   number of many beam copies to use for precession
    !         real(kind=real64),intent(in),optional       ::      theta                   !   precession angle
    !         type(IntegrateManyBeams)                    ::      this

    !         real(kind=real64),dimension(3,3)        ::      dfg_bar             !   average deformation gradient, used to fix correct positions of the g-vectors
    !         real(kind=real64),dimension(3,3)        ::      RR                  !   precession rotation matrix
    !         integer                                 ::      mm,nG
    !         real(kind=real64),dimension(3)          ::      kk                  !   electgron beam vector
    !         real(kind=real64),dimension(3)          ::      xyz_offset    
    !         real(kind=real64),dimension(3,3)        ::      a_super          !  periodic supercell lattice vectors
    !         real(kind=real64)                       ::      dphi
    !         integer,dimension(2,3)                  ::      bb
    !         this = IntegrateManyBeams_null()

    !         call setGvectors( this,a0_in,latticename,rho_in=2.0d0 )  
    !         call readXifile(this,xifilename)
    !         call readInputXyzFile( this,filename,dfg_bar,xyz_offset,a_super )    
    !         this%a = a0_in(1)/2
    !         this%sigma = a0_in(1)/2

    !         this%as = AtomSpace_ctor(this%a,xyz_offset,a_super)
    !         call setR(this%as,RotationMatrix_identity)
    !         call setThickness(this%as)
    !         this%omega0 = getOmega0(this%latt) * a0_in(1)**3
    
    !         if (rank==0) call report(this%as)
    !         if (rank==0) print *,""
    !         nG = getn(this%gv)
    !         call setR(this%gv,dfg_bar)    
    !         this%nPrec = max(1,nPrec)
    !         allocate(this%mb(this%nPrec))
    !         kk = (/0,0,1/) * ME * velocity( V*1000 )/ HBAR     
    !         this%mb(1) = ManyBeam_ctor(kk,this%gv,this%xi)
    !         if (this%nPrec > 1) then
    !             !   precession
    !             dphi = 2*PI / (this%nPrec-1)
    !             do mm = 2,this%nPrec
    !                 RR = RotationMatrix_ctor( (/0.0d0,1.0d0,0.0d0/),theta )                      !   rotation about y axis by angle theta
    !                 RR = matmul( RotationMatrix_ctor( (/0.0d0,0.0d0,1.0d0/),(mm-2)*dphi ),RR )   !   rotation about z axis by angle phi 
    !                 this%mb(mm) = ManyBeam_ctor( matmul(RR,kk),this%gv,this%xi )
    !             end do
    !         end if
    !         call setImagingSpace(this)                        
    !         call computePhaseFields(this)         
            
    !     !---    allocate memory for solution
    !         call getBounds(this%is,bb)
    !         do mm = 1,this%nPrec
    !             allocate(this%slice(mm))
    !             allocate(this%slice(mm)%phi(0:nG,bb(1,1)-1:bb(2,1)+1,bb(1,2)-1:bb(2,2)+1))
    !             this%slice(mm)%phi = 0
    !         end do
    
    !         return
    !     end function IntegrateManyBeams_ctor0


        function IntegrateManyBeams_ctor1( latticename,a0_in,filename,T,V,nPrec,theta,columnar,lossy,hkl_in ) result(this)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^                   
    !*      default constructor
    !        
            character(len=*),intent(inout)              ::      latticename             !   eg "bcc" the name of the lattice type 
            real(kind=real64),dimension(3),intent(inout)::      a0_in                   !   lattice parameter(s) - may be set to defaults if a0_in = 0
            character(len=*),intent(in)                 ::      filename                !   atom position file
            real(kind=real64),intent(in)                ::      T                       !   temperature (K)
            real(kind=real64),intent(in)                ::      V                       !   accelerating voltage (kV)
            integer,intent(in)                          ::      nPrec                   !   number of many beam copies to use for precession
            real(kind=real64),intent(in)                ::      theta                   !   precession angle
            logical,intent(in),optional                 ::      columnar                !   use columnar approximation
            logical,intent(in),optional                 ::      lossy                   !   use imaginary parts of crystal structure factor
            integer,dimension(:),intent(in),optional    ::      hkl_in                  !   two-beam mode
            type(IntegrateManyBeams)                    ::      this

            type(Gvectors)                                  ::      gv2                 !   expanded gvector set with all g" = g - g'
            character(len=XYZFILE_ATOMNAMELENGTH)           ::      element             !   name of element
            real(kind=real64),dimension(3,3)                ::      dfg_bar             !   average deformation gradient, used to fix correct positions of the g-vectors
            real(kind=real64),dimension(3)                  ::      xyz_offset          !   offset of atom positions read from xyz file
            real(kind=real64),dimension(3,3)                ::      a_super             !  periodic supercell lattice vectors
            real(kind=real64),dimension(3,3)                ::      RR                  !   precession rotation matrix

            real(kind=real64),dimension(3)                  ::      kk                  !   electron beam vector
            real(kind=real64)                               ::      dphi
            integer                                         ::      nG,nG2              !   number of g-vectors tracked
            integer,dimension(:,:),allocatable              ::      hkl                 !   reflections
            complex(kind=real64),dimension(:),allocatable   ::      Vg                  !   crystal structure factor
            integer                                         ::      ii,mm
            integer,dimension(2,3)                          ::      bb
           ! real(kind=real64)                               ::      LL
           ! real(kind=real64),dimension(:),allocatable      ::      Ig,Igm
             

            this = IntegrateManyBeams_null()
            if (present(columnar)) this%columnar = columnar
            if (present(lossy)) this%lossy = lossy

            if (rank==0) then
                print *,""
                print *,"reading input atomic positions file"
                print *,"^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^"
            end if
    

            call readInputXyzFile( this,filename,dfg_bar,xyz_offset,a_super,element )    
            
            ii = whichElement(element)
            if (ii == 0) then
                if (rank==0) print *,"Lib_IntegrateManyBeams::IntegrateManyBeams_ctor1 warning - do not recognise element """//trim(element)//""" - expect problems"
            end if

            if (any(a0_in<=0)) then            
                a0_in = getLatticeConstant(element)
                if (whichElement(element) == 0) then
                    call errorExit("error - can't determine lattice parameter from element name """//trim(element)//"""")
                else
                    if (rank==0) write(*,fmt='(a,3f12.5,a)') " Lib_IntegrateManyBeams::IntegrateManyBeams_ctor1 info - setting lattice parameter ",a0_in," from element name"
                end if
            end if

 
            if (latticename == UNKNOWN_LATTICE) then
                latticename = getLatticeName(element)
                if (latticename == UNKNOWN_LATTICE) then
                    call errorExit("error - can't determine lattice type from element name """//trim(element)//"""")
                else
                    if (rank==0) print *,"Lib_IntegrateManyBeams::IntegrateManyBeams_ctor1 info - setting lattice type """//trim(latticename)//""" from element name"
                end if
            end if

            
            if (rank==0) then
                print *,""
                print *,"setting g-vectors"
                print *,"^^^^^^^^^^^^^^^^^"
            end if
    

            if (present(hkl_in)) then
                call setGvectors( this,a0_in,latticename,V,hkl_in=hkl_in )  
            else
                call setGvectors( this,a0_in,latticename,V,rho_in=2.0d0 )  
            end if
            this%omega0 = getOmega0(this%latt) * a0_in(1)**3



            if (rank==0) then
                print *,""
                print *,"setting AtomSpace"
                print *,"^^^^^^^^^^^^^^^^^"
            end if

            this%a = a0_in(1)/4
            this%sigma = a0_in(1)/2
            this%as = AtomSpace_ctor(this%a,xyz_offset,a_super)
            call setR(this%as,RotationMatrix_identity)          
            call setThickness(this%as)
    
            if (rank==0) call report(this%as)            
            if (rank==0) print *,""




            if (rank==0) then
                print *,""
                print *,"computing crystal structure factors"
                print *,"^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^"
            end if


            !if (rank==0) print *,"Lib_IntegrateManyBeams::IntegrateManyBeams_ctor1 info - computing crystal structure factors"
            nG = getn(gv2)
            gv2 = Gvectors_ctor(this%gv,doubleSet=.true.)
            !call report(gv2)
            nG2 = getn(gv2)
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
                call computeCrystalStructureFactors( element,hkl, T,V*1000, Vg )       !    note computeCrystalStructureFactors wants accelerator voltage in V not kV.
            end if
#ifdef MPI
            call MPI_BCAST( Vg,nG2,MPI_DOUBLE_COMPLEX,0,MPI_COMM_WORLD,ii)
#endif        
            do ii = 1,nG2                
                this%xi(ii) = PI * HBAR * velocity( V*1000 ) / Vg(ii)
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
                print *,""
            end if




            if (rank==0) then
                print *,""
                print *,"initialising many-beam calculations"
                print *,"^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^"
            end if


            call setR(this%gv,dfg_bar)    
            call setR(gv2,dfg_bar)    
            this%nPrec = max(1,nPrec)
            allocate(this%mb(this%nPrec))
            kk = (/0,0,1/) * ME * velocity( V*1000.0d0 )/ HBAR     
            this%mb(1) = ManyBeam_ctor(kk,this%gv,gv2,this%xi)
            if (this%nPrec > 2) then
                !   precession
                dphi = 2*PI / (this%nPrec-1)
                do mm = 2,this%nPrec
                    RR = RotationMatrix_ctor( (/0.0d0,1.0d0,0.0d0/),theta )                      !   rotation about y axis by angle theta
                    RR = matmul( RotationMatrix_ctor( (/0.0d0,0.0d0,1.0d0/),(mm-2)*dphi ),RR )   !   rotation about z axis by angle phi 
                    this%mb(mm) = ManyBeam_ctor( matmul(RR,kk),this%gv,gv2,this%xi )
                end do
            end if


        !     if (rank==0) then
        !         print *,""
        !         print *,"perfect lattice intensity"
        !         print *,"^^^^^^^^^^^^^^^^^^^^^^^^^"
        !     end if
            
        ! !---    find the perfect lattice intensities
        !     LL = getThickness(this%as)
        !     !print *,"IntegrateManyBeams_ctor1 info - foil thickness ",LL
        !     nG = getn(this%gv)
        !     allocate(Ig(0:nG))
        !     allocate(Igm(0:nG))
        !     Ig = 0
        !     do mm = 1,this%nPrec
        !         call perfectLatticeIntensity(this%mb(mm), LL, Igm)
        !         Ig = Ig + Igm
        !     end do
        !     Ig = Ig / this%nPrec
        !     if (rank==0) then
        !         if (isMillerBravais(gv2)) then
        !             write(*,fmt='(a16,a16)') "reflection","Intensity"
        !         else
        !             write(*,fmt='(a12,a16)') "reflection","Intensity"
        !         end if
        !         do ii = 0,nG
        !             if (isMillerBravais(gv2)) then
        !                 write(*,fmt='(4i4)',advance="no") gethjkl(this%gv,ii)
        !             else
        !                 write(*,fmt='(3i4)',advance="no") gethkl(this%gv,ii)
        !             end if
        !             write(*,fmt='(f16.8)',advance="yes") Ig(ii)
        !         end do
        !         write(*,fmt='(a,f12.8)') " IntegrateManyBeams_ctor1 info - sum(I_g) (perfect lattice) ",sum(Ig)
        !     end if
            

            if (rank==0) then
                print *,""
                print *,"computing phase fields"
                print *,"^^^^^^^^^^^^^^^^^^^^^^"
            end if


            call setImagingSpace(this)                        
            call computePhaseFields(this)     
            

            if (rank==0) then
                print *,""
                print *,"allocating memory for solution"
                print *,"^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^"
            end if
            

        !---    allocate memory for solution
            nG = getn(this%gv)
            call getBounds(this%is,bb)
            allocate(this%slice(this%nPrec))
            do mm = 1,this%nPrec                
                allocate(this%slice(mm)%phi(0:nG,bb(1,1)-1:bb(2,1)+1,bb(1,2)-1:bb(2,2)+1))
                this%slice(mm)%phi = 0
            end do            
            if (rank==0) then
                print *,"Lib_IntegrateManyBeams::IntegrateManyBeams_ctor1 info - memory for phi ", size(this%slice(1)%phi)*16.0*this%nPrec/(1024*1024)," Mb"
                print *,""
            end if
             
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

            if (associated(this%grad_arg_x)) deallocate(this%grad_arg_x)
            
            if (associated(this%rho)) deallocate(this%rho)
      
            if (this%nAtoms>0) then
                deallocate(this%r)
            end if
            
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

        subroutine perfectLatticeIntensity0(this,hkl,I_g)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      compute the perfect lattice intensity estimate at reflection hkl
            type(IntegrateManyBeams),intent(in)         ::      this
            integer,dimension(:),intent(in)             ::      hkl    
            real(kind=real64),intent(out)               ::      I_g

            real(kind=real64)           ::      LL
            integer                     ::      ii,mm,nG,kk
            real(kind=real64),dimension(:),allocatable      ::      Ig,Igm

            ! if (rank==0) then
            !     print *,""
            !     print *,"perfect lattice intensity"
            !     print *,"^^^^^^^^^^^^^^^^^^^^^^^^^"
            ! end if
            
        !---    find the perfect lattice intensities
            LL = getThickness(this%as)
            if(rank==0) print *,"perfectLatticeIntensity0 info - foil thickness ",LL
            nG = getn(this%gv)
            allocate(Ig(0:nG))
            allocate(Igm(0:nG))
            Ig = 0
            do mm = 1,this%nPrec
                call perfectLatticeIntensity(this%mb(mm), LL, Igm)
                Ig = Ig + Igm
            end do
            Ig = Ig / this%nPrec
            if (rank==0) then
                if (isMillerBravais(this%gv)) then
                    write(*,fmt='(a16,a16)') "reflection","Intensity"
                    do ii = 0,nG
                        write(*,fmt='(4i4,f16.8)') gethjkl(this%gv,ii),Ig(ii)
                    end do
                else
                    write(*,fmt='(a12,a16)') "reflection","Intensity"
                    do ii = 0,nG
                        write(*,fmt='(3i4,f16.8)') gethkl(this%gv,ii),Ig(ii)
                    end do    
                end if
                ! do ii = 0,nG
                !     if (isMillerBravais(this%gv)) then
                !         write(*,fmt='(4i4)',advance="no") gethjkl(this%gv,ii)
                !     else
                !         write(*,fmt='(3i4)',advance="no") gethkl(this%gv,ii)
                !     end if
                !     write(*,fmt='(f16.8)',advance="yes") Ig(ii)
                ! end do
                write(*,fmt='(a,f12.8)') " Lib_IntegrateManyBeams::perfectLatticeIntensity0 info - sum(I_g) (perfect lattice) ",sum(Ig)
            end if

            kk = whichg(this%gv,hkl)
            I_g = Ig(kk)

            return
        end subroutine perfectLatticeIntensity0
            


        subroutine applyFoilTilt(this,R)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      apply a foil tilt R
    !*      this will change the G-vectors and the AtomSpace
    !*      and needs to be propagated to the ManyBeams
            type(IntegrateManyBeams),intent(inout)          ::      this
            real(kind=real64),dimension(3,3),intent(in)     ::      R
            real(kind=real64),dimension(3,3)     ::      R0,RR

            integer             ::      mm

            R0 = getR(this%gv)          !   old tilt
            RR = matmul(R,R0)           !   new tilt

            call setR(this%gv,RR)       !   update g-vectors
            call setR(this%as,RR)       !   update atom space
            do mm = 1,this%nPrec
                call setOrientationDependence(this%mb(mm))          !   mb has a pointer to gv
            end do
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
            allocate(I_tmp(0:nx-1,0:ny-1))
            I_tmp = I
            call MPI_REDUCE( I_tmp,I,nx*ny,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,ierror )
            deallocate(I_tmp)
#endif        

            I = I / this%nPrec

            return
        end subroutine getIntensity0


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


            real(kind=real64),dimension(3,3,nTrial) ::      RR
            real(kind=real64)                       ::      L_R
            integer                                 ::      ii,nG,kk
            real(kind=real64),dimension(:),allocatable    ::      Ig_R
            !integer                                 ::      besti
            real(kind=real64)                       ::      qq,Ik_R 

            integer,dimension(:),allocatable            ::      besti,besti_tmp
            real(kind=real64),dimension(:),allocatable  ::      bestq,bestq_tmp
        !---    which reflection are we looking at? Do we have it in our set?
            nG = getn(this%gv)
            allocate(Ig_R(0:nG))
            kk = whichg(this%gv,hkl)
            if (kk == LIB_GVECTORS_UNSET) call errorExit("selectFoilTilt error - could not find desired reflection in set")

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
                call perfectLatticeIntensity( this%mb(1), L_R, RR(:,:,ii), Ig_R )
                if (bright) then
                    Ik_R = Ig_R(kk)
                    Ig_R(kk) = 0.0d0                
                    qq = Ik_R - maxval(Ig_R(1:nG))              !   high if selected reflection is bright                    
                else
                    Ik_R = Ig_R(kk)
                    qq = 1 - Ik_R                               !   high if selected reflection is dark
                end if
                !print *,"rank ",rank," trial ",ii," quality ",qq,RR(:,3,ii),L_R
                if (qq > bestq(rank)) then
                    bestq(rank) = qq
                    besti(rank) = ii                    
                end if
            end do
            !print *,"rank ",rank," bestq ",bestq
#ifdef MPI
            allocate(besti_tmp(0:nProcs-1))
            besti_tmp = besti
            call MPI_ALLREDUCE(besti_tmp,besti,nProcs,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,ii )
            deallocate(besti_tmp)
            allocate(bestq_tmp(0:nProcs-1))
            bestq_tmp = bestq
            call MPI_ALLREDUCE(bestq_tmp,bestq,nProcs,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ii )
            deallocate(bestq_tmp)
#endif
            !if (rank==0) print *,"bestq ",bestq
            ii = maxloc(bestq,dim=1)-1      !   which rank found the best quality?
            qq = bestq(ii)
           ! if (rank == 0) print *,"besti ",qq,ii,besti(ii)
            ii = besti(ii)                  !   which trial was it?
            R(:,3) = RR(:,3,ii)
            R(:,1) = (/1,0,0/)
            call completeBasis( R(:,3),R(:,1),R(:,2) )
            L_R = getThickness( this%as,R )
            if (rank==0) print *,R,L_R
            call perfectLatticeIntensity( this%mb(1), L_R, R, Ig_R )
            if (bright) then
                Ik_R = Ig_R(kk)
                Ig_R(kk) = 0.0d0                
                qq = Ik_R - maxval(Ig_R(1:nG))              !   high if selected reflection is bright                    
            else
                Ik_R = Ig_R(kk)
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



        subroutine integrate0(this)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      propagate many beam equations in z through the system
            type(IntegrateManyBeams),intent(inout)     ::      this

            type(RK4)           ::      integrator
            complex(kind=real64),dimension(:),pointer       ::      rk_phip
            complex(kind=real64),dimension(:),pointer       ::      rk_dphip
            complex(kind=real64),dimension(:,:,:),pointer   ::      dphidz
            
            integer,dimension(2,3)                          ::      bb
            integer             ::      iz,rk_step
            integer             ::      mm
            integer             ::      nx,ny,ng
            integer             ::      my_nx,my_ny
            !integer             ::      ix,iy
            integer             ::      dz
            !logical,dimension(:,:),allocatable              ::      mask

        !---    allocate memory and temporaries needed for RK4 integration
            nx = getNx(this%is)
            ny = getNy(this%is)
            nG = getn(this%gv)
            call getBounds(this%is,bb)
            bb(1,:) = bb(1,:) - 1           !   buffer pixel
            bb(2,:) = bb(2,:) + 1           !   buffer pixel
            my_nx = bb(2,1)+1-bb(1,1)
            my_ny = bb(2,2)+1-bb(1,2)

            integrator = RK4_ctor( n=(nG+1) * my_nx * my_ny, deltaz=2*geta(this%is) )     !   dz = 2a because the RK4 integrator needs to evaluate at half steps
            if (rank==0) call report(integrator)
            if (rank==0) print *,"Lib_IntegrateManyBeams::integrate0 info - bounds: nx,ny = ",nx,ny," lbx,ubx,lby,uby = ",bb(1:2,1),bb(1:2,2)," my_nx,my_ny = ",my_nx,my_ny
            !if ( (bb(2,1)+1-bb(1,1)/=nx) .or. (bb(2,2)+1-bb(1,2)/=ny) ) then                
            !    call errorExit("error")
            !end if

        !---    get pointers to current state and derivative
            call getphip_dphip(integrator , rk_phip,rk_dphip)
           ! print *,"a"
            allocate(dphidz(0:nG,bb(1,1):bb(2,1),bb(1,2):bb(2,2)))
            !print *,"b"
            dphidz = 0      
            rk_phip = 0
            rk_dphip = 0

        !---    set boundary conditions at z=0
            do mm = 1,this%nPrec
                this%slice(mm)%phi(0   ,:,:) = 1.0d0
                this%slice(mm)%phi(1:nG,:,:) = 0.0d0
            end do
            ! do mm = 1,this%nPrec
            !     print *,"c ",mm,size(this%slice(mm)%phi,dim=1),size(this%slice(mm)%phi,dim=2),size(this%slice(mm)%phi,dim=3)
            ! end do
            !allocate(mask())

!            rk_phip = pack(this%phi,.true.)
             

            
        !---    integrate!   
            do iz = 0,getNz(this%is)-3,2

                do mm = 1,this%nPrec
                    !print *,"d",iz,mm
                    ! print *,"size(rk_phip) ",size(rk_phip),(nG+1) * my_nx * my_ny,size(this%slice(mm)%phi)
                    rk_phip = pack(this%slice(mm)%phi,.true.)
                    !print *,"e"
                    call update(integrator)
                    !print *,"f"
                    dz = 0
                    do rk_step = 1,4
                        ! print *,"integrate0 slice ",rank
                    !   ensure that phi is distributed correctly 
                        this%slice(mm)%phi(:,bb(1,1):bb(2,1),bb(1,2):bb(2,2)) = reshape( rk_phip, (/nG+1,my_nx,my_ny/) )
                        ! print *,"integrate0 sendrecv ",rank
                        if (.not. this%columnar) call sendrecv( this%is,this%slice(mm)%phi )
                        ! print *,"integrate0 sendrecv done ",rank
                    !   compute dphidz
                        ! print *,"integrate0 sendrecv ",rank," ass ",associated(this%slice(mm)%phi),associated(dphidz)
                        call finddPhidz( this%mb(mm) , this%slice(mm)%phi , iz+dz, dphidz, this%columnar )
                        rk_dphip = pack(dphidz,.true.)
                                            
                        ! print *,"integrate0 update ",rank
                    !   update accumulators in RK4
                        call update(integrator,rk_step,dz)
                       
                    end do

                    !if (iz == getNz(this%is)-2) then
                        !---    extract result
                        this%slice(mm)%phi(0:nG,bb(1,1):bb(2,1),bb(1,2):bb(2,2)) = reshape( rk_phip, (/nG+1,my_nx,my_ny/) )
                    !end if

                end do

                !stop

            end do

        

            call delete(integrator)
            return

        end subroutine integrate0

!-------

        subroutine setGvectors( this,a0_in,latticename,V,rho_in,hkl_in )
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      given the name of the lattice, and a cutoff radius rho for 
    !*      the reflections, construct a set of g-vectors.

            type(IntegrateManyBeams),intent(inout)          ::      this
            real(kind=real64),dimension(3),intent(in)       ::      a0_in           !   lattice parameter(s) - expect 3 for hcp
            character(len=*),intent(in)                     ::      latticename
            real(kind=real64),intent(in),optional           ::      V               !   electron voltage (keV)
            real(kind=real64),intent(in),optional           ::      rho_in          !   compute all g-vectors within given radius
            integer,dimension(3),intent(in),optional        ::      hkl_in          !   compute a single reflection

            real(kind=real64),dimension(3,3)        ::      a_conventional_cell
            integer                                 ::      nn
            integer,dimension(:,:),allocatable      ::      hkl,hjkl
            logical,dimension(:),allocatable        ::      keep
            integer,dimension(3)                    ::      hkl_swp
            real(kind=real64)                       ::      a0,sg
            real(kind=real64),dimension(3)          ::      kk,gg
            integer                                 ::      ii,jj
            logical                 ::      done

            if (rank==0) then
                if (present(rho_in)) then
                    print *,"Lib_IntegrateManyBeams::setGvectors info - many beam mode, radius ",rho_in
                else if (present(hkl_in)) then
                    print *,"Lib_IntegrateManyBeams::setGvectors info - two beam mode, hkl = ",hkl_in
                end if
            end if
            this%latt = Lattice_ctor(latticename)
            
            if (getLatticeType(this%latt) == LATTICE_HCP) then
                a0 = a0_in(1)
                call setCoverA(this%latt,a0_in(3)/a0_in(1))                 
            else
                a0 = a0_in(1)
            end if
            if (rank==0) call report(this%latt)

            a_conventional_cell = a0*getConventionalCell(this%latt)

            if (present(rho_in)) then
            !---    find all reflections within certain radius
                call permittedReflections( this%latt,nn,hkl,rho_in )            
            else if (present(hkl_in)) then
                nn = 2
                allocate(hkl(3,nn))
                hkl(:,1) = 0
                hkl(:,2) = hkl_in(:)
            end if
 

        !---    check the deviation parameter of the reflections
            kk = (/0,0,1/) * ME * velocity( V*1000.0d0 )/ HBAR     
            allocate(keep(nn))
            do ii = 1,nn
                gg = getg( a_conventional_cell,hkl(:,ii) )
                sg = deviationParameter( kk,gg )
                
                keep(ii) = (abs(2*PI*sg*a0) < 0.5d0) 
                if (rank == 0) then
                    if (keep(ii)) then
                        write(*,fmt='(a,3i4,a,f12.5,a)') " Lib_IntegrateManyBeams::setGvectors info - permitted reflection ",hkl(:,ii)," sg = ",sg
                    else
                        write(*,fmt='(a,3i4,a,f12.5,a)') " Lib_IntegrateManyBeams::setGvectors info - permitted reflection ",hkl(:,ii)," sg = ",sg," blocked "                    
                    end if
                end if
            end do
            jj = 0
            do ii = 1,nn
                if (keep(ii)) then
                    jj = jj + 1
                    hkl(:,jj) = hkl(:,ii)
                end if
            end do
            nn = count(keep)
            

 

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
            allocate(this%gv)
            if (getLatticeType(this%latt) == LATTICE_HCP) then
                allocate(hjkl(4,nn))
                do ii = 1,nn
                    hjkl(:,ii) = nint( MillerToMillerBravais_plane(real(hkl(:,ii),kind=real64)) )
                end do
                this%gv = Gvectors_ctor(a_conventional_cell,hjkl(:,2:nn))
            else
                this%gv = Gvectors_ctor(a_conventional_cell,hkl(:,2:nn))
            end if
 
 
            if (rank==0) print *,""
            return
        end subroutine setGvectors
 
 



    !     subroutine readXifile(this,xifilename)
    ! !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    ! !*      when the g-vectors are set, can read in the xi data
    !         type(IntegrateManyBeams),intent(inout)          ::      this
    !         character(len=*),intent(in)                     ::      xifilename

    !         integer,dimension(:,:),allocatable                ::      hkl
    !         integer             ::              ii,nG
    !         logical             ::              ok

    !         if (rank==0) print *,"Lib_IntegrateManyBeams::readXifile info"
    !         nG = getn(this%gv)
    !         allocate(this%xi(0:nG))
    !         allocate(hkl(3,0:nG))
    !         hkl(:,0) = 0
    !         do ii = 1,nG
    !             hkl(:,ii) =  gethkl(this%gv,ii)
    !         end do
    !         inquire(file=trim(xifilename),exist=ok)
    !         if (.not. ok) call errorExit("Lib_IntegrateManyBeams::readXifile ERROR - file not found """//trim(xifilename)//"""")
    !         call readExtinctionDistance( xifilename,hkl,getSymmetry(this%latt),this%xi, ok )
    !         if (.not. ok) call errorExit("Lib_IntegrateManyBeams::readXifile ERROR - error reading file """//trim(xifilename)//"""")

    !         do ii = 0,nG
    !             this%xi(ii) = complex( -abs(real(this%xi(ii))) , -abs(aimag(this%xi(ii))) )
    !         end do

    !         if (rank==0) print *,""
                
    !         return
    !     end subroutine readXifile




        subroutine readInputXyzFile( this,filename,dfg_bar,xyz_offset,a_super,element )
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      reads in the atom positions ( shared by all processes )
    !*      and establishes the AtomSpace
    !*      returns the average deformation gradient, offset of atoms in file, and periodic supercell lattice vectors
     
            type(IntegrateManyBeams),intent(inout)              ::      this
            character(len=*),intent(in)                         ::      filename
            real(kind=real64),dimension(3),intent(out)          ::      xyz_offset  !   atom position offset read from file
            real(kind=real64),dimension(3,3),intent(out)        ::      a_super
            real(kind=real64),dimension(3,3),intent(inout)      ::      dfg_bar
            character(len=XYZFILE_ATOMNAMELENGTH),intent(out),optional   ::      element

            type(XYZFile)       ::      xyz
            logical             ::      ok,ok_tmp
            !real(kind=real64),dimension(3,3)        ::  a_super
            !real(kind=real64),dimension(3)          ::  xyz_offset
            
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
                    this%r = this%r + 1.0d-8
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
                                    
                        call computeAvgDefGrad(eps_bar,dfg_bar,weight)
                        do ii = 1,this%nAtoms
                            if (  nint( this%r(4,ii) ) == kk) call computeAvgDefGrad( this%r(5:13,ii),1.0d0, eps_bar,dfg_bar,weight)
                        end do
                        call computeAvgDefGrad(eps_bar,dfg_bar,weight , dfg) ; dfg_bar = dfg
                        print *,"d2bi info - average deformation gradient, strain, rot read from extended .xyz file "

                        call DefGradToStrainAndRotMat(dfg_bar,eps_bar,rot_bar)
        
                    else        !   don't have columns data
                        if (rank==0) then
                            print *,"Lib_IntegrateManyBeams::readInputXyzFile info - don't have deformation gradient data in input file"
                            print *,"   assuming lattice orientation aligned with xyz coordinates"
                        end if
                        dfg_bar = reshape( (/1,0,0,0,1,0,0,0,1/),(/3,3/) )
                        eps_bar = 0
                        rot_bar = reshape( (/1,0,0,0,1,0,0,0,1/),(/3,3/) )
                    
                    end if
                    if (rank==0) then
                        write (*,fmt='(a36,a,a36,a,a36)') "average def grad","    ","average strain","    ","average rotation"
                        write (*,fmt='(3f12.6,a,3f12.6,a,3f12.6)') dfg_bar(1,:),"    ",eps_bar(1,:),"    ",rot_bar(1,:)
                        write (*,fmt='(3f12.6,a,3f12.6,a,3f12.6)') dfg_bar(2,:),"    ",eps_bar(2,:),"    ",rot_bar(2,:)
                        write (*,fmt='(3f12.6,a,3f12.6,a,3f12.6)') dfg_bar(3,:),"    ",eps_bar(3,:),"    ",rot_bar(3,:)
                    end if   

                !---    attempt to id the most populous element
                    if (present(element)) then
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
                        element = getAtomName(xyz,mostPopType)
                        print *,"Lib_IntegrateManyBeams::readInputXyzFile info - most common atom type """//trim(element)//""""
                    end if

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
            if (present(element)) call MPI_BCAST(element,XYZFILE_ATOMNAMELENGTH,MPI_CHARACTER,0,MPI_COMM_WORLD,ii)
            if (rank/=0) allocate(this%r(3,this%nAtoms))
            call MPI_BCAST(this%r(1:3,1:this%nAtoms),3*this%nAtoms,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ii)                
            call MPI_BCAST(a_super,9,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ii)
            call MPI_BCAST(xyz_offset,3,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ii)
            call MPI_BCAST(dfg_bar,9,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ii)
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

            type(IntegrateManyBeams),intent(inout)          ::      this

            integer             ::     Nx,Ny,Nz
            real(kind=real64)   ::      theta           !   max angle required for beam dispersion
            integer             ::      ii,jj,mm

            real(kind=real64),dimension(:,:),pointer        ::      rt_tmp
            real(kind=real64),dimension(3,9)                ::      xtp
            integer                                         ::      np
            integer                                         ::      mynAtoms
            if (rank==0) print *,"Lib_IntegrateManyBeams::setImagingSpace"

            theta = 0
            do mm = 1,this%nPrec
                do ii = 1,getn(this%gv)
                    theta = max(theta,getAngle(this%mb(mm),ii))
                   ! print *,"Lib_IntegrateManyBeams::setImagingSpace info - angle prec,g ",mm,ii," = ",getAngle(this%mb(mm),ii)
                end do
            end do
            
            call suggestImagingSpace(this%as,Nx,Ny,Nz,theta)
            if (rank==0) print *,"Lib_IntegrateManyBeams::setImagingSpace info - theta = ",theta," Nx,Ny,Nz = ",Nx,Ny,Nz
            call setDelta(this%as,Nx,Ny,Nz)

            this%is = ImagingSpace_ctor(this%a,this%sigma, getdelta(this%as),Nx,Ny,Nz)
            if (rank==0) call report(this%is) 
             

        !---    can now make my transformation on the atom positions, and store all the periodic copies I need.
        !       count the number of atoms after periodic copies
            mynAtoms = 0
            do ii = 1,this%nAtoms
                call periodicCopies(this%as,this%is,this%r(:,ii),np,xtp)
                do jj = 1,np
                    if (inMyCell(this%is,xtp(:,jj),buffered=.true.)) mynAtoms = mynAtoms + 1
                end do
            end do
            print *,"Lib_IntegrateManyBeams::setImagingSpace info - rank ",rank," mynAtoms = ",mynAtoms,"/",this%nAtoms

            allocate(rt_tmp(3,mynAtoms))
            mynAtoms = 0
            do ii = 1,this%nAtoms
                call periodicCopies(this%as,this%is,this%r(:,ii),np,xtp)
                do jj = 1,np
                    if (inMyCell(this%is,xtp(:,jj),buffered=.true.)) then
                        mynAtoms = mynAtoms + 1
                        rt_tmp(1:3,mynAtoms) = xtp(1:3,jj)
                    end if
                end do
            end do

            deallocate(this%r)
            this%nAtoms = mynAtoms
            this%r => rt_tmp
            if (rank==0) print *,""

            return
        end subroutine setImagingSpace


        subroutine computePhaseFields(this)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      once the imaging space is set, I can compute the phase factors in my region
    !*      note that this subroutine does a lot of the big memory allocations
            type(IntegrateManyBeams),intent(inout)          ::      this
            integer,dimension(2,3)          ::      bb
            integer                         ::      nG
            real(kind=real64),dimension(:,:),allocatable        ::      gg          !   (3,nG)      g-vectors
            !integer,dimension(:),allocatable        ::      gindx
            integer                                 ::      ii,mm,nGcalc
             
            real                            ::      mem

            if (rank==0) print *,"Lib_IntegrateManyBeams::computePhaseFields"

            call getBounds(this%is,bb)
            nG = getn(this%gv)
            allocate(this%rho   (     bb(1,1):bb(2,1),bb(1,2):bb(2,2),bb(1,3):bb(2,3)))        
            mem = (bb(2,1)+1-bb(1,1))*(bb(2,2)+1-bb(1,2)) 

            if (rank==0) then
                print *,"Lib_IntegrateManyBeams::computePhaseFields info - memory for x,rho ", (4*nG+0.5)*mem*16.0/(1024*1024)," Mb"
                print *,"Lib_IntegrateManyBeams::computePhaseFields info - bb = ",bb," nG = ",nG
            end if


        !---    we don't need to compute phase factors for all the g-vectors.
        !       g=0 gives a constant, x(g=0) = exp[ -i g.u ] = 1
        !       and x(-g) =  exp[ -i (-g).u ] = x(g)*
            nGcalc = nPositiveg(this%gv)
            allocate(gg(3,nGcalc))

            nGcalc = 0
            do ii = 1,nG
                if (isPositiveg(this%gv,ii)) then
                    nGcalc = nGcalc + 1
                    gg(:,nGcalc) = getG(this%gv,ii)
                    if (rank==0) write(*,fmt='(3(a,3f10.5))') " Lib_IntegrateManyBeams::computePhaseFields info - calculating x_(g=",getG(this%gv,ii),")"
                else
                    !   already have computed -g, so don't need to compute g
                    if (rank==0) write(*,fmt='(3(a,3f10.5))') " Lib_IntegrateManyBeams::computePhaseFields info - skipping calculation of x_(g=",getG(this%gv,ii),") = x*_(g=",getG(this%gv,getMinusg(this%gv,ii)),")"
                end if                
            end do

            allocate(this%grad_arg_x(3,nGcalc,bb(1,1):bb(2,1),bb(1,2):bb(2,2),bb(1,3):bb(2,3)))

            call computePhaseFactor( this%nAtoms,this%r,gg, this%is,this%grad_arg_x,this%rho )

            this%rho = this%rho/this%a**3
            if (rank==0) print *,"Lib_IntegrateManyBeams::computePhaseFields info - minmaxavg rho      ",minval(abs(this%rho)),maxval(abs(this%rho)),sum(abs(this%rho))/(size(this%rho))," (1/A^3) "
            this%rho = densityScale( this%rho, this%omega0 )
            if (rank==0) print *,"Lib_IntegrateManyBeams::computePhaseFields info - minmaxavg rho'     ",minval(abs(this%rho)),maxval(abs(this%rho)),sum(abs(this%rho))/(size(this%rho))

        !---    unpack the computed g-vectors




            if (rank==0) print *,"Lib_IntegrateManyBeams::computePhaseFields info - minmaxavg |grad_arg_x| ",minval(abs(this%grad_arg_x)),maxval(abs(this%grad_arg_x)),sum(abs(this%grad_arg_x))/(size(this%grad_arg_x))
            
            do mm = 1,this%nPrec
                call setPhaseFactors( this%mb(mm),this%a, this%grad_arg_x,this%rho )
            end do
            

            return
        end subroutine computePhaseFields


        pure logical function sanityCheck0(this)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      returns true if all the pointers are set
            type(IntegrateManyBeams),intent(in)           ::      this            
            integer         ::      mm

            sanityCheck0 = (getn(this%gv)>0)

            !sanityCheck0 = sanityCheck0 .and. associated(this%x)
            sanityCheck0 = sanityCheck0 .and. associated(this%grad_arg_x)
            sanityCheck0 = sanityCheck0 .and. associated(this%rho)

            do mm = 1,this%nPrec
                sanityCheck0 = sanityCheck0 .and. sanityCheck(this%mb(mm))
            end do
            

            return
        end function sanityCheck0

        pure real(kind=real64) function getSigma0(this)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      returns width of blurring used to construct phase fields
            type(IntegrateManyBeams),intent(in)           ::      this            
            getSigma0 = this%sigma
            return
        end function getSigma0

        pure real(kind=real64) function getA0(this)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      returns pixel size
            type(IntegrateManyBeams),intent(in)           ::      this            
            getA0 = this%a
            return
        end function getA0




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
