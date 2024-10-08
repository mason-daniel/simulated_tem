
    module Lib_IntegrateManyBeams
!---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
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
        use Lib_CrystalStructureFactor
#ifdef MPI
        use mpi_f08
#endif
        implicit none
        private

        real(kind=real64),parameter                 ::      PI = 3.14159265390d0

        integer,private                 ::      rank = 0, nProcs = 1

        public          ::      Lib_IntegrateManyBeams_init_MPI

        public          ::      IntegrateManyBeams_ctor
        public          ::      report
        public          ::      delete



        type,public     ::      IntegrateManyBeams
            private
            type(Lattice)                                       ::      latt
            type(AtomSpace)                                     ::      as
            type(ImagingSpace)                                  ::      is
            type(Gvectors)                                      ::      gv
            integer                                             ::      nPrec            !   how many sets used for precession averaging?
            type(ManyBeam),dimension(:),pointer                 ::      mb
            real(kind=real64)                                   ::      a , sigma           !   cell side , imaging blur

            complex(kind=real64),dimension(:),pointer           ::      xi                  !   (0:nG) complex extimnction distances
            

            complex(kind=real64),dimension(:,:,:,:),pointer     ::      x                   !   (nG,lbx:ubx,lby:uby,lbz:ubz) pointer to chunk of phase factors
            complex(kind=real64),dimension(:,:,:,:,:),pointer   ::      grad_x              !   (3,nG,lbx:ubx,lby:uby,lbz:ubz) pointer to chunk of phase factor gradientss            
            real(kind=real64),dimension(:,:,:),pointer          ::      rho                 !   (lbx:ubx,lby:uby,lbz:ubz) pointer to chunk of atom densities

            integer                                             ::      nAtoms
            real(kind=real64),dimension(:,:),pointer            ::      r                   !   (3,nAtoms) positions of atoms
        end type


        interface       IntegrateManyBeams_ctor
            module procedure            IntegrateManyBeams_null
            module procedure            IntegrateManyBeams_ctor0
            module procedure            IntegrateManyBeams_ctor1
        end interface

        interface       report
            module procedure            report0
        end interface

        interface       delete
            module procedure            delete0
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
            this%gv = Gvectors_ctor()
            this%nPrec = 0
            nullify(this%mb)
            this%a = 1.0d0
            this%sigma = 1.0d0
            nullify(this%xi)
            nullify(this%x)
            nullify(this%grad_x)
            nullify(this%rho)
            this%nAtoms = 0
            nullify(this%r)
            
            return
        end function IntegrateManyBeams_null


        function IntegrateManyBeams_ctor0( xifilename, latticename,a0_in,filename,V,nPrec,theta ) result(this)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^                      
    !*      default constructor
            character(len=*),intent(in)                 ::      xifilename              !   where to find crystal potential data
            character(len=*),intent(in)                 ::      latticename             !   eg "bcc" the name of the lattice type 
            real(kind=real64),dimension(3),intent(in)   ::      a0_in                   !   lattice parameter(s)
            character(len=*),intent(in)                 ::      filename                !   atom position file
            real(kind=real64),intent(in)                ::      V                       !   accelerating voltage (keV)
            integer,intent(in)                          ::      nPrec                   !   number of many beam copies to use for precession
            real(kind=real64),intent(in),optional       ::      theta                   !   precession angle
            type(IntegrateManyBeams)                    ::      this

            real(kind=real64),dimension(3,3)        ::      dfg_bar             !   average deformation gradient, used to fix correct positions of the g-vectors
            real(kind=real64),dimension(3,3)        ::      RR                  !   precession rotation matrix
            integer                                 ::      mm
            real(kind=real64),dimension(3)          ::      kk                  !   electgron beam vector
            real(kind=real64)                       ::      dphi
            this = IntegrateManyBeams_null()

            call setLattice( this,a0_in,latticename,rho_in=2.0d0 )  
            call readXifile(this,xifilename)
            call readInputXyzFile( this,filename,dfg_bar )    
            call setR(this%gv,dfg_bar)    
            this%nPrec = max(1,nPrec)
            allocate(this%mb(this%nPrec))
            kk = (/0,0,1/) * ME * velocity( V/1000 )/ HBAR     
            this%mb(1) = ManyBeam_ctor(kk,this%gv,this%xi)
            if (this%nPrec > 1) then
                !   precession
                dphi = 2*PI / (this%nPrec-1)
                do mm = 2,this%nPrec
                    RR = RotationMatrix_ctor( (/0.0d0,1.0d0,0.0d0/),theta )                      !   rotation about y axis by angle theta
                    RR = matmul( RotationMatrix_ctor( (/0.0d0,0.0d0,1.0d0/),(mm-2)*dphi ),RR )   !   rotation about z axis by angle phi 
                    this%mb(mm) = ManyBeam_ctor( matmul(RR,kk),this%gv,this%xi )
                end do
            end if
            call setImagingSpace(this)                        
            call computePhaseFields(this)                     
 
            return
        end function IntegrateManyBeams_ctor0


        function IntegrateManyBeams_ctor1( element, latticename,a0_in,filename,T,V,nPrec,theta ) result(this)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^                     
    !*      default constructor
            character(len=*),intent(in)                 ::      element                 !   name of element
            character(len=*),intent(in)                 ::      latticename             !   eg "bcc" the name of the lattice type 
            real(kind=real64),dimension(3),intent(in)   ::      a0_in                   !   lattice parameter(s)
            character(len=*),intent(in)                 ::      filename                !   atom position file
            real(kind=real64),intent(in)                ::      T                       !   temperature (K)
            real(kind=real64),intent(in)                ::      V                       !   accelerating voltage (kV)
            integer,intent(in)                          ::      nPrec                   !   number of many beam copies to use for precession
            real(kind=real64),intent(in)                ::      theta                   !   precession angle
            type(IntegrateManyBeams)                    ::      this

            real(kind=real64),dimension(3,3)        ::      dfg_bar             !   average deformation gradient, used to fix correct positions of the g-vectors
            real(kind=real64),dimension(3,3)        ::      RR                  !   precession rotation matrix
            integer                                 ::      mm
            real(kind=real64),dimension(3)          ::      kk                  !   electgron beam vector
            real(kind=real64)                       ::      dphi
            integer                                 ::      nG                  !   number of g-vectors tracked
            integer,dimension(:,:),allocatable      ::      hkl                 !   reflections
            complex(kind=real64),dimension(:),allocatable   ::      Vg          !   crystal structure factor
            integer                                 ::      ii

            this = IntegrateManyBeams_null()

            call setLattice( this,a0_in,latticename,rho_in=2.0d0 )  



            if (rank==0) print *,"Lib_IntegrateManyBeams::IntegrateManyBeams_ctor1 info - computing crystal structure factors"
            nG = getn(this%gv)
            allocate(this%xi(0:nG))
            allocate(hkl(3,0:nG))
            hkl(:,0) = 0
            do ii = 1,nG
                hkl(:,ii) =  gethkl(this%gv,ii)
            end do
            allocate(Vg(0:nG))
            call computeCrystalStructureFactors( element,hkl, T,V*1000.0d0, Vg )       !    note computeCrystalStructureFactors wants accelerator voltage in V not kV.
            do ii = 0,nG
                this%xi(ii) = PI * HBAR * velocity( V ) / Vg(ii)
            end do
            if (rank==0) print *,""


            !call readXifile(this,xifilename)
            call readInputXyzFile( this,filename,dfg_bar )    
            call setR(this%gv,dfg_bar)    
            this%nPrec = max(1,nPrec)
            allocate(this%mb(this%nPrec))
            kk = (/0,0,1/) * ME * velocity( V/1000 )/ HBAR     
            this%mb(1) = ManyBeam_ctor(kk,this%gv,this%xi)
            if (this%nPrec > 1) then
                !   precession
                dphi = 2*PI / (this%nPrec-1)
                do mm = 2,this%nPrec
                    RR = RotationMatrix_ctor( (/0.0d0,1.0d0,0.0d0/),theta )                      !   rotation about y axis by angle theta
                    RR = matmul( RotationMatrix_ctor( (/0.0d0,0.0d0,1.0d0/),(mm-2)*dphi ),RR )   !   rotation about z axis by angle phi 
                    this%mb(mm) = ManyBeam_ctor( matmul(RR,kk),this%gv,this%xi )
                end do
            end if
            call setImagingSpace(this)                        
            call computePhaseFields(this)                     
 
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
            deallocate(this%xi)
            do mm = 1,this%nPrec
                call delete(this%mb(mm))
            end do
            deallocate(this%mb)
            if (associated(this%x)) then
                deallocate(this%x)
                deallocate(this%grad_x)
                deallocate(this%rho)
            end if
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
                call report(this%mb(mm))
            end do
            return
        end subroutine report0

    !---

!-------

        subroutine setLattice( this,a0_in,latticename,rho_in )
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ 
    !*      given the name of the lattice, and a cutoff radius rho for 
    !*      the reflections, construct a set of g-vectors.

            type(IntegrateManyBeams),intent(inout)          ::      this
            real(kind=real64),dimension(3),intent(in)       ::      a0_in           !   lattice parameter(s) - expect 3 for hcp
            character(len=*),intent(in)                     ::      latticename
            real(kind=real64),intent(in)                    ::      rho_in

            real(kind=real64),dimension(3,3)        ::      a_conventional_cell
            integer                                 ::      nn
            integer,dimension(:,:),allocatable      ::      hkl,hjkl
            real(kind=real64)                       ::      a0
            integer                                 ::      ii

            if (rank==0) print *,"Lib_IntegrateManyBeams::setLattice info"
            this%latt = Lattice_ctor(latticename)
            
            if (getLatticeType(this%latt) == LATTICE_HCP) then
                a0 = a0_in(1)
                call setCoverA(this%latt,a0_in(3)/a0_in(1))                 
            else
                a0 = a0_in(1)
            end if
            if (rank==0) call report(this%latt)

            a_conventional_cell = a0*getConventionalCell(this%latt)
            call permittedReflections( this%latt,nn,hkl,rho_in )

    !---    construct a g-vector object
            if (getLatticeType(this%latt) == LATTICE_HCP) then
                allocate(hjkl(4,nn))
                do ii = 1,nn
                    hjkl(:,ii) = nint( MillerToMillerBravais_plane(real(hkl(:,ii),kind=real64)) )
                end do
                this%gv = Gvectors_ctor(a_conventional_cell,hjkl)
            else
                this%gv = Gvectors_ctor(a_conventional_cell,hkl)
            end if
 
            if (rank==0) call report(this%gv)
 
            this%a = a0/2
            this%sigma = a0/2
            if (rank==0) print *,""
            return
        end subroutine setLattice
 




        subroutine readXifile(this,xifilename)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      when the g-vectors are set, can read in the xi data
            type(IntegrateManyBeams),intent(inout)          ::      this
            character(len=*),intent(in)                     ::      xifilename

            integer,dimension(:,:),allocatable                ::      hkl
            integer             ::              ii,nG
            logical             ::              ok

            if (rank==0) print *,"Lib_IntegrateManyBeams::readXifile info"
            nG = getn(this%gv)
            allocate(this%xi(0:nG))
            allocate(hkl(3,0:nG))
            hkl(:,0) = 0
            do ii = 1,nG
                hkl(:,ii) =  gethkl(this%gv,ii)
            end do
            inquire(file=trim(xifilename),exist=ok)
            if (.not. ok) call errorExit("Lib_IntegrateManyBeams::readXifile ERROR - file not found """//trim(xifilename)//"""")
            call readExtinctionDistance( xifilename,hkl,getSymmetry(this%latt),this%xi, ok )
            if (.not. ok) call errorExit("Lib_IntegrateManyBeams::readXifile ERROR - error reading file """//trim(xifilename)//"""")

            do ii = 0,nG
                this%xi(ii) = complex( -abs(real(this%xi(ii))) , -abs(aimag(this%xi(ii))) )
            end do

            if (rank==0) print *,""
                
            return
        end subroutine readXifile




        subroutine readInputXyzFile( this,filename,dfg_bar,element )
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      reads in the atom positions ( shared by all processes )
    !*      and establishes the AtomSpace
    !*      returns the average deformation gradient
     
            type(IntegrateManyBeams),intent(inout)              ::      this
            character(len=*),intent(in)                         ::      filename
            real(kind=real64),dimension(3,3),intent(inout)      ::      dfg_bar
            character(len=XYZFILE_ATOMNAMELENGTH),intent(out),optional   ::      element

            type(XYZFile)       ::      xyz
            logical             ::      ok
            real(kind=real64),dimension(3,3)        ::  a_super
            real(kind=real64),dimension(3)          ::  xyz_offset
            
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
        
                        if (rank==0) then
                            write (*,fmt='(a36,a,a36,a,a36)') "average def grad","    ","average strain","    ","average rotation"
                            write (*,fmt='(3f12.6,a,3f12.6,a,3f12.6)') dfg_bar(1,:),"    ",eps_bar(1,:),"    ",rot_bar(1,:)
                            write (*,fmt='(3f12.6,a,3f12.6,a,3f12.6)') dfg_bar(2,:),"    ",eps_bar(2,:),"    ",rot_bar(2,:)
                            write (*,fmt='(3f12.6,a,3f12.6,a,3f12.6)') dfg_bar(3,:),"    ",eps_bar(3,:),"    ",rot_bar(3,:)
                        end if   
                    else        !   don't have columns data
                        print *,"Lib_IntegrateManyBeams::readInputXyzFile info - don't have deformation gradient data in input file"
                        !dfg_bar = reshape( (/1,0,0,0,1,0,0,0,1/),(/3,3/))
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
            call MPI_BCAST(dfg_bar,9,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ii)
#else
            if (.not. ok) call errorExit("readInputXyzFile ERROR - could not read file """//trim(filename)//"""")
#endif  
            xyz_offset = 0.0d0

            this%as = AtomSpace_ctor(this%a,xyz_offset,a_super)
            call setR(this%as,RotationMatrix_identity)
            call setThickness(this%as)
    
            if (rank==0) call report(this%as)
            if (rank==0) print *,""


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
            type(IntegrateManyBeams),intent(inout)          ::      this
            integer,dimension(2,3)          ::      bb
            integer                         ::      nG
            real(kind=real64),dimension(:,:),allocatable        ::      gg          !   (3,nG)      g-vectors
            integer                                 ::      ii

            if (rank==0) print *,"Lib_IntegrateManyBeams::computePhaseFields"

            call getBounds(this%is,bb)
            nG = getn(this%gv)
            allocate(this%x(nG,bb(1,1):bb(2,1),bb(1,2):bb(2,2),bb(1,3):bb(2,3)))
            allocate(this%rho(bb(1,1):bb(2,1),bb(1,2):bb(2,2),bb(1,3):bb(2,3)))        
            allocate(this%grad_x(3,nG,bb(1,1):bb(2,1),bb(1,2):bb(2,2),bb(1,3):bb(2,3)))
            allocate(gg(3,nG))
            if (rank==0) print *,"Lib_IntegrateManyBeams::computePhaseFields info - bb = ",bb," nG = ",nG
            do ii = 1,nG
                gg(:,ii) = getG(this%gv,ii)
            end do
            call computePhaseFactor( this%nAtoms,this%r,gg, this%is,this%x,this%grad_x,this%rho )
            if (rank==0) print *,""
            return
        end subroutine computePhaseFields

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
