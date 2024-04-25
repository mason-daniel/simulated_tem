
    module VoidIsosurfaces
!---^^^^^^^^^^^^^^^^^^^^^^
        use Lib_TopologicallyCorrectClustering
        use Lib_MarchingCubes
        use Lib_SimpleProgressBar        
        use Lib_LinkCell3D
        use Lib_Lattices
        use Lib_SimpleSupercells
        use Lib_DeformationGradients
        ! use Lib_Voxelise
        use iso_fortran_env
        implicit none
        private
        

    !---
    
        public      ::      VoidIsosurface_ctor
        public      ::      delete
        public      ::      report
     
        public      ::      setAtomSublattices 
        public      ::      setDistanceThresh
        !public      ::      voxeliseDistance
        public      ::      setMaxVol
        !public      ::      setDefGrad
        public      ::      getPsi
        public      ::      getDistanceField
        public      ::      getDefGradField
        public      ::      findPhaseField
        public      ::      imposeDistanceThreshold
        public      ::      countVoids
        public      ::      getLC3D
        
    !---
    
        logical,public          ::      VoidIsosurface_dbg = .false.
                   
    !---
    
        type,public     ::      VoidIsosurface
            private 
            type(LinkCell3D)                                ::      lc3d                !   stores the positions of atoms as link-cell list
            type(SimpleSupercell)                           ::      super               !   the periodic cell determining the phase field
            type(Lattice)                                   ::      latt                !   lattice type
            real(kind=real64)                               ::      a0                  !   indicative lengthscale            
            real(kind=real64)                               ::      minVol              !   minimum void size to accept 
            real(kind=real64)                               ::      maxVol              !   maximum void size to accept
            real(kind=real64)                               ::      isolevel            !   minimum void size to accept 
            logical                                         ::      hasDistanceThresh
            real(kind=real64)                               ::      distanceThresh 
            real(kind=real32),dimension(:,:,:),pointer      ::      psi                 !   phase field
            real(kind=real32),dimension(:,:,:,:),pointer    ::      defGrad             !   deformation gradient with resolution of phase field
            real(kind=real32),dimension(:,:,:),pointer      ::      distanceField       !   distance field with resolutoin of phase field            
            integer,dimension(:),pointer                    ::      atomSublattice     
         end type VoidIsosurface
        
    !---
    
        interface VoidIsosurface_ctor
            module procedure    VoidIsosurface_null
            module procedure    VoidIsosurface_ctor0
        end interface
                
        interface delete
            module procedure    delete0
        end interface
        
        interface report
            module procedure    report0
        end interface
        
    contains
!---^^^^^^^^

        function VoidIsosurface_null() result(this)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
            type(VoidIsosurface)           ::      this
            this%lc3d = LinkCell3D_ctor()
            this%super = SimpleSupercell_ctor()
            this%latt = Lattice_ctor()
            this%a0 = 1.0d0
            this%minVol = 0.0d0
            this%maxVol = huge(1.0)
            this%isolevel = 1.05d0
            this%hasDistanceThresh = .false.
            this%distanceThresh = 0.0d0
             
            nullify( this%psi )
            nullify( this%defGrad )
            nullify( this%distanceField )             
            nullify( this%atomSublattice )
            return
        end function VoidIsosurface_null
                         
        function VoidIsosurface_ctor0(x,super,latt,a0,minVol,isolevel) result(this)        
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ 
            real(kind=real64),dimension(:,:),intent(in)         ::      x
            type(SimpleSupercell),intent(in)    ::      super
            type(Lattice),intent(in)            ::      latt
            real(kind=real64),intent(in)        ::      a0,minVol,isolevel
            type(VoidIsosurface)                ::      this
            
            integer                     ::      Nx,Ny,Nz
            integer                     ::      ix,iy,iz , ii  , nAtoms,nMax
            real(kind=real64),dimension(3,3)    ::      superA,cellA
            real                        ::      mm
            
            this = VoidIsosurface_null()
             
            this%super = super
            this%latt = latt
            this%a0 = a0
            this%minVol = minVol
            
            this%isolevel = isolevel
             
            
        !---    estimate maximum atoms per cell - 1% chance of being larger than this number
            nAtoms = size(x,dim=2)
            nMax = 16*getnMotif(this%latt)
            
        !---    construct link cell list    
            superA = getSuperA(this%super)        
            Nx = max(3,floor( norm2(superA(:,1))/(2*this%a0) ))
            Ny = max(3,floor( norm2(superA(:,2))/(2*this%a0) ))
            Nz = max(3,floor( norm2(superA(:,3))/(2*this%a0) ))
            cellA(1:3,1) = superA(1:3,1) / Nx
            cellA(1:3,2) = superA(1:3,2) / Ny
            cellA(1:3,3) = superA(1:3,3) / Nz
            this%lc3d = LinkCell3D_ctor(cellA,nMax,Nx,Ny,Nz)
            do ii = 1,nAtoms
                call add(this%lc3d,ii,x(:,ii))
            end do
            
            
            
            Nx = getNx(this%super) 
            Ny = getNy(this%super) 
            Nz = getNz(this%super) 
            
            print *,"VoidIsosurfaces::VoidIsosurface_ctor0() info - allocating memory for ",Nx*Ny*Nz," cells for deformation gradient and phase field"
            mm = 10*4
            mm = mm*Nx*Ny*Nz/(1024.0*1024.0)
            print *,"VoidIsosurfaces::VoidIsosurface_ctor0() info - memory cost ",mm," Mb"
            
            allocate(this%defGrad(9,0:Nx-1,0:Ny-1,0:Nz-1))                        
            allocate(this%psi(0:Nx-1,0:Ny-1,0:Nz-1))
            do iz = 0,Nz-1
                do iy = 0,Ny-1
                    do ix = 0,Nx-1
                        this%defGrad(1:9,ix,iy,iz) = (/ 1,0,0,0,1,0,0,0,1 /)
                    end do
                end do
            end do      
            
            allocate(this%atomSublattice(getNPoints(this%lc3d)))
            this%atomSublattice = 1
            
            
            return
        end function VoidIsosurface_ctor0
                       
    !---
    
        subroutine delete0(this)
    !---^^^^^^^^^^^^^^^^^^^^^^^^
            type(VoidIsosurface),intent(inout)    ::      this
            if (this%minVol == 0) return                      
            
            if (this%hasDistanceThresh) deallocate(this%distanceField)
            deallocate(this%psi)             
            deallocate(this%defGrad)
            deallocate(this%atomSublattice)
            call delete(this%lc3d)
            this = VoidIsosurface_null()
            return
        end subroutine delete0
        
    !---
    
        subroutine report0(this,u,o)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^
            type(VoidIsosurface),intent(in)    ::      this
            integer,intent(in),optional     ::      u,o
            integer     ::      uu,oo
            uu = 6 ; if (present(u)) uu = u
            oo = 0 ; if (present(o)) oo = o
            write(unit=uu,fmt='(4(a,f12.5))') repeat(" ",oo)//"VoidIsosurface [a0,isolevel,minVol = ",this%a0,",",this%isolevel,",",this%minVol,"]"   
            write(unit=uu,fmt='(a)') repeat(" ",oo+2)//"LinkCell3d list storing atom position"
            call report(this%lc3d,uu,oo+4)
            write(unit=uu,fmt='(a)') repeat(" ",oo+2)//"Supercell for voxel positions"
            call report(this%super,uu,oo+4)         
            write(unit=uu,fmt='(a)') repeat(" ",oo+2)//"Lattice for expected atom positions"
            call report(this%latt,uu,oo+4)         
            return
        end subroutine report0
    
    !---
        
        function getLC3D( this ) result(lc3d)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      return link cell list, often more useful than atom positions alone
    !*      note: returns a shallow copy of lc3d, so copies pointers
            type(VoidIsosurface),intent(inout)              ::      this
            type(LinkCell3d)                                ::      lc3d
            lc3d = this%lc3d
            return
        end function getLC3D
    
    ! !---
    !     subroutine voxeliseDistance(this,x,distancePerAtom)
    ! !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ 
    !         type(VoidIsosurface),intent(inout)              ::      this
    !         real(kind=real64),dimension(:,:),intent(in)     ::      x
    !         real(kind=real64),dimension(:),intent(in)       ::      distancePerAtom
       
    !         real(kind=real64),dimension(3,3)        ::       superA
    !         real(kind=real64)           ::      dd
             
    !         superA = getSuperA(this%super)
            
    !         print *,"VoidIsosurfaces::voxeliseDistance() info - voxelising"
             
    !         call voxeliseAtomicData( x, distancePerAtom, this%a0/2, superA, this%distanceField, fillin=.true. )
            
    !         if (VoidIsosurface_dbg) then
    !             print *,"VoidIsosurfaces::voxeliseDistance() dbg - minmaxavg distance field input  "   &
    !                         ,minval(distancePerAtom),maxval(distancePerAtom),sum(distancePerAtom)/size(distancePerAtom)
    !             dd = 1.0d0/(size(this%distanceField,dim=1)*size(this%distanceField,dim=2)*size(this%distanceField,dim=3))
    !             print *,"VoidIsosurfaces::voxeliseDistance() dbg - minmaxavg distance field output "   &
    !                         ,minval(this%distanceField),maxval(this%distanceField),sum(this%distanceField)*dd
    !         end if
    !         return
    !     end subroutine voxeliseDistance
    
    
        subroutine setDistanceThresh(this,distanceThresh)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ 
            type(VoidIsosurface),intent(inout)              ::      this
            real(kind=real64),intent(in)                    ::      distanceThresh
            
            integer                     ::      Nx,Ny,Nz  
              
            this%hasDistanceThresh = .true.
            this%distanceThresh = distanceThresh
             
            
            Nx = getNx(this%super) 
            Ny = getNy(this%super) 
            Nz = getNz(this%super) 
            
            allocate(this%distanceField(0:Nx-1,0:Ny-1,0:Nz-1))     
            this%distanceField = 0.0d0
            
            return
        end subroutine setDistanceThresh
    
        subroutine setMaxVol(this,maxVol)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^        
            type(VoidIsosurface),intent(inout)                  ::      this
            real(kind=real64),intent(in)                        ::      maxVol
            
            this%maxVol = maxVol
            return
        end subroutine setMaxVol
            
            
                    
    !     subroutine setDefGrad(this,x,defGradPerAtom,singleDefGrad)
    ! !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    ! !*      compute the deformation gradient on the oriented voxel grid, given the def grad on the atoms
    ! !*      if singleDefGrad then can assume there ss only one deformation gradient everywhere.
    !         type(VoidIsosurface),intent(inout)                  ::      this
    !         real(kind=real64),dimension(:,:),intent(in)         ::      x               ! (3,nAtoms)
    !         real(kind=real64),dimension(:,:),intent(in)         ::      defGradPerAtom  ! (9,nAtoms) or (9,1) if singleDefGrad
            
    !         logical,intent(in)                                  ::      singleDefGrad
            
    !         real(kind=real64),dimension(3,3)                    ::      superA,TT,eps,rot
    !         real(kind=real32),dimension(:,:),allocatable        ::      epsOnAtoms
    !         real(kind=real32),dimension(:,:,:,:),allocatable    ::      epsOnVox
    !         integer                     ::      Nx,Ny,Nz , nAtoms
    !         integer                     ::      ii,ix,iy,iz
            
    !         superA = getsuperA(this%super)
            
    !         Nx = getNx(this%super) 
    !         Ny = getNy(this%super) 
    !         Nz = getNz(this%super) 
    !         nAtoms = size(x,dim=2)
    !         !
            
            
    !         if (singleDefGrad) then
    !             print *,"VoidIsosurfaces::setDefGrad() info - deformation gradient taken to be the same everywhere"
    !             TT(1:3,1:3) = reshape( defGradPerAtom(1:9,1),(/3,3/) )
    !             !TT(1:3,1:3) = reshape( (/1,0,0,0,1,0,0,0,1/) , (/3,3/) )
    !             do iz = 0,Nz-1
    !                 do iy = 0,Ny-1
    !                     do ix = 0,Nx-1
    !                         !call progressBar( ix+1 + Nx*(iy + Ny*iz)+1,Nx*Ny*Nz )         
    !                         this%defGrad(1:3,ix,iy,iz) = real( TT(1:3,1),kind=real32 )
    !                         this%defGrad(4:6,ix,iy,iz) = real( TT(1:3,2),kind=real32 )
    !                         this%defGrad(7:9,ix,iy,iz) = real( TT(1:3,3),kind=real32 )
    !                     end do
    !                 end do
    !             end do    
    !             return 
    !         end if
            
            
            
            
            
    !     !---    find the voxelised strain voxelised deformation gradient
    !         print *,"VoidIsosurfaces::setDefGrad() info - voxelising strain elements"
    !         allocate(epsOnVox(6,0:Nx-1,0:Ny-1,0:Nz-1))            
    !         allocate(epsOnAtoms(6,nAtoms))
    !         do ii = 1,nAtoms
    !             TT(1:3,1:3) = reshape( defGradPerAtom(1:9,ii),(/3,3/) )
    !             call DefGradToStrain( TT , eps )
    !             epsOnAtoms(1,ii) = real( eps(1,1),kind=real32 )
    !             epsOnAtoms(2,ii) = real( eps(2,2),kind=real32 )
    !             epsOnAtoms(3,ii) = real( eps(3,3),kind=real32 )
    !             epsOnAtoms(4,ii) = real( eps(1,2),kind=real32 )
    !             epsOnAtoms(5,ii) = real( eps(2,3),kind=real32 )
    !             epsOnAtoms(6,ii) = real( eps(3,1),kind=real32 )
    !         end do
    !         call voxeliseAtomicData( x, epsOnAtoms, this%a0/2, superA, epsOnVox, fillin=.true. )
    !         deallocate(epsOnAtoms)
    !         eps(1,1) = sum( epsOnVox(1,:,:,:) )
    !         eps(2,2) = sum( epsOnVox(2,:,:,:) )
    !         eps(3,3) = sum( epsOnVox(3,:,:,:) )
    !         eps(1,2) = sum( epsOnVox(4,:,:,:) )
    !         eps(2,3) = sum( epsOnVox(5,:,:,:) )
    !         eps(3,1) = sum( epsOnVox(6,:,:,:) )
    !         eps = eps / (Nx*Ny*Nz)
 
            
    !     !---    find the voxelised deformation gradient
    !         print *,"VoidIsosurfaces::setDefGrad() info - voxelising deformation gradient"      
    !         call voxeliseAtomicData( x, defGradPerAtom, this%a0/2, superA, this%defGrad, fillin=.true. )
            
    !     !---    use the voxelised deformation gradient to find rotation matrices 
    !         print *,"VoidIsosurfaces::setDefGrad() info - using interpolated strains and rotation matrices to construct deformation gradients"
    !         do iz = 0,Nz-1
    !             do iy = 0,Ny-1
    !                 do ix = 0,Nx-1
    !                     call progressBar( ix+1 + Nx*(iy + Ny*iz)+1,Nx*Ny*Nz )         
    !                 !---    interpolated strain
                        
    !                     eps(1,1) = epsOnVox(1,ix,iy,iz) 
    !                     eps(2,2) = epsOnVox(2,ix,iy,iz) 
    !                     eps(3,3) = epsOnVox(3,ix,iy,iz) 
    !                     eps(1,2) = epsOnVox(4,ix,iy,iz) 
    !                     eps(2,3) = epsOnVox(5,ix,iy,iz) 
    !                     eps(3,1) = epsOnVox(6,ix,iy,iz) 
    !                     eps(2,1) = eps(1,2)
    !                     eps(3,2) = eps(2,3)
    !                     eps(1,3) = eps(3,1)
    !                 !---    interpolated rot mat
    !                     TT(1:3,1:3) = reshape( this%defGrad(1:9,ix,iy,iz),(/3,3/) )
    !                     call DefGradToRotMat( TT , rot )  
                        
    !                 !---    construct interpolated def grad
    !                     call StrainAndRotMatToDefGrad( eps,rot , TT )
    !                     this%defGrad(1:3,ix,iy,iz) = real( TT(1:3,1),kind=real32 )
    !                     this%defGrad(4:6,ix,iy,iz) = real( TT(1:3,2),kind=real32 )
    !                     this%defGrad(7:9,ix,iy,iz) = real( TT(1:3,3),kind=real32 )
    !                 end do
    !             end do
    !         end do
             
            
            
    !         return
    !     end subroutine setDefGrad
    
    
    !---        
        
        subroutine getPsi(this,psip)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^ 
            type(VoidIsosurface),intent(inout)                      ::      this      
            real(kind=real32),dimension(:,:,:),pointer,intent(out)  ::      psip
            psip => this%psi
            return
        end subroutine getPsi
        
        subroutine getDistanceField(this,dp)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
            type(VoidIsosurface),intent(inout)                      ::      this      
            real(kind=real32),dimension(:,:,:),pointer,intent(out)  ::      dp
            dp => this%distanceField
            return
        end subroutine getDistanceField
    
        subroutine getDefGradField(this,dgp)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
            type(VoidIsosurface),intent(inout)                          ::      this      
            real(kind=real32),dimension(:,:,:,:),pointer,intent(out)    ::      dgp
            dgp => this%defGrad
            return
        end subroutine getDefGradField
    
    !---
    
        subroutine setAtomSublattices(this,s)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      set the most probable sublattice for each atom
            type(VoidIsosurface),intent(inout)                      ::      this    
            integer,dimension(:),intent(in)     ::      s
            integer                             ::      nAtoms,ii
            nAtoms = size(s)
            this%atomSublattice(1:nAtoms) = s(1:nAtoms) 
            if (maxval(s(1:nAtoms))>getnMotif(this%latt)) then
                print *,"VoidIsosurfaces::setAtomSublattices error - motif number in input file does not correspond to expected motif number in lattice"
                stop
            end if
            do ii = 1,getnMotif(this%latt)
                print *,"VoidIsosurfaces::setAtomSublattices info - motif ",ii," count ",count(this%atomSublattice(1:nAtoms)==ii)
            end do
            return
        end subroutine setAtomSublattices
            
    
    
    !---
    
        subroutine findPhaseField(this)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      compute the phase field
            type(VoidIsosurface),intent(inout)                      ::      this    
        
            integer                     ::      Nx,Ny,Nz , ix,iy,iz 
            integer                     ::      mm,nn,nMotif
            integer                     ::      kk,subl
                     
            real(kind=real64),dimension(:,:),allocatable        ::      dx          !   vector from voxel point to atom
             
            integer,dimension(:),allocatable                    ::      id
            
            real(kind=real64),dimension(:,:,:),allocatable      ::      dr0
            real(kind=real64),dimension(3,3)                    ::      TT
            real(kind=real64),dimension(3)                      ::      xx          !   position of voxel point   
            real(kind=real64)                                   ::      rc          !   cutoff range for looking for atoms from voxel point
            
            real(kind=real64),dimension(3)                      ::      dr
            real(kind=real64)                                   ::      dd,psi,psimax,psimaxmin,ia0
            
            
            
            Nx = getNx(this%super) 
            Ny = getNy(this%super) 
            Nz = getNz(this%super) 
            
            nn = getnNeighMax(this%lc3d)           !   max number of atoms near a voxel point            
            allocate(dx(3,nn))
            allocate(id(nn))
             
            
            nn = getNneighbours(this%latt)      !   max number of reference lattice points near a reference lattice point
            nMotif = getnMotif(this%latt)       !   number of sublattices
            
            allocate(dr0(3,nn,nMotif))
            do kk = 1,nMotif
                nn = getNneighbours(this%latt,kk)
                dr0(1:3,1:nn,kk) =  getNeighbours(this%latt,kk)
            end do
            
            
            
            rc = 1.5*this%a0 
            ia0 = 1/this%a0
            print *,"VoidIsosurfaces::findPhaseField info - taking characteristic length a0 = ",this%a0," cutoff ",rc
            
            
                      !      xx = (/ 0.500000 , 0.866025 , 0.816497 /) ! position of atom 1
                      !      dr = realSpaceToCell( this%super,xx )
                      !      ix = floor(dr(1))
                      !      iy = floor(dr(2))
                      !      iz = floor(dr(3))
                      !      print *,"xx = ",xx," cell ",ix,iy,iz
                      !      call neighbourList( this%lc3d,xx,rc, nn,id,dx ) 
                      !   
                      !      
                      !  !---    find the deformation gradient at this point
                      !      TT(1:3,1:3) = reshape( this%defGrad( 1:9 ,ix,iy,iz),(/3,3/) )
                      !      
                      !      print *,"TT = ",TT
                      !      
                      !  !---    find the expected neighbour positions for this atom on its sublattice, and from this the expected normals
                      !      psimaxmin = huge(1.0)
                      !      print *,"nneigh ",nn," id ",id(1:nn)
                      !      !dx = dx*ia0
                      !      do kk = 1,nn
                      !      
                      !      !   for atom k in the neighbour list find the sublattice, and find the shortest normal projection
                      !          subl = this%atomSublattice(id(kk))
                      !          print *,"neigh ",kk," subl ",subl," at ",dx(1:3,kk)*ia0
                      !      end do
                      !      
                      !      do kk = 1,nn
                      !      
                      !      !   for atom k in the neighbour list find the sublattice, and find the shortest normal projection
                      !          subl = this%atomSublattice(id(kk))
                      !          print *,"neigh ",kk," subl ",subl," at ",dx(1:3,kk)
                      !          psimax = 0
                      !          
                      !          do mm = 1,getNneighbours(this%latt,subl)
                      !              dr = dr0(1:3,mm,subl)                                       !   expected position of neighbour
                      !              dr = TT(1:3,1)*dr(1) + TT(1:3,2)*dr(2) + TT(1:3,3)*dr(3)    !   rotated/strained position of neighbour
                      !              psi = dr(1)*dx(1,kk) + dr(2)*dx(2,kk) + dr(3)*dx(3,kk)      !   projection of voxel on expected neighbour location
                      !              !if (psi<=0) cycle
                      !              dd = dr(1)*dr(1) + dr(2)*dr(2) + dr(3)*dr(3)                                                      
                      !              psi = -psi/max(1.0d-8,dd)                                    !   scale so projection = 1/2 at surface of WS cell
                      !              psimax = max(psi,psimax)                                    !   find furthest projection distance from atom to voxel
                      !              write(*,fmt='(a,i4,a,3f12.6,a,3f12.6,4(a,f12.6))') "vec ",mm," dr0 ",dr0(1:3,mm,subl)," T dr0 ",dr(1:3)," dr.dx ",dot_product( dr,dx(:,kk) )," dd ",dd," psi ",psi
                      !          end do
                      !          
                      !          psimaxmin = min(psimaxmin,psimax)                               !   find furthest distance from any atom to voxel
                      !          print *,"psimax psimaxmin ",psimax,psimaxmin*2*ia0
                      !         
                      !      end do
                      !      !stop
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
!            this%psi = real(this%a0,kind=real32)         !   assume that we are miles from any atom, unless proved otherwise.
            this%psi = 2.0                                !   assume that we are miles from any atom, unless proved otherwise.
            do iz = 0,Nz-1
                do iy = 0,Ny-1
                    do ix = 0,Nx-1
                        call progressBar( ix + Nx*(iy + Ny*iz)+1,Nx*Ny*Nz )                    

                    !---    find a position on the phase field grid            
                        xx = cellToRealSpace( this%super,ix,iy,iz )
                        
                    !---    now find the observed neighbour list around the actual point on voxel grid
                        call neighbourList( this%lc3d,xx,rc, nn,id,dx ) 
                        
                        if (nn < 1) cycle       !   no atoms here, phase field must be default 
                        
                        
                    !---    find the deformation gradient at this point
                        TT(1:3,1:3) = reshape( this%defGrad( 1:9 ,ix,iy,iz),(/3,3/) )
                        
                        !if (ix*ix + iy*iy + iz*iz == 0) print *,"0,0,0 ",xx,nn,TT
                        
                    !---    find the expected neighbour positions for this atom on its sublattice, and from this the expected normals
                        psimaxmin = huge(1.0)
                        do kk = 1,nn
                            
                        !   for atom k in the neighbour list find the sublattice, and find the shortest normal projection
                            subl = this%atomSublattice(id(kk))
                            psimax = 0
                            
                            do mm = 1,getNneighbours(this%latt,subl)
                                dr = dr0(1:3,mm,subl)                                       !   expected position of neighbour
                                dr = TT(1:3,1)*dr(1) + TT(1:3,2)*dr(2) + TT(1:3,3)*dr(3)    !   rotated/strained position of neighbour
                                psi = dr(1)*dx(1,kk) + dr(2)*dx(2,kk) + dr(3)*dx(3,kk)      !   projection of voxel on expected neighbour location
                                if (psi>0) cycle
                                dd = dr(1)*dr(1) + dr(2)*dr(2) + dr(3)*dr(3)                                                      
                                psi = -psi/max(1.0d-8,dd)                                    !   scale so projection = 1/2 at surface of WS cell
                                psimax = max(psi,psimax)                                    !   find furthest projection distance from atom to voxel
                                !if (ix*ix + iy*iy + iz*iz == 0) print *,"neigh ",kk,nn," vec ",ll,mm," psi ",psi,psimax
                            end do
                            
                            psimaxmin = min(psimaxmin,psimax)                               !   find furthest distance from any atom to voxel
                        end do
                        
                        !if (ix*ix + iy*iy + iz*iz == 0) print *,"psi ",psimaxmin,psimaxmin*ia0
                        this%psi(ix,iy,iz) = real( psimaxmin*2*ia0 ,kind=real32 )           !   scale correctly with lattice parameter.
                        !this%psi(ix,iy,iz) = real( psimaxmin*2  ,kind=real32 )           !   scale correctly with lattice parameter.
                        
                        
                 !   !---    debugging hcp                        
                 !       if (this%psi(ix,iy,iz)>1.25) then
                 !           print *,"this%psi(ix,iy,iz) = ",this%psi(ix,iy,iz)," at ",ix,iy,iz,xx
                 !           print *,"nneigh ",nn," id ",id(1:nn)
                 !           !dx = dx*ia0
                 !           do kk = 1,nn
                 !           
                 !           !   for atom k in the neighbour list find the sublattice, and find the shortest normal projection
                 !               subl = this%atomSublattice(id(kk))
                 !               print *,"neigh ",kk," subl ",subl," at ",dx(1:3,kk)
                 !           end do
                 !           
                 !           do kk = 1,nn
                 !           
                 !           !   for atom k in the neighbour list find the sublattice, and find the shortest normal projection
                 !               subl = this%atomSublattice(id(kk))
                 !               print *,"neigh ",kk," subl ",subl," at ",dx(1:3,kk)
                 !               psimax = 0
                 !               
                 !               do mm = 1,getNneighbours(this%latt,subl)
                 !                   dr = dr0(1:3,mm,subl)                                       !   expected position of neighbour
                 !                   dr = TT(1:3,1)*dr(1) + TT(1:3,2)*dr(2) + TT(1:3,3)*dr(3)    !   rotated/strained position of neighbour
                 !                   psi = dr(1)*dx(1,kk) + dr(2)*dx(2,kk) + dr(3)*dx(3,kk)      !   projection of voxel on expected neighbour location
                 !                   !if (psi<=0) cycle
                 !                   dd = dr(1)*dr(1) + dr(2)*dr(2) + dr(3)*dr(3)                                                      
                 !                   psi = -psi/max(1.0d-8,dd)                                    !   scale so projection = 1/2 at surface of WS cell
                 !                   psimax = max(psi,psimax)                                    !   find furthest projection distance from atom to voxel
                 !                   write(*,fmt='(a,i4,a,3f12.6,a,3f12.6,4(a,f12.6))') "vec ",mm," dr0 ",dr0(1:3,mm,subl)," T dr0 ",dr(1:3)," dr.dx ",dot_product( dr,dx(:,kk) )," dd ",dd," psi ",psi
                 !               end do
                 !               
                 !               psimaxmin = min(psimaxmin,psimax)                               !   find furthest distance from any atom to voxel
                 !               print *,"psimax psimaxmin ",psimax,psimaxmin*2*ia0
                 !              
                 !           end do
                 !           stop
                 !       end if
                    
                    end do
                end do
            end do
            
            dd = Nx*Ny*Nz
            print *,"VoidIsosurfaces::findPhaseField() dbg - minmaxavg psi ",minval(this%psi),maxval(this%psi),sum(this%psi)/dd
            print *,"VoidIsosurfaces::findPhaseField() dbg - psi>1    ",count(this%psi>1.001),"/",nint(dd)," = ",count(this%psi>1.001)/(0.01*dd),"%"                  
            
            
            return
        end subroutine findPhaseField
        
        subroutine imposeDistanceThreshold(this)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      impose a distance threshold ( if one is set ) on the phase field
            type(VoidIsosurface),intent(inout)                      ::      this    
        
            integer                     ::      nExcluded
            real(kind=real64)           ::      dd
            integer,dimension(0:20)     ::      hist
            real(kind=real64)           ::      dmin,dmax,id
            integer                     ::      ix,iy,iz,jj
            if (.not. this%hasDistanceThresh) return        !   nothing to do
            
                        
            print *,"VoidIsosurfaces::imposeDistanceThreshold() info - applying distance threshold to phase field "            
        
            dd = 0.01*size(this%psi,dim=1)*size(this%psi,dim=2)*size(this%psi,dim=3)
            
        !---    find an indicative histogram of distances
            dmin = minval(this%distanceField)
            dmax = maxval(this%distanceField)    
            id = (size(hist)-1)/(max(1.0d-8,dmax-dmin))
            hist = 0
            do iz = 0,size(this%distanceField,dim=3)-1
                do iy = 0,size(this%distanceField,dim=2)-1
                    do ix = 0,size(this%distanceField,dim=1)-1
                        jj = int( (this%distanceField(ix,iy,iz)-dmin)*id )
                        hist(jj) = hist(jj) + 1
                    end do
                end do
            end do
            print *,"VoidIsosurfaces::imposeDistanceThreshold() info - distance histogram"
            write(*,fmt='(a6,a12,a12,a12)') "bin "," D "," count "," %" 
            do jj = 0,size(hist)-1
                write(*,fmt='(i6,f12.5,i12,f12.3)') jj,dmin + jj*(dmax-dmin)/(size(hist)-1),hist(jj),hist(jj)/dd
            end do
            
        !---    apply distance thresh
            if (this%distanceThresh<0) then
                print *,"VoidIsosurfaces::imposeDistanceThreshold() info - grain interior calc D > ",-this%distanceThresh
                where (this%distanceField > this%distanceThresh)
                    this%psi = min(1.0,this%psi)
                end where
                nExcluded = count(this%distanceField > this%distanceThresh)
                print *,"VoidIsosurfaces::imposeDistanceThreshold() info - excluded voxels ",nExcluded," = ",nExcluded/dd,"%"
            else
                print *,"VoidIsosurfaces::imposeDistanceThreshold() info - grain boundary calc D < ",this%distanceThresh
                where (this%distanceField <-this%distanceThresh)
                    this%psi = min(1.0,this%psi)
                end where                   
                nExcluded = count(this%distanceField <-this%distanceThresh)
                print *,"VoidIsosurfaces::imposeDistanceThreshold() info - excluded voxels ",nExcluded," = ",nExcluded/dd,"%"
            end if                
               
            return
        end subroutine imposeDistanceThreshold
    
    !---
    
    
            
        subroutine countVoids( this,area,volume,nVoids,nVacs , rescaleVacVol)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      given the phase field, compute the area and volume of those regions > isolevel = 1
    !*      by finding the area and volume at this%isolevel, together with derivatives, and extrapolating to isolevel = 1
            type(VoidIsosurface),intent(inout)          ::      this      
            real(kind=real64),intent(out)               ::      area,volume            
            integer,intent(out)                         ::      nVoids,nVacs
            logical,intent(in)                          ::      rescaleVacVol
            
            integer                                     ::      pp  
            real(kind=real64)                           ::      DELTA = 0.0001d0
            real(kind=real64)                           ::      darea,dvolume,omega0,volumePerMonoVacancy,volumeMin,volumeMax!,volumeEst,dV!,d1,d2
            
            integer,dimension(:,:,:),allocatable            ::      indx
            integer                                         ::      nn,kk!,mm
            real(kind=real64),dimension(:,:),allocatable    ::      areap,volumep
            logical                                         ::      ok
            
            real(kind=real64),dimension(:),allocatable      ::      avgVolPerCluster,avgVol2PerCluster
            integer(kind=real64),dimension(:),allocatable   ::      nCluster
            integer                                         ::      nLargestCluster,nAbsLargestCluster,clusterSize
            
            logical,dimension(:,:),allocatable              ::      hasx,hasy,hasz
            integer                                         ::      ix,iy,iz,nNodes
            integer                                         ::      Nx,Ny,Nz
            
            integer                                         ::      alphaExcluded,tinyExcluded
            real(kind=real64),dimension(:,:),allocatable    ::      centreOfmass
            real(kind=real32),dimension(:,:,:),pointer      ::      psi_crop
            
        !---    compute area and volume at range of levels close to this%isolevel          
            print *,"VoidIsosurfaces::countVoids() info - finding volume and area as function of isolevel ",this%isolevel
            
            
            
            Nx = getNx(this%super) 
            Ny = getNy(this%super) 
            Nz = getNz(this%super) 
            
            
        !---    find unique clusters    
            DELTA = (this%isolevel - 1)*0.01d0
            allocate( indx(0:Nx-1,0:Ny-1,0:Nz-1) )
            call topologicallyCorrectCluster( this%psi-this%isolevel+2*DELTA , indx,nn )
            print *,"VoidIsosurfaces::countVoids() info - found ",nn," regions of interest"
            
           ! print *,"indx(1:5,248,7) ",indx(1:5,248,7)
           ! print *,"indx(1:5,249,7) ",indx(1:5,249,7)
           ! print *,"indx(1:5,250,7) ",indx(1:5,250,7)
           ! print *,"indx(1:5,251,7) ",indx(1:5,251,7)
           ! print *,"indx(1:5,252,7) ",indx(1:5,252,7)
            
            
      
        !---    find extent of each cluster      
            print *,"VoidIsosurfaces::countVoids() dbg - allocating memory"
            allocate( hasx(0:Nx-1,0:nn) ) ; hasx = .false.
            allocate( hasy(0:Ny-1,0:nn) ) ; hasy = .false.
            allocate( hasz(0:Nz-1,0:nn) ) ; hasz = .false.
            allocate(centreOfmass(3,0:nn)) ; centreOfMass = 0
            allocate(nCluster(0:nn)) ; nCluster = 0
            do iz = 0,Nz-1
                do iy = 0,Ny-1
                    do ix = 0,Nx-1
                        kk = indx(ix,iy,iz)
                        hasx(ix,kk) = .true.
                        hasy(iy,kk) = .true.
                        hasz(iz,kk) = .true.
                        centreOfMass(:,kk) = centreOfMass(:,kk) + (/ ix,iy,iz/)
                        nCluster(kk) = nCluster(kk) + 1
                    end do
                end do
            end do
            do kk = 0,nn
                centreOfMass(:,kk) = centreOfMass(:,kk) / max(1,nCluster(kk) )
            end do
            
        !---    find the area and volume of each cluster
            print *,"VoidIsosurfaces::countVoids() dbg - find the area and volume of each cluster"
            allocate( areap(-2:2,nn) ) ; areap = 0
            allocate( volumep(-2:2,0:nn) ) ; volumep = 0     
            nullify(psi_crop)
        !   while we do this, we will estimate the vacancy count of the largest cluster
            omega0 = getOmega0(this%latt)*this%a0**3
            nLargestCluster = 0            
            nAbsLargestCluster = 0
            nNodes = Nx*Ny*Nz
             
                      
            do kk = 1,nn
                
            !   crop part of the phase field containing cluster k
                call cropPhaseField( this,this%isolevel-2*DELTA,indx,kk,ok ,hasx(0:,kk),hasy(0:,kk),hasz(0:,kk),psi_crop )
                
                !if ( (kk-309)*(kk-4272) == 0) print *,"crop ",kk,ok,count(psi_crop>this%isolevel) 
                
                             
                if (.not. ok)  cycle        !   crop returns "fail" if the number of lit nodes is too small.
                 
               
            !   make an initial guess of cluster volume based on nodes above threshold
                volume = count(psi_crop>this%isolevel)                               
           
                
             !   write(*,fmt='(a,i8,a,3f12.4,a,i8,a,f10.5,a)') &
             !       "isosurface ",kk," at ",centreOfmass(:,kk)," has ",nCluster(kk)," nodes = ",volume*100.0/(getNx(this%bsuper)*getNy(this%bsuper)*getNz(this%bsuper)),"% "
                
                
            !   find volume and area using marching cubes at a few levels around threshold
                do pp = -2,2
                    call findVolumeAndArea3( this,this%isolevel+pp*DELTA,areap(pp,kk) ,volumep(pp,kk),psi_crop  ,dbg = .false.)        
                end do 
             
            !   extrapolate volume and area at isolevel = 1
            
                call extrapolate( areap(:,kk),volumep(:,kk) ,this%isolevel, DELTA , darea,dvolume )
                areap(0,kk) = darea
                volumep(0,kk) = dvolume
                clusterSize = nint(dvolume/omega0)
                nLargestCluster = max(nLargestCluster,clusterSize)                    
                nAbsLargestCluster = max(nAbsLargestCluster,abs(clusterSize))
                                
            end do                                
            
            deallocate(nCluster)                        
        !   did we find anything at all?
            if (nAbsLargestCluster == 0) then
                this%psi = min( 1.0,this%psi )
                volume = 0
                area = 0
                nVoids = 0
                nVacs = 0
                return        !   nope. 
            end if
            
            
            !stop
            
        !   now get a good first estimate for the size of each cluster
            print *,"VoidIsosurfaces::countVoids() info - largest cluster ",nLargestCluster," assuming omega0 ",omega0
            

            allocate(avgVolPerCluster(0:2*nLargestCluster))
            allocate(avgVol2PerCluster(0:2*nLargestCluster))
            allocate(nCluster(0:2*nLargestCluster))
            avgVolPerCluster = 0
            nCluster = 0
            do kk = 1,nn
                dvolume = volumep(0,kk) 
                if (dvolume<0) cycle
                clusterSize = nint(dvolume/omega0)
                nCluster(clusterSize) = nCluster(clusterSize) + 1
                avgVolPerCluster(clusterSize) = avgVolPerCluster(clusterSize) + dvolume
            end do
             
                                           
        !---    use these sizes to find an average volume per monovacancy  
            volumePerMonoVacancy = omega0
            if (rescaleVacVol) then
                do clusterSize = 1,nLargestCluster
                    if (nCluster(clusterSize)>0) then
                        volumePerMonoVacancy = avgVolPerCluster(clusterSize)/(clusterSize*nCluster(clusterSize))
                        exit
                    end if
                end do
            end if
            
        !---    and use the average volume per monovacancy to decide which voids to cut, and the vacancy count of each            
            nCluster = 0
            avgVolPerCluster = 0
            avgVol2PerCluster = 0
            nVoids = 0
            nVacs = 0
            area = 0
            volume = 0
            nLargestCluster = 0
            volumeMin = volumePerMonoVacancy*this%minVol
            volumeMax = volumePerMonoVacancy*this%maxVol
            do kk = 1,nn
                dvolume = volumep(0,kk) 
                darea = areap(0,kk)
                if ( (darea > 0) .and. (dvolume > volumeMin) .and. (dvolume < volumeMax) ) then
                    nVoids = nVoids + 1
                        
                    clusterSize = nint( dvolume/volumePerMonoVacancy )   
                                        
                    !clusterSize = nint( (dvolume/volumePerMonoVacancy) )     !   assume volume of cluster V ~ V_1 N , not ~ V_1 N^(2/3). This is better for homogeneous defects.
                    nCluster(clusterSize) = nCluster(clusterSize) + 1
                    avgVolPerCluster(clusterSize) = avgVolPerCluster(clusterSize) + dvolume
                    avgVol2PerCluster(clusterSize) = avgVol2PerCluster(clusterSize) + dvolume*dvolume
                    nVacs = nVacs + clusterSize
                    area = area + darea
                    volume = volume + dvolume
                    nLargestCluster = max(nLargestCluster,clusterSize)
                    print *,"void ",nVoids," area,volume ",darea,dvolume," nVacs ",clusterSize
                else if ( (darea < 0) .and. (abs(dvolume) > volumeMin) .and. (abs(dvolume) < volumeMax) ) then
                    print *,"alpha",kk," area,volume ",darea,dvolume," atoms ",nint( abs(dvolume)/volumePerMonoVacancy )                   
                    volumep(0,kk) = 0.0d0       
                    areap(0,kk) = 0.0d0
                else
                    volumep(0,kk) = -1.0d0       
                    areap(0,kk) = 0.0d0
                end if
            end do
            
        !---    output stats based on revised volume per vacancy
            print *,"VoidIsosurfaces::countVoids() info - largest cluster ",nLargestCluster," assuming vac vol ",volumePerMonoVacancy            
            write (*,fmt='(2a8,2a16)') "size","count","<vol>","+/-"
            do clusterSize = 0,nLargestCluster
                if (nCluster(clusterSize)>0) then
                    avgVolPerCluster(clusterSize) = avgVolPerCluster(clusterSize)/nCluster(clusterSize)
                    if (nCluster(clusterSize)>1) then
                        avgVol2PerCluster(clusterSize) = avgVol2PerCluster(clusterSize)/nCluster(clusterSize)
                        avgVol2PerCluster(clusterSize) = sqrt( avgVol2PerCluster(clusterSize) - avgVolPerCluster(clusterSize)*avgVolPerCluster(clusterSize) )
                        avgVol2PerCluster(clusterSize) = avgVol2PerCluster(clusterSize) * sqrt( 1.0d0*nCluster(clusterSize)/(nCluster(clusterSize) - 1.5d0) )
                    else
                        avgVol2PerCluster(clusterSize) = 0
                    end if
                    write (*,fmt='(2i8,2f16.5)') clusterSize,nCluster(clusterSize),avgVolPerCluster(clusterSize),avgVol2PerCluster(clusterSize)                     
                end if
            end do            
            
            
        !---    finally tidy up the phase field
            alphaExcluded = 0
            tinyExcluded  = 0
        
            do iz = 0,Nz-1
                do iy = 0,Ny-1           
                    do ix = 0,Nx-1
                        kk = indx(ix,iy,iz)
                        if (volumep(0,kk) == 0) then
                            this%psi(ix,iy,iz) = min(1.000012345,this%psi(ix,iy,iz))
                            alphaExcluded = alphaExcluded + 1
                        else if (volumep(0,kk) < 0 ) then
                            this%psi(ix,iy,iz) = min(1.000098765,this%psi(ix,iy,iz))
                            tinyExcluded = tinyExcluded + 1
                        end if
                    end do
                end do
            end do
            
            print *,"VoidIsosurfaces::countVoids() dbg - clusters excluded (alpha) ",alphaExcluded
            print *,"VoidIsosurfaces::countVoids() dbg - clusters excluded (small) ",tinyExcluded
             

            deallocate( indx )
            deallocate( hasx )
            deallocate( hasy )
            deallocate( hasz )
            return
        end subroutine countVoids
            
        subroutine extrapolate( areap,volumep ,isolevel, delta , area,volume )
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      use a series of measurements at isolevel +/- delta to extrapolate the area/volume at isolevel = 1
            real(kind=real64),dimension(-2:2),intent(in)        ::  areap,volumep
            real(kind=real64),intent(in)                        ::  isolevel,delta
            real(kind=real64),intent(out)                       ::  area,volume
            real(kind=real64)                       ::      darea,dvolume,d1,d2
            real(kind=real64)                       ::      diso,dvolmax
            
        !---    find the distance we need to extrapolate, and the maximum permitted change in volume            
            diso = 1-isolevel
            dvolmax = abs( 3*volumep(0)*(isolevel**3-1) )
            
            
        !---    extrapolate using central difference            
            d1 = ( volumep(-2) -8*volumep(-1) + 8*volumep(1) - volumep(2) )/(12*delta)
            d2 = (-volumep(-2) +16*volumep(-1) -30*volumep(0) +16*volumep(1) - volumep(2) )/(12*delta*delta)
            dvolume = diso*(d1 + d2*diso/2)
            
           ! if (abs(volumep(0))>1e5) print *,"diso,dvolmax,volumep(0),dvolume ",diso,dvolmax,volumep(0),dvolume
           ! 
           ! print *,"diso*d1     ",diso*d1
           ! print *,"diso^2*d1/2 ",diso*diso*d2/2
           ! 
           ! 
           ! print *,"d1",d1,( 3*volumep(0) - 4*volumep(-1) + volumep(-2) )/(2*delta),(-3*volumep(0) + 4*volumep(1) - volumep(2) )/(2*delta),(-volumep(-1) + volumep(1))/(2*delta)
           ! print *,"d2",d2,(volumep(-2) -2*volumep(-1) + volumep(0) )/(delta*delta), (volumep(0) -2*volumep(1) + volumep(2) )/(delta*delta),(volumep(-1) -2*volumep(0) + volumep(1) )/(delta*delta)
            
            if (isExtrapolationGood( d1,d2,diso,dvolmax )) then
              !  if (abs(volumep(0))>1e5) print *,"full range extrap good v',v""",d1,d2
                dvolume = diso*(d1 + d2*diso/2)           
                d1 = ( areap(-2) -8*areap(-1) + 8*areap(1) - areap(2) )/(12*delta)
                d2 = (-areap(-2) +16*areap(-1) -30*areap(0) +16*areap(1) - areap(2) )/(12*delta*delta)
              !  if (abs(volumep(0))>1e5)print *,"full range extrap good a',a""",d1,d2
                darea = diso*(d1 + d2*diso/2)
                
            else 
            !---    extrapolate using lower isolevels ( closer to 1 )
                d1 = ( 3*volumep(0) - 4*volumep(-1) + volumep(-2) )/(2*delta)
                d2 = (volumep(-2) -2*volumep(-1) + volumep(0) )/(delta*delta)
                if (isExtrapolationGood( d1,d2,diso,dvolmax )) then
               !     if (abs(volumep(0))>1e5)print *,"lower range extrap good v',v""",d1,d2
                    dvolume = diso*(d1 + d2*diso/2)           
                    d1 = ( 3*areap(0) - 4*areap(-1) + areap(-2) )/(2*delta)
                    d2 = (areap(-2) -2*areap(-1) + areap(0) )/(delta*delta)
               !     if (abs(volumep(0))>1e5)print *,"lower range extrap good a',a""",d1,d2
                    darea = diso*(d1 + d2*diso/2)
                else
                !---    extrapolate using upper isolevels ( further from 1 )
                    d1 = (-3*volumep(0) + 4*volumep(1) - volumep(2) )/(2*delta)
                    d2 = (volumep(0) -2*volumep(1) + volumep(2) )/(delta*delta)
                    if (isExtrapolationGood( d1,d2,diso,dvolmax )) then
                !        if (abs(volumep(0))>1e5)print *,"upper range extrap good v',v""",d1,d2
                        dvolume = diso*(d1 + d2*diso/2)           
                        d1 = (-3*areap(0) + 4*areap(1) - areap(2) )/(2*delta)
                        d2 = (areap(0) -2*areap(1) + areap(2) )/(delta*delta)
                !        if (abs(volumep(0))>1e5)print *,"upper range extrap good a',a""",d1,d2
                        darea = diso*(d1 + d2*diso/2)
                    else
                    !---    extrapolate using central isolevels
                        d1 = (-volumep(-1) + volumep(1))/(2*delta)
                        d2 = (volumep(-1) -2*volumep(0) + volumep(1) )/(delta*delta)
                        if (isExtrapolationGood( d1,d2,diso,dvolmax )) then
                 !           if (abs(volumep(0))>1e5)print *,"central range extrap good v',v""",d1,d2
                            dvolume = diso*(d1 + d2*diso/2)           
                            d1 = (-areap(-1) + areap(1))/(2*delta)
                            d2 = (areap(-1) -2*areap(0) + areap(1) )/(delta*delta)
                 !           if (abs(volumep(0))>1e5)print *,"central range extrap good a',a""",d1,d2
                            darea = diso*(d1 + d2*diso/2)
                        else
                           ! if (VIS_DBG) print *,"VoidIsosurfaces::findVolumeAndArea() warning - extrapolation using derivatives failed, using linear extrapolation"
                            
                            !dvolume = volumep(0)*(isolevel**3-1) 
                            !darea = areap(0)*(isolevel**2-1)
                            d1 = (-volumep(-1) + volumep(1))/(2*delta)
                 !           if (abs(volumep(0))>1e5)print *,"linear v'",d1
                            dvolume = d1*diso
                            
                            d2 = volumep(0)*(isolevel**3-1) 
                            if (abs(dvolume)>abs(d2)) dvolume = d2
                            
                            d1 = (-areap(-1) + areap(1))/(2*delta)
                   !         if (abs(volumep(0))>1e5)print *,"linear a'",d1
                            darea = d1*diso
                            d2 = areap(0)*(isolevel**2-1) 
                            if (abs(darea)>abs(d2)) darea = d2
                            
                            
                        end if                        
                    end if
                end if
            end if            
          !  if (abs(volumep(0))>1e5) print *,"vol,area ",volumep(0),areap(0)," dvol,darea ",dvolume,darea
            volume = volumep(0) + dvolume
            area = areap(0) + darea
           
            return
            
        contains
    !---^^^^^^^^
    
            function isExtrapolationGood( d1,d2,diso,dvolmax ) result(is)
        !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
        !*      assume volume = vol0 + diso*(d1 + d2*diso/2) 
        !*      if |diso*(d1 + d2*diso/2)| > dvolmax then the expected difference is too large
        !*      also if |d2*diso| > 10*d1 then expected difference is too large
                real(kind=real64),intent(in)            ::      d1,d2,diso,dvolmax
                logical                                 ::      is
                real(kind=real64)       ::      dvol
                dvol = diso*(d1 + d2*diso/2)
                is = ( abs(dvol) < dvolmax )
                is = is .and. ( abs(d2*diso) < abs(10*d1) )
                return
            end function isExtrapolationGood
            
            
        end subroutine extrapolate
        
        
    
        
        subroutine cropPhaseField( this,isolevel,indx,surf,ok ,hasx,hasy,hasz , psi_crop)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      given the phase field crop those regions > isolevel in the masked region only
    !*      note isolevel /= this%isolevel - we may be doing a numerical derivative
    !*      return ok = false if volume < 3 nodes, in which case also fill in psi = 1  for those nodes. 
    
            type(VoidIsosurface),intent(inout)          ::      this
            real(kind=real64),intent(in)                ::      isolevel
            integer,dimension(0:,0:,0:),intent(in)      ::      indx
            integer,intent(in)                          ::      surf
            logical,intent(out)                         ::      ok
            real(kind=real32),dimension(:,:,:),pointer  ::      psi_crop
            
            logical,dimension(0:),intent(in)            ::      hasx,hasy,hasz
             
            integer             ::      Mx,My,Mz,Mx1,My1,Mz1,Mx2,My2,Mz2 
            integer             ::      ix,iy,iz , jx,jy,jz
                                             
            
            Mx = getNx(this%super)
            My = getNy(this%super)
            Mz = getNz(this%super)

            if (associated(psi_crop)) deallocate(psi_crop)
            
            call findCoverage( hasx,mx1,mx2 )
            call findCoverage( hasy,my1,my2 )
            call findCoverage( hasz,mz1,mz2 )
            
            if ( (mx2<mx1).or.(my2<my1).or.(mz2<mz1) ) then
                ok = .false.
                nullify(psi_crop)
                return
            end if 
            
            ok = .true.
            allocate(psi_crop(0:mx2-mx1,0:my2-my1,0:mz2-mz1))
            
            do iz = 0,mz2-mz1
                jz = mod( iz + mz1,Mz )
                do iy = 0,my2-my1
                    jy = mod( iy + my1,My )
                    do ix = 0,mx2-mx1
                        jx = mod( ix + mx1,Mx )
                        if ( indx(jx,jy,jz)*(indx(jx,jy,jz)-surf)==0 ) then  
                            !   in the cluster we want, or in no cluster at all
                            psi_crop(ix,iy,iz) = this%psi(jx,jy,jz)
                        else
                            !   in a different cluster. Need to set below isolevel in order to not draw a surface here.
                            psi_crop(ix,iy,iz) = min(this%psi(jx,jy,jz),1.0)   
                        end if
                    end do
                end do
            end do
            
            
            if (count(psi_crop>isolevel)<3) then
                !   must be noise. 
                do iz = 0,mz2-mz1
                    jz = mod( iz + mz1,Mz )
                    do iy = 0,my2-my1
                        jy = mod( iy + my1,My )
                        do ix = 0,mx2-mx1
                            jx = mod( ix + mx1,Mx )
                            if ( indx(jx,jy,jz)==surf ) then  
                                !   in the cluster we want, but we have established this cluster is noise. 
                                !   though it doesn't make a difference to the void counting, it is nice to set under threshold, so that the .chgcar doesn't contain it.
                                this%psi(jx,jy,jz) = min(psi_crop(ix,iy,iz),1.0) 
                            end if
                        end do
                    end do
                end do                
                ok = .false.
            end if 
            
            
            
            
            
            return
            
            
        contains
    !---^^^^^^^^
            
            subroutine findCoverage( has,m1,m2 )
        !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
        !*      the array is true from from mod( m1,Mx ) to mod( m2,Mx ) with m2>m1
        !*      m1 and m2 include a buffer of one where necessary
        
        
                logical,dimension(0:),intent(in)        ::      has
                integer,intent(out)                     ::      m1,m2
                
                integer         ::      ii,mm,nn
                
                mm = size(has)
                nn = count(has)
                
                if (nn >= mm-2) then
                    !   lit pretty much all the way across. No advantage in a buffer.
                    m1 = 0 ; m2 = mm-1
                else if (has(0)) then
                    !   lit on the left boundary. Does it stretch across the periodic repeat?
                    if (has(mm-1)) then
                        !   starts midway through cell and crosses boundary
                        do ii = mm-1,1,-1
                            if (has(ii)) then
                                m1 = ii-1 
                            else
                                exit
                            end if
                        end do
                        do ii = m1+1,m1+mm-1                            
                            if (has( mod( ii,mm ) )) then
                                m2 = ii+1 
                            else
                                exit
                            end if
                        end do
                    else
                        !   starts at 0 and ends midway through cell. Because of pbc, have to start on other side.
                        m1 = mm-1
                        do ii = m1+1,m1+mm-1
                            if (has( mod( ii,mm ) )) then
                                m2 = ii+1 
                            else
                                exit
                            end if
                        end do
                    end if
                else if (has(mm-1)) then
                    !   starts midway through cell and ends at mm-1
                    do ii = 0,mm-1
                        if (has(ii)) then
                            m1 = ii-1 ; exit
                        end if
                    end do
                    m2 = mm   
                else if (nn==0) then   
                !   unlit all the way across???
                    m1 = mm-1 ; m2 = 0                                                        
                else 
                    !   starts and ends midway through cell
                    do ii = 1,mm-1
                        if (has(ii)) then
                            m1 = ii-1 ; exit
                        end if
                    end do
                    do ii = mm-2,m1,-1
                        if (has(ii)) then
                            m2 = ii+1 ; exit
                        end if
                    end do                
                end if
                                
                return
            end subroutine findCoverage
                       
        end subroutine cropPhaseField
          
        
        
          
        subroutine findVolumeAndArea3( this,isolevel,area,volume,psi_crop , dbg )
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      given the phase field, compute the area and volume of those regions > isolevel in the masked region only
    !*      note isolevel /= this%isolevel - we may be doing a numerical derivative
    !*      discard any volumes under the minimum
            type(VoidIsosurface),intent(inout)          ::      this
            real(kind=real64),intent(in)                ::      isolevel
            real(kind=real64),intent(out)               ::      area,volume
        
            logical,intent(in)      ::      dbg
            real(kind=real32),dimension(:,:,:),pointer  ::      psi_crop
            integer             ::      ix,iy,iz
                                             
        
            real(kind=real64),dimension(:,:,:),pointer      ::      tris,tris_tmp
            integer                                         ::      nTri 
            real(kind=real64),dimension(8)                  ::      vv
            type(Gridcell)                                  ::      gg
            integer                                         ::      nTriang
            real(kind=real64),dimension(:,:,:),allocatable  ::      triang
            real(kind=real64),dimension(3)                  ::      xx
            integer             ::      jj,kk
            logical             ::      split  
            real(kind=real64),dimension(3,3)        ::      a_cell
            integer             ::      Mx,My,Mz          
            
            a_cell = getA(this%super)
            Mx = size(psi_crop,dim=1)
            My = size(psi_crop,dim=2)
            Mz = size(psi_crop,dim=3)
            
            allocate(triang(3,3,1000))
                
            if (dbg) print *,"Mx,My,Mz ",Mx,My,Mz
            if (dbg) print *,"lbound,ubound",lbound(psi_crop,dim=1),ubound(psi_crop,dim=1)
            if (dbg) then
                do iz = 0,Mz-1
                    print *,"iz = ",iz
                    do iy = 0,My-1
                        write(*,fmt='(100f6.3)') psi_crop(:  ,iy  ,iz  )
                    end do
                end do
            end if    
                         
        !---    marching cubes in portion of full phase field
            allocate(tris(3,3,100)) ; nTri = 0 ; tris = 0
            do iz = 0,Mz-2              !   -2 because the cell counts go 0:M-1, and we always generate a cell 0:1
                do iy = 0,My-2
                    do ix = 0,Mx-2
                        
                        vv(1) = psi_crop(ix  ,iy  ,iz  )
                        vv(2) = psi_crop(ix+1,iy  ,iz  )
                        vv(3) = psi_crop(ix+1,iy+1,iz  )
                        vv(4) = psi_crop(ix  ,iy+1,iz  )
                        vv(5) = psi_crop(ix  ,iy  ,iz+1)
                        vv(6) = psi_crop(ix+1,iy  ,iz+1)
                        vv(7) = psi_crop(ix+1,iy+1,iz+1)
                        vv(8) = psi_crop(ix  ,iy+1,iz+1)
                                            
                        
                        gg = Gridcell_ctor(vv) ; nTriang = 0
                        call Polygonize( gg,isolevel, triang, nTriang,split )                       
                        if (nTriang == 0) cycle
                        
                    !---    make triangles into real space positions
                        do kk = 1,nTriang
                        
                            do jj = 1,3               ! vertex j
                                xx = triang(1:3,jj,kk) + (/ix,iy,iz/)
                                
                                triang(1:3,jj,kk) = a_cell(1:3,1)*xx(1) + a_cell(1:3,2)*xx(2) + a_cell(1:3,3)*xx(3) 
                                
                            end do
                                                                      
                        end do        
                        
                        if (dbg) print *,ix,iy,iz," nTriang ",nTriang,nTri
                        
                    !---    add to collection
                        kk = size(tris,dim=3)
                        if (nTri + nTriang > kk) then
                            allocate(tris_tmp( 3,3,max(kk*2,nTri+nTriang) ))
                            tris_tmp(1:3,1:3,1:kk) = tris(1:3,1:3,1:kk) 
                            deallocate(tris)
                            tris => tris_tmp
                        end if
 
                        tris(1:3,1:3,nTri+1:nTri+nTriang) = triang(1:3,1:3,1:nTriang)
                        nTri = nTri + nTriang
                        
                    end do     
                end do
            end do

            area = 0 ; volume = 0
            !if (count( psi_crop > isolevel) > 1e6 ) print *,"findVolumeAndArea3 very large ",Mx,My,Mz
            !if (count( psi_crop > isolevel) > 1e6 ) print *,"findVolumeAndArea3 very large ",nTri
            
            if (dbg) print *,"<x> ",sum(tris(1,:,1:nTri))/3*nTri,sum(tris(2,:,1:nTri))/3*nTri,sum(tris(3,:,1:nTri))/3*nTri
            
            
            if (nTri > 0) call areaAndVolume2( tris(1:3,1:3,1:nTri),area,volume  )       
            deallocate(tris)                   
            
            
            return
        end subroutine findVolumeAndArea3
        
        
        subroutine areaAndVolume2( tris,a,v )
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
            real(kind=real64),dimension(:,:,:),intent(in)           ::      tris
            real(kind=real64),intent(out)                           ::      a,v
         
            
            real(kind=real64),dimension(3)      ::      x12,x13,vecta
            real(kind=real64)                   ::      dd
            integer     ::      kk,nTri
            
            
            
                 !   debug
               !  integer             ::      nAp,nVp 
               !  real(kind=real64)   ::      Ap,Vp
               !   nAp = 0 ; nVp = 0 ; Ap = 0 ; Vp = 0
                
            nTri = size(tris,dim=3)
            a = 0.0d0 ; v = 0.0d0
            do kk = 1,nTri                
                                    
            !---    find offset from vertex 1
                x12(1:3) = tris(1:3,2,kk) - tris(1:3,1,kk)
                x13(1:3) = tris(1:3,3,kk) - tris(1:3,1,kk)
                
            !---    find cross product
                vecta(1) = x12(2)*x13(3) - x12(3)*x13(2)
                vecta(2) = x12(3)*x13(1) - x12(1)*x13(3)
                vecta(3) = x12(1)*x13(2) - x12(2)*x13(1)
                
            !---    find area
                dd = norm2(vecta) / 2
                a = a + dd
                
              !  if (sum(vecta)>0) then
              !      nAp = nAp + 1
              !      Ap = Ap + dd
              !  end if
                
                
            !---    find volume                        
                dd = ( tris(1,1,kk)*vecta(1) + tris(2,1,kk)*vecta(2) + tris(3,1,kk)*vecta(3) )/6       
                v = v + dd
                
              !  if (dd>0) then
              !      nVp = nVp + 1
              !      Vp = Vp + dd
              !  end if
                
            end do
            
            
            
            !if (nTri > 1e6) print *,"areaAndVolume2 very large ",nTri," a ",nAp,Ap,a," v ",nVp,Vp,v
            return
        end subroutine areaAndVolume2
        
        
        
        
        
        
        
        
        
    end module VoidIsosurfaces
    
    
!!   gfortran -ffree-line-length-256 -Og -g VoidIsosurfaces.f90 -o testVoidIsosurfaces.exe
!
!    
!    program testVoidIsosurfaces
!!---^^^^^^^^^^^^^^^^^^^^^^^^
!        use iso_fortran_env
!        use VoidIsosurfaces
!        implicit none
!        
!        type(VoidIsosurface)           ::      this
!        
!        this = VoidIsosurface_ctor()
!        
!        call report(this)
!        
!        call delete(this)
!        
!        print *,""
!        print *,"done"
!        print *,""
!        
!    end program testVoidIsosurfaces
!    
!        
    
        
        
    