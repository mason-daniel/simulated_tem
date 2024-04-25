
    module simpleMD
!---^^^^^^^^^^^^^^^
  
!*      very simple molecular dynamics  


        use Lib_LinkCell3d          !   add,move and neighbour list functions
        use Lib_SimpleSupercells    !   defines pbc
        use EAM_PotentialTables     !   defines potential function
        use Lib_RandomSeed          !   gaussian variates for Langevin
        use Lib_ConjugateGradient
        use iso_fortran_env

        implicit none
        private

    !---

        real(kind=real64),parameter     ::      KB = 0.00008617d0       !   Boltzmann const eV/K


    !---

        public      ::      MD_ctor
        public      ::      delete
        public      ::      report

    !---

        public      ::      integrateVV,integrateRK4
        public      ::      potentialEnergyArray
        public      ::      electronDensityArray
        public      ::      embeddingEnergyArray
        public      ::      findKineticEnergy
        public      ::      totalEnergy
        public      ::      kickAtom
        public      ::      setDamping
        public      ::      setPositions
        public      ::      wignerSeitz
        public      ::      countPointDefectMovements
        public      ::      countPointDefects

        public      ::      computeEnergyDerivatives
        public      ::      computeSaddlePointEnergy
        public      ::      CGrelax    , CGrelax_old
        public      ::      elasticConstants
        public      ::      dipoleTensor
        
 !       public      ::      tidyShortBonds

    !---

        logical,public          ::      SIMPLEMD_DBG = .false.
        logical,public          ::      SIMPLEMD_KICKRELAX_MOVEALL = .false.
    !---


        type,public     ::      MD
            !private
            type(LinkCell3d)                            ::      lc3d
            real(kind=real64)                           ::      dt              !   timestep (fs)
            real(kind=real64),dimension(:),pointer      ::      mass            !   inverse mass ( eV (A/fs)^-2 )
            real(kind=real64),dimension(:),pointer      ::      imass           !   inverse mass ( eV (A/fs)^-2 )

            type(EAM_Alloy)                             ::      EAM             !   potential
            integer                                     ::      nAtoms
            real(kind=real64),dimension(:,:),pointer    ::      x,v             !   pointer to positions,velocities held externally
            integer,dimension(:),pointer                ::      at              !   atom types
            real(kind=real64),dimension(:,:),pointer    ::      f               !   forces
            real(kind=real64)                           ::      beta,T,sqrt2kTbeta

            integer,dimension(:,:),pointer                  ::      indx
            real(kind=real32),dimension(:,:,:,:),pointer    ::      hess

        end type



    !---

        interface         MD_ctor
            module procedure    MD_null
            module procedure    MD_ctor0
        end interface

        interface         delete
            module procedure    delete0
        end interface

        interface         report
            module procedure    report0
            module procedure    report1
        end interface



    contains
!---^^^^^^^^

        function MD_null() result(this)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
            type(MD)            ::      this
            this%lc3d = LinkCell3d_ctor()
            this%EAM = EAM_Alloy_ctor()
            this%nAtoms = 0
            nullify(this%x)
            nullify(this%f)
            nullify(this%v)
            this%dt = 0.0d0
            this%beta = 0.0d0
            this%T = 0.0d0
            this%sqrt2kTbeta = 0.0d0
            nullify(this%at)
            nullify(this%mass)
            nullify(this%imass)
            nullify(this%indx)
            nullify(this%hess)
            return
        end function MD_null


        function MD_ctor0(x,at,v,a_super,eam,dt) result(this)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
            real(kind=real64),dimension(:,:),pointer        ::      x
            integer,dimension(:),pointer                    ::      at
            real(kind=real64),dimension(:,:),pointer        ::      v
            real(kind=real64),dimension(3,3),intent(in)     ::      a_super
            type(EAM_Alloy),intent(in)                      ::      eam
            real(kind=real64),intent(in)                    ::      dt
            type(MD)                                        ::      this

            type(SimpleSupercell)               ::      super
            real(kind=real64),dimension(3,3)    ::      a_cell
            integer         ::      ii
            integer         ::      nMax,Nx,Ny,Nz

            this = MD_null()


            this%nAtoms = size(x,dim=2)
            this%eam = eam

            this%at => at
            allocate(this%mass(this%nAtoms))
            allocate(this%imass(this%nAtoms))
            do ii = 1,this%nAtoms
                this%mass(ii) = getMass(this%eam,this%at(ii))
                this%imass(ii) = 1/this%mass(ii) 
            end do

            Nx = max(3,floor( norm2(a_super(:,1))/getCutoff(this%eam) ))
            Ny = max(3,floor( norm2(a_super(:,2))/getCutoff(this%eam) ))
            Nz = max(3,floor( norm2(a_super(:,3))/getCutoff(this%eam) ))
            a_cell(:,1) = a_super(:,1) / Nx
            a_cell(:,2) = a_super(:,2) / Ny
            a_cell(:,3) = a_super(:,3) / Nz
            super = SimpleSupercell_ctor(a_cell,Nx,Ny,Nz)

            nMax = estimateMaxCountPerCell(Nx*Ny*Nz,this%nAtoms)

            this%lc3d = LinkCell3D_ctor(nMax,super)


            this%x => x
            this%v => v

            call add(this%lc3d,this%x)



            allocate(this%f(3,this%nAtoms))

            this%dt = dt
            this%beta = 0.0d0
            this%T = 0.0d0
            this%sqrt2kTbeta = 0.0d0

            return
        end function MD_ctor0


    !---


        subroutine delete0(this)
    !---^^^^^^^^^^^^^^^^^^^^^^^^
            type(MD),intent(inout)    ::      this
!            integer     ::      ix,iy,iz
            call delete(this%lc3d)
            call delete(this%EAM)
            if (this%nAtoms>0) then
                deallocate(this%f)
                deallocate(this%mass)
                deallocate(this%imass)
            end if
            if (associated(this%indx)) then
                deallocate(this%indx)
                deallocate(this%hess)
            end if
            this = MD_null()
            return
        end subroutine delete0

    !---

        subroutine report0(this,u,o)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^
            type(MD),intent(in)             ::      this
            integer,intent(in),optional     ::      u,o
            integer                 ::      uu,oo
            real(kind=real64)       ::      pe,ke,Tp,Tk
            uu = 6 ; if (present(u)) uu = u
            oo = 0 ; if (present(o)) oo = o
            call findPotentialEnergy(this,pe)
            call findKineticEnergy(this,ke)
            Tp = pe / (this%nAtoms*1.5*KB)
            Tk = ke / (this%nAtoms*1.5*KB)
            write(unit=uu,fmt='(10(a,f16.6))') repeat(" ",oo)//"MD [PE,KE,total = ",pe,",",ke,",",pe+ke," T(PE),T(KE) ",Tp,",",Tk," ]"

            return
        end subroutine report0


        subroutine report1(this,verbose,u,o)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
            type(MD),intent(in)             ::      this
            integer,intent(in),optional     ::      u,o
            logical,intent(in)              ::      verbose
            integer                 ::      uu,oo
!            real(kind=real64)       ::      pe,ke
            uu = 6 ; if (present(u)) uu = u
            oo = 0 ; if (present(o)) oo = o
            call report(this,uu,oo)
            if (.not. verbose) return
            write(unit=uu,fmt='(10(a,f16.6))') repeat(" ",oo+4)//"dt,beta,T = ",this%dt,",",this%beta,",",this%T
            call report(this%lc3d,uu,oo+4)
            call report(this%EAM,uu,oo+4)
            return
        end subroutine report1

    !---

        subroutine setDamping( this,beta,T )
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
            type(MD),intent(inout)          ::      this
            real(kind=real64),intent(in)    ::      beta,T
            this%beta = beta
            this%T = T
            this%sqrt2kTbeta =  sqrt( 2 * KB * this%T * this%beta / this%dt )
            return
        end subroutine setDamping

        subroutine setPositions(this,x )
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
            type(MD),intent(inout)                      ::      this
            real(kind=real64),dimension(:,:),intent(in) ::      x
            call clear(this%lc3d)
            this%x(1:3,1:this%nAtoms) = x(1:3,1:this%nAtoms)
            call add(this%lc3d,this%x)
            return
        end subroutine setPositions

    !---

        subroutine findForce(this)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      computes the force vector for every atom in the system.
    !*      note: damping only works for atoms of type 1

            type(MD),intent(inout)      ::      this

            integer             ::      ix,iy,iz, ik,nn,ii,jk
            integer             ::      nNeigh,nNeigh_trim
            integer             ::      Nx,Ny,Nz

            real(kind=real64),dimension(:,:),allocatable        ::      dx
            integer,dimension(:),allocatable                    ::      neigh,neigh_trim , indx,at
            real(kind=real64)                                   ::      dd,rc2
            real(kind=real64),dimension(:,:),allocatable        ::      dx_trim
            real(kind=real64),dimension(:),allocatable          ::      dr_trim


!            real(kind=real64)       ::  e1,e2,rc

             real(kind=real64)       ::  dy1,dy2,dy3 , dy01,dy02,dy03



            nn = getnNeighMax(this%lc3d)
            allocate(dx(3,nn))
            allocate(neigh(nn))
            allocate(at(0:nn))
            allocate(indx(getNmax(this%lc3d)))
            allocate(dx_trim(3,nn))
            allocate(dr_trim(nn))
            allocate(neigh_trim(nn))

            Nx = getNx(this%lc3d)
            Ny = getNy(this%lc3d)
            Nz = getNz(this%lc3d)

            !rc = getCutoff(this%eam)
            rc2 = getCutoff(this%eam)*getCutoff(this%eam)

            if (this%beta>0) then
                this%f(1:3,1:this%nAtoms) = reshape( gaussianVariate(3*this%nAtoms), (/3,this%nAtoms/) )
                do ii = 1,this%nAtoms
                    if (this%at(ii) == 1) this%f(1:3,ii) = this%sqrt2kTbeta*this%f(1:3,ii) - this%beta*this%v(1:3,ii)
                end do
            else
                this%f = 0.0d0
            end if

        ! !---    strategy is to loop through each of the atoms, computing each in turn
        !     do ii = 1,this%nAtoms
        !
        !         call neighbourList( this%lc3d,this%x(:,ii), getCutoff(this%eam), nNeigh,neigh,dx )
        !
        !     !---    remove self, and make vectors unit length
        !         nNeigh_trim = 0
        !         do jj = 1,nNeigh
        !             if (neigh(jj)==ii) cycle
        !             dy(1:3) = dx(1:3,jj)
        !             dd = dy(1)*dy(1) + dy(2)*dy(2) + dy(3)*dy(3)
        !             dd = sqrt(dd)
        !             nNeigh_trim = nNeigh_trim + 1
        !             neigh_trim(nNeigh_trim) = neigh(jj)
        !             dr_trim(nNeigh_trim) = dd
        !             dd = 1/dd
        !             dx_trim(1:3,nNeigh_trim) = dy(1:3)*dd
        !         end do
        !
        !         call addForce( this%eam,ii,nNeigh_trim,dr_trim,dx_trim,neigh_trim,this%f )
        !
        !     end do
        !    return


        !---    strategy is to loop through each of the link-cells, computing each in turn
        !   e1 = 0
        !   e2 = 0
            do iz = 0,Nz-1
                do iy = 0,Ny-1
                    do ix = 0,Nx-1

                    !---    find the complete neighbour list- all atoms in all 26+1 nearest cells
                        call neighbourList( this%lc3d,ix,iy,iz, nNeigh,neigh,dx )




                    !---    find the subset of atoms in the neighbour list which are in the central cell, and where they appear in the neighbour list.
                        nn = getNpoints( this%lc3d,ix,iy,iz )
                        do ik = 1,nn
                            ii = getId(this%lc3d,ix,iy,iz,ik)
                            do jk = 1,nNeigh
                                if (neigh(jk) == ii) then
                                    indx(ik) = jk
                                    exit
                                end if
                            end do
                        end do

                    !---    now for each atom in the central cell I can find the subset within cutoff range


                        do ik = 1,nn

                            ii = indx(ik)
                            dy01 = dx(1,ii)
                            dy02 = dx(2,ii)
                            dy03 = dx(3,ii)

                            nNeigh_trim = 0
                            at(0) = this%at(ii)
                            do jk = 1,nNeigh
                                if (ii==jk) cycle

                                dy1 = dx(1,jk) - dy01
                                dd = dy1*dy1
                                dy2 = dx(2,jk) - dy02
                                dd = dd + dy2*dy2
                                dy3 = dx(3,jk) - dy03
                                dd = dd + dy3*dy3

                                if (dd<rc2) then

                                    dd = sqrt(dd)

                                    nNeigh_trim = nNeigh_trim + 1
                                    dr_trim(nNeigh_trim) = dd
                                    neigh_trim(nNeigh_trim) = neigh(jk)

                                    dd = 1/dd
                                    dx_trim(1,nNeigh_trim) = dy1*dd
                                    dx_trim(2,nNeigh_trim) = dy2*dd
                                    dx_trim(3,nNeigh_trim) = dy3*dd

                                    at(nNeigh_trim) = this%at(neigh(jk))

                                end if

                            end do


                            call addForce( this%eam,neigh(ii),nNeigh_trim,at,dr_trim,dx_trim,neigh_trim,this%f )


                           ! if (ii==1) then
                           !     print *,"atom ",ii," neigh ",nNeigh_trim
                           !     print *,dr_trim(1:nNeigh_trim)
                           !     print *,neigh_trim(1:nNeigh_trim)
                           ! end if

                         !  e1 = e1 +     nNeigh_trim

                        end do

                    end do
                end do
            end do

            !print *,"avg neighbours ",e1/this%nAtoms," avg distance ",e2/e1


        !     print *,"force on atom 1 ",this%f(:,1)
        !
        !     do ii = 1,3
        !         call cut(this%lc3d,1,this%x(:,1))
        !         this%x(ii,1) = this%x(ii,1) + 1.0d-6
        !         call add(this%lc3d,1,this%x(:,1))
        !         call findPotentialEnergy( this,e1 )
        !         call cut(this%lc3d,1,this%x(:,1))
        !         this%x(ii,1) = this%x(ii,1) - 2.0d-6
        !         call add(this%lc3d,1,this%x(:,1))
        !         call findPotentialEnergy( this,e2 )
        !         print *,"force on atom 1 ",(e2-e1)/(2*1.0d-6)
        !         call cut(this%lc3d,1,this%x(:,1))
        !         this%x(ii,1) = this%x(ii,1) + 2.0d-6
        !         call add(this%lc3d,1,this%x(:,1))
        !     end do
        !
        !
        !
        !    stop


!             if (this%beta>0) then
!                 do ii = 1,this%nAtoms
!                     dy(1:3) = this%sqrt2kTbeta * gaussianVariate(3)
!                     dy(1:3) = dy(1:3) - this%beta*this%v(1:3,ii)
!                     this%f(1:3,ii) = this%f(1:3,ii) + dy(1:3)
!                 end do
!             end if

            return
        end subroutine findForce



        subroutine findPotentialEnergy(this,e)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      computes the potential energy for every atom in the system.
            type(MD),intent(in)             ::      this
            real(kind=real64),intent(out)   ::      e

            integer             ::      ix,iy,iz, ik,nn,ii,jk
            integer             ::      nNeigh,nNeigh_trim
            integer             ::      Nx,Ny,Nz

            real(kind=real64),dimension(:,:),allocatable        ::      dx
            integer,dimension(:),allocatable                    ::      neigh , indx , at    ! , neigh_trim
            real(kind=real64),dimension(3)                      ::      dy,dy0
            real(kind=real64)                                   ::      dd,rc2
            real(kind=real64),dimension(:),allocatable          ::      dr_trim



            nn = getnNeighMax(this%lc3d)
            allocate(dx(3,nn))
            allocate(neigh(nn))
            allocate(at(0:nn))
            allocate(indx(getNmax(this%lc3d)))
            allocate(dr_trim(nn))
            !allocate(neigh_trim(nn))

            Nx = getNx(this%lc3d)
            Ny = getNy(this%lc3d)
            Nz = getNz(this%lc3d)

            rc2 = getCutoff(this%eam)*getCutoff(this%eam)
            e = 0.0d0


        !---    strategy is to loop through each of the link-cells, computing each in turn
            do iz = 0,Nz-1
                do iy = 0,Ny-1
                    do ix = 0,Nx-1

                    !---    find the complete neighbour list- all atoms in all 26+1 nearest cells
                        call neighbourList( this%lc3d,ix,iy,iz, nNeigh,neigh,dx )

                    !---    find the subset of atoms in the neighbour list which are in the central cell, and where they appear in the neighbour list.
                        nn = getNpoints( this%lc3d,ix,iy,iz )
                        do ik = 1,nn
                            ii = getId(this%lc3d,ix,iy,iz,ik)
                            do jk = 1,nNeigh
                                if (neigh(jk) == ii) then
                                    indx(ik) = jk
                                    exit
                                end if
                            end do
                        end do

                    !---    now for each atom in the central cell I can find the subset within cutoff range
                        do ik = 1,nn
                            dy0(1:3) = -dx(1:3,indx(ik))
                            ii = getId(this%lc3d,ix,iy,iz,ik)
                            at(0) = this%at(ii)
                            nNeigh_trim = 0
                            do jk = 1,nNeigh
                                dy(1:3) = dx(1:3,jk) + dy0(1:3)
                                dd = dy(1)*dy(1) + dy(2)*dy(2) + dy(3)*dy(3)
                                if ((dd-1.0d-8)*(rc2-dd) > 0.0d0) then
                                    dd = sqrt(dd)
                                    nNeigh_trim = nNeigh_trim + 1
                                    dr_trim(nNeigh_trim) = dd
                                    at(nNeigh_trim) = this%at(neigh(jk))
                                    !neigh_trim(nNeigh_trim) = neigh(jk)
                                end if
                            end do
                            !write(*,fmt='(a,i6,a,4i6,a,f16.8,a,i6,a,1000i6)') "atom ",getId(this%lc3d,ix,iy,iz,ik)," at ",ix,iy,iz,ik," e ",getEnergy( this%eam,nNeigh_trim,dr_trim )," neigh ",nNeigh_trim,":",neigh_trim(1:nNeigh_trim)
                            e = e + getEnergy( this%eam,nNeigh_trim,at,dr_trim )
                        end do

                    end do
                end do
            end do
            !print *,"total energy ",e
            !stop



            return
        end subroutine findPotentialEnergy


        subroutine potentialEnergyArray(this,e)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      computes the potential energy for every atom in the system.
            type(MD),intent(in)             ::      this
            real(kind=real64),dimension(:),intent(out)   ::      e

            integer             ::      ix,iy,iz, ik,nn,ii,jk
            integer             ::      nNeigh,nNeigh_trim
            integer             ::      Nx,Ny,Nz

            real(kind=real64),dimension(:,:),allocatable        ::      dx
            integer,dimension(:),allocatable                    ::      neigh , indx , at    ! , neigh_trim
            real(kind=real64),dimension(3)                      ::      dy,dy0
            real(kind=real64)                                   ::      dd,rc2
            real(kind=real64),dimension(:),allocatable          ::      dr_trim



            nn = getnNeighMax(this%lc3d)
            allocate(dx(3,nn))
            allocate(neigh(nn))
            allocate(at(0:nn))
            allocate(indx(getNmax(this%lc3d)))
            allocate(dr_trim(nn))
            !allocate(neigh_trim(nn))

            Nx = getNx(this%lc3d)
            Ny = getNy(this%lc3d)
            Nz = getNz(this%lc3d)

            rc2 = getCutoff(this%eam)*getCutoff(this%eam)
            e = 0.0d0


        !---    strategy is to loop through each of the link-cells, computing each in turn
            do iz = 0,Nz-1
                do iy = 0,Ny-1
                    do ix = 0,Nx-1

                    !---    find the complete neighbour list- all atoms in all 26+1 nearest cells
                        call neighbourList( this%lc3d,ix,iy,iz, nNeigh,neigh,dx )

                    !---    find the subset of atoms in the neighbour list which are in the central cell, and where they appear in the neighbour list.
                        nn = getNpoints( this%lc3d,ix,iy,iz )
                        do ik = 1,nn
                            ii = getId(this%lc3d,ix,iy,iz,ik)
                            do jk = 1,nNeigh
                                if (neigh(jk) == ii) then
                                    indx(ik) = jk
                                    exit
                                end if
                            end do
                        end do

                    !---    now for each atom in the central cell I can find the subset within cutoff range
                        do ik = 1,nn
                            !ii = getId(this%lc3d,ix,iy,iz,ik)
                            dy0(1:3) = -dx(1:3,indx(ik))
                            nNeigh_trim = 0
                            at(0) = this%at( neigh(indx(ik)) )
                            do jk = 1,nNeigh
                                dy(1:3) = dx(1:3,jk) + dy0(1:3)
                                dd = dy(1)*dy(1) + dy(2)*dy(2) + dy(3)*dy(3)
                                if ((dd-1.0d-8)*(rc2-dd) > 0.0d0) then
                                    dd = sqrt(dd)
                                    nNeigh_trim = nNeigh_trim + 1
                                    dr_trim(nNeigh_trim) = dd
                                    at(nNeigh_trim) = this%at( neigh(jk) )
                                    !neigh_trim(nNeigh_trim) = neigh(jk)
                                end if
                            end do
                            !write(*,fmt='(a,i6,a,4i6,a,f16.8,a,i6,a,1000i6)') "atom ",getId(this%lc3d,ix,iy,iz,ik)," at ",ix,iy,iz,ik," e ",getEnergy( this%eam,nNeigh_trim,dr_trim )," neigh ",nNeigh_trim,":",neigh_trim(1:nNeigh_trim)
                            ii = getId(this%lc3d,ix,iy,iz,ik)
                            e(ii) = getEnergy( this%eam,nNeigh_trim,at,dr_trim )
                        end do

                    end do
                end do
            end do
            !print *,"total energy ",e
            !stop



            return
        end subroutine potentialEnergyArray


        subroutine electronDensityArray(this,e)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      computes the electron density energy for every atom in the system.
            type(MD),intent(in)                         ::      this
            real(kind=real64),dimension(:),intent(out)   ::      e

            integer             ::      ix,iy,iz, ik,nn,ii,jk
            integer             ::      nNeigh,nNeigh_trim
            integer             ::      Nx,Ny,Nz

            real(kind=real64),dimension(:,:),allocatable        ::      dx
            integer,dimension(:),allocatable                    ::      neigh , indx , at    ! , neigh_trim
            real(kind=real64),dimension(3)                      ::      dy,dy0
            real(kind=real64)                                   ::      dd,rc2
            real(kind=real64),dimension(:),allocatable          ::      dr_trim



            nn = getnNeighMax(this%lc3d)
            allocate(dx(3,nn))
            allocate(neigh(nn))
            allocate(at(0:nn))
            allocate(indx(getNmax(this%lc3d)))
            allocate(dr_trim(nn))
            !allocate(neigh_trim(nn))

            Nx = getNx(this%lc3d)
            Ny = getNy(this%lc3d)
            Nz = getNz(this%lc3d)

            rc2 = getCutoff(this%eam)*getCutoff(this%eam)
            e = 0.0d0


        !---    strategy is to loop through each of the link-cells, computing each in turn
            do iz = 0,Nz-1
                do iy = 0,Ny-1
                    do ix = 0,Nx-1

                    !---    find the complete neighbour list- all atoms in all 26+1 nearest cells
                        call neighbourList( this%lc3d,ix,iy,iz, nNeigh,neigh,dx )

                    !---    find the subset of atoms in the neighbour list which are in the central cell, and where they appear in the neighbour list.
                        nn = getNpoints( this%lc3d,ix,iy,iz )
                        do ik = 1,nn
                            ii = getId(this%lc3d,ix,iy,iz,ik)
                            do jk = 1,nNeigh
                                if (neigh(jk) == ii) then
                                    indx(ik) = jk
                                    exit
                                end if
                            end do
                        end do

                    !---    now for each atom in the central cell I can find the subset within cutoff range
                        do ik = 1,nn
                            !ii = getId(this%lc3d,ix,iy,iz,ik)
                            dy0(1:3) = -dx(1:3,indx(ik))
                            nNeigh_trim = 0
                            at(0) = this%at( neigh(indx(ik)) )
                            do jk = 1,nNeigh
                                dy(1:3) = dx(1:3,jk) + dy0(1:3)
                                dd = dy(1)*dy(1) + dy(2)*dy(2) + dy(3)*dy(3)
                                if ((dd-1.0d-8)*(rc2-dd) > 0.0d0) then
                                    dd = sqrt(dd)
                                    nNeigh_trim = nNeigh_trim + 1
                                    dr_trim(nNeigh_trim) = dd
                                    at(nNeigh_trim) = this%at( neigh(jk) )
                                    
                                end if
                            end do
                            !write(*,fmt='(a,i6,a,4i6,a,f16.8,a,i6,a,1000i6)') "atom ",getId(this%lc3d,ix,iy,iz,ik)," at ",ix,iy,iz,ik," e ",getEnergy( this%eam,nNeigh_trim,dr_trim )," neigh ",nNeigh_trim,":",neigh_trim(1:nNeigh_trim)
                            ii = getId(this%lc3d,ix,iy,iz,ik)
                            e(ii) = getElectronDensity( this%eam,nNeigh_trim,at,dr_trim )
                        end do

                    end do
                end do
            end do
            !print *,"total energy ",e
            !stop



            return
        end subroutine electronDensityArray
        

        subroutine embeddingEnergyArray(this,e)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      computes the electron density energy for every atom in the system.
            type(MD),intent(in)                             ::      this
            real(kind=real64),dimension(:),intent(out)      ::      e

            integer             ::      ix,iy,iz, ik,nn,ii,jk
            integer             ::      nNeigh,nNeigh_trim
            integer             ::      Nx,Ny,Nz

            real(kind=real64),dimension(:,:),allocatable        ::      dx
            integer,dimension(:),allocatable                    ::      neigh , indx , at    ! , neigh_trim
            real(kind=real64),dimension(3)                      ::      dy,dy0
            real(kind=real64)                                   ::      dd,rc2
            real(kind=real64),dimension(:),allocatable          ::      dr_trim



            nn = getnNeighMax(this%lc3d)
            allocate(dx(3,nn))
            allocate(neigh(nn))
            allocate(at(0:nn))
            allocate(indx(getNmax(this%lc3d)))
            allocate(dr_trim(nn))
            !allocate(neigh_trim(nn))

            Nx = getNx(this%lc3d)
            Ny = getNy(this%lc3d)
            Nz = getNz(this%lc3d)

            rc2 = getCutoff(this%eam)*getCutoff(this%eam)
            e = 0.0d0


        !---    strategy is to loop through each of the link-cells, computing each in turn
            do iz = 0,Nz-1
                do iy = 0,Ny-1
                    do ix = 0,Nx-1

                    !---    find the complete neighbour list- all atoms in all 26+1 nearest cells
                        call neighbourList( this%lc3d,ix,iy,iz, nNeigh,neigh,dx )

                    !---    find the subset of atoms in the neighbour list which are in the central cell, and where they appear in the neighbour list.
                        nn = getNpoints( this%lc3d,ix,iy,iz )
                        do ik = 1,nn
                            ii = getId(this%lc3d,ix,iy,iz,ik)
                            do jk = 1,nNeigh
                                if (neigh(jk) == ii) then
                                    indx(ik) = jk
                                    exit
                                end if
                            end do
                        end do

                    !---    now for each atom in the central cell I can find the subset within cutoff range
                        do ik = 1,nn
                            !ii = getId(this%lc3d,ix,iy,iz,ik)
                            dy0(1:3) = -dx(1:3,indx(ik))
                            nNeigh_trim = 0
                            at(0) = this%at( neigh(indx(ik)) )
                            do jk = 1,nNeigh
                                dy(1:3) = dx(1:3,jk) + dy0(1:3)
                                dd = dy(1)*dy(1) + dy(2)*dy(2) + dy(3)*dy(3)
                                if ((dd-1.0d-8)*(rc2-dd) > 0.0d0) then
                                    dd = sqrt(dd)
                                    nNeigh_trim = nNeigh_trim + 1
                                    dr_trim(nNeigh_trim) = dd
                                    at(nNeigh_trim) = this%at( neigh(jk) )
                                    
                                end if
                            end do
                            !write(*,fmt='(a,i6,a,4i6,a,f16.8,a,i6,a,1000i6)') "atom ",getId(this%lc3d,ix,iy,iz,ik)," at ",ix,iy,iz,ik," e ",getEnergy( this%eam,nNeigh_trim,dr_trim )," neigh ",nNeigh_trim,":",neigh_trim(1:nNeigh_trim)
                            ii = getId(this%lc3d,ix,iy,iz,ik)
                            e(ii) = getEmbeddingEnergy( this%eam,nNeigh_trim,at,dr_trim )
                        end do

                    end do
                end do
            end do
            !print *,"total energy ",e
            !stop



            return
        end subroutine embeddingEnergyArray
        

        subroutine totalEnergy(this,pe,ke)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      computes the potential energy for every atom in the system.
            type(MD),intent(in)             ::      this
            real(kind=real64),intent(out)   ::      pe,ke
            call findPotentialEnergy(this,pe)
            call findKineticEnergy(this,ke)
            return
        end subroutine totalEnergy


        subroutine findKineticEnergy(this,e)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      computes the potential energy for every atom in the system.
            type(MD),intent(in)             ::      this
            real(kind=real64),intent(out)   ::      e

            integer             ::      ii
            real(kind=real64)   ::      vv


            e = 0.0d0
            do ii = 1,this%nAtoms
                vv = this%v(1,ii)*this%v(1,ii) + this%v(2,ii)*this%v(2,ii) + this%v(3,ii)*this%v(3,ii)
                e = e + this%mass(ii)*vv
            end do

            e = e / 2

            return
        end subroutine findKineticEnergy

        subroutine integrateVV( this,initialise )
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      perform a single velocity verlet integration timestep
            type(MD),intent(inout)          ::      this
            logical,intent(in),optional     ::      initialise

!            real(kind=real64),dimension(3)  ::      zz
!            logical                         ::      ok
            integer                         ::      ii
!            integer                 ::      ix,iy,iz,ik,nn
!            type(SimpleSupercell)   ::      super
!            integer                 ::      Nx,Ny,Nz

!            real(kind=real64)   ::      pe,ke



          !  super = getSuper(this%lc3d)
          !  Nx = getNx(super)
          !  Ny = getNy(super)
          !  Nz = getNz(super)

        !---    compute force at first timestep?
             if (present(initialise)) then
                 if (initialise) then
                   call findForce(this)
                 end if
             end if

         !---   compute velocity at half timestep
            do ii = 1,this%nAtoms
                this%v(1:3,ii) = this%v(1:3,ii) + (this%dt*this%imass(ii)/2)*this%f(1:3,ii)
            end do
         !---   compute position at full timestep using estimated velocity
             this%x(1:3,1:this%nAtoms) = this%x(1:3,1:this%nAtoms) + (this%dt)*this%v(1:3,1:this%nAtoms)
         !---   recompute force
             call clear(this%lc3d)
             call add(this%lc3d,this%x)
             call findForce(this)
         !--    compute velocity at full timestep
            do ii = 1,this%nAtoms
                this%v(1:3,ii) = this%v(1:3,ii) + (this%dt*this%imass(ii)/2)*this%f(1:3,ii)
            end do

            
!             call totalEnergy(this,pe,ke)
!             print *,"pe,ke ",pe,ke
            return
        end subroutine integrateVV

        subroutine integrateRK4( this )
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      perform a single RK4 integration timestep
            type(MD),intent(inout)          ::      this

!            real(kind=real64),dimension(3)  ::      zz
!            logical                         ::      ok
            integer                         ::      ii
!            integer                 ::      ix,iy,iz,ik,nn
!            type(SimpleSupercell)   ::      super
!            integer                 ::      Nx,Ny,Nz

            real(kind=real64),dimension(:,:),allocatable,save       ::      x0,v0,x1,v1
            logical,save                                            ::      firstcall = .true.

            if (firstcall) then
                allocate(x0(3,this%nAtoms))
                allocate(x1(3,this%nAtoms))
                allocate(v0(3,this%nAtoms))
                allocate(v1(3,this%nAtoms))
                firstcall = .false.
            end if



        !
        !    super = getSuper(this%lc3d)
        !    Nx = getNx(super)
        !    Ny = getNy(super)
        !    Nz = getNz(super)
        !    mass = 1/this%imass

        !---    store position and velocity at beginning of step- these will be reused to make successive approximations.
            x0(1:3,1:this%nAtoms) = this%x(1:3,1:this%nAtoms)
            v0(1:3,1:this%nAtoms) = this%v(1:3,1:this%nAtoms)



        !---    compute velocities and forces at beginning of step
            call clear(this%lc3d)
            call add(this%lc3d,this%x)
            call findForce(this)
        !---    store the partial position and velocity
            x1(1:3,1:this%nAtoms) = x0(1:3,1:this%nAtoms) + (this%dt/6)           *this%v(1:3,1:this%nAtoms)
            do ii = 1,this%nAtoms
                v1(1:3,ii) = v0(1:3,ii) + (this%imass(ii)*this%dt/6)*this%f(1:3,ii)
            end do



        !---    find the estimated position and velocity at middle of step using input velocity and force
            this%x(1:3,1:this%nAtoms) = x0(1:3,1:this%nAtoms) + (this%dt/2)           *this%v(1:3,1:this%nAtoms)
            do ii = 1,this%nAtoms
                this%v(1:3,ii) = v0(1:3,ii) + (this%imass(ii)*this%dt/2)*this%f(1:3,ii)
            end do
            !this%v(1:3,1:this%nAtoms) = v0(1:3,1:this%nAtoms) + (this%imass*this%dt/2)*this%f(1:3,1:this%nAtoms)
            call clear(this%lc3d)
            call add(this%lc3d,this%x)
            call findForce(this)
        !---    update the partial position and velocity
            x1(1:3,1:this%nAtoms) = x1(1:3,1:this%nAtoms) + (this%dt/3)           *this%v(1:3,1:this%nAtoms)
            !v1(1:3,1:this%nAtoms) = v1(1:3,1:this%nAtoms) + (this%imass*this%dt/3)*this%f(1:3,1:this%nAtoms)
            do ii = 1,this%nAtoms
                v1(1:3,ii) = v1(1:3,ii) + (this%imass(ii)*this%dt/3)*this%f(1:3,ii)
            end do


        !---    recalculate the estimated position and velocity at middle of step using newly estimated velocity and force
            this%x(1:3,1:this%nAtoms) = x0(1:3,1:this%nAtoms) + (this%dt/2)           *this%v(1:3,1:this%nAtoms)
            !this%v(1:3,1:this%nAtoms) = v0(1:3,1:this%nAtoms) + (this%imass*this%dt/2)*this%f(1:3,1:this%nAtoms)
            do ii = 1,this%nAtoms
                this%v(1:3,ii) = v0(1:3,ii) + (this%imass(ii)*this%dt/2)*this%f(1:3,ii)
            end do
            call clear(this%lc3d)
            call add(this%lc3d,this%x)
            call findForce(this)
         !---    update the partial position and velocity
            x1(1:3,1:this%nAtoms) = x1(1:3,1:this%nAtoms) + (this%dt/3)           *this%v(1:3,1:this%nAtoms)
            !v1(1:3,1:this%nAtoms) = v1(1:3,1:this%nAtoms) + (this%imass*this%dt/3)*this%f(1:3,1:this%nAtoms)
            do ii = 1,this%nAtoms
                v1(1:3,ii) = v1(1:3,ii) + (this%imass(ii)*this%dt/3)*this%f(1:3,ii)
            end do



        !---    find the estimated position and velocity at end of step
            this%x(1:3,1:this%nAtoms) = x0(1:3,1:this%nAtoms) + (this%dt)           *this%v(1:3,1:this%nAtoms)
            !this%v(1:3,1:this%nAtoms) = v0(1:3,1:this%nAtoms) + (this%imass*this%dt)*this%f(1:3,1:this%nAtoms)
            do ii = 1,this%nAtoms
                this%v(1:3,ii) = v0(1:3,ii) + (this%imass(ii)*this%dt)*this%f(1:3,ii)
            end do
            call clear(this%lc3d)
            call add(this%lc3d,this%x)
            call findForce(this)

        !---    find the updated position and velocity at the end of the step
            this%x(1:3,1:this%nAtoms) = x1(1:3,1:this%nAtoms) + (this%dt/6)           *this%v(1:3,1:this%nAtoms)
            !this%v(1:3,1:this%nAtoms) = v1(4:6,1:this%nAtoms) + (this%imass*this%dt/6)*this%f(1:3,1:this%nAtoms)
            do ii = 1,this%nAtoms
                this%v(1:3,ii) = v1(1:3,ii) + (this%imass(ii)*this%dt/6)*this%f(1:3,ii)
            end do
            call clear(this%lc3d)
            call add(this%lc3d,this%x)

            return
        end subroutine integrateRK4




        subroutine kickAtom( this,x_old,atom,R0,Ekick,deltax,nTrials , ok )
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ 
    !*      attempt to kick atom with energy Ekick
    !*      atoms initially in relaxed positions 
    !*      compute compensation when kicked atom has travelled a distance deltax
    
            type(MD),intent(inout)          ::      this
            real(kind=real64),dimension(:,:),intent(in)     ::      x_old 
            integer,intent(in)              ::      atom
            real(kind=real64),intent(in)    ::      R0
            real(kind=real64),intent(in)    ::      Ekick
            real(kind=real64),intent(in)    ::      deltax
            integer,intent(in)              ::      nTrials
            logical,intent(out)             ::      ok

            
            
            !integer,parameter               ::      NTRIALS = 8
            

        !---    physically relevant properties
            integer                         ::      nSphere,nCore   !   number of atoms in sphere, number in core ( moving ) region
            real(kind=real64)               ::      rc              !   sphere cutoff radius, atom density
            real(kind=real64)               ::      max_disp2       !   largest displacement (squared) of any atom in sphere
            integer,dimension(:),allocatable,save               ::      indx            !   (1:nAtoms) = 1:nCore where in the core region is an atom
            integer,dimension(:),allocatable,save               ::      id_sphere       !   (1:nSphere) = 1:nAtoms = which atom is which in sphere region
            integer,dimension(:),allocatable,save               ::      at_sphere       !   (1:nSphere) = type of atom
            logical,save                                        ::      firstCall = .true.
            real(kind=real64),dimension(:),allocatable,save     ::      dr2_sphere      !   (1:nNeigh)          distance squared to a neighbour atom
            real(kind=real64),dimension(:,:),allocatable,save   ::      x_sphere       !   (3,nSphere)
            real(kind=real64),dimension(:,:),allocatable,save   ::      x0_sphere       !   (3,nSphere)
!            real(kind=real64),dimension(:),allocatable,save     ::      modx_sphere    !   (NREPLICAS)
            real(kind=real64),dimension(:,:,:),allocatable,save   ::      dx_sphere       !   (3,nSphere,NTRIALS)
            real(kind=real64),dimension(:,:),allocatable,save   ::      md_sphere       !   (3,nSphere)
            real(kind=real64),dimension(:,:),allocatable,save   ::      v_sphere        !   (3,nSphere)
            real(kind=real64),dimension(:),allocatable,save     ::      m_sphere        !   (nSphere)
            real(kind=real64),dimension(:,:),allocatable,save   ::      ff,f0
            real(kind=real32),dimension(:,:),allocatable,save   ::      hess    !   (1:3*nCore+3,1:3*nCore+3) Hessian in core region
            
!            real(kind=real32),dimension(:,:,:),allocatable,save     ::      hessEinstein
            
            real(kind=real64),dimension(:),allocatable                    ::      ee              !   energy of trials

            real(kind=real64),dimension(7),parameter          ::      ALPHA = (/ -0.1d0,0.0d0,0.1d0,0.3d0,0.6d0,1.0d0,2.0d0 /)
            real(kind=real64),dimension(7)                    ::      eeee            !   energy of trial displacements

!            real(kind=real64)               ::      dE_in       !   used as a sanity check 
!            real(kind=real64)               ::      neb_k = 5.0d0   !   nudged elastic band spring constant.
            real(kind=real64)               ::      KE_old,KE_new,PE_old,PE_new
            real(kind=real64)               ::      deltat
!            integer                         ::      saddle_sphere
            integer                         ::      atom_sphere,atom_core
!            integer                         ::      best_displacement
!            real(kind=real64),dimension(:,:),allocatable,save   ::      xhat_sphere
            
        !---    computational workspace
            real(kind=real32),dimension(:),allocatable,save     ::      uu,work !   1:3*nCore+3 RHS of linear equation , LAPACK workspace
            integer,dimension(:),allocatable,save               ::      ipiv
            real(kind=real64),dimension(:),allocatable          ::      dr2     !   (1:nNeigh)          distance squared to a neighbour atom
            real(kind=real64),dimension(:,:),allocatable        ::      dx      !   (1:3,1:nNeigh)      vector separation to neighbour atom
            integer,dimension(:),allocatable                    ::      id      !   (1:nNeigh) =1:nAtoms index of atoms in neighbour list
            integer,dimension(:),allocatable                    ::      at      !   (1:nNeigh) = type of atom

            integer                                             ::      trial

            real(kind=real64),dimension(3)              ::      e_trial!,alpha_trial
            

        !---    general variables
            real(kind=real64)               ::      aa,bb,dd!,x2
            real(kind=real64),dimension(3)  ::      xx  , eta
            integer                         ::      ii,jj,kk, nn
!            real(kind=real32),dimension(3)      ::      h1,h2,h3,gg

!            real(kind=real64),dimension(3,3)    ::      a_cell  !,ia
!            logical,dimension(:),allocatable    ::      hasx

            type(SimpleSupercell)   ::      super

            !  
!             integer                 ::      Nx,Ny,Nz !, ix,iy,iz
!            character(len=256)      ::      dummy
!            character(len=5)        ::      aaaa
!            logical                 ::      ok
!           integer,dimension(3)    ::      nnnn

            
        !---    initialise            
            super = getSuper(this%lc3d)
            rc = getCutoff(this%eam)
            
            
            if (SIMPLEMD_DBG) print *,""
            if (SIMPLEMD_DBG) write(*,fmt='(a,f12.4,a,f12.4,a,i8,a,3f12.4)') "SimpleMD::kickAtom() info - R0 = ",R0," Ekick = ",Ekick," atom = ",atom," at ",x_old(1:3,atom)
            
            
            
        !---    quick escape if R0<0
            if (R0<0) then
                !   am not allowing compensation            
                !   want 1/2 m v^2 = Ekick
                
            !---    select random velocity direction as unit vector on sphere
                do
                    call random_number(eta)
                    eta(1:3) = 2*eta(1:3) - 1
                    dd = eta(1)*eta(1) + eta(2)*eta(2) + eta(3)*eta(3)
                    if (dd*(1-dd)>0) exit
                end do
                
                dd = sqrt( 2 * Ekick * this%imass( atom ) / dd )
                
                this%v(1:3,atom) = dd*eta(1:3)
                ok = .true.
                return
            end if
                
            
            
            
            
            
            
            
            
            
            
            
            
            
        !---    find the number of atoms in the sphere
!             call clear(this%lc3d)
!             call add(this%lc3d,this%x(1:3,1:this%nAtoms))
            !call setPositions(this,x_old(1:3,1:this%nAtoms))
            !call neighbourList_long( this%lc3d,x_old(1:3,atom),R0+2*rc, nSphere )                !   just returns number of atoms in sphere R0+2*rc
            call setPositions(this,this%x(1:3,1:this%nAtoms))
            call neighbourList_long( this%lc3d,this%x(1:3,atom),R0+2*rc, nSphere )                !   just returns number of atoms in sphere R0+2*rc


        !---    step 0: allocate memory to work on the atoms in the local sphere
            if (firstCall) then
                allocate(indx(nSphere))
                allocate(dr2_sphere(nSphere))
                allocate(id_sphere(nSphere))
                allocate(x_sphere(3,nSphere))
                allocate(x0_sphere(3,nSphere))
                allocate(dx_sphere(3,nSphere,NTRIALS))
                allocate(md_sphere(3,nSphere))
                allocate(v_sphere(3,nSphere))
                allocate(m_sphere(nSphere))
                allocate(at_sphere(nSphere))
            else
                if (nSphere > size(indx)) then
                    deallocate(indx)
                    deallocate(dr2_sphere)
                    deallocate(id_sphere)
                    deallocate(x_sphere)
                    deallocate(x0_sphere)
                    deallocate(dx_sphere)
                    deallocate(md_sphere)
                    deallocate(v_sphere)
                    deallocate(m_sphere)
                    deallocate(at_sphere)

                    allocate(indx(nSphere))
                    allocate(dr2_sphere(nSphere))
                    allocate(id_sphere(nSphere))
                    allocate(x_sphere(3,nSphere))
                    allocate(x0_sphere(3,nSphere))
                    allocate(dx_sphere(3,nSphere,NTRIALS))
                    allocate(md_sphere(3,nSphere))
                    allocate(v_sphere(3,nSphere))
                    allocate(m_sphere(nSphere))
                    allocate(at_sphere(nSphere))
                end if
            end if



            !call neighbourList_long( this%lc3d,x_old(1:3,atom),R0+2*rc, nSphere,id_sphere,dr2_sphere )
            call neighbourList_long( this%lc3d,this%x(1:3,atom),R0+2*rc, nSphere,id_sphere,dr2_sphere )
            nCore = count( dr2_sphere(1:nSphere) <= R0*R0 )        !   number which will become core
            if (firstCall) then
                allocate(ff(3,(nCore+1)))
                allocate(f0(3,(nCore+1)))
                allocate(hess(3*(nCore+1),3*(nCore+1)))
                allocate(work(66*3*(nCore+1)))
                allocate(ipiv(3*(nCore+1)))
                allocate(uu(3*(nCore+1)))
                !allocate(hessEinstein(3,3,(nCore+1)))
            else
                if (3*nCore+3 > size(ipiv)) then
                    deallocate(ff)
                    deallocate(f0)
                    deallocate(hess)
                    deallocate(work)
                    deallocate(ipiv)
                    deallocate(uu)
                    !deallocate(hessEinstein)

                    allocate(ff(3,(nCore+1)))
                    allocate(f0(3,(nCore+1)))
                    allocate(hess(3*(nCore+1),3*(nCore+1)))
                    allocate(work(66*3*(nCore+1)))
                    allocate(ipiv(3*(nCore+1)))
                    allocate(uu(3*(nCore+1)))
                    !allocate(hessEinstein(3,3,(nCore+1)))
                end if
            end if
            firstCall = .false.

            
            

        !---    find which atoms are in the core region.
        !       store the types of the atoms in the sphere region
             
            nn = 0
            do kk = 1,nSphere
                ii = id_sphere(kk)
                
                at_sphere(kk) = this%at(ii)
                if ( dr2_sphere(kk) <= R0*R0 ) then
                    nn = nn + 1
                    indx(kk) = nn
                    
                else if ( dr2_sphere(kk) <= (R0+rc)*(R0+rc) ) then
                    indx(kk) = nCore+1                        !   label index = nCore+1 for atoms which are not in the core region, but which directly contribute to force/hessian
                else
                    indx(kk) = -(nCore+1)                     !     label index = -(nCore+1) for atoms which do not directly contribute to force/hessian
                end if
            end do
                
            !print *,"central atom ",atom
            if (SIMPLEMD_DBG) then
                print *,"rCore ",R0," rInner ",R0+rc," rOuter ",R0+2*rc
                print *,"nCore ",nCore," nInner ",count(indx(1:nSphere)==(nCore+1))," nOuter ",count(indx(1:nSphere)<0)," total ",nSphere
            end if
            
            


!        !---    centre the sphere of atoms
!            a_cell = getA(this%lc3d)
!            call inverse3Mat(a_cell,ia)
!            nnnn(1) = getNx(this%lc3d)
!            nnnn(2) = getNy(this%lc3d)
!            nnnn(3) = getNz(this%lc3d)
!            allocate(hasx(0:maxval(nnnn)-1))  
!            offset = 0.0d0 ; ok = .true.
!            do jj = 1,3
!                hasx = .false.
!                do kk = 1,nSphere
!                    ii = id_sphere(kk)
!                    ii = floor( ia(jj,1)*this%x(1,ii) + ia(jj,2)*this%x(2,ii) + ia(jj,3)*this%x(3,ii) )
!                    ii = mod( ii + nnnn(jj)*8 , nnnn(jj) )
!                    hasx(ii) = .true.
!                end do
!               ! write (*,fmt='(100l2)') hasx(0:nnnn(jj)-1)
!                if (hasx(0) .and. hasx(nnnn(jj)-1)) then        !   spans x pbc
!                    do ii = nnnn(jj)-2,1,-1
!                        if (.not. hasx(ii)) then                !   cell ii is empty. Shift back to here
!                            offset(1:3) = offset(1:3) + ii*a_cell(1:3,jj)
!                            ok = .false.
!                            exit
!                        end if
!                    end do
!                end if
!            end do
!            deallocate(hasx)
             


        !---    construct a sphere of atoms to work on in first periodic cell BEFORE kick
        !   x2 = 0.0d0
        !   do kk = 1,nSphere
        !       ii = id_sphere(kk)                      !   = 1:nAtoms
        !       !xx(1:3) = x_old(1:3,ii) - offset(1:3)                 
        !       !x0_sphere(1:3,kk) = wrapPBC( super,xx(1:3) ) 
        !       
        !       !x0_sphere(1:3,kk) = x_old(1:3,ii)
        !       
        !       m_sphere(kk) = this%mass(ii)
        !       v_sphere(1:3,kk) = this%v(1:3,ii)   
        !       if (ii==atom) then
        !           atom_sphere = kk
        !           atom_core = indx(kk)
        !       end if
        !       
        !   end do
        !   if (SIMPLEMD_DBG) print *,"kicked atom ",atom," in sphere ",atom_sphere," in core ",atom_core

        !---    compute the "initial" energy in the core region of the replicas.
        !       also compute initial force  
            nn = getnNeighMax(this%lc3d)                !   maximum number of neighbours in neighbour list
            allocate(dx(3,nn))
            allocate(id(nn))
            allocate(at(0:nn))
            allocate(dr2(nn))
            allocate(ee(nTrials))
            
            
            
            call clear( this%lc3d )
            do kk = 1,nSphere
                ii = id_sphere(kk)                      !   = 1:nAtoms
                x_sphere(1:3,kk) = this%x(1:3,ii)
            end do
            
            call add( this%lc3d,x_sphere(1:3,1:nSphere) )
            f0 = 0.0d0
             

            do kk = 1,nSphere
                if (indx(kk)<0) cycle                   !   only interested in energy in core + inner
                call neighbourList( this%lc3d, x_sphere(1:3,kk) ,rc, nn,id,dx,dr2 )
            !   find self as a neighbour, and remove from neighbour list
                jj = minloc(dr2(1:nn),dim=1)            !   closest neighbour to atom ii is itself.
                id(jj) = id(nn)                         !   replace neighbour list entry jj with the last entry
                dr2(jj) = dr2(nn)
                dx(1:3,jj) = dx(1:3,nn)
                nn = nn - 1

            !   check which core atom type referenced
                at(0) = at_sphere(kk)
                do jj = 1,nn
                    at(jj) = at_sphere( id(jj) )
                    dd = sqrt(dr2(jj))
                    dr2(jj) = dd
                    dx(1:3,jj) = dx(1:3,jj)/dd
                    id(jj) = abs( indx( id(jj) ) )      !   abs because atoms outside R0+rc have indx = -nCore-1          
 
                end do
 
                call addForce( this%eam, indx(kk) ,nn,at,dr2,dx,id, f0 )
                  
            end do
            
         
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            do kk = 1,nSphere
                ii = id_sphere(kk)                      !   = 1:nAtoms
                 x0_sphere(1:3,kk) = x_old(1:3,ii)
                
                m_sphere(kk) = this%mass(ii)
                v_sphere(1:3,kk) = this%v(1:3,ii)   
                if (ii==atom) then
                    atom_sphere = kk
                    atom_core = indx(kk)
                end if
                
            end do
            if (SIMPLEMD_DBG) print *,"kicked atom ",atom," in sphere ",atom_sphere," in core ",atom_core

            
            
            PE_old = 0.0d0
            KE_old = 0.0d0
            call clear( this%lc3d )
            call add( this%lc3d,x0_sphere(1:3,1:nSphere) )
            !f0 = 0.0d0
             

            do kk = 1,nSphere
                if (indx(kk)<0) cycle                   !   only interested in energy in core + inner

                 call neighbourList( this%lc3d, x0_sphere(1:3,kk) ,rc, nn,id,dr2 )
                !call neighbourList( this%lc3d, x0_sphere(1:3,kk) ,rc, nn,id,dx,dr2 )
            !   find self as a neighbour, and remove from neighbour list
                jj = minloc(dr2(1:nn),dim=1)            !   closest neighbour to atom ii is itself.
                id(jj) = id(nn)                         !   replace neighbour list entry jj with the last entry
                dr2(jj) = dr2(nn)
                !dx(1:3,jj) = dx(1:3,nn)
                nn = nn - 1

            !   check which core atom type referenced
                at(0) = at_sphere(kk)
                do jj = 1,nn
                    at(jj) = at_sphere( id(jj) )
                    dd = sqrt(dr2(jj))
                    dr2(jj) = dd
                    !dx(1:3,jj) = dx(1:3,jj)/dd
                    !id(jj) = abs( indx( id(jj) ) )      !   abs because atoms outside R0+rc have indx = -nCore-1          
 
                end do

                PE_old = PE_old + getEnergy(this%eam, nn,at,dr2 )
                !call addForce( this%eam, indx(kk) ,nn,at,dr2,dx,id, f0 )
                 
                KE_old = KE_old + m_sphere(kk)*(v_sphere(1,kk)*v_sphere(1,kk) + v_sphere(2,kk)*v_sphere(2,kk) + v_sphere(3,kk)*v_sphere(3,kk))
            end do
            KE_old = KE_old / 2

            
        !---    compute the offset time          
            deltat = sqrt( this%mass(atom)/( 2*Ekick ) ) * deltax
            if (SIMPLEMD_DBG) print *,"offset displacment, time ",deltax,deltat
            
            dd = 0.0d0
            do kk = 1,nSphere          
                jj = indx(kk) 
                if (abs(jj)>nCore) cycle
                dd = dd + dot_product(f0(1:3,jj),f0(1:3,jj))
            end do
            if (SIMPLEMD_DBG) print *,"initial force <|f|^2> ",dd/nCore," |f(at)|dt^2/(2m) = ",norm2(f0(1:3,atom_core))*deltat*deltat/(2*m_sphere(atom_sphere))
            
            
            
            
        !!---    divide force through by t*t/mass to get displacement units.
        !    do kk = 1,nSphere          
        !        jj = indx(kk) 
        !        if (jj<0) cycle
        !        f0(1:3,jj) = f0(1:3,jj)*deltat*deltat/(2*m_sphere(kk))
        !    end do
            
            
            !print *,"step ",0," replica ",1," energy ",ee(1)

            
            
            dx_sphere = 0.0d0
            ee = huge(1.0)
            do trial = 1,NTRIALS
            
                
            !---    select random velocity direction as unit vector on sphere
                do
                    call random_number(eta)
                    eta(1:3) = 2*eta(1:3) - 1
                    dd = eta(1)*eta(1) + eta(2)*eta(2) + eta(3)*eta(3)
                    if (dd*(1-dd)>0) exit
                end do
                
!!---------------*********    DEBUG FUDGE - vacancy migration in bcc is 1/2[111], right?  ***************
!                select case(trial)
!                    case(1) ; eta = (/ 1, 1, 1/)                   
!                    case(2) ; eta = (/-1, 1, 1/)                   
!                    case(3) ; eta = (/ 1,-1, 1/)                   
!                    case(4) ; eta = (/-1,-1, 1/)                   
!                    case(5) ; eta = (/ 1, 1,-1/)                   
!                    case(6) ; eta = (/-1, 1,-1/)                   
!                    case(7) ; eta = (/ 1,-1,-1/)                   
!                    case(8) ; eta = (/-1,-1,-1/)                   
!                end select 
!                dd = 3
!!---------------******************************************   
                
                eta = eta/sqrt(dd)
                if (SIMPLEMD_DBG) print *,"trial ",trial," direction ",eta," v.eta dt ",dot_product(v_sphere(1:3,atom_sphere),eta)*deltat
                
                x_sphere  = x0_sphere
            !---    place atoms at displaced position.
            !    do kk = 1,nSphere          
            !        !x_sphere(1:3,kk) = x0_sphere(1:3,kk)
            !        jj = indx(kk) 
            !        if (abs(jj)>nCore) cycle
            !        x_sphere(1:3,kk) = x_sphere(1:3,kk) + v_sphere(1:3,kk)*deltat + f0(1:3,jj)*deltat*deltat/(2*m_sphere(kk))     !   note f0 dt^2/(2 m) in length units.
            !    end do
                x_sphere(1:3,atom_sphere) = x_sphere(1:3,atom_sphere) + deltax*eta(1:3) 
                
                
            !       compute energy again displaced but _before compensation_
            !       also compute force and hessian
                
                call clear( this%lc3d )
                call add( this%lc3d,x_sphere(1:3,1:nSphere) )
                PE_new = 0.0d0
                ff = 0.0d0
                hess = 0.0d0
                do kk = 1,nSphere
                    if (indx(kk)<0) cycle                   !   only interested in energy in core + inner
    
                    call neighbourList( this%lc3d, x_sphere(1:3,kk) ,rc, nn,id,dx,dr2 )
    
                !   find self as a neighbour, and remove from neighbour list
                    jj = minloc(dr2(1:nn),dim=1)            !   closest neighbour to atom ii is itself.
                    id(jj) = id(nn)                         !   replace neighbour list entry jj with the last entry
                    dr2(jj) = dr2(nn)
                    dx(1:3,jj) = dx(1:3,nn)
                    nn = nn - 1
    
                !   check which core atom type referenced
                    at(0) = at_sphere(kk)
                    do jj = 1,nn
                        at(jj) = at_sphere( id(jj) )
                        dd = sqrt(dr2(jj))
                        dr2(jj) = dd
                        dx(1:3,jj) = dx(1:3,jj)/dd
                        id(jj) = abs( indx( id(jj) ) )      !   abs because atoms outside R0+rc have indx = -nCore-1                    
                    end do
    
                    PE_new = PE_new + getEnergy(this%eam, nn,at,dr2 )
                    
                    call addForce( this%eam, indx(kk) ,nn,at,dr2,dx,id, ff )
                    call addHessian( this%eam, indx(kk), nn,at,dr2,dx,id, hess )
                    
                     
                end do
                KE_new = KE_old - m_sphere(atom_sphere)*(v_sphere(1,atom_sphere)*v_sphere(1,atom_sphere) + v_sphere(2,atom_sphere)*v_sphere(2,atom_sphere) + v_sphere(3,atom_sphere)*v_sphere(3,atom_sphere))/2 + Ekick
                !print *,"step ",0," replica ",2," energy ",ee(2)
    
    
                if (SIMPLEMD_DBG) print *,"SimpleMD::kickAtom() info - PE before compensation ",PE_old,"->",PE_new," dPE = ",PE_new-PE_old
                if (SIMPLEMD_DBG) print *,"SimpleMD::kickAtom() info - KE before compensation ",KE_old,"->",KE_new," dKE = ",KE_new-KE_old
                
                
                
                
            !---    add a lagrange multiplier to the hessian preventing the displaced atom moving back along eta, and pack force into a 3N vector
                hess( :,3*nCore+1:3*nCore+3 ) = 0.0
               ! hess( 3*nCore+1,: ) = 0.0
                hess( atom_core*3-2:atom_core*3 , 3*nCore+1 ) = real( eta(1:3),kind=real32 )    
               ! hess( 3*nCore+1,atom_core*3-2:atom_core*3 ) = real( eta(1:3),kind=real32 )    
               
               do kk = 1,nSphere
                    jj = abs(indx(kk))              !   where in the core is this sphere atom?
                    if (jj>nCore) cycle             !   only interested in core atoms
                    hess( jj*3-2:jj*3,3*nCore+2 ) = real( f0( 1:3,jj ),kind=real32 )
                end do

               
               
                uu = 0.0d0                
                do kk = 1,nSphere
                    jj = abs(indx(kk))              !   where in the core is this sphere atom?
                    if (jj>nCore) cycle             !   only interested in core atoms
                    uu( jj*3-2:jj*3 ) = real( ff(1:3,jj),kind=real32 )  
                end do
                if (SIMPLEMD_DBG) then
                    dd = dot_product(uu,uu)
                    print *,"SimpleMD::kickAtom() info - |f|^2 ",dd," <|f|^2> ",dd/nCore
                end if
    !                                
                
            !---    solve     
                !call SSYSV( "U",3*nCore+1, 1, hess(1:3*nCore+1,1:3*nCore+1), 3*nCore+1,ipiv(1:3*nCore+1), uu(1:3*nCore+1), 3*nCore+1, work, size(work),ii )       
                call SSYSV( "U",3*nCore+2, 1, hess(1:3*nCore+2,1:3*nCore+2), 3*nCore+2,ipiv(1:3*nCore+2), uu(1:3*nCore+2), 3*nCore+2, work, size(work),ii )       
                if (ii/=0) then
                    !   error in the linear algebra. Not seen one in the wild...
                    print *,"simpleMD::kickAtom error - SSYSV returns ii = ",ii
                    cycle
                end if
                
                
            !---    unpack the solution vector. Check the projection of the solution on eta      
!                dx_sphere = 0.0d0       
                if (SIMPLEMD_DBG) then
                    dd = dot_product( uu( atom_core*3-2:atom_core*3 ),eta )
                    print *,"SimpleMD::kickAtom() info - u.eta = ",dd
                end if
                do kk = 1,nSphere                                                             
                    jj = abs(indx(kk))              !   where in the core is this sphere atom?
                    if (jj>nCore) cycle                                                       
                    dx_sphere( 1:3,kk,trial ) = uu( jj*3-2:jj*3 )                 
                end do
                
                
            !---    carefully move atoms along the solution vector to find optimum displacement.     
                call computeEnergyFnVecLength( this,nSphere,x_sphere,dx_sphere(:,:,trial),at_sphere,ALPHA,eeee,rc/4,indx>0 )
                
            !       where is the lowest energy?
                ii = minloc( eeee,dim=1 )                                   
                if ( (ii-1)*(size(eeee)-ii) == 0 ) then             
                    !   at one or other end
                    dd = ALPHA(ii)
                else
                    !   compute the best ( quadratic approx ) displacement
                    dd = ALPHA(ii-1)*(eeee(ii)-eeee(ii+1)) + ALPHA(ii)*(eeee(ii+1)-eeee(ii-1)) + ALPHA(ii+1)*(eeee(ii-1)-eeee(ii))
                    if (abs(dd)>1.0d-16) dd = 1/(2*dd)
                    dd = ( ALPHA(ii-1)*ALPHA(ii-1)*(eeee(ii)-eeee(ii+1)) + ALPHA(ii)*ALPHA(ii)*(eeee(ii+1)-eeee(ii-1)) + ALPHA(ii+1)*ALPHA(ii+1)*(eeee(ii-1)-eeee(ii)) )*dd
                end if    
                        
                
            
            !---    sanity check with my assumed best displacement dd near minloc ii.
            
                if (SIMPLEMD_DBG) write (*,fmt='(a,2i4,100f16.6)') "SimpleMD::kickAtom() info - ",minloc( eeee,dim=1 ) ,ii,eeee,dd
                call computeEnergyFnVecLength( this,nSphere,x_sphere,dx_sphere(:,:,trial),at_sphere,(/dd/),e_trial,rc/4,indx>0 )   
                !print *,"sanity check 1 - displacement ",dd," energy ",e_trial(1)
                if (e_trial(1) > PE_new) then    
                    !   this is actually quite a likely scenario. If you fire an atom at the neighbours it is going to be horrible.                                     
                    if (SIMPLEMD_DBG) print *,"sanity check failed at ",dd," energy ",e_trial(1)," > ",PE_new
                    cycle
                end if
                PE_new = e_trial(1)
                 
                max_disp2 = 0.0d0
                do kk = 1,nSphere
                    max_disp2 = max( max_disp2,norm2( dx_sphere(1:3,kk,trial) ) )
                end do
                !if (SIMPLEMD_DBG) write (*,fmt='(a,2i4,100f16.6)') "SimpleMD::kickAtom() info - ",minloc( eeee,dim=1 ) ,ii,eeee,dd,max_disp2
                if (SIMPLEMD_DBG) print *,"SimpleMD::kickAtom() info - PE after compensation ",PE_old,"->",PE_new," dPE = ",PE_new-PE_old
                
                ee(trial) = PE_new
                
            end do !    trial
           
            if (SIMPLEMD_DBG) then
                do trial = 1,NTRIALS
                    if (ee(trial)<huge(1.0)/2) write (*,fmt='(i8,f16.4)') trial,ee(trial)-PE_old
                end do     
            end if
            trial = minloc(ee,dim=1)
            PE_new = ee(trial)
            ok = (PE_new < huge(1.0)/2)
            if (.not. ok) return
            
           
            
            
        !---    at this point, I have a vector, dx_sphere, describing the best direction to push the atoms to get them over a barrier
        !       ( assuming there is a barrier there ), starting from CG relaxed positions.
        !       What I haven't accounted for yet is the current MD positions or velocities.
        !       want 1/2 m (v + dv)^2 = 1/2 m v^2 + Ekick
        !       so compute my velocity scaling factor   dd
        !       
            aa = 0.0d0  
            bb = 0.0d0
            dx_sphere(1:3,atom_sphere,trial) = dx_sphere(1:3,atom_sphere,trial) + deltax*eta(1:3)     !   return original kick
            do kk = 1,nSphere
                jj = abs(indx(kk))              !   where in the core is this sphere atom?
                if (jj>nCore) cycle           
                  
                ii = id_sphere(kk)                         
                aa = aa + m_sphere(kk)*( dx_sphere(1,kk,trial)*dx_sphere(1,kk,trial) + dx_sphere(2,kk,trial)*dx_sphere(2,kk,trial) + dx_sphere(3,kk,trial)*dx_sphere(3,kk,trial) )
                bb = bb + m_sphere(kk)*( v_sphere(1,kk)*dx_sphere(1,kk,trial) + v_sphere(2,kk)*dx_sphere(2,kk,trial) + v_sphere(3,kk)*dx_sphere(3,kk,trial) )
            end do
            aa = aa/2
            
            dd = bb*bb + 4*aa*Ekick
             
            dd = sqrt(dd)
                
            if ( (dd - bb)*aa > 0 ) then
                dd = (dd - bb)/(2*aa) 
            else
                dd = (-dd - bb)/(2*aa) 
            end if
            
            if (dd<0) print *,"SimpleMD::kickAtom() warning - kicking atom in the wrong direction! dd = ",dd
            
            !print *,"aa,bb,dd ",aa,bb,dd
            
            
            KE_new = 0.0d0
            do kk = 1,nSphere
                if (indx(kk)<0) cycle                   !   only interested in energy in core + inner
                ii = id_sphere(kk)
                xx(1:3) = v_sphere(1:3,kk) + dd * dx_sphere(1:3,kk,trial)
                this%v(1:3,ii) = xx(1:3)
                KE_new = KE_new + m_sphere(kk)*dot_product(xx,xx)
            end do
            KE_new = KE_new/2
            
            
            if (SIMPLEMD_DBG) print *,"SimpleMD::kickAtom() info - PE after compensation ",PE_old,"->",PE_new," dPE = ",PE_new-PE_old
            if (SIMPLEMD_DBG) print *,"SimpleMD::kickAtom() info - KE after compensation ",KE_old,"->",KE_new," dKE = ",KE_new-KE_old
            
            
            write(*,fmt='(a,f12.4,a,i8,a,3f12.4,a,4f12.4)') "kickAtom ",PE_new-PE_old," atom ",atom,"( "//trim(getName(this%eam,this%at(atom)))//" ) x ",x_old(1:3,atom)," v ",this%v(1:3,atom),this%mass(atom)*dot_product(this%v(1:3,atom),this%v(1:3,atom))/2
            return
        end subroutine kickAtom






        recursive subroutine kickAtom_(this,atom,R0,Ekick,deltax , ok ,attempt)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      take an atom and give it a high kinetic energy in a random direction
    !*      compute the response of the neighbour atoms a short time later,
    !*      and find the best local relaxation mode to accommodate this fast moving atom
    !*      subject to some constraints.
    !*      The result is modifiying the velocity field to include a kicked atom and the local response.
    !*      the total additional energy added to the simulation should be all kinetic, and equal to Ekick.
    !*
    !*      notes: subroutine computes a list of atoms in the local region id_core(1:nCore)
    !*      and a list of atoms in the local region + buffer id_sphere(1:nSphere)
    !*      so I need to (eg) compute force(1:3,1:nCore)
    !*      To save unnecessary conditionals, I refer to a dummy atom outside the core with the index nCore+1
    !*      so will occasionally write to force(1:3,nCore+1) - note this changes alloc ranges
    !*
            type(MD),intent(inout)          ::      this
            integer,intent(in)              ::      atom                !   number of the atom which is to be kicked
            real(kind=real64),intent(in)    ::      R0                  !   radius ( in Angstrom ) of local region which will respond. -ve for no compensation.
            real(kind=real64),intent(in)    ::      Ekick               !   kinetic energy ( eV ) added
            real(kind=real64),intent(in)    ::      deltax              !   displacement at which local relaxation should be computed- note this is immediately turned into a time.
            logical,intent(out)             ::      ok                  !   success? Many factors might make it impossible to construct a good kick and response.
            integer,intent(in),optional     ::      attempt


        !---    physical quantities
            integer                         ::      nCore                   !   number of atoms given velocity total ( = local region )
            integer                         ::      nSphere                 !   number of atoms that need consideration for force calculation ( = local region + buffer )
            real(kind=real64)               ::      tt                      !   time at which local response is computed
            real(kind=real64)               ::      vv                      !   magnitude of the velocity of the kicked atom
            real(kind=real64)               ::      rc                      !   cut-off radius for potential
            real(kind=real64),dimension(3)  ::      eta                 !   random unit vector direction of kicked atom


            real(kind=real64),dimension(:,:),allocatable    ::      f0      !   (1:3,1:nCore+1) force on atoms in core region before the kick
            real(kind=real64),dimension(:,:),allocatable    ::      f1      !   (1:3,1:nCore+1) force on atoms in core region a time t after the kick
            real(kind=real32),dimension(:,:),allocatable    ::      hess    !   (1:3*nCore+3,1:3*nCore+3) Hessian in core region a time t after the kick
            real(kind=real64),dimension(:,:),allocatable    ::      v0      !   (1:3,1:nCore+1) velocity on atoms in core region before the kick
            real(kind=real64),dimension(:,:),allocatable    ::      x0      !   (1:3,1:nCore+1) position of atoms in core region before the kick
            integer                                         ::      atomindx    !   where is the kicked atom in the local region? indx(atom) = atomindx

            real(kind=real64)               ::      e0,e1,e2                   !   potential energy of atoms before and at time t


        !---    computational workspace
            real(kind=real64),dimension(:),allocatable      ::      dr2     !   (1:nNeigh)          distance squared to a neighbour atom
            real(kind=real64),dimension(:,:),allocatable    ::      dx      !   (1:3,1:nNeigh)      vector separation to neighbour atom
            integer,dimension(:),allocatable                ::      id              !   =1:nAtoms index of atoms in neighbour list
            integer,dimension(:),allocatable                ::      id_sphere       !   =1:nAtomsindex of atoms in local+buffer region
            integer,dimension(:),allocatable                ::      indx    !   (1:nAtoms) = 1:nCore where in the core region is an atom?
            real(kind=real64)                               ::      udotu,v0dotv0,udotv0,vpdotu,vpdotv0,udotf0,f0dotf0,f1dotf1,f2dotf2,vpdotvp
            integer,dimension(:),allocatable                ::      at
            real(kind=real32),dimension(:),allocatable      ::      uu,work !   1:3*nCore+3 RHS of linear equation , LAPACK workspace
            integer,dimension(:),allocatable                ::      ipiv

            integer,parameter                               ::      MAX_ATTEMPTS = 8

        !---    dummy
            real(kind=real64)               ::      dd,alpha,beta
            integer                         ::      ii,jj,kk,nn
            real(kind=real64),dimension(3)  ::      yy

            if (SIMPLEMD_DBG) print *,"simpleMD::kickAtom info - considering new kick for atom ",atom
            
            
        !---    find the properties of the kicked atom: velocity magnitude,direction, and time
        !       select random velocity direction as unit vector on sphere
            do
                call random_number(eta)
                eta(1:3) = 2*eta(1:3) - 1
                dd = eta(1)*eta(1) + eta(2)*eta(2) + eta(3)*eta(3)
                if (dd*(1-dd)>0) exit
            end do
            dd = 1/sqrt(dd)
            eta(1:3) = eta(1:3) * dd
        !       determine velocity



        !       quick escape for a simple case: do we have no local region compensation?
            if (R0 <= 0) then
                v0dotv0 = dot_product( this%v(1:3,atom),this%v(1:3,atom) )
                vpdotv0 = dot_product( eta(1:3),this%v(1:3,atom) )
                dd = 2*Ekick*this%imass(atom) + vpdotv0*vpdotv0 - v0dotv0
                if (dd<0) then
                    !   can't quadratic solve 
                    vv = sqrt( 2*Ekick*this%imass(atom) )
                    this%v(1:3,atom) = vv*eta(1:3)
                else
                    dd = sqrt(dd)
                    if (vpdotv0<0) then
                        vv = -dd - vpdotv0
                    else
                        vv = dd - vpdotv0
                    end if
                    this%v(1:3,atom) = this%v(1:3,atom) + vv*eta(1:3)
                end if
                if (SIMPLEMD_DBG) print *,"simpleMD::kickAtom info - no compensation return"
                ok = .true.
                return
            else

                vv = sqrt( 2*Ekick*this%imass(atom) )
            end if

        !       compute time at which compensation is calculated

            tt = deltax/vv
            if (SIMPLEMD_DBG) &
            print *,"simpleMD::kickAtom info - time to reach distance ",deltax," is ",tt





        !---    find the properties of the core region before the atom is kicked
            rc = getCutoff(this%eam)
            call neighbourList_long( this%lc3d,this%x(1:3,atom),R0+2*rc, nSphere )        !   just returns number of atoms in sphere R0+rc
            allocate(dr2(nSphere))
            allocate(id_sphere(nSphere))
            call neighbourList_long( this%lc3d,this%x(1:3,atom),R0+2*rc, nSphere,id_sphere,dr2 )


        !---    count and identify the atoms in the core region. compute "before" positions, forces, velocities.
            nCore = count( dr2(1:nSphere) <= R0*R0 )
            allocate(x0(3,nCore))
            allocate(v0(3,nCore))
            allocate(indx(this%nAtoms))
            indx(:) = nCore+1
            nCore = 0
            do kk = 1,nSphere
                if ( dr2(kk) <= R0*R0 ) then
                    ii = id_sphere(kk)
                    nCore = nCore + 1
                    if (ii==atom) atomindx = nCore
                    indx(ii) = nCore
                    v0(1:3,nCore) = this%v(1:3,ii)
                    x0(1:3,nCore) = this%x(1:3,ii)
                end if
            end do
            deallocate(dr2)

            if (SIMPLEMD_DBG) &
            print *,"simpleMD::kickAtom info - atom kicked is ",atom," at ",x0(1:3,atomindx)



            if (SIMPLEMD_DBG) then
                dd = 0.0d0 ; nn = 0.0
                do kk = 1,nSphere
                    ii = id_sphere(kk)
                    jj = indx(ii)
                    if (jj <= nCore) then
                        alpha = this%v(1,ii)*this%v(1,ii) + this%v(2,ii)*this%v(2,ii) + this%v(3,ii)*this%v(3,ii)
                        alpha = alpha*this%mass(ii)
                        if (alpha>dd) then
                            dd = alpha
                            nn = ii
                        end if
                        !dd = max(dd,alpha)
                    end if
                end do
                dd = dd/2
                print *,"simpleMD::kickAtom info - max KE before kick ",dd," atom ",nn
            end if





        !---    now compute the force and energy before the kick
            nn = getnNeighMax(this%lc3d)                !   maximum number of neighbours in neighbour list
            allocate(dx(3,nn))
            allocate(dr2(nn))
            allocate(id(nn))
            allocate(at(0:nn))
            allocate(f0(3,nCore+1))
            e0 = 0.0d0
            f0 = 0.0d0
            do kk = 1,nSphere
                ii = id_sphere(kk)                      !   = 1:nAtoms , the number of the atom
                call neighbourList( this%lc3d, this%x(1:3,ii) ,rc, nn,id,dx,dr2 )       !   returns index, vector separation and distance squared to all neighbours. Number of neighbours = nn.

            !   find self as a neighbour, and remove from neighbour list
                jj = minloc(dr2(1:nn),dim=1)            !   closest neighbour to atom ii is itself.
                id(jj) = id(nn)                         !   replace neighbour list entry jj with the last entry
                dx(1:3,jj) = dx(1:3,nn)
                dr2(jj) = dr2(nn)
                nn = nn - 1
                at(0) = this%at(ii)

            !   check which core atom is referenced, which neighbour
                do jj = 1,nn
                    at(jj) = this%at(id(jj))
                    id(jj) = indx( id(jj) )             !   either a core atom or to nCore+1

                !   convert square distance and vector to distance and normal vector - this is what my force routine wants as input.
                    dd = sqrt(dr2(jj))
                    dr2(jj) = dd
                    dd = 1/dd
                    dx(1:3,jj) = dx(1:3,jj) * dd
                end do

                call addForce( this%eam, indx(ii) ,nn,at,dr2(1:nn),dx(1:3,1:nn),id(1:nn), f0 )
                e0 = e0 + getEnergy(this%eam, nn,at,dr2(1:nn) )
            end do


            if ( .not. SIMPLEMD_KICKRELAX_MOVEALL ) then
            !---    propagate the central atom only
                call cut(this%lc3d,atom,x0(1:3,atomIndx))
                this%x(1:3,atom) = x0(1:3,atomIndx) + vv*eta(1:3)*tt
                call add(this%lc3d,atom,this%x(1:3,atom))
            else
                do kk = 1,nSphere
                    ii = id_sphere(kk)
                    jj = indx(ii)
                    if (jj <= nCore) then
                    !   only cut and replace atoms in core region
                        call cut(this%lc3d,ii,this%x(1:3,ii) )
                        yy(1:3) = x0(1:3,jj) + v0(1:3,jj)*tt + f0(1:3,jj)*tt*tt*this%imass(ii)/2
                        if (ii == atom) yy(1:3) = yy(1:3) + vv*eta(1:3)*tt
                        this%x(1:3,ii) = yy(1:3)
                        call add( this%lc3d,ii,this%x(1:3,ii) )
                    end if
                end do
            end if

        !---    recompute the energy, forces, and the hessian at this new displaced position.
            allocate(f1(3,nCore+1))
            allocate(hess(3*(nCore+1),3*(nCore+1)))
            f1 = 0.0d0
            hess = 0.0d0
            e1 = 0.0d0
            do kk = 1,nSphere
                ii = id_sphere(kk)                      !   = 1:nAtoms , the number of the atom
                call neighbourList( this%lc3d, this%x(1:3,ii) ,rc, nn,id,dx,dr2 )

            !   find self as a neighbour, and remove from neighbour list
                jj = minloc(dr2(1:nn),dim=1)            !   closest neighbour to atom ii is itself.
                id(jj) = id(nn)                         !   replace neighbour list entry jj with the last entry
                dx(1:3,jj) = dx(1:3,nn)
                dr2(jj) = dr2(nn)
                nn = nn - 1
                at(0) = this%at(ii)

            !   check which core atom is referenced, which neighbour
                do jj = 1,nn
                    at(jj) = this%at(id(jj))
                    id(jj) = indx( id(jj) )             !   either pointing to a core atom or to nCore+1

                !   convert square distance and vector to distance and normal vector
                    dd = sqrt(dr2(jj))
                    dr2(jj) = dd
                    dd = 1/dd
                    dx(1:3,jj) = dx(1:3,jj) * dd
                end do

                call addForce  ( this%eam, indx(ii), nn,at,dr2(1:nn),dx(1:3,1:nn),id(1:nn), f1 )
                call addHessian( this%eam, indx(ii), nn,at,dr2(1:nn),dx(1:3,1:nn),id(1:nn), hess )
                e1 = e1 + getEnergy(this%eam, nn,at,dr2(1:nn) )
            end do

            f0dotf0 = 0.0d0
            f1dotf1 = 0.0d0
            do ii = 1,nCore
                f0dotf0 = f0dotf0 + dot_product( f0(1:3,ii),f0(1:3,ii) )
                f1dotf1 = f1dotf1 + dot_product( f1(1:3,ii),f1(1:3,ii) )
            end do
            if (SIMPLEMD_DBG) &
            print *,"simpleMD::kickAtom info - rms force before, after ",sqrt(f0dotf0)/(3*nCore),sqrt(f1dotf1)/(3*nCore)


!
            if ( .not. SIMPLEMD_KICKRELAX_MOVEALL ) then
            !---    don't forget to restore the central atom
                call cut(this%lc3d,atom,this%x(1:3,atom))
                this%x(1:3,atom) = x0(1:3,atomIndx)
                call add(this%lc3d,atom,this%x(1:3,atom))
            else
                do kk = 1,nSphere
                    ii = id_sphere(kk)
                    jj = indx(ii)
                    if (jj <= nCore) then
                    !   only cut and replace atoms in core region
                        call cut(this%lc3d,ii,this%x(1:3,ii) )
                        this%x(1:3,ii) = x0(1:3,jj)
                        call add( this%lc3d,ii,this%x(1:3,ii) )
                    end if
                end do
            end if


        !---    now we have a chance of a short-cut. If e1>>e0, then this means
        !       that displacing the atoms to time t has a catastrophic effect on the potential energy.
        !       so reject move if too large
            if (SIMPLEMD_DBG) &
            print *,"simpleMD::kickAtom info - potential energy rise due to atom move  ",e1-e0
            if (e1 - e0 > 10*Ekick) then
            !   potential energy rise is 10x Ekick before compensating. This is never going to work!
                ok = .false.
                deallocate(id_sphere)
                deallocate(indx     )
                deallocate(x0       )
                deallocate(v0       )
                deallocate(f0       )
                deallocate(f1       )
                deallocate(hess     )
                deallocate(dx       )
                deallocate(dr2      )
                deallocate(id       )
                if (present(attempt)) then
                    if (attempt==MAX_ATTEMPTS) then
                        if (SIMPLEMD_DBG) print *,"simpleMD::kickAtom info - high energy rise detected too many times. Return failure"
                        return
                    end if
                    call kickAtom_(this,atom,R0,Ekick,deltax , ok ,attempt+1)
                else
                    call kickAtom_(this,atom,R0,Ekick,deltax , ok ,attempt=1)
                end if                
                return
            end if




        !---    we can now solve the constrained minimisation problem.
        !   Add Lagrange multiplier constraints to the hessian.
            allocate(uu(3*nCore+3))                                     !   solution vector = compensating velocities in core
            uu                                        = 0.0d0
            hess(:,3*nCore+1: )                       = 0.0d0           !   set to zero last three columns of hessian as they might have junk in them.
            hess(3*atomindx-2:3*atomindx,3*nCore+1)   = real( eta(1:3) , kind=real32 )      !   set solution vector u.eta = 0
            dd = 1/tt
            do ii = 1,nCore
                hess(ii*3-2:ii*3,3*nCore+2)           = real( f0(1:3,ii) , kind=real32 )    !   set solution vector u.f0 = 0
                uu(ii*3-2:ii*3)                       = real( f1(1:3,ii)*dd , kind=real32 ) !   set RHS = f1/tt
            end do
            hess(3*nCore+3,3*nCore+3) = 1.0d0

        !---    solve minimisation with LAPACK
            allocate(work(66*(3*nCore+3)))
            allocate(ipiv(3*nCore+3))
            call SSYSV( "U",3*nCore+3,1,hess,3*nCore+3,ipiv,uu,3*nCore+3,work,size(work),ii )      !   note: LDA = 3*nCore+3
            deallocate(work)
            deallocate(ipiv)

            if (ii/=0) then
                !   error in the linear algebra. Not seen one in the wild...
                if (SIMPLEMD_DBG) &
                print *,"simpleMD::kickAtom error - SSYSV returns ii = ",ii
                ok = .false.
                deallocate(id_sphere)
                deallocate(indx     )
                deallocate(x0       )
                deallocate(v0       )
                deallocate(f0       )
                deallocate(f1       )
                deallocate(uu       )
                deallocate(hess     )
                deallocate(dx       )
                deallocate(dr2      )
                deallocate(id       )
                return
                 
            end if




!



            uu(3*nCore+1:) = 0.0d0                                 !    unnecessary tidy-up step
            v0dotv0 = 0.0d0
            udotu   = 0.0d0
            udotv0  = 0.0d0
            udotf0  = 0.0d0
            vpdotu  = vv* dot_product( eta(1:3),uu(3*atomindx-2:3*atomindx) )
            vpdotv0 = vv* dot_product( eta(1:3),this%v(1:3,atom) )
            do ii = 1,nCore
                yy(1:3) = uu(3*ii-2:3*ii)
                udotu   = udotu   + dot_product( yy(1:3),yy(1:3) )
                udotv0  = udotv0  + dot_product( yy(1:3),v0(1:3,ii) )
                v0dotv0 = v0dotv0 + dot_product( v0(1:3,ii),v0(1:3,ii) )
                udotf0  = udotf0  + dot_product( yy(1:3),f0(1:3,ii) )
            end do

            if (SIMPLEMD_DBG) then
                print *,"simpleMD::kickAtom info - 1/2 m v0.v0  = ",v0dotv0 /(2*this%imass(atomindx))
                print *,"simpleMD::kickAtom info - 1/2 m u .u   = ",udotu /(2*this%imass(atomindx))
                print *,"simpleMD::kickAtom info - 1/2 m u .v0  = ",udotv0/(2*this%imass(atomindx))
                print *,"simpleMD::kickAtom info - 1/2 m v*.u   = ",vpdotu/(2*this%imass(atomindx))
                print *,"simpleMD::kickAtom info - 1/2 m v*.v0  = ",vpdotv0/(2*this%imass(atomindx))
                print *,"simpleMD::kickAtom info - 1/2 m v*.v*  = ",vv*vv/(2*this%imass(atomindx))
                print *,"simpleMD::kickAtom info -       u .f0  = ",udotf0


            end if

        !   add the central kick atom to the vector u so I'm only dealing with a single vector.
            uu(3*atomindx-2:3*atomindx) = real( uu(3*atomindx-2:3*atomindx) + vv*eta(1:3) , kind=real32 )

        !   this updates the dot products slightly
            udotu  = udotu  + 2*vpdotu + vv*vv                       !  note that vpdotu should be zero.
            udotv0 = udotv0 + vpdotv0

        !   now the energy rise is  1/2 m u.(u+2v0) > Ekick
        !   rescale u = alpha u + beta v0
            if (v0dotv0 == 0) then
                !   zero temperature case is special, and easy.
                beta = 0.0d0        !   (unnecessary)
                alpha = sqrt( 2*Ekick*this%imass(atom) / udotu )
                if (SIMPLEMD_DBG) &
                print *,"simpleMD::kickAtom info - changing KE of kicked atom, alpha,beta = ",alpha,beta
            else
                !   finite temperature case: may be possible to just reduce KE in local region




                dd = udotv0*udotv0 + v0dotv0*( 2*Ekick*this%imass(atom) - udotu + v0dotv0 )

                !print *,"simpleMD::kickAtom info udotv0,v0dotv0, 2*Ekick*this%imass,udotu,dd ",  udotv0,v0dotv0, 2*Ekick*this%imass,udotu,dd
                if (dd >= 0) then
                    !   yes, can compensate by reducing energy in local region
                    alpha = ( -sqrt( dd ) - udotv0 )/v0dotv0
                    beta = ( sqrt( dd ) - udotv0  )/v0dotv0
                    if ((alpha<1).and.(alpha > beta)) beta = alpha
                    alpha = 1

                    if (SIMPLEMD_DBG) &
                    print *,"simpleMD::kickAtom info - rescaling KE in local region, alpha,beta = ",alpha,beta

                    if ( beta*(1-beta)<0 ) then
                        beta = 0.0d0
                        alpha = sqrt( ( 2*Ekick*this%imass(atom) + v0dotv0 )/udotu )
                        if (SIMPLEMD_DBG) &
                        print *,"simpleMD::kickAtom info - minimising KE in local region, alpha,beta = ",alpha,beta
                    end if
                else
                    !   no, can't compensate completely
                    beta = 0.0d0
                    alpha = sqrt( ( 2*Ekick*this%imass(atom) + v0dotv0 )/udotu )

                !    dd    = (udotu*v0dotv0 - udotv0*udotv0)
                !    alpha = sqrt( ( 2*Ekick*this%imass + v0dotv0 )*v0dotv0 / dd )
                !    beta = -1.0d0 - alpha*udotv0/v0dotv0
                    if (SIMPLEMD_DBG) &
                    print *,"simpleMD::kickAtom info - minimising KE in local region, alpha,beta = ",alpha,beta
                end if
            end if

            if ((SIMPLEMD_DBG).and. (max( abs(alpha),abs(beta) )>1)) print *,"simpleMD::kickAtom warning - KE rescaling, alpha,beta = ",alpha,beta



            do ii = 1,nCore
                uu(3*ii-2:3*ii) = real( alpha*uu(3*ii-2:3*ii) + beta*v0(1:3,ii),kind=real32 )
            end do

            if (SIMPLEMD_DBG) then
                udotu   = 0.0d0
                do ii = 1,nCore
                    yy(1:3) = uu(3*ii-2:3*ii)
                    udotu   = udotu   + dot_product( yy(1:3),yy(1:3) )
                end do
                print *,"simpleMD::kickAtom info - 1/2 m u'.u'  = ",udotu /(2*this%imass(atom))
            end if

            if (SIMPLEMD_DBG) then
        !---    now the energy injected in by this KE boost should be 1/2 m u^2 - 1/2 m v0^2 = Ekick
        !       I actually want it to be Ekick.
        !       compute dd = u.( u + 2 v0 )
        !       check the above actually works...
                dd = 0.0d0
                do ii = 1,nCore
                    yy(1:3) = uu(3*ii-2:3*ii)
                    dd = dd + yy(1)*yy(1)
                    dd = dd + yy(2)*yy(2)
                    dd = dd + yy(3)*yy(3)
                    yy(1:3) = v0(1:3,ii)
                    dd = dd - yy(1)*yy(1)
                    dd = dd - yy(2)*yy(2)
                    dd = dd - yy(3)*yy(3)
                end do
                print *,"simpleMD::kickAtom info - kinetic energy rise due to compensation ",dd/(2*this%imass(atom))
           end if



           if (SIMPLEMD_DBG) then
        !---    the potential energy at time t should be lower if the compensation is added.


            !---    move atoms to new displacements including compensation
                do kk = 1,nSphere
                    ii = id_sphere(kk)
                    jj = indx(ii)
                    if (jj <= nCore) then
                    !   only cut and replace atoms in core region
                        call cut(this%lc3d,ii,this%x(1:3,ii) )
                        yy(1:3) = ( uu(3*jj-2:3*jj) - v0(1:3,jj) )*tt            !   delta v dt
                        this%x(1:3,ii) = x0(1:3,jj) + yy(1:3)
                        call add( this%lc3d,ii,this%x(1:3,ii) )
                    end if
                end do

                e2 = 0.0d0
                f0 = 0.0d0
                do kk = 1,nSphere
                    ii = id_sphere(kk)                      !   = 1:nAtoms , the number of the atom
                    call neighbourList( this%lc3d, this%x(1:3,ii) ,rc, nn,id,dx,dr2 )       !   returns index, vector separation and distance squared to all neighbours. Number of neighbours = nn.

                !   find self as a neighbour, and remove from neighbour list
                    jj = minloc(dr2(1:nn),dim=1)            !   closest neighbour to atom ii is itself.
                    id(jj) = id(nn)                         !   replace neighbour list entry jj with the last entry
                    dx(1:3,jj) = dx(1:3,nn)
                    dr2(jj) = dr2(nn)
                    nn = nn - 1
                    at(0) = this%at(ii)

                !   check which core atom is referenced, which neighbour
                    do jj = 1,nn
                        at(jj) = this%at( id(jj) )
                        id(jj) = indx( id(jj) )             !   either a core atom or to nCore+1

                    !   convert square distance and vector to distance and normal vector - this is what my force routine wants as input.
                        dd = sqrt(dr2(jj))
                        dr2(jj) = dd
                        dd = 1/dd
                        dx(1:3,jj) = dx(1:3,jj) * dd
                    end do

                    call addForce( this%eam, indx(ii) ,nn,at,dr2(1:nn),dx(1:3,1:nn),id(1:nn), f0 )
                    e2 = e2 + getEnergy(this%eam, nn,at,dr2(1:nn) )
                end do
                f2dotf2 = 0.0d0
                do ii = 1,nCore
                    f2dotf2 = f2dotf2 + dot_product( f0(1:3,ii),f0(1:3,ii) )
                end do
                if (SIMPLEMD_DBG) then
                    print *,"simpleMD::kickAtom info - rms force before, after, compensated  ",sqrt(f0dotf0)/(3*nCore),sqrt(f1dotf1)/(3*nCore),sqrt(f2dotf2)/(3*nCore)
                    print *,"simpleMD::kickAtom info - pot energy before, after, compensated ",e0,e1,e2
                    print *,"simpleMD::kickAtom info - pot energy delta after, compensated   ",(e1-e0),(e2-e0)
                end if


            !---    don't forget to restore the positions to their original state of the atoms
                do kk = 1,nSphere
                    ii = id_sphere(kk)
                    jj = indx(ii)
                    if (jj <= nCore) then
                    !   only cut and replace atoms in core region
                        call cut(this%lc3d,ii,this%x(1:3,ii) )
                        this%x(1:3,ii) = x0(1:3,jj)
                        call add( this%lc3d,ii,this%x(1:3,ii) )
                    end if
                end do

                 if (e2 > e1) then
                     if (SIMPLEMD_DBG) print *,"simpleMD::kickAtom info- test failed e2>e1"
!                     ok = .false.
!                     deallocate(id_sphere)
!                     deallocate(indx     )
!                     deallocate(x0       )
!                     deallocate(v0       )
!                     deallocate(f0       )
!                     deallocate(f1       )
!                     deallocate(uu       )
!                     deallocate(hess     )
!                     deallocate(dx       )
!                     deallocate(dr2      )
!                     deallocate(id       )
! !                     if (present(attempt)) then
! !                         if (attempt==MAX_ATTEMPTS) return
! !                         call kickAtom(this,atom,R0,Ekick,deltax , ok ,attempt+1)
! !                     else
! !                         call kickAtom(this,atom,R0,Ekick,deltax , ok ,attempt=1)
! !                     end if
!                     return
                  end if



               end if
!



        !---    test for unreasonable KE
            dd = 0.0d0
            do kk = 1,nSphere
                ii = id_sphere(kk)
                jj = indx(ii)
                if (jj <= nCore) then
                    alpha = this%mass(ii)*dot_product( uu( 3*jj-2:3*jj ) , uu( 3*jj-2:3*jj ) )
                    dd = max(dd,alpha)
                end if
            end do
            if (SIMPLEMD_DBG) print *,"simpleMD::kickAtom info - KE rescaling gives max energy ",dd/2

            if (dd>2*Ekick) then
                ok = .false.
                deallocate(id_sphere)
                deallocate(indx     )
                deallocate(x0       )
                deallocate(v0       )
                deallocate(f0       )
                deallocate(f1       )
                deallocate(uu       )
                deallocate(hess     )
                deallocate(dx       )
                deallocate(dr2      )
                deallocate(id       )

                if (present(attempt)) then
                    if (attempt==MAX_ATTEMPTS) then
                        if (SIMPLEMD_DBG) print *,"simpleMD::kickAtom info - too many attempts, return fail"
                        return
                    end if
                    if (SIMPLEMD_DBG) print *,"simpleMD::kickAtom info - try again"
                    call kickAtom_(this,atom,R0,Ekick,deltax , ok ,attempt+1)
                else
                    if (SIMPLEMD_DBG) print *,"simpleMD::kickAtom info - try again"
                    call kickAtom_(this,atom,R0,Ekick,deltax , ok ,attempt=1)
                end if
                return
            end if



            if (SIMPLEMD_DBG) then
                vpdotvp = dot_product( this%v(1:3,atom),this%v(1:3,atom) )
                print *,"simpleMD::kickAtom info - kinetic energy of kicked atom ",this%mass(atomindx)*vpdotvp/2
                if ( vpdotvp > 4*Ekick*this%imass(atomindx) ) print *, "simpleMD::kickAtom error - kinetic energy of kicked atom > 2 Ekick "
            end if



        !---    done, successfully.
            ok = .true.


        !---    add the velocities to the atoms

            do kk = 1,nSphere
                ii = id_sphere(kk)
                jj = indx(ii)
                if (jj <= nCore) then
                !   change velocity of atoms in core region
                    this%v(1:3,ii) =  uu( 3*jj-2:3*jj )
                end if
            end do

            if (SIMPLEMD_DBG) then
                print *,"simpleMD::kickAtom info - success"
                print *,""
            end if
            
            return
        end subroutine kickAtom_


        subroutine computeSaddlePointEnergy( this,x_old,x_new,atom,R0,e_old, nreplicas, dE , leaveInSaddle,ee_out,outfile,disp_scale )
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      attempt to compute the energy of a saddle point
    !*      if the atoms have moved from x_old to their current positions at x_new.
    !*      note that these are distinct from their positions within MD%x(:,:)
    !*      We are looking for a local saddle crossing in the region of atom
    !*      with a relaxation sphere radius R0. Return saddle barrier height dE.
    !*      WARNING - this leaves the link cell list in a disordered state...
    !*      it is expected that the caller will deal with the mess.
    !*      on input dE is the highest energy the saddle can reasonably be expected to reach.

            type(MD),intent(inout)          ::      this
            real(kind=real64),dimension(:,:),intent(in)     ::      x_old,x_new
            integer,intent(in)              ::      atom
            real(kind=real64),intent(in)    ::      e_old,R0
            real(kind=real64),intent(inout) ::      dE
            integer,intent(in)              ::      nReplicas
            logical,intent(in),optional     ::      leaveInSaddle
            real(kind=real64),dimension(:),intent(out),optional     ::      ee_out
            character(len=*),intent(in),optional                    ::      outfile
            real(kind=real64),intent(in),optional                   ::      disp_scale

            !integer                                             ::      NREPLICAS = 7   !   including "old" and "new"

            integer,parameter               ::      NSTEPS = 30


        !---    physically relevant properties
            integer                         ::      nSphere,nCore   !   number of atoms in sphere, number in core ( moving ) region
            real(kind=real64)               ::      rc              !   sphere cutoff radius, atom density
            real(kind=real64)               ::      max_disp2       !   largest displacement (squared) of any atom in sphere
            integer,dimension(:),allocatable,save               ::      indx            !   (1:nAtoms) = 1:nCore where in the core region is an atom
            integer,dimension(:),allocatable,save               ::      id_sphere       !   (1:nSphere) = 1:nAtoms = which atom is which in sphere region
            integer,dimension(:),allocatable,save               ::      at_sphere       !   (1:nSphere) = type of atom
            logical,save                                        ::      firstCall = .true.
            real(kind=real64),dimension(:),allocatable,save     ::      dr2_sphere      !   (1:nNeigh)          distance squared to a neighbour atom
            real(kind=real64),dimension(:,:,:),allocatable,save ::      x_replica       !   (3,nSphere,NREPLICAS)
            real(kind=real64),dimension(:),allocatable,save     ::      modx_replica    !   (NREPLICAS)
            real(kind=real64),dimension(:,:,:),allocatable,save ::      dx_replica       !   (3,nSphere,2:NREPLICAS-1)
            real(kind=real64),dimension(:,:),allocatable,save   ::      ff
            real(kind=real32),dimension(:,:),allocatable,save   ::      hess    !   (1:3*nCore+3,1:3*nCore+3) Hessian in core region
            
            real(kind=real32),dimension(:,:,:),allocatable,save     ::      hessEinstein
            
            real(kind=real64),dimension(:),allocatable          ::      ee,eold              !   energy of replicas

            real(kind=real64),dimension(7),parameter          ::      ALPHA = (/ -0.1d0,0.0d0,0.1d0,0.3d0,0.6d0,1.0d0,2.0d0 /)
            real(kind=real64),dimension(7)                    ::      eeee            !   energy of trial displacements

            real(kind=real64)               ::      dE_in       !   used as a sanity check 
!            real(kind=real64)               ::      neb_k = 5.0d0   !   nudged elastic band spring constant.
!            real(kind=real64)               ::      modx_bar,modx
            integer                         ::      saddle_replica
            real(kind=real64)               ::      displacement_scale
!            integer                         ::      best_displacement
            real(kind=real64),dimension(:,:),allocatable,save   ::      xhat_replica
            
        !---    computational workspace
            real(kind=real32),dimension(:),allocatable,save     ::      uu,work !   1:3*nCore+3 RHS of linear equation , LAPACK workspace
            integer,dimension(:),allocatable,save               ::      ipiv
            real(kind=real64),dimension(:),allocatable          ::      dr2     !   (1:nNeigh)          distance squared to a neighbour atom
            real(kind=real64),dimension(:,:),allocatable        ::      dx      !   (1:3,1:nNeigh)      vector separation to neighbour atom
            integer,dimension(:),allocatable                    ::      id      !   (1:nNeigh) =1:nAtoms index of atoms in neighbour list
            integer,dimension(:),allocatable                    ::      at      !   (1:nNeigh) = type of atom
            logical,dimension(:),allocatable                    ::      sane
            integer                                             ::      step

            real(kind=real64),dimension(3)              ::      e_trial,alpha_trial
            

        !---    general variables
            real(kind=real64)               ::      dd,x2
            real(kind=real64),dimension(3)  ::      xx,  offset 
            integer                         ::      ii,jj,kk,mm,nn
!            real(kind=real32),dimension(3)      ::      h1,h2,h3,gg

            real(kind=real64),dimension(3,3)    ::      aa,ia
            logical,dimension(:),allocatable    ::      hasx

            type(SimpleSupercell)   ::      super

            !  
!             integer                 ::      Nx,Ny,Nz !, ix,iy,iz
            character(len=256)      ::      dummy
            character(len=5)        ::      aaaa
            logical                 ::      ok
            integer,dimension(3)    ::      nnnn

            super = getSuper(this%lc3d)
            rc = getCutoff(this%eam)
            dE_in = dE
            dE = 0.0d0
            displacement_scale = rc/4
            if (present(disp_scale)) displacement_scale = disp_scale
            
            if (SIMPLEMD_DBG) print *,""
            if (SIMPLEMD_DBG) print *,"SimpleMD::computeSaddlePointEnergy() info - R0 = ",R0," dE_in = ",dE_in," atom = ",atom
            
            
            
        !---    find the number of atoms in the sphere
            call clear(this%lc3d)
!            do ii = 1,this%nAtoms
!                xx(1:3) = minimumImage( super,x_new(1:3,ii) - x_old(1:3,ii) )
!                this%x(1:3,ii) = x_old(1:3,ii) + xx(1:3)/2
!                call add(this%lc3d,ii,this%x(1:3,ii))
!            end do
!            call neighbourList_long( this%lc3d,this%x(1:3,atom),R0+2*rc, nSphere )                !   just returns number of atoms in sphere R0+2*rc
            call add(this%lc3d,x_old)
            call neighbourList_long( this%lc3d,x_old(1:3,atom),R0+2*rc, nSphere )                !   just returns number of atoms in sphere R0+2*rc

        !---    step 0: allocate memory to work on the atoms in the local sphere
            if (firstCall) then
                allocate(indx(nSphere))
                allocate(dr2_sphere(nSphere))
                allocate(id_sphere(nSphere))
                allocate(x_replica(3,nSphere,NREPLICAS))
                allocate(modx_replica(NREPLICAS))
                allocate(dx_replica(3,nSphere,NREPLICAS))
                allocate(xhat_replica(3,nSphere))
                allocate(at_sphere(nSphere))
            else
                if (nSphere > size(indx)) then
                    deallocate(indx)
                    deallocate(dr2_sphere)
                    deallocate(id_sphere)
                    deallocate(x_replica)
                    deallocate(modx_replica)
                    deallocate(dx_replica)
                    deallocate(xhat_replica)
                    deallocate(at_sphere)

                    allocate(indx(nSphere))
                    allocate(dr2_sphere(nSphere))
                    allocate(id_sphere(nSphere))
                    allocate(x_replica(3,nSphere,NREPLICAS))
                    allocate(modx_replica(NREPLICAS))
                    allocate(dx_replica(3,nSphere,NREPLICAS))
                    allocate(xhat_replica(3,nSphere))
                    allocate(at_sphere(nSphere))
                end if
            end if



            !call neighbourList_long( this%lc3d,this%x(1:3,atom),R0+2*rc, nSphere,id_sphere,dr2_sphere )
            call neighbourList_long( this%lc3d,x_old(1:3,atom),R0+2*rc, nSphere,id_sphere,dr2_sphere )
            nCore = count( dr2_sphere(1:nSphere) <= R0*R0 )        !   number which will become core
            if (firstCall) then
                allocate(ff(3,(nCore+1)))
                allocate(hess(3*nCore+3,3*nCore+3))
                allocate(work(66*(3*nCore+3)))
                allocate(ipiv(3*nCore+3))
                allocate(uu(3*nCore+3))
                allocate(hessEinstein(3,3,(nCore+1)))
            else
                if (3*nCore+3 > size(ipiv)) then
                    deallocate(ff)
                    deallocate(hess)
                    deallocate(work)
                    deallocate(ipiv)
                    deallocate(uu)
                    deallocate(hessEinstein)

                    allocate(ff(3,(nCore+1)))
                    allocate(hess(3*nCore+3,3*nCore+3))
                    allocate(work(66*(3*nCore+3)))
                    allocate(ipiv(3*nCore+3))
                    allocate(uu(3*nCore+3))
                    allocate(hessEinstein(3,3,(nCore+1)))
                end if
            end if
            firstCall = .false.

            
            

        !---    Check what the maximum displacement is of an atom in the core region.
        !       and escape if not large enough
        !       find which atoms are in the core region.
        !       place atoms in sphere region into their "old" positions.
        !       store the types of the atoms in the sphere region
             
            nn = 0
            max_disp2 = 0 ; jj = 0
            do kk = 1,nSphere
                ii = id_sphere(kk)
                at_sphere(kk) = this%at(ii)
                if ( dr2_sphere(kk) <= R0*R0 ) then
                    nn = nn + 1
                    indx(kk) = nn
                    xx(1:3) = minimumImage( super,x_new(1:3,ii) - x_old(1:3,ii) )
                    dd = xx(1)*xx(1) + xx(2)*xx(2) + xx(3)*xx(3)
                    if (dd>max_disp2) then
                        max_disp2 = dd
                        jj = ii
                    end if
                    !max_disp2 = max( max_disp2 , xx(1)*xx(1) + xx(2)*xx(2) + xx(3)*xx(3) )
                else if ( dr2_sphere(kk) <= (R0+rc)*(R0+rc) ) then
                    indx(kk) = nCore+1                        !   label index = nCore+1 for atoms which are not in the core region, but which directly contribute to force/hessian
                else
                    indx(kk) = -(nCore+1)                     !     label index = -(nCore+1) for atoms which do not directly contribute to force/hessian
                end if
            end do
            max_disp2 = sqrt(max_disp2)
            
            if ( (.not. SIMPLEMD_DBG).and.(max_disp2 < displacement_scale/2) ) return           !   saddle point barrier energy is returned as zero.
                
            !print *,"central atom ",atom
            if (SIMPLEMD_DBG) then
                print *,"max disp ",max_disp2," for atom ",jj
                print *,"old ",x_old(:,jj)," new ",x_new(:,jj)," A"
                print *,"rCore ",R0," rInner ",R0+rc," rOuter ",R0+2*rc
                print *,"nCore ",nn,nCore," nInner ",count(indx(1:nSphere)==(nCore+1))," nOuter ",count(indx(1:nSphere)<0)," total ",nSphere
            end if
            
            


        !---    centre the replicas
            
            allocate(ee(NREPLICAS))
            allocate(eold(NREPLICAS)) 
            allocate(sane(NREPLICAS))
             
            aa = getA(this%lc3d)
            call inverse3Mat(aa,ia)
            nnnn(1) = getNx(this%lc3d)
            nnnn(2) = getNy(this%lc3d)
            nnnn(3) = getNz(this%lc3d)
            allocate(hasx(0:maxval(nnnn)-1))  
            offset = 0.0d0 ; ok = .true.
            do jj = 1,3
                hasx = .false.
                do kk = 1,nSphere
                    ii = id_sphere(kk)
                    !ii = floor( ia(jj,1)*this%x(1,ii) + ia(jj,2)*this%x(2,ii) + ia(jj,3)*this%x(3,ii) )
                    ii = floor( ia(jj,1)*x_old(1,ii) + ia(jj,2)*x_old(2,ii) + ia(jj,3)*x_old(3,ii) )
                    ii = mod( ii + nnnn(jj)*8 , nnnn(jj) )
                    hasx(ii) = .true.
                end do
               ! write (*,fmt='(100l2)') hasx(0:nnnn(jj)-1)
                if (hasx(0) .and. hasx(nnnn(jj)-1)) then        !   spans x pbc
                    do ii = nnnn(jj)-2,1,-1
                        if (.not. hasx(ii)) then                !   cell ii is empty. Shift back to here
                            offset(1:3) = offset(1:3) + ii*aa(1:3,jj)
                            ok = .false.
                            exit
                        end if
                    end do
                end if
            end do
            deallocate(hasx)
           ! print *,"offset ",offset
            
           ! do kk = 1,nSphere
           !     ii = id_sphere(kk)  
           !     this%x(1:3,ii) = this%x(1:3,ii) - offset(1:3)
           ! end do
            


        !---    construct a set of replicas by making a evenly spaced string through the saddle
        !       note that after this step, all the x_replica(:) are within the bounds of the first periodic cell, and no longer need minimum image.
            x2 = 0.0d0
            do kk = 1,nSphere
                ii = id_sphere(kk)                      !   = 1:nAtoms
                
                do mm = 1,NREPLICAS
                    dd = real(mm-1,kind=real64)/(NREPLICAS-1)                  !    scales evenly from 0:1
                    xx(1:3) = minimumImage( super,x_new(1:3,ii) - x_old(1:3,ii) )
                    xhat_replica(1:3,kk) = xx(1:3)
                    xx(1:3) = x_old(1:3,ii) - offset(1:3) + xx(1:3)*dd 
                    x_replica(1:3,kk,mm) = wrapPBC( super,xx(1:3) )
                end do
                x2 = x2 + xhat_replica(1,kk)*xhat_replica(1,kk) + xhat_replica(2,kk)*xhat_replica(2,kk) + xhat_replica(3,kk)*xhat_replica(3,kk)
            end do
            x2 = 1/sqrt(x2)
            xhat_replica = xhat_replica*x2

            
!            if (SIMPLEMD_DBG) then
!                do mm = 1,NREPLICAS
!                    do jj = 1,3
!                        print *,"minmaxavg ",jj,mm,minval( x_replica( jj,1:nSphere,mm ) ),minval( x_replica( jj,1:nSphere,mm ) ),sum( x_replica( jj,1:nSphere,mm ) )/nSphere
!                    end do
!                end do
!            end if









        !---    compute the "initial" energy in the core region of the replicas.

            nn = getnNeighMax(this%lc3d)                !   maximum number of neighbours in neighbour list
            allocate(dx(3,nn))
            allocate(id(nn))
            allocate(at(0:nn))
            allocate(dr2(nn))

            ee = 0.0d0
            do mm = 1,NREPLICAS

                call clear( this%lc3d )
                call add( this%lc3d,x_replica(1:3,1:nSphere,mm) )


                do kk = 1,nSphere
                    if (indx(kk)<0) cycle                   !   only interested in energy in core + inner

                    call neighbourList( this%lc3d, x_replica(1:3,kk,mm) ,rc, nn,id,dr2 )

                !   find self as a neighbour, and remove from neighbour list
                    jj = minloc(dr2(1:nn),dim=1)            !   closest neighbour to atom ii is itself.
                    id(jj) = id(nn)                         !   replace neighbour list entry jj with the last entry
                    dr2(jj) = dr2(nn)
                    nn = nn - 1

                !   check which core atom type referenced
                    at(0) = at_sphere(kk)
                    do jj = 1,nn
                        at(jj) = at_sphere( id(jj) )
                        dr2(jj) = sqrt(dr2(jj))
                    end do

                    ee(mm) = ee(mm) + getEnergy(this%eam, nn,at,dr2 )
                end do

                if (SIMPLEMD_DBG) print *,"step ",0," replica ",mm," energy ",ee(mm)

            end do




        !---    quick escape for no saddle
            saddle_replica = maxloc( ee,dim=1 )
            if ( saddle_replica == 1 ) then
                !   no saddle can be found if the first (relaxed) end has a larger energy than the (unrelaxed) replicas
                dE = 0.0d0
                return            
            end if

            

            eold = ee(saddle_replica)
            do step = 1,NSTEPS
            
            !---    find distance between replicas
                modx_replica = 0
                do mm = 2,NREPLICAS 
                    do kk = 1,nSphere       !   note compute contribution from all atoms in sphere, even though some are outside core, they can still contribute inside core
                        !if (indx(kk)<0) cycle                  !   only interested in energy in core.
                        xx(1:3) = x_replica(1:3,kk,mm) - x_replica(1:3,kk,1)                           
                        modx_replica(mm) = modx_replica(mm) + dot_product( xx,xhat_replica(1:3,kk) )
                    end do
                   ! modx_replica(mm) = sqrt(modx_replica(mm))
                                        
                end do                 
               if (SIMPLEMD_DBG) then
                    do mm = 1,NREPLICAS 
                        print *,"step ",step," replica ",mm," |x| ",modx_replica(mm) 
                    end do
               end if 
                
            
            !---    relax the interior replicas
                dx_replica = 0.0d0
                do mm = 2,NREPLICAS-1

                    if ( (mm==1).or.(mm==NREPLICAS) ) then
                        !   test if ends are converged. They should be easier than the saddle
                        if ( eold(mm)-ee(mm) < 1.0d-4 ) cycle
                    end if
                
                
                    call clear( this%lc3d )
                    call add( this%lc3d,x_replica(1:3,1:nSphere,mm) )

 
                    
                !---    compute force and hessian in replica
                    ff = 0.0d0
                    hess = 0.0d0
                    uu = 0.0d0
                    
                    
                !---    add neb spring
                !    if ( (mm-1)*(mm-NREPLICAS)*(mm-saddle_replica) /= 0 ) then
                !        alpha_trial(1) = ( 1 - modx_bar/modx_replica(mm) )
                !        alpha_trial(2) = ( 1 - modx_bar/modx_replica(mm-1) )
                !        
                !        print *,"neb spring ",modx_replica(mm-1),modx_bar,modx_replica(mm)," alpha ",alpha_trial(1:2)
                !        
                !        do kk = 1,nSphere       !   note compute contribution from all atoms in sphere, even though some are outside core, they can still contribute inside core
                !            if (indx(kk)<0) cycle                  !   only interested in energy in core.
                !            
                !            xx(1:3) = x_replica(1:3,kk,mm+1) - x_replica(1:3,kk,mm)         !   vector linking replicas
                !            xx = minimumImage( super,xx ) * alpha_trial(1)          
                !                                     
                !            x2(1:3) = x_replica(1:3,kk,mm) - x_replica(1:3,kk,mm-1)         !   vector linking replicas
                !            x2 = minimumImage( super,x2 ) * alpha_trial(2)
                !             
                !            
                !            ff(1:3,indx(kk)) = neb_k*( xx(1:3) + x2(1:3) )
                !        end do    
                !    end if                                      


                    do kk = 1,nSphere       !   note compute contribution from all atoms in sphere, even though some are outside core, they can still contribute inside core
                        if (indx(kk)<0) cycle                  !   only interested in energy in core.

                        call neighbourList( this%lc3d, x_replica(1:3,kk,mm) ,rc, nn,id,dx,dr2 )

                    !   find self as a neighbour, and remove from neighbour list
                        jj = minloc(dr2(1:nn),dim=1)            !   closest neighbour to atom ii is itself.
                        id(jj) = id(nn)                         !   replace neighbour list entry jj with the last entry
                        dr2(jj) = dr2(nn)
                        dx(1:3,jj) = dx(1:3,nn)
                        nn = nn - 1

                    !   check which core atom type referenced
                        at(0) = at_sphere(kk)
                        do jj = 1,nn
                            at(jj) = at_sphere( id(jj) )
                            dd = sqrt(dr2(jj))
                            dr2(jj) = dd
                            dx(1:3,jj) = dx(1:3,jj)/dd
                            id(jj) = abs( indx( id(jj) ) )      !   abs because atoms outside R0+rc have indx = -nCore-1
                        end do

                        call addForce( this%eam, indx(kk) ,nn,at,dr2,dx,id, ff )
                        call addHessian( this%eam, indx(kk), nn,at,dr2,dx,id, hess )

                    end do
                     
                    
                    
                     
                    

                    !print *,"minimise"
                    if ( (mm-1)*(mm-NREPLICAS) == 0 ) then
                    !   straight relaxation of endpoints
                    
                        do jj = 1,nCore                                                        
                            uu( jj*3-2:jj*3 ) = real( ff(1:3,jj),kind=real32 )
                        end do
    
                    !       solve with LAPACK
                        call SSYSV( "U",3*nCore, 1, hess(1:3*nCore,1:3*nCore), 3*nCore,ipiv, uu(1:3*nCore), 3*nCore, work, size(work),ii )      !   note: LDA = 3*nCore
                        if (ii/=0) then
                            !   error in the linear algebra. Not seen one in the wild...
                            print *,"simpleMD::computeSaddlePointEnergy error - SSYSV returns ii = ",ii
                            dE = 0.0d0
                            return
                        end if
                    
                    else                    
                        
                    !---    now we want to minimise
                    !           E = E0 + E'.dx + 1/2 dx. E" dx
                    !       subject to the constraint
                    !           dx.v = 0
                    !       with v = direction along replicas
                    !           solution
                    !           H' dx = f'
                    !       with top left H' = E"           /   top f' = -E'
                    !       and bottom row/right column = v /   bottom f' = 0
    
                    !       will construct upper triangle H' and a single row vector u to capture the force
                    
                        
                        
                        dd = 0.0d0 
                        hess( :, 3*nCore+1 ) = 0.0d0
                        do kk = 1,nSphere
                            jj = abs(indx(kk))              !   where in the core is this sphere atom?
                            if (jj>nCore) cycle             !   only interested in core atoms
                            hess( jj*3-2:jj*3 , 3*nCore+1 ) = real( xhat_replica(1:3,kk),kind=real32 )    
                            dd = dd + dot_product( ff(1:3,jj),xhat_replica(1:3,kk) )
                        end do
                        if (SIMPLEMD_DBG) print *,"f.xhat ",dd
                        uu( 3*nCore+1: ) = 0.0d0
                        
                        do kk = 1,nSphere
                            jj = abs(indx(kk))              !   where in the core is this sphere atom?
                            if (jj>nCore) cycle  
                            uu( jj*3-2:jj*3 ) = real( ff(1:3,jj) - dd*xhat_replica(1:3,kk) ,kind=real32 )                             
                        end do
                         
                        if (SIMPLEMD_DBG) then
                            dd = 0.0d0 
                            do kk = 1,nSphere
                                jj = abs(indx(kk))              !   where in the core is this sphere atom?
                                if (jj>nCore) cycle             !   only interested in core atoms
                                dd = dd + dot_product( uu( jj*3-2:jj*3 ),xhat_replica(1:3,kk) )
                            end do
                            print *,"u.xhat ",dd
                        end if
                        
                        
                    !       solve with LAPACK
                         call SSYSV( "U",3*nCore+1, 1, hess(1:3*nCore+1,1:3*nCore+1), 3*nCore+1,ipiv(1:3*nCore+1), uu(1:3*nCore+1), 3*nCore+1, work, size(work),ii )      !   note: LDA = 3*nCore+1
                       ! call SSYSV( "U",3*nCore+3, 1, hess(1:3*nCore+3,1:3*nCore+3), 3*nCore+3,ipiv, uu, 3*nCore+3, work, size(work),ii )      !   note: LDA = 3*nCore+1
                        !call SSYSV( "U",3*nCore+1, 1, hess, 3*nCore+1,ipiv, uu, 3*nCore+1, work, size(work),ii )      !   note: LDA = 3*nCore
                        if (ii/=0) then
                            !   error in the linear algebra. Not seen one in the wild...
                            print *,"simpleMD::computeSaddlePointEnergy error - SSYSV returns ii = ",ii
                            dE = 0.0d0
                            return
                        end if
                        
                        if (SIMPLEMD_DBG) then
                            dd = 0.0d0  
                            do kk = 1,nSphere
                                jj = abs(indx(kk))              !   where in the core is this sphere atom?
                                if (jj>nCore) cycle             !   only interested in core+inner atoms
                                dd = dd + dot_product( uu( jj*3-2:jj*3 ),xhat_replica(1:3,kk)  )
                            end do
                            print *,"u.xhat ",dd
                        end if
                        
                        
                        
                      !  do kk = 1,nSphere
                      !      jj = abs(indx(kk))              !   where in the core is this sphere atom?
                      !      if (jj>nCore) cycle  
                      !      uu( jj*3-2:jj*3 ) = uu( jj*3-2:jj*3 ) - dd*xhat_replica(1:3,kk)
                      !  end do
                      !    
                      !  
                      !
                      !  dd = 0.0d0
                      !  do kk = 1,nSphere
                      !      jj = abs(indx(kk))              !   where in the core is this sphere atom?
                      !      if (jj>nCore) cycle             !   only interested in core+inner atoms
                      !      dd = dd + dot_product( uu( jj*3-2:jj*3 ),xhat_replica(1:3,kk) )
                      !  end do
                      !  
                      !  print *," u.dx = ",dd
                        
                    end if                    
                    
                    !print *,"unpack"
                !       unpack solution vector into dx_replica

                    do kk = 1,nSphere
                        jj = abs(indx(kk))              !   where in the core is this sphere atom?
                        if (jj>nCore) cycle             !   only move core atoms
                        dx_replica(1:3,kk,mm) = uu(3*jj-2:3*jj)
                    end do
                     
                    
                !       carefully move in direction of solution vector                    
                    call computeEnergyFnVecLength( this,nSphere,x_replica(:,:,mm),dx_replica(:,:,mm),at_sphere,ALPHA,eeee,displacement_scale,indx>0,xhat_replica )


                !       where is the lowest energy?
                    ii = minloc( eeee,dim=1 )                                   
                    if ( (ii-1)*(size(eeee)-ii) == 0 ) then             ! at one or other end
                        dd = ALPHA(ii)
                    else
                        !   compute the best ( quadratic approx ) displacement
                        dd = ALPHA(ii-1)*(eeee(ii)-eeee(ii+1)) + ALPHA(ii)*(eeee(ii+1)-eeee(ii-1)) + ALPHA(ii+1)*(eeee(ii-1)-eeee(ii))
                        if (abs(dd)>1.0d-16) dd = 1/(2*dd)
                        dd = ( ALPHA(ii-1)*ALPHA(ii-1)*(eeee(ii)-eeee(ii+1)) + ALPHA(ii)*ALPHA(ii)*(eeee(ii+1)-eeee(ii-1)) + ALPHA(ii+1)*ALPHA(ii+1)*(eeee(ii-1)-eeee(ii)) )*dd
                    end if    
                    
                    
                !---    sanity check with my assumed best displacement dd near minloc ii.
                
                    if (SIMPLEMD_DBG) write (*,fmt='(a,3i4,100f16.6)') "replica",mm,minloc( eeee,dim=1 ) ,ii,eeee,dd
                    ok = .true.
                    call computeEnergyFnVecLength( this,nSphere,x_replica(:,:,mm),dx_replica(:,:,mm),at_sphere,(/dd/),e_trial,displacement_scale,indx>0,xhat_replica )   
                    !print *,"sanity check 1 - displacement ",dd," energy ",e_trial(1)
                    if (e_trial(1) > eold(mm)) then                                         
                        if (SIMPLEMD_DBG) print *,"sanity check 1 failed at ",dd," energy ",e_trial(1)," > ",eold(mm)
                        if ( (ii-1)*(NREPLICAS-ii) == 0 ) then
                            !   can't go any further here. set displacement to zero
                            ok = .false.                       
                        else
                            !   have more data to improve displacement guess. Try again.
                            if (dd < ALPHA(ii)) then
                                alpha_trial(1:3) = (/ ALPHA(ii-1),dd,ALPHA(ii) /)
                                e_trial(1:3) = (/ eeee(ii-1),e_trial(1),eeee(ii) /)
                            else
                                alpha_trial(1:3) = (/ ALPHA(ii),dd,ALPHA(ii+1) /)
                                e_trial(1:3) = (/ eeee(ii),e_trial(1),eeee(ii+1) /)
                            end if
                            if (SIMPLEMD_DBG) write (*,fmt='(a,3i4,100f16.6)') "replica",mm,minloc( e_trial,dim=1 ) ,2,e_trial,dd
                            dd = alpha_trial(1)*(e_trial(2)-e_trial(3)) + alpha_trial(2)*(e_trial(3)-e_trial(1)) + alpha_trial(3)*(e_trial(1)-e_trial(2))
                            if (abs(dd)>1.0d-16) dd = 1/(2*dd)
                            dd = ( alpha_trial(1)*alpha_trial(1)*(e_trial(2)-e_trial(3)) + alpha_trial(2)*alpha_trial(2)*(e_trial(3)-e_trial(1)) + alpha_trial(3)*alpha_trial(3)*(e_trial(1)-e_trial(2)) )*dd
                            call computeEnergyFnVecLength( this,nSphere,x_replica(:,:,mm),dx_replica(:,:,mm),at_sphere,(/dd/),e_trial,displacement_scale,indx>0,xhat_replica )   
                            !print *,"sanity check 2 - displacement ",dd," energy ",e_trial(1)
                            if (e_trial(1) > eold(mm)) then      
                                if (SIMPLEMD_DBG) print *,"sanity check 2 failed at ",dd," energy ",e_trial(1)," > ",eold(mm)
                                ok = .false.
                                sane(mm) = .false.
                                
                            else if (e_trial(1) > eeee(ii)) then
                                if (SIMPLEMD_DBG) print *,"sanity check 2 could be improved at ",dd," energy ",e_trial(1)," > ",eeee(ii)     
                                dd = ALPHA(ii)                       
                                e_trial(1) = eeee(ii)
                                sane(mm) = .true.
                            else  
                                if (SIMPLEMD_DBG) print *,"sanity check 2 passed at ",dd," energy ",e_trial(1)," < ",eold(mm)
                                sane(mm) = .true.
                            end if
                        end if
                    else if (e_trial(1) > eeee(ii)) then
                        if (SIMPLEMD_DBG) print *,"sanity check 1 could be improved at ",dd," energy ",e_trial(1)," > ",eeee(ii)     
                        dd = ALPHA(ii)                       
                        e_trial(1) = eeee(ii)
                        sane(mm) = .true.
                    else
                        if (SIMPLEMD_DBG) print *,"sanity check 1 passed at ",dd," energy ",e_trial(1)," < ",eold(mm)
                        sane(mm) = .true.
                    end if        
                    
                    if (ok) then
                        call scaleDisplacementVector( dx_replica(:,:,mm),dd,displacement_scale,xhat_replica )
                        !x_replica(1:3,1:nSphere,mm) = x_replica(1:3,1:nSphere,mm) + dx_replica(1:3,1:nSphere,mm)
                        ee(mm) = e_trial(1)
                    else
                        !call random_number( dx_replica(:,:,mm) )
                        !dx_replica(:,:,mm) = (dx_replica(:,:,mm) - 0.5d0)
                        !call scaleDisplacementVector( dx_replica(:,:,mm),0.001d0,rc/4,xhat_replica )
                        dx_replica(:,:,mm) = 0.0d0
                        ee(mm) = eold(mm)
                    end if
                            
                     dd = 0.0d0  
                        do kk = 1,nSphere
                        jj = abs(indx(kk))              !   where in the core is this sphere atom?
                        if (jj>nCore) cycle             !   only interested in core+inner atoms
                        dd = dd + dot_product( dx_replica(1:3,kk,mm),xhat_replica(1:3,kk)  )
                    end do
                   
                    if (SIMPLEMD_DBG) print *,"dx.xhat ",dd       
                            
                            
                                                 
                    !   save the displacement needed
                    
                     
                    max_disp2 = 0.0d0
                    do kk = 1,nSphere
                        max_disp2 = max( max_disp2,norm2( dx_replica(1:3,kk,mm) ) )
                    end do
                    if (SIMPLEMD_DBG) write (*,fmt='(a,3i4,100f16.6)') "replica",mm,minloc( eeee,dim=1 ) ,ii,eeee,dd   , max_disp2

                end do      !   mm = 2,REPLICAS-1

                
                do mm = 2,NREPLICAS-1 
                    if (.not. sane(mm)) then
                        if (sane(mm-1) .and. sane(mm+1)) then
                            if (SIMPLEMD_DBG) print *,"sanitising replica ",mm," by redrawing as midpoint"
                       
                            do kk = 1,nSphere
                                x_replica(1:3,kk,mm) = ( x_replica(1:3,kk,mm+1)+x_replica(1:3,kk,mm-1) )/2
                            end do
                            dx_replica(1:3,1:nSphere,mm) = 0.0d0
                        else
                            if (SIMPLEMD_DBG) print *,"attempting to sanitise replica ",mm," by jiggling"
                            call random_number( dx_replica(:,:,mm) )
                            dx_replica(:,:,mm) = (dx_replica(:,:,mm) - 0.5d0)
                            call scaleDisplacementVector( dx_replica(:,:,mm),0.001d0,displacement_scale,xhat_replica )
                        end if
                    end if
                end do
                
                
                

            !---    move the replicas using the stored displacements
                !modx = 0.0d0
                if (SIMPLEMD_DBG) print *,"step ",step," replica ",1," disp ",0.0d0," energy ",ee(1)
                do mm = 2,NREPLICAS-1 
!
                   x_replica(1:3,1:nSphere,mm) = x_replica(1:3,1:nSphere,mm) + dx_replica(1:3,1:nSphere,mm)
                   call clear( this%lc3d )
                   call add( this%lc3d,x_replica(1:3,1:nSphere,mm) )

                   ee(mm) = 0.0d0
                   do kk = 1,nSphere
                       if (indx(kk)<0) cycle                  !   only interested in energy in core.

                       call neighbourList( this%lc3d, x_replica(1:3,kk,mm) ,rc, nn,id,dr2 )

                   !   find self as a neighbour, and remove from neighbour list
                       jj = minloc(dr2(1:nn),dim=1)            !   closest neighbour to atom ii is itself.
                       id(jj) = id(nn)                         !   replace neighbour list entry jj with the last entry
                       dr2(jj) = dr2(nn)
                       nn = nn - 1

                   !   check which core atom type referenced
                       at(0) = at_sphere(kk)
                       do jj = 1,nn
                           at(jj) = at_sphere( id(jj) )
                           dr2(jj) = sqrt(dr2(jj))
                       end do

                       ee(mm) = ee(mm) + getEnergy(this%eam, nn,at,dr2 )
                   end do
                   
                   !if (ee(mm) > eold(mm)) then
                   !     if (SIMPLEMD_DBG) print *,"sanity check 3 failed energy ",ee(mm)," > ",eold(mm)
                   !     x_replica(1:3,1:nSphere,mm) = x_replica(1:3,1:nSphere,mm) - dx_replica(1:3,1:nSphere,mm)
                   !     ee(mm) = eold(mm)
                   !end if
                   ! 
                    !modx = modx + modx_replica(mm-1)/(modx_bar*NREPLICAS)
                    if (SIMPLEMD_DBG) print *,"step ",step," replica ",mm," disp ",modx_replica(mm)," energy ",ee(mm)," sane ",sane(mm)

                end do     
                !modx = modx + modx_replica(mm)/(modx_bar*NREPLICAS)
                if (SIMPLEMD_DBG) print *,"step ",step," replica ",NREPLICAS," disp ",modx_replica(NREPLICAS)," energy ",ee(NREPLICAS)

                mm = maxloc( eold-ee,dim=1)
                dE = eold(mm) - ee(mm)
                saddle_replica = maxloc(ee,dim=1) 
                if (SIMPLEMD_DBG) print *,"step ",step," de ",dE," at ",mm," saddle = ",maxval(ee)-ee(1)," at ",saddle_replica
                if (SIMPLEMD_DBG) print *,""
                if (dE < 1.0d-4) exit                       !   converged energies
                
                if ( max( eold(1)-ee(1),eold(NREPLICAS)-ee(NREPLICAS) )<1.0d-4 ) then
                    !   both ends are converged.                    
                     
                    if ( (saddle_replica-1)*(saddle_replica-NREPLICAS) == 0 ) exit      !   no barrier peak found
                end if
                    
                eold = ee
                
                
            end do      !   step

!        !---    unpack saddle position
!            x_replica(1:3,1:nSphere,mm) = x_replica(1:3,1:nSphere,mm) - dx_replica(1:3,1:nSphere,mm)
!            call clear(this%lc3d)
!            do kk = 1,nSphere
!                ii = id_sphere(kk)
!                xx(1:3) = x_replica(1:3,kk,saddle_replica) - x_replica(1:3,kk,1)
!                this%x(1:3,ii) = wrapPBC( super,x_old(1:3,ii) + xx(1:3) )
!            end do
!            call add(this%lc3d,this%x(1:3,1:this%nAtoms))
            
            if (present(outfile)) then
                do mm = 1,NREPLICAS
                    write(aaaa,fmt='(i5)') mm ; aaaa = adjustl(aaaa) ; aaaa = repeat("0",5-len_trim(aaaa)) // trim(aaaa)
                    open(unit=880,file=trim(outfile)//"."//trim(aaaa)//".xyz",action="write")
                         write(unit=880,fmt='(i8)') nSphere
                         write(dummy,fmt='(9f12.4)') getSuperA( super )
                         write(unit=880,fmt='(a)') "Lattice="""//trim(dummy)//""" Properties=Species:S:1:Pos:R:3:Energy:R:1:n:I:1"
                         call clear( this%lc3d )                          
                         call add( this%lc3d,x_replica(1:3,1:nSphere,mm) )
                         
                         do kk = 1,nSphere
                            call neighbourList( this%lc3d, x_replica(1:3,kk,mm) ,rc, nn,id,dr2 )
                
                        !   find self as a neighbour, and remove from neighbour list
                            jj = minloc(dr2(1:nn),dim=1)            !   closest neighbour to atom ii is itself.
                            id(jj) = id(nn)                         !   replace neighbour list entry jj with the last entry
                            dr2(jj) = dr2(nn)
                            nn = nn - 1
                
                        !   check which core atom type referenced
                            at(0) = at_sphere(kk)
                            do jj = 1,nn
                                at(jj) = at_sphere( id(jj) )
                                dr2(jj) = sqrt(dr2(jj))
                            end do
                            dd = getEnergy(this%eam, nn,at,dr2 )
                
                            write(unit=880,fmt='(a,4f16.6,i8)') getName(this%eam,at(0)),x_replica(1:3,kk,mm),dd,indx(kk)
                
                         end do
                    close(unit=880)
                end do
            end if



             
            if (saddle_replica<NREPLICAS) then
                dE = maxval(ee) - ee(1)
            else
                dE = 0
            end if
            write(*,fmt='(a,i8,a,3f12.4,a,f16.4)') "simpleMD::computeSaddlePointEnergy info - atom ",atom,"( "//trim(getName(this%eam,this%at(atom)))//" ) at ",x_old(1:3,atom)," crosses saddle barrier ",dE
            
            !if (dE>2*dE_in) stop

            if (present(leaveInSaddle)) then     
                if (leaveInSaddle) then           
                    
                    this%x(1:3,1:this%nAtoms) = x_old(1:3,1:this%nAtoms)
                    do kk = 1,nSphere
                        ii = id_sphere(kk)
                        xx(1:3) = x_replica(1:3,kk,saddle_replica) - x_replica(1:3,kk,1)
                        this%x(1:3,ii) = wrapPBC( super,x_old(1:3,ii) + xx(1:3) )
                    end do
                    call clear(this%lc3d)
                    call add(this%lc3d,this%x(1:3,1:this%nAtoms))
                end if
            end if    
            
            if (present(ee_out)) then
                ee_out(1:NREPLICAS) = ee(1:NREPLICAS)
            end if
            
            
            return





            return
        end subroutine computeSaddlePointEnergy



        
        
        
        
        
        
        
        
        
        

        subroutine computeSaddlePointEnergy_old( this,x_old,x_new,atom,R0,e_old, dE )
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      attempt to compute the energy of a saddle point
    !*      if the atoms have moved from x_old to their current positions
    !*      and we are looking for a local saddle crossing in the region of atom
    !*      with a relaxation sphere radius R0

            type(MD),intent(inout)          ::      this
            real(kind=real64),dimension(:,:),intent(in)     ::      x_old,x_new
            integer,intent(in)              ::      atom
            real(kind=real64),intent(in)    ::      e_old,R0
            real(kind=real64),intent(inout)   ::      dE



        !---    physical quantities
            integer                         ::      nCore                   !   number of atoms given velocity total ( = local region )
            integer                         ::      nSphere                 !   number of atoms that need consideration for force calculation ( = local region + buffer )
            real(kind=real64)               ::      rc                      !   cut-off radius for potential


            real(kind=real64),dimension(:,:),allocatable    ::      f0!,f1   !   (1:3,1:nCore+1) force on atoms in core region
            real(kind=real32),dimension(:,:),allocatable    ::      hess    !   (1:3*nCore+3,1:3*nCore+3) Hessian in core region
            !real(kind=real64),dimension(:,:),allocatable    ::      x_stored     !   (1:3,1:nAtoms) position of atoms

            integer                                         ::      atomindx    !   where is the kicked atom in the local region? indx(atom) = atomindx

            !real(kind=real64)               ::      e0,e1,e2,e3,e4,eHarm                  !   potential energy of atoms before and after


            real(kind=real64)               ::      e_before,e_saddle_unrelaxed,e_saddle_anharmonic,e_saddle_harmonic         !,e_after

        !---    computational workspace
            real(kind=real64),dimension(:),allocatable      ::      dr2     !   (1:nNeigh)          distance squared to a neighbour atom
            real(kind=real64),dimension(:,:),allocatable    ::      dx      !   (1:3,1:nNeigh)      vector separation to neighbour atom
            integer,dimension(:),allocatable                ::      id              !   =1:nAtoms index of atoms in neighbour list
            integer,dimension(:),allocatable                ::      id_sphere       !   =1:nAtomsindex of atoms in local+buffer region
            integer,dimension(:),allocatable                ::      indx    !   (1:nAtoms) = 1:nCore where in the core region is an atom
            integer,dimension(:),allocatable                ::      at
!           real(kind=real64)                               ::      udotu,v0dotv0,udotv0,vpdotu,vpdotv0,udotf0,f0dotf0,f1dotf1,f2dotf2,vpdotvp

            real(kind=real32),dimension(:),allocatable      ::      uu,work !   1:3*nCore+3 RHS of linear equation , LAPACK workspace
            integer,dimension(:),allocatable                ::      ipiv

        !---    dummy
            real(kind=real64),dimension(3)  ::      deltax
            real(kind=real64)               ::      dd , max_disp2
            integer                         ::      ii,jj,kk,nn



        !---    find the properties of the core region
            rc = getCutoff(this%eam)
            deltax = (x_old(1:3,atom) + x_new(1:3,atom))/2
            call neighbourList_long( this%lc3d,deltax,R0+rc, nSphere )        !   just returns number of atoms in sphere R0+rc
            allocate(dr2(nSphere))
            allocate(id_sphere(nSphere))
            call neighbourList_long( this%lc3d,deltax,R0+rc, nSphere,id_sphere,dr2 )


        !---    count and identify the atoms in the core region. compute "before" positions, forces, velocities.
            nCore = count( dr2(1:nSphere) <= R0*R0 )
            if (SIMPLEMD_DBG) write(*,fmt='(a,f12.4,4(a,i8))') "simpleMD::computeSaddlePointEnergy info - sphere range ",R0," contains ",nCore,"/",nSphere," atoms centred on ",atom







        !---    store the position of the atoms currently in the sphere. Index the local atoms
            allocate(indx(this%nAtoms))
            indx(:) = nCore+1           !   denote "outside the core region" as nCore + 1.
            nCore = 0
            do kk = 1,nSphere
                if ( dr2(kk) < R0*R0 ) then
                    ii = id_sphere(kk)
                    nCore = nCore + 1
                    if (ii==atom) atomindx = nCore
                    indx(ii) = nCore
                end if
            end do
            deallocate(dr2)


        !---    did the atoms actually go anywhere in the core region?

            max_disp2 = 0.0d0
            do kk = 1,nSphere
                ii = id_sphere(kk)
                jj = indx(ii)
                if (jj <= nCore) then
                    deltax(1:3) = x_new(1:3,ii) - x_old(1:3,ii)
                    dd = deltax(1)*deltax(1) + deltax(2)*deltax(2) + deltax(3)*deltax(3)
                    max_disp2 = max(dd,max_disp2)
                end if
            end do
            if (SIMPLEMD_DBG) print *,"simpleMD::computeSaddlePointEnergy info - max atom displacement ",sqrt(max_disp2)
            if (max_disp2 < 0.25d0) then
                dE = 0.0d0      !   saddle point barrier energy is returned
                return
            end if


            nn = getnNeighMax(this%lc3d)                !   maximum number of neighbours in neighbour list
            allocate(dx(3,nn))
            allocate(id(nn))
            allocate(at(0:nn))
            allocate(dr2(nn))

            allocate(f0(3,nCore+1))
            !allocate(f1(3,nCore+1))
            allocate(hess(3*(nCore+1),3*(nCore+1)))
            allocate(work(66*(3*nCore)))
            allocate(ipiv(3*nCore))
            allocate(uu(3*nCore))

        !---    compute energy within the sphere ( = core + rc ) before the transition

            do kk = 1,nSphere
                ii = id_sphere(kk)
                call cut(this%lc3d,ii,this%x(1:3,ii) )
                this%x(1:3,ii) = x_old(1:3,ii)
                call add( this%lc3d,ii,this%x(1:3,ii) )
            end do

            e_before = 0.0d0
            do kk = 1,nSphere
                ii = id_sphere(kk)                      !   = 1:nAtoms , the true number of the atom
                if (indx(ii)>nCore) cycle               !   only interested in energy in core.
                call neighbourList( this%lc3d, this%x(1:3,ii) ,rc, nn,id,dr2 )       !   returns index, vector separation and distance squared to all neighbours. Number of neighbours = nn.

            !   find self as a neighbour, and remove from neighbour list
                jj = minloc(dr2(1:nn),dim=1)            !   closest neighbour to atom ii is itself.
                id(jj) = id(nn)                         !   replace neighbour list entry jj with the last entry
                dx(1:3,jj) = dx(1:3,nn)
                dr2(jj) = dr2(nn)
                nn = nn - 1

                at(0) = this%at(ii)
            !   check which core atom is referenced, which neighbour
                do jj = 1,nn
                    at(jj) = this%at( id(jj) )
                    dr2(jj) = sqrt(dr2(jj))
                end do

                !print *,"before atom in sphere ",kk," in core ",indx(ii)," in .xyz ",ii," neighbours ",nn," energy ", getEnergy(this%eam, nn,at,dr2 )
                e_before = e_before + getEnergy(this%eam, nn,at,dr2 )
            end do
            !if (SIMPLEMD_DBG) &
            print *,"energy before in core ",e_before







            if (SIMPLEMD_DBG) &
            print *,"simpleMD::computeSaddlePointEnergy info - atom at centre is ",atom," at ",this%x(1:3,atom)," = ",atomindx,"/",nCore





        !---    now place atoms into positions close to the saddle point
            do kk = 1,nSphere
                ii = id_sphere(kk)
            !   only cut and replace atoms in core region
                call cut(this%lc3d,ii,this%x(1:3,ii) )
                this%x(1:3,ii) = ( x_old(1:3,ii) + x_new(1:3,ii) )/2
                call add( this%lc3d,ii,this%x(1:3,ii) )
            end do





        !---    compute energy, force, hessian within the sphere ( = core + rc ) at the saddle point
            e_saddle_unrelaxed = 0.0d0
            f0 = 0.0d0
            hess = 0.0d0
            do kk = 1,nSphere
                ii = id_sphere(kk)                      !   = 1:nAtoms , the true number of the atom
                call neighbourList( this%lc3d, this%x(1:3,ii) ,rc, nn,id,dx,dr2 )       !   returns index, vector separation and distance squared to all neighbours. Number of neighbours = nn.

            !   find self as a neighbour, and remove from neighbour list
                jj = minloc(dr2(1:nn),dim=1)            !   closest neighbour to atom ii is itself.
                id(jj) = id(nn)                         !   replace neighbour list entry jj with the last entry
                dx(1:3,jj) = dx(1:3,nn)
                dr2(jj) = dr2(nn)
                nn = nn - 1

                at(0) = this%at(ii)
            !   check which core atom is referenced, which neighbour
                do jj = 1,nn
                    at(jj) = this%at( id(jj) )
                    id(jj) = indx( id(jj) )             !   either a core atom or to nCore+1

                !   convert square distance and vector to distance and normal vector - this is what my force routine wants as input.
                    dd = sqrt(dr2(jj))
                    dr2(jj) = dd
                    dd = 1/dd
                    dx(1:3,jj) = dx(1:3,jj) * dd

                end do

                call addForce( this%eam, indx(ii) ,nn,at,dr2,dx,id, f0 )
                call addHessian( this%eam, indx(ii), nn,at,dr2,dx,id, hess )

                if (indx(ii) <= nCore) then
                    !print *,"unrelx atom in sphere ",kk," in core ",indx(ii)," in .xyz ",ii," neighbours ",nn," energy ", getEnergy(this%eam, nn,at,dr2 )
                    e_saddle_unrelaxed = e_saddle_unrelaxed + getEnergy(this%eam, nn,at,dr2 )
                end if

            end do
            !if (SIMPLEMD_DBG) &
            print *,"unrelaxed saddle ",e_saddle_unrelaxed







        !---    now compute the _relaxed_ atom positions near to the old position
        !---    solve harmonic minimisation with LAPACK

            do jj = 1,nCore         !   make a 32 bit copy of the forces
                uu(3*jj-2:3*jj) = real(f0(1:3,jj),kind=real32)
            end do
            call SSYSV( "U",3*nCore,1,hess(1:3*nCore,1:3*nCore),3*nCore,ipiv,uu,3*nCore,work,size(work),ii )      !   note: LDA = 3*nCore
            if (ii/=0) then
                !   error in the linear algebra. Not seen one in the wild...
                if (SIMPLEMD_DBG) &
                print *,"simpleMD::computeSaddlePointEnergy error - SSYSV returns ii = ",ii
                dE = 0.0d0
                return
            end if
            e_saddle_harmonic = 0.0d0
            do jj = 1,nCore
                e_saddle_harmonic = e_saddle_harmonic + dot_product( f0(1:3,jj) , uu(3*jj-2:3*jj) )
                f0(1:3,jj) = uu(3*jj-2:3*jj)
            end do
            e_saddle_harmonic = -e_saddle_harmonic/2
            !if (SIMPLEMD_DBG) &
            print *,"harmonic saddle ",e_saddle_unrelaxed,"+",e_saddle_harmonic," = ",e_saddle_unrelaxed+e_saddle_harmonic

        !---    I want a good calculation of the saddle energy too...
        !   place atoms at saddle point
            do kk = 1,nSphere
                ii = id_sphere(kk)
                jj = indx(ii)
                if (jj <= nCore) then           !   only cut and replace atoms in core region. Note that in the sphere the atoms are midway between old & new.
                    call cut(this%lc3d,ii,this%x(1:3,ii) )
                    this%x(1:3,ii) = this%x(1:3,ii) + f0(1:3,jj)
                    call add( this%lc3d,ii,this%x(1:3,ii) )
                end if
            end do


       !---     compute anharmonic energy at saddle position in core region
            e_saddle_anharmonic = 0.0d0
            do kk = 1,nSphere
                ii = id_sphere(kk)                      !   = 1:nAtoms , the true number of the atom
                if (indx(ii)>nCore) cycle               !   only interested in energy in core.
                call neighbourList( this%lc3d, this%x(1:3,ii) ,rc, nn,id,dr2 )       !   returns index, vector separation and distance squared to all neighbours. Number of neighbours = nn.

            !   find self as a neighbour, and remove from neighbour list
                jj = minloc(dr2(1:nn),dim=1)            !   closest neighbour to atom ii is itself.
                id(jj) = id(nn)                         !   replace neighbour list entry jj with the last entry
                dr2(jj) = dr2(nn)
                nn = nn - 1

                at(0) = this%at(ii)
            !   check which core atom is referenced, which neighbour
                do jj = 1,nn
                    at(jj) = this%at( id(jj) )
                    dr2(jj) = sqrt( dr2(jj) )
                end do
                !print *,"anharm atom in sphere ",kk," in core ",indx(ii)," in .xyz ",ii," neighbours ",nn," energy ", getEnergy(this%eam, nn,at,dr2 )
                e_saddle_anharmonic = e_saddle_anharmonic + getEnergy(this%eam, nn,at,dr2 )
            end do
            !if (SIMPLEMD_DBG) &
            print *,"anharm saddle ",e_saddle_anharmonic




            !if (SIMPLEMD_DBG) &
            print *,"computeSaddlePointEnergy info - relaxation of saddle state ",e_before," -> ",e_saddle_unrelaxed,e_saddle_anharmonic," de = ",  &
                            e_saddle_unrelaxed-e_before,e_saddle_anharmonic-e_before,e_saddle_unrelaxed-e_before+e_saddle_harmonic

            dE = e_saddle_unrelaxed+e_saddle_harmonic-e_before
            !if (SIMPLEMD_DBG) &
            print *,"computeSaddlePointEnergy info - saddle point barrier energy ",dE




        !---    return atoms to original positions
 
             do kk = 1,nSphere
                ii = id_sphere(kk)
                call cut(this%lc3d,ii,this%x(1:3,ii) )
                this%x(1:3,ii) = x_new(1:3,ii)
                call add( this%lc3d,ii,this%x(1:3,ii) )
            end do


            return
        end subroutine computeSaddlePointEnergy_old







!         subroutine kickAtom_Old(this,atom,R0,Ekick,deltax , ok)
!     !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
!             type(MD),intent(inout)          ::      this
!             integer,intent(in)              ::      atom
!             real(kind=real64),intent(in)    ::      R0
!             real(kind=real64),intent(in)    ::      Ekick
!             real(kind=real64),intent(in)    ::      deltax
!             logical,intent(out)             ::      ok
!
!
!             integer                                         ::      nSphere,nCore
!             real(kind=real64),dimension(:),allocatable      ::      dr,dr2
!             integer,dimension(:),allocatable                ::      id,id2,indx
!             real(kind=real64),dimension(:,:),allocatable    ::      force,dx,f1
!             integer                                         ::      coreatom
!
!             integer             ::      k2,ii,jj,kk,nn
!
!             real(kind=real64)   ::      e0,e1,rc,dd , d1,d2,d3 , tt
!             real(kind=real64)   ::      ke_before,ke_after
!
!             real(kind=real64),dimension(3)          ::      vv,x0
!             real(kind=real32),dimension(:),allocatable          ::      bb,work
!             integer,dimension(:),allocatable                    ::      ipiv
!             real(kind=real32),dimension(:,:),allocatable        ::      hess
!             real(kind=real64),dimension(3)                      ::      ubar
!
!
!             rc = getCutoff(this%eam)
!             nn = getnNeighMax(this%lc3d)
!             allocate(dx(3,nn))
!             allocate(dr(nn))
!             allocate(id(nn))
!             ke_before = dot_product(this%v(1:3,atom),this%v(1:3,atom))/(2*this%imass)
!             ok = .false.
!
!
!             call neighbourList_long( this%lc3d,this%x(1:3,atom),R0+rc, nSphere )
!             allocate(dr2(nSphere))
!             allocate(id2(nSphere))
!             call neighbourList_long( this%lc3d,this%x(1:3,atom),R0+rc, nSphere,id2,dr2 )
!
!           !  print *,"number of atoms in sphere ",nSphere
!
!             allocate(indx(this%nAtoms)) ; indx(:) = -1
!             !allocate(x0(3,nSphere))                         !   store a copy of the old atom positions in order to recover the old state.
!
!         !---    identify the atoms in the core region
!             nCore = 0
!             do k2 = 1,nSphere
!                 if ( dr2(k2) < R0*R0 ) then
!                     ii = id2(k2)
!                     nCore = nCore + 1
!                     indx(ii) = nCore
!                     if (ii==atom) coreatom = nCore
!                    ! x0(1:3,nCore) = this%x(1:3,ii)
!                 end if
!             end do
!             where (indx==-1)
!                 indx = nCore+1
!             end where
!            ! print *,"core atom count ",nCore," central atom ",coreatom
!
!
!
!         !---    compute energy in sphere region
!             allocate(force(3,nCore+1))
!             force = 0.0d0
!             e0 = 0.0d0
!             do k2 = 1,nSphere
!
!                 call neighbourList( this%lc3d,this%x(1:3,id2(k2)),rc, nn,id,dx,dr )
!
!             !   find self as a neighbour, and remove from neighbour list
!                 jj = minloc(dr(1:nn),dim=1)
!                 do kk = jj+1,nn
!                     id(kk-1) = id(kk) ; dx(:,kk-1) = dx(:,kk) ; dr(kk-1) = dr(kk)
!                 end do
!                 nn = nn - 1
!
!             !   check which core atom is referenced, which neighbour
!                 do kk = 1,nn
!                     id(kk) = indx( id(kk) )             !   either pointing to a core atom or to nCore+1
!
!                 !   convert square distance and vector to distance and normal vector
!                     dr(kk) = sqrt(dr(kk))
!                     dx(:,kk) = dx(:,kk) / dr(kk)
!                 end do
!
!                 call addForce( this%eam, indx(k2) ,nn,dr(1:nn),dx(1:3,1:nn),id(1:nn), force )
!                 e0 = e0 + getEnergy(this%eam, nn,dr(1:nn) )
!             end do
!            ! print *,"sphere energy ",e0
!
!
!
!
!         !---    select random velocity vector with Ekick energy
!             do
!                 call random_number(vv)
!                 vv = 2*vv - 1
!                 dd = vv(1)*vv(1) + vv(2)*vv(2) + vv(3)*vv(3)
!                 if (dd*(1-dd)>0) exit
!             end do
!             dd = sqrt( 2*Ekick*this%imass/dd )
!             vv = vv * dd
!         !   print *,"kick energy ",Ekick," velocity ",vv," KE ",dot_product(vv,vv)/(2*this%imass)
!
!         !---    displace central atom by distance deltax
!
!
!             tt = deltax/norm2(vv)
!             x0(1:3) = this%x(1:3,atom)
!             call cut(this%lc3d,atom,this%x(1:3,atom))
!             this%x(1:3,atom) = x0(1:3)+vv(1:3)*tt
!             call add(this%lc3d,atom,this%x(1:3,atom))
!
!
!
!
!
!
!         !---    compute hessian and force with central atom displaced
!
!
!
!
!             allocate(hess(3*(nCore+1),3*(nCore+1)))
!             allocate(f1(3,nCore+1))
!             f1 = 0
!             e1 = 0
!             hess = 0
!             do k2 = 1,nSphere
!
!                 call neighbourList( this%lc3d,this%x(1:3,id2(k2)),rc, nn,id,dx,dr )
!
!             !   find self as a neighbour, and remove from neighbour list
!                 jj = minloc(dr(1:nn),dim=1)
!                 do kk = jj+1,nn
!                     id(kk-1) = id(kk) ; dx(:,kk-1) = dx(:,kk) ; dr(kk-1) = dr(kk)
!                 end do
!                 nn = nn - 1
!
!             !   check which core atom is referenced, which neighbour
!                 do kk = 1,nn
!                     id(kk) = indx( id(kk) )             !   either pointing to a core atom or to nCore+1
!                     dr(kk) = sqrt(dr(kk))
!                     dx(:,kk) = dx(:,kk) / dr(kk)
!                 end do
!
!                 call addForce(this%eam, indx(k2),nn,dr(1:nn),dx(1:3,1:nn),id(1:nn), f1 )
!                 call addHessian(this%eam, indx(k2),nn,dr(1:nn),dx(1:3,1:nn),id(1:nn), hess )
!                 e1 = e1 + getEnergy(this%eam, nn,dr(1:nn) )
!             end do
!        !    print *,"potential energy before ",e0
!        !    print *,"potential energy after  ",e1
!        !    print *,"number of atoms in sphere ",nSphere
!        !    print *,"number of atoms in core   ",nCore
!
!             if (e1-e0 > 10*Ekick) then
!
!             !---    replace atom in original position
!                 call cut(this%lc3d,atom,this%x(1:3,atom))
!                 this%x(1:3,atom) = x0(1:3)
!                 call add(this%lc3d,atom,this%x(1:3,atom))
!          !       print *,""
!                 return
!             end if
!
!
!
!     !---    construct Lagrange multipler problem hess uu = bb
!             allocate(bb(3*nCore+2))
!             allocate(work(66*(3*nCore+2)))
!             allocate(ipiv(3*nCore+2))
!             bb(:) = 0.0d0
!             hess(:,3*nCore+1:3*nCore+3) = 0.0d0
!             do k2 = 1,nCore
!                 bb(k2*3-2:k2*3) =  real( f1(:,k2)/tt,kind=real32 )
!                 hess(3*k2-2:3*k2,3*nCore+2) = real( force(1:3,k2),kind=real32 )
!             end do
!             hess(3*coreAtom-2:3*coreAtom,3*nCore+1) = real( vv(1:3),kind=real32 )
!
!             call SSYSV( "U",3*nCore+2,1,hess,3*nCore+3,ipiv,bb,3*nCore+2,work,size(work),ii )      !   note: LDA = 3*nCore+3
!          !   print *,"SSYSV returns info = ",ii
!
!         !---    check constraints
!             bb(3*nCore+1:) = 0.0d0
!             ke_after = sum( bb*bb )/(2*this%imass)
!       !      print *,"ke in u ",ke_after," eV"
!
!          !  if (ke_after > Ekick/2) then
!          !  !---    replace atom in original position and return
!          !      call cut(this%lc3d,atom,this%x(1:3,atom))
!          !      this%x(1:3,atom) = x0(1:3)
!          !      call add(this%lc3d,atom,this%x(1:3,atom))
!        ! !       print *,""
!          !      return
!          !  end if
!          !
!
!             ubar = 0.0d0 ; d1 = 0.0d0 ; d2 = 0.0d0 ; d3 = 0.0d0
!             do k2 = 1,nSphere
!                 ii = id2(k2)
!                 jj = indx(ii)
!                 if (jj <= nCore) d1 = d1 + dot_product( bb(jj*3-2:jj*3),this%v(1:3,ii) )
!             end do
!             do k2 = 1,nCore
!                 ubar = ubar + bb(k2*3-2:k2*3)
!                 d2 = d2 + dot_product( bb(k2*3-2:k2*3),f1(:,k2) )
!                 d3 = d3 + dot_product( bb(k2*3-2:k2*3),force(:,k2) )
!             end do
!             ubar = ubar / nCore
!         !    print *,"<u>     ",ubar
!         !    print *,"u.v     ",d1
!         !    print *,"u.f     ",d2
!         !    print *,"u.f0    ",d3
!
!
!
!             bb(coreatom*3-2:coreatom*3) = real( bb(coreatom*3-2:coreatom*3) + vv(1:3),kind=real32 )
!             ke_before = 0.0d0
!             ke_after = 0.0d0
!             do k2 = 1,nSphere
!                 ii = id2(k2)
!                 jj = indx(ii)
!                 if (jj <= nCore) then
!                     ke_before = ke_before + this%v(1,ii)*this%v(1,ii) + this%v(2,ii)*this%v(2,ii) + this%v(3,ii)*this%v(3,ii)
!                     this%v(1:3,ii) = this%v(1:3,ii) + bb(jj*3-2:jj*3)
!                     ke_after = ke_after + this%v(1,ii)*this%v(1,ii) + this%v(2,ii)*this%v(2,ii) + this%v(3,ii)*this%v(3,ii)
!                 end if
!             end do
!             ke_before = ke_before /(2*this%imass)
!             ke_after = ke_after /(2*this%imass)
!       !      print *,"ke_before in core ",ke_before
!       !      print *,"ke_after in core  ",ke_after
!             if (ke_after - ke_before > 2*Ekick) then
!             !---    replace atom in original position and return
!                 call cut(this%lc3d,atom,this%x(1:3,atom))
!                 this%x(1:3,atom) = x0(1:3)
!                 call add(this%lc3d,atom,this%x(1:3,atom))
!                 ! print *,""
!                 return
!             end if
!
!
!         !---    replace atom in original position
!             call cut(this%lc3d,atom,this%x(1:3,atom))
!             this%x(1:3,atom) = x0(1:3)
!             call add(this%lc3d,atom,this%x(1:3,atom))
!             print *,"kick atom ",atom," at ",x0," direction ",vv
!       !      print *,""
!             ok = .true.
!             return
!         end subroutine kickAtom_old



        subroutine wignerSeitz( this,x_ref, occupation )
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      compute a wigner-seitz analysis using reference positions
            type(MD),intent(in)                             ::      this
            real(kind=real64),dimension(:,:),intent(in)     ::      x_ref
            integer,dimension(:),intent(out)                ::      occupation

            type(LinkCell3d),save                   ::      lc3d_ref
            logical,save                            ::      firstcall = .true.



            integer             ::      nSites_ref
            integer             ::      ii,jj
            real(kind=real64)   ::      dr2

            if (firstcall) then
                call clone( lc3d_ref,this%lc3d )    !   deep copy with allocate of approximately correct storage size
                call clear( lc3d_ref )              !   empty atoms in MD
                call add( lc3d_ref,x_ref )              !   add reference atoms instead
                firstcall = .false.
            end if

            nSites_ref = size(x_ref,dim=2)
            occupation = 0
            do ii = 1,this%nAtoms
                call nearestNeighbour( lc3d_ref,this%x(1:3,ii),jj,dr2,self=.true. )
                occupation(jj) = occupation(jj) + 1
            end do


            return
        end subroutine wignerSeitz

        subroutine countPointDefectMovements( occupation1,occupation2 , nv,ni )
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      given occupation1 at time 1 and occupation2 at time 2,
    !*      return the number of vacancy and interstitial movements.
            integer,dimension(:),intent(in)             ::      occupation1
            integer,dimension(:),intent(in)             ::      occupation2
            integer,intent(out)                         ::      nv,ni

            integer     ::      ii
            integer     ::      nSites_ref

            nSites_ref = size(occupation1)

            nv = 0 ; ni = 0
            do ii = 1,nSites_ref

                if (occupation1(ii) == occupation2(ii)) cycle       !   no change.

                if (occupation1(ii) == 0) then
                !   before has a vacancy here
                    nv = nv + 1
                else if (occupation1(ii) > max(1,occupation2(ii)) ) then
                !   before has one or more interstitials
                    ni = ni + occupation1(ii)-occupation2(ii)
                else
                !   this would be double counting movements.
                end if

            end do

            return
        end subroutine countPointDefectMovements


        subroutine countPointDefects( occupation , nv,ni )
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      given occupation
    !*      return the number of vacancy and interstitial
            integer,dimension(:),intent(in)             ::      occupation
            integer,intent(out)                         ::      nv,ni

            integer     ::      ii
            integer     ::      nSites_ref

            nSites_ref = size(occupation)

            nv = 0 ; ni = 0
            do ii = 1,nSites_ref

                if (occupation(ii) == 0) then
                !   has a vacancy here
                    nv = nv + 1
                else
                    ni = ni + occupation(ii)-1
                end if

            end do

            return
        end subroutine countPointDefects

        subroutine computeEnergyDerivatives( this,e,findF,findH )
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
            type(MD),intent(inout)                                          ::      this
            real(kind=real64),intent(out)                                   ::      e
            logical,intent(in)                                              ::      findF,findH

            integer             ::      ii,nn,kk

            real(kind=real64),dimension(:,:),allocatable        ::      dx
            real(kind=real64),dimension(:),allocatable          ::      dr
            integer,dimension(:),allocatable                    ::      id,at

            integer             ::      LL                      !   bandwidth
            real(kind=real64)   ::      vv,rc

            rc = getCutoff(this%EAM)                            !   potential cutoff
            vv = (4*3.141592654d0/3) * (8*rc*rc*rc)             !   sphere radius with 2x potential cutoff
            LL = estimateMaxCountPerVolume(this%lc3d,vv)        !   99.9% confident will have less than equal to this number atoms in sphere

            if (findH .and. .not. associated(this%indx)) then
                allocate(this%hess(1:3,1:3,1:LL,1:this%nAtoms))
                allocate(this%indx(0:LL,1:this%nAtoms))
            end if

            rc = getCutoff(this%eam)
            nn = getnNeighMax(this%lc3d)
            allocate(at(0:nn))
            allocate(id(nn))
            if (findF.or.findH) then
                allocate(dx(3,nn))


            end if
            allocate(dr(nn))

            if (findH) then
                this%indx(0,1:this%nAtoms) = 0                           !   this is the number of neighbours. Don't need to set what neighbour it is or what the element is otherwise
                this%hess(:,:,:,:) = 0.0


            !---    test: can I make things faster by storing the short range connectivity?
                do ii = 1,this%nAtoms
                    call neighbourList( this%lc3d,this%x(1:3,ii),rc, nn,id )
                    this%indx(0,ii) = nn
                    this%indx(1:nn,ii) = id(1:nn)
                end do

            end if





            if (findF) this%f(:,:) = 0.0d0
            e = 0.0d0

            do ii = 1,this%nAtoms

                if (findF.or.findH) then
                    call neighbourList( this%lc3d,this%x(1:3,ii),rc, nn,id,dx,dr )

                !   find self as a neighbour, and remove from neighbour list
                    kk = minloc(dr(1:nn),dim=1)            !   closest neighbour to atom ii is itself.
                    id(kk) = id(nn)                         !   replace neighbour list entry jj with the last entry
                    dx(1:3,kk) = dx(1:3,nn)
                    dr(kk) = dr(nn)
                    nn = nn - 1
                    at(0) = this%at(ii)
                    do kk = 1,nn
                        dr(kk) = sqrt(dr(kk))
                        dx(1:3,kk) = dx(1:3,kk) / dr(kk)
                        at(kk) = this%at(id(kk))
                    end do

                    if (findF) call addForce( this%eam, ii,nn,at,dr,dx,id, this%f )
                    if (findH) call addHessian( this%eam, ii,nn,at,dr,dx,id, this%hess,this%indx )
                else

                    call neighbourList( this%lc3d,this%x(1:3,ii),rc, nn,id,dr )


                !   find self as a neighbour, and remove from neighbour list
                    kk = minloc(dr(1:nn),dim=1)            !   closest neighbour to atom ii is itself.
                    dr(kk) = dr(nn)
                    nn = nn - 1
                    at(0) = this%at(ii)
                    do kk = 1,nn
                        dr(kk) = sqrt(dr(kk))
                        at(kk) = this%at(id(kk))
                    end do
                end if


                e = e + getEnergy( this%eam, nn,at,dr )

            end do

            return
        end subroutine computeEnergyDerivatives

        subroutine computeEnergyDerivativesEinstein( this,nAtoms,x,t,e,f,h, core )
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      compute force ( in this%f ) and einstein hessian
            type(MD),intent(inout)                                          ::      this
            integer,intent(in)                                              ::      nAtoms
            real(kind=real64),dimension(:,:),intent(in)                     ::      x
            integer,dimension(:),intent(in)                                 ::      t
            real(kind=real64),intent(out)                                   ::      e
            real(kind=real64),dimension(:,:),intent(out)                    ::      f
            real(kind=real32),dimension(:,:,:),intent(inout)                ::      h
            logical,dimension(:),intent(in),optional                        ::      core
            integer             ::      ii,nn,kk

            real(kind=real64),dimension(:,:),allocatable        ::      dx
            real(kind=real64),dimension(:),allocatable          ::      dr
            integer,dimension(:),allocatable                    ::      id,at
            real(kind=real64)                                   ::      rc


            rc = getCutoff(this%eam)
            nn = getnNeighMax(this%lc3d)
            allocate(at(0:nn))
            allocate(id(nn))
            allocate(dx(3,nn))
            allocate(dr(nn))

            f = 0.0d0
            h = 0.0d0
            e = 0.0d0

            do ii = 1,nAtoms
            
                if (present(core)) then
                    if (.not. core(ii)) cycle
                end if


                call neighbourList( this%lc3d,x(1:3,ii),rc, nn,id,dx,dr )

            !   find self as a neighbour, and remove from neighbour list
                do kk = 1,nn
                    if (id(kk)==ii) then
                        id(kk) = id(nn)                         !   replace neighbour list entry jj with the last entry
                        dx(1:3,kk) = dx(1:3,nn)
                        dr(kk) = dr(nn)        
                        exit
                    end if
                end do
!                kk = minloc(dr(1:nn),dim=1)            !   closest neighbour to atom ii is itself.
                
                nn = nn - 1
                at(0) = t(ii)
                do kk = 1,nn
                    dr(kk) = sqrt(dr(kk))
                    dx(1:3,kk) = dx(1:3,kk) / dr(kk)
                    at(kk) = t(id(kk))
                end do

                call addForce( this%eam, ii,nn,at,dr,dx,id, f )
                call HessianEinstein( this%eam, ii,nn,at,dr,dx,id, h )
                e = e + getEnergy( this%eam, nn,at,dr )

            end do

            return
        end subroutine computeEnergyDerivativesEinstein


        subroutine computeEnergyFnVecLength( this,nAtoms,x,y,t,alpha,e ,dmax, core,yhat)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      compute the energy at positions x + alpha y
    !*      with alpha being an array of lengths
    !*      under the assumption that the neighbour list is unchanged for all possible alpha
            type(MD),intent(in)                                             ::      this
            integer,intent(in)                                              ::      nAtoms
            real(kind=real64),dimension(:,:),intent(in)                     ::      x,y
            integer,dimension(:),intent(in)                                 ::      t
            real(kind=real64),dimension(:),intent(in)                       ::      alpha
            real(kind=real64),dimension(:),intent(out)                      ::      e
            real(kind=real64),intent(in)                                    ::      dmax      !   max displacment scale
            logical,dimension(:),intent(in),optional                        ::      core
            real(kind=real64),dimension(:,:),intent(in),optional            ::      yhat
            integer             ::      ii,nn,kk , jj

            real(kind=real64),dimension(:,:),allocatable        ::      dx0,dy,dyhat
            real(kind=real64),dimension(:),allocatable          ::      dr,beta
            real(kind=real64),dimension(3)                      ::      dx 
            integer,dimension(:),allocatable                    ::      id,at
            real(kind=real64)                                   ::      rc,d2r,idmax,dd 
            real(kind=real64),dimension(:),allocatable          ::      modx

            rc = getCutoff(this%eam) * 1.1d0        !   a little extra buffer zone is prudent
            nn = getnNeighMax(this%lc3d)
            allocate(at(0:nn))
            allocate(id(nn))
            allocate(dx0(3,nn))
            allocate(dy(3,0:nn))
            allocate(dr(nn))
            allocate(modx(0:nn))
            
            if (present(yhat)) then
            !   compute beta(jj) = F[alpha(jj)*y].yhat 
                allocate(beta(size(alpha)))
                allocate(dyhat(3,nn))
                beta = 0.0d0
                do ii = 1,nAtoms    
                
                    if (present(core)) then
                        if (.not. core(ii)) cycle
                    end if
                    
                    dx = y(1:3,ii)
                    dd = norm2(dx)      
                    do jj = 1,size(alpha)
                        dx = y(1:3,ii) /( 1 + alpha(jj)*dd )     !   scale displacement vector to range dmax                            
                        beta(jj) = beta(jj) + dx(1)*yhat(1,ii) + dx(2)*yhat(2,ii) + dx(3)*yhat(3,ii)
                    end do
                    
                end do
                do jj = 1,size(alpha)
                    beta(jj) = beta(jj)*alpha(jj)
                end do
            end if

            idmax = 1/dmax
            e = 0.0d0

            
            do ii = 1,nAtoms

                if (present(core)) then
                    if (.not. core(ii)) cycle
                end if


                call neighbourList( this%lc3d,x(1:3,ii),rc, nn,id,dx0 )      !   neighbour list without displacement

            !   find self as a neighbour, and remove from neighbour list
                do kk = 1,nn
                     d2r = dx0(1,kk)*dx0(1,kk) + dx0(2,kk)*dx0(2,kk) + dx0(3,kk)*dx0(3,kk)
                     if (d2r < 1.0d-8) then
                         id(kk) = id(nn)
                         dx0(1:3,kk) = dx0(1:3,nn)
                         nn = nn - 1
                         exit
                     end if
                end do

                if (present(yhat)) then
                !   store arrays for atom i
                    at(0) = t(ii)
                    dy(1:3,0) = y(1:3,ii)
                    modx(0) = norm2(dy(1:3,0))*idmax                    
                    do kk = 1,nn
                        jj = id(kk)
                        at(kk) = t(jj)
                        dy(1:3,kk) = y(1:3,jj)  
                        modx(kk) = norm2(dy(1:3,kk) )*idmax
                        dyhat(1:3,kk) = yhat(1:3,jj)  
                    end do
                else
                !   store arrays for atom i
                    at(0) = t(ii)
                    dy(1:3,0) = y(1:3,ii)
                    modx(0) = norm2(dy(1:3,0))*idmax                    
                    do kk = 1,nn
                        jj = id(kk)
                        at(kk) = t(jj)
                        dy(1:3,kk) = y(1:3,jj)  
                        modx(kk) = norm2(dy(1:3,kk) )*idmax                         
                    end do
                end if
                
            !   now consider displacements
                do jj = 1,size(alpha)

                    if (present(yhat)) then
                    !---    find displaced positions of atoms
                        do kk = 1,nn     
                        
                            dx = alpha(jj)*dy(1:3,kk) / ( 1 + alpha(jj)*modx(kk) )      & 
                               - alpha(jj)*dy(1:3,0 ) / ( 1 + alpha(jj)*modx(0 ) )      &
                               + dx0(1:3,kk)                                            &
                               - beta(jj)*dyhat(1:3,kk)
                               
                            dr(kk) = norm2(dx)
                        end do
                        
                    else    
                    !---    find displaced positions of atoms
                        do kk = 1,nn                            
                            dx = alpha(jj)*dy(1:3,kk) / ( 1 + alpha(jj)*modx(kk) )      & 
                               - alpha(jj)*dy(1:3,0 ) / ( 1 + alpha(jj)*modx(0 ) )      &
                               + dx0(1:3,kk)                                            
                               
                            dr(kk) = norm2(dx)
                        end do
                    end if

                !---    compute energy
                    e(jj) = e(jj) + getEnergy( this%eam, nn,at,dr )


                end do

            end do

            return
        end subroutine computeEnergyFnVecLength


        subroutine CGrelax( this,ebest,disp_scale )
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      relax the atoms to a local minimum

            type(MD),intent(inout)          ::      this
            real(kind=real64),intent(out)   ::      ebest
            real(kind=real64),intent(in),optional   ::      disp_scale
            
            integer             ::      step  , loop
            real(kind=real32),dimension(:,:,:),allocatable,save     ::      hessEinstein
            logical,save                                            ::      firstCall = .true.
            real(kind=real32)                   ::      dd
            real(kind=real32),dimension(3)      ::      h1,h2,h3,gg
            real(kind=real64),dimension(0:6),parameter      ::      ALPHA = (/ 0.0d0,-0.1d0,0.1d0,0.3d0,0.6d0,1.0d0,2.0d0 /)
            real(kind=real64),dimension(0:6)                ::      ee
            real(kind=real64)   ::      f2,eps,eold,rc,cc
            real(kind=real64)   ::      displacement_scale
            real(kind=real64),parameter     ::      MAXEPS = 1.0d-6
            integer             ::      ii

        !---    CG works well if we are near a harmonic minimum, and badly if we're not.
        !       start with a couple of steepestish descentish steps, using Einstein oscillator approx to turn force into displacement
            if (firstCall) then
                allocate(hessEinstein(3,3,this%nAtoms))
                firstCall = .false.

            end if
            rc = getCutoff(this%eam)
            displacement_scale = rc/4
            if (present(disp_scale)) displacement_scale = disp_scale
            eps = MAXEPS
            eold = huge(1.0)
            do loop = 1,10
              do step = 1,3
              
              
                !print *,"loop,step ",loop,step," SD atom 17039,17539 at ",this%x(1:3,17039) ,this%x(1:3,17539)

                call computeEnergyDerivativesEinstein( this,this%nAtoms,this%x,this%at,ee(0),this%f,hessEinstein )
                !   now E = E0 + E' dx + 1/2 dx. E" dx
                !         = E0 + E' dx + 1/2 dx. (E_1)" dx + 1/2 dx. (E_2)" dx          , where E_1 is Einstein part
                !   dE/d dx = E' + [ (E_1)" + (E_2)" ) dx
                !         = E' + (E_1)"  dx                                             , in approximation E_1 >> E_2
                !         = 0,
                !         E_1 dx = f
                !
                !   note that I can solve this by finding the cross products of the 3x3 matrix E_1(i) . Don't need lapack :)
                do ii = 1,this%nAtoms
                    h1(1:3) = hessEinstein(1:3,1,ii)
                    h2(1:3) = hessEinstein(1:3,2,ii)
                    h3(1:3) = hessEinstein(1:3,3,ii)
                    gg(1:3) = crossProduct( h2,h3 )
                    dd = h1(1)*gg(1) + h1(2)*gg(2) + h1(3)*gg(3)
                    if (abs(dd)>1.0d-16) dd = 1/dd

                    this%v(1,ii) = dd*( gg(1)*this%f(1,ii) + gg(2)*this%f(2,ii) + gg(3)*this%f(3,ii) )
                    gg(1:3) = crossProduct( h3,h1 )
                    this%v(2,ii) = dd*( gg(1)*this%f(1,ii) + gg(2)*this%f(2,ii) + gg(3)*this%f(3,ii) )
                    gg(1:3) = crossProduct( h1,h2 )
                    this%v(3,ii) = dd*( gg(1)*this%f(1,ii) + gg(2)*this%f(2,ii) + gg(3)*this%f(3,ii) )
                end do




                !   now find the energies as a function of the magnitude of this displacement vector yy
!                call scaleDisplacementVector( this%v,1.0d0,rc/4 )
                call computeEnergyFnVecLength( this,this%nAtoms,this%x,this%v,this%at,alpha(1:),ee(1:),displacement_scale )


                !   where is the lowest energy?
                ii = minloc( ee,dim=1 ) - 1             !   note: -1 because want answer 0:

                !   compute the best ( quadratic approx ) displacement
                ii = min( size(ALPHA)-2,max(1,ii) )          !   because I want to look at range i-1:i+1 for quadratic
                cc = ALPHA(ii-1)*(ee(ii)-ee(ii+1)) + ALPHA(ii)*(ee(ii+1)-ee(ii-1)) + ALPHA(ii+1)*(ee(ii-1)-ee(ii))
                if (abs(cc)>1.0d-16) cc = 1/(2*cc)
                cc = ( ALPHA(ii-1)*ALPHA(ii-1)*(ee(ii)-ee(ii+1)) + ALPHA(ii)*ALPHA(ii)*(ee(ii+1)-ee(ii-1)) + ALPHA(ii+1)*ALPHA(ii+1)*(ee(ii-1)-ee(ii)) )*cc

                !
                
                !write (*,fmt='(a2,i4,100f16.6)') "SD",step,ee(1),ee(0),ee(2:),dd
                call clear(this%lc3d)
                
              !  do ii = 1,this%nAtoms
              !      gg(1:3) = dd*this%v(1:3,ii)
              !      f2 = norm2(gg)
              !      f2 = 1/(1 + 4*f2/rc)
              !      this%x(1:3,ii) = this%x(1:3,ii) + f2*gg(1:3)
              !  end do
              
                                                                                           
                call scaleDisplacementVector( this%v,cc,displacement_scale )
                this%x(1:3,1:this%nAtoms) = this%x(1:3,1:this%nAtoms) + this%v(1:3,1:this%nAtoms)
                                 
                !print *,"loop,step ",loop,step," SD atom 17039,17539 v ",this%v(1:3,17039) ,this%v(1:3,17539)
                                                                          
                !this%x(1:3,1:this%nAtoms) = this%x(1:3,1:this%nAtoms) + dd*this%v(1:3,1:this%nAtoms)
                call add(this%lc3d,this%x)
                
                
                

            end do

        !---    now try some real CG steps

            do step = 1,1

                call computeEnergyDerivatives( this,ee(0),findF=.true.,findH=.true. )
                f2 = dot_product( this%f(1,:),this%f(1,:) )+dot_product( this%f(2,:),this%f(2,:) )+dot_product( this%f(3,:),this%f(3,:) )

                !print *,"simpleMD::CGrelax info - step ",step," energy ",ee ," |f^2| ",f2
                this%v(1:3,1:this%nAtoms) = 0.0d0

            !---    solve E(x+v) = E(x) - f.v + 1/2 v.Dv
            !                D v = f
                eps = min(eps,MAXEPS)

                call conjgrad( this%hess,this%indx, this%v,this%f , eps , restartable = .true.)
!                 print *,"simpleMD::CGrelax info - step ",ii," energy ",eold," -> ",ee," err ",eps
!                print *,dot_product( this%v(1,:),this%v(1,:) )+dot_product( this%v(2,:),this%v(2,:) )+dot_product( this%v(3,:),this%v(3,:) )    &
!                       ,dot_product( this%f(1,:),this%f(1,:) )+dot_product( this%f(2,:),this%f(2,:) )+dot_product( this%f(3,:),this%f(3,:) )
!
                !call scaleDisplacementVector( this%v,1.0d0,rc/4 )   
                call computeEnergyFnVecLength( this,this%nAtoms,this%x,this%v,this%at,alpha(1:),ee(1:),displacement_scale )

                !   where is the lowest energy?
                ii = minloc( ee,dim=1 ) - 1             !   note: -1 because want answer 0:
                ebest = ee(ii)

                !   compute the best ( quadratic approx ) displacement
                ii = min( size(ALPHA)-2,max(1,ii) )          !   because I want to look at range i-1:i+1 for quadratic
                dd = real( ALPHA(ii-1)*(ee(ii)-ee(ii+1)) + ALPHA(ii)*(ee(ii+1)-ee(ii-1)) + ALPHA(ii+1)*(ee(ii-1)-ee(ii)) , kind=real32 ) 
                if (abs(dd)>1.0d-16) dd = 1/(2*dd)
                dd = real( ALPHA(ii-1)*ALPHA(ii-1)*(ee(ii)-ee(ii+1)) + ALPHA(ii)*ALPHA(ii)*(ee(ii+1)-ee(ii-1)) + ALPHA(ii+1)*ALPHA(ii+1)*(ee(ii-1)-ee(ii)) , kind=real32 )*dd

                !write (*,fmt='(a2,i4,100f16.6)') "CG",step,ee(1),ee(0),ee(2:),dd , eps,f2
                
                !print *,"CG atom 2612 at ",this%x(1:3,2612)," -> ",this%x(1:3,2612) + dd*this%v(1:3,2612)
                
                call clear(this%lc3d)
                !this%x(1:3,1:this%nAtoms) = this%x(1:3,1:this%nAtoms) + dd*this%v(1:3,1:this%nAtoms)
                
                
             !  do ii = 1,this%nAtoms
             !      gg(1:3) = dd*this%v(1:3,ii)
             !      f2 = norm2(gg)
             !      f2 = 1/(1 + 4*f2/rc)
             !      this%x(1:3,ii) = this%x(1:3,ii) + f2*gg(1:3)
             !  end do
                call scaleDisplacementVector( this%v,real(dd,kind=real64),displacement_scale )                                                     
                this%x(1:3,1:this%nAtoms) = this%x(1:3,1:this%nAtoms) + this%v(1:3,1:this%nAtoms)
                                                                                                                 
                call add(this%lc3d,this%x)



            end do

            if ( abs(ebest - eold) < 1.0d-3 ) exit
            eold = ebest

            end do


            this%v(1:3,1:this%nAtoms) = 0.0d0

            return

        contains
    !---^^^^^^^^

            pure function crossProduct( a,b ) result( c )
        !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
                real(kind=real32),dimension(3),intent(in)       ::      a,b
                real(kind=real32),dimension(3)                  ::      c
                c(1) = a(2)*b(3) - a(3)*b(2)
                c(2) = a(3)*b(1) - a(1)*b(3)
                c(3) = a(1)*b(2) - a(2)*b(1)
                return
            end function crossProduct

        end subroutine CGrelax



        function elasticConstants( this ) result( C )
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      construct the  elastic_constants tensor
    !*      assumes the neighbour list is in a good state

            type(MD),intent(in)                     ::      this
            real(kind=real64),dimension(6,6)        ::      C
            real(kind=real64),dimension(3,3,3,3)    ::      elastic_constants
            integer                 ::      ii,kk,nn , aa,bb
            real(kind=real64)       ::      vol
            
            
            real(kind=real64),dimension(:,:),allocatable        ::      dx 
            real(kind=real64),dimension(:),allocatable          ::      dr 
            integer,dimension(:),allocatable                    ::      id,at
            real(kind=real64)                                   ::      rc 
 

            rc = getCutoff(this%eam) 
            nn = getnNeighMax(this%lc3d)
            allocate(at(0:nn))
            allocate(id(nn))
            allocate(dx(3,nn))
            allocate(dr(nn))

            
            

            vol = superCellVolume(getSuper(this%lc3d))

            elastic_constants = 0.0d0

            do ii = 1,this%nAtoms

                call neighbourList( this%lc3d,this%x(1:3,ii),rc, nn,id,dx,dr )

            !   find self as a neighbour, and remove from neighbour list
                do kk = 1,nn
                    if (id(kk)==ii) then
                        id(kk) = id(nn)                         !   replace neighbour list entry jj with the last entry
                        dx(1:3,kk) = dx(1:3,nn)
                        dr(kk) = dr(nn)        
                        exit
                    end if
                end do
!                kk = minloc(dr(1:nn),dim=1)            !   closest neighbour to atom ii is itself.
                
                nn = nn - 1
                at(0) = this%at(ii)
                do kk = 1,nn
                    dr(kk) = sqrt(dr(kk))
                    dx(1:3,kk) = dx(1:3,kk) / dr(kk)
                    at(kk) = this%at(id(kk))
                end do

                
                call addelastic_constants( this%eam,nn,at,dr,dx,elastic_constants )

            end do


            C(1,1) = elastic_constants(1,1,1,1)
            C(2,1) = elastic_constants(2,2,1,1)
            C(3,1) = elastic_constants(3,3,1,1)
            C(4,1) = elastic_constants(2,3,1,1)
            C(5,1) = elastic_constants(3,1,1,1)
            C(6,1) = elastic_constants(1,2,1,1)

            C(2,2) = elastic_constants(2,2,2,2)
            C(3,2) = elastic_constants(3,3,2,2)
            C(4,2) = elastic_constants(2,3,2,2)
            C(5,2) = elastic_constants(3,1,2,2)
            C(6,2) = elastic_constants(1,2,2,2)

            C(3,3) = elastic_constants(3,3,3,3)
            C(4,3) = elastic_constants(2,3,3,3)
            C(5,3) = elastic_constants(3,1,3,3)
            C(6,3) = elastic_constants(1,2,3,3)

            C(4,4) = elastic_constants(2,3,2,3)
            C(5,4) = elastic_constants(3,1,2,3)
            C(6,4) = elastic_constants(1,2,2,3)

            C(5,5) = elastic_constants(3,1,3,1)
            C(6,5) = elastic_constants(3,1,1,2)

            C(6,6) = elastic_constants(1,2,1,2)

            do aa = 1,5
                do bb = aa+1,6
                    C(aa,bb) = C(bb,aa)
                end do
            end do

            C = C / vol


            return
        end function elasticConstants




        function dipoleTensor( this ) result( P )
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      construct the  dipole tensor
    !*      assumes the neighbour list is in a good state

            type(MD),intent(in)                     ::      this
            real(kind=real64),dimension(3,3)        ::      P
            integer                 ::      ii,kk,nn  
            !real(kind=real64)       ::      vol
            
            
            real(kind=real64),dimension(:,:),allocatable        ::      dx 
            real(kind=real64),dimension(:),allocatable          ::      dr 
            integer,dimension(:),allocatable                    ::      id,at
            real(kind=real64)                                   ::      rc 
 

            rc = getCutoff(this%eam) 
            nn = getnNeighMax(this%lc3d)
            allocate(at(0:nn))
            allocate(id(nn))
            allocate(dx(3,nn))
            allocate(dr(nn))

             

            !vol = superCellVolume(getSuper(this%lc3d))

            P = 0.0d0

            do ii = 1,this%nAtoms

                call neighbourList( this%lc3d,this%x(1:3,ii),rc, nn,id,dx,dr )

            !   find self as a neighbour, and remove from neighbour list
                do kk = 1,nn
                    if (id(kk)==ii) then
                        id(kk) = id(nn)                         !   replace neighbour list entry jj with the last entry
                        dx(1:3,kk) = dx(1:3,nn)
                        dr(kk) = dr(nn)        
                        exit
                    end if
                end do

                                
                nn = nn - 1
                at(0) = this%at(ii)
                do kk = 1,nn
                    dr(kk) = sqrt(dr(kk))
                    dx(1:3,kk) = dx(1:3,kk) / dr(kk)
                    at(kk) = this%at(id(kk))
                end do

                call addVirial( this%eam,nn,at,dr,dx, P )

            end do
 

            return
        end function dipoleTensor




        subroutine CGrelax_old( this,ebest )
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
            type(MD),intent(inout)          ::      this
            real(kind=real64),intent(out)   ::      ebest

            real(kind=real64),parameter     ::      MAXEPS = 1.0d-4
            real(kind=real64)               ::      eps,ee,eold,elast,f2,f2old
            integer                         ::      ii,jj,jbest

            real(kind=real64),dimension(:,:),allocatable    ::      x0

            real(kind=real64)     ::      JIGGLE

            elast = huge(1.0)/2
            ebest = elast
            eold = huge(1.0)
            allocate(x0(3,this%nAtoms))
            f2 = 0
            f2old = -1
            JIGGLE = 0.001d0

            do ii = 1,32


                call computeEnergyDerivatives( this,ee,findF=.true.,findH=.true. )
                f2 = dot_product( this%f(1,:),this%f(1,:) )+dot_product( this%f(2,:),this%f(2,:) )+dot_product( this%f(3,:),this%f(3,:) )

                print *,"simpleMD::CGrelax info - step ",ii," energy ",ee ," |f^2| " &
                     ,dot_product( this%f(1,:),this%f(1,:) )+dot_product( this%f(2,:),this%f(2,:) )+dot_product( this%f(3,:),this%f(3,:) )
                this%v(1:3,1:this%nAtoms) = 0.0d0

            !---    solve E(x+v) = E(x) - f.v + 1/2 v.Dv
            !                D v = f
                eps = min(eps,MAXEPS)

                call conjgrad( this%hess,this%indx, this%v,this%f , eps , restartable = .true.)
!                 print *,"simpleMD::CGrelax info - step ",ii," energy ",eold," -> ",ee," err ",eps
!                print *,dot_product( this%v(1,:),this%v(1,:) )+dot_product( this%v(2,:),this%v(2,:) )+dot_product( this%v(3,:),this%v(3,:) )    &
!                       ,dot_product( this%f(1,:),this%f(1,:) )+dot_product( this%f(2,:),this%f(2,:) )+dot_product( this%f(3,:),this%f(3,:) )
!
                x0(1:3,1:this%nAtoms) = this%x(1:3,1:this%nAtoms)
                jbest = 0 ; ebest = ee
                do jj = 1,10
                    this%x(1:3,1:this%nAtoms) = x0(1:3,1:this%nAtoms) + jj*this%v(1:3,1:this%nAtoms)/5
                    call clear(this%lc3d)
                    call add(this%lc3d,this%x)
                    call computeEnergyDerivatives( this,ee,findF=.false.,findH=.false. )
!                     print *,"jj = ",jj,jj*0.2d0,ee
                    if (ee<ebest) then
                        jbest = jj ; ebest = ee
                    else
                        exit
                    end if
                end do

                if (ebest < elast) then
                    this%x(1:3,1:this%nAtoms) = x0(1:3,1:this%nAtoms) + jbest*this%v(1:3,1:this%nAtoms)/5
                else
                    jbest = 0 ; ebest = ee
                    do jj = 1,10
                        this%x(1:3,1:this%nAtoms) = x0(1:3,1:this%nAtoms) + jj*this%f(1:3,1:this%nAtoms)/20
                        call clear(this%lc3d)
                        call add(this%lc3d,this%x)
                        call computeEnergyDerivatives( this,ee,findF=.false.,findH=.false. )
    !                     print *,"jj = ",jj,jj*0.2d0,ee
                        if (ee<ebest) then
                            jbest = jj ; ebest = ee
                        else
                            exit
                        end if
                    end do

                    this%x(1:3,1:this%nAtoms) = x0(1:3,1:this%nAtoms) + jbest*this%f(1:3,1:this%nAtoms)/20

                end if




                call clear(this%lc3d)
                call add(this%lc3d,this%x)

                if ( abs(ebest - eold) < 1.0d-3 ) exit
                eold = elast
                elast = ebest
                if ( f2 > f2old*0.9d0 ) then
                    call random_number(x0)
                    x0 = x0*JIGGLE - JIGGLE/2
                    this%x(1:3,1:this%nAtoms) = this%x(1:3,1:this%nAtoms) + x0(1:3,1:this%nAtoms)
                end if
                f2old = f2
            end do
            this%v(1:3,1:this%nAtoms) = 0.0d0

            return
        end subroutine CGrelax_old

        
        
        

        pure subroutine scaleDisplacementVector( dx,alpha,rc,yhat )
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      scale a displacement vector so its magnitude is no larger than rc
    !*      using the simple switch F[x] = x/(1+x)
     
    
            real(kind=real64),dimension(:,:),intent(inout)          ::      dx
            real(kind=real64),intent(in)                            ::      alpha,rc
            real(kind=real64),dimension(:,:),intent(in),optional    ::      yhat
            real(kind=real64),dimension(3)      ::      xx
            real(kind=real64)                   ::      modx,irc,beta
            integer                             ::      ii
            
            irc = 1/rc
            
            if (present(yhat)) then         
                beta = 0.0d0       
                do ii = 1,size(dx,dim=2)
                    xx(1:3) = alpha*dx(1:3,ii)
                    modx = norm2(xx)
                    modx = 1/( 1 + modx*irc )
                    beta = beta + modx*( xx(1)*yhat(1,ii) + xx(2)*yhat(2,ii) + xx(3)*yhat(3,ii) )
                end do                
                do ii = 1,size(dx,dim=2)
                    xx(1:3) = alpha*dx(1:3,ii) 
                    modx = norm2(xx)
                    modx = 1/( 1 + modx*irc )
                    dx(1:3,ii) = modx*xx(1:3) - beta*yhat(1:3,ii)
                end do    
            else
                do ii = 1,size(dx,dim=2)
                    xx(1:3) = alpha*dx(1:3,ii)
                    modx = norm2(xx)
                    modx = 1/( 1 + modx*irc )
                    dx(1:3,ii) = modx*xx(1:3)
                end do                
            end if
            return
        end subroutine scaleDisplacementVector

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

 


    end module simpleMD

