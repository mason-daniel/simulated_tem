
    program test_Lib_ImagingSpace
!---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

        use Lib_ImagingSpace
        use Lib_AtomSpace
        use Lib_RotationMatrices
        use Lib_XYZFiles
#ifdef MPI
        use mpi_f08
#endif
        use iso_fortran_env
        implicit none


        type(ImagingSpace)  ::          is
        integer,parameter                               ::      MX = 50,MY = 40, MZ = 30
        integer,parameter                               ::      NATOMS = MX*MY*MZ*2
        real(kind=real64),dimension(3,3)                ::      a_super
        real(kind=real64),dimension(3)                  ::      xyzoffset
        real(kind=real64),dimension(3,NATOMS)           ::      x           !   atom positions
        real(kind=real64)                               ::      a           !   grid spacing
        real(kind=real64)                               ::      a0          !   atom lattice const
        real(kind=real64)                               ::      sigma       !   blurring 

        type(AtomSpace)         ::      as
        real(kind=real64),dimension(3,3)        ::      R,eps              !   rotation matrix,strain matrix
        integer                 ::      Nx,Ny,Nz
        integer             ::          mynAtoms
        real(kind=real64),dimension(:,:),pointer        ::      myxt

        integer             ::          ii,jj,np
        integer             ::          ix,iy,iz
        real(kind=real64),dimension(3,9)                ::      xtp  
        real(kind=real64),dimension(:,:),pointer        ::      myxt_tmp
        
        logical             ::          ok
        integer             ::          rank = 0 ,nprocs = 1
        type(XYZFile)       ::          xyz
        character(len=256)  ::          filename
        
#ifdef MPI
        call MPI_INIT(ii)
        call MPI_COMM_SIZE(MPI_COMM_WORLD, nProcs, ii)
        call MPI_COMM_RANK(MPI_COMM_WORLD, rank, ii)        
#endif
        if (rank==0) print *,"test_Lib_ImagingSpace running on ",nProcs,"procs"
        call Lib_ImagingSpace_init_MPI()


        
        
        
        a0 = 3.165d0        !   atom lattice constant
        a = a0/4            !   grid spacing
        sigma = a0/2        !   blurring
        a_super = 0
        a_super(1,1) = MX
        a_super(2,2) = MY
        a_super(3,3) = MZ
        xyzoffset = 0                                                           !   ... with zero offset ...
                                                                      

    !---    place atoms in bcc fractional positions (0:1)
        ii = 0
        do iz = 0,MZ-1
            do iy = 0,MY-1
                do ix = 0,MX-1
                    ii = ii + 1 ; x(:,ii) = (/ real(ix,kind=real64)/MX,real(iy,kind=real64)/MY,real(iz,kind=real64)/MZ /)      
                    ii = ii + 1 ; x(:,ii) = (/ real(ix+0.5d0,kind=real64)/MX,real(iy+0.5d0,kind=real64)/MY,real(iz+0.5d0,kind=real64)/MZ /)    
                end do
            end do
        end do    
         
        a_super = a_super * a0
        eps = 0
        R = RotationMatrix_identity


        ! R = RotationMatrix_ctor( (/0.0d0,1.0d0,0.0d0/), 30*3.14159265930d0/180 )       
        ! eps(1,2) = 0.05d0 ; eps(2,1) = eps(1,2)
        ! eps(1,3) = 0.05d0 ; eps(3,1) = eps(1,3)
        ! eps(2,3) = 0.05d0 ; eps(3,2) = eps(2,3)



        a_super = a_super + matmul(eps,a_super)
        x = rotateVector(a_super,x)

        
        as = AtomSpace_ctor(a,xyzoffset,a_super)
        call setR(as,R)
        call setThickness(as)
        call suggestImagingSpace(as,Nx,Ny,Nz) 
        call setDelta(as,Nx,Ny,Nz)
        if (rank==0) call report(as) 
        if (rank==0) print *,""


        is = ImagingSpace_ctor(a,sigma, getdelta(as),Nx,Ny,Nz)
        if (rank==0) call report(is) 
        if (rank==0) print *,""

        allocate(myxt(3,ceiling(real(NATOMS)/nprocs)))
        mynAtoms = 0
        
        do ii = 1,NATOMS
            call periodicCopies(as,is,x(:,ii),np,xtp)
            do jj = 1,np
                if (inMyCell(is,xtp(:,jj),buffered=.true.)) then
                    mynAtoms = mynAtoms + 1
                    if (mynAtoms > size(myxt,dim=2)) then
                        allocate(myxt_tmp(3,ceiling(size(myxt,dim=2)*1.5)))
                        myxt_tmp(1:3,1:size(myxt,dim=2)) = myxt(:,:) 
                        deallocate(myxt)
                        myxt => myxt_tmp
                    end if
                    myxt(1:3,mynAtoms) = xtp(:,jj)
                end if
            end do
        end do
        
        print *,"rank ",rank," has ",mynAtoms," atoms "

        if (nProcs <= 4) then
            write(filename,fmt='(i6)') rank ; filename = "test.p"//trim(adjustl(filename))//".xyz"
            xyz = XYZFile_ctor(filename)
            call setAtomNames(xyz,(/"int","buf"/))
            call setnAtoms(xyz,mynAtoms)
            call setColumn_Description(xyz,Nx,Ny,Nz,RotationMatrix_identity)
            call setnHeaderLines(xyz,0)
            call setnColumns(xyz,3)
            do ii = 1,mynAtoms
                if (inMyCell(is,myxt(:,ii),buffered=.false.)) then
                    call setAtomType(xyz,ii,1)
                else
                    call setAtomType(xyz,ii,2)
                end if
                call setColumns(xyz,ii,myxt(:,ii))
            end do
            if (rank==0) call report(xyz)
            call output(xyz)
        end if

        !ok = (all(myxt(3,1:mynAtoms)>=0)).and.(all(myxt(3,1:mynAtoms)<=Nz))
        ok = .true.

        call MPI_FINALIZE(ii)
        if (rank==0) then
            print *,""
            if (ok) print *,"PASS"
            print *,"done"
            print *,""
        end if


    end program test_Lib_ImagingSpace

