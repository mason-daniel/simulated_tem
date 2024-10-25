
    program test_Lib_ComputePhaseFactor
!---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
!*      program to test correct working of phase factor code

        use Lib_ComputePhaseFactor
        use Lib_ImagingSpace
        use Lib_AtomSpace
        use Lib_RotationMatrices
        use Lib_XYZFiles
        use Lib_Callipers
#ifdef MPI
        use mpi_f08
#endif
        use iso_fortran_env
        implicit none

        integer,parameter                               ::      MX = 20,MY = 30, MZ = 30
        integer,parameter                               ::      NATOMS = MX*MY*MZ*2
        integer,parameter                               ::      NG = 4
        real(kind=real64),parameter                     ::      PI = 3.14159265390d0


        real(kind=real64),dimension(3,3)                ::      a_super
        real(kind=real64),dimension(3)                  ::      xyzoffset
        real(kind=real64),dimension(3,NATOMS)           ::      r           !   atom positions
        real(kind=real64)                               ::      a           !   grid spacing
        real(kind=real64)                               ::      a0          !   atom lattice constant
        real(kind=real64)                               ::      sigma       !   blurring
        type(AtomSpace)         ::      as
        type(ImagingSpace)  ::          is
        real(kind=real64),dimension(3,3)        ::      RR,eps              !   rotation matrix,strain matrix
        integer                 ::      Nx,Ny,Nz , Nx0,Ny0,Nz0 , dNxy,dNz
        integer             ::          mynAtoms
        real(kind=real64),dimension(:,:),pointer        ::      myrt

        integer             ::          ii,jj,np
        real(kind=real64),dimension(3,9)                ::      rtp  
        real(kind=real64),dimension(:,:),pointer        ::      myrt_tmp

        integer,dimension(2,3)                          ::      bb
        real(kind=real64),dimension(3,NG)               ::      g
        complex(kind=real64),dimension(:,:,:,:),pointer      ::      x       !   (1:nGvec,lbx:ubx,lby:uby,lbz:ubz)
        real(kind=real64),dimension(:,:,:,:,:),pointer    ::      grad_arg_x       !   (3,1:nGvec,lbx:ubx,lby:uby,lbz:ubz)
        real(kind=real64),dimension(:,:,:),pointer      ::      rho
        real(kind=real64),dimension(12)                  ::      dat
        integer                                         ::      ix,iy,iz

        
        logical             ::          ok
        integer             ::          rank = 0 ,nprocs = 1
        type(XYZFile)       ::          xyz
        character(len=256)  ::          filename
        type(Callipers)     ::          timer
        real(kind=real64)   ::          t_phaseField,t_phaseField_tmp
        
#ifdef MPI
        call MPI_INIT(ii)
        call MPI_COMM_SIZE(MPI_COMM_WORLD, nProcs, ii)
        call MPI_COMM_RANK(MPI_COMM_WORLD, rank, ii)        
#endif
        if (rank==0) print *,"test_Lib_ImagingSpace running on ",nProcs,"procs"
        call Lib_ImagingSpace_init_MPI()


        
        a0 = 3.165d0        !   atom lattice constant
        a =  a0/2           !   grid spacing
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
                    ii = ii + 1 ; r(:,ii) = (/ real(ix,kind=real64)/MX,real(iy,kind=real64)/MY,real(iz,kind=real64)/MZ /)      
                    ii = ii + 1 ; r(:,ii) = (/ real(ix+0.5d0,kind=real64)/MX,real(iy+0.5d0,kind=real64)/MY,real(iz+0.5d0,kind=real64)/MZ /)    
                end do
            end do
        end do    
         
        a_super = a_super * a0

        RR = RotationMatrix_identity
        eps = 0
        
        ! more complex test
        eps(1,2) = 0.05d0 ; eps(2,1) = eps(1,2)
        eps(1,3) = 0.05d0 ; eps(3,1) = eps(1,3)
        eps(2,3) = 0.05d0 ; eps(3,2) = eps(2,3)        
        RR = RotationMatrix_ctor( (/0.0d0,1.0d0,0.0d0/), 15*3.14159265930d0/180 )      

        a_super = a_super + matmul(eps,a_super)

    !---    
        r = rotateVector(a_super,r)

        as = AtomSpace_ctor(a,xyzoffset,a_super)
        call setR(as,RR)
        call setThickness(as)
        call suggestImagingSpace(as,Nx0,Ny0,Nz0)                    !   imaging space required if no dispersion angle        
        call suggestImagingSpace(as,Nx,Ny,Nz , theta = 0.25d0)      !   imaging space required with dispersion angle
        dNxy = (Nx-Nx0)/2
        dNz = (Nz-Nz0)/2
        call setDelta(as,Nx,Ny,Nz)
        if (rank==0) call report(as) 
        if (rank==0) print *,""


        is = ImagingSpace_ctor(a,sigma,getdelta(as),Nx,Ny,Nz)
        if (rank==0) call report(is) 
        if (rank==0) print *,""

        allocate(myrt(3,ceiling(real(NATOMS)/nprocs)))
        mynAtoms = 0
        
        do ii = 1,NATOMS
            call periodicCopies(as,is,r(:,ii),np,rtp)
            do jj = 1,np
                if (inMyCell(is,rtp(:,jj),buffered=.true.)) then
                    mynAtoms = mynAtoms + 1
                    if (mynAtoms > size(myrt,dim=2)) then
                        allocate(myrt_tmp(3,ceiling(size(myrt,dim=2)*1.5)))
                        myrt_tmp(1:3,1:size(myrt,dim=2)) = myrt(:,:) 
                        deallocate(myrt)
                        myrt => myrt_tmp
                    end if
                    myrt(1:3,mynAtoms) = rtp(:,jj)
                end if
            end do
        end do
        
        print *,"rank ",rank," has ",mynAtoms," atoms "
 
        g(:,1) = 2*PI*(/ -2.0d0,0.0d0,0.0d0 /)/a
        g(:,2) = 2*PI*(/ 2.0d0,0.0d0,0.0d0 /)/a
        g(:,3) = 2*PI*(/ 1.0d0,1.0d0,0.0d0 /)/a
        g(:,4) = 2*PI*(/ 1.0d0,-1.0d0,0.0d0 /)/a
        g = rotateVector(RR,g)
    
    !---    find region I am responsible for, _excluding_ buffers
        call getBounds(is,bb)

        !allocate(x(NG,bb(1,1):bb(2,1),bb(1,2):bb(2,2),bb(1,3):bb(2,3)))
        allocate(rho(bb(1,1):bb(2,1),bb(1,2):bb(2,2),bb(1,3):bb(2,3)))        
        allocate(grad_arg_x(3,NG,bb(1,1):bb(2,1),bb(1,2):bb(2,2),bb(1,3):bb(2,3)))
        allocate(x(NG,bb(1,1):bb(2,1),bb(1,2):bb(2,2),bb(1,3):bb(2,3)))        
        !print *,"rank ",rank," start computePhaseFactor"

        Lib_ComputePhaseFactor_DBG = .true.
        timer = Callipers_ctor()
        call computePhaseFactor( myNatoms,myrt,g, is,getdelta(as) ,grad_arg_x,rho,x )
        print *,"bounds x ",lbound(grad_arg_x,dim=3),lbound(rho,dim=1),lbound(x,dim=2) , ":" , ubound(grad_arg_x,dim=3),ubound(rho,dim=1),ubound(x,dim=2)
        print *,"bounds y ",lbound(grad_arg_x,dim=4),lbound(rho,dim=2),lbound(x,dim=3) , ":" , ubound(grad_arg_x,dim=4),ubound(rho,dim=2),ubound(x,dim=3)
        print *,"bounds z ",lbound(grad_arg_x,dim=5),lbound(rho,dim=3),lbound(x,dim=4) , ":" , ubound(grad_arg_x,dim=5),ubound(rho,dim=3),ubound(x,dim=4)

        call pause(timer)
        t_phaseField = elapsed(timer)
#ifdef MPI
        t_phaseField_tmp = t_phaseField
        call MPI_REDUCE(t_phaseField_tmp,t_phaseField,1,MPI_DOUBLE_PRECISION,MPI_MAX,0,MPI_COMM_WORLD,ii)
        if (rank==0) print *,"computePhaseFactor time(max) ",t_phaseField 
        call MPI_REDUCE(t_phaseField_tmp,t_phaseField,1,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,ii)
        if (rank==0) print *,"computePhaseFactor time(avg) ",t_phaseField/nProcs
#else
        print *,"computePhaseFactor time ",t_phaseField 
#endif

        if (rank==0) then
            print *,"derivative test ",bb
            print *,"grad_arg_x(g=1,1,1,1)  ",grad_arg_x(:,1,1,1,1)
            !print *,"numerical         ",( x(1,2,1,1)-x(1,0,1,1) )/(2*a),( x(1,1,2,1)-x(1,1,0,1) )/(2*a),( x(1,1,1,2)-x(1,0,1,0) )/(2*a)
        end if

        !print *,"rank ",rank," done computePhaseFactor"

        if (nProcs <= 4) then
            write(filename,fmt='(i6)') rank ; filename = "test.p"//trim(adjustl(filename))//".xyz"
            xyz = XYZFile_ctor(filename)
            call setAtomNames(xyz,(/"int","buf","g.r"/))
            call setnAtoms(xyz,mynAtoms+size(rho,dim=1)*size(rho,dim=2)*size(rho,dim=3))
            call setColumn_Description(xyz,Nx,Ny,Nz,RotationMatrix_identity,":grad_arg(x):R:3:rho:R:1:Re(x):R:1:Im(x):R:1" )
            call setnHeaderLines(xyz,0)
            call setnColumns(xyz,9)
            dat = 0
            do ii = 1,mynAtoms
                if (inMyCell(is,myrt(:,ii),buffered=.false.)) then
                    call setAtomType(xyz,ii,1)
                else
                    call setAtomType(xyz,ii,2)
                end if                
                dat(1:3) = myrt(1:3,ii) 
                call setColumns(xyz,ii,dat)
            end do
            ii = mynAtoms
            do iz = bb(1,3),bb(2,3)
                do iy = bb(1,2),bb(2,2)
                    do ix = bb(1,1),bb(2,1)
                        ii = ii + 1
                        call setAtomType(xyz,ii,3)
                        dat(1:3) = (/ ix,iy,iz /) + 0.5d0              !   +0.5 because x nodes are at (1/2,1/2,1/2) positions
                        dat(4:6) = grad_arg_x(:,1,ix,iy,iz)
                        !dat(5) = aimag( x(2,ix,iy,iz) )
                        dat(7) = rho(ix,iy,iz)
                        dat(8) = real(x(1,ix,iy,iz))
                        dat(9) = aimag(x(1,ix,iy,iz))
                        call setColumns(xyz,ii,dat(1:9))
                    end do
                end do
            end do

            if (rank==0) call report(xyz)
            call output(xyz)


        end if


        ok = .true.

        call MPI_FINALIZE(ii)
        if (rank==0) then
            print *,""
            if (ok) print *,"PASS"
            print *,"done"
            print *,""
        end if

    end program test_Lib_ComputePhaseFactor