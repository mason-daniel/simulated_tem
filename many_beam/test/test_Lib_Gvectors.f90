
    program test_Lib_Gvectors
!---^^^^^^^^^^^^^^^^^^^^^^^^^
!*      simple test program to illustrate correct function of g-vector code
        use Lib_Gvectors
        use Lib_Lattices
        use Lib_RotationMatrices
        use iso_fortran_env
        implicit none


        type(Lattice)                               ::      latt
        real(kind=real64),dimension(3,3)            ::      a_conventional_cell
        integer                                     ::      n
        integer,dimension(:,:),allocatable          ::      hkl,hjkl


        !real(kind=real64),dimension(3,3)            ::      R
        type(Gvectors)                              ::      gv

        integer             ::      ii

    !---    find the g-vectors for bcc
        latt = Lattice_ctor("bcc")
        call report(latt)
        a_conventional_cell = getConventionalCell(latt)
        call permittedReflections( latt,n,hkl,rho_in=2.0d0 )
        print *,""

    !---    construct a g-vector object
        gv = Gvectors_ctor(a_conventional_cell,hkl)
        call report(gv)
        print *,""

        !R = RotationMatrix_ctor( (/0.0d0,1.0d0,0.0d0/), 30*3.14159265930d0/180 )      
        !call report(R)
        !call setR(gv,R)
        !call report(gv)
        !print *,""

        latt = Lattice_ctor("hcp")
        call report(latt)
        a_conventional_cell = getConventionalCell(latt)
        call permittedReflections( latt,n,hkl,rho_in=1.5d0 )
        print *,""

    !---    construct a g-vector object
        allocate(hjkl(4,n))
        do ii = 1,n
            hjkl(:,ii) = nint( MillerToMillerBravais_plane(real(hkl(:,ii),kind=real64)) )
        end do


        gv = Gvectors_ctor(a_conventional_cell,hjkl)
        call report(gv)
        print *,""
 


        print *,""
        print *,"PASS"
        print *,"done"
        print *,""


    end program test_Lib_Gvectors
