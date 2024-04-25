
    program strainXYZ
!---^^^^^^^^^^^^^^^^^
!*      quick code to strain an .xyz box.
!*
!*      Daniel Mason
!*      (c) UKAEA May 2023  
!*

        use Lib_CommandLineArguments
        use Lib_XYZFiles
        use iso_fortran_env
        implicit none

        type(CommandLineArguments)      ::      cla

    !---    input parameters
        character(len=256)              ::      filename = ""           !   input filename
        character(len=256)              ::      outfile = ""            !   output filename
        character(len=8)                ::      VERSION = "1.0.0"



    !---    physically meaningful parameters deduced from input file
        real(kind=real64),dimension(3,3)    ::      a_super , ia_super
        type(XYZFile)                       ::      xyz
        real(kind=real64),dimension(3,3)    ::      strain = 0
        real(kind=real64),dimension(:,:),pointer    ::      colp
        integer                             ::      nAtoms
        real(kind=real64),dimension(3,3)    ::      stretch

    !---    dummy parameter

        logical                             ::      ok
        real(kind=real64),dimension(9)      ::      edat = 0
        real(kind=real64),dimension(3)      ::      xx
        integer,dimension(3)                ::      yy
        integer                             ::      nn,ii
        character(len=1024)                 ::      dummy,coldesc

    !---    read command line arguments
        cla = CommandLineArguments_ctor(10)

        call setProgramDescription( cla, "strainXYZ.exe <filename>" )
        call setProgramVersion( cla, VERSION )

        call get( cla,"f",filename ,LIB_CLA_REQUIRED,"          input .xyz filename" )
        call get( cla,"o",outfile ,LIB_CLA_REQUIRED,"          output .xyz filename" )
        nn = 0
        call get( cla,"e",edat,nn ,LIB_CLA_OPTIONAL,"       strain [e_hydro] or [ exx,eyy,ezz ] or [ exx,eyy,ezz,exy,eyz,ezx ] or [exx,eyx,ezx,exy,...,ezz]" )
        if (nn==1) then
            strain(1,1) = edat(1)
            strain(2,2) = edat(1)
            strain(3,3) = edat(1)
        else if (nn==3) then
            strain(1,1) = edat(1)
            strain(2,2) = edat(2)
            strain(3,3) = edat(3)
        else if (nn==6) then
            strain(1,1) = edat(1)
            strain(2,2) = edat(2)
            strain(3,3) = edat(3)
            strain(1,2) = edat(4) ; strain(2,1) = edat(4)
            strain(2,3) = edat(5) ; strain(3,2) = edat(5)
            strain(3,1) = edat(6) ; strain(1,3) = edat(6)
        else
            strain(1:3,1:3) = reshape( edat,(/3,3/) )
        end if

        call report(cla)
        if (.not. allRequiredArgumentsSet(cla)) stop
        if (hasHelpArgument(cla)) stop
        call delete(cla)


        print *,""
        print *,"strainXYZ.exe v"//trim(VERSION)
        print *,"^^^^^^^^^^^^^^^"//repeat("^",len_trim(VERSION))
        print *,""
        print *,"strain"
        do ii = 1,3
            write (*,fmt='(3f16.8)') strain(ii,:)
        end do
        print *,""


    !---
        print *,""
        print *,"reading input file"
        print *,"^^^^^^^^^^^^^^^^^^"
        print *,""
        xyz = XYZFile_ctor(filename)
        call readHeader(xyz,ok)
        call input(xyz,verbose=.true.)
        call getSupercell(xyz,a_super,ok)
        if (ok) then
            print *,"strainXYZ.exe info - supercell read from file"
            print *,a_super(1,:)
            print *,a_super(2,:)
            print *,a_super(3,:)
        else
            print *,"strainXYZ.exe warning - supercell not read from file"
        end if
        call getColumnsp(xyz,colp)
        nAtoms = getNatoms(xyz)
        call report(xyz)
    !---


    !---    convert strain to stretch
        stretch(:,:) = strain(:,:)
        stretch(1,1) = stretch(1,1) + 1.0d0
        stretch(2,2) = stretch(2,2) + 1.0d0
        stretch(3,3) = stretch(3,3) + 1.0d0


        call inverse3Mat(a_super,ia_super)

    !   change supercell params.
        a_super = matmul( stretch,a_super )
        dummy = getColumn_description(xyz)
        nn = index( dummy,"Properties" )
        call setColumn_description(xyz,a_super)
        if (nn>0) then
            coldesc = getColumn_description(xyz)
            ii = index( dummy,"Properties" )
            coldesc = coldesc(1:ii-1)//dummy(nn:)
            call setColumn_description(xyz,coldesc)
        end if


    !   strain atoms
        do ii = 1,nAtoms
        !   strained atom position
            xx(1:3) = stretch(1:3,1)*colp(1,ii) + stretch(1:3,2)*colp(2,ii) + stretch(1:3,3)*colp(3,ii)
        !   which supercell periodic repeat?
            yy(1:3) = floor( ia_super(1:3,1)*xx(1) + ia_super(1:3,2)*xx(2) + ia_super(1:3,3)*xx(3) )
        !   remove supercell offset
            colp(1:3,ii) = xx(1:3) - a_super(1:3,1)*yy(1) - a_super(1:3,2)*yy(2) - a_super(1:3,3)*yy(3)
        end do


    !---    output result
        call setFilename(xyz,outfile)
        call setLammpsFormat( xyz,outfile( len_trim(outfile)-6:len_trim(outfile) )==".lammps" )
        call report(xyz)
        call output(xyz)


    !---    bye bye
        call delete(xyz)
        print *,""
        print *,"done"
        print *,""



    contains
!---^^^^^^^^


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



    end program strainXYZ