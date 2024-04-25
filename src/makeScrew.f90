
    program makeScrew
!---^^^^^^^^^^^^^^^^^^
        use Lib_CommandLineArguments
        use Lib_XYZFiles
        use iso_fortran_env
        implicit none
         
        type(CommandLineArguments)      ::      cla

    !---    input parameters        
        character(len=256)              ::      filename = ""           !    filename
        real(kind=real64)               ::      b 
        real(kind=real64),dimension(2)      ::      offset = 0.0d0
        
    !---    physically meaningful parameters deduced from input file       
        real(kind=real64),dimension(3,3)    ::      a_super , ia_super
        type(XYZFile)                       ::      xyz 
        integer                             ::      nAtoms
        
        
        
    !---    dummy parameter     
        real(kind=real64),dimension(:,:),pointer    ::      colp
        logical                 ::      ok
        integer                 ::      ii
        real(kind=real64),dimension(3)      ::      xx,yy
        real(kind=real64)       ::      zz,uu,vv
        character(len=256)      ::  dummy
        
    !---    read command line arguments
        cla = CommandLineArguments_ctor(10)  
         
        call setProgramDescription( cla, "makeScrew.exe -f filename -b burgers [-x offset]" )
        call setProgramVersion( cla, "0.1.1" )               
        
          
        call get( cla,"f",filename ,LIB_CLA_REQUIRED,"          input filename" )           
        
        call get( cla,"b",b ,LIB_CLA_REQUIRED,"          burgers vector length (A)" )           
        ii = 2 ; call get( cla,"x",offset ,ii,LIB_CLA_OPTIONAL,"          (x,y) offset of centre of displacement (A)" )     
        
            
        call report(cla)
        if (.not. allRequiredArgumentsSet(cla)) stop
        if (hasHelpArgument(cla)) stop
        call delete(cla)
            print *,""
        print *,"reading input file"
        print *,"^^^^^^^^^^^^^^^^^^"
        print *,""
        xyz = XYZFile_ctor(filename)
        call readHeader(xyz,ok)
        call input(xyz,verbose=.true.)
        call getSupercell(xyz,a_super,ok) 
        if (ok) then
            print *,"makeScrew.exe info - supercell read from file"
            print *,a_super(1,:)
            print *,a_super(2,:)
            print *,a_super(3,:)
        else
            print *,"makeScrew.exe warning - supercell not read from file"
        end if
        call report(xyz)        
        
    !---    
        call setNColumns(xyz,4)
        call getColumnsp(xyz,colp)
        nAtoms = getnAtoms(xyz)
        call inverse3Mat(a_super,ia_super)
        
        do ii = 1,nAtoms
            xx = colp(1:3,ii) 
            xx(1:2) = xx(1:2) - offset(1:2)
            
        !---    make the cut parallel to the z- axis
        !       with atoms y<L/2 , x<L/2 shift down and atoms x>L/2 shift up            
            yy = ia_super(1:3,1)*xx(1) + ia_super(1:3,2)*xx(2) + ia_super(1:3,3)*xx(3)  
            
            uu = (yy(1) + yy(2))/2              !   distance along [110] with 0 at 0,0 and 1 at 1,1
            vv = (yy(1) - yy(2))                !   distance along [1-10] with -1 at 0,1 and +1 at 1,0   
            
            
            if (vv>0) then               
                zz = sin( 3.141592654d0*uu )         
                zz = (vv-1)*zz*zz                  
            else
                zz = sin( 3.141592654d0*uu )      
                zz = (1+vv)*zz*zz
            end if
            colp(3,ii) = colp(3,ii) + b*zz/2
            colp(4,ii) = zz     
            
        end do
                
        
    !---
        call setFilename(xyz,trim(removeSuffix(filename))//".screw.xyz")
        call setNHeaderLines(xyz,0)
        write(dummy,fmt='(9f12.5)') a_super
        dummy = "Lattice="""//trim(dummy)//""" Properties=species:S:1:pos:R:3:z:R:1"
        call setColumn_description(xyz,dummy)        
        call output(xyz)            
        
    !---    
        print *,""
        print *,"done"
        print *,""
        
    contains
!---^^^^^^^^


        function removeSuffix( file_in ) result( file_out )
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      return filename with suffix removed - eg "directory/foo.bar.txt" becomes "directory/foo.bar"
    
            character(len=*),intent(in)     ::      file_in
            character(len=len(file_in))     ::      file_out
            
            integer     ::      ii
            file_out = ""
            ii = index( file_in,".",back=.true. )
            if (ii>0) then  
                file_out(1:ii-1) = file_in(1:ii-1)        
            else
                file_out = file_in
            end if
            return
        end function removeSuffix        
         
    

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
            
    end program makeScrew    