

    module Lib_XYZFiles
!---^^^^^^^^^^^^^^^^^^^^
!*
!*      Define the XYZ file format
!*
!*      Author      :   Daniel Mason
!*      Version     :   1.0
!*      Revision    :   July 2018
!*
!*
!*-----
!*
!       MIT License
!
!       Copyright (c) 2018 Daniel Robert Mason, CCFE, UKAEA, Abingdon, Oxfordshire, OX14 3DB, UK
!
!
!       Permission is hereby granted, free of charge, to any person obtaining a copy
!       of this software and associated documentation files (the "Software"), to deal
!       in the Software without restriction, including without limitation the rights
!       to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
!       copies of the Software, and to permit persons to whom the Software is
!       furnished to do so, subject to the following conditions:
!
!       The above copyright notice and this permission notice shall be included in all
!       copies or substantial portions of the Software.
!
!       THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
!       IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
!       FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
!       AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
!       LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
!       OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
!       SOFTWARE.
!
!*
!*-----
!*

        use NBAX_StringTokenizers
        use iso_fortran_env
        implicit none
        private


    !---

        public      ::      XYZFile_ctor
        public      ::      report
        public      ::      delete
        public      ::      clone       !   deep copy with allocation clone(this,that) this = that

        public      ::      readHeader
        public      ::      input
        public      ::      output,outputQhull

        public      ::      clear
        public      ::      getnAtoms
        public      ::      getnAtomNames
        public      ::      getAtomName
        public      ::      getAtomTypeFromName
        public      ::      getBoxRepeats
        public      ::      getnColumns
        public      ::      setnColumns
        public      ::      getColumn_description
        public      ::      setColumn_description
        public      ::      setFilename
        public      ::      getFilename
        public      ::      getNumberedFilename
        public      ::      setFilen


        public      ::      setnAtoms
        public      ::      setnAtomNames
        public      ::      setAtomNames
        public      ::      setAtomTypes
        public      ::      getAtomTypes
        public      ::      setAtomType
        public      ::      getAtomType
        public      ::      setColumns
        public      ::      getColumns
        public      ::      getColumnsp
        public      ::      getTypesp
        public      ::      getHeaderLines
        public      ::      setHeaderLines
        public      ::      getCol_dataType
        public      ::      setCol_dataType
        public      ::      getnHeaderLines
        public      ::      setnHeaderLines
        public      ::      displace

        public      ::      findOffset,adjustOffset

        public      ::      convertCSV
        public      ::      getMinmax
        public      ::      voxelisedCount
        public      ::      randomiseAtomTypes
        public      ::      getSupercell
        public      ::      getOrigin
        public      ::      getPBC

        public      ::      shuffle
        
        
        public      ::      idByMass
        public      ::      getMass
        public      ::      isLammpsFormat
        public      ::      setLammpsFormat
        
    !---

        integer,public,parameter        ::      XYZFILE_FILENAMELENGTH = 2048
        integer,public,parameter        ::      XYZFILE_COMMENTLINELENGTH = 256
        integer,public,parameter        ::      XYZFILE_ATOMNAMELENGTH = 6
        character(len=1),dimension(3),parameter     ::  XYZFILE_COMMENT = (/ "#","%","!" /)     !   line treated as comment if starts with this character
        
        
        integer,private,parameter       ::      XYZFILE_NATOMSIDBYMASS = 54
        character(len=2),dimension(XYZFILE_NATOMSIDBYMASS),private,parameter        ::    XYZFILE_ATOMSIDBYMASS_NAME = (/   & 
                        "H ","D ","T ","He"                                         &
                       ,'Li','Be','C ','O ','Na','Mg','Al','Si','K '                    &
                       ,'Ca','Ti','V ','Cr','Mn','Fe','Co','Ni','Cu','Zn','Ga','Ge'     &
                       ,'As','Rb','Sr','Y ','Zr','Nb','Mo','Ru','Rh','Pd','Ag','Cd'     &
                       ,'Sn','Cs','Ba','La','Ce','Gd','Hf','Ta','W ','Re','Os','Ir'     &
                       ,'Pt','Au','Tl','Pb','Th' /)  
               
        real(kind=real64),dimension(XYZFILE_NATOMSIDBYMASS),private,parameter       ::    XYZFILE_ATOMSIDBYMASS_MASS = (/   &
                        1.0d0 , 2.0d0 , 3.0d0 , 4.003d0                                     &
                        ,7.0d0,9.012183d0,12.011d0,15.999d0,22.9897693d0,24.305d0,26.981538d0,28.085d0,39.0983d0                                 &
                        ,40.08d0,47.867d0,50.9415d0,51.996d0,54.93804d0,55.84d0,58.93319d0,58.693d0,63.55d0,65.4d0,69.723d0,72.63d0              &
                        ,74.92159d0,85.468d0,87.62d0,88.90584d0,91.22d0,92.90637d0,95.95d0,101.1d0,102.9055d0,106.42d0,107.868d0,112.41d0        &
                        ,118.71d0,132.9054520d0,137.33d0,138.9055d0,140.116d0,157.2d0,178.49d0,180.9479d0,183.84d0,186.207d0,190.2d0,192.22d0    &
                        ,195.08d0,196.96657d0,204.383d0,207d0,232.038d0  /)
        
        logical,public      ::      LIBXYZ_DBG = .false.                                                                            
        integer,public,parameter      ::      LIBXYZ_COL_STRING = 0
        integer,public,parameter      ::      LIBXYZ_COL_FLOAT = 1
        integer,public,parameter      ::      LIBXYZ_COL_UINT = 2 
        integer,public,parameter      ::      LIBXYZ_COL_INT = 3
        
        character(len=5),dimension(0:3),private,parameter  ::  LIBXYZ_COL_NAME = (/ "str  ","float","uint ","int  " /)                                                                                     
        
        
    !---

        type,public     ::         XYZFile
            private
            integer                 ::      filen                               !   filenumber -1 for none, or 00000,00001,...
            integer                 ::      nHeaderLines
            integer                 ::      nAtoms
            integer                 ::      nColumns
            integer                 ::      nAtomNames
            character(len=XYZFILE_FILENAMELENGTH)                               ::      filename
            character(len=XYZFILE_COMMENTLINELENGTH),dimension(:),pointer       ::      header
            character(len=XYZFILE_COMMENTLINELENGTH)                            ::      column_description
            character(len=XYZFILE_ATOMNAMELENGTH),dimension(:),pointer          ::      atomName        !   (1:nAtomNames)
            integer,dimension(:),pointer                                        ::      atom            !   (1:nAtoms)
            real(kind=real64),dimension(:,:),pointer                            ::      dat             !   (1:nColumns,1:nAtoms)
            integer,dimension(:),pointer                                        ::      col_dataType
            logical                 ::      lammpsFormat
        end type

    !---




        interface   XYZFile_ctor
            module procedure        XYZFile_null
            module procedure        XYZFile_ctor0
            module procedure        XYZFile_ctor1
        end interface


        interface   report
            module procedure        report0
        end interface

        interface   delete
            module procedure        delete0
        end interface

        interface   input
            module procedure        input0
        end interface

        interface   output
            module procedure        output0
            module procedure        output1
        !    module procedure        output2
        end interface

        interface   clear
            module procedure        clear0
            module procedure        clear1
        end interface

        interface   clone
            module procedure        clone0
        end interface

        interface   getFilename
            module procedure        getFilename0
        end interface




        interface   getnAtoms
            module procedure        getnAtoms0
            module procedure        getnAtoms1
        end interface

        interface   getnAtomNames
            module procedure        getnAtomNames0
        end interface

        interface   getAtomName
            module procedure        getAtomName0
        end interface

        interface   getnColumns
            module procedure        getnColumns0
        end interface

        interface   setnColumns
            module procedure        setnColumns0
        end interface
        
        interface   setColumn_Description
            module procedure        setColumn_Description0
            module procedure        setColumn_Description1
            module procedure        setColumn_Description2
            module procedure        setColumn_Description3
        end interface
        
        

        interface   setnHeaderLines
            module procedure        setnHeaderLines0
        end interface

        interface   setnAtomNames
            module procedure        setnAtomNames0
        end interface

        interface   getnHeaderLines
            module procedure        getnHeaderLines0
        end interface

        interface   setAtomTypes
            module procedure        setAtomTypes0
            module procedure        setAtomTypes1
        end interface

        interface   getAtomTypes
            module procedure        getAtomTypes0
        end interface

        interface   setAtomType
            module procedure        setAtomType0
        end interface

        interface   getAtomType
            module procedure        getAtomType0
        end interface

        interface   setColumns
            module procedure        setColumns0
            module procedure        setColumns1
            module procedure        setColumns2
        end interface

        interface   getColumns
            module procedure        getColumns0
            module procedure        getColumns1
        end interface

        interface   getHeaderLines
            module procedure        getHeaderLines0
            module procedure        getHeaderLines1
        end interface

        interface   setHeaderLines
            module procedure        setHeaderLines0
            module procedure        setHeaderLines0a
            module procedure        setHeaderLines1
            module procedure        setHeaderLines2
        end interface


        interface   getCol_dataType
            module procedure        getCol_dataType0
            module procedure        getCol_dataType1
        end interface

        interface   setCol_dataType
            module procedure        setCol_dataType0
            module procedure        setCol_dataType1
            module procedure        setCol_dataType2
        end interface


        interface   setAtomNames
            module procedure        setAtomNames0
            module procedure        setAtomNames0a
        end interface

        interface   displace
            module procedure        displace0
        end interface

        interface   shuffle
            module procedure        shuffle0
        end interface


        interface   adjustOffset
            module procedure        adjustOffset0
            module procedure        adjustOffset1
        end interface
        
        interface   voxelisedCount
            module procedure        voxelisedCount0
            module procedure        voxelisedCount1
            module procedure        voxelisedCount2
            module procedure        voxelisedCount1a
            module procedure        voxelisedCount2a
        end interface


        interface   getSupercell
            module procedure        getSupercell0
            module procedure        getSupercell1
        end interface
        
        interface   getOrigin
            module procedure        getOrigin0
        end interface
        
        interface   getPBC
            module procedure        getPBC0
        end interface
        
        interface   getMass
            module procedure        getMass0
        end interface
        
        
    !---

    contains
!---^^^^^^^^

        function XYZFile_null() result(this)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
            type(XYZFile)       ::      this
            this%nAtoms = 0
            this%nAtomNames = 0
            this%nColumns = 0
            this%nHeaderLines = 0
            this%filename = ""
            this%column_description = ""
            this%filen = -1
            nullify(this%header)
            nullify(this%atomName)
            nullify(this%atom)
            nullify(this%dat)
            nullify(this%col_dataType)
            this%lammpsFormat = .false.
            return
        end function XYZFile_null

        function XYZFile_ctor0(nAtoms,nAtomNames,nColumns,nHeaderLines) result(this)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
            integer,intent(in)      ::      nAtoms,nAtomNames,nColumns,nHeaderLines
            type(XYZFile)       ::      this
            this = XYZFile_null()
            call ensureAdequateStorage(this,nAtoms,nAtomNames,nColumns,nHeaderLines)
            this%filename = ""
            this%column_description = ""
            this%lammpsFormat = .false.
            return
        end function XYZFile_ctor0

        function XYZFile_ctor1(filename) result(this)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
            character(len=*),intent(in)     ::      filename
            type(XYZFile)       ::      this
            this = XYZFile_null()
            this%filename = trim(filename)
            this%lammpsFormat = .false.
            return
        end function XYZFile_ctor1

    !---


        subroutine report0(this,u,o)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^
            type(XYZFile),intent(in)     ::      this
            integer,intent(in),optional             ::      u,o
            integer     ::      uu,oo
            integer     ::      ii,jj
            real(kind=real64),dimension(3,this%nColumns,0:this%nAtomNames)  ::      analysis_dat
            integer,dimension(0:this%nAtomNames)                            ::      atom_count
            
            uu = 6 ; if (present(u)) uu = u
            oo = 0 ; if (present(o)) oo = o
            write (unit=uu,fmt='(a)') repeat(" ",oo)//"XYZFile["//trim(this%filename)//"]"
            write (unit=uu,fmt='(a,i12)') repeat(" ",oo+2)//"header lines: ",this%nHeaderLines
            write (unit=uu,fmt='(a,i12)') repeat(" ",oo+2)//"data columns: ",this%nColumns
            write (unit=uu,fmt='(a,i12)') repeat(" ",oo+2)//"atoms       : ",this%nAtoms
            if (this%nAtomNames>0) then
                write (unit=uu,fmt='(a)',advance="no") repeat(" ",oo+2)//"atom names  : "
                do ii = 1,this%nAtomNames-1
                    write (unit=uu,fmt='(a)',advance="no") this%atomName(ii)//","
                end do
                write (unit=uu,fmt='(a)',advance="yes") this%atomName(this%nAtomNames)
                if (this%nAtoms>0) then
                    atom_count = 0
                    analysis_dat(1,:,:) = 0.0d0         !   will be average 
                    analysis_dat(2,:,:) = huge(1.0)     !   minimum 
                    analysis_dat(3,:,:) = -huge(1.0)    !   maximum
                    do ii = 1,this%nAtoms
                        jj = this%atom(ii)
                        atom_count(jj) = atom_count(jj) + 1
                        analysis_dat(1,1:this%nColumns,jj) = analysis_dat(1,1:this%nColumns,jj) + this%dat( 1:this%nColumns,ii )
                        analysis_dat(2,1:this%nColumns,jj) = min(analysis_dat(2,1:this%nColumns,jj),this%dat( 1:this%nColumns,ii ) )
                        analysis_dat(3,1:this%nColumns,jj) = max(analysis_dat(3,1:this%nColumns,jj),this%dat( 1:this%nColumns,ii ) )
                    end do
                    
                    
                    write (unit=uu,fmt='(a)') repeat(" ",oo+2)//"atom counts"
                    do jj = 1,this%nAtomNames
                        write (unit=uu,fmt='(a,i8,a,f7.3,a)' ) repeat(" ",oo+4)//this%atomName(jj),atom_count(jj)," (",atom_count(jj)*100.0d0/this%nAtoms," %)"
                    end do               
                    
                             
                    ! !print *,"this%atom ",minval(this%atom(1:this%nAtoms)),maxval(this%atom(1:this%nAtoms))
                    ! do ii = 1,this%nColumns
                    !     print *,"col ",ii,minval(this%dat(ii,1:this%nAtoms)),maxval(this%dat(ii,1:this%nAtoms)),minloc(this%dat(ii,1:this%nAtoms),dim=1),maxloc(this%dat(ii,1:this%nAtoms),dim=1)
                    ! end do
                    
                    
                    write (unit=uu,fmt='(a,a8,a6,a8,3a16)') repeat(" ",oo+1),"column","type","(atom)","average","minimum","maximum"                 
                    do jj = 1,this%nAtomNames
                        if (atom_count(jj)==0) cycle
                        analysis_dat(1,1:this%nColumns,jj) = analysis_dat(1,1:this%nColumns,jj) / atom_count(jj)
                        write (unit=uu,fmt='(a,a8,a6,a8,3g16.6)') repeat(" ",oo+2),"x ",LIBXYZ_COL_NAME(this%col_dataType(1)),this%atomName(jj),analysis_dat(1:3,1,jj)
                        write (unit=uu,fmt='(a,a8,a6,a8,3g16.6)') repeat(" ",oo+2),"y ",LIBXYZ_COL_NAME(this%col_dataType(2)),this%atomName(jj),analysis_dat(1:3,2,jj)
                        write (unit=uu,fmt='(a,a8,a6,a8,3g16.6)') repeat(" ",oo+2),"z ",LIBXYZ_COL_NAME(this%col_dataType(3)),this%atomName(jj),analysis_dat(1:3,3,jj)
                        do ii = 4,this%nColumns
                            if (this%col_dataType(ii) == 1) then
                                write (unit=uu,fmt='(a,i7,a7,a8,3g16.6)') repeat(" ",oo+2),ii," "//LIBXYZ_COL_NAME(this%col_dataType(ii)),this%atomName(jj),analysis_dat(1:3,ii,jj)
                            else
                                write (unit=uu,fmt='(a,i7,a7,a8,g16.6,2i16)') repeat(" ",oo+2),ii," "//LIBXYZ_COL_NAME(this%col_dataType(ii)),this%atomName(jj),analysis_dat(1,ii,jj),nint( analysis_dat(2:3,ii,jj) )
                            end if                            
                        end do
                    end do

                end if
            end if
            
!             if (this%nAtoms>0) then
!             
!             
!                 do ii = 1,this%nColumns
!                     avg(ii) = sum(this%dat(ii,1:this%nAtoms))
!                 end do
!                 avg(1:this%nColumns) = avg(1:this%nColumns)/this%nAtoms
!                     
!             
!                 write (unit=uu,fmt='(a)') repeat(" ",oo+2)//"  column     minimum      (atom)             maximum     (atom)             average"                 
!                 n1 = int(minloc(this%dat(1,1:this%nAtoms),dim=1)) ; n2 = int(maxloc(this%dat(1,1:this%nAtoms),dim=1))            
!                 write (unit=uu,fmt='(a,a8,3(g16.6,a8,i8))') repeat(" ",oo+2),"x",this%dat(1,n1)," "//this%atomName(this%atom(n1)),n1 ,this%dat(1,n2)," "//this%atomName(this%atom(n2)),n2,avg(1)
!                 n1 = int(minloc(this%dat(2,1:this%nAtoms),dim=1)) ; n2 = int(maxloc(this%dat(2,1:this%nAtoms),dim=1))                
!                 write (unit=uu,fmt='(a,a8,3(g16.6,a8,i8))') repeat(" ",oo+2),"y",this%dat(2,n1)," "//this%atomName(this%atom(n1)),n1 ,this%dat(2,n2)," "//this%atomName(this%atom(n2)),n2,avg(2)
!                 n1 = int(minloc(this%dat(3,1:this%nAtoms),dim=1)) ; n2 = int(maxloc(this%dat(3,1:this%nAtoms),dim=1))                
!                 write (unit=uu,fmt='(a,a8,3(g16.6,a8,i8))') repeat(" ",oo+2),"z",this%dat(3,n1)," "//this%atomName(this%atom(n1)),n1 ,this%dat(3,n2)," "//this%atomName(this%atom(n2)),n2,avg(3)
!                 do ii = 4,this%nColumns
!                     n1 = int(minloc(this%dat(ii,1:this%nAtoms),dim=1)) ; n2 = int(maxloc(this%dat(ii,1:this%nAtoms),dim=1))                
!                     write (unit=uu,fmt='(a,i8,3(g16.6,a8,i8))') repeat(" ",oo+2),ii,this%dat(ii,n1)," "//this%atomName(this%atom(n1)),n1 ,this%dat(ii,n2)," "//this%atomName(this%atom(n2)),n2,avg(ii)
!                 end do                    
!             end if
            return
        end subroutine report0


        subroutine delete0(this)
    !---^^^^^^^^^^^^^^^^^^^^^^^^
            type(XYZFile),intent(inout)  ::      this
            
            if (this%nHeaderLines >0 ) deallocate(this%header)
            if (this%nAtomNames >0 )   deallocate(this%atomName)
            if (this%nAtoms >0 )       deallocate(this%atom)
            if (this%nColumns >0 )     then
                deallocate(this%dat)
                deallocate(this%col_DataType)
            end if
            this = XYZFile_null()
            return
        end subroutine delete0

    !---

        subroutine ensureAdequateStorage(this,nAtoms,nAtomNames,nColumns,nHeaderLines)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
            type(XYZFile),intent(inout)     ::      this
            integer,intent(in)              ::      nAtoms,nAtomNames,nColumns,nHeaderLines
            call setNColumns( this,nColumns,nAtoms )
            call setNHeaderLines( this,nHeaderLines )
            call setNAtomNames( this,nAtomNames )
            return
        end subroutine ensureAdequateStorage




    !---

        subroutine readHeader(this,ok)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      read header. if necessary, reallocate memory
            type(XYZFile),intent(inout)     ::      this
            logical,intent(out)             ::      ok

            integer                         ::      uu
            type(StringTokenizer)           ::      st
            integer                 ::      ii,ioerr
            integer                 ::      nAtoms,nColumns,nHeaderLines
            character(len=XYZFILE_COMMENTLINELENGTH)    ::  dummy
            character(len=4096)     ::      longDummy
            character(len=1)        ::      letter

            real(kind=real64),dimension(1000)   ::  dat


           !ii = index(trim(this%filename),".lammps",back=.true.)
           !if ( (ii>0).and.(ii >= len_trim(this%filename)-7) ) then                
           !    call readHeader_lammps(this,ok)
           !    return
           !end if            
            
            
            inquire(file=trim(this%filename),exist=ok)
            if (.not. ok) return
            
            uu = findGoodUnit()
            open(unit=uu,file=trim(this%filename),action="read")
                read(unit=uu,fmt='(a)',iostat=ioerr) dummy
                if (ioerr/=0) then
                    print *,"Lib_XYZFiles::readHeader error - could not read first line of file """//trim(this%filename)//""""
                    stop
                end if
                dummy = adjustl(dummy)
                if (dummy(1:6) == "LAMMPS") then
                    close(unit=uu)
                    call readHeader_lammps(this,ok)      
                    this%lammpsFormat = .true.                      
                    return
                end if
            close(unit=uu)
             
            
            

            nHeaderLines = 0
            open(unit=uu,file=trim(this%filename),action="read")
                do ii = 1,1000      !   don't loop forever...
                    read(unit=uu,fmt='(a)',iostat=ioerr) dummy

                    if (ioerr /= 0) then
                        print *,"Lib_XYZFiles::readHeader error - could not read header line ",ii
                        stop
                    end if

                    dummy = adjustl(dummy) ; letter = dummy(1:1)
                    if (any( letter == XYZFILE_COMMENT(:) )) then
                        !   a comment line
                        nHeaderLines = nHeaderLines + 1
                    else
                        !   not a comment line. Must be the atom count line
                        exit
                    end if
                end do

                st = StringTokenizer_ctor( dummy, " " )
                call nextToken(st,dummy)
                call parse( dummy,nAtoms,ok )
                if (.not. ok) then
                    print *,"Lib_XYZFiles::readHeader error - could not parse number of atoms from line "//trim(dummy)
                    stop
                end if

                read(unit=uu,fmt='(a)',iostat=ioerr) this%column_description
                if (ioerr /= 0) then
                    print *,"Lib_XYZFiles::readHeader error - could not read column description line"
                    stop
                end if


                read(unit=uu,fmt='(a)',iostat=ioerr) longDummy
                if (ioerr /= 0) then
                    print *,"Lib_XYZFiles::readHeader error - could not read first data line"
                    stop
                end if

                st = StringTokenizer_ctor( longDummy, " " )
                call nextToken(st,dummy)
                call parse(getRemainingTokens(st),dat,nColumns)
                if (nColumns < 3) then
                    print *,"Lib_XYZFiles::readHeader error - expected at least 3 columns for position x,y,z. Found ",nColumns
                    stop
                end if

            close(unit=uu)


            call ensureAdequateStorage(this,nAtoms,0,nColumns,nHeaderLines)

            this%nAtoms = nAtoms
            this%nAtomNames = 0
            this%nColumns = nColumns
            this%nHeaderLines = nHeaderLines
            
            
            open(unit=uu,file=trim(this%filename),action="read")
                do ii = 1,nHeaderLines
                    read(unit=uu,fmt='(a)',iostat=ioerr) this%header(ii)
                end do
            close(unit=uu)
            
            
            

            ok = .true.

            return
        end subroutine readHeader


        subroutine readHeader_lammps(this,ok)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      read header, assuming .lammps format file. if necessary, reallocate memory
            type(XYZFile),intent(inout)     ::      this
            logical,intent(out)             ::      ok

            integer                         ::      uu
            integer                 ::      ii,ioerr,jj,kk
            integer                 ::      nAtoms,nColumns,nHeaderLines,nAtomNames
            character(len=XYZFILE_COMMENTLINELENGTH)    ::  dummy
            character(len=XYZFILE_ATOMNAMELENGTH),dimension(:),allocatable       ::      newatomName
            real(kind=real64)       ::      mm


!             ii = index(trim(this%filename),".lammps",back=.true.)
!             if ( ii==0 ) then
!                 print *,"Lib_XYZFiles::readHeader_lammps error - file does not have expected .lammps suffix """//trim(this%filename)//"""" 
!                 stop
!             end if
            uu = findGoodUnit()
            inquire(file=trim(this%filename),exist=ok)
            if (.not. ok) return

            nHeaderLines = 0
            open(unit=uu,file=trim(this%filename),action="read")
                do ii = 1,1000      !   don't loop forever...
                    read(unit=uu,fmt='(a)',iostat=ioerr) dummy

                    if (ioerr /= 0) then
                        print *,"Lib_XYZFiles::readHeader_lammps error - could not read header line ",ii
                        stop
                    end if

                    dummy = adjustl(dummy) 
                    
                    if (index(trim(dummy),"atoms")/=0) then
                        !   this is the expected position of the atom count line.
                        call parse( dummy,nAtoms,ok )
                        if (.not. ok) then
                            print *,"Lib_XYZFiles::readHeader_lammps error - expected line to have number of atoms, actually read """//trim(dummy)//""""
                            stop
                        else
                            print *,"Lib_XYZFiles::readHeader_lammps info - read atom count ",nAtoms
                        end if
                    end if
                    
                    if (index(dummy,"atom types")/=0) then
                        !   this line carries number of atom types...
                        call parse( dummy,nAtomNames,ok )
                        if (.not. ok) then
                            print *,"Lib_XYZFiles::readHeader_lammps error - expected line to have number of atom types, actually read """//trim(dummy)//""""
                            stop
                        else
                            print *,"Lib_XYZFiles::readHeader_lammps info - read number of atom types ",nAtomNames
                        end if
                        allocate(newatomName(nAtomNames))
                    end if
                    
                    if (trim(dummy(1:6))=="Masses") then
                        !   next line should be blank, then mass of each atom type
                        read(unit=uu,fmt='(a)',iostat=ioerr) dummy
                        do jj = 1,nAtomNames
                            read(unit=uu,fmt='(a)',iostat=ioerr) dummy
                            read(dummy,fmt=*,iostat=ioerr) kk,mm
                            if (ioerr/=0) then
                                print *,"Lib_XYZFiles::readHeader_lammps error - expected to read mass of atom type ",jj," could not parse """//trim(dummy)//""""
                                stop
                            end if
                            newatomName(kk) = idByMass(mm)
                            if (trim(newatomName(kk))=="unset") then
                                write(dummy,fmt='(i2)') ii ; dummy = adjustl(dummy) ; dummy = "atom"//repeat("0",2-len_trim(dummy))//trim(dummy)
                                newatomName(kk) = trim(dummy)
                            else
                                print *,"Lib_XYZFiles::readHeader_lammps info - identified atom """//trim(newatomName(kk))//""" by mass ",mm
                            end if                            
                        end do    
                        nHeaderLines = nHeaderLines + 1 + nAtomNames
                    end if                 
                    
                    
                    if (trim(dummy(1:5))=="Atoms") then
                        !   next line should be blank, then the data...
                        read(unit=uu,fmt='(a)',iostat=ioerr) dummy
                        nHeaderLines = nHeaderLines + 2
                        exit
                    end if
                    
                    nHeaderLines = nHeaderLines + 1
                    
                end do


                this%column_description = "atom position_x position_y position_z"
                nColumns = 3    

            close(unit=uu)


            call ensureAdequateStorage(this,nAtoms,nAtomNames,nColumns,nHeaderLines)

            
            open(unit=uu,file=trim(this%filename),action="read")
                do ii = 1,nHeaderLines
                    read(unit=uu,fmt='(a)',iostat=ioerr) dummy 
                    this%header(ii) = "# "//trim(dummy)
                end do
            close(unit=uu)
            
            
!             
!             
!                 do ii = 1,nAtomNames
!                     write(dummy,fmt='(i2)') ii ; dummy = adjustl(dummy) ; dummy = "atom"//repeat("0",2-len_trim(dummy))//trim(dummy)
!                     newatomName(ii) = trim(dummy)
!                 end do
                call setAtomNames(this,newAtomName)
                
            deallocate(newatomName)
            

            ok = .true.

            return
        end subroutine readHeader_lammps
        
        
    !---

        subroutine input0(this,verbose)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
             type(XYZFile),intent(inout)     ::      this
            logical,intent(in),optional     ::      verbose
            integer                                     ::      uu
            integer                                     ::      ii,jj,nn,ioerr,nAtomNames
            character(len=XYZFILE_ATOMNAMELENGTH)       ::      atom
            logical                                     ::      ok,allokdat,allokatm,findAtomNames
            character(len=XYZFILE_ATOMNAMELENGTH),dimension(1000)       ::      atom_names

            logical         ::      op
            character(len=XYZFILE_FILENAMELENGTH)       ::      f
            
            
!            character(len=256)    ::      dummy
!            integer     ::      kk
            
                     
            if (this%lammpsFormat) then
                call input_lammps(this,verbose)
                return
            end if            
                        
            
            
            f = getNumberedFilename(this)


            op = .false. ; if (present(verbose)) op = verbose

            uu = findGoodUnit()
            findAtomNames = (this%nAtomNames == 0)
            nAtomNames = 0
            allokdat = .true.
            allokatm = .true.
            
            open(unit=uu,file=trim(f),action="read")


 
                do ii = 1,this%nHeaderLines
                    read(unit=uu,fmt='(a)',iostat=ioerr) this%header(ii)
                end do
                read(unit=uu,fmt=*,iostat=ioerr) nn

                if (nn > this%nAtoms) then
                    print *,"Lib_XYZFiles::input0 error - read number of atoms ",nn," more than expected ",this%nAtoms
                    stop
                end if
                read(unit=uu,fmt='(a)',iostat=ioerr) this%column_description
                if (op) write(unit=*,fmt='(a)') trim(this%column_description)

                
                
                do ii = 1,nn                 
                
                    read(unit=uu,fmt=*,iostat=ioerr) atom,this%dat(1:this%nColumns,ii)         !   spectacularly slow
                    
                    
                    if ( (ii==1).and.(ioerr/=0) ) then
                        print *,"Lib_XYZFiles::input0 error - attempting to parse first line into ""species"" then ",this%nColumns," columns failed"
                        stop
                    end if
                     
                    !read(unit=uu,fmt='(a)',iostat=ioerr) dummy       !   spectacularly slow
                    !read(unit=uu,fmt='(a)') dummy       !   spectacularly slow
                    !kk = index(dummy," ")
                    !atom = dummy(1:kk)
                    !call parse( dummy(kk:),this%dat(1:this%nColumns,ii),kk )
                    
                    
                    
                    allokdat = allokdat .and. (ioerr==0)
                    
                    if (findAtomNames) then
                        ok = .false.
                        do jj = 1,nAtomNames
                            if (trim(atom) == trim(atom_names(jj))) then
                                this%atom(ii) = jj
                                ok = .true.
                                exit
                            end if
                        end do
                        if (.not. ok) then
                            nAtomNames = nAtomNames + 1
                            if (nAtomNames > size(atom_names)) then
                                print *,"Lib_XYZFiles::input0 error - attempting to add atom type number ",nAtomNames
                                print *,"                     note: expect .xyz file format ""species"" position_x position_y position_z [extra columns]"
                                
                                stop
                            end if
                            atom_names(nAtomNames) = trim(atom)
                            this%atom(ii) = nAtomNames
                        end if
                    else
                        this%atom(ii) = getAtomTypeFromName( this,atom )
                        allokatm = allokatm .and. (this%atom(ii)>0)
                    end if

                    if (op) then
                        if (ii<4) then
                            write(unit=*,fmt='(a8,1000f16.6)') atom,this%dat(1:this%nColumns,ii)
                        else if ( (ii==4).and.(nn>5) ) then
                            write(unit=*,fmt='(a)') "... reading ..."
                        else if ( (ii>nn-3).and.(nn>5) ) then
                            write(unit=*,fmt='(a8,1000f16.6)') atom,this%dat(1:this%nColumns,ii)
                        end if
                    end if


                end do
            close(unit=uu)

            if (.not. allokatm) then
                print *,"Lib_XYZFiles::input0 error - could not read atom names correctly"
                stop
            end if


            if (.not. allokdat) then
                print *,"Lib_XYZFiles::input0 error - could not read data correctly"
                stop
            end if

            if (findAtomNames) then
                if (this%nAtomNames>0) deallocate(this%atomName)
                allocate(this%atomName(nAtomNames))
                this%nAtomNames = nAtomNames
                this%atomName(1:this%nAtomNames) = atom_names(1:this%nAtomNames)
            end if


            
            
            
            return
        end subroutine input0

    !---        

        subroutine getOrigin0(this,origin,ok)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      looks for attribute 'Origin="ox oy oz"' in the column description
            type(XYZFile),intent(in)                    ::      this
            real(kind=real64),dimension(3),intent(out)  ::      origin
            logical,intent(out)                         ::      ok

            character(len=XYZFILE_COMMENTLINELENGTH)    ::  dummy
            integer                         ::      ii

            ok = .false.
            origin = 0.0d0      !   assume no information in file

            ii = index(trim(convertUpperCase(this%column_description)),"ORIGIN=")
            if (ii/=0) then
                dummy = this%column_description(ii+8:)
                call parse(dummy,origin,ii)
                ok = (ii==3)
            end if

            if (.not. ok) then
                ii = index(trim(convertUpperCase(this%column_description)),"BOXSIZE")       !   check for parcas format
                if (ii/=0) then
                    dummy = this%column_description(ii+8:)
                    call parse( dummy,origin,ii )
                    ok = (ii==3)
                    if (ok) then
                        origin = -origin/2                      !   parcas has atoms from -A/2 to A/2
                    else
                        origin = 0
                    end if
                end if            
            end if
            return
        end subroutine getOrigin0


        subroutine getPBC0(this,pbc,ok)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ 
    !*      looks for attribute 'pbc="px py pz"' in the column description
    !*      with px = "true","T",".true." etc
            type(XYZFile),intent(in)            ::      this
            logical,dimension(3),intent(out)    ::      pbc
            logical,intent(out)                 ::      ok

            character(len=XYZFILE_COMMENTLINELENGTH)    ::  dummy
            integer                         ::      ii
            ok = .false.
            pbc = .true.            !   default if no information in file

            ii = index(trim(convertUpperCase(this%column_description)),"PBC=")
            if (ii/=0) then
                dummy = this%column_description(ii+5:)
                call parse(dummy,pbc,ii)
                ok = (ii==3)
            end if
            return
        end subroutine getPBC0


        subroutine getSupercell1(this,a,ok)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
            type(XYZFile),intent(inout)     ::      this
            real(kind=real64),dimension(3,3),intent(out)    ::      a            
            logical,intent(out)             ::      ok        
            integer         ::      nx,ny,nz
            call getSupercell0(this,a,nx,ny,nz,ok)
            return
        end subroutine getSupercell1
        
        subroutine getSupercell0(this,a,nx,ny,nz,ok)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      the supercell size may be stored with the .xyz file.
    !*      look for it, and return it if possible
            type(XYZFile),intent(inout)     ::      this
            real(kind=real64),dimension(3,3),intent(out)    ::      a
            integer,intent(out)             ::      nx,ny,nz
            logical,intent(out)             ::      ok
            
            integer                                     ::      uu
            integer                                     ::      ii,jj,kk,ioerr
            
            real(kind=real64),dimension(3)              ::      xmin,xmax , offdiag
            real(kind=real64),dimension(9)              ::      aaaa
 

            character(len=XYZFILE_COMMENTLINELENGTH)    ::  dummy,dummy2

            ok = .false.
            a = reshape( (/1,0,0,0,1,0,0,0,1/),(/3,3/) )
            nx = 1 ; ny = 1 ; nz = 1
            if (this%lammpsFormat) then
                print *,"Lib_XYZFiles::getSupercell0 info - file interpretted as .lammps format" 
                uu = findGoodUnit()
                xmax = 0.0d0 ; xmin = 0.0d0 ; offdiag = 0.0d0  
                open(unit=uu,file=trim(this%filename),action="read")
    
                    do ii = 1,this%nHeaderLines
                        read(unit=uu,fmt='(a)',iostat=ioerr) dummy
    
                        jj = index(dummy,"xlo") 
                        if (jj/=0) read(dummy,fmt=*) xmin(1),xmax(1)
                        jj = index(dummy,"ylo") 
                        if (jj/=0) read(dummy,fmt=*) xmin(2),xmax(2)
                        jj = index(dummy,"zlo") 
                        if (jj/=0) read(dummy,fmt=*) xmin(3),xmax(3)
                                                        
                        
                        jj = index(dummy,"xy xz yz") 
                        if (jj/=0) read(dummy,fmt=*) offdiag(1:3)                                       
                              
                    end do
                    if (norm2(xmax)>norm2(xmin)) then
                        xmax = xmax - xmin
                        a(1,1) = xmax(1) ; a(2,2) = xmax(2) ; a(3,3) = xmax(3)
                        a(1,2) = offdiag(1) ; a(1,3) = offdiag(2) ; a(2,3) = offdiag(3) 
                        ok = .true.
                        
                    end if
                close(unit=uu)            
                return    
            end if    
            
            ii = index(trim(convertUpperCase(this%column_description)),"LATTICE")
            if (ii/=0) then
                print *,"Lib_XYZFiles::getSupercell0 info - file interpretted as extended-xyz format" 
               ! uu = findGoodUnit()
               ! open(unit=uu,file=trim(this%filename),action="read")
                    dummy2 = this%column_description(ii+9:)
                    call parse( dummy2,aaaa,jj )
                    if (jj==3) then
                        a(1,1) = aaaa(1) ; a(2,2) = aaaa(2) ; a(3,3) = aaaa(3)                        
                        ok = .true.
                    else if (jj==9) then
                        a(1:3,1:3) = reshape(aaaa,(/3,3/))
                        ok = .true.
                    else
                        print *,"Lib_XYZFiles::getSupercell0 ERROR - can't parse """//trim(dummy2)//""""
                    end if
              !  close(unit=uu)            
                return    
            end if    
                
            ii = index(trim(convertUpperCase(this%column_description)),"BOXSIZE")
            if (ii/=0) then
                print *,"Lib_XYZFiles::getSupercell0 info - file interpretted as PARCAS-xyz format" 
               ! uu = findGoodUnit()
               ! open(unit=uu,file=trim(this%filename),action="read")
                    dummy2 = this%column_description(ii+8:)
                    call parse( dummy2,aaaa,jj )
                    if (jj>=3) then
                        a(1,1) = aaaa(1) ; a(2,2) = aaaa(2) ; a(3,3) = aaaa(3)                        
                        ok = .true.
                    else
                        print *,"Lib_XYZFiles::getSupercell0 ERROR - can't parse """//trim(dummy2)//""""
                    end if
                !close(unit=uu)            
                return    
            end if    
            
            
            print *,"Lib_XYZFiles::getSupercell0 WARNING - unknown file format - hunting for text ""supercell"" " 
            uu = findGoodUnit()
            open(unit=uu,file=trim(this%filename),action="read")

                do ii = 1,this%nHeaderLines
                    read(unit=uu,fmt='(a)',iostat=ioerr) dummy

                    jj = max( index(dummy,"supercell"), index(dummy,"super cell") )
                    if (jj/=0) then
                        dummy = dummy(jj+9:)
                        call parse( dummy,aaaa,jj )
                        if (jj==3) then
                            xmax(1:3) = aaaa(1:3)
                            if (LIBXYZ_DBG) write(*,fmt='(a,3i6)') " Lib_XYZFiles::getSupercell0 info - parsing """//trim(dummy)//""" as cell repeats ",nint(xmax)
                            nx = nint(xmax(1)) ; ny = nint(xmax(2)) ; nz = nint(xmax(3))
                        else
                            xmax = 1
                            print *,"Lib_XYZFiles::getSupercell0 warning  - can't parse """//trim(dummy)//""" as cell repeats. Assuming 1x1x1 unit cells."

                        end if
                        
                        do kk = 1,3
                            read(unit=uu,fmt='(a)',iostat=ioerr) dummy
                            call parse( dummy(3:),aaaa,jj )
                            if (jj==3) then
                                a(kk,1:3) = aaaa(1:3)
                                if (LIBXYZ_DBG) write(*,fmt='(a,9f12.5)') " Lib_XYZFiles::getSupercell0 info - parsing """//trim(dummy(3:))//""" as lattice vectors ",aaaa(1:3)
                                if (kk==3) then 
                                    ok = .true.
                                    a(:,1) = a(:,1)*nint(xmax(1))
                                    a(:,2) = a(:,2)*nint(xmax(2))
                                    a(:,3) = a(:,3)*nint(xmax(3))                                    
                                    close(unit=uu)     
                                    return
                                end if
                            else if (jj==9) then
                                if (LIBXYZ_DBG) write(*,fmt='(a,9f12.5)') " Lib_XYZFiles::getSupercell0 info - parsing """//trim(dummy(3:))//""" as lattice vectors ",aaaa(1:9)
                                a(:,:) = reshape( aaaa(1:9),(/3,3/) )
                                a(:,1) = a(:,1)*nint(xmax(1))
                                a(:,2) = a(:,2)*nint(xmax(2))
                                a(:,3) = a(:,3)*nint(xmax(3))
                                ok = .true.
                                close(unit=uu)     
                                return
                            else 
                                print *,"Lib_XYZFiles::getSupercell0 ERROR - can't parse """//trim(dummy)//""" as cell vectors."
                                ok = .false.            
                                close(unit=uu)                             
                                return
                            end if
                        end do
                            
                        exit
                    end if                        
                                            
                    
                end do
                print *,"Lib_XYZFiles::getSupercell0 WARNING - unknown file format could not read supercell"                 
            close(unit=uu)            
            return    
        end subroutine getSupercell0
                
          

        subroutine input_lammps(this,verbose)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
            type(XYZFile),intent(inout)     ::      this
            logical,intent(in),optional     ::      verbose
            integer                                     ::      uu
            integer                                     ::      ii,jj,nn,ioerr
            
            logical                                     ::      allokdat
            real(kind=real64),dimension(3)              ::      xmin,xmax,xx

            logical         ::      op
            character(len=XYZFILE_COMMENTLINELENGTH)    ::  dummy
                     
            
            


            op = .false. ; if (present(verbose)) op = verbose

            uu = findGoodUnit()
            allokdat = .true.
            
            xmax = 0.0d0 ; xmin = 0.0d0
            
            inquire(file=trim(this%filename),exist=allokdat)
            if (.not. allokdat) then
                print *,"Lib_XYZFiles::input_lammps error - file """//trim(this%filename)//""" not found"
                return
            end if
            open(unit=uu,file=trim(this%filename),action="read")

                do ii = 1,this%nHeaderLines
                    read(unit=uu,fmt='(a)',iostat=ioerr) dummy

                    jj = index(dummy,"xlo") 
                    if (jj/=0) read(dummy,fmt=*) xmin(1),xmax(1)
                    jj = index(dummy,"ylo") 
                    if (jj/=0) read(dummy,fmt=*) xmin(2),xmax(2)
                    jj = index(dummy,"zlo") 
                    if (jj/=0) read(dummy,fmt=*) xmin(3),xmax(3)
                                      
                   
                                                             
                end do
                write(*,fmt='(a,3f12.6,a,3f12.6)') " Lib_XYZFiles::input_lammps info - read box size ",xmin,":",xmax
                xmax = xmax - xmin
                
                
                

                do ii = 1,this%nAtoms
                
!                     box = 0
!                     read(unit=uu,fmt=*,iostat=ioerr) jj,nn,xx,box                   
!                     xx(1:3) = xx(1:3) - xmin(1:3) + box(1:3)*xmax(1:3)
                    
                    read(unit=uu,fmt=*,iostat=ioerr) jj,nn,xx                   
                    if (allokdat .and.(ioerr/=0)) then  
                        print *,"Lib_XYZFiles::input_lammps error - line ",ii," ioerr ",ioerr," read ",jj,nn,xx                        
                    end if
                    allokdat = allokdat .and. (ioerr==0)
                    !xx(1:3) = xx(1:3) - xmin(1:3)
                    
                    this%dat(1:3,ii) = xx(1:3)
                    this%atom(ii) = nn
!                 
!                     read(unit=uu,fmt=*,iostat=ioerr) jj,nn,this%dat(1:this%nColumns,ii)
!                     allokdat = allokdat .and. (ioerr==0)
!                     this%atom(ii) = nn
                    
                    if (op) then
                        if (ii<4) then
                            write(unit=*,fmt='(a8,1000f16.6)',iostat=ioerr) this%atomName(nn),this%dat(1:this%nColumns,ii)
                        else if ( (ii==4).and.(this%nAtoms>5) ) then
                            write(unit=*,fmt='(a)',iostat=ioerr) "..."
                        else if ( (ii>this%nAtoms-3).and.(this%nAtoms>5) ) then
                            write(unit=*,fmt='(a8,1000f16.6)',iostat=ioerr) this%atomName(nn),this%dat(1:this%nColumns,ii)
                        end if
                    end if


                end do
            close(unit=uu)

            if (.not. allokdat) then
                print *,"Lib_XYZFiles::input_lammps warning - may not have read data correctly"                
            end if


            return
        end subroutine input_lammps

        
    !---                

        subroutine output0(this,mask)
    !---^^^^^^^^^^^^^^^^^^^^^^^^
            type(XYZFile),intent(in)     ::      this
            logical,dimension(:),intent(in),optional          ::      mask

            integer                 ::      uu
            integer                 ::      ii,jj,ioerr,mm
            character(len=XYZFILE_FILENAMELENGTH)       ::      f
            
            character(len=256)                              ::      dummy                         
            integer,dimension(:),allocatable                ::      nn 
            real(kind=real64),dimension(1:this%nColumns)    ::      datmax 
            real(kind=real64),dimension(9)                  ::      a_super
            logical                                         ::      triclinic
            integer,dimension(0:8)                          ::      col_type
            type(StringTokenizer)                           ::      stok
           
            uu = findGoodUnit()
            f = getNumberedFilename(this)
            open(unit=uu,file=trim(f)//".tmp",action="write")
            
                if (this%lammpsFormat) then
            
                    write(unit=uu,fmt='(a)') "LAMMPS data file in write_data format written by UKAEA Lib_XYZFiles library"
                    write(unit=uu,fmt='(a)') ""
                    if (present(mask)) then
                        write(unit=uu,fmt='(i10,a)') count(mask)," atoms"
                    else
                        write(unit=uu,fmt='(i10,a)') this%nAtoms," atoms" 
                    end if
                    write(unit=uu,fmt='(i4,a)') this%nAtomNames," atom types"
                    write(unit=uu,fmt='(a)') ""
                    
                    stok = StringTokenizer_ctor(this%column_description,"'""")
                    call nextToken( stok,dummy )
                    call nextToken( stok,dummy )
                    call parse( dummy,a_super,ii )
                    if (ii<9) then
                        print *,"Lib_XYZFiles::output0 error - trying to write lammps write_data() format"
                        print *,"Lib_XYZFiles::output0 error - can't parse supercell from column_description """//trim(this%column_description)//""""
                        stop
                    end if
                    triclinic = (maxval(abs( (/a_super(2),a_super(3),a_super(4),a_super(6),a_super(7),a_super(8)/) )) > 1.0d-8)

                    if (triclinic) then
                        print *,"Lib_XYZFiles::output0 error - trying to write lammps write_data() format"
                        print *,"Lib_XYZFiles::output0 error - haven't coded triclinic cell yet"
                        stop
                    else    
                        write (unit=uu,fmt='(2g24.16,a)') 0.0d0,a_super(1)," xlo xhi"
                        write (unit=uu,fmt='(2g24.16,a)') 0.0d0,a_super(5)," ylo yhi"
                        write (unit=uu,fmt='(2g24.16,a)') 0.0d0,a_super(9)," zlo zhi"
                    end if                
                    write(unit=uu,fmt='(a)') ""
                    write(unit=uu,fmt='(a)') "Masses"
                    write(unit=uu,fmt='(a)') "" 
                    do ii = 1,this%nAtomNames   
                        write(unit=uu,fmt='(i4,f10.4)') ii,getMass( getAtomName(this,ii) )
                    end do
                    write(unit=uu,fmt='(a)') "" 
                    write(unit=uu,fmt='(a)') "Atoms # atomic" 
                    write(unit=uu,fmt='(a)') "" 
                    dummy = ""
                    a_super = 0.0d0
                    col_type = (/ LIBXYZ_COL_STRING,LIBXYZ_COL_UINT,LIBXYZ_COL_UINT,LIBXYZ_COL_FLOAT,LIBXYZ_COL_FLOAT,LIBXYZ_COL_FLOAT,LIBXYZ_COL_UINT,LIBXYZ_COL_UINT,LIBXYZ_COL_UINT /)                    
                    allocate(nn(0:8))
                    nn(0) = 0
                    nn(1) = ceiling(log10(1.000001d0*this%nAtoms))
                    nn(2) = ceiling(log10(1.000001d0*this%nAtomNames)) 
                    datmax = 1.000001d0
                    do ii = 1,this%nAtoms
                        if (present(mask)) then
                            if (.not. mask(ii)) cycle
                        end if
                        do jj = 1,3
                            datmax(jj) = max(datmax(jj),abs(this%dat(jj,ii)))
                        end do
                    end do    
                    do jj = 1,3
                        nn(2+jj) = ceiling(log10(1.000001d0*datmax(jj)))                
                    end do     
                    nn(6:8) = 1
                    
                    mm = 0
                    do ii = 1,this%nAtoms
                        if (present(mask)) then
                            if (.not. mask(ii)) cycle
                        end if
                        mm = mm + 1
                        a_super(1:5) = (/   real(mm,kind=real64),              &
                                            real(this%atom(ii),kind=real64),   &
                                            this%dat(1,ii),this%dat(2,ii),this%dat(3,ii) /) 
                                     
                        dummy = opLine( "",a_super(1:8),nn , col_type ) 
                        write(unit=uu,fmt='(a)') trim(dummy)
                    end do
                    write(unit=uu,fmt='(a)') "" 
                    write(unit=uu,fmt='(a)') "Velocities" 
                    write(unit=uu,fmt='(a)') "" 
                    dummy = ""
                    if (present(mask)) then
                        do ii = 1,count(mask)
                            call utoa(ii,nn(1),dummy,jj)
                            dummy(jj+2:jj+6) = "0 0 0" 
                            write(unit=uu,fmt='(a)') trim(dummy)
                        end do

                    else
                        do ii = 1,this%nAtoms
                            call utoa(ii,nn(1),dummy,jj)
                            dummy(jj+2:jj+6) = "0 0 0" 
                            write(unit=uu,fmt='(a)') trim(dummy)
                        end do
                    end if
                else
                    allocate(nn(0:this%nColumns))
                    nn = 0
                    do ii = 1,this%nAtomNames
                        nn(0) = max(nn(0),len_trim(this%atomName(ii)))
                    end do
                    datmax = 1.000001d0
                    do ii = 1,this%nAtoms
                        if (present(mask)) then
                            if (.not. mask(ii)) cycle
                        end if
                        do jj = 1,this%nColumns
                            datmax(jj) = max(datmax(jj),abs(this%dat(jj,ii)))
                        end do
                    end do
                    do jj = 1,this%nColumns
                        nn(jj) = ceiling(log10(1.000001d0*datmax(jj)))                
                    end do             
                
                    if (this%nHeaderLines>0) then
                        do ii = 1,this%nHeaderLines
                            write(unit=uu,fmt='(a)',iostat=ioerr) trim(this%header(ii))
                        end do
                    end if
                    write(unit=uu,fmt='(i12)',iostat=ioerr) this%nAtoms
    
                    write(unit=uu,fmt='(a)',iostat=ioerr) trim(this%column_description)
                    dummy = ""
                    do ii = 1,this%nAtoms
                        if (present(mask)) then
                            if (.not. mask(ii)) cycle
                        end if
                        dummy = opLine(this%atomName(this%atom(ii)),this%dat(1:this%nColumns,ii),nn , this%col_dataType(0:this%nColumns)) 
                        write(unit=uu,fmt='(a)') trim(dummy)
                       
                    end do
                    
                end if
                
            close(unit=uu)
            call system("mv -fv "//trim(f)//".tmp "//trim(f))
            return
             
!            
!                
        end subroutine output0

        subroutine output1(this,ovito)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
            type(XYZFile),intent(in)    ::      this
            logical,intent(in)          ::      ovito

            integer                 ::      uu
            integer                 ::      ii,jj,ioerr
            character(len=XYZFILE_FILENAMELENGTH)       ::      f
            
            
            character(len=256)                              ::      dummy                         
            integer,dimension(0:this%nColumns)              ::      nn 
            real(kind=real64),dimension(1:this%nColumns)    ::      datmax 
             
            if (.not. ovito) then
                call output0(this)
                return
            end if
           
            nn = 0
            do ii = 1,this%nAtomNames
                nn(0) = max(nn(0),len_trim(this%atomName(ii)))
            end do
            datmax = 1.000001d0
            do ii = 1,this%nAtoms
                do jj = 1,this%nColumns
                    datmax(jj) = max(datmax(jj),abs(this%dat(jj,ii)))
                end do
            end do
            do jj = 1,this%nColumns
                nn(jj) = ceiling(log10(1.000001d0*datmax(jj)))                
            end do         
           
            
            uu = findGoodUnit()
            f = getNumberedFilename(this,ovito)
            open(unit=uu,file=trim(f)//".tmp",action="write")
                write(unit=uu,fmt=*,iostat=ioerr) this%nAtoms

                write(unit=uu,fmt='(a)',iostat=ioerr) trim(this%column_description)

                do ii = 1,this%nAtoms
!                    write(unit=uu,fmt='(a,1000f16.6)',iostat=ioerr) adjustl(trim(this%atomName(this%atom(ii)))),this%dat(1:this%nColumns,ii)
                  dummy = opLine(this%atomName(this%atom(ii)),this%dat(1:this%nColumns,ii),nn) 
                  write(unit=uu,fmt='(a)') trim(dummy)
                end do
            close(unit=uu)
            call system("mv "//trim(f)//".tmp "//trim(f))
            return
            
 
                        
        end subroutine output1

    !     subroutine output2(this,mask)
    ! !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !         type(XYZFile),intent(in)     ::      this
    !         logical,dimension(:),intent(in)          ::      mask

    !         integer                 ::      uu
    !         integer                 ::      ii,jj,ioerr
    !         character(len=XYZFILE_FILENAMELENGTH)       ::      f
            
            
    !         character(len=256)                              ::      dummy                         
    !         integer,dimension(0:this%nColumns)              ::      nn 
    !         real(kind=real64),dimension(1:this%nColumns)    ::      datmax 
              
    !         nn = 0
    !         do ii = 1,this%nAtomNames
    !             nn(0) = max(nn(0),len_trim(this%atomName(ii)))
    !         end do
    !         datmax = 1.000001d0
    !         do ii = 1,this%nAtoms
    !             do jj = 1,this%nColumns
    !                 datmax(jj) = max(datmax(jj),abs(this%dat(jj,ii)))
    !             end do
    !         end do
    !         do jj = 1,this%nColumns
    !             nn(jj) = ceiling(log10(1.000001d0*datmax(jj)))                
    !         end do 
            
    !         uu = findGoodUnit()
    !         f = getNumberedFilename(this)
    !         open(unit=uu,file=trim(f)//".tmp",action="write")
    !             write(unit=uu,fmt=*,iostat=ioerr) count(mask)

    !             write(unit=uu,fmt='(a)',iostat=ioerr) trim(this%column_description)

    !             do ii = 1,this%nAtoms
    !                 if (.not. mask(ii)) cycle
    !                 dummy = opLine(this%atomName(this%atom(ii)),this%dat(1:this%nColumns,ii),nn, this%col_dataType(0:this%nColumns))  
    !                 write(unit=uu,fmt='(a)') trim(dummy)
                    
    !             end do
    !         close(unit=uu)
    !         call system("mv "//trim(f)//".tmp "//trim(f))
    !         return
             
    
                        
    !     end subroutine output2
        

        pure function opLine(t,dat,n,dt) result(line)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      generate a single output line eg
    !*      "atom" 123.456 123.456 123.456 654.321 1.00 1.00 1.00
            character(len=*),intent(in)                     ::      t
            real(kind=real64),dimension(:),intent(in)       ::      dat
            integer,dimension(0:),intent(in)                ::      n
            integer,dimension(0:),intent(in),optional       ::      dt
            character(len=256)                              ::      line
            integer     ::      ii,kk,mm
            
            line = t
            !call ftoa( dat,n(1:),6,line(n(0)+1:),kk )   
            
            if (present(dt)) then
                kk = n(0) 
                do ii = 1,size(dat)
                    if (dt(ii) == LIBXYZ_COL_FLOAT) then
                        call ftoa(dat(ii),n(ii),6,line(kk+1:),mm) 
                        kk = kk + mm
                    else if (dt(ii) == LIBXYZ_COL_UINT) then
                        call utoa(nint(abs(dat(ii))),n(ii),line(kk+2:),mm) 
                        kk = kk + mm + 1
                    else if (dt(ii) == LIBXYZ_COL_INT) then
                        call itoa(nint(dat(ii)),n(ii),line(kk+2:),mm) 
                        kk = kk + mm + 1
                    end if
                    
                end do                
            else                    
                kk = n(0) 
                do ii = 1,size(dat)
                    call ftoa(dat(ii),n(ii),6,line(kk+1:),mm) 
                    kk = kk + mm
                end do
            end if                    
            
            
            return
        end function opLine        

        subroutine outputQhull(this)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      output a file in qhull readable format
            type(XYZFile),intent(in)     ::      this
            integer                 ::      uu
            integer                 ::      ii,ioerr
            character(len=XYZFILE_FILENAMELENGTH)       ::      f
            uu = findGoodUnit()
            f = getNumberedFilename(this)
            open(unit=uu,file=trim(f)//".tmp",action="write")
                if (this%nHeaderLines>0) then
                    do ii = 1,this%nHeaderLines
                        write(unit=uu,fmt='(a)',iostat=ioerr) trim(this%header(ii))
                    end do
                end if
                write(unit=uu,fmt='(i12)',iostat=ioerr) 3
                write(unit=uu,fmt='(i12)',iostat=ioerr) this%nAtoms
                do ii = 1,this%nAtoms
                    write(unit=uu,fmt='(3f16.8)',iostat=ioerr) this%dat(1:3,ii)
                end do
            close(unit=uu)
            call system("mv "//trim(f)//".tmp "//trim(f)//".qhin" )
            return
        end subroutine outputQhull
    !---

        subroutine clear0(this)
    !---^^^^^^^^^^^^^^^^^^^^^^^
    !*      cut all atoms
            type(XYZFile),intent(inout)     ::      this
            this%nAtoms = 0
            return
        end subroutine clear0


        subroutine clear1(this,i)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      cut atom number i
            type(XYZFile),intent(inout)     ::      this
            integer,intent(in)              ::      i
            if (real(i)*(this%nAtoms+1-i)<=0) return
            this%atom(i) = this%atom(this%nAtoms)
            this%dat(1:this%nColumns,i) = this%dat(1:this%nColumns,this%nAtoms)
            this%nAtoms = this%nAtoms - 1
            return
        end subroutine clear1

    !---

        subroutine clone0(this,that)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      this = that deep copy with allocation
            type(XYZFile),intent(inout)     ::      this
            type(XYZFile),intent(in)     ::      that
            call ensureAdequateStorage(this,that%nAtoms,that%nAtomNames,that%nColumns,that%nHeaderLines)
            this%filen                              = that%filen
            this%nHeaderLines                       = that%nHeaderLines
            this%nAtoms                             = that%nAtoms
            this%nColumns                           = that%nColumns
            this%nAtomNames                         = that%nAtomNames
            this%filename                           = that%filename
            this%header(1:this%nHeaderLines)        = that%header(1:this%nHeaderLines)
            this%column_description                 = that%column_description
            this%atomName(1:this%nAtomNames)        = that%atomName(1:this%nAtomNames)
            this%atom(1:this%nAtoms)                = that%atom(1:this%nAtoms)
            this%dat(1:this%nColumns,1:this%nAtoms) = that%dat(1:this%nColumns,1:this%nAtoms)
            return
        end subroutine clone0










    !---

        function getNumberedFilename(this,ovito) result(f)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
            type(XYZFile),intent(in)    ::      this
            logical,intent(in),optional ::      ovito
            character(len=XYZFILE_FILENAMELENGTH)       ::      f
            character(len=5)            ::      aaaa
            integer                     ::      lastdot
            logical     ::  ovitofilename
            if (this%lammpsFormat) then
                f = this%filename
                return
            end if
            ovitofilename = .false. ; if (present(ovito)) ovitofilename = .true.
            lastdot = index(this%filename,".",back=.true.)
            if (this%filen<0) then
                f = this%filename
                if (ovitofilename) f = trim(f)//".ovito"
                if (lastdot == 0) then
                    f = trim(f)//".xyz"
                else
                    if (this%filename(len_trim(f)-3:len_trim(f))/=".xyz") f = trim(f)//".xyz"
                end if
                
            else
                write(aaaa,fmt='(i5)') this%filen ; aaaa = adjustl(aaaa) ; aaaa = repeat("0",5-len_trim(aaaa))//trim(aaaa)
                if (lastdot > 0) then
                    if (ovitofilename) then
                        f = this%filename(1:lastdot)//"ovito."//aaaa//".xyz"
                    else
                        f = this%filename(1:lastdot)//aaaa//".xyz"
                    end if
                else
                    if (ovitofilename) then
                        f = trim(this%filename)//".ovito."//aaaa//".xyz"
                    else
                        f = trim(this%filename)//"."//aaaa//".xyz"
                    end if
                end if
            end if
            return
        end function getNumberedFilename


        pure function getFilename0(this) result(f)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
            type(XYZFile),intent(in)    ::      this
            character(len=XYZFILE_FILENAMELENGTH)       ::      f
            f = this%filename
            return
        end function getFilename0

    !---


        pure function getnAtoms0(this) result(nAtoms)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
            type(XYZFile),intent(in)         ::      this
            integer                             ::      nAtoms
            nAtoms = this%nAtoms
            return
        end function getnAtoms0

        pure function getnAtoms1(this,bytype) result(nAtoms)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
            type(XYZFile),intent(in)        ::      this
            logical,intent(in)              ::      bytype
            integer,dimension(this%nAtomNames)     ::      nAtoms
            integer     ::      ii
            if (.not. bytype) then
                nAtoms = 0
                nAtoms(1) = this%nAtoms
            else
                do ii = 1,this%nAtomNames
                    nAtoms(ii) = count( this%atom(1:this%nAtoms) == ii )
                end do
            end if
            return
        end function getnAtoms1


        pure function getnColumns0(this) result(nColumns)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
            type(XYZFile),intent(in)         ::      this
            integer                             ::      nColumns
            nColumns = this%nColumns
            return
        end function getnColumns0



        pure function getnHeaderLines0(this) result(nHeaderLines)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
            type(XYZFile),intent(in)         ::      this
            integer                          ::      nHeaderLines
            nHeaderLines = this%nHeaderLines
            return
        end function getnHeaderLines0

    !---

        pure function getnAtomNames0(this) result(nAtomNames)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
            type(XYZFile),intent(in)         ::      this
            integer                             ::      nAtomNames
            nAtomNames = this%nAtomNames
            return
        end function getnAtomNames0

        pure function getAtomName0(this,i) result(atom)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      return the name associated with type i
            type(XYZFile),intent(in)        ::      this
            integer,intent(in)              ::      i
            character(len=XYZFILE_ATOMNAMELENGTH)       ::      atom
            atom = ""
            if ( i*(this%nAtomNames+1-i) > 0 ) atom = this%atomName(i)
            return
        end function getAtomName0

        pure function getAtomTypeFromName(this,atom) result(i)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      do you recognise this name? return a type. 0 if not recognised
            type(XYZFile),intent(in)        ::      this
            character(len=*),intent(in)     ::      atom
            integer                         ::      i
            integer         ::      jj

            i = 0
            do jj = 1,this%nAtomNames
                if (trim(atom) == trim(this%atomName(jj))) then
                    i = jj
                    return
                end if
            end do

            return
        end function getAtomTypeFromName

        pure function getBoxRepeats( this,a0 ) result(Nx)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      assuming a cubic cellside a0, how many box repeats?
            type(XYZFile),intent(in)            ::      this
            real(kind=real64),intent(in)        ::      a0
            integer,dimension(3)                ::      Nx

            real(kind=real64),dimension(3)      ::      xmin,xmax
            call getMinMax( this,xmin,xmax )
            xmax(1:3) = (xmax(1:3) - xmin(1:3) + 0.25d0)
            Nx(1:3) = nint( xmax(1:3)/a0 )
            return
        end function getBoxRepeats

        subroutine setnAtomNames0(this,natomNames) 
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
            type(XYZFile),intent(inout)       ::      this
            integer,intent(in)                ::      natomNames
            character(len=XYZFILE_ATOMNAMELENGTH),dimension(:),pointer       ::      newatomName            
            if (associated(this%atomName)) then
                if (natomNames > this%natomNames) then
                    if (this%natomNames > 0) then
                        allocate(newatomName(natomNames))
                        newatomName(1:this%natomNames) = this%atomName(1:this%natomNames)
                        newatomName(this%natomNames+1) = ""
                        deallocate(this%atomName)
                        this%atomName => newatomName
                    else
                        deallocate(this%atomName)
                        allocate(this%atomName(natomNames))
                        this%atomName = ""
                    end if
                end if
            else
                allocate(this%atomName(natomNames))
                this%atomName = ""
            end if
            
            this%natomNames = natomNames
            return
        end subroutine setnAtomNames0


        subroutine setnHeaderLines0(this,nHeaderLines)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
            type(XYZFile),intent(inout)       ::      this
            integer,intent(in)                ::      nHeaderLines
            character(len=XYZFILE_COMMENTLINELENGTH),dimension(:),pointer       ::      newHeader
            if (associated(this%header)) then
                if (nHeaderLines > size(this%header)) then
                    if (this%nHeaderLines > 0) then
                        allocate(newHeader(nHeaderLines))
                        newHeader(1:this%nHeaderLines) = this%header(1:this%nHeaderLines)
                        newHeader(this%nHeaderLines+1) = ""
                        deallocate(this%header)
                        this%header => newHeader
                    else
                        deallocate(this%header)
                        allocate(this%header(nHeaderLines))
                        this%header = ""
                    end if
                end if
            else
                allocate(this%header(nHeaderLines))
                this%header = ""
            end if
            this%nHeaderLines = nHeaderLines
            return
        end subroutine setnHeaderLines0

        subroutine setnColumns0(this,nColumns,nAtoms)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
            type(XYZFile),intent(inout)       ::      this
            integer,intent(in)                ::      nColumns
            integer,intent(in),optional       ::      nAtoms
            real(kind=real64),dimension(:,:),pointer        ::      dat_tmp
            integer,dimension(:),pointer                    ::      atom_tmp,dt_tmp
            
            
            real            ::      mm
            
            if (present(nAtoms)) then
                !print *,"Lib_XYZFiles::setnColumns0 info - allocating ",nColumns," columns for ",nAtoms," atoms"
                mm = nColumns*8 + 4
                mm = mm * nAtoms / (1024.0*1024.0)
                !print *,"Lib_XYZFiles::setnColumns0 info - memory usage ",mm," Mb"
                
                
                if (associated(this%dat)) then
                    if ( (nColumns > size(this%dat,dim=1)).or.(nAtoms > size(this%dat,dim=2)) ) then
                        if (this%nColumns*this%nAtoms > 0) then
                            allocate(dat_tmp(nColumns,nAtoms))
                            dat_tmp(1:this%nColumns,1:this%nAtoms) = this%dat(1:this%nColumns,1:this%nAtoms)
                            dat_tmp(this%nColumns+1:,:) = 0.0d0
                            dat_tmp(:,this%nAtoms+1:) = 0.0d0
                            deallocate(this%dat)
                            this%dat => dat_tmp
                            
                            allocate(dt_tmp(0:nColumns))
                            dt_tmp(0:this%nColumns) = this%col_dataType(0:this%nColumns)
                            dt_tmp(this%nColumns+1:) = LIBXYZ_COL_FLOAT
                            if (associated(this%col_dataType))  deallocate(this%col_dataType)
                            this%col_dataType => dt_tmp
                            
                        else
                            deallocate(this%dat)
                            allocate(this%dat(nColumns,nAtoms))
                            this%dat = 0.0d0
                            
                            if (associated(this%col_dataType))  deallocate(this%col_dataType)
                            allocate(this%col_dataType(0:nColumns))
                            this%col_dataType(0) = LIBXYZ_COL_STRING
                            this%col_dataType(1:) = LIBXYZ_COL_FLOAT                            
                        end if
                    end if
                else
                    allocate(this%dat(nColumns,nAtoms))
                    this%dat = 0.0d0
                    allocate(this%col_dataType(0:nColumns))
                    this%col_dataType(0) = LIBXYZ_COL_STRING
                    this%col_dataType(1:) = LIBXYZ_COL_FLOAT          
                end if

                if (associated(this%atom)) then
                    if (nAtoms>size(this%atom)) then
                        if (this%nAtoms>0) then
                            allocate(atom_tmp(nAtoms))
                            atom_tmp(1:this%nAtoms) = this%atom(1:this%nAtoms)
                            atom_tmp(this%nAtoms+1:) = 0
                            deallocate(this%atom)
                            this%atom => atom_tmp
                        else
                            deallocate(this%atom)
                            allocate(this%atom(nAtoms))
                            this%atom = 0
                        end if
                    end if
                else
                    allocate(this%atom(nAtoms))
                    this%atom = 0
                end if


                this%nColumns = nColumns
                this%nAtoms = nAtoms
            else
                print *,"Lib_XYZFiles::setnColumns0 info - allocating ",nColumns," columns for ",this%nAtoms," atoms"
                mm = nColumns*8 + 4
                mm = mm * this%nAtoms / (1024.0*1025.0)
                print *,"Lib_XYZFiles::setnColumns0 info - memory usage ",mm," Mb"
                
            
                if (associated(this%dat)) then
                    if (nColumns > size(this%dat,dim=1)) then
                        if (this%nColumns > 0) then
                            allocate(dat_tmp(nColumns,this%nAtoms))
                            dat_tmp(1:this%nColumns,1:this%nAtoms) = this%dat(1:this%nColumns,1:this%nAtoms)
                            dat_tmp(this%nColumns+1:,:) = 0.0d0
                            deallocate(this%dat)
                            this%dat => dat_tmp
                            
                            allocate(dt_tmp(0:nColumns))
                            dt_tmp(0:this%nColumns) = this%col_dataType(0:this%nColumns)
                            dt_tmp(this%nColumns+1:) = LIBXYZ_COL_FLOAT
                            if (associated(this%col_dataType))  deallocate(this%col_dataType)
                            this%col_dataType => dt_tmp
                            
                            
                        else
                            deallocate(this%dat)
                            allocate(this%dat(nColumns,this%nAtoms))
                            this%dat = 0.0d0
                            
                            if (associated(this%col_dataType))  deallocate(this%col_dataType)
                            allocate(this%col_dataType(0:nColumns))
                            this%col_dataType(0) = LIBXYZ_COL_STRING
                            this%col_dataType(1:) = LIBXYZ_COL_FLOAT        
                        end if
                    end if
                else
                    allocate(this%dat(nColumns,this%nAtoms))
                    this%dat = 0.0d0
                    allocate(this%col_dataType(0:nColumns))
                    this%col_dataType(0) = LIBXYZ_COL_STRING
                    this%col_dataType(1:) = LIBXYZ_COL_FLOAT        
                end if
                this%nColumns = nColumns
            end if
            
            if (this%nColumns>3) this%lammpsFormat = .false.
            return
        end subroutine setnColumns0
        
        
        function getColumn_Description(this) result(column_description)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
            type(XYZFile),intent(in)                        ::      this
            character(len=XYZFILE_COMMENTLINELENGTH)        ::      column_description
            column_description = this%column_description 
            return
        end function getColumn_Description


        subroutine setColumn_Description0(this,column_description)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      set a general column descriptor line
            type(XYZFile),intent(inout)     ::      this
            character(len=*),intent(in)     ::      column_description
            this%column_description = trim(column_description)
            return
        end subroutine setColumn_Description0


        subroutine setColumn_Description1(this,a,desc)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      set an extended .xyz column description line using lattice vectors a
            type(XYZFile),intent(inout)     ::      this
            real(kind=real64),dimension(3,3),intent(in)     ::      a
            character(len=*),intent(in),optional            ::      desc
            character(len=XYZFILE_COMMENTLINELENGTH)        ::      dummy
            write(dummy,fmt='(9f16.8)') a           ;   dummy = adjustl(dummy)
            this%column_description = "Lattice="""//trim(dummy)//""" Properties=species:S:1:pos:R:3"
            if (present(desc)) this%column_description = trim(this%column_description)//trim(desc)
            return
        end subroutine setColumn_Description1        
        
        subroutine setColumn_Description2(this,Nx,Ny,Nz,a,desc)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      set an extended .xyz column description line using unit cell vectors a
            type(XYZFile),intent(inout)     ::      this
            integer,intent(in)              ::      Nx,Ny,Nz
            real(kind=real64),dimension(3,3),intent(in)    ::      a
            character(len=*),intent(in),optional            ::      desc
            character(len=XYZFILE_COMMENTLINELENGTH)        ::      dummy
            write(dummy,fmt='(9f16.8)') a(1:3,1)*Nx,a(1:3,2)*Ny,a(1:3,3)*Nz     ;   dummy = adjustl(dummy)
            this%column_description = "Lattice="""//trim(dummy)//""" Properties=species:S:1:pos:R:3"
            if (present(desc)) this%column_description = trim(this%column_description)//trim(desc)
            return
        end subroutine setColumn_Description2        
        
        subroutine setColumn_Description3(this,a,origin,pbc,desc)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      set an extended .xyz column description line using lattice vectors a, a box origin and pbc
            type(XYZFile),intent(inout)     ::      this
            real(kind=real64),dimension(3,3),intent(in)     ::      a
            real(kind=real64),dimension(3),intent(in)       ::      origin
            logical,dimension(3),intent(in),optional        ::      pbc
            character(len=*),intent(in),optional            ::      desc
            character(len=XYZFILE_COMMENTLINELENGTH)        ::      dummy,dummy2,dummy3
            write(dummy,fmt='(9f16.8)') a                       ;   dummy = adjustl(dummy)
            write(dummy2,fmt='(3f16.8)') origin                 ;   dummy2 = adjustl(dummy2)
            dummy3 = "T T T"
            if (present(pbc)) write(dummy3,fmt='(3l2)') pbc     ;   dummy3 = adjustl(dummy3)
            this%column_description = "Lattice="""//trim(dummy)//""" Origin="""//trim(dummy2)//""" pbc="""//trim(dummy3)//""" Properties=species:S:1:pos:R:3"
            if (present(desc)) this%column_description = trim(this%column_description)//trim(desc)
            return
        end subroutine setColumn_Description3        
        
        
        subroutine setFilename(this,filename)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
            type(XYZFile),intent(inout)     ::      this
            character(len=*),intent(in)     ::      filename
            this%filename = trim(filename)
            return
        end subroutine setFilename


        subroutine setFilen(this,filen)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
            type(XYZFile),intent(inout)     ::      this
            integer,intent(in)              ::      filen
            this%filen = filen
            return
        end subroutine setFilen



        subroutine setAtomNames0(this,atom_names)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
            type(XYZFile),intent(inout)       ::      this
            character(len=*),dimension(:),intent(in)       ::      atom_names
            integer     ::      nAtomNames
            nAtomNames = size(atom_names)
            if (this%nAtomNames>0) deallocate(this%atomName)
            allocate(this%atomName(nAtomNames))
            this%nAtomNames = nAtomNames
            this%atomName(1:this%nAtomNames) = adjustl( atom_names(1:this%nAtomNames) )
            return
        end subroutine setAtomNames0

        subroutine setAtomNames0a(this,atom_names)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
            type(XYZFile),intent(inout)       ::      this
            character(len=*),intent(in)       ::      atom_names
            if (this%nAtomNames>0) deallocate(this%atomName)
            allocate(this%atomName(1))
            this%nAtomNames = 1
            this%atomName(1) = adjustl( atom_names )
            return
        end subroutine setAtomNames0a




        subroutine setnAtoms(this,nAtoms)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
            type(XYZFile),intent(inout)     ::      this
            integer,intent(in)              ::      nAtoms
            integer,dimension(:),pointer                    ::      atom_tmp
            real(kind=real64),dimension(:,:),pointer        ::      dat_tmp
            if (associated(this%dat)) then
                if (nAtoms>size(this%dat,dim=2)) then
                    if (this%nAtoms>0) then
                        allocate(dat_tmp(this%nColumns,nAtoms))
                        dat_tmp(1:this%nColumns,1:this%nAtoms) = this%dat(1:this%nColumns,1:this%nAtoms)
                        dat_tmp(:,this%nAtoms+1:) = 0.0d0
                        deallocate(this%dat)
                        this%dat => dat_tmp
                    else
                        deallocate(this%dat)
                        allocate(this%dat(this%nColumns,nAtoms))
                        this%dat = 0.0d0
                    end if
                end if
            else
                allocate(this%dat(this%nColumns,nAtoms))
                this%dat = 0.0d0
            end if
            if (associated(this%atom)) then
                if (nAtoms>size(this%atom)) then
                    if (this%nAtoms>0) then
                        allocate(atom_tmp(nAtoms))
                        atom_tmp(1:this%nAtoms) = this%atom(1:this%nAtoms)
                        atom_tmp(this%nAtoms+1:) = 0
                        deallocate(this%atom)
                        this%atom => atom_tmp
                    else
                        deallocate(this%atom)
                        allocate(this%atom(nAtoms))
                        this%atom = 0
                    end if
                end if
            else
                allocate(this%atom(nAtoms))
                this%atom = 0
            end if
            this%nAtoms = nAtoms
            return
        end subroutine setnAtoms

    !---

        pure subroutine setAtomTypes0(this,t)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
            type(XYZFile),intent(inout)      ::      this
            integer,dimension(:),intent(in)  ::      t
            this%atom(1:this%nAtoms) = t(1:this%nAtoms)
            return
        end subroutine setAtomTypes0

        pure subroutine setAtomTypes1(this,t)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
            type(XYZFile),intent(inout)      ::      this
            integer,intent(in)               ::      t
            this%atom(1:this%nAtoms) = t
            return
        end subroutine setAtomTypes1


        pure subroutine getAtomTypes0(this,t)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
            type(XYZFile),intent(in)            ::      this
            integer,dimension(:),intent(out)    ::      t
            t(1:this%nAtoms) = this%atom(1:this%nAtoms)
            return
        end subroutine getAtomTypes0

        pure function getAtomType0(this,i) result(t)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
            type(XYZFile),intent(in)            ::      this
            integer,intent(in)                  ::      i
            integer                             ::      t
            t = this%atom(i)
            return
        end function getAtomType0

        pure subroutine setAtomType0(this,i,t)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
            type(XYZFile),intent(inout)         ::      this
            integer,intent(in)                  ::      i
            integer,intent(in)                  ::      t
            this%atom(i) = t
            return
        end subroutine setAtomType0

    !---


        pure subroutine setColumns0(this,c)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      can be used to set first 3 cols ( position )
            type(XYZFile),intent(inout)      ::      this
            real(kind=real64),dimension(:,:),intent(in)  ::      c
            integer     ::  nc,na
            nc = min(this%nColumns,size(c,dim=1))
            na = min(this%nAtoms,size(c,dim=2))
            this%dat(1:nc,1:na) = c(1:nc,1:na)
            return
        end subroutine setColumns0


        pure subroutine setColumns1(this,i,c)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      can be used to set first 3 cols ( position )
            type(XYZFile),intent(inout)                 ::      this
            integer,intent(in)                          ::      i
            real(kind=real64),dimension(:),intent(in)   ::      c
            integer     ::  nc
            if (real(i)*(this%nAtoms+1-i)<=0) return
            nc = min(this%nColumns,size(c,dim=1))
            this%dat(1:nc,i) = c(1:nc)
            return
        end subroutine setColumns1
        
        pure subroutine setColumns2(this,i,j,c)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      can be used to set jth col
            type(XYZFile),intent(inout)                 ::      this
            integer,intent(in)                          ::      i,j
            real(kind=real64),intent(in)   ::      c
     
            this%dat(j,i) = c
            return
        end subroutine setColumns2

        pure subroutine getColumns0(this,c)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      can be used to get first 3 cols ( position )
            type(XYZFile),intent(in)                        ::      this
            real(kind=real64),dimension(:,:),intent(out)    ::      c
            integer     ::  nc,na
            nc = min(this%nColumns,size(c,dim=1))
            na = min(this%nAtoms,size(c,dim=2))
            c(1:nc,1:na) = this%dat(1:nc,1:na)
            return
        end subroutine getColumns0

        pure subroutine getColumns1(this,i,c)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      can be used to get first 3 cols ( position )
            type(XYZFile),intent(inout)                 ::      this
            integer,intent(in)                          ::      i
            real(kind=real64),dimension(:),intent(out)  ::      c
            integer     ::  nc
            if (real(i)*(this%nAtoms+1-i)<=0) return
            nc = min(this%nColumns,size(c,dim=1))
            c(1:nc) = this%dat(1:nc,i)
            return
        end subroutine getColumns1
        
        subroutine getColumnsp(this,cp)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      pointer to columns
            type(XYZFile),intent(inout)                 ::      this
            real(kind=real64),dimension(:,:),pointer    ::      cp
            cp => this%dat
            return
        end subroutine getColumnsp


        subroutine getTypesp(this,tp)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      pointer to atom types
            type(XYZFile),intent(inout)                 ::      this
            integer,dimension(:),pointer                ::      tp
            tp => this%atom
            return
        end subroutine getTypesp
    
    !---
    
        

        pure function getCol_dataType0(this,i) result(dt)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      return data type of column i
            type(XYZFile),intent(in)                        ::      this
            integer,intent(in)                              ::      i
            integer                                         ::      dt
            dt = LIBXYZ_COL_FLOAT               
            if ( i*(this%nColumns - i) >= 0) dt = this%col_dataType(i)
            return
        end function getCol_dataType0


        pure function getCol_dataType1(this) result(dt)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ 
    !*      return data type of columns
            type(XYZFile),intent(in)                        ::      this
            integer,dimension(0:this%nColumns)              ::      dt
            dt(0:this%nColumns) = this%col_dataType(0:this%nColumns)
            return
        end function getCol_dataType1


        pure subroutine setCol_dataType0(this,i,dt)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      set data type of column i
            type(XYZFile),intent(inout)                    ::      this
            integer,intent(in)                             ::      i
            integer,intent(in)                             ::      dt
            if ( i*(this%nColumns - i) >= 0) this%col_dataType(i) = dt
            return
        end subroutine setCol_dataType0


        pure subroutine setCol_dataType1(this,dt)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      set data type of columns
            type(XYZFile),intent(inout)                    ::      this
            integer,dimension(0:this%nColumns),intent(in)  ::      dt
            this%col_dataType(0:this%nColumns) = dt(0:this%nColumns)
            return
        end subroutine setCol_dataType1

        pure subroutine setCol_dataType2(this,dt)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      set data type of all columns
            type(XYZFile),intent(inout)                     ::      this
            integer,intent(in)                              ::      dt
            this%col_dataType(0) = LIBXYZ_COL_STRING
            this%col_dataType(0:this%nColumns) = dt 
            return
        end subroutine setCol_dataType2



    !---

        pure subroutine getHeaderLines0(this,h)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      return all header lines
            type(XYZFile),intent(in)                        ::      this
            character(len=*),dimension(:),intent(out)       ::      h
            h = ""
            if (this%nHeaderLines==0) return
            h(1:this%nHeaderLines) = this%header(1:this%nHeaderLines)
            return
        end subroutine getHeaderLines0

        pure subroutine getHeaderLines1(this,i,h)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      return one header line
            type(XYZFile),intent(in)                        ::      this
            integer,intent(in)                              ::      i
            character(len=*),intent(out)                    ::      h
            h = ""
            if (i*(this%nHeaderLines+1-i)<=0) return
            h = this%header(i)
            return
        end subroutine getHeaderLines1


    !---

        pure subroutine setHeaderLines0(this,h)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      set all header lines
            type(XYZFile),intent(inout)                     ::      this
            character(len=*),dimension(:),intent(in)        ::      h
            integer     ::  nn,ii
            if (this%nHeaderLines==0) return
            nn = min(this%nHeaderLines,size(h))
            do ii = 1,nn
                this%header(ii) = h(ii)
            end do
            return
        end subroutine setHeaderLines0

        pure subroutine setHeaderLines0a(this,h)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      set all header lines
            type(XYZFile),intent(inout)                     ::      this
            character(len=*),intent(in)                     ::      h
            if (this%nHeaderLines==0) return
            this%header(1) = h
            if (this%nHeaderLines>1) this%header(2:) = ""
            return
        end subroutine setHeaderLines0a


        pure subroutine setHeaderLines1(this,i,h)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      set one header line
            type(XYZFile),intent(inout)                     ::      this
            integer,intent(in)                              ::      i
            character(len=*),intent(in)                     ::      h
            if (i*(this%nHeaderLines+1-i)<=0) return
            this%header(i) = h
            return
        end subroutine setHeaderLines1


        subroutine setHeaderLines2(this,nx,ny,nz,a)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      append additional header lines to report supercell data
            type(XYZFile),intent(inout)                     ::      this
            integer,intent(in)                              ::      nx,ny,nz
            real(kind=real64),dimension(3,3),intent(in)     ::      a
            
            character(len=XYZFILE_COMMENTLINELENGTH)    ::      dummy
            
            integer     ::      jj,kk,nn
            real(kind=real64),dimension(9)      ::          aaaa
            logical     ::      ok
            logical     ::      supercellAsOneLine = .false.
            
        
            do nn = 1,this%nHeaderLines
                dummy = this%header(nn)
                jj = index(dummy,"supercell") 
                if (jj/=0) then
                    if (LIBXYZ_DBG) print *,"Lib_XYZFiles::setHeaderLines2 info - found existing supercell data? checking ..." 
                    dummy = dummy(jj+9:)
                    ok = .false.                                                    
                    call parse( dummy,aaaa,jj )
                    if (jj==3) then
                        if (LIBXYZ_DBG) print *,"Lib_XYZFiles::setHeaderLines2 info - parsing """//trim(dummy)//""" as cell repeats "
                    else
                        print *,"Lib_XYZFiles::setHeaderLines2 warning  - can't parse """//trim(dummy)//""" as cell repeats. Assuming 1x1x1 unit cells..."
                    end if
                    
                    do kk = 1,3
                        if (nn+kk>this%nHeaderLines) exit                                                        
                        dummy = this%header(nn+kk) 
                        call parse( dummy(3:),aaaa,jj )
                        if (jj==3) then
                            if (LIBXYZ_DBG) write(*,fmt='(a)') " Lib_XYZFiles::setHeaderLines2 info - parsing """//trim(dummy(3:))//""" as lattice vectors "
                            if (kk==3) then 
                                ok = .true.  
                                supercellAsOneLine = .false.                                  
                            end if
                        else if (jj==9) then
                            if (LIBXYZ_DBG) write(*,fmt='(a)') " Lib_XYZFiles::setHeaderLines2 info - parsing """//trim(dummy(3:))//""" as lattice vectors "
                            supercellAsOneLine = .true.
                            ok = .true.
                            exit
                        else 
                            print *,"Lib_XYZFiles::setHeaderLines2 error - can't parse """//trim(dummy)//""" as cell vectors."
                            ok = .false.                                    
                            exit
                        end if
                    end do
                    
                    if (ok) then
                        if (LIBXYZ_DBG) print *,"Lib_XYZFiles::setHeaderLines2 info - verified existing supercell data, replacing"
                        write(dummy,fmt='(4(a,i4))') "# supercell ",nx,",",ny,",",nz," repeats"
                        this%header(nn) = trim(dummy)
                        if (supercellAsOneLine) then
                            write(dummy,fmt='(a2,9f14.8)') "# ",a(:,:) ; this%header(nn+1) = trim(dummy)
                        else
                            write(dummy,fmt='(a2,3f14.8)') "# ",a(1,:) ; this%header(nn+1) = trim(dummy)
                            write(dummy,fmt='(a2,3f14.8)') "# ",a(2,:) ; this%header(nn+2) = trim(dummy)
                            write(dummy,fmt='(a2,3f14.8)') "# ",a(3,:) ; this%header(nn+3) = trim(dummy)       
                        end if
                        if (LIBXYZ_DBG) then
                            print *,"Lib_XYZFiles::setHeaderLines2 info - replacement lines"
                            do kk = nn,nn+3
                                print *,trim(this%header(kk))
                            end do
                        end if
                        return          
                    end if
                    
                end if                        
            end do            
            
            nn = getNHeaderLines(this)
            call setNHeaderLines(this,nn+4)
            write(dummy,fmt='(4(a,i4))') "# supercell = ",nx,",",ny,",",nz," repeats"
            this%header(nn+1) = trim(dummy)
            write(dummy,fmt='(a2,3f14.8)') "# ",a(1,:) ; this%header(nn+2) = trim(dummy)
            write(dummy,fmt='(a2,3f14.8)') "# ",a(2,:) ; this%header(nn+3) = trim(dummy)
            write(dummy,fmt='(a2,3f14.8)') "# ",a(3,:) ; this%header(nn+4) = trim(dummy)           

            print *,"Lib_XYZFiles::setHeaderLines2 info - new lines"
            do kk = nn+1,nn+4
                print *,trim(this%header(kk))
            end do

            
            return
        end subroutine setHeaderLines2

    !---

        pure subroutine displace0(this,dx)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      displace all atoms by fixed vector dx
            type(XYZFile),intent(inout)                     ::      this
            real(kind=real64),dimension(3),intent(in)       ::      dx
            integer     ::      ii
            do ii = 1,this%nAtoms
                this%dat(1:3,ii) = this%dat(1:3,ii) + dx(1:3)
            end do
            return
        end subroutine displace0

        subroutine shuffle0(this,indx)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      shuffle the ordering of the atoms
    !*      this%dat(1:this%nColumns,ii) = old_dat(1:this%nColumns,indx(ii))
            type(XYZFile),intent(inout)                     ::      this
            integer,dimension(:),intent(in)                 ::      indx
            real(kind=real64),dimension(this%nColumns,this%nAtoms)  ::      old_dat
            integer,dimension(this%nAtoms)  ::      old_atom
            integer             ::      ii
            old_dat(1:this%nColumns,1:this%nAtoms) = this%dat(1:this%nColumns,1:this%nAtoms)
            old_atom(1:this%nAtoms) = this%atom(1:this%nAtoms)
            do ii = 1,this%nAtoms
                this%dat(1:this%nColumns,ii) = old_dat(1:this%nColumns,indx(ii))
                this%atom(ii) = old_atom(indx(ii))
            end do            
            return
        end subroutine shuffle0

    !---
        function findOffset(this,a0,offset) result(xx)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      from the positions x,
    !*      compute the offset to line up left front bottom with cell(1/8,1/8,1/8) or (offset,offset,offset)
    !*      but do not move the atoms
            type(XYZFile),intent(in)            ::      this
            real(kind=real64),intent(in)        ::      a0      !   indicative unit cell side
            real(kind=real64),intent(in),optional   ::  offset
            real(kind=real64),dimension(3)      ::      minx,maxx,xx,xoffbar
            integer                             ::      ii
            integer                             ::      nsx,nsy,nsz
         
           
            minx(1:3) = huge(1.0)
            maxx(1:3) = -huge(1.0)
            do ii = 1,this%nAtoms
                xx(1:3) = this%dat(1:3,ii)
                minx(1:3) = min(xx(1:3),minx(1:3))
                maxx(1:3) = max(xx(1:3),maxx(1:3))
            end do

        !   find atoms on surface, find their average offset
            xoffbar(1:3) = 0.0d0 ; nsx = 0; nsy = 0; nsz = 0;
            do ii = 1,this%nAtoms
                xx(1:3) = this%dat(1:3,ii) - minx(1:3)
                if (4*xx(1) <= a0) then
                    xoffbar(1) = xoffbar(1) + xx(1) ; nsx = nsx + 1
                end if
                if (4*xx(2) <= a0) then
                    xoffbar(2) = xoffbar(2) + xx(2) ; nsy = nsy + 1
                end if
                if (4*xx(3) <= a0) then
                    xoffbar(3) = xoffbar(3) + xx(3) ; nsz = nsz + 1
                end if
            end do
            xoffbar(1) = xoffbar(1) / max(nsx,1)
            xoffbar(2) = xoffbar(2) / max(nsy,1)
            xoffbar(3) = xoffbar(3) / max(nsz,1)
            xx(1:3) = - minx(1:3) - xoffbar(1:3)
            if (present(offset)) then
                xx(1:3) = xx(1:3) + a0*offset
            else
                xx(1:3) = xx(1:3) + a0/8
            end if
        !   now have offset required to line up atoms on 0,0,0 boundary

            return
        end function findOffset


        subroutine adjustOffset0(this,a0,offset,verbose)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      from the positions x,
    !*      compute the offset to line up left front bottom with cell(1/8,1/8,1/8) or (offset,offset,offset)

            type(XYZFile),intent(inout)         ::      this
            real(kind=real64),intent(in)        ::      a0      !   indicative unit cell side
            real(kind=real64),intent(in),optional   ::  offset
            logical,intent(in),optional         ::      verbose
            real(kind=real64),dimension(3)      ::      minx,maxx,xx,xoffbar
            integer                             ::      ii
            integer                             ::      nsx,nsy,nsz
            logical                             ::      op

            op = .false. ; if (present(verbose)) op = verbose
            minx(1:3) = huge(1.0)
            maxx(1:3) = -huge(1.0)
            do ii = 1,this%nAtoms
                xx(1:3) = this%dat(1:3,ii)
                minx(1:3) = min(xx(1:3),minx(1:3))
                maxx(1:3) = max(xx(1:3),maxx(1:3))
            end do

            if (op) then
                print *,"Lib_XYZFiles::findOffset info - indicative lattice parameter ",a0
                write (*,fmt='(a12,3f16.8)') "  minx   = ",minx(1:3)
                write (*,fmt='(a12,3f16.8)') "  maxx   = ",maxx(1:3)
            end if

        !   find atoms on surface, find their average offset
            xoffbar(1:3) = 0.0d0 ; nsx = 0; nsy = 0; nsz = 0;
            do ii = 1,this%nAtoms
                xx(1:3) = this%dat(1:3,ii) - minx(1:3)
                if (4*xx(1) <= a0) then
                    xoffbar(1) = xoffbar(1) + xx(1) ; nsx = nsx + 1
                end if
                if (4*xx(2) <= a0) then
                    xoffbar(2) = xoffbar(2) + xx(2) ; nsy = nsy + 1
                end if
                if (4*xx(3) <= a0) then
                    xoffbar(3) = xoffbar(3) + xx(3) ; nsz = nsz + 1
                end if
            end do
            xoffbar(1) = xoffbar(1) / max(nsx,1)
            xoffbar(2) = xoffbar(2) / max(nsy,1)
            xoffbar(3) = xoffbar(3) / max(nsz,1)
            xx(1:3) = - minx(1:3) - xoffbar(1:3)
            if (present(offset)) then
                xx(1:3) = xx(1:3) + a0*offset
            else
                xx(1:3) = xx(1:3) + a0/8
            end if
        !   now have offset required to line up atoms on 0,0,0 boundary

            if (op) then
                write (*,fmt='(a12,3i6)')    "  surf # = ",nsx,nsy,nsz
                write (*,fmt='(a12,3f16.8)') "  offset = ",xoffbar(1:3)
                write (*,fmt='(a12,3f16.8)') "  adjust = ",xx(1:3)
            end if

            do ii = 1,this%nAtoms
                this%dat(1:3,ii) = this%dat(1:3,ii) + xx(1:3)
            end do

            if (op) then
                minx(1:3) = huge(1.0)
                maxx(1:3) = -huge(1.0)
                do ii = 1,this%nAtoms
                    xx(1:3) = this%dat(1:3,ii)
                    minx(1:3) = min(xx(1:3),minx(1:3))
                    maxx(1:3) = max(xx(1:3),maxx(1:3))
                end do
                print *,"Lib_XYZFiles::findOffset info - after"
                write (*,fmt='(a12,3f16.8)') "  minx   = ",minx(1:3)
                write (*,fmt='(a12,3f16.8)') "  maxx   = ",maxx(1:3)
            end if


            return
        end subroutine adjustOffset0


        subroutine adjustOffset1(this,a0,offset,verbose)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      from the positions x,
    !*      compute the offset to line up left front bottom with cell(1/8,1/8,1/8) or (offset,offset,offset)

            type(XYZFile),intent(inout)         ::      this
            real(kind=real64),dimension(3),intent(in)        ::      a0      !   indicative unit cell sides
            real(kind=real64),intent(in),optional   ::  offset
            logical,intent(in),optional         ::      verbose
            real(kind=real64),dimension(3)      ::      minx,maxx,xx,xoffbar
            integer                             ::      ii
            integer                             ::      nsx,nsy,nsz
            logical                             ::      op

            op = .false. ; if (present(verbose)) op = verbose
            minx(1:3) = huge(1.0)
            maxx(1:3) = -huge(1.0)
            do ii = 1,this%nAtoms
                xx(1:3) = this%dat(1:3,ii)
                minx(1:3) = min(xx(1:3),minx(1:3))
                maxx(1:3) = max(xx(1:3),maxx(1:3))
            end do

            if (op) then
                print *,"Lib_XYZFiles::findOffset info - indicative lattice parameters ",a0
                write (*,fmt='(a12,3f16.8)') "  minx   = ",minx(1:3)
                write (*,fmt='(a12,3f16.8)') "  maxx   = ",maxx(1:3)
            end if

        !   find atoms on surface, find their average offset
            xoffbar(1:3) = 0.0d0 ; nsx = 0; nsy = 0; nsz = 0;
            do ii = 1,this%nAtoms
                xx(1:3) = this%dat(1:3,ii) - minx(1:3)
                if (4*xx(1) <= a0(1)) then
                    xoffbar(1) = xoffbar(1) + xx(1) ; nsx = nsx + 1
                end if
                if (4*xx(2) <= a0(2)) then
                    xoffbar(2) = xoffbar(2) + xx(2) ; nsy = nsy + 1
                end if
                if (4*xx(3) <= a0(3)) then
                    xoffbar(3) = xoffbar(3) + xx(3) ; nsz = nsz + 1
                end if
            end do
            xoffbar(1) = xoffbar(1) / max(nsx,1)
            xoffbar(2) = xoffbar(2) / max(nsy,1)
            xoffbar(3) = xoffbar(3) / max(nsz,1)
            xx(1:3) = - minx(1:3) - xoffbar(1:3)
            if (present(offset)) then
                xx(1:3) = xx(1:3) + a0(1:3)*offset
            else
                xx(1:3) = xx(1:3) + a0(1:3)/8
            end if
        !   now have offset required to line up atoms on 0,0,0 boundary

            if (op) then
                write (*,fmt='(a12,3i6)')    "  surf # = ",nsx,nsy,nsz
                write (*,fmt='(a12,3f16.8)') "  offset = ",xoffbar(1:3)
                write (*,fmt='(a12,3f16.8)') "  adjust = ",xx(1:3)
            end if

            do ii = 1,this%nAtoms
                this%dat(1:3,ii) = this%dat(1:3,ii) + xx(1:3)
            end do

            if (op) then
                minx(1:3) = huge(1.0)
                maxx(1:3) = -huge(1.0)
                do ii = 1,this%nAtoms
                    xx(1:3) = this%dat(1:3,ii)
                    minx(1:3) = min(xx(1:3),minx(1:3))
                    maxx(1:3) = max(xx(1:3),maxx(1:3))
                end do
                print *,"Lib_XYZFiles::findOffset info - after"
                write (*,fmt='(a12,3f16.8)') "  minx   = ",minx(1:3)
                write (*,fmt='(a12,3f16.8)') "  maxx   = ",maxx(1:3)
            end if


            return
        end subroutine adjustOffset1
    !---

        subroutine convertCSV( this,verbose )
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      assume the data is in the format
    !*          dummyline
    !*          x1,y1,z1,name1
    !*          x2,y2,z2,name2
    !*          ...
    !*      read and convert format
            type(XYZFile),intent(inout)                 ::      this
            logical,intent(in),optional                 ::      verbose
            integer                                     ::      uu
            integer                                     ::      ii,jj,nn,ioerr,nAtomNames
            character(len=XYZFILE_ATOMNAMELENGTH)       ::      atom
            logical                                     ::      ok
            character(len=XYZFILE_ATOMNAMELENGTH),dimension(1000)       ::      atom_names

            logical                                     ::      op
            real(kind=real64),dimension(3)              ::      xx



            op = .false. ; if (present(verbose)) op = verbose

            uu = findGoodUnit()

            nAtomNames = 0

            open(unit=uu,file=trim(this%filename),action="read")

                read(uu,fmt='(a)') atom         !   dummyline
                nn = 0
                do
                    read(uu,fmt=*,iostat=ii) xx,atom
                    if (ii == 0) then
                        nn = nn + 1
                        ok = .false.
                        do jj = 1,nAtomNames
                            if (trim(atom) == trim(atom_names(jj))) then
                                ok = .true.
                                exit
                            end if
                        end do
                        if (.not. ok) then
                            nAtomNames = nAtomNames + 1
                            atom_names(nAtomNames) = trim(atom)
                        end if
                    else
                        exit
                    end if
                end do

                if (op) print *,"Lib_XYZFiles::convertCSV info - read ",nn," atom data lines"

            close(unit=uu)

            call ensureAdequateStorage(this,nn,nAtomNames,nColumns=3,nHeaderLines=1)

            this%nAtomNames = nAtomNames
            this%atomName(1:this%nAtomNames) = atom_names(1:this%nAtomNames)

            this%header(1) = "# converted from .csv format data"
            this%column_description = "atom    position_x      position_y      position_z"

            open(unit=uu,file=trim(this%filename),action="read")

                read(uu,fmt='(a)') atom         !   dummyline
                do ii = 1,nn
                    read(uu,fmt=*) this%dat(1:3,ii),atom

                    this%atom(ii) = getAtomTypeFromName( this,atom )

                    if (op) then
                        if (ii<4) then
                            write(unit=*,fmt='(a8,1000f16.6)',iostat=ioerr) atom,this%dat(1:this%nColumns,ii)
                        else if ( (ii==4).and.(nn>5) ) then
                            write(unit=*,fmt='(a)',iostat=ioerr) "..."
                        else if ( (ii>nn-3).and.(nn>5) ) then
                            write(unit=*,fmt='(a8,1000f16.6)',iostat=ioerr) atom,this%dat(1:this%nColumns,ii)
                        end if
                    end if

                end do

            close(unit=uu)


            return
        end subroutine convertCSV

    !---

        pure subroutine getMinMax( this,xmin,xmax )
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      returns the minimum and maximum data values of the atoms
            type(XYZFile),intent(in)                    ::      this
            real(kind=real64),dimension(:),intent(out)  ::      xmin,xmax
            integer     ::      nc,ii
            nc = min(this%nColumns,size(xmin))

            if (this%nAtoms>0) then
                do ii = 1,nc
                    xmin(ii) = minval( this%dat(ii,1:this%natoms) )
                    xmax(ii) = maxval( this%dat(ii,1:this%natoms) )
                end do
            else
                xmin(1:nc) = 0.0d0
                xmax(1:nc) = 0.0d0
            end if
            return
        end subroutine getMinMax

    !---

        subroutine randomiseAtomTypes( this )
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      randomise the atom types, but not their positions
            type(XYZFile),intent(inout)     ::      this

            integer,dimension(this%nAtomNames)          ::      nType
            integer         ::      ii,jj,kk,mm,nn
            integer,dimension(this%nAtoms)              ::      indx
            real(kind=real64)       ::      xx

        !---    first count the atoms
            nType = 0       !   count of atoms of each type
            do ii = 1,this%nAtoms
                indx(ii) = ii
                jj = this%atom(ii)
                nType(jj) = nType(jj) + 1
            end do


        !---    fill up the atom occupations one by one
            nn = this%nAtoms        !   number of slots left to fill
            do jj = 1,this%nAtomNames
                do kk = 1,nType(jj)
                    call random_number(xx)
                    mm = int(xx*nn) + 1             !   number from 1 to number of unassigned atoms
                    this%atom( indx(mm) ) = jj      !   set the unassigned atom to type jj
                    indx(mm) = indx(nn)             !   fill in the hole in the unassigned atom list
                    nn = nn - 1                     !   reduce the number of unassigned atoms.
                end do
            end do

            return
        end subroutine randomiseAtomTypes



        subroutine setLammpsFormat(this,lammps)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      set output mode to lammps format
            type(XYZFile),intent(inout)     ::      this
            logical,intent(in)              ::      lammps
            this%lammpsFormat = lammps
            return
        end subroutine setLammpsFormat


        pure logical function isLammpsFormat(this)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      set output mode to lammps format
            type(XYZFile),intent(in)     ::      this
            isLammpsFormat = this%lammpsFormat 
            return
        end function isLammpsFormat
        
        
        
        



    !---

        subroutine voxelisedCount0( this, atom, a,sig, v,vtot )
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      construct voxelised data v(:,:,:)
    !*      for the concentration of atoms of type "atom"
    !*      where each voxel has length a
    !*      and each atom is smeared out with a gaussian blur sig
    !*      also return ( optionally ) the total count of all atoms
    !*      voxels are 0:a , a:2a, 2a:3a ...
    !*      so that nodes are a/2 , 3a/2, 5a/2 ...
    !*      set atom="" to count all atoms
    !*      set sig=0 for no smearing
            type(XYZFile),intent(in)        ::      this
            character(len=*),intent(in)     ::      atom
            real(kind=real64),intent(in)    ::      a,sig
            real(kind=real64),dimension(0:,0:,0:),intent(out)           ::      v
            real(kind=real64),dimension(0:,0:,0:),intent(out),optional  ::      vtot

            real(kind=real64),dimension(3)      ::      xx,yy
            integer                             ::      ii,ix,iy,iz, jx,jy,jz ,nx,ny,nz , nd
            integer                             ::      ixm,ixp,iym,iyp,izm,izp
            integer                             ::      tt
            real(kind=real64)                   ::      ia,i2s2,ff
            real(kind=real64),dimension(:,:,:),allocatable  ::      gg

            integer                         ::      nSig


            tt = getAtomTypeFromName( this,atom )
            nx = size(v,dim=1)          !   note v range expected (0:nx-1)
            ny = size(v,dim=2)
            nz = size(v,dim=3)
            ia = 1/a
            v  = 0.0d0
            if (present(vtot)) vtot = 0.0d0

            nSig = 3
            if ( (sig>0).and.(sig<0.5) ) nSig = 4

            if (sig*nSig*ia<1) then

            !   no gaussian smearing.
                do ii = 1,this%nAtoms
                    xx(1:3) = this%dat(1:3,ii)*ia
                    ix = floor(xx(1)) ; iy = floor(xx(2)) ; iz = floor(xx(3))
                    if ( ix*(nx-1-ix) < 0 ) cycle
                    if ( iy*(ny-1-iy) < 0 ) cycle
                    if ( iz*(nz-1-iz) < 0 ) cycle
                    if (tt*(this%atom(ii)-tt) == 0) v(ix,iy,iz) = v(ix,iy,iz) + 1
                    if (present(vtot)) vtot(ix,iy,iz) = vtot(ix,iy,iz) + 1
                end do

            else

            !   gaussian smearing
                nd = ceiling( nSig*sig*ia )          !   max spread of voxels to cover
                i2s2 = 1/(2*sig*sig)                !   1/(2 sigma^2)
                allocate(gg(-nd:nd,-nd:nd,-nd:nd))
                print *,"Lib_XYZFiles::voxelisedCount0 info - gaussian kernel range ",-nd,":",nd


                do ii = 1,this%nAtoms

                !---    where to put this gaussian kernel?
                    xx(1:3) = this%dat(1:3,ii)*ia - 0.50d0                              !

                    ix = nint(xx(1)) ; iy = nint(xx(2)) ; iz = nint(xx(3))      !   central voxel to consider - will get gg(0,0,0)
                    !   boundary issues
                    ixm = max(0,ix-nd) ; ixp = min( nx-1,ix+nd )
                    iym = max(0,iy-nd) ; iyp = min( ny-1,iy+nd )
                    izm = max(0,iz-nd) ; izp = min( nz-1,iz+nd )

                    if ( ixp*(nx-1-ixm) < 0 ) cycle                     !   if ix<1-nd or ix>nx-1+nd out of range
                    if ( iyp*(ny-1-iym) < 0 ) cycle
                    if ( izp*(nz-1-izm) < 0 ) cycle

                !---    construct gaussian kernel for this atom
                    xx(1:3) = xx(1:3) - (/ ix,iy,iz /)      !   offset to nearest node
                    gg = 0.0d0
                    do jz = -nd,nd
                        yy(3) = xx(3) - jz
                        do jy = -nd,nd
                            yy(2) = xx(2) - jy
                            do jx = -nd,nd
                                yy(1) = xx(1) - jx
                                ff = yy(1)*yy(1) + yy(2)*yy(2) + yy(3)*yy(3)
                                ff = ff*i2s2
!                                if (ii==1) print *,jx,jy,jz,yy,ff,exp(-ff)
                                if (ff <= nSig*nSig*0.51d0) gg( jx,jy,jz ) = exp( -ff )
                            end do
                        end do
                    end do
!                     ff = 1/sum(gg)
!                     gg(-nd:nd,-nd:nd,-nd:nd) = gg(-nd:nd,-nd:nd,-nd:nd)*ff
                    ff = sum( gg(ixm-ix:ixp-ix,iym-iy:iyp-iy,izm-iz:izp-iz) )
                    ff = 1/ff
                    gg(-nd:nd,-nd:nd,-nd:nd) = gg(-nd:nd,-nd:nd,-nd:nd)*ff

                    if (tt*(this%atom(ii)-tt) == 0) v(ixm:ixp,iym:iyp,izm:izp) = v(ixm:ixp,iym:iyp,izm:izp) + gg(ixm-ix:ixp-ix,iym-iy:iyp-iy,izm-iz:izp-iz)
                    if (present(vtot)) vtot(ixm:ixp,iym:iyp,izm:izp) = vtot(ixm:ixp,iym:iyp,izm:izp) + gg(ixm-ix:ixp-ix,iym-iy:iyp-iy,izm-iz:izp-iz)

                end do

                deallocate(gg)
            end if
            print *,"Lib_XYZFiles::voxelisedCount0 info - count atom type "//trim(atom)//" ",count(this%atom(1:this%nAtoms)==tt)," -> ",sum(v)

            return
        end subroutine voxelisedCount0


        subroutine voxelisedCount1( this, a,sig, v,uniform )
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      construct voxelised data v(1:ntypes,0:nx-1,0:ny-1,0:nz-1)
    !*      for the concentration of atoms of all types
    !*      where each voxel has length a
    !*      and each atom is smeared out with a gaussian blur sig
    !*      voxels are 0:a , a:2a, 2a:3a ...
    !*      so that nodes are a/2 , 3a/2, 5a/2 ...
    !*      set sig=0 for no smearing
            type(XYZFile),intent(in)        ::      this
            real(kind=real64),intent(in)    ::      a,sig
            real(kind=real64),dimension(1:,0:,0:,0:),intent(out)           ::      v
            logical,intent(in),optional     ::      uniform

            real(kind=real64),dimension(3)      ::      xx,yy
            integer                             ::      ii,ix,iy,iz, jx,jy,jz ,nx,ny,nz , nd , tt
            integer                             ::      ixm,ixp,iym,iyp,izm,izp
            real(kind=real64)                   ::      ia,i2s2,ff
            real(kind=real64),dimension(:,:,:),allocatable  ::      gg

            real(kind=real64),dimension(this%nAtomNames)    ::      cc

            logical     ::      uni
            integer                         ::      nSig
            uni = .false. ; if (present(uniform)) uni = uniform


            if (uni) then
            !   find concentration of atoms
                do tt = 1,this%nAtomNames
                    cc(tt) = count(this%atom(1:this%nAtoms)==tt)
                end do
                ff = sum(cc) ; ff = 1/ff
                cc(:) = cc(:)*ff
            end if

            nx = size(v,dim=2)          !   note v range expected (0:nx-1)
            ny = size(v,dim=3)
            nz = size(v,dim=4)
            ia = 1/a
            v  = 0.0d0

            nSig = 3
            if ( (sig>0).and.(sig*ia<0.5) ) nSig = 4
            if (sig*nSig*ia<1) then

            !   no gaussian smearing.
                do ii = 1,this%nAtoms
                    xx(1:3) = this%dat(1:3,ii)*ia
                    ix = floor(xx(1)) ; iy = floor(xx(2)) ; iz = floor(xx(3))
                    if ( ix*(nx-1-ix) < 0 ) cycle
                    if ( iy*(ny-1-iy) < 0 ) cycle
                    if ( iz*(nz-1-iz) < 0 ) cycle

                    if (uni) then
                        do tt = 1,this%nAtomNames
                            v(tt,ix,iy,iz) = v(tt,ix,iy,iz) + cc(tt)
                        end do
                    else
                        tt = this%atom(ii)
                        v(tt,ix,iy,iz) = v(tt,ix,iy,iz) + 1
                    end if
                end do

            else

            !   gaussian smearing
                nd = ceiling( nSig*sig*ia )          !   max spread of voxels to cover
                i2s2 = 1/(2*sig*sig)                !   1/(2 sigma^2)
                allocate(gg(-nd:nd,-nd:nd,-nd:nd))
                print *,"Lib_XYZFiles::voxelisedCount1 info - gaussian kernel range ",-nd,":",nd


                do ii = 1,this%nAtoms



                !---    where to put this gaussian kernel?
                    xx(1:3) = this%dat(1:3,ii)*ia - 0.50d0                              !

                    ix = nint(xx(1)) ; iy = nint(xx(2)) ; iz = nint(xx(3))      !   central voxel to consider - will get gg(0,0,0)
                    !   boundary issues
                    ixm = max(0,ix-nd) ; ixp = min( nx-1,ix+nd )
                    iym = max(0,iy-nd) ; iyp = min( ny-1,iy+nd )
                    izm = max(0,iz-nd) ; izp = min( nz-1,iz+nd )

                    if ( ixp*(nx-1-ixm) < 0 ) cycle                     !   if ix<1-nd or ix>nx-1+nd out of range
                    if ( iyp*(ny-1-iym) < 0 ) cycle
                    if ( izp*(nz-1-izm) < 0 ) cycle

                !---    construct gaussian kernel for this atom
                    xx(1:3) = xx(1:3) - (/ ix,iy,iz /)      !   offset to nearest node
                    gg = 0.0d0
                    do jz = -nd,nd
                        yy(3) = xx(3) - jz
                        do jy = -nd,nd
                            yy(2) = xx(2) - jy
                            do jx = -nd,nd
                                yy(1) = xx(1) - jx
                                ff = yy(1)*yy(1) + yy(2)*yy(2) + yy(3)*yy(3)
                                ff = ff*i2s2
!                                if (ii==1) print *,jx,jy,jz,yy,ff,exp(-ff)
                                if (ff <= nSig*nSig*0.51d0 ) gg( jx,jy,jz ) = exp( -ff )
                            end do
                        end do
                    end do
                    ff = sum( gg(ixm-ix:ixp-ix,iym-iy:iyp-iy,izm-iz:izp-iz) )
                    ff = 1/ff
                    gg(-nd:nd,-nd:nd,-nd:nd) = gg(-nd:nd,-nd:nd,-nd:nd)*ff

                    if (uni) then
                        do tt = 1,this%nAtomNames
                            v(tt,ixm:ixp,iym:iyp,izm:izp) = v(tt,ixm:ixp,iym:iyp,izm:izp) + cc(tt)*gg(ixm-ix:ixp-ix,iym-iy:iyp-iy,izm-iz:izp-iz)
                        end do
                    else
                        tt = this%atom(ii)
                        v(tt,ixm:ixp,iym:iyp,izm:izp) = v(tt,ixm:ixp,iym:iyp,izm:izp) + gg(ixm-ix:ixp-ix,iym-iy:iyp-iy,izm-iz:izp-iz)
                    end if

                end do

                deallocate(gg)
            end if
            do tt = 1,this%nAtomNames
                print *,"count "//this%atomName(tt)//" ",count(this%atom(1:this%nAtoms)==tt)," -> ",sum(v(tt,:,:,:))
            end do
            return
        end subroutine voxelisedCount1

        subroutine voxelisedCount2( this,col, a,sig, v )
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      construct voxelised data v(1:ntypes,0:nx-1,0:ny-1,0:nz-1)
    !*      for the set of columns col(1),col(2),...,col(nc)
    !*      where each voxel has length a
    !*      and each atom is smeared out with a gaussian blur sig
    !*      returns total count in column nc+1 if possible
    !*      voxels are 0:a , a:2a, 2a:3a ...
    !*      so that nodes are a/2 , 3a/2, 5a/2 ...
    !*      set sig=0 for no smearing
            type(XYZFile),intent(in)        ::      this
            integer,dimension(:),intent(in) ::      col
            real(kind=real64),intent(in)    ::      a,sig
            real(kind=real64),dimension(1:,0:,0:,0:),intent(out)           ::      v


            real(kind=real64),dimension(3)      ::      xx,yy
            integer                             ::      nc
            integer                             ::      ii,ix,iy,iz, jx,jy,jz ,nx,ny,nz , nd , tt
            integer                             ::      ixm,ixp,iym,iyp,izm,izp
            real(kind=real64)                   ::      ia,i2s2,ff
            real(kind=real64),dimension(:,:,:),allocatable  ::      gg



            real(kind=real64),dimension(0:size(v,dim=2)-1,0:size(v,dim=2)-1,0:size(v,dim=2)-1)           ::      vcount
            integer                         ::      nSig


            nc = size(col)
            do tt = 1,nc
                if (col(tt) > this%nColumns) stop "Lib_XYZFiles::voxelisedCount2 error - attempt to count column out of data range"
            end do


            nx = size(v,dim=2)          !   note v range expected (0:nx-1)
            ny = size(v,dim=3)
            nz = size(v,dim=4)
            ia = 1/a
            v  = 0.0d0
            vcount = 0

            nSig = 3
            if ( (sig>0).and.(sig*ia<0.5) ) nSig = 4
            if (sig*nSig*ia<1) then

            !   no gaussian smearing.
                do ii = 1,this%nAtoms
                    xx(1:3) = this%dat(1:3,ii)*ia
                    ix = floor(xx(1)) ; iy = floor(xx(2)) ; iz = floor(xx(3))
                    if ( ix*(nx-1-ix) < 0 ) cycle
                    if ( iy*(ny-1-iy) < 0 ) cycle
                    if ( iz*(nz-1-iz) < 0 ) cycle

                    do tt = 1,nc
                        v(tt,ix,iy,iz) = v(tt,ix,iy,iz) + this%dat(col(tt),ii)
                    end do
                    vcount(ix,iy,iz) = vcount(ix,iy,iz) + 1


                end do

            else

            !   gaussian smearing
                nd = ceiling( nSig*sig*ia )          !   max spread of voxels to cover
                i2s2 = 1/(2*sig*sig)                !   1/(2 sigma^2)
                allocate(gg(-nd:nd,-nd:nd,-nd:nd))
                print *,"Lib_XYZFiles::voxelisedCount1 info - gaussian kernel range ",-nd,":",nd


                do ii = 1,this%nAtoms



                !---    where to put this gaussian kernel?
                    xx(1:3) = this%dat(1:3,ii)*ia - 0.50d0                              !

                    ix = nint(xx(1)) ; iy = nint(xx(2)) ; iz = nint(xx(3))      !   central voxel to consider - will get gg(0,0,0)
                    !   boundary issues
                    ixm = max(0,ix-nd) ; ixp = min( nx-1,ix+nd )
                    iym = max(0,iy-nd) ; iyp = min( ny-1,iy+nd )
                    izm = max(0,iz-nd) ; izp = min( nz-1,iz+nd )

                    if ( ixp*(nx-1-ixm) < 0 ) cycle                     !   if ix<1-nd or ix>nx-1+nd out of range
                    if ( iyp*(ny-1-iym) < 0 ) cycle
                    if ( izp*(nz-1-izm) < 0 ) cycle

                !---    construct gaussian kernel for this atom
                    xx(1:3) = xx(1:3) - (/ ix,iy,iz /)      !   offset to nearest node
                    gg = 0.0d0
                    do jz = -nd,nd
                        yy(3) = xx(3) - jz
                        do jy = -nd,nd
                            yy(2) = xx(2) - jy
                            do jx = -nd,nd
                                yy(1) = xx(1) - jx
                                ff = yy(1)*yy(1) + yy(2)*yy(2) + yy(3)*yy(3)
                                ff = ff*i2s2
!                                if (ii==1) print *,jx,jy,jz,yy,ff,exp(-ff)
                                if (ff <= nSig*nSig*0.51d0) gg( jx,jy,jz ) = exp( -ff )
                            end do
                        end do
                    end do
!                     ff = 1/sum(gg)
!                     gg(-nd:nd,-nd:nd,-nd:nd) = gg(-nd:nd,-nd:nd,-nd:nd)*ff
                    ff = sum( gg(ixm-ix:ixp-ix,iym-iy:iyp-iy,izm-iz:izp-iz) )
                    ff = 1/ff
                    gg(-nd:nd,-nd:nd,-nd:nd) = gg(-nd:nd,-nd:nd,-nd:nd)*ff



                    do tt = 1,nc
                        v(tt,ixm:ixp,iym:iyp,izm:izp) = v(tt,ixm:ixp,iym:iyp,izm:izp) + gg(ixm-ix:ixp-ix,iym-iy:iyp-iy,izm-iz:izp-iz)*this%dat(col(tt),ii)
                    end do
                    vcount(ixm:ixp,iym:iyp,izm:izp) = vcount(ixm:ixp,iym:iyp,izm:izp) + gg(ixm-ix:ixp-ix,iym-iy:iyp-iy,izm-iz:izp-iz)

                end do

                deallocate(gg)
            end if

            do iz = 0,nz-1
                do iy = 0,ny-1
                    do ix = 0,nx-1
                        ff = vcount(ix,iy,iz)
                        if (ff>0) v(1:nc,ix,iy,iz) = v(1:nc,ix,iy,iz) / ff
                    end do
                end do
            end do

            if (size(v,dim=1)>nc) v(nc+1,0:,0:,0:) = vcount(0:,0:,0:)


            return
        end subroutine voxelisedCount2

        subroutine voxelisedCount1a( this, a,sig, v,uniform )
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      construct voxelised data v(1:ntypes,0:nx-1,0:ny-1,0:nz-1)
    !*      for the concentration of atoms of all types
    !*      where each voxel has length a
    !*      and each atom is smeared out with a gaussian blur sig
    !*      voxels are 0:a , a:2a, 2a:3a ...
    !*      so that nodes are a/2 , 3a/2, 5a/2 ...
    !*      set sig=0 for no smearing
            type(XYZFile),intent(in)        ::      this
            real(kind=real64),dimension(3),intent(in)    ::      a
            real(kind=real64),intent(in)    ::      sig
            real(kind=real64),dimension(1:,0:,0:,0:),intent(out)           ::      v
            logical,intent(in),optional     ::      uniform

            real(kind=real64),dimension(3)      ::      xx,yy,ia
            integer                             ::      ii,ix,iy,iz, jx,jy,jz ,nx,ny,nz , nd , tt
            integer                             ::      ixm,ixp,iym,iyp,izm,izp
            real(kind=real64)                   ::      i2s2,ff,iamax
            real(kind=real64),dimension(:,:,:),allocatable  ::      gg

            real(kind=real64),dimension(this%nAtomNames)    ::      cc

            logical     ::      rnd
            integer                         ::      nSig
            rnd = .false. ; if (present(uniform)) rnd = uniform


            if (rnd) then
            !   find concentration of atoms
                do tt = 1,this%nAtomNames
                    cc(tt) = count(this%atom(1:this%nAtoms)==tt)
                end do
                ff = sum(cc) ; ff = 1/ff
                cc(:) = cc(:)*ff
            end if

            nx = size(v,dim=2)          !   note v range expected (0:nx-1)
            ny = size(v,dim=3)
            nz = size(v,dim=4)
            ia(1:3) = 1/a(1:3)
            iamax = maxval(ia)
            v  = 0.0d0

            nSig = 3
            if ( (sig>0).and.(sig*iamax<0.5) ) nSig = 4
            if (sig*nSig*iamax<1) then

            !   no gaussian smearing.
                do ii = 1,this%nAtoms
                    xx(1:3) = this%dat(1:3,ii)*ia(1:3)
                    ix = floor(xx(1)) ; iy = floor(xx(2)) ; iz = floor(xx(3))
                    if ( ix*(nx-1-ix) < 0 ) cycle
                    if ( iy*(ny-1-iy) < 0 ) cycle
                    if ( iz*(nz-1-iz) < 0 ) cycle

                    if (rnd) then
                        do tt = 1,this%nAtomNames
                            v(tt,ix,iy,iz) = v(tt,ix,iy,iz) + cc(tt)
                        end do
                    else
                        tt = this%atom(ii)
                        v(tt,ix,iy,iz) = v(tt,ix,iy,iz) + 1
                    end if
                end do

            else

            !   gaussian smearing
                nd = ceiling( nSig*sig*iamax )          !   max spread of voxels to cover
                i2s2 = 1/(2*sig*sig)                !   1/(2 sigma^2)
                allocate(gg(-nd:nd,-nd:nd,-nd:nd))
                print *,"Lib_XYZFiles::voxelisedCount1 info - gaussian kernel range ",-nd,":",nd


                do ii = 1,this%nAtoms



                !---    where to put this gaussian kernel?
                    xx(1:3) = this%dat(1:3,ii)*ia(1:3) - 0.50d0                              !

                    ix = nint(xx(1)) ; iy = nint(xx(2)) ; iz = nint(xx(3))      !   central voxel to consider - will get gg(0,0,0)
                    !   boundary issues
                    ixm = max(0,ix-nd) ; ixp = min( nx-1,ix+nd )
                    iym = max(0,iy-nd) ; iyp = min( ny-1,iy+nd )
                    izm = max(0,iz-nd) ; izp = min( nz-1,iz+nd )

                    if ( ixp*(nx-1-ixm) < 0 ) cycle                     !   if ix<1-nd or ix>nx-1+nd out of range
                    if ( iyp*(ny-1-iym) < 0 ) cycle
                    if ( izp*(nz-1-izm) < 0 ) cycle

                !---    construct gaussian kernel for this atom
                    xx(1:3) = xx(1:3) - (/ ix,iy,iz /)      !   offset to nearest node
                    gg = 0.0d0
                    do jz = -nd,nd
                        yy(3) = xx(3) - jz
                        do jy = -nd,nd
                            yy(2) = xx(2) - jy
                            do jx = -nd,nd
                                yy(1) = xx(1) - jx
                                ff = yy(1)*yy(1) + yy(2)*yy(2) + yy(3)*yy(3)
                                ff = ff*i2s2
!                                if (ii==1) print *,jx,jy,jz,yy,ff,exp(-ff)
                                if (ff <= nSig*nSig*0.51d0) gg( jx,jy,jz ) = exp( -ff )
                            end do
                        end do
                    end do
!                     ff = 1/sum(gg)
!                     gg(-nd:nd,-nd:nd,-nd:nd) = gg(-nd:nd,-nd:nd,-nd:nd)*ff
                    ff = sum( gg(ixm-ix:ixp-ix,iym-iy:iyp-iy,izm-iz:izp-iz) )
                    ff = 1/ff
                    gg(-nd:nd,-nd:nd,-nd:nd) = gg(-nd:nd,-nd:nd,-nd:nd)*ff

                    if (rnd) then
                        do tt = 1,this%nAtomNames
                            v(tt,ixm:ixp,iym:iyp,izm:izp) = v(tt,ixm:ixp,iym:iyp,izm:izp) + cc(tt)*gg(ixm-ix:ixp-ix,iym-iy:iyp-iy,izm-iz:izp-iz)
                        end do
                    else
                        tt = this%atom(ii)
                        v(tt,ixm:ixp,iym:iyp,izm:izp) = v(tt,ixm:ixp,iym:iyp,izm:izp) + gg(ixm-ix:ixp-ix,iym-iy:iyp-iy,izm-iz:izp-iz)
                    end if

                end do

                deallocate(gg)
            end if
            do tt = 1,this%nAtomNames
                print *,"count "//this%atomName(tt)//" ",count(this%atom(1:this%nAtoms)==tt)," -> ",sum(v(tt,:,:,:))
            end do
            return
        end subroutine voxelisedCount1a

        subroutine voxelisedCount2a( this,col, a,sig, v )
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      construct voxelised data v(1:ntypes,0:nx-1,0:ny-1,0:nz-1)
    !*      for the set of columns col(1),col(2),...,col(nc)
    !*      where each voxel has length a
    !*      and each atom is smeared out with a gaussian blur sig
    !*      returns total count in column nc+1 if possible
    !*      voxels are 0:a , a:2a, 2a:3a ...
    !*      so that nodes are a/2 , 3a/2, 5a/2 ...
    !*      set sig=0 for no smearing
            type(XYZFile),intent(in)        ::      this
            integer,dimension(:),intent(in) ::      col
            real(kind=real64),dimension(3),intent(in)    ::      a
            real(kind=real64),intent(in)    ::      sig
            real(kind=real64),dimension(1:,0:,0:,0:),intent(out)           ::      v


            real(kind=real64),dimension(3)      ::      xx,yy,ia
            integer                             ::      nc
            integer                             ::      ii,ix,iy,iz, jx,jy,jz ,nx,ny,nz , nd , tt
            integer                             ::      ixm,ixp,iym,iyp,izm,izp
            real(kind=real64)                   ::      i2s2,ff,iamax

            real(kind=real64),dimension(:,:,:),allocatable  ::      gg

            integer                             ::      nSig

            real(kind=real64),dimension(0:size(v,dim=2)-1,0:size(v,dim=2)-1,0:size(v,dim=2)-1)           ::      vcount



            nc = size(col)
            do tt = 1,nc
                if (col(tt) > this%nColumns) stop "Lib_XYZFiles::voxelisedCount2 error - attempt to count column out of data range"
            end do


            nx = size(v,dim=2)          !   note v range expected (0:nx-1)
            ny = size(v,dim=3)
            nz = size(v,dim=4)
            ia(1:3) = 1/a(1:3)
            iamax = maxval(ia)
            v  = 0.0d0
            vcount = 0

            nSig = 3
            if ( (sig>0).and.(sig*iamax<0.5) ) nSig = 4
            if (sig*nSig*iamax<1) then

            !   no gaussian smearing.
                do ii = 1,this%nAtoms
                    xx(1:3) = this%dat(1:3,ii)*ia(1:3)
                    ix = floor(xx(1)) ; iy = floor(xx(2)) ; iz = floor(xx(3))
                    if ( ix*(nx-1-ix) < 0 ) cycle
                    if ( iy*(ny-1-iy) < 0 ) cycle
                    if ( iz*(nz-1-iz) < 0 ) cycle

                    do tt = 1,nc
                        v(tt,ix,iy,iz) = v(tt,ix,iy,iz) + this%dat(col(tt),ii)
                    end do
                    vcount(ix,iy,iz) = vcount(ix,iy,iz) + 1


                end do

            else

            !   gaussian smearing
                nd = ceiling( nSig*sig*iamax )          !   max spread of voxels to cover
                i2s2 = 1/(2*sig*sig)                !   1/(2 sigma^2)
                allocate(gg(-nd:nd,-nd:nd,-nd:nd))
                print *,"Lib_XYZFiles::voxelisedCount1 info - gaussian kernel range ",-nd,":",nd


                do ii = 1,this%nAtoms



                !---    where to put this gaussian kernel?
                    xx(1:3) = this%dat(1:3,ii)*ia(1:3) - 0.50d0                              !

                    ix = nint(xx(1)) ; iy = nint(xx(2)) ; iz = nint(xx(3))      !   central voxel to consider - will get gg(0,0,0)
                    !   boundary issues
                    ixm = max(0,ix-nd) ; ixp = min( nx-1,ix+nd )
                    iym = max(0,iy-nd) ; iyp = min( ny-1,iy+nd )
                    izm = max(0,iz-nd) ; izp = min( nz-1,iz+nd )

                    if ( ixp*(nx-1-ixm) < 0 ) cycle                     !   if ix<1-nd or ix>nx-1+nd out of range
                    if ( iyp*(ny-1-iym) < 0 ) cycle
                    if ( izp*(nz-1-izm) < 0 ) cycle

                !---    construct gaussian kernel for this atom
                    xx(1:3) = xx(1:3) - (/ ix,iy,iz /)      !   offset to nearest node
                    gg = 0.0d0
                    do jz = -nd,nd
                        yy(3) = xx(3) - jz
                        do jy = -nd,nd
                            yy(2) = xx(2) - jy
                            do jx = -nd,nd
                                yy(1) = xx(1) - jx
                                ff = yy(1)*yy(1) + yy(2)*yy(2) + yy(3)*yy(3)
                                ff = ff*i2s2
!                                if (ii==1) print *,jx,jy,jz,yy,ff,exp(-ff)
                                if (ff <= nSig*nSig*0.51d0) gg( jx,jy,jz ) = exp( -ff )
                            end do
                        end do
                    end do
!                     ff = 1/sum(gg)
!                     gg(-nd:nd,-nd:nd,-nd:nd) = gg(-nd:nd,-nd:nd,-nd:nd)*ff
                    ff = sum( gg(ixm-ix:ixp-ix,iym-iy:iyp-iy,izm-iz:izp-iz) )
                    ff = 1/ff
                    gg(-nd:nd,-nd:nd,-nd:nd) = gg(-nd:nd,-nd:nd,-nd:nd)*ff



                    do tt = 1,nc
                        v(tt,ixm:ixp,iym:iyp,izm:izp) = v(tt,ixm:ixp,iym:iyp,izm:izp) + gg(ixm-ix:ixp-ix,iym-iy:iyp-iy,izm-iz:izp-iz)*this%dat(col(tt),ii)
                    end do
                    vcount(ixm:ixp,iym:iyp,izm:izp) = vcount(ixm:ixp,iym:iyp,izm:izp) + gg(ixm-ix:ixp-ix,iym-iy:iyp-iy,izm-iz:izp-iz)

                end do

                deallocate(gg)
            end if

            do iz = 0,nz-1
                do iy = 0,ny-1
                    do ix = 0,nx-1
                        ff = vcount(ix,iy,iz)
                        if (ff>0) v(1:nc,ix,iy,iz) = v(1:nc,ix,iy,iz) / ff
                    end do
                end do
            end do

            if (size(v,dim=1)>nc) v(nc+1,0:,0:,0:) = vcount(0:,0:,0:)


            return
        end subroutine voxelisedCount2a

    !---
    
        function idByMass(m) result(name)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      attempt to id an atom by mass. 
            real(kind=real64),intent(in)            ::      m
            character(len=XYZFILE_ATOMNAMELENGTH)   ::      name 
            integer         ::      ii
            name = "unset"
            do ii = 1,XYZFILE_NATOMSIDBYMASS
                if ( abs(m - XYZFILE_ATOMSIDBYMASS_MASS(ii))<0.50d0 ) then   
                    !   close enough.
                    name = XYZFILE_ATOMSIDBYMASS_NAME(ii)
                    return
                end if
            end do
            return
        end function idByMass
            
        function getMass0(name) result(m)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      attempt to find mass from id 
            character(len=XYZFILE_ATOMNAMELENGTH),intent(in)    ::      name 
            real(kind=real64)                                   ::      m
            integer         ::      ii
            m = 1.0d0       !   don't know.
            do ii = 1,XYZFILE_NATOMSIDBYMASS
                if (trim(name)==trim(XYZFILE_ATOMSIDBYMASS_NAME(ii))) then
                    m = XYZFILE_ATOMSIDBYMASS_MASS(ii)
                    return
                end if
            end do
            return
        end function getMass0
        
        
!-------

        function findGoodUnit() result(u)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      find a good (unopened) unit to read/write
            integer             ::      u
            logical             ::      lod
            u = 400
            do
                inquire (unit = u,opened = lod)
                if (.not. lod) exit
                u = u + 1
            end do
            return
        end function findGoodUnit





    end module Lib_XYZFiles