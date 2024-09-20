
    module Lib_CommandLineArguments
!---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
!*      find key-pair values from command-line arguments 
!*      in the form "-key value[,val2,val3]"
!*
!*      usage:
!*          type(CommandLineArguments)      ::      cla
!*
!*          cla = CommandLineArguments_ctor(10)         !   allocates space for up to 10 arguments
!*
!*          !   optionally add a number of subcategories to organise options
!*          call setCategories(cla,(/ "category 1","category 2" /))
!*
!*          !   optionally request that the first few arguments are handled elsewhere, and not of form "-key value" 
!*          !   (eg might want ./a.exe <filename> [-option1 ... ] - in which case ignore the first argument with)
!*          call setIgnore( cla,1 )
!*
!*          !   optionally might want to set the program name and some description          
!*          call setProgramDescription( cla, "./a.exe <filename> [-options ...] \n does something awesome to <filename>" )
!*          call setProgramVersion( cla, "1.0 beta" )
!*
!*
!*          !   add an option, return its value
!*          call get( cla,"key",value [,n_value] ,optional,"description" [,category] )  
!*                       "key" should be used without -sign in front, ie key="m" means expect command line arg "-m"  
!*                        optional = LIB_CLA_OPTIONAL/LIB_CLA_REQUIRED
!*                        category required if set in previous step. 
!*                        n_value returns number of array arguments returned. If on input /= 0 then is required minimum number.
!*
!*          !   report the options available with 
!*          call report(cla [,u] [,o])              !   gives full output iff there is an error, or "-h" is a cla
!*
!*          !   optionally eg exit
!*          if (.not. allRequiredArgumentsSet(cla)) stop
!*          if (hasHelpArgument(cla)) stop
!*
!*          call delete(cla)
!*
!*
!*      Daniel Mason, UKAEA
!*      April 2022
!*
!*

        use iso_fortran_env
        use NBAX_StringTokenizers
        implicit none
        private
        
    !---
    
        logical,public,parameter        ::      LIB_CLA_OPTIONAL = .true.
        logical,public,parameter        ::      LIB_CLA_REQUIRED = .false.
        
        integer(kind=int64),private,parameter               ::      BADF00D = int( z'BADF00D',kind=int64 )
        real(kind=real64),public,parameter                  ::      LIB_CLA_NODEFAULT_R = transfer( (BADF00D+ishft(BADF00D,32_int64)),1.0d0 )
        integer,public,parameter                            ::      LIB_CLA_NODEFAULT_I = int(BADF00D,kind=int32)
        !character(len=7),public,parameter                   ::      LIB_CLA_NODEFAULT_C = "BADF00D"    
    
        integer,private,parameter       ::      LIB_CLA_BOOL  = 1       !   defines an argument to be of type logical
        integer,private,parameter       ::      LIB_CLA_INT   = 2
        integer,private,parameter       ::      LIB_CLA_REAL  = 3
        integer,private,parameter       ::      LIB_CLA_CHAR  = 4
        integer,private,parameter       ::      LIB_CLA_BOOLA = 5
        integer,private,parameter       ::      LIB_CLA_INTA  = 6       !   defines an argument to be of type integer array                                                   
        integer,private,parameter       ::      LIB_CLA_REALA = 7                                                       
        integer,private,parameter       ::      LIB_CLA_CHARA = 8  
                                                             
        
        
        integer,private,parameter       ::      LIB_CLA_KEYLEN = 32
        integer,private,parameter       ::      LIB_CLA_CATLEN = 128
        integer,private,parameter       ::      LIB_CLA_DESCLEN = 4096
        
    !---
    
        public      ::      commandLineArguments_ctor
        public      ::      delete
        public      ::      report
        
        public      ::      reportProgramDescription
        public      ::      setCategories
        public      ::      setIgnore
        public      ::      setProgramDescription
        public      ::      setProgramVersion
        
        public      ::      get
        public      ::      allRequiredArgumentsSet
        public      ::      hasArgument
        public      ::      hasHelpArgument
        
    !---
    
        type,public         ::      CommandLineArguments
            private           
            
            integer         ::      nArgs                   !   number of arguments registered
            integer         ::      nArgs0                  !   length of arrays = max number which can be registered                               
            
            integer         ::      ignore                  !   ignore the first few CLA - they are handled elsewhere
            integer         ::      nCategories             !   allows for sorting of arguments into categories.
        
            logical,dimension(:),pointer        ::      arg_optional
            integer,dimension(:),pointer        ::      arg_type
            integer,dimension(:),pointer        ::      arg_cat
            
            character(len=LIB_CLA_DESCLEN)      ::      program_description
            character(len=32)                   ::      program_version
            
            character(len=LIB_CLA_KEYLEN),dimension(:),pointer     ::      arg_key
            character(len=LIB_CLA_CATLEN),dimension(:),pointer     ::      cat_description
            character(len=LIB_CLA_DESCLEN),dimension(:),pointer    ::      arg_description
            
            logical,dimension(:),pointer        ::      arg_optional_set        
        end type 
        
    !---
    
        interface    commandLineArguments_ctor
            module procedure       commandLineArguments_null
            module procedure       commandLineArguments_ctor0
        end interface
    
        interface    report
            module procedure       report0
        end interface
    
        interface    delete
            module procedure       delete0
        end interface
    
        interface    get
            module procedure       get1         !   bool
            module procedure       get2         !   int
            module procedure       get3         !   real
            module procedure       get4         !   char
            module procedure       get5         !   int a
            module procedure       get6         !   real a
            module procedure       get7         !   char a
            module procedure       get8         !   bool a
        end interface
        
        
    contains
!---^^^^^^^^


        function commandLineArguments_null() result(this)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
            type(commandLineArguments)          ::      this
            
            this%nArgs = 0
            this%nArgs0 = 0
            this%nCategories = 0
            this%ignore = 0
            this%program_description = ""
            this%program_version = "unspecified"
            nullify(this%arg_optional)
            nullify(this%arg_type)
            nullify(this%arg_optional_set)
            nullify(this%arg_cat)    
            nullify(this%arg_key)    
            nullify(this%cat_description)    
            nullify(this%arg_description)
            
            return
        end function commandLineArguments_null    

        function commandLineArguments_ctor0(nArgs) result(this)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
            integer,intent(in)                  ::      nArgs
            type(commandLineArguments)          ::      this
            this = commandLineArguments_null()
            this%nArgs0 = nArgs
            allocate(this%arg_optional(this%nArgs0))
            allocate(this%arg_type(this%nArgs0))
            allocate(this%arg_optional_set(this%nArgs0))
            allocate(this%arg_cat(this%nArgs0))    
            allocate(this%arg_key(this%nArgs0))    
            allocate(this%cat_description(this%nArgs0))   
            allocate(this%arg_description(this%nArgs0))        
            this%arg_optional_set = .false.
            
            return
        end function commandLineArguments_ctor0    
        
    !---
    
        subroutine delete0(this)
    !---^^^^^^^^^^^^^^^^^^^^^^^^
            type(commandLineArguments),intent(inout)     ::      this
            if (this%nArgs0 == 0) return
            deallocate(this%arg_optional)
            deallocate(this%arg_type)
            deallocate(this%arg_optional_set)
            deallocate(this%arg_cat)    
            deallocate(this%arg_key)    
            deallocate(this%cat_description)
            deallocate(this%arg_description)            
            this = commandLineArguments_null()
            return
        end subroutine delete0
        
        
    !---
    
        function hasArgument(this,arg_key) result(is)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ 
            type(commandLineArguments),intent(in)       ::      this
            character(len=*),intent(in)                 ::      arg_key
            logical                                     ::      is
            integer         ::      kk,nn   
            character(len=256)  ::      cla         
            is = .false.
            
            kk = this%ignore 
            nn = command_argument_count()
            do kk = this%ignore+1,nn
                call get_command_argument( kk,cla )
                if (trim(cla) == "-"//trim(arg_key)) then          !   have identified a key
                    is = .true.
                    exit
                else if (trim(cla) == "-no"//trim(arg_key)) then          !   have identified a -nokey
                    is = .true.
                    exit
                end if
            end do
                        
            return
        end function hasArgument
        
        function allRequiredArgumentsSet(this) result(is)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
            type(commandLineArguments),intent(in)       ::      this
            logical                                     ::      is
            integer         ::      jj            
            is = .true.
            do jj = 1,this%nArgs
                if ( (.not. this%arg_optional(jj)).and.(.not. this%arg_optional_set(jj) ) ) is = .false.
            end do
            return
        end function allRequiredArgumentsSet
    
        function hasHelpArgument(this) result(has)      !   note: dummy argument not used.
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
            type(commandLineArguments),intent(in)   ::      this        
            logical                                 ::      has
            integer             ::      kk
            character(len=256)  ::      cla
            has = .false.
            do kk = 1,command_argument_count()
                call get_command_argument( kk,cla )
                if ( (trim(cla)=="-h").or.(trim(cla)=="-help").or.(trim(cla)=="--h").or.(trim(cla)=="--help") ) then 
                    has = .true.
                    exit
                end if
            end do
            return
        end function hasHelpArgument
        
        subroutine reportProgramDescription(this,u,o)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
            type(commandLineArguments),intent(in)       ::      this
            integer,intent(in),optional                 ::      u,o
            integer                                     ::      uu,oo
                        
            character(len=LIB_CLA_DESCLEN),dimension(:),pointer       ::      desc
            integer                                     ::      nDesc
            integer     ::      kk
            
            uu = 6 ; if (present(u)) uu = u
            oo = 0 ; if (present(o)) oo = o 
            
            if ( (len_trim(this%program_description)==0) .and. (this%nArgs0 == 0) ) then
                write (unit=uu,fmt='(a)') repeat(" ",oo)//"CommandLineArguments[not set]"
                return
            end if
            
            if ( len_trim(this%program_description)>0 ) then                
 
                call breakUpLine( this%program_description,nDesc,desc )                      
                write (unit=uu,fmt='(a)',advance="yes") trim(desc(1))
                do kk = 2,nDesc
                    write (unit=uu,fmt='(a)',advance="yes") trim(desc(kk))
                end do                
                deallocate(desc)                                            
                write (unit=uu,fmt='(a)') ""        
            end if
            
            return
        end subroutine reportProgramDescription
            
        
        subroutine breakUpLine( desc_in,nDesc_lines,desc_out )
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      break up the input line "desc_in" into several output lines
            character(len=*),intent(in)                                         ::      desc_in
            integer,intent(out)                                                 ::      nDesc_lines
            character(len=LIB_CLA_DESCLEN),dimension(:),pointer,intent(out)     ::      desc_out
            character(len=LIB_CLA_DESCLEN)      ::      desc0,desc1
            integer             ::      kk
            
            desc0 = desc_in
            nDesc_lines = 0
            do 
                kk = index( desc0,"\n" )
                if (kk == 0) then
                    nDesc_lines = nDesc_lines + 1
                    exit                                   
                else
                    nDesc_lines = nDesc_lines + 1
                    desc0 = desc0(kk+2:)
                end if
            end do     
            desc0 = desc_in
            allocate(desc_out(nDesc_lines))  
            nDesc_lines = 0    
            do 
                kk = index( desc0,"\n" )
                if (kk == 0) then
                    nDesc_lines = nDesc_lines + 1
                    desc_out(nDesc_lines) = desc0  
                    exit                                   
                else
                    nDesc_lines = nDesc_lines + 1
                    desc1 = desc0(1:kk-1)
                    desc_out(nDesc_lines) = desc1
                    desc0 = desc0(kk+2:)
                end if
            end do      
            !print *,"...done" 
            return
        end subroutine breakUpLine    
    
        subroutine report0(this,u,o)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^
            type(commandLineArguments),intent(inout)     ::      this
            integer,intent(in),optional         ::      u,o
            integer     ::      uu,oo
            
            character(len=LIB_CLA_DESCLEN),dimension(:),pointer         ::      desc
            integer                                     ::      nDesc
            
            logical     ::      ok,ok_tmp
            integer     ::      ii,jj,nn,kk
            character(len=256)  ::      cla
            
            logical,dimension(2)    ::      loga
            integer                 ::      nloga
            character(len=256)      ::      cla_tmp
            uu = 6 ; if (present(u)) uu = u
            oo = 0 ; if (present(o)) oo = o 
            
            
            
            if ( this%nArgs0 == 0 ) return      
            
            nn = command_argument_count() 
            ok = .not. hasHelpArgument(this)
!             
            do jj = this%ignore+1,this%nArgs
                if ( (.not. this%arg_optional(jj)).and.(.not. this%arg_optional_set(jj) ) ) ok = .false.
            end do
            

            kk = this%ignore 
            do
                kk = kk + 1
                if (kk>nn) exit
                call get_command_argument( kk,cla )
                if ( (trim(cla)=="-h").or.(trim(cla)=="-help").or.(trim(cla)=="--h").or.(trim(cla)=="--help") ) cycle   !   always accept help agument.
                if ( (trim(cla)=="-v").or.(trim(cla)=="-version").or.(trim(cla)=="--v").or.(trim(cla)=="--version") ) then
                    call reportProgramDescription(this,uu,oo)
                    if (len_trim(this%program_version) > 0) then
                        write (unit=uu,fmt='(a)') repeat(" ",oo+4)//"version "//trim(this%program_version)     
                    else
                        write (unit=uu,fmt='(a)') repeat(" ",oo+4)//"version number not set"
                    end if  
                    stop    
                end if
                ok_tmp = .false.
                do jj = 1,this%nArgs
                    if ( ( (trim(cla) == "-"//trim(this%arg_key(jj))).or.(trim(cla) == "-no"//trim(this%arg_key(jj))) ) .and. (this%arg_type(jj)==LIB_CLA_BOOLA) )  then          !   have identified a "-key T,F" or "-key" or "-nokey T,F" or "-nokey"
                        ok_tmp = .true.         !   this is a valid argument
                        if (kk<nn) then
                            !   check if next argument is "T,F"
                            call get_command_argument( kk+1,cla_tmp )
                            call parse( cla_tmp,loga,nloga )
                            if (nloga>0) kk = kk + 1                            !   skip next argument, it is the value. 
                        end if
                    else if (trim(cla) == "-"//trim(this%arg_key(jj))) then          !   have identified a key
                        if ( this%arg_type(jj)/=LIB_CLA_BOOL ) kk = kk + 1      !   skip next argument, it is the value.
                        ok_tmp = (kk<=nn)       !   ok if value exists.
                        exit
                    else if (( this%arg_type(jj)==LIB_CLA_BOOL ).and. (trim(cla) == "-no"//trim(this%arg_key(jj)))) then          !   have identified a -nokey
                        ok_tmp = (kk<=nn)       !   ok if value exists.
                        exit
                    end if
                end do
                !print *,"arg ",jj,kk," """//trim(cla)//""" ",ok_tmp
                if (.not. ok_tmp) then
                    print *,"Lib_CommandLineArguments::report0 warning - option """//trim(cla)//""" not recognised, skipping..."                                        
                    print *,""
                end if
            end do
                        
            
            

            if (ok) return           
            
            write (unit=uu,fmt='(a)') repeat(" ",oo)//"usage"
            write (unit=uu,fmt='(a)') repeat(" ",oo)//repeat("^",5)
                
            if (this%nCategories>0) then
                do ii = 1,this%nCategories 
                    write (unit=uu,fmt='(a)') repeat(" ",oo+4)//trim(this%cat_description(ii))
                    do jj = 1,this%nArgs
                        if (this%arg_cat(jj) /= ii) cycle
                        
                        write (unit=uu,fmt='(a)',advance="no") repeat(" ",oo+8)
                        if (this%arg_optional(jj)) write (unit=uu,fmt='(a)',advance="no") "["
                        write (unit=uu,fmt='(a)',advance="no") "-"//trim(this%arg_key(jj))
                        
                        if (this%arg_optional(jj)) then
                            select case(this%arg_type(jj))  
                                case(LIB_CLA_BOOL )
                                    write (unit=uu,fmt='(a)',advance="no") " ]              "
                                case(LIB_CLA_INT  )                                        
                                    write (unit=uu,fmt='(a)',advance="no") " <int> ]        "
                                case(LIB_CLA_REAL )                                        
                                    write (unit=uu,fmt='(a)',advance="no") " <float> ]      "
                                case(LIB_CLA_CHAR )                                        
                                    write (unit=uu,fmt='(a)',advance="no") " <char> ]       "
                                case(LIB_CLA_BOOLA)                                        
                                    write (unit=uu,fmt='(a)',advance="no") " <bool_array> ] "
                                case(LIB_CLA_INTA )                                        
                                    write (unit=uu,fmt='(a)',advance="no") " <int_array> ]  "
                                case(LIB_CLA_REALA)                                        
                                    write (unit=uu,fmt='(a)',advance="no") " <float_array> ]"
                                case(LIB_CLA_CHARA)                                        
                                    write (unit=uu,fmt='(a)',advance="no") " <char_array> ] "
                                
                            end select
                        else
                            select case(this%arg_type(jj))  
                                case(LIB_CLA_BOOL )
                                    write (unit=uu,fmt='(a)',advance="no") "               "
                                case(LIB_CLA_INT  )                                        
                                    write (unit=uu,fmt='(a)',advance="no") " <int>         "
                                case(LIB_CLA_REAL )                                        
                                    write (unit=uu,fmt='(a)',advance="no") " <float>       "
                                case(LIB_CLA_CHAR )                                        
                                    write (unit=uu,fmt='(a)',advance="no") " <char>        "
                                case(LIB_CLA_BOOLA)                                        
                                    write (unit=uu,fmt='(a)',advance="no") " <bool_array>  "
                                case(LIB_CLA_INTA )                                        
                                    write (unit=uu,fmt='(a)',advance="no") " <int_array>   "
                                case(LIB_CLA_REALA)                                        
                                    write (unit=uu,fmt='(a)',advance="no") " <float_array> "
                                case(LIB_CLA_CHARA)                                        
                                    write (unit=uu,fmt='(a)',advance="no") " <char_array> ] "
                                
                            end select
                        end if                    
                       
                        call breakUpLine( this%arg_description(jj),nDesc,desc )

                        
                                                
                        write (unit=uu,fmt='(a)',advance="yes") trim(desc(1))
                        do kk = 2,nDesc
                            write (unit=uu,fmt='(a)',advance="yes") repeat(" ",oo+32)//trim(desc(kk))
                        end do                                            
                        deallocate(desc)                               
                    end do
                    write (unit=uu,fmt='(a)') ""
                end do
                
            else

                do jj = 1,this%nArgs
                        
                    write (unit=uu,fmt='(a)',advance="no") repeat(" ",oo+8)
                    if (this%arg_optional(jj)) write (unit=uu,fmt='(a)',advance="no") "["
                    write (unit=uu,fmt='(a)',advance="no") " -"//trim(this%arg_key(jj))
                    
                    if (this%arg_optional(jj)) then
                        select case(this%arg_type(jj))  
                            case(LIB_CLA_BOOL )
                                write (unit=uu,fmt='(a)',advance="no") " ]              "
                            case(LIB_CLA_INT  )                                        
                                write (unit=uu,fmt='(a)',advance="no") " <int> ]        "
                            case(LIB_CLA_REAL )                                        
                                write (unit=uu,fmt='(a)',advance="no") " <float> ]      "
                            case(LIB_CLA_CHAR )                                        
                                write (unit=uu,fmt='(a)',advance="no") " <char> ]       "
                            case(LIB_CLA_BOOLA)                                        
                                write (unit=uu,fmt='(a)',advance="no") " <bool_array> ] "
                            case(LIB_CLA_INTA )                                        
                                write (unit=uu,fmt='(a)',advance="no") " <int_array> ]  "
                            case(LIB_CLA_REALA)                                        
                                write (unit=uu,fmt='(a)',advance="no") " <float_array> ]"
                            case(LIB_CLA_CHARA)                                        
                                write (unit=uu,fmt='(a)',advance="no") " <char_array> ] "
                            
                        end select
                    else
                        select case(this%arg_type(jj))  
                            case(LIB_CLA_BOOL )
                                write (unit=uu,fmt='(a)',advance="no") "               "
                            case(LIB_CLA_INT  )                                        
                                write (unit=uu,fmt='(a)',advance="no") " <int>         "
                            case(LIB_CLA_REAL )                                        
                                write (unit=uu,fmt='(a)',advance="no") " <float>       "
                            case(LIB_CLA_CHAR )                                        
                                write (unit=uu,fmt='(a)',advance="no") " <char>        "
                            case(LIB_CLA_BOOLA)                                        
                                write (unit=uu,fmt='(a)',advance="no") " <bool_array>  "
                            case(LIB_CLA_INTA )                                        
                                write (unit=uu,fmt='(a)',advance="no") " <int_array>   "
                            case(LIB_CLA_REALA)                                        
                                write (unit=uu,fmt='(a)',advance="no") " <float_array> "
                            case(LIB_CLA_CHARA)                                        
                                write (unit=uu,fmt='(a)',advance="no") " <char_array> ] "
                            
                        end select
                    end if                    
!                     
                    call breakUpLine( this%arg_description(jj),nDesc,desc )             
                    write (unit=uu,fmt='(a)',advance="yes") trim(desc(1))
                    do kk = 2,nDesc
                        write (unit=uu,fmt='(a)',advance="yes") repeat(" ",oo+32)//trim(desc(kk))
                    end do                                                    
                    deallocate(desc) 
                end do
                write (unit=uu,fmt='(a)') ""

            
            end if
 
            
            return
        end subroutine report0
             
    !---
    
        subroutine setIgnore( this,ignore )
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
            type(CommandLineArguments),intent(inout)        ::      this
            integer,intent(in)          ::      ignore
            this%ignore = ignore
            return
        end subroutine setIgnore
    
        subroutine setCategories( this,cat_description )
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
            type(CommandLineArguments),intent(inout)        ::      this
            character(len=*),dimension(:),intent(in)     ::      cat_description
            integer     ::  ii
            this%nCAtegories = size(cat_description)
            do ii = 1,this%nCategories
                this%cat_description(ii) = cat_description(ii)
            end do
            return
        end subroutine setCategories
            
        
        
        subroutine setProgramDescription( this,program_description )
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
            type(CommandLineArguments),intent(inout)        ::      this
            character(len=*),intent(in)     ::      program_description
            this%program_description = program_description            
            return
        end subroutine setProgramDescription
            
        subroutine setProgramVersion( this,program_version )
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
            type(CommandLineArguments),intent(inout)        ::      this
            character(len=*),intent(in)     ::      program_version
            this%program_version = program_version            
            return
        end subroutine setProgramVersion
            
    !---
    
    
        subroutine get00( this,arg_key,arg_optional,arg_description,arg_cat )
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
            type(CommandLineArguments),intent(inout)        ::      this
            logical,intent(in)              ::      arg_optional
            character(len=*),intent(in)     ::      arg_key
            character(len=*),intent(in)     ::      arg_description
            integer,intent(in),optional     ::      arg_cat
            
            this%nArgs = this%nArgs + 1
            if (this%nArgs > this%nArgs0) then
                print *,"Lib_CommandLineArguments::get00 error - attempt to specify too many arguments"
                stop
            end if
            
            if (this%nCategories == 0) then
                this%arg_cat(this%nArgs) = 0
            else 
                if (present(arg_cat)) then
                    this%arg_cat(this%nArgs) = arg_cat
                else
                    print *,"Lib_CommandLineArguments::get00 error - category required "
                    stop
                end if
            end if
            this%arg_key(this%nArgs) = arg_key
            this%arg_optional(this%nArgs) = arg_optional
            this%arg_description(this%nArgs) = arg_description
            return
        end subroutine get00
            
            

            
        subroutine get1( this,arg_key,arg_value,arg_optional,arg_description,arg_cat )
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
            type(CommandLineArguments),intent(inout)        ::      this
            character(len=*),intent(in)     ::      arg_key
            logical,intent(inout)           ::      arg_value
            logical,intent(in)              ::      arg_optional
            character(len=*),intent(in)     ::      arg_description
            integer,intent(in),optional     ::      arg_cat
            character(len=LIB_CLA_DESCLEN)  ::      arg_default
            integer         ::      kk,nn
            character(len=256)              ::      cla
            
            if (present(arg_cat)) then
                call get00( this,arg_key,arg_optional,arg_description,arg_cat )
            else
                call get00( this,arg_key,arg_optional,arg_description )
            end if
            
            this%arg_type(this%nArgs) = LIB_CLA_BOOL
                        
            write (arg_default,fmt='(l2)') arg_value
            arg_default =  "[ default "//trim(adjustl(arg_default))//" ]"           
            kk = index(this%arg_description(this%nArgs),"\n")
            if (kk>0) then
                this%arg_description(this%nArgs) = trim(this%arg_description(this%nArgs))//"\n"//trim(arg_default)
            else
                this%arg_description(this%nArgs) = trim(this%arg_description(this%nArgs))//" "//trim(arg_default)
            end if
            
            nn = command_argument_count()
            this%arg_optional_set(this%nArgs) = .false.
            
            do kk = this%ignore+1,nn
                call get_command_argument( kk,cla )
                if (trim(cla) == "-"//trim(this%arg_key(this%nArgs))) then
                    arg_value = .true.
!                    print *,"Lib_CommandLineArguments::get1 info - set "//trim(this%arg_key(this%nArgs))//" to true"
                    this%arg_optional_set(this%nArgs) = .true.
                    exit
                else if (trim(cla) == "-no"//trim(this%arg_key(this%nArgs))) then
                    arg_value = .false.
!                    print *,"Lib_CommandLineArguments::get1 info - set "//trim(this%arg_key(this%nArgs))//" to false"
                    this%arg_optional_set(this%nArgs) = .true.
                    exit
                end if
            end do
            
            
            return
        end subroutine get1
            
            
        subroutine get2( this,arg_key,arg_value,arg_optional,arg_description,arg_cat )
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
            type(CommandLineArguments),intent(inout)        ::      this
            character(len=*),intent(in)     ::      arg_key
            integer,intent(inout)           ::      arg_value
            logical,intent(in)              ::      arg_optional
            character(len=*),intent(in)     ::      arg_description
            integer,intent(in),optional     ::      arg_cat
            character(len=LIB_CLA_DESCLEN)  ::      arg_default
            integer         ::      kk,nn
            character(len=256)              ::      cla
            logical         ::      ok
            
            if (present(arg_cat)) then
                call get00( this,arg_key,arg_optional,arg_description,arg_cat )
            else
                call get00( this,arg_key,arg_optional,arg_description )
            end if
            
            this%arg_type(this%nArgs) = LIB_CLA_INT
            
             
            if (arg_value/=LIB_CLA_NODEFAULT_I) then
                write (arg_default,fmt='(i8)') arg_value
                arg_default =  "[ default "//trim(adjustl(arg_default))//" ]"                       
                kk = index(this%arg_description(this%nArgs),"\n")
                if (kk>0) then
                    this%arg_description(this%nArgs) = trim(this%arg_description(this%nArgs))//"\n"//trim(arg_default)
                else
                    this%arg_description(this%nArgs) = trim(this%arg_description(this%nArgs))//" "//trim(arg_default)
                end if
            end if
                        
            nn = command_argument_count()
            this%arg_optional_set(this%nArgs) = .false.
            do kk = this%ignore+1,nn-1
                call get_command_argument( kk,cla )
                if (trim(cla) == "-"//this%arg_key(this%nArgs)) then
                    call get_command_argument( kk+1,cla )
                    call parse( cla,arg_value,ok )
!                    print *,"Lib_CommandLineArguments::get2 info - set "//trim(this%arg_key(this%nArgs))//" to ",arg_value
                    this%arg_optional_set(this%nArgs) = .true.                    
                    if (.not. ok) then
                        print *,"Lib_CommandLineArguments::get2 error - could not parse integer value for -"//trim(this%arg_key(this%nArgs))//" from """//trim(cla)//""""
                        this%arg_optional_set(this%nArgs) = .false.                    
                    end if
                    exit
                end if
            end do
            
            
            return
        end subroutine get2
            
        subroutine get3( this,arg_key,arg_value,arg_optional,arg_description,arg_cat )
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
            type(CommandLineArguments),intent(inout)        ::      this
            character(len=*),intent(in)     ::      arg_key
            real(kind=real64),intent(inout) ::      arg_value
            logical,intent(in)              ::      arg_optional
            character(len=*),intent(in)     ::      arg_description
            integer,intent(in),optional     ::      arg_cat
            character(len=LIB_CLA_DESCLEN)  ::      arg_default
            integer         ::      kk,nn
            character(len=256)              ::      cla
            logical         ::      ok
            if (present(arg_cat)) then
                call get00( this,arg_key,arg_optional,arg_description,arg_cat )
            else
                call get00( this,arg_key,arg_optional,arg_description )
            end if
            
            this%arg_type(this%nArgs) = LIB_CLA_REAL
            
             
            if (arg_value/=LIB_CLA_NODEFAULT_R) then
                write (arg_default,fmt='(g16.4)') arg_value
                arg_default =  "[ default "//trim(adjustl(arg_default))//" ]"           
                kk = index(this%arg_description(this%nArgs),"\n")
                if (kk>0) then
                    this%arg_description(this%nArgs) = trim(this%arg_description(this%nArgs))//"\n"//trim(arg_default)
                else
                    this%arg_description(this%nArgs) = trim(this%arg_description(this%nArgs))//" "//trim(arg_default)
                end if
            end if
            
            nn = command_argument_count()
            this%arg_optional_set(this%nArgs) = .false.
            
            do kk = this%ignore+1,nn-1
                call get_command_argument( kk,cla )
                if (trim(cla) == "-"//this%arg_key(this%nArgs)) then
                    call get_command_argument( kk+1,cla )
                    call parse( cla,arg_value,ok )
!                     print *,"Lib_CommandLineArguments::get3 info - set "//trim(this%arg_key(this%nArgs))//" to ",arg_value
                    this%arg_optional_set(this%nArgs) = .true.
                    if (.not. ok) then
                        print *,"Lib_CommandLineArguments::get3 error - could not parse real value for -"//trim(this%arg_key(this%nArgs))//" from """//trim(cla)//""""
                        this%arg_optional_set(this%nArgs) = .false.                        
                    end if
                    exit
                end if
            end do
            
            
            return
        end subroutine get3
            
        subroutine get4( this,arg_key,arg_value,arg_optional,arg_description,arg_cat )
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
            type(CommandLineArguments),intent(inout)        ::      this
            character(len=*),intent(in)     ::      arg_key
            character(len=*),intent(inout)  ::      arg_value
            logical,intent(in)              ::      arg_optional
            character(len=*),intent(in)     ::      arg_description
            integer,intent(in),optional     ::      arg_cat
            character(len=LIB_CLA_DESCLEN)  ::      arg_default
            integer         ::      kk,nn
            character(len=256)              ::      cla
            
            if (present(arg_cat)) then
                call get00( this,arg_key,arg_optional,arg_description,arg_cat )
            else
                call get00( this,arg_key,arg_optional,arg_description )
            end if
            
            this%arg_type(this%nArgs) = LIB_CLA_CHAR
            
            if (len_trim(arg_value)>0) then
                arg_default = "[ default """//trim(adjustl(arg_value))//""" ]"
                kk = index(this%arg_description(this%nArgs),"\n")
                if (kk>0) then
                    this%arg_description(this%nArgs) = trim(this%arg_description(this%nArgs))//"\n"//trim(arg_default)
                else
                    this%arg_description(this%nArgs) = trim(this%arg_description(this%nArgs))//" "//trim(arg_default)
                end if
            end if
            
                        
            nn = command_argument_count()
            this%arg_optional_set(this%nArgs) = .false.
            
            do kk = this%ignore+1,nn 
                call get_command_argument( kk,cla )
                if (trim(cla) == "-"//this%arg_key(this%nArgs)) then
                    call get_command_argument( kk+1,arg_value )
!                     print *,"Lib_CommandLineArguments::get4 info - set "//trim(this%arg_key(this%nArgs))//" to """,trim(arg_value),""""
                    this%arg_optional_set(this%nArgs) = .true.  
                    exit
                end if
            end do
            
            
            return
        end subroutine get4
            
            
        subroutine get5( this,arg_key,arg_value,n_value,arg_optional,arg_description,arg_cat )
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
            type(CommandLineArguments),intent(inout)        ::      this
            character(len=*),intent(in)     ::      arg_key
            integer,dimension(:),intent(inout)           ::      arg_value
            integer,intent(inout)           ::      n_value
            logical,intent(in)              ::      arg_optional
            character(len=*),intent(in)     ::      arg_description
            integer,intent(in),optional     ::      arg_cat
            character(len=LIB_CLA_DESCLEN)  ::      arg_default,arg_default_tmp
            integer         ::      kk,nn,n_val
            character(len=256)              ::      cla
             
            
            if (present(arg_cat)) then
                call get00( this,arg_key,arg_optional,arg_description,arg_cat )
            else
                call get00( this,arg_key,arg_optional,arg_description )
            end if
            
            this%arg_type(this%nArgs) = LIB_CLA_INTA
            
            if (arg_value(1)/=LIB_CLA_NODEFAULT_I) then
                arg_default = ""
                if (size(arg_value)>=1) then
                    write (arg_default,fmt='(i8)') arg_value(1)
                end if
                if (size(arg_value)>=2) then    
                    write (arg_default_tmp,fmt='(i8)') arg_value(2)
                    arg_default = trim(arg_default)//","//trim(adjustl(arg_default_tmp)) 
                end if   
                if (size(arg_value)>=3) then    
                    write (arg_default_tmp,fmt='(i8)') arg_value(3)
                    arg_default =  trim(arg_default)//","//trim(adjustl(arg_default_tmp)) 
                end if   
                if (size(arg_value)>=4) then                   
                    arg_default =  trim(arg_default)//",..."
                end if   
                if (len_trim(arg_default)>0) then
                    arg_default =  "[ default "//trim(adjustl(arg_default))//" ]"    
                    kk = index(this%arg_description(this%nArgs),"\n")
                    if (kk>0) then
                        this%arg_description(this%nArgs) = trim(this%arg_description(this%nArgs))//"\n"//trim(arg_default)
                    else
                        this%arg_description(this%nArgs) = trim(this%arg_description(this%nArgs))//" "//trim(arg_default)
                    end if
                end if
            end if
                        
            nn = command_argument_count()
            this%arg_optional_set(this%nArgs) = .false.
            do kk = this%ignore+1,nn-1
                call get_command_argument( kk,cla )
                if (trim(cla) == "-"//this%arg_key(this%nArgs)) then
                    call get_command_argument( kk+1,cla )
                    call parse( cla,arg_value,n_val )
!                     print *,"Lib_CommandLineArguments::get5 info - set "//trim(this%arg_key(this%nArgs))//" to ",arg_value(1:n_val)

                    this%arg_optional_set(this%nArgs) = .true.  
                    if (n_val < n_value) then
                        print *,"Lib_CommandLineArguments::get5 error - could not parse enough integer values for -"//trim(this%arg_key(this%nArgs))//" from """//trim(cla)//""""                        
                        this%arg_optional_set(this%nArgs) = .false.
                    end if
                    n_value = n_val              
                    exit
                end if
            end do
            
            
            return
        end subroutine get5
            
        subroutine get6( this,arg_key,arg_value,n_value,arg_optional,arg_description,arg_cat )
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
            type(CommandLineArguments),intent(inout)        ::      this
            character(len=*),intent(in)     ::      arg_key
            real(kind=real64),dimension(:),intent(inout) ::      arg_value
            integer,intent(inout)           ::      n_value
            logical,intent(in)              ::      arg_optional
            character(len=*),intent(in)     ::      arg_description
            integer,intent(in),optional     ::      arg_cat
            character(len=LIB_CLA_DESCLEN)  ::      arg_default,arg_default_tmp
            integer         ::      kk,nn,n_val
            character(len=256)              ::      cla
            
            if (present(arg_cat)) then
                call get00( this,arg_key,arg_optional,arg_description,arg_cat )
            else
                call get00( this,arg_key,arg_optional,arg_description )
            end if
            
            this%arg_type(this%nArgs) = LIB_CLA_REALA
 
            if (arg_value(1)/=LIB_CLA_NODEFAULT_R) then
                arg_default = ""
                if (size(arg_value)>=1) then
                    write (arg_default,fmt='(f16.4)') arg_value(1)
                end if
                if (size(arg_value)>=2) then    
                    write (arg_default_tmp,fmt='(f16.4)') arg_value(2)
                    arg_default = trim(arg_default)//","//trim(adjustl(arg_default_tmp)) 
                end if   
                if (size(arg_value)>=3) then    
                    write (arg_default_tmp,fmt='(f16.4)') arg_value(3)
                    arg_default =  trim(arg_default)//","//trim(adjustl(arg_default_tmp)) 
                end if   
                if (size(arg_value)>=4) then                   
                    arg_default =  trim(arg_default)//",..."
                end if   
                if (len_trim(arg_default)>0) then
                    arg_default =  "[ default "//trim(adjustl(arg_default))//" ]"    
                    kk = index(this%arg_description(this%nArgs),"\n")
                    if (kk>0) then
                        this%arg_description(this%nArgs) = trim(this%arg_description(this%nArgs))//"\n"//trim(arg_default)
                    else
                        this%arg_description(this%nArgs) = trim(this%arg_description(this%nArgs))//" "//trim(arg_default)
                    end if
                end if            
            end if   
           
            nn = command_argument_count()
            this%arg_optional_set(this%nArgs) = .false.
            
            do kk = this%ignore+1,nn-1
                call get_command_argument( kk,cla )
                if (trim(cla) == "-"//this%arg_key(this%nArgs)) then
                    call get_command_argument( kk+1,cla )
                    call parse( cla,arg_value,n_val )
                    this%arg_optional_set(this%nArgs) = .true.  
!                     print *,"Lib_CommandLineArguments::get6 info - set "//trim(this%arg_key(this%nArgs))//" to ",arg_value(1:n_val)
                    if (n_val < n_value) then
                        print *,"Lib_CommandLineArguments::get6 error - could not parse enough real values for -"//trim(this%arg_key(this%nArgs))//" from """//trim(cla)//""""                        
                        this%arg_optional_set(this%nArgs) = .false.
                    end if
                    n_value = n_val              
                    exit
                end if
            end do
            
            
            return
        end subroutine get6
            
        subroutine get7( this,arg_key,arg_value,n_value,arg_optional,arg_description,arg_cat )
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
            type(CommandLineArguments),intent(inout)        ::      this
            character(len=*),intent(in)     ::      arg_key
            character(len=*),dimension(:),intent(inout) ::      arg_value
            integer,intent(inout)           ::      n_value
            logical,intent(in)              ::      arg_optional
            character(len=*),intent(in)     ::      arg_description
            integer,intent(in),optional     ::      arg_cat
            character(len=LIB_CLA_DESCLEN)  ::      arg_default
            integer         ::      kk,nn,n_val
            character(len=256)              ::      cla
            
            if (present(arg_cat)) then
                call get00( this,arg_key,arg_optional,arg_description,arg_cat )
            else
                call get00( this,arg_key,arg_optional,arg_description )
            end if
            
            this%arg_type(this%nArgs) = LIB_CLA_CHARA
 
            if (len(arg_value(1))>0) then
                arg_default = ""
                if (size(arg_value)>=1) then
                    arg_default = arg_value(1)
                end if
                if (size(arg_value)>=2) then    
                    arg_default = trim(arg_default)//","//arg_value(2)
                end if   
                if (size(arg_value)>=3) then    
                    arg_default =  trim(arg_default)//","//arg_value(3)
                end if   
                if (size(arg_value)>=4) then                   
                    arg_default =  trim(arg_default)//",..."
                end if   
                if (len_trim(arg_default)>0) then
                    arg_default =  "[ default "//trim(adjustl(arg_default))//" ]"    
                    kk = index(this%arg_description(this%nArgs),"\n")
                    if (kk>0) then
                        this%arg_description(this%nArgs) = trim(this%arg_description(this%nArgs))//"\n"//trim(arg_default)
                    else
                        this%arg_description(this%nArgs) = trim(this%arg_description(this%nArgs))//" "//trim(arg_default)
                    end if
                end if            
            end if   
           
            nn = command_argument_count()
            this%arg_optional_set(this%nArgs) = .false.
            
            do kk = this%ignore+1,nn-1
                call get_command_argument( kk,cla )
                if (trim(cla) == "-"//this%arg_key(this%nArgs)) then
                    call get_command_argument( kk+1,cla )
                    call parse( cla,arg_value,n_val )
                    this%arg_optional_set(this%nArgs) = .true.  
                    if (n_val < n_value) then
                        print *,"Lib_CommandLineArguments::get7 error - could not parse enough char values for -"//trim(this%arg_key(this%nArgs))//" from """//trim(cla)//""""                        
                        this%arg_optional_set(this%nArgs) = .false.
                    end if
                    n_value = n_val              
                    exit
                end if
            end do
            
            
            return
        end subroutine get7
        

        subroutine get8( this,arg_key,arg_value,n_value,arg_optional,arg_description,arg_cat )
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
            type(CommandLineArguments),intent(inout)        ::      this
            character(len=*),intent(in)     ::      arg_key
            logical,dimension(:),intent(inout) ::      arg_value
            integer,intent(inout)           ::      n_value
            logical,intent(in)              ::      arg_optional
            character(len=*),intent(in)     ::      arg_description
            integer,intent(in),optional     ::      arg_cat
            character(len=LIB_CLA_DESCLEN)  ::      arg_default,arg_default_tmp
            integer         ::      kk,nn,n_val
            character(len=256)              ::      cla
            
            if (present(arg_cat)) then
                call get00( this,arg_key,arg_optional,arg_description,arg_cat )
            else
                call get00( this,arg_key,arg_optional,arg_description )
            end if
            
            this%arg_type(this%nArgs) = LIB_CLA_BOOLA
 
            arg_default = ""
            if (size(arg_value)>=1) then
                write (arg_default,fmt='(l2)') arg_value(1)
            end if
            if (size(arg_value)>=2) then    
                write (arg_default_tmp,fmt='(l2)') arg_value(2)
                arg_default = trim(arg_default)//","//trim(adjustl(arg_default_tmp)) 
            end if   
            if (size(arg_value)>=3) then    
                write (arg_default_tmp,fmt='(l2)') arg_value(3)
                arg_default =  trim(arg_default)//","//trim(adjustl(arg_default_tmp)) 
            end if   
            if (size(arg_value)>=4) then                   
                arg_default =  trim(arg_default)//",..."
            end if   
            if (len_trim(arg_default)>0) then
                arg_default =  "[ default "//trim(adjustl(arg_default))//" ]"    
                kk = index(this%arg_description(this%nArgs),"\n")
                if (kk>0) then
                    this%arg_description(this%nArgs) = trim(this%arg_description(this%nArgs))//"\n"//trim(arg_default)
                else
                    this%arg_description(this%nArgs) = trim(this%arg_description(this%nArgs))//" "//trim(arg_default)
                end if
            end if            
           
            nn = command_argument_count()
            this%arg_optional_set(this%nArgs) = .false.
            
            do kk = this%ignore+1,nn-1
                call get_command_argument( kk,cla )
                if (trim(cla) == "-"//this%arg_key(this%nArgs)) then
                    call get_command_argument( kk+1,cla )
                    call parse( cla,arg_value,n_val )
                    
                    if (n_val == 0) then
                        !   have entered a single logical value "-fruit" rather than an array "-fruit T,T,F,T"
                        arg_value(1) = .true.
                        n_val = 1
                    else
                        this%arg_optional_set(this%nArgs) = .true.  
                        if (n_val < n_value) then
                            print *,"Lib_CommandLineArguments::get8 error - could not parse enough logical values for -"//trim(this%arg_key(this%nArgs))//" from """//trim(cla)//""""                        
                            this%arg_optional_set(this%nArgs) = .false.
                        end if                                 
                    end if
                    n_value = n_val    
                    exit
                else if (trim(cla) == "-no"//this%arg_key(this%nArgs)) then
                    call get_command_argument( kk+1,cla )
                    call parse( cla,arg_value,n_val )                    
                    if (n_val == 0) then
                        !   have entered a single logical value "-nofruit" rather than an array "-nofruit T,T,F,T"
                        arg_value(1) = .false.
                        n_val = 1
                    else
                        this%arg_optional_set(this%nArgs) = .true.  
                        arg_value = .not. arg_value
                        if (n_val < n_value) then
                            print *,"Lib_CommandLineArguments::get8 error - could not parse enough logical values for -"//trim(this%arg_key(this%nArgs))//" from """//trim(cla)//""""                        
                            this%arg_optional_set(this%nArgs) = .false.
                        end if                                 
                    end if
                    n_value = n_val    
                    exit                    
                end if
            end do
            
            
            return
        end subroutine get8
                    
             
            
    end module Lib_CommandLineArguments

!     
! !   gfortran -ffree-line-length-256 ${MYF90LIB}/NBAX_StringTokenizers.f90 Lib_CommandLineArguments.f90 -Og -g 
!     
!     program testLib_CommandLineArguments
! !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
!         use iso_fortran_env
!         use Lib_CommandLineArguments
!         implicit none
!         
!         
!         type(CommandLineArguments)      ::      cla
!         
!         character(len=256)      ::  filename = "test.dat"
!         integer                 ::  vari = 12345
!         real(kind=real64)       ::  varr = 1.0d0
!         logical                 ::  varl = .true.
!         integer,dimension(3)    ::  varia = (/1,2,3/)
!         real(kind=real64),dimension(4)    ::  varra = (/3.0d0,2.0d0,1.0d0,0.0d0/)
!         integer                 ::  nvaria = 3 , nvarra = 0
!         
!         cla = CommandLineArguments_ctor(10)
!         call setCategories(cla,(/"category 1","category 2"/))
!         
!         call get(cla,"f",filename,LIB_CLA_REQUIRED,"input filename",1)
!         call get(cla,"i",vari,LIB_CLA_REQUIRED,"integer variable",1)
!         call get(cla,"r",varr,LIB_CLA_REQUIRED,"real variable",1)
!         call get(cla,"l",varl,LIB_CLA_OPTIONAL,"option",2)
!         call get(cla,"int_array",varia,nvaria,LIB_CLA_OPTIONAL,"integer array ",2)
!         call get(cla,"real_array",varra,nvarra,LIB_CLA_OPTIONAL,"real array",2)
!         
!         
!         call report(cla)
!         
!         
!         print *,"filename = """//trim(filename)//""""
!         print *,"vari     = ",vari
!         print *,"varr     = ",varr
!         print *,"varl     = ",varl
!         print *,"varia    = ",varia
!         print *,"varra    = ",varra
!         
!         if (.not. allRequiredArgumentsSet(cla)) print *,"required arguments unset"
!         
!         call delete(cla)
!         print *,""
!         print *,"done"
!         print *,""
!     end program testLib_CommandLineArguments
!         
!         
        
        
        
        
        
        
        
        