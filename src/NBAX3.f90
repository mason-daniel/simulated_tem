
    module Lib_NBAX
!---^^^^^^^^^^^^^^^
!*      version 2.0

!*
!*    This is the Fortran version of my NoBrainer Api for Xml
!*    - a very simple set of subroutines to read and write xml
!*    from disk. You don't need an extra parser, a document object
!*    model or any additional rubbish to use this: instead my
!*    idea is to get started with *writing* good xml immediately
!*    and to be able to read in what you wrote at least.
!
!*    In the future someone may decide to implement a *real(kind=real64)* xml
!*    library. If so, then all the data written is accessible, you
!*    merely gain a little extra read functionality.
!
!*    for my purposes XML data structure is as follows
!*
!*v         < (name) [attribute-list] >
!*v             [text]
!*v             [ < (child1) ... > ]
!*v             [ < (child2) ... > ]
!*v         < (/name) >
!*
!*    or
!*v         < (name) [attribute-list] />
!*
!

!
!*    where attribute list is a comma or space delineated list of
!*    simple keypairs ( key="value" )
!


!*    A test program is attached.
!
!*    Author      Daniel Mason
!*    Version     2.01                  !   2.01 includes Xinsert support and is thread safe
!*    Revision    Aug 2010
!


        use iso_fortran_env
        use NBAX_StringTokenizers

        implicit none
        private


!-------        standard member functions
        public      ::      NBAX_ctor
        public      ::      Keypair_ctor
        public      ::      delete
        public      ::      report
!
! !-------        operators
        public      ::      operator(==)
        public              ::      assignment(=)
!
! !-------        io
        public      ::      input
        public      ::      output
!
! !-------        accessors/mutators
        public      ::      setName
        public      ::      getName
        public      ::      getKey
        public      ::      getValue
        public      ::      toString
        public      ::      KeyPairfromString
        public      ::      KeyfromString
        public      ::      extractKeyPair
        public      ::      ValuefromString
        public      ::      isDefined
        public      ::      hasAttributes
        public      ::      getNAttributes
        public      ::      getAttribute
        public      ::      getAttributeValue
        public      ::      changeAttributeValue
        public      ::      addAttribute
        public      ::      hasChildren
        public      ::      getNChildren
        public      ::      getChild
        public      ::      addChild
        public      ::      hasText
        public      ::      getNText
        public      ::      getText
        public      ::      addText
        public      ::      cutText
        public      ::      cutChildren
        public      ::      cutAttributes
        
        



!-------        public types

        integer,public,parameter                ::      KEYPAIR_MAX_STRING = 256
        character(len=1),public,parameter       ::      SEPARATOR="/"

        type,public         ::      KeyPair
            private
            character(len=KEYPAIR_MAX_STRING)   ::      key
            character(len=KEYPAIR_MAX_STRING)   ::      value
        end type

        type,public :: NBAX
            private
            integer                                     ::  nContent,N_CONTENT
            integer                                     ::  nChildren,N_CHILDREN
            integer                                     ::  nTextLines,N_TEXTLINES
            integer                                     ::  nAttributes,N_ATTRIBUTES
            integer,dimension(:),pointer                ::  content_type
            character(len=2048)                         ::  name
            character(len=2048),dimension(:),pointer    ::  textLines
            type(KeyPair),dimension(:),pointer          ::  attributes
            type(NBAX),dimension(:),pointer             ::  children
        end type


!-------        interface specifications

        interface NBAX_ctor
            module procedure    NBAX_Constructor0
            module procedure    NBAX_Constructor1
        end interface
        
        interface   Keypair_ctor
            module procedure    KeyPairNull
            module procedure    KeyPairConstructor1
            module procedure    KeyPairConstructor2
            module procedure    KeyPairConstructor3
            module procedure    KeyPairConstructor4
            module procedure    KeyPairConstructor5
        end interface
        
!
        interface delete
            module procedure    delete1
            module procedure    deleteKeyPair
        end interface

!-------

        interface operator(==)
            module procedure    NBAXequalsA
            module procedure        KPequalsKP
            module procedure        KPequalsA            
        end interface

!-------

        interface input
            module procedure    inputNBAX0
            module procedure    inputNBAX
        end interface

        interface output
            module procedure    outputNBAX0
            module procedure    outputNBAX1
            module procedure    outputNBAX2
        end interface

!-------

        interface setName
            module procedure    setNameNBAX
        end interface

        interface getName
            module procedure    getNameNBAX
        end interface

!-------

        interface getChild
            module procedure    getChildNBAX0
            module procedure    getChildNBAX1
            module procedure    getChildNBAX2
        end interface
!
        interface getNChildren
            module procedure    getNChildren0
            module procedure    getNChildren1
        end interface


        interface getAttributeValue
            module procedure    getValueNBAX0
            module procedure    getValueNBAX1
            module procedure    getValueNBAX2
            !module procedure    getValueNBAX3
            module procedure    getValueNBAX4
        end interface

        interface changeAttributeValue
            module procedure    changeValueNBAX0
            module procedure    changeValueNBAX1
            module procedure    changeValueNBAX2
            !module procedure    changeValueNBAX3
            module procedure    changeValueNBAX4
        end interface


        interface addAttribute
            module procedure    addAttributeNBAX0
            module procedure    addAttributeNBAX1
            module procedure    addAttributeNBAX2
            module procedure    addAttributeNBAX3
            module procedure    addAttributeNBAX4
            module procedure    addAttributeNBAX5
        end interface
        

!-------

        interface   getKey
            module procedure    getKey1
        end interface

        interface   getValue
            module procedure    getValue1
            module procedure    getValue2
        end interface

        interface   toString
            module procedure    toStringKP
        end interface

        interface assignment(=)
            module procedure        assignStringFromKP
        end interface

        interface isDefined
            module procedure        isDefinedKP
        end interface


        interface report
            module procedure        outputKP
        end interface
        
        
        
        

!-------        PRIVATE DATA MEMBERS
        integer,private,parameter               ::  NBAX_TEXT               = 1
        integer,private,parameter               ::  NBAX_OPEN               = 2
        integer,private,parameter               ::  NBAX_OPENCLOSE          = 3
        integer,private,parameter               ::  NBAX_CLOSE              = 4
        integer,private,parameter               ::  NBAX_OPENCOMMENT        = 5
        integer,private,parameter               ::  NBAX_OPENCLOSECOMMENT   = 6
        integer,private,parameter               ::  NBAX_CLOSECOMMENT       = 7
        integer,private,parameter               ::  NBAX_XMLSPEC            = 8
        integer,private,parameter               ::  NBAX_OPENDOCTYPE        = 9
        integer,private,parameter               ::  NBAX_OPENCLOSEDOCTYPE   = 10
        integer,private,parameter               ::  NBAX_CLOSEDOCTYPE       = 11

        integer,private,parameter               ::  NBAX_STARTNATT          = 4
        integer,private,parameter               ::  NBAX_STARTNTXT          = 4
        integer,private,parameter               ::  NBAX_STARTNCHL          = 4
        
        character(len=KEYPAIR_MAX_STRING),private,parameter ::  KP_TRUE = "true";
        character(len=KEYPAIR_MAX_STRING),private,parameter ::  KP_FALSE = "false";
        

    contains
!---^^^^^^^^

        recursive function NBAX_null(N_CH,N_TE,N_AT) result(this)
!-------^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
        !*  Constructs an empty xml wrapper with adequate storage
            integer,intent(in)              ::      N_CH,N_TE,N_AT
            type(NBAX)                      ::      this
            integer                         ::      ii
            this%name = ""
            this%N_CONTENT = N_CH + N_TE
            this%N_CHILDREN = N_CH
            this%N_TEXTLINES = N_TE
            this%N_ATTRIBUTES = N_AT
            if (this%N_CONTENT>0) then
                allocate(this%content_type(this%N_CONTENT))
                do ii = 1,this%N_CONTENT
                    this%content_type(ii) = 0
                end do
            else
                nullify(this%content_type)
            end if
            if (this%N_CHILDREN>0) then
                allocate(this%children(this%N_CHILDREN))
                do ii = 1,this%N_CHILDREN
                    this%children(ii) = NBAX_null(0,0,0)
                end do
            else
                nullify(this%children)
            end if
            if (this%N_TEXTLINES>0) then
                allocate(this%textLines(this%N_TEXTLINES))
                do ii = 1,this%N_TEXTLINES
                    this%textLines(ii) = ""
                end do
            else
                nullify(this%textLines)
            end if
            if (this%N_ATTRIBUTES>0) then
                allocate(this%attributes(this%N_ATTRIBUTES))
                do ii = 1,this%N_ATTRIBUTES
                    this%attributes(ii) = KeyPair_ctor()
                end do
            else
                nullify(this%attributes)
            end if
            this%nTextLines = 0
            this%nChildren = 0
            this%nContent = 0
            this%nAttributes = 0
            this%name = ""
            return
        end function NBAX_null


        function NBAX_constructor0() result(this)
!-------^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
        !*  Constructs an empty xml wrapper.
            type(NBAX)                      ::      this
            this = NBAX_null( 0,0,0 )
            return
        end function NBAX_constructor0



        function NBAX_constructor1(a) result(this)
!-------^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
        !*  Constructs an xml wrapper named a.
            character(len=*),intent(in)     ::      a
            type(NBAX)                      ::      this
!             print *,"NBAX_constructor1 ",trim(a)
            this = NBAX_null( NBAX_STARTNCHL,NBAX_STARTNTXT,NBAX_STARTNATT )
            this%name = trim(a)
            return
        end function NBAX_constructor1


!-------

        recursive subroutine delete1(this)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
            type(NBAX),intent(inout)         ::      this
            integer         ::      ii

!D$         print *,"          deleting"
!D$         call outputNBAX1(this,6,10)

            if (this%N_ATTRIBUTES > 0) then
                do ii = 1,this%nAttributes
                    call delete( this%attributes(ii) )
                end do
                deallocate(this%attributes)
            end if

            if (this%N_TEXTLINES > 0) deallocate(this%textLines)

            if (this%N_CHILDREN > 0) then
                do ii = 1,this%nChildren
                    call delete1( this%children(ii) )
                end do
                deallocate(this%children)
            end if

            if (this%N_CONTENT > 0) deallocate(this%content_type)

            this = NBAX_null( 0,0,0 )

            return
        end subroutine delete1

    !---


!-------

        subroutine reallocAttributes(this,n)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      ensure there is sufficient memory for n attributes
            type(NBAX),intent(inout)                ::      this
            integer,intent(in)                      ::      n
            type(KeyPair),dimension(:),pointer      ::      kp
            integer             ::      ii,mm
            if (n<=this%N_ATTRIBUTES) return
            mm = max( (3*this%N_ATTRIBUTES)/2,n )
            mm = max ( mm,NBAX_STARTNATT )
!             print *,"reallocAttributes ",mm
            allocate(kp(mm))
!             print *,"done"
            if (this%nAttributes > 0) then
                do ii = 1,this%nAttributes
                    kp(ii) = this%attributes(ii)
                end do
            end if
            do ii = this%nAttributes+1,mm
                kp(ii) = KeyPair_ctor()
            end do
            if (this%N_ATTRIBUTES > 0) deallocate(this%attributes)
            this%attributes => kp
            this%N_ATTRIBUTES = mm
            return
        end subroutine reallocAttributes

        subroutine reallocTextLines(this,n)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      ensure there is sufficient memory for n textLines
            type(NBAX),intent(inout)                    ::      this
            integer,intent(in)                          ::      n
            character(len=2048),dimension(:),pointer    ::      txt
            integer,dimension(:),pointer                ::      cnt
            integer             ::      ii,mm
            if (n<=this%N_TEXTLINES) return
            mm = max( (3*this%N_TEXTLINES)/2,n )
            mm = max ( mm,NBAX_STARTNTXT )
!             print *,"reallocTextLines ",mm
            allocate(txt(mm)) 
!             print *,"done"
            if (this%nTextlines > 0) then
                do ii = 1,this%nTextlines
                    txt(ii) = this%textlines(ii)
                end do
            end if
            do ii = this%nTextlines+1,mm
                txt(ii) = ""
            end do
            if (this%N_TEXTLINES > 0) deallocate(this%textlines)
            this%textlines => txt
            this%N_TEXTLINES = mm

            allocate(cnt(mm+this%N_CHILDREN))
            if (this%nContent > 0) then
                do ii = 1,this%nContent
                    cnt(ii) = this%content_type(ii)
                end do
            end if
            do ii = this%nContent+1,mm
                cnt(ii) = 0
            end do
            if (this%N_CONTENT > 0) deallocate(this%content_type)
            this%content_type => cnt
            this%N_CONTENT = this%N_TEXTLINES + this%N_CHILDREN
            return
        end subroutine reallocTextLines

        subroutine reallocChildren(this,n)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      ensure there is sufficient memory for n children
            type(NBAX),intent(inout)                    ::      this
            integer,intent(in)                          ::      n
            type(NBAX),dimension(:),pointer             ::      chl
            integer,dimension(:),pointer                ::      cnt
            integer             ::      ii,mm
            if (n<=this%N_CHILDREN) return
            mm = max( (3*this%N_CHILDREN)/2,n )
            mm = max ( mm,NBAX_STARTNCHL )
!             print *,"reallocChildren ",mm
            allocate(chl(mm))
!             print *,"done"
            if (this%nChildren > 0) then
                do ii = 1,this%nChildren
                    chl(ii) = this%children(ii)
                end do
            end if
            do ii = this%nChildren+1,mm
                chl(ii) = NBAX_null(0,0,0)
            end do
            if (this%N_CHILDREN > 0) deallocate(this%children)
            this%children => chl
            this%N_CHILDREN = mm

            allocate(cnt(mm+this%N_TEXTLINES))
            if (this%nContent > 0) then
                do ii = 1,this%nContent
                    cnt(ii) = this%content_type(ii)
                end do
            end if
            do ii = this%nContent+1,mm
                cnt(ii) = 0
            end do
            if (this%N_CONTENT > 0) deallocate(this%content_type)
            this%content_type => cnt
            this%N_CONTENT = this%N_TEXTLINES + this%N_CHILDREN
            return
        end subroutine reallocChildren







!-------


        subroutine addAttributeNBAX0(this,kk)
!-------^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
        !*  Adds attribute kk to this xml wrapper.
            type(NBAX),intent(inout)        ::      this
            type(KeyPair),intent(in)        ::      kk
            call reallocAttributes(this,this%nAttributes+1)
            this%nAttributes = this%nAttributes + 1
            this%attributes(this%nAttributes) = kk
            return
        end subroutine addAttributeNBAX0

        subroutine addAttributeNBAX1(this,key,value)
!-------^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
        !*  Adds an attribute"key"="value" to this xml wrapper.
            type(NBAX),intent(inout)        ::      this
            character(len=*),intent(in)     ::      key,value
            type(KeyPair)                   ::      kk
            kk = Keypair_ctor(key,Value)
            call addAttributeNBAX0(this,kk)
            return
        end subroutine addAttributeNBAX1

        subroutine addAttributeNBAX2(this,key,ii)
!-------^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
        !*  Adds an attribute"key"=ii to this xml wrapper.
            type(NBAX),intent(inout)        ::      this
            character(len=*),intent(in)     ::      key
            integer,intent(in)              ::      ii
            type(KeyPair)                   ::      kk
            kk = Keypair_ctor(key,ii)
            call addAttributeNBAX0(this,kk)
            return
        end subroutine addAttributeNBAX2

        subroutine addAttributeNBAX3(this,key,xx)
!-------^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
        !*  Adds an attribute"key"=xx to this xml wrapper.
            type(NBAX),intent(inout)        ::      this
            character(len=*),intent(in)     ::      key
            real(kind=real64),intent(in)                 ::      xx
            type(KeyPair)                   ::      kk
            kk = Keypair_ctor(key,xx)
            call addAttributeNBAX0(this,kk)
            return
        end subroutine addAttributeNBAX3

        subroutine addAttributeNBAX4(this,key,xx)
!-------^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
        !*  Adds an attribute"key"=xx to this xml wrapper.
            type(NBAX),intent(inout)        ::      this
            character(len=*),intent(in)     ::      key
            complex,intent(in)              ::      xx
            type(KeyPair)                   ::      kk
            kk = Keypair_ctor(key,xx)
            call addAttributeNBAX0(this,kk)
            return
        end subroutine addAttributeNBAX4

        subroutine addAttributeNBAX5(this,key,ll)
!-------^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
        !*  Adds an attribute"key"=ll to this xml wrapper.
            type(NBAX),intent(inout)        ::      this
            character(len=*),intent(in)     ::      key
            logical,intent(in)              ::      ll
            type(KeyPair)                   ::      kk
            if (ll) then
                kk = Keypair_ctor(key,"true")
            else
                kk = Keypair_ctor(key,"false")
            end if
            call addAttributeNBAX0(this,kk)
            return
        end subroutine addAttributeNBAX5

!-------

        subroutine addText(this,aa)
!-------^^^^^^^^^^^^^^^^^^^^^^^^^^^
        !*  Adds a line of text to this xml wrapper.
            type(NBAX),intent(inout)        ::      this
            character(len=*),intent(in)     ::      aa
            call reallocTextLines(this,this%nTextLines+1)
            this%nTextLines = this%nTextLines+1
            this%textLines(this%nTextLines) = trim(aa)
            this%nContent = this%nContent+1
            this%content_type(this%nContent) = -this%nTextLines
            return
        end subroutine addText

!-------

        subroutine addChild(this,child)
!-------^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
        !*  Adds a child to this xml wrapper.
            type(NBAX),intent(inout)        ::      this
            type(NBAX),intent(in)           ::      child
            call reallocChildren(this,this%nChildren+1)
            this%nChildren = this%nChildren+1
            this%children(this%nChildren) = child
            this%nContent = this%nContent+1
            this%content_type(this%nContent) = this%nChildren
            return
        end subroutine addChild

!-------

        subroutine outputNBAX0(this,u)
!-------^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
        !*  Dumps xml to unit u.
            type(NBAX),intent(in)           ::      this
            integer,intent(in),optional     ::      u
            integer     ::      uu,nn
            uu = 6
            if (present(u)) uu = u
            nn = 0
            call outputNBAX1(this,uu,nn)
            return
        end subroutine outputNBAX0

        subroutine outputNBAX2(this,u,topLevel)
!-------^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
        !*  Dumps xml to unit u.
            type(NBAX),intent(in)           ::      this
            integer,intent(in)              ::      u
            logical,intent(in)              ::      topLevel
            integer         ::      ioerr
            if (topLevel) then
                write(unit=u,fmt='(a)',iostat=ioerr) "<?xml version=""1.0"" standalone=""no""?>"
                if (ioerr == 0) call outputNBAX1(this,u,0)
            else
                call outputNBAX1(this,u,0)
            end if
            return
        end subroutine outputNBAX2


        recursive subroutine outputNBAX1(this,uu,nn)
!-------^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
            type(NBAX),intent(in)           ::      this
            integer,intent(in)              ::      uu,nn
            integer                     ::      ii,jj,ioerr
            character(len=2048)         ::      aa,bb
!            call nest(uu,nn)
            bb ="<"//trim(this%name)
            do ii = 1,this%nAttributes
                aa = toString(this%attributes(ii))
                bb = trim(bb)//" "//trim(aa)
            end do
            write(unit=uu,fmt="(a)",iostat=ioerr,advance="no") repeat(" ",nn)
            if (this%nContent > 0) then
                write(unit=uu,fmt="(a)",iostat=ioerr) trim(bb)//">"
                do ii = 1,this%nContent
                    jj = this%content_type(ii)
                    if (jj < 0) then
!                        call nest(uu,nn+1)
                        write(unit=uu,fmt="(a)",iostat=ioerr) repeat(" ",nn+2)//trim(this%textLines(-jj))
                    else if (jj > 0) then
                        call outputNBAX1(this%children(jj),uu,nn+2)
                    end if
                end do
!                call nest(uu,nn)
                write(unit=uu,fmt="(a)",iostat=ioerr) repeat(" ",nn)//"</"//trim(this%name)//">"
            else
                write(unit=uu,fmt="(a)",iostat=ioerr) trim(bb)//"/>"
            end if
            return
        end subroutine outputNBAX1

        

        
!         subroutine nest(uu,ii)
! !-------^^^^^^^^^^^^^^^^^^^^^^
!             integer,intent(in)      ::      uu,ii
!             integer     ::      jj,ioerr
!             do jj = 1,ii
!                 write(unit=uu,fmt='(a)',advance="no",iostat=ioerr)" "
!             end do
!             return
!         end subroutine nest

!-------

        pure function hasChildren(this) result(has)
!-------^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
        !*  Returns true if node has children.
            type(NBAX),intent(in)           ::      this
            logical     ::      has
            has = this%nChildren > 0
            return
        end function hasChildren

        pure function hasAttributes(this) result(has)
!-------^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
        !*  Returns true if this xml wrapper has any attributes.
            type(NBAX),intent(in)           ::      this
            logical     ::      has
            has = this%nAttributes > 0
            return
        end function hasAttributes


        pure function hasText(this) result(has)
!-------^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
        !*  Returns true if this xml wrapper has any text.
            type(NBAX),intent(in)           ::      this
            logical     ::      has
            has = this%nTextLines > 0
            return
        end function hasText



!-------



        pure function getNAttributes(this) result(ii)
!-------^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
        !*  Returns the number of attributes.
            type(NBAX),intent(in)           ::      this
            integer     ::      ii
            ii = this%nAttributes
            return
        end function getNAttributes

        pure function getNChildren0(this) result(ii)
!-------^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
        !*  Returns the number of children.
            type(NBAX),intent(in)           ::      this
            integer     ::      ii
            ii = this%nChildren
            return
        end function getNChildren0


        function getNChildren1(this,aa) result(ii)
!-------^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
        !*  Returns the number of children called aa
            type(NBAX),intent(in)           ::      this
            character(len=*),intent(in)     ::      aa
            integer                         ::      ii,jj
            ii = 0
            do jj = 1,this%nChildren
                if (this%children(jj) == aa) ii = ii + 1
            end do
            return
        end function getNChildren1


        pure function getNText(this) result(ii)
!-------^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
        !*  Returns the number of lines of text.
            type(NBAX),intent(in)           ::      this
            integer     ::      ii
            ii = this%nTextLines
            return
        end function getNText

!-------

    !---

        pure function NBAXequalsA(this,aa) result(is)
!-------^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
        !*  Equality operator tests whether the name of this xml is aa.
            character(len=*),intent(in)     ::      aa
            type(NBAX),intent(in)           ::      this
            logical             ::      is
            is = (trim(this%name) == trim(aa))
            return
        end function NBAXequalsA

!-------

        subroutine getChildNBAX0(this,ii,m)
!-------^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
        !*  Returns the iith child.
            type(NBAX),intent(in)           ::      this
            integer,intent(in)              ::      ii
            type(NBAX),pointer              ::      m
            if (ii>this%nChildren) then
                nullify(m)
            else
                m => this%children(ii)
            end if
            return
        end subroutine getChildNBAX0

        subroutine getChildNBAX1(this,aa,m,ok)
!-------^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
        !*  Looks for a child with the given name = aa.
            type(NBAX),intent(in)           ::      this
            character(len=*),intent(in)     ::      aa
            type(NBAX),pointer              ::      m
            logical,intent(out)             ::      ok
            integer      ::      jj
            ok = .false.
            nullify(m)
            do jj = 1,this%nChildren
                if (this%children(jj) == aa) then
                    m => this%children(jj)
                    ok = .true.
                    return
                end if
            end do
            return
        end subroutine getChildNBAX1

        subroutine getChildNBAX2(this,aa,i,m,ok)
!-------^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
        !*  Looks for the ith child with the given name = aa.
            type(NBAX),intent(in)           ::      this
            character(len=*),intent(in)     ::      aa
            integer,intent(in)              ::      i
            type(NBAX),pointer              ::      m
            logical,intent(out)             ::      ok
            integer      ::      jj,ll
            ok = .false.
            nullify(m)
            ll = 0
            do jj = 1,this%nChildren
                if (this%children(jj) == aa) then
                    ll = ll + 1
                    if (ll == i) then
                        m => this%children(jj)
                        ok = .true.
                        return
                    end if
                end if
            end do
            return
        end subroutine getChildNBAX2


        subroutine getAttribute(this,ii,k)
!-------^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
        !*  Returns attribute number ii.
            type(NBAX),intent(in)           ::      this
            integer,intent(in)              ::      ii
            type(KeyPair),intent(out)       ::      k
            if (ii>this%nAttributes) then
                k = KeyPair_ctor()
                return
            end if
            k = this%attributes(ii)
            return
        end subroutine getAttribute

        subroutine getText(this,ii,aa,ok)
!-------^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
        !*  Returns the iith line of text.
            type(NBAX),intent(in)           ::      this
            integer,intent(in)              ::      ii
            character(len=*),intent(out)    ::      aa
            logical,intent(out)             ::      ok
            ok = .false.
            aa = ""
            if (ii > this%nTextLines) return
            aa = this%textLines(ii)
            ok = .true.
            return
        end subroutine getText


!-------

        subroutine getValueNBAX0(this,key,ii,ok)
!-------^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
        !*  Looks for an attribute with the given key.
        !*  Tries to return an integer representation of the value.
            type(NBAX),intent(in)           ::      this
            character(len=*),intent(in)     ::      key
            integer,intent(inout)           ::      ii
            logical,intent(out)             ::      ok
            integer                         ::      jj
            ok = .false.
            do jj = 1,this%nAttributes
                if (this%attributes(jj) == key) then
                    call parse( getValue(this%attributes(jj)),ii,ok )
                    return
                end if
            end do
            return
        end subroutine getValueNBAX0

        subroutine getValueNBAX1(this,key,xx,ok)
!-------^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
        !*  Looks for an attribute with the given key.
        !*  Tries to return a real(kind=real64) representation of the value.
            type(NBAX),intent(in)           ::      this
            character(len=*),intent(in)     ::      key
            real(kind=real64),intent(inout)              ::      xx
            logical,intent(out)             ::      ok
            integer                         ::      jj
            ok = .false.
            do jj = 1,this%nAttributes
                if (this%attributes(jj) == key) then
                    call parse( getValue(this%attributes(jj)),xx,ok )
                    return
                end if
            end do
            return
        end subroutine getValueNBAX1

        subroutine getValueNBAX2(this,key,aa,ok)
!-------^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
        !*  Looks for an attribute with the given key.
        !*  Tries to return the value.
            type(NBAX),intent(in)           ::      this
            character(len=*),intent(in)     ::      key
            character(len=*),intent(out)    ::      aa
            logical,intent(out)             ::      ok
            integer                         ::      jj
            ok = .false.
            do jj = 1,this%nAttributes
                if (this%attributes(jj) == key) then
                    aa = getValue(this%attributes(jj))
                    ok = .true.
                    return
                end if
            end do
            return
        end subroutine getValueNBAX2

!        subroutine getValueNBAX3(this,key,xx,ok)
!!-------^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
!        !*  Looks for an attribute with the given key.
!        !*  Tries to return a complex representation of the value.
!            type(NBAX),intent(in)           ::      this
!            character(len=*),intent(in)     ::      key
!            complex,intent(inout)           ::      xx
!            logical,intent(out)             ::      ok
!            integer                         ::      jj
!            ok = .false.
!            do jj = 1,this%nAttributes
!                if (this%attributes(jj) == key) then
!                    call parse( getValue(this%attributes(jj)),xx,ok )
!                    return
!                end if
!            end do
!            return
!        end subroutine getValueNBAX3

        subroutine getValueNBAX4(this,key,ll,ok)
!-------^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
        !*  Looks for an attribute with the given key.
        !*  Tries to return a logical representation of the value.
            type(NBAX),intent(in)           ::      this
            character(len=*),intent(in)     ::      key
            logical,intent(inout)           ::      ll
            logical,intent(out)             ::      ok
            integer                         ::      jj
            ok = .false.
            do jj = 1,this%nAttributes
                if (this%attributes(jj) == key) then
                    call parse( getValue(this%attributes(jj)),ll,ok )
                    return
                end if
            end do
            return
        end subroutine getValueNBAX4

    !---

        subroutine changeValueNBAX0(this,key,ii,ok)
!-------^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
        !*  Looks for an attribute with the given key.
        !*  tries to change an integer representation of the value.
            type(NBAX),intent(inout)        ::      this
            character(len=*),intent(in)     ::      key
            integer,intent(in)              ::      ii
            logical,intent(out)             ::      ok
            integer                         ::      jj
            integer                         ::      iiii
            ok = .false.
            do jj = 1,this%nAttributes
                if (this%attributes(jj) == key) then
                    call parse( getValue(this%attributes(jj)),iiii,ok )
                    if (ok) this%attributes(jj) = KeyPair_ctor( trim(getKey(this%attributes(jj))),ii )
                    return
                end if
            end do
            if (.not. ok) call addAttribute(this,key,ii)
            return
        end subroutine changeValueNBAX0

        subroutine changeValueNBAX1(this,key,xx,ok)
!-------^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
        !*  Looks for an attribute with the given key.
        !*  tries to change a real(kind=real64) representation of the value.
            type(NBAX),intent(inout)        ::      this
            character(len=*),intent(in)     ::      key
            real(kind=real64),intent(in)                 ::      xx
            logical,intent(out)             ::      ok
            integer                         ::      jj
            real(kind=real64)                            ::      xxxx
            ok = .false.
            do jj = 1,this%nAttributes
                if (this%attributes(jj) == key) then
                    call parse( getValue(this%attributes(jj)),xxxx,ok )
                    if (ok) this%attributes(jj) = KeyPair_ctor( trim(getKey(this%attributes(jj))),xx )
                    return
                end if
            end do
            if (.not. ok) call addAttribute(this,key,xx)
            return
        end subroutine changeValueNBAX1

        subroutine changeValueNBAX2(this,key,aa,ok)
!-------^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
        !*  Looks for an attribute with the given key.
        !*  tries to change the value.
            type(NBAX),intent(inout)        ::      this
            character(len=*),intent(in)     ::      key
            character(len=*),intent(in)     ::      aa
            logical,intent(out)             ::      ok
            integer                         ::      jj
            ok = .false.
            do jj = 1,this%nAttributes
                if (this%attributes(jj) == key) then
                    this%attributes(jj) = KeyPair_ctor( trim(getKey(this%attributes(jj))),trim(aa) )
                    ok = .true.
                    return
                end if
            end do
            if (.not. ok) call addAttribute(this,key,aa)
            return
        end subroutine changeValueNBAX2

!        subroutine changeValueNBAX3(this,key,xx,ok)
!!-------^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
!        !*  Looks for an attribute with the given key.
!        !*  tries to change a complex representation of the value.
!            type(NBAX),intent(inout)        ::      this
!            character(len=*),intent(in)     ::      key
!            complex,intent(in)              ::      xx
!            logical,intent(out)             ::      ok
!            integer                         ::      jj
!            complex                         ::      xxxx
!            ok = .false.
!            do jj = 1,this%nAttributes
!                if (this%attributes(jj) == key) then
!                    call parse( getValue(this%attributes(jj)),xxxx,ok )
!                    if (ok) this%attributes(jj) = KeyPair_ctor( trim(getKey(this%attributes(jj))),xxxx )
!                    return
!                end if
!            end do
!            if (.not. ok) call addAttribute(this,key,xx)
!            return
!        end subroutine changeValueNBAX3

        subroutine changeValueNBAX4(this,key,ll,ok)
!-------^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
        !*  Looks for an attribute with the given key.
        !*  tries to change a logical representation of the value.
            type(NBAX),intent(inout)        ::      this
            character(len=*),intent(in)     ::      key
            logical,intent(in)              ::      ll
            logical,intent(out)             ::      ok
            integer                         ::      jj
            logical                         ::      llll
            ok = .false.
            do jj = 1,this%nAttributes
                if (this%attributes(jj) == key) then
                    call parse( getValue(this%attributes(jj)),llll,ok )
                    if (ok) then
                        if (ll) then
                            this%attributes(jj) = KeyPair_ctor( trim(getKey(this%attributes(jj))),"true" )
                        else
                            this%attributes(jj) = KeyPair_ctor( trim(getKey(this%attributes(jj))),"false" )
                        end if
                    end if
                    return
                end if
            end do
            if (.not. ok) call addAttribute(this,key,ll)
            return
        end subroutine changeValueNBAX4



!-------

        subroutine inputNBAX0(this,filename,ok)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ 
            type(NBAX),intent(inout)        ::      this
            character(len=256),intent(in)   ::      filename
            logical,intent(out)             ::      ok
            integer             ::      uu
            this = NBAX_ctor()
            inquire(file=trim(filename),exist=ok)
            if (.not. ok) then
                print *,"NBAX3::inputNBAX0() error - file """//trim(filename)//""" not found"
                stop
            end if
            uu = 400 ; call findGoodUnit_NBAX(uu)
            open(unit=uu,file=trim(filename),action="read")
            call inputNBAX(this,uu,ok)
            close(unit=uu)
            return
        end subroutine inputNBAX0
            
            

        recursive subroutine inputNBAX(this,uu,loadok,oo)
!-------^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
        !*  Attempts to read xml from unit uu.
            type(NBAX),intent(inout)        ::      this
            integer,intent(in)              ::      uu
            logical,intent(inout)           ::      loadok
            logical,intent(in),optional     ::      oo
            logical         ::      opened,ok
            character(len=2048)   ::      ss,sss
            integer         ::      ioerr,lt,uuu
            type(NBAX)            ::      xml_tmp
            type(StringTokenizer)           ::      st
            type(KeyPair)   ::      kp
            if (present(oo)) then
                opened = oo
            else
                opened = .false.
                loadok = .false.
            end if
            do
                ss =""
                read(unit=uu,fmt='(a2048)',iostat = ioerr) ss
                if ( ioerr/=0 ) return
                call cutComments(ss)
                lt = lineType(ss)
                !write (*,fmt='(a,i4,l4,a,l4)')"       read line",lt,opened,""//trim(ss)
                st = StringTokenizer_ctor(ss,"< >"//TAB_CHARACTER)
                select case(lt)
                    case (NBAX_OPEN:NBAX_OPENCLOSE)
                        call nextToken(st,sss)
                        if (opened) then
                            xml_tmp = NBAX_ctor(trim(sss))      !   this is the line which allocates memory for xml_tmp
                            sss = getRemainingTokens(st)
                            ok = .true.
                            do
                                call extractKeypair(sss,kp,ok)
                                if (ok) then
                                    call addAttribute(xml_tmp,kp)
                                     ! print *,'addAttribute ',trim(toString(kp))
                                else
                                    exit
                                end if
                            end do
                            if (lt == NBAX_OPEN) then
                                ! print *,trim(this%name)," opened, opening child",trim(xml_tmp%name)
                                call input(xml_tmp,uu,loadok,.true.)
                            end if
                            if (xml_tmp=="xi:include") then
                                call getAttributeValue(xml_tmp,"href",ss,ok)
                                call replaceEnv(ss,sss)
                                if (ok) then
                                    ! print *,trim(this%name)," - special instruction insert file"//trim(sss)
                                    if (.not. hasChildren(xml_tmp)) then
                                        inquire( file=trim(sss),exist=ok )
                                        if (ok) then
!                                             call delete(xml_tmp%attributes,.true.)
!                                             call delete(xml_tmp%text,.true.)
                                            uuu = uu + 1
                                            call findGoodUnit_NBAX(uuu)
                                            open(unit=uuu,file=trim(sss),action="read")
                                                call inputNBAX( xml_tmp,uuu,loadok,.false.)
                                            close(unit=uuu)
                                            ! print *,"   adding ",trim(xml_tmp%name)," to parent ",trim(this%name)
                                            call addChild(this,xml_tmp)
                                        else
                                            write(unit=0,fmt='(a)')"NBAX::inputNBAX error - malformed xi:include, file not found " &
                                                //trim(sss)
                                            call delete(xml_tmp)
                                        end if
                                    else
                                        write(unit=0,fmt='(a)')"NBAX::inputNBAX error - malformed xi:include, children found?"
                                        call delete(xml_tmp)
                                    end if
                                else
                                    write(unit=0,fmt='(a)')"NBAX::inputNBAX error - malformed xi:include, missing href"
                                    call delete(xml_tmp)
                                end if
                            else
                                ! print *,"adding",trim(xml_tmp%name)," to parent",trim(this%name)
                                call addChild(this,xml_tmp)
                            end if
                        else
                            this%name = sss
                            sss = getRemainingTokens(st)
                            ok = .true.
                            do
                                call extractKeypair(sss,kp,ok)
                                if (ok) then
                                    call addAttribute(this,kp)
                                    !  print *,'addAttribute ',trim(toString(kp))
                                else
                                    exit
                                end if
                            end do
                            if (lt == NBAX_OPEN) then
                                ! print *,"reading xml",trim(this%name)
                                opened = .true.
                            end if
                        end if
                        loadok = .true.
                    case (NBAX_OPENCOMMENT)
                        do
                            ss =""
                            read(unit=uu,fmt='(a)',iostat = ioerr) ss
                            if ( ioerr/=0 ) return
                            lt = lineType(ss)
                            if (lt==NBAX_CLOSECOMMENT) exit
                        end do
                    case (NBAX_OPENDOCTYPE)
                        do
                            ss =""
                            read(unit=uu,fmt='(a)',iostat = ioerr) ss
                            if ( ioerr/=0 ) return
                            lt = lineType(ss)
                            if (lt==NBAX_CLOSEDOCTYPE) exit
                        end do
                    case (NBAX_TEXT)
                        if (opened) then
                            sss = adjustl(getRemainingTokens(st))
                            if (sss /="") then
                                !print *,'addText ',trim(sss)
                                call addText(this,sss)
                            end if
                        end if
                    case (NBAX_CLOSE)
                        loadok = .true.
                        ! print *,"xml closed",trim(this%name)
                        opened = .false.
                        ! print *,"NBAX_CLOSE"
                        return
                end select
            end do
            return
        end subroutine inputNBAX

        subroutine cutComments(aa)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      remove internal comments from the line
    !*       eg <Herring attrib="1" <!-- a herring --> /> =  <Herring attrib="1"/>
            character(len=*),intent(inout)  ::      aa
            integer         ::      lt,ioc1,ioc2
            do
                lt = lineType(aa)
                if (lt == NBAX_OPENCLOSECOMMENT) then
                !   there is an internal comment to remove
                    ioc1 = index(aa,"<!--")
                    if (ioc1 <= 0) return
                    ioc2 = index(aa,"-->")
                    if( ioc1 > 1) then
                        if (ioc2+3 < len(aa)) then
                            aa = aa(1:ioc1-1)//aa(ioc2+3:)
                        else
                            aa = aa(1:ioc1-1)
                        end if
                    else
                        if (ioc2 < len(aa)) then
                            aa = aa(ioc2:)
                        else
                            aa =""
                        end if
                    end if
                else
                    return
                end if
            end do
            return
        end subroutine cutComments




        pure function lineType(aa) result(ii)
!-------^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
            character(len=*),intent(in)     ::      aa
            integer                         ::      ii
            integer         ::  io
            ii = NBAX_TEXT

            if (index(aa,"<!--") > 0) then
                if (index(aa,"-->") > 0) then
                    ii = NBAX_OPENCLOSECOMMENT
                else
                    ii = NBAX_OPENCOMMENT
                end if
                return
            else
                if (index(aa,"-->") > 0) then
                    ii = NBAX_CLOSECOMMENT
                    return
                end if
            end if

            io = index(aa,"<?xml")
            if (io > 0) then
                ii = NBAX_XMLSPEC
                return
            end if

            io = index(aa,"<!DOCTYPE")
            if (io > 0) then
                if (index(aa,">") > 0) then
                    ii = NBAX_OPENCLOSEDOCTYPE
                else
                    ii = NBAX_OPENDOCTYPE
                end if
                return
            end if

            io = index(aa,"]>")
            if (io > 0) then
                ii = NBAX_CLOSEDOCTYPE
                return
            end if


            io = index(aa,"<")
            if ( io > 0 ) then
                if (index(aa,"/>") > 0) then
                    ii = NBAX_OPENCLOSE
                else if (index(aa,"</") == io) then
                    ii = NBAX_CLOSE
                else
                    ii = NBAX_OPEN
                end if
            end if
            return
        end function lineType


        subroutine replaceEnv(strin,strout)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*
            character(len=*),intent(in)         ::      strin
            character(len=*),intent(out)        ::      strout
            integer             ::      i1,i2,i3,jj
            character(len=2048)  ::      env
            strout = strin
            jj = 1
            do
                i1 = index( strout(jj:),"${" )
                i2 = 0 ; i3 = len_trim(strout)
                if (i1 > 0) then
                    i2 = i1 + jj + 1
                    i3 = index( strout(i2:),"}" ) + i2 - 2
                else
                    i1 = index( strout(jj:),"%" )
                    i2 = i1 + jj
                    if (i1 > 0) i3 = index( strout(i2:),"%" ) + i2 - 2
                end if
                if (i1 == 0) exit
!                call getenv_f( strout( i2:i3 ),env )
                call get_environment_variable( strout( i2:i3 ),env )
                env = trim(env)//strout(i3+2:)
                if (i1 > 1) then
                    strout = strout(1:i1+jj-2) // trim(env)
                else
                    strout = trim(env)
                end if
                jj = i2 + 1
            end do

            return
        end subroutine replaceEnv


        subroutine findGoodUnit_NBAX(uu)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      find a good (unopened) unit to read/write
            integer,intent(inout)       ::      uu
            integer             ::      i
            logical             ::      lod
            ! find a good unit number ( i )
            i = uu
            do
                inquire (opened = lod,unit = i)
                if (.not. lod) exit
                i = i + 1
            end do
            uu = i
            return
        end subroutine findGoodUnit_NBAX

!-------

        pure subroutine setNameNBAX(this,a)
!-------^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
        !*  Sets the name of this xml wrapper to a.
            character(len=*),intent(in)     ::      a
            type(NBAX),intent(inout)        ::      this
            this%name = a
            return
        end subroutine setNameNBAX

        pure function getNameNBAX(this) result(a)
!-------^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
        !*  Returns the name of this xml wrapper.
            type(NBAX),intent(in)           ::      this
            character(len=len(this%name))   ::      a
            a = this%name
            return
        end function getNameNBAX

!-------

        subroutine cutChildren(this)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^
            type(NBAX),intent(inout)        ::      this
            integer                             ::  ii,nn
            integer,dimension(:),allocatable    ::  cc
            do ii = 1,this%nChildren
                call delete(this%children(ii))
            end do
            this%nChildren = 0
            if (this%N_CONTENT>0) then
                allocate(cc(this%N_CONTENT))
                nn = 0
                cc = 0
                do ii = 1,this%nContent
                    if (this%content_type(ii)<=0) then
                        nn = nn + 1
                        cc(nn) = this%content_type(ii)
                    end if
                end do
                this%nContent = nn
                this%content_type = 0
                this%content_type(1:nn) = cc(1:nn)
                deallocate(cc)
            end if
            return
        end subroutine cutChildren

        subroutine cutText(this)
    !---^^^^^^^^^^^^^^^^^^^^^^^^
            type(NBAX),intent(inout)        ::      this
            integer                             ::  ii,nn
            integer,dimension(:),allocatable    ::  cc
            this%nTextLines = 0
            if (this%N_CONTENT>0) then
                allocate(cc(this%N_CONTENT))
                nn = 0
                cc = 0
                do ii = 1,this%nContent
                    if (this%content_type(ii)>=0) then
                        nn = nn + 1
                        cc(nn) = this%content_type(ii)
                    end if
                end do
                this%nContent = nn
                this%content_type = 0
                this%content_type(1:nn) = cc(1:nn)
                deallocate(cc)
            end if
            return
        end subroutine cutText

        subroutine cutAttributes(this)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
            type(NBAX),intent(inout)            ::      this
            integer                             ::      ii
            do ii = 1,this%nAttributes
                call delete(this%attributes(ii))
            end do
            this%nAttributes = 0
            return
        end subroutine cutAttributes


        pure function KeyPairNull() result(this)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
        !*      Produces an empty keypair.
            type(KeyPair)                   ::      this
            this%key = ""
            this%value = KP_FALSE
            return
        end function KeyPairNull

        pure function KeyPairConstructor1(k,v) result(this)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
        !*      Produces a standard keypair k = v.
            character(len=*),intent(in)     ::      k,v
            type(KeyPair)                   ::      this
            integer         ::      i
            i = min( len_trim(k),KEYPAIR_MAX_STRING )
            this%key = k(1:i)
            i = min( len_trim(v),KEYPAIR_MAX_STRING )
            this%value = v(1:i)
            return
        end function KeyPairConstructor1

        pure function KeyPairConstructor2(k) result(this)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
        !*      Produces a simple keypair k = "true".
            character(len=*),intent(in)     ::      k
            type(KeyPair)                   ::      this
            integer         ::      i
            i = min( len_trim(k),KEYPAIR_MAX_STRING )
            this%key = k(1:i)
            this%value = KP_TRUE
            return
        end function KeyPairConstructor2


        function KeyPairConstructor3(k,ii) result(this)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
        !*  Utility constructor produces k = <text representation of ii>.
            character(len=*),intent(in)     ::      k
            integer,intent(in)              ::      ii
            type(KeyPair)                   ::      this
            integer         ::      i
            i = min( len_trim(k),KEYPAIR_MAX_STRING )
            this%key = k(1:i)
            write(this%value,fmt=*) ii
            this%value = adjustl(this%value)
            return
        end function KeyPairConstructor3


        function KeyPairConstructor4(k,xx) result(this)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
        !*  Utility constructor produces k = <text representation of xx>.
            character(len=*),intent(in)     ::      k
            real(kind=real64),intent(in)                 ::      xx
            type(KeyPair)                   ::      this
            integer         ::      i
            i = min( len_trim(k),KEYPAIR_MAX_STRING )
            this%key = k(1:i)
            write(this%value,fmt=*) xx
            this%value = adjustl(this%value)
            return
        end function KeyPairConstructor4

        function KeyPairConstructor5(k,xx) result(this)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
        !*  Utility constructor produces k = <text representation of xx>.
            character(len=*),intent(in)     ::      k
            complex,intent(in)              ::      xx
            type(KeyPair)                   ::      this
            integer         ::      i
            i = min( len_trim(k),KEYPAIR_MAX_STRING )
            this%key = k(1:i)
            write(this%value,fmt=*) xx
            this%value = adjustl(this%value)
            return
        end function KeyPairConstructor5

!-------

        subroutine deleteKeyPair(this)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
        !*  No dynamic memory- this does nothing.
            type(KeyPair),intent(inout)     ::  this
            this%key = ""
            this%value = ""
            return
        end subroutine deleteKeyPair

!-------

        pure function getKey1(this) result(k)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
        !*  Returns the "key" half of the pair.
            type(KeyPair),intent(in)            ::      this
            character(len=KEYPAIR_MAX_STRING)   ::      k
            k = this%key
            return
        end function getKey1

        pure function getValue1(this) result(v)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
        !*  Returns the "value" half of the pair.
            type(KeyPair),intent(in)            ::      this
            character(len=KEYPAIR_MAX_STRING)   ::      v
            v = this%value
            return
        end function getValue1

        pure function getValue2(kp,k) result(v)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
        !*  Searches through the array kp for one with key = k
        !*  Returns the "value" half of the pair.
            type(KeyPair),dimension(:),intent(in)            ::      kp
            character(len=*),intent(in)         ::      k
            character(len=KEYPAIR_MAX_STRING)   ::      v
            integer             ::      i
            do i = 1,size(kp)
                if (trim(k) == trim(kp(i)%key)) then
                    v = kp(i)%value
                    return
                end if
            end do
            v = KP_FALSE
            return
        end function getValue2

!-------

        pure function toStringKP(this) result(s)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
        !*  Returns a standard string representation of a keypair.
        !*  In the form "<key>"="<value>"
            type(KeyPair),intent(in)                ::      this
            character(len=2*KEYPAIR_MAX_STRING+1)   ::      s
            s = ""
            s = trim(this%key)//'="'//trim(this%value)//'"'
            return
        end function toStringKP

!-------

        subroutine extractKeyPair(ss,this,ok)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
        !*  Try to extract a keypair from a standard input character string.
        !*  There is a little flexibility, as this method uses StringTokenizers
        !*  but try to input in the standard form "<key>"="<value>".
        !*  On output ss has the remainder of the input string, so this subroutine
        !*  can be used to grab all the keypairs from a single input line.
            character(len=*),intent(inout)      ::      ss
            type(KeyPair),intent(out)           ::      this
            logical,intent(out)                 ::      ok
            type(StringTokenizer)               ::      st
            character(len=STRINGTOKENIZER_MAX_STRING)   ::      aa,bb,sss

            ok = .false.
            st = StringTokenizer_ctor(ss,' ="')
            aa = ""
            bb = ""
            call nextToken(st,aa)
!            print *,"extractKeyPair::nextToken ",trim(aa)
            sss = getRemainingTokens(st)
!             print *,"extractKeyPair::getRemainingTokens ",trim(sss)
            if (index(ss,"=")>0) then
                if (index(sss,'"')>0) then
                    st = StringTokenizer_ctor(sss,'">')
                else
                    st = StringTokenizer_ctor(sss," ,;>")
                end if
                call nextToken(st,bb)
!                 print *,"extractKeyPair::nextToken ",trim(bb)
                call trimPrefix(st,", ;")
                ss = getRemainingTokens(st)
!                 print *,"extractKeyPair::getRemainingTokens ",trim(ss)
                ok = .true.
            end if
            this%key = trim(aa)
            this%value = trim(bb)
            return
        end subroutine extractKeyPair

        function KeyPairfromString(s) result(this)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
        !*  This is a simple minded extractor.
            character(len=*),intent(in)         ::      s
            type(KeyPair)                       ::      this
            type(StringTokenizer)               ::      st
            character(len=KEYPAIR_MAX_STRING)   ::      k,v
            st = StringTokenizer_ctor(s,'"= ')
            k = ""
            v = KP_TRUE
            if (hasMoreTokens(st)) call nextToken(st,k)
            if (hasMoreTokens(st)) call nextToken(st,v)
            this = KeyPairConstructor1(k,v)
            return
        end function KeyPairfromString


        function KeyfromString(s) result(k)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
        !*  Returns the "key" half from an appropriate string.
            character(len=*),intent(in)         ::      s
            type(StringTokenizer)               ::      st
            character(len=KEYPAIR_MAX_STRING)   ::      k
            st = StringTokenizer_ctor(s,'"= ')
            k = ""
            if (hasMoreTokens(st)) call nextToken(st,k)
            return
        end function KeyfromString

!-------

        function ValuefromString(s) result(v)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
        !*  Returns the "value" half from an appropriate string.
            character(len=*),intent(in)         ::      s
            type(StringTokenizer)               ::      st
            character(len=KEYPAIR_MAX_STRING)   ::      k,v
            st = StringTokenizer_ctor(s,'"= ')
            k = ""
            v = KP_TRUE
            if (hasMoreTokens(st)) call nextToken(st,k)
            if (hasMoreTokens(st)) call nextToken(st,v)
            return
        end function ValuefromString

!-------

        elemental subroutine assignStringFromKP(s,this)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
        !*  Assignment operator s = this.
            character(len=2*KEYPAIR_MAX_STRING+1),intent(out)       ::      s
            type(keyPair),intent(in)                            ::      this
            s = toStringKP(this)
            return
        end subroutine assignStringFromKP

!-------

        pure function isDefinedKP(this) result(is)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
        !*  Returns true if the value half is not false.
            type(KeyPair),pointer       ::      this
            logical                     ::      is
            is = associated(this)
            if (is) then
                is = .not.( trim(this%value)==trim(KP_FALSE) )
            end if
            return
        end function isDefinedKP

        elemental function KPequalsKP(kp1,kp2) result(is)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
        !*  Equality operator kp1==kp2.
            type(KeyPair),intent(in)        ::      kp1,kp2
            logical                         ::      is
            is = ( ( trim(kp1%key) == trim(kp2%key) ).and.          &
                    ( trim(kp1%value) == trim(kp2%value) ) )
            return
        end function KPequalsKP

        elemental function KPequalsA(this,a) result(is)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
        !*  Equlity operator checks if the key is a.
            type(KeyPair),intent(in)        ::      this
            character(len=*),intent(in)     ::      a
            logical                         ::      is
            is = ( trim(this%key) == trim(a) )
            return
        end function KPequalsA

!-------

        subroutine outputKP(this,u)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^
        !*  Dumps a simple representation to unit u.
            type(KeyPair),intent(in)        ::      this
            integer,intent(in),optional     ::      u
            integer             ::      uu
            uu = 6
            if (present(u)) uu = u
            write(unit=uu,fmt='(a)') trim( toString(this) )
            return
        end subroutine outputKP


    end module Lib_NBAX




!
!
! ! !------------------------------------------------------------------------------
!
!     program testNBAX
! !---^^^^^^^^^^^^^^^^
!
!         use Lib_NBAX
! !         use NBAX_KeyPairs
!         implicit none
!
!         type(NBAX),pointer      ::  m1,m2,m3,m4,m5 , mm
!         logical                 ::  ok
!         real(kind=real64)                    ::  xx
!
! !-------
!         print *,""
!         print *,"Test XML construction"
!         allocate(m1,m2,m3,m4,m5)
!         m1 = NBAX_ctor("test_outer1")
!         m2 = NBAX_ctor("test_inner2")
!         m3 = NBAX_ctor("test_inner3")
!         m4 = NBAX_ctor("test_inner4")
!         m5 = NBAX_ctor("xi:include")
!
!         print *,""
!         print *,"addAttribute"
!         call addAttribute( m1,"outermost","true" )
!         call addAttribute( m2,"real(kind=real64)_value",0.123 )
!         call addAttribute( m3,"attr1",1 )
!         call addAttribute( m3,"attr2","2" )
!         call addAttribute( m5,"href","insert.xml" )
!
!         print *,""
!         print *,"addText"
!         call addText( m1,"text add demo" )
!         call addText( m4,"Mary had a little lamb" )
!         call addText( m4,"her fleece was white as snow" )
!
!         print *,""
!         print *,"addChild"
!         call addChild(m2,m4)
!         call addChild(m1,m2)
!         call addChild(m1,m3)
!         call addChild(m1,m5)
!
!         call addText( m1,"more text content demo" )
!
!
! ! !-------
!         print *,""
!         print *,"Test XML output"
!         call output(m1)
!
! !-------
!         print *,""
!         print *,"Test XML data extraction"
!         call getChild(m1,"test_inner2",mm,ok)
!         if (ok) then
!             call getAttributeValue(mm,"real(kind=real64)_value",xx,ok)
!             if (ok) print *,"test_inner2::real(kind=real64)_value =",xx
!         end if
!
! !-------
!         print *,""
!         print *,"Test XML write to disk"
!         open(unit=20,file="test.xml",action="write",status="unknown")
!             call output(m1,20,.true.)
!         close(unit=20)
!
! !-------
!         print *,""
!         print *,"Test XML delete"
!         call delete(m1)
!
!         print *,""
!         print *,"Test XML output"
!         call output(m1)
!
!
! !-------
!         print *,""
!         print *,"Test XML read from disk"
!         open(unit=20,file="test.xml",action="read",status="old")
!             call input(m1,20,ok)
!         close(unit=20)
!         print *,""
!         print *,"Test XML output"
!         if (ok) call output(m1)
!
! !-------
!         print *,""
!         print *,"Test cut children from inner2, cut attributes from outer1"
!         call getChild(m1,"test_inner2",mm,ok)
!         call cutChildren( mm )
!         call cutAttributes( m1 )
!         print *,""
!         print *,"Test XML output"
!         call output(m1)
!
!         print *,""
!          print *,"Tests complete"
!          call delete(m1)
!          deallocate(m1,m2,m3,m4,m5)
!
!     !---
!         print *,""
!         print *,"done"
!         print *,""
!     end program testNBAX
