
    module NBAX_StringTokenizers
!---^^^^^^^^^^^^^^^^^^^^^^^^^^^
!
!*      A simple implementation of a StringTokenizer
!*      as an analogue to the one found in Java.
!*      Also contains some parsing routines.
!
!*      author      Daniel Mason
!*      version     1.0
!*      revision    May 2005
!

!*      Usage:
!*      A character string "Hello, world. How are you today?"
!*      can be parsed into its constituent tokens if the delimiters
!*      between tokens are given. For instance I might choose space and
!*      punctuation. The delimiters are then del = " .,?!:;".
!*
!*      character(len=256)  :: string = "Hello, world. How are you today?"
!*      character(len=16)   :: del = " .,?!:;"
!*      character(len=64)   :: token
!*      type(StringTokenizer) :: st
!*      st = StringTokenizer_ctor(string,del)
!*      do
!*          if (.not. hasMoreTokens(st)) exit
!*          call nextToken(st,token)
!*          print *,token
!*      end do
!* 
!*      produces
!*
!*      Hello
!*      world
!*      how
!*      are
!*      you
!*      today
!*
!*
        use iso_fortran_env
        implicit none
        private

        public              ::      StringTokenizer_ctor            !   generates an object with information about how to break up the string
        public              ::      nextToken                       !   returns the next token in a string
        public              ::      hasMoreTokens                   !   true if there are any tokens left
        public              ::      getRemainingTokens              !   returns the rest of the strin
        public              ::      trimPrefix                      !   removes any junk made up of delimiters from the start of the string
        public              ::      parse                           !   converts string to number 
        public              ::      adjustl                         !   removes delimiters from start of string
        public              ::      ftoa,utoa,itoa                  !   converts number to string
        public              ::      convertUpperCase
        public              ::      convertLowerCase

        integer,public,parameter            ::      STRINGTOKENIZER_MAX_STRING = 1024
        integer,public,parameter            ::      STRINGTOKENIZER_MAX_TOKENS = 8

        character(len=1),public,parameter   ::      TAB_CHARACTER = char(9)
        character(len=1),public,parameter   ::      CR_CHARACTER = char(13)
        character(len=1),public,parameter   ::      NULL_CHARACTER = char(0)


!-------    TYPE DEFINITIONS

        type,public         ::  StringTokenizer
            private
            character(len=STRINGTOKENIZER_MAX_STRING)   ::      buffer
            integer                                     ::      nChars
            character(len=STRINGTOKENIZER_MAX_TOKENS)   ::      chr
        end type

!-------    INTERFACE

        interface   StringTokenizer_ctor
            module procedure    StringTokenizerConstructor
        end interface


        interface   Parse
            module procedure            ParseIntArray1
            module procedure            ParseRealArray1
            module procedure            ParseLogicalArray
            module procedure            ParseCharacterArray
            module procedure            ParseInt
            module procedure            ParseReal
!            module procedure            ParseComplex
            module procedure            ParseLogical
        end interface

        interface   adjustl
            module procedure            adjustl0
        end interface
        
        interface   ftoa
            module procedure            ftoa0       !   float to string
            module procedure            fatoa1      !   float array to string, fixed width
            module procedure            fatoa1a     !   float rank 2 array to string, fixed width
            module procedure            fatoa2      !   float array to string, variable width
            module procedure            fatoa3      !   float array to string, variable width     
        end interface

        
        interface   itoa
            module procedure            itoa0       !   int to string
            module procedure            iatoa1      !   int array to string 
        end interface


    contains
!---^^^^^^^^

        pure function StringTokenizerConstructor(s1,s2) result(st)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
        !*  Constructs a string tokenizer from the string s1 with
        !*  delimiters s2.
            character(len=*),intent(in)         ::      s1
            character(len=*),intent(in)         ::      s2
            type(StringTokenizer)               ::      st
            integer             ::      i
            i = min(STRINGTOKENIZER_MAX_STRING,len_trim(s1))
            st%buffer = ""
            st%buffer(1:i) = s1
            st%nChars = min(STRINGTOKENIZER_MAX_TOKENS,len(s2))
            do i = 1,st%nChars
                st%chr(i:i) = s2(i:i)
            end do
            return
        end function StringTokenizerConstructor

!-------

        pure function hasMoreTokens(st) result(has)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
        !*  Returns true if st has more tokens.
            type(StringTokenizer),intent(in)    ::      st
            logical                             ::      has
            has = (len_trim(st%buffer)>0)
            return
        end function hasMoreTokens

        pure function getRemainingTokens(st) result(ss)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
        !*  Returns a character string containing all remaining tokens, discarding the first i
            type(StringTokenizer),intent(in)    ::      st
            character(len=STRINGTOKENIZER_MAX_STRING)   ::  ss
            ss = trim(st%buffer)
            return
        end function getRemainingTokens


        function adjustl0(st) result(ss)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
        !*  Returns a character string containing all remaining tokens, discarding the first i
            type(StringTokenizer),intent(inout)         ::      st
            character(len=STRINGTOKENIZER_MAX_STRING)   ::      ss
            integer             ::      ii,jj
            logical             ::      done
            do ii = 1,len_trim(st%buffer)
                done = .true.
                do jj = 1,st%nChars
                    if (st%buffer(ii:ii) == st%chr(jj:jj) ) done = .false.
                end do
                if (done) then
                    st%buffer = st%buffer(ii:)
                    exit
                end if
            end do
            ss = trim(st%buffer)
            return
        end function adjustl0


!-------

        pure function getIndex(st,s) result(jj)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

            type(StringTokenizer),intent(in)    ::      st
            character(len=*),intent(in)         ::      s
            character(len=1)                    ::      c
            integer     ::      jj
            integer     ::      i,kk
            jj = len(s)+1
            do i = 1,st%nChars
                c = st%chr(i:i)
                kk = index( s,c )
                if ( kk>0 ) then
                    jj = min(kk,jj)
                end if
!                print *,"test "//s//" with "//c//" ",kk,jj
            end do
            return
        end function getIndex

        recursive subroutine nextToken(st,t)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
        !*  Returns t as the next token present in st.
            type(StringTokenizer),intent(inout)             ::      st
            character(len=*),intent(inout)                  ::      t
            integer         ::      jj,kk
            integer         ::      i
            t = ""
            kk = len_trim(st%buffer)
            if (kk==0) return
            jj = getIndex(st,st%buffer(1:kk))
            i = min(len(t),jj-1)
            t = st%buffer(1:i)
            st%buffer(1:kk) = st%buffer(jj+1:kk)//repeat(" ",jj)
            if (i==0) call nextToken(st,t)
            return
        end subroutine nextToken

        recursive subroutine trimPrefix(st,delim)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
        !*  Trims the start of the string held in st of all the
        !*  given delimiters.
            type(StringTokenizer),intent(inout)             ::      st
            character(len=*),intent(in)                     ::      delim
            integer         ::      jj,kk
            kk = len_trim(st%buffer)
            if (kk==0) return
            do jj = 1,len(delim)
                if (st%buffer(1:1) == delim(jj:jj)) then
                    st%buffer(1:kk) = st%buffer(2:kk)//" "
                    call trimPrefix(st,delim)
                end if
            end do
            return
        end subroutine trimPrefix



        pure function convertUpperCase(text_in) result(text_out)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      convert string to all upper case
            character(len=*),intent(in)                         ::      text_in
            character(len=len(text_in))                         ::      text_out
            integer             ::      ii,ascii
            integer,parameter   ::      ASCII_a = iachar("a")
            integer,parameter   ::      ASCII_z = iachar("z")
            integer,parameter   ::      ASCII_CASE_OFFSET = iachar("A") - iachar("a")
            text_out = text_in
            do ii = 1,len_trim(text_in)
                ascii = iachar( text_in(ii:ii) )        !   convert character i to its ascii representation
                if ( (ascii>=ASCII_a).and.(ascii<=ASCII_z) ) text_out(ii:ii) = achar( ascii + ASCII_CASE_OFFSET )
            end do
            return
        end function convertUpperCase


        pure function convertLowerCase(text_in) result(text_out)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      convert string to all lower case
            character(len=*),intent(in)                         ::      text_in
            character(len=len(text_in))                         ::      text_out
            integer             ::      ii,ascii
            integer,parameter   ::      ASCII_A = iachar("A")
            integer,parameter   ::      ASCII_Z = iachar("Z")
            integer,parameter   ::      ASCII_CASE_OFFSET = iachar("a") - iachar("A")
            text_out = text_in
            do ii = 1,len_trim(text_in)
                ascii = iachar( text_in(ii:ii) )        !   convert character i to its ascii representation
                if ( (ascii>=ASCII_A).and.(ascii<=ASCII_Z) ) text_out(ii:ii) = achar( ascii + ASCII_CASE_OFFSET )
            end do
            return
        end function convertLowerCase
!-------
! 
!         recursive subroutine ParseRealArray(dummy,ra,n,m)
!     !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
!     !*      Parse the character string dummy into an array of n real(kind=real64)s.
!     !*      m is used internally.
!             character(len=*),intent(in)         ::      dummy
!             real(kind=real64),dimension(:),intent(inout)     ::      ra
!             integer,intent(out)                 ::      n
!             integer,intent(in),optional         ::      m       ! number of filled entries in ia
!             character(len=len(dummy))           ::      dummy2,dummy3
!             integer                             ::      u
!             real(kind=real64)                                ::      rr
!             integer                             ::      ioerr,jj,n0
!             if (present(m)) then
!                 n = m
!             else
!                 n = 0
!             end if
!             u = findGoodUnit()
!             open(unit=u,status="scratch")
!                 write(unit=u,fmt="(a)") dummy
!             rewind(unit=u,iostat=ioerr)
!                 if( ioerr==0) read(unit=u,fmt=*,iostat=ioerr) rr
!             close(unit=u)
!             if (ioerr == 0) then
!                 ! have successfully found another real to read
!                 n = n + 1
!                 ra(n) = rr
!                 dummy2 = adjustl(dummy)
!                 jj = index(dummy2," ")
!                 if ((jj>0).and.(jj<len_trim(dummy2))) then
!                     ! still more to parse
!                     dummy3 = dummy2(jj:) // repeat( " ",jj-1 )
!                     n0 = n
!                     call ParseRealArray(dummy3,ra,n,n0)
!                 end if
!             end if
!             return
!         end subroutine ParseRealArray

        subroutine ParseRealArray1(dummy,ra,n)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      Parse the character string dummy into an array of n reals
    !*      New version
            character(len=*),intent(in)         ::      dummy
            real(kind=real64),dimension(:),intent(inout)     ::      ra
            integer,intent(out)                 ::      n

            type(StringTokenizer)       ::      stok
            character(len=32)           ::      tok
            logical                     ::      ok
            integer                     ::      ii
            real(kind=real64)                        ::      rr
            stok = StringTokenizer_ctor(trim(dummy)," ,;"//TAB_CHARACTER//CR_CHARACTER//NULL_CHARACTER)
            n = 0
            do ii = 1,1000000
                if (.not. hasMoreTokens(stok)) exit
                call nextToken(stok,tok)
                call ParseReal(tok,rr,ok)
                if (ok) then
                    n = n + 1
                    ra(n) = rr
                    if (n==size(ra)) return
                end if
            end do
            return
        end subroutine ParseRealArray1



        subroutine ParseReal(text,r,ok)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      read the next real(kind=real64) number from the string "text"
    !*      cut this from text on exit
            character(len=*),intent(in)         ::      text
            real(kind=real64),intent(inout)                    ::      r
            logical,intent(out),optional        ::      ok
            character(len=64)           ::      str
            character(len=32)           ::      aaaa,bbbb,cccc
            real(kind=real64)                        ::      xaaa,xbbb
            integer                     ::      iccc
            integer                     ::      ii,jj
            logical                     ::      inNum
            real(kind=real64)                        ::      xx,rr
            character(len=1),dimension(17),parameter    ::      NUMBERS =           &
                    (/ "0","1","2","3","4","5","6","7","8","9",                     &
                       ".","+","-",                                                 &
                       "e","E","d","D"                                              /)
            character(len=1)            ::      cc


        !---    extract "possible" real(kind=real64) from text
            jj = 0
            inNum = .false.
            str = ""
            rr = 0.0
            if (present(ok)) ok = .false.
            do ii = 1,len(text)
                cc = text(ii:ii)
                if (.not. inNum) then
                    if (cc==" ") cycle     !   not hit anything reasonable yet
                end if
                if (any(cc==NUMBERS)) then
                    inNum = .true.
                    jj = jj + 1
                    str(jj:jj) = cc
                    if (jj == len(str)) exit
                else
                    exit
                end if
            end do
            if (jj == 0) return     !   failed to find anything
            rr = 0.0
        !---    convert "D"s to standard "e" form
            do ii = 1,jj
                cc = str(ii:ii)
                select case(cc)
                    case ("d")
                        cc = "e"
                    case ("D")
                        cc = "e"
                    case ("E")
                        cc = "e"
                end select
                str(ii:ii) = cc
            end do

        !---    attempt to parse into format
        !               aaaa.bbbb e cccc
        !       options are  1   aaaa
        !                    2   aaaa e cccc
        !                    3   .bbbb
        !                    4   .bbbb e cccc
        !                    5   aaaa.bbbb
        !                    6   aaaa.bbbb e cccc
            aaaa = "0"
            bbbb = "0"
            cccc = "0"
            ii = index(str,".")
            jj = index(str,"e")
            if (ii == 1) then       !   3,4
                if (jj>1) then      !   4
                    bbbb = str(2:jj-1)
                    cccc = str(jj+1:len_trim(str))
                else                !   3
                    bbbb = str(2:len_trim(str))
                end if
            else if (ii>1) then     !   5,6
                aaaa = str(1:ii-1)
                if (jj>1) then      !   6
                    bbbb = str(ii+1:jj-1)
                    cccc = str(jj+1:len_trim(str))
                else                !   5
                    bbbb = str(ii+1:len_trim(str))
                end if
            else                    !   1,2
                if (jj>1) then      !   2
                    aaaa = str(1:jj-1)
                    cccc = str(jj+1:len_trim(str))
                else                !   1
                    aaaa = str(1:len_trim(str))
                end if
            end if

        !---    convert to integers
            xaaa = 0.0d0
            do ii = 1,len_trim(aaaa)
                cc = aaaa(ii:ii)
                jj = ( iachar(aaaa(ii:ii)) - iachar("0") )
                if ((jj>=0).and.(jj<=9)) xaaa = xaaa*10 + jj
            end do
            iccc = 0
            do ii = 1,len_trim(cccc)
                jj = ( iachar(cccc(ii:ii)) - iachar("0") )
                if ((jj>=0).and.(jj<=9)) iccc = iccc*10 + jj
            end do
            if (cccc(1:1)=="-") iccc = -iccc
            xx = 1.0d0
            xbbb = 0.0d0
            do ii = 1,len_trim(bbbb)
                xx = xx * 0.1d0
                xbbb = xbbb + xx*( iachar(bbbb(ii:ii)) - iachar("0") )
            end do

        !---    put it all together
            rr = ( xaaa + xbbb ) * (10.0d0**iccc)
            if (aaaa(1:1)=="-") rr = -rr
            if (abs(rr)==0) rr = 0
            if (present(ok)) ok = .true.
            r = rr
            return
        end subroutine ParseReal



!        subroutine parseComplex(dummy,rr,ok)
!    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
!    !*      Parse dummy into a complex variable.
!            character(len=*),intent(in)         ::      dummy
!            complex,intent(out)                 ::      rr
!            logical,intent(out),optional        ::      ok
!            integer                             ::      u
!            integer                             ::      ioerr
!            u = findGoodUnit()
!            rr = 0.0
!            open(unit=u,status="scratch")
!                write(unit=u,fmt="(a)") dummy
!            rewind(unit=u,iostat=ioerr)
!                if( ioerr==0) read(unit=u,fmt=*,iostat=ioerr) rr
!            close(unit=u)
!            if (present(ok)) ok = (ioerr == 0)
!            return
!        end subroutine parseComplex

!-------


! 
!         recursive subroutine parseIntArray(dummy,ia,n,m)
!     !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
!     !*      Parse the character string dummy into an array of n integers.
!     !*      m is used internally.
!             character(len=*),intent(in)         ::      dummy
!             integer,dimension(:),intent(inout)  ::      ia
!             integer,intent(out)                 ::      n
!             integer,intent(in),optional         ::      m       ! number of filled entries in ia
!             character(len=len(dummy))           ::      dummy2,dummy3
!             integer                             ::      u
!             integer                             ::      ii,ioerr,jj,nn
!             if (present(m)) then
!                 n = m
!             else
!                 n = 0
!             end if
!             u = findGoodUnit()
!             open(unit=u,status="scratch")
!                 write(unit=u,fmt="(a)") dummy
!             rewind(unit=u,iostat=ioerr)
!                 if( ioerr==0) read(unit=u,fmt=*,iostat=ioerr) ii
!             close(unit=u)
!             if (ioerr == 0) then
!                 ! have successfully found another int to read
!                 n = n + 1
!                 ia(n) = ii
!                 dummy2 = adjustl(dummy)
!                 jj = index(dummy2," ")
!                 if ((jj>0).and.(jj<len_trim(dummy2))) then
!                     ! still more to parse
!                     dummy3 = dummy2(jj:) // repeat( " ",jj-1 )
!                     nn = n ; call parseIntArray(dummy3,ia,n,nn)
!                 end if
!             end if
!             return
!         end subroutine parseIntArray


        subroutine parseIntArray1(dummy,ia,n)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      Parse the character string dummy into an array of n integers.
    !*      New version
            character(len=*),intent(in)         ::      dummy
            integer,dimension(:),intent(inout)  ::      ia
            integer,intent(out)                 ::      n

            type(StringTokenizer)       ::      stok
            character(len=32)           ::      tok
            logical                     ::      ok
            integer                     ::      ii,jj
            stok = StringTokenizer_ctor(dummy," ,;"//TAB_CHARACTER//CR_CHARACTER//NULL_CHARACTER)
            n = 0
            do ii = 1,1000000
                if (.not. hasMoreTokens(stok)) exit
                call nextToken(stok,tok)
                call parseInt(tok,jj,ok)
                if (ok) then
                    n = n + 1
                    ia(n) = jj
                    if (n==size(ia)) return
                end if
            end do
            return
        end subroutine parseIntArray1


!         subroutine parseInt0(dummy,ii)
!     !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
!     !*      Parse dummy into a integer variable. Original version
!             character(len=*),intent(in)         ::      dummy
!             integer,intent(out)                 ::      ii
!             integer                             ::      u
!             integer                             ::      ioerr
!             u = findGoodUnit()
!             ii = 0
!             open(unit=u,status="scratch")
!                 write(unit=u,fmt="(a)") dummy
!             rewind(unit=u,iostat=ioerr)
!                 if( ioerr==0) read(unit=u,fmt=*,iostat=ioerr) ii
!             close(unit=u)
!             return
!         end subroutine parseInt0

        subroutine parseInt(dummy,i,ok)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      Parse dummy into a integer variable. New version
            character(len=*),intent(in)         ::      dummy
            integer,intent(inout)               ::      i
            logical,intent(out),optional        ::      ok
            integer                 ::      ii,jj
            integer                 ::      letter
            logical                 ::      numberStarted,negative
            numberStarted = .false.
            negative = .false.
            ii = 0
            if (present(ok)) ok = .false.
            do jj = 1,len_trim(dummy)
                letter = iachar( dummy(jj:jj) )
                if ( (letter>=iachar("0")).and.(letter<=iachar("9")) ) then
                    ii = ii*10 + (letter - iachar("0"))
                    numberStarted = .true.
                else
                    if (numberStarted) then
                        exit
                    else if (letter==iachar("-")) then
                        negative = .true.
                    end if
                end if
            end do
            if (negative) ii = -ii
            if (present(ok)) ok = numberStarted
            if (numberStarted) i = ii
            return
        end subroutine parseInt



!-------


        subroutine parseLogicalArray(dummy,la,n)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      Parse the character string dummy into an array of n logicals.
            character(len=*),intent(in)         ::      dummy
            logical,dimension(:),intent(inout)  ::      la
            integer,intent(out)                 ::      n
            character(len=len(dummy))           ::      dummy2
            type(StringTokenizer)               ::      st
            logical                             ::      ok
            st = StringTokenizer_ctor(trim(dummy)," ,")
            n = 0
            do
                if (.not.hasMoreTokens(st)) exit
                call nextToken(st,dummy2)
                n = n + 1
                call parseLogical(dummy2,la(n),ok)
                if (.not. ok) then
                    n = n - 1
                    exit
                end if                
                if (n==size(la)) exit
            end do
            return
        end subroutine parseLogicalArray

        subroutine parseCharacterArray(dummy,ca,n)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      Parse the character string dummy into an array of n character strings
            character(len=*),intent(in)         ::      dummy
            character(len=*),dimension(:),intent(inout)  ::      ca
            integer,intent(out)                 ::      n
            character(len=len(dummy))           ::      dummy2
            type(StringTokenizer)               ::      st
            st = StringTokenizer_ctor(trim(dummy)," ,")
            n = 0
            do
                if (.not.hasMoreTokens(st)) exit
                call nextToken(st,dummy2)
                n = n + 1
                ca(n) = dummy2
                if (n==size(ca)) exit
            end do
            return
        end subroutine parseCharacterArray



        subroutine parseLogical(dummy,ll,ok)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      Parse dummy into a Logical variable.
            character(len=*),intent(in)         ::      dummy
            Logical,intent(out)                 ::      ll
            logical,intent(out),optional        ::      ok
            if (len_trim(dummy)==0) then
                ll = .false.
                if (present(ok)) ok = .false.
            end if
            ll = ( (trim(dummy)=="t") .or. (trim(dummy)=="true") .or. (trim(dummy)==".true.")    &
                  .or. (trim(dummy)=="T") .or. (trim(dummy)=="TRUE") .or. (trim(dummy)==".TRUE.") )
            if (present(ok)) then
                if (ll) then
                    ok = .true.
                else
                    ok = ( (trim(dummy)=="f") .or. (trim(dummy)=="false") .or. (trim(dummy)==".false.")    &
                          .or. (trim(dummy)=="F") .or. (trim(dummy)=="FALSE") .or. (trim(dummy)==".FALSE.") )
                end if
            end if

            return
        end subroutine parseLogical

    !---
    
    
            pure subroutine fatoa1(dat,n,m,a,k) 
        !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
        !*      return the floating point array dat in format f nnnn.mmmm
        !*      with n digits before decimal place, no leading zeroes
        !*      and m digits after decimal place.
        !*      n,m <= 9
        !*      on output k is length of string
        !*      do not call with n=0 or m=0
                real(kind=real64),dimension(:),intent(in)        ::      dat
                integer,intent(in)                  ::      n,m
                character(len=*),intent(out)        ::      a
                integer,intent(out)                 ::      k
                integer         ::      ii,kk
                k = 0
                do ii = 1,size(dat)
                    call ftoa0(dat(ii),n,m,a(k+1:),kk) 
                    k = k + kk
                end do
                return
            end subroutine fatoa1
            

            pure subroutine fatoa1a(dat,n,m,a,k) 
        !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
        !*      return the floating point array dat in format f nnnn.mmmm
        !*      with n digits before decimal place, no leading zeroes
        !*      and m digits after decimal place.
        !*      n,m <= 9
        !*      on output k is length of string
        !*      do not call with n=0 or m=0
                real(kind=real64),dimension(:,:),intent(in)        ::      dat
                integer,intent(in)                  ::      n,m
                character(len=*),intent(out)        ::      a
                integer,intent(out)                 ::      k
                integer         ::      ii,jj,kk
                k = 0
                do ii = 1,size(dat,dim=2)
                    do jj = 1,size(dat,dim=1)
                        call ftoa0(dat(jj,ii),n,m,a(k+1:),kk) 
                        k = k + kk
                    end do
                end do
                return
            end subroutine fatoa1a
            

            pure subroutine fatoa2(dat,n,m,a,k)  
        !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
        !*      return the floating point array dat in format f nnnn.mmmm
        !*      with n digits before decimal place, no leading zeroes
        !*      and m digits after decimal place.
        !*      n,m <= 9
        !*      on output k is length of string
        !*      do not call with n=0 or m=0
                real(kind=real64),dimension(:),intent(in)        ::      dat
                integer,dimension(:),intent(in)                  ::      n
                integer,intent(in)                  ::      m
                character(len=*),intent(out)        ::      a
                integer,intent(out)                 ::      k
                integer         ::      ii,kk
                k = 0
                do ii = 1,size(dat)
                    call ftoa0(dat(ii),n(ii),m,a(k+1:),kk) 
                    k = k + kk
                end do
                return
            end subroutine fatoa2
            

            pure subroutine fatoa3(dat,n,m,a,k) 
        !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
        !*      return the floating point array dat in format f nnnn.mmmm
        !*      with n digits before decimal place, no leading zeroes
        !*      and m digits after decimal place.
        !*      n,m <= 9
        !*      on output k is length of string
        !*      do not call with n=0 or m=0
                real(kind=real64),dimension(:),intent(in)        ::      dat
                integer,dimension(:),intent(in)                  ::      n,m
                character(len=*),intent(out)        ::      a
                integer,intent(out)                 ::      k
                integer         ::      ii,kk
                k = 0
                do ii = 1,size(dat)
                    call ftoa0(dat(ii),n(ii),m(ii),a(k+1:),kk) 
                    k = k + kk 
                end do
                return
            end subroutine fatoa3                        
            
            pure subroutine ftoa0(dat,n,m,a,k) 
        !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
        !*      return the floating point number dat in format f nnnn.mmmm
        !*      with n digits before decimal place, no leading zeroes
        !*      and m digits after decimal place.
        !*      n,m <= 9
        !*      on output k is length of string
        !*      do not call with n=0 or m=0
                real(kind=real64),intent(in)        ::      dat
                integer,intent(in)                  ::      n,m
                character(len=*),intent(out)        ::      a
                integer,intent(out)                 ::      k
                
                integer             ::      uu,vv,m9,n9
                real(kind=real64)   ::      absdat
                
                
                m9 = min(m,9)
                n9 = min(n,9)
                
                a = ""
                if (dat<0) a(2:2) = "-"
                absdat = abs(dat)                
                uu = int( absdat )      !   unsigned integer before decimal place
                
                absdat = absdat - uu
                vv = nint( absdat*10**m9 )
                if (vv == 10**m9 ) then          !   FP number was actually 123.999995 and should be rounded up to 124.00000
                    call utoa(uu+1,n9,a(3:),k)
                    a(3+k:3+k+m9) = "."//repeat("0",m9)
                else
                    call utoa(uu,n9,a(3:),k)
                    a(3+k:3+k) = "."
                    call utoalz(vv,m9,a(4+k:))
                end if
                
                k = 3+k+m9
                !print *,"ftoa ",dat,""""//a(1:k)//""""
                
                return
            end subroutine ftoa0
                
            
            pure subroutine utoa(u,n,a,k)
        !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
        !*      unsigned integer u with known fixed maximum digits n to string a
        !*      eg utoa(1234,4) = "1234", utoa(123,6) = "123"
        !*      on output k is the length of the string
        !*      do not call with n = 0, or u < 0 or you'll get garbage.
                integer,intent(in)              ::      u
                integer,intent(in)              ::      n
                character(len=*),intent(out)    ::      a
                integer,intent(out)             ::      k
                integer         ::      jj,uu
                logical         ::      ok
                
                a = "" ; k = 0
                select case(n)
                    case (1)
                        k = 1
                        a = achar( u + 48 )
                    case (2)                         
                        jj = u/10
                        if (jj>0) then
                            k = 2
                            a(1:1) = achar( jj + 48 )
                            a(2:2) = achar( u - jj*10 + 48 )
                        else
                            k = 1
                            a = achar( u + 48 )
                        end if
                    case (3)                            
                        k = 1 ; ok = .false.
                        uu = u ; jj = uu/100                        
                        a(k:k) = achar( jj + 48 )
                        ok = (jj>0) ; if (ok) k = k + 1
                        uu = uu - jj*100 ; jj = uu/10
                        a(k:k) = achar( jj + 48 )
                        ok = ok .or. (jj>0) ; if (ok) k = k + 1
                        uu = uu - jj*10                       
                        a(k:k) = achar( uu + 48 )
                    case (4)
                        k = 1 ; ok = .false. 
                        uu = u ; jj = uu/1000                       
                        a(k:k) = achar( jj + 48 )
                        ok = (jj>0) ; if (ok) k = k + 1
                        uu = uu - jj*1000 ; jj = uu/100
                        a(k:k) = achar( jj + 48 )
                        ok = ok .or. (jj>0) ; if (ok) k = k + 1
                        uu = uu - jj*100 ; jj = uu/10                      
                        a(k:k) = achar( jj + 48 )
                        ok = ok .or. (jj>0) ; if (ok) k = k + 1
                        uu = uu - jj*10   
                        a(k:k) = achar( uu + 48 )
                    case (5)
                        k = 1 ; ok = .false. 
                        uu = u ; jj = uu/10000                       
                        a(k:k) = achar( jj + 48 )
                        ok = (jj>0) ; if (ok) k = k + 1
                        uu = uu - jj*10000 ; jj = uu/1000
                        a(k:k) = achar( jj + 48 )
                        ok = ok .or. (jj>0) ; if (ok) k = k + 1
                        uu = uu - jj*1000 ; jj = uu/100
                        a(k:k) = achar( jj + 48 )
                        ok = ok .or. (jj>0) ; if (ok) k = k + 1
                        uu = uu - jj*100 ; jj = uu/10                      
                        a(k:k) = achar( jj + 48 )
                        ok = ok .or. (jj>0) ; if (ok) k = k + 1
                        uu = uu - jj*10   
                        a(k:k) = achar( uu + 48 )                   
                    case (6)
                        k = 1 ; ok = .false. 
                        uu = u ; jj = uu/100000                       
                        a(k:k) = achar( jj + 48 )
                        ok = (jj>0) ; if (ok) k = k + 1
                        uu = uu - jj*100000 ; jj = uu/10000
                        a(k:k) = achar( jj + 48 )
                        ok = ok .or. (jj>0) ; if (ok) k = k + 1
                        uu = uu - jj*10000 ; jj = uu/1000
                        a(k:k) = achar( jj + 48 )
                        ok = ok .or. (jj>0) ; if (ok) k = k + 1
                        uu = uu - jj*1000 ; jj = uu/100
                        a(k:k) = achar( jj + 48 )
                        ok = ok .or. (jj>0) ; if (ok) k = k + 1
                        uu = uu - jj*100 ; jj = uu/10                      
                        a(k:k) = achar( jj + 48 )
                        ok = ok .or. (jj>0) ; if (ok) k = k + 1
                        uu = uu - jj*10   
                        a(k:k) = achar( uu + 48 )                                                 
                    case (7)
                        k = 1 ; ok = .false. 
                        uu = u ; jj = uu/1000000                       
                        a(k:k) = achar( jj + 48 )
                        ok = (jj>0) ; if (ok) k = k + 1
                        uu = uu - jj*1000000 ; jj = uu/100000
                        a(k:k) = achar( jj + 48 )
                        ok = ok .or. (jj>0) ; if (ok) k = k + 1
                        uu = uu - jj*100000 ; jj = uu/10000
                        a(k:k) = achar( jj + 48 )
                        ok = ok .or. (jj>0) ; if (ok) k = k + 1
                        uu = uu - jj*10000 ; jj = uu/1000
                        a(k:k) = achar( jj + 48 )
                        ok = ok .or. (jj>0) ; if (ok) k = k + 1
                        uu = uu - jj*1000 ; jj = uu/100
                        a(k:k) = achar( jj + 48 )
                        ok = ok .or. (jj>0) ; if (ok) k = k + 1
                        uu = uu - jj*100 ; jj = uu/10                      
                        a(k:k) = achar( jj + 48 )
                        ok = ok .or. (jj>0) ; if (ok) k = k + 1
                        uu = uu - jj*10   
                        a(k:k) = achar( uu + 48 )                                                 
                    case (8)
                        k = 1 ; ok = .false. 
                        uu = u ; jj = uu/10000000                       
                        a(k:k) = achar( jj + 48 )
                        ok = (jj>0) ; if (ok) k = k + 1
                        uu = uu - jj*10000000 ; jj = uu/1000000
                        a(k:k) = achar( jj + 48 )
                        ok = ok .or. (jj>0) ; if (ok) k = k + 1
                        uu = uu - jj*1000000 ; jj = uu/100000
                        a(k:k) = achar( jj + 48 )
                        ok = ok .or. (jj>0) ; if (ok) k = k + 1
                        uu = uu - jj*100000 ; jj = uu/10000
                        a(k:k) = achar( jj + 48 )
                        ok = ok .or. (jj>0) ; if (ok) k = k + 1
                        uu = uu - jj*10000 ; jj = uu/1000
                        a(k:k) = achar( jj + 48 )
                        ok = ok .or. (jj>0) ; if (ok) k = k + 1
                        uu = uu - jj*1000 ; jj = uu/100
                        a(k:k) = achar( jj + 48 )
                        ok = ok .or. (jj>0) ; if (ok) k = k + 1
                        uu = uu - jj*100 ; jj = uu/10                      
                        a(k:k) = achar( jj + 48 )
                        ok = ok .or. (jj>0) ; if (ok) k = k + 1
                        uu = uu - jj*10   
                        a(k:k) = achar( uu + 48 )                                                 
                    case (9)
                        k = 1 ; ok = .false. 
                        uu = u ; jj = uu/100000000                       
                        a(k:k) = achar( jj + 48 )
                        ok = (jj>0) ; if (ok) k = k + 1
                        uu = uu - jj*100000000 ; jj = uu/10000000
                        a(k:k) = achar( jj + 48 )
                        ok = ok .or. (jj>0) ; if (ok) k = k + 1
                        uu = uu - jj*10000000 ; jj = uu/1000000
                        a(k:k) = achar( jj + 48 )
                        ok = ok .or. (jj>0) ; if (ok) k = k + 1
                        uu = uu - jj*1000000 ; jj = uu/100000
                        a(k:k) = achar( jj + 48 )
                        ok = ok .or. (jj>0) ; if (ok) k = k + 1
                        uu = uu - jj*100000 ; jj = uu/10000
                        a(k:k) = achar( jj + 48 )
                        ok = ok .or. (jj>0) ; if (ok) k = k + 1
                        uu = uu - jj*10000 ; jj = uu/1000
                        a(k:k) = achar( jj + 48 )
                        ok = ok .or. (jj>0) ; if (ok) k = k + 1
                        uu = uu - jj*1000 ; jj = uu/100
                        a(k:k) = achar( jj + 48 )
                        ok = ok .or. (jj>0) ; if (ok) k = k + 1
                        uu = uu - jj*100 ; jj = uu/10                      
                        a(k:k) = achar( jj + 48 )
                        ok = ok .or. (jj>0) ; if (ok) k = k + 1
                        uu = uu - jj*10   
                        a(k:k) = achar( uu + 48 ) 
                end select
                
                 
                return
            end subroutine utoa                    
                
            pure subroutine utoalz(u,n,a)
        !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
        !*      unsigned integer u with known fixed maximum digits n to string a
        !*      with leading zeroes
        !*      eg utoa(1234,4) = "1234", utoa(123,6) = "000123"   
        !*      do not call with n = 0, or u < 0 or you'll get garbage.
                integer,intent(in)              ::      u
                integer,intent(in)              ::      n
                character(len=*),intent(out)    ::      a                
                integer         ::      jj,uu
                
                a = "" 
                select case(n)
                    case (1)
                        a = achar( u + 48 )
                    case (2)                         
                        jj = u/10
                        a(1:1) = achar( jj + 48 )
                        a(2:2) = achar( u - jj*10 + 48 )
                    case (3)                            
                        uu = u ; jj = uu/100                        
                        a(1:1) = achar( jj + 48 )
                        uu = uu - jj*100 ; jj = uu/10
                        a(2:2) = achar( jj + 48 )
                        uu = uu - jj*10                       
                        a(3:3) = achar( uu + 48 )
                    case (4)
                        uu = u ; jj = uu/1000                        
                        a(1:1) = achar( jj + 48 )
                        uu = uu - jj*1000 ; jj = uu/100
                        a(2:2) = achar( jj + 48 )
                        uu = uu - jj*100 ; jj = uu/10                       
                        a(3:3) = achar( jj + 48 )
                        uu = uu - jj*10                       
                        a(4:4) = achar( uu + 48 )
                    case (5)
                        uu = u ; jj = uu/10000                        
                        a(1:1) = achar( jj + 48 )
                        uu = uu - jj*10000 ; jj = uu/1000
                        a(2:2) = achar( jj + 48 )
                        uu = uu - jj*1000 ; jj = uu/100                       
                        a(3:3) = achar( jj + 48 )
                        uu = uu - jj*100 ; jj = uu/10                       
                        a(4:4) = achar( jj + 48 )
                        uu = uu - jj*10                       
                        a(5:5) = achar( uu + 48 )
                    case (6)
                        uu = u ; jj = uu/100000                        
                        a(1:1) = achar( jj + 48 )
                        uu = uu - jj*100000 ; jj = uu/10000
                        a(2:2) = achar( jj + 48 )
                        uu = uu - jj*10000 ; jj = uu/1000                       
                        a(3:3) = achar( jj + 48 )
                        uu = uu - jj*1000 ; jj = uu/100                       
                        a(4:4) = achar( jj + 48 )
                        uu = uu - jj*100 ; jj = uu/10                       
                        a(5:5) = achar( jj + 48 )
                        uu = uu - jj*10                       
                        a(6:6) = achar( uu + 48 )
                    case (7)
                        uu = u ; jj = uu/1000000                        
                        a(1:1) = achar( jj + 48 )
                        uu = uu - jj*1000000 ; jj = uu/100000
                        a(2:2) = achar( jj + 48 )
                        uu = uu - jj*100000 ; jj = uu/10000                       
                        a(3:3) = achar( jj + 48 )
                        uu = uu - jj*10000 ; jj = uu/1000                      
                        a(4:4) = achar( jj + 48 )
                        uu = uu - jj*1000 ; jj = uu/100                       
                        a(5:5) = achar( jj + 48 )
                        uu = uu - jj*100 ; jj = uu/10                       
                        a(6:6) = achar( jj + 48 )
                        uu = uu - jj*10                       
                        a(7:7) = achar( uu + 48 )
                    case (8)
                        uu = u ; jj = uu/10000000                        
                        a(1:1) = achar( jj + 48 )
                        uu = uu - jj*10000000 ; jj = uu/1000000
                        a(2:2) = achar( jj + 48 )
                        uu = uu - jj*1000000 ; jj = uu/100000                       
                        a(3:3) = achar( jj + 48 )
                        uu = uu - jj*100000 ; jj = uu/10000                       
                        a(4:4) = achar( jj + 48 )
                        uu = uu - jj*10000 ; jj = uu/1000                       
                        a(5:5) = achar( jj + 48 )
                        uu = uu - jj*1000 ; jj = uu/100                       
                        a(6:6) = achar( jj + 48 )
                        uu = uu - jj*100 ; jj = uu/10                       
                        a(7:7) = achar( jj + 48 )
                        uu = uu - jj*10                       
                        a(8:8) = achar( uu + 48 )
                    case (9)
                        uu = u ; jj = uu/100000000                        
                        a(1:1) = achar( jj + 48 )
                        uu = uu - jj*100000000 ; jj = uu/10000000
                        a(2:2) = achar( jj + 48 )
                        uu = uu - jj*10000000 ; jj = uu/1000000                       
                        a(3:3) = achar( jj + 48 )
                        uu = uu - jj*1000000 ; jj = uu/100000                       
                        a(4:4) = achar( jj + 48 )
                        uu = uu - jj*100000 ; jj = uu/10000                       
                        a(5:5) = achar( jj + 48 )
                        uu = uu - jj*10000 ; jj = uu/1000                       
                        a(6:6) = achar( jj + 48 )
                        uu = uu - jj*1000 ; jj = uu/100                       
                        a(7:7) = achar( jj + 48 )
                        uu = uu - jj*100 ; jj = uu/10                       
                        a(8:8) = achar( jj + 48 )
                        uu = uu - jj*10                       
                        a(9:9) = achar( uu + 48 )
                end select
                return
            end subroutine utoalz                   
                
            pure subroutine itoa0(u,n,a,k)
        !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
        !*      signed integer u with known fixed maximum digits n to string a
        !*      eg itoa(1234,4) = "1234", itoa(-123,6) = "-123"
        !*      on output k is the length of the string
        !*      do not call with n = 0, or u < 0 or you'll get garbage.
                integer,intent(in)              ::      u
                integer,intent(in)              ::      n
                character(len=*),intent(out)    ::      a
                integer,intent(out)             ::      k
                call utoa(abs(u),min(n,9),a,k)
                if (u<0) then
                    a = "-"//trim(a)
                    k = k+1
                end if                
                return
            end subroutine itoa0                    
                

            pure subroutine iatoa1(dat,n,a,k)  
        !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
        !*      return the int array dat with known fixed maximum digits n to string a
        !*      on output k is length of string
        !*      do not call with n=0  
                integer,dimension(:),intent(in)     ::      dat
                integer,intent(in)                  ::      n
                character(len=*),intent(out)        ::      a
                integer,intent(out)                 ::      k
                integer         ::      ii,kk
                k = 0
                do ii = 1,size(dat)
                    call itoa0(dat(ii),n,a(k+1:),kk) 
                    !print *,"itoa0 ",ii,n,k,dat(ii)," """//trim(a(k+1:))//""" ",kk
                    k = k + kk + 1
                    a(k:k) = " "
                end do
                return
            end subroutine iatoa1
                        
        
        
!-------

!         function findGoodUnit() result(i)
! !-------^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
! !*          find a good (unopened) unit to read/write
!             integer             ::      i
!             logical             ::      lod
!             ! find a good unit number ( i )
!             i = 20
!             do
!                 inquire (opened = lod,unit = i)
!                 if (.not. lod) exit
!                 i = i + 1
!             end do
!             return
!         end function findGoodUnit


!-------
    end module NBAX_StringTokenizers


