
    module Lib_Gvectors
!---^^^^^^^^^^^^^^^^^^^
!*      short module defining a set of g-vectors for imaging
!*      The g-vectors are defined by integer reflections [hkl] 
!*      and reciprocal lattice vectors B
!*          g = B [hkl]
!*      where B is defined from _conventional_ lattice vectors A 
!*          B = 2 pi A^-1
!*      Note that if the conventional lattice vectors have a rotation or strain, then so will the g-vectors
!*
!*      Daniel Mason
!*      (c) UKAEA September 2024
!*

        use iso_fortran_env
        implicit none
        private

        integer,public,parameter                        ::      LIB_GVECTORS_UNSET = 0
        real(kind=real64),parameter                     ::      PI = 3.14159265390d0
    
        public      ::      Gvectors_ctor  
        public      ::      report          
        public      ::      delete          
        public      ::      clone           !   deep copy this = that

        public      ::      sanityCheck     !   tests whether there is a -g for every g

        public      ::      getn  
        public      ::      getg
        public      ::      gethkl
        public      ::      gethjkl
        public      ::      whichg
        public      ::      getMinusg
        public      ::      isPositiveg             !   returns true if this g-vector has no negative, or is [000], or the negative is later in the set
        public      ::      nPositiveg              !   counts number of "positive" g-vectors in the set
        public      ::      rotate
        public      ::      isMillerBravais


        type,public     ::      Gvectors
            private
            integer                                     ::      n = 0           !   number of g-vectors in set
            real(kind=real64),dimension(3,3)            ::      B               !   (3,3)   reciprocal lattice vectors
            logical                                     ::      MillerBravais
            integer,dimension(:,:),pointer              ::      hkl             !   (3,n)   Miller indices of g vectors. Always set, as we need reflection g = B [hkl]
            integer,dimension(:,:),pointer              ::      hjkl            !   (4,n)   Miller-Bravais indices of g vectors. Only set if needed
            integer,dimension(:),pointer                ::      minusg          !   (n)     which g-vector (if any) in the list is my -g?  g(:,minusg(ii)) = - g(:,ii).  minusg(ii) = 0 if not in set.
        end type
        
    
        interface   Gvectors_ctor
            module procedure        Gvectors_null
            module procedure        Gvectors_ctor0
            module procedure        Gvectors_ctor1
        end interface
    
        interface   report
            module procedure        report0
        end interface
        
        interface   delete
            module procedure        delete0
        end interface

        interface   clone
            module procedure        clone0
        end interface

        interface   getn
            module procedure        getn0
        end interface

        interface   rotate
            module procedure        rotate0
        end interface

        interface   sanityCheck
            module procedure        sanityCheck0
        end interface 

        interface   getg
            module procedure        getg0
            module procedure        getg1
            module procedure        getg2
            module procedure        getg3
        end interface
        
        interface   gethkl
            module procedure        gethkl0
            module procedure        gethkl1
            module procedure        gethkl2
        end interface
        
        interface   gethjkl
            module procedure        gethjkl0
            module procedure        gethjkl1
        end interface

        interface   whichg
            module procedure        whichg0
            module procedure        whichg1
        end interface

    contains
!---^^^^^^^^

        
        function Gvectors_null() result(this)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^    
            type(Gvectors)                     ::      this
            this%n = 0
            this%B = reshape( (/1,0,0,0,1,0,0,0,1/),(/3,3/) )
            this%MillerBravais = .false.
            nullify(this%hkl)
            nullify(this%hjkl)
        !    nullify(this%g)
            nullify(this%minusg)
            return
        end function Gvectors_null
        

        function Gvectors_ctor0(a_cell_conventional,hkl) result(this)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ 
            type(Gvectors)                                  ::      this
            real(kind=real64),dimension(3,3),intent(in)     ::      a_cell_conventional     !   conventional unit cell
            integer,dimension(:,:),intent(in)               ::      hkl                     !   (1:3,1:n) or (1:4,1:n) Miller indices or Miller-Bravais indices.            
        !    real(kind=real64),dimension(3,3),intent(in),optional        ::      R

            integer         ::      ii
           

            this%n = size(hkl,dim=2)            
            
            call inverse3Mat(transpose(a_cell_conventional),this%B) 
            this%B = 2*PI*this%B  
            allocate(this%hkl(3,this%n))
            if (size(hkl,dim=1) == 4) then
                this%MillerBravais = .true.
                allocate(this%hjkl(4,this%n))
                this%hjkl(1:4,1:this%n) = hkl(1:4,1:this%n)
                do ii = 1,this%n
                    this%hkl(1:3,ii) = MillerBravaisToMiller_plane( this%hjkl(1:4,ii) )
                end do
            else
                this%MillerBravais = .false.
                nullify(this%hjkl)
                this%hkl(1:3,1:this%n) = hkl(1:3,1:this%n)
            end if
 

            allocate(this%minusg(this%n))             
            do ii = 1,this%n
                this%minusg(ii) = whichg( this,-getG(this,ii) )
            end do

            return
        end function Gvectors_ctor0



        function Gvectors_ctor1(that,doubleSet) result(this)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      clone, or produce a new set of g-vectors containing all g" = g - g'
            type(Gvectors),intent(in)                       ::      that
            logical,intent(in),optional                     ::      doubleSet
            type(Gvectors)                                  ::      this

            integer         ::      ii,jj,kk,nn
            integer,dimension(:,:),pointer                  ::      hkl_tmp
            integer,dimension(:),pointer                    ::      minusg_tmp
            real(kind=real64),dimension(3)                  ::      gg

        !---    start by cloning. Then add more vectors if necessary
            this%n = 0
            call clone( this,that )

            if (present(doubleSet)) then
                if (doubleSet) then
                    !print *,"doubling set ",that%n
                    do ii = 0,that%n
                        do jj = 0,that%n
                            if (ii == 0) then
                                if (jj == 0) then
                                    gg = 0
                                else
                                    gg = - getg(that,jj)
                                end if
                            else
                                if (jj == 0) then
                                    gg = getg(that,ii) 
                                else
                                    gg = getg(that,ii) - getg(that,jj)
                                end if
                            end if
                            
                            kk = whichg(this,gg) 
                            !print *,"Gvectors_ctor1 info - searching for ",ii,"-",jj," = ",gg,kk
                            if (kk == LIB_GVECTORS_UNSET) then   !   can't find g" = g - g'

                            !---   if no space left, then need to increase alloc                    
                                nn = size(this%minusg)
                                if (this%n == nn) then
                                    nn = nn * 2

                                    allocate(minusg_tmp(nn))
                                    minusg_tmp = LIB_GVECTORS_UNSET
                                    minusg_tmp(1:this%n) = this%minusg(1:this%n)
                                    deallocate(this%minusg)
                                    this%minusg => minusg_tmp

                                    if (this%MillerBravais) then
                                        allocate(hkl_tmp(4,nn))                                        
                                        hkl_tmp(1:4,1:this%n) = this%hjkl(1:4,1:this%n)
                                        deallocate(this%hjkl)
                                        this%hjkl => hkl_tmp
                                    end if

                                    allocate(hkl_tmp(3,nn))
                                    hkl_tmp(1:3,1:this%n) = this%hkl(1:3,1:this%n)
                                    deallocate(this%hkl)
                                    this%hkl => hkl_tmp    
                                    
                                end if

                            !---    add new g-vector to set. Note add minusg at the end
                                this%n = this%n + 1
                                this%hkl(1:3,this%n) = getHkl(this,gg)
                                if (this%MillerBravais) then
                                    this%hjkl(1:4,this%n) = MillerToMillerBravais_plane( this%hkl(1:3,this%n) )
                                end if

                            end if
                        end do
                    end do

                    do ii = 1,this%n
                        this%minusg(ii) = whichg( this,-getg(this,ii) )
                    end do
 
                    return
                end if
            end if

            call clone( this,that )

            return
        end function Gvectors_ctor1



        subroutine delete0(this)
    !---^^^^^^^^^^^^^^^^^^^^^^^^
            type(Gvectors),intent(inout)       ::      this
            if (this%n==0) return
            if (this%MillerBravais) deallocate(this%hjkl)
            deallocate(this%hkl)
            deallocate(this%minusg)
            this = Gvectors_null()
            return
        end subroutine delete0

        
        subroutine report0(this,u,o)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^
            type(Gvectors),intent(in)          ::      this
            integer,intent(in),optional         ::      u,o
            integer     ::      uu,oo
            integer     ::      ii
            uu = 6 ; if (present(u)) uu = u
            oo = 0 ; if (present(o)) oo = o
            write(unit=uu,fmt='(a,i6,a)') repeat(" ",oo)//"Gvectors [n=",this%n,"]"
            write(unit=uu,fmt='(a,a36,a6,a36)')      repeat(" ",oo+4),"reciprocal lattice vector"  
            do ii = 1,3
                write(unit=uu,fmt='(a,3f12.6,a6,3f12.6)')      repeat(" ",oo+4),this%B(ii,:) 
            end do
            if (this%MillerBravais) then
                write(unit=uu,fmt='(a,a16,a36,a12)') repeat(" ",oo+4),"Miller-Bravais","g-vector","|g|"
                do ii = 1,this%n
                    write(unit=uu,fmt='(a,4i4,5f12.6)') repeat(" ",oo+4),this%hjkl(:,ii),getg(this,ii),norm2(getg(this,ii))
                end do
            else
                write(unit=uu,fmt='(a,a12,a36,a12)') repeat(" ",oo+4),"Miller index","g-vector","|g|"
                do ii = 1,this%n
                    write(unit=uu,fmt='(a,3i4,5f12.6)') repeat(" ",oo+4),this%hkl(:,ii),getg(this,ii),norm2(getg(this,ii))
                end do
            end if
            
            return
        end subroutine report0

         
        subroutine clone0(this,that)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      make a deep copy this = that
            type(Gvectors),intent(inout)        ::      this
            type(Gvectors),intent(in)           ::      that

            if (this%n /= that%n) call delete(this)
            this%n = that%n
            this%B = that%B
            allocate(this%hkl(3,this%n))
            this%hkl(:,:) = that%hkl(:,:)
            this%MillerBravais = that%MillerBravais
            if (this%MillerBravais) then
                allocate(this%hjkl(4,this%n))
                this%hjkl(:,:) = that%hjkl(:,:)
            end if
            allocate(this%minusg(this%n))
            this%minusg(:) = that%minusg(:)
            return
        end subroutine clone0



!-------        

        pure subroutine rotate0(this,R)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      The conventional unit cell has been rotated by R.
    !*      apply this rotation to the lattice vectors, so getg() will show the correct new direction
    !*      CAUTION - you _can_ undo this with rotate(this,transpose(R)), but this will slowly accumulate errors 
            type(Gvectors),intent(inout)                    ::      this
            real(kind=real64),dimension(3,3),intent(in)     ::      R
            this%B = matmul( R, this%B )
            return
        end subroutine rotate0


        pure logical function isMillerBravais(this)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      returns true if reflections are stored as 4 component Miller-Briavais form 
    !*      (ie is hcp)
            type(Gvectors),intent(in)                       ::      this
            isMillerBravais = this%MillerBravais
            return
        end function isMillerBravais
 

        pure logical function sanityCheck0(this)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      returns true if there is a -g for every g
            type(Gvectors),intent(in)              ::      this            
            integer                             ::      ii
            real(kind=real64),dimension(3)      ::      gg

            sanityCheck0 = (this%n>0)
            do ii = 1,this%n
                gg = getg(this,ii)
                sanityCheck0 = sanityCheck0 .and. ( whichg(this,-gg)>0 )
            end do
            return
        end function sanityCheck0


        
        pure integer function getn0(this)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      return number of g-vectors
            type(Gvectors),intent(in)              ::      this            
            getn0 = this%n
            return
        end function getn0
        
        pure function getg0(this,i) result(g)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      return g-vector i
            type(Gvectors),intent(in)               ::      this    
            integer,intent(in)                      ::      i 
            real(kind=real64),dimension(3)          ::      g
            g = 0
            if ( (i>0).and.(i<=this%n) ) g = this%B(:,1)*this%hkl(1,i) + this%B(:,2)*this%hkl(2,i) + this%B(:,3)*this%hkl(3,i)
            return
        end function getg0

        pure function getg1(this) result(g)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      return all g-vectors
            type(Gvectors),intent(in)               ::      this    
            real(kind=real64),dimension(3,this%n)   ::      g
            integer             ::      ii
            do ii = 1,this%n
                g(:,ii) = getg0(this,ii)
            end do
            return
        end function getg1


        pure function getg2(this,hkl_in) result(g)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      return g-vector from the reflection
            type(Gvectors),intent(in)               ::      this    
            integer,dimension(:),intent(in)         ::      hkl_in
            real(kind=real64),dimension(3)          ::      g

            integer,dimension(3)            ::      hkl

            if (this%MillerBravais) then
                hkl = MillerBravaisToMiller_plane( hkl_in ) 
            else
                hkl = hkl_in
            end if

            g = this%B(:,1)*hkl(1) + this%B(:,2)*hkl(2) + this%B(:,3)*hkl(3)
            return
        end function getg2



        pure function getg3(a_cell_conventional,hkl_in) result(g)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      return g-vector from the reflection without constructing a Gvectors object
            real(kind=real64),dimension(3,3),intent(in)     ::      a_cell_conventional 
            integer,dimension(:),intent(in)                 ::      hkl_in
            real(kind=real64),dimension(3)                  ::      g
            real(kind=real64),dimension(3,3)    ::      BB 
            integer,dimension(3)                ::      hkl

            if (size(hkl_in)==4) then
                hkl = MillerBravaisToMiller_plane( hkl_in ) 
            else
                hkl = hkl_in
            end if

            call inverse3Mat(transpose(a_cell_conventional),BB) 
            BB = 2*PI*BB  

            g = BB(:,1)*hkl(1) + BB(:,2)*hkl(2) + BB(:,3)*hkl(3)

            return
        end function getg3


 
        pure function gethkl0(this,i) result(hkl)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      return g-vector i
            type(Gvectors),intent(in)               ::      this    
            integer,intent(in)                      ::      i 
            integer,dimension(3)                    ::      hkl
            hkl = 0
            if ( (i>0).and.(i<=this%n) ) hkl = this%hkl(:,i)
            return
        end function gethkl0

        pure function gethkl1(this) result(hkl)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      return g-vectors
            type(Gvectors),intent(in)               ::      this    
            integer,dimension(3,this%n)             ::      hkl

            hkl(:,:) = this%hkl(:,:)
            return
        end function gethkl1

        pure function gethkl2(this,g) result(hkl)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      return reflection given g-vector. Inverse function to getg2
    !*      note: always returns [hkl], not miller bravais [hjkl]
            type(Gvectors),intent(in)                   ::      this    
            real(kind=real64),dimension(3),intent(in)   ::      g
            integer,dimension(3)                        ::      hkl
            real(kind=real64),dimension(3)      ::      xx
            real(kind=real64),dimension(3,3)    ::      iB
        !---    find inverse reciprocal lattice vectors
            call inverse3Mat(this%B,iB)
        !---    find real reflection
            xx(:) = iB(:,1)*g(1) + iB(:,2)*g(2) + iB(:,3)*g(3)
        !---    find integer reflection
            hkl(:) = nint( xx )
            return
        end function gethkl2


        pure function gethjkl0(this,i) result(hjkl)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      return g-vector i
            type(Gvectors),intent(in)               ::      this    
            integer,intent(in)                      ::      i 
            integer,dimension(4)                    ::      hjkl
            hjkl = 0
            if ( (this%MillerBravais).and.(i>0).and.(i<=this%n) ) hjkl = this%hjkl(:,i) 
            return
        end function gethjkl0

        pure function gethjkl1(this) result(hjkl)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      return g-vectors
            type(Gvectors),intent(in)               ::      this    
            integer,dimension(4,this%n)             ::      hjkl
            integer     ::  ii
            hjkl = 0
            if (this%MillerBravais) then
                do ii = 1,this%n
                    hjkl(:,ii) = gethjkl0(this,ii)
                end do
            end if
            return
        end function gethjkl1        

        pure integer function whichg0(this,g) 
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      inverse function of getg
    !*      returns the g vector index closest to vector g
            type(Gvectors),intent(in)                       ::      this    
            real(kind=real64),dimension(3),intent(in)       ::      g
            integer             ::      ii
            real(kind=real64)   ::      dd,dmin
            whichg0 = LIB_GVECTORS_UNSET
            dmin = 1.0d-8           !   don't accept any old rubbish, has to be quite close.
            do ii = 1,this%n
                dd = norm2( g(:) - getg0(this,ii) )
                if (dd<dmin) then
                    dmin = dd
                    whichg0 = ii
                end if
            end do
            return
        end function whichg0

        pure integer function whichg1(this,hkl) 
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      inverse function of getg
    !*      returns the g vector index closest to vector hkl
            type(Gvectors),intent(in)                       ::      this    
            integer,dimension(:),intent(in)                 ::      hkl
            integer             ::      ii
            whichg1 = LIB_GVECTORS_UNSET
            if (this%MillerBravais) then
                do ii = 1,this%n
                    if (all( hkl(1:4) == this%hkl(1:4,ii) )) then
                        whichg1 = ii
                        return
                    end if
                end do
            else
                do ii = 1,this%n
                    if (all( hkl(1:3) == this%hkl(1:3,ii) )) then
                        whichg1 = ii
                        return
                    end if
                end do
            end if
            return
        end function whichg1

        pure integer function getMinusg(this,i) 
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      returns the g vector index corresponding to -g, or zero if not in set
            type(Gvectors),intent(in)                       ::      this    
            integer,intent(in)                              ::      i
            getMinusg = this%minusg(i)
            return
        end function getMinusg

        pure logical function isPositiveg(this,i)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      returns true if g-vector is the "+ve" version
    !*      this is the case if no reverse exists, or if the index of -g is higher than the index of g.
            type(Gvectors),intent(in)                       ::      this    
            integer,intent(in)                              ::      i

            isPositiveg = (this%minusg(i)==LIB_GVECTORS_UNSET) .or. (this%minusg(i)>i)
            return
        end function isPositiveg

        
        pure integer function nPositiveg(this)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      returns the count of g-vector which are the "+ve" versions
    !*      this is the case if no reverse exists, if g=000, or if the index of -g is higher than the index of g.
            type(Gvectors),intent(in)                       ::      this    
            
            integer         ::      ii
            nPositiveg = 0
            do ii = 1,this%n
                if (isPositiveg(this,ii)) nPositiveg = nPositiveg + 1
            end do
         
            return
        end function nPositiveg

!-------
        
        
        pure function MillerBravaisToMiller_plane( hkjl ) result( hkl )
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      HCP specific routine
    !
    !       if a plane is defined by  d = [hkl].[uvw]
    !       then  
    !           d = [hkl]. ( 1  0 -1  0 ) [uvtw]
    !                      ( 0  1 -1  0 )
    !                      ( 0  0  0  1 )
    !        so
    !           [hkil] = (  1  0  0 ) [hkl]
    !                    (  0  1  0 )
    !                    ( -1 -1  0 )
    !                    (  0  0  1 )
    !       left-inverse to give
    !                                    
    !           [hkl]  = (  1  0  0  0 ) [hkil]
    !                    (  0  1  0  0 )
    !                    (  0  0  0  1 )
    
            integer,dimension(4),intent(in)   ::      hkjl
            integer,dimension(3)              ::      hkl
 
            hkl(1) = hkjl(1) 
            hkl(2) = hkjl(2)  
            hkl(3) = hkjl(4)     
                        
            return
        end function MillerBravaisToMiller_plane               
        
        
        pure function MillerToMillerBravais_plane( hkl ) result( hkil )
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      HCP specific routine
    
            integer,dimension(3),intent(in)   ::      hkl
            integer,dimension(4)              ::      hkil
 
            hkil(1) = hkl(1)
            hkil(2) = hkl(2)
            hkil(3) = (-hkl(1) - hkl(2))
            hkil(4) = hkl(3)     
            
            
            return
        end function MillerToMillerBravais_plane                     
                    
                    
        pure subroutine inverse3Mat(M,N)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      utility monadic operator
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
    !*      utility returns the determinant of M
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

    end module Lib_Gvectors
