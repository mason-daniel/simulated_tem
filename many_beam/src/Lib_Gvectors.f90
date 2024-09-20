
    module Lib_Gvectors
!---^^^^^^^^^^^^^^^^^^^
!*      short module defining a set of g-vectors for imaging
!*
!*      Daniel Mason
!*      (c) UKAEA September 2024
!*

        use iso_fortran_env
        implicit none
        private

        real(kind=real64),parameter                     ::      PI = 3.14159265390d0
    
        public      ::      Gvectors_ctor  
        public      ::      report          
        public      ::      delete          
        public      ::      clone           !   deep copy this = that

        public      ::      sanityCheck     !   tests whether there is a -g for every g

        public      ::      setR
        public      ::      getn  
        public      ::      getg
        public      ::      gethkl
        public      ::      whichg
        


        type,public     ::      Gvectors
            private
            integer                                     ::      n               !   number of g-vectors in set
            real(kind=real64),dimension(3,3)            ::      B               !   (3,3)   reciprocal lattice vectors
            logical                                     ::      MillerBravais
            integer,dimension(:,:),pointer              ::      hkl             !   (3,n)   Miller indices of g vectors
            integer,dimension(:,:),pointer              ::      hjkl            !   (4,n)   Miller-Bravais indices of g vectors
            real(kind=real64),dimension(3,3)            ::      R               !   rotation of foil 
            real(kind=real64),dimension(:,:),pointer    ::      g               !   (3,n)   g-vectors expressed in lab frame, g = R B hkl
        end type
        
    
        interface   Gvectors_ctor
            module procedure        Gvectors_null
            module procedure        Gvectors_ctor0
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

        interface   sanityCheck
            module procedure        sanityCheck0
        end interface
        
        interface   setR
            module procedure        setR0
        end interface

        interface   getg
            module procedure        getg0
            module procedure        getg1
        end interface
        
        interface   gethkl
            module procedure        gethkl0
            module procedure        gethkl1
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
            this%R = reshape( (/1,0,0,0,1,0,0,0,1/),(/3,3/) )
            this%B = reshape( (/1,0,0,0,1,0,0,0,1/),(/3,3/) )
            this%MillerBravais = .false.
            nullify(this%hkl)
            nullify(this%hjkl)
            nullify(this%g)
            return
        end function Gvectors_null
        

        function Gvectors_ctor0(a_cell_conventional,hkl,R) result(this)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ 
            type(Gvectors)                                  ::      this
            real(kind=real64),dimension(3,3),intent(in)     ::      a_cell_conventional     !   conventional unit cell
            integer,dimension(:,:),intent(in)               ::      hkl                     !   (1:3,1:n) or (1:4,1:n) Miller indices or Miller-Bravais indices.            
            real(kind=real64),dimension(3,3),intent(in),optional        ::      R

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

            allocate(this%g(3,this%n))            
            
            if (present(R)) then
                this%R = R
            else
                this%R = reshape( (/1,0,0,0,1,0,0,0,1/),(/3,3/) )
            end if
            call setR(this,this%R)

            return
        end function Gvectors_ctor0



        subroutine delete0(this)
    !---^^^^^^^^^^^^^^^^^^^^^^^^
            type(Gvectors),intent(inout)       ::      this
            if (this%n==0) return
            if (this%MillerBravais) deallocate(this%hjkl)
            deallocate(this%hkl)
            deallocate(this%g)
            this = Gvectors_null()
            return
        end subroutine delete0

        
        subroutine report0(this,u,o)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^
            type(Gvectors),intent(in)          ::      this
            integer,intent(in),optional         ::      u,o
            integer     ::      uu,oo
            integer     ::      ii
            real(kind=real64),dimension(3,3)        ::      RB
            uu = 6 ; if (present(u)) uu = u
            oo = 0 ; if (present(o)) oo = o
            write(unit=uu,fmt='(a,i6,a)') repeat(" ",oo)//"Gvectors [n=",this%n,"]"
            write(unit=uu,fmt='(a,a36,a6,a36)')      repeat(" ",oo+4),"reciprocal lattice vector (crystal)"," ","reciprocal lattice vector (lab)"
            RB = matmul(this%R,this%B)
            do ii = 1,3
                write(unit=uu,fmt='(a,3f12.6,a6,3f12.6)')      repeat(" ",oo+4),this%B(ii,:)," ",RB(ii,:)
            end do
            if (this%MillerBravais) then

                write(unit=uu,fmt='(a,a16,a36,a16)') repeat(" ",oo+4),"Miller-Bravais","g-vector (lab frame)","|g|"
                do ii = 1,this%n
                    write(unit=uu,fmt='(a,4i4,5f12.6)') repeat(" ",oo+4),this%hjkl(:,ii),this%g(:,ii),norm2(this%g(:,ii))
                end do
            else
                write(unit=uu,fmt='(a,a12,a36,a16)') repeat(" ",oo+4),"Miller index","g-vector (lab frame)","|g|"
                do ii = 1,this%n
                    write(unit=uu,fmt='(a,3i4,5f12.6)') repeat(" ",oo+4),this%hkl(:,ii),this%g(:,ii),norm2(this%g(:,ii))
                end do
            end if
            
            return
        end subroutine report0

         
        subroutine clone0(this,that)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^
            type(Gvectors),intent(inout)        ::      this
            type(Gvectors),intent(in)           ::      that

            if (this%n /= that%n) call delete(this)
            this%n = that%n
            this%B = that%B
            this%R = that%R
            allocate(this%hkl(3,this%n))
            this%hkl(:,:) = that%hkl(:,:)
            this%MillerBravais = that%MillerBravais
            if (this%MillerBravais) then
                allocate(this%hjkl(4,this%n))
                this%hjkl(:,:) = that%hjkl(:,:)
            end if
            allocate(this%g(3,this%n))
            this%g(:,:) = that%g(:,:)
            return
        end subroutine clone0

!-------        


        pure subroutine setR0(this,R)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      sets the rotation matrix, and the g-vectors after rotation
            type(Gvectors),intent(inout)                       ::      this
            real(kind=real64),dimension(3,3),intent(in)         ::      R
            integer             ::      ii
            this%R = R

            do ii = 1,this%n
                this%g(:,ii) = this%B(:,1)*this%hkl(1,ii) + this%B(:,2)*this%hkl(2,ii) + this%B(:,3)*this%hkl(3,ii)
                this%g(:,ii) = this%R(:,1)*this%g(1,ii) + this%R(:,2)*this%g(2,ii) + this%R(:,3)*this%g(3,ii)
            end do
            

            return
        end subroutine setR0

        pure logical function sanityCheck0(this)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      returns true if there is a -g for every g
            type(Gvectors),intent(in)              ::      this            
            integer                             ::      ii
            real(kind=real64),dimension(3)      ::      gg

            sanityCheck0 = (this%n>0)
            do ii = 1,this%n
                gg = this%g(:,ii)
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

            g = this%g(:,i)
            return
        end function getg0

        pure function getg1(this) result(g)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      return all g-vectors
            type(Gvectors),intent(in)               ::      this    
            real(kind=real64),dimension(3,this%n)   ::      g

            g = this%g
            return
        end function getg1


        pure function gethkl0(this,i) result(hkl)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      return g-vector i
            type(Gvectors),intent(in)               ::      this    
            integer,intent(in)                      ::      i 
            integer,dimension(3)          ::      hkl

            hkl = this%hkl(:,i)
            return
        end function gethkl0

        pure function gethkl1(this) result(hkl)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      return g-vectors
            type(Gvectors),intent(in)               ::      this    
            integer,dimension(3,this%n)          ::      hkl

            hkl(:,:) = this%hkl(:,:)
            return
        end function gethkl1

        pure integer function whichg0(this,g) 
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      inverse function of getg
    !*      returns the g vector index closest to vector g
            type(Gvectors),intent(in)                       ::      this    
            real(kind=real64),dimension(3),intent(in)       ::      g
            integer             ::      ii
            real(kind=real64)   ::      dd,dmin
            whichg0 = 0
            dmin = huge(1.0)
            do ii = 1,this%n
                dd = norm2( g(:) - this%g(:,ii) )
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
            integer,dimension(3),intent(in)                 ::      hkl
            integer             ::      ii
            whichg1 = 0
            do ii = 1,this%n
                if (all( hkl(1:3) == this%hkl(1:3,ii) )) then
                    whichg1 = ii
                    return
                end if
            end do
            return
        end function whichg1

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