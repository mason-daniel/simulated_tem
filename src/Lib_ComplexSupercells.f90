
    module Lib_ComplexSupercells
!---^^^^^^^^^^^^^^^^^^^^^^^^^^^^
!*      provide bare functionality for a simple periodic 3d supercell with repeat vectors A
!*      constructed from a set of primitive unit cells b
!*      where A = b N   and N is an integer matrix.

        use iso_fortran_env
        implicit none
        private
        
    !---
    
        public      ::      ComplexSupercell_ctor
        public      ::      delete
        public      ::      report
     
    !---
    
        public      ::      cellToRealSpace
        public      ::      realSpaceToCell
        public      ::      getNnodes 
        public      ::      getN
        public      ::      getA,getb,getConventionalCell,estimateUnitCell
        public      ::      getLatticeVector,getSuperVector
        public      ::      getSuperCellSideLength
        public      ::      getCellSideLength  
        public      ::      wrapPBC
        public      ::      getNode
        public      ::      pointInSupercell
        public      ::      unitCellVolume
        public      ::      supercellVolume
        public      ::      suggestSupercellOrientedWithN
        public      ::      completeBasis
        public      ::      inverse3Mat
        public      ::      determinant3Mat
        
         
    !---    it is convenient for me to index lattice points to nodes using a compressed 30 bit format. This kind indicates a bitmap is used.
        integer,private,parameter       ::      int60 = int64      
        
    !---
    
        logical,public          ::      ComplexSupercell_dbg = .false.
                   
    !---
    
        type,public     ::      ComplexSupercell
            private
            integer,dimension(3,3)                      ::      N           !   number of repeats
            integer                                     ::      nNodes      !   number of unit cells in supercell
            real(kind=real64),dimension(3,3)            ::      A           !   supercell vectors - vector 1 is A(1:3,1) etc
            real(kind=real64),dimension(3,3)            ::      iA          !   inverse supercell vectors             
            real(kind=real64),dimension(3,3)            ::      b           !   unit cell vectors - vector 1 is b(1:3,1) etc
            real(kind=real64),dimension(3,3)            ::      ib          !   inverse unit cell vectors             
            integer,dimension(3,3)                      ::      iNxdetN     !   inverse repeats x determinant of N.
            integer(kind=int60),dimension(:),pointer    ::      indx        !   which unit cell does each node point to ( 30 bit bitmap )
        end type ComplexSupercell
        
    !---
    
        interface ComplexSupercell_ctor
            module procedure    ComplexSupercell_null
            module procedure    ComplexSupercell_ctor0
            module procedure    ComplexSupercell_ctor1
        end interface
                
        interface delete
            module procedure    delete0
        end interface
        
        interface report
            module procedure    report0
        end interface
        

        interface   getN
            module procedure        getN0
        end interface
        
        interface   getA
            module procedure        getA0
        end interface
        
        interface   getConventionalCell
            module procedure        getConventionalCell0
        end interface
        
        interface   estimateUnitCell
            module procedure        estimateUnitCell0
        end interface
 
        interface   getb
            module procedure        getb0
        end interface
        
        interface   getNnodes
            module procedure        getNnodes0
        end interface   
                 
        interface   getNode
            module procedure        getNode0
        end interface   

        interface   realSpaceToCell
            module procedure        realSpaceToCell0
        end interface

        interface   cellToRealSpace
            module procedure        cellToRealSpace0
            module procedure        cellToRealSpace1
        end interface

        interface   unitCellVolume
            module procedure        unitCellVolume0
        end interface
        interface   superCellVolume
            module procedure        superCellVolume0
        end interface
         
        
        interface   wrapPBC
            module procedure        wrapPBC1
        end interface
         
        interface   pointInSupercell
            module procedure        pointInSupercell0
        end interface
!        

        interface  determinant3Mat 
            module procedure        d_determinant3Mat
            module procedure        i_determinant3Mat
        end interface
        
        interface  inverse3Mat 
            module procedure        d_inverse3Mat
            module procedure        i_inverse3Matxdet
        end interface
                           
        interface   suggestSupercellOrientedWithN
            module procedure        suggestSupercellOrientedWithN0
        end interface  
        
        interface   getCellSideLength
            module procedure    getCellSideLength0
        end interface
        
        interface   getSupercellSideLength
            module procedure    getSupercellSideLength0
        end interface
        
        interface   getLatticeVector
            module procedure    getLatticeVector0
        end interface
        
        interface   getSuperVector
            module procedure    getSuperVector0
        end interface
        
    contains
!---^^^^^^^^

        function ComplexSupercell_null() result(this)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      provide empty null instance with no allocated storage
            type(ComplexSupercell)              ::      this
            real(kind=real64),dimension(3,3)    ::      identity3 = reshape( (/1.0d0,0.0d0,0.0d0 , 0.0d0,1.0d0,0.0d0 , 0.0d0,0.0d0,1.0d0/),(/3,3/) )
            this%A = identity3
            this%iA = identity3
            this%b = identity3
            this%ib = identity3
            this%N = nint(identity3)
            this%iNxdetN = nint(identity3)
            this%nNodes = 1
            nullify(this%indx)
            return
        end function ComplexSupercell_null
        
        function ComplexSupercell_ctor1(A,N) result(this)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      construct a complex supercell given the size of the simulation box
    !*      and the primitive cell repeats A = b N
            type(ComplexSupercell)                          ::      this
            real(kind=real64),dimension(3,3),intent(in)     ::      A
            integer,dimension(3,3),intent(in)               ::      N
            integer,dimension(3,3)              ::      iNxdetN
            real(kind=real64),dimension(3,3)    ::      bb
            integer                             ::      detN
            call inverse3Mat( N,iNxdetN )
            detN = determinant3Mat( N ) 
            if (detN == 0) stop "Lib_ComplexSupercells::ComplexSupercell_ctor1 ERROR - attempt to define supercell with det(N) = 0"
            bb = matmul( A, real( iNxdetN,kind=real64) )
            bb = bb / detN
            this =  ComplexSupercell_ctor0(A,bb) 
            return
        end function ComplexSupercell_ctor1
            
                         
        function ComplexSupercell_ctor0(A,b , noindx) result(this)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      construct a commensurate supercell with periodic repeat vectors A
    !*      and primitive unit cell b.
    !*      optionally don't complete the indexing ( Warning!! This will produce an incomplete construction!!! )
            type(ComplexSupercell)                          ::      this
            real(kind=real64),dimension(3,3),intent(in)     ::      A
            real(kind=real64),dimension(3,3),intent(in)     ::      b
            logical,intent(in),optional                     ::      noindx
            integer         ::      i1,i2,i3
            integer         ::      l1,l2,l3,u1,u2,u3 , nn , m1,m2,m3  !, bb
            
            real(kind=real64),dimension(3)      ::      xx
            
            this%A = A
            call inverse3Mat( this%A,this%iA ) 
            call inverse3Mat( b,this%ib )                           !   not the finalised ib
            this%N = nint( matmul( this%ib,this%A ) )               !   commensurate number of primitive cell repeats
            
            
            this%nNodes = determinant3Mat( this%N )                        
            call inverse3Mat( this%N,this%iNxdetN )
                             
            this%b = matmul( this%A,this%iNxdetN )/this%nNodes      !   commensurate unit cell close to b
            call inverse3Mat( this%b,this%ib )
             
            
            if (ComplexSupercell_dbg) then
                print *,"Lib_ComplexSupercells::ComplexSupercell_ctor0() info - check node periodicity"
                xx(1:3) = this%b(1:3,1)*this%N(1,1) + this%b(1:3,2)*this%N(2,1) + this%b(1:3,3)*this%N(3,1)
                print *,"   direction 1: ",xx,wrapPBC(this,xx,realSpace=.true.)            
                xx(1:3) = this%b(1:3,1)*this%N(1,2) + this%b(1:3,2)*this%N(2,2) + this%b(1:3,3)*this%N(3,2)
                print *,"   direction 2: ",xx,wrapPBC(this,xx,realSpace=.true.)
                xx(1:3) = this%b(1:3,1)*this%N(1,3) + this%b(1:3,2)*this%N(2,3) + this%b(1:3,3)*this%N(3,3)
                print *,"   direction 3: ",xx,wrapPBC(this,xx,realSpace=.true.)
            end if
            
            
            if (present(noindx)) then
                nullify(this%indx)
                if (noindx) return
            end if
            
            allocate( this%indx( this%nNodes ) )
            
            l1 = sum( min(-1,this%N(1,:)) )
            l2 = sum( min(-1,this%N(2,:)) )
            l3 = sum( min(-1,this%N(3,:)) )
            
            
            u1 = sum( max(0,this%N(1,:)+1) )
            u2 = sum( max(0,this%N(2,:)+1) )
            u3 = sum( max(0,this%N(3,:)+1) )
            nn = 0           
           
            do i3 = l3,u3
                do i2 = l2,u2
                    do i1 = l1,u1
                   
                        m1 = this%iNxdetN(1,1)*i1 + this%iNxdetN(1,2)*i2 + this%iNxdetN(1,3)*i3
                        if ( (m1<0).or.(m1>=this%nNodes) ) cycle             !   outside sim cell
                        m2 = this%iNxdetN(2,1)*i1 + this%iNxdetN(2,2)*i2 + this%iNxdetN(2,3)*i3
                        if ( (m2<0).or.(m2>=this%nNodes) ) cycle
                        m3 = this%iNxdetN(3,1)*i1 + this%iNxdetN(3,2)*i2 + this%iNxdetN(3,3)*i3 
                        if ( (m3<0).or.(m3>=this%nNodes) ) cycle
                        
                        nn = nn + 1
                        !bb = signedTripletToThirtyBit( i1,i2,i3 )   
                        this%indx( nn ) = signedTripletToSixtyBit( i1,i2,i3 )   
                        
                    end do
                end do
            end do
            if (ComplexSupercell_dbg) then
                print *,"Lib_ComplexSupercells::ComplexSupercell_ctor0() info - nodes expected ",this%nNodes," nodes found ",nn
                do i1 = 1,nn-1
                    if (any(this%indx(i1+1:)==this%indx(i1))) print *,"repeated index ",i1,this%indx(i1)
                end do
            end if
            
            !print *,"Lib_ComplexSupercells::ComplexSupercell_ctor0() info - nodes expected ",this%nNodes," nodes found ",nn
            
            
            return
        end function ComplexSupercell_ctor0
                       
    !---
    
        subroutine delete0(this)
    !---^^^^^^^^^^^^^^^^^^^^^^^^
    !*      deallocate dynamic memory
            type(ComplexSupercell),intent(inout)    ::      this
            if (associated(this%indx)) deallocate(this%indx)
            this = ComplexSupercell_null()
            return
        end subroutine delete0
        
    !--- 
    
        subroutine report0(this,u,o,verbose)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      make a simple dump of the supercell (default to screen) 
    !*      optionally verbose
    !*      optionally with a left margin of o spaces.
            type(ComplexSupercell),intent(in)    ::      this
            integer,intent(in),optional         ::      u,o
            logical,intent(in),optional         ::      verbose

            integer         ::      uu,oo
            
            real(kind=real64),parameter     ::      DEGINRAD = 180.0d0/3.141592654d0
            uu = 6 ; if (present(u)) uu = u
            oo = 0 ; if (present(o)) oo = o


            write(unit=uu,fmt='(a,i8,a)')   repeat(" ",oo)//"ComplexSupercell [nodes = ",this%nNodes,"]"
            write(unit=uu,fmt='(a)')        repeat(" ",oo+2)//"periodic repeat vectors A "
            write(unit=uu,fmt='(a,3f14.8)') repeat(" ",oo+2),this%A(1,1),this%A(1,2),this%A(1,3)
            write(unit=uu,fmt='(a,3f14.8)') repeat(" ",oo+2),this%A(2,1),this%A(2,2),this%A(2,3)
            write(unit=uu,fmt='(a,3f14.8)') repeat(" ",oo+2),this%A(3,1),this%A(3,2),this%A(3,3)
            write(unit=uu,fmt='(a)')        repeat(" ",oo+2)//"primitive unit cell vectors b"
            write(unit=uu,fmt='(a,3f14.8)') repeat(" ",oo+2),this%b(1,1),this%b(1,2),this%b(1,3)
            write(unit=uu,fmt='(a,3f14.8)') repeat(" ",oo+2),this%b(2,1),this%b(2,2),this%b(2,3)
            write(unit=uu,fmt='(a,3f14.8)') repeat(" ",oo+2),this%b(3,1),this%b(3,2),this%b(3,3)
            write(unit=uu,fmt='(a)')        repeat(" ",oo+2)//"unit cell vector repeats A = b n"
            write(unit=uu,fmt='(a,3i14)')   repeat(" ",oo+2),this%n(1,1),this%n(1,2),this%n(1,3)
            write(unit=uu,fmt='(a,3i14)')   repeat(" ",oo+2),this%n(2,1),this%n(2,2),this%n(2,3)
            write(unit=uu,fmt='(a,3i14)')   repeat(" ",oo+2),this%n(3,1),this%n(3,2),this%n(3,3)
            if (present(verbose)) then
                if (verbose) then
                    write(unit=uu,fmt='(a)')        repeat(" ",oo+2)//"repeat lengths"
                    write(unit=uu,fmt='(a,3f14.8)') repeat(" ",oo+2),norm2(this%A(:,1)),norm2(this%A(:,2)),norm2(this%A(:,3))
                    write(unit=uu,fmt='(a)')        repeat(" ",oo+2)//"repeat angles (deg)"
                    write(unit=uu,fmt='(a,3f14.8)') repeat(" ",oo+2),DEGINRAD*acos( dot_product(this%A(:,2),this%A(:,3))/(norm2(this%A(:,2))*norm2(this%A(:,3))) )          &
                                                                    ,DEGINRAD*acos( dot_product(this%A(:,3),this%A(:,1))/(norm2(this%A(:,3))*norm2(this%A(:,1))) )          &
                                                                    ,DEGINRAD*acos( dot_product(this%A(:,1),this%A(:,2))/(norm2(this%A(:,1))*norm2(this%A(:,2))) )
                    write(unit=uu,fmt='(a)')        repeat(" ",oo+2)//"unit cell angles (deg)"
                    write(unit=uu,fmt='(a,3f14.8)') repeat(" ",oo+2),DEGINRAD*acos( dot_product(this%b(:,2),this%b(:,3))/(norm2(this%b(:,2))*norm2(this%b(:,3))) )          &
                                                                    ,DEGINRAD*acos( dot_product(this%b(:,3),this%b(:,1))/(norm2(this%b(:,3))*norm2(this%b(:,1))) )          &
                                                                    ,DEGINRAD*acos( dot_product(this%b(:,1),this%b(:,2))/(norm2(this%b(:,1))*norm2(this%b(:,2))) )
                        
                end if
            end if
            return
        end subroutine report0

    
    !---
    

        pure function unitCellVolume0( this ) result(V)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      return volume of the unit cell
            type(ComplexSupercell),intent(in)                ::      this
            real(kind=real64)                               ::      V
            V = determinant3Mat(this%b)
            return
        end function unitCellVolume0

        pure function superCellVolume0( this ) result(V)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      return volume of the super cell
            type(ComplexSupercell),intent(in)                ::      this
            real(kind=real64)                               ::      V
            V = determinant3Mat(this%A)
            return
        end function superCellVolume0

    !---
                
    
        pure function cellToRealSpace0( this,x,wrap ) result( y )
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      given a point as a multiple of unit cells, return the position in real space
            type(ComplexSupercell),intent(in)               ::      this
            real(kind=real64),dimension(3),intent(in)       ::      x
            real(kind=real64),dimension(3)                  ::      y
            logical,intent(in),optional                     ::      wrap
            
            y(1:3) = this%b(1:3,1)*x(1) + this%b(1:3,2)*x(2) + this%b(1:3,3)*x(3)
            if (present(wrap)) then
                if (wrap) y = wrapPBC(this,y,realSpace=.true.) 
            end if
            return
        end function cellToRealSpace0
    
        pure function cellToRealSpace1( this,c,wrap ) result( y )
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      given a point as a multiple of unit cells, return the position in real space
            type(ComplexSupercell),intent(in)               ::      this
            integer,dimension(3),intent(in)                 ::      c
            real(kind=real64),dimension(3)                  ::      y
            logical,intent(in),optional                     ::      wrap
            
            y(1:3) = this%b(1:3,1)*c(1) + this%b(1:3,2)*c(2) + this%b(1:3,3)*c(3)
            if (present(wrap)) then
                if (wrap) y = wrapPBC(this,y,realSpace=.true.) 
            end if
            return
        end function cellToRealSpace1
        
        pure function realSpaceToCell0( this,y,wrap ) result( x )
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ 
    !*      given a point in real space, return the position in unit cell space.    
            type(ComplexSupercell),intent(in)               ::      this
            real(kind=real64),dimension(:),intent(in)       ::      y
            real(kind=real64),dimension(3)                  ::      x
            logical,intent(in),optional                     ::      wrap
            x(1:3) = this%ib(1:3,1)*y(1) + this%ib(1:3,2)*y(2) + this%ib(1:3,3)*y(3)
            if (present(wrap)) then
                if (wrap) x = wrapPBC(this,x,realSpace=.false.) 
            end if
            return
        end function realSpaceToCell0
         
!         
            
    !---
        
        pure function getb0(this) result(b)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      return unit cell
            type(ComplexSupercell),intent(in)          ::      this
            real(kind=real64),dimension(3,3)    ::      b
            b = this%b
            return
        end function getb0
        
        pure function getN0(this) result(n)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      return primitive cell repeats
            type(ComplexSupercell),intent(in)          ::      this
            integer,dimension(3,3)                      ::      n
            n(1:3,1) = this%N(1:3,1) 
            n(1:3,2) = this%N(1:3,2) 
            n(1:3,3) = this%N(1:3,3) 
            return
        end function getN0

        pure function getA0(this) result(a)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      return super cell
            type(ComplexSupercell),intent(in)          ::      this
            real(kind=real64),dimension(3,3)    ::      a
            a(1:3,1) = this%A(1:3,1) 
            a(1:3,2) = this%A(1:3,2) 
            a(1:3,3) = this%A(1:3,3) 
            return
        end function getA0
        
        pure function getConventionalCell0(this,primitive) result(c)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      return a conventional unit cell given the primitive unit cell
            type(ComplexSupercell),intent(in)           ::      this
            real(kind=real64),dimension(3,3),intent(in) ::      primitive
            real(kind=real64),dimension(3,3)            ::      c
            call inverse3Mat(primitive,c)
            c = matmul( this%b,c )
             
            return
        end function getConventionalCell0
        
    
        pure function estimateUnitCell0(this,Nx,primitive) result(b)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      return a primitive unit cell given a number of cell repeats in each dimension
    !*      b = c p where c(:,1) = A(:,1)/Nx(1) etc
    !*      note: does not set this%b
            type(ComplexSupercell),intent(in)           ::      this
            integer,dimension(3),intent(in)             ::      Nx
            real(kind=real64),dimension(3,3),intent(in) ::      primitive
            real(kind=real64),dimension(3,3)            ::      b
            
            b(:,1) = this%A(:,1)/Nx(1)      !   this is actually a conventional cell ...
            b(:,2) = this%A(:,2)/Nx(2)
            b(:,3) = this%A(:,3)/Nx(3)
            b = matmul( b,primitive )       !   ... which I multiply to get a unit cell.
             
            return
        end function estimateUnitCell0
        
    
    
        
    !---
    
    
        
        pure function getSuperCellSideLength0(this,i) result(a)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      return super cell side i length
            type(ComplexSupercell),intent(in)        ::      this
            integer,intent(in)                      ::      i
            real(kind=real64)                       ::      a
            if (i==1) then
                a = norm2(this%A(1:3,1)) 
            else if (i==2) then
                a = norm2(this%A(1:3,2)) 
            else
                a = norm2(this%A(1:3,3))                                
            end if
            return
        end function getSuperCellSideLength0
        
        pure function getCellSideLength0(this,i) result(b)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      return unit cell side i length
            type(ComplexSupercell),intent(in)        ::      this
            integer,intent(in)                      ::      i
            real(kind=real64)                       ::      b
            b = norm2(this%b(1:3,i))
            return
        end function getCellSideLength0
        
        pure function getLatticeVector0(this,i) result(b)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      return lattice vector i
            type(ComplexSupercell),intent(in)        ::      this
            integer,intent(in)                      ::      i
            real(kind=real64),dimension(3)          ::      b
            if (i==1) then
                b = this%b(1:3,1)
            else if (i==2) then
                b = this%b(1:3,2)
            else
                b = this%b(1:3,3)                              
            end if
            return
        end function getLatticeVector0
        
        pure function getSuperVector0(this,i) result(a)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      return supercell vector i
            type(ComplexSupercell),intent(in)        ::      this
            integer,intent(in)                      ::      i
            real(kind=real64),dimension(3)          ::      a
            if (i==1) then
                a = this%A(1:3,1)
            else if (i==2) then
                a = this%A(1:3,2)
            else
                a = this%A(1:3,3)                               
            end if
            return
        end function getSuperVector0
        
        
    !---

        pure function getNnodes0(this) result(n)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      returns the number of lattice points in the supercell
            type(ComplexSupercell),intent(in)      ::      this
            integer                         ::      n
            n = this%nNodes
            return
        end function getNnodes0

!         
!         pure function getNode0(this,n) result(x)
!     !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
!     !*      given the number of the node n, return the realspace position of the node 
!             type(ComplexSupercell),intent(in)       ::      this
!             integer,intent(in)                      ::      n
!             real(kind=real64),dimension(3)          ::      x
!             integer         ::      ii,bb
!             
!             bb = this%indx(n)
!             
!             ii = tenBitToSignedInt( iand( bb,1023 ) )      
!             x(1:3) = this%b(1:3,1)*ii     
!             ii = tenBitToSignedInt( iand( ishft(bb,-10),1023 ) )
!             x(1:3) = x(1:3) + this%b(1:3,2)*ii   
!             ii = tenBitToSignedInt( iand( ishft(bb,-20),1023 ) ) 
!             x(1:3) = x(1:3) + this%b(1:3,3)*ii   
!             return
!         end function getNode0


        
        pure function getNode0(this,n) result(x)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      given the number of the node n, return the realspace position of the node 
            type(ComplexSupercell),intent(in)       ::      this
            integer,intent(in)                      ::      n
            real(kind=real64),dimension(3)          ::      x
            integer                 ::      ii 
            integer(kind=int60)     ::      bb
            bb = this%indx(n)
            
            ii = twentyBitToSignedInt( iand( bb,1048575_int60 ) )      
            x(1:3) = this%b(1:3,1)*ii     
            ii = twentyBitToSignedInt( iand( ishft(bb,-21_int60),1048575_int60 ) )
            x(1:3) = x(1:3) + this%b(1:3,2)*ii   
            ii = twentyBitToSignedInt( iand( ishft(bb,-42_int60),1048575_int60 ) ) 
            x(1:3) = x(1:3) + this%b(1:3,3)*ii   
            return
        end function getNode0
        
      
    !---
    
        pure function pointInSupercell0( this,x ) result(is)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      returns true if wrapping x back into the supercell, does not actually move it
    !*      cf wrapPBS which moves the point but doesn't tell you.
            type(ComplexSupercell),intent(in)               ::      this
            real(kind=real64),dimension(3),intent(in)       ::      x
            logical                                         ::      is
    
            integer             ::      ix,iy,iz
                        
            ix = floor( this%iA(1,1)*x(1) + this%iA(1,2)*x(2) + this%iA(1,3)*x(3) )
            iy = floor( this%iA(2,1)*x(1) + this%iA(2,2)*x(2) + this%iA(2,3)*x(3) )
            iz = floor( this%iA(3,1)*x(1) + this%iA(3,2)*x(2) + this%iA(3,3)*x(3) )
            is = (ix*ix + iy*iy + iz*iz == 0)
            return
        end function pointInSupercell0
            
        pure function wrapPBC1( this,x,realspace ) result( xp )
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      wraps a position x back into the periodic supercell
            type(ComplexSupercell),intent(in)               ::      this
            real(kind=real64),dimension(3),intent(in)       ::      x
            logical,intent(in)                              ::      realspace
            real(kind=real64),dimension(3)                  ::      xp
            integer             ::      ix,iy,iz
            real(kind=real64)   ::      idetN
            
            if (realspace) then
                ix = floor( this%iA(1,1)*x(1) + this%iA(1,2)*x(2) + this%iA(1,3)*x(3) )
                iy = floor( this%iA(2,1)*x(1) + this%iA(2,2)*x(2) + this%iA(2,3)*x(3) )
                iz = floor( this%iA(3,1)*x(1) + this%iA(3,2)*x(2) + this%iA(3,3)*x(3) )
                xp(1:3) = x(1:3) - this%A(1:3,1)*ix - this%A(1:3,2)*iy - this%A(1:3,3)*iz
            else
                idetN = 1.0d0/this%nNodes
                ix = floor( (this%iNxdetN(1,1)*x(1) + this%iNxdetN(1,2)*x(2) + this%iNxdetN(1,3)*x(3))*idetN )
                iy = floor( (this%iNxdetN(2,1)*x(1) + this%iNxdetN(2,2)*x(2) + this%iNxdetN(2,3)*x(3))*idetN )
                iz = floor( (this%iNxdetN(3,1)*x(1) + this%iNxdetN(3,2)*x(2) + this%iNxdetN(3,3)*x(3))*idetN )               
                xp(1:3) = x(1:3) - this%N(1:3,1)*ix - this%N(1:3,2)*iy - this%N(1:3,3)*iz                
            end if
            
            return
        end function wrapPBC1
        
        

    !---        
        
         subroutine suggestSupercellOrientedWithN0( this, n ,a0, that , mx,my,mz)
    !----^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      suggest a new unit cell with the same periodic supercell volume
    !*      but orthorhombic and oriented with z-axis along n
    !*      Aim for unit cell side lengths of a0, but optionally allow (mx) subdivisions per side   
            type(ComplexSupercell),intent(in)               ::      this 
            real(kind=real64),dimension(3),intent(in)       ::      n
            real(kind=real64),intent(in)                    ::      a0
            integer,intent(in),optional                     ::      mx,my,mz
            type(ComplexSupercell),intent(out)              ::      that
        
            real(kind=real64),dimension(3)      ::      dx1,dx2,dx3
            real(kind=real64),dimension(3,3)    ::      bb


        !---    suggest three axes
            call suggestPlane( this%A ,n, dx1,dx2,dx3 )
                        
        !---    get the cubic cell
            bb(1:3,1) = a0*dx1(1:3)/norm2(dx1)
            bb(1:3,2) = a0*dx2(1:3)/norm2(dx2)
            bb(1:3,3) = a0*dx3(1:3)/norm2(dx3)
            
        !---    optionally add subdivisions per axis
            if (present(mx)) then
                that = ComplexSupercell_ctor( this%A,bb,noindx=.true. )   !   partial construction- fixes a but does not allocate indx.
                bb = that%b
                bb(:,1) = bb(:,1)/mx
                bb(:,2) = bb(:,2)/my
                bb(:,3) = bb(:,3)/mz
            end if
            
        !---    construct commensurate unit cell
            that = ComplexSupercell_ctor( this%A,bb )
            
! 
            return
        end subroutine suggestSupercellOrientedWithN0
                
            

        subroutine suggestPlane( a_super,n, dx1,dx2,dx3 )
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      given that we want to sample points on a plane y.n = 0
    !*      suggest three vectors dx1,dx2,dx3
    !*      which cover the supercell
            real(kind=real64),dimension(3,3),intent(in)     ::      a_super
            real(kind=real64),dimension(3),intent(in)       ::      n
            real(kind=real64),dimension(3),intent(out)      ::      dx1,dx2,dx3

            real(kind=real64),dimension(3)      ::      nn,pp
            real(kind=real64)                   ::      d1,d2,d3,vol
            
        !---    which face of the supercell is closest to n?
        !   use this to decide dx1 and dx2.
        
        !   find projection of given normal on each face normal
            nn = cross_product( a_super(1:3,1),a_super(1:3,2) )
            nn = nn/norm2(nn)
            pp(1) = dot_product( nn(1:3),n(1:3) )
            nn = cross_product( a_super(1:3,2),a_super(1:3,3) )
            nn = nn/norm2(nn)
            pp(2) = dot_product( nn(1:3),n(1:3) )
            nn = cross_product( a_super(1:3,3),a_super(1:3,1) )
            nn = nn/norm2(nn)
            pp(3) = dot_product( nn(1:3),n(1:3) )
            
            
        !   which has greatest projection?
            pp = abs(pp)            
            nn(1:3) = n(1:3) / norm2(n(1:3))
            
           ! print *,"Lib_ComplexSupercells::suggestPlane() info - n ",nn," p ",pp
            
            
            if (pp(1)>max(pp(2),pp(3))) then
                !   n has greatest projection on a1,a2 face
                dx1 = a_super(1:3,1)                !   hint to push vector 1 along a1
                call completeBasis( nn,dx1,dx2 )                                
                !   find projection along a1,a2 vectors
                d1 = dot_product( dx1,a_super(1:3,1) )
                d2 = dot_product( dx2,a_super(1:3,2) )                

            else if (pp(2)>max(pp(3),pp(1))) then
                !   n has greatest projection on a2,a3 face
                dx1 = a_super(1:3,2)                !   hint to push vector 1 along a2
                call completeBasis( nn,dx1,dx2 )                
                !   find projection along a1,a2 vectors
                d1 = dot_product( dx1,a_super(1:3,2) )
                d2 = dot_product( dx2,a_super(1:3,3) )                
            
            else !if (pp(3)>max(pp(1),pp(2))) then
                !   n has greatest projection on a3,a1 face
                dx1 = a_super(1:3,3)                !   hint to push vector 1 along a3
                call completeBasis( nn,dx1,dx2 )                
                !   find projection along a1,a2 vectors
                d1 = dot_product( dx1,a_super(1:3,3) )
                d2 = dot_product( dx2,a_super(1:3,1) )                
            
            
            end if

        !---    can now complete the job by ensuring resultant box has same volume as initial box            
            vol = determinant3Mat(a_super)
            d3 = vol / (d1*d2)
            
            dx1 = dx1*d1
            dx2 = dx2*d2
            dx3 = nn*d3
            
            

            return
        end subroutine suggestPlane

    !---

        pure subroutine d_inverse3Mat(M,N)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
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
        end subroutine d_inverse3Mat
        
        pure subroutine i_inverse3Matxdet(M,N)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      monadic operator
    !*      finds the inverse of an integer three matrix multiplied by its determinant 
    !*      ( aka the matrix of cofactors ) 
            integer,dimension(3,3),intent(in)       ::  M
            integer,dimension(3,3),intent(out)      ::  N

            N(1,1)   =  M(2,2)*M(3,3) - M(2,3)*M(3,2)   
            N(2,1)   =  M(2,3)*M(3,1) - M(2,1)*M(3,3)   
            N(3,1)   =  M(2,1)*M(3,2) - M(2,2)*M(3,1)   
                                                        
            N(1,2)   =  M(1,3)*M(3,2) - M(1,2)*M(3,3)   
            N(2,2)   =  M(1,1)*M(3,3) - M(1,3)*M(3,1)   
            N(3,2)   =  M(1,2)*M(3,1) - M(1,1)*M(3,2)   
                                                        
            N(1,3)   =  M(1,2)*M(2,3) - M(1,3)*M(2,2)   
            N(2,3)   =  M(1,3)*M(2,1) - M(1,1)*M(2,3)   
            N(3,3)   =  M(1,1)*M(2,2) - M(1,2)*M(2,1)   

            return
        end subroutine i_inverse3Matxdet

        pure function d_determinant3Mat(M) result(d)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      returns the determinant of real matrix M
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
        end function d_determinant3Mat
        
        pure function i_determinant3Mat(M) result(d)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      returns the determinant of integer matrix M
            integer,dimension(3,3),intent(in)      ::      M
            integer                                ::      d
            d     = M(1,1)*( M(2,2)*M(3,3) - M(2,3)*M(3,2) )        &
                  + M(1,2)*( M(2,3)*M(3,1) - M(2,1)*M(3,3) )        &
                  + M(1,3)*( M(2,1)*M(3,2) - M(2,2)*M(3,1) )
            return
        end function i_determinant3Mat

    !---
    
            
        subroutine completeBasis( z,x,y )
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      given the z-axis vector, complete the basis to provide a triplet x,y,z
    !*      if on input x is set, then use this as a hint to attempt to place x along this direction
            real(kind=real64),dimension(3),intent(in)       ::      z
            real(kind=real64),dimension(3),intent(inout)    ::      x
            real(kind=real64),dimension(3),intent(out)      ::      y

            real(kind=real64),dimension(3)      ::      zz        !   normalised
            real(kind=real64)                   ::      xxxx,zzzz  !   vector lengths
            real(kind=real64)                   ::      zdotx

        !---    check we have a normalised z as input
            zzzz = norm2(z)
            if (zzzz == 0) then
                !   can't do anything with zero input vector
                x = (/ 1,0,0 /)
                y = (/ 0,1,0 /)
                return
            end if
            zz = z/zzzz

        !---    check for a sensible hint for the x-direction
            xxxx = norm2(x)
            zdotx = dot_product( zz,x )
            if ( (xxxx < 0.001d0 ).or.(abs(zdotx) >= xxxx*0.999d0) ) then
                if (ComplexSupercell_dbg) print *,"Lib_ComplexSupercells::completeBasis info - haven't got a good hint for x-vector"
                !   haven't got a good hint. Make random hint
                if (abs(zz(3))>max(abs(zz(1)),abs(zz(2)))) then
                    !   z points along 3-axis
                    x(1) = 1.0d0
                    x(3) = 0.0d0
                    if (abs(zz(1))>0) then
                        x(2) = - zz(2)/zz(1)
                    else
                        x(2) = 0.0d0
                    end if
                    zdotx = dot_product( zz,x )
                else if (abs(zz(2))>max(abs(zz(1)),abs(zz(3)))) then
                    !   z points along 2-axis
                    x(3) = 1.0d0
                    x(2) = 0.0d0
                    if (abs(zz(3))>0) then
                        x(1) = - zz(1)/zz(3)
                    else
                        x(1) = 0.0d0
                    end if
                    zdotx = dot_product( zz,x )
                else
                    !   z points along 1-axis
                    x(2) = 1.0d0
                    x(1) = 0.0d0
                    if (abs(zz(2))>0) then
                        x(3) = - zz(3)/zz(2)
                    else
                        x(3) = 0.0d0
                    end if
                    zdotx = dot_product( zz,x )
                end if
                !print *,"suggest ",x,norm2(x),zdotx
            end if


        !---    remove projection of z on x-direction
            x = x - zz*zdotx
            xxxx = norm2(x)
            x = x/xxxx

        !---    construct y-direction
            y = cross_product( zz,x )

            if (ComplexSupercell_dbg) then
                print *,"Lib_ComplexSupercells::completeBasis info -"
                print *,"    x :",x
                print *,"    y :",y
                print *,"    z :",z
            end if
            return
        end subroutine completeBasis
        

        pure function cross_product(x,y) result(z)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      ronseal
            real(kind=real64),dimension(3),intent(in)       ::      x,y
            real(kind=real64),dimension(3)                  ::      z
            z(1) = x(2)*y(3) - x(3)*y(2)
            z(2) = x(3)*y(1) - x(1)*y(3)
            z(3) = x(1)*y(2) - x(2)*y(1)
            return
        end function cross_product

!         
!         pure function signedTripletToThirtyBit( n1,n2,n3 ) result( b )
!     !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
!     !*      convert the signed triplet of integers ( n1,n2,n3 ) = (-512:511) 
!     !*      into a single 30 bit number
!             integer,intent(in)                     ::      n1,n2,n3
!             integer(kind=int60)                    ::      b
!             
!             
!             b =        signedIntToTenBit( n1 )            &
!               + ishft( signedIntToTenBit( n2 ),10 )       &
!               + ishft( signedIntToTenBit( n3 ),20 )
!               
!             return
!         end function signedTripletToThirtyBit
!             
!          
!         
!         elemental function tenBitToSignedInt( b ) result ( i )
!     !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^    
!     !*      convert a ten bit unsigned integer (0:1023) to a signed integer (-512:511)
!             integer,intent(in)                  ::      b
!             integer                             ::      i
!             i = iand( b,511 ) - iand( b,512 )
!             return
!         end function tenBitToSignedInt
!             
!         elemental function signedIntToTenBit( i ) result ( b )
!     !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^    
!     !*      convert a signed integer (-512:511) to a ten bit unsigned integer (0:1023) 
!             integer,intent(in)                  ::      i
!             integer                ::      b
!             if (i<0) then
!                 b = iand( i + 1024 , 1023 )
!             else
!                 b = iand( i , 1023 )                
!             end if
!             return
!         end function signedIntToTenBit
!         

        pure function signedTripletToSixtyBit( n1,n2,n3 ) result( b )
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      convert the signed triplet of integers ( n1,n2,n3 ) = (-524288:524287)   
    !*      into a single 30 bit number
            integer,intent(in)                     ::      n1,n2,n3
            integer(kind=int60)                    ::      b
            
            
            b =        signedIntToTwentyBit( n1 )            &
              + ishft( signedIntToTwentyBit( n2 ),21_int60 )       &
              + ishft( signedIntToTwentyBit( n3 ),42_int60 )
              
            return
        end function signedTripletToSixtyBit
            
        
!         
!         pure function sixtyBitToSignedTriplet( b ) result( n )
!     !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
!     !*      convert a single 30 bit number
!     !*      int a signed triplet of integers ( n1,n2,n3 ) = (-524288:524287)     
!             integer(kind=int60),intent(in)                      ::      b
!             integer,dimension(3)                                ::      n
!             
!             n(1) = twentyBitToSignedInt( iand( b,524287_int60 ) )
!             n(2) = twentyBitToSignedInt( iand( ishft(b,-20_int60),524287_int60 ) )
!             n(3) = twentyBitToSignedInt( iand( ishft(b,-40_int60),524287_int60 ) )
!             
!             return
!         end function sixtyBitToSignedTriplet
!             
        
        
        elemental function twentyBitToSignedInt( b ) result ( i )
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^    
    !*      convert a 20 bit unsigned integer (0:1048575) to a signed integer (-524288:524287)
            integer(kind=int60),intent(in)                  ::      b
            integer                                         ::      i
            i = int( iand( b,524287_int60 ) - iand( b,524288_int60 ),kind=int32 )
            return
        end function twentyBitToSignedInt
            
        elemental function signedIntToTwentyBit( i ) result ( b )
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^    
    !*      convert a signed integer (-524288:524287) to a 20 bit unsigned integer (0:1048575) 
            integer,intent(in)                  ::      i
            integer(kind=int60)                 ::      b
            if (i<0) then
                b = iand( int(i,kind=int60) + 1048576_int60 , 1048575_int60 )
            else
                b = iand( int(i,kind=int60) , 1048575_int60 )                
            end if
            return
        end function signedIntToTwentyBit
                    
            
                
!         pure function asBinary1(h) result(txt)
!     !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
!     !*      returns a 32-bit integer as a binary string
!             integer(kind=int32),intent(in)      ::      h
!             character(len=32)                   ::      txt
!             integer     ::      ii
!             txt = ""
!             do ii = 1,32
!                 if (btest(h,ii-1)) then
!                     txt(ii:ii) = "1"
!                 else
!                     txt(ii:ii) = "0"
!                 end if
!             end do                
!             return
!         end function asBinary1  
             
    end module Lib_ComplexSupercells
    
