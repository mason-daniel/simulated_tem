
    module Lib_TriangulatedSurface
!---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
!*      a module defining a surface constructed from interlocking triangles
!*      Note: a TriangulatedSurface has publicly accessible fields
!*      as it is likely to be subclassed.

        use iso_fortran_env
        implicit none
        private

    !---

        public      ::      TriangulatedSurface_ctor
        public      ::      GeodesicGrid
        public      ::      Cube,Octahedron,Dodecahedron,Icosahedron,TruncatedIcosahedron
        public      ::      report
        public      ::      delete
        public      ::      getn_triangle,getn_node
        public      ::      getArea
        public      ::      distanceToSurface
        public      ::      within
        public      ::      scaleRadius
        public      ::      getTriangle
        public      ::      getTriangleInfo
        public      ::      getTriangleNodes
        public      ::      getNode,setNode
        public      ::      reallocMemory
        public      ::      quarterTriangle,thirdTriangle
!         public      ::      findBoundingSpheres
        public      ::      generateNeighbourList
        public      ::      getN_neigh
        public      ::      pointOnSurface
        public      ::      getNeigh
        public      ::      getNorm
        public      ::      findVolume,findArea
        public      ::      findVolumeAndArea
        public      ::      findAreaCorrection
        public      ::      findNormals
        public      ::      addTriangle

    !---

        integer,parameter,public            ::      ICOSAHEDRONSHAPE = 0
        integer,parameter,public            ::      DODECAHEDRONSHAPE = 1
        integer,parameter,public            ::      TRUNCATEDICOSAHEDRONSHAPE = 2
        integer,parameter,public            ::      OCTAHEDRONSHAPE = 3
        integer,parameter,public            ::      CUBESHAPE = 4
        integer,parameter,public            ::      WIGNERCELLSHAPE = 5
        integer,parameter,public            ::      USER = 6

    !---


        type,public     ::      TriangulatedSurface
            integer                             ::      original    !   original shape used for mesh
            integer                             ::      n_node      !   number of nodes
            integer                             ::      n_triangle  !   number of nodes
            integer                             ::      N_neigh     !   max number neighbours for a node
            integer                             ::      N_array     !   length of array
            integer,dimension(:,:),pointer      ::      triangle    !   (3,n_triangle) array defines triangles from triplet of nodes.
            real(kind=real64),dimension(:,:),pointer         ::      node        !   (3,n_node) array defines nodes
            integer,dimension(:,:),pointer      ::      neigh       !   (0:N_neigh,n_node) array defines triangles to which a node maps
            real(kind=real64)                                ::      rmin2,rmax2
        end type

    !---

        interface report
            module procedure        report1
            module procedure        report2
        end interface

        interface delete
            module procedure        delete1
        end interface


        interface TriangulatedSurface_ctor
            module procedure        TriangulatedSurface_ctor0
            module procedure        TriangulatedSurface_ctor1
            module procedure        TriangulatedSurface_ctor2
            module procedure        TriangulatedSurface_null
        end interface

        interface   getArea
            module procedure        getArea1
        end interface

        interface   within
            module procedure        within0
            module procedure        within1
        end interface

        interface   distanceToSurface
            module procedure        distanceToSurface0
        end interface

        interface   reallocMemory
            module procedure        reallocMemory0
        end interface

        interface   getN_neigh
            module procedure        getN_neigh0
        end interface

        interface   getNeigh
            module procedure        getNeigh0
        end interface

        interface   thirdTriangle
            module procedure        thirdTriangle0
            module procedure        thirdTriangle1
        end interface


        interface   quarterTriangle
            module procedure        quarterTriangle0
            module procedure        quarterTriangle1
            module procedure        quarterTriangle2
            module procedure        quarterTriangle3
        end interface

        interface   findVolume
            module procedure        findVolume0
        end interface

        interface   findArea
            module procedure        findArea0
            module procedure        findArea1
            module procedure        findArea2
        end interface

        interface   findNormals
            module procedure        findNormals0
            module procedure        findNormals1
        end interface


        interface   findVolumeAndArea
            module procedure        findVolumeAndArea0
            module procedure        findVolumeAndArea1
            module procedure        findVolumeAndArea2
            module procedure        findVolumeAndArea3
        end interface

        interface   getn_triangle
            module procedure        getn_triangle0
        end interface

        interface   getn_node
            module procedure        getn_node0
        end interface

    contains
!---^^^^^^^^


        function TriangulatedSurface_ctor2( that,n ) result(this)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      produce a new triangulated surface with n nodes/triangles
    !*      populate the first part with "that"
            type(TriangulatedSurface),intent(in)        ::      that
            integer,intent(in)                          ::      n
            type(TriangulatedSurface)                   ::      this
            integer         ::      ii
            this = TriangulatedSurface_ctor0(n)
            do ii = 1,min(that%N_array,this%N_array)
                this%node(:,ii) = that%node(:,ii)
                this%triangle(:,ii) = that%triangle(:,ii)
            end do
            this%n_node = that%n_node
            this%n_triangle = that%n_triangle
            call findBoundingSpheres(this)
            call generateNeighbourList(this)
            return
        end function TriangulatedSurface_ctor2


        function TriangulatedSurface_ctor1( node,tri,usershape ) result(this)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
            real(kind=real64),dimension(:,:),intent(in)              ::      node
            integer,dimension(:,:),intent(in)           ::      tri
            integer,intent(in),optional                 ::      usershape
            type(TriangulatedSurface)                   ::      this
            this = TriangulatedSurface_null( )
            this%n_node = size(node,dim=2)
            this%n_triangle = size(tri,dim=2) 
            this%N_array = max(this%n_node,this%n_triangle)
            allocate(this%node(3,this%N_array))
            allocate(this%triangle(3,this%N_array))
            this%node(1:3,1:this%n_node) = node
            this%triangle(1:3,1:this%n_triangle) = tri
            call findBoundingSpheres(this)
            call generateNeighbourList(this)
            if (present(usershape)) this%original = usershape
            return
        end function TriangulatedSurface_ctor1

        function TriangulatedSurface_ctor0( n ) result(this)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
            integer,intent(in)              ::      n
            type(TriangulatedSurface)                   ::      this
            this = TriangulatedSurface_null( )
            this%N_array = n
            allocate(this%node(3,n))
            allocate(this%triangle(3,n))
            this%node = 0.0
            this%triangle = 0
            this%N_neigh = 0
            return
        end function TriangulatedSurface_ctor0




        function TriangulatedSurface_null( ) result(this)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
            type(TriangulatedSurface)                   ::      this
            this%N_array = 0
            this%n_triangle = 0
            this%n_node = 0
            this%N_neigh = 0
            nullify(this%node)
            nullify(this%triangle)
            nullify(this%neigh)
            this%rmin2 = 0.0
            this%rmax2 = huge(1.0)
            this%original = USER
            return
        end function TriangulatedSurface_null

    !---


        function GeodesicGrid( n,startingshape ) result( this )
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      construct a geodesic grid with approximately n triangles
    !*      ( optionally starting from a truncated icosahedron )
            integer,intent(in)              ::      n
            integer,intent(in),optional     ::      startingshape
            type(TriangulatedSurface)       ::      this
            type(TriangulatedSurface)       ::      start
            integer     ::      nn,aa
            integer     ::      ii
            integer                 ::      i1,i2,i3,i4,i5
            integer                 ::      starting
            starting = -1
            if (present(startingshape)) starting = startingshape
            select case(starting)
                case (-1)   !   decide which start is the best
                    i1 = 20 * 4**(nint( log(n/20.0)/log(4.0) ))
                    i2 = 60 * 4**(nint( log(n/60.0)/log(4.0) ))
                    i3 = 180 * 4**(nint( log(n/180.0)/log(4.0) ))
                    i4 = 8 * 4**(nint( log(n/8.0)/log(4.0) ))
                    i5 = 12 * 4**(nint( log(n/12.0)/log(4.0) ))
                    print *,i1,i2,i3,i4,i5
                    if ( abs(i1-n) < minval( (/abs(i2-n),abs(i3-n),abs(i4-n),abs(i5-n)/) ) ) then
                        start = Icosahedron()
                        starting = ICOSAHEDRONSHAPE
                    else if ( abs(i2-n) < minval( (/abs(i1-n),abs(i3-n),abs(i4-n),abs(i5-n)/) )  ) then 
                        start = Dodecahedron()
                        starting = DODECAHEDRONSHAPE
                    else if ( abs(i3-n) < minval( (/abs(i1-n),abs(i2-n),abs(i4-n),abs(i5-n)/) )  ) then
                        start = TruncatedIcosahedron()
                        starting =TRUNCATEDICOSAHEDRONSHAPE
                    else if ( abs(i4-n) < minval( (/abs(i1-n),abs(i2-n),abs(i3-n),abs(i5-n)/) )  ) then
                        start = Octahedron()
                        starting = OCTAHEDRONSHAPE
                    else
                        start = Cube()
                        starting = CUBESHAPE
                    end if
                case(ICOSAHEDRONSHAPE)
                    start = Icosahedron()
                case(DODECAHEDRONSHAPE)
                    start = Dodecahedron()
                case(TRUNCATEDICOSAHEDRONSHAPE)
                    start = TruncatedIcosahedron()
                case(OCTAHEDRONSHAPE)
                    start = Octahedron()
                case(CUBESHAPE)
                    start = Cube()
                case default
                    stop "Lib_TriangulatedSurface::GeodesicGrid error - unknown starting shape"
            end select
        !           n ~ 20 * 4^a
        !       log (n/20) / log (4) ~ a
            aa = nint( log(n/real(start%n_triangle))/log(4.0) )
            nn = start%n_triangle * 4**aa
            this = TriangulatedSurface_ctor0( nn )
            this%original = starting


            this%n_node = start%n_node
            this%n_triangle = start%n_triangle
            this%node(1:3,1:this%n_node) = start%node(1:3,1:this%n_node)
            this%triangle(1:3,1:this%n_triangle) = start%triangle(1:3,1:start%n_triangle)
!           print *,aa,nn
            do ii = 1,aa
                call refineMesh4( this )
                call inflate(this)
            end do
            if (aa==0) call inflate(this)

            call delete(start)
            call findBoundingSpheres(this)
            call generateNeighbourList(this)
            this%N_array = max(this%n_node,this%n_triangle)
            return

        contains
    !---^^^^^^^^

            pure subroutine inflate(this)
        !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
                type(TriangulatedSurface),intent(inout)     ::      this
                integer             ::      ii
                real(kind=real64),dimension(3)   ::      xx
                real(kind=real64)                ::      dd
                do ii = 1,this%n_node
                    xx = this%node(:,ii)
                    dd = 1.0/sqrt( xx(1)*xx(1) + xx(2)*xx(2) + xx(3)*xx(3) )
                    this%node(:,ii) = this%node(:,ii)*dd
                end do
                return
            end subroutine inflate

        end function GeodesicGrid

!-------

        subroutine reallocMemory0( this,n )
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      increase size of the storage in this surface to accomodate n triangles
            type(TriangulatedSurface),intent(inout)         ::      this
            integer,intent(in)                              ::      n

            real(kind=real64),dimension(:,:),pointer         ::  ttp
            integer,dimension(:,:),pointer      ::  tip,tnp
            if (this%N_array >= n) return
            this%N_array = n
            allocate(ttp(3,this%N_array))
            allocate(tip(3,this%N_array))
            ttp(:,1:this%n_node) = this%node(:,1:this%n_node)
            tip(:,1:this%n_triangle) = this%triangle(:,1:this%n_triangle)
            deallocate(this%node)
            deallocate(this%triangle)
            this%node => ttp
            this%triangle => tip
            if (hasNeighbourList(this)) then
                allocate(tnp(0:this%N_neigh,this%N_array))
                tnp(:,1:this%n_node) = this%neigh(:,1:this%n_node)
                deallocate(this%neigh)
                this%neigh => tnp
            end if
            return
        end subroutine reallocMemory0


!-------

        subroutine refineMesh4( this )
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*  a triangle is easily split into 4 by halving the side length
    !*              1                           1
    !*             / \                         / \
    !*            /   \             ->        4   5
    !*           /     \                     /     \
    !*          2-------3                   2---6---3
            type(TriangulatedSurface),intent(inout)     ::      this
!           integer             ::      i1,i2,i3,i4,i5,i6
            integer             ::      jj
!           real(kind=real64),dimension(3)   ::      x1,x2,x3,x4,x5,x6
            call reallocMemory( this,this%n_triangle*4 )

            do jj = 1,this%n_triangle
                call quarterTriangle( this,jj )
!               i1 = this%triangle(1,jj)
!               i2 = this%triangle(2,jj)
!               i3 = this%triangle(3,jj)
!               x1 = this%node(:,i1)
!               x2 = this%node(:,i2)
!               x3 = this%node(:,i3)
!               x4 = 0.5*(x1+x2)
!               x5 = 0.5*(x1+x3)
!               x6 = 0.5*(x2+x3)
!               call addTriangle( this,x4,x2,x6 )
!               call addTriangle( this,x4,x6,x5 )
!               i4 = this%triangle(1,this%n_triangle)
!               i5 = this%triangle(3,this%n_triangle)
!               i6 = this%triangle(2,this%n_triangle)
!               call addTriangle( this,x5,x6,x3 )
!               this%triangle( :,jj ) = (/ i1,i4,i5 /)
            end do
            return
        end subroutine refineMesh4

        subroutine quarterTriangle0( this,i )
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      assumes sufficient memory allocated
    !*              1              1
    !*             / \            / \
    !*            /   \          4---5
    !*           /     \        / \ / \
    !*          2-------3      2---6---3
            type(TriangulatedSurface),intent(inout)     ::      this
            integer,intent(in)                          ::      i
            integer             ::      i1,i2,i3,i4,i5
            real(kind=real64),dimension(3)   ::      x1,x2,x3,x4,x5,x6
            i1 = this%triangle(1,i)
            i2 = this%triangle(2,i)
            i3 = this%triangle(3,i)
            x1 = this%node(:,i1)
            x2 = this%node(:,i2)
            x3 = this%node(:,i3)
            x4 = 0.5*(x1+x2)
            x5 = 0.5*(x1+x3)
            x6 = 0.5*(x2+x3)
            call addTriangle( this,x4,x2,x6 )
            call addTriangle( this,x4,x6,x5 )
            i4 = this%triangle(1,this%n_triangle)
            i5 = this%triangle(3,this%n_triangle)
            call addTriangle( this,x5,x6,x3 )
            this%triangle( :,i ) = (/ i1,i4,i5 /)
            return
        end subroutine quarterTriangle0


        subroutine quarterTriangle1( this,n,quarterList )
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      assumes sufficient memory allocated
            type(TriangulatedSurface),intent(inout)     ::      this
            integer,intent(in)                          ::      n
            integer,dimension(:),intent(in)             ::      quarterList
            integer     ::      ii,jj,kk,ll
            integer     ::      i1,i2,i3
            integer     ::      nhalf
            integer,dimension(3,this%n_triangle)        ::  oldTriangle
            integer,dimension(2,this%n_triangle)        ::  halfList
            oldTriangle = this%triangle(1:3,1:this%n_triangle)
        !---    first quarter all the triangles
            do ii = 1,n
                jj = quarterList(ii)
!               xx = this%node(:,oldTriangle(1,jj)) + this%node(:,oldTriangle(2,jj)) + this%node(:,oldTriangle(3,jj))
!               if (all(xx>0)) &
!               write (*,fmt='(a,5i8)') "quarter ",ii,jj,oldTriangle(:,jj)
                call quarterTriangle0(this,jj)
            end do
        !---    next check whether this leaves a node half way across any triangle
            halfList = 0
            nHalf = 0
            do ii = 1,n
                jj = quarterList(ii)
                i1 = oldTriangle(1,jj)
                i2 = oldTriangle(2,jj)
                i3 = oldTriangle(3,jj)
!               xx = this%node(:,oldTriangle(1,jj)) + this%node(:,oldTriangle(2,jj)) + this%node(:,oldTriangle(3,jj))
                do kk = 1,this%n_triangle
                    if (kk==jj) cycle
!                     is = inHalfList( halfList(:,1:nHalf),kk )
                    ll = triangleSharesSide( this,kk,i1,i2 )
                    if (ll /= 0) then
                        nHalf = nHalf + 1
                        halfList(:,nHalf) = (/ kk,ll /)
!               if (all(xx>0)) &
!                       write (*,fmt='(a,7i8)') "half ",ii,jj,oldTriangle(:,jj),kk,ll
                    end if
                    ll = triangleSharesSide( this,kk,i2,i3 )
                    if (ll /= 0) then
                        nHalf = nHalf + 1
                        halfList(:,nHalf) = (/ kk,ll /)
!               if (all(xx>0)) &
!                       write (*,fmt='(a,7i8)') "half ",ii,jj,oldTriangle(:,jj),kk,ll
                    end if
                    ll = triangleSharesSide( this,kk,i3,i1 )
                    if (ll /= 0) then
                        nHalf = nHalf + 1
                        halfList(:,nHalf) = (/ kk,ll /)
!               if (all(xx>0)) &
!                       write (*,fmt='(a,7i8)') "half ",ii,jj,oldTriangle(:,jj),kk,ll
                    end if
                end do
            end do

!           print *,"quarterTriangle1 ",this%n_node,this%n_triangle,this%N_array
            do ii = 1,nHalf
                call halfTriangle(this,halfList(1,ii),halfList(2,ii))
            end do
!             print *,"quarterTriangle1 ",this%n_node,this%n_triangle,this%N_array

            call generateNeighbourList(this)

            return

        contains
    !---^^^^^^^^

            pure function inHalfList( halfList,kk ) result(is)
        !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
        !*      returns true if triangle kk is in the halflist ( do not bother about sides )
                integer,dimension(:,:),intent(in)       ::      halflist
                integer,intent(in)                      ::      kk
                logical                                 ::      is
                integer     ::      ii
                is = .false.
                do ii = 1,size(halflist,dim=2)
                    if (kk == halflist(2,ii)) then
                        is = .true.
                        exit
                    end if
                end do
                return
            end function inHalfList


        end subroutine quarterTriangle1


        subroutine quarterTriangle2( this )
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      assumes sufficient memory allocated. Quarters all triangles
            type(TriangulatedSurface),intent(inout)     ::      this
            integer     ::      ii,nn
        !---    quarter all the triangles
            nn = this%n_triangle
            do ii = 1,nn
                call quarterTriangle0(this,ii)
            end do
            call generateNeighbourList(this)

            return
        end subroutine quarterTriangle2


        subroutine quarterTriangle3( this,complete )
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      Repeatedly quarters all triangles until memory is filled
            type(TriangulatedSurface),intent(inout)     ::      this
            logical,intent(in)                          ::      complete
            integer     ::      ii,nn,mm
            do
            !---    quarter all the triangles
                nn = this%n_triangle
                mm = this%n_node
                if (max(mm,nn)*5 > this%N_array) exit       !   4x for periodic repeat, 5x for safety.
                do ii = 1,nn
                    call quarterTriangle0(this,ii)
                end do
                call generateNeighbourList(this)
            end do
            return
        end subroutine quarterTriangle3



        subroutine thirdTriangle0( this,i )
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      assumes sufficient memory allocated
            type(TriangulatedSurface),intent(inout)     ::      this
            integer,intent(in)                          ::      i
            integer             ::      i1,i2,i3,i4
            real(kind=real64),dimension(3)   ::      x1,x2,x3,x4
            i1 = this%triangle(1,i)
            i2 = this%triangle(2,i)
            i3 = this%triangle(3,i)
            x1 = this%node(:,i1)
            x2 = this%node(:,i2)
            x3 = this%node(:,i3)
            x4 = 0.3333333333333333*(x1+x2+x3)
            call addTriangle( this,x1,x2,x4 )
            i4 = this%triangle(3,this%n_triangle)
            call addTriangle( this,x4,x2,x3 )
            this%triangle( :,i ) = (/ i1,i4,i3 /)
            return
        end subroutine thirdTriangle0


        subroutine thirdTriangle1( this,n,thirdList )
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      assumes sufficient memory allocated
            type(TriangulatedSurface),intent(inout)     ::      this
            integer,intent(in)                          ::      n
            integer,dimension(:),intent(in)             ::      thirdList
            integer     ::      ii,jj
        !---    first third all the triangles
            do ii = 1,n
                jj = thirdList(ii)
                call thirdTriangle0(this,jj)
            end do
        !---    update neighbour list
            call generateNeighbourList(this)

            return
        end subroutine thirdTriangle1



        function triangleSharesSide( this,j,i1,i2 ) result(k)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      returns non-zero if triangle j shares nodes i1,i2
    !*      k = 1 if i1-i2 are side 1-2
    !*      k = 2 if i1-i2 are side 2-3
    !*      k = 3 if i1-i2 are side 3-1
            type(TriangulatedSurface),intent(inout)     ::      this
            integer,intent(in)                          ::      j
            integer,intent(in)                          ::      i1,i2
            integer                                     ::      k
            k = 0
            if (this%triangle(1,j)==i1) then
                if (this%triangle(2,j)==i2) then
                    k = 1
                else if (this%triangle(3,j)==i2) then
                    k = 3
                end if
            else if (this%triangle(2,j)==i1) then
                if (this%triangle(1,j)==i2) then
                    k = 1
                else if (this%triangle(3,j)==i2) then
                    k = 2
                end if
            else if (this%triangle(3,j)==i1) then
                if (this%triangle(2,j)==i2) then
                    k = 2
                else if (this%triangle(1,j)==i2) then
                    k = 3
                end if
            end if
            return
        end function triangleSharesSide

        subroutine halfTriangle( this,i,j )
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      assumes sufficient memory allocated
    !*      j = 1,split side 1-2                        1
    !*      j = 2,split side 2-3                       / \
    !*      j = 3,split side 3-1                      2---3

            type(TriangulatedSurface),intent(inout)     ::      this
            integer,intent(in)                          ::      i,j
            integer             ::      i1,i2,i3,i4
            real(kind=real64),dimension(3)   ::      x1,x2,x3,x4
            i1 = this%triangle(1,i)
            i2 = this%triangle(2,i)
            i3 = this%triangle(3,i)
            x1 = this%node(:,i1)
            x2 = this%node(:,i2)
            x3 = this%node(:,i3)
            select case (j)
                case(1)
                    x4 = 0.5*(x1+x2)
                    call addTriangle( this,x1,x4,x3 )
                    i4 = this%triangle(2,this%n_triangle)
                    this%triangle( :,i ) = (/ i4,i2,i3 /)
                case(2)
                    x4 = 0.5*(x2+x3)
                    call addTriangle( this,x1,x2,x4 )
                    i4 = this%triangle(3,this%n_triangle)
                    this%triangle( :,i ) = (/ i1,i4,i3 /)
                case(3)
                    x4 = 0.5*(x3+x1)
                    call addTriangle( this,x1,x2,x4 )
                    i4 = this%triangle(3,this%n_triangle)
                    this%triangle( :,i ) = (/ i4,i2,i3 /)
            end select
            return
        end subroutine halfTriangle
!-------

        subroutine addTriangle(this,x1,x2,x3)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
            type(TriangulatedSurface),intent(inout)     ::      this
            real(kind=real64),dimension(3),intent(in)                ::      x1,x2,x3
            integer                 ::      ii
            integer                 ::      i1,i2,i3
            real(kind=real64),dimension(3)       ::      xx
            real(kind=real64)                    ::      dd
            i1 = 0 ; i2 = 0 ; i3 = 0

            do ii = 1,this%n_node
                xx = x1 - this%node(:,ii)
                dd = xx(1)*xx(1) + xx(2)*xx(2) + xx(3)*xx(3)
                if (dd < 1.0d-8) then
                    i1 = ii
                    exit
                end if
            end do
            do ii = 1,this%n_node
                xx = x2 - this%node(:,ii)
                dd = xx(1)*xx(1) + xx(2)*xx(2) + xx(3)*xx(3)
                if (dd < 1.0d-8) then
                    i2 = ii
                    exit
                end if
            end do
            do ii = 1,this%n_node
                xx = x3 - this%node(:,ii)
                dd = xx(1)*xx(1) + xx(2)*xx(2) + xx(3)*xx(3)
                if (dd < 1.0d-8) then
                    i3 = ii
                    exit
                end if
            end do

            if (i1 == 0) then
                this%n_node = this%n_node + 1
                i1 = this%n_node
                this%node(:,i1) = x1
!               if (all(x1>0)) &
!               write (*,fmt='(a,i8,3f12.5)') "add node ",i1,x1
            end if
            if (i2 == 0) then
                this%n_node = this%n_node + 1
                i2 = this%n_node
                this%node(:,i2) = x2
!               if (all(x2>0)) &
!               write (*,fmt='(a,i8,3f12.5)') "add node ",i2,x2
            end if
            if (i3 == 0) then
                this%n_node = this%n_node + 1
                i3 = this%n_node
                this%node(:,i3) = x3
!               if (all(x3>0)) &
!               write (*,fmt='(a,i8,3f12.5)') "add node ",i3,x3
            end if

            this%n_triangle = this%n_triangle + 1
            this%triangle(:,this%n_triangle) = (/ i1,i2,i3 /)
            return
        end subroutine addTriangle

!-------

        function Icosahedron() result( this )
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      construct a icosahedron
            type(TriangulatedSurface)       ::      this

            real(kind=real64),parameter      ::      SIDE = 1.0514622242382672120513381696958d0    !   side length required for sphere radius 1
            real(kind=real64),parameter      ::      GOLD = 1.6180339887498948482045868343656d0    !   golden ratio
            real(kind=real64),dimension(3,12)    ::      vertex =        &
                    reshape( (/ 0.d0,1.d0,GOLD , 0.d0,1.d0,-GOLD , 0.d0,-1.d0,GOLD , 0.d0,-1.d0,-GOLD ,     &
                                1.d0,GOLD,0.d0 , 1.d0,-GOLD,0.d0 , -1.d0,GOLD,0.d0 , -1.d0,-GOLD,0.d0 ,     &
                                GOLD,0.d0,1.d0 , GOLD,0.d0,-1.d0 , -GOLD,0.d0,1.d0 , -GOLD,0.d0,-1.d0 /)    &
                                , (/ 3,12 /) )*SIDE/2

            integer,dimension(3,20) ::      triangle =      &
                    reshape( (/  1, 3 , 9 ,   1, 11, 3,    1, 5 , 7 ,   1, 9 , 5 ,                 &
                                 1, 7 , 11,   2, 10, 4,    2, 4 , 12,   2, 7 , 5 ,                 &
                                 2, 5 , 10,   2, 12, 7,    3, 8 , 6 ,   3, 6 , 9 ,                 &
                                 3, 11, 8,    4, 6 , 8 ,   4, 10, 6,    4, 8 , 12,                 &
                                 5, 9 , 10,   6, 10, 9,    7, 12, 11,   8, 11, 12     /) , (/ 3,20 /) )

            this = TriangulatedSurface_ctor(20)
            this%node(1:3,1:12) = vertex
            this%triangle(1:3,1:20) = triangle
            this%n_node = 12
            this%n_triangle = 20
            this%N_array = max(this%n_node,this%n_triangle)
            return
        end function Icosahedron


        function Dodecahedron() result( this )
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      construct a dodecahedron
            type(TriangulatedSurface)       ::      this

            real(kind=real64),parameter      ::      SIDE = 0.57735026918962576450914878050196d0   !   side length required for sphere radius 1
            real(kind=real64),parameter      ::      GOLD = 1.6180339887498948482045868343656d0    !   golden ratio
            real(kind=real64),parameter      ::      IGOL = 1.d0/GOLD                             !   (inverse) golden ratio
            real(kind=real64),dimension(3,32)    ::      vertex =        &
                    reshape( (/     1.d0,  1.d0,  1.d0  ,  1.d0,  1.d0, -1.d0  ,  1.d0, -1.d0,  1.d0  ,  1.d0, -1.d0, -1.d0  ,  &
                                   -1.d0,  1.d0,  1.d0  , -1.d0,  1.d0, -1.d0  , -1.d0, -1.d0,  1.d0  , -1.d0, -1.d0, -1.d0  ,  &
                                    0.d0, IGOL, GOLD  ,  0.d0, IGOL,-GOLD  ,  0.d0,-IGOL, GOLD  ,  0.d0,-IGOL,-GOLD  ,  &
                                   IGOL, GOLD,  0.d0  , IGOL,-GOLD,  0.d0  ,-IGOL, GOLD,  0.d0  ,-IGOL,-GOLD,  0.d0  ,  &
                                   GOLD,  0.d0, IGOL  , GOLD,  0.d0,-IGOL  ,-GOLD,  0.d0, IGOL  ,-GOLD,  0.d0,-IGOL  ,  &
                                   0.d0,0.d0,0.d0 , 0.d0,0.d0,0.d0 , 0.d0,0.d0,0.d0 , 0.d0,0.d0,0.d0 ,  &
                                   0.d0,0.d0,0.d0 , 0.d0,0.d0,0.d0 , 0.d0,0.d0,0.d0 , 0.d0,0.d0,0.d0 ,  &
                                   0.d0,0.d0,0.d0 , 0.d0,0.d0,0.d0 , 0.d0,0.d0,0.d0 , 0.d0,0.d0,0.d0    &
                                /), (/ 3,32 /) )*SIDE
            integer,dimension(5,12) ::      pentagon = reshape( (/  1, 9, 5,15,13   ,           &
                                                                    1, 9,11, 3,17   ,           &
                                                                    9, 5,19, 7,11   ,           &
                                                                    5,15, 6,20,19   ,           &
                                                                   15,13, 2,10, 6   ,           &
                                                                   13, 1,17,18, 2   ,           &
                                                                   14,16, 8,12, 4   ,           &
                                                                   14,16, 7,11, 3   ,           &
                                                                   16, 8,20,19, 7   ,           &
                                                                    8,12,10, 6,20   ,           &
                                                                   12, 4,18, 2,10   ,           &
                                                                    4,14, 3,17,18   /),(/5,12/) )


            integer,dimension(3,60) ::      triangle =      &
                    reshape( (/      1,21, 9,    1, 9,22,    1,13,21,    1,26,13,       &
                                     1,22,17,    1,17,26,    2,10,25,    2,31,10,       &
                                     2,25,13,    2,13,26,    2,26,18,    2,18,31,       &
                                     3,22,11,    3,11,28,    3,28,14,    3,14,32,       &
                                     3,17,22,    3,32,17,    4,27,12,    4,12,31,       &
                                     4,14,27,    4,32,14,    4,31,18,    4,18,32,       &
                                     5, 9,21,    5,23, 9,    5,21,15,    5,15,24,       &
                                     5,19,23,    5,24,19,    6,25,10,    6,10,30,       &
                                     6,24,15,    6,15,25,    6,20,24,    6,30,20,       &
                                     7,11,23,    7,28,11,    7,16,28,    7,29,16,       &
                                     7,23,19,    7,19,29,    8,12,27,    8,30,12,       &
                                     8,27,16,    8,16,29,    8,29,20,    8,20,30,       &
                                     9,11,22,    9,23,11,   10,12,30,   10,31,12,       &
                                    13,15,21,   13,25,15,   14,16,27,   14,28,16,       &
                                    17,18,26,   17,32,18,   19,24,20,   19,20,29        &
                        /) , (/ 3,60 /) )
            integer     ::      ii,jj
            real(kind=real64),dimension(3)   ::      xx
            this = TriangulatedSurface_ctor(60)
            do jj = 1,12
                xx = 0.0
                do ii = 1,5
                    xx = xx + vertex(:,pentagon(ii,jj))
                end do
                vertex(:,20+jj) = xx*0.2
            end do
            this%node(1:3,1:32) = vertex
            this%triangle(1:3,1:60) = triangle
            this%n_node = 32
            this%n_triangle = 60
            this%N_array = max(this%n_node,this%n_triangle)
            return
        end function Dodecahedron



        function TruncatedIcosahedron() result( this )
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      construct a truncated icosahedron (buckyball)
    !*      a buckyball has 20 hexagonal faces and 12 pentagons
    !*      - to triangulate gives 180 triangles and 92 vertices
    !*      to construct add nodes 1/3 and 2/3 along each edge
    !*                 1                          1
    !*                / \                        / \
    !*               /   \          ->          4---5
    !*              /     \                    / \ / \
    !*             /       \                  6--10---7
    !*            /         \                / \ / \ / \
    !*           2-----------3              2---8---9---3
    !*      ... turning each triangle into 9
    !*      Then bring back nodes 1,2,3 into the pentagonal planes

            type(TriangulatedSurface)       ::      this
            type(TriangulatedSurface)       ::      icos
            real(kind=real64),dimension(3)       ::      x1,x2,x3,x4,x5,x6,x7,x8,x9,x10
            integer                 ::      i1,i2,i3,i4,i5,i6,i7,i8,i9
            integer,dimension(5)    ::      pentagon
            integer                 ::      ii,np
            real(kind=real64)                    ::      dd
            icos = Icosahedron()
            this = TriangulatedSurface_ctor(180)
        !---    put in all the vertices
            this%n_node = 12
            this%node(:,1:12) = icos%node(:,1:12)
            this%n_triangle = 20
            this%triangle(:,1:20) = icos%triangle(:,1:20)

        !---    increase the number of triangles 9-fold
            np = 0                      !   number of points on the pentagon surrounding vertex 1 found
            pentagon = 0
            do ii = 1,icos%n_triangle
                i1 = icos%triangle(1,ii)
                i2 = icos%triangle(2,ii)
                i3 = icos%triangle(3,ii)


                x1 = icos%node(:,i1)
                x2 = icos%node(:,i2)
                x3 = icos%node(:,i3)
                x4 = (2*x1 + x2)/3
                x5 = (2*x1 + x3)/3
                x6 = (2*x2 + x1)/3
                x7 = (2*x3 + x1)/3
                x8 = (2*x2 + x3)/3
                x9 = (2*x3 + x2)/3
                x10 = (x1 + x2 + x3)/3

                call addTriangle( this,x4,x6,x10 )
                call addTriangle( this,x4,x10,x5 )
                i4 = this%triangle(1,this%n_triangle)
                i5 = this%triangle(3,this%n_triangle)
                call addTriangle( this,x5,x10,x7 )
                call addTriangle( this,x6,x8,x10 )
                call addTriangle( this,x6,x2,x8 )
                i6 = this%triangle(1,this%n_triangle)
                i8 = this%triangle(3,this%n_triangle)
                call addTriangle( this,x10,x8,x9 )
                call addTriangle( this,x7,x10,x9 )
                i7 = this%triangle(1,this%n_triangle)
                i9 = this%triangle(3,this%n_triangle)
                call addTriangle( this,x7,x9,x3 )
                this%triangle( :,ii ) = (/ i1,i4,i5 /)

            !---    look out for the 5 triangles surrounding original vertex 1
                if (i1 == 1) then
                    if (.not. any(pentagon==i4)) then
                        np = np + 1
                        pentagon(np) = i4
                    end if
                    if (.not. any(pentagon==i5)) then
                        np = np + 1
                        pentagon(np) = i5
                    end if
                else if (i2 == 1) then
                    if (.not. any(pentagon==i6)) then
                        np = np + 1
                        pentagon(np) = i6
                    end if
                    if (.not. any(pentagon==i8)) then
                        np = np + 1
                        pentagon(np) = i8
                    end if
                else if (i3 == 1) then
                    if (.not. any(pentagon==i7)) then
                        np = np + 1
                        pentagon(np) = i7
                    end if
                    if (.not. any(pentagon==i9)) then
                        np = np + 1
                        pentagon(np) = i9
                    end if
                end if

            end do

        !---    put back the original vertices
            x1 = this%node(:,pentagon(1))
            x2 = this%node(:,pentagon(2))
            x3 = this%node(:,pentagon(3))
            x4 = this%node(:,pentagon(4))
            x5 = this%node(:,pentagon(5))
            x6 = (x1+x2+x3+x4+x5)*0.2
            dd = x6(1)*x6(1) + x6(2)*x6(2) + x6(3)*x6(3)    !   correct distance squared to pentagon face
            dd = sqrt(dd)
            do ii = 1,12
                this%node(:,ii) = this%node(:,ii)*dd        !   note: icosahedron has vertices at radius 1
            end do

        !---    now inflate vertices to radius 1
            dd = x1(1)*x1(1) + x1(2)*x1(2) + x1(3)*x1(3)    !   distance squared to a new vertex
            dd = 1/sqrt(dd)
            do ii = 1,92
                this%node(:,ii) = this%node(:,ii)*dd
            end do
            this%N_array = max(this%n_node,this%n_triangle)
            return
        end function TruncatedIcosahedron


        function Octahedron() result( this )
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      construct an octahedron
            type(TriangulatedSurface)       ::      this

            real(kind=real64),dimension(3,6) ::      vertex =        &
                    reshape( (/     0.0,0.0,1.0  ,  1.0,0.0,0.0 , 0.0,1.0,0.0 ,         &
                                    -1.,0.0,0.0  ,  0.0,-1.,0.0 , 0.0,0.0,-1.           &
                                /), (/ 3,6 /) )


            integer,dimension(3,8)  ::      triangle =      &
                    reshape( (/      1,2,3  ,  1,3,4  , 1,4,5  ,  1,5,2  ,  &
                                     2,6,3  ,  3,6,4  , 4,6,5  ,  5,6,2     &
                        /) , (/ 3,8 /) )
            this = TriangulatedSurface_ctor(8)
            this%node(1:3,1:6) = vertex
            this%triangle(1:3,1:8) = triangle
            this%n_node = 6
            this%n_triangle = 8
            this%N_array = max(this%n_node,this%n_triangle)
            return
        end function Octahedron

        function Cube() result( this )
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      construct a Cube
            type(TriangulatedSurface)       ::      this
            real(kind=real64),parameter      ::      SIDE = 0.57735026918962576450914878050196
            real(kind=real64),dimension(3,8) ::      vertex =        &
                    reshape( (/     1,1,1 , 1,1,-1 , 1,-1,1 , 1,-1,-1   ,   &
                                   -1,1,1 ,-1,1,-1 ,-1,-1,1 ,-1,-1,-1       &
                              /)*SIDE, (/ 3,8 /) )


            integer,dimension(3,12) ::      triangle =      &
                    reshape( (/     3,2,1, 4,2,3, 3,7,4, 4,7,8      ,   &
                                    1,5,7, 1,7,3, 7,5,6, 7,6,8      ,   &
                                    5,1,2, 5,2,6, 6,2,8, 2,4,8          &
                        /) , (/ 3,12 /) )
            this = TriangulatedSurface_ctor(12)
            this%node(1:3,1:8) = vertex
            this%triangle(1:3,1:12) = triangle
            this%n_node = 8
            this%n_triangle = 12
            this%N_array = max(this%n_node,this%n_triangle)
            return
        end function Cube


!-------

        subroutine report1(this,u)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^
            type(TriangulatedSurface),intent(in)            ::      this
            integer,intent(in),optional     ::      u
            integer         ::      uu
            uu = 6
            if (present(u)) uu = u
            select case(this%original)
                case(ICOSAHEDRONSHAPE)
                    write (unit=uu,fmt='(a)')        "TriangulatedSurface   (icosahedron)"
                case(DODECAHEDRONSHAPE)
                    write (unit=uu,fmt='(a)')        "TriangulatedSurface   (dodecahedron)"
                case(TRUNCATEDICOSAHEDRONSHAPE)
                    write (unit=uu,fmt='(a)')        "TriangulatedSurface   (truncated icosahedron)"
                case(OCTAHEDRONSHAPE)
                    write (unit=uu,fmt='(a)')        "TriangulatedSurface   (octahedron)"
                case(CUBESHAPE)
                    write (unit=uu,fmt='(a)')        "TriangulatedSurface   (cube)"
                case(WIGNERCELLSHAPE)
                    write (unit=uu,fmt='(a)')        "TriangulatedSurface   (WignerCell)"
                case(USER)
                    write (unit=uu,fmt='(a)')        "TriangulatedSurface   (user)"
                case default
                    write (unit=uu,fmt='(a)')        "TriangulatedSurface   "
            end select
            write (unit=uu,fmt='(a,2i8)')    "    nodes,triangles : ",this%n_node,this%n_triangle
            return
        end subroutine report1

        subroutine report2(this,u,advfile)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
            type(TriangulatedSurface),intent(in)    ::      this
            integer,intent(in)                      ::      u
            logical,intent(in)                      ::      advfile
            integer                 ::      ii
            real(kind=real64),dimension(3,3)     ::      tri

            do ii = 1,this%n_triangle
                tri = getTriangle( this,ii )
                write (unit=u,fmt='(a5,100f16.8)')   "tri  ",                                &
                                tri(1:3,1), tri(1:3,2), tri(1:3,3),  tri(3,1:3)
            end do
            return
        end subroutine report2



        subroutine delete1(this)
    !---^^^^^^^^^^^^^^^^^^^^^^^^
            type(TriangulatedSurface),intent(inout)         ::      this
            if (this%n_array == 0) return
            deallocate(this%node)
            deallocate(this%triangle)
            if (this%N_neigh > 0) deallocate(this%neigh)
            this = TriangulatedSurface_null( )
            return
        end subroutine delete1
        


!-------

        pure function getn_triangle0(this) result(n)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
            type(TriangulatedSurface),intent(in)    ::      this
            integer                                 ::      n
            n = this%n_triangle
            return
        end function getn_triangle0

        pure function getn_node0(this) result(n)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
            type(TriangulatedSurface),intent(in)    ::      this
            integer                                 ::      n
            n = this%n_node
            return
        end function getn_node0

        pure function getN_neigh0(this,i) result(n)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      returns the number of triangles to which node i belongs
            type(TriangulatedSurface),intent(in)    ::      this
            integer,intent(in)                      ::      i
            integer                                 ::      n
            n = this%neigh(0,i)
            return
        end function getN_neigh0

        pure function getNeigh0(this,j,i) result(n)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      returns the number of the jth triangle to which node i belongs
            type(TriangulatedSurface),intent(in)    ::      this
            integer,intent(in)                      ::      i
            integer,intent(in)                      ::      j
            integer                                 ::      n
            n = this%neigh(j,i)
            return
        end function getNeigh0


    !---

        pure function getTriangle(this,i) result( tri )
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
            type(TriangulatedSurface),intent(in)        ::      this
            integer,intent(in)                          ::      i
            real(kind=real64),dimension(3,3)                         ::      tri
            tri(:,1) = this%node( :,this%triangle(1,i) )
            tri(:,2) = this%node( :,this%triangle(2,i) )
            tri(:,3) = this%node( :,this%triangle(3,i) )
            return
        end function getTriangle

        pure subroutine getTriangleInfo(this,i,tri,n )
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
            type(TriangulatedSurface),intent(in)        ::      this
            integer,intent(in)                          ::      i
            real(kind=real64),dimension(3,3),intent(out)             ::      tri
            integer,dimension(3),intent(out)            ::      n
            n = this%triangle(:,i)
            tri(:,1) = this%node( :,n(1) )
            tri(:,2) = this%node( :,n(2) )
            tri(:,3) = this%node( :,n(3) )
            return
        end subroutine getTriangleInfo

        pure function getTriangleNodes(this,i) result(n)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
            type(TriangulatedSurface),intent(in)        ::      this
            integer,intent(in)                          ::      i
            integer,dimension(3)                        ::      n
            n = this%triangle(:,i)
            return
        end function getTriangleNodes

        pure function getNode(this,i) result( x )
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
            type(TriangulatedSurface),intent(in)        ::      this
            integer,intent(in)                          ::      i
            real(kind=real64),dimension(3)                           ::      x
            x(:) = this%node(:,i)
            return
        end function getNode

        pure subroutine setNode(this,i,x )
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
            type(TriangulatedSurface),intent(inout)     ::      this
            integer,intent(in)                          ::      i
            real(kind=real64),dimension(3),intent(in)                ::      x
            this%node(:,i) = x(:)
            return
        end subroutine setNode



    !---


        pure function crossProduct(a,b) result(c)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      ronseal
            real(kind=real64),dimension(3),intent(in)        ::      a,b
            real(kind=real64),dimension(3)                   ::      c
            c(1) = a(2)*b(3) - a(3)*b(2)
            c(2) = a(3)*b(1) - a(1)*b(3)
            c(3) = a(1)*b(2) - a(2)*b(1)
            return
        end function crossProduct


        pure function getArea1(this,ii,x0) result(a)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      returns (signed) area of triangle ii
            type(TriangulatedSurface),intent(in)    ::      this
            integer,intent(in)                      ::      ii
            real(kind=real64),dimension(3),intent(in),optional   ::      x0
            real(kind=real64)                                ::      a
            real(kind=real64),dimension(3)       ::      nn,cc
        !---    area of a triangle = 1/2 |a x b|
        !---    use nn,cc as dummy storage
            cc = this%node( 1:3,this%triangle(2,ii) )-this%node( 1:3,this%triangle(1,ii) )
            nn = this%node( 1:3,this%triangle(3,ii) )-this%node( 1:3,this%triangle(1,ii) )
        !-- nn  = (x2-x1) x (x3-x1)
            nn = crossProduct( cc,nn )
            a = nn(1)*nn(1) + nn(2)*nn(2) + nn(3)*nn(3)
            if (a > 1.0d-16) then
                a = 0.5*sqrt(a)
            else
                a = 0.0
            end if
            if (present(x0)) then
                cc = ( this%node( 1:3,this%triangle(1,ii) )         &
                    + this%node( 1:3,this%triangle(2,ii) )      &
                    + this%node( 1:3,this%triangle(3,ii) ) )*0.3333333333333333
                if ( cc(1)*nn(1) + cc(2)*nn(2) + cc(3)*nn(3)            &
                    < x0(1)*nn(1) + x0(2)*nn(2) + x0(3)*nn(3) ) a = -a
            end if
            return
        end function getArea1


        pure function getNorm(this,ii,x0) result(n)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      returns outward normal of triangle ii
            type(TriangulatedSurface),intent(in)    ::      this
            integer,intent(in)                      ::      ii
            real(kind=real64),dimension(3),intent(in),optional   ::      x0
            real(kind=real64),dimension(3)       ::      n
            real(kind=real64),dimension(3)       ::      cc
            real(kind=real64)                    ::      dd
        !---    area of a triangle = 1/2 |a x b|
        !---    use nn,cc as dummy storage
            cc = this%node( 1:3,this%triangle(2,ii) )-this%node( 1:3,this%triangle(1,ii) )
            n = this%node( 1:3,this%triangle(3,ii) )-this%node( 1:3,this%triangle(1,ii) )
        !-- n  = (x2-x1) x (x3-x1)
            n = crossProduct( cc,n )
            if (present(x0)) then
                cc = ( this%node( 1:3,this%triangle(1,ii) )         &
                    + this%node( 1:3,this%triangle(2,ii) )      &
                    + this%node( 1:3,this%triangle(3,ii) ) )*0.3333333333333333
                if ( cc(1)*n(1) + cc(2)*n(2) + cc(3)*n(3)           &
                    < x0(1)*n(1) + x0(2)*n(2) + x0(3)*n(3) ) n = -n
            end if
            dd = 1.0/sqrt( n(1)*n(1) + n(2)*n(2) + n(3)*n(3) )
            n = n * dd
            return
        end function getNorm


!-------

        subroutine generateNeighbourList(this,full)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      find which triangles each node maps to
    !*      if full = true then can't assume each node is distinct- first go through the list culling
    !*      those nodes which are duplicates
            type(TriangulatedSurface),intent(inout)     ::      this
            logical,intent(in),optional                 ::      full
            integer,dimension(0:min(this%n_triangle,100),this%n_node)       ::      neigh
            integer                 ::      ii,jj,kk,mm,nn
            real(kind=real64),dimension(3)       ::      xi,dx
            real(kind=real64)                    ::      dd
            integer                 ::      found
            if (present(full)) then
                if (full) then
                    nn = 0              !   count of distinct nodes
                    do ii = 1,this%n_node
                        xi = this%node(:,ii)
                        found = 0
                        do jj = 1,nn
                            dx = this%node(:,jj) - xi
                            dd = (dx(1)*dx(1) + dx(2)*dx(2) + dx(3)*dx(3))
                            if (dd < 1.0d-8) then
                                !   node ii is in fact the same point as node jj
                                !   do not add ii to the list of distinct nodes
                                found = jj
                                exit
                            end if
                        end do
                        if (found == 0) then
                            !   add node ii to the list of distinct nodes
                            nn = nn + 1
                            this%node(:,nn) = xi
                            found = nn
                        end if
                        !   now make sure all the triangles which thought they include node "ii"
                        !   now know they actually include node "found"
                        do jj = 1,this%n_triangle
                            do kk = 1,3
                                if (this%triangle(kk,jj) == ii) this%triangle(kk,jj)=found
                            end do
                        end do
                    end do
                    print *,"reduced node set from ",this%n_node," to ",nn,maxval(this%triangle(:,1:this%n_triangle))
                    this%n_node = nn                    
                    do jj = 1,this%n_triangle
                        do kk = 1,3
                            if (this%triangle(kk,jj) > nn) print *,"error corner ",kk,jj,this%triangle(kk,jj)
                        end do
                    end do

!                    call expandNeighbourList(this,nn*2)

                end if
            end if




            neigh = 0
            do ii = 1,this%n_triangle               !   for each triangle
                do jj = 1,3
                    kk = this%triangle(jj,ii)       !   find the node on a corner
                    mm = neigh(0,kk) + 1
                    neigh(mm,kk) = ii               !   add triangle ii to the list
                    neigh(0,kk) = mm
                end do
            end do
            mm = maxval(neigh(0,:))
!           print *,"maxval ",mm
            call expandNeighbourList(this,mm)
            this%neigh = 0
            do ii = 1,this%n_triangle
                do jj = 1,3
                    kk = this%triangle(jj,ii)
                    mm = this%neigh(0,kk) + 1
                    this%neigh(mm,kk) = ii
                    this%neigh(0,kk) = mm
                end do
            end do
            return
        end subroutine generateNeighbourList

    !---

        pure function hasNeighbourList(this) result(has)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      returns true if an up-to-date neighbour list is available
            type(TriangulatedSurface),intent(in)        ::      this
            logical                                     ::      has
            has = (this%N_neigh > 0)
            return
        end function hasNeighbourList


    !---

        subroutine removeNeighbourList(this)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      for one reason or other, the neighbour list is either out of date or not needed
    !*      deallocate the memory
            type(TriangulatedSurface),intent(inout)     ::      this
            if (.not. hasNeighbourList(this)) return
            this%N_neigh = 0
            deallocate(this%neigh)
            nullify(this%neigh)
            return
        end subroutine removeNeighbourList

        subroutine expandNeighbourList(this,n)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      ensure there is sufficient memory for n neighbours in list
            type(TriangulatedSurface),intent(inout)     ::      this
            integer,intent(in)                          ::      n
            integer,dimension(:,:),pointer      ::  tnp
            if (this%N_neigh >= n) return
            if (.not. hasNeighbourList(this)) then
                allocate(this%neigh(0:n,this%N_array))
            else
                allocate(tnp(0:n,this%N_array))
                tnp = 0
                tnp(0:this%N_neigh,1:this%n_node) = this%neigh(0:this%N_neigh,1:this%n_node)
                deallocate(this%neigh)
                this%neigh => tnp
            end if
            this%N_neigh = n
            return
        end subroutine expandNeighbourList

!-------

        subroutine intersectLineTriangle( y , a,b , lambda , ok , positive )
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      if the line x = a + lambda b intersects the triangle
    !*      described by the three points y then return the lambda value at which it does
    !*      otherwise return ok = .false.
            real(kind=real64),dimension(3,3),intent(in)         ::      y                       !   positions of vertices of the triangle
            real(kind=real64),dimension(3),intent(in)           ::      a,b                     !   describing the line                                            
            real(kind=real64),intent(out)                       ::      lambda                  !   point of intersection                        
            logical,intent(out)                                 ::      ok                      !   does it intersect?
            logical,intent(in),optional                         ::      positive                !   only accept if the triangle normal is strictly pointing in same direction as line
            
            real(kind=real64),dimension(3)       ::      y31,y21,yp1 !   vector from 1-3 , from 1-2, from 1-p
            real(kind=real64),dimension(3)       ::      nn          !   unnormalised normal to triangle
            real(kind=real64)                    ::      dd,ee
            real(kind=real64)                    ::      v00,v01,v02,v12,v11
            lambda = 0.0
            y21 = y(1:3,2) - y(1:3,1)
            y31 = y(1:3,3) - y(1:3,1)
        !---    equation of the plane is x.n = d , so a.n + lambda b.n = d , lambda = ( x.n - a.n )/b.n
        
        !---    construct unnormalised normal to plane of triangle
            nn = crossProduct( y21,y31 )
            
        !---    check for outward pointing normal 
            yp1(1:3) = y(1:3,1)+y(1:3,2)+y(1:3,3) - 3*a(1:3) !   3x centroid of triangle wrt a
            if (nn(1)*yp1(1) + nn(2)*yp1(2) + nn(3)*yp1(3)  <0) nn = - nn

            ee = b(1)*nn(1) + b(2)*nn(2) + b(3)*nn(3)
            ok = (abs(ee) > 1.0d-12)
            if (present(positive)) then
                if (positive) ok = ee > 1.0d-12
            end if
            if (.not. ok) return
            dd = ( y(1,1)*nn(1) + y(2,1)*nn(2) + y(3,1)*nn(3) )

        !---    find distance of point a to plane
            lambda = ( dd - a(1)*nn(1) - a(2)*nn(2) - a(3)*nn(3) ) / ee

        !---    point on plane is p = a + lambda b = y1 + dd y31 + ee y21
            yp1(1:3) = a + lambda * b - y(:,1)

        !---    find those distances dd and ee.
            v00 = y31(1)*y31(1) + y31(2)*y31(2) + y31(3)*y31(3)
            v11 = y21(1)*y21(1) + y21(2)*y21(2) + y21(3)*y21(3)
            v01 = y31(1)*y21(1) + y31(2)*y21(2) + y31(3)*y21(3)
            v02 = y31(1)*yp1(1) + y31(2)*yp1(2) + y31(3)*yp1(3)
            v12 = y21(1)*yp1(1) + y21(2)*yp1(2) + y21(3)*yp1(3)

            dd = v12*v00 - v02*v01
            ee = v02*v11 - v12*v01

        !---    to be inside triangle need 0 <= dd,ee,dd+ee <= v00*v11 - v01*v01
            ok = ( (dd>-1.0d-12).and.(ee>-1.0d-12) )
            ok = ok .and. ( dd+ee-1.0d-12 < v00*v11 - v01*v01 )
            ok = ok .and. ( dd-1.0d-12 < v00*v11 - v01*v01 ) .and. ( ee-1.0d-12 < v00*v11 - v01*v01 )
            
            return
        end subroutine intersectLineTriangle

!-------

        function pointOnSurface(this,x) result(i)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      is point x on the surface?
    !*      if so return the triangle number i, if not return 0
            type(TriangulatedSurface),intent(in)        ::      this
            real(kind=real64),dimension(3),intent(in)                ::      x
            integer                                     ::      i
            integer                 ::      ii
            real(kind=real64),dimension(3)       ::      norm
            logical                 ::      ok
            real(kind=real64)                    ::      lambda
            i = 0
            do ii = 1,this%n_triangle
                norm = getNorm(this,ii)
                call intersectLineTriangle( getTriangle(this,ii) , x,norm , lambda , ok , positive=.true. )
                if (ok) then
                    if (lambda < 1.0d-12) then
                        i = ii
                        return
                    end if
                end if
            end do
            return
        end function pointOnSurface



!-------


        function within0(this,x) result(is)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      returns true if the point x is inside the surface
    !*      ( note- this is pretty meaningless unless the surface is closed... )
            type(TriangulatedSurface),intent(in)        ::      this
            real(kind=real64),dimension(3),intent(in)                ::      x
            logical                                     ::      is
            real(kind=real64)        ::      dd
            dd = x(1)*x(1) + x(2)*x(2) + x(3)*x(3)
            is = (dd <= this%rmin2)
            if (is) return
            is = (dd <= this%rmax2)
            if (.not. is) return
            is = within1(this,x,raytrace = .true.)
            return
        end function within0


        function within1(this,x,raytrace) result(is)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      returns true if the point x is inside the surface
    !*      ( note- this is pretty meaningless unless the surface is closed... )
            type(TriangulatedSurface),intent(in)        ::      this
            real(kind=real64),dimension(3),intent(in)                ::      x
            logical,intent(in)                          ::      raytrace
            logical                                     ::      is
            integer             ::      ii
            real(kind=real64),dimension(3)   ::      bb
            real(kind=real64)                ::      lambda,lmin
            logical             ::      ok
        !---    to be inside the surface you must
        !       be able to head out along an arbitrary direction b
        !       and hit the back side of a triangle _before_ hitting a front side.
        !       ie the minimum distance along must be +ve.
            bb = (/ 1.0,0.0,0.0 /)
            lmin = huge(1.0)
            do ii = 1,this%n_triangle
                call intersectLineTriangle( getTriangle(this,ii), x,bb , lambda , ok , positive=.true.)
                if (ok) lmin = min(lmin,lambda)
            end do
            is = (lmin < 0.5*huge(1.0))
            return
        end function within1

    !---

        subroutine distanceToSurface0(this,a,b,lambda,ok)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      if you move along the line x = a + lambda b
    !*      compute the minimum +ve lambda to hit the surface
    !*      or return ok = .false. if not successful
            type(TriangulatedSurface),intent(in)        ::      this
            real(kind=real64),dimension(3),intent(in)                ::      a,b
            real(kind=real64),intent(out)                            ::      lambda
            logical,intent(out)                         ::      ok
            real(kind=real64)        ::      dd
            integer     ::      ii
            lambda = 0.0
            dd = a(1)*a(1) + a(2)*a(2) + a(3)*a(3)
            ok = (dd <= this%rmax2)
!D$         if (.not. within(this,a)) &
!D$             write(*,fmt='(a,6f16.8)') "Lib_TriangulatedSurface::distanceToSurface0 warning - point outside cell",a,this%rmax2,dd
            if (.not. ok) return
            dd = huge(1.0)
            do ii = 1,this%n_triangle
                call intersectLineTriangle( getTriangle(this,ii), a,b , lambda , ok , positive=.true.)
                if (ok) then
                    dd = min(dd,lambda)
!                   print *,"intersectLineTriangle ",ii,dd,lambda
                end if
            end do
            ok = (dd < 0.5*huge(1.0))
            lambda = dd
            return
        end subroutine distanceToSurface0


!-------

        pure subroutine findBoundingSpheres(this)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      find the radius of the spheres circumscribing and inscribing the surface
            type(TriangulatedSurface),intent(inout)     ::      this
            real(kind=real64)        ::      dd
            integer     ::      ii

            this%rmin2 = huge(1.0)
            this%rmax2 = 0.0
            do ii = 1,this%n_node
                dd = this%node(1,ii)*this%node(1,ii) + this%node(2,ii)*this%node(2,ii)  &
                   + this%node(3,ii)*this%node(3,ii)
                this%rmin2 = min(this%rmin2,dd)
                this%rmax2 = max(this%rmax2,dd)
            end do
            return
        end subroutine findBoundingSpheres

        subroutine scaleRadius(this,x)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      scale the radius of the surface by a constant factor x
    !*      useful if we have started with a unit sphere from GeodesicGrid
            type(TriangulatedSurface),intent(inout)     ::      this
            real(kind=real64),intent(in)                             ::      x
            this%node = this%node * x
            return
        end subroutine scaleRadius

!-------

        pure subroutine findNormalAndAreaAndCentroid( tri,n,A,c, x0 )
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      find a normal vector to the triangle and area and centre of mass simultaneously
    !*      If x0 is present, normal points away from it
            real(kind=real64),dimension(3,3),intent(in)          ::      tri
            real(kind=real64),dimension(3),intent(out)           ::      n,c
            real(kind=real64),intent(out)                        ::      A
            real(kind=real64),dimension(3),intent(in),optional   ::      x0
            real(kind=real64)            ::      nn

            n = crossProduct( tri(1:3,2)-tri(1:3,1),tri(1:3,3)-tri(1:3,1) )
            nn   = n(1)*n(1) + n(2)*n(2) + n(3)*n(3)
            if (nn > 1.0d-16) then
                nn = sqrt(nn)
                n = n / nn
                if (present(x0)) then
                    if ( n(1)*( tri(1,1)+tri(1,2)+tri(1,3)-3*x0(1) )                    &
                       + n(2)*( tri(2,1)+tri(2,2)+tri(2,3)-3*x0(2) )                    &
                       + n(3)*( tri(3,1)+tri(3,2)+tri(3,3)-3*x0(3) ) < 0 ) n = -n
                end if
                A = 0.5*nn
            else
                n = 0.0
                A = 0.0
            end if
            c(:) = ( tri(1:3,1) + tri(1:3,2) + tri(1:3,3) ) / 3.0
            return
        end subroutine findNormalAndAreaAndCentroid

        pure subroutine findNormalAndSignedAreaAndCentroid( tri,n,A,c, x0 )
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      find a normal vector to the triangle and area and centre of mass simultaneously
    !*      If x0 is present, normal points away from it
            real(kind=real64),dimension(3,3),intent(in)          ::      tri
            real(kind=real64),dimension(3),intent(out)           ::      n,c
            real(kind=real64),intent(out)                        ::      A
            real(kind=real64),dimension(3),intent(in)   ::      x0
            real(kind=real64)            ::      nn

            n = crossProduct( tri(1:3,2)-tri(1:3,1),tri(1:3,3)-tri(1:3,1) )
            nn   = n(1)*n(1) + n(2)*n(2) + n(3)*n(3)
            if (nn > 1.0d-16) then
                nn = sqrt(nn)
                n = n / nn
                if ( n(1)*( tri(1,1)+tri(1,2)+tri(1,3)-3*x0(1) )                    &
                    + n(2)*( tri(2,1)+tri(2,2)+tri(2,3)-3*x0(2) )                    &
                    + n(3)*( tri(3,1)+tri(3,2)+tri(3,3)-3*x0(3) ) < 0 ) then
                    A = -0.5*nn
                else
                    A = 0.5*nn
                end if
            else
                n = 0.0
                A = 0.0
            end if
            c(:) = ( tri(1:3,1) + tri(1:3,2) + tri(1:3,3) ) / 3.0
            return
        end subroutine findNormalAndSignedAreaAndCentroid

        pure function findVolume0(this,x0) result(v)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      find the volume of a triangulated surface enclosing x0 using
    !*      the divergence theorem.
    !*          F(x,y,z)  = (x+y+z)/3
    !*          div F = 1
    !*      int div F dV = V = int F.n dS
            type(TriangulatedSurface),intent(in)        ::      this
            real(kind=real64),dimension(3),intent(in),optional    ::      x0
            real(kind=real64)                            ::      v
            integer                 ::      ii
            real(kind=real64),dimension(3,3)     ::      tri
            real(kind=real64),dimension(3)       ::      nn

            v = 0.0
            if (present(x0)) then
                do ii = 1,this%n_triangle
                    tri = getTriangle(this,ii)
                    nn = crossProduct( tri(1:3,2)-x0,tri(1:3,3)-x0 )
                    v = v + (tri(1,1)-x0(1))*nn(1) + (tri(2,1)-x0(2))*nn(2) + (tri(3,1)-x0(3))*nn(3)
                end do
            else
                do ii = 1,this%n_triangle
                    tri = getTriangle(this,ii)
                    nn = crossProduct( tri(1:3,2),tri(1:3,3) )
                    v = v + tri(1,1)*nn(1) + tri(2,1)*nn(2) + tri(3,1)*nn(3)
                end do
            end if
            v = v / 6.0

            return
        end function findVolume0

        pure function findArea0(this,x0,onSurface) result(a)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      find the area of a triangulated surface enclosing x0
            type(TriangulatedSurface),intent(in)        ::      this
            real(kind=real64),dimension(3),intent(in)                ::      x0
            logical,dimension(:),intent(in),optional    ::  onSurface
            real(kind=real64)                            ::      a
            integer                 ::      ii
            real(kind=real64),dimension(3,3)     ::      tri
            real(kind=real64),dimension(3)       ::      nn,cc
            integer,dimension(3)    ::      trin
            real(kind=real64)                    ::      aa
            logical,dimension(3)    ::      on

            a = 0.0
            if (present(onSurface)) then
                do ii = 1,this%n_triangle
                    call getTriangleInfo(this,ii,tri,trin)
                    on = (/ onSurface(trin(1)), onSurface(trin(2)), onSurface(trin(3)) /)
                    call findNormalAndSignedAreaAndCentroid( tri,nn,aa,cc, x0 )
                    a = a + aa*count(on)/3.0
                end do
            else
                do ii = 1,this%n_triangle
                    tri = getTriangle(this,ii)
                    call findNormalAndSignedAreaAndCentroid( tri,nn,aa,cc, x0 )
                    a = a + aa
                end do
            end if
            return
        end function findArea0

        pure function findArea1(this,x0,onSurfaceNode,onSurfaceTri) result(a)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      find the area of a triangulated surface enclosing x0
            type(TriangulatedSurface),intent(in)        ::      this
            real(kind=real64),dimension(3),intent(in)                ::      x0
            logical,dimension(:),intent(in)             ::  onSurfaceNode
            logical,dimension(:),intent(in)             ::  onSurfaceTri
            real(kind=real64)                            ::      a
            integer                 ::      ii
            real(kind=real64),dimension(3,3)     ::      tri
            real(kind=real64),dimension(3)       ::      nn,cc
            integer,dimension(3)    ::      trin
            real(kind=real64)                    ::      aa
            logical,dimension(3)    ::      on

            a = 0.0
            do ii = 1,this%n_triangle
                call getTriangleInfo(this,ii,tri,trin)
                on = (/ onSurfaceNode(trin(1)), onSurfaceNode(trin(2)), onSurfaceNode(trin(3)) /)
                if (all(on).and. .not. onSurfaceTri(ii)) cycle
                call findNormalAndSignedAreaAndCentroid( tri,nn,aa,cc, x0 )
                a = a + aa*count(on)/3.0
            end do
            return
        end function findArea1

!-------


        function findNormals0( this,signed ) result( norm )
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
            type(TriangulatedSurface),intent(in)        ::      this
            logical,intent(in),optional                 ::      signed
            real(kind=real64),dimension(3,this%n_node)               ::      norm
            integer             ::      ii
            real(kind=real64),dimension(3)   ::      nn,cc
            real(kind=real64)                ::      dd
            if (present(signed)) then
                if (signed) then
                    norm = 0.0
                    do ii = 1,this%n_triangle
                        cc = this%node( 1:3,this%triangle(2,ii) )-this%node( 1:3,this%triangle(1,ii) )
                        nn = this%node( 1:3,this%triangle(3,ii) )-this%node( 1:3,this%triangle(1,ii) )
                        nn = crossProduct( cc,nn )
                        dd = nn(1)*nn(1) + nn(2)*nn(2) + nn(3)*nn(3)
                        dd = 1.0
                        if (dd>1.0d-16) then
                            norm(:,this%triangle(1,ii)) = norm(:,this%triangle(1,ii)) + nn(:)*dd
                            norm(:,this%triangle(2,ii)) = norm(:,this%triangle(2,ii)) + nn(:)*dd
                            norm(:,this%triangle(3,ii)) = norm(:,this%triangle(3,ii)) + nn(:)*dd
                        end if
                    end do
                    do ii = 1,this%n_node
                        nn = norm(:,ii)
                        dd = nn(1)*nn(1) + nn(2)*nn(2) + nn(3)*nn(3)
                        if (dd > 1.0d-16) then
                            dd = 1.0/sqrt(dd)
                            norm(:,ii) = norm(:,ii)*dd
                        else
                            cc = this%node( :,ii )
                            dd = cc(1)*cc(1) + cc(2)*cc(2) + cc(3)*cc(3)
                            if (dd > 1.0d-16) then
                                dd = 1.0/sqrt(dd)
                                norm(:,ii) = cc * dd
                            else
                                norm(:,ii) = (/ 0.0,0.0,1.0 /)
                            end if
                        end if
                    end do

                    return
                end if
            end if
            norm = 0.0
            if (maxval(this%triangle(:,:)) > this%n_node) & 
            print *,"findNormals0 this%n_triangle,size(this%node,dim=2),size(this%triangle,dim=2),maxval(this%triangle(:,:)),this%n_node "      &
                                 ,this%n_triangle,size(this%node,dim=2),size(this%triangle,dim=2),maxval(this%triangle(:,:)),this%n_node
            do ii = 1,this%n_triangle
                cc = this%node( 1:3,this%triangle(2,ii) )-this%node( 1:3,this%triangle(1,ii) )
                nn = this%node( 1:3,this%triangle(3,ii) )-this%node( 1:3,this%triangle(1,ii) )
                nn = crossProduct( cc,nn )
                cc = ( this%node( :,this%triangle(1,ii) )       &
                     + this%node( :,this%triangle(2,ii) )           &
                     + this%node( :,this%triangle(3,ii) ) )*0.3333333333333333
                if ( cc(1)*nn(1) + cc(2)*nn(2) + cc(3)*nn(3) < 0.0 ) nn = -nn
                dd = nn(1)*nn(1) + nn(2)*nn(2) + nn(3)*nn(3)
                dd = 1.0
                if (dd>1.0d-16) then
                    norm(:,this%triangle(1,ii)) = norm(:,this%triangle(1,ii)) + nn(:)*dd
                    norm(:,this%triangle(2,ii)) = norm(:,this%triangle(2,ii)) + nn(:)*dd
                    norm(:,this%triangle(3,ii)) = norm(:,this%triangle(3,ii)) + nn(:)*dd
                end if
            end do
            do ii = 1,this%n_node
                nn = norm(:,ii)
                dd = nn(1)*nn(1) + nn(2)*nn(2) + nn(3)*nn(3)
                if (dd > 1.0d-16) then
                    dd = 1.0/sqrt(dd)
                    norm(:,ii) = norm(:,ii) * dd
                else
                    cc = this%node( :,ii )
                    dd = cc(1)*cc(1) + cc(2)*cc(2) + cc(3)*cc(3)
                    if (dd > 1.0d-16) then
                        dd = 1.0/sqrt(dd)
                        norm(:,ii) = cc * dd
                    else
                        norm(:,ii) = (/ 0.0,0.0,1.0 /)
                    end if
                end if
            end do
            return
        end function findNormals0


        function findNormals1( this,i ) result( norm )
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      find normal of node i only
            type(TriangulatedSurface),intent(in)        ::      this
            integer,intent(in)                          ::      i
            real(kind=real64),dimension(3)                           ::      norm
            integer             ::      jj,kk
            real(kind=real64),dimension(3)   ::      nn,cc
            real(kind=real64)                ::      dd
            norm = 0.0

            do jj = 1,this%neigh(0,i)
                kk = this%neigh(jj,i)
                cc = this%node( 1:3,this%triangle(2,kk) )-this%node( 1:3,this%triangle(1,kk) )
                nn = this%node( 1:3,this%triangle(3,kk) )-this%node( 1:3,this%triangle(1,kk) )
                nn = crossProduct( cc,nn )
                cc = ( this%node( :,this%triangle(1,kk) )       &
                     + this%node( :,this%triangle(2,kk) )           &
                     + this%node( :,this%triangle(3,kk) ) )*0.3333333333333333
                if ( cc(1)*nn(1) + cc(2)*nn(2) + cc(3)*nn(3) < 0.0 ) nn = -nn
                dd = nn(1)*nn(1) + nn(2)*nn(2) + nn(3)*nn(3)
                dd = 1.0
                if (dd>1.0d-16) then
                    norm(:) = norm(:) + nn(:)*dd
                    norm(:) = norm(:) + nn(:)*dd
                    norm(:) = norm(:) + nn(:)*dd
                end if
            end do

            dd = norm(1)*norm(1) + norm(2)*norm(2) + norm(3)*norm(3)
            if (dd>1.0d-12) then
                norm = norm / sqrt(dd)
            else
                norm = (/ 0.0,0.0,1.0 /)
            end if
            return
        end function findNormals1


!-------

        subroutine findVolumeAndArea0(this,x0,v,a,onSurface)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      find the volume of a triangulated surface enclosing x0 using
    !*      the divergence theorem.
    !*          F(x,y,z)  = (x+y+z)/3
    !*          div F = 1
    !*      int div F dV = V = int F.n dS
            type(TriangulatedSurface),intent(in)        ::      this
            real(kind=real64),dimension(3),intent(in)    ::      x0
            real(kind=real64),intent(out)                ::      V,A
            logical,dimension(:),intent(in),optional    ::      onSurface
            integer                 ::      ii
            real(kind=real64),dimension(3,3)     ::      tri,trinorm
            real(kind=real64),dimension(3)       ::      nn,cc
            real(kind=real64)                    ::      aa,dV,dA,vv
            real(kind=real64),dimension(3,this%n_node)   ::  norm
            logical                 ::      ok
            logical,dimension(3)    ::      on
            V = 0.0
            A = 0.0
            norm = findNormals( this )
            if (present(onSurface)) then
                do ii = 1,this%n_triangle
                    on = (/ onSurface(this%triangle(1,ii)), onSurface(this%triangle(2,ii)), onSurface(this%triangle(3,ii)) /)
                    tri(:,1) = this%node( :,this%triangle(1,ii) )
                    tri(:,2) = this%node( :,this%triangle(2,ii) )
                    tri(:,3) = this%node( :,this%triangle(3,ii) )
                    trinorm(:,1) = norm( :,this%triangle(1,ii) )
                    trinorm(:,2) = norm( :,this%triangle(2,ii) )
                    trinorm(:,3) = norm( :,this%triangle(3,ii) )

                    call findNormalAndSignedAreaAndCentroid( tri,nn,aa,cc, x0 )
                    vv = abs(aa)*( nn(1)*cc(1) + nn(2)*cc(2) + nn(3)*cc(3) )/3.0
                    call sphereToTriangle( tri,trinorm, dV,dA, ok )
                    if (ok) then
                        if (abs(dV) < vv) V = V + dV!*count(on)*0.333333333333333

!                           A = A + dA*count(on)*count(on)*0.1111111111111111
!                           V = V + dV*count(on)*0.333333333333333
!                       end if
                        if (abs(dA) < aa) A = A + dA*count(on)*0.333333333333333
!                       if ( dV<0) write (*,fmt='(i6,3l2,100f20.12)') ii,on,vv,dV,aa,dA,cc  ,   &
!                           dot_product( abs(cc),(/ 1.0,1.0,1.0 /) )
                    end if
                    A = A + aa*count(on)*0.333333333333333
                    V = V + vv
                end do
            else
                do ii = 1,this%n_triangle
                    tri(:,1) = this%node( :,this%triangle(1,ii) )
                    tri(:,2) = this%node( :,this%triangle(2,ii) )
                    tri(:,3) = this%node( :,this%triangle(3,ii) )
                    trinorm(:,1) = norm( :,this%triangle(1,ii) )
                    trinorm(:,2) = norm( :,this%triangle(2,ii) )
                    trinorm(:,3) = norm( :,this%triangle(3,ii) )
                    call findNormalAndSignedAreaAndCentroid( tri,nn,aa,cc, x0 )                                                                
                    vv = abs(aa)*( nn(1)*cc(1) + nn(2)*cc(2) + nn(3)*cc(3) )/3.0
                    call sphereToTriangle( tri,trinorm, dV,dA, ok )
                    if (ok) then
                        A = A + dA
                        V = V + dV
                    end if
                    A = A + aa
                    V = V + vv
                end do
            end if
            return
        end subroutine findVolumeAndArea0

        function findArea2( tri ) result(aa)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
            real(kind=real64),dimension(3,3),intent(in)      ::      tri
            real(kind=real64)                                ::      aa
            real(kind=real64),dimension(3)       ::      nn

            nn  = crossProduct( tri(1:3,2)-tri(1:3,1),tri(1:3,3)-tri(1:3,1) )
            aa  = sqrt( nn(1)*nn(1) + nn(2)*nn(2) + nn(3)*nn(3) )*0.5

            return
        end function findArea2


        subroutine findAreaCorrection( tri,trinorm,onSurface,x0, aa,dA )
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
            real(kind=real64),dimension(3,3),intent(in)      ::      tri,trinorm
            logical,dimension(3),intent(in)     ::      onSurface
            real(kind=real64),dimension(3),intent(in)        ::      x0
            real(kind=real64),intent(out)                    ::      aa,dA
            real(kind=real64),dimension(3)       ::      nn,cc
            real(kind=real64)                    ::      vv,dV
            logical                 ::      ok

            call findNormalAndSignedAreaAndCentroid( tri,nn,aa,cc, x0 )
            vv = abs(aa)*( nn(1)*cc(1) + nn(2)*cc(2) + nn(3)*cc(3) )/3.0
            call sphereToTriangle( tri,trinorm, dV,dA, ok )
            if (ok) then
                if ( (abs(dV) < vv).and.(abs(dA) < aa) ) then
                    dA = dA*count(onSurface)*count(onSurface)*0.1111111111111111
                end if
            end if
            aa = aa*count(onSurface)*0.333333333333333
            return
        end subroutine findAreaCorrection

        subroutine findVolumeAndArea1(this,norm,x0,v,a,onSurface)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      find the volume of a triangulated surface enclosing x0 using
    !*      the divergence theorem.
    !*          F(x,y,z)  = (x+y+z)/3
    !*          div F = 1
    !*      int div F dV = V = int F.n dS
            type(TriangulatedSurface),intent(in)        ::      this
            real(kind=real64),dimension(3,this%n_node),intent(in)    ::      norm
            real(kind=real64),dimension(3),intent(in)    ::      x0
            real(kind=real64),intent(out)                ::      V,A
            logical,dimension(:),intent(in),optional    ::      onSurface
            integer                 ::      ii
            real(kind=real64),dimension(3,3)     ::      tri,trinorm
            real(kind=real64),dimension(3)       ::      nn,cc
            real(kind=real64)                    ::      aa,dV,dA,vv , xx
            logical                 ::      ok
            logical,dimension(3)    ::      on
            V = 0.0
            A = 0.0
!             norm = findNormals( this )
            if (present(onSurface)) then
                do ii = 1,this%n_triangle
                    on = (/ onSurface(this%triangle(1,ii)), onSurface(this%triangle(2,ii)), onSurface(this%triangle(3,ii)) /)
                    tri(:,1) = this%node( :,this%triangle(1,ii) )
                    tri(:,2) = this%node( :,this%triangle(2,ii) )
                    tri(:,3) = this%node( :,this%triangle(3,ii) )
                    trinorm(:,1) = norm( :,this%triangle(1,ii) )
                    trinorm(:,2) = norm( :,this%triangle(2,ii) )
                    trinorm(:,3) = norm( :,this%triangle(3,ii) )

                    call findNormalAndSignedAreaAndCentroid( tri,nn,aa,cc, x0 )
                    xx = 0.0
!                     nnnn = crossProduct( tri(:,2)-tri(:,1), tri(:,3)-tri(:,1) )
!                     xx = cc(1)*( trinorm(1,1)+trinorm(1,2)+trinorm(1,3) ) + cc(2)*( trinorm(2,1)+trinorm(2,2)+trinorm(2,3) ) + cc(3)*( trinorm(3,1)+trinorm(3,2)+trinorm(3,3) )
!                     if (xx<0) then
!                       aa = - aa; nn = - nn
!                     end if
                    aa = abs(aa)
                    vv = aa*( nn(1)*cc(1) + nn(2)*cc(2) + nn(3)*cc(3) )/3.0
                    call sphereToTriangle( tri,trinorm, dV,dA, ok )
!                    print *,"findVolumeAndArea1 ",ii,xx,aa,vv,dA*count(on)*count(on)*0.1111111111111111,dV*count(on)*0.333333333333333,ok
!                     dA = 0.0 ; dV = 0.0
                    if (ok) then
                        if ( (abs(dV) < abs(vv)).and.(abs(dA) < abs(aa)) ) then
                            A = A + dA*count(on)*count(on)*0.1111111111111111
                            V = V + dV*count(on)*0.333333333333333
                        end if
!                       if ( dV<0) write (*,fmt='(i6,3l2,100f20.12)') ii,on,vv,dV,aa,dA,cc  ,   &
!                           dot_product( abs(cc),(/ 1.0,1.0,1.0 /) )
                    end if
                    A = A + aa*count(on)*0.333333333333333
                    V = V + vv
                end do
            else
                do ii = 1,this%n_triangle
                    tri(:,1) = this%node( :,this%triangle(1,ii) )
                    tri(:,2) = this%node( :,this%triangle(2,ii) )
                    tri(:,3) = this%node( :,this%triangle(3,ii) )
                    trinorm(:,1) = norm( :,this%triangle(1,ii) )
                    trinorm(:,2) = norm( :,this%triangle(2,ii) )
                    trinorm(:,3) = norm( :,this%triangle(3,ii) )
                    call findNormalAndSignedAreaAndCentroid( tri,nn,aa,cc, x0 )
                    aa = abs(aa)
                    vv = aa*( nn(1)*cc(1) + nn(2)*cc(2) + nn(3)*cc(3) )/3.0
                    call sphereToTriangle( tri,trinorm, dV,dA, ok )
                    if (ok) then
                        A = A + dA
                        V = V + dV
                    end if
                    A = A + aa
                    V = V + vv
                end do
            end if
            V = abs(V) ; A = abs(A)
            
            
            return
        end subroutine findVolumeAndArea1

        subroutine findVolumeAndArea2(this,norm,x0,v,a,onSurfaceNode,onSurfaceTri)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      find the volume of a triangulated surface enclosing x0 using
    !*      the divergence theorem.
    !*          F(x,y,z)  = (x+y+z)/3
    !*          div F = 1
    !*      int div F dV = V = int F.n dS
            type(TriangulatedSurface),intent(in)        ::      this
            real(kind=real64),dimension(3,this%n_node),intent(in)    ::      norm
            real(kind=real64),dimension(3),intent(in)    ::      x0
            real(kind=real64),intent(out)                ::      V,A
            logical,dimension(:),intent(in) ::      onSurfaceNode
            logical,dimension(:),intent(in) ::      onSurfaceTri
            integer                 ::      ii
            real(kind=real64),dimension(3,3)     ::      tri,trinorm
            real(kind=real64),dimension(3)       ::      nn,cc
            real(kind=real64)                    ::      aa,dV,dA,vv
            logical                 ::      ok
            logical,dimension(3)    ::      on
            V = 0.0
            A = 0.0
            do ii = 1,this%n_triangle
                on = (/ onSurfaceNode(this%triangle(1,ii)),     &
                        onSurfaceNode(this%triangle(2,ii)),     &
                        onSurfaceNode(this%triangle(3,ii)) /)
                tri(:,1) = this%node( :,this%triangle(1,ii) )
                tri(:,2) = this%node( :,this%triangle(2,ii) )
                tri(:,3) = this%node( :,this%triangle(3,ii) )
                trinorm(:,1) = norm( :,this%triangle(1,ii) )
                trinorm(:,2) = norm( :,this%triangle(2,ii) )
                trinorm(:,3) = norm( :,this%triangle(3,ii) )

                call findNormalAndSignedAreaAndCentroid( tri,nn,aa,cc, x0 )
                vv = abs(aa)*( nn(1)*cc(1) + nn(2)*cc(2) + nn(3)*cc(3) )/3.0
                call sphereToTriangle( tri,trinorm, dV,dA, ok )
                if (ok) then
                    if ( (abs(dV) < vv).and.(abs(dA) < aa) ) then
                        if (onSurfaceTri(ii)) A = A + dA*count(on)*count(on)*0.1111111111111111
                        V = V + dV*count(on)*0.333333333333333
                    end if
                end if
                if (onSurfaceTri(ii)) A = A + aa*count(on)*0.333333333333333
                V = V + vv
            end do
            V = abs(V) ; A = abs(A)
            return
        end subroutine findVolumeAndArea2



        subroutine findVolumeAndArea3(this,x0,v,a,onSurfaceNorm,onSurfaceTri)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      find the volume of a triangulated surface enclosing x0 using
    !*      the divergence theorem.
    !*          F(x,y,z)  = (x+y+z)/3
    !*          div F = 1
    !*      int div F dV = V = int F.n dS
            type(TriangulatedSurface),intent(in)        ::      this
            real(kind=real64),dimension(3),intent(in)    ::      x0
            real(kind=real64),intent(out)                ::      V,A
            logical,dimension(:),intent(in) ::      onSurfaceNorm
            logical,dimension(:),intent(in) ::      onSurfaceTri
            integer                 ::      ii
            real(kind=real64),dimension(3,3)     ::      tri,trinorm
            real(kind=real64),dimension(3)       ::      nn,cc
            real(kind=real64)                    ::      aa,dV,dA,vv
            real(kind=real64),dimension(3,this%n_node)   ::  norm
            logical                 ::      ok
            logical,dimension(3)    ::      on
            V = 0.0
            A = 0.0
            norm = findNormals( this )
            do ii = 1,this%n_triangle
                on = (/ onSurfaceNorm(this%triangle(1,ii)),     &
                        onSurfaceNorm(this%triangle(2,ii)),     &
                        onSurfaceNorm(this%triangle(3,ii)) /)
                tri(:,1) = this%node( :,this%triangle(1,ii) )
                tri(:,2) = this%node( :,this%triangle(2,ii) )
                tri(:,3) = this%node( :,this%triangle(3,ii) )
                trinorm(:,1) = norm( :,this%triangle(1,ii) )
                trinorm(:,2) = norm( :,this%triangle(2,ii) )
                trinorm(:,3) = norm( :,this%triangle(3,ii) )

                call findNormalAndSignedAreaAndCentroid( tri,nn,aa,cc, x0 )
                vv = abs(aa)*( nn(1)*cc(1) + nn(2)*cc(2) + nn(3)*cc(3) )/3.0
                call sphereToTriangle( tri,trinorm, dV,dA, ok )
                if (ok) then
                    if (abs(dV) < vv) V = V + dV
                    if ( onSurfaceTri(ii).and.(abs(dA) < aa) ) A = A + dA*count(on)*0.333333333333333
                end if
                if ( onSurfaceTri(ii) ) A = A + aa*count(on)*0.333333333333333
                V = V + vv
            end do
            V = abs(V) ; A = abs(A)
            return
        end subroutine findVolumeAndArea3



!-------

        subroutine sphereToTriangle( a,n, dV,dA, ok )
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      construct the best fit sphere (x-c)^2 = r^2
    !*      to the points a1,a2,a3 which are associated with the normals
    !*      n1,n2,n3
    !*      find omega, the solid angle subtended. Returns ok = .false. if
    !*      the normals are parallel etc the sphere cannot be found.
    !*      returns the volume and area corrections

            real(kind=real64),dimension(3,3),intent(in)          ::      a,n
            real(kind=real64),intent(out)                        ::      dV,dA
            logical,intent(out)                     ::      ok
            real(kind=real64),dimension(3)       ::      cc,ca1,ca2,ca3
            real(kind=real64),dimension(3)       ::      ca2xca3
            real(kind=real64)                    ::      cc1,cc2,cc3,d1,d2,d3,V,rr,omega

            call midpoint3skew( a(:,1),a(:,2),a(:,3),n(:,1),n(:,2),n(:,3), cc , ok )    !   find where the normals "intersect"
            if (.not. ok) then
                dV = 0.0
                dA = 0.0
                return
            end if

        !---
            ca1 = a(:,1) - cc(:)
            ca2 = a(:,2) - cc(:)
            ca3 = a(:,3) - cc(:)
            ca2xca3 = crossProduct( ca2,ca3 )

            cc1 = ca1(1)*ca1(1) + ca1(2)*ca1(2) + ca1(3)*ca1(3)
            cc2 = ca2(1)*ca2(1) + ca2(2)*ca2(2) + ca2(3)*ca2(3)
            cc3 = ca3(1)*ca3(1) + ca3(2)*ca3(2) + ca3(3)*ca3(3)

        !---    find radius of sphere centre cc which (roughly) passes through a1,a2,a3
            rr = cc1 + cc2 + cc3
            rr = sqrt(rr/3.0)                   !   rms radius

!***********        WHAT IF RADII ARE OPPOSITE SIGNS?

        !---    is sphere above or below triangle?
            ca1 = ( a(:,1) + a(:,2) + a(:,3) )*0.3333333333333333 - cc
            if (ca2xca3(1)*ca1(1) + ca2xca3(2)*ca1(2) + ca2xca3(3)*ca1(3) < 0) rr = -rr

        !---    now have a best-guess sphere. The trouble is that it doesn't actually
        !   pass through the triangle points.
        !   so construct a new triangle which does map exactly to the sphere
            ca1 = n(:,1)
            ca2 = n(:,2)
            ca3 = n(:,3)


            ca2xca3 = crossProduct( ca2,ca3 )


        !---    for solid angle see Oosterom and Strackee ( via wikipedia's page on solid angles )
            d3 = ca1(1)*ca2(1) + ca1(2)*ca2(2) + ca1(3)*ca2(3)
            d1 = ca3(1)*ca2(1) + ca3(2)*ca2(2) + ca3(3)*ca2(3)
            d2 = ca3(1)*ca1(1) + ca3(2)*ca1(2) + ca3(3)*ca1(3)

            V = abs( (ca1(1)*ca2xca3(1) + ca1(2)*ca2xca3(2) + ca1(3)*ca2xca3(3) ) ) ! volume of tetrahedron * 6/(r^3)
            omega = 1 + d1 + d2 + d3

            if (abs(omega)<1.0d-14) then
                if (V < 1.0d-14) then
                    omega = 0.0
                else
                    omega = 3.1415926535897932384626433832795
                end if
            else
                if (omega < 0) then
                    omega = atan( abs(V)/omega ) + 3.1415926535897932384626433832795
                else
                    omega = atan( abs(V)/omega )
                end if
            end if
            omega = 2*omega


        !---    with the radius and solid angle I can compute the volume and area corrections easily
        !       volume correction ( unsigned )
            dV = ( omega - V*0.5) * abs(rr*rr*rr)/3.0

        !       find normal to triangle
            ca2xca3 = crossProduct( ca2-ca1,ca3-ca1 )
            ca1 = cc + rr*( ca1(:) + ca2(:) + ca3(:) )*0.3333333333333333
            d3 = ( ca1(1)*ca2xca3(1) + ca1(2)*ca2xca3(2) + ca1(3)*ca2xca3(3) )
            if (d3 < 0) ca2xca3 = - ca2xca3

        !       which side of triangle does centre of sphere lie?
            ca1 = ca1 - cc
            d2 = ( ca1(1)*ca2xca3(1) + ca1(2)*ca2xca3(2) + ca1(3)*ca2xca3(3) )
            if (d2 < 0) then
                dV = -dV
            end if

            d1 = ( ca2xca3(1)*ca2xca3(1) + ca2xca3(2)*ca2xca3(2) + ca2xca3(3)*ca2xca3(3) )
            d1 = 0.5*sqrt(d1)
    !       area is always increased by taking sphere approx.
            dA = max(0.0,omega - d1)*rr*rr
!
!               write (*,fmt='(100(a,3f16.8))') "S2T a1 ",a(:,1)," a2 ",a(:,2)," a3 ",a(:,3),               &
!                                                  " aa1 ",aa(:,1)," aa2 ",aa(:,2)," aa3 ",aa(:,3),         &
!                                                  " n1 ",n(:,1)," n2 ",n(:,2)," n3 ",n(:,3),               &
!                                                  " c ",cc," r,da,dv ",rr,dA,dV," omega,A,V ",omega,d1*rr*rr,V*abs(rr*rr*rr)/6.0,  &
!                                                  " omega calc ",V,1+d1+d2+d3,V/(1+d1+d2+d3)
!
!           if (abs(rr)>10*sqrt(d1)*rr*rr) then
! ! !               print *,"too flat"
!               dA = 0.0
!               dV = 0.0
!           end if
            return
        end subroutine sphereToTriangle


        subroutine midpoint3skew( a1,a2,a3, b1,b2,b3, c , ok )
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      find the point c which is closest to the three lines
    !*          a1 + lambda1 b1 , a2 + lambda2 b2 , a3 + lambda3 b3

            real(kind=real64),dimension(3),intent(in)            ::      a1,a2,a3
            real(kind=real64),dimension(3),intent(in)            ::      b1,b2,b3
            real(kind=real64),dimension(3),intent(out)           ::      c
            logical,intent(out)                     ::      ok
            real(kind=real64),dimension(3,6)     ::      hex
            logical,dimension(3)    ::      skew
            real(kind=real64),dimension(3)       ::      d
            real(kind=real64)        ::      b1b2

            ok = .true.
            hex = 0.0
            b1b2 = b1(1)*b2(1) + b1(2)*b2(2) + b1(3)*b2(3)
            skew(1) = (1 - b1b2*b1b2 > 0.00001)
            if (skew(1)) then
                c = b1(:) - b2(:)*b1b2
                d = a2(:) - a1(:)
                hex(:,1) = a1(:) + b1(:)*( c(1)*d(1) + c(2)*d(2) + c(3)*d(3) )/( 1 - b1b2*b1b2 )
                c = b1(:)*b1b2 - b2(:)
                hex(:,2) = a2(:) + b2(:)*( c(1)*d(1) + c(2)*d(2) + c(3)*d(3) )/( 1 - b1b2*b1b2 )
            end if

            b1b2 = b2(1)*b3(1) + b2(2)*b3(2) + b2(3)*b3(3)
            skew(2) = (1 - b1b2*b1b2 > 0.00001)
            if (skew(2)) then
                c = b2(:) - b3(:)*b1b2
                d = a3(:) - a2(:)
                hex(:,3) = a2(:) + b2(:)*( c(1)*d(1) + c(2)*d(2) + c(3)*d(3) )/( 1 - b1b2*b1b2 )
                c = b2(:)*b1b2 - b3(:)
                hex(:,4) = a3(:) + b3(:)*( c(1)*d(1) + c(2)*d(2) + c(3)*d(3) )/( 1 - b1b2*b1b2 )
            end if

            b1b2 = b3(1)*b1(1) + b3(2)*b1(2) + b3(3)*b1(3)
            skew(3) = (1 - b1b2*b1b2 > 0.00001)
            if (skew(3)) then
                c = b3(:) - b1(:)*b1b2
                d = a1(:) - a3(:)
                hex(:,5) = a3(:) + b3(:)*( c(1)*d(1) + c(2)*d(2) + c(3)*d(3) )/( 1 - b1b2*b1b2 )
                c = b3(:)*b1b2 - b1(:)
                hex(:,6) = a1(:) + b1(:)*( c(1)*d(1) + c(2)*d(2) + c(3)*d(3) )/( 1 - b1b2*b1b2 )
            end if
!***********        WHAT IF RADII ARE OPPOSITE SIGNS?

!           write (*,fmt='(3l6,7(a,3f12.5))') skew," h1 ",hex(:,1)," h2 ",hex(:,2)," h3 ",hex(:,3)  &
!                                                 ," h4 ",hex(:,4)," h5 ",hex(:,5)," h6 ",hex(:,6)  &
!                                                 ," dots ",b1(1)*b2(1) + b1(2)*b2(2) + b1(3)*b2(3) &
!                                                   ,b2(1)*b3(1) + b2(2)*b3(2) + b2(3)*b3(3)        &
!                                                   ,b3(1)*b1(1) + b3(2)*b1(2) + b3(3)*b1(3)

            if (skew(1)) then
                if (skew(2)) then
                    if (skew(3)) then
                        c = ( hex(:,1) + hex(:,2) + hex(:,3) + hex(:,4) + hex(:,5) + hex(:,6) )*0.16666666666667
                    else
                        c = ( hex(:,1) + hex(:,2) + hex(:,3) + hex(:,4) )*0.25
                    end if
                else
                    if (skew(3)) then
                        c = ( hex(:,1) + hex(:,2) + hex(:,5) + hex(:,6) )*0.25
                    else
                        c = ( hex(:,1) + hex(:,2) )*0.5
                    end if
                end if
            else
                if (skew(2)) then
                    if (skew(3)) then
                        c = ( hex(:,3) + hex(:,4) + hex(:,5) + hex(:,6) )*0.25
                    else
                        c = ( hex(:,3) + hex(:,4) )*0.5
                    end if
                else
                    if (skew(3)) then
                        c = ( hex(:,5) + hex(:,6) )*0.5
                    else
                        c = ( a1(:) + a2(:) + a3(:) ) * 0.3333333333333333
                        ok = .false.
                    end if
                end if
            end if
            return
        end subroutine midpoint3skew

!-------



    end module Lib_TriangulatedSurface

! !
!     program testTriangulatedSurface
! !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
! !*        a quickie to test the functioning of    Lib_TriangulatedSurface
!         use Lib_TriangulatedSurface
!         use iso_fortran_env
!         implicit none
! 
! 
!         type(TriangulatedSurface)           ::      t
!         integer                             ::      n,ii
!         real(kind=real64),dimension(3)      ::      xx
!         
!    !---
!         print *,"enter number of points on sphere"
!         read(*,fmt='(i12)') n
!         t = GeodesicGrid( n )
!         call report(t)
!         print *,""
!         
!    !---
! 
!         open(unit=400,file="sphere.adv",action="write")
!              write (unit=400,fmt='(a)')  "# Surface"
!              write (unit=400,fmt='(i6,a)')   getn_triangle(t)," # number of tiles"
!              write (unit=400,fmt='(a)')  "atom position_x position_y position_z"
!              call report( t,400,.true. )
!         close(unit=400)
!         print *,""
!    !---
!         
!    !--- count number of nodes in 1/48th part
!         n = 0
!         do ii = 1,getn_node(t)
!             xx = getNode(t,ii)
!             if (PointInTriangle(xx, (/1.d0,0.d0,0.d0/),(/0.707107d0,0.707107d0,0.d0/),(/0.577350d0,0.577350d0,0.577350d0/))) then
!                 n = n + 1
!                 write (*,fmt='(a6,3f12.8)') "Cu ",xx
!             else
!                 write (*,fmt='(a6,3f12.8)') "Du ",xx
!             end if
!         end do
!         print *,"number ",n
!         print *,""
!    
!    !---
!    
!         call delete(t)
!    !---
!         print *,""
!         print *,"done"
!         print *,""
!         
!         
!     contains
! !---^^^^^^^^
! 
!         pure function SameSide(p1,p2, a,b) result(is)
!     !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
!     !*      tests if p1 is the same side of the line a-b as p2.
!             real(kind=real64),dimension(3),intent(in)       ::      p1,p2                                        
!             real(kind=real64),dimension(3),intent(in)       ::      a,b
!             logical                                         ::      is          
!             real(kind=real64),dimension(3)                  ::      cp1,cp2
!             cp1 = CrossProduct(b-a, p1-a)
!             cp2 = CrossProduct(b-a, p2-a)
!             is = ( cp1(1)*cp2(1) + cp1(2)*cp2(2) + cp1(3)*cp2(3) >= 0 )
!             return
!         end function SameSide
!         
!         pure function PointInTriangle(p, a,b,c) result(is)
!     !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
!     !*      returns true if point p is within the triangle a-b-c
!             real(kind=real64),dimension(3),intent(in)       ::      p                                   
!             real(kind=real64),dimension(3),intent(in)       ::      a,b,c
!             logical                                         ::      is          
!             
!             is = ( SameSide(p,a, b,c) .and. SameSide(p,b, a,c) .and. SameSide(p,c, a,b) )
!             return
!         end function PointInTriangle
!     
!         pure function CrossProduct(a,b) result(c)
!     !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
!             real(kind=real64),dimension(3),intent(in)       ::      a,b
!             real(kind=real64),dimension(3)                  ::      c
!             c(1) = a(2)*b(3) - a(3)*b(2)
!             c(2) = a(3)*b(1) - a(1)*b(3)
!             c(3) = a(1)*b(2) - a(2)*b(1)
!             return
!         end function CrossProduct
!             
!         
!         
! 
!     end program testTriangulatedSurface
