
    module Lib_Graphs
!---^^^^^^^^^^^^^^^^^^
!*      A data type capable of storing a small unweighted graph of up to 2^15 vertices and 2^15 edges

!       constuct a simple graph from the upper triangle matrix of edges
!       
!       eg  edge = transpose( reshape( (/   X,T,T,F,F ,         &
!                                           X,X,T,T,F ,         &                               
!                                           X,X,X,F,F ,         &                               
!                                           X,X,X,X,T ,         &                               
!                                           X,X,X,X,X       /),(/NV,NV/) ) )                       
!       g = graph_ctor( vertex,edge ) 
!
!       then you have
!           subgraph(g1,g2)
!           g1 == g2
!           hasVertex(g,v)
!           hasEdge(g,v1,v2)

!
!       Current version only supports undirected graphs
!       but there is flexibility in the storage format of edges to extend to directed if necessary.
!


        use Lib_Quicksort
        use iso_fortran_env
        implicit none
        
        private
        
    !---
    
        public      ::      Graph_ctor
        public      ::      delete
        public      ::      report
        public      ::      clone
        
        public      ::      operator(==)        
        public      ::      getNV,getNE
        public      ::      subgraph     
        public      ::      hasVertex 
        public      ::      hasEdge     
        public      ::      getEdge   
        public      ::      getVertex  
        public      ::      edgeCode
        public      ::      addVertex
        public      ::      addEdge
        public      ::      whichEdge
            
    !---
    
        logical,public          ::      Graph_dbg = .false.
                   
        integer,private,parameter   ::          MASK = 32767
        
    !---
    
        type,public     ::      Graph
            private
            integer                         ::      nVE             !   number of vertices in low 15 bits, number of edges in high 15 bits.
            integer,dimension(:),pointer    ::      vertex          !   a single 32 bit integer index for each vertex
            integer,dimension(:),pointer    ::      edge            !   low 15 bits vertex from, high 15 bits vertex to. 
        end type Graph
        
    !---
    
        interface Graph_ctor
            module procedure    Graph_null
            module procedure    Graph_ctor0
            module procedure    Graph_ctor1
        end interface
                
        interface delete
            module procedure    delete0
        end interface
        
        interface report
            module procedure    report0
        end interface
        
        interface clone
            module procedure    clone0
        end interface
        
        interface operator(==)
            module procedure    equals
        end interface
        
        interface getEdge
            module procedure    getEdge1            
        end interface
        
        interface whichEdge
            module procedure    whichEdge1
            module procedure    whichEdge2            
        end interface
        
    contains
!---^^^^^^^^

        function Graph_null() result(this)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
            type(Graph)           ::      this
            this%nVE = 0
            nullify(this%vertex)
            nullify(this%edge)
            return
        end function Graph_null
                         
        function Graph_ctor0(v,e) result(this)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      construct a simple undirected graph with an index for each vertex
    !*      and the upper triangle of e stores whether vertex i connects to vertex j
            integer,dimension(:),intent(in)     ::      v
            logical,dimension(:,:),intent(in)   ::      e
            type(Graph)                 ::      this
            
            integer         ::      ii,jj,nV,nE,indxi,indxj
            integer,dimension(size(v))  ::      indx
            nV = size(v)
            
        !---    count the number of edges
            nE = 0
            do jj = 2,nV
                nE = nE + count( e(1:jj-1,jj) )
            end do
            !print *,"nE = ",nE,count(e)
            
            this%nVE = nV + ishft(nE,16)
            
            allocate(this%vertex(nV))
            allocate(this%edge(nE))
            
        !---    add the vertices, sorting first
            do ii = 1,nV                 
                indx(ii) = ii
            end do
            call quicksort( v,indx )    !   now v(indx(2)) > v(indx(1))
            do ii = 1,nV
                jj = indx(ii)
                this%vertex(ii) = v(jj)
                !print *,"vertex ",ii," is ",jj,v(jj)
            end do
            
        !---    add the edges
            nE = 0
            do jj = 1,nV
                indxj = indx(jj)
                do ii = 1,nV
                    indxi = indx(ii)
                    
                    if ( indxi < indxj ) then
                        !print *,"considering i<j ",this%vertex(ii),this%vertex(jj),indxi,indxj,e(indxi,indxj)
                        if (e(indxi,indxj)) then
                            nE = nE + 1
                            this%edge(nE) = edgeCode(ii,jj)
                        end if
                    !else if ( indxj < indxi ) then
                    !   ! print *,"considering j<i ",this%vertex(ii),this%vertex(jj),indxj,indxi,e(indxj,indxi)
                    !   if (e(indxj,indxi)) then
                    !       nE = nE + 1
                    !       this%edge(nE) = edge(jj,ii)
                    !   end if
                    !
                    end if                  
                end do
            end do   
            
            return
        end function Graph_ctor0
                       
        function Graph_ctor1(v,e) result(this)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      construct an null graph but with storage
            integer,intent(in)     ::      v,e
            type(Graph)                 ::      this
            
            this%nVE = 0
            allocate(this%vertex(v))
            allocate(this%edge(e)) 
            this%vertex = 0
            this%edge = 0
            return
        end function Graph_ctor1
                       
    !---
    
        subroutine delete0(this)
    !---^^^^^^^^^^^^^^^^^^^^^^^^
            type(Graph),intent(inout)    ::      this
            if (this%nVE > 0) then
                deallocate(this%vertex)
                deallocate(this%edge)
            end if
            this = Graph_null()
            return
        end subroutine delete0
        
    !---
    
        subroutine report0(this,u,o)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^
            type(Graph),intent(in)          ::      this
            integer,intent(in),optional     ::      u,o
            integer     ::      uu,oo
            integer     ::      ii,jj,nV,nE
            uu = 6 ; if (present(u)) uu = u
            oo = 0 ; if (present(o)) oo = o
            
            nV = getnV(this)
            nE = getnE(this)
            write(unit=uu,fmt='(4(a,i5))') repeat(" ",oo)//"Graph [nV,nE = ",nV,",",nE,"]"            
            if (nV<=1) return
            
            
            if (maxval(this%vertex(1:nV)) < 100) then
                write(unit=uu,fmt='(a,10000i3)') repeat(" ",oo+8),this%vertex(1:nV)
                do ii = 1,nV
                    write(unit=uu,fmt='(a,i3)',advance="no") repeat(" ",oo+5),this%vertex(ii)
                    do jj = 1,nV
                        if (ii==jj) then                        
                            write(unit=uu,fmt='(a3)',advance="no") "   "                
                        else if (any( undirectedEdgeEquals( this%edge(1:nE), edgeCode(ii,jj)))) then
                            write(unit=uu,fmt='(a3)',advance="no") "  *"                
                        else
                            write(unit=uu,fmt='(a3)',advance="no") "  ."                
                        end if
                    end do
                    write(unit=uu,fmt='(a)',advance="yes") ""
                end do
            else if (maxval(this%vertex(1:nV)) < 1000) then
                write(unit=uu,fmt='(a,10000i4)') repeat(" ",oo+6),this%vertex(1:nV)
                do ii = 1,nV
                    write(unit=uu,fmt='(a,i4)',advance="no") repeat(" ",oo+2),this%vertex(ii)
                    do jj = 1,nV                        
                        if (ii==jj) then                        
                            write(unit=uu,fmt='(a4)',advance="no") "    "                
                        else if (any( undirectedEdgeEquals( this%edge(1:nE), edgeCode(ii,jj) ))) then
                            write(unit=uu,fmt='(a4)',advance="no") "  **"                   
                        else
                            write(unit=uu,fmt='(a4)',advance="no") "  .."                   
                        end if
                    end do
                    write(unit=uu,fmt='(a)',advance="yes") ""
                end do
            else                
                write(unit=uu,fmt='(a,10000i6)') repeat(" ",oo+8),this%vertex(1:nV)
                do ii = 1,nV
                    write(unit=uu,fmt='(a,i6)',advance="no") repeat(" ",oo+2),this%vertex(ii)
                    do jj = 1,nV                        
                        if (ii==jj) then                        
                            write(unit=uu,fmt='(a6)',advance="no") "      "                
                        else if (any( undirectedEdgeEquals( this%edge(1:nE), edgeCode(ii,jj)))) then
                            write(unit=uu,fmt='(a6)',advance="no") "    **"                 
                        else
                            write(unit=uu,fmt='(a6)',advance="no") "    .."                 
                        end if
                    end do
                    write(unit=uu,fmt='(a)',advance="yes") ""
                end do
            end if
            return
        end subroutine report0
    
    !---
    
        subroutine clone0(this,that) 
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      deep copy with allocation this = that
            type(Graph),intent(inout)   ::      this
            type(Graph),intent(in)      ::      that
            
            integer         ::      thisnV,thatnV,thisnE,thatnE
                      
            thisnV = 0  
            if (associated(this%vertex))  thisnV = size(this%vertex)
            thisnE = 0
            if (associated(this%edge))  thisnE = size(this%edge)
             
            thatnV = getNV( that )
            thatnE = getNE( that )
            
            if (thisnV == 0) then   
                !   fresh allocation
                allocate(this%vertex(thatnV))                
            else if (thisnV < thatnV) then
                !   reallocation
                deallocate(this%vertex)
                allocate(this%vertex(thatnV))                
            end if
            this%vertex(1:thatnV) = that%vertex(1:thatnV)
             
            if (thisnE == 0) then   
                !   fresh allocation
                allocate(this%edge(thatnE))                
            else if (thisnE < thatnE) then
                !   reallocation
                deallocate(this%edge)
                allocate(this%edge(thatne))                
            end if
            this%edge(1:thatnE) = that%edge(1:thatnE)
            
            this%nVE = that%nVE
            return
        end subroutine clone0
            
    !---    
    
     
        subroutine addVertex( this,v )
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      if vertex v has not previously been seen, add it to the graph
            type(Graph),intent(inout)       ::      this
            integer,intent(in)              ::      v
            integer,dimension(:),pointer    ::      vertex_tmp
                
            integer         ::      nV
            nV = getNV( this )
            if (any(this%vertex(1:nV) == v)) return
            
            if (.not. associated(this%vertex)) then
                allocate(this%vertex(1:10))
            else if (nV == size(this%vertex)) then
                allocate(vertex_tmp(1:2*nV))
                vertex_tmp(1:nV) = this%vertex(1:nV)
                deallocate(this%vertex)
                this%vertex => vertex_tmp
            end if
            this%nVE = this%nVE + 1     !   add one vertex, don't change edges.
            this%vertex(nV+1) = v
            
            return
        end subroutine addVertex
        
        
        subroutine addEdge( this,g1,g2 )
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      if edge g1-g2 has not previously been seen, add it to the graph
            type(Graph),intent(inout)       ::      this
            integer,intent(in)              ::      g1,g2
            integer,dimension(:),pointer    ::      edge_tmp
                
            integer         ::      v1,v2,ii,nE,ee
             
            
        !---    check the vertices are present. Add if not.
            v1 = whichVertex( this,g1 )
            if (v1 == 0) then
                call addVertex( this,g1 )
                v1 = getNV(this)
            end if
            v2 = whichVertex( this,g2 )
            if (v2 == 0) then
                call addVertex( this,g2 )
                v2 = getNV(this)
            end if
            
            
        !---    have we already got this edge?
            nE = getNE( this )
            ee = edgeCode(v1,v2)             
            do ii = 1,nE
                if (undirectedEdgeEquals( this%edge(ii),ee   )) return
            end do
            
            
        !---    we need to add it    
            if (.not. associated(this%edge)) then
                allocate(this%edge(1:10))
            else if (nE == size(this%edge)) then
                allocate(edge_tmp(1:2*nE))
                edge_tmp(1:nE) = this%edge(1:nE)
                deallocate(this%edge)
                this%edge => edge_tmp
            end if
            this%nVE = this%nVE + ishft(1,16)     !   add one edge, don't change vertices
            this%edge(nE+1) = ee
            
            
            return
        end subroutine addEdge
    !---
    
        elemental function equals( this,that ) result(is)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      note that two graphs can only be the same if their vertices are the same.
    !*      and that vertices are ordered, so that the arrays must be identical
    !*      but this also means that the edges must have the same values. 
    !*      If you use the constructor properly, they should also be in the same order.
            type(Graph),intent(in)      ::      this,that
            logical                     ::      is
            integer             ::      ii
            
            is = (this%nVE == that%nVE)
            if (.not. is) return
            do ii = 1,getNV( this )
                is = is .and. ( this%vertex(ii) == that%vertex(ii) )
            end do
            if (.not. is) return
            do ii = 1,getNE( this )
                is = is .and. undirectedEdgeEquals( this%edge(ii),that%edge(ii) )
            end do
            return
        end function equals
        
        
          function subgraph(this,that) result(is)   
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      returns true if that is a subgraph of this
            type(Graph),intent(in)      ::      this,that
            logical                     ::      is
            
            integer         ::      thisnV,thisnE,thatnV,thatnE,ee
            integer         ::      ii,jj,v1,v2
            
        
            
            thisnV = getNV( this )
            thisnE = getNE( this )
            thatnV = getNV( that )
            thatnE = getNE( that )
            
            
            is = .false.
            if ( (thisnV<thatnV).or.(thisnE<thatnE) ) return
                 
            
            jj = 1                  !   search position in this list
            do ii = 1,thatnV        !   search for vertex i in this list
                do
                    if (this%vertex(jj) == that%vertex(ii)) exit        !   found it.
                    jj = jj + 1
                    if (jj>thisnV) return                               !   have got to the end of the list, and didn't find it.
                end do
            end do
                            
            
            jj = 1                  !   search position in this list
            do ii = 1,thatnE        !   search for edge i in that list
                call getEdge( that,ii,v1,v2 )       !   vertices for that edge
                ee = findEdge( this,v1,v2 )         !   corresponding edge in this
                if (ee == 0) return                 !   edge not found in this
             
            end do
                
            is = .true.
                
            return
        end function subgraph
            
    !---
    
        elemental function edgeCode( v1,v2 ) result(e)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
            integer,intent(in)              ::      v1,v2
            integer                         ::      e
            e = v1 + ishft(v2,16)
            return
        end function edgeCode
        
    !---
    
        elemental function hasVertex( this,v ) result(is)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^                           
            type(Graph),intent(in)      ::      this
            integer,intent(in)          ::      v
            logical                     ::      is
            integer         ::      ii
            is = .false.
            do ii = 1,getNv(this)
                if (v == this%vertex(ii)) then
                    is = .true.
                    exit
                else if (v < this%vertex(ii)) then
                    exit
                end if
            end do
            return
        end function hasVertex    
         
        elemental function whichVertex( this,v ) result(i)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      return the vertex number 1:nV matching v. return 0 if not found.                        
            type(Graph),intent(in)      ::      this
            integer,intent(in)          ::      v
            integer                     ::      i 
            integer         ::      ii
            i = 0
            do ii = 1,getNv(this)
                if (v == this%vertex(ii)) then
                    i = ii
                    exit
                else if (v < this%vertex(ii)) then
                    exit
                end if
            end do
            return
        end function whichVertex     
        
        
        elemental function hasEdge( this,v1,v2 ) result(is)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^                         
            type(Graph),intent(in)      ::      this
            integer,intent(in)          ::      v1,v2
            logical                     ::      is
            integer         ::      i1,i2,ii,ee
            
            i1 = whichVertex( this,v1 )
            i2 = whichVertex( this,v2 )
            is = .false.
            
            if (i1*i2 == 0) return      !   one or other vertex not present, so edge can't be.
            
            ee = edgeCode(i1,i2)
            do ii = 1,getNe(this)
                if (undirectedEdgeEquals( this%edge(ii),ee   )) then
                    is = .true.
                    exit
                end if
            end do
            
            return
        end function hasEdge                            
                                            
                               
        pure subroutine getEdge1( this,i,v1,v2 ) 
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      returns the vertices associated with edge i
            type(Graph),intent(in)      ::      this
            integer,intent(in)          ::      i
            integer,intent(out)         ::      v1,v2
            
            v1 = iand(MASK,this%edge(i))
            v2 = ishft(this%edge(i),-16)
            
            v1 = this%vertex(v1)
            v2 = this%vertex(v2)
            return
        end subroutine getEdge1                        
                          
        elemental function getVertex( this,i ) result(v) 
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      returns the ith vertex
            type(Graph),intent(in)      ::      this
            integer,intent(in)          ::      i
            integer                     ::      v
            
            v = this%vertex(i)
            return
        end function getVertex                       
                          
        
        pure function findEdge( this,v1,v2 ) result(e) 
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      returns the edge associated with vertices v1,v2.
    !*      returns 0 if not found.
            type(Graph),intent(in)      ::      this             
            integer,intent(in)          ::      v1,v2
            integer                     ::      e
            integer         ::      i1,i2,ii
            
            e = 0
            i1 = whichVertex( this,v1 )
            i2 = whichVertex( this,v2 )          
            if (i1*i2 == 0) return      !   one or other vertex not present, so edge can't be.
            
            !   lets assume it is present, and check for it.
            e = edgeCode(i1,i2)
            do ii = 1,getNe(this)
                if (undirectedEdgeEquals( this%edge(ii),e )) return             
            end do
            
            !   failed to find it.
            e = 0
            return
        end function findEdge                      
                                    
        pure function whichEdge1( this,e ) result(i) 
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      returns the edge 
    !*      returns 0 if not found.
            type(Graph),intent(in)      ::      this             
            integer,intent(in)          ::      e
            integer                     ::      i 
             
            do i = 1,getNe(this)
                if (undirectedEdgeEquals( this%edge(i),e )) return             
            end do
            i = 0
            return
        end function whichEdge1                    
                                           
        pure function whichEdge2( this,v1,v2 ) result(i) 
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      returns the edge 
    !*      returns 0 if not found.
            type(Graph),intent(in)      ::      this             
            integer,intent(in)          ::      v1,v2
            integer                     ::      i
            
            integer         ::      i1,i2,ee
            
            i = 0
            i1 = whichVertex( this,v1 )
            i2 = whichVertex( this,v2 )          
            if (i1*i2 == 0) return      !   one or other vertex not present, so edge can't be.
            
            !   lets assume it is present, and check for it.
            ee = edgeCode(i1,i2)
            do i = 1,getNe(this)
                if (undirectedEdgeEquals( this%edge(i),ee )) return             
            end do
            
            !   failed to find it.
            i = 0
            return
        end function whichEdge2                              
                               
        
    !---
    
        elemental function getNV( this ) result(nV)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ 
            type(Graph),intent(in)      ::      this
            integer                     ::      nV
            nV = iand( MASK,this%nVE )
            return
        end function getNV
        
        elemental function getNE( this ) result(nE)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ 
            type(Graph),intent(in)      ::      this
            integer                     ::      nE
            nE = iand( MASK,ishft(this%nVE,-16) )
            return
        end function getNE
        
        elemental function reverseEdge( e ) result( f )
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ 
            integer,intent(in)      ::      e
            integer                 ::      f            
            f = ishft( iand( MASK,e ),16 ) + ishft( e,-16 )
            return
        end function reverseEdge
        
        elemental function undirectedEdgeEquals( e,f ) result( is )
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ 
            integer,intent(in)      ::      e,f     
            logical                 ::      is
            is = (e == f) .or. ( e == reverseEdge( f ) ) 
            return
        end function undirectedEdgeEquals
        
    end module Lib_Graphs
    
    