
    module GrainNodes
!---^^^^^^^^^^^^^^^^^
        use Lib_Graphs
        use Lib_FitTanhFunction
        use Lib_MarchingCubes
        use SolveTripleJunction
        use iso_fortran_env
        implicit none
        private
        include "permutations.h"
        
    !---
    
        public      ::      GrainNode_ctor
        public      ::      delete
        public      ::      report
     
        public      ::      getPhi
        public      ::      cubeProperties
        public      ::      getnGrainBoundaries
        public      ::      getnTripleJunctions
        public      ::      getnQuadJunctions
        public      ::      getn5Junctions
        public      ::      getn6Junctions
        
        
    !---
    
        logical,public          ::      GrainNode_dbg = .false.
                   
         
           
        type(Graph),public                              ::      fullGraph     
        integer,public                                  ::      nClique3 = 0
        type(Graph),dimension(:),pointer,public         ::      clique3  
        integer,public                                  ::      nClique4 = 0
        type(Graph),dimension(:),pointer,private        ::      clique4  
        integer,public                                  ::      nClique5 = 0
        type(Graph),dimension(:),pointer,private        ::      clique5  
        integer,public                                  ::      nClique6 = 0
        type(Graph),dimension(:),pointer,private        ::      clique6  
        
        
        
    !---
    
        type,public     ::      GrainNode
            private
            integer             ::      nGrains
            integer             ::      nGrainBoundaries
            integer             ::      nTripleJunctions
            integer             ::      nQuadJunctions
            integer             ::      n5Junctions
            integer             ::      n6Junctions
            type(Graph)         ::      g            
            real(kind=real64),dimension(:),pointer      ::      phi     !   (1:nGrainBoundaries) = distance from grain boundary with -1 at from grain and +1 at to grain and 0 at GB      
            integer             ::      grain                           !   if there are no grain boundaries, then this will store the actual grain this node is in. 
            integer,dimension(:),pointer                ::      tripleJunction
            integer,dimension(:),pointer                ::      quadJunction
            integer,dimension(:),pointer                ::      fiveJunction
            integer,dimension(:),pointer                ::      sixJunction
        end type GrainNode
        
    !---
    
        interface GrainNode_ctor
            module procedure    GrainNode_null
            module procedure    GrainNode_ctor0
        end interface
                
        interface delete
            module procedure    delete0
        end interface
        
        interface report
            module procedure    report0
        end interface
        
        interface   getPhi
            module procedure    getPhi0
        end interface
        
    contains
!---^^^^^^^^

        function GrainNode_null() result(this)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
            type(GrainNode)           ::      this
            this%nGrains = 0
            this%nGrainBoundaries = 0
            this%nTripleJunctions = 0
            this%nQuadJunctions   = 0
            this%n5Junctions      = 0
            this%n6Junctions      = 0
            nullify(this%phi)
            nullify(this%tripleJunction)
            nullify(this%quadJunction  )
            nullify(this%fiveJunction  )
            nullify(this%sixJunction   )
            this%g = Graph_ctor()
            return
        end function GrainNode_null
                         
        function GrainNode_ctor0(g,x,rc) result(this)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
            type(GrainNode)                 ::      this
            integer,dimension(:),intent(in)                 ::      g
            real(kind=real64),dimension(:,:),intent(in)     ::      x
            real(kind=real64),intent(in)                    ::      rc
            
            integer,dimension(0:size(g))        ::      vertex
            logical,dimension(:,:),allocatable  ::      edge
            integer,dimension(:),allocatable    ::      indx
            real(kind=real64),dimension(3)      ::      dx
            real(kind=real64)                   ::      dd,ddsum,ddmax
            integer             ::      ii,jj,kk,nn
            integer             ::      v1,v2,v3,v4,v5,v6
            integer             ::      g1,g2
            
            real(kind=real64),dimension(:,:),allocatable  ::      yy
            real(kind=real64),dimension(:),allocatable    ::      ff
            
            type(Graph)                         ::      clique
            integer,dimension(:,:),pointer      ::      list
            logical                             ::      ok
            type(Graph),dimension(:),pointer    ::      clique_tmp
            
            integer,dimension(1000)              ::      junction
            
            
        !---    I have a set of atoms at relative positions to the node x
        !       and grain indices g.
        !       use this to construct a graph with edges where atoms of different grains are within rc.
            
            nn = size(g)
            allocate(indx(0:maxval(g))) ; indx = 0
            this = GrainNode_null()
            
        !   first count different grains
            this%nGrains = 0
            vertex = 0
            do ii = 1,nn
                if (g(ii) == 0) cycle                               !   not interested in "unset" grain
                if (any(vertex(1:this%nGrains)==g(ii))) cycle       !   have already seen this grain
                this%nGrains = this%nGrains + 1
                vertex(this%nGrains) = g(ii)
                indx(g(ii)) = this%nGrains                
            end do
            
            if (this%nGrains>8) then
                print *,"GrainNodes::GrainNode_ctor0 warning - large number of grains detected at one node nGrains=",this%nGrains                
                !GrainNode_dbg = .true.
            end if
            
            
            if (GrainNode_dbg) print *,"GrainNodes::GrainNode_ctor0 dbg - number of atoms,highest grain index ",nn,maxval(g) 
            if (GrainNode_dbg) print *,"GrainNodes::GrainNode_ctor0 dbg - number different grains ",this%nGrains," grains ",vertex(1:this%nGrains)
            
            
            
        !   now count edges    
            allocate(edge(this%nGrains,this%nGrains)) ; edge = .false.
            do ii = 1,nn 
                if (g(ii) == 0) cycle                   !   not interested in "unset" grain
                 
                v1 = indx(g(ii))
                do jj = 1,nn
                    
                    !if (g(jj)*(g(ii)-g(jj))==0) cycle   !   not interested in "unset" grain or same grain as i    
                    v2 = indx(g(jj))
                    if (v1>=v2) cycle                    !   only need to complete top half of edge matrix - note this also removes "unset" grain or same grain as i                   
                     
                    !   sort by distance                                
                      dx = x(1:3,jj) - x(1:3,ii)
                      dd = maxval(abs(dx))
                      if (dd<=2*rc) edge( v1,v2 ) = .true.
                    
                    ! edge( v1,v2 ) = .true.
                    
                end do
            end do
            if (GrainNode_dbg) print *,"GrainNodes::GrainNode_ctor0 dbg - number grain boundaries ",count(edge)
            
            
        !   can now construct graph for this node.
            this%g = Graph_ctor( vertex(1:this%nGrains),edge )            
            deallocate(edge)
            
            if (GrainNode_dbg) print *,"GrainNodes::GrainNode_ctor0 dbg - partially constructed grain node- graph only"
            if (GrainNode_dbg) call report(this%g) 
            
            
        !   establish most likely grain type.
            if (this%nGrains == 0) then
                this%grain = 0
            else
            
                !   which grain has highest weight on the node?
                ddmax = 0.0d0 ; kk = 0
                do ii = 1,this%nGrains
                    ddsum = 0.0d0
                    do jj = 1,nn
                        if (g(jj) == vertex(ii)) then   
                            dd = x(1,jj)*x(1,jj) + x(2,jj)*x(2,jj) + x(3,jj)*x(3,jj)
                            dd = dd / (2*rc*rc/4)       !   use rc/2 as weighting length
                            dd = exp( - dd )            !   gaussian weight
                            ddsum = ddsum + dd
                        end if
                    end do
                    if (GrainNode_dbg) print *,"GrainNodes::GrainNode_ctor0 dbg - grain ",vertex(ii)," count ",count(g(1:nn)==vertex(ii))," weight ",ddsum 
                    if (ddsum > ddmax) then
                        kk = ii ; ddmax = ddsum
                    end if
                end do
                this%grain = vertex(kk)
                    
                
            
            end if                 
            if (GrainNode_dbg) print *,"GrainNodes::GrainNode_ctor0 dbg - this grain ",this%grain
            
            
        !   can find grain boundary count = number of edges, and can find projections of each grain at the node
            this%nGrainBoundaries = getNe(this%g)    
            allocate(this%phi(this%nGrainBoundaries))
            this%phi = -1
            allocate(yy(3,nn))
            allocate(ff(nn))
            do kk = 1,this%nGrainBoundaries
                call getEdge( this%g,kk,g1,g2 )         !   consider the atoms in grains g1 and g2...
                call addEdge( fullGraph,g1,g2 )         !   make sure they are in the full graph too.
                
                !print *,"find projection edge ",kk," grains ",g1,g2
                jj = 0
                do ii = 1,nn
!                   if (g(ii) == 0) then
!                         jj = jj + 1
!                         ff(jj) = 0                        !   set grain boundary to value 0
!                         yy(1:3,jj) = x(1:3,ii)                            
!                     else 
                    if (g(ii) == g1) then
                        jj = jj + 1
                        ff(jj) = +1                     !   set grain v1 to value +1
                        yy(1:3,jj) = x(1:3,ii)
                    else if (g(ii) == g2) then
                        jj = jj + 1
                        ff(jj) = -1                     !   set grain v2 to value -1
                        yy(1:3,jj) = x(1:3,ii)
                    end if
                end do
                
                !call fitTanhFunction( yy(1:3,1:jj),ff(1:jj), dx,dd )    !   find a logistic function to interpolate between grains
                call fitTanhFunction( yy(1:3,1:jj),ff(1:jj), dx,dd,rc )    !   find a logistic function to interpolate between grains
                    
                if (GrainNode_dbg) print *,"GrainNodes::GrainNode_ctor0 dbg - edge ",g1,g2," has ",jj," atoms. Phi = ",tanh(dd) 
                
                
                this%phi(kk) = tanh(dd)                                 !   this is the value of the field at the node.
                  
            end do   
            deallocate(yy)
            deallocate(ff)
            
            
                  
               
        !   now can look for triple point junctions - these are three-vertex cliques in the graph
            this%nTripleJunctions = 0
            if ( (this%nGrains >=3).and.(this%nGrainBoundaries >=3) ) then
                !   possibly. Lets look.    
                allocate(edge(3,3)) ; edge = .true.
                call permutations( this%nGrains,3, list )
                do ii = 1,size(list,dim=2)
                !---    construct a subgraph with 3 vertices permuted from the list of grains, fully connected
                    v1 = getVertex(this%g,list(1,ii))
                    v2 = getVertex(this%g,list(2,ii))
                    v3 = getVertex(this%g,list(3,ii))
                    clique = Graph_ctor( (/v1,v2,v3/) , edge )
                    
                !---    is this connected subgraph within the graph on the node? 
                    ok = subgraph( this%g,clique )
                    if (ok) then
                        this%nTripleJunctions = this%nTripleJunctions + 1
                        
                    !---    is this a new clique?
                        ok = .false.
                        do jj = 1,nClique3
                            if (clique3(jj) == clique) then
                                ok = .true.
                                junction(this%nTripleJunctions) = jj
                                exit
                            end if
                        end do
                        
                        if (.not. ok) then  !   it is a new one
                        !---    reallocate space if necessary
                            if (nClique3 == 0) then
                                allocate(clique3(10))                            
                            else if (nClique3 == size(clique3)) then
                                allocate(clique_tmp(nClique3*2))
                                do jj = 1,nClique3
                                    clique_tmp(jj) = Graph_ctor()
                                    call clone( clique_tmp(jj),clique3(jj) )
                                    call delete(clique3(jj))
                                end do
                                deallocate(clique3)
                                clique3 => clique_tmp                            
                            end if
                                
                        !---    add triple junction to list
                            nClique3 = nClique3 + 1
                            clique3(nClique3) = Graph_ctor()
                            call clone( clique3(nClique3),clique )
                            junction(this%nTripleJunctions) = nClique3
                            
                            if (GrainNode_dbg) then
                                print *,"GrainNodes::GrainNode_ctor0() info - added new triple junction ",nClique3
                                call report(clique3(nClique3))
                            end if
                            
                        end if                        
                        
                    end if
                    call delete(clique)
                end do
                deallocate(edge)
            end if
            
            if (this%nTripleJunctions>0) then
                allocate(this%tripleJunction(this%nTripleJunctions))
                this%tripleJunction(1:this%nTripleJunctions) = junction(1:this%nTripleJunctions)
            end if
            
               
        !   now can look for four point junctions - these are four-vertex cliques in the graph

        !   now can look for quad point junctions - these are four-vertex cliques in the graph
            this%nQuadJunctions = 0
            if ( (this%nGrains >=4).and.(this%nGrainBoundaries >=4) ) then
                !   possibly. Lets look.    
                allocate(edge(4,4)) ; edge = .true.
                call permutations( this%nGrains,4, list )
                 
                do ii = 1,size(list,dim=2)
                !---    construct a subgraph with 3 vertices permuted from the list of grains, fully connected
                    v1 = getVertex(this%g,list(1,ii))
                    v2 = getVertex(this%g,list(2,ii))
                    v3 = getVertex(this%g,list(3,ii))
                    v4 = getVertex(this%g,list(4,ii))
                    clique = Graph_ctor( (/v1,v2,v3,v4/) , edge )
                    
                !---    is this connected subgraph within the graph on the node? 
                    ok = subgraph( this%g,clique )
                    if (ok) then
                        this%nQuadJunctions = this%nQuadJunctions + 1
                        
                    !---    is this a new clique?
                        ok = .false.
                        do jj = 1,nClique4
                            if (clique4(jj) == clique) then
                                ok = .true.
                                junction(this%nQuadJunctions) = jj
                                exit
                            end if
                        end do
                        
                        if (.not. ok) then  !   it is a new one
                        !---    reallocate space if necessary
                            if (nClique4 == 0) then
                                allocate(clique4(10))                            
                            else if (nClique4 == size(clique4)) then
                                allocate(clique_tmp(nClique4*2))
                                do jj = 1,nClique4
                                    clique_tmp(jj) = Graph_ctor()
                                    call clone( clique_tmp(jj),clique4(jj) )
                                    call delete(clique4(jj))
                                end do
                                deallocate(clique4)
                                clique4 => clique_tmp                            
                            end if
                                
                        !---    add quad junction to list
                            nClique4 = nClique4 + 1
                            clique4(nClique4) = Graph_ctor()
                            call clone( clique4(nClique4),clique )
                            junction(this%nQuadJunctions) = nClique4
                            
                            if (GrainNode_dbg) then
                                print *,"GrainNodes::GrainNode_ctor0() info - added new quad junction ",nClique4
                                call report(clique4(nClique4))
                            end if
                            
                        end if                        
                        
                    end if
                    call delete(clique)
                end do
                deallocate(edge)
            end if
            
            if (this%nQuadJunctions>0) then
                allocate(this%quadJunction(this%nQuadJunctions))
                this%quadJunction(1:this%nQuadJunctions) = junction(1:this%nQuadJunctions)
            end if
                        
            !return
            
            
        !   now can look for five point junctions - these are five-vertex cliques in the graph
            this%n5Junctions = 0
            if ( (this%nGrains >=5).and.(this%nGrainBoundaries >=5) ) then
                !   possibly. Lets look.    
                allocate(edge(5,5)) ; edge = .true.
                call permutations( this%nGrains,5, list )
                 
                do ii = 1,size(list,dim=2)
                !---    construct a subgraph with 4 vertices permuted from the list of grains, fully connected
                    v1 = getVertex(this%g,list(1,ii))
                    v2 = getVertex(this%g,list(2,ii))
                    v3 = getVertex(this%g,list(3,ii))
                    v4 = getVertex(this%g,list(4,ii))
                    v5 = getVertex(this%g,list(5,ii))
                    clique = Graph_ctor( (/v1,v2,v3,v4,v5/) , edge )
                    
                !---    is this connected subgraph within the graph on the node? 
                    ok = subgraph( this%g,clique )
                    if (ok) then
                        this%n5Junctions = this%n5Junctions + 1
                        
                    !---    is this a new clique?
                        ok = .false.
                        do jj = 1,nClique5
                            if (clique5(jj) == clique) then
                                ok = .true.
                                junction(this%n5Junctions) = jj
                                exit
                            end if
                        end do
                        
                        if (.not. ok) then  !   it is a new one
                        !---    reallocate space if necessary
                            if (nClique5 == 0) then
                                allocate(clique5(10))                            
                            else if (nClique5 == size(clique5)) then
                                allocate(clique_tmp(nClique5*2))
                                do jj = 1,nClique5
                                    clique_tmp(jj) = Graph_ctor()
                                    call clone( clique_tmp(jj),clique5(jj) )
                                    call delete(clique5(jj))
                                end do
                                deallocate(clique5)
                                clique5 => clique_tmp                            
                            end if
                                
                        !---    add quad junction to list
                            nClique5 = nClique5 + 1
                            clique5(nClique5) = Graph_ctor()
                            call clone( clique5(nClique5),clique )
                            junction(this%n5Junctions) = nClique5
                            
                            if (GrainNode_dbg) then
                                print *,"GrainNodes::GrainNode_ctor0() info - added new 5 junction ",nClique5
                                call report(clique5(nClique5))
                            end if
                            
                        end if                        
                        
                    end if
                    call delete(clique)
                end do
                deallocate(edge)
            end if
            
            if (this%n5Junctions>0) then
                allocate(this%fiveJunction(this%n5Junctions))
                this%fiveJunction(1:this%n5Junctions) = junction(1:this%n5Junctions)
            end if
                        
            
        !   now can look for six point junctions - these are six-vertex cliques in the graph
            this%n6Junctions = 0
            if ( (this%nGrains >=6).and.(this%nGrainBoundaries >=6) ) then
                !   possibly. Lets look.    
                allocate(edge(6,6)) ; edge = .true.
                call permutations( this%nGrains,6, list )
                !if (.not. associated(list)) return
                do ii = 1,size(list,dim=2)
                !---    construct a subgraph with 4 vertices permuted from the list of grains, fully connected
                    v1 = getVertex(this%g,list(1,ii))
                    v2 = getVertex(this%g,list(2,ii))
                    v3 = getVertex(this%g,list(3,ii))
                    v4 = getVertex(this%g,list(4,ii))
                    v5 = getVertex(this%g,list(5,ii))
                    v6 = getVertex(this%g,list(6,ii))
                    clique = Graph_ctor( (/v1,v2,v3,v4,v5,v6/) , edge )
                    
                !---    is this connected subgraph within the graph on the node? 
                    ok = subgraph( this%g,clique )
                    if (ok) then
                        this%n6Junctions = this%n6Junctions + 1
                        
                    !---    is this a new clique?
                        ok = .false.
                        do jj = 1,nClique6
                            if (clique6(jj) == clique) then
                                ok = .true.
                                junction(this%n6Junctions) = jj
                                exit
                            end if
                        end do
                        
                        if (.not. ok) then  !   it is a new one
                        !---    reallocate space if necessary
                            if (nClique6 == 0) then
                                allocate(clique6(10))                            
                            else if (nClique6 == size(clique6)) then
                                allocate(clique_tmp(nClique6*2))
                                do jj = 1,nClique6
                                    clique_tmp(jj) = Graph_ctor()
                                    call clone( clique_tmp(jj),clique6(jj) )
                                    call delete(clique6(jj))
                                end do
                                deallocate(clique6)
                                clique6 => clique_tmp                            
                            end if
                                
                        !---    add 6 junction to list
                            nClique6 = nClique6 + 1
                            clique6(nClique6) = Graph_ctor()
                            call clone( clique6(nClique6),clique )
                            junction(this%n6Junctions) = nClique6
                            
                            if (GrainNode_dbg) then
                                print *,"GrainNodes::GrainNode_ctor0() info - added new 6 junction ",nClique6
                                call report(clique6(nClique6))
                            end if
                            
                        end if                        
                        
                    end if
                    call delete(clique)
                end do
                deallocate(edge)
            end if
            
            if (this%n6Junctions>0) then
                allocate(this%sixJunction(this%n6Junctions))
                this%sixJunction(1:this%n6Junctions) = junction(1:this%n6Junctions)
            end if
                        
            
            
            return
        end function GrainNode_ctor0
                       
    !---
    
        subroutine delete0(this)
    !---^^^^^^^^^^^^^^^^^^^^^^^^
            type(GrainNode),intent(inout)    ::      this
            if (this%nGrains == 0) return
            deallocate(this%phi)
            call delete(this%g)
            
            if (this%nTripleJunctions>0) deallocate(this%tripleJunction)
            if (this%nQuadJunctions>0)   deallocate(this%quadJunction)
            if (this%n5Junctions>0)      deallocate(this%fiveJunction)
            if (this%n6Junctions>0)      deallocate(this%sixJunction)
                
            
            this = GrainNode_null()
            return
        end subroutine delete0
        
    !---
    
        subroutine report0(this,u,o)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^
            type(GrainNode),intent(in)    ::      this
            integer,intent(in),optional     ::      u,o
            integer     ::      uu,oo
            integer     ::      kk,g1,g2
            
            uu = 6 ; if (present(u)) uu = u
            oo = 0 ; if (present(o)) oo = o
            write(unit=uu,fmt='(8(a,i3))') repeat(" ",oo)//"GrainNode [grainID,grains,GB,TJ,QJ,5J,6J = ",this%grain,","       &
                    ,this%nGrains,",",this%nGrainBoundaries,",",this%nTripleJunctions,",",this%nQuadJunctions,",",this%n5Junctions,",",this%n6Junctions,"]"            
            call report(this%g,uu,oo+2)
            do kk = 1,this%nGrainBoundaries
                call getEdge( this%g,kk,g1,g2 )         !   consider the atoms in grains g1 and g2...
                write(unit=uu,fmt='(a,i6,a,i6,a,f10.5)') repeat(" ",oo+2)//"GB  ",g1,"-",g2," phi=",this%phi(kk)  
            end do
             
            do kk = 1,getnV(this%g)
                g1 = getVertex( this%g,kk )
                write(unit=uu,fmt='(a,i6 ,f10.5)') repeat(" ",oo+2)//"phi ",g1,getPhi( this,g1 )
            end do
            return
        end subroutine report0
     
    !---
    
        elemental function getPhi0(this,g) result(phi)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      return the phase field value for this node point for grain g
    !*      returns -1 for inside grain, +1 for outside grain and 0 for grain boundary.
            type(GrainNode),intent(in)     ::      this
            integer,intent(in)             ::      g
            real(kind=real64)              ::      phi
            
            integer         ::      kk,g1,g2
            logical         ::      ok 
             
            
            if (hasVertex( this%g,g )) then
                !   this node at least recognises the grain
                phi = +1.0d0 
                !   are we on the grain boundary? 
                if (this%nGrainBoundaries > 0) then
                    ok = .false.                !   recognise grain, but do we recognise GB?          
                    do kk = 1,this%nGrainBoundaries
                        call getEdge( this%g,kk,g1,g2 )
                        if (g1 == g) then
                            phi = min( phi,this%phi(kk) )
                            ok = .true.
                        else if (g2 == g) then
                            phi = min( phi,-this%phi(kk) )
                            ok = .true.
                        end if
                    end do
                    if (.not. ok) phi = -1.0d0      !   nope. There might be one or two atoms of grain g nearby, but not in range to pick up the grain boundary.
                else
                    if (g /= this%grain) then
                        phi = -1.0d0
                    end if
                end if
            else
                !   this node doesn't recognise the grain - must be way outside
                phi = -1.0d0
            end if
            return
        end function getPhi0
            
    !---
    
        subroutine cubeProperties( node,i0,a , Garea,Gvol,TJlen,GBarea,nTJline,TJline,TJlineId,nQJpoint , nTriOut,TriGrain,TriOut)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      given eight nodes with node 1 at (cell) offset i0 and x,y,z axis given by a(:,1),a(:,2),a(:,3) respectively
    !*          8-------7
    !*         /|      /|
    !*        / |     / |       z
    !*       5-------6  |       | y
    !*       |  4----|--3       |/
    !*       | /     | /        *---x
    !*       |/      |/
    !*       1-------2
    !*   
    !*      find the properties of grains, grain boundaries etc.
    !*      note that the offset is needed to accumulate area,volume properly
    !*      Note too that the cumulative volume, area is only sensible if the surface is closed- 
    !*      ... which may not be true for a periodic boundary. 
    !*      Make sure you the input offset starts at zero, and continues to the maximum level needed
    !*      rather than having some elements at i0 = (0,0,0) and some at (N,N,N)
       
            type(GrainNode),dimension(8),intent(in)             ::      node
            integer,dimension(3),intent(in)                     ::      i0
            real(kind=real64),dimension(3,3),intent(in)         ::      a
            real(kind=real64),dimension(:),intent(out)          ::      Garea    
            real(kind=real64),dimension(:),intent(out)          ::      Gvol  
            real(kind=real64),dimension(:),intent(out)          ::      TJlen 
            real(kind=real64),dimension(:),intent(out)          ::      GBarea  
            integer,intent(out)                                 ::      nTJline
            real(kind=real64),dimension(:,:,:),intent(out)      ::      TJline              !   (3,2,nLineMax)
            integer,dimension(:),intent(out)                    ::      TJlineId
            integer,intent(out)                                 ::      nQJpoint
            
            integer,intent(out),optional                                            ::      nTriOut        !   (1:nGrains)
            integer,dimension(:),allocatable,intent(out),optional                   ::      TriGrain        !   (1:nTriOut)
            real(kind=real64),dimension(:,:,:),allocatable,optional,intent(out)     ::      TriOut        !   (1:3 xyz,1:3 vertex, 1:nTriOut)
                                                                                                     
            integer         ::      n5Jpoint,n6Jpoint
            
            integer         ::      ii,gg,kk,jj
            logical         ::      ok
            integer         ::      nTriang,nGrains
            real(kind=real64),dimension(:,:,:),allocatable  ::      triang
            type(GridCell)                                  ::      grid
            real(kind=real64),dimension(:,:),allocatable    ::      phi
            real(kind=real64),dimension(3)                  ::      xx
            logical,dimension(:),allocatable         ::      hasG
            integer                     ::      nTJ!,nQJ,n5J,n6J
            integer,dimension(1000)      ::      nC3,nC4,nC5,nC6
            
            
            integer                             ::      v1,v2,v3 !, e1,e2,e3 , g1,g2,g3
            !real(kind=real64),dimension(3,2,100)  ::      xLine
            integer                             ::      nc,nTJL
            real(kind=real64),dimension(8)      ::      aa,bb,cc
            real(kind=real64)                   ::      totalA
            real(kind=real64),dimension(3,2,1000)            ::      TJL
              
            !real(kind=real64)                   ::      vv
        !---    first find all the grains and all the grain boundaries detected here.
           
            nGrains = size(Garea)
            Garea = 0.0d0
            Gvol = 0.0d0
            TJlen = 0.0d0
            GBarea = 0.0d0 
            nQJpoint = 0
            allocate(phi(8,nGrains)) ; phi = -1.0d0
            allocate(hasG(nGrains))  ; hasG = .false.
            allocate(triang(3,3,1000))

            do ii = 1,8
                do gg = 1,nGrains
                    phi(ii,gg) = getPhi( node(ii),gg )
                    hasG(gg) = hasG(gg) .or. (phi(ii,gg)>-1)
                end do
            end do
            
            if (present(nTriOut))  nTriOut = 0
            
            
            do gg = 1,nGrains
                if (.not. hasG(gg)) cycle
                 
                grid = GridCell_ctor(phi(:,gg)) ; nTriang = 0
                call Polygonize( grid,0.0d0, triang, nTriang,ok ) 
                if (nTriang==0) cycle           !   no surface here
               
            !---    convert cell space position of triangles into real space positions.
                do kk = 1,nTriang
                    do ii = 1,3               ! vertex j
                        xx = triang(1:3,ii,kk) + i0 + 0.5d0
                        triang(1:3,ii,kk) = a(1:3,1)*xx(1) + a(1:3,2)*xx(2) + a(1:3,3)*xx(3)                         
                    end do
                end do
                
                call areaAndVolume( triang(:,:,1:nTriang),Garea(gg),Gvol(gg) )                
                
                if (present(nTriOut)) nTriOut = nTriOut + nTriang
                
            end do
            
            
            
            if (present(nTriOut)) then
                allocate( TriOut(3,3,nTriOut) )                
                allocate( TriGrain(nTriOut) )    
                nTriOut = 0            
                do gg = 1,nGrains
                    if (.not. hasG(gg)) cycle                     
                    grid = GridCell_ctor(phi(:,gg)) ; nTriang = 0
                    call Polygonize( grid,0.0d0, triang, nTriang,ok ) 
                !---    convert cell space position of triangles into real space positions.
                    do kk = 1,nTriang
                        TriGrain(nTriOut+kk) = gg
                        do ii = 1,3               ! vertex i, triangle k
                            xx = triang(1:3,ii,kk) + i0 + 0.5d0
                            TriOut(1:3,ii,nTriOut+kk) = a(1:3,1)*xx(1) + a(1:3,2)*xx(2) + a(1:3,3)*xx(3)                         
                        end do
                    end do        
                    nTriOut = nTriOut + nTriang            
                end do
            end if    
            deallocate(phi)       
            deallocate(hasG)
            kk = 0
            do ii = 1,8
                kk = kk + getNe(node(ii)%g)                
            end do      
            totalA = sum(Garea)/max(1,kk)  
            do ii = 1,8
                do jj = 1,getNe(node(ii)%g)                
                    call getEdge( node(ii)%g ,jj , v1,v2 )
                    kk = whichEdge(fullGraph,v1,v2 )                        
                    GBarea(kk) = GBarea(kk) + totalA
                end do
            end do
                
            
            
        !---    find six junction points             
            call findMatching6J( n6Jpoint,nC6 )
            
            
            
        !---    find five junction points             
            call findMatching5J( n5Jpoint,nC5 )
            
            
            
            
        !---    find quad junction points             
            call findMatchingQJ( nQJpoint,nC4 )
            
            
            
        !---    find triple junctions
            TJlen = 0
            call findMatchingTJ( nTJ,nC3 )
            
            
           ! GrainNode_dbg = GrainNode_dbg .or. (n6Jpoint>0)
             
             
            if (GrainNode_dbg) then 
                print *,"debug node at ",i0
                write (*,fmt='(a)',advance="no") "node 1 " ; call report(node(1 ))
                write (*,fmt='(a)',advance="no") "node 2 " ; call report(node(2 ))
                write (*,fmt='(a)',advance="no") "node 3 " ; call report(node(3 ))
                write (*,fmt='(a)',advance="no") "node 4 " ; call report(node(4 ))
                write (*,fmt='(a)',advance="no") "node 5 " ; call report(node(5 ))
                write (*,fmt='(a)',advance="no") "node 6 " ; call report(node(6 ))
                write (*,fmt='(a)',advance="no") "node 7 " ; call report(node(7 ))
                write (*,fmt='(a)',advance="no") "node 8 " ; call report(node(8 ))
                  
            end if
            nTJline = 0
            do kk = 1,nTJ
                
                nc = nC3(kk)
            
            
                if (GrainNode_dbg) print *,"findMatchingTJ (#",nc,") at ",i0
                v1 = getVertex( clique3(nc),1 )
                v2 = getVertex( clique3(nc),2 )
                v3 = getVertex( clique3(nc),3 )
                
                
                
                aa = -2 ; bb = -2 ; cc = -2
                
                
                !---    send distance across grain boundaries
             !  do ii = 1,8
             !      
             !      do jj = 1,node(ii)%nGrainBoundaries
             !          call getEdge( node(ii)%g,jj,g1,g2 )
             !          if ( (g1-v1)*(g1-v1) + (g2-v2)*(g2-v2) == 0) then
             !              aa(ii) = node(ii)%phi(jj)
             !          else if ( (g1-v2)*(g1-v2) + (g2-v1)*(g2-v1) == 0) then
             !              aa(ii) = -node(ii)%phi(jj)
             !          else if ( (g1-v2)*(g1-v2) + (g2-v3)*(g2-v3) == 0) then
             !              bb(ii) = node(ii)%phi(jj)
             !          else if ( (g1-v3)*(g1-v3) + (g2-v2)*(g2-v2) == 0) then
             !              bb(ii) = -node(ii)%phi(jj)
             !          else if ( (g1-v3)*(g1-v3) + (g2-v1)*(g2-v1) == 0) then
             !              cc(ii) = node(ii)%phi(jj)
             !          else if ( (g1-v1)*(g1-v1) + (g2-v3)*(g2-v3) == 0) then
             !              cc(ii) = -node(ii)%phi(jj)
             !          end if
             !      end do
             !  end do
               
                 !---     send value of grain
                
                  aa(1:8) = getPhi(node(1:8),v1) 
                  bb(1:8) = getPhi(node(1:8),v2) 
                  cc(1:8) = getPhi(node(1:8),v3) 
                 
                 
                if (GrainNode_dbg) then
                    if (any( (/ minval(aa),minval(bb),minval(cc) /) < -1 )) then
                        print *,"TJ vertices ",v1,v2,v3
                        print *,"error minval a,b,c"
                        print *,aa
                        print *,bb
                        print *,cc
                        stop
                    end if
                end if                    
               
                !
                if (GrainNode_dbg) then
                    !print *,"vertex data ",i0
                    print *,"TJ ",kk,"/",nTJ," (",v1,",",v2,",",v3,")"
                    write(*,fmt='(a8,3i10)') " ",v1,v2,v3
                    write(*,fmt='(a8,3f10.5)') "node 1 ",aa(1),bb(1),cc(1)
                    write(*,fmt='(a8,3f10.5)') "node 2 ",aa(2),bb(2),cc(2)
                    write(*,fmt='(a8,3f10.5)') "node 3 ",aa(3),bb(3),cc(3)
                    write(*,fmt='(a8,3f10.5)') "node 4 ",aa(4),bb(4),cc(4)
                    write(*,fmt='(a8,3f10.5)') "node 5 ",aa(5),bb(5),cc(5)
                    write(*,fmt='(a8,3f10.5)') "node 6 ",aa(6),bb(6),cc(6)
                    write(*,fmt='(a8,3f10.5)') "node 7 ",aa(7),bb(7),cc(7)
                    write(*,fmt='(a8,3f10.5)') "node 8 ",aa(8),bb(8),cc(8)
                    
                    print *,"4,5,6 junction ",nQJpoint ,n5Jpoint ,n6Jpoint
                    
                end if
                
                 dbg_SolveTripleJunction = GrainNode_dbg!
                call SolveTripleJunctionLine( aa,bb,cc , nTJL,TJL )
                 
                    
                if (GrainNode_dbg .and. (nTJline==0)) print *,"no TJ lines found"  
                
                do ii = 1,nTJL
                  !  if (GrainNode_dbg) write(*,fmt='(a,i4,a,6f10.4,a,f10.4)') "TJ line ",nTJline+ii," x ",TJL(:,1,ii),TJL(:,2,ii)," len ",norm2( TJL(:,2,ii)-TJL(:,1,ii) )
                    xx = TJL(:,1,ii) + i0 + 0.5d0
                    TJline(:,1,nTJLine+ii) = a(1:3,1)*xx(1) + a(1:3,2)*xx(2) + a(1:3,3)*xx(3)  
                    xx = TJL(:,2,ii) + i0 + 0.5d0
                    TJline(:,2,nTJLine+ii) = a(1:3,1)*xx(1) + a(1:3,2)*xx(2) + a(1:3,3)*xx(3)  
                    
                    xx = TJline(:,2,nTJLine+ii) - TJline(:,1,nTJLine+ii) 
                    TJlen(nc) = TJlen(nc) + norm2( xx )
                    
                     if (GrainNode_dbg) write(*,fmt='(2(a,i4),a,6f10.4,a,f10.4)') "TJ line (",kk,")",nTJline+ii," x ",TJL(:,1,ii),TJL(:,2,ii)," len(real) ",norm2( xx )
                    
                  !  do jj = 1,10
                  !      write(*,fmt='(a,3f16.8,i8)') "Tj ",xLine(:,1,ii) *(10-jj)*0.1d0 + xLine(:,2,ii) *jj*0.1d0,nc
                  !  end do
                  
                  TJlineId(nTJline+ii) = nc
                end do                 
                nTJLine = nTJLine + nTJL   
                
            end do
            !if (GrainNode_dbg .and. (n6Jpoint>0)) stop
            
            return
            
        contains
    !---^^^^^^^^
    
            subroutine findMatchingTJ( nTJ,nC3 )
        !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ 
        !*      try to find a matching triple point junction in nodes 
        !*      return nC3 = 0 if one can't be found             
                integer,intent(out)                         ::      nTJ
                integer,dimension(:),intent(out)            ::      nC3
                
                integer         ::      k1,k2, nc   !,k3,k4,k5,k6,k7,k8
                
                !logical         ::      ok
                
                nTJ = 0
                nC3 = 0
                
                
                if (product( node(:)%nTripleJunctions ) == 0) return                                 !    non-zero only if all nodes have a triple junction
                
                do k1 = 1,8
                    do k2 = 1,node(k1)%nTripleJunctions
                        nc = node(k1)%tripleJunction(k2)            !   this is the one I'm trying to match
                        if (any(NC3(1:nTJ)==nc)) cycle
                        nTJ = nTJ+1
                        NC3(nTJ) = nc
                    end do
                end do
                
                return
            end subroutine findMatchingTJ
    
!             subroutine findMatchingTJ( nTJ,nC3 )
!         !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ 
!         !*      try to find a matching triple point junction in nodes 
!         !*      return nC3 = 0 if one can't be found             
!               integer,intent(out)                         ::      nTJ
!                 integer,dimension(:),intent(out)            ::      nC3
!                 
!                 integer         ::      k1,k2,k3,k4,k5,k6,k7,k8, nc
!                 
!                 logical         ::      ok
!                 
!                 nTJ = 0
!                 nC3 = 0
!                 
!                 
!                 if (product( node(:)%nTripleJunctions ) == 0) return                                 !    non-zero only if all nodes have a triple junction
!                 
!                 
!                 do k1 = 1,node(1)%nTripleJunctions
!                     nc = node(1)%tripleJunction(k1)            !   this is the one I'm trying to match
!                     ok = .false.
!                     do k2 = 1,node(2)%nTripleJunctions
!                         if ( nc == node(2)%tripleJunction(k2) ) then
!                             do k3 = 1,node(3)%nTripleJunctions
!                                 if ( nc == node(3)%tripleJunction(k3) ) then
!                                     do k4 = 1,node(4)%nTripleJunctions
!                                         if ( nc == node(4)%tripleJunction(k4) ) then
!                                             do k5 = 1,node(5)%nTripleJunctions
!                                                 if ( nc == node(5)%tripleJunction(k5) ) then
!                                                     do k6 = 1,node(6)%nTripleJunctions
!                                                         if ( nc == node(6)%tripleJunction(k6) ) then
!                                                             do k7 = 1,node(7)%nTripleJunctions
!                                                                 if ( nc == node(7)%tripleJunction(k7) ) then
!                                                                     do k8 = 1,node(8)%nTripleJunctions
!                                                                         if ( nc == node(8)%tripleJunction(k8) ) then
!                                                                             nTJ = nTJ + 1
!                                                                             nC3(nTJ) = nc
!                                                                             ok = .true.
!                                                                         end if
!                                                                         if (ok) exit
!                                                                     end do
!                                                                 end if
!                                                                 if (ok) exit
!                                                             end do
!                                                         end if
!                                                         if (ok) exit
!                                                     end do
!                                                 end if
!                                                 if (ok) exit                                                                    
!                                             end do
!                                         end if
!                                         if (ok) exit
!                                     end do
!                                 end if
!                                 if (ok) exit
!                             end do
!                         end if
!                         if (ok) exit
!                     end do
!                 end do
!                 return
!             end subroutine findMatchingTJ
!     
            subroutine findMatchingQJ( nQJ,nC4 )
        !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ 
        !*      try to find a matching quad point junction in nodes 
        !*      return nC4 = 0 if one can't be found             
                integer,intent(out)                         ::      nQJ
                integer,dimension(:),intent(out)            ::      nC4
                
                integer         ::      k1,k2,k3,k4,k5,k6,k7,k8, nc
                
                logical         ::      ok
                
                nQJ = 0
                nC4 = 0
                
                
                if (product( node(:)%nQuadJunctions ) == 0) return                                 !    non-zero only if all nodes have a triple junction
                
                
                do k1 = 1,node(1)%nQuadJunctions
                    nc = node(1)%quadJunction(k1)            !   this is the one I'm trying to match
                    ok = .false.
                    do k2 = 1,node(2)%nQuadJunctions
                        if ( nc == node(2)%quadJunction(k2) ) then
                            do k3 = 1,node(3)%nQuadJunctions
                                if ( nc == node(3)%quadJunction(k3) ) then
                                    do k4 = 1,node(4)%nQuadJunctions
                                        if ( nc == node(4)%quadJunction(k4) ) then
                                            do k5 = 1,node(5)%nQuadJunctions
                                                if ( nc == node(5)%quadJunction(k5) ) then
                                                    do k6 = 1,node(6)%nQuadJunctions
                                                        if ( nc == node(6)%quadJunction(k6) ) then
                                                            do k7 = 1,node(7)%nQuadJunctions
                                                                if ( nc == node(7)%quadJunction(k7) ) then
                                                                    do k8 = 1,node(8)%nQuadJunctions
                                                                        if ( nc == node(8)%quadJunction(k8) ) then
                                                                            nQJ = nQJ + 1
                                                                            nC4(nQJ) = nc
                                                                            ok = .true.
                                                                        end if
                                                                        if (ok) exit
                                                                    end do
                                                                end if
                                                                if (ok) exit
                                                            end do
                                                        end if
                                                        if (ok) exit
                                                    end do
                                                end if
                                                if (ok) exit                                                                    
                                            end do
                                        end if
                                        if (ok) exit
                                    end do
                                end if
                                if (ok) exit
                            end do
                        end if
                        if (ok) exit
                    end do
                end do
                return
            end subroutine findMatchingQJ
            
            subroutine findMatching5J( n5J,NC5 )
        !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ 
        !*      try to find a matching five point junction in nodes 
        !*      return NC5 = 0 if one can't be found             
                 integer,intent(out)                        ::      n5J
                integer,dimension(:),intent(out)            ::      nC5

                
                integer         ::      k1,k2,k3,k4,k5,k6,k7,k8, nc
          
                logical         ::      ok
                
                n5J = 0
                NC5 = 0
                
                
                if (product( node(:)%n5Junctions ) == 0) return                                 !    non-zero only if all nodes have a triple junction
                
                
                do k1 = 1,node(1)%n5Junctions
                    nc = node(1)%FiveJunction(k1)            !   this is the one I'm trying to match
                    ok = .false.
                    do k2 = 1,node(2)%n5Junctions
                        if ( nc == node(2)%FiveJunction(k2) ) then
                            do k3 = 1,node(3)%n5Junctions
                                if ( nc == node(3)%FiveJunction(k3) ) then
                                    do k4 = 1,node(4)%n5Junctions
                                        if ( nc == node(4)%FiveJunction(k4) ) then
                                            do k5 = 1,node(5)%n5Junctions
                                                if ( nc == node(5)%FiveJunction(k5) ) then
                                                    do k6 = 1,node(6)%n5Junctions
                                                        if ( nc == node(6)%FiveJunction(k6) ) then
                                                            do k7 = 1,node(7)%n5Junctions
                                                                if ( nc == node(7)%FiveJunction(k7) ) then
                                                                    do k8 = 1,node(8)%n5Junctions
                                                                        if ( nc == node(8)%FiveJunction(k8) ) then
                                                                            n5J = n5J + 1
                                                                            NC5(n5J) = nc
                                                                            ok = .true.
                                                                        end if
                                                                        if (ok) exit
                                                                    end do
                                                                end if
                                                                if (ok) exit
                                                            end do
                                                        end if
                                                        if (ok) exit
                                                    end do
                                                end if
                                                if (ok) exit                                                                    
                                            end do
                                        end if
                                        if (ok) exit
                                    end do
                                end if
                                if (ok) exit
                            end do
                        end if
                        if (ok) exit
                    end do
                end do
                return
            end subroutine findMatching5J
            
            subroutine findMatching6J( n6J,NC6 )
        !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ 
        !*      try to find a matching six point junction in nodes 
        !*      return NC6 = 0 if one can't be found             
                integer,intent(out)                         ::      n6J
                integer,dimension(:),intent(out)            ::      nC6

                
                integer         ::      k1,k2,k3,k4,k5,k6,k7,k8, nc
                
                logical         ::      ok
                
                n6J = 0
                NC6 = 0
                
                
                if (product( node(:)%n6Junctions ) == 0) return                                 !    non-zero only if all nodes have a triple junction
                
                
                do k1 = 1,node(1)%n6Junctions
                    nc = node(1)%sixJunction(k1)            !   this is the one I'm trying to match
                    ok = .false.
                    do k2 = 1,node(2)%n6Junctions
                        if ( nc == node(2)%sixJunction(k2) ) then
                            do k3 = 1,node(3)%n6Junctions
                                if ( nc == node(3)%sixJunction(k3) ) then
                                    do k4 = 1,node(4)%n6Junctions
                                        if ( nc == node(4)%sixJunction(k4) ) then
                                            do k5 = 1,node(5)%n6Junctions
                                                if ( nc == node(5)%sixJunction(k5) ) then
                                                    do k6 = 1,node(6)%n6Junctions
                                                        if ( nc == node(6)%sixJunction(k6) ) then
                                                            do k7 = 1,node(7)%n6Junctions
                                                                if ( nc == node(7)%sixJunction(k7) ) then
                                                                    do k8 = 1,node(8)%n6Junctions
                                                                        if ( nc == node(8)%sixJunction(k8) ) then
                                                                            n6J = n6J + 1
                                                                            NC6(n6J) = nc
                                                                            ok = .true.
                                                                        end if
                                                                        if (ok) exit
                                                                    end do
                                                                end if
                                                                if (ok) exit
                                                            end do
                                                        end if
                                                        if (ok) exit
                                                    end do
                                                end if
                                                if (ok) exit                                                                    
                                            end do
                                        end if
                                        if (ok) exit
                                    end do
                                end if
                                if (ok) exit
                            end do
                        end if
                        if (ok) exit
                    end do
                end do
                return
            end subroutine findMatching6J
            
            
        end subroutine cubeProperties
        
    !---
        
            
        elemental function getnGrainBoundaries(this) result(n)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^  
            type(GrainNode),intent(in)      ::      this
            integer                         ::      n
            n = this%nGrainBoundaries
            return
        end function getnGrainBoundaries
        
        elemental function getnTripleJunctions(this) result(n)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^  
            type(GrainNode),intent(in)      ::      this
            integer                         ::      n
            n = this%nTripleJunctions
            return
        end function getnTripleJunctions
    
        elemental function getnQuadJunctions(this) result(n)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^  
            type(GrainNode),intent(in)      ::      this
            integer                         ::      n
            n = this%nQuadJunctions
            return
        end function getnQuadJunctions
    
        elemental function getn5Junctions(this) result(n)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^  
            type(GrainNode),intent(in)      ::      this
            integer                         ::      n
            n = this%n5Junctions
            return
        end function getn5Junctions
    
        elemental function getn6Junctions(this) result(n)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^  
            type(GrainNode),intent(in)      ::      this
            integer                         ::      n
            n = this%n6Junctions
            return
        end function getn6Junctions
    
    
    
    
    
    !---
              
    
    
        subroutine areaAndVolume( tris,a,v )
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
            real(kind=real64),dimension(:,:,:),intent(in)           ::      tris
            real(kind=real64),intent(out)                           ::      a,v
            
            real(kind=real64),dimension(3)      ::      x12,x13,vecta
            real(kind=real64)                   ::      dd
            integer     ::      kk,nTri
            
            nTri = size(tris,dim=3)
            a = 0.0d0 ; v = 0.0d0
            do kk = 1,nTri                
                                    
            !---    find offset from vertex 1
                x12(1:3) = tris(1:3,2,kk) - tris(1:3,1,kk)
                x13(1:3) = tris(1:3,3,kk) - tris(1:3,1,kk)
                
            !---    find cross product
                vecta(1) = x12(2)*x13(3) - x12(3)*x13(2)
                vecta(2) = x12(3)*x13(1) - x12(1)*x13(3)
                vecta(3) = x12(1)*x13(2) - x12(2)*x13(1)
                
            !---    find area
                dd = norm2(vecta) / 2
                a = a + dd
                
                
                
            !---    find volume                        
                dd = ( tris(1,1,kk)*vecta(1) + tris(2,1,kk)*vecta(2) + tris(3,1,kk)*vecta(3) )/6       
                v = v + dd
                
            end do

            return
        end subroutine areaAndVolume
    
    !---
    
        subroutine permutations( n,k, list )
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      find all permutations of k from n , (k,n>3)
            integer,intent(in)                                ::      n,k
            integer,dimension(:,:),pointer,intent(out)        ::      list
            
            !print *,"permutations ",n,k
            nullify(list)
            select case( n )
                case(2)
                    select case(k)
                        case(2) ; list => C22
                    end select
                case(3)
                    select case(k)
                        case(2) ; list => C23
                        case(3) ; list => C33
                    end select
                case(4)
                    select case(k)
                        case(2) ; list => C24
                        case(3) ; list => C34
                        case(4) ; list => C44
                    end select
                case(5)
                    select case(k)
                        case(2) ; list => C25
                        case(3) ; list => C35
                        case(4) ; list => C45
                        case(5) ; list => C55
                    end select                
                case(6)
                    select case(k)
                        case(2) ; list => C26
                        case(3) ; list => C36
                        case(4) ; list => C46
                        case(5) ; list => C56
                        case(6) ; list => C66
                    end select     
                case(7)
                    select case(k)
                        case(2) ; list => C27
                        case(3) ; list => C37
                        case(4) ; list => C47
                        case(5) ; list => C57
                        case(6) ; list => C67
                        case(7) ; list => C77
                    end select     
                case(8)
                    select case(k)
                        case(2) ; list => C28
                        case(3) ; list => C38
                        case(4) ; list => C48
                        case(5) ; list => C58
                        case(6) ; list => C68
                        case(7) ; list => C78
                        case(8) ; list => C88
                    end select                
                case(9)
                    select case(k)
                        case(2) ; list => C29
                        case(3) ; list => C39
                        case(4) ; list => C49      
                        case(5) ; list => C59                        
                        case(6) ; list => C69                        
                    end select                   
                case(10)
                    select case(k)
                        case(2) ; list => C210
                        case(3) ; list => C310
                        case(4) ; list => C410     
                        case(5) ; list => C510     
                        case(6) ; list => C610     
                    end select                                                         
            end select

            if (.not. associated(list)) then
                print *,"GrainNodes::permutations() warning - have not coded nCk for n=",n," k=",k
                print *,"                                     this probably violates the assumptions about number of grains at one node"
                print *,"                                     consider reducing node spacing and check input files"
                stop
            end if
                         
            return 
        end subroutine permutations                 
            
    end module GrainNodes
    
    
!   gfortran -ffree-line-length-256 ${MYF90LIB}/Lib_Quicksort.f90  Lib_Graphs.f90 Lib_FitTanhFunction.f90 ${MYF90LIB}/NBAX_StringTokenizers.f90 ${MYF90LIB}/Lib_MarchingCubes.f90 ${MYF90LIB}/Lib_ReadVTK.f90 SolveTripleJunction.f90 GrainNodes.f90  -o testGrainNodes.exe -llapack  -fbacktrace -fbounds-check -g -Wall -ffpe-trap=invalid,zero,overflow -Wno-unused-function -DDEBUG -Wno-unused-dummy-argument 

!    
!    program testGrainNodes
!!---^^^^^^^^^^^^^^^^^^^^^^^^
!        use iso_fortran_env
!        use GrainNodes
!        use Lib_ReadVTK
!        implicit none
!        
!        real(kind=real64)                   ::      grainDev = 0.01 
!        real(kind=real64)                   ::      rc 
!        integer                             ::      N = 10000     !   number of atoms
!        integer                             ::      Nx = 8
!        
!         
!        
!        integer,dimension(3),parameter      ::      GRAIN = (/ 5,10,15 /)
!        
!        
!        integer             ::      ii,nn
!        real(kind=real64)   ::      d1,d2,d3
!        
!
!        integer                         ::        ix,iy,iz  , jx,jy,jz
!        real(kind=real64),dimension(3)  ::      dx
!        character(len=5)                ::      aaaa
!        character(len=256)              ::      dummy
!        
!        
!        
!        
!        
!        real(kind=real64),dimension(:,:),allocatable    ::      x           !   position of atoms
!        integer,dimension(:),allocatable                ::      g           !   grain number of atoms        
!        real(kind=real64),dimension(:,:),allocatable    ::      yy,xx          !   position of atoms
!        integer,dimension(:),allocatable                ::      gg          !   grain number of atoms
!        logical,dimension(:),allocatable                ::      mask
!        type(GrainNode),dimension(:,:,:),allocatable    ::      this
!        real(kind=real64),dimension(:,:,:,:),allocatable::      phi
!        
!        type(GrainNode),dimension(8)                    ::      node
!        real(kind=real64),dimension(3,3)                ::      aa
!                
!        !integer,dimension(:),allocatable                  ::      grainsFound
!        real(kind=real64),dimension(:),allocatable        ::      area,darea,volume,dvolume,TJlen,dTJlen
!        
!        
!        call init_random_seed(12345)
!        
!        call get_command_argument(1,dummy)
!        read(dummy,fmt=*) N
!        call get_command_argument(2,dummy)
!        read(dummy,fmt=*) Nx
!        
!        allocate(x(3,N))
!        allocate(xx(3,N))
!        allocate(yy(3,N))
!        allocate(g(N))
!        allocate(gg(N))
!        allocate(mask(N))
!        allocate(this(-1:Nx,-1:Nx,-1:Nx))
!        allocate(phi(0:Nx-1,0:Nx-1,0:Nx-1,size(GRAIN)))
!        
!        do iz = -1,Nx
!            do iy = -1,Nx
!                do ix = -1,Nx
!                    this(ix,iy,iz) = GrainNode_ctor()
!                end do
!            end do
!        end do
!        
!    !---    place atoms in random positions 0:1
!        call random_number(x) 
!       
!       
!    !---    determine grain number by considering projection on three lobes
!    !       [001],[ (sqrt 3)/2,0,-1/2 ], [ -(sqrt 3)/2,0,-1/2 ]
!        do ii = 1,N
!            dx = x(1:3,ii) - (/0.5d0,0.5d0,0.5d0/)
!            d1 = dx(3) 
!            d2 = ( dx(1)*sqrt(3.0) - dx(3) )/2
!            d3 = (-dx(1)*sqrt(3.0) - dx(3) )/2
!            
!            if (d1 > max(d2,d3)+grainDev) then
!                g(ii) = GRAIN(1)
!            else if (d2 > max(d1,d3)+grainDev) then
!                g(ii) = GRAIN(2)
!            else if (d3 > max(d1,d2)+grainDev) then
!                g(ii) = GRAIN(3)
!            else
!                g(ii) = 0
!            end if
!            
!        end do
!        open(unit=400,file="test.xyz",action="write")
!            write(400,fmt='(i8)') N
!            write(400,fmt='(a)') "Lattice=""1.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 1.0"" Properties=species:S:1:pos:R:3:grain:I:1"
!
!            do ii = 1,N
!                write(400,fmt='(a,3f16.8,i6)') "Du ",x(:,ii),g(ii)
!            end do
!        close(unit=400)
!        
!        !rc = 1.0d0/Nx
!        rc = sqrt(3.0d0)/Nx
!        
!        
!        do iz = 0,Nx-1
!            do iy = 0,Nx-1
!                do ix = 0,Nx-1
!                                
!                    dx = ((/ ix,iy,iz /)+0.5d0)/Nx
!                    call neighbourList( x,dx,rc,mask,yy )
!                    nn = 0
!                    do ii = 1,N
!                        if (mask(ii)) then
!                            nn = nn + 1
!                            xx(:,nn) = yy(:,ii)
!                            gg(nn) = g(ii)
!                        end if
!                    end do
!                    
!           !         print *,"grain ",ix,iy,iz," at ",dx," nn = ",nn
!                    this(ix,iy,iz) = GrainNode_ctor(gg(1:nn),xx(:,1:nn),rc)
!                    if (count(gg(1:nn)/=0)<5) print *,"testGrainNodes.exe warning - node has very few (",count(gg(1:nn)/=0),") grain id atoms"
!           !         call report(this(ix,iy,iz))
!        
!                end do
!            end do
!        end do
!        
!               
!        do iz = 0,Nx-1
!            do iy = 0,Nx-1
!                do ix = 0,Nx-1
!                     do ii = 1,size(GRAIN) !   grain index        
!                        phi(ix,iy,iz,ii) = getPhi( this(ix,iy,iz),GRAIN(ii) )
!                    end do
!                    dx = ((/ ix,iy,iz /)+0.5d0)/Nx
!                    !print *,"grain ",ix,iy,iz," at ",dx," phi = ",phi(ix,iy,iz,:)
!                end do
!            end do
!        end do    
!        do ii = 1,size(GRAIN)
!            write(aaaa,fmt='(i5)') ii ; aaaa = adjustl(aaaa) ; aaaa = repeat("0",5-len_trim(aaaa))//trim(aaaa)
!            call writeChgcar( "test."//aaaa//".chgcar",Nx,Nx,Nx,(/0.0d0,0.0d0,0.0d0/),(/1.0d0,1.0d0,1.0d0/), phi(:,:,:,ii) )        
!        end do    
!        
!       ! print *,"grain 1"        
!       ! do iz = 0,Nx-1
!       !     write(*,fmt='(1000f9.4)')  phi(:,Nx/2,iz,1) 
!       ! end do
!       ! print *,"grain 2"        
!       ! do iz = 0,Nx-1
!       !     write(*,fmt='(1000f9.4)')  phi(:,Nx/2,iz,2) 
!       ! end do
!       ! print *,"grain 3"        
!       ! do iz = 0,Nx-1
!       !     write(*,fmt='(1000f9.4)')  phi(:,Nx/2,iz,3) 
!       ! end do
!       ! print *,"sum"        
!       ! do iz = 0,Nx-1
!       !     write(*,fmt='(1000f9.4)')  phi(:,Nx/2,iz,1)+phi(:,Nx/2,iz,2)+phi(:,Nx/2,iz,3) 
!       ! end do
!        !call delete(this)
!        
!        
!        aa = 0.0d0
!        aa(1,1) = 1.0d0/Nx
!        aa(2,2) = 1.0d0/Nx
!        aa(3,3) = 1.0d0/Nx
!        
!        nn = maxval(GRAIN)
!        allocate( area(nn) ) ; area = 0
!        allocate( darea(nn) )
!        allocate( volume(nn) ) ; volume = 0
!        allocate( dvolume(nn) )
!        nn = nClique3
!        allocate( TJlen(nn) ) ; TJlen = 0
!        allocate( dTJlen(nn) )
!        
!        do iz = -1,Nx-1
!            do iy = -1,Nx-1
!                do ix = -1,Nx-1
!                
!                !---    find nodes for volume calculation - includes empty cells outside box
!                    node(1) = this(ix  ,iy  ,iz  )
!                    node(2) = this(ix+1,iy  ,iz  )
!                    node(3) = this(ix+1,iy+1,iz  )
!                    node(4) = this(ix  ,iy+1,iz  )
!                    node(5) = this(ix  ,iy  ,iz+1)
!                    node(6) = this(ix+1,iy  ,iz+1)
!                    node(7) = this(ix+1,iy+1,iz+1)
!                    node(8) = this(ix  ,iy+1,iz+1)
!                                                  
!                    call cubeProperties( node,(/ix,iy,iz/),aa ,darea,dvolume,dTJlen )                    
!                    volume = volume + dvolume
!                    !area = area + darea 
!                
!                    if ( minval( (/ix,iy,iz,Nx-2-ix,Nx-2-iy,Nx-2-iz/) ) >= 0 ) then
!                        area = area + darea
!                        TJlen = TJlen + dTJlen
!                    !    print *,ix,iy,iz,darea(5:15:5)
!                    else if ( minval( (/ix,iy,iz,Nx-1-ix,Nx-1-iy,Nx-1-iz/) ) == 0 ) then
!                        node(1) = this(    ix         ,    iy         ,    iz         )
!                        node(2) = this(mod(ix+1+Nx,Nx),    iy         ,    iz         )
!                        node(3) = this(mod(ix+1+Nx,Nx),mod(iy+1+Nx,Nx),    iz         )
!                        node(4) = this(    ix         ,mod(iy+1+Nx,Nx),    iz         )
!                        node(5) = this(    ix         ,    iy         ,mod(iz+1+Nx,Nx))
!                        node(6) = this(mod(ix+1+Nx,Nx),    iy         ,mod(iz+1+Nx,Nx))
!                        node(7) = this(mod(ix+1+Nx,Nx),mod(iy+1+Nx,Nx),mod(iz+1+Nx,Nx))
!                        node(8) = this(    ix         ,mod(iy+1+Nx,Nx),mod(iz+1+Nx,Nx))
!                        call cubeProperties( node,(/ix,iy,iz/),aa ,darea,dvolume,dTJlen )               
!                        area = area + darea
!                        TJlen = TJlen + dTJlen   
!                    !    print *,ix,iy,iz,darea(5:15:5)                 
!                    end if
!              
!                end do
!            end do
!        end do    
!        
!        do ii = 1,size(area)
!            nn = count(g==ii)
!            if (nn>0) print *,"grain ",ii," count ",nn," area ",area(ii)," volume ",volume(ii)
!        end do
!        
!        do ii = 1,nClique3           
!            print *,"Triple Junction ",ii," length ",TJlen(ii)
!        end do
!        
!        print *,""
!        print *,"done"
!        print *,""
!        
!    contains
!!---^^^^^^^^
!
!        subroutine neighbourList( x,dx,rc,mask,y )
!    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
!            real(kind=real64),dimension(:,:),intent(in)     ::  x
!            real(kind=real64),dimension(3),intent(in)       ::  dx
!            real(kind=real64),intent(in)                    ::  rc
!            logical,dimension(:),intent(out)                ::  mask
!            real(kind=real64),dimension(:,:),intent(out)    ::  y
!            integer             ::      ii
!            real(kind=real64),dimension(3)  ::      delta
!            do ii = 1,size(x,dim=2)
!                delta = x(1:3,ii) - dx(1:3)
!                if (delta(1)>0.5) delta(1) = delta(1) - 1
!                if (delta(1)<-.5) delta(1) = delta(1) + 1
!                if (delta(2)>0.5) delta(2) = delta(2) - 1
!                if (delta(2)<-.5) delta(2) = delta(2) + 1
!                if (delta(3)>0.5) delta(3) = delta(3) - 1
!                if (delta(3)<-.5) delta(3) = delta(3) + 1
!                !mask(ii) = all(abs(delta)<=rc)
!                mask(ii) = norm2(delta)<=rc
!                y(1:3,ii) = delta
!            end do
!            return
!        end subroutine neighbourList
!        
!        
!
!        subroutine init_random_seed(seed_in)
!    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
!            integer,intent(in),optional         ::      seed_in
!            integer                             ::      ii, nn, clock
!            integer, dimension(:), allocatable  ::      seed
!
!            call random_seed(size = nn)
!            allocate(seed(nn))
!
!            if (present(seed_in)) then
!                seed = seed_in + 37 * (/ (ii - 1, ii = 1, nn) /)            
!            else
!                call system_clock(count=clock)
!                seed = clock + 37 * (/ (ii - 1, ii = 1, nn) /)
!            end if
!            
!            call random_seed(put = seed)
!
!            deallocate(seed)
!            return
!        end subroutine init_random_seed
!        
!        
!    end program testGrainNodes
!    
!        
!   
       
        
    