

    
    program testLib_Graphs
!---^^^^^^^^^^^^^^^^^^^^^^^
!*      simple program to test functioning of Lib_DeformationGradients
!*      
!*          successful result        
!*
!*      Graph [nV,nE =     5,    5]
!*                1  2  4  5 12
!*             1     .  *  *  .
!*             2  .     .  *  *
!*             4  *  .     *  .
!*             5  *  *  *     .
!*            12  .  *  .  .
!*       hasVertex(g1,5)   T
!*       hasVertex(g1,15)  F
!*       hasEdge(g1,1,5)   T
!*       hasEdge(g1,1,12)  F
!*      Graph [nV,nE =     3,    3]
!*                1  4  5
!*             1     *  *
!*             4  *     *
!*             5  *  *
!*      Graph [nV,nE =     2,    1]
!*                2 12
!*             2     *
!*            12  *
!*       g1==g1  T
!*       g1==g2  F
!*       g1==g3  F
!*       g2==g2  T
!*       g2==g3  F
!*       g3==g3  T
!*       subgraph(g1,g1)  T
!*       subgraph(g1,g2)  T
!*       subgraph(g1,g3)  T
!*       subgraph(g2,g2)  T
!*       subgraph(g2,g3)  F
!*       subgraph(g3,g3)  T
!*      
!*       done
!*      
!*      
!*



        use iso_fortran_env 
        use Lib_ColouredTerminal
        use Lib_Graphs
        implicit none
        
        type(Graph)            ::      g1,g2,g3
        
        
        logical,parameter                      ::      T = .true. , F = .false. , X = .false.
        integer,parameter                      ::      NV = 5
        integer,dimension(NV),parameter        ::      v1 = (/ 1,5,4,2,12 /)
        logical,dimension(NV,NV),parameter     ::      e1 = transpose( reshape( (/     X,T,T,F,F ,         &
                                                                                       X,X,T,T,F ,         &
                                                                                       X,X,X,F,F ,         &
                                                                                       X,X,X,X,T ,         &
                                                                                       X,X,X,X,X       /),(/NV,NV/) ) )
        
        integer,dimension(3),parameter         ::      v2 = (/ 4,1,5 /)
        logical,dimension(3,3),parameter       ::      e2 = transpose( reshape( (/     X,T,T ,     &
                                                                                       X,X,T ,      &
                                                                                       X,X,X     /),(/3,3/) ) )
        
        integer,dimension(2),parameter         ::      v3 = (/ 2,12 /)
        logical,dimension(2,2),parameter       ::      e3 = transpose( reshape( (/     X,T ,       &
                                                                                       X,X        /),(/2,2/) ) )
                                                                                       
                                                                                       


        character(len=256),dimension(16) ::      output
        character(len=*),dimension(16),parameter   ::      output0 = (/ " hasVertex(g1,5)   T         ",    &
                                                                        " hasVertex(g1,15)  F         ",    &
                                                                        " hasEdge(g1,1,5)   T         ",    &
                                                                        " hasEdge(g1,1,12)  F         ",    &
                                                                        " g1==g1  T                   ",    &
                                                                        " g1==g2  F                   ",    &
                                                                        " g1==g3  F                   ",    &
                                                                        " g2==g2  T                   ",    &
                                                                        " g2==g3  F                   ",    &
                                                                        " g3==g3  T                   ",    &
                                                                        " subgraph(g1,g1)  T          ",    &
                                                                        " subgraph(g1,g2)  T          ",    &
                                                                        " subgraph(g1,g3)  T          ",    &
                                                                        " subgraph(g2,g2)  T          ",    &
                                                                        " subgraph(g2,g3)  F          ",    &
                                                                        " subgraph(g3,g3)  T          "         /)
                                                                        
                                                                        
                                                                        
                                                                        
        logical                     ::  ok 
        integer                     ::      ii                                                             
                                                                   
                                                      
               
       
                                                                                              
                                                                                       
                                                                                       
                                                                                       
                                                                                       
                                                                                       
                                                                                       
                                                                                       
                                                                                       
                                                                                       
                                                                                       
                                                                                       
        
        g1 = Graph_ctor( v1,e1 )        
        call report(g1)

        write(output(1),fmt='(a,l4)') "hasVertex(g1,5)  ",hasVertex(g1,5)       
        write(output(2),fmt='(a,l4)') "hasVertex(g1,15) ",hasVertex(g1,15)       
        
        write(output(3),fmt='(a,l4)') "hasEdge(g1,1,5)  ",hasEdge(g1,1,5)   
        write(output(4),fmt='(a,l4)') "hasEdge(g1,1,12) ",hasEdge(g1,1,12)   
        
        g2 = Graph_ctor( v2,e2 )        
        call report(g2)
        
        
        g3 = Graph_ctor( v3,e3 )        
        call report(g3)
        write(output(5) ,fmt='(a,l4)')  "g1==g1 ",g1==g1
        write(output(6) ,fmt='(a,l4)')  "g1==g2 ",g1==g2
        write(output(7) ,fmt='(a,l4)')  "g1==g3 ",g1==g3
        write(output(8) ,fmt='(a,l4)')  "g2==g2 ",g2==g2
        write(output(9) ,fmt='(a,l4)')  "g2==g3 ",g2==g3
        write(output(10),fmt='(a,l4)')  "g3==g3 ",g3==g3
        
        write(output(11),fmt='(a,l4)') "subgraph(g1,g1) ",subgraph(g1,g1)
        write(output(12),fmt='(a,l4)') "subgraph(g1,g2) ",subgraph(g1,g2)
        write(output(13),fmt='(a,l4)') "subgraph(g1,g3) ",subgraph(g1,g3)
        write(output(14),fmt='(a,l4)') "subgraph(g2,g2) ",subgraph(g2,g2)
        write(output(15),fmt='(a,l4)') "subgraph(g2,g3) ",subgraph(g2,g3)
        write(output(16),fmt='(a,l4)') "subgraph(g3,g3) ",subgraph(g3,g3)
        
        
          
    !---    here is the simple test: does the output look like the stored output?
        ok = .true.
        do ii = 1,size(output)            
            if ( trim(cutSpaces(output(ii))) == trim(cutSpaces(output0(ii))) ) then
                write (*,fmt='(a)') trim(cutSpaces(output(ii)))
            else
                ok = .false.
                write (*,fmt='(a)') trim(cutSpaces(output0(ii)))//"    "//colour(RED,trim(cutSpaces(output(ii))))
            end if
        end do   
        
       
    !---    output the result "PASS" or "FAIL"    
        if (ok) then
            print *,colour(LIGHT_GREEN,"PASS")
        else
            print *,colour(RED,"FAIL")
        end if
        
        
        call delete(g1)
        call delete(g2)
        
        print *,""
        print *,"done"
        print *,""
        
    end program testLib_Graphs
   
        
    
         
         
     