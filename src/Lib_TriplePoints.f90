
    module Lib_TriplePoints
!---^^^^^^^^^^^^^^^^^^^^^^^
!*      given a distribution of grains
!*      it is important to be able to find the area of intersection, 
!*      plus triple point lines, and higher index points of coincidence
!*      This module returns the triple point lines.
!*

!       consider a cubic cell with 8 nodes, each of which is given a grain index    
!
!
!             ________________
!            7               8       grain = (/ a,b,c,d,e,f,g,h /)
!           /.              /|
!          / .             / |
!         5______________6   |     
!        |   .           |   |
!        |   .           |   |
!        |   .           |   |
!        |   3...........|...4
!        |  .            |  / 
!        | .             | /
!        1_______________2/
!         
!
!       A triple line bounding three grains must pass through one face, and end on a high point
!       or must pass through two faces.
!
!       So there are only seven points of interest in the cube - the centres of the faces and the centre of the cube.
!       A face is intesting if one of the four triangles made from its corners has three dissimilar grains.

        use iso_fortran_env
        implicit none
        private

    contains
!---^^^^^^^^

            
        subroutine analyseCube( ff )
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^                                                                                     
            integer,dimension(8),intent(in)     ::      ff                                                               
            integer,dimension(4,6)      ::      face                                                               
!            integer,dimension(3)        ::      tri                                                                
            integer,dimension(8)        ::      cc                                                                 
            integer             ::      jj,kk                                              
 !           logical             ::      centreInteresting                                                           
                
            integer,dimension(4,6),parameter        ::      FACEPT = reshape((/ 1,2,3,4     &                                  
                                                                               ,1,2,5,6     &
                                                                               ,1,3,5,8     &
                                                                               ,2,4,6,8     &
                                                                               ,3,4,7,8     &
                                                                               ,5,6,7,8 /),(/4,6/) )
                    
!             integer,dimension(3,8),parameter        ::      TRIPT =  reshape( (/ 4,6,7      &
!                                                                                 ,3,5,8      &
!                                                                                 ,2,5,8      &
!                                                                                 ,1,6,7      &
!                                                                                 ,2,3,8      &
!                                                                                 ,1,4,7      &
!                                                                                 ,1,4,6      &
!                                                                                 ,2,3,5 /),(/3,8/) ) 
            
                      
        !---    analyse the cube, look for anything of interest                                                    
            cc = ff                                                                                                
            call bubblesortAndTrim( cc )                                                                           
            if (cc(1) == 0) return          !   nothing to do, there is only 1 or 2 distinct numbers in ff.        
                                                                                                                   
                                                                                                                   
        !---    OK, I know there is something interesting here. But what, and where?
            do kk = 1,6
                do jj = 1,4
                    face(jj,kk) = ff( FACEPT(jj,kk))
                end do
                call bubblesortAndTrim( face(:,kk) )
            end do
            
             
!              
!         !---    interesting points appear where there are non-zero values.
!             if (centreInteresting) then
!                 !   all roads lead in and out of the centre.
!                 
!             else
!                 !   all roads lead from one face to another.
!                 
!             end if
            return
       end subroutine analyseCube                   
            

       pure subroutine bubblesortAndTrim( array )
   !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
   !*       sort the values of array into order, then delete duplicates
   !*       if there are fewer than 3 remaining numbers, delete all
   !*       eg (/1,2,4,3/) -> (/1,2,3,4/)
   !*          (/1,3,2,1/) -> (/1,1,2,3/) -> (/1,2,3,0/)
   !*          (/1,1,1,2/) -> (/1,1,1,2/) -> (/1,2,0,0/) -> (/0,0,0,0/)
   
            integer,dimension(:),intent(inout)      ::      array
            integer                         ::      kk,jj,aa,nn
            logical                         ::      done
            nn = size(array)
            do 
                done = .true.
                do kk = 1,nn-1                    
                    if (array(kk)>array(kk+1)) then
                        aa = array(kk) ; array(kk) = array(kk+1) ; array(kk+1) = aa ; done = .false.
                    end if 
                end do
                if (done) exit
            end do
            do kk = 1,nn-1
                if (array(kk) == array(kk+1)) then
                    do jj = kk+1,nn-1
                        array(jj) = array(jj+1)
                    end do
                    array(nn) = 0
                end do
            end do 
            if (array(3) == 0) array(1:2) = 0
            return
        end subroutine bubblesortAndTrim
                        
                    
                    



    end module Lib_TriplePoints

