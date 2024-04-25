
    module Lib_MitchellsBestCandidate
!---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
!*      provides functionality for Mitchells Best Candidate algorithm
!*      used to produce "blue" noise
!*      1) Start with N>=1 fixed points
!*      2) Choose M candidate points at random
!*      3) Using a distance metric, pick the candidate which is furthest from all others. Reject the others
!*      4) Add this best candidate to the fixed point list and return to 1.
!*     
!*      You should increase M as a function of N.
!*      This algorithm is written in a very general way, so uses an external function 
!*          d = distance( x1,x2 )
!*      or
!*          d = minDistance( x1(:,:),x2 )
!*      for the metric, and
!*          x = trial(d) to generate a candidate, where d is dimensionality
!*      These will have to be aware of the boundary conditions.

        use Lib_SimpleProgressBar
        use iso_fortran_env
        implicit none
        private
         
        public      ::      bestCandidate_pairDist
        public      ::      bestCandidate_minDist
             
    
    contains
!---^^^^^^^^

        subroutine bestCandidate_pairDist( x_fixed, a, distance, trial, x_best )
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      select the best candidate from M = a N + 1 trials.
            real(kind=real64),dimension(:,:),intent(in)             ::      x_fixed
            integer,intent(in)                                      ::      a
            real(kind=real64),dimension(:),intent(out)              ::      x_best
            real(kind=real64),dimension(size(x_fixed,dim=1))        ::      x_trial
            
                
         interface
         
             function distance(x1,x2) result(d)
                 use iso_fortran_env
                 real(kind=real64),dimension(:),intent(in)       ::      x1,x2
                 real(kind=real64)                               ::      d
             end function distance
         
             function trial(d) result(x)
                 use iso_fortran_env
                 integer,intent(in)                              ::      d
                 real(kind=real64),dimension(d)                  ::      x
             end function trial
         
         end interface
        
        
            
            integer             ::      ii,jj,nn,mm,dd
            real(kind=real64)   ::      dij,dmin,dbest
            
            dd = size(x_fixed,dim=1)        !   dimensionality
            nn = size(x_fixed,dim=2)        !   number of fixed points
            
            if (nn==0) then
                x_best(1:dd) = trial(dd)
                return
            end if
            
        !---    generate trial points            
            mm = a*nn + 1                   !   number of trial points
            
            dbest = 0.0d0
            do jj = 1,mm
                x_trial(1:dd) = trial(dd)
                dmin = huge(1.0)
                do ii = 1,nn
                    dij = distance( x_fixed(:,ii),x_trial(:) )
                    dmin = min(dmin,dij)
                end do
                
                if (dmin > dbest) then
                    dbest = dmin
                    x_best(1:dd) = x_trial(1:dd)
                end if
            end do
          
            return
        end subroutine bestCandidate_pairDist
            
       
        subroutine bestCandidate_minDist( x_fixed, a, minDistance, trial, x_best , d )
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      select the best candidate from M = a N + 1 trials.
            real(kind=real64),dimension(:,:),intent(in)             ::      x_fixed
            integer,intent(in)                                      ::      a
            real(kind=real64),dimension(:),intent(out)              ::      x_best
            real(kind=real64),dimension(size(x_fixed,dim=1))        ::      x_trial
            real(kind=real64),intent(out),optional                  ::      d
                
         interface
         
             function minDistance(x1,x2) result(d)
                 use iso_fortran_env
                 real(kind=real64),dimension(:,:),intent(in)     ::      x1
                 real(kind=real64),dimension(:),intent(in)       ::      x2
                 real(kind=real64)                               ::      d
             end function minDistance
         
             function trial(d) result(x)
                 use iso_fortran_env
                 integer,intent(in)                              ::      d
                 real(kind=real64),dimension(d)                  ::      x
             end function trial
         
         end interface
        
        
            
            integer             ::      jj,nn,mm,dd
            real(kind=real64)   ::      dmin,dbest
            
            dd = size(x_fixed,dim=1)        !   dimensionality
            nn = size(x_fixed,dim=2)        !   number of fixed points
            
            if (nn==0) then
                x_best(1:dd) = trial(dd)
                return
            end if
            
        !---    generate trial points            
            mm = a*nn + 1                   !   number of trial points
            
            dbest = 0.0d0
            do jj = 1,mm
                x_trial(1:dd) = trial(dd)
                dmin = minDistance( x_fixed(:,:),x_trial(:) )
                
                if (dmin > dbest) then
                    dbest = dmin
                    x_best(1:dd) = x_trial(1:dd)
                end if
                
            end do
            !print *,"bestCandidate_minDist accepts min dist ",dbest
          
            if (present(d)) d = dbest
            return
        end subroutine bestCandidate_minDist
            
           


    end module Lib_MitchellsBestCandidate
    