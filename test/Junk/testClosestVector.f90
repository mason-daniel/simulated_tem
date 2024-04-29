
        
!   gfortran -ffree-line-length-256 ${MYF90LIB}/Lib_Callipers.f90 ${MYF90LIB}/Lib_RotationMatrices.f90 Lib_ClosestVector.f90       
        
    program testLib_ClosestVector
!---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
!*      test against bcc lattice
        use iso_fortran_env
        use Lib_ClosestVector       
        use Lib_RotationMatrices 
        use Lib_Callipers        
        implicit none
        
        integer,parameter                       ::      NRNDTRIALS = 10000000
        integer,parameter                       ::      NTRIALS = 1000000 
        integer,parameter                       ::      NSETS = 100
        
        real(kind=real64),parameter             ::      RMAX = 1.7d0
        real(kind=real64),dimension(3,1000)     ::      vv
        real(kind=real64),dimension(:,:,:),allocatable     ::      vset
        real(kind=real64),dimension(:,:),allocatable     ::      wset
        real(kind=real64),dimension(3)          ::      xx
        integer                                 ::      ix,iy,iz,ik,mm
        real(kind=real64),dimension(3,3)        ::      RR
        
        type(ClosestVector)                     ::      cc
        type(Callipers)                         ::      tt
        real(kind=real64)                       ::      t1,t2,dd,dx
        
        
        integer,dimension(NSETS)                ::      kk
        real(kind=real64),dimension(NSETS)      ::      d2
        
    !---    find a set of target vectors in a sphere         
        mm = 0                      !   number of target vectors
        do iz = -ceiling(RMAX),ceiling(RMAX)
            do iy = -ceiling(RMAX),ceiling(RMAX)
                do ix = -ceiling(RMAX),ceiling(RMAX)
                    do ik = 0,1
                        xx = (/ix,iy,iz/) + 0.5d0*ik        !   bcc lattice 
                        if (norm2(xx) <= RMAX+1.0d-12) then
                            mm = mm + 1
                            vv(1:3,mm) = xx(1:3)
                            print *,"v ",mm,xx , norm2(xx)
                        end if
                    end do
                end do
            end do
        end do
        
        
    !---    rotate target vectors randomly
        allocate(vset(3,mm,NSETS))
        do ik = 1,NSETS
            call random_number(xx) ; xx = xx*3.141592654d0      !   xx is Euler angles
            RR = RotationMatrix_ctor( xx(1),xx(2),xx(3) )       !   xx(1) is rotation about x axis. But that doesn't matter much here.        
            vset(1:3,1:mm,ik) = rotateVector(RR,vv(1:3,1:mm)) 
            call zeroCentreOfMass( vset(1:3,1:mm,ik) )          !   zeroCentreOfMass shouldn't do anything here, but is good practice generally
        end do
        
        
    !---    now construct a closestVector object
        
        cc = ClosestVector_ctor( vset )       
        call report(cc)
        
        
    !---    find time taken to generate random numbers in a sphere.
        
        tt = Callipers_ctor() 
        do ik = 1,NRNDTRIALS
            xx = rndVecInSphere()*RMAX*1.1d0          
        end do
        t1 = elapsed(tt)/NRNDTRIALS       
 
        print *,"time for random vector ",t1
       
    !---    find time taken for nearest vector
        dx = 0
        tt = Callipers_ctor() 
        do ik = 1,NTRIALS
        
            xx = rndVecInSphere()*RMAX*1.1d0
            
            call closestTargetVector(cc,xx,1,ix,dd)
            dx = dx + dd             
        end do
        t2 = elapsed(tt) 
        t2 = t2/NTRIALS - t1
        print *,"time for closest vector ",t2
        
        print *,"mean square distance ",dx/NTRIALS
        
        
        print *,""
        print *,"complex algorithm"
        cc = ClosestVector_ctor( vset,10 )       
        call report(cc)
        dx = 0
        tt = Callipers_ctor() 
        do ik = 1,NTRIALS
        
            xx = rndVecInSphere()*RMAX*1.1d0
            
            call closestTargetVector(cc,xx,kk,d2)
            !
        end do
        t2 = elapsed(tt) 
        t2 = t2/NTRIALS  - t1
        print *,"time for closest vector ",t2/NSETS
        
        tt = Callipers_ctor() 
        dx = 0.0d0
        do ik = 1,NTRIALS
            dx = dx + 1.0d0
        end do
        t2 = elapsed(tt) /NTRIALS 
        print *,"dx ",dx
        print *,"time for add ",t2
       
        
        
!         do ik = 1,NTRIALS
!         
!             xx = rndVecInSphere()*RMAX*1.1d0
!             call closestTargetVector(cc,xx,1,iy,dd )
!             call closestTargetVector(cc,xx,1,ix,dx)
!             
!             if (ix/=iy) then
!                 print *,"error ",iy,dd,ix,dx 
!                
!                 print *,""
!                 stop
!             end if
!             
!                       
!         end do
        
        
        
        allocate(wset(3,mm))
        dd = 0
        tt = Callipers_ctor()
        do ik = 1,NRNDTRIALS
        
            call random_number(dd)
            ix = floor(dd*NSETS) + 1
            wset(1:3,1:mm) = vset(1:3,1:mm,ix)
            dd = dd + sum(wset) 
                         
        end do
        t1 = elapsed(tt)/NRNDTRIALS
        print *,dd,t1
        
        
        dx = 0
        tt = Callipers_ctor()
        do ik = 1,(NTRIALS/100)
        
            call random_number(dd)
            ix = floor(dd*NSETS) + 1
            wset(1:3,1:mm) = vset(1:3,1:mm,ix)
            call closestTargetVectorSet(cc,wset,0.01d0,iy,dd,iz )
               
            dx = dx + abs(ix-iy)
                      
        end do
        t2 = elapsed(tt)/(NTRIALS/100) - t1
        print *,dx,t2
        
        
        
        
        
        
        
        
        
        
        
        call delete(cc)
        print *,""
        print *,"done"
        print *,""
        
    contains
!---^^^^^^^^

        function rndVecInSphere() result(xx)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
            real(kind=real64),dimension(3)  ::      xx
            real(kind=real64)               ::      dd
            do 
                call random_number(xx)
                xx = 2*xx - 1
                dd = xx(1)*xx(1) + xx(2)*xx(2) + xx(3)*xx(3)
                if (dd <= 1.0d0) exit

            end do 
            return
        end function rndVecInSphere
            
        
    end program testLib_ClosestVector        
        
        


    
    