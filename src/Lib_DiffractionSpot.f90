
    module Lib_DiffractionSpot  
!---^^^^^^^^^^^^^^^^^^^^^^^^^^
!*      specialised fourier transform code
!*      to find the first moment of the intensity of a diffraction spot
!*          <k> = int k |F(k)|^2 d3k / int |F(k)|^2 d3k 
!*      where the integral is over a sphere centred on k0
!*      and the fourier transform is over discrete atomic positions
!*          F(k) = 1/sqrt(N) sum_j Exp[ i k.xj ]
!*
!*      Can be run in parallel if compiled with -D MPI
!*  


        use iso_fortran_env
        implicit none
        private
#ifdef MPI 
        include 'mpif.h'
#endif        
        

        
        public      ::      firstMoment
        
    contains
!---^^^^^^^^


        subroutine firstMoment( x, k0,rk,nk , kbar , distributedAtoms)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      compute the first moment in sphere range |k-k0| < |rk| using nk divisions per axis
    !*      note: can run in parallel
    !*      note: two parallel input options - either x(1:3,1:nAtoms) is available on all processes, 
    !*      or x(1:3,1:nLocal) is available on each process, with sum nLocal = nAtoms
            real(kind=real64),dimension(:,:),intent(in)         ::      x           !   atom positions
            real(kind=real64),dimension(3),intent(in)           ::      k0          !   centre of spot
            real(kind=real64),intent(in)                        ::      rk          !   radius of spot
            integer,intent(in)                                  ::      nk          !   subdivisions
            real(kind=real64),dimension(3),intent(out)          ::      kbar        !   first moment
            logical,intent(in),optional                         ::      distributedAtoms
            
            
            complex(kind=real32),dimension(:,:,:,:),pointer   ::      FT   !   Fourier transform at discrete k points
            complex(kind=real32),dimension(0:nk)              ::      G1   !   Fourier transform at discrete k points along 1,2,3 directions
            complex(kind=real32),dimension(0:nk)              ::      G2   !   Fourier transform at discrete k points along 1,2,3 directions
            complex(kind=real32),dimension(0:nk)              ::      G3   !   Fourier transform at discrete k points along 1,2,3 directions
            real(kind=real32)           ::      dk                  !   increment of k along 1,2,3 directions
            integer                     ::      ii,nAtoms           !   atom counters
            real(kind=real32)           ::      x1,x2,x3            !   position of atom in 1,2,3 directions
            real(kind=real32),dimension(8)           ::      FTmod2              !   modulus squared of FT
            
            complex(kind=real32)        ::      eik0x,eik0k3x,eik0mk3x,eik0k3k2x,eik0k3mk2x,eik0mk3k2x,eik0mk3mk2x
            real(kind=real32)           ::      kdotx,cc,ss,k1p,k1m,k2p,k2m,k3p,k3m,ftmodtot
            integer                     ::      jj,kk,ll,mk2
            
            complex(kind=real32),dimension(:,:),pointer        ::      FT_fetch

            integer,dimension(4)        ::      imax
            
        !---    parallel version
            integer                     ::      rank = 0
#ifdef MPI 
            logical                     ::      distribute
            integer                     ::      nproc, ierror
            complex(kind=real32),dimension(:,:,:,:),pointer      ::      FT_tmp   !   Fourier transform at discrete k points parallel copy for reduce op
            real(kind=real32),dimension(4)                       ::      zz,zz_tmp
            call MPI_COMM_SIZE(MPI_COMM_WORLD, nproc, ierror)
            call MPI_COMM_RANK(MPI_COMM_WORLD, rank, ierror)
            distribute = .false.
            if (present(distributedAtoms)) distribute = distributedAtoms
            if (rank==0) print *,"Lib_DiffractionSpot::firstMoment info - running with ",nproc," processes"     
            if (rank==0) print *,"Lib_DiffractionSpot::firstMoment info - distributed atoms? ",distribute     
            allocate(FT_tmp(8,0:nk,0:nk,0:nk))        
            FT_tmp = 0.0
#endif
   
            
            
            
            allocate(FT(8,0:nk,0:nk,0:nk))          !   8 for +/- k in each direction 
                                                    !   FT(1,j,k,l) = F( j, k, l)
                                                    !   FT(2,j,k,l) = F(-j, k, l)
                                                    !   FT(3,j,k,l) = F( j,-k, l)
                                                    !   FT(4,j,k,l) = F(-j,-k, l)
                                                    !   FT(5,j,k,l) = F( j, k,-l)
                                                    !   FT(6,j,k,l) = F(-j, k,-l)
                                                    !   FT(7,j,k,l) = F( j,-k,-l)
                                                    !   FT(8,j,k,l) = F(-j,-k,-l)
                                                    !   will construct all 8, so I have F(0,k,l) and F(-0,k,l)
                                                    !   then will set the F(-0,k,l) etc to zero
                                                   
                                                      
                                                    
            FT = 0.0
            
            G1(0) = cmplx( 1.0,0.0,kind=real32 )     
            G2(0) = cmplx( 1.0,0.0,kind=real32 )     
            G3(0) = cmplx( 1.0,0.0,kind=real32 )     
            
            dk = real(rk/nk,kind=real32)      
            nAtoms = size(x,dim=2)
            
            print *,"searching ",nk,rk," = ",k0 - nk*dk,":",k0 + nk*dk
            
            do ii = 1,nAtoms
#ifdef MPI                
                if ((.not. distribute) .and.  mod(ii-1,nproc)/=rank) cycle
#endif
                x1 = real(x(1,ii),kind=real32)
                x2 = real(x(2,ii),kind=real32) 
                x3 = real(x(3,ii),kind=real32)      !   extract position of atoms
                
            !---    compute Ft at k0
                kdotx = real( k0(1)*x1 + k0(2)*x2 + k0(3)*x3 , kind=real32 )
                eik0x = exp( cmplx(0.0,kdotx) ) 
                              
            !---    compute FT at regular spaced intervals
                do jj = 1,nk
                    kdotx = x1*dk*jj
                    G1(jj) = exp( cmplx(0.0,kdotx) )
                end do
                do jj = 1,nk
                    kdotx = x2*dk*jj
                    G2(jj) = exp( cmplx(0.0,kdotx) )
                end do                    
                do jj = 1,nk
                    kdotx = x3*dk*jj
                    G3(jj) = exp( cmplx(0.0,kdotx) )
                end do
                
            !---    combine to find FT in sphere. Note I am constructing both F(0,k,l) and F(-0,k,l) here.
                do ll = 0,nk
                    eik0k3x = eik0x * G3(ll)
                    eik0mk3x = eik0x * conjg( G3(ll) )
                    do kk = 0,nk 
                        eik0k3k2x   = eik0k3x  * G2(kk)
                        eik0k3mk2x  = eik0k3x  * conjg(G2(kk))
                        eik0mk3k2x  = eik0mk3x * G2(kk)
                        eik0mk3mk2x = eik0mk3x * conjg(G2(kk))
                        mk2 = nk*nk-kk*kk-ll*ll
                        
#ifdef MPI                         
                        FT_Fetch => FT_tmp(:,:,kk,ll)
#else
                        FT_Fetch => FT(:,:,kk,ll)
#endif                        
                        do jj = 0,nk
                            if (jj*jj > mk2) exit    !   strictly compute in sphere                 
                            
                            FT_Fetch(1,jj+1) = FT_fetch(1,jj+1) + eik0k3k2x  *G1(jj)
                            FT_Fetch(2,jj+1) = FT_fetch(2,jj+1) + eik0k3k2x  *conjg(G1(jj))
                            FT_Fetch(3,jj+1) = FT_fetch(3,jj+1) + eik0k3mk2x *G1(jj)       
                            FT_Fetch(4,jj+1) = FT_fetch(4,jj+1) + eik0k3mk2x *conjg(G1(jj))
                            FT_Fetch(5,jj+1) = FT_fetch(5,jj+1) + eik0mk3k2x *G1(jj)       
                            FT_Fetch(6,jj+1) = FT_fetch(6,jj+1) + eik0mk3k2x *conjg(G1(jj))
                            FT_Fetch(7,jj+1) = FT_fetch(7,jj+1) + eik0mk3mk2x*G1(jj)       
                            FT_Fetch(8,jj+1) = FT_fetch(8,jj+1) + eik0mk3mk2x*conjg(G1(jj))
                            
                        end do
                    end do
                end do                                    
                
            end do
            
#ifdef MPI 
            do ll = 0,nk
                kk = mod(ll,nproc)
                call MPI_REDUCE( FT_tmp(:,:,:,ll),FT(:,:,:,ll),8*(nk+1)**2,MPI_COMPLEX,MPI_SUM,kk,MPI_COMM_WORLD,ierror ) 
                if (rank==kk) then
                    FT(2,0,:,ll) = 0.0
                    FT(4,0,:,ll) = 0.0
                    FT(6,0,:,ll) = 0.0
                    FT(8,0,:,ll) = 0.0
                    FT(3,:,0,ll) = 0.0
                    FT(4,:,0,ll) = 0.0
                    FT(7,:,0,ll) = 0.0
                    FT(8,:,0,ll) = 0.0
                    if (ll==0) then
                        FT(5,:,:,0) = 0.0
                        FT(6,:,:,0) = 0.0
                        FT(7,:,:,0) = 0.0
                        FT(8,:,:,0) = 0.0
                    end if
                    !print *,"ll,kk,rank ",ll,kk,rank," FT(1:8,0,0,ll) ",FT(1:8,0,0,ll)
                 end if
            end do  
!           call MPI_ALLREDUCE( FT_tmp,FT,8*(nk+1)**3,MPI_COMPLEX,MPI_SUM,MPI_COMM_WORLD,ierror )
!#endif            
#else
           FT(2,0,:,:) = 0.0
           FT(4,0,:,:) = 0.0
           FT(6,0,:,:) = 0.0
           FT(8,0,:,:) = 0.0
           FT(3,:,0,:) = 0.0
           FT(4,:,0,:) = 0.0
           FT(7,:,0,:) = 0.0
           FT(8,:,0,:) = 0.0
           FT(5,:,:,0) = 0.0
           FT(6,:,:,0) = 0.0
           FT(7,:,:,0) = 0.0
           FT(8,:,:,0) = 0.0
#endif            
        
            
           ! imax = maxloc( abs(FT*conjg(FT)) ) ; imax(2:4) = imax(2:4)-1 
           ! select case(imax(1))
           !     case(1)
           !         jj = imax(2) ; kk = imax(3) ; ll = imax(4)
           !     case(2)                                
           !         jj =-imax(2) ; kk = imax(3) ; ll = imax(4)
           !     case(3)                                
           !         jj = imax(2) ; kk =-imax(3) ; ll = imax(4)
           !     case(4)                                
           !         jj =-imax(2) ; kk =-imax(3) ; ll = imax(4)
           !     case(5)                                
           !         jj = imax(2) ; kk = imax(3) ; ll =-imax(4)
           !     case(6)                                
           !         jj =-imax(2) ; kk = imax(3) ; ll =-imax(4)
           !     case(7)                                
           !         jj = imax(2) ; kk =-imax(3) ; ll =-imax(4)
           !     case(8)                                
           !         jj =-imax(2) ; kk =-imax(3) ; ll =-imax(4)
           ! end select
           ! kbar = k0 + (/jj,kk,ll/)*dk
           ! eik0x =   FT(imax(1),imax(2),imax(3),imax(4))
           ! eik0k3x = directFT( x, kbar )
           ! print *,"max val FT ",imax,kbar,abs(eik0x*conjg(eik0x)),abs(eik0k3x*conjg(eik0k3x))
                            
            
        !---    have now computed FT(0:nk), can compute first moment
            kbar = 0.0
            ftmodtot = 0.0
            do ll = 0,nk
#ifdef MPI
                if (mod(ll,nproc)/=rank) cycle
#endif                
                k3p = real( k0(3) + ll*dk , kind=real32)
                k3m = real( k0(3) - ll*dk , kind=real32)   
                do kk = 0,nk
                    k2p = real( k0(2) + kk*dk , kind=real32)  
                    k2m = real( k0(2) - kk*dk , kind=real32)  
                    mk2 = nk*nk-kk*kk-ll*ll    
                    do jj = 0,nk
                        if (jj*jj > mk2) exit    !   strictly compute in sphere                 
                        k1p = real( k0(1) + jj*dk , kind=real32)  
                        k1m = real( k0(1) - jj*dk , kind=real32)      
                    
                    !---    compute modulus squared of FT
                        x1 = 0.0
                        do ii = 1,8
                            cc = real (FT(ii,jj,kk,ll))
                            ss = aimag(FT(ii,jj,kk,ll))
                            x2 = cc*cc + ss*ss
                            FTmod2(ii) = x2
                            x1 = x1 + x2
                        end do 
                        ftmodtot = ftmodtot + x1
                            
                    !---    add to mean                    
                         kbar(1) = kbar(1) + k1p*( FTmod2(1)+FTmod2(3)+FTmod2(5)+FTmod2(7) )                 &
                                           + k1m*( FTmod2(2)+FTmod2(4)+FTmod2(6)+FTmod2(8) )
                         kbar(2) = kbar(2) + k2p*( FTmod2(1)+FTmod2(2)+FTmod2(5)+FTmod2(6) )                 &
                                           + k2m*( FTmod2(3)+FTmod2(4)+FTmod2(7)+FTmod2(8) )
                         kbar(3) = kbar(3) + k3p*( FTmod2(1)+FTmod2(2)+FTmod2(3)+FTmod2(4) )                 &
                                           + k3m*( FTmod2(5)+FTmod2(6)+FTmod2(7)+FTmod2(8) )
                    end do
                end do
            end do      
            

#ifdef MPI 
            zz(1:3) = kbar(1:3)
            zz(4) = ftmodtot
            !call MPI_REDUCE( zz,zz_tmp,4,MPI_REAL,MPI_SUM,0,MPI_COMM_WORLD,ierror )
            call MPI_ALLREDUCE( zz,zz_tmp,4,MPI_REAL,MPI_SUM,MPI_COMM_WORLD,ierror )
            kbar(1:3) = zz_tmp(1:3)
            ftmodtot = zz_tmp(4)
#endif                        
            
            !print *," ftmodtot ",ftmodtot," on rank ",rank
            
            
            ftmodtot = 1.0/ftmodtot
            kbar = kbar * ftmodtot
             
            
!!           do ll = -nk,nk
!            ll=3
!               do kk = -nk,nk
!                   do jj = -nk,nk
!                       if (jj<0) then
!                           if (kk<0) then
!                                if (ll<0) then
!                                   eik0x = FT(8,-jj,-kk,-ll)
!                                else
!                                   eik0x = FT(4,-jj,-kk,ll)
!                                end if
!                           else
!                               if (ll<0) then
!                                   eik0x = FT(6,-jj,kk,-ll)                                 
!                               else
!                                   eik0x = FT(2,-jj,kk,ll)
!                               end if
!                           end if       
!                       else
!                           if (kk<0) then
!                                if (ll<0) then
!                                   eik0x = FT(7,jj,-kk,-ll)                                                                  
!                                else
!                                   eik0x = FT(3,jj,-kk,ll)                                                                  
!                                end if
!                           else
!                               if (ll<0) then
!                                   eik0x = FT(5,jj,kk,-ll)                                                                  
!                               else
!                                   eik0x = FT(1,jj,kk,ll)                                                                  
!                               end if
!                           end if             
!                       end if 
!                       write (*,fmt='(f10.2,a)',advance="No") abs(eik0x*conjg(eik0x))," "
!                   end do
!                   print *,""
!               end do
!               print *,""
!               print *,""
!!           end do

            deallocate(FT)
            
#ifdef MPI             
            deallocate(FT_tmp)
#endif            
            
            !eik0k3x = directFT( x, kbar )
            !print *,"max val FT ",kbar,abs(eik0k3x*conjg(eik0k3x))
            !stop
            return
        end subroutine firstMoment
            
            
            

        pure function directFT( x, k0 ) result ( FT )
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      compute the direct fourier transform at point k0
             real(kind=real64),dimension(:,:),intent(in)         ::      x           !   atom positions
            real(kind=real64),dimension(3),intent(in)           ::      k0          !   centre of spot
            complex(kind=real64)                    ::      FT   !   Fourier transform at discrete k points
            
            
            integer                     ::      ii,nAtoms           !   atom counters
            real(kind=real64)           ::      x1,x2,x3            !   position of atom in 1,2,3 directions
            
            real(kind=real64)           ::      kdotx
                
            FT = 0.0
            
            nAtoms = size(x,dim=2)
            
            do ii = 1,nAtoms

                x1 = x(1,ii)
                x2 = x(2,ii)
                x3 = x(3,ii)     !   extract position of atoms
                
            !---    compute Ft at k0
                kdotx =  k0(1)*x1 + k0(2)*x2 + k0(3)*x3 
                FT = FT + exp( cmplx(0.0d0,kdotx) )               
            end do
            
            return
        end function directFT
            
            
    end module Lib_DiffractionSpot        
    
!   gfortran -ffree-line-length-256 -Og -pg -g -p ${MYF90LIB}/Lib_Callipers.f90 Lib_DiffractionSpot.f90 -o Lib_DiffractionSpot.exe
!
!    program testLib_DiffractionSpot
!!---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
!        use Lib_DiffractionSpot
!        use Lib_Callipers
!        use iso_fortran_env
!        implicit none
!#ifdef MPI 
!        include 'mpif.h'
!#endif        
!        
!        real(kind=real64),parameter     ::      PI = 3.141592654d0
!        integer       ::      Nx = 40
!        real(kind=real64),dimension(:,:),allocatable        ::      xx
!        integer                 ::      ix,iy,iz,ik,ii
!        
!        real(kind=real64),dimension(3)      ::      k0,kbar
!        real(kind=real64)                   ::      rk
!        integer                             ::      nk
!        character(len=256)      ::      dummy
!        type(Callipers)         ::      tt
!        
!        integer                 ::      rank = 0
!        integer                 ::      nAtoms
!#ifdef MPI 
!           integer                      ::      nproc, ierror
!           call MPI_INIT(ierror)
!           call MPI_COMM_SIZE(MPI_COMM_WORLD, nproc, ierror)
!           call MPI_COMM_RANK(MPI_COMM_WORLD, rank, ierror)
!           if (rank==0) print *,"Lib_DiffractionSpot.exe info - running with ",nproc," processes"     
!#endif        
!        
!        if (rank==0) then
!            print *,"test program usage:"
!            print *,"   ./testLib_DiffractionSpot.exe k0x,k0y,k0z rk nk nx"
!        end if
!        
!        call get_command_argument(1,dummy)        
!        read(dummy,fmt=*)   k0
!        call get_command_argument(2,dummy)        
!        read(dummy,fmt=*)   rk
!        call get_command_argument(3,dummy)        
!        read(dummy,fmt=*)   nk
!        
!        call get_command_argument(4,dummy)        
!        read(dummy,fmt=*)   nx
!        
!        if (rank==0) print *,"read k0,rk,nk,nx = ",k0,rk,nk,nx 
!        
!        
!        
!        
!        k0 = 2*PI*k0
!        rk = rk*2*PI
!
!        allocate(xx(3,2*Nx*Nx*Nx))
!        if (rank==0) print *,"nAtoms = ",2*Nx*Nx*Nx
!        !call random_number(xx) ; xx = (2*xx - 1)*0.1d0
!        nAtoms = 0
!        ii = 0
!        do iz = 0,Nx-1
!            do iy = 0,Nx-1
!                do ix = 0,Nx-1
!                    do ik = 0,1
!                        ii = ii + 1
!#ifdef MPI 
!                        if (mod(ii-1,nProc)/=rank) cycle      !   don't put this atom on the process                  
!#endif
!                        nAtoms = nAtoms + 1
!                        xx(1:3,nAtoms) = xx(1:3,nAtoms) + (/ ix,iy,iz /) + 0.5d0*ik
!                    end do
!                end do
!            end do
!        end do
!        print *,"nAtoms (rank=",rank,") = ",nAtoms
!        
!        !k0 = 2*PI*(/ 1,0,0 /)
!        !rk = 2*PI*0.2d0
!        !nk = 6
!        tt = Callipers_ctor()
!        call firstMoment( xx(:,1:nAtoms), k0,rk,nk , kbar , distributedAtoms = .true.)
!        call pause(tt)
!        
!#ifdef MPI             
!            call MPI_FINALIZE(ierror)
!#endif            
!        if (rank==0) then
!            print *,"k0 = ",k0
!            print *,"rk,nk = ",rk,nk
!            
!            print *,"kbar = ",kbar
!            
!            print *,"kbar/2pi = ",kbar/(2*PI)
!            
!            print *,"elapsed ",elapsed(tt)
!            
!            print *,""
!            print *,"done"
!            print *,""
!        end if    
!        
!        
!            
!    end program testLib_DiffractionSpot    