
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

        use OMP_LIB
        use iso_fortran_env
        implicit none
        private
        
        integer,public      ::      DIFFSPOTS_M = 1

        
        public      ::      firstMoment,findMaximum
        
        interface           intensity
            !module procedure    intensity0
            !module procedure    intensity1
            module procedure    intensity2
        end interface
        
    contains
!---^^^^^^^^


        subroutine firstMoment( x, k0,rk,nk , kbar )
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      compute the first moment in sphere range |k-k0| < |rk| using nk divisions per axis
    !*      note: can run in parallel
    !*      note: two parallel input options - either x(1:3,1:nAtoms) is available on all processes, 
    !*      or x(1:3,1:nLocal) is available on each process, with sum nLocal = nAtoms
            real(kind=real64),dimension(:,:),intent(in)         ::      x           !   atom positions
            real(kind=real64),dimension(3),intent(in)           ::      k0          !   centre of spot
            real(kind=real64),intent(in)                        ::      rk          !   radius of spot
            integer,intent(in)                                  ::      nk          !   subdivisions
            real(kind=real64),dimension(3),intent(out)          ::      kbar        !   first moment
           
            complex(kind=real64),dimension(:,:,:,:,:),allocatable   ::      FT   !   Fourier transform at discrete k points
            complex(kind=real64),dimension(0:nk)              ::      G1,G2,G3  !   Fourier transform at discrete k points along 1,2,3 directions
            real(kind=real64)           ::      dk                  !   increment of k along 1,2,3 directions
            integer                     ::      ii,nAtoms           !   atom counters
            real(kind=real64)           ::      x1,x2,x3            !   position of atom in 1,2,3 directions
            real(kind=real64),dimension(8)           ::      FTmod2              !   modulus squared of FT
            
            complex(kind=real64)        ::      eik0x,eiknx,eiknk3x,eiknkk3x,eiknk3k2x,eiknk3mk2x,eiknkk3k2x,eiknkk3mk2x
            real(kind=real64)           ::      kdotx,k1p,k1m,k2p,k2m,k3p,k3m,ftmodtot
            integer                     ::      jj,kk,ll,mk2,nn
            complex(kind=real64),dimension(:,:,:,:,:),pointer    ::      FT_tmp
            complex(kind=real64),dimension(:,:),pointer          ::      FT_fetch
            real(kind=real64),dimension(3)      ::      kbar_part
            real(kind=real64)                   ::      ftmod_part
            
            
            !real(kind=real64),dimension(-nk:nk,3)      ::      ibar 
            !real(kind=real64),dimension(0:nk)          ::      sumi
            
            
            real(kind=real64)           ::      k0dotx            
            
            
            allocate(FT(8,0:nk,0:nk,0:nk,DIFFSPOTS_M))          !   8 for +/- k in each direction 
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
            
            dk = real(rk/nk,kind=real64)      
            nAtoms = size(x,dim=2)
            
            !print *,"nAtoms, sum x ",nAtoms, sum(x(1,:)), sum(x(2,:)), sum(x(3,:)) ," dk,nk ",dk,nk
            
            
!$OMP PARALLEL  PRIVATE(FT_tmp,FT_fetch,x1,x2,x3,eik0x,eiknx,ii,jj,kk,ll,nn,mk2,G1,G2,G3,eiknk3x,eiknkk3x,eiknk3k2x,eiknkk3k2x,eiknk3mk2x,eiknkk3mk2x) 
!!!!  !$OMP PARALLEL  PRIVATE(FT_tmp,FT_fetch,x1,x2,x3,k0dotx,kdotx,ii,jj,kk,ll,nn,mk2) 
            allocate(FT_tmp(8,0:nk,0:nk,0:nk,DIFFSPOTS_M))
            FT_tmp = 0.0 
           G1(0) = cmplx( 1.0d0,0.0d0,kind=real64 )     
           G2(0) = cmplx( 1.0d0,0.0d0,kind=real64 )     
           G3(0) = cmplx( 1.0d0,0.0d0,kind=real64 )     
            

    !$OMP DO                            
            do ii = 1,nAtoms

                x1 = x(1,ii) 
                x2 = x(2,ii)  
                x3 = x(3,ii)       !   extract position of atoms
                
            !---    compute Ft at k0
                kdotx = real( k0(1)*x1 + k0(2)*x2 + k0(3)*x3 , kind=real64 )
                eik0x = exp( cmplx(0.0d0,k0(1)*x1 + k0(2)*x2 + k0(3)*x3,kind=real64) ) 
                        
                k0dotx = k0(1)*x1 + k0(2)*x2 + k0(3)*x3
                
                
            !---    combine to find FT in sphere. Note I am constructing both F(0,k,l) and F(-0,k,l) here.
            
                eiknx = cmplx( 1.0d0,0.0d0,kind=real64 )  
                do nn = 1,DIFFSPOTS_M     !   look at k0,2k0,3k0,4k0 spots...
                    eiknx = eiknx*eik0x
            
                    
                      
                !---    compute FT at regular spaced intervals
                   do jj = 1,nk
                       G1(jj) = exp( cmplx(0.0d0,nn*jj*dk*x1,kind=real64) )
                       G2(jj) = exp( cmplx(0.0d0,nn*jj*dk*x2,kind=real64) )
                       G3(jj) = exp( cmplx(0.0d0,nn*jj*dk*x3,kind=real64) )
                   end do
                    
                    
                    
                    
                do ll = 0,nk
                    eiknk3x  = eik0x * G3(ll)
                    eiknkk3x = eik0x * conjg( G3(ll) )
                    do kk = 0,nk 
                        eiknk3k2x   = eiknk3x  * G2(kk)
                        eiknk3mk2x  = eiknk3x  * conjg(G2(kk))
                        eiknkk3k2x  = eiknkk3x * G2(kk)
                        eiknkk3mk2x = eiknkk3x * conjg(G2(kk))
                        mk2 = nk*nk-kk*kk-ll*ll
                        
                        

                        FT_Fetch => FT_tmp(:,:,kk,ll,nn)
                  
                        do jj = 0,nk
                            if (jj*jj > mk2) exit    !   strictly compute in sphere                 
                                    
                          !  kdotx = k0dotx + ( jj*x1+kk*x2+ll*x3)*dk
                          !  FT_Fetch(1,jj+1) = FT_fetch(1,jj+1) + exp( cmplx(0.0d0,nn*kdotx,kind=real64) )
                          !  
                          !  kdotx = k0dotx + (-jj*x1+kk*x2+ll*x3)*dk
                          !  FT_Fetch(2,jj+1) = FT_fetch(2,jj+1) + exp( cmplx(0.0d0,nn*kdotx,kind=real64) )
                          !  
                          !  kdotx = k0dotx + ( jj*x1-kk*x2+ll*x3)*dk
                          !  FT_Fetch(3,jj+1) = FT_fetch(3,jj+1) + exp( cmplx(0.0d0,nn*kdotx,kind=real64) )
                          !  
                          !  kdotx = k0dotx + (-jj*x1-kk*x2+ll*x3)*dk
                          !  FT_Fetch(4,jj+1) = FT_fetch(4,jj+1) + exp( cmplx(0.0d0,nn*kdotx,kind=real64) )
                          !  
                          !  kdotx = k0dotx + ( jj*x1+kk*x2-ll*x3)*dk
                          !  FT_Fetch(5,jj+1) = FT_fetch(5,jj+1) + exp( cmplx(0.0d0,nn*kdotx,kind=real64) )
                          !  
                          !  kdotx = k0dotx + (-jj*x1+kk*x2-ll*x3)*dk
                          !  FT_Fetch(6,jj+1) = FT_fetch(6,jj+1) + exp( cmplx(0.0d0,nn*kdotx,kind=real64) )
                          !  
                          !  kdotx = k0dotx + ( jj*x1-kk*x2-ll*x3)*dk
                          !  FT_Fetch(7,jj+1) = FT_fetch(7,jj+1) + exp( cmplx(0.0d0,nn*kdotx,kind=real64) )
                          !  
                          !  kdotx = k0dotx + (-jj*x1-kk*x2-ll*x3)*dk
                          !  FT_Fetch(8,jj+1) = FT_fetch(8,jj+1) + exp( cmplx(0.0d0,nn*kdotx,kind=real64) )
                                       
                            FT_Fetch(1,jj+1) = FT_fetch(1,jj+1) + eiknk3k2x  *G1(jj)
                            FT_Fetch(2,jj+1) = FT_fetch(2,jj+1) + eiknk3k2x  *conjg(G1(jj))
                            FT_Fetch(3,jj+1) = FT_fetch(3,jj+1) + eiknk3mk2x *G1(jj)       
                            FT_Fetch(4,jj+1) = FT_fetch(4,jj+1) + eiknk3mk2x *conjg(G1(jj))
                            FT_Fetch(5,jj+1) = FT_fetch(5,jj+1) + eiknkk3k2x *G1(jj)       
                            FT_Fetch(6,jj+1) = FT_fetch(6,jj+1) + eiknkk3k2x *conjg(G1(jj))
                            FT_Fetch(7,jj+1) = FT_fetch(7,jj+1) + eiknkk3mk2x*G1(jj)       
                            FT_Fetch(8,jj+1) = FT_fetch(8,jj+1) + eiknkk3mk2x*conjg(G1(jj))
                          
                        end do
                    end do
                end do
                
                end do                                    
                
            end do
      !$OMP END DO         
      
      !$OMP CRITICAL
            FT(:,:,:,:,:) = FT(:,:,:,:,:) + FT_tmp(:,:,:,:,:)
      !$OMP END CRITICAL      

       
       
       deallocate(FT_tmp)
!$OMP END PARALLEL           
        
!$OMP PARALLEL SHARED(FT)
    !$OMP DO
        do jj = 0,nk
            FT(2,0,:,jj,:) = 0.0d0
            FT(4,0,:,jj,:) = 0.0d0
            FT(6,0,:,jj,:) = 0.0d0
            FT(8,0,:,jj,:) = 0.0d0
            FT(3,:,0,jj,:) = 0.0d0
            FT(4,:,0,jj,:) = 0.0d0
            FT(7,:,0,jj,:) = 0.0d0
            FT(8,:,0,jj,:) = 0.0d0
            FT(5,:,jj,0,:) = 0.0d0
            FT(6,:,jj,0,:) = 0.0d0
            FT(7,:,jj,0,:) = 0.0d0
            FT(8,:,jj,0,:) = 0.0d0
        end do
    !$OMP END DO
!$OMP END PARALLEL   



        !---    have now computed FT(0:nk), can compute first moment
            kbar = 0.0d0
            ftmodtot = 0.0d0
!$OMP PARALLEL SHARED(kbar,ftmodtot),PRIVATE(ii,jj,kk,ll,mk2,nn, k3m,k3p,k2m,k2p,k1m,k1p, FTmod2,kbar_part,ftmod_part)   
            kbar_part(1:3) = 0.0d0
            ftmod_part = 0.0d0
    !$OMP DO            
            do ll = 0,nk
              
                k3p = k0(3) + ll*dk ; k3m = k0(3) - ll*dk
                do kk = 0,nk
                    k2p = k0(2) + kk*dk ; k2m = k0(2) - kk*dk 
                    mk2 = nk*nk-kk*kk-ll*ll    
                    do jj = 0,nk
                        if (jj*jj > mk2) exit    !   strictly compute in sphere                 
                        k1p = k0(1) + jj*dk ; k1m = k0(1) - jj*dk  
                    
                    !---    compute modulus squared of FT    
                    !
                    !    FTmod2 = 0.0d0
                    !    do nn = 1,DIFFSPOTS_M                   
                    !        do ii = 1,8
                    !            cc = real (FT(ii,jj,kk,ll,nn))
                    !            ss = aimag(FT(ii,jj,kk,ll,nn))
                    !            FTmod2(ii) = FTmod2(ii) + cc*cc + ss*ss
                    !        end do 
                    !    end do
                         FTmod2(1:8) = intensity( FT(:,jj,kk,ll,:) )                               
                        
                            
                    !---    add to mean                    
                         kbar_part(1) = kbar_part(1) + k1p*( FTmod2(1)+FTmod2(3)+FTmod2(5)+FTmod2(7) )                 &
                                                     + k1m*( FTmod2(2)+FTmod2(4)+FTmod2(6)+FTmod2(8) )
                         kbar_part(2) = kbar_part(2) + k2p*( FTmod2(1)+FTmod2(2)+FTmod2(5)+FTmod2(6) )                 &
                                                     + k2m*( FTmod2(3)+FTmod2(4)+FTmod2(7)+FTmod2(8) )
                         kbar_part(3) = kbar_part(3) + k3p*( FTmod2(1)+FTmod2(2)+FTmod2(3)+FTmod2(4) )                 &
                                                     + k3m*( FTmod2(5)+FTmod2(6)+FTmod2(7)+FTmod2(8) )
                                                     
                         ftmod_part = ftmod_part + sum(FTmod2) 
                         
                    end do
                end do
            end do      
    !$OMP END DO             
    
    !$OMP CRITICAL
                kbar(1:3) = kbar(1:3) + kbar_part(1:3)
                ftmodtot = ftmodtot + ftmod_part
    !$OMP END CRITICAL      

            
            
!$OMP END PARALLEL            

            !print *," ftmodtot ",ftmodtot
            
            ftmodtot = 1.0/ftmodtot
            kbar = kbar * ftmodtot
            
            
!             ibar = 0.0d0
!             sumi = 0.0d0
!             
!                 do ll = 0,nk
!                     do kk = 0,nk
!                         do jj = 0,nk
!                             if (jj*jj+kk*kk+ll*ll>nk*nk) exit 
!                             FTmod2 = 0.0d0
!                             do nn = 1,DIFFSPOTS_M
!                                 do ii = 1,8
!                                     cc = real (FT(ii,jj,kk,ll,nn))
!                                     ss = aimag(FT(ii,jj,kk,ll,nn))
!                                     FTmod2(ii) = FTmod2(ii) + cc*cc + ss*ss
!                                 end do
!                             end do
!                             ibar(jj,1) = ibar(jj,1)   + FTmod2(1)+FTmod2(3)+FTmod2(5)+FTmod2(7)
!                             ibar(-jj,1) = ibar(-jj,1) + FTmod2(2)+FTmod2(4)+FTmod2(6)+FTmod2(8)
!                             ibar(kk,2) = ibar(kk,2)   + FTmod2(1)+FTmod2(2)+FTmod2(5)+FTmod2(6)
!                             ibar(-kk,2) = ibar(-jj,2) + FTmod2(3)+FTmod2(4)+FTmod2(7)+FTmod2(8)
!                             ibar(ll,3) = ibar(ll,3)   + FTmod2(1)+FTmod2(2)+FTmod2(3)+FTmod2(4)
!                             ibar(-ll,3) = ibar(-ll,3) + FTmod2(5)+FTmod2(6)+FTmod2(7)+FTmod2(8)
!                             sumi(jj) = sumi(jj) + 1.0d0
!                         end do
!                     end do
!                 end do             
!                 do jj = 0,nk
!                     cc = 1/(natoms*sumi(jj))
!                     ibar(jj,:) = ibar(jj,:)*cc
!                     if (jj>0) ibar(-jj,:) = ibar(-jj,:)*cc
!                     
!                 end do
!                 do jj = -nk,nk
!                     write (*,fmt='(3f12.5,a,3f16.5)') k0+jj*dk ,"    ", ibar(jj,:)
!                 end do
!                 stop
            
!            do ll = -nk,nk
!                do kk = -nk,nk
!                    do jj = -nk,nk
!                        if (jj<0) then
!                            if (kk<0) then
!                                 if (ll<0) then
!                                    eik0x = FT(8,-jj,-kk,-ll)
!                                 else
!                                    eik0x = FT(4,-jj,-kk,ll)
!                                 end if
!                            else
!                                if (ll<0) then
!                                    eik0x = FT(6,-jj,kk,-ll)                                 
!                                else
!                                    eik0x = FT(2,-jj,kk,ll)
!                                end if
!                            end if       
!                        else
!                            if (kk<0) then
!                                 if (ll<0) then
!                                    eik0x = FT(7,jj,-kk,-ll)                                                                  
!                                 else
!                                    eik0x = FT(3,jj,-kk,ll)                                                                  
!                                 end if
!                            else
!                                if (ll<0) then
!                                    eik0x = FT(5,jj,kk,-ll)                                                                  
!                                else
!                                    eik0x = FT(1,jj,kk,ll)                                                                  
!                                end if
!                            end if             
!                        end if 
!                        write (*,fmt='(f10.2,a)',advance="No") abs(eik0x*conjg(eik0x))," "
!                    end do
!                    print *,""
!                end do
!                print *,""
!                print *,""
!            end do

            deallocate(FT)
                 
            
            
            return
        end subroutine firstMoment
            
            
        subroutine findMaximum( x, k0,rk,nk , kmax,ok )
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      compute the first maximum in sphere range |k-k0| < |rk| using nk divisions per axis
    !*      note: can run in parallel
    !*      note: two parallel input options - either x(1:3,1:nAtoms) is available on all processes, 
    !*      or x(1:3,1:nLocal) is available on each process, with sum nLocal = nAtoms
            real(kind=real64),dimension(:,:),intent(in)         ::      x           !   atom positions
            real(kind=real64),dimension(3),intent(in)           ::      k0          !   centre of spot
            real(kind=real64),intent(in)                        ::      rk          !   radius of spot
            integer,intent(in)                                  ::      nk          !   subdivisions
            real(kind=real64),dimension(3),intent(out)          ::      kmax        !   maximum
            logical,intent(out)                                 ::      ok          !   returns false if max on border of sphere
           
            complex(kind=real64),dimension(:,:,:,:,:),allocatable   ::      FT   !   Fourier transform at discrete k points
            complex(kind=real64),dimension(0:nk)              ::      G1   !   Fourier transform at discrete k points along 1,2,3 directions
            complex(kind=real64),dimension(0:nk)              ::      G2   !   Fourier transform at discrete k points along 1,2,3 directions
            complex(kind=real64),dimension(0:nk)              ::      G3   !   Fourier transform at discrete k points along 1,2,3 directions
            real(kind=real64)           ::      dk                  !   increment of k along 1,2,3 directions
            integer                     ::      ii,nAtoms           !   atom counters
            real(kind=real64)           ::      x1,x2,x3            !   position of atom in 1,2,3 directions
            real(kind=real64),dimension(8)           ::      FTmod2              !   modulus squared of FT
            
            complex(kind=real64)        ::      eik0x,eiknx,eiknk3x,eiknkk3x,eiknk3k2x,eiknk3mk2x,eiknkk3k2x,eiknkk3mk2x
!            real(kind=real64)           ::      cc,ss   !,ftmodtot,kdotx,k1p,k1m,k2p,k2m,k3p,k3m
            integer                     ::      jj,kk,ll,mk2,nn
            complex(kind=real64),dimension(:,:,:,:,:),pointer    ::      FT_tmp
            complex(kind=real64),dimension(:,:),pointer          ::      FT_fetch
            
            integer,dimension(4)        ::      jmax
            real(kind=real64)           ::      imax
             
!            integer     ::      dj1,dj2,dj3
            
            allocate(FT(8,0:nk,0:nk,0:nk,DIFFSPOTS_M))          !   8 for +/- k in each direction 
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
            
            dk = real(rk/nk,kind=real64)      
            nAtoms = size(x,dim=2)
            
            !print *,"nAtoms, sum x ",nAtoms, sum(x(1,:)), sum(x(2,:)), sum(x(3,:)) ," dk,nk ",dk,nk," spots ", DIFFSPOTS_M
            
            
!$OMP PARALLEL  PRIVATE(FT_tmp,FT_fetch,x1,x2,x3,eik0x,eiknx,ii,jj,kk,ll,nn,mk2,G1,G2,G3,eiknk3x,eiknkk3x,eiknk3k2x,eiknkk3k2x,eiknk3mk2x,eiknkk3mk2x) 
            allocate(FT_tmp(8,0:nk,0:nk,0:nk,DIFFSPOTS_M))
            FT_tmp = 0.0 
            G1(0) = cmplx( 1.0d0,0.0d0,kind=real64 )     
            G2(0) = cmplx( 1.0d0,0.0d0,kind=real64 )     
            G3(0) = cmplx( 1.0d0,0.0d0,kind=real64 )     
             

    !$OMP DO                            
            do ii = 1,nAtoms

                x1 = x(1,ii) 
                x2 = x(2,ii)  
                x3 = x(3,ii)       !   extract position of atoms
                
            !---    compute Ft at k0
                !kdotx = real( k0(1)*x1 + k0(2)*x2 + k0(3)*x3 , kind=real64 )
                eik0x = exp( cmplx(0.0d0,k0(1)*x1 + k0(2)*x2 + k0(3)*x3,kind=real64) ) 
                    
            !---    combine to find FT in sphere. Note I am constructing both F(0,k,l) and F(-0,k,l) here.
            
                eiknx = cmplx( 1.0d0,0.0d0,kind=real64 )  
                do nn = 1,DIFFSPOTS_M     !   look at k0,2k0,3k0,4k0 spots...
                    eiknx = eiknx*eik0x
                    
                              
            !---    compute FT at regular spaced intervals
                do jj = 1,nk
                    G1(jj) = exp( cmplx(0.0d0,nn*jj*dk*x1,kind=real64) )
                    G2(jj) = exp( cmplx(0.0d0,nn*jj*dk*x2,kind=real64) )
                    G3(jj) = exp( cmplx(0.0d0,nn*jj*dk*x3,kind=real64) )
                end do
                
            
                do ll = 0,nk
                    eiknk3x  = eik0x * G3(ll)
                    eiknkk3x = eik0x * conjg( G3(ll) )
                    do kk = 0,nk 
                        eiknk3k2x   = eiknk3x  * G2(kk)
                        eiknk3mk2x  = eiknk3x  * conjg(G2(kk))
                        eiknkk3k2x  = eiknkk3x * G2(kk)
                        eiknkk3mk2x = eiknkk3x * conjg(G2(kk))
                        mk2 = nk*nk-kk*kk-ll*ll
                        

                        FT_Fetch => FT_tmp(:,:,kk,ll,nn)
                  
                        do jj = 0,nk
                            if (jj*jj > mk2) exit    !   strictly compute in sphere                 
                                                        
                            FT_Fetch(1,jj+1) = FT_fetch(1,jj+1) + eiknk3k2x  *G1(jj)
                            FT_Fetch(2,jj+1) = FT_fetch(2,jj+1) + eiknk3k2x  *conjg(G1(jj))
                            FT_Fetch(3,jj+1) = FT_fetch(3,jj+1) + eiknk3mk2x *G1(jj)       
                            FT_Fetch(4,jj+1) = FT_fetch(4,jj+1) + eiknk3mk2x *conjg(G1(jj))
                            FT_Fetch(5,jj+1) = FT_fetch(5,jj+1) + eiknkk3k2x *G1(jj)       
                            FT_Fetch(6,jj+1) = FT_fetch(6,jj+1) + eiknkk3k2x *conjg(G1(jj))
                            FT_Fetch(7,jj+1) = FT_fetch(7,jj+1) + eiknkk3mk2x*G1(jj)       
                            FT_Fetch(8,jj+1) = FT_fetch(8,jj+1) + eiknkk3mk2x*conjg(G1(jj))
                          
                        end do
                    end do
                end do
                
                end do                                    
                
            end do
      !$OMP END DO         
      
      !$OMP CRITICAL
            FT(:,:,:,:,:) = FT(:,:,:,:,:) + FT_tmp(:,:,:,:,:)
      !$OMP END CRITICAL      

       
       
       deallocate(FT_tmp)
!$OMP END PARALLEL           
        
!$OMP PARALLEL SHARED(FT)
    !$OMP DO
        do jj = 0,nk
            FT(2,0,:,jj,:) = 0.0d0
            FT(4,0,:,jj,:) = 0.0d0
            FT(6,0,:,jj,:) = 0.0d0
            FT(8,0,:,jj,:) = 0.0d0
            FT(3,:,0,jj,:) = 0.0d0
            FT(4,:,0,jj,:) = 0.0d0
            FT(7,:,0,jj,:) = 0.0d0
            FT(8,:,0,jj,:) = 0.0d0
            FT(5,:,jj,0,:) = 0.0d0
            FT(6,:,jj,0,:) = 0.0d0
            FT(7,:,jj,0,:) = 0.0d0
            FT(8,:,jj,0,:) = 0.0d0
        end do
    !OMP END DO
!$OMP END PARALLEL   



        !---    have now computed FT(0:nk), can compute maximum
             imax = 0.0d0 ; jmax = 0
            
            do ll = 0,nk
                do kk = 0,nk
                    mk2 = nk*nk-kk*kk-ll*ll    
                    do jj = 0,nk
                        if (jj*jj > mk2) exit    !   strictly compute in sphere                 
                        
                    !---    compute modulus squared of FT    
                    !   FTmod2 = 0.0d0
                    !   do nn = 1,DIFFSPOTS_M                   
                    !       do ii = 1,8
                    !           cc = real (FT(ii,jj,kk,ll,nn))
                    !           ss = aimag(FT(ii,jj,kk,ll,nn))
                    !           FTmod2(ii) = FTmod2(ii) + cc*cc + ss*ss
                    !       end do 
                    !   end do
                        FTmod2(1:8) = intensity( FT(1:8,jj,kk,ll,:) )                  
                                                    
                    !---    is this a maximum 
                        do ii = 1,8 
                            if (FTmod2(ii)>imax) then
                                jmax(1:4) = (/ ii,jj,kk,ll /)
                                imax = FTmod2(ii)
                                !print *," max ",ii,jj,kk,ll,imax,FT(ii,jj,kk,ll,:)                   
                                                    
                            end if
                        end do
                        
                    end do
                end do
            end do     
            
            
      !     print *,jmax,intensity( FT(jmax(1),jmax(2),jmax(3),jmax(4),:) )
      !     
      !     do ll = 0,nk
      !         do kk = 0,nk
      !             mk2 = nk*nk-kk*kk-ll*ll    
      !             do jj = 0,nk
      !                 if (jj*jj > mk2) exit    !   strictly compute in sphere                 
      !                 FT(:,jj,kk,ll,1) = intensity( FT(:,jj,kk,ll,:) )
      !             end do
      !         end do
      !     end do
      !     jmax = maxloc( real(FT(:,:,:,:,1)) )
      !     jmax(2:4) = jmax(2:4)-1
      !     
      !     
      !     print *,jmax,intensity( FT(jmax(1),jmax(2),jmax(3),jmax(4),:) )
      !     
            
            
            deallocate(FT)
                                    
            select case(jmax(1))
                case(1)
                    jj = jmax(2) ; kk = jmax(3) ; ll = jmax(4)
                case(2)
                    jj =-jmax(2) ; kk = jmax(3) ; ll = jmax(4)
                case(3)
                    jj = jmax(2) ; kk =-jmax(3) ; ll = jmax(4)
                case(4)
                    jj =-jmax(2) ; kk =-jmax(3) ; ll = jmax(4)
                case(5)
                    jj = jmax(2) ; kk = jmax(3) ; ll =-jmax(4)
                case(6)
                    jj =-jmax(2) ; kk = jmax(3) ; ll =-jmax(4)
                case(7)
                    jj = jmax(2) ; kk =-jmax(3) ; ll =-jmax(4)
                case(8)
                    jj =-jmax(2) ; kk =-jmax(3) ; ll =-jmax(4)
            end select
            
            kmax = k0 + dk*(/jj,kk,ll/)  
            
            ok = ( (jj+1)*(jj+1) + kk*kk + ll*ll < nk*nk )
            
            
            
        !   write (*,fmt='(a,3f12.5,a,f12.5,a,3f12.5,a,3i4,a,i4)') "search ",k0,"+/-",rk," maximum at ",kmax," = ",jj,kk,ll,"/",nk
        !   do dj3 = -2,2
        !       do dj2 = -2,2
        !           do dj1 = -2,2
        !               write (*,fmt='(f16.4)',advance="no") intensity(directFT( x, kmax+(/dj1,dj2,dj3/)*dk ))
        !           end do
        !           print *,""
        !       end do
        !       print *,""
        !   end do 
        !    
        !   
        !   stop
            
           ! kmax = k0  
           ! eik0x =  directFT( x, kmax )              
           ! print *,"direct FT ",kmax,abs(eik0x*conjg(eik0x)) 
           ! 
             
           ! eik0x =  directFT( x, kmax )              
           ! print *,"direct FT ",kmax,abs(eik0x*conjg(eik0x)),imax/DIFFSPOTS_M
            
            !kmax = (/ 4.31613 , -1.99226 , 0.96605 /)
            !eik0x =  directFT( x, kmax )              
            !print *,"direct FT ",kmax,abs(eik0x*conjg(eik0x)) 
            
            
            return
        end subroutine findMaximum
        
!         pure function intensity0( FT ) result( aFT )
!     !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
!     !*      take a sum over the Fourier transform spots and find intensity
!             complex(kind=real64),intent(in)        ::      FT        
!             real(kind=real64)                      ::      aFT
!             real(kind=real64)       ::      cc,ss
!             cc = real (FT)
!             ss = aimag(FT)
!             aFT =  cc*cc + ss*ss
!             return
!         end function intensity0

            
!         pure function intensity1( FT ) result( aFT )
!     !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
!     !*      take a sum over the Fourier transform spots and find intensity
!             complex(kind=real64),dimension(:),intent(in)        ::      FT        
!             real(kind=real64)                                   ::      aFT
!             integer                 ::      nn
!             real(kind=real64)       ::      cc,ss
!             aFT = 0.0d0
!             do nn = 1,size(FT)
!                 cc = real (FT(nn))
!                 ss = aimag(FT(nn))
!                 aFT = aFT + cc*cc + ss*ss
!             end do 
!             return
!         end function intensity1

        pure function intensity2( FT ) result( aFT )
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      take a sum over the Fourier transform spots and find intensity
            complex(kind=real64),dimension(:,:),intent(in)      ::      FT        
            real(kind=real64),dimension(size(FT,dim=1))         ::      aFT
            integer                 ::      nn,ii
            real(kind=real64)       ::      cc,ss
            aFT = 0.0d0
            do ii = 1,size(FT,dim=1)
                do nn = 1,size(FT,dim=2)
                    cc = real (FT(ii,nn))
                    ss = aimag(FT(ii,nn))
                    aFT(ii) = aFT(ii) + cc*cc + ss*ss
                end do 
            end do
            return            
        end function intensity2

        

!         function directFT( x, k0 ) result ( FT )
!     !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
!     !*      compute the direct fourier transform at point k0
!              real(kind=real64),dimension(:,:),intent(in)         ::      x           !   atom positions
!             real(kind=real64),dimension(3),intent(in)           ::      k0          !   centre of spot
!             complex(kind=real64)                    ::      FT   !   Fourier transform at discrete k points
!             
!             
!             integer                     ::      ii                  !   atom counters
!             real(kind=real64)           ::      x1,x2,x3            !   position of atom in 1,2,3 directions
!             complex(kind=real64)        ::      FT_part
!             real(kind=real64)           ::      kdotx
!                 
!             FT = 0.0d0             
!         
! !$OMP PARALLEL  PRIVATE(FT_part,x1,x2,x3,kdotx,ii),SHARED(FT)
!             FT_part = 0.0d0
!     !$OMP DO              
!             do ii = 1,size(x,dim=2)
! 
!                 x1 = x(1,ii)
!                 x2 = x(2,ii)
!                 x3 = x(3,ii)     !   extract position of atoms
!                 
!             !---    compute Ft at k0
!                 kdotx =  k0(1)*x1 + k0(2)*x2 + k0(3)*x3 
!                 FT_part = FT_part + exp( cmplx(0.0d0,kdotx,kind=real64) )               
!             end do
!     !$OMP END DO
!     
!     !$OMP CRITICAL
!             FT = FT + FT_part
!     !$OMP END CRITICAL      
! 
! !$OMP END PARALLEL        
!             
!             return
!         end function directFT
!             
    end module Lib_DiffractionSpot        