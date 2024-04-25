
    module Lib_FindLocalStrains
!---^^^^^^^^^^^^^^^^^^^^^^^^^^^
!*      self contained code to find local strains 
!*      from a link-cell list

        use Lib_LinkCell3d
        use Lib_SimpleProgressBar
        
        use Lib_SimpleSupercells
        use iso_fortran_env
        implicit none
        private
        
    !---
    
        public      ::      iterativelyComputeLocalStrains
        public      ::      iterativelyComputeLocalDisplacements
        public      ::      computeHomogeneousStrain
        public      ::      computeLocalStrainOnGrid
        public      ::      computeLocalDisplacementOnGrid
        public      ::      refineLocalStrainOnGrid
        public      ::      refineLocalDisplacementOnGrid
        
    !---
    
        logical,public              ::      FLS_DBG = .false.
    
        external    ::      DGESV
        external    ::      DSYSV
        
         
    contains
!---^^^^^^^^

        subroutine iterativelyComputeLocalStrains( lc3d,r0,a0,T,delta )
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      homogeneous strain computed and stored in element (0,0,0)
            type(LinkCell3d),intent(in)                              ::      lc3d
            real(kind=real64),dimension(:,:),intent(in)              ::      r0
            real(kind=real64),intent(in)                             ::      a0
            real(kind=real64),dimension(:,:,0:,0:,0:),intent(inout)  ::      T
            real(kind=real64),dimension(:,0:,0:,0:),intent(inout)    ::      delta
  

            integer         ::      Nx,Ny,Nz
            integer         ::      Mx,My,Mz,Lx,Ly,Lz
            integer         ::      ix,iy,iz
            real(kind=real64),dimension(:,:,:,:,:),allocatable  ::      T_tmp,uu
            real(kind=real64),dimension(:,:,:,:),allocatable    ::      delta_tmp
            real(kind=real64),dimension(3,3)        ::      Tsum,udotu
            real(kind=real64),dimension(3)          ::      aa
            real(kind=real64)           ::      ss
            
            integer         ::      step
            
        !---    find the size of the grids
            Nx = size(T,dim=3)
            Ny = size(T,dim=4)
            Nz = size(T,dim=5)
                
            allocate(T_tmp(3,3,0:Nx-1,0:Ny-1,0:Nz-1))
            allocate(delta_tmp(3,0:Nx-1,0:Ny-1,0:Nz-1))
            

            if (FLS_DBG) then
                write(*,fmt='(7(a,i6))') "Lib_FindLocalStrains::iterativelyComputeLocalStrains() info - grid ",Nx,"x",Ny,"x",Nz
            end if                                
        !---    start with 2x2x2 and end with powers of 2
            Mx = 1 ; My = 1 ; Mz = 1
            
         
            do
            !---    make a copy of the strains at previous resolution
                T_tmp = T ; delta_tmp = delta
            !---    interpolate a guess for the new grid
                Lx = min(4*Mx,Nx) ; Ly = min(4*Mx,Ny) ; Lz = min(4*Mx,Nz)
                call refineLocalStrainOnGrid( T_tmp(:,:,0:Mx-1,0:My-1,0:Mz-1),delta_tmp(:,0:Mx-1,0:My-1,0:Mz-1),    &
                                              T(:,:,0:Lx-1,0:Ly-1,0:Lz-1),delta(:,0:Lx-1,0:Ly-1,0:Lz-1) )
            !---    find strains at new resolution                                              
                call computeLocalStrainOnGrid( lc3d, r0,a0, T(:,:,0:Lx-1,0:Ly-1,0:Lz-1),delta(:,0:Lx-1,0:Ly-1,0:Lz-1),noSmooth = .false.  )                                              

                             
                                                              
            !---    all done?
                if (all( (/Lx-Nx,Ly-Ny,Lz-Nz/) == 0)) exit      !   done at full resolution                              
                Mx = Lx ; My = Ly ; Mz = Lz                                                                              
                
                
            end do
           
              return
            
            
            
            
            
        !--- now compute consistent delta
            
             !delta = 0.0d0          !   need to start with delta(1:3,0,0,0) = 0
              
            
        !---    T is computed as a best fit. But I want my displacements to be periodic
        !       since the displacement is defined by
        !       u(x) = r'(x) - r(x) = T(x) r(x) - r(x) = (T(x)-I) r(x)
        !       first remove the identity, and make the integral zero 
        !           T'(x) = T(x) - J
        !       such that int T'(x) dV = 0
        !           
            Tsum = 0.0d0
            do iz = 0,Nz-1
                do iy = 0,Ny-1
                    do ix = 0,Nx-1
                        Tsum(1:3,1:3) = Tsum(1:3,1:3) + T(1:3,1:3,ix,iy,iz) 
                    end do
                end do
            end do
            Tsum = Tsum/(Nx*Ny*Nz)
            do iz = 0,Nz-1
                do iy = 0,Ny-1
                    do ix = 0,Nx-1
                        T_tmp(1:3,1:3,ix,iy,iz) = T(1:3,1:3,ix,iy,iz) - Tsum(1:3,1:3)
                    end do
                end do
            end do
            
        !---    now each time I do a row or column of delta 
        !   I want to make sure that delta(Nx) = delta(0)
        !   I can do this by finding T" = T - K
        !   with  int T" dx = 0
        !   , or more cheaply I can work out the expected "extra" displacement  delta(Nx) - delta(0)        
        !   and subtract this from the row.
            aa(1) = getSuperCellSideLength( getSuper(lc3d),1 )/Nx
            aa(2) = getSuperCellSideLength( getSuper(lc3d),2 )/Ny
            aa(3) = getSuperCellSideLength( getSuper(lc3d),3 )/Nz
            
            print *,"integration steps ",aa
            
        !---    find the x,y and z lines     
            allocate(uu(1:3,0:Nx-1,0:Ny-1,0:Nz-1,3))
            
            
            do step = 1,10
            
                call  integrateStrain_array( T_tmp,delta,aa(1),1,uu(:,:,:,:,1) )
                call  integrateStrain_array( T_tmp,delta,aa(2),2,uu(:,:,:,:,2) )
                call  integrateStrain_array( T_tmp,delta,aa(3),3,uu(:,:,:,:,3) )
                
                print *,"minmax delta ",minval( delta(1,:,:,:) ),minval( delta(2,:,:,:) ),minval( delta(3,:,:,:) ),",",maxval( delta(1,:,:,:) ),maxval( delta(2,:,:,:) ),maxval( delta(3,:,:,:) )
                print *,"minmax u 1   ",minval( uu(1,:,:,:,1)  ),minval( uu(2,:,:,:,1)  ),minval( uu(3,:,:,:,1)  ),",",maxval( uu(1,:,:,:,1)  ),maxval( uu(2,:,:,:,1)  ),maxval( uu(3,:,:,:,1)  )
                print *,"minmax u 2   ",minval( uu(1,:,:,:,2)  ),minval( uu(2,:,:,:,2)  ),minval( uu(3,:,:,:,2)  ),",",maxval( uu(1,:,:,:,2)  ),maxval( uu(2,:,:,:,2)  ),maxval( uu(3,:,:,:,2)  )
                print *,"minmax u 3   ",minval( uu(1,:,:,:,3)  ),minval( uu(2,:,:,:,3)  ),minval( uu(3,:,:,:,3)  ),",",maxval( uu(1,:,:,:,3)  ),maxval( uu(2,:,:,:,3)  ),maxval( uu(3,:,:,:,3)  )
                
                udotu = 0.0d0
                do iz = 0,Nz-1
                    do iy = 0,Ny-1
                        do ix = 0,Nx-1 
                            udotu(1,1) = udotu(1,1) + uu(1,ix,iy,iz,1)*uu(1,ix,iy,iz,1) + uu(2,ix,iy,iz,1)*uu(2,ix,iy,iz,1) + uu(3,ix,iy,iz,1)*uu(3,ix,iy,iz,1)
                            udotu(2,1) = udotu(2,1) + uu(1,ix,iy,iz,2)*uu(1,ix,iy,iz,1) + uu(2,ix,iy,iz,2)*uu(2,ix,iy,iz,1) + uu(3,ix,iy,iz,2)*uu(3,ix,iy,iz,1)
                            udotu(3,1) = udotu(3,1) + uu(1,ix,iy,iz,3)*uu(1,ix,iy,iz,1) + uu(2,ix,iy,iz,3)*uu(2,ix,iy,iz,1) + uu(3,ix,iy,iz,3)*uu(3,ix,iy,iz,1)
                            udotu(2,2) = udotu(2,2) + uu(1,ix,iy,iz,2)*uu(1,ix,iy,iz,2) + uu(2,ix,iy,iz,2)*uu(2,ix,iy,iz,2) + uu(3,ix,iy,iz,2)*uu(3,ix,iy,iz,2)
                            udotu(3,2) = udotu(3,2) + uu(1,ix,iy,iz,3)*uu(1,ix,iy,iz,2) + uu(2,ix,iy,iz,3)*uu(2,ix,iy,iz,2) + uu(3,ix,iy,iz,3)*uu(3,ix,iy,iz,2)
                            udotu(3,3) = udotu(3,3) + uu(1,ix,iy,iz,3)*uu(1,ix,iy,iz,3) + uu(2,ix,iy,iz,3)*uu(2,ix,iy,iz,3) + uu(3,ix,iy,iz,3)*uu(3,ix,iy,iz,3)
                        end do
                    end do
                end do
                udotu(1,1) = sqrt( udotu(1,1) )
                udotu(2,2) = sqrt( udotu(2,2) )
                udotu(3,3) = sqrt( udotu(3,3) )
                print *,"|u1| = ",udotu(1,1)
                print *,"|u2| = ",udotu(2,2)
                print *,"|u3| = ",udotu(3,3)
                
                udotu(2,1) = udotu(2,1)/max(1.0d-16,udotu(1,1)*udotu(2,2))
                udotu(3,2) = udotu(3,2)/max(1.0d-16,udotu(2,2)*udotu(3,3))
                udotu(3,1) = udotu(3,1)/max(1.0d-16,udotu(1,1)*udotu(3,3))
                            
                print *,"dot_product u1.u2 = ",udotu(2,1) 
                print *,"dot_product u2.u3 = ",udotu(3,2) 
                print *,"dot_product u3.u1 = ",udotu(3,1) 
                
                
                do ix = 0,Nx-1 
                    write (*,fmt='(i4,100f12.5)') ix,delta(1:3,ix,1,1),uu(1:3,ix,1,1,1),uu(1:3,ix,iy,iz,2),uu(1:3,ix,1,1,3)               
                end do
                print *,""
                
                
                
                
                
                
                
                
                
                
                ss = (1 - maxval( abs( (/udotu(2,1),udotu(3,2),udotu(3,1)/) ) ))/6
               ! ss = 1.0d0/3
                do iz = 0,Nz-1        
                    do iy = 0,Ny-1    
                        do ix = 0,Nx-1
                            delta(1:3,ix,iy,iz) = delta(1:3,ix,iy,iz) + ss*( uu(1:3,ix,iy,iz,1) + uu(1:3,ix,iy,iz,2) + uu(1:3,ix,iy,iz,3) )
                        end do
                    end do
                end do
                   
            end do        
                    
!                
!             call integrateStrain_line( T_tmp(1:3,1,:,0,0),delta(1:3,:,0,0),aa(1) )        
!             call integrateStrain_line( T_tmp(1:3,2,0,:,0),delta(1:3,0,:,0),aa(2) )        
!             call integrateStrain_line( T_tmp(1:3,3,0,0,:),delta(1:3,0,0,:),aa(3) )        
!          
!              
!              do iy = 1,Ny-1
!                  do ix = 1,Nx-1
!                      delta(1:3,ix,iy,0) = ( delta(1:3,ix-1,iy,0) + delta(1:3,ix,iy-1,0) )/2      &
!                              + ( T_tmp(1:3,1,ix-1,iy,0) + T_tmp(1:3,1,ix,iy,0) + T_tmp(1:3,2,ix,iy-1,0) + T_tmp(1:3,2,ix,iy,0)  )/4 
!                  end do
!              end do
!              
!              do iz = 1,Nz-1
!                  do iy = 1,Ny-1
!                      delta(1:3,0,iy,iz) = ( delta(1:3,0,iy-1,iz) + delta(1:3,0,iy,iz-1) )/2      &
!                              + ( T_tmp(1:3,2,0,iy-1,iz) + T_tmp(1:3,2,0,iy,iz) + T_tmp(1:3,3,0,iy,iz-1) + T_tmp(1:3,3,0,iy,iz))/4 
!                      !print *,iy,iz,delta(1:3,0,iy,iz)
!                  end do
!              end do
!              
!              do iz = 1,Nz-1
!                  do ix = 1,Nx-1
!                      delta(1:3,ix,0,iz) = ( delta(1:3,ix-1,0,iz) + delta(1:3,ix,0,iz-1) )/2      &
!                              + ( T_tmp(1:3,1,ix-1,0,iz) + T_tmp(1:3,1,ix,0,iz) + T_tmp(1:3,3,ix,0,iz-1) + T_tmp(1:3,3,ix,0,iz)  )/4
!                  end do
!              end do
!              
!              
!              do iz = 1,Nz-1
!                  do iy = 1,Ny-1
!                      do ix = 1,Nx-1
!                          delta(1:3,ix,iy,iz) = ( delta(1:3,ix-1,iy,iz) + delta(1:3,ix,iy-1,iz) + delta(1:3,ix,iy,iz-1) )/3    &
!                                      + ( T_tmp(1:3,1,ix-1,iy,iz) + T_tmp(1:3,1,ix,iy,iz) + T_tmp(1:3,2,ix,iy-1,iz) + T_tmp(1:3,2,ix,iy,iz) + T_tmp(1:3,3,ix,iy,iz-1) + T_tmp(1:3,3,ix,iy,iz)  )/6
!                      end do
!                  end do
!              end do
!              
!              dsum = 0.0d0
!              do iz = 0,Nz-1
!                  do iy = 0,Ny-1
!                      do ix = 0,Nx-1
!                         dsum = dsum + delta(1:3,ix,iy,iz)
!                      end do
!                  end do
!              end do
!              dsum = dsum/(Nx*Ny*Nz)
!              do iz = 0,Nz-1
!                  do iy = 0,Ny-1
!                      do ix = 0,Nx-1
!                         delta(1:3,ix,iy,iz) = delta(1:3,ix,iy,iz) - dsum 
!                      end do
!                  end do
!              end do
!              
!              
!              delta = delta * a0
!              
!              print *,"self-consistent delta"
!              print *,minval(delta(1,:,:,:)),maxval(delta(1,:,:,:)),sum(delta(1,:,:,:))/(Nx*Ny*Nz)
!              print *,minval(delta(2,:,:,:)),maxval(delta(2,:,:,:)),sum(delta(2,:,:,:))/(Nx*Ny*Nz)
!              print *,minval(delta(3,:,:,:)),maxval(delta(3,:,:,:)),sum(delta(3,:,:,:))/(Nx*Ny*Nz)
!              
!             
!             
!             
!             
            
            return
            
        contains
    !---^^^^^^^^    
            
            subroutine integrateStrain_line( T,delta,dx,u )            
        !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
        !*      given the array of values of strain T(1:3,1:N)
        !*      update the displacements delta(1:3,1:N)
        !*      this is done by solving for T(:,ii) = ( delta(:,ii+1) - delta(:,ii-1) )/(2 dx)
        !*      note that delta is inout
                real(kind=real64),dimension(:,:),intent(in)        ::      T            !   (3,nn)  taken from (1:3,1,1:nn,ix,iy)
                real(kind=real64),dimension(:,:),intent(in)        ::      delta        !   (3,nn)  taken from (1:3,1:nn,ix,iy)
                real(kind=real64),intent(in)                       ::      dx
                real(kind=real64),dimension(:,:),intent(out)       ::      u            !   (3,nn)  taken from (1:3,1:nn,ix,iy)
                
                real(kind=real64)            ::      a1,a2,a3,v1,v2,v3
                integer     ::     ii,nn
                                
                nn = size(T,dim=2)                
            !---    forward predict u in a line
                u(1:3,1) = 0.0d0
                do ii = 2,nn
                    u(1:3,ii) = ( delta(:,ii-1) - delta(:,ii) + u(:,ii-1)  ) + ( T(:,ii-1) + T(:,ii) )*dx/2
                end do
                
                if (FLS_DBG) then
                print *,""
                    print *,"integrateStrain_line"
                    print *,""
                do ii = 1,nn
                    write(*,fmt='(i4,100f12.4)') ii,delta(1:3,ii),(delta(1:3,mod(ii+nn+1,nn)+1)-delta(1:3,mod(ii+nn-1,nn)+1))/(2*dx),T(1:3,ii)
                end do
                end if
                
                
                
            !---    make a guess for delta(n+1) + u(n+1)
                v1 = ( delta(1,nn) + u(1,nn)  ) + ( T(1,nn) + T(1,1) )*dx/2
                v2 = ( delta(2,nn) + u(2,nn)  ) + ( T(2,nn) + T(2,1) )*dx/2
                v3 = ( delta(3,nn) + u(3,nn)  ) + ( T(3,nn) + T(3,1) )*dx/2
                
            !---    ... but it could well be that u(i) ~ a i dx. Find the linear term
                a1 = ( v1 - u(1,1) )/(nn*dx)
                a2 = ( v2 - u(2,1) )/(nn*dx)
                a3 = ( v3 - u(3,1) )/(nn*dx)
                
!                if (FLS_DBG) write(*,fmt='(a4,100f12.4)') "v,a",v1,v2,v3,a1,a2,a3
                
            !---    remove linear term                
                do ii = 1,nn
                    u(1,ii) = u(1,ii) - ii*a1*dx
                    u(2,ii) = u(2,ii) - ii*a2*dx
                    u(3,ii) = u(3,ii) - ii*a3*dx                    
                end do
                
!                 if (FLS_DBG) then
!                     do ii = 1,nn
!                         write(*,fmt='(i4,100f12.4)') ii,u(1:3,ii)
!                     end do                
!                 end if
            !---    compute constant term
                v1 = ( sum(delta(1,1:nn))+sum(u(1,1:nn)) ) / nn
                v2 = ( sum(delta(2,1:nn))+sum(u(2,1:nn)) ) / nn
                v3 = ( sum(delta(3,1:nn))+sum(u(3,1:nn)) ) / nn
                
!                if (FLS_DBG) print *, "v",v1,v2,v3 
                
            !---    remove constant term                
                do ii = 1,nn
                    u(1,ii) = u(1,ii) - v1
                    u(2,ii) = u(2,ii) - v2
                    u(3,ii) = u(3,ii) - v3                  
                end do            
                
                if (FLS_DBG) then
                    
                    do ii = 1,nn
                        write(*,fmt='(i4,100f12.4)') ii,u(1:3,ii),(u(1:3,mod(ii+nn+1,nn)+1)-u(1:3,mod(ii+nn-1,nn)+1))/(2*dx),T(1:3,ii)
                    end do    
                    print *,""
                    print *,""
                    print *,""
                end if            
                !stop
                
                
              ! 
              ! off1 = u(1,nn) + delta(1,nn) - delta(1,1)
              ! off2 = u(2,nn) + delta(2,nn) - delta(2,1)
              ! off3 = u(3,nn) + delta(3,nn) - delta(3,1)
              ! off1 = off1/nn
              ! off2 = off2/nn
              ! off3 = off3/nn
              ! u(1,1:nn) = u(1,1:nn) - off1
              ! u(2,1:nn) = u(2,1:nn) - off2
              ! u(3,1:nn) = u(3,1:nn) - off3
                
              ! if (any(abs(u(1,:))>1)) then
              !     print *,"integrateStrain_line"
              !      nn = size(T,dim=2)                
              !         u(1:3,1) = 0.0d0
              !         do ii = 2,nn
              !             u(1:3,ii) = ( delta(:,ii-1) + u(:,ii-1) - delta(:,ii) ) + ( T(:,ii-1) + T(:,ii) )*dx/2
              !             write(*,fmt='(i4,100f12.4)') ii,( delta(:,ii-1) - delta(:,ii) ),( T(:,ii-1) + T(:,ii) )*dx/2,u(1:3,ii)
              !         end do
              !         a1 = u(1,nn) + delta(1,nn) - delta(1,1)
              !         a2 = u(2,nn) + delta(2,nn) - delta(2,1)
              !         a3 = u(3,nn) + delta(3,nn) - delta(3,1)
              !         a1 = off1/nn
              !         a2 = off2/nn
              !         a3 = off3/nn
              !         write(*,fmt='(a4,100f12.4)') "off ",off1,off2,off3
              !         u(1,1:nn) = u(1,1:nn) - off1
              !         u(2,1:nn) = u(2,1:nn) - off2
              !         u(3,1:nn) = u(3,1:nn) - off3
              !         stop
              !   end if
                
                
                return
            end subroutine integrateStrain_line
            
              
            subroutine integrateStrain_array( T,delta,dx,d,u )            
        !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
        !*      given the array of values of strain T(1:3,1:N,:,:)
        !*      update the displacements delta(1:3,1:N)
        !*      this is done by solving for T(:,ii) = ( delta(:,ii+1) - delta(:,ii-1) + u(:,ii+1) - u(:,ii-1) )/(2 dx)
        
                real(kind=real64),dimension(:,:,:,:,:),intent(in)   ::      T           !   (3,3,Nx,Ny,Nz)
                real(kind=real64),dimension(:,:,:,:),intent(in)     ::      delta       !   (3,Nx,Ny,Nz)
                real(kind=real64),dimension(:,:,:,:),intent(out)    ::      u           !   (3,Nx,Ny,Nz)
                integer,intent(in)                                  ::      d
                real(kind=real64),intent(in)                        ::      dx
                
                real(kind=real64),dimension(:,:),allocatable        ::      aa
                real(kind=real64),dimension(:,:),allocatable        ::      bb
                integer,dimension(:),allocatable                    ::      ipiv
                integer     ::      ii,nn
                integer     ::      n1,n2,i1,i2,mm
                
                real(kind=real64),dimension(:,:),allocatable   ::      Ttmp           !   (3,nn)
                real(kind=real64),dimension(3)                 ::      Tbar           !    
                
                                
            !---    determine integration direction and orthogonal array directions                
                select case(d) 
                    case(1)
                        nn = size(T,dim=3)
                        n1 = size(T,dim=4)
                        n2 = size(T,dim=5)                                
                    case(2)
                        nn = size(T,dim=4)
                        n1 = size(T,dim=3)
                        n2 = size(T,dim=5)
                    case(3)
                        nn = size(T,dim=5)
                        n1 = size(T,dim=3)
                        n2 = size(T,dim=4)
                end select
                
                
                
                
               
              
              
               select case(d) 
                  case(1)
                      do i2 = 1,n2
                          do i1 = 1,n1
                            FLS_DBG = ( (i1-1)*(i1-1) + (i2-1)*(i2-1) )==0
                              call integrateStrain_line( T(1:3,1,1:nn,i1,i2),delta(1:3,1:nn,i1,i2),dx,u(1:3,1:nn,i1,i2) )            
                          end do
                      end do
                              
                  case(2)
                      do i2 = 1,n2
                          do i1 = 1,n1
                              call integrateStrain_line( T(1:3,2,i1,1:nn,i2),delta(1:3,i1,1:nn,i2),dx,u(1:3,i1,1:nn,i2) )            
                          end do                                                                                      
                      end do                                                                                          
                                                                                                                      
                  case(3)                                                                                             
                      do i2 = 1,n2                                                                                    
                          do i1 = 1,n1                                                                                
                              call integrateStrain_line( T(1:3,3,i1,i2,1:nn),delta(1:3,i1,i2,1:nn),dx,u(1:3,i1,i2,1:nn) )            
                          end do
                      end do
                              
              end select
              
              return
               
                
                
                
                
                allocate(Ttmp(3,nn)) 
                
                
                
                
                
                
                
            !---    define tridiagonal (+pbc) matrix (-1 0 1)
            
                 
                allocate( aa(nn+1,nn+1) )      ; aa = 0  
                do ii = 1,nn
                    aa(ii, mod(ii+nn-2,nn)+1 ) = 1.0d0/6
                    aa(ii, mod(ii+nn-1,nn)+1 ) = -8.0d0/6
                    aa(ii, mod(ii+nn+1,nn)+1 ) = 8.0d0/6
                    aa(ii, mod(ii+nn+2,nn)+1 ) = -1.0d0/6                       
                end do
                
                 
                aa(1:nn,nn+1) = 1
                aa(nn+1,1:nn) = 1
               
                
                allocate(ipiv(nn+1))
                
               ! 
               !
               ! call DGETRF(nn+1, nn+1, aa, nn+1, ipiv, ii)
               ! print *,"DGETRF returns ",ii
               ! 
               ! stop
                
                
            !---    define right hand sides
                print *,"integrateStrain_array d,nn,n1,n2 ",d,nn,n1,n2
                allocate( bb(nn+1,3*n1*n2) ) ; bb = 0
                select case(d) 
                    case(1)
                        do i2 = 1,n2
                            do i1 = 1,n1
                            
                               Ttmp(1:3,1:nn) =  T(1:3,1,1:nn,i1,i2)
                               Tbar(1) = sum(Ttmp(1,1:nn))/nn
                               Tbar(2) = sum(Ttmp(2,1:nn))/nn
                               Tbar(3) = sum(Ttmp(3,1:nn))/nn
                               Ttmp(1,1:nn) = Ttmp(1,1:nn) - Tbar(1)
                               Ttmp(2,1:nn) = Ttmp(2,1:nn) - Tbar(2)
                               Ttmp(3,1:nn) = Ttmp(3,1:nn) - Tbar(3)
                            
                               mm = 3*(i1-1) + (3*n1)*(i2-1) 
                                bb(1,mm+1:mm+3) = 2*dx* Ttmp(1:3,1) - (delta(1:3,2,i1,i2)-delta(1:3,nn,i1,i2)) 
                                do ii = 2,nn-1
                                    bb(ii,mm+1:mm+3) = 2*dx* Ttmp(1:3,ii) - (delta(1:3,ii+1,i1,i2)-delta(1:3,ii-1,i1,i2)) 
                                end do
                                bb(nn,mm+1:mm+3) = 2*dx* Ttmp(1:3,nn) - (delta(1:3,1,i1,i2)-delta(1:3,nn-1,i1,i2)) 
                                                               
                            
                               ! mm = 3*(i1-1) + (3*n1)*(i2-1) 
                               ! bb(1,mm+1:mm+3) = 2*dx* T(1:3,1,1,i1,i2) - (delta(1:3,2,i1,i2)-delta(1:3,nn,i1,i2)) 
                               ! do ii = 2,nn-1
                               !     bb(ii,mm+1:mm+3) = 2*dx* T(1:3,1,ii,i1,i2) - (delta(1:3,ii+1,i1,i2)-delta(1:3,ii-1,i1,i2)) 
                               ! end do
                               ! bb(nn,mm+1:mm+3) = 2*dx* T(1:3,1,nn,i1,i2) - (delta(1:3,1,i1,i2)-delta(1:3,nn-1,i1,i2)) 
                            end do
                        end do    
                    case(2)
                        do i2 = 1,n2
                            do i1 = 1,n1
                                mm = 3*(i1-1) + (3*n1)*(i2-1) 
                                bb(1,mm+1:mm+3) = 2*dx* T(1:3,2,i1,1,i2) - (delta(1:3,i1,2,i2)-delta(1:3,i1,nn,i2)) 
                                do ii = 2,nn-1
                                    bb(ii,mm+1:mm+3) = 2*dx* T(1:3,2,i1,ii,i2) - (delta(1:3,i1,ii+1,i2)-delta(1:3,i1,ii-1,i2)) 
                                end do
                                bb(nn,mm+1:mm+3) = 2*dx* T(1:3,2,i1,nn,i2) - (delta(1:3,i1,1,i2)-delta(1:3,i1,nn-1,i2))
                            end do
                        end do    
                    case(3)
                        do i2 = 1,n2
                            do i1 = 1,n1
                                mm = 3*(i1-1) + (3*n1)*(i2-1) 
                                bb(1,mm+1:mm+3) = 2*dx* T(1:3,3,i1,i2,1) - (delta(1:3,i1,i2,2)-delta(1:3,i1,i2,nn)) 
                                do ii = 2,nn-1
                                    bb(ii,mm+1:mm+3) = 2*dx* T(1:3,3,i1,i2,ii) - (delta(1:3,i1,i2,ii+1)-delta(1:3,i1,i2,ii-1)) 
                                end do
                                bb(nn,mm+1:mm+3) = 2*dx* T(1:3,3,i1,i2,nn) - (delta(1:3,i1,i2,1)-delta(1:3,i1,i2,nn-1)) 
                            end do
                        end do    
                end select
                
               ! do ii = 1,nn+1
               !     write (*,fmt='(1000f8.3)',advance="no") aa(ii,:)
               !     write (*,fmt='(a,1000f8.3)',advance="yes") " : ",bb(ii,1:3)
               ! end do
                
               ! do ii = 1,nn+1
               !     write (*,fmt='(3f6.2,a,100f6.2)') bb(ii,1:3),",",aa(ii,:)
               ! end do
                
!                call DGESV( nn+1,3*n1*n2,aa,nn+1,ipiv,bb,nn+1,ii )      !   note NRHS = 3, there are 3 independent right hand sides for each displacment.
                call DGESV( nn+1,3*n1*n2,aa,nn+1,ipiv,bb,nn+1,ii )      !   note NRHS = 3, there are 3 independent right hand sides for each displacment.
                
               ! print *,"DGESV returns ",ii
                
                
                
                if (ii/=0) stop "DGESV fail"
                
              !  do ii = 1,nn+1                        
              !      write (*,fmt='(3f6.2,a,100f6.2)') bb(ii,1:3),",",aa(ii,:)
              !  end do                                
                                                      
                
                
                
                
                select case(d) 
                    case(1)
                        do i2 = 1,n2
                            do i1 = 1,n1
                                mm = 3*(i1-1) + (3*n1)*(i2-1) 
                                do ii = 1,nn
                                    u(1:3,ii,i1,i2) = bb(ii,mm+1:mm+3)
                                end do
                                !print *,"d=1 ",i1,i2,sum(u(1,:,i1,i2)),sum(u(2,:,i1,i2)),sum(u(3,:,i1,i2)),maxval(abs(u(1,:,i1,i2))),maxval(abs(u(2,:,i1,i2))),maxval(abs(u(3,:,i1,i2)))
                            end do
                        end do    
                    case(2)
                        do i2 = 1,n2
                            do i1 = 1,n1
                                mm = 3*(i1-1) + (3*n1)*(i2-1) 
                                do ii = 1,nn
                                    u(1:3,i1,ii,i2) = bb(ii,mm+1:mm+3)
                                end do
                            end do
                        end do    
                    case(3)
                        do i2 = 1,n2
                            do i1 = 1,n1
                                mm = 3*(i1-1) + (3*n1)*(i2-1) 
                                do ii = 1,nn
                                    u(1:3,i1,i2,ii) = bb(ii,mm+1:mm+3)
                                end do
                            end do
                        end do    
                end select
                
                
                
                do ii = 1,nn
                    write(*,fmt='(i4,100f12.4)') ii,u(1:3,ii,1,1),(u(1:3,mod(ii+nn+1,nn)+1,1,1)-u(1:3,mod(ii+nn-1,nn)+1,1,1))/(2*dx),T(1:3,1,ii,1,1)
                end do                
                stop
                
                
                
                
                
                
                
                
                
                !stop
                return
            end subroutine integrateStrain_array
            
        end subroutine iterativelyComputeLocalStrains

        subroutine iterativelyComputeLocalDisplacements( lc3d,r0,a0,delta )
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      homogeneous strain computed and stored in element (0,0,0)
            type(LinkCell3d),intent(in)                              ::      lc3d
            real(kind=real64),dimension(:,:),intent(in)              ::      r0
            real(kind=real64),intent(in)                             ::      a0
            real(kind=real64),dimension(:,0:,0:,0:),intent(inout)    ::      delta

            integer         ::      Nx,Ny,Nz
            integer         ::      Mx,My,Mz,Lx,Ly,Lz

            real(kind=real64),dimension(:,:,:,:),allocatable    ::      delta_tmp
            
            
        !---    find the size of the grids
            Nx = size(delta,dim=2)
            Ny = size(delta,dim=3)
            Nz = size(delta,dim=4)
                            
            allocate(delta_tmp(3,0:Nx-1,0:Ny-1,0:Nz-1))
            

            if (FLS_DBG) then
                write(*,fmt='(7(a,i6))') "Lib_FindLocalStrains::iterativelyComputeLocalDisplacements() info - grid ",Nx,"x",Ny,"x",Nz
            end if                                
        !---    start with 2x2x2 and end with powers of 2
            Mx = 1 ; My = 1 ; Mz = 1
            
        
            do
            
            !---    make a copy of the strains at previous resolution
                delta_tmp = delta
            !---    interpolate a guess for the new grid
                Lx = min(4*Mx,Nx) ; Ly = min(4*Mx,Ny) ; Lz = min(4*Mx,Nz)
                call refineLocalDisplacementOnGrid( delta_tmp(:,0:Mx-1,0:My-1,0:Mz-1),delta(:,0:Lx-1,0:Ly-1,0:Lz-1) )
            !---    find strains at new resolution                                              
                call computeLocalDisplacementOnGrid( lc3d, r0,a0, delta(:,0:Lx-1,0:Ly-1,0:Lz-1) )                                              

                             
                                                              
            !---    all done?
                if (all( (/Lx-Nx,Ly-Ny,Lz-Nz/) == 0)) exit      !   done at full resolution                              
                Mx = Lx ; My = Ly ; Mz = Lz                                                                              
                
                
            end do
            
            return
        end subroutine iterativelyComputeLocalDisplacements

        subroutine computeHomogeneousStrain( lc3d, r0,a0, T,delta )
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
            type(LinkCell3d),intent(in)                     ::      lc3d
            real(kind=real64),dimension(:,:),intent(in)     ::      r0
            real(kind=real64),intent(in)                    ::      a0
            real(kind=real64),dimension(3,3),intent(out)    ::      T
            real(kind=real64),dimension(3),intent(out)      ::      delta

            real(kind=real64),dimension(3,3,0:0,0:0,0:0)    ::      TT
            real(kind=real64),dimension(3,0:0,0:0,0:0)      ::      dd
    
            TT(1:3,1:3,0,0,0) = reshape( (/1,0,0, 0,1,0, 0,0,1/) , (/3,3/) )
            dd(1:3,0,0,0) = 0
            
            call computeLocalStrainOnGrid( lc3d, r0,a0, TT,dd,noSmooth = .false. )
            
            T(1:3,1:3) = TT(1:3,1:3,0,0,0)
            delta(1:3) = dd(1:3,0,0,0)
            
            if (FLS_DBG) then
                print *,"Lib_FindLocalStrains::computeHomogeneousStrain() info - "
                write(*,fmt='(3f14.8,a4,f14.8)') T(1,1:3),"   ",delta(1)
                write(*,fmt='(3f14.8,a4,f14.8)') T(2,1:3),"   ",delta(2)
                write(*,fmt='(3f14.8,a4,f14.8)') T(3,1:3),"   ",delta(3)
            end if
            
            return
        end subroutine computeHomogeneousStrain
            
            
            
        subroutine refineLocalStrainOnGrid( T_in,delta_in, T_out,delta_out )
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      have strain computed on a grid.
    !*      now interpolate onto a new grid
    
            real(kind=real64),dimension(:,:,0:,0:,0:),intent(inout)  ::      T_in
            real(kind=real64),dimension(:,0:,0:,0:),intent(inout)    ::      delta_in
            
            real(kind=real64),dimension(:,:,0:,0:,0:),intent(inout)  ::      T_out
            real(kind=real64),dimension(:,0:,0:,0:),intent(inout)    ::      delta_out
            
            integer         ::      Nx,Ny,Nz
            integer         ::      Mx,My,Mz
            
            integer             ::      ix,iy,iz , jx,jy,jz
            integer             ::      kx,ky,kz
            real(kind=real64)   ::      ux,uy,uz
            
            real(kind=real64),dimension(3,3)        ::      TT
            real(kind=real64),dimension(3)          ::      dd
                        
            
        !---    find the size of the grids
            Nx = size(T_in,dim=3)
            Ny = size(T_in,dim=4)
            Nz = size(T_in,dim=5)
            Mx = size(T_out,dim=3)
            My = size(T_out,dim=4)
            Mz = size(T_out,dim=5)

            if (FLS_DBG) then
                write(*,fmt='(7(a,i6))') "Lib_FindLocalStrains::refineLocalStrainOnGrid() info - refine ",Nx,"x",Ny,"x",Nz," to ",Mx,"x",My,"x",Mz
            end if
                        
        !---    quick exit?
            if (maxval( (/Nx,Ny,Nz/) )==1) then
                !   single grid point on input lattice
                TT = T_in(1:3,1:3,0,0,0)
                dd = delta_in(1:3,0,0,0)
                do kz = 0,Mz-1
                    do ky = 0,My-1
                        do kx = 0,Mx-1
                            T_out(1:3,1:3,kx,ky,kz) = TT
                            delta_out(1:3,kx,ky,kz) = dd
                        end do
                    end do
                end do
                return
            end if            
            
        !---    loop through each point on the output grid
            do kz = 0,Mz-1
                do ky = 0,My-1
                    do kx = 0,Mx-1  
                        !   position of this point is (kx+1/2,ky+1/2,kz+1/2)
                        
                        !   which half grid point is this nearest to in input lattice?
                        !   want (ix+1/2)/Nx <= (kx+1/2)/Mx <= (ix+3/2)/Nx
                         
                        
                        ux = Nx*real(kx,kind=real64)/Mx                       !   cell position of output grid point in input grid
                        uy = Ny*real(ky,kind=real64)/My
                        uz = Nz*real(kz,kind=real64)/Mz
                        
                        ix = floor( ux ) 
                        iy = floor( uy ) 
                        iz = floor( uz ) 
            
                        ux = ux - ix          !   distance from input grid = 0 if aligned LHS and =1 if aligned RHS
                        uy = uy - iy
                        uz = uz - iz
                        
                        ix = mod( ix+Nx,Nx ) ; jx = mod( ix+Nx+1,Nx )
                        iy = mod( iy+Ny,Ny ) ; jy = mod( iy+Ny+1,Ny )
                        iz = mod( iz+Nz,Nz ) ; jz = mod( iz+Nz+1,Nz )
                        
                       ! write(*,fmt='(a,3i4,a,3i4,a,3i4,a,3f12.5)')  "refine ",kx,ky,kz," from ",ix,iy,iz," : ",jx,jy,jz," u ",ux,uy,uz
                              
                        TT = T_in(1:3,1:3,ix,iy,iz)*(1-ux)*(1-uy)*(1-uz)         & 
                           + T_in(1:3,1:3,jx,iy,iz)*(  ux)*(1-uy)*(1-uz)         &
                           + T_in(1:3,1:3,ix,jy,iz)*(1-ux)*(  uy)*(1-uz)         &
                           + T_in(1:3,1:3,jx,jy,iz)*(  ux)*(  uy)*(1-uz)         &
                           + T_in(1:3,1:3,ix,iy,jz)*(1-ux)*(1-uy)*(  uz)         &
                           + T_in(1:3,1:3,jx,iy,jz)*(  ux)*(1-uy)*(  uz)         &
                           + T_in(1:3,1:3,ix,jy,jz)*(1-ux)*(  uy)*(  uz)         &
                           + T_in(1:3,1:3,jx,jy,jz)*(  ux)*(  uy)*(  uz) 
                                                   
                        dd = delta_in(1:3,ix,iy,iz)*(1-ux)*(1-uy)*(1-uz)         & 
                           + delta_in(1:3,jx,iy,iz)*(  ux)*(1-uy)*(1-uz)         &
                           + delta_in(1:3,ix,jy,iz)*(1-ux)*(  uy)*(1-uz)         &
                           + delta_in(1:3,jx,jy,iz)*(  ux)*(  uy)*(1-uz)         &
                           + delta_in(1:3,ix,iy,jz)*(1-ux)*(1-uy)*(  uz)         &
                           + delta_in(1:3,jx,iy,jz)*(  ux)*(1-uy)*(  uz)         &
                           + delta_in(1:3,ix,jy,jz)*(1-ux)*(  uy)*(  uz)         &
                           + delta_in(1:3,jx,jy,jz)*(  ux)*(  uy)*(  uz) 
                           
                        T_out(1:3,1:3,kx,ky,kz) = TT
                        delta_out(1:3,kx,ky,kz) = dd
                        
                    end do
                end do
            end do
            
           ! stop
            return
        end subroutine refineLocalStrainOnGrid
                        
         
            
        subroutine refineLocalDisplacementOnGrid( delta_in, delta_out )
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ 
    !*      have strain computed on a grid.
    !*      now interpolate onto a new grid    
            real(kind=real64),dimension(:,0:,0:,0:),intent(inout)    ::      delta_in            
            real(kind=real64),dimension(:,0:,0:,0:),intent(inout)    ::      delta_out
            
            integer         ::      Nx,Ny,Nz
            integer         ::      Mx,My,Mz
            
            integer             ::      ix,iy,iz , jx,jy,jz
            integer             ::      kx,ky,kz
            real(kind=real64)   ::      ux,uy,uz
            
            
            real(kind=real64),dimension(3)          ::      dd
                        
            
        !---    find the size of the grids
            Nx = size(delta_in,dim=2)
            Ny = size(delta_in,dim=3)
            Nz = size(delta_in,dim=4)
            Mx = size(delta_out,dim=2)
            My = size(delta_out,dim=3)
            Mz = size(delta_out,dim=4)

            if (FLS_DBG) then
                write(*,fmt='(7(a,i6))') "Lib_FindLocalStrains::refineLocalDisplacementOnGrid() info - refine ",Nx,"x",Ny,"x",Nz," to ",Mx,"x",My,"x",Mz
            end if
                        
        !---    quick exit?
            if (maxval( (/Nx,Ny,Nz/) )==1) then
                !   single grid point on input lattice
                dd = delta_in(1:3,0,0,0)
                do kz = 0,Mz-1
                    do ky = 0,My-1
                        do kx = 0,Mx-1
                            delta_out(1:3,kx,ky,kz) = dd
                        end do
                    end do
                end do
                return
            end if            
            
        !---    loop through each point on the output grid
            do kz = 0,Mz-1
                do ky = 0,My-1
                    do kx = 0,Mx-1  
                        !   position of this point is (kx+1/2,ky+1/2,kz+1/2)
                        
                        !   which half grid point is this nearest to in input lattice?
                        
                        ux = Nx*real(kx,kind=real64)/Mx                       !   cell position of output grid point in input grid
                        uy = Ny*real(ky,kind=real64)/My
                        uz = Nz*real(kz,kind=real64)/Mz
                        
                        ix = floor( ux ) 
                        iy = floor( uy ) 
                        iz = floor( uz ) 
            
                        ux = ux - ix          !   distance from input grid = 0 if aligned LHS and =1 if aligned RHS
                        uy = uy - iy
                        uz = uz - iz
                        
                        ix = mod( ix+Nx,Nx ) ; jx = mod( ix+Nx+1,Nx )
                        iy = mod( iy+Ny,Ny ) ; jy = mod( iy+Ny+1,Ny )
                        iz = mod( iz+Nz,Nz ) ; jz = mod( iz+Nz+1,Nz )
                        
                       ! write(*,fmt='(a,3i4,a,3i4,a,3i4,a,3f12.5)')  "refine ",kx,ky,kz," from ",ix,iy,iz," : ",jx,jy,jz," u ",ux,uy,uz
                                                 
                        dd = delta_in(1:3,ix,iy,iz)*(1-ux)*(1-uy)*(1-uz)         & 
                           + delta_in(1:3,jx,iy,iz)*(  ux)*(1-uy)*(1-uz)         &
                           + delta_in(1:3,ix,jy,iz)*(1-ux)*(  uy)*(1-uz)         &
                           + delta_in(1:3,jx,jy,iz)*(  ux)*(  uy)*(1-uz)         &
                           + delta_in(1:3,ix,iy,jz)*(1-ux)*(1-uy)*(  uz)         &
                           + delta_in(1:3,jx,iy,jz)*(  ux)*(1-uy)*(  uz)         &
                           + delta_in(1:3,ix,jy,jz)*(1-ux)*(  uy)*(  uz)         &
                           + delta_in(1:3,jx,jy,jz)*(  ux)*(  uy)*(  uz) 
                           
                      
                        delta_out(1:3,kx,ky,kz) = dd
                        
                    end do
                end do
            end do
            
           ! stop
            return
        end subroutine refineLocalDisplacementOnGrid
                        

        subroutine computeLocalStrainOnGrid( lc3d, r0,a0, T,delta , nosmooth )
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      compute the strain at evenly spaced points, within a width sigma
    !*      given the expected separations r0
    !*      on input, T, delta are expected roto-strain transformation and displacement
    !*      on exit, T, delta are better fitted 
    !*          sense:  r = T r0 + delta
    !*      values are computed at helf-grid points 1/2,1/2,1/2  not 0,0,0
    
            type(LinkCell3d),intent(in)                     ::      lc3d
            real(kind=real64),dimension(:,:),intent(in)     ::      r0
            real(kind=real64),intent(in)        ::      a0   
            real(kind=real64),dimension(:,:,0:,0:,0:),intent(inout)  ::      T
            real(kind=real64),dimension(:,0:,0:,0:),intent(inout)    ::      delta
            logical,intent(in)                                      ::      nosmooth
            
            real(kind=real64),dimension(3,3)    ::      aa,TT,T_lcl
            real(kind=real64),dimension(:,:),pointer        ::      dx_all,dx_all_tmp
            real(kind=real64),dimension(:,:),allocatable    ::      dx,pos
            integer,dimension(:),allocatable                ::      id
            real(kind=real64),dimension(:),pointer          ::      weight,weight_tmp
            real(kind=real64),dimension(:),allocatable      ::      dr2
            real(kind=real64),dimension(3,size(r0,dim=2))   ::      r0_lcl
            type(SimpleSupercell)           ::      super
            integer             ::      Nx,Ny,Nz, ix,iy,iz,ik, ii  !,mx,my,mz
            real(kind=real64)   ::      sigma,i2s2,ww
            real(kind=real64),dimension(3)  ::      xx,yy,dd,d_lcl
            integer             ::      nAtoms,nn,nn_all,nn2,nFails,nR0
            logical             ::      ok
            
        !---    find the size of the grid
            Nx = size(T,dim=3)
            Ny = size(T,dim=4)
            Nz = size(T,dim=5)

            
            
              !nn = Nx*Ny*Nz
              ! 
              !print *,"Lib_FindLocalStrains::computeLocalStrainOnGrid() dbg - input"
              !print *,"Lib_FindLocalStrains::computeLocalStrainOnGrid() info - maximum values"
              !write(*,fmt='(3f14.8,a4,f14.8)') maxval(T(1,1,:,:,:)),maxval(T(1,2,:,:,:)),maxval(T(1,3,:,:,:)),"   ",maxval(delta(1,:,:,:))
              !write(*,fmt='(3f14.8,a4,f14.8)') maxval(T(2,1,:,:,:)),maxval(T(2,2,:,:,:)),maxval(T(2,3,:,:,:)),"   ",maxval(delta(2,:,:,:))
              !write(*,fmt='(3f14.8,a4,f14.8)') maxval(T(3,1,:,:,:)),maxval(T(3,2,:,:,:)),maxval(T(3,3,:,:,:)),"   ",maxval(delta(3,:,:,:))
              !print *,"Lib_FindLocalStrains::computeLocalStrainOnGrid() info - minimum values"
              !write(*,fmt='(3f14.8,a4,f14.8)') minval(T(1,1,:,:,:)),minval(T(1,2,:,:,:)),minval(T(1,3,:,:,:)),"   ",minval(delta(1,:,:,:))
              !write(*,fmt='(3f14.8,a4,f14.8)') minval(T(2,1,:,:,:)),minval(T(2,2,:,:,:)),minval(T(2,3,:,:,:)),"   ",minval(delta(2,:,:,:))
              !write(*,fmt='(3f14.8,a4,f14.8)') minval(T(3,1,:,:,:)),minval(T(3,2,:,:,:)),minval(T(3,3,:,:,:)),"   ",minval(delta(3,:,:,:))
              !print *,"Lib_FindLocalStrains::computeLocalStrainOnGrid() info - average values"
              !write(*,fmt='(3f14.8,a4,f14.8)') sum(T(1,1,:,:,:))/nn,sum(T(1,2,:,:,:))/nn,sum(T(1,3,:,:,:))/nn,"   ",sum(delta(1,:,:,:))/nn
              !write(*,fmt='(3f14.8,a4,f14.8)') sum(T(2,1,:,:,:))/nn,sum(T(2,2,:,:,:))/nn,sum(T(2,3,:,:,:))/nn,"   ",sum(delta(2,:,:,:))/nn
              !write(*,fmt='(3f14.8,a4,f14.8)') sum(T(3,1,:,:,:))/nn,sum(T(3,2,:,:,:))/nn,sum(T(3,3,:,:,:))/nn,"   ",sum(delta(3,:,:,:))/nn
              !
            
            
            
            
            
            
            
            
            
            
            
            
            nR0 = size(r0,dim=2)
            
            if (FLS_DBG) then
                print *,"Lib_FindLocalStrains::computeLocalStrainOnGrid() info - cell size hint ",a0
                print *,"Lib_FindLocalStrains::computeLocalStrainOnGrid() info - input expected points ",nR0
                call report( lc3d )
            end if
            
        !---    find number of atoms in link cell. Make initial allocation of arrays
            nAtoms = getNPoints( lc3d )
            allocate(pos(3,nAtoms))            
            
            do iz = 0,getNz(getSuper(lc3d))-1
                do iy = 0,getNy(getSuper(lc3d))-1
                    do ix = 0,getNx(getSuper(lc3d))-1
                        do ik = 1,getNPoints(lc3d,ix,iy,iz)
                            xx = getPosition( lc3d,ix,iy,iz,ik,includeCellOffset=.true. )
                            nn = getID( lc3d,ix,iy,iz,ik )
                            pos(1:3,nn) = xx(1:3)
                        end do
                    end do
                end do
            end do
            !print *,"find number of atoms in link cell ",nn,nAtoms
            if (FLS_DBG) then
                print *,"Lib_FindLocalStrains::computeLocalStrainOnGrid() info - minmaxavg positions"
                print *,minval(pos(1,:)),maxval(pos(1,:)),sum(pos(1,:))/nAtoms
                print *,minval(pos(2,:)),maxval(pos(2,:)),sum(pos(2,:))/nAtoms
                print *,minval(pos(3,:)),maxval(pos(3,:)),sum(pos(3,:))/nAtoms
            end if
            
            
        !---                
            
                        
            allocate(dx(3,getnNeighMax(lc3d)))
            allocate(dx_all(3,nAtoms))
            allocate(weight(nAtoms))
                        
            if (maxval( (/Nx,Ny,Nz/) )==1) then
                print *,"Lib_FindLocalStrains::computeLocalStrainOnGrid() info - homogeneous strain calc"
            !   only a single grid point. Add all atoms with equal weight
                nn_all = 0
                do ii = 1,nAtoms
                    call progressBar(ii,nAtoms)
                    yy = pos(1:3,ii)
                    call neighbourList( lc3d,yy,1.9d0*a0, nn,dx )
                    if (nn_all + nn>size(weight)) then
                        allocate(dx_all_tmp(3,max(2*nn_all,nn_all+2*nn)))
                        allocate(weight_tmp(max(2*nn_all,nn_all+nn)))
                        dx_all_tmp(1:3,1:nn_all) = dx_all(1:3,1:nn_all)
                        weight_tmp(1:nn_all) = weight(1:nn_all)
                        deallocate(dx_all)
                        deallocate(weight)
                        dx_all => dx_all_tmp
                        weight => weight_tmp
                    end if                                   
                    dx_all(1:3,nn_all+1:nn_all+nn) = dx(1:3,1:nn)
                    weight(nn_all+1:nn_all+nn) = 1.0d0
                    nn_all = nn_all + nn
                end do
                
            !---    expect to see r' = T r0 + d                    
                T_lcl(1:3,1:3) = T(1:3,1:3,0,0,0)
                d_lcl(1:3) = delta(1:3,0,0,0)
                do ii = 1,nR0
                    r0_lcl( 1:3,ii ) = T_lcl(1:3,1)*r0(1,ii) + T_lcl(1:3,2)*r0(2,ii) + T_lcl(1:3,3)*r0(3,ii) + d_lcl(1:3)
                end do                
                call localStrain( r0_lcl,dx_all(1:3,1:nn_all),a0, TT,dd ,ok, .false., weight )   
            !   add local disp / strain to result             
                dd(1:3) = T_lcl(1:3,1)*dd(1) + T_lcl(1:3,2)*dd(2) + T_lcl(1:3,3)*dd(3) + d_lcl(1:3)
                TT(1:3,1:3) = matmul( TT,T_lcl )
                
                if (ok) then
                    T(1:3,1:3,0,0,0) = TT(1:3,1:3)
                    delta(1:3,0,0,0) = dd(1:3)
                else                
                    T(1:3,1:3,0,0,0) = reshape( (/1,0,0 , 0,1,0 , 0,0,1/) , (/3,3/) )
                    delta(1:3,0,0,0) = 0.0d0
                    if (FLS_DBG) print *,"Lib_FindLocalStrains::computeLocalStrainOnGrid() warning - single point strain calc failed"
                end if
                
            else                        
                        
            !---    find a unit cell for grid
                
                aa = getSuperA( getSuper(lc3d) )
                aa(1:3,1) = aa(1:3,1)/Nx
                aa(1:3,2) = aa(1:3,2)/Ny
                aa(1:3,3) = aa(1:3,3)/Nz
                super = SimpleSupercell_ctor( aa,Nx,Ny,Nz )
                if (noSmooth) then
                    sigma = 0
                else
                    sigma = minval( (/getCellSideLength(super,1),getCellSideLength(super,2),getCellSideLength(super,3)/) )
                end if 
                
                !print *,"unit cell for grid",sigma
                !call report(super)
                print *,"Lib_FindLocalStrains::computeLocalStrainOnGrid() info - grid size ",Nx,"x",Ny,"x",Nz," smoothing ",sigma/a0   
                
                
            !---    loop through each grid point, and use weights
                allocate(id(nAtoms))
                allocate(dr2(nAtoms))            
                i2s2 = 1/(2*sigma*sigma/2)
                !FLS_DBG = .true.
                nFails = 0
                do iz = 0,Nz-1
                    do iy = 0,Ny-1
                        do ix = 0,Nx-1
                        
                            call progressBar(ix+Nx*(iy+Ny*iz),Nx*Ny*Nz-1)
                        
                            xx = (/ix,iy,iz/) !+ 0.5d0                 !   half grid point
                            xx = cellToRealSpace( super,xx )          !   position of point in real space on my grid
                            
                            
                            !FLS_DBG = ( (ix==0).and.(iy==0).and.(iz==Nx/2) ).or.( (ix==0).and.(iy==0).and.(iz==0) )
                            !LinkCell3D_dbg = FLS_DBG
                            
                        !---    find the long range atom neighbours of this point in real space
                            call neighbourList_long( lc3d,xx,max(a0,sigma), nn2,id,dr2 )
                            if (FLS_DBG) print *,"Lib_FindLocalStrains::computeLocalStrainOnGrid() info - neighbourList_long range ",max(a0,sigma)," returns ",nn2," neighbours"
                            
                          
                            if (noSmooth) then
                                nn_all = 0
                                ii = minloc(dr2(1:nn2),dim=1)
                                yy = pos(1:3,id(ii))
                                call neighbourList( lc3d,yy,1.9d0*a0, nn,dx )
                                if (nn_all + nn>size(weight)) then
                                    allocate(dx_all_tmp(3,max(2*nn_all,nn_all+2*nn)))
                                    allocate(weight_tmp(max(2*nn_all,nn_all+nn)))
                                    dx_all_tmp(1:3,1:nn_all) = dx_all(1:3,1:nn_all)
                                    weight_tmp(1:nn_all) = weight(1:nn_all)
                                    deallocate(dx_all)
                                    deallocate(weight)
                                    dx_all => dx_all_tmp
                                    weight => weight_tmp
                                end if
                                dx_all(1:3,nn_all+1:nn_all+nn) = dx(1:3,1:nn)
                                ww = exp( - dr2(ii)*i2s2 )
                                weight(nn_all+1:nn_all+nn) = ww
                                nn_all = nn_all + nn                            
                            else                                    
                            !---    find the short range atom neighbours of each one of these 
                                nn_all = 0
                                do ii = 1,nn2
                                    yy = pos(1:3,id(ii))
                                    call neighbourList( lc3d,yy,1.9d0*a0, nn,dx )
                                    if (nn_all + nn>size(weight)) then
                                        allocate(dx_all_tmp(3,max(2*nn_all,nn_all+2*nn)))
                                        allocate(weight_tmp(max(2*nn_all,nn_all+nn)))
                                        dx_all_tmp(1:3,1:nn_all) = dx_all(1:3,1:nn_all)
                                        weight_tmp(1:nn_all) = weight(1:nn_all)
                                        deallocate(dx_all)
                                        deallocate(weight)
                                        dx_all => dx_all_tmp
                                        weight => weight_tmp
                                    end if
                                    dx_all(1:3,nn_all+1:nn_all+nn) = dx(1:3,1:nn)
                                    ww = exp( - dr2(ii)*i2s2 )
                                    weight(nn_all+1:nn_all+nn) = ww
                                    nn_all = nn_all + nn
                                end do
                            end if                            
                            
                            
                            
                            !call localStrain( r0,dx_all(1:3,1:nn_all),a0, TT,dd ,ok, .false., weight )
                        !---    expect to see r' = T r0 + d                    
                            T_lcl(1:3,1:3) = T(1:3,1:3,ix,iy,iz)
                            d_lcl(1:3) = delta(1:3,ix,iy,iz)
                            do ii = 1,nR0
                                r0_lcl( 1:3,ii ) = T_lcl(1:3,1)*r0(1,ii) + T_lcl(1:3,2)*r0(2,ii) + T_lcl(1:3,3)*r0(3,ii) + d_lcl(1:3)
                            end do                
                            
                            
                            
                            !print *,"d_lcl ",d_lcl
                            !print *,"nn    ",nn2,nn_all
                            
                            
                            
                            
                            ! weight = 1
                            
                            call localStrain( r0_lcl,dx_all(1:3,1:nn_all),a0, TT,dd ,ok, .false., weight )   
                        !   add local disp / strain to result             
                            dd(1:3) = T_lcl(1:3,1)*dd(1) + T_lcl(1:3,2)*dd(2) + T_lcl(1:3,3)*dd(3) + d_lcl(1:3)
                            TT(1:3,1:3) = matmul( TT,T_lcl )
                            
                                                        
!                          if ( (ix == Nx/2).and.(iy == Ny/2).and.(iz == Nz/2).and.(Nx>=32) ) then
!                              print *,"dbg local strain"
!                              print *,ix,iy,iz
!                              do ii = 1,nn_all
!                                  print *,ii,weight(ii),dx_all(1:3,ii)
!                              end do                            
!                              print *,"delta ",dd(1:3)
!                          end if
!                       
                            !write (*,fmt='(a,3i4,3f12.5,a,2i8,a,3f12.5,a,100f12.5)') "grid ",ix,iy,iz,xx," nn    ",nn2,nn_all," dd    ",dd," tt    ",TT                            
                            
                            if (ok) then
                                T(1:3,1:3,ix,iy,iz) = TT(1:3,1:3)
                                delta(1:3,ix,iy,iz) = dd(1:3)
                            else                
                                T(1:3,1:3,ix,iy,iz) = reshape( (/1,0,0 , 0,1,0 , 0,0,1/) , (/3,3/) )
                                delta(1:3,ix,iy,iz) = 0.0d0
                                nFails = nFails + 1
                            end if
                            
                            !if ((ix==0).and.(iy==0).and.(iz==1)) stop
                            
                        end do
                    end do
                end do
                
                nn = Nx*Ny*Nz
                print *,"Lib_FindLocalStrains::computeLocalStrainOnGrid() info - nFails = ",nFails,"/",nn
                
                print *,"Lib_FindLocalStrains::computeLocalStrainOnGrid() info - maximum values"
                write(*,fmt='(3f14.8,a4,f14.8)') maxval(T(1,1,:,:,:)),maxval(T(1,2,:,:,:)),maxval(T(1,3,:,:,:)),"   ",maxval(delta(1,:,:,:))
                write(*,fmt='(3f14.8,a4,f14.8)') maxval(T(2,1,:,:,:)),maxval(T(2,2,:,:,:)),maxval(T(2,3,:,:,:)),"   ",maxval(delta(2,:,:,:))
                write(*,fmt='(3f14.8,a4,f14.8)') maxval(T(3,1,:,:,:)),maxval(T(3,2,:,:,:)),maxval(T(3,3,:,:,:)),"   ",maxval(delta(3,:,:,:))
                print *,"Lib_FindLocalStrains::computeLocalStrainOnGrid() info - minimum values"
                write(*,fmt='(3f14.8,a4,f14.8)') minval(T(1,1,:,:,:)),minval(T(1,2,:,:,:)),minval(T(1,3,:,:,:)),"   ",minval(delta(1,:,:,:))
                write(*,fmt='(3f14.8,a4,f14.8)') minval(T(2,1,:,:,:)),minval(T(2,2,:,:,:)),minval(T(2,3,:,:,:)),"   ",minval(delta(2,:,:,:))
                write(*,fmt='(3f14.8,a4,f14.8)') minval(T(3,1,:,:,:)),minval(T(3,2,:,:,:)),minval(T(3,3,:,:,:)),"   ",minval(delta(3,:,:,:))
                print *,"Lib_FindLocalStrains::computeLocalStrainOnGrid() info - average values"
                write(*,fmt='(3f14.8,a4,f14.8)') sum(T(1,1,:,:,:))/nn,sum(T(1,2,:,:,:))/nn,sum(T(1,3,:,:,:))/nn,"   ",sum(delta(1,:,:,:))/nn
                write(*,fmt='(3f14.8,a4,f14.8)') sum(T(2,1,:,:,:))/nn,sum(T(2,2,:,:,:))/nn,sum(T(2,3,:,:,:))/nn,"   ",sum(delta(2,:,:,:))/nn
                write(*,fmt='(3f14.8,a4,f14.8)') sum(T(3,1,:,:,:))/nn,sum(T(3,2,:,:,:))/nn,sum(T(3,3,:,:,:))/nn,"   ",sum(delta(3,:,:,:))/nn
                
             end if
                                
            deallocate(weight)
            deallocate(dx_all)
                
            return
        end subroutine computeLocalStrainOnGrid
    
        subroutine computeLocalDisplacementOnGrid( lc3d, r0,a0, delta )
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      compute the strain at evenly spaced points, within a width sigma
    !*      given the expected separations r0
    !*      on input, T, delta are expected roto-strain transformation and displacement
    !*      on exit, T, delta are better fitted 
    !*          sense:  r = T r0 + delta
    !*      values are computed at helf-grid points 1/2,1/2,1/2  not 0,0,0
    
            type(LinkCell3d),intent(in)                     ::      lc3d
            real(kind=real64),dimension(:,:),intent(in)     ::      r0
            real(kind=real64),intent(in)        ::      a0   
            real(kind=real64),dimension(:,0:,0:,0:),intent(inout)    ::      delta
            
            real(kind=real64),dimension(3,3)    ::      aa , TT
            real(kind=real64),dimension(:,:),pointer        ::      dx,dx_all,dx_all_tmp,pos
            integer,dimension(:),allocatable                ::      id
            real(kind=real64),dimension(:),pointer          ::      weight,weight_tmp
            real(kind=real64),dimension(:),allocatable      ::      dr2
            real(kind=real64),dimension(3,size(r0,dim=2))   ::      r0_lcl
            type(SimpleSupercell)           ::      super
            integer             ::      Nx,Ny,Nz, ix,iy,iz,ik, ii  !,mx,my,mz
            real(kind=real64)   ::      sigma,i2s2,ww
            real(kind=real64),dimension(3)  ::      xx,yy,dd,d_lcl
            integer             ::      nAtoms,nn,nn_all,nn2,nFails,nR0
            logical             ::      ok
            
        !---    find the size of the grid
            Nx = size(delta,dim=2)
            Ny = size(delta,dim=3)
            Nz = size(delta,dim=4)

            nR0 = size(r0,dim=2)
            print *,"Lib_FindLocalStrains::computeLocalDisplacementOnGrid() info - grid size ",Nx,"x",Ny,"x",Nz    
            if (FLS_DBG) then
                print *,"Lib_FindLocalStrains::computeLocalDisplacementOnGrid() info - cell size hint ",a0
                print *,"Lib_FindLocalStrains::computeLocalDisplacementOnGrid() info - input expected points ",nR0
                call report( lc3d )
            end if
            
        !---    find number of atoms in link cell. Make initial allocation of arrays
            nAtoms = getNPoints( lc3d )
            allocate(pos(3,nAtoms))            
            
            do iz = 0,getNz(getSuper(lc3d))-1
                do iy = 0,getNy(getSuper(lc3d))-1
                    do ix = 0,getNx(getSuper(lc3d))-1
                        do ik = 1,getNPoints(lc3d,ix,iy,iz)
                            xx = getPosition( lc3d,ix,iy,iz,ik,includeCellOffset=.true. )
                            nn = getID( lc3d,ix,iy,iz,ik )
                            pos(1:3,nn) = xx(1:3)
                        end do
                    end do
                end do
            end do
            !print *,"find number of atoms in link cell ",nn,nAtoms
            if (FLS_DBG) then
                print *,"Lib_FindLocalStrains::computeLocalDisplacementOnGrid() info - minmaxavg positions"
                print *,minval(pos(1,:)),maxval(pos(1,:)),sum(pos(1,:))/nAtoms
                print *,minval(pos(2,:)),maxval(pos(2,:)),sum(pos(2,:))/nAtoms
                print *,minval(pos(3,:)),maxval(pos(3,:)),sum(pos(3,:))/nAtoms
            end if
            
            
        !---                
            
                        
            allocate(dx(3,getnNeighMax(lc3d)))
            allocate(dx_all(3,nAtoms))
            allocate(weight(nAtoms))
                        
            if (maxval( (/Nx,Ny,Nz/) )==1) then
            
            !   only a single grid point. Add all atoms with equal weight
                nn_all = 0
                do ii = 1,nAtoms
                    yy = pos(1:3,ii)
                    call neighbourList( lc3d,yy,1.9d0*a0, nn,dx )
                    if (nn_all + nn>size(weight)) then
                        allocate(dx_all_tmp(3,max(2*nn_all,nn_all+2*nn)))
                        allocate(weight_tmp(max(2*nn_all,nn_all+nn)))
                        dx_all_tmp(1:3,1:nn_all) = dx_all(1:3,1:nn_all)
                        weight_tmp(1:nn_all) = weight(1:nn_all)
                        deallocate(dx_all)
                        deallocate(weight)
                        dx_all => dx_all_tmp
                        weight => weight_tmp
                    end if                                   
                    dx_all(1:3,nn_all+1:nn_all+nn) = dx(1:3,1:nn)
                    weight(nn_all+1:nn_all+nn) = 1.0d0
                    nn_all = nn_all + nn
                end do
                
            !---    expect to see r' = T r0 + d                    
                d_lcl(1:3) = delta(1:3,0,0,0)
                do ii = 1,nR0
                    r0_lcl( 1:3,ii ) = d_lcl(1:3)
                end do                
                call localStrain( r0_lcl,dx_all(1:3,1:nn_all),a0, TT,dd ,ok, .true., weight )   
            !   add local disp / strain to result             
                dd(1:3) = dd(1:3) + d_lcl(1:3)                 
                
                if (ok) then
                    delta(1:3,0,0,0) = dd(1:3)
                else                
                    delta(1:3,0,0,0) = 0.0d0
                    if (FLS_DBG) print *,"Lib_FindLocalStrains::computeLocalDisplacementOnGrid() warning - single point strain calc failed"
                end if
                
            else                        
                        
            !---    find a unit cell for grid
            
                
                aa = getSuperA( getSuper(lc3d) )
                aa(1:3,1) = aa(1:3,1)/Nx
                aa(1:3,2) = aa(1:3,2)/Ny
                aa(1:3,3) = aa(1:3,3)/Nz
                super = SimpleSupercell_ctor( aa,Nx,Ny,Nz )
                sigma = minval( (/getCellSideLength(super,1),getCellSideLength(super,2),getCellSideLength(super,3)/) )
                
                !print *,"unit cell for grid",sigma
                !call report(super)
                
                
            !---    loop through each grid point, and use weights
                allocate(id(nAtoms))
                allocate(dr2(nAtoms))            
                i2s2 = 1/(2*sigma*sigma/2)
                !FLS_DBG = .true.
                nFails = 0
                do iz = 0,Nz-1
                    do iy = 0,Ny-1
                        do ix = 0,Nx-1
                            xx = (/ix,iy,iz/)                  
                            xx = cellToRealSpace( super,xx )          !   position of point in real space on my grid
                            
                            
                            call progressBar(ix+1+Nx*(iy+Ny*iz),Nx*Ny*Nz)
                            
                            !FLS_DBG = ( (ix==0).and.(iy==0).and.(iz==Nx/2) ).or.( (ix==0).and.(iy==0).and.(iz==0) )
                            !LinkCell3D_dbg = FLS_DBG
                            
                        !---    find the long range atom neighbours of this point in real space
                            call neighbourList_long( lc3d,xx,max(a0,sigma), nn2,id,dr2 )
                            if (FLS_DBG) print *,"Lib_FindLocalStrains::computeLocalDisplacementOnGrid() info - neighbourList_long range ",1.5*sigma," returns ",nn2," neighbours"
                            
                        !---    find the short range atom neighbours of each one of these 
                            nn_all = 0
                            do ii = 1,nn2
                                yy = pos(1:3,id(ii))
                                call neighbourList( lc3d,yy,1.9d0*a0, nn,dx )
                                if (nn_all + nn>size(weight)) then
                                    allocate(dx_all_tmp(3,max(2*nn_all,nn_all+2*nn)))
                                    allocate(weight_tmp(max(2*nn_all,nn_all+nn)))
                                    dx_all_tmp(1:3,1:nn_all) = dx_all(1:3,1:nn_all)
                                    weight_tmp(1:nn_all) = weight(1:nn_all)
                                    deallocate(dx_all)
                                    deallocate(weight)
                                    dx_all => dx_all_tmp
                                    weight => weight_tmp
                                end if
                                dx_all(1:3,nn_all+1:nn_all+nn) = dx(1:3,1:nn)
                                ww = exp( - dr2(ii)*i2s2 )
                                weight(nn_all+1:nn_all+nn) = ww
                                nn_all = nn_all + nn
                            end do
                            
                            !call localStrain( r0,dx_all(1:3,1:nn_all),a0, TT,dd ,ok, .false., weight )
                        !---    expect to see r' = T r0 + d                    
                            d_lcl(1:3) = delta(1:3,ix,iy,iz)
                            do ii = 1,nR0
                                r0_lcl( 1:3,ii ) = r0(1:3,ii) + d_lcl(1:3)
                            end do                
                                                        
                            
                            call localStrain( r0_lcl,dx_all(1:3,1:nn_all),a0, TT,dd ,ok, .true., weight )   
                        !   add local disp / strain to result             
                            dd(1:3) = dd(1:3) + d_lcl(1:3)
                            
                            !write (*,fmt='(a,3i4,3f12.5,a,2i8,a,3f12.5,a,100f12.5)') "grid ",ix,iy,iz,xx," nn    ",nn2,nn_all," dd    ",dd," tt    ",TT                            
                            
                            if (ok) then
                                delta(1:3,ix,iy,iz) = dd(1:3)
                            else                
                                delta(1:3,ix,iy,iz) = 0.0d0
                                nFails = nFails + 1
                            end if
                            
                            !if ((ix==0).and.(iy==0).and.(iz==1)) stop
                            
                        end do
                    end do
                end do
                
                nn = Nx*Ny*Nz
                print *,"Lib_FindLocalStrains::computeLocalDisplacementOnGrid() info - nFails = ",nFails,"/",nn
                
                print *,"Lib_FindLocalStrains::computeLocalDisplacementOnGrid() info - maximum values"
                write(*,fmt='(3f14.8,a4,f14.8)') maxval(delta(1,:,:,:)) ,maxval(delta(2,:,:,:)),maxval(delta(3,:,:,:))
                print *,"Lib_FindLocalStrains::computeLocalDisplacementOnGrid() info - minimum values"
                write(*,fmt='(3f14.8,a4,f14.8)') minval(delta(1,:,:,:)),minval(delta(2,:,:,:)),minval(delta(3,:,:,:))
                print *,"Lib_FindLocalStrains::computeLocalDisplacementOnGrid() info - average values"
                write(*,fmt='(3f14.8,a4,f14.8)') sum(delta(1,:,:,:))/nn,sum(delta(2,:,:,:))/nn,sum(delta(3,:,:,:))/nn
                
             end if
                                
                
            return
        end subroutine computeLocalDisplacementOnGrid
    
            
        subroutine localStrain( r0,r,a0, T,delta ,ok, displacementOnly, weight )
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      expect to see atoms at separations r0
    !*      actually see atoms at separations r
    !*      return the affine transformation T,delta
    !*      to best fit 
    !*          assume r = T r0 + delta
    !*      with T = (1 + eps) rot is the complete basis transform
    !*      each expected atom r0 can be matched to multiple observed atoms r
    !*      if displacementOnly then don't bother computing the rotostrain transform T - return identity
    !*      use weight to determine the relative importance of each of the input atom separations
    
            real(kind=real64),dimension(:,:),intent(in)     ::      r0
            real(kind=real64),dimension(:,:),intent(in)     ::      r
            real(kind=real64),intent(in)                    ::      a0
            real(kind=real64),dimension(3,3),intent(out)    ::      T
            real(kind=real64),dimension(3),intent(out)      ::      delta
            logical,intent(out)                             ::      ok
            logical,intent(in)                              ::      displacementOnly
            real(kind=real64),dimension(:),intent(in)       ::      weight
            
            integer                                 ::      n0,nn,nPaired           !   count of expected atoms, count of input atoms, number of pairs made
            integer                                 ::      ii,jj,kk,mm
            integer,dimension(:,:),allocatable      ::      indx            !   (2,1:n)       r(1,indx(i)) matches to r0(2,indx(i))
            
            real(kind=real64)               ::      dd,dc,dmin,ww,wsum
            real(kind=real64),dimension(3)  ::      dx
            
            real(kind=real64),dimension(12,12)      ::      AA
            real(kind=real64),dimension(12)         ::      BB
            integer,dimension(12)                   ::      ipiv
            real(kind=real64),dimension(12*66)      ::      work

        !---    check for sensible input
            n0 = size(r0,dim=2)
            nn = size(r,dim=2)            
            T = reshape( (/1,0,0, 0,1,0, 0,0,1/) , (/3,3/) )
            delta = 0.0
            ok = .true.
            if (nn<4) then
            !   no chance of finding a correct anything with 3 input points.
                ok = .false.
                return
            end if
                
        !---    friendly help wanted?                        
             if (FLS_DBG) then
                 delta = 0 ; wsum = 0.0d0
                 do ii = 1,nn
                     delta = delta + r(1:3,ii)
                     wsum = wsum + weight(ii)
                 end do
                 print *,"observed mean position ",nn,delta/nn
                 print *,"weighted mean position ",wsum,delta/wsum
             
                 delta = 0
                 do ii = 1,n0
                     delta = delta + r0(1:3,ii)
                 end do
                 print *,"expected mean position ",n0,delta/n0
                 print *,"average obs/exp        ",real(nn)/n0
                 print *,"hint a0                ",a0                
             end if
            
            
        !---    first task is to order r the same as r0.            
        !       find the distance between each pair of observed/expected            
        !       use this to find best match
        !       we can also compute the mean separation while we are doing this task.
            dc = a0*a0/16               !   make the characteristic max separation a0/4.            
            allocate(indx(2,nn)) 
            delta = 0.0d0
            nPaired = 0                      !   number of matches made
            wsum = 0.0d0
            do ii = 1,nn
            
            !---    find closest match in expected input within dc
                dmin = dc ; mm = 0      !   mm will be the match index in expected input    
                do jj = 1,n0
                    dx(1:3) = r(1:3,ii) - r0(1:3,jj) 
                    dd = dx(1)*dx(1) + dx(2)*dx(2) + dx(3)*dx(3)
                    if ( dd <= dmin )then
                        dmin = dd ; mm = jj
                    end if
                end do
                
            !---    add match to indexing
                if ( mm /= 0) then
                    nPaired = nPaired + 1
                    indx(1:2,nPaired) = (/ii,mm/) 
                    delta(1:3) = delta(1:3) + weight(ii)*( r(1:3,ii) - r0(1:3,mm) )
                    wsum = wsum + weight(ii)
                end if
                
            end do
            
            if (FLS_DBG) then
                print *,"pairs per observed     ",real(nPaired)/nn
                print *,"pairs per expected     ",real(nPaired)/n0
                print *,"avg weight per pair    ",wsum / nPaired      
                print *,"mean offset            ",delta / max(1.0d-16,wsum)
            end if
            
            
        !---    catastrophic fail?
            if (nPaired*wsum == 0) then
                ok = .false.
                return
            end if
            
            delta = delta / wsum                
                             
        !---    are we done?
            if (displacementOnly) return
            
        !---    now the easy bit - linear least squares
        !       S = sum   w_i ( (T (r0_j+delta') + delta) - r_i )^2
        !       where sum is over pairs and observed i matched to expected j
        !       delta' is the best fit to centre of mass computed above
        !
        !       dS/d T_ab    = 2 sum w_i ( (T_ac (r0_jc+delta'_c) + delta_a) - r_ia ) (r0_jb+delta'_b)        
        !       dS/d delta_a = 2 sum w_i ( (T_ab (r0_jb+delta'_b) + delta_a) - r_ia ) 
        !
        !
                                    
            AA = 0.0d0
            BB = 0.0d0      !   ordering T11,T21,T31,T12,T22,T32,T13,T23,T33,d1,d2,d3
            
                      
            AA(10,10) = wsum
            AA(11,11) = wsum
            AA(12,12) = wsum                
            
            do kk = 1,nPaired           !   loop over matched pairs
                ii = indx(1,kk)         !   observed atom
                jj = indx(2,kk)         !   expected atom
                ww = weight(ii)
                                
            !---    dS/dT  ( lower corner ) 
                dd = ww*(r0(1,jj)+delta(1))*(r0(1,jj)+delta(1)) ; AA(1,1) = AA(1,1) + dd ; AA(2,2) = AA(2,2) + dd ; AA(3,3) = AA(3,3) + dd                 
                dd = ww*(r0(2,jj)+delta(2))*(r0(1,jj)+delta(1)) ; AA(4,1) = AA(4,1) + dd ; AA(5,2) = AA(5,2) + dd ; AA(6,3) = AA(6,3) + dd      
                dd = ww*(r0(3,jj)+delta(3))*(r0(1,jj)+delta(1)) ; AA(7,1) = AA(7,1) + dd ; AA(8,2) = AA(8,2) + dd ; AA(9,3) = AA(9,3) + dd  
                dd = ww*(r0(2,jj)+delta(2))*(r0(2,jj)+delta(2)) ; AA(4,4) = AA(4,4) + dd ; AA(5,5) = AA(5,5) + dd ; AA(6,6) = AA(6,6) + dd 
                dd = ww*(r0(3,jj)+delta(3))*(r0(2,jj)+delta(2)) ; AA(7,4) = AA(7,4) + dd ; AA(8,5) = AA(8,5) + dd ; AA(9,6) = AA(9,6) + dd 
                dd = ww*(r0(3,jj)+delta(3))*(r0(3,jj)+delta(3)) ; AA(7,7) = AA(7,7) + dd ; AA(8,8) = AA(8,8) + dd ; AA(9,9) = AA(9,9) + dd 
                         
            !---    dS/ddelta                                       
                dd = ww*(r0(1,jj)+delta(1)) ; BB(1:3) = BB(1:3) + r(1:3,ii)*dd ; AA(10,1) = AA(10,1) + dd ; AA(11,2) = AA(11,2) + dd ; AA(12,3) = AA(12,3) + dd                
                dd = ww*(r0(2,jj)+delta(2)) ; BB(4:6) = BB(4:6) + r(1:3,ii)*dd ; AA(10,4) = AA(10,4) + dd ; AA(11,5) = AA(11,5) + dd ; AA(12,6) = AA(12,6) + dd  
                dd = ww*(r0(3,jj)+delta(3)) ; BB(7:9) = BB(7:9) + r(1:3,ii)*dd ; AA(10,7) = AA(10,7) + dd ; AA(11,6) = AA(11,6) + dd ; AA(12,9) = AA(12,9) + dd  
                                
                BB(10:12) = BB(10:12) + ww*r(1:3,ii)
                                            
            end do
           
             if (FLS_DBG) then
                 do ii = 1,12
                     do jj = ii+1,12
                         AA(ii,jj) = AA(jj,ii)
                     end do
                     write (*,fmt='(g14.6,a,12g14.6)') BB(ii),"    ",AA(ii,:)
                 end do
            end if
            
            call DSYSV( "L",12,1,AA(1:12,1:12),12,ipiv,BB(1:12),12,work,size(work),ii )                    
                    
            T(1:3,1) = BB(1:3)
            T(1:3,2) = BB(4:6)
            T(1:3,3) = BB(7:9)
            delta(1:3) = T(1:3,1)*delta(1) + T(1:3,2)*delta(2) + T(1:3,3)*delta(3)  + BB(10:12)
            
            
            ok = (ii==0)
           
            
             if (FLS_DBG) then
                 print *,"dsysv ",ok
                 do ii = 1,3
                     write (*,fmt='(3g14.6,a,g14.6,a,g14.6)') T(ii,:),"    ",delta(ii)
                 end do                            
             end if
            
            
            return
        end subroutine localStrain
            
            
            



    end module Lib_FindLocalStrains
         
    