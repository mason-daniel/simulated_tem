
    module Lib_TopologicallyCorrectClustering
!---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
!*      given a voxelised field on a roughly cubic lattice,
!*      return clusters for voxels over 0

        use iso_fortran_env
        implicit none
        private
        
        
        public      ::      topologicallyCorrectCluster
        public      ::      findClusterExtents
        public      ::      cropField
         
        
        logical,public      ::      DBG_TCC = .false.
        
    contains
!---^^^^^^^^

        subroutine topologicallyCorrectCluster( f , indx,nClust , invertLarge)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      assumes input f is periodic
    !*      if this is not desired, just add a plane of -huge(1.0) and the clusters will be split
    !*      if invertLarge, then make large 
            real(kind=real64),dimension(0:,0:,0:),intent(in)        ::      f
            integer,dimension(0:,0:,0:),intent(out)                 ::      indx
            integer,intent(out)                                     ::      nClust
            logical,intent(in),optional                             ::      invertLarge
            
            integer             ::      Nx,Ny,Nz
            
            integer             ::      ix,iy,iz , jx,jy,jz ,ii, i1,i2 ,kx,ky,kz
            
                        
            integer             ::      highIndx
            logical             ::      ok
            integer,dimension(:),allocatable        ::      reindx
            integer,dimension(:,:,:),allocatable    ::      indx_tmp
            Nx = size(f,dim=1)
            Ny = size(f,dim=2)
            Nz = size(f,dim=3)        
            
            
        !---    first step: identify those pixels above threshold. Don't worry about clustering at the moment...
            !print *,"lit voxels ",count(f>=0)
            nClust = 0      
            highIndx = Nx*Ny*Nz+1
            indx = highIndx
            do iz = 0,Nz-1
                do iy = 0,Ny-1
                    do ix = 0,Nx-1
                        if (f(ix,iy,iz)>=0) then
                            nClust = nClust + 1
                            indx(ix,iy,iz) = nClust
                            !print *,"lit ",nClust,ix,iy,iz
                        end if
                    end do
                end do
            end do
            
            if (DBG_TCC) print *,"Lib_topologicallyCorrectClustering::topologicallyCorrectCluster() info - number of lit voxels ",nClust,"/",highIndx-1
            
                    
        !---    second step: combine clusters             
        
        !   first do a quick and dirty scan of nearest neighbours only. This will produce too many clusters, but they won't be wrong. This is good for merging large volumes.
        !   then do the full topologically correct scan.
            print *,"Lib_topologicallyCorrectClustering::topologicallyCorrectCluster() info - quick forwards cluster pass " 
            do ii = 1,highIndx
                ok = .true.
                
                do iz = 0,Nz-1
                    jz = mod( iz + 1,Nz )
                    kz = mod( iz + Nz - 1,Nz )
                    do iy = 0,Ny-1
                        jy = mod( iy + 1,Ny )
                        ky = mod( iy + Ny - 1,Ny )
                        do ix = 0,Nx-1
                            jx = mod( ix + 1,Nx )
                            kx = mod( ix + Nx - 1,Nx )
                            
                            i1 = indx(ix,iy,iz)
                            if (i1==highIndx) cycle        !   voxel ix,iy,iz is not lit. Do not connect anything to it
                                        
                            
                        !---    nearest neighbours                                                
                            i2 = indx(jx,iy,iz)
                            if ((i2-highIndx)*(i2-i1)/=0) then
                                i1 = min(i1,i2)
                                indx(ix,iy,iz) = i1
                                indx(jx,iy,iz) = i1
                                ok = .false.
                            end if                        
                                                              
                            i2 = indx(ix,jy,iz)
                            if ((i2-highIndx)*(i2-i1)/=0) then
                                i1 = min(i1,i2)
                                indx(ix,iy,iz) = i1
                                indx(ix,jy,iz) = i1
                                ok = .false.
                            end if                        
                                                
                            i2 = indx(ix,iy,jz)
                            if ((i2-highIndx)*(i2-i1)/=0) then
                                i1 = min(i1,i2)
                                indx(ix,iy,iz) = i1
                                indx(ix,iy,jz) = i1
                                ok = .false.
                            end if          
                                
                        end do
                    end do
                end do
                if (ok) exit
            end do

        
        
            print *,"Lib_topologicallyCorrectClustering::topologicallyCorrectCluster() info - topologically correct cluster pass " 
            allocate(indx_tmp(0:Nx-1,0:Ny-1,0:Nz-1))            
            do ii = 1,highIndx
                ok = .true.
                !print *,"Lib_topologicallyCorrectClustering::topologicallyCorrectCluster() info - topologically correct cluster pass ",ii
                if (DBG_TCC) call output()

                indx_tmp = indx
                do iz = 0,Nz-1
                    jz = mod( iz + 1,Nz )
                    kz = mod( iz + Nz - 1,Nz )
                    do iy = 0,Ny-1
                        jy = mod( iy + 1,Ny )
                        ky = mod( iy + Ny - 1,Ny )
                        do ix = 0,Nx-1
                            jx = mod( ix + 1,Nx )
                            kx = mod( ix + Nx - 1,Nx )
                             
                            i1 = indx(ix,iy,iz)
                            if (i1==highIndx) cycle        !   voxel ix,iy,iz is not lit. Do not connect anything to it
                                        
                            
                        !---    nearest neighbours                                                
                            i2 = indx(jx,iy,iz)
                            if ((i2-highIndx)*(i2-i1)/=0) then
                                i1 = min(i1,i2)
                                indx_tmp(ix,iy,iz) = i1
                                indx_tmp(jx,iy,iz) = i1
                                ok = .false.
                            end if                        
                                                              
                            i2 = indx(ix,jy,iz)
                            if ((i2-highIndx)*(i2-i1)/=0) then
                                i1 = min(i1,i2)
                                indx_tmp(ix,iy,iz) = i1
                                indx_tmp(ix,jy,iz) = i1
                                ok = .false.
                            end if                        
                                                
                            i2 = indx(ix,iy,jz)
                            if ((i2-highIndx)*(i2-i1)/=0) then
                                i1 = min(i1,i2)
                                indx_tmp(ix,iy,iz) = i1
                                indx_tmp(ix,iy,jz) = i1
                                ok = .false.
                            end if          
                                           
                       
                        !---    second nearest neighbours                                                
                            i2 = indx(jx,jy,iz)
                            if ((i2-highIndx)*(i2-i1)/=0) then
                                if (f(ix,iy,iz)+f(jx,iy,iz)+f(ix,jy,iz)+f(jx,jy,iz) >= 0) then
                                    i1 = min(i1,i2)
                                    indx_tmp(ix,iy,iz) = i1
                                    indx_tmp(jx,jy,iz) = i1
                                    ok = .false.
                                end if
                            end if                        
                                                              
                            i2 = indx(jx,iy,jz)
                            if ((i2-highIndx)*(i2-i1)/=0) then
                                if (f(ix,iy,iz)+f(jx,iy,iz)+f(ix,iy,jz)+f(jx,iy,jz) >= 0) then
                                    i1 = min(i1,i2)
                                    indx_tmp(ix,iy,iz) = i1
                                    indx_tmp(jx,iy,jz) = i1
                                    ok = .false.
                                end if
                            end if                            
                                        
                            i2 = indx(ix,jy,jz)
                            if ((i2-highIndx)*(i2-i1)/=0) then
                                if (f(ix,iy,iz)+f(ix,jy,iz)+f(ix,iy,jz)+f(ix,jy,jz) >= 0) then
                                    i1 = min(i1,i2)
                                    indx_tmp(ix,iy,iz) = i1
                                    indx_tmp(ix,jy,jz) = i1
                                    ok = .false.
                                end if
                            end if
                            
                            i2 = indx(jx,ky,iz)
                            if ((i2-highIndx)*(i2-i1)/=0) then
                                if (f(ix,iy,iz)+f(jx,iy,iz)+f(ix,ky,iz)+f(jx,ky,iz) >= 0) then
                                    i1 = min(i1,i2)
                                    indx_tmp(ix,iy,iz) = i1
                                    indx_tmp(jx,ky,iz) = i1
                                    ok = .false.
                                end if
                            end if                        
                                                              
                            i2 = indx(jx,iy,kz)
                            if ((i2-highIndx)*(i2-i1)/=0) then
                                if (f(ix,iy,iz)+f(jx,iy,iz)+f(ix,iy,kz)+f(jx,iy,kz) >= 0) then
                                    i1 = min(i1,i2)
                                    indx_tmp(ix,iy,iz) = i1
                                    indx_tmp(jx,iy,kz) = i1
                                    ok = .false.
                                end if
                            end if                            
                                        
                            i2 = indx(ix,jy,kz)
                            if ((i2-highIndx)*(i2-i1)/=0) then
                                if (f(ix,iy,iz)+f(ix,jy,iz)+f(ix,iy,kz)+f(ix,jy,kz) >= 0) then
                                    i1 = min(i1,i2)
                                    indx_tmp(ix,iy,iz) = i1
                                    indx_tmp(ix,jy,kz) = i1
                                    ok = .false.
                                end if
                            end if
                            
                           
                        
                        !---    third nearest neighbours                                                
                            i2 = indx(jx,jy,jz)
                            if ((i2-highIndx)*(i2-i1)/=0) then
                                if (f(ix,iy,iz)+f(jx,iy,iz)+f(ix,jy,iz)+f(jx,jy,iz)+f(ix,iy,jz)+f(jx,iy,jz)+f(ix,jy,jz)+f(jx,jy,jz) >= 0) then
                                    i1 = min(i1,i2)
                                    indx_tmp(ix,iy,iz) = i1
                                    indx_tmp(jx,jy,jz) = i1
                                    ok = .false.
                                end if
                            end if                        
                                                   
                            i2 = indx(kx,jy,jz)
                            if ((i2-highIndx)*(i2-i1)/=0) then
                                if (f(ix,iy,iz)+f(kx,iy,iz)+f(ix,jy,iz)+f(kx,jy,iz)+f(ix,iy,jz)+f(kx,iy,jz)+f(ix,jy,jz)+f(kx,jy,jz) >= 0) then
                                    i1 = min(i1,i2)
                                    indx_tmp(ix,iy,iz) = i1
                                    indx_tmp(kx,jy,jz) = i1
                                    ok = .false.
                                end if
                            end if      
                            
                            i2 = indx(jx,ky,jz)
                            if ((i2-highIndx)*(i2-i1)/=0) then
                                if (f(ix,iy,iz)+f(jx,iy,iz)+f(ix,ky,iz)+f(jx,ky,iz)+f(ix,iy,jz)+f(jx,iy,jz)+f(ix,ky,jz)+f(jx,ky,jz) >= 0) then
                                    i1 = min(i1,i2)
                                    indx_tmp(ix,iy,iz) = i1
                                    indx_tmp(jx,ky,jz) = i1
                                    ok = .false.
                                end if
                            end if                        
                                      
                            i2 = indx(jx,jy,kz)
                            if ((i2-highIndx)*(i2-i1)/=0) then
                                if (f(ix,iy,iz)+f(jx,iy,iz)+f(ix,jy,iz)+f(jx,jy,iz)+f(ix,iy,kz)+f(jx,iy,kz)+f(ix,jy,kz)+f(jx,jy,kz) >= 0) then
                                    i1 = min(i1,i2)
                                    indx_tmp(ix,iy,iz) = i1
                                    indx_tmp(jx,jy,kz) = i1
                                    ok = .false.
                                end if
                            end if             
                                                  
                                      
                            
                        
                        end do
                    end do
                end do
                indx = indx_tmp
                if (ok) exit
            end do
            deallocate(indx_tmp)
                        
        !---    final step: renumber the clusters           
        
        !---    count the clusters by looking for different cluster indices  
            allocate(reindx(0:highIndx))
            reindx = -1 ; reindx(highIndx) = 0
            nClust = 0
            do iz = 0,Nz-1
                do iy = 0,Ny-1
                    do ix = 0,Nx-1            
                        i1 = indx(ix,iy,iz)                        
                        if (reindx(i1) == -1) then       !   ... which won't be true if indx(ix,iy,iz) = highIndx
                            nClust = nClust + 1
                            reindx(i1) = nClust                             
                        end if
                    end do
                end do
            end do
            
       !--- now change each cluster to correct index
            do iz = 0,Nz-1
                do iy = 0,Ny-1
                    do ix = 0,Nx-1            
                        i1 = indx(ix,iy,iz) 
                        indx(ix,iy,iz) = reindx(i1)                      
                    end do
                end do
            end do
             
            !if (DBG_TCC) call output()
            if (DBG_TCC) print *,"Lib_topologicallyCorrectClustering::topologicallyCorrectCluster() info - number of unknown,unlit,lit voxels ",count(indx<0),count(indx==0),count(indx>0)
            
            
            
        !---    count the size of each cluster
            if (present(invertLarge)) then
                if (invertLarge) then
                    deallocate(reindx)
                    allocate(reindx(0:nClust))
                    reindx = 0
                    do iz = 0,Nz-1
                        do iy = 0,Ny-1
                            do ix = 0,Nx-1
                                i1 = indx(ix,iy,iz) 
                                reindx(i1) = reindx(i1) + 1
                            end do
                        end do
                    end do
                    do ix = 1,nClust
                        if (reindx(ix)*2>(highIndx-1)) then
                            print *,"Lib_topologicallyCorrectClustering::topologicallyCorrectCluster() info - cluster ",ix," has more than half the nodes ",reindx(ix)," = ",100.0*reindx(ix)/(highIndx-1),"%"
                            call invertBigCluster( indx,ix )
                            exit
                        end if
                    end do
                end if
            end if    
            
            
            return
            
        contains
    !---^^^^^^^^
                
            subroutine output()
        !---^^^^^^^^^^^^^^^^^^^
                integer         ::      ix,iy,iz
                do iz = 0,Nz-1                    
                    do iy = 0,Ny-1
                        do ix = 0,Nx-1
                            if ((indx(ix,iy,iz) == highIndx).or.(indx(ix,iy,iz) == 0)) then
                                write(*,fmt='(a4)',advance="no") " "
                            else
                                write(*,fmt='(i4)',advance="no") indx(ix,iy,iz)
                            end if
                        end do
                        write(*,fmt='(a)',advance="yes") ""
                    end do                                
                    write(*,fmt='(a)',advance="yes") ""
                end do                                
                return
            end subroutine output
                        
        end subroutine topologicallyCorrectCluster                    

        subroutine invertBigCluster( indx,k )
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      cluster k has more than half of the nodes.
    !*      invert it, so that the nodes labelled "0" become k and vice versa.
    !*      to stay topologically correct, add a buffer around each "0"
            integer,dimension(0:,0:,0:),intent(inout)       ::      indx
            integer,intent(in)                              ::      k
            
            integer             ::      Nx,Ny,Nz
            integer             ::      ix,iy,iz
            integer             ::      jx,jy,jz,kx,ky,kz
            integer,dimension(:,:,:),allocatable        ::      indx_tmp
            
            Nx = size(indx,dim=1)
            Ny = size(indx,dim=2)
            Nz = size(indx,dim=3)
            allocate(indx_tmp(0:Nx-1,0:Ny-1,0:Nz-1))
             
        !---    copy the original map, but reset regions labelled "k" to zero.
            indx_tmp(0:Nx-1,0:Ny-1,0:Nz-1) = indx(0:Nx-1,0:Ny-1,0:Nz-1)
            where (indx_tmp == k)
                indx_tmp = 0
            end where
         
            if (DBG_TCC) then
                print *,"Lib_topologicallyCorrectClustering::invertBigCluster() info - original cluster ",0," has ",count(indx==0)," nodes"
                print *,"Lib_topologicallyCorrectClustering::invertBigCluster() info - original cluster ",k," has ",count(indx==k)," nodes"
            end if
            do iz = 0,Nz-1
                do iy = 0,Ny-1
                    do ix = 0,Nx-1
                    
                    !   for each point in the original index, look for regions of type "0". These will be flipped to "k" in the new map                    
                        if (indx(ix,iy,iz) == 0) then
                        
                        !   to preserve the topology, need to flip the 27 cubes surrounding, if possible.
                            do jz = -1,1
                                kz = mod(iz+jz+Nz,Nz)
                                do jy = -1,1
                                    ky = mod(iy+jy+Ny,Ny)
                                    do jx = -1,1
                                        kx = mod(ix+jx+Nx,Nx)  
                                        if ( indx(kx,ky,kz)*(indx(kx,ky,kz)-k) == 0 ) then  !   original map has 0 or k here
                                            indx_tmp(kx,ky,kz) = k       
                                        end if
                                    end do
                                end do
                            end do
                        
                        end if
                    end do
                end do
            end do
            indx(0:Nx-1,0:Ny-1,0:Nz-1) = indx_tmp(0:Nx-1,0:Ny-1,0:Nz-1)
            
            if (DBG_TCC) then
                print *,"Lib_topologicallyCorrectClustering::invertBigCluster() info - new cluster ",0," has ",count(indx==0)," nodes"
                print *,"Lib_topologicallyCorrectClustering::invertBigCluster() info - new cluster ",k," has ",count(indx==k)," nodes"
            end if
            
            return
            
        end subroutine invertBigCluster
        
        
        
        
        subroutine findClusterExtents( indx,nClust, hasx,hasy,hasz )
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^    
    !---    find extent of each cluster - note that hasx is allocated here, and will need deallocating elsewhere.    
            integer,dimension(0:,0:,0:),intent(in)                  ::      indx
            integer,intent(in)                                      ::      nClust
            logical,dimension(:,:),allocatable,intent(out)          ::      hasx,hasy,hasz
            
            integer         ::      ix,iy,iz,kk
            integer         ::      Nx,Ny,Nz
            
            Nx = size(indx,dim=1)
            Ny = size(indx,dim=2)
            Nz = size(indx,dim=3)
            
              
            allocate( hasx(0:Nx-1,0:nClust) ) ; hasx = .false.
            allocate( hasy(0:Ny-1,0:nClust) ) ; hasy = .false.
            allocate( hasz(0:Nz-1,0:nClust) ) ; hasz = .false.
            do iz = 0,Nz-1
                do iy = 0,Ny-1
                    do ix = 0,Nx-1
                        kk = indx(ix,iy,iz)
                        hasx(ix,kk) = .true.
                        hasy(iy,kk) = .true.
                        hasz(iz,kk) = .true.
                    end do
                end do
            end do
            return
        end subroutine findClusterExtents
        
        
        
        
        
        subroutine cropField( f,indx,clust,ok ,hasx,hasy,hasz , fcrop)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ 
    !*      given the phase field crop those regions > 0 in the masked region only
    !*      return ok = false if no volume
    !*      note: allocates fcrop - this will need deallocating elsewhere
    
            real(kind=real64),dimension(0:,0:,0:),intent(in)          ::      f
            integer,dimension(0:,0:,0:),intent(in)      ::      indx
            integer,intent(in)                          ::      clust
            logical,intent(out)                         ::      ok
            
            logical,dimension(0:),intent(in)            ::      hasx,hasy,hasz
            
            real(kind=real64),dimension(:,:,:),allocatable,intent(out)          ::      fcrop
            
            integer             ::      Mx,My,Mz,Mx1,My1,Mz1,Mx2,My2,Mz2 
            integer             ::      ix,iy,iz , jx,jy,jz
                                             
            
            Mx = size(f,dim=1)
            My = size(f,dim=2)
            Mz = size(f,dim=3)

           
            call findCoverage( hasx,mx1,mx2 )
            call findCoverage( hasy,my1,my2 )
            call findCoverage( hasz,mz1,mz2 )
            
            if ( (mx2<mx1).or.(my2<my1).or.(mz2<mz1) ) then
                ok = .false.
                return
            end if
            
            
            
            ok = .true.
            allocate(fcrop(0:mx2-mx1,0:my2-my1,0:mz2-mz1))
            
            do iz = 0,mz2-mz1
                jz = mod( iz + mz1,Mz )
                do iy = 0,my2-my1
                    jy = mod( iy + my1,My )
                    do ix = 0,mx2-mx1
                        jx = mod( ix + mx1,Mx )
                        if ( indx(jx,jy,jz)*(indx(jx,jy,jz)-clust)==0 ) then  
                            !   inside requested cluster, or outside all clusters
                            fcrop(ix,iy,iz) = f(jx,jy,jz)
                        else
                            !   inside _the wrong cluster_
                            fcrop(ix,iy,iz) = - 0.01d0
                        end if
                    end do
                end do
            end do
            return
            
            
        contains
    !---^^^^^^^^
            
            subroutine findCoverage( has,m1,m2 )
        !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
        !*      the array is true from from mod( m1,Mx ) to mod( m2,Mx ) with m2>m1
        !*      m1 and m2 include a buffer of one where necessary
        
        
                logical,dimension(0:),intent(in)        ::      has
                integer,intent(out)                     ::      m1,m2
                
                integer         ::      ii,mm,nn
                
                mm = size(has)
                nn = count(has)
                
                if (nn >= mm-2) then
                    !   lit pretty much all the way across. No advantage in a buffer.
                    m1 = 0 ; m2 = mm-1
                else if (has(0)) then
                    !   lit on the left boundary. Does it stretch across the periodic repeat?
                    if (has(mm-1)) then
                        !   starts midway through cell and crosses boundary
                        do ii = mm-1,1,-1
                            if (has(ii)) then
                                m1 = ii-1 
                            else
                                exit
                            end if
                        end do
                        do ii = m1+1,m1+mm-1                            
                            if (has( mod( ii,mm ) )) then
                                m2 = ii+1 
                            else
                                exit
                            end if
                        end do
                    else
                        !   starts at 0 and ends midway through cell. Because of pbc, have to start on other side.
                        m1 = mm-1
                        do ii = m1+1,m1+mm-1
                            if (has( mod( ii,mm ) )) then
                                m2 = ii+1 
                            else
                                exit
                            end if
                        end do
                    end if
                else if (has(mm-1)) then
                    !   starts midway through cell and ends at mm-1
                    do ii = 0,mm-1
                        if (has(ii)) then
                            m1 = ii-1 ; exit
                        end if
                    end do
                    m2 = mm   
                else if (nn==0) then   
                !   unlit all the way across???
                    m1 = mm-1 ; m2 = 0                                                        
                else 
                    !   starts and ends midway through cell
                    do ii = 1,mm-1
                        if (has(ii)) then
                            m1 = ii-1 ; exit
                        end if
                    end do
                    do ii = mm-2,m1,-1
                        if (has(ii)) then
                            m2 = ii+1 ; exit
                        end if
                    end do                
                end if
                                
                return
            end subroutine findCoverage
                       
        end subroutine cropField
                  

    end module Lib_TopologicallyCorrectClustering   