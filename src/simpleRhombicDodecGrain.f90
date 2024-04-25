
    !   gfortran -ffree-line-length-256 src/simpleRhombicDodecGrain.f90 -o bin/simpleRhombicDodecGrain.exe

    program simpleRhombicDodecGrain
!---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

        use iso_fortran_env
        implicit none
        
        
        
        integer                     ::      N       !   number of atoms to place
        integer                     ::      Nx      !   number of fcc cells per side
        real(kind=real64)           ::      a       !   supercell side length
        real(kind=real64)           ::      eps     !   grain boundary parameter
        character(len=256)          ::      filename
        
        
        
        
        
        real(kind=real64),dimension(:,:),allocatable        ::      x_seed
        real(kind=real64),dimension(:,:),allocatable        ::      x_atom
        integer,dimension(:),allocatable                    ::      g_atom
        
        real(kind=real64),dimension(3,4),parameter  ::      fcc_motif = reshape( (/     0.25,0.25,0.25,                 &
                                                                                        0.75,0.75,0.25,                 &
                                                                                        0.75,0.25,0.75,                 &
                                                                                        0.25,0.75,0.75      /),(/3,4/) )
        
        
        integer             ::      ii,ik,ix,iy,iz
        integer             ::      jj,j1,j2,kk
        integer             ::      M                   !   number of seeds
        real(kind=real64)   ::      vv,dd,d1,d2
        real(kind=real64),dimension(3)  ::      dx       
        character(len=256)          ::      dummy
        
        print *,"usage: simpleRhombicDodecGrain.exe N Nx eps"
        print *,"N = number of atoms"
        print *,"with Nx = number of fcc cell sides ( M = 4 Nx^3 )"
        print *,"and eps = distance parameter for grain boundary"
        print *,""
        
        call get_command_argument( 1,dummy)
        read(dummy,fmt=*) N ; print *,"read   N = ",N
        call get_command_argument( 2,dummy)
        read(dummy,fmt=*) Nx ; print *,"read  Nx = ",Nx
        call get_command_argument( 3,dummy)
        read(dummy,fmt=*) eps ; print *,"read eps = ",eps
        
        M = Nx*Nx*Nx*4      
       ! eps = 0.01
       ! N = 10000
        allocate(x_seed(3,M))
        allocate(x_atom(3,N))
        allocate(g_atom(N))
        
        
    !---    place seeds in fcc lattice
        ii = 0
        do iz = 0,Nx-1
            do iy = 0,Nx-1
                do ix = 0,Nx-1
                    do ik = 1,4
                        ii = ii + 1
                        x_seed(:,ii) = (/ix,iy,iz/) + fcc_motif(:,ik)
                    end do
                end do
            end do
        end do
        x_seed = x_seed / Nx
       
    !---    add atoms
        call random_number(x_atom)
        g_atom = 0
        
    !---    which grain?
        do ii = 1,N
            d1 = huge(1.0) ; j1 = 0         !   nearest neighbour
            d2 = huge(1.0) ; j2 = 0         !   second nearest
            do jj = 1,M
                dx = x_atom(:,ii) - x_seed(:,jj)
                do kk = 1,3
                    if (2*dx(kk)<-1) then
                        dx(kk) = dx(kk) + 1
                    else if (2*dx(kk)>1) then
                        dx(kk) = dx(kk) - 1
                    end if
                end do
                dd = dx(1)*dx(1) + dx(2)*dx(2) + dx(3)*dx(3) 
                if (dd < d1) then
                    d2 = d1 ; j2 = j1
                    d1 = dd ; j1 = jj
                else if (dd < d2) then
                    d2 = dd ; j2 = jj
                end if
            end do
            
            
            if (d1<(1-eps)*(1-eps)*d2) then
                !   closest is much closer than second closest
                g_atom(ii) = j1
            end if
        end do
                    
        do jj = 0,M
            print *,"grain ",jj," count ",count(g_atom(:) == jj)
        end do 
        
    !---    scale lattice to have rhomb dodec side a0 = 1
        ! vol = 16/9 sqrt(3) a0^3 = a^3/M
        !   a = a0 ( 16 M sqrt(3) / 9 )^1/3
        a = ( 16 * M * sqrt(3.0d0) / 9 )**(1/3.0)        
        x_seed = x_seed * a
        x_atom = x_atom * a
        
        write(dummy,fmt='(i7)') N ; dummy = adjustl(dummy) ; dummy = repeat("0",7-len_trim(dummy))//trim(dummy)
        filename = "test/dodec."//trim(dummy)
        write(dummy,fmt='(i5)') M ; dummy = adjustl(dummy) ; dummy = repeat("0",5-len_trim(dummy))//trim(dummy)
        filename = trim(filename)//"."//trim(dummy)
        write(dummy,fmt='(f5.3)') eps ; dummy = adjustl(dummy) ; dummy = repeat("0",5-len_trim(dummy))//trim(dummy)
        filename = trim(filename)//"."//trim(dummy)
        
        open(unit=850,file=trim(filename)//".tmp",action="write")
            write(850,fmt='(i8)') N
            write(850,fmt='(4(a,f12.8))') "Lattice=""",a," 0.0 0.0 0.0 ",a," 0.0 0.0 0.0 ",a,""" Properties=species:S:1:pos:R:3:grain:I:1"
            do ii = 1,N
                write(850,fmt='(a,3f12.8,i6)') "Du ",x_atom(:,ii),g_atom(ii)
            end do
        close(unit=850)  
        call system("mv "//trim(filename)//".tmp "//trim(filename)//".xyz")
             
        print *,""
        print *,"done"
        print *,""
        
    end program simpleRhombicDodecGrain
        
        
        