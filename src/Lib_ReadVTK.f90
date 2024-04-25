
    module Lib_ReadVTK
!---^^^^^^^^^^^^^^^^^^
!*      read/write vtk files
!*      as format 
!*          # vtk DataFile Version 3.0
!*          Saved using Lib_ReadVTK
!*          ASCII
!*          
!*          DATASET RECTILINEAR_GRID
!*          DIMENSIONS  nx ny nz                             !  number of NODE points along x,y,z directions
!*          X_COORDINATES nx  float
!*          real(xx(1:nx),kind=real32)                       !  positions of NODE points 
!*          Y_COORDINATES ny float                           !  the nodes are at the centre of the voxel cubes
!*          real(yy(1:ny),kind=real32)                       !  so note that the dimension of the supercell extends outside the node points
!*          Z_COORDINATES nz float
!*          real(zz(1:nz),kind=real32)
!*          POINT_DATA nx*ny*nz
!*          SCALARS miscellaneous float
!*          LOOKUP_TABLE default
!*          real(dat(:,:,:),kind=real32)                     !  dat(1,1,1) dat(2,1,1) dat(3,1,1) ... dat(1,2,1) dat(2,2,1) ... dat(nx,ny,nz)



!*          
    
        use iso_fortran_env
        use NBAX_StringTokenizers
        implicit none
        private
        
    !---
    
        public      ::      readVTKheader
        public      ::      readVTK
        public      ::      readCSVheader
        public      ::      readCSV
        public      ::      writeVTK
        public      ::      writeChgcar
    
    !---
    
        interface   readVTK
            module procedure    readVTK0
            module procedure    readVTK1
        end interface
        
        interface   writeVTK
            module procedure    writeVTK0
            module procedure    writeVTK1
        end interface
        
        interface   report
            module procedure    report0
            module procedure    report1
        end interface
        
        
        interface   writeChgcar
            module procedure    writeChgcar1
            module procedure    writeChgcar2
            module procedure    writeChgcar2_32
        end interface
    
    contains
!---^^^^^^^^

        subroutine readVTKheader( filename,nx,ny,nz )
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^           
            character(len=*),intent(in)     ::      filename
            integer,intent(out)             ::      nx,ny,nz
            character(len=256)      ::      dummy
            integer                 ::      ii
            
            open(unit=700,file=trim(filename),action="read")
                do ii = 1,5
                    read(700,fmt='(a)') dummy
                end do
                read(700,fmt=*) dummy,nx,ny,nz
            close(unit=700)
            print *,"Lib_ReadVTK::readVTKheader info - nx,ny,nz = ",nx,ny,nz
            return
         end subroutine readVTKheader
            
            
        subroutine readCSVheader( filename,nx,ny,nz )
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^       
    !*      figure out from the CSV file what the number of voxels are
            character(len=*),intent(in)     ::      filename
            integer,intent(out)             ::      nx,ny,nz
            character(len=256)              ::      dummy
            integer                         ::      ii,nn
            real(kind=real64),dimension(3)  ::      xx
            real(kind=real64),dimension(1000,3) ::      xknots
            integer,dimension(3)                ::      nknots
            logical                             ::      ok
            
!            print *,"readCSVheader"
            open(unit=700,file=trim(filename),action="read")
                read(700,fmt='(a)') dummy
                nKnots = 0
                xKnots = 0.0d0
                nn = 0
                do 
                    read(700,fmt=*,iostat=ii) xx(1:3)
                    if (ii/=0) exit
                    nn = nn + 1
                    do ii = 1,3
                        ok = any( abs(xx(ii)-xKnots(1:nKnots(ii),ii))<1.0d-6 )
                        if (.not. ok) then
                            nKnots(ii) = nKnots(ii) + 1
                            xKnots(nKnots(ii),ii) = xx(ii)
!                           print *,"new knot ",ii,nKnots(ii),xx(ii)
                        end if                              
                    end do
                end do                  
            close(unit=700)
            nx = nKnots(1)
            ny = nKnots(2)
            nz = nKnots(3)
            
            print *,"Lib_ReadVTK::readCSVheader info - nx,ny,nz = ",nx,ny,nz," lines = ",nn
            return
         end subroutine readCSVheader
            

            
         subroutine report0( filename,nx,ny,nz,xmin,xmax, f )
     !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
     !*      note that xmin,xmax are the CELL boundaries, not the limits of the knot points.
            character(len=*),intent(in)                     ::      filename
            integer,intent(in)                              ::      nx,ny,nz 
            real(kind=real64),dimension(3),intent(in)       ::      xmin,xmax
            real(kind=real64),dimension(:,:,:),intent(in)   ::  f
            real(kind=real64)       ::      fsum
            integer                 ::      nnonzero
            real(kind=real64),dimension(3)                  ::      deltax
            print *,"vtk file "//trim(filename)
            deltax(1) = (xmax(1)-xmin(1))/nx
            deltax(2) = (xmax(2)-xmin(2))/ny
            deltax(3) = (xmax(3)-xmin(3))/nz
            write(*,fmt='(a,f16.8,a,f16.8,a,i8,a,f16.8)') "  knots x : ",xmin(1)+deltax(1)/2,":",xmax(1)-deltax(1)/2," nx = ",nx," deltax = ",deltax(1)
            write(*,fmt='(a,f16.8,a,f16.8,a,i8,a,f16.8)') "  knots y : ",xmin(2)+deltax(2)/2,":",xmax(2)-deltax(2)/2," ny = ",ny," deltay = ",deltax(2)
            write(*,fmt='(a,f16.8,a,f16.8,a,i8,a,f16.8)') "  knots z : ",xmin(3)+deltax(3)/2,":",xmax(3)-deltax(3)/2," nz = ",nz," deltaz = ",deltax(3)
            write(*,fmt='(a,f16.8,a,f16.8,a,i8,a,f16.8)') "  cell x  : ",xmin(1)            ,":",xmax(1)            
            write(*,fmt='(a,f16.8,a,f16.8,a,i8,a,f16.8)') "  cell y  : ",xmin(2)            ,":",xmax(2)            
            write(*,fmt='(a,f16.8,a,f16.8,a,i8,a,f16.8)') "  cell z  : ",xmin(3)            ,":",xmax(3)            
            write(*,fmt='(a,f16.8,a,f16.8)')              "  f       : ",minval(f)          ,":",maxval(f)
            fsum = real(sum(f),kind=real64)
            nnonzero = count(f>0)
            write(*,fmt='(a,f16.4)') " sum  = ",fsum                       
            write(*,fmt='(a,f16.4)') " mean = ",fsum/(nx*ny*nz)
            write(*,fmt='(a,f16.4)') " mean (non-zero cells) = ",fsum/max(nnonzero,1)
            write(*,fmt='(a,i16)  ') " number non zero cells = ",nnonzero
            print *,""     
            return
        end subroutine report0
         
            
         subroutine report1( filename,nx,ny,nz,xmin,xmax, f )
     !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
     !*      note that xmin,xmax are the CELL boundaries, not the limits of the knot points.
            character(len=*),intent(in)                     ::      filename
            integer,intent(in)                              ::      nx,ny,nz 
            real(kind=real64),dimension(3),intent(in)       ::      xmin,xmax
            integer,dimension(:,:,:),intent(in)             ::      f
            real(kind=real64)       ::      fsum
            integer                 ::      nnonzero
            real(kind=real64),dimension(3)                  ::      deltax
            print *,"vtk file "//trim(filename)
            deltax(1) = (xmax(1)-xmin(1))/nx
            deltax(2) = (xmax(2)-xmin(2))/ny
            deltax(3) = (xmax(3)-xmin(3))/nz
            write(*,fmt='(a,f16.8,a,f16.8,a,i8,a,f16.8)') "  knots x : ",xmin(1)+deltax(1)/2,":",xmax(1)-deltax(1)/2," nx = ",nx," deltax = ",deltax(1)
            write(*,fmt='(a,f16.8,a,f16.8,a,i8,a,f16.8)') "  knots y : ",xmin(2)+deltax(2)/2,":",xmax(2)-deltax(2)/2," ny = ",ny," deltay = ",deltax(2)
            write(*,fmt='(a,f16.8,a,f16.8,a,i8,a,f16.8)') "  knots z : ",xmin(3)+deltax(3)/2,":",xmax(3)-deltax(3)/2," nz = ",nz," deltaz = ",deltax(3)
            write(*,fmt='(a,f16.8,a,f16.8,a,i8,a,f16.8)') "  cell x  : ",xmin(1)            ,":",xmax(1)            
            write(*,fmt='(a,f16.8,a,f16.8,a,i8,a,f16.8)') "  cell y  : ",xmin(2)            ,":",xmax(2)            
            write(*,fmt='(a,f16.8,a,f16.8,a,i8,a,f16.8)') "  cell z  : ",xmin(3)            ,":",xmax(3)            
            write(*,fmt='(a,i16,a,i16)')                  "  f       : ",minval(f)          ,":",maxval(f)
            fsum = real(sum(f),kind=real64)
            nnonzero = count(f>0)
            write(*,fmt='(a,f16.4)') " sum  = ",fsum
            write(*,fmt='(a,f16.4)') " mean = ",fsum/(nx*ny*nz)
            write(*,fmt='(a,f16.4)') " mean = ",fsum/(nx*ny*nz)
            write(*,fmt='(a,f16.4)') " mean (non-zero cells) = ",fsum/max(nnonzero,1)
            write(*,fmt='(a,i16)  ') " number non zero cells = ",nnonzero
            print *,""     
            return
        end subroutine report1
         
    !---        
        
        subroutine readCSV( filename,nx,ny,nz,xmin,xmax, f , column)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^      
    !*      read csv data in format x,y,z,c1,c2,c3...
    !*      note that xmin,xmax are the CELL boundaries, not the limits of the knot points.
            character(len=*),intent(in)                 ::      filename
            integer,intent(in)                          ::      nx,ny,nz 
            real(kind=real64),dimension(3),intent(out)  ::      xmin,xmax
            real(kind=real64),dimension(:,:,:),intent(inout)    ::  f
            integer,intent(in),optional                         ::  column
            
            real(kind=real64),dimension(20)     ::      xx
            
            character(len=256)      ::       dummy
            integer                 ::      ii,jj,col,kk
            logical                 ::      ok
            
            real(kind=real64),dimension(1000,3) ::      xknots
            integer,dimension(3)                ::      nknots,ll
            real(kind=real64),dimension(3)      ::      deltax
            xmin = 0
            xmax = 0
            col = 1 ; if (present(column)) col = column
            nKnots = 0
            xKnots = 0.0d0
            open(unit=700,file=trim(filename),action="read")
                read(700,fmt='(a)') dummy
                do ii = 1,nx*ny*nz
                    read(700,fmt=*) xx(1:3+col)
                    do jj = 1,3
                        ok = .false.
                        do kk = 1,nKnots(jj)
                            if ( abs(xx(jj)-xKnots(kk,jj))<1.0d-6 ) then
                                ll(jj) = kk
                                ok = .true.
                                exit
                            end if
                        end do
                        if (.not. ok) then
                            nKnots(jj) = nKnots(jj) + 1
                            xKnots(nKnots(jj),jj) = xx(jj)
                            ll(jj) = nKnots(jj)                     
                        end if                              
                    end do
                    f( ll(1),ll(2),ll(3) ) = xx(3+col)
                end do
            close(unit=700)
            
        !---    find extent of the NODE positions            
            xmin = minval(xKnots)            
            xmax = maxval(xKnots)            
            deltax(1) = ( xmax(1) - xmin(1) )/(nx-1)
            deltax(2) = ( xmax(2) - xmin(2) )/(ny-1)
            deltax(3) = ( xmax(3) - xmin(3) )/(nz-1)
        !----   find extent of the CELL             
            xmin(1:3) = xmin(1:3) - deltax(1:3)/2         
            xmax(1:3) = xmax(1:3) + deltax(1:3)/2        
            
                
            call report( filename,nx,ny,nz,xmin,xmax, f )
           
            
            
            
            return
         end subroutine readCSV
         
    
    !---         
         
        subroutine readVTK0( filename,nx,ny,nz,xmin,xmax, f )
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*       
    !*      note that xmin,xmax are the CELL boundaries, not the limits of the knot points.    
            character(len=*),intent(in)                 ::      filename
            integer,intent(in)                          ::      nx,ny,nz 
            real(kind=real64),dimension(3),intent(out)  ::      xmin,xmax
            real(kind=real64),dimension(:,:,:),intent(inout)    ::  f
            real(kind=real64),dimension(max(max(nx,ny),nz))     ::  xx
            character(len=256)      ::       dummy
            integer                 ::      ii,ix,iy,iz
            real(kind=real64),dimension(3)      ::      deltax
!            print *,"Lib_ReadVTK::readVTK info - nx,ny,nz = ",nx,ny,nz
            xmin = 0
            xmax = 0
            open(unit=700,file=trim(filename),action="read")
                do ii = 1,6
                    read(700,fmt='(a)') dummy
                end do
                read(700,fmt='(a)') dummy
                read(700,fmt=*) xx(1:nx) 
                xmin(1) = minval( xx(1:nx) )
                xmax(1) = maxval( xx(1:nx) )                
                read(700,fmt='(a)') dummy
                read(700,fmt=*) xx(1:ny)
                xmin(2) = minval( xx(1:ny) )
                xmax(2) = maxval( xx(1:ny) )                
                read(700,fmt='(a)') dummy
                read(700,fmt=*) xx(1:nz)
                xmin(3) = minval( xx(1:nz) )
                xmax(3) = maxval( xx(1:nz) )                
                read(700,fmt='(a)') dummy
                read(700,fmt='(a)') dummy
                read(700,fmt='(a)') dummy
                do iz = 1,nz
                    do iy = 1,ny
                        do ix = 1,nx
                            read(700,fmt=*) f(ix,iy,iz)
                        end do
                    end do
                end do
            close(unit=700)
            
        !---    find extent of the NODE positions            
            deltax(1) = ( xmax(1) - xmin(1) )/(nx-1)
            deltax(2) = ( xmax(2) - xmin(2) )/(ny-1)
            deltax(3) = ( xmax(3) - xmin(3) )/(nz-1)
        !----   find extent of the CELL             
            xmin(1:3) = xmin(1:3) - deltax(1:3)/2         
            xmax(1:3) = xmax(1:3) + deltax(1:3)/2        
            
                
            call report( filename,nx,ny,nz,xmin,xmax, f )
           
            
            
            
            return
         end subroutine readVTK0
         

    
        subroutine readVTK1( filename,nx,ny,nz,xmin,xmax, f )
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      note that xmin,xmax are the CELL boundaries, not the limits of the knot points.
            character(len=*),intent(in)                 ::      filename
            integer,intent(in)                          ::      nx,ny,nz 
            real(kind=real64),dimension(3),intent(out)  ::      xmin,xmax
            integer,dimension(:,:,:),intent(inout)      ::  f
            real(kind=real64),dimension(max(max(nx,ny),nz))     ::  xx
            
            character(len=256)      ::       dummy
            integer                 ::      ii,ix,iy,iz
            real(kind=real64),dimension(3)      ::      deltax
!            print *,"Lib_ReadVTK::readVTK info - nx,ny,nz = ",nx,ny,nz
            xmin = 0
            xmax = 0
            open(unit=700,file=trim(filename),action="read")
                do ii = 1,6
                    read(700,fmt='(a)') dummy
                end do
                read(700,fmt='(a)') dummy
                read(700,fmt=*) xx(1:nx) 
                xmin(1) = minval( xx(1:nx) )
                xmax(1) = maxval( xx(1:nx) )                
                read(700,fmt='(a)') dummy
                read(700,fmt=*) xx(1:ny)
                xmin(2) = minval( xx(1:ny) )
                xmax(2) = maxval( xx(1:ny) )                
                read(700,fmt='(a)') dummy
                read(700,fmt=*) xx(1:nz)
                xmin(3) = minval( xx(1:nz) )
                xmax(3) = maxval( xx(1:nz) )                
                read(700,fmt='(a)') dummy
                read(700,fmt='(a)') dummy
                read(700,fmt='(a)') dummy
                do iz = 1,nz
                    do iy = 1,ny
                        do ix = 1,nx
                            read(700,fmt=*) f(ix,iy,iz)
                        end do
                    end do
                end do
            close(unit=700)
                
        !---    find extent of the NODE positions            
            deltax(1) = ( xmax(1) - xmin(1) )/(nx-1)
            deltax(2) = ( xmax(2) - xmin(2) )/(ny-1)
            deltax(3) = ( xmax(3) - xmin(3) )/(nz-1)
        !----   find extent of the CELL             
            xmin(1:3) = xmin(1:3) - deltax(1:3)/2         
            xmax(1:3) = xmax(1:3) + deltax(1:3)/2        
            
            
            call report( filename,nx,ny,nz,xmin,xmax, f )
           
            
            
            
            return
         end subroutine readVTK1         
         
         
         
    
        subroutine writeVTK0( filename,nx,ny,nz,xmin,xmax, f )
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      note that xmin,xmax are the CELL boundaries, not the limits of the knot points.
            character(len=*),intent(in)                         ::      filename
            integer,intent(in)                                  ::      nx,ny,nz 
            real(kind=real64),dimension(3),intent(in)           ::      xmin,xmax
            real(kind=real64),dimension(:,:,:),intent(in)       ::  f
            real(kind=real64),dimension(max(max(nx,ny),nz))     ::  xx

            integer                 ::      ii,ix,iy,iz
            real(kind=real64),dimension(3)      ::      deltax

           
            
            
            
            call report( filename,nx,ny,nz,xmin,xmax, f )
           

        !---    find extent of the CELL positions            
            deltax(1) = ( xmax(1) - xmin(1) )/nx
            deltax(2) = ( xmax(2) - xmin(2) )/ny
            deltax(3) = ( xmax(3) - xmin(3) )/nz            
              
                        
            open(unit=700,file=trim(filename),action="write")
                write (700,fmt='(a)') "# vtk DataFile Version 3.0"
                write (700,fmt='(a)') "Saved using Lib_ReadVTK"
                write (700,fmt='(a)') "ASCII"
                write (700,fmt='(a)') ""
                write (700,fmt='(a)') "DATASET RECTILINEAR_GRID"
                write (700,fmt='(a,3i6)') "DIMENSIONS ",nx,ny,nz
                
                write (700,fmt='(a,i6,a)') "X_COORDINATES ",nx," float"
                do ii = 1,nx
                    xx(ii) = xmin(1) + (ii-0.5d0)*deltax(1)
                end do              
                write (700,fmt=*) real(xx(1:nx),kind=real32)
                write (700,fmt='(a,i6,a)') "Y_COORDINATES ",ny," float"
                do ii = 1,ny
                    xx(ii) = xmin(2) + (ii-0.5d0)*deltax(2)
                end do              
                write (700,fmt=*) real(xx(1:ny),kind=real32)
                write (700,fmt='(a,i6,a)') "Z_COORDINATES ",nz," float"
                do ii = 1,nz
                    xx(ii) = xmin(3)  + (ii-0.5d0)*deltax(3)
                end do              
                write (700,fmt=*) real(xx(1:nz),kind=real32)
                
                write (700,fmt='(a,i12)') "POINT_DATA ",nx*ny*nz
                write (700,fmt='(a,i12)') "SCALARS miscellaneous float"
                write (700,fmt='(a,i12)') "LOOKUP_TABLE default"

                do iz = 1,nz
                    do iy = 1,ny
                        do ix = 1,nx
                            write(700,fmt=*) real(f(ix,iy,iz),kind=real32)
                        end do
                    end do
                end do
            close(unit=700)
                
           
            
            
            
            return
            
            
            
            
         end subroutine writeVTK0
         
    
        subroutine writeVTK1( filename,nx,ny,nz,xmin,xmax, f )
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      note that xmin,xmax are the CELL boundaries, not the limits of the knot points.
            character(len=*),intent(in)                         ::      filename
            integer,intent(in)                                  ::      nx,ny,nz 
            real(kind=real64),dimension(3),intent(in)           ::      xmin,xmax
            integer,dimension(:,:,:),intent(in)                 ::      f
            real(kind=real64),dimension(max(max(nx,ny),nz))     ::      xx
            

            integer                 ::      ii,ix,iy,iz
            real(kind=real64),dimension(3)      ::      deltax
!            print *,"Lib_ReadVTK::writeVTK info - nx,ny,nz = ",nx,ny,nz

            call report( filename,nx,ny,nz,xmin,xmax, f )
           
            
            
        !---    find extent of the CELL positions            
            deltax(1) = ( xmax(1) - xmin(1) )/nx
            deltax(2) = ( xmax(2) - xmin(2) )/ny
            deltax(3) = ( xmax(3) - xmin(3) )/nz            
              
                        
            open(unit=700,file=trim(filename),action="write")
                write (700,fmt='(a)') "# vtk DataFile Version 3.0"
                write (700,fmt='(a)') "Saved using Lib_ReadVTK"
                write (700,fmt='(a)') "ASCII"
                write (700,fmt='(a)') ""
                write (700,fmt='(a)') "DATASET RECTILINEAR_GRID"
                write (700,fmt='(a,3i6)') "DIMENSIONS ",nx,ny,nz
                
                write (700,fmt='(a,i6,a)') "X_COORDINATES ",nx," float"
                do ii = 1,nx
                    xx(ii) = xmin(1) + (ii-0.5d0)*deltax(1)
                end do              
                write (700,fmt=*) real(xx(1:nx),kind=real32)
                write (700,fmt='(a,i6,a)') "Y_COORDINATES ",ny," float"
                do ii = 1,ny
                    xx(ii) = xmin(2) + (ii-0.5d0)*deltax(2)
                end do              
                write (700,fmt=*) real(xx(1:ny),kind=real32)
                write (700,fmt='(a,i6,a)') "Z_COORDINATES ",nz," float"
                do ii = 1,nz
                    xx(ii) = xmin(3) + (ii-0.5d0)*deltax(3)
                end do              
                write (700,fmt=*) real(xx(1:nz),kind=real32)
                
                write (700,fmt='(a,i12)') "POINT_DATA ",nx*ny*nz
                write (700,fmt='(a,i12)') "SCALARS miscellaneous int"
                write (700,fmt='(a,i12)') "LOOKUP_TABLE default"

                do iz = 1,nz
                    do iy = 1,ny
                        do ix = 1,nx
                            write(700,fmt=*) f(ix,iy,iz)
                        end do
                    end do
                end do
            close(unit=700)
                
           
            
            
            
            return
         end subroutine writeVTK1         
            


        subroutine writeChgcar1( filename,nx,ny,nz,xmin,xmax, f )
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      note that xmin,xmax are the CELL boundaries, not the limits of the knot points.
            character(len=*),intent(in)                         ::      filename
            integer,intent(in)                                  ::      nx,ny,nz 
            real(kind=real64),dimension(3),intent(in)           ::      xmin,xmax
            real(kind=real64),dimension(:,:,:),intent(in)       ::  f
            real(kind=real64)       ::      dd
            
            integer                 ::      ix,iy,iz
            
            character(len=256)                              ::      dummy                         
            integer                                         ::      nn 
            real(kind=real64)                               ::      datmax 
            
            datmax = maxval(abs(f))
            nn = ceiling(log10(max(1.0d0,datmax)))
                                 
            
            
            call report( filename,nx,ny,nz,xmin,xmax, f )
            open(unit=700,file=trim(filename)//".tmp",action="write")
            write (unit=700,fmt='(a)') "Lib_ReadVTK::writeChgcar()"
            write (unit=700,fmt='(f16.8)') xmax(1)-xmin(1)
            write (unit=700,fmt='(3f12.5)') 1.0d0,0.0d0,0.0d0
            write (unit=700,fmt='(3f12.5)') 0.0d0,(xmax(2)-xmin(2))/(xmax(1)-xmin(1)),0.0d0
            write (unit=700,fmt='(3f12.5)') 0.0d0,0.0d0,(xmax(3)-xmin(3))/(xmax(1)-xmin(1))
            write (unit=700,fmt='(a)')   "Du"
            write (unit=700,fmt='(a)')   "1"
            write (unit=700,fmt='(a)')   "direct"
            write (unit=700,fmt='(a)')   "  0.000000  0.000000  0.000000"
            write (unit=700,fmt='(a)')
            write (unit=700,fmt='(3i5)')   nx,ny,nz
            dd = (xmax(1)-xmin(1))*(xmax(2)-xmin(2))*(xmax(3)-xmin(3))
            
            datmax = maxval( abs( f(:,:,:) ) )* dd
            nn = ceiling(log10(max(1.0001d0,datmax)))
            
            do iz = 1,nz
                do iy = 1,ny
                    do ix = 1,nx,12
                        dummy = opLine(f(ix:min(ix+11,nx),iy,iz) * dd,nn)
                        write (unit=700,fmt='(a)') trim(dummy)
                    end do
                end do
            end do
           
            !write (unit=700,fmt=*)   f(:,:,:) * dd

            close(unit=700)
            call system( "mv "//trim(filename)//".tmp "//trim(filename) )
            return            
            
            
        end subroutine writeChgcar1

        subroutine writeChgcar2( filename,nx,ny,nz,a, f )
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ 
    !*      note that xmin,xmax are the CELL boundaries, not the limits of the knot points.
            character(len=*),intent(in)                         ::      filename
            integer,intent(in)                                  ::      nx,ny,nz 
            real(kind=real64),dimension(3,3),intent(in)         ::      a
            real(kind=real64),dimension(:,:,:),intent(in)       ::      f
            real(kind=real64)       ::      dd
            
            integer                 ::      ix,iy,iz
            
            character(len=256)                              ::      dummy                         
            integer                                         ::      nn 
            real(kind=real64)                               ::      datmax , vol
            
            datmax = maxval(abs(f)) 
            nn = ceiling(log10(max(1.0d0,datmax)))
                                 
            
            
            !call report( filename,nx,ny,nz,xmin,xmax, f )
            open(unit=700,file=trim(filename)//".tmp",action="write")
            !write (unit=700,fmt='(a)') "Lib_ReadVTK::writeChgcar()"
            !write (unit=700,fmt='(f16.8)') xmax(1)-xmin(1)
            !write (unit=700,fmt='(3f12.5)') 1.0d0,0.0d0,0.0d0
            !write (unit=700,fmt='(3f12.5)') 0.0d0,(xmax(2)-xmin(2))/(xmax(1)-xmin(1)),0.0d0
            !write (unit=700,fmt='(3f12.5)') 0.0d0,0.0d0,(xmax(3)-xmin(3))/(xmax(1)-xmin(1))
            
            vol = determinant3Mat(a)
            dd = vol**(1.0d0/3)
            
            
            write (unit=700,fmt='(a)') "Lib_ReadVTK::writeChgcar()"
            write (unit=700,fmt='(f16.8)') dd
            write (unit=700,fmt='(3f12.5)') a(1:3,1)/dd
            write (unit=700,fmt='(3f12.5)') a(1:3,2)/dd
            write (unit=700,fmt='(3f12.5)') a(1:3,3)/dd
            
           
            
            
            write (unit=700,fmt='(a)')   "Du"
            write (unit=700,fmt='(a)')   "1"
            write (unit=700,fmt='(a)')   "direct"
            write (unit=700,fmt='(a)')   "  0.000000  0.000000  0.000000"
            write (unit=700,fmt='(a)')
            write (unit=700,fmt='(3i5)')   nx,ny,nz
            
            
            
           ! dd = (xmax(1)-xmin(1))*(xmax(2)-xmin(2))*(xmax(3)-xmin(3))
            
            datmax = maxval( abs( f(:,:,:) ) )* vol
            nn = ceiling(log10(max(1.0001d0,datmax)))
            
            do iz = 1,nz
                do iy = 1,ny
                    do ix = 1,nx,12
                        dummy = opLine(f(ix:min(ix+11,nx),iy,iz) * vol,nn)
                        write (unit=700,fmt='(a)') trim(dummy)
                    end do
                end do
            end do
           
            !write (unit=700,fmt=*)   f(:,:,:) * dd

            close(unit=700)
            call system( "mv "//trim(filename)//".tmp "//trim(filename) )
            return            
            
            
        end subroutine writeChgcar2

         
        subroutine writeChgcar2_32( filename,nx,ny,nz,a, f )
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      note that xmin,xmax are the CELL boundaries, not the limits of the knot points.
            character(len=*),intent(in)                         ::      filename
            integer,intent(in)                                  ::      nx,ny,nz 
            real(kind=real64),dimension(3,3),intent(in)         ::      a
            real(kind=real32),dimension(:,:,:),intent(in)       ::      f
            real(kind=real64)       ::      dd 
            
            integer                 ::      ix,iy,iz
            
            character(len=256)                              ::      dummy                         
            integer                                         ::      nn 
            real(kind=real64)                               ::      datmax , vol
            
            datmax = maxval(abs(f)) 
            nn = ceiling(log10(max(1.0d0,datmax)))
                                 
            
            
            !call report( filename,nx,ny,nz,xmin,xmax, f )
            open(unit=700,file=trim(filename)//".tmp",action="write")
            !write (unit=700,fmt='(a)') "Lib_ReadVTK::writeChgcar()"
            !write (unit=700,fmt='(f16.8)') xmax(1)-xmin(1)
            !write (unit=700,fmt='(3f12.5)') 1.0d0,0.0d0,0.0d0
            !write (unit=700,fmt='(3f12.5)') 0.0d0,(xmax(2)-xmin(2))/(xmax(1)-xmin(1)),0.0d0
            !write (unit=700,fmt='(3f12.5)') 0.0d0,0.0d0,(xmax(3)-xmin(3))/(xmax(1)-xmin(1))
            
            vol = determinant3Mat(a)
            dd = vol**(1.0d0/3)
            
            
            write (unit=700,fmt='(a)') "Lib_ReadVTK::writeChgcar()"
            write (unit=700,fmt='(f16.8)') dd
            write (unit=700,fmt='(3f12.5)') a(1:3,1)/dd
            write (unit=700,fmt='(3f12.5)') a(1:3,2)/dd
            write (unit=700,fmt='(3f12.5)') a(1:3,3)/dd
            
           
            
            
            write (unit=700,fmt='(a)')   "Du"
            write (unit=700,fmt='(a)')   "1"
            write (unit=700,fmt='(a)')   "direct"
            write (unit=700,fmt='(a)')   "  0.000000  0.000000  0.000000"
            write (unit=700,fmt='(a)')
            write (unit=700,fmt='(3i5)')   nx,ny,nz
            
            
            
           ! dd = (xmax(1)-xmin(1))*(xmax(2)-xmin(2))*(xmax(3)-xmin(3))
            
            datmax = maxval( abs( f(:,:,:) ) )* vol
            nn = ceiling(log10(max(1.0001d0,datmax)))
            
            do iz = 1,nz
                do iy = 1,ny
                    do ix = 1,nx,12
                        dummy = opLine(f(ix:min(ix+11,nx),iy,iz) * vol,nn)
                        write (unit=700,fmt='(a)') trim(dummy)
                    end do
                end do
            end do
           
            !write (unit=700,fmt=*)   f(:,:,:) * dd

            close(unit=700)
            call system( "mv "//trim(filename)//".tmp "//trim(filename) )
            return            
            
            
        end subroutine writeChgcar2_32

         
        pure function determinant3Mat(M) result(d)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      returns the determinant of M
            real(kind=real64),dimension(3,3),intent(in)      ::      M
            real(kind=real64)                                ::      d
            real(kind=real64),dimension(9)       ::      dd
            dd(1) = M(1,1)*( M(2,2)*M(3,3) - M(2,3)*M(3,2) )
            dd(2) = M(1,2)*( M(2,3)*M(3,1) - M(2,1)*M(3,3) )
            dd(3) = M(1,3)*( M(2,1)*M(3,2) - M(2,2)*M(3,1) )
            dd(4) = M(2,1)*( M(3,2)*M(1,3) - M(3,3)*M(1,2) )
            dd(5) = M(2,2)*( M(3,3)*M(1,1) - M(3,1)*M(1,3) )
            dd(6) = M(2,3)*( M(3,1)*M(1,2) - M(3,2)*M(1,1) )
            dd(7) = M(3,1)*( M(1,2)*M(2,3) - M(1,3)*M(2,2) )
            dd(8) = M(3,2)*( M(1,3)*M(2,1) - M(1,1)*M(2,3) )
            dd(9) = M(3,3)*( M(1,1)*M(2,2) - M(1,2)*M(2,1) )
            d = (1.0d0/3.0d0) * sum(dd)
            return
        end function determinant3Mat
         
        pure function opLine(dat,n) result(line)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      generate a single output line eg
    !*      123.456 123.456 123.456 654.321 1.00 1.00 1.00
            real(kind=real64),dimension(:),intent(in)       ::      dat
            integer,intent(in)                              ::      n
            character(len=240)                              ::      line
            integer     ::      kk
            call ftoa( dat,n,8,line,kk )   
            return
        end function opLine
    
         
         
    end module Lib_ReadVTK         
         