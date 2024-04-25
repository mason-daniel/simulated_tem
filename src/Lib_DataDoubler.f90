
    module Lib_DataDoubler
!---^^^^^^^^^^^^^^^^^^^^^^
!*      module to reduce data set in size by 1/2 and then reexpand it
        use Lib_Splines1   
        use iso_fortran_env
        implicit none 
        private
        
        public      ::      halfLine,doubleLine
        public      ::      halfLine32,doubleLine32
        public      ::      halfImage,doubleImage
        public      ::      halfImage32,doubleImage32
        public      ::      halfField,doubleField
        public      ::      halfField32,doubleField32
        
        public      ::      halfData,doubleData
        
        interface       halfData
            module procedure    halfLine   
            module procedure    halfImage   
            module procedure    halfField   
            module procedure    halfLine32    
            module procedure    halfImage32   
            module procedure    halfField32   
        end interface  
        
        interface       doubleData
            module procedure    doubleLine   
            module procedure    doubleImage   
            module procedure    doubleField   
            module procedure    doubleLine32    
            module procedure    doubleImage32   
            module procedure    doubleField32   
        end interface  
        

    contains
!---^^^^^^^^

        subroutine halfLine( n,y2,y1, pbc )
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      take the data from y2(1:n) and interpolate onto the line y1(1:(n+1)/2)
            integer,intent(in)                              ::      n
            real(kind=real64),dimension(:),intent(in)       ::      y2
            real(kind=real64),dimension(:),intent(out)      ::      y1
            logical,intent(in)                              ::      pbc
    
            integer             ::      ii,mm
            type(Spline1)       ::      ss
            real(kind=real64)   ::      xx
            
            ss = Spline1_ctor(y2,pbc)
             
            mm = (n+1)/2
            y1 = 0.0d0
            if (pbc) then
                do ii = 1,mm
                    xx = real((ii-1)*n,kind=real64)/mm
                    !print *,"half point ",ii," at double point ",xx+1
                    y1(ii) = Splint(ss,xx)                
                end do
            else
                do ii = 1,mm
                    xx = real((ii-0.5d0)*(n-1),kind=real64)/mm
                    !print *,"half point ",ii," at double point ",xx+1
                    y1(ii) = Splint(ss,xx)                
                end do
            end if                
            return
        end subroutine halfLine 
            
        subroutine doubleLine( n,y1,y2, pbc )
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      take the data from y1(1:(n+1)/2) and interpolate onto the line y2(1:n)
            integer,intent(in)                              ::      n
            real(kind=real64),dimension(:),intent(in)       ::      y1
            real(kind=real64),dimension(:),intent(out)      ::      y2
            logical,intent(in)                              ::      pbc
    
            integer             ::      ii,mm
            type(Spline1)       ::      ss
            real(kind=real64)   ::      xx
            
            
            mm = (n+1)/2
            ss = Spline1_ctor(y1(1:mm),pbc)
            !call report(ss,6,.true.)
            if (pbc) then
                do ii = 1,n
                    xx = real((ii-1)*mm,kind=real64)/n
                    !print *,"double point ",ii," at half point ",xx+1
                    y2(ii) = Splint(ss,xx)                
                end do
            else
                do ii = 1,n
                    xx = real((ii-1)*mm,kind=real64)/(n-1)-0.5d0
                    !print *,"double point ",ii," at half point ",xx+1
                    y2(ii) = Splint(ss,xx)                
                end do
            end if                
            return
        end subroutine doubleLine 

        subroutine halfLine32( n,y2,y1, pbc )
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      take the data from y2(1:n) and interpolate onto the line y1(1:(n+1)/2)
            integer,intent(in)                              ::      n
            real(kind=real32),dimension(:),intent(in)       ::      y2
            real(kind=real32),dimension(:),intent(out)      ::      y1
            logical,intent(in)                              ::      pbc
    
            integer             ::      ii,mm
            type(Spline1)       ::      ss
            real(kind=real64)   ::      xx
            real(kind=real64),dimension(size(y2))           ::      y264
            
            y264 = y2
            ss = Spline1_ctor(y264,pbc)
             
            mm = (n+1)/2
            y1 = 0.0d0
            if (pbc) then
                do ii = 1,mm
                    xx = real((ii-1)*n,kind=real64)/mm
                    !print *,"half point ",ii," at double point ",xx+1
                    y1(ii) = real( Splint(ss,xx),kind=real32 )                
                end do
            else
                do ii = 1,mm
                    xx = real((ii-0.5d0)*(n-1),kind=real64)/mm
                    !print *,"half point ",ii," at double point ",xx+1
                    y1(ii) = real( Splint(ss,xx),kind=real32 )                
                end do
            end if                
            return
        end subroutine halfLine32 
            
        subroutine doubleLine32( n,y1,y2, pbc )
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      take the data from y1(1:(n+1)/2) and interpolate onto the line y2(1:n)
            integer,intent(in)                              ::      n
            real(kind=real32),dimension(:),intent(in)       ::      y1
            real(kind=real32),dimension(:),intent(out)      ::      y2
            logical,intent(in)                              ::      pbc
    
            integer             ::      ii,mm
            type(Spline1)       ::      ss
            real(kind=real64)   ::      xx
            real(kind=real64),dimension(size(y1))       ::      y164
            
            y164 = y1
            mm = (n+1)/2
            ss = Spline1_ctor(y164(1:mm),pbc)
            !call report(ss,6,.true.)
            if (pbc) then
                do ii = 1,n
                    xx = real((ii-1)*mm,kind=real64)/n
                    !print *,"double point ",ii," at half point ",xx+1
                    y2(ii) = real( Splint(ss,xx),kind=real32 )
                end do
            else
                do ii = 1,n
                    xx = real((ii-1)*mm,kind=real64)/(n-1)-0.5d0
                    !print *,"double point ",ii," at half point ",xx+1
                    y2(ii) = real( Splint(ss,xx),kind=real32 )                
                end do
            end if                
            return
        end subroutine doubleLine32 
            
        subroutine halfImage( nx,ny,y2,y1, pbc )
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      take the data from y2(1:n) and interpolate onto the line y1(1:(n+1)/2)
            integer,intent(in)                              ::      nx,ny
            real(kind=real64),dimension(:,:),intent(in)     ::      y2
            real(kind=real64),dimension(:,:),intent(out)    ::      y1
            logical,intent(in)                              ::      pbc
            
            real(kind=real64),dimension(:,:),allocatable    ::      y2_tmp
    
            integer             ::      ii,mx
            mx = (nx+1)/2
            allocate(y2_tmp(mx,ny))
            do ii = 1,ny
                call halfLine( nx,y2(:,ii),y2_tmp(:,ii),pbc )
            end do
            do ii = 1,mx
                call halfLine( ny,y2_tmp(ii,:),y1(ii,:),pbc )
            end do
            deallocate(y2_tmp)
            
             
            return
        end subroutine halfImage
            
        subroutine doubleImage( nx,ny,y1,y2, pbc )
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      take the data from y1(1:(n+1)/2) and interpolate onto the line y2(1:n)
            integer,intent(in)                              ::      nx,ny
            real(kind=real64),dimension(:,:),intent(in)     ::      y1
            real(kind=real64),dimension(:,:),intent(out)    ::      y2
            logical,intent(in)                              ::      pbc
    
            real(kind=real64),dimension(:,:),allocatable    ::      y1_tmp
    
            integer             ::      ii,my
            my = (ny+1)/2
            allocate(y1_tmp(nx,my))
            do ii = 1,my
                call doubleLine( nx,y1(:,ii),y1_tmp(:,ii),pbc )
            end do
            do ii = 1,nx
                call doubleLine( ny,y1_tmp(ii,:),y2(ii,:),pbc )
            end do
            deallocate(y1_tmp)
            
            return
        end subroutine doubleImage
           
        subroutine halfImage32( nx,ny,y2,y1, pbc )
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      take the data from y2(1:n) and interpolate onto the line y1(1:(n+1)/2)
            integer,intent(in)                              ::      nx,ny
            real(kind=real32),dimension(:,:),intent(in)     ::      y2
            real(kind=real32),dimension(:,:),intent(out)    ::      y1
            logical,intent(in)                              ::      pbc
            
            real(kind=real32),dimension(:,:),allocatable    ::      y2_tmp
    
            integer             ::      ii,mx
            mx = (nx+1)/2
            allocate(y2_tmp(mx,ny))
            do ii = 1,ny
                call halfLine32( nx,y2(:,ii),y2_tmp(:,ii),pbc )
            end do
            do ii = 1,mx
                call halfLine32( ny,y2_tmp(ii,:),y1(ii,:),pbc )
            end do
            deallocate(y2_tmp)
            
             
            return
        end subroutine halfImage32
            
        subroutine doubleImage32( nx,ny,y1,y2, pbc )
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      take the data from y1(1:(n+1)/2) and interpolate onto the line y2(1:n)
            integer,intent(in)                              ::      nx,ny
            real(kind=real32),dimension(:,:),intent(in)     ::      y1
            real(kind=real32),dimension(:,:),intent(out)    ::      y2
            logical,intent(in)                              ::      pbc
    
            real(kind=real32),dimension(:,:),allocatable    ::      y1_tmp
    
            integer             ::      ii,my
            my = (ny+1)/2
            allocate(y1_tmp(nx,my))
            do ii = 1,my
                call doubleLine32( nx,y1(:,ii),y1_tmp(:,ii),pbc )
            end do
            do ii = 1,nx
                call doubleLine32( ny,y1_tmp(ii,:),y2(ii,:),pbc )
            end do
            deallocate(y1_tmp)
            
            return
        end subroutine doubleImage32
        
         
        
        subroutine halfField( nx,ny,nz,y2,y1, pbc )
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      take the data from y2(1:n) and interpolate onto the line y1(1:(n+1)/2)
            integer,intent(in)                                  ::      nx,ny,nz
            real(kind=real64),dimension(:,:,:),intent(in)       ::      y2
            real(kind=real64),dimension(:,:,:),intent(out)      ::      y1
            logical,intent(in)                                  ::      pbc
            
            real(kind=real64),dimension(:,:,:),allocatable      ::      y2_tmp
    
            integer             ::      ii,jj,mx,my
                        
            mx = (nx+1)/2
            my = (ny+1)/2
            allocate(y2_tmp(mx,my,nz))
            do ii = 1,nz
                call halfImage( nx,ny,y2(:,:,ii),y2_tmp(:,:,ii),pbc )
            end do
            do jj = 1,my
                do ii = 1,mx
                    call halfLine( nz,y2_tmp(ii,jj,:),y1(ii,jj,:),pbc )                   
                end do
            end do
            deallocate(y2_tmp)
            
             
            return
        end subroutine halfField
            
        subroutine doubleField( nx,ny,nz,y1,y2, pbc )
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      take the data from y1(1:(n+1)/2) and interpolate onto the line y2(1:n)
            integer,intent(in)                              ::      nx,ny,nz
            real(kind=real64),dimension(:,:,:),intent(in)   ::      y1
            real(kind=real64),dimension(:,:,:),intent(out)  ::      y2
            logical,intent(in)                              ::      pbc
    
            real(kind=real64),dimension(:,:,:),allocatable    ::      y1_tmp
    
            integer             ::      ii,jj,mz 
             
            mz = (nz+1)/2
            allocate(y1_tmp(nx,ny,mz))
            do ii = 1,mz
                call doubleImage( nx,ny,y1(:,:,ii),y1_tmp(:,:,ii),pbc )
            end do           
            do jj = 1,ny                              
                do ii = 1,nx
                    call doubleLine( nz,y1_tmp(ii,jj,:),y2(ii,jj,:),pbc )
                end do
            end do
            deallocate(y1_tmp)
            
            return
        end subroutine doubleField
        
        
        subroutine halfField32( nx,ny,nz,y2,y1, pbc )
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      take the data from y2(1:n) and interpolate onto the line y1(1:(n+1)/2)
            integer,intent(in)                                  ::      nx,ny,nz
            real(kind=real32),dimension(:,:,:),intent(in)       ::      y2
            real(kind=real32),dimension(:,:,:),intent(out)      ::      y1
            logical,intent(in)                                  ::      pbc
            
            real(kind=real32),dimension(:,:,:),allocatable      ::      y2_tmp
    
            integer             ::      ii,jj,mx,my
                        
            mx = (nx+1)/2
            my = (ny+1)/2
            allocate(y2_tmp(mx,my,nz))
            do ii = 1,nz
                call halfImage32( nx,ny,y2(:,:,ii),y2_tmp(:,:,ii),pbc )
            end do
            do jj = 1,my
                do ii = 1,mx
                    call halfLine32( nz,y2_tmp(ii,jj,:),y1(ii,jj,:),pbc )                   
                end do            
            end do
            deallocate(y2_tmp)
            
             
            return
        end subroutine halfField32
            
        subroutine doubleField32( nx,ny,nz,y1,y2, pbc )
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      take the data from y1(1:(n+1)/2) and interpolate onto the line y2(1:n)
            integer,intent(in)                              ::      nx,ny,nz
            real(kind=real32),dimension(:,:,:),intent(in)   ::      y1
            real(kind=real32),dimension(:,:,:),intent(out)  ::      y2
            logical,intent(in)                              ::      pbc
    
            real(kind=real32),dimension(:,:,:),allocatable    ::      y1_tmp
    
            integer             ::      ii,jj,mz 
             
            mz = (nz+1)/2
            allocate(y1_tmp(nx,ny,mz))
            do ii = 1,mz
                call doubleImage32( nx,ny,y1(:,:,ii),y1_tmp(:,:,ii),pbc )
            end do           
            do jj = 1,ny                              
                do ii = 1,nx
                    call doubleLine32( nz,y1_tmp(ii,jj,:),y2(ii,jj,:),pbc )
                end do
            end do
            deallocate(y1_tmp)
            
            return
        end subroutine doubleField32
        
        
    end module Lib_DataDoubler  