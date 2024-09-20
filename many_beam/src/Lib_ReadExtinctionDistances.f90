
    module Lib_ReadExtinctionDistances
!---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
!*      a short module to find an extinction distance in a table 
!*      The extinction distances are expected to be in the format 
!*      # crystalStructureFactors.exe v.1.0
!*       element          g        T (K)      V (V)         Re(U_g)         Im(U_g)            xi_g
!*            FE   0   0   0     300.000  150000.000     29.10544372      0.41419664    135.09100206
!*            FE   1   1   0     300.000  150000.000     11.85052544      0.37953687    331.65335318
!*      ...


        use iso_fortran_env
        implicit none
        private
               
        
        public      ::      readExtinctionDistance
       
        
        integer,public,parameter       ::      LATTICE_SYM_CUSTOM      = 0
        integer,public,parameter       ::      LATTICE_SYM_CUBIC       = 1
        integer,public,parameter       ::      LATTICE_SYM_HEXAGONAL   = 2
        
        interface           readExtinctionDistance
            module procedure    readExtinctionDistance0
            module procedure    readExtinctionDistance1
        end interface
         
    contains
!---^^^^^^^^

        subroutine readExtinctionDistance0( filename,hkl,symm,xi_g, ok )
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      attempt to read in xi_g for reflection [hkl]
            character(len=*),intent(in)             ::      filename
            integer,dimension(:,:),intent(in)       ::      hkl         !   3,nReflection
            integer,intent(in)                      ::      symm
            real(kind=real64),dimension(:),intent(out)           ::      xi_g
            logical,intent(out)                     ::      ok
            
            integer                         ::      ii,ioerr,nReflection,jj
            character(len=256)              ::      dummy
            integer,dimension(3,size(hkl,dim=2))    ::      hkl_s
            logical,dimension(size(hkl,dim=2))      ::      found
            integer,dimension(3)            ::      idat 
            real(kind=real64),dimension(5)  ::      rdat
            
            
            ok = .false.
            xi_g = huge(1.0)        !   most likely reason for failure is [hkl] too large, and xi_g -> infinity.
            
            inquire(file=trim(filename),exist=ok)
            if (.not. ok) then
                print *,"Lib_ReadExtinctionDistances::readExtinctionDistance0 error - file not found """//trim(filename)//""""
                return
            end if
            
            
        !---    not all hkl reflections will be listed.
            nReflection = size(hkl,dim=2)
            do jj = 1,nReflection
                hkl_s(:,jj) = symmetry( hkl(:,jj),symm )     !   symmetry-related reflection to search
            end do        
            
            
            found = .false.
            open(unit=850,file=trim(filename),action="read")
            
            !---    first two lines of file are expected to be comments
                read(unit=850,fmt='(a)',iostat=ioerr) dummy
                read(unit=850,fmt='(a)',iostat=ioerr) dummy
                if (ioerr/=0) then
                    print *,"Lib_ReadExtinctionDistances::readExtinctionDistance0 error - can't read header """//trim(filename)//""""
                    return
                end if
                
            !---    find the line we are looking for.
                do ii = 1,1000      !   don't search forever... if the file is too long, its bound to be an error!
                    read(unit=850,fmt=*,iostat=ioerr) dummy,idat,rdat

                    if (ioerr/=0) then
                        print *,"Lib_ReadExtinctionDistances::readExtinctionDistance0 error - end of file reached/incompatible format """//trim(filename)//""""
                        exit
                    end if
                    
                    idat = symmetry(idat,symm)           !   make standard form
                    do jj = 1,nReflection
                        if (all(idat(1:3) == hkl_s(1:3,jj))) then
                            found(jj) = .true.
                            xi_g(jj) = rdat(5)
                            
                        end if
                    end do
                    if (all(found)) exit
                         
                end do
                
                
            close(unit=850)
            
            ok = all(found)
            
            if (.not. ok) then
                do jj = 1,nReflection
                    if (.not. found(jj)) write(*,fmt='(a,3i4,a)') " Lib_ReadExtinctionDistances::readExtinctionDistance0 warning - reflection ",hkl(:,jj)," not found"
                end do
            end if
            
            
            
            return            
        end subroutine readExtinctionDistance0          
        

        subroutine readExtinctionDistance1( filename,hkl,symm,xi_g, ok )
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      attempt to read in xi_g for reflection [hkl]
            character(len=*),intent(in)             ::      filename
            integer,dimension(:,:),intent(in)       ::      hkl         !   3,nReflection
            integer,intent(in)                      ::      symm
            complex(kind=real64),dimension(:),intent(out)           ::      xi_g
            logical,intent(out)                     ::      ok
            
            integer                         ::      ii,ioerr,nReflection,jj
            character(len=256)              ::      dummy
            integer,dimension(3,size(hkl,dim=2))    ::      hkl_s
            logical,dimension(size(hkl,dim=2))      ::      found
            integer,dimension(3)            ::      idat 
            real(kind=real64),dimension(5)  ::      rdat
            
            
            ok = .false.
            xi_g = huge(1.0)        !   most likely reason for failure is [hkl] too large, and xi_g -> infinity.
            
            inquire(file=trim(filename),exist=ok)
            if (.not. ok) then
                print *,"Lib_ReadExtinctionDistances::readExtinctionDistance1 error - file not found """//trim(filename)//""""
                return
            end if
            
            
        !---    not all hkl reflections will be listed.
            nReflection = size(hkl,dim=2)
            do jj = 1,nReflection
                hkl_s(:,jj) = symmetry( hkl(:,jj),symm )     !   symmetry-related reflection to search
            end do        
            
            
            found = .false.
            open(unit=850,file=trim(filename),action="read")
            
            !---    first two lines of file are expected to be comments
                read(unit=850,fmt='(a)',iostat=ioerr) dummy
                read(unit=850,fmt='(a)',iostat=ioerr) dummy
                if (ioerr/=0) then
                    print *,"Lib_ReadExtinctionDistances::readExtinctionDistance1 error - can't read header """//trim(filename)//""""
                    return
                end if
                
            !---    find the line we are looking for.
                do ii = 1,1000      !   don't search forever... if the file is too long, its bound to be an error!
                    read(unit=850,fmt=*,iostat=ioerr) dummy,idat,rdat
                     if (ioerr/=0) then
                        print *,"Lib_ReadExtinctionDistances::readExtinctionDistance1 error - end of file reached/incompatible format """//trim(filename)//""""
                        exit
                    end if
                    
                    idat = symmetry(idat,symm)           !   make standard form
                    do jj = 1,nReflection
                        if (all(idat(1:3) == hkl_s(1:3,jj))) then
                            found(jj) = .true.
                            xi_g(jj) = complex( rdat(3),rdat(4) )                            
                        end if
                    end do
                    if (all(found)) exit
                         
                end do
                
                
            close(unit=850)
            
            ok = all(found)
            
            if (.not. ok) then
                do jj = 1,nReflection
                    if (.not. found(jj)) write(*,fmt='(a,3i4,a)') " Lib_ReadExtinctionDistances::readExtinctionDistance1 warning - reflection ",hkl(:,jj)," not found"
                end do
            end if
            
            
            
            return            
        end subroutine readExtinctionDistance1             
            
        pure function symmetry( hkl_in,symm ) result( hkl_out )
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      find a standardised symmetry-related reflection
            integer,dimension(3),intent(in)     ::      hkl_in
            integer,intent(in)                  ::      symm
            integer,dimension(3)                ::      hkl_out
            
            
            integer         ::      hh,kk,ll!,qq
            
            select case( symm )
            
                case(LATTICE_SYM_CUBIC)
                    !   h,k,l   -> |h| |k| |l|  with |h| >= |k| >= |l|
                    
                    hh = abs(hkl_in(1))
                    kk = abs(hkl_in(2))
                    ll = abs(hkl_in(3))
                
                    if (hh>=max(kk,ll)) then
                        if (kk>=ll) then
                            hkl_out(1:3) = (/ hh,kk,ll /)
                        else
                            hkl_out(1:3) = (/ hh,ll,kk /)
                        end if
                    
                    else if (kk>=max(hh,ll)) then
                        if (hh>=ll) then
                            hkl_out(1:3) = (/ kk,hh,ll /)
                        else
                            hkl_out(1:3) = (/ kk,ll,hh /)
                        end if
                    
                    else  
                        if (hh>=kk) then
                            hkl_out(1:3) = (/ ll,hh,kk /)
                        else
                            hkl_out(1:3) = (/ ll,kk,hh /)
                        end if
                    end if

                case(LATTICE_SYM_HEXAGONAL)
                    !   h k l -> |h| k |l| with |h| >= |k|
                 
                    hh = abs(hkl_in(1))
                    kk = abs(hkl_in(2))
                    ll = abs(hkl_in(3))
                
                    
                    if ( hkl_in(1)*hkl_in(2) >= 0 ) then
                        !   h k l or -h -k l or -h 0 l or h 0 l 
                        if (hh>=kk) then
                            hkl_out(1:3) = (/ hh,kk,ll /)
                        else
                            hkl_out(1:3) = (/kk,hh,ll /)
                        end if
                    else
                        !   h -k l or -h k l
                        if (hh>=kk) then
                            hkl_out(1:3) = (/ hh,-kk,ll /)
                        else
                            hkl_out(1:3) = (/kk,-hh,ll /)
                        end if                    
                    end if
               
            end select                    
                                
            return
        end function symmetry
            
    end module Lib_ReadExtinctionDistances

    
!!   gfortran -ffree-line-length-256 Lib_ReadExtinctionDistances.f90 -o testLib_ReadExtinctionDistances.exe    
!    
!    program testLib_ReadExtinctionDistances
!!---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
!
!        use Lib_ReadExtinctionDistances
!        use iso_fortran_env
!        implicit none
!        
!        
!        character(len=256)      ::      filename,dummy
!        integer,dimension(3,1)  ::      hkl
!        real(kind=real64),dimension(1)       ::      xi_g
!        logical                 ::      ok
!        
!        print *,"usage: testLib_ReadExtinctionDistances.exe filename h,k,l"
!        call get_command_argument(1,filename)
!        call get_command_argument(2,dummy)
!        read(dummy,fmt=*)   hkl
!        
!        
!        print *,"filename """//trim(filename)//""""
!        print *,"reflection ",hkl
!        
!        
!        call readExtinctionDistance( filename,hkl,xi_g, ok )
!        
!        print *,"extinction distance xi_g = ",xi_g," success ",ok
!        
!        print *,""
!        print *,"done"
!        print *,""
!      
!    end program testLib_ReadExtinctionDistances
        
        
        
        