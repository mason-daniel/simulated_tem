
    program test_Lib_AtomSpace
!---^^^^^^^^^^^^^^^^^^^^^^^^^^
!*      short test program to confirm correct functioning of Lib_AtomSpace

        use Lib_AtomSpace
        
        use Lib_RotationMatrices
        use iso_fortran_env
        implicit none


        integer,parameter                               ::      NATOMS = 10000
        real(kind=real64),dimension(3,3)                ::      a_super
        real(kind=real64),dimension(3)                  ::      xyzoffset
        real(kind=real64),dimension(3,NATOMS)           ::      x           !   atom positions
        real(kind=real64)                               ::      a           !   grid spacing
        type(AtomSpace)         ::      as
        real(kind=real64),dimension(3,3)        ::      R               !   rotation matrix
        real(kind=real64),dimension(3,3)        ::      eps             !   strain matrix
        
       ! integer                 ::      ii
        integer                 ::      Nx,Ny,Nz
        
        logical                 ::      ok


        ok = .true.

        !a_super = 10.d0 *RotationMatrix_identity                                !   start with a simple cubic cell...
        a_super = 0
        a_super(1,1) = 12.0d0
        a_super(2,2) = 10.0d0
        a_super(3,3) =  9.0d0
        xyzoffset = 0                                                           !   ... with zero offset ...
        a = 1.0d0                                                               !   unit grid spacing
        call random_number(x)                                                   !   ... randomly placed atoms ...
        x = rotateVector(a_super,x)

        print *,""
        print *,"simple test case: R = I  "
        R = RotationMatrix_identity
        as = AtomSpace_ctor(a,xyzoffset,a_super)
        call setR(as,R)
        call setThickness(as)
        call suggestImagingSpace(as,Nx,Ny,Nz) 
        print *,"suggestImagingSpace ",Nx,Ny,Nz
        call setDelta(as,Nx,Ny,Nz)
        call report(as) 
        call test(as,ok)
        print *,""
        print *,""


    !   non-trivial rotation matrix - rotation of 30 deg about y-axis
        print *,"test case 1:   R = rotation 30 deg about y-axis"
        R = RotationMatrix_ctor( (/0.0d0,1.0d0,0.0d0/), 30*3.14159265930d0/180 )        
         
        as = AtomSpace_ctor(a,xyzoffset,a_super)
        call setR(as,R)
        call setThickness(as)
        call suggestImagingSpace(as,Nx,Ny,Nz) 
        print *,"suggestImagingSpace ",Nx,Ny,Nz
        call setDelta(as,Nx,Ny,Nz)
        call report(as) 
        call test(as,ok)
        print *,""
        print *,""

    !   strain the supercell 
        print *,"test case 2: A 5% shear strain, R = rotation 30 deg about y-axis "
        eps = 0
        eps(1,2) = 0.05d0 ; eps(2,1) = eps(1,2)
        eps(1,3) = 0.05d0 ; eps(3,1) = eps(1,3)
        eps(2,3) = 0.05d0 ; eps(3,2) = eps(2,3)
        a_super = a_super + matmul(eps,a_super)
        call random_number(x)                                                   !   ... randomly placed atoms ...        
        x = rotateVector(a_super,x)
        as = AtomSpace_ctor(a,xyzoffset,a_super)
        call setR(as,R)
        call setThickness(as,x)
        call suggestImagingSpace(as,Nx,Ny,Nz) 
        print *,"suggestImagingSpace ",Nx,Ny,Nz
        call setDelta(as,Nx,Ny,Nz)
        call report(as) 
        call test(as,ok)
        print *,""
        print *,""



    !   strain the supercell, introduce an explicit surface
        print *,"test case 3: A 10% shear strain,  R = rotation 30 deg about y-axis, surface at +10%,+90% z."       
        a_super = a_super + matmul(eps,a_super)
        call random_number(x)                                                   !   ... randomly placed atoms ...
        x(3,:) = x(3,:)*0.8d0 + 0.1d0                                           !   ... explicit surface
        x = rotateVector(a_super,x)
        as = AtomSpace_ctor(a,xyzoffset,a_super)
        call setR(as,R)
        call setThickness(as,x)
        call suggestImagingSpace(as,Nx,Ny,Nz) 
        print *,"suggestImagingSpace ",Nx,Ny,Nz
        call setDelta(as,Nx,Ny,Nz)
        call report(as) 
        call test(as,ok)
        print *,""
        print *,""

        print *,""
        if (ok) print *,"PASS"
        print *,"done"
        print *,""


    contains
!---^^^^^^^^
    
        subroutine test(as,ok)
    !---^^^^^^^^^^^^^^^^^^^^^^
            type(AtomSpace),intent(in)          ::      as
            logical,intent(inout)               ::      ok
            
            real(kind=real64),dimension(3)          ::      xt,delta , minx,maxx,avgx
            real(kind=real64),dimension(3,9)        ::      xtp
            integer                                 ::      np
            !integer                 ::      Nx,Ny,Nz
            integer                 ::      ii,nn,jj
            real(kind=real64),parameter     ::      TOL = 1.0d-6

            
            delta(:) = ( a_super(:,1) + a_super(:,2) + a_super(:,3) )/2 + xyzoffset(:)
            write(*,fmt='(3(a,3f12.6))') "atom at ",delta," transforms to ",toImagingSpace(as,delta)," cf ",(/Nx,Ny,Nz/)*0.5d0
            ok = ok .and. maxval( abs( toImagingSpace(as,delta) - (/Nx,Ny,Nz/)*0.5d0 ) )<1.0d-3                
            minx = huge(1.0)
            maxx = -huge(1.0)
            avgx = 0
            do ii = 1,NATOMS
                xt(:) = toImagingSpace( as, x(:,ii) )
                minx = min( xt(:),minx )
                maxx = max( xt(:),maxx )
                avgx = avgx + xt(:)
            end do
            print *,"bounds without periodic copies"
            avgx = avgx / NATOMS
            do ii = 1,3
                print *,"min/max/avg ",ii,minx(ii),maxx(ii),avgx(ii)
            end do

            
            nn = 0
            do ii = 1,NATOMS
                call periodicCopies(as,real(x(:,ii),kind=real32),Nx,Ny,0,np,xtp)
                do jj = 1,np
                    nn = nn + 1
                    minx = min( xtp(:,jj),minx )
                    maxx = max( xtp(:,jj),maxx )
                    avgx = avgx + xtp(:,jj)
                end do
            end do
            print *,"bounds with periodic copies, no buffer"
            avgx = avgx / nn
            do ii = 1,3
                print *,"min/max/avg ",ii,minx(ii),maxx(ii),avgx(ii)
            end do
            ok = ok .and. ( (minx(1)>=-TOL).and.(minx(1)<=Nx+TOL) )
            ok = ok .and. ( (minx(2)>=-TOL).and.(minx(2)<=Ny+TOL) )
            ok = ok .and. ( (minx(3)>=-TOL).and.(minx(3)<=Nz+TOL) )
            
            return
        end subroutine test

    end program test_Lib_AtomSpace


