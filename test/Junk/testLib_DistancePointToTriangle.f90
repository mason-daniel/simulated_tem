     
    program testLib_DistancePointToTriangle
!---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
!*      simple program to test functioning of Lib_DeformationGradients
!*      
!*          successful result        
!*
!*      distance (0,0,0)-tri     -11.54700538
!*       !! note: generates a file "testDistancePointToTriangle.xyz" to demonstrate correct working
!*
!*                 26222.659
!*
!*       done
!*


        use iso_fortran_env
        use Lib_ColouredTerminal
        use Lib_DistancePointToTriangle
        implicit none
        
        integer,parameter      ::      Nx = 20
        real(kind=real64),dimension(3)      ::      p1,p2,p3
        
        real(kind=real64),dimension(3,0:Nx-1)      ::      xx
        integer             ::       ix,iy,iz
        real(kind=real64)   ::      ddsum
        real(kind=real64),dimension(0:Nx-1)    ::      dd
        
        
        
        character(len=256),dimension(2) ::      output
        character(len=*),dimension(2),parameter   ::      output0 = (/  "distance (0,0,0)-tri     -11.54700538 ",         &
                                                                        "26222.659                             "  /)      
        integer                     ::      ii                                                             
        logical                     ::      ok      
        
        
        
        p1 = (/ Nx/3.0d0,Nx/3.0d0,Nx/3.0d0 /)
        p2 = p1 + (/ Nx/3.0d0,0.0d0,0.0d0 /)
        p3 = p1 + (/ 0.0d0,Nx/3.0d0,0.0d0 /)
        
        
        write(output(1),fmt='(a,f16.8)') "distance (0,0,0)-tri " , distancePointToTriangle( (/0.0d0,0.0d0,0.0d0/),p1,p2,p3 )
         

!        print *,"note: generates a file ""Output/testDistancePointToTriangle.xyz"" to demonstrate correct working"
!        print *,""
!       open(unit=777,file="Output/testDistancePointToTriangle.xyz",action="write")
!       write(unit=777,fmt='(i10)') Nx**3
!       write(unit=777,fmt='(a,9i4,a)') "Lattice=""",Nx,0,0,0,Nx,0,0,0,Nx,""" Properties=species:S:1:pos:R:3:dist:R:1"
!       
        ddsum = 0
        do iz = 0,Nx-1
            do iy = 0,Nx-1
               
                do ix = 0,Nx-1
                    xx(1:3,ix) = (/ix,iy,iz/) 
                end do     
                
                dd = distancePointToTriangle( xx,p1,p2,p3 )
                ddsum = ddsum + sum(dd)
                   
                do ix = 0,Nx-1 
                 
                   write(unit=777,fmt='(a4,4f10.4)') "D ",xx(:,ix),dd(ix)
                end do
                
            end do
        end do
!        close(unit=777)
        write(output(2),fmt='(g24.8)') ddsum
        
        
        
       
    !---    here is the simple test: does the output look like the stored output?
        ok = .true.
        do ii = 1,size(output)            
            if ( trim(cutSpaces(output(ii))) == trim(cutSpaces(output0(ii))) ) then
                write (*,fmt='(a)') trim(cutSpaces(output(ii)))
            else
                ok = .false.
                write (*,fmt='(a)') trim(cutSpaces(output0(ii)))//"    "//colour(RED,trim(cutSpaces(output(ii))))
            end if
        end do   
       
       
    !---    output the result "PASS" or "FAIL"    
        if (ok) then
            print *,colour(LIGHT_GREEN,"PASS")
        else
            print *,colour(RED,"FAIL")
        end if
        
        print *,""
        print *,"done"
        print *,""  
        
    end program testLib_DistancePointToTriangle