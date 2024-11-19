
!   gfortran -ffree-line-length-256 src/Lib_Lattices.f90 src/testLib_Lattices.f90 -o Test/testLib_Lattices.exe

    
    program testLib_Lattices
!---^^^^^^^^^^^^^^^^^^^^^^^^
!*      simple program to test correct functioning of Lib_Lattices 
!*
!*          correct functioning
!*
!*              $ ./Test/testLib_Lattices.exe
!*              Lattice [fcc   ,cubic       ,nMotif=   1,nNeigh(max)=  19]
!*                      neighbours, motif  1
!*                     0.000   0.000   0.000
!*                     0.000  -0.500  -0.500
!*                     0.500   0.000  -0.500
!*                     0.500  -0.500   0.000
!*                    -0.500   0.000  -0.500
!*                     0.000   0.500  -0.500
!*                    -0.500  -0.500   0.000
!*                     0.500   0.500   0.000
!*                     0.000  -0.500   0.500
!*                     0.500   0.000   0.500
!*                    -0.500   0.500   0.000
!*                    -0.500   0.000   0.500
!*                     0.000   0.500   0.500
!*                     0.000   0.000  -1.000
!*                     0.000  -1.000   0.000
!*                     1.000   0.000   0.000
!*                    -1.000   0.000   0.000
!*                     0.000   1.000   0.000
!*                     0.000   0.000   1.000
!*              getOmega0                  0.25000
!*              getnMotif                  1
!*              getConventionalnMotif      4
!*              diff spot      1     2.00000     0.00000     0.00000
!*              diff spot      2     0.00000     2.00000     0.00000
!*              diff spot      3     0.00000     0.00000     2.00000
!*              
!*               primitive cell
!*                    0.00000000      0.50000000      0.50000000
!*                    0.50000000      0.00000000      0.50000000
!*                    0.50000000      0.50000000      0.00000000
!*               conventional cell
!*                    1.00000000      0.00000000      0.00000000
!*                    0.00000000      1.00000000      0.00000000
!*                    0.00000000      0.00000000      1.00000000
!*              
!*               done
!*              

        use iso_fortran_env
        use Lib_ColouredTerminal
        use Lib_Lattices
        implicit none
        
        type(Lattice)           ::      latt
         
        real(kind=real64),dimension(3,3)        ::      bb
        
        
        


        character(len=256),dimension(14)            ::      output
        character(len=*),dimension(14),parameter   ::      output0 = (/ " getOmega0                 0.69282                   ",    &
                                                                        " getnMotif                  2                        ",    &
                                                                        " getConventionalnMotif      2                        ",    &
                                                                        " diff spot      1      1.00000 1.00000 0.00000       ",    &
                                                                        " diff spot      2     -1.00000 1.00000 0.00000       ",    &
                                                                        " diff spot      3      0.00000 0.00000 1.00000       ",    &                                                                        
                                                                        "  primitive cell                                     ",    &
                                                                        "       0.50000000  0.50000000 0.00000000             ",    &
                                                                        "       -0.86602540 0.86602540 0.00000000             ",    &
                                                                        "       0.00000000 0.00000000 1.60000000              ",    &
                                                                        "  conventional cell                                  ",    &
                                                                        "        0.50000000 0.50000000 0.00000000             ",    &
                                                                        "       -0.86602540 0.86602540 0.00000000             ",    &
                                                                        "        0.00000000 0.00000000 1.60000000             "  /)
                                                                        
                                                                        
                                                                        
                                                                        
        logical                     ::      ok 
        integer                     ::      ii                                                             
                                                                   
                                                      
        
        
        
        
        
        
        
        
        
        
        latt = Lattice_ctor("hcp")
        call setCoverA(latt,1.60d0)
        call report(latt)
        
        write(output(1),fmt='(a,f12.5)') "getOmega0             ",getOmega0(latt)
        write(output(2),fmt='(a,i6)')    "getnMotif             ",getnMotif(latt)
        write(output(3),fmt='(a,i6)')    "getConventionalnMotif ",getConventionalnMotif(latt)
        do ii = 1,getNDiffractionSpots(latt)
            write(output(3+ii),fmt='(a,i6,3f12.5)') "diff spot ",ii,getDiffractionSpot(latt,ii)
        end do
        print *,""
        bb = getPrimitiveCell(latt)
        output(7)="primitive cell"
        write(output(8) ,fmt='(3f16.8)') bb(1,:)
        write(output(9) ,fmt='(3f16.8)') bb(2,:)
        write(output(10),fmt='(3f16.8)') bb(3,:)
        bb = getConventionalCell(latt)
        output(11)="conventional cell"
        write(output(12) ,fmt='(3f16.8)') bb(1,:)
        write(output(13),fmt='(3f16.8)') bb(2,:)
        write(output(14),fmt='(3f16.8)') bb(3,:)
        
        print *,""
        
        call delete(latt)
         
        
        
        
        
        
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
        
        
    !---    list permitted reflections
        latt = Lattice_ctor("fcc")        
        call listPermittedReflections( latt,rho_in=20.0d0 )
        
        print *,""
        print *,"done"
        print *,""
        
    end program testLib_Lattices
    
        
   
        
        
    