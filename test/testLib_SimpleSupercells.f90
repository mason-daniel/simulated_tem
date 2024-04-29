
    program testLib_SimpleSupercells
!---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
!*  
!*      simple program to test functioning of Lib_SimpleSupercells
!*      
!*          successful result        
!*          
!*              Supercell[Nx,Ny,Nz =    6,   7,   8 ]
!*                lattice vectors
!*                    1.00000000    0.00000000    0.00000000
!*                    0.00000000    1.10000000    0.00000000
!*                    0.00000000    0.00000000    1.20000000
!*              cell space                  3.05000     5.55455    -1.02500
!*              real space                  3.05000     6.11000    -1.23000
!*              real space (pbc)            3.05000     6.11000     8.37000
!*              minimum image               3.05000    -1.59000    -1.23000
!*              pbc/min image in cell?      T     F
!*              
!*               align supercell with this vector
!*              Supercell[Nx,Ny,Nz =    5,   6,  11 ]
!*                lattice vectors
!*                    1.10433751   -0.00000000    0.31801890
!*                   -0.19163864    0.83720265    0.63708048
!*                   -0.26252298   -0.61114793    0.87272727
!*              old cell/super vols         1.32000   443.52000
!*              old cell/super vols         1.34400   443.52000
!*              
!*               done
!*          
!*          
        use Lib_ColouredTerminal   
        use Lib_SimpleSupercells
        use iso_fortran_env
        implicit none
        
        type(SimpleSupercell)                   ::      super , super2
        
        real(kind=real64),dimension(3,3)        ::      a_cell
        integer                                 ::      Nx,Ny,Nz
        real(kind=real64),dimension(3)          ::      uu
      

        character(len=256),dimension(10)           ::      output
        character(len=*),dimension(10),parameter   ::      output0 = (/ "cell space                  3.05000     5.55455    -1.02500 ", &
                                                                        "real space                  3.05000     6.11000    -1.23000 ", &
                                                                        "real space (pbc)            3.05000     6.11000     8.37000 ", &
                                                                        "minimum image               3.05000    -1.59000    -1.23000 ", &
                                                                        "pbc/min image in cell?      T     F                         ", &
                                                                        " 1.104338 0.000000 0.318019                                 ", &   
                                                                        " -0.191639 0.837203 0.637080                                ", &
                                                                        "-0.262523 -0.611148 0.872727                                ", &
                                                                        "old cell/super vols         1.32000   443.52000             ",  &
                                                                        "old cell/super vols         1.34400   443.52000             "   /)
                                                                        
        character(len=*),dimension(10),parameter   ::      output0_alternate = (/ "cell space                  3.05000     5.55455    -1.02500 ", &     !   negative zero output is acceptable
                                                                        "real space                  3.05000     6.11000    -1.23000 ", &
                                                                        "real space (pbc)            3.05000     6.11000     8.37000 ", &
                                                                        "minimum image               3.05000    -1.59000    -1.23000 ", &
                                                                        "pbc/min image in cell?      T     F                         ", &
                                                                        " 1.104338 -0.000000 0.318019                                ", &   
                                                                        " -0.191639 0.837203 0.637080                                ", &
                                                                        "-0.262523 -0.611148 0.872727                                ", &
                                                                        "old cell/super vols         1.32000   443.52000             ",  &
                                                                        "old cell/super vols         1.34400   443.52000             "   /)                                                                        
                                                                         
        logical                     ::      ok 
        integer                     ::      ii          
        
        
        
        
        
        
        
        
        
    !---    create a supercell
        a_cell = reshape( (/1.0d0,0.0d0,0.0d0 , 0.0d0,1.1d0,0.0d0 , 0.0d0,0.0d0,1.2d0 /) , (/3,3/) )        
        Nx = 6; Ny = 7; Nz = 8     
        super = SimpleSupercell_ctor(a_cell,Nx,Ny,Nz)        
        call report(super)
        
    !---    choose a random point within the cell. Convert real to cell space and back  
        uu = (/ 3.05d0,6.11d0,-1.23d0 /)
        uu = realSpaceToCell( super,uu )
        write (output(1),fmt='(a,3f12.5)') "cell space             ",uu
        uu = cellToRealSpace( super,uu )
        
        write (output(2),fmt='(a,3f12.5)') "real space             ",uu
        write (output(3),fmt='(a,3f12.5)') "real space (pbc)       ",wrapPBC(super,uu)
        write (output(4),fmt='(a,3f12.5)') "minimum image          ",minimumImage(super,uu)
        write (output(5),fmt='(a,2l6)')    "pbc/min image in cell? ",pointInSupercell(super,wrapPBC(super,uu)),pointInSupercell(super,minimumImage(super,uu)) 
        
        
    !---    try for a supercell aligned with this vector
        uu = wrapPBC(super,uu)
        print *,""
        print *,"align supercell with this vector"
        call suggestSupercellOrientedWithN( super, uu , super2 )
        call report(super2)
        
        a_cell = getA(super2)
        ! print *,"a_cell(1,2)",a_cell(1,2), (a_cell(1,2)==0), (a_cell(1,2)==-0)
        ! where (a_cell==-0.0d0)
        !     a_cell=0.0d0
        ! endwhere
        write (output(6),fmt='(3f16.6)')  a_cell(1,1:3)
        write (output(7),fmt='(3f16.6)')  a_cell(2,1:3)
        write (output(8),fmt='(3f16.6)')  a_cell(3,1:3)
        write (output(9),fmt='(a,3f12.5)') "old cell/super vols    ",unitCellVolume(super),superCellVolume(super)
        write (output(10),fmt='(a,3f12.5)') "old cell/super vols    ",unitCellVolume(super2),superCellVolume(super2)
        
        
        
        call delete(super)
        
        
        
    !---    here is the simple test: does the output look like the stored output?
        ok = .true.
        do ii = 1,size(output)            
            if ( ( trim(cutSpaces(output(ii))) == trim(cutSpaces(output0(ii))) ).or.( trim(cutSpaces(output(ii))) == trim(cutSpaces(output0_alternate(ii))) ) ) then
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
        
    end program testLib_SimpleSupercells
        