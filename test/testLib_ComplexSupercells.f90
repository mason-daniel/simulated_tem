     
    program testLib_ComplexSupercells
!---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^   
!*  
!*      simple program to test functioning of Lib_ComplexSupercells
!*      
!*          successful result        
!*          
!*              ComplexSupercell [nodes =  1049760]
!*                periodic repeat vectors A
!*                  256.91000366    0.00000000    0.00000000
!*                    0.00000000  253.85949707    0.00000000
!*                    0.00000000    0.00000000  255.62690735
!*                primitive unit cell vectors b
!*                   -1.58586422    1.58586422    1.58586422
!*                    1.58662186   -1.58662186    1.58662186
!*                    1.57794387    1.57794387   -1.57794387
!*                unit cell vector repeats A = b n
!*                             0            80            81
!*                            81             0            81
!*                            81            80             0
!*              Supercell volume    16671742.6003
!*              unit cell volume          15.8815
!*              supercell sides          256.9100        253.8595        255.6269
!*              unit cell sides            2.7427          2.7427          2.7427
!*              point in supercell              10.0000         10.0000        -10.0000   F
!*              wrapped point in super          10.0000         10.0000        245.6269   T
!*              unit cell with z axis along (/1,2,3/)
!*                        3.0546          0.0029          0.8479
!*                       -0.4751          2.6229          1.6843
!*                       -0.7122         -1.7482          2.5252
!*              
!*               done
!*              
!*          
        use iso_fortran_env
        use Lib_ColouredTerminal
        use Lib_ComplexSupercells         
        implicit none
        
        
        real(kind=real64),parameter             ::      a0 = 3.1652d0
        real(kind=real64),dimension(3,3)        ::      A = reshape( (/256.91000688,0.0,0.0 , 0.0,253.85950143,0.0 , 0.0,0.0,255.62690141/) ,(/3,3/) )  
        
        real(kind=real64),dimension(3,3)        ::      b = reshape( (/-1.0,1.0,1.0 , 1.0,-1.0,1.0 , 1.0,1.0,-1.0/) ,(/3,3/) ) * a0/2
      
         
        type(ComplexSupercell)          ::      super , super2
        integer                         ::      ii 
        
        real(kind=real64),dimension(3)  ::      xx
        
        character(len=256),dimension(9) ::      output
        character(len=*),dimension(9),parameter   ::      output0 = (/  "Supercell volume    16671742.6003                                          ",         &
                                                                        "unit cell volume          15.8815                                          ",         &
                                                                        "supercell sides          256.9100        253.8595        255.6269          ",         &
                                                                        "unit cell sides            2.7427          2.7427          2.7427          ",         &
                                                                        "point in supercell              10.0000         10.0000        -10.0000   F",         &
                                                                        "wrapped point in super          10.0000         10.0000        245.6269   T",         &
                                                                        "          3.0546          0.0029          0.8479                           ",         &
                                                                        "         -0.4751          2.6229          1.6843                           ",         &
                                                                        "         -0.7122         -1.7482          2.5252                           "          &
                                                                     /)      
                                                                     
        logical                     ::      ok                                                             
                                                                     
                                                                     
                                                                     
                                                                     
                                                                                      
        super = ComplexSupercell_ctor( A,b )
        call report(super) 
        
        write (output(1),fmt='(a,f16.4)')  "Supercell volume ",supercellVolume(super)
        write (output(2),fmt='(a,f16.4)')  "unit cell volume ",unitcellVolume(super)
        write (output(3),fmt='(a,3f16.4)') "supercell sides  ",( getSuperCellSideLength(super,ii),ii=1,3 )
        write (output(4),fmt='(a,3f16.4)') "unit cell sides  ",( getCellSideLength(super,ii),ii=1,3 )

        
        xx = (/10.0d0,10.0d0,-10.0d0/)
        write (output(5),fmt='(a,3f16.4,l4)') "point in supercell     ",xx,pointInSupercell(super,xx)
        xx = wrapPBC(super,xx,realspace=.true.)
        write (output(6),fmt='(a,3f16.4,l4)') "wrapped point in super ",xx,pointInSupercell(super,xx)
                
         
        
    !---    suggest a second supercell with one axis pointing along 1,2,3
        xx = (/ 1.0d0,2.0d0,3.0d0 /)
        call suggestSupercellOrientedWithN( super, xx ,a0, super2 )
         
        b = getb(super2)
        write (*,fmt='(a)') "unit cell with z axis along (/1,2,3/)"
        write (output(7),fmt='(3f16.4)') b(1,:)
        write (output(8),fmt='(3f16.4)') b(2,:)
        write (output(9),fmt='(3f16.4)') b(3,:)                
        
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
        
    end program testLib_ComplexSupercells
        
        