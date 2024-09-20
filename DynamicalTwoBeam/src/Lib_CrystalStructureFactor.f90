
    module Lib_CrystalStructureFactor
!---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ 
!*      very short helper module describing the crystal structure factor
!*      in terms of Doyle-Turner ( expansion of gaussian terms ) coefficients
!*
!*      L.M. PENG,  G. REN, S. L. DUDAREV AND M. J. WHELAN
!*      Debye-Waller Factors and Absorptive Scattering Factors of Elemental Crystals
!*      Acta Cryst. (1996). A52, 456-470
!*
!*      call with
!*           Vg = crystalStructureFactor( a0, lattice, g, a,b,dwf )
!*      with
!*          a0        lattice parameter  
!*          lattice = "bcc" or "fcc" 
!*          g         g vector (reduced units eg [110])   
!*          a(:,:),b(:,:) being complex Doyle Turner coefficients for each motif point
!*          dwf       Debye Waller factor (A^-2)
  

!*      units eV/fs/A
!*
        use Lib_Elements
        use Lib_Lattices
        use iso_fortran_env
        implicit none
        private
        
        
        real(kind=real64),parameter,private     ::      ME =  0.056856301036d0      !   electron mass in weird units eV A^-2 fs^2
        real(kind=real64),parameter,private     ::      PI =  3.14159265359d0      
        real(kind=real64),parameter,private     ::      HBAR = 0.6582119569d0       !   Dirac constant eV.fs
            
       
        real(kind=real64),dimension(3,2),parameter,private      ::  MOTIF_BCC = reshape( (/ 0.0d0,0.0d0,0.0d0 , 0.5d0,0.5d0,0.5d0 /) , (/3,2/) )
        real(kind=real64),dimension(3,4),parameter,private      ::  MOTIF_FCC = reshape( (/ 0.0d0,0.0d0,0.0d0 , 0.5d0,0.5d0,0.0d0 , 0.5d0,0.0d0,0.5d0 , 0.0d0,0.5d0,0.5d0/) , (/3,4/) )
        real(kind=real64),dimension(3,2),parameter,private      ::  MOTIF_HCP = reshape( (/ 0.0d0,0.0d0,0.0d0 , 1.0d0/3,2.0d0/3,0.5d0 /),(/3,2/) )  
              
        real(kind=real64),dimension(3,3),parameter,private      ::  CONVENTIONAL_BCC = reshape( (/  1,0,0 , 0,1,0  , 0,0,1  /) , (/3,3/) )
        real(kind=real64),dimension(3,3),parameter,private      ::  CONVENTIONAL_FCC = reshape( (/  1,0,0 , 0,1,0  , 0,0,1  /) , (/3,3/) )
        real(kind=real64),dimension(3,3),parameter,private      ::  CONVENTIONAL_HCP = reshape( (/  0.5d0,-sqrt(0.75d0),0.0d0 , 0.5d0,sqrt(0.75d0),0.0d0 , 0.0d0,0.0d0,sqrt(8.0d0/3)  /) , (/3,3/) )  

    !---    
    
        public      ::      crystalStructureFactor
                                                                                                                                                                                
        interface           crystalStructureFactor
        !    module procedure            crystalStructureFactor0
            module procedure            crystalStructureFactor1
        end interface    
        
    contains
!---^^^^^^^^

        pure function realDoyleTurnerCoeff( a,b, dwf, s ) result( f )
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      find the Doyle-Turner coefficient
    !*          f = sum_i a_i Exp[ - (b_i + dwf) s^2 ]
    !*      for real input coefficients a_i,b_i,dwf.
    !*      note that s = g/(4 pi)
    
            real(kind=real64),dimension(:),intent(in)       ::      a,b
            real(kind=real64),intent(in)                    ::      dwf
            real(kind=real64),intent(in)                    ::      s
            
            real(kind=real64)                               ::      f
            
            integer             ::      ii
            real(kind=real64)   ::      xx
            
            f = 0.0d0
            do ii = 1,size(a)
                xx = (b(ii) + dwf)*s*s                
                f  = f + a(ii)*exp( - xx )
            end do
        
            return
        end function realDoyleTurnerCoeff
                    
            
        pure function complexDoyleTurnerCoeff( a,b, dwf, s ) result( f )            
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^        
    !*      find the Doyle-Turner coefficient
    !*          f = sum_i Re(a)_i Exp[ - (Re(a)_i + dwf) s^2 ] 
    !*            + i sum_i Im(a)_i Exp[ - (Im(a)_i + dwf/2 ) s^2 ] 
    !*      for complex input coefficients a_i,b_i 
    !*      and real Debye-Waller factor dwf
    !*      note that s = g/(4 pi)
    
            complex(kind=real64),dimension(:),intent(in)    ::      a,b
            real(kind=real64),intent(in)                    ::      dwf
            real(kind=real64),intent(in)                    ::      s
            
            complex(kind=real64)                            ::      f
            
            real(kind=real64)               ::      rDTc,iDTc
            
            rDTc = realDoyleTurnerCoeff( real(a) ,real(b) , dwf  , s )
            iDTc = realDoyleTurnerCoeff( aimag(a),aimag(b), dwf/2, s )
            f = complex( rDTc,iDTc )
               
        
            return
        end function complexDoyleTurnerCoeff
                    
        
    !     function crystalStructureFactor0( omega , a,b,dwf,g, r ) result (Vg)
    ! !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    ! !*      assuming a unit cell volume omega  (A^3) and atoms at r
    ! !*      and g-vector g
    ! !*      find the crystal structure factor from the Doyle Turner Coefficients
    ! !*      ( same for all atoms )
                    
    !         real(kind=real64),intent(in)                    ::      omega 
    !         complex(kind=real64),dimension(:),intent(in)    ::      a,b
    !         real(kind=real64),intent(in)                    ::      dwf
    !         real(kind=real64),dimension(3),intent(in)       ::      g
    !         real(kind=real64),dimension(:,:),intent(in)     ::      r
    !         complex(kind=real64)                            ::      Vg
            
    !         real(kind=real64)           ::      ss,gdotr
    !         complex(kind=real64)        ::      cDTc,expgr 
    !         integer                     ::      ii
            
    !         ss = norm2(g)/(4*PI)           
    !         cDTc = complexDoyleTurnerCoeff( a,b, dwf, ss )
            
    !         expgr =  0.0d0
    !         do ii = 1,size(r,dim=2)
    !             gdotr = g(1)*r(1,ii) + g(2)*r(2,ii) + g(3)*r(3,ii)
    !             expgr = expgr + exp( complex( 0.0d0,-gdotr ) )
    !             !print *,"Lib_CrystalStructureFactor::crystalStructureFactor0 info - g.r ",gdotr," g.r/pi ",gdotr/PI
    !         end do
            
    !         !print *,"Lib_CrystalStructureFactor::crystalStructureFactor0 info - sum exp[ -g.r ] ",expgr
    !         ss = (HBAR*HBAR/(2*ME))*(4*PI/Omega)
                        
    !         Vg = ss * expgr * cDTc
            
    !         return
    !     end function crystalStructureFactor0
        
        
        
        function crystalStructureFactor1( a0, latticename, hkl, a,b,dwf ) result (Vg)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      find the crystal structure factor from the Doyle Turner Coefficients
    !*      on input g is in reduced coordinates eg [110] and converted to A^-1 
                    
            real(kind=real64),dimension(3),intent(in)       ::      a0
            character(len=*),intent(in)                     ::      latticename
            complex(kind=real64),dimension(:,:),intent(in)  ::      a,b         !   (5,nMotif)
            real(kind=real64),intent(in)                    ::      dwf
            integer,dimension(3),intent(in)                 ::      hkl
            complex(kind=real64)                            ::      Vg
            
            type(Lattice)                                   ::      latt
            real(kind=real64)                               ::      omega 
            real(kind=real64)                               ::      ss , gdotr
            real(kind=real64),dimension(3)                  ::      rr
            real(kind=real64),dimension(3)                  ::      gg
            real(kind=real64),dimension(3,3)                ::      bb
            complex(kind=real64)        ::      cDTc,expgr 
            integer         ::      nMotif
            integer         ::      ii
            
            if (getLatticeType(latticename) == LATTICE_CUSTOM) then
                print *,"Lib_CrystalStructureFactor::crystalStructureFactor1 error - did not recognise lattice """//trim(latticename)//""""
                print *,"   only coded"
                call listAvailableLattices()
                stop
            end if            
            latt = Lattice_ctor(latticename)
            if (getLatticeType(latt) == LATTICE_HCP) call setCoverA(latt,a0(3)/a0(1))

        !---    find the reciprocal lattice vectors             
            bb = getReciprocalLatticeVectors(latt) / a0(1)
        
        !---    find the reflection considered
            gg = bb(:,1)*hkl(1) + bb(:,2)*hkl(2) + bb(:,3)*hkl(3) 

        !---    find the s parameter, so that V_g = V_g^0 Exp[ -B s^2 ] where B is the Debye-Waller factor   
        !       Acta Cryst. (1996). A52, 456-470 Eqn 1     
            ss = norm2(gg)/(4*PI)           
                
        !---    find the motif of the lattice, sum the crystal structure factor
        !       Acta Cryst. (1996). A52, 456-470 Eqn 13    
            nMotif = getNmotif(latt)
            Vg=0
            do ii = 1,nMotif
                rr(:) = a0(1) * getMotif(latt,ii)
                gdotr = dot_product( gg,rr )
                cDTc = complexDoyleTurnerCoeff( a(:,ii),b(:,ii), dwf, ss )          
                expgr = complex( cos(gdotr),-sin(gdotr) )                           !   Exp( -i g.r )
                Vg = Vg + expgr * cDTc
            end do



        !---    find the volume per atom and scale Vg
            omega = getOmega0(latt) * a0(1)**3

            Vg = Vg * (HBAR*HBAR/(2*ME))*(4*PI/Omega)


            print *,"Lib_CrystalStructureFactor::crystalStructureFactor1 info"
            ! print *,"   DT coeff a      ",a
        !    ! print *,"   DT coeff b      ",b
            print *,"   g-vector [hkl]  ",hkl
            print *,"   g-vector (A^-1) ",gg 
            print *,"   omega    (A^3)  ",omega 
            
           

        ! !---    find motif and volume of unit cell and g-vector in A^-1 based on lattice            
        !     select case( trim(lattice) )
        !         case( "bcc" )
        !             aa = CONVENTIONAL_BCC*a0
        !             call inverse3Mat(aa,ia)
        !             allocate(rr(3,size(MOTIF_BCC,dim=2)))
        !             rr(:,:) = MOTIF_BCC(:,:)
        !         case( "fcc" )
        !             aa = CONVENTIONAL_FCC*a0
        !             call inverse3Mat(aa,ia)
        !             allocate(rr(3,size(MOTIF_FCC,dim=2)))
        !             rr(:,:) = MOTIF_FCC(:,:)
        !         case( "hcp" )
        !             aa = CONVENTIONAL_HCP*a0
        !             call inverse3Mat(aa,ia)
        !             allocate(rr(3,size(MOTIF_HCP,dim=2)))
        !             rr(:,:) = MOTIF_HCP(:,:)
        !         case default
        !             print *,"Lib_CrystalStructureFactor::crystalStructureFactor1 error - did not recognise lattice """//trim(lattice)//""""
        !             Vg = 0.0d0
        !             return
        !     end select 
            
        !     do ii = 1,3
        !        rr(:,ii) = aa(:,1)*rr(1,ii) + aa(:,2)*rr(2,ii) + aa(:,3)*rr(3,ii)  
        !     end do
             
        !     ia = 2*PI*transpose(ia)
        !  !   print *,"Lib_CrystalStructureFactor::crystalStructureFactor1 info - lattice points "     
        !  !   write (*,fmt='(100f12.6)') rr(1,:)
        !  !   write (*,fmt='(100f12.6)') rr(2,:)
        !  !   write (*,fmt='(100f12.6)') rr(3,:)       
        !  !   
        !  !   print *,"Lib_CrystalStructureFactor::crystalStructureFactor1 info - lattice vectors "            
        !  !   write (*,fmt='(3f12.6)') aa(1,:)
        !  !   write (*,fmt='(3f12.6)') aa(2,:)
        !  !   write (*,fmt='(3f12.6)') aa(3,:)
        !  !   
        !  !   print *,"Lib_CrystalStructureFactor::crystalStructureFactor1 info - reciprocal lattice vectors "            
        !  !   write (*,fmt='(3f12.6)') ia(1,:)/(2*PI)
        !  !   write (*,fmt='(3f12.6)') ia(2,:)/(2*PI)
        !  !   write (*,fmt='(3f12.6)') ia(3,:)/(2*PI)
        !  !   
        !  !   print *,"dot prods "
        !  !   write (*,fmt='(3f12.6)') dot_product( aa(:,1),ia(:,1) )/(2*PI),dot_product( aa(:,1),ia(:,2) )/(2*PI),dot_product( aa(:,1),ia(:,3) )/(2*PI)
        !  !   write (*,fmt='(3f12.6)') dot_product( aa(:,2),ia(:,1) )/(2*PI),dot_product( aa(:,2),ia(:,2) )/(2*PI),dot_product( aa(:,2),ia(:,3) )/(2*PI)
        !  !   write (*,fmt='(3f12.6)') dot_product( aa(:,3),ia(:,1) )/(2*PI),dot_product( aa(:,3),ia(:,2) )/(2*PI),dot_product( aa(:,3),ia(:,3) )/(2*PI)
            
        !     omega  = determinant3Mat(aa)
        !     gg(1:3) = ia(1:3,1)*g(1) + ia(1:3,2)*g(2) + ia(1:3,3)*g(3) 
        !     gg(1:3) = gg(1:3)
            
        !     print *,"Lib_CrystalStructureFactor::crystalStructureFactor1 info"
        !    ! print *,"   DT coeff a      ",a
        !    ! print *,"   DT coeff b      ",b
        !     print *,"   g-vector [hkl]  ",g
        !     print *,"   g-vector (A^-1) ",gg/(2*PI)
        !     print *,"   omega    (A^3)  ",omega 
            
        !     Vg = crystalStructureFactor0( omega , a,b,dwf,gg, rr ) 
            
            return
        end function crystalStructureFactor1
        
        
        

        pure subroutine inverse3Mat(M,N)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      monadic operator
    !*      finds the inverse of a general three matrix
            real(kind=real64),dimension(3,3),intent(in)       ::  M
            real(kind=real64),dimension(3,3),intent(out)      ::  N
            real(kind=real64)            ::      idd
            real(kind=real64),dimension(3,3),parameter        ::      &
            IDENTITY3MAT = reshape( (/ 1.0d0,0.0d0,0.0d0 , 0.0d0,1.0d0,0.0d0 , 0.0d0,0.0d0,1.0d0 /) &
                                   ,(/ 3,3 /) )

            idd = determinant3Mat(M)
            if (abs(idd) < tiny(1.0d0)) then
                N = IDENTITY3MAT
                return
            end if
            idd = 1.0/idd

            N(1,1)   = ( M(2,2)*M(3,3) - M(2,3)*M(3,2) ) * idd
            N(2,1)   = ( M(2,3)*M(3,1) - M(2,1)*M(3,3) ) * idd
            N(3,1)   = ( M(2,1)*M(3,2) - M(2,2)*M(3,1) ) * idd

            N(1,2)   = ( M(1,3)*M(3,2) - M(1,2)*M(3,3) ) * idd
            N(2,2)   = ( M(1,1)*M(3,3) - M(1,3)*M(3,1) ) * idd
            N(3,2)   = ( M(1,2)*M(3,1) - M(1,1)*M(3,2) ) * idd

            N(1,3)   = ( M(1,2)*M(2,3) - M(1,3)*M(2,2) ) * idd
            N(2,3)   = ( M(1,3)*M(2,1) - M(1,1)*M(2,3) ) * idd
            N(3,3)   = ( M(1,1)*M(2,2) - M(1,2)*M(2,1) ) * idd

            return
        end subroutine inverse3Mat

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

        
        
    end module Lib_CrystalStructureFactor     
    
    
