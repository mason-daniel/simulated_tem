
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
        use Lib_RelativisticElectrons
        use Lib_inputReal
        use Lib_inputTDS
        use Lib_FitDoyleTurner
        use iso_fortran_env
        implicit none
        private

        character(len=8),public,parameter       ::      LIB_CRYSTALSTRUCTUREFACTOR_VERSION = "0.0.1"        
        real(kind=real64),parameter,private     ::      PI =  3.14159265359d0      
            
        
    !---    
    
        public      ::      computeCrystalStructureFactors
        public      ::      crystalStructureFactor
  
                                                                                                                                                                                
        interface           crystalStructureFactor
            module procedure            crystalStructureFactor1
        end interface    
        
    contains
!---^^^^^^^^
 
            
        subroutine computeCrystalStructureFactors( element, hkl, T,V ,Vg , a,b)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      find the crystal structure factor from the element name and reflection at the given temperature and accelerating voltage
                    
            character(len=*),intent(in)                     ::      element
            integer,dimension(:,:),intent(in)               ::      hkl         !   (3,nReflections)
            real(kind=real64),intent(in)                    ::      T       
            real(kind=real64),intent(in)                    ::      V
            complex(kind=real64),dimension(:),intent(out)   ::      Vg          !   (1:nReflections)
            complex(kind=real64),dimension(:),intent(out),optional      ::      a,b

            character(len=LATTNAME_LEN)                     ::      latticename
            real(kind=real64)                               ::      dwf
            real(kind=real64),dimension(3)                  ::      a0
            type(Lattice)                                   ::      latt
            integer                                         ::      nMotif,nReflections,ii
            

            integer                                         ::      nS
            real(kind=real64)                               ::      smin,smax,ds        
            real(kind=real64),dimension(DT_COEFFS)          ::      ar,br,ai,bi             !   Doyle-Turner coefficients
            real(kind=real64),dimension(:),allocatable      ::      ff,ss           !   Scattering amplitudes
            real(kind=real64)                               ::      mser,msei
            complex(kind=real64),dimension(:,:),allocatable ::      aa,bb

        !---    find the lattice
            latticename =  getLatticeName(element) 
            if (getLatticeType(latticename) == LATTICE_CUSTOM) then
                print *,"Lib_CrystalStructureFactor::crystalStructureFactor0 error - did not recognise lattice """//trim(latticename)//""""
                print *,"   only coded"
                call listAvailableLattices()
                stop
            end if            
            latt = Lattice_ctor(latticename)
            nMotif = getConventionalnMotif(latt)
            allocate(aa(DT_COEFFS,nMotif))
            allocate(bb(DT_COEFFS,nMotif))         
            print *,"Lib_CrystalStructureFactor::crystalStructureFactor0 info -"
            call report(latt)

        !---    Debye Waller factor and lattice constant
            dwf = DebyeWallerFactor(element,T,latticename)
            a0 = getLatticeConstant(element)          
            write(*,fmt='(a,f12.6,a)')  " Lib_CrystalStructureFactor::crystalStructureFactor0 info - Debye-Waller factor ",dwf," (A^2)"
            write(*,fmt='(a,3f12.6,a)') " Lib_CrystalStructureFactor::crystalStructureFactor0 info - lattice parameters  ",a0," (A)"

        !---    find real DT coeffs
            smin = LIB_INPUTREAL_SMIN
            smax = LIB_INPUTREAL_SMAX
            nS = LIB_INPUTREAL_NS
            ds = (smax-smin)/nS
            allocate(ss(0:ns))
            allocate(ff(0:nS))
            do ii = 0,nS
                ss(ii) = smin + ii*ds
            end do
            call getScatteringAmplitude(element,ss,ff)
            ar = getrealDoyleTurner_a(element)
            br = getrealDoyleTurner_b(element)
            call fitDoyleTurner( ss,ff,ar,br, mser )
            deallocate(ss)
            deallocate(ff)
           ! print *,"Lib_CrystalStructureFactor::crystalStructureFactor0 info - Doyle-Turner coeffs for real part of scattering "
            ! write(*,fmt='(a6,2a16)') "coeff"," a_r(ii)"," b_r(ii)"
            ! do ii = 1,DT_COEFFS
            !     write(*,fmt='(i6,2f16.8))') ii,ar(ii),br(ii)
            ! end do

            
        !---    find imag DT coeffs
            smin = LIB_INPUTTDS_SMIN
            smax = LIB_INPUTTDS_SMAX
            nS = LIB_INPUTTDS_NS
            ds = (smax-smin)/nS
            allocate(ss(0:ns))
            allocate(ff(0:nS))
            call getImagScatteringAmplitude(element,V,dwf,ss,ff)
            ai = getscatDoyleTurner_a(element,T)
            bi = getscatDoyleTurner_b(element,T)
            call fitDoyleTurner( ss,ff,ai,bi, msei )
            deallocate(ss)
            deallocate(ff)
            print *,"Lib_CrystalStructureFactor::crystalStructureFactor0 info - Doyle-Turner coeffs "
            write(*,fmt='(a6,4a16)') "coeff"," a_r(ii)"," b_r(ii)"," a_i(ii)"," b_i(ii)"
            do ii = 1,DT_COEFFS
                write(*,fmt='(i6,4f16.8)') ii,ar(ii),br(ii),ai(ii),bi(ii)
            end do
            write(*,fmt='(3(a,f12.8))')" Lib_CrystalStructureFactor::crystalStructureFactor0 info - mean square error ",mser," (real) ",msei," (imag)"
            
        !---    construct complex DT coeffs
            
            do ii = 1,DT_COEFFS
                aa(ii,:) = complex( ar(ii),ai(ii) )
                bb(ii,:) = complex( br(ii),bi(ii) )
            end do

            nReflections = size(hkl,dim=2)
            if (size(hkl,dim=1)==3) then
                !   hkl given as Miller indices
                do ii = 1,nReflections
                    Vg(ii) = crystalStructureFactor1( a0, latticename, hkl(:,ii), aa,bb,dwf )
                end do
            else
                !   hkl given as Miller-Bravais indices
                do ii = 1,nReflections
                    Vg(ii) = crystalStructureFactor1( a0, latticename, MillerBravaisToMiller_plane(hkl(:,ii)), aa,bb,dwf )
                end do
            end if


            if (present(a)) then
                do ii = 1,DT_COEFFS
                    a(ii) = complex( ar(ii),ai(ii) )
                end do
            end if
            if (present(b)) then
                do ii = 1,DT_COEFFS
                    b(ii) = complex( br(ii),bi(ii) )
                end do
            end if

            return
        end subroutine computeCrystalStructureFactors

        function crystalStructureFactor1( a0, latticename, hkl, a,b,dwf ) result (Vg)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      find the crystal structure factor from the Doyle Turner Coefficients
    !*      on input g is in reduced coordinates eg [110] and converted to A^-1 
                    
            real(kind=real64),dimension(3),intent(in)       ::      a0          !   for cubic crystals only use a0(1). For hcp will also set c/a ratio.
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
            complex(kind=real64)        ::      cDT,expgr 
            integer         ::      nMotif,ii
             
            
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
            nMotif = getConventionalnMotif(latt)
            Vg = 0
            do ii = 1,nMotif
                rr(:) = a0(1) * getConventionalMotif( latt,ii,removeOffset=.true. )
                gdotr = dot_product( gg,rr )
                !print *,"motif point ",ii,rr," g.r ",gdotr," exp[ -i g.r ]",exp( complex(0,-gdotr) )
                cDT = complexDoyleTurnerSum( a(:,ii),b(:,ii), dwf, ss )          
                expgr = complex( cos(gdotr),-sin(gdotr) )                           !   Exp( -i g.r )
                Vg = Vg + expgr * cDT
            end do



        !---    find the volume per conventional cell and scale Vg
            omega = getOmega0(latt) * a0(1)**3 * getConventionalnMotif( latt )

            Vg = Vg * (HBAR*HBAR/(2*ME))*(4*PI/Omega)


            ! print *,"Lib_CrystalStructureFactor::crystalStructureFactor1 info"
            ! print *,"   g-vector [hkl]  ",hkl
            ! print *,"   g-vector (A^-1) ",gg 
            ! print *,"   omega    (A^3)  ",omega 
            
           
 
            return
        end function crystalStructureFactor1
        
         
                    
            
        pure function complexDoyleTurnerSum( a,b, dwf, s ) result( f )            
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^        
    !*      find the Doyle-Turner sum with complex coefficients
    !*          f = sum_i Re(a)_i Exp[ - (Re(a)_i + dwf) s^2 ] 
    !*            + i sum_i Im(a)_i Exp[ - (Im(a)_i + dwf/2 ) s^2 ] 
    !*      for complex input coefficients a_i,b_i 
    !*      and real Debye-Waller factor dwf
    !*      note that s = g/(4 pi)
    
            complex(kind=real64),dimension(:),intent(in)    ::      a,b     !   complex Doyle-Turner coefficient
            real(kind=real64),intent(in)                    ::      dwf     !   Debye-Waller factor
            real(kind=real64),intent(in)                    ::      s       !   scattering angle s = g/4pi = sin(theta)/lambda
            
            complex(kind=real64)                            ::      f
            
            real(kind=real64)               ::      rDTc,iDTc
            
            rDTc = DoyleTurnerSum( s, real(a) ,real(b) + dwf    )
            iDTc = DoyleTurnerSum( s, aimag(a),aimag(b) + dwf/2 )
            f = complex( rDTc,iDTc )
        
            return
        end function complexDoyleTurnerSum
                    
         
        
    end module Lib_CrystalStructureFactor     
    
    
