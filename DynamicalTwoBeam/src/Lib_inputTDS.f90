! ******************************************************************************
! 
!                   FORTRAN SOURCE FILE       inputTDS.f
!
! *******************************************************************************
!                                                        by S. L. Dudarev, 
!                                                        Department of Materials,
!                                                               University of Oxford
!                                                                Parks Road
!                                                                OXFORD OX1 3PH
!                                                                UK
!
!                        present address:                EURATOM/UKAEA Fusion Association
!                                                        Culham Science Centre
!                                                        Oxfordshire OX14 3DB
!                                                        UK
!
! application of this program has been described in the following publications:
!
!  S L Dudarev and M J Whelan, Surface Science 310 (1994) 373 
!  S L Dudarev, D D Vvedensky and M J Whelan, Phys. Review B50 (1994) 14525 
!  S L Dudarev, L M Peng and M J Whelan, Surface Science 330 (1995) 86
!  see also
!  L M Peng, S L Dudarev and M J Whelan, High-Energy Electron Diffraction and Microscopy
! (Oxford University Press, 2004) ISBN 0-19-850074-2
! ******************************************************************************** 
!   Edit by D.R. Mason 28/01/21
!       changed some obsolete numbered do loops
!       removed obsolete pause command
!       tweaked filenames to allow directory paths
!       no functional changes

!   DRM Sept 2024
!   change to f08 module

    module Lib_inputTDS
!---^^^^^^^^^^^^^^^^^^^
        use iso_fortran_env
        use Lib_Elements
        use Lib_SimpleProgressBar
        use Lib_inputREAL
        implicit none
        private

        real(kind=real64),public,parameter     ::      LIB_INPUTTDS_SMIN = 0.0d0          !   minimum scattering s = q/4pi
        real(kind=real64),public,parameter     ::      LIB_INPUTTDS_SMAX = 15.0d0         !   maximum scattering s = q/4pi
        integer,public,parameter               ::      LIB_INPUTTDS_NS = 150              !   number of s intervals
        real(kind=real64),parameter ::      ABS_ERR_TOL = 1.0d-7         !   absolute error tolerance
        real(kind=real64),parameter ::      REL_ERR_TOL = 1.0d-5         !   relative error tolerance

    !---    Gaussian quadrature with n=16 points - note that we use abscissa +/- x_i and weights w_i
    !       see eg https://pomax.github.io/bezierinfo/legendre-gauss.html
        real(kind=real64),private,dimension(16),parameter     ::          GQ_abscissa =      &
                        (/  -.98940093499164993D0, -.94457502307323258D0,       &
                            -.86563120238783174D0, -.75540440835500303D0,       &
                            -.61787624440264375D0, -.45801677765722739D0,       &
                            -.28160355077925891D0, -.09501250983763744D0,       &
                            0.09501250983763744D0, 0.28160355077925891D0,       &
                            0.45801677765722739D0, 0.61787624440264375D0,       &
                            0.75540440835500303D0, 0.86563120238783174D0,       &
                            0.94457502307323258D0, 0.98940093499164993D0        /)

        real(kind=real64),private,dimension(16),parameter     ::          GQ_weight =      &
                        (/  0.02715245941175409D0, 0.06225352393864789D0,       &
                            0.09515851168249278D0, 0.12462897125553387D0,       &
                            0.14959598881657673D0, 0.16915651939500254D0,       &
                            0.18260341504492359D0, 0.18945061045506850D0,       &
                            0.18945061045506850D0, 0.18260341504492359D0,       &
                            0.16915651939500254D0, 0.14959598881657673D0,       &
                            0.12462897125553387D0, 0.09515851168249278D0,       &
                            0.06225352393864789D0, 0.02715245941175409D0    /)


        public      ::      getImagScatteringAmplitude

    contains
!---^^^^^^^^

        subroutine getImagScatteringAmplitude(el,V,bdw,fimag)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ 
            character(len=*),intent(in)                     ::      el      !   element name
            real(kind=real64),intent(in)                    ::      V       !   accelerating voltage (V)
            real(kind=real64),intent(in)                    ::      bdw     !   debye-waller factor (A^2)
            real(kind=real64),dimension(0:),intent(out)     ::      fimag
            real(kind=real64)           ::      alamda,coeff
            
            real(kind=real64),parameter ::      A = -20.0d0         !   estimated lower bound of integral
            real(kind=real64),parameter ::      B = 20.0d0          !   estimated upper bound of integral

            real(kind=real64)           ::      smax,smin,ds,ss     !   range, interval, s parameter
            integer                     ::      nS                  !   number of s params
            integer                     ::      ee                  !   element index
            integer                     ::      ii

            real(kind=real64)           ::      sum
 
            alamda=12.264306d0/dsqrt(V)/dsqrt(1.d0+0.97846707d-06*V)
            coeff=2.d0*alamda*(1.d0+V/0.511d06)
            ee = whichElement(el)
            if (ee==0) then 
                print *,"Lib_InputTDS::getImagScatteringAmplitude error - did not recognise element """//trim(el)//""""
                fimag = 0
                return
            end if

            smin = LIB_INPUTTDS_SMIN
            smax = LIB_INPUTTDS_SMAX
            nS = LIB_INPUTTDS_NS
            ds = (smax-smin)/nS

            do ii = 0,nS

                call progressBar(ii+1,nS+1)
                ss = smin + ds*ii           !   s = g/4pi = sin(theta)/lambda 
                                            !   where theta is the scattering angle and lambda the electron wavelength

        ! two-dimensional integration using gaussian quadratures
                call gauss2( ee,bdw,ss,A,B,sum )

                fimag(ii)=sum*coeff*dexp(-0.5d0*ss*ss*bdw)
    
            end do
            return
        end subroutine getImagScatteringAmplitude


! ***************************************************************************

!     subroutine fx(x,w)
! !---^^^^^^^^^^^^^^^^^^
!     !implicit real*8 (a-h,o-z)
!         real(kind=real64),intent(in)        ::      x
!         real(kind=real64),intent(out)       ::      w
!         w=1.d0
!         return
!     end subroutine fx

! ***************************************************************************

    real(kind=real64) function fxy(e,B,s,x,y)
!---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
!*      given the Debye-Waller factor B and element index e
!*      and scattering parameter s = q/4pi
!*      compute the value of the function fxy at point x,y
        integer,intent(in)                  ::      e
        real(kind=real64),intent(in)        ::      B
        real(kind=real64),intent(in)        ::      s
        real(kind=real64),intent(in)        ::      x
        real(kind=real64),intent(in)        ::      y
         
        real(kind=real64)       ::      x1,x2
        real(kind=real64)       ::      xq,yq
        real(kind=real64)       ::      ar1,ar2,are,ex
        x1=x-0.5d0*s
        x2=x+0.5d0*s
        xq=x*x
        yq=y*y
        ar1=sqrt(yq+x1*x1)
        ar2=sqrt(yq+x2*x2)
        are=2*B*(s*s/4-xq-yq)
        ex=exp(are)

        fxy = scatteringAmplitude(e,ar1)*scatteringAmplitude(e,ar2)*(1.d0-ex)

        return
    end function fxy


!     subroutine fxy(el,bdw,s,x,y,w)
! !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
! !*      given the Debye-Waller factor bdw and element name
! !*      and scattering parameter s = q/4pi
!         character(len=*),intent(in)         ::      el
!         real(kind=real64),intent(in)        ::      bdw
!         real(kind=real64),intent(in)        ::      s
!         real(kind=real64),intent(in)        ::      x
!         real(kind=real64),intent(in)        ::      y
!         real(kind=real64),intent(out)       ::      w
!     !   implicit real*8 (a-h,o-z)
!     !character*2 elname,elnames(50)
!     !common s,bdw,elname,elnames
!         real(kind=real64)       ::      x1,x2
!         real(kind=real64)       ::      xq,yq
!         real(kind=real64)       ::      ar1,ar2,are,ex
!         x1=x-0.5d0*s
!         x2=x+0.5d0*s
!         xq=x*x
!         yq=y*y
!         ar1=dsqrt(yq+x1*x1)
!         ar2=dsqrt(yq+x2*x2)
!         are=2*bdw*(s*s/4-xq-yq)
!         ex=dexp(are)

!         w = scatteringAmplitude(el,ar1)*scatteringAmplitude(el,ar2)*(1.d0-ex)

!     ! if(elname.eq.elnames(1)) W=FLI(ar1)*FLI(ar2)*(1.d0-ex)
!     ! if(elname.eq.elnames(2)) W=FBE(ar1)*FBE(ar2)*(1.d0-ex)
!     ! if(elname.eq.elnames(3)) W=FC(ar1)*FC(ar2)*(1.d0-ex)
!     ! if(elname.eq.elnames(4)) W=FO(ar1)*FO(ar2)*(1.d0-ex)
!     ! if(elname.eq.elnames(5)) W=FNA(ar1)*FNA(ar2)*(1.d0-ex)
!     ! if(elname.eq.elnames(6)) W=FMG(ar1)*FMG(ar2)*(1.d0-ex)
!     ! if(elname.eq.elnames(7)) W=FAL(ar1)*FAL(ar2)*(1.d0-ex)
!     ! if(elname.eq.elnames(8)) W=FSI(ar1)*FSI(ar2)*(1.d0-ex)
!     ! if(elname.eq.elnames(9)) W=FK(ar1)*FK(ar2)*(1.d0-ex)
!     ! if(elname.eq.elnames(10)) W=FCA(ar1)*FCA(ar2)*(1.d0-ex)
!     ! if(elname.eq.elnames(11)) W=FTI(ar1)*FTI(ar2)*(1.d0-ex)
!     ! if(elname.eq.elnames(12)) W=FV(ar1)*FV(ar2)*(1.d0-ex)
!     ! if(elname.eq.elnames(13)) W=FCR(ar1)*FCR(ar2)*(1.d0-ex)
!     ! if(elname.eq.elnames(14)) W=FMN(ar1)*FMN(ar2)*(1.d0-ex)
!     ! if(elname.eq.elnames(15)) W=FFE(ar1)*FFE(ar2)*(1.d0-ex)
!     ! if(elname.eq.elnames(16)) W=FCO(ar1)*FCO(ar2)*(1.d0-ex)
!     ! if(elname.eq.elnames(17)) W=FNI(ar1)*FNI(ar2)*(1.d0-ex)
!     ! if(elname.eq.elnames(18)) W=FCU(ar1)*FCU(ar2)*(1.d0-ex)
!     ! if(elname.eq.elnames(19)) W=FZN(ar1)*FZN(ar2)*(1.d0-ex)
!     ! if(elname.eq.elnames(20)) W=FGA(ar1)*FGA(ar2)*(1.d0-ex)
!     ! if(elname.eq.elnames(21)) W=FGE(ar1)*FGE(ar2)*(1.d0-ex)
!     ! if(elname.eq.elnames(22)) W=FAS(ar1)*FAS(ar2)*(1.d0-ex)
!     ! if(elname.eq.elnames(23)) W=FRB(ar1)*FRB(ar2)*(1.d0-ex)
!     ! if(elname.eq.elnames(24)) W=FSR(ar1)*FSR(ar2)*(1.d0-ex)
!     ! if(elname.eq.elnames(25)) W=FY(ar1)*FY(ar2)*(1.d0-ex)
!     ! if(elname.eq.elnames(26)) W=FZR(ar1)*FZR(ar2)*(1.d0-ex)
!     ! if(elname.eq.elnames(27)) W=FNB(ar1)*FNB(ar2)*(1.d0-ex)
!     ! if(elname.eq.elnames(28)) W=FMO(ar1)*FMO(ar2)*(1.d0-ex)
!     ! if(elname.eq.elnames(29)) W=FRU(ar1)*FRU(ar2)*(1.d0-ex)
!     ! if(elname.eq.elnames(30)) W=FRH(ar1)*FRH(ar2)*(1.d0-ex)    
!     ! if(elname.eq.elnames(31)) W=FPD(ar1)*FPD(ar2)*(1.d0-ex)
!     ! if(elname.eq.elnames(32)) W=FAG(ar1)*FAG(ar2)*(1.d0-ex)
!     ! if(elname.eq.elnames(33)) W=FCD(ar1)*FCD(ar2)*(1.d0-ex)
!     ! if(elname.eq.elnames(34)) W=FSN(ar1)*FSN(ar2)*(1.d0-ex)
!     ! if(elname.eq.elnames(35)) W=FCS(ar1)*FCS(ar2)*(1.d0-ex)
!     ! if(elname.eq.elnames(36)) W=FBA(ar1)*FBA(ar2)*(1.d0-ex)
!     ! if(elname.eq.elnames(37)) W=FLA(ar1)*FLA(ar2)*(1.d0-ex)
!     ! if(elname.eq.elnames(38)) W=FCE(ar1)*FCE(ar2)*(1.d0-ex)
!     ! if(elname.eq.elnames(39)) W=FGD(ar1)*FGD(ar2)*(1.d0-ex)
!     ! if(elname.eq.elnames(40)) W=FHF(ar1)*FHF(ar2)*(1.d0-ex)    
!     ! if(elname.eq.elnames(41)) W=FTA(ar1)*FTA(ar2)*(1.d0-ex)
!     ! if(elname.eq.elnames(42)) W=FW(ar1)*FW(ar2)*(1.d0-ex)
!     ! if(elname.eq.elnames(43)) W=FRE(ar1)*FRE(ar2)*(1.d0-ex)
!     ! if(elname.eq.elnames(44)) W=FOS(ar1)*FOS(ar2)*(1.d0-ex)
!     ! if(elname.eq.elnames(45)) W=FIR(ar1)*FIR(ar2)*(1.d0-ex)
!     ! if(elname.eq.elnames(46)) W=FPT(ar1)*FPT(ar2)*(1.d0-ex)
!     ! if(elname.eq.elnames(47)) W=FAU(ar1)*FAU(ar2)*(1.d0-ex)
!     ! if(elname.eq.elnames(48)) W=FTL(ar1)*FTL(ar2)*(1.d0-ex)
!     ! if(elname.eq.elnames(49)) W=FPB(ar1)*FPB(ar2)*(1.d0-ex)    
!     ! if(elname.eq.elnames(50)) W=FTH(ar1)*FTH(ar2)*(1.d0-ex)    

!         return
!     end subroutine fxy



! ***************************************************************************

    subroutine ylim(x,ylim1,ylim2)
!---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
        real(kind=real64),intent(in)        ::      x
        real(kind=real64),intent(out)       ::      ylim1,ylim2
    !   implicit real*8 (a-h,o-z)
        ylim1=-dsqrt(400.d0-x*x)
        ylim2=-ylim1
        return
    end subroutine ylim
    
    SUBROUTINE GAUSS2(e,bdw,s,A,B,SUM )
!---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ 
        integer,intent(in)                  ::      e
        real(kind=real64),intent(in)        ::      s
        real(kind=real64),intent(in)        ::      bdw
        real(kind=real64),intent(in)        ::      A,B                         !   initial extimate for integral limits
        !real(kind=real64),intent(in)        ::      REL_ERR_TOL,ABS_ERR_TOL     !   relative, absolute error tolerance
         
        real(kind=real64),intent(out)       ::      SUM      
        

    !    real(kind=real64),dimension(8)       ::      xx,ww
        real(kind=real64),dimension(80)      ::      aa,bb,w
        real(kind=real64)       ::      A0,B0,A1,B1
        real(kind=real64)       ::      SSSS,SA,SB,S1,ZZ,SS,ZM  !,EA,ER
        real(kind=real64)       ::      W2!,W1,WW1,WW2
        real(kind=real64)       ::      X,YLIM1,YLIM2
        integer                 ::      L,ISUM,K

!    IMPLICIT REAL*8(A-H,O-Z)
!    EXTERNAL GAUSS1,YLIM,FX,FXY
    !DIMENSION XX(8),WW(8),AA(80),BB(80),W(80)
        ! DATA XX/ 0.09501250983763744D0, 0.28160355077925891D0,      &
        !         0.45801677765722739D0, 0.61787624440264375D0,      &
        !         0.75540440835500303D0, 0.86563120238783174D0,      &
        !         0.94457502307323258D0, 0.98940093499164993D0/
        ! DATA WW/ 0.18945061045506850D0, 0.18260341504492359D0,      &
        !         0.16915651939500254D0, 0.14959598881657673D0,      &
        !         0.12462897125553387D0, 0.09515851168249278D0,      &
        !         0.06225352393864789D0, 0.02715245941175409D0/
    SUM=0.D0
    ! ER=REL_ERR_TOL*0.1D0
    ! EA=AE*0.01D0
    L=1
    A0=A
    B0=B
    AA(1)=A0
    BB(1)=B0
    ISUM=1
    GO TO 100
  10    W(1)=SSSS
  25    A0=AA(L)
    B0=0.5D0*(BB(L)+A0)
    SS=W(L)
    ISUM=2
    GO TO 100
  20    A1=A0
    B1=B0
    S1=SSSS
    B0=B0+B0-A0
    A0=B1
    ISUM=3
    GO TO 100 
  30    SA=SSSS+S1
    SB=DABS(SA-SS)
    IF(SB.LT.REL_ERR_TOL*DABS(SA)) GO TO 40
    IF(SB.LT.REL_ERR_TOL*DABS(SUM)) GO TO 40
    IF(SB.LT.ABS_ERR_TOL) GO TO 40
    AA(L)=A1
    BB(L)=B1
    W(L)=S1
    L=L+1
    IF(L.GT.80) GO TO 50
    AA(L)=A0
    BB(L)=B0
    W(L)=SSSS
    GO TO 25
  40    SUM=SUM+SA
    L=L-1
    IF(L.EQ.0) RETURN
    GO TO 25
  50    WRITE(6,*) ' number of iterations = 80 = limiting value'
    STOP
 100    ZZ=0.5D0*(B0-A0)
    ZM=ZZ+A0
    SSSS=0.D0
    DO K=1,size(GQ_abscissa)
        X=ZM+ZZ*GQ_abscissa(K) !   XX(K)
        CALL YLIM(X,YLIM1,YLIM2)
!        CALL FX(X,WW2)     !   just sets ww2=1
        !CALL GaussianQuadratureIntegration_y(e,bdw,s,YLIM1,YLIM2,W2,X)
        call GG_Integration_y(YLIM1,YLIM2,fxy,e,bdw,s,W2,x)
        SSSS=SSSS+GQ_weight(K)*W2
!         X=ZM+ZM-X
!         CALL YLIM(X,YLIM1,YLIM2)
! !        CALL FX(X,WW1)     !   just sets ww1=1
!         CALL GaussianQuadratureIntegration_y(el,bdw,s,YLIM1,YLIM2,REL_ERR_TOL,AE,W1,X)
!        SSSS=SSSS+WW(K)*(W1*WW1+W2*WW2)
    !    SSSS=SSSS+WW(K)*(W1+W2)     !   equivalent
        ! SSSS=SSSS+GaussianQuadrature_w(K)*(W1+W2)
    end do
    
    SSSS=SSSS*ZZ
    GO TO (10,20,30),ISUM
    END SUBROUTINE GAUSS2

!------------------------------------------------------------------------------------

    subroutine GaussianQuadratureIntegration_y(e,bdw,s,A,B,integral,x)
!---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
!*      given the Debye-Waller factor bdw and element index e
!*      and scattering parameter s = q/4pi
!*      use Gaussian quadrature to integrate over y at x position X1
        integer,intent(in)                  ::      e
        real(kind=real64),intent(in)        ::      bdw
        real(kind=real64),intent(in)        ::      s
        real(kind=real64),intent(in)        ::      A,B         !   integration ranges
        !real(kind=real64),intent(in)        ::      REL_ERR_TOL,ABS_ERR_TOL       !   relative, absolute error tolerances
        real(kind=real64),intent(out)       ::      integral    !   the result, integral_a^b FXY(X,y) dy
        real(kind=real64),intent(in)        ::      x          !  x position

        ! real(kind=real64),dimension(8)       ::      xx,ww
        real(kind=real64),dimension(20)      ::      aa,bb,w
        real(kind=real64)       ::      A0,B0,A1,B1
        real(kind=real64)       ::      SSSS,SA,SB,S1,S2,S_total,SS!,EA,ER
        integer                 ::      ii 

 
        integral=0.D0       !   this is the result of the integration
        ii=1
        A0=A
        B0=B
        AA(1)=A0
        BB(1)=B0
        S_total = GaussianQuadrature( A0,B0,fxy,e,bdw,s,x )
        W(1)=S_total
        do

        !   compute integral on left half of the interval AA(ii):BB(ii)
            A0=AA(ii)
            B0=(BB(ii)+AA(ii))/2
            SS=W(ii)
            S1=GaussianQuadrature( A0,B0,fxy,e,bdw,s,x )
            A1=A0
            B1=B0
            
        !   compute integral on right half of the interval AA(ii):BB(ii)            
            A0=(BB(ii)+AA(ii))/2
            B0=BB(ii)
            S2 = GaussianQuadrature( A0,B0,fxy,e,bdw,s,x )
            S_total=S1+S2
            SB=abs(S_total-SS)

            if (SB < maxval( (/REL_ERR_TOL*abs(S_total),REL_ERR_TOL*abs(integral),ABS_ERR_TOL/) )) then
                integral=integral+S_total
                ii=ii-1
                if (ii == 0) return
                cycle
            end if

            AA(ii)=A1
            BB(ii)=B1
            W(ii)=S1
            ii=ii+1
            if (ii > 20) exit

            AA(ii)=A0
            BB(ii)=B0
            W(ii)=S2
        end do
        write(6,*) ' number of iterations = 20 = limiting value'
        return
    end subroutine GaussianQuadratureIntegration_y

    
    recursive subroutine GG_Integration_y(ymin,ymax,func,e,B,s,integral,x,level)
!---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
!*      given the Debye-Waller factor bdw and element index e
!*      and scattering parameter s = q/4pi
!*      use Gaussian quadrature to integrate over y at x position x
!*
!*      This subroutine works by comparing an input estimate to the integral 
!*          integral = int_ymin^ymax func(x,y) dy
!*      generated by a 16 point Gaussian Quadrature
!*      with an improved estimate made by splitting the range into two, ie two x 16 point Gaussian Quadratures
!*          integral' = int_ymin^ybar func(x,y) dy + int_ybar^ymax func(x,y) dy 
!*      with ybar = (ymin+ymax)/2      
!*      If the error is within tolerance, return. Otherwise converge the left and right hand sides recursively.
!*      
!*
!*      The recursive operation is handled by the optional parameter level
!*      If level is absent, then this must be the first call, we don't have an estimate for integral over (a:b).
!*      otherwise we have made an estimate, and found it to be inaccurate.
!*      The next call then increases level by one.
!*      Each recursive call bisects the interval, so quit when level is 20.
!*
        real(kind=real64),intent(in)        ::      ymin,ymax   !   integration ranges
        real(kind=real64),external          ::      func        !   function to integrate.
        integer,intent(in)                  ::      e           !   element index (See Lib_Elements)
        real(kind=real64),intent(in)        ::      B           !   Debye-Waller B parameter (See Lib_Elements)
        real(kind=real64),intent(in)        ::      s           !   scattering parameter s = q/4pi
        real(kind=real64),intent(inout)     ::      integral    !   the result, integral_ymin^ymax func(x,y) dy
        real(kind=real64),intent(in)        ::      x           !   x position
        integer,intent(in),optional         ::      level       !   controls number of recursion loops

      
        real(kind=real64)       ::      ybar
        real(kind=real64)       ::      integral_left,integral_right,integral_diff
 
        
    !---    decide what calculation is required        
        if (.not. present(level)) then
            !   first call to the subroutine, so I don't have an estimate for integral over whole interval = int_ymin,ymax func(x,y) dy
            integral = GaussianQuadrature( ymin,ymax,func,e,B,s,x )
        else if (level == 20) then
            !   I have now bisected the interval 20 times. Probably something is going wrong.
            print *,"Lib_InputTDS::GG_Integration_y warning - number of iterations = 20 = limiting value"
            return 
        else 
        !   we should bisect the interval, check if we have convergence now. I do have a previous estimate for integral over whole interval.
            !print *,"Lib_InputTDS::GG_Integration_y info - ",level,integral
        end if


    !---    bisect the interval, and compute the left and right integrals. 
    !       This gives a better estimate for integral over whole interval.
        ybar = (ymin + ymax)/2                                                  !   bisection point
        integral_left  = GaussianQuadrature( ymin,ybar,func,e,B,s,x )
        integral_right = GaussianQuadrature( ybar,ymax,func,e,B,s,x )



    !---    have we converged yet?        
        integral_diff  = integral - ( integral_left + integral_right )
        if (integral_diff >= max( ABS_ERR_TOL, REL_ERR_TOL*max(abs(integral),abs(integral_left + integral_right) ))) then
        !   no, we have not yet converged. Better try to converge the left and right sides.
            if (present(level)) then
                call GG_Integration_y(ymin,ybar,func,e,B,s,integral_left,x,level+1)
                call GG_Integration_y(ybar,ymax,func,e,B,s,integral_right,x,level+1)
            else
                !   don't have a value for level, so can't ask for level+1. I'll assume input level was in fact zero.
                call GG_Integration_y(ymin,ybar,func,e,B,s,integral_left,x,1)
                call GG_Integration_y(ybar,ymax,func,e,B,s,integral_right,x,1)
            end if
        end if


    !---    at this point we must have an accurate answer for integral_left and integral_right. 
    !       (or have run out of time)
    !       So return my best guess for int_ymin,ymax func(x,y) dy
        integral = integral_left + integral_right           !   updated guess

        return
    end subroutine GG_Integration_y

 

    real(kind=real64) function GaussianQuadrature( ymin,ymax,func,e,B,s,x )
!---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
!*      compute an approximation of the integral func(x,y) dy in the interval (ymin,ymax)
!*      using gaussian quadrature. Note that the GQ abscissae with +/- x_i have the same weight.
        real(kind=real64),intent(in)                ::      ymin,ymax   !   integration ranges
        real(kind=real64),external                  ::      func        !   function to integrate.
        integer,intent(in)                          ::      e           !   element index
        real(kind=real64),intent(in)                ::      B           !   Debye-Waller factor B
        real(kind=real64),intent(in)                ::      s           !   scattering s = q/4pi
        real(kind=real64),intent(in)                ::      x           !   fixed x position for integral over dy

        integer             ::      kk
        real(kind=real64)   ::      yy
        GaussianQuadrature = 0.0d0
        do kk = 1,size(GQ_abscissa)
            yy = (ymin+ymax)/2 + GQ_abscissa(kk)*(ymax-ymin)/2             
            GaussianQuadrature = GaussianQuadrature + func(e,B,s,x,yy) * GQ_weight(kk)             
        end do
        GaussianQuadrature = GaussianQuadrature * (ymax-ymin)/2
        return
    end function GaussianQuadrature

! --------------------------------------------------------------------------

    
end module Lib_inputTDS