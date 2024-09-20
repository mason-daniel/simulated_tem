c ******************************************************************************
c 
c                   FORTRAN SOURCE FILE       inputTDS.f
c
c *******************************************************************************
c                                                        by S. L. Dudarev, 
c                                                        Department of Materials,
c                      					                     University of Oxford
c							                                    Parks Road
c							                                    OXFORD OX1 3PH
c							                                    UK
c
c                        present address:                EURATOM/UKAEA Fusion Association
c                                                        Culham Science Centre
c                                                        Oxfordshire OX14 3DB
c                                                        UK
c
c application of this program has been described in the following publications:
c
c  S L Dudarev and M J Whelan, Surface Science 310 (1994) 373 
c  S L Dudarev, D D Vvedensky and M J Whelan, Phys. Review B50 (1994) 14525 
c  S L Dudarev, L M Peng and M J Whelan, Surface Science 330 (1995) 86
c  see also
c  L M Peng, S L Dudarev and M J Whelan, High-Energy Electron Diffraction and Microscopy
c (Oxford University Press, 2004) ISBN 0-19-850074-2
c ******************************************************************************** 
!   Edit by D.R. Mason 28/01/21
!       changed some obsolete numbered do loops
!       removed obsolete pause command
!       tweaked filenames to allow directory paths
!       no functional changes

	IMPLICIT REAL*8 (a-h,o-z)
	EXTERNAL fx,fxy,ylim,gauss1
	PARAMETER (pi=3.14159265359d0)
	character*2 elnames(50),elname
	common s,bdw,elname,elnames

	character*60 filename
	
	data elnames /'LI','BE','C ','O ','NA','MG','AL','SI','K ',
     ,'CA','TI','V ','CR','MN','FE','CO','NI','CU','ZN','GA','GE',
     ,'AS','RB','SR','Y ','ZR','NB','MO','RU','RH','PD','AG','CD',
     ,'SN','CS','BA','LA','CE','GD','HF','TA','W ','RE','OS','IR',
     ,'PT','AU','TL','PB','TH'/	
     

 	print *,"inputTDS.f"
 	print *,"^^^^^^^^^^"
 	print *,""     
     
 1	write(6,*) 
     *'element name, upper case, 2 characters, left justified ='
	read(5,*) elname

	iflag=0
	do i=1,50
    	if(elname.eq.elnames(i)) iflag=1
	end do
	if(iflag.eq.0) go to 1
	
	write(6,*) 'output filename (no more than 60 characters) ='
	read(5,fmt='(a)') filename
	
	write(6,*) 'kinetic energy (eV), E='
	read(5,*) E
	alamda=12.264306d0/dsqrt(E)/dsqrt(1.d0+0.97846707d-06*E)
	coeff=2.d0*alamda*(1.d0+E/0.511d06)

	write(6,*) 'mean square thermal displacement (A^2), U2='
	read(5,*) U2
	
	print *,"inputTDS.f info - element """,elname,""" KE ",E," U2 ",U2," filename """,trim(filename),""""

	
	BDW=8.d0*pi**2*u2               !   this is the Debye Waller factor B  Acta Cryst. (1996). A52, 456-470 eqn 3.
	
	OPEN(33,FILE=trim(filename),action="write")

	ae=1.d-6
	re=1.d-4
	a=-20.d0
	b=20.d0

	smin=0.d0
	smax=15.d0
	ds=0.1d0
	
	ns=(smax-smin)/ds+1

	do 3 is=1,ns
	
 		s=smin + ds*(is-1)          !   s = g/4pi = sin(theta)/lambda 
 		                            !   where theta is the scattering angle and lambda the electron wavelength

		ae=1.d-6
		re=1.d-4
		a=-20.d0
		b=20.d0

! two-dimensional integration using gaussian quadratures
		call gauss2(a,b,re,ae,sum,ylim,fx,fxy)

		fimag=sum*coeff*dexp(-0.5d0*s*s*bdw)
	
		write(33,99) is,s,fimag

 3	continue

 	close(35)

  99	format(1x,i3,1x,2(1x,e12.6))

  
 	print *,""
 	print *,"done"
 	print *,""
  
	stop
	end


c ***************************************************************************

	subroutine fx(x,w)
	implicit real*8 (a-h,o-z)
	w=1.d0
	return
	end

c ***************************************************************************

	subroutine fxy(x,y,w)
	implicit real*8 (a-h,o-z)
	character*2 elname,elnames(50)
	common s,bdw,elname,elnames

	x1=x-0.5d0*s
	x2=x1+s
	xq=x*x
	yq=y*y
	ar1=dsqrt(yq+x1*x1)
	ar2=dsqrt(yq+x2*x2)
	are=2.d0*bdw*(0.25d0*s*s-xq-yq)
	ex=dexp(are)

	if(elname.eq.elnames(1)) W=FLI(ar1)*FLI(ar2)*(1.d0-ex)
	if(elname.eq.elnames(2)) W=FBE(ar1)*FBE(ar2)*(1.d0-ex)
	if(elname.eq.elnames(3)) W=FC(ar1)*FC(ar2)*(1.d0-ex)
	if(elname.eq.elnames(4)) W=FO(ar1)*FO(ar2)*(1.d0-ex)
	if(elname.eq.elnames(5)) W=FNA(ar1)*FNA(ar2)*(1.d0-ex)
	if(elname.eq.elnames(6)) W=FMG(ar1)*FMG(ar2)*(1.d0-ex)
	if(elname.eq.elnames(7)) W=FAL(ar1)*FAL(ar2)*(1.d0-ex)
	if(elname.eq.elnames(8)) W=FSI(ar1)*FSI(ar2)*(1.d0-ex)
	if(elname.eq.elnames(9)) W=FK(ar1)*FK(ar2)*(1.d0-ex)
	if(elname.eq.elnames(10)) W=FCA(ar1)*FCA(ar2)*(1.d0-ex)
	if(elname.eq.elnames(11)) W=FTI(ar1)*FTI(ar2)*(1.d0-ex)
	if(elname.eq.elnames(12)) W=FV(ar1)*FV(ar2)*(1.d0-ex)
	if(elname.eq.elnames(13)) W=FCR(ar1)*FCR(ar2)*(1.d0-ex)
	if(elname.eq.elnames(14)) W=FMN(ar1)*FMN(ar2)*(1.d0-ex)
	if(elname.eq.elnames(15)) W=FFE(ar1)*FFE(ar2)*(1.d0-ex)
	if(elname.eq.elnames(16)) W=FCO(ar1)*FCO(ar2)*(1.d0-ex)
	if(elname.eq.elnames(17)) W=FNI(ar1)*FNI(ar2)*(1.d0-ex)
	if(elname.eq.elnames(18)) W=FCU(ar1)*FCU(ar2)*(1.d0-ex)
	if(elname.eq.elnames(19)) W=FZN(ar1)*FZN(ar2)*(1.d0-ex)
	if(elname.eq.elnames(20)) W=FGA(ar1)*FGA(ar2)*(1.d0-ex)
	if(elname.eq.elnames(21)) W=FGE(ar1)*FGE(ar2)*(1.d0-ex)
	if(elname.eq.elnames(22)) W=FAS(ar1)*FAS(ar2)*(1.d0-ex)
	if(elname.eq.elnames(23)) W=FRB(ar1)*FRB(ar2)*(1.d0-ex)
	if(elname.eq.elnames(24)) W=FSR(ar1)*FSR(ar2)*(1.d0-ex)
	if(elname.eq.elnames(25)) W=FY(ar1)*FY(ar2)*(1.d0-ex)
	if(elname.eq.elnames(26)) W=FZR(ar1)*FZR(ar2)*(1.d0-ex)
	if(elname.eq.elnames(27)) W=FNB(ar1)*FNB(ar2)*(1.d0-ex)
	if(elname.eq.elnames(28)) W=FMO(ar1)*FMO(ar2)*(1.d0-ex)
	if(elname.eq.elnames(29)) W=FRU(ar1)*FRU(ar2)*(1.d0-ex)
	if(elname.eq.elnames(30)) W=FRH(ar1)*FRH(ar2)*(1.d0-ex)	
	if(elname.eq.elnames(31)) W=FPD(ar1)*FPD(ar2)*(1.d0-ex)
	if(elname.eq.elnames(32)) W=FAG(ar1)*FAG(ar2)*(1.d0-ex)
	if(elname.eq.elnames(33)) W=FCD(ar1)*FCD(ar2)*(1.d0-ex)
	if(elname.eq.elnames(34)) W=FSN(ar1)*FSN(ar2)*(1.d0-ex)
	if(elname.eq.elnames(35)) W=FCS(ar1)*FCS(ar2)*(1.d0-ex)
	if(elname.eq.elnames(36)) W=FBA(ar1)*FBA(ar2)*(1.d0-ex)
	if(elname.eq.elnames(37)) W=FLA(ar1)*FLA(ar2)*(1.d0-ex)
	if(elname.eq.elnames(38)) W=FCE(ar1)*FCE(ar2)*(1.d0-ex)
	if(elname.eq.elnames(39)) W=FGD(ar1)*FGD(ar2)*(1.d0-ex)
	if(elname.eq.elnames(40)) W=FHF(ar1)*FHF(ar2)*(1.d0-ex)	
	if(elname.eq.elnames(41)) W=FTA(ar1)*FTA(ar2)*(1.d0-ex)
	if(elname.eq.elnames(42)) W=FW(ar1)*FW(ar2)*(1.d0-ex)
	if(elname.eq.elnames(43)) W=FRE(ar1)*FRE(ar2)*(1.d0-ex)
	if(elname.eq.elnames(44)) W=FOS(ar1)*FOS(ar2)*(1.d0-ex)
	if(elname.eq.elnames(45)) W=FIR(ar1)*FIR(ar2)*(1.d0-ex)
	if(elname.eq.elnames(46)) W=FPT(ar1)*FPT(ar2)*(1.d0-ex)
	if(elname.eq.elnames(47)) W=FAU(ar1)*FAU(ar2)*(1.d0-ex)
	if(elname.eq.elnames(48)) W=FTL(ar1)*FTL(ar2)*(1.d0-ex)
	if(elname.eq.elnames(49)) W=FPB(ar1)*FPB(ar2)*(1.d0-ex)	
	if(elname.eq.elnames(50)) W=FTH(ar1)*FTH(ar2)*(1.d0-ex)	

	return
	end

c ***************************************************************************

	subroutine ylim(x,ylim1,ylim2)
	implicit real*8 (a-h,o-z)
	ylim1=-dsqrt(400.d0-x*x)
	ylim2=-ylim1
	return
	end
	
	SUBROUTINE GAUSS2(A,B,RE,AE,SUM,YLIM,FX,FXY)
	IMPLICIT REAL*8(A-H,O-Z)
	EXTERNAL GAUSS1,YLIM,FX,FXY
	DIMENSION XX(8),WW(8),AA(80),BB(80),W(80)
	DATA XX/ 0.09501250983763744D0, 0.28160355077925891D0,
     ,		 0.45801677765722739D0, 0.61787624440264375D0,
     ,		 0.75540440835500303D0, 0.86563120238783174D0,
     ,		 0.94457502307323258D0, 0.98940093499164993D0/
	DATA WW/ 0.18945061045506850D0, 0.18260341504492359D0,
     ,		 0.16915651939500254D0, 0.14959598881657673D0,
     ,		 0.12462897125553387D0, 0.09515851168249278D0,
     ,		 0.06225352393864789D0, 0.02715245941175409D0/
	SUM=0.D0
	ER=RE*0.1D0
	EA=AE*0.01D0
	L=1
	A0=A
	B0=B
	AA(1)=A0
	BB(1)=B0
	ISUM=1
	GO TO 100
  10	W(1)=S
  25	A0=AA(L)
	B0=0.5D0*(BB(L)+A0)
	SS=W(L)
	ISUM=2
	GO TO 100
  20	A1=A0
	B1=B0
	S1=S
	B0=B0+B0-A0
	A0=B1
	ISUM=3
	GO TO 100
  30	SA=S+S1
	SB=DABS(SA-SS)
	IF(SB.LT.ER*DABS(SA)) GO TO 40
	IF(SB.LT.ER*DABS(SUM)) GO TO 40
	IF(SB.LT.EA) GO TO 40
	AA(L)=A1
	BB(L)=B1
	W(L)=S1
	L=L+1
	IF(L.GT.80) GO TO 50
	AA(L)=A0
	BB(L)=B0
	W(L)=S
	GO TO 25
  40	SUM=SUM+SA
	L=L-1
	IF(L.EQ.0) RETURN
	GO TO 25
  50	WRITE(6,*) ' number of iterations = 80 = limiting value'
	STOP
 100	ZZ=0.5D0*(B0-A0)
	ZM=ZZ+A0
	S=0.D0
	DO K=1,8
    	X=ZM+ZZ*XX(K)
    	CALL YLIM(X,YLIM1,YLIM2)
    	CALL FX(X,WW2)
    	CALL GAUSS1(YLIM1,YLIM2,RE,AE,W2,X,FXY)
    	X=ZM+ZM-X
    	CALL YLIM(X,YLIM1,YLIM2)
    	CALL FX(X,WW1)
    	CALL GAUSS1(YLIM1,YLIM2,RE,AE,W1,X,FXY)
       	S=S+WW(K)*(W1*WW1+W2*WW2)
	end do
	S=S*ZZ
	GO TO (10,20,30),ISUM
	END

c------------------------------------------------------------------------------------

	SUBROUTINE GAUSS1(A,B,RE,AE,SUM,X1,FXY)
	IMPLICIT REAL*8(A-H,O-Z)
	EXTERNAL FXY
	DIMENSION XX(8),WW(8),AA(20),BB(20),W(20)
	DATA XX/ 0.09501250983763744D0, 0.28160355077925891D0,
     ,		 0.45801677765722739D0, 0.61787624440264375D0,
     ,		 0.75540440835500303D0, 0.86563120238783174D0,
     ,		 0.94457502307323258D0, 0.98940093499164993D0/
	DATA WW/ 0.18945061045506850D0, 0.18260341504492359D0,
     ,		 0.16915651939500254D0, 0.14959598881657673D0,
     ,		 0.12462897125553387D0, 0.09515851168249278D0,
     ,		 0.06225352393864789D0, 0.02715245941175409D0/
	SUM=0.D0
	ER=RE*0.1D0
	EA=AE*0.01D0
	L=1
	A0=A
	B0=B
	AA(1)=A0
	BB(1)=B0
	ISUM=1
	GO TO 100
  10	W(1)=S
  25	A0=AA(L)
	B0=0.5D0*(BB(L)+A0)
	SS=W(L)
	ISUM=2
	GO TO 100
  20	A1=A0
	B1=B0
	S1=S
	B0=B0+B0-A0
	A0=B1
	ISUM=3
	GO TO 100
  30	SA=S+S1
	SB=DABS(SA-SS)
	IF(SB.LT.ER*DABS(SA)) GO TO 40
	IF(SB.LT.ER*DABS(SUM)) GO TO 40
	IF(SB.LT.EA) GO TO 40
	AA(L)=A1
	BB(L)=B1
	W(L)=S1
	L=L+1
	IF(L.GT.20) GO TO 50
	AA(L)=A0
	BB(L)=B0
	W(L)=S
	GO TO 25
  40	SUM=SUM+SA
	L=L-1
	IF(L.EQ.0) RETURN
	GO TO 25
  50	WRITE(6,*) ' number of iterations = 20 = limiting value'
	RETURN
 100	ZZ=0.5D0*(B0-A0)
	ZM=ZZ+A0
	S=0.D0
	DO K=1,8
    	X=ZM+ZZ*XX(K)
    	CALL FXY(X1,X,W2)
    	X=ZM+ZM-X
    	CALL FXY(X1,X,W1)
    	S=S+WW(K)*(W1+W2)
	end do
	S=S*ZZ
	GO TO (10,20,30),ISUM
	END

c --------------------------------------------------------------------------

	REAL*8 FUNCTION FAG(X)
C Program for evaluating the scattering amplitude for SILVER atom, inter-
C polation of data from International Tables for X-ray Crystallography, 1974
C argument X equals to q/4pi within the notations of Landau & Lifshits
	IMPLICIT REAL*8(A-H,O-Z)
	DIMENSION A(56),XX(56)
	PARAMETER (cf=41.78214d0)
	DATA Z/47.d0/
       DATA A/8.671D0,8.654D0,8.599D0,8.510D0,8.391D0,8.244D0,
     ,8.075D0,7.888D0,7.689D0,7.480D0,7.267D0,7.052D0,6.837D0,
     ,6.625D0,6.418D0,6.215D0,6.018D0,5.827D0,5.643D0,5.464D0,
     ,5.293D0,4.967D0,4.665D0,4.522D0,4.384D0,4.122D0,3.878D0,
     ,3.651D0,3.440D0,3.339D0,3.242D0,3.058D0,2.886D0,2.726D0,
     ,2.576D0,2.505D0,2.436D0,2.306D0,2.185D0,1.915D0,1.688D0,
     ,1.497d0,1.335D0,1.082D0,0.897D0,0.758D0,0.652D0,0.568D0,
     ,0.500D0,0.444D0,0.397D0,0.357D0,0.321D0,0.291D0,0.265D0,
     ,0.241D0/
C
       DATA XX/0.D0,0.01D0,0.02D0,0.03D0,0.04d0,0.05d0,0.06d0,
     ,0.07d0,0.08d0,0.09d0,0.10d0,0.11d0,0.12d0,0.13d0,0.14d0,
     ,0.15d0,0.16d0,0.17d0,0.18d0,0.19d0,0.20d0,0.22d0,0.24d0,
     ,0.25D0,0.26d0,0.28d0,0.30d0,0.32d0,0.34d0,0.35d0,0.36d0,
     ,0.38d0,0.40d0,0.42d0,0.44d0,0.45d0,0.46d0,0.48d0,0.50d0,
     ,0.55d0,0.60d0,0.65d0,0.70d0,0.80d0,0.90d0,1.00d0,1.10d0,
     ,1.20d0,1.30d0,1.40d0,1.50d0,1.60d0,1.70d0,1.80d0,1.90d0,
     ,2.00d0/
        X02=Z/(cf*A(56))-XX(56)**2
	IF(X.GT.XX(1)) GO TO 1
	FAG=A(1)
	RETURN
 1	IF(X.LT.XX(56)) GO TO 2
C for large arguments the screened Coulomb law is assumed, 
C F proportional to 1./(X**2+X02)
	FAG=(Z/cf)/(X**2+X02)
	RETURN
 2	IF(X.GT.XX(2)) GO TO 21
 	FAG=A(1)-(A(1)-A(2))*X*X/XX(2)**2
 	RETURN
 21	N=1
	GO TO 4
 3	N=N+1
 4	IF(X.LE.XX(N)) GO TO 5
	GO TO 3
 5	FAG=A(N-1)+(X-XX(N-1))*(A(N)-A(N-1))/(XX(N)-XX(N-1))
	RETURN
	END
	
	REAL*8 FUNCTION FAL(X)
C Program for evaluating the scattering amplitude for ALUMINIUM atom, inter-
C polation of data from International Tables for X-ray Crystallography, 1974
C argument X equals to q/4pi within the notations of Landau & Lifshits
	IMPLICIT REAL*8(A-H,O-Z)
	DIMENSION A(56),XX(56)
	PARAMETER (CF=41.78214D0)
	DATA Z/13.D0/
       DATA A/5.889D0,5.867D0,5.800D0,5.692D0,5.547D0,5.371D0,
     ,5.170D0,4.949D0,4.717D0,4.478D0,4.237D0,3.999D0,3.767D0,
     ,3.544D0,3.330D0,3.128D0,2.938D0,2.760D0,2.595D0,2.441D0,
     ,2.299D0,2.046D0,1.832D0,1.737D0,1.650D0,1.495D0,1.363D0,
     ,1.251D0,1.154D0,1.111D0,1.070D0,0.997D0,0.932D0,0.875D0,
     ,0.825D0,0.801D0,0.779D0,0.737D0,0.700D0,0.618D0,0.551D0,
     ,0.494d0,0.445D0,0.366D0,0.304D0,0.255D0,0.217D0,0.185D0,
     ,0.160D0,0.139D0,0.123D0,0.109D0,0.096D0,0.087D0,0.078D0,
     ,0.070D0/
C
       DATA XX/0.D0,0.01D0,0.02D0,0.03D0,0.04d0,0.05d0,0.06d0,
     ,0.07d0,0.08d0,0.09d0,0.10d0,0.11d0,0.12d0,0.13d0,0.14d0,
     ,0.15d0,0.16d0,0.17d0,0.18d0,0.19d0,0.20d0,0.22d0,0.24d0,
     ,0.25D0,0.26d0,0.28d0,0.30d0,0.32d0,0.34d0,0.35d0,0.36d0,
     ,0.38d0,0.40d0,0.42d0,0.44d0,0.45d0,0.46d0,0.48d0,0.50d0,
     ,0.55d0,0.60d0,0.65d0,0.70d0,0.80d0,0.90d0,1.00d0,1.10d0,
     ,1.20d0,1.30d0,1.40d0,1.50d0,1.60d0,1.70d0,1.80d0,1.90d0,
     ,2.00d0/
         X02=Z/(CF*A(56))-XX(56)**2
	IF(X.GT.XX(1)) GO TO 1
	FAL=A(1)
	RETURN
 1	IF(X.LT.XX(56)) GO TO 2
C for large arguments the screened Coulomb law is assumed
	FAL=(Z/CF)/(X**2+X02)
	RETURN
 2	IF(X.GT.XX(2)) GO TO 21
 	FAL=A(1)-(A(1)-A(2))*X*X/XX(2)**2
 	RETURN
 21	N=1
	GO TO 4
 3	N=N+1
 4	IF(X.LE.XX(N)) GO TO 5
	GO TO 3
 5	FAL=A(N-1)+(X-XX(N-1))*(A(N)-A(N-1))/(XX(N)-XX(N-1))
	RETURN
	END

	REAL*8 FUNCTION FAS(X)
C Program for evaluating the scattering amplitude for ARSENIC atom, inter-
C polation of data from International Tables for X-ray Crystallography, 1974
C argument X equals to q/4pi within the notations of Landau & Lifshits
	IMPLICIT REAL*8(A-H,O-Z)
	DIMENSION A(56),XX(56)
	PARAMETER (CF=41.78214D0)
	DATA Z/33.D0/
       DATA A/7.320D0,7.306D0,7.260D0,7.184D0,7.081D0,6.953D0,
     ,6.803D0,6.634D0,6.449D0,6.253D0,6.048D0,5.838D0,5.625D0,
     ,5.411D0,5.200D0,4.992D0,4.789D0,4.593D0,4.404D0,4.222D0,
     ,4.048D0,3.724D0,3.433D0,3.299D0,3.172D0,2.940D0,2.733D0,
     ,2.548D0,2.384D0,2.308D0,2.237D0,2.105D0,1.986D0,1.878D0,
     ,1.780D0,1.734D0,1.691D0,1.608D0,1.533D0,1.367D0,1.228D0,
     ,1.108d0,1.004D0,0.832D0,0.697D0,0.589D0,0.502D0,0.431D0,
     ,0.374D0,0.327D0,0.288D0,0.256D0,0.229D0,0.206D0,0.187D0,
     ,0.170D0/
C
       DATA XX/0.D0,0.01D0,0.02D0,0.03D0,0.04d0,0.05d0,0.06d0,
     ,0.07d0,0.08d0,0.09d0,0.10d0,0.11d0,0.12d0,0.13d0,0.14d0,
     ,0.15d0,0.16d0,0.17d0,0.18d0,0.19d0,0.20d0,0.22d0,0.24d0,
     ,0.25D0,0.26d0,0.28d0,0.30d0,0.32d0,0.34d0,0.35d0,0.36d0,
     ,0.38d0,0.40d0,0.42d0,0.44d0,0.45d0,0.46d0,0.48d0,0.50d0,
     ,0.55d0,0.60d0,0.65d0,0.70d0,0.80d0,0.90d0,1.00d0,1.10d0,
     ,1.20d0,1.30d0,1.40d0,1.50d0,1.60d0,1.70d0,1.80d0,1.90d0,
     ,2.00d0/
         X02=Z/(CF*A(56))-XX(56)**2
	IF(X.GT.XX(1)) GO TO 1
	FAS=A(1)
	RETURN
 1	IF(X.LT.XX(56)) GO TO 2
C for large arguments the screened Coulomb law is assumed
	FAS=(Z/CF)/(X**2+X02)
	RETURN
 2	IF(X.GT.XX(2)) GO TO 21
 	FAS=A(1)-(A(1)-A(2))*X*X/XX(2)**2
 	RETURN
 21	N=1
	GO TO 4
 3	N=N+1
 4	IF(X.LE.XX(N)) GO TO 5
	GO TO 3
 5	FAS=A(N-1)+(X-XX(N-1))*(A(N)-A(N-1))/(XX(N)-XX(N-1))
	RETURN
	END

	REAL*8 FUNCTION FAU(X)
C Program for evaluating the scattering amplitude for GOLD atom, inter-
C polation of data from International Tables for X-ray Crystallography, 1974
C argument X equals to q/4pi within the notations of Landau & Lifshits
	IMPLICIT REAL*8(A-H,O-Z)
	DIMENSION A(56),XX(56)
	PARAMETER (CF=41.78214D0)
	DATA Z/79.D0/
       DATA A/10.573D0,10.559D0,10.511D0,10.434D0,10.328D0,10.195D0,
     ,10.040D0,9.865D0,9.673D0,9.467D0,9.251D0,9.028D0,8.799D0,
     ,8.568D0,8.337D0,8.106D0,7.877D0,7.652D0,7.431D0,7.214D0,
     ,7.003D0,6.598D0,6.216D0,6.035D0,5.859D0,5.525D0,5.214D0,
     ,4.924D0,4.654D0,4.526D0,4.403D0,4.169D0,3.952D0,3.750D0,
     ,3.562D0,3.472D0,3.386D0,3.223D0,3.070D0,2.732D0,2.446D0,
     ,2.203d0,1.995D0,1.661D0,1.407D0,1.208D0,1.048D0,0.918D0,
     ,0.809D0,0.717D0,0.639D0,0.572D0,0.514D0,0.465D0,0.422D0,
     ,0.384D0/
C
       DATA XX/0.D0,0.01D0,0.02D0,0.03D0,0.04d0,0.05d0,0.06d0,
     ,0.07d0,0.08d0,0.09d0,0.10d0,0.11d0,0.12d0,0.13d0,0.14d0,
     ,0.15d0,0.16d0,0.17d0,0.18d0,0.19d0,0.20d0,0.22d0,0.24d0,
     ,0.25D0,0.26d0,0.28d0,0.30d0,0.32d0,0.34d0,0.35d0,0.36d0,
     ,0.38d0,0.40d0,0.42d0,0.44d0,0.45d0,0.46d0,0.48d0,0.50d0,
     ,0.55d0,0.60d0,0.65d0,0.70d0,0.80d0,0.90d0,1.00d0,1.10d0,
     ,1.20d0,1.30d0,1.40d0,1.50d0,1.60d0,1.70d0,1.80d0,1.90d0,
     ,2.00d0/
        X02=Z/(CF*A(56))-XX(56)**2
	IF(X.GT.XX(1)) GO TO 1
	FAU=A(1)
	RETURN
 1	IF(X.LT.XX(56)) GO TO 2
C for large arguments the screened Coulomb law is assumed
	FAU=(Z/cf)/(X**2+X02)
	RETURN
 2	IF(X.GT.XX(2)) GO TO 21
 	FAU=A(1)-(A(1)-A(2))*X*X/XX(2)**2
 	RETURN
 21	N=1
	GO TO 4
 3	N=N+1
 4	IF(X.LE.XX(N)) GO TO 5
	GO TO 3
 5	FAU=A(N-1)+(X-XX(N-1))*(A(N)-A(N-1))/(XX(N)-XX(N-1))
	RETURN
	END
c
	REAL*8 FUNCTION FBA(X)
C Program for evaluating the scattering amplitude for BARIUM atom, inter-
C polation of data from International Tables for X-ray Crystallography, 1974
C argument X equals to q/4pi within the notations of Landau & Lifshits
	IMPLICIT REAL*8(A-H,O-Z)
	DIMENSION A(56),XX(56)
	PARAMETER (CF=41.78214D0)
	DATA Z/56.D0/
       DATA A/18.267D0,18.157D0,17.828D0,17.309D0,16.636D0,15.854D0,
     ,15.008D0,14.138D0,13.278D0,12.451D0,11.675D0,10.958D0,10.302D0,
     ,9.707D0,9.168D0,8.682D0,8.241D0,7.840D0,7.474D0,7.139D0, 
     ,6.829D0,6.275D0,5.791D0,5.570D0,5.361D0,4.975D0,4.628D0,
     ,4.313D0,4.028D0,3.895D0,3.769D0,3.533D0,3.318D0,3.123D0,
     ,2.944D0,2.861D0,2.781D0,2.631D0,2.494D0,2.197D0,1.951d0,
     ,1.745D0,1.570D0,1.288D0,1.073D0,0.904D0,0.772D0,0.666D0,
     ,0.580D0,0.511D0,0.456D0,0.411D0,0.367D0,0.337D0,0.304D0,
     ,0.277D0/
C
       DATA XX/0.D0,0.01D0,0.02D0,0.03D0,0.04d0,0.05d0,0.06d0,
     ,0.07d0,0.08d0,0.09d0,0.10d0,0.11d0,0.12d0,0.13d0,0.14d0,
     ,0.15d0,0.16d0,0.17d0,0.18d0,0.19d0,0.20d0,0.22d0,0.24d0,
     ,0.25D0,0.26d0,0.28d0,0.30d0,0.32d0,0.34d0,0.35d0,0.36d0,
     ,0.38d0,0.40d0,0.42d0,0.44d0,0.45d0,0.46d0,0.48d0,0.50d0,
     ,0.55d0,0.60d0,0.65d0,0.70d0,0.80d0,0.90d0,1.00d0,1.10d0,
     ,1.20d0,1.30d0,1.40d0,1.50d0,1.60d0,1.70d0,1.80d0,1.90d0,
     ,2.00d0/
        X02=Z/(CF*A(56))-XX(56)**2
	IF(X.GT.XX(1)) GO TO 1
	FBA=A(1)
	RETURN
 1	IF(X.LT.XX(56)) GO TO 2
C for large arguments the screened Coulomb law is assumed
	FBA=(Z/CF)/(X**2+X02)
	RETURN
 2	IF(X.GT.XX(2)) GO TO 21
 	FBA=A(1)-(A(1)-A(2))*X*X/XX(2)**2
 	RETURN
 21	N=1
	GO TO 4
 3	N=N+1
 4	IF(X.LE.XX(N)) GO TO 5
	GO TO 3
 5	FBA=A(N-1)+(X-XX(N-1))*(A(N)-A(N-1))/(XX(N)-XX(N-1))
	RETURN
	END

	REAL*8 FUNCTION FBE(X)
C Program for evaluating the scattering amplitude for  a BERYLLIUM atom, inter-
C polation of data from International Tables for X-ray Crystallography, 1974
C argument X equals to q/4pi within the notations of Landau & Lifshits
	IMPLICIT REAL*8(A-H,O-Z)
	DIMENSION A(56),XX(56)
	PARAMETER (CF=41.78214D0)
	DATA Z/4.D0/
       DATA A/3.052D0,3.042D0,3.011D0,2.961D0,2.892D0,2.807D0,
     ,2.710D0,2.601D0,2.484D0,2.362D0,2.237D0,2.111D0,1.987D0,
     ,1.865D0,1.748D0,1.635D0,1.528D0,1.427D0,1.332D0,1.243D0,
     ,1.161D0,1.013D0,0.887D0,0.832D0,0.781D0,0.690D0,0.614D0,
     ,0.549D0,0.494D0,0.469D0,0.446D0,0.406D0,0.371D0,0.341D0,
     ,0.314D0,0.302D0,0.291D0,0.271D0,0.253D0,0.215D0,0.186D0,
     ,0.164d0,0.145D0,0.117D0,0.096D0,0.081D0,0.069D0,0.059D0,
     ,0.051D0,0.045D0,0.040D0,0.035D0,0.031D0,0.028D0,0.026D0,
     ,0.023D0/
C
       DATA XX/0.D0,0.01D0,0.02D0,0.03D0,0.04d0,0.05d0,0.06d0,
     ,0.07d0,0.08d0,0.09d0,0.10d0,0.11d0,0.12d0,0.13d0,0.14d0,
     ,0.15d0,0.16d0,0.17d0,0.18d0,0.19d0,0.20d0,0.22d0,0.24d0,
     ,0.25D0,0.26d0,0.28d0,0.30d0,0.32d0,0.34d0,0.35d0,0.36d0,
     ,0.38d0,0.40d0,0.42d0,0.44d0,0.45d0,0.46d0,0.48d0,0.50d0,
     ,0.55d0,0.60d0,0.65d0,0.70d0,0.80d0,0.90d0,1.00d0,1.10d0,
     ,1.20d0,1.30d0,1.40d0,1.50d0,1.60d0,1.70d0,1.80d0,1.90d0,
     ,2.00d0/
         X02=Z/(CF*A(56))-XX(56)**2
	IF(X.GT.XX(1)) GO TO 1
	FBE=A(1)
	RETURN
 1	IF(X.LT.XX(56)) GO TO 2
C for large arguments the screened Coulomb law is assumed
	FBE=(Z/CF)/(X**2+X02)
	RETURN
 2	IF(X.GT.XX(2)) GO TO 21
 	FBE=A(1)-(A(1)-A(2))*X*X/XX(2)**2
 	RETURN
 21	N=1
	GO TO 4
 3	N=N+1
 4	IF(X.LE.XX(N)) GO TO 5
	GO TO 3
 5	FBE=A(N-1)+(X-XX(N-1))*(A(N)-A(N-1))/(XX(N)-XX(N-1))
	RETURN
	END

	REAL*8 FUNCTION FC(X)
C Program for evaluating the scattering amplitude for  CARBON atom, inter-
C polation of data from International Tables for X-ray Crystallography, 1974
C argument X equals to q/4pi within the notations of Landau & Lifshits
	IMPLICIT REAL*8(A-H,O-Z)
	DIMENSION A(56),XX(56)
	PARAMETER (CF=41.78214D0)
	DATA Z/6.D0/
       DATA A/2.509D0,2.505D0,2.492D0,2.471D0,2.442D0,2.406D0,
     ,2.363D0,2.313D0,2.259D0,2.200D0,2.138D0,2.072D0,2.005D0,
     ,1.936D0,1.866D0,1.796D0,1.727D0,1.658D0,1.591D0,1.524D0,
     ,1.460D0,1.337D0,1.222D0,1.168D0,1.117D0,1.020D0,0.932D0,
     ,0.853D0,0.781D0,0.748D0,0.717D0,0.658D0,0.606D0,0.559D0,
     ,0.517D0,0.497D0,0.479D0,0.444D0,0.413D0,0.348D0,0.297D0,
     ,0.256d0,0.223D0,0.175D0,0.141D0,0.117D0,0.099D0,0.085D0,
     ,0.073D0,0.064D0,0.057D0,0.051D0,0.045D0,0.041D0,0.037D0,
     ,0.034D0/
C
       DATA XX/0.D0,0.01D0,0.02D0,0.03D0,0.04d0,0.05d0,0.06d0,
     ,0.07d0,0.08d0,0.09d0,0.10d0,0.11d0,0.12d0,0.13d0,0.14d0,
     ,0.15d0,0.16d0,0.17d0,0.18d0,0.19d0,0.20d0,0.22d0,0.24d0,
     ,0.25D0,0.26d0,0.28d0,0.30d0,0.32d0,0.34d0,0.35d0,0.36d0,
     ,0.38d0,0.40d0,0.42d0,0.44d0,0.45d0,0.46d0,0.48d0,0.50d0,
     ,0.55d0,0.60d0,0.65d0,0.70d0,0.80d0,0.90d0,1.00d0,1.10d0,
     ,1.20d0,1.30d0,1.40d0,1.50d0,1.60d0,1.70d0,1.80d0,1.90d0,
     ,2.00d0/
         X02=Z/(CF*A(56))-XX(56)**2
	IF(X.GT.XX(1)) GO TO 1
	FC=A(1)
	RETURN
 1	IF(X.LT.XX(56)) GO TO 2
C for large arguments the screened Coulomb law is assumed
	FC=(Z/CF)/(X**2+X02)
	RETURN
 2	IF(X.GT.XX(2)) GO TO 21
 	FC=A(1)-(A(1)-A(2))*X*X/XX(2)**2
 	RETURN
 21	N=1
	GO TO 4
 3	N=N+1
 4	IF(X.LE.XX(N)) GO TO 5
	GO TO 3
 5	FC=A(N-1)+(X-XX(N-1))*(A(N)-A(N-1))/(XX(N)-XX(N-1))
	RETURN
	END

	REAL*8 FUNCTION FCA(X)
C Program for evaluating the scattering amplitude for CALCIUM atom, inter-
C polation of data from International Tables for X-ray Crystallography, 1974
C argument X equals to q/4pi within the notations of Landau & Lifshits
	IMPLICIT REAL*8(A-H,O-Z)
	DIMENSION A(56),XX(56)
	PARAMETER (CF=41.78214D0)
	DATA Z/20.D0/
       DATA A/9.913D0,9.860D0,9.699D0,9.442D0,9.104D0,8.703D0,
     ,8.258D0,7.789D0,7.312D0,6.841D0,6.388D0,5.959D0,5.560D0,
     ,5.192D0,4.855D0,4.550D0,4.273D0,4.023D0,3.797D0,3.593D0,
     ,3.408D0,3.086D0,2.815D0,2.695D0,2.584D0,2.383D0,2.206D0,
     ,2.048D0,1.905D0,1.838D0,1.775D0,1.657D0,1.548D0,1.449D0,
     ,1.357D0,1.314D0,1.272D0,1.194D0,1.123D0,0.966D0,0.838D0,
     ,0.733d0,0.647D0,0.515D0,0.422D0,0.354D0,0.303D0,0.262D0,
     ,0.230D0,0.202D0,0.181D0,0.162D0,0.144D0,0.132D0,0.118D0,
     ,0.107D0/
C
       DATA XX/0.D0,0.01D0,0.02D0,0.03D0,0.04d0,0.05d0,0.06d0,
     ,0.07d0,0.08d0,0.09d0,0.10d0,0.11d0,0.12d0,0.13d0,0.14d0,
     ,0.15d0,0.16d0,0.17d0,0.18d0,0.19d0,0.20d0,0.22d0,0.24d0,
     ,0.25D0,0.26d0,0.28d0,0.30d0,0.32d0,0.34d0,0.35d0,0.36d0,
     ,0.38d0,0.40d0,0.42d0,0.44d0,0.45d0,0.46d0,0.48d0,0.50d0,
     ,0.55d0,0.60d0,0.65d0,0.70d0,0.80d0,0.90d0,1.00d0,1.10d0,
     ,1.20d0,1.30d0,1.40d0,1.50d0,1.60d0,1.70d0,1.80d0,1.90d0,
     ,2.00d0/
         X02=Z/(CF*A(56))-XX(56)**2
	IF(X.GT.XX(1)) GO TO 1
	FCA=A(1)
	RETURN
 1	IF(X.LT.XX(56)) GO TO 2
C for large arguments the screened Coulomb law is assumed
	FCA=(Z/CF)/(X**2+X02)
	RETURN
 2	IF(X.GT.XX(2)) GO TO 21
 	FCA=A(1)-(A(1)-A(2))*X*X/XX(2)**2
 	RETURN
 21	N=1
	GO TO 4
 3	N=N+1
 4	IF(X.LE.XX(N)) GO TO 5
	GO TO 3
 5	FCA=A(N-1)+(X-XX(N-1))*(A(N)-A(N-1))/(XX(N)-XX(N-1))
	RETURN
	END

	REAL*8 FUNCTION FCD(X)
C Program for evaluating the scattering amplitude for  a CADMIUM atom, inter-
C polation of data from International Tables for X-ray Crystallography, 1974
C argument X equals to q/4pi within the notations of Landau & Lifshits
	IMPLICIT REAL*8(A-H,O-Z)
	DIMENSION A(56),XX(56)
	PARAMETER (CF=41.78214D0)
	DATA Z/48.D0/
       DATA A/9.232D0,9.213D0,9.153D0,9.057D0,8.926D0,8.764D0,
     ,8.577D0,8.369D0,8.144D0,7.909D0,7.666D0,7.421D0,7.176D0,
     ,6.933D0,6.695D0,6.464D0,6.240D0,6.024D0,5.817D0,5.618D0,
     ,5.427D0,5.070D0,4.745D0,4.592D0,4.447D0,4.173D0,3.922D0,
     ,3.690D0,3.476D0,3.375D0,3.278D0,3.093D0,2.922D0,2.762D0,
     ,2.613D0,2.542D0,2.474D0,2.344D0,2.223D0,1.953D0,1.724D0,
     ,1.531d0,1.366D0,1.107D0,0.916D0,0.773D0,0.664D0,0.578D0,
     ,0.508D0,0.451D0,0.403D0,0.362D0,0.327D0,0.297D0,0.270D0,
     ,0.246D0/
C
       DATA XX/0.D0,0.01D0,0.02D0,0.03D0,0.04d0,0.05d0,0.06d0,
     ,0.07d0,0.08d0,0.09d0,0.10d0,0.11d0,0.12d0,0.13d0,0.14d0,
     ,0.15d0,0.16d0,0.17d0,0.18d0,0.19d0,0.20d0,0.22d0,0.24d0,
     ,0.25D0,0.26d0,0.28d0,0.30d0,0.32d0,0.34d0,0.35d0,0.36d0,
     ,0.38d0,0.40d0,0.42d0,0.44d0,0.45d0,0.46d0,0.48d0,0.50d0,
     ,0.55d0,0.60d0,0.65d0,0.70d0,0.80d0,0.90d0,1.00d0,1.10d0,
     ,1.20d0,1.30d0,1.40d0,1.50d0,1.60d0,1.70d0,1.80d0,1.90d0,
     ,2.00d0/
         X02=Z/(CF*A(56))-XX(56)**2
	IF(X.GT.XX(1)) GO TO 1
	FCD=A(1)
	RETURN
 1	IF(X.LT.XX(56)) GO TO 2
C for large arguments the screened Coulomb law is assumed
	FCD=(Z/CF)/(X**2+X02)
	RETURN
 2	IF(X.GT.XX(2)) GO TO 21
 	FCD=A(1)-(A(1)-A(2))*X*X/XX(2)**2
 	RETURN
 21	N=1
	GO TO 4
 3	N=N+1
 4	IF(X.LE.XX(N)) GO TO 5
	GO TO 3
 5	FCD=A(N-1)+(X-XX(N-1))*(A(N)-A(N-1))/(XX(N)-XX(N-1))
	RETURN
	END

	REAL*8 FUNCTION FCE(X)
C Program for evaluating the scattering amplitude for CERIUM atom, inter-
C polation of data from International Tables for X-ray Crystallography, 1974
C argument X equals to q/4pi within the notations of Landau & Lifshits
	IMPLICIT REAL*8(A-H,O-Z)
	DIMENSION A(53),XX(53)
	PARAMETER (CF=41.78214D0)
	DATA Z/58.D0/
       DATA A/17.378D0,16.10D0,15.46D0,14.77D0,
     ,14.03D0,13.29D0,12.56D0,11.85D0,11.19D0,10.561D0,9.981D0,
     ,9.448D0,8.958D0,8.507D0,8.094D0,7.714D0,7.365D0,7.041D0, 
     ,6.462D0,5.957D0,5.728D0,5.512D0,5.115D0,4.759D0,4.438D0,
     ,4.146D0,4.010D0,3.881D0,3.640D0,3.420D0,3.219D0,3.035D0,
     ,2.949D0,2.866D0,2.712D0,2.570D0,2.262D0,2.008D0,1.796d0,
     ,1.617D0,1.329D0,1.109D0,0.936D0,0.799D0,0.690D0,0.602D0,
     ,0.530D0,0.470D0,0.421D0,0.380D0,0.345D0,0.314D0,0.288D0/
C
       DATA XX/0.D0,0.04d0,0.05d0,0.06d0,
     ,0.07d0,0.08d0,0.09d0,0.10d0,0.11d0,0.12d0,0.13d0,0.14d0,
     ,0.15d0,0.16d0,0.17d0,0.18d0,0.19d0,0.20d0,0.22d0,0.24d0,
     ,0.25D0,0.26d0,0.28d0,0.30d0,0.32d0,0.34d0,0.35d0,0.36d0,
     ,0.38d0,0.40d0,0.42d0,0.44d0,0.45d0,0.46d0,0.48d0,0.50d0,
     ,0.55d0,0.60d0,0.65d0,0.70d0,0.80d0,0.90d0,1.00d0,1.10d0,
     ,1.20d0,1.30d0,1.40d0,1.50d0,1.60d0,1.70d0,1.80d0,1.90d0,
     ,2.00d0/
        X02=Z/(CF*A(53))-XX(53)**2
	IF(X.GT.XX(1)) GO TO 1
	FCE=A(1)
	RETURN
 1	IF(X.LT.XX(53)) GO TO 2

C for large arguments the screened Coulomb law is assumed

	FCE=(Z/CF)/(X**2+X02)
	RETURN
 2	IF(X.GT.XX(2)) GO TO 21
 	FCE=A(1)-(A(1)-A(2))*X*X/XX(2)**2
 	RETURN
 21	N=1
	GO TO 4
 3	N=N+1
 4	IF(X.LE.XX(N)) GO TO 5
	GO TO 3
 5	FCE=A(N-1)+(X-XX(N-1))*(A(N)-A(N-1))/(XX(N)-XX(N-1))
	RETURN
	END

	REAL*8 FUNCTION FCO(X)
C Program for evaluating the scattering amplitude for COBALT atom, inter-
C polation of data from International Tables for X-ray Crystallography, 1974
C argument X equals to q/4pi within the notations of Landau & Lifshits
	IMPLICIT REAL*8(A-H,O-Z)
	DIMENSION A(56),XX(56)
	PARAMETER (CF=41.78214D0)
	DATA Z/27.D0/
       DATA A/6.854D0,6.836D0,6.779D0,6.687D0,6.562D0,6.410D0,
     ,6.234D0,6.040D0,5.834D0,5.619D0,5.401D0,5.182D0,4.964D0,
     ,4.758D0,4.555D0,4.361D0,4.177D0,4.002D0,3.836D0,3.689D0,
     ,3.534D0,3.267D0,3.032D0,2.924D0,2.823D0,2.637D0,2.471D0,
     ,2.321D0,2.184D0,2.121D0,2.060D0,1.946D0,1.841D0,1.743D0,
     ,1.653D0,1.610D0,1.569D0,1.490D0,1.416D0,1.251D0,1.110D0,
     ,0.988d0,0.883D0,0.712D0,0.583D0,0.485D0,0.409D0,0.350D0,
     ,0.303D0,0.265D0,0.235D0,0.209D0,0.188D0,0.170D0,0.154D0,
     ,0.141D0/
C
       DATA XX/0.D0,0.01D0,0.02D0,0.03D0,0.04d0,0.05d0,0.06d0,
     ,0.07d0,0.08d0,0.09d0,0.10d0,0.11d0,0.12d0,0.13d0,0.14d0,
     ,0.15d0,0.16d0,0.17d0,0.18d0,0.19d0,0.20d0,0.22d0,0.24d0,
     ,0.25D0,0.26d0,0.28d0,0.30d0,0.32d0,0.34d0,0.35d0,0.36d0,
     ,0.38d0,0.40d0,0.42d0,0.44d0,0.45d0,0.46d0,0.48d0,0.50d0,
     ,0.55d0,0.60d0,0.65d0,0.70d0,0.80d0,0.90d0,1.00d0,1.10d0,
     ,1.20d0,1.30d0,1.40d0,1.50d0,1.60d0,1.70d0,1.80d0,1.90d0,
     ,2.00d0/
         X02=Z/(CF*A(56))-XX(56)**2
	IF(X.GT.XX(1)) GO TO 1
	FCO=A(1)
	RETURN
 1	IF(X.LT.XX(56)) GO TO 2
C for large arguments the screened Coulomb law is assumed
	FCO=(Z/CF)/(X**2+X02)
	RETURN
 2	IF(X.GT.XX(2)) GO TO 21
 	FCO=A(1)-(A(1)-A(2))*X*X/XX(2)**2
 	RETURN
 21	N=1
	GO TO 4
 3	N=N+1
 4	IF(X.LE.XX(N)) GO TO 5
	GO TO 3
 5	FCO=A(N-1)+(X-XX(N-1))*(A(N)-A(N-1))/(XX(N)-XX(N-1))
	RETURN
	END

	REAL*8 FUNCTION FCR(X)
C Program for evaluating the scattering amplitude for CHROMIUM atom, inter-
C polation of data from International Tables for X-ray Crystallography, 1974
C argument X equals to q/4pi within the notations of Landau & Lifshits
	IMPLICIT REAL*8(A-H,O-Z)
	DIMENSION A(56),XX(56)
	PARAMETER (CF=41.78214D0)
	DATA Z/24.D0/
       DATA A/6.969D0,6.945D0,6.875D0,6.762D0,6.610D0,6.427D0,
     ,6.221D0,5.997D0,5.764D0,5.527D0,5.291D0,5.061D0,4.838D0,
     ,4.625D0,4.423D0,4.231D0,4.051D0,3.882D0,3.723D0,3.574D0,
     ,3.434D0,3.179D0,2.953D0,2.849D0,2.750D0,2.568D0,2.403D0,
     ,2.252D0,2.114D0,2.049D0,1.987D0,1.870D0,1.761D0,1.660D0,
     ,1.567D0,1.523D0,1.480D0,1.399D0,1.323D0,1.155D0,1.014D0,
     ,0.894d0,0.792D0,0.631D0,0.514D0,0.427D0,0.361D0,0.310D0,
     ,0.269D0,0.237D0,0.210D0,0.188D0,0.169D0,0.154D0,0.139D0,
     ,0.123D0/
C
       DATA XX/0.D0,0.01D0,0.02D0,0.03D0,0.04d0,0.05d0,0.06d0,
     ,0.07d0,0.08d0,0.09d0,0.10d0,0.11d0,0.12d0,0.13d0,0.14d0,
     ,0.15d0,0.16d0,0.17d0,0.18d0,0.19d0,0.20d0,0.22d0,0.24d0,
     ,0.25D0,0.26d0,0.28d0,0.30d0,0.32d0,0.34d0,0.35d0,0.36d0,
     ,0.38d0,0.40d0,0.42d0,0.44d0,0.45d0,0.46d0,0.48d0,0.50d0,
     ,0.55d0,0.60d0,0.65d0,0.70d0,0.80d0,0.90d0,1.00d0,1.10d0,
     ,1.20d0,1.30d0,1.40d0,1.50d0,1.60d0,1.70d0,1.80d0,1.90d0,
     ,2.00d0/
         X02=Z/(CF*A(56))-XX(56)**2
	IF(X.GT.XX(1)) GO TO 1
	FCR=A(1)
	RETURN
 1	IF(X.LT.XX(56)) GO TO 2
C for large arguments the screened Coulomb law is assumed
	FCR=(Z/CF)/(X**2+X02)
	RETURN
 2	IF(X.GT.XX(2)) GO TO 21
 	FCR=A(1)-(A(1)-A(2))*X*X/XX(2)**2
 	RETURN
 21	N=1
	GO TO 4
 3	N=N+1
 4	IF(X.LE.XX(N)) GO TO 5
	GO TO 3
 5	FCR=A(N-1)+(X-XX(N-1))*(A(N)-A(N-1))/(XX(N)-XX(N-1))
	RETURN
	END

	REAL*8 FUNCTION FCS(X)
C Program for evaluating the scattering amplitude for CESIUM atom, inter-
C polation of data from International Tables for X-ray Crystallography, 1974
C argument X equals to q/4pi within the notations of Landau & Lifshits
	IMPLICIT REAL*8(A-H,O-Z)
	DIMENSION A(56),XX(56)
	PARAMETER (CF=41.78214D0)
	DATA Z/55.D0/
       DATA A/16.508D0,16.391D0,16.050D0,15.521D0,14.855D0,14.106D0,
     ,13.326D0,12.556D0,11.823D0,11.145D0,10.525D0,9.965D0,9.458D0,
     ,9.000D0,8.583D0,8.201D0,7.848D0,7.519D0,7.212D0,6.922D0, 
     ,6.649D0,6.143D0,5.684D0,5.471D0,5.268D0,4.890D0,4.547D0,
     ,4.235D0,3.953D0,3.822D0,3.697D0,3.465D0,3.255D0,3.064D0,
     ,2.890D0,2.809D0,2.731D0,2.586D0,2.453D0,2.163D0,1.923d0,
     ,1.721D0,1.548D0,1.269D0,1.055D0,0.888D0,0.758D0,0.654D0,
     ,0.570D0,0.502D0,0.450D0,0.405D0,0.363D0,0.332D0,0.299D0,
     ,0.272D0/
C
       DATA XX/0.D0,0.01D0,0.02D0,0.03D0,0.04d0,0.05d0,0.06d0,
     ,0.07d0,0.08d0,0.09d0,0.10d0,0.11d0,0.12d0,0.13d0,0.14d0,
     ,0.15d0,0.16d0,0.17d0,0.18d0,0.19d0,0.20d0,0.22d0,0.24d0,
     ,0.25D0,0.26d0,0.28d0,0.30d0,0.32d0,0.34d0,0.35d0,0.36d0,
     ,0.38d0,0.40d0,0.42d0,0.44d0,0.45d0,0.46d0,0.48d0,0.50d0,
     ,0.55d0,0.60d0,0.65d0,0.70d0,0.80d0,0.90d0,1.00d0,1.10d0,
     ,1.20d0,1.30d0,1.40d0,1.50d0,1.60d0,1.70d0,1.80d0,1.90d0,
     ,2.00d0/
        X02=Z/(CF*A(56))-XX(56)**2
	IF(X.GT.XX(1)) GO TO 1
	FCS=A(1)
	RETURN
 1	IF(X.LT.XX(56)) GO TO 2
C for large arguments the screened Coulomb law is assumed
	FCS=(Z/CF)/(X**2+X02)
	RETURN
 2	IF(X.GT.XX(2)) GO TO 21
 	FCS=A(1)-(A(1)-A(2))*X*X/XX(2)**2
 	RETURN
 21	N=1
	GO TO 4
 3	N=N+1
 4	IF(X.LE.XX(N)) GO TO 5
	GO TO 3
 5	FCS=A(N-1)+(X-XX(N-1))*(A(N)-A(N-1))/(XX(N)-XX(N-1))
	RETURN
	END

	REAL*8 FUNCTION FCU(X)
C Program for evaluating the scattering amplitude for a COPPER atom, inter-
C polation of data from International Tables for X-ray Crystallography, 1974
C argument X equals to q/4pi within the notations of Landau & Lifshits
	IMPLICIT REAL*8(A-H,O-Z)
	DIMENSION A(56),XX(56)
	PARAMETER (CF=41.78214D0)
	DATA Z/29.D0/
       DATA A/5.600D0,5.587D0,5.547D0,5.482D0,5.395D0,5.287D0,
     ,5.165D0,5.029D0,4.886D0,4.737D0,4.585D0,4.434D0,4.285D0,
     ,4.139D0,3.998D0,3.862D0,3.731D0,3.607D0,3.488D0,3.375D0,
     ,3.276D0,3.067D0,2.885D0,2.800D0,2.719D0,2.568D0,2.428D0,
     ,2.299D0,2.180D0,2.123D0,2.069D0,1.965D0,1.868D0,1.777D0,
     ,1.691D0,1.651D0,1.611D0,1.535D0,1.464D0,1.303D0,1.163D0,
     ,1.041d0,0.935D0,0.761D0,0.626D0,0.523D0,0.442D0,0.378D0,
     ,0.327D0,0.285D0,0.252D0,0.224D0,0.201D0,0.182D0,0.165D0,
     ,0.150D0/
C
       DATA XX/0.D0,0.01D0,0.02D0,0.03D0,0.04d0,0.05d0,0.06d0,
     ,0.07d0,0.08d0,0.09d0,0.10d0,0.11d0,0.12d0,0.13d0,0.14d0,
     ,0.15d0,0.16d0,0.17d0,0.18d0,0.19d0,0.20d0,0.22d0,0.24d0,
     ,0.25D0,0.26d0,0.28d0,0.30d0,0.32d0,0.34d0,0.35d0,0.36d0,
     ,0.38d0,0.40d0,0.42d0,0.44d0,0.45d0,0.46d0,0.48d0,0.50d0,
     ,0.55d0,0.60d0,0.65d0,0.70d0,0.80d0,0.90d0,1.00d0,1.10d0,
     ,1.20d0,1.30d0,1.40d0,1.50d0,1.60d0,1.70d0,1.80d0,1.90d0,
     ,2.00d0/
         X02=Z/(CF*A(56))-XX(56)**2
	IF(X.GT.XX(1)) GO TO 1
	FCU=A(1)
	RETURN
 1	IF(X.LT.XX(56)) GO TO 2
C for large arguments the screened Coulomb law is assumed
	FCU=(Z/CF)/(X**2+X02)
	RETURN
 2	IF(X.GT.XX(2)) GO TO 21
 	FCU=A(1)-(A(1)-A(2))*X*X/XX(2)**2
 	RETURN
 21	N=1
	GO TO 4
 3	N=N+1
 4	IF(X.LE.XX(N)) GO TO 5
	GO TO 3
 5	FCU=A(N-1)+(X-XX(N-1))*(A(N)-A(N-1))/(XX(N)-XX(N-1))
	RETURN
	END

	REAL*8 FUNCTION FFE(X)
C Program for evaluating the scattering amplitude for IRON atom, inter-
C polation of data from International Tables for X-ray Crystallography, 1974
C argument X equals to q/4pi within the notations of Landau & Lifshits
	IMPLICIT REAL*8(A-H,O-Z)
	DIMENSION A(56),XX(56)
	PARAMETER (CF=41.78214D0)
	DATA Z/26.D0/
       DATA A/7.165D0,7.145D0,7.081D0,6.978D0,6.839D0,6.669D0,
     ,6.474D0,6.260D0,6.032D0,5.796D0,5.558D0,5.320D0,5.087D0,
     ,4.861D0,4.644D0,4.436D0,4.240D0,4.054D0,3.880D0,3.716D0,
     ,3.562D0,3.284D0,3.039D0,2.928D0,2.824D0,2.632D0,2.461D0,
     ,2.308D0,2.168D0,2.104D0,2.042D0,1.925D0,1.818D0,1.719D0,
     ,1.628D0,1.584D0,1.542D0,1.462D0,1.388D0,1.222D0,1.080D0,
     ,0.959d0,0.854D0,0.686D0,0.561D0,0.466D0,0.393D0,0.336D0,
     ,0.291D0,0.255D0,0.226D0,0.202D0,0.181D0,0.165D0,0.149D0,
     ,0.136D0/
C
       DATA XX/0.D0,0.01D0,0.02D0,0.03D0,0.04d0,0.05d0,0.06d0,
     ,0.07d0,0.08d0,0.09d0,0.10d0,0.11d0,0.12d0,0.13d0,0.14d0,
     ,0.15d0,0.16d0,0.17d0,0.18d0,0.19d0,0.20d0,0.22d0,0.24d0,
     ,0.25D0,0.26d0,0.28d0,0.30d0,0.32d0,0.34d0,0.35d0,0.36d0,
     ,0.38d0,0.40d0,0.42d0,0.44d0,0.45d0,0.46d0,0.48d0,0.50d0,
     ,0.55d0,0.60d0,0.65d0,0.70d0,0.80d0,0.90d0,1.00d0,1.10d0,
     ,1.20d0,1.30d0,1.40d0,1.50d0,1.60d0,1.70d0,1.80d0,1.90d0,
     ,2.00d0/
         X02=Z/(CF*A(56))-XX(56)**2
	IF(X.GT.XX(1)) GO TO 1
	FFE=A(1)
	RETURN
 1	IF(X.LT.XX(56)) GO TO 2
C for large arguments the screened Coulomb law is assumed
	FFE=(Z/CF)/(X**2+X02)
	RETURN
 2	IF(X.GT.XX(2)) GO TO 21
 	FFE=A(1)-(A(1)-A(2))*X*X/XX(2)**2
 	RETURN
 21	N=1
	GO TO 4
 3	N=N+1
 4	IF(X.LE.XX(N)) GO TO 5
	GO TO 3
 5	FFE=A(N-1)+(X-XX(N-1))*(A(N)-A(N-1))/(XX(N)-XX(N-1))
	RETURN
	END

	REAL*8 FUNCTION FGA(X)
C Program for evaluating the scattering amplitude for GALLIUM atom, inter-
C polation of data from International Tables for X-ray Crystallography, 1974
C argument X equals to q/4pi within the notations of Landau & Lifshits
	IMPLICIT REAL*8(A-H,O-Z)
	DIMENSION A(56),XX(56)
	PARAMETER (CF=41.78214D0)
	DATA Z/31.D0/
       DATA A/7.108D0,7.088D0,7.027D0,6.927D0,6.792D0,6.629D0,
     ,6.441D0,6.236D0,6.017D0,5.792D0,5.564D0,5.337D0,5.113D0,
     ,4.896D0,4.686D0,4.486D0,4.295D0,4.114D0,3.942D0,3.781D0,
     ,3.629D0,3.352D0,3.108D0,2.997D0,2.892D0,2.701D0,2.531D0,
     ,2.379D0,2.242D0,2.179D0,2.119D0,2.006D0,1.903D0,1.808D0,
     ,1.720D0,1.679D0,1.639D0,1.563D0,1.492D0,1.334D0,1.197D0,
     ,1.078d0,0.973D0,0.800D0,0.665D0,0.558D0,0.473D0,0.405D0,
     ,0.350D0,0.306D0,0.270D0,0.240D0,0.215D0,0.194D0,0.175D0,
     ,0.160D0/
C
       DATA XX/0.D0,0.01D0,0.02D0,0.03D0,0.04d0,0.05d0,0.06d0,
     ,0.07d0,0.08d0,0.09d0,0.10d0,0.11d0,0.12d0,0.13d0,0.14d0,
     ,0.15d0,0.16d0,0.17d0,0.18d0,0.19d0,0.20d0,0.22d0,0.24d0,
     ,0.25D0,0.26d0,0.28d0,0.30d0,0.32d0,0.34d0,0.35d0,0.36d0,
     ,0.38d0,0.40d0,0.42d0,0.44d0,0.45d0,0.46d0,0.48d0,0.50d0,
     ,0.55d0,0.60d0,0.65d0,0.70d0,0.80d0,0.90d0,1.00d0,1.10d0,
     ,1.20d0,1.30d0,1.40d0,1.50d0,1.60d0,1.70d0,1.80d0,1.90d0,
     ,2.00d0/
         X02=Z/(CF*A(56))-XX(56)**2
	IF(X.GT.XX(1)) GO TO 1
	FGA=A(1)
	RETURN
 1	IF(X.LT.XX(56)) GO TO 2
C for large arguments the screened Coulomb law is assumed
	FGA=(Z/CF)/(X**2+X02)
	RETURN
 2	IF(X.GT.XX(2)) GO TO 21
 	FGA=A(1)-(A(1)-A(2))*X*X/XX(2)**2
 	RETURN
 21	N=1
	GO TO 4
 3	N=N+1
 4	IF(X.LE.XX(N)) GO TO 5
	GO TO 3
 5	FGA=A(N-1)+(X-XX(N-1))*(A(N)-A(N-1))/(XX(N)-XX(N-1))
	RETURN
	END

	REAL*8 FUNCTION FGD(X)
C Program for evaluating the scattering amplitude for a GADOLINIUM atom, inter-
C polation of data from International Tables for X-ray Crystallography, 1974
C argument X equals to q/4pi within the notations of Landau & Lifshits
	IMPLICIT REAL*8(A-H,O-Z)
	DIMENSION A(53),XX(53)
	PARAMETER (CF=41.78214D0)
	DATA Z/64.D0/
       DATA A/15.266D0,14.30D0,13.81D0,13.27D0,
     ,12.70D0,12.11D0,11.52D0,10.95D0,10.39D0,9.871D0,9.382D0,
     ,8.926D0,8.505D0,8.114D0,7.754D0,7.422D0,7.114D0,6.828D0, 
     ,6.316D0,5.868D0,5.664D0,5.472D0,5.117D0,4.796D0,4.504D0,
     ,4.238D0,4.113D0,3.993D0,3.767D0,3.559D0,3.367D0,3.189D0,
     ,3.105D0,3.025D0,2.872D0,2.730D0,2.419D0,2.158D0,1.937d0,
     ,1.749D0,1.446D0,1.213D0,1.030D0,0.883D0,0.763D0,0.666D0,
     ,0.585D0,0.519D0,0.463D0,0.416D0,0.377D0,0.343D0,0.313D0/
C
       DATA XX/0.D0,0.04d0,0.05d0,0.06d0,
     ,0.07d0,0.08d0,0.09d0,0.10d0,0.11d0,0.12d0,0.13d0,0.14d0,
     ,0.15d0,0.16d0,0.17d0,0.18d0,0.19d0,0.20d0,0.22d0,0.24d0,
     ,0.25D0,0.26d0,0.28d0,0.30d0,0.32d0,0.34d0,0.35d0,0.36d0,
     ,0.38d0,0.40d0,0.42d0,0.44d0,0.45d0,0.46d0,0.48d0,0.50d0,
     ,0.55d0,0.60d0,0.65d0,0.70d0,0.80d0,0.90d0,1.00d0,1.10d0,
     ,1.20d0,1.30d0,1.40d0,1.50d0,1.60d0,1.70d0,1.80d0,1.90d0,
     ,2.00d0/
        X02=Z/(CF*A(53))-XX(53)**2
	IF(X.GT.XX(1)) GO TO 1
	FGD=A(1)
	RETURN
 1	IF(X.LT.XX(53)) GO TO 2

C for large arguments the screened Coulomb law is assumed

	FGD=(Z/CF)/(X**2+X02)
	RETURN
 2	IF(X.GT.XX(2)) GO TO 21
 	FGD=A(1)-(A(1)-A(2))*X*X/XX(2)**2
 	RETURN
 21	N=1
	GO TO 4
 3	N=N+1
 4	IF(X.LE.XX(N)) GO TO 5
	GO TO 3
 5	FGD=A(N-1)+(X-XX(N-1))*(A(N)-A(N-1))/(XX(N)-XX(N-1))
	RETURN
	END
	REAL*8 FUNCTION FGE(X)
C Program for evaluating the scattering amplitude for  SILICON atom, inter-
C polation of data from International Tables for X-ray Crystallography, 1974
C argument X equals to q/4pi within the notations of Landau & Lifshits
	IMPLICIT REAL*8(A-H,O-Z)
	DIMENSION A(56),XX(56)
	PARAMETER (CF=41.78214D0)
	DATA Z/32.D0/
       DATA A/7.378D0,7.359D0,7.303D0,7.211D0,7.088D0,6.935D0,
     ,6.759D0,6.562D0,6.351D0,6.129D0,5.902D0,5.672D0,5.442D0,
     ,5.217D0,4.996D0,4.783D0,4.578D0,4.382D0,4.195D0,4.017D0,
     ,3.849D0,3.541D0,3.268D0,3.143D0,3.026D0,2.813D0,2.623D0,
     ,2.455D0,2.304D0,2.235D0,2.169D0,2.048D0,1.938D0,1.837D0,
     ,1.745D0,1.702D0,1.661D0,1.583D0,1.510D0,1.349D0,1.212D0,
     ,1.093d0,0.989D0,0.817D0,0.681D0,0.574D0,0.488D0,0.418D0,
     ,0.362D0,0.317D0,0.279D0,0.248D0,0.222D0,0.200D0,0.181D0,
     ,0.165D0/
C
       DATA XX/0.D0,0.01D0,0.02D0,0.03D0,0.04d0,0.05d0,0.06d0,
     ,0.07d0,0.08d0,0.09d0,0.10d0,0.11d0,0.12d0,0.13d0,0.14d0,
     ,0.15d0,0.16d0,0.17d0,0.18d0,0.19d0,0.20d0,0.22d0,0.24d0,
     ,0.25D0,0.26d0,0.28d0,0.30d0,0.32d0,0.34d0,0.35d0,0.36d0,
     ,0.38d0,0.40d0,0.42d0,0.44d0,0.45d0,0.46d0,0.48d0,0.50d0,
     ,0.55d0,0.60d0,0.65d0,0.70d0,0.80d0,0.90d0,1.00d0,1.10d0,
     ,1.20d0,1.30d0,1.40d0,1.50d0,1.60d0,1.70d0,1.80d0,1.90d0,
     ,2.00d0/
         X02=Z/(CF*A(56))-XX(56)**2
	IF(X.GT.XX(1)) GO TO 1
	FGE=A(1)
	RETURN
 1	IF(X.LT.XX(56)) GO TO 2
C for large arguments the screened Coulomb law is assumed
	FGE=(Z/CF)/(X**2+X02)
	RETURN
 2	IF(X.GT.XX(2)) GO TO 21
 	FGE=A(1)-(A(1)-A(2))*X*X/XX(2)**2
 	RETURN
 21	N=1
	GO TO 4
 3	N=N+1
 4	IF(X.LE.XX(N)) GO TO 5
	GO TO 3
 5	FGE=A(N-1)+(X-XX(N-1))*(A(N)-A(N-1))/(XX(N)-XX(N-1))
	RETURN
	END
c **************************************************************************

	REAL*8 FUNCTION FHF(S)
C Program for evaluation of scattering amplitude for a HAFNIUM atom, inter-
C polation of data from International Tables for X-ray Crystallography, 1974
C argument S equals to q/(4*pi) within the notation of Landau & Lifshits
	IMPLICIT REAL*8(A-H,O-Z)
	DIMENSION A(53),SS(53)
	PARAMETER (CF=41.78214D0)
	DATA Z/72.D0/
       DATA A/13.177D0,12.55D0,12.23D0,
     ,11.85D0,11.45D0,11.02D0,10.59D0,10.16D0,9.73D0,9.308D0,
     ,8.907D0,8.525D0,8.163D0,7.822D0,7.502D0,7.202D0,6.922D0,
     ,6.660D0,6.185D0,5.768D0,5.578D0,5.399D0,5.069D0,4.772D0,
     ,4.503D0,4.258D0,4.143D0,4.033D0,3.825D0,3.632D0,3.454D0,
     ,3.288D0,3.209D0,3.133D0,2.988D0,2.853D0,2.551D0,2.294D0,
     ,2.073D0,1.882D0,1.571D0,1.329D0,1.138D0,0.982D0,0.854D0,
     ,0.748D0,0.660D0,0.585D0,0.522D0,0.469D0,0.423D0,0.384D0,
     ,0.350D0/
C
       DATA SS/0.D0,0.04d0,0.05d0,0.06d0,
     ,0.07d0,0.08d0,0.09d0,0.10d0,0.11d0,0.12d0,0.13d0,0.14d0,
     ,0.15d0,0.16d0,0.17d0,0.18d0,0.19d0,0.20d0,0.22d0,0.24d0,
     ,0.25D0,0.26d0,0.28d0,0.30d0,0.32d0,0.34d0,0.35d0,0.36d0,
     ,0.38d0,0.40d0,0.42d0,0.44d0,0.45d0,0.46d0,0.48d0,0.50d0,
     ,0.55d0,0.60d0,0.65d0,0.70d0,0.80d0,0.90d0,1.00d0,1.10d0,
     ,1.20d0,1.30d0,1.40d0,1.50d0,1.60d0,1.70d0,1.80d0,1.90d0,
     ,2.00d0/
        S02=Z/(CF*A(53))-SS(53)**2
	IF(S.GT.SS(1)) GO TO 1
	FHF=A(1)
	RETURN
 1	IF(S.LT.SS(53)) GO TO 2
C for large arguments the screened Coulomb law is assumed
	FHF=(Z/CF)/(S**2+S02)
	RETURN
 2	IF(S.GT.SS(2)) GO TO 21	
        FHF=A(1)-(A(1)-A(2))*S*S/SS(2)**2
        RETURN
 21     N=1
	GO TO 4
 3	N=N+1
 4	IF(S.LE.SS(N)) GO TO 5
	GO TO 3
 5	FHF=A(N-1)+(S-SS(N-1))*(A(N)-A(N-1))/(SS(N)-SS(N-1))
	RETURN
	END
	


c **************************************************************************

	REAL*8 FUNCTION FIR(S)
C Program for evaluation of scattering amplitude for an IRIDIUM atom, inter-
C polation of data from International Tables for X-ray Crystallography, 1974
C argument S equals to q/(4*pi) within the notation of Landau & Lifshits
	IMPLICIT REAL*8(A-H,O-Z)
	DIMENSION A(53),SS(53)
	PARAMETER (CF=41.78214D0)
	DATA Z/77.D0/
       DATA A/11.718D0,11.37D0,11.18D0,
     ,10.96D0,10.72D0,10.45D0,10.17D0,9.88D0,9.58D0,9.281D0,
     ,8.982D0,8.686D0,8.395D0,8.111D0,7.836D0,7.570D0,
     ,7.313D0,7.067D0,6.604D0,6.180D0,5.982D0,5.792D0,5.437D0,
     ,5.113D0,4.816D0,4.543D0,4.415D0,4.293D0,4.061D0,3.848D0,
     ,3.651D0,3.468D0,3.382D0,3.299D0,3.141D0,2.994D0,2.669D0,
     ,2.394D0,2.160D0,1.959D0,1.634D0,1.385D0,1.189D0,1.031D0,
     ,0.901D0,0.793D0,0.702D0,0.624D0,0.558D0,0.502D0,0.453D0,
     ,0.411D0,0.374d0/
C
       DATA SS/0.D0,0.04d0,0.05d0,0.06d0,
     ,0.07d0,0.08d0,0.09d0,0.10d0,0.11d0,0.12d0,0.13d0,0.14d0,
     ,0.15d0,0.16d0,0.17d0,0.18d0,0.19d0,0.20d0,0.22d0,0.24d0,
     ,0.25D0,0.26d0,0.28d0,0.30d0,0.32d0,0.34d0,0.35d0,0.36d0,
     ,0.38d0,0.40d0,0.42d0,0.44d0,0.45d0,0.46d0,0.48d0,0.50d0,
     ,0.55d0,0.60d0,0.65d0,0.70d0,0.80d0,0.90d0,1.00d0,1.10d0,
     ,1.20d0,1.30d0,1.40d0,1.50d0,1.60d0,1.70d0,1.80d0,1.90d0,
     ,2.00d0/
        S02=Z/(CF*A(53))-SS(53)**2
	IF(S.GT.SS(1)) GO TO 1
	FIR=A(1)
	RETURN
 1	IF(S.LT.SS(53)) GO TO 2
C for large arguments the screened Coulomb law is assumed
	FIR=(Z/CF)/(S**2+S02)
	RETURN
 2	IF(S.GT.SS(2)) GO TO 21	
        FIR=A(1)-(A(1)-A(2))*S*S/SS(2)**2
        RETURN
 21     N=1
	GO TO 4
 3	N=N+1
 4	IF(S.LE.SS(N)) GO TO 5
	GO TO 3
 5	FIR=A(N-1)+(S-SS(N-1))*(A(N)-A(N-1))/(SS(N)-SS(N-1))
	RETURN
	END
	


	REAL*8 FUNCTION FK(X)
C Program for evaluating the scattering amplitude for  POTASSIUM atom, inter-
C polation of data from International Tables for X-ray Crystallography, 1974
C argument X equals to q/4pi within the notations of Landau & Lifshits
	IMPLICIT REAL*8(A-H,O-Z)
	DIMENSION A(56),XX(56)
	PARAMETER (CF=41.78214D0)
	DATA Z/19.D0/
       DATA A/8.984D0,8.921D0,8.731D0,8.434D0,8.054D0,7.619D0,
     ,7.157D0,6.691D0,6.239D0,5.815D0,5.426D0,5.073D0,4.756D0,
     ,4.474D0,4.222D0,3.997D0,3.795D0,3.612D0,3.446D0,3.295D0,
     ,3.154D0,2.902D0,2.680D0,2.578D0,2.481D0,2.299D0,2.134D0,
     ,1.982D0,1.842D0,1.776D0,1.714D0,1.595D0,1.487D0,1.387D0,
     ,1.295D0,1.252D0,1.211D0,1.134D0,1.064D0,0.912D0,0.790D0,
     ,0.690d0,0.609D0,0.488D0,0.402D0,0.339D0,0.290D0,0.252D0,
     ,0.220D0,0.194D0,0.174D0,0.156D0,0.138D0,0.127D0,0.112D0,
     ,0.101D0/
C
       DATA XX/0.D0,0.01D0,0.02D0,0.03D0,0.04d0,0.05d0,0.06d0,
     ,0.07d0,0.08d0,0.09d0,0.10d0,0.11d0,0.12d0,0.13d0,0.14d0,
     ,0.15d0,0.16d0,0.17d0,0.18d0,0.19d0,0.20d0,0.22d0,0.24d0,
     ,0.25D0,0.26d0,0.28d0,0.30d0,0.32d0,0.34d0,0.35d0,0.36d0,
     ,0.38d0,0.40d0,0.42d0,0.44d0,0.45d0,0.46d0,0.48d0,0.50d0,
     ,0.55d0,0.60d0,0.65d0,0.70d0,0.80d0,0.90d0,1.00d0,1.10d0,
     ,1.20d0,1.30d0,1.40d0,1.50d0,1.60d0,1.70d0,1.80d0,1.90d0,
     ,2.00d0/
         X02=Z/(CF*A(56))-XX(56)**2
	IF(X.GT.XX(1)) GO TO 1
	FK=A(1)
	RETURN
 1	IF(X.LT.XX(56)) GO TO 2
C for large arguments the screened Coulomb law is assumed
	FK=(Z/CF)/(X**2+X02)
	RETURN
 2	IF(X.GT.XX(2)) GO TO 21
 	FK=A(1)-(A(1)-A(2))*X*X/XX(2)**2
 	RETURN
 21	N=1
	GO TO 4
 3	N=N+1
 4	IF(X.LE.XX(N)) GO TO 5
	GO TO 3
 5	FK=A(N-1)+(X-XX(N-1))*(A(N)-A(N-1))/(XX(N)-XX(N-1))
	RETURN
	END

	REAL*8 FUNCTION FLA(X)
C Program for evaluating the scattering amplitude for LANTAN atom, inter-
C polation of data from International Tables for X-ray Crystallography, 1974
C argument X equals to q/4pi within the notations of Landau & Lifshits
	IMPLICIT REAL*8(A-H,O-Z)
	DIMENSION A(53),XX(53)
	PARAMETER (CF=41.78214D0)
	DATA Z/57.D0/
       DATA A/17.805D0,16.45D0,15.79D0,15.05D0,
     ,14.28D0,13.51D0,12.74D0,12.01D0,11.32D0,10.671D0,10.072D0,
     ,9.522D0,9.017D0,8.555D0,8.131D0,7.742D0,7.384D0,7.053D0, 
     ,6.462D0,5.948D0,5.714D0,5.495D0,5.092D0,4.730D0,4.405D0,
     ,4.111D0,3.974D0,3.844D0,3.602D0,3.381D0,3.180D0,2.997D0,
     ,2.911D0,2.829D0,2.676D0,2.535D0,2.230D0,1.979D0,1.770d0,
     ,1.592D0,1.308D0,1.090D0,0.920D0,0.785D0,0.678D0,0.591D0,
     ,0.521D0,0.463D0,0.415D0,0.374D0,0.340D0,0.310D0,0.284D0/
C
       DATA XX/0.D0,0.04d0,0.05d0,0.06d0,
     ,0.07d0,0.08d0,0.09d0,0.10d0,0.11d0,0.12d0,0.13d0,0.14d0,
     ,0.15d0,0.16d0,0.17d0,0.18d0,0.19d0,0.20d0,0.22d0,0.24d0,
     ,0.25D0,0.26d0,0.28d0,0.30d0,0.32d0,0.34d0,0.35d0,0.36d0,
     ,0.38d0,0.40d0,0.42d0,0.44d0,0.45d0,0.46d0,0.48d0,0.50d0,
     ,0.55d0,0.60d0,0.65d0,0.70d0,0.80d0,0.90d0,1.00d0,1.10d0,
     ,1.20d0,1.30d0,1.40d0,1.50d0,1.60d0,1.70d0,1.80d0,1.90d0,
     ,2.00d0/
        X02=Z/(CF*A(53))-XX(53)**2
	IF(X.GT.XX(1)) GO TO 1
	FLA=A(1)
	RETURN
 1	IF(X.LT.XX(53)) GO TO 2
C for large arguments the screened Coulomb law is assumed
	FLA=(Z/CF)/(X**2+X02)
	RETURN
 2	IF(X.GT.XX(2)) GO TO 21
 	FLA=A(1)-(A(1)-A(2))*X*X/XX(2)**2
 	RETURN
 21	N=1
	GO TO 4
 3	N=N+1
 4	IF(X.LE.XX(N)) GO TO 5
	GO TO 3
 5	FLA=A(N-1)+(X-XX(N-1))*(A(N)-A(N-1))/(XX(N)-XX(N-1))
	RETURN
	END
 
 	REAL*8 FUNCTION FLI(X)
C Program for evaluating the scattering amplitude for  LITHIUM atom, inter-
C polation of data from International Tables for X-ray Crystallography, 1974
C argument X equals to q/4pi within the notations of Landau & Lifshits
	IMPLICIT REAL*8(A-H,O-Z)
	DIMENSION A(56),XX(56)
	PARAMETER (CF=41.78214D0)
	DATA Z/3.D0/
       DATA A/3.286D0,3.265D0,3.200D0,3.097D0,2.961D0,2.800D0,
     ,2.622D0,2.435D0,2.245D0,2.058D0,1.879D0,1.710D0,1.554D0,
     ,1.411D0,1.282D0,1.166D0,1.063D0,0.971D0,0.889D0,0.817D0,
     ,0.753D0,0.646D0,0.562D0,0.526D0,0.494D0,0.440D0,0.396D0,
     ,0.359D0,0.328D0,0.314D0,0.301D0,0.279D0,0.259D0,0.241D0,
     ,0.226D0,0.219D0,0.212D0,0.200D0,0.188D0,0.164D0,0.145D0,
     ,0.128d0,0.115D0,0.093D0,0.077D0,0.064D0,0.054D0,0.046D0,
     ,0.040D0,0.035D0,0.031D0,0.028D0,0.024D0,0.022D0,0.019D0,
     ,0.017D0/
C
       DATA XX/0.D0,0.01D0,0.02D0,0.03D0,0.04d0,0.05d0,0.06d0,
     ,0.07d0,0.08d0,0.09d0,0.10d0,0.11d0,0.12d0,0.13d0,0.14d0,
     ,0.15d0,0.16d0,0.17d0,0.18d0,0.19d0,0.20d0,0.22d0,0.24d0,
     ,0.25D0,0.26d0,0.28d0,0.30d0,0.32d0,0.34d0,0.35d0,0.36d0,
     ,0.38d0,0.40d0,0.42d0,0.44d0,0.45d0,0.46d0,0.48d0,0.50d0,
     ,0.55d0,0.60d0,0.65d0,0.70d0,0.80d0,0.90d0,1.00d0,1.10d0,
     ,1.20d0,1.30d0,1.40d0,1.50d0,1.60d0,1.70d0,1.80d0,1.90d0,
     ,2.00d0/
         X02=Z/(CF*A(56))-XX(56)**2
	IF(X.GT.XX(1)) GO TO 1
	FLI=A(1)
	RETURN
 1	IF(X.LT.XX(56)) GO TO 2
C for large arguments the screened Coulomb law is assumed
	FLI=(Z/CF)/(X**2+X02)
	RETURN
 2	IF(X.GT.XX(2)) GO TO 21
 	FLI=A(1)-(A(1)-A(2))*X*X/XX(2)**2
 	RETURN
 21	N=1
	GO TO 4
 3	N=N+1
 4	IF(X.LE.XX(N)) GO TO 5
	GO TO 3
 5	FLI=A(N-1)+(X-XX(N-1))*(A(N)-A(N-1))/(XX(N)-XX(N-1))
	RETURN
	END
	
	REAL*8 FUNCTION FMG(X)
C Program for evaluating the scattering amplitude for a MAGNESIUM atom, inter-
C polation of data from International Tables for X-ray Crystallography, 1974
C argument X equals to q/4pi within the notations of Landau & Lifshits
	IMPLICIT REAL*8(A-H,O-Z)
	DIMENSION A(56),XX(56)
	PARAMETER (CF=41.78214D0)
	DATA Z/12.D0/
       DATA A/5.207D0,5.187D0,5.124D0,5.022D0,4.884D0,4.717D0,
     ,4.527D0,4.320D0,4.102D0,3.879D0,3.656D0,3.437D0,3.226D0,
     ,3.025D0,2.835D0,2.657D0,2.492D0,2.340D0,2.199D0,2.071D0,
     ,1.953D0,1.748D0,1.577D0,1.502D0,1.434D0,1.313D0,1.211D0,
     ,1.123D0,1.047D0,1.013D0,0.980D0,0.921D0,0.868D0,0.821D0,
     ,0.777D0,0.757D0,0.738D0,0.701D0,0.667D0,0.592D0,0.528D0,
     ,0.473d0,0.425D0,0.347D0,0.286D0,0.239D0,0.202D0,0.172D0,
     ,0.148D0,0.129D0,0.113D0,0.100D0,0.089D0,0.080D0,0.072D0,
     ,0.065D0/
C
       DATA XX/0.D0,0.01D0,0.02D0,0.03D0,0.04d0,0.05d0,0.06d0,
     ,0.07d0,0.08d0,0.09d0,0.10d0,0.11d0,0.12d0,0.13d0,0.14d0,
     ,0.15d0,0.16d0,0.17d0,0.18d0,0.19d0,0.20d0,0.22d0,0.24d0,
     ,0.25D0,0.26d0,0.28d0,0.30d0,0.32d0,0.34d0,0.35d0,0.36d0,
     ,0.38d0,0.40d0,0.42d0,0.44d0,0.45d0,0.46d0,0.48d0,0.50d0,
     ,0.55d0,0.60d0,0.65d0,0.70d0,0.80d0,0.90d0,1.00d0,1.10d0,
     ,1.20d0,1.30d0,1.40d0,1.50d0,1.60d0,1.70d0,1.80d0,1.90d0,
     ,2.00d0/
         X02=Z/(CF*A(56))-XX(56)**2
	IF(X.GT.XX(1)) GO TO 1
	FMG=A(1)
	RETURN
 1	IF(X.LT.XX(56)) GO TO 2
C for large arguments the screened Coulomb law is assumed
	FMG=(Z/CF)/(X**2+X02)
	RETURN
 2	IF(X.GT.XX(2)) GO TO 21
 	FMG=A(1)-(A(1)-A(2))*X*X/XX(2)**2
 	RETURN
 21	N=1
	GO TO 4
 3	N=N+1
 4	IF(X.LE.XX(N)) GO TO 5
	GO TO 3
 5	FMG=A(N-1)+(X-XX(N-1))*(A(N)-A(N-1))/(XX(N)-XX(N-1))
	RETURN
	END

	REAL*8 FUNCTION FMN(X)
C Program for evaluating the scattering amplitude for a MANGANESE atom, inter-
C polation of data from International Tables for X-ray Crystallography, 1974
C argument X equals to q/4pi within the notations of Landau & Lifshits
	IMPLICIT REAL*8(A-H,O-Z)
	DIMENSION A(56),XX(56)
	PARAMETER (CF=41.78214D0)
	DATA Z/25.D0/
       DATA A/7.506D0,7.484D0,7.412D0,7.296D0,7.140D0,6.949D0,
     ,6.732D0,6.493D0,6.241D0,5.981D0,5.719D0,5.459D0,5.206D0,
     ,4.962D0,4.728D0,4.506D0,4.297D0,4.100D0,3.916D0,3.743D0,
     ,3.583D0,3.292D0,3.039D0,2.924D0,2.817D0,2.620D0,2.445D0,
     ,2.288D0,2.146D0,2.080D0,2.017D0,1.899D0,1.790D0,1.690D0,
     ,1.597D0,1.553D0,1.511D0,1.431D0,1.356D0,1.189D0,1.047D0,
     ,0.927d0,0.824D0,0.659D0,0.538D0,0.446D0,0.377D0,0.323D0,
     ,0.280D0,0.246D0,0.218D0,0.195D0,0.175D0,0.159D0,0.144D0,
     ,0.132D0/
C
       DATA XX/0.D0,0.01D0,0.02D0,0.03D0,0.04d0,0.05d0,0.06d0,
     ,0.07d0,0.08d0,0.09d0,0.10d0,0.11d0,0.12d0,0.13d0,0.14d0,
     ,0.15d0,0.16d0,0.17d0,0.18d0,0.19d0,0.20d0,0.22d0,0.24d0,
     ,0.25D0,0.26d0,0.28d0,0.30d0,0.32d0,0.34d0,0.35d0,0.36d0,
     ,0.38d0,0.40d0,0.42d0,0.44d0,0.45d0,0.46d0,0.48d0,0.50d0,
     ,0.55d0,0.60d0,0.65d0,0.70d0,0.80d0,0.90d0,1.00d0,1.10d0,
     ,1.20d0,1.30d0,1.40d0,1.50d0,1.60d0,1.70d0,1.80d0,1.90d0,
     ,2.00d0/
         X02=Z/(CF*A(56))-XX(56)**2
	IF(X.GT.XX(1)) GO TO 1
	FMN=A(1)
	RETURN
 1	IF(X.LT.XX(56)) GO TO 2
C for large arguments the screened Coulomb law is assumed
	FMN=(Z/CF)/(X**2+X02)
	RETURN
 2	IF(X.GT.XX(2)) GO TO 21
 	FMN=A(1)-(A(1)-A(2))*X*X/XX(2)**2
 	RETURN
 21	N=1
	GO TO 4
 3	N=N+1
 4	IF(X.LE.XX(N)) GO TO 5
	GO TO 3
 5	FMN=A(N-1)+(X-XX(N-1))*(A(N)-A(N-1))/(XX(N)-XX(N-1))
	RETURN
	END
 
 	REAL*8 FUNCTION FMO(X)
C Program for evaluating the scattering amplitude for MOLIBDENUM atom, inter-
C polation of data from International Tables for X-ray Crystallography, 1974
C argument X equals to q/4pi within the notations of Landau & Lifshits
	IMPLICIT REAL*8(A-H,O-Z)
	DIMENSION A(56),XX(56)
	PARAMETER (CF=41.78214D0)
	DATA Z/42.D0/
       DATA A/10.260D0,10.230D0,10.138D0,9.989D0,9.790D0,9.548D0,
     ,9.272D0,8.972D0,8.655D0,8.330D0,8.004D0,7.680D0,7.364D0,
     ,7.058D0,6.763D0,6.481D0,6.213D0,5.957D0,5.715D0,5.486D0,
     ,5.269D0,4.868D0,4.507D0,4.341D0,4.182D0,3.888D0,3.622D0,
     ,3.379D0,3.158D0,3.054D0,2.955D0,2.770D0,2.600D0,2.444D0,
     ,2.300D0,2.233D0,2.168D0,2.046D0,1.934D0,1.690D0,1.490D0,
     ,1.324d0,1.185D0,0.971D0,0.814D0,0.695D0,0.601D0,0.525D0,
     ,0.462D0,0.408D0,0.364D0,0.325D0,0.291D0,0.263D0,0.238D0,
     ,0.216D0/
C
       DATA XX/0.D0,0.01D0,0.02D0,0.03D0,0.04d0,0.05d0,0.06d0,
     ,0.07d0,0.08d0,0.09d0,0.10d0,0.11d0,0.12d0,0.13d0,0.14d0,
     ,0.15d0,0.16d0,0.17d0,0.18d0,0.19d0,0.20d0,0.22d0,0.24d0,
     ,0.25D0,0.26d0,0.28d0,0.30d0,0.32d0,0.34d0,0.35d0,0.36d0,
     ,0.38d0,0.40d0,0.42d0,0.44d0,0.45d0,0.46d0,0.48d0,0.50d0,
     ,0.55d0,0.60d0,0.65d0,0.70d0,0.80d0,0.90d0,1.00d0,1.10d0,
     ,1.20d0,1.30d0,1.40d0,1.50d0,1.60d0,1.70d0,1.80d0,1.90d0,
     ,2.00d0/
        X02=Z/(CF*A(56))-XX(56)**2
	IF(X.GT.XX(1)) GO TO 1
	FMO=A(1)
	RETURN
 1	IF(X.LT.XX(56)) GO TO 2
C for large arguments the screened Coulomb law is assumed
	FMO=(Z/CF)/(X**2+X02)
	RETURN
 2	IF(X.GT.XX(2)) GO TO 21
 	FMO=A(1)-(A(1)-A(2))*X*X/XX(2)**2
 	RETURN
 21	N=1
	GO TO 4
 3	N=N+1
 4	IF(X.LE.XX(N)) GO TO 5
	GO TO 3
 5	FMO=A(N-1)+(X-XX(N-1))*(A(N)-A(N-1))/(XX(N)-XX(N-1))
	RETURN
	END
	REAL*8 FUNCTION FNA(X)
C Program for evaluating the scattering amplitude for  SODIUM atom, inter-
C polation of data from International Tables for X-ray Crystallography, 1974
C argument X equals to q/4pi within the notations of Landau & Lifshits
	IMPLICIT REAL*8(A-H,O-Z)
	DIMENSION A(56),XX(56)
	PARAMETER (CF=41.78214D0)
	DATA Z/11.D0/
       DATA A/4.778D0,4.749D0,4.663D0,4.527D0,4.348D0,4.138D0,
     ,3.908D0,3.667D0,3.425D0,3.190D0,2.967D0,2.759D0,2.569D0,
     ,2.395D0,2.239D0,2.099D0,1.974D0,1.863D0,1.763D0,1.674D0,
     ,1.594D0,1.458D0,1.344D0,1.295D0,1.249D0,1.167D0,1.095D0,
     ,1.031D0,0.973D0,0.946D0,0.921D0,0.872D0,0.827D0,0.785D0,
     ,0.746D0,0.727D0,0.709D0,0.675D0,0.642D0,0.569D0,0.505D0,
     ,0.450d0,0.403D0,0.325D0,0.266D0,0.221D0,0.185D0,0.158D0,
     ,0.135D0,0.117D0,0.103D0,0.092D0,0.081D0,0.073D0,0.065D0,
     ,0.059D0/
C
       DATA XX/0.D0,0.01D0,0.02D0,0.03D0,0.04d0,0.05d0,0.06d0,
     ,0.07d0,0.08d0,0.09d0,0.10d0,0.11d0,0.12d0,0.13d0,0.14d0,
     ,0.15d0,0.16d0,0.17d0,0.18d0,0.19d0,0.20d0,0.22d0,0.24d0,
     ,0.25D0,0.26d0,0.28d0,0.30d0,0.32d0,0.34d0,0.35d0,0.36d0,
     ,0.38d0,0.40d0,0.42d0,0.44d0,0.45d0,0.46d0,0.48d0,0.50d0,
     ,0.55d0,0.60d0,0.65d0,0.70d0,0.80d0,0.90d0,1.00d0,1.10d0,
     ,1.20d0,1.30d0,1.40d0,1.50d0,1.60d0,1.70d0,1.80d0,1.90d0,
     ,2.00d0/
         X02=Z/(CF*A(56))-XX(56)**2
	IF(X.GT.XX(1)) GO TO 1
	FNA=A(1)
	RETURN
 1	IF(X.LT.XX(56)) GO TO 2
C for large arguments the screened Coulomb law is assumed
	FNA=(Z/CF)/(X**2+X02)
	RETURN
 2	IF(X.GT.XX(2)) GO TO 21
 	FNA=A(1)-(A(1)-A(2))*X*X/XX(2)**2
 	RETURN
 21	N=1
	GO TO 4
 3	N=N+1
 4	IF(X.LE.XX(N)) GO TO 5
	GO TO 3
 5	FNA=A(N-1)+(X-XX(N-1))*(A(N)-A(N-1))/(XX(N)-XX(N-1))
	RETURN
	END
c **************************************************************************

	REAL*8 FUNCTION FNB(S)
C Program for evaluation of scattering amplitude for a NIOBIUM atom, inter-
C polation of data from International Tables for X-ray Crystallography, 1974
C argument S equals to q/(4*pi) within the notation of Landau & Lifshits
	IMPLICIT REAL*8(A-H,O-Z)
	DIMENSION A(53),SS(53)
	PARAMETER (CF=41.78214D0)
	DATA Z/41.D0/
       DATA A/10.679D0,10.13D0,9.86D0,
     ,9.54D0,9.20D0,8.85D0, 8.49D0, 8.12D0, 7.77D0, 7.421D0,
     ,7.090d0,6.772d0,6.472d0,6.187d0,5.918d0,5.665d0,5.427d0,
     ,5.203D0,4.792D0,4.426D0,4.258D0,4.099D0,3.804D0,3.539D0,
     ,3.298D0,3.080D0,2.978D0,2.880D0,2.698D0,2.531D0,2.379D0,
     ,2.239D0,2.173D0,2.110D0,1.991D0,1.883D0,1.646D0,1.452D0,
     ,1.292D0,1.159D0,0.952D0,0.800D0,0.684D0,0.591D0,0.516D0,
     ,0.454D0,0.401D0,0.356D0,0.318D0,0.285D0,0.257D0,0.233D0,
     ,0.211D0/
C
       DATA SS/0.D0,0.04d0,0.05d0,0.06d0,
     ,0.07d0,0.08d0,0.09d0,0.10d0,0.11d0,0.12d0,0.13d0,0.14d0,
     ,0.15d0,0.16d0,0.17d0,0.18d0,0.19d0,0.20d0,0.22d0,0.24d0,
     ,0.25D0,0.26d0,0.28d0,0.30d0,0.32d0,0.34d0,0.35d0,0.36d0,
     ,0.38d0,0.40d0,0.42d0,0.44d0,0.45d0,0.46d0,0.48d0,0.50d0,
     ,0.55d0,0.60d0,0.65d0,0.70d0,0.80d0,0.90d0,1.00d0,1.10d0,
     ,1.20d0,1.30d0,1.40d0,1.50d0,1.60d0,1.70d0,1.80d0,1.90d0,
     ,2.00d0/
        S02=Z/(CF*A(53))-SS(53)**2
	IF(S.GT.SS(1)) GO TO 1
	FNB=A(1)
	RETURN
 1	IF(S.LT.SS(53)) GO TO 2
C for large arguments the screened Coulomb law is assumed
	FNB=(Z/CF)/(S**2+S02)
	RETURN
 2	IF(S.GT.SS(2)) GO TO 21	
        FNB=A(1)-(A(1)-A(2))*S*S/SS(2)**2
        RETURN
 21     N=1
	GO TO 4
 3	N=N+1
 4	IF(S.LE.SS(N)) GO TO 5
	GO TO 3
 5	FNB=A(N-1)+(S-SS(N-1))*(A(N)-A(N-1))/(SS(N)-SS(N-1))
	RETURN
	END
	


	REAL*8 FUNCTION FNI(X)
C Program for evaluating the scattering amplitude for NICKEL atom, inter-
C polation of data from International Tables for X-ray Crystallography, 1974
C argument X equals to q/4pi within the notations of Landau & Lifshits
	IMPLICIT REAL*8(A-H,O-Z)
	DIMENSION A(56),XX(56)
	PARAMETER (CF=41.78214D0)
	DATA Z/28.D0/
       DATA A/6.569D0,6.552D0,6.501D0,6.418D0,6.306D0,6.169D0,
     ,6.010D0,5.834D0,5.646D0,5.449D0,5.249D0,5.048D0,4.848D0,
     ,4.654D0,4.565D0,4.283D0,4.110D0,3.944D0,3.788D0,3.640D0,
     ,3.500D0,3.245D0,3.018D0,2.914D0,2.816D0,2.636D0,2.474D0,
     ,2.328D0,2.195D0,2.133D0,2.073D0,1.962D0,1.858D0,1.763D0,
     ,1.674D0,1.631D0,1.591D0,1.513D0,1.440D0,1.277D0,1.136D0,
     ,1.015d0,0.909D0,0.737D0,0.605D0,0.504D0,0.425D0,0.364D0,
     ,0.315D0,0.275D0,0.243D0,0.217D0,0.194D0,0.176D0,0.160D0,
     ,0.146D0/
C
       DATA XX/0.D0,0.01D0,0.02D0,0.03D0,0.04d0,0.05d0,0.06d0,
     ,0.07d0,0.08d0,0.09d0,0.10d0,0.11d0,0.12d0,0.13d0,0.14d0,
     ,0.15d0,0.16d0,0.17d0,0.18d0,0.19d0,0.20d0,0.22d0,0.24d0,
     ,0.25D0,0.26d0,0.28d0,0.30d0,0.32d0,0.34d0,0.35d0,0.36d0,
     ,0.38d0,0.40d0,0.42d0,0.44d0,0.45d0,0.46d0,0.48d0,0.50d0,
     ,0.55d0,0.60d0,0.65d0,0.70d0,0.80d0,0.90d0,1.00d0,1.10d0,
     ,1.20d0,1.30d0,1.40d0,1.50d0,1.60d0,1.70d0,1.80d0,1.90d0,
     ,2.00d0/
         X02=Z/(CF*A(56))-XX(56)**2
	IF(X.GT.XX(1)) GO TO 1
	FNI=A(1)
	RETURN
 1	IF(X.LT.XX(56)) GO TO 2
C for large arguments the screened Coulomb law is assumed
	FNI=(Z/CF)/(X**2+X02)
	RETURN
 2	IF(X.GT.XX(2)) GO TO 21
 	FNI=A(1)-(A(1)-A(2))*X*X/XX(2)**2
 	RETURN
 21	N=1
	GO TO 4
 3	N=N+1
 4	IF(X.LE.XX(N)) GO TO 5
	GO TO 3
 5	FNI=A(N-1)+(X-XX(N-1))*(A(N)-A(N-1))/(XX(N)-XX(N-1))
	RETURN
	END
	REAL*8 FUNCTION FO(X)
C Program for evaluating the scattering amplitude for an OXIGEN atom, inter-
C polation of data from International Tables for X-ray Crystallography, 1974
C argument X equals to q/4pi within the notations of Landau & Lifshits
	IMPLICIT REAL*8(A-H,O-Z)
	DIMENSION A(56),XX(56)
	PARAMETER (CF=41.78214D0)
	DATA Z/8.D0/
       DATA A/1.983D0,1.982D0,1.976D0,1.966D0,1.953D0,1.937D0,
     ,1.917D0,1.893D0,1.867D0,1.839D0,1.808D0,1.774D0,1.739D0,
     ,1.702D0,1.664D0,1.625D0,1.585D0,1.545D0,1.504D0,1.463D0,
     ,1.422D0,1.341D0,1.261D0,1.222D0,1.184D0,1.110D0,1.040D0,
     ,0.974D0,0.911D0,0.881D0,0.853D0,0.798D0,0.747D0,0.700D0,
     ,0.656D0,0.635D0,0.615D0,0.577D0,0.542D0,0.466D0,0.403D0,
     ,0.350d0,0.307D0,0.241D0,0.193D0,0.159D0,0.133D0,0.113D0,
     ,0.097D0,0.085D0,0.074D0,0.066D0,0.059D0,0.053D0,0.048D0,
     ,0.044D0/
C
       DATA XX/0.D0,0.01D0,0.02D0,0.03D0,0.04d0,0.05d0,0.06d0,
     ,0.07d0,0.08d0,0.09d0,0.10d0,0.11d0,0.12d0,0.13d0,0.14d0,
     ,0.15d0,0.16d0,0.17d0,0.18d0,0.19d0,0.20d0,0.22d0,0.24d0,
     ,0.25D0,0.26d0,0.28d0,0.30d0,0.32d0,0.34d0,0.35d0,0.36d0,
     ,0.38d0,0.40d0,0.42d0,0.44d0,0.45d0,0.46d0,0.48d0,0.50d0,
     ,0.55d0,0.60d0,0.65d0,0.70d0,0.80d0,0.90d0,1.00d0,1.10d0,
     ,1.20d0,1.30d0,1.40d0,1.50d0,1.60d0,1.70d0,1.80d0,1.90d0,
     ,2.00d0/
         X02=Z/(CF*A(56))-XX(56)**2
	IF(X.GT.XX(1)) GO TO 1
	FO=A(1)
	RETURN
 1	IF(X.LT.XX(56)) GO TO 2
C for large arguments the screened Coulomb law is assumed
	FO=(Z/CF)/(X**2+X02)
	RETURN
 2	IF(X.GT.XX(2)) GO TO 21
 	FO=A(1)-(A(1)-A(2))*X*X/XX(2)**2
 	RETURN
 21	N=1
	GO TO 4
 3	N=N+1
 4	IF(X.LE.XX(N)) GO TO 5
	GO TO 3
 5	FO=A(N-1)+(X-XX(N-1))*(A(N)-A(N-1))/(XX(N)-XX(N-1))
	RETURN
	END
c **************************************************************************

	REAL*8 FUNCTION FOS(S)
C Program for evaluation of scattering amplitude for a OSMIUM atom, inter-
C polation of data from International Tables for X-ray Crystallography, 1974
C argument S equals to q/(4*pi) within the notation of Landau & Lifshits
	IMPLICIT REAL*8(A-H,O-Z)
	DIMENSION A(53),SS(53)
	PARAMETER (CF=41.78214D0)
	DATA Z/76.D0/
       DATA A/11.987d0,11.59D0,11.39D0,
     ,11.15D0,10.88D0,10.59D0,10.29D0,9.98D0,9.65D0,9.334D0,
     ,9.016D0,8.702D0,8.396D0,8.099D0,7.813D0,7.537D0,7.272D0,
     ,7.019D0,6.547D0,6.117D0,5.917D0,5.727D0,5.372D0,5.049D0,
     ,4.755D0,4.485D0,4.359D0,4.237D0,4.010D0,3.800D0,3.606D0,
     ,3.427D0,3.342D0,3.260D0,3.105D0,2.961D0,2.641D0,2.371D0,
     ,2.140D0,1.942D0,1.621D0,1.374D0,1.179D0,1.022D0,0.892D0,
     ,0.784D0,0.694D0,0.617D0,0.551D0,0.495D0,0.447D0,0.406D0,
     ,0.370D0/
C
       DATA SS/0.D0,0.04d0,0.05d0,0.06d0,
     ,0.07d0,0.08d0,0.09d0,0.10d0,0.11d0,0.12d0,0.13d0,0.14d0,
     ,0.15d0,0.16d0,0.17d0,0.18d0,0.19d0,0.20d0,0.22d0,0.24d0,
     ,0.25D0,0.26d0,0.28d0,0.30d0,0.32d0,0.34d0,0.35d0,0.36d0,
     ,0.38d0,0.40d0,0.42d0,0.44d0,0.45d0,0.46d0,0.48d0,0.50d0,
     ,0.55d0,0.60d0,0.65d0,0.70d0,0.80d0,0.90d0,1.00d0,1.10d0,
     ,1.20d0,1.30d0,1.40d0,1.50d0,1.60d0,1.70d0,1.80d0,1.90d0,
     ,2.00d0/
        S02=Z/(CF*A(53))-SS(53)**2
	IF(S.GT.SS(1)) GO TO 1
	FOS=A(1)
	RETURN
 1	IF(S.LT.SS(53)) GO TO 2
C for large arguments the screened Coulomb law is assumed
	FOS=(Z/CF)/(S**2+S02)
	RETURN
 2	IF(S.GT.SS(2)) GO TO 21	
        FOS=A(1)-(A(1)-A(2))*S*S/SS(2)**2
        RETURN
 21     N=1
	GO TO 4
 3	N=N+1
 4	IF(S.LE.SS(N)) GO TO 5
	GO TO 3
 5	FOS=A(N-1)+(S-SS(N-1))*(A(N)-A(N-1))/(SS(N)-SS(N-1))
	RETURN
	END
	



	REAL*8 FUNCTION FPB(X)
C Program for evaluating the scattering amplitude for a LEAD atom, inter-
C polation of data from International Tables for X-ray Crystallography, 1974
C argument X equals to q/4pi within the notations of Landau & Lifshits
	IMPLICIT REAL*8(A-H,O-Z)
	DIMENSION A(56),XX(56)
	PARAMETER (CF=41.78214D0)
	DATA Z/82.D0/
       DATA A/12.597D0,12.573D0,12.494D0,12.366D0,12.193D0,11.979D0,
     ,11.730D0,11.454D0,11.155D0,10.840D0,10.516D0,10.186D0,9.855D0,
     ,9.527D0,9.203D0,8.888D0,8.581D0,8.285D0,7.999D0,7.724D0,
     ,7.461D0,6.969D0,6.520D0,6.310D0,6.110D0,5.736D0,5.395D0,
     ,5.083D0,4.797D0,4.662D0,4.533D0,4.290D0,4.066D0,3.858D0,
     ,3.665D0,3.573D0,3.485D0,3.318D0,3.162D0,2.816D0,2.522D0,
     ,2.271d0,2.055D0,1.708D0,1.444D0,1.239D0,1.075D0,0.942D0,
     ,0.831D0,0.738D0,0.659D0,0.591D0,0.533D0,0.482D0,0.438D0,
     ,0.399D0/
C
       DATA XX/0.D0,0.01D0,0.02D0,0.03D0,0.04d0,0.05d0,0.06d0,
     ,0.07d0,0.08d0,0.09d0,0.10d0,0.11d0,0.12d0,0.13d0,0.14d0,
     ,0.15d0,0.16d0,0.17d0,0.18d0,0.19d0,0.20d0,0.22d0,0.24d0,
     ,0.25D0,0.26d0,0.28d0,0.30d0,0.32d0,0.34d0,0.35d0,0.36d0,
     ,0.38d0,0.40d0,0.42d0,0.44d0,0.45d0,0.46d0,0.48d0,0.50d0,
     ,0.55d0,0.60d0,0.65d0,0.70d0,0.80d0,0.90d0,1.00d0,1.10d0,
     ,1.20d0,1.30d0,1.40d0,1.50d0,1.60d0,1.70d0,1.80d0,1.90d0,
     ,2.00d0/
        X02=Z/(CF*A(56))-XX(56)**2
	IF(X.GT.XX(1)) GO TO 1
	FPB=A(1)
	RETURN
 1	IF(X.LT.XX(56)) GO TO 2
C for large arguments the screened Coulomb law is assumed
	FPB=(Z/cf)/(X**2+X02)
	RETURN
 2	IF(X.GT.XX(2)) GO TO 21
 	FPB=A(1)-(A(1)-A(2))*X*X/XX(2)**2
 	RETURN
 21	N=1
	GO TO 4
 3	N=N+1
 4	IF(X.LE.XX(N)) GO TO 5
	GO TO 3
 5	FPB=A(N-1)+(X-XX(N-1))*(A(N)-A(N-1))/(XX(N)-XX(N-1))
	RETURN
	END
c
c **************************************************************************

	REAL*8 FUNCTION FPD(S)
C Program for evaluation of scattering amplitude for a PALLADIUM atom, inter-
C polation of data from International Tables for X-ray Crystallography, 1974
C argument S equals to q/(4*pi) within the notation of Landau & Lifshits
	IMPLICIT REAL*8(A-H,O-Z)
	DIMENSION A(53),SS(53)
	PARAMETER (CF=41.78214D0)
	DATA Z/46.D0/
       DATA A/7.583D0,7.43D0,7.35D0,
     ,7.26D0,7.16D0,7.03D0, 6.91D0, 6.77D0, 6.62D0, 6.474D0,
     ,6.319d0,6.162d0,6.003d0,5.843d0,5.684d0,5.526d0,5.369d0,
     ,5.214D0,4.913D0,4.626D0,4.487D0,4.352D0,4.093D0,3.850D0,
     ,3.622D0,3.408D0,3.306D0,3.208D0,3.022D0,2.848D0,2.686D0,
     ,2.535D0,2.464D0,2.395D0,2.264D0,2.143D0,1.875D0,1.650D0,
     ,1.462D0,1.304D0,1.058D0,0.879D0,0.745D0,0.641D0,0.559D0,
     ,0.493D0,0.437D0,0.391D0,0.351D0,0.316D0,0.286D0,0.260D0,
     ,0.237D0/
C
       DATA SS/0.D0,0.04d0,0.05d0,0.06d0,
     ,0.07d0,0.08d0,0.09d0,0.10d0,0.11d0,0.12d0,0.13d0,0.14d0,
     ,0.15d0,0.16d0,0.17d0,0.18d0,0.19d0,0.20d0,0.22d0,0.24d0,
     ,0.25D0,0.26d0,0.28d0,0.30d0,0.32d0,0.34d0,0.35d0,0.36d0,
     ,0.38d0,0.40d0,0.42d0,0.44d0,0.45d0,0.46d0,0.48d0,0.50d0,
     ,0.55d0,0.60d0,0.65d0,0.70d0,0.80d0,0.90d0,1.00d0,1.10d0,
     ,1.20d0,1.30d0,1.40d0,1.50d0,1.60d0,1.70d0,1.80d0,1.90d0,
     ,2.00d0/
        S02=Z/(CF*A(53))-SS(53)**2
	IF(S.GT.SS(1)) GO TO 1
	FPD=A(1)
	RETURN
 1	IF(S.LT.SS(53)) GO TO 2
C for large arguments the screened Coulomb law is assumed
	FPD=(Z/CF)/(S**2+S02)
	RETURN
 2	IF(S.GT.SS(2)) GO TO 21	
        FPD=A(1)-(A(1)-A(2))*S*S/SS(2)**2
        RETURN
 21     N=1
	GO TO 4
 3	N=N+1
 4	IF(S.LE.SS(N)) GO TO 5
	GO TO 3
 5	FPD=A(N-1)+(S-SS(N-1))*(A(N)-A(N-1))/(SS(N)-SS(N-1))
	RETURN
	END
	


c **************************************************************************

	REAL*8 FUNCTION FPT(S)
C Program for evaluation of scattering amplitude for a PLATINUM atom, inter-
C polation of data from International Tables for X-ray Crystallography, 1974
C argument S equals to q/(4*pi) within the notation of Landau & Lifshits
	IMPLICIT REAL*8(A-H,O-Z)
	DIMENSION A(53),SS(53)
	PARAMETER (CF=41.78214D0)
	DATA Z/78.D0/
       DATA A/10.813D0,10.55D0,10.40D0,
     ,10.23D0,10.03D0,9.82D0, 9.60D0, 9.37D0, 9.13D0, 8.882D0,
     ,8.636D0,8.389D0,8.145D0,7.904D0,7.667D0,7.436D0,7.210D0,
     ,6.991D0,6.572D0,6.181D0,5.995D0,5.817D0,5.478D0,5.164D0,
     ,4.873D0,4.603D0,4.475D0,4.352D0,4.120D0,3.905D0,3.704D0,
     ,3.518D0,3.430D0,3.345D0,3.184D0,3.034D0,2.701D0,2.420D0,
     ,2.181D0,1.976D0,1.647D0,1.396D0,1.198D0,1.040D0,0.909D0,
     ,0.801D0,0.709D0,0.632D0,0.565D0,0.508D0,0.459D0,0.416D0,
     ,0.379D0/
C
       DATA SS/0.D0,0.04d0,0.05d0,0.06d0,
     ,0.07d0,0.08d0,0.09d0,0.10d0,0.11d0,0.12d0,0.13d0,0.14d0,
     ,0.15d0,0.16d0,0.17d0,0.18d0,0.19d0,0.20d0,0.22d0,0.24d0,
     ,0.25D0,0.26d0,0.28d0,0.30d0,0.32d0,0.34d0,0.35d0,0.36d0,
     ,0.38d0,0.40d0,0.42d0,0.44d0,0.45d0,0.46d0,0.48d0,0.50d0,
     ,0.55d0,0.60d0,0.65d0,0.70d0,0.80d0,0.90d0,1.00d0,1.10d0,
     ,1.20d0,1.30d0,1.40d0,1.50d0,1.60d0,1.70d0,1.80d0,1.90d0,
     ,2.00d0/
        S02=Z/(CF*A(53))-SS(53)**2
	IF(S.GT.SS(1)) GO TO 1
	FPT=A(1)
	RETURN
 1	IF(S.LT.SS(53)) GO TO 2
C for large arguments the screened Coulomb law is assumed
	FPT=(Z/CF)/(S**2+S02)
	RETURN
 2	IF(S.GT.SS(2)) GO TO 21	
        FPT=A(1)-(A(1)-A(2))*S*S/SS(2)**2
        RETURN
 21     N=1
	GO TO 4
 3	N=N+1
 4	IF(S.LE.SS(N)) GO TO 5
	GO TO 3
 5	FPT=A(N-1)+(S-SS(N-1))*(A(N)-A(N-1))/(SS(N)-SS(N-1))
	RETURN
	END
	


	REAL*8 FUNCTION FRB(X)
C Program for evaluating the scattering amplitude for RUBIDIUM atom, inter-
C polation of data from International Tables for X-ray Crystallography, 1974
C argument X equals to q/4pi within the notations of Landau & Lifshits
	IMPLICIT REAL*8(A-H,O-Z)
	DIMENSION A(56),XX(56)
	PARAMETER (CF=41.78214D0)
	DATA Z/37.D0/
       DATA A/11.778D0,11.699D0,11.460D0,11.088D0,10.613D0,10.073D0,
     ,9.504D0,8.934D0,8.385D0,7.872D0,7.402D0,6.976D0,6.593D0,
     ,6.248D0,5.938D0,5.658D0,5.403D0,5.170D0,4.954D0,4.754D0,
     ,4.566D0,4.224D0,3.916D0,3.773D0,3.636D0,3.382D0,3.149D0,
     ,2.936D0,2.742D0,2.651D0,2.564D0,2.402D0,2.254D0,2.119D0,
     ,1.995D0,1.938D0,1.883D0,1.780D0,1.686D0,1.483D0,1.319D0,
     ,1.182d0,1.068D0,0.887D0,0.749D0,0.640D0,0.551D0,0.478D0,
     ,0.417D0,0.365D0,0.325D0,0.290D0,0.257D0,0.233D0,0.208D0,
     ,0.188D0/
C
       DATA XX/0.D0,0.01D0,0.02D0,0.03D0,0.04d0,0.05d0,0.06d0,
     ,0.07d0,0.08d0,0.09d0,0.10d0,0.11d0,0.12d0,0.13d0,0.14d0,
     ,0.15d0,0.16d0,0.17d0,0.18d0,0.19d0,0.20d0,0.22d0,0.24d0,
     ,0.25D0,0.26d0,0.28d0,0.30d0,0.32d0,0.34d0,0.35d0,0.36d0,
     ,0.38d0,0.40d0,0.42d0,0.44d0,0.45d0,0.46d0,0.48d0,0.50d0,
     ,0.55d0,0.60d0,0.65d0,0.70d0,0.80d0,0.90d0,1.00d0,1.10d0,
     ,1.20d0,1.30d0,1.40d0,1.50d0,1.60d0,1.70d0,1.80d0,1.90d0,
     ,2.00d0/
         X02=Z/(CF*A(56))-XX(56)**2
	IF(X.GT.XX(1)) GO TO 1
	FRB=A(1)
	RETURN
 1	IF(X.LT.XX(56)) GO TO 2
C for large arguments the screened Coulomb law is assumed
	FRB=(Z/CF)/(X**2+X02)
	RETURN
 2	IF(X.GT.XX(2)) GO TO 21
 	FRB=A(1)-(A(1)-A(2))*X*X/XX(2)**2
 	RETURN
 21	N=1
	GO TO 4
 3	N=N+1
 4	IF(X.LE.XX(N)) GO TO 5
	GO TO 3
 5	FRB=A(N-1)+(X-XX(N-1))*(A(N)-A(N-1))/(XX(N)-XX(N-1))
	RETURN
	END
c **************************************************************************

	REAL*8 FUNCTION FRE(S)
C Program for evaluation of scattering amplitude for a HAFNIUM atom, inter-
C polation of data from International Tables for X-ray Crystallography, 1974
C argument S equals to q/(4*pi) within the notation of Landau & Lifshits
	IMPLICIT REAL*8(A-H,O-Z)
	DIMENSION A(53),SS(53)
	PARAMETER (CF=41.78214D0)
	DATA Z/75.D0/
       DATA A/12.263d0,11.83D0,11.60D0,
     ,11.34D0,11.04D0,10.73D0,10.40D0,10.05D0,9.71D0,9.366D0,
     ,9.028D0,8.697D0,8.376D0,8.067D0,7.769D0,7.485D0,7.213D0,
     ,6.954D0,6.475D0,6.043D0,5.843D0,5.653D0,5.301D0,4.981D0,
     ,4.691D0,4.425D0,4.301D0,4.182D0,3.959D0,3.753D0,3.563D0,
     ,3.387D0,3.304D0,3.224D0,3.072D0,2.930D0,2.616D0,2.349D0,
     ,2.122D0,1.926D0,1.608D0,1.363D0,1.169D0,1.012D0,0.883D0,
     ,0.776D0,0.685D0,0.609D0,0.544D0,0.489D0,0.441D0,0.400D0,
     ,0.365D0/
C
       DATA SS/0.D0,0.04d0,0.05d0,0.06d0,
     ,0.07d0,0.08d0,0.09d0,0.10d0,0.11d0,0.12d0,0.13d0,0.14d0,
     ,0.15d0,0.16d0,0.17d0,0.18d0,0.19d0,0.20d0,0.22d0,0.24d0,
     ,0.25D0,0.26d0,0.28d0,0.30d0,0.32d0,0.34d0,0.35d0,0.36d0,
     ,0.38d0,0.40d0,0.42d0,0.44d0,0.45d0,0.46d0,0.48d0,0.50d0,
     ,0.55d0,0.60d0,0.65d0,0.70d0,0.80d0,0.90d0,1.00d0,1.10d0,
     ,1.20d0,1.30d0,1.40d0,1.50d0,1.60d0,1.70d0,1.80d0,1.90d0,
     ,2.00d0/
        S02=Z/(CF*A(53))-SS(53)**2
	IF(S.GT.SS(1)) GO TO 1
	FRE=A(1)
	RETURN
 1	IF(S.LT.SS(53)) GO TO 2
C for large arguments the screened Coulomb law is assumed
	FRE=(Z/CF)/(S**2+S02)
	RETURN
 2	IF(S.GT.SS(2)) GO TO 21	
        FRE=A(1)-(A(1)-A(2))*S*S/SS(2)**2
        RETURN
 21     N=1
	GO TO 4
 3	N=N+1
 4	IF(S.LE.SS(N)) GO TO 5
	GO TO 3
 5	FRE=A(N-1)+(S-SS(N-1))*(A(N)-A(N-1))/(SS(N)-SS(N-1))
	RETURN
	END
	


c **************************************************************************

	REAL*8 FUNCTION FRH(S)
C Program for evaluation of scattering amplitude for a RH atom, inter-
C polation of data from International Tables for X-ray Crystallography, 1974
C argument S equals to q/(4*pi) within the notation of Landau & Lifshits
	IMPLICIT REAL*8(A-H,O-Z)
	DIMENSION A(52),SS(52)
	PARAMETER (CF=41.78214D0)
	DATA Z/45.D0/
       DATA A/9.242D0,8.90D0,8.73D0,
     ,8.53D0,8.31D0, 7.83D0, 7.58D0, 7.33D0, 7.079D0,
     ,6.836d0,6.598d0,6.366d0,6.143d0,5.929d0,5.722d0,5.524d0,
     ,5.334D0,4.976D0,4.648D0,4.493D0,4.345D0,4.066D0,3.809D0,
     ,3.572D0,3.353D0,3.249D0,3.150D0,2.962D0,2.788D0,2.626D0,
     ,2.477D0,2.406D0,2.338D0,2.210D0,2.090D0,1.828D0,1.609D0,
     ,1.426D0,1.273D0,1.035D0,0.861D0,0.731D0,0.631D0,0.551D0,
     ,0.485D0,0.431D0,0.384D0,0.345D0,0.310D0,0.281D0,0.255D0,
     ,0.232D0/
C
       DATA SS/0.D0,0.04d0,0.05d0,0.06d0,
     ,0.07d0,0.09d0,0.10d0,0.11d0,0.12d0,0.13d0,0.14d0,
     ,0.15d0,0.16d0,0.17d0,0.18d0,0.19d0,0.20d0,0.22d0,0.24d0,
     ,0.25D0,0.26d0,0.28d0,0.30d0,0.32d0,0.34d0,0.35d0,0.36d0,
     ,0.38d0,0.40d0,0.42d0,0.44d0,0.45d0,0.46d0,0.48d0,0.50d0,
     ,0.55d0,0.60d0,0.65d0,0.70d0,0.80d0,0.90d0,1.00d0,1.10d0,
     ,1.20d0,1.30d0,1.40d0,1.50d0,1.60d0,1.70d0,1.80d0,1.90d0,
     ,2.00d0/
        S02=Z/(CF*A(52))-SS(52)**2
	IF(S.GT.SS(1)) GO TO 1
	FRH=A(1)
	RETURN
 1	IF(S.LT.SS(52)) GO TO 2
C for large arguments the screened Coulomb law is assumed
	FRH=(Z/CF)/(S**2+S02)
	RETURN
 2	IF(S.GT.SS(2)) GO TO 21	
        FRH=A(1)-(A(1)-A(2))*S*S/SS(2)**2
        RETURN
 21     N=1
	GO TO 4
 3	N=N+1
 4	IF(S.LE.SS(N)) GO TO 5
	GO TO 3
 5	FRH=A(N-1)+(S-SS(N-1))*(A(N)-A(N-1))/(SS(N)-SS(N-1))
	RETURN
	END
	


c **************************************************************************

	REAL*8 FUNCTION FRU(S)
C Program for evaluation of scattering amplitude for a RUTENIUM atom, inter-
C polation of data from International Tables for X-ray Crystallography, 1974
C argument S equals to q/(4*pi) within the notation of Landau & Lifshits
	IMPLICIT REAL*8(A-H,O-Z)
	DIMENSION A(53),SS(53)
	PARAMETER (CF=41.78214D0)
	DATA Z/44.D0/
       DATA A/9.558D0,9.18D0,8.99D0,
     ,8.77D0,8.53D0,8.27D0,8.00D0,7.73D0,7.46D0,7.190D0,
     ,6.928d0,6.672d0,6.426d0,6.188d0,5.960d0,5.741d0,
     ,5.533D0,5.332D0,4.959D0,4.618D0,4.459D0,4.306D0,4.021D0,
     ,3.759D0,3.518D0,3.296D0,3.192D0,3.092D0,2.904D0,2.730D0,
     ,2.570D0,2.421D0,2.351D0,2.284D0,2.157D0,2.040D0,1.782D0,
     ,1.569D0,1.391D0,1.243D0,1.013D0,0.845D0,0.719D0,0.621D0,
     ,0.542D0,0.478D0,0.423D0,0.378D0,0.338D0,0.304D0,0.275D0,
     ,0.249D0,0.227D0/
C
       DATA SS/0.D0,0.04d0,0.05d0,0.06d0,
     ,0.07d0,0.08d0,0.09d0,0.10d0,0.11d0,0.12d0,0.13d0,0.14d0,
     ,0.15d0,0.16d0,0.17d0,0.18d0,0.19d0,0.20d0,0.22d0,0.24d0,
     ,0.25D0,0.26d0,0.28d0,0.30d0,0.32d0,0.34d0,0.35d0,0.36d0,
     ,0.38d0,0.40d0,0.42d0,0.44d0,0.45d0,0.46d0,0.48d0,0.50d0,
     ,0.55d0,0.60d0,0.65d0,0.70d0,0.80d0,0.90d0,1.00d0,1.10d0,
     ,1.20d0,1.30d0,1.40d0,1.50d0,1.60d0,1.70d0,1.80d0,1.90d0,
     ,2.00d0/
        S02=Z/(CF*A(53))-SS(53)**2
	IF(S.GT.SS(1)) GO TO 1
	FRU=A(1)
	RETURN
 1	IF(S.LT.SS(53)) GO TO 2
C for large arguments the screened Coulomb law is assumed
	FRU=(Z/CF)/(S**2+S02)
	RETURN
 2	IF(S.GT.SS(2)) GO TO 21	
        FRU=A(1)-(A(1)-A(2))*S*S/SS(2)**2
        RETURN
 21     N=1
	GO TO 4
 3	N=N+1
 4	IF(S.LE.SS(N)) GO TO 5
	GO TO 3
 5	FRU=A(N-1)+(S-SS(N-1))*(A(N)-A(N-1))/(SS(N)-SS(N-1))
	RETURN
	END
	


	REAL*8 FUNCTION FSI(X)
C Program for evaluating the scattering amplitude for  SILICON atom, inter-
C polation of data from International Tables for X-ray Crystallography, 1974
C argument X equals to q/4pi within the notations of Landau & Lifshits
	IMPLICIT REAL*8(A-H,O-Z)
	DIMENSION A(56),XX(56)
	PARAMETER (CF=41.78214D0)
	DATA Z/14.D0/
       DATA A/5.828D0,5.810D0,5.759D0,5.675D0,5.561D0,5.421D0,
     ,5.258D0,5.077D0,4.882D0,4.677D0,4.467D0,4.255D0,4.043D0,
     ,3.835D0,3.632D0,3.437D0,3.249D0,3.070D0,2.900D0,2.740D0,
     ,2.589D0,2.315D0,2.076D0,1.969D0,1.869D0,1.689D0,1.534D0,
     ,1.400D0,1.284D0,1.231D0,1.182D0,1.094D0,1.017D0,0.949D0,
     ,0.888D0,0.861D0,0.834D0,0.786D0,0.743D0,0.651D0,0.578D0,
     ,0.517d0,0.465D0,0.383D0,0.320D0,0.270D0,0.231D0,0.198D0,
     ,0.172D0,0.150D0,0.132D0,0.117D0,0.104D0,0.093D0,0.084D0,
     ,0.076D0/
C
       DATA XX/0.D0,0.01D0,0.02D0,0.03D0,0.04d0,0.05d0,0.06d0,
     ,0.07d0,0.08d0,0.09d0,0.10d0,0.11d0,0.12d0,0.13d0,0.14d0,
     ,0.15d0,0.16d0,0.17d0,0.18d0,0.19d0,0.20d0,0.22d0,0.24d0,
     ,0.25D0,0.26d0,0.28d0,0.30d0,0.32d0,0.34d0,0.35d0,0.36d0,
     ,0.38d0,0.40d0,0.42d0,0.44d0,0.45d0,0.46d0,0.48d0,0.50d0,
     ,0.55d0,0.60d0,0.65d0,0.70d0,0.80d0,0.90d0,1.00d0,1.10d0,
     ,1.20d0,1.30d0,1.40d0,1.50d0,1.60d0,1.70d0,1.80d0,1.90d0,
     ,2.00d0/
         X02=Z/(CF*A(56))-XX(56)**2
	IF(X.GT.XX(1)) GO TO 1
	FSI=A(1)
	RETURN
 1	IF(X.LT.XX(56)) GO TO 2
C for large arguments the screened Coulomb law is assumed
	FSI=(Z/CF)/(X**2+X02)
	RETURN
 2	IF(X.GT.XX(2)) GO TO 21
 	FSI=A(1)-(A(1)-A(2))*X*X/XX(2)**2
 	RETURN
 21	N=1
	GO TO 4
 3	N=N+1
 4	IF(X.LE.XX(N)) GO TO 5
	GO TO 3
 5	FSI=A(N-1)+(X-XX(N-1))*(A(N)-A(N-1))/(XX(N)-XX(N-1))
	RETURN
	END

	REAL*8 FUNCTION FSN(X)
C Program for evaluating the scattering amplitude for a SN atom, inter-
C polation of data from International Tables for X-ray Crystallography, 1974
C argument X equals to q/4pi within the notations of Landau & Lifshits
	IMPLICIT REAL*8(A-H,O-Z)
	DIMENSION A(56),XX(56)
	PARAMETER (CF=41.78214D0)
	DATA Z/50.D0/
       DATA A/10.859D0,10.833D0,10.750D0,10.615D0,10.433D0,10.209D0,
     ,9.950D0,9.664D0,9.357D0,9.037D0,8.709D0,8.380D0,8.053D0,
     ,7.732D0,7.419D0,7.118D0,6.829D0,6.552D0,6.289D0,6.039D0,
     ,5.803D0,5.368D0,4.979D0,4.801D0,4.633D0,4.323D0,4.044D0,
     ,3.792D0,3.564D0,3.458D0,3.356D0,3.165D0,2.990D0,2.828D0,
     ,2.678D0,2.608D0,2.539D0,2.409D0,2.288D0,2.019D0,1.790D0,
     ,1.574d0,1.426D0,1.157D0,0.956D0,0.805D0,0.688D0,0.597D0,
     ,0.525D0,0.465D0,0.416D0,0.374D0,0.337D0,0.307D0,0.279D0,
     ,0.255D0/
C
       DATA XX/0.D0,0.01D0,0.02D0,0.03D0,0.04d0,0.05d0,0.06d0,
     ,0.07d0,0.08d0,0.09d0,0.10d0,0.11d0,0.12d0,0.13d0,0.14d0,
     ,0.15d0,0.16d0,0.17d0,0.18d0,0.19d0,0.20d0,0.22d0,0.24d0,
     ,0.25D0,0.26d0,0.28d0,0.30d0,0.32d0,0.34d0,0.35d0,0.36d0,
     ,0.38d0,0.40d0,0.42d0,0.44d0,0.45d0,0.46d0,0.48d0,0.50d0,
     ,0.55d0,0.60d0,0.65d0,0.70d0,0.80d0,0.90d0,1.00d0,1.10d0,
     ,1.20d0,1.30d0,1.40d0,1.50d0,1.60d0,1.70d0,1.80d0,1.90d0,
     ,2.00d0/
         X02=Z/(CF*A(56))-XX(56)**2
	IF(X.GT.XX(1)) GO TO 1
	FSN=A(1)
	RETURN
 1	IF(X.LT.XX(56)) GO TO 2
C for large arguments the screened Coulomb law is assumed
	FSN=(Z/CF)/(X**2+X02)
	RETURN
 2	IF(X.GT.XX(2)) GO TO 21
 	FSN=A(1)-(A(1)-A(2))*X*X/XX(2)**2
 	RETURN
 21	N=1
	GO TO 4
 3	N=N+1
 4	IF(X.LE.XX(N)) GO TO 5
	GO TO 3
 5	FSN=A(N-1)+(X-XX(N-1))*(A(N)-A(N-1))/(XX(N)-XX(N-1))
	RETURN
	END

	REAL*8 FUNCTION FSR(X)
C Program for evaluating the scattering amplitude for a STRONCIUM atom, inter-
C polation of data from International Tables for X-ray Crystallography, 1974
C argument X equals to q/4pi within the notations of Landau & Lifshits
	IMPLICIT REAL*8(A-H,O-Z)
	DIMENSION A(56),XX(56)
	PARAMETER (CF=41.78214D0)
	DATA Z/38.D0/
       DATA A/13.109D0,13.035D0,12.816D0,12.468D0,12.013D0,11.476D0,
     ,10.888D0,10.273D0,9.655D0,9.052D0,8.478D0,7.940D0,7.443D0,
     ,6.988D0,6.575D0,6.200D0,5.862D0,5.555D0,5.278D0,5.025D0,
     ,4.794D0,4.387D0,4.039D0,3.882D0,3.735D0,3.465D0,3.224D0,
     ,3.007D0,2.810D0,2.718D0,2.630D0,2.466D0,2.315D0,2.178D0,
     ,2.052D0,1.993D0,1.936D0,1.830D0,1.733D0,1.522D0,1.350D0,
     ,1.208d0,1.089D0,0.902D0,0.762D0,0.651D0,0.562D0,0.488D0,
     ,0.427D0,0.375D0,0.333D0,0.297D0,0.264D0,0.239D0,0.214D0,
     ,0.194D0/
C
       DATA XX/0.D0,0.01D0,0.02D0,0.03D0,0.04d0,0.05d0,0.06d0,
     ,0.07d0,0.08d0,0.09d0,0.10d0,0.11d0,0.12d0,0.13d0,0.14d0,
     ,0.15d0,0.16d0,0.17d0,0.18d0,0.19d0,0.20d0,0.22d0,0.24d0,
     ,0.25D0,0.26d0,0.28d0,0.30d0,0.32d0,0.34d0,0.35d0,0.36d0,
     ,0.38d0,0.40d0,0.42d0,0.44d0,0.45d0,0.46d0,0.48d0,0.50d0,
     ,0.55d0,0.60d0,0.65d0,0.70d0,0.80d0,0.90d0,1.00d0,1.10d0,
     ,1.20d0,1.30d0,1.40d0,1.50d0,1.60d0,1.70d0,1.80d0,1.90d0,
     ,2.00d0/
         X02=Z/(CF*A(56))-XX(56)**2
	IF(X.GT.XX(1)) GO TO 1
	FSR=A(1)
	RETURN
 1	IF(X.LT.XX(56)) GO TO 2
C for large arguments the screened Coulomb law is assumed
	FSR=(Z/CF)/(X**2+X02)
	RETURN
 2	IF(X.GT.XX(2)) GO TO 21
 	FSR=A(1)-(A(1)-A(2))*X*X/XX(2)**2
 	RETURN
 21	N=1
	GO TO 4
 3	N=N+1
 4	IF(X.LE.XX(N)) GO TO 5
	GO TO 3
 5	FSR=A(N-1)+(X-XX(N-1))*(A(N)-A(N-1))/(XX(N)-XX(N-1))
	RETURN
	END
c **************************************************************************

	REAL*8 FUNCTION FTA(S)
C Program for evaluation of scattering amplitude for a TALLIUM atom, inter-
C polation of data from International Tables for X-ray Crystallography, 1974
C argument S equals to q/(4*pi) within the notation of Landau & Lifshits
	IMPLICIT REAL*8(A-H,O-Z)
	DIMENSION A(53),SS(53)
	PARAMETER (CF=41.78214D0)
	DATA Z/73.D0/
       DATA A/12.856D0,12.31D0,12.01D0,
     ,11.69D0,11.33D0,10.95D0,10.55D0,10.15D0,9.75D0,9.363D0,
     ,8.982D0,8.616D0,8.266D0,7.933D0,7.617D0,7.321D0,7.040D0,
     ,6.776D0,6.295D0,5.867D0,5.672D0,5.487D0,5.147D0,4.840D0,
     ,4.563D0,4.310D0,4.191D0,4.078D0,3.865D0,3.668D0,3.486D0,
     ,3.317D0,3.237D0,3.159D0,3.013D0,2.876D0,2.571D0,2.311D0,
     ,2.089D0,1.896D0,1.583D0,1.341D0,1.148D0,0.993D0,0.864D0,
     ,0.758D0,0.668D0,0.593D0,0.530D0,0.475D0,0.429D0,0.389D0,
     ,0.355D0/
C
       DATA SS/0.D0,0.04d0,0.05d0,0.06d0,
     ,0.07d0,0.08d0,0.09d0,0.10d0,0.11d0,0.12d0,0.13d0,0.14d0,
     ,0.15d0,0.16d0,0.17d0,0.18d0,0.19d0,0.20d0,0.22d0,0.24d0,
     ,0.25D0,0.26d0,0.28d0,0.30d0,0.32d0,0.34d0,0.35d0,0.36d0,
     ,0.38d0,0.40d0,0.42d0,0.44d0,0.45d0,0.46d0,0.48d0,0.50d0,
     ,0.55d0,0.60d0,0.65d0,0.70d0,0.80d0,0.90d0,1.00d0,1.10d0,
     ,1.20d0,1.30d0,1.40d0,1.50d0,1.60d0,1.70d0,1.80d0,1.90d0,
     ,2.00d0/
        S02=Z/(CF*A(53))-SS(53)**2
	IF(S.GT.SS(1)) GO TO 1
	FTA=A(1)
	RETURN
 1	IF(S.LT.SS(53)) GO TO 2
C for large arguments the screened Coulomb law is assumed
	FTA=(Z/CF)/(S**2+S02)
	RETURN
 2	IF(S.GT.SS(2)) GO TO 21	
        FTA=A(1)-(A(1)-A(2))*S*S/SS(2)**2
        RETURN
 21     N=1
	GO TO 4
 3	N=N+1
 4	IF(S.LE.SS(N)) GO TO 5
	GO TO 3
 5	FTA=A(N-1)+(S-SS(N-1))*(A(N)-A(N-1))/(SS(N)-SS(N-1))
	RETURN
	END
	


c **************************************************************************

	REAL*8 FUNCTION FTH(S)
C Program for evaluation of scattering amplitude for a THORIUM atom, inter-
C polation of data from International Tables for X-ray Crystallography, 1974
C argument S equals to q/(4*pi) within the notation of Landau & Lifshits
	IMPLICIT REAL*8(A-H,O-Z)
	DIMENSION A(53),SS(53)
	PARAMETER (CF=41.78214D0)
	DATA Z/90.D0/
       DATA A/20.115D0,18.92D0,18.33D0,
     ,17.66D0,16.93D0,16.19D0,15.43D0,14.68D0,13.95D0,13.255D0,
     ,12.594D0,11.972D0,11.388D0,10.845D0,10.339D0,9.868D0,9.430D0,
     ,9.022D0,8.287D0,7.645D0,7.353D0,7.079D0,6.578D0,6.129D0,
     ,5.727D0,5.364D0,5.196D0,5.036D0,4.768D0,4.466D0,4.218D0,
     ,3.992D0,3.885D0,3.784D0,3.592D0,3.416D0,3.032D0,2.712D0,
     ,2.442D0,2.212D0,1.840D0,1.554D0,1.330D0,1.152D0,1.007D0,
     ,0.889D0,0.791D0,0.708D0,0.637D0,0.576D0,0.523D0,0.477D0,
     ,0.436D0/
C
       DATA SS/0.D0,0.04d0,0.05d0,0.06d0,
     ,0.07d0,0.08d0,0.09d0,0.10d0,0.11d0,0.12d0,0.13d0,0.14d0,
     ,0.15d0,0.16d0,0.17d0,0.18d0,0.19d0,0.20d0,0.22d0,0.24d0,
     ,0.25D0,0.26d0,0.28d0,0.30d0,0.32d0,0.34d0,0.35d0,0.36d0,
     ,0.38d0,0.40d0,0.42d0,0.44d0,0.45d0,0.46d0,0.48d0,0.50d0,
     ,0.55d0,0.60d0,0.65d0,0.70d0,0.80d0,0.90d0,1.00d0,1.10d0,
     ,1.20d0,1.30d0,1.40d0,1.50d0,1.60d0,1.70d0,1.80d0,1.90d0,
     ,2.00d0/
        S02=Z/(CF*A(53))-SS(53)**2
	IF(S.GT.SS(1)) GO TO 1
	FTH=A(1)
	RETURN
 1	IF(S.LT.SS(53)) GO TO 2
C for large arguments the screened Coulomb law is assumed
	FTH=(Z/CF)/(S**2+S02)
	RETURN
 2	IF(S.GT.SS(2)) GO TO 21	
        FTH=A(1)-(A(1)-A(2))*S*S/SS(2)**2
        RETURN
 21     N=1
	GO TO 4
 3	N=N+1
 4	IF(S.LE.SS(N)) GO TO 5
	GO TO 3
 5	FTH=A(N-1)+(S-SS(N-1))*(A(N)-A(N-1))/(SS(N)-SS(N-1))
	RETURN
	END
	


	REAL*8 FUNCTION FTI(X)
C Program for evaluating the scattering amplitude for a TITANIUM atom, inter-
C polation of data from International Tables for X-ray Crystallography, 1974
C argument X equals to q/4pi within the notations of Landau & Lifshits
	IMPLICIT REAL*8(A-H,O-Z)
	DIMENSION A(56),XX(56)
	PARAMETER (CF=41.78214D0)
	DATA Z/22.D0/
       DATA A/8.776D0,8.740D0,8.631D0,8.455D0,8.220D0,7.937D0,
     ,7.618D0,7.274D0,6.917D0,6.556D0,6.199D0,5.853D0,5.522D0,
     ,5.209D0,4.916D0,4.643D0,4.390D0,4.157D0,3.943D0,3.745D0,
     ,3.564D0,3.242D0,2.967D0,2.844D0,2.730D0,2.523D0,2.341D0,
     ,2.178D0,2.032D0,1.964D0,1.899D0,1.778D0,1.668D0,1.566D0,
     ,1.472D0,1.428D0,1.385D0,1.305D0,1.230D0,1.067D0,0.930D0,
     ,0.816d0,0.721D0,0.573D0,0.467D0,0.389D0,0.330D0,0.285D0,
     ,0.249D0,0.219D0,0.195D0,0.175D0,0.157D0,0.143D0,0.129D0,
     ,0.117D0/
C
       DATA XX/0.D0,0.01D0,0.02D0,0.03D0,0.04d0,0.05d0,0.06d0,
     ,0.07d0,0.08d0,0.09d0,0.10d0,0.11d0,0.12d0,0.13d0,0.14d0,
     ,0.15d0,0.16d0,0.17d0,0.18d0,0.19d0,0.20d0,0.22d0,0.24d0,
     ,0.25D0,0.26d0,0.28d0,0.30d0,0.32d0,0.34d0,0.35d0,0.36d0,
     ,0.38d0,0.40d0,0.42d0,0.44d0,0.45d0,0.46d0,0.48d0,0.50d0,
     ,0.55d0,0.60d0,0.65d0,0.70d0,0.80d0,0.90d0,1.00d0,1.10d0,
     ,1.20d0,1.30d0,1.40d0,1.50d0,1.60d0,1.70d0,1.80d0,1.90d0,
     ,2.00d0/
         X02=Z/(CF*A(56))-XX(56)**2
	IF(X.GT.XX(1)) GO TO 1
	FTI=A(1)
	RETURN
 1	IF(X.LT.XX(56)) GO TO 2
C for large arguments the screened Coulomb law is assumed
	FTI=(Z/CF)/(X**2+X02)
	RETURN
 2	IF(X.GT.XX(2)) GO TO 21
 	FTI=A(1)-(A(1)-A(2))*X*X/XX(2)**2
 	RETURN
 21	N=1
	GO TO 4
 3	N=N+1
 4	IF(X.LE.XX(N)) GO TO 5
	GO TO 3
 5	FTI=A(N-1)+(X-XX(N-1))*(A(N)-A(N-1))/(XX(N)-XX(N-1))
	RETURN
	END
c **************************************************************************

	REAL*8 FUNCTION FTL(S)
C Program for evaluation of scattering amplitude for a TALLIUM atom, inter-
C polation of data from International Tables for X-ray Crystallography, 1974
C argument S equals to q/(4*pi) within the notation of Landau & Lifshits
	IMPLICIT REAL*8(A-H,O-Z)
	DIMENSION A(53),SS(53)
	PARAMETER (CF=41.78214D0)
	DATA Z/81.D0/
       DATA A/12.109D0,11.71D0,11.51D0,
     ,11.27D0,11.00D0,10.72D0,10.42D0,10.12D0, 9.81D0,9.500D0,
     ,9.195D0,8.896D0,8.603D0,8.320D0,8.046D0,7.781D0,7.526D0,
     ,7.282D0,6.822D0,6.399D0,6.201D0,6.011D0,5.654D0,5.327D0,
     ,5.025D0,4.746D0,4.614D0,4.488D0,4.249D0,4.028D0,3.823D0,
     ,3.632D0,3.541D0,3.454D0,3.288D0,3.133D0,2.789D0,2.497D0,
     ,2.248D0,2.035D0,1.692D0,1.431D0,1.228D0,1.066D0,0.934D0,
     ,0.824D0,0.731D0,0.653D0,0.585D0,0.527D0,0.476D0,0.432D0,
     ,0.394D0/
C
       DATA SS/0.D0,0.04d0,0.05d0,0.06d0,
     ,0.07d0,0.08d0,0.09d0,0.10d0,0.11d0,0.12d0,0.13d0,0.14d0,
     ,0.15d0,0.16d0,0.17d0,0.18d0,0.19d0,0.20d0,0.22d0,0.24d0,
     ,0.25D0,0.26d0,0.28d0,0.30d0,0.32d0,0.34d0,0.35d0,0.36d0,
     ,0.38d0,0.40d0,0.42d0,0.44d0,0.45d0,0.46d0,0.48d0,0.50d0,
     ,0.55d0,0.60d0,0.65d0,0.70d0,0.80d0,0.90d0,1.00d0,1.10d0,
     ,1.20d0,1.30d0,1.40d0,1.50d0,1.60d0,1.70d0,1.80d0,1.90d0,
     ,2.00d0/
        S02=Z/(CF*A(53))-SS(53)**2
	IF(S.GT.SS(1)) GO TO 1
	FTL=A(1)
	RETURN
 1	IF(S.LT.SS(53)) GO TO 2
C for large arguments the screened Coulomb law is assumed
	FTL=(Z/CF)/(S**2+S02)
	RETURN
 2	IF(S.GT.SS(2)) GO TO 21	
        FTL=A(1)-(A(1)-A(2))*S*S/SS(2)**2
        RETURN
 21     N=1
	GO TO 4
 3	N=N+1
 4	IF(S.LE.SS(N)) GO TO 5
	GO TO 3
 5	FTL=A(N-1)+(S-SS(N-1))*(A(N)-A(N-1))/(SS(N)-SS(N-1))
	RETURN
	END
	


	REAL*8 FUNCTION FV(X)
C Program for evaluating the scattering amplitude for VANADIUM atom, inter-
C polation of data from International Tables for X-ray Crystallography, 1974
C argument X equals to q/4pi within the notations of Landau & Lifshits
	IMPLICIT REAL*8(A-H,O-Z)
	DIMENSION A(56),XX(56)
	PARAMETER (CF=41.78214D0)
	DATA Z/23.D0/
       DATA A/8.305D0,8.274D0,8.180D0,8.029D0,7.826D0,7.581D0,
     ,7.303D0,7.002D0,6.686D0,6.365D0,6.045D0,5.732D0,5.430D0,
     ,5.142D0,4.871D0,4.616D0,4.378D0,4.158D0,3.953D0,3.763D0,
     ,3.588D0,3.276D0,3.006D0,2.885D0,2.772D0,2.568D0,2.386D0,
     ,2.225D0,2.079D0,2.011D0,1.947D0,1.826D0,1.716D0,1.614D0,
     ,1.520D0,1.476D0,1.433D0,1.352D0,1.277D0,1.111D0,0.973D0,
     ,0.856d0,0.757D0,0.602D0,0.490D0,0.408D0,0.345D0,0.297D0,
     ,0.259D0,0.228D0,0.203D0,0.182D0,0.163D0,0.148D0,0.134D0,
     ,0.122D0/
C
       DATA XX/0.D0,0.01D0,0.02D0,0.03D0,0.04d0,0.05d0,0.06d0,
     ,0.07d0,0.08d0,0.09d0,0.10d0,0.11d0,0.12d0,0.13d0,0.14d0,
     ,0.15d0,0.16d0,0.17d0,0.18d0,0.19d0,0.20d0,0.22d0,0.24d0,
     ,0.25D0,0.26d0,0.28d0,0.30d0,0.32d0,0.34d0,0.35d0,0.36d0,
     ,0.38d0,0.40d0,0.42d0,0.44d0,0.45d0,0.46d0,0.48d0,0.50d0,
     ,0.55d0,0.60d0,0.65d0,0.70d0,0.80d0,0.90d0,1.00d0,1.10d0,
     ,1.20d0,1.30d0,1.40d0,1.50d0,1.60d0,1.70d0,1.80d0,1.90d0,
     ,2.00d0/
         X02=Z/(CF*A(56))-XX(56)**2
	IF(X.GT.XX(1)) GO TO 1
	FV=A(1)
	RETURN
 1	IF(X.LT.XX(56)) GO TO 2
C for large arguments the screened Coulomb law is assumed
	FV=(Z/CF)/(X**2+X02)
	RETURN
 2	IF(X.GT.XX(2)) GO TO 21
 	FV=A(1)-(A(1)-A(2))*X*X/XX(2)**2
 	RETURN
 21	N=1
	GO TO 4
 3	N=N+1
 4	IF(X.LE.XX(N)) GO TO 5
	GO TO 3
 5	FV=A(N-1)+(X-XX(N-1))*(A(N)-A(N-1))/(XX(N)-XX(N-1))
	RETURN
	END
c **************************************************************************

	REAL*8 FUNCTION FW(S)
C Program for evaluation of scattering amplitude for a TUNGSTEN atom, inter-
C polation of data from International Tables for X-ray Crystallography, 1974
C argument S equals to q/(4*pi) within the notation of Landau & Lifshits
	IMPLICIT REAL*8(A-H,O-Z)
	DIMENSION A(53),SS(53)
	PARAMETER (CF=41.78214D0)
	DATA Z/74.D0/
       DATA A/12.543D0,12.06D0,11.80D0,
     ,11.51D0,11.18D0,10.83D0,10.47D0,10.10D0,9.74D0,9.369D0,
     ,9.011D0,8.663D0,8.327D0,8.006D0,7.699D0,7.408D0,
     ,7.132D0,6.870D0,6.388D0,5.957D0,5.759D0,5.571D0,5.224D0,
     ,4.910D0,4.626D0,4.366D0,4.245D0,4.129D0,3.910D0,3.709D0,
     ,3.523D0,3.350D0,3.269D0,3.190D0,3.041D0,2.902D0,2.592D0,
     ,2.330D0,2.105D0,1.911D0,1.596D0,1.352D0,1.159D0,1.003D0,
     ,0.874D0,0.767D0,0.677D0,0.601D0,0.537D0,0.482D0,0.435D0,
     ,0.395D0,0.360d0/
C
       DATA SS/0.D0,0.04d0,0.05d0,0.06d0,
     ,0.07d0,0.08d0,0.09d0,0.10d0,0.11d0,0.12d0,0.13d0,0.14d0,
     ,0.15d0,0.16d0,0.17d0,0.18d0,0.19d0,0.20d0,0.22d0,0.24d0,
     ,0.25D0,0.26d0,0.28d0,0.30d0,0.32d0,0.34d0,0.35d0,0.36d0,
     ,0.38d0,0.40d0,0.42d0,0.44d0,0.45d0,0.46d0,0.48d0,0.50d0,
     ,0.55d0,0.60d0,0.65d0,0.70d0,0.80d0,0.90d0,1.00d0,1.10d0,
     ,1.20d0,1.30d0,1.40d0,1.50d0,1.60d0,1.70d0,1.80d0,1.90d0,
     ,2.00d0/
        S02=Z/(CF*A(53))-SS(53)**2
	IF(S.GT.SS(1)) GO TO 1
	FW=A(1)
	RETURN
 1	IF(S.LT.SS(53)) GO TO 2
C for large arguments the screened Coulomb law is assumed
	FW=(Z/CF)/(S**2+S02)
	RETURN
 2	IF(S.GT.SS(2)) GO TO 21	
        FW=A(1)-(A(1)-A(2))*S*S/SS(2)**2
        RETURN
 21     N=1
	GO TO 4
 3	N=N+1
 4	IF(S.LE.SS(N)) GO TO 5
	GO TO 3
 5	FW=A(N-1)+(S-SS(N-1))*(A(N)-A(N-1))/(SS(N)-SS(N-1))
	RETURN
	END
	


c **************************************************************************

	REAL*8 FUNCTION FY(S)
C Program for evaluation of scattering amplitude for a YTTRIUM atom, inter-
C polation of data from International Tables for X-ray Crystallography, 1974
C argument S equals to q/(4*pi) within the notation of Landau & Lifshits
	IMPLICIT REAL*8(A-H,O-Z)
	DIMENSION A(53),SS(53)
	PARAMETER (CF=41.78214D0)
	DATA Z/39.D0/
       DATA A/12.674D0,11.79D0,11.34D0,
     ,10.84D0,10.31D0,9.77D0, 9.23D0, 8.70D0, 8.20D0, 7.722D0,
     ,7.278d0,6.865d0,6.485d0,6.136d0,5.816d0,5.523d0,5.254d0,
     ,5.008D0,4.570D0,4.195D0,4.027D0,3.869D0,3.583D0,3.329D0,
     ,3.101D0,2.895D0,2.799D0,2.708D0,2.538D0,2.383D0,2.241D0,
     ,2.111D0,2.049D0,1.991D0,1.881D0,1.780D0,1.562D0,1.383D0,
     ,1.235D0,1.111D0,0.918D0,0.774D0,0.662D0,0.572D0,0.498D0,
     ,0.436D0,0.384D0,0.341D0,0.303D0,0.272D0,0.244D0,0.221D0,
     ,0.201D0/
C
       DATA SS/0.D0,0.04d0,0.05d0,0.06d0,
     ,0.07d0,0.08d0,0.09d0,0.10d0,0.11d0,0.12d0,0.13d0,0.14d0,
     ,0.15d0,0.16d0,0.17d0,0.18d0,0.19d0,0.20d0,0.22d0,0.24d0,
     ,0.25D0,0.26d0,0.28d0,0.30d0,0.32d0,0.34d0,0.35d0,0.36d0,
     ,0.38d0,0.40d0,0.42d0,0.44d0,0.45d0,0.46d0,0.48d0,0.50d0,
     ,0.55d0,0.60d0,0.65d0,0.70d0,0.80d0,0.90d0,1.00d0,1.10d0,
     ,1.20d0,1.30d0,1.40d0,1.50d0,1.60d0,1.70d0,1.80d0,1.90d0,
     ,2.00d0/
        S02=Z/(CF*A(53))-SS(53)**2
	IF(S.GT.SS(1)) GO TO 1
	FY=A(1)
	RETURN
 1	IF(S.LT.SS(53)) GO TO 2
C for large arguments the screened Coulomb law is assumed
	FY=(Z/CF)/(S**2+S02)
	RETURN
 2	IF(S.GT.SS(2)) GO TO 21	
        FY=A(1)-(A(1)-A(2))*S*S/SS(2)**2
        RETURN
 21     N=1
	GO TO 4
 3	N=N+1
 4	IF(S.LE.SS(N)) GO TO 5
	GO TO 3
 5	FY=A(N-1)+(S-SS(N-1))*(A(N)-A(N-1))/(SS(N)-SS(N-1))
	RETURN
	END
	


	REAL*8 FUNCTION FZN(X)
C Program for evaluating the scattering amplitude for a ZINK atom, inter-
C polation of data from International Tables for X-ray Crystallography, 1974
C argument X equals to q/4pi within the notations of Landau & Lifshits
	IMPLICIT REAL*8(A-H,O-Z)
	DIMENSION A(56),XX(56)
	PARAMETER (CF=41.78214D0)
	DATA Z/30.D0/
       DATA A/6.065D0,6.051D0,6.009D0,5.941D0,5.849D0,5.735D0,
     ,5.603D0,5.457D0,5.299D0,5.133D0,4.962D0,4.790D0,4.618D0,
     ,4.449D0,4.283D0,4.123D0,3.969D0,3.822D0,3.681D0,3.547D0,
     ,3.421D0,3.186D0,2.977D0,2.880D0,2.789D0,2.620D0,2.468D0,
     ,2.329D0,2.203D0,2.144D0,2.087D0,1.980D0,1.882D0,1.790D0,
     ,1.704D0,1.663D0,1.624D0,1.549d0,1.478D0,1.319D0,1.181D0,
     ,1.061d0,0.955D0,0.781D0,0.646D0,0.541D0,0.457D0,0.391D0,
     ,0.339D0,0.296D0,0.261D0,0.232D0,0.208D0,0.188D0,0.170D0,
     ,0.155D0/
C
       DATA XX/0.D0,0.01D0,0.02D0,0.03D0,0.04d0,0.05d0,0.06d0,
     ,0.07d0,0.08d0,0.09d0,0.10d0,0.11d0,0.12d0,0.13d0,0.14d0,
     ,0.15d0,0.16d0,0.17d0,0.18d0,0.19d0,0.20d0,0.22d0,0.24d0,
     ,0.25D0,0.26d0,0.28d0,0.30d0,0.32d0,0.34d0,0.35d0,0.36d0,
     ,0.38d0,0.40d0,0.42d0,0.44d0,0.45d0,0.46d0,0.48d0,0.50d0,
     ,0.55d0,0.60d0,0.65d0,0.70d0,0.80d0,0.90d0,1.00d0,1.10d0,
     ,1.20d0,1.30d0,1.40d0,1.50d0,1.60d0,1.70d0,1.80d0,1.90d0,
     ,2.00d0/
         X02=Z/(CF*A(56))-XX(56)**2
	IF(X.GT.XX(1)) GO TO 1
	FZN=A(1)
	RETURN
 1	IF(X.LT.XX(56)) GO TO 2
C for large arguments the screened Coulomb law is assumed
	FZN=(Z/CF)/(X**2+X02)
	RETURN
 2	IF(X.GT.XX(2)) GO TO 21
 	FZN=A(1)-(A(1)-A(2))*X*X/XX(2)**2
 	RETURN
 21	N=1
	GO TO 4
 3	N=N+1
 4	IF(X.LE.XX(N)) GO TO 5
	GO TO 3
 5	FZN=A(N-1)+(X-XX(N-1))*(A(N)-A(N-1))/(XX(N)-XX(N-1))
	RETURN
	END
c **************************************************************************

	REAL*8 FUNCTION FZR(S)
C Program for evaluation of scattering amplitude for a ZIRCONIUM atom, inter-
C polation of data from International Tables for X-ray Crystallography, 1974
C argument S equals to q/(4*pi) within the notation of Landau & Lifshits
	IMPLICIT REAL*8(A-H,O-Z)
	DIMENSION A(53),SS(53)
	PARAMETER (CF=41.78214D0)
	DATA Z/40.D0/
       DATA A/12.166D0,11.41D0,11.04D0,
     ,10.62D0,10.15D0,9.68D0, 9.20D0, 8.72D0, 8.26D0, 7.818D0,
     ,7.400d0,7.007d0,6.640d0,6.299d0,5.983d0,5.689d0,
     ,5.419D0,5.168D0,4.721D0,4.333D0,4.158D0,3.995D0,3.697D0,
     ,3.433D0,3.196D0,2.982D0,2.883D0,2.789D0,2.613D0,2.452D0,
     ,2.305D0,2.171D0,2.108D0,2.047D0,1.934D0,1.829D0,1.603D0,
     ,1.417D0,1.263D0,1.135D0,0.935D0,0.787D0,0.673D0,0.582D0,
     ,0.507D0,0.445D0,0.393D0,0.349D0,0.311D0,0.278D0,0.251D0,
     ,0.227D0,0.206D0/
C
       DATA SS/0.D0,0.04d0,0.05d0,0.06d0,
     ,0.07d0,0.08d0,0.09d0,0.10d0,0.11d0,0.12d0,0.13d0,0.14d0,
     ,0.15d0,0.16d0,0.17d0,0.18d0,0.19d0,0.20d0,0.22d0,0.24d0,
     ,0.25D0,0.26d0,0.28d0,0.30d0,0.32d0,0.34d0,0.35d0,0.36d0,
     ,0.38d0,0.40d0,0.42d0,0.44d0,0.45d0,0.46d0,0.48d0,0.50d0,
     ,0.55d0,0.60d0,0.65d0,0.70d0,0.80d0,0.90d0,1.00d0,1.10d0,
     ,1.20d0,1.30d0,1.40d0,1.50d0,1.60d0,1.70d0,1.80d0,1.90d0,
     ,2.00d0/
        S02=Z/(CF*A(53))-SS(53)**2
	IF(S.GT.SS(1)) GO TO 1
	FZR=A(1)
	RETURN
 1	IF(S.LT.SS(53)) GO TO 2
C for large arguments the screened Coulomb law is assumed
	FZR=(Z/CF)/(S**2+S02)
	RETURN
 2	IF(S.GT.SS(2)) GO TO 21	
        FZR=A(1)-(A(1)-A(2))*S*S/SS(2)**2
        RETURN
 21     N=1
	GO TO 4
 3	N=N+1
 4	IF(S.LE.SS(N)) GO TO 5
	GO TO 3
 5	FZR=A(N-1)+(S-SS(N-1))*(A(N)-A(N-1))/(SS(N)-SS(N-1))
	RETURN
	END
	


