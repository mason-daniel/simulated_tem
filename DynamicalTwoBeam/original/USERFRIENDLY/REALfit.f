
c============================================================================= 
c                                  REALfit.f
c=============================================================================
c
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


	implicit real*8 (a-h,o-z)
	parameter (ndata=501,ma=10)
	REAL*8 s(ndata),f(ndata),sig(ndata),a(ma)
	REAL*8 covar(ma,ma),alpha(ma,ma)
	INTEGER ia(ma)
	CHARACTER*1 answer
	CHARACTER*60 filename1,filename2
	EXTERNAL funcs

	nca=ma

	write(6,*) 
     *'INPUT filename (=the output filename for inputREAL.exe)='
	read(5,*) filename1
	write(6,*) 'OUTPUT filename (not more than 60 characters)='
	read(5,*) filename2

	open(unit=101,file=filename1)
	open(unit=102,file=filename2)

	do 20 j=1,ndata
	read(101,*) k,s(k),f(k)
 20	continue
	close(101)
	
	do 2 i=1,ndata
	sig(i)=1.d-01
 2	continue

	write(6,*) 'input the initial 5 values of a(i)='
	read(5,*) a(1),a(3),a(5),a(7),a(9)
	write(6,*) 'input the initial 5 values of b(i)='
	read(5,*) a(2),a(4),a(6),a(8),a(10)

	do 3 i=1,ma
	ia(i)=1
 3	continue
	alamda=-1.d0
 4	call mrqmin(s,f,sig,ndata,a,ia,ma,covar,alpha,nca,
     *          chisq,funcs,alamda)
        write(6,*) 'chi^2=',(chisq/float(ndata)),'  alamda=',alamda
        write(6,*) 'continue? (y/n)'
        read(5,*) answer
	if(answer.eq.'y') go to 4
	write(*,*) 'writing to the file ',filename2

        write(102,98) a(1),a(3),a(5),a(7),a(9)
        write(102,*)
        write(102,99) a(2),a(4),a(6),a(8),a(10)
        write(102,*)
        write(102,*) (chisq/float(ndata))

	close(101)
	close(102)
  
  98	format(1x,6Ha(i)= ,5(e12.5,1x))
  99	format(1x,6Hb(i)= ,5(e12.5,1x))
  
 	STOP
 	END

c =========================================================================
 	
 	SUBROUTINE funcs(ss,a,ff,dfda,ma)
 	INTEGER ma
 	REAL*8 ss,ff,a(ma),dfda(ma)
 	REAL*8 arg,ex
 	ff=0.d0
	do 11 i=1,ma-1,2
	arg=a(i+1)*ss*ss
	ex=dexp(-arg)
	ff=ff+a(i)*ex
	dfda(i)=ex
	dfda(i+1)=-a(i)*ss*ss*ex
 11	continue
 	return
 	end 	

c ==========================================================================
	SUBROUTINE mrqmin(x,y,sig,ndata,a,ia,ma,covar,alpha,nca,
     *		chisq,funcs,alamda)

     	INTEGER ma,nca,ndata,ia(ma),MMAX
     	REAL*8 alamda,chisq,funcs,a(ma),alpha(nca,nca),covar(nca,nca),
     *		sig(ndata),x(ndata),y(ndata)
c
     	PARAMETER (MMAX=20)
	EXTERNAL funcs
	INTEGER j,k,l,m,mfit
	REAL*8 ochisq,atry(MMAX),beta(MMAX),da(MMAX)
	SAVE ochisq,atry,beta,da,mfit
	
c Initialization

	if(alamda.lt.0.d0) then
c
c determination of how many parameters are to be determined
c
        	mfit=0
		do 11 j=1,ma
			if(ia(j).ne.0) mfit=mfit+1
  11		continue
		alamda=0.001d0
	call mrqcof(x,y,sig,ndata,a,ia,ma,alpha,beta,nca,
     *    chisq,funcs)
		ochisq=chisq
		do 12 j=1,ma
			atry(j)=a(j)
  12		continue
	endif
	j=0
c
c After linearized fitting matrix, by augmenting diagonal elements
c
	do 14 l=1,ma
	  if(ia(l).ne.0) then
	  j=j+1
	  k=0
	  do 13 m=1,ma
	     if(ia(m).ne.0) then
	        k=k+1
	        covar(j,k)=alpha(j,k)
	     endif
  13	  continue
	  covar(j,j)=alpha(j,j)*(1.d0+alamda)
	  da(j)=beta(j)
	  endif
  14	continue
c
c Matrix solution. Once converged, evaluate covariance matrix.
c
	call gaussj(covar,mfit,nca,da,1,1)
	if(alamda.eq.0.d0) then
	  call covsrt(covar,nca,ma,ia,mfit)
	  return
	endif
	j=0
c 
c Did the trial succeed?
c
	do 15 l=1,ma
	   if(ia(l).ne.0) then
	      j=j+1
	      atry(l)=a(l)+da(j)
	   endif
  15	continue
  
  
  	do 51 iw=1,ma
  	write(6,*) a(iw),atry(iw)
  51	continue
  
  
	call mrqcof(x,y,sig,ndata,atry,ia,ma,covar,da,nca,chisq,funcs)
c
c Success, accept the new solution
c
	if(chisq.lt.ochisq) then
	   alamda=0.1d0*alamda
	   ochisq=chisq
	   j=0
	   do 17 l=1,ma
	      if(ia(l).ne.0) then
	         j=j+1
	         k=0
	         do 16 m=1,ma
	            if(ia(m).ne.0) then
	              k=k+1
	              alpha(j,k)=covar(j,k)
	            endif
  16	         continue
	         beta(j)=da(j)
	         a(l)=atry(l)
	      endif
  17	   continue
	 else
c
c False, increase alamda and return
c	 
	   alamda=10.d0*alamda	
	   chisq=ochisq
	endif
	return
	END
c
c **************************************************************************   
c 
	SUBROUTINE mrqcof(x,y,sig,ndata,a,ia,ma,alpha,beta,nalp,
     *	    chisq,funcs)
        INTEGER ma,nalp,ndata,ia(ma),MMAX
        REAL*8 chisq,a(ma),alpha(nalp,nalp),beta(ma),sig(ndata),
     *     x(ndata),y(ndata)
     	EXTERNAL funcs
     	PARAMETER (MMAX=20)
	INTEGER mfit,i,j,k,l,m
	REAL*8 dy,sig2i,wt,ymod,dyda(MMAX)     	
	mfit=0
	do 11 j=1,ma
	   if(ia(j).ne.0) mfit=mfit+1
 11	continue
c
c Initialize (symmetric) alpha, beta
c
	do 13 j=1,mfit
	   do 12 k=1,j
	      alpha(j,k)=0.d0
 12	   continue
 	   beta(j)=0.d0
 13	continue
 	chisq=0.d0
c
c Summation loop over all data
c
	do 16 i=1,ndata
	    call funcs(x(i),a,ymod,dyda,ma)
	    sig2i=1.d0/(sig(i)*sig(i))
	    dy=y(i)-ymod
	    j=0
	    do 15 l=1,ma
	      if(ia(l).ne.0) then
	         j=j+1
	         wt=dyda(l)*sig2i
	         k=0
	         do 14  m=1,l
	           if(ia(m).ne.0) then
	              k=k+1
	              alpha(j,k)=alpha(j,k)+wt*dyda(m)
	           endif
 14		 continue
 		 beta(j)=beta(j)+dy*wt
 	      endif
  15	    continue
c
c and find `chisq'
c
	chisq=chisq+dy*dy*sig2i
  16	continue
c
c Fill in the symmetric side
c
	do 18 j=2,mfit
	   do 17 k=1,j-1
	       alpha(k,j)=alpha(j,k)
 17	   continue
 18	continue
 	return
 	END
c
c ***************************************************************************	       	 	    	            	 	   	  
c
	SUBROUTINE gaussj(a,n,np,b,m,mp)
	INTEGER m,mp,n,np,NMAX
	REAL*8 a(np,np),b(np,mp)
	PARAMETER (NMAX=50)
	INTEGER i,icol,irow,j,k,l,ll,indxc(NMAX),indxr(NMAX),
     *      ipiv(NMAX)
     	REAL*8 big,dum,pivinv
     	do 11 j=1,n
     	   ipiv(j)=0
  11	continue
c
c This is the main loop over the columns to be reduced
c
	do 22 i=1,n
	     big=0.d0
c
c This is the outer loop of the search  for a pivot element
c
	     do 13 j=1,n
	        if(ipiv(j).ne.1) then
	            do 12 k=1,n
	               if (ipiv(k).eq.0) then
	                 if(dabs(a(j,k)).ge.big) then
	                    big=dabs(a(j,k))
	                    irow=j
	                    icol=k
	                 endif
	               else if (ipiv(k).gt.1) then
	                 pause 'singular matrix in gaussj (1)'
	               endif
 12		    continue
 		endif
 13	     continue
 	     ipiv(icol)=ipiv(icol)+1
c
c We now have the pivot element, so we interchange rows, if needed,
c to put the pivot element on the diagonal. The columns are not physically
c interchange, only relabeled: `indxc(i)', the column of the i'th pivot
c element, is the i'th column that is reduced, while `indxr(i)' is the
c row in which that pivot element was originally located. If `indxr(i) \ne
c indxc(i)' there is an implied column interchange. With this form of
c bookkeeping, the solution b's will end up in the correct order, and
c the inverse matrix will be scrambled by columns.
c
	     if(irow.ne.icol) then
	        do 14 l=1,n
	           dum=a(irow,l)
	           a(irow,l)=a(icol,l)
	           a(icol,l)=dum
 14		continue
 		do 15 l=1,m
 		   dum=b(irow,l)
 		   b(irow,l)=b(icol,l)
 		   b(icol,l)=dum
 15		continue
             endif
c
c We are now ready to devide the pivot row by the pivot element,
c located at irow and icol
c
	     indxr(i)=irow
	     indxc(i)=icol
	     if(a(icol,icol).eq.0.d0) pause
     * 'singular matrix in gaussj (2)'
             pivinv=1.d0/a(icol,icol)
             a(icol,icol)=1.d0
             do 16 l=1,n
                a(icol,l)=a(icol,l)*pivinv
 16	     continue
 	     do 17 l=1,m
 	         b(icol,l)=b(icol,l)*pivinv
 17	     continue
c
c Next, we reduce the rows... except for the pivot one, of course
c
 	     do 21 ll=1,n
 	        if(ll.ne.icol)then
 	           dum=a(ll,icol)
 	           a(ll,icol)=0.d0
 	           do 18 l=1,n
 	               a(ll,l)=a(ll,l)-a(icol,l)*dum
 18		   continue
 		   do 19 l=1,m
 		       b(ll,l)=b(ll,l)-b(icol,l)*dum
 19		   continue
 		endif
 21	     continue
 22	continue
c
c This was the end of the main loop over columns of the reduction.
c It only remains to unscramble the solution in view of the column
c interchanges. We do this by interchanging pairs of columns in the
c reverse order that the permutation was built up.
c
	do 24 l=n,1,-1
	   if(indxr(l).ne.indxc(l)) then
	     do 23 k=1,n
	        dum=a(k,indxr(l))
	        a(k,indxr(l))=a(k,indxc(l))
	        a(k,indxc(l))=dum
 23	     continue
           endif
 24	continue
 	return
 	END
c
c ************************************************************************** 	           	     		 		        	           
c
	SUBROUTINE covsrt(covar,npc,ma,ia,mfit)
	INTEGER ma,mfit,npc,ia(ma)
	REAL*8 covar(npc,npc)
c
c Expand in storage the covariance matrix `covar', so as to take
c into account parameters that are being held fixed. (For the latter,
c return zero covariances).
c
	INTEGER i,j,k
	REAL*8 swap
	do 12 i=mfit+1,ma
	   do 11 j=1,i
	      covar(i,j)=0.d0
	      covar(j,i)=0.d0
 11	   continue
 12	continue
 	k=mfit
 	do 15 j=ma,1,-1
 	   if(ia(j).ne.0) then
 	      do 13 i=1,ma
 	          swap=covar(i,k)
 	          covar(i,k)=covar(i,j)
 	          covar(i,j)=swap
 13	      continue
 	      do 14 i=1,ma
 	          swap=covar(k,i)
 	          covar(k,i)=covar(j,i)
 	          covar(j,i)=swap
 14	      continue
 	      k=k-1
 	   endif
 15	continue
 	return
 	END
c
c *************************************************************************** 	 	    	           	          	      
