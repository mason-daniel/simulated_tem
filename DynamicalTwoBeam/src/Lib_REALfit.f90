
!============================================================================= 
!                                  REALfit.f
!=============================================================================
!
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
!       added explicit interface to subroutine funcs
!       removed obsolete pause command
!       tweaked filenames to allow directory paths
!       no functional changes

!   Edit DRM 17/09/24
!       made into f08 compatible module

    module Lib_REALfit
!---^^^^^^^^^^^^^^^^^^
!*      f08 module containing functionality of old REALfit.f code
!*      fits the Doyle-Truner coefficients a and b to the scattering amplitudes.
!*
        use iso_fortran_env
        use Lib_inputReal
        implicit none
        private
        
        

        public      ::      fitDoyleTurnerCoefficients

  


    contains
!---^^^^^^^^
    
        subroutine fitDoyleTurnerCoefficients( a,b, f, chisq_bar )
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      given the scattering amplitudes f(s), refine and improve the Doyle-Turner coefficients a and b
    !*      return the average chi-square value
            real(kind=real64),dimension(5),intent(inout)        ::      a           
            real(kind=real64),dimension(5),intent(inout)        ::      b
            real(kind=real64),dimension(0:),intent(in)          ::      f           !   (0:LIB_INPUTREAL_NS) scattering amplitudes
            real(kind=real64),intent(out)                       ::      chisq_bar

            real(kind=real64),dimension(:),allocatable      ::      ss
            real(kind=real64),dimension(:),allocatable      ::      sig

            integer                 ::      nS,mab,ncab
            real(kind=real64)       ::      smin,smax,ds        
            integer                 ::      ii

            real(kind=real64),dimension(10)                     ::      ab
            real(kind=real64),dimension(10,10)                  ::      covar
            real(kind=real64),dimension(10,10)                  ::      alpha
            integer,dimension(10)                               ::      iab

            real(kind=real64)       ::      alamda
            real(kind=real64)       ::      chisq,oldchisq
            real(kind=real64)       ::      TOL = 1.0d-12

        !---    set the range of s = 1/4pi
            smin = LIB_INPUTREAL_SMIN
            smax = LIB_INPUTREAL_SMAX
            nS = LIB_INPUTREAL_NS
            ds = (smax-smin)/nS
            allocate(ss(0:nS))
            allocate(sig(0:nS))
            do ii = 0,nS
                ss(ii) = sMin + ds*ii
                sig(ii) = 0.1d0
            end do

        !---    pack the a and b D-T coefficients into a single array.
            mab = 10
            ncab = mab
            ab(1:9:2) = a(1:5)
            ab(2:10:2) = b(1:5)
            iab(1:10) = 1

            alamda=-1.0d0

            oldchisq = huge(1.0)
            do ii = 1,1000              !   don't loop forever
                call mrqmin(ss,f,sig,(nS+1),ab,iab,mab,covar,alpha,ncab,chisq,alamda)
                !print *,ii,chisq,oldchisq,abs(chisq-oldchisq)<TOL*nS
                if (abs(chisq-oldchisq)<TOL*nS) exit
                oldchisq = chisq
            end do

                    
        !---    unpack the a and b D-T coefficients
            a(1:5) = ab(1:9:2) 
            b(1:5) = ab(2:10:2) 
            chisq_bar = chisq/nS
            
            return
        end subroutine fitDoyleTurnerCoefficients
            



!     implicit real*8 (a-h,o-z)
!     parameter (ndata=501,ma=10)
!     real(kind=real64) s(ndata),f(ndata),sig(ndata),a(ma)
!     real(kind=real64) covar(ma,ma),alpha(ma,ma)
!     INTEGER ia(ma) 
!     CHARACTER*1 answer
!     CHARACTER*60 filename1,filename2
    
    
!     nca=ma

!     write(6,*) 
!      *'INPUT filename (=the output filename for inputREAL.exe)='
!     read(5,fmt='(a)') filename1
!     write(6,*) 'OUTPUT filename (not more than 60 characters)='
!     read(5,fmt='(a)') filename2

!     print *,"writing to file """//trim(filename2)//""""
!     open(unit=101,file=trim(filename1),action="read")
!     open(unit=102,file=trim(filename2),action="write")

!     do 20 j=1,ndata
!     read(101,*) k,s(k),f(k)
!  20    continue
!     close(101)
    
!     do 2 i=1,ndata
!     sig(i)=1.d-01
!  2    continue

!     write(6,*) 'input the initial 5 values of a(i)='
!     read(5,*) a(1),a(3),a(5),a(7),a(9)
!     write(6,*) 'input the initial 5 values of b(i)='
!     read(5,*) a(2),a(4),a(6),a(8),a(10)

!     do 3 i=1,ma
!     ia(i)=1
!  3    continue
!     alamda=-1.d0
!     do iiii = 1,100
!     call mrqmin(s,f,sig,ndata,a,ia,ma,covar,alpha,nca,
!      *          chisq,funcs,alamda)
! !         write(6,*) 'continue? (y/n)'
! !         read(5,*) answer
! !    if(answer.eq.'y') go to 4
!     end do  
!         write(6,*) 'chi^2=',(chisq/float(ndata)),'  alamda=',alamda
!     write(*,*) 'writing to the file ',filename2

!         write(102,98) a(1),a(3),a(5),a(7),a(9)
!         write(102,*)
!         write(102,99) a(2),a(4),a(6),a(8),a(10)
!         write(102,*)
!         write(102,*) (chisq/float(ndata))

!     close(101)
!     close(102)
  
!   98    format(1x,6Ha(i)= ,5(e12.5,1x))
!   99    format(1x,6Hb(i)= ,5(e12.5,1x))
  
!      STOP
!      END

! =========================================================================
     
        pure SUBROUTINE funcs(ss,a,ff,dfda,ma)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
            INTEGER,intent(in)              ::      ma
            real(kind=real64),intent(in)    ::      ss,a(ma)
            real(kind=real64),intent(out)   ::      ff,dfda(ma)
            real(kind=real64)      ::      arg,ex
            integer                ::      ii
            ff=0.d0
            do ii=1,ma-1,2
                arg=a(ii+1)*ss*ss
                ex=dexp(-arg)
                ff=ff+a(ii)*ex
                dfda(ii)=ex
                dfda(ii+1)=-a(ii)*ss*ss*ex
            end do
            return
        end SUBROUTINE funcs    

! ==========================================================================
    SUBROUTINE mrqmin(x,y,sig,ndata,a,ia,ma,covar,alpha,nca, chisq,alamda)

         INTEGER ma,nca,ndata,ia(ma),MMAX
!          real(kind=real64) alamda,chisq,funcs,a(ma),alpha(nca,nca),covar(nca,nca),
!      *        sig(ndata),x(ndata),y(ndata)
         real(kind=real64) alamda,chisq,a(ma),alpha(nca,nca),covar(nca,nca), sig(ndata),x(ndata),y(ndata)
!       
         PARAMETER (MMAX=20)
   !EXTERNAL funcs
    INTEGER j,k,l,m,mfit!,iw
    real(kind=real64) ochisq,atry(MMAX),beta(MMAX),da(MMAX)
    SAVE ochisq,atry,beta,da,mfit
    
! Initialization

    if(alamda.lt.0.d0) then
!
! determination of how many parameters are to be determined
!
            mfit=0
        do 11 j=1,ma
            if(ia(j).ne.0) mfit=mfit+1
  11        continue
        alamda=0.001d0
    call mrqcof(x,y,sig,ndata,a,ia,ma,alpha,beta,nca, chisq)
        ochisq=chisq
        do 12 j=1,ma
            atry(j)=a(j)
  12        continue
    endif
    j=0
!
! After linearized fitting matrix, by augmenting diagonal elements
!
    do 14 l=1,ma
      if(ia(l).ne.0) then
      j=j+1
      k=0
      do 13 m=1,ma
         if(ia(m).ne.0) then
            k=k+1
            covar(j,k)=alpha(j,k)
         endif
  13      continue
      covar(j,j)=alpha(j,j)*(1.d0+alamda)
      da(j)=beta(j)
      endif
  14    continue
!
! Matrix solution. Once converged, evaluate covariance matrix.
!
    call gaussj(covar,mfit,nca,da,1,1)
    if(alamda.eq.0.d0) then
      call covsrt(covar,nca,ma,ia,mfit)
      return
    endif
    j=0
! 
! Did the trial succeed?
!
    do 15 l=1,ma
       if(ia(l).ne.0) then
          j=j+1
          atry(l)=a(l)+da(j)
       endif
  15    continue
  
  
!       do 51 iw=1,ma
!       write(6,*) a(iw),atry(iw)
!   51    continue
  
  
    call mrqcof(x,y,sig,ndata,atry,ia,ma,covar,da,nca,chisq)
!
! Success, accept the new solution
!
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
  16             continue
             beta(j)=da(j)
             a(l)=atry(l)
          endif
  17       continue
     else
!
! False, increase alamda and return
!     
       alamda=10.d0*alamda    
       chisq=ochisq
    endif
    return
    END
!
! **************************************************************************   
! 
    SUBROUTINE mrqcof(x,y,sig,ndata,a,ia,ma,alpha,beta,nalp, chisq)
        INTEGER ma,nalp,ndata,ia(ma),MMAX
        real(kind=real64) chisq,a(ma),alpha(nalp,nalp),beta(ma),sig(ndata), x(ndata),y(ndata)
        ! EXTERNAL funcs
         PARAMETER (MMAX=20)
    INTEGER mfit,i,j,k,l,m
    real(kind=real64) dy,sig2i,wt,ymod,dyda(MMAX)         
    mfit=0
    do 11 j=1,ma
       if(ia(j).ne.0) mfit=mfit+1
 11    continue
!
! Initialize (symmetric) alpha, beta
!
    do 13 j=1,mfit
       do 12 k=1,j
          alpha(j,k)=0.d0
 12       continue
        beta(j)=0.d0
 13    continue
     chisq=0.d0
!
! Summation loop over all data
!
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
 14         continue
          beta(j)=beta(j)+dy*wt
           endif
  15        continue
!
! and find `chisq'
!
    chisq=chisq+dy*dy*sig2i
  16    continue
!
! Fill in the symmetric side
!
    do 18 j=2,mfit
       do 17 k=1,j-1
           alpha(k,j)=alpha(j,k)
 17       continue
 18    continue
     return
     END
!
! ***************************************************************************                                                          
!
    SUBROUTINE gaussj(a,n,np,b,m,mp)
    INTEGER m,mp,n,np,NMAX
    real(kind=real64) a(np,np),b(np,mp)
    PARAMETER (NMAX=50)
    INTEGER i,icol,irow,j,k,l,ll,indxc(NMAX),indxr(NMAX),ipiv(NMAX)


        
         real(kind=real64) big,dum,pivinv
         do 11 j=1,n
            ipiv(j)=0
  11    continue
!
! This is the main loop over the columns to be reduced
!
    do 22 i=1,n
         big=0.d0
!
! This is the outer loop of the search  for a pivot element
!
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
                     stop 'singular matrix in gaussj (1)'
                   endif
 12            continue
         endif
 13         continue
          ipiv(icol)=ipiv(icol)+1
!
! We now have the pivot element, so we interchange rows, if needed,
! to put the pivot element on the diagonal. The columns are not physically
! interchange, only relabeled: `indxc(i)', the column of the i'th pivot
! element, is the i'th column that is reduced, while `indxr(i)' is the
! row in which that pivot element was originally located. If `indxr(i) \ne
! indxc(i)' there is an implied column interchange. With this form of
! bookkeeping, the solution b's will end up in the correct order, and
! the inverse matrix will be scrambled by columns.
!
         if(irow.ne.icol) then
            do 14 l=1,n
               dum=a(irow,l)
               a(irow,l)=a(icol,l)
               a(icol,l)=dum
 14        continue
         do 15 l=1,m
            dum=b(irow,l)
            b(irow,l)=b(icol,l)
            b(icol,l)=dum
 15        continue
             endif
!
! We are now ready to devide the pivot row by the pivot element,
! located at irow and icol
!
         indxr(i)=irow
         indxc(i)=icol
         if(a(icol,icol).eq.0.d0) stop 'singular matrix in gaussj (2)'
             pivinv=1.d0/a(icol,icol)
             a(icol,icol)=1.d0
             do 16 l=1,n
                a(icol,l)=a(icol,l)*pivinv
 16         continue
          do 17 l=1,m
              b(icol,l)=b(icol,l)*pivinv
 17         continue
!
! Next, we reduce the rows... except for the pivot one, of course
!
          do 21 ll=1,n
             if(ll.ne.icol)then
                dum=a(ll,icol)
                a(ll,icol)=0.d0
                do 18 l=1,n
                    a(ll,l)=a(ll,l)-a(icol,l)*dum
 18           continue
            do 19 l=1,m
                b(ll,l)=b(ll,l)-b(icol,l)*dum
 19           continue
         endif
 21         continue
 22    continue
!
! This was the end of the main loop over columns of the reduction.
! It only remains to unscramble the solution in view of the column
! interchanges. We do this by interchanging pairs of columns in the
! reverse order that the permutation was built up.
!
    do 24 l=n,1,-1
       if(indxr(l).ne.indxc(l)) then
         do 23 k=1,n
            dum=a(k,indxr(l))
            a(k,indxr(l))=a(k,indxc(l))
            a(k,indxc(l))=dum
 23         continue
           endif
 24    continue
     return
     END
!
! **************************************************************************                                                                 
!
    SUBROUTINE covsrt(covar,npc,ma,ia,mfit)
    INTEGER ma,mfit,npc,ia(ma)
    real(kind=real64) covar(npc,npc)
!
! Expand in storage the covariance matrix `covar', so as to take
! into account parameters that are being held fixed. (For the latter,
! return zero covariances).
!
    INTEGER i,j,k
    real(kind=real64) swap
    do 12 i=mfit+1,ma
       do 11 j=1,i
          covar(i,j)=0.d0
          covar(j,i)=0.d0
 11       continue
 12    continue
     k=mfit
     do 15 j=ma,1,-1
        if(ia(j).ne.0) then
           do 13 i=1,ma
               swap=covar(i,k)
               covar(i,k)=covar(i,j)
               covar(i,j)=swap
 13          continue
           do 14 i=1,ma
               swap=covar(k,i)
               covar(k,i)=covar(j,i)
               covar(j,i)=swap
 14          continue
           k=k-1
        endif
 15    continue
     return
     END
!
! ***************************************************************************                                                     

    end module Lib_REALfit