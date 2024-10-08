!============================================================================= 
!                                  TDSfit.f
!=============================================================================
!
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
!

!   DRM Sept 2024
!   change to f08 module

    module Lib_TDSfit
!---^^^^^^^^^^^^^^^^^
        use Lib_inputTDS
        use iso_fortran_env
        implicit none
        private

        public          ::      TDSfit


    contains
!---^^^^^^^^

        subroutine TDSfit(a,b , f,chisq_bar )
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

            real(kind=real64),dimension(:),intent(inout)        ::      a           
            real(kind=real64),dimension(:),intent(inout)        ::      b
            real(kind=real64),dimension(0:),intent(in)          ::      f           !   (0:LIB_INPUTREAL_NS) scattering amplitudes
            real(kind=real64),intent(out)                       ::      chisq_bar

            real(kind=real64),dimension(2*size(a))              ::      ab
            integer                                             ::      mab
            integer                                             ::      na
            integer,dimension(2*size(a))                        ::      iab
            real(kind=real64),dimension(2*size(a),2*size(a))    ::      covar
            real(kind=real64),dimension(2*size(a),2*size(a))    ::      alpha
            real(kind=real64),dimension(:),allocatable          ::      ss
            real(kind=real64),dimension(:),allocatable          ::      sig

            real(kind=real64)       ::      alamda
            real(kind=real64)       ::      chisq,oldchisq
            real(kind=real64)       ::      TOL = 1.0d-12

            real(kind=real64)       ::      smin,smax,ds        
            integer                 ::      nS,ii


        !---    set the range of s = 1/4pi
            smin = LIB_INPUTTDS_SMIN
            smax = LIB_INPUTTDS_SMAX
            nS = LIB_INPUTTDS_NS
            ds = (smax-smin)/nS
            allocate(ss(0:nS))
            allocate(sig(0:nS))
            do ii = 0,nS
                ss(ii) = sMin + ds*ii
                sig(ii) = 0.1d0
            end do


        !---    pack input arguments a,b into a single array
            na = size(a)            !   number of Doyle-turner coefficients (5)
            mab = 2*na              !   number of coeffs in combined array (10)
            ab(1:mab-1:2) = a(1:na)
            ab(2:mab:2)   = b(1:na)

            iab = 1
            alamda = -1.0d0


            oldchisq = huge(1.0)
            do ii = 1,1000          !   don't loop forever
                call mrqmin(ss,f,sig,(nS+1),ab,iab,mab,covar,alpha,mab,chisq,alamda)
                write(6,*) 'chi^2=',(chisq/float(nS+1)),'  alamda=',alamda
                if (abs(chisq-oldchisq)<TOL*nS) exit
                oldchisq = chisq
            end do

                    
        !---    unpack the a and b D-T coefficients
            a(1:5) = ab(1:9:2) 
            b(1:5) = ab(2:10:2) 
            chisq_bar = chisq/nS
            

            return
        end subroutine TDSfit
            


!---
     
        pure subroutine funcs(s,a,ff,dfda,ma)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
            integer,intent(in)                              ::  ma
            real(kind=real64),intent(in)                    ::  s
            real(kind=real64),dimension(:),intent(in)       ::  a           !   (1:ma)
            real(kind=real64),intent(out)                   ::  ff
            real(kind=real64),dimension(:),intent(out)      ::  dfda        !   (1:ma)

            real(kind=real64)       ::      xx
            integer                 ::      ii

            ff=0.d0           
            do ii = 1,ma-1,2
                xx = a(ii+1)*s*s
                xx = exp(-xx)
                ff = ff + a(ii)*xx
                dfda(ii) = xx
                dfda(ii+1) = -a(ii)*s*s*xx
            end do
            return
        end subroutine funcs
 

! ==========================================================================
!
    SUBROUTINE mrqmin(x,y,sig,ndata,a,ia,ma,covar,alpha,nca,chisq,alamda)
         INTEGER ma,nca,ndata,ia(ma),MMAX
!          REAL*8 alamda,chisq,funcs,a(ma),alpha(nca,nca),covar(nca,nca),
!      *        sig(ndata),x(ndata),y(ndata)
         REAL*8 alamda,chisq,a(ma),alpha(nca,nca),covar(nca,nca),sig(ndata),x(ndata),y(ndata)
!
         PARAMETER (MMAX=20)
    ! EXTERNAL funcs
    INTEGER j,k,l,m,mfit,iw
    REAL*8 ochisq,atry(MMAX),beta(MMAX),da(MMAX)
    SAVE ochisq,atry,beta,da,mfit
    
! Initialization
        print *,"mrqmin ",a,ma,ndata
    if(alamda.lt.0.d0) then
!
! determination of how many parameters are to be determined
!
            mfit=0
        do 11 j=1,ma
            if(ia(j).ne.0) mfit=mfit+1
  11        continue
        alamda=0.001d0
    call mrqcof(x,y,sig,a,ia,alpha,beta,chisq )
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
      call covsrt(covar,ma,ia,mfit)
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
  
  
      do 51 iw=1,ma
      write(6,*) a(iw),atry(iw)
  51    continue
  
  
    call mrqcof(x,y,sig,atry,ia,covar,da,chisq )
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
        subroutine mrqcof(x,y,sig,a,ia,alpha,beta,chisq)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^            
            integer,dimension(:),intent(in)                 ::      ia          !   (1:ma)
            !INTEGER ma,nalp,ndata,ia(ma) 
            real(kind=real64),dimension(:),intent(in)       ::      a           !   (1:ma)
            real(kind=real64),dimension(:),intent(in)       ::      x,y,sig     !   (1:ndata)
            real(kind=real64),dimension(:,:),intent(out)    ::      alpha       !   (1:nalp,1:nalp)
            real(kind=real64),dimension(:),intent(out)      ::      beta        !   (1:nalp)
            real(kind=real64),intent(out)                   ::      chisq
            !REAL*8 chisq,a(ma),alpha(nalp,nalp),beta(ma),sig(ndata),x(ndata),y(ndata)
         
            integer,parameter       ::      MMAX = 20
            integer                 ::      ma,ndata
            integer                 ::      mfit,ii,jj,kk,ll,mm
            real(kind=real64)       ::      dy,isig2,wt,ymod
            real(kind=real64),dimension(MMAX)       ::      dyda 
            
        !---    establish size of problem
            ma = size(ia)
            ndata = size(x)
            mfit = count(ia/=0)

        !   Initialize (symmetric) alpha, beta
            do jj = 1,mfit
                alpha(jj,1:jj) = 0.0d0
                beta(jj) = 0.0d0
            end do
            chisq = 0.0d0

            print *,"mrqcof ma,ndata,nfit ",ma,ndata,mfit
        !   Summation loop over all data            
            do ii = 1,ndata
                call funcs(x(ii),a,ymod,dyda,ma)
                isig2 = 1/(sig(ii)*sig(ii))
                dy = y(ii) - ymod
                jj = 0
                do ll = 1,ma
                    if (ia(ll)/=0) then
                        jj = jj + 1
                        wt = dyda(ll) * isig2
                        kk = 0
                        do mm = 1,ll
                            if (ia(mm).ne.0) then
                                kk = kk+1
                                alpha(jj,kk) = alpha(jj,kk) + wt*dyda(mm)
                            end if
                        end do
                        beta(jj) = beta(jj) + dy*wt
                    end if                      
                end do  !   ll
                chisq = chisq + dy*dy*isig2
            end do  !   ii

        !---    symmetrise
            do jj = 2,mfit
                do kk = 1,jj-1
                    alpha(kk,jj) = alpha(jj,kk)
                end do
            end do


! ! Summation loop over all data
! !
!     do 16 i=1,ndata
!         call funcs(x(i),a,ymod,dyda,ma)
!         sig2i=1.d0/(sig(i)*sig(i))
!         dy=y(i)-ymod
!         j=0
!         do 15 l=1,ma
!           if(ia(l).ne.0) then
!              j=j+1
!              wt=dyda(l)*sig2i
!              k=0
!              do 14  m=1,l
!                if(ia(m).ne.0) then
!                   k=k+1
!                   alpha(j,k)=alpha(j,k)+wt*dyda(m)
!                endif
!  14         continue
!           beta(j)=beta(j)+dy*wt
!            endif
!   15        continue
! !
! ! and find `chisq'
! !
!     chisq=chisq+dy*dy*sig2i
!   16    continue
! !
! ! Fill in the symmetric side
! !
!     do 18 j=2,mfit
!        do 17 k=1,j-1
!            alpha(k,j)=alpha(j,k)
!  17       continue
!  18    continue
            return
        end subroutine mrqcof
!
! ***************************************************************************                                                          
!
    SUBROUTINE gaussj(a,n,np,b,m,mp)
    INTEGER m,mp,n,np,NMAX
    REAL*8 a(np,np),b(np,mp)
    PARAMETER (NMAX=50)
!
    INTEGER i,icol,irow,j,k,l,ll,indxc(NMAX),indxr(NMAX),ipiv(NMAX)
         REAL*8 big,dum,pivinv


         irow = 0; icol = 0      !   to avoid maybe uninitialized warning ( can't be 

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
        subroutine covsrt(covar,ma,ia,mfit)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
            real(kind=real64),dimension(:,:),intent(inout)  ::      covar
            integer,intent(in)                              ::      ma,mfit
            integer,dimension(:),intent(in)                 ::      ia
!
            integer             ::      ii,jj,kk
            real(kind=real64)   ::      swap
    
            do ii = mfit+1,ma
                covar(ii,1:ii) = 0.d00
                covar(1:ii,ii) = 0.0d0
            end do

            kk = mfit
            do jj = ma,1,-1
                if (ia(jj)/=0) then
                    do ii = 1,ma
                        swap = covar(ii,kk)
                        covar(ii,kk) = covar(ii,jj)
                        covar(ii,jj) = swap
                    end do
                    do ii = 1,ma
                        swap = covar(kk,ii)
                        covar(kk,ii) = covar(jj,ii)
                        covar(jj,ii) = swap
                    end do
                    kk = kk - 1
                end if
            end do

            return
        end subroutine covsrt


!         c
!         c ************************************************************************** 	           	     		 		        	           
!         c
!             SUBROUTINE covsrt(covar,npc,ma,ia,mfit)
!             INTEGER ma,mfit,npc,ia(ma)
!             REAL*8 covar(npc,npc)
!         c
!             INTEGER i,j,k
!             REAL*8 swap
!             do 12 i=mfit+1,ma
!                do 11 j=1,i
!                   covar(i,j)=0.d0
!                   covar(j,i)=0.d0
!          11	   continue
!          12	continue
!              k=mfit
!              do 15 j=ma,1,-1
!                 if(ia(j).ne.0) then
!                    do 13 i=1,ma
!                        swap=covar(i,k)
!                        covar(i,k)=covar(i,j)
!                        covar(i,j)=swap
!          13	      continue
!                    do 14 i=1,ma
!                        swap=covar(k,i)
!                        covar(k,i)=covar(j,i)
!                        covar(j,i)=swap
!          14	      continue
!                    k=k-1
!                 endif
!          15	continue
!              return
!              END
!         c
! ! ***************************************************************************                                                     


    end module Lib_TDSfit



