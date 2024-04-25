
    module Lib_RAPFanalysis
!---^^^^^^^^^^^^^^^^^^^^^^^
        use iso_fortran_env
        implicit none
        private
        
    !---
     
        
        real(kind=real64),private,parameter ::      PI = 3.1415926535898d0  
        integer,parameter                   ::      NDIV = 40000
        real(kind=real64),parameter         ::      LOGSONWMIN = -50
        real(kind=real64),parameter         ::      LOGSONWMAX =  20
        real(kind=real64),parameter         ::      QONWMIN = 0
        real(kind=real64),parameter         ::      QONWMAX = 400
        real(kind=real64),parameter         ::      DLOGSONW = (LOGSONWMAX-LOGSONWMIN)/NDIV
        real(kind=real64),parameter         ::      DQONW = (QONWMAX-QONWMIN)/NDIV
            
    !---
    
        logical,public          ::      RAPFanalysis_dbg = .false.
                           
    !---
        public      ::      computeMomentsRPS
        public      ::      findMeanAndStdDev
     
    contains
!---^^^^^^^^

        pure real(kind=real64) function safeExp( x ) 
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
            real(kind=real64),intent(in)        ::      x
            if (x<-200.0d0) then
                safeExp = 0.d0
            else if (x>200.0d0) then
                safeExp = exp( 200.0d0 )
            else
                safeExp = exp( x )
            end if
            return
        end function safeExp

        function lognormal( s,logsonw,w ) result( p )
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      returns the lognormal distribution with unit mean
    !*      p = 1/sqrt(2 pi) 1/s 1/w Esp[ - (log s)^2 / 2 w^2 ]
            real(kind=real64),intent(in)        ::      s,logsonw,w
            real(kind=real64)                   ::      p
            real(kind=real64),parameter         ::      ISQRT2PI = 1/sqrt(2*PI)
            
            !logsonw = log(s)/w    
            if (s<1.0d-16) then
                p = 0
            else if (s>1.0d16) then
                p = 0
            else
                p = safeExp( - logsonw*logsonw/2 )
                p = ISQRT2PI * p / (w*s)
            end if
            return
        end function lognormal
            

        function rps( q,w ) result( c )
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      compute the radial power function assuming equal intensity Gaussian spots distributed with 
    !*      according to a lognormal distribution with width w and unit mean
            real(kind=real64),intent(in)        ::      q,w
            real(kind=real64)                   ::      c
             
            
            real(kind=real64)           ::      s1,s2,pp,gg,c1,c2
            
            real(kind=real64)           ::      logsonw  , sump,p1
            
            
            integer             ::      ii
            
            !aa = 4*PI*PI*q*q 
            
        !---    test: s = w
            !c = 2*PI*q*exp( -q*q*w*w )
            !return    
            
            
            
            
            s2 = safeExp( w * LOGSONWMIN )      !   smallest value of s to consider
            pp = lognormal( s2,LOGSONWMIN,w )            
            gg = safeExp( - q*q * s2*s2 )
            c2 = gg*pp
                        
            c = 0
            sump = 0.0d0
            do ii = 1,NDIV
                s1 = s2 ; c1 = c2 ; p1 = pp
                logsonw = LOGSONWMIN + ii*DLOGSONW                
                s2 = safeExp( w * logsonw )     
                                  
                pp = lognormal( s2,logsonw,w )
                gg = safeExp( - q*q * s2*s2 )
                c2 = gg*pp
                         
                c  = c + (s2-s1)*(c2+c1)             !   trapezium rule ( note s 2 )
                
                sump = sump + (s2-s1)*(pp+p1) 
            end do
                 
            c = c * PI*q !/ sump                      !   normalise ( remembering / 2 )
            
            !if (abs(sump-1)>0.001d0) print *,"sump ",sump
            !print *,"rps q,w ",q,w,c
            return
        end function rps
            
        subroutine computeMomentsRPS( w, m1,m2,m3 )
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ 
    !*      compute the first and second moments of the rps
    !*      assuming gaussian blobs and lognormal distribution
            real(kind=real64),intent(in)        ::      w
            real(kind=real64),intent(out)       ::      m1,m2,m3
             
            
            
            real(kind=real64)           ::      qq,cc,m0
            integer                     ::      ii
            
            qq   = w*QONWMIN        !  = 0
            cc   = rps( qq,w )      !  = 0
            
            m0 = 0 ; m1 = 0 ; m2 = 0 ; m3 = 0
           ! print *,"computeMomentsRPS info - qmin,w,qmax = ",w*QONWMIN,w,w*( QONWMIN + NDIV*DQONW )
            do ii = 0,NDIV-1
                qq = w*( QONWMIN + ii*DQONW )
                
                !print *,ii,qq,cc
                
                m0  = m0 + cc                !   computed at qq 
                m1  = m1 + qq*cc
                m2  = m2 + qq*qq*cc
                m3  = m3 + qq*qq*qq*cc

                qq  = qq + w*DQONW/2
                cc  = rps( qq,w )           !   computed at qq + 1/2
                m0  = m0 + 4*cc              !   Simpsons rule
                m1  = m1 + 4*qq*cc
                m2  = m2 + 4*qq*qq*cc
                m3  = m3 + 4*qq*qq*qq*cc
                
                qq  = qq + w*DQONW/2
                cc  = rps( qq,w )           !   computed at qq + 1
                m0  = m0 + cc      
                m1  = m1 + qq*cc   
                m2  = m2 + qq*qq*cc
                m3  = m3 + qq*qq*qq*cc
            end do
            
        !---    normalise integrals
            
            !print *,"computeMomentsRPS info - m0,m1,m2 ",(/m0,m1,m2/)*DQONW*w/6
            m1 = m1/m0
            m2 = m2/m0
            m3 = m3/m0
            
            !stop
            
            return
        end subroutine computeMomentsRPS

 
       subroutine findMeanAndStdDev( q,q2, rbar,rstd,r2bar )
   !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
   !*      given an input first and second moment <q>,<q^2>
   !*      find the most probable width of a log-normal distribution
   !*      then find the true size scale
   !*      then use the properties of the lognormal to find the real space size and std dev and mean square size.
           real(kind=real64),intent(in)                ::      q,q2
           real(kind=real64),intent(out)               ::      rbar,rstd
           real(kind=real64),intent(out),optional      ::      r2bar
           
           real(kind=real64),dimension(33),parameter  ::      W  = (/     &
                                     1.0000000000000000D-002  ,     &
                                     1.1547819846894580D-002  ,     &
                                     1.3335214321633230D-002  ,     &
                                     1.5399265260594912D-002  ,     &
                                     1.7782794100389226D-002  ,     &
                                     2.0535250264571453D-002  ,     &
                                     2.3713737056616540D-002  ,     &
                                     2.7384196342643607D-002  ,     &
                                     3.1622776601683784D-002  ,     &
                                     3.6517412725483755D-002  ,     &
                                     4.2169650342858217D-002  ,     &
                                     4.8696752516586297D-002  ,     &
                                     5.6234132519034884D-002  ,     &
                                     6.4938163157621118D-002  ,     &
                                     7.4989420933245551D-002  ,     &
                                     8.6596432336006529D-002  ,     &
                                     9.9999999999999978D-002  ,     &
                                    0.11547819846894578d0     ,     &
                                    0.13335214321633240d0     ,     &
                                    0.15399265260594916d0     ,     &
                                    0.17782794100389226d0     ,     &
                                    0.20535250264571459d0     ,     &
                                    0.23713737056616546d0     ,     &
                                    0.27384196342643607d0     ,     &
                                    0.31622776601683789d0     ,     &
                                    0.36517412725483767d0     ,     &
                                    0.42169650342858223d0     ,     &
                                    0.48696752516586306d0     ,     &
                                    0.56234132519034907d0     ,     &
                                    0.64938163157621132d0     ,     &
                                    0.74989420933245576d0     ,     &
                                    0.86596432336006535d0     ,     &
                                     1.0000000000000000d0     /)
           
           real(kind=real64),dimension(33),parameter  ::      MUQ   = (/     &
                                0.88644812551211793d0,        &     
                                0.88652242299810435d0,        &     
                                0.88662100281460343d0,        &     
                                0.88675247502872256d0,        &     
                                0.88692782637369461d0,        &     
                                0.88716171510236153d0,        &     
                                0.88747370670236447d0,        &     
                                0.88788992490665752d0,        &     
                                0.88844526453484574d0,        &     
                                0.88918636229947001d0,        &     
                                0.89017559397952462d0,        &     
                                0.89149646815843631d0,        &     
                                0.89326093219080460d0,        &     
                                0.89561931759101820d0,        &     
                                0.89877396588497571d0,        &     
                                0.90299805237847786d0,        &     
                                0.90866186688258588d0,        &     
                                0.91626999641251461d0,        &     
                                0.92651480639662909d0,        &     
                                0.94035492519035191d0,        &     
                                0.95913322199353113d0,        &     
                                0.98475922257449888d0,        &     
                                 1.0200004756934731d0,        &     
                                 1.0689654533160899d0,        &     
                                 1.1379378972338530d0,        &     
                                 1.2368874715107214d0,        &     
                                 1.3823536845773430d0,        &     
                                 1.6032954844162379d0,        &     
                                 1.9538214590194602d0,        &     
                                 2.5432912691230105d0,        &     
                                 3.6149216141347091d0,        &     
                                 5.7769064246748929d0,        &     
                                 10.763779613679084d0       /)
           
           
           real(kind=real64),dimension(33),parameter  ::      MU2Q2  = (/     &
                                 1.0005982786888210d0           ,           &   
                                 1.0008004199223273d0           ,           &   
                                 1.0010675370417095d0           ,           &   
                                 1.0014238369180208d0           ,           &   
                                 1.0018991677350668d0           ,           &   
                                 1.0025333826248637d0           ,           &   
                                 1.0033797464581300d0           ,           &   
                                 1.0045095025980137d0           ,           &   
                                 1.0060180360540680d0           ,           &   
                                 1.0080332231628431d0           ,           &   
                                 1.0107268004421270d0           ,           &   
                                 1.0143299454527257d0           ,           &   
                                 1.0191548098015206d0           ,           &   
                                 1.0256245972817493d0           ,           &   
                                 1.0343161456625924d0           ,           &   
                                 1.0460212203306549d0           ,           &   
                                 1.0618365465453545d0           ,           &   
                                 1.0832992936457495d0           ,           &   
                                 1.1125968245810001d0           ,           &   
                                 1.1529022079011693d0           ,           &   
                                 1.2089311953225719d0           ,           &   
                                 1.2879063328008604d0           ,           &   
                                 1.4013061908449240d0           ,           &   
                                 1.5682126409730404d0           ,           &   
                                 1.8221188003893858d0           ,           &   
                                 2.2257921156511808d0           ,           &   
                                 2.9065524275197601d0           ,           &   
                                 4.1488211074550758d0           ,           &   
                                 6.6683109420152613d0           ,           &   
                                 12.555751081114192d0           ,           &   
                                 29.195302173790651d0           ,           &   
                                 89.739991842172046d0           ,           &   
                                 384.04821468487228d0   /)   
                                  
            real(kind=real64),dimension(size(MUQ)),parameter  ::      RATIO  = MU2Q2/(MUQ*MUQ)                   
                                  
                                  
            real(kind=real64)       ::      q2onqq
            real(kind=real64)       ::      mu,ww,xx
            integer                 ::      ii
            
            q2onqq = q2/(q*q)
            if (RAPFanalysis_dbg) print *,"Lib_RAPFanalysis::findMeanAndStdDev dbg - q2onqq ",4/PI,q2onqq,RATIO(size(MUQ))
            if (q2onqq*PI<4) then
                ww = 0.0d0           !   essentially zero width in this distribution.                   
                mu = sqrt(PI)/(2*q)                
            else if (q2onqq>=RATIO(size(MUQ))) then
                !   the width of the distribution is greater than can be handled by a gaussian blob
                !   this could mean that the underlying objects are disks rather than smooth.
                !   set a very simple radius instead
                
                rbar = sqrt(PI)/(2*q)
                ww = sqrt( max(0.0d0,q2 - q*q) )
                ww = abs( 1/(q-ww) - 1/(q+ww) )
                rstd = sqrt(PI)*ww
                if (present(r2bar)) r2bar = rbar*rbar + rstd*rstd
                return
            else            
                do ii = 1,size(MUQ)-1
                    !print *,ii,RATIO(ii),q2onqq,RATIO(ii+1)
                    if ( (RATIO(ii)-q2onqq)*(RATIO(ii+1)-q2onqq) <= 0) then
                        !   sign changes in this interval = have found solution
                        xx = (RATIO(ii+1)-q2onqq)/(RATIO(ii+1)-RATIO(ii))       !   fraction of distance across interval
                        !print *,"interval ",ii," + ",xx
                        ww = W(ii) + (W(ii+1)-W(ii))*xx              
                        !print *,"w ",ww
                        mu = MUQ(ii) + (MUQ(ii+1)-MUQ(ii))*xx                   !   good value for mu q
                        !print *,"mu q ",mu 
                        if (RAPFanalysis_dbg) print *,"Lib_RAPFanalysis::findMeanAndStdDev dbg - ratio interval ",ii,"+",xx," muq = ",mu," w = ",ww
                        mu = mu / q 
                        !print *,"mu ",mu
                        exit
                    end if
                end do
            end if  
            
            !if (RAPFanalysis_dbg) &          
            print *,"Lib_RAPFanalysis::findMeanAndStdDev info - <q^2>/<q>^2 = ",q2onqq," lognormal mu = ",mu," w = ",ww
             
            
            rbar = mu * exp( ww*ww/2 )
            rstd = sqrt( max(0.0d0,exp( ww*ww ) - 1 ) )*rbar
            if (present(r2bar)) r2bar = mu*mu * exp( 2*ww*ww )
            return
        end subroutine findMeanAndStdDev
        
    end module Lib_RAPFanalysis
    
    
!   gfortran -ffree-line-length-256 -Og -g Lib_RAPFanalysis.f90 -o testRAPFanalysis.exe

   
!   program testRAPFanalysiss
!!---^^^^^^^^^^^^^^^^^^^^^^^^
!       use iso_fortran_env
!       use Lib_RAPFanalysis
!       implicit none
!         
!       real(kind=real64)        ::      ww
!       real(kind=real64)        ::      m1,m2,m3
!       real(kind=real64)        ::      dbar,dstd
!       integer          ::      ii
!       
!       call findMeanAndStdDev( 7.6617853637458599d-002,8.6021500016066205d-003 , dbar,dstd )
!       print *,"d = ",dbar," +/- ",dstd
!       
!       do ii = -32,0
!            ww = 10.0d0**(ii/16.0d0)
!            call computeMomentsRPS( ww, m1,m2,m3 )
!            print *,"ww ,m1,m2 ",ww,m1,m2,m3
!       end do
!       
!       
!      
!       print *,""
!       print *,"done"
!       print *,""
!       
!   end program testRAPFanalysiss
!   
!!        
    
        
        
    