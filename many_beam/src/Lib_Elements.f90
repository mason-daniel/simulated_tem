
    module Lib_Elements
!---^^^^^^^^^^^^^^^^^^^
!*      A simple module which stores some basic facts about elements
!*
!*      Daniel Mason
!*      (c) UKAEA Sept 2024
!*
!*      version history
!*      0.0.1           Sept 2024       First version
!*

        use iso_fortran_env
        implicit none
        private


        
    !---    define the elements recognised        
        integer,private,parameter                               ::      NELNAMES = 50
        integer,public,parameter                                ::      ELNAME_LEN = 2
        integer,public,parameter                                ::      LATTNAME_LEN = 3
        integer,public,parameter                                ::      DT_COEFFS = 5

        character(len=ELNAME_LEN),dimension(NELNAMES),parameter  ::      ELNAMES = (/   &
                        'Li','Be','C ','O ','Na','Mg','Al','Si','K '                    &
                       ,'Ca','Ti','V ','Cr','Mn','Fe','Co','Ni','Cu','Zn','Ga','Ge'     &
                       ,'As','Rb','Sr','Y ','Zr','Nb','Mo','Ru','Rh','Pd','Ag','Cd'     &
                       ,'Sn','Cs','Ba','La','Ce','Gd','Hf','Ta','W ','Re','Os','Ir'     &
                       ,'Pt','Au','Tl','Pb','Th' /)        

        character(len=LATTNAME_LEN),dimension(NELNAMES),parameter  ::      ELLATTICE = (/           &
                        'bcc','hcp','dia','unk','bcc','hcp','fcc','dia','bcc'                       &
                       ,'fcc','hcp','bcc','bcc','bcc','bcc','hcp','fcc','fcc','fcc','ort','dia'     &
                       ,'rho','bcc','fcc','hcp','hcp','bcc','bcc','hcp','fcc','fcc','fcc','hcp'     &
                       ,'tet','bcc','bcc','dhc','dhc','hcp','hcp','bcc','bcc','hcp','hcp','fcc'     &
                       ,'fcc','fcc','hcp','fcc','fcc' /)        
                
        include "realDoyleTurnerCoefficients.h"
        include "absorbativeScatteringFactors.h"

    !   https://pubchem.ncbi.nlm.nih.gov/ptable/atomic-mass/                       
        real(kind=real64),dimension(NELNAMES),parameter             ::      ELMASS = (/                                                         &
                        7.0d0,9.012183d0,12.011d0,15.999d0,22.9897693d0,24.305d0,26.981538d0,28.085d0,39.0983d0                                 &
                       ,40.08d0,47.867d0,50.9415d0,51.996d0,54.93804d0,55.84d0,58.93319d0,58.693d0,63.55d0,65.4d0,69.723d0,72.63d0              &
                       ,74.92159d0,85.468d0,87.62d0,88.90584d0,91.22d0,92.90637d0,95.95d0,101.1d0,102.9055d0,106.42d0,107.868d0,112.41d0        &
                       ,118.71d0,132.9054520d0,137.33d0,138.9055d0,140.116d0,157.2d0,178.49d0,180.9479d0,183.84d0,186.207d0,190.2d0,192.22d0    &
                       ,195.08d0,196.96657d0,204.383d0,207d0,232.038d0  /)

    !---    lattice parameter (pm)       
    !   https://periodictable.com/Properties/A/LatticeConstants.html                                   
        real(kind=real64),dimension(NELNAMES),parameter             ::      ELA0 = (/                                                           &
                        351.00d0,228.58d0,246.40d0,540.30d0,429.06d0,320.94d0,404.95d0,543.09d0,532.80d0                                        &
                       ,558.84d0,295.08d0,303.00d0,291.00d0,891.25d0,286.65d0,250.71d0,352.40d0,361.49d0,266.49d0,451.97d0,565.75d0             &
                       ,375.98d0,558.50d0,608.49d0,364.74d0,323.20d0,330.04d0,314.70d0,270.59d0,380.34d0,389.07d0,408.53d0,297.94d0             &
                       ,583.18d0,614.10d0,502.80d0,377.20d0,362.00d0,363.60d0,319.64d0,330.13d0,316.52d0,276.10d0,273.44d0,383.90d0             &
                       ,392.42d0,407.82d0,345.66d0,495.08d0,508.42d0                  /)

        real(kind=real64),dimension(NELNAMES),parameter             ::      ELB0 = (/                                                           &
                        351.00d0,228.58d0,246.40d0,342.90d0,429.06d0,320.94d0,404.95d0,543.09d0,532.80d0                                        &
                       ,558.84d0,295.08d0,303.00d0,291.00d0,891.25d0,286.65d0,250.71d0,352.40d0,361.49d0,266.49d0,766.33d0,565.75d0             &
                       ,375.98d0,558.50d0,608.49d0,364.74d0,323.20d0,330.04d0,314.70d0,270.59d0,380.34d0,389.07d0,408.53d0,297.94d0             &
                       ,583.18d0,614.10d0,502.80d0,377.20d0,362.00d0,363.60d0,319.64d0,330.13d0,316.52d0,276.10d0,273.44d0,383.90d0             &
                       ,392.42d0,407.82d0,345.66d0,495.08d0,508.42d0                  /)

        real(kind=real64),dimension(NELNAMES),parameter             ::      ELC0 = (/                                                           &
                        351.00d0,358.43d0,671.10d0,508.60d0,429.06d0,521.08d0,404.95d0,543.09d0,532.80d0                                        &
                       ,558.84d0,468.55d0,303.00d0,291.00d0,891.25d0,286.65d0,406.95d0,352.40d0,361.49d0,494.68d0,452.60d0,565.75d0             &
                      ,1054.75d0,558.50d0,608.49d0,573.06d0,514.70d0,330.04d0,314.70d0,428.15d0,380.34d0,389.07d0,408.53d0,561.86d0             &
                       ,318.19d0,614.10d0,502.80d0,1214.4d0,599.00d0,578.26d0,505.11d0,330.13d0,316.52d0,445.60d0,431.73d0,383.90d0             &
                       ,392.42d0,407.82d0,552.48d0,495.08d0,508.42d0                  /)


        public          ::      numberOfElements        
        public          ::      getElementName        
        public          ::      getLatticeName
        public          ::      whichElement
        public          ::      DebyeWallerFactor
        public          ::      getLatticeConstant
        public          ::      getElementMass
        public          ::      getrealDoyleTurner_a
        public          ::      getrealDoyleTurner_b
        public          ::      getscatDoyleTurner_a
        public          ::      getscatDoyleTurner_b


        interface       getElementMass
            module procedure    getElementMass0
        end interface



        interface getLatticeName
            module procedure        getLatticeName0
        end interface
        


    contains
!---^^^^^^^^
    
        pure integer function numberOfElements()
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
            numberOfElements = NELNAMES
            return
        end function numberOfElements

        pure character(len=ELNAME_LEN) function getElementName(i)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
            integer,intent(in)              ::      i
            getElementName = ""
            if ( i*(NELNAMES+1-i) > 0 ) getElementName = ELNAMES(i)
            return
        end function getElementName
            


    
        pure integer function whichElement(el)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      return the element code for the element.
    !*      note: this is not the atomic number or mass, but the position in the ELNAMES table
            character(len=*),intent(in)         ::      el
            character(len=len(el))              ::      el_uc
            integer         ::      ii
            el_uc = convertUpperCase(el) 
            do ii = 1,NELNAMES
                !print *,ii,""""//trim(el_uc)//""" """//trim(convertUpperCase(ELNAMES(ii)) )//""" ",(trim(el_uc) == trim(convertUpperCase(ELNAMES(ii)) ))
                if (trim(el_uc) == trim(convertUpperCase(ELNAMES(ii)) )) then
                    whichElement = ii
                    return
                end if
            end do
            whichElement = 0
            return
        end function whichElement





        pure character(len=LATTNAME_LEN) function getLatticeName0(el) 
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      return name of lattice recorded for element
            character(len=*),intent(in)         ::      el
            integer         ::      ii
            ii = whichElement(el)
            if (ii /= 0) then
                getLatticeName0 = ELLATTICE(ii)
            else
                getLatticeName0 = "unk"
            end if
            return
        end function getLatticeName0


        pure real(kind=real64) function getElementMass0(el) 
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      return mass recorded for element
            character(len=*),intent(in)         ::      el
            integer         ::      ii
            ii = whichElement(el)
            if (ii /= 0) then
                getElementMass0 = ELMASS(ii)
            else
                getElementMass0 = 0.0d0
            end if
            return
        end function getElementMass0



        pure function getLatticeConstant(el) result(a0)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      return lattice constants for element in A
            character(len=*),intent(in)         ::      el
            real(kind=real64),dimension(3)      ::      a0            
            integer         ::      ii
            ii = whichElement(el)
            if (ii /= 0) then
                a0 = (/ ELA0(ii),ELB0(ii),ELC0(ii) /) / 100.0d0
            else
                a0 = 1.0d0
            end if
            return
        end function getLatticeConstant


        real(kind=real64) function DebyeWallerFactor(el,T,latt)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      try to return the debye waller factor for this element at temperature T
    !*      optionally force the crystal structure to be latt

            character(len=*),intent(in)                 ::      el
            real(kind=real64),intent(in)                ::      T
            character(len=*),intent(in),optional        ::      latt

        !---    elements recognised
            integer,parameter       ::      NELEMENTS = 42
            character(len=5),dimension(NELEMENTS),parameter     ::      ELEMENT = (/                &
                                                                "LI   ","BE   ","C    ","NA   ",    &     
                                                                "MG   ","A1   ","SI   ","K    ",    &  
                                                                "CA   ","CABCC","SC   ","TI   ",    &  
                                                                "V    ","CR   ","FE   ","FEFCC",    &
                                                                "NI   ","CU   ","ZN   ","GE   ",    &  
                                                                "KR   ","RB   ","SR   ","Y    ",    &  
                                                                "ZR   ","NB   ","MO   ","PD   ",    &
                                                                "AG   ","SN   ","XE   ","CS   ",    & 
                                                                "BA   ","LA   ","TB   ","HO   ",    &  
                                                                "TA   ","W    ","PT   ","AU   ",    &
                                                                "PB   ","TH   "                  /)
            
        !---    temperature range tabulated
            integer,parameter       ::      NTEMPERATURES = 17
            real(kind=real64),dimension(NTEMPERATURES),parameter           ::      TEMPERATURE = (/            &
                                                                85.0d0 ,                           &
                                                                90.0d0 ,                           &
                                                                110.0d0 ,                           &
                                                                120.0d0 ,                           &
                                                                130.0d0 ,                           &
                                                                140.0d0 ,                           &
                                                                160.0d0 ,                           &
                                                                170.0d0 ,                           &
                                                                180.0d0 ,                           &
                                                                190.0d0 ,                           &
                                                                210.0d0 ,                           &
                                                                220.0d0 ,                           &
                                                                230.0d0 ,                           &
                                                                240.0d0 ,                           &
                                                                260.0d0 ,                           &
                                                                270.0d0 ,                           &
                                                                280.0d0     /)                                                            
            
                
            real(kind=real64),dimension(NTEMPERATURES,NELEMENTS),parameter  ::      DWF_TABLE = reshape( (/          &                                                            
                    1.9009d0, 1.4780d0, 1.8063d0, 1.9705d0, 2.1347d0, 2.2989d0, 2.6272d0, 2.7913d0, 2.9554d0, 3.1195d0, 3.4477d0, 3.6117d0, 3.7758d0, 3.9398d0, 4.2678d0, 4.4318d0, 4.5950d0,       &
                    0.3637d0, 0.3652d0, 0.3719d0, 0.3758d0, 0.3800d0, 0.3846d0, 0.3947d0, 0.4002d0, 0.4061d0, 0.4123d0, 0.3032d0, 0.3176d0, 0.3321d0, 0.3465d0, 0.3754d0, 0.3898d0, 0.4040d0,       &
                    0.1299d0, 0.1300d0, 0.1307d0, 0.1311d0, 0.1315d0, 0.1320d0, 0.1330d0, 0.1336d0, 0.1342d0, 0.1349d0, 0.1363d0, 0.1370d0, 0.1378d0, 0.1386d0, 0.1403d0, 0.1412d0, 0.1420d0,       &
                    1.9377d0, 2.0516d0, 2.5072d0, 2.7348d0, 2.9624d0, 3.1900d0, 3.6448d0, 3.8721d0, 4.0993d0, 4.3264d0, 4.7802d0, 5.0070d0, 5.2336d0, 5.4601d0, 5.9128d0, 6.1389d0, 6.3640d0,       &
                    0.5263d0, 0.5572d0, 0.6810d0, 0.7429d0, 0.8048d0, 0.8667d0, 0.9904d0, 1.0523d0, 1.1141d0, 1.1759d0, 1.2996d0, 1.3614d0, 1.4232d0, 1.4849d0, 1.6085d0, 1.6702d0, 1.7310d0,       &
                    0.3301d0, 0.3374d0, 0.2932d0, 0.3198d0, 0.3465d0, 0.3731d0, 0.4264d0, 0.4531d0, 0.4797d0, 0.5064d0, 0.5596d0, 0.5863d0, 0.6129d0, 0.6395d0, 0.6928d0, 0.7194d0, 0.7460d0,       &
                    0.2170d0, 0.2201d0, 0.2342d0, 0.2423d0, 0.2511d0, 0.2607d0, 0.2648d0, 0.2814d0, 0.2979d0, 0.3145d0, 0.3476d0, 0.3641d0, 0.3807d0, 0.3972d0, 0.4303d0, 0.4469d0, 0.4630d0,       &
                    3.1319d0, 3.3159d0, 4.0510d0, 4.4182d0, 4.7852d0, 5.1518d0, 5.8840d0, 6.2496d0, 6.6148d0, 6.9796d0, 7.7077d0, 8.0711d0, 8.4340d0, 8.7963d0, 9.5192d0, 9.8798d0,10.2390d0,       &
                    0.5805d0, 0.6146d0, 0.7511d0, 0.8194d0, 0.8876d0, 0.9558d0, 1.0922d0, 1.1603d0, 1.2285d0, 1.2966d0, 1.4328d0, 1.5008d0, 1.5688d0, 1.6368d0, 1.7728d0, 1.8407d0, 1.9080d0,       &
                    0.7556d0, 0.8000d0, 0.9777d0, 1.0665d0, 1.1553d0, 1.2441d0, 1.4216d0, 1.5104d0, 1.5991d0, 1.6878d0, 1.8651d0, 1.9537d0, 2.0423d0, 2.1308d0, 2.3078d0, 2.3963d0, 2.4840d0,       &
                    0.2165d0, 0.2292d0, 0.2801d0, 0.3056d0, 0.3310d0, 0.3565d0, 0.4074d0, 0.4328d0, 0.4583d0, 0.4837d0, 0.5346d0, 0.5600d0, 0.5854d0, 0.6109d0, 0.6617d0, 0.6871d0, 0.7120d0,       &
                    0.1491d0, 0.1579d0, 0.1930d0, 0.2105d0, 0.2281d0, 0.2456d0, 0.2807d0, 0.2982d0, 0.3158d0, 0.3333d0, 0.3684d0, 0.3859d0, 0.4034d0, 0.4210d0, 0.4560d0, 0.4735d0, 0.4910d0,       &
                    0.1660d0, 0.1758d0, 0.2148d0, 0.2343d0, 0.2539d0, 0.2734d0, 0.3124d0, 0.3320d0, 0.3515d0, 0.3710d0, 0.4100d0, 0.4295d0, 0.4490d0, 0.4685d0, 0.5075d0, 0.5270d0, 0.5460d0,       &
                    0.1259d0, 0.1274d0, 0.0943d0, 0.1028d0, 0.1114d0, 0.1200d0, 0.1371d0, 0.1457d0, 0.1542d0, 0.1628d0, 0.1799d0, 0.1885d0, 0.1970d0, 0.2056d0, 0.2227d0, 0.2313d0, 0.2390d0,       &
                    0.1461d0, 0.1493d0, 0.1221d0, 0.1332d0, 0.1443d0, 0.1554d0, 0.1775d0, 0.1886d0, 0.1997d0, 0.2108d0, 0.2330d0, 0.2441d0, 0.2552d0, 0.2663d0, 0.2884d0, 0.2995d0, 0.3100d0,       &
                    0.1619d0, 0.1715d0, 0.2095d0, 0.2286d0, 0.2476d0, 0.2667d0, 0.3047d0, 0.3238d0, 0.3428d0, 0.3618d0, 0.3999d0, 0.4189d0, 0.4379d0, 0.4570d0, 0.4950d0, 0.5140d0, 0.5330d0,       &
                    0.1502d0, 0.1094d0, 0.1337d0, 0.1458d0, 0.1580d0, 0.1701d0, 0.1944d0, 0.2066d0, 0.2187d0, 0.2309d0, 0.2551d0, 0.2673d0, 0.2794d0, 0.2916d0, 0.3158d0, 0.3280d0, 0.3400d0,       &
                    0.1598d0, 0.1692d0, 0.2068d0, 0.2256d0, 0.2444d0, 0.2632d0, 0.3008d0, 0.3196d0, 0.3384d0, 0.3571d0, 0.3947d0, 0.4135d0, 0.4322d0, 0.4510d0, 0.4886d0, 0.5073d0, 0.5260d0,       &
                    0.3317d0, 0.3512d0, 0.4292d0, 0.4682d0, 0.5072d0, 0.5462d0, 0.6242d0, 0.6632d0, 0.7022d0, 0.7411d0, 0.8190d0, 0.8580d0, 0.8969d0, 0.9359d0, 1.0137d0, 1.0526d0, 1.0910d0,       &
                    0.1670d0, 0.1844d0, 0.2254d0, 0.2458d0, 0.2663d0, 0.2868d0, 0.3278d0, 0.3483d0, 0.3687d0, 0.3892d0, 0.4302d0, 0.4506d0, 0.4711d0, 0.4915d0, 0.5325d0, 0.5529d0, 0.5730d0,       &
                    2.8102d0, 2.9748d0, 3.6321d0, 0.0000d0, 0.0000d0, 0.0000d0, 0.0000d0, 0.0000d0, 0.0000d0, 0.0000d0, 0.0000d0, 0.0000d0, 0.0000d0, 0.0000d0, 0.0000d0, 0.0000d0, 0.0000d0,       &
                    3.8986d0, 4.1269d0, 5.0388d0, 5.4935d0, 5.9475d0, 6.4004d0, 7.3033d0, 7.7531d0, 8.2016d0, 8.6488d0, 9.5392d0, 9.9821d0,10.4235d0,10.8632d0,11.7373d0,12.1717d0,12.6040d0,       &
                    1.1061d0, 1.1710d0, 1.4308d0, 1.5606d0, 1.6903d0, 1.8199d0, 2.0788d0, 2.2082d0, 2.3374d0, 2.4665d0, 2.7243d0, 2.8530d0, 2.9816d0, 3.1100d0, 3.3664d0, 3.4944d0, 3.6220d0,       &
                    0.2486d0, 0.2633d0, 0.3217d0, 0.3510d0, 0.3802d0, 0.4094d0, 0.4678d0, 0.4970d0, 0.5262d0, 0.5554d0, 0.6137d0, 0.6429d0, 0.6720d0, 0.7011d0, 0.7594d0, 0.7885d0, 0.8170d0,       &
                    0.1652d0, 0.1749d0, 0.2137d0, 0.2331d0, 0.2526d0, 0.2720d0, 0.3108d0, 0.3302d0, 0.3496d0, 0.3690d0, 0.4078d0, 0.4272d0, 0.4466d0, 0.4660d0, 0.5047d0, 0.5241d0, 0.5430d0,       &
                    0.1306d0, 0.1383d0, 0.1690d0, 0.1843d0, 0.1997d0, 0.2150d0, 0.2457d0, 0.2611d0, 0.2764d0, 0.2918d0, 0.3225d0, 0.3378d0, 0.3531d0, 0.3685d0, 0.3991d0, 0.4144d0, 0.4290d0,       &
                    0.0626d0, 0.0662d0, 0.0809d0, 0.0883d0, 0.0957d0, 0.1030d0, 0.1177d0, 0.1251d0, 0.1324d0, 0.1398d0, 0.1545d0, 0.1619d0, 0.1692d0, 0.1766d0, 0.1912d0, 0.1986d0, 0.2050d0,       &
                    0.1294d0, 0.1370d0, 0.1674d0, 0.1827d0, 0.1979d0, 0.2131d0, 0.2435d0, 0.2587d0, 0.2739d0, 0.2891d0, 0.3196d0, 0.3347d0, 0.3499d0, 0.3651d0, 0.3955d0, 0.4107d0, 0.4250d0,       &
                    0.2134d0, 0.2259d0, 0.2761d0, 0.3012d0, 0.3262d0, 0.3513d0, 0.4015d0, 0.4265d0, 0.4516d0, 0.4766d0, 0.5267d0, 0.5517d0, 0.5767d0, 0.6017d0, 0.6517d0, 0.6767d0, 0.7017d0,       &
                    0.3290d0, 0.3483d0, 0.4257d0, 0.4644d0, 0.5031d0, 0.5418d0, 0.6191d0, 0.6577d0, 0.6964d0, 0.7350d0, 0.8123d0, 0.8509d0, 0.8895d0, 0.9281d0, 1.0053d0, 1.0439d0, 1.0825d0,       &
                    2.5987d0, 2.7507d0, 3.3574d0, 3.6598d0, 3.9614d0, 4.2622d0, 4.8611d0, 0.0000d0, 0.0000d0, 0.0000d0, 0.0000d0, 0.0000d0, 0.0000d0, 0.0000d0, 0.0000d0, 0.0000d0, 0.0000d0,       &
                    5.2797d0, 5.5878d0, 6.8158d0, 7.4268d0, 8.0355d0, 8.6418d0, 9.8464d0,10.4443d0,11.0390d0,11.6304d0,12.8022d0,13.3823d0,13.9583d0,14.5300d0,15.6596d0,16.2173d0,16.7698d0,       &
                    0.8961d0, 0.9488d0, 1.1591d0, 1.2641d0, 1.3691d0, 1.4739d0, 1.6833d0, 1.7878d0, 1.8922d0, 1.9965d0, 2.2046d0, 2.3085d0, 2.4121d0, 2.5157d0, 2.7222d0, 2.8251d0, 2.9279d0,       &
                    0.5419d0, 0.5738d0, 0.7010d0, 0.7646d0, 0.8281d0, 0.8916d0, 1.0185d0, 1.0818d0, 1.1451d0, 1.2083d0, 1.3345d0, 1.3975d0, 1.4605d0, 1.5233d0, 1.6488d0, 1.7113d0, 1.7739d0,       &
                    0.2944d0, 0.3117d0, 0.3809d0, 0.4155d0, 0.4500d0, 0.4846d0, 0.5536d0, 0.5881d0, 0.6226d0, 0.6570d0, 0.7259d0, 0.7603d0, 0.7946d0, 0.8290d0, 0.8975d0, 0.9318d0, 0.9660d0,       &
                    0.2444d0, 0.2588d0, 0.3162d0, 0.3449d0, 0.3737d0, 0.4023d0, 0.4597d0, 0.4883d0, 0.5170d0, 0.5456d0, 0.6028d0, 0.6313d0, 0.6599d0, 0.6884d0, 0.7454d0, 0.7739d0, 0.8023d0,       &
                    0.0936d0, 0.0991d0, 0.1211d0, 0.1321d0, 0.1431d0, 0.1541d0, 0.1761d0, 0.1871d0, 0.1981d0, 0.2090d0, 0.2310d0, 0.2420d0, 0.2529d0, 0.2639d0, 0.2858d0, 0.2968d0, 0.3078d0,       &
                    0.0464d0, 0.0491d0, 0.0600d0, 0.0654d0, 0.0709d0, 0.0763d0, 0.0872d0, 0.0927d0, 0.0981d0, 0.1036d0, 0.1145d0, 0.1199d0, 0.1254d0, 0.1308d0, 0.1417d0, 0.1471d0, 0.1526d0,       &
                    0.1082d0, 0.1145d0, 0.1400d0, 0.1527d0, 0.1654d0, 0.1781d0, 0.2035d0, 0.2162d0, 0.2289d0, 0.2416d0, 0.2670d0, 0.2797d0, 0.2924d0, 0.3051d0, 0.3304d0, 0.3431d0, 0.3557d0,       &
                    0.1802d0, 0.1908d0, 0.2332d0, 0.2544d0, 0.2755d0, 0.2967d0, 0.3390d0, 0.3602d0, 0.3813d0, 0.4025d0, 0.4448d0, 0.4659d0, 0.4870d0, 0.5081d0, 0.5503d0, 0.5714d0, 0.5925d0,       &
                    0.6228d0, 0.6593d0, 0.8055d0, 0.8784d0, 0.9514d0, 1.0242d0, 1.1697d0, 1.2423d0, 1.3148d0, 1.3873d0, 1.5318d0, 1.6040d0, 1.6760d0, 1.7478d0, 1.8912d0, 1.9627d0, 2.0341d0,       &
                    0.2125d0, 0.2250d0, 0.2750d0, 0.2999d0, 0.3249d0, 0.3498d0, 0.3997d0, 0.4246d0, 0.4495d0, 0.4744d0, 0.5241d0, 0.5490d0, 0.5738d0, 0.5986d0, 0.6482d0, 0.6730d0, 0.6977d0        &
                                                                                            /),(/ NTEMPERATURES,NELEMENTS /) )
                                                                
            
            character(len=len(el)+3)              ::      el_uc                                                                                
            integer                 ::      ii,ee   
            real(kind=real64)       ::      xx
                                                 
            if (present(latt)) then
                if ( (trim(convertUpperCase(el))=="FE") .and. (trim(convertUpperCase(latt))=="BCC") ) then
                    el_uc = "FE"
                else if ( (trim(convertUpperCase(el))=="CA") .and. (trim(convertUpperCase(latt))=="FCC") ) then
                    el_uc = "CA"
                else
                    el_uc = trim(convertUpperCase(el)) 
                end if
            else
                el_uc = trim(convertUpperCase(el))
            end if
            ee = 0
            do ii = 1,NELEMENTS
                if (trim(el_uc) == trim(convertUpperCase(ELEMENT(ii)) )) then
                    ee = ii
                    exit
                end if
            end do

            if (ee==0) then
                print *,"Lib_Elements::DebyeWallerFactor error - did not recognise element """//trim(el_uc)//""""
                print *,"   tabulated only for"
                print *,(""""//trim(ELEMENT(ii))//""" ",ii=1,NELEMENTS)
            end if

            DebyeWallerFactor = 0    
        !---    are we under , inside, over the known temperature block?
            if (T < TEMPERATURE(1)) then
                !   under
                DebyeWallerFactor = 0                   !   could use Debye function to find this?
            else if (T >= TEMPERATURE(NTEMPERATURES)) then
                !   over
                DebyeWallerFactor = DWF_TABLE(NTEMPERATURES,ee) * T / TEMPERATURE(NTEMPERATURES)           
            else 
                !   inside
                do ii = 1,NTEMPERATURES-1
                    if ( (T>=TEMPERATURE(ii)).and.(T<TEMPERATURE(ii+1)) ) then
                        !   bounded by ii,ii+1
                        xx = ( T - TEMPERATURE(ii) )/( TEMPERATURE(ii+1) - TEMPERATURE(ii) )        !   how far across interval are we?    
                        DebyeWallerFactor = DWF_TABLE(ii,ee)*(1-xx) + DWF_TABLE(ii+1,ee)*xx
                        exit
                    end if 
                end do
            end if

            return
        end function DebyeWallerFactor



        function getrealDoyleTurner_a(el) result(a)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      return the real Doyle-Turner coefficients a(1:5) for element el
            character(len=*),intent(in)         ::      el
            real(kind=real64),dimension(5)      ::      a 
            integer         ::      ee

            a = 0
            ee = whichElement(el)
            if (ee>0) a(1:5) = realDoyleTurner_a(1:5,ee)

            return
        end function getrealDoyleTurner_a


        function getrealDoyleTurner_b(el) result(b)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      return the imaginary Doyle-Turner coefficients b(1:5) for element el
            character(len=*),intent(in)         ::      el
            real(kind=real64),dimension(5)      ::      b

            integer         ::      ee

            b = 0
            ee = whichElement(el)
            if (ee>0) b(1:5) = realDoyleTurner_b(1:5,ee)
            

            return
        end function getrealDoyleTurner_b        


        
        function getscatDoyleTurner_a(el,T,latt) result(a)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      return the absorbative scattering factor Doyle-Turner coefficients a(1:5) for element el
    !*      optionally force the crystal structure to be latt

            character(len=*),intent(in)                 ::      el
            real(kind=real64),intent(in)                ::      T
            character(len=*),intent(in),optional        ::      latt
            real(kind=real64),dimension(5)              ::      a 

            character(len=len(el)+3)                    ::      el_uc           
            real(kind=real64),dimension(:,:),pointer    ::      absScatFact
            real(kind=real64),dimension(:,:),pointer    ::      DT_a
            real(kind=real64),dimension(:),pointer      ::      Temperature
            integer,parameter                           ::      TCOL = 1
            integer,parameter                           ::      ACOL = 3
            integer                                     ::      ii,nTemp
            real(kind=real64)                           ::      xx

        !---    find the name of the scattering table in conventional format            
            if (present(latt)) then
                if ( (trim(convertUpperCase(el))=="FE") .and. (trim(convertUpperCase(latt))=="BCC") ) then
                    el_uc = "FE"
                else if ( (trim(convertUpperCase(el))=="CA") .and. (trim(convertUpperCase(latt))=="FCC") ) then
                    el_uc = "CA"
                else
                    el_uc = trim(convertUpperCase(el))//trim(convertUpperCase(latt))
                end if
            else
                el_uc = trim(convertUpperCase(el))
            end if

        !---    find the scattering table
            
            select case(trim(el_uc))
                case("LI")    ;  absScatFact => AbsScatFact_Li
                case("BE")    ;  absScatFact => AbsScatFact_Be
                case("C")     ;  absScatFact => AbsScatFact_C 
                case("NA")    ;  absScatFact => AbsScatFact_Na
                case("MG")    ;  absScatFact => AbsScatFact_Mg
                case("AL")    ;  absScatFact => AbsScatFact_Al
                case("SI")    ;  absScatFact => AbsScatFact_Si
                case("AR")    ;  absScatFact => AbsScatFact_Ar
                case("K")     ;  absScatFact => AbsScatFact_K
                case("CAFCC") ;  absScatFact => AbsScatFact_CaFcc
                case("CA")    ;  absScatFact => AbsScatFact_CaFcc
                case("CABCC") ;  absScatFact => AbsScatFact_CaBcc
                case("SC")    ;  absScatFact => AbsScatFact_Sc
                case("TI")    ;  absScatFact => AbsScatFact_Ti
                case("V")     ;  absScatFact => AbsScatFact_V
                case("CR")    ;  absScatFact => AbsScatFact_Cr
                case("FEBCC") ;  absScatFact => AbsScatFact_FeBcc
                case("FE")    ;  absScatFact => AbsScatFact_FeBcc
                case("FEFCC") ;  absScatFact => AbsScatFact_FeFcc
                case("NI")    ;  absScatFact => AbsScatFact_Ni
                case("CU")    ;  absScatFact => AbsScatFact_Cu
                case("ZN")    ;  absScatFact => AbsScatFact_Zn
                case("GE")    ;  absScatFact => AbsScatFact_Ge
                case("KR")    ;  absScatFact => AbsScatFact_Kr
                case("RB")    ;  absScatFact => AbsScatFact_Rb
                case("SR")    ;  absScatFact => AbsScatFact_Sr
                case("Y")     ;  absScatFact => AbsScatFact_Y
                case("ZR")    ;  absScatFact => AbsScatFact_Zr
                case("NB")    ;  absScatFact => AbsScatFact_Nb
                case("MO")    ;  absScatFact => AbsScatFact_Mo
                case("PD")    ;  absScatFact => AbsScatFact_Pd
                case("AG")    ;  absScatFact => AbsScatFact_Ag
                case("SN")    ;  absScatFact => AbsScatFact_Sn
                case("CS")    ;  absScatFact => AbsScatFact_Cs
                case("BA")    ;  absScatFact => AbsScatFact_Ba
                case("LA")    ;  absScatFact => AbsScatFact_La
                case("TB")    ;  absScatFact => AbsScatFact_Tb
                case("HO")    ;  absScatFact => AbsScatFact_Ho
                case("TA")    ;  absScatFact => AbsScatFact_Ta
                case("W")     ;  absScatFact => AbsScatFact_W
                case("PT")    ;  absScatFact => AbsScatFact_Pt
                case("AU")    ;  absScatFact => AbsScatFact_Au
                case("PB")    ;  absScatFact => AbsScatFact_Pb
                case("TH")    ;  absScatFact => AbsScatFact_Th
                case default
                    print *,"Lib_Elements::getscatDoyleTurner_a ERROR - did not recognise element """//trim(el_uc)//""""
                    a = 0
                    return                    
            end select


            Temperature => absScatFact(:,TCOL)
            DT_a => absScatFact(:,ACOL:ACOL+4)


        !---    interpolate correct temperature
            nTemp = size(Temperature,dim=1)
            a(1:5) = DT_a(1,:)
        !---    are we under , inside, over the known temperature block?
            if (T < Temperature(1)) then
                !   under
                a(1:5) = DT_a(1,:)
            else if (T >= Temperature(nTemp)) then
                !   over
                a(1:5) = DT_a(nTemp,:)          
            else 
                !   inside
                do ii = 1,nTemp-1
                    if ( (T>=Temperature(ii)).and.(T<Temperature(ii+1)) ) then
                        !   bounded by ii,ii+1
                        xx = ( T - Temperature(ii) )/( Temperature(ii+1) - Temperature(ii) )        !   how far across interval are we?    
                        a(1:5) = DT_a(ii,:)*(1-xx) + DT_a(ii+1,:)*xx
                        !print *,"T,Temperature(ii),xx,DT_a(ii,:)",T,Temperature(ii),xx,DT_a(ii,:)
                        exit
                    end if 
                end do
            end if


            return
        end function getscatDoyleTurner_a

        function getscatDoyleTurner_b(el,T,latt) result(b)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      return the absorbative scattering factor imaginary Doyle-Turner coefficients a(1:5) for element el
    !*      optionally force the crystal structure to be latt

            character(len=*),intent(in)                 ::      el
            real(kind=real64),intent(in)                ::      T
            character(len=*),intent(in),optional        ::      latt
            real(kind=real64),dimension(5)              ::      b 

            character(len=len(el)+3)                    ::      el_uc           
            real(kind=real64),dimension(:,:),pointer    ::      absScatFact
            real(kind=real64),dimension(:,:),pointer    ::      DT_b
            real(kind=real64),dimension(:),pointer      ::      Temperature
            integer,parameter                           ::      TCOL = 1
            integer,parameter                           ::      BCOL = 8
            integer                                     ::      ii,nTemp
            real(kind=real64)                           ::      xx

        !---    find the name of the scattering table in conventional format            
            if (present(latt)) then
                if ( (trim(convertUpperCase(el))=="FE") .and. (trim(convertUpperCase(latt))=="BCC") ) then
                    el_uc = "FE"
                else if ( (trim(convertUpperCase(el))=="CA") .and. (trim(convertUpperCase(latt))=="FCC") ) then
                    el_uc = "CA"
                else
                    el_uc = trim(convertUpperCase(el))//trim(convertUpperCase(latt))
                end if
            else
                el_uc = trim(convertUpperCase(el))
            end if

        !---    find the scattering table
            select case(trim(el_uc))
                case("LI")    ;  absScatFact => AbsScatFact_Li
                case("BE")    ;  absScatFact => AbsScatFact_Be
                case("C")     ;  absScatFact => AbsScatFact_C 
                case("NA")    ;  absScatFact => AbsScatFact_Na
                case("MG")    ;  absScatFact => AbsScatFact_Mg
                case("AL")    ;  absScatFact => AbsScatFact_Al
                case("SI")    ;  absScatFact => AbsScatFact_Si
                case("AR")    ;  absScatFact => AbsScatFact_Ar
                case("K")     ;  absScatFact => AbsScatFact_K
                case("CAFCC") ;  absScatFact => AbsScatFact_CaFcc
                case("CA")    ;  absScatFact => AbsScatFact_CaFcc
                case("CABCC") ;  absScatFact => AbsScatFact_CaBcc
                case("SC")    ;  absScatFact => AbsScatFact_Sc
                case("TI")    ;  absScatFact => AbsScatFact_Ti
                case("V")     ;  absScatFact => AbsScatFact_V
                case("CR")    ;  absScatFact => AbsScatFact_Cr
                case("FEBCC") ;  absScatFact => AbsScatFact_FeBcc
                case("FE")    ;  absScatFact => AbsScatFact_FeBcc
                case("FEFCC") ;  absScatFact => AbsScatFact_FeFcc
                case("NI")    ;  absScatFact => AbsScatFact_Ni
                case("CU")    ;  absScatFact => AbsScatFact_Cu
                case("ZN")    ;  absScatFact => AbsScatFact_Zn
                case("GE")    ;  absScatFact => AbsScatFact_Ge
                case("KR")    ;  absScatFact => AbsScatFact_Kr
                case("RB")    ;  absScatFact => AbsScatFact_Rb
                case("SR")    ;  absScatFact => AbsScatFact_Sr
                case("Y")     ;  absScatFact => AbsScatFact_Y
                case("ZR")    ;  absScatFact => AbsScatFact_Zr
                case("NB")    ;  absScatFact => AbsScatFact_Nb
                case("MO")    ;  absScatFact => AbsScatFact_Mo
                case("PD")    ;  absScatFact => AbsScatFact_Pd
                case("AG")    ;  absScatFact => AbsScatFact_Ag
                case("SN")    ;  absScatFact => AbsScatFact_Sn
                case("CS")    ;  absScatFact => AbsScatFact_Cs
                case("BA")    ;  absScatFact => AbsScatFact_Ba
                case("LA")    ;  absScatFact => AbsScatFact_La
                case("TB")    ;  absScatFact => AbsScatFact_Tb
                case("HO")    ;  absScatFact => AbsScatFact_Ho
                case("TA")    ;  absScatFact => AbsScatFact_Ta
                case("W")     ;  absScatFact => AbsScatFact_W
                case("PT")    ;  absScatFact => AbsScatFact_Pt
                case("AU")    ;  absScatFact => AbsScatFact_Au
                case("PB")    ;  absScatFact => AbsScatFact_Pb
                case("TH")    ;  absScatFact => AbsScatFact_Th
                case default
                    print *,"Lib_Elements::getscatDoyleTurner_b ERROR - did not recognise element """//trim(el_uc)//""""
                    b = 0
                    return                    
            end select


            Temperature => absScatFact(:,TCOL)
            DT_b => absScatFact(:,BCOL:BCOL+4)


        !---    interpolate correct temperature
            nTemp = size(Temperature,dim=1)
        !---    are we under , inside, over the known temperature block?
            b(1:5) = DT_b(1,:)
            if (T < Temperature(1)) then
                !   under
                b(1:5) = DT_b(1,:)
            else if (T >= Temperature(nTemp)) then
                !   over
                b(1:5) = DT_b(nTemp,:)          
            else 
                !   inside
                do ii = 1,nTemp-1
                    if ( (T>=Temperature(ii)).and.(T<Temperature(ii+1)) ) then
                        !   bounded by ii,ii+1
                        xx = ( T - Temperature(ii) )/( Temperature(ii+1) - Temperature(ii) )        !   how far across interval are we?    
                        b(1:5) = DT_b(ii,:)*(1-xx) + DT_b(ii+1,:)*xx
                        exit
                    end if 
                end do
            end if


            return
        end function getscatDoyleTurner_b


        pure function convertUpperCase(text_in) result(text_out)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      convert string to all upper case
            character(len=*),intent(in)                         ::      text_in
            character(len=len(text_in))                         ::      text_out
            integer             ::      ii
            character(len=1)    ::      letter
            text_out = text_in
            do ii = 1,len_trim(text_in)
                letter = text_in(ii:ii)
                if (iachar(letter)>96) text_out(ii:ii) = achar( iachar(letter)-32 )
            end do
            return
        end function convertUpperCase
                 
        
                



    end module Lib_Elements