
    
    program crystalStructureFactors
!---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
        use iso_fortran_env
        use Lib_CommandLineArguments
        use Lib_RelativisticElectrons
        use Lib_CrystalStructureFactor
        implicit none

    !---    version history
    !       v.0.0.1                 Working code used for two beam parameterization
            
        
        
        
    !---    define the elements recognised        
        integer,parameter                               ::      NELNAMES = 50
        character(len=2),dimension(NELNAMES),parameter  ::      ELNAMES = (/            &
                        'LI','BE','C ','O ','NA','MG','AL','SI','K '                    &
                       ,'CA','TI','V ','CR','MN','FE','CO','NI','CU','ZN','GA','GE'     &
                       ,'AS','RB','SR','Y ','ZR','NB','MO','RU','RH','PD','AG','CD'     &
                       ,'SN','CS','BA','LA','CE','GD','HF','TA','W ','RE','OS','IR'     &
                       ,'PT','AU','TL','PB','TH' /)        
        character(len=3),dimension(NELNAMES),parameter  ::      ELLATTICE = (/                      &
                        'bcc','hcp','dia','unk','bcc','hcp','fcc','dia','bcc'                       &
                       ,'fcc','hcp','bcc','bcc','bcc','bcc','hcp','fcc','fcc','fcc','ort','dia'     &
                       ,'rho','bcc','fcc','hcp','hcp','bcc','bcc','hcp','fcc','fcc','fcc','hcp'     &
                       ,'tet','bcc','bcc','dhc','dhc','hcp','hcp','bcc','bcc','hcp','hcp','fcc'     &
                       ,'fcc','fcc','hcp','fcc','fcc' /)        
        
        !real(kind=real64),parameter         ::      HBAR = 0.6582119569d0       !   Dirac constant eV.fs
            
        real(kind=real64),parameter         ::      PI = 3.14159265390d0        !   Ronseal
        character(len=8),parameter          ::      VERSION = "0.0.1"
                       
                               
        type(CommandLineArguments)          ::      cla
    
        character(len=8)                    ::      el = ""
        real(kind=real64)                   ::      a0 = 1.0d0
        real(kind=real64),dimension(3)      ::      g = (/ 2,0,0 /)
        real(kind=real64)                   ::      V = 150000d0
        real(kind=real64)                   ::      T = 300d0
                
        complex(kind=real64)                ::      U0,Ug        
        real(kind=real64)                   ::      xi0,xig
        character(len=256)                  ::      outfile = "Output/extinctionDistances.dat"          
        
        character(len=2)        ::      elu
        character(len=3)        ::      lattice
        integer                 ::      ii
        character(len=256)      ::      syscall,dummy,fitfile,dummy0,dummyg
        
        real(kind=real64),dimension(5)      ::      DTcoeffar,DTcoeffai,DTcoeffbr,DTcoeffbi
        complex(kind=real64),dimension(5)   ::      aa,bb
        real(kind=real64)                   ::      dwf
        
        logical                 ::      ok,okg
        
        cla = CommandLineArguments_ctor(20) 
        call setProgramDescription( cla, "crystalStructureFactor" )
        call setProgramVersion( cla, "1.1" )
        
        call get( cla,"e",el ,LIB_CLA_REQUIRED,"           element" ) 
        call get( cla,"a0",a0 ,LIB_CLA_REQUIRED,"          unit cell dimension (A)" ) 
        call get( cla,"V",V ,LIB_CLA_REQUIRED,"           electron accelerating voltage (eV)" ) 
        call get( cla,"T",T ,LIB_CLA_REQUIRED,"           temperature (K)" ) 
        ii=3 ; call get( cla,"g",g,ii,LIB_CLA_REQUIRED,"           g-vector (reduced units eg 1,1,0)" )         
        call get( cla,"o",outfile,LIB_CLA_OPTIONAL,"      output extinction distances file")
        
        call report(cla)
        if (.not. allRequiredArgumentsSet(cla)) stop
        if (hasHelpArgument(cla)) stop
        call delete(cla)
        
        
        print *,""
        print *,"crystalStructureFactor v."//trim(VERSION)
        print *,"^^^^^^^^^^^^^^^^^^^^^^^^^"//repeat("^",len_trim(VERSION))
        print *,""
        print *,"   element     """//trim(el)//""""
        print *,"   a0           ",a0," (A)"
        print *,"   V            ",V ," (eV)"
        print *,"   T            ",T ," (K)"
        print *,"   g            ",g
        print *,"   outfile     """//trim(outfile)//""""        
        print *,""
        
        
        
        
    !---    first convert element name to upper case            
        elu(1:2) = adjustl(el)
        if (iachar(elu(1:1))>96) elu(1:1) = achar( iachar(elu(1:1))-32 )
        if (iachar(elu(2:2))>96) elu(2:2) = achar( iachar(elu(2:2))-32 )

        
    !---    is this an acceptable name?
        if (.not. any(ELNAMES == elu) ) then
            print *,"crystalStructureFactor error - did not recognise element name """//trim(el)//""""
            stop
        end if
        
    !---    what is the lattice?
        do ii = 1,NELNAMES
            if (ELNAMES(ii) == elu) then
                lattice = ELLATTICE(ii)
                exit
            end if
        end do
        
        
        fitfile = "Output/"//trim(elu)//".DT.dat"
        print *,"crystalStructureFactor info - writing to file """//trim(fitfile)//""""
        
     
        
        
        ok = .false. ; okg = .false.
        write(dummy0,fmt='(a8,3i4,2f12.3,3f20.6)') elu,0,0,0,T,V ,real(U0),aimag(U0),xi0 
        write(dummyg,fmt='(a8,3i4,2f12.3,3f20.6)') elu,nint(g),T,V ,real(Ug),aimag(Ug),xig 
        open(unit=710,file=trim(outfile),action="read")
            do 
                read(unit=710,fmt='(a44)',iostat=ii) dummy(1:44)
                if (ii/=0) exit
                if (dummy(1:44) == dummy0(1:44)) ok = .true.
                if (dummy(1:44) == dummyg(1:44)) okg = .true.
                if (ok .and. okg) exit                
            end do
        close(unit=710)
        if (.not. ok) print *,"adding result U_0"
        if (.not. okg) print *,"adding result U_g"
        if (ok .and. okg) then
            print *,"results for U_0 and U_g exist already"
            stop
        end if
                
        
        
        
        
        
        
        
        
        
    !---    find the imaginary Doyle-Turner coefficients.            
        syscall = "./bin/run.TDSfit.sh "//trim(elu) 
        write(dummy,fmt='(g16.3)') V ; syscall = trim(syscall)//" "//trim(adjustl(dummy))
        write(dummy,fmt='(g16.3)') T ; syscall = trim(syscall)//" "//trim(adjustl(dummy))//" > tmp"  
        print *,trim(syscall)  
        call system( syscall )    
        
    !   then extract the real Doyle-Turner coefficients from "Fit/element.fit.dat"
        syscall = "grep 'a' Fit/"//trim(elu)//".fit.dat | awk '{print $2 "" "" $3 "" "" $4 "" "" $5 "" "" $6}' > "//trim(fitfile)        
        print *,trim(syscall)  
        call system( syscall )    
        syscall = "grep 'b' Fit/"//trim(elu)//".fit.dat | awk '{print $2 "" "" $3 "" "" $4 "" "" $5 "" "" $6}' >> "//trim(fitfile)        
        print *,trim(syscall)  
        call system( syscall )    
        
    !   and extract the imaginary Doyle-Turner coefficients from "Fit/element.VVVV.TTTT.tds.fit.dat"
        write(syscall,fmt='(g16.3)') V  
        write(dummy,fmt='(g16.3)') T  
        dummy = "Fit/"//trim(elu)//".V"//trim(adjustl(syscall))//".T"//trim(adjustl(dummy))//".tds.fit.dat"
        syscall = "grep 'a' "//trim(dummy)//" | awk '{print $2 "" "" $3 "" "" $4 "" "" $5 "" "" $6}' >> "//trim(fitfile)
        print *,trim(syscall)  
        call system( syscall )    
        syscall = "grep 'b' "//trim(dummy)//" | awk '{print $2 "" "" $3 "" "" $4 "" "" $5 "" "" $6}' >> "//trim(fitfile)  
        print *,trim(syscall)  
        call system( syscall )    
        
        
    !   finally find the Debye Waller factor
        write(dummy,fmt='(g16.3)') T        
        syscall = "./bin/DebyeWallerFactor.exe "//trim(elu)//" "//trim(adjustl(dummy))//" | grep 'B =' | awk '{print $3}' >> "//trim(fitfile)
        print *,trim(syscall)          
        call system( syscall )  
        
        
    !   finally read back in the complex Doyle-Turner coefficients and Debye Waller factor
               
        open(unit=701,file=trim(fitfile),action="read")
            read(unit=701,fmt=*) DTcoeffar
            read(unit=701,fmt=*) DTcoeffbr
            read(unit=701,fmt=*) DTcoeffai
            read(unit=701,fmt=*) DTcoeffbi
            read(unit=701,fmt=*) dwf
            do ii = 1,5
                aa(ii) = complex( DTcoeffar(ii),DTcoeffai(ii) )
                bb(ii) = complex( DTcoeffbr(ii),DTcoeffbi(ii) )
            end do                
        close(unit=701)     
        
    !---    and return the answer
        print *,"crystalStructureFactor info - DWF = ",dwf
        write (*,fmt='(4a16)') "Re(a)","Im(a)","Re(b)","Im(b)"
        do ii = 1,5
            write (*,fmt='(4f16.8)') real(aa(ii)),aimag(aa(ii)),real(bb(ii)),aimag(bb(ii)) 
        end do
        
        
        print *,""
        U0 = crystalStructureFactor( a0, lattice, (/0.0d0,0.0d0,0.0d0/), aa,bb,dwf )
        xi0 = PI * HBAR * velocity( V )/abs(U0)
        write (*,fmt='(5(a,f16.6))') "crystalStructureFactor U_0 = ",real(U0)," +i ",aimag(U0)," |U_0| = ",abs(U0)," eV , xi_0 ",xi0," A"
        Ug = crystalStructureFactor( a0, lattice, g, aa,bb,dwf )
        xig = PI * HBAR * velocity( V )/abs(Ug)
        write (*,fmt='(5(a,f16.6))') "crystalStructureFactor U_g = ",real(Ug)," +i ",aimag(Ug)," |U_g| = ",abs(Ug)," eV , xi_g ",xig," A"
        print *,""
        
        print *,"adding result to data table """//trim(outfile)//""""
        inquire(file=trim(outfile),exist=ok)
        if (.not. ok) then
            open(unit=800,file=trim(outfile),action="write")
                write(unit=800,fmt='(a)') "# crystalStructureFactors.exe v.1.0"
                write(unit=800,fmt='(a8,3a12,3a16)') "element","g ","T (K)","V (V)","Re(U_g)","Im(U_g)","xi_g"            
            close(unit=800)
        end if
        ok = .false. ; okg = .false.
        write(dummy0,fmt='(a8,3i4,2f12.3,3f16.6)') elu,0,0,0,T,V ,real(U0),aimag(U0),xi0 
        write(dummyg,fmt='(a8,3i4,2f12.3,3f16.6)') elu,nint(g),T,V ,real(Ug),aimag(Ug),xig 
        open(unit=710,file=trim(outfile),action="read")
            do 
                read(unit=710,fmt='(a44)',iostat=ii) dummy(1:44)
                if (ii/=0) exit
                if (dummy(1:44) == dummy0(1:44)) ok = .true.
                if (dummy(1:44) == dummyg(1:44)) okg = .true.
                if (ok .and. okg) exit                
            end do
        close(unit=710)
        if (.not. ok) print *,"adding result U_0"
        if (.not. okg) print *,"adding result U_g"
        if (.not. (ok .and. okg)) then
            open(unit=810,file=trim(outfile),action="write",position="append")
                if (.not. ok) write(unit=810,fmt='(a)') trim(dummy0)
                if (.not. okg) write(unit=810,fmt='(a)') trim(dummyg) 
            close(unit=810)
        else
            print *,"results for U_0 and U_g exist already"
        end if
        
        
                
        print *,""
        print *,"done"
        print *,""
        
    end program crystalStructureFactors         