
    program test_Lib_CrystalStructureFactor
!---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
        use Lib_ColouredTerminal
        use Lib_CommandLineArguments
        use Lib_Elements        
        use Lib_CrystalStructureFactor
        use Lib_RelativisticElectrons
        use Lib_InputTDS
        use Lib_InputReal
        use Lib_FitDoyleTurner
        use iso_fortran_env
        implicit none


        real(kind=real64),parameter     ::      PI =  3.14159265359d0      

        logical                 ::      ok
        character(len=8)        ::      element = "Fe"
        real(kind=real64)       ::      V = 200000
        real(kind=real64)       ::      T = 300
        integer,dimension(:,:),allocatable      ::      hkl
        integer,dimension(12000) ::      dat
        integer                 ::      nReflections
        character(len=8)        ::      latticename
        integer                 ::      ii,nS
        complex(kind=real64),dimension(:),allocatable   ::      Vg          !   (1:nReflections)
        complex(kind=real64),dimension(:),allocatable   ::      xig         !   (1:nReflections)
        complex(kind=real64),dimension(DT_COEFFS)       ::      aa,bb       !   fitted complex D-T coefficients
        type(CommandLineArguments)      ::      cla
        real(kind=real64)       ::      dwf,ds,ss,f_elastic,f_inelastic
        

        cla = CommandLineArguments_ctor(10)
        call setProgramDescription( cla, "test_Lib_CrystalStructureFactor" )
        call get( cla,"e",element  ,LIB_CLA_OPTIONAL,"        element" )  
        call get( cla,"V",V        ,LIB_CLA_OPTIONAL,"        accelerating voltage (V)" )  
        call get( cla,"T",T        ,LIB_CLA_OPTIONAL,"        temperature (K)" )  
        nReflections = 0
        call get( cla,"hkl",dat,nReflections,LIB_CLA_REQUIRED,"     reflections" )  

        latticename = getLatticeName(element) 
        if (latticename == "hcp") then
            if (mod(nReflections,4)/=0) stop "test_Lib_CrystalStructureFactor error - expect sets of hikl Miller-Bravais indices"
            nReflections = nReflections/4
            allocate(hkl(4,nReflections))
            hkl = reshape( dat(1:4*nReflections),(/4,nReflections/) )

        !---    check for [0000] - always needed
            ok = .false.
            do ii = 1,nReflections
                if (all(hkl(:,ii)==0)) ok = .true.
            end do
            if (.not. ok) then
                deallocate(hkl)
                allocate(hkl(4,nReflections+1))
                hkl(:,1) = 0
                hkl(:,2:nReflections+1) = reshape( dat(1:4*nReflections),(/4,nReflections/) )
                nReflections = nReflections+1
            end if
        else
            if (mod(nReflections,3)/=0) stop "test_Lib_CrystalStructureFactor error - expect sets of hkl Miller indices"
            nReflections = nReflections/3
            allocate(hkl(3,nReflections))
            hkl = reshape( dat(1:3*nReflections),(/3,nReflections/) )
            
        !---    check for [000] - always needed
            ok = .false.
            do ii = 1,nReflections
                if (all(hkl(:,ii)==0)) ok = .true.
            end do
            if (.not. ok) then
                deallocate(hkl)
                allocate(hkl(3,nReflections+1))
                hkl(:,1) = 0
                hkl(:,2:nReflections+1) = reshape( dat(1:3*nReflections),(/3,nReflections/) )
                nReflections = nReflections+1
            end if
        end if

        call report(cla)
        if (.not. allRequiredArgumentsSet(cla)) stop
        if (hasHelpArgument(cla)) stop
        call delete(cla)

        
    !---    compute crystal structure factors
        allocate(Vg(1:nReflections)) ; Vg = 0
        ! print *,""""//element//""""
        ! print *,hkl
        ! print *,T
        ! print *,V
        
        call computeCrystalStructureFactors( element,hkl(:,1:nReflections), T,V, Vg , aa,bb)
        allocate(xig(1:nReflections))
        do ii = 1,nReflections
            xig(ii) = PI * HBAR * velocity( V ) / Vg(ii)
        end do

    !---    output result
        if (latticename == "hcp") then
            write (*,fmt='(a8,4a16)') "hkil","Re(V_g)","Im(V_g)","Re(xi_g)","Im(xi_g)"
        else
             write (*,fmt='(a6,4a16)') "hkl","Re(V_g)","Im(V_g)","Re(xi_g)","Im(xi_g)"
        end if
        do ii = 1,nReflections
            if (latticename == "hcp") then
                write (*,fmt='(4i2,4f16.8)') hkl(:,ii) ,  real(Vg(ii)),aimag(Vg(ii)),  real(xig(ii)),aimag(xig(ii))
            else
                write (*,fmt='(3i2,4f16.8)') hkl(:,ii) ,  real(Vg(ii)),aimag(Vg(ii)),  real(xig(ii)),aimag(xig(ii))
            end if
        end do

        print *,""
        ns = 100
        ds = max(LIB_INPUTREAL_SMAX,LIB_INPUTTDS_SMAX)/ns
        dwf = DebyeWallerFactor(element,T,latticename) 
        write(*,fmt='(100a16)') "s=q/(4pi)","f_el(s)","f_inel(s)","f_DT^real(s)","f_DT^imag(s)"
        do ii = 0,ns
            ss = ii*ds
            f_elastic = scatteringAmplitude(element,ss)
            call getImagScatteringAmplitude(element,V,dwf,ss,f_inelastic)
            write(*,fmt='(100f16.8)') ss,f_elastic,f_inelastic,DoyleTurnerSum(ss,real(aa),real(bb)),DoyleTurnerSum(ss,aimag(aa),aimag(bb))
        end do
        print *,""


        ok = .true.
        
        if (ok) then
            print *,colour(GREEN,"PASS")
        else
            print *,colour(RED,"FAIL")
        end if


    end program test_Lib_CrystalStructureFactor