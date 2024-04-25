
    module EAM_PotentialTables
!---^^^^^^^^^^^^^^^^^^^^^^^^^^
!*      Read in embedded atom potential tables from disk, and construct an interpolating quintic spline table for them
!*          E_i = sum_j V(r_ij) + F[ sum_j phi(r_ij) ]
!*      the functions V,phi,F are tabulated

!*      the name of the potential should be of the form "A" or "A:B"


        use NBAX_StringTokenizers
        use Lib_NBAX
        use Lib_QuinticSplinesEven
        use iso_fortran_env
        implicit none
        private

    !---

        public      ::      EAM_PotentialTable_ctor
        public      ::      EAM_Alloy_ctor
        public      ::      report
        public      ::      delete


        public      ::      getE0
        public      ::      getCutoff
        public      ::      getPhi
        public      ::      getV
        public      ::      getF
        public      ::      getdPhidr
        public      ::      getdVdr
        public      ::      getdFdrho
        public      ::      getd2Phidr2
        public      ::      getd2Vdr2
        public      ::      getd2Fdrho2
        public      ::      getNTypes

        public      ::      getPhi_array_inrange_sum
        public      ::      getV_array_inrange_sum
        public      ::      getdPhidr_array_inrange
        public      ::      getdVdr_array_inrange

        public      ::      getEnergy
        public      ::      getElectronDensity
        public      ::      getEmbeddingEnergy
        public      ::      addForce
        public      ::      addHessian
        public      ::      HessianEinstein
        public      ::      hessianToDynamicalMatrix
        public      ::      addVirial

        public      ::      addTrS        
        public      ::      addElastic_constants        
        public      ::      addBulkModulus
       

        public      ::      outputAsXML
        public      ::      inputFromXML

        public      ::      setName
        public      ::      getName
        public      ::      getMass

    !---

        integer,public,parameter                ::      EAM_NAME_LEN = 16                
        real(kind=real64),public,parameter      ::      EAM_MASS_UNIT = 103.6427d0      !   1 Da in eV (fs/A)^2
    !---

        type,public     ::      EAM_PotentialTable
            private
            character(len=EAM_NAME_LEN) ::      name
            real(kind=real64)           ::      E0           !   energy scale ( equilibrium cohesive energy )
            real(kind=real64)           ::      mass         !
            real(kind=real64)           ::      r_phi_max
            real(kind=real64)           ::      r_V_max
            type(QuinticSplineEven)     ::      q_phi
            type(QuinticSplineEven)     ::      q_V
            type(QuinticSplineEven)     ::      q_F
        end type EAM_PotentialTable

        type,public     ::      EAM_Alloy
        !   
            private
            integer                                             ::      n               !   number of components in the alloy
             
            character(len=EAM_NAME_LEN),dimension(:),pointer ::      name
            real(kind=real64),dimension(:),pointer           ::      E0           !   energy scale ( equilibrium cohesive energy )
            real(kind=real64),dimension(:),pointer           ::      mass         !
            real(kind=real64),dimension(:),pointer           ::      r_phi_max
            real(kind=real64),dimension(:,:),pointer         ::      r_V_max        !   lower triangle filled
            type(QuinticSplineEven),dimension(:),pointer     ::      q_phi
            type(QuinticSplineEven),dimension(:,:),pointer   ::      q_V            !   lower triangle filled
            type(QuinticSplineEven),dimension(:),pointer     ::      q_F
             
            
        end type EAM_alloy

    !---


        interface   EAM_PotentialTable_ctor
            module procedure    EAM_PotentialTable_null
            module procedure    EAM_PotentialTable_ctor0
        end interface


        interface   EAM_Alloy_ctor
            module procedure    EAM_Alloy_null
            module procedure    EAM_Alloy_ctor0
        end interface

        interface   report
            module procedure    report0
            module procedure    report0a
            module procedure    report1
        end interface

        interface   delete
            module procedure    delete0
            module procedure    delete1
        end interface
        
        interface   getE0
            module procedure    getE00
            module procedure    getE01
        end interface
        interface   getCutoff
            module procedure    getCutoff0
            module procedure    getCutoff1
            module procedure    getCutoff2
        end interface

        interface   getPhi
            module procedure    getPhi0
            module procedure    getPhi1
        end interface
        
        interface   getV
            module procedure    getV0
            module procedure    getV1
        end interface
        
        interface   getF
            module procedure    getF0
            module procedure    getF1
        end interface


        interface   getdPhidr
            module procedure    getdPhidr0
            module procedure    getdPhidr1
        end interface
        interface   getdVdr
            module procedure    getdVdr0
            module procedure    getdVdr1
        end interface
        interface   getdFdrho
            module procedure    getdFdrho0
            module procedure    getdFdrho1
        end interface

        interface   getd2Phidr2
            module procedure    getd2Phidr20
            module procedure    getd2Phidr21
        end interface
        interface   getd2Vdr2
            module procedure    getd2Vdr20
            module procedure    getd2Vdr21
        end interface
        interface   getd2Fdrho2
            module procedure    getd2Fdrho20
            module procedure    getd2Fdrho21
        end interface



        interface   getPhi_array_inrange_sum
            module procedure    getPhi_array_inrange_sum0
        end interface
        interface   getV_array_inrange_sum
            module procedure    getV_array_inrange_sum0
        end interface


        interface   getdPhidr_array_inrange
            module procedure    getdPhidr_array_inrange0
        end interface
        interface   getdVdr_array_inrange
            module procedure    getdVdr_array_inrange0
        end interface
        interface   getdPhidr_array_inrange_sum
            module procedure    getdPhidr_array_inrange_sum0
        end interface

        interface   getEnergy
            module procedure    getEnergy0
            module procedure    getEnergy1
        end interface

        interface   getElectronDensity
            module procedure    getElectronDensity0
            module procedure    getElectronDensity1
        end interface
        
        interface   getEmbeddingEnergy
            module procedure    getEmbeddingEnergy0
            module procedure    getEmbeddingEnergy1
        end interface
        
        interface   addForce
            module procedure    addForce0
            module procedure    addForce0a
            module procedure    addForce0b
            module procedure    addForce1
            !module procedure    addForce1a
        end interface

        interface   addHessian
            module procedure    addHessian0
            module procedure    addHessian0a
            module procedure    addHessian1
            module procedure    addHessian1a
        end interface
        
        interface   HessianEinstein
            module procedure    HessianEinstein1a
        end interface
        
        interface   hessianToDynamicalMatrix 
            module procedure    hessianToDynamicalMatrix0
            module procedure    hessianToDynamicalMatrix1
        end interface
        

        interface   addVirial
            module procedure    addVirial0
            module procedure    addVirial1
        end interface



        interface   inputFromXML
            module procedure    inputFromXML0
            module procedure    inputFromXML1
        end interface

        interface   outputAsXML
            module procedure    outputAsXML0
!            module procedure    outputAsXML1
        end interface


        interface   getName
            module procedure    getName0
            module procedure    getName1
        end interface

        interface   setName
            module procedure    setName0
!             module procedure    setName1 
        end interface

        interface   getMass
            module procedure    getMass0
            module procedure    getMass1
        end interface

    !---

    contains
!---^^^^^^^^

        function EAM_PotentialTable_null() result(this)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
            type(EAM_PotentialTable)            ::      this

            this%E0 = 0.0d0
            this%r_phi_max = 0.0d0
            this%r_V_max = 0.0d0
            this%q_phi = QuinticSplineEven_ctor()
            this%q_V = QuinticSplineEven_ctor()
            this%q_F = QuinticSplineEven_ctor()
            this%name = ""
            return
        end function EAM_PotentialTable_null


        function EAM_PotentialTable_ctor0(E0,mass,q_phi,q_V,q_F,name) result(this)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
            real(kind=real64),intent(in)            ::      E0,mass       !   energy scale,mass
            type(QuinticSplineEven),intent(in)      ::      q_phi
            type(QuinticSplineEven),intent(in)      ::      q_V
            type(QuinticSplineEven),intent(in)      ::      q_F
            character(len=*),intent(in),optional    ::      name
            type(EAM_PotentialTable)            ::      this
            this%E0 = E0
            this%mass = mass
            this%q_phi = q_phi
            this%q_V = q_V
            this%q_F = q_F
            this%r_phi_max = getXmax(this%q_phi)
            this%r_V_max = getXmax(this%q_V)
            this%name = "" ; if (present(name)) this%name = name
            return
        end function EAM_PotentialTable_ctor0

    !---

        function EAM_Alloy_null() result(this)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
            type(EAM_Alloy)            ::      this
            this%n = 0
            
            nullify(this%name)
            nullify(this%E0)
            nullify(this%mass)
            nullify(this%r_phi_max)
            nullify(this%r_V_max)
            nullify(this%q_phi)
            nullify(this%q_v)
            nullify(this%q_f)
            
             
            return
        end function EAM_Alloy_null


        function EAM_Alloy_ctor0(E0,name,mass,q_phi,q_V,q_F) result(this)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
            real(kind=real64),dimension(:),intent(in)            ::      E0          !   energy scale,mass
            real(kind=real64),dimension(:),intent(in)            ::      mass        !   mass
            character(len=*),dimension(:),intent(in)             ::      name
            type(QuinticSplineEven),dimension(:),intent(in)      ::      q_phi
            type(QuinticSplineEven),dimension(:,:),intent(in)    ::      q_V
            type(QuinticSplineEven),dimension(:),intent(in)      ::      q_F
            
            type(EAM_Alloy)             ::      this
            integer                     ::      ii,jj
            this%n = size(E0)
            
            allocate(this%name(this%n))
            allocate(this%E0(this%n))
            allocate(this%mass(this%n))
            allocate(this%r_phi_max(this%n))
            allocate(this%r_V_max(this%n,this%n))
            allocate(this%q_phi(this%n))
            allocate(this%q_v(this%n,this%n))
            allocate(this%q_f(this%n))
            
            
            
            
            do ii = 1,this%n
                this%name(ii) = name(ii)
                this%E0(ii)   = E0(ii)
                this%mass(ii) = mass(ii)
                this%q_phi(ii) = QuinticSplineEven_ctor()
                call clone( this%q_phi(ii),q_phi(ii) )
                this%q_f(ii) = QuinticSplineEven_ctor()
                call clone( this%q_f(ii),q_f(ii) )
                this%r_phi_max(ii) = getXmax(this%q_phi(ii))
                do jj = 1,this%n
                    this%q_v(jj,ii) = QuinticSplineEven_ctor()
                    call clone( this%q_v(jj,ii),q_v(jj,ii) )
                     
                    this%r_V_max(jj,ii) = getXmax(this%q_V(jj,ii))
                end do
            end do                
            
            
            return
        end function EAM_Alloy_ctor0

    !---

        subroutine delete0(this)
    !---^^^^^^^^^^^^^^^^^^^^^^^^
            type(EAM_PotentialTable),intent(inout)          ::      this
            if (this%r_phi_max == 0.0d0) return
            call delete(this%q_phi)
            call delete(this%q_V)
            call delete(this%q_F)
            this = EAM_PotentialTable_null()
        end subroutine delete0

        subroutine delete1(this)
    !---^^^^^^^^^^^^^^^^^^^^^^^^
            type(EAM_Alloy),intent(inout)          ::      this
            integer         ::          ii,jj
            if (this%n == 0) return
            
            deallocate(this%name)
            deallocate(this%E0)
            deallocate(this%mass)
            deallocate(this%r_phi_max)
            deallocate(this%r_V_max)
             
            
            do ii = 1,this%n                
                call delete(this%q_phi(ii))
                call delete(this%q_f(ii))                
                do jj = 1,this%n
                    call delete(this%q_v(jj,ii))
                end do
            end do
            
            deallocate(this%q_phi)
            deallocate(this%q_f)                
            deallocate(this%q_v)                
            
            this = EAM_Alloy_null()
            return
        end subroutine delete1

    !---


        subroutine report0( this,u,o )
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
            type(EAM_PotentialTable),intent(in)     ::      this
            integer,intent(in),optional             ::      u,o
            integer             ::      uu,oo
            integer             ::      ii
            real(kind=real64),parameter ::  RMAX = 5.0d0,RHOMAX = 50.0d0
            integer,parameter   ::      NTABLE = 50
            real(kind=real64)   ::      rr,rho,V,dV,d2V,phi,dphi,d2phi,F,dF,d2F
            uu = 6 ; if (present(u)) uu = u
            oo = 0 ; if (present(o)) oo = o




            write (unit=uu,fmt='(a,3f12.5,a)') repeat(" ",oo)//"EAM_PotentialTable ("//trim(this%name)//")"
!             write (unit=uu,fmt='(a,3f12.5,a)') repeat(" ",oo+2)//"A0,E0,rc= ",this%A0,this%E0,getCutoff(this)
            write (unit=uu,fmt='(a,3f12.5,a)') repeat(" ",oo+2)//"E0,rc= ",this%E0,getCutoff(this)
            write (unit=uu,fmt='(a,3f12.5,a)') repeat(" ",oo+2)//"mass (Da,eV (fs/A)^2= ",getMass(this)/EAM_MASS_UNIT,getMass(this)
! 
!             write (unit=uu,fmt='(a,f16.8,a)') repeat(" ",oo+2)//"electronic density function ( r_phi_max = ",this%r_phi_max,")"
!             call report(this%q_phi,uu,oo+2)
!             write (unit=uu,fmt='(a,f16.8,a)') repeat(" ",oo+2)//"repulsive function ( r_V_max = ",this%r_V_max,")"
!             call report(this%q_V,uu,oo+2)
!             write (unit=uu,fmt='(a,i6)') repeat(" ",oo+2)//"embedding function"
!             call report(this%q_F,uu,oo+2)

            write(unit=uu,fmt='(a,100a16)') repeat(" ",oo+2),"rr","V","dV","d2V","phi","dphi","d2phi","rho","F","dF","d2F"
            do ii = 0,NTABLE
                rr = max(1.0d-8,ii*getCutoff(this)/NTABLE)
                call qsplint_inrange( this%q_V,rr,V,dV,d2V )
                call qsplint_inrange( this%q_phi,rr,phi,dphi,d2phi )
                rho = ii*getXmax(this%q_F)/NTABLE
                call qsplint( this%q_F,rho,F,dF,d2F )
                write(unit=uu,fmt='(a,100g16.8)') repeat(" ",oo+2),rr,V,dV,d2V,phi,dphi,d2phi,rho,F,dF,d2F
            end do
            
            
            
    !        print *,""
    !        do ii = 0,100
    !            rr = 1.0d0 + ii*1.0d-2
    !            call qsplint_inrange( this%q_V,rr,V,dV,d2V ) 
    !            write(unit=uu,fmt='(a,100f24.12)') repeat(" ",oo+2),rr,V,dV,d2V
    !        end do

            return
        end subroutine report0

        subroutine report0a( this,verbose,u,o )
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
            type(EAM_PotentialTable),intent(in)     ::      this
            logical,intent(in)                      ::      verbose
            integer,intent(in),optional             ::      u,o
            integer             ::      uu,oo
!            integer             ::      ii
            real(kind=real64),parameter ::  RMAX = 5.0d0,RHOMAX = 50.0d0
            integer,parameter   ::      NTABLE = 50
!            real(kind=real64)   ::      rr,rho  !   ,V,phi,F  !   ,dV,d2V,dphi,d2phi,dF,d2F
            uu = 6 ; if (present(u)) uu = u
            oo = 0 ; if (present(o)) oo = o
            
            call report0( this,uu,oo )
            if (.not. verbose) return
                        
            write (unit=uu,fmt='(a,f16.8,a)') repeat(" ",oo+2)//"electronic density function ( r_phi_max = ",this%r_phi_max,")"
            call report(this%q_phi,uu,oo+2)
            write (unit=uu,fmt='(a,f16.8,a)') repeat(" ",oo+2)//"repulsive function ( r_V_max = ",this%r_V_max,")"
            call report(this%q_V,uu,oo+2)
            write (unit=uu,fmt='(a,i6)') repeat(" ",oo+2)//"embedding function"
            call report(this%q_F,uu,oo+2)
            

            return
        end subroutine report0a


        subroutine report1( this,u,o )
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
            type(EAM_Alloy),intent(in)              ::      this
            integer,intent(in),optional             ::      u,o
            integer             ::      uu,oo
            integer             ::      jj
            uu = 6 ; if (present(u)) uu = u
            oo = 0 ; if (present(o)) oo = o
            write (unit=uu,fmt='(a)',advance="no") repeat(" ",oo)//"EAM_Alloy["
            do jj = 1,this%n-1                
                write (unit=uu,fmt='(a)',advance="no") trim(this%name(jj))//":"
            end do
            write (unit=uu,fmt='(a)',advance="yes") trim(this%name(this%n))//"]"             
            return
        end subroutine report1
    !---
!
! 

    !---

        pure function getNTypes(this) result(n)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
            type(EAM_Alloy),intent(in)              ::      this
            integer                                 ::      n
            n = this%n
            return
        end function getNTypes


        pure function getE00(this) result(E0)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
            type(EAM_PotentialTable),intent(in)     ::      this
            real(kind=real64)                       ::      E0
            E0 = this%E0
            return
        end function getE00

        pure function getCutoff0(this) result(rc)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
            type(EAM_PotentialTable),intent(in)     ::      this
            real(kind=real64)                       ::      rc
            rc = max(this%r_V_max,this%r_phi_max)
            return
        end function getCutoff0

    !---

        pure function getE01(this,ti) result(E0)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
            type(EAM_Alloy),intent(in)              ::      this
            integer,intent(in)                      ::      ti
            real(kind=real64)                       ::      E0
            E0 = this%E0(ti)
            return
        end function getE01

        pure function getCutoff1(this,ti) result(rc)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
            type(EAM_Alloy),intent(in)              ::      this
            integer,intent(in)                      ::      ti
            real(kind=real64)                       ::      rc
            integer         ::      tj
            rc = this%r_phi_max(ti)
            do tj = 1,ti-1
                rc = max(rc,this%r_V_max(tj,ti))
            end do
            return
        end function getCutoff1

        pure function getCutoff2(this) result(rc)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
            type(EAM_Alloy),intent(in)              ::      this

            real(kind=real64)                       ::      rc
            integer         ::      ti
            rc = 0.0d0
            do ti = 1,this%n
                rc = max(rc,getCutoff1(this,ti))
            end do
            return
        end function getCutoff2
    !---

        pure function getPhi0(this,r) result(phi)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
            type(EAM_PotentialTable),intent(in)     ::      this
            real(kind=real64),intent(in)            ::      r
            real(kind=real64)                       ::      phi
            phi = 0.00d0
            if (r<this%r_phi_max) phi = splint(this%q_phi,r)
            
            return
        end function getPhi0
        
        
        pure function getPhi1(this,ti,r) result(phi)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
            type(EAM_Alloy),intent(in)              ::      this
            integer,intent(in)                      ::      ti
            real(kind=real64),intent(in)            ::      r
            real(kind=real64)                       ::      phi
            phi = 0
            if (r<this%r_phi_max(ti)) phi = splint(this%q_phi(ti),r)
            return
        end function getPhi1

        pure function getV0(this,r) result(V)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
            type(EAM_PotentialTable),intent(in)     ::      this
            real(kind=real64),intent(in)            ::      r
            real(kind=real64)                       ::      V
            V = 0.00d0
            if (r<this%r_V_max)  V = splint(this%q_V,r)
            return
        end function getV0

        pure function getV1(this,ti,tj,r) result(V)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      note: lower triangle filled
            type(EAM_Alloy),intent(in)              ::      this
            integer,intent(in)                      ::      ti,tj
            real(kind=real64),intent(in)            ::      r
            real(kind=real64)                       ::      V
            integer         ::      tt1,tt2
             
            tt1 = max(ti,tj)
            tt2 = min(ti,tj)
            V = 0 
            if (r<this%r_V_max(tt1,tt2)) V = splint(this%q_V(tt1,tt2),r)
            return
        end function getV1        

        pure function getF0(this,rho) result(F)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
            type(EAM_PotentialTable),intent(in)     ::      this
            real(kind=real64),intent(in)            ::      rho
            real(kind=real64)                       ::      F
            F = splint(this%q_F,rho)
            return
        end function getF0



        pure function getF1(this,ti,rho) result(F)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
            type(EAM_Alloy),intent(in)              ::      this
            integer,intent(in)                      ::      ti
            real(kind=real64),intent(in)            ::      rho
            real(kind=real64)                       ::      F
            F = splint(this%q_F(ti),rho)
            return
        end function getF1

    !---

        pure subroutine getdPhidr0(this,r, phi,dphi)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
            type(EAM_PotentialTable),intent(in)     ::      this
            real(kind=real64),intent(in)            ::      r
            real(kind=real64),intent(out)           ::      phi,dphi
            if (r>=this%r_phi_max) then
                phi = 0.00d0
                dphi = 0.0d0
            else
                call qsplint(this%q_phi,r,phi,dphi)
            end if
            return
        end subroutine getdPhidr0

        pure subroutine getdVdr0(this,r, V,dV)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
            type(EAM_PotentialTable),intent(in)     ::      this
            real(kind=real64),intent(in)            ::      r
            real(kind=real64),intent(out)           ::      V,dV
            if (r>=this%r_V_max) then
                V = 0.00d0
                dV = 0.0d0
            else
                call qsplint(this%q_V,r,V,dV)
            end if
            return
        end subroutine getdVdr0

        pure subroutine getdFdrho0(this,rho, F,dF)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
            type(EAM_PotentialTable),intent(in)     ::      this
            real(kind=real64),intent(in)            ::      rho
            real(kind=real64),intent(out)           ::      F,dF
            call qsplint(this%q_F,rho,F,dF)
            return
        end subroutine getdFdrho0

    !---

        pure subroutine getdPhidr1(this,ti,r, phi,dphi)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
            type(EAM_Alloy),intent(in)              ::      this
            integer,intent(in)                      ::      ti
            real(kind=real64),intent(in)            ::      r
            real(kind=real64),intent(out)           ::      phi,dphi
            if (r>=this%r_phi_max(ti)) then
                phi = 0.00d0
                dphi = 0.0d0
            else
                call qsplint(this%q_phi(ti),r,phi,dphi)
            end if
            return
        end subroutine getdPhidr1

        pure subroutine getdVdr1(this,ti,tj,r, V,dV)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
            type(EAM_Alloy),intent(in)              ::      this
            integer,intent(in)                      ::      ti,tj
            real(kind=real64),intent(in)            ::      r
            real(kind=real64),intent(out)           ::      V,dV
            integer         ::      tt1,tt2
            tt1 = ti ; tt2 = tj
            if (tj>ti) then
                tt1 = ti ; tt2 = tj
            end if
            if (r>=this%r_V_max(tt1,tt2)) then
                V = 0.00d0
                dV = 0.0d0
            else
                call qsplint(this%q_V(tt1,tt2),r,V,dV)
            end if
            return
        end subroutine getdVdr1

        pure subroutine getdFdrho1(this,ti,rho, F,dF)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
            type(EAM_Alloy),intent(in)              ::      this
            integer,intent(in)                      ::      ti
            real(kind=real64),intent(in)            ::      rho
            real(kind=real64),intent(out)           ::      F,dF
            call qsplint(this%q_F(ti),rho,F,dF)
            return
        end subroutine getdFdrho1

    !---

        pure subroutine getd2Phidr20(this,r, phi,dphi,d2phi)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
            type(EAM_PotentialTable),intent(in)     ::      this
            real(kind=real64),intent(in)            ::      r
            real(kind=real64),intent(out)           ::      phi,dphi,d2phi
            if (r>=this%r_phi_max) then
                phi = 0.00d0
                dphi = 0.0d0
                d2phi = 0.0d0
            else
                call qsplint(this%q_phi,r,phi,dphi,d2phi)
            end if
            return
        end subroutine getd2Phidr20

        pure subroutine getd2Vdr20(this,r, V,dV,d2V)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
            type(EAM_PotentialTable),intent(in)     ::      this
            real(kind=real64),intent(in)            ::      r
            real(kind=real64),intent(out)           ::      V,dV,d2V
            if (r>=this%r_V_max) then
                V = 0.00d0
                dV = 0.0d0
                d2V = 0.0d0
            else
                call qsplint(this%q_V,r,V,dV,d2V)
            end if
            return
        end subroutine getd2Vdr20

        pure subroutine getd2Fdrho20(this,rho, F,dF,d2F)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
            type(EAM_PotentialTable),intent(in)     ::      this
            real(kind=real64),intent(in)            ::      rho
            real(kind=real64),intent(out)           ::      F,dF,d2F
            call qsplint(this%q_F,rho,F,dF,d2F)
            return
        end subroutine getd2Fdrho20

    !---

        pure subroutine getd2Phidr21(this,ti,r, phi,dphi,d2phi)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
            type(EAM_Alloy),intent(in)              ::      this
            integer,intent(in)                      ::      ti
            real(kind=real64),intent(in)            ::      r
            real(kind=real64),intent(out)           ::      phi,dphi,d2phi
            call qsplint(this%q_phi(ti),r,phi,dphi,d2phi)
            return
        end subroutine getd2Phidr21

        pure subroutine getd2Vdr21(this,ti,tj,r, V,dV,d2V)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
            type(EAM_Alloy),intent(in)              ::      this
            integer,intent(in)                      ::      ti,tj
            real(kind=real64),intent(in)            ::      r
            real(kind=real64),intent(out)           ::      V,dV,d2V
            integer         ::      tt1,tt2
            tt1 = ti ; tt2 = tj
            if (tj>ti) then
                tt1 = ti ; tt2 = tj
            end if
            if (r>=this%r_V_max(tt1,tt2)) then
                V = 0.00d0
                dV = 0.0d0
                d2V = 0.0d0
            else
                call qsplint(this%q_V(tt1,tt2),r,V,dV,d2V)
            end if
            return
        end subroutine getd2Vdr21

        pure subroutine getd2Fdrho21(this,ti,rho, F,dF,d2F)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
            type(EAM_Alloy),intent(in)              ::      this
            integer,intent(in)                      ::      ti
            real(kind=real64),intent(in)            ::      rho
            real(kind=real64),intent(out)           ::      F,dF,d2F
            call qsplint(this%q_F(ti),rho,F,dF,d2F)
            return
        end subroutine getd2Fdrho21

    !---

        pure function getPhi_array_inrange_sum0(this,n,r) result(rho)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      on exit returns the sum of phi
            type(EAM_PotentialTable),intent(in)                 ::      this
            integer,intent(in)                                  ::      n
            real(kind=real64),dimension(:),intent(in)           ::      r
            real(kind=real64)                                   ::      rho
            call qsplint_array_inrange_sum(this%q_phi,n,r,rho)
            return
        end function getPhi_array_inrange_sum0

        pure function getV_array_inrange_sum0(this,n,r) result(V)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
            type(EAM_PotentialTable),intent(in)                 ::      this
            integer,intent(in)                                  ::      n
            real(kind=real64),dimension(:),intent(in)           ::      r
            real(kind=real64)                                   ::      V
            call qsplint_array_inrange_sum(this%q_V,n,r,V)
            return
        end function getV_array_inrange_sum0


        pure subroutine getdPhidr_array_inrange0(this,n,r, phi,dphi)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      on entry r(:) is the length of 3 vector
    !*      on exit dphi(1:nn) = d phi/dx
            type(EAM_PotentialTable),intent(in)             ::      this
            integer,intent(in)                              ::      n
            real(kind=real64),dimension(:),intent(in)       ::      r
            real(kind=real64),dimension(:),intent(out)      ::      phi
            real(kind=real64),dimension(:),intent(out)      ::      dphi
            call qsplint_array_inrange(this%q_phi,n,r(:),phi(:),dphi(:))
            return
        end subroutine getdPhidr_array_inrange0


        pure subroutine getdPhidr_array_inrange_sum0(this,n,r, rho,dphi)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      on entry  r(:) being length of 3 vector
    !*      on exit dphi(1:3,1:nn) = d phi/dx and rho = sum phi
            type(EAM_PotentialTable),intent(in)             ::      this
            integer,intent(in)                              ::      n
            real(kind=real64),dimension(:),intent(in)       ::      r
            real(kind=real64),intent(out)                   ::      rho
            real(kind=real64),dimension(:),intent(out)      ::      dphi
            call qsplint_array_inrange_sum(this%q_phi,n,r(:),rho,dphi(:))
            return
        end subroutine getdPhidr_array_inrange_sum0


        pure subroutine getdVdr_array_inrange0(this,n,r, dV)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      on entry r(:) being length of 3 vector
            type(EAM_PotentialTable),intent(in)             ::      this
            integer,intent(in)                              ::      n
            real(kind=real64),dimension(:),intent(in)       ::      r
            real(kind=real64),dimension(:),intent(out)      ::      dV
            call qsplint_array_inrange(this%q_V,n,r(:),dV(:),derivative = .true.)
            return
        end subroutine getdVdr_array_inrange0

!-------

        pure function getEnergy0(this,n,r) result(E)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      note that this works for a single atom species only- otherwise the types need to be taken into account
            type(EAM_PotentialTable),intent(in)             ::      this
            integer,intent(in)                              ::      n
            real(kind=real64),dimension(:),intent(in)       ::      r
            real(kind=real64)                               ::      E
            real(kind=real64)           ::      V,rho
            E = 0.0d0 ; if (n==0) return
            V = getV_array_inrange_sum0(this,n,r(:))                !   compute V = sum( V(r) )
            rho = getPhi_array_inrange_sum0(this,n,r(:))            !   rho = sum( phi(r) )
            E = V/2 + getF(this,rho) - this%E0                      !   E = V/2 + F[rho]
            
            return
        end function getEnergy0

        pure function getEnergy1(this,n,t,r) result(E)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      t(0) is type of central atom, t(1:) are types of neighbours
            type(EAM_Alloy),intent(in)                      ::      this
            integer,intent(in)                              ::      n
            integer,dimension(0:),intent(in)                ::      t
            real(kind=real64),dimension(:),intent(in)       ::      r
            real(kind=real64)                               ::      E
            real(kind=real64)           ::      V,rho
            integer                     ::      ii
            E = 0.0d0 ; if (n==0) return
            V = 0.0d0 ; rho = 0.0d0
            do ii = 1,n
                V = V + getV( this,t(0),t(ii),r(ii) )
                rho = rho + getPhi( this,t(ii),r(ii) )
            end do
            E = V/2 + getF(this,t(0),rho) - this%E0(t(0))                      !   E = V/2 + F[rho]
            return
        end function getEnergy1


        pure function getElectronDensity0(this,n,r) result(rho)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      note that this works for a single atom species only- otherwise the types need to be taken into account
            type(EAM_PotentialTable),intent(in)             ::      this
            integer,intent(in)                              ::      n
            real(kind=real64),dimension(:),intent(in)       ::      r
            real(kind=real64)                               ::      rho
            rho = 0.0d0 ; if (n==0) return
            rho = getPhi_array_inrange_sum0(this,n,r(:))            !   rho = sum( phi(r) )
            return
        end function getElectronDensity0

        pure function getElectronDensity1(this,n,t,r) result(rho)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      t(0) is type of central atom, t(1:) are types of neighbours
            type(EAM_Alloy),intent(in)                      ::      this
            integer,intent(in)                              ::      n
            integer,dimension(0:),intent(in)                ::      t
            real(kind=real64),dimension(:),intent(in)       ::      r
            real(kind=real64)                               ::      rho
            integer                     ::      ii            
            rho = 0.0d0 ; if (n==0) return
            do ii = 1,n
                rho = rho + getPhi( this,t(ii),r(ii) )
            end do
            return
        end function getElectronDensity1


        pure function getEmbeddingEnergy0(this,n,r) result(E)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      note that this works for a single atom species only- otherwise the types need to be taken into account
            type(EAM_PotentialTable),intent(in)             ::      this
            integer,intent(in)                              ::      n
            real(kind=real64),dimension(:),intent(in)       ::      r
            real(kind=real64)                               ::      E
            real(kind=real64)           ::      rho
            E = 0.0d0 ; if (n==0) return
            rho = getPhi_array_inrange_sum0(this,n,r(:))            !   rho = sum( phi(r) )
            E = getF(this,rho)
            
            return
        end function getEmbeddingEnergy0

        real(kind=real64) function getEmbeddingEnergy1(this,n,t,r) !result(E)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      t(0) is type of central atom, t(1:) are types of neighbours
            type(EAM_Alloy),intent(in)                      ::      this
            integer,intent(in)                              ::      n
            integer,dimension(0:),intent(in)                ::      t
            real(kind=real64),dimension(:),intent(in)       ::      r
            !real(kind=real64)                               ::      E
            real(kind=real64)           ::      rho
            integer                     ::      ii
            getEmbeddingEnergy1 = 0.0d0 ; if (n==0) return
            rho = 0.0d0
            do ii = 1,n                
                rho = rho + getPhi( this,t(ii),r(ii) )
            end do
            !E = getF(this,t(0),rho) 
            getEmbeddingEnergy1 = getF(this,t(0),rho) 
            return
        end function getEmbeddingEnergy1


    !---

        subroutine addForce0(this,i,n,r,x,neigh,force)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ 
    !*      on entry x(1:3,1:n) with r(:) being length of 3 vector and x(1:3,:) of unit length
    !*      neigh(1:n) is the number of each neighbouring atom
    !*      add force to atom i
    !*      note that this works for a single atom species only- otherwise the types need to be taken into account
            type(EAM_PotentialTable),intent(in)             ::      this
            integer,intent(in)                              ::      i,n
            real(kind=real64),dimension(:),intent(in)       ::      r
            real(kind=real64),dimension(:,:),intent(in)     ::      x
            integer,dimension(:),intent(in)                 ::      neigh
            real(kind=real64),dimension(:,:),intent(inout)  ::      force
            real(kind=real64),dimension(n)                  ::      dV,dPhi       

            real(kind=real64)               ::      rho,dFdrho
            integer                         ::      kk,jj
            real(kind=real64)               ::      dEdr
            real(kind=real64),dimension(3)  ::      dEdx

            if (n==0) return
            call getdVdr_array_inrange0(this,n,r , dV)                !   compute dV/dr for each distance r
            call getdPhidr_array_inrange_sum0(this,n,r, rho,dPhi)     !   compute dphi/dr for each distance r
            call qsplint(this%q_F,rho,dFdrho)                                       !   compute dF/drho

            do kk = 1,n
                dEdr = dV(kk)/2 + dFdrho*dPhi(kk)                     ! = 1/2 dV/dr + dF/drho dphi/dr
                dEdx(1:3) =  dEdr * x(1:3,kk)                                      ! now = 1/2 dV/dr dr/dx + dF/drho dphi/dr dr/dx = 1/2 dV/dx + dF/dx
                force(1:3,i) = force(1:3,i) + dEdx(1:3)                             ! sign convention: x(:,kk) = x(:,jj)-x(:,ii)
                jj = neigh(kk)
                force(1:3,jj) = force(1:3,jj) - dEdx(1:3)
            end do

            return
        end subroutine addForce0

    !---
        subroutine addForce0a(this,i,n,r,x,neigh,force,workspace)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      on entry x(1:3,1:n) with r(:) being length of 3 vector and x(1:3,:) of unit length
    !*      neigh(1:n) is the number of each neighbouring atom
    !*      add force to atom i
    !*      note that this works for a single atom species only- otherwise the types need to be taken into account
            type(EAM_PotentialTable),intent(in)             ::      this
            integer,intent(in)                              ::      i,n
            real(kind=real64),dimension(:),intent(in)       ::      r
            real(kind=real64),dimension(:,:),intent(in)     ::      x
            integer,dimension(:),intent(in)                 ::      neigh
            real(kind=real64),dimension(:,:),intent(inout)  ::      force
            real(kind=real64),dimension(:),intent(inout)    ::      workspace       !   (1:2*n) allocating workspace memory externally saves time

            real(kind=real64)               ::      rho,dFdrho
            integer                         ::      kk,jj
            real(kind=real64)               ::      dEdr
            real(kind=real64),dimension(3)  ::      dEdx

            if (n==0) return
            call getdVdr_array_inrange0(this,n,r(:), workspace(1:n))                !   compute dV/dr for each distance r
            call getdPhidr_array_inrange_sum0(this,n,r(:), rho,workspace(n+1:))     !   compute dphi/dr for each distance r
            call qsplint(this%q_F,rho,dFdrho)                                       !   compute dF/drho

            do kk = 1,n
                dEdr = workspace(kk)/2 + dFdrho*workspace(n+kk)                     ! = 1/2 dV/dr + dF/drho dphi/dr
                dEdx(1:3) = dEdr * x(1:3,kk)                                      ! now = 1/2 dV/dr dr/dx + dF/drho dphi/dr dr/dx = 1/2 dV/dx + dF/dx
                force(1:3,i) = force(1:3,i) + dEdx(1:3)                             ! sign convention: x(:,kk) = x(:,jj)-x(:,ii)
                jj = neigh(kk)
                force(1:3,jj) = force(1:3,jj) - dEdx(1:3)
            end do

            return
        end subroutine addForce0a

        subroutine addForce0b(this,i,n,r,x,neigh,force)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ 
    !*      on entry x(1:3,1:n) with r(:) being length of 3 vector and x(1:3,:) of unit length
    !*      neigh(1:n) is the number of each neighbouring atom
    !*      add force to atom i
    !*      note that this works for a single atom species only- otherwise the types need to be taken into account
            type(EAM_PotentialTable),intent(in)             ::      this
            integer,intent(in)                              ::      i,n
            real(kind=real64),dimension(:),intent(in)       ::      r
            real(kind=real64),dimension(:,:),intent(in)     ::      x
            integer,dimension(:),intent(in)                 ::      neigh
            real(kind=real64),dimension(:),intent(inout)    ::      force
            real(kind=real64),dimension(n)                  ::      dV,dPhi       

            real(kind=real64)               ::      rho,dFdrho
            integer                         ::      kk,jj
            real(kind=real64)               ::      dEdr
            real(kind=real64),dimension(3)  ::      dEdx

            if (n==0) return
            call getdVdr_array_inrange0(this,n,r , dV)                !   compute dV/dr for each distance r
            call getdPhidr_array_inrange_sum0(this,n,r, rho,dPhi)     !   compute dphi/dr for each distance r
            call qsplint(this%q_F,rho,dFdrho)                         !   compute dF/drho

            do kk = 1,n
                dEdr = dV(kk)/2 + dFdrho*dPhi(kk)                     ! = 1/2 dV/dr + dF/drho dphi/dr
                dEdx(1:3) = dEdr * x(1:3,kk)                          ! now = 1/2 dV/dr dr/dx + dF/drho dphi/dr dr/dx = 1/2 dV/dx + dF/dx
                force(i*3-2:i*3) = force(i*3-2:i*3) + dEdx(1:3)       ! sign convention: x(:,kk) = x(:,jj)-x(:,ii)
                jj = neigh(kk)
                force(jj*3-2:jj*3) = force(jj*3-2:jj*3) - dEdx(1:3)
            end do

            return
        end subroutine addForce0b

    !---
    !---

        subroutine addForce1(this,i,n,t,r,x,neigh,force)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ 
    !*      on entry x(1:3,1:n) with r(:) being length of 3 vector and x(1:3,:) of unit length
    !*      t(0) is type of central atom, t(1:) are types of neighbours
    !*      neigh(1:n) is the number of each neighbouring atom
    !*      add force to atom i
            type(EAM_Alloy),intent(in)                      ::      this
            integer,intent(in)                              ::      i,n
            integer,dimension(0:),intent(in)                ::      t
            real(kind=real64),dimension(:),intent(in)       ::      r
            real(kind=real64),dimension(:,:),intent(in)     ::      x
            integer,dimension(:),intent(in)                 ::      neigh
            real(kind=real64),dimension(:,:),intent(inout)  ::      force
            real(kind=real64),dimension(2*n)                ::      workspace 
            real(kind=real64)               ::      vv,phi,rho,ff,dFdrho
            integer                         ::      kk,jj
            real(kind=real64)               ::      dEdr
            real(kind=real64),dimension(3)  ::      dEdx

            if (n==0) return

            rho = 0.0d0
            do kk = 1,n
                call getdVdr( this,t(0),t(kk),r(kk),vv,workspace(kk) )
                call getdphidr( this,t(kk),r(kk),phi,workspace(n+kk) )
                rho = rho + phi
            end do
            call getdFdrho(this,t(0),rho,ff,dfdrho)

            do kk = 1,n
                dEdr = workspace(kk)/2 + dFdrho*workspace(n+kk)                     ! = 1/2 dV/dr + dF/drho dphi/dr
                dEdx(1:3) = - dEdr * x(1:3,kk)                                      ! now = 1/2 dV/dr dr/dx + dF/drho dphi/dr dr/dx = 1/2 dV/dx + dF/dx
                force(1:3,i) = force(1:3,i) - dEdx(1:3)                             ! sign convention: x(:,kk) = x(:,jj)-x(:,ii)
                jj = neigh(kk)
                force(1:3,jj) = force(1:3,jj) + dEdx(1:3)
            end do

            return
        end subroutine addForce1

!         subroutine addForce1a(this,i,n,t,r,x,neigh,force,workspace)
!     !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
!     !*      on entry x(1:3,1:n) with r(:) being length of 3 vector and x(1:3,:) of unit length
!     !*      t(0) is type of central atom, t(1:) are types of neighbours
!     !*      neigh(1:n) is the number of each neighbouring atom
!     !*      add force to atom i
!     !*      note that this works for a single atom species only- otherwise the types need to be taken into account
!             type(EAM_Alloy),intent(in)                      ::      this
!             integer,intent(in)                              ::      i,n
!             integer,dimension(0:),intent(in)                ::      t
!             real(kind=real64),dimension(:),intent(in)       ::      r
!             real(kind=real64),dimension(:,:),intent(in)     ::      x
!             integer,dimension(:),intent(in)                 ::      neigh
!             real(kind=real64),dimension(:,:),intent(inout)  ::      force
!             real(kind=real64),dimension(:),intent(inout)    ::      workspace       !   (1:2*n) allocating workspace memory externally saves time
! 
!             real(kind=real64)               ::      vv,phi,rho,ff,dFdrho
!             integer                         ::      kk,jj
!             real(kind=real64)               ::      dEdr
!             real(kind=real64),dimension(3)  ::      dEdx
! 
!             if (n==0) return
! 
!             rho = 0.0d0
!             do kk = 1,n
!                 call getdVdr( this,t(0),t(kk),r(kk),vv,workspace(kk) )
!                 call getdphidr( this,t(0),t(kk),r(kk),phi,workspace(n+kk) )
!                 rho = rho + phi
!             end do
!             call getdFdrho(this,t(0),rho,ff,dfdrho)
! 
!             do kk = 1,n
!                 dEdr = workspace(kk)/2 + dFdrho*workspace(n+kk)                     ! = 1/2 dV/dr + dF/drho dphi/dr
!                 dEdx(1:3) = - dEdr * x(1:3,kk)                                      ! now = 1/2 dV/dr dr/dx + dF/drho dphi/dr dr/dx = 1/2 dV/dx + dF/dx
!                 force(1:3,i) = force(1:3,i) - dEdx(1:3)                             ! sign convention: x(:,kk) = x(:,jj)-x(:,ii)
!                 jj = neigh(kk)
!                 force(1:3,jj) = force(1:3,jj) + dEdx(1:3)
!             end do
! 
!             return
!         end subroutine addForce1a
! 

!         subroutine addForce2(this,i,n,t,r,x,neigh,force)
!     !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
!     !*      on entry x(1:3,1:n) with r(:) being length of 3 vector and x(1:3,:) of unit length
!     !*      t(0) is type of central atom, t(1:) are types of neighbours
!     !*      neigh(1:n) is the number of each neighbouring atom
!     !*      add force to atom i
!     !*      note that this works for a single atom species only- otherwise the types need to be taken into account
!             type(EAM_Alloy),intent(in)                      ::      this
!             integer,intent(in)                              ::      i,n
!             integer,dimension(0:),intent(in)                ::      t
!             real(kind=real64),dimension(:),intent(in)       ::      r
!             real(kind=real64),dimension(:,:),intent(in)     ::      x
!             integer,dimension(:),intent(in)                 ::      neigh
!             real(kind=real64),dimension(:,:),intent(inout)  ::      force
!             
!             real(kind=real64),dimension(2*n)    ::      workspace       
! 
!             real(kind=real64)               ::      vv,phi,rho,ff,dFdrho
!             integer                         ::      kk,jj
!             real(kind=real64)               ::      dEdr
!             real(kind=real64),dimension(3)  ::      dEdx
! 
!             if (n==0) return
! 
!             rho = 0.0d0
!             do kk = 1,n
!                 call getdVdr( this,t(0),t(kk),r(kk),vv,workspace(kk) )
!                 call getdphidr( this,t(0),t(kk),r(kk),phi,workspace(n+kk) )
!                 rho = rho + phi
!             end do
!             call getdFdrho(this,t(0),rho,ff,dfdrho)
! 
!             do kk = 1,n
!                 dEdr = workspace(kk)/2 + dFdrho*workspace(n+kk)                     ! = 1/2 dV/dr + dF/drho dphi/dr
!                 dEdx(1:3) = - dEdr * x(1:3,kk)                                      ! now = 1/2 dV/dr dr/dx + dF/drho dphi/dr dr/dx = 1/2 dV/dx + dF/dx
!                 force(1:3,i) = force(1:3,i) - dEdx(1:3)                             ! sign convention: x(:,kk) = x(:,jj)-x(:,ii)
!                 jj = neigh(kk)
!                 force(1:3,jj) = force(1:3,jj) + dEdx(1:3)
!             end do
! 
!             return
!         end subroutine addForce2
! 
!     !---
!     
    

        subroutine addVirial0(this,n,t,r,x,virial)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      on entry x(1:3,1:n) with r(:) being length of 3 vector and x(1:3,:) of unit length
    !*      t(0) is type of central atom, t(1:) are types of neighbours
    !*      neigh(1:n) is the number of each neighbouring atom
    !*      add virial stress for atom 
            type(EAM_Alloy),intent(in)                      ::      this
            integer,intent(in)                              ::      n
            integer,dimension(0:),intent(in)                ::      t
            real(kind=real64),dimension(:),intent(in)       ::      r
            real(kind=real64),dimension(:,:),intent(in)     ::      x
            real(kind=real64),dimension(3,3),intent(inout)  ::      virial
            
            real(kind=real64),dimension(2*n)    ::      workspace       

            real(kind=real64)               ::      vv,phi,rho,ff,dFdrho
            integer                         ::      kk
            real(kind=real64)               ::      dEdr
            if (n==0) return

            rho = 0.0d0
            do kk = 1,n
                call getdVdr( this,t(0),t(kk),r(kk),vv,workspace(kk) )
                call getdphidr( this,t(kk),r(kk),phi,workspace(n+kk) )
                rho = rho + phi
            end do
            call getdFdrho(this,t(0),rho,ff,dfdrho)

            do kk = 1,n
                dEdr = workspace(kk)/2 + dFdrho*workspace(n+kk)                     ! = 1/2 dV/dr + dF/drho dphi/dr
                virial(1:3,1:3) = virial(1:3,1:3) - dEdr * drdstrain( x(1:3,kk),r(kk) )
            end do

            return
        end subroutine addVirial0
        
        subroutine addVirial1(this,i,n,t,r,x,neigh,virial)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      on entry x(1:3,1:n) with r(:) being length of 3 vector and x(1:3,:) of unit length
    !*      t(0) is type of central atom, t(1:) are types of neighbours
    !*      neigh(1:n) is the number of each neighbouring atom
    !*      add virial stress to atoms
            type(EAM_Alloy),intent(in)                      ::      this
            integer,intent(in)                              ::      i,n
            integer,dimension(0:),intent(in)                ::      t
            real(kind=real64),dimension(:),intent(in)       ::      r
            real(kind=real64),dimension(:,:),intent(in)     ::      x
            integer,dimension(:),intent(in)                 ::      neigh
            real(kind=real64),dimension(:,:,:),intent(inout)::      virial
            
            real(kind=real64),dimension(2*n)    ::      workspace       

            real(kind=real64)                   ::      vv,phi,rho,ff,dFdrho
            integer                             ::      kk,jj
            real(kind=real64)                   ::      dEdr
            real(kind=real64),dimension(3,3)    ::      stress
            
            if (n==0) return

            rho = 0.0d0
            do kk = 1,n
                call getdVdr( this,t(0),t(kk),r(kk),vv,workspace(kk) )
                call getdphidr( this,t(kk),r(kk),phi,workspace(n+kk) )
                rho = rho + phi
            end do
            call getdFdrho(this,t(0),rho,ff,dfdrho)

            do kk = 1,n
                dEdr = workspace(kk)/2 + dFdrho*workspace(n+kk)                     ! = 1/2 dV/dr + dF/drho dphi/dr
                stress(1:3,1:3) = - dEdr * drdstrain( x(1:3,kk),r(kk) )/2
                virial(1:3,1:3,i) = virial(1:3,1:3,i) + stress(1:3,1:3)
                jj = neigh(kk)
                virial(1:3,1:3,jj) = virial(1:3,1:3,jj) + stress(1:3,1:3)
            end do

            return
        end subroutine addVirial1
        

!         subroutine add3x3Hess(this,i,n,t,r,x,neigh,hess)
!     !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
!     !*      on entry x(1:3,1:n) with r(:) being length of 3 vector and x(1:3,:) of unit length
!     !*      t(0) is type of central atom, t(1:) are types of neighbours
!     !*      neigh(1:n) is the number of each neighbouring atom
!     !*      add 3x3 block of Hessian diagonal elements for atom
!             type(EAM_Alloy),intent(in)                      ::      this
!             integer,intent(in)                              ::      i,n
!             integer,dimension(0:),intent(in)                ::      t
!             real(kind=real64),dimension(:),intent(in)       ::      r
!             real(kind=real64),dimension(:,:),intent(in)     ::      x
!             integer,dimension(:),intent(in)                 ::      neigh
!             real(kind=real32),dimension(:,:,:),intent(inout)::      hess
!             
!             real(kind=real64)                   ::  imodrij,rho
!             integer                             ::  jj,kk
!             real(kind=real64)                   ::  VV,dVdr,d2Vdr2
!             real(kind=real64)                   ::  phi,dphidr,d2phidr2
!             real(kind=real64)                   ::  FF,dFdrho,d2Fdrho2
!             real(kind=real64),dimension(3)      ::  dr1
!             real(kind=real64),dimension(3,3)    ::  hs,dr2
!             real(kind=real64),dimension(3,0:n)  ::  fi
! 
! !
!             if (n==0) return
!             rho = 0.0d0
!             do kk = 1,n
!                 rho = rho + getPhi( this,t(kk),r(kk) )
!             end do
!             call getd2Fdrho2(this,t(0),rho,FF,dFdrho,d2Fdrho2)
! 
!             fi = 0.0d0
!             do kk = 1,n
! 
!                 call getd2Vdr2( this,t(0),t(kk),r(kk),VV,dVdr,d2Vdr2)
!                 call getd2Phidr2( this,t(kk),r(kk),phi,dphidr,d2phidr2)
! !
!                 imodrij = 1.0d0/r(kk)
!                 dr1(1:3) = -x(1:3,kk)
! 
! 
!                 fi(1:3,0) = fi(1:3,0) + dphidr*dr1(1:3)   ! = sum_j dphi_{ij}/dx_i = drho_i/dx_i
!                 fi(1:3,kk) =  - dphidr*dr1(1:3)         ! = dphi_{ij}/dx_j = drho_i/dx_j
! 
!                 dr2(1:3,1:3) = d2rijbydxidxi( x(1:3,kk),imodrij )
! 
!                 hs(1:3,1) = ( (d2Vdr2/2+dFdrho*d2phidr2)*dr1(1:3)*dr1(1) + (dVdr/2+dFdrho*dphidr)*dr2(1:3,1) )
!                 hs(1:3,2) = ( (d2Vdr2/2+dFdrho*d2phidr2)*dr1(1:3)*dr1(2) + (dVdr/2+dFdrho*dphidr)*dr2(1:3,2) )
!                 hs(1:3,3) = ( (d2Vdr2/2+dFdrho*d2phidr2)*dr1(1:3)*dr1(3) + (dVdr/2+dFdrho*dphidr)*dr2(1:3,3) )
!                 jj = neigh(kk)
! 
!                 hess(1:3,1:3,i) = hess(1:3,1:3,i)   + real(hs(1:3,1:3),kind=real32)
!                 hess(1:3,1:3,jj) = hess(1:3,1:3,jj) + real(hs(1:3,1:3),kind=real32)               
! 
!             end do
! !             stop
! 
!             hs(1:3,1) = d2Fdrho2*fi(1:3,0)*fi(1,0)          !   = d2E/drho_i2 * drho_i/dx_i * drho_i/dx_i
!             hs(1:3,2) = d2Fdrho2*fi(1:3,0)*fi(2,0)
!             hs(1:3,3) = d2Fdrho2*fi(1:3,0)*fi(3,0)
!             hess(1:3,1:3,i) = hess(1:3,1:3,i) + real(hs(1:3,1:3),kind=real32)   
! 
!             do kk = 1,n
!                 jj = neigh(kk)
! !
! !                hs(1:3,1) = d2Fdrho2*fi(1:3,0)*fi(1,kk)        !    = d2E/drho_i2 * drho_i/dx_i * drho_i/dx_j
! !                hs(1:3,2) = d2Fdrho2*fi(1:3,0)*fi(2,kk)
! !                hs(1:3,3) = d2Fdrho2*fi(1:3,0)*fi(3,kk)
! !                call add( hess,i,jj,hs )
! !                call add( hess,jj,i,transpose(hs) )
! !
! !                do ll = 1,n
! !                    mm = neigh(ll)
! !                    if (mm/=jj) cycle
! !                    hs(1:3,1) = d2Fdrho2*fi(1:3,ll)*fi(1,kk)      !    = d2E/drho_i2 * drho_i/dx_m * drho_i/dx_j
! !                    hs(1:3,2) = d2Fdrho2*fi(1:3,ll)*fi(2,kk)
! !                    hs(1:3,3) = d2Fdrho2*fi(1:3,ll)*fi(3,kk)
!                     hs(1:3,1) = d2Fdrho2*fi(1:3,kk)*fi(1,kk)      !    = d2E/drho_i2 * drho_i/dx_m * drho_i/dx_j
!                     hs(1:3,2) = d2Fdrho2*fi(1:3,kk)*fi(2,kk)
!                     hs(1:3,3) = d2Fdrho2*fi(1:3,kk)*fi(3,kk)
! 
! !                    call add( hess,mm,jj,hs )                 !     note (mm,jj) and (jj,mm) are both visited in this loop
!                     hess(1:3,1:3,jj) = hess(1:3,1:3,jj) + real(hs(1:3,1:3),kind=real32)   
! 
! !                end do
!             end do
! 
!             return
!         end subroutine add3x3Hess
        
        subroutine addTrS(this,i,n,t,r,x,neigh,trS)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      on entry x(1:3,1:n) with r(:) being length of 3 vector and x(1:3,:) of unit length
    !*      t(0) is type of central atom, t(1:) are types of neighbours
    !*      neigh(1:n) is the number of each neighbouring atom
    !*      add trace of stress for atom i
            type(EAM_Alloy),intent(in)                      ::      this
            integer,intent(in)                              ::      i,n
            integer,dimension(0:),intent(in)                ::      t
            real(kind=real64),dimension(:),intent(in)       ::      r
            real(kind=real64),dimension(:,:),intent(in)     ::      x
            integer,dimension(:),intent(in)                 ::      neigh
            real(kind=real64),dimension(:),intent(inout)    ::      trS 
            
            real(kind=real64),dimension(2*n)    ::      workspace       

            real(kind=real64)                   ::      vv,phi,rho,ff,dFdrho
            integer                             ::      kk
            real(kind=real64)                   ::      dEdr,TrSon2
            real(kind=real64),dimension(3,3)    ::      stress
            if (n==0) return

            rho = 0.0d0
            do kk = 1,n
                call getdVdr( this,t(0),t(kk),r(kk),vv,workspace(kk) )
                call getdphidr( this,t(kk),r(kk),phi,workspace(n+kk) )
                rho = rho + phi
            end do
            call getdFdrho(this,t(0),rho,ff,dfdrho)

            do kk = 1,n
                dEdr = workspace(kk)/2 + dFdrho*workspace(n+kk)                     ! = 1/2 dV/dr + dF/drho dphi/dr
                stress(1:3,1:3) = - dEdr * drdstrain( x(1:3,kk),r(kk) )
                TrSon2 = ( stress(1,1) + stress(2,2) + stress(3,3) )/2
                trS(i) = trS(i) + TrSon2
                trS( neigh(kk) )  = trS( neigh(kk) ) + TrSon2
            end do

            return
        end subroutine addTrS
        
        
        
        
        

        subroutine addElastic_constants(this,n,t,r,x,c)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      on entry x(1:3,1:n) with r(:) being length of 3 vector and x(1:3,:) of unit length
    !*      t(0) is type of central atom, t(1:) are types of neighbours
    !*      add elastic constants tensor for atom 
    !*      note that this works for a single atom species only- otherwise the types need to be taken into account
            type(EAM_Alloy),intent(in)                      ::      this
            integer,intent(in)                              ::      n
            integer,dimension(0:),intent(in)                ::      t
            real(kind=real64),dimension(:),intent(in)       ::      r
            real(kind=real64),dimension(:,:),intent(in)     ::      x
            real(kind=real64),dimension(3,3,3,3),intent(inout)  ::      c
            
            real(kind=real64)               ::      VV,dVdr,d2Vdr2
            real(kind=real64)               ::      phi,dphidr,d2phidr2
            real(kind=real64)               ::      FF,dFdrho,d2Fdrho2
            real(kind=real64)               ::      dedr,d2edr2
            integer                         ::      kk
            real(kind=real64)               ::      rho
            real(kind=real64),dimension(3,3)    ::      gamma,drde
            integer                         ::      aa,bb
            if (n==0) return

            rho = 0.0d0
            do kk = 1,n
                rho = rho + getPhi( this,t(kk),r(kk) )
            end do
            call getd2Fdrho2(this,t(0),rho,FF,dFdrho,d2Fdrho2)

            gamma = 0.0d0
            do kk = 1,n

                call getd2Vdr2( this,t(0),t(kk),r(kk),VV,dVdr,d2Vdr2)
                call getd2Phidr2( this,t(kk),r(kk),phi,dphidr,d2phidr2)
            
                dedr = (dVdr/2 + dFdrho*dphidr)
                d2edr2 = (d2Vdr2/2 + dFdrho*d2phidr2)
                
                c(1:3,1:3,1:3,1:3) = c(1:3,1:3,1:3,1:3) + dedr * d2rdstrain2( x(1:3,kk),r(kk) )
                drde(1:3,1:3) = drdstrain( x(1:3,kk),r(kk) )
                gamma(1:3,1:3) = gamma(1:3,1:3) + dphidr*drde(1:3,1:3)
                
                do bb = 1,3
                    do aa = 1,3
                        c(aa,bb,1:3,1:3) = c(aa,bb,1:3,1:3) + d2edr2*drde(1:3,1:3)*drde(aa,bb)
                    end do
                end do
                
            end do
            
            do bb = 1,3
                do aa = 1,3
                    c(aa,bb,1:3,1:3) = c(aa,bb,1:3,1:3) + d2Fdrho2*gamma(1:3,1:3)*gamma(aa,bb)
                end do
            end do
            
            return
        end subroutine addElastic_constants
        
        pure function drdstrain( x,r ) result(dr)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      x_ij is the unit vector, magnitude r
    !*      return derivative of r_ij wrt strain
            real(kind=real64),dimension(3),intent(in)       ::      x
            real(kind=real64),intent(in)                    ::      r
            real(kind=real64),dimension(3,3)                ::      dr
            
            dr(1:3,1) = x(1:3)*x(1)*r               
            dr(1:3,2) = x(1:3)*x(2)*r
            dr(1:3,3) = x(1:3)*x(3)*r
            return
        end function drdstrain
                                   
        function d2rdstrain2( x,r ) result(d2r)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      x_ij is the unit vector, magnitude r
    !*      return 2nd derivative of r_ij wrt strain
            real(kind=real64),dimension(3),intent(in)       ::      x
            real(kind=real64),intent(in)                    ::      r
            real(kind=real64),dimension(3,3,3,3)            ::      d2r
            integer                     ::      aa,bb,cc,dd
                                                                
            do dd = 1,3
                do cc = 1,3
                    do bb = 1,3
                        do aa = 1,3
                            d2r(aa,bb,cc,dd) = - x(aa)*x(bb)*x(cc)*x(dd)
                            if (aa==cc) d2r(aa,bb,cc,dd) = d2r(aa,bb,cc,dd) + x(bb)*x(dd)
                            if (bb==cc) d2r(aa,bb,cc,dd) = d2r(aa,bb,cc,dd) + x(aa)*x(dd)
                            if (aa==dd) d2r(aa,bb,cc,dd) = d2r(aa,bb,cc,dd) + x(bb)*x(cc)
                            if (bb==dd) d2r(aa,bb,cc,dd) = d2r(aa,bb,cc,dd) + x(aa)*x(cc)
                        end do
                    end do
                end do
            end do
            d2r(1:3,1:3,1:3,1:3) = d2r(1:3,1:3,1:3,1:3)*r
            return
        end function d2rdstrain2
        
    !---        

        subroutine addBulkModulus(this,n,t,r,b)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      on entry r(:) being length of 3 vector 
    !*      t(0) is type of central atom, t(1:) are types of neighbours
    !*      add bulk modulus for an atom 
    !*      note: returns B*V - divide by volume for correct bulk mod.
            type(EAM_Alloy),intent(in)                      ::      this
            integer,intent(in)                              ::      n
            integer,dimension(0:),intent(in)                ::      t
            real(kind=real64),dimension(:),intent(in)       ::      r
            real(kind=real64),intent(inout)                 ::      b
            
            real(kind=real64)               ::      VV,dVdr,d2Vdr2
            real(kind=real64)               ::      phi,dphidr,d2phidr2
            real(kind=real64)               ::      FF,dFdrho,d2Fdrho2
            integer                         ::      kk
            real(kind=real64)               ::      rho,aa,bb,cc
            
            if (n==0) return

            rho = 0.0d0
            do kk = 1,n
                rho = rho + getPhi( this,t(kk),r(kk) )
            end do
            call getd2Fdrho2(this,t(0),rho,FF,dFdrho,d2Fdrho2)

            aa = 0.0d0 ; bb = 0.0d0 ; cc = 0.0d0
            do kk = 1,n

                call getd2Vdr2( this,t(0),t(kk),r(kk),VV,dVdr,d2Vdr2)
                call getd2Phidr2( this,t(kk),r(kk),phi,dphidr,d2phidr2)
                
                aa = aa + r(kk)*( r(kk)*d2Vdr2 - 2*dVdr )
                bb = bb + dFdrho*r(kk)*( r(kk)*d2phidr2 - 2*dphidr )
                cc = cc + r(kk)*dphidr
                            
            end do
                        
            b = b + (aa/2 + bb + d2Fdrho2*cc*cc)/9
            
            
            return
        end subroutine addBulkModulus
                
        
     
    !---

        subroutine hessianToDynamicalMatrix0(this,hess)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      divide element i,j by sqrt(mi mj)
            type(EAM_PotentialTable), intent(in)            ::  this
            real(kind=real32),dimension(:,:),intent(inout)                      ::  hess
            
             
            real(kind=real64)       ::      mi 
            real(kind=real32)       ::      isqmimj              
            
            mi = getMass(this)
            isqmimj = real(1/mi,kind=real32)
            hess = hess * isqmimj
             
            return
        end subroutine hessianToDynamicalMatrix0
            
            

        subroutine hessianToDynamicalMatrix1(this,t,hess)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      divide element i,j by sqrt(mi mj)
            type(EAM_Alloy), intent(in)            ::      this
            integer,dimension(:),intent(in)                 ::      t
            real(kind=real32),dimension(:,:),intent(inout)                      ::  hess
            
             
            integer                 ::      ii,jj,ti,tj
            real(kind=real64),dimension(this%n)       ::      isqmi
            real(kind=real64)       ::      mi
            real(kind=real32)       ::      isqmimj 
            integer         ::      nn
            
            do ii = 1,this%n
                mi = getMass(this,ii)
                isqmi(ii) = real(1/sqrt(mi),kind=real32)
            end do
             
            nn = size(t)
             
            do jj = 1,nn
                tj = t(jj)
                do ii = 1,nn
                    ti = t(ii)
                    isqmimj = real(isqmi(ti)*isqmi(tj),kind=real32)
                    hess(3*ii-2:3*ii,3*jj-2:3*jj) = hess(3*ii-2:3*ii,3*jj-2:3*jj) * isqmimj
                end do
            end do
            return
        end subroutine hessianToDynamicalMatrix1
            
            

        subroutine addHessian0(this, i,n,r,x,neigh, hess )
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      Calculates the Hessian ( force constant matrix) due a single atom i given its neighbour list
    !*      on input x are unit vectors with length r
    !*      note that this works for a single atom species only- otherwise the types need to be taken into account
    !*      note this is the dense version with a full 3Nx3N matrix
    !*      System potential energy is 
    !*             E = sum_i  1/2 sum_j V(r_ij) + F[ sum_j phi(r_ij) ] 
    
            type(EAM_PotentialTable), intent(in)            ::  this
            integer,intent(in)                              ::  i,n            !   atom i is being computed. n neighbours
            real(kind=real64),dimension(:),intent(in)       ::  r              !   (1:n) distance of atom i to each neighbour
            real(kind=real64),dimension(:,:),intent(in)     ::  x              !   (1:3,1:n) unit vector to each neighbour
            integer,dimension(:),intent(in)                 ::  neigh          !   (1:n) index of each neighbour
            real(kind=real32),dimension(:,:),intent(inout)  ::  hess           !   (1:3N,1:3N) dense matrix representing Hessian
 
            
            real(kind=real64)                   ::  imodrij,rho
            integer                             ::  jj,kk,ll,mm
            real(kind=real64)                   ::  VV,dVdr,d2Vdr2
            real(kind=real64)                   ::  phi,dphidr,d2phidr2
            real(kind=real64)                   ::  FF,dFdrho,d2Fdrho2
            real(kind=real64),dimension(3)      ::  dr1
            real(kind=real64),dimension(3,3)    ::  hs,dr2
            real(kind=real64),dimension(3,0:n)  ::  fi
 
        !---   quick return if no atoms in the neighbour list
            if (n==0) return
            
                                    
        !---   compute rho = sum_i phi( r_i )    
            rho = getPhi_array_inrange_sum0(this,n,r(:))             
            
        !---   compute F[rho], F'[rho], F"[rho]
            call qsplint(this%q_F,rho,FF,dFdrho,d2Fdrho2)
            
           ! print *,"DBG r       ",r(1:n)
           ! print *,"DBG rho     ",rho
           ! print *,"DBG F,F',F''",FF,dFdrho,d2Fdrho2
            
                         
 
            fi = 0.0d0         ! = sum_j dphi_{ij}/dx_i = drho_i/dx_i
            do kk = 1,n
            
            !---   compute V(r_k) and V' and V"
                call qsplint_inrange(this%q_V,r(kk),VV,dVdr,d2Vdr2)
            
            !---   compute phi(r_k) and phi' and phi"
                call qsplint_inrange(this%q_phi,r(kk),phi,dphidr,d2phidr2)
 
                !print *,"DBG k,r(k),V,phi ",kk,r(kk),VV,dVdr,d2Vdr2,phi,dphidr,d2phidr2
                
                
                imodrij = 1.0d0/r(kk)                  
                dr1(1:3) = -x(1:3,kk)
 
                fi(1:3,0) = fi(1:3,0) + dphidr*dr1(1:3)         ! = sum_j dphi_{ij}/dx_i = drho_i/dx_i
                fi(1:3,kk) =  - dphidr*dr1(1:3)                 ! = dphi_{ij}/dx_j = drho_i/dx_j
 
                dr2(1:3,1:3) = d2rijbydxidxi( x(1:3,kk),imodrij )
 
!                hs(1:3,1) = ( (d2Vdr2+ 2*dFdrho*d2phidr2)*dr1(1:3)*dr1(1) + (dVdr+ 2*dFdrho*dphidr)*dr2(1:3,1) )
!                hs(1:3,2) = ( (d2Vdr2+ 2*dFdrho*d2phidr2)*dr1(1:3)*dr1(2) + (dVdr+ 2*dFdrho*dphidr)*dr2(1:3,2) )
!                hs(1:3,3) = ( (d2Vdr2+ 2*dFdrho*d2phidr2)*dr1(1:3)*dr1(3) + (dVdr+ 2*dFdrho*dphidr)*dr2(1:3,3) )
                hs(1:3,1) = ( (d2Vdr2/2+dFdrho*d2phidr2)*dr1(1:3)*dr1(1) + (dVdr/2+dFdrho*dphidr)*dr2(1:3,1) )
                hs(1:3,2) = ( (d2Vdr2/2+dFdrho*d2phidr2)*dr1(1:3)*dr1(2) + (dVdr/2+dFdrho*dphidr)*dr2(1:3,2) )
                hs(1:3,3) = ( (d2Vdr2/2+dFdrho*d2phidr2)*dr1(1:3)*dr1(3) + (dVdr/2+dFdrho*dphidr)*dr2(1:3,3) )
                jj = neigh(kk)
 
                hess(3*i -2:3*i ,3*i -2:3*i)   = hess(3*i -2:3*i ,3*i -2:3*i) + real(hs(1:3,1:3),kind=real32)
                hess(3*jj-2:3*jj,3*jj-2:3*jj)  = hess(3*jj-2:3*jj,3*jj-2:3*jj) + real(hs(1:3,1:3),kind=real32)
                hess(3*i -2:3*i ,3*jj-2:3*jj)  = hess(3*i -2:3*i ,3*jj-2:3*jj) - real(hs(1:3,1:3),kind=real32)
                hess(3*jj-2:3*jj,3*i -2:3*i )  = hess(3*jj-2:3*jj,3*i -2:3*i ) - real(hs(1:3,1:3),kind=real32)
                 
                !call add( hess,i,i,hs )
                !call add( hess,jj,jj,hs )
                !call add( hess,i,jj,-hs )
                !call add( hess,jj,i,-hs )       
            end do
 
            !stop
            
             
             !   print *,"addHessian0 ",i,rho,FF,dFdrho,d2Fdrho2,fi(:,0)
             !   stop
         
            
            hs(1:3,1) =  d2Fdrho2*fi(1:3,0)*fi(1,0)          !   = d2E/drho_i2 * drho_i/dx_i * drho_i/dx_i
            hs(1:3,2) =  d2Fdrho2*fi(1:3,0)*fi(2,0)
            hs(1:3,3) =  d2Fdrho2*fi(1:3,0)*fi(3,0)
            
            hess(3*i -2:3*i ,3*i -2:3*i )  = hess(3*i -2:3*i ,3*i -2:3*i )   + real(hs(1:3,1:3),kind=real32)
  
                 !call add( hess,i,i,hs )
                 ! 
 
            do kk = 1,n
                jj = neigh(kk)
 
                hs(1:3,1) = d2Fdrho2*fi(1:3,0)*fi(1,kk)        !    = d2E/drho_i2 * drho_i/dx_i * drho_i/dx_j
                hs(1:3,2) = d2Fdrho2*fi(1:3,0)*fi(2,kk)
                hs(1:3,3) = d2Fdrho2*fi(1:3,0)*fi(3,kk)
                
                hess(3*i -2:3*i ,3*jj-2:3*jj)  = hess(3*i -2:3*i ,3*jj-2:3*jj) + real(hs(1:3,1:3),kind=real32)
                hess(3*jj-2:3*jj,3*i -2:3*i )  = hess(3*jj-2:3*jj,3*i -2:3*i ) + transpose( real(hs(1:3,1:3),kind=real32) )
                    
                 !call add( hess,i,jj,hs )
                 !call add( hess,jj,i,transpose(hs) )
                do ll = 1,n
                    mm = neigh(ll)
                    hs(1:3,1) = d2Fdrho2*fi(1:3,ll)*fi(1,kk)      !    = d2E/drho_i2 * drho_i/dx_m * drho_i/dx_j
                    hs(1:3,2) = d2Fdrho2*fi(1:3,ll)*fi(2,kk)
                    hs(1:3,3) = d2Fdrho2*fi(1:3,ll)*fi(3,kk)
 
                    hess(3*mm-2:3*mm,3*jj-2:3*jj)  = hess(3*mm-2:3*mm,3*jj-2:3*jj) + real(hs(1:3,1:3),kind=real32)                      
  !    call add( hess,mm,jj,hs )                 !     note (mm,jj) and (jj,mm) are both visited in this loop
 
                end do
            end do
 
            
            
            return
        end subroutine addHessian0
 
     !---
         

            subroutine addHessian0a(this, i,n,r,x,neigh, hess,indx )
        !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
        !*      Calculates the Hessian ( force constant matrix) due a single atom i given its neighbour list
        !*      on input x are unit vectors with length r
        !*      note that this works for a single atom species only- otherwise the types need to be taken into account
        !*      note this is the dense version with a full 3Nx3N matrix
        !*      System potential energy is 
        !*             E = sum_i  1/2 sum_j V(r_ij) + F[ sum_j phi(r_ij) ] 
        
                type(EAM_PotentialTable), intent(in)            ::  this
                integer,intent(in)                              ::  i,n            !   atom i is being computed. n neighbours
                real(kind=real64),dimension(:),intent(in)       ::  r              !   (1:n) distance of atom i to each neighbour
                real(kind=real64),dimension(:,:),intent(in)     ::  x              !   (1:3,1:n) unit vector to each neighbour
                integer,dimension(:),intent(in)                 ::  neigh          !   (1:n) index of each neighbour
                real(kind=real32),dimension(:,:,:,:),pointer,intent(inout)      ::      hess        !   (1:3*bandwidth, 1:3*nAtoms)
                integer,dimension(:,:),pointer,intent(inout)                    ::      indx        !   (0:bandwidth,1:nAtoms)           with 0 = number of neighbours, 1: = neighbour index
     
                
                real(kind=real64)                   ::  imodrij,rho
                integer                             ::  jj,kk,ll,mm
                real(kind=real64)                   ::  VV,dVdr,d2Vdr2
                real(kind=real64)                   ::  phi,dphidr,d2phidr2
                real(kind=real64)                   ::  FF,dFdrho,d2Fdrho2
                real(kind=real64),dimension(3)      ::  dr1
                real(kind=real64),dimension(3,3)    ::  dr2
                real(kind=real32),dimension(3,3)    ::  hs
                real(kind=real64),dimension(3,0:n)  ::  fi
     
                real(kind=real32),dimension(3,3,0:n)  ::  hsi
                
                
            !---   quick return if no atoms in the neighbour list
                if (n==0) return
                
                                        
            !---   compute rho = sum_i phi( r_i )    
                rho = getPhi_array_inrange_sum0(this,n,r(:))             
                
            !---   compute F[rho], F'[rho], F"[rho]
                call qsplint(this%q_F,rho,FF,dFdrho,d2Fdrho2)
                 
                             
                hsi = 0.0
     
                fi = 0.0d0         ! = sum_j dphi_{ij}/dx_i = drho_i/dx_i
                do kk = 1,n
                
                !---   compute V(r_k) and V' and V"
                    call qsplint_inrange(this%q_V,r(kk),VV,dVdr,d2Vdr2)
                
                !---   compute phi(r_k) and phi' and phi"
                    call qsplint_inrange(this%q_phi,r(kk),phi,dphidr,d2phidr2)
     
                    !print *,"DBG k,r(k),V,phi ",kk,r(kk),VV,dVdr,d2Vdr2,phi,dphidr,d2phidr2
                    
                    
                    imodrij = 1.0d0/r(kk)                  
                    dr1(1:3) = -x(1:3,kk)
     
                    fi(1:3,0) = fi(1:3,0) + dphidr*dr1(1:3)         ! = sum_j dphi_{ij}/dx_i = drho_i/dx_i
                    fi(1:3,kk) =  - dphidr*dr1(1:3)                 ! = dphi_{ij}/dx_j = drho_i/dx_j
     
                    dr2(1:3,1:3) = d2rijbydxidxi( x(1:3,kk),imodrij )
      
                    hs(1:3,1) = real( (d2Vdr2/2+dFdrho*d2phidr2)*dr1(1:3)*dr1(1) + (dVdr/2+dFdrho*dphidr)*dr2(1:3,1) , kind=real32 )
                    hs(1:3,2) = real( (d2Vdr2/2+dFdrho*d2phidr2)*dr1(1:3)*dr1(2) + (dVdr/2+dFdrho*dphidr)*dr2(1:3,2) , kind=real32 )
                    hs(1:3,3) = real( (d2Vdr2/2+dFdrho*d2phidr2)*dr1(1:3)*dr1(3) + (dVdr/2+dFdrho*dphidr)*dr2(1:3,3) , kind=real32 )
                    jj = neigh(kk)
      
                     
                 
                    call addSparseHessian3x3( hess,indx,jj,jj,hs )
                   
                    call addSparseHessian3x3( hess,indx,jj,i,-hs ) 
                    
                    hsi(1:3,1:3,0) = hsi(1:3,1:3,0)   + hs(1:3,1:3)
                    hsi(1:3,1:3,kk) = hsi(1:3,1:3,kk) - hs(1:3,1:3)
                          
                end do
      
             
                
                hs(1:3,1) = real( d2Fdrho2*fi(1:3,0)*fi(1,0), kind=real32 )           !   = d2E/drho_i2 * drho_i/dx_i * drho_i/dx_i
                hs(1:3,2) = real( d2Fdrho2*fi(1:3,0)*fi(2,0), kind=real32 ) 
                hs(1:3,3) = real( d2Fdrho2*fi(1:3,0)*fi(3,0), kind=real32 ) 
                 
                hsi(1:3,1:3,0) = hsi(1:3,1:3,0)   + hs(1:3,1:3) 
     
                do kk = 1,n
                    jj = neigh(kk)
     
                    hs(1:3,1) = real( d2Fdrho2*fi(1:3,0)*fi(1,kk), kind=real32 )         !    = d2E/drho_i2 * drho_i/dx_i * drho_i/dx_j
                    hs(1:3,2) = real( d2Fdrho2*fi(1:3,0)*fi(2,kk), kind=real32 ) 
                    hs(1:3,3) = real( d2Fdrho2*fi(1:3,0)*fi(3,kk), kind=real32 ) 
                     
                    hsi(1:3,1:3,kk) = hsi(1:3,1:3,kk)   + hs(1:3,1:3) 
                    call addSparseHessian3x3( hess,indx,jj,i,transpose(hs) )
                    do ll = 1,n
                        mm = neigh(ll)
                        hs(1:3,1) = real( d2Fdrho2*fi(1:3,ll)*fi(1,kk) , kind=real32 )     !    = d2E/drho_i2 * drho_i/dx_m * drho_i/dx_j
                        hs(1:3,2) = real( d2Fdrho2*fi(1:3,ll)*fi(2,kk) , kind=real32 )
                        hs(1:3,3) = real( d2Fdrho2*fi(1:3,ll)*fi(3,kk) , kind=real32 )
          
                        call addSparseHessian3x3( hess,indx, mm,jj,hs )                 !     note (mm,jj) and (jj,mm) are both visited in this loop
     
                    end do
                end do
     
                call addSparseHessian3x3xm( hess,indx,i,n,neigh,hsi )                    
                
                
                return
            end subroutine addHessian0a
         
            subroutine reallocateHessian( H,indx,n )
        !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
        !*      reallocate the hessian
                    
                real(kind=real32),dimension(:,:,:,:),pointer      ::      H           !   (1:3*bandwidth, 1:3*nAtoms)
                integer,dimension(:,:),pointer                    ::      indx        !   (0:bandwidth,1:nAtoms)           with 0 = number of neighbours, 1: = neighbour index
                integer,intent(in),optional                       ::      n   
                
                integer             ::      nAtoms
                
                real(kind=real32),dimension(:,:,:,:),pointer      ::      H_tmp
                integer,dimension(:,:),pointer                    ::      indx_tmp
                
                integer             ::      LL,LL_tmp      !   bandwidth
                
                LL = size(indx,dim=1)-1
                nAtoms = size(indx,dim=2)
                if (present(n)) then
                    LL_tmp = max( int(LL*1.25),LL+n )
                else
                    LL_tmp = max( int(LL*1.25),LL+10 )
                end if
                print *,"reallocateHessian ",LL,"->",LL_tmp
                allocate(H_tmp(3,3,LL_tmp,nAtoms))
                allocate(indx_tmp(0:LL_tmp,1:nAtoms))
                H_tmp(1:3,1:3,1:LL,1:nAtoms) = H(1:3,1:3,1:LL,1:nAtoms)
                indx_tmp(0:LL,1:nAtoms)  = indx(0:LL,1:nAtoms)                       
                deallocate(H)
                deallocate(indx)
                H => H_tmp
                indx => indx_tmp
                return
            end subroutine reallocateHessian                
                    
                    
                 

            subroutine addSparseHessian3x3( hess,indx,i,j,h3x3 )
        !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ 
        !*      add matrix 3x3 block to atom position i,j
                real(kind=real32),dimension(:,:,:,:),pointer    ::      hess        !   1:3*bandwidth,1:3*nAtoms
                integer,dimension(:,:),pointer                  ::      indx        !   note: 0:bandwidth,1:nAtoms
                integer,intent(in)                              ::      i,j
                real(kind=real32),dimension(3,3),intent(in)     ::      h3x3

                integer             ::      kk,nn
                
                
                nn = indx(0,i)      !   number of neighbours for atom i : note indx(1:)!
                
                do kk = 1,nn
                    if (indx(kk,i) == j) then
                        !   add to an existing entry
                        hess(1:3,1:3,kk,i) = hess(1:3,1:3,kk,i) + h3x3(1:3,1:3)
                        return
                    end if
                end do
               
            !---    new neighbour
                if (nn+1 >= size(indx,dim=1)) then              !   because indx(0:L) and so if n+1>L then already out of bounds    
                !   requires a reallocation
                    call reallocateHessian( hess,indx )
                end if

                
                nn = nn + 1
                indx(0,i) = nn
                indx(nn,i) = j
                hess(1:3,1:3,nn,i) = h3x3(1:3,1:3)

                                
                return
            end subroutine addSparseHessian3x3
                
           subroutine addSparseHessian3x3xm( hess,indx,i,m,neigh,h3x3xm )
        !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ 
        !*      add matrix 3x3 block to atom position i 
                real(kind=real32),dimension(:,:,:,:),pointer    ::      hess        !   1:3*bandwidth,1:3*nAtoms
                integer,dimension(:,:),pointer                  ::      indx        !   note: 0:bandwidth,1:nAtoms
                integer,intent(in)                              ::      i,m
                integer,dimension(:),intent(in)                 ::      neigh
                real(kind=real32),dimension(:,:,0:),intent(in)  ::      h3x3xm

                integer                 ::      kk,jj,nn,mm
                logical,dimension(0:m)  ::      done
                
!                if (lbound(indx,dim=1)==1) stop "addSparseHessian3x3 error - passed pointer with lbound = 1, expected lbound = 0"
                
                
                
                nn = indx(0,i)      !   number of neighbours for atom i : note indx(1:)!
                                
            !---    find which members of the neighbour list are already in the matrix, and add those entries
                done = .false.
                do kk = 1,nn
                    jj = indx(kk,i)
                    if (jj == i) then
                        !   add self to an existing entry
                        hess(1:3,1:3,kk,i) = hess(1:3,1:3,kk,i) + h3x3xm(1:3,1:3,0)
                        done(0) = .true.
                    else
                        do mm = 1,m
                            if (jj == neigh(mm)) then
                                hess(1:3,1:3,kk,i) = hess(1:3,1:3,kk,i) + h3x3xm(1:3,1:3,mm)
                                done(mm) = .true.
                                exit
                            end if  
                        end do                                                         
                    end if
                end do
                
            !---    any more entries to add?
                kk = (m+1) - count( done )            !   number of entries needed to add to matrix
                if (nn + kk >= size(indx,dim=1)) call reallocateHessian( hess,indx,kk )
                
                if (.not. done(0)) then                        
                    nn = nn + 1 
                    indx(0,i) = nn
                    indx(nn,i) = i
                    hess(1:3,1:3,nn,i) = h3x3xm(1:3,1:3,0)
                end if
                
                do mm = 1,m
                    if (done(mm)) cycle                                                           
                        
                    nn = nn + 1
                    indx(0,i) = nn
                    indx(nn,i) = neigh(mm)
                    hess(1:3,1:3,nn,i) = h3x3xm(1:3,1:3,mm)
    
                end do
                                    
                return
            end subroutine addSparseHessian3x3xm
                
                
     !---
 
         subroutine addHessian1(this, i,n,t,r,x,neigh, hess )
     !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
     !*      Calculates the Hessian ( force constant matrix) due a single atom i given its neighbour list
     !*      on input x are unit vectors with length r
     !*      t(0) is type of central atom, t(1:) are types of neighbours
     !*      note that this works for a single atom species only- otherwise the types need to be taken into account
 
             type(EAM_Alloy), intent(in)                     ::  this
             integer,intent(in)                              ::  i,n
             integer,dimension(0:),intent(in)                ::      t
             real(kind=real64),dimension(:),intent(in)       ::  r
             real(kind=real64),dimension(:,:),intent(in)     ::  x
             integer,dimension(:),intent(in)                 ::  neigh
             !type(Matrix),intent(inout)                      ::  hess
             real(kind=real32),dimension(:,:),intent(inout)                      ::  hess
             
 !             integer                             ::  nn
             real(kind=real64)                   ::  imodrij,rho
             integer                             ::  jj,kk,ll,mm
             real(kind=real64)                   ::  VV,dVdr,d2Vdr2
             real(kind=real64)                   ::  phi,dphidr,d2phidr2
             real(kind=real64)                   ::  FF,dFdrho,d2Fdrho2
             real(kind=real64),dimension(3)      ::  dr1
             real(kind=real64),dimension(3,3)    ::  hs,dr2
             real(kind=real64),dimension(3,0:n)  ::  fi
 
 
 
 !
             if (n==0) return
             rho = 0.0d0
             do kk = 1,n
                 rho = rho + getPhi( this,t(kk),r(kk) )
             end do
             call getd2Fdrho2(this,t(0),rho,FF,dFdrho,d2Fdrho2)
 
             fi = 0.0d0
             do kk = 1,n
 
                 call getd2Vdr2( this,t(0),t(kk),r(kk),VV,dVdr,d2Vdr2)
                 call getd2Phidr2( this,t(kk),r(kk),phi,dphidr,d2phidr2)
 !
                 imodrij = 1.0d0/r(kk)
                 dr1(1:3) = -x(1:3,kk)
 
 
                 fi(1:3,0) = fi(1:3,0) + dphidr*dr1(1:3)   ! = sum_j dphi_{ij}/dx_i = drho_i/dx_i
                 fi(1:3,kk) =  - dphidr*dr1(1:3)         ! = dphi_{ij}/dx_j = drho_i/dx_j
 
                 dr2(1:3,1:3) = d2rijbydxidxi( x(1:3,kk),imodrij )
 
                 hs(1:3,1) = ( (d2Vdr2/2+dFdrho*d2phidr2)*dr1(1:3)*dr1(1) + (dVdr/2+dFdrho*dphidr)*dr2(1:3,1) )
                 hs(1:3,2) = ( (d2Vdr2/2+dFdrho*d2phidr2)*dr1(1:3)*dr1(2) + (dVdr/2+dFdrho*dphidr)*dr2(1:3,2) )
                 hs(1:3,3) = ( (d2Vdr2/2+dFdrho*d2phidr2)*dr1(1:3)*dr1(3) + (dVdr/2+dFdrho*dphidr)*dr2(1:3,3) )
                 jj = neigh(kk)
 
 !                print *,"addHess(A) ",i,jj,modrij
 
 
                 hess(3*i -2:3*i ,3*i -2:3*i )  = hess(3*i -2:3*i ,3*i -2:3*i ) + real(hs(1:3,1:3),kind=real32)
                 hess(3*jj-2:3*jj,3*jj-2:3*jj)  = hess(3*jj-2:3*jj,3*jj-2:3*jj) + real(hs(1:3,1:3),kind=real32)
                 hess(3*i -2:3*i ,3*jj-2:3*jj)  = hess(3*i -2:3*i ,3*jj-2:3*jj) - real(hs(1:3,1:3),kind=real32)
                 hess(3*jj-2:3*jj,3*i -2:3*i )  = hess(3*jj-2:3*jj,3*i -2:3*i ) - real(hs(1:3,1:3),kind=real32)
                 
 
                !call add( hess,i,i,hs )
                !call add( hess,jj,jj,hs )
                !call add( hess,i,jj,-hs )
                !call add( hess,jj,i,-hs )               !   note hs = transpose(hs)
 
 !                 if (i==1) print *,r(kk),dFdrho,dphidr,dr2(:,1),dphidr*dr2(:,1),dFdrho*dphidr*dr2(:,1)
 !                 if (i==1) write (*,fmt='(a,100f16.8)') "        ",r(kk),dVdr,d2Vdr2,dFdrho,d2Fdrho2,dphidr,d2phidr2,dr1,dr2,hs
 
             end do
 !             stop
 
             hs(1:3,1) = d2Fdrho2*fi(1:3,0)*fi(1,0)          !   = d2E/drho_i2 * drho_i/dx_i * drho_i/dx_i
             hs(1:3,2) = d2Fdrho2*fi(1:3,0)*fi(2,0)
             hs(1:3,3) = d2Fdrho2*fi(1:3,0)*fi(3,0)
             
             hess(3*i -2:3*i ,3*i -2:3*i )  = hess(3*i -2:3*i ,3*i -2:3*i )   + real(hs(1:3,1:3),kind=real32)
 
            !call add( hess,i,i,hs )
 
             do kk = 1,n
                 jj = neigh(kk)
 
                 hs(1:3,1) = d2Fdrho2*fi(1:3,0)*fi(1,kk)        !    = d2E/drho_i2 * drho_i/dx_i * drho_i/dx_j
                 hs(1:3,2) = d2Fdrho2*fi(1:3,0)*fi(2,kk)
                 hs(1:3,3) = d2Fdrho2*fi(1:3,0)*fi(3,kk)
                 
                 hess(3*i -2:3*i ,3*jj-2:3*jj)  = hess(3*i -2:3*i ,3*jj-2:3*jj) + real(hs(1:3,1:3),kind=real32)
                 hess(3*jj-2:3*jj,3*i -2:3*i )  = hess(3*jj-2:3*jj,3*i -2:3*i ) + transpose( real(hs(1:3,1:3),kind=real32) )
                 
                 
                 !call add( hess,i,jj,hs )
                 !call add( hess,jj,i,transpose(hs) )
 
                 do ll = 1,n
                     mm = neigh(ll)
                     hs(1:3,1) = d2Fdrho2*fi(1:3,ll)*fi(1,kk)      !    = d2E/drho_i2 * drho_i/dx_m * drho_i/dx_j
                     hs(1:3,2) = d2Fdrho2*fi(1:3,ll)*fi(2,kk)
                     hs(1:3,3) = d2Fdrho2*fi(1:3,ll)*fi(3,kk)
 
 !                     dr1 = rij(:,kk) - rij(:,ll)
 !                     print *,"addHess(B) ",jj,mm,sqrt(dot_product(dr1,dr1))
 
                    hess(3*mm-2:3*mm,3*jj-2:3*jj)  = hess(3*mm-2:3*mm,3*jj-2:3*jj) + real(hs(1:3,1:3),kind=real32)
                     
                 !    call add( hess,mm,jj,hs )                 !     note (mm,jj) and (jj,mm) are both visited in this loop
 
 
                 end do
             end do
 
             return
         end subroutine addHessian1

         


            subroutine addHessian1a(this, i,n,t,r,x,neigh, hess,indx )
        !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
        !*      Calculates the Hessian ( force constant matrix) due a single atom i given its neighbour list
        !*      on input x are unit vectors with length r
        !*      note that this works for a single atom species only- otherwise the types need to be taken into account
        !*      note this is the dense version with a full 3Nx3N matrix
        !*      System potential energy is 
        !*             E = sum_i  1/2 sum_j V(r_ij) + F[ sum_j phi(r_ij) ] 
        
         
                type(EAM_Alloy), intent(in)                     ::  this
                integer,intent(in)                              ::  i,n            !   atom i is being computed. n neighbours
                integer,dimension(0:),intent(in)                ::      t
                real(kind=real64),dimension(:),intent(in)       ::  r              !   (1:n) distance of atom i to each neighbour
                real(kind=real64),dimension(:,:),intent(in)     ::  x              !   (1:3,1:n) unit vector to each neighbour
                integer,dimension(:),intent(in)                 ::  neigh          !   (1:n) index of each neighbour
                real(kind=real32),dimension(:,:,:,:),pointer,intent(inout)      ::      hess        !   (1:3*bandwidth, 1:3*nAtoms)
                integer,dimension(:,:),pointer,intent(inout)                    ::      indx        !   (0:bandwidth,1:nAtoms)           with 0 = number of neighbours, 1: = neighbour index
     
                
                real(kind=real64)                   ::  imodrij,rho
                integer                             ::  jj,kk,ll,mm
                real(kind=real64)                   ::  VV,dVdr,d2Vdr2
                real(kind=real64)                   ::  phi,dphidr,d2phidr2
                real(kind=real64)                   ::  FF,dFdrho,d2Fdrho2
                real(kind=real64),dimension(3)      ::  dr1
                real(kind=real64),dimension(3,3)    ::  dr2
                real(kind=real32),dimension(3,3)    ::  hs
                real(kind=real64),dimension(3,0:n)  ::  fi
     
                real(kind=real32),dimension(3,3,0:n)  ::  hsi
                
                
                
        !---   quick return if no atoms in the neighbour list
            if (n==0) return
            
            rho = 0.0d0
            do kk = 1,n
                rho = rho + getPhi( this,t(kk),r(kk) )
            end do
            call getd2Fdrho2(this,t(0),rho,FF,dFdrho,d2Fdrho2)
  
 
                  
                hsi = 0.0
     
                fi = 0.0d0         ! = sum_j dphi_{ij}/dx_i = drho_i/dx_i
                do kk = 1,n
 
                    call getd2Vdr2( this,t(0),t(kk),r(kk),VV,dVdr,d2Vdr2)                                
                    call getd2Phidr2( this,t(kk),r(kk),phi,dphidr,d2phidr2)             
 
                    imodrij = 1.0d0/r(kk)                  
                    dr1(1:3) = -x(1:3,kk)
     
                    fi(1:3,0) = fi(1:3,0) + dphidr*dr1(1:3)         ! = sum_j dphi_{ij}/dx_i = drho_i/dx_i
                    fi(1:3,kk) =  - dphidr*dr1(1:3)                 ! = dphi_{ij}/dx_j = drho_i/dx_j
     
                    dr2(1:3,1:3) = d2rijbydxidxi( x(1:3,kk),imodrij )
     
    !                
                    hs(1:3,1) = real( (d2Vdr2/2+dFdrho*d2phidr2)*dr1(1:3)*dr1(1) + (dVdr/2+dFdrho*dphidr)*dr2(1:3,1) , kind=real32)
                    hs(1:3,2) = real( (d2Vdr2/2+dFdrho*d2phidr2)*dr1(1:3)*dr1(2) + (dVdr/2+dFdrho*dphidr)*dr2(1:3,2) , kind=real32)
                    hs(1:3,3) = real( (d2Vdr2/2+dFdrho*d2phidr2)*dr1(1:3)*dr1(3) + (dVdr/2+dFdrho*dphidr)*dr2(1:3,3) , kind=real32)
                    jj = neigh(kk)
                    
                    
                 
                    call addSparseHessian3x3( hess,indx,jj,jj,hs )
                    call addSparseHessian3x3( hess,indx,jj,i,-hs )
                 
                                                
                    
                    hsi(1:3,1:3,0) = hsi(1:3,1:3,0)   + hs(1:3,1:3)
                    hsi(1:3,1:3,kk) = hsi(1:3,1:3,kk) - hs(1:3,1:3)
                          
                end do
      
             
                
                hs(1:3,1) = real( d2Fdrho2*fi(1:3,0)*fi(1,0) , kind=real32 )         !   = d2E/drho_i2 * drho_i/dx_i * drho_i/dx_i
                hs(1:3,2) = real( d2Fdrho2*fi(1:3,0)*fi(2,0) , kind=real32 )
                hs(1:3,3) = real( d2Fdrho2*fi(1:3,0)*fi(3,0) , kind=real32 )
                 
                hsi(1:3,1:3,0) = hsi(1:3,1:3,0) + hs(1:3,1:3) 
     
                do kk = 1,n
                    jj = neigh(kk)
     
                    hs(1:3,1) = real( d2Fdrho2*fi(1:3,0)*fi(1,kk) , kind=real32 )        !    = d2E/drho_i2 * drho_i/dx_i * drho_i/dx_j
                    hs(1:3,2) = real( d2Fdrho2*fi(1:3,0)*fi(2,kk) , kind=real32 )
                    hs(1:3,3) = real( d2Fdrho2*fi(1:3,0)*fi(3,kk) , kind=real32 )
                     
                    hsi(1:3,1:3,kk) = hsi(1:3,1:3,kk)   + hs(1:3,1:3) 
                    call addSparseHessian3x3( hess,indx,jj,i,transpose(hs) )
                    do ll = 1,n
                        mm = neigh(ll)
                        hs(1:3,1) = real( d2Fdrho2*fi(1:3,ll)*fi(1,kk) , kind=real32 )      !    = d2E/drho_i2 * drho_i/dx_m * drho_i/dx_j
                        hs(1:3,2) = real( d2Fdrho2*fi(1:3,ll)*fi(2,kk) , kind=real32 )
                        hs(1:3,3) = real( d2Fdrho2*fi(1:3,ll)*fi(3,kk) , kind=real32 )
      
                        call addSparseHessian3x3( hess,indx, mm,jj,hs )                 !     note (mm,jj) and (jj,mm) are both visited in this loop
     
                    end do
                end do
     
                call addSparseHessian3x3xm( hess,indx,i,n,neigh,hsi )                    
                
                
                return
            end subroutine addHessian1a         
         
    !---
        subroutine hessianEinstein1a(this, i,n,t,r,x,neigh, hess )
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      Calculates the on-site (Einstein oscillator) Hessian ( force constant matrix) at a single atom i given its neighbour list
    !*      on input x are unit vectors with length r
    !*      note that this works for a single atom species only- otherwise the types need to be taken into account
    !*      note this is the dense version with a full 3Nx3N matrix
    !*      System potential energy is 
    !*             E = sum_i  1/2 sum_j V(r_ij) + F[ sum_j phi(r_ij) ] 
    
    
    
            type(EAM_Alloy), intent(in)                     ::  this
            integer,intent(in)                              ::  i,n            !   atom i is being computed. n neighbours
            integer,dimension(0:),intent(in)                ::  t               !   t(0) is type of atom i
            real(kind=real64),dimension(:),intent(in)       ::  r              !   (1:n) distance of atom i to each neighbour
            real(kind=real64),dimension(:,:),intent(in)     ::  x              !   (1:3,1:n) unit vector to each neighbour
            integer,dimension(:),intent(in)                 ::  neigh          !   (1:n) index of each neighbour
            
            real(kind=real32),dimension(:,:,:),intent(out)    ::  hess            !   3x3xn
            
            real(kind=real64)                   ::  imodrij,rho
            integer                             ::  jj,kk
            real(kind=real64)                   ::  VV,dVdr,d2Vdr2
            real(kind=real64)                   ::  phi,dphidr,d2phidr2
            real(kind=real64)                   ::  FF,dFdrho,d2Fdrho2
            real(kind=real64),dimension(3)      ::  dr1
            real(kind=real64),dimension(3,3)    ::  dr2
            real(kind=real32),dimension(3,3)    ::  hs,hsi
            real(kind=real64),dimension(3,0:n)  ::  fi
  
                
                
        !---   quick return if no atoms in the neighbour list
            if (n==0) return
            
            rho = 0.0d0
            do kk = 1,n
                rho = rho + getPhi( this,t(kk),r(kk) )
            end do
            call getd2Fdrho2(this,t(0),rho,FF,dFdrho,d2Fdrho2)
  
 
              
            hsi = 0.0 
            fi = 0.0d0         ! = sum_j dphi_{ij}/dx_i = drho_i/dx_i
            
            do kk = 1,n

                call getd2Vdr2( this,t(0),t(kk),r(kk),VV,dVdr,d2Vdr2)                                
                call getd2Phidr2( this,t(kk),r(kk),phi,dphidr,d2phidr2)             

                imodrij = 1.0d0/r(kk)                  
                dr1(1:3) = -x(1:3,kk)
 
                fi(1:3,0) = fi(1:3,0) + dphidr*dr1(1:3)         ! = sum_j dphi_{ij}/dx_i = drho_i/dx_i
                fi(1:3,kk) =  - dphidr*dr1(1:3)                 ! = dphi_{ij}/dx_j = drho_i/dx_j
 
                dr2(1:3,1:3) = d2rijbydxidxi( x(1:3,kk),imodrij )
 
!                
                hs(1:3,1) = real( (d2Vdr2/2+dFdrho*d2phidr2)*dr1(1:3)*dr1(1) + (dVdr/2+dFdrho*dphidr)*dr2(1:3,1) , kind=real32)
                hs(1:3,2) = real( (d2Vdr2/2+dFdrho*d2phidr2)*dr1(1:3)*dr1(2) + (dVdr/2+dFdrho*dphidr)*dr2(1:3,2) , kind=real32)
                hs(1:3,3) = real( (d2Vdr2/2+dFdrho*d2phidr2)*dr1(1:3)*dr1(3) + (dVdr/2+dFdrho*dphidr)*dr2(1:3,3) , kind=real32)
                
                hsi(1:3,1:3) = hsi(1:3,1:3)   + hs(1:3,1:3)
                      
            end do
  
         
            
            hs(1:3,1) = real( d2Fdrho2*fi(1:3,0)*fi(1,0) , kind=real32 )         !   = d2E/drho_i2 * drho_i/dx_i * drho_i/dx_i
            hs(1:3,2) = real( d2Fdrho2*fi(1:3,0)*fi(2,0) , kind=real32 )
            hs(1:3,3) = real( d2Fdrho2*fi(1:3,0)*fi(3,0) , kind=real32 )
             
            hsi(1:3,1:3) =  hsi(1:3,1:3) + hs(1:3,1:3) 
   
            hess(1:3,1:3,i) =  hess(1:3,1:3,i) + hsi(1:3,1:3)  
            

            do kk = 1,n
                jj = neigh(kk)
 
!                hs(1:3,1) = d2Fdrho2*fi(1:3,0)*fi(1,kk)        !    = d2E/drho_i2 * drho_i/dx_i * drho_i/dx_j
!                hs(1:3,2) = d2Fdrho2*fi(1:3,0)*fi(2,kk)
!                hs(1:3,3) = d2Fdrho2*fi(1:3,0)*fi(3,kk)
!                 
!                hsi(1:3,1:3,kk) = hsi(1:3,1:3,kk)   + hs(1:3,1:3) 
!                call addSparseHessian3x3( hess,indx,jj,i,transpose(hs) )
!                do ll = 1,n
!                    mm = neigh(ll)
                    hs(1:3,1) = real( d2Fdrho2*fi(1:3,kk)*fi(1,kk) , kind=real32 )     !    = d2E/drho_i2 * drho_i/dx_m * drho_i/dx_j
                    hs(1:3,2) = real( d2Fdrho2*fi(1:3,kk)*fi(2,kk) , kind=real32 )
                    hs(1:3,3) = real( d2Fdrho2*fi(1:3,kk)*fi(3,kk) , kind=real32 )
  
                 hess(1:3,1:3,jj) =  hess(1:3,1:3,jj) + hs(1:3,1:3)  
                       
                 !call addSparseHessian3x3( hess,indx, mm,jj,hs )                 !     note (mm,jj) and (jj,mm) are both visited in this loop
 
                !end do
            end do
 
            !call addSparseHessian3x3xm( hess,indx,i,n,neigh,hsi )                    
                        
            
            
            return
        end subroutine hessianEinstein1a         
         


        pure function d2rijbydxidxi( xij,imodrij ) result( d2r )
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      compute d2/d x_i d x_i r_ij = d/x_i (x_i - x_j)/rij
    !*                                       = I/rij - (x_i - x_j)(x_i - x_j)/rij^3
    !*                                       = - d2rijbydxidxj
    !*      on entry xij = rij/|rij|
            real(kind=real64),dimension(3),intent(in)        ::      xij
            real(kind=real64),intent(in)                     ::      imodrij
            real(kind=real64),dimension(3,3)                 ::      d2r
            d2r(:,1) = (/ imodrij,0.0d0,0.0d0 /) - xij(:)*xij(1)*imodrij
            d2r(:,2) = (/ 0.0d0,imodrij,0.0d0 /) - xij(:)*xij(2)*imodrij
            d2r(:,3) = (/ 0.0d0,0.0d0,imodrij /) - xij(:)*xij(3)*imodrij
            return
        end function d2rijbydxidxi


!-------

        subroutine inputFromXML0(xml,this,ok)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      read in a single element potential from xml
    !*                                                  
    !*          <EAM_PotentialTable E0="in ev" name="element" mass="in daltons">
    !*             <QuinticSplineEven type="F"/> 
    !*             <QuinticSplineEven type="PHI"/>
    !*             <QuinticSplineEven type="V"/>
    !*          </EAM_PotentialTable>
    !*
    !*
    
            type(EAM_PotentialTable),intent(out)    ::      this
            type(NBAX),intent(in)                   ::      xml
            logical,intent(out)                     ::      ok

            type(NBAX),pointer                      ::      xmlp
            character(len=256)                      ::      qspline_type
            logical                                 ::      allok
            integer                                 ::      ii
            type(QuinticSplineEven)     ::      q_phi,q_V,q_F
            real(kind=real64)           ::      E0,mass
            character(len=EAM_NAME_LEN) ::      name

            if (.not. (xml == "EAM_PotentialTable")) then
                print *,"EAM_PotentialTables::inputFromXML0 warning - xml is not an ""EAM_PotentialTable"""
                ok = .false.
                return
            end if

            allok = .true.
            E0 = 0.0d0 ; call getAttributeValue(xml,"E0",E0,ok) ; allok = allok .and. ok
            if (.not.ok) write(unit=0,fmt='(a,i6)') "EAM_PotentialTables::inputFromXML0 error - could not read attribute ""E0"""
            mass = 1.0d0 ; call getAttributeValue(xml,"mass",mass,ok)
            name = "" ; call getAttributeValue(xml,"name",name,ok) ; allok = allok .and. ok
            if (.not.ok) write(unit=0,fmt='(a,i6)') "EAM_PotentialTables::inputFromXML0 error - could not read attribute ""name"""
            do ii = 1,3
                call getChild(xml,"QuinticSplineEven",ii,xmlp,ok) ; allok = allok .and. ok
                if (.not. ok) write(unit=0,fmt='(a,i6)') "EAM_PotentialTables::inputFromXML0 error - could not read child ""QuinticSplineEven"" ",ii
                call getAttributeValue(xmlp,"type",qspline_type,ok) ; allok = allok .and. ok
                if (.not. ok) write(unit=0,fmt='(a,i6)') "EAM_PotentialTables::inputFromXML0 error - could not read attribute ""type"" for child ",ii
                if ( (trim(qspline_type)=="phi").or.(trim(qspline_type)=="PHI") ) then
                    print *,"EAM_PotentialTables::inputFromXML0 info - reading phi()"
                    call inputFromXML(xmlp,q_phi,ok) ; allok = allok .and. ok
                else if ( (trim(qspline_type)=="v").or.(trim(qspline_type)=="V") ) then
                    print *,"EAM_PotentialTables::inputFromXML0 info - reading V()"
                    call inputFromXML(xmlp,q_V,ok) ; allok = allok .and. ok
                else if ( (trim(qspline_type)=="f").or.(trim(qspline_type)=="F") ) then
                    print *,"EAM_PotentialTables::inputFromXML0 info - reading F()"
                    call inputFromXML(xmlp,q_F,ok) ; allok = allok .and. ok
                else                    
                    allok = .false.
                    write(unit=0,fmt='(a,i6,a)') "EAM_PotentialTables::inputFromXML0 error - expected attribute ""type"" for child ",ii," to be one of (phi,V,F) , read "//trim(qspline_type)
                end if
                if (.not.ok) write(unit=0,fmt='(a,i6)') "EAM_PotentialTables::inputFromXML0 error - could not read child QuinticSplineEven::"//trim(qspline_type)
            end do

            ! DEBUG
            ! call clear(q_V)
            
            
            this = EAM_PotentialTable_ctor( E0,mass*EAM_MASS_UNIT,q_phi,q_V,q_F,name )
            
            print *,"EAM_PotentialTables::inputFromXML0 info - EAM_PotentialTable loaded from xml ",allok
            ok = allok
            if (.not. ok) then
                print *,"EAM_PotentialTables::inputFromXML0 WARNING - unable to read EAM_PotentialTable from xml"
            end if
            
            
            
            
            return
        end subroutine inputFromXML0


        subroutine outputAsXML0(this,xml)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
            type(EAM_PotentialTable),intent(in)     ::      this
            type(NBAX),intent(out)                  ::      xml

            type(NBAX)                              ::      xmlp1,xmlp2,xmlp3

            
            xml = NBAX_ctor("EAM_PotentialTable")
!             call addAttribute(xml,"A0",this%A0)
            call addAttribute(xml,"E0",this%E0)
            call addAttribute(xml,"mass",this%mass/EAM_MASS_UNIT)
            call addAttribute(xml,"name",trim(this%name))


            xmlp1 = NBAX_ctor("QuinticSplineEven")
            call outputAsXML(this%q_phi,xmlp1)
            call addAttribute(xmlp1,"type","phi")
            call addChild(xml,xmlp1)

            xmlp2 = NBAX_ctor("QuinticSplineEven")
            call outputAsXML(this%q_V,xmlp2)
            call addAttribute(xmlp2,"type","V")
            call addChild(xml,xmlp2)

            xmlp3 = NBAX_ctor("QuinticSplineEven")
            call outputAsXML(this%q_F,xmlp3)
            call addAttribute(xmlp3,"type","F")
            call addChild(xml,xmlp3)

            return
        end subroutine outputAsXML0

    !---


        subroutine inputFromXML1(xml,this,ok)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*          reads in an alloy potential from xml
    !*          <EAM_Alloy n="number of elements" name="element1,element2,element3..." E0="eV1,eV2,eV3..." mass="m1,m2,m3..." >
    !*              <QuinticSplineEven type="F"/>  
    !*              <QuinticSplineEven type="PHI"/>
    !*              <QuinticSplineEven type="F"/>  
    !*              <QuinticSplineEven type="PHI"/>
    !*              <QuinticSplineEven type="F"/>  
    !*              <QuinticSplineEven type="PHI"/>
    !*              ...
    !*              <QuinticSplineEven type="V"/>          !   11
    !*              <QuinticSplineEven type="V"/>          !   21
    !*              <QuinticSplineEven type="V"/>          !   22
    !*              <QuinticSplineEven type="V"/>          !   31
    !*              <QuinticSplineEven type="V"/>          !   32
    !*              <QuinticSplineEven type="V"/>          !   33
    !*              ... 
    !*          </EAM_Alloy>
    !*
    !*
    
            type(EAM_Alloy),intent(out)             ::      this
            type(NBAX),intent(in)                   ::      xml
            logical,intent(out)                     ::      ok

            type(NBAX),pointer                      ::      xmlp
            logical                                 ::      allok
            integer                                 ::      kk,nn
 
            
            real(kind=real64),dimension(:),allocatable          ::      E0          !   energy scale,mass
            real(kind=real64),dimension(:),allocatable          ::      mass        !   mass
            character(len=16),dimension(:),allocatable          ::      name
            type(QuinticSplineEven),dimension(:),allocatable    ::      q_phi
            type(QuinticSplineEven),dimension(:,:),allocatable  ::      q_V
            type(QuinticSplineEven),dimension(:),allocatable    ::      q_F
            character(len=256)                      ::      qspline_type

            character(len=256)          ::      value
             
            integer                     ::      ii,jj,nPhi,nF
           
            
            if (.not. (xml == "EAM_Alloy")) then
                print *,"EAM_PotentialTables::inputFromXML1 warning - xml is not an ""EAM_Alloy"""
                ok = .false.
                return
            end if

            call getAttributeValue(xml,"n",nn,ok)
            if (.not. ok) stop "EAM_PotentialTables::inputFromXML1 error - could not read attribute ""n"""
            
            
            allocate(name(nn))
            allocate(E0(nn))
            allocate(mass(nn))
            allocate(q_phi(nn))
            allocate(q_v(nn,nn))
            allocate(q_f(nn))

            
            call getAttributeValue(xml,"E0",value,ok)
            if (.not. ok) stop "EAM_PotentialTables::inputFromXML1 error - could not read attribute ""E0"""
            ii = nn ; call parse(value,E0,ii)
            
            call getAttributeValue(xml,"name",value,ok)
            if (.not. ok) stop "EAM_PotentialTables::inputFromXML1 error - could not read attribute ""name"""
            ii = nn ; call parse(value,name,ii)
            
            call getAttributeValue(xml,"mass",value,ok)
            if (.not. ok) stop "EAM_PotentialTables::inputFromXML1 error - could not read attribute ""mass"""
            ii = nn ; call parse(value,mass,ii)
            
            if (getNChildren(xml,"QuinticSplineEven") /= 2*nn + (nn*(nn+1))/2 )  stop "EAM_PotentialTables::inputFromXML1 error - not enough splines"
            nPhi = 0 ; nF = 0
            
            do ii = 1,nn
                do jj = 1,nn
                    q_v(jj,ii) = QuinticSplineEven_ctor()
                end do
            end do
            
            ii = 0 ; jj = 0; nPhi = 0 ; nF = 0 ; allok = .true.
           
            do kk = 1,2*nn + (nn*(nn+1))/2
                call getChild(xml,"QuinticSplineEven",kk,xmlp,ok)      
                call getAttributeValue(xmlp,"comment",value,ok)                 
                call getAttributeValue(xmlp,"type",qspline_type,ok) ; allok = allok .and. ok                
                if (.not. ok) write(unit=0,fmt='(a,i6)') "EAM_PotentialTables::inputFromXML1 error - could not read attribute ""type"" for child ",kk
                
                if ( (trim(qspline_type)=="phi").or.(trim(qspline_type)=="PHI") ) then
                    nPhi = nPhi + 1
                    write (*,fmt='(a,i2,a)') "EAM_PotentialTables::inputFromXML1 info - reading phi(",nPhi,") """//trim(value)//""""
                    call inputFromXML(xmlp,q_phi(nPhi),ok) ; allok = allok .and. ok
                else if ( (trim(qspline_type)=="f").or.(trim(qspline_type)=="F") ) then
                    nF = nF + 1
                     write (*,fmt='(a,i2,a)') "EAM_PotentialTables::inputFromXML1 info - reading F(",nF,") """//trim(value)//""""
                    call inputFromXML(xmlp,q_F(nF),ok) ; allok = allok .and. ok
                else if ( (trim(qspline_type)=="v").or.(trim(qspline_type)=="V") ) then
                    jj = jj + 1
                    if (jj>ii) then
                        jj = 1
                        ii = ii + 1
                    end if
                    write (*,fmt='(a,i2,a,i2,a)')  "EAM_PotentialTables::inputFromXML1 info - reading V(",ii,",",jj,") """//trim(value)//""""
                    call inputFromXML(xmlp,q_V(ii,jj),ok) ; allok = allok .and. ok
                    if (ii /= jj) call clone( q_v(jj,ii),q_v(ii,jj) )
                else                    
                    allok = .false.
                    write(unit=0,fmt='(a,i6,a)') "EAM_PotentialTables::inputFromXML1 error - expected attribute ""type"" for child ",kk," to be one of (phi,V,F) , read "//trim(qspline_type)
                end if
                if (.not.ok) write(unit=0,fmt='(a,i6)') "EAM_PotentialTables::inputFromXML1 error - could not read child QuinticSplineEven::"//trim(qspline_type)
            end do
                
            ok = ok .and. allok
            if (.not. ok) print *,"EAM_PotentialTables::inputFromXML1 WARNING - could not read potential"
              
            this = EAM_Alloy_ctor0( E0,name,mass*EAM_MASS_UNIT,q_phi,q_V,q_f )     
            
            
            
            do ii = 1,nn
                call delete(q_phi(ii))
                call delete(q_F(ii))
                do jj = 1,nn
                    call delete(q_V(jj,ii))
                end do
            end do
            deallocate(name)
            deallocate(E0)
            deallocate(mass)
            deallocate(q_phi)
            deallocate(q_v)
            deallocate(q_f)
            
            

            
             
            return
        end subroutine inputFromXML1


!       subroutine outputAsXML1(this,xml)
!   !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
!           type(EAM_Alloy),intent(in)              ::      this
!           type(NBAX),intent(out)                  ::      xml
!
!           type(NBAX)                              ::      xmlp
!           character(len=1024)                     ::      names
!           integer                 ::      ii,jj
!
!           xml = NBAX_ctor("EAM_Alloy")
!           call addAttribute(xml,"n",this%n)
!           names = ""
!           do ii = 1,this%n
!               names = trim(names)//" "//getName(this,ii)
!           end do
!           call addAttribute(xml,"names",trim(names))
!           do ii = 1,this%n
!               do jj = ii,this%n
!                   call outputAsXML(this%eam(ii,jj),xmlp)
!                   call addChild(xml,xmlp)
!               end do
!           end do
!
!           return
!       end subroutine outputAsXML1



    !---

        pure function getName0(this) result(name)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
            type(EAM_PotentialTable),intent(in)         ::      this
            character(len=EAM_NAME_LEN)                 ::      name
            name = this%name
            return
        end function getName0


        function getName1(this,ti) result(name)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      return the first half of the name
            type(EAM_Alloy),intent(in)                  ::      this
            integer,intent(in)                          ::      ti
            character(len=EAM_NAME_LEN)                 ::      name
            name = this%name(ti)
            return
        end function getName1


        pure subroutine setName0(this,name)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
            type(EAM_PotentialTable),intent(inout)      ::      this
            character(len=*),intent(in)      ::      name
            this%name = name
            return
        end subroutine setName0
! 
!         pure subroutine setName1(this,name1,name2)
!     !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
!             type(EAM_PotentialTable),intent(inout)      ::      this
!             character(len=*),intent(in)      ::      name1,name2
!             this%name = trim(name1(1:16))//":"//trim(name2(1:16))
!             return
!         end subroutine setName1


        pure function getMass0(this) result(mass)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
            type(EAM_PotentialTable),intent(in)         ::      this
            real(kind=real64)                           ::      mass
            mass = this%mass
            return
        end function getMass0


        pure function getMass1(this,ti) result(mass)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
            type(EAM_Alloy),intent(in)                  ::      this
            integer,intent(in)                          ::      ti
            real(kind=real64)                           ::      mass
            mass = this%mass(ti)
            return
        end function getMass1


    end module EAM_PotentialTables

!
!   program testEAM_PotentialTables
! !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
!
!       use Lib_NBAX
!       use EAM_PotentialTables
!       use Lib_Timers
!       use iso_fortran_env
!       implicit none
!
!       type(NBAX)                  ::      xml
!       type(EAM_PotentialTable)    ::      p
!
!       integer                     ::      ii,jj,kk
!       real(kind=real64)           ::      rr,phi,dphi,d2phi,V,dV,d2V,rho,F,dF,d2F
!       logical                     ::      ok
!
!
!       integer,parameter           ::      NN = 5
!       integer,parameter           ::      NATOMS = 2*NN*NN*NN
!       integer         ::      ix,iy,iz,ik
!       real(kind=real64),dimension(3,NATOMS)       ::  x
!       type(Timer)                                 ::      tt
!       integer                     ::      NX
!
!       real(kind=real64),dimension(1:14,NATOMS)    ::  rij
!       real(kind=real64),dimension(1:3,1:14,NATOMS)    ::  xij
!       integer,dimension(0:14,NATOMS)  ::  neigh
!       real(kind=real64),dimension(3)          ::  dx
!       real(kind=real64),dimension(1:3,NATOMS) ::  force
!       real(kind=real64),dimension(1:100)  ::  work
!
!       xml = NBAX_ctor()
!       open(unit=400,file="Data/Potentials/EAM_W.xml",action="read")
!           call input(xml,400,ok)
!           print *,"input from disk ",ok
!           call inputFromXML(xml,p,ok)
!           print *,"input from xml ",ok
!       close(unit=400)
!       call delete(xml)
!
!       call report(p)
!
!       do ii = 0,50
!           rr = ii*5.0d0/50
!           call getd2Phidr2(p,rr, phi,dphi,d2phi)
!           call getd2Vdr2(p,rr, V,dV,d2V)
!
!           rho = ii*50.0d0/50
!           call getd2Fdrho2(p,rho, F,dF,d2F)
!           write (*,fmt='(i6,100g24.16)') ii,rr,phi,dphi,d2phi,V,dV,d2V,rho,F,dF,d2F
!       end do
!
!
!       call outputAsXML(p,xml)
!       open(unit=500,file="test.xml",action="write")
!           call output(xml,500,.true.)
!       close(unit=500)
!
!   !---    create chunk of crystal
!       ii = 0
!       do iz = 0,NN-1
!           do iy = 0,NN-1
!               do ix = 0,NN-1
!                   do ik = 0,1
!                       ii = ii + 1
!                       x(:,ii) = (/ ix,iy,iz /) + 0.5d0*ik
!                   end do
!               end do
!           end do
!       end do
!
!   !---    construct neighbour list
!       neigh = 0
!       rij = 0.0d0
!       do ii = 1,NATOMS-1
!           do jj = ii+1,NATOMS
!               dx = x(:,jj) - x(:,ii)
!               do kk = 1,3
!                   if (2*dx(kk) > NN) dx(kk) = dx(kk)-NN
!                   if (2*dx(kk) <-NN) dx(kk) = dx(kk)+NN
!               end do
!               rr = dot_product( dx,dx )
!               if (rr < 1.5d0) then
!                   rr = sqrt(rr) ; dx = dx/rr
!                   kk = neigh(0,ii) + 1 ; neigh(0,ii) = kk ; neigh(kk,ii) = jj ; xij(1:3,kk,ii) = dx ; rij(kk,ii) = rr
!                   kk = neigh(0,jj) + 1 ; neigh(0,jj) = kk ; neigh(kk,jj) = ii ; xij(1:3,kk,jj) = -dx ; rij(kk,jj) = rr
!               end if
!           end do
!       end do
!   !---    scale to latt param
!       rij = rij * 3.1652d0
!
!
!       call gaugeTransform(p,rij(:,1))
!       call report(p)
!
!
!
!
!       dF = 0.0d0
!       NX = 10000
!       tt = Timer_ctor()
!       do jj = 1,NX
!           do ii = 1,NATOMS
!               kk = neigh(0,ii)
!               F = getEnergy( p,rij(1:kk,ii) )
!               dF = dF + F
!           end do
!       end do
!       print *,"elapsed (energy) ",real(elapsed(tt))/(NX*NATOMS)
!       print *,"energy (avoid opt)",dF/(NX*NATOMS)
!
!
!       dF = 0.0d0
!       NX = 10000
!       tt = Timer_ctor()
!       do jj = 1,NX
!           force = 0.0d0
!           do ii = 1,NATOMS
!               kk = neigh(0,ii)
!               call addForce( p,ii,kk,rij(:,ii),xij(:,:,ii),neigh(1:,ii), force ,work)
!           end do
!           do ii = 1,3
!               dx(ii) = sum(force(ii,:))
!           end do
!       end do
!       print *,"elapsed (force) ",real(elapsed(tt))/(NX*NATOMS)
!       print *,"force (avoid opt)",dx/(NX)
!
!
!       call delete(xml)
!
!
!
!       print *,""
!       print *,"done"
!       print *,""
!
!
!   end program testEAM_PotentialTables