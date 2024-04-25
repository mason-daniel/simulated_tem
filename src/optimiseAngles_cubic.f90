!   gfortran -O2 -ffree-line-length-256 ${MYF90LIB}/NBAX_StringTokenizers.f90 ${MYF90LIB}/Lib_CommandLineArguments.f90 ${MYF90LIB}/Lib_RotationMatrices.f90 optimiseAngles_cubic.f90 -o optimiseAngles.exe 

    program optimiseAngles
!---^^^^^^^^^^^^^^^^^^^^^^
        use iso_fortran_env
        use Lib_RotationMatrices
        use Lib_CommandLineArguments
        implicit none
        
        type(CommandLineArguments)      ::      cla       
        integer                         ::      NJ = 3
        integer                         ::      freq = 1
        character(len=256)              ::      outfile = "test"
        character(len=256)              ::      infile = ""
        integer                         ::      maxSteps = 1000
        real(kind=real64)               ::      kick = 0.0d0
        integer                         ::      firstnOut = 0
        
        real(kind=real64),parameter     ::      beta = 0.01d0
        integer,parameter               ::      NS = 1000 
        real(kind=real64)               ::      DT 
        
        
        character(len=4)        ::      aaaa , bbbb
        real(kind=real64),dimension(:,:,:),allocatable      ::      R_in
        real(kind=real64),dimension(:,:,:),allocatable      ::      R
        real(kind=real64),dimension(:,:),allocatable        ::      v
        integer                 ::      ii,jj,nj_in
        real(kind=real64)       ::      ee,tt , elast , eold
        logical                 ::      ok
        integer                 ::      nout,firsti
        real(kind=real64),dimension(3)  ::      zeta
        
        
        cla = CommandLineArguments_ctor(10)
        call setProgramDescription( cla, "optimiseAngles.exe" )
        call setProgramVersion( cla, "1.0" )
        
        call get( cla,"n",NJ ,LIB_CLA_REQUIRED,"number of jacks" )          
        call get( cla,"o",outfile ,LIB_CLA_OPTIONAL," output file prefix - will write to PREFIX.NJACKS.STEP.xyz" )  
        call get( cla,"i",infile ,LIB_CLA_OPTIONAL," input file" )  
        call get( cla,"f",freq ,LIB_CLA_OPTIONAL," output frequency ( 0 for off )" )          
        call get( cla,"steps",maxSteps ,LIB_CLA_OPTIONAL," maximum steps " )          
        call get( cla,"kick",kick ,LIB_CLA_OPTIONAL," random kick Euler angles by maximum this " )  
        call get( cla,"nOut",firstnOut ,LIB_CLA_OPTIONAL," first output number" )  
        
        
        call report(cla)
        if (.not. allRequiredArgumentsSet(cla)) stop
        if (hasHelpArgument(cla)) stop
        call delete(cla)
        
                 
       allocate(R(3,3,NJ)) 
       allocate(v(3,NJ)  ) 
       write(bbbb,fmt='(i4)') NJ ; bbbb = adjustl(bbbb) ; bbbb = repeat("0",4-len_trim(bbbb))//trim(bbbb)
       
       
        if (len_trim(infile)/=0) then
            inquire(file=trim(infile),exist=ok)
            if (.not. ok) stop "input file not found"
            open(unit=500,file=trim(infile),action="read")
                call inputXYZ(500,R_in)                 
            close(unit=500)
            nj_in = size(R_in,dim=3)
            if (nj_in < NJ) then
                R(1:3,1:3,1:nj_in) = R_in(1:3,1:3,1:nj_in)
                call initRotMats(  R(1:3,1:3,nj_in+1:NJ)  )   
            else
            !   should really randomly select jacks to drop. This will do to avoid a fail.
                R(1:3,1:3,1:NJ) = R_in(1:3,1:3,1:NJ)
            end if
            deallocate(R_in)
            
        else
            R(1:3,1,1) = (/ 1.0d0,0.0d0,0.0d0 /)
            R(1:3,2,1) = (/ 0.0d0,1.0d0,0.0d0 /)
            R(1:3,3,1) = (/ 0.0d0,0.0d0,1.0d0 /)
            call initRotMats(  R(1:3,1:3,2:NJ)  )   
        end if
        
        do ii = 1,size(R,dim=3)
            call tidyRotationMatrix ( R(1:3,1:3,ii) )
        end do
        
        !call outputXyz(6,R)
        v = 0.0d0
        
        if ( (len_trim(infile)/=0).and.(nj_in<NJ) ) then    
            elast = huge(1.0)        
            do ii = 1,100
                DT = 10.0**(log10(real(ii))-3*log(real(NJ)))
                call dampedDynamicsJacks( R,v,beta,NS,DT,nj_in )         
                
                ee = totalEnergy(R)           
                tt = tt + DT
                if (ee>elast) v = v/2
                elast = ee
                    
                print *,"init ",ii," time ",tt," energy ",ee  , sum(v*v)         
                  
            end do
        end if    
        
        if (kick>0) then
            print *,"kicking Euler angles ",kick
            do ii = 2,NJ
                call random_number(zeta) 
                zeta = (2*zeta-1.0d0)*kick                
                R(1:3,1:3,ii) = matmul( rotmat(zeta(1),zeta(2),zeta(3)),R(1:3,1:3,ii) )
            end do
        end if
        
        
        tt = 0.0d0
        elast = huge(1.0) ; eold = huge(1.0)
        nout = firstnOut
        if ( (freq > 0).and.(nout==0) ) then            
            open(unit=600,file=trim(outfile)//bbbb//".0000.xyz",action="write")
                call outputXyz(600,R)
            close(unit=600)
        end if          
        
        if (freq>0) firsti = firstnOut * freq
        do ii = firsti+1,firsti+maxSteps
            DT =  10.0**(log10(real(ii))-3*log(real(NJ)))
             
            if (freq > 0) then
                if (mod(ii,freq)==0) then
                    nout = nout + 1
                    write(aaaa,fmt='(i4)') nout ; aaaa = adjustl(aaaa) ; aaaa = repeat("0",4-len_trim(aaaa))//trim(aaaa)
                    open(unit=600,file=trim(outfile)//bbbb//"."//aaaa//".xyz",action="write")
                        call outputXyz(600,R)
                    close(unit=600)
                    
                end if
            end if            
            call dampedDynamicsJacks( R,v,beta,NS,DT,1 )
            
            
            ee = totalEnergy(R)        
            tt = tt + DT
            if (ee>elast+1.0d-6) then
                v = v/2
                eold = huge(1.0)
            end if
            elast = ee
                
            print *,"step ",ii," time ",tt," energy ",ee  , sum(v*v)!,(eold - ee)         
            
            if (mod(ii,10)==0) then
                if (eold - ee < 1.0d-6) then
                    exit
                else
                    eold = ee
                end if
            end if
            
            
            
        end do    
    !---    always output final step
    
        do ii = 1,NJ
            call tidyRotationMatrix ( R(1:3,1:3,ii) )
        end do
        
        write(aaaa,fmt='(i4)') nout + 1; aaaa = adjustl(aaaa) ; aaaa = repeat("0",4-len_trim(aaaa))//trim(aaaa)
        open(unit=600,file=trim(outfile)//bbbb//"."//aaaa//".xyz",action="write")
            call outputXyz(600,R)               
        close(unit=600)

         
        open(unit=601,file=trim(outfile)//bbbb//".dat",action="write")
            write(unit=601,fmt=*) NJ,ee
            do ii = 1,NJ
                write(unit=601,fmt=*) R(:,:,ii)
            end do              
        close(unit=601)        
       ! print *,""
       ! call outputXyz(6,R)
        
        print *,""
        print *,"done"
        print *,""
        
    contains
!---^^^^^^^^        

        function totalEnergy(R) result(e)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
            real(kind=real64),dimension(:,:,:),intent(in)       ::      R
            real(kind=real64)                                   ::      e
            integer         ::      NJ,ii
            NJ = size(R,dim=3)
            e = 0.0d0
            do ii = 1,NJ
                e = e + jacksEnergy( R(1:3,1:3,ii),ii,R )
            end do    
            e = e/(2*NJ*(NJ-1))             
            return
        end function totalEnergy
            
    

            
        subroutine outputXyz(u,R)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^
            integer,intent(in)                                  ::      u
            real(Kind=real64),dimension(:,:,:),intent(in)       ::      R
            integer     ::      ii
            write(unit=u,fmt='(i12)') 6*NJ
            !write(unit=u,fmt='(a)')  "Boxsize 1.0 1.0 1.0"
            !write(unit=u,fmt='(a)')  "Lattice=""2.0 0.0 0.0 0.0 2.0 0.0 0.0 0.0 2.0"" Properties=species:S:1:pos:R:3:conc:R:1"
            write(unit=u,fmt='(a)')  "atom position_x position_y position_z charge"
            do ii = 1,NJ
                write(unit=u,fmt='(a,4f12.5)') "J ",R(1:3,1,ii) , real(ii)
                write(unit=u,fmt='(a,4f12.5)') "J ",R(1:3,2,ii) , real(ii)
                write(unit=u,fmt='(a,4f12.5)') "J ",R(1:3,3,ii) , real(ii)
                write(unit=u,fmt='(a,4f12.5)') "J ",-R(1:3,1,ii), real(ii)
                write(unit=u,fmt='(a,4f12.5)') "J ",-R(1:3,2,ii), real(ii)
                write(unit=u,fmt='(a,4f12.5)') "J ",-R(1:3,3,ii), real(ii)
            end do
            write(600,fmt='(a)') ""
            write(600,fmt='(a)') "# energy ",totalEnergy(R)    
            return
        end subroutine outputXyz
        
        subroutine inputXyz(u,R_in)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^
            integer,intent(in)                                          ::      u
            real(kind=real64),dimension(:,:,:),allocatable,intent(out)  ::      R_in
            integer             ::      ii,NJ
            character(len=256)  ::      dummy
            read(unit=u,fmt='(i12)') ii
            NJ = ii/6            
            allocate(R_in(3,3,NJ))
            read(unit=u,fmt='(a)') dummy
            do ii = 1,NJ
                read(unit=u,fmt=*) dummy(1:2),R_in(1:3,1,ii)  
                read(unit=u,fmt=*) dummy(1:2),R_in(1:3,2,ii)  
                read(unit=u,fmt=*) dummy(1:2),R_in(1:3,3,ii)  
                read(unit=u,fmt='(a)') dummy
                read(unit=u,fmt='(a)') dummy
                read(unit=u,fmt='(a)') dummy
            end do
            return
        end subroutine inputXyz
                       
        
        subroutine dampedDynamicsJacks( R,v,beta,n,dt,nStatic )
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ 
            real(kind=real64),dimension(:,:,:),intent(inout)        ::      R
            real(kind=real64),dimension(3,size(R,dim=3)),intent(inout)    ::      v
            real(kind=real64),intent(in)                            ::      beta
            integer,intent(in)                                      ::      n
            real(kind=real64),intent(in)                            ::      dt
            integer,intent(in)                                      ::      nStatic
            real(kind=real64),dimension(3,size(R,dim=3))    ::      ff
            real(kind=real64),dimension(3,3)                ::      rot,oldR
            
            real(kind=real64)   ::      thetax,thetay,thetaz , ee
            integer             ::      ii,jj
            
            do jj = 1,n
                call jackForce( R, ff, nStatic )
                ff = ff - beta*v 
            
                do ii = nStatic+1,size(R,dim=3)
                    thetax = v(1,ii)*dt
                    thetay = v(2,ii)*dt
                    thetaz = v(3,ii)*dt
            
                    rot = rotmat(thetax,thetay,thetaz)        
                 !   print *,determinant3Mat( R(1:3,1:3,ii) )
                    !R(1:3,1:3,ii) = matmul( rot,R(1:3,1:3,ii) )
                    oldR(1:3,1:3) = R(1:3,1:3,ii)
                    R(1:3,1,ii) = rot(1:3,1)*oldR(1,1) + rot(1:3,2)*oldR(2,1) + rot(1:3,3)*oldR(3,1)
                    R(1:3,2,ii) = rot(1:3,1)*oldR(1,2) + rot(1:3,2)*oldR(2,2) + rot(1:3,3)*oldR(3,2)
                    R(1:3,3,ii) = rot(1:3,1)*oldR(1,3) + rot(1:3,2)*oldR(2,3) + rot(1:3,3)*oldR(3,3)
                    
                 !   print *,determinant3Mat( R(1:3,1:3,ii) )
                    
                    v(1:3,ii) = v(1:3,ii) + ff(1:3,ii)*dt
                end do
                
            end do   
            
            
            
            return
        end subroutine dampedDynamicsJacks
                
                
    
            pure function determinant3Mat(M) result(d)
        !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
        !*      returns the determinant of M
                real(kind=real64),dimension(3,3),intent(in)      ::      M
                real(kind=real64)                                ::      d
                real(kind=real64),dimension(9)       ::      dd
                dd(1) = M(1,1)*( M(2,2)*M(3,3) - M(2,3)*M(3,2) )
                dd(2) = M(1,2)*( M(2,3)*M(3,1) - M(2,1)*M(3,3) )
                dd(3) = M(1,3)*( M(2,1)*M(3,2) - M(2,2)*M(3,1) )
                dd(4) = M(2,1)*( M(3,2)*M(1,3) - M(3,3)*M(1,2) )
                dd(5) = M(2,2)*( M(3,3)*M(1,1) - M(3,1)*M(1,3) )
                dd(6) = M(2,3)*( M(3,1)*M(1,2) - M(3,2)*M(1,1) )
                dd(7) = M(3,1)*( M(1,2)*M(2,3) - M(1,3)*M(2,2) )
                dd(8) = M(3,2)*( M(1,3)*M(2,1) - M(1,1)*M(2,3) )
                dd(9) = M(3,3)*( M(1,1)*M(2,2) - M(1,2)*M(2,1) )
                d = (1.0d0/3.0d0) * sum(dd)
                return
            end function determinant3Mat
            
        
        subroutine initRotMats( R )
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      set up randomly oriented rotation matrices
            real(kind=real64),dimension(:,:,:),intent(inout)       ::      R
            integer                             ::      ii
              
            do ii = 1,size(R,dim=3)
                R(1:3,1,ii) = rndVecOnSphere() 
                do 
                    R(1:3,3,ii) = rndVecOnSphere() 
                    if (abs( dot_product( R(1:3,3,ii),R(1:3,1,ii) ) )<0.9 ) exit
                end do
                call completeBasis( R(1:3,3,ii),R(1:3,1,ii),R(1:3,2,ii) )
           
            end do
            !stop
            return            
            
        end subroutine initRotMats
            
        
        
            function rndVecOnSphere() result(x)
        !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
                real(kind=real64),dimension(3)  ::      x
                real(kind=real64)               ::      dd
                do 
                    call random_number(x)
                    x = 2*x - 1
                    dd = x(1)*x(1) + x(2)*x(2) + x(3)*x(3)
                    if (dd*(1-dd) > 0.0d0) exit
    
                end do 
                x = x/sqrt(dd)
                
                return
            end function rndVecOnSphere 
        
        pure function rotmat(thetax,thetay,thetaz) result(R)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
            real(kind=real64),intent(in)        ::      thetax,thetay,thetaz
            real(kind=real64),dimension(3,3)    ::      R
            real(kind=real64)       ::      cx,sx,cy,sy,cz,sz
            sx = sin(thetax)
            cx = cos(thetax)       
            sy = sin(thetay)
            cy = cos(thetay)       
            sz = sin(thetaz)
            cz = cos(thetaz)       
            R(1:3,1) = (/ cy*cz , cz*sx*sy - cx*sz , cx*cz*sy + sx*sz /)
            R(1:3,2) = (/ cy*sz , sx*sy*sz + cx*cz , cx*sy*sz - cz*sx /)
            R(1:3,3) = (/ -sy   , cy*sx            , cx*cy            /) 
            return
        end function rotmat
        
        
        subroutine jackForce( R, f ,nStatic)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ 
    !*      compute the force with respect to infinitessimal rotations about the x,y,z axes
            real(kind=real64),dimension(:,:,:),intent(in)       ::      R
            real(kind=real64),dimension(:,:),intent(out)        ::      f       !   (1:3,N)
            integer,intent(in)                                  ::      nStatic
            real(kind=real64),parameter                 ::      DELTA = 1.0d-6
            
            real(kind=real64),dimension(3,3),parameter  ::  RX = reshape( (/        &
                            1.0d0,      0.0d0,      0.0d0,                          &
                            0.0d0,      cos(DELTA), sin(DELTA),                     &
                            0.0d0,     -sin(DELTA), cos(DELTA)      /),(/3,3/) )
            
            real(kind=real64),dimension(3,3),parameter  ::  RY = reshape( (/        &
                            cos(DELTA), 0.0d0,      -sin(DELTA),                    &
                            0.0d0,      1.0d0,      0.0d0,                          &
                            sin(DELTA), 0.0d0,      cos(DELTA)      /),(/3,3/) )
                            
            real(kind=real64),dimension(3,3),parameter  ::  RZ = reshape( (/        &
                            cos(DELTA), sin(DELTA), 0.0d0,                          &
                           -sin(DELTA), cos(DELTA), 0.0d0,                          &
                            0.0d0,      0.0d0,      1.0d0           /),(/3,3/) )
                            
                            
            real(kind=real64),dimension(3,3),parameter  ::  RXT = transpose(RX)
            real(kind=real64),dimension(3,3),parameter  ::  RYT = transpose(RY)                        
            real(kind=real64),dimension(3,3),parameter  ::  RZT = transpose(RZ)
                            
            integer         ::      nn,ii,jj
            real(kind=real64),dimension(3,3)        ::      Ri
            nn = size(R,dim=3)
            
            f = 0.0d0
            do ii = nStatic+1,nn
                Ri = matmul( RXT,R(:,:,ii) )        
                f(1,ii) = f(1,ii) - jacksEnergy(Ri,ii,R)
                Ri = matmul( RX,R(:,:,ii) )        
                f(1,ii) = f(1,ii) + jacksEnergy(Ri,ii,R)
                
                Ri = matmul( RYT,R(:,:,ii) )        
                f(2,ii) = f(2,ii) - jacksEnergy(Ri,ii,R)
                Ri = matmul( RY,R(:,:,ii) )        
                f(2,ii) = f(2,ii) + jacksEnergy(Ri,ii,R)
                
                Ri = matmul( RZT,R(:,:,ii) )        
                f(3,ii) = f(3,ii) - jacksEnergy(Ri,ii,R)
                Ri = matmul( RZ,R(:,:,ii) )        
                f(3,ii) = f(3,ii) + jacksEnergy(Ri,ii,R)
            end do
            return
        end subroutine jackForce     
     
                
        pure function jacksEnergy(Ri,i,R) result(e)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !        return energy of jack i if has rotation Ri
            real(kind=real64),dimension(3,3),intent(in)         ::      Ri
            real(kind=real64),dimension(:,:,:),intent(in)       ::      R
            integer,intent(in)                                  ::      i
            real(kind=real64)                                   ::      e
            
            integer         ::      jj
            e = 0.0d0
            do jj = 1,size(R,dim=3)
                if (i == jj) cycle         
                e = e + jackEnergy(Ri,R(:,:,jj))
            end do
            return
        end function jacksEnergy
    
                
        pure function jackEnergy(R1,R2) result(e)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      consider two six-armed jacks, oriented with the rotation matrices R1,R2
    !*      
    !*              o
    !*              | o
    !*              |/
    !*        o-----+-----o
    !*             /|
    !*            o |
    !*              o
    !*      compute energy between the beads on the end of the arms
    !*
    
            real(kind=real64),dimension(3,3),intent(in)         ::      R1,R2
    
            real(kind=real64)                                   ::      e
             
            integer         ::      ii,jj
            real(kind=real64),dimension(3)       ::      x1,x2 
            real(kind=real64)        ::      dx,dy,dz,dd
            
            e = 0.0d0
            
            do ii = 1,3
                do jj = 1,3
                    dx = R2(1,ii) - R1(1,jj)
                    dy = R2(2,ii) - R1(2,jj)
                    dz = R2(3,ii) - R1(3,jj)
                    dd = dx*dx + dy*dy + dz*dz ; e = e + 2/(dd*dd*dd)        
                                                                                     
                    dx = R2(1,ii) + R1(1,jj)
                    dy = R2(2,ii) + R1(2,jj)
                    dz = R2(3,ii) + R1(3,jj)                                       
                    dd = dx*dx + dy*dy + dz*dz ; e = e + 2/(dd*dd*dd)          
                end do
            end do                          
            
            
            
            
!             do ii = 1,6
!                 select case(ii)
!                     case (1) ; x1(1:3) = R1(1:3,1)
!                     case (2) ; x1(1:3) = R1(1:3,2)
!                     case (3) ; x1(1:3) = R1(1:3,3)
!                     case (4) ; x1(1:3) = -R1(1:3,1)
!                     case (5) ; x1(1:3) = -R1(1:3,2)
!                     case (6) ; x1(1:3) = -R1(1:3,3)
!                 end select       
!                 do jj = 1,6
!                     select case(jj)
!                         case (1) ; x2(1:3) = R2(1:3,1)
!                         case (2) ; x2(1:3) = R2(1:3,2)
!                         case (3) ; x2(1:3) = R2(1:3,3)
!                         case (4) ; x2(1:3) = -R2(1:3,1)
!                         case (5) ; x2(1:3) = -R2(1:3,2)
!                         case (6) ; x2(1:3) = -R2(1:3,3)
!                     end select        
!                     e = e + efunc(x1,x2)
!                 end do
!             end do
            return
        end function jackEnergy
            
            
            
        pure real(kind=real64) function efunc(x1,x2) 
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      inverse 6th power repulsive function
            real(kind=real64),dimension(3),intent(in)       ::      x1,x2
            
            real(kind=real64)               ::      dx,dd
            real(kind=real64),parameter     ::      d2min = 1.0d-12
            
            dx = x2(1) - x1(1) ; dd = dx*dx
            dx = x2(2) - x1(2) ; dd = dd + dx*dx
            dx = x2(3) - x1(3) ; dd = dd + dx*dx
            
            dd = max(d2min,dd)
            efunc = 1/(dd*dd*dd)
            
            return
        end function efunc  
            
            

    end program optimiseAngles 
        
                