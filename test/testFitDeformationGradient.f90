    program testfitDeformationGradient
!---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
        use Lib_RandomSeed
        use fitDeformationGradient
        use Lib_RotationMatrices
        use Lib_Quaternions
        use Lib_DeformationGradients
        use iso_fortran_env
        implicit none
        
        character(len=256)      ::      filename = "Result/test0040.dat"
        
        real(kind=real64),dimension(3,39),parameter       ::      v0 = reshape(            &
                                                                (/  0.0d0,0.0d0,0.0d0,      &
                                                                    0.5d0,0.5d0,0.5d0,      &
                                                                    -.5d0,0.5d0,0.5d0,      &
                                                                    0.5d0,-.5d0,0.5d0,      &
                                                                    0.5d0,0.5d0,-.5d0,      &
                                                                    -.5d0,-.5d0,0.5d0,      &
                                                                    -.5d0,0.5d0,-.5d0,      &
                                                                    0.5d0,-.5d0,-.5d0,      &                                                                                                                                                                        
                                                                    -.5d0,-.5d0,-.5d0,      &
                                                                    1.0d0,0.0d0,0.0d0,      &
                                                                    -1.d0,0.0d0,0.0d0,      &
                                                                    0.0d0,1.0d0,0.0d0,      &
                                                                    0.0d0,-1.d0,0.0d0,      &
                                                                    0.0d0,0.0d0,1.0d0,      &
                                                                    0.0d0,0.0d0,-1.d0,      &
                                                                    1.5d0,0.5d0,0.5d0,      &
                                                                    1.5d0,-.5d0,0.5d0,      &                                                                    
                                                                    1.5d0,0.5d0,-.5d0,      &                                                                    
                                                                    1.5d0,-.5d0,-.5d0,      &                                                                    
                                                                   -1.5d0,0.5d0,0.5d0,      &
                                                                   -1.5d0,-.5d0,0.5d0,      &                                                                    
                                                                   -1.5d0,0.5d0,-.5d0,      &                                                                    
                                                                   -1.5d0,-.5d0,-.5d0,      &                                                                    
                                                                   0.5d0, 1.5d0,0.5d0,      &
                                                                   -.5d0, 1.5d0,0.5d0,      &                                                                    
                                                                   0.5d0, 1.5d0,-.5d0,      &                                                                    
                                                                   -.5d0, 1.5d0,-.5d0,      &                                                                    
                                                                   0.5d0,-1.5d0,0.5d0,      &
                                                                   -.5d0,-1.5d0,0.5d0,      &                                                                    
                                                                   0.5d0,-1.5d0,-.5d0,      &                                                                    
                                                                   -.5d0,-1.5d0,-.5d0,      &                                                                    
                                                                   0.5d0,0.5d0, 1.5d0,      &
                                                                   -.5d0,0.5d0, 1.5d0,      &                                                                    
                                                                   0.5d0,-.5d0, 1.5d0,      &                                                                    
                                                                   -.5d0,-.5d0, 1.5d0,      &                                                                    
                                                                   0.5d0,0.5d0,-1.5d0,      &
                                                                   -.5d0,0.5d0,-1.5d0,      &                                                                    
                                                                   0.5d0,-.5d0,-1.5d0,      &                                                                    
                                                                   -.5d0,-.5d0,-1.5d0   /),(/3,39/) )                                                                    
                                                                    
                                                                    
                                                                    
                                                                    
                                                                    
                                                                    
                                                                    
                                                                    
                                                                   
        real(kind=real64),dimension(3,3)    ::      RR,AA,T0,TT
        real(kind=real64),dimension(3,size(v0,dim=2))   ::      vv      
        real(kind=real64),dimension(3)      ::      delta                                                     
        type(FitDG)         ::      fdg
        integer             ::      ii
        integer             ::      steps
        integer,parameter   ::      NSTEPS = 1000 
        real(kind=real64)   ::      e1,e2 ,ee
        character(len=256)  ::      dummy
        real(kind=real64)   ::      maxstrain,maxtrans,maxdisp
        !type(Quaternion)    ::      qq
        !real(kind=real64),dimension(4)  ::      qq
        
        print *,"usage: testfitDeformationGradient.exe maxstrain  maxtrans  maxdisp"
        call get_command_argument( 1,dummy ) ; read(dummy,fmt=*) maxstrain
        call get_command_argument( 2,dummy ) ; read(dummy,fmt=*) maxtrans
        call get_command_argument( 3,dummy ) ; read(dummy,fmt=*) maxdisp
        
        
        
         call init_random_seed(12345)
        
       ! do ii = 1,100
       !     RR = rndRotMat()
       !     call report(RR)
       !     qq = RotMatToQuaternion(RR)
       !     write(*,fmt='(4f16.8)') qq
       !     
       !     TT = quaternionToRotMat(qq)
       !     call report(TT)
       !     
       !     e1 = RotMatDistance(RR,TT)
       !     print *,e1
       !     print *,""
       !     print *,""
       ! end do
       ! stop
        
        
         
        
        fdg =    FitDG_Ctor(v0,filename)
        
        call report(fdg)
        
        
        e1 = 0.0d0
        e2 = 0.0d0
        do steps = 1,NSTEPS
        !    print *,""
        !    print *,""
        !    print *,""
        !    print *,""
        
        !---    generate a random deformation gradient
            RR = rndRotMat() ; call tidyRotationMatrix(RR)
            call random_number(AA) ; AA = 2*(AA-1) * maxstrain ; AA = (AA + transpose(AA))/2             
            call StrainAndRotMatToDefGrad( AA,RR , T0 )
             
        !---    rotate and add a random offset    
            call random_number(delta) ; delta = (2*delta-1) * maxtrans
            do ii = 1,size(vv,dim=2)
                vv(:,ii) = rotateVector(T0,v0(:,ii)) + delta(:)
            end do            
            do ii = 1,size(vv,dim=2)
                call random_number(delta) ; delta = (2*delta-1) * maxdisp
                vv(:,ii) = vv(:,ii) + delta(:)
            end do
            
            
        !---    fit     
            
            call fit( fdg, vv, TT ,  0.5d0 , ee)
            e1 = e1 + ee
            
            ee = distance( T0,TT,QUATERNION_SYMMETRY_CUBIC )
            e2 = e2 + ee
            
         !   if (ee>1) then
         !       
         !       print *,""
         !       print *,"random rotation"
         !       call report(RR)
         !       print *,"random strain"
         !       write(*,fmt='(3f12.6)') AA(1,:)
         !       write(*,fmt='(3f12.6)') AA(2,:)
         !       write(*,fmt='(3f12.6)') AA(3,:)        
         !       print *,"random translation "
         !       print *,delta
         !       
         !       print *,"random deformation gradient"
         !       print *,T0(1,:)
         !       print *,T0(2,:)
         !       print *,T0(3,:)
         !       print *,""
!        !        do ii = 1,getnRotationMatrices(fdg)                
!        !            print *,"quat dist ",ii,RotMatDistance(getRotationMatrix(fdg,ii),RR)
!        !        end do
         !       print *,"output def grad"
         !       print *,TT(1,:)
         !       print *,TT(2,:)
         !       print *,TT(3,:)
         !       print *,"err ",ee
         !       
         !       print *,""
         !       open(unit=700,file="test.xyz",action="write")
         !       write (unit=700,fmt=*) "78"
         !       !write (unit=700,fmt=*) "Lattice=""2.0 0.0 0.0 0.0 2.0 0.0 0.0 0.0 2.0"" Properties=species:S:1:pos:R:3:conc:R:1"
         !       write (unit=700,fmt=*) "atom position_x position_y position_z charge"
         !       do ii = 1,39
         !           write (unit=700,fmt='(a2,4f12.5)') "I ",vv(:,ii),1.0d0
         !       end do
         !       do ii = 1,39
         !           write (unit=700,fmt='(a2,4f12.5)') "I ",matmul( TT(:,:),v0(:,ii) ),2.0d0
         !       end do
         !       close(unit=700)
         !       stop
         !   end if
            
        end do
        
        e1 = sqrt( (e1/NSTEPS)/size(vv,dim=2) )
        e2 = sqrt( (e2/NSTEPS)/3 )
        		
        print *,maxstrain,maxtrans,maxdisp,"pos err ",e1,"def grad err ",e2 
            
        
        
                                                                            
        call delete(fdg)
        print *,""
        print *,"done"
        print *,""                                                            
    
            
        
        
    end program testfitDeformationGradient
    