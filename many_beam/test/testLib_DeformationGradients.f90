     
    
    program testLib_DeformationGradients
!---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
!*  
!*      simple program to test functioning of Lib_DeformationGradients
!*      
!*          successful result        
!*          
!*           rotation matrix
!*               -0.810008   -0.558865   -0.177643
!*                0.298518   -0.132221   -0.945201
!*                0.504751   -0.818650    0.273932
!*           strain matrix
!*               -0.093257   -0.057282    0.040720
!*               -0.057282   -0.056864    0.043367
!*                0.040720    0.043367   -0.078460
!*           deformation gradient
!*               -0.731015   -0.532508   -0.095779
!*                0.349832   -0.128191   -0.869397
!*                0.445111   -0.782909    0.204215
!*           backwards
!*           rotation matrix
!*               -0.810008   -0.558865   -0.177643
!*                0.298518   -0.132221   -0.945201
!*                0.504751   -0.818650    0.273932
!*           strain matrix
!*               -0.093257   -0.057282    0.040720
!*               -0.057282   -0.056864    0.043367
!*                0.040720    0.043367   -0.078460
!*          strain invariants         -0.22858171      0.01026152     -0.00009126
!*          
!*           done
!*          
!*              

        use iso_fortran_env
        use Lib_RotationMatrices
        use Lib_ColouredTerminal
        use Lib_DeformationGradients
        implicit none
         
        
        
        
        real(kind=real64),dimension(3,3)   ::   TT
        
        real(kind=real64),dimension(3,3)   ::   RR  = reshape( (/-0.810007854270d0, 0.298518395941d0,0.504751466869d0,-0.558865047959d0,-0.132220996548d0,-0.818649782411d0,-0.177643277924d0,-0.945200706388d0, 0.273931543367d0 /),(/3,3/) )
        real(kind=real64),dimension(3,3)   ::   eps = reshape( (/ -0.093257239753d0,-0.057282484545d0,0.040719920472d0,-0.057282484545d0,-0.056864344007d0, 0.043366567941d0, 0.040719920472d0, 0.043366567941d0,-0.078460122236d0 /),(/3,3/) )

        real(kind=real64)           ::      i1,i2,i3
        
        
        character(len=256),dimension(4) ::      output
        character(len=*),dimension(4),parameter   ::      output0 = (/  "    -0.093257   -0.057282    0.040720                                ",         &
                                                                        "    -0.057282   -0.056864    0.043367                                ",         &
                                                                        "     0.040720    0.043367   -0.078460                                ",         &
                                                                        "strain invariants         -0.22858171      0.01026152     -0.00009126"          &
                                                                        /)      
        integer                     ::      ii                                                             
        logical                     ::      ok                                                             
                                                      
        
        
    !---    generate test
       ! RR = rndRotMat()
       ! print *,"random rotation matrix used "
       ! write(*,fmt='(9f20.12)') RR 
       ! call random_number(eps)
       ! eps = (eps*2 - 1) * 0.1d0  
       ! eps = ( eps + transpose(eps) )/2
       ! write(*,fmt='(9f20.12)') eps 
    !---        
        
        print *,"rotation matrix"
        call opMat(RR)
        
        print *,"strain matrix"
        call opMat(eps)
        
        call StrainAndRotMatToDefGrad( eps,RR , TT )        
        print *,"deformation gradient"
        call opMat(TT)
        
        eps = 0 ; RR = 0
        call DefGradToStrainAndRotMat(TT, eps,RR)
        print *,"backwards"
        
        print *,"rotation matrix"
        call opMat(RR)
        
        print *,"strain matrix"
        call opMat(eps)
        
        call strainInvariants(eps,i1,i2,i3)
        write(*,fmt='(a,3f16.8)') "strain invariants    ",i1,i2,i3
        
    !---    check result
    !    write(*,fmt='(a,3f16.8)') "von Mises = sqrt(I1^2 - 3 I2) = ",sqrt(i1*i1-3*i2)
    !    write(*,fmt='(a,3f16.8)') "cf Tr,von Mises,det  ",eps(1,1)+eps(2,2)+eps(3,3),vonMises3Mat(eps),determinant3Mat(eps)
    
        write( output(1),fmt='(3f12.6 )') eps(1,1:3)
        write( output(2),fmt='(3f12.6 )') eps(2,1:3)
        write( output(3),fmt='(3f12.6 )') eps(3,1:3)
        write( output(4),fmt='(a,3f16.8)') "strain invariants    ",i1,i2,i3
            
       
    !---    here is the simple test: does the output look like the stored output?
        ok = .true.
        do ii = 1,size(output)            
            if ( trim(cutSpaces(output(ii))) == trim(cutSpaces(output0(ii))) ) then
                write (*,fmt='(a)') trim(cutSpaces(output(ii)))
            else
                ok = .false.
                write (*,fmt='(a)') trim(cutSpaces(output0(ii)))//"    "//colour(RED,trim(cutSpaces(output(ii))))
            end if
        end do   
        
       
    !---    output the result "PASS" or "FAIL"    
        if (ok) then
            print *,colour(LIGHT_GREEN,"PASS")
        else
            print *,colour(RED,"FAIL")
        end if
        
        print *,""
        print *,"done"
        print *,""
        
        
      contains
  !---^^^^^^^^
  
!          pure function determinant3Mat(M) result(d)
!      !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
!     !*      returns the determinant of M
!             real(kind=real64),dimension(3,3),intent(in)      ::      M
!             real(kind=real64)                                ::      d
!             real(kind=real64),dimension(9)       ::      dd
!             dd(1) = M(1,1)*( M(2,2)*M(3,3) - M(2,3)*M(3,2) )
!             dd(2) = M(1,2)*( M(2,3)*M(3,1) - M(2,1)*M(3,3) )
!             dd(3) = M(1,3)*( M(2,1)*M(3,2) - M(2,2)*M(3,1) )
!             dd(4) = M(2,1)*( M(3,2)*M(1,3) - M(3,3)*M(1,2) )
!             dd(5) = M(2,2)*( M(3,3)*M(1,1) - M(3,1)*M(1,3) )
!             dd(6) = M(2,3)*( M(3,1)*M(1,2) - M(3,2)*M(1,1) )
!             dd(7) = M(3,1)*( M(1,2)*M(2,3) - M(1,3)*M(2,2) )
!             dd(8) = M(3,2)*( M(1,3)*M(2,1) - M(1,1)*M(2,3) )
!             dd(9) = M(3,3)*( M(1,1)*M(2,2) - M(1,2)*M(2,1) )
!             d = (1.0d0/3.0d0) * sum(dd)
!             return
!         end function determinant3Mat
!         
!         
!         pure function vonMises3Mat(M) result(d)
!     !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
!     !*      returns the von mises invariant of M
!             real(kind=real64),dimension(3,3),intent(in)      ::      M
!             real(kind=real64)                                ::      d
! 
!             d =     M(2,3)*M(3,2)
!             d = d + M(3,1)*M(1,3)
!             d = d + M(1,2)*M(2,1)
!             
!                         
!             d = 6*d + (M(1,1) - M(2,2))**2
!             d =   d + (M(2,2) - M(3,3))**2
!             d =   d + (M(3,3) - M(1,1))**2
!              
!             d = sqrt(max(0.0d0,d/2))                !   max to avoid -0.0 rounding error
!             
!             
!             return
!         end function vonMises3Mat
        
        
        subroutine opMat(this,u,o)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      simple output. Defaults to unit 6 = screen.
    !*      optional argument o determines left hand margin
            real(kind=real64),dimension(3,3),intent(in)    ::      this
            integer,intent(in),optional     ::      u,o
            integer     ::      uu,oo
            uu = 6 ; if (present(u)) uu = u
            oo = 0 ; if (present(o)) oo = o
            write(unit=uu,fmt='(a,3f12.6)') repeat(" ",oo+2),this(1,1:3)
            write(unit=uu,fmt='(a,3f12.6)') repeat(" ",oo+2),this(2,1:3)
            write(unit=uu,fmt='(a,3f12.6)') repeat(" ",oo+2),this(3,1:3)
            return
        end subroutine opMat
    
    
        
    end program testLib_DeformationGradients
