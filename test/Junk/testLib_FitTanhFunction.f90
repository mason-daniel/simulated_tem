    
    program testLib_FitTanhFunction
!---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
!*  
!*      simple program to test functioning of Lib_DeformationGradients
!*      
!*          successful result            
!*      
!*       solution b,b~,b-b~
!*           -0.78296868     -0.77953408     -0.00343460
!*           -0.74386581     -0.73980624     -0.00405957
!*           -0.87656189     -0.87705720      0.00049531
!*       solution c,c~,c-c~
!*            0.31963714      0.32280668     -0.00316953
!*       solution |b|,|b~|
!*            1.39094826      1.38716113
!*       solution b.b~
!*            0.99999633
!*       solution with target vector length
!*       solution |b|,|b~|
!*            1.39094826      1.09659421
!*       solution b.b~
!*            0.99997361
!*       
!*       done
!*       

!*      







       use iso_fortran_env
       use Lib_FitTanhFunction
       use Lib_ColouredTerminal
       implicit none
       
       integer,parameter       ::      d = 3
       integer,parameter       ::      n = 20
       
       real(kind=real64),dimension(d)      ::      b = (/ -0.7829686807d0, -0.7438658125d0, -0.8765618915d0 /)
       real(kind=real64)                   ::      c =  0.3196371424d0
       
       real(kind=real64),dimension(d)      ::      b_out       
       real(kind=real64)                   ::      c_out
       
       real(kind=real64),dimension(d,n)    ::      x = reshape( (/    0.8235128657d0, 0.1950620945d0, 0.7831443263d0,               &
                                                                     -0.7907728630d0,-0.0215808419d0,-0.6239430759d0,               &
                                                                     -0.3558610493d0, 0.9162936023d0, 0.4633187763d0,               &
                                                                      0.7979425051d0,-0.2932464982d0,-0.1787896825d0,               &
                                                                     -0.9739186154d0, 0.8377647627d0,-0.7767146954d0,               &
                                                                      0.0153054681d0,-0.1576081955d0, 0.1100798451d0,               &
                                                                     -0.0276193427d0,-0.0092659799d0, 0.1355033687d0,               &
                                                                      0.5104462035d0,-0.5541775385d0, 0.7692343908d0,               &
                                                                      0.3266051025d0, 0.3228215467d0, 0.9393348078d0,               &
                                                                      0.3884169796d0,-0.3738354128d0,-0.6930719825d0,               &
                                                                     -0.5710550598d0,-0.5457554733d0,-0.3386480622d0,               &
                                                                      0.0079896445d0,-0.8078890519d0,-0.7472727827d0,               &
                                                                     -0.6489579109d0, 0.2043826115d0, 0.3522026367d0,               &
                                                                      0.3247775299d0,-0.4970638046d0,-0.2634439352d0,               &
                                                                      0.8483528617d0,-0.3702038542d0,-0.2081967089d0,               &
                                                                      0.9932846647d0, 0.4338593807d0,-0.3928083854d0,               &
                                                                      0.9483526196d0, 0.1432590122d0,-0.5343774553d0,               &
                                                                     -0.4974771350d0,-0.5091789089d0, 0.6982024250d0,               &
                                                                      0.5019900072d0,-0.2803025843d0,-0.7373132798d0,               &
                                                                      0.5792579242d0,-0.7753215218d0, 0.4346848297d0    /),(/d,n/) )
       
       
       real(kind=real64),dimension(n)      ::      f = (/ -0.8222197440d0 ,     &
                                                           0.9071461090d0 ,     &
                                                          -0.4490032259d0 ,     &
                                                           0.0780673817d0 ,     &
                                                           0.8180390509d0 ,     &
                                                           0.3127426478d0 ,     &
                                                           0.2346310912d0 ,     &
                                                          -0.3224951584d0 ,     &
                                                          -0.7628814952d0 ,     &
                                                           0.7206844215d0 ,     &
                                                           0.8938379761d0 ,     &
                                                           0.9187024788d0 ,     &
                                                           0.3561895199d0 ,     &
                                                           0.5817015842d0 ,     &
                                                           0.1212550724d0 ,     &
                                                          -0.4065057985d0 ,     &
                                                          -0.0554514414d0 ,     &
                                                           0.4367851406d0 ,     &
                                                           0.6512540780d0 ,     &
                                                           0.0625137288d0 /)
       
       
        character(len=256),dimension(15)            ::      output
        character(len=*),dimension(15),parameter   ::      output0 = (/ "solution b,b~,b-b~                             ",         &
                                                                        "    -0.78297        -0.77953        -0.00343   ",         &
                                                                        "    -0.74387        -0.73981        -0.00406   ",         &
                                                                        "    -0.87656        -0.87706         0.00050   ",         &
                                                                        "solution c,c~,c-c~                             ",         &
                                                                        "     0.31964         0.32281        -0.00317   ",         &
                                                                        "solution |b|,|b~|                              ",         &                                                            
                                                                        "     1.39095         1.38716                   ",         &
                                                                        "solution b.b~                                  ",         &
                                                                        "     1.00000                                   ",         &
                                                                        "solution with target vector length             ",         &
                                                                        "solution |b|,|b~|                              ",         &
                                                                        "     1.39095         1.09659                   ",         &
                                                                        "solution b.b~                                  ",         &                                                                           
                                                                        "     0.99997                                   "    /)
        
         
        integer                     ::      ii                                                             
        logical                     ::      ok                                                             
                                                      
               
                                                                  
                                                           
                                                           
                                                           
                                                           
                                                           
                                                           
                                                           
                                                           
                                                           
   
!   !---     generate test
!   !---    determine arbitrary vector and offset
!       call random_number( b ) ; b = 2*( 2*b - 1 )
!       call random_number( c ) ; c = 2*c - 1
!       write(*,fmt='(a,100f16.10)') "b ",b
!       write(*,fmt='(a,100f16.10)') "c ",c
!       
!   !---    generate some points in 3d 
!       call random_number( x ) ; x = 2*x - 1
!                     
!   !---    set data to approximate tanh function
!        call random_number(f)
!        f = (2*f-1) * 0.01d0
!       do ii = 1,n
!           f(ii) = f(ii) + tanh( dot_product(b,x(:,ii)) + c )
!           write(*,fmt='(i6,100f16.10)') ii,x(:,ii),f(ii)
!       end do
!       
       
       
       call fitTanhFunction( x,f, b_out,c_out )
       

       output(1) = "solution b,b~,b-b~"
       do ii = 1,d
           write (output(1+ii),fmt='(100f16.5)') b(ii),b_out(ii),b(ii)-b_out(ii)
       end do
       output(5) = "solution c,c~,c-c~"
       write (output(6),fmt='(100f16.5)') c,c_out,c-c_out
       output(7) = "solution |b|,|b~|"
       write (output(8),fmt='(100f16.5)') norm2( b ),norm2(b_out )
       
       output(9) = "solution b.b~"
       write (output(10),fmt='(100f16.5)') dot_product( b,b_out )/( norm2(b)*norm2(b_out) )
       
       call fitTanhFunction( x,f, b_out,c_out,rc=sqrt(d*1.0d0) )
       

       output(11) ="solution with target vector length"
       output(12) ="solution |b|,|b~|"
       write (output(13),fmt='(100f16.5)') norm2( b ),norm2(b_out )
       
       output(14) ="solution b.b~"
       write (output(15),fmt='(100f16.5)') dot_product( b,b_out )/( norm2(b)*norm2(b_out) )
       
       
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
            
       
   end program testLib_FitTanhFunction