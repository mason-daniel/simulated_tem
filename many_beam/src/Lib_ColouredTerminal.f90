
    module Lib_ColouredTerminal
!---^^^^^^^^^^^^^^^^^^^^^^^^^^^
!*      Persuades fortran to output coloured text on the terminal.
!*      Note: can't guarantee that this will work with all compilers and all terminals.
!*      
!*      Daniel Mason
!*      (c) UKAEA
!*      May 2023
!*
!*
!*      Basic operation:
!*
!*          print *,"this is in "//colour(RED,"red")
!*      

        use iso_fortran_env
        implicit none
        private



        integer,public,parameter        ::      DEFAULT = 0
        integer,public,parameter        ::      DARK_GREY = 1
        integer,public,parameter        ::      PEACH = 2
        integer,public,parameter        ::      LIGHT_GREEN = 3
        integer,public,parameter        ::      LIGHT_YELLOW = 4
        integer,public,parameter        ::      LIGHT_BLUE = 5
        integer,public,parameter        ::      PINK = 6
        integer,public,parameter        ::      LIGHT_AQUA = 7
        integer,public,parameter        ::      PEARL_WHITE = 8
        integer,public,parameter        ::      BLACK = 9
        integer,public,parameter        ::      RED = 10
        integer,public,parameter        ::      GREEN = 11
        integer,public,parameter        ::      YELLOW = 12
        integer,public,parameter        ::      BLUE = 13
        integer,public,parameter        ::      PURPLE = 14
        integer,public,parameter        ::      AQUA = 15
        

        character(len=2),dimension(15),private,parameter                                &
                                        ::      COLOUR_CHARCODES = (/                   &
                                        "90","91","92","93","94","95","96","97",        &
                                        "30","31","32","33","34","35","36" /)    


        character(len=1),public,parameter   ::      NL = achar(10)
        character(len=1),public,parameter   ::      BS = achar(8)
                                        
                                        


        public          ::      colour
        public          ::      cutSpaces
        public          ::      compareResult
        
        
        interface       compareResult
            module procedure        compareResult_s
            module procedure        compareResult_i
            module procedure        compareResult_r32
            module procedure        compareResult_r64
            module procedure        compareResult_im
            module procedure        compareResult_r32m
            module procedure        compareResult_r64m
        end interface
        
    contains
!---^^^^^^^^

        pure function colour(c,text) result(ctext)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^  
    !*      produce a text string "text" decorated with terminal markup so it appears in the requested colour.
            integer,intent(in)                  ::      c
            character(len=*),intent(in)         ::      text
            character(len=len_trim(text)+9)    ::      ctext
            
            if ( (c-1)*(size(COLOUR_CHARCODES)-c) >= 0 ) then
                ctext = achar(27)//"["//COLOUR_CHARCODES(c)//"m"//trim(text)//achar(27)//"[0m"
            else
                ctext = trim(text)
            end if
            
            return
        end function colour
        
        
        pure function cutSpaces(text) result(ctext)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      remove leading and repeated spaces
            character(len=*),intent(in)         ::      text
            character(len=len(text))            ::      ctext   
            integer                 ::      ii,jj
            logical                 ::      lastSpace
            character(len=1)        ::      letter
            lastSpace = .true.
            jj = 0
            ctext = ""
            do ii = 1,len_trim(text) 
                letter = text(ii:ii)
                if (letter == " ") then
                    if (lastSpace) cycle
                    jj = jj + 1
                    ctext(jj:jj) = letter                   
                    lastSpace = .true.
                else
                    jj = jj + 1
                    ctext(jj:jj) = letter                    
                    lastSpace = .false.
                end if
            end do
            return
        end function cutSpaces
            
        
        
        pure logical function compareResult_s(groundTruth,test)  
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      returns true if test = groundTruth
            character(len=*),dimension(:),intent(in)         ::      groundTruth
            character(len=*),dimension(:),intent(in)         ::      test
            integer         ::      ii
            compareResult_s = .true.
            do ii = 1,size(groundTruth)
                compareResult_s = compareResult_s .and. ( trim(cutSpaces(groundTruth(ii))) == trim(cutSpaces(test(ii))) )
            end do
            return
        end function compareResult_s
        
        pure logical function compareResult_i(groundTruth,test)  
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      returns true if test = groundTruth
            integer,dimension(:),intent(in)         ::      groundTruth
            integer,dimension(:),intent(in)         ::      test
            integer         ::      ii
            compareResult_i = .true.
            do ii = 1,size(groundTruth)
                compareResult_i = compareResult_i .and. ( groundTruth(ii) == test(ii) )
            end do
            return
        end function compareResult_i
        
        pure logical function compareResult_r32(groundTruth,test,tol)  
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      returns true if test = groundTruth
            real(kind=real32),dimension(:),intent(in)       ::      groundTruth
            real(kind=real32),dimension(:),intent(in)       ::      test
            real(kind=real32),intent(in)                    ::      tol
            integer         ::      ii
            compareResult_r32 = .true.
            do ii = 1,size(groundTruth)
                compareResult_r32 = compareResult_r32 .and. ( abs( groundTruth(ii) - test(ii) ) <= abs(groundTruth(ii)*tol) )
            end do
            return
        end function compareResult_r32
        
        pure logical function compareResult_r64(groundTruth,test,tol)  
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      returns true if test = groundTruth
            real(kind=real64),dimension(:),intent(in)       ::      groundTruth
            real(kind=real64),dimension(:),intent(in)       ::      test
            real(kind=real64),intent(in)                    ::      tol
            integer         ::      ii
            compareResult_r64 = .true.
            do ii = 1,size(groundTruth)
                compareResult_r64 = compareResult_r64 .and. ( abs( groundTruth(ii) - test(ii) ) <= abs(groundTruth(ii)*tol) )
            end do
            return
        end function compareResult_r64                                  
        
        pure logical function compareResult_im(groundTruth,test)  
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      returns true if test = groundTruth
            integer,dimension(:,:),intent(in)         ::      groundTruth
            integer,dimension(:,:),intent(in)         ::      test
            integer         ::      ii,jj
            compareResult_im = .true.
            do ii = 1,size(groundTruth,dim=2)
                do jj = 1,size(groundTruth,dim=1)
                    compareResult_im = compareResult_im .and. ( groundTruth(jj,ii) == test(jj,ii) )
                end do
            end do
            return
        end function compareResult_im
         
         
        pure logical function compareResult_r32m(groundTruth,test,tol)  
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      returns true if test = groundTruth
            real(kind=real32),dimension(:,:),intent(in)       ::      groundTruth
            real(kind=real32),dimension(:,:),intent(in)       ::      test
            real(kind=real32),intent(in)                    ::      tol
            integer         ::      ii,jj
            compareResult_r32m = .true.
            do ii = 1,size(groundTruth,dim=2)
                do jj = 1,size(groundTruth,dim=1)
                    compareResult_r32m = compareResult_r32m .and. ( abs( groundTruth(jj,ii) - test(jj,ii) ) <= abs(groundTruth(jj,ii)*tol) )
                end do
            end do
            return
        end function compareResult_r32m
        
        
        pure logical function compareResult_r64m(groundTruth,test,tol)  
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      returns true if test = groundTruth
            real(kind=real64),dimension(:,:),intent(in)     ::      groundTruth
            real(kind=real64),dimension(:,:),intent(in)     ::      test
            real(kind=real64),intent(in)                    ::      tol
            integer         ::      ii,jj
            compareResult_r64m = .true.
            do ii = 1,size(groundTruth,dim=2)
                do jj = 1,size(groundTruth,dim=1)
                    compareResult_r64m = compareResult_r64m .and. ( abs( groundTruth(jj,ii) - test(jj,ii) ) <= abs(groundTruth(jj,ii)*tol) )
                end do
            end do
            return
        end function compareResult_r64m
            
    end module Lib_ColouredTerminal
    
    

           
        
        
