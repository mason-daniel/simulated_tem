
    module Lib_ColourScale
!---^^^^^^^^^^^^^^^^^^^^^^
!*      simple module for converting a value (0:1) into a colour scale 24-bit rgb
        use iso_fortran_env
        implicit none
        private
        
    !---        

        public      ::      getRGB,getRGB_double 
        public      ::      listAvailableColourScales      
        public      ::      getColourScale  
        public      ::      getColourScaleName
        
    !---        
        
        integer,public,parameter        ::      NCOLOURSCALES = 6
        integer,public,parameter        ::      COLOURSCALE_UNSET     = -1
        integer,public,parameter        ::      COLOURSCALE_GREYSCALE = 0
        integer,public,parameter        ::      COLOURSCALE_REDBLUE   = 1
        integer,public,parameter        ::      COLOURSCALE_BLUERED   = 2
        integer,public,parameter        ::      COLOURSCALE_VIRIDIS   = 3
        integer,public,parameter        ::      COLOURSCALE_MAGMA     = 4
        integer,public,parameter        ::      COLOURSCALE_BGR       = 5
        
        integer,private,parameter       ::      NKNOTS = 9
        
        integer,dimension(3,0:NKNOTS-1,0:NCOLOURSCALES-1),private,parameter   ::      COLOURSCALE_KNOTS = reshape( (/       &
                            000,000,000 , 031,031,031 , 063,063,063 , 095,095,095 , 127,127,127 , 159,159,159 , 191,191,191 , 223,223,223 , 255,255,255 ,           &
                            255,000,000 , 255,063,063 , 255,127,127 , 255,191,191 , 255,255,255 , 191,191,255 , 127,127,255 , 063,063,255 , 000,000,255 ,           &
                            000,000,255 , 063,063,255 , 127,127,255 , 191,191,255 , 255,255,255 , 255,191,191 , 255,127,127 , 255,063,063 , 255,000,000 ,           &
                            068,021,084 , 062,053,112 , 057,086,140 , 044,118,139 , 031,150,139 , 073,173,112 , 115,208,085 , 184,219,061 , 253,231,037 ,           &
                            000,000,000 , 042,011,062 , 085,023,125 , 134,044,123 , 184,055,121 , 218,096,109 , 252,138,098 , 252,195,144 , 253,253,191 ,           &
                            000,000,127 , 000,000,255 , 015,111,239 , 031,223,223 , 000,255,000 , 255,255,000 , 255,127,000 , 255,000,000 , 127,000,000             &
                             /) , (/3,NKNOTS,NCOLOURSCALES/) )
                            
        character(len=10),dimension(0:NCOLOURSCALES-1),public,parameter   ::      COLOURSCALE_NAME =  (/    &
                        "greyscale ",                                                                       &
                        "redblue   ",                                                                       &
                        "bluered   ",                                                                       &
                        "viridis   ",                                                                       &
                        "magma     ",                                                                       &
                        "bgr       "     /)
                            
                            
       
    !---        
        
        interface   getRGB
            module procedure    getRGB0
        end interface                            
        
    !---        
    
    contains
!---^^^^^^^^

        subroutine listAvailableColourScales()
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
            integer     ::      ii
            print *,"Lib_ColourScale::listAvailableColourScales() info"
            do ii = 0,NCOLOURSCALES-1
                write(*,fmt='(a)') """"//trim(COLOURSCALE_NAME(ii))//""""
            end do
            return
        end subroutine listAvailableColourScales
        

        function getColourScale( cs ) result(colourscale)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      convert character string to integer
            character(len=*),intent(in)     ::      cs
            integer                         ::     colourscale
            integer         ::      ii
            colourscale = COLOURSCALE_UNSET
            do ii = 0,NCOLOURSCALES-1
                if (trim(cs)==trim(COLOURSCALE_NAME(ii))) then
                    colourscale = ii
                    return
                end if
            end do
            
            return
        end function getColourScale
        

        function getColourScaleName( cs ) result(colourscale)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      convert integer to character string 
            integer,intent(in)    ::      cs
            character(len=10)     ::      colourscale
            colourscale = "unset"
            if ( cs*(NCOLOURSCALES-1 - cs) >= 0) colourscale = COLOURSCALE_NAME(cs)
            
            return
        end function getColourScaleName        

        function getRGB0(colourscale,x) result( rgb )
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
            integer,intent(in)                      ::      colourscale
            real(kind=real64),intent(in)            ::      x
            integer                                 ::      rgb
            
            integer         ::      rr,gg,bb
            integer         ::      ii
            real(kind=real64)       ::      xx
            
            ii = floor( x*NKNOTS )
            if (ii<0) then
                rr = COLOURSCALE_KNOTS(1,0,colourscale)
                gg = COLOURSCALE_KNOTS(2,0,colourscale)
                bb = COLOURSCALE_KNOTS(3,0,colourscale)
            else if (ii>=NKNOTS-1) then                               
                rr = COLOURSCALE_KNOTS(1,NKNOTS-1,colourscale)
                gg = COLOURSCALE_KNOTS(2,NKNOTS-1,colourscale)
                bb = COLOURSCALE_KNOTS(3,NKNOTS-1,colourscale)
            else
                xx = x*NKNOTS - ii
                rr = int( COLOURSCALE_KNOTS(1,ii,colourscale)*(1-xx) + COLOURSCALE_KNOTS(1,ii+1,colourscale)*xx )
                gg = int( COLOURSCALE_KNOTS(2,ii,colourscale)*(1-xx) + COLOURSCALE_KNOTS(2,ii+1,colourscale)*xx )
                bb = int( COLOURSCALE_KNOTS(3,ii,colourscale)*(1-xx) + COLOURSCALE_KNOTS(3,ii+1,colourscale)*xx )
            end if
            
            
            
            rgb = ishft(rr,16) + ishft(gg,8) + bb
            !print *," getRGB0 ",colourscale,x,rr,gg,bb,rgb
            return
        end function getRGB0
                
                
        function getRGB_double(colourscale,x) result( rgb )
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
            integer,intent(in)                      ::      colourscale
            real(kind=real64),intent(in)            ::      x
            real(kind=real64),dimension(3)          ::      rgb
            
             
            integer                 ::      ii
            real(kind=real64)       ::      xx
            
            ii = floor( x*NKNOTS )
            if (ii<0) then
                rgb(1:3) = COLOURSCALE_KNOTS(1:3,0,colourscale)
            else if (ii>=NKNOTS-1) then                               
                rgb(1:3) = COLOURSCALE_KNOTS(1:3,NKNOTS-1,colourscale)
            else
                xx = x*NKNOTS - ii
                rgb(1:3) =  COLOURSCALE_KNOTS(1:3,ii,colourscale)*(1-xx) + COLOURSCALE_KNOTS(1:3,ii+1,colourscale)*xx 
            end if
            rgb(1:3) = rgb(1:3) / 256.0d0
            
            return
        end function getRGB_double
                
                
                
                
        


    end module Lib_ColourScale
        