
    module Lib_Filenames
!---^^^^^^^^^^^^^^^^^^^^ 
!*      
!*      short module to handle file numbering and suffixes
!*  
!*      Daniel Mason, UKAEA
!*      April 2022
!*
!*
    
        use NBAX_StringTokenizers
        use iso_fortran_env
        implicit none
        private
        
    !---
        
    
        public          ::      getSuffix
        public          ::      removeSuffix
        public          ::      numberFile
        
        public          ::      versionNumber
        public          ::      getVersionNumber
        public          ::      compareVersionNumber
        
    !--
        
        interface   numberFile
            module procedure    numberFile1 
            module procedure    numberFile2
        end interface
        
    contains
!---^^^^^^^^
        
        function getSuffix( file_in ) result( suffix )
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ 
    !*      return filename suffix - eg "directory/foo.bar.txt" becomes "txt"
    
            character(len=*),intent(in)     ::      file_in
            character(len=len(file_in))     ::      suffix
            
            integer     ::      ii
            ii = index( file_in,".",back=.true. )
            
            suffix = ""
            if (ii>index( file_in,"/",back=.true. )) then                   !   last "." after last "/"
                suffix = file_in(ii+1:)          
            end if
            return
        end function getSuffix
        
        
        function removeSuffix( file_in ) result( file_out )
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      return filename with suffix removed - eg "directory/foo.bar.txt" becomes "directory/foo.bar"
    
            character(len=*),intent(in)     ::      file_in
            character(len=len(file_in))     ::      file_out
            
            integer     ::      ii
            file_out = ""
            ii = index( file_in,".",back=.true. )
            if (ii>0) then  
                file_out(1:ii-1) = file_in(1:ii-1)        
            else
                file_out = file_in
            end if
            return
        end function removeSuffix
        
        function numberFile1( file_prefix,n ) result( file_out )
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      return filename with number added
    !*      eg "directory/foo.bar" + [n=123] + "png" returns "directory/foo.bar.00123.png"
    !*      if suffix not supplied then it is automatically found
    !*      eg "directory/foo.bar" + [n=123]         returns "directory/foo.00123.bar"
    
            character(len=*),intent(in)             ::      file_prefix
            integer,intent(in)                      ::      n
            
            character(len=len(file_prefix)+7)       ::      file_out
            character(len=5)                    ::      aaaa
            
            integer     ::      nn
            nn = mod( n+900000,100000 )    !   ensure a number between 0:99999
            write (aaaa,fmt='(i5)') nn 
            aaaa = adjustl(aaaa)
            aaaa = repeat("0",5-len_trim(aaaa))//trim(aaaa)
            
            file_out = getSuffix( file_prefix )
            if (len_trim(file_out)>0) then
                file_out = trim(removeSuffix(file_prefix))//"."//aaaa//"."//trim(file_out)
            else
                file_out = trim(file_prefix)//"."//aaaa
            end if
            return
           
        end function numberFile1
        
  
        function numberFile2( file_prefix,n,file_suffix ) result( file_out )
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      return filename with number added
    !*      eg "directory/foo.bar" + [n=123] + "png" returns "directory/foo.bar.00123.png"
    !*      if suffix not supplied then it is automatically found
    !*      eg "directory/foo.bar" + [n=123]         returns "directory/foo.00123.bar"
    
            character(len=*),intent(in)             ::      file_prefix
            integer,intent(in)                      ::      n
            character(len=*),intent(in)             ::      file_suffix
            
            character(len=len(file_prefix)+len(file_suffix)+7)       ::      file_out
            character(len=5)                    ::      aaaa
            
            integer     ::      nn
            nn = mod( n,100000 )    !   ensure a number between 0:99999
            write (aaaa,fmt='(i5)') nn 
            aaaa = adjustl(aaaa)
            aaaa = repeat("0",5-len_trim(aaaa))//trim(aaaa)
            
            if (len_trim(file_suffix)>0) then
                file_out = trim(file_prefix)//"."//aaaa//"."//trim(file_suffix)
            else
                file_out = trim(file_prefix)//"."//aaaa
            end if
            
            return
           
        end function numberFile2
        
        pure function versionNumber(major,minor,patch) result(txt)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      write txt = "major"."minor"."patch"
    !*      assuming each value <= 999
            integer,intent(in)                  ::      major,minor,patch    
            character(len=11)                   ::      txt
            character(len=3)                    ::      val
            write(val,fmt='(i3)') major
            txt = trim(adjustl(val)) 
            write(val,fmt='(i3)') minor
            txt = trim(txt)//"."//trim(adjustl(val))
            write(val,fmt='(i3)') patch
            txt = trim(txt)//"."//trim(adjustl(val))
            return
        end function versionNumber
            
        subroutine getVersionNumber(txt,major,minor,patch)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      attempt to break up txt = "1.0.13" into major 1 minor 0 and patch 13 
    !*      defaults to 0 0 0 if there's an error
            character(len=*),intent(in)         ::      txt
            integer,intent(out)                 ::      major,minor,patch
            type(StringTokenizer)       ::      stok
            character(len=32)           ::      tok
            stok = StringTokenizer_ctor(txt," .,;"//TAB_CHARACTER//CR_CHARACTER//NULL_CHARACTER)
            major = 0
            minor = 0
            patch = 0
            call nextToken( stok,tok )
            call parse( tok,major )
            call nextToken( stok,tok )
            call parse( tok,minor )
            call nextToken( stok,tok )
            call parse( tok,patch )
            return
        end subroutine getVersionNumber    
            
        logical function compareVersionNumber(v1,v2)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      return true if version v1 is >= version v2
            character(len=*),intent(in)         ::      v1,v2
             
            integer                     ::      major1,minor1,patch1
            integer                     ::      major2,minor2,patch2
            
            call getVersionNumber(v1,major1,minor1,patch1)
            call getVersionNumber(v2,major2,minor2,patch2)
            compareVersionNumber = major1>=major2
            if (compareVersionNumber) compareVersionNumber = minor1>=minor2
            if (compareVersionNumber) compareVersionNumber = patch1>=patch2
          
           !print *,"Lib_Filenames::compareVersionNumber info - """//trim(v1)//""" = ",major1,minor1,patch1
           !print *,"                                           """//trim(v2)//""" = ",major2,minor2,patch2
           !print *,"                                           result ",compareVersionNumber
            
            return
        end function compareVersionNumber    
        
    end module Lib_Filenames
    
    