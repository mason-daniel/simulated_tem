
    module Lib_quicksort
!---^^^^^^^^^^^^^^^^^^^^
!*      simple implementation of the quicksort algorithm
!*          call quicksort(array [,back])
!*      shuffles the elements of array into size order, from smallest to largest ( or if back then largest to smallest )
!*          call quicksort(array,indx [,back])
!*      shuffles the indexing such that array(indx(i)) >= array(indx(i-1)), doesn't touch array


        use iso_fortran_env
        implicit none
        private
        
    !---        
        
        public      ::      quickSort
        
    !---        
    
        interface       quickSort
            module procedure    quickSort_i32       
            module procedure    quickSort_i64 
            module procedure    quickSort_r32       
            module procedure    quickSort_r64 
            module procedure    quickSort_indx_i32       
            module procedure    quickSort_indx_i64       
            module procedure    quickSort_indx_r32       
            module procedure    quickSort_indx_r64   
        end interface
        
    contains
!---^^^^^^^^                
    

        subroutine quickSort_i32( array,back )
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^       
            integer(kind=int32),dimension(:),intent(inout)          ::      array
            logical,intent(in),optional                             ::      back
            logical     ::      bb            
            integer     ::      nn
            integer(kind=int32) ::  aa
            nn = size(array)
            bb = .false. ; if (present(back)) bb = back
            if (nn<=1) then
                return
            else if (nn==2) then    
                aa = array(1)
                if (bb) then
                    if (array(2)>aa) then
                        array(1) = array(2) 
                        array(2) = aa
                    end if
                else
                    if (aa>array(2)) then
                        array(1) = array(2) 
                        array(2) = aa
                    end if                
                end if
            else
                call quickSort_i32a( array,bb,1,nn )
            end if
            return
        end subroutine quickSort_i32
                
        subroutine quickSort_i64( array,back )
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^       
            integer(kind=int64),dimension(:),intent(inout)          ::      array
            logical,intent(in),optional                             ::      back
            logical     ::      bb
            integer     ::      nn
            integer(kind=int64) ::  aa
            nn = size(array)
            bb = .false. ; if (present(back)) bb = back
            if (nn<=1) then
                return
            else if (nn==2) then    
                aa = array(1)
                if (bb) then
                    if (array(2)>aa) then
                        array(1) = array(2) 
                        array(2) = aa
                    end if
                else
                    if (aa>array(2)) then
                        array(1) = array(2) 
                        array(2) = aa
                    end if                
                end if
            else
                call quickSort_i64a( array,bb,1,nn )
            end if
            return
        end subroutine quickSort_i64
        
        subroutine quickSort_r32( array,back )
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^       
            real(kind=real32),dimension(:),intent(inout)            ::      array
            logical,intent(in),optional                             ::      back
            logical     ::      bb
            integer     ::      nn
            real(kind=real32) ::  aa
            nn = size(array)            
            bb = .false. ; if (present(back)) bb = back
            if (nn<=1) then
                return
            else if (nn==2) then    
                aa = array(1)
                if (bb) then
                    if (array(2)>aa) then
                        array(1) = array(2) 
                        array(2) = aa
                    end if
                else
                    if (aa>array(2)) then
                        array(1) = array(2) 
                        array(2) = aa
                    end if                
                end if
            else
                call quickSort_r32a( array,bb,1,nn )
            end if
            return
        end subroutine quickSort_r32
        
        subroutine quickSort_r64( array,back )
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^       
            real(kind=real64),dimension(:),intent(inout)            ::      array
            logical,intent(in),optional                             ::      back
            logical     ::      bb
            integer     ::      nn
            real(kind=real64) ::  aa
            nn = size(array)            
            bb = .false. ; if (present(back)) bb = back
            if (nn<=1) then
                return
            else if (nn==2) then    
                aa = array(1)
                if (bb) then
                    if (array(2)>aa) then
                        array(1) = array(2) 
                        array(2) = aa
                    end if
                else
                    if (aa>array(2)) then
                        array(1) = array(2) 
                        array(2) = aa
                    end if                
                end if
            else
                call quickSort_r64a( array,bb,1,nn )
            end if
            return
        end subroutine quickSort_r64     
           
        
        subroutine quickSort_indx_i32( array,indx,back )
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^  
            integer(kind=int32),dimension(:),intent(in)             ::      array
            integer,dimension(:),intent(inout)                      ::      indx            
            logical,intent(in),optional                             ::      back
            logical             ::      bb
            integer             ::      nn,ii
            nn = size(array)
            
            bb = .false. ; if (present(back)) bb = back
            if (nn<=1) then
                return
            else if (nn==2) then    
                ii = indx(1)
                if (bb) then
                    if (array(2)>array(1)) then
                        indx(1) = indx(2)
                        indx(2) = ii
                    end if
                else
                    if (array(1)>array(2)) then
                        indx(1) = indx(2) 
                        indx(2) = ii     
                    end if                
                end if
            else
                call quickSort_indx_i32a( array,indx,bb,1,nn )
            end if
                
            return
        end subroutine quickSort_indx_i32
        
        subroutine quickSort_indx_i64( array,indx,back )
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^  
            integer(kind=int64),dimension(:),intent(in)             ::      array
            integer,dimension(:),intent(inout)                      ::      indx            
            logical,intent(in),optional                             ::      back
            logical             ::      bb
            integer             ::      nn,ii
            nn = size(array)
            
            bb = .false. ; if (present(back)) bb = back
            if (nn<=1) then
                return
            else if (nn==2) then    
                ii = indx(1)
                if (bb) then
                    if (array(2)>array(1)) then
                        indx(1) = indx(2)
                        indx(2) = ii
                    end if
                else
                    if (array(1)>array(2)) then
                        indx(1) = indx(2) 
                        indx(2) = ii     
                    end if                
                end if
            else
                call quickSort_indx_i64a( array,indx,bb,1,nn )
            end if
            return
        end subroutine quickSort_indx_i64        
        
        subroutine quickSort_indx_r32( array,indx,back )
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^  
            real(kind=real32),dimension(:),intent(in)               ::      array
            integer,dimension(:),intent(inout)                      ::      indx            
            logical,intent(in),optional                             ::      back
            logical             ::      bb
            integer             ::      nn,ii
            nn = size(array)
            bb = .false. ; if (present(back)) bb = back
            if (nn<=1) then
                return
            else if (nn==2) then    
                ii = indx(1)
                if (bb) then
                    if (array(2)>array(1)) then
                        indx(1) = indx(2)
                        indx(2) = ii
                    end if
                else
                    if (array(1)>array(2)) then
                        indx(1) = indx(2) 
                        indx(2) = ii     
                    end if                
                end if
            else
                call quickSort_indx_r32a( array,indx,bb,1,nn )
            end if
            return
        end subroutine quickSort_indx_r32
        
        
        subroutine quickSort_indx_r64( array,indx,back )
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^  
            real(kind=real64),dimension(:),intent(in)               ::      array
            integer,dimension(:),intent(inout)                      ::      indx            
            logical,intent(in),optional                             ::      back
            logical             ::      bb
            integer             ::      nn,ii
            nn = size(array)
            
            bb = .false. ; if (present(back)) bb = back
            if (nn<=1) then
                return
            else if (nn==2) then    
                ii = indx(1)
                if (bb) then
                    if (array(2)>array(1)) then
                        indx(1) = indx(2)
                        indx(2) = ii
                    end if
                else
                    if (array(1)>array(2)) then
                        indx(1) = indx(2) 
                        indx(2) = ii     
                    end if                
                end if
            else
                call quickSort_indx_r64a( array,indx,bb,1,nn )
            end if
            return
        end subroutine quickSort_indx_r64
        
    !---
        
        recursive subroutine quickSort_i32a( array,back,left,right )
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^       
            integer(kind=int32),dimension(:),intent(inout)          ::      array
            logical,intent(in)                                      ::      back
            integer,intent(in)                                      ::      left,right
            integer             ::      ii
                        
            call partition_i32( array,back,left,right,ii )
            if (left < ii-1) call quickSort_i32a( array,back,left,ii-1 )
            if (ii < right)   call quickSort_i32a( array,back,ii,right )
            return
        end subroutine quickSort_i32a
        
        recursive subroutine quickSort_i64a( array,back,left,right )
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^       
            integer(kind=int64),dimension(:),intent(inout)          ::      array
            logical,intent(in)                                      ::      back
            integer,intent(in)                                      ::      left,right
            integer             ::      ii
            call partition_i64( array,back,left,right,ii )
            if (left < ii-1) call quickSort_i64a( array,back,left,ii-1 )
            if (ii < right)   call quickSort_i64a( array,back,ii,right )
            return
        end subroutine quickSort_i64a
        
        recursive subroutine quickSort_r32a( array,back,left,right )
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^       
            real(kind=real32),dimension(:),intent(inout)            ::      array
            logical,intent(in)                                      ::      back
            integer,intent(in)                                      ::      left,right
            integer             ::      ii
            call partition_r32( array,back,left,right,ii )
            if (left < ii-1) call quickSort_r32a( array,back,left,ii-1 )
            if (ii < right)   call quickSort_r32a( array,back,ii,right )
            return
        end subroutine quickSort_r32a
        
        recursive subroutine quickSort_r64a( array,back,left,right )
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^       
            real(kind=real64),dimension(:),intent(inout)            ::      array
            logical,intent(in)                                      ::      back
            integer,intent(in)                                      ::      left,right
            integer             ::      ii
            call partition_r64( array,back,left,right,ii )
            if (left < ii-1) call quickSort_r64a( array,back,left,ii-1 )
            if (ii < right)  call quickSort_r64a( array,back,ii,right )
            return
        end subroutine quickSort_r64a
                      
    !---                                
        
        recursive subroutine quickSort_indx_i32a( array,indx,back,left,right )
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^   
            integer(kind=int32),dimension(:),intent(in)             ::      array
            integer,dimension(:),intent(inout)                      ::      indx            
            logical,intent(in)                                      ::      back
            integer,intent(in)                                      ::      left,right
            integer             ::      ii
!             print *,"quickSort_indx_i32a ",back,left,right,ii
            call partition_indx_i32( array,indx,back,left,right,ii )
            if (left < ii-1) call quickSort_indx_i32a( array,indx,back,left,ii-1 )
            if (ii < right)  call quickSort_indx_i32a( array,indx,back,ii,right )
            return
        end subroutine quickSort_indx_i32a
               
        recursive subroutine quickSort_indx_i64a( array,indx,back,left,right )
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^   
            integer(kind=int64),dimension(:),intent(in)             ::      array
            integer,dimension(:),intent(inout)                      ::      indx            
            logical,intent(in)                                      ::      back
            integer,intent(in)                                      ::      left,right
            integer             ::      ii
            call partition_indx_i64( array,indx,back,left,right,ii )
            if (left < ii-1) call quickSort_indx_i64a( array,indx,back,left,ii-1 )
            if (ii < right)  call quickSort_indx_i64a( array,indx,back,ii,right )
            return
        end subroutine quickSort_indx_i64a     
               
        recursive subroutine quickSort_indx_r32a( array,indx,back,left,right )
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^   
            real(kind=real32),dimension(:),intent(in)               ::      array
            integer,dimension(:),intent(inout)                      ::      indx            
            logical,intent(in)                                      ::      back
            integer,intent(in)                                      ::      left,right
            integer             ::      ii
            call partition_indx_r32( array,indx,back,left,right,ii )
            if (left < ii-1) call quickSort_indx_r32a( array,indx,back,left,ii-1 )
            if (ii < right)  call quickSort_indx_r32a( array,indx,back,ii,right )
            return
        end subroutine quickSort_indx_r32a
                
        recursive subroutine quickSort_indx_r64a( array,indx,back,left,right )
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^   
            real(kind=real64),dimension(:),intent(in)               ::      array
            integer,dimension(:),intent(inout)                      ::      indx            
            logical,intent(in)                                      ::      back
            integer,intent(in)                                      ::      left,right
            integer             ::      ii
            call partition_indx_r64( array,indx,back,left,right,ii )
            if (left < ii-1) call quickSort_indx_r64a( array,indx,back,left,ii-1 )
            if (ii < right)  call quickSort_indx_r64a( array,indx,back,ii,right )
            return
        end subroutine quickSort_indx_r64a
        
    !---           
       
        subroutine partition_i32( array,back,left,right,i )
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
            integer(kind=int32),dimension(:),intent(inout)  ::      array
            logical,intent(in)                              ::      back
            integer,intent(in)                              ::      left,right
            integer,intent(out)                             ::      i
            integer                     ::      j
            integer(kind=int32)         ::      tmp,pivot
            
            i = left ; j = right
            pivot = array( (i+j)/2 )
            if (back) then
                do
                    if (array(i)>pivot) then
                        i = i + 1 ; cycle
                    end if
                    if (array(j)<pivot) then
                        j = j - 1 ; cycle
                    end if                       
                    if (i <= j) then
                        tmp = array(i) ; array(i) = array(j) ; array(j) = tmp
                        i = i + 1 ; j = j - 1
                    end if                    
                    if (i>j) exit
                end do            
            else
                do
                    if (array(i)<pivot) then
                        i = i + 1 ; cycle
                    end if
                    if (array(j)>pivot) then
                        j = j - 1 ; cycle
                    end if                       
                    if (i <= j) then
                        tmp = array(i) ; array(i) = array(j) ; array(j) = tmp
                        i = i + 1 ; j = j - 1
                    end if                    
                    if (i>j) exit
                end do
            end if
            return          
        end subroutine partition_i32
        
        subroutine partition_i64( array,back,left,right,i )
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
            integer(kind=int64),dimension(:),intent(inout)  ::      array
            logical,intent(in)                              ::      back
            integer,intent(in)                              ::      left,right
            integer,intent(out)                             ::      i
            integer                     ::      j
            integer(kind=int64)         ::      tmp,pivot
            
            i = left ; j = right
            pivot = array( (i+j)/2 )
            if (back) then
                do
                    if (array(i)>pivot) then
                        i = i + 1 ; cycle
                    end if
                    if (array(j)<pivot) then
                        j = j - 1 ; cycle
                    end if                       
                    if (i <= j) then
                        tmp = array(i) ; array(i) = array(j) ; array(j) = tmp
                        i = i + 1 ; j = j - 1
                    end if                    
                    if (i>j) exit
                end do            
            else
                do
                    if (array(i)<pivot) then
                        i = i + 1 ; cycle
                    end if
                    if (array(j)>pivot) then
                        j = j - 1 ; cycle
                    end if                       
                    if (i <= j) then
                        tmp = array(i) ; array(i) = array(j) ; array(j) = tmp
                        i = i + 1 ; j = j - 1
                    end if                    
                    if (i>j) exit
                end do
            end if
            return          
        end subroutine partition_i64
        
        
        subroutine partition_r32( array,back,left,right,i )
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
            real(kind=real32),dimension(:),intent(inout)    ::      array
            logical,intent(in)                              ::      back
            integer,intent(in)                              ::      left,right
            integer,intent(out)                             ::      i
            integer                     ::      j
            real(kind=real32)           ::      tmp,pivot
            
            i = left ; j = right
            pivot = array( (i+j)/2 )
            if (back) then
                do
                    if (array(i)>pivot) then
                        i = i + 1 ; cycle
                    end if
                    if (array(j)<pivot) then
                        j = j - 1 ; cycle
                    end if                       
                    if (i <= j) then
                        tmp = array(i) ; array(i) = array(j) ; array(j) = tmp
                        i = i + 1 ; j = j - 1
                    end if                    
                    if (i>j) exit
                end do            
            else
                do
                    if (array(i)<pivot) then
                        i = i + 1 ; cycle
                    end if
                    if (array(j)>pivot) then
                        j = j - 1 ; cycle
                    end if                       
                    if (i <= j) then
                        tmp = array(i) ; array(i) = array(j) ; array(j) = tmp
                        i = i + 1 ; j = j - 1
                    end if                    
                    if (i>j) exit
                end do
            end if
            return          
        end subroutine partition_r32
        
        
        subroutine partition_r64( array,back,left,right,i )
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
            real(kind=real64),dimension(:),intent(inout)    ::      array
            logical,intent(in)                              ::      back
            integer,intent(in)                              ::      left,right
            integer,intent(out)                             ::      i
            integer                     ::      j
            real(kind=real64)           ::      tmp,pivot
            
            i = left ; j = right
            pivot = array( (i+j)/2 )
            if (back) then
                do
                    if (array(i)>pivot) then
                        i = i + 1 ; cycle
                    end if
                    if (array(j)<pivot) then
                        j = j - 1 ; cycle
                    end if                       
                    if (i <= j) then
                        tmp = array(i) ; array(i) = array(j) ; array(j) = tmp
                        i = i + 1 ; j = j - 1
                    end if                    
                    if (i>j) exit
                end do            
            else
                do
                    if (array(i)<pivot) then
                        i = i + 1 ; cycle
                    end if
                    if (array(j)>pivot) then
                        j = j - 1 ; cycle
                    end if                       
                    if (i <= j) then
                        tmp = array(i) ; array(i) = array(j) ; array(j) = tmp
                        i = i + 1 ; j = j - 1
                    end if                    
                    if (i>j) exit
                end do
            end if
            return          
        end subroutine partition_r64
        
    !---        
            
        subroutine partition_indx_i32( array,indx,back,left,right,i )
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
            integer(kind=int32),dimension(:),intent(in)     ::      array
            integer,dimension(:),intent(inout)              ::      indx
            logical,intent(in)                              ::      back
            integer,intent(in)                              ::      left,right
            integer,intent(out)                             ::      i
            integer                     ::      tmp,j
            integer(kind=int32)         ::      pivot
            
            i = left ; j = right
            pivot = array( indx((i+j)/2) )
            if (back) then
                do
                    if (array( indx(i) )>pivot) then
                        i = i + 1 ; cycle
                    end if
                    if (array( indx(j) )<pivot) then
                        j = j - 1 ; cycle
                    end if
                        
                    if (i <= j) then
                        tmp = indx(i) ; indx(i) = indx(j) ; indx(j) = tmp
                        i = i + 1 ; j = j - 1
                    end if
                    
                    if (i>j) exit
                end do            
            else
                do
!                 
!                     print *,"partition_indx_i32 i,j,pivot ",i,j,pivot
!                     print *,"indx(i,j) ",indx(i),indx(j)
!                     print *,"array(indx(i,j)) ",array( indx(i) ),array( indx(j) )
!                 
                    if (array( indx(i) )<pivot) then
                        i = i + 1 ; cycle
                    end if
                    if (array( indx(j) )>pivot) then
                        j = j - 1 ; cycle
                    end if
                        
                    if (i <= j) then
                        tmp = indx(i) ; indx(i) = indx(j) ; indx(j) = tmp
                        i = i + 1 ; j = j - 1
                    end if
                    
                    if (i>j) exit
                end do
            end if
            return
        end subroutine partition_indx_i32        
        

        subroutine partition_indx_i64( array,indx,back,left,right,i )
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
            integer(kind=int64),dimension(:),intent(in)     ::      array
            integer,dimension(:),intent(inout)              ::      indx
            logical,intent(in)                              ::      back
            integer,intent(in)                              ::      left,right
            integer,intent(out)                             ::      i
            integer                     ::      tmp,j
            integer(kind=int64)         ::      pivot
            
            i = left ; j = right
            pivot = array( indx((i+j)/2) )
            if (back) then
                do
                    if (array( indx(i) )>pivot) then
                        i = i + 1 ; cycle
                    end if
                    if (array( indx(j) )<pivot) then
                        j = j - 1 ; cycle
                    end if
                        
                    if (i <= j) then
                        tmp = indx(i) ; indx(i) = indx(j) ; indx(j) = tmp
                        i = i + 1 ; j = j - 1
                    end if
                    
                    if (i>j) exit
                end do            
            else
                do
                    if (array( indx(i) )<pivot) then
                        i = i + 1 ; cycle
                    end if
                    if (array( indx(j) )>pivot) then
                        j = j - 1 ; cycle
                    end if
                        
                    if (i <= j) then
                        tmp = indx(i) ; indx(i) = indx(j) ; indx(j) = tmp
                        i = i + 1 ; j = j - 1
                    end if
                    
                    if (i>j) exit
                end do
            end if
            return
        end subroutine partition_indx_i64        
        
        subroutine partition_indx_r32( array,indx,back,left,right,i )
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
            real(kind=real32),dimension(:),intent(in)       ::      array
            integer,dimension(:),intent(inout)              ::      indx
            logical,intent(in)                              ::      back
            integer,intent(in)                              ::      left,right
            integer,intent(out)                             ::      i
            integer                     ::      tmp,j
            real(kind=real32)           ::      pivot
            
            i = left ; j = right
            pivot = array( indx((i+j)/2) )
            if (back) then
                do
                    if (array( indx(i) )>pivot) then
                        i = i + 1 ; cycle
                    end if
                    if (array( indx(j) )<pivot) then
                        j = j - 1 ; cycle
                    end if
                        
                    if (i <= j) then
                        tmp = indx(i) ; indx(i) = indx(j) ; indx(j) = tmp
                        i = i + 1 ; j = j - 1
                    end if
                    
                    if (i>j) exit
                end do            
            else
                do
                    if (array( indx(i) )<pivot) then
                        i = i + 1 ; cycle
                    end if
                    if (array( indx(j) )>pivot) then
                        j = j - 1 ; cycle
                    end if
                        
                    if (i <= j) then
                        tmp = indx(i) ; indx(i) = indx(j) ; indx(j) = tmp
                        i = i + 1 ; j = j - 1
                    end if
                    
                    if (i>j) exit
                end do
            end if
            return
        end subroutine partition_indx_r32               

        subroutine partition_indx_r64( array,indx,back,left,right,i )
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
            real(kind=real64),dimension(:),intent(in)       ::      array
            integer,dimension(:),intent(inout)              ::      indx
            logical,intent(in)                              ::      back
            integer,intent(in)                              ::      left,right
            integer,intent(out)                             ::      i
            integer                     ::      tmp,j
            real(kind=real64)           ::      pivot
            
            i = left ; j = right
            pivot = array( indx((i+j)/2) )
            if (back) then
                do
                    if (array( indx(i) )>pivot) then
                        i = i + 1 ; cycle
                    end if
                    if (array( indx(j) )<pivot) then
                        j = j - 1 ; cycle
                    end if
                        
                    if (i <= j) then
                        tmp = indx(i) ; indx(i) = indx(j) ; indx(j) = tmp
                        i = i + 1 ; j = j - 1
                    end if
                    
                    if (i>j) exit
                end do            
            else
                do
                    if (array( indx(i) )<pivot) then
                        i = i + 1 ; cycle
                    end if
                    if (array( indx(j) )>pivot) then
                        j = j - 1 ; cycle
                    end if
                        
                    if (i <= j) then
                        tmp = indx(i) ; indx(i) = indx(j) ; indx(j) = tmp
                        i = i + 1 ; j = j - 1
                    end if
                    
                    if (i>j) exit
                end do
            end if
            return
        end subroutine partition_indx_r64        
        
                
    end module Lib_quicksort                
    
    
!     
!     program testLib_quicksort
! !---^^^^^^^^^^^^^^^^^^^^^^^^^
! 
!         use Lib_quicksort
!         implicit none
!         
!         integer,dimension(12)           ::      array_1 = (/ 4,5,78,21,4,7,9,23,5,6,2,1 /)
!         real(kind=8),dimension(12)      ::      array_2 = (/ 4.0d0,5.1d0,78.2d0,21.2d0,4.3d0,7.4d0,9.5d0,23.1d0,5.1d0,6.8d0,2.9d0,1.0d0 /)
!         integer,dimension(size(array_2))           ::      indx_2 
!         integer         ::      ii
!         
!         print *,"sort int"
!         print *,array_1
!         call quicksort(array_1)
!         print *,array_1
! !         
!         print *,""
!         print *,"sort real"
!         print *,array_2
!         call quicksort(array_2 )
!         print *,array_2
!         
!         
!         
!         
!         print *,""
!         print *,"sort index backwards"
! !         array_1 = (/ 4,5,78,21,4,7,9,23,5,6,2,1 /)
!         do ii = 1,size(array_2)
!             indx_2(ii) = ii
!         end do
!         call quicksort(array_2,indx_2,back=.true.)
!         do ii = 1,size(array_2)
!             print *,ii,array_2(ii),indx_2(ii),array_2(indx_2(ii))
!         end do
!         
!         print *,""    
!         print *,"done"            
!         print *,""    
! 
!    end program testLib_quicksort               
!
