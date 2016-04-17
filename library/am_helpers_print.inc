    pure function  centertitle(title,length) result(title_centered)
        !
        implicit none
        !
        character(len=*), intent(in) :: title
        integer, intent(in) :: length
        character(len=length) :: title_centered
        integer :: blanks
        integer :: i 
        !
        title_centered = title
        blanks = 0 
        do  i = length,1,-1                  ! starting on the right side 
        if( title_centered(i:i) .ne. ' ' ) goto 33    ! check for trailing blanks 
        blanks = blanks + 1                  ! count the blanks 
        enddo 
        33  continue 
        if ( blanks .gt. 1 ) then 
        blanks = blanks/2                      ! cut half of the blanks 
        do i = length , 1 , -1 
          if ( i-blanks .gt. 0 ) title_centered(i:i) = title_centered(i-blanks:i-blanks) 
          if ( i-blanks .le. 0 ) title_centered(i:i) = ' '     ! blank front half 
        enddo
        endif
    end function   centertitle

    subroutine     am_print_title(title)
        !
        implicit none
        !
        character(len=*), intent(in) :: title
        integer :: fid
        !
        fid = 6
        !
        write(unit=fid,fmt='(a,a,a)') "+", repeat("-",column_width-2), "+"
        write(unit=fid,fmt='(a,a,a)') "|",  centertitle(title,column_width-2), "|" 
        write(unit=fid,fmt='(a,a,a)') "+", repeat("-",column_width-2), "+"
        !
    end subroutine am_print_title
    
    subroutine     am_print_sparse_double(name,A,flags,permission,filename)
        !
        implicit none
        !
        real(dp), intent(in) :: A(:,:)
        integer :: i, j
        !
#include 'am_print_header.inc'
        !
        if (fid.eq.6) write(unit=fid,fmt='(a,a," = ")') emph, trim(name)
        write(unit=fid,fmt='(5x,a,a,a)') ,'+', repeat('-',size(A,2)), '+'
        !
        do i = 1, size(A,1)
            write(unit=fid,fmt='(5x)',advance='no')
            write(unit=fid,fmt='(a1)',advance='no') '|'
            do j = 1, size(A,2)
                if(abs(A(i,j)).gt.tiny) then
                    write(*,'(a1)',advance='no') 'x'
                else
                    write(*,'(a1)',advance='no') ' '
                endif
            enddo
            write(unit=fid,fmt='(a1)',advance='no') '|'
            write(unit=fid,fmt=*)
        enddo
        !
        write(unit=fid,fmt='(5x,a,a,a)') ,'+', repeat('-',size(A,2)), '+'
        !
#include 'am_print_footer.inc'
        !
    end subroutine am_print_sparse_double

    subroutine     am_print_sparse_integer(name,A,flags,permission,filename)
        !
        implicit none
        !
        integer, intent(in) :: A(:,:)
        integer :: i, j
        !
#include 'am_print_header.inc'
        !
        if (fid.eq.6) write(unit=fid,fmt='(a,a," = ")') emph, trim(name)
        write(unit=fid,fmt='(5x,a,a,a)') ,'+', repeat('-',size(A,2)), '+'
        !
        do i = 1, size(A,1)
            write(unit=fid,fmt='(5x)',advance='no')
            write(unit=fid,fmt='(a1)',advance='no') '|'
            do j = 1, size(A,2)
                if(A(i,j).ne.0) then
                    write(*,'(a1)',advance='no') 'x'
                else
                    write(*,'(a1)',advance='no') ' '
                endif
            enddo
            write(unit=fid,fmt='(a1)',advance='no') '|'
            write(unit=fid,fmt=*)
        enddo
        !
        write(unit=fid,fmt='(5x,a,a,a)') ,'+', repeat('-',size(A,2)), '+'
        !
#include 'am_print_footer.inc'
        !
    end subroutine am_print_sparse_integer

    subroutine     am_print_logical(name,A,flags,permission,filename)
        !
        implicit none
        !
        logical, intent(in) :: A
        !
#include 'am_print_header.inc'
        !
        write(unit=fid,fmt='(a,a," = ",l)') emph, name, A
        !
#include 'am_print_footer.inc'
        !
    end subroutine am_print_logical

    subroutine     am_print_double(name,A,flags,permission,filename)
        !
        implicit none
        !
        real(dp), intent(in) :: A
        !
#include 'am_print_header.inc'
        !
        write(unit=fid,fmt='(a,a," = ",g15.5)') emph, name, dbl2char(A)
        !
#include 'am_print_footer.inc'
        !
    end subroutine am_print_double

    subroutine     am_print_integer(name,A,flags,permission,filename)
        !
        implicit none
        !
        integer, intent(in) :: A
        !
#include 'am_print_header.inc'
        !
        if (fid.eq.6) write(unit=fid,fmt='(a,a," = ",a)') emph, name, int2char(A)
        !
#include 'am_print_footer.inc'
        !
    end subroutine am_print_integer

    subroutine     am_print_string(name,A,flags,permission,filename)
        !
        implicit none
        !
        character(len=*), intent(in) :: A
        !
#include 'am_print_header.inc'
        !
        if (fid.eq.6) write(unit=fid,fmt='(a,a," = ",a)') emph, name, A
        !
#include 'am_print_footer.inc'
        !
    end subroutine am_print_string

    subroutine     am_print_logical_matrix(name,A,flags,permission,filename)
        !
        implicit none
        !
        logical, intent(in) :: A(:,:)
        integer :: i, j
        !
#include 'am_print_header.inc'
        !
        if (fid.eq.6) write(unit=fid,fmt='(a,a," = ")') emph, trim(name)
        do i = 1, size(A,1)
            write(unit=fid,fmt="(5x,300l)") ( A(i,j), j = 1, size(A,2) )
        enddo
        !
#include 'am_print_footer.inc'
        !
    end subroutine am_print_logical_matrix

    subroutine     am_print_double_matrix(name,A,flags,permission,filename)
        !
        implicit none
        !
        real(dp), intent(in) :: A(:,:)
        integer :: i, j
        !
#include 'am_print_header.inc'
        !
        if (fid.eq.6) write(unit=fid,fmt='(a,a," = ")') emph, trim(name)
        do i = 1, size(A,1)
            write(unit=fid,fmt='(5x)',advance='no')
            do j = 1, size(A,2)
            write(unit=fid,fmt='(f15.5)',advance='no') A(i,j)
            enddo
            write(unit=fid,fmt=*)
        enddo
        !
#include 'am_print_footer.inc'
        !
    end subroutine am_print_double_matrix

    subroutine     am_print_complex_matrix(name,A,flags,permission,filename)
        !
        implicit none
        !
        complex(dp), intent(in) :: A(:,:)
        integer :: i, j
        !
#include 'am_print_header.inc'
        !
        if (fid.eq.6) write(unit=fid,fmt='(a,a," = ",a)') emph, name
        do i = 1, size(A,1)
            write(unit=fid,fmt="(5x,100f15.5)") ( A(i,j), j = 1, size(A,2) )
        enddo
        !
#include 'am_print_footer.inc'
        !
    end subroutine am_print_complex_matrix

    subroutine     am_print_int_matrix(name,A,flags,permission,filename)
        !
        implicit none
        !
        integer, intent(in) :: A(:,:)
        integer :: i, j
        !
#include 'am_print_header.inc'
        !
        if (fid.eq.6) write(unit=fid,fmt='(a,a," = ",a)') emph, name
        do i = 1, size(A,1)
            write(unit=fid,fmt="(5x,300i4)") ( A(i,j), j = 1, size(A,2) )
        enddo
        !
#include 'am_print_footer.inc'
        !
    end subroutine am_print_int_matrix

    subroutine     am_print_double_vector(name,A,flags,permission,filename)
        !
        implicit none
        !
        real(dp), intent(in) :: A(:)
        integer :: i
        !
#include 'am_print_header.inc'
        !
        if (fid.eq.6) write(unit=fid,fmt='(a,a," = ",a)') emph, name
        write(unit=fid,fmt="(5x,100f15.5)") ( A(i), i = 1, size(A,1) )
        !
#include 'am_print_footer.inc'
        !
    end subroutine am_print_double_vector

    subroutine     am_print_int_vector(name,A,flags,permission,filename)
        !
        implicit none
        !
        integer, intent(in) :: A(:)
        integer :: i
        !
#include 'am_print_header.inc'
        !
        if (fid.eq.6) write(unit=fid,fmt='(a,a," = ",a)') emph, name
        write(unit=fid,fmt="(5x,300i4)") ( A(i), i = 1, size(A,1) )
        !
#include 'am_print_footer.inc'
        !
    end subroutine am_print_int_vector




