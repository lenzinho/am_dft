module am_stdout

    use am_constants
    use am_matlab

    implicit none

    interface am_print
        ! function overloading must start with the same name as the interface, i.e. "am_print" in this case
        module procedure &
            am_print_str,       &
            am_print_int,       &
            am_print_int_vec,   &
            am_print_int_mat,   &
            am_print_dble,      &
            am_print_dble_vec,  &
            am_print_dble_mat,  &
            am_print_cmplx_vec, &
            am_print_cmplx_mat, &
            am_print_bool,      &
            am_print_bool_mat
    end interface ! am_print
    
    interface am_print_sparse
        module procedure am_print_sparse_double, am_print_sparse_integer
    end interface ! am_print_sparse

contains

    !
    ! convert data types  
    !

    function      int2char(int,iopt_sign) result(str)
        !
        implicit none
        !
        integer, intent(in) :: int
        character(len=*), intent(in), optional :: iopt_sign
        character(100) :: str
        !
        if (.not.present(iopt_sign)) then
            write(str,'(i)') int
            str = trim( adjustl( str ) )
        else
            if (index(iopt_sign,'SP').ge.1) then 
                write(str,'(SP,i)') int
                str = trim( adjustl( str ) )
            else
                call am_print('ERROR','Unknown sign requested.')
                stop
            endif
        endif
        !
    end function  int2char
    
    function      int2charSP(int) result(str)
        !
        implicit none
        !
        integer, intent(in) :: int
        character(100) :: str
        !
        write(str,'(SP,i)') int
        !
    end function  int2charSP

    function      dbl2char(dbl,iopt_len) result(str)
        !
        implicit none
        !
        real(dp), intent(in) :: dbl
        integer, intent(in), optional :: iopt_len
        character(:), allocatable :: str
        character(100) :: buffer
        !
        if (present(iopt_len)) then
            allocate(character(iopt_len) :: str)
            write(buffer,'(f)') nint(dbl*(10.0_dp**iopt_len))/(10.0_dp**iopt_len)
        else
            allocate(character(100) :: str)
            write(buffer,'(f)') dbl
        endif
        !
        write(buffer,'(f)') dbl
        buffer = trim( adjustl( buffer ) )
        str = buffer(1:len(str))
        !
    end function  dbl2char

    function      dbl2charSP(dbl,iopt_len) result(str)
        !
        implicit none
        !
        real(dp), intent(in) :: dbl
        integer, intent(in), optional :: iopt_len
        character(:), allocatable :: str
        character(100) :: buffer
        !
        if (present(iopt_len)) then
            allocate(character(iopt_len) :: str)
            write(buffer,'(SP,f)') nint(dbl*(10.0_dp**iopt_len))/(10.0_dp**iopt_len)
        else
            allocate(character(100) :: str)
            write(buffer,'(SP,f)') dbl
        endif
        !
        buffer = trim( adjustl( buffer ) )
        str = buffer(1:len(str))
        !
    end function  dbl2charSP

    function      dbl2frac(dbl) result(nd)
        !
        implicit none
        !
        real(dp), intent(in) :: dbl
        integer :: n, c, d
        integer :: nd(2)
        !
        d = 100
        n = nint(dbl*d)
        c = gcd(n,d)
        n = n/c
        d = d/c
        nd(1) = n ! numerator
        nd(2) = d ! denominator
        !
    end function  dbl2frac

    function      print_pretty(dbl) result(str)
        !
        ! call it with trim
        !
        implicit none
        !
        real(dp), intent(in) :: dbl
        character(len=column_width) :: str
        integer :: num_den(2)
        !
        !
        ! convert double to 0
        if (abs(dbl).lt.tiny) then
            str='0'
            return
        endif
        ! convert double to integer
        if (abs(nint(dbl)-dbl).lt.tiny) then
            str=int2char(nint(dbl))
            return
        endif
        ! convert double to fraction
        num_den = dbl2frac(dbl)
        if (abs(num_den(2)).gt.tiny) then
        if (abs(num_den(1)/real(num_den(2),dp)-dbl).lt.tiny) then
            write(str,'(a,a,a)') trim(int2char(num_den(1))), '/', trim(int2char(num_den(2)))
            return
        endif
        endif
        ! do nothing
        write(str,*) dbl
        !
    end function  print_pretty

    subroutine     am_print_two_matrices_side_by_side(name, Atitle, Btitle, A, B , iopt_emph, iopt_fid , iopt_teaser )
        !
        implicit none
        !
        character(*), intent(in) :: name
        character(*), intent(in) :: Atitle
        character(*), intent(in) :: Btitle
        real(dp)    , intent(in) :: A(:,:)
        real(dp)    , intent(in) :: B(:,:)
        logical     , intent(in), optional :: iopt_teaser
        integer     , intent(in), optional :: iopt_fid
        character(5), intent(in), optional :: iopt_emph
        integer, parameter :: maximum_lines_to_print = 10
        logical :: teaser
        character(5) :: emph
        integer :: fid
        integer :: i, j, m, n
        !
        fid=6
        if ( present(iopt_fid)  ) fid = iopt_fid
        emph = "     "
        if ( present(iopt_emph) ) emph = iopt_emph
        teaser = .false.
        if ( present(iopt_teaser) ) teaser = iopt_teaser
        !
        !
        !
        if (size(A,1).ne.size(B,1)) then
            call am_print('ERROR','Size mismatch: A and B.',flags='E')
            stop
        endif
        if (size(A,2).ne.size(B,2)) then
            call am_print('ERROR','Size mismatch: A and B.',flags='E')
            stop
        endif
        !
        m=size(A,1)
        n=size(A,2)
        !
        if (teaser) m = minval([maximum_lines_to_print,m])
        !
        write(fid,'(a,a," = ")') emph, name
        !
        write(fid,'(5x,a)',advance='no') centertitle(Atitle,n*15)
        write(fid,'(5x,a)',advance='no') centertitle(Btitle,n*15)
        write(fid,*)
        !
        write(fid,'(5x,a)',advance='no') centertitle(repeat('-',n*15),n*15)
        write(fid,'(5x,a)',advance='no') centertitle(repeat('-',n*15),n*15)
        write(fid,*)
        !
        do i = 1, m
            !
            write(fid,'(5x)',advance='no')
            do j = 1,n
                write(fid,'(f15.8)',advance='no') A(i,j)
            enddo
            !
            write(fid,'(5x)',advance='no')
            do j = 1,n
                write(fid,'(f15.8)',advance='no') B(i,j)
            enddo
            !
            write(fid,*)
        enddo
        !
        if (size(A,1).gt.maximum_lines_to_print) then
            write(fid,'(5x,a)',advance='no') centertitle(flare,n*15)
            write(fid,'(5x,a)',advance='no') centertitle(flare,n*15)
            write(fid,*)
        endif
        !
    end subroutine am_print_two_matrices_side_by_side

    !
    ! generic math functions 
    !

    pure function regspace(d1,d2,d) result(y)
        !
        implicit none
        !
        real(dp), intent(in) :: d1
        real(dp), intent(in) :: d2
        real(dp), intent(in) :: d
        real(dp), allocatable :: y(:)
        integer :: n
        !
        n=nint((d2-d1)/d)
        !
        allocate(y(n))
        y = d1+d*[0:n:1]
        !
    end function  regspace

    !
    ! print related stuff
    !
    
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

    subroutine     am_print_bool(name,A,flags,permission,filename)
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
    end subroutine am_print_bool

    subroutine     am_print_dble(name,A,flags,permission,filename)
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
    end subroutine am_print_dble

    subroutine     am_print_int(name,A,flags,permission,filename)
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
    end subroutine am_print_int

    subroutine     am_print_str(name,A,flags,permission,filename)
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
    end subroutine am_print_str

    subroutine     am_print_bool_mat(name,A,flags,permission,filename)
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
    end subroutine am_print_bool_mat

    subroutine     am_print_dble_mat(name,A,flags,permission,filename)
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
    end subroutine am_print_dble_mat

    subroutine     am_print_cmplx_mat(name,A,flags,permission,filename)
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
        do j = 1, size(A,2)
            write(unit=fid,fmt='(5x,2f10.2,"i ")',advance='no') real(A(i,j)),aimag(A(i,j))
        enddo
        write(unit=fid,fmt=*)
        enddo
        !
#include 'am_print_footer.inc'
        !
    end subroutine am_print_cmplx_mat

    subroutine     am_print_int_mat(name,A,flags,permission,filename)
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
    end subroutine am_print_int_mat

    subroutine     am_print_dble_vec(name,A,flags,permission,filename)
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
    end subroutine am_print_dble_vec

    subroutine     am_print_int_vec(name,A,flags,permission,filename)
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
    end subroutine am_print_int_vec

    subroutine     am_print_cmplx_vec(name,A,flags,permission,filename)
        !
        implicit none
        !
        complex(dp), intent(in) :: A(:)
        integer :: i
        !
#include 'am_print_header.inc'
        !
        if (fid.eq.6) write(unit=fid,fmt='(a,a," = ",a)') emph, name
        do i = 1, size(A,1)
            write(unit=fid,fmt='(5x,2f10.2,"i ")',advance='no') real(A(i)),aimag(A(i))
            write(unit=fid,fmt=*)
        enddo
        !
#include 'am_print_footer.inc'
        !
    end subroutine am_print_cmplx_vec





    end module







