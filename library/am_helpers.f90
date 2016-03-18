    module am_helpers
    !
    use am_constants
    !
    implicit none
    !
    interface unique
        module procedure unique_integer, unique_columns_integer, unique_columns_double, unique_matrices_double
    end interface
    !
    interface am_print
        ! function overloading must start with the same name as the interface, i.e. "am_print" in this case
        module procedure am_print_double, &
            am_print_integer, &
            am_print_string, &
            am_print_logical, &
            am_print_double_matrix, &
            am_print_complex_matrix, &
            am_print_int_matrix, &
            am_print_logical_matrix, &
            am_print_double_vector, &
            am_print_int_vector
    end interface ! am_print
    !
    interface print_sparse
        module procedure print_sparse_double, print_sparse_integer
    end interface ! print_sparse
    !
    interface linspace
        module procedure linspace_double, linspace_integer
    end interface ! linspace
    !
    interface inv
        module procedure inv_3x3_dbl
    end interface ! inv
    !
    interface det
        module procedure det_3x3_dbl
    end interface ! det
    !
    interface rref
        module procedure rref_dbl, rref_int
    end interface ! rref

    contains

    !
    ! PRINT TO STDOUT
    !

#include "am_helpers_print.f90"
 !
    ! CONVERT DATA TYPES
    !

    function      int2char(int) result(string)
        !
        implicit none
        !
        integer :: int
        character(100) :: string
        !
        write(string,'(i)') int
        string = trim( adjustl( string ) )
    end function  int2char
    
    function      dbl2char(dbl) result(string)
        !
        implicit none
        !
        real(dp) :: dbl
        character(100) :: string
        !
        write(string,'(f)') dbl
        string = trim( adjustl( string ) )
        !
    end function  dbl2char

    function      dbl2charSP(dbl,iopt_len) result(string)
        !
        implicit none
        !
        real(dp), intent(in) :: dbl
        integer, intent(in), optional :: iopt_len
        character(100) :: string
        !
        write(string,'(SP,f)') dbl
        string = trim( adjustl( string ) )
        if (present(iopt_len)) string = string(1:iopt_len)
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

    function      print_pretty(dbl) result(string)
        !
        ! call it with trim
        !
        implicit none
        !
        real(dp), intent(in) :: dbl
        character(len=column_width) :: string
        integer :: num_den(2)
        !
        !
        ! convert double to 0
        if (abs(dbl).lt.tiny) then
            string='0'
            return
        endif
        ! convert double to integer
        if (abs(nint(dbl)-dbl).lt.tiny) then
            string=int2char(nint(dbl))
            return
        endif
        ! convert double to fraction
        num_den = dbl2frac(dbl)
        if (abs(num_den(2)).gt.tiny) then
        if (abs(num_den(1)/real(num_den(2),dp)-dbl).lt.tiny) then
            write(string,'(a,a,a)') trim(int2char(num_den(1))), '/', trim(int2char(num_den(2)))
            return
        endif
        endif
        ! do nothing
        write(string,*) dbl
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
            call am_print('ERROR','Size mismatch: A and B.',' >>> ')
            stop
        endif
        if (size(A,2).ne.size(B,2)) then
            call am_print('ERROR','Size mismatch: A and B.',' >>> ')
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
            write(fid,'(5x,a)',advance='no') centertitle(' ... ',n*15)
            write(fid,'(5x,a)',advance='no') centertitle(' ... ',n*15)
            write(fid,*)
        endif
        !
    end subroutine am_print_two_matrices_side_by_side

    !
    ! OPEN FILE
    !

    integer function unit_number(filename)
        !> Borrowed from Olle
        character(len=*),intent(in) :: filename
        !
        integer :: unit, i=100
        logical :: file_open

        inquire(file=filename,number=unit,opened=file_open)
        if ( file_open .eqv. .true. ) then
            unit_number=unit
        else
            do
                inquire(unit=i,opened=file_open)
                if ( file_open .eqv. .false. ) then
                    unit_number=i
                    exit
                endif
                i=i+1
                if ( i == 1000 ) then
                    write(*,*) 'No available units between 100 and 1000. Are you sure you want that'
                    write(*,*) 'many open files?'
                endif
            enddo
        endif
    end function unit_number

    integer function open_file(rw,filename)
        !> Borrowed from Olle
        implicit none
        !
        character(len=*),intent(in) :: filename
        character(len=*),intent(in) :: rw !> read/write
        integer :: unit, status
        !
        status=0
        !
        unit=unit_number(filename)
        select case(trim(rw))
        case('r')
            open(unit=unit, file=filename, status='old', action='read',iostat=status)
            open_file=unit
        case('w')
            open(unit=unit, file=filename, status='replace', action='write',iostat=status)
            open_file=unit
            case default
            write(*,*) "Please open the file with either 'r' or 'w', not something in between"
            stop
        end select
        !
        if ( status .ne. 0 ) then
            write(*,*) 'Could not open ',filename
            stop
        endif
    end function open_file 

    !
    ! progress bar
    !

    subroutine     progressbar(message,j,totn)
        !> Borrowed from Olle; original from: https://software.intel.com/en-us/forums/intel-fortran-compiler-for-linux-and-mac-os-x/topic/270155
        integer, intent(in) :: j !> current index
        integer, intent(in) :: totn !> max indes
        character(len=*), intent(in) :: message
        !
        real(dp) :: f0
        integer :: i,k
        character(len=50) :: bar
        character(len=40) :: msg
        !
        ! updates the fraction of calculation done
        !
        bar="     % |                                        |"
        msg="                                        "
        ! write the percentage in the bar
        if ( totn .gt. 1 ) then
            f0=100.0_dp*(j-1.0_dp)/(totn*1.0_dp-1)
        else
            f0=100.0_dp
        endif
        write(unit=bar(1:5),fmt="(F5.1)") f0
        !fill in the bar
        k=ceiling(40*f0/100)
        do i = 1, k
            bar(8+i:8+i)="="
        enddo
        ! add the message
        do i=1,min(len_trim(message),40)
            msg(i:i)=message(i:i)
        enddo
        !
        ! Print the progress bar. Different compilers like it differently.
        !
        !#if ifortprogressbar
        ! write(unit=6,fmt="(1X,a1,a1,a40,a50)") '+',char(13),msg,bar
        ! unit 6 is the screen, stdout
        write(unit=6,fmt="(1X,a1,a1,a40,a50)",advance='no') ' ',char(13),msg,bar
        ! if ( j .eq. totn ) write(*,*) " "
        !#elif gfortranprogressbar
        !    if ( j .lt. totn ) then
        !        write(*,'(1X,a,a,a)',advance='no') char(13),msg,bar
        !    else
        !        write(*,'(1X,a,a,a)') char(13),msg,bar
        !    endif
        !#else
        !    ! Boring version.
        !    if ( totn .gt. 100 ) then
        !        if ( mod(j,totn/100) .eq. 0 ) write(*,*) msg,real(f0),'%'
        !    else
        !        write(*,*) msg,real(f0),'%'
        !    endif
        !#endif
    end subroutine progressbar

    !
    ! MATRIX FUNCTIONS
    !

    pure function eye(n)
        !> nxn identity matrix
        implicit none
        !
        integer, intent(in) :: n
        real(dp), dimension(n,n) :: eye
        integer :: i
        !
        eye=0.0_dp
        do i = 1, n
            eye(i,i) = 1.0_dp
        enddo
    end function  eye

    pure function mesh_grid(n) result(grid_points)
        !> 
        !> Generates a mesh of lattice vectors (fractional coordinates) around the origin. 
        !> n(1), n(2), n(3) specifies the number of mesh points away from 0, i.e. [-n:1:n]
        !> To obtain a mesh of reciprocal lattice vectors: matmul(recbas,grid_points)
        !> To obtain a mesh of real-space lattice vectors: matmul(bas,grid_points)
        !>
        implicit none
        !
        integer , intent(in) :: n(3)
        real(dp), allocatable :: grid_points(:,:) !> voronoi points (27=3^3)
        integer :: i1,i2,i3,j
        !
        allocate(grid_points(3,product(2*n+1)))
        !
        j=0
        do i1 = -n(1),n(1)
            do i2 = -n(2),n(2)
                do i3 = -n(3),n(3)
                    j=j+1
                    grid_points(1:3,j)=[i1,i2,i3]
                enddo
            enddo
        enddo
    end function  mesh_grid

    pure function det_3x3_dbl(a) result(det)
        !
        implicit none
        !
        real(dp), intent(in) :: a(3,3)
        real(dp) :: det
        !
        det = a(1,1)*a(2,2)*a(3,3) &
            - a(1,1)*a(2,3)*a(3,2) &
            - a(1,2)*a(2,1)*a(3,3) &
            + a(1,2)*a(2,3)*a(3,1) &
            + a(1,3)*a(2,1)*a(3,2) &
            - a(1,3)*a(2,2)*a(3,1)
        !
    end function  det_3x3_dbl
    
    pure function inv_3x3_dbl(a) result(ainv)
        !
        implicit none
        !
        real(dp), intent(in) :: a(3,3)
        real(dp) :: cofactor(3,3)
        real(dp) :: determinant
        real(dp) :: ainv(3,3)
        !
        determinant = &
              + a(1,1)*a(2,2)*a(3,3) &
              - a(1,1)*a(2,3)*a(3,2) &
              - a(1,2)*a(2,1)*a(3,3) &
              + a(1,2)*a(2,3)*a(3,1) &
              + a(1,3)*a(2,1)*a(3,2) &
              - a(1,3)*a(2,2)*a(3,1)
        !
        cofactor(1,1) = +(a(2,2)*a(3,3)-a(2,3)*a(3,2))
        cofactor(1,2) = -(a(2,1)*a(3,3)-a(2,3)*a(3,1))
        cofactor(1,3) = +(a(2,1)*a(3,2)-a(2,2)*a(3,1))
        cofactor(2,1) = -(a(1,2)*a(3,3)-a(1,3)*a(3,2))
        cofactor(2,2) = +(a(1,1)*a(3,3)-a(1,3)*a(3,1))
        cofactor(2,3) = -(a(1,1)*a(3,2)-a(1,2)*a(3,1))
        cofactor(3,1) = +(a(1,2)*a(2,3)-a(1,3)*a(2,2))
        cofactor(3,2) = -(a(1,1)*a(2,3)-a(1,3)*a(2,1))
        cofactor(3,3) = +(a(1,1)*a(2,2)-a(1,2)*a(2,1))
        !
        ainv = transpose(cofactor) / determinant
        !
    end function  inv_3x3_dbl

    pure subroutine rref_dbl(A)
        !
        implicit none
        !
        real(dp), dimension(:,:), intent(inout) :: A
        integer :: pivot, norow, nocolumn
        integer :: r, i
        real, dimension(:), allocatable :: trow
        pivot = 1
        norow = size(A, 1)
        nocolumn = size(A, 2)
        !
        allocate(trow(nocolumn))
        !
        do r = 1, norow
            if ( nocolumn <= pivot ) exit
            i = r
            do while ( A(i, pivot) == 0 )
                i = i + 1
                if ( norow == i ) then
                    i = r
                    pivot = pivot + 1
                    if ( nocolumn == pivot ) return
                end if
            end do
            trow = A(i, :)
            A(i, :) = A(r, :)
            A(r, :) = trow
            A(r, :) = A(r, :) / A(r, pivot)
            do i = 1, norow
                if ( i /= r ) A(i, :) = A(i, :) - A(r, :) * A(i, pivot)
            end do
            pivot = pivot + 1
        end do
        deallocate(trow)
    end subroutine  rref_dbl

    pure subroutine rref_int(A)
        !
        implicit none
        !
        integer, dimension(:,:), intent(inout) :: A
        integer :: pivot, norow, nocolumn
        integer :: r, i
        integer, allocatable :: trow(:)
        pivot = 1
        norow = size(A,1)
        nocolumn = size(A,2)
        !
        allocate(trow(nocolumn))
        !
        do r = 1, norow
            if ( nocolumn <= pivot ) exit
            i = r
            do while ( A(i, pivot) == 0 )
                i = i + 1
                if ( norow == i ) then
                    i = r
                    pivot = pivot + 1
                    if ( nocolumn == pivot ) return
                end if
            end do
            trow = A(i, :)
            A(i, :) = A(r, :)
            A(r, :) = trow
            A(r, :) = A(r, :) / A(r, pivot)
            do i = 1, norow
                if ( i /= r ) A(i, :) = A(i, :) - A(r, :) * A(i, pivot)
            end do
            pivot = pivot + 1
        end do
        deallocate(trow)
    end subroutine  rref_int

    !
    ! generic math functions
    !

#include "am_helpers_unique.f90"

    pure function issubset(A,B,iopt_sym_prec)
        !> returns true if true if B(:) is a subset of A(:,i) within numerical precision for any i
        implicit none
        !
        real(dp), intent(in) :: A(:,:)
        real(dp), intent(in) :: B(:)
        logical :: issubset
        integer :: i
        real(dp), intent(in), optional :: iopt_sym_prec
        real(dp) :: sym_prec
        !
        if (present(iopt_sym_prec)) then
            sym_prec = iopt_sym_prec
        else
            sym_prec = tiny
        endif
        !
        issubset = .false.
        !
        do i = 1,size(A,2)
            if ( all(abs(A(:,i)-B).lt.sym_prec) ) then
                issubset = .true.
                return
            endif
        enddo
        !
    end function  issubset

    pure function linspace_double(d1,d2,n) result(y)
        !
        implicit none
        !
        real(dp), intent(in) :: d1
        real(dp), intent(in) :: d2
        integer, intent(in)  :: n
        real(dp) :: d
        real(dp) :: y(n)
        integer :: i
        !
        d = (d2-d1)/(n-1.0_dp)
        do i = 1,n
            y(i) = real(d1,dp) + (i-1.0_dp)*d;
        enddo
        !
    end function  linspace_double

    pure function linspace_integer(d1,d2,n) result(y)
        !
        implicit none
        !
        integer, intent(in) :: d1
        integer, intent(in) :: d2
        integer, intent(in)  :: n
        real(dp) :: d
        real(dp) :: y(n)
        integer :: i
        !
        d = (d2-d1)/(n-1.0_dp)
        do i = 1,n
            y(i) = real(d1,dp) + (i-1.0_dp)*d;
        enddo
        !
    end function  linspace_integer

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

    pure function cross_product(a,b) result(c)
        !
        implicit none
        !
        real(dp) :: c(3)
        real(dp), INTENT(IN) :: a(3), b(3)
        !
        c(1) = a(2) * b(3) - a(3) * b(2)
        c(2) = a(3) * b(1) - a(1) * b(3)
        c(3) = a(1) * b(2) - a(2) * b(1)
        !
    end function  cross_product

    pure function rotation_from_vectors(a,b) result(R)
        ! based on Rodrigues' Rotation Formula
        ! "A Mathematical Introduction to Robotic Manipulation", Richard M. Murray, Zexiang Li, S. Shankar Sastry, pp. 26-28
        implicit none
        !
        real(dp), intent(in) :: a(3)
        real(dp), intent(in) :: b(3)
        real(dp) :: R(3,3)
        real(dp) :: s ! sine of angle
        real(dp) :: c ! cosine of angle
        real(dp) :: v(3) ! cross product of a and b
        real(dp) :: vx(3,3) ! skew symmetric cross product of v
        ! 
        v = cross_product(a,b)
        s = norm2(v)
        c = dot_product(a,b)
        vx(1:3,1) = [0.0_dp,v(3),-v(2)]
        vx(1:3,2) = [-v(3),0.0_dp,v(1)]
        vx(1:3,3) = [v(2),-v(1),0.0_dp]
        !
        R = eye(3) + vx + matmul(vx,vx)*(1.0_dp - c)/(s**2)
        !
    end function  rotation_from_vectors
    
    pure function are_equal(a)
        !
        ! are all elements of a equal to each other
        !
        implicit none
        !
        real(dp), intent(in) :: a(:)
        logical :: are_equal
        !
        if ( all(abs(a(1)-a).lt.tiny) ) then
            are_equal = .true.
        else
            are_equal = .false.
        endif
        !
    end function  are_equal

    pure function are_different(a)
        !
        ! are all elements of a different from each other?
        !
        implicit none
        !
        real(dp), intent(in) :: a(:)
        logical :: are_different
        integer :: i, j, n
        !
        n = size(a)
        !
        if (n.eq.1) then
            are_different = .false.
            return
        endif    
        !
        are_different = .true.
        do i = 1, n
        do j = 1, i
            if (i.ne.j) then
            if (abs(a(i)-a(j)).lt.tiny) then
                are_different = .false.
                return
            endif
            endif
        enddo
        enddo        
    end function  are_different

    pure function lcm(a,b)
        !
        implicit none
        !
        integer, intent(in) :: a,b
        integer :: lcm
        !
        lcm = a*b / gcd(a,b)
        !
    end function  lcm
 
    pure function gcd(a,b)
        !
        implicit none
        !
        integer, intent(in) :: a, b
        integer :: t, bc, ac
        integer :: gcd
        !
        bc = b
        ac = a
        !
        do while (bc/=0)
            t = bc
            bc = mod(ac,bc)
            ac = t
        end do
        gcd = abs(ac)
        !
    end function  gcd

    !
    ! things that should've existed in fortran already
    !
    
    function      exists(fname)
        !
        implicit none
        !
        character(max_argument_length) :: fname
        logical :: exists
        !
        inquire(file=trim(fname),exist=exists)
        !
    end function  exists

    pure function strsplit(str,delimiter) result(word)
        !
        character(500), intent(in) :: str
        character(1), intent(in) :: delimiter
        character(len=:), allocatable :: word(:) ! character(500) :: word(500)
        integer :: i, j, k
        !
        j = 0; k = 0
        do i = 1, len(str)
            if (i .gt. k) then
            if ( str(i:i) .ne. delimiter ) then
                ! determine the next position of delimiter
                k = i+index(str(i:),delimiter)-1
                ! increase word count
                j=j+1
            endif
            endif
        enddo
        !
        allocate(character(500) :: word(j))
        !
        j = 0; k = 0
        do i = 1, len(str)
            if (i .gt. k) then
            if ( str(i:i) .ne. delimiter ) then
                ! determine the next position of delimiter
                k = i+index(str(i:),delimiter)-1
                ! increase word count
                j=j+1
                word(j) = str(i:k)
                ! write(*,"(I,A)") j, word(j)
            endif
            endif
        enddo
        !
    end function  strsplit
    

    end module








