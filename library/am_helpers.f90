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
    interface am_print_sparse
        module procedure am_print_sparse_double, am_print_sparse_integer
    end interface ! am_print_sparse
    !
    interface linspace
        module procedure linspace_double, linspace_integer
    end interface ! linspace
    !
    interface issubset
        module procedure issubset_vec, issubset_mat
    end interface ! issubset
    !
    interface diag
        module procedure diag1, diag2
    end interface ! diag
    !
    contains

    !
    ! include files
    !

#include "am_helpers_print.inc"

#include "am_helpers_unique.inc"

#include "am_helpers_matlab.inc"

    !
    ! CONVERT DATA TYPES  
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
        else
            allocate(character(100) :: str)
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
        else
            allocate(character(100) :: str)
        endif
        !
        write(buffer,'(SP,f)') dbl
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
            write(fid,'(5x,a)',advance='no') centertitle(' ... ',n*15)
            write(fid,'(5x,a)',advance='no') centertitle(' ... ',n*15)
            write(fid,*)
        endif
        !
    end subroutine am_print_two_matrices_side_by_side

    !
    ! matrix functions
    !

    pure function  mesh_grid(n) result(grid_points)
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
    end function   mesh_grid

    !
    ! generic math functions 
    !

    pure function issubset_vec(group,element,iopt_sym_prec) result(issubset)
        !> returns true if true if element(:) is a subset of group(:,i) within numerical precision for any i
        implicit none
        !
        real(dp), intent(in) :: group(:,:)
        real(dp), intent(in) :: element(:)
        real(dp), intent(in), optional :: iopt_sym_prec
        real(dp) :: sym_prec
        logical :: issubset
        integer :: i
        !
        if (present(iopt_sym_prec)) then
            sym_prec = iopt_sym_prec
        else
            sym_prec = tiny
        endif
        !
        issubset = .false.
        !
        do i = 1,size(group,2)
            if ( all(abs(group(:,i)-element).lt.sym_prec) ) then
                issubset = .true.
                return
            endif
        enddo
        !
    end function  issubset_vec

    pure function issubset_mat(group,element,iopt_sym_prec) result(issubset)
        !> returns true if true if element(:,:) is a subset of group(:,:,i) within numerical precision for any i
        implicit none
        !
        real(dp), intent(in) :: group(:,:,:)
        real(dp), intent(in) :: element(:,:)
        real(dp), intent(in), optional :: iopt_sym_prec
        real(dp) :: sym_prec
        integer :: i
        logical :: issubset
        !
        if (present(iopt_sym_prec)) then
            sym_prec = iopt_sym_prec
        else
            sym_prec = tiny
        endif
        !
        issubset = .false.
        !
        do i = 1,size(group,3)
            if ( all(abs(group(:,:,i)-element(:,:)).lt.sym_prec) ) then
                issubset = .true.
                return
            endif
        enddo
        !
    end function  issubset_mat

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


    end module








