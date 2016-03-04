    module am_helpers
    !
    use am_constants
    !
    implicit none
    !
    interface unique
        module procedure unique_integer, unique_double
    end interface
    !
    interface am_print
        ! function overloading must start with the same name as the interface, i.e. "am_print" in this case
        module procedure am_print_double, am_print_integer, am_print_string, am_print_logical, &
            am_print_double_matrix, am_print_int_matrix, am_print_logical_matrix, &
            am_print_double_vector, am_print_int_vector
    end interface
    !
    contains


    !
    ! ERROR HANDELING
    !

    subroutine errore(file,line,calling_routine,message)
        !
        implicit none
        !
        character(len=*), intent(in), optional :: file
        integer, intent(in), optional :: line
        character(len=*), intent(in) :: calling_routine, message
        !
        write(*,'(/,x,79("%"))')
        write(*,'(5x,a10)') "ERROR! "
        write(*,'(5x,a10,a)') "routine: ", trim(calling_routine)
        if (present(file)) write(*,'(5x,a10,a)') "file: ", file
        if (present(line)) write(*,'(5x,a10,a)') "line: ", int2char(line)
        write(*,'(5x,a10,a)') "message: ", trim(message)
        write(*,'(1x,79("%"),/)')
        !
        stop
        !
        return
    end subroutine errore

    !
    ! PRINT TO STDOUT
    !

    pure function centertitle(title,length) result(title_centered)
        ! ----------------------------------------------------------------------+ 
        !      this subroutine centers a title by elimating 1/2 trailing " "s.  | 
        ! ----------------------------------------------------------------------+ 
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
    end function centertitle
    subroutine am_print_title(title)
        !
        implicit none
        !
        character(len=*), intent(in) :: title
        integer :: fid
        !
        fid = 6
        !
        write(fid,'(a,a,a)') "+", repeat("-",column_width), "+"
        write(fid,'(a,a,a)') "|",  centertitle(title,column_width), "|" 
        write(fid,'(a,a,a)') "+", repeat("-",column_width), "+"
        !
    end subroutine am_print_title
    subroutine am_print_logical( name, bool , in_emph, in_fid )
        !
        implicit none
        !
        character(len=*), intent(in) :: name
        logical, intent(in) :: bool
        character(len=5), optional, intent(in) :: in_emph
        character(len=5) :: emph
        integer, optional, intent(in) :: in_fid
        integer :: fid
        !
        fid = 6
        if ( present(in_fid)  ) fid = in_fid
        emph = "     "
        if ( present(in_emph) ) emph = in_emph
        !
        write(fid,'(a,a," = ",l)') emph, name, bool
    end subroutine am_print_logical
    subroutine am_print_double( name, dbl , in_emph, in_fid )
        !
        implicit none
        !
        character(len=*), intent(in) :: name
        real(dp) :: dbl
        character(len=5), optional, intent(in) :: in_emph
        character(len=5) :: emph
        integer, optional, intent(in) :: in_fid
        integer :: fid
        !
        fid = 6
        if ( present(in_fid)  ) fid = in_fid
        emph = "     "
        if ( present(in_emph) ) emph = in_emph
        !
        write(fid,'(a,a," = ",g15.5)') emph, name, dbl2char(dbl)
    end subroutine am_print_double
    subroutine am_print_integer( name, i , in_emph, in_fid )
        !
        implicit none
        !
        character(len=*), intent(in) :: name
        integer, intent(in) :: i
        character(len=5), optional, intent(in) :: in_emph
        character(len=5) :: emph
        integer, optional, intent(in) :: in_fid
        integer :: fid
        !
        fid = 6
        if ( present(in_fid)  ) fid = in_fid
        emph = "     "
        if ( present(in_emph) ) emph = in_emph
        !
        write(fid,'(a,a," = ",a)') emph, name, int2char(i)
    end subroutine am_print_integer
    subroutine am_print_string( name, str , in_emph, in_fid )
        !
        implicit none
        !
        character(len=*), intent(in) :: name, str
        character(len=5), optional, intent(in) :: in_emph
        character(len=5) :: emph
        integer, optional, intent(in) :: in_fid
        integer :: fid
        !
        fid = 6
        if ( present(in_fid)  ) fid = in_fid
        emph = "     "
        if ( present(in_emph) ) emph = in_emph
        !
        write(fid,'(a,a," = ",a)') emph, name, str
    end subroutine am_print_string
    subroutine am_print_logical_matrix( name, A , in_emph, in_fid )
        !
        implicit none
        !
        character(len=*), intent(in) :: name
        logical, intent(in) :: A(:,:)
        integer :: i,j,m,n
        character(len=5), optional, intent(in) :: in_emph
        character(len=5) :: emph
        integer, optional, intent(in) :: in_fid
        integer :: fid
        !
        fid = 6
        if ( present(in_fid)  ) fid = in_fid
        emph = "     "
        if ( present(in_emph) ) emph = in_emph
        !
        m=size(A,1)
        n=size(A,2)
        write(*,'(a,a," = ",a)') emph, name
        do i = 1, m
            write(fid,"(5x,300l)") ( A(i,j), j = 1, n )
        enddo
    end subroutine am_print_logical_matrix
    subroutine am_print_double_matrix( name, A , in_emph, in_fid )
        !
        implicit none
        !
        character(len=*), intent(in) :: name
        real(dp), intent(in) :: A(:,:)
        integer :: i,j,m,n
        character(len=5), optional, intent(in) :: in_emph
        character(len=5) :: emph
        integer, optional, intent(in) :: in_fid
        integer :: fid
        !
        fid=6
        if ( present(in_fid)  ) fid = in_fid
        emph = "     "
        if ( present(in_emph) ) emph = in_emph
        !
        m=size(A,1)
        n=size(A,2)
        write(fid,'(a,a," = ",a)') emph, name
        do i = 1, m
            write(fid,"(5x,100g15.5)") ( A(i,j), j = 1, n )
        enddo
    end subroutine am_print_double_matrix
    subroutine am_print_int_matrix( name, A , in_emph, in_fid )
        !
        implicit none
        !
        character(len=*), intent(in) :: name
        integer, intent(in) :: A(:,:)
        integer :: i,j,m,n
        character(len=5), optional, intent(in) :: in_emph
        character(len=5) :: emph
        integer, optional, intent(in) :: in_fid
        integer :: fid
        !
        fid = 6
        if ( present(in_fid)  ) fid = in_fid
        emph = "     "
        if ( present(in_emph) ) emph = in_emph
        !
        m=size(A,1)
        n=size(A,2)
        write(fid,'(a,a," = ",a)') emph, name
        do i = 1, m
            write(fid,"(5x,300i4)") ( A(i,j), j = 1, n )
        enddo
    end subroutine am_print_int_matrix
    subroutine am_print_double_vector( name, A , in_emph, in_fid )
        !
        implicit none
        !
        character(len=*), intent(in) :: name
        real(dp), intent(in) :: A(:)
        integer :: i,m
        character(len=5), optional, intent(in) :: in_emph
        character(len=5) :: emph
        integer, optional, intent(in) :: in_fid
        integer :: fid
        !
        fid = 6
        if ( present(in_fid)  ) fid = in_fid
        emph = "     "
        if ( present(in_emph) ) emph = in_emph
        !
        m=size(A)
        write(fid,'(a,a," = ",a)') emph, name
        write(fid,"(5x,100g15.5)") ( A(i), i = 1, m )
    end subroutine am_print_double_vector
    subroutine am_print_int_vector( name, A , in_emph, in_fid )
        !
        implicit none
        !
        character(len=*), intent(in) :: name
        integer, intent(in) :: A(:)
        integer :: j,m
        character(len=5), optional, intent(in) :: in_emph
        character(len=5) :: emph
        integer, optional, intent(in) :: in_fid
        integer :: fid
        !
        fid = 6
        if ( present(in_fid)  ) fid = in_fid
        emph = "     "
        if ( present(in_emph) ) emph = in_emph
        !
        m=size(A)
        write(fid,'(a,a," = ",a)') emph, name
        write(fid,"(5x,300i4)") ( A(j), j = 1, m )
    end subroutine am_print_int_vector

    !
    ! CONVERT DATA TYPES
    !
    character(100) function int2char(int) result(string)
    !
    implicit none
    !
    integer :: int
    !
    write(string,'(i)') int
    string = trim( adjustl( string ) )
    end function int2char
    character(100) function dbl2char(dbl) result(string)
    !
    implicit none
    !
    real(dp) :: dbl
    !
    write(string,'(f)') dbl
    string = trim( adjustl( string ) )
    end function dbl2char

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
    end function
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
        end function
    subroutine progressbar(message,j,totn)
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
    end function eye
    ! generic math functions
    pure function unique_double(A,iopt_tiny) result(B)
        !> returns unique columns of double matrix A(:,i) within numerical precision
        implicit none
        !
        real(dp), intent(in) :: A(:,:)
        real(dp) :: storage(size(A,1),size(A,2))
        real(dp), allocatable :: B(:,:)
        integer :: i1,i2,k
        real(dp), intent(in), optional :: iopt_tiny
        real(dp) :: tiny_prec
        !
        if ( present(iopt_tiny) ) then
            tiny_prec = iopt_tiny
        else
            tiny_prec = tiny
        endif
        !
        k=1
        storage(:,1) = A(:,1)
        !
        try_loop : do i1 = 1,size(A,2)
            do i2 = 1,k
                if ( all(abs(A(:,i1)-storage(:,i2)).lt.tiny_prec) ) cycle try_loop
            enddo
            k = k + 1
            storage(:,k) = A(:,i1)
        enddo try_loop
        !
        allocate(B(size(A,1),k))
        !
        do i1 = 1,k
            B(:,i1) = storage(:,i1)
        enddo
    end function unique_double
    pure function unique_integer(A) result(B)
        !> returns unique columns of double matrix A(:,i) within numerical precision
        implicit none
        !
        integer, intent(in) :: A(:,:)
        integer :: storage(size(A,1),size(A,2))
        integer, allocatable :: B(:,:)
        integer :: i1,i2,k
        !
        k=1
        storage(:,1) = A(:,1)
        !
        try_loop : do i1 = 1,size(A,2)
            do i2 = 1,k
                if ( all(abs(A(:,i1)-storage(:,i2)).eq.0) ) cycle try_loop
            enddo
            k = k + 1
            storage(:,k) = A(:,i1)
        enddo try_loop
        !
        allocate(B(size(A,1),k))
        !
        do i1 = 1,k
            B(:,i1) = storage(:,i1)
        enddo
    end function unique_integer
    pure function issubset(A,B)
        !> returns true if true if B(:) is a subset of A(:,i) within numerical precision for any i
        implicit none
        !
        real(dp), intent(in) :: A(:,:)
        real(dp), intent(in) :: B(:)
        logical :: issubset
        integer :: i
        !
        issubset = .false.
        !
        do i = 1,size(A,2)
            if ( all(abs(A(:,i)-B).lt.tiny) ) then
                issubset = .true.
                exit
            endif
        enddo
        !
    end function issubset
    pure function linspace(d1,d2,n) result(y)
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
        d = (d2-d1)/(n-1)
        do i = 1,n
            y(i) = d1 + (i-1.0_dp)*d;
        enddo
        !
    end function linspace
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
    END FUNCTION cross_product
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
    end function
    
    
    ! helper functions
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
    end function strsplit
    

    end module








