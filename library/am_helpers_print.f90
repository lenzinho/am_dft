    subroutine     print_sparse_double(title,A,iopt_emph,iopt_fid)
        !
        implicit none
        !
        character(len=*), intent(in) :: title
        real(dp), intent(in) :: A(:,:)
        character(len=5), optional, intent(in) :: iopt_emph
        character(len=5) :: emph
        integer, optional, intent(in) :: iopt_fid
        integer :: i, j
        integer :: fid
        !
        if ( present(iopt_emph) ) then
            emph = iopt_emph
        else
            emph = "     "
        endif
        if ( present(iopt_fid)  ) then
            fid = iopt_fid
        else
            fid = 6
        endif
        !
        write(fid,'(a,a," = ")') emph, trim(title)
        write(fid,'(5x,a,a,a)') ,'+', repeat('-',size(A,2)), '+'
        !
        do i = 1, size(A,1)
            write(fid,'(5x)',advance='no')
            write(fid,'(a1)',advance='no') '|'
            do j = 1, size(A,2)
                if(abs(A(i,j)).gt.tiny) then
                    write(*,'(a1)',advance='no') 'x'
                else
                    write(*,'(a1)',advance='no') ' '
                endif
            enddo
            write(fid,'(a1)',advance='no') '|'
            write(fid,*)
        enddo
        !
        write(fid,'(5x,a,a,a)') ,'+', repeat('-',size(A,2)), '+'
        !
    end subroutine print_sparse_double

    subroutine     print_sparse_integer(title,A,iopt_emph,iopt_fid)
        !
        implicit none
        !
        character(len=*), intent(in) :: title
        integer, intent(in) :: A(:,:)
        character(len=5), optional, intent(in) :: iopt_emph
        character(len=5) :: emph
        integer, optional, intent(in) :: iopt_fid
        integer :: i, j
        integer :: fid
        !
        if ( present(iopt_emph) ) then
            emph = iopt_emph
        else
            emph = "     "
        endif
        if ( present(iopt_fid)  ) then
            fid = iopt_fid
        else
            fid = 6
        endif
        !
        write(fid,'(a,a," = ")') emph, trim(title)
        write(fid,'(5x,a,a,a)') ,'+', repeat('-',size(A,2)), '+'
        !
        do i = 1, size(A,1)
            write(fid,'(5x)',advance='no')
            write(fid,'(a1)',advance='no') '|'
            do j = 1, size(A,2)
                if(A(i,j).ne.0) then
                    write(*,'(a1)',advance='no') 'x'
                else
                    write(*,'(a1)',advance='no') ' '
                endif
            enddo
            write(fid,'(a1)',advance='no') '|'
            write(fid,*)
        enddo
        !
        write(fid,'(5x,a,a,a)') ,'+', repeat('-',size(A,2)), '+'
        !
    end subroutine print_sparse_integer

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
        write(fid,'(a,a,a)') "+", repeat("-",column_width-2), "+"
        write(fid,'(a,a,a)') "|",  centertitle(title,column_width-2), "|" 
        write(fid,'(a,a,a)') "+", repeat("-",column_width-2), "+"
        !
    end subroutine am_print_title

    subroutine     am_print_logical( name, bool , iopt_emph, iopt_fid )
        !
        implicit none
        !
        character(len=*), intent(in) :: name
        logical, intent(in) :: bool
        character(len=5), optional, intent(in) :: iopt_emph
        character(len=5) :: emph
        integer, optional, intent(in) :: iopt_fid
        integer :: fid
        !
        fid = 6
        if ( present(iopt_fid)  ) fid = iopt_fid
        emph = "     "
        if ( present(iopt_emph) ) emph = iopt_emph
        !
        write(fid,'(a,a," = ",l)') emph, name, bool
    end subroutine am_print_logical

    subroutine     am_print_double( name, dbl , iopt_emph, iopt_fid )
        !
        implicit none
        !
        character(len=*), intent(in) :: name
        real(dp) :: dbl
        character(len=5), optional, intent(in) :: iopt_emph
        character(len=5) :: emph
        integer, optional, intent(in) :: iopt_fid
        integer :: fid
        !
        fid = 6
        if ( present(iopt_fid)  ) fid = iopt_fid
        emph = "     "
        if ( present(iopt_emph) ) emph = iopt_emph
        !
        write(fid,'(a,a," = ",g15.5)') emph, name, dbl2char(dbl)
    end subroutine am_print_double

    subroutine     am_print_integer( name, i , iopt_emph, iopt_fid )
        !
        implicit none
        !
        character(len=*), intent(in) :: name
        integer, intent(in) :: i
        character(len=5), optional, intent(in) :: iopt_emph
        character(len=5) :: emph
        integer, optional, intent(in) :: iopt_fid
        integer :: fid
        !
        fid = 6
        if ( present(iopt_fid)  ) fid = iopt_fid
        emph = "     "
        if ( present(iopt_emph) ) emph = iopt_emph
        !
        write(fid,'(a,a," = ",a)') emph, name, int2char(i)
    end subroutine am_print_integer

    subroutine     am_print_string( name, str , iopt_emph, iopt_fid )
        !
        implicit none
        !
        character(len=*), intent(in) :: name, str
        character(len=5), optional, intent(in) :: iopt_emph
        character(len=5) :: emph
        integer, optional, intent(in) :: iopt_fid
        integer :: fid
        !
        fid = 6
        if ( present(iopt_fid)  ) fid = iopt_fid
        emph = "     "
        if ( present(iopt_emph) ) emph = iopt_emph
        !
        write(fid,'(a,a," = ",a)') emph, name, str
    end subroutine am_print_string

    subroutine     am_print_logical_matrix( name, A , iopt_emph, iopt_fid )
        !
        implicit none
        !
        character(len=*), intent(in) :: name
        logical, intent(in) :: A(:,:)
        integer :: i,j,m,n
        character(len=5), optional, intent(in) :: iopt_emph
        character(len=5) :: emph
        integer, optional, intent(in) :: iopt_fid
        integer :: fid
        !
        fid = 6
        if ( present(iopt_fid)  ) fid = iopt_fid
        emph = "     "
        if ( present(iopt_emph) ) emph = iopt_emph
        !
        m=size(A,1)
        n=size(A,2)
        write(*,'(a,a," = ",a)') emph, name
        do i = 1, m
            write(fid,"(5x,300l)") ( A(i,j), j = 1, n )
        enddo
    end subroutine am_print_logical_matrix

    subroutine     am_print_double_matrix(iopt_name, A, iopt_emph, iopt_fname)
        !
        implicit none
        !
        character(*), intent(in), optional :: iopt_name
        character(5), intent(in), optional :: iopt_emph
        character(*), intent(in), optional :: iopt_fname
        real(dp), intent(in) :: A(:,:)
        integer :: i, j
        integer :: fid
        !
        if (present(iopt_fname)) then
            fid = 1
            open(unit=fid,file=trim(iopt_fname),status='replace',action='write')
        else
            fid=6
        endif
        !
        if (present(iopt_emph)) then
            write(unit=fid,fmt='(a)',advance='no') iopt_emph
        else
            write(unit=fid,fmt='(a)',advance='no') '     '
        endif
        if (present(iopt_name)) then
            write(unit=fid,fmt='(a)',advance='no') iopt_name
            write(unit=fid,fmt='(a)',advance='no') ' = '
        endif
        write(unit=fid,fmt='(a)')
        !
        do i = 1, size(A,1)
            write(unit=fid,fmt='(5x)',advance='no')
            do j = 1, size(A,2)
            write(unit=fid,fmt='(f15.5)',advance='no') A(i,j)
            enddo
            write(unit=fid,fmt=*)
        enddo
        !
        if (present(iopt_fname)) close(fid)
        !
    end subroutine am_print_double_matrix

    subroutine     am_print_complex_matrix( name, A , iopt_emph, iopt_fid )
        !
        implicit none
        !
        character(len=*), intent(in) :: name
        complex(dp), intent(in) :: A(:,:)
        integer :: i,j,m,n
        character(len=5), optional, intent(in) :: iopt_emph
        character(len=5) :: emph
        integer, optional, intent(in) :: iopt_fid
        integer :: fid
        !
        fid=6
        if ( present(iopt_fid)  ) fid = iopt_fid
        emph = "     "
        if ( present(iopt_emph) ) emph = iopt_emph
        !
        m=size(A,1)
        n=size(A,2)
        write(fid,'(a,a," = ",a)') emph, name
        do i = 1, m
            write(fid,"(5x,100f15.5)") ( A(i,j), j = 1, n )
        enddo
    end subroutine am_print_complex_matrix

    subroutine     am_print_int_matrix( name, A , iopt_emph, iopt_fid )
        !
        implicit none
        !
        character(len=*), intent(in) :: name
        integer, intent(in) :: A(:,:)
        integer :: i,j,m,n
        character(len=5), optional, intent(in) :: iopt_emph
        character(len=5) :: emph
        integer, optional, intent(in) :: iopt_fid
        integer :: fid
        !
        fid = 6
        if ( present(iopt_fid)  ) fid = iopt_fid
        emph = "     "
        if ( present(iopt_emph) ) emph = iopt_emph
        !
        m=size(A,1)
        n=size(A,2)
        write(fid,'(a,a," = ",a)') emph, name
        do i = 1, m
            write(fid,"(5x,300i4)") ( A(i,j), j = 1, n )
        enddo
    end subroutine am_print_int_matrix

    subroutine     am_print_double_vector( name, A , iopt_emph, iopt_fid )
        !
        implicit none
        !
        character(len=*), intent(in) :: name
        real(dp), intent(in) :: A(:)
        integer :: i,m
        character(len=5), optional, intent(in) :: iopt_emph
        character(len=5) :: emph
        integer, optional, intent(in) :: iopt_fid
        integer :: fid
        !
        fid = 6
        if ( present(iopt_fid)  ) fid = iopt_fid
        emph = "     "
        if ( present(iopt_emph) ) emph = iopt_emph
        !
        m=size(A)
        write(fid,'(a,a," = ",a)') emph, name
        write(fid,"(5x,100f15.5)") ( A(i), i = 1, m )
    end subroutine am_print_double_vector

    subroutine     am_print_int_vector( name, A , iopt_emph, iopt_fid )
        !
        implicit none
        !
        character(len=*), intent(in) :: name
        integer, intent(in) :: A(:)
        integer :: j,m
        character(len=5), optional, intent(in) :: iopt_emph
        character(len=5) :: emph
        integer, optional, intent(in) :: iopt_fid
        integer :: fid
        !
        fid = 6
        if ( present(iopt_fid)  ) fid = iopt_fid
        emph = "     "
        if ( present(iopt_emph) ) emph = iopt_emph
        !
        m=size(A)
        write(fid,'(a,a," = ",a)') emph, name
        write(fid,"(5x,300i4)") ( A(j), j = 1, m )
    end subroutine am_print_int_vector