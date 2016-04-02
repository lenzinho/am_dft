        character(*), intent(in) :: name
        character(*), intent(in), optional :: flags
        character(*), intent(in), optional :: filename
        character(*), intent(in), optional :: permission
        character(len=5) :: emph
        integer :: fid
        !
        ! defaults
        !
        emph = ' ... '
        if (present(flags)) then
            if (index(flags,'W')) emph = ' --> '
            if (index(flags,'E')) emph = ' >>> '
            if (index(flags,'N')) emph = '     '
        endif
        !
        if (present(filename)) then
            if (present(permission)) then
                fid = fopen(filename=filename,permission=permission)
            else
                fid = fopen(filename=filename,permission='w')
            endif
        else
            fid = 6
        endif
        !