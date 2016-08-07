module am_options
    !
    use am_constants 
    use am_matlab
    !
    private
    !
    type, public :: am_class_options
        integer  :: verbosity
        real(dp) :: prec
        integer  :: n(3) ! monkhorst pack mesh grid dimensions
        real(dp) :: s(3) ! monkhorst pack mesh grid shift
        character(max_argument_length) :: flags
        character(max_argument_length) :: ibzkpt
        character(max_argument_length) :: procar
        character(max_argument_length) :: prjcar
        character(max_argument_length) :: eigenval
        character(max_argument_length) :: doscar
        character(max_argument_length) :: poscar
        character(max_argument_length) :: wan_amn
        character(max_argument_length) :: wan_mmn
        character(max_argument_length) :: wan_win
        character(max_argument_length) :: wan_eig
        integer :: supercell(3)
        integer :: deformation_code
        real(dp):: maxstrain
        integer :: nstrains
        ! TB
        real(dp):: pair_cutoff
        integer :: skip_band
        contains
        procedure :: set_default
        procedure :: parse_control_file
        procedure :: parse_command_line_uc
    end type am_class_options

    type     am_class_file
        integer :: nlines
        character(:), allocatable :: line(:)
        contains
        procedure :: load
        procedure :: get_keyword
    end type am_class_file

contains

    subroutine     set_default(opts)
        !
        implicit none
        !
        class(am_class_options), intent(inout) :: opts
        !
        opts%verbosity = 1
        opts%prec     = 1.0D-5           ! this is ultimately limited by the input
        opts%n        = [9,9,9]          ! monkhorst pack mesh grid dimensions
        opts%s        = real([0,0,0],dp) ! monkhorst pack mesh grid shift
        opts%flags    = ''
        opts%ibzkpt   = 'IBZKPT'
        opts%procar   = 'PROCAR'
        opts%prjcar   = 'PRJCAR'
        opts%doscar   = 'DOSCAR'
        opts%eigenval = 'EIGENVAL'
        opts%poscar   = 'POSCAR'
        opts%wan_amn  = 'wannier90.amn'
        opts%wan_mmn  = 'wannier90.mmn'
        opts%wan_win  = 'wannier90.win'
        opts%wan_eig  = 'wannier90.eig'
        opts%supercell= [1,1,1] ! no supercell
        opts%deformation_code = -1000
        opts%maxstrain        = -1000
        opts%nstrains         = -1000
        opts%pair_cutoff      =  1000
        opts%skip_band        = 0
        !
    end subroutine set_default

    ! control file

    subroutine     load(file,fname)
        !
        implicit none
        !
        class(am_class_file), intent(out) :: file
        character(*), intent(in) :: fname
        character(maximum_buffer_size) :: buffer ! read buffer
        character :: tab
        integer :: pos
        integer :: nlines
        integer :: fid, istat
        integer :: i
        ! 
        tab = char(9)
        ! initialize fid
        fid = 1
        ! get number of lines
        open(unit=fid,file=trim(fname),status='old',action='read',iostat=istat)
            if (istat.ne. 0) then 
                call template_infile_control()
                stop 'ERROR [param_in_file]: problem reading infile.control. Creating template.control'
            endif
            nlines=0
            file_loop : do
                ! read line
                read(unit=fid,fmt='(a)',iostat=istat) buffer
                if (istat.eq.-1) exit file_loop ! end of file
                if (istat.ne. 0) stop 'ERROR [param_in_file]: problem reading buffer'
                ! increase line count
                nlines=nlines+1
            enddo file_loop
        close(fid)
        ! allocate space
        allocate(character(maximum_buffer_size)::file%line(nlines))
        ! initialize as empty
        file%line = ''
        ! read lines
        open(unit=fid,file='infile.control',status='old',action='read')
        do i = 1, nlines
            ! read line
            read(unit=fid,fmt='(a)',iostat=istat) buffer
            ! convert tabs into spaces
            pos = index(buffer,tab)
            do while (pos.ne.0)
               buffer(pos:pos) = ' '
               pos = index(buffer,tab)
            enddo
            ! check for comment
            buffer=adjustl(buffer)
            if (.not.(buffer(1:1).eq.'!'.or.buffer(1:1).eq.'#')) then
                ! check that the line has something
                if(len_trim(buffer).gt.0) then
                    file%line(i) = buffer
                endif
            endif
            !
        enddo
        ! trim
        file%line = trim_null(file%line)
        file%nlines = size(file%line)
        contains
        subroutine     template_infile_control()
            !
            implicit none
            !
            integer :: fid
            !
            fid = 1
            open(unit=fid,file='template.control',status='replace',action='write')
                write(unit=fid,fmt='(a)') 'verbosity   = 1'
                write(unit=fid,fmt='(a)') 'prec        = 1E-8'
                write(unit=fid,fmt='(a)') 'pair_cutoff = 3.2'
                write(unit=fid,fmt='(a)') 'fit_tb      = eigenval'
                write(unit=fid,fmt='(a)') 'skip_bands  = 5'
                write(unit=fid,fmt='(a)') '# vasp outputs'
                write(unit=fid,fmt='(a)') 'ibzkpt      = ./vasp/IBZKPT'
                write(unit=fid,fmt='(a)') 'doscar      = ./vasp/DOSCAR'
                write(unit=fid,fmt='(a)') 'prjcar      = ./vasp/PRJCAR'
                write(unit=fid,fmt='(a)') 'poscar      = ./vasp/POSCAR'
                write(unit=fid,fmt='(a)') 'procar      = ./vasp/PROCAR'
                write(unit=fid,fmt='(a)') 'eigenval    = ./vasp/EIGENVAL'
            close(fid)
        end subroutine template_infile_control
    end subroutine load

    subroutine     get_keyword(file,keyword,isfound,c_value,l_value,i_value,r_value)
        ! note: this routine mangles file object, copy file if it needs to be read later
        implicit none
        !
        class(am_class_file),  intent(in) :: file
        character(*), intent(in)  :: keyword
        logical     , intent(out) :: isfound
        character(*),optional, intent(inout) :: c_value
        logical     ,optional, intent(inout) :: l_value
        integer     ,optional, intent(inout) :: i_value
        real(dp)    ,optional, intent(inout) :: r_value
        integer :: kl, in,loop,itmp
        character(len=maximum_buffer_size) :: buffer
        !
        kl=len_trim(keyword)

        isfound=.false.
        ! this loop produces buffer at the end
        do loop = 1, file%nlines
           in=index(file%line(loop),trim(keyword))
           if (in==0 .or. in>1 ) cycle
           itmp=in+len(trim(keyword))
           if (file%line(loop)(itmp:itmp)/='=' &
                .and. file%line(loop)(itmp:itmp)/=':' &
                .and. file%line(loop)(itmp:itmp)/=' ') cycle
           if (isfound) then
              stop 'ERROR [parse_control]: a keyword occurs more than once'
           endif
           isfound=.true.
           buffer=file%line(loop)(kl+1:)
           ! file%line(loop)(:) = ' '
           buffer=adjustl(buffer)
           if ( buffer(1:1)=='=' .or. buffer(1:1)==':') then
              buffer=buffer(2:)
              buffer=adjustl(buffer)
           endif
        enddo
        ! buffer is now parsed
        if (isfound) then
           if ( present(c_value) ) c_value=buffer
           if ( present(l_value) ) then
              if (index(buffer,'t') > 0) then
                 l_value=.true.
              elseif (index(buffer,'f') > 0) then
                 l_value=.false.
              else
                 stop 'ERROR [parse_control]: problem reading logical keyword'
              endif
           endif
           if ( present(i_value) ) read(buffer,*) i_value
           if ( present(r_value) ) read(buffer,*) r_value
        end if
        !
    end subroutine get_keyword

    ! program-specific parsing

    subroutine      parse_command_line_uc(opts)
        !
        implicit none
        !
        class(am_class_options), intent(inout) :: opts
        character(max_argument_length) :: argument
        integer :: narg
        integer :: i, j
        !
        call opts%set_default
        !
        narg=command_argument_count()
        !
        i=0
        do while ( i < narg )
            i=i+1
            call get_command_argument(i,argument)
            select case(argument)
                case('-h')
                    call help_command_line_uc
                    stop
                case('-verbosity')
                    i=i+1
                    call get_command_argument(i,argument)
                    read(argument,*) opts%verbosity
                case('-prec')
                    i=i+1
                    call get_command_argument(i,argument)
                    read(argument,*) opts%prec
                case('-poscar')
                    i=i+1
                    call get_command_argument(i,opts%poscar)
                    opts%poscar = trim_null(opts%poscar)
                case('-conventional') 
                    opts%flags = trim(opts%flags)//' conv'
                case('-primitive') 
                    opts%flags = trim(opts%flags)//' prim'
                case('-supercell')
                    opts%flags = trim(opts%flags)//' supercell'
                    do j = 1,3
                        i=i+1
                        call get_command_argument(i,argument)
                        read(argument,*) opts%supercell(j)
                    enddo
                case('-strain')
                    opts%flags = trim(opts%flags)//' def'
                    i=i+1
                    call get_command_argument(i,argument)
                    read(argument,*) opts%deformation_code
                    i=i+1
                    call get_command_argument(i,argument)
                    read(argument,*) opts%maxstrain
                    i=i+1
                    call get_command_argument(i,argument)
                    read(argument,*) opts%nstrains
                case default
                    write(0,*) trim(argument),' is not a recognized option, these are:'
                    call help_command_line_uc
                    stop
            end select
        enddo
        contains
        subroutine      help_command_line_uc
            write(*,'(5x,a)') '-h'
            write(*,'(5x,a)') '     Print this message and exit.'
            write(*,'(5x,a)') ''
            write(*,'(5x,a)') '-verbosity <int>'
            write(*,'(5x,a)') '     Controls stdout: (0) supress all output, (1) default, (2) verbose.'
            write(*,'(5x,a)') ''
            write(*,'(5x,a)') '-primitive'
            write(*,'(5x,a)') '     Converts cell to primitive cell.'
            write(*,'(5x,a)') ''
            write(*,'(5x,a)') '-conventional'
            write(*,'(5x,a)') '     Converts cell to conventional setting.'
            write(*,'(5x,a)') ''
            write(*,'(5x,a)') '-supercell <a:int> <b:int> <c:int>'
            write(*,'(5x,a)') '     Creates supercell of dimesnions a,b,c based on input structure.'
            write(*,'(5x,a)') ''
            write(*,'(5x,a)') '-strain <defcode:int> <eta:dbl> <strains:int>'
            write(*,'(5x,a)') '     Applies elastic deformation to cell. Deformation code species type of deformation. eta'
            write(*,'(5x,a)') '     `specifies the size of the deformation. Strains specifies how many states to consider.'
            write(*,'(5x,a)') ''
            write(*,'(5x,a)') '     `defcode  strain                                description'
            write(*,'(5x,a)') '     `-------  ------------------------------------  --------------------------------------'
            write(*,'(5x,a)') '     `    0    [ eta,  eta,  eta,    0,    0,    0]  volume strain '
            write(*,'(5x,a)') '     `    1    [ eta,    0,    0,    0,    0,    0]  linear strain along x '
            write(*,'(5x,a)') '     `    2    [   0,  eta,    0,    0,    0,    0]  linear strain along y '
            write(*,'(5x,a)') '     `    3    [   0,    0,  eta,    0,    0,    0]  linear strain along z '
            write(*,'(5x,a)') '     `    4    [   0,    0,    0,  eta,    0,    0]  yz shear strain'
            write(*,'(5x,a)') '     `    5    [   0,    0,    0,    0,  eta,    0]  xz shear strain'
            write(*,'(5x,a)') '     `    6    [   0,    0,    0,    0,    0,  eta]  xy shear strain 2C44'
            write(*,'(5x,a)') '     `    7    [   0,    0,    0,  eta,  eta,  eta]  shear strain along (111)'
            write(*,'(5x,a)') '     `    8    [ eta,  eta,    0,    0,    0,    0]  xy in-plane strain '
            write(*,'(5x,a)') '     `    9    [ eta, -eta,    0,    0,    0,    0]  xy in-plane shear strain'
            write(*,'(5x,a)') '     `   10    [ eta,  eta,  eta,  eta,  eta,  eta]  global strain'
            write(*,'(5x,a)') '     `   11    [  1e,   2e,   3e,   4e,   5e,   6e]  universal strain 1'
            write(*,'(5x,a)') '     `   12    [ -2e,   1e,   4e,  -3e,   6e,  -5e]  universal strain 2'
            write(*,'(5x,a)') '     `   13    [  3e,  -5e,  -1e,   6e,   2e,  -4e]  universal strain 3'
            write(*,'(5x,a)') '     `   14    [ -4e,  -6e,   5e,   1e,  -3e,   2e]  universal strain 4'
            write(*,'(5x,a)') '     `   15    [  5e,   4e,   6e,  -2e,  -1e,  -3e]  universal strain 5'
            write(*,'(5x,a)') '     `   16    [ -6e,   3e,  -2e,   5e,  -4e,   1e]  universal strain 6'
            write(*,'(5x,a)') '     `   20    [   0,    0,    0,  eta,  eta,    0]  yz zx shear'
            write(*,'(5x,a)') '         40    [ eta,  eta,  -2e,    0,    0,    0]  tetragonal for 3/2(c11-c12)'
            write(*,'(5x,a)') '         41    [ eta, -eta,    0,    0,    0,    0]  orthorhombic for (c11-c12)'
            write(*,'(5x,a)') ''
        end subroutine help_command_line_uc
    end subroutine  parse_command_line_uc

    ! used tb

    subroutine      parse_control_file(opts)
        !
        implicit none
        !
        class(am_class_options), intent(out) :: opts
        type(am_class_file) :: file
        character(maximum_buffer_size) :: buffer
        logical :: isfound
        !
        ! load control file
        call file%load('infile.control')
        ! set defaults
        call opts%set_default
        ! global
        call file%get_keyword(keyword='verbosity'       ,isfound=isfound, i_value=opts%verbosity)
        call file%get_keyword(keyword='prec'            ,isfound=isfound, r_value=opts%prec)
        call file%get_keyword(keyword='restart'         ,isfound=isfound, c_value=buffer)
        if (isfound) opts%flags = trim(opts%flags)//'restart'
        ! vasp
        call file%get_keyword(keyword='ibzkpt'          ,isfound=isfound, c_value=opts%ibzkpt)
        call file%get_keyword(keyword='doscar'          ,isfound=isfound, c_value=opts%doscar)
        call file%get_keyword(keyword='prjcar'          ,isfound=isfound, c_value=opts%prjcar)
        call file%get_keyword(keyword='poscar'          ,isfound=isfound, c_value=opts%poscar)
        call file%get_keyword(keyword='procar'          ,isfound=isfound, c_value=opts%procar)
        call file%get_keyword(keyword='eigenval'        ,isfound=isfound, c_value=opts%eigenval)
        ! wannier
        call file%get_keyword(keyword='wan_amn'         ,isfound=isfound, c_value=opts%wan_amn)
        call file%get_keyword(keyword='wan_mmn'         ,isfound=isfound, c_value=opts%wan_mmn)
        call file%get_keyword(keyword='wan_win'         ,isfound=isfound, c_value=opts%wan_win)
        call file%get_keyword(keyword='wan_eig'         ,isfound=isfound, c_value=opts%wan_eig)
        ! tb
        call file%get_keyword(keyword='pair_cutoff'     ,isfound=isfound, r_value=opts%pair_cutoff)
        call file%get_keyword(keyword='fit_tb'          ,isfound=isfound, c_value=buffer)
        if     (isfound) then
            opts%flags = trim(opts%flags)//' fit_tb:'
            if     (index(buffer,'eigenval')) then
                opts%flags=trim(opts%flags)//'eigenval'
            elseif (index(buffer,'procar')) then
                opts%flags=trim(opts%flags)//'procar'
            else
                stop 'ERROR [parse_control_file]: fit_tb /= eigenval or procar'
            endif
            call file%get_keyword(keyword='skip_band'   ,isfound=isfound, i_value=opts%skip_band)
            if (.not.isfound) then
                stop 'ERROR [parse_control_file]: skip_band is required for fitting'
            endif
        endif
    end subroutine  parse_control_file


end module
