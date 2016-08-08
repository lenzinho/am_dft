#:include "fypp_macros.fpp"
module am_options
    !
    use am_constants 
    use am_matlab
    !
    private
    !
    type, public :: am_class_options
        ! generic
        integer  :: verbosity
        real(dp) :: prec
        character(max_argument_length) :: flags
        ! uc
        integer :: supercell(3)
        integer :: deformation_code
        real(dp):: maxstrain
        integer :: nstrains
        ! tbvsk
        real(dp):: pair_cutoff
        ! tbfit
        integer :: skip
        !
        integer  :: n(3) ! monkhorst pack mesh grid dimensions
        real(dp) :: s(3) ! monkhorst pack mesh grid shift
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
        contains
        procedure :: set_default
        procedure :: parse_command_line_sym
        procedure :: parse_command_line_tbfit
        procedure :: parse_command_line_tbvsk
        procedure :: parse_command_line_uc
    end type am_class_options

contains

    subroutine     set_default(opts)
        !
        implicit none
        !
        class(am_class_options), intent(inout) :: opts
        !
        opts%verbosity        = 1
        opts%prec             = 1.0D-5           ! this is ultimately limited by the input
        opts%flags            = ''
        ! uc
        opts%supercell        = [1,1,1] ! no supercell
        opts%deformation_code = -1000
        opts%maxstrain        = -1000
        opts%nstrains         = -1000
        ! tbvsk
        opts%pair_cutoff      = 0.0_dp
        ! tbfit
        opts%skip             = 0
        !
        opts%n                = [9,9,9]          ! monkhorst pack mesh grid dimensions
        opts%s                = real([0,0,0],dp) ! monkhorst pack mesh grid shift
        opts%ibzkpt           = 'IBZKPT'
        opts%procar           = 'PROCAR'
        opts%prjcar           = 'PRJCAR'
        opts%doscar           = 'DOSCAR'
        opts%eigenval         = 'EIGENVAL'
        opts%poscar           = 'POSCAR'
        opts%wan_amn          = 'wannier90.amn'
        opts%wan_mmn          = 'wannier90.mmn'
        opts%wan_win          = 'wannier90.win'
        opts%wan_eig          = 'wannier90.eig'
    end subroutine set_default

    subroutine     parse_command_line_uc(opts)
        !
        implicit none
        !
        class(am_class_options), intent(inout) :: opts
        character(max_argument_length) :: argument
        integer :: narg
        integer :: i, j
        integer :: iostat
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
                    write(*,'(5x,a)') '-h'
                    write(*,'(5x,a)') '     Prints this message and exits.'
                    write(*,'(5x,a)') ''
                    write(*,'(5x,a)') '-verbosity <int>'
                    write(*,'(5x,a)') '     Controls stdout: (0) supress all output, (1) default, (2) verbose.'
                    write(*,'(5x,a)') ''
                    write(*,'(5x,a)') '-prec <dbl>'
                    write(*,'(5x,a)') '     Numerical precision.'
                    write(*,'(5x,a)') ''
                    write(*,'(5x,a)') '-primitive'
                    write(*,'(5x,a)') '     Produces outfile.primitive'
                    write(*,'(5x,a)') ''
                    write(*,'(5x,a)') '-conventional'
                    write(*,'(5x,a)') '     Produces outfile.conventional'
                    write(*,'(5x,a)') ''
                    write(*,'(5x,a)') '-supercell <a:int> <b:int> <c:int>'
                    write(*,'(5x,a)') '     Produces outfile.supercell containing cell expanded [a,b,c] times along fractional directions.'
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
                    stop
                case('-verbosity')
                    i=i+1
                    call get_command_argument(i,argument)
                    read(argument,*,iostat=iostat) opts%verbosity
                    if (iostat.ne.0) stop 'ERROR [parse_command_line_uc]: iostat /= 0'
                case('-prec')
                    i=i+1
                    call get_command_argument(i,argument)
                    read(argument,*,iostat=iostat) opts%prec
                    if (iostat.ne.0) stop 'ERROR [parse_command_line_uc]: iostat /= 0'
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
                        read(argument,*,iostat=iostat) opts%supercell(j)
                        if (iostat.ne.0) stop 'ERROR [parse_command_line_uc]: iostat /= 0'
                    enddo
                case('-strain')
                    opts%flags = trim(opts%flags)//' def'
                    i=i+1
                    call get_command_argument(i,argument)
                    read(argument,*,iostat=iostat) opts%deformation_code
                    if (iostat.ne.0) stop 'ERROR [parse_command_line_uc]: iostat /= 0'
                    i=i+1
                    call get_command_argument(i,argument)
                    read(argument,*,iostat=iostat) opts%maxstrain
                    if (iostat.ne.0) stop 'ERROR [parse_command_line_uc]: iostat /= 0'
                    i=i+1
                    call get_command_argument(i,argument)
                    read(argument,*,iostat=iostat) opts%nstrains
                    if (iostat.ne.0) stop 'ERROR [parse_command_line_uc]: iostat /= 0'
                case default
                    stop 'ERROR [parse_command_line_uc]: unrecognized option'
            end select
        enddo
    end subroutine parse_command_line_uc

    subroutine     parse_command_line_sym(opts)
        !
        implicit none
        !
        class(am_class_options), intent(inout) :: opts
        character(max_argument_length) :: argument
        integer :: narg
        integer :: i
        integer :: iostat
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
                    write(*,'(5x,a)') '-h'
                    write(*,'(5x,a)') '     Prints this message and exits.'
                    write(*,'(5x,a)') ''
                    write(*,'(5x,a)') '-verbosity <int>'
                    write(*,'(5x,a)') '     Controls stdout: (0) supress all output, (1) default, (2) verbose.'
                    write(*,'(5x,a)') ''
                    write(*,'(5x,a)') '-prec <dbl>'
                    write(*,'(5x,a)') '     Numerical precision.'
                    write(*,'(5x,a)') ''
                    stop
                case('-verbosity')
                    i=i+1
                    call get_command_argument(i,argument)
                    read(argument,*,iostat=iostat) opts%verbosity
                    if (iostat.ne.0) stop 'ERROR [parse_command_line_sym]: iostat /= 0'
                case('-prec')
                    i=i+1
                    call get_command_argument(i,argument)
                    read(argument,*,iostat=iostat) opts%prec
                    if (iostat.ne.0) stop 'ERROR [parse_command_line_sym]: iostat /= 0'
                case default
                    stop 'ERROR [parse_command_line_sym]: unrecognized option'
            end select
        enddo
    end subroutine parse_command_line_sym

    subroutine     parse_command_line_tbfit(opts)
        !
        implicit none
        !
        class(am_class_options), intent(inout) :: opts
        character(max_argument_length) :: argument
        integer :: narg
        integer :: i
        integer :: iostat
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
                    write(*,'(5x,a)') '-h'
                    write(*,'(5x,a)') '     Prints this message and exits.'
                    write(*,'(5x,a)') ''
                    write(*,'(5x,a)') '-verbosity <int>'
                    write(*,'(5x,a)') '     Controls stdout: (0) supress all output, (1) default, (2) verbose.'
                    write(*,'(5x,a)') ''
                    write(*,'(5x,a)') '-prec <dbl>'
                    write(*,'(5x,a)') '     Numerical precision.'
                    write(*,'(5x,a)') ''
                    write(*,'(5x,a)') '-file <eigenval/procar>'
                    write(*,'(5x,a)') '     File containing dispersion to fit.'
                    write(*,'(5x,a)') ''
                    write(*,'(5x,a)') '-shift <fermi/minimum/none>'
                    write(*,'(5x,a)') '     Shifts energies of DFT bands used in tight binding fit.'
                    write(*,'(5x,a)') '     fermi   - sets EF = 0, requires DOSCAR and tetrahedron integration, preferably.'
                    write(*,'(5x,a)') '     minimum - sets lowest energy of first tight binding band at zero.'
                    write(*,'(5x,a)') '     none    - (default) does nothing.'
                    write(*,'(5x,a)') ''
                    write(*,'(5x,a)') '-skip <int>'
                    write(*,'(5x,a)') '     Number of bands skipped by fit.'
                    write(*,'(5x,a)') ''
                    stop
                case('-verbosity')
                    i=i+1
                    call get_command_argument(i,argument)
                    read(argument,*,iostat=iostat) opts%verbosity
                    if (iostat.ne.0) stop 'ERROR [parse_command_line_tbfit]: iostat /= 0'
                case('-prec')
                    i=i+1
                    call get_command_argument(i,argument)
                    read(argument,*,iostat=iostat) opts%prec
                    if (iostat.ne.0) stop 'ERROR [parse_command_line_tbfit]: iostat /= 0'
                case('-file')
                    i=i+1
                    call get_command_argument(i,argument)
                    if     (index(argument,'eigenval')) then
                        opts%flags=trim(opts%flags)//'eigenval'
                    elseif (index(argument,'procar')) then
                        opts%flags=trim(opts%flags)//'procar'
                    else
                        stop 'ERROR [parse_command_line_tbfit]: file /= eigenval or procar'
                    endif
                case('-shift')
                    i=i+1
                    call get_command_argument(i,argument)
                    if     (index(argument,'fermi')) then
                        opts%flags=trim(opts%flags)//'fermi'
                    elseif (index(argument,'minimum')) then
                        opts%flags=trim(opts%flags)//'minimum'
                    elseif (index(argument,'none')) then
                        opts%flags=trim(opts%flags)//'none'
                    else
                        stop 'ERROR [parse_command_line_tbfit]: shift /= fermi, minimum, or none'
                    endif
                case('-skip')
                    i=i+1
                    call get_command_argument(i,argument)
                    read(argument,*,iostat=iostat) opts%skip
                    if (iostat.ne.0) stop 'ERROR [parse_command_line_tbfit]: iostat /= 0'
                    opts%flags=trim(opts%flags)//'skip'
                case default
                    stop 'ERROR [parse_command_line_tbfit]: unrecognized option'
            end select
        enddo
        !
        if (index(opts%flags,'skip').ne.0) then
            stop 'ERROR [parse_command_line_tbfit]: skip must be specified! See tbfit -h.'
        endif
        !
    end subroutine parse_command_line_tbfit

    subroutine     parse_command_line_tbvsk(opts)
        !
        implicit none
        !
        class(am_class_options), intent(inout) :: opts
        character(max_argument_length) :: argument
        integer :: narg
        integer :: i
        integer :: iostat
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
                    write(*,'(5x,a)') '-h'
                    write(*,'(5x,a)') '     Prints this message and exits.'
                    write(*,'(5x,a)') ''
                    write(*,'(5x,a)') '-verbosity <int>'
                    write(*,'(5x,a)') '     Controls stdout: (0) supress all output, (1) default, (2) verbose.'
                    write(*,'(5x,a)') ''
                    write(*,'(5x,a)') '-prec <dbl>'
                    write(*,'(5x,a)') '     Numerical precision.'
                    write(*,'(5x,a)') ''
                    write(*,'(5x,a)') '-pair_cutoff <int>'
                    write(*,'(5x,a)') '     Real-space cutoff for two-center tight-binding integrals.'
                    write(*,'(5x,a)') ''
                    stop
                case('-verbosity')
                    i=i+1
                    call get_command_argument(i,argument)
                    read(argument,*,iostat=iostat) opts%verbosity
                    if (iostat.ne.0) stop 'ERROR [parse_command_line_uc]: iostat /= 0'
                case('-prec')
                    i=i+1
                    call get_command_argument(i,argument)
                    read(argument,*,iostat=iostat) opts%prec
                    if (iostat.ne.0) stop 'ERROR [parse_command_line_uc]: iostat /= 0'
                case('-pair_cutoff')
                    i=i+1
                    call get_command_argument(i,argument)
                    read(argument,*,iostat=iostat) opts%pair_cutoff
                    if (iostat.ne.0) stop 'ERROR [parse_command_line_uc]: iostat /= 0'
                case default
                    stop 'ERROR [parse_command_line_tbvsk]: unrecognized option'
            end select
        enddo
    end subroutine parse_command_line_tbvsk

end module


















