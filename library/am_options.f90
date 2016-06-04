module am_options
    !
    use am_constants 
    !
    type am_class_options
        integer  :: verbosity
        real(dp) :: prec
        integer  :: n(3) ! monkhorst pack mesh grid dimensions
        real(dp) :: s(3) ! monkhorst pack mesh grid shift
        character(max_argument_length) :: flags
        character(max_argument_length) :: tbf
        character(max_argument_length) :: ibzkpt
        character(max_argument_length) :: procar
        character(max_argument_length) :: prjcar
        character(max_argument_length) :: eigenval
        character(max_argument_length) :: poscar
        character(max_argument_length) :: wan_amn
        character(max_argument_length) :: wan_mmn
        character(max_argument_length) :: wan_win
        character(max_argument_length) :: wan_eig
        integer :: supercell(3)
        logical :: primitive
        logical :: conventional
        logical :: deformation
        integer :: deformation_code
        real(dp):: maxstrain
        integer :: nstrains
        real(dp):: pair_cutoff

    contains
        procedure :: defaults
        procedure :: get
    end type am_class_options    
    !
    contains

    pure subroutine defaults(opts)
        !
        implicit none
        !
        class(am_class_options), intent(inout) :: opts
        !
        opts%verbosity = 1
        opts%prec     = 1.0D-5  ! this is ultimately limited by the input
        opts%n        = [9,9,9] ! monkhorst pack mesh grid dimensions
        opts%s        = real([0,0,0],dp) ! monkhorst pack mesh grid shift
        opts%flags    = ''
        opts%ibzkpt   = 'IBZKPT'
        opts%procar   = 'PROCAR'
        opts%prjcar   = 'PRJCAR'
        opts%eigenval = 'EIGENVAL'
        opts%poscar   = 'POSCAR'
        opts%wan_amn  = 'wannier90.amn'
        opts%wan_mmn  = 'wannier90.mmn'
        opts%wan_win  = 'wannier90.win'
        opts%wan_eig  = 'wannier90.eig'
        opts%tbf      = 'infile.tightbinding'
        opts%supercell        = [1,1,1] ! no supercell
        opts%deformation_code = -1000
        opts%maxstrain        = -1000
        opts%nstrains         = -1000
        opts%pair_cutoff      =  1000
        !
    end subroutine  defaults

    subroutine      get(opts)
        !
        implicit none
        !
        class(am_class_options), intent(inout) :: opts
        character(max_argument_length) :: argument
        integer :: narg
        integer :: i, j
        !
        call opts%defaults
        !
        narg=command_argument_count()
        !
        i=0
        do while ( i < narg )
            i=i+1
            call get_command_argument(i,argument)
            select case(argument)
                case('-h')
                    call am_printhelp
                    stop
                case('-verbosity')
                    i=i+1
                    call get_command_argument(i,argument)
                    read(argument,*) opts%verbosity
                case('-baseline')
                    opts%flags = trim(opts%flags)//' baseline'
                case('-query')
                    opts%flags = trim(opts%flags)//' query'
                case('-prec')
                    i=i+1
                    call get_command_argument(i,argument)
                    read(argument,*) opts%prec
                ! ************************************************
                ! vasp input parameters
                ! 
                case('-ibzkpt')
                    i=i+1
                    call get_command_argument(i,opts%ibzkpt)
                    opts%ibzkpt = trim(opts%ibzkpt)
                case('-procar')
                    i=i+1
                    call get_command_argument(i,opts%procar)
                    opts%procar = trim(opts%procar)
                case('-prjcar')
                    i=i+1
                    call get_command_argument(i,opts%prjcar)
                    opts%prjcar = trim(opts%prjcar)
                case('-eigenval')
                    i=i+1
                    call get_command_argument(i,opts%eigenval)
                    opts%eigenval = trim(opts%eigenval)
                case('-poscar')
                    i=i+1
                    call get_command_argument(i,opts%poscar)
                    opts%poscar = trim(opts%poscar)
                ! ************************************************
                ! prog_tb parameters
                case('-pair_cutoff')
                    i=i+1
                    call get_command_argument(i,argument)
                    read(argument,*) opts%pair_cutoff
                case('-tbfit')
                    opts%flags = trim(opts%flags)//' tbfit:'
                    i=i+1
                    call get_command_argument(i,argument)
                    if     (index(argument,'eigenval').ne.0) then
                        opts%flags = trim(opts%flags)//'eigenval'
                    elseif (index(argument,'procar').ne.0) then
                        opts%flags = trim(opts%flags)//'procar'
                    else
                        stop '-tbfit /= eigenval or procar!'
                    endif
                ! ************************************************
                ! prog_uc parameters
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
                !
                ! monkhorst pack mesh
                !
                case('-s')
                    do j = 1,3
                        i=i+1
                        call get_command_argument(i,argument)
                        read(argument,*) opts%s(j)
                    enddo
                case('-n')
                    do j = 1,3
                        i=i+1
                        call get_command_argument(i,argument)
                        read(argument,*) opts%n(j)
                    enddo
                !
                ! wannier
                !
                case('-wan_amn')
                    i=i+1
                    call get_command_argument(i,opts%wan_amn)
                    opts%wan_amn=trim(opts%wan_amn)
                case('-wan_mmn')
                    i=i+1
                    call get_command_argument(i,opts%wan_mmn)
                    opts%wan_mmn=trim(opts%wan_mmn)
                case('-wan_win')
                    i=i+1
                    call get_command_argument(i,opts%wan_win)
                    opts%wan_win=trim(opts%wan_win)
                case('-wan_eig')
                    i=i+1
                    call get_command_argument(i,opts%wan_eig)
                    opts%wan_eig=trim(opts%wan_eig)
                case default
                    write(0,*) trim(argument),' is not a recognized option, these are:'
                    call am_printhelp
                    stop
            end select
        enddo
    end subroutine  get

    subroutine      am_printhelp
        write(*,'(5x,a)') '-h'
        write(*,'(5x,a)') '     Print this message and exit.'
        write(*,'(5x,a)') ''
        write(*,'(5x,a)') '-verbosity <int>'
        write(*,'(5x,a)') '     Controls output: (0) supress output, (1) normal, (2) verbose.'
        write(*,'(5x,a)') ''
        write(*,'(5x,a)') '-baseline'
        write(*,'(5x,a)') '     Generates baseline output for benchmarking.'
        write(*,'(5x,a)') ''
        write(*,'(5x,a)') '-query'
        write(*,'(5x,a)') '     Query structure.'
        write(*,'(5x,a)') ''
        write(*,'(5x,a)') '-prec <dbl>'
        write(*,'(5x,a)') '     Controls precision used for symmetry determination and calculation. Default = tiny.'
        write(*,'(5x,a)') '' ! PROG_TB STUFF BEGINS HERE
        write(*,'(5x,a)') '-pair_cutoff <dbl>'
        write(*,'(5x,a)') '     Realspace cutoff for second order force constants.'
        write(*,'(5x,a)') '' ! PROG_UC STUFF BEGINS HERE
        write(*,'(5x,a)') '-primitive'
        write(*,'(5x,a)') '     Converts cell to primitive cell.'
        write(*,'(5x,a)') ''
        write(*,'(5x,a)') '-conventional'
        write(*,'(5x,a)') '     Converts cell to conventional setting.'
        write(*,'(5x,a)') ''
        write(*,'(5x,a)') '-strain <defcode:int> <eta:dbl> <strains:int>'
        write(*,'(5x,a)') '     Applies elastic deformation to cell. Deformation code species type of deformation. eta'
        write(*,'(5x,a)') '     specifies the size of the deformation. Strains specifies how many states to consider.'
        write(*,'(5x,a)') ''
        write(*,'(5x,a)') 'defcode                                   strain                       description'
        write(*,'(5x,a)') ' ------  ---------------------------------------  --------------------------------'
        write(*,'(5x,a)') '      0     [ eta,  eta,  eta,    0,    0,    0]      volume strain '
        write(*,'(5x,a)') '      1     [ eta,    0,    0,    0,    0,    0]      linear strain along x '
        write(*,'(5x,a)') '      2     [   0,  eta,    0,    0,    0,    0]      linear strain along y '
        write(*,'(5x,a)') '      3     [   0,    0,  eta,    0,    0,    0]      linear strain along z '
        write(*,'(5x,a)') '      4     [   0,    0,    0,  eta,    0,    0]      yz shear strain'
        write(*,'(5x,a)') '      5     [   0,    0,    0,    0,  eta,    0]      xz shear strain'
        write(*,'(5x,a)') '      6     [   0,    0,    0,    0,    0,  eta]      xy shear strain 2C44'
        write(*,'(5x,a)') '      7     [   0,    0,    0,  eta,  eta,  eta]      shear strain along (111)'
        write(*,'(5x,a)') '      8     [ eta,  eta,    0,    0,    0,    0]      xy in-plane strain '
        write(*,'(5x,a)') '      9     [ eta, -eta,    0,    0,    0,    0]      xy in-plane shear strain'
        write(*,'(5x,a)') '     10     [ eta,  eta,  eta,  eta,  eta,  eta]      global strain'
        write(*,'(5x,a)') '     11     [  1e,   2e,   3e,   4e,   5e,   6e]      universal strain 1'
        write(*,'(5x,a)') '     12     [ -2e,   1e,   4e,  -3e,   6e,  -5e]      universal strain 2'
        write(*,'(5x,a)') '     13     [  3e,  -5e,  -1e,   6e,   2e,  -4e]      universal strain 3'
        write(*,'(5x,a)') '     14     [ -4e,  -6e,   5e,   1e,  -3e,   2e]      universal strain 4'
        write(*,'(5x,a)') '     15     [  5e,   4e,   6e,  -2e,  -1e,  -3e]      universal strain 5'
        write(*,'(5x,a)') '     16     [ -6e,   3e,  -2e,   5e,  -4e,   1e]      universal strain 6'
        write(*,'(5x,a)') '     20     [   0,    0,    0,  eta,  eta,    0]      yz zx shear'
        write(*,'(5x,a)') '     40     [ eta,  eta,  -2e,    0,    0,    0]      tetragonal for 3/2(c11-c12)'
        write(*,'(5x,a)') '     41     [ eta, -eta,    0,    0,    0,    0]      orthorhombic for (c11-c12)'
        write(*,'(5x,a)') '-supercell <a:int> <b:int> <c:int>'
        write(*,'(5x,a)') '     Creates supercell of dimesnions a,b,c based on input structure.'
        write(*,'(5x,a)') ''


!                 !
!                 ! monkhorst pack mesh
!                 !
!                 case('-s')
!                     do j = 1,3
!                         i=i+1
!                         call get_command_argument(i,argument)
!                         read(argument,*) opts%s(j)
!                     enddo
!                 case('-n')
!                     do j = 1,3
!                         i=i+1
!                         call get_command_argument(i,argument)
!                         read(argument,*) opts%n(j)
!                     enddo
!                 !
!                 ! vasp
!                 ! 
!                 case('-ibzkpt')
!                     i=i+1
!                     call get_command_argument(i,opts%ibzkpt)
!                     opts%ibzkpt = trim(opts%ibzkpt)
!                 case('-procar')
!                     i=i+1
!                     call get_command_argument(i,opts%procar)
!                     opts%procar = trim(opts%procar)
!                 case('-prjcar')
!                     i=i+1
!                     call get_command_argument(i,opts%prjcar)
!                     opts%prjcar = trim(opts%prjcar)
!                 case('-eigenval')
!                     i=i+1
!                     call get_command_argument(i,opts%eigenval)
!                     opts%eigenval = trim(opts%eigenval)
!                 case('-poscar')
!                     i=i+1
!                     call get_command_argument(i,opts%poscar)
!                     opts%poscar = trim(opts%poscar)
!                 !
!                 ! wannier
!                 !
!                 case('-wan_amn')
!                     i=i+1
!                     call get_command_argument(i,opts%wan_amn)
!                     opts%wan_amn=trim(opts%wan_amn)
!                 case('-wan_mmn')
!                     i=i+1
!                     call get_command_argument(i,opts%wan_mmn)
!                     opts%wan_mmn=trim(opts%wan_mmn)
!                 case('-wan_win')
!                     i=i+1
!                     call get_command_argument(i,opts%wan_win)
!                     opts%wan_win=trim(opts%wan_win)
!                 case('-wan_eig')
!                     i=i+1
!                     call get_command_argument(i,opts%wan_eig)
!                     opts%wan_eig=trim(opts%wan_eig)
!                 case default
!                     write(0,*) trim(argument),' is not a recognized option, these are:'
!                     call am_printhelp
!                     stop
    end subroutine  am_printhelp

end module
