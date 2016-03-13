module am_options
    !
    use am_constants 
    !
    type am_class_options
        integer  :: verbosity
        real(dp) :: sym_prec
        integer  :: n(3) ! monkhorst pack mesh grid dimensions
        real(dp) :: s(3) ! monkhorst pack mesh grid shift
        character(max_argument_length) :: ibzkpt
        character(max_argument_length) :: procar
        character(max_argument_length) :: prjcar
        character(max_argument_length) :: eigenval
        character(max_argument_length) :: poscar
        character(max_argument_length) :: wan_amn
        character(max_argument_length) :: wan_mmn
        character(max_argument_length) :: wan_win
        character(max_argument_length) :: wan_eig

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
        opts%sym_prec = tiny ! seems to work.
        opts%n        = [9,9,9] ! monkhorst pack mesh grid dimensions
        opts%s        = real([0,0,0],dp) ! monkhorst pack mesh grid shift
        opts%ibzkpt   = 'IBZKPT'
        opts%procar   = 'PROCAR'
        opts%prjcar   = 'PRJCAR'
        opts%eigenval = 'EIGENVAL'
        opts%poscar   = 'POSCAR'
        opts%wan_amn  = 'wannier90.amn'
        opts%wan_mmn  = 'wannier90.mmn'
        opts%wan_win  = 'wannier90.win'
        opts%wan_eig  = 'wannier90.eig'
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
                case('-sym_prec')
                    i=i+1
                    call get_command_argument(i,argument)
                    read(argument,*) opts%sym_prec
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
                ! vasp
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
        write(*,*) ' -h'
        write(*,*) '        Print this message and exit.'
        write(*,*) ''
    end subroutine  am_printhelp

end module
