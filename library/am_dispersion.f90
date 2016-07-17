module am_dispersion

    use dispmodule
    use am_brillouin_zone
    use am_unit_cell
    use am_constants
    use am_stdout
    use am_vasp_io
    use am_mkl
    use am_options
    use am_histogram
    use am_shells
    use am_symmetry_rep
    use am_tight_binding

    implicit none

    private

    type, public :: am_class_dispersion
        integer :: nbands
        integer :: nspins
        real(dp), allocatable :: E(:,:) ! E(nbands,nkpts) energies
    contains
        procedure :: debug_dump => debug_dump_dr
    end type am_class_dispersion

    type, public, extends(am_class_dispersion) :: am_class_dispersion_dft
        contains
        procedure :: load
    end type am_class_dispersion_dft

    type, public, extends(am_class_dispersion) :: am_class_dispersion_tb
        complex(dp), allocatable :: C(:,:,:) ! tight binding coefficients
        contains
        procedure :: get_dispersion
    end type am_class_dispersion_tb

contains

    subroutine     debug_dump_dr(dr,fname)
        !
        implicit none
        !
        class(am_class_dispersion), intent(in) :: dr
        character(*), intent(in) :: fname
        !
                                                call dump(A=dr%nbands            ,fname=trim(fname)//'.nbands'           )
        if (allocated(dr%E                   )) call dump(A=dr%E                 ,fname=trim(fname)//'.E'                )
        ! ! dump dr specific data                                                                                        
        ! select type (dr)                                                                                         
        ! class is (am_class_dispersion)                                                                                    
        ! class is (am_class_seitz_group)                                                                           
        ! class default
        !     stop 'ERROR [debug_dump]: invalid group class'
        ! end select

        !
    end subroutine debug_dump_dr

    subroutine     load(dr,opts,flags)
        !
        implicit none
        !
        class(am_class_dispersion_dft), intent(out) :: dr
        type(am_class_options)        , intent(in)  :: opts
        character(*)                  , intent(in)  :: flags
        !
        if (opts%verbosity.ge.1) call print_title('Electronic band dispersion (DFT)')
        !
        if     (index(flags,'eigenval')) then
            write(*,'(a,a,a)') flare, 'file = ', 'eigenval'
            call read_eigenval(E=dr%E, nspins=dr%nspins, nbands=dr%nbands, iopt_verbosity=0, iopt_filename=opts%eigenval) ! , lmproj=dr%lmproj)
        elseif (index(flags,'procar')) then
            write(*,'(a,a,a)') flare, 'file = ', 'procar'
            write(*,*) 'PROCAR NOT YET SUPPORTED.'
            stop
            !     call read_procar(E=dr%E, iopt_filename=opts%procar, iopt_verbosity=0, &
            !         nbands=dr%nbands, nions=dr%nions, nspins=dr%nspins, norbitals=dr%norbitals, orbitals=dr%orbitals, lmproj=dr%lmproj)
        else
            stop 'ERROR: eigenval/procar flag required.'
        endif
        !
        if (opts%verbosity.ge.1) then
            write(*,'(a,a,a)') flare, 'bands = ', tostring(dr%nbands)
            write(*,'(a,a,a)') flare, 'spins = ', tostring(dr%nspins)
            write(*,'(a,a)  ') flare, 'band energy statistics = '
            write(*,'(5x,a,a)')     'lowest = ',  tostring(minval(pack(dr%E,.true.)))
            write(*,'(5x,a,a)')     'highest = ', tostring(maxval(pack(dr%E,.true.)))
            write(*,'(5x,a,a)')     'mean = ',    tostring(sum(pack(dr%E,.true.))/size(pack(dr%E,.true.)))
            ! make nice dos histogram
            write(*,'(a,a)') flare, 'density of states'
            call plot_histogram(binned_data=histogram(x=pack(dr%E,.true.),m=98))
        endif
        !
        if (dr%nbands.eq.0) stop 'ERROR: nbands = 0'
        !
    end subroutine load

    subroutine     get_dispersion(dr,tb,tbpg,bz,pp)
        !
        implicit none
        !
        class(am_class_dispersion_tb),intent(out) :: dr
        type(am_class_tb_group)      , intent(in) :: tbpg ! point group in tight binding representation
        type(am_class_tight_binding) , intent(in) :: tb
        type(am_class_bz)            , intent(in) :: bz
        type(am_class_prim_pair)     , intent(in) :: pp
        complex(dp), allocatable :: H(:,:), V(:,:)
        real(dp)   , allocatable :: D(:)
        integer :: i
        !
        ! tight binding coefficients
        allocate(dr%C(tbpg%nbases,tbpg%nbases,bz%nkpts))
        ! eigenvalues
        allocate(dr%E(tbpg%nbases,bz%nkpts))
        ! loop over kpoints
        do i = 1, bz%nkpts
            ! construct hamiltonian
            H = tb%get_hamiltonian(tbpg=tbpg, pp=pp, kpt=bz%kpt_cart(:,i))
            ! diagonalize hamiltonian
            call am_zheev(A=H,V=V,D=D)
            dr%C(:,:,i) = V
            dr%E(:,i) = D
        enddo
        !
    end subroutine get_dispersion

!     subroutine     write_dispersion(tb)
!         !
!         implicit none
!         !
! !         ! write it in fractional coordinates so that  

! !     end subroutine



end module am_dispersion