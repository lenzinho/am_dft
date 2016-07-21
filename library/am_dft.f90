module am_dft


    use am_dispersion

    private

    type, extends(am_class_dispersion) :: am_class_dispersion_dft
        contains
        procedure :: load
    end type am_class_dispersion_dft

    type, public :: am_class_dft
    	type(am_class_dispersion_dft) :: dr
		type(am_class_dispersion_dft) :: dos
	    contains
    end type am_class_dft


contains




            call dr_dft%read_doscar(natoms,efermi,nedos,dos,E,iopt_filename,iopt_verbosity)




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


end module am_dft