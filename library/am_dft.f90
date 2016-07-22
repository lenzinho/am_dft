module am_dft

    use dispmodule
    use am_options
	use am_vasp_io
    use am_constants
    use am_histogram
    use am_dispersion
    use am_dos

    private

    type, public :: am_class_dft
    	type(am_class_dispersion) :: dr
		type(am_class_dos) 		  :: dos
		contains
			procedure :: load
    end type am_class_dft

contains

    subroutine     load(dft,opts,flags)
        ! flags = dispersion: fermi/(eigenval,procar)
        implicit none
        !
        class(am_class_dft)   ,intent(out) :: dft
        type(am_class_options), intent(in) :: opts
        character(*)          , intent(in) :: flags
        real(dp) :: Ef
        !
        if (opts%verbosity.ge.1) call print_title('DFT')
        !
		if     (index(flags,'dispersion').ne.0) then
			! load dispersion relations
	        if     (index(flags,'eigenval').ne.0) then
	        	! load dispersion form eigenval
	            write(*,'(a,a,a)') flare, 'file = ', trim(opts%eigenval)
	            call read_eigenval(E=dft%dr%E, nspins=dft%dr%nspins, nbands=dft%dr%nbands, iopt_filename=opts%eigenval, iopt_verbosity=0) ! , lmproj=dft%dr%lmproj)
	        elseif (index(flags,'procar').ne.0) then
	        	! load dispersion from procar
	            write(*,'(a,a,a)') flare, 'file = ', trim(opts%procar)
	            stop ' ERROR [load]: PROCAR NOT YET SUPPORTED.'
	            !     call read_procar(E=dft%dr%E, iopt_filename=opts%procar, iopt_verbosity=0, &
	            !         nbands=dft%dr%nbands, nions=dft%dr%nions, nspins=dft%dr%nspins, norbitals=dft%dr%norbitals, orbitals=dft%dr%orbitals, lmproj=dft%dr%lmproj)
	        else
	            stop 'ERROR [load]: eigenval or procar flag required'
	        endif
	        ! check if flags is asking to set fermi level at zero
	        if     (index(flags,'fermi').ne.0) then
	        	call read_doscar(efermi=Ef, nedos=dft%dos%nEs, dos=dft%dos%D, E=dft%dos%E, iopt_filename=opts%doscar, iopt_verbosity=0)
	        	! shift energies
	        	dft%dos%E = dft%dos%E - Ef
	        	dft%dr%E  = dft%dr%E  - Ef
	        endif
	        ! print stdout
	        if (opts%verbosity.ge.1) then
	            write(*,'(a,a,a)') flare, 'bands = ', tostring(dft%dr%nbands)
	            write(*,'(a,a,a)') flare, 'spins = ', tostring(dft%dr%nspins)
	            write(*,'(a,a)  ') flare, 'band energies = '
	            write(*,'(5x,a,a)')      'lowest = ', tostring(minval(pack(dft%dr%E,.true.)))
	            write(*,'(5x,a,a)')     'highest = ', tostring(maxval(pack(dft%dr%E,.true.)))
	            write(*,'(5x,a,a)')        'mean = ', tostring(sum(pack(dft%dr%E,.true.))/size(pack(dft%dr%E,.true.)))
	            if (index(flags,'fermi')) write(*,'(5x,a)') 'Note, energies have been shifted, putting Ef at zero. Originally, Ef = '//tostring(Ef)
	            ! make nice dos histogram
	            write(*,'(a,a)') flare, 'density of states'
	            call plot_histogram(binned_data=histogram(x=pack(dft%dr%E,.true.),m=98))
	        endif
	        !
	        if (dft%dr%nbands.eq.0) stop 'ERROR: nbands = 0'

	    endif
        !
    end subroutine load



end module am_dft