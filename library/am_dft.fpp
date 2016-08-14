#:include "fypp_macros.fpp"
module am_dft

    use dispmodule
    use am_constants
    use am_options
	use am_vasp_io
    use am_histogram
    use am_unit_cell
    use am_brillouin_zone
    use am_dispersion
    use am_dos

    private

    type, public :: am_class_dft
    	type(am_class_unit_cell)  :: uc
		type(am_class_bz) 		  :: bz
    	type(am_class_dispersion) :: dr
		type(am_class_dos) 		  :: dos
		contains
			procedure :: load => load_dft
    end type am_class_dft

contains

	!

    subroutine     load_dft(dft,opts,flags)
        ! flags = dispersion: fermi/(eigenval,procar)
        implicit none
        !
        class(am_class_dft)   ,intent(out) :: dft
        type(am_class_options), intent(in) :: opts
        character(*)          , intent(in) :: flags
        real(dp), allocatable :: kpt_recp(:,:) ! kpoint coordinates
        real(dp), allocatable :: w(:)     ! normalized weights
        real(dp) :: Ef
        real(dp) :: offset
        !
        if (opts%verbosity.ge.1) call print_title('DFT')
        !
        offset = 0.0_dp
        !
		if     (index(flags,'dispersion').ne.0) then
	        ! load poscar used in DFT
	        call dft%uc%load_poscar(opts)
			! load dispersion relations
	        if     (index(flags,'eigenval').ne.0) then
	        	! load dispersion form eigenval
	            if (opts%verbosity.ge.1) write(*,'(a,a,a)') flare, 'file = ', trim(opts%eigenval)
	            ! load band energies
	            call read_eigenval(E=dft%dr%E, nbands=dft%dr%nbands, iopt_filename=opts%eigenval, iopt_verbosity=0) ! , lmproj=dft%dr%lmproj)
	            ! load kpoints from eigenval
	            call read_eigenval(kpt=kpt_recp, w=w, iopt_filename=opts%eigenval, iopt_verbosity=0)
	            ! check weights
	            if ((sum(w)-1.0_dp).gt.tiny) stop 'ERROR [load_dft]: eigenval kpoint weights does not sum to one'
	            ! create instance
	            call dft%bz%create_bz(kpt_recp=kpt_recp, bas=dft%uc%bas, w=w, prec=opts%prec)
	            !
	        elseif (index(flags,'procar').ne.0) then
	        	! load dispersion from procar
	            if (opts%verbosity.ge.1) write(*,'(a,a,a)') flare, 'file = ', trim(opts%procar)
	            stop ' ERROR [load]: PROCAR NOT YET SUPPORTED.'
	            !     call read_procar(E=dft%dr%E, iopt_filename=opts%procar, iopt_verbosity=0, &
	            !         nbands=dft%dr%nbands, nions=dft%dr%nions, nspins=dft%dr%nspins, norbitals=dft%dr%norbitals, orbitals=dft%dr%orbitals, lmproj=dft%dr%lmproj)
	            ! check weights
	            if ((sum(w)-1.0_dp).gt.tiny) stop 'ERROR [load_dft]: procar kpoint weights does not sum to one'
	            ! load kpoints
	            call read_procar(kpt=kpt_recp, w=w, iopt_filename=opts%procar, iopt_verbosity=0)
	            ! create instance
	            call dft%bz%create_bz(kpt_recp=kpt_recp, bas=dft%uc%bas, w=w, prec=opts%prec)
	            !
	        else
	            stop 'ERROR [load]: eigenval or procar flag required'
	        endif
	        ! check for flags asking for shift in energy
	        if     (index(flags,'minimum').ne.0) then
	        	! shift lowest band energy to zero; E(nbands,nkpts)
	        	offset = minval(dft%dr%E((opts%skip+1):,:))
				dft%dr%E  = dft%dr%E  - offset
        	elseif (index(flags,'fermi').ne.0) then
        		! shift fermi level to zero
	        	call read_doscar(efermi=Ef, nedos=dft%dos%nEs, dos=dft%dos%D, E=dft%dos%E, iopt_filename=opts%doscar, iopt_verbosity=0)
	        	offset = Ef
	        	! shift energies
	        	dft%dos%E = dft%dos%E - offset
	        	dft%dr%E  = dft%dr%E  - offset
	        endif
	        ! stdout
	        if (opts%verbosity.ge.1) then
	        	write(*,'(a,a,a)') flare, 'bands = ', tostring(dft%dr%nbands)
	            write(*,'(a,a)  ') flare, 'band energies = '
	            write(*,'(5x,a,a)')      'lowest = ', tostring(minval(pack(dft%dr%E,.true.)))
	            write(*,'(5x,a,a)')     'highest = ', tostring(maxval(pack(dft%dr%E,.true.)))
	            write(*,'(5x,a,a)')        'mean = ', tostring(sum(pack(dft%dr%E,.true.))/size(pack(dft%dr%E,.true.)))
	            write(*,'(a,a)') flare, 'Energies shift = '//tostring(offset)
	            if (index(flags,'fermi'))   write(*,'(a,a)') flare, 'Fermi energy = '//tostring(Ef)
	            ! make nice dos histogram
	            write(*,'(a,a)') flare, 'density of states'
	            call plot_histogram(binned_data=histogram(x=pack(dft%dr%E,.true.),m=98))
	        endif
	        !
	        if (dft%dr%nbands.eq.0) stop 'ERROR: nbands = 0' 
	    else
	    	stop 'ERROR [dft:load]: unknown flag'
	    endif
        !
    end subroutine load_dft

end module am_dft








