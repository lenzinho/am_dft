! #:include "fypp_macros.fpp"
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

    type, public :: am_class_wavecar
    	integer :: n(3) ! total divisions along each primitive direction
    	integer :: n_fft(3) ! utilized for creation of fft mesh x = [-n(1):n(1)], etc, followed by meshgrid(z,y,x)
    	integer :: nrpts
    	integer :: nkpts
    	integer :: nbands
    	integer :: nspins
    	real(dp)   , allocatable :: rpt_frac(:,:) ! rpt(3,nrpts)
        complex(sp), allocatable :: psi(:,:,:,:) ! psi(nrpts,kpt,band,spin)
    end type am_class_wavecar

    type, public, extends(am_class_dft) :: am_class_vasp
		type(am_class_wavecar)    :: wc
	end type am_class_vasp

contains

    subroutine     load_dft(dft,opts,flags)
        ! flags = dispersion: fermi/(eigenval,procar,wavecar)
        ! 		  wavefunction
       	!  		  structure 
        implicit none
        !
        class(am_class_dft)   ,intent(out) :: dft
        type(am_class_options), intent(in) :: opts
        character(*)          , intent(in) :: flags
        real(dp), allocatable :: kpt_recp(:,:) ! kpoint coordinates
        real(dp), allocatable :: w(:)     ! normalized weights
        real(dp) :: Ef
        real(dp) :: offset
        integer  :: i,j,k
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
	            call read_eigenval(E=dft%dr%E, nbands=dft%dr%nbands, fname=opts%eigenval, verbosity=0) ! , lmproj=dft%dr%lmproj)
	            ! load kpoints from eigenval
	            call read_eigenval(kpt=kpt_recp, w=w, fname=opts%eigenval, verbosity=0)
	            ! check weights
	            if ((sum(w)-1.0_dp).gt.tiny) stop 'ERROR [load_dft]: eigenval kpoint weights does not sum to one'
	            ! create instance
	            call dft%bz%create_bz(kpt_recp=kpt_recp, bas=dft%uc%bas, prec=opts%prec)
	            !
	        elseif (index(flags,'procar').ne.0) then
	        	! load dispersion from procar
	            if (opts%verbosity.ge.1) write(*,'(a,a,a)') flare, 'file = ', trim(opts%procar)
	            stop ' ERROR [load]: PROCAR NOT YET SUPPORTED.'
	            !     call read_procar(E=dft%dr%E, fname=opts%procar, verbosity=0, &
	            !         nbands=dft%dr%nbands, nions=dft%dr%nions, nspins=dft%dr%nspins, norbitals=dft%dr%norbitals, orbitals=dft%dr%orbitals, lmproj=dft%dr%lmproj)
	            ! check weights
	            if ((sum(w)-1.0_dp).gt.tiny) stop 'ERROR [load_dft]: procar kpoint weights does not sum to one'
	            ! load kpoints
	            call read_procar(kpt=kpt_recp, w=w, fname=opts%procar, verbosity=0)
	            ! create instance
	            call dft%bz%create_bz(kpt_recp=kpt_recp, bas=dft%uc%bas, prec=opts%prec)
	            !
            elseif  (index(flags,'wavecar').ne.0) then
                ! make sure it is a vasp structure
                select type (dft)
                type is (am_class_vasp)
	            	! query wavecar
			        call query_wavecar(nkpts=dft%bz%nkpts, nbands=dft%dr%nbands, fname=opts%wavecar, verbosity=0)
			        ! allocate kpoint space
			        allocate(kpt_recp(3,dft%bz%nkpts))
			        ! read kpoints
			        do i = 1, dft%bz%nkpts
						call read_wavecar(kpt_id=i, band_id=1, spin_id=1, kpt=kpt_recp(1:3,i), fname=opts%wavecar, verbosity=0)
			        enddo
		            ! create bz instance
		            ! NOTE: BAS USED FOR BZ CREATION IS READ FROM POSCAR RATHER THAN FROM WAVECAR
		            call dft%bz%create_bz(kpt_recp=kpt_recp, bas=dft%uc%bas, prec=opts%prec)
			        ! allocate space for dispersion
			        allocate(dft%dr%E(dft%dr%nbands,dft%bz%nkpts))
			        ! read dispersion
			        do i = 1, dft%bz%nkpts
						call read_wavecar(kpt_id=i, band_id=1, spin_id=1, E=dft%dr%E(:,i), fname=opts%wavecar, verbosity=0)
			        enddo
                class default
                    ! do nothing
                    stop 'ERROR [load]: type must be vasp in order to read wavecar'
                end select
	        else
	            stop 'ERROR [load]: eigenval, procar, or wavecar flag required'
	        endif
	        ! check for flags asking for shift in energy
	        if     (index(flags,'shift:minimum').ne.0) then
	        	! shift lowest band energy to zero; E(nbands,nkpts)
	        	offset = minval(dft%dr%E((opts%skip+1):,:))
				dft%dr%E  = dft%dr%E  - offset
        	elseif (index(flags,'shift:fermi').ne.0) then
        		! shift fermi level to zero
	        	call read_doscar(efermi=Ef, nedos=dft%dos%nEs, dos=dft%dos%D, E=dft%dos%E, fname=opts%doscar, verbosity=0)
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
	    ! <\LOAD DISPERSION>
        elseif (index(flags,'wavefunction').ne.0) then
		! <LOAD WAVEFUNCTION>
            ! make sure it is a vasp structure
            select type (dft)
            type is (am_class_vasp)
            	! DOES NOT LOAD ALL WAVEFUNCTIONS. ONLY THE ONES ASKED FOR BY OPTS
		        ! transfer options
		        dft%wc%nkpts  = size(opts%kpt_id)
		        dft%wc%nbands = size(opts%band_id)
		        dft%wc%nspins = size(opts%spin_id)
		        ! query workspace
		        call query_wavecar(n_fft=dft%wc%n_fft,fname=opts%wavecar,verbosity=opts%verbosity)
				! get mesh dimensions
				dft%wc%n = dft%wc%n_fft*2+1
		        ! number of real space points
		        dft%wc%nrpts = product(dft%wc%n)
			    ! generate real-space grid in fractional coordinates between [0,1)
			    dft%wc%rpt_frac = meshgrid( v1=[0:2*dft%wc%n_fft(1)]/real(1+2*dft%wc%n_fft(1),dp), &
									      & v2=[0:2*dft%wc%n_fft(2)]/real(1+2*dft%wc%n_fft(2),dp), &
									      & v3=[0:2*dft%wc%n_fft(3)]/real(1+2*dft%wc%n_fft(3),dp))
		        ! allocate flattened psi wavefunction
		        allocate(dft%wc%psi(dft%wc%nrpts, dft%wc%nkpts, dft%wc%nbands, dft%wc%nspins))
		        ! read band
		        do i = 1, dft%wc%nkpts
		    	do j = 1, dft%wc%nbands
				do k = 1, dft%wc%nspins
					call read_wavecar(rpt_frac=dft%wc%rpt_frac, kpt_id=opts%kpt_id(i), band_id=opts%band_id(j), &
						spin_id=opts%spin_id(k), psi=dft%wc%psi(:,i,j,k), fname=opts%wavecar, verbosity=0)
		        enddo
		        enddo
		        enddo
            class default
                ! do nothing
                stop 'ERROR [load]: type must be vasp in order to read wavecar'
            end select
		! <\LOAD WAVEFUNCTION>
        elseif (index(flags,'structure').ne.0) then
		! <LOAD POSCAR>
	        ! load poscar used in DFT
	        call dft%uc%load_poscar(opts)
        ! <\LOAD POSCAR>
	    else
	    	stop 'ERROR [dft:load]: unknown flag'
	    endif
        !
    end subroutine load_dft

end module am_dft








