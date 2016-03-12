module am_wannier
	!
	use am_constants
	use am_helpers
	!
	implicit none
	
	type am_class_wannier
		integer :: nwanniers
		integer :: nkpts
		integer :: nbands
		complex(dp), allocatable :: U(:,:,:) ! u_matrix(nbands,nbands,nkpts)
		real(dp), allocatable :: E(:,:) ! E(nbands,nkpts)
	contains 
		procedure :: load_amn
		procedure :: load_eig
	end type am_class_wannier

contains

    subroutine     load_amn(wan,opts)
        ! 
        ! Reads wannier amn file.
        !
        use am_options
        use am_vasp_io
        use am_brillouin_zone
        !
        implicit none
        !
        class(am_class_wannier), intent(inout) :: wan
        type(am_class_options), intent(in) :: opts
        !
        call read_amn(U=wan%U,&
        	nbands=wan%nbands,&
        	nkpts=wan%nkpts,&
        	nwanniers=wan%nwanniers,&
        	iopt_filename=opts%wan_amn,&
        	iopt_verbosity=opts%verbosity)
        !
    end subroutine load_amn

    subroutine     load_eig(wan,opts)
        ! 
        ! Reads wannier amn file.
        !
        use am_options
        use am_vasp_io
        use am_brillouin_zone
        !
        implicit none
        !
        class(am_class_wannier), intent(inout) :: wan
        type(am_class_options), intent(in) :: opts
        !
        call read_eig(E=wan%E,&
        	nbands=wan%nbands,&
        	nkpts=wan%nkpts,&
        	iopt_filename=opts%wan_eig,&
        	iopt_verbosity=opts%verbosity)
        !
    end subroutine load_eig


end module