module am_wannier
	!
	use am_unit_cell
	use am_constants
	use am_helpers
	!
	implicit none
	
	type am_class_wannier
		integer :: nwanniers
		integer :: nkpts
		integer :: nbands
		integer :: nrpts
		real(dp)   , allocatable :: rpt(:,:)   ! rpt(3,nrpts)
		complex(dp), allocatable :: U(:,:,:)   ! U(nbands,nbands,nkpts)
		real(dp)   , allocatable :: E(:,:)     ! E(nbands,nkpts)
		complex(dp), allocatable :: H_k(:,:,:) ! H(nbands,nbands,nkpts) reciprocal space hamiltonian
		complex(dp), allocatable :: H_r(:,:,:) ! H(nwanniers,nwanniers,nkpts) real space hamiltonian
	contains 
		procedure :: load_amn
		procedure :: load_eig
		!
		procedure :: get_hamiltonian
		procedure :: wannier_development_space
	end type am_class_wannier

contains

	!
	! create real space mesh
	!

	!	function


	!
    ! functions which operate on wannier object
    !

    subroutine     load_amn(wan,opts)
        ! 
        ! Reads wannier amn file.
        !
        use am_options
        use am_vasp_io
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

    subroutine     get_hamiltonian(wan)
    	!
    	class(am_class_wannier), intent(inout) :: wan
    	!
    	integer :: i
    	integer :: k
  		real(dp), allocatable :: E(:,:)
    	!
    	allocate(E(wan%nbands,wan%nbands))
    	!
		allocate(wan%H_k(wan%nbands,wan%nbands,wan%nkpts))
    	!
    	do k = 1, wan%nkpts
    		!
    		! hamiltonian in diagonalized bloch-state basis
    		!
	    	E=0.0_dp
	    	do i = 1,wan%nbands
		    	E(i,i) = wan%E(i,k)
	    	enddo
	    	!
	    	! hamiltonian in rotated wannier basis
	    	!
			wan%H_k(:,:,k) = matmul(transpose(conjg(wan%U(:,:,k))),matmul(E,wan%U(:,:,k)))
			!
		enddo
		!
	end subroutine get_hamiltonian

	!
	! functions which operate on hamiltonian
	!

	function 	   fft_kernel(K,R,direction) result(kernel)
		!
		! Eq. 25 PRB 65, 035109
		! Forward f(K) = Kernel * f(R) : K -> R : exp(-i*2*pi)
		!         mx1  =   mxn     nx1
		!
		! Eq. 26 PRB 65, 035109
		! Reverse f(R) = Kernel * f(K) : R -> K : exp( i*2*pi)
		!         mx1  =   mxn     nx1
		!
		implicit none
		!
    	real(dp), intent(in) :: K(:,:)
		real(dp), intent(in) :: R(:,:)
    	character(7), intent(in) :: direction ! forward/reverse
    	integer :: i,j,nx,nk
    	real(dp), allocatable :: kernel(:,:)
    	!
    	nx = size(R,2)
    	nk = size(K,2)
    	!
    	select case (direction)
    	case('forward')
    		! Eq. 25 PRB 65, 035109 : k -> r
	    	allocate(kernel(nx,nk))
			do i = 1, nx
			do j = 1, nk
			  	kernel(i,j) = exp( -cmplx_i*twopi*dot_product(R(:,i),K(:,j)) )
			enddo
			enddo
			kernel = kernel/real(nk,dp)
		case('reverse')
			! Eq. 26 PRB 65, 035109 : r -> k
			allocate(kernel(nk,nx))
			do j = 1, nk
			do i = 1, nx
			  	kernel(j,i) = exp(  cmplx_i*twopi*dot_product(R(:,i),K(:,j)) )
			enddo
			enddo
			! notice, no normalization factor in revese operation!
		case default
			call am_print('ERROR','FFT option unknown.',' >>> ')
			stop
		end select
	end function   fft_kernel

	function       fft_hamiltonian(H_k,kpt,rpt) result(H_r)
    	!
    	! Eq. 25 PRB 65, 035109 : k -> r
    	!
    	real(dp), intent(in) :: kpt(:,:)
    	real(dp), intent(in) :: rpt(:,:)
    	complex(dp), intent(in) :: H_k(:,:,:)
    	complex(dp), allocatable :: H_r(:,:,:)
    	integer :: m,n,nm,nn,nr
    	complex(dp), allocatable :: kernel(:,:)
    	!
    	nn = size(H_k,1)
    	nm = size(H_k,2)
    	nr = size(rpt,2)
    	!
    	allocate(H_r(nn,nm,nr))
    	!
    	kernel = fft_kernel(K=kpt,R=rpt,direction='forward')
    	!
		do m=1,nm
		do n=1,nn
		  	H_r(n,m,:) = matmul(kernel,H_k(n,m,:))
		enddo
		enddo
		!
	end function   fft_hamiltonian

	function       ifft_hamiltonian(H_r,kpt,rpt) result(H_k)
    	!
    	! Eq. 26 PRB 65, 035109 : r -> k
    	!
 		real(dp), intent(in) :: kpt(:,:)
    	real(dp), intent(in) :: rpt(:,:)
    	complex(dp), intent(in) :: H_r(:,:,:)
    	complex(dp), allocatable :: H_k(:,:,:)
    	integer :: m,n,nm,nn,nk
    	complex(dp), allocatable :: kernel(:,:)
    	!
    	nn = size(H_r,1)
    	nm = size(H_r,2)
    	nk = size(kpt,2)
    	!
    	allocate(H_k(nn,nm,nk))
    	!
    	kernel = fft_kernel(K=kpt,R=rpt,direction='reverse')
    	!
		do m=1,nm
		do n=1,nn
		  	H_k(n,m,:) = matmul(kernel,H_r(n,m,:))
		enddo
		enddo
		!
	end function   ifft_hamiltonian

	!
	! other
	!

	subroutine     wannier_development_space(wan,uc,bz,opts)
    	!
    	use am_brillouin_zone
    	!
    	implicit none
    	!
    	class(am_class_wannier), intent(inout) :: wan
    	type(am_class_bz), intent(in) :: bz
		type(am_class_unit_cell), intent(in) :: uc
		type(am_class_options) :: opts
		real(dp), allocatable :: grid_points(:,:)
    	!
	    call wan%load_amn(opts=opts)
	    !
	    call wan%load_eig(opts=opts)
	    !
	    call wan%get_hamiltonian()
	    !
	    grid_points = matmul( uc%bas, real(mesh_grid([2,2,2]),dp) )
	    !
	    call am_print('H',wan%H_k(:,:,1))
	    !
	    wan%H_r = fft_hamiltonian(H_k=wan%H_k,kpt=bz%kpt,rpt=grid_points)
		!
		wan%H_k = ifft_hamiltonian(H_r=wan%H_r,kpt=bz%kpt,rpt=grid_points)
		!
		call am_print('H',wan%H_k(:,:,1))

		!
		! test kernel (should get delta function.)
! 		!
! 		call am_print('test',&
! 		matmul( fft_kernel(X=grid_points,K=bz%kpt(:,1:10),direction='reverse'), &
! 			    fft_kernel(X=grid_points,K=bz%kpt(:,1:10),direction='forward') ) )



	end subroutine wannier_development_space








end module





