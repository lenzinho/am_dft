module am_wannier
	!
	use am_unit_cell
	use am_constants
	use am_stdout
	use am_matlab
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
		procedure :: fourier_interpolation
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
    		! hamiltonian in diagonalized bloch-state basis
	    	E=0.0_dp
	    	do i = 1,wan%nbands
		    	E(i,i) = wan%E(i,k)
	    	enddo
	    	! hamiltonian in rotated wannier basis
			wan%H_k(:,:,k) = matmul(transpose(conjg(wan%U(:,:,k))),matmul(E,wan%U(:,:,k)))
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
    	real(dp), allocatable :: RdotK(:,:)
    	!
    	nx = size(R,2)
    	nk = size(K,2)
    	!
    	allocate(RdotK(nx,nk))
    	allocate(kernel(nx,nk))
		do i = 1, nx
		do j = 1, nk
		  	RdotK(i,j) = dot_product(R(:,i),K(:,j))
		enddo
		enddo
		!
    	select case (direction)
    	case('forward')
    		! Eq. 25 PRB 65, 035109 : k -> r
	    	kernel = exp(-cmplx_i*twopi*RdotK)
			kernel = kernel/real(nk,dp)
		case('reverse')
			! Eq. 26 PRB 65, 035109 : r -> k
			kernel = exp(+cmplx_i*twopi*transpose(RdotK))
			! notice, no normalization factor in revese operation!
		case default
			call am_print('ERROR','FFT option unknown.',flags='E')
			stop
		end select
		!
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

	function       get_fft_mesh(kpt) result(fourier_points)
		!
		implicit none
		!
		real(dp) :: kpt(:,:)
		real(dp), allocatable :: fourier_points(:,:)
		logical :: mask(size(kpt,2))
		integer :: n(3)
		integer :: i
		!
		do i = 1,3
			mask = (abs(kpt(i,:)).gt.tiny)
			n(i) = nint(1/(minval(abs(pack(kpt(i,:),mask)))))
		enddo
		!
	    fourier_points = meshgrid( [0:n(1)]-floor(n(1)/2.0_dp),[0:n(1)]-floor(n(2)/2.0_dp),[0:n(1)]-floor(n(3)/2.0_dp) )
		!
	end function   get_fft_mesh

	!
	! other
	!

	subroutine     fourier_interpolation(wan,uc,bz,opts)
    	!
    	use am_brillouin_zone
    	!
    	implicit none
    	!
    	class(am_class_wannier), intent(inout) :: wan
    	type(am_class_bz), intent(inout) :: bz
		type(am_class_unit_cell), intent(in) :: uc
		type(am_class_options) :: opts
		real(dp), allocatable :: fourier_points(:,:)
		!
	    ! should use a gamma-centered odd mesh
	    ! gamma provides the dc component; the pairs of - and + provide the cosine and sine components.
	    call am_print('1.0_dp/minval(bz%kpt(1,:))',nint(1.0_dp/minval(bz%kpt(1,:))))
	    call am_print('nint(1.0_dp/minval(bz%kpt(1,:)))',nint(1.0_dp/minval(bz%kpt(1,:))))
	    call am_print('nint(1.0_dp/minval(bz%kpt(1,:)))/2',nint(1.0_dp/minval(bz%kpt(1,:)))/2)
	    !
	    fourier_points = get_fft_mesh(bz%kpt)
	    !
	    call am_print('fourier_points',transpose(fourier_points))
	    !
	    call am_print('H',wan%H_k(:,:,1))
	    !
	    wan%H_r = fft_hamiltonian(H_k=wan%H_k,kpt=bz%kpt,rpt=fourier_points)
		!
		wan%H_k = ifft_hamiltonian(H_r=wan%H_r,kpt=bz%kpt,rpt=fourier_points)
		!
		call am_print('H',wan%H_k(:,:,1))

		!
		! test kernel (should get delta function.)
		! clear; clc
		! % THIS IS THE FFT ALGORITHM:
		! % X MUST BE DEFINED BETWEEN 0 and 1 WITHOUT DOUBLE COUNTING EDGES
		! % i.e. X must be defined within the brillouin zone, no double counting boundaries
		!
		! % points inside the BZ (not carefully count edge so that no points are included twice)
		! % best to use cartesian coordinates here. figure out how to do it. cart
		! % will let the bz assume the nice wigner seitz cell shape.
		! x = linspace(-0.5,0.5,100); x=x(1:end-1);
		! k(1,:)=x(:); k(3,:)=0;
		!
		! % energy values to interpolate
		! E(:,1) = sin(x*2*pi)+sin(x*2*pi*4);
		!
		! plot(x,E,'o-')
		!
		! % Y MUST BE DEFINED AT INTEGER VALUES and MUST BE SYMMETRY AT 0.
		! % These are real space lattice sites centered, use wigner-seitz BZ and
		! % wigner-seitz real-space superlattice. Number of real space points
		! % determine the number of fourier components. 
		! % FIRST TRY RE-IMPLEMENTING BOLTZTRAP. FOURIER INTERPOLATION.
		!
		! % fourier components corresponging to real-space primitive lattice sites
		! n = 20;
		! y=[-n:1:n];
		! r(1,:)=y(:); r(3,:)=0;
		!
		! % start main script here
		!
		! nk=size(k,2);
		! nr=size(r,2);
		!
		! % F converts k -> r
		! for i = 1:nk
		! for j = 1:nr
		!     F(j,i) = exp( -sqrt(-1)*2*pi*dot(k(:,i),r(:,j)) );
		! end
		! end
		! F = F/nk; % NOTICE ONLY FORWARD TRANSFER HAS NORMALIZATION. See Eq. 25.
		!
		! % WHEN RESAMPLED. THE K NORMALIZATION FACTOR NEEDS TO BE THE SAME.
		!
		! % RESAMPLE K HERE
		! x = linspace(0,1,100); x=x(1:end-1);
		! clear k; k(1,:)=x(:); k(3,:)=0; nk=size(k,2);
		!
		! % I converts r -> k
		! for i = 1:nk
		! for j = 1:nr
		!     I(i,j) = exp( +sqrt(-1)*2*pi*dot(k(:,i),r(:,j)) );
		! end
		! end
		! % INVERSE TRANSFORM DOES NOT HAVE NORMALIZATION.
		!
		! % Fourier components are obtained from:
		! fourer_components(1:nr) = F*E;
		!
		! % norm(I*(F*E)-E) % ~ 1E-15
		! hold on;
		! plot(x,I*(F*E),'*')
		!
		! hold off;
		!
		! % surf(abs(F*I),'edgecolor','none'); view([0 0 1])
		! %%
		!




! 		!
! 		call am_print('test',&
! 		matmul( fft_kernel(X=grid_points,K=bz%kpt(:,1:10),direction='reverse'), &
! 			    fft_kernel(X=grid_points,K=bz%kpt(:,1:10),direction='forward') ) )



	end subroutine fourier_interpolation








end module





