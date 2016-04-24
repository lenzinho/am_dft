module am_tight_binding

    use am_constants
    use am_helpers
    use am_matlab
    use am_options
    use am_shells
    use am_irre_cell
    use am_atom
    use am_symmetry
    use am_mkl

	implicit none

	private

! 	public :: tight_binding_main
! 	public :: test_orbitals
	public :: test_functions

contains

	pure function  slater_koster(l1,m1,l2,m2,R) result(sk)
		!
		! A. V. Podolskiy and P. Vogl, Phys. Rev. B 69, 233101 (2004). Eq 23
		!
		implicit none
		!
		integer , intent(in) :: l1,m1,l2,m2
		real(dp), intent(in) :: R(3)
		real(dp), allocatable :: sk(:)
		integer  :: mp
		real(dp) :: V(3)
		real(dp) :: L, N ! directional cosines
		real(dp) :: gamma
		real(dp) :: prefac
		!
		! V = [l,m,n] directional cosines
		V = R/norm2(R)
		L = V(1)
		N = V(3)
		gamma = -acos(L/sqrt(1.0_dp-N**2))
		!
		! One slater koster term per value of m between 0 and min(l1,l2)
		! 
		! sk(m+1) = ( min(l1,l2), max(l1,l2), m )
		!
		! For l1=0 and l2=1 ...
		! sk(0+1) = (0,1,0) = (spσ)
		! For l1=1 and l2=1 ...
		! sk(0+1) = (1,1,0) = (ppσ)
		! sk(1+1) = (1,1,1) = (ppπ)
		! For l1=1 and l2=2 ...
		! sk(0+1) = (1,2,0) = (pdσ)
		! sk(1+1) = (1,2,1) = (pdπ)
		! For l1=3 and l2=2 ...
		! sk(0+1) = (2,3,0) = (dfσ)
		! sk(1+1) = (2,3,1) = (dfπ)
		! sk(2+1) = (2,3,2) = (dfδ)
		!
		! l=0 (s), l=1 (p), l=2 (d), l=3 (f), l=4 (g)
		!
		allocate(sk(min(l1,l2)))
		!
		prefac=(-1.0_dp)**((l1-l2+abs(l1-l2))/2.0_dp)
		!
		! add m = 0 term first
		sk(1) = prefac*2.0_dp*Am(m1,gamma)*Am(m2,gamma)*dlmmp(l1,abs(m1),0,N)*dlmmp(l1,abs(m2),0,N)
		! 
		! add remaining terms
		do mp = 1, min(l1,l2)
			sk(mp+1) = prefac * &
				& (   slmmp(l1,m1,abs(mp),N,gamma)*slmmp(l2,m2,abs(mp),N,gamma) &
				&   + tlmmp(l1,m1,abs(mp),N,gamma)*tlmmp(l2,m2,abs(mp),N,gamma) )
		enddo
		!
		contains
		pure function dlmmp(l,m,mp,N)
			!
			! A. V. Podolskiy and P. Vogl, Phys. Rev. B 69, 233101 (2004). Eq 7
			!
			implicit none
			!
			integer , intent(in) :: l,m,mp
			real(dp), intent(in) :: N
			integer  :: t
			real(dp) :: dlmmp
			real(dp) :: fac1,fac2,fac3,fac4,fac5
			!
			fac1 = ((1.0_dp+N)/2.0_dp)**l
			fac2 = ((1.0_dp-N)/(1.0_dp+N))**((m-mp)/2.0_dp)
			fac3 = sqrt(real(factorial(l+mp)*factorial(l-mp)*factorial(l+m)*factorial(l-m),dp))
			fac4 = 0.0_dp
			do t = 0, (2*l+1)
				fac4 = fac4 + (-1.0_dp)**t/(factorial(l+mp-t)*factorial(l-m-t)*factorial(t)*factorial(t+m-mp))
			enddo
			fac5 = ((1.0_dp-N)/(1.0_dp+N))**t
			!
			dlmmp = fac1*fac2*fac3*fac4*fac5
			!
		end function  dlmmp
		pure function Am(m,gamma)
			!
			! A. V. Podolskiy and P. Vogl, Phys. Rev. B 69, 233101 (2004). Eq 21
			!
			implicit none
			!
			integer , intent(in) :: m
			real(dp), intent(in) :: gamma
			real(dp) :: Am
			!
			if (m.ne.0) then
				Am = (-1.0_dp)*abs(m)*( heavi(m)*cos(abs(m)*gamma) - heavi(-m)*sin(abs(m)*gamma) )
			else
				Am = 1/sqrt(2.0_dp)
			endif
			!
		end function  Am
		pure function Bm(m,gamma)
			!
			! A. V. Podolskiy and P. Vogl, Phys. Rev. B 69, 233101 (2004). Eq 21
			!
			implicit none
			!
			integer , intent(in) :: m
			real(dp), intent(in) :: gamma
			real(dp) :: Bm
			!
			Bm = (-1.0_dp)*abs(m)*( heavi(m)*cos(abs(m)*gamma) + heavi(-m)*sin(abs(m)*gamma) )
			!
		end function  Bm
		pure function slmmp(l,m,mp,N,gamma)
			!
			! A. V. Podolskiy and P. Vogl, Phys. Rev. B 69, 233101 (2004). Eq 24
			!
			implicit none
			!
			integer , intent(in) :: l,m,mp
			real(dp), intent(in) :: N, gamma
			real(dp) :: slmmp
			!
			slmmp = Am(m,gamma) * ( (-1.0_dp)**abs(mp)*dlmmp(l,abs(m),abs(mp),N) + dlmmp(l,abs(m),-abs(mp),N) )
			!
		end function  slmmp
		pure function tlmmp(l,m,mp,N,gamma)
			!
			! A. V. Podolskiy and P. Vogl, Phys. Rev. B 69, 233101 (2004). Eq 25
			!
			implicit none
			!
			integer , intent(in) :: l,m,mp
			real(dp), intent(in) :: N, gamma
			real(dp) :: tlmmp
			!
			if (m.eq.0) then
				tlmmp = 0.0_dp
			else
				tlmmp = Bm(m,gamma)*( (-1.0_dp)**abs(mp)*dlmmp(l,abs(m),abs(mp),N) - dlmmp(l,abs(m),-abs(mp),N) )
			endif
			!
		end function  tlmmp
			!
	end function   slater_koster

    function       spherical_harmonics(l,m,th,phi) result(Ylm)
        !
        ! computes spherical harmonics. th and phi in radians. th and phi should be the same size.
        !
        ! W. H. Press, B. P. Flannery, S. A. Teukolsky, and W. T. Vetterling, Numerical Recipes in Fortran 77: The Art
        ! of Scientific Computing, 2 edition (Cambridge University Press, Cambridge England ; New York, 1992), p 246.
        ! 
        implicit none
        !
        integer , intent(in)  :: l, m
        real(dp), intent(in) :: th(:), phi(:)
        complex(dp), allocatable :: Ylm(:)
        !
        allocate(Ylm(size(th)))
        !
        Ylm = sqrt( real(2*l+1,dp)/fourpi * factorial(l-m)/real(factorial(l+m),dp) ) * legendre(l,m,cos(th)) * exp(cmplx_i*m*phi)
        !
        Ylm = Ylm/am_dznrm2(Ylm)
        !
    end function   spherical_harmonics

    subroutine 	   setup_query_points_for_wavefunction(r,th,phi,flags)
    	!
    	implicit none
    	!
        character(*), intent(in), optional  :: flags
		real(dp), intent(out), allocatable :: r(:)   ! real space spherical coordinates
        real(dp), intent(out), allocatable :: th(:)  ! real space spherical coordinates
        real(dp), intent(out), allocatable :: phi(:) ! real space spherical coordinates
        real(dp), allocatable :: sph(:,:)
		integer :: ndims_r
		integer :: ndims_th
		integer :: ndims_phi ! dimensions for spheical grid on which wave unction is evaluated
        logical :: ismesh ! evaluates a mesh on the surface on the sphere rather than considerng phi and th independent
        logical :: issurf ! disregards radial contribution... will not produce wavefucnctions that are orthogonal in n quantum number
        !
        ismesh = .false.
        issurf = .false.
        if (present(flags)) then
            if (index(flags,'mesh')) then
                ismesh = .true.    
            endif
            if (index(flags,'surf')) then
                issurf = .true.
            endif
        endif
		!
		! make mesh for integration (because each component is independent, rather than a line is probed for each)
		! r = [0,inf); th = [0, pi]; phi= [0,2pi)
        ! 1000 already seems good enough...
        ndims_r = 1000
        ndims_th = 100
        ndims_phi = 100
        !
        if (issurf) then
            allocate(r(1))
            r = 10.0_dp
        else
            allocate(r(ndims_r))
            r = linspace(0.0_dp,100.0_dp,ndims_r) + 1.0D-14
        endif
        allocate(th(ndims_th))
        th = pi*[0:ndims_th]/real(ndims_th,dp)
        allocate(phi(ndims_phi-1))
        phi = 2*pi*[1:ndims_phi]/real(ndims_phi,dp)
        !
        if (ismesh) then
            !
            sph = meshgrid(th,phi)
            deallocate(th)
            deallocate(phi)
            allocate(th,source=sph(1,:))
            allocate(phi,source=sph(2,:))
            !
        endif
        !
    end subroutine setup_query_points_for_wavefunction

	subroutine     atomic_orbital(n,l,m,r,th,phi,psi_r,psi_th,psi_phi)
		!
		! returns real orthonormal wavefunctions
		! requires volume element such that when transpose(Y)_{n,l,m})*Y_{n',l',m'} = delta(n,n') delta(l,l') delta(m,m')
		! uses tesseral harmonics (wavefunction is real)
		!
		implicit none
		!
		integer , intent(in) :: n
		integer , intent(in) :: l
		integer , intent(in) :: m
	    real(dp), intent(in) :: r(:)   ! real space spherical coordinates
        real(dp), intent(in) :: th(:)  ! real space spherical coordinates
        real(dp), intent(in) :: phi(:) ! real space spherical coordinates
		real(dp), intent(out) :: psi_r(:)
        real(dp), intent(out) :: psi_th(:)
        real(dp), intent(out) :: psi_phi(:)
        !
		psi_r   = r_nl(n=n,l=l,r=r) * r
        psi_th  = legendre(l=l,m=m,x=cos(th)) * sqrt(abs(sin(th)))
        psi_phi = phi_m(m=m,phi=phi)
        !
        psi_r   = psi_r/am_dnrm2(psi_r)
        psi_th  = psi_th/am_dnrm2(psi_th)
        psi_phi = psi_phi/am_dnrm2(psi_phi)
        !
		contains
	    function     phi_m(m,phi)
	    	!
	        ! L. Pauling, E. B. W. Jr, and Physics, Introduction to Quantum Mechanics with Applications to Chemistry, 3134th
	        ! edition (Dover Publications, New York, N.Y, 1985). Eq 21.3 p 132, see also Table 21-1 p 133.
	        !
	        implicit none
	        !
	        integer , intent(in) :: m
	        real(dp), intent(in) :: phi(:)
	        real(dp), allocatable :: phi_m(:)
	        !
	        allocate(phi_m,mold=phi)
	        !
	        if     (m.lt.0) then; phi_m = sin(m*phi)
	    	elseif (m.gt.0) then; phi_m = cos(m*phi)
			elseif (m.eq.0) then; phi_m = 1.0_dp
			endif
			!
	    end function phi_m
		function     r_nl(n,l,r) result(y)
			!
			! Unnormalized radial hyrdogen wavefunction 
			!
			! R. Shankar, Principles of Quantum Mechanics, Softcover reprint of the original 1st ed. 1980 edition (Springer,
			! 2013), p 356. Griffiths, David J. (2005), Introduction to Quantum Mechanics, 2nd Edition; Pearson Education -
			! Problem 4.12.
			! L. Pauling, E. B. W. Jr, and Physics, Introduction to Quantum Mechanics with Applications to Chemistry, 3134th
			! edition (Dover Publications, New York, N.Y, 1985), p 131.
			!
			!             Laguerre L_n       functions orthogonal over R = [0,inf) with weighing function        exp(-r)
			! Generalized Laguerre L_{n}^{k} functions orthogonal over R = [0,inf) with weighing function r**k * exp(-r) 
			! Divide weight factor by two because the wave functions will multiply each other to produce unity.
			!
			! The radial components are built out of Laguerre polynomials, whose orthogonality only holds when leaving the
			! secondary index l fixed. <i,j|i,j> = 1, <i,k|j,k> = 0, and <i,j|k,l> for i/=j,k/=l is not 0 or 1.
			!
			real(dp), intent(in) :: r(:)
			integer , intent(in) :: n,l
			real(dp), allocatable :: y(:)
			integer :: k
			!
			allocate(y(size(r)))
			!
			! associated Laguerre polynomial L_k^p(x)(k,p,x) of degree k and order p
			y = laguerre(k=(n-l-1),p=(2*l+1),x=(r/n)) * (r/n)**l * exp(-r/(2*n))
			!
        end function r_nl
	end subroutine atomic_orbital

	subroutine     get_Ly_eigenstates(l, My, D)
		!
		! Calculates the eigenvalues D and eigevectors My of Ly operator in the Lz basis using raising/lowering operators.
		!
		! R. M. Martin, Electronic Structure: Basic Theory and Practical Methods, 1 edition
		! (Cambridge University Press, Cambridge, UK ; New York, 2008), p 573.
		!
		! Romero Nichols "Density Functional Study of Fullerene-Based Solids: Crystal Structure,
		! Doping, and Electron- Phonon Interaction", Ph.D. Thesis UIUC, p 96.
		!
		! C. Cohen-Tannoudji, B. Diu, and F. Laloe, Quantum Mechanics, 1 edition (Wiley-VCH, New
		! York; Paris, 1992), p 666.
		!
		! J. J. Sakurai, Modern Quantum Mechanics, Revised edition (Addison Wesley, Reading, Mass,
		! 1993). p 207
		!
		implicit none
		!
		integer    , intent(in) :: l ! l is the orbital quantum number,  
		complex(dp), intent(out), allocatable :: My(:,:) ! My is returned as a (l**2+1) x (l**2+1) matrix containing all possible values of m : |m| .le. l
		complex(dp), intent(out), allocatable :: D(:) ! eigenvalues
		complex(dp), allocatable :: Ly(:,:) !
		integer :: k, nstates, m
		! 
		! number of states (m=-l,-l+1,-l+2,...,0,...,l-2,l-1,l)
		nstates = l**2 + 1
		!
		! indices of My refer to value of m
		allocate(My(nstates,nstates))
		!
		! initialize banded Ly matrix
		allocate(Ly(2,nstates))
		Ly = 0
		!
		! Construct banded Ly matrix in the Lz basis from raising and lower operators
		! Laloe, p 666; Martin, p 573, Eq N4; Sakurai, p 207. Here, hbar is set to 1.
		k = 0
		do m = -l, l-1
			k=k+1
			Ly(2,k) = -0.5_dp * cmplx_i * sqrt( real(l*(l+1)-m*(m+1),dp) )
		enddo
		!
		call am_print('Ly',Ly)
		write(*,*) 'asdfasdfsd'
		!
		! Diagonalize the Ly matrix
		call am_zhbev(A_b=Ly,V=My,D=D)
		!
	end Subroutine get_Ly_eigenstates

	subroutine     test_functions()
		!
		implicit none
		integer :: l
		complex(dp), allocatable :: My(:,:)
		complex(dp), allocatable :: D(:)
		!
		l = 1
		!
		call get_Ly_eigenstates(l, My, D)
		!
		call am_print('My',My)
		!
		call am_print('D',D)
		!
	end subroutine test_functions

! 	subroutine     test_orbitals()
! 		!
! 		implicit none
! 		!
! 	    real(dp),allocatable :: r(:)   ! real space spherical coordinates
!         real(dp),allocatable :: th(:)  ! real space spherical coordinates
!         real(dp),allocatable :: phi(:) ! real space spherical coordinates
! 		integer ,allocatable :: state(:,:)
! 		real(dp),allocatable :: psi_r(:,:)
!         real(dp),allocatable :: psi_th(:,:)
!         real(dp),allocatable :: psi_phi(:,:)
!         real(dp),allocatable :: overlap(:,:)
!         real(dp),allocatable :: overlap_r(:,:)
!         real(dp),allocatable :: overlap_th(:,:)
!         real(dp),allocatable :: overlap_phi(:,:)
!         complex(dp), allocatable :: My(:,:)
! 		integer :: n,l,m
! 		integer :: nstates
! 		integer :: i
! 		!
! 		!
! 		allocate(state(3,100))
! 		nstates=0
! 		do n = 1,2
! 		do l = 0,n-1
! 		do m = -l,l
! 			nstates=nstates+1
! 			state(:,nstates)=[n,l,m]
! 		enddo
! 		enddo
! 		enddo 
! 		call am_print('state',state(:,1:nstates))
! 		!
! 		! setup points on which wave function will be evaluated
! 		call setup_query_points_for_wavefunction(r=r,th=th,phi=phi,flags='mesh,surf')
! 		!
! 		! evalute wave functions
! 		allocate(psi_r(size(r),nstates))
!         allocate(psi_th(size(th),nstates))
!         allocate(psi_phi(size(phi),nstates))
! 		do i = 1,nstates
! 			!
! 			n = state(1,i)
! 			l = state(2,i)
! 			m = state(3,i)
! 			!
! 			call atomic_orbital(n=n,l=l,m=m,r=r,th=th,phi=phi,psi_r=psi_r(:,i),psi_th=psi_th(:,i),psi_phi=psi_phi(:,i))
!             !
!         enddo
! 		!
! 		! calculate overlaps
!         overlap_r=am_dgemm(psi_r,psi_r,flags='AT')
!         overlap_th=am_dgemm(psi_th,psi_th,flags='AT')
!         overlap_phi=am_dgemm(psi_phi,psi_phi,flags='AT')
! 		call am_print('r overlap',overlap_r)
!         call am_print('th overlap',overlap_th)
!         call am_print('phi overlap',overlap_phi)
!         !
!         ! calculated final overlap
!         allocate(overlap(nstates,nstates))
!         overlap = overlap_r * overlap_th * overlap_phi
!         call am_print('full overlap',overlap)
!         !
!         call am_print('real(My)',Real(My))
!         call am_print('complex(My)',Aimag(My))
!         !
! 	end subroutine test_orbitals

	subroutine     test_spherical_harmonic_expansion()
		!
		implicit none
		!
	    real(dp),allocatable :: r(:)   ! real space spherical coordinates
        real(dp),allocatable :: th(:)  ! real space spherical coordinates
        real(dp),allocatable :: phi(:) ! real space spherical coordinates
		integer ,allocatable :: state(:,:)
		real(dp),allocatable :: psi_r(:,:)
        real(dp),allocatable :: psi_th(:,:)
        real(dp),allocatable :: psi_phi(:,:)
        real(dp),allocatable :: overlap(:,:)
        real(dp),allocatable :: overlap_r(:,:)
        real(dp),allocatable :: overlap_th(:,:)
        real(dp),allocatable :: overlap_phi(:,:)
        complex(dp), allocatable :: My(:,:)
		integer :: n,l,m
		integer :: nstates
		integer :: i
		!
		!
		allocate(state(3,100))
		nstates=0
		do n = 3,3
		do l = 0,n-1
		do m = -l,l
			nstates=nstates+1
			state(:,nstates)=[n,l,m]
		enddo
		enddo
		enddo 
		call am_print('state',state(:,1:nstates))
		!
		! setup points on which wave function will be evaluated
		call setup_query_points_for_wavefunction(r=r,th=th,phi=phi,flags='mesh,surf')
		!
		! evalute wave functions
		allocate(psi_r(size(r),nstates))
        allocate(psi_th(size(th),nstates))
        allocate(psi_phi(size(phi),nstates))
		do i = 1,nstates
			!
			n = state(1,i)
			l = state(2,i)
			m = state(3,i)
			!
			call atomic_orbital(n=n,l=l,m=m,r=r,th=th,phi=phi,psi_r=psi_r(:,i),psi_th=psi_th(:,i),psi_phi=psi_phi(:,i))
            !
        enddo
		!
		! calculate overlaps
        overlap_r=am_dgemm(psi_r,psi_r,flags='AT')
        overlap_th=am_dgemm(psi_th,psi_th,flags='AT')
        overlap_phi=am_dgemm(psi_phi,psi_phi,flags='AT')
		call am_print('r overlap',overlap_r)
        call am_print('th overlap',overlap_th)
        call am_print('phi overlap',overlap_phi)
        !
        ! calculated final overlap
        allocate(overlap(nstates,nstates))
        overlap = overlap_r * overlap_th * overlap_phi
        call am_print('full overlap',overlap)
        !
        call am_print('real(My)',Real(My))
        call am_print('complex(My)',Aimag(My))
        !
	end subroutine test_spherical_harmonic_expansion

!  	subroutine     tight_binding_main(irrpair,ic,opts)
!  		!
!  		implicit none
!         !
!         class(am_class_irre_cell):: ic ! irreducible cell
!         type(am_class_pair_shell):: irrpair ! irreducible pairs
!         type(am_class_options)   :: opts
!         real(dp),allocatable :: ck1(:,:) ! atom 1: ck(i,k) are the k expansion coefficients for cartesian atomic orbital i quantize with respect to the molecular coordinate system (bond vector aligned along z)
!         real(dp),allocatable :: ck2(:,:) ! atom 2: ck(i,k) are the k expansion coefficients for cartesian atomic orbital i quantize with respect to the molecular coordinate system (bond vector aligned along z)
!         real(dp),allocatable :: dmat(:,:) ! delta matrix
!         integer :: i,j,k
!         !
!         if (opts%verbosity.ge.1) call am_print_title('Quantizing irreducible orbitals with respect to irreducible bonds (pair shells)')
!         !
!         ! generate orbitals or load them from input
!         ! note, this procedure below does not necessarily generate all m degenerates states for a given l.
!         ! this, however is a necessary condition when quantizing the orbital with respect to the bond vector
!         ! need to modify it below. requested input should be s p d f and their combinations s,sp,spd, spds* , etc.
!         allocate(ic%atom(ic%natoms))
!         do i = 1, ic%natoms
!             ! last number in nlms array is the state index starting with H 1s
!             call ic%atom(i)%generate_states(Z=ic%Z(i),n_valence_states=4,n_excited_states=4)
!         enddo
!         !
!         !
!         do k = 1, irrpair%nshells
!             i = irrpair%shell(k)%i
!             j = irrpair%shell(k)%j
!             !
!             if (opts%verbosity.ge.1) then
!                 call am_print('quantizing state',unique(ic%atom(i)%state(1:3,:)))
!             endif

!             ! get expansion coefficients from state = [n,l,m,s,#]
!             ! the unique() ignores the spin quantum number
!             ck1 = quantize_to_bond(v=irrpair%shell(k)%tau(1:3,1), state=unique(ic%atom(i)%state(1:3,:)))
!             ck2 = quantize_to_bond(v=irrpair%shell(k)%tau(1:3,1), state=unique(ic%atom(j)%state(1:3,:)))
!             !
!             ! form outer product of ck1(i,:) and ck2(j,:) for all i and j
            
!             stop
!         enddo

!     end subroutine tight_binding_main
 !






	! stuff that is broken
	! stuff that is broken
	! stuff that is broken
	! stuff that is broken




! 	pure function spdf2l(spdf) result(l)
! 		!
! 		character(1), intent(in) :: spdf
! 		integer :: l
! 		!
! 		if     (spdf(1:1).eq.'s') then; l = 0
! 		elseif (spdf(1:1).eq.'p') then; l = 1
! 		elseif (spdf(1:1).eq.'d') then; l = 2
! 		elseif (spdf(1:1).eq.'f') then; l = 3
! 		endif
! 	end function  spdf2l





! 	function       quantize_to_bond(v,l,m) result(overlap)
! 		!
!         ! quantize all states on a given atom with respect to vector v
!         ! 
! 		implicit none
! 		!
! 		real(dp), intent(in) :: v(3)   ! bond vector
! 		integer , intent(in) :: l(:)   ! l and m of the state that is to be quantized with respect to v
! 		integer , intent(in) :: m(:)   ! l(:) and m(:) vectors are used to probe multiple states at once.
! 	    real(dp),allocatable :: r(:)   ! real space spherical coordinates
!         real(dp),allocatable :: th(:)  ! real space spherical coordinates
!         real(dp),allocatable :: phi(:) ! real space spherical coordinates
! 		real(dp),allocatable :: psi1_r(:,:)
!         real(dp),allocatable :: psi1_th(:,:)
!         real(dp),allocatable :: psi1_phi(:,:)
! 		real(dp),allocatable :: psi2_r(:,:)
!         real(dp),allocatable :: psi2_th(:,:)
!         real(dp),allocatable :: psi2_phi(:,:)
!         real(dp),allocatable :: overlap(:,:) ! overlap corresponds to expansion coefficients
!         real(dp),allocatable :: overlap_r(:,:)
!         real(dp),allocatable :: overlap_th(:,:)
!         real(dp),allocatable :: overlap_phi(:,:)
!         real(dp):: dcosines(3)
!         real(dp):: offset_th
!         real(dp):: offset_phi
! 		integer :: np,lp,mp
! 		integer :: nstates
! 		integer :: i
!         !
! 		! setup points on which wave function will be evaluated
! 		call setup_query_points_for_wavefunction(r=r,th=th,phi=phi,flags='mesh,surf')
! 		!
! 		! setup s,p,d,f states in cartesian coordinates
! 		! l = 0 , states = 1; for n = 1, total states = 1
! 		! l = 1 , states = 3; for n = 2, total states = 4
! 		! l = 2 , states = 5; for n = 3, total states = 9
! 		! l = 3 , states = 7; for n = 4, total states = 16
! 		np = 4
! 		nstates = np**2
! 		allocate(state(3,nstates))
! 		do lp = 0,np-1
! 		do mp = -lp,lp
! 			state(:,nstates)=[np,lp,mp]
! 		enddo
! 		enddo
! 		!
! 		! setup atomic orbitals in cartesian coordinates
! 		allocate(psi1_r(size(r),1))
!         allocate(psi1_th(size(th),1))
!         allocate(psi1_phi(size(phi),1))
! 		do i = 1,nstates
! 			!
! 			np = state(1,i)
! 			lp = state(2,i)
! 			mp = state(3,i)
! 			!
! 			call atomic_orbital(n=np,l=lp,m=mp,r=r,th=th,phi=phi, &
!                 psi_r  =psi1_r(:,i), &
!                 psi_th =psi1_th(:,i),&
!                 psi_phi=psi1_phi(:,i))
!             !
!         enddo
!         !
!         ! rotate coordinate system ...
!         ! z axis is at th = 0. this operation below brings vector v parallel with z axis, pointing up.
!         ! if v=[0,0,0] no quantization is needed. keep vector aligned as it already is, i.e. z pointing along z.
!         if (norm2(v).gt.tiny) then
!             dcosines = v/norm2(v)
!         else
!             dcosines = real([0,0,1],dp)
!         endif
!         !
!         offset_th = acos(dcosines(3))
!         offset_phi = atan(dcosines(2)/(dcosines(1)+1.0D-14))
!         !
!         th  = th - offset_th
!         phi = phi - offset_phi
!         !
!         call am_print('v',v)
!         call am_print('directional cosines',dcosines)
!         call am_print('offset theta',offset_th)
!         call am_print('offset phi',offset_phi)
!         !
!         ! determine atomic orbitals with respect to bond coordinate system
!         ! define new z axis to be parallel with bond. Recall R*f(x) = f(inv(R)*x))
! 		nstates = size(l)
! 		allocate(psi2_r(size(r),nstates))
!         allocate(psi2_th(size(th),nstates))
!         allocate(psi2_phi(size(phi),nstates))
! 		do i = 1, nstates
! 			!
! 			lp = l(i)
! 			mp = m(i)
!             !
! 			call atomic_orbital(n=np,l=l,m=m,r=r,th=th,phi=phi,&
!                 psi_r  =psi2_r(:,i), &
!                 psi_th =psi2_th(:,i),&
!                 psi_phi=psi2_phi(:,i))
!             !
!         enddo
!         !
! 		! calculate overlaps
!         !overlap_r=matmul(transpose(psi1_r),psi2_r)
!         !overlap_th=matmul(transpose(psi1_th),psi2_th)
!         !overlap_phi=matmul(transpose(psi1_phi),psi2_phi)
!         overlap_r=am_dgemm(  psi1_r,  psi2_r,  flags='AT')
!         overlap_th=am_dgemm( psi1_th, psi2_th, flags='AT')
!         overlap_phi=am_dgemm(psi1_phi,psi2_phi,flags='AT')
!         !
! 		call am_print('r overlap',overlap_r)
!         call am_print('th overlap',overlap_th)
!         call am_print('phi overlap',overlap_phi)
!         !
!         ! calculated final overlap
!         allocate(overlap(nstates,nstates))
!         overlap = overlap_r * overlap_th * overlap_phi
!         !
!         call am_print('full overlap',overlap)
! 	end function   quantize_to_bond


end module am_tight_binding











