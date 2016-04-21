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

	! public :: tight_binding_main
	public :: test_orbitals

contains

	pure function slater_koster(l1,m1,l2,m2,R) result(sk)
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
		! One slater koster term per value of m betwen 0 and min(l1,l2)
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
	end function  slater_koster

    function     spherical_harmonics(l,m,theta,phi) result(Ylm)
        !
        ! computes spherical harmonics. Theta and phi in radians. Theta and phi should be the same size.
        !
        ! W. H. Press, B. P. Flannery, S. A. Teukolsky, and W. T. Vetterling, Numerical Recipes in Fortran 77: The Art
        ! of Scientific Computing, 2 edition (Cambridge University Press, Cambridge England ; New York, 1992), p 246.
        ! 
        implicit none
        !
        integer , intent(in)  :: l, m
        real(dp), intent(in) :: theta(:), phi(:)
        complex(dp), allocatable :: Ylm(:)
        !
        allocate(Ylm(size(theta)))
        !
        Ylm = sqrt( real(2*l+1,dp)/fourpi * factorial(l-m)/real(factorial(l+m),dp) ) * legendre(l,m,cos(theta)) * exp(cmplx_i*m*phi)
        !
        Ylm = Ylm/am_dznrm2(Ylm)
        !
    end function spherical_harmonics

    function     tesseral_harmonics(l,m,theta,phi) result(Ylm)
        ! L. Pauling, E. B. W. Jr, and Physics, Introduction to Quantum Mechanics with Applications to Chemistry, 3134th
        ! edition (Dover Publications, New York, N.Y, 1985). Eq 21.3 p 132, see also Table 21-1 p 133.
        !
        implicit none
        !
        integer , intent(in)  :: l, m
        real(dp), intent(in) :: theta(:), phi(:)
        real(dp), allocatable :: Ylm(:)
        !
        allocate(Ylm(size(theta)))
        !
        if     (m.lt.0) then; Ylm = legendre(l,abs(m),cos(theta)) *  sin(m*phi)
    	elseif (m.gt.0) then; Ylm = legendre(l,abs(m),cos(theta)) *  cos(m*phi)
		elseif (m.eq.0) then; Ylm = legendre(l,abs(m),cos(theta))
		endif
        !
        ! normalize wave function
        ! Ylm = Ylm/am_dnrm2(Ylm)
        !
    end function tesseral_harmonics

	function     Rnl(r,n,l) result(y)
		! Unnormalized radial hyrdogen wavefunction 
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
		! associated Laguerre polynomial L_k^p(x) of degree k and order p
		! the last r is part of the volume integration factor for spherical integration
		y = r**(l/2) * exp(-r/2) * laguerre(k=(n-l-1),p=(2*l+1),x=(r/n))
		! 
		! FROM LINUS PAULING'S BOOK:
		! y=0
		! do k = 0, (n-l-1)
		! y = y + (-1)**(k+1) * factorial(n+l)**2/(factorial(n-l-k-1)*factorial(2*l+1+k)*factorial(k)) * (r/n)**k
		! enddo
		! y = y * exp(-r/(2*n)) ! * (r/n)**l
		!
	end function Rnl

	subroutine    test_orbitals()
		!
		implicit none
		!
	    real(dp),allocatable :: sph(:,:)  ! real space spherical coordinates sph(1:3,i) = [r,th,phi] 
		integer ,allocatable :: nlm(:,:)
		real(dp),allocatable :: Y(:,:)
		real(dp),allocatable :: dV(:) ! integration volume element
		integer :: n,l,m
		integer :: k, i
		!
		!
		n = 3 ! 180/4
		! r     = [0,inf)
		! theta = [0, pi]
		! phi   = [0,2pi)
		sph = meshgrid(linspace(0.0_dp,10000.0_dp, 2 ), pi*[0:n]/real(n,dp), [pi*[1:(2*n)]/real(2*n,dp)] )
		! volume element: in reality dV = sqrt(dV) because it will multiply both bra and ket wave functions to form recover its power in the integral
		dV  = sph(1,:) * aimag(exp(cmplx_i*sph(2,:)/2))
		!
		allocate(nlm(3,100))
		k=0
		do n = 3,3
		do l = 0,n-1
		do m = -l,l
			k=k+1
			nlm(:,k)=[n,l,m]
		enddo
		enddo
		enddo 
		call am_print('nlm',nlm(:,1:k))
		!
		!
		allocate(Y(size(sph,2),k))
		Y = 0
		!
		do i = 1,k
			n = nlm(1,i)
			l = nlm(2,i)
			m = nlm(3,i)
			Y(:,i) = tesseral_harmonics(l=l,m=m,theta=sph(2,:),phi=sph(3,:)) * dV
			! Y(:,i) = Rnl(n=n,l=l,r=sph(1,:))
			Y(:,i) = Y(:,i)/am_dnrm2(Y(:,i))
		enddo
		!
		! Y = get_orbitals(x=xyz(1,:),y=xyz(2,:),z=xyz(3,:),nlm=nlm(:,1:k))
		!
		call am_print('Y^T*Y',matmul(transpose((Y)),Y))




	end subroutine test_orbitals


! 	subroutine tight_binding_main(ic,shell,pg,opts)
! 		!
! 		implicit none
! 		!
! 		class(am_class_irre_cell):: ic
! 		type(am_class_symmetry)  :: pg ! point group
! 	    type(am_class_options)   :: opts
! 	    type(am_class_pair_shell):: shell
!     	real(dp),allocatable :: r(:)      ! radial coordinates 
! 	    real(dp),allocatable :: sph(:,:)  ! real space spherical coordinates sph(1:3,i) = [r,th,phi]
! 	    real(dp),allocatable :: xyz(:,:)  ! real space cartsian  coordinates xyz(1:3,i) = [x,y,z]
! 	    real(dp),allocatable :: xyzp(:,:) ! real space cartesian coordinates xyz(1:3,i) = [x,y,z] rotated 
! 	    real(dp),allocatable :: Y(:,:)    ! Y(r,lm combined indices)
! 	    integer ,allocatable :: nlms(:,:) ! quantum numbers
! 	    real(dp),allocatable :: s(:)
! 	    real(dp) :: Rot(3,3)
! 	    integer :: i, j, k
! 	    integer :: N ! spherical coordinates partitions
! 		!
! 		! pair_prototypes = shell%get_pair_prototypes(ic=ic,pg=pg,opts=opts)
! 		!
! 		if (opts%verbosity.ge.1) call am_print_title('Determining tight-binding parameters')
! 		!
! 		! setup real space coordinates
! 		! Note: d_z^2 requires a very large range in order to be orthogonal with s orbital
! 		!
!         ! something is wrong with spherical integration...
!         !
! 		! setup spherical coordinates
! ! 		N = 2**4
! !         write(*,*) N
! !         ! r, 0 <= theta <= pi , 0 <= phi < 2*pi
! ! 		sph = meshgrid(linspace(0.0_dp,100000.0_dp,3), pi*[0:n]/real(n,dp), twopi*[0:(n-1)]/real(n,dp))
! ! 		! convert to cartesian coordinates
! ! 		allocate(xyz,mold=sph)
! ! 		xyz(1,:)=sph(1,:)*cos(sph(2,:))*sin(sph(3,:))
! ! 		xyz(2,:)=sph(1,:)*sin(sph(2,:))*sin(sph(3,:))
! ! 		xyz(3,:)=sph(1,:)*cos(sph(3,:))
!         !
!        	N = 100
!   		xyz = meshgrid(linspace(-N,N,N),linspace(-N,N,N),linspace(-N,N,N))
! 		!
! 		if (opts%verbosity.ge.1) call am_print('number of irreducible atoms',ic%natoms)
! 		!
! 		!
! 		i = 1
! 		nlms = get_states(Z=ic%Z(i),nvalences=10,nexciteds=4)
! 		!
! 		! get cubic harmonics
! 		! get_orbitals(nlm,x,y,z)
! 		Y = get_orbitals(x=xyz(1,:),y=xyz(2,:),z=xyz(3,:),nlm=unique(nlms(1:3,:)))
! 		!
! 		call am_print('orbitals',unique(nlms(1:3,:)))
! 		! call am_print('Y',Y)
! 		call am_print('Y^T*Y',matmul(transpose(Y),Y))
! 		!
        
        
!         stop
        
        
        
! 		! confirm that states are orthogonal
! 		if (any( abs(am_dgemm(A=Y,B=Y,flags='AT')-eye(k)).gt.tiny )) then
! 			call am_print('ERROR','Cubic harmonics not orthonormal.')
! 			call am_print('Y^T*Y', am_dgemm(A=Y,B=Y,flags='AT') )
! 			call am_print('log(abs(Y^T*Y-I))', log(abs(am_dgemm(A=Y,B=Y,flags='AT')-eye(k))) )
! 			stop
! 		endif

! ! 		! get atomic orbitals of intrest (eventually this will be an input)
! ! 		! nvo = number of valence orbital
! ! 		! neo = number of excited orbitals
! ! 		allocate(ic%atom(ic%natoms))
! ! 		do i = 1, ic%natoms
! ! 			call ic%atom(i)%generate_orbitals(Z=ic%Z(i),opts=opts,nvalences=4,nexciteds=4)
! ! 			! last number in nlms array is the state index starting with H 1s
! ! 		enddo


	

! ! 		do i = 1, 1 ! ic%natoms
! ! 			lm = unique(ic%atom(i)%)
		
! ! 		!
! ! 		! quantize each orbital with respect to bond axes
! ! 		allocate(r_rot,mold=r)
! ! 		do i = 1, 1  !, ic%natoms
! ! 			! get rotation matrix which aligns Z axis parallel to bond
! ! 			Rot = rotmat(a=real([0,0,1],dp),b=[1,1,1]/real(4,dp))
! ! 			r_rot = matmul(Rot,r)

! ! 		enddo
! ! 		stop





! 	end subroutine tight_binding_main




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



! 	function      get_orbitals(nlm,x,y,z) result(orb)
! 		!
! 		implicit none
! 		!
! 		real(dp), intent(in) :: x(:)
! 		real(dp), intent(in) :: y(:)
! 		real(dp), intent(in) :: z(:)
! 		integer , intent(in) :: nlm(:,:) ! nlm(:,nstates) [n,l,m] for each state
! 		integer :: n,l,m
! 		real(dp), allocatable :: orb(:,:)
! 		real(dp), allocatable :: r(:)
! 		integer :: i
! 		integer :: nstates
! 		!
! 		nstates=size(nlm,2)
! 		!
! 		allocate(r,source=sqrt( x**2 + y**2 + z**2 + 1.0D-14 ))
! 		!
! 		allocate(orb(size(x),nstates))
! 		orb = 0
! 		!
! 		do i = 1, nstates
! 			n = nlm(1,i)
! 			l = nlm(2,i)
! 			m = nlm(3,i)
! 			! psi = (cubic harmonics) * (radial part)
! 			orb(:,i) = Xlc(l=l,m=m,x=x,y=y,z=z,r=r) * Rnl(n=n,l=l,r=r) 
! ! 			orb(:,i) =  Rnl(n=n,l=l,r=r) 
! 			!
! 			! normalize
! 			orb(:,i) = orb(:,i)/am_dnrm2(orb(:,i))
! 		enddo
! 		!
! 	end function  get_orbitals



end module am_tight_binding











