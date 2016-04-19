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

	public :: tight_binding_main

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

	pure function heavi(m)
		!
		! A. V. Podolskiy and P. Vogl, Phys. Rev. B 69, 233101 (2004). Eq 15
		!
		implicit none
		!
		integer, intent(in) :: m
		real(dp) :: heavi
		!
		if (m.ge.0) then
			heavi = 1.0_dp
		else
			heavi = 0.0_dp
		endif
		!
	end function  heavi

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
	end function  slater_koster

    pure function Ylm(l,m,theta,phi)
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
        Ylm = Ylm/sqrt(dot_product(conjg(Ylm),Ylm))
        !
    end function  Ylm

	function      Xlc(l,m,x,y,z,r)
		! unnormalized cubic harmonics
		! R. McWeeny, Symmetry: An Introduction to Group Theory and Its Applications (Courier Corporation, 2002). p 263.
        ! See also https://en.wikipedia.org/wiki/Table_of_spherical_harmonics
		!
		real(dp), intent(in) :: x(:)
		real(dp), intent(in) :: y(:)
		real(dp), intent(in) :: z(:)
		real(dp), intent(in) :: r(:)
		integer , intent(in) :: l
		integer , intent(in) :: m
		real(dp), allocatable :: Xlc(:)
		! real(dp), allocatable :: r(:)
		!
		allocate(Xlc(size(x)))
		!
		! allocate(r,source=sqrt( x**2 + y**2 + z**2 + 1.0D-14 ))
		!
		if     (l.eq. 0) then ! s
							       Xlc = 1.0_dp
		elseif (l.eq. 1) then ! p
			if     (m.eq.-1) then; Xlc = z/r ! z
			elseif (m.eq. 0) then; Xlc = x/r ! x
			elseif (m.eq. 1) then; Xlc = y/r ! y
			endif
        elseif (l.eq. 2) then ! d
            ! 2 -2 does seem to orthogonal to 0,0
			if     (m.eq.-2) then; Xlc = (3*z**2)/r**2-1.0_dp ! ~ z^2 ! E_g
			elseif (m.eq.-1) then; Xlc = (x**2-y**2)/r**2 ! x^2-y^2   ! E_g
			elseif (m.eq. 0) then; Xlc = (x*y)/r**2 ! xy ! d_2tg
			elseif (m.eq. 1) then; Xlc = (z*x)/r**2 ! zx ! d_2tg
			elseif (m.eq. 2) then; Xlc = (y*z)/r**2 ! yz ! d_2tg
			endif
		else 
			call am_print('ERROR','d and f states not yet supported.')
			stop
		! These are not orthogonal for some reason?
		! 		elseif (l.eq. 3) then ! f
		! 			if     (m.eq.-3) then; Xlc = z*(x**2-y**2)/r**3 ! 
		! 			elseif (m.eq.-2) then; Xlc = x*(y**2-z**2)/r**3 ! 
		! 			elseif (m.eq.-1) then; Xlc = y*(z**2-x**2)/r**3 ! 
		! 			elseif (m.eq. 0) then; Xlc = x*y*z/r**3
		! 			elseif (m.eq. 1) then; Xlc = z*(5*z**2 - 3*r**2)/r**3 ! z^2
		! 			elseif (m.eq. 2) then; Xlc = x*(5*x**2 - 3*r**2)/r**3 ! x^2
		! 			elseif (m.eq. 3) then; Xlc = y*(5*y**2 - 3*r**2)/r**3 ! y^2
		! 			endfi
		endif
		!
	end function  Xlc

	function      radial(n,l,r) result(y)
		! unnormalized radial hyrdogen wavefunction
		! R. Shankar, Principles of Quantum Mechanics, Softcover reprint of the original 1st ed. 1980 edition (Springer, 2013), p 356.
		! Griffiths, David J. (2005), Introduction to Quantum Mechanics, 2nd Edition; Pearson Education - Problem 4.12.
		!
		real(dp), intent(in) :: r(:)
		integer , intent(in) :: n,l
		real(dp), allocatable :: y(:)
		real(dp) :: ao, s
		!
		allocate(y(size(r)))
		!
		! y = rho**l * exp(-rho/2) * laguerre(n=n,l=l,x=rho) ! laguerre(p=(n+l),q=(2*l+1),x=rho)
		! 
		! bohr radius
		ao = 1
		! scaling factor
		s = 1/factorial(n - l + 2*l) * sqrt( (2/(ao*n))**3 * factorial(n-l-1)/(2*n*factorial(n+l)))
		! radial wave function
	 	y = s * (2*r/(ao*n)) * exp(-r/(ao*n)) * laguerre(n-l-1, 2*l+1, 2*r/(ao*n))
	 	!
	end function  radial

	function      get_orbitals(nlm,x,y,z) result(orb)
		!
		implicit none
		!
		real(dp), intent(in) :: x(:)
		real(dp), intent(in) :: y(:)
		real(dp), intent(in) :: z(:)
		integer , intent(in) :: nlm(:,:) ! nlm(:,nstates) [n,l,m] for each state
		integer :: n,l,m
		real(dp), allocatable :: orb(:,:)
		real(dp), allocatable :: r(:)
		integer :: i
		integer :: nstates
		!
		nstates=size(nlm,2)
		!
		allocate(r,source=sqrt( x**2 + y**2 + z**2 + 1.0D-14 ))
		!
		allocate(orb(size(x),nstates))
		orb = 0
		!
		do i = 1, nstates
			n = nlm(1,i)
			l = nlm(2,i)
			m = nlm(3,i)
			! psi = (cubic harmonics) * (radial part)
            orb(:,i) = radial(n=n,l=l,r=r)
			! orb(:,i) = Xlc(l=l,m=m,x=x,y=y,z=z,r=r) * Rnl(n=n,l=l,r=r) 
 			! orb(:,i) =  Rnl(n=n,l=l,r=r) 
			!
			! normalize
			orb(:,i) = orb(:,i)/am_dnrm2(orb(:,i))
		enddo
		!
	end function  get_orbitals

	pure function spdf2l(spdf) result(l)
		!
		character(1), intent(in) :: spdf
		integer :: l
		!
		if     (spdf(1:1).eq.'s') then; l = 0
		elseif (spdf(1:1).eq.'p') then; l = 1
		elseif (spdf(1:1).eq.'d') then; l = 2
		elseif (spdf(1:1).eq.'f') then; l = 3
		endif
	end function  spdf2l

	subroutine tight_binding_main(ic,shell,pg,opts)
		!
		implicit none
		!
		class(am_class_irre_cell):: ic
		type(am_class_symmetry)  :: pg ! point group
	    type(am_class_options)   :: opts
	    type(am_class_pair_shell):: shell
    	real(dp),allocatable :: r(:)      ! radial coordinates 
	    real(dp),allocatable :: sph(:,:)  ! real space spherical coordinates sph(1:3,i) = [r,th,phi]
	    real(dp),allocatable :: xyz(:,:)  ! real space cartsian  coordinates xyz(1:3,i) = [x,y,z]
	    real(dp),allocatable :: xyzp(:,:) ! real space cartesian coordinates xyz(1:3,i) = [x,y,z] rotated 
	    real(dp),allocatable :: Y(:,:)    ! Y(r,lm combined indices)
	    integer ,allocatable :: nlms(:,:) ! quantum numbers
	    real(dp),allocatable :: s(:)
	    real(dp) :: Rot(3,3)
	    integer :: i, j, k
	    integer :: N ! spherical coordinates partitions
		!
		! pair_prototypes = shell%get_pair_prototypes(ic=ic,pg=pg,opts=opts)
		!
		if (opts%verbosity.ge.1) call am_print_title('Determining tight-binding parameters')
		!
		! setup real space coordinates
		! Note: d_z^2 requires a very large range in order to be orthogonal with s orbital
		!
        ! something is wrong with spherical integration...
        !
		! setup spherical coordinates
! 		N = 2**4
!         write(*,*) N
!         ! r, 0 <= theta <= pi , 0 <= phi < 2*pi
! 		sph = meshgrid(linspace(0.0_dp,100000.0_dp,3), pi*[0:n]/real(n,dp), twopi*[0:(n-1)]/real(n,dp))
! 		! convert to cartesian coordinates
! 		allocate(xyz,mold=sph)
! 		xyz(1,:)=sph(1,:)*cos(sph(2,:))*sin(sph(3,:))
! 		xyz(2,:)=sph(1,:)*sin(sph(2,:))*sin(sph(3,:))
! 		xyz(3,:)=sph(1,:)*cos(sph(3,:))
        !
       	N = 100
  		xyz = meshgrid(linspace(-N,N,N),linspace(-N,N,N),linspace(-N,N,N))
		!
		if (opts%verbosity.ge.1) call am_print('number of irreducible atoms',ic%natoms)
		!
		!
		i = 1
		nlms = get_states(Z=ic%Z(i),nvalences=10,nexciteds=4)
		!
		! get cubic harmonics
		! get_orbitals(nlm,x,y,z)
		Y = get_orbitals(x=xyz(1,:),y=xyz(2,:),z=xyz(3,:),nlm=unique(nlms(1:3,:)))
		!
		call am_print('orbitals',unique(nlms(1:3,:)))
		! call am_print('Y',Y)
		call am_print('Y^T*Y',matmul(transpose(Y),Y))
		!
        
        
        stop
        
        
        
		! confirm that states are orthogonal
		if (any( abs(am_dgemm(A=Y,B=Y,flags='AT')-eye(k)).gt.tiny )) then
			call am_print('ERROR','Cubic harmonics not orthonormal.')
			call am_print('Y^T*Y', am_dgemm(A=Y,B=Y,flags='AT') )
			call am_print('log(abs(Y^T*Y-I))', log(abs(am_dgemm(A=Y,B=Y,flags='AT')-eye(k))) )
			stop
		endif

! 		! get atomic orbitals of intrest (eventually this will be an input)
! 		! nvo = number of valence orbital
! 		! neo = number of excited orbitals
! 		allocate(ic%atom(ic%natoms))
! 		do i = 1, ic%natoms
! 			call ic%atom(i)%generate_orbitals(Z=ic%Z(i),opts=opts,nvalences=4,nexciteds=4)
! 			! last number in nlms array is the state index starting with H 1s
! 		enddo


	

! 		do i = 1, 1 ! ic%natoms
! 			lm = unique(ic%atom(i)%)
		
! 		!
! 		! quantize each orbital with respect to bond axes
! 		allocate(r_rot,mold=r)
! 		do i = 1, 1  !, ic%natoms
! 			! get rotation matrix which aligns Z axis parallel to bond
! 			Rot = rotmat(a=real([0,0,1],dp),b=[1,1,1]/real(4,dp))
! 			r_rot = matmul(Rot,r)

! 		enddo
! 		stop





	end subroutine tight_binding_main










end module am_tight_binding











