module am_tight_binding

    use am_constants
    use am_helpers
    use am_matlab
    use am_options
    use am_shells
    use am_irre_cell
    use am_atom
    use am_symmetry

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

	pure function two_center_integral(l1,l2) result(sk)
		!
		! Numerical alternative to slater_koster
		!
		implicit none
		!
		integer, intent(in) :: l1
		integer, intent(in) :: l2
		real(dp), allocatable :: sk(:)
		!

	end function  two_center_integral

	subroutine tight_binding_main(ic,shell,pg,opts)
		!
		implicit none
		!
		class(am_class_irre_cell):: ic
		type(am_class_symmetry)  :: pg ! point group
	    type(am_class_options)   :: opts
	    type(am_class_pair_shell):: shell
	    integer, allocatable :: nvs(:)    ! number of valence states per atom
	    integer, allocatable :: nes(:)    ! number of excited states per atom
	    integer, allocatable :: nlms(:,:) ! quantum numbers
	    integer :: i, k
		!
		! pair_prototypes = shell%get_pair_prototypes(ic=ic,pg=pg,opts=opts)
		!
		if (opts%verbosity.ge.1) call am_print_title('Determining tight-binding parameters')
		!
		if (opts%verbosity.ge.1) call am_print('number of irreducible atoms',ic%natoms)
		!
		! these numbers do not take spin into account!
		allocate(nvs(ic%natoms))
		nvs=8
		allocate(nes(ic%natoms))
		nes=8
		!
		do i = 1, ic%natoms
			nlms = atm_state(Z=ic%Z(i),nvs=nvs(i),nes=nes(i),verbosity=opts%verbosity)
			! last number in nlms array is the state index starting with H 1s
			call am_print('nlms',transpose(nlms))
			!
			if (opts%verbosity.ge.1) then
				write(*,'(5x,a)',advance='no') 'atom '//trim(int2char(i))//', '//trim(atm_symb(ic%Z(i)))//': '
				do k = 1, (nvs(i)+nes(i))
					! quantum number n
					write(*,'(i1)',advance='no') nlms(1,k)
					select case(nlms(2,k))
					case(0); write(*,'(a)',advance='no') 's'
					case(1); write(*,'(a)',advance='no') 'p'
					case(2); write(*,'(a)',advance='no') 'd'
					case(3); write(*,'(a)',advance='no') 'f'
					case(4); write(*,'(a)',advance='no') 'g'
					end select
					write(*,'(a)',advance='no') trim(int2char( count( (nlms(1,:).eq.nlms(1,k)) &
														        .and. (nlms(2,:).eq.nlms(2,k)) &
														        .and. (nlms(5,:).le.ic%Z(i)) )))//' '
				enddo
				write(*,*)
			endif
			!
		enddo


	end subroutine tight_binding_main










end module am_tight_binding











