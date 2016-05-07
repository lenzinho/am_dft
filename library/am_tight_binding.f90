module am_tight_binding

    use am_constants
    use am_stdout
    use am_matlab
    use am_options
    use am_shells
    use am_irre_cell
    use am_symmetry
    use am_mkl

	implicit none

	private
    
    public :: test_rotation

	type, public :: am_class_tb
		!
		! matrix elements
		integer :: nVsks ! number of irreducible matrix elements
		real(dp), allocatable :: Vsk(:) ! irreducible matrix elements
		integer , allocatable :: Vsk_label(:,:) ! irreducible matrix element labels [i,j,li,lj,m]
        ! i irreducible cell atom, jth neighbor shell around i, li orbital on i, lj orbital on j, m overlap
		!
	contains
		procedure :: template_matrix_elements
		procedure :: print_matrix_elements ! to stdout
        procedure :: write_matrix_elements ! to file
		procedure :: read_matrix_elements  ! from fil
        
	end type am_class_tb

contains

	subroutine     template_matrix_elements(tb,irrpair,ic,opts)
		!
    	! useful if no states on atoms are defined. Gives each atom all the possible states
    	! correspoding to the L shell (n=2) [ s p   states (1+3=4   states) ]
        ! correspoding to the M shell (n=3) [ s p d states (1+3+5=9 states) ]
    	!
        ! there will be 1 matrix element for every irreducible pair of orbitals
		!
		implicit none
        !
        class(am_class_tb)        , intent(inout) :: tb
        class(am_class_pair_shell), intent(inout) :: irrpair
        class(am_class_irre_cell) , intent(inout) :: ic
        type(am_class_options)    , intent(in) :: opts
        integer, allocatable :: Vsk_label(:,:)
        integer, allocatable :: Vsk_label_unique(:,:)
        logical :: is_spin_polarized
    	integer :: nVsks
    	integer :: n,l,m,s
    	integer :: i,j,k,ii,jj,kk
        integer, allocatable :: nlsi(:,:)
        integer, allocatable :: nlsj(:,:)
        !
        is_spin_polarized = .false.
        !
        if (opts%verbosity.ge.1) call am_print_title('Creating orbital basis functions from template')
        !
    	! create basis on each irreducible atoms
		allocate(ic%atom(ic%natoms))
		do i = 1, ic%natoms
			n = 2
			! number of orbitals
            ic%atom(i)%norbitals = n**2
            ! factor of two here accounts for spin
            if (is_spin_polarized) ic%atom(i)%norbitals = 2*ic%atom(i)%norbitals
            !
			allocate(ic%atom(i)%orbital(4, ic%atom(i)%norbitals))
			k=0
			do l = 0, (n-1)
			do m = -l, l
                if (is_spin_polarized) then
			        do s = -1,1,2
				        k = k+1
				        ic%atom(i)%orbital(1,k) = n
				        ic%atom(i)%orbital(2,k) = l
				        ic%atom(i)%orbital(3,k) = m
				        ic%atom(i)%orbital(4,k) = s
                    enddo
                else
				        k = k+1
				        ic%atom(i)%orbital(1,k) = n
				        ic%atom(i)%orbital(2,k) = l
				        ic%atom(i)%orbital(3,k) = m
				        ic%atom(i)%orbital(4,k) = 0
                endif
			enddo
			enddo
		enddo
		if (opts%verbosity.ge.1) call print_oribtal_basis(irrpair,ic)
        !
        ! quick qestimate to determine how much space to allocate for irreducible matrix elements
        nVsks = 0
        do k = 1, irrpair%nshells
            i = irrpair%shell(k)%i
            j = irrpair%shell(k)%j
            nlsi = unique(ic%atom(i)%orbital([1,2,4],:))
            nlsj = unique(ic%atom(j)%orbital([1,2,4],:))
            do ii = 1, size(nlsi,2)
            do jj = 1, size(nlsj,2)
                  nVsks = nVsks + abs(min( nlsi(2,ii), nlsj(2,jj) )) + 1
            enddo
            enddo
        enddo
		!
    	! now, really determine the number of irreducible matrix elements
        ! to determine number of slater koster matrix elements per irreducible pairs, use the facts that:
        ! 1) (l,l',m) = (l,l',-m)
        ! 2) (l,l',m) = (-1)^(l-m) * (l',l,m)
        allocate(Vsk_label(6,nVsks))
        kk=0
        do k = 1, irrpair%nshells
            ! j irreducible atom shell entered on atom irreducible atom i
            i = irrpair%shell(k)%i
            j = irrpair%shell(k)%j
            !
            nlsi = unique(ic%atom(i)%orbital([1,2,4],:))
            nlsj = unique(ic%atom(j)%orbital([1,2,4],:))
            !
            ! loop over all combination of orbitals (outer product of li x lj)
            do ii = 1, size(nlsi,2)
            do jj = 1, size(nlsj,2)
                ! because of cylindrical symmetry, only the matrix elements for which m = m' remain in each li and lj 
                ! also, -m = m, so the irreducible tight binding parameters include only m = 0 through min(li,lj)
                do m = 0, abs(min( nlsi(2,ii) , nlsj(2,jj) ))
                    !
                    kk=kk+1
                    Vsk_label(1,kk) = k
                    ! sites
                    Vsk_label(2,kk) = irrpair%shell(k)%i
                    Vsk_label(3,kk) = irrpair%shell(k)%j
                    ! azimuthal quantum number
                    Vsk_label(4,kk) = nlsi(2,ii)
                    Vsk_label(5,kk) = nlsj(2,jj)
                    ! (flip the indices to get sp ~= ps )
                    if (Vsk_label(4,kk).gt.Vsk_label(5,kk)) Vsk_label([4,5],kk) = Vsk_label([5,4],kk)
                    ! magnetic quantum number
                    Vsk_label(6,kk) = m
                enddo
            enddo
            enddo
        enddo
        !
        ! get unique (irreducible matrix elements)
        Vsk_label_unique = unique(Vsk_label)
        !
        tb%nVsks = size(Vsk_label_unique,2)
        !
        allocate(tb%Vsk(nVsks))
        allocate(tb%Vsk_label(6,nVsks))
        do k = 1, tb%nVsks
            ! copy irreducible labels
            tb%Vsk_label(:,k) = Vsk_label_unique(:,k)
            ! add a placeholder slater koster parameter to initialize
            tb%Vsk(k) = 1.0_dp
        enddo
        !
		if (opts%verbosity.ge.1) call tb%print_matrix_elements
        !
        contains
	    subroutine     print_oribtal_basis(irrpair,ic)
		    !
		    implicit none
            !
            type(am_class_pair_shell), intent(inout) :: irrpair
            type(am_class_irre_cell) , intent(inout) :: ic
    	    integer :: n,l,m,s
    	    integer :: i,j
            !
		    do i = 1, ic%natoms
			    write(*,'(a5,a)',advance='no') ' ... ', 'irreducible atom '//trim(int2char(i))//' contributes '//trim(int2char(ic%atom(i)%norbitals))//' orbitals (n,l,m,s):'
			    !
			    do j = 1, ic%atom(i)%norbitals
				    if (mod(j,6).eq.1) then
					    write(*,*)
					    write(*,'(5x)',advance='no')
				    endif
				    n = ic%atom(i)%orbital(1,j)
				    l = ic%atom(i)%orbital(2,j)
				    m = ic%atom(i)%orbital(3,j)
				    s = ic%atom(i)%orbital(4,j)
				    write(*,'(a2,i2,a1,i2,a1,i2,a1,i2,a2)',advance='no') ' (', n, ',', l, ',', m, ',', s, ') '
			    enddo
			    write(*,*)
		    enddo
	    end subroutine print_oribtal_basis
    end subroutine template_matrix_elements

    subroutine     print_matrix_elements(tb)
        !
        use am_atom, only : l2spdf, m2spdf
        !
        implicit none
        !
        class(am_class_tb), intent(in) :: tb
        integer :: k
        !
        ! call am_print_title('Tight binding irreducible matrix elements')
        !
        call am_print('irreducible matrix elements',tb%nVsks)
        !
        write(*,'(a5,7a6,a10)') '     ', '#', 'shell', 'i', 'j', 'li', 'lj', 'm', 'Vsk'
        write(*,'(a5,a42,a10)') '     ', repeat(' -----',7), ' '//repeat('-',9)
        !
        do k = 1, tb%nVsks
            write(*,'(5x, i6, 3i6, 3a6, f10.3)') &
                k, &
                tb%Vsk_label(1,k), &               ! shell
                tb%Vsk_label(2,k), &               ! site i
                tb%Vsk_label(3,k), &               ! site j
                trim(l2spdf(tb%Vsk_label(4,k))), & ! li (s,p,d,f)
                trim(l2spdf(tb%Vsk_label(5,k))), & ! lj (s,p,d,f)
                trim(m2spdf(tb%Vsk_label(6,k))), & ! m  (sigma, pi, delta, phi)
                tb%Vsk(k)                          ! matrix element
        enddo
        !
    end subroutine print_matrix_elements

    subroutine     write_matrix_elements(tb)
    	!
    	implicit none
    	!
    	class(am_class_tb), intent(inout) :: tb
    	integer :: fid, k
    	!
        fid = 1
        open(unit=fid,file='outfile.tightbinding',status='replace',action='write')
	        !
	        write(fid,'(i6)') tb%nVsks
            do k = 1, tb%nVsks
            write(fid,'(3i6, 3i6, f)') &
                tb%Vsk_label(1,k), & ! shell
                tb%Vsk_label(2,k), & ! site i
                tb%Vsk_label(3,k), & ! site j
                tb%Vsk_label(4,k), & ! li (s,p,d,f)
                tb%Vsk_label(5,k), & ! lj (s,p,d,f)
                tb%Vsk_label(6,k), & ! m  (sigma, pi, delta, phi)
                tb%Vsk(k) ! matrix element
            enddo
			!
		close(fid)
		!
	end subroutine write_matrix_elements

    subroutine     read_matrix_elements(tb,opts)
    	!
    	implicit none
    	!
    	class(am_class_tb), intent(inout) :: tb
        type(am_class_options), intent(in) :: opts
    	integer :: fid, k, i
    	real(dp) :: b
        character(maximum_buffer_size) :: buffer ! read buffer
        character(len=:), allocatable :: word(:) ! read buffer
    	!
    	if (opts%verbosity.ge.1) call am_print_title('Reading tight binding matrix elements')
    	!
    	if (opts%verbosity.ge.1) call am_print('input file',trim(opts%tbf))
    	!
        fid = 1
        open(unit=fid,file=trim(opts%tbf),status="old",action='read')
	        !
            read(unit=fid,fmt='(a)') buffer
            word = strsplit(buffer,delimiter=' ')
            read(word(1),*) tb%nVsks
	        !
	        if (allocated(tb%Vsk_label)) deallocate(tb%Vsk_label)
	        if (allocated(tb%Vsk)) 		 deallocate(tb%Vsk)
	        !
	        allocate(tb%Vsk_label(6,tb%nVsks))
	        allocate(tb%Vsk(tb%nVsks))
            !
            do k = 1, tb%nVsks
                read(unit=fid,fmt='(a)') buffer
                word = strsplit(buffer,delimiter=' ')
            
                read(word(1),*) tb%Vsk_label(1,k) ! shell
                read(word(2),*) tb%Vsk_label(2,k) ! site i
                read(word(3),*) tb%Vsk_label(3,k) ! site j
                read(word(4),*) tb%Vsk_label(4,k) ! li (s,p,d,f)
                read(word(5),*) tb%Vsk_label(5,k) ! lj (s,p,d,f)
                read(word(6),*) tb%Vsk_label(6,k) ! m  (sigma, pi, delta, phi)
                read(word(7),*) tb%Vsk(k)         ! matrix element
                !
            enddo
            !
            call tb%print_matrix_elements
            !
            !
        close(fid)
		!
	end subroutine read_matrix_elements

	subroutine     get_Ly_in_Lz_basis(l, V, D)
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
		complex(dp), intent(out), allocatable :: V(:,:)  ! My is returned as a (l**2+1) x (l**2+1) matrix containing all possible values of m : |m| .le. l
		real(dp)   , intent(out), allocatable :: D(:)    ! because Ly operator is Hermitian, eigenvalues are real
		complex(dp), allocatable :: Ly(:,:)
		integer :: k, nstates, m
		! 
		! number of states (m=-l,-l+1,-l+2,...,0,...,l-2,l-1,l)
		nstates = l*2 + 1
		!
		! indices of V refer to value of m
		allocate(V(nstates,nstates))
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
		! Diagonalize the Ly matrix
		call am_zhbev(A_b=Ly,V=V,D=D)
		!
		! Equivalent to above, but using full matrix instead
		! complex(dp), allocatable :: Lyf(:,:)
		! Lyf = diag(Ly(2,1:k),-1)
		! Lyf = Lyf + adjoint(Lyf)
		! call am_print('Lyf',Lyf)
		! call am_zheev(A=Lyf,V=My,D=D)
		! call am_print('zheev : V',My)
		! call am_print('zheev : D',D)
		!
    end subroutine get_Ly_in_Lz_basis

    function       spherical_to_tesseral(l) result(B)
        ! R. R. Sharma, Phys. Rev. B. 19, 2813 (1979).
        !     [1    0    0    0     0     0    0    0     0  ]
        !     [0   0.7   0   0.7i   0     0    0    0     0  ]
        !     [0    0    1    0     0     0    0    0     0  ]
        !     [0  -0.7   0   0.7i   0     0    0    0     0  ]
        ! B = [0    0    0    0    0.7    0    0    0    0.7i]
        !     [0    0    0    0     0    0.7   0   0.7i   0  ]
        !     [0    0    0    0     0     0    1    0     0  ]
        !     [0    0    0    0     0   -0.7   0   0.7i   0  ]
        !     [0    0    0    0    0.7    0    0    0   -0.7i]
        implicit none
        !
		integer, intent(in) :: l ! l is the orbital quantum number
		integer  :: nstates
        real(dp) :: invsqrt2
        complex(dp), allocatable :: B(:,:)
        integer :: m,mp
        !
        invsqrt2 = 0.70710678118655_dp
        !
        nstates = l*2 + 1
        allocate(B(-l:l,-l:l))
        B = 0.0_dp
        !
		do m = -l, l
        do mp= -l, l
            ! Real orbitals which are given by Y_{l0}.
            if (m.eq.mp) then
            if (m.eq.0 ) then
                B(m,mp) = 1.0_dp
            endif
            endif
            ! Real orbitals which are given by (Y_{l-m} + (-1)^m Y_{lm}^\ast)/sqrt(2)
			if (mp.lt.0) then
                if (mp.eq.-m) then
                    B(m,mp) = CMPLX(invsqrt2, 0.0_dp)
                endif
                if (mp.eq.m) then
                    B(m,mp) = CMPLX(((-1)**m)*invsqrt2, 0.0_dp)
                endif
            endif
            ! Real orbitals which are given by i*(Y_{l-m} - (-1)^m Y_{lm}^\ast)/sqrt(2)
			if (mp.gt.0) then
                if (mp.eq.-m) then
                    B(m,mp) = CMPLX(0.0_dp, invsqrt2)
                endif
                if (mp.eq.m) then
                    B(m,mp) = CMPLX(0.0_dp, -((-1)**m)*invsqrt2)
                endif
            endif
        enddo
        enddo
        !
    end function   spherical_to_tesseral
    
 	subroutine     test_rotation
 		!
 		implicit none
 		!
        integer :: l
 		complex(dp), allocatable :: V(:,:)
 		real(dp)   , allocatable :: D(:)
        real(dp)   :: th, phi, dcosines(3), vec(3)
        complex(dp), allocatable :: th_mat(:,:)
        complex(dp), allocatable :: phi_mat(:,:)
        complex(dp), allocatable :: U(:,:)
        complex(dp), allocatable :: B(:,:)
        real(dp)   , allocatable :: R(:,:) !  the rotation matrix
        real(dp) :: aa(4)
        !
        !
        l = 1
        !
        ! 0 s, 1 p, 3 d, 4 f
        call get_Ly_in_Lz_basis(l=l, V=V, D=D)
        !
        call am_print('V',V)
        call am_print('D',D)
        !
        th = 30.0*pi/180.0_dp
        phi= 00.0*pi/180.0_dp
 		call am_print('th',th)
        call am_print('phi',phi)
        !
        th_mat  = diag(exp(cmplx_i*th*D))
        phi_mat = diag(exp(cmplx_i*phi*D))
 		call am_print('th_mat',th_mat)
        call am_print('phi_mat',phi_mat)
        !
        U = matmul(matmul(V,th_mat),matmul(adjoint(V),phi_mat))
        call am_print('U',U)
        !
        ! check that U is unitary
        call am_print('det(U)',abs(det(U)))
        if (abs(abs(det(U))-1.0_dp).gt.tiny) then
            call am_print('det(U)',abs(det(U)))
            call am_print('ERROR','U is not unitary.')
            stop
        endif
        !
        ! convert from spherical to tesseral harmonics 
        B = spherical_to_tesseral(l)
        call am_print('B',B)
        call am_print('B\dagB',matmul(adjoint(B),B))
        !
        ! check that B is unitary
        call am_print('det(B)',abs(det(B)))
        if (abs(abs(det(B))-1.0_dp).gt.tiny) then
            call am_print('det(B)',abs(det(U)))
            call am_print('ERROR','B is not unitary.')
            stop
        endif
        !
        !
        R = matmul(adjoint(B),matmul(U,B))
        call am_print('R',R)
        call am_print('det(R)',abs(det(R)))
        if (abs(abs(det(U))-1.0_dp).gt.tiny) then
            call am_print('det(R)',abs(det(U)))
            call am_print('ERROR','R is not unitary.')
            stop
        endif
        !
        aa = R2axis_angle(R)
        call am_print('aa',aa)
        
        
        
        !
        ! convert from spherical to tesseral harmonics 
        !
        !
        vec = [0,0,1]
        !
 		if (norm2(vec).gt.tiny) then
 		    dcosines = vec/norm2(vec)
 		else
 		    dcosines = real([0,0,1],dp)
        endif
 		! th is polar angle (measured from z axis)
        ! phi is azimulthal angle (around z axis)
 		th = acos(dcosines(3))
 		phi = atan(dcosines(2)/(dcosines(1)+1.0D-14))
 		
        !
        call am_print('vec',vec)
        call am_print('dcosines',dcosines)
 		call am_print('th',th)
        call am_print('phi',phi)
 		!
        stop
 	end subroutine test_rotation
    
    !subroutine     build_K(tb,uc,K)
    !    !
    !    ! k-independent part of the Hamiltonian
    !    implicit none
    !    !
    !    class(am_class_tb), intent(in) :: tb
    !    class(am_class_pair_shell), intent(in) :: irrpair
    !    real(dp), allocatable, intent(out) :: K(:,:)
    !    integer :: kdim
    !    integer :: k
    !    integer :: i,j,alpha,beta,ii,jj,iuc,juc
    !    !
    !    ! determine size of K
    !    kdim = 0
    !    do i = 1, uc%natoms
    !        kdim = kdim + ic%atom( uc%ic_identifier(i) )%norbitals
    !    enddo
    !    allocate(K(kdim,kdim))
    !    K = 0.0_dp
    !    !
    !    ! build K_{i,alpha,j,beta} i,j are primitive atom indices; alpha, beta are orbital indices
    !    ! in reality, make K a matrix rather than a tensor by using a compound index ii = (i,alpha) and jj = (j,beta), i.e. K_{ii,jj}
    !    ii=0
    !    do iuc = 1, uc%natoms
    !        i = uc%ic_identifier(iuc)
    !        do alpha = 1, ic%atom( i )%norbitals
    !            ii=ii+1
    !            !
    !            !
    !            jj=0
    !            do juc = 1, uc%natoms
    !                j = uc%ic_identifier(juc)
    !                do beta = 1, ic%atom( j )%norbitals
    !                    jj=jj+1
    !                    !
    !                    ! orbitals quantum numbers [n,l,m,s]
    !                    li = ic%atom(i)%orbital(2,alpha)
    !                    lj = ic%atom(j)%orbital(2,beta)
    !                    mi = ic%atom(i)%orbital(3,alpha)
    !                    mj = ic%atom(j)%orbital(3,beta)
    !                    !
    !                    ! this is not right because the index j below is shell
    !                    ! here, it is being used as a primitive cell index.
    !                    ! NOT CORRECT!
    !                    ! NOT CORRECT!
    !                    ! NOT CORRECT!
    !                    do k = 1, tb%nVsks
    !                        if (mi.eq.mj) then
    !                        if (i .eq.Vsk_label(1,k)) then
    !                        if (j .eq.Vsk_label(2,k)) then
    !                        if (li.eq.Vsk_label(3,k)) then
    !                        if (lj.eq.Vsk_label(4,k)) then
    !                        if (m .eq.Vsk_label(5,k)) then
    !                        K(ii,jj) = Vsk(k)
    !                        endif
    !                        endif
    !                        endif
    !                        endif
    !                        endif
    !                        endif
    !                    enddo
    !                    ! Vsk_label(:,:) ! irreducible matrix element labels [i,j,li,lj,m]
    !                   !  K(ii,jj) =
    !
    !                enddo
    !            enddo
    !        enddo
    !    enddo
    !    
    !    
    !        
    !        
    !        !i  = tb%Vsk_label(1,k)
    !        !j  = tb%Vsk_label(2,k)
    !        !li = tb%Vsk_label(3,k)
    !        !lj = tb%Vsk_label(4,k)
    !        !m  = tb%Vsk_label(5,k)
    !
    !    enddo
    !    !
    !
    !end subroutine build_K


	!
	! stuff below is not in use
	!

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











