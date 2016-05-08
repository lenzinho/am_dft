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
    
    public :: test_SO3

	type, public :: am_class_tb
		!
		! matrix elements
		integer :: nVsks ! number of irreducible matrix elements
		real(dp), allocatable :: Vsk(:) ! irreducible matrix elements
		integer , allocatable :: Vsk_label(:,:) ! irreducible matrix element labels [k,i,j,li,lj,m]
        ! i irreducible cell atom, jth neighbor shell around i, li orbital on i, lj orbital on j, m overlap
		!
	contains
		procedure :: template_matrix_elements
		procedure :: print_matrix_elements ! to stdout
        procedure :: write_matrix_elements ! to file
		procedure :: read_matrix_elements  ! from file
        
	end type am_class_tb

contains

	subroutine     template_matrix_elements(tb,ip,ic,opts)
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
        class(am_class_pair_shell), intent(inout) :: ip
        class(am_class_irre_cell) , intent(inout) :: ic
        type(am_class_options)    , intent(in) :: opts
        integer, allocatable :: Vsk_label(:,:)
        integer, allocatable :: Vsk_label_unique(:,:)
        logical :: is_spin_polarized
    	integer :: nVsks
        integer, allocatable :: slist(:)
        integer :: smax
    	integer :: n,l,m,s
    	integer :: i,j,k,ii,jj,kk
        integer, allocatable :: nlsi(:,:)
        integer, allocatable :: nlsj(:,:)
        !
        is_spin_polarized = .false.
        if (is_spin_polarized) then
            allocate(slist,source=[-1,1])
            smax = size(slist)
        else
            allocate(slist,source=[0])
            smax = size(slist)
        endif
        !
        if (opts%verbosity.ge.1) call am_print_title('Creating orbital basis functions from template')
        !
    	! create basis on each irreducible atoms
		allocate(ic%atom(ic%natoms))
		do i = 1, ic%natoms
			n = 2
			! number of orbitals
            ic%atom(i)%norbitals = n**2
            ! multiply by two for spin polarized case
            if (is_spin_polarized) ic%atom(i)%norbitals = 2*ic%atom(i)%norbitals
            ! allocate space for orbitals
			allocate(ic%atom(i)%orbital(4, ic%atom(i)%norbitals))
            ! generate orbitals
			k=0
			do l = 0, (n-1)
			do m = -l, l
		        do s = 1, smax
			        k = k+1
			        ic%atom(i)%orbital(1,k) = n
			        ic%atom(i)%orbital(2,k) = l
			        ic%atom(i)%orbital(3,k) = m
			        ic%atom(i)%orbital(4,k) = slist(s)
                enddo
			enddo
			enddo
		enddo
		if (opts%verbosity.ge.1) call print_oribtal_basis(ip,ic)
        !
        ! quick estimate to determine how much space to allocate for irreducible matrix elements
        nVsks = 0
        do k = 1, ip%nshells
            i = ip%shell(k)%i
            j = ip%shell(k)%j
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
        do k = 1, ip%nshells
            ! j irreducible atom shell entered on atom irreducible atom i
            i = ip%shell(k)%i
            j = ip%shell(k)%j
            ! make a list of orbital indices which omits degeneracies in m the magnetic quantum nmber
            nlsi = unique(ic%atom(i)%orbital([1,2,4],:))
            nlsj = unique(ic%atom(j)%orbital([1,2,4],:))
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
                    Vsk_label(2,kk) = ip%shell(k)%i
                    Vsk_label(3,kk) = ip%shell(k)%j
                    ! azimuthal quantum number
                    Vsk_label(4,kk) = nlsi(2,ii)
                    Vsk_label(5,kk) = nlsj(2,jj)
                    ! if irreducible atoms are identical flip the indices to get sp ~= ps
                    if (Vsk_label(2,kk).eq.Vsk_label(3,kk)) then
                    if (Vsk_label(4,kk).gt.Vsk_label(5,kk)) Vsk_label([4,5],kk) = Vsk_label([5,4],kk)
                    endif
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
	    subroutine     print_oribtal_basis(ip,ic)
		    !
		    implicit none
            !
            type(am_class_pair_shell), intent(inout) :: ip
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
        write(*,'(5x,a)') 'Definitions:'
        write(*,'(5x,a)') 'i , j : irreducible atoms indicies'
        write(*,'(5x,a)') 'li,lj : orbitals on irreducible atoms'
        write(*,'(5x,a)') 'm     : type of overlap'
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
        character(maximum_buffer_size) :: buffer ! read buffer
        character(len=:), allocatable :: word(:) ! read buffer
        integer :: fid, k
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

    function       SO3(l,th,phi) result(R)
        !
        ! Definitions:
        ! l       angular quantum number
        ! theta   polar angle (rotate away from Z axis, i.e. around Y axis)
        ! phi     azimuthal angle (rotate around Z axis)
        !
        ! For l = 1 (p orbitals transform like x, y, z), R is equivalent to a rotation around the Y
        ! axis, followed by a rotation around the Z axis. That is,
        !
        !         R = Z*Y with
        !
        !         Y = axis_angle2rot([0.0_dp, 1.0_dp, 0.0_dp, theta_phi(1)])
        !         Z = axis_angle2rot([0.0_dp, 0.0_dp, 1.0_dp, theta_phi(2)])
        !
        ! The operation
        !       
        !         K = Z' * Y * Z
        ! 
        !         takes directional cosines of {th,phi} to [0 0 1]. That is,
        !
        !             [ sin(th)cos(phi) ]   [ 0 ]
        !         K * [ sin(th)sin(phi) ] = [ 0 ]
        !             [     cos(th)     ]   [ 1 ]
        !
        ! Equivalently, using the transpose distributive properties, (A*B)' = B'*A', and the fact that K is unitary, inv(K) = K':
        !
        !         K' = Z' * Y' * Z
        !
        !              [ 0 ]   [ sin(th)cos(phi) ]
        !         K' * [ 0 ] = [ sin(th)sin(phi) ]
        !              [ 1 ]   [     cos(th)     ]
        !
        ! Quick matlab script (note different sign conventions for theta/phi in vrrotvec2mat vs. that used here)
        !
        !     v = [3 2 1]';
        !     dcosines = v/norm(v);
        !     theta_phi(1) = acos(dcosines(3));
        !     theta_phi(2) = atan(dcosines(2)/(dcosines(1)+1E-14));
        !     Y=vrrotvec2mat([0 1 0 -theta_phi(1)]);
        !     Z=vrrotvec2mat([0 0 1 -theta_phi(2)]);
        !     R=(Z*Y)
        !     K=(Z'*Y*Z)
        !     K*dcosines
        !
        implicit none
        !
        integer , intent(in) :: l
        real(dp), intent(in) :: th
        real(dp), intent(in) :: phi
        integer :: n
        real(dp)   , allocatable :: R(:,:) ! tesseral harmonics rotation matrix, (2*l+1) x (2*l+1) irrep of rotation in 3 dimensional space 
        complex(dp), allocatable :: V(:,:) ! eigenvector of Ly in Lz basis
        real(dp)   , allocatable :: D(:)   ! eigenvalues of Ly in Lz basis = -l, -l+1, ... -1, 0, 1, ... l-1, l
        complex(dp), allocatable :: th_mat(:,:)  ! matrix corresponding to exp(-i*theta*Ly)
        complex(dp), allocatable :: phi_mat(:,:) ! matrix corresponding to exp(-i* phi *Ly)
        complex(dp), allocatable :: U(:,:) ! rotation matrix in complex coordinates (spherical harmonics)
        complex(dp), allocatable :: B(:,:) ! similarity transform which maps complex spherical harmonics onto real tesseral harmonics
        !
        ! get dimensions of matrices
        n = 2*l + 1
        !
        ! get Ly eigenvector and eigenvalus in Lz basis
        call get_Ly_in_Lz(l=l, V=V, D=D)
        !
        ! setup rotation matrices
        th_mat  = diag(exp(cmplx_i* th*D)) ! away from Z
        phi_mat = diag(exp(cmplx_i*phi*D)) ! around Z
        !
        ! construct U (rotation in complex coordinates) by first rotating AWAY FROM Z then AROUND Z
        allocate(U(n,n))
        U = matmul(matmul(V,phi_mat),matmul(adjoint(V),th_mat))
        !
        ! check that U is unitary
        if (abs(abs(det(U))-1.0_dp).gt.tiny) then
            call am_print('det(U)',abs(det(U)))
            call am_print('U\dag * U',matmul(adjoint(U),U))
            call am_print('ERROR','U is not unitary.',flags='E')
            stop
        endif
        !
        ! determine similarity transform to convert spherical into tesseral harmonics (complex to real)
        B = spherical2tesseral(l)
        ! B IS INTRODUCING ALOT OF NUMERICAL NOISE. TO DO: FIGURE OUT WHY... LOOK AT DET(B), not very close to 1
        !         call am_print('det(B)',abs(abs(det(B))-1.0_dp))
        !         call am_print('det(B)',abs(abs(det(B))-1.0_dp).gt.tiny)
        !         call am_print('det(B)',tiny)
        !
        ! check that B is unitary
        if ( abs(abs(det(B))-1.0_dp) .gt. tiny ) then
            call am_print('det(B)',abs(det(U)))
            call am_print('ERROR','B is not unitary.',flags='E')
            stop
        endif
        !
        ! convert U to real coordinates
        allocate(R(n,n))
        R = matmul(adjoint(B),matmul(U,B))
        !
        ! check that R is unitary
        if (abs(abs(det(R))-1.0_dp).gt.tiny) then
            call am_print('det(R)',abs(det(R)))
            call am_print('R\dag * R',matmul(transpose(R),R))
            call am_print('ERROR','R is not unitary.',flags='E')
            stop
        endif
        !
        contains
        subroutine     get_Ly_in_Lz(l, V, D)
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
        end subroutine get_Ly_in_Lz
        function       spherical2tesseral(l) result(B)
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
            invsqrt2 = 0.707106781186547_dp
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
        end function   spherical2tesseral
    end function   SO3

 	subroutine     test_SO3
 		!
 		implicit none
 		!
        integer  :: l
        real(dp) :: theta_phi(2), dcosines(3), vec(3)
        real(dp), allocatable :: R(:,:) !  the rotation matrix
        real(dp) :: Z(3,3), Y(3,3)
        !
        vec = real([2.0,2.0,1.0],dp)
        !
        dcosines = vec2dcosines(vec)
        !
        theta_phi = dcosines2thetaphi(dcosines)
        !
        call am_print('vec',vec)
        call am_print('dcosines',dcosines)
        call am_print('th & phi',theta_phi*180.0_dp/pi)
        !
        l = 1
        !
        R = SO3(l=l,th=theta_phi(1),phi=theta_phi(2))
        call am_print('R',R)
        !
        ! R is equivalent to first rotating around Y (away from Z) and then around Z
        Y = axis_angle2rot([0.0_dp, 1.0_dp, 0.0_dp, theta_phi(1)]) ! rot around Y, away from Z
        Z = axis_angle2rot([0.0_dp, 0.0_dp, 1.0_dp, theta_phi(2)]) ! rot around Z, in XY plane
        call am_print('Z*Y',matmul(Z,Y))
        !
        ! should equal zero
        call am_print('(Z*Y-R) * 1E5',(matmul(Z,Y)-R)*1.0D5)
        !
        ! should equal [0,0,1]
        call am_print('[001]', matmul(matmul(transpose(Z),matmul(Y,Z)),dcosines) )
        !
        stop
        !
 	end subroutine test_SO3
    
!     subroutine     build_K(tb,uc,ic,pp,ip,K)
!        !
!        ! k-independent part of the Hamiltonian
!        implicit none
!        !
!        class(am_class_tb), intent(in) :: tb ! tight binding parametrs
!        class(am_class_unit_cell) , intent(in) :: uc ! irreducible cell
!        class(am_class_irre_cell) , intent(in) :: ic ! irreducible cell
!        class(am_class_pair_shell), intent(in) :: pp ! primitive pairs
!        class(am_class_pair_shell), intent(in) :: ip ! irreducible pairs
!        real(dp), allocatable, intent(out) :: K(:,:,:) ! matrix elements for each shell
!        integer :: kdim
!        integer :: k
!        integer :: i,j,alpha,beta,ii,jj,iuc,juc
!        !
!        ! determine size of K
!        kdim = 0
!        do i = 1, uc%natoms
!            kdim = kdim + ic%atom( uc%ic_identifier(i) )%norbitals
!        enddo
!        allocate(K(kdim,kdim))
!        K = 0.0_dp
!        !
!        ! build K_{i,alpha,j,beta} i,j are primitive atom indices; alpha, beta are orbital indices
!        ! in reality, make K a matrix rather than a tensor by using a compound index ii = (i,alpha) and jj = (j,beta), i.e. K_{ii,jj}
!        ii=0
!        do iuc = 1, uc%natoms
!            i = uc%ic_identifier(iuc)
!            do alpha = 1, ic%atom( i )%norbitals
!                ii=ii+1
!                !

!            enddo
!        enddo
       
       
           
           
!        enddo
!        !
    
!     end subroutine build_K


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

!   subroutine     test_orbitals()
!       !
!       implicit none
!       !
!       real(dp),allocatable :: r(:)   ! real space spherical coordinates
!         real(dp),allocatable :: th(:)  ! real space spherical coordinates
!         real(dp),allocatable :: phi(:) ! real space spherical coordinates
!       integer ,allocatable :: state(:,:)
!       real(dp),allocatable :: psi_r(:,:)
!         real(dp),allocatable :: psi_th(:,:)
!         real(dp),allocatable :: psi_phi(:,:)
!         real(dp),allocatable :: overlap(:,:)
!         real(dp),allocatable :: overlap_r(:,:)
!         real(dp),allocatable :: overlap_th(:,:)
!         real(dp),allocatable :: overlap_phi(:,:)
!         complex(dp), allocatable :: My(:,:)
!       integer :: n,l,m
!       integer :: nstates
!       integer :: i
!       !
!       !
!       allocate(state(3,100))
!       nstates=0
!       do n = 1,2
!       do l = 0,n-1
!       do m = -l,l
!           nstates=nstates+1
!           state(:,nstates)=[n,l,m]
!       enddo
!       enddo
!       enddo 
!       call am_print('state',state(:,1:nstates))
!       !
!       ! setup points on which wave function will be evaluated
!       call setup_query_points_for_wavefunction(r=r,th=th,phi=phi,flags='mesh,surf')
!       !
!       ! evalute wave functions
!       allocate(psi_r(size(r),nstates))
!         allocate(psi_th(size(th),nstates))
!         allocate(psi_phi(size(phi),nstates))
!       do i = 1,nstates
!           !
!           n = state(1,i)
!           l = state(2,i)
!           m = state(3,i)
!           !
!           call atomic_orbital(n=n,l=l,m=m,r=r,th=th,phi=phi,psi_r=psi_r(:,i),psi_th=psi_th(:,i),psi_phi=psi_phi(:,i))
!             !
!         enddo
!       !
!       ! calculate overlaps
!         overlap_r=am_dgemm(psi_r,psi_r,flags='AT')
!         overlap_th=am_dgemm(psi_th,psi_th,flags='AT')
!         overlap_phi=am_dgemm(psi_phi,psi_phi,flags='AT')
!       call am_print('r overlap',overlap_r)
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
!   end subroutine test_orbitals







end module am_tight_binding











