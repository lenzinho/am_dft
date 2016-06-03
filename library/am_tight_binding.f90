module am_tight_binding

    use am_constants
    use am_stdout
    use am_matlab
    use am_options
    use am_shells
    use am_irre_cell
    use am_prim_cell
    use am_symmetry
    use am_symmetry_rep
    use am_symmetry_relations
    use am_mkl
    use dispmodule

	implicit none

	private

    public :: test_hamiltonian, set_Vsk_dummies ! for testing

    type, public, extends(am_class_tensor) :: am_class_tens_tb
        ! <INHERITED>
        ! character(100)        :: property       ! name of property
        ! character(100)        :: flags          ! axial/polar/i==j/i/=j
        ! integer               :: rank           ! tensor rank
        ! integer , allocatable :: dims(:)        ! tensor dimensions
        ! real(dp), allocatable :: relations(:,:) ! relations connecting tensor elements
        ! real(dp), allocatable :: V(:)           ! the value of the tensor, use reshape(V,dims)
        ! </INHERITED>
    end type am_class_tens_tb

    type, public :: am_class_tight_binding
        integer :: nshells                     ! how many shells irreducible atoms
        type(am_class_tens_tb), allocatable :: tbvsk(:)  ! tbvsk(nshells)
        integer, allocatable :: nVs            ! number of irreducible matrix element values
        integer, allocatable :: V(:)           ! their values
        integer, allocatable :: V_ind(:,:)     ! their indices V_ind( [i,j], nVs); j = compound [alpha x beta index, see initialize_tb subroutine]
        contains
        procedure :: initialize_tb
        procedure :: get_hamiltonian
        procedure :: write_tb_explicit
        procedure :: export_tb2matlab
    end type am_class_tight_binding

contains

    subroutine     initialize_tb(tb,pg,ip,ic,pc,opts)
        !
        implicit none
        !
        class(am_class_tight_binding), intent(out) :: tb   ! tight binding parameters
        type(am_class_point_group)        , intent(in)  :: pg   ! seitz point group
        type(am_class_prim_cell)          , intent(in)  :: pc   ! primitive cell
        type(am_class_irre_cell)          , intent(in)  :: ic   ! irreducible cell
        type(am_class_irre_pair)          , intent(in)  :: ip   ! irreducible pairs
        type(am_class_options)            , intent(in)  :: opts
        logical, allocatable :: is_independent(:)
        integer :: sub(2)
        integer :: nterms
        integer :: i,j,k
        integer :: a,b,c
        integer :: m,n,o
        !
        !
        if (opts%verbosity.ge.1) call print_title('Symmetry-adapted tight-binding parameters')
        !
        tb%nshells = ip%nshells
        allocate(tb%tbvsk(tb%nshells))
        !
        do k = 1, tb%nshells
            ! get rank, dims, flags, property
            call initialize_tbvsk(tbvsk=tb%tbvsk(k), pc=pc, ic=ic, shell=ip%shell(k))
            ! get relations
            call get_shell_relations(tbvsk=tb%tbvsk(k), pg=pg, ic=ic, shell=ip%shell(k), opts=opts)
        enddo
        ! get number of independent (irreducible) matrix elements
        tb%nVs = 0
        do i = 1, tb%nshells
            tb%nVs = tb%nVs + count(get_independent(tb%tbvsk(i)%relations)) 
        enddo
        ! allocate space for independent (irreducible) matrix elements V
        allocate(tb%V(tb%nVs))
        tb%V = 0
        ! indices : V_ind( [i, subd2ind(dims=tb%tbvsk(i)%dims,sub=[alpha,beta])], nVs)
        allocate(tb%V_ind(3,tb%nVs))
        tb%V_ind = 0
        k=0
        do i = 1, tb%nshells
            is_independent = get_independent(tb%tbvsk(i)%relations)
            do j = 1, product(tb%tbvsk(i)%dims)
            if (is_independent(j)) then
                k=k+1
                sub = ind2sub(dims=tb%tbvsk(i)%dims, ind=j)
                tb%V_ind(1,k) = i      ! shell
                tb%V_ind(2,k) = sub(1) ! alpha
                tb%V_ind(3,k) = sub(2) ! beta
            endif
            enddo
        enddo
        !
        if (opts%verbosity.ge.1) then
            ! print statistics about parameters
            write(*,'(a5,2a6,3a14)') ' ... ', 'shell', 'terms', 'null', 'dependent', 'independent'
            write(*,'(5x,a)') repeat(' '//repeat('-',5),2)//repeat(' '//repeat('-',13),3)
            ! print table
            a=0;b=0;c=0;nterms=0
            do i = 1, tb%nshells
                m = count(get_null(tb%tbvsk(i)%relations))
                n = count(get_depenent(tb%tbvsk(i)%relations)) 
                o = count(get_independent(tb%tbvsk(i)%relations))
                nterms = m+n+o
                write(*,'(5x,i6)'        ,advance='no') i
                write(*,'(i6)'           ,advance='no') nterms
                write(*,'(i5,a2,f5.1,a2)',advance='no') m, '(', (m*100_dp)/real(nterms,dp) , '%)'
                write(*,'(i5,a2,f5.1,a2)',advance='no') n, '(', (n*100_dp)/real(nterms,dp) , '%)'
                write(*,'(i5,a2,f5.1,a2)',advance='no') o, '(', (o*100_dp)/real(nterms,dp) , '%)'
                write(*,*)
                a = a + m
                b = b + n
                c = c + o
            enddo
            ! print total
            write(*,'(5x,a6)'        ,advance='no') 'total'
            write(*,'(i6)'           ,advance='no') (a+b+c)
            write(*,'(i5,a2,f5.1,a2)',advance='no') a, '(', (a*100_dp)/real(a+b+c,dp) , '%)'
            write(*,'(i5,a2,f5.1,a2)',advance='no') b, '(', (b*100_dp)/real(a+b+c,dp) , '%)'
            write(*,'(i5,a2,f5.1,a2)',advance='no') c, '(', (c*100_dp)/real(a+b+c,dp) , '%)'
            write(*,*)
            ! print symmetry relations
            do i = 1, tb%nshells
                write(*,'(a5,a)') ' ... ', 'shell '//trim(int2char(i))//' irreducible symmetry relations:'
                call print_relations(relations=tb%tbvsk(i)%relations, dims=tb%tbvsk(i)%dims, flags='print:dependent,independent')
            enddo
        endif
        !
        contains
        subroutine     initialize_tbvsk(tbvsk,pc,ic,shell)
            !
            implicit none
            !
            type(am_class_tens_tb)  , intent(out) :: tbvsk
            type(am_class_prim_cell), intent(in)  :: pc
            type(am_class_irre_cell), intent(in)  :: ic
            type(am_shell_cell)     , intent(in)  :: shell
            integer :: i
            !
            tbvsk%property = 'tight binding'
            tbvsk%rank  = 2
            ! detetmine if irreducible atoms are the same
            if (shell%i.eq.shell%j) then
                tbvsk%flags = 'i==j'
            else
                tbvsk%flags = 'i/=j'
            endif
            ! set dimensions of matrix elements
            allocate(tbvsk%dims(tbvsk%rank))
            tbvsk%dims = 0 
            ! dimension of Hamiltonian subsection corresponding to primitive atom m (irreducible atoms i)
            do i = 1, ic%atom(shell%i)%nazimuthals
                tbvsk%dims(1) = tbvsk%dims(1) + ic%atom(shell%i)%azimuthal(i)*2+1
            enddo
            ! dimension of Hamiltonian subsection corresponding to primitive atom n (irreducible atoms j)
            do i = 1, ic%atom(shell%j)%nazimuthals
                tbvsk%dims(2) = tbvsk%dims(2) + ic%atom(shell%j)%azimuthal(i)*2+1
            enddo
            ! allocate space for matrix elements
            allocate(tbvsk%V(tbvsk%dims(1)*tbvsk%dims(2)))
            !
        end subroutine initialize_tbvsk
        subroutine     get_shell_relations(tbvsk,pg,ic,shell,opts)
                !
                implicit none
                !
                type(am_class_tens_tb)    , intent(inout) :: tbvsk
                type(am_class_point_group), intent(in) :: pg   ! seitz point group
                type(am_class_irre_cell)  , intent(in) :: ic
                type(am_shell_cell)       , intent(in) :: shell
                type(am_class_options)    , intent(in) :: opts
                type(am_class_flat_group)  :: flat_pg
                type(am_class_flat_group)  :: flat_ig
                type(am_class_point_group) :: stab
                !
                ! determine intrinsic symmetries (if both irreducible atoms are of the same irreducible type)
                ! Interchange of indices: (l,l',m) = (-1)^(l+l') * (l',l,m), due to parity of wavefunction. (s,d are even under inversion, p,f are odd)
                ! E. Scheer, Molecular Electronics: An Introduction to Theory and Experiment, p 245. Also see R. Martin.
                ! NOTE: FOR SOME REASON, MUST CALL atom_m with shell%j and atom_n with shell%i (INDICES FLIPPED) - probably has to do with reshape and column-major ordering
                call flat_ig%get_flat_intrinsic_group(tens=tbvsk, atom_m=ic%atom(shell%j), atom_n=ic%atom(shell%i) )
                ! determine stabilizers relations
                call stab%get_stabilizer_group(pg=pg, v=shell%tau_cart(1:3,1), opts=opts, flags='cart')
                ! get stabilizer symmetries in the flattened hamiltonin basis
                call flat_pg%get_flat_point_group(tens=tbvsk, pg=stab, atom_m=ic%atom(shell%j), atom_n=ic%atom(shell%i))
                ! get combined relations
                tbvsk%relations = combine_relations(relationsA=flat_pg%relations, relationsB=flat_ig%relations)
                ! 
        end subroutine get_shell_relations
    end subroutine initialize_tb

    function       get_hamiltonian(tb,tbpg,ic,ip,pp,kpt,iopt_mask) result(H)
        ! 
        ! Get tight binding Hamiltonian at kpt.
        ! The only reason ip is included is for the mask.
        ! 
        implicit none
        !
        class(am_class_tight_binding), intent(in) :: tb     ! tight binding matrix elements
        type(am_class_tb_group)      , intent(in) :: tbpg   ! point group in tight binding representation
        type(am_class_irre_cell)     , intent(in) :: ic     ! irreducible cell
        type(am_class_prim_pair)     , intent(in) :: pp     ! primitive pairs
        type(am_class_irre_pair)     , intent(in) :: ip     ! irreducible pairs
        real(dp)                     , intent(in) :: kpt(3) ! fractional
        logical, optional            , intent(in) :: iopt_mask(:)
        logical    , allocatable :: mask(:)
        integer    , allocatable :: S(:),E(:)
        complex(dp), allocatable, target :: Hsub_target(:,:)
        real(dp)   , allocatable, target :: pg_target(:,:,:)
        complex(dp), allocatable :: H(:,:)
        complex(dp), pointer :: Hsub(:,:)
        real(dp), pointer :: Dm(:,:), Dn(:,:)
        integer :: m ! primitive atom 1 index 
        integer :: n ! primitive atom 2 index
        integer :: i ! irreducible atom 1 index 
        integer :: j ! irreducible atom 2 index
        integer :: k ! shell (primitive)
        integer :: l ! shell (irreducible)
        integer :: p ! atoms
        !
        ! mask irreducible pair shells 
        !      -----------
        if (present(iopt_mask)) then
            allocate(mask,source=iopt_mask)
        else
            allocate(mask(ip%nshells))
            mask = .true.
        endif
        !
        ! allocate space for vectors demarking start and end of Hamiltonian subsection
        allocate(S, source=tbpg%H_start)
        allocate(E, source=tbpg%H_end)
        ! allocate workspace for H subsection (initialized later)
        allocate(Hsub_target(tbpg%nbases,tbpg%nbases))
        ! allocate workspace for H subsection (initialized later)
        allocate(pg_target, source=tbpg%sym)
        ! allocate and initialize
        allocate(H(tbpg%nbases,tbpg%nbases))
        H = cmplx(0,0,dp)
        ! construct Hamiltonian
        do l = 1, ip%nshells
        if ( mask(l) ) then
            do k = 1, pp%nshells
            if (abs(pp%ip_id(k)).eq.l) then
                ! primitive atom indicies
                m = pp%shell(k)%m
                n = pp%shell(k)%n
                ! irreducible atom indicies
                i = pp%shell(k)%i
                j = pp%shell(k)%j
                ! compute bloch sum by loop over atoms in shell
                do p = 1, pp%shell(k)%natoms
                    ! set pointers
                    Hsub => Hsub_target(S(m):E(m), S(n):E(n))
                    Dm   => pg_target(S(m):E(m), S(m):E(m), pp%shell(k)%pg_id(p) )
                    Dn   => pg_target(S(n):E(n), S(n):E(n), pp%shell(k)%pg_id(p) )
                    ! get matrix elements (initialize Hsub)
                    Hsub = get_Vsk(tb=tb, ip_id=pp%ip_id(k), atom_m=ic%atom(i), atom_n=ic%atom(j))
                    ! rotate matrix elements as needed to get from the tau_frac(:,1) => tau_frac(:,x)
                    Hsub = matmul(matmul(transpose(Dm), Hsub), Dn)
                    ! multiply exponential factor from Bloch sum
                    Hsub = Hsub * exp(-itwopi*dot_product(pp%shell(k)%tau_cart(1:3,p), kpt)) ! kpt [cart]
                    ! this pair's contribution to the Hamiltonian
                    H(S(m):E(m), S(n):E(n)) = H(S(m):E(m), S(n):E(n)) + Hsub
                enddo
            endif
            enddo
        endif
        enddo
    end function   get_hamiltonian

    subroutine     test_hamiltonian(tb,tbpg,ic,ip,pp)
        ! makes sure Hamiltonian at Gamma commutes with all point symmetry operations
        implicit none
        !
        type(am_class_tight_binding), intent(in) :: tb   ! tight binding matrix elements
        type(am_class_tb_group)          , intent(in) :: tbpg ! point group in tight binding representation
        type(am_class_irre_cell)         , intent(in) :: ic   ! irreducible cell
        type(am_class_prim_pair)         , intent(in) :: pp   ! primitive pairs
        type(am_class_irre_pair)         , intent(in) :: ip   ! irreducible pairs
        complex(dp), allocatable :: H(:,:)
        real(dp)   , allocatable :: R(:,:)
        logical    , allocatable :: mask(:) ! mask the irreducible shells which are considered in construction of the hamiltonin
        integer :: i, j
        !
        ! allocate space for mask
        allocate(mask(ip%nshells))
        ! allocate space for symmetry in tb basis
        allocate(R(tbpg%nbases,tbpg%nbases))
        !
        call print_title('Checking Hamiltonian at Gamma')
        !
        ! loop over irreducible shells
        do j = 2,2 !1, ip%nshells
            ! ignore all irreducible shells, except j
            mask    = .false.
            mask(j) = .true.
            !
            H = tb%get_hamiltonian(tbpg=tbpg, ic=ic, ip=ip, pp=pp, kpt=real([0,0,0],dp),iopt_mask=mask)
            !
            ! check that H is hermitian
            if ( .not. isequal(H,adjoint(H)) ) then
                call disp('H',H)
                stop 'H is not Hermitian.'
            endif
            !
            ! check that H commutes with all point symmetries
            do i = 1, tbpg%nsyms
                R = tbpg%sym(:,:,i)
                if (.not.isequal(matmul(H,R),matmul(R,H))) then
                    call disp('H',H,style='above')
                    call disp('R',R,style='above')
                    call disp('[H,R]',matmul(H,R)-matmul(R,H),style='above')
                    call disp("R'*H*R - H",matmul(matmul(transpose(R),H),R)-H,style='above')
                    stop 'H does not commute with symmetry.'
                endif
            enddo
            !
            write(*,'(a5,a)') ' ... ' ,'shell '//tostring(j)//': OK!'
        enddo
    end subroutine test_hamiltonian

    function       get_Vsk(tb,atom_m,atom_n,ip_id) result(V)
        ! irreducible pair (ip_id), can be negative (corresponds to pair n-m rathet than m-n, on the
        ! opposite [upper/lower] side of the Hamiltonian), in which case adjoint of V is returned
        use am_atom
        !
        implicit none
        !
        class(am_class_tight_binding), intent(in) :: tb
        type(am_class_atom), intent(in) :: atom_m
        type(am_class_atom), intent(in) :: atom_n
        integer , intent(in) :: ip_id
        real(dp), allocatable :: V(:,:)
        real(dp), allocatable :: S(:,:)
        integer :: id
        integer :: m,n
        !
        ! note the absolute value
        id = abs(ip_id)
        ! get primitive atom indices
        m = tb%tbvsk(id)%dims(1)
        n = tb%tbvsk(id)%dims(2)
        ! if irreducible pair id is negative, it means the pair was flipped
        if (ip_id.lt.0) then
            allocate(V(n,m))
            allocate(S(n,m))
            V = adjoint( reshape(tb%tbvsk(id)%V, [m,n]) )
            ! atom_m and atom_n are taken from pp on input, as a result they are already in the correct order (no need to flip them)
            ! m and n are taken from ip on input, thus, they need to be flipped.
            S = transp_parity_sign(atom_m=atom_m, atom_n=atom_n)
            ! if (.not.isequal(shape(S),shape(V))) stop 'S * V : dimension mismatch'
            ! Adding the sign flip here seems to cause H to not be Hermitian...
            ! Looking at Chadi/Cohen and Vogl. There doesn't seem to be a sign flip.
            ! The sign flip is probably only used in determining the irreducible matrix elements.
            V = V ! * S
        else
            if (.not.isequal(tb%tbvsk(id)%dims,[atom_m%norbitals,atom_n%norbitals])) stop 'dims /= norbitals'
            allocate(V(m,n))
            V = reshape(tb%tbvsk(id)%V, [m,n])
        endif
        !
    end function   get_Vsk

    subroutine     set_Vsk_dummies(tb,flags)
        ! flags = seq/rand
        implicit none
        !
        class(am_class_tight_binding), intent(inout) :: tb
        character(*), intent(in) :: flags
        integer :: i, j, k, alpha, beta
        !
        ! set irreducible matrix elements
        do i = 1, tb%nVs
            if      (index(flags, 'seq').ne.0) then
                tb%V(i) = i
            elseif  (index(flags,'rand').ne.0) then
                stop 'rand does not seem to be working here'
                tb%V(i) = rand()
            else
                stop 'Unknown flag set_Vsk_dummies'
            endif
        enddo
        !
        ! clear irreducible matrix elements in each shell
        do k = 1, tb%nshells
            tb%tbvsk(k)%V = 0
        enddo
        !
        ! transfer irreducible matrix elements to each shell
        do i = 1, tb%nVs
            k     = tb%V_ind(1,i)
            alpha = tb%V_ind(2,i)
            beta  = tb%V_ind(3,i)
            j     = sub2ind(dims=tb%tbvsk(k)%dims, sub=[alpha,beta])
            !
            tb%tbvsk(k)%V(j) = tb%V(i)
        enddo
        !
        ! once the irreducible matrix elements have been copied, symmetrize
        do k = 1, tb%nshells
            tb%tbvsk(k)%V(:) = matmul(tb%tbvsk(k)%relations,tb%tbvsk(k)%V(:))
        enddo
        ! things are looking good up to this point.
        !  ... tb%tbvsk(k)%V(:) =
        !              1.00000        0.00000        0.00000        0.00000
        !              0.00000        2.00000        0.00000        0.00000
        !              0.00000        0.00000        2.00000        0.00000
        !              0.00000        0.00000        0.00000        2.00000
        !  ... tb%tbvsk(k)%V(:) =
        !              3.00000        4.00000        4.00000        4.00000
        !             -4.00000        6.00000        5.00000        5.00000
        !             -4.00000        5.00000        6.00000        5.00000
        !             -4.00000        5.00000        5.00000        6.00000
    end subroutine set_Vsk_dummies

    subroutine     write_irreducible_Vsk(tb)
        !
        implicit none
        !
        class(am_class_tight_binding), intent(in) :: tb
        integer :: k ! shell index
        integer :: alpha
        integer :: beta
        integer :: i
        !
        do i = 1, tb%nVs
            ! get shell index
            k     = tb%V_ind(1,i)
            ! get orbital indices
            alpha = tb%V_ind(2,i)
            beta  = tb%V_ind(3,i)
            ! write 
            write(*,*) k, alpha, beta, tb%V(i)
        enddo
    end subroutine write_irreducible_Vsk

    subroutine     write_tb_explicit(tb,tbpg,ic,pp)
        !
        implicit none
        !
        class(am_class_tight_binding), intent(in) :: tb    ! irreducible tight binding matrix elements
        type(am_class_tb_group)           , intent(in) :: tbpg  ! point group in tight binding representation
        type(am_class_irre_cell)          , intent(in) :: ic    ! irreducible cell
        type(am_class_prim_pair)          , intent(in) :: pp    ! primitive pairs
        real(dp), allocatable, target :: Hsub_target(:,:)
        real(dp), allocatable, target :: pg_target(:,:,:)
        integer , allocatable :: S(:) ,E(:)
        real(dp), pointer :: Dm(:,:), Dn(:,:)
        real(dp), pointer :: Hsub(:,:)
        integer :: fid
        integer :: k,p, i,j, m,n
        !
        ! allocate space for vectors demarking start and end of Hamiltonian subsection
        allocate(Hsub_target(tbpg%nbases,tbpg%nbases))
        allocate(pg_target, source=tbpg%sym)
        allocate(S, source=tbpg%H_start)
        allocate(E, source=tbpg%H_end)
        !
        call tostring_set(sep=' ')
        !
        fid = 1
        open(unit=fid,file='outfile.tb_matrix_elements',status='replace',action='write')
            ! print stuff
            write(fid,'(a,a)')                 tostring(pp%nshells), ' primitive shells'
            write(fid,'(a,a)')                          tostring(S), ' subregion start'
            write(fid,'(a,a)')                          tostring(E), ' subregion end'
            do k = 1, pp%nshells
            ! abbreviations
            m = pp%shell(k)%m
            n = pp%shell(k)%n
            i = pp%shell(k)%i
            j = pp%shell(k)%j
            ! print stuff
            write(fid,'(a,a)')                       repeat('=',79), repeat('=',79)
            write(fid,'(a,a)')              centertitle(tostring(k)//' shell',158)
            write(fid,'(a,a)')                          tostring(m), ' primitive m'
            write(fid,'(a,a)')                          tostring(n), ' primitive n'
            write(fid,'(a,a)')         tostring(pp%shell(k)%natoms), ' atoms in shell'
            do p = 1, pp%shell(k)%natoms
            ! abbreviations
            Dm   => pg_target(S(m):E(m), S(m):E(m), pp%shell(k)%pg_id(p) )
            Dn   => pg_target(S(n):E(n), S(n):E(n), pp%shell(k)%pg_id(p) )
            Hsub => Hsub_target(S(m):E(m), S(n):E(n))
            ! print stuff
            Hsub = get_Vsk(tb=tb, ip_id=pp%ip_id(k), atom_m=ic%atom(i), atom_n=ic%atom(j))
            call disp(unit=fid, fmt='f18.10',X=matmul(matmul(transpose(Dm), Hsub), Dn),title=tostring(pp%shell(k)%tau_cart(1:3,p),fmt='f10.5')//' [cart.]',style='above')
            enddo
            enddo
        close(fid)
        !
    end subroutine write_tb_explicit

    subroutine     export_tb2matlab(tb,tbpg,ip,pp)
        ! 
        ! Get tight binding Hamiltonian at kpt.
        ! The only reason ip is included is for the mask.
        ! 
        implicit none
        !
        class(am_class_tight_binding), intent(in) :: tb     ! tight binding matrix elements
        type(am_class_tb_group)      , intent(in) :: tbpg   ! point group in tight binding representation
        type(am_class_prim_pair)     , intent(in) :: pp     ! primitive pairs
        type(am_class_irre_pair)     , intent(in) :: ip     ! irreducible pairs
        integer, allocatable :: Sv(:) ! start and end if irreducible vector
        integer, allocatable :: Ev(:) ! start and end if irreducible vector
        integer, allocatable :: S(:)  ! start and end of hamiltonian subsection
        integer, allocatable :: E(:)  ! start and end of hamiltonian subsection
        integer :: i,j,l,k,p,m,n, fid
        ! parameters
        character(100) :: H_fnc_name
        character(100) :: V_fnc_name
        character(100) :: str_H
        character(100) :: str_Dm
        character(100) :: str_Dn
        character(100) :: str_tau
        character(100) :: str_Vab
        !
        H_fnc_name = 'getH'
        V_fnc_name = 'getV'
        !
        ! create abbreviations
        allocate(S, source=tbpg%H_start)
        allocate(E, source=tbpg%H_end)
        !
        ! export symmetry relations for each irreducible pair
        !
        allocate(Sv(ip%nshells))
        allocate(Ev(ip%nshells))
        j = 1
        do i = 1, ip%nshells
            Sv(i) = j
            call export_relations2matlab(relations=tb%tbvsk(i)%relations, dims=tb%tbvsk(i)%dims, fnc_name=trim(V_fnc_name)//tostring(i))
            j = j + count(get_independent(tb%tbvsk(i)%relations)) - 1
            Ev(i) = j
        enddo
        !
        fid = 1
        open(unit=fid,file=trim(H_fnc_name)//'.m',status='replace',action='write')
            !
            write(fid,'(a,a,a)') 'function [H] = ', trim(H_fnc_name), '(pg,v,kpt)'
            write(fid,'(a)') 'i2pi = 2*sqrt(-1)*pi;'
            write(fid,'(a)') 'H(1:'//tostring(E(ip%nshells))//',1:'//tostring(E(ip%nshells))//') = 0;'
                !
                ! construct Hamiltonian
                do l = 1, ip%nshells
                do k = 1, pp%nshells
                if (abs(pp%ip_id(k)).eq.l) then
                    !
                    write(fid,'(a)') '% irreducible shell '//tostring(l)
                    ! compute bloch sum by loop over atoms in shell
                    do p = 1, pp%shell(k)%natoms
                        ! primitive atom indicies
                        m = pp%shell(k)%m
                        n = pp%shell(k)%n
                        str_H  =      'H('//tostring(S(m))//':'//tostring(E(m))//','//tostring(S(n))//':'//tostring(E(n))//')'
                        str_Dm = 'pg.sym('//tostring(S(m))//':'//tostring(E(m))//','//tostring(S(m))//':'//tostring(E(m))//','//tostring(pp%shell(k)%pg_id(p))//')'
                        str_Dn = 'pg.sym('//tostring(S(n))//':'//tostring(E(n))//','//tostring(S(n))//':'//tostring(E(n))//','//tostring(pp%shell(k)%pg_id(p))//')'
                        str_tau= '['//tostring(pp%shell(k)%tau_cart(1:3,p),fmt='SP,f10.5')//']'
                        str_Vab= trim(V_fnc_name)//tostring(abs(pp%ip_id(k)))//'(v('//tostring(Sv(abs(pp%ip_id(k))))//':'//tostring(Ev(abs(pp%ip_id(k))))//'))'
                        if (pp%ip_id(k).lt.0) str_Vab=trim(str_Vab)//"'" ! transpose
                        !
                        write(fid,'(a)',advance='no') trim(str_H)//' = '//trim(str_H)//' + '
                        write(fid,'(a)',advance='no') trim(str_Dm)//"' * "//trim(str_Vab)//" * "//trim(str_Dn)
                        write(fid,'(a)',advance='no') ' * exp(i2pi*dot(kpt,'//trim(str_tau)//'));'
                        write(fid,*)
                        !
                    enddo
                endif
                enddo
                enddo
                !
            write(fid,'(a)') 'end'
        close(fid)
    end subroutine export_tb2matlab








































end module am_tight_binding











