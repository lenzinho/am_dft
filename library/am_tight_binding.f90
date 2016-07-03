module am_tight_binding

    use dispmodule
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
    use am_dispersion
    use am_brillouin_zone

    implicit none

    private

    type, public :: am_class_tight_binding
        integer :: nshells                     ! how many shells irreducible atoms
        type(am_class_tensor), allocatable :: tens(:)  ! tens(nshells)
        integer , allocatable :: nVs            ! number of irreducible matrix element values
        real(dp), allocatable :: V(:)           ! their values
        integer , allocatable :: V_ind(:,:)     ! their indices V_ind( [i,j], nVs); j = compound [alpha x beta index, see initialize_tb subroutine]
        contains
        procedure :: set_Vsk
        procedure :: initialize_tb
        procedure :: get_hamiltonian
        procedure :: test_hamiltonian
        procedure :: export_to_matlab
        procedure :: read_matrix_elements_irr
        procedure :: write_matrix_elements_irr
        procedure :: write_matrix_elements_full
        procedure :: optimize_matrix_elements
    end type am_class_tight_binding

    type, private :: am_class_tb_optimizer
        integer :: maxiter      ! maximum number of iterations  
        integer :: nbands
        integer :: nkpts
        integer :: nxs
        integer :: nrs
        integer               :: skip_band
        integer , allocatable :: selector_shell(:)
        integer , allocatable :: selector_kpoint(:)
        real(dp), allocatable :: x(:) 
        real(dp), allocatable :: r(:)
        real(dp)              :: rms
    end type am_class_tb_optimizer

    type, public, extends(am_class_dispersion) :: am_class_dispersion_tb
        complex(dp), allocatable :: C(:,:,:) ! tight binding coefficients
        contains
        procedure :: get_dispersion
    end type am_class_dispersion_tb

contains

    subroutine     initialize_tb(tb,pg,ip,ic,pc,opts)
        !
        implicit none
        !
        class(am_class_tight_binding), intent(out) :: tb   ! tight binding parameters
        type(am_class_point_group)   , intent(in)  :: pg   ! seitz point group
        type(am_class_prim_cell)     , intent(in)  :: pc   ! primitive cell
        type(am_class_irre_cell)     , intent(in)  :: ic   ! irreducible cell
        type(am_class_irre_pair)     , intent(in)  :: ip   ! irreducible pairs
        type(am_class_options)       , intent(in)  :: opts
        logical, allocatable :: is_independent(:)
        integer :: sub(2)
        integer :: i,j,k
        !
        !
        if (opts%verbosity.ge.1) call print_title('Symmetry-adapted tight-binding parameters')
        !
        tb%nshells = ip%nshells
        allocate(tb%tens(tb%nshells))
        ! get rank, dims, flags, property, and symmetry relations
        do k = 1, tb%nshells
            call tb%tens(k)%symmetrize(pg=pg,pc=pc,ic=ic,shell=ip%shell(k),opts=opts,property='irreducible tight binding shell '//tostring(k))
        enddo
        ! get number of independent (irreducible) matrix elements
        tb%nVs = 0
        do i = 1, tb%nshells
            tb%nVs = tb%nVs + count(get_independent(tb%tens(i)%relations)) 
        enddo
        ! allocate space for independent (irreducible) matrix elements V
        allocate(tb%V(tb%nVs))
        tb%V = 0
        ! indices : V_ind( [i, subd2ind(dims=tb%tens(i)%dims,sub=[alpha,beta])], nVs)
        allocate(tb%V_ind(3,tb%nVs))
        tb%V_ind = 0
        k=0
        do i = 1, tb%nshells
            is_independent = get_independent(tb%tens(i)%relations)
            do j = 1, product(tb%tens(i)%dims)
            if (is_independent(j)) then
                k=k+1
                sub = ind2sub(dims=tb%tens(i)%dims, ind=j)
                tb%V_ind(1,k) = i      ! shell
                tb%V_ind(2,k) = sub(1) ! alpha
                tb%V_ind(3,k) = sub(2) ! beta
            endif
            enddo
        enddo
        ! write template for irreducible matrix elements
        if (fexists('infile.tb_matrix_elements_irreducible').ne.0) then 
            ! print title
            call print_title('Input tight-binding matrix elements')
            ! read matrix elements
            call tb%read_matrix_elements_irr()
            ! display input file on stdout
            call disp(x=[0],zeroas=' ',advance='no')
            call disp(x=tb%V,title='V',style='underline',advance='no')
            call disp(x=tb%V_ind(1,:),title='shell',style='underline',advance='no')
            call disp(x=tb%V_ind(2,:),title='alpha',style='underline',advance='no')
            call disp(x=tb%V_ind(3,:),title='beta',style='underline',advance='yes')
        else
            ! print title
            write(*,'(a,a)') flare, 'infile.tb_matrix_elements_irreducible not found. Using template instead.'
            !
            call tb%set_Vsk(flags='seq')
            !
            call tb%write_matrix_elements_irr()
        endif
        !
        contains
        subroutine     initialize_tens(tens,pc,ic,shell)
            !
            implicit none
            !
            type(am_class_tensor) , intent(out) :: tens
            type(am_class_prim_cell), intent(in)  :: pc
            type(am_class_irre_cell), intent(in)  :: ic
            type(am_shell_cell)     , intent(in)  :: shell
            integer :: i
            !
            tens%property = 'tight binding'
            tens%rank  = 2
            ! detetmine if irreducible atoms are the same
            if (shell%i.eq.shell%j) then
                tens%flags = 'i==j'
            else
                tens%flags = 'i/=j'
            endif
            ! set dimensions of matrix elements
            allocate(tens%dims(tens%rank))
            tens%dims = 0 
            ! dimension of Hamiltonian subsection corresponding to primitive atom m (irreducible atoms i)
            do i = 1, ic%atom(shell%i)%nazimuthals
                tens%dims(1) = tens%dims(1) + ic%atom(shell%i)%azimuthal(i)*2+1
            enddo
            ! dimension of Hamiltonian subsection corresponding to primitive atom n (irreducible atoms j)
            do i = 1, ic%atom(shell%j)%nazimuthals
                tens%dims(2) = tens%dims(2) + ic%atom(shell%j)%azimuthal(i)*2+1
            enddo
            ! allocate space for matrix elements
            allocate(tens%V(tens%dims(1)*tens%dims(2)))
            !
        end subroutine initialize_tens
        subroutine     get_shell_relations(tens,pg,ic,shell,opts)
                !
                implicit none
                !
                type(am_class_tensor)   , intent(inout) :: tens
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
                call flat_ig%get_flat_intrinsic_group(tens=tens, atom_m=ic%atom(shell%j), atom_n=ic%atom(shell%i) )
                ! determine stabilizers relations
                call stab%get_stabilizer_group(pg=pg, v=shell%tau_cart(1:3,1), opts=opts, flags='cart')
                ! get stabilizer symmetries in the flattened hamiltonin basis
                call flat_pg%get_flat_point_group(tens=tens, pg=stab, atom_m=ic%atom(shell%j), atom_n=ic%atom(shell%i))
                ! get combined relations
                tens%relations = combine_relations(relationsA=flat_pg%relations, relationsB=flat_ig%relations)
                ! correct rounding error
                call correct_rounding_error(tens%relations)
                ! 
        end subroutine get_shell_relations
    end subroutine initialize_tb

    function       get_hamiltonian(tb,tbpg,pp,kpt,selector_shell) result(H)
        ! 
        ! Get tight binding Hamiltonian at kpt.
        ! 
        implicit none
        !
        class(am_class_tight_binding), intent(in) :: tb     ! tight binding matrix elements
        type(am_class_tb_group)      , intent(in) :: tbpg   ! point group in tight binding representation
        type(am_class_prim_pair)     , intent(in) :: pp     ! primitive pairs
        real(dp)                     , intent(in) :: kpt(3) ! cart
        integer, optional            , intent(in) :: selector_shell(:)
        integer    , allocatable :: selector_shell_internal(:)
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
        integer :: ip_nshells
        ! get number of irreducile shell
        ip_nshells = maxval(abs(pp%ip_id(:)))
        ! set selector
        if (present(selector_shell)) then
            allocate(selector_shell_internal,source=selector_shell)
        else
            allocate(selector_shell_internal,source=[1:ip_nshells])
        endif
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
        do l = 1, ip_nshells
        if ( any(selector_shell_internal.eq.l) ) then
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
                    Dm   =>   pg_target(S(m):E(m), S(m):E(m), pp%shell(k)%pg_id(p))
                    Dn   =>   pg_target(S(n):E(n), S(n):E(n), pp%shell(k)%pg_id(p))
                    ! get matrix elements (initialize Hsub)
                    Hsub = get_Vsk(tb=tb, ip_id=pp%ip_id(k))
                    ! rotate matrix elements as needed to get from the tau_frac(:,1) => tau_frac(:,x)
                    Hsub = matmul(matmul(transpose(Dm), Hsub), Dn)
                    ! multiply exponential factor from Bloch sum
                    Hsub = Hsub * exp(-itwopi*dot_product(pp%shell(k)%tau_cart(1:3,p), kpt)) ! kpt [rec. cart.]
                    ! this pair's contribution to the Hamiltonian
                    H(S(m):E(m), S(n):E(n)) = H(S(m):E(m), S(n):E(n)) + Hsub
                enddo
            endif
            enddo
            !
            if (debug) then
            if (.not.ishermitian(H)) stop 'ERROR [get_hamiltonian]: H is not Hermitian!'
            endif
        endif
        enddo
    end function   get_hamiltonian

    subroutine     test_hamiltonian(tb,tbpg,pp)
        ! makes sure Hamiltonian at Gamma commutes with all point symmetry operations
        implicit none
        !
        class(am_class_tight_binding), intent(inout) :: tb   ! tight binding matrix elements
        type(am_class_tb_group)      , intent(in)    :: tbpg ! point group in tight binding representation
        type(am_class_prim_pair)     , intent(in)    :: pp   ! primitive pairs
        complex(dp), allocatable :: H(:,:)
        real(dp)   , allocatable :: R(:,:)
        integer :: i, j, k
        integer :: ip_nshells
        real(dp):: kpt(3,1)
        !
        call print_title('Checking model Hamiltonian at Gamma')
        !
        ! get number of irreducile shell
        ip_nshells = maxval(abs(pp%ip_id(:)))
        ! allocate space for symmetry in tb basis
        allocate(R(tbpg%nbases,tbpg%nbases))
        ! initilize tb matrix elements
        call tb%set_Vsk(flags='zero')
        ! loop over irreducible shells
        do j = 1, ip_nshells
            !
            H = tb%get_hamiltonian(tbpg=tbpg, pp=pp, kpt=real([0,0,0],dp),selector_shell=[j])
            !
            ! check that H commutes with all point symmetries
            do i = 1, tbpg%nsyms
                R = tbpg%sym(:,:,i)
                if (.not.isequal(matmul(H,R),matmul(R,H))) then
                    call disp('H',H,style='above')
                    call disp('R',R,style='above')
                    call disp('[H,R]',matmul(H,R)-matmul(R,H),style='above')
                    call disp("R'*H*R - H",matmul(matmul(transpose(R),H),R)-H,style='above')
                    stop 'ERROR [test_hamiltonian]: H does not commute with point symmetries.'
                endif
            enddo
            !
            write(*,'(a,a)') flare ,'shell '//tostring(j)//': OK!'
        enddo
        !
        call print_title('Checking model Hamiltonian at arbitrary k-points')
        !
        do k = 1, 10
            kpt = reshape([rand(),rand(),rand()],[3,1])
            ! kpt = matmul(uc%recbas,kpt), ideally... it should be in cart.
            write(*,'(a,a)') flare, 'k-point '//tostring(k)//' = '//tostring(pack(kpt,.true.),fmt='f18.10')
            do j = 1, ip_nshells
                write(*,'(a,a)',advance='no') flare, 'shell '//tostring(j)//' '
                ! ignore all irreducible shells, except j
                H = tb%get_hamiltonian(tbpg=tbpg, pp=pp, kpt=real([0,0,0],dp),selector_shell=[j])
                !
                if (.not.ishermitian(H)) then
                    stop 'ERROR [test_hamiltonian]: H is not Hermitian'
                endif
                write(*,*) 'OK!'
            enddo
        enddo
        ! reset tb matrix elments to zero
        call tb%set_Vsk(flags='zero')
    end subroutine test_hamiltonian

    function       get_Vsk(tb,ip_id) result(V)
        ! irreducible pair (ip_id), can be negative (corresponds to pair n-m rathet than m-n, on the
        ! opposite [upper/lower] side of the Hamiltonian), in which case adjoint of V is returned
        !
        implicit none
        !
        class(am_class_tight_binding), intent(in) :: tb
        integer , intent(in) :: ip_id
        real(dp), allocatable :: V(:,:)
        integer :: id
        integer :: m,n
        !
        ! note the absolute value
        id = abs(ip_id)
        ! get primitive atom indices
        m = tb%tens(id)%dims(1)
        n = tb%tens(id)%dims(2)
        ! if irreducible pair id is negative, it means the pair was flipped
        if (ip_id.lt.0) then
            allocate(V, source=adjoint(reshape(tb%tens(id)%V,[m,n])) )
        else
            allocate(V, source=        reshape(tb%tens(id)%V,[m,n])  )
        endif
        !
    end function   get_Vsk

    subroutine     set_Vsk(tb,V,flags)
        ! flags = seq/rand
        implicit none
        !
        class(am_class_tight_binding), intent(inout) :: tb
        real(dp)    , optional       , intent(in)    :: V(:)
        character(*), optional       , intent(in)    :: flags
        integer :: i, j, k, alpha, beta
        !
        ! set irreducible matrix elements
        if      (index(flags,   'seq').ne.0) then
            do i = 1, tb%nVs
                tb%V(i) = i
            enddo
        elseif  (index(flags,  'rand').ne.0) then
            stop 'rand does not seem to be working here'
            do i = 1, tb%nVs
                tb%V(i) = rand()
            enddo
        elseif  (index(flags,  'zero').ne.0) then
            do i = 1, tb%nVs
                tb%V(i) = 0
            enddo
        elseif  (present(V)                ) then
            if (size(V).ne.tb%nVs) stop 'ERROR: V /= nVs dimension mismatch'
            do i = 1, tb%nVs
                tb%V(i) = V(i)
            enddo
        else
            stop 'Unknown flag set_Vsk_dummies'
        endif
        ! clear irreducible matrix elements in each shell
        do k = 1, tb%nshells
            tb%tens(k)%V = 0
        enddo
        ! transfer irreducible matrix elements to each shell
        do i = 1, tb%nVs
            k     = tb%V_ind(1,i)
            alpha = tb%V_ind(2,i)
            beta  = tb%V_ind(3,i)
            j     = sub2ind(dims=tb%tens(k)%dims, sub=[alpha,beta])
            !
            tb%tens(k)%V(j) = tb%V(i)
        enddo
        ! once the irreducible matrix elements have been copied, symmetrize
        do k = 1, tb%nshells
            tb%tens(k)%V(:) = matmul(tb%tens(k)%relations,tb%tens(k)%V(:))
        enddo
        ! things are looking good up to this point.
        !  ... tb%tens(k)%V(:) =
        !              1.00000        0.00000        0.00000        0.00000
        !              0.00000        2.00000        0.00000        0.00000
        !              0.00000        0.00000        2.00000        0.00000
        !              0.00000        0.00000        0.00000        2.00000
        !  ... tb%tens(k)%V(:) =
        !              3.00000        4.00000        4.00000        4.00000
        !             -4.00000        6.00000        5.00000        5.00000
        !             -4.00000        5.00000        6.00000        5.00000
        !             -4.00000        5.00000        5.00000        6.00000
    end subroutine set_Vsk

    subroutine     write_matrix_elements_irr(tb)
        !
        implicit none
        !
        class(am_class_tight_binding), intent(in) :: tb
        integer :: k ! shell index
        integer :: alpha
        integer :: beta
        integer :: i
        integer :: fid
        ! create tb dir
        call execute_command_line ('mkdir -p '//trim(outfile_dir_tb))
        ! write point group
        fid = 1
        open(unit=fid,file=trim(outfile_dir_tb)//'/'//'outfile.tb_matrix_elements_irreducible',status='replace',action='write')
            !
            write(fid,'(a,a16,3a6)') '#', 'V', 'shell', 'alpha', 'beta'
            do i = 1, tb%nVs
                ! get shell index
                k     = tb%V_ind(1,i)
                ! get orbital indices
                alpha = tb%V_ind(2,i)
                beta  = tb%V_ind(3,i)
                ! write 
                write(fid,'(i5,f16.8,3i6)') i, tb%V(i), k, alpha, beta
            enddo
            !
        close(fid)
    end subroutine write_matrix_elements_irr

    subroutine     read_matrix_elements_irr(tb)
        !
        implicit none
        !
        class(am_class_tight_binding), intent(inout) :: tb
        character(maximum_buffer_size) :: buffer ! read buffer
        character(len=:), allocatable :: word(:) ! read buffer
        integer :: fid
        real(dp), allocatable :: V(:)
        integer :: k
        integer :: alpha
        integer :: beta
        integer :: i
        integer :: j
        !
        allocate(V(tb%nVs))
        !
        fid = 1
        open(unit=fid,file='infile.tb_matrix_elements_irreducible',status="old",action='read')
            ! skip header
            read(fid,*)
            ! read matrix elements
            do j = 1, tb%nVs
                read(unit=fid,fmt='(a)') buffer
                word = strsplit(buffer,delimiter=' ')
                !
                read(word(1),*) i
                read(word(2),*) V(i)
                read(word(3),*) k
                read(word(4),*) alpha
                read(word(5),*) beta
                !
                if (.not.isequal(    k,tb%V_ind(1,i))) stop 'ERROR [read_matrix_elements_irr]: k /= V_ind(1,i)'
                if (.not.isequal(alpha,tb%V_ind(2,i))) stop 'ERROR [read_matrix_elements_irr]: alpha /= V_ind(2,i)'
                if (.not.isequal( beta,tb%V_ind(3,i))) stop 'ERROR [read_matrix_elements_irr]: beta /= V_ind(3,i)'
            enddo
            call tb%set_Vsk(V)
        close(fid)
    end subroutine read_matrix_elements_irr

    subroutine     write_matrix_elements_full(tb,tbpg,pp)
        !
        implicit none
        !
        class(am_class_tight_binding)     , intent(in) :: tb    ! irreducible tight binding matrix elements
        type(am_class_tb_group)           , intent(in) :: tbpg  ! point group in tight binding representation
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
            ! create directory
            call execute_command_line ('mkdir -p '//trim(outfile_dir_tb))
            ! export symmetry
            fid = 1
            open(unit=fid,file=trim(outfile_dir_tb)//'/'//'outfile.tb_matrix_elements',status='replace',action='write')
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
                Hsub = get_Vsk(tb=tb, ip_id=pp%ip_id(k)) 
                ! print
                call disp(unit=fid, fmt='f18.10',x=matmul(matmul(transpose(Dm), Hsub), Dn),title=tostring(pp%shell(k)%tau_cart(1:3,p),fmt='f10.5')//' [cart.]',style='above')
                enddo
                enddo
            close(fid)
        ! return delimiter to ","
        call tostring_set(sep=',')
        !
    end subroutine write_matrix_elements_full

    subroutine     export_to_matlab(tb,ip,tbpg,pp,flags)
        ! flags = symbolic/numerical cart/frac
        implicit none
        !
        class(am_class_tight_binding), intent(in) :: tb     ! tight binding matrix elements
        type(am_class_tb_group)      , intent(in) :: tbpg   ! point group in tight binding representation
        type(am_class_prim_pair)     , intent(in) :: pp     ! primitive pairs
        type(am_class_irre_pair)     , intent(in) :: ip     ! irreducible pairs
        character(*)                 , intent(in) :: flags  ! symbolic/numerical cart/frac
        integer, allocatable :: Sv(:) ! start and end if irreducible vector
        integer, allocatable :: Ev(:) ! start and end if irreducible vector
        integer, allocatable :: S(:)  ! start and end of hamiltonian subsection
        integer, allocatable :: E(:)  ! start and end of hamiltonian subsection
        integer :: i,j,l,k,p,m,n,end,fid
        ! parameters
        character(100) :: H_fnc_name
        character(100) :: V_fnc_name
        character(100) :: str_H
        character(100) :: str_Dm
        character(100) :: str_Dn
        character(100) :: str_tau
        character(100) :: str_Vab
        ! tau substitution
        ! integer :: ntaus
        ! real(dp), allocatable :: taus(:)
        !
        ! set function names
        if     (index(flags,'symbolic').ne.0) then
            H_fnc_name = 'get_H_symbolic'
        elseif (index(flags, 'numeric').ne.0) then
            H_fnc_name = 'get_H_numeric'
        else
            stop 'ERROR [export_to_matlab]: flags != symbolic/numeric'
        endif
        !
        if     (index(flags,'cart').ne.0) then
            H_fnc_name=trim(H_fnc_name)//'_cart'
        elseif (index(flags,'frac').ne.0) then
            H_fnc_name=trim(H_fnc_name)//'_frac'
        else
            stop 'ERROR [export_to_matlab]: flags != cart/frac'
        endif
        !
        ! get number of unique coordinates
        ! if     (index(flags, 'symbolic').ne.0) then
        ! allocate(taus( pp%nshells * maxval(pp%shell(:)%natoms)))
        ! ntaus = 0
        ! do i = 1, pp%nshells
        ! do j = 1, pp%shell(k)%natoms
        ! do k = 1, 3
        !     if (abs(pp%shell(i)%tau_frac(k,j))).gt.opts%prec) then
        !     if (.not.issubset(taus(1:k),pp%shell(i)%tau_frac(k,j))) then
        !         ntaus = ntaus+1
        !         tau(ntaus) = pp%shell(i)%tau_frac(k,j)
        !     endif
        !     endif
        ! enddo
        ! enddo
        ! enddo
        ! endif
        !
        V_fnc_name = 'getV'
        ! create abbreviations
        allocate(S, source=tbpg%H_start)
        allocate(E, source=tbpg%H_end)
        end = size(tbpg%H_end)
        ! allocate start and end vectors
        allocate(Sv(ip%nshells))
        allocate(Ev(ip%nshells))
        j = 0
        do i = 1, ip%nshells
            j = j + 1
            Sv(i) = j
            j = j + count(get_independent(tb%tens(i)%relations)) - 1
            Ev(i) = j
        enddo
        ! export hamiltonian
        fid = 1
        ! create directory
        call execute_command_line ('mkdir -p '//trim(outfile_dir_tb))
        call execute_command_line( 'mkdir -p '//trim(outfile_dir_tb)//'/matlab')
        ! export file
        open(unit=fid,file=trim(outfile_dir_tb)//'/matlab/'//trim(H_fnc_name)//'.m',status='replace',action='write')
            write(fid,'(a,a,a)') 'function [H] = ', trim(H_fnc_name), '(pg,v,kpt,shell)'
            write(fid,'(a)') 'i2pi = 2*sqrt(-1)*pi;'
            write(fid,'(a)') 'H(1:'//tostring(E(end))//',1:'//tostring(E(end))//') = 0;'
            if (index(flags,'symb').ne.0) write(fid,'(a)') 'H = sym(H);'
                ! construct Hamiltonian
                do l = 1, ip%nshells
                write(fid,'(a)') 'if any(shell=='//tostring(l)//')'
                    do k = 1, pp%nshells
                    if (abs(pp%ip_id(k)).eq.l) then
                        ! compute bloch sum by loop over atoms in shell
                        do p = 1, pp%shell(k)%natoms
                            ! determine how to write stuff
                            m = pp%shell(k)%m
                            n = pp%shell(k)%n
                            str_H  =      'H('//tostring(S(m))//':'//tostring(E(m))//','//tostring(S(n))//':'//tostring(E(n))//')'
                            str_Dm = 'pg.sym('//tostring(S(m))//':'//tostring(E(m))//','//tostring(S(m))//':'//tostring(E(m))//','//tostring(pp%shell(k)%pg_id(p))//')'
                            str_Dn = 'pg.sym('//tostring(S(n))//':'//tostring(E(n))//','//tostring(S(n))//':'//tostring(E(n))//','//tostring(pp%shell(k)%pg_id(p))//')'
                            if     (index(flags,'cart').ne.0) then
                                str_tau= '['//tostring(pp%shell(k)%tau_cart(1:3,p),fmt='SP,f10.5')//']'
                            elseif (index(flags,'frac').ne.0) then
                                str_tau= '['//tostring(pp%shell(k)%tau_frac(1:3,p),fmt='SP,f10.5')//']'
                            endif
                            str_Vab= trim(V_fnc_name)//tostring(l)//'(v('//tostring(Sv(l))//':'//tostring(Ev(l))//'))'
                            if (pp%ip_id(k).lt.0) str_Vab=trim(str_Vab)//"'" ! adjoint
                            ! write stuff
                            write(fid,'(a)',advance='no') trim(str_H)       // ' = ' // trim(str_H)   // ' + '
                            write(fid,'(a)',advance='no') trim(str_Dm)//"'" // ' * ' // trim(str_Vab) // ' * ' // trim(str_Dn)
                            write(fid,'(a)',advance='no')     ' * exp(i2pi*dot(kpt,' // trim(str_tau) // '));'
                            write(fid,*)
                        enddo
                    endif
                    enddo
                write(fid,'(a)') 'end'
                enddo
            write(fid,'(a)') 'end'
            ! export symmetry relations for each irreducible pair (append to the same file)
            do i = 1, ip%nshells
                call export_relations2matlab(relations=tb%tens(i)%relations, dims=tb%tens(i)%dims, fnc_name=trim(V_fnc_name)//tostring(i), fid_append=fid, flags='append')
            enddo
        close(fid)
    end subroutine export_to_matlab

    subroutine     optimize_matrix_elements(tb,tbpg,bz,dr,ip,pp,opts)
        !
        implicit none
        !   
        class(am_class_tight_binding) , intent(inout) :: tb
        type(am_class_tb_group)       , intent(in) :: tbpg
        type(am_class_bz)             , intent(in) :: bz
        type(am_class_dispersion_vasp), intent(in) :: dr
        type(am_class_irre_pair)      , intent(in) :: ip
        type(am_class_prim_pair)      , intent(in) :: pp
        type(am_class_options)        , intent(in) :: opts
        type(am_class_tb_optimizer) :: ft
        !
        if (opts%verbosity.ge.1) call print_title('Optimizing matrix elements')
        !
        ! initialize optimizer
        ft%maxiter    = 10
        ft%nbands     = maxval(tbpg%H_end)
        ft%nkpts      = bz%nkpts
        ft%nxs        = tb%nVs
        ft%skip_band  = opts%skip_band
        ft%nrs        = bz%nkpts * ft%nbands
        allocate(ft%selector_shell,  source=[1:ip%nshells])
        allocate(ft%selector_kpoint, source=[1:bz%nkpts])
        allocate(ft%r(ft%nrs))
        ft%r = 0
        allocate(ft%x(ft%nxs))
        ft%x = 0
        !
        if (debug) then
        write(*,'(a,a,a)') flare, 'irreducible matrix elements = '//tostring(ft%nxs)
        write(*,'(a,a,a)') flare, 'k-points = '                   //tostring(ft%nkpts)
        write(*,'(a,a,a)') flare, 'bands = '                      //tostring(ft%nbands)
        write(*,'(a,a,a)') flare, 'number of bands to skip = '    //tostring(ft%skip_band)
        write(*,'(a,a,a)') flare, 'residual vector length = '     //tostring(ft%nrs)
        write(*,'(a,a,a)') flare, 'max iterations = '             //tostring(ft%maxiter)
        endif
        !
        call perform_optimization(ft=ft, tb=tb, tbpg=tbpg, bz=bz, dr=dr, ip=ip, pp=pp)
        !
        ! call dump(ft%x,'outfile.debug.V')
        !
        contains
        subroutine     perform_optimization(ft,tb,tbpg,bz,dr,pp,ip)
            !
            use mkl_rci
            use mkl_rci_type
            !
            implicit none
            !
            class(am_class_tb_optimizer)  , intent(inout) :: ft
            type(am_class_tight_binding)  , intent(inout) :: tb
            type(am_class_tb_group)       , intent(in) :: tbpg
            type(am_class_bz)             , intent(in) :: bz
            type(am_class_dispersion_vasp), intent(in) :: dr
            type(am_class_irre_pair)      , intent(in) :: ip
            type(am_class_prim_pair)      , intent(in) :: pp
            real(dp), allocatable :: x(:)
            real(dp), allocatable :: FVEC(:)
            real(dp), allocatable :: FJAC(:,:)
            type(handle_tr) :: handle
            real(dp) :: eps(6)
            integer  :: info(6)
            integer  :: rci_request
            integer  :: successful
            integer  :: i
            !
            ! default internal paramters
            eps(1:6) = 0.01
            !
            ! initialize fitting parameters (irreducible matrix elements) 
            allocate(x(ft%nxs))
            x = tb%V 
            ! initialize residual vector
            allocate(FVEC(ft%nrs))
            FVEC = 0.0_dp
            ! initialize jacobian matrix
            allocate(FJAC(ft%nrs,ft%nxs))
            FJAC = 0
            ! handle    in/out: tr solver handle
            ! n         in:     number of function variables
            ! m         in:     dimension of function value
            ! x         in:     solution vector. contains values x for f(x)
            ! eps       in:     precisions for stop-criteria
            ! iter1     in:     maximum number of iterations
            ! iter2     in:     maximum number of iterations of calculation of trial-step
            ! rs        in:     initial step bound
            if (dtrnlsp_init(handle=handle, n=ft%nxs, m=ft%nrs, x=x, eps=eps, iter1=ft%maxiter, iter2=100, rs=100.0D0) /= tr_success) then
                print *, 'error in dtrnlsp_init'
                call mkl_free_buffers
                stop 1
            end if
            ! n         in:     length of x
            ! m         in:     length of f(x)
            ! fjac      in:     [n*m] jacobian matrix
            ! fvec      in:     [m] function values at x ( fvec(i) = y_i - f_i(x) )
            ! eps       in:     precisions for stop-criteria
            ! info      out:    result of input checking
            if (dtrnlsp_check(handle=handle, n=ft%nxs, m=ft%nrs, FJAC=FJAC, FVEC=FVEC, eps=eps, info=info) /= tr_success) then
                print *, ' error in dtrnlspbc_init'
                call mkl_free_buffers
                stop 1
            else
                if ( info(1) /= 0 .or. info(2) /= 0 .or. info(3) /= 0 .or. info(4) /= 0 ) then
                    print *, 'input parameters are not valid'
                    call mkl_free_buffers
                    stop 1
                end if
            end if
            ! set initial rci variables
            rci_request = 0
            successful = 0
            i = 0
            ! write header
            write(*,'(5x,a8,a10,a10,a)') 'iter', 'rms', 'log(d rms)',centertitle('parameters',ft%nxs*10)
            write(*,'(5x,a8,a10,a10,a)') ' '//repeat('-',8-1), ' '//repeat('-',10-1), ' '//repeat('-',10-1), ' '//repeat('-',ft%nxs*10-1)
            ! enter optimization loop
            do while (successful == 0)
                ! handle        in/out: tr solver handle
                ! fvec          in:     vector
                ! fjac          in:     jacobi matrix
                ! rci_request   in/out: return number which denote next step for performing
                if (dtrnlsp_solve(handle=handle, FVEC=FVEC, FJAC=FJAC, rci_request=rci_request) /= tr_success) then
                    print *, 'error in dtrnlsp_solve'
                    call mkl_free_buffers
                    stop 1
                end if
                select case (rci_request)
                case (-1, -2, -3, -4, -5, -6)
                    successful = 1
                case (1)
                    ! update x in tb model
                    call tb%set_Vsk(V=x)
                    ! recalculate function
                    call compute_residual(ft=ft, tb=tb, tbpg=tbpg, bz=bz, dr=dr, pp=pp, R=FVEC)
                    ! save result
                    ft%r   = FVEC 
                    ft%rms = sqrt(sum(FVEC**2)/ft%nrs)
                    ! increase counter
                    i = i + 1
                    ! print results
                    if (i.eq.1) then
                        write(*,'(5x,i8,f10.2,10x,SP,1000f10.2)')   i, ft%rms, x
                    else
                        write(*,'(5x,i8,f10.2,f10.2,SP,1000f10.2)') i, ft%rms, log10(abs(norm(1.0D-14 + x-ft%x ))), x 
                    endif
                    ! save results
                    ft%x = x
                case (2)
                    ! update x in tb model
                    call tb%set_Vsk(V=x)
                    ! compute jacobian matrix (uses central difference)
                    call compute_jacobian(ft=ft, tb=tb, tbpg=tbpg, bz=bz, dr=dr,  pp=pp, FJAC=FJAC)
                end select
            end do
            ! clean up
            if (dtrnlsp_delete(handle) /= tr_success) then
                print *, 'error in dtrnlsp_delete'
                call mkl_free_buffers
                stop 1
            else
                call mkl_free_buffers
            endif
        end subroutine perform_optimization
        subroutine     compute_residual(ft,tb,tbpg,bz,dr,pp,R)
            !
            implicit none
            !
            class(am_class_tb_optimizer)  , intent(in) :: ft
            type(am_class_tight_binding)  , intent(in) :: tb
            type(am_class_tb_group)       , intent(in) :: tbpg
            type(am_class_bz)             , intent(in) :: bz
            type(am_class_dispersion_vasp), intent(in) :: dr 
            type(am_class_prim_pair)      , intent(in) :: pp
            real(dp), allocatable         , intent(out):: R(:)
            type(am_class_dispersion_tb) :: dr_tb 
            integer, allocatable :: inds(:)
            integer :: i
            ! index vector
            allocate(inds(ft%nbands))
            ! create residual vector
            allocate(R(ft%nrs))
            R = 0
            ! get tb dispersion 
            call dr_tb%get_dispersion(tb=tb,tbpg=tbpg,bz=bz,pp=pp)
            ! loop over kpoints
            do i = 1, bz%nkpts
                ! get indices
                inds = [1:ft%nbands] + (i-1)*ft%nbands
                ! calculate residual vector indices
                R(inds) = dr%E([1:ft%nbands]+ft%skip_band, i ) - dr_tb%E(:,i)
            enddo
            !
        end subroutine compute_residual
        subroutine     compute_jacobian(ft,tb,tbpg,bz,dr,pp,FJAC)
            !
            use mkl_rci
            !
            implicit none 
            !
            class(am_class_tb_optimizer)  , intent(in) :: ft
            type(am_class_tight_binding)  , intent(inout) :: tb
            type(am_class_tb_group)       , intent(in) :: tbpg
            type(am_class_bz)             , intent(in) :: bz
            type(am_class_dispersion_vasp), intent(in) :: dr
            type(am_class_prim_pair)      , intent(in) :: pp
            real(dp), allocatable         , intent(out):: FJAC(:,:) ! fjac(m,n) jacobian matrix
            real(dp), allocatable :: f1(:)     ! f1(m)     residual vector
            real(dp), allocatable :: f2(:)     ! f1(m)     residual vector
            real(dp), allocatable :: x(:)      ! x(n)      parameter vector
            real(dp) :: eps_jac
            integer*8:: handle 
            integer  :: successful
            integer  :: rci_request
            ! set epsilon for jacobian calc
            eps_jac = 1.d-5
            ! allocate space
            allocate(  f1(ft%nrs))
            allocate(  f2(ft%nrs))
            allocate(FJAC(ft%nrs,ft%nxs))
            ! initialize x
            allocate(x, source=tb%V)
            ! begin jacobian computation
            if (djacobi_init(handle=handle, n=ft%nxs, m=ft%nrs, x=x, FJAC=FJAC, eps=eps_jac) .ne. tr_success) then
                print *, '#fail: error in djacobi_init' 
                call mkl_free_buffers
                stop 1
            end if
            rci_request = 0
            successful  = 0
            do while (successful.eq.0)
                if (djacobi_solve(handle=handle, f1=f1, f2=f2, rci_request=rci_request) .ne. tr_success) then
                    print *, '#fail: error in djacobi_solve'
                    call mkl_free_buffers
                    stop 1
                end if
                if (rci_request .eq. 1) then
                    ! update x in TB model
                    call tb%set_Vsk(V=x)
                    ! compute jacobian
                    call compute_residual(ft=ft,tb=tb,tbpg=tbpg,bz=bz,dr=dr,pp=pp,R=f1)
                else if (rci_request .eq. 2) then
                    ! update x in TB model
                    call tb%set_Vsk(V=x)
                    ! compute jacobian
                    call compute_residual(ft=ft,tb=tb,tbpg=tbpg,bz=bz,dr=dr,pp=pp,R=f2)
                else if (rci_request .eq. 0) then
                    successful = 1 
                end if
            end do
            if (djacobi_delete(handle) .ne. tr_success) then
                print *, '#fail: error in djacobi_delete' 
                call mkl_free_buffers
                stop 1
            else
                call mkl_free_buffers
            end if
        end subroutine compute_jacobian
    end subroutine optimize_matrix_elements

    subroutine     get_dispersion(dr,tb,tbpg,bz,pp)
        !
        implicit none
        !
        class(am_class_dispersion_tb),intent(out) :: dr
        type(am_class_tb_group)      , intent(in) :: tbpg ! point group in tight binding representation
        type(am_class_tight_binding) , intent(in) :: tb
        type(am_class_bz)            , intent(in) :: bz
        type(am_class_prim_pair)     , intent(in) :: pp
        complex(dp), allocatable :: H(:,:), V(:,:)
        real(dp)   , allocatable :: D(:)
        integer :: i
        !
        ! tight binding coefficients
        allocate(dr%C(tbpg%nbases,tbpg%nbases,bz%nkpts))
        ! eigenvalues
        allocate(dr%E(tbpg%nbases,bz%nkpts))
        ! loop over kpoints
        do i = 1, bz%nkpts
            ! construct hamiltonian
            H = tb%get_hamiltonian(tbpg=tbpg, pp=pp, kpt=bz%kpt_cart(:,i))
            ! diagonalize hamiltonian
            call am_zheev(A=H,V=V,D=D)
            dr%C(:,:,i) = V
            dr%E(:,i) = D
        enddo
        !
    end subroutine get_dispersion




































end module am_tight_binding











