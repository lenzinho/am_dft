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
    use am_brillouin_zone

    implicit none

    private

    type, public :: am_class_tight_binding
        integer :: nshells                             ! how many shells irreducible atoms
        integer , allocatable :: nVs                   ! number of irreducible matrix element values
        real(dp), allocatable :: V(:)                  ! their values
        integer , allocatable :: V_ind(:,:)            ! their indices
        type(am_class_tensor), allocatable :: tens(:)  ! tens(nshells)
        ! type(am_class_hamiltonian)
        contains
        procedure :: set_Vsk
        procedure :: initialize_tb
        procedure :: get_hamiltonian
        procedure :: test_hamiltonian
        procedure :: export_to_matlab
        procedure :: read_matrix_elements_irr
        procedure :: write_matrix_elements_irr
        procedure :: write_matrix_elements_full
    end type am_class_tight_binding

    type, public, extends(am_class_symrep_group) :: am_class_tb_group
        integer, allocatable :: S(:) ! Hstart(m) start of hamiltonian section corresponding to atom m
        integer, allocatable :: E(:) ! Hend(m)   end of hamiltonian section corresponding to atom m
        contains
        procedure :: get_tight_binding_point_group
    end type am_class_tb_group

contains

    subroutine     get_tight_binding_point_group(tbpg,pg,pc,ic)
        !
        implicit none
        !
        class(am_class_tb_group)   , intent(out) :: tbpg ! point symmetries in tight binding basis
        class(am_class_point_group), intent(in)  :: pg   ! seitz point group (rev stab rot groups as well)
        type(am_class_prim_cell)   , intent(in)  :: pc   ! primitive cell
        type(am_class_irre_cell)   , intent(in)  :: ic   ! irreducible cell
        integer  :: i
        !
        ! get number of bases functions in representation (see below) ...
        tbpg%nbases = 0
        ! ... and also determine subsections of rotations in Hamiltonian basis corresponding to each atom
        allocate(tbpg%S(pc%natoms))
        allocate(tbpg%E(pc%natoms))
        do i = 1, pc%natoms
            tbpg%S(i)   = tbpg%nbases + 1
            tbpg%nbases = tbpg%nbases + ic%atom(pc%ic_id(i))%norbitals
            tbpg%E(i)   = tbpg%nbases
        enddo
        ! number of symmetries
        tbpg%nsyms = pg%nsyms
        ! generate intrinsic symmetries
        allocate(tbpg%sym(tbpg%nbases,tbpg%nbases,tbpg%nsyms))
        tbpg%sym = 0
        ! determine rotation in the hamiltonian basis
        ! Nye, J.F. "Physical properties of crystals: their representation by tensors and matrices". p 133 Eq 7
        do i = 1, pg%nsyms
            ! convert rotation R to symmetry representation in tight-binding basis (direct sum of wigner matrices)
            tbpg%sym(:,:,i) = ps2D(S=tbpg%S, E=tbpg%E, R_cart=pg%seitz_cart(1:3,1:3,i), pc=pc, ic=ic)
        enddo
        ! check that identity is first
        if (.not.isequal(tbpg%sym(:,:,1),eye(tbpg%nbases))) stop 'ERROR [get_tight_binding_point_group]: Identity is not first.'
        ! correct basic rounding errors
        call correct_rounding_error(tbpg%sym)
        ! copy symmetry ids
        allocate(tbpg%ps_id, source=pg%ps_id)
        ! copy classes
        allocate(tbpg%cc%id, source=pg%cc%id)
        ! check that multiplication table is identical to pg%mt
        call tbpg%get_multiplication_table()
        if (.not.isequal(tbpg%mt%multab,pg%mt%multab)) stop 'ERROR [get_tight_binding_point_group]: Multiplication table mismatch.'
        ! get conjugacy clases
        call tbpg%get_conjugacy_classes()
        ! get character table
        call tbpg%get_character_table()
        ! create tb dir
        call execute_command_line('mkdir -p '//trim(outfile_dir_tb))
        ! write point group
        call dump(A=tbpg%sym,fname=trim(outfile_dir_tb)//'/outfile.tbpg.sym')
        ! dump debug files
        if (debug) then
        call execute_command_line('mkdir -p '//trim(debug_dir)//'/tbpg')
        call tbpg%debug_dump(fname=            trim(debug_dir)//'/tbpg'//'/outfile.tbpg')
        endif
        !
        contains
        function       ps2D(S,E,R_cart,pc,ic) result(Dsum)
            ! produces rotation which commutes with the entire Hamiltonian (useful building Hamiltonian and probably later for kpoints stuff too)
            implicit none
            !
            integer , intent(in) :: S(:)
            integer , intent(in) :: E(:)
            real(dp), intent(in) :: R_cart(3,3)
            type(am_class_prim_cell), intent(in) :: pc ! primitive cell
            type(am_class_irre_cell), intent(in) :: ic ! irreducible cell
            real(dp), allocatable :: Dsum(:,:)
            integer :: nbases
            integer :: i
            ! get hamil
            nbases = maxval(E)
            ! allocate and initialize space for rotations in Hamiltonin bais
            allocate(Dsum(nbases,nbases))
            Dsum = 0.0_dp
            ! construct rotation in the Hamiltonian basis
            do i = 1, pc%natoms
                Dsum(S(i):E(i), S(i):E(i)) = ps2tb(R_cart=R_cart, atom=ic%atom(pc%ic_id(i)) )
            enddo
            !
        end function   ps2D
    end subroutine get_tight_binding_point_group

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
        allocate(S, source=tbpg%S)
        allocate(E, source=tbpg%E)
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
        allocate(S, source=tbpg%S)
        allocate(E, source=tbpg%E)
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
        integer, allocatable :: Sv(:) ! start and end of irreducible vector
        integer, allocatable :: Ev(:) ! start and end of irreducible vector
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
        allocate(S, source=tbpg%S)
        allocate(E, source=tbpg%E)
        end = size(tbpg%E)
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




































end module am_tight_binding











