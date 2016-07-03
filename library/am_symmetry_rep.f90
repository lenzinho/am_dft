module am_symmetry_rep

    use am_constants
    use am_stdout
    use am_options
    use am_symmetry
    use am_matlab
    use am_atom
    use am_prim_cell
    use am_irre_cell
    use am_shells
    use am_symmetry_relations
    use dispmodule

    implicit none

    private

    public :: ps2tb ! used in am_tight_binding

    public :: print_relations, combine_relations, dump_relations

    public :: get_null, get_independent, get_depenent

    public :: transp_parity_sign

    type, public, extends(am_class_symrep_group) :: am_class_flat_group
        real(dp), allocatable :: relations(:,:)
        contains
        procedure :: get_flat_intrinsic_group
        procedure :: get_flat_point_group
        procedure :: get_direct_product
        procedure :: get_relations
    end type am_class_flat_group

    type, public :: am_class_tensor
        character(100)        :: property       ! name of propertty
        character(100)        :: flags          ! axial/polar/tight
        integer               :: rank           ! tensor rank
        integer , allocatable :: dims(:)        ! tensor dimensions
        real(dp), allocatable :: relations(:,:) ! relations connecting tensor elements
        real(dp), allocatable :: V(:)           ! the value of the tensor, use reshape(V,dims)
        contains
        procedure :: symmetrize
    end type am_class_tensor

	contains

    subroutine     get_flat_intrinsic_group(flat_ig,tens,atom_m,atom_n)
        !
        implicit none
        !
        class(am_class_flat_group), intent(out) :: flat_ig ! intrinsic symmetry group
        class(am_class_tensor)    , intent(in)  :: tens
        type(am_class_atom)       , intent(in), optional :: atom_m ! only required if property = tb
        type(am_class_atom)       , intent(in), optional :: atom_n ! only required if property = tb
        real(dp), allocatable :: T(:,:) ! transpositional operator building block
        integer :: ndims
        integer :: k
        !
        ! basic checks
        if (index(tens%property, 'tight').ne.0) then
        if (index(tens%flags,'i==j'     ).ne.0) then
            if (.not.present(atom_m)) stop 'ERROR [get_flat_intrinsic_group]: atom_m required'
            if (.not.present(atom_n)) stop 'ERROR [get_flat_intrinsic_group]: atom_n required'
        endif
        endif
        !
        ! number of spatial dimensions
        ndims = 3
        ! get transpositional operator
        T = transp_operator(ndims)
        ! number of bases functions in representation
        flat_ig%nbases = product(tens%dims)
        ! generate intrinsic symmetries
        k=0
        if     (index(tens%property,'thermoelectricity').ne.0 &
         & .or. index(tens%property,'tightbinding'     ).ne.0) then
            !
            flat_ig%nsyms = 1
            allocate(flat_ig%sym(flat_ig%nbases,flat_ig%nbases,flat_ig%nsyms))
            k=k+1; flat_ig%sym(:,:,k) = eye(flat_ig%nbases)    ! E
            !
        elseif (index(tens%property,'conductivity'     ).ne.0 &
         & .or. index(tens%property,'resistivity'      ).ne.0 &
         & .or. index(tens%property,'voigt'            ).ne.0) then
            !
            flat_ig%nsyms = 2
            allocate(flat_ig%sym(flat_ig%nbases,flat_ig%nbases,flat_ig%nsyms))
            k=k+1; flat_ig%sym(:,:,k) = eye(flat_ig%nbases)    ! E
            k=k+1; flat_ig%sym(:,:,k) = T                      ! s_ij = s_ji
            !
        elseif (index(tens%property,'piezoelectricity' ).ne.0) then
            !
            flat_ig%nsyms = 2
            allocate(flat_ig%sym(flat_ig%nbases,flat_ig%nbases,flat_ig%nsyms))
            k=k+1; flat_ig%sym(:,:,k) = eye(flat_ig%nbases)    ! E
            k=k+1; flat_ig%sym(:,:,k) = kron(T,eye(ndims))     ! d_ijk = d_ikj
            !
        elseif (index(tens%property,'elasticity'       ).ne.0) then
            !
            flat_ig%nsyms = 4
            allocate(flat_ig%sym(flat_ig%nbases,flat_ig%nbases,flat_ig%nsyms))
            k=k+1; flat_ig%sym(:,:,k) = eye(flat_ig%nbases)    ! E
            k=k+1; flat_ig%sym(:,:,k) = kron(eye(ndims**2),T)  ! cijkl = cjikl
            k=k+1; flat_ig%sym(:,:,k) = kron(T,eye(ndims**2))  ! cijkl = cjilk
            k=k+1; flat_ig%sym(:,:,k) = kron(T,T)              ! cijkl = cjilk
            !
        elseif (index(tens%property,'tight'            ).ne.0) then
        if     (index(tens%flags   ,'i==j'             ).ne.0) then
            ! atoms correspond to the same irreducible atom
            flat_ig%nsyms = 2
            allocate(flat_ig%sym(flat_ig%nbases,flat_ig%nbases,flat_ig%nsyms))
            k=k+1; flat_ig%sym(:,:,k) = eye(flat_ig%nbases)    ! E
            k=k+1; flat_ig%sym(:,:,k) = orbital_parity(atom_m,atom_n) ! (l,l',m) = (-1)^(l+l') (l',l,m)
        elseif (index(tens%flags   ,'i/=j'             ).ne.0) then
            ! atoms correspond to different irreducible atoms
            flat_ig%nsyms = 1
            allocate(flat_ig%sym(flat_ig%nbases,flat_ig%nbases,flat_ig%nsyms))
            k=k+1; flat_ig%sym(:,:,k) = eye(flat_ig%nbases)    ! E
        endif
        else
            stop 'ERROR [get_flat_intrinsic_group]: undefined propery'
        endif
        ! sometimes for 1-dimensional reps, multiple symmetries can be the same; a unique is added
        ! below to handle such cases. if ignored, it will cause problems later on when attempting to
        ! determine the multiplication table
        flat_ig%sym = unique(flat_ig%sym)
        flat_ig%nsyms = size(flat_ig%sym,3)
        !
        allocate(flat_ig%ps_id(flat_ig%nsyms))
        flat_ig%ps_id    = default_ps_id_value
        flat_ig%ps_id(1) = 1
        ! get multiplication table
        call flat_ig%get_multiplication_table()
        ! get conjugacy classes
        call flat_ig%get_conjugacy_classes()
        ! get character table
        call flat_ig%get_character_table()
        ! get relations
        flat_ig%relations = flat_ig%get_relations()
        !
    end subroutine get_flat_intrinsic_group

    subroutine     get_flat_point_group(flat_pg,tens,pg,atom_m,atom_n)
        !
        class(am_class_flat_group), intent(out):: flat_pg ! flat point group
        class(am_class_tensor)    , intent(in) :: tens    ! tensor
        class(am_class_point_group),intent(in) :: pg      ! seitz point group (rev stab rot groups as well)
        type(am_class_atom)       , intent(in), optional :: atom_m ! only required if property = tb
        type(am_class_atom)       , intent(in), optional :: atom_n ! only required if property = tb
        integer , allocatable :: uinds(:)
        real(dp), allocatable :: R_cart(:,:)
        integer  :: i
        !
        ! basic checks
        if (index(tens%property,'tight').ne.0) then
            if (.not.present(atom_m)) stop 'ERROR [get_flat_point_group]: atom_m required'
            if (.not.present(atom_n)) stop 'ERROR [get_flat_point_group]: atom_n required'
        endif
        !
        ! number of bases functions in representation
        flat_pg%nbases = product(tens%dims)
        ! number of symmetries
        flat_pg%nsyms = pg%nsyms
        ! allocate space for R_cart
        allocate(R_cart(pg%nbases,pg%nbases))
        ! generate intrinsic symmetries
        allocate(flat_pg%sym(flat_pg%nbases,flat_pg%nbases,flat_pg%nsyms))
        flat_pg%sym = 0
        ! Nye, J.F. "Physical properties of crystals: their representation by tensors and matrices". p 133 Eq 7
        do i = 1, pg%nsyms
            R_cart = pg%seitz_cart(1:3,1:3,i)
            ! determine rotation in the basis
            if     (index(tens%flags   ,'axial').ne.0) then; flat_pg%sym(:,:,i) = kron_pow(R_cart, tens%rank) * det(R_cart)
            elseif (index(tens%flags   ,'polar').ne.0) then; flat_pg%sym(:,:,i) = kron_pow(R_cart, tens%rank)
            elseif (index(tens%property,'tight').ne.0) then; flat_pg%sym(:,:,i) = kron(ps2tb(R_cart,atom_m), ps2tb(R_cart,atom_n))
            else
                stop 'ERROR [get_flat_point_group]: unknown flag'
            endif
        enddo
        ! correct basic rounding error
        where (abs(flat_pg%sym).lt.tiny) flat_pg%sym = 0
        ! get unique symmetry indices
        uinds = unique_inds(flat_pg%sym)
        ! update nsyms
        flat_pg%nsyms = size(uinds)
        ! copy unique symmetries
        flat_pg%sym = reallocate(flat_pg%sym(:,:,uinds))
        ! copy symmetry ids
        allocate(flat_pg%ps_id, source=pg%ps_id(uinds))
        ! get multiplication table
        call flat_pg%get_multiplication_table()
        ! get conjugacy classes
        call flat_pg%get_conjugacy_classes()
        ! get character table
        call flat_pg%get_character_table()
        ! get relations
        flat_pg%relations = flat_pg%get_relations()
        ! check that identity is first
        if (.not.isequal(flat_pg%sym(:,:,1),eye(flat_pg%nbases))) then
            stop 'ERROR [get_flat_point_group]: identity is not first.'
        endif
        !
    end subroutine get_flat_point_group

    function       get_relations(flat) result(relations)
        ! 
        implicit none
        !
        class(am_class_flat_group), intent(inout) :: flat
        real(dp), allocatable :: relations(:,:)     !
        integer , allocatable :: class_member(:,:)  ! 
        real(dp), allocatable :: A(:,:)             ! A(2*flat%nbases,2*flat%nbases) - augmented matrix equation
        real(dp), allocatable :: LHS(:,:)           ! LHS(flat%nbases,flat%nbases)   - left hand side of augmented matrix equation (should be identity after reducing to row echlon form)
        real(dp), allocatable :: RHS(:,:)           ! RHS(flat%nbases,flat%nbases)   - right hand side of augmented matrix equation
        integer , allocatable :: indices(:)         ! used for clarity
        integer :: i, j                             ! loop variables
        !
        !
        ! get conjugacy class members
        class_member = id_member(id=flat%cc%id)
        ! intialize indices
        allocate(indices,source=[1:flat%nbases])
        ! initialize augmented workspace matrix A
        allocate(A(3*flat%nbases,2*flat%nbases))
        A = 0
        ! use LU factorization to incorporate, one symmetry at a time, the effect of all symmetries on the flattened basis
        ! enough to loop over class representatives, which are group generators
        do j = 1, size(class_member,1) ! loop over classes
            ! get index of class representative
            i=class_member(j,1)
            ! construct slice of A
            A(2*flat%nbases+indices,0*flat%nbases+indices) = flat%sym(:,:,i)
            A(2*flat%nbases+indices,1*flat%nbases+indices) = eye(flat%nbases)
            ! incorporate symmetry via lu factorization (equivalent to applying rref)
            call lu(A)
        enddo
        ! Apply Gram-Schmidt orthogonalization to obtain A in reduced row echelon form
        A = rref(A)
        ! correct basic rounding error
        where (abs(nint(A)-A).lt.tiny) A = nint(A)
        ! correct basic rounding errors
        ! At this point, A = [ LHS | RHS ], in which LHS = E, identity matrix; A completely specifies all relationships between variables: LHS = RHS.
        allocate(LHS(flat%nbases,flat%nbases))
        allocate(RHS(flat%nbases,flat%nbases))
        LHS = A(0*flat%nbases+indices,0*flat%nbases+indices)
        RHS = A(0*flat%nbases+indices,1*flat%nbases+indices)
        ! checks
        if (.not.isequal(LHS,eye(flat%nbases))) then
            stop 'Failed to reduce matrix to row echlon form.'
        endif
        if (count(get_null(RHS))+count(get_independent(RHS))+count(get_depenent(RHS)).ne.flat%nbases) then
            stop 'Number of null, independent, and dependent terms do not sum to the number of terms.'
        endif
        !
        allocate(relations, source=RHS)
        !
    end function   get_relations

    ! operates on tens

    subroutine     symmetrize(tens,pg,opts,property,pc,ic,shell)
        !                                           ----------- tight binding input
        !
        implicit none
        !
        class(am_class_tensor),      intent(out):: tens
        class(am_class_point_group), intent(in) :: pg
        type(am_class_options),      intent(in) :: opts
        character(*),                intent(in) :: property
        type(am_class_prim_cell)   , optional, intent(in) :: pc   ! primitive cell
        type(am_class_irre_cell)   , optional, intent(in) :: ic    ! required for tight binding input
        type(am_shell_cell)        , optional, intent(in) :: shell ! required for tight binding input
        type(am_class_flat_group) :: flat_ig
        type(am_class_flat_group) :: flat_pg
        character(:), allocatable :: str
        integer :: m, n, o, nterms
        !
        if (opts%verbosity.ge.1) call print_title('Symmeterized '//trim(property))
        ! check intput for tight binding
        if    (index(property,'tight binding').ne.0) then
            ! tb matrix element
            if (.not.present(pc))    stop 'ERROR [symmetrize]: missing pc for tight-binding intput'
            if (.not.present(ic))    stop 'ERROR [symmetrize]: missing ic for tight-binding intput'
            if (.not.present(shell)) stop 'ERROR [symmetrize]: missing shell for tight-binding intput'
            ! initialize
            call initialize_tensor(tens=tens, pc=pc, ic=ic, shell=shell, property=property)
            ! determine intrinsic symmetries (if both irreducible atoms are of the same irreducible type)
            ! Interchange of indices: (l,l',m) = (-1)^(l+l') * (l',l,m), due to parity of wavefunction. (s,d are even under inversion, p,f are odd)
            ! E. Scheer, Molecular Electronics: An Introduction to Theory and Experiment, p 245. Also see R. Martin.
            ! NOTE: FOR SOME REASON, MUST CALL atom_m with shell%j and atom_n with shell%i (INDICES FLIPPED) - probably has to do with reshape and column-major ordering
            call flat_ig%get_flat_intrinsic_group(tens=tens, atom_m=ic%atom(shell%j), atom_n=ic%atom(shell%i) )
            ! get stabilizer symmetries in the flattened hamiltonin basis
            call flat_pg%get_flat_point_group(tens=tens, pg=shell%stab, atom_m=ic%atom(shell%j), atom_n=ic%atom(shell%i))
            ! combined relations
            tens%relations = combine_relations(flat_ig%relations, flat_pg%relations)
            ! correct rounding error
            call correct_rounding_error(tens%relations)
            ! print stabilzier group parameters
            if (opts%verbosity.ge.1) then
                write(*,'(a,a)') flare, 'atomic pair = '//trim(atm_symb(pc%Z( shell%i )))//'-'//trim(atm_symb(pc%Z( shell%j )))
                write(*,'(a,a)') flare, 'stabilizer group = '//trim(decode_pointgroup(point_group_schoenflies( shell%stab%ps_id )))
                write(*,'(a,a)') flare, 'orbit = '
                call disp_indent()
                call disp(X=[1:shell%natoms]         ,title='#'     ,style='underline',advance='no')
                call disp(X=transpose(shell%tau_cart),title='atomic basis [cart]',style='underline',advance='no')
                call disp(X=transpose(shell%tau_frac),title='atomic basis [frac]',style='underline',advance='yes')
            endif
        else
            ! initialize
            call initialize_tensor(tens=tens, property=property)
            ! get intrinsic symmetries
            call flat_ig%get_flat_intrinsic_group(tens=tens)
            ! get point symmetries
            call flat_pg%get_flat_point_group(tens=tens, pg=pg)
            ! combined relations
            tens%relations = combine_relations(flat_ig%relations, flat_pg%relations)
            ! correct rounding error
            call correct_rounding_error(tens%relations)
            ! dump
            if (debug) then
                allocate(str, source = strrep(property,' ','_') )
                call execute_command_line ('mkdir -p '//trim(outfile_dir_sym)//'/debug/'//str)
                call flat_pg%debug_dump(fname=          trim(outfile_dir_sym)//'/debug/'//str//'/outfile.'//str)
            endif
        endif
        ! print relations
        if (opts%verbosity.ge.1) then
            ! print statistics about parameters
            write(*,'(a,a)') flare, 'matrix elements = '
            write(*,'(5x,2a6,3a14)') 'shell', 'terms', 'null', 'dependent', 'independent'
            write(*,'(5x,a)') repeat(' '//repeat('-',5),2)//repeat(' '//repeat('-',13),3)
            ! intrinsic symmetries
            m = count(get_null(flat_ig%relations))
            n = count(get_depenent(flat_ig%relations)) 
            o = count(get_independent(flat_ig%relations))
            nterms = m+n+o
            write(*,'(5x,a6)'        ,advance='no') 'intr.'
            write(*,'(i6)'           ,advance='no') nterms
            write(*,'(i5,a2,f5.1,a2)',advance='no') m, '(', (m*100_dp)/real(nterms,dp) , '%)'
            write(*,'(i5,a2,f5.1,a2)',advance='no') n, '(', (n*100_dp)/real(nterms,dp) , '%)'
            write(*,'(i5,a2,f5.1,a2)',advance='no') o, '(', (o*100_dp)/real(nterms,dp) , '%)'
            write(*,*)
            ! point symmetries
            m = count(get_null(flat_pg%relations))
            n = count(get_depenent(flat_pg%relations)) 
            o = count(get_independent(flat_pg%relations))
            nterms = m+n+o
            write(*,'(5x,a6)'        ,advance='no') 'pnt.'
            write(*,'(i6)'           ,advance='no') nterms
            write(*,'(i5,a2,f5.1,a2)',advance='no') m, '(', (m*100_dp)/real(nterms,dp) , '%)'
            write(*,'(i5,a2,f5.1,a2)',advance='no') n, '(', (n*100_dp)/real(nterms,dp) , '%)'
            write(*,'(i5,a2,f5.1,a2)',advance='no') o, '(', (o*100_dp)/real(nterms,dp) , '%)'
            write(*,*)
            ! total
            m = count(get_null(tens%relations))
            n = count(get_depenent(tens%relations)) 
            o = count(get_independent(tens%relations))
            nterms = m+n+o
            write(*,'(5x,a6)'        ,advance='no') 'total'
            write(*,'(i6)'           ,advance='no') nterms
            write(*,'(i5,a2,f5.1,a2)',advance='no') m, '(', (m*100_dp)/real(nterms,dp) , '%)'
            write(*,'(i5,a2,f5.1,a2)',advance='no') n, '(', (n*100_dp)/real(nterms,dp) , '%)'
            write(*,'(i5,a2,f5.1,a2)',advance='no') o, '(', (o*100_dp)/real(nterms,dp) , '%)'
            write(*,*)
            ! print symmetry relations
            if (m.ne.nterms) then
                write(*,'(a,a)') flare, 'irreducible symmetry relations:'
                call print_relations(relations=tens%relations, dims=tens%dims, flags='print:dependent,independent')
            endif
            ! show block structure
            write(*,'(a,a)') flare, 'irreducible block structure [flattened point group, ignores intrinsic symmetries]:'
            call print_blocks(sym=flat_pg%sym,block_proj=flat_pg%ct%block_proj)
        endif
        !
        contains
        subroutine     initialize_tensor(tens,property,pc,ic,shell)
            !
            implicit none
            !
            class(am_class_tensor)  , intent(out) :: tens
            character(*)            , intent(in)  :: property
            type(am_class_prim_cell), intent(in), optional :: pc    ! required for tb
            type(am_class_irre_cell), intent(in), optional :: ic    ! required for tb
            type(am_shell_cell)     , intent(in), optional :: shell ! required for tb
            integer :: i
            !
            !----------------------------------------------------------------------------------------------------------------------
            !
            ! Tensors transferm differently whether they are axial or polar:
            ! POLAR ( odd under inv.):  T_{i,j,k,l,...} =        sum_{ip,jp,kp,lp,...} R_{i,ip} R_{j,jp} R_{k,kp} R_{l,lp} T_{ip,jp,kp,lp,...}
            ! AXIAL (even under inv.):  T_{i,j,k,l,...} = det(R) sum_{ip,jp,kp,lp,...} R_{i,ip} R_{j,jp} R_{k,kp} R_{l,lp} T_{ip,jp,kp,lp,...}
            ! Thus, all axial tensors of even rank and polar tensors of odd rank are null are null.  Wooten p 485. Eq. 13.21. 
            !
            ! Onsager’s Principle requires that the electric resistivity and thermal conductivity tensors be symmetric.
            ! This does not hold for the Seebeck and Peltier (thermoelectric) tensors which relate two different flows. Thus
            ! there are, at most, nine independent parameters rather than six. [Newnham "Properties of Materials"]
            !
            tens%property = property
            !------------------------------------------------- FIRST-RANK TENSORS --------------------------------------------------
            if     (index(tens%property,'pyroelectricity')          .ne.0) then; tens%rank = 1                       !    ! P_{i}     = p_{i} \Delta T
            !------------------------------------------------- SECOND-RANK TENSORS -------------------------------------------------
            elseif (index(tens%property,'dielectric')               .ne.0) then; tens%rank = 2; tens%flags = 'polar' 
            elseif (index(tens%property,'electrical susceptibility').ne.0) then; tens%rank = 2; tens%flags = 'polar' ! S  ! P_{i}     = \alpha_{ij}  E_{j}
            elseif (index(tens%property,'magnetic susceptibility')  .ne.0) then; tens%rank = 2; tens%flags = 'axial' ! S  ! M_{i}     = \mu_{ij}     H_{j}
            elseif (index(tens%property,'magneto electric')         .ne.0) then; tens%rank = 2; tens%flags = 'axial' 
            elseif (index(tens%property,'thermal expansion')        .ne.0) then; tens%rank = 2; tens%flags = 'polar' ! S  ! \eps_{ij} = \alpha_{ij}  \Delta T
            elseif (index(tens%property,'electrical conductivity')  .ne.0) then; tens%rank = 2; tens%flags = 'polar' ! S  ! J_{i}     = \sigma_{ij}  E_{i}
            elseif (index(tens%property,'electrical resistivity')   .ne.0) then; tens%rank = 2; tens%flags = 'polar' ! S  ! E_{i}     = \rho_{ij}    J_{j}
            elseif (index(tens%property,'thermal conductivity')     .ne.0) then; tens%rank = 2; tens%flags = 'polar' ! S  ! q_{i}     = \kappa_{ij}  \frac{\partial T}/{\partial r_{j}}
            elseif (index(tens%property,'thermoelectricity')        .ne.0) then; tens%rank = 2; tens%flags = 'polar' ! N  ! 
            elseif (index(tens%property,'seebeck')                  .ne.0) then; tens%rank = 2; tens%flags = 'polar' ! N  ! E_{i}     = \beta_{ij}   \frac{\partial T}/{\partial r_{j}}
            elseif (index(tens%property,'peltier')                  .ne.0) then; tens%rank = 2; tens%flags = 'polar' ! N  ! q_{i}     = \pi_{ij}     J_{j}
            !------------------------------------------------- THIRD-RANK TENSORS -------------------------------------------------
            elseif (index(tens%property,'hall')                     .ne.0) then; tens%rank = 3;                      !    ! E_{i}     = h_{ijk}      J_{j} H_{k} - has two polar componnts {ij} and one axial component {k}
            elseif (index(tens%property,'piezoelectricity')         .ne.0) then; tens%rank = 3; tens%flags = 'polar' !    ! P_{i}     = d_{ijk}      \sigma_{jk}
            elseif (index(tens%property,'piezomagnetic')            .ne.0) then; tens%rank = 3; tens%flags = 'axial' !    ! M_{i}     = Q_{ijk}      \sigma_{jk}
            !------------------------------------------------- FOURTH-RANK TENSORS ------------------------------------------------
            elseif (index(tens%property,'elasticity')               .ne.0) then; tens%rank = 4; tens%flags = 'polar' !    ! 
            elseif (index(tens%property,'piezo optic')              .ne.0) then; tens%rank = 4                       !    ! 
            elseif (index(tens%property,'kerr')                     .ne.0) then; tens%rank = 4                       !    ! 
            elseif (index(tens%property,'electrostriction')         .ne.0) then; tens%rank = 4                       !    ! 
            !------------------------------------------------- SXITH-RANK TENSORS -------------------------------------------------
            elseif (index(tens%property,'third-order elasticity')   .ne.0) then; tens%rank = 6; tens%flags = 'polar' !    ! 
            !----------------------------------------------------------------------------------------------------------------------
            elseif (index(tens%property,'tight binding')            .ne.0) then
                ! set rank
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
            else
                stop 'ERROR [initialize_tensor]: unknown property'
            endif
            !----------------------------------------------------------------------------------------------------------------------
            !
            ! if property is macroscopic enter if condition
            ! if tight binding, skip it.
            if (index(tens%property,'tight binding').eq.0) then
                allocate(tens%dims(tens%rank))
                tens%dims = 3
            endif
            !
        end subroutine initialize_tensor
    end subroutine symmetrize

    ! other stuff

    subroutine     get_direct_product(C,A,B)
        !
        implicit none
        !
        class(am_class_flat_group) :: C
        type(am_class_flat_group)  :: A
        type(am_class_flat_group)  :: B
        real(dp), allocatable :: wkr(:,:,:)
        real(dp), allocatable :: try(:,:)
        integer :: i, j, k
        !
        if (A%nbases/=B%nbases) stop 'dimension mismatch'
        !
        C%nbases = A%nbases
        !
        allocate(C%ps_id(A%nsyms*B%nsyms))
        C%ps_id = default_ps_id_value ! default
        !
        allocate(try(C%nbases,C%nbases))
        allocate(wkr(C%nbases,C%nbases,A%nsyms*B%nsyms))
        !
        k = 0
        do i = 1, A%nsyms
        do j = 1, B%nsyms
            !
            k = k + 1
            wkr(:,:,k) = matmul(A%sym(:,:,i),B%sym(:,:,j))
            !
            if(mod(k,10).eq.0) call show_progress(iteration=k, maximum=A%nsyms*B%nsyms)
            !
            ! ps_id has very little meaning here
            if (A%ps_id(i).ne.0) C%ps_id(k) = A%ps_id(i)
            if (B%ps_id(j).ne.0) C%ps_id(k) = B%ps_id(j)
            !
        end do
        end do
        !
        allocate(C%sym,source=unique(wkr))
        C%nsyms = size(C%sym,3)
        !
        ! get multiplication table (and conjugacy classes in the process)
        call C%get_multiplication_table()
        !
        ! sort symmetries based on parameters
        call C%sort_symmetries(criterion=real(C%ps_id,dp), flags='acsend')
        call C%sort_symmetries(criterion=real(C%cc%id,dp), flags='ascend')
        !
        ! get character table
        call C%get_character_table()
        !
    end subroutine get_direct_product

    ! produces rotation which operates on a subsection of the Hamiltonian (useful for determining symmetry relations)

    function       ps2tb(R_cart,atom) result(H)
        ! construct rotation representation which transforms all orbitals on input atom
        implicit none
        !
        real(dp), intent(in) :: R_cart(3,3)
        type(am_class_atom), intent(in) :: atom
        real(dp), allocatable :: H(:,:)
        integer , allocatable :: S(:), E(:)
        integer :: Hdim ! hamiltonian dimensions
        integer :: i
        !
        ! allocate space for defining subsections of the rotation in Hamiltonian basis
        allocate(S(atom%nazimuthals))
        allocate(E(atom%nazimuthals))
        ! determine subsections: 1 per primitive atom per azimuthal quantum number per magnetic number per spin
        ! ignoring spin contributions at this stage
        Hdim = 0
        do i = 1, atom%nazimuthals
            S(i) = Hdim + 1
            Hdim = Hdim + atom%azimuthal(i)*2+1
            E(i) = Hdim
        enddo
        ! allocate and initialize space for rotations in Hamiltonin bais
        allocate(H(Hdim,Hdim))
        H = 0.0_dp
        ! construct rotation in the Hamiltonian basis
        do i = 1, atom%nazimuthals
            H(S(i):E(i), S(i):E(i)) = rot2irrep(l=atom%azimuthal(i), R=R_cart)
        enddo
        !
    end function   ps2tb

    ! super operators

    pure function  transp_operator(n) result(M_superop)
        ! transposition superoperator (c_ij -> c_ji), flat rep
        implicit none
        !
        integer, intent(in) :: n
        integer, allocatable :: M_superop(:,:)
        integer  :: i, j
        !
        allocate(M_superop(n**2,n**2))
        M_superop=0
        !
        do i = 1, n
        do j = 1, n
           M_superop(i+n*(j-1),j+n*(i-1)) = 1
        enddo
        enddo
        !
    end function   transp_operator

    pure function  transp_parity_sign(atom_m,atom_n) result(S)
        ! sign under : (l,l',m) = (-1)^(l+l') (l',l,m)
        implicit none
        !
        type(am_class_atom), intent(in) :: atom_m
        type(am_class_atom), intent(in) :: atom_n
        integer, allocatable :: S(:,:) ! sign
        integer :: alpha,beta
        !
        allocate(S(atom_m%norbitals, atom_n%norbitals))
        S = 0
        do alpha = 1, atom_m%norbitals
        do beta  = 1, atom_n%norbitals
            S(alpha,beta) = (-1)**( atom_m%orbital(2,alpha) + atom_n%orbital(2,beta) )
        enddo
        enddo
        !
    end function   transp_parity_sign

    function       orbital_parity(atom_m,atom_n) result(L_superop)
        !
        ! Transpose superoperator: (l,l',m) = (-1)^(l+l') (l',l,m)
        !
        ! For example: 
        !        p-1 p 0 p+1 
        !                    d-2 d-1 d 0 d+1 d+2
        !  p-1  [  1,  2,  3, -4, -5, -6, -7, -8]                               [  1,  2,  3, -4, -5, -6, -7, -8]
        !  p 0  [  9, 10, 11,-12,-13,-14,-15,-16]             ---->             [  9, 10, 11,-12,-13,-14,-15,-16]
        !  p+1  [ 17, 18, 19,-20,-21,-22,-23,-24]         transposition         [ 17, 18, 19,-20,-21,-22,-23,-24]
        !  d-2  [-25,-26,-27, 28, 29, 30, 31, 32]         and sign flip         [-25,-26,-27, 28, 29, 30, 31, 32]
        !  d-1  [-33,-34,-35, 36, 37, 38, 39, 40]                               [-33,-34,-35, 36, 37, 38, 39, 40]
        !  d 0  [-41,-42,-43, 44, 45, 46, 47, 48]   (l,l',m) =                  [-41,-42,-43, 44, 45, 46, 47, 48]
        !  d+1  [-49,-50,-51, 52, 53, 54, 55, 56]      (-1)^(l+l') * (l',l,m)   [-49,-50,-51, 52, 53, 54, 55, 56]
        !  d+2  [-57,-58,-59, 60, 61, 62, 63, 64]                               [-57,-58,-59, 60, 61, 62, 63, 64]
        !
        implicit none
        !
        type(am_class_atom), intent(in) :: atom_m
        type(am_class_atom), intent(in) :: atom_n
        integer, allocatable :: L_superop(:,:)
        integer, allocatable :: A(:,:)
        integer, allocatable :: A_flat(:)
        integer, allocatable :: Ap(:,:)
        integer, allocatable :: Ap_flat(:)
        integer :: nflats
        integer :: i
        !
        nflats = atom_m%norbitals * atom_n%norbitals
        !
        allocate(A_flat(nflats))
        A_flat=[1:nflats]
        !
        allocate(A(atom_m%norbitals, atom_n%norbitals))
        A = reshape(A_flat,[atom_m%norbitals, atom_n%norbitals])
        !
        allocate(Ap(atom_n%norbitals, atom_m%norbitals))
        Ap = transpose(A) * transp_parity_sign(atom_m,atom_n)
        !
        allocate(Ap_flat(nflats))
        Ap_flat = reshape(Ap,[nflats])
        !
        ! construct flat superoperator
        allocate(L_superop(nflats,nflats))
        L_superop = 0
        do i = 1, nflats
            L_superop( A_flat(i) ,  abs(Ap_flat(i))  ) = sign( 1, Ap_flat(i) )
        enddo
        ! check to make sure things are correct
        if (.not.isequal(Ap_flat,matmul(L_superop,A_flat))) stop 'A /= Ap'
        !
    end function   orbital_parity

end module am_symmetry_rep









