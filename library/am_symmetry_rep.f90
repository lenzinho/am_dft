module am_symmetry_rep

    use am_constants
    use am_stdout
    use am_options
    use am_symmetry
    use am_matlab
    use am_atom
    use am_prim_cell
    use am_irre_cell
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

    type, public, extends(am_class_symrep_group) :: am_class_regular_group
	    contains
        procedure :: get_regular_group
    end type am_class_regular_group

    type, public, extends(am_class_symrep_group) :: am_class_tb_group
        integer, allocatable :: H_start(:) ! Hstart(m) start of hamiltonian section corresponding to atom m
        integer, allocatable :: H_end(:)   ! Hend(m)   end of hamiltonian section corresponding to atom m
        contains
        procedure :: get_tight_binding_point_group
    end type am_class_tb_group

    type, public :: am_class_tensor
        character(100)        :: property       ! name of propertty
        character(100)        :: flags          ! axial/polar/tight
        integer               :: rank           ! tensor rank
        integer , allocatable :: dims(:)        ! tensor dimensions
        real(dp), allocatable :: relations(:,:) ! relations connecting tensor elements
        real(dp), allocatable :: V(:)           ! the value of the tensor, use reshape(V,dims)
                                                ! to apply symmetry operations: V_symmetrized = matmul(relations,V); computes dependent from independent parameters
        contains
    end type am_class_tensor

    type, public, extends(am_class_tensor) :: am_class_property
        contains
        procedure :: get_property
    end type am_class_property

	contains

	! procedures which create representations

	subroutine     get_regular_group(rr,sg)
		!
		implicit none
		!
		class(am_class_regular_group), intent(out) :: rr
        class(am_class_seitz_group)  , intent(in) :: sg
		!
		rr%sym   = rep_regular(seitz_frac=sg%seitz_frac)
		!
		rr%nbases= size(rr%sym,1)
        !
		rr%nsyms = size(rr%sym,3)
		!
        contains
        function       rep_regular(seitz_frac) result(reg_rep)
            !
            ! Generates regular representation for each operator
            !
            ! For details, see page 73-74, especially Eqs. 4.18 Symmetry and Condensed Matter Physics: A Computational
            ! Approach. 1 edition. Cambridge, UK; New York: Cambridge  University Press, 2008.
            !
            ! Requires identity to be the first element of group. regular representation is built from multiplication table
            ! by making sure that the identity appears along the diagonal. this sometimes causes the row to be permuted with
            ! respect to the columns. so, first find the permutation necessary to bring the identities to the center and
            ! then build regular representations from this new multiplication table.
            !
            implicit none
            !
            real(dp), intent(in) :: seitz_frac(:,:,:)
            integer, allocatable :: reg_rep(:,:,:)
            integer, allocatable :: multab(:,:)
            integer :: n
            integer :: i
            !
            ! obtain multiplication table
            multab = get_multab(sym=seitz_frac,flags='seitz,sort')
            !
            ! construct regular representation from multiplication table
            n = size(seitz_frac,3)
            allocate(reg_rep(n,n,n))
            reg_rep = 0
            do i = 1, n
                where (multab.eq.i) reg_rep(:,:,i) = 1
            enddo
            !
        end function   rep_regular
	end subroutine get_regular_group

    ! get point group in the hamilonian tight binding basis 

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
        allocate(tbpg%H_start(pc%natoms))
        allocate(tbpg%H_end(pc%natoms))
        do i = 1, pc%natoms
            tbpg%H_start(i) = tbpg%nbases + 1
            tbpg%nbases     = tbpg%nbases + ic%atom(pc%ic_id(i))%norbitals
            tbpg%H_end(i)   = tbpg%nbases
        enddo
        ! number of symmetries
        tbpg%nsyms = pg%nsyms
        ! generate intrinsic symmetries
        allocate(tbpg%sym(tbpg%nbases,tbpg%nbases,tbpg%nsyms))
        tbpg%sym = 0
        ! determine rotation in the hamiltonian basis
        ! Nye, J.F. "Physical properties of crystals: their representation by tensors and matrices". p 133 Eq 7
        do i = 1, pg%nsyms
            tbpg%sym(:,:,i) = ps2tb_H(R_cart=pg%seitz_cart(1:3,1:3,i), pc=pc, ic=ic)
        enddo
        ! copy symmetry ids
        allocate(tbpg%ps_id, source=pg%ps_id)
        ! copy classes
        allocate(tbpg%cc%id, source=pg%cc%id)
        ! copy multiplication table
        tbpg%mt = pg%mt
        ! copy character table
        tbpg%ct = pg%ct
        ! check that identity is first
        if (.not.isequal(tbpg%sym(:,:,1),eye(tbpg%nbases))) stop 'tbpg: Identity is not first.'
        !
        contains
        function       ps2tb_H(R_cart,pc,ic) result(H)
            ! produces rotation which commutes with the entire Hamiltonian (useful building Hamiltonian and probably later for kpoints stuff too)
            implicit none
            !
            real(dp), intent(in) :: R_cart(3,3)
            type(am_class_prim_cell), intent(in) :: pc ! primitive cell
            type(am_class_irre_cell), intent(in) :: ic ! irreducible cell
            real(dp), allocatable :: H(:,:)
            integer , allocatable :: H_start(:), H_end(:)
            integer :: Hdim ! hamiltonian dimensions
            integer :: i
            !
            ! allocate space for defining subsections of rotation in Hamiltonian basis
            allocate(H_start(pc%natoms))
            allocate(H_end(pc%natoms))
            ! determine subsections of rotations in Hamiltonian basis corresponding to each atom
            Hdim = 0
            do i = 1, pc%natoms
                H_start(i) = Hdim + 1
                Hdim = Hdim + ic%atom(pc%ic_id(i))%norbitals
                H_end(i) = Hdim
            enddo
            ! allocate and initialize space for rotations in Hamiltonin bais
            allocate(H(Hdim,Hdim))
            H = 0.0_dp
            ! construct rotation in the Hamiltonian basis
            do i = 1, pc%natoms
                H(H_start(i):H_end(i), H_start(i):H_end(i)) = ps2tb(R_cart=R_cart, atom=ic%atom(pc%ic_id(i)) )
            enddo
            !
        end function   ps2tb_H
    end subroutine get_tight_binding_point_group

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
            if (.not.present(atom_m)) stop 'atom_m is required for the creation of tb flat intrinsic symmetry group'
            if (.not.present(atom_n)) stop 'atom_n is required for the creation of tb flat intrinsic symmetry group'
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
            stop 'Undefined propery encountered.'
        endif
        !
        allocate(flat_ig%ps_id(flat_ig%nsyms))
        flat_ig%ps_id    = default_ps_id_value
        flat_ig%ps_id(1) = 1
        !
        ! get multiplication table (determines conjugacy classes in the process)
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
        real(dp), allocatable :: R_cart(:,:)
        integer  :: i
        !
        ! basic checks
        if (index(tens%property,'tight').ne.0) then
            if (.not.present(atom_m)) stop 'atom_m is required for the creation of tb flat pg'
            if (.not.present(atom_n)) stop 'atom_n is required for the creation of tb flat pg'
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
            else;
                stop 'tens%flags unknown'
            endif
        enddo
        !
        ! copy symmetry ids
        allocate(flat_pg%ps_id, source=pg%ps_id)
        ! copy classes
        allocate(flat_pg%cc%id, source=pg%cc%id)
        ! copy multiplication table
        flat_pg%mt = pg%mt
        ! copy character table
        flat_pg%ct = pg%ct
        ! get relations
        flat_pg%relations = flat_pg%get_relations()
        !
        if (.not.isequal(flat_pg%sym(:,:,1),eye(flat_pg%nbases))) then
            call am_print('flat_pg%sym(:,:,1)',flat_pg%sym(:,:,1))
            stop 'flat_pg: Identity is not first.'
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
        call rref(A)
        ! correct basic rounding error
        where (abs(nint(A)-A).lt.tiny) A = nint(A)
        ! At this point, A = [ LHS | RHS ], in which LHS = E, identity matrix; A completely specifies all relationships between variables: LHS = RHS.
        allocate(LHS(flat%nbases,flat%nbases))
        allocate(RHS(flat%nbases,flat%nbases))
        LHS = A(0*flat%nbases+indices,0*flat%nbases+indices)
        RHS = A(0*flat%nbases+indices,1*flat%nbases+indices)
        !
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

    ! operates on prop

    subroutine     get_property(prop,pg,opts,property)
        !
        implicit none
        !
        class(am_class_property),    intent(out):: prop
        class(am_class_point_group), intent(in) :: pg
        type(am_class_options),      intent(in) :: opts
        character(*),                intent(in) :: property
        type(am_class_flat_group) :: flat_ig
        type(am_class_flat_group) :: flat_pg
        integer :: m, n, o, nterms
        !
        if (opts%verbosity.ge.1) call print_title('Determining group of '//trim(property))
        !
        ! initialize
        call initialize_property(prop=prop, property=property)
        ! get intrinsic symmetries
        call flat_ig%get_flat_intrinsic_group(tens=prop)
        ! get point symmetries
        call flat_pg%get_flat_point_group(tens=prop, pg=pg)
        ! combined relations
        prop%relations = combine_relations(flat_ig%relations, flat_pg%relations)
        !
        ! print relations
        if (opts%verbosity.ge.1) then
            ! print statistics about parameters
            write(*,'(a5,2a6,3a14)') ' ... ', 'shell', 'terms', 'null', 'dependent', 'independent'
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
            m = count(get_null(prop%relations))
            n = count(get_depenent(prop%relations)) 
            o = count(get_independent(prop%relations))
            nterms = m+n+o
            write(*,'(5x,a6)'        ,advance='no') 'total'
            write(*,'(i6)'           ,advance='no') nterms
            write(*,'(i5,a2,f5.1,a2)',advance='no') m, '(', (m*100_dp)/real(nterms,dp) , '%)'
            write(*,'(i5,a2,f5.1,a2)',advance='no') n, '(', (n*100_dp)/real(nterms,dp) , '%)'
            write(*,'(i5,a2,f5.1,a2)',advance='no') o, '(', (o*100_dp)/real(nterms,dp) , '%)'
            write(*,*)
            ! print symmetry relations
            if (m.ne.nterms) then
                write(*,'(a5,a)') ' ... ', 'irreducible symmetry relations:'
                call print_relations(relations=prop%relations, dims=prop%dims, flags='print:dependent,independent')
            endif
        endif
        !
        contains
        subroutine     initialize_property(prop,property)
            !
            implicit none
            !
            class(am_class_property), intent(out) :: prop
            character(*)            , intent(in)  :: property
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
            prop%property = property
            !------------------------------------------------- FIRST-RANK TENSORS --------------------------------------------------
            if     (index(prop%property,'pyroelectricity')          .ne.0) then; prop%rank = 1                       !    ! P_{i}     = p_{i} \Delta T
            !------------------------------------------------- SECOND-RANK TENSORS -------------------------------------------------
            elseif (index(prop%property,'dielectric')               .ne.0) then; prop%rank = 2; prop%flags = 'polar' 
            elseif (index(prop%property,'electrical susceptibility').ne.0) then; prop%rank = 2; prop%flags = 'polar' ! S  ! P_{i}     = \alpha_{ij}  E_{j}
            elseif (index(prop%property,'magnetic susceptibility')  .ne.0) then; prop%rank = 2; prop%flags = 'axial' ! S  ! M_{i}     = \mu_{ij}     H_{j}
            elseif (index(prop%property,'magneto-electric')         .ne.0) then; prop%rank = 2; prop%flags = 'axial' 
            elseif (index(prop%property,'thermal expansion')        .ne.0) then; prop%rank = 2; prop%flags = 'polar' ! S  ! \eps_{ij} = \alpha_{ij}  \Delta T
            elseif (index(prop%property,'electrical conductivity')  .ne.0) then; prop%rank = 2; prop%flags = 'polar' ! S  ! J_{i}     = \sigma_{ij}  E_{i}
            elseif (index(prop%property,'electrical resistivity')   .ne.0) then; prop%rank = 2; prop%flags = 'polar' ! S  ! E_{i}     = \rho_{ij}    J_{j}
            elseif (index(prop%property,'thermal conductivity')     .ne.0) then; prop%rank = 2; prop%flags = 'polar' ! S  ! q_{i}     = \kappa_{ij}  \frac{\partial T}/{\partial r_{j}}
            elseif (index(prop%property,'thermoelectricity')        .ne.0) then; prop%rank = 2; prop%flags = 'polar' ! N  ! 
            elseif (index(prop%property,'seebeck')                  .ne.0) then; prop%rank = 2; prop%flags = 'polar' ! N  ! E_{i}     = \beta_{ij}   \frac{\partial T}/{\partial r_{j}}
            elseif (index(prop%property,'peltier')                  .ne.0) then; prop%rank = 2; prop%flags = 'polar' ! N  ! q_{i}     = \pi_{ij}     J_{j}
            !------------------------------------------------- THIRD-RANK TENSORS -------------------------------------------------
            elseif (index(prop%property,'hall')                     .ne.0) then; prop%rank = 3;                      !    ! E_{i}     = h_{ijk}      J_{j} H_{k} - has two polar componnts {ij} and one axial component {k}
            elseif (index(prop%property,'piezoelectricity')         .ne.0) then; prop%rank = 3; prop%flags = 'polar' !    ! P_{i}     = d_{ijk}      \sigma_{jk}
            elseif (index(prop%property,'piezomagnetic')            .ne.0) then; prop%rank = 3; prop%flags = 'axial' !    ! M_{i}     = Q_{ijk}      \sigma_{jk}
            !------------------------------------------------- FOURTH-RANK TENSORS ------------------------------------------------
            elseif (index(prop%property,'elasticity')               .ne.0) then; prop%rank = 4; prop%flags = 'polar' !    ! 
            elseif (index(prop%property,'piezo-optic')              .ne.0) then; prop%rank = 4                       !    ! 
            elseif (index(prop%property,'kerr')                     .ne.0) then; prop%rank = 4                       !    ! 
            elseif (index(prop%property,'electrostriction')         .ne.0) then; prop%rank = 4                       !    ! 
            !------------------------------------------------- SXITH-RANK TENSORS -------------------------------------------------
            elseif (index(prop%property,'third-order elasticity')   .ne.0) then; prop%rank = 6; prop%flags = 'polar' !    ! 
            !----------------------------------------------------------------------------------------------------------------------
            else
                stop 'Unknown property.'
            endif
            !----------------------------------------------------------------------------------------------------------------------
            !
            allocate(prop%dims(prop%rank))
            prop%dims = 3
            !
        end subroutine initialize_property
    end subroutine get_property

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

    ! write point symmetry operation in tight binding basis ()

    ! produces rotation which operates on a subsection of the Hamiltonian (useful for determining symmetry relations)

    function       ps2tb(R_cart,atom) result(H)
        ! construct rotation representation which transforms all orbitals on input atom
        implicit none
        !
        real(dp), intent(in) :: R_cart(3,3)
        type(am_class_atom), intent(in) :: atom
        real(dp), allocatable :: H(:,:)
        integer , allocatable :: H_start(:), H_end(:)
        integer :: Hdim ! hamiltonian dimensions
        integer :: i
        !
        ! allocate space for defining subsections of the rotation in Hamiltonian basis
        allocate(H_start(atom%nazimuthals))
        allocate(H_end(atom%nazimuthals))
        ! determine subsections: 1 per primitive atom per azimuthal quantum number per magnetic number per spin
        ! ignoring spin contributions at this stage
        Hdim = 0
        do i = 1, atom%nazimuthals
            H_start(i) = Hdim + 1
            Hdim = Hdim + atom%azimuthal(i)*2+1
            H_end(i) = Hdim
        enddo
        ! allocate and initialize space for rotations in Hamiltonin bais
        allocate(H(Hdim,Hdim))
        H = 0.0_dp
        ! construct rotation in the Hamiltonian basis
        do i = 1, atom%nazimuthals
            H(H_start(i):H_end(i), H_start(i):H_end(i)) = rot2irrep(l=atom%azimuthal(i), R=R_cart)
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
        ! transpose superoperator: (l,l',m) = (-1)^(l+l') (l',l,m)
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









