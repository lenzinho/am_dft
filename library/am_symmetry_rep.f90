module am_symmetry_rep

    use am_constants
    use am_stdout
    use am_mkl
    use am_options
    use am_symmetry
    use am_matlab
    use am_atom
    use am_prim_cell
    use am_irre_cell
    use am_symmetry_tables , only : get_multab

    implicit none

    private

    public :: ps2tb ! used in am_tight_binding

    public :: print_relations, combine_relations, dump_relations

    public :: get_null, get_independent, get_depenent

    type, public, extends(am_class_group) :: am_class_flat_group
        real(dp), allocatable :: relations(:,:)
        contains
        procedure :: get_flat_intrinsic_group
        procedure :: get_flat_point_group
        procedure :: get_direct_product
        procedure :: get_relations
    end type am_class_flat_group

    type, public, extends(am_class_group) :: am_class_regular_group
	    contains
        procedure :: create => create_regular
    end type am_class_regular_group

    type, public :: am_class_tensor
        character(100)        :: property       ! name of propertty
        character(100)        :: flags          ! axial/polar
        integer               :: rank           ! tensor rank
        integer , allocatable :: dims(:)        ! tensor dimensions
        real(dp), allocatable :: relations(:,:) ! relations connecting tensor elements
        real(dp), allocatable :: V(:)           ! the value of the tensor, use reshape(V,dims)
                                                ! to apply symmetry operations: V_symmetrized = matmul(relations,V); computes dependent from independent parameters
    end type am_class_tensor

    type, public, extends(am_class_tensor) :: am_class_property
        type(am_class_flat_group) :: flat_ig        ! flattened intrinsic symmetry group
        type(am_class_flat_group) :: flat_pg        ! flattened point group
        contains
        procedure :: get_property
    end type am_class_property

	contains

	! procedures which create representations

	subroutine     create_regular(rr,seitz)
		!
		implicit none
		!
		class(am_class_regular_group), intent(out) :: rr
        real(dp), intent(in) :: seitz(:,:,:)
		!
		rr%sym   = rep_regular(seitz=seitz)
		!
		rr%nbases= size(rr%sym,1)
        !
		rr%nsyms = size(rr%sym,3)
		!
        contains
        function       rep_regular(seitz) result(reg_rep)
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
            real(dp), intent(in) :: seitz(:,:,:)
            integer, allocatable :: reg_rep(:,:,:)
            integer, allocatable :: multab(:,:)
            integer :: n
            integer :: i
            !
            ! obtain multiplication table
            multab = get_multab(sym=seitz,flags='seitz,sort')
            !
            ! construct regular representation from multiplication table
            n = size(seitz,3)
            allocate(reg_rep(n,n,n))
            reg_rep = 0
            do i = 1, n
                where (multab.eq.i) reg_rep(:,:,i) = 1
            enddo
            !
        end function   rep_regular
	end subroutine create_regular

    ! operators on tensors

    subroutine     get_flat_intrinsic_group(flat_ig,tens)
        !
        implicit none
        !
        class(am_class_flat_group), intent(out) :: flat_ig ! intrinsic symmetry group
        class(am_class_tensor)    , intent(in)  :: tens
        real(dp), allocatable :: T(:,:) ! transpositional operator building block
        integer :: ndims
        integer :: k
        !
        ! number of spatial dimensions
        ndims = 3
        ! get transpositional operator
        T = transp_operator(ndims)
        ! number of bases functions in representation
        flat_ig%nbases = ndims**tens%rank
        ! generate intrinsic symmetries
        k=0
        if     (index(tens%property,'thermoelectricity').ne.0 &
         & .or. index(tens%property,'tightbinding'     ).ne.0) then
            !
            flat_ig%nsyms = 1
            allocate(flat_ig%sym(flat_ig%nbases,flat_ig%nbases,flat_ig%nsyms))
            k=k+1; flat_ig%sym(:,:,k) = eye(flat_ig%nbases)        ! E
            !
        elseif (index(tens%property,'conductivity'     ).ne.0 &
         & .or. index(tens%property,'resistivity'      ).ne.0 &
         & .or. index(tens%property,'voigt'            ).ne.0) then
            !
            flat_ig%nsyms = 2
            allocate(flat_ig%sym(flat_ig%nbases,flat_ig%nbases,flat_ig%nsyms))
            k=k+1; flat_ig%sym(:,:,k) = eye(flat_ig%nbases)        ! E
            k=k+1; flat_ig%sym(:,:,k) = T                      ! s_ij = s_ji
            !
        elseif (index(tens%property,'piezoelectricity' ).ne.0) then
            !
            flat_ig%nsyms = 2
            allocate(flat_ig%sym(flat_ig%nbases,flat_ig%nbases,flat_ig%nsyms))
            k=k+1; flat_ig%sym(:,:,k) = eye(flat_ig%nbases)        ! E
            k=k+1; flat_ig%sym(:,:,k) = kron(T,eye(ndims))     ! d_ijk = d_ikj
            !
        elseif (index(tens%property,'elasticity'       ).ne.0) then
            !
            flat_ig%nsyms = 4
            allocate(flat_ig%sym(flat_ig%nbases,flat_ig%nbases,flat_ig%nsyms))
            k=k+1; flat_ig%sym(:,:,k) = eye(flat_ig%nbases)        ! E
            k=k+1; flat_ig%sym(:,:,k) = kron(eye(ndims**2),T)  ! cijkl = cjikl
            k=k+1; flat_ig%sym(:,:,k) = kron(T,eye(ndims**2))  ! cijkl = cjilk
            k=k+1; flat_ig%sym(:,:,k) = kron(T,T)              ! cijkl = cjilk
            !
        else
            stop 'Undefined propery encountered.'
        endif
        !
        allocate(flat_ig%ps_id(flat_ig%nsyms))
        flat_ig%ps_id    = default_ps_id_value
        flat_ig%ps_id(1) = 1
        !
        ! get multiplication table (determines determine conjugacy classes in the process)
        call flat_ig%get_multiplication_table()
        !
        call flat_ig%sort_symmetries(criterion=real(flat_ig%ps_id,dp), flags='acsend')
        call flat_ig%sort_symmetries(criterion=real(flat_ig%cc%id,dp), flags='ascend')
        !
        ! get character table
        call flat_ig%get_character_table()
        !
        contains
        pure function transp_operator(n) result(M)
            ! transposition operator (c_ij -> c_ji) in the flattened representation
            implicit none
            !
            integer, intent(in) :: n
            integer, allocatable :: M(:,:)
            integer  :: i, j
            !
            allocate(M(n**2,n**2))
            M=0
            !
            do i = 1, n
            do j = 1, n
               M(i+n*(j-1),j+n*(i-1)) = 1
            enddo
            enddo
        end function  transp_operator
    end subroutine get_flat_intrinsic_group

    subroutine     get_flat_point_group(flat_pg,tens,pg,pc,atom_m,atom_n)
        !
        class(am_class_flat_group), intent(out):: flat_pg  ! flat point group
        class(am_class_tensor)    , intent(in) :: tens ! tensor
        class(am_class_group)     , intent(in) :: pg   ! seitz point group (rev stab rot groups as well)
        type(am_class_prim_cell)  , intent(in) :: pc   ! primitive cell
        type(am_class_atom)       , intent(in), optional :: atom_m ! only required if property = tb
        type(am_class_atom)       , intent(in), optional :: atom_n ! only required if property = tb
        real(dp), allocatable :: R_cart(:,:)
        integer  :: i
        logical  :: is_seitz
        !
        ! basic checks
        if (index(tens%flags,'tight').ne.0) then
            if (.not.present(atom_m)) stop 'atom_m is required for the creation of tb flat pg'
            if (.not.present(atom_n)) stop 'atom_n is required for the creation of tb flat pg'
        endif
        !
        ! set type
        select type (pg)
        class is (am_class_seitz_group)
            is_seitz = .true.
        class default
            is_seitz = .false.
        end select
        !
        ! number of bases functions in representation
        flat_pg%nbases = product(tens%dims)
        ! number of symmetries
        flat_pg%nsyms = pg%nsyms
        ! generate intrinsic symmetries
        allocate(flat_pg%sym(flat_pg%nbases,flat_pg%nbases,flat_pg%nsyms))
        flat_pg%sym = 0
        ! Nye, J.F. "Physical properties of crystals: their representation by tensors and matrices". p 133 Eq 7
        do i = 1, pg%nsyms
            ! convert to cartesian
            if (is_seitz) then
                R_cart = ps_frac2cart(R_frac=pg%sym(1:3,1:3,i),bas=pc%bas)
            else
                R_cart = pg%sym(:,:,i)
            endif
            ! determine rotation in the basis
            if     (index(tens%flags,'axial').ne.0) then; flat_pg%sym(:,:,i) = kron_pow(R_cart, tens%rank) * det(R_cart)
            elseif (index(tens%flags,'polar').ne.0) then; flat_pg%sym(:,:,i) = kron_pow(R_cart, tens%rank)
            elseif (index(tens%flags,'tight').ne.0) then; flat_pg%sym(:,:,i) = kron(ps2tb(R_cart,atom_m), ps2tb(R_cart,atom_n))
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
        !
        if (.not.isequal(flat_pg%sym(:,:,1),eye(flat_pg%nbases))) stop 'flat_pg: Identity is not first.'
        !
    end subroutine get_flat_point_group

    ! operates on prop

    subroutine     get_property(prop,pc,pg,opts,property)
        !
        implicit none
        !
        class(am_class_property),    intent(out):: prop
        class(am_class_prim_cell),   intent(in) :: pc
        class(am_class_point_group), intent(in) :: pg
        type(am_class_options),      intent(in) :: opts
        character(*),                intent(in) :: property
        !
        if (opts%verbosity.ge.1) call am_print_title('Determining group of '//trim(property))
        !
        ! initialize
        call initialize_property(prop=prop, property=property)
        !
        ! get intrinsic symmetries
        call prop%flat_ig%get_flat_intrinsic_group(tens=prop)
        ! write number of intrinsic symmetries
        if (opts%verbosity.ge.1) call am_print('intrinsict symmetries', prop%flat_ig%nsyms)
        ! get symmetry relations
        prop%flat_ig%relations = prop%flat_ig%get_relations()
        ! print relations
        call print_relations(relations=prop%flat_ig%relations,flags='')
        !
        ! get point symmetries
        call prop%flat_pg%get_flat_point_group(tens=prop, pg=pg, pc=pc)
        ! write number of point symmetries
        if (opts%verbosity.ge.1) call am_print('point symmetries', prop%flat_pg%nsyms)
        ! get symmetry relations
        prop%flat_pg%relations = prop%flat_pg%get_relations()
        ! print relations
        call print_relations(relations=prop%flat_pg%relations,flags='')
        !
        ! combine point and intrisnci symmetries
        if (opts%verbosity.ge.1) call am_print('total symmetries', prop%flat_pg%nsyms * prop%flat_ig%nsyms )
        ! combined relations
        prop%relations = combine_relations(prop%flat_ig%relations, prop%flat_pg%relations)
        ! print relations
        call print_relations(relations=prop%relations, dims=prop%dims, flags='show:dependent,independent')
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
        call C%sort_symmetries(criterion=real(C%ps_id,dp)   , flags='acsend')
        call C%sort_symmetries(criterion=real(C%cc%id,dp), flags='ascend')
        !
        ! get character table
        call C%get_character_table()
        !
    end subroutine get_direct_product

    ! operates on relations

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
        class_member = member(id=flat%cc%id)
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

    function       combine_relations(relationsA,relationsB,relationsC) result(relations)
        !
        ! combins up to three relations
        !
        implicit none
        !
        real(dp), intent(in) :: relationsA(:,:)
        real(dp), intent(in) :: relationsB(:,:)
        real(dp), intent(in), optional :: relationsC(:,:)
        real(dp), allocatable :: relations(:,:)     !
        real(dp), allocatable :: A(:,:)             ! A(2*flat%nbases,2*flat%nbases) - augmented matrix equation
        real(dp), allocatable :: LHS(:,:)           ! LHS(flat%nbases,flat%nbases)   - left hand side of augmented matrix equation (should be identity after reducing to row echlon form)
        real(dp), allocatable :: RHS(:,:)           ! RHS(flat%nbases,flat%nbases)   - right hand side of augmented matrix equation
        integer , allocatable :: indices(:)         ! used for clarity
        integer :: nbases
        !
        !
        if (size(relationsA,1).ne.size(relationsA,2)) stop 'Dimension mismatch: A1 vs A2'
        if (size(relationsA,1).ne.size(relationsB,1)) stop 'Dimension mismatch: A1 vs B1'
        if (size(relationsA,1).ne.size(relationsB,2)) stop 'Dimension mismatch: A1 vs B2'
        if (present(relationsC)) then
        if (size(relationsA,1).ne.size(relationsC,1)) stop 'Dimension mismatch: A1 vs C1'
        if (size(relationsA,1).ne.size(relationsC,2)) stop 'Dimension mismatch: A1 vs C2'
        endif
        !
        nbases = size(relationsA,1)
        !
        ! get indices
        allocate(indices,source=[1:nbases])
        ! initialize augmented workspace matrix A
        allocate(A(3*nbases,2*nbases))
        A = 0
        ! construct slice of A
        A(0*nbases+indices,1*nbases+indices) = relationsA
        A(0*nbases+indices,0*nbases+indices) = eye(nbases)
        A(1*nbases+indices,1*nbases+indices) = relationsB
        A(1*nbases+indices,0*nbases+indices) = eye(nbases)
        if (present(relationsC)) then
        A(2*nbases+indices,1*nbases+indices) = relationsC
        A(2*nbases+indices,0*nbases+indices) = eye(nbases)
        endif
        ! incorporate symmetry via lu factorization (equivalent to applying rref)
        call lu(A)
        ! Apply Gram-Schmidt orthogonalization to obtain A in reduced row echelon form
        call rref(A)
        ! correct basic rounding error
        where (abs(nint(A)-A).lt.tiny) A = nint(A)
        ! At this point, A = [ LHS | RHS ], in which LHS = E, identity matrix; A completely specifies all relationships between variables: LHS = RHS.
        allocate(LHS(nbases,nbases))
        allocate(RHS(nbases,nbases))
        LHS = A(0*nbases+indices,0*nbases+indices)
        RHS = A(0*nbases+indices,1*nbases+indices)
        !
        ! checks
        if (.not.isequal(LHS,eye(nbases))) then
            stop 'Failed to reduce matrix to row echlon form.'
        endif
        if (count(get_null(RHS))+count(get_independent(RHS))+count(get_depenent(RHS)).ne.nbases) then
            stop 'Number of null, independent, and dependent terms do not sum to the number of terms.'
        endif
        !
        allocate(relations, source=RHS)
        !
    end function   combine_relations

    pure function  get_null(relations) result(is_null)
        !
        implicit none
        !
        real(dp), intent(in) :: relations(:,:)
        logical , allocatable :: is_null(:)
        !
        ! null terms (equal zero)
        allocate(is_null, source=all(abs(relations).lt.tiny,2))
        !
    end function   get_null

    pure function  get_independent(relations) result(is_independent)
        !
        implicit none
        !
        real(dp), intent(in) :: relations(:,:)
        logical , allocatable :: is_independent(:)
        !
        ! independent terms (equal themselves and nothing else)
        allocate(is_independent, source=all(abs(relations-eye(size(relations,1))).lt.tiny,2))
        !
    end function   get_independent

    pure function  get_depenent(relations) result(is_dependent)
        !
        implicit none
        !
        real(dp), intent(in) :: relations(:,:)
        logical , allocatable :: is_dependent(:)
        !
        ! dependent terms (can be written via independent terms)
        allocate(is_dependent, source=any(abs(relations).gt.tiny,2))
        is_dependent = (is_dependent.and..not.get_independent(relations))
        !
    end function   get_depenent

    subroutine     print_relations(relations,dims,flags)
        !
        ! flags - null, dependent, independent, header
        !       - bra, ket
        !
        implicit none
        !
        real(dp), intent(in) :: relations(:,:)
        integer , optional, intent(in) :: dims(:)
        character(*), intent(in) :: flags
        integer :: i, j 
        logical, allocatable :: is_null(:)
        logical, allocatable :: is_independent(:)
        logical, allocatable :: is_dependent(:)
        character(:), allocatable :: str
        integer :: nterms
        !
        ! get terms
        nterms = size(relations,1)
        is_null = get_null(relations)
        is_independent = get_independent(relations)
        is_dependent = get_depenent(relations)
        !
        ! print header
        if (index(flags,'header').ne.0) then
            write(*,'(a5,a,a)',advance='no') ' ... ', trim(int2char(nterms)), ' terms = '
            write(*,'(i4,a,f5.1,a)',advance='no') count(is_null)       , ' null ('       , count(is_null)       /real(nterms,dp)*100.0_dp , '%) '
            write(*,'(i4,a,f5.1,a)',advance='no') count(is_dependent)  , ' dependent ('  , count(is_dependent)  /real(nterms,dp)*100.0_dp , '%) '
            write(*,'(i4,a,f5.1,a)',advance='no') count(is_independent), ' independent (', count(is_independent)/real(nterms,dp)*100.0_dp , '%) '
            writE(*,*)
        endif
        !
        ! print irreducible symmetry relations
        if ((index(flags,'independent').ne.0).or.(index(flags,'dependent').ne.0).or.(index(flags,'null').ne.0)) then
            !
            if (.not.present(dims)) stop 'dims is required for printing relations'
            !
            ! write the independent terms (equal only to themselves)
            if (index(flags,'independent').ne.0) then
            do i = 1, nterms
                if (is_independent(i)) then
                    str = get_basis_label(dims=dims, ind=i ,flags=flags)
                    write(*,'(5x,a,a,a)',advance='no') trim(str),' = ', trim(str)
                    write(*,*)
                endif
            enddo
            endif
            !
            ! write the dependent terms
            if (index(flags,'dependent').ne.0) then
            do i = 1,nterms
                if (is_dependent(i)) then
                    str = get_basis_label(dims=dims, ind=i ,flags=flags)
                    write(*,'(5x,a,a)',advance='no') trim(str), ' = '
                    do j = 1,nterms
                        if (abs(relations(i,j)).gt.tiny) then
                            str = get_basis_label(dims=dims, ind=j ,flags=flags)
                            write(*,'(a,a,a)',advance='no') trim(dbl2charSP(relations(i,j),7)), '*', trim(str)
                        endif
                    enddo
                    write(*,*)
                endif
            enddo
            endif
            !
            ! write null terms
            if (index(flags,'null').ne.0) then
            do i = 1, nterms
                if (is_null(i)) then
                    str = get_basis_label(dims=dims, ind=i ,flags=flags)
                    write(*,'(5x,a,a)',advance='no') trim(str),' = 0'
                    write(*,*)
                endif
            enddo
            endif
        endif
        !
        contains
            function     get_basis_label(dims,ind,flags) result(str)
                !
                implicit none
                !
                integer     , intent(in) :: dims(:)
                integer     , intent(in) :: ind
                character(*), intent(in) :: flags
                character(:), allocatable :: str
                integer     , allocatable :: sub(:)
                integer :: i, n
                !
                allocate(character(100) :: str)
                !
                n = size(dims)
                !
                sub = ind2sub(dims=dims,ind=ind)
                !
                i = 1
                if     (index(flags,'bra').ne.0) then
                    str = '<'//trim(int2char(sub(i)))
                elseif (index(flags,'ket').ne.0) then
                    str = '|'//trim(int2char(sub(i)))
                else
                    str = 'a('//trim(int2char(sub(i)))
                endif
                !
                if (n.ge.2) then
                do i = 2, n
                    str = trim(str)//','//trim(int2char(sub(i)))
                enddo
                endif
                !
                if     (index(flags,'bra').ne.0) then
                    str = trim(str)//'|'
                elseif (index(flags,'ket').ne.0) then
                    str = trim(str)//'>'
                else
                    str = trim(str)//')'
                endif
                !
            end function get_basis_label
    end subroutine print_relations

    subroutine     dump_relations(relations,iopt_filename)
        !
        implicit none
        !
        real(dp), intent(in) :: relations(:,:)
        character(*), intent(in), optional :: iopt_filename
        character(100) :: fname
        integer :: fid
        integer :: m,n
        integer :: j,k
        !
        ! set default
        fname = 'dump.relations'
        if (present(iopt_filename)) fname = iopt_filename
        !
        m = size(relations,1)
        n = size(relations,2)
        !
        fid = 1
        open(unit=fid,file=trim(fname),status="replace",action='write')
            !
            do j = 1, m
            do k = 1, n
                write(fid,"(f)",advance='no') relations(j,k)
            enddo
            write(fid,*)
            enddo
            !
        close(fid)
        !
    end subroutine dump_relations

    ! write point symmetry operation in tight binding basis ()

    ! produces rotation which operates on a subsection of the Hamiltonian (useful for determining symmetry relations)

    function       ps2tb(R,atom) result(H)
        ! construct rotation representation which transforms all orbitals on input atom
        implicit none
        !
        real(dp), intent(in) :: R(3,3)
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
            H(H_start(i):H_end(i), H_start(i):H_end(i)) = rot2irrep(l=atom%azimuthal(i), R=R)
        enddo
        !
    end function   ps2tb

    ! produces rotation which commutes with the entire Hamiltonian (will be useful later when going to kpoints)

    function       ps2tb_H(R,pc,ic) result(H)
        !
        implicit none
        !
        real(dp), intent(in) :: R(3,3)
        type(am_class_prim_cell), intent(in) :: pc ! primitive cell
        type(am_class_irre_cell), intent(in) :: ic ! irreducible cell
        real(dp), allocatable :: H(:,:)
        integer , allocatable :: H_start(:), H_end(:)
        integer :: Hdim ! hamiltonian dimensions
        integer :: maxnazimuthals
        integer :: i,j,k
        !
        ! get largest azimuthal quantum number
        maxnazimuthals = 0
        do i = 1, ic%natoms
            if (maxnazimuthals.lt.ic%atom(i)%nazimuthals) maxnazimuthals = ic%atom(i)%nazimuthals
        enddo
        ! allocate space for defining subsections of the rotation in Hamiltonian basis
        i = maxnazimuthals * pc%natoms
        allocate(H_start(i))
        allocate(H_end(i))
        ! determine subsections: 1 per primitive atom per azimuthal quantum number per magnetic number per spin
        ! ignoring spin contributions at this stage
        Hdim = 0
        k = 0
        do i = 1, pc%natoms
        do j = 1, ic%atom(pc%ic_id(i))%nazimuthals
            k=k+1
            H_start(k) = Hdim + 1
            Hdim = Hdim + ic%atom(pc%ic_id(i))%azimuthal(j)*2+1
            H_end(k) = Hdim
        enddo
        enddo
        ! allocate and initialize space for rotations in Hamiltonin bais
        allocate(H(Hdim,Hdim))
        H = 0.0_dp
        ! construct rotation in the Hamiltonian basis
        k=0
        do i = 1, pc%natoms
        do j = 1, ic%atom(pc%ic_id(i))%nazimuthals
            k=k+1
            H(H_start(k):H_end(k), H_start(k):H_end(k)) = rot2irrep(l=ic%atom(pc%ic_id(i))%azimuthal(j), R=R)
        enddo
        enddo
        !
    end function   ps2tb_H

end module am_symmetry_rep









