module am_symmetry_rep

    use am_constants
    use am_stdout
    use am_mkl
    use am_options
    use am_symmetry
    use am_matlab
    use am_unit_cell

    implicit none

    private

    type, public, extends(am_class_group) :: am_class_flattened_group
        real(dp), allocatable :: relations(:,:)
        contains
        procedure :: get_flat_intrinsic_group
        procedure :: get_flat_point_group
        procedure :: get_direct_product
        procedure :: get_relations
    end type am_class_flattened_group

    type, public, extends(am_class_group) :: am_class_regular_group
	    contains
        procedure :: create => create_regular
    end type am_class_regular_group

    type, public :: am_class_property
        character(50) :: property
        character(5)  :: axial_polar
        integer       :: rank ! tensor rank
        integer , allocatable :: dims(:) ! tensor dimensions
        real(dp), allocatable :: relations(:,:)
        ! direct product of intrinsic and point groups takes too long (merging relations instead)
        type(am_class_flattened_group) :: fig  ! flattened intrinsic symmetry group
        type(am_class_flattened_group) :: fpg  ! flattened point group
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
            multab = get_multiplication_table(sym=seitz,flags='seitz,sort')
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

    ! operates on prop

    subroutine     get_property(prop,uc,pg,opts,property)
        !
        implicit none
        !
        class(am_class_property),    intent(out):: prop
        class(am_class_unit_cell),   intent(in) :: uc
        class(am_class_point_group), intent(in) :: pg
        type(am_class_options),      intent(in) :: opts
        character(*),                intent(in) :: property
        !
        if (opts%verbosity.ge.1) call am_print_title('Determining group of '//trim(property))
        !
        call initialize_property(prop=prop, property=property)
        !
        call prop%fig%get_flat_intrinsic_group(prop=prop)
        if (opts%verbosity.ge.1) call am_print('intrinsict symmetries', prop%fig%nsyms)
        !
        prop%fig%relations = prop%fig%get_relations()
        call print_relations(relations=prop%fig%relations,flags='')
        !
        call prop%fpg%get_flat_point_group(prop=prop, pg=pg, uc=uc)
        if (opts%verbosity.ge.1) call am_print('point symmetries', prop%fpg%nsyms)
        !
        prop%fpg%relations = prop%fpg%get_relations()
        call print_relations(relations=prop%fpg%relations,flags='')
        !
        if (opts%verbosity.ge.1) call am_print('total symmetries', prop%fpg%nsyms * prop%fig%nsyms )
        !
        prop%relations = combine_relations(prop%fig%relations, prop%fpg%relations)
        call print_relations(relations=prop%relations, dims=prop%dims, flags='show:dependent,independent')
        !
        ! call prop%flat%get_direct_product(A=prop%fig,B=prop%fpg)
        ! !
        ! if (opts%verbosity.ge.1) call am_print('total symmetries (direct product: intrinsic * point)', prop%flat%nsyms)
        ! !
        ! if (opts%verbosity.ge.1) write(*,'(a5,a)') ' ... ', 'getting symmetry relations'
        ! !
        ! call prop%flat%get_relations()
        ! !
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
            ! POLAR :  T_{i,j,k,l,...} =        sum_{ip,jp,kp,lp,...} R_{i,ip} R_{j,jp} R_{k,kp} R_{l,lp} T_{ip,jp,kp,lp,...}
            ! AXIAL :  T_{i,j,k,l,...} = det(R) sum_{ip,jp,kp,lp,...} R_{i,ip} R_{j,jp} R_{k,kp} R_{l,lp} T_{ip,jp,kp,lp,...}
            ! Thus, all axial tensors of even rank and polar tensors of odd rank are null are null.  Wooten p 485. Eq. 13.21. 
            !
            ! Onsager’s Principle requires that the electric resistivity and thermal conductivity tensors be symmetric.
            ! This does not hold for the Seebeck and Peltier (thermoelectric) tensors which relate two different flows. Thus
            ! there are, at most, nine independent parameters rather than six. [Newnham "Properties of Materials"]
            !
            prop%property = property
            !------------------------------------------------- FIRST-RANK TENSORS --------------------------------------------------
            if     (index(prop%property,'pyroelectricity')          .ne.0) then; prop%rank = 1                             !    ! P_{i}     = p_{i} \Delta T
            !------------------------------------------------- SECOND-RANK TENSORS -------------------------------------------------
            elseif (index(prop%property,'electrical susceptibility').ne.0) then; prop%rank = 2; prop%axial_polar = 'polar' ! S  ! P_{i}     = \alpha_{ij}  E_{j}
            elseif (index(prop%property,'magnetic susceptibility')  .ne.0) then; prop%rank = 2; prop%axial_polar = 'axial' ! S  ! M_{i}     = \mu_{ij}     H_{j}
            elseif (index(prop%property,'magneto-electric')         .ne.0) then; prop%rank = 2; prop%axial_polar = 'axial' 
            elseif (index(prop%property,'thermal expansion')        .ne.0) then; prop%rank = 2; prop%axial_polar = 'polar' ! S  ! \eps_{ij} = \alpha_{ij}  \Delta T
            elseif (index(prop%property,'electrical conductivity')  .ne.0) then; prop%rank = 2; prop%axial_polar = 'polar' ! S  ! J_{i}     = \sigma_{ij}  E_{i}
            elseif (index(prop%property,'electrical resistivity')   .ne.0) then; prop%rank = 2; prop%axial_polar = 'polar' ! S  ! E_{i}     = \rho_{ij}    J_{j}
            elseif (index(prop%property,'thermal conductivity')     .ne.0) then; prop%rank = 2; prop%axial_polar = 'polar' ! S  ! q_{i}     = \kappa_{ij}  \frac{\partial T}/{\partial r_{j}}
            elseif (index(prop%property,'thermoelectricity')        .ne.0) then; prop%rank = 2; prop%axial_polar = 'polar' ! N  ! 
            elseif (index(prop%property,'seebeck')                  .ne.0) then; prop%rank = 2; prop%axial_polar = 'polar' ! N  ! E_{i}     = \beta_{ij}   \frac{\partial T}/{\partial r_{j}}
            elseif (index(prop%property,'peltier')                  .ne.0) then; prop%rank = 2; prop%axial_polar = 'polar' ! N  ! q_{i}     = \pi_{ij}     J_{j}
            !------------------------------------------------- THIRD-RANK TENSORS -------------------------------------------------
            elseif (index(prop%property,'hall')                     .ne.0) then; prop%rank = 3;                            !    ! E_{i}     = h_{ijk}      J_{j} H_{k} 
            elseif (index(prop%property,'piezoelectricity')         .ne.0) then; prop%rank = 3; prop%axial_polar = 'polar' !    ! P_{i}     = d_{ijk}      \sigma_{jk}
            elseif (index(prop%property,'piezomagnetic')            .ne.0) then; prop%rank = 3; prop%axial_polar = 'axial' !    ! M_{i}     = Q_{ijk}      \sigma_{jk}
            !------------------------------------------------- FOURTH-RANK TENSORS ------------------------------------------------
            elseif (index(prop%property,'elasticity')               .ne.0) then; prop%rank = 4; prop%axial_polar = 'polar' !    ! 
            elseif (index(prop%property,'piezo-optic')              .ne.0) then; prop%rank = 4                             !    ! 
            elseif (index(prop%property,'kerr')                     .ne.0) then; prop%rank = 4                             !    ! 
            elseif (index(prop%property,'electrostriction')         .ne.0) then; prop%rank = 4                             !    ! 
            !------------------------------------------------- SXITH-RANK TENSORS -------------------------------------------------
            elseif (index(prop%property,'third-order elasticity')   .ne.0) then; prop%rank = 6; prop%axial_polar = 'polar' !    ! 
            !----------------------------------------------------------------------------------------------------------------------
            else
                stop 'Unknown property.'
            endif
            !----------------------------------------------------------------------------------------------------------------------
            allocate(prop%dims(prop%rank))
            prop%dims = 3
            !
        end subroutine initialize_property
    end subroutine get_property

    subroutine     get_flat_intrinsic_group(fig,prop)
        !
        implicit none
        !
        class(am_class_flattened_group), intent(out) :: fig ! intrinsic symmetry group
        class(am_class_property),        intent(in)  :: prop
        real(dp), allocatable :: T(:,:) ! transpositional operator building block
        integer :: ndims
        integer :: k
        !
        ! number of spatial dimensions
        ndims = 3
        ! get transpositional operator
        T = transp_operator(ndims)
        ! number of bases functions in representation
        fig%nbases = ndims**prop%rank
        ! generate intrinsic symmetries
        k=0
        if     (index(prop%property,'conductivity'     ).ne.0 &
         & .or. index(prop%property,'resistivity'      ).ne.0 &
         & .or. index(prop%property,'voigt'            ).ne.0) then
            !
            fig%nsyms = 2
            allocate(fig%sym(fig%nbases,fig%nbases,fig%nsyms))
            k=k+1; fig%sym(:,:,k) = eye(fig%nbases)        ! E
            k=k+1; fig%sym(:,:,k) = T                      ! s_ij = s_ji
            !
        elseif (index(prop%property,'piezoelectricity' ).ne.0) then
            !
            fig%nsyms = 2
            allocate(fig%sym(fig%nbases,fig%nbases,fig%nsyms))
            k=k+1; fig%sym(:,:,k) = eye(fig%nbases)        ! E
            k=k+1; fig%sym(:,:,k) = kron(T,eye(ndims))     ! d_ijk = d_ikj
            !
        elseif (index(prop%property,'elasticity'       ).ne.0) then
            !
            fig%nsyms = 4
            allocate(fig%sym(fig%nbases,fig%nbases,fig%nsyms))
            k=k+1; fig%sym(:,:,k) = eye(fig%nbases)        ! E
            k=k+1; fig%sym(:,:,k) = kron(eye(ndims**2),T)  ! cijkl = cjikl
            k=k+1; fig%sym(:,:,k) = kron(T,eye(ndims**2))  ! cijkl = cjilk
            k=k+1; fig%sym(:,:,k) = kron(T,T)              ! cijkl = cjilk
            !
        elseif (index(prop%property,'thermoelectricity').ne.0) then
            fig%nsyms = 1
            allocate(fig%sym(fig%nbases,fig%nbases,fig%nsyms))
            k=k+1; fig%sym(:,:,k) = eye(fig%nbases)        ! E
        else
            stop 'Undefined propery encountered.'
        endif
        !
        allocate(fig%ps_id(fig%nsyms))
        fig%ps_id    = default_ps_id_value
        fig%ps_id(1) = 1
        !
        ! get multiplication table
        fig%multab   = get_multiplication_table(sym=fig%sym,flags='')
        !
        ! determine conjugacy classes (needs identity first, inversion second)
        fig%class_id = get_conjugacy_classes(multab=fig%multab, ps_id=fig%ps_id)
        !
        ! sort symmetries based on parameters
        call fig%sort_symmetries(criterion=fig%sym(1,4,:)        , flags='ascend')
        call fig%sort_symmetries(criterion=fig%sym(2,4,:)        , flags='ascend')
        call fig%sort_symmetries(criterion=fig%sym(3,4,:)        , flags='ascend')
        call fig%sort_symmetries(criterion=real(fig%ps_id,dp)    , flags='acsend')
        call fig%sort_symmetries(criterion=real(fig%class_id,dp) , flags='ascend')
        !
        ! get character table
        fig%chartab = get_character_table(multab=fig%multab,ps_id=fig%ps_id)
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

    subroutine     get_flat_point_group(fpg,prop,pg,uc)
        !
        class(am_class_flattened_group), intent(out):: fpg ! flattened point group
        class(am_class_property),        intent(in) :: prop! properties
        type(am_class_point_group)     , intent(in) :: pg  ! seitz point group
        type(am_class_unit_cell)       , intent(in) :: uc  ! unit cell
        integer :: i
        !
        !
        ! number of bases functions in representation
        fpg%nbases = product(prop%dims)
        ! number of symmetries
        fpg%nsyms = pg%nsyms
        ! generate intrinsic symmetries
        allocate(fpg%sym(fpg%nbases,fpg%nbases,fpg%nsyms))
        fpg%sym = 0
        ! Nye, J.F. "Physical properties of crystals: their representation by tensors and matrices". p 133 Eq 7
        do i = 1, pg%nsyms
            !
            if     (index(prop%axial_polar,'axial').ne.0) then
                fpg%sym(:,:,i) = kron_pow(ps_frac2cart(R_frac=pg%sym(1:3,1:3,i),bas=uc%bas),prop%rank) * det(pg%sym(1:3,1:3,i))
            elseif (index(prop%axial_polar,'polar').ne.0) then
                fpg%sym(:,:,i) = kron_pow(ps_frac2cart(R_frac=pg%sym(1:3,1:3,i),bas=uc%bas),prop%rank)
            endif
            ! correct basic rounding error
            where (abs(fpg%sym(:,:,i)).lt.tiny)  fpg%sym(:,:,i) = 0
            where (abs(fpg%sym(:,:,i)-nint(fpg%sym(:,:,i))).lt.tiny) fpg%sym(:,:,i) = nint(fpg%sym(:,:,i))
        enddo
        !
        ! copy symmetry ids
        allocate(fpg%ps_id, source=pg%ps_id)
        ! copy classes
        allocate(fpg%class_id, source=pg%class_id)
        ! copy multiplication table
        allocate(fpg%multab, source=pg%multab)
        ! copy character table
        allocate(fpg%chartab, source=pg%chartab)
        !
        if (.not.isequal(fpg%sym(:,:,1),eye(fpg%nbases))) stop 'FPG: Identity is not first.'
        !
    end subroutine get_flat_point_group

    subroutine     get_direct_product(C,A,B)
        !
        implicit none
        !
        class(am_class_flattened_group) :: C
        type(am_class_flattened_group)  :: A
        type(am_class_flattened_group)  :: B
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
        ! get multiplication table
        C%multab = get_multiplication_table(sym=C%sym,flags='prog')
        ! determine conjugacy classes (needs identity first, inversion second)
        C%class_id = get_conjugacy_classes(multab=C%multab,ps_id=C%ps_id)
        !
        ! sort symmetries based on parameters
        call C%sort_symmetries(criterion=real(C%ps_id,dp)   , flags='acsend')
        call C%sort_symmetries(criterion=real(C%class_id,dp), flags='ascend')
        !
        ! get character table
        C%chartab = get_character_table(multab=C%multab, ps_id=C%ps_id)
        !
    end subroutine get_direct_product

    ! operates on relations

    function       get_relations(flat) result(relations)
        ! 
        implicit none
        !
        class(am_class_flattened_group), intent(inout) :: flat
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
        class_member = member(id=flat%class_id)
        ! intialize indices
        allocate(indices,source=[1:flat%nbases])
        ! initialize augmented workspace matrix A
        allocate(A(3*flat%nbases,2*flat%nbases))
        A = 0
        ! use LU factorization to incorporate, one symmetry at a time, the effect of all symmetries on the flattened basis
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

    function       combine_relations(relationsA,relationsB) result(relationsC)
        !
        implicit none
        !
        real(dp), intent(in)  :: relationsA(:,:)    !
        real(dp), intent(in)  :: relationsB(:,:)    !
        real(dp), allocatable :: relationsC(:,:)    !
        integer :: nbases
        real(dp), allocatable :: A(:,:)             ! A(2*flat%nbases,2*flat%nbases) - augmented matrix equation
        real(dp), allocatable :: LHS(:,:)           ! LHS(flat%nbases,flat%nbases)   - left hand side of augmented matrix equation (should be identity after reducing to row echlon form)
        real(dp), allocatable :: RHS(:,:)           ! RHS(flat%nbases,flat%nbases)   - right hand side of augmented matrix equation
        integer , allocatable :: indices(:)         ! used for clarity
        !
        !
        if (size(relationsA,1).ne.size(relationsB,1)) stop 'Dimension mismatch: A1 vs B1'
        if (size(relationsA,2).ne.size(relationsB,2)) stop 'Dimension mismatch: A2 vs B2'
        if (size(relationsA,1).ne.size(relationsA,2)) stop 'Dimension mismatch: A1 vs A2'
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
        ! incorporate symmetry via lu factorization (equivalent to applying rref)
        call lu(A)
        ! Apply Gram-Schmidt orthogonalization to obtain A in reduced row echelon form
        call rref(A)
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
        allocate(relationsC, source=RHS)
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
        ! flags - null, dependent, independent
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
        write(*,'(a5,a,a)',advance='no') ' ... ', trim(int2char(nterms)), ' terms = '
        write(*,'(i4,a,f5.1,a)',advance='no') count(is_null)       , ' null ('       , count(is_null)       /real(nterms,dp)*100.0_dp , '%) '
        write(*,'(i4,a,f5.1,a)',advance='no') count(is_independent), ' independent (', count(is_independent)/real(nterms,dp)*100.0_dp , '%) '
        write(*,'(i4,a,f5.1,a)',advance='no') count(is_dependent)  , ' dependent ('  , count(is_dependent)  /real(nterms,dp)*100.0_dp , '%) '
        writE(*,*)
        !
        ! irreducible symmetry relations
        if ((index(flags,'independent').ne.0).or.(index(flags,'dependent').ne.0).or.(index(flags,'null').ne.0)) then
            !
            if (.not.present(dims)) stop 'dims is required for priting relations'
            !
            write(*,'(a5,a)') ' ... ', 'irreducible symmetry relations'
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
                    str = get_basis_label(dims=dims, ind=i ,flags='')
                    write(*,'(5x,a,a)',advance='no') trim(str), ' = '
                    do j = 1,nterms
                        if (abs(relations(i,j)).gt.tiny) then
                            str = get_basis_label(dims=dims, ind=j ,flags=flags)
                            write(*,'(a,a)',advance='no') trim(dbl2charSP(relations(i,j),7)), trim(str)
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
                if     (index(flags,'bra')) then
                    str = '<'//trim(int2char(sub(i)))
                elseif (index(flags,'ket')) then
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
                if     (index(flags,'bra')) then
                    str = trim(str)//'|'
                elseif (index(flags,'ket')) then
                    str = trim(str)//'>'
                else
                    str = trim(str)//')'
                endif
                !
            end function get_basis_label
    end subroutine print_relations


end module am_symmetry_rep









