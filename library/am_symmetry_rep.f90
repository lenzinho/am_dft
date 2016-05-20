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
        procedure :: get_flattened_intrinsic_symmetry_group
        procedure :: get_flattened_point_group
        procedure :: get_direct_product
        procedure :: get_symmetry_relations
    end type am_class_flattened_group

    type, public, extends(am_class_group) :: am_class_regular_group
	    contains
        procedure :: create => create_regular
    end type am_class_regular_group

    type, public :: am_class_property
        character(50) :: property
        character(5)  :: axial_polar
        integer       :: tensor_rank
        type(am_class_flattened_group) :: fig  ! flattened intrinsic symmetry group
        type(am_class_flattened_group) :: fpg  ! flattened point group
        type(am_class_flattened_group) :: flat ! flattened group (direct product: fig * pg)
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

    !

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
        call prop%fig%get_flattened_intrinsic_symmetry_group(prop=prop)
        if (opts%verbosity.ge.1) call am_print('intrinsic symmetries', prop%fig%nsyms)
        !
        call prop%fpg%get_flattened_point_group(prop=prop, pg=pg, uc=uc)
        if (opts%verbosity.ge.1) call am_print('point symmetries', prop%fpg%nsyms)
        !
        call prop%flat%get_direct_product(A=prop%fig,B=prop%fpg)
        if (opts%verbosity.ge.1) call am_print('total symmetries (direct product: intrinsic * point)', prop%flat%nsyms)
        !
        if (opts%verbosity.ge.1) write(*,'(a5,a)') ' ... ', 'getting symmetry relations'
        call prop%flat%get_symmetry_relations()
        !
        call print_relations(relations=prop%flat%relations)
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
            if     (index(prop%property,'pyroelectricity')          .ne.0) then; prop%tensor_rank = 1                             !    ! P_{i}     = p_{i} \Delta T
            !------------------------------------------------- SECOND-RANK TENSORS -------------------------------------------------
            elseif (index(prop%property,'electrical susceptibility').ne.0) then; prop%tensor_rank = 2; prop%axial_polar = 'polar' ! S  ! P_{i}     = \alpha_{ij}  E_{j}
            elseif (index(prop%property,'magnetic susceptibility')  .ne.0) then; prop%tensor_rank = 2; prop%axial_polar = 'axial' ! S  ! M_{i}     = \mu_{ij}     H_{j}
            elseif (index(prop%property,'magneto-electric')         .ne.0) then; prop%tensor_rank = 2; prop%axial_polar = 'axial' 
            elseif (index(prop%property,'thermal expansion')        .ne.0) then; prop%tensor_rank = 2; prop%axial_polar = 'polar' ! S  ! \eps_{ij} = \alpha_{ij}  \Delta T
            elseif (index(prop%property,'electrical conductivity')  .ne.0) then; prop%tensor_rank = 2; prop%axial_polar = 'polar' ! S  ! J_{i}     = \sigma_{ij}  E_{i}
            elseif (index(prop%property,'electrical resistivity')   .ne.0) then; prop%tensor_rank = 2; prop%axial_polar = 'polar' ! S  ! E_{i}     = \rho_{ij}    J_{j}
            elseif (index(prop%property,'thermal conductivity')     .ne.0) then; prop%tensor_rank = 2; prop%axial_polar = 'polar' ! S  ! q_{i}     = \kappa_{ij}  \frac{\partial T}/{\partial r_{j}}
            elseif (index(prop%property,'thermoelectricity')        .ne.0) then; prop%tensor_rank = 2; prop%axial_polar = 'polar' ! N  ! 
            elseif (index(prop%property,'seebeck')                  .ne.0) then; prop%tensor_rank = 2; prop%axial_polar = 'polar' ! N  ! E_{i}     = \beta_{ij}   \frac{\partial T}/{\partial r_{j}}
            elseif (index(prop%property,'peltier')                  .ne.0) then; prop%tensor_rank = 2; prop%axial_polar = 'polar' ! N  ! q_{i}     = \pi_{ij}     J_{j}
            !------------------------------------------------- THIRD-RANK TENSORS -------------------------------------------------
            elseif (index(prop%property,'hall')                     .ne.0) then; prop%tensor_rank = 3;                            !    ! E_{i}     = h_{ijk}      J_{j} H_{k} 
            elseif (index(prop%property,'piezoelectricity')         .ne.0) then; prop%tensor_rank = 3; prop%axial_polar = 'polar' !    ! P_{i}     = d_{ijk}      \sigma_{jk}
            elseif (index(prop%property,'piezomagnetic')            .ne.0) then; prop%tensor_rank = 3; prop%axial_polar = 'axial' !    ! M_{i}     = Q_{ijk}      \sigma_{jk}
            !------------------------------------------------- FOURTH-RANK TENSORS ------------------------------------------------
            elseif (index(prop%property,'elasticity')               .ne.0) then; prop%tensor_rank = 4; prop%axial_polar = 'polar' !    ! 
            elseif (index(prop%property,'piezo-optic')              .ne.0) then; prop%tensor_rank = 4                             !    ! 
            elseif (index(prop%property,'kerr')                     .ne.0) then; prop%tensor_rank = 4                             !    ! 
            elseif (index(prop%property,'electrostriction')         .ne.0) then; prop%tensor_rank = 4                             !    ! 
            !------------------------------------------------- SXITH-RANK TENSORS -------------------------------------------------
            elseif (index(prop%property,'third-order elasticity')   .ne.0) then; prop%tensor_rank = 6; prop%axial_polar = 'polar' !    ! 
            !----------------------------------------------------------------------------------------------------------------------
            else
                stop 'Unknown property.'
            endif
            !----------------------------------------------------------------------------------------------------------------------
        end subroutine initialize_property
    end subroutine get_property

    subroutine     get_flattened_intrinsic_symmetry_group(fig,prop)
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
        !
        ! get transpositional operator
        T = transp_operator(ndims)
        !
        ! number of bases functions in representation
        fig%nbases = ndims**prop%tensor_rank
        !
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
    end subroutine get_flattened_intrinsic_symmetry_group

    subroutine     get_flattened_point_group(fpg,prop,pg,uc)
        !
        class(am_class_flattened_group), intent(out):: fpg ! flattened point group
        class(am_class_property),        intent(in) :: prop! properties
        type(am_class_point_group)     , intent(in) :: pg  ! seitz point group
        type(am_class_unit_cell)       , intent(in) :: uc  ! unit cell
        integer :: ndims       ! number of spatial dimensions (ndims = 3, usually)
        integer :: i
        !
        !
        ! number of spatial dimensions
        ndims = 3
        ! number of bases functions in representation
        fpg%nbases = ndims**prop%tensor_rank
        ! number of symmetries
        fpg%nsyms = pg%nsyms
        ! generate intrinsic symmetries
        allocate(fpg%sym(fpg%nbases,fpg%nbases,fpg%nsyms))
        fpg%sym = 0
        ! Nye, J.F. "Physical properties of crystals: their representation by tensors and matrices". p 133 Eq 7
        do i = 1, pg%nsyms
            if     (index(prop%axial_polar,'axial').ne.0) then
                fpg%sym(:,:,i) = kron_pow(ps_frac2cart(R_frac=pg%sym(1:3,1:3,i),bas=uc%bas),prop%tensor_rank) * det(pg%sym(1:3,1:3,i))
            elseif (index(prop%axial_polar,'polar').ne.0) then
                fpg%sym(:,:,i) = kron_pow(ps_frac2cart(R_frac=pg%sym(1:3,1:3,i),bas=uc%bas),prop%tensor_rank)
            endif
            ! correct basic rounding error
            where (abs(fpg%sym(:,:,i)).lt.tiny) fpg%sym(:,:,i) = 0
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
    end subroutine get_flattened_point_group

    subroutine     get_direct_product(C,A,B)
        !
        use am_progress_bar
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
            if(mod(k,10).eq.0) call progress_bar(iteration=k, maximum=A%nsyms*B%nsyms)
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

    subroutine     get_symmetry_relations(flat)
        ! At this point A is an augmented matrix of the form
        !
        !        [ A1 , B1 ]
        !   A =  [ A2 , B2 ]
        !        [ A3 , I  ]
        !
        ! in which the square matrix [A1,B1;A2,B2] is an augmented upper triangular matrix corresponding to 2*flat%nbases
        ! coupled equations which describe how the terms are interrelated. The last augmented matrix [A3,I]
        ! describes how the tensor rotates under the most recently considered symmetry operation R. In order
        ! incorporate the effect of symmetry operation R on the set of linearly coupled 2*flat%nbases equations, LU
        ! factorization is performed. By doing so, the augmented matrix A is converted into its factored upper
        ! triangular matrix U, which, by definition, only occupies rows corresponding to the smallest dimension of
        ! the matrix [either row or column, i.e. min(m,n)]. For this case, it occupies the top 2x2 augmented blocks.
        ! This is also why the symmetry equations are added as components in the [A3,I] block -- so they don't
        ! overlap with U.
        !
        implicit none
        !
        class(am_class_flattened_group), intent(inout) :: flat
        integer , allocatable :: class_member(:,:)  !
        integer , allocatable :: class_nelements(:) ! 
        real(dp), allocatable :: A(:,:)             ! A(2*flat%nbases,2*flat%nbases) - augmented matrix equation
        real(dp), allocatable :: LHS(:,:)           ! LHS(flat%nbases,flat%nbases)   - left hand side of augmented matrix equation (should be identity after reducing to row echlon form)
        real(dp), allocatable :: RHS(:,:)           ! RHS(flat%nbases,flat%nbases)   - right hand side of augmented matrix equation
        integer , allocatable :: indices(:)         ! used for clarity
        integer :: nequations                       ! equation counter, for fun
        integer :: i, j                             ! loop variables
        !
        !
        ! get conjugacy class members
        class_member = member(id=flat%class_id)
        ! get nmber of elements in each class
        class_nelements = nelements(id=flat%class_id)
        ! intialize indices
        allocate(indices,source=[1:flat%nbases])
        ! initialize equation counter
        nequations = 0
        ! initialize augmented workspace matrix A
        allocate(A(3*flat%nbases,2*flat%nbases))
        A = 0
        ! use LU factorization to incorporate, one symmetry at a time, the effect of all symmetries on the flattened basis
        do j = 1, size(class_member,1) ! loop over classes
            ! track number of symmetry equations for fun
            nequations = nequations + flat%nbases*class_nelements(j)
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
        if (count(get_zeros(RHS))+count(get_indep(RHS))+count(get_depen(RHS)).ne.flat%nbases) then
            stop 'Number of null, independent, and dependent terms do not sum to the number of terms.'
        endif
        !
        allocate(flat%relations, source=RHS)
        !
    end subroutine get_symmetry_relations

    function       get_zeros(relations) result(is_zero)
        !
        implicit none
        !
        real(dp), intent(in) :: relations(:,:)
        logical , allocatable :: is_zero(:)
        !
        ! null terms (equal zero)
        allocate(is_zero, source=all(abs(relations).lt.tiny,2))
        !
    end function   get_zeros

    function       get_indep(relations) result(is_independent)
        !
        implicit none
        !
        real(dp), intent(in) :: relations(:,:)
        logical , allocatable :: is_independent(:)
        !
        ! independent terms (equal themselves and nothing else)
        allocate(is_independent, source=all(abs(relations-eye(size(relations,1))).lt.tiny,2))
        !
    end function   get_indep

    function       get_depen(relations) result(is_dependent)
        !
        implicit none
        !
        real(dp), intent(in) :: relations(:,:)
        logical , allocatable :: is_dependent(:)
        !
        ! dependent terms (can be written via independent terms)
        allocate(is_dependent, source=any(abs(relations).gt.tiny,2))
        is_dependent = (is_dependent.and..not.get_indep(relations))
        !
    end function   get_depen

    subroutine     print_relations(relations)
        !
        implicit none
        !
        real(dp), intent(in) :: relations(:,:)
        integer :: i, j 
        logical, allocatable :: is_zero(:)
        logical, allocatable :: is_independent(:)
        logical, allocatable :: is_dependent(:)
        integer :: nterms
        !
        ! get terms
        nterms = size(relations,1)
        ! null terms (equal zero)
        is_zero = get_zeros(relations)
        call am_print('null terms',count(is_zero),' ... ')
        ! independent terms (equal themselves and nothing else)
        is_independent = get_indep(relations)
        call am_print('independent terms',count(is_independent),' ... ')
        ! dependent terms (can be written via independent terms)
        is_dependent = get_depen(relations)
        call am_print('number of dependent terms',count(is_dependent),' ... ')
        ! irreducible symmetry relations
        write(*,'(a5,a)') ' ... ', 'irreducible symmetry relations (null terms omitted)'
        ! write the independent terms (equal only to themselves)
        do i = 1, nterms
            if (is_independent(i)) then
                write(*,'(5x,a1,a,a1,a1,a)',advance='no') 'a',trim(int2char(i)),'=', 'a', trim(int2char(i))
                write(*,*)
            endif
        enddo
        ! write the dependent terms
        do i = 1,nterms
            if (is_dependent(i)) then
                write(*,'(5x,a1,a,a1)',advance='no') 'a',trim(int2char(i)),'='
                do j = 1,nterms
                    if (abs(relations(i,j)).gt.tiny) then
                        write(*,'(a,a,a)',advance='no') trim(dbl2charSP(relations(i,j),7)), '*a', trim(int2char(j))
                    endif
                enddo
                write(*,*)
            endif
        enddo
        !
    end subroutine print_relations


end module am_symmetry_rep









