module am_symmetry_rep

    use am_constants
    use am_stdout
    use am_mkl
    use am_options
    use am_symmetry
    use am_matlab

    implicit none

    private

    type, public, extends(am_class_group) :: am_class_flattened_group
        contains
        procedure :: get_intrinsic_symmetry_group
    end type am_class_flattened_group

    type, public, extends(am_class_group) :: am_class_regular_group
	    contains
        procedure :: create => create_regular
    end type am_class_regular_group

    type, public :: am_class_property
        character(50) :: property
        integer :: tensor_rank
        integer :: axial_polar
        type(am_class_flattened_group) :: ig ! intrinsic symmetries in flattened basis
        type(am_class_flattened_group) :: pg ! point symmetries in flattened basis
        type(am_class_flattened_group) :: fg ! all symmetries in flattened basis (direct product: ig * pg)
        contains
        procedure :: initialize_property
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
		rr%ndims = size(rr%sym,1)
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
            multab = get_multiplication_table(seitz=seitz,flags='seitz,sort')
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

    subroutine     get_property(prop,property,opts)
        !
        implicit none
        !
        class(am_class_property), intent(out) :: prop
        character(*)            , intent(in)  :: property
        type(am_class_options)  , intent(in)  :: opts
        !
        if (opts%verbosity.ge.1) call am_print_title('Determining group of '//trim(property))
        !
        call prop%initialize_property(property=property)
        !
        call ig%get_intrinsic_symmetry_group(ndims=ndims, nbases=nbases, property=property)
        !
    end subroutine get_property

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

    subroutine     get_flattened_group(grp,pg,uc,opts,property)
        !
        class(am_class_flattened_group), intent(out):: grp
        type(am_class_point_group)     , intent(in) :: pg
        type(am_class_unit_cell)       , intent(in) :: uc
        type(am_class_options)         , intent(in) :: opts
        character(*)                   , intent(in) :: property
        type(am_class_flattened_group) :: ig ! intrinsic symmetry group
        integer  :: ndims       ! number of spatial dimensions (ndims = 3, usually)
        integer  :: nbases      ! number of basis functions
        integer  :: tensor_rank ! tensor_rank rank of tensor matrix M, determined automatically by code
        real(dp), allocatable :: wrk(:,:,:) ! workspace for symmetries
        real(dp),allocatable :: sort_parameter(:)
        integer :: max_syms
        character(10) :: flags
        integer :: i, j
        !
        !
        ! number of spatial dimensions
        ndims = 3
        ! number of bases functions in representation
        nbases = ndims**tensor_rank
        !
        ! get intrinsic symmety group
        call ig%get_intrinsic_symmetry_group(ndims=ndims, nbases=nbases, property=property)
        !
        ! get point symmetries in the flattened basis (convert to cartesian coordinates in the process)
        allocate(wrk(nbases,nbases,pg%nsyms))
        do i = 1, pg%nsyms
            ! Nye, J.F. "Physical properties of crystals: their representation by tensors and matrices". p 133 Eq 7
            if     (index(flags,'axial').ne.0) then
                wrk(:,:,i) = kron_pow(ps_frac2cart(R_frac=pg%sym(1:3,1:3,i),bas=uc%bas),tensor_rank) * det(pg%sym(1:3,1:3,i))
            elseif (index(flags,'polar').ne.0) then
                wrk(:,:,i) = kron_pow(ps_frac2cart(R_frac=pg%sym(1:3,1:3,i),bas=uc%bas),tensor_rank)
            endif
        enddo
        !
        ! merge point symmetries with intrinsic symmetries
        grp%sym = direct_product(symA=wrk(:,:,1:pg%nsyms),symB=igs%sym(:,:,1:S_nsyms))
        ! get number of symmetries
        grp%nsyms = size(grp%sym,3)
        ! get size of bases
        grp%nbases = nbases
        ! get multiplication table
        grp%multab = get_multiplication_table(seitz=grp%sym, flags='')
        ! transfer ps_id
        write(*,*) 'Need a more inteligent way to do ps_id here.'
        allocate(grp%ps_id(grp%nsyms))
        grp%ps_id = 100
        do j = 1, grp%nsyms
        search : do i = 1, pg%nsyms
            if (isequal(grp%sym(:,:,j),wrk(:,:,i))) then
                grp%ps_id(j) = pg%ps_id(i)
                exit search
            endif
        enddo search
        enddo
        ! determine conjugacy classes (needs identity first, inversion second)
        grp%class_id = get_conjugacy_classes(multab=grp%multab, ps_id=grp%ps_id)
        !
        ! sort symmetries based on parameters
        allocate(sort_parameter(grp%nsyms))
        ! 
        do i = 1, grp%nsyms; sort_parameter(i) = det(grp%sym(:,:,i)); enddo
        call sort_symmetries(sg=grp, criterion=sort_parameter, flags='acsend')
        do i = 1, grp%nsyms; sort_parameter(i) = trace(grp%sym(:,:,i)); enddo
        call sort_symmetries(sg=grp, criterion=sort_parameter, flags='acsend')
        call sort_symmetries(sg=sg,  criterion=real(sg%class_id,dp), flags='ascend')
        !
        ! get character table
        sg%chartab = get_character_table(multab=sg%multab,ps_id=sg%ps_id)
        !
        contains
        function       direct_product(symA,symB) result(symC)
            !
            implicit none
            !
            real(dp), intent(in) :: symA(:,:,:)
            real(dp), intent(in) :: symB(:,:,:)
            real(dp), allocatable:: symC(:,:,:)
            real(dp), allocatable:: try(:,:)
            integer :: nsymsA
            integer :: nsymsB
            integer :: nsymsC
            integer :: i,j,k
            !
            if (size(symsA,1)/=size(symsB,1)) stop 'row mismatch'
            if (size(symsA,2)/=size(symsB,2)) stop 'column mismatch'
            !
            nbases = size(symsA,1)
            !
            nsymsA = size(symsA,3)
            nsymsB = size(symsB,3)
            nsymsC = nsymsA*nsymsB
            !
            allocate(try(nbases,nbases))
            !
            allocate(symC(nbases,nbases,nsymsC))
            symC = 0
            !
            k = 0
            do i = 1, nsymsA
            do j = 1, nsymsB
                try = matmul(symA(:,:,i),symB(:,:,j))
                if (.not.issubset(A=symC,B=try)) then
                    k = k + 1
                    symC(:,:,k) = try
                endif
            end do
            end do
            !
        end function   direct_product
        subroutine     multiply_to_closure(sym)
            !
            implicit none
            !
            real(dp), intent(inout) :: sym(:,:,:)
            real(dp), intent(in) :: closed_grp(:,:,:)
            integer :: max_iterations
            integer :: i
            !
            max_iterations = 10
            !
            do i = 1, max_iterations
                closed_grp = direct_product(symA=sym,symB=sym)
                if (size(sym,3).eq.size(closed_grp,3)) then
                    exit
                else
                    deallocate(sym)
                    allocate(sym,source=closed_grp)
                    if (i.eq.max_iterations) stop 'Too many intrisinc symmetries.'
                endif
            end do
            !
            !
            end subroutine multiply_to_closure
    end subroutine get_flattened_group

    subroutine     get_intrinsic_symmetry_group(ig,ndims,nbases,property)
        !
        implicit none
        !
        class(am_class_flattened_group), intent(out) :: ig ! intrinsic symmetry group
        integer     , intent(in) :: ndims
        integer     , intent(in) :: nbases
        character(*), intent(in) :: property
        real(dp), allocatable :: T(:,:)   ! transpositional operator building block
        integer :: k
        integer :: max_syms
        !
        ! generate intrinsic symmetry group
        T = transp_operator(ndims)
        !
        ! allocate a little space just for generators below
        ! eventually this space will expand in multiply_to_closure
        max_syms = 100
        allocate(ig%sym(nbases,nbases,max_syms))
        !
        k=0
        if     (index(property,'conductivity'    ).ne.0 &
         & .or. index(property,'resistivity'     ).ne.0 &
         & .or. index(property,'voigt'           ).ne.0) then
            k=k+1; ig%sym(:,:,k) = T                      ! s_ij  = s_ji
        elseif (index(property,'piezoelectricity').ne.0) then
            k=k+1; ig%sym(:,:,k) = kron(T,eye(ndims))     ! d_ijk = d_ikj
        elseif (index(property,'elasticity'      ).ne.0) then
            k=k+1; ig%sym(:,:,k) = eye(nbases)            ! cijkl = cijkl
            k=k+1; ig%sym(:,:,k) = kron(eye(ndims**2),T)  ! cijkl = cjikl
            k=k+1; ig%sym(:,:,k) = kron(T,eye(ndims**2))  ! cijkl = cjilk
            k=k+1; ig%sym(:,:,k) = kron(T,T)              ! cijkl = cjilk
        endif
        ! multiply to make ensure group closure
        call multiply_to_closure(sym=ig%sym)
        !
        ! get number of intrisic symmetries and the size of the basis
        ig%nsyms = size(ig%sym,3)
        ig%nbases= nbases
        !
        ! label symmetries by id
        allocate(ig%ps_id(ig%nsyms)))
        ig%ps_id = 10
        do i = 1, ig%nsyms
            if (isequal(ig%sym(:,:,i),eye(ig%nbases))) then
                ig%ps_id(i) = 1
                exit
            endif
        enddo
        !
    end function   get_intrinsic_symmetry_group

end module am_symmetry_rep









