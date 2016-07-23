module am_symmetry

    use am_constants
    use am_stdout
    use am_unit_cell
    use am_options
    use am_mkl
    use am_symmetry_tables
    use am_matlab
    use dispmodule

    implicit none

    public

    type, public :: am_class_group
        integer :: nsyms                              ! number of symmetries
        integer :: nbases                             ! number of bases functions (dimensions of rep)
        real(dp), allocatable :: sym(:,:,:)           ! symmetry elements 
                                                      ! (for seitz groups: regular rep; for reppresentation groups: rep engendered by the basis)
        integer , allocatable :: ps_id(:)             ! integer which identifies point symmetries (see decode_pointsymmetry)
        type(am_class_conjugacy_class)      :: cc     ! conjugacy classes
        type(am_class_character_table)      :: ct     ! character table
        type(am_class_multiplication_table) :: mt     ! multiplication table
        integer , allocatable :: supergroup_id(:)     ! if this is a subgroup, supergroup_id identifis the symmetry in the encompassing group (used for stabilizer group)
        integer , allocatable :: ax_id(:)
        contains
        procedure :: sort_symmetries
        procedure :: get_multiplication_table
        procedure :: get_conjugacy_classes
        procedure :: get_character_table
        procedure :: print_character_table
        procedure :: debug_dump => debug_dump_grp
    end type am_class_group

    type, public, extends(am_class_group)       :: am_class_representation_group
        ! generic symmetry group
    end type am_class_representation_group

    type, public, extends(am_class_group)       :: am_class_seitz_group
        real(dp) :: bas(3,3)                          ! basis in which [frac.] are defined
        real(dp) :: recbas(3,3)                       ! basis in which [recfrac] are defined
        real(dp), allocatable :: seitz_cart(:,:,:)    ! symmetry elements [cart.]
        real(dp), allocatable :: seitz_frac(:,:,:)    ! symmetry elements [frac, direct]; acts on atoms
        real(dp), allocatable :: seitz_recfrac(:,:,:) ! symmetry elements [frac, reciprocal]; acts on kpoints
        contains
        procedure :: create_seitz_group
        procedure :: write_action_table
    end type am_class_seitz_group

    type, public, extends(am_class_seitz_group) :: am_class_point_group
        integer :: pg_code
        character(:), allocatable :: pg_name
        contains
        procedure :: get_point_group
        procedure :: get_rotational_group
        procedure :: get_stabilizer_group
        procedure :: get_reversal_group
    end type am_class_point_group

    type, public, extends(am_class_seitz_group) :: am_class_space_group
        contains
        procedure :: get_space_group
    end type am_class_space_group

contains

    ! operate on group representation

    subroutine     sort_symmetries(grp,criterion,flags)
        ! flags = 'ascend'/'descend'
        use am_rank_and_sort
        !
        implicit none
        !
        class(am_class_group), intent(inout) :: grp
        real(dp)    , intent(in) :: criterion(:)
        character(*), intent(in) :: flags
        integer     ,allocatable :: inds(:)
        integer     ,allocatable ::rinds(:) ! reverse inds
        !
        allocate(inds(grp%nsyms))
        if     (index(flags,'descend').ne.0) then
            call rank(-criterion,inds)
        elseif (index(flags,'ascend' ).ne.0) then
            call rank(+criterion,inds)
        else
            stop 'ERROR [sort_symmetries]: ascend/descend?'
        endif
        ! get reverse indices
        allocate(rinds(0:grp%nsyms))
        rinds(0)    = 0
        rinds(inds) = [1:grp%nsyms]
        ! sort symmetries
        select type (grp)
        class is (am_class_representation_group)
            ! [undefined units]  
            grp%sym = grp%sym(:,:,inds)
        class is (am_class_seitz_group)
            ! symmerties [frac.]
            grp%seitz_frac = grp%seitz_frac(:,:,inds)
            ! symmerties [cart.]
            grp%seitz_cart = grp%seitz_cart(:,:,inds)
            ! symmetries [recfrac]
            grp%seitz_recfrac = grp%seitz_recfrac(:,:,inds)
            ! regular representation
            if (allocated(grp%sym)) grp%sym = grp%sym(inds,inds,inds)
        class default
            stop 'ERROR [sort_symmetries]: class unknown'
        end select
        ! reshape(pack( .true.),shape())
        if (allocated(grp%ps_id))             grp%ps_id             =                            grp%ps_id(inds)
        if (allocated(grp%mt%multab))         grp%mt%multab         =         reshape(rinds(pack(grp%mt%multab(inds,inds)           , .true.)),shape(grp%mt%multab))
        if (allocated(grp%mt%corder))         grp%mt%corder         =                            grp%mt%corder(inds)        
        if (allocated(grp%mt%commutator_id))  grp%mt%commutator_id  = sort_nz(reshape(rinds(pack(grp%mt%commutator_id(inds,:)       , .true.)),shape(grp%mt%commutator_id)))
        if (allocated(grp%cc%id))             grp%cc%id             =                            grp%cc%id(inds)
        if (allocated(grp%cc%member))         grp%cc%member         =         reshape(rinds(pack(grp%cc%member                      , .true.)),shape(grp%cc%member))
        ! sort class matrices only for seitz group
        select type (grp)
        class is (am_class_seitz_group)
        if (allocated(grp%cc%matrices))       grp%cc%matrices       =                            grp%cc%matrices(inds,inds,:)
        end select
        ! generators need to be reevaluated. no other way. because generators are not unique.
        if (allocated(grp%mt%gen_id)) then
            ! symmetry generators
            grp%mt%gen_id = get_generator_id(multab=grp%mt%multab)
        endif
        if (allocated(grp%mt%gen)) then
            ! group generators
            deallocate(grp%mt%gen)
            allocate(grp%mt%gen, source=pack(grp%mt%gen_id,grp%mt%gen_id.ne.0) )
            grp%mt%gen = nint(sort(real(unique(grp%mt%gen),dp)))
        endif
        contains
        function     sort_nz(A) result(A_sorted)
            ! sort nonzero values in each row (puts zeros at the end of rows)
            implicit none
            !
            integer, intent(in) :: A(:,:)
            integer,allocatable :: A_sorted(:,:)
            logical,allocatable :: mask(:)
            integer :: m, n, i
            integer,allocatable :: inds(:)
            !
            m = size(A,1)
            n = size(A,2)
            !
            allocate(mask(n))
            !
            allocate(A_sorted(m,n))
            A_sorted = 0
            !
            allocate(inds(n))
            inds = [1:n]
            !
            do i = 1, m
                where (A(i,:).ne.0)
                    mask = .true.
                elsewhere
                    mask = .false.
                endwhere
                A_sorted(i,1:count(mask)) = sort( real(A(i, pack(inds,mask)) ,dp) )
            enddo
            !
        end function sort_nz
    end subroutine sort_symmetries

    subroutine     debug_dump_grp(grp,fname)
        !
        implicit none
        !
        class(am_class_group), intent(in) :: grp
        character(*), intent(in) :: fname
        !
                                              call dump(A=grp%nsyms            ,fname=trim(fname)//'.nsyms'            )
                                              call dump(A=grp%nbases           ,fname=trim(fname)//'.nbases'           )
        if (allocated(grp%sym              )) call dump(A=grp%sym              ,fname=trim(fname)//'.sym'              )
                                              call dump(A=grp%ps_id            ,fname=trim(fname)//'.ps_id'            )
        if (allocated(grp%supergroup_id    )) call dump(A=grp%supergroup_id    ,fname=trim(fname)//'.supergroup_id'    )
        if (allocated(grp%mt%multab        )) call dump(A=grp%mt%multab        ,fname=trim(fname)//'.mt.multab'        )
        if (allocated(grp%mt%gen_id        )) call dump(A=grp%mt%gen_id        ,fname=trim(fname)//'.mt.gen_id'        )
        if (allocated(grp%mt%gen           )) call dump(A=grp%mt%gen           ,fname=trim(fname)//'.mt.gen'           )
        if (allocated(grp%mt%inv_id        )) call dump(A=grp%mt%inv_id        ,fname=trim(fname)//'.mt.inv_id'        )
        if (allocated(grp%mt%corder        )) call dump(A=grp%mt%corder        ,fname=trim(fname)//'.mt.order'         )
        if (allocated(grp%mt%commutator_id )) call dump(A=grp%mt%commutator_id ,fname=trim(fname)//'.mt.commutator_id' )
                                              call dump(A=grp%cc%nclasses      ,fname=trim(fname)//'.cc.nclasses'      )
        if (allocated(grp%cc%id            )) call dump(A=grp%cc%id            ,fname=trim(fname)//'.cc.id'            )
        if (allocated(grp%cc%nelements     )) call dump(A=grp%cc%nelements     ,fname=trim(fname)//'.cc.nelements'     )
        if (allocated(grp%cc%member        )) call dump(A=grp%cc%member        ,fname=trim(fname)//'.cc.member'        )
        if (allocated(grp%cc%matrices      )) call dump(A=grp%cc%matrices      ,fname=trim(fname)//'.cc.matrices'      )
                                              call dump(A=grp%ct%nirreps       ,fname=trim(fname)//'.ct.nirreps'       )
        if (allocated(grp%ct%chartab       )) call dump(A=grp%ct%chartab       ,fname=trim(fname)//'.ct.chartab'       )
        if (allocated(grp%ct%irrep_dim     )) call dump(A=grp%ct%irrep_dim     ,fname=trim(fname)//'.ct.irrep_dim'     )
        if (allocated(grp%ct%irrep_decomp  )) call dump(A=grp%ct%irrep_decomp  ,fname=trim(fname)//'.ct.irrep_decomp'  )
        if (allocated(grp%ct%irrep_id      )) call dump(A=grp%ct%irrep_id      ,fname=trim(fname)//'.ct.irrep_id'      )
        if (allocated(grp%ct%irrep_proj    )) call dump(A=grp%ct%irrep_proj    ,fname=trim(fname)//'.ct.irrep_proj'    )
        if (allocated(grp%ct%irrep_proj_V  )) call dump(A=grp%ct%irrep_proj_V  ,fname=trim(fname)//'.ct.irrep_proj_V'  )
        if (allocated(grp%ct%block_proj    )) call dump(A=grp%ct%block_proj    ,fname=trim(fname)//'.ct.block_proj'    )
        ! dump group specific data                                                                                        
        select type (grp)                                                                                         
        class is (am_class_representation_group)                                                                                          
        if (allocated(grp%sym              )) call dump(A=grp%sym              ,fname=trim(fname)//'.seitz_cart'       )
        class is (am_class_seitz_group)                                                                           
        if (allocated(grp%seitz_cart       )) call dump(A=grp%seitz_cart       ,fname=trim(fname)//'.seitz_cart'       )
        if (allocated(grp%seitz_frac       )) call dump(A=grp%seitz_frac       ,fname=trim(fname)//'.seitz_frac'       )
        class default
            stop 'ERROR [debug_dump]: invalid group class'
        end select
    end subroutine debug_dump_grp

    subroutine     get_multiplication_table(grp)
        !
        implicit none
        !
        class(am_class_group), intent(inout) :: grp
        !
        ! get multiplication table
        select type (grp)
        class is (am_class_representation_group)
            ! [unitless]
            grp%mt%multab = get_multab(sym=grp%sym, flags='') !//'prog')
        class is (am_class_seitz_group)
            ! symmerties [frac.]
            grp%mt%multab = get_multab(sym=grp%seitz_frac, flags='seitz') !//'prog')
            ! get regular representation
            grp%sym       = get_regular_rep(multab=grp%mt%multab)
        class default
            stop 'ERROR [get_multiplication_table]: class unknown'
        end select
        ! get inverse indices
        grp%mt%inv_id = get_inverse_indices(multab=grp%mt%multab)
        ! get cyclic order
        grp%mt%corder = get_cyclic_order(multab=grp%mt%multab)
        ! get generator id
        grp%mt%gen_id = get_generator_id(multab=grp%mt%multab)
        ! get generators
        allocate(grp%mt%gen, source=pack(grp%mt%gen_id,grp%mt%gen_id.ne.0) )
        grp%mt%gen = nint(sort(real(unique(grp%mt%gen),dp)))
        ! identifty commutators
        grp%mt%commutator_id = get_commutator_id(multab=grp%mt%multab)
        !
    end subroutine get_multiplication_table

    subroutine     get_conjugacy_classes(grp)
        !
        implicit none
        !
        class(am_class_group), intent(inout) :: grp
        ! check
        if (.not.allocated(grp%mt%multab)) stop 'ERROR [get_conjugacy_classes]: multiplication table required'
        if (.not.allocated(grp%ps_id    )) stop 'ERROR [get_conjugacy_classes]: point symmetry identities required'
        if (.not.allocated(grp%mt%inv_id)) stop 'ERROR [get_conjugacy_classes]: inverse indices look-up table required'
        ! get conjugacy classes (ps_id is required for sorting the classes: proper before improper)
        grp%cc%id               = get_class_id(multab=grp%mt%multab, inv_id=grp%mt%inv_id, ps_id=grp%ps_id)
        ! number of classes
        grp%cc%nclasses         = maxval(grp%cc%id)
        ! get number of elements in each class
        grp%cc%nelements        = id_nelements(id=grp%cc%id)
        ! get members
        grp%cc%member           = id_member(id=grp%cc%id)
        ! get class matrices
        grp%cc%matrices         = get_class_matrices(sym=grp%sym, class_nelements=grp%cc%nelements, class_member=grp%cc%member)
    end subroutine get_conjugacy_classes

    subroutine     get_character_table(grp)
        !
        implicit none
        !
        class(am_class_group), intent(inout) :: grp
        complex(dp), allocatable :: chi(:)
        ! 
        ! check input
        if (.not.allocated(grp%mt%multab)) stop 'ERROR [get_character_table]: multiplication table is required'
        if (.not.allocated(grp%sym))       stop 'ERROR [get_character_table]: sym is required' ! for seitz group, regular rep can be saved as grp%sym.
        if (.not.allocated(grp%cc%member)) stop 'ERROR [get_character_table]: conjugacy classes are required'
        ! set number of irreps
        grp%ct%nirreps      = grp%cc%nclasses
        ! get character table (ps_id is required for sorting the irreps: proper before improper)
        grp%ct%chartab      = get_chartab(            multab=grp%mt%multab, class_nelements=grp%cc%nelements, class_member=grp%cc%member, ps_id=grp%ps_id)
        ! get muliken labels for irreps
        grp%ct%irrep_label  = get_irrep_label(      chartab=grp%ct%chartab,                                           class_id=grp%cc%id, ps_id=grp%ps_id,flags='number',ax_id=grp%ax_id)
        ! get irrep dimensions
        grp%ct%irrep_dim    = get_irrep_dimension(  chartab=grp%ct%chartab, class_nelements=grp%cc%nelements, class_member=grp%cc%member, ps_id=grp%ps_id)
        ! get rep character
        chi                 = get_rep_characters(                                                             class_member=grp%cc%member,     sym=grp%sym)
        ! get number of times each irrep appears in the decomposition of the representation
        grp%ct%irrep_decomp = get_irrep_decomp(     chartab=grp%ct%chartab, class_nelements=grp%cc%nelements,                                     chi=chi)
        ! identify which columns of proj_V and block_proj belong to which irrep
        grp%ct%irrep_id     = get_irrep_id(                                       irrep_dim=grp%ct%irrep_dim,            irrep_decomp=grp%ct%irrep_decomp)
        ! get irrep projection operators
        grp%ct%irrep_proj   = get_irrep_proj(       chartab=grp%ct%chartab,       irrep_dim=grp%ct%irrep_dim,         class_id=grp%cc%id,     sym=grp%sym)
        ! get irrep projection eigenvectors
        grp%ct%irrep_proj_V = get_irrep_proj_V(   irrep_id=grp%ct%irrep_id,     irrep_proj=grp%ct%irrep_proj)
        ! convert irrep projection eigenvectors to produce block-diagnaolized matrices
        grp%ct%block_proj   = get_block_proj(     irrep_id=grp%ct%irrep_id, irrep_proj_V=grp%ct%irrep_proj_V,                                 sym=grp%sym)
        !
    end subroutine get_character_table

    subroutine     print_character_table(grp)
        !
        use am_atom, only : l2spdf
        !
        implicit none
        !
        class(am_class_group), intent(in) :: grp
        type(am_class_chartab_printer) :: printer
        character(:), allocatable :: rep_label(:)
        integer :: i,j
        !
        ! print character table
        write(*,'(a,a)') flare, 'character table'
        ! initilize table
        call printer%printer_initialize()
        ! print header
        call printer%printer_print_chartab_header(class_nelements=grp%cc%nelements,class_member=grp%cc%member,ps_id=grp%ps_id)
        ! print irrep characters (character table)
        do i = 1, grp%ct%nirreps
            write(*,'(5x,i2,a8)',advance='no') i, grp%ct%irrep_label(i)
            do j = 1, grp%cc%nclasses
                write(*,printer%fmts(3),advance='no') printer%printer_chi2symb(chi=grp%ct%chartab(i,j))
            enddo
            write(*,*)
        enddo
        write(*,*)
        call printer%printer_print_definitions()
        ! ------------------------------------------------------------------------------------------
        ! print irrep decompositions
        write(*,'(a,a)') flare, 'decomposition of irrep products:'
        ! allocate space for reduction coefficient
        call print_product_decomp(chartab=grp%ct%chartab,class_nelements=grp%cc%nelements,&
            irrep_label=grp%ct%irrep_label,rep_label=grp%ct%irrep_label,chi=grp%ct%chartab)
        ! print group-specific things
        select type (grp)
        class is (am_class_point_group)
            ! ------------------------------------------------------------------------------------------
            ! initialize decomposition
            write(*,'(a,a)') flare, 'decomposition of symmeterized irrep products:'
            ! create rep label
            if (allocated(rep_label)) deallocate(rep_label)
            allocate(character(100) :: rep_label(grp%ct%nirreps))
            do i = 1, grp%ct%nirreps
                rep_label(i) = '['//trim(grp%ct%irrep_label(i))//' * '//trim(grp%ct%irrep_label(i))//']'
            enddo
            ! print decomposition
            call print_rep_decomp(chartab=grp%ct%chartab,class_nelements=grp%cc%nelements,&
                irrep_label=grp%ct%irrep_label,rep_label=rep_label,&
                chi=  get_chi_symmeterized_irrep(multab=grp%mt%multab,chartab=grp%ct%chartab,&
                              class_id=grp%cc%id,class_member=grp%cc%member,flags='plus') )
            ! ------------------------------------------------------------------------------------------
            ! initialize decomposition
            write(*,'(a,a)') flare, 'decomposition of antisymmeterized irrep products:'
            ! create rep label
            if (allocated(rep_label)) deallocate(rep_label)
            allocate(character(100) :: rep_label(grp%ct%nirreps))
            do i = 1, grp%ct%nirreps
                rep_label(i) = '{'//trim(grp%ct%irrep_label(i))//' * '//trim(grp%ct%irrep_label(i))//'}'
            enddo
            ! print decomposition
            call print_rep_decomp(chartab=grp%ct%chartab,class_nelements=grp%cc%nelements,&
                irrep_label=grp%ct%irrep_label,rep_label=rep_label,&
                chi=  get_chi_symmeterized_irrep(multab=grp%mt%multab,chartab=grp%ct%chartab,&
                              class_id=grp%cc%id,class_member=grp%cc%member,flags='minus') )
            ! ------------------------------------------------------------------------------------------
            ! initialize decomposition
            write(*,'(a,a)') flare, 'decompositions of subduced wigner matrices (orbital):'
            ! create rep label
            if (allocated(rep_label)) deallocate(rep_label)
            allocate(rep_label, source=['s','p','d','f'])
            ! print decomposition
            call print_product_decomp(chartab=grp%ct%chartab,class_nelements=grp%cc%nelements,&
                irrep_label=grp%ct%irrep_label,rep_label=rep_label,&
                chi=  get_chi_wigner(seitz_cart=grp%seitz_cart,class_member=grp%cc%member)  )
            ! ------------------------------------------------------------------------------------------
            ! initialize decomposition
            if (allocated(grp%ct%subduced_chi)) then
                ! initialize decomposition
                write(*,'(a,a)') flare, 'decompositions of subduced supergroup irrep products:'
                ! print decomposition
                call print_product_decomp(chartab=grp%ct%chartab,class_nelements=grp%cc%nelements,&
                    irrep_label=grp%ct%irrep_label,rep_label=grp%ct%subduced_label,chi=grp%ct%subduced_chi)
                write(*,'(a,a)') flare, 'decompositions of subduced supergroup irreps:'
                ! print decomposition
                call print_rep_decomp(chartab=grp%ct%chartab,class_nelements=grp%cc%nelements,&
                    irrep_label=grp%ct%irrep_label,rep_label=grp%ct%subduced_label,&
                    chi=  grp%ct%subduced_chi  )
            endif
            ! ------------------------------------------------------------------------------------------

            ! ADD VECTOR PRINT, RAMAN and IR selection rules.


            ! Irreps Decompositions
            ! Oh(m-3m)    A1g A1u A2g A2u Eu  Eg  T2u T2g T1u T1g
            ! V   ·   ·   ·   ·   ·   ·   ·   ·   1   ·
            ! [V2]    1   ·   ·   ·   ·   1   ·   1   ·   ·
            ! [V3]    ·   ·   ·   1   ·   ·   1   ·   2   ·
            ! [V4]    2   ·   ·   ·   ·   2   ·   2   ·   1
            ! A   ·   ·   ·   ·   ·   ·   ·   ·   ·   1
            ! [A2]    1   ·   ·   ·   ·   1   ·   1   ·   ·
            ! [A3]    ·   ·   1   ·   ·   ·   ·   1   ·   2
            ! [A4]    2   ·   ·   ·   ·   2   ·   2   ·   1
            ! [V2]xV  ·   ·   ·   1   1   ·   2   ·   3   ·
            ! [[V2]2] 3   ·   ·   ·   ·   3   ·   3   ·   1
            ! {V2}    ·   ·   ·   ·   ·   ·   ·   ·   ·   1
            ! {A2}    ·   ·   ·   ·   ·   ·   ·   ·   ·   1
            ! {[V2]2} ·   ·   1   ·   ·   1   ·   2   ·   2

            ! V ≡ the vector representation 
            ! A ≡ the axial representation 


            ! IR Selection Rules
            ! IR  A1g A1u A2g A2u Eu  Eg  T2u T2g T1u T1g
            ! A1g ·   ·   ·   ·   ·   ·   ·   ·   x   ·
            ! A1u ·   ·   ·   ·   ·   ·   ·   ·   ·   x
            ! A2g ·   ·   ·   ·   ·   ·   x   ·   ·   ·
            ! A2u ·   ·   ·   ·   ·   ·   ·   x   ·   ·
            ! Eu  ·   ·   ·   ·   ·   ·   ·   x   ·   x
            ! Eg  ·   ·   ·   ·   ·   ·   x   ·   x   ·
            ! T2u ·   ·   x   ·   ·   x   ·   x   ·   x
            ! T2g ·   ·   ·   x   x   ·   x   ·   x   ·
            ! T1u x   ·   ·   ·   ·   x   ·   x   ·   x
            ! T1g ·   x   ·   ·   x   ·   x   ·   x   ·

            ! [ Note: x means allowed ]


            ! Raman Selection Rules
            ! Raman   A1g A1u A2g A2u Eu  Eg  T2u T2g T1u T1g
            ! A1g x   ·   ·   ·   ·   x   ·   x   ·   ·
            ! A1u ·   x   ·   ·   x   ·   x   ·   ·   ·
            ! A2g ·   ·   x   ·   ·   x   ·   ·   ·   x
            ! A2u ·   ·   ·   x   x   ·   ·   ·   x   ·
            ! Eu  ·   x   ·   x   x   ·   x   ·   x   ·
            ! Eg  x   ·   x   ·   ·   x   ·   x   ·   x
            ! T2u ·   x   ·   ·   x   ·   x   ·   x   ·
            ! T2g x   ·   ·   ·   ·   x   ·   x   ·   x
            ! T1u ·   ·   ·   x   x   ·   x   ·   x   ·
            ! T1g ·   ·   x   ·   ·   x   ·   x   ·   x

            ! [ Note: x means allowed ]

        class is (am_class_space_group) 
            ! do nothing.
        class is (am_class_representation_group)
            ! do nothing.
        class default
            stop 'ERROR [print_character_table]: class unknown'
        end select
        !
        contains
        function       get_chi_wigner(seitz_cart,class_member) result(chi)
            !
            implicit none
            !
            real(dp), intent(in) :: seitz_cart(:,:,:)
            integer , intent(in) :: class_member(:,:)
            complex(dp),allocatable :: chi(:,:)
            integer :: nclasses
            ! initialize
            nclasses = size(class_member,1)
            allocate(chi(4,nclasses))
            ! representation based on s orbitals
            chi(1,:) = get_orb_characters(l=0, R=seitz_cart(1:3,1:3,:), class_member=class_member)
            ! representation based on p orbitals
            chi(2,:) = get_orb_characters(l=1, R=seitz_cart(1:3,1:3,:), class_member=class_member)
            ! representation based on d orbitals
            chi(3,:) = get_orb_characters(l=2, R=seitz_cart(1:3,1:3,:), class_member=class_member)
            ! representation based on f orbitals
            chi(4,:) = get_orb_characters(l=3, R=seitz_cart(1:3,1:3,:), class_member=class_member)
        end function   get_chi_wigner
        function       get_orb_characters(l,R,class_member) result(orbchi)
            !
            implicit none
            !
            integer,  intent(in) :: l
            real(dp), intent(in) :: R(:,:,:)
            integer,  intent(in) :: class_member(:,:)
            complex(dp), allocatable :: orbchi(:)
            integer :: nclasses
            integer :: i
            !
            nclasses = size(class_member,1)
            !
            allocate(orbchi(nclasses))
            do i = 1, nclasses
                orbchi(i) = get_orbchi(l=l,R=R(:,:,class_member(i,1)))
            enddo
            !
        end function   get_orb_characters
        function       get_orbchi(l,R) result(chi)
            !
            ! Character of (2l+1)-dimensional orthogonal group O(3) irrep parameterized by 3x3 rotoinversion matrix R
            !       
            ! T. Wolfram and Ş. Ellialtıoğlu, Applications of Group Theory to Atoms, Molecules, 
            ! and Solids, 1 edition (Cambridge University Press, Cambridge, 2014), p 74, Eq. 3.20.
            !
            ! M. Tinkham, Group Theory and Quantum Mechanics (McGraw-Hill, 1964), p 66 Eq 4.6 
            ! and p 100 bottom.
            !
            implicit none
            !
            integer,  intent(in) :: l ! dimension of irrep
            real(dp), intent(in) :: R(3,3) ! rotation angle
            real(dp) :: aa(4)
            real(dp) :: d
            complex(dp) :: chi
            !
            ! if R is a rotoinversion, get the angle and axis of the rotational part only (without the inversion)
            d = R(1,1)*R(2,2)*R(3,3)-R(1,1)*R(2,3)*R(3,2)-R(1,2)*R(2,1)*R(3,3)+R(1,2)*R(2,3)*R(3,1)+R(1,3)*R(2,1)*R(3,2)-R(1,3)*R(2,2)*R(3,1)
            ! convert to proper rotation
            aa = rot2axis_angle(d*R)
            if (abs(aa(4)).lt.tiny) then
                ! limiting case is obtained by expanding sin( (j+1/2)*x ) / sin( x/2 ) at zero
                chi = (2*l+1)
            else
                ! general case
                chi = sin( (l+0.5_dp)*aa(4) ) / sin( aa(4)*0.5_dp )
            endif
            ! convert back to improper rotation
            chi = chi * sign(1.0_dp,d)**(l)
            !
        end function   get_orbchi
        function       get_chi_symmeterized_irrep(multab,chartab,class_id,class_member,flags) result(chi)
            !
            implicit none
            !
            integer    , intent(in) :: multab(:,:)
            complex(dp), intent(in) :: chartab(:,:)
            integer    , intent(in) :: class_id(:)
            integer    , intent(in) :: class_member(:,:)
            character(*),intent(in) :: flags
            complex(dp),allocatable :: chi(:,:)
            integer :: nirreps, nclasses
            integer :: class_rep
            integer :: i,j,ii
            !
            nirreps = size(chartab,1)
            nclasses= size(chartab,2)
            !
            allocate(chi(nirreps,nclasses))
            !
            do j = 1, nirreps
            do i = 1, nclasses
                ! get class rep
                class_rep = class_member(i,1)
                ! get class of rep squared
                ii = class_id(multab(class_rep,class_rep))
                ! get antisymmetric character
                if     (index(flags,'plus').ne.0) then
                    ! symmetric
                    chi(j,i) = 0.5_dp * (chartab(j,i)**2 + chartab(j,ii) )
                elseif (index(flags,'minus').ne.0) then
                    ! antisymmetric
                    chi(j,i) = 0.5_dp * (chartab(j,i)**2 - chartab(j,ii) )
                else 
                    stop 'ERROR [get_chi_symmeterized_irrep]: unknown flag'
                endif
            enddo
            enddo
            !
        end function   get_chi_symmeterized_irrep
        subroutine     print_product_decomp(chartab,class_nelements,irrep_label,rep_label,chi)
            !
            implicit none
            !
            complex(dp) , intent(in) :: chartab(:,:) ! class constant chartab which becomes character chartab (can be complex!)
            integer     , intent(in) :: class_nelements(:)
            character(*), intent(in) :: irrep_label(:)
            character(*), intent(in) :: rep_label(:)
            complex(dp) , intent(in) :: chi(:,:)
            integer :: i, j, k
            integer :: nirreps
            integer :: nclasses
            integer :: nreps
            character(:), allocatable :: str(:)
            integer, allocatable :: beta(:,:,:)
            integer :: maxlen
            ! get dimensions
            nirreps  = size(chartab,1)
            nclasses = nirreps
            nreps    = size(chi,1)
            ! allocate space for string
            allocate(character(150)::str(nreps))
            ! allocate space for reduction cofficients
            allocate(beta(nclasses,nreps,nreps))
            ! initialize
            beta = 0
            ! get reduction coefficients
            do i = 1, nreps
            do j = 1, i
                beta(:,i,j) = get_reduction_coefficient(chartab=chartab,class_nelements=class_nelements,chi_rep=chi(i,:)*chi(j,:))
            enddo
            enddo
            ! get string and begin printing
            call disp_indent()
            call disp(X=rep_label,title='*',style='underline',advance='no',fmt='10a')
            do i = 1, nreps
                do j = 1, nreps
                    if (all(beta(:,i,j).eq.0)) then
                        str(j) = '.'
                    else
                        str(j) = ''
                        do k = 1, nirreps
                        if (beta(k,i,j).ne.0) then
                            if     (beta(k,i,j).eq.+1) then
                                if (len_trim(str(j)).eq.0) then
                                    ! do nothing
                                else
                                    str(j) = trim(str(j))//'+'
                                endif
                            else
                                str(j) = trim(str(j))//trim(tostring(X=beta(k,i,j),fmt='SP,i5'))
                            endif
                            str(j) = trim(str(j))//trim(irrep_label(k))
                        endif
                        enddo
                    endif
                enddo
                ! get max length
                maxlen = len_trim(rep_label(i))
                do j = 1, nreps
                    if (len_trim(str(j)).gt.maxlen) maxlen = len_trim(str(j))
                enddo
                ! print column of table
                call disp(X=trim_null(str(1:nreps)),title=trim(rep_label(i)),style='underline',advance='no',fmt='a'//tostring(maxlen))
            enddo
            ! write table
            call disp(X=0,zeroas=' ',advance='yes')
            !
        end subroutine print_product_decomp
        subroutine     print_rep_decomp(chartab,class_nelements,irrep_label,rep_label,chi)
            !
            implicit none
            !
            complex(dp) , intent(in) :: chartab(:,:) ! class constant chartab which becomes character chartab (can be complex!)
            integer     , intent(in) :: class_nelements(:)
            character(*), intent(in) :: irrep_label(:)
            character(*), intent(in) :: rep_label(:)
            complex(dp) , intent(in) :: chi(:,:)
            integer :: i, j, k
            integer :: nirreps
            integer :: nclasses
            integer :: nreps
            character(:), allocatable :: str
            integer, allocatable :: beta(:,:)
            ! get dimensions
            nirreps = size(chartab,1)
            nclasses = nirreps
            nreps = size(chi,1)
            ! allocate space for string
            allocate(character(150)::str)
            ! allocate space for reduction cofficients
            allocate(beta(nreps,nreps))
            ! initialize
            beta = 0
            ! print rep decompositions
            do i = 1, nreps
                beta(:,i) = get_reduction_coefficient(chartab=chartab,class_nelements=class_nelements,chi_rep=chi(i,:))
            enddo
            ! print
            do i = 1, nreps
            if (any(beta(:,i).ne.0)) then
                write(*,'(5x,a,a)',advance='no') trim(rep_label(i)), ' = '
                str = ''
                do j = 1, nirreps
                if (beta(j,i).ne.0) then
                    if     (beta(j,i).eq.1) then
                        if (len_trim(str).eq.0) then
                            ! do nothing
                        else
                            str = '+'
                        endif
                    else
                        str = trim(tostring(X=beta(j,i),fmt='SP,i5'))
                    endif
                    write(*,'(a)',advance='no') trim(str)//trim(irrep_label(j))
                endif
                enddo
                write(*,*)
            endif
            enddo
            if (all(beta.eq.0)) then
                write(*,'(5x,a)') 'all reduction coefficients are zero'
            endif
        end subroutine print_rep_decomp
    end subroutine print_character_table

    ! high level routines which operate on sg/pg

    subroutine     create_seitz_group(sg,bas,seitz_frac)
        !
        implicit none
        !
        class(am_class_seitz_group), intent(out):: sg
        real(dp), intent(in) :: bas(3,3)
        real(dp), intent(in) :: seitz_frac(:,:,:)
        integer  :: i
        !
        ! basis in which [frac.] are defined
        sg%bas = bas
        ! get reciprocal basis in which [recfrac] are defined
        sg%recbas = inv(bas)
        ! get number of symmetries
        sg%nsyms = size(seitz_frac,3)
        ! number of "bases functions", doesn't have any real meaning here.
        ! using nsyms since regular rep is used in decomposing seitz groups
        sg%nbases = sg%nsyms
        ! get symmetries [frac.]
        allocate(sg%seitz_frac,source=seitz_frac)
        ! put identity first
        call put_identity_first(sg%seitz_frac)
        ! get symmetries [cart.]
        allocate(sg%seitz_cart(4,4,sg%nsyms))
        do i = 1, sg%nsyms
            sg%seitz_cart(:,:,i) = seitz_frac2cart(seitz_frac=sg%seitz_frac(:,:,i), bas=sg%bas, recbas=sg%recbas)
        enddo
        ! get symmetries [rec. frac.]
        allocate(sg%seitz_recfrac(4,4,sg%nsyms))
        do i = 1, sg%nsyms
            sg%seitz_recfrac(:,:,i) = seitz_cart2recfrac(seitz_cart=sg%seitz_cart(:,:,i), bas=sg%bas, recbas=sg%recbas)
        enddo
        ! correct rounding error in frac
        where(abs(sg%seitz_frac).lt.tiny) sg%seitz_frac = 0
        ! correct rounding error in cart
        call correct_rounding_error(sg%seitz_cart)
        ! identify point symmetries
        sg%ps_id = get_ps_id(R=sg%seitz_cart)
        ! get ax_id
        sg%ax_id = get_ax_id(R=sg%seitz_cart)
        ! sort symmetries based on parameters 
        call sg%sort_symmetries(criterion=sg%seitz_frac(1,4,:)      , flags='ascend')
        call sg%sort_symmetries(criterion=sg%seitz_frac(2,4,:)      , flags='ascend')
        call sg%sort_symmetries(criterion=sg%seitz_frac(3,4,:)      , flags='ascend')
        call sg%sort_symmetries(criterion=real(sg%ps_id,dp)         , flags='ascend')
        call sg%sort_symmetries(criterion=cc(sg%seitz_frac,sg%ps_id), flags='ascend') ! cc_id wrapper
        ! get multiplication table
        call sg%get_multiplication_table()
        ! get conjugacy classes
        call sg%get_conjugacy_classes()
        ! get character table
        call sg%get_character_table()
        ! get class-specifc things
        select type (sg)
        class is (am_class_point_group)
            ! identify point group
            sg%pg_code = get_pg_code(sg%ps_id)
            ! get name of point group
            sg%pg_name = trim_null(get_pg_name(sg%pg_code))
        class is (am_class_space_group)
            ! do nothing
        class default
            stop 'ERROR [create_seitz_group]: class unknown'
        end select
        !
        contains
        function       seitz_frac2cart(seitz_frac,bas,recbas) result(seitz_cart)
            !
            implicit none
            !
            real(dp), intent(in) :: seitz_frac(4,4)
            real(dp), intent(in) :: bas(3,3)
            real(dp), intent(in) :: recbas(3,3)
            real(dp) :: seitz_cart(4,4)
            !
            seitz_cart(1:3,1:3) = matmul(bas,matmul(seitz_frac(1:3,1:3),recbas))
            seitz_cart(1:3,4)   = matmul(bas,seitz_frac(1:3,4))
            seitz_cart(4,4)     = 1.0_dp
            seitz_cart(4,1:3)   = 0.0_dp
            !
        end function   seitz_frac2cart
        function       seitz_cart2recfrac(seitz_cart,bas,recbas) result(seitz_recfrac)
            !
            implicit none
            !
            real(dp), intent(in) :: seitz_cart(4,4)
            real(dp), intent(in) :: bas(3,3)
            real(dp), intent(in) :: recbas(3,3)
            real(dp) :: seitz_recfrac(4,4)
            !
            seitz_recfrac(1:3,1:3) = matmul(bas,matmul(seitz_cart(1:3,1:3),recbas))
            seitz_recfrac(1:3,4)   = matmul(bas,seitz_cart(1:3,4))
            seitz_recfrac(4,4)     = 1.0_dp
            seitz_recfrac(4,1:3)   = 0.0_dp
            !
        end function   seitz_cart2recfrac
        function       cc(seitz_frac,ps_id) result(cc_id)
            !
            implicit none
            !
            real(dp), intent(in) :: seitz_frac(:,:,:)
            integer , intent(in) :: ps_id(:)
            real(dp),allocatable :: cc_id(:)
            integer ,allocatable :: multab(:,:)
            !
            multab = get_multab(sym=seitz_frac, flags='seitz')
            cc_id  = real(get_class_id(multab=multab, inv_id=get_inverse_indices(multab), ps_id=ps_id),dp)
        end function   cc
    end subroutine create_seitz_group

    subroutine     get_space_group(sg,pc,opts)
        !
        implicit none
        !
        class(am_class_space_group), intent(inout) :: sg
        class(am_class_unit_cell)  , intent(in) :: pc ! space groups are tabulated in the litearture for conventional cells, but primitive or arbitrary cell works just as well.
        type(am_class_options)     , intent(in) :: opts
        real(dp)    , allocatable :: seitz_frac(:,:,:)
        character(20) :: str
        integer :: mpl
        integer :: i
        !
        if (opts%verbosity.ge.1) call print_title('Space group symmetries')
        !
        ! determine space symmetries from atomic basis [frac.]
        seitz_frac = space_symmetries_from_basis(bas=pc%bas, tau=pc%tau_frac,Z=pc%Z, prec=opts%prec, verbosity=opts%verbosity)
        ! put identity first
        call put_identity_first(seitz=seitz_frac)
        ! create space group instance
        call sg%create_seitz_group(seitz_frac=seitz_frac, bas=pc%bas)
        ! print stdout
        if (opts%verbosity.ge.1) then
            ! print global information
            write(*,'(a,a,a)') flare, 'space symmetries = ' , tostring(sg%nsyms)
            write(*,'(a,a,a)') flare, 'conjugacy classes = ', tostring(sg%cc%nclasses)
            write(*,'(a,a,a)') flare, 'group generators = ' , tostring(sg%mt%gen)
            write(*,'(a,a)')   flare, 'seitz symmetries [frac/cart] = '
            ! print seitz operators
            mpl=2
            do i = 1,sg%nsyms
                str = tostring(i)//': '//get_ps_name(sg%ps_id(i))
                if     ((mod(i,mpl).eq.0).or.(i.eq.sg%nsyms)) then
                    call disp_indent()
                    call disp(title=trim(str)//' [frac]', X=sg%seitz_frac(:,:,i) ,style='underline',fmt='f5.2',zeroas='0',advance='no' , trim='no')
                    call disp(title=trim(str)//' [cart]', X=sg%seitz_cart(:,:,i) ,style='underline',fmt='f5.2',zeroas='0',advance='yes', trim='no')
                else
                    call disp_indent()
                    call disp(title=trim(str)//' [frac]', X=sg%seitz_frac(:,:,i) ,style='underline',fmt='f5.2',zeroas='0',advance='no' , trim='no')
                    call disp(title=trim(str)//' [cart]', X=sg%seitz_cart(:,:,i) ,style='underline',fmt='f5.2',zeroas='0',advance='no ', trim='no')
                endif
            enddo
            ! print multiplication table
            call print_title('Space group multiplication table')
            call disp_indent()
            call disp(X=sg%mt%multab,advance='yes',trim='yes')
        endif
        ! dump debugging
        if (debug) then
        call execute_command_line('mkdir -p '//trim(debug_dir)//'/space_group/')
        call sg%debug_dump(fname=              trim(debug_dir)//'/space_group/outfile.sg')
        endif
        ! write action table
        call execute_command_line('mkdir -p '//trim(outfile_dir_sym))
        call sg%write_action_table(uc=pc,fname=trim(outfile_dir_sym)//'/'//'outfile.space_group_action',opts=opts)
        !
        contains
        function        space_symmetries_from_basis(bas,tau,Z,prec,verbosity) result(seitz)
            !
            !The program first finds the primitive unit cell by looking for
            !additional lattice vectors within the unit cell given in the input.
            !
            ! It chooses one of the atoms and tests each vector from that atom to
            !every other atom of the same type in the unit cell. If that vector
            !will take us from each atom in the unit cell to another atom of the
            !same type (not necessarily in the same unit cell), then we have
            !found a new lattice vector.
            !
            ! This process is repeated until we have
            !found a complete list of lattice vectors. From that list, we find three
            !which can serve as primitive basis vectors ti of the lattice, i.e. every
            !lattice vector can be expressed as an integer linear combination of
            !these basis vectors.
            !
            use am_rank_and_sort
            !
            implicit none
            ! subroutine i/o
            real(dp), intent(in) :: bas(3,3)
            real(dp), intent(in) :: tau(:,:)
            integer , intent(in) :: Z(:)
            real(dp), intent(in) :: prec
            integer , intent(in) :: verbosity
            real(dp), allocatable :: seitz(:,:,:)
            real(dp), allocatable :: try(:,:)
            real(dp), allocatable :: T(:,:)
            real(dp), allocatable :: R(:,:,:)
            real(dp), allocatable :: wrkspace(:,:,:)
            real(dp) :: shift(3)
            real(dp) :: T_shifted(3)
            integer :: nTs
            integer :: nRs
            integer :: i, j, k
            !
            R = lattice_symmetries(bas=bas, prec=prec)
            nRs = size(R,3)
            !
            T = translations_from_basis(tau=tau, Z=Z, prec=prec, flags='zero,relax')
            nTs = size(T,2)
            !
            allocate(wrkspace(4,4,nTs*nRs))
            wrkspace = 0.0_dp
            k=0
            do i = 1, nRs
            do j = 1, nTs
                ! Shift tau to put the first atom at the origin. This is important because the
                ! rotational part of the symmetry operation must leave at least one atom in it's
                ! original location. If the atom already has a displacement, then it will be rotated
                ! to another position.
                if (.not.any(all(abs(tau).lt.prec,2))) then
                    shift = tau(:,1)
                    ! tau' = R*tau + T
                    ! tau' - s = R*(tau-s) + T
                    ! tau' = R*tau + (R*S+T-s) = R*tau + T_shifted ! sign is wrong here for some reason.
                    ! thus, T_shifted = T - R*s + s
                    T_shifted = T(:,j) - matmul(R(:,:,i),shift) + shift
                    ! reduce to primitive
                    T_shifted = modulo(T_shifted+opts%prec,1.0_dp)-opts%prec
                else
                    T_shifted = T(:,j)
                endif
                !
                try = eye(4)
                try(1:3,1:3) = R(1:3,1:3,i)
                try(1:3,4) = T_shifted
                !
                if ( is_symmetry_valid(tau=tau, Z=Z, seitz=try, prec=opts%prec)) then
                    k = k + 1
                    wrkspace(:,:,k)=try
                endif
            enddo
            enddo
            !
            allocate(seitz(4,4,k))
            seitz = wrkspace(1:4,1:4,1:k)
            !
            ! print stdout
            if (verbosity.ge.2) then
                call am_print('possible (im-)propers',nRs,' ... ')
                call am_print('possible translations',nTs,' ... ')
                call am_print_two_matrices_side_by_side(name='translations',&
                    Atitle='fractional',A=transpose(T),&
                    Btitle='cartesian' ,B=transpose(matmul(bas,T)),&
                    iopt_emph=' ... ',iopt_teaser=.true.)
            endif
        end function    space_symmetries_from_basis
        function        lattice_symmetries(bas,prec) result(R_frac)
            !>
            !> Given a metric tensor, this function returns the (arthimetic) a-holohodry:
            !> the group of point symmetries R (fractional) that are compatible with the
            !> metric tensor).
            !>
            !> Symmetry and Condensed Matter Physics: A Computational Approach. 1 edition.
            !> Cambridge, UK; New York: Cambridge University Press, 2008. pages 278-279,
            !> Eq. 10.60 and 10.61
            !>
            !> The Mathematical Theory of Symmetry in Solids: Representation Theory
            !> for Point Groups and Space Groups. 1 edition. Oxford?: New York:
            !> Oxford University Press, 2010. page 134, table 3.9
            !>
            !> Structure of Materials: An Introduction to Crystallography, Diffraction
            !> and Symmetry. 2 edition. New York: Cambridge University Press, 2012.
            !>
            !> page 169:
            !> "If [point symmetries] are described with respect to the 'primed' basis
            !>  vectors (fractional coordinates), then once again all [point symmetries]
            !>  can be represented by matrices which contain only -1, 0, +1"
            !>
            !> page 166:
            !> "The inverse of a unitary matrix (rotation) is equal to its transpose
            !>  (this is only true in cartesian coordinates.)"
            !>
            implicit none
            ! subroutine i/o
            real(dp), intent(in) :: bas(3,3)
            real(dp), intent(in) :: prec
            real(dp), allocatable :: R_frac(:,:,:) !> point symmetries (fractional)
            !
            real(dp) :: id(3,3)
            real(dp) :: recbas(3,3)
            real(dp) :: metric(3,3)
            real(dp) :: buffer(3,3,48) ! buffer for point symmetries
            integer  :: k ! point symmetry counter
            real(dp) :: o(3,3)  ! point symmetry in fractional
            integer  :: i11, i12, i13, i21, i22, i23, i31, i32, i33 ! used to generate unitary rotational matrices
            !
            id = eye(3)
            !
            recbas = inv(bas)
            !
            metric = matmul(transpose(bas),bas)
            !
            buffer = 0.0_dp
            !
            k = 0
            !
            do i11 = -1,1
            do i12 = -1,1
            do i13 = -1,1
                do i21 = -1,1
                do i22 = -1,1
                do i23 = -1,1
                    do i31 = -1,1
                    do i32 = -1,1
                    do i33 = -1,1
                        o(1,1:3)=real([i11,i12,i13],dp)
                        o(2,1:3)=real([i21,i22,i23],dp)
                        o(3,1:3)=real([i31,i32,i33],dp)
                        ! Check that metric is left unchanged by symmetry opreation in fractional coordinates.
                        ! i.e. that metric tensor commutes with point symmetry.
                        if ( all( abs(matmul(transpose(o),matmul(metric,o))-metric) .lt. prec ) ) then
                            k = k + 1
                            ! store point symmetry in fractional coordinates
                            buffer(:,:,k) = o
                        endif
                        !
                    enddo
                    enddo
                    enddo
                enddo
                enddo
                enddo
            enddo
            enddo
            enddo
            !
            allocate(R_frac,source=buffer(:,:,1:k))
            !
        end function    lattice_symmetries
    end subroutine get_space_group

    subroutine     get_point_group(pg,pc,sg,opts)
        !
        implicit none
        !
        class(am_class_point_group), intent(inout) :: pg
        type(am_class_space_group) , intent(in) :: sg
        class(am_class_unit_cell)  , intent(in) :: pc  ! accepts unit cell, primitive cell, conventional cell.
        type(am_class_options)     , intent(in) :: opts
        real(dp), allocatable  :: wrk_frac(:,:,:)
        character(20) :: str
        integer :: mpl
        integer :: i
        !
        if (opts%verbosity.ge.1) call print_title('Point group symmetries')
        !
        ! create point group instance
        allocate(wrk_frac,source=sg%seitz_frac)
        wrk_frac(1:3,4,:) = 0
        call pg%create_seitz_group(seitz_frac=unique(wrk_frac), bas=pc%bas)
        !
        if (opts%verbosity.ge.1) then
            write(*,'(a,a,a)') flare, 'point group = '      , get_pg_name(get_pg_code(pg%ps_id))
            write(*,'(a,a,a)') flare, 'point symmetries = ' , tostring(pg%nsyms)
            write(*,'(a,a,a)') flare, 'conjugacy classes = ', tostring(pg%cc%nclasses)
            write(*,'(a,a)')   flare, 'group generators = ' , tostring(pg%mt%gen)
            write(*,'(a,a)')   flare, 'point symmetries [frac/cart] = '
            mpl=2 ! matrices per line
            do i = 1,pg%nsyms
                str = tostring(i)//': '//get_ps_name(pg%ps_id(i))
                if ((mod(i,mpl).eq.0).or.(i.eq.sg%nsyms)) then
                    call disp_indent()
                    call disp(title=trim(str)//' [frac]', X=pg%seitz_frac(:,:,i) ,style='underline',fmt='f5.2',zeroas='0',advance='no' , trim='no')
                    call disp(title=trim(str)//' [cart]', X=pg%seitz_cart(:,:,i) ,style='underline',fmt='f5.2',zeroas='0',advance='yes', trim='no')
                else
                    call disp_indent()
                    call disp(title=trim(str)//' [frac]', X=pg%seitz_frac(:,:,i) ,style='underline',fmt='f5.2',zeroas='0',advance='no' , trim='no')
                    call disp(title=trim(str)//' [cart]', X=pg%seitz_cart(:,:,i) ,style='underline',fmt='f5.2',zeroas='0',advance='no ', trim='no')
                endif
            enddo
            ! print multiplication table
            call print_title('Point group multiplication table')
            call disp(spread(' ',1,2),orient='row',advance='no',trim='no') ! add some space
            call disp(X=pg%mt%multab,advance='yes',trim='yes')
            ! print character table
            call print_title('Point group character properties')
            call pg%print_character_table()
        endif
        ! dump debugging
        if (debug) then
        call execute_command_line ('mkdir -p '//trim(debug_dir)//'/point_group')
        call pg%debug_dump(fname=               trim(debug_dir)//'/point_group/outfile.pg')
        endif
        ! write action table
        call pg%write_action_table(uc=pc,fname=trim(outfile_dir_sym)//'/'//'outfile.point_group_action',opts=opts)
        !
    end subroutine get_point_group

    subroutine     get_rotational_group(rotg,uc,pg,opts,flags)
        !
        ! Get point symmetries which are compatible with the atomic basis (essentially keeps
        ! rotational part of space symmetries which are able to map all atoms onto each other)
        !
        ! This is also useful for determining the point symmetries compatible with kpoint mesh
        !
        ! Name "point group" - NOTE: not really a point group because translational components of
        ! space symmetries were ignored in the genereation process this would be the name of the
        ! group, if the group were symmorphic (in which case it will match the real point group
        ! name)...
        !
        implicit none
        ! subroutine i/o
        class(am_class_point_group),intent(out):: rotg
        class(am_class_unit_cell),  intent(in) :: uc
        type(am_class_point_group), intent(in) :: pg
        type(am_class_options),     intent(in) :: opts
        character(*)              , intent(in) :: flags
        integer, allocatable   :: indicies(:)
        logical, allocatable   :: mask(:)
        integer :: i,j
        !
        allocate(mask(pg%nsyms))
        mask = .false. 
        !
        if      (index(flags,'cart').ne.0) then
            do i = 1, pg%nsyms
                if ( is_symmetry_valid(tau=uc%tau_cart, Z=uc%Z, seitz=pg%seitz_cart, prec=opts%prec, flags='cart')) then
                    mask(i) = .true.
                endif
            enddo
        else
            stop 'Note: rotational group only makes sense for [cart], for which symmetry operations are strictly unitary.'
        endif
        ! create indices
        allocate(indicies(pg%nsyms))
        indicies=[1:pg%nsyms]
        !
        ! create group
        call rotg%create_seitz_group(seitz_frac=pg%seitz_frac(:,:,pack(indicies,mask)), bas=pg%bas)
        ! save id's of corresponding symmetries in the encompasing group
        allocate(rotg%supergroup_id(rotg%nsyms))
        do i = 1, rotg%nsyms
        search : do j = 1, pg%nsyms
            if (isequal(rotg%seitz_cart(:,:,i),pg%seitz_cart(:,:,j))) then
                rotg%supergroup_id(i) = j
                exit search
            endif
        enddo search
        enddo
        ! get subduced chi
        rotg%ct%subduced_chi = get_subduced_chi(supergroup_chartab=pg%ct%chartab,&
            supergroup_class_id=pg%cc%id,supergroup_id=rotg%supergroup_id)
        ! get subduced irrep labels
        allocate(character(200)::rotg%ct%subduced_label(pg%ct%nirreps))
        do i = 1, pg%ct%nirreps
            rotg%ct%subduced_label(i) = trim(pg%ct%irrep_label(i))//'('//trim(pg%pg_name)//')'
        enddo
        rotg%ct%subduced_label = trim_null(rotg%ct%subduced_label)
        !
    end subroutine get_rotational_group

    subroutine     get_stabilizer_group(stab,pg,v,opts,flags)
        !
        implicit none
        !
        class(am_class_point_group), intent(out):: stab ! stabilizer group associated with vector v
        type(am_class_point_group) , intent(in) :: pg
        real(dp)                   , intent(in) :: v(3) ! vector which is stabilized (should be same units as R, which is fractional)
        type(am_class_options)     , intent(in) :: opts
        character(*)               , intent(in) :: flags
        integer, allocatable :: indicies(:)
        logical, allocatable :: mask(:)
        integer :: i,j
        !
        allocate(mask(pg%nsyms))
        mask = .false.
        !
        if      (index(flags,'frac').ne.0) then
            do i = 1, pg%nsyms
                if (is_symmetry_valid(tau=reshape(v,[3,1]), Z=[1], seitz=pg%seitz_frac(:,:,i), prec=opts%prec, flags=flags//'exact')) then
                    mask(i) = .true.
                endif
            enddo
        elseif  (index(flags,'cart').ne.0) then
            do i = 1, pg%nsyms
                if (is_symmetry_valid(tau=reshape(v,[3,1]), Z=[1], seitz=pg%seitz_cart(:,:,i), prec=opts%prec, flags=flags//'exact')) then
                    mask(i) = .true.
                endif
            enddo
        else
            stop 'flags must contain frac/cart'
        endif
        !
        allocate(indicies(pg%nsyms))
        indicies=[1:pg%nsyms]
        ! create group
        call stab%create_seitz_group(seitz_frac=pg%seitz_frac(:,:,pack(indicies,mask)), bas=pg%bas)
        ! save id's of corresponding symmetries in the encompasing group
        allocate(stab%supergroup_id(stab%nsyms))
        do i = 1, stab%nsyms
        search : do j = 1, pg%nsyms
            if (isequal(stab%seitz_cart(:,:,i),pg%seitz_cart(:,:,j))) then
                stab%supergroup_id(i) = j
                exit search
            endif
        enddo search
        enddo
        ! get subduced chi
        stab%ct%subduced_chi = get_subduced_chi(supergroup_chartab=pg%ct%chartab,&
            supergroup_class_id=pg%cc%id,supergroup_id=stab%supergroup_id)
        ! get subduced irrep labels
        allocate(character(200)::stab%ct%subduced_label(pg%ct%nirreps))
        do i = 1, pg%ct%nirreps
            stab%ct%subduced_label(i) = trim(pg%ct%irrep_label(i))//'('//trim(pg%pg_name)//')'
        enddo
        stab%ct%subduced_label = trim_null(stab%ct%subduced_label)
        !
    end subroutine get_stabilizer_group
    
    subroutine     get_reversal_group(revg,pg,v,opts,flags)
        !
        ! there is a bug in the reversal group. Si-Si first nearest neighbor shells on different primitive cell atoms does not match.
        !
        !  shell Zi-Zj   i-j   m-n    m   stab.    rot.    rev. |v(cart)|           v(cart)             |v(frac)|           v(frac)
        !  ----- ----- ----- ----- ---- ------- ------- ------- --------- ----------------------------- --------- -----------------------------
        !      1 Si-Si   1-1   1-1    1     o_h     o_h     t_d     0.000    -0.000    -0.000    -0.000     0.000    -0.000    -0.000    -0.000
        !      2 Si-Si   1-1   1-2    4    c_3v     o_h    d_3d     2.351     1.357     1.357     1.357     0.217     0.125     0.125     0.125
        !  primitive atom 2 at 0.12,0.12,0.12 (frac) has 2 nearest-neighbor shells
        !  shell Zi-Zj   i-j   m-n    m   stab.    rot.    rev. |v(cart)|           v(cart)             |v(frac)|           v(frac)
        !  ----- ----- ----- ----- ---- ------- ------- ------- --------- ----------------------------- --------- -----------------------------
        !      3 Si-Si   1-1   2-2    1     o_h     o_h     t_d     0.000     0.000     0.000     0.000     0.000     0.000     0.000     0.000
        !      4 Si-Si   1-1   2-1    4    c_3v     o_h    c_3v     2.351    -1.357    -1.357    -1.357     0.217    -0.125    -0.125    -0.125
        !
        implicit none
        !
        class(am_class_point_group), intent(out):: revg ! reversal group associated with vector v, i.e. the space symmetries with bring v to [0,0,0], essentially flipping the atoms at the corners of the vector
        type(am_class_point_group) , intent(in) :: pg ! space group
        real(dp)                   , intent(in) :: v(3)
        type(am_class_options)     , intent(in) :: opts
        character(*)               , intent(in) :: flags
        integer, allocatable :: indicies(:)
        logical, allocatable :: mask(:)
        real(dp),allocatable :: seitz_frac(:,:,:)
        integer :: i,j
        !
        allocate(mask(pg%nsyms))
        mask = .false.
        !
        if      (index(flags,'frac').ne.0) then
            do i = 1, pg%nsyms
                if (is_symmetry_valid(tau=reshape(v,[3,1]), Z=[1], seitz=pg%seitz_frac(:,:,i), prec=opts%prec, flags=flags//'reversal')) then
                    mask(i) = .true.
                endif
            enddo
        elseif  (index(flags,'cart').ne.0) then
            do i = 1, pg%nsyms
                if (is_symmetry_valid(tau=reshape(v,[3,1]), Z=[1], seitz=pg%seitz_cart(:,:,i), prec=opts%prec, flags=flags//'reversal')) then
                    mask(i) = .true.
                endif
            enddo
        else
            stop 'flags must contain frac/cart'
        endif
        !
        allocate(indicies(pg%nsyms))
        indicies=[1:pg%nsyms]
        !
        allocate(seitz_frac, source=pg%seitz_frac(:,:,pack(indicies,mask)))
        seitz_frac(1:3,4,:) = 0 ! set all translations to zero
        !
        ! create group
        call revg%create_seitz_group(seitz_frac=unique(seitz_frac), bas=pg%bas)
        ! save id's of corresponding symmetries in the encompasing group
        allocate(revg%supergroup_id(revg%nsyms))
        do i = 1, revg%nsyms
        search : do j = 1, pg%nsyms
            if (isequal(revg%seitz_cart(:,:,i),pg%seitz_cart(:,:,j))) then
                revg%supergroup_id(i) = j
                exit search
            endif
        enddo search
        enddo
        ! get subduced chi
        revg%ct%subduced_chi = get_subduced_chi(supergroup_chartab=pg%ct%chartab,&
            supergroup_class_id=pg%cc%id,supergroup_id=revg%supergroup_id)
        ! get subduced irrep labels
        allocate(character(200)::revg%ct%subduced_label(pg%ct%nirreps))
        do i = 1, pg%ct%nirreps
            revg%ct%subduced_label(i) = trim(pg%ct%irrep_label(i))//'('//trim(pg%pg_name)//')'
        enddo
        revg%ct%subduced_label = trim_null(revg%ct%subduced_label)
        !
    end subroutine get_reversal_group

    subroutine     write_action_table(sg,uc,fname,opts)
        !
        implicit none
        !
        class(am_class_seitz_group), intent(in) :: sg
        class(am_class_unit_cell),   intent(in) :: uc
        character(*),                intent(in) :: fname
        type(am_class_options),      intent(in) :: opts
        real(dp), allocatable :: R(:,:)   ! used for basis functions output
        character(:), allocatable :: xyz(:) ! basis functions output
        integer, allocatable :: PM(:,:)
        integer, allocatable :: P(:,:,:)
        real(dp) :: aa(4)  ! axis and angle 
        integer  :: fid
        integer  :: i,j,k ! loop variables
        character(10) :: buffer
        character(1)  :: s ! sign 
        character(10) :: fmt1
        character(10) :: fmt2
        character(10) :: fmt3
        character(10) :: fmt4
        character(10) :: fmt5
        character(10) :: fmt6
        !
        fmt1 = '(a5)'
        fmt2 = '(i5)'
        fmt3 = '(a10)'
        fmt4 = '(f10.2)'
        fmt5 = '(i10)'
        fmt6 = '(l10)'
        !
        fid = 1
        open(unit=fid,file=trim(fname),status="replace",action='write')
            !
            ! SYMMETRY NUMBER
            write(fid,'(5x)',advance='no')
            write(fid,fmt1,advance='no') '#'
            do i = 1, sg%nsyms
                write(fid,fmt5,advance='no') i
            enddo
            write(fid,*)
            ! CLASS IDENTIFIED
            write(fid,'(5x)',advance='no')
            write(fid,fmt1,advance='no') 'class'
            do i = 1, sg%nsyms
                write(fid,fmt5,advance='no') sg%cc%id(i)
            enddo
            write(fid,*)
            ! HEADER / SEPERATOR
            write(fid,'(5x)',advance='no')
            write(fid,'(5x)',advance='no')
            do i = 1, sg%nsyms
                write(fid,fmt3,advance='no') ' '//repeat('-',9)
            enddo
            write(fid,*)
            ! POINT-SYMMETRY IDENTIFIED
            write(fid,'(5x)',advance='no')
            write(fid,fmt1,advance='no') 'schf.'
            do i = 1, sg%nsyms
                write(fid,fmt3,advance='no') trim(get_ps_name(sg%ps_id(i)))
            enddo
            write(fid,*)
            ! POINT SYMMETRY COMPONENTS
            do i = 1,3
            do j = 1,3
                write(fid,'(5x)',advance='no')
                write(fid,fmt1,advance='no') adjustr('R'//trim(int2char(i))//trim(int2char(j)))
                do k = 1, sg%nsyms
                    write(fid,fmt4,advance='no') sg%seitz_cart(i,j,k)
                enddo
                write(fid,*)
            enddo
            enddo
            ! HEADER / SEPERATOR
            write(fid,'(5x)',advance='no')
            write(fid,fmt1,advance='no') 'ax/th'
            do i = 1, sg%nsyms
                write(fid,fmt3,advance='no') ' '//repeat('-',9)
            enddo
            write(fid,*)
            ! POINT-SYMMETRY ROTATION AXIS
            allocate(character(1)::xyz(3))
            i=0
            i=i+1; xyz(i)='x'
            i=i+1; xyz(i)='y'
            i=i+1; xyz(i)='z'
            do j = 1,3
            write(fid,'(5x)',advance='no')
            write(fid,fmt1,advance='no') 'ax_'//trim(xyz(j))
            do i = 1, sg%nsyms
                aa = rot2axis_angle(R=sg%seitz_cart(1:3,1:3,i))
                write(fid,fmt4,advance='no') aa(j)
            enddo
            write(fid,*)
            enddo
            deallocate(xyz)
            ! POINT-SYMMETRY ROTATION ANGLE
            write(fid,'(5x)',advance='no')
            write(fid,fmt1,advance='no') 'R_th'
            do i = 1, sg%nsyms
                aa = rot2axis_angle(R=sg%seitz_cart(1:3,1:3,i))
                write(fid,fmt4,advance='no') aa(4)*180.0_dp/pi
            enddo
            write(fid,*)
            ! POINT-SYMMETRY ROTOINVERSION VS INVERSION
            write(fid,'(5x)',advance='no')
            write(fid,fmt1,advance='no') 'inv'
            do i = 1, sg%nsyms
                write(fid,fmt6,advance='no') (abs(det(sg%seitz_cart(1:3,1:3,i))+1).lt.tiny)
            enddo
            write(fid,*)
            ! HEADER / SEPERATOR
            write(fid,'(5x)',advance='no')
            write(fid,fmt1,advance='no') 'trans'
            do i = 1, sg%nsyms
                write(fid,fmt3,advance='no') ' '//repeat('-',9)
            enddo
            write(fid,*)
            ! TRANSLATIONAL COMPONENTS (FRAC)
            do j = 1,3
                write(fid,'(5x)',advance='no')
                write(fid,fmt1,advance='no') adjustr('T'//trim(int2char(j)))
                do i = 1, sg%nsyms
                    write(fid,fmt4,advance='no') sg%seitz_cart(j,4,i)
                enddo
                write(fid,*)
            enddo
            ! ACTION chartab (AXES)
            write(fid,'(5x)',advance='no')
            write(fid,fmt1,advance='no') 'axes'
            do i = 1, sg%nsyms
                write(fid,fmt3,advance='no') ' '//repeat('-',9)
            enddo
            write(fid,*)
            allocate(character(1)::xyz(3))
            i=0
            i=i+1; xyz(i)='x'
            i=i+1; xyz(i)='y'
            i=i+1; xyz(i)='z'
            do j = 1,3
                write(fid,'(5x)',advance='no')
                write(fid,fmt1,advance='no') trim(xyz(j))
                do i = 1, sg%nsyms
                    ! inverse here because f(Rr) = R^-1 * f(r)
                    ! for fractional R; its elements will always be an integer...
                    R=inv(sg%seitz_cart(1:3,1:3,i))
                    ! if statement are for the cases in which output needs to be trimmed
                    if (count(abs(R(j,:)).gt.tiny).gt.3) then
                        buffer=' ... ' 
                    else
                        buffer=''
                        do k = 1,3
                        if (abs(R(j,k)).gt.tiny) then
                            ! s is only the sign
                            s = trim(int2char(nint(R(j,k)),'SP'))
                            buffer = trim(buffer)//trim(s)//trim(xyz(k))
                            !
                        endif
                        enddo
                    endif
                    write(fid,fmt3,advance='no') trim(buffer)
                enddo
                write(fid,*)
            enddo
            deallocate(xyz)
            ! ACTION chartab (BASIS FUNCTIONS FOR P ORBITALS)
            write(fid,'(5x)',advance='no')
            write(fid,fmt1,advance='no') 'p orbs.'
            do i = 1, sg%nsyms
                write(fid,fmt3,advance='no') ' '//repeat('-',9)
            enddo
            write(fid,*)
            allocate(character(1)::xyz(3))
            i=0
            i=i+1; xyz(i)='x'
            i=i+1; xyz(i)='y'
            i=i+1; xyz(i)='z'
            do j = 1, size(xyz)
                write(fid,'(5x)',advance='no')
                write(fid,fmt1,advance='no') trim(xyz(j))
                do i = 1, sg%nsyms
                    !
                    R = rot2O3(l=1,R=sg%seitz_cart(1:3,1:3,i))
                    ! if statement are for the cases in which output needs to be trimmed
                    if (count(abs(R(j,:)).gt.tiny).gt.3) then
                        buffer=' ... ' 
                    else
                        buffer=''
                        do k = 1, size(xyz)
                            if (abs(R(j,k)).gt.tiny) then
                                ! s is only the sign
                                s = trim(int2char(nint(R(j,k)),'SP'))
                                buffer = trim(buffer)//trim(s)//trim(xyz(k))
                                !
                            endif
                        enddo
                    endif
                    write(fid,fmt3,advance='no') trim(buffer)
                enddo
                write(fid,*)
            enddo
            deallocate(xyz)
            ! ACTION chartab (BASIS FUNCTIONS FOR D ORBITALS)
            write(fid,'(5x)',advance='no')
            write(fid,fmt1,advance='no') 'd orbs.'
            do i = 1, sg%nsyms
                write(fid,fmt3,advance='no') ' '//repeat('-',9)
            enddo
            write(fid,*)
            allocate(character(2)::xyz(5))
            i=0
            i=i+1; xyz(i)='xx'
            i=i+1; xyz(i)='xy'
            i=i+1; xyz(i)='yy'
            i=i+1; xyz(i)='yz'
            i=i+1; xyz(i)='zz'
            do j = 1,5
                write(fid,'(5x)',advance='no')
                write(fid,fmt1,advance='no') trim(xyz(j))
                do i = 1, sg%nsyms
                    !
                    R = rot2O3(l=2,R=sg%seitz_cart(1:3,1:3,i))
                    !
                    if (count(abs(R(j,:)).gt.tiny).gt.3) then
                        buffer=' ... ' 
                    else
                        buffer=''
                        do k = 1,5
                            if (abs(R(j,k)).gt.tiny) then
                                ! s is only the sign
                                s = trim(int2char(nint(R(j,k)),'SP'))
                                buffer = trim(buffer)//trim(s)//trim(xyz(k))
                                !
                            endif
                        enddo
                    endif
                    write(fid,fmt3,advance='no') trim(buffer)
                enddo
                write(fid,*)
            enddo
            deallocate(xyz)
            ! ACTION chartab (BASIS FUNCTIONS FOR D ORBITALS)
            write(fid,'(5x)',advance='no')
            write(fid,fmt1,advance='no') 'f orbs.'
            do i = 1, sg%nsyms
                write(fid,fmt3,advance='no') ' '//repeat('-',9)
            enddo
            write(fid,*)
            allocate(character(1)::xyz(7))
            i=0
            i=i+1; xyz(i)='1'
            i=i+1; xyz(i)='2'
            i=i+1; xyz(i)='3'
            i=i+1; xyz(i)='4'
            i=i+1; xyz(i)='5'
            i=i+1; xyz(i)='6'
            i=i+1; xyz(i)='7'
            do j = 1,7
                write(fid,'(5x)',advance='no')
                write(fid,fmt1,advance='no') trim(xyz(j))
                do i = 1, sg%nsyms
                    !
                    R = rot2O3(l=3,R=sg%seitz_cart(1:3,1:3,i))
                    !
                    if (count(abs(R(j,:)).gt.tiny).gt.3) then
                        buffer=' ... ' 
                    else
                        buffer=''
                        do k = 1,7
                            if (abs(R(j,k)).gt.tiny) then
                                ! s is only the sign
                                s = trim(int2char(nint(R(j,k)),'SP'))
                                buffer = trim(buffer)//trim(s)//trim(xyz(k))
                                !
                            endif
                        enddo
                    endif
                    write(fid,fmt3,advance='no') trim(buffer)
                enddo
                write(fid,*)
            enddo
            deallocate(xyz)
            ! ACTION chartab (ATOM PERMUTATION)
            select type (sg)
            type is (am_class_space_group)
                !
                P  = permutation_rep(seitz=sg%seitz_frac,tau=uc%tau_frac,flags='',prec=opts%prec)
                PM = permutation_map(P)
                !
                write(fid,'(5x)',advance='no')
                write(fid,fmt1,advance='no') 'atoms'
                do i = 1, sg%nsyms
                    write(fid,fmt3,advance='no') ' '//repeat('-',9)
                enddo
                write(fid,*)
                !
                do j = 1, uc%natoms
                    write(fid,'(5x)',advance='no')
                    write(fid,fmt2,advance='no') j
                    do i = 1, sg%nsyms
                        write(fid,fmt5,advance='no') PM(j,i)
                    enddo
                    write(fid,*)
                enddo
            end select
            !
            write(fid,*)
            write(fid,'(a)') 'p,d,f orbital expansion coefficients are omited for clarity. Also omitted are orbitals which decompose into more than three terms.'
            write(fid,'(a)') 'Cartesian coordinates are used throughout this output.'
            !
        ! close file
        close(fid)
    end subroutine write_action_table

    ! decoding/naming functions 

    function       get_ps_id(R) result(ps_id)
        !
        ! Point symmetries in fractional coordinates so that they are nice integers which can be easily classified.
        !  (should not matter... trace and det are invariant under similarity transforms...)
        !    element:         e    i  c_2  c_3  c_4  c_6  s_2  s_6  s_4  s_3
        !    trace:          +3   -3   -1    0   +1   +2   +1    0   -1   -2
        !    determinant:    +1   -1   +1   +1   +1   +1   -1   -1   -1   -1
        !
        implicit none
        !
        real(dp), intent(in) :: R(:,:,:)
        integer, allocatable :: ps_id(:)
        integer :: nsyms
        integer :: i
        real(dp) :: tr
        real(dp) :: d
        !
        nsyms = size(R,3)
        !
        allocate(ps_id(nsyms))
        ps_id = 0
        !
        do i = 1, nsyms
            ! get trace and determinant (fractional)
            tr = trace(R(1:3,1:3,i))
            d  = det(R(1:3,1:3,i))
            ! The Mathematical Theory of Symmetry in Solids: Representation Theory for
            ! Point Groups and Space Groups. 1 edition. Oxford?: New York: Oxford University
            ! Press, 2010. page 138, chartab 3.8.
            if     ( (abs(tr - 3).lt.tiny) .and. (abs(d - 1).lt.tiny) ) then; ps_id(i) = 1  ! 'e'
            elseif ( (abs(tr + 1).lt.tiny) .and. (abs(d - 1).lt.tiny) ) then; ps_id(i) = 2  ! 'c_2'
            elseif ( (abs(tr - 0).lt.tiny) .and. (abs(d - 1).lt.tiny) ) then; ps_id(i) = 3  ! 'c_3'
            elseif ( (abs(tr - 1).lt.tiny) .and. (abs(d - 1).lt.tiny) ) then; ps_id(i) = 4  ! 'c_4'
            elseif ( (abs(tr - 2).lt.tiny) .and. (abs(d - 1).lt.tiny) ) then; ps_id(i) = 5  ! 'c_6'
            elseif ( (abs(tr + 3).lt.tiny) .and. (abs(d + 1).lt.tiny) ) then; ps_id(i) = 6  ! 'i'
            elseif ( (abs(tr - 1).lt.tiny) .and. (abs(d + 1).lt.tiny) ) then; ps_id(i) = 7  ! 's_2'
            elseif ( (abs(tr - 0).lt.tiny) .and. (abs(d + 1).lt.tiny) ) then; ps_id(i) = 8  ! 's_6'
            elseif ( (abs(tr + 1).lt.tiny) .and. (abs(d + 1).lt.tiny) ) then; ps_id(i) = 9  ! 's_4'
            elseif ( (abs(tr + 2).lt.tiny) .and. (abs(d + 1).lt.tiny) ) then; ps_id(i) = 10 ! 's_3'
            endif
            if (ps_id(i).eq.0) stop 'Unable to identify point symmetry.'
        enddo
    end function   get_ps_id

    function       get_ax_name(ax_id) result(ax_name)
        ! 0 principal rotation axis, +1 parallel-rotation, 2 perp.-rotation, 4 reflection parallel (sv), 5 reflection perpendicular (sh), 6 reflection diagonal (sd)
        implicit none
        !
        integer ,intent(in) :: ax_id
        character(:), allocatable :: ax_name
        !
        allocate(character(10)::ax_name)
        !
        select case(ax_id)
        case(1); ax_name = 'P'
        case(2); ax_name = '='
        case(3); ax_name = 'T'
        case(4); ax_name = 's_v'
        case(5); ax_name = 's_h'
        case(6); ax_name = 's_d'
        case default
            ax_name = ''
        end select
        !
    end function   get_ax_name

    function       get_pg_name(pg_code) result(pg_name)
        ! Hermann-Mauguin     Shubnikov[1]        Schoenflies     Orbifold    Coxeter 
        ! 1                   1\                  C1              11          [ ]+    
        ! -1                  \tilde{2}           S2 (=Ci)        ×           [2+,2+] 
        ! 2                   2\                  C2              22          [2]+    
        ! m                   m\                  C1h             *           [ ]     
        ! 2/m                 2:m\                C2h             2*          [2,2+]  
        ! 222                 2:2\                D2              222         [2,2]+  
        ! mm2                 2 \cdot m\          C2v             *22         [2]     
        ! mmm                 m \cdot 2:m\        D2h             *222        [2,2]   
        ! 4                   4\                  C4              44          [4]+    
        ! -4                  \tilde{4}           S4              2×          [2+,4+] 
        ! 4/m                 4:m\                C4h             4*          [2,4+]  
        ! 422                 4:2\                D4              422         [4,2]+  
        ! 4mm                 4 \cdot m\          C4v             *44         [4]     
        ! 42m                 \tilde{4}\cdot m    D2d             2*2         [2+,4]  
        ! 4/mmm               m \cdot 4:m\        D4h             *422        [4,2]   
        ! 3                   3\                  C3              33          [3]+    
        ! -3                  \tilde{6}           S6 (=C3i)       3×          [2+,6+] 
        ! 32                  3:2\                D3              322         [3,2]+  
        ! 3m                  3 \cdot m\          C3v             *33         [3]     
        ! -3m                 \tilde{6}\cdot m    D3d             2*3         [2+,6]  
        ! 6                   6\                  C6              66          [6]+    
        ! -6                  3:m\                C3h             3*          [2,3+]  
        ! 6/m                 6:m\                C6h             6*          [2,6+]  
        ! 622                 6:2\                D6              622         [6,2]+  
        ! 6mm                 6 \cdot m\          C6v             *66         [6]     
        ! 6m2                 m \cdot 3:m\        D3h             *322        [3,2]   
        ! 6/mmm               m \cdot 6:m\        D6h             *622        [6,2]   
        ! 23                  3/2\                T               332         [3,3]+  
        ! m3                  \tilde{6}/2         Th              3*2         [3+,4]  
        ! 432                 3/4\                O               432         [4,3]+  
        ! 43m                 3/\tilde{4}         Td              *332        [3,3]   
        ! m3m                 \tilde{6}/4         Oh              *432        [4,3]   
        !
        implicit none
        !
        integer, intent(in) :: pg_code
        character(string_length_schoenflies) :: pg_name
        !
        select case (pg_code)
                case(1);  pg_name = 'c_1'
                case(2);  pg_name = 's_2'
                case(3);  pg_name = 'c_2'
                case(4);  pg_name = 'c_1h'
                case(5);  pg_name = 'c_2h'
                case(6);  pg_name = 'd_2'
                case(7);  pg_name = 'c_2v'
                case(8);  pg_name = 'd_2h'
                case(9);  pg_name = 'c_3'
                case(10); pg_name = 's_6'
                case(11); pg_name = 'd_3'
                case(12); pg_name = 'c_3v'
                case(13); pg_name = 'd_3d'
                case(14); pg_name = 'c_4'
                case(15); pg_name = 's_4'
                case(16); pg_name = 'c_4h'
                case(17); pg_name = 'd_4'
                case(18); pg_name = 'c_4v'
                case(19); pg_name = 'd_2d'
                case(20); pg_name = 'd_4h'
                case(21); pg_name = 'c_6'
                case(22); pg_name = 'c_3h'
                case(23); pg_name = 'c_6h'
                case(24); pg_name = 'd_6'
                case(25); pg_name = 'c_6v'
                case(26); pg_name = 'd_3h'
                case(27); pg_name = 'd_6h'
                case(28); pg_name = 't'
                case(29); pg_name = 't_h'
                case(30); pg_name = 'o'
                case(31); pg_name = 't_d'
                case(32); pg_name = 'o_h'
            case default
                call am_print('ERROR','Schoenflies code for point-group unknown.',flags='E')
                call am_print('pg_code',pg_code)
                stop
            end select
    end function   get_pg_name

    subroutine     put_identity_first(seitz)
        !
        implicit none
        !
        real(dp), intent(inout) :: seitz(:,:,:)
        real(dp) :: id(4,4)
        integer  :: nsyms
        integer  :: i
        !
        id = eye(4)
        nsyms = size(seitz,3)
        !
        search : do i = 1, nsyms
            if (isequal(seitz(:,:,i),id)) then
                seitz(:,:,i) = seitz(:,:,1)
                seitz(:,:,1) = id
                exit search
            endif
        enddo search
        ! 
    end subroutine put_identity_first

    function       get_ax_id(R) result(ax_id)
        ! ax_id: 1  principal axis
        !        2  rotation   parallel
        !        3  rotation   perpendicular
        !        4  reflection parallel (sv)
        !        5  reflection perpendicular (sh)
        !        6  reflection diagonal (sd)
        implicit none
        !
        real(dp), intent(in) :: R(:,:,:)
        integer , allocatable :: ax_id(:)
        integer , allocatable :: ps_ns_rank(:)
        real(dp), allocatable :: aa(:,:)
        real(dp), allocatable :: v(:,:) ! null space vectors
        integer , allocatable :: ps_id(:)
        real(dp) :: c(3) ! cross product
        real(dp) :: d    ! dot product
        integer  :: nsyms
        real(dp) :: highsym
        integer  :: highsym_i
        real(dp) :: zero(3)
        real(dp) :: identity(3,3)
        integer  :: i,j,k
        ! get dimensions
        nsyms = size(R,3)
        ! get ps_id
        ps_id = get_ps_id(R=R)
        ! allocate space
        allocate(ax_id(nsyms))
        ax_id = 0
        ! initialize zero
        zero = 0
        ! initialize identity
        identity = eye(3)
        ! identify rotations and reflections
        allocate(ps_ns_rank(nsyms))
        do i = 1, nsyms
            ps_ns_rank(i) = size(null_svd(R(1:3,1:3,i)-identity),2)
            if (ps_ns_rank(i).eq.0) then
            ps_ns_rank(i) =-size(null_svd(R(1:3,1:3,i)+identity),2)
            endif
        enddo
        ! among the rotations find the highest symmetry (principal) axis
        allocate(aa(4,nsyms))
        aa = 0
        do i = 1, nsyms
        if (ps_ns_rank(i).eq.1) then
            aa(1:4,i) = rot2axis_angle(R(1:3,1:3,i))
        endif
        enddo
        ! find highest symmetry axis, if multiple reflctions are the same choose the on alined with Z
        highsym   = 3.15_dp ! slightly larger than pi
        highsym_i = 1
        search : do i = 1, nsyms
        if (ps_ns_rank(i).eq.1) then
            if (abs(aa(4,i)).le.highsym+tiny) then
                highsym   = abs(aa(4,i))
                highsym_i = i
                ! if multiple reflctions are the same choose the on alined with Z
                if (isequal(aa(1:3,i),real([0,0,1],dp))) then
                    highsym = highsym - tiny
                endif
            endif
        endif
        enddo search
        ! set it as the principal axis
        ax_id(highsym_i) = 1
        ! see which rotations are perpendicular to/parallel with the principal symmetry axis
        do i = 1, nsyms
        if (i.ne.highsym_i)     then
        if (ps_ns_rank(i).eq.1) then
            d = dot_product(aa(1:3,i),aa(1:3,highsym_i))
            if     (    abs(d)-0.0_dp .lt.tiny) then
                ! find perpendicular axes
                ax_id(i) = 2
            elseif (abs(abs(d)-1.0_dp).lt.tiny) then
                ! find parallel axes
                ax_id(i) = 3
            endif
        endif
        endif
        enddo
        ! see which planes of reflections contain the principal axis
        do i = 1, nsyms
        if (ps_ns_rank(i).eq.2) then
            v = null_svd(R(1:3,1:3,i)-identity)
            ! define reflection plane as the normal vector
            c = cross_product(v(:,1),v(:,2))
            ! get dot product with highest symmetry axis
            d = dot_product(c,aa(1:3,highsym_i))
            ! note that perp and para are switched below to agree with convetion.
            ! (using the normal vector to determine whether the reflection plane contains or is perpendicular to the principal axis)
            if     (    abs(d)-0.0_dp .lt.tiny) then
                ! If plane contains the principle rotation axis (i.e., parallel), it is a vertical plane (sv)
                ax_id(i) = 4
            elseif (abs(abs(d)-1.0_dp).lt.tiny) then
                ! If plane is perpendicular to the principle rotation axis, it is a horizontal plane (sh)
                ax_id(i) = 5
            endif
            ! if it is still equal to zero, check whether it is a sd ( requires the reflection normal to bisect the angle between two C2 axes)
            if (ax_id(i).eq.0) then
                sd_search : do j = 1, nsyms
                if (ps_ns_rank(j).eq.1) then
                do k = j, nsyms
                if (j.ne.k) then
                if (ps_ns_rank(k).eq.1) then
                    if ( abs(dot_product(aa(1:3,j),c)-dot_product(aa(1:3,k),c)).lt.tiny ) then
                        ax_id(i) = 6
                        exit sd_search
                    endif
                endif
                endif
                enddo
                endif
                enddo sd_search
            endif
        endif
        enddo
    end function   get_ax_id

end module am_symmetry


