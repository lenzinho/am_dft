module am_symmetry_tables

    use am_constants
    use am_matlab
    use am_mkl
    use am_stdout
    use dispmodule

    implicit none

    public

    type, public :: am_class_conjugacy_class
        integer                  :: nclasses                ! nclasses number of classes
        integer    , allocatable :: id(:)                   ! class_id(nsyms) identify conjugacy classes
        integer    , allocatable :: nelements(:)            ! class_nelements(nclasses) get number of elements in each class
        integer    , allocatable :: member(:,:)             ! members(nclass,maxval(class_nelements)) record indicies of each class_member element for each class 
        integer    , allocatable :: representative(:)       ! class_nelements(nclasses) get number of elements in each class
        integer    , allocatable :: class_matrices(:,:,:)   ! class_nelements(nclasses) get number of elements in each class
        contains
    end type am_class_conjugacy_class

    type, public :: am_class_character_table
        integer                  :: nirreps
        complex(dp), allocatable :: chartab(:,:)            ! character table ( irreps * classes )
        character(:),allocatable :: irrep_label(:)          ! irrep label
        integer    , allocatable :: irrep_dim(:)            ! irrep dimension
        complex(dp), allocatable :: irrep_proj(:,:,:)       !
        complex(dp), allocatable :: irrep_proj_V(:,:)       ! eigenvectors of projection
        complex(dp), allocatable :: wigner_proj(:,:,:)      !
    end type am_class_character_table

    type, public :: am_class_multiplication_table
        integer    , allocatable :: multab(:,:)             ! multiplication table
        integer    , allocatable :: inv_id(:)               ! index of inverse element
        integer    , allocatable :: corder(:)               ! cyclic order of the element
        integer    , allocatable :: gen_id(:,:)             ! gen_id(i,:) identifies generators which produce element i
        integer    , allocatable :: gen(:)                  ! complete list of group generators
        integer    , allocatable :: commutator_id(:,:)      ! get_commutator_id(i,:) list of group symmetries which commute with symmetry i 
        integer    , allocatable :: rr(:,:,:)               ! regular representation matrices
    end type am_class_multiplication_table

contains

    ! multiplication table

    function       get_multab(sym,flags) result(multab)
        !
        ! flag options:
        !  seitz   - reduces seitz(1:3,4) to between 0 and 1
        !  prog    - shows progress bar
        !
        implicit none
        !
        real(dp)    , intent(in) :: sym(:,:,:) ! list of 2D reps...
        character(*), intent(in) :: flags
        integer , allocatable :: multab(:,:)
        real(dp), allocatable :: W(:,:) ! workspace
        integer :: ref ! used for checking for closure
        integer :: nsyms, nbases
        integer :: i, j ! loop variabls
        integer :: k, kmax ! for progress bar
        ! get dimensions
        nbases=size(sym,1)
        nsyms=size(sym,3)
        ! check that first element is the identity
        if (.not.isequal(sym(:,:,1),eye(nbases))) stop 'ERROR [get_multab]: idensity not first'
        ! allocate space
        allocate(W(nbases,nbases))
        allocate(multab(nsyms,nsyms))
        !
        k=0
        kmax=nsyms**2
        do i = 1, nsyms
        do j = 1, nsyms
            ! show progress bar
            if (index(flags,'prog').ne.0) then
                k=k+1
                call show_progress(iteration=k, maximum=kmax)
            endif
            ! multiply the two sym operators
            W = matmul(sym(:,:,i),sym(:,:,j))
            ! if 'seitz' flagged, reduce translational part to primitive cell 
            if (index(flags,'seitz')) then
                W(1:3,4) = modulo(W(1:3,4)+tiny,1.0_dp)-tiny
            endif
            ! get matching symmetry
            multab(i,j) = get_matching_element_index_in_list(list=sym,elem=W)
        enddo
        enddo
        ! quick consitency check for closure. not comprehensive
        ref = nsyms*(nsyms+1.0_dp)/2.0_dp
        do i = 1, nsyms
            if (sum(multab(:,i)).ne.ref) stop 'ERROR [get_multab]: multab is incorrect'
        enddo
        !
        contains
        function       get_matching_element_index_in_list(list,elem) result(i)
            !
            implicit none
            real(dp), intent(in) :: list(:,:,:)
            real(dp), intent(in) :: elem(:,:)
            integer :: i
            !
            do i = 1,size(list,3)
                if (all(abs(list(:,:,i)-elem).lt.tiny)) return
            enddo
        end function   get_matching_element_index_in_list
    end function   get_multab

    pure function  get_inverse_indices(multab) result(inv_id)
        ! returns the inv_id of the inverse elements given the multiplication chartab.
        ! requires identity as first element
        implicit none
        !
        integer, intent(in) :: multab(:,:)
        integer, allocatable :: inv_id(:)
        integer, allocatable :: P(:,:)
        !
        inv_id = [1:size(multab,1)]
        !
        allocate(P,source=multab)
        where (P.ne.1) P = 0
        !
        inv_id = matmul(P,inv_id)
        !
    end function   get_inverse_indices

    function       get_cyclic_order(multab) result(corder)
        !
        implicit none
        !
        integer, intent(in) :: multab(:,:)
        integer,allocatable :: corder(:)
        integer :: power
        integer :: i, j
        integer :: nsyms
        !
        nsyms = size(multab,2)
        !
        allocate(corder(nsyms))
        corder = 0
        do i = 1, nsyms
            power = i
            ! apply operation until identity (1) is returned
            do j = 1, 100
                if (power.eq.1) then
                    corder(i) = j
                    exit
                endif
                power = multab(i,power)
            enddo
            if (corder(i).eq.0) stop 'ERROR [get_cyclic_order]: unable to determine symmetry order'
        enddo
        !
    end function   get_cyclic_order

    function       get_regular_rep(multab) result(rr)
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
        integer, intent(in)  :: multab(:,:)
        integer, allocatable :: rr(:,:,:)
        integer, allocatable :: sortmat(:,:)
        integer, allocatable :: inds(:)
        integer, allocatable :: multab_sorted(:,:)
        integer :: nsyms
        integer :: i
        ! get number of symmetries (and size of basis)
        nsyms = size(multab,1)
        ! initializ indices
        allocate(inds, source=[1:nsyms])
        ! initialize sortmat
        allocate(sortmat(nsyms,nsyms))
        where (multab.eq.1)
            sortmat = 1
        elsewhere              
            sortmat = 0
        endwhere
        ! re-order rows so that identities are along diagonal
        allocate(multab_sorted(nsyms,nsyms))
        inds = matmul(sortmat(:,:),inds)
        multab_sorted = multab(inds,:)
        ! construct regular representation from multiplication table
        allocate(rr(nsyms,nsyms,nsyms))
        do i = 1, nsyms
            where (multab_sorted.eq.i) 
                rr(:,:,i) = 1
            elsewhere
                rr(:,:,i) = 0
            endwhere
        enddo
    end function   get_regular_rep

    function       get_generator_id(multab) result(gen_id)
        !
        implicit none
        !
        integer, intent(in) :: multab(:,:) ! multiplication chartab
        integer,allocatable, target :: gen_id(:,:) ! which generators does the symmetry belong to
        integer, pointer :: pntr_i(:)
        integer, pointer :: pntr_j(:)
        integer, target  :: fundamental(2)
        integer,allocatable :: gen_u(:)
        integer :: i, j
        integer :: m, n
        integer :: nsyms
        !
        nsyms = size(multab,2)
        !
        allocate(gen_id(nsyms,nsyms))
        gen_id = 0
        !
        do i = 1, nsyms
            ! loop over elements
        do j = 1, i
            ! check if generator for the element found
            if (gen_id(multab(i,j),1).eq.0) then
                fundamental(1:2) = [i,j]
                ! check if generator for generator found
                if (gen_id(i,1).ne.0) then
                    m = count(gen_id(i,:).ne.0)
                    pntr_i => gen_id(i,1:m)
                else
                    m = 1
                    pntr_i => fundamental(1:1)
                endif
                if (gen_id(i,1).ne.0) then
                    n = count(gen_id(j,:).ne.0)
                    pntr_j => gen_id(j,1:n)
                else
                    n = 1
                    pntr_j => fundamental(2:2)
                endif
                ! 
                gen_u = unique([pntr_i,pntr_j])
                ! save generator
                gen_id(multab(i,j),1:size(gen_u)) = nint(sort(real(gen_u,dp)))
            endif
        enddo
        enddo
        !
    end function   get_generator_id

    function       get_commutator_id(multab) result(commutator_id)
        !
        implicit none
        !
        integer, intent(in) :: multab(:,:) ! multiplication chartab
        integer,allocatable :: commutator_id(:,:) ! symmetries which commute with symmetry i
        integer :: i, j, k
        integer :: nsyms
        !
        nsyms = size(multab,2)
        !
        allocate(commutator_id(nsyms,nsyms))
        commutator_id = 0
        !
        do i = 1, nsyms
            k = 0
            do j = 1, nsyms
            if (multab(i,j).eq.multab(j,i)) then
                k = k + 1
                commutator_id(i,k) = j
            endif
            enddo
        enddo
        !
    end function   get_commutator_id

    ! conjugacy classes

    function       get_class_id(multab,inv_id,ps_id) result(class_id)
        !
        ! for AX = XB, if elements A and B are conjugate pairs for some other element X in the group, then they are in the same class
        !
        implicit none
        !
        integer, intent(in) :: multab(:,:) ! multiplication chartab
        integer, intent(in) :: ps_id(:)    ! sort classes basd on point symmetry id: puts identity first and orders classes containing inversion nicely in the second half
        integer, intent(in) :: inv_id(:)   ! inv_id(nsyms) index of the inverse of symmetry element i
        integer, allocatable :: conjugates(:)
        integer, allocatable :: class_id(:)
        integer :: i, j, k
        integer :: nsyms
        integer :: nclasses
        !
        nsyms = size(multab,2)
        ! allocate space for conjugacy class
        allocate(class_id(nsyms))
        class_id = 0
        ! allocate space for conjugate elements
        allocate(conjugates(nsyms))
        ! determine conjugacy classes
        k = 0
        do i = 1, nsyms
        if (class_id(i).eq.0) then
            k=k+1
            ! conjugate each element with all other group elements
            ! A = X(j) * B * X(j)^-1
            do j = 1, nsyms
                conjugates(j) = multab(j,multab(i,inv_id(j)))
            enddo
            ! for each subgroup element created by conjugation find the corresponding index of the element in the group
            ! in order to save the class class_id number
            do j = 1, nsyms
                class_id( conjugates(j) ) = k
            enddo
        endif
        enddo
        nclasses = k
        ! relabel based on how many elements each class has
        call relabel_based_on_occurances(class_id)
        ! relabel classes based on point symmetry id of representative class element
        call relabel_based_on_ps_id(class_id,ps_id)
        ! make sure the identity is in the first class
        ! class_id(i=1) is the class of the first element (the identity); swap it's location with whaterver elements are in the first class
        where (class_id.eq.1) class_id = class_id(1)
        class_id(1)=1
        ! check that classes are disjoint and complete
        if (any(class_id.eq.0)) stop 'Not every element in the group has been assigned a conjugacy class.'
        !
        contains        
        subroutine     relabel_based_on_occurances(class_id)
            !
            use am_rank_and_sort
            !
            implicit none
            !
            integer , intent(inout) :: class_id(:)
            integer , allocatable :: A_sorted(:)
            integer , allocatable :: occurances(:) 
            integer , allocatable :: reverse_sort(:)
            integer , allocatable :: sorted_indices(:)
            integer , allocatable :: list_relabled(:)
            integer :: nsyms, i, j
            !
            nsyms = size(class_id,1)
            !
            allocate(occurances(nsyms))
            do i = 1, nsyms
                occurances(i) = count(class_id(i).eq.class_id)
            enddo
            !
            ! this quick and dirty procedure lifts degeneracies
            allocate(sorted_indices(nsyms))
            allocate(A_sorted(nsyms))
            call rank((1+maxval(class_id))*occurances+class_id,sorted_indices)
            A_sorted = class_id(sorted_indices)
            !
            allocate(list_relabled(nsyms))
            list_relabled = 0
            !
            j=1
            list_relabled(1) = 1
            do i = 2, nsyms
                if (A_sorted(i).ne.A_sorted(i-1)) then
                    j=j+1
                endif
                list_relabled(i) = j
            enddo
            ! return everything to the original order at call
            allocate(reverse_sort(nsyms))
            reverse_sort(sorted_indices)=[1:nsyms]
            class_id=list_relabled(reverse_sort)
        end subroutine relabel_based_on_occurances
        subroutine     relabel_based_on_ps_id(class_id,ps_id)
            !
            use am_rank_and_sort
            !
            implicit none
            !
            integer, intent(inout) :: class_id(:)
            integer, intent(in)  :: ps_id(:)
            integer :: nclasses
            integer, allocatable :: s(:)
            integer, allocatable :: inds(:)
            integer, allocatable :: rinds(:)
            integer, allocatable :: class_member(:,:)
            integer :: nsyms, i, j
            ! get dimensions
            nsyms = size(class_id)
            nclasses = maxval(class_id)
            ! create inds
            allocate(inds,source=[1:nclasses])
            ! identify class members 
            ! members(nclass,maxval(class_nelements))
            class_member = id_member(class_id)
            ! for each class, identify the representative element
            allocate(s(nclasses))
            do i = 1, nclasses
                s(i) = ps_id( class_member(i,1) )
            enddo
            ! sort classes based on representative element id
            call rank(s,inds)
            allocate(rinds(nclasses))
            rinds(inds) = [1:nclasses] 
            ! at this point inds contains the new labeling for class_id
            ! [1:nclasses] => inds; relabel according to new labeling
            do i = 1, nsyms
                class_id(i) = rinds( class_id(i) )
            enddo
        end subroutine relabel_based_on_ps_id
    end function   get_class_id

    function       get_class_matrices(rr,nclasses,class_nelements,class_member) result(class_matrices)
        !
        implicit none
        !
        integer, intent(in) :: rr(:,:,:)
        integer, intent(in) :: nclasses
        integer, intent(in) :: class_nelements(:)
        integer, intent(in) :: class_member(:,:)
        integer,allocatable :: class_matrices(:,:,:)
        integer :: i,j
        integer :: nbases
        !
        ! get number of things
        nbases  = size(rr,1)
        ! class matrices
        allocate(class_matrices(nbases,nbases,nclasses))
        class_matrices = 0
        do i = 1, nclasses
        do j = 1, class_nelements(i)
            class_matrices(:,:,i) = class_matrices(:,:,i) + rr(:,:,class_member(i,j))
        enddo
        enddo
        ! check that their sum equals ones() matrix
        if (.not.isequal( sum(class_matrices,3), nint(ones(nbases)))) then
            stop 'ERROR [get_class_matrices]: class matrices do not sum to ones matrix'
        endif
    end function   get_class_matrices

    ! character table

    function       get_chartab(multab,nclasses,class_nelements,class_member,ps_id) result(chartab)
        !
        use am_rank_and_sort
        !
        implicit none
        !
        integer, intent(in) :: multab(:,:)
        integer, intent(in) :: ps_id(:)
        integer, intent(in) :: nclasses
        integer, intent(in) :: class_nelements(:) ! number of elements in each class
        integer, intent(in) :: class_member(:,:) ! members(nclass,maxval(class_nelements))
        complex(dp), allocatable :: chartab(:,:) ! class constant chartab which becomes character chartab (can be complex!)
        real(dp)   , allocatable :: Re(:,:)
        real(dp)   , allocatable :: Im(:,:)
        integer :: i, j, k
        integer :: nsyms
        integer :: nirreps  ! just to make things clearer; always equal to nclasses
        integer, allocatable :: H(:,:,:) ! class multiplication 
        integer, allocatable :: irrep_dim(:) ! sym dimensions
        integer, allocatable :: indices(:)
        real(dp) :: wrk
        !
        nsyms = size(multab,2)
        ! get class coefficients (thenumber of times class k appears in the pdocut o class j and k )
        ! Wooten p 35, Eq 2.10; p 40, Example 2.8; p 89, Example 4.6
        allocate(H(nclasses,nclasses,nclasses))
        do i = 1,nclasses
        do j = 1,nclasses
        do k = 1,nclasses
            H(j,k,i) = count( pack( multab(class_member(i,1:class_nelements(i)),class_member(j,1:class_nelements(j))) , .true. ) .eq. class_member(k,1) )
        enddo
        enddo
        enddo
        ! chartab(nirreps,nclasses) determine eigenvector which simultaneously diagonalizes i Hjk(:,:) matrices
        chartab = eigenanalysis(cmplx(H,0,dp))
        ! get dimensions of representations using Eq. 4.53, p 91 of Wooten
        ! Note: The absolute square is the result of the conjugate multiplication
        nirreps = nclasses
        allocate(irrep_dim(nirreps))
        do i = 1, nirreps
            ! sum on classes
            wrk = 0
            do j = 1, nclasses
                wrk = wrk + abs(chartab(i,j))**2/real(class_nelements(j),dp)
            enddo
            ! compute irrep dimension
            irrep_dim(i) = nint((nsyms/wrk)**0.5_dp)
        enddo
        ! convert class constant chartab into character chartab
        do i = 1, nirreps
        do j = 1, nclasses
            chartab(i,j) = irrep_dim(i)/real(class_nelements(j),dp)*chartab(i,j)
        enddo
        enddo
        ! correct basic rounding error
        allocate(Re, source= real(chartab))
        allocate(Im, source=aimag(chartab))
        where (abs(nint(Re)-Re).lt.tiny) Re = nint(Re)
        where (abs(nint(Im)-Im).lt.tiny) Im = nint(Im)
        chartab = cmplx(Re,Im,dp)
        ! initialize indices sort sorting
        allocate(indices(nirreps))
        ! sort irreps based on character of class containing identity
        ! find class containing inversion (1 is the ps_id for identity)
        id_search : do i = 1, nclasses
        do j = 1, class_nelements(i)
            if (ps_id(class_member(i,j)).eq.1) then
                ! sort based on inversion-containing class character
                call rank(nint(real(chartab(:,i))),indices)
                ! indices = indices(nirreps:1:-1)
                chartab = chartab(indices,:)
                ! exit loop
                exit id_search
            endif
        enddo
        enddo id_search
        ! sort irreps based on character of class containing inversion
        ! find class containing inversion (6 is the ps_id for inversion)
        inv_search : do i = 1, nclasses
        do j = 1, class_nelements(i)
            if (ps_id(class_member(i,j)).eq.6) then
                ! sort based on inversion-containing class character
                call rank(-sign(1,nint(real(chartab(:,i)))),indices)
                ! indices = indices(nirreps:1:-1)
                chartab = chartab(indices,:)
                ! exit loop
                exit inv_search
            endif
        enddo
        enddo inv_search
        ! check that the number of elements ri in class Ci is a divisor of theorder of the group
        ! Symmetry and Condensed Matter Physics: A Computational Approach. 1 edition. Cambridge, UK?; New York: Cambridge University Press, 2008. p 35.
        do i = 1, nclasses
            if (modulo(nsyms,class_nelements(i)).ne.0) then
                stop 'Number of elements in class is not a divisor of the group order.'
            endif
        enddo
        !
        contains
        function       eigenanalysis(A) result(D)
            ! determines characters which simultaneously diagonalizes n square matrices A(:,:,i)
            ! only possible if A(:,:,i) commute with each other.
            implicit none
            !
            complex(dp), intent(in)  :: A(:,:,:)
            complex(dp), allocatable :: V(:,:) ! these two should match in this particular situation (character chartab determination from class multiplication )
            complex(dp), allocatable :: D(:,:)
            integer :: m, n, i
            !
            ! get dimensions
            n = size(A,1)
            m = size(A,3)
            ! get similarity transform which simultaneosuly diagolnaizes matrices A
            V = eigenspace(A)
            ! check that matrices have been diagonalized properly and save eigenvalues
            allocate(D(n,n))
            do i = 1, m
                D(:,i) = diag( matmul(inv(V),matmul(A(:,:,i),V)) )
            enddo
            !
            ! should enforce check here...
            !
            ! ! normalize each eigencolumn to the element of the column (that is, the first element, i.e. row, corresponds to the class containing the identity element)
            ! do i = 1, m
            !     V(:,i)=V(:,i)
            ! enddo
            ! V = transpose(V)
            !
            ! at this point V and D should be identical. The elements are to within +/- sign.
            ! D is more robust 
            !
        end function   eigenanalysis
    end function   get_chartab

    function       get_muliken_label(chartab,nclasses,class_nelements,class_member,ps_id) result(irrep_label)
        !
        implicit none
        !
        complex(dp), intent(in) :: chartab(:,:) ! class constant chartab which becomes character chartab (can be complex!)
        integer, intent(in) :: nclasses
        integer, intent(in) :: class_nelements(:)
        integer, intent(in) :: class_member(:,:)
        integer, intent(in) :: ps_id(:)
        character(:), allocatable :: irrep_label(:)
        integer :: i, j, k
        integer :: nirreps
        integer :: class_containing_identity
        integer :: class_containing_inversion
        integer :: class_containing_highsym
        integer, allocatable :: highsymlist(:)
        !
        ! they are always identical...
        nirreps = nclasses
        !
        ! create irrep label
        allocate(character(7) :: irrep_label(nirreps))
        !
        class_containing_identity  = find_class_containing_ps(nclasses=nclasses, id=1, ps_id=ps_id, class_nelements=class_nelements, class_member=class_member)
        class_containing_inversion = find_class_containing_ps(nclasses=nclasses, id=6, ps_id=ps_id, class_nelements=class_nelements, class_member=class_member)
        allocate(highsymlist, source=[5,10,4,9,3,8,2,7,1,6])
        do i = 1,10
            class_containing_highsym = find_class_containing_ps(nclasses=nclasses, id=highsymlist(i), ps_id=ps_id, class_nelements=class_nelements, class_member=class_member)
            if (class_containing_highsym.ne.0) exit
        enddo
        !
        do j = 1, nirreps
            ! 0) nullify irrep label
            irrep_label(j) = ''
            ! 1) find class containing identity (ps_id=1)
            i = class_containing_identity
            select case (nint(real(chartab(j,i))))
                case (1)
                    k = class_containing_highsym
                    if (real(chartab(j,k)).ge.0) then 
                        irrep_label(j) = trim(irrep_label(j))//'A'
                    else
                        irrep_label(j) = trim(irrep_label(j))//'B'
                    endif
                case (2); irrep_label(j) = trim(irrep_label(j))//'E'
                case (3); irrep_label(j) = trim(irrep_label(j))//'T'
            end select
            ! 2) find class containing inversion (ps_id=6)
            i = class_containing_inversion
            if (i.ne.0) then
            select case (sign(1,nint(real(chartab(j,i)))))
                case (+1); irrep_label(j) = trim(irrep_label(j))//'_g'
                case (-1); irrep_label(j) = trim(irrep_label(j))//'_u'
            end select
            endif
            !
        enddo
        !
        contains
            pure function find_class_containing_ps(nclasses, id, ps_id, class_nelements, class_member) result(i)
                !
                implicit none
                !
                integer, intent(in) :: id
                integer, intent(in) :: nclasses
                integer, intent(in) :: ps_id(:)
                integer, intent(in) :: class_nelements(:)
                integer, intent(in) :: class_member(:,:)
                integer :: i, j
                !
                do i = 1, nclasses
                do j = 1, class_nelements(i)
                    if (ps_id(class_member(i,j)).eq.id) then
                        return
                    endif
                enddo
                enddo
                ! if not found
                i = 0
            end function  find_class_containing_ps
    end function   get_muliken_label

    subroutine     print_chartab(chartab,nclasses,class_nelements,class_member,irrep_label,ps_id,rep_chi,rep_label)
        !
        implicit none
        !
        complex(dp) , intent(in) :: chartab(:,:) ! class constant chartab which becomes character chartab (can be complex!)
        integer     , intent(in) :: nclasses
        integer     , intent(in) :: class_nelements(:)
        integer     , intent(in) :: class_member(:,:)
        character(*), intent(in) :: irrep_label(:)
        integer     , intent(in) :: ps_id(:)
        complex(dp) , intent(in), optional :: rep_chi(:,:)
        character(*), intent(in), optional :: rep_label(:)
        integer :: i, j, k
        integer :: nirreps
         ! rep to irrep decompostion (used only if rep_chi nd rep_label are present)
        integer :: nreps
        integer :: nsyms
        integer, allocatable :: beta(:,:)
        complex(dp) :: try
        ! complex to symbol
        complex(dp), allocatable :: s(:)
        character(:), allocatable :: str
        integer :: char_start
        ! print formats
        character(50) :: fmts(5)
        !
        ! they are always identical...
        nirreps = nclasses
        !
        ! left-side headers (sum to 10)
        fmts(2) = '(5x,a10)'
        fmts(4) = '(5x,a6,i4)'
        ! chartab (sum to 6)
        fmts(1) = '(i6)'
        fmts(3) = '(a6)'
        fmts(5) = '(f6.2)'
        !
        call print_chartab_header(nclasses,class_nelements,class_member,ps_id,fmts)
        !
        allocate(s(100)) ! value of symbolic output, 30 possible options
        allocate(character(4)::str)
        char_start = 96 ! 97 = a
        k=0
        !
        do i = 1, nirreps
            write(*,'(5x,i2,a8)',advance='no') i, irrep_label(i)
            do j = 1, nclasses
                call print_chartab_symb(Z=chartab(i,j),char_start=char_start,s=s,k=k,fmts=fmts,str=str)
                write(*,fmts(3),advance='no') str
            enddo
            write(*,*)
        enddo
        ! 
        if (present(rep_chi).and.present(rep_label)) then
            if (size(rep_chi,2).ne.nclasses) stop 'dimensions of rep character does not match number of classes'
            !
            ! rep characters
            call print_bar(fmts=fmts, nclasses=nclasses)
            nreps = size(rep_chi,1)
            do i = 1, nreps
                !
                write(*,'(5x,i2,a8)',advance='no') i, rep_label(i)
                do j = 1, nclasses
                    call print_chartab_symb(Z=rep_chi(i,j),char_start=char_start,s=s,k=k,fmts=fmts,str=str)
                    write(*,fmts(3),advance='no') str
                enddo
                write(*,*)
            enddo
        endif
        !
        ! rep decompositions
        if (present(rep_chi).and.present(rep_label)) then
            call print_bar(fmts=fmts, nclasses=nclasses)
            !
            write(*,'(5x,a)') 'Decompositions:'
            ! initialize
            allocate(beta(nreps,nirreps))
            beta = 0
            nsyms = sum(class_nelements)
            ! calculate
            do i = 1, nreps
            do j = 1, nirreps
                try = dot(rep_chi(i,:),chartab(j,:)*class_nelements)/real(nsyms,dp)
                if (abs(aimag(try)).gt.tiny) then
                    stop 'irrep decomposition coefficient is not real'
                endif
                if (abs(nint(real(try))-real(try)).gt.tiny) then
                    stop 'irrep decomposition coefficient is not an integer'
                endif
                beta(i,j) = nint(real(try))
            enddo
            enddo
            ! print
            do i = 1, nreps
            if (any(beta(i,:).ne.0)) then
                write(*,'(5x,a9,a)',advance='no') trim(rep_label(i)), ' = '
                do j = 1, nirreps
                if (beta(i,j).ne.0) then
                    write(*,'(a)',advance='no') trim(int2char(beta(i,j),'SP'))//'*'//trim(irrep_label(j))//'('//trim(int2char(j))//') '
                endif
                enddo
                write(*,*)
            endif
            enddo
        endif
        !
        !
        !
        !
        !
        !
        if (k.ne.0) then
            write(*,'(5x,a)') 'Definitions:'
            do i = 1, k
                write(*,'(5x)',advance='no')
                write(*,'(a)' ,advance='no') char(char_start+k)//' = '
                write(*,'(a)' ,advance='no') trim(dbl2char(real(s(i)),8))//trim(dbl2charSP(aimag(s(i)),9))//'i'
                write(*,*)
            enddo
            !
        endif
        !
        contains
        subroutine     print_chartab_header(nclasses,class_nelements,class_member,ps_id,fmts)
            !
            implicit none
            !
            integer     , intent(in) :: nclasses
            integer     , intent(in) :: class_nelements(:)
            integer     , intent(in) :: class_member(:,:)
            integer     , intent(in) :: ps_id(:)
            character(*), intent(in) :: fmts(:)
            integer :: i
            !
            ! start printing
            write(*,fmts(2),advance='no') 'class'
            do i = 1, nclasses
                write(*,fmts(1),advance='no') i
            enddo
            write(*,*)
            !
            write(*,fmts(2),advance='no') 'elements'
            do i = 1, nclasses
                write(*,fmts(1),advance='no') class_nelements(i)
            enddo
            write(*,*)
            !
            write(*,fmts(2),advance='no') 'repr.'
            do i = 1, nclasses
                write(*,fmts(3),advance='no') trim(decode_pointsymmetry(ps_id(class_member(i,1))))
            enddo
            write(*,*)
            !
            call print_bar(fmts=fmts, nclasses=nclasses)
            !
        end subroutine print_chartab_header
        subroutine     print_chartab_symb(Z,char_start,s,k,fmts,str)
            !
            implicit none
            !
            complex(dp), intent(in) :: Z
            integer    , intent(in) :: char_start
            complex(dp), intent(inout) :: s(:)
            integer    , intent(inout) :: k
            character(*), intent(in) :: fmts(:)
            character(:), allocatable, intent(out) :: str
            !
            integer :: kk, k_exp
            complex(dp) :: s_exp
            real(dp) :: Zr,Zi
            logical :: strmatch
            !
            strmatch = .false.
            !
            Zr= real(Z)
            Zi= aimag(Z)
            !
            if ( isint(Zr) .and. iszero(Zi) ) then
                ! no imaginary, integer real
                strmatch = .true.
                str = trim(int2char(nint(Zr)))
                !
            elseif ( iszero(Zr) .and. isint(Zi) ) then
                ! no real, imaginary integer
                strmatch = .true.
                str = trim(int2char(nint(Zi)))//'i'
                !
            elseif ( (.not. isint(Zr)) .and. (.not. isint(Zi)) ) then
                ! complex number
                !
                if (k.ge.1) then
                search : do kk = 1, k
                    !
                    s_exp = cmplx(0.0_dp,0.0_dp)
                    k_exp = 0 
                    do while ( .not. isequal(s_exp,cmplx(1,0,dp)) )
                        ! do a full loop. complex numbers form a cyclic abelian group.
                        ! exponentiate it until it loops back to one, the identity
                        k_exp = k_exp+1
                        s_exp = s(k)**k_exp
                        ! check positive
                        if ( isequal(Z,s_exp) ) then
                            !
                            strmatch = .true.
                            if (k_exp.eq.1) then
                                str = char(char_start+k)
                            else
                                str = char(char_start+k)//trim(int2char(k_exp))
                            endif
                            exit search
                            !
                        endif
                        ! check negative
                        if ( isequal(Z,-s_exp) ) then
                            !
                            strmatch = .true.
                            if (k_exp.eq.1) then
                                str = '-'//char(char_start+k)
                            else
                                str = '-'//char(char_start+k)//trim(int2char(k_exp))
                            endif
                            exit search
                            !
                        endif
                    enddo
                enddo search
                endif
                !
                ! if match is not found assign a new character to variable
                if (.not.strmatch) then
                    !
                    strmatch = .true.
                    k = k + 1
                    s(k) = Z
                    str = char(char_start+k)
                    !
                endif
            endif
            !
            if (.not.strmatch) then
                write(str,fmts(5)) real(Z)
            endif
            !
        end subroutine print_chartab_symb
        subroutine     print_bar(fmts,nclasses)
            !
            implicit none
            !
            integer     , intent(in) :: nclasses
            character(*), intent(in) :: fmts(:)
            integer :: i
            !
            write(*,fmts(2),advance='no') repeat('-',10)
            do i = 1, nclasses
                write(*,fmts(3),advance='no') repeat('-',6)
            enddo
            write(*,*)
        end subroutine print_bar
    end subroutine print_chartab

    function       get_irrep_dimension(chartab,nclasses,class_nelements,class_member,ps_id) result(irrep_dim)
        !
        implicit none
        !
        complex(dp), intent(in) :: chartab(:,:)
        integer    , intent(in) :: nclasses
        integer    , intent(in) :: class_nelements(:)
        integer    , intent(in) :: class_member(:,:)
        integer    , intent(in) :: ps_id(:)
        integer    , allocatable :: irrep_dim(:)
        integer :: i ! class
        integer :: j ! irrep
        integer :: nirreps
        ! find class containing identity symmetry
        id_search : do i = 1, nclasses
        do j = 1, class_nelements(i)
            if (ps_id(class_member(i,j)).eq.1) then
                exit id_search
            endif
        enddo
        enddo id_search
        ! get dimensions
        nirreps = nclasses
        allocate(irrep_dim(nirreps))
        do j = 1, nirreps
            irrep_dim(j) = nint(real(chartab(j,i)))
        enddo
    end function   get_irrep_dimension

    ! regular representation

    function       get_irrep_proj(rr,chartab,irrep_dim,class_matrices) result(irrep_proj)
        !
        implicit none
        !
        integer    , intent(in) :: rr(:,:,:)
        complex(dp), intent(in) :: chartab(:,:)
        integer    , intent(in) :: irrep_dim(:)
        integer    , intent(in) :: class_matrices(:,:,:)
        complex(dp),allocatable :: irrep_proj(:,:,:)
        integer :: nirreps, nsyms, nbases, nclasses
        integer :: i,j
        !
        ! get number of things
        nbases  = size(rr,1)
        nsyms   = nbases
        nclasses= size(class_matrices,3)
        nirreps = nclasses
        ! allocate space for irrep projection operator
        allocate(irrep_proj(nbases,nbases,nirreps))
        irrep_proj = 0
        ! loop over irreos
        do j = 1, nirreps
            ! loop over symmetries (by looping over class and members)
            do i = 1, nclasses
                irrep_proj(:,:,j) = irrep_proj(:,:,j) + class_matrices(:,:,i) * chartab(j,i)
            enddo
            irrep_proj(:,:,j) = irrep_proj(:,:,j) * irrep_dim(j)/real(nsyms,dp)
        enddo
        ! check that each irrep projection operator is hermitian
        do j = 1, nirreps
            if (.not.isequal(adjoint(irrep_proj(:,:,j)),irrep_proj(:,:,j))) then
                stop 'ERROR [get_irrep_projection]: irrep_proj is not hermitian'
            endif
        enddo
        ! check that each irrep projection operator is idempotent (A^2 = A)
        do j = 1, nirreps
            if (.not.isequal(matmul(irrep_proj(:,:,j),irrep_proj(:,:,j)),irrep_proj(:,:,j))) then
                stop 'ERROR [get_irrep_projection]: irrep_proj is not idempotent'
            endif
        enddo
    end function   get_irrep_proj

    function       get_irrep_proj_V(irrep_proj,irrep_dim) result(irrep_proj_V)
        !
        implicit none
        !
        complex(dp), intent(in) :: irrep_proj(:,:,:) ! can be complex because of complex characters
        integer    , intent(in) :: irrep_dim(:)
        complex(dp),allocatable :: irrep_proj_V(:,:)
        integer,    allocatable :: S_sqr(:)
        integer,    allocatable :: E_sqr(:)
        integer,    allocatable :: inds(:)
        complex(dp),allocatable :: D(:)
        complex(dp),allocatable :: V(:,:)
        logical,    allocatable :: mask(:)
        integer :: i, n, nirreps, nbases
        ! get number of irreps
        nbases  = size(irrep_proj,1)
        nirreps = size(irrep_proj,3)
        ! allocate/initialize stuff
        allocate(inds(nbases))
        inds = [1:nbases]
        allocate(mask(nbases))
        allocate(D(nbases))
        allocate(V(nbases,nbases))
        allocate(S_sqr(nirreps))
        allocate(E_sqr(nirreps))
        allocate(irrep_proj_V(nbases,nbases))
        !
        n = 0
        do i = 1, nirreps
            ! get eigenvectors/eigenvalues
            call am_zgeev(A=irrep_proj(:,:,i), VR=V, D=D)
            ! create mask to exclude nullspace
            where (abs(D).gt.tiny)
                mask = .true.
            elsewhere
                mask = .false.
            endwhere
            ! save indics
            n = n + 1
            S_sqr(i) = n
            n = n + count(mask) - 1
            E_sqr(i) = n
            ! save projection
            irrep_proj_V(:,S_sqr(i):E_sqr(i)) = orth_svd( V(:,pack(inds,mask)) )
        enddo
        ! check
        do i = 1, nirreps
        if ((E_sqr(i)-S_sqr(i)+1).ne.irrep_dim(i)**2) stop 'ERROR [get_irrep_proj_V]: eigenspace mismatch'
        enddo
        !
    end function   get_irrep_proj_V

    function       get_irrep_diag(rr,chartab,irrep_proj,irrep_dim,class_member) result(irrep_diag)
        !
        implicit none
        !
        integer    , intent(in) :: rr(:,:,:)
        complex(dp), intent(in) :: chartab(:,:)
        complex(dp), intent(in) :: irrep_proj(:,:,:) ! projection
        integer    , intent(in) :: irrep_dim(:)
        integer    , intent(in) :: class_member(:,:)
        complex(dp),allocatable :: irrep_diag(:,:)
        complex(dp),allocatable :: phi(:,:)
        complex(dp),allocatable :: wrk(:,:)
        complex(dp),allocatable :: tr(:,:)
        complex(dp),allocatable :: Q(:)
        integer    ,allocatable :: S(:)
        integer    ,allocatable :: E(:)
        integer :: nirreps, nclasses, nbases, nsyms
        integer :: i, j, k, n
        ! get dimensions
        nsyms   = size(rr,3)
        nbases  = size(rr,1)
        nirreps = size(irrep_proj,3)
        nclasses= nirreps
        ! allocate phi
        allocate(phi(nbases,sum(irrep_dim)))
        phi = 0
        ! allocate local character table
        allocate(tr(nirreps,nsyms))
        ! allocate irrep_diag
        allocate(irrep_diag(sum(irrep_dim),nsyms))
        irrep_diag = 0
        ! allocate start/end
        allocate(S(nirreps))
        allocate(E(nirreps))
        ! allocate workspace
        allocate(wrk(nbases,nsyms))
        ! allocate random vector
        allocate(Q(nbases))
        ! Q = sqrt(real([1:nbases],dp))
        Q = pack( rand(nbases,1), .true.)
        ! 
        n = 0
        do j = 1, nirreps
            ! start/end
            n = n + 1
            S(j) = n
            n = n + irrep_dim(j) - 1 
            E(j) = n
            ! project a random vector onto the j-th irreducible subspace
            wrk(:,1) = matmul(irrep_proj(:,:,j), Q )
            ! generate partner symmetry functions
            if (irrep_dim(j).ge.2) then
                k = 2
                do i = 2, nsyms
                    wrk(:,k) = matmul(rr(:,:,i), wrk(:,1))
                    if (k.eq.rank_svd(wrk(:,1:k))) then
                    if (k.eq.irrep_dim(j)) then
                        exit
                    else
                        k = k + 1
                    endif
                    endif
                enddo
            endif
            ! save phi
            phi(:,S(j):E(j)) = orth_svd( wrk(:,1:irrep_dim(j)) )
        enddo
        ! get diagonal symmetry elements
        do i = 1, nsyms
        do j = 1, nirreps
            do k = 1, irrep_dim(j)
                irrep_diag(S(j)+k-1,i) = dot_product(matmul( conjg(phi(:,S(j)+k-1)), rr(:,:,i)), phi(:,S(j)+k-1))
            enddo
            tr(j,i) = sum( irrep_diag(S(j):E(j),i) )
        enddo
        enddo
        ! check that trace matches value in character table
        do j = 1, nirreps
        do i = 1, nclasses
            if ( abs(chartab(j,i)-tr(j,class_member(i,1))) .gt. tiny ) then
                stop 'ERROR [get_irrep_diag]: character table mismatch'
            endif
        enddo
        enddo
    end function   get_irrep_diag

    function       get_wigner_proj(rr,irrep_diag,irrep_dim,class_matrices) result(wigner_proj)
        !
        implicit none
        !
        integer    , intent(in) :: rr(:,:,:)
        complex(dp), intent(in) :: irrep_diag(:,:)
        integer    , intent(in) :: irrep_dim(:)
        integer    , intent(in) :: class_matrices(:,:,:)
        complex(dp),allocatable :: wigner_proj(:,:,:,:) ! wigner_proj(nbases,nbases, maxval(irrep_dim), nirreps )
        integer :: nirreps, nsyms, nbases, nclasses
        integer :: i,j,k,n
        !
        ! get number of things
        nbases  = size(rr,1)
        nsyms   = nbases
        nclasses= size(class_matrices,3)
        nirreps = nclasses
        ! allocate space for wigner projection operator
        allocate(wigner_proj(nbases,nbases,maxval(irrep_dim),nirreps))
        wigner_proj = 0
        ! loop over irreps
        n = 0
        do j = 1, nirreps
            ! loop over dimensions of irreps
            do k = 1, irrep_dim(j)
                ! n is a compound irrep and irrep dim index
                n = n + 1
                ! wooten eq 6.24: sum over symmetries
                do i = 1, nsyms
                wigner_proj(:,:,k,j) = wigner_proj(:,:,k,j) + rr(:,:,i) * irrep_diag(n,i)
                enddo
                wigner_proj(:,:,k,j) = wigner_proj(:,:,k,j) * irrep_dim(j)/real(nsyms,dp)
            enddo
        enddo
        ! check that each wigner projection operator is hermitian
        do j = 1, nirreps
        do k = 1, irrep_dim(j)
            if (.not.isequal(adjoint(wigner_proj(:,:,k,j)),wigner_proj(:,:,k,j))) then
                stop 'ERROR [get_irrep_projection]: irrep_proj is not hermitian'
            endif
        enddo
        enddo
        ! check that each wigner projection operator is idempotent (A^2 = A)
        do j = 1, nirreps
        do k = 1, irrep_dim(j)
        if (.not.isequal(matmul(wigner_proj(:,:,k,j),wigner_proj(:,:,k,j)),wigner_proj(:,:,k,j))) then
            call disp(X=wigner_proj(:,:,k,j),title=tostring(k)//tostring(j),style='underline')
            stop 'ERROR [get_irrep_projection]: irrep_proj is not idempotent'
        endif
        enddo
        enddo
        !
    end function   get_wigner_proj

    function       get_block_transform(rr,irrep_dim,wigner_proj) result(phi)
        !
        implicit none
        !
        integer    , intent(in) :: rr(:,:,:)
        integer    , intent(in) :: irrep_dim(:)
        complex(dp), intent(in) :: wigner_proj(:,:,:,:) ! wigner_proj(nbases,nbases, maxval(irrep_dim), nirreps )
        complex(dp),allocatable :: phi(:,:)
        complex(dp),allocatable :: rr_block(:,:,:)
        complex(dp),allocatable :: Q(:)
        integer :: nirreps, nbases, nsyms
        integer :: sumdim
        integer :: j,k,n, S,E
        ! get dimensions
        nsyms   = size(rr,3)
        nbases  = size(rr,1)
        nirreps = size(wigner_proj,4)
        sumdim  = sum(irrep_dim)
        ! allocate phi
        allocate(phi(nbases,sumdim))
        phi = 0
        ! get a random vector
        allocate(Q(nbases))
        Q = pack( rand(nbases,1), .true. )
        !wigner_proj
        n = 0
        do j = 1, nirreps
            n = n + 1
            S = n
            n = n + irrep_dim(j) - 1
            E = n
            do k = 1, irrep_dim(j)
                phi(:,(S-1)+k) = matmul(wigner_proj(:,:,k,j),Q)
                ! call am_zgeev(A=wigner_proj(:,:,k,j), VR=V, D=Q)
                ! call disp(X=Q,orient='row')
                ! phi(:,S-1+k) = V(:,1)
            enddo
            phi(:,S:E) = orth_svd(phi(:,S:E))
            ! phi(:,S:E) = phi(:,S:E)
        enddo
        !
        allocate(rr_block(sumdim,sumdim,nsyms))
        do j = 1, nsyms
            rr_block(:,:,j) = matmul(matmul(adjoint(phi),rr(:,:,j)), phi)
        enddo
        !
        call disp(rr_block(:,:,5))

        call dump(A=rr_block, fname=trim(outfile_dir_sym)//'/debug'//'/outfile.rr_block')
        call dump(A=phi     , fname=trim(outfile_dir_sym)//'/debug'//'/outfile.phi')
        call dump(A=wigner_proj(:,:,:,5), fname=trim(outfile_dir_sym)//'/debug'//'/outfile.wigner_proj')





        stop
        !
    end function   get_block_transform

    function       get_cycle_structure_str(rr,order) result(cycle_structure_str)
        !
        implicit none
        !
        integer, intent(in) :: rr(:,:,:)
        integer, intent(in) :: order(:)
        character(300), allocatable :: cycle_structure_str(:)
        integer,allocatable :: R(:,:)
        integer,allocatable :: P(:,:) ! power of rr
        logical,allocatable :: T(:)
        integer,allocatable :: inds(:)
        integer :: nsyms
        integer :: nbases
        integer :: i,j 
        !
        nbases = size(rr,1)
        nsyms = size(rr,3)
        !
        allocate(cycle_structure_str(nsyms))
        allocate(R(nbases,nbases))
        allocate(P(nbases,nbases))
        allocate(T(nbases))
        allocate(inds(nbases))
        inds = [1:nbases]
        !
        do i = 1, nsyms
            R = 0
            P = nint(eye(nbases))
            do j = 1, order(i)
                P = matmul(rr(:,:,i),P)
                R = R + P
            enddo
            ! reduce to row echelon form
            P = nint(rref(real(R,dp)))
            ! initialize cycle strucutre string
            cycle_structure_str(i)=' '
            ! write cycle structure
            do j = 1, nbases
                where (P(j,:).eq.1) 
                    T = .true.
                elsewhere
                    T = .false.
                endwhere
                if (count(T).ge.1) then
                    cycle_structure_str(i) = trim(cycle_structure_str(i))//' ('//tostring(pack(inds,T))//')'
                endif
            enddo
        enddo
        !
    end function   get_cycle_structure_str

!     function       get_block_similarity_transform(rr,irrep_dim,irrep_proj,nclasses,commutator_id,class_member) result(phi)
!         !
!         use am_rank_and_sort
!         !
!         implicit none
!         !
!         integer    , intent(in) :: rr(:,:,:)
!         integer    , intent(in) :: irrep_dim(:)
!         complex(dp), intent(in) :: irrep_proj(:,:,:) ! can be complex because of complex characters
!         integer    , intent(in) :: nclasses
!         integer    , intent(in) :: class_member(:,:)
!         integer    , intent(in) :: commutator_id(:,:) ! symmetries which commute with symmetry i
!         complex(dp),allocatable :: phi(:,:)
!         integer    ,allocatable :: S(:)
!         integer    ,allocatable :: E(:) 
!         complex(dp),allocatable :: V(:,:,:) ! eigenvectors
!         complex(dp),allocatable :: D(:,:) ! eigenvalues
!         complex(dp),allocatable :: zeros(:)
!         complex(dp),allocatable :: try(:)
!         complex(dp),allocatable :: V_least(:,:)
!         integer    ,allocatable :: commutator_least(:)
!         ! complex(dp),allocatable :: V_proj(:,:)
!         ! integer,    allocatable :: S_sqr(:)
!         ! integer,    allocatable :: E_sqr(:)
!         complex(dp),allocatable :: V_least_proj(:,:) ! eigenvectors
!         logical :: isorthogonal
!         integer :: nbases
!         integer :: nsyms
!         integer :: nphis
!         integer :: nirreps
!         integer :: i, i_phi, j, k
!         integer :: verbosity
!         ! verbosity for debugging
!         verbosity = 1
!         ! get bases dimensions
!         nbases = size(rr,1)
!         nsyms  = size(rr,3)
!         nirreps= size(irrep_proj,3)
!         nphis  = sum(irrep_dim)
!         ! allocate space
!         allocate(phi(nbases,nphis))
!         allocate(S(nirreps))
!         allocate(E(nirreps))
!         ! work space
!         allocate(try(nbases))
!         ! allocate zero vector
!         allocate(zeros(nbases))
!         zeros = 0
!         ! get symmetries eigenvectors and eigenvalues
!         allocate(D(nbases,nsyms))
!         allocate(V(nbases,nbases,nsyms))
!         do i = 1, nsyms
!             call am_dgeev(A=real(rr(:,:,i),dp),VR=V(:,:,i),D=D(:,i))
!         enddo
!         ! debug
!         if (debug) then
!         call execute_command_line('mkdir -p '//trim(outfile_dir_sym)//'/debug')
!         call dump(A=class_member       , fname=trim(outfile_dir_sym)//'/debug'//'/outfile.class_member' )
!         call dump(A=rr                 , fname=trim(outfile_dir_sym)//'/debug'//'/outfile.rr'           )
!         call dump(A=irrep_proj         , fname=trim(outfile_dir_sym)//'/debug'//'/outfile.irrep_proj'   )
!         call dump(A=irrep_dim          , fname=trim(outfile_dir_sym)//'/debug'//'/outfile.irrep_dim'    )
!         endif
!         ! get similarity transform which block diagonalizes irrep projections
!         ! call get_irrep_proj_eigvec(irrep_proj=irrep_proj, V_proj=V_proj,S=S_sqr,E=E_sqr)
!         ! call dump(A=V_proj, fname=trim(outfile_dir_sym)//'/outfile.V_proj')
!         ! initialize i_phi
!         i_phi = 0
!         ! loop over irrep projections
!         do i = 1, nirreps
!             if (verbosity.ge.1) write(*,'(a,a)',advance='no') flare, 'irrep '//tostring(i)//': '
!             ! save indices
!             i_phi= i_phi + 1
!             S(i) = i_phi
!             E(i) = S(i) + irrep_dim(i) - 1
!             if     (irrep_dim(i).eq.1) then
!                 if (verbosity.ge.1) write(*,'(a)',advance='no') '1D '
!                 ! project a random vector on the irreducible subspace
!                 ! since the subspace is 1 dimension, the randomness goes away.
!                 phi(:,i_phi) = matmul(irrep_proj(:,:,i), rand(nbases) )
!             elseif (irrep_dim(i).ge.2) then
!                 ! if the irrep is two-dimensional or larger...
!                 if (verbosity.ge.1) write(*,'(a)',advance='no') tostring(irrep_dim(i))//'D '
!                 ! search for an eigenvector of a regular rep symmetry element that has unique (nondegenerate eigenvalues)
!                 if    (successfully_get_least_degenerate_eigenspace(V=V,D=D,irrep_proj_i=irrep_proj(:,:,i),V_least=V_least)) then
!                     ! project the nondegenerate vector onto the irreducible subspace
!                     phi(:,i_phi) = matmul(irrep_proj(:,:,i), V_least(:,1) )
!                 elseif (successfully_restrict_eigenspace(V=V,D=D,irrep_dim_i=irrep_dim(i), &
!                     irrep_proj_i=irrep_proj(:,:,i),V_least=V_least,V_least_proj=V_least_proj)) then
!                     ! use the lowest degenerate eigenvalue to search for a starting eigensubspace
!                     phi(:,i_phi) = matmul(irrep_proj(:,:,i), V_least_proj(:,1) )
!                 else
!                     stop 'need to program this still... '
!                 endif
!                 ! search for partners
!                 search : do k = 1, nclasses
!                     ! exit search when all partners symmetry functions are found 
!                     if (i_phi.eq.E(i)) then
!                         exit search
!                     endif
!                     ! try partner function
!                     try = matmul(rr(:,:,class_member(k,1)), phi(:,S(i)))
!                     ! make sure: non-zero
!                     if (.not.isequal(try,zeros)) then
!                         ! make sure: orthogonal to other vectors in the irrep basis
!                         ! since in the regular rep symmetries have matrix elements (0,1), 
!                         ! it suffices to check whether the vectors are not the same.
!                         ! ideal, a cross product would be taken.
!                         isorthogonal = .true.
!                         search_orthogonal : do j = S(i), i_phi
!                             if (isequal(phi(:,j),try)) then
!                                 isorthogonal = .false.
!                                 exit search_orthogonal
!                             endif
!                         enddo search_orthogonal
!                         if (isorthogonal) then
!                             i_phi = i_phi + 1
!                             phi(:,i_phi) = try/norm(try)
!                         endif 
!                     endif
!                 enddo search
!             else
!                 stop 'ERROR [get_block_similarity_transform]: incorrect irrep_dim'
!             endif
!             ! orthonormalize the irreducible subspace
!             phi(:,S(i):E(i)) = orthonormalize( phi(:,S(i):E(i)) )
!             ! line break
!             if (verbosity.ge.1) write(*,*)
!         enddo
!         if (i_phi.ne.nphis) stop 'ERROR [get_block_similarity_transform]: i_phi /= nirreps'

!     end function   get_block_similarity_transform

    ! <AUX:get_block_similarity_transform>

    function       successfully_restrict_eigenspace(V,D,irrep_dim_i,irrep_proj_i,V_least,V_least_proj) result(is_successful)
        !
        implicit none
        !
        complex(dp), intent(in) :: V(:,:,:) ! eigenvectors
        complex(dp), intent(in) :: D(:,:)   ! eigenvalues
        integer    , intent(in) :: irrep_dim_i
        complex(dp), intent(in) :: irrep_proj_i(:,:)
        complex(dp), intent(in) :: V_least(:,:)
        complex(dp),allocatable, intent(out) :: V_least_proj(:,:) ! eigenvectors
        complex(dp),allocatable :: V_i(:,:)
        complex(dp),allocatable :: ns(:,:)  ! eigenvectors
        integer    ,allocatable :: inds(:)
        logical :: is_successful
        integer :: nsyms, nbases
        integer :: ns_i, ns_j
        integer :: i
        ! set is_successful (is_successful = .true. if dimension of V_least_proj is fully reduced to irrep_dim )
        is_successful = .false.
        ! get dimensions
        nbases = size(V,1)
        nsyms  = size(V,3)
        ! allocate inds
        allocate(inds(nbases))
        inds = [1:nbases]
        ! project onto irreducible subspace and orthonormalize
        V_least_proj = orth_svd( matmul(irrep_proj_i, V_least) )

! something is wrong here... V_least_proj is empty?
            call disp(V_least_proj)
            call disp(irrep_proj_i)
            call disp(V_least)
            call disp(matmul(irrep_proj_i, V_least))

        ! initialize null space dimension counter
        ns_i = nbases
        ! attempt to reduce projected subspace by intsersecting it with eigenspace of regular reps
        reduction_loop : do i = 1, nsyms
            ! get eigenvectors which nonzero eigenvalues
            V_i = orth_svd( V(:, pack(inds,abs(D(:,i)).gt.tiny), i) )
            ! get subspace intersection
            ns  = subspace_intersection(V_i, V_least_proj)
            ! get dimension of intersection
            ns_j = size(ns,2)
            ! check that something was obtained form the intserection
            if (ns_j.ne.0) then
                ! check that the intersection reduced the dimension of the space
                if (ns_j.lt.ns_i) then
                    ! update subspace
                    V_least_proj = orth_svd(ns)
                    ! update subspace dimension
                    ns_i = ns_j
                    ! end if subspace is completely reduced (this is the smallest it will ever get)
                    if (ns_i.eq.irrep_dim_i) then
                        is_successful = .true.
                        exit reduction_loop
                    endif
                endif
            endif
        enddo reduction_loop
        ! 
    end function   successfully_restrict_eigenspace
    
    function       successfully_get_least_degenerate_eigenspace(V,D,irrep_proj_i, V_least) result(is_successful)
        !
        implicit none
        !
        complex(dp), intent(in) :: V(:,:,:) ! eigenvectors
        complex(dp), intent(in) :: D(:,:)   ! eigenvalues
        complex(dp), intent(in) :: irrep_proj_i(:,:) ! projection
        complex(dp),allocatable, intent(out) :: V_least(:,:)
        logical :: is_successful
        complex(dp),allocatable :: V_try(:,:)
        complex(dp),allocatable :: V_j(:,:) ! eigenvectors
        complex(dp),allocatable :: D_j(:)   ! eigenvalues
        complex(dp),allocatable :: D_j_unique(:) ! eigenvalues
        logical    ,allocatable :: mask(:)
        integer    ,allocatable :: inds(:)
        integer :: n_D_j_unique_s
        integer :: nsyms, nbases
        logical :: pass
        integer :: degen_i ! counter
        integer :: degen_j ! counter
        integer :: j, k
        ! successful if 1D space is returned
        is_successful = .false.
        ! get dimnsions
        nbases = size(V,1)
        nsyms  = size(V,3)
        ! allocate stuff
        allocate(V_j(nbases,nbases))
        allocate(D_j(nbases))
        ! allocate masks
        allocate(mask(nbases))
        ! allocate inds
        allocate(inds(nbases))
        inds = [1:nbases]
        ! initialize degeneracy counter
        degen_i = nbases
        ! loop over symmetries
        do j = 1, nsyms
            ! get eigenvector/values
            V_j = V(:,:,j)
            D_j = D(:,j)
            ! filter eigenvectors belonging to least degenerate eigenvalue
            D_j_unique = unique(D_j)
            ! get number of unique eigenvalues
            n_D_j_unique_s = size(D_j_unique)
            ! loop over unique eigenvalues
            do k = 1, n_D_j_unique_s
                ! mask eigenvalues which match the k-th unique value
                where (abs(D_j-D_j_unique(k)).lt.tiny)
                    mask = .true.
                elsewhere
                    mask = .false.
                endwhere
                ! get numer of eigenvalues matches
                degen_j = count(mask)
                ! check if this eigenvalue has fewer matchs than previous ones
                if (degen_j.lt.degen_i) then
                    ! get orthonormalized columns
                    V_try = orth_svd( V(:,pack(inds,mask),j) )
                    ! check that this eigenvectors which are not in the null space of the projection oprators
                    pass = .true.
                    ! seems impossible to find a set of vectors that are not simultaenosuly in the null space of at least one irrep.
                    ! do i = i, nirreps
                    if ( any(all( abs(matmul(irrep_proj_i,V_try)).lt.tiny ,1) ) ) then
                        pass = .false.
                        exit
                    endif
                    ! enddo
                    !
                    if (pass) then
                        ! updat degeneracy counter to reflect smaller eigenspace
                        degen_i = degen_j
                        ! update eigenspace
                        if (allocated(V_least)) deallocate(V_least)
                        allocate(V_least, source = V_try)
                        ! if a 1D space is found return (it doesn't get much smaller than this!)
                        if (degen_j.eq.1) then 
                            is_successful = .true.
                            return
                        endif
                    endif
                endif
            enddo
        enddo
        !
        if (.not.pass) then
            stop 'ERROR [successfully_get_least_degenerate_eigenspace]: failed to find basis'
        endif
    end function   successfully_get_least_degenerate_eigenspace


    ! </AUX:get_block_similarity_transform>

    ! identifier functions which operate on identifiers

    function       id_member(id) result(class_member)
        !
        implicit none
        !
        integer, intent(in) :: id(:)
        integer, allocatable :: class_nelements(:) ! number of elements in each class
        integer, allocatable :: class_member(:,:) ! members(nclass,maxval(class_nelements))
        integer :: i, j, k
        integer :: nsyms
        integer :: nclasses
        !
        nsyms = size(id,1)
        !
        nclasses = maxval(id)
        !
        class_nelements = id_nelements(id)
        !
        allocate(class_member(nclasses,maxval(class_nelements)))
        class_member = 0
        do i = 1, nclasses
            k=0
            do j = 1, nsyms
                if (id(j).eq.i) then
                    k=k+1
                    class_member(i,k) = j
                endif
            enddo
        enddo
    end function   id_member

    function       id_representative(id) result(class_representative)
        !
        implicit none
        !
        integer, intent(in) :: id(:)
        integer, allocatable :: class_member(:,:) ! members(nclass,maxval(class_nelements))
        integer, allocatable :: class_representative(:)
        !
        class_member = id_member(id)
        !
        allocate(class_representative, source=class_member(:,1))
        !
    end function   id_representative

    function       id_nelements(id) result(class_nelements)
        !
        implicit none
        !
        integer, intent(in)  :: id(:)
        integer, allocatable :: class_nelements(:)
        integer :: nclasses
        integer :: i
        !
        nclasses = maxval(id)
        !
        allocate(class_nelements(nclasses))
        !
        do i = 1, nclasses
            class_nelements(i) = count(i.eq.id)
        enddo
    end function   id_nelements

    ! get point symmetry name

    function       decode_pointsymmetry(ps_id) result(schoenflies)
        !
        implicit none
        !
        integer, intent(in) :: ps_id
        character(len=string_length_schoenflies) :: schoenflies
        !
        select case (ps_id)
        case(1);  schoenflies = 'e'
        case(2);  schoenflies = 'c_2'
        case(3);  schoenflies = 'c_3'
        case(4);  schoenflies = 'c_4'
        case(5);  schoenflies = 'c_6'
        case(6);  schoenflies = 'i'
        case(7);  schoenflies = 's_2'
        case(8);  schoenflies = 's_6'
        case(9);  schoenflies = 's_4'
        case(10); schoenflies = 's_3'
        case default
            stop 'Schoenflies code unknown.'
        end select
    end function   decode_pointsymmetry


end module am_symmetry_tables
