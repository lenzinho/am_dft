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
        real(dp)   , allocatable :: matrices(:,:,:)         ! 
    end type am_class_conjugacy_class

    type, public :: am_class_character_table
        integer                  :: nirreps
        complex(dp), allocatable :: chartab(:,:)            ! character table ( irreps * classes )
        character(:),allocatable :: irrep_label(:)          ! irrep label
        integer    , allocatable :: irrep_dim(:)            ! irrep dimension
        integer    , allocatable :: irrep_id(:)             ! identifies the irrep to which the column of  
        integer    , allocatable :: irrep_decomp(:)         ! number of times irrep appears in the deocmposition of rep
        complex(dp), allocatable :: irrep_proj(:,:,:)       ! 
        complex(dp), allocatable :: irrep_proj_V(:,:)       ! eigenvectors of projection
        complex(dp), allocatable :: block_proj(:,:)         ! eigenvectors of projection
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

    type, private :: am_class_chartab_printer
        integer :: k
        ! complex to symbol
        complex(dp), allocatable :: s(:)
        character(:), allocatable :: str
        integer :: char_start
        ! print formats
        character(50) :: fmts(5)
    contains
        procedure :: printer_initialize
        procedure :: printer_print_chartab_header
        procedure :: printer_chi2symb
        procedure :: printer_print_bar
        procedure :: printer_print_definitions
    end type am_class_chartab_printer

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
        integer :: nsyms, nbases
        integer :: i, j, h ! loop variabls
        integer :: k, kmax ! for progress bar
        ! get dimensions
        nbases  = size(sym,1)
        nsyms   = size(sym,3)
        ! check that first element is the identity
        if (.not.isequal(sym(:,:,1),eye(nbases))) stop 'ERROR [get_multab]: identity not first'
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
            if (index(flags,'seitz').ne.0) then
                W(1:3,4) = modulo(W(1:3,4)+tiny,1.0_dp)-tiny
            endif
            ! get matching symmetry
            do h = 1, nsyms
                if (isequal(sym(:,:,h),W)) then
                    multab(i,j) = h
                    exit
                endif
            enddo
        enddo
        enddo
        ! quick consitency check for closure. not comprehensive
        do i = 1, nsyms
            if (2*sum(multab(:,i)).ne.nsyms*(nsyms+1)) then
                stop 'ERROR [get_multab]: multab is incorrect'
            endif
        enddo
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

    function       get_regular_rep(multab) result(sym)
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
        integer , intent(in)  :: multab(:,:)
        real(dp), allocatable :: sym(:,:,:)
        integer , allocatable :: sortmat(:,:)
        integer , allocatable :: inds(:)
        integer , allocatable :: multab_sorted(:,:)
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
        allocate(sym(nsyms,nsyms,nsyms))
        do i = 1, nsyms
            where (multab_sorted.eq.i) 
                sym(:,:,i) = 1
            elsewhere
                sym(:,:,i) = 0
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

    ! subgroup

!     function       get_subgroup_member(multab)
!         !
!         implicit none
!         !
!         integer, intent(in) :: multab(:,:)
!         integer :: nsyms
!         !
!         nsyms = size(multab,1)
!         ! get all cyclic subgroups
        


!     end function

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

    ! character table

    function       get_chartab(multab,class_nelements,class_member,ps_id) result(chartab)
        !
        use am_rank_and_sort
        !
        implicit none
        !
        integer, intent(in) :: multab(:,:)
        integer, intent(in) :: ps_id(:)
        integer, intent(in) :: class_nelements(:) ! number of elements in each class
        integer, intent(in) :: class_member(:,:) ! members(nclass,maxval(class_nelements))
        complex(dp), allocatable :: chartab(:,:) ! class constant chartab which becomes character chartab (can be complex!)
        real(dp), allocatable :: Re(:,:)
        real(dp), allocatable :: Im(:,:)
        integer , allocatable :: H(:,:,:) ! class multiplication 
        integer , allocatable :: irrep_dim(:) ! sym dimensions
        integer , allocatable :: indices(:)
        integer  :: nsyms
        integer  :: nirreps  ! just to make things clearer; always equal to nclasses
        integer  :: nclasses
        real(dp) :: wrk
        integer  :: i, j, k
        ! get dimensions
        nsyms    = size(multab,2)
        nclasses = size(class_member,1)
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
                stop 'ERROR [get_chartab]: class order is not a factor of group order'
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

    function       get_muliken_label(chartab,class_nelements,class_member,ps_id) result(irrep_label)
        !
        implicit none
        !
        complex(dp), intent(in) :: chartab(:,:) ! class constant chartab which becomes character chartab (can be complex!)
        integer, intent(in) :: class_nelements(:)
        integer, intent(in) :: class_member(:,:)
        integer, intent(in) :: ps_id(:)
        character(:), allocatable :: irrep_label(:)
        integer :: i, j, k
        integer :: nirreps
        integer :: nclasses
        integer :: class_containing_identity
        integer :: class_containing_inversion
        integer :: class_containing_highsym
        integer, allocatable :: highsymlist(:)
        ! they are always identical...
        nirreps  = size(chartab,1)
        nclasses = size(chartab,2)
        ! create irrep label
        allocate(character(7) :: irrep_label(nirreps))
        ! get clases containing stuff
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

    subroutine     print_chartab(chartab,class_nelements,class_member,irrep_label,ps_id,chi_rep,rep_label)
        !
        implicit none
        !
        complex(dp) , intent(in) :: chartab(:,:) ! class constant chartab which becomes character chartab (can be complex!)
        integer     , intent(in) :: class_nelements(:)
        integer     , intent(in) :: class_member(:,:)
        character(*), intent(in) :: irrep_label(:)
        integer     , intent(in) :: ps_id(:)
        complex(dp) , intent(in), optional :: chi_rep(:,:)
        character(*), intent(in), optional :: rep_label(:)
        type(am_class_chartab_printer) :: printer
        integer :: i, j
        integer :: nirreps
        integer :: nclasses
         ! rep to irrep decompostion (used only if chi_rep nd rep_label are present)
        integer :: nreps
        integer, allocatable :: beta(:,:)
        !
        ! get dimensions
        nirreps  = size(chartab,1)
        nclasses = size(chartab,2)
        !
        call printer%printer_initialize()
        !
        call printer%printer_print_chartab_header(nclasses=nclasses,&
            class_nelements=class_nelements,class_member=class_member,ps_id=ps_id)
        !
        ! print irrep characters (character table)
        do i = 1, nirreps
            write(*,'(5x,i2,a8)',advance='no') i, irrep_label(i)
            do j = 1, nclasses
                write(*,printer%fmts(3),advance='no') printer%printer_chi2symb(chi=chartab(i,j))
            enddo
            write(*,*)
        enddo
        ! print rep characters
        if (present(chi_rep).and.present(rep_label)) then
            ! check for error
            if (size(chi_rep,2).ne.nclasses) stop 'ERROR [print_chartab]: dimensions of rep character does not match number of classes'
            ! rep characters
            call printer%printer_print_bar(nclasses=nclasses)
            nreps = size(chi_rep,1)
            ! allocate beta (reduction coefficients)
            allocate(beta(nreps,nirreps))
            ! print rep characters
            do i = 1, nreps
                write(*,'(5x,i2,a8)',advance='no') i, rep_label(i)
                do j = 1, nclasses
                    write(*,printer%fmts(3),advance='no') printer%printer_chi2symb(chi=chi_rep(i,j))
                enddo
                write(*,*)
            enddo
        endif
        ! print rep decompositions
        if (present(chi_rep).and.present(rep_label)) then
            call printer%printer_print_bar(nclasses=nclasses)
            write(*,'(5x,a)') 'Decompositions:'
            do i = 1, nreps
                beta(i,:) = get_reduction_coefficient(chartab=chartab,class_nelements=class_nelements,chi_rep=chi_rep(i,:))
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
        ! print definitions
        call printer%printer_print_definitions()
        !
    end subroutine print_chartab

!     function      get_aux_reduction_matrix(chartab,)
!         !
!         ! matlab: ( c_4v.ct.chartab .* spread(c_4v.cc.nelements,1,c_4v.ct.nirreps))'/c_4v.nsyms
!         implicit none
!         !

!         !
!     end function  get_aux_reduction_matrix


    function       get_reduction_coefficient(chartab,class_nelements,chi_rep) result(beta)
        !
        implicit none
        !
        complex(dp), intent(in) :: chartab(:,:)
        complex(dp), intent(in) :: chi_rep(:)
        integer, intent(in) :: class_nelements(:)
        integer,allocatable :: beta(:)
        complex(dp) :: try
        integer :: nsyms
        integer :: nirreps
        integer :: j
        !
        nirreps = size(chartab,1)
        nsyms   = sum(class_nelements)
        ! initialize
        allocate(beta(nirreps))
        ! calculate
        do j = 1, nirreps
            try = dot(chi_rep,chartab(j,:)*class_nelements)/real(nsyms,dp)
            if (abs(aimag(try)).gt.tiny)                stop 'irrep reduction coefficient is not real'
            if (abs(nint(real(try))-real(try)).gt.tiny) stop 'irrep reduction coefficient is not an integer'
            ! get save integer beta
            beta(j) = nint(real(try))
        enddo
    end function   get_reduction_coefficient

    function       get_chi_subduced(supergroup_chartab,supergroup_class_id,supergroup_id) result(chi_subduced)
        !
        ! matlab:
        ! o_h.ct.chartab(:,unique(o_h.cc.id(c_4v.supergroup_id)))
        implicit none
        !
        complex(dp), intent(in) :: supergroup_chartab(:,:)
        integer    , intent(in) :: supergroup_class_id(:)
        integer    , intent(in) :: supergroup_id(:)
        integer    ,allocatable :: chi_subduced(:,:)
        !
        allocate(chi_subduced,source=supergroup_chartab(:,unique(supergroup_class_id(supergroup_id))))
        !
    end function   get_chi_subduced

    function       get_irrep_dimension(chartab,class_nelements,class_member,ps_id) result(irrep_dim)
        !
        implicit none
        !
        complex(dp), intent(in) :: chartab(:,:)
        integer    , intent(in) :: class_nelements(:)
        integer    , intent(in) :: class_member(:,:)
        integer    , intent(in) :: ps_id(:)
        integer    , allocatable :: irrep_dim(:)
        integer :: i ! class
        integer :: j ! irrep
        integer :: nirreps
        integer :: nclasses
        ! get dimensions
        nirreps  = size(chartab,1)
        nclasses = size(chartab,1)
        !
        ! find class containing identity symmetry
        id_search : do i = 1, nclasses
        do j = 1, class_nelements(i)
            if (ps_id(class_member(i,j)).eq.1) then
                exit id_search
            endif
        enddo
        enddo id_search
        !
        allocate(irrep_dim(nirreps))
        do j = 1, nirreps
            irrep_dim(j) = nint(real(chartab(j,i)))
        enddo
    end function   get_irrep_dimension

    ! irrep construction representation

    function       get_rep_characters(sym,class_member) result(chi)
        !
        implicit none
        !
        real(dp), intent(in) :: sym(:,:,:)
        integer , intent(in) :: class_member(:,:)
        complex(dp), allocatable :: chi(:)
        integer :: nclasses
        integer :: i
        !
        nclasses = size(class_member,1)
        !
        allocate(chi(nclasses))
        do i = 1, nclasses
            chi(i) = trace( sym(:,:,class_member(i,1)) )
        enddo
        !
    end function   get_rep_characters

    function       get_irrep_decomp(chi,chartab,class_nelements) result(beta)
        !
        implicit none
        !
        complex(dp), intent(in) :: chi(:) ! of rep's class reps.
        complex(dp), intent(in) :: chartab(:,:)
        integer    , intent(in) :: class_nelements(:)
        integer    ,allocatable :: beta(:)
        integer :: nsyms
        integer :: nirreps
        complex(dp) :: wrk
        integer :: j
        ! get dimensions
        nsyms   = sum(class_nelements)
        nirreps = size(chartab,1)
        ! allocate space
        allocate(beta(nirreps))
        ! get decomposition coefficients
        do j = 1, nirreps
            wrk = dot(chi,chartab(j,:)*class_nelements)/real(nsyms,dp)
            if (abs(aimag(wrk)).gt.tiny)                stop 'ERROR [get_irrep_decomp]: irrep_decomp coefficient is not real'
            if (abs(nint(real(wrk))-real(wrk)).gt.tiny) stop 'ERROR [get_irrep_decomp]: irrep_decomp coefficient is not an integer'
            beta(j) = nint(real(wrk))
        enddo
    end function   get_irrep_decomp

    function       get_irrep_id(irrep_dim,irrep_decomp) result(irrep_id)
        !
        implicit none
        !
        integer, intent(in) :: irrep_decomp(:)
        integer, intent(in) :: irrep_dim(:)
        integer,allocatable :: irrep_id(:)
        integer :: nirreps
        integer :: i,j,k,n
        ! get dimesnions
        nirreps = size(irrep_dim)
        ! get indices
        n = sum( irrep_dim * irrep_decomp )
        ! allocate space
        allocate(irrep_id(n))
        !  save n's
        n = 0
        do i = 1, nirreps
        do j = 1, irrep_dim(i)
        do k = 1, irrep_decomp(i)
            n = n + 1
            irrep_id(n) = i
        enddo
        enddo
        enddo
        !
    end function   get_irrep_id

    function       get_class_matrices(sym,class_nelements,class_member) result(class_matrices)
        !
        implicit none
        !
        real(dp), intent(in) :: sym(:,:,:)
        integer , intent(in) :: class_nelements(:)
        integer , intent(in) :: class_member(:,:)
        real(dp),allocatable :: class_matrices(:,:,:)
        integer :: i,j
        integer :: nbases
        integer :: nclasses
        ! get number of things
        nbases  = size(sym,1)
        nclasses = size(class_member,1)
        ! class matrices
        allocate(class_matrices(nbases,nbases,nclasses))
        class_matrices = 0
        do i = 1, nclasses
        do j = 1, class_nelements(i)
            class_matrices(:,:,i) = class_matrices(:,:,i) + sym(:,:,class_member(i,j))
        enddo
        enddo
        !
    end function   get_class_matrices

    function       get_irrep_proj(sym,chartab,class_id,irrep_dim) result(irrep_proj)
        !
        implicit none
        !
        real(dp)   , intent(in) :: sym(:,:,:)
        complex(dp), intent(in) :: chartab(:,:)
        integer    , intent(in) :: class_id(:)
        integer    , intent(in) :: irrep_dim(:)
        complex(dp),allocatable :: irrep_proj(:,:,:)
        integer :: nbases
        integer :: nsyms
        integer :: nirreps
        integer :: i,j
        ! get dimensions
        nbases  = size(sym,1)
        nsyms   = size(sym,3)
        nirreps = size(irrep_dim)
        ! allocate space for irrep projection operator
        allocate(irrep_proj(nbases,nbases,nirreps))
        irrep_proj = 0
        ! loop over irreos
        do j = 1, nirreps
            ! sum over symmetries
            do i = 1, nsyms
                irrep_proj(:,:,j) = irrep_proj(:,:,j) + chartab(j, class_id(i) ) * sym(:,:,i)
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

    function       get_irrep_proj_V(irrep_proj,irrep_id) result(irrep_proj_V)
        !
        implicit none
        !
        complex(dp), intent(in) :: irrep_proj(:,:,:) ! can be complex because of complex characters
        integer    , intent(in) :: irrep_id(:)
        complex(dp),allocatable :: irrep_proj_V(:,:)
        real(dp)   ,allocatable :: D(:)
        complex(dp),allocatable :: V(:,:)
        integer :: nbases
        integer :: nirreps
        integer :: i
        ! get dimensions
        nbases  = size(irrep_proj,1)
        nirreps = size(irrep_proj,3)
        ! allocate space
        allocate(D(nbases))
        allocate(V(nbases,nbases))
        ! get irrep_proj_V
        allocate(irrep_proj_V(nbases,nbases))
        irrep_proj_V = 0
        !
        do i = 1, nirreps
        if (any(irrep_id.eq.i)) then
            ! get eigenvectors/eigenvalues (HERMITIAN)
            call am_zheev(A=irrep_proj(:,:,i), V=V, D=D)
            ! save projection
            irrep_proj_V(:, selector(irrep_id.eq.i) ) = orth_svd( V(:, selector(abs(D).gt.tiny)) )
        endif
        enddo
        !
    end function   get_irrep_proj_V

    function       get_block_proj(sym,irrep_proj_V,irrep_id) result(block_proj)  
        !
        implicit none
        !
        real(dp)   , intent(in) :: sym(:,:,:)
        complex(dp), intent(in) :: irrep_proj_V(:,:)
        integer    , intent(in) :: irrep_id(:)
        complex(dp),allocatable :: block_proj(:,:)
        complex(dp),allocatable :: V(:,:)
        complex(dp),allocatable :: H(:,:)
        complex(dp),allocatable :: M(:,:,:)
        real(dp)   ,allocatable :: C(:,:)
        real(dp)   ,allocatable :: D(:)
        real(dp)   ,allocatable :: D_unique(:)
        integer    ,allocatable :: inds(:)
        integer    ,allocatable :: inds2(:)
        integer    ,allocatable :: blocks(:)
        logical :: isdecomposed
        integer :: nbases
        integer :: nsyms
        integer :: nirreps
        integer :: ninds
        integer :: i,j,k,n
        ! get dimensions
        nbases  = size(sym,1)
        nsyms   = size(sym,3)
        nirreps = maxval(irrep_id)
        ! allocate H
        allocate(H(nbases,nbases))
        ! allocate block diagonalized regular rep
        allocate(M(nbases,nbases,nsyms))
        ! allocate block structure
        allocate(C(nbases,nbases))
        ! initialize decomposition loop ...
        ! ... by using eigenvectors of irrep projection operator to perform preliminar block decomposition 
        allocate(block_proj,source=irrep_proj_V)
        if (size(irrep_id).ne.nbases) stop 'ERROR [get_block_proj]: irrep_id /= nbases'
        allocate(blocks,source=irrep_id)
        ninds = nirreps
        do i = 1, nsyms
            M(:,:,i) = matmul(matmul( adjoint(irrep_proj_V), sym(:,:,i)), irrep_proj_V )
        enddo
        ! ALTERNATIVE : initialize decomposition loop ... (no initial projection) DOESN'T WORK. H IS NOT HERMITIAN FOR SOME REASON...
        !  block_proj = eye(nbases)
        !  allocate(blocks(nbases))
        !  blocks = 1
        !  ninds = 1
        !   do i = 1, nsyms
        !      M(:,:,i) = sym(:,:,i)
        !   enddo
        ! should be done by 10 loops
        decompose_loop : do k = 1, 10
            isdecomposed = .true.
            ! loop over block structures
            do j = 1, ninds
                ! call disp(X=blocks,orient='row',title=tostring(j)//' ')
                ! get blocks
                inds = selector(blocks.eq.j)
                ! get size of block diagonalized basis
                n = count(blocks.eq.j)
                ! get dixon H
                H = 0
                H(1:n,1:n) = get_dixon_H( M(inds,inds,:) )
                ! exit loop if identity is returned
                if (.not.isequal(H(1:n,1:n), H(1,1)*eye(n) )) then
                    ! check if atleast one block was decomposed
                    isdecomposed = .false.
                    ! check that H is hermitian
                    if (.not.isequal(adjoint(H(1:n,1:n)),H(1:n,1:n))) stop 'ERROR [get_block_proj]: H is not Hermitian'
                    ! get eigenvectors using complex Hermitian routine 
                    call am_zheev(A=H(1:n,1:n),V=V,D=D)
                    ! get unique eigenvalues
                    D_unique = unique(D)
                    ! orthonormalize the eigenbasis; each unique eigenvalue at a time
                    do i = 1, size(D_unique)
                        ! get basis corresponding to unique eigenvalue
                        inds2 = selector( abs(D-D_unique(i)).lt.tiny )
                        ! orthonormalize
                        V(:,inds2) = orth_svd( V(:,inds2) )
                    enddo
                    ! update M
                    do i = 1, nsyms
                        M(inds,inds,i) = matmul(matmul(adjoint(V), M(inds,inds,i)), V)
                    enddo
                    ! update block_proj
                    block_proj(:,inds) = matmul(block_proj(:,inds), V)
                endif
            enddo
            ! get block structure
            where ( sum(abs(M),3).gt.tiny )
                C = 1
            elsewhere
                C = 0
            endwhere
            ! reduce block structure to row echelon form
            C = rref(C)
            ! get number of indices
            ninds = count( any(abs(C).gt.tiny,2) )
            do i = 1, ninds
                blocks( selector(abs(C(i,:)).gt.tiny) ) = i
            enddo
            if (count(any(abs(C).gt.tiny,1)).ne.nbases) stop 'ERROR [get_block_proj]: not all rows have been assigned to a block'
            ! check if finished
            if (isdecomposed) exit decompose_loop
        enddo decompose_loop
        if (.not.isdecomposed) stop 'ERROR [get_block_proj]: irrep decomposition taking longer than 10 iterations...'
        ! sort projection based on block structure
        k = 0
        get_blocks : do i = 1, nbases
        do j = 1, nbases
        if (abs(C(i,j)).gt.tiny) then
            k = k + 1
            blocks(k) = j
            ! Note: C is not in RREF, hence the exit here.
            if (k.eq.nbases) exit get_blocks
        endif
        enddo
        enddo get_blocks
        ! sort block_proj
        block_proj = block_proj(:,blocks)
        contains
        function       get_dixon_H(rr) result(H)
            ! REF: 
            ! MATHEMATICS OF COMPUTATION, VOLUME 24, NUMBER 111, JULY, 1970
            ! Computing Irreducible Representations of Groups
            ! By John D. Dixon
            implicit none
            !
            complex(dp), intent(in) :: rr(:,:,:)
            complex(dp),allocatable :: Hrs(:,:)
            complex(dp),allocatable :: H(:,:)
            logical :: isreducible
            integer :: nbases, nsyms
            integer :: r,s,i
            !
            nbases=size(rr,1)
            nsyms =size(rr,3)
            !
            isreducible = .false.
            !
            allocate(Hrs(nbases,nbases))
            allocate(H(nbases,nbases))
            !
            search : do r = 1, nbases
            do s = 1, nbases
                ! construct Hrs
                Hrs(1:nbases,1:nbases) = 0
                if     (r.eq.s) then
                    Hrs(r,r) = 1.0_dp
                elseif (r.gt.s) then
                    Hrs(r,s) = 1.0_dp
                    Hrs(s,r) = 1.0_dp
                elseif (r.lt.s) then
                    Hrs(r,s) = cmplx_i
                    Hrs(s,r) =-cmplx_i
                endif
                ! construct H
                H(1:nbases,1:nbases) = 0
                do i = 1, nsyms
                H = H + matmul(matmul(adjoint(rr(:,:,i)), Hrs), rr(:,:,i))
                enddo
                H = H/real(nsyms,dp)
                ! if H is not a scalar, then the rep is reducible
                if (.not.isequal(H(1,1)*eye(nbases),H)) then
                    isreducible = .true.
                    exit search
                endif
            enddo
            enddo search
            ! if it isn't reducible return identity
            if (.not.isreducible) H = eye(nbases)
            !
            ! At the end of the previous post we saw that in order to decompose a representation (V,ρ), it is
            ! enough to find a non-scalar matrix T that commutes with ρ(g) for every g ∈ G . This first step finds
            ! a Hermitian non-scalar H that commutes with (G) (if there is one to be found). Let Ers denote the
            ! n×nn×n matrix with a 11 in the (r,s)th entry and zeros everywhere else. Here nn is the dimension of
            ! VV in the representation (V,ρ). Define:
            !
            ! H_{rs} = 
            ! \begin{cases}
            !            E_{rr}  &\text{if } r = s \\
            !   E_{rs} + E_{sr}  &\text{if } r > s \\
            ! i(E_{rs} - E_{sr}) &\text{if } r < s,
            ! \end{cases}
            ! 
            ! then the set of matrices HrsHrs forms a Hermitian basis for the n×n matrices over ℂ. Now for each r,s 
            ! compute the sum: H=1|G|∑g∈Gρ(g)∗Hrsρ(g). Observe that HH has the following properties:
            ! 1) it is hermitian
            ! 2) it commutes with ρ(g)ρ(g) for all g∈Gg∈G
            !
            ! If ρρ is irreducible, then HH is a scalar matrix for all r,sr,s. Otherwise, it turns out that there will 
            ! be some r,s such that HH is non-scalar (this is due to the fact that the Hrs matrices form a basis of 
            ! the n×n matrices.
            !
        end function   get_dixon_H
    end function   get_block_proj

    subroutine     print_blocks(sym,block_proj)
        !
        implicit none
        !
        real(dp)   , intent(in) :: sym(:,:,:)
        complex(dp), intent(in) :: block_proj(:,:)
        real(dp)   ,allocatable :: M(:,:,:)
        integer :: i
        !
        allocate(M,mold=sym)
        !
        do i = 1, size(sym,3)
            M(:,:,i) = matmul(matmul(adjoint(block_proj), sym(:,:,i)), block_proj)
        enddo
        !
        call spy( sum(abs(M),3).gt.tiny )
        !
    end subroutine print_blocks





!     function       get_irrep_diag(sym,chartab,irrep_proj,irrep_dim,class_member) result(irrep_diag)
!         !
!         implicit none
!         !
!         real(dp)   , intent(in) :: sym(:,:,:)
!         complex(dp), intent(in) :: chartab(:,:)
!         complex(dp), intent(in) :: irrep_proj(:,:,:) ! projection
!         integer    , intent(in) :: irrep_dim(:)
!         integer    , intent(in) :: class_member(:,:)
!         complex(dp),allocatable :: irrep_diag(:,:)
!         complex(dp),allocatable :: phi(:,:)
!         complex(dp),allocatable :: wrk(:,:)
!         complex(dp),allocatable :: tr(:,:)
!         complex(dp),allocatable :: Q(:)
!         integer    ,allocatable :: S(:)
!         integer    ,allocatable :: E(:)
!         integer :: nirreps, nclasses, nbases, nsyms
!         integer :: i, j, k, n
!         ! get dimensions
!         nsyms   = size(rr,3)
!         nbases  = size(rr,1)
!         nirreps = size(irrep_proj,3)
!         nclasses= nirreps
!         ! allocate phi
!         allocate(phi(nbases,sum(irrep_dim)))
!         phi = 0
!         ! allocate local character table
!         allocate(tr(nirreps,nsyms))
!         ! allocate irrep_diag
!         allocate(irrep_diag(sum(irrep_dim),nsyms))
!         irrep_diag = 0
!         ! allocate start/end
!         allocate(S(nirreps))
!         allocate(E(nirreps))
!         ! allocate workspace
!         allocate(wrk(nbases,nsyms))
!         ! allocate random vector
!         allocate(Q(nbases))
!         ! Q = sqrt(real([1:nbases],dp))
!         Q = pack( rand(nbases,1), .true.)
!         ! 
!         n = 0
!         do j = 1, nirreps
!             ! start/end
!             n = n + 1
!             S(j) = n
!             n = n + irrep_dim(j) - 1 
!             E(j) = n
!             ! project a random vector onto the j-th irreducible subspace
!             wrk(:,1) = matmul(irrep_proj(:,:,j), Q )
!             ! generate partner symmetry functions
!             if (irrep_dim(j).ge.2) then
!                 k = 2
!                 do i = 2, nsyms
!                     wrk(:,k) = matmul(sym(:,:,i), wrk(:,1))
!                     if (k.eq.rank_svd(wrk(:,1:k))) then
!                     if (k.eq.irrep_dim(j)) then
!                         exit
!                     else
!                         k = k + 1
!                     endif
!                     endif
!                 enddo
!             endif
!             ! save phi
!             phi(:,S(j):E(j)) = orth_svd( wrk(:,1:irrep_dim(j)) )
!         enddo
!         ! get diagonal symmetry elements
!         do i = 1, nsyms
!         do j = 1, nirreps
!             do k = 1, irrep_dim(j)
!                 irrep_diag(S(j)+k-1,i) = dot_product(matmul( conjg(phi(:,S(j)+k-1)), sym(:,:,i)), phi(:,S(j)+k-1))
!             enddo
!             tr(j,i) = sum( irrep_diag(S(j):E(j),i) )
!         enddo
!         enddo
!         ! check that trace matches value in character table
!         do j = 1, nirreps
!         do i = 1, nclasses
!             if ( abs(chartab(j,i)-tr(j,class_member(i,1))) .gt. tiny ) then
!                 stop 'ERROR [get_irrep_diag]: character table mismatch'
!             endif
!         enddo
!         enddo
!     end function   get_irrep_diag

!     function       get_wigner_proj(rr,irrep_diag,irrep_dim,class_matrices) result(wigner_proj)
!         !
!         implicit none
!         !
!         integer    , intent(in) :: rr(:,:,:)
!         complex(dp), intent(in) :: irrep_diag(:,:)
!         integer    , intent(in) :: irrep_dim(:)
!         integer    , intent(in) :: class_matrices(:,:,:)
!         complex(dp),allocatable :: wigner_proj(:,:,:,:) ! wigner_proj(nbases,nbases, maxval(irrep_dim), nirreps )
!         integer :: nirreps, nsyms, nbases, nclasses
!         integer :: i,j,k,n
!         !
!         ! get number of things
!         nbases  = size(rr,1)
!         nsyms   = nbases
!         nclasses= size(class_matrices,3)
!         nirreps = nclasses
!         ! allocate space for wigner projection operator
!         allocate(wigner_proj(nbases,nbases,maxval(irrep_dim),nirreps))
!         wigner_proj = 0
!         ! loop over irreps
!         n = 0
!         do j = 1, nirreps
!             ! loop over dimensions of irreps
!             do k = 1, irrep_dim(j)
!                 ! n is a compound irrep and irrep dim index
!                 n = n + 1
!                 ! wooten eq 6.24: sum over symmetries
!                 do i = 1, nsyms
!                 wigner_proj(:,:,k,j) = wigner_proj(:,:,k,j) + rr(:,:,i) * irrep_diag(n,i)
!                 enddo
!                 wigner_proj(:,:,k,j) = wigner_proj(:,:,k,j) * irrep_dim(j)/real(nsyms,dp)
!             enddo
!         enddo
!         ! check that each wigner projection operator is hermitian
!         do j = 1, nirreps
!         do k = 1, irrep_dim(j)
!             if (.not.isequal(adjoint(wigner_proj(:,:,k,j)),wigner_proj(:,:,k,j))) then
!                 stop 'ERROR [get_irrep_projection]: irrep_proj is not hermitian'
!             endif
!         enddo
!         enddo
!         ! check that each wigner projection operator is idempotent (A^2 = A)
!         do j = 1, nirreps
!         do k = 1, irrep_dim(j)
!         if (.not.isequal(matmul(wigner_proj(:,:,k,j),wigner_proj(:,:,k,j)),wigner_proj(:,:,k,j))) then
!             call disp(X=wigner_proj(:,:,k,j),title=tostring(k)//tostring(j),style='underline')
!             stop 'ERROR [get_irrep_projection]: irrep_proj is not idempotent'
!         endif
!         enddo
!         enddo
!         !
!     end function   get_wigner_proj

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

    ! printer for chartab

    subroutine     printer_print_definitions(printer)
        !
        implicit none
        !
        class(am_class_chartab_printer), intent(in) :: printer
        integer :: i
        !
        if (printer%k.ne.0) then
            write(*,'(5x,a)') 'Definitions:'
            do i = 1, printer%k
                write(*,'(5x,a)') char(printer%char_start+printer%k)//' = '//tostring(printer%s(i))
            enddo
        endif
    end subroutine printer_print_definitions

    subroutine     printer_print_chartab_header(printer,nclasses,class_nelements,class_member,ps_id)
        !
        implicit none
        !
        class(am_class_chartab_printer), intent(in) :: printer
        integer     , intent(in) :: nclasses
        integer     , intent(in) :: class_nelements(:)
        integer     , intent(in) :: class_member(:,:)
        integer     , intent(in) :: ps_id(:)
        integer :: i
        !
        ! start printing
        write(*,printer%fmts(2),advance='no') 'class'
        do i = 1, nclasses
            write(*,printer%fmts(1),advance='no') i
        enddo
        write(*,*)
        !
        write(*,printer%fmts(2),advance='no') 'elements'
        do i = 1, nclasses
            write(*,printer%fmts(1),advance='no') class_nelements(i)
        enddo
        write(*,*)
        !
        write(*,printer%fmts(2),advance='no') 'repr.'
        do i = 1, nclasses
            write(*,printer%fmts(3),advance='no') trim(decode_pointsymmetry(ps_id(class_member(i,1))))
        enddo
        write(*,*)
        !
        call printer%printer_print_bar(nclasses=nclasses)
        !
    end subroutine printer_print_chartab_header

    function       printer_chi2symb(printer,chi) result(str)
        !
        implicit none
        !
        class(am_class_chartab_printer), intent(inout) :: printer
        complex(dp), intent(in) :: chi !  value o
        character(:), allocatable :: str
        !
        integer :: kk, k_exp
        complex(dp) :: s_exp
        real(dp) :: Zr,Zi
        logical :: strmatch
        !
        strmatch = .false.
        !
        Zr= real(chi)
        Zi= aimag(chi)
        !
        if ( isint(Zr) .and. iszero(Zi) ) then
            ! no imaginary, integer real
            strmatch = .true.
            str = trim(int2char(nint(Zr)))
        elseif ( iszero(Zr) .and. isint(Zi) ) then
            ! no real, imaginary integer
            strmatch = .true.
            str = trim(int2char(nint(Zi)))//'i'
        elseif ( (.not. isint(Zr)) .and. (.not. isint(Zi)) ) then
            ! complex number
            if (printer%k.ge.1) then
            search : do kk = 1, printer%k
                s_exp = cmplx(0.0_dp,0.0_dp)
                k_exp = 0 
                do while ( .not. isequal(s_exp,cmplx(1,0,dp)) )
                    ! do a full loop. complex numbers form a cyclic abelian group.
                    ! exponentiate it until it loops back to one, the identity
                    k_exp = k_exp+1
                    s_exp = printer%s(printer%k)**k_exp
                    ! check positive
                    if ( isequal(chi,s_exp) ) then
                        strmatch = .true.
                        if (k_exp.eq.1) then
                            str = char(printer%char_start+printer%k)
                        else
                            str = char(printer%char_start+printer%k)//trim(int2char(k_exp))
                        endif
                        exit search
                    endif
                    ! check negative
                    if ( isequal(chi,-s_exp) ) then
                        strmatch = .true.
                        if (k_exp.eq.1) then
                            str = '-'//char(printer%char_start+printer%k)
                        else
                            str = '-'//char(printer%char_start+printer%k)//trim(int2char(k_exp))
                        endif
                        exit search
                    endif
                enddo
            enddo search
            endif
            ! if match is not found assign a new character to variable
            if (.not.strmatch) then
                strmatch = .true.
                printer%k = printer%k + 1
                printer%s(printer%k) = chi
                str = char(printer%char_start+printer%k)
            endif
        endif
        ! if nothing found...
        if (.not.strmatch) then
            write(str,printer%fmts(5)) real(chi)
        endif
        !
    end function   printer_chi2symb

    subroutine     printer_initialize(printer)
        !
        implicit none
        !
        class(am_class_chartab_printer), intent(out) :: printer
        !
        ! left-side headers (sum to 10)
        printer%fmts(2) = '(5x,a10)'
        printer%fmts(4) = '(5x,a6,i4)'
        ! chartab (sum to 6)
        printer%fmts(1) = '(i6)'
        printer%fmts(3) = '(a6)'
        printer%fmts(5) = '(f6.2)'
        !
        allocate(printer%s(100)) ! value of symbolic output, 30 possible options
        allocate(character(4) :: printer%str)
        printer%char_start = 96 ! 97 = a
        printer%k=0
        !
    end subroutine printer_initialize

    subroutine     printer_print_bar(printer,nclasses)
        !
        implicit none
        !
        class(am_class_chartab_printer), intent(in) :: printer
        integer, intent(in) :: nclasses
        integer :: i
        !
        write(*,printer%fmts(2),advance='no') repeat('-',10)
        do i = 1, nclasses
            write(*,printer%fmts(3),advance='no') repeat('-',6)
        enddo
        write(*,*)
    end subroutine printer_print_bar

end module am_symmetry_tables
