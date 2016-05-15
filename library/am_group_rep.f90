module am_group_rep

    use am_constants
    use am_stdout
    use am_mkl
    use am_options

    implicit none

    public

    type am_class_rep
        integer :: nsyms ! dimensions of reps
        integer :: ndims ! dimensions of reps
        real(dp), allocatable :: sym(:,:,:)
    end type am_class_rep

    type, extends(am_class_rep) :: am_class_rep_reg
	    contains
        procedure :: create => create_regular
    end type am_class_rep_reg

    type, extends(am_class_rep) :: am_class_rep_perm
    	!
    	! basis of permutation rep are atomic positions, thus:
    	! ndims = natoms
    	!
    	integer, allocatable :: PM(:,:) ! permutation map
    	!
		contains
	    	!
	        procedure :: create => create_permutation
       		!
    end type am_class_rep_perm

	contains

	! procedures which create representations

	subroutine     create_regular(rr,R,T)
		!
		implicit none
		!
		class(am_class_rep_reg), intent(out) :: rr
        real(dp), intent(in) :: R(:,:,:)
        real(dp), intent(in), optional :: T(:,:)
		!
		if (present(T)) then
			rr%sym = rep_regular(sym=get_seitz(R=R,T=T),flags='seitz')
		else
			rr%sym = rep_regular(sym=R,flags='')
		endif
		!
		rr%ndims = size(rr%sym,1)
		rr%nsyms = size(rr%sym,3)
		!
	end subroutine create_regular

	subroutine     create_permutation(pr,R,T,tau,sym_prec)
		!
		! pr%create_perm(R=sg%R,T=sg%R,tau=uc%tau,sym_prec=opts%sym_prec)
		!
		implicit none
		!
		class(am_class_rep_perm), intent(out) :: pr
        real(dp), intent(in) :: R(:,:,:)
        real(dp), intent(in), optional :: T(:,:)
        real(dp), intent(in) :: tau(:,:)
        real(dp), intent(in) :: sym_prec
        integer :: i,j
		!
		if (present(T)) then
			pr%sym = rep_permutation(seitz=get_seitz(R=R,T=T),tau=tau,flags='',sym_prec=sym_prec)
		else
			pr%sym = rep_permutation(seitz=get_seitz(R=R),tau=tau,flags='relax_pbc',sym_prec=sym_prec)
		endif
        !
        pr%ndims = size(pr%sym,1)
        pr%nsyms = size(pr%sym,3)
        !
        ! check that each column and row of the sym sums to 1; i.e. sym(:,:,i) is orthonormal for all i.
        ! this check is equivalent to verifying whther all atoms have been permute onto another other atom
        do i = 1, pr%nsyms
        do j = 1, pr%ndims
            !
            if (sum(pr%sym(:,j,i)).ne.1) then
               call am_print('ERROR','Permutation matrix representation has a column which does not sum to 1.')
               call am_print('i',i)
               call am_print_sparse('spy(P_i)',pr%sym(:,:,i))
               call am_print('pr%sym',pr%sym(:,:,i))
               stop
            endif
            !
            if (sum(pr%sym(j,:,i)).ne.1) then
               call am_print('ERROR','Permutation matrix representation has a row which does not sum to 1.')
               call am_print('i',i)
               call am_print_sparse('spy(P_i)',pr%sym(:,:,i))
               call am_print('pr%sym',pr%sym(:,:,i))
               stop
            endif
            !
        enddo
        enddo
		!
		! get permutation map
		pr%PM = get_map(pr=pr)
		!
		contains
		    function       get_map(pr) result(PM)
		        !
		        ! PM(uc%natoms,sg%nsyms) permutation map; shows how atoms are permuted by each space symmetry operation
		        !
		        implicit none
		        !
		        class(am_class_rep_perm), intent(in) :: pr
		        integer, allocatable :: PM(:,:)
		        integer, allocatable :: ind(:)
		        integer :: i
		        !
		        ! allocate
		        allocate(PM(pr%ndims,pr%nsyms))
		        ! initialize
		        PM=0.0_dp
		        ! get indices
		        ind = [1:pr%ndims]
		        ! construct map
		        do i = 1, pr%nsyms
		            PM(:,i) = matmul(pr%sym(:,:,i),ind)
		        enddo
		        !
		    end function   get_map
	end subroutine create_permutation

    ! engines that create representations

    function       rep_regular(sym,flags) result(reg_rep)
        !>
        !> Generates regular representation for each operator
        !>
        !> For details, see page 73-74, especially Eqs. 4.18 Symmetry and Condensed Matter Physics: A Computational
        !> Approach. 1 edition. Cambridge, UK; New York: Cambridge  University Press, 2008.
        !>
        !> Requires identity to be the first element of group. regular representation is built from multiplication table
        !> by making sure that the identity appears along the diagonal. this sometimes causes the row to be permuted with
        !> respect to the columns. so, first find the permutation necessary to bring the identities to the center and
        !> then build regular representations from this new multiplication table.
        !>
        implicit none
        !
        real(dp), intent(in) :: sym(:,:,:)
        character(*), intent(in) :: flags ! can be nothing flags='' at call
        integer, allocatable :: reg_rep(:,:,:)
        integer, allocatable :: cayley_table(:,:)
        integer :: n
        integer :: i
        !
        ! obtain cayley table
        cayley_table = get_cayley_table(sym=sym,flags=flags//'reg')
        !
        ! construct regular representation from cayley table
        n = size(sym,3)
        allocate(reg_rep(n,n,n))
        reg_rep = 0
        do i = 1, n
            where (cayley_table.eq.i) reg_rep(:,:,i) = 1
        enddo
        !
    end function   rep_regular

    function       rep_permutation(seitz,tau,flags,sym_prec) result(sym)
        !
        ! find permutation representation; i.e. which atoms are connected by space symmetry oprations R, T.
        ! also works to find which kpoint are connected by point group operations
        !
        implicit none
        !
        real(dp), intent(in) :: seitz(:,:,:)
        real(dp), intent(in) :: tau(:,:)
        character(*), intent(in) :: flags
        real(dp), intent(in) :: sym_prec
        integer, allocatable :: sym(:,:,:)
        real(dp) :: tau_rot(3)
        integer :: i,j,k
        integer :: nsyms
        integer :: ntaus ! number of tau points tau(1:3,ntaus)
        !
        ntaus = size(tau,2)
        !
        nsyms = size(seitz,3)
        !
        allocate(sym(ntaus,ntaus,nsyms))
        sym = 0
        !
        do i = 1, nsyms
            ! determine the permutations of atomic indicies which results from each space symmetry operation
            do j = 1,ntaus
                ! apply rotational component
                tau_rot = matmul(seitz(1:3,1:3,i),tau(:,j))
                ! apply translational component
                tau_rot = tau_rot + seitz(1:3,4,i)
                ! reduce rotated+translated point to unit cell
                if (index(flags,'relax_pbc').eq.0) then
                    tau_rot = modulo(tau_rot+sym_prec,1.0_dp)-sym_prec
                endif
                ! find matching atom
                search : do k = 1, ntaus
                    if (all(abs(tau_rot-tau(:,k)).lt.sym_prec)) then
                    sym(j,k,i) = 1
                    exit search ! break loop
                    endif
                enddo search
            enddo
        enddo
        !
    end function   rep_permutation

    ! procedures which operate on matrix-rep symmetry elements

    subroutine     get_irreps(rr)
    	!
    	implicit none
    	!
    	class(am_class_rep) :: rr ! any type of representation... regular, permutation, except for seitz
        !
    end subroutine get_irreps

    function       get_order_sym(sym) result(i)
        !
        implicit none
        !
        real(dp), intent(in)  :: sym(:,:)
        real(dp), allocatable :: identity(:,:)
        real(dp), allocatable :: sym_pow(:,:)
        integer :: i
        !
        allocate(identity,mold=sym)
        identity = eye(size(sym,1))
        !
        allocate(sym_pow,source=sym)
        !
        i = 0
        do while ( all(abs(identity-sym_pow).lt.tiny) .or. (i.eq.50) )
            i=i+1
            sym_pow = matmul(sym,sym_pow)
        enddo
        !
        if (i.eq.50) stop 'Error: unable to determine symmetry order.'
        !
    end function   get_order_sym

    function       get_cayley_table(sym,flags) result(cayley_table)
        !
        ! flag options:
        !  reg    -  sorts rows to put identity along diagonals (useful for constructing regular representation)
        !  seitz  -  reduces translational part to primitive cell
        !
        implicit none
        !
        real(dp), intent(in) :: sym(:,:,:) ! list of 2D reps...
        character(*), intent(in), optional :: flags
        integer, allocatable :: cayley_table(:,:)
        integer, allocatable :: sortmat(:,:)
        integer, allocatable :: indices(:)
        real(dp) :: W(size(sym,1),size(sym,2)) ! workspace
        integer  :: n
        integer  :: i, j
        !
        ! before doing anything, confirm that first element is the identity
        if (any(abs(sym(:,:,1)-eye(size(sym,1))).gt.tiny)) then
            call am_print('ERROR','First element of sym is not the identity.')
            call am_print('first element of sym',sym(:,:,1))
            stop
        endif
        !
        n=size(sym,3)
        !
        allocate(cayley_table(n,n))
        !
        do i = 1, n
        do j = 1, n
            !
            W = matmul(sym(:,:,i),sym(:,:,j))
            !
            ! if 'seitz' flagged, reduce translational part to primitive cell 
            if (present(flags)) then
            if (index(flags,'seitz').ne.0) then
                W(1:3,4) = modulo(W(1:3,4)+tiny,1.0_dp)-tiny
            endif
            endif
            !
            cayley_table(i,j) = get_matching_element_index_in_list(list=sym,elem=W)
        enddo
        enddo
        !
        ! if 'reg' flagged, re-order rows so that identities are along diagonal
        if (present(flags)) then
        if (index(flags,'reg').ne.0) then
            !
            allocate(sortmat(n,n))
            allocate(indices(n))
            indices = [1:n]
            sortmat = 0
            where (cayley_table.eq.1) sortmat(:,:) = 1
            indices = matmul(sortmat(:,:),indices)
            cayley_table = cayley_table(indices,:)
        endif
        endif
        !
        ! quick consitency check. not comprehensive
        do i = 1, n
            if (sum(cayley_table(:,i)).ne.sum(cayley_table(:,1))) then
                call am_print('ERROR','Cayley table sums along rows/columns are not consistent.')
                stop
            endif
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
    end function   get_cayley_table

    ! procedures which utilize cayley table

    pure function  get_inverse_indices(cayley_table) result(ind)
        ! returns the ind of the inverse elements given the cayley table.
        ! requires identity as first element
        implicit none
        !
        integer, intent(in) :: cayley_table(:,:)
        integer, allocatable :: ind(:)
        integer, allocatable :: P(:,:)
        !
        ind = [1:size(cayley_table,1)]
        !
        allocate(P,source=cayley_table)
        where (P.ne.1) P = 0
        !
        ind = matmul(P,ind)
        !
    end function   get_inverse_indices

    function       get_conjugacy_classes(sym,flags) result(cc_id)
        !
        ! for AX = XB, if elements A and B are conjugate pairs for some other element X in the group, then they are in the same class
        !
        use am_rank_and_sort
        !
        implicit none
        !
        real(dp), intent(in) :: sym(:,:,:)
        character(*), intent(in), optional :: flags
        integer, allocatable :: cayley_table(:,:) ! cayley table
        integer, allocatable :: sinv(:) ! sinv(nsyms) index of the inverse of symmetry element i
        integer, allocatable :: cc_id(:)
        integer, allocatable :: celem(:)
        ! integer, allocatable :: class_member(:,:) ! used to sort based on determinant
        ! integer, allocatable :: indices(:)        ! used to sort based on determinant
        ! integer, allocatable :: relabel_indices(:)! used to sort based on determinant
        ! real(dp),allocatable :: ccdet(:)          ! used to sort based on determinant
        integer :: i, j, k
        integer :: nsyms
        integer :: nclasses
        !
        nsyms = size(sym,3)
        !
        ! get cayley table
        if (present(flags)) then
            cayley_table = get_cayley_table(sym=sym,flags=flags)
        else
            cayley_table = get_cayley_table(sym=sym)
        endif
        !
        ! allocate space for conjugacy class
        allocate(cc_id(nsyms))
        cc_id = 0
        !
        ! get element inverse
        sinv = get_inverse_indices(cayley_table)
        !
        ! allocate space for conjugate elements
        allocate(celem(nsyms))
        !
        ! determine conjugacy classes
        k = 0
        do i = 1, nsyms
        if (cc_id(i).eq.0) then
            k=k+1
            ! conjugate each element with all other group elements
            ! A = X(j) * B * X(j)^-1
            do j = 1, nsyms
                celem(j) = cayley_table(j,cayley_table(i,sinv(j)))
            enddo
            ! for each subgroup element created by conjugation find the corresponding index of the element in the group
            ! in order to save the class cc_id number
            do j = 1, nsyms
                cc_id( celem(j) ) = k
            enddo
            !
        endif
        enddo
        nclasses = k
        !
        ! relabel classes based on number of elements in each class
        call relabel_based_on_occurances(cc_id)
        !
        ! relabel classes based on det of class representative
        ! class_member = member(cc_id)
        ! allocate(ccdet(nclasses)) 
        ! do i = 1,nclasses
        !     ccdet(i) = det( sym(:,:,class_member(i,1)) )
        ! enddo
        ! !
        ! allocate(indices(nclasses))
        ! call rank(ccdet,indices)
        ! !
        ! allocate(relabel_indices(nclasses))
        ! do i = 1, nclasses
        !     relabel_indices(indices(i)) = i
        ! enddo
        ! do i = 1, nsyms
        !     cc_id(i) = relabel_indices( cc_id(i) )
        ! enddo
        !
        ! make sure the identity is in the first class
        ! cc_id(i=1) is the class of the first element (the identity); swap it's location with whaterver elements are in the first class
        where (cc_id.eq.1) cc_id = cc_id(1)
        cc_id(1)=1
        !
        ! check that classes are disjoint and complete
        if (any(cc_id.eq.0)) then
            call am_print('ERROR','Not every element in the group has been asigned a conjugacy class.',flags='E')
            stop
        endif
        !
        contains
        subroutine     relabel_based_on_occurances(list)
            !
            use am_rank_and_sort
            !
            implicit none
            !
            integer , intent(inout) :: list(:)
            integer , allocatable :: A_sorted(:)
            integer , allocatable :: occurances(:) 
            integer , allocatable :: reverse_sort(:)
            integer , allocatable :: sorted_indices(:)
            integer , allocatable :: list_relabled(:)
            integer :: nsyms, i, j
            !
            nsyms = size(list,1)
            !
            allocate(occurances(nsyms))
            do i = 1, nsyms
                occurances(i) = count(list(i).eq.list)
            enddo
            !
            ! this quick and dirty procedure lifts degeneracies
            !
            allocate(sorted_indices(nsyms))
            allocate(A_sorted(nsyms))
            call rank((1+maxval(list))*occurances+list,sorted_indices)
            A_sorted = list(sorted_indices)
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
            !
            ! return everything to the original order at call
            !
            allocate(reverse_sort(nsyms))
            reverse_sort(sorted_indices)=[1:nsyms]
            list=list_relabled(reverse_sort)
        end subroutine relabel_based_on_occurances
    end function   get_conjugacy_classes

    function       get_coset(g_id,H_id,cayley_table,flags) result(ind)
        !
        ! In mathematics, if G is a group, and H is a subgroup of G, and g is an element of G, then:
        ! gH = { gh : h an element of H } is the left coset of H in G with respect to g, and
        ! Hg = { hg : h an element of H } is the right coset of H in G with respect to g.
        ! Wikipedia.
        !
        implicit none
        !
        integer, intent(in) :: g_id
        integer, intent(in) :: H_id(:)
        integer, intent(in) :: cayley_table(:,:)
        character(*), intent(in) :: flags ! left/right
        integer, allocatable :: ind(:)
        integer :: i, m
        !
        m = size(H_id)
        !
        allocate(ind(m))
        !
        do i = 1, m
            ! left/right coset, coonjucate 
            if     (index(flags,'right')) then; ind(i) = cayley_table(H_id(i),g_id)
            elseif (index(flags,'left' )) then; ind(i) = cayley_table(g_id,H_id(i))
            else
                call am_print('ERROR','Flag '//trim(flags)//' not valid.',flags='E')
                stop
            endif
        enddo
        !
    end function   get_coset

    ! functions which operate on conjugate classes ids

    function       member(id) result(class_member)
        !
        implicit none
        !
        integer, allocatable, intent(in) :: id(:)
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
        class_nelements = nelements(id)
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
    end function   member

    function       nelements(id) result(class_nelements)
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
    end function   nelements

    ! other

    function       get_seitz(R,T) result(sym)
        !
        implicit none
        !
        real(dp), intent(in) :: R(:,:,:)
        real(dp), intent(in), optional :: T(:,:)
        real(dp), allocatable :: sym(:,:,:)
        integer :: i, n
        !
        n = size(R,3)
        !
        allocate(sym(4,4,n))
        sym = 0
        !
        do i = 1, n
            sym(1:3,1:3,i) = R(:,:,i)
            sym(  4,  4,:) = 1.0_dp
        enddo
        !
        if (present(T)) then
        do i = 1, n
            sym(1:3,4,i) = T(:,i)
        enddo
        endif
        !
    end function   get_seitz

end module am_group_rep









