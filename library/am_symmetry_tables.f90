module am_symmetry_tables

	use am_constants
	use am_matlab
	use am_mkl
	use am_stdout

	implicit none

	public

	type, public :: am_class_conjugacy_class
		integer :: nclasses							! nclasses number of classes
		integer, allocatable :: id(:)				! class_id(nsyms) identify conjugacy classes
		integer, allocatable :: nelements(:)		! class_nelements(nclasses) get number of elements in each class
		integer, allocatable :: member(:,:)			! members(nclass,maxval(class_nelements)) record indicies of each class_member element for each class 
		integer, allocatable :: representative(:)	! class_nelements(nclasses) get number of elements in each class
		contains
	end type am_class_conjugacy_class

	type, public :: am_class_multiplication_table
		integer, allocatable :: multab(:,:)
		integer, allocatable :: inv_id(:)
		integer, allocatable :: order(:)
		contains
	end type am_class_multiplication_table

contains

	! multiplication table

    function       get_multab(sym,flags) result(multab)
        !
        ! flag options:
        !  sort    - sorts rows to put identity along diagonals (useful for constructing regular representation)
        !  seitz   - reduces seitz(1:3,4) to between 0 and 1
        !  prog    - shows progress bar
        !
        implicit none
        !
        real(dp), intent(in) :: sym(:,:,:) ! list of 2D reps...
        character(*), intent(in) :: flags
        integer, allocatable :: multab(:,:)
        integer, allocatable :: sortmat(:,:)
        integer, allocatable :: indices(:)
        real(dp) :: W(size(sym,1),size(sym,2)) ! workspace
        integer :: ref ! used for checking for closure
        integer :: n
        integer :: i, j ! loop variabls
        integer :: k, kmax ! for progress bar
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
        allocate(multab(n,n))
        !
        k=0
        kmax=n**2
        do i = 1, n
        do j = 1, n
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
            ! get matching lment
            multab(i,j) = get_matching_element_index_in_list(list=sym,elem=W)
        enddo
        enddo
        !
        ! if 'sort' flagged, re-order rows so that identities are along diagonal
        if (index(flags,'sort').ne.0) then
            !
            allocate(sortmat(n,n))
            allocate(indices(n))
            indices = [1:n]
            sortmat = 0
            where (multab.eq.1) sortmat(:,:) = 1
            indices = matmul(sortmat(:,:),indices)
            multab = multab(indices,:)
        endif
        !
        ! quick consitency check for closure. not comprehensive
        ref = sum(multab(:,1))
        do i = 1, n
            if (sum(multab(:,i)).ne.ref) then
                call am_print('ERROR','multiplication chartab sums along rows/columns are not consistent.')
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

    function       get_cyclic_order(multab) result(order)
        !
        implicit none
        !
        integer, intent(in) :: multab(:,:)
        integer, allocatable :: order(:)
        integer :: power
        integer :: i, j
        integer :: nsyms
        !
        nsyms = size(multab,2)
        !
        allocate(order(nsyms))
        !
        do i = 1, nsyms
            ! initialize
            j = 0
            power = i
            ! apply operation until identity (1) is returned
            do while (power.ne.1)
                !
                j=j+1
                power = multab(i,power)
                !
                if (j.eq.50) stop 'Error: unable to determine symmetry order.'
            enddo
            !
            order(i) = j
        enddo
        !
    end function   get_cyclic_order

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
        !
        ! allocate space for conjugacy class
        allocate(class_id(nsyms))
        class_id = 0
        !
        ! allocate space for conjugate elements
        allocate(conjugates(nsyms))
        !
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
            !
        endif
        enddo
        nclasses = k
        !
        ! relabel based on how many elements each class has
        call relabel_based_on_occurances(class_id)
        ! relabel classes based on point symmetry id of representative class element
        call relabel_based_on_ps_id(class_id,ps_id)
        !
        ! make sure the identity is in the first class
        ! class_id(i=1) is the class of the first element (the identity); swap it's location with whaterver elements are in the first class
        where (class_id.eq.1) class_id = class_id(1)
        class_id(1)=1
        !
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
            !
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
            !
            ! return everything to the original order at call
            !
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
            !
            nsyms = size(class_id)
            !
            nclasses = maxval(class_id)
            !
            allocate(inds,source=[1:nclasses])
            !
            ! identify class members 
            ! members(nclass,maxval(class_nelements))
            class_member = member(class_id)
            !
            ! for each class, identify the representative element
            allocate(s(nclasses))
            do i = 1, nclasses
                s(i) = ps_id( class_member(i,1) )
            enddo
            ! 
            ! sort classes based on representative element id
            call rank(s,inds)
            !
            !
            allocate(rinds(nclasses))
            rinds(inds) = [1:nclasses] 
            !
            ! at this point inds contains the new labeling for class_id
            ! [1:nclasses] => inds; relabel according to new labeling
            do i = 1, nsyms
                class_id(i) = rinds( class_id(i) )
            enddo
            !
        end subroutine relabel_based_on_ps_id
    end function   get_class_id

    ! identifier functions which operate on identifiers

    function       member(id) result(class_member)
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

    function       representative(id) result(class_representative)
        !
        implicit none
        !
        integer, intent(in) :: id(:)
        integer, allocatable :: class_member(:,:) ! members(nclass,maxval(class_nelements))
        integer, allocatable :: class_representative(:)
        !
        class_member = member(id)
        !
        allocate(class_representative, source=class_member(:,1))
        !
    end function   representative

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