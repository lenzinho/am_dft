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

    type, public :: am_class_character_table
        complex(dp), allocatable :: chartab(:,:)    ! character table ( irreps * classes )
        character(:),allocatable :: irrep_label(:)  ! irrep label
    end type am_class_character_table

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
            class_member = id_member(class_id)
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
        integer :: i, j, k
        integer :: nsyms
        integer :: nirreps  ! just to make things clearer; always equal to nclasses
        integer, allocatable :: H(:,:,:) ! class multiplication 
        integer, allocatable :: irrep_dim(:) ! sym dimensions
        integer, allocatable :: indices(:)
        real(dp) :: wrk
        !
        nsyms = size(multab,2)
        !
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
        !
        ! chartab(nirreps,nclasses) determine eigenvector which simultaneously diagonalizes i Hjk(:,:) matrices
        chartab = eigenanalysis(cmplx(H,0,dp))
        !
        ! get dimensions of representations using Eq. 4.53, p 91 of Wooten
        ! Note: The absolute square is the result of the conjugate multiplication
        nirreps = nclasses
        allocate(irrep_dim(nirreps))
        !
        do i = 1, nirreps
            ! sum on classes
            wrk = 0
            do j = 1, nclasses
                wrk = wrk + abs(chartab(i,j))**2/real(class_nelements(j),dp)
            enddo
            ! compute irrep dimension
            irrep_dim(i) = nint((nsyms/wrk)**0.5_dp)
        enddo
        !
        ! convert class constant chartab into character chartab
        do i = 1, nirreps
        do j = 1, nclasses
            chartab(i,j) = irrep_dim(i)/real(class_nelements(j),dp)*chartab(i,j)
        enddo
        enddo
        !
        ! initialize indices sort sorting
        allocate(indices(nirreps))
        !
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
        !
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
        !
        !
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
            !
            ! determines characters which simultaneously diagonalizes i square matrices A(:,:,i)
            ! only possible if A(:,:,i) commute with each other. Assumes all dimensions of A are the same.
            ! 
            !
            implicit none
            !
            complex(dp), intent(in)  :: A(:,:,:)
            complex(dp), allocatable :: V(:,:) ! these two should match in this particular situation (character chartab determination from class multiplication )
            complex(dp), allocatable :: D(:,:) ! these two should match in this particular situation (character chartab determination from class multiplication )
            complex(dp), allocatable :: M(:,:)
            integer :: n, i
            !
            n = size(A,3)
            !
            if (n.eq.0) stop 'eigenanalysis: n = 0'
            !
            ! matrix pencil approach to lift the degeneracy of eigenvalues
            allocate(M(n,n))
            M = 0
            do i = 1,n
                M = M + rand() * A(:,:,i)
            enddo
            !
            ! get left and right eigenvectors: transpose(VL)*A*VR = D
            call am_zgeev(A=M,VR=V)
            !
            ! check that matrices have been diagonalized properly and save eigenvalues
            allocate(d(n,n))
            do i = 1, n
                !
                M = matmul(inv(V),matmul(A(:,:,i),V))
                !
                if (any( abs(diag(diag(M))-M).gt. (tiny*n**2) )) then
                    ! tiny may be too small a criterion here... should probably use a value which scales with the size of the matrix.
                    stop 'Unable to perform simultaneous matrix diagonalization. Check that whether they commute.'
                endif
                !
                ! save diagonal elements as row vectors
                D(:,i) = diag( M )
                !
            enddo
            !
            ! should enforce check here...
            !
            ! ! normalize each eigencolumn to the element of the column (that is, the first element, i.e. row, corresponds to the class containing the identity element)
            ! do i = 1, n
            !     V(:,i)=V(:,i)
            ! enddo
            ! V = transpose(V)
            !
            ! at this point V and D should be identical. The elements are to within +/- sign.
            ! D is more robust 
            !
        end function   eigenanalysis
    end function   get_chartab

    function       get_muliken(chartab,nclasses,class_nelements,class_member,ps_id) result(irrep_label)
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
    end function   get_muliken

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