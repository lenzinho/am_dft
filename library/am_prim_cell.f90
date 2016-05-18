module am_prim_cell

    use am_constants
    use am_stdout
    use am_options
    use am_unit_cell
    use am_matlab
    use am_mkl , only : det, inv

    implicit none

    private

    type, public, extends(am_class_unit_cell) :: am_class_prim_cell
    contains
        procedure :: get_primitive
    end type am_class_prim_cell

contains

    subroutine     get_primitive(pc,uc,opts)
        !
        use am_rank_and_sort, only : rank
        !
        implicit none
        !
        class(am_class_prim_cell), intent(inout) :: pc
        class(am_class_unit_cell), intent(inout) :: uc
        type(am_class_options), intent(in) :: opts
        integer , allocatable :: indices(:)
        real(dp), allocatable :: T(:,:)
        real(dp) :: uc2prim(3,3)
        logical  :: isfirst
        integer  :: nTs
        integer  :: i, j, k 
        !
        if (opts%verbosity.ge.1) call am_print_title('Reducing to primitive cell')
        !
        ! get basis translations which could serve as primitive cell vectors
        T = translations_from_basis(tau=uc%tau, Z=uc%Z, prec=opts%prec, flags='prim')
        T = matmul(uc%bas,T)
        nTs = size(T,2)
        if (opts%verbosity.ge.2) call am_print("lattice translations found",nTs," ... ")
        !
        ! sort primitive vectors based on magnitude (smallest last)
        allocate(indices(nTs))
        call rank(norm2(T,1),indices)
        T=T(1:3, indices(nTs:1:-1) )
        if (opts%verbosity.ge.2) then
            call am_print_two_matrices_side_by_side(name='primitive lattice vector candidates',&
                Atitle='fractional (original)',A=transpose(matmul(inv(uc%bas),T)),&
                Btitle='cartesian' ,B=transpose(T),&
                iopt_emph=' ... ',iopt_teaser=.true.)
        endif
        !
        ! select three primitive vectors which yield the smallest cell volume by ...
        ! ... selecting smallest vector as first vector
        i = nTs
        ! ... selecting smallest pritimive vector which is noncollinear to the first (non-zero cross product) as second vector
        do j = nTs,1,-1
            if ( any(cross_product(T(1:3,i),T(1:3,j)) .gt. tiny) ) exit
        enddo
        if (any([i,j].eq.0)) then
            call am_print('ERROR','No pair of non-collinear lattice vectors found',flags='E')
            stop
        endif
        ! ... selecting smallest primitive vector which produces a non-zero cell volume as third vector
        do k = nTs,1,-1
            ! adding regularization to prevent division by zero
            if ( abs( det( T(1:3,[i,j,k])+eye(3)*1.0D-14 ) ) .gt. tiny ) exit
        enddo
        if (any([i,j,k].eq.0)) then
            call am_print('ERROR','No primitive basis found',flags='E')
            stop
        endif
        pc%bas = T(1:3,[i,j,k])
        ! output to stdout
        if (opts%verbosity.ge.2) call am_print("volume",det(pc%bas)," ... ")
        ! output to stdout
        if (opts%verbosity.ge.1) then
            call am_print_two_matrices_side_by_side(name='original basis',&
                Atitle='fractional (primitive)',A=matmul(inv(pc%bas),uc%bas),&
                Btitle='cartesian'             ,B=uc%bas,&
                iopt_emph=' ... ',iopt_teaser=.true.)
            call am_print_two_matrices_side_by_side(name='primitive basis',&
                Atitle='fractional (original)' ,A=matmul(inv(uc%bas),pc%bas),&
                Btitle='cartesian'             ,B=pc%bas,&
                iopt_emph=' ... ',iopt_teaser=.true.)
        endif
        !
        !
        ! reduce atoms to primitive cell by ...
        if (opts%verbosity.ge.2) call am_print('atoms in original cell',uc%natoms,' ... ')
        ! ... allocating enough space
        allocate(pc%tau,source=uc%tau)
        ! ... converting all atoms to primitive fractional coordinates
        uc2prim = matmul(inv(pc%bas),uc%bas)
        do i = 1, uc%natoms
            pc%tau(1:3,i) = matmul(uc2prim, pc%tau(1:3,i))
        enddo
        ! ... reducing all atoms to primitive cell
        do i = 1,uc%natoms
            pc%tau(1:3,i) = modulo(pc%tau(1:3,i)+opts%prec,1.0_dp)-opts%prec
        enddo
        ! ... getting unique values
        pc%tau = unique(pc%tau,opts%prec)
        pc%natoms = size(pc%tau,2)
        if (opts%verbosity.ge.1) call am_print('primitive cell atoms',pc%natoms,' ... ')
        if (opts%verbosity.ge.1) then
            call am_print_two_matrices_side_by_side(name='reduced atomic basis',&
                Atitle='fractional (primitive)',A=transpose(pc%tau),&
                Btitle='cartesian'             ,B=transpose(matmul(pc%bas,pc%tau)),&
                iopt_emph=' ... ',iopt_teaser=.true.)
        endif
        !
        !
        ! transfer atomic species by comparing atomic coordinates
        allocate(pc%pc_id(pc%natoms),source=[1:pc%natoms])
        allocate(pc%uc_id(pc%natoms)); pc%uc_id = 0
        allocate(uc%pc_id(uc%natoms)); uc%pc_id=0
        if (.not.allocated(uc%uc_id)) then
            allocate(uc%uc_id,source=[1:uc%natoms])
        endif
        !
        allocate(pc%Z(pc%natoms))
        do i = 1, pc%natoms
            pc%Z(i) = 0
            isfirst = .true.
            search_for_atom : do j = 1, uc%natoms
                if ( all(abs(pc%tau(1:3,i)-modulo( matmul(uc2prim,uc%tau(1:3,j))+opts%prec,1.0_dp)+opts%prec).lt.opts%prec) ) then
                    uc%pc_id(j) = i
                    pc%Z(i) = uc%Z(j)
                    !
                    if (isfirst) then
                        isfirst = .false.
                        pc%uc_id(i) = j
                    endif
                endif
            enddo search_for_atom
            if (pc%Z(i).eq.0) then
                call am_print('ERROR','Unable to identify atom in reduced cell',flags='E')
                call am_print('atom #',i)
                call am_print('atomic coordinate of unidentified atoms',pc%tau(1:3,i))
                call am_print('atom types',uc%Z(:))
                call am_print('atom original',transpose(uc%tau))
                call am_print('atom reduced',transpose(pc%tau))
                stop
            endif
        enddo
        !
        if (opts%verbosity.ge.1) then
            write(*,'(a5,a)',advance='no') ' ... ', 'atomic mapping (to primitive cell)'
            do i = 1, uc%natoms
                if (modulo(i,11).eq.1) then
                    write(*,*)
                    write(*,'(5x)',advance='no')
                endif
                write(*,'(a8)',advance='no') trim(int2char(i))//'->'//trim(int2char(uc%pc_id(i)))
            enddo
            write(*,*)
            !
            write(*,'(a5,a)',advance='no') ' ... ', 'atomic mapping (from primitive cell)'
            do i = 1, pc%natoms
                if (modulo(i,11).eq.1) then
                    write(*,*)
                    write(*,'(5x)',advance='no')
                endif
                write(*,'(a8)',advance='no') trim(int2char(i))//'->'//trim(int2char(pc%uc_id(i)))
            enddo
            write(*,*)
        endif
        !
    end subroutine get_primitive

end module am_prim_cell