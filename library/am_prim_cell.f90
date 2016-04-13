module am_prim_cell

    use am_constants
    use am_helpers
    use am_options
    use am_unit_cell
    use am_symmetry
    use am_mkl

    implicit none

    private

    type, public, extends(am_class_unit_cell) :: am_class_prim_cell
    contains
        procedure :: get_primitive
    end type am_class_prim_cell

contains

    subroutine     get_primitive(prim,uc,sg,opts)
        !
        use am_rank_and_sort, only : rank
        !
        implicit none
        !
        class(am_class_prim_cell), intent(inout) :: prim
        class(am_class_unit_cell), intent(in) :: uc
        type(am_class_options), intent(in) :: opts
        type(am_class_symmetry), optional, intent(in) :: sg ! not used for anything
        integer , allocatable :: indices(:)
        real(dp), allocatable :: T(:,:)
        real(dp) :: uc2prim(3,3) 
        real(dp) :: tau_ref(3) ! currently not in use, this shifts the basis, putting one atom at the origin
        integer  :: nTs
        integer  :: i, j, k 
        !
        if ( opts%verbosity .ge. 1) call am_print_title('Reducing to primitive cell')
        !
        !
        ! get basis translations which could serve as primitive cell vectors
        T = translations_from_basis(tau=uc%tau,atype=uc%atype,iopt_include=trim('prim'),iopt_sym_prec=opts%sym_prec)
        T = matmul(uc%bas,T)
        nTs = size(T,2)
        if (opts%verbosity.ge.1) call am_print("possible primitive lattice translations found",nTs," ... ")
        !
        !
        ! sort primitive vectors based on magnitude (smallest last)
        allocate(indices(nTs))
        call rank(norm2(T,1),indices)
        T=T(1:3, indices(nTs:1:-1) )
        if (opts%verbosity.ge.1) then
            call am_print_two_matrices_side_by_side(name='possible primitive lattice vectors',&
                Atitle='fractional (original)',A=transpose(matmul(inv(uc%bas),T)),&
                Btitle='cartesian' ,B=transpose(T),&
                iopt_emph=' ... ',iopt_teaser=.true.)
        endif
        !
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
        prim%bas = T(1:3,[i,j,k])
        ! output to stdout
        if (opts%verbosity.ge.1) then
            call am_print("volume",det(prim%bas)," ... ")
            call am_print_two_matrices_side_by_side(name='original basis',&
                Atitle='fractional (primitive)',A=matmul(inv(prim%bas),uc%bas),&
                Btitle='cartesian'             ,B=uc%bas,&
                iopt_emph=' ... ',iopt_teaser=.true.)
            call am_print_two_matrices_side_by_side(name='primitive basis',&
                Atitle='fractional (original)' ,A=matmul(inv(uc%bas),prim%bas),&
                Btitle='cartesian'             ,B=prim%bas,&
                iopt_emph=' ... ',iopt_teaser=.true.)
        endif
        !
        !
        ! reduce atoms to primitive cell by ...
        if (opts%verbosity.ge.1) call am_print('number of atoms in original cell',uc%natoms,' ... ')
        ! ... allocating enough space
        allocate(prim%tau(3,uc%natoms))
        ! ... setting first atom at the origin (no need to do this...)
        ! tau_ref = uc%tau(1:3,1)
        tau_ref = real([0,0,0],dp)
        do i = 1, uc%natoms
            prim%tau(1:3,i) = uc%tau(1:3,i) - tau_ref
        enddo
        ! ... converting all atoms to primitive fractional coordinates
        uc2prim = matmul(inv(prim%bas),uc%bas)
        do i = 1, uc%natoms
            prim%tau(1:3,i) = matmul(uc2prim, prim%tau(1:3,i))
        enddo
        ! ... reducing all atoms to primitive cell
        do i = 1,uc%natoms
            prim%tau(1:3,i) = modulo(prim%tau(1:3,i)+opts%sym_prec,1.0_dp)-opts%sym_prec
        enddo
        ! ... getting unique values
        prim%tau = unique(prim%tau,opts%sym_prec)
        prim%natoms = size(prim%tau,2)
        if (opts%verbosity.ge.1) call am_print('number of atoms in primitive cell',prim%natoms,' ... ')
        if (opts%verbosity.ge.1) then
            call am_print_two_matrices_side_by_side(name='reduced atomic basis',&
                Atitle='fractional (primitive)',A=transpose(prim%tau),&
                Btitle='cartesian'             ,B=transpose(matmul(prim%bas,prim%tau)),&
                iopt_emph=' ... ',iopt_teaser=.true.)
        endif
        !
        !
        ! transfer atomic species by comparing atomic coordinates
        allocate(prim%atype(prim%natoms))
        do i = 1, prim%natoms
            prim%atype(i) = 0
            search_for_atom : do j = 1, uc%natoms
                if ( all(abs(prim%tau(1:3,i)-modulo( matmul(uc2prim,uc%tau(1:3,j)-tau_ref)+opts%sym_prec,1.0_dp)+opts%sym_prec).lt.opts%sym_prec) ) then
                    prim%atype(i) = uc%atype(j)
                    exit search_for_atom
                endif
            enddo search_for_atom
            if (prim%atype(i).eq.0) then
                call am_print('ERROR','Unable to identify atom in reduced cell',flags='E')
                call am_print('atom #',i)
                call am_print('coordinate of reference atom',tau_ref)
                call am_print('atomic coordinate of unidentified atoms',prim%tau(1:3,i))
                call am_print('atom types',uc%atype(:))
                call am_print('atom original',transpose(uc%tau))
                call am_print('atom reduced',transpose(prim%tau))
                stop
            endif
        enddo
        ! 
        !
        ! lastly, transfer other stuff
        allocate(character(500)::prim%symb(uc%nspecies))
        prim%symb = uc%symb
        prim%nspecies = uc%nspecies
        !
    end subroutine get_primitive

end module am_prim_cell