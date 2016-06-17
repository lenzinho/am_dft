module am_prim_cell

    use am_constants
    use am_options
    use am_unit_cell
    use am_matlab
    use am_mkl, only : det, inv
    use dispmodule

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
        type(am_class_options)   , intent(in) :: opts
        character(:), allocatable :: str(:)
        real(dp) :: uc2prim(3,3)
        integer  :: i
        !
        if (opts%verbosity.ge.1) call print_title('Primitive cell')
        !
        ! get primitive cell basis from translations which leave unit cell invariant
        pc%bas = get_primitive_basis(bas=uc%bas,tau_frac=uc%tau_frac,Z=uc%Z,prec=opts%prec)
        ! get reciprocal basis
        pc%recbas = inv(pc%bas)
        ! get basis transformation (uc fractional to pc fractional)
        uc2prim = matmul(pc%recbas,uc%bas)
        ! reduce atoms to primitive cell
        pc%tau_frac = reduce_atoms_to_primitive_cell(uc2prim=uc2prim, uc_tau_frac=uc%tau_frac, prec=opts%prec)
        ! get cartesian coordiantes
        allocate(pc%tau_cart, source=matmul(pc%bas,pc%tau_frac))
        ! get numer of atoms
        pc%natoms = size(pc%tau_frac,2)
        ! get map: pc -> pc
        allocate(pc%pc_id(pc%natoms),source=[1:pc%natoms])
        ! get map: uc -> uc
        if (.not.allocated(uc%uc_id)) stop 'uc%uc_id is not allocated'
        ! allocate(uc%uc_id(uc%natoms),source=[1:uc%natoms])
        ! get map: uc -> pc
        uc%pc_id = get_uc2pc_map(uc2prim=uc2prim, pc_natoms=pc%natoms, pc_tau_frac=pc%tau_frac, uc_natoms=uc%natoms, uc_tau_frac=uc%tau_frac, prec=opts%prec)
        ! get map: pc -> uc (same algorithm as used in irreducible pair shell)
        allocate(pc%uc_id(pc%natoms))
        pc%uc_id = 0
        do i = 1, uc%natoms
        if (pc%uc_id(uc%pc_id(i)).eq.0) then
            pc%uc_id(uc%pc_id(i)) = i
        endif
        enddo
        ! transfer Z
        allocate(pc%Z(pc%natoms))
        pc%Z = uc%Z(pc%uc_id)
        !
        ! print stdout
        if (opts%verbosity.ge.1) then
            allocate(character(2) :: str(pc%natoms))
            !
            write(*,'(a,a)') flare, 'original basis ='
            call disp_indent()
            call disp(title='fractional (primitive)', X=matmul(pc%recbas,uc%bas), style='underline', fmt='f13.8',advance='no' , trim = 'no')
            call disp(title='cartesian'             , X=uc%bas                  , style='underline', fmt='f13.8',advance='yes', trim = 'no')
            !
            write(*,'(a,a)') flare, 'primitive basis ='
            call disp_indent()
            call disp(title='fractional (input)'    , X=matmul(uc%recbas,pc%bas), style='underline', fmt='f13.8',advance='no' , trim = 'no')
            call disp(title='cartesian'             , X=pc%bas                  , style='underline', fmt='f13.8',advance='yes', trim = 'no')
            !
            write(*,'(a,a,a)') flare, 'primitive cell atoms = ', tostring(pc%natoms)
            do i = 1, pc%natoms
                str(i) = atm_symb(uc%Z(pc%uc_id(i)))
            enddo
            call disp_indent()
            call disp(title='#'                     , X=[1:pc%natoms]           , style='underline', fmt='i7'   ,advance='no' , trim = 'yes')
            call disp(title='Z'                     , X=str                     , style='underline', fmt='a3'   ,advance='no' , trim = 'no')
            call disp(title='fractional (primitive)', X=transpose(pc%tau_frac)  , style='underline', fmt='f11.8',advance='no' , trim = 'no')
            call disp(title='cartesian'             , X=transpose(pc%tau_cart)  , style='underline', fmt='f11.8',advance='yes', trim = 'no')
            !
            write(*,'(a,a)',advance='no') flare, 'atomic mapping (to primitive: uc->pc)'
            call id_print_map(uc%pc_id)
            !
            write(*,'(a,a)',advance='no') flare, 'atomic mapping (to input: pc->uc)'
            call id_print_map(pc%uc_id)
            !
        endif
        !
        contains
        function       get_primitive_basis(bas,tau_frac,Z,prec) result(pc_bas)
            !
            implicit none
            !
            real(dp), intent(in) :: bas(3,3)
            real(dp), intent(in) :: tau_frac(:,:)
            integer , intent(in) :: Z(:)
            real(dp), intent(in) :: prec
            real(dp) :: pc_bas(3,3)
            integer , allocatable :: indices(:)
            real(dp), allocatable :: T(:,:)
            integer :: nTs
            integer :: i,j,k
            !
            ! get basis translations which could serve as primitive cell vectors
            T = translations_from_basis(tau=tau_frac, Z=Z, prec=prec, flags='prim')
            T = matmul(bas,T)
            nTs = size(T,2)
            ! sort primitive vectors based on magnitude (smallest last)
            allocate(indices(nTs))
            call rank(norm2(T,1),indices)
            T=T(1:3, indices(nTs:1:-1))
            ! select three primitive vectors which yield the smallest cell volume by ...
            ! ... selecting smallest vector as first vector
            i = nTs
            ! ... selecting smallest pritimive vector which is noncollinear to the first (non-zero cross product) as second vector
            do j = nTs,1,-1
                if ( any(cross_product(T(1:3,i),T(1:3,j)) .gt. tiny) ) exit
                enddo
            if (any([i,j].eq.0)) stop 'ERROR [get_primitive_basis]: No pair of non-collinear lattice vectors found'
            ! ... selecting smallest primitive vector which produces a non-zero cell volume as third vector
            do k = nTs,1,-1
                ! regularization added to prevent division by zero
                if ( abs( det( T(1:3,[i,j,k])+eye(3)*1.0D-14 ) ) .gt. tiny ) exit
            enddo
            if (any([i,j,k].eq.0)) stop 'ERROR [get_primitive_basis]: No primitive basis found'
            !
            pc_bas = T(1:3,[i,j,k])
        end function   get_primitive_basis
        function       reduce_atoms_to_primitive_cell(uc2prim,uc_tau_frac,prec) result(pc_tau_frac)
            !
            implicit none
            !
            real(dp), intent(in) :: uc2prim(3,3)
            real(dp), intent(in) :: uc_tau_frac(:,:)
            real(dp), intent(in) :: prec
            real(dp), allocatable :: pc_tau_frac(:,:)
            !
            ! allocating space
            allocate(pc_tau_frac,source=uc_tau_frac)
            ! converting atoms to pc fractional coordinates
            pc_tau_frac = matmul(uc2prim, pc_tau_frac)
            ! reduce atoms to pc
            pc_tau_frac = modulo(pc_tau_frac+prec,1.0_dp)-prec
            ! get unique values
            pc_tau_frac = unique(pc_tau_frac,prec)
            !
        end function   reduce_atoms_to_primitive_cell
        function       get_uc2pc_map(uc2prim,pc_natoms,pc_tau_frac,uc_natoms,uc_tau_frac,prec) result(uc_pc_id)
            !
            implicit none
            !
            real(dp), intent(in) :: uc2prim(3,3)
            integer , intent(in) :: pc_natoms
            integer , intent(in) :: uc_natoms
            real(dp), intent(in) :: uc_tau_frac(:,:)
            real(dp), intent(in) :: pc_tau_frac(:,:)
            real(dp), intent(in) :: prec
            integer :: uc_pc_id(uc_natoms)
            real(dp):: uc_tau_prim_frac(3,uc_natoms)
            integer :: i,j
            !
            ! intiailize
            uc_pc_id=0
            ! get uc atoms in primitive fractional coordinates
            uc_tau_prim_frac = matmul(uc2prim,uc_tau_frac)
            ! reduce to primitive cell between [0,1)
            uc_tau_prim_frac = modulo(uc_tau_prim_frac+prec,1.0_dp)-prec
            ! make map
            do i = 1, pc_natoms
            do j = 1, uc_natoms
                if ( isequal(pc_tau_frac(:,i),uc_tau_prim_frac(:,j)) ) uc_pc_id(j) = i
            enddo
            enddo
            !
            if (any(uc_pc_id.eq.0)) stop 'Unit cell -> primitive cell map failed.'
            !
        end function   get_uc2pc_map
    end subroutine get_primitive

end module am_prim_cell