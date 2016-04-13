module am_irre_cell

    use am_constants
    use am_helpers
    use am_options
    use am_symmetry
    use am_unit_cell
    use am_prim_cell
    use am_mkl

    implicit none

    private

    type, public, extends(am_class_prim_cell) :: am_class_irre_cell
    contains
        procedure :: get_irreducible
    end type am_class_irre_cell

contains

    subroutine     get_irreducible(ic,pc,sg,opts)
        !
        implicit none
        !
        class(am_class_irre_cell), intent(inout) :: ic ! irreducible cell
        class(am_class_unit_cell), intent(in) :: pc
        type(am_class_symmetry), intent(in) :: sg
        type(am_class_options), intent(in) :: opts
        type(am_class_options) :: notalk
        integer, allocatable :: PM(:,:)
        logical, allocatable :: mask(:)
        integer, allocatable :: ind(:)
        integer :: i,j,k
        !
        if (opts%verbosity.ge.1) call am_print_title('Reducing to irreducible cell')
        !
        notalk = opts
        notalk%verbosity=0
        !
        if (opts%verbosity.ge.1) call am_print('number of atoms in primitive cell',pc%natoms,' ... ')
        !
        call sg%symmetry_action(uc=pc,flags='',oopt_PM=PM,opts=notalk)
        !
        allocate(mask(pc%natoms))
        allocate(ind(pc%natoms))
        mask = .true.
        !
        k=0
        do i = 1, pc%natoms
        if (mask(i)) then
            k=k+1
            ind(k)=i
            ! PM(uc%natoms,sg%nsyms) shows how atoms are permuted by each space symmetry operation
            do j = 1,sg%nsyms
                mask(PM(i,j))=.false.
            enddo
        endif
        enddo
        !
        ic%bas = pc%bas
        ic%natoms = k
        ic%nspecies = pc%nspecies
        !
        allocate(ic%symb ,source=pc%symb)
        allocate(ic%tau  ,source=pc%tau(:,ind(1:k)))
        allocate(ic%atype,source=pc%atype(ind(1:k)))
        !
        if (opts%verbosity.ge.1) call am_print('number of atoms in irreducible cell',ic%natoms,' ... ')
        if (opts%verbosity.ge.1) then
            call am_print_two_matrices_side_by_side(name='irreducible atomic basis',&
                Atitle='fractional',A=transpose(ic%tau),&
                Btitle='cartesian' ,B=transpose(matmul(ic%bas,ic%tau)),&
                iopt_emph=' ... ',iopt_teaser=.true.)
        endif
        !
    end subroutine get_irreducible

end module am_irre_cell