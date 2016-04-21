module am_irre_cell

    use am_constants
    use am_helpers
    use am_options
    use am_symmetry
    use am_unit_cell
    use am_prim_cell
    use am_mkl
    use am_atom

    implicit none

    private
    
    type, public, extends(am_class_prim_cell) :: am_class_irre_cell
        !
    contains
        procedure :: get_irreducible
    end type am_class_irre_cell

contains

    
    
    subroutine     get_irreducible(ic,pc,uc,sg,opts)
        !
        implicit none
        !
        class(am_class_irre_cell), intent(inout) :: ic ! irreducible cell
        class(am_class_unit_cell), intent(inout) :: pc
        class(am_class_unit_cell), intent(inout), optional :: uc
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
        ! PM(uc%natoms,sg%nsyms) shows how atoms are permuted by each space symmetry operation
        call sg%symmetry_action(uc=pc,flags='',oopt_PM=PM,opts=notalk)
        !
        allocate(mask(pc%natoms))
        allocate(ind(pc%natoms))
        mask = .true.
        !
        !
        ! determine irreducible atoms.
        k=0
        do i = 1, pc%natoms
        if (mask(i)) then
            k=k+1
            ind(k)=i
            do j = 1,sg%nsyms
                mask(PM(i,j))=.false.
            enddo
        endif
        enddo
        !
        ic%bas = pc%bas
        ic%natoms = k
        !
        allocate(ic%tau,source=pc%tau(:,ind(1:k)))
        allocate(ic%Z,source=pc%Z(ind(1:k)))
        !
        if (opts%verbosity.ge.1) call am_print('number of atoms in irreducible cell',ic%natoms,' ... ')
        if (opts%verbosity.ge.1) then
            call am_print_two_matrices_side_by_side(name='irreducible atomic basis',&
                Atitle='fractional',A=transpose(ic%tau),&
                Btitle='cartesian' ,B=transpose(matmul(ic%bas,ic%tau)),&
                iopt_emph=' ... ',iopt_teaser=.true.)
        endif
        !
        ! make maps
        allocate(ic%uc_identifier,source=pc%uc_identifier(ind(1:k)))
        allocate(ic%pc_identifier,source=ind(1:k))
        allocate(ic%ic_identifier,source=[1:ic%natoms])
        !
        allocate(pc%ic_identifier(pc%natoms))
        do i = 1, ic%natoms
            where (PM(i,:).eq.i) pc%ic_identifier = i
        enddo
        !
        if (present(uc)) then
            allocate(uc%ic_identifier(pc%natoms))
            do i = 1, uc%natoms
                uc%ic_identifier(i) = pc%ic_identifier(uc%pc_identifier(i))
            enddo
        endif
        !
        if (opts%verbosity.ge.1) then
            write(*,'(5x,a)',advance='no') 'atomic mapping (irreducible to original cell)'
            do i = 1, ic%natoms
                if (modulo(i,10).eq.1) then
                    write(*,*)
                    write(*,'(5x)',advance='no')
                endif
                write(*,'(a8)',advance='no') trim(int2char(i))//'->'//trim(int2char(ic%uc_identifier(i)))
            enddo
            write(*,*)
            !
            write(*,'(5x,a)',advance='no') 'atomic mapping (irreducible to primitive cell)'
            do i = 1, ic%natoms
                if (modulo(i,10).eq.1) then
                    write(*,*)
                    write(*,'(5x)',advance='no')
                endif
                write(*,'(a8)',advance='no') trim(int2char(i))//'->'//trim(int2char(ic%pc_identifier(i)))
            enddo
            write(*,*)
            !
            write(*,'(5x,a)',advance='no') 'atomic mapping (primitive to irreducible cell)'
            do i = 1, pc%natoms
                if (modulo(i,10).eq.1) then
                    write(*,*)
                    write(*,'(5x)',advance='no')
                endif
                write(*,'(a8)',advance='no') trim(int2char(i))//'->'//trim(int2char(pc%ic_identifier(i)))
            enddo
            write(*,*)
            !
            if (present(uc)) then
            write(*,'(5x,a)',advance='no') 'atomic mapping (original to irreducible cell)'
            do i = 1, uc%natoms
                if (modulo(i,10).eq.1) then
                    write(*,*)
                    write(*,'(5x)',advance='no')
                endif
                write(*,'(a8)',advance='no') trim(int2char(i))//'->'//trim(int2char(uc%ic_identifier(i)))
            enddo
            write(*,*)
            endif
            !
        endif
        !     
    end subroutine get_irreducible

end module am_irre_cell