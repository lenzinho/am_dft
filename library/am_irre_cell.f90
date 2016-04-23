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
        class(am_class_atom), allocatable :: atom(:)
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
        ! determine irreducible atoms, i.e. get all atoms, in increments, which have not been already mapped onto by a space symmetry operation
        ! thus, primitive cell atom i corresponds to irreducible atom k ...
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
        !     !
        !     ! this mapping is broken.
        !     ! 
        !     ! make maps
        !     ! Note: ind(1:k) are indices primitive cell atoms correspond to irreducible cell atoms
        !     ! maps irreducible cell atom onto -> irreducible cell
        !     allocate(ic%ic_identifier,source=[1:ic%natoms]) 
        !     ! maps irreducible cell atom onto -> primitive cell
        !     allocate(ic%pc_identifier,source=ind(1:k))
        !     ! maps irreducible cell atom onto -> unit cell
        !     allocate(ic%uc_identifier,source=pc%uc_identifier(ind(1:k)))
        !     ! maps primitive cell atom onto -> irreducible cell
        !     allocate(pc%ic_identifier(pc%natoms))
        !     do i = 1, ic%natoms
        !         ! PM(uc%natoms,sg%nsyms) shows how atoms are permuted by each space symmetry operation
        !         ! PM(1,:) shows all atoms onto which atom 1 is mapped by all space symmetry operations
        !         do j = 1, pc%natoms
        !         do k = 1, sg%nsyms
        !         if (PM(j,k).eq.i) pc%ic_identifier(j) = i
        !         enddo
        !         enddo
        !     enddo
        !     ! maps unit cell atom onto -> irreducible cell
        !     if (present(uc)) then
        !         allocate(uc%ic_identifier(uc%natoms))
        !         do i = 1, uc%natoms
        !             ! j = uc%pc_identifier(i) shows to which primitive atom j, unit cell atom i is associated with  
        !             ! now find... to which irreducible atom k, primitive cell atom j is associated with... k = 
        !         do j = 1, ic%natoms
        !             uc%ic_identifier(i) = pc%uc_identifier(ic%pc_identifier(j))
        !         enddo
        !         enddo
        !     endif
        !     !
        !     if (opts%verbosity.ge.1) then
        !         write(*,'(a5,a)',advance='no') ' ... ', 'atomic mapping (irreducible to original cell)'
        !         do i = 1, ic%natoms
        !             if (modulo(i,10).eq.1) then
        !                 write(*,*)
        !                 write(*,'(5x)',advance='no')
        !             endif
        !             write(*,'(a8)',advance='no') trim(int2char(i))//'->'//trim(int2char(ic%uc_identifier(i)))
        !         enddo
        !         write(*,*)
        !         !
        !         write(*,'(a5,a)',advance='no') ' ... ', 'atomic mapping (irreducible to primitive cell)'
        !         do i = 1, ic%natoms
        !             if (modulo(i,10).eq.1) then
        !                 write(*,*)
        !                 write(*,'(5x)',advance='no')
        !             endif
        !             write(*,'(a8)',advance='no') trim(int2char(i))//'->'//trim(int2char(ic%pc_identifier(i)))
        !         enddo
        !         write(*,*)
        !         !
        !         write(*,'(a5,a)',advance='no') ' ... ', 'atomic mapping (primitive to irreducible cell)'
        !         do i = 1, pc%natoms
        !             if (modulo(i,10).eq.1) then
        !                 write(*,*)
        !                 write(*,'(5x)',advance='no')
        !             endif
        !             write(*,'(a8)',advance='no') trim(int2char(i))//'->'//trim(int2char(pc%ic_identifier(i)))
        !         enddo
        !         write(*,*)
        !         !
        !         if (present(uc)) then
        !         write(*,'(a5,a)',advance='no') ' ... ', 'atomic mapping (original to irreducible cell)'
        !         do i = 1, uc%natoms
        !             if (modulo(i,10).eq.1) then
        !                 write(*,*)
        !                 write(*,'(5x)',advance='no')
        !             endif
        !             write(*,'(a8)',advance='no') trim(int2char(i))//'->'//trim(int2char(uc%ic_identifier(i)))
        !         enddo
        !         write(*,*)
        !         endif
        !         !
        !     endif
        !     
    end subroutine get_irreducible

end module am_irre_cell