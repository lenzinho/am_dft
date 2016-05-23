module am_irre_cell

    use am_constants
    use am_stdout
    use am_options
    use am_symmetry
    use am_unit_cell
    use am_prim_cell
    use am_mkl
    use am_atom

    implicit none

    private
    
    type, public, extends(am_class_prim_cell) :: am_class_irre_cell
        ! atom is used by tight binding
        class(am_class_atom), allocatable :: atom(:)
        !
    contains
        procedure :: get_irreducible
        procedure :: initialize_orbitals
        procedure :: print_orbital_basis
    end type am_class_irre_cell

contains
    
    subroutine     get_irreducible(ic,pc,uc,sg,opts)
        !
        implicit none
        !
        class(am_class_irre_cell), intent(out)   :: ic ! irreducible cell
        class(am_class_prim_cell), intent(inout) :: pc
        class(am_class_unit_cell), intent(inout), optional :: uc ! if present, mapping from ic <-> uc is obtained
        type(am_class_space_group),intent(in) :: sg
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
        if (opts%verbosity.ge.1) call am_print('primitive cell atoms',pc%natoms,' ... ')
        !
        ! PM(uc%natoms,sg%nsyms) shows how atoms are permuted by each space symmetry operation
        PM = permutation_map( permutation_rep(seitz=sg%sym, tau=pc%tau, flags='', prec=opts%prec) )
        !
        allocate(mask(pc%natoms))
        allocate(ind(pc%natoms))
        mask = .true.
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
        if (opts%verbosity.ge.1) call am_print('irreducible cell atoms',ic%natoms,' ... ')
        if (opts%verbosity.ge.1) then
            call am_print_two_matrices_side_by_side(name='irreducible atomic basis',&
                Atitle='fractional',A=transpose(ic%tau),&
                Btitle='cartesian' ,B=transpose(matmul(ic%bas,ic%tau)),&
            iopt_emph=' ... ',iopt_teaser=.true.)
        endif
        !
        ! <MAP>
            ! Note: ind(1:k) are indices on primitive cell atoms correspond to irreducible cell atoms
            ! maps irreducible cell atom onto -> irreducible cell
            allocate(ic%ic_id,source=[1:ic%natoms])
            ! maps irreducible cell atom onto -> primitive cell
            allocate(ic%pc_id,source=ind(1:k))
            ! maps irreducible cell atom onto -> unit cell
            allocate(ic%uc_id,source=pc%uc_id(ind(1:k)))
            ! maps primitive cell atom onto -> irreducible cell
            allocate(pc%ic_id(pc%natoms))
            pc%ic_id = 0
            do i = 1, ic%natoms
                ! PM(uc%natoms,sg%nsyms) shows how atoms are permuted by each space symmetry operation
                ! PM(1,:) shows all atoms onto which atom 1 is mapped by all space symmetry operations
                do j = 1, pc%natoms
                    search : do k = 1, sg%nsyms
                        if (PM(j,k).eq.ic%pc_id(i)) then
                            pc%ic_id(j) = i
                            exit search
                        endif
                    enddo search
                enddo
            enddo
            if (any(pc%ic_id.eq.0)) stop 'ERROR: prim->ic mapping failed.'
            ! maps (input) unit cell atom onto -> irreducible cell
            if (present(uc)) then
                allocate(uc%ic_id(uc%natoms))
                uc%ic_id = 0
                do i = 1, uc%natoms
                    ! j = uc%pc_id(i) shows to which primitive atom j, unit cell atom i is associated with
                    ! now find... to which irreducible atom k, primitive cell atom j is associated with... k =
                    uc%ic_id(i) = pc%ic_id(uc%pc_id(i))
                enddo
            endif
            !
            if (opts%verbosity.ge.1) then
                write(*,'(a5,a)',advance='no') ' ... ', 'atomic mapping (to input: irr->input)'
                do i = 1, ic%natoms
                    if (modulo(i,10).eq.1) then
                        write(*,*)
                        write(*,'(5x)',advance='no')
                    endif
                    write(*,'(a8)',advance='no') trim(int2char(i))//'->'//trim(int2char(ic%uc_id(i)))
                enddo
                write(*,*)
                !
                write(*,'(a5,a)',advance='no') ' ... ', 'atomic mapping (to primitive: irr->prim)'
                do i = 1, ic%natoms
                    if (modulo(i,10).eq.1) then
                        write(*,*)
                        write(*,'(5x)',advance='no')
                    endif
                    write(*,'(a8)',advance='no') trim(int2char(i))//'->'//trim(int2char(ic%pc_id(i)))
                enddo
                write(*,*)
                !
                write(*,'(a5,a)',advance='no') ' ... ', 'atomic mapping (from primitive: prim->irr)'
                do i = 1, pc%natoms
                    if (modulo(i,10).eq.1) then
                        write(*,*)
                        write(*,'(5x)',advance='no')
                    endif
                    write(*,'(a8)',advance='no') trim(int2char(i))//'->'//trim(int2char(pc%ic_id(i)))
                enddo
                write(*,*)
                !
                if (present(uc)) then
                    write(*,'(a5,a)',advance='no') ' ... ', 'atomic mapping (from input: input->irr)'
                    do i = 1, uc%natoms
                        if (modulo(i,10).eq.1) then
                            write(*,*)
                            write(*,'(5x)',advance='no')
                        endif
                        write(*,'(a8)',advance='no') trim(int2char(i))//'->'//trim(int2char(uc%ic_id(i)))
                    enddo
                    write(*,*)
                endif
                !
            endif
        ! </MAP>
        !
    end subroutine get_irreducible

    subroutine     initialize_orbitals(ic,opts)
        !
        ! Useful if no states on atoms are defined. Gives each atom all the possible states
        !
        ! Correspoding to the L shell (n=2) [ s p   states (1+3=4   states) ]
        ! Correspoding to the M shell (n=3) [ s p d states (1+3+5=9 states) ]
        !
        ! There will be 1 matrix element for every irreducible pair of orbitals
        !
        implicit none
        !
        class(am_class_irre_cell), intent(inout) :: ic
        type(am_class_options)   , intent(in) :: opts
        integer :: i
        !
        ! create basis on each irreducible atoms
        allocate(ic%atom(ic%natoms))
        do i = 1, ic%natoms
            call ic%atom(i)%gen_orbitals(orbital_flags='2s,2p')
        enddo
        if (opts%verbosity.ge.1) call ic%print_orbital_basis()
        !
    end subroutine initialize_orbitals

    subroutine     print_orbital_basis(ic)
        !
        implicit none
        !
        class(am_class_irre_cell) , intent(inout) :: ic
        integer :: n,l,m,s
        integer :: i,j
        !
        do i = 1, ic%natoms
            write(*,'(a5,a)',advance='no') ' ... ', 'irreducible atom '//trim(int2char(i))//' contributes '//trim(int2char(ic%atom(i)%norbitals))//' orbitals (n,l,m,s):'
            !
            do j = 1, ic%atom(i)%norbitals
                if (mod(j,6).eq.1) then
                    write(*,*)
                    write(*,'(5x)',advance='no')
                endif
                n = ic%atom(i)%orbital(1,j)
                l = ic%atom(i)%orbital(2,j)
                m = ic%atom(i)%orbital(3,j)
                s = ic%atom(i)%orbital(4,j)
                write(*,'(a2,i2,a1,i2,a1,i2,a1,i2,a2)',advance='no') ' (', n, ',', l, ',', m, ',', s, ') '
            enddo
            write(*,*)
        enddo
    end subroutine print_orbital_basis

    function       get_Vsk_ind(ic,i,alpha,j,beta,k) result(ind)
        ! returns matrix element indices given irreducible atoms (i,j), orbitals (alpha,beta), anda shell k
        ! ind = [i,ni,li,mi,si,j,nj,lj,mj,sj,k]
        implicit none
        !
        class(am_class_irre_cell), intent(in) :: ic
        integer, intent(in) :: i ! irreducible atom index
        integer, intent(in) :: j ! irreducible atom index
        integer, intent(in) :: alpha ! orbital index
        integer, intent(in) :: beta  ! orbital index
        integer, intent(in) :: k ! shell index
        integer :: ind(11)
        !
        ind(1)  = i
        ind(2)  = ic%atom(i)%orbital(1,alpha)
        ind(3)  = ic%atom(i)%orbital(2,alpha)
        ind(4)  = ic%atom(i)%orbital(3,alpha)
        ind(5)  = ic%atom(j)%orbital(4,alpha)
        ind(6)  = j
        ind(7)  = ic%atom(j)%orbital(1,beta)
        ind(8)  = ic%atom(j)%orbital(2,beta)
        ind(9)  = ic%atom(j)%orbital(3,beta)
        ind(10) = ic%atom(j)%orbital(4,beta)
        ind(11) = k
        !
        ! impose symmetry requirements
        !
        ! 1) (l,l',m) = (l,l',-m)
        ind(4)=abs(ind(4))
        ind(9)=abs(ind(9))
        !
        ! 2) (l,l',m) = (-1)^(l-m) * (l',l,m)
        ! if irreducible atoms are identical, flipped orbitals are related
        if (ind(1).eq.ind(6)) then !  i = j
            if (ind(2).gt.ind(7))  ind = flip_orbitals(ind) ! n > n'
            if (ind(5).gt.ind(10)) ind = flip_orbitals(ind) ! s > s'
            if (ind(4).gt.ind(9))  ind = flip_orbitals(ind) ! l > l'
        endif
        !
        contains
        function       flip_orbitals(ind) result(ind_flipped)
            ! returns the corresponding matrix element given the indices
            ! convert
            ! ind = [i,ni,li,mi,si, j,nj,lj,mj,sj, k]
            ! to
            ! ind = [j,nj,lj,mj,sj, i,ni,li,mi,si, k]
            implicit none
            !
            integer, intent(in) :: ind(1:11)
            integer :: ind_flipped(1:11)
            !
            ind_flipped([1,2,3,4,5, 6,7,8,9,10, 11]) = ind([6,7,8,9,10, 1,2,3,4,5, 11])
            !
        end function   flip_orbitals
    end function   get_Vsk_ind

end module am_irre_cell