module am_irre_cell

    use am_constants
    use am_matlab
    use am_stdout
    use am_options
    use am_symmetry
    use am_unit_cell
    use am_prim_cell
    use am_mkl
    use am_atom
    use dispmodule

    implicit none

    private
    
    type, public, extends(am_class_prim_cell) :: am_class_irre_cell
        ! atom is used by tight binding
        class(am_class_atom), allocatable :: atom(:)
        !
    contains
        procedure :: get_irreducible
        procedure :: initialize_orbitals
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
        type(am_class_options)   , intent(in) :: opts
        character(:), allocatable :: str(:)
        integer, allocatable :: PM(:,:)
        logical, allocatable :: mask(:)
        integer, allocatable :: ind(:)
        integer :: i,j,k
        !
        if (opts%verbosity.ge.1) call print_title('Reducing to irreducible cell')
        !
        ! get primitive basis
        ic%bas = pc%bas
        ! get reciprocal basis
        ic%recbas = pc%recbas
        ! get permutation map which shows how atoms are permuted by each space symmetry operation, PM(pc%natoms,sg%nsyms) 
        PM = permutation_map( permutation_rep(seitz=sg%seitz_frac, tau=pc%tau_frac, flags='', prec=opts%prec) )
        ! determine irreducible atoms, i.e. get all atoms, in increments, which have not been already mapped onto by a space symmetry operation
        ! thus, primitive cell atom i corresponds to irreducible atom k ...
        allocate(mask(pc%natoms))
        mask = .true.
        ! ind(1:k) are indices on primitive atoms corresponding to irreducible atoms
        allocate(ind(pc%natoms)) 
        ind  = 0
        !
        k=0
        do i = 1, pc%natoms
        if (mask(i)) then
            k=k+1
            ind(k)=i
            do j = 1, sg%nsyms
                ! mask keeps track of pc atoms which were already mapped onto irreducible atoms
                mask(PM(i,j))=.false.
            enddo
        endif
        enddo
        ! get number of atoms
        ic%natoms = k
        ! transfer irreducible atoms (frac)
        allocate(ic%tau_frac,source=pc%tau_frac(:,ind(1:k)))
        ! transfer irreducible atoms (cart)
        allocate(ic%tau_cart,source=pc%tau_cart(:,ind(1:k)))
        ! transfer Z
        allocate(ic%Z,source=pc%Z(ind(1:k)))
        ! map irreducible atom -> irreducible atom
        allocate(ic%ic_id,source=[1:ic%natoms])
        ! map irreducible atom onto -> primitive atom
        allocate(ic%pc_id,source=ind(1:k))
        ! map irreducible atom onto -> unit atom
        allocate(ic%uc_id,source=pc%uc_id(ind(1:k)))
        ! map primitive atom onto -> irreducible atom
        allocate(pc%ic_id(pc%natoms))
        pc%ic_id = 0
        do i = 1, ic%natoms
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
        if (any(pc%ic_id.eq.0)) stop 'ERROR: pc->ic mapping failed.'
        ! maps (input) unit cell atom onto -> irreducible atom
        if (present(uc)) then
            allocate(uc%ic_id(uc%natoms))
            uc%ic_id = 0
            do i = 1, uc%natoms
                uc%ic_id(i) = pc%ic_id(uc%pc_id(i))
            enddo
        endif
        ! print stdout
        if (opts%verbosity.ge.1) then
            allocate(character(4) :: str(ic%natoms))
            !
            write(*,'(a5,a,a)') ' ... ', 'input atoms = ', tostring(uc%natoms)
            write(*,'(a5,a,a)') ' ... ', 'primitive atoms = ', tostring(pc%natoms)
            write(*,'(a5,a,a)') ' ... ', 'irreducible atoms = ', tostring(ic%natoms)
            do i = 1, ic%natoms
                str(i) = atm_symb(ic%Z(i))
            enddo
            call disp(title=' '                     , X=zeros([1,1])            , style='underline', fmt='i2'   ,advance='no' , zeroas=' ' )
            call disp(title=' id '                  , X=[1:ic%natoms]           , style='underline', fmt='i7'   ,advance='no' , trim = 'yes')
            call disp(title='atom'                  , X=str                     , style='underline', fmt='a3'   ,advance='no' , trim = 'no')
            call disp(title='fractional (primitive)', X=transpose(ic%tau_frac)  , style='underline', fmt='f11.8',advance='no' , trim = 'no')
            call disp(title='cartesian'             , X=transpose(ic%tau_cart)  , style='underline', fmt='f11.8',advance='yes', trim = 'no')
            !
            write(*,'(a5,a)',advance='no') ' ... ', 'atomic mapping (to input: ic->uc)'
            call id_print_map(ic%uc_id)
            !
            write(*,'(a5,a)',advance='no') ' ... ', 'atomic mapping (to primitive: ic->pc)'
            call id_print_map(ic%pc_id)
            !
            write(*,'(a5,a)',advance='no') ' ... ', 'atomic mapping (from primitive: pc->ic)'
            call id_print_map(pc%ic_id)
            !
            if (present(uc)) then
            write(*,'(a5,a)',advance='no') ' ... ', 'atomic mapping (from input: uc->ic)'
            call id_print_map(uc%ic_id)
            endif
            !
        endif
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
        character(100) :: fname
        !
        call print_title('Atomic orbitals')
        !
        fname = 'infile.tb_orbitals'
        if (fexists(fname)) then
            ! read it
            call read_orbital_basis(ic,fname,opts%verbosity)
        else
            ! create basis on each irreducible atoms
            allocate(ic%atom(ic%natoms))
            do i = 1, ic%natoms
                call ic%atom(i)%gen_orbitals(orbital_flags='2s,2p')
            enddo
            call write_orbital_basis(ic)
        endif
        !
        if (opts%verbosity.ge.1) call print_orbital_basis(ic)
        !
        contains
        subroutine     write_orbital_basis(ic)
            !
            implicit none
            !
            class(am_class_irre_cell), intent(inout) :: ic
            integer :: i, j
            integer :: fid
            !
            fid = 1
            open(unit=fid,file='outfile.tb_orbitals',status='replace',action='write')
                ! spin polarized?
                write(fid,'(a)') 'spin: off'
                ! irreducible atoms
                write(fid,'(a,a)') 'irreducible atoms: ', tostring(ic%natoms)
                ! loop over irreducible atoms
                do i = 1, ic%natoms
                    write(fid,'(i10)',advance='no') i
                    ! loop over azimuthal quantum numbers
                    do j = 1, ic%atom(i)%nazimuthals
                        write(fid,'(a)',advance='no') ' '//trim(ic%atom(i)%orbname(j))
                    enddo
                    write(fid,*)
                enddo
            close(fid)
            !
        end subroutine write_orbital_basis
        subroutine     read_orbital_basis(ic,fname,verbosity)
            !
            implicit none
            !
            class(am_class_irre_cell), intent(inout) :: ic
            character(*), intent(in) :: fname
            integer, intent(in) :: verbosity
            integer :: i,j,k
            integer :: fid
            integer :: iostat
            character(maximum_buffer_size) :: buffer ! read buffer
            character(len=:), allocatable :: word(:) ! read buffer
            integer :: ic_natoms, nazimuthals
            character(len=100) :: orbital_flags
            character(len=100) :: spin_polarized_flag
            !
            fid = 1
            open(unit=fid,file=trim(fname),status="old",action='read')
                !
                write(*,'(a5,a,a)') ' ... ', 'atomic orbitals read = ', trim(fname)
                ! spin polarized?
                read(unit=fid,fmt='(a)') buffer
                word = strsplit(buffer,delimiter=' ')
                read(word(2),*) spin_polarized_flag
                ! irreducible atoms
                read(unit=fid,fmt='(a)') buffer
                word = strsplit(buffer,delimiter=' ')
                read(word(3),*) ic_natoms
                write(*,'(a5,a,a)') ' ... ', 'irreducible atoms = ', tostring(ic_natoms)
                ! set irreducible atoms
                if (ic%natoms.ne.ic_natoms) stop 'number of irreducible atoms input does not match internally calculated.'
                allocate(ic%atom(ic_natoms))
                ! loop over irreducible atoms
                do j = 1, ic_natoms
                    read(unit=fid,fmt='(a)') buffer
                    word = strsplit(buffer,delimiter=' ')
                    read(word(1),*) i
                    !
                    nazimuthals = size(word) - 1
                    write(*,'(a5,a,a,a)', advance='no') ' ... ', 'atom ', tostring(j), ' azimuthals ('//tostring(nazimuthals)//'):'
                    if (nazimuthals.le.0) stop 'number of azimuthals < 0'
                    !
                    orbital_flags=trim(spin_polarized_flag)
                    do k = 1, nazimuthals
                        orbital_flags=trim(orbital_flags)//','//trim(word(k+1))
                        write(*,'(a)',advance='no') ' '//trim(word(k+1))
                    enddo
                    write(*,*)
                    !
                    call ic%atom(i)%gen_orbitals(orbital_flags=orbital_flags)
                    !
                enddo
            close(fid)
            !
        end subroutine read_orbital_basis
        subroutine     print_orbital_basis(ic)
            !
            implicit none
            !
            class(am_class_irre_cell) , intent(inout) :: ic
            integer :: n,l,m,s
            integer :: i,j,a
            !
            do i = 1, ic%natoms
                write(*,'(a5,a)') ' ... ', 'irreducible atom '//tostring(i)//' contributes '//tostring((ic%atom(i)%norbitals))//' orbitals |n,l,m,s> :'
                do a = 1, ic%atom(i)%nazimuthals
                    write(*,'(5x,a,a)') trim(ic%atom(i)%orbname(a)), ':'
                do j = 1, ic%atom(i)%norbitals
                    n = ic%atom(i)%orbital(1,j)
                    l = ic%atom(i)%orbital(2,j)
                    m = ic%atom(i)%orbital(3,j)
                    s = ic%atom(i)%orbital(4,j)
                    if (ic%atom(i)%azimuthal(a).eq.l) write(*,'(5x,a2,i2,a1,i2,a1,i2,a1,i2,a2)') ' |', n, ',', l, ',', m, ',', s, '> '
                enddo
                enddo
            enddo
        end subroutine print_orbital_basis
    end subroutine initialize_orbitals

end module am_irre_cell