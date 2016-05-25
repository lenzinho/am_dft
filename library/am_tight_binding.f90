module am_tight_binding

    use am_constants
    use am_stdout
    use am_matlab
    use am_options
    use am_shells
    use am_irre_cell
    use am_prim_cell
    use am_symmetry
    use am_symmetry_rep
    use am_mkl

	implicit none

	private

    type, public, extends(am_class_group) :: am_class_tb_pg
    end type am_class_tb_pg

    type, public, extends(am_class_tensor) :: am_class_tens_tb
        ! <INHERITED: am_class_tensor>
        ! character(100)        :: property       ! name of propertty
        ! character(100)        :: flags          ! axial/polar
        ! integer               :: rank           ! tensor rank
        ! integer , allocatable :: dims(:)        ! tensor dimensions
        ! real(dp), allocatable :: relations(:,:) ! relations connecting tensor elements
        ! </INHERITED: am_class_tensor>
    end type am_class_tens_tb

    type, public :: am_class_tight_binding
        integer :: nshells ! how many shells irreducible atoms
        type(am_class_tens_tb), allocatable :: tbvsk(:)  ! tbvsk(nshells)
        !
        integer, allocatable :: niVs            ! number of irreducible matrix element values
        integer, allocatable :: iV(:)           ! their values
        integer, allocatable :: iV_indices(:,:) ! their indices iV_indices( [i,j], niVs); j = compound [alpha x beta index, see set Vsk subroutine]
        !
        type(am_class_hamiltonian_group) :: ham_pg
        !
        contains
            procedure :: initialize_tb
    end type am_class_tight_binding

contains

    subroutine     initialize_tb(tb,sg,pg,ip,ic,pc,opts)
        !
        implicit none
        !
        class(am_class_tight_binding), intent(out) :: tb   ! tight binding parameters
        type(am_class_space_group)   , intent(in)  :: sg   ! seitz space group
        type(am_class_point_group)   , intent(in)  :: pg   ! seitz point group
        type(am_class_prim_cell)     , intent(in)  :: pc   ! primitive cell
        type(am_class_irre_cell)     , intent(in)  :: ic ! irreducible cell
        type(am_class_pair_shell)    , intent(in)  :: ip   ! irreducible pairs
        type(am_class_options)       , intent(in)  :: opts
        logical, allocatable :: is_independent(:)
        integer :: nterms
        integer :: i,j,k
        integer :: a,b,c
        integer :: m,n,o
        !
        !
        if (opts%verbosity.ge.1) call am_print_title('Imposing symmetry constraints on tight-binding model')
        !
        call tb%ham_pg%get_hamiltonian_point_group(pg=pg, pc=pc, ic=ic)
        !
        tb%nshells = ip%nshells
        allocate(tb%tbvsk(tb%nshells))
        !
        do k = 1, tb%nshells
            call get_shell_relations(tbvsk=tb%tbvsk(k), sg=sg, pg=pg, pc=pc, ic=ic, shell=ip%shell(k), opts=opts)
        enddo
        ! get number of independent (irreducible) matrix elements
        tb%niVs = 0
        do i = 1, tb%nshells
            tb%niVs = tb%niVs + count(get_independent(tb%tbvsk(i)%relations)) 
        enddo
        ! allocate space for iV independent (irreducible) matrix elements
        allocate(tb%iV(tb%niVs))
        tb%iV = 0
        ! indices : iV_indices( [i, subd2ind(dims=tb%tbvsk(i)%dims,sub=[alpha,beta])], niVs)
        allocate(tb%iV_indices(2,tb%niVs))
        tb%iV_indices = 0
        k=0
        do i = 1, tb%nshells
            is_independent = get_independent(tb%tbvsk(i)%relations)
            do j = 1, product(tb%tbvsk(i)%dims)
                k=k+1
                if (is_independent(j)) then
                    tb%iV_indices(1,k) = i
                    tb%iV_indices(2,k) = j
                endif
            enddo
        enddo
        !
        if (opts%verbosity.ge.1) then
            ! print statistics about parameters
            write(*,'(a5,2a6,3a14)') ' ... ', 'shell', 'terms', 'null', 'dependent', 'independent'
            write(*,'(5x,a)') repeat(' '//repeat('-',5),2)//repeat(' '//repeat('-',13),3)
            ! print table
            a=0;b=0;c=0;nterms=0
            do i = 1, tb%nshells
                m = count(get_null(tb%tbvsk(i)%relations))
                n = count(get_depenent(tb%tbvsk(i)%relations)) 
                o = count(get_independent(tb%tbvsk(i)%relations))
                nterms = m+n+o
                write(*,'(5x,i6)'        ,advance='no') i
                write(*,'(i6)'           ,advance='no') nterms
                write(*,'(i5,a2,f5.1,a2)',advance='no') m, '(', (m*100_dp)/real(nterms,dp) , '%)'
                write(*,'(i5,a2,f5.1,a2)',advance='no') n, '(', (n*100_dp)/real(nterms,dp) , '%)'
                write(*,'(i5,a2,f5.1,a2)',advance='no') o, '(', (o*100_dp)/real(nterms,dp) , '%)'
                write(*,*)
                a = a + m
                b = b + n
                c = c + o
            enddo
            ! print total
            write(*,'(5x,a6)'        ,advance='no') 'total'
            write(*,'(i6)'           ,advance='no') (a+b+c)
            write(*,'(i5,a2,f5.1,a2)',advance='no') a, '(', (a*100_dp)/real(a+b+c,dp) , '%)'
            write(*,'(i5,a2,f5.1,a2)',advance='no') b, '(', (b*100_dp)/real(a+b+c,dp) , '%)'
            write(*,'(i5,a2,f5.1,a2)',advance='no') c, '(', (c*100_dp)/real(a+b+c,dp) , '%)'
            write(*,*)
            ! print symmetry relations
            do i = 1, tb%nshells
                write(*,'(a5,a)') ' ... ', 'shell '//trim(int2char(i))//' irreducible symmetry relations:'
                call print_relations(relations=tb%tbvsk(i)%relations, dims=tb%tbvsk(i)%dims, flags='print:dependent,independent')
            enddo
        endif
        !
        contains
        subroutine     initialize_tbvsk(tbvsk,pc,ic,shell)
            !
            implicit none
            !
            type(am_class_tens_tb)  , intent(out) :: tbvsk
            type(am_class_prim_cell), intent(in)  :: pc
            type(am_class_irre_cell), intent(in)  :: ic
            type(am_shell_cell)     , intent(in)  :: shell
            integer :: i
            !
            tbvsk%property = 'tight binding'
            tbvsk%rank  = 2
            ! detetmine if irreducible atoms are the same
            if (shell%i.eq.shell%j) then
                tbvsk%flags = 'symmetric'
            else
                tbvsk%flags = 'asymmetric'
            endif
            ! set dimensions of matrix elements
            allocate(tbvsk%dims(tbvsk%rank))
            tbvsk%dims = 0 
            ! dimension of Hamiltonian subsection corresponding to primitive atom m (irreducible atoms i)
            do i = 1, ic%atom(shell%i)%nazimuthals
                tbvsk%dims(1) = tbvsk%dims(1) + ic%atom(shell%i)%azimuthal(i)*2+1
            enddo
            ! dimension of Hamiltonian subsection corresponding to primitive atom n (irreducible atoms j)
            do i = 1, ic%atom(shell%j)%nazimuthals
                tbvsk%dims(2) = tbvsk%dims(2) + ic%atom(shell%j)%azimuthal(i)*2+1
            enddo
            ! allocate space for matrix elements
            allocate(tbvsk%V(tbvsk%dims(1)*tbvsk%dims(2)))
            ! allocate space for matrix descriptors
        end subroutine initialize_tbvsk
        subroutine     get_shell_relations(tbvsk,sg,pg,pc,ic,shell,opts)
                !
                implicit none
                !
                type(am_class_tens_tb)      , intent(inout) :: tbvsk
                type(am_class_space_group), intent(in) :: sg   ! seitz space group
                type(am_class_point_group), intent(in) :: pg   ! seitz point group
                type(am_class_prim_cell)  , intent(in) :: pc
                type(am_class_irre_cell)  , intent(in) :: ic
                type(am_shell_cell)       , intent(in) :: shell
                type(am_class_options)    , intent(in) :: opts
                type(am_class_flat_group)  :: flat_pg
                type(am_class_flat_group)  :: flat_ig
                type(am_class_point_group) :: stab
                !
                ! get rank, dims, flags, property
                call initialize_tbvsk(tbvsk=tbvsk, pc=pc, ic=ic, shell=shell)
                !
                ! determine intrinsic symmetries (if both irreducible atoms are of the same irreducible type)
                ! Interchange of indices: (l,l',m) = (-1)^(l+l') * (l',l,m), due to parity of wavefunction. (s,d are even under inversion, p,f are odd)
                ! E. Scheer, Molecular Electronics: An Introduction to Theory and Experiment, p 245.
                call flat_ig%get_flat_intrinsic_group(tens=tbvsk, atom_m=ic%atom(shell%i), atom_n=ic%atom(shell%j) )
                ! 
                ! determine stabilizers relations
                call stab%get_stabilizer_group(pg=pg, v=shell%tau(1:3,1), opts=opts)
                ! get stabilizer symmetries in the flattened hamiltonin basis
                call flat_pg%get_flat_point_group(tens=tbvsk, pg=stab, pc=pc, atom_m=ic%atom(shell%i), atom_n=ic%atom(shell%j))
                !
                ! get combined relations
                tbvsk%relations = combine_relations(relationsA=flat_pg%relations, relationsB=flat_ig%relations)
                !
        end subroutine get_shell_relations
    end subroutine initialize_tb

    function       get_Hamiltonian(tb,ip,ic,pp,pc,kpt) result(H)
        ! 
        ! Get tight binding Hamiltonian at kpt.
        ! 
        implicit none
        !
        class(am_class_tight_binding), intent(out):: tb ! tight binding matrix elements
        class(am_class_pair_shell)   , intent(in) :: ip ! irreducible pairs
        type(am_class_irre_cell)     , intent(in) :: ic ! irreducible cell
        type(am_class_pair_shell)    , intent(in) :: pp ! primitive pairs
        type(am_class_prim_cell)     , intent(in) :: pc ! primitive cell
        real(dp)                     , intent(in) :: kpt(3) ! fractional
        complex(dp), allocatable :: H(:,:)
        integer    , allocatable :: H_start(:), H_end(:)
        complex(dp) :: E ! exponential factor in bloch sum
        integer :: Hdim ! hamiltonian dimensions
        integer :: m ! primitive atom 1 index 
        integer :: n ! primitive atom 2 index
        integer :: i,j,k
        !
        ! allocate space for Hamiltonian
        ! determine Hamiltonian subsections corresponding to primitive atoms
        Hdim = 0
        allocate(H_start(pc%natoms))
        allocate(H_end(pc%natoms))
        do i = 1, pc%natoms
            H_start(i) = Hdim + 1
            Hdim = Hdim + ic%atom( pc%ic_id(i) )%norbitals
            H_end(i) = Hdim
        enddo
        !
        ! allocate and initialize
        allocate(H(Hdim,Hdim))
        H = cmplx(0,0,dp)
        !
        ! construct Hamiltonian
        do k = 1, pp%nshells
            ! primitive atom indicies
            m = pp%shell(k)%m
            n = pp%shell(k)%n
            ! irreducible atom indicies
            i = pp%shell(k)%i
            j = pp%shell(k)%j
            ! compute bloch sum by loop over atoms in shell
            E = cmplx(0,0,dp)
            do j = 1, pp%shell(k)%natoms
                ! get vector connecting primitive atom to atom j in shell k
                R = pp%shell(k)%tau(1:3,j)
                !
                ! NEED TO ROTATE V_sk
                H(H_start(m):H_end(m), H_start(n):H_end(n)) = get_Vsk(tb=tb, i=pp%ip_id(k))
                !
                if (pp%shell(k)%ip_id.lt.0) then
                ! if ip_id < 0 -> the irreducible pair has its indices flipped with respect to the primitive pair
                ! as a result, the matrix elements need to be transposed and the signed adjusted to reflect oribtal parity
                H(H_start(m):H_end(m), H_start(n):H_end(n)) = &
                H(H_start(m):H_end(m), H_start(n):H_end(n)) * transp_parity_sign(atom_m=ic%atom(i), atom_n=ic%atom(j))
                endif
                ! multiply exponential factor from Bloch sum
                H(H_start(m):H_end(m), H_start(n):H_end(n)) = &
                H(H_start(m):H_end(m), H_start(n):H_end(n)) * exp(-cmplx_i*dot_product(R,kpt))
                !
            enddo
            !
            call am_print('H'//trim(int2char(k)),H)
            !
        enddo
    end function   get_Hamiltonian

    ! functions which operate on V

    function       get_Vsk(tb,i) result(V)
        !
        ! irreducible pair (i), can be negative (corresponds to pair n-m rathet than m-n, on the
        ! opposite [upper/lower] side of the Hamiltonian), in which case adjoint of V is returned
        !
        implicit none
        !
        class(am_class_tight_binding), intent(in) :: tb
        integer , intent(in) :: i
        real(dp), allocatable :: V(:,:)
        integer :: m,n
        !
        m = tb%tbvsk(i)%dims(1)
        n = tb%tbvsk(i)%dims(2)
        !
        if (i.gt.0) then
            allocate(V(m,n))
            V =          reshape(tb%tbvsk(i)%V, [m,n])
        else
            allocate(V(tb%tbvsk(i)%dims(2), tb%tbvsk(i)%dims(1)))
            V = adjoint( reshape(tb%tbvsk(i)%V, [m,n]) )
        endif
        !
    end function   get_Vsk

    subroutine     set_Vsk(tb,i,alpha,beta,V)
        !
        implicit none
        !
        class(am_class_tight_binding), intent(inout) :: tb
        integer, intent(in) :: i           ! irreducible pair (i)
        integer, intent(in) :: alpha, beta ! orbitals on atoms alpha,beta (matrix coordinates)
        real(dp),intent(in) :: V           ! matrix element
        integer :: ind
        !
        ind = sub2ind(dims=tb%tbvsk(i)%dims,sub=[alpha,beta])
        !
        tb%tbvsk(i)%V(ind) = V
        !
    end subroutine set_Vsk


!     subroutine     apply_relations(Vsk,is_independent,relations)
!         !
!         implicit none
!         !
!         real(dp), intent(in) :: Vsk(:)
!         real(dp), intent(in) :: relations(:,:)
!         !
!         !
!         tb%tbvsk(i)%Vsk(ind) = V
!         !
!         matmul(relations,V)
!         !
!     end subroutine apply_relations

!     subroutine     print_Vsk(tb)
!         !
!         use am_atom, only : l2spdf, m2spdf
!         !
!         implicit none
!         !
!         class(am_class_tb), intent(in) :: tb
!         integer :: k
!         !
!         !
!         call am_print('irreducible matrix elements',tb%nVsks)
!         !
!         write(*,'(5x,12a6,a20)') '#','i','n_i','l_i','m_i','s_i','j','n_j','l_j','m_j','s_j','shell','Vsk'
!         write(*,'(5x,a72,a20)') repeat(' -----',12), ' '//repeat('-',19)
!         !
!         do k = 1, tb%nVsks
!             write(*,'(5x, 3i6, 2a6, 3i6, 2a6, 2i6, f20.8)') &
!                         k, &
!                         tb%ind( 1,k),  & ! atom 1     
!                         tb%ind( 2,k),  & ! n   1
!             trim(l2spdf(tb%ind( 3,k))),& ! l   1 (s,p,d,f)
!             trim(m2spdf(tb%ind( 4,k))),& ! m   1 (sigma, pi, delta, phi)
!                         tb%ind( 5,k),  & ! s   1
!                         tb%ind( 6,k),  & ! atom 2
!                         tb%ind( 7,k),  & ! n   2
!             trim(l2spdf(tb%ind( 8,k))),& ! l   2 (s,p,d,f)
!             trim(m2spdf(tb%ind( 9,k))),& ! m   2 (sigma, pi, delta, phi)
!                         tb%ind(10,k),  & ! s   2
!                         tb%ind(11,k),  & ! shell
!                         tb%Vsk(k)        ! matrix element
!         enddo
!         !
!         write(*,'(5x,a)') 'Definitions:'
!         write(*,'(5x,a)') 'i , j : irreducible atoms indicies'
!         write(*,'(5x,a)') 'li,lj : orbitals on irreducible atoms'
!         write(*,'(5x,a)') 'm     : type of overlap'
!         !
!     end subroutine print_Vsk

!     subroutine     write_Vsk(tb)
!       !
!         implicit none
!         !
!         class(am_class_tb), intent(inout) :: tb
!         integer :: fid, k
!         !
!         fid = 1
!         open(unit=fid,file='outfile.tightbinding',status='replace',action='write')
!         !
!         call am_print('irreducible matrix elements',tb%nVsks)
!             !
!             write(fid,'(i6)') tb%nVsks
!             !
!             write(fid,'(5x,11a3,a20)') 'i','n','l','m','s','j','n','l','m','s','sh','Vsk'
!             !
!             do k = 1, tb%nVsks
!                 write(fid,'(5x, 11i3, f20.8)') &
!                             tb%ind( 1,k), & ! atom 1     
!                             tb%ind( 2,k), & ! n   1
!                             tb%ind( 3,k), & ! l   1 (s,p,d,f)
!                             tb%ind( 4,k), & ! m   1 (sigma, pi, delta, phi)
!                             tb%ind( 5,k), & ! s   1
!                             tb%ind( 6,k), & ! atom 2
!                             tb%ind( 7,k), & ! n   2
!                             tb%ind( 8,k), & ! l   2 (s,p,d,f)
!                             tb%ind( 9,k), & ! m   2 (sigma, pi, delta, phi)
!                             tb%ind(10,k), & ! s   2
!                             tb%ind(11,k), & ! shell
!                             tb%Vsk(k)       ! matrix element
!             enddo
!             !
!         close(fid)
!         !
!     end subroutine write_Vsk

!     subroutine     read_Vsk(tb,opts)
!       !
!       implicit none
!       !
!       class(am_class_tb), intent(inout) :: tb
!         type(am_class_options), intent(in) :: opts
!         character(maximum_buffer_size) :: buffer ! read buffer
!         character(len=:), allocatable :: word(:) ! read buffer
!         integer :: fid, k
!       !
!       if (opts%verbosity.ge.1) call am_print_title('Reading tight binding matrix elements')
!       !
!       if (opts%verbosity.ge.1) call am_print('input file',trim(opts%tbf))
!       !
!         fid = 1
!         open(unit=fid,file=trim(opts%tbf),status="old",action='read')
!             !
!             read(unit=fid,fmt='(a)') buffer
!             word = strsplit(buffer,delimiter=' ')
!             read(word(1),*) tb%nVsks
!             !
!             allocate(tb%ind(11,tb%nVsks))
!             allocate(tb%Vsk(tb%nVsks))
!             !
!             do k = 1, tb%nVsks
!                 read(unit=fid,fmt='(a)') buffer
!                 word = strsplit(buffer,delimiter=' ')
!                 !
!                 read(word( 1),*) tb%ind( 1,k) ! atom 1     
!                 read(word( 2),*) tb%ind( 2,k) ! n   1
!                 read(word( 3),*) tb%ind( 3,k) ! l   1 (s,p,d,f)
!                 read(word( 4),*) tb%ind( 4,k) ! m   1 (sigma, pi, delta, phi)
!                 read(word( 5),*) tb%ind( 5,k) ! s   1
!                 read(word( 6),*) tb%ind( 6,k) ! atom 2
!                 read(word( 7),*) tb%ind( 7,k) ! n   2
!                 read(word( 8),*) tb%ind( 8,k) ! l   2 (s,p,d,f)
!                 read(word( 9),*) tb%ind( 9,k) ! m   2 (sigma, pi, delta, phi)
!                 read(word(10),*) tb%ind(10,k) ! s   2
!                 read(word(11),*) tb%ind(11,k) ! shell
!                 read(word(12),*) tb%Vsk(k)    ! matrix element
!                 !
!             enddo
!             !
!             call tb%print_Vsk
!             !
!         close(fid)
!       !
!     end subroutine read_Vsk











































end module am_tight_binding











