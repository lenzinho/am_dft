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
        contains
    end type am_class_tb_pg

    type, public, extends(am_class_tensor) :: am_class_tbvsk
    end type am_class_tbvsk

    type, public :: am_class_tight_binding
        integer :: nshells ! how many shells irreducible atoms
        type(am_class_tbvsk), allocatable :: tbvsk(:)  ! tbvsk(nshells)
    end type am_class_tight_binding

! 	type, public :: am_class_tb
! 		!
! 		! matrix elements
!         integer :: nshells ! how many shells irreducible atoms
!         !
! 		integer :: nVsks ! number of irreducible matrix elements
! 		real(dp), allocatable :: Vsk(:)   ! irreducible matrix elements
! 		integer , allocatable :: ind(:,:) ! ind([i,n,l,m,s,j,n,l,m,s,k],:)
!         ! i irreducible cell atom, jth neighbor shell around i, li orbital on i, lj orbital on j, m overlap
! 		!
! 	contains
!         procedure :: get_Vsk   ! gets irreducible matrix element given overlap indicies 
!         procedure :: set_Vsk   ! sets irreducible matrix element given overlap indicies
!         procedure :: build_Vsk ! 
!         procedure :: print_Vsk ! to stdout
!         procedure :: write_Vsk ! to file
!    		procedure :: read_Vsk  ! from file
! 	end type am_class_tb

contains

    subroutine     initialize_tb(tb,pg,ip,ic,pc,opts)
        !
        implicit none
        !
        class(am_class_tight_binding), intent(out):: tb
        type(am_class_point_group)   , intent(in) :: pg ! seitz point group
        type(am_class_prim_cell)     , intent(in) :: pc ! primitive cell
        type(am_class_irre_cell)     , intent(inout) :: ic ! irreducible cell
        class(am_class_pair_shell)   , intent(in) :: ip ! irreducible pairs
        type(am_class_options)       , intent(in) :: opts
        type(am_class_options)                    :: notalk
        type(am_class_point_group)                :: stab
        integer :: i,j,k
        !
        notalk = opts
        notalk%verbosity = 0
        !
        ! ADD OPTION TO READ ORBITALS HERE
        call ic%initialize_orbitals(opts=opts)
        ! call ic%read_orbitals(opts=opts)
        !
        !
        tb%nshells = ip%nshells
        !
        allocate(tb%tbvsk(tb%nshells))
        do k = 1, ip%nshells
            ! get rank, dims, flags, property
            call initialize_tbvsk(tbvsk=tb%tbvsk(i), pc=pc, ic=ic)
            ! determine stabilizers of a prototypical bond in shell (vector v)
            call stab%get_stabilizer_group(pg=pg, v=ip%shell(k)%tau(1:3,1), opts=notalk)
            !
            call tb%tbvsk(i)%fpg%get_flat_point_group(tens=tb%tbvsk(i), pg=pg, pc=pc)
            ! call revg%get_reversal_group(sg=sg, v=ip%shell(k)%tau(1:3,1), opts=notalk)
        enddo
        !
        contains
        subroutine     initialize_tbvsk(tbvsk,pc,ic)
            !
            implicit none
            !
            class(am_class_tbvsk), intent(out):: tbvsk
            class(am_class_prim_cell), intent(in) :: pc
            class(am_class_irre_cell), intent(in) :: ic
            !
            tbvsk%property = 'tb_matrix_elements'
            tbvsk%rank  = 2
            tbvsk%flags = 'tight' 
            tbvsk%dims  = shape(ps2tb(R=eye(3),pc=pc,ic=ic))
            !
        end subroutine initialize_tbvsk
    end subroutine initialize_tb

































! subroutine     get_tb_relations(tb,ip,ic,opts)
!     !
!     implicit none
!     !
!     class(am_class_tb)        , intent(inout) :: tb
!     class(am_class_pair_shell), intent(inout) :: ip
!     class(am_class_irre_cell) , intent(inout) :: ic
!     type(am_class_options)    , intent(in) :: opts
!     integer, allocatable :: ind(:,:)
!     integer :: nVsks
!     integer :: label_dim
!     integer :: i,j,alpha,beta,k,kk
!     !
!     ! get total matrix elements
!     nVsks = 0
!     do k = 1, ip%nshells
!         nVsks = nVsks + ic%atom( ip%shell(k)%i )%norbitals * ic%atom( ip%shell(k)%j )%norbitals
!     enddo
!     call am_print('total matrix elements',nVsks)
!     !
!     ! determine irreducible matrix elements
!     label_dim = 11
!     allocate(ind(label_dim,nVsks))
!     kk=0
!     do k = 1, ip%nshells
!         ! irreducible atoms i and j
!         i = ip%shell(k)%i
!         j = ip%shell(k)%j
!         ! loop over all combination of orbitals (outer product of li x lj)
!         do alpha = 1, ic%atom(i)%norbitals
!         do beta  = 1, ic%atom(j)%norbitals
!             ! matrix element vanishes unless m = m'
!             if (ic%atom(i)%orbital(3,alpha).eq.ic%atom(j)%orbital(3,beta)) then
!                 kk=kk+1
!                 ! (l,l',m) = (l,l',-m) and (l,l',m) = (-1)^(l-m) * (l',l,m) symmetries enforced in get_Vsk_ind
!                 ind(:,kk) = ic%get_Vsk_ind(i=i,alpha=alpha,j=j,beta=beta,k=k)
!                 !
!             endif
!         enddo
!         enddo
!     enddo
!     ! get number of matrix elements equal to zero (those for which m/=m')
!     call am_print('null matrix elements',nVsks-kk)
!     ! determine irrducible matrix elements
!     tb%ind = unique(ind(:,1:kk))
!     ! allocate space and initialize irrducible matrix elements as 1.0
!     tb%nVsks = size(tb%ind,2)
!     allocate(tb%Vsk(tb%nVsks))
!     tb%Vsk = 1.0_dp
!     ! write matrix elements to stdout
!     if (opts%verbosity.ge.1) call tb%print_Vsk()
!     !
! end subroutine initialize_orbitals

! ! functions which operate on Vsk

! function       get_Vsk(tb,ind) result(Vsk)
!     ! returns the corresponding matrix element given the indices
!     ! ind = [i,ni,li,mi,si,j,nj,lj,mj,sj,k]
!     implicit none
!     !
!     class(am_class_tb), intent(in) :: tb
!     integer, intent(in) :: ind(11)
!     real(dp) :: Vsk
!     integer  :: i,j
!     !
!     search : do j = 1, tb%nVsks
!         !
!         do i = 1, 11
!             if (tb%ind(i,j).ne.ind(i)) then
!                 cycle search
!             endif
!         enddo
!         !
!         Vsk = tb%Vsk(j)
!         !
!         return
!         !
!     enddo search
!     !
!     ! if nothing matches return zero
!     Vsk = 0.0_dp
!     !
! end function   get_Vsk

! subroutine     set_Vsk(tb,ind,Vsk)
!     ! set the corresponding matrix element given the indices
!     ! ind = [i,ni,li,mi,si,j,nj,lj,mj,sj,k]
!     implicit none
!     !
!     class(am_class_tb), intent(inout) :: tb
!     integer, intent(in) :: ind(1:11)
!     real(dp) :: Vsk
!     integer :: i,j
!     !
!     search : do j = 1, tb%nVsks
!         !
!         do i = 1, 11
!             if (tb%ind(i,j).ne.ind(i)) then
!                 cycle search
!             endif
!         enddo
!         !
!         tb%Vsk(j) = Vsk
!         !
!         return
!         !
!     enddo search
!     !
!     ! if nothing matches return error
!     call am_print('ERROR','Cannot find where to set Vsk',flags='E')
!     stop
!     !
! end subroutine set_Vsk


! subroutine     print_Vsk(tb)
!     !
!     use am_atom, only : l2spdf, m2spdf
!     !
!     implicit none
!     !
!     class(am_class_tb), intent(in) :: tb
!     integer :: k
!     !
!     !
!     call am_print('irreducible matrix elements',tb%nVsks)
!     !
!     write(*,'(5x,12a6,a20)') '#','i','n_i','l_i','m_i','s_i','j','n_j','l_j','m_j','s_j','shell','Vsk'
!     write(*,'(5x,a72,a20)') repeat(' -----',12), ' '//repeat('-',19)
!     !
!     do k = 1, tb%nVsks
!         write(*,'(5x, 3i6, 2a6, 3i6, 2a6, 2i6, f20.8)') &
!                     k, &
!                     tb%ind( 1,k),  & ! atom 1     
!                     tb%ind( 2,k),  & ! n   1
!         trim(l2spdf(tb%ind( 3,k))),& ! l   1 (s,p,d,f)
!         trim(m2spdf(tb%ind( 4,k))),& ! m   1 (sigma, pi, delta, phi)
!                     tb%ind( 5,k),  & ! s   1
!                     tb%ind( 6,k),  & ! atom 2
!                     tb%ind( 7,k),  & ! n   2
!         trim(l2spdf(tb%ind( 8,k))),& ! l   2 (s,p,d,f)
!         trim(m2spdf(tb%ind( 9,k))),& ! m   2 (sigma, pi, delta, phi)
!                     tb%ind(10,k),  & ! s   2
!                     tb%ind(11,k),  & ! shell
!                     tb%Vsk(k)        ! matrix element
!     enddo
!     !
!     write(*,'(5x,a)') 'Definitions:'
!     write(*,'(5x,a)') 'i , j : irreducible atoms indicies'
!     write(*,'(5x,a)') 'li,lj : orbitals on irreducible atoms'
!     write(*,'(5x,a)') 'm     : type of overlap'
!     !
! end subroutine print_Vsk

! subroutine     build_Vsk(tb,ic,ip)
!     !
!     ! Get k-independent part of the Hamiltonian, ip%shell(k)%Vsk
!     !
!     implicit none
!     !
!     class(am_class_tb),         intent(in) :: tb ! tight binding parametrs
!     class(am_class_irre_cell) , intent(in) :: ic ! irreducible cell
!     class(am_class_pair_shell), intent(inout) :: ip ! irreducible pairs
!     integer :: k ! shell index
!     integer :: i, j ! i and j irreducible atoms 
!     integer :: Ki,Kj ! number of orbitals on i and j
!     integer :: alpha, beta ! orbital indices on i and j
!     !
!     ! loop over irreducible pairs
!     do k = 1, ip%nshells
!         ! irreducible atom indices
!         i = ip%shell( k )%i
!         j = ip%shell( k )%j
!         ! get number of orbitals (size of matrix)
!         Ki = ic%atom( i )%norbitals
!         Kj = ic%atom( j )%norbitals
!         ! allocate space
!         allocate( ip%shell( k )%Vsk(Ki,Kj) )
!         ! loop over orbitals
!         do alpha = 1, Ki
!         do beta  = 1, Kj
!             !
!             ip%shell( k )%Vsk( alpha, beta ) = tb%get_Vsk(ind=ic%get_Vsk_ind(i=i,alpha=alpha,j=j,beta=beta,k=k))
!             !
!         enddo
!         enddo

!         call am_print('K'//trim(int2char(k)), ip%shell( k )%Vsk )
!     enddo
! end subroutine build_Vsk

! subroutine     write_Vsk(tb)
! 	!
!     implicit none
!     !
!     class(am_class_tb), intent(inout) :: tb
!     integer :: fid, k
!     !
!     fid = 1
!     open(unit=fid,file='outfile.tightbinding',status='replace',action='write')
!     !
!     call am_print('irreducible matrix elements',tb%nVsks)
!         !
!         write(fid,'(i6)') tb%nVsks
!         !
!         write(fid,'(5x,11a3,a20)') 'i','n','l','m','s','j','n','l','m','s','sh','Vsk'
!         !
!         do k = 1, tb%nVsks
!             write(fid,'(5x, 11i3, f20.8)') &
!                         tb%ind( 1,k), & ! atom 1     
!                         tb%ind( 2,k), & ! n   1
!                         tb%ind( 3,k), & ! l   1 (s,p,d,f)
!                         tb%ind( 4,k), & ! m   1 (sigma, pi, delta, phi)
!                         tb%ind( 5,k), & ! s   1
!                         tb%ind( 6,k), & ! atom 2
!                         tb%ind( 7,k), & ! n   2
!                         tb%ind( 8,k), & ! l   2 (s,p,d,f)
!                         tb%ind( 9,k), & ! m   2 (sigma, pi, delta, phi)
!                         tb%ind(10,k), & ! s   2
!                         tb%ind(11,k), & ! shell
!                         tb%Vsk(k)       ! matrix element
!         enddo
!         !
!     close(fid)
!     !
! end subroutine write_Vsk

! subroutine     read_Vsk(tb,opts)
! 	!
! 	implicit none
! 	!
! 	class(am_class_tb), intent(inout) :: tb
!     type(am_class_options), intent(in) :: opts
!     character(maximum_buffer_size) :: buffer ! read buffer
!     character(len=:), allocatable :: word(:) ! read buffer
!     integer :: fid, k
! 	!
! 	if (opts%verbosity.ge.1) call am_print_title('Reading tight binding matrix elements')
! 	!
! 	if (opts%verbosity.ge.1) call am_print('input file',trim(opts%tbf))
! 	!
!     fid = 1
!     open(unit=fid,file=trim(opts%tbf),status="old",action='read')
!         !
!         read(unit=fid,fmt='(a)') buffer
!         word = strsplit(buffer,delimiter=' ')
!         read(word(1),*) tb%nVsks
!         !
!         allocate(tb%ind(11,tb%nVsks))
!         allocate(tb%Vsk(tb%nVsks))
!         !
!         do k = 1, tb%nVsks
!             read(unit=fid,fmt='(a)') buffer
!             word = strsplit(buffer,delimiter=' ')
!             !
!             read(word( 1),*) tb%ind( 1,k) ! atom 1     
!             read(word( 2),*) tb%ind( 2,k) ! n   1
!             read(word( 3),*) tb%ind( 3,k) ! l   1 (s,p,d,f)
!             read(word( 4),*) tb%ind( 4,k) ! m   1 (sigma, pi, delta, phi)
!             read(word( 5),*) tb%ind( 5,k) ! s   1
!             read(word( 6),*) tb%ind( 6,k) ! atom 2
!             read(word( 7),*) tb%ind( 7,k) ! n   2
!             read(word( 8),*) tb%ind( 8,k) ! l   2 (s,p,d,f)
!             read(word( 9),*) tb%ind( 9,k) ! m   2 (sigma, pi, delta, phi)
!             read(word(10),*) tb%ind(10,k) ! s   2
!             read(word(11),*) tb%ind(11,k) ! shell
!             read(word(12),*) tb%Vsk(k)    ! matrix element
!             !
!         enddo
!         !
!         call tb%print_Vsk
!         !
!     close(fid)
! 	!
! end subroutine read_Vsk

! ! functions used to setup hamiltonian

! function       get_Hamiltonian(ip,ic,pp,pc,kpt) result(H)
!     ! 
!     ! Get tight binding Hamiltonian at kpt.
!     ! 
!     implicit none
!     !
!     class(am_class_pair_shell), intent(in) :: ip ! irreducible pairs
!     type(am_class_irre_cell)  , intent(in) :: ic ! irreducible cell
!     type(am_class_pair_shell) , intent(in) :: pp ! primitive pairs
!     type(am_class_prim_cell)  , intent(in) :: pc ! primitive cell
!     real(dp)                  , intent(in) :: kpt(3) ! fractional
!     complex(dp), allocatable, target :: H(:,:)
!     integer, allocatable :: H_start(:), H_end(:)
!     complex(dp), pointer :: H_sub(:,:) 
!     integer :: Hdim ! hamiltonian dimensions
!     real(dp), allocatable :: T(:,:) ! similarity transform usd to get molecular axis paralle to z
!     real(dp) :: R(3) ! vector connecting atom to shell
!     real(dp) :: E ! exponential factor in bloch sum
!     integer  :: k ! shell index
!     integer  :: m ! primitive atom 1 index 
!     integer  :: n ! primitive atom 2 index
!     integer  :: i,j ! loop variable
!     !
!     ! allocate space for Hamiltonian and determine the subsections of the Hamiltonian corresponding to each primitive atom
!     Hdim = 0
!     allocate(H_start(pc%natoms))
!     allocate(H_end(pc%natoms))
!     do i = 1, pc%natoms
!         H_start(i) = Hdim + 1
!         Hdim = Hdim + ic%atom( pc%ic_id(i) )%norbitals
!         H_end(i) = Hdim
!     enddo
!     !
!     ! allocate and initialize
!     allocate(H(Hdim,Hdim))
!     H = 0.0_dp
!     !
!     ! construct Hamiltonian
!     do k = 1, pp%nshells
!         ! primitive atom indicies
!         m = pp%shell( k )%m
!         n = pp%shell( k )%n
!         ! set pointer to subsection of hamiltonian
!         H_sub => H(H_start(m):H_end(m), H_start(n):H_end(n))
!         ! loop over atoms in shell
!         do j = 1, pp%shell( k )%natoms
!             ! get vector connecting primitive atom to atom j in shell k
!             R = pp%shell(k)%tau(1:3,j)
!             ! get rotation that will transform matrix elements


!             ! WORKING ON STUFF HERE.... 
!             ! T = O3_rotations(azimuthal=[ic%atom(m)%azimuthal,ic%atom(n)%azimuthal],th=th,phi=phi,is_spin_polarized=is_spin_polarized)



!             ! get exponent factor
!             E = exp(-cmplx_i*dot_product(R,kpt))
!             ! multiply everything and add to the sum
!             H_sub = H_sub + matmul(transpose(T), matmul(ip%shell(pp%ip_id(k))%Vsk, T)) * E
!             !
!         enddo
!         !
!         call am_print('H'//trim(int2char(k)),H)
!         !
!     enddo
! end function   get_Hamiltonian

! function       O3_rotations(azimuthal,th,phi,is_spin_polarized) result(R)
!     !
!     implicit none
!     !
!     integer , intent(in) :: azimuthal(:)
!     real(dp), intent(in) :: th,phi
!     logical , intent(in) :: is_spin_polarized
!     real(dp), allocatable :: R(:,:)
!     integer , allocatable :: l_start(:)
!     integer , allocatable :: l_end(:)
!     integer :: i, n, m
!     !
!     m = size(azimuthal)
!     !
!     allocate(l_start(m))
!     allocate(l_end(m))
!     !
!     n = 0
!     do i = 1,m
!         l_start(i) = n + 1
!         n = n + 2*azimuthal(i)+1
!         l_end(i) = n
!     enddo
!     !
!     allocate(R(n,n))
!     R = 0.0_dp
!     !
!     ! construct a direct sum of rotations
!     do i = 1, m
!         R(l_start(i):l_end(i),l_start(i):l_end(i)) = euler2SO3(l=azimuthal(i),euler=[0.0_dp,th,phi])
!     enddo
!     !
!     if (is_spin_polarized) then
!         ! need to add something to account for the different spins. perhaps a kroner product of the R matrix obtained above and the corresponding rotations using the pauling spin matrices?
!         call am_print('ERROR','spin polarized is not yet supported.')
!         stop
!     endif
!     !
! end function   O3_rotations











end module am_tight_binding











