module am_shells

    use am_constants
    use am_helpers
    use am_unit_cell
    use am_prim_cell
    use am_irre_cell
    use am_options
    use am_mkl
    use am_atom
    use am_symmetry

    implicit none

    private

    public :: am_class_pair_shell

    type am_class_pair_shell
        !
        integer, allocatable :: nshells(:) ! how many pairs atom i has
        type(am_class_unit_cell), allocatable :: shell(:,:) ! shell(ic%natoms,npairs(natoms)) - pair%shell(i,j) pair between irreducible atom i and its shell j
        !
    contains
        !
        procedure :: get_pair_shells
        procedure :: get_pair_prototypes
        !
    end type am_class_pair_shell

contains

    function       get_pair_prototypes(pair,pg,ic,opts) result(pair_prototypes)
        !
        implicit none
        !
        class(am_class_pair_shell), intent(inout) :: pair
        type(am_class_symmetry)   , intent(in) :: pg
        class(am_class_unit_cell) , intent(in) :: ic
        type(am_class_options)    , intent(in) :: opts
        type(am_class_options)  :: notalk
        type(am_class_symmetry) :: vg ! stabilizer of vector v
        integer , allocatable :: pair_prototypes(:,:) ! unique indicies i and j (see below)
        integer , allocatable :: ij(:,:) ! indicies i and j (see below)
        integer , allocatable :: p(:,:) ! Z of each atom in pair
        integer , allocatable :: s(:) ! stabilizer point group identifier
        real(dp), allocatable :: d(:) ! distances of atoms
        integer , allocatable :: ind(:) ! used for multiple things
        real(dp) :: v(3)
        integer :: totalshells
        integer :: i,j,k
        ! just for stdout
        integer :: totalshells_u
        integer , allocatable :: ind_u(:)
        !
        !
        if (opts%verbosity.ge.1) call am_print_title('Determining irreducible bond prototypes')
        !
        notalk = opts
        notalk%verbosity = 0
        !
        totalshells = sum(pair%nshells)
        !
        if (opts%verbosity.ge.1) call am_print('possible bond prototypes',totalshells)
        !
        allocate(d(totalshells))
        allocate(s(totalshells))
        allocate(p(2,totalshells))
        allocate(ij(2,totalshells))
        !
        k=0
        do i = 1, ic%natoms
        do j = 1, pair%nshells(i)
            k=k+1
            ! get bond vector
            v = pair%shell(i,j)%tau(:,1)
            ! record distances
            d(k) = norm2(matmul(pair%shell(i,j)%bas,v))
            ! record indices
            ij(1:2,k) = [i,j]
            ! record sorted pair
            p(1:2,k) = [ic%Z(i), pair%shell(i,j)%Z(1)]
            if (p(1,k).gt.p(2,k)) p([1,2],k) = p([2,1],k)
            ! record stabilzier group
            call vg%get_stabilizer_group(pg=pg, v=v, opts=notalk)
            s(k) = vg%pg_identifier
            !
            ! call pair%shell(i,j)%write_poscar(file_output_poscar='shell_'//trim(int2char(i))//'_'//trim(int2char(j)))
        enddo
        enddo
        if (k.ne.totalshells) then
            call am_print('ERROR','k.ne.totalshells',flags='E')
            stop
        endif
        !
        allocate(ind(totalshells))
        !
        ! compare to figure out which pairs are unique
        ind = 0
        do i = 1, totalshells
            if (ind(i).eq.0) then
                do j = 1, totalshells
                    ! check that bond length is the same
                    if ( abs(d(i)-d(j)).lt.opts%sym_prec) then
                    ! check that atoms are the same type
                    if ( all(p(1:2,i).eq.p(1:2,j)) ) then
                    ! check that stabilizers are the same
                    if ( s(i).eq.s(j) ) then
                        ind(j) = i
                    endif
                    endif
                    endif
                enddo
            endif
        enddo
        !
        ! output
        allocate(pair_prototypes,source=ij(1:2,unique(ind)))
        !
        !
        if (opts%verbosity.ge.1) then
            ind_u = unique(ind)
            totalshells_u = size(ind_u)
            call am_print('unique bond prototypes ('//trim(int2char(totalshells_u))//')',ij(:,ind_u))
            !
            write(*,'(5x)', advance='no')
            do k = 1,totalshells_u
                i = ij(1,ind_u(k))
                j = ij(2,ind_u(k))
                write(*,'(a,"-",a)',advance='no') "("//trim(atm_symb(ic%Z(i))), trim(atm_symb(pair%shell(i,j)%Z(1)))//") "
            enddo
            write(*,*)
        endif
    end function   get_pair_prototypes

    subroutine     get_pair_shells(pair,ic,pg,sg,pair_cutoff,uc,opts)
        !
        implicit none
        !
        class(am_class_pair_shell)   , intent(inout) :: pair
        class(am_class_unit_cell), intent(in) :: ic ! irreducible cell, only determine pairs for which atleast one atom is in the primitive cell
        type(am_class_symmetry) , intent(in) :: sg ! space group
        type(am_class_symmetry) , intent(in) :: pg ! point group
        type(am_class_unit_cell), intent(in) :: uc ! unit cell
        type(am_class_options)  , intent(in) :: opts
        real(dp), intent(inout) :: pair_cutoff
        !
        type(am_class_symmetry) :: rg ! local point groupas seen by rotating the shell
        type(am_class_symmetry) :: revg ! reversal group of a typical bond in the shell v
        type(am_class_symmetry) :: vg ! stabilizer of a typical bond in the shell v
        type(am_class_unit_cell) :: sphere ! sphere containing atoms up to a cutoff
        type(am_class_options) :: notalk ! supress verbosity
        integer , allocatable :: pair_nelements(:)
        integer , allocatable :: pair_member(:,:)
        integer , allocatable :: ind_u(:)
        integer , allocatable :: ind(:)
        integer  :: npairs !  number of pairs
        integer  :: maxpairs
        real(dp) :: D(3)
        integer  :: k,i,j
        !
        ! print title
        if (opts%verbosity.ge.1) call am_print_title('Determining nearest-neighbor atomic pairs')
        !
        ! set notalk option
        notalk = opts 
        notalk%verbosity = 0
        !
        ! set pair cutoff radius (smaller than half the smallest cell dimension, lapger than the smallest distance between atoms)
        pair_cutoff = minval([norm2(uc%bas(:,:),1)/real(2,dp), pair_cutoff])
        call am_print('pair cutoff radius',pair_cutoff,' ... ')
        !
        ! get maxmimum number of pair shells
        maxpairs=0
        do i = 1, ic%natoms
            ! get center of sphere in fractional supercell coordinates
            D = matmul(matmul(inv(uc%bas),ic%bas),ic%tau(:,i))
            ! create sphere
            sphere = create_sphere(uc=uc, sphere_center=D, pair_cutoff=pair_cutoff, opts=opts )
            ! get number of pairs
            npairs = maxval(identify_pairs(sphere=sphere,pg=pg))
            if (npairs.gt.maxpairs) then
                maxpairs = npairs
            endif
        enddo
        !
        ! allocate pair space, pair pair centered on atom ic%natoms 
        allocate( pair%nshells(ic%natoms) )
        allocate( pair%shell(ic%natoms,maxpairs) )
        ! 
        do i = 1, ic%natoms
            !
            ! make a sphere containing atoms a maximum distance of a choosen atom; translate all atoms around the sphere.
            ! this distance is the real-space cut-off radius used in the construction of second-order force constants
            ! its value cannot exceed half the smallest dimension of the supercell
            !
            ! get center of sphere in fractional supercell coordinates
            D = matmul(matmul(inv(uc%bas),ic%bas),ic%tau(:,i))
            !
            ! create sphere
            sphere = create_sphere(uc=uc, sphere_center=D, pair_cutoff=pair_cutoff, opts=opts )
            !
            ! get pair identifier
            ind_u = identify_pairs(sphere=sphere,pg=pg)
            ! get number of pairs
            npairs = maxval(ind_u)
            ! get pair members [pair representative is given by pair_member(:,1)]
            pair_member = member(ind_u)
            ! get number of pair elements (essentially the oribital weights)
            pair_nelements = nelements(ind_u)
            !
            if (opts%verbosity.ge.1) then
                !
                write(*,'(" ... irreducible atom ",a," at "   ,a,",",a,",",a,  " (frac) has ",a," nearest-neighbor pairs")') &
                    & trim(int2char(i)), (trim(dbl2char(D(k),4)),k=1,3), trim(int2char(npairs))
                !
                write(*,'(5x)' ,advance='no')
                write(*,'(a5)' ,advance='no') 'shell'
                write(*,'(a6)' ,advance='no') 'i-j'
                write(*,'(a5)' ,advance='no') 'm'
                write(*,'(a8)' ,advance='no') 'stab.'
                write(*,'(a8)' ,advance='no') 'rot.'
                write(*,'(a8)' ,advance='no') 'rev.'
                write(*,'(a10)',advance='no') '|v(cart)|'
                write(*,'(a30)',advance='no') centertitle('v(cart)',30)
                write(*,'(a30)',advance='no') centertitle('v(frac)',30)
                write(*,*)
                write(*,'(5x)' ,advance='no')
                write(*,'(a5)' ,advance='no')      repeat('-',5)
                write(*,'(a6)' ,advance='no') ' '//repeat('-',5)
                write(*,'(a5)' ,advance='no') ' '//repeat('-',4)
                write(*,'(a8)', advance='no') ' '//repeat('-',7)
                write(*,'(a8)', advance='no') ' '//repeat('-',7)
                write(*,'(a8)', advance='no') ' '//repeat('-',7)
                write(*,'(a10)',advance='no') ' '//repeat('-',9)
                write(*,'(a30)',advance='no') ' '//repeat('-',29)
                write(*,'(a30)',advance='no') ' '//repeat('-',29)
                write(*,*)
                !
            endif
            !
            ! each atom has npairs, with a representative for each; save their positions
            ! representative is given by the first element of the pair_member(npair,1)
            pair%nshells(i) = npairs
            do j = 1, npairs
                !
                if (allocated(ind)) deallocate(ind)
                allocate(ind,source=pair_member(j,1:pair_nelements(j)))
                !
                call pair%shell(i,j)%initialize(bas=ic%bas,tau=sphere%tau(:,ind),Z=sphere%Z(ind))
                !
                ! print to stdout
                if (opts%verbosity.ge.1) then
                    ! determine rotations which leaves shell invariant
                    call rg%get_rotational_group(pg=pg, uc=pair%shell(i,j), opts=notalk)
                    ! determine stabilizers of a prototypical bond in shell (vector v)
                    call vg%get_stabilizer_group(pg=pg, v=pair%shell(i,j)%tau(1:3,1), opts=notalk)
                    ! determine reversal group, space symmetries, which interchange the position of atoms at the edges of the bond
                    call revg%get_reversal_group(sg=sg,v=pair%shell(i,j)%tau(1:3,1),opts=notalk)
                    write(*,'(5x)'    ,advance='no')
                    write(*,'(i5)'    ,advance='no') j
                    write(*,'(a6)'    ,advance='no') trim(atm_symb(ic%Z(i)))//'-'//trim(atm_symb(pair%shell(i,j)%Z(1)))
                    write(*,'(i5)'    ,advance='no') pair%shell(i,j)%natoms
                    write(*,'(a8)'    ,advance='no') trim(decode_pointgroup( vg%pg_identifier ))
                    write(*,'(a8)'    ,advance='no') trim(decode_pointgroup( rg%pg_identifier ))
                    write(*,'(a8)'    ,advance='no') trim(decode_pointgroup( revg%pg_identifier ))
                    write(*,'(f10.3)' ,advance='no') norm2(matmul(pair%shell(i,j)%bas,pair%shell(i,j)%tau(1:3,1)))
                    write(*,'(3f10.3)',advance='no') matmul(pair%shell(i,j)%bas,pair%shell(i,j)%tau(1:3,1))
                    write(*,'(3f10.3)',advance='no') pair%shell(i,j)%tau(1:3,1)
                    write(*,*)
                endif
                !
            enddo
            !
        enddo ! primitive cell atoms
        !
        contains
        function       create_sphere(uc,sphere_center,pair_cutoff,opts) result(sphere)
            !
            implicit none
            !
            type(am_class_unit_cell), intent(in) :: uc
            type(am_class_options)  , intent(in) :: opts
            type(am_class_unit_cell) :: sphere
            real(dp), intent(in)  :: sphere_center(3)
            real(dp), intent(in)  :: pair_cutoff
            integer , allocatable :: atoms_inside(:)
            real(dp), allocatable :: grid_points(:,:) ! used for wigner-seitz reduction
            type(am_class_options) :: notalk ! supress verbosity
            real(dp) :: bas(3,3), D(3)
            logical :: check_center
            integer :: i,j
            !
            ! set notalk option
            notalk = opts 
            notalk%verbosity = 0
            !
            ! generate grid points for wigner_seitz reduction
            D = real([-1:1],dp)
            grid_points = meshgrid(D,D,D)
            grid_points = matmul(uc%bas,grid_points)
            !
            ! create sphere instance
            bas = 2.0_dp*eye(3)
            call sphere%get_supercell(uc=uc, bscfp=bas, opts=notalk)
            D = matmul(inv(bas),sphere_center)

            ! call sphere%copy(uc=uc)
            !
            ! elements with incomplete orbits should be ignored, since they do not have enough information to build full pairs
            allocate(atoms_inside(sphere%natoms))
            atoms_inside = 0
            !
            ! filter sphere keeping only atoms iniside pair cutoff raidus
            j=0
            do i = 1, sphere%natoms
                ! turn sphere into a block with select atom at the origin
                sphere%tau(:,i) = sphere%tau(:,i) - D
                ! translate atoms to be as close to the origin as possible
                ! sphere%tau(:,i) = reduce_to_wigner_seitz(pnt=sphere%tau(:,i),grid_points=grid_points,bas=uc%bas,sym_prec=opts%sym_prec)
                sphere%tau(:,i) = modulo(sphere%tau(:,i) + 0.5_dp + opts%sym_prec, 1.0_dp) - 0.5_dp - opts%sym_prec
                ! take note of points within the predetermied pair cutoff radius
                if (norm2(matmul(sphere%bas,sphere%tau(:,i))).le.pair_cutoff + opts%sym_prec) then
                    j=j+1
                    atoms_inside(j) = i
                endif
            enddo
            call sphere%filter(indices=atoms_inside(1:j))
            !
            ! sphere is in fractional. which means, if a supercell was created by expanding the basis by a factor of 2, the fractional distance between atoms shrunk by half.
            sphere%tau = matmul(bas,sphere%tau)
            sphere%bas = matmul(inv(bas),sphere%bas)
            !
            !
            ! check that there is one atom at the origin
            check_center = .false.
            do i = 1,sphere%natoms
                if (all(abs(sphere%tau(:,i)).lt.tiny)) then
                    check_center = .true.
                    exit
                endif
            enddo
            if (check_center.eq..false.) then
                call am_print('ERROR','No atom at the origin.',flags='E')
                call am_print('sphere_center',sphere_center)
                call am_print('sphere%tau',transpose(sphere%tau))
                stop
            endif
            !
        end function   create_sphere
        function       identify_pairs(sphere,pg) result(ind_u)
            !
            use am_rank_and_sort
            !
            implicit none
            !
            type(am_class_unit_cell), intent(in) :: sphere
            type(am_class_symmetry) , intent(in) :: pg ! rotational group (space symmetries which have translational part set to zero and are still compatbile with the atomic basis)
            integer , allocatable :: PM(:,:) ! PM(uc%natoms,sg%nsyms) shows how atoms are permuted by each space symmetry operation
            real(dp), allocatable :: d(:)    ! d(npairs) array containing distances between atoms
            integer , allocatable :: indices(:)
            integer , allocatable :: ind_u(:)
            integer :: i, jj, j, k
            !
            ! write action file and get permutation map PM(sphere%natoms,sg%nsyms) which shows how space symmetries permute atomic positions
            call pg%symmetry_action(uc=sphere,oopt_PM=PM,flags='relax_pbc',opts=notalk)
            !
            ! get distance of atoms
            allocate(d(sphere%natoms))
            do i = 1, sphere%natoms
                d(i) = norm2(matmul(sphere%bas,sphere%tau(:,i)))
            enddo
            !
            ! rank atoms according to their distances
            allocate(indices(sphere%natoms))
            call rank(d,indices)
            !
            ! get pairs starting with closest atoms first
            allocate(ind_u(sphere%natoms))
            ind_u=0
            k=0
            do jj = 1, sphere%natoms
                i = indices(jj)
                if (ind_u(i).eq.0) then
                    k=k+1
                    do j = 1, pg%nsyms
                    if ((PM(i,j)).ne.0) then
                        ind_u(PM(i,j)) = k
                    endif
                    enddo
                endif
            enddo
            !
        end function   identify_pairs
    end subroutine get_pair_shells

    pure function  transform_2nd_order_force_constants(M,R,PM) result(RM)
        !
        implicit none
        !
        real(dp), intent(in) :: M(:,:,:,:) ! M(3,3,uc%natoms,uc%natoms) 
        real(dp), intent(in) :: R(3,3)  ! rotational part of space symmetry
        integer , intent(in) :: PM(:) ! ! PM(uc%natoms) permutation map contains a table of indices describing atoms space symmetry map other atoms to
        real(dp) :: RM(size(M,1),size(M,2),size(M,3),size(M,4))
        integer  :: i, j, alpha, beta, gamma, delta
        integer  :: natoms
        !
        natoms = size(PM,1)
        !
        ! second order force constants transform like 
        ! p 678 "Symmetry and condensed matter physics" Wooten
        !
        RM = 0
        !
        do i = 1, natoms
        do j = 1, natoms
        do beta  = 1,3
        do alpha = 1,3
            ! gamma -> alpha
            ! delta -> beta
            do gamma = 1,3
            do delta = 1,3
                RM(alpha,beta,PM(i),PM(j)) = RM(alpha,beta,PM(i),PM(j)) + R(alpha,gamma)*R(beta,delta)*M(gamma,delta,i,j)
            enddo
            enddo
        enddo
        enddo
        enddo
        enddo
        !
    end function   transform_2nd_order_force_constants

end module am_shells


