module am_shells

    use am_constants
    use am_stdout
    use am_unit_cell
    use am_prim_cell
    use am_irre_cell
    use am_options
    use am_mkl
    use am_atom
    use am_symmetry

    implicit none

    private

    type, public, extends(am_class_unit_cell) :: am_shell_cell
        !
        ! identifies irreducible atoms
        integer :: i
        integer :: j
        ! identifies primitive atoms
        integer :: m
        integer :: n
        ! tight binding matrix elements
        real(dp), allocatable :: Vsk(:,:)
        !
        ! <Inherited from am_class_unit_cell>
        ! real(dp) :: bas(3,3) !> column vectors a(1:3,i), a(1:3,j), a(1:3,k)
        ! integer  :: natoms   !> number of atoms
        ! real(dp), allocatable :: tau(:,:) !> tau(3,natoms) fractional atomic coordinates
        ! integer , allocatable :: Z(:) !> identifies type of element
        ! integer , allocatable :: uc_id(:) ! identifies corresponding atom in unit cell
        ! integer , allocatable :: pc_id(:) ! identifies corresponding atom in primitive cell
        ! integer , allocatable :: ic_id(:) ! identifies corresponding atom in irreducible cell
        ! </Inherited from am_class_unit_cell>
        !
    end type am_shell_cell

    type, public :: am_class_pair_shell
        ! 
        integer :: nshells ! how many shells irreducible atoms
        type(am_shell_cell), allocatable :: shell(:) ! shell(k)
        !
        integer, allocatable :: ip_id(:) ! identifies irreducible pair
        integer, allocatable :: pp_id(:) ! identifies primitive cell pair
        !
    contains
        !
        procedure :: get_primitive
        procedure :: get_irreducible
        !
    end type am_class_pair_shell

contains
    
    subroutine     get_primitive(pp,pc,pg,sg,pair_cutoff,uc,opts)
        !
        implicit none
        !
        class(am_class_pair_shell), intent(inout) :: pp ! primitive pairs
        class(am_class_unit_cell) , intent(in) :: pc ! primitive cell
        type(am_class_space_group), intent(in) :: sg ! space group
        type(am_class_point_group), intent(in) :: pg ! point group
        type(am_class_unit_cell)  , intent(in) :: uc ! unit cell
        type(am_class_options)    , intent(in) :: opts
        real(dp), intent(inout) :: pair_cutoff
        !
        character(11) :: pctype
        type(am_class_point_group) :: rotg ! local point groupas seen by rotating the shell
        type(am_class_point_group) :: revg ! reversal group of a typical bond in the shell v
        type(am_class_point_group) :: stab ! stabilizer of a typical bond in the shell v
        type(am_class_unit_cell) :: sphere ! sphere containing atoms up to a cutoff
        type(am_class_options) :: notalk ! supress verbosity
        integer , allocatable :: pair_nelements(:)
        integer , allocatable :: pair_member(:,:)
        integer , allocatable :: ind_u(:)
        integer , allocatable :: ind(:)
        integer  :: npairs !  number of pairs
        integer  :: nshells
        real(dp) :: D(3)
        integer  :: k,i,j
        !
        select type (pc)
        type is (am_class_prim_cell); pctype = 'primitive'
        type is (am_class_irre_cell); pctype = 'irreducible'
        class default
            stop 'pc type is not valid.'
        end select
        !
        ! print title
        if (opts%verbosity.ge.1) call am_print_title('Determining '//trim(pctype)//' nearest-neighbor pairs')
        !
        ! set notalk option
        notalk = opts 
        notalk%verbosity = 0
        !
        ! set pair cutoff radius (smaller than half the smallest cell dimension, larger than the smallest distance between atoms)
        pair_cutoff = minval([norm2(uc%bas(:,:),1)/real(2,dp), pair_cutoff])
        call am_print('pair cutoff radius',pair_cutoff,' ... ')
        !
        ! get maxmimum number of pair shells
        nshells=0
        do i = 1, pc%natoms
            ! get center of sphere in fractional supercell coordinates
            D = matmul(matmul(inv(uc%bas),pc%bas),pc%tau(:,i))
            ! create sphere
            sphere = create_sphere(uc=uc, sphere_center=D, pair_cutoff=pair_cutoff, opts=opts )
            ! get number of pairs
            nshells = nshells + maxval(identify_pairs(sphere=sphere,pg=pg))
        enddo
        !
        ! allocate pair space, pair centered on atom pc%natoms 
        pp%nshells = nshells
        allocate( pp%shell(pp%nshells) )
        ! 
        k = 0 ! shell index
        do i = 1, pc%natoms
            !
            ! make a sphere containing atoms a maximum distance of a choosen atom; translate all atoms around the sphere.
            ! this distance is the real-space cut-off radius used in the construction of second-order force constants
            ! its value cannot exceed half the smallest dimension of the supercell
            !
            ! get center of sphere in fractional supercell coordinates
            D = matmul(matmul(inv(uc%bas),pc%bas),pc%tau(:,i))
            !
            ! create sphere
            sphere = create_sphere(uc=uc, sphere_center=D, pair_cutoff=pair_cutoff, opts=opts )
            !
            ! identify atoms in the whole sphere that can be mapped onto each other by point symmetry operations and return the id of the pair for each atom (including the center atom)
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
                write(*,'(" ... ",a," atom ",a," at "   ,a,",",a,",",a,  " (frac) has ",a," nearest-neighbor shells")') &
                    & trim(pctype), trim(int2char(i)), (trim(dbl2char(D(j),4)),j=1,3), trim(int2char(npairs))
                !
                write(*,'(5x)' ,advance='no')
                write(*,'(a5)' ,advance='no') 'shell'
                write(*,'(a6)' ,advance='no') 'Zi-Zj'
                write(*,'(a6)' ,advance='no') 'i-j'
                write(*,'(a6)' ,advance='no') 'm-n'
                write(*,'(a5)' ,advance='no') 'm'
                write(*,'(a8)' ,advance='no') 'stab.'
                write(*,'(a8)' ,advance='no') 'rot.'
                write(*,'(a8)' ,advance='no') 'rev.'
                write(*,'(a10)',advance='no') '|v(cart)|'
                write(*,'(a30)',advance='no') centertitle('v(cart)',30)
                write(*,'(a10)',advance='no') '|v(frac)|'
                write(*,'(a30)',advance='no') centertitle('v(frac)',30)
                write(*,*)
                write(*,'(5x)' ,advance='no')
                write(*,'(a5)' ,advance='no')      repeat('-',5)
                write(*,'(a6)' ,advance='no') ' '//repeat('-',5)
                write(*,'(a6)' ,advance='no') ' '//repeat('-',5)
                write(*,'(a6)' ,advance='no') ' '//repeat('-',5)
                write(*,'(a5)' ,advance='no') ' '//repeat('-',4)
                write(*,'(a8)', advance='no') ' '//repeat('-',7)
                write(*,'(a8)', advance='no') ' '//repeat('-',7)
                write(*,'(a8)', advance='no') ' '//repeat('-',7)
                write(*,'(a10)',advance='no') ' '//repeat('-',9)
                write(*,'(a30)',advance='no') ' '//repeat('-',29)
                write(*,'(a10)',advance='no') ' '//repeat('-',9)
                write(*,'(a30)',advance='no') ' '//repeat('-',29)
                write(*,*)
                !
            endif
            !
            ! each atom has npairs, with a representative for each; save their positions
            ! representative is given by the first element of the pair_member(npair,1)
            do j = 1, npairs
                !
                if (allocated(ind)) deallocate(ind)
                allocate(ind,source=pair_member(j,1:pair_nelements(j)))
                !
                ! create the jth shell centered on primitive atom i
                k=k+1
                call pp%shell(k)%copy(uc=sphere)
                call pp%shell(k)%filter(indices=ind)
                !
                pp%shell(k)%i = pc%ic_id(i)
                pp%shell(k)%j = pp%shell(k)%ic_id(1)
                !
                pp%shell(k)%m = pc%pc_id(i)
                pp%shell(k)%n = pp%shell(k)%pc_id(1)
                !
                ! print to stdout
                if (opts%verbosity.ge.1) then
                    ! determine stabilizers of a prototypical bond in shell (vector v)
                    call stab%get_stabilizer_group(pg=pg, v=pp%shell(k)%tau(1:3,1), opts=notalk)
                    ! determine rotations which leaves shell invariant
                    call rotg%get_rotational_group(pg=pg, uc=pp%shell(k), opts=notalk)
                    ! determine reversal group, space symmetries, which interchange the position of atoms at the edges of the bond
                    call revg%get_reversal_group(sg=sg, v=pp%shell(k)%tau(1:3,1), opts=notalk)
                    write(*,'(5x)'    ,advance='no')
                    write(*,'(i5)'    ,advance='no') j
                    write(*,'(a6)'    ,advance='no') trim(atm_symb(pc%Z( pp%shell(k)%i )))//'-'//trim(atm_symb(pc%Z( pp%shell(k)%j )))
                    write(*,'(a6)'    ,advance='no') trim(int2char(pp%shell(k)%i))//'-'//trim(int2char(pp%shell(k)%j))
                    write(*,'(a6)'    ,advance='no') trim(int2char(pp%shell(k)%m))//'-'//trim(int2char(pp%shell(k)%n))
                    write(*,'(i5)'    ,advance='no') pp%shell(k)%natoms
                    write(*,'(a8)'    ,advance='no') trim(decode_pointgroup( stab%pg_id ))
                    write(*,'(a8)'    ,advance='no') trim(decode_pointgroup( rotg%pg_id ))
                    write(*,'(a8)'    ,advance='no') trim(decode_pointgroup( revg%pg_id ))
                    write(*,'(f10.3)' ,advance='no') norm2(matmul(pp%shell(k)%bas,pp%shell(k)%tau(1:3,1)))
                    write(*,'(3f10.3)',advance='no') matmul(pp%shell(k)%bas,pp%shell(k)%tau(1:3,1))
                    write(*,'(f10.3)' ,advance='no') norm2(pp%shell(k)%tau(1:3,1))
                    write(*,'(3f10.3)',advance='no') pp%shell(k)%tau(1:3,1)
                    write(*,*)
                endif
                !
            enddo
            !
        enddo ! primitive cell atoms
        ! write some definitions just to remember
        write(*,'(a5,a)') ' ... ', 'Definitions:'
        write(*,'(5x,a)') ' - i, j             : irreducible indicies'
        write(*,'(5x,a)') ' - m, n             : primitive indicies'
        write(*,'(5x,a)') ' - multiplicity (m) : pair apperences'
        write(*,'(5x,a)') ' - stabilizer (stab): point symmetries which leave a representative bond in shell invariant'
        write(*,'(5x,a)') ' - rotational (rot) : point symmetries which leave shell invariant'
        write(*,'(5x,a)') ' - reversal   (rev) : (im)proper rotational parts of space symmetries which flip bond endpoints'
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
                ! sphere%tau(:,i) = reduce_to_wigner_seitz(pnt=sphere%tau(:,i),grid_points=grid_points,bas=uc%bas,prec=opts%prec)
                sphere%tau(:,i) = modulo(sphere%tau(:,i) + 0.5_dp + opts%prec, 1.0_dp) - 0.5_dp - opts%prec
                ! take note of points within the predetermied pair cutoff radius
                ! right now, pair cutoff is in fractional. should use cartesian coordinats.
                if (norm2(matmul(sphere%bas,sphere%tau(:,i))).le.pair_cutoff + opts%prec) then
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
            type(am_class_unit_cell)  , intent(in) :: sphere
            type(am_class_point_group), intent(in) :: pg ! rotational group (space symmetries which have translational part set to zero and are still compatbile with the atomic basis)
            integer , allocatable :: P(:,:,:)
            integer , allocatable :: PM(:,:) ! PM(uc%natoms,sg%nsyms) shows how atoms are permuted by each space symmetry operation
            real(dp), allocatable :: d(:)    ! d(npairs) array containing distances between atoms
            integer , allocatable :: indices(:)
            integer , allocatable :: ind_u(:)
            integer :: i, jj, j, k
            !
            ! PM(uc%natoms,sg%nsyms) shows how atoms are permuted by each space symmetry operation
            P  = permutation_rep(seitz=pg%sym, tau=sphere%tau, flags='relax_pbc', prec=opts%prec)
            PM = permutation_map( P )
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
    end subroutine get_primitive

    subroutine     get_irreducible(ip,pp,pg,sg,ic,opts)
        !
        implicit none
        !
        class(am_class_pair_shell), intent(out) :: ip
        class(am_class_pair_shell), intent(inout) :: pp
        type(am_class_space_group), intent(in) :: sg ! space group
        type(am_class_point_group), intent(in) :: pg ! point group
        class(am_class_unit_cell) , intent(in) :: ic
        type(am_class_options)    , intent(in) :: opts
        type(am_class_point_group) :: stab
        type(am_class_point_group) :: revg
        type(am_class_point_group) :: rotg
        type(am_class_options)     :: notalk
        integer , allocatable :: ip_id(:) ! can have negative, it means bond was flipped!
        integer , allocatable :: ip_id_unique(:) ! strictly positive
        integer , allocatable :: multiplicity(:)
        integer :: i,j,k
        !
        if (opts%verbosity.ge.1) call am_print_title('Determining irreducible nearest-neighbor pairs')
        !
        ! supress output from some subroutines
        notalk = opts
        notalk%verbosity = 0
        !
        ! output shells
        if (opts%verbosity.ge.1) call am_print('pair shells',pp%nshells)
        !
        ! determine irreducible pair shells
        ip_id = identify_irreducible(pp=pp,pg=pg,ic=ic,opts=opts)
        ! call am_print('ip_id',ip_id)
        !
        ! allocate space
        ip_id_unique = unique(abs(ip_id))
        ! call am_print('ip_id_unique',ip_id_unique)
        ip%nshells = size(ip_id_unique)
        if (opts%verbosity.ge.1) call am_print('irreducible pair shells',ip%nshells)
        !
        ! create ip instance
        allocate(ip%shell(ip%nshells))
        do i = 1,ip%nshells
            ip%shell(i) = pp%shell(ip_id_unique(i))
        enddo
        !
        allocate(multiplicity(ip%nshells))
        multiplicity=0
        do k = 1, ip%nshells
        do j = 1, pp%nshells
            if (ip_id_unique(k).eq.abs(ip_id(j))) multiplicity(k) = multiplicity(k) + ip%shell(k)%natoms
        enddo
        enddo
        !
        ! write to stdout
        if (opts%verbosity.ge.1) then
            !
            write(*,'(5x)' ,advance='no')
            write(*,'(a5)' ,advance='no') 'shell'
            write(*,'(a6)' ,advance='no') 'Zi-Zj'
            write(*,'(a6)' ,advance='no') 'i-j'
            write(*,'(a6)' ,advance='no') 'm-n'
            write(*,'(a5)' ,advance='no') 'm'
            write(*,'(a8)' ,advance='no') 'stab.'
            write(*,'(a8)' ,advance='no') 'rot.'
            write(*,'(a8)' ,advance='no') 'rev.'
            write(*,'(a10)',advance='no') '|v(cart)|'
            write(*,'(a30)',advance='no') centertitle('v(cart)',30)
            write(*,'(a10)',advance='no') '|v(frac)|'
            write(*,'(a30)',advance='no') centertitle('v(frac)',30)
            write(*,*)
            write(*,'(5x)' ,advance='no')
            write(*,'(a5)' ,advance='no')      repeat('-',5)
            write(*,'(a6)' ,advance='no') ' '//repeat('-',5)
            write(*,'(a6)' ,advance='no') ' '//repeat('-',5)
            write(*,'(a6)' ,advance='no') ' '//repeat('-',5)
            write(*,'(a5)' ,advance='no') ' '//repeat('-',4)
            write(*,'(a8)', advance='no') ' '//repeat('-',7)
            write(*,'(a8)', advance='no') ' '//repeat('-',7)
            write(*,'(a8)', advance='no') ' '//repeat('-',7)
            write(*,'(a10)',advance='no') ' '//repeat('-',9)
            write(*,'(a30)',advance='no') ' '//repeat('-',29)
            write(*,'(a10)',advance='no') ' '//repeat('-',9)
            write(*,'(a30)',advance='no') ' '//repeat('-',29)
            write(*,*)
            !
        endif
        do k = 1, ip%nshells
            !
            ! print to stdout
            if (opts%verbosity.ge.1) then
                i = ip%shell(k)%i
                j = ip%shell(k)%j
                ! determine stabilizers of a prototypical bond in shell (vector v)
                call stab%get_stabilizer_group(pg=pg, v=ip%shell(k)%tau(1:3,1), opts=notalk)
                ! determine rotations which leaves shell invariant
                call rotg%get_rotational_group(pg=pg, uc=ip%shell(k), opts=notalk)
                ! determine reversal group, space symmetries, which interchange the position of atoms at the edges of the bond
                call revg%get_reversal_group(sg=sg, v=ip%shell(k)%tau(1:3,1), opts=notalk)
                write(*,'(5x)'    ,advance='no')
                write(*,'(i5)'    ,advance='no') k ! shell
                write(*,'(a6)'    ,advance='no') trim(atm_symb(ic%Z( ip%shell(k)%i )))//'-'//trim(atm_symb(ic%Z( ip%shell(k)%j )))
                write(*,'(a6)'    ,advance='no') trim(int2char(ip%shell(k)%i))//'-'//trim(int2char(ip%shell(k)%j))
                write(*,'(a6)'    ,advance='no') trim(int2char(pp%shell(k)%m))//'-'//trim(int2char(pp%shell(k)%n))
                write(*,'(i5)'    ,advance='no') multiplicity(k)
                write(*,'(a8)'    ,advance='no') trim(decode_pointgroup( stab%pg_id ))
                write(*,'(a8)'    ,advance='no') trim(decode_pointgroup( rotg%pg_id ))
                write(*,'(a8)'    ,advance='no') trim(decode_pointgroup( revg%pg_id ))
                write(*,'(f10.3)' ,advance='no') norm2(matmul(ip%shell(k)%bas,ip%shell(k)%tau(1:3,1)))
                write(*,'(3f10.3)',advance='no') matmul(ip%shell(k)%bas,ip%shell(k)%tau(1:3,1))
                write(*,'(f10.3)' ,advance='no') norm2(ip%shell(k)%tau(1:3,1))
                write(*,'(3f10.3)',advance='no') ip%shell(k)%tau(1:3,1)
                write(*,*)
            endif
        enddo
        !
        ! create map
        allocate(pp%pp_id  , source = [1:pp%nshells])
        allocate(pp%ip_id  , source = ip_id)
        allocate(ip%ip_id , source = [1:ip%nshells])
        allocate(ip%pp_id(ip%nshells))
        do i = 1, ip%nshells
            ip%pp_id(i) = ip_id_unique(i)
        enddo
        !
        if (opts%verbosity.ge.1) then
            write(*,'(a5,a)',advance='no') ' ... ', 'pair mapping (from primitive: prim->irr)'
            do i = 1, pp%nshells
                if (modulo(i,11).eq.1) then
                    write(*,*)
                    write(*,'(5x)',advance='no')
                endif
                write(*,'(a8)',advance='no') trim(int2char(i))//'->'//trim(int2char(pp%ip_id(i)))
            enddo
            write(*,*)
            !
            write(*,'(a5,a)',advance='no') ' ... ', 'pair mapping (to primitive: irr->prim)'
            do i = 1, ip%nshells
                if (modulo(i,11).eq.1) then
                    write(*,*)
                    write(*,'(5x)',advance='no')
                endif
                write(*,'(a8)',advance='no') trim(int2char(i))//'->'//trim(int2char(ip%pp_id(i)))
            enddo
            write(*,*)
            !
        endif
        !
        if (opts%verbosity.ge.1) then
            write(*,'(a5,a)') ' ... ', 'Definitions:'
            write(*,'(5x,a)') ' - i, j             : irreducible indicies'
            write(*,'(5x,a)') ' - m, n             : primitive indicies'
            write(*,'(5x,a)') ' - multiplicity (m) : pair apperences'
            write(*,'(5x,a)') ' - stabilizer (stab): point symmetries which leave a representative bond in shell invariant'
            write(*,'(5x,a)') ' - rotational (rot) : point symmetries which leave shell invariant'
            write(*,'(5x,a)') ' - reversal   (rev) : (im)proper rotational parts of space symmetries which flip bond endpoints'
            write(*,'(5x,a)') ' - negative mappings indicate that endpoints have been flipped'
        endif
        !
        contains
        function       identify_irreducible(pp,pg,ic,opts) result(ip_id)
            !
            ! negative ip_id is returned for bonds that are flipped
            !
            implicit none
            !
            class(am_class_pair_shell), intent(in) :: pp ! primitive pairs
            type(am_class_point_group), intent(in) :: pg ! point group
            class(am_class_unit_cell) , intent(in) :: ic ! irreducible cell
            type(am_class_options)    , intent(in) :: opts
            type(am_class_point_group) :: stab ! stabilizer of vector v
            type(am_class_options) :: notalk
            integer , allocatable :: Z(:,:) ! Z of each atom in pair
            integer , allocatable :: s(:) ! stabilizer point group id
            real(dp), allocatable :: d(:) ! distances of atoms
            integer , allocatable :: ip_id(:) ! can have negative, it means bond was flipped!
            logical , allocatable :: isflipped(:)
            real(dp) :: v(3)
            integer :: i,j,k
            !
            ! supress output from some subroutines
            notalk = opts
            notalk%verbosity = 0
            !
            ! set total number of pairs as the sum of the number of shells on each atom 
            ! note that totalshells is not the absolute number of pairs, it is the number of shells on all atoms; i.e. the number of unique pairs all atoms have by themselves without referencing other atoms.
            !
            ! allocate space 
            allocate(d(pp%nshells))   ! distance between pair of atoms (having atoms the same distance apart is a prerequesit for the bond to be the same)
            allocate(s(pp%nshells))   ! stabilizer group id (having the same stabilier group is a prerequesit for the bond to be the same)
            allocate(Z(2,pp%nshells)) ! atomic number of elements in pair (having the same types of atoms is a prerequesit for the bond to be the same)
            allocate(isflipped(pp%nshells)) ! indicates which pairs were flipped in the comparison, returns negative ip_id if the corresponding bond was flipped
            !
            isflipped = .false.
            do k = 1, pp%nshells
                ! get bond vector
                v = pp%shell(k)%tau(:,1)
                ! record distances
                d(k) = norm2(matmul(pp%shell(k)%bas,v))
                ! record sorted pair
                Z(1:2,k) = [ic%Z(pp%shell(k)%i), pp%shell(k)%Z(1)]
                if (Z(1,k).gt.Z(2,k)) then
                    isflipped(k) = .true.
                    Z([1,2],k) = Z([2,1],k)
                endif
                ! record stabilzier group
                call stab%get_stabilizer_group(pg=pg, v=v, opts=notalk)
                s(k) = stab%pg_id
                ! possible write shell for debuging
                ! call pair%shell(i,j)%write_poscar(file_output_poscar='shell_'//trim(int2char(i))//'_'//trim(int2char(j)))
            enddo
            !
            ! compare all the prerequisists listed above (distances, atom types, and bond stabilizer) to figure out which pairs are irreducible
            allocate(ip_id(pp%nshells))
            ip_id = 0
            do i = 1, pp%nshells
                if (ip_id(i).eq.0) then
                    do j = 1, pp%nshells
                        ! check that bond length is the same
                        if ( abs(d(i)-d(j)).lt.opts%prec) then
                        ! check that atoms are the same type
                        if ( all(Z(1:2,i).eq.Z(1:2,j)) ) then
                        ! check that stabilizers are the same
                        if ( s(i).eq.s(j) ) then
                            ip_id(j) = i
                        endif
                        endif
                        endif
                    enddo
                endif
            enddo
            !
            ! indicate which atoms have been flipped with negative sign
            do k = 1, pp%nshells
                if (isflipped(k)) ip_id(k) = -ip_id(k)
            enddo
            !
        end function   identify_irreducible
    end subroutine get_irreducible

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


    
!    subroutine     get_irreducible_force_constants(pair,ic,pg,sg,opts)
!        !
!        implicit none
!        class(am_class_pair_shell), intent(inout) :: pair
!        class(am_class_unit_cell) , intent(in) :: ic
!        type(am_class_seitz_group)   , intent(in) :: pg
!        type(am_class_seitz_group)   , intent(in) :: sg
!        type(am_class_options)    , intent(in) :: opts
!        integer , allocatable :: ip_id(:) ! unique indicies i and j (see below)
!        !
!        !
!        ip_id = pair%identify_irreducible(ic=ic,pg=pg,opts=opts)
!        
!!        do i = 1, pair%
! !       call get_symmetry_adapted_tensor(pg=pg,uc=uc,opts=opts,property='reversal',relations)
!        !
!        
!
!        !
!        !! output
!        !allocate(irrep_pairs,source=ij(1:2,unique(ip_id)))
!        !!
!        !!
!        !if (opts%verbosity.ge.1) then
!        !    ind_u = unique(ip_id)
!        !    npairs_u = size(ind_u)
!        !    call am_print('irreducible pair ('//trim(int2char(npairs_u))//')',ij(:,ind_u))
!        !    !
!        !    write(*,'(5x)', advance='no')
!        !    do k = 1,npairs_u
!        !        i = ij(1,ind_u(k))
!        !        j = ij(2,ind_u(k))
!        !        write(*,'(a,"-",a)',advance='no') "("//trim(atm_symb(ic%Z(i))), trim(atm_symb(pair%shell(i,j)%Z(1)))//") "
!        !    enddo
!        !    write(*,*)
!        !endif
!    end subroutine get_irreducible_force_constants
    
end module am_shells


