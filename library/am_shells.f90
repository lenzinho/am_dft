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

    public :: am_class_pair_shell

    type, public, extends(am_class_unit_cell) :: am_shell_cell
        ! 
        ! j th shell centered on irreducible atom i 
        integer :: i
        integer :: j
        ! <Inherited from am_class_unit_cell>
        ! real(dp) :: bas(3,3) !> column vectors a(1:3,i), a(1:3,j), a(1:3,k)
        ! integer  :: natoms   !> number of atoms
        ! real(dp), allocatable :: tau(:,:) !> tau(3,natoms) fractional atomic coordinates
        ! integer , allocatable :: Z(:) !> identifies type of element
        ! integer , allocatable :: uc_identifier(:) ! identifies corresponding atom in unit cell
        ! integer , allocatable :: pc_identifier(:) ! identifies corresponding atom in primitive cell
        ! integer , allocatable :: ic_identifier(:) ! identifies corresponding atom in irreducible cell
        ! </Inherited from am_class_unit_cell>
    end type am_shell_cell
    
    type am_class_pair_shell
        ! 
        integer, allocatable :: nshells ! how many shells irreducible atoms
        type(am_shell_cell), allocatable :: shell(:) ! shell(k)
        !
    contains
        !
        procedure :: get_pair_shells
        procedure :: get_irreducible_pair_shells
        procedure :: identify_irreducible_pair_shells
        ! procedure :: get_irreducible_force_constants
        !
    end type am_class_pair_shell

contains
    
    subroutine     get_irreducible_pair_shells(irrpair,pair,pg,sg,ic,opts)
        !
        implicit none
        !
        class(am_class_pair_shell), intent(inout) :: irrpair
        class(am_class_pair_shell), intent(in) :: pair
        type(am_class_symmetry)   , intent(in) :: pg
        type(am_class_symmetry)   , intent(in) :: sg
        class(am_class_unit_cell) , intent(in) :: ic
        type(am_class_options)    , intent(in) :: opts
        integer , allocatable :: irrpair_identifier(:) ! can have negative, it means bond was flipped!
        integer , allocatable :: irrpair_identifier_unique(:) ! strictly positive
        type(am_class_symmetry) :: stab
        type(am_class_symmetry) :: revg
        type(am_class_symmetry) :: rotg
        type(am_class_options)  :: notalk
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
        if (opts%verbosity.ge.1) call am_print('pair shells',pair%nshells)
        !
        ! determine irreducible pair shells
        irrpair_identifier = identify_irreducible_pair_shells(pair=pair,pg=pg,ic=ic,opts=opts)
        ! call am_print('irrpair_identifier',irrpair_identifier)
        !
        ! allocate space
        irrpair_identifier_unique = unique(abs(irrpair_identifier))
        ! call am_print('irrpair_identifier_unique',irrpair_identifier_unique)
        irrpair%nshells = size(irrpair_identifier_unique)
        if (opts%verbosity.ge.1) call am_print('irreducible pair shells',irrpair%nshells)
        !
        ! create irrpair instance
        allocate(irrpair%shell(irrpair%nshells))
        do i = 1,irrpair%nshells
            irrpair%shell(i) = pair%shell(irrpair_identifier_unique(i))
        enddo
        !
        allocate(multiplicity(irrpair%nshells))
        multiplicity=0
        do k = 1, irrpair%nshells
        do j = 1, pair%nshells
            if (irrpair_identifier_unique(k).eq.abs(irrpair_identifier(j))) multiplicity(k) = multiplicity(k) + irrpair%shell(k)%natoms
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
        do k = 1, irrpair%nshells
            !
            ! print to stdout
            if (opts%verbosity.ge.1) then
                i = irrpair%shell(k)%i
                j = irrpair%shell(k)%j
                ! determine stabilizers of a prototypical bond in shell (vector v)
                call stab%get_stabilizer_group(pg=pg, v=irrpair%shell(k)%tau(1:3,1), opts=notalk)
                ! determine rotations which leaves shell invariant
                call rotg%get_rotational_group(pg=pg, uc=irrpair%shell(k), opts=notalk)
                ! determine reversal group, space symmetries, which interchange the position of atoms at the edges of the bond
                call revg%get_reversal_group(sg=sg, v=irrpair%shell(k)%tau(1:3,1), opts=notalk)
                write(*,'(5x)'    ,advance='no')
                write(*,'(i5)'    ,advance='no') k ! shell
                write(*,'(a6)'    ,advance='no') trim(atm_symb(ic%Z( irrpair%shell(k)%i )))//'-'//trim(atm_symb(ic%Z( irrpair%shell(k)%j )))
                write(*,'(a6)'    ,advance='no') trim(int2char(irrpair%shell(k)%i))//'-'//trim(int2char(irrpair%shell(k)%j))
                write(*,'(i5)'    ,advance='no') multiplicity(k)
                write(*,'(a8)'    ,advance='no') trim(decode_pointgroup( stab%pg_identifier ))
                write(*,'(a8)'    ,advance='no') trim(decode_pointgroup( rotg%pg_identifier ))
                write(*,'(a8)'    ,advance='no') trim(decode_pointgroup( revg%pg_identifier ))
                write(*,'(f10.3)' ,advance='no') norm2(matmul(irrpair%shell(k)%bas,irrpair%shell(k)%tau(1:3,1)))
                write(*,'(3f10.3)',advance='no') matmul(irrpair%shell(k)%bas,irrpair%shell(k)%tau(1:3,1))
                write(*,'(3f10.3)',advance='no') irrpair%shell(k)%tau(1:3,1)
                write(*,*)
            endif
        enddo
        write(*,'(5x,a)') 'Definitions:'
        write(*,'(5x,a)') 'multiplicity (m): number of times pair appear'
        !
    end subroutine get_irreducible_pair_shells
    
    function       identify_irreducible_pair_shells(pair,pg,ic,opts) result(irrpair_identifier)
        !
        ! negative irrpair_identifier is returned for bonds that are flipped
        !
        implicit none
        !
        class(am_class_pair_shell), intent(in) :: pair
        type(am_class_symmetry)   , intent(in) :: pg
        class(am_class_unit_cell) , intent(in) :: ic
        type(am_class_options)    , intent(in) :: opts
        type(am_class_options)  :: notalk
        type(am_class_symmetry) :: stab ! stabilizer of vector v
        integer , allocatable :: Z(:,:) ! Z of each atom in pair
        integer , allocatable :: s(:) ! stabilizer point group identifier
        real(dp), allocatable :: d(:) ! distances of atoms
        integer , allocatable :: irrpair_identifier(:) ! can have negative, it means bond was flipped!
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
        allocate(d(pair%nshells))   ! distance between pair of atoms (having atoms the same distance apart is a prerequesit for the bond to be the same)
        allocate(s(pair%nshells))   ! stabilizer group identifier (having the same stabilier group is a prerequesit for the bond to be the same)
        allocate(Z(2,pair%nshells)) ! atomic number of elements in pair (having the same types of atoms is a prerequesit for the bond to be the same)
        allocate(isflipped(pair%nshells)) ! indicates which pairs were flipped in the comparison, returns negative irrpair_identifier if the corresponding bond was flipped
        !
        isflipped = .false.
        do k = 1, pair%nshells
            ! get bond vector
            v = pair%shell(k)%tau(:,1)
            ! record distances
            d(k) = norm2(matmul(pair%shell(k)%bas,v))
            ! record sorted pair
            Z(1:2,k) = [ic%Z(pair%shell(k)%i), pair%shell(k)%Z(1)]
            if (Z(1,k).gt.Z(2,k)) then
                isflipped(k) = .true.
                Z([1,2],k) = Z([2,1],k)
            endif
            ! record stabilzier group
            call stab%get_stabilizer_group(pg=pg, v=v, opts=notalk)
            s(k) = stab%pg_identifier
            ! possible write shell for debuging
            ! call pair%shell(i,j)%write_poscar(file_output_poscar='shell_'//trim(int2char(i))//'_'//trim(int2char(j)))
        enddo
        !
        ! compare all the prerequisists listed above (distances, atom types, and bond stabilizer) to figure out which pairs are irreducible
        allocate(irrpair_identifier(pair%nshells))
        irrpair_identifier = 0
        do i = 1, pair%nshells
            if (irrpair_identifier(i).eq.0) then
                do j = 1, pair%nshells
                    ! check that bond length is the same
                    if ( abs(d(i)-d(j)).lt.opts%sym_prec) then
                    ! check that atoms are the same type
                    if ( all(Z(1:2,i).eq.Z(1:2,j)) ) then
                    ! check that stabilizers are the same
                    if ( s(i).eq.s(j) ) then
                        irrpair_identifier(j) = i
                    endif
                    endif
                    endif
                enddo
            endif
        enddo
        !
        ! indicate which atoms have been flipped with negative sign
        do k = 1, pair%nshells
            if (isflipped(k)) irrpair_identifier(k) = -irrpair_identifier(k)
        enddo
        !
    end function   identify_irreducible_pair_shells

    subroutine     get_pair_shells(pair,pc,pg,sg,pair_cutoff,uc,opts)
        !
        implicit none
        !
        class(am_class_pair_shell), intent(inout) :: pair
        class(am_class_unit_cell), intent(in) :: pc ! irreducible cell or primitive cell
        type(am_class_symmetry) , intent(in) :: sg ! space group
        type(am_class_symmetry) , intent(in) :: pg ! point group
        type(am_class_unit_cell), intent(in) :: uc ! unit cell
        type(am_class_options)  , intent(in) :: opts
        real(dp), intent(inout) :: pair_cutoff
        !
        character(11) :: pctype
        type(am_class_symmetry) :: rotg ! local point groupas seen by rotating the shell
        type(am_class_symmetry) :: revg ! reversal group of a typical bond in the shell v
        type(am_class_symmetry) :: stab ! stabilizer of a typical bond in the shell v
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
        ! allocate pair space, pair pair centered on atom pc%natoms 
        pair%nshells = nshells
        allocate( pair%shell(pair%nshells) )
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
            ! identify atoms in the whole sphere that can be mapped onto each other by point symmetry operations and return the identifier of the pair for each atom (including the center atom)
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
            do j = 1, npairs
                !
                if (allocated(ind)) deallocate(ind)
                allocate(ind,source=pair_member(j,1:pair_nelements(j)))
                !
                ! create the jth shell centered on primitive atom i
                k=k+1
                call pair%shell(k)%copy(uc=sphere)
                call pair%shell(k)%filter(indices=ind)
                !
                pair%shell(k)%i = pc%ic_identifier(i)
                pair%shell(k)%j = pair%shell(k)%ic_identifier(1)
                !
                ! print to stdout
                if (opts%verbosity.ge.1) then
                    ! determine stabilizers of a prototypical bond in shell (vector v)
                    call stab%get_stabilizer_group(pg=pg, v=pair%shell(k)%tau(1:3,1), opts=notalk)
                    ! determine rotations which leaves shell invariant
                    call rotg%get_rotational_group(pg=pg, uc=pair%shell(k), opts=notalk)
                    ! determine reversal group, space symmetries, which interchange the position of atoms at the edges of the bond
                    call revg%get_reversal_group(sg=sg, v=pair%shell(k)%tau(1:3,1), opts=notalk)
                    write(*,'(5x)'    ,advance='no')
                    write(*,'(i5)'    ,advance='no') j
                    write(*,'(a6)'    ,advance='no') trim(atm_symb(pc%Z( pair%shell(k)%i )))//'-'//trim(atm_symb(pc%Z( pair%shell(k)%j )))
                    write(*,'(a6)'    ,advance='no') trim(int2char(pair%shell(k)%i))//'-'//trim(int2char(pair%shell(k)%j))
                    write(*,'(i5)'    ,advance='no') pair%shell(k)%natoms
                    write(*,'(a8)'    ,advance='no') trim(decode_pointgroup( stab%pg_identifier ))
                    write(*,'(a8)'    ,advance='no') trim(decode_pointgroup( rotg%pg_identifier ))
                    write(*,'(a8)'    ,advance='no') trim(decode_pointgroup( revg%pg_identifier ))
                    write(*,'(f10.3)' ,advance='no') norm2(matmul(pair%shell(k)%bas,pair%shell(k)%tau(1:3,1)))
                    write(*,'(3f10.3)',advance='no') matmul(pair%shell(k)%bas,pair%shell(k)%tau(1:3,1))
                    write(*,'(3f10.3)',advance='no') pair%shell(k)%tau(1:3,1)
                    write(*,*)
                endif
                !
            enddo
            !
        enddo ! primitive cell atoms
        ! write some definitions just to remember
        write(*,'(5x,a)') 'Definitions:'
        write(*,'(5x,a)') 'stabilizer (stab): point symmetries which leave a prototypical bond in shell invariant'
        write(*,'(5x,a)') 'rotational (rot) : point symmetries which leave shell invariant'
        write(*,'(5x,a)') 'reversal   (rev) : (im)proper rotational parts of space symmetries which flip bond endpoints'
        write(*,'(5x,a)') 'i-j              : irreducible atom indicies'
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
                ! sphere%tau(:,i) = reduce_to_wigner_seitz(pnt=sphere%tau(:,i),grid_points=grid_points,bas=uc%bas,sym_prec=opts%sym_prec)
                sphere%tau(:,i) = modulo(sphere%tau(:,i) + 0.5_dp + opts%sym_prec, 1.0_dp) - 0.5_dp - opts%sym_prec
                ! take note of points within the predetermied pair cutoff radius
                ! right now, pair cutoff is in fractional. should use cartesian coordinats.
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


    
!    subroutine     get_irreducible_force_constants(pair,ic,pg,sg,opts)
!        !
!        implicit none
!        class(am_class_pair_shell), intent(inout) :: pair
!        class(am_class_unit_cell) , intent(in) :: ic
!        type(am_class_symmetry)   , intent(in) :: pg
!        type(am_class_symmetry)   , intent(in) :: sg
!        type(am_class_options)    , intent(in) :: opts
!        integer , allocatable :: irrpair_identifier(:) ! unique indicies i and j (see below)
!        !
!        !
!        irrpair_identifier = pair%identify_irreducible_pair_shells(ic=ic,pg=pg,opts=opts)
!        
!!        do i = 1, pair%
! !       call get_symmetry_adapted_tensor(pg=pg,uc=uc,opts=opts,property='reversal',relations)
!        !
!        
!
!        !
!        !! output
!        !allocate(irrep_pairs,source=ij(1:2,unique(irrpair_identifier)))
!        !!
!        !!
!        !if (opts%verbosity.ge.1) then
!        !    ind_u = unique(irrpair_identifier)
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


