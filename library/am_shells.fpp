#:include "fypp_macros.fpp"
module am_shells

    use am_constants
    use am_unit_cell
    use am_irre_cell
    use am_options
    use am_mkl
    use am_atom
    use am_symmetry
    use dispmodule

    implicit none

    private

    type, public, extends(am_class_unit_cell) :: am_shell_cell
        real(dp):: center(3) ! center of shell [cart]
        integer :: i ! identifies irreducible atoms (center)
        integer :: j ! identifies irreducible atoms (shell)
        integer :: m ! identifies primitive atoms (center)
        integer :: n ! identifies primitive atoms (shell)
        integer, allocatable :: pg_id(:)   ! identifies the symmetry in the point group which takes atom tau_frac(:,1) to atom(:,i)
        type(am_class_point_group) :: rotg ! local point group as seen by rotating the shell
        type(am_class_point_group) :: revg ! reversal group of a typical bond in the shell v
        type(am_class_point_group) :: stab ! stabilizer of a typical bond in the shell v
    end type am_shell_cell

    type, public :: am_class_pair_shell
        integer :: nshells ! how many shells irreducible atoms
        type(am_shell_cell), allocatable :: shell(:) ! shell(k)
        integer, allocatable :: ip_id(:) ! identifies irreducible pair
        integer, allocatable :: pp_id(:) ! identifies primitive cell pair
    end type am_class_pair_shell

    type, public, extends(am_class_pair_shell) :: am_class_prim_pair
        contains
        procedure :: get_primitive
    end type am_class_prim_pair

    type, public, extends(am_class_pair_shell) :: am_class_irre_pair
        contains
        procedure :: get_irreducible
    end type am_class_irre_pair

contains

    subroutine     get_primitive(pp,pc,pg,pair_cutoff,opts)
        !
        use am_symmetry_tables
        !
        implicit none
        !
        class(am_class_prim_pair) , intent(inout) :: pp ! primitive pairs
        type(am_class_prim_cell)  , intent(in) :: pc ! primitive cell
        type(am_class_point_group), intent(in) :: pg ! point group
        type(am_class_options)    , intent(in) :: opts
        real(dp), intent(inout) :: pair_cutoff
        !
        type(am_class_unit_cell) :: sphere ! sphere containing atoms up to a cutoff
        integer , allocatable :: shell_nelements(:)
        integer , allocatable :: shell_member(:,:)
        integer , allocatable :: shell_id(:)
        integer , allocatable :: ind(:)
        integer  :: nshells ! number of shells
        integer  :: k,i,j,m,n
        integer  :: s
        !
        ! print title
        if (opts%verbosity.ge.1) call print_title('Primitive neighbor pairs')
        !
        ! get maxmimum number of pair shells
        nshells=0
        do i = 1, pc%natoms
            ! create sphere
            sphere = create_sphere(pc=pc, sphere_center=pc%tau_cart(:,i), pair_cutoff=pair_cutoff, opts=opts )
            ! get number of pairs
            nshells = nshells + maxval(identify_shells(sphere=sphere,pg=pg))
        enddo
        !
        ! allocate pair space, pair centered on atom pc%natoms 
        pp%nshells = nshells
        allocate( pp%shell(pp%nshells) )
        ! 
        k = 0 ! shell index
        do i = 1, pc%natoms
          	! create sphere containing atoms a maximum distance of a choosen atom (atoms are translated to as close to sphere center as possible)
            sphere = create_sphere(pc=pc, sphere_center=pc%tau_cart(:,i), pair_cutoff=pair_cutoff, opts=opts )
            ! identify atoms in the whole sphere that can be mapped onto each other by point symmetry operations and return the id of the pair for each atom (including the center atom)
            shell_id = identify_shells(sphere=sphere,pg=pg)
            ! get number of pairs
            nshells = maxval(shell_id)
            ! get number of pair elements (essentially the orbital weights)
            shell_nelements = id_nelements(shell_id)
            ! get pair members [pair representative is given by shell_member(:,1)]
            shell_member = id_member(shell_id)
            do j = 1, nshells
                if (allocated(ind)) deallocate(ind)
                allocate(ind,source=shell_member(j,1:shell_nelements(j)))
                ! create the jth shell centered on primitive atom i
                k=k+1
                call pp%shell(k)%copy(uc=sphere)
                call pp%shell(k)%filter(ind=ind)
                ! record of shell center
                pp%shell(k)%center = pc%tau_cart(:,i)
                ! get irreducible atom indices
                pp%shell(k)%i = pc%ic_id(i)
                pp%shell(k)%j = pp%shell(k)%ic_id(1)
                ! get primitive atom indices
                pp%shell(k)%m = pc%pc_id(i)
                pp%shell(k)%n = pp%shell(k)%pc_id(1)
                ! fix sign before sorting to handle nonsymorphic space group cases
                s = 1
                if (pp%shell(k)%m .lt. pp%shell(k)%n) s = -s
                ! this sort below is absolutely necessary, othewise there will be problems later when pg_id is called when building the tb hamiltonian
                call pp%shell(k)%sort_atoms(criterion=s*real(nint(pp%shell(k)%tau_cart(1,:)*1.0D3),dp),flags='descend')
                call pp%shell(k)%sort_atoms(criterion=s*real(nint(pp%shell(k)%tau_cart(2,:)*1.0D3),dp),flags='descend')
                call pp%shell(k)%sort_atoms(criterion=s*real(nint(pp%shell(k)%tau_cart(3,:)*1.0D3),dp),flags='descend')
                ! take note of point symmetry which takes atom tau_frac(:,m) to atom tau_frac(:,1)
                ! should use cart here because a unitary transformation is necessary
                ! NOTE: THIS SECTION ON FINDING THE POINT SYMMETRY MUST STRICTLY COME AFTER THE SORT
                ! ABOVE, BECAUSE THE SORT ABOVE DOES NOT REARRANGE pg_id!
                allocate(pp%shell(k)%pg_id(pp%shell(k)%natoms))
                pp%shell(k)%pg_id = 0
                do m = 1, pp%shell(k)%natoms
                    search : do n = 1, pg%nsyms
                        if (isequal(pp%shell(k)%tau_cart(:,1), matmul(pg%seitz_cart(1:3,1:3,n),pp%shell(k)%tau_cart(:,m)))) then
                            pp%shell(k)%pg_id(m) = n
                            exit search
                        endif
                    enddo search
                enddo
                if (any(pp%shell(k)%pg_id.eq.0)) stop 'ERROR: Unable to identify a symmetry operation in the shell (at least one pg_id(:) = 0)!'
                ! determine stabilizers of a prototypical bond in shell (vector v)
                call pp%shell(k)%stab%get_stabilizer_group(pg=pg, v=pp%shell(k)%tau_cart(1:3,1), opts=opts, flags='cart')
                ! determine rotations which leaves shell invariant
                call pp%shell(k)%rotg%get_rotational_group(pg=pg, uc=pp%shell(k), opts=opts, flags='cart')
                ! determine reversal group, space symmetries, which interchange the position of atoms at the edges of the bond
                call pp%shell(k)%revg%get_reversal_group(pg=pg  , v=pp%shell(k)%tau_cart(1:3,1), opts=opts, flags='cart')
            enddo
        enddo
        ! print stuff
        if (opts%verbosity.ge.1) then
            write(*,'(5a,a,a)') flare, 'pair cutoff radius = ', tostring(pair_cutoff)
            ! write the number of shells each primitive cell atom has
            do i = 1, pc%natoms
                write(*,'(" ... ",a," atom ",a," at ",a," (cart) has ",a," nearest-neighbor shells")') &
                    & 'primitive', tostring(i), tostring(pc%tau_cart(:,i),'f4.2'), tostring(nshells)
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
                do k = 1, pp%nshells
                if (i.eq.pp%shell(k)%m) then
                    write(*,'(5x)'    ,advance='no')
                    write(*,'(i5)'    ,advance='no') k
                    write(*,'(a6)'    ,advance='no') trim(atm_symb(pc%Z( pp%shell(k)%i )))//'-'//trim(atm_symb(pc%Z( pp%shell(k)%j )))
                    write(*,'(a6)'    ,advance='no') trim(int2char(pp%shell(k)%i))//'-'//trim(int2char(pp%shell(k)%j))
                    write(*,'(a6)'    ,advance='no') trim(int2char(pp%shell(k)%m))//'-'//trim(int2char(pp%shell(k)%n))
                    write(*,'(i5)'    ,advance='no') pp%shell(k)%natoms
                    write(*,'(a8)'    ,advance='no') trim(get_pg_name(get_pg_code( pp%shell(k)%stab%ps_id )))
                    write(*,'(a8)'    ,advance='no') trim(get_pg_name(get_pg_code( pp%shell(k)%rotg%ps_id )))
                    write(*,'(a8)'    ,advance='no') trim(get_pg_name(get_pg_code( pp%shell(k)%revg%ps_id )))
                    write(*,'(f10.3)' ,advance='no') norm2(pp%shell(k)%tau_cart(1:3,1))
                    write(*,'(3f10.3)',advance='no') 	   pp%shell(k)%tau_cart(1:3,1)
                    write(*,'(f10.3)' ,advance='no') norm2(pp%shell(k)%tau_frac(1:3,1))
                    write(*,'(3f10.3)',advance='no') 	   pp%shell(k)%tau_frac(1:3,1)
                    write(*,*)
                endif
                enddo
            enddo
            write(*,'(a,a)') flare, 'Definitions:'
            write(*,'(5x,a)') ' - i, j             : irreducible indicies'
            write(*,'(5x,a)') ' - m, n             : primitive indicies'
            write(*,'(5x,a)') ' - multiplicity (m) : pair apperences'
            write(*,'(5x,a)') ' - stabilizer (stab): point symmetries which leave a representative bond in shell invariant'
            write(*,'(5x,a)') ' - rotational (rot) : point symmetries which leave shell invariant'
            write(*,'(5x,a)') ' - reversal   (rev) : (im)proper rotational parts of space symmetries which flip bond endpoints'
        endif
        !
        contains
        function       create_sphere(pc,sphere_center,pair_cutoff,opts) result(sphere)
            !
            implicit none
            !
            type(am_class_prim_cell), intent(in) :: pc
            type(am_class_options)  , intent(in) :: opts
            type(am_class_unit_cell) :: sphere
            real(dp), intent(in)  :: sphere_center(3)
            real(dp), intent(in)  :: pair_cutoff
            integer , allocatable :: atoms_inside(:)
            type(am_class_options) :: notalk ! supress verbosity
        	integer , allocatable :: ind(:)
            real(dp) :: bscfp(3,3), inv_bscfp(3,3)
            real(dp) :: grid_points(3,27) ! voronoi points
            logical  :: check_center
            integer  :: i,j
            !
            ! set notalk option
            notalk = opts 
            notalk%verbosity = 0
            ! create sphere instance based on a supercell
            bscfp     = eye(3) * ceiling(2.0_dp*pair_cutoff/minval(norm2(pc%bas(:,:),1)))
            inv_bscfp = inv(bscfp)
            call sphere%get_supercell(uc=pc, bscfp=bscfp, opts=notalk)
            ! generate voronoi points [cart]
            grid_points = meshgrid([-1:1],[-1:1],[-1:1])
            grid_points = matmul(sphere%bas,grid_points)
            ! filter sphere keeping only atoms inside pair cutoff radius (elements with incomplete
            ! symmetry orbits are ignored, since they do not have enough information to build full pairs)
            allocate(atoms_inside(sphere%natoms))
            atoms_inside = 0
            j=0
            do i = 1, sphere%natoms
                ! turn sphere into a block with select atom at the origin 
                sphere%tau_cart(:,i) = sphere%tau_cart(:,i) - sphere_center
                ! translate atoms to be as close to the origin as possible
                sphere%tau_cart(:,i) = reduce_to_wigner_seitz(tau=sphere%tau_cart(:,i), grid_points=grid_points)
                ! also apply same tranformations to [frac], shift to origin
				sphere%tau_frac(:,i) = sphere%tau_frac(:,i) - matmul(matmul(inv_bscfp,pc%recbas),sphere_center)
				! also apply translation to get [frac] in the range [0,1)
            	sphere%tau_frac(:,i) = modulo(sphere%tau_frac(:,i) + 0.5_dp + opts%prec, 1.0_dp) - opts%prec - 0.5_dp
                ! take note of points within the predetermined pair cutoff radius [cart.]
                if (norm2(sphere%tau_cart(:,i)).le.pair_cutoff) then
                    j=j+1
                    atoms_inside(j) = i
                endif
            enddo
            ! revert basis back to primitive cell basis and reciprocal basis
            ! (IMPORTANT: this basis update must come after grid_points generation; otherwise problem later with identify_shell)
            sphere%bas    = matmul(inv_bscfp,sphere%bas)
            sphere%recbas = inv(sphere%bas)
            ! revert tau_frac back to primitive fractional (rather than supercell fractional)
            ! sphere is in fractional. which means, if a supercell was created by expanding the basis by 
            ! a factor of 2, the fractional distance between atoms shrunk by half.
            sphere%tau_frac = matmul(bscfp,sphere%tau_frac)
            ! correct rounding error
            where (abs(sphere%tau_cart).lt.opts%prec) sphere%tau_cart = 0
        	where (abs(sphere%tau_frac).lt.opts%prec) sphere%tau_frac = 0
            ! filter atoms outside sphere (keep ones within cutoff radius)
            call sphere%filter(ind=atoms_inside(1:j))
            ! sort atoms by increasing distance from [0,0,0]
            call sphere%sort_atoms(criterion=norm2(sphere%tau_cart,1),flags='ascend')
            ! check that there is one atom at the origin
            check_center = .false.
            do i = 1,sphere%natoms
                if (isequal(sphere%tau_frac(:,i),real([0,0,0],dp))) then
                    check_center = .true.
                    exit
                endif
            enddo
            if (check_center.eq..false.) stop 'No atom [frac] at the origin.'
            !
            do i = 1,sphere%natoms
                if (isequal(sphere%tau_frac(:,i),real([0,0,0],dp))) then
                    check_center = .true.
                    exit
                endif
            enddo
            if (check_center.eq..false.) stop 'No atom [cart] at the origin.'
            !
        end function   create_sphere
        function       identify_shells(sphere,pg) result(shell_id)
            !
            implicit none
            !
            type(am_class_unit_cell)  , intent(in) :: sphere
            type(am_class_point_group), intent(in) :: pg
            integer , allocatable :: PM(:,:)
            integer , allocatable :: shell_id(:)
            integer  :: i,j,k
            !
            ! PM(uc%natoms,sg%nsyms) shows how atoms are permuted by each space symmetry operation
            ! (IMPORTANT: Must use [cart]! The rotation operation needs to be strictly unitary.)
            ! (NOTE 	: If rotational group, rather than point group were used, can probably perform check by ommiting skip_check.)
            PM = permutation_map( permutation_rep(seitz=pg%seitz_cart, tau=sphere%tau_cart, flags='relax_pbc,skip_check', prec=opts%prec) )
            ! get pairs starting with closest atoms first
            allocate(shell_id(sphere%natoms))
            shell_id=0
            k=0
            do i = 1, sphere%natoms
                if (shell_id(i).eq.0) then
                    k=k+1
                    do j = 1, pg%nsyms
                    if ((PM(i,j)).ne.0) then
                        shell_id(PM(i,j)) = k
                    endif
                    enddo
                endif
            enddo
            !
            !  Permutation matrix for Si (first atom is in a shell by itself, second three atoms are in a shell by themselves):
			!  atom 1 [ 0.00000, 0.00000, 0.00000]:  1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1
			!  atom 2 [ 1.35750, 1.35750, 1.35750]:  2   4   3   5   0   0   0   0   0   0   5   5   3   2   4   2   4   3   0   0   0   0   0   0   0   0   0   0   5   2   2   4   2   3   0   0   0   0   0   0   0   0   5   5   4   3   4   3
			!  atom 3 [-1.35750,-1.35750, 1.35750]:  3   5   2   4   0   0   0   0   0   0   2   3   5   5   3   4   2   4   0   0   0   0   0   0   0   0   0   0   3   5   3   3   4   2   0   0   0   0   0   0   0   0   4   2   5   5   2   4
			!  atom 4 [ 1.35750,-1.35750,-1.35750]:  4   2   5   3   0   0   0   0   0   0   4   2   4   3   5   5   3   2   0   0   0   0   0   0   0   0   0   0   4   4   5   2   3   4   0   0   0   0   0   0   0   0   2   3   3   2   5   5
			!  atom 5 [-1.35750, 1.35750,-1.35750]:  5   3   4   2   0   0   0   0   0   0   3   4   2   4   2   3   5   5   0   0   0   0   0   0   0   0   0   0   2   3   4   5   5   5   0   0   0   0   0   0   0   0   3   4   2   4   3   2
			!
        end function   identify_shells
    end subroutine get_primitive

    subroutine     get_irreducible(ip,pp,pg,ic,opts)
        !
        use am_rank_and_sort
        !
        implicit none
        !
        class(am_class_irre_pair) , intent(out) :: ip
        type(am_class_prim_pair)  , intent(inout) :: pp
        type(am_class_point_group), intent(in) :: pg
        type(am_class_irre_cell)  , intent(in) :: ic
        type(am_class_options)    , intent(in) :: opts
        integer , allocatable :: ip_id(:) ! can have negative, it means bond was flipped!
        real(dp), allocatable :: d(:)
        integer , allocatable :: ind(:), rind(:)
        integer , allocatable :: pg_id_unique(:)
        character(20) :: str
        integer :: mpl
        integer :: i,j,k
        !
        if (opts%verbosity.ge.1) call print_title('Irreducible neighbor pairs')
        !
        ! determine irreducible pair shells
        ip_id = identify_irreducible(pp=pp, pg=pg, ic=ic, opts=opts)
        ! get number of shells
        ip%nshells = size(unique(abs(ip_id)))
        ! create maps irreducible pair to/from primitive pair
        allocate(pp%pp_id, source = [1:pp%nshells])
        allocate(pp%ip_id, source = ip_id)
        allocate(ip%ip_id, source = [1:ip%nshells])
        allocate(ip%pp_id(ip%nshells))
        ip%pp_id = 0
        do i = 1, pp%nshells
        if (ip%pp_id(abs(ip_id(i))).eq.0) then
            ip%pp_id(abs(ip_id(i))) = i
        endif
        enddo
        ! create ip instance
        allocate(ip%shell(ip%nshells))
        do i = 1,ip%nshells
            ip%shell(i) = pp%shell( ip%pp_id(i) )
        enddo
        ! sort irreducible pair shells based on distance 
        allocate(d(ip%nshells))
        do i = 1, ip%nshells
            d(i) = norm2(ip%shell(i)%tau_cart(:,1))
        enddo
        allocate( ind(ip%nshells))
        call rank(real(nint(d*1.0D5)),ind)
        ip%shell = ip%shell(ind)
        ip%pp_id = ip%pp_id(ind)
        allocate(rind(ip%nshells))
        rind(ind) = [1:ip%nshells]
        do i = 1, pp%nshells
            pp%ip_id(i) = rind(abs(pp%ip_id(i))) * sign(1,pp%ip_id(i))
        enddo
        ! print stdout
        if (opts%verbosity.ge.1) then
            write(*,'(a,a,a)') flare, 'primitive pair shells = '  , tostring(pp%nshells)
            write(*,'(a,a,a)') flare, 'irreducible pair shells = ', tostring(ip%nshells)
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
            do k = 1, ip%nshells
                write(*,'(5x)'    ,advance='no')
                write(*,'(i5)'    ,advance='no') k ! shell
                write(*,'(a6)'    ,advance='no') trim(atm_symb(ic%Z( ip%shell(k)%i )))//'-'//trim(atm_symb(ic%Z( ip%shell(k)%j )))
                write(*,'(a6)'    ,advance='no') tostring(ip%shell(k)%i)//'-'//tostring(ip%shell(k)%j)
                write(*,'(a6)'    ,advance='no') tostring(ip%shell(k)%m)//'-'//tostring(ip%shell(k)%n)
                write(*,'(i5)'    ,advance='no') ip%shell(k)%natoms
                write(*,'(a8)'    ,advance='no') trim(get_pg_name(get_pg_code( ip%shell(k)%stab%ps_id )))
                write(*,'(a8)'    ,advance='no') trim(get_pg_name(get_pg_code( ip%shell(k)%rotg%ps_id )))
                write(*,'(a8)'    ,advance='no') trim(get_pg_name(get_pg_code( ip%shell(k)%revg%ps_id )))
                write(*,'(f10.3)' ,advance='no') norm2(ip%shell(k)%tau_cart(1:3,1))
                write(*,'(3f10.3)',advance='no')       ip%shell(k)%tau_cart(1:3,1)
                write(*,'(f10.3)' ,advance='no') norm2(ip%shell(k)%tau_frac(1:3,1))
                write(*,'(3f10.3)',advance='no')       ip%shell(k)%tau_frac(1:3,1)
                write(*,*)
            enddo
            ! write maps
            write(*,'(a,a)',advance='no') flare, 'pair mapping (to irreducible: pp->ip)'
            call id_print_map(pp%ip_id)
            !
            write(*,'(a,a)',advance='no') flare, 'pair mapping (to primitive: ip->pp)'
            call id_print_map(ip%pp_id)
            !
            ! write definitions
            write(*,'(a,a)') flare, 'Definitions:'
            write(*,'(5x,a)') ' - i, j             : irreducible indicies'
            write(*,'(5x,a)') ' - m, n             : primitive indicies'
            write(*,'(5x,a)') ' - multiplicity (m) : pair apperences'
            write(*,'(5x,a)') ' - stabilizer (stab): point symmetries which leave a representative bond in shell invariant'
            write(*,'(5x,a)') ' - rotational (rot) : point symmetries which leave shell invariant'
            write(*,'(5x,a)') ' - reversal   (rev) : (im)proper rotational parts of space symmetries which flip bond endpoints'
            write(*,'(5x,a)') ' - negative mappings indicate that endpoints have been flipped'
        endif
        !
        ! print the character tables of all subgroups which reprsent bond stabilizers
        if (opts%verbosity.ge.1) then
            ! 
            call print_title('Stabilizer subgroup character tables')
            ! get unique point group identifiers
            allocate(pg_id_unique(ip%nshells))
            do k = 1, ip%nshells
            pg_id_unique(k) = get_pg_code( ip%shell(k)%stab%ps_id )
            enddo
            pg_id_unique = unique(pg_id_unique)
            ! loop over unique point subgroups
            do k = 1, size(pg_id_unique)
            do i = 1, ip%nshells
                if (pg_id_unique(k).eq. get_pg_code(ip%shell(i)%stab%ps_id) ) then
                    ! print point group
                    call print_title(trim(get_pg_name(pg_id_unique(k))))
                    ! left coset expansion
                    write(*,'(a,a)') flare, 'left coset representatives '//tostring(ip%shell(i)%stab%nlcrs)
                    mpl=2 ! matrices per line
                    do j = 1,ip%shell(i)%stab%nlcrs
                        str = tostring(j)
                        if ((mod(j,mpl).eq.0).or.(j.eq.ip%shell(i)%stab%nlcrs)) then
                            call disp_indent()
                            call disp(title=trim(str)//' [frac]', X=pg%seitz_frac(:,:,ip%shell(i)%stab%lcr_id(j)) ,style='underline',fmt='f5.2',zeroas='0',advance='no' , trim='no')
                            call disp(title=trim(str)//' [cart]', X=pg%seitz_cart(:,:,ip%shell(i)%stab%lcr_id(j)) ,style='underline',fmt='f5.2',zeroas='0',advance='yes', trim='no')
                        else
                            call disp_indent()
                            call disp(title=trim(str)//' [frac]', X=pg%seitz_frac(:,:,ip%shell(i)%stab%lcr_id(j)) ,style='underline',fmt='f5.2',zeroas='0',advance='no' , trim='no')
                            call disp(title=trim(str)//' [cart]', X=pg%seitz_cart(:,:,ip%shell(i)%stab%lcr_id(j)) ,style='underline',fmt='f5.2',zeroas='0',advance='no ', trim='no')
                        endif
                    enddo
                    ! print character table
                    call ip%shell(i)%stab%print_character_table()
                    ! dump debugging
                    if (debug) then
!                     call execute_command_line('mkdir -p  '//trim(debug_dir)//'/subgroups/')
!                     call ip%shell(i)%stab%debug_dump(fname= trim(debug_dir)//'/subgroups/outfile.'//trim(get_pg_name(pg_id_unique(k))))
                    endif
                    ! break loop
                    exit
                    !
                endif
            enddo
            enddo
            !
        endif
        !
        contains
        function       identify_irreducible(pp,pg,ic,opts) result(ip_id)
            !
            ! negative ip_id is returned for bonds that are flipped
            !
            implicit none
            !
            class(am_class_prim_pair) , intent(in) :: pp ! primitive pairs
            type(am_class_point_group), intent(in) :: pg ! point group
            type(am_class_irre_cell)  , intent(in) :: ic ! irreducible cell
            type(am_class_options)    , intent(in) :: opts
            integer , allocatable :: Z(:,:) ! Z of each atom in pair
            integer , allocatable :: M(:,:) !primitive atom indices
            real(dp), allocatable :: n(:) ! number of atoms in shell
            integer , allocatable :: s(:) ! stabilizer point group id
            real(dp), allocatable :: d(:) ! distances of atoms
            integer , allocatable :: ip_id(:) ! can have negative, it means bond was flipped!
            logical , allocatable :: isflipped(:)
            logical :: isfound
            real(dp) :: v(3)
            integer :: i,j,k
            !
            ! set total number of pairs as the sum of the number of shells on each atom 
            ! note that totalshells is not the absolute number of pairs, it is the number of shells on all atoms; i.e. the number of unique pairs all atoms have by themselves without referencing other atoms.
            !
            ! allocate space 
            allocate(d(pp%nshells))         ! distance between pair of atoms (having atoms the same distance apart is a prerequesit for the bond to be the same)
            allocate(s(pp%nshells))         ! rotational group should be the same
            allocate(n(pp%nshells))         ! number of atoms in shell
            allocate(Z(2,pp%nshells))       ! atomic number of elements in pair (having the same types of atoms is a prerequesit for the bond to be the same)
            allocate(M(2,pp%nshells))       ! primitive atom indices
            allocate(isflipped(pp%nshells)) ! indicates which pairs were flipped in the comparison, returns negative ip_id if the corresponding bond was flipped
            !
            isflipped = .false.
            do k = 1, pp%nshells
                ! record distances
                d(k) = norm2(pp%shell(k)%tau_cart(:,1))
                ! record number of atoms in shell
                n(k) = pp%shell(k)%natoms
                ! record sorted pair atomic number
                Z(1:2,k) = [ic%Z(pp%shell(k)%i), pp%shell(k)%Z(1)]
                if (Z(1,k).gt.Z(2,k)) then
                    Z([1,2],k) = Z([2,1],k)
                endif
                ! record primitive pair indices 
                M(1:2,k) = [pp%shell(k)%m, pp%shell(k)%n]
                ! flip pairs based on primitive cell index
                if (M(1,k).gt.M(2,k)) then
                    isflipped(k) = .true.
                    M([1,2],k) = M([2,1],k)
                endif
                ! record stabilzier group
                s(k) = get_pg_code( pp%shell(k)%rotg%ps_id )
            enddo
            !
            ! compare all the prerequisists listed above (distances, atom types, and bond stabilizer) to figure out which pairs are irreducible
            allocate(ip_id(pp%nshells))
            k = 1
            ip_id = 0
            do i = 1, pp%nshells
                ! check that shell has not already been assigned
                if (ip_id(i).eq.0) then
                    ! isfound is used to count the number of irreducible shells
                    isfound = .false.
                    do j = 1, pp%nshells
                        ! check that number of atoms in shell are the same
                        if ( abs(n(i)-n(j)).lt.opts%prec) then
                        ! check that bond length is the same
                        if ( abs(d(i)-d(j)).lt.opts%prec) then
                        ! check that atoms are the same type
                        if ( all(Z(1:2,i).eq.Z(1:2,j)) ) then
                        ! check that stabilizers are the same
                        if ( s(i).eq.s(j) ) then
                        ! if everything matches, the shells are the same
                        ! mark all shells that are the same before augmenting irrep shell counter
                            isfound = .true.
                            ip_id(j) = k
                        endif
                        endif
                        endif
                        endif
                    enddo
                    ! augment irrep shell counter
                    if (isfound) then
                        k = k+1
                    endif
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

end module am_shells


