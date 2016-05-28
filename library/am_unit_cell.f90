module am_unit_cell

    use am_constants
    use am_stdout
    use am_options
    use am_mkl
    use am_atom

    implicit none

    private

    public :: lattice_symmetries
    public :: space_symmetries_from_basis
    public :: translations_from_basis
    public :: is_symmetry_valid
    public :: deform

    type, public :: am_class_unit_cell
        real(dp) :: bas(3,3)    ! basis vectors a(1:3,i), a(1:3,j), a(1:3,k)
        real(dp) :: recbas(3,3) ! reciprocal basis vectors a(1:3,i), a(1:3,j), a(1:3,k)
        integer  :: natoms      ! number of atoms
        real(dp), allocatable :: tau_frac(:,:) ! tau_frac(3,natoms) fractional atomic coordinates
        real(dp), allocatable :: tau_cart(:,:) ! tau_frac(3,natoms) cartesean  atomic coordinates
        integer , allocatable :: Z(:)     ! identifies type of element
        integer , allocatable :: uc_id(:) ! identifies corresponding atom in unit cell
        integer , allocatable :: pc_id(:) ! identifies corresponding atom in primitive cell
        integer , allocatable :: ic_id(:) ! identifies corresponding atom in irreducible cell
    contains
        procedure :: load_poscar
        procedure :: write_poscar
        procedure :: get_supercell
        procedure :: copy
        procedure :: filter
        procedure :: initialize
    end type am_class_unit_cell


contains

    !
    ! symmetry related functions for creating space and point groups
    !

    function       lattice_symmetries(bas,prec) result(R_frac)
        !>
        !> Given a metric tensor, this function returns the (arthimetic) a-holohodry:
        !> the group of point symmetries R (fractional) that are compatible with the
        !> metric tensor).
        !>
        !> Symmetry and Condensed Matter Physics: A Computational Approach. 1 edition.
        !> Cambridge, UK; New York: Cambridge University Press, 2008. pages 278-279,
        !> Eq. 10.60 and 10.61
        !>
        !> The Mathematical Theory of Symmetry in Solids: Representation Theory
        !> for Point Groups and Space Groups. 1 edition. Oxford?: New York:
        !> Oxford University Press, 2010. page 134, table 3.9
        !>
        !> Structure of Materials: An Introduction to Crystallography, Diffraction
        !> and Symmetry. 2 edition. New York: Cambridge University Press, 2012.
        !>
        !> page 169:
        !> "If [point symmetries] are described with respect to the 'primed' basis
        !>  vectors (fractional coordinates), then once again all [point symmetries]
        !>  can be represented by matrices which contain only -1, 0, +1"
        !>
        !> page 166:
        !> "The inverse of a unitary matrix (rotation) is equal to its transpose
        !>  (this is only true in cartesian coordinates.)"
        !>
        implicit none
        ! subroutine i/o
        real(dp), intent(in) :: bas(3,3)
        real(dp), intent(in) :: prec
        real(dp), allocatable :: R_frac(:,:,:) !> point symmetries (fractional)
        !
        real(dp) :: id(3,3)
        real(dp) :: recbas(3,3)
        real(dp) :: metric(3,3)
        real(dp) :: buffer(3,3,48) ! buffer for point symmetries
        integer  :: k ! point symmetry counter
        real(dp) :: o(3,3)  ! point symmetry in fractional
        integer  :: i11, i12, i13, i21, i22, i23, i31, i32, i33 ! used to generate unitary rotational matrices
        !
        id = eye(3)
        !
        recbas = inv(bas)
        !
        metric = matmul(transpose(bas),bas)
        !
        buffer = 0.0_dp
        !
        k = 0
        !
        do i11 = -1,1
        do i12 = -1,1
        do i13 = -1,1
            do i21 = -1,1
            do i22 = -1,1
            do i23 = -1,1
                do i31 = -1,1
                do i32 = -1,1
                do i33 = -1,1
                    !
                    o(1,1:3)=real([i11,i12,i13],dp)
                    o(2,1:3)=real([i21,i22,i23],dp)
                    o(3,1:3)=real([i31,i32,i33],dp)
                    !
                    ! Check that metric is left unchanged by symmetry opreation in fractional coordinates.
                    ! i.e. that metric tensor commutes with point symmetry.
                    if ( all( abs(matmul(transpose(o),matmul(metric,o))-metric) .lt. prec ) ) then
                        k = k + 1
                        ! store point symmetry in fractional coordinates
                        buffer(:,:,k) = o
                    endif
                    !
                enddo
                enddo
                enddo
            enddo
            enddo
            enddo
        enddo
        enddo
        enddo
        !
        allocate(R_frac,source=buffer(:,:,1:k))
        !
    end function   lattice_symmetries

    pure function  translations_from_basis(tau_frac,Z,prec,flags) result(T_frac)
        !
        !> flags string can contain 'prim', 'zero', 'relax' and any combination of each
        !> if it has 'zero'  : add [0,0,0] to the T vectors returned
        !> if it has 'prim'  : add primitive basis to the T vectors returned
        !> if it has 'relax' : relax symmetry check and return all translations connecting reference atom to every other atom
        !> the addition of 'relax' is useful for determining the space group
        !> the omition of 'relax', which is the default procedure, is useful for reducing an arbitrary cell to the primitive cell
        !
        implicit none
        ! subroutine i/o
        integer , intent(in) :: Z(:)     !> Z(natoms) list identifing type of atom
        real(dp), intent(in) :: tau_frac(:,:) !> tau_frac(3,natoms) fractional atomic coordinates
        real(dp), intent(in) :: prec
        character(len=*), intent(in), optional :: flags 
        real(dp), allocatable :: T_frac(:,:)   !> T(3,nTs) translation which leaves basis invariant
        real(dp), allocatable :: wrk(:,:)    ! wrk(1:3,nTs) list of possible lattice vectors that are symmetry-compatible volume
        real(dp) :: seitz(4,4)
        real(dp) :: wrkr(3)
        integer  :: natoms ! natoms number of atoms
        integer  :: nTs ! nTs number of primitive lattice vector candidates
        logical  :: relax_symmetry
        integer  :: i, j
        !
        !
        relax_symmetry = .false.
        if (present(flags)) then
            if (index(flags,'relax').ne.0) relax_symmetry = .true.
        endif 
        !
        natoms = size(tau_frac,2)
        allocate(wrk(3,natoms-1+3)) 
        wrk = 0.0_dp
        ! choose reference atom (ideally choose an atom corresponding to the species with the fewest number of atoms in unit cell)
        i = 1
        ! search for lattice vectors using, as translational components, vectors connecting the choosen atom to other atoms of the same species
        nTs = 0
        do j = 1,natoms
            if ( i .ne. j) then
            if ( Z(i) .eq. Z(j) ) then
                ! shift to put reference atom at zero.
                wrkr(1:3) = tau_frac(1:3,j) - tau_frac(1:3,i)
                ! wrkr = modulo(wrkr+prec,1.0_dp)-prec ! added this in for testing.
                !
                if (relax_symmetry) then
                    nTs = nTs + 1
                    wrk(1:3,nTs) = wrkr
                else
                    !
                    seitz = eye(4)
                    seitz(1:3,4) = wrkr
                    !
                    if ( is_symmetry_valid(seitz_frac=seitz, Z=Z, tau_frac=tau_frac, prec=prec) ) then
                        nTs = nTs + 1
                        wrk(1:3,nTs) = wrkr
                    endif
                endif
            endif
            endif
        enddo
        if (present(flags)) then
            if (index(flags,"prim").ne.0) then
                ! add native basis to vectors
                nTs = nTs+1; wrk(1:3,nTs) = real([1,0,0],dp)
                nTs = nTs+1; wrk(1:3,nTs) = real([0,1,0],dp)
                nTs = nTs+1; wrk(1:3,nTs) = real([0,0,1],dp)
            endif
            if (index(flags,"zero").ne.0) then
                ! add origin as a possible translation
                nTs = nTs+1; wrk(1:3,nTs) = real([0,0,0],dp)
            endif
        endif
        ! allocate output
        allocate(T_frac(3,nTs))
        T_frac = wrk(1:3,1:nTs)
        !
    end function   translations_from_basis

    pure function  is_symmetry_valid(tau_frac,Z,seitz_frac,prec,flags)
        !
        ! check whether symmetry operation is valid
        !
        implicit none
        ! function i/o
        real(dp), intent(in) :: tau_frac(:,:) !> tau_frac(3,natoms) fractional atomic basis
        integer , intent(in) :: Z(:)     !> Z(natoms) list identify type of atom
        real(dp), intent(in) :: seitz_frac(4,4)
        real(dp), intent(in) :: prec
        character(*), intent(in), optional :: flags
        !
        real(dp), allocatable :: tau_frac_internal(:,:)
        real(dp), allocatable :: tau_frac_ref(:,:) ! what to compare to
        logical  :: isexact    ! if present, do not apply mod. ; useful for determining stabilizers
        logical  :: iszero     ! if present, returns that the symmetry is valid if it reduces the tau_frac to [0,0,0]. useful for determing which symmetries flip bonds
        real(dp) :: tau_frac_rot(3) ! rotated 
        logical  :: is_symmetry_valid
        integer  :: natoms
        integer  :: i, j, m
        logical  :: overlap_found
        !
        natoms = size(tau_frac,2)
        !
        !
        ! load flag options
        isexact = .false.
        iszero  = .false.
        if (present(flags)) then
        if (index(flags,'exact').ne.0) isexact = .true.
        if (index(flags,'zero').ne.0)  iszero  = .true.
        endif
        !
        ! start main part of function here
        allocate(tau_frac_internal, source=tau_frac)
        if (.not.isexact) then
        do i = 1, natoms
            tau_frac_internal(:,i) = modulo(tau_frac_internal(:,i)+prec,1.0_dp)-prec
        enddo
        endif
        !
        ! tau_frac ref is what the rotated vector is compared to
        if (iszero) then
            tau_frac_ref(:,i) = 0.0_dp
        else
            ! default behavior: just compare it to what it was originally before it was rotated
            allocate(tau_frac_ref,source=tau_frac_internal)
        endif
        !
        m = 0
        do i = 1, natoms
            ! apply symmetry operation
            tau_frac_rot(1:3) = matmul(seitz_frac(1:3,1:3),tau_frac_internal(1:3,i)) + seitz_frac(1:3,4)
            ! reduce rotated+translated point to unit cell
            if (.not.isexact) tau_frac_rot(1:3) = modulo(tau_frac_rot(1:3)+prec,1.0_dp)-prec
            ! check that newly created point matches something already present
            overlap_found = .false.
            check_overlap : do j = 1,natoms
                if (Z(i) .eq. Z(j)) then
                    if (all(abs(tau_frac_rot(1:3)-tau_frac_internal(1:3,j)).lt.prec)) then
                        m = m + 1
                        overlap_found = .true.
                        exit check_overlap
                    endif
                endif
            enddo check_overlap
            !
            if ( .not. overlap_found ) then
                is_symmetry_valid = .false.
                return
            endif
            !
        enddo
        !
        if (m .eq. natoms) then
            is_symmetry_valid = .true.
            return
        endif
    end function   is_symmetry_valid

    function       space_symmetries_from_basis(uc,opts) result(seitz_frac)
        !
        !The program first finds the primitive unit cell by looking for
        !additional lattice vectors within the unit cell given in the input.
        !
        ! It chooses one of the atoms and tests each vector from that atom to
        !every other atom of the same type in the unit cell. If that vector
        !will take us from each atom in the unit cell to another atom of the
        !same type (not necessarily in the same unit cell), then we have
        !found a new lattice vector.
        !
        ! This process is repeated until we have
        !found a complete list of lattice vectors. From that list, we find three
        !which can serve as primitive basis vectors ti of the lattice, i.e. every
        !lattice vector can be expressed as an integer linear combination of
        !these basis vectors.
        !
        use am_rank_and_sort
        !
        implicit none
        ! subroutine i/o
        type(am_class_unit_cell), intent(in) :: uc
        type(am_class_options), intent(in) :: opts
        real(dp), allocatable :: seitz_frac(:,:,:)
        real(dp), allocatable :: seitz_probe(:,:)
        real(dp), allocatable :: T(:,:)
        real(dp), allocatable :: R(:,:,:)
        real(dp), allocatable :: wrk_frac(:,:,:)
        real(dp) :: shift(3)
        real(dp) :: T_shifted(3)
        integer :: nTs
        integer :: nRs
        integer :: i, j, m
        !
        R = lattice_symmetries(bas=uc%bas,prec=opts%prec)
        nRs = size(R,3)
        !
        if (opts%verbosity.ge.1) call am_print('possible (im-)proper rotations',nRs,' ... ')
        !
        T = translations_from_basis(tau_frac=uc%tau_frac, Z=uc%Z, prec=opts%prec, flags='zero,relax')
        nTs = size(T,2)
        !
        if (opts%verbosity.ge.1) call am_print('possible translations',nTs,' ... ')
        !
        if (opts%verbosity.ge.1) then
            call am_print_two_matrices_side_by_side(name='translations',&
                Atitle='fractional',A=transpose(T),&
                Btitle='cartesian' ,B=transpose(matmul(uc%bas,T)),&
                iopt_emph=' ... ',iopt_teaser=.true.)
        endif 
        !
        allocate(wrk_frac(4,4,nTs*nRs))
        wrk_frac = 0.0_dp
        m=0
        do i = 1, nRs
        do j = 1, nTs
            ! Shift tau_frac to put the first atom at the origin. This is important because the
            ! rotational part of the symmetry operation must leave at least one atom in it's
            ! original location. If the atom already has a displacement, then it will be rotated
            ! to another position.
            if (.not.any(all(abs(uc%tau_frac).lt.opts%prec,2))) then
                shift = uc%tau_frac(:,1)
                ! tau_frac' = R*tau_frac + T
                ! tau_frac' - s = R*(tau_frac-s) + T
                ! tau_frac' = R*tau_frac + (R*S+T-s) = R*tau_frac + T_shifted ! sign is wrong here for some reason.
                ! thus, T_shifted = T - R*s + s
                T_shifted = T(:,j) - matmul(R(:,:,i),shift) + shift
                ! reduce to primitive
                T_shifted = modulo(T_shifted+opts%prec,1.0_dp)-opts%prec
            else
                T_shifted = T(:,j)
            endif
            !
            seitz_probe = eye(4)
            seitz_probe(1:3,1:3) = R(1:3,1:3,i)
            seitz_probe(1:3,4) = T_shifted
            !
            if ( is_symmetry_valid(tau_frac=uc%tau_frac, Z=uc%Z, seitz_frac=seitz_probe, prec=opts%prec)) then
                m = m + 1
                wrk_frac(:,:,m)=seitz_probe
            endif
        enddo
        enddo
        !
        allocate(seitz_frac(4,4,m))
        seitz_frac = wrk_frac(1:4,1:4,1:m)
        !
    end function   space_symmetries_from_basis

    !
    ! functions which operate on uc or similar derived types
    !

    subroutine      load_poscar(uc,opts)
        ! 
        ! Reads the poscar file, look below: code is short and self explanatory.
        !
        use am_vasp_io
        !
        implicit none
        !
        class(am_class_unit_cell), intent(inout) :: uc
        type(am_class_options), intent(in) :: opts
        ! ommiting these parameters in favor of Z
        integer :: nspecies
        character(:), allocatable :: symb(:)
        integer, allocatable :: atype(:)
        integer :: i
        !
        call read_poscar(&
            bas=uc%bas,&
            natoms=uc%natoms,&
            nspecies=nspecies,&
            symb=symb,&
            tau_frac=uc%tau_frac,&
            atype=atype,&
            iopt_filename=opts%poscar,&
            iopt_verbosity=opts%verbosity)
        !
        allocate(uc%Z(uc%natoms))
        do i = 1, uc%natoms
            uc%Z(i) = atm_Z(trim(symb(atype(i))))
        enddo
        !
        allocate(uc%tau_cart(3,uc%natoms))
        uc%tau_cart = matmul(uc%bas,uc%tau_frac)
        !
        uc%recbas=inv(uc%bas)
        !
    end subroutine  load_poscar

    subroutine      write_poscar(uc,file_output_poscar)
        ! 
        ! Reads the poscar file, look below: code is short and self explanatory.
        !
        use am_vasp_io, write_poscar_internal => write_poscar
        !
        implicit none
        !
        class(am_class_unit_cell), intent(in) :: uc
        character(*), intent(in) :: file_output_poscar
        !
        integer, allocatable :: unique_Z(:)
        integer :: nspecies
        character(:), allocatable :: symb(:)
        integer, allocatable :: atype(:)
        integer :: i
        !
        if (uc%natoms.lt.1) then
            call am_print('ERROR','Number of atoms < 1.',flags='E')
            stop
        endif
        !
        ! nspecies
        unique_Z = unique(uc%Z)
        nspecies = size(unique_Z)
        if (nspecies.lt.1) then
            call am_print('ERROR','Number of atomic species < 1.',flags='E')
            stop
        endif
        !
        ! symbs
        allocate(character(500)::symb(nspecies))
        do i = 1, nspecies
            symb(i) = atm_symb(unique_Z(i))
        enddo
        !
        ! atype
        allocate(atype(uc%natoms))
        do i = 1, nspecies
            where(uc%Z(:).eq.unique_Z(i)) atype = i
        enddo
        !
        call write_poscar_internal(&
            bas=uc%bas,&
            natoms=uc%natoms,&
            nspecies=nspecies,&
            symb=symb,&
            tau_frac=uc%tau_frac,&
            atype=atype,&
            iopt_filename=file_output_poscar)
        !
    end subroutine  write_poscar

    subroutine      get_supercell(sc,uc,bscfp,opts)
        !
        ! bscfp = supercell basis in fractional primitive lattice coordinates
        ! for primitive -> conventional cell, use: bscfp = inv(centering_matrix)
        !
        implicit none
        !
        class(am_class_unit_cell), intent(out) :: sc
        class(am_class_unit_cell), intent(in)  :: uc
        real(dp), intent(in) :: bscfp(3,3)
        real(dp) :: inv_bscfp(3,3)
        real(dp) :: tau_frac_wrk(3)
        integer :: i1, i2, i3, j, m, i
        type(am_class_options), intent(in) :: opts
        !
        if (opts%verbosity.ge.1) call am_print_title('Expanding to supercell')
        !
        inv_bscfp = inv(bscfp)
        !
        sc%bas    = matmul(uc%bas,bscfp)
        sc%recbas = inv(sc%bas)
        !
        if ((mod(det(bscfp)+opts%prec,1.0_dp)-opts%prec).gt.opts%prec) then
            call am_print('ERROR','The determinant of the supercell basis in primitive fractional coordinates must be an integer.',flags='E')
            call am_print('bscfp',bscfp)
            call am_print('det(bscfp)',det(bscfp))
            call am_print('check',modulo(det(bscfp),1.0_dp).lt.tiny)
            write(*,*) opts%prec
            write(*,*) mod(det(bscfp),1.0_dp)
            write(*,*) mod(det(bscfp)+opts%prec,1.0_dp)-opts%prec
            stop
        endif
        !
        sc%natoms = nint(det(bscfp)*uc%natoms)
        !
        allocate(sc%tau_frac(3,sc%natoms))
        allocate(sc%Z(sc%natoms))
        allocate(sc%uc_id(sc%natoms)) ! map
        !
        sc%tau_frac = 1.0D5
        m=0
        map : do i1 = 0, nint(sum(abs(bscfp(:,1)))) ! sc%natoms
              do i2 = 0, nint(sum(abs(bscfp(:,2)))) ! sc%natoms
              do i3 = 0, nint(sum(abs(bscfp(:,3)))) ! sc%natoms
                  do j = 1, uc%natoms
                      !
                      tau_frac_wrk = uc%tau_frac(1:3,j)
                      ! atomic basis in primitive fractional
                      tau_frac_wrk = tau_frac_wrk+real([i1,i2,i3],dp)
                      ! convert to supercell fractional
                      tau_frac_wrk = matmul(inv_bscfp,tau_frac_wrk)
                      ! reduce to primitive supercell
                      tau_frac_wrk = modulo(tau_frac_wrk+opts%prec,1.0_dp) - opts%prec
                      !
                      if (.not.issubset(sc%tau_frac,tau_frac_wrk,opts%prec)) then
                          m = m+1 
                          sc%tau_frac(1:3,m) = tau_frac_wrk
                          sc%Z(m) = uc%Z(j)
                          sc%uc_id(m) = j
                      endif
                      if (m.eq.sc%natoms) then
                          exit map
                      endif
                  enddo
              enddo
              enddo
        enddo map
        ! correct basic rounding error
        where(abs(sc%tau_frac).lt.opts%prec) sc%tau_frac=0.0_dp 
        ! get cartesaian atomic basis
        allocate(sc%tau_cart(3,uc%natoms))
        sc%tau_cart = matmul(sc%bas,sc%tau_frac)
        ! print stdout
        !
        if (opts%verbosity.ge.1) then
            !
            call am_print_two_matrices_side_by_side(name='primitive basis',&
                Atitle='fractional (supercell)',A=inv_bscfp,&
                Btitle='cartesian'             ,B=uc%bas,&
                iopt_emph=' ... ',iopt_teaser=.true.)
            !
            call am_print_two_matrices_side_by_side(name='supercell basis',&
                Atitle='fractional (primitive)',A=bscfp,&
                Btitle='cartesian'             ,B=sc%bas,&
                iopt_emph=' ... ',iopt_teaser=.true.)
            !
            call am_print('number of atoms (supercell)',sc%natoms,' ... ')
            !
            call am_print('number of atoms (primitive)',uc%natoms,' ... ')
            !
            call am_print_two_matrices_side_by_side(name='primitive atomic basis',&
                Atitle='fractional (primitive)',A=transpose(uc%tau_frac),&
                Btitle='cartesian'             ,B=transpose(uc%tau_cart),&
                iopt_emph=' ... ',iopt_teaser=.true.)
            !
            call am_print_two_matrices_side_by_side(name='supercell atomic basis',&
                Atitle='fractional (supercell)',A=transpose(sc%tau_frac),&
                Btitle='cartesian'             ,B=transpose(sc%tau_cart),&
                iopt_emph=' ... ',iopt_teaser=.true.)
        endif
        !
        ! <MAP>
        ! mapping of each atom in the sphere to the corresponding unit cell atom was already performed in the creation of the supercell
        ! now map each supercell atom the corresponding primitive atom if the mapping from uc->pc was already established
        if (allocated(uc%pc_id)) then
            allocate(sc%pc_id(sc%natoms))
            do i = 1, sc%natoms
                sc%pc_id(i) = uc%pc_id( sc%uc_id(i) )
            enddo
            !
            if (opts%verbosity.ge.1) then
                write(*,'(a5,a)',advance='no') ' ... ', 'atomic mapping (to primitive cell: sc->pc)'
                do i = 1, sc%natoms
                    if (modulo(i,10).eq.1) then
                        write(*,*)
                        write(*,'(5x)',advance='no')
                    endif
                    write(*,'(a8)',advance='no') trim(int2char(i))//'->'//trim(int2char(sc%pc_id(i)))
                enddo
                write(*,*)
            endif
        endif
        ! now map each supercell atom the corresponding irreducible cell atom if the mapping from uc->ic was already established
        if (allocated(uc%ic_id)) then
            allocate(sc%ic_id(sc%natoms))
            do i = 1, sc%natoms
                sc%ic_id(i) = uc%ic_id( sc%uc_id(i) )
            enddo
            !
            if (opts%verbosity.ge.1) then
                write(*,'(a5,a)',advance='no') ' ... ', 'atomic mapping (to irreducible cell: sc->ic)'
                do i = 1, sc%natoms
                    if (modulo(i,10).eq.1) then
                        write(*,*)
                        write(*,'(5x)',advance='no')
                    endif
                    write(*,'(a8)',advance='no') trim(int2char(i))//'->'//trim(int2char(sc%ic_id(i)))
                enddo
                write(*,*)
            endif
        endif
        ! </MAP>
    end subroutine  get_supercell

    pure subroutine copy(cp,uc)
        !
        implicit none
        !
        class(am_class_unit_cell), intent(inout) :: cp
        type(am_class_unit_cell) , intent(in) :: uc
        !
        cp%natoms = uc%natoms
        cp%bas    = uc%bas
        cp%recbas = uc%recbas
        allocate(cp%tau_frac, source=uc%tau_frac)
        allocate(cp%tau_cart, source=uc%tau_cart)
        allocate(cp%Z       , source=uc%Z)
        if (allocated(uc%uc_id)) allocate(cp%uc_id, source=uc%uc_id)
        if (allocated(uc%pc_id)) allocate(cp%pc_id, source=uc%pc_id)
        if (allocated(uc%ic_id)) allocate(cp%ic_id, source=uc%ic_id)
        !
    end subroutine  copy
   
    subroutine      filter(uc,ind)
        !
        implicit none
        !
        class(am_class_unit_cell), intent(inout) :: uc
        integer , intent(in)  :: ind(:)
        !
        ! count how many atoms should be passed
        uc%natoms = size(ind)
        ! filter stuff
        if (allocated(uc%tau_frac)) call dvfilter(values=uc%tau_frac, ind=ind)
        if (allocated(uc%tau_cart)) call dvfilter(values=uc%tau_cart, ind=ind)
        if (allocated(uc%Z))        call  ifilter(values=uc%Z       , ind=ind)
        if (allocated(uc%uc_id))    call  ifilter(values=uc%uc_id   , ind=ind)
        if (allocated(uc%pc_id))    call  ifilter(values=uc%pc_id   , ind=ind)
        if (allocated(uc%ic_id))    call  ifilter(values=uc%ic_id   , ind=ind)
        !
        contains
        subroutine     ifilter(values,ind)
            !
            implicit none
            !
            integer , allocatable, intent(inout) :: values(:)
            integer , intent(in)    :: ind(:)
            integer , allocatable   :: Z(:) ! temporary variable
            integer :: i, ninds
            !
            ninds = size(ind)
            ! setup temporary variable
            allocate(Z,source=values)
            ! clear original variable
            deallocate(values)
            ! reallocate space
            allocate(values(ninds))
            ! transfer content through filter
            do i = 1, ninds
                values(i) = Z(ind(i))
            enddo
            ! deallocate
            deallocate(Z)
        end subroutine ifilter
        subroutine     dvfilter(values,ind)
            !
            implicit none
            !
            real(dp), allocatable, intent(inout) :: values(:,:)
            integer , intent(in)    :: ind(:)
            real(dp), allocatable   :: Z(:,:) ! temporary variable
            integer :: i, ninds, d1
            !
            ninds = size(ind)
            d1    = size(values,1)
            ! setup temporary variable
            allocate(Z,source=values)
            ! clear original variable
            deallocate(values)
            ! reallocate space
            allocate(values(d1,ninds))
            ! transfer content through filter
            do i = 1, ninds
                values(:,i) = Z(:,ind(i))
            enddo
            ! deallocate
            deallocate(Z)
        end subroutine dvfilter
    end subroutine  filter

    subroutine      initialize(uc,bas,tau_frac,Z,uc_id,pc_id,ic_id)
        !
        implicit none
        !
        class(am_class_unit_cell) , intent(out) :: uc
        real(dp), intent(in) :: bas(3,3)
        real(dp), intent(in) :: tau_frac(:,:)
        integer , intent(in) :: Z(:)
        integer , intent(in), optional :: uc_id(:)
        integer , intent(in), optional :: pc_id(:)
        integer , intent(in), optional :: ic_id(:)
        !
        uc%bas    = bas
        uc%recbas = inv(bas)
        !
        uc%natoms = size(tau_frac,2)
        !
        allocate(uc%tau_frac,source=tau_frac)
        allocate(uc%tau_cart,source=matmul(bas,tau_frac))
        !
        allocate(uc%Z,source=Z)
        !
        if (present(uc_id)) allocate(uc%uc_id,source=uc_id)
        if (present(pc_id)) allocate(uc%pc_id,source=pc_id)
        if (present(ic_id)) allocate(uc%ic_id,source=ic_id)
        !
    end subroutine  initialize

    !
    ! deform (strain related stuff)
    !

    subroutine      deform(def,uc,opts)
        !
        implicit none
        !
        type(am_class_unit_cell), allocatable, intent(out) :: def(:)
        type(am_class_unit_cell), intent(in) :: uc
        type(am_class_options)  , intent(in) :: opts
        real(dp), allocatable :: bas_def(:,:,:)
        integer :: i
        !
        if (opts%verbosity.ge.1) call am_print_title('Elastically deforming structure')
        !
        bas_def = apply_elastic_deformation(bas=uc%bas,deformation_code=opts%deformation_code,strain_max=opts%maxstrain,nstrains=opts%nstrains,iopt_verbosity=opts%verbosity)
        !
        allocate(def(opts%nstrains))
        !
        do i = 1,opts%nstrains
            !
            call def(i)%copy(uc)
            !
            def(i)%bas = bas_def(:,:,i)
            !
        enddo
        !
    end subroutine  deform

    function        apply_elastic_deformation(bas,deformation_code,strain_max,nstrains,iopt_verbosity) result(bas_def)
        !
        !
        ! needs editing to conform. make this a "pure function" , push all outputs except errors into its caller.
        !
        ! creates at bas_def(3,3,nTs) matrix containing the nTs deformed primitive basises.
        !
        ! bas_def = apply_elastic_deformation(bas=uc%bas,deformation_code=6,strain_max=0.1_dp,nstrains=10,iopt_verbosity=verbosity)
        !
        implicit none
        !
        real(dp), intent(in) :: bas(3,3)
        integer , intent(in) :: deformation_code
        real(dp), intent(in) :: strain_max
        integer , intent(in) :: nstrains
        real(dp), allocatable :: bas_def(:,:,:)
        !
        character(len=6) :: dc ! deformation code string
        integer :: i, j ! loop variable
        real(dp) :: eta_step ! the deformation step size
        real(dp) :: eta ! the magnitude of the deformation
        real(dp) :: e(6) ! strain voigt
        real(dp) :: eta_matrix(3,3) ! deformation matrix (lagrangian)
        real(dp) :: eps_matrix(3,3) ! deformation matrix (true)
        real(dp), allocatable :: def_mat(:,:,:) ! deformation matrix
        !
        integer, intent(in), optional :: iopt_verbosity
        integer :: verbosity
        if ( present(iopt_verbosity) ) then
            verbosity = iopt_verbosity
        else
            verbosity = 1
        endif
        !
        if (nstrains .le. 1) then
            call am_print('ERROR','More than one strain point is required.',flags='E')
            stop
        endif
        !
        if (strain_max .gt. 0.1) then
            call am_print('ERROR','Strain value too large.',flags='E')
            stop
        endif
        !
        dc = 'XXXXXX'
        select case(deformation_code)
            case(0);  dc='EEE000'
            case(1);  dc='E00000'
            case(2);  dc='0E0000'
            case(3);  dc='00E000'
            case(4);  dc='000E00'
            case(5);  dc='0000E0'
            case(6);  dc='00000E'
            case(7);  dc='000EEE'
            case(8);  dc='EE0000'
            case(9);  dc='Ee0000'
            case(10); dc='EEEEEE'
            case(11); dc='123456'
            case(12); dc='d14t6c'
            case(13); dc='3cu62q'
            case(14); dc='qs51t2'
            case(15); dc='546dut'
            case(16); dc='s3d541'
            case(20); dc='000EE0'
            case(40); dc='EEd000'
            case(41); dc='Ee0000'
            case default
            call am_print('ERROR','Deformation code not found. Select one of the following:',flags='E')
            write(*,'(5x,a)') "  0 =>  ( eta,  eta,  eta,    0,    0,    0)  | volume strain "
            write(*,'(5x,a)') "  1 =>  ( eta,    0,    0,    0,    0,    0)  | linear strain along x "
            write(*,'(5x,a)') "  2 =>  (   0,  eta,    0,    0,    0,    0)  | linear strain along y "
            write(*,'(5x,a)') "  3 =>  (   0,    0,  eta,    0,    0,    0)  | linear strain along z "
            write(*,'(5x,a)') "  4 =>  (   0,    0,    0,  eta,    0,    0)  | yz shear strain"
            write(*,'(5x,a)') "  5 =>  (   0,    0,    0,    0,  eta,    0)  | xz shear strain"
            write(*,'(5x,a)') "  6 =>  (   0,    0,    0,    0,    0,  eta)  | xy shear strain 2C44"
            write(*,'(5x,a)') "  7 =>  (   0,    0,    0,  eta,  eta,  eta)  | shear strain along (111)"
            write(*,'(5x,a)') "  8 =>  ( eta,  eta,    0,    0,    0,    0)  | xy in-plane strain "
            write(*,'(5x,a)') "  9 =>  ( eta, -eta,    0,    0,    0,    0)  | xy in-plane shear strain"
            write(*,'(5x,a)') " 10 =>  ( eta,  eta,  eta,  eta,  eta,  eta)  | global strain"
            write(*,'(5x,a)') " 11 =>  (  1e,   2e,   3e,   4e,   5e,   6e)  | universal strain 1"
            write(*,'(5x,a)') " 12 =>  ( -2e,   1e,   4e,  -3e,   6e,  -5e)  | universal strain 2"
            write(*,'(5x,a)') " 13 =>  (  3e,  -5e,  -1e,   6e,   2e,  -4e)  | universal strain 3"
            write(*,'(5x,a)') " 14 =>  ( -4e,  -6e,   5e,   1e,  -3e,   2e)  | universal strain 4"
            write(*,'(5x,a)') " 15 =>  (  5e,   4e,   6e,  -2e,  -1e,  -3e)  | universal strain 5"
            write(*,'(5x,a)') " 16 =>  ( -6e,   3e,  -2e,   5e,  -4e,   1e)  | universal strain 6"
            write(*,'(5x,a)') " 20 =>  (   0,    0,    0,  eta,  eta,    0)  | yz zx shear"
            write(*,'(5x,a)') " 40 =>  ( eta,  eta,  -2e,    0,    0,    0)  | tetragonal for 3/2(c11-c12)"
            write(*,'(5x,a)') " 41 =>  ( eta, -eta,    0,    0,    0,    0)  | orthorhombic for (c11-c12)"
            stop
        end select
        !
        eta_step = 2.0_dp*strain_max/(nstrains-1.0_dp)
        !
        allocate(bas_def(3,3,nstrains))
        allocate(def_mat(3,3,nstrains))
        !
        if (verbosity.ge.1) call am_print('Number of deformations',nstrains,' ... ')
        !
        do i = 1,nstrains
            ! set deformation magnitude
            eta = (i-1)*eta_step - strain_max
            if (abs(eta) < 1.0E-6_dp) eta=1.0E-6_dp ! why is this here?
            ! convert string to deformation
            e(1:6)=0.0_dp
            do j = 1,6
                select case (dc(j:j))
                    case('E'); e(j)= eta
                    case('e'); e(j)=-eta
                    case('0'); e(j)=0.0_dp
                    case('1'); e(j)= 1.0_dp*eta
                    case('2'); e(j)= 2.0_dp*eta
                    case('3'); e(j)= 3.0_dp*eta
                    case('4'); e(j)= 4.0_dp*eta
                    case('5'); e(j)= 5.0_dp*eta
                    case('6'); e(j)= 6.0_dp*eta
                    case('u'); e(j)=-1.0_dp*eta
                    case('d'); e(j)=-2.0_dp*eta
                    case('t'); e(j)=-3.0_dp*eta
                    case('q'); e(j)=-4.0_dp*eta
                    case('c'); e(j)=-5.0_dp*eta
                    case('s'); e(j)=-6.0_dp*eta
                    case default
                        call am_print('ERROR','Deformation not allowed',flags='E')
                        stop
                end select
            enddo
            ! build deformation matrix
            eta_matrix(1:3,1) = [e(1)       ,e(6)/2.0_dp,e(5)/2.0_dp]
            eta_matrix(1:3,2) = [e(6)/2.0_dp,e(2)       ,e(4)/2.0_dp]
            eta_matrix(1:3,3) = [e(5)/2.0_dp,e(4)/2.0_dp,e(3)       ]
            !
            ! Solve for eps (true strain) from eta (lagrangian strain) by iteration
            !
            ! "ElaStic: A tool for calculating second-order elastic constants from first principles"
            ! Eq. 2 : tau = det(1+eps)(1+eps)^-1 . sigma . (1+eps)^-1 ; dot indicates tensor product
            ! Eq. 1 : eta = eps + 1/2 * eps*eps
            ! Also see Wallace "Thermodynamics of crystals"
            !
            ! % Matlab script
            ! % maximum difference is 1E-17 for e < 0.1, which is more than enough strain for elastic constants
            ! % if more strain is desired, then a differen method for solving for eps is required.
            ! for j = 1:200
            ! e = rand(6,1)*0.1;
            ! eta = [e(1),e(6)/2,e(5)/2;
            !        e(6)/2,e(2),e(4)/2;
            !        e(5)/2,e(4)/2,e(3)];
            ! eps = eta;
            ! for i = 1:300
            !     % eta = eps + 1/2eps*eps => rearrange:  eps = eta - 1/2eps*eps
            !     eps = eta - 1/2*eps*eps;
            ! end
            ! hold on; plot(j,max((eps+ 1/2*eps*eps-eta)),'s')
            ! end
            !
            eps_matrix = eta_matrix
            solve_for_eps : do j = 1,200
                eps_matrix = eta_matrix - 0.5_dp*matmul(eps_matrix,eps_matrix)
                if ( all(abs(eps_matrix).gt.10.0_dp) ) then
                    ! if this is the case then a differen method for solving for eps is required than iteration.
                    call am_print('ERROR','eps too large; possibly because it is diverging with each iteration: deformation may be too large.',flags='E')
                    call am_print('eps (true strain matrix)',eps_matrix)
                    stop
                endif
                if ( all(abs(eta-eps_matrix).lt.1E-15_dp) ) exit solve_for_eps
            enddo solve_for_eps
            !
            def_mat(1:3,1:3,i) = eye(3) + eps_matrix
            bas_def(1:3,1:3,i)=matmul(def_mat(1:3,1:3,i),bas)
            !
        enddo
        !
        if (verbosity.ge.1) then
            write(*,'(a5,a)') ' ... ', 'Deformation matrices (true strain, Voigt notation)'
            write(*,'(a5,2x,6a16)') ' ... ', 'e1', 'e2', 'e3', 'e4', 'e5', 'e6'
            do i = 1,nstrains
                write(*,'(5x,i2,6f16.10)') i,def_mat(1,1,i),def_mat(2,2,i),def_mat(3,3,i),def_mat(2,3,i),def_mat(1,3,i),def_mat(1,2,i)
            enddo
        endif
    end function    apply_elastic_deformation

end module am_unit_cell




