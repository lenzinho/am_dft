module am_vasp_io
    !
    use dispmodule
    use am_constants
    use am_stdout
    !
    implicit none
    !
    private
    !
    public :: read_ibzkpt
    public :: read_procar
    public :: read_prjcar
    public :: read_eigenval
    public :: read_poscar
    public :: read_amn
    public :: read_eig
    public :: write_poscar
    !
    contains

    !
    ! vasp
    !

    subroutine     read_ibzkpt(nkpts,kpt,w,ntets,vtet,tet,wtet,iopt_filename,iopt_verbosity) 
        !>
        !> Reads IBZKPT into variables.
        !>
        !> nkpts number of kpoints
        !> kpt(3,nkpts) kpoint vectors
        !> w(nkpts) weights (normalizd)
        !> ntets number of tetrahedra
        !> vtet(ntets) tetrahedron volume
        !> tet(4,ntets) tetrahedron connection table
        !> wtet(ntets) weights
        !>
        implicit none
        !
        character(len=*), optional, intent(in) :: iopt_filename
        character(max_argument_length) :: fname
        integer, optional, intent(in) :: iopt_verbosity
        integer :: verbosity
        integer ,              intent(out), optional :: nkpts   !> nkpts number of kpoints
        real(dp), allocatable, intent(out), optional :: kpt(:,:)!> kpt(3,nkpts) kpoint vectors
        real(dp), allocatable, intent(out), optional :: w(:)    !> weights normalized
        integer ,              intent(out), optional :: ntets   !> ntets number of tetrahedra
        real(dp), allocatable, intent(out), optional :: vtet(:) !> vtet(ntets) tetrahedron volume
        integer , allocatable, intent(out), optional :: tet(:,:)!> tet(4,ntets) tetrahedron connection table
        real(dp), allocatable, intent(out), optional :: wtet(:) !> wtet(ntets) weights
        integer  :: nkpts_internal, ntets_internal
        real(dp) :: vtet_tmp
        integer, allocatable :: w_int(:) ! w_int(nkpts) integer weights
        !
        integer :: iostat
        integer :: fid ! file id
        character(maximum_buffer_size) :: buffer ! read buffer
        character(len=:), allocatable :: word(:) ! read buffer
        integer :: i ! loop variables
        !
        fname = "IBZKPT"
        if ( present(iopt_filename) ) fname = iopt_filename
        verbosity = 1
        if ( present(iopt_verbosity) ) verbosity = iopt_verbosity
        !
        if (verbosity.ge.1) call print_title('Reading IBZKPT')
        !
        fid = 1
        open(unit=fid,file=trim(fname),status="old",action='read')
            !
            if (verbosity.ge.1) call am_print('actual file read is',trim(fname),' ... ')
            !
            ! (LINE 1) Automatically generated mesh
            read(fid,'(a)') buffer
            ! (LINE 2)      816
            read(unit=fid,fmt='(a)') buffer
            word = strsplit(buffer,delimiter=' ')
            !> nkpts number of kpoints
            read(word(1),*) nkpts_internal
            if (present(nkpts)) nkpts = nkpts_internal
            if (verbosity.ge.1) call am_print('number of kpoints',nkpts_internal,' ... ')
            !
            if (present(kpt)) allocate(kpt(3,nkpts_internal))
            allocate(w_int(nkpts_internal))
            ! (LINE 3) Reciprocal lattice
            read(fid,'(a)')
            ! (LINE 4) 0.00000000000000    0.00000000000000    0.00000000000000             1
            do i = 1, nkpts_internal
                read(unit=fid,fmt='(a)') buffer
                word = strsplit(buffer,delimiter=' ')
                !> kpt(3,nkpts) kpoint vectors
                if (present(kpt)) then
                    read(word(1),*) kpt(1,i)
                    read(word(2),*) kpt(2,i)
                    read(word(3),*) kpt(3,i)
                endif
                !> w_int(nkpts) weights
                read(word(4),*) w_int(i)
            enddo
            !> normalized weights
            if (present(w)) then
                allocate(w(nkpts_internal))
                w = w_int/real(sum(w_int),dp)
            endif
            ! (LINE 820) Tetrahedra
            if (present(vtet).or.present(tet).or.present(wtet)) then
                !
                read(unit=fid,fmt='(a)',iostat=iostat) buffer
                if ( iostat .eq. 0 ) then
                if ( buffer(1:1) .eq. "T" ) then
                    ! (LINE 821) 3976    0.00000559453079
                    read(unit=fid,fmt='(a)') buffer
                    word = strsplit(buffer,delimiter=' ')
                    !> ntets number of tetrahedra
                    read(word(1),*) ntets_internal
                    if (present(ntets)) ntets = ntets_internal
                    if (verbosity.ge.1) call am_print('number of tetrahedra',ntets_internal,' ... ')
                    !> vtet(ntets) tetrahedron volume
                    read(word(2),*) vtet_tmp
                    if (present(vtet)) then
                        allocate(vtet(ntets_internal))
                        vtet = vtet_tmp
                    endif
                    !
                    if (verbosity.ge.1) call am_print('tetrahedron volume (equal for all tetrahedra)',vtet_tmp,' ... ')
                    !
                    if (present(tet))  allocate(tet(4,ntets_internal))
                    if (present(wtet)) allocate(wtet(ntets_internal))
                    ! (LINE 822)  24         1         2         2        17
                    do i = 1, ntets_internal
                        read(unit=fid,fmt='(a)') buffer
                        word = strsplit(buffer,delimiter=' ')
                        !> wtet(ntets) weights
                        if (present(wtet)) then
                            read(word(1),*) wtet(i)
                        endif
                        !> tet(4,ntets) tetrahedron connection table
                        if (present(tet)) then
                            read(word(2),*) tet(1,i)
                            read(word(3),*) tet(2,i)
                            read(word(4),*) tet(3,i)
                            read(word(5),*) tet(4,i)
                        endif
                    enddo
                    if (verbosity.ge.1) then
                        if (present(tet)) then
                            if (i.gt.10) call am_print('irreducible k-point indicies of first ten tetrahedra',transpose(tet(:,1:10)),' ... ')
                            if (i.le.10) call am_print('irreducible k-point indicies tetrahedra',transpose(tet),' ... ')
                        endif
                    endif
                endif
                endif
                !
            endif
        close(fid)
    end subroutine read_ibzkpt 

    subroutine     read_procar(nkpts,nbands,nions,norbitals,nspins,E,occ,kpt,w,orbitals,lmproj,iopt_filename,iopt_verbosity)
        !>
        !> Reads PROCAR into file.
        !>
        !> nkpts number of kpoints
        !> nbands number of bands
        !> nions number of ions
        !> norbitals number of orbitals
        !> nspins number of spins
        !> E(nbands,nkpts) energies
        !> occ(nbands,nkpts) occupancies
        !> kpt(3,nkpts) kpoint
        !> w(nkpts) weights (normalized)
        !> orbitals(norbitals) names of orbitals
        !> lmproj(nspins,norbitals,nions,nbands,nkpts)
        !>
        implicit none
        !
        character(*), optional :: iopt_filename
        character(max_argument_length) :: fname
        integer, optional :: iopt_verbosity
        integer :: verbosity
        integer              , intent(out), optional :: nkpts     !> nkpts number of kpoints
        integer              , intent(out), optional :: nbands    !> nbands number of bands
        integer              , intent(out), optional :: nions     !> nions number of ions
        integer              , intent(out), optional :: norbitals !> norbitals number of orbitals
        integer              , intent(out), optional :: nspins    !> nspins number of spins
        real(dp), allocatable, intent(out), optional :: E(:,:)    !> E(nbands,nkpts) energies
        real(dp), allocatable, intent(out), optional :: occ(:,:)  !> occ(nbands,nkpts) occupancies
        real(dp), allocatable, intent(out), optional :: kpt(:,:)  !> kpt(3,nkpts) kpoint
        real(dp), allocatable, intent(out), optional :: w(:)      !> w(nkpts) weights (normalized)
        character(:), allocatable, intent(out), optional :: orbitals(:)     !> orbitals(norbitals) names of orbitals
        real(dp), allocatable, intent(out), optional :: lmproj(:,:,:,:,:)   !> lmproj(nspins,norbitals,nions,nbands,nkpts)
        ! <INTERNAL>
        integer :: internal_nkpts     !> nkpts number of kpoints
        integer :: internal_nbands    !> nbands number of bands
        integer :: internal_nions     !> nions number of ions
        integer :: internal_norbitals !> norbitals number of orbitals
        integer :: internal_nspins    !> nspins number of spins
        integer :: fid ! file id
        character(maximum_buffer_size) :: buffer ! read buffer
        character(len=:), allocatable :: word(:) ! read buffer
        integer :: i, j, l, m ! loop variables
        !
        fname = "PROCAR"
        if ( present(iopt_filename) ) fname = iopt_filename
        verbosity = 1
        if ( present(iopt_verbosity) ) verbosity = iopt_verbosity
        !
        internal_nspins = 1 ! not yet implemented!
        if (present(nspins)) nspins = internal_nspins
        !
        !
        !
        if (verbosity.ge.1) call print_title('Reading PROCAR')
        !
        fid = 1
        open(unit=fid,file=trim(fname),status="old",action='read')
            !
            if (verbosity.ge.1) call am_print('actual file read is',trim(fname),' ... ')
            ! (LINE 1) PROCAR lm decomposed
            read(fid,'(a)')
            ! (LINE 2) # of kpt-points:  160         # of bands:   8         # of ions:   2
            read(unit=fid,fmt='(a)') buffer
            word = strsplit(buffer,delimiter=' ')
            !> nkpts number of kpoints
            read(word(4),*)  internal_nkpts
            if (present(nkpts)) nkpts = internal_nkpts
            !> nbands number of bands
            read(word(8),*)  internal_nbands
            if (present(nbands)) nbands = internal_nbands
            !> nions number of ions
            read(word(12),*) internal_nions
            if (present(nions)) nions = internal_nions
            !
            if (verbosity.ge.1) call am_print('number of kpoints',internal_nkpts,' ... ')
            if (verbosity.ge.1) call am_print('number of bands',internal_nbands,' ... ')
            if (verbosity.ge.1) call am_print('number of ions',internal_nions,' ... ')
            !
            !> E(nbands,nkpts) energies
            !> occ(nbands,nkpts) occupancies
            !> kpt(3,nkpts) kpoint
            !> w(nkpts) normalized weights
            if (present(E))   allocate(E(internal_nbands,internal_nkpts))
            if (present(occ)) allocate(occ(internal_nbands,internal_nkpts))
            if (present(kpt)) allocate(kpt(3,internal_nkpts))
            if (present(w))   allocate(w(internal_nkpts))
            !
            do i = 1,internal_nkpts
                ! (LINE 3,62) skip line
                read(unit=fid,fmt=*)
                ! (LINE 4,63) kpt-point    1 :    0.50000000 0.50000000 0.50000000     weight = 0.00625000
                read(unit=fid,fmt='(a)') buffer
                word = strsplit(buffer,delimiter=' ')
                !> kpt(3,nkpts) kpoint
                !> w(nkpts) normalized weights
                if (present(kpt)) then
                    read(word(4),*) kpt(1,i)
                    read(word(5),*) kpt(2,i)
                    read(word(6),*) kpt(3,i)
                endif
                if (present(w)) then
                    read(word(9),*) w(i)
                endif
                !
                do j = 1, internal_nbands
                    ! (LINE 5,64) skip line
                    read(unit=fid,fmt=*)
                    ! (LINE 6,65) band   1 # energy   -3.91737559 # occ.  2.00000000
                    read(unit=fid,fmt='(a)') buffer
                    word = strsplit(buffer,delimiter=' ')
                    !> E(nbands,nkpts) energies
                    !> occ(nbands,nkpts) occupancies
                    if (present(E))   read(word(5),*) E(j,i)
                    if (present(occ)) read(word(8),*) occ(j,i)
                    ! (LINE 7,66) skip line
                    read(unit=fid,fmt=*)
                    ! (LINE 8,67) ion      s     py     pz     px    dxy    dyz    dz2    dxz    dx2    tot
                    read(unit=fid,fmt='(a)') buffer
                    if ( (i .eq. 1 ) .and. (j .eq. 1) ) then
                        word = strsplit(buffer,delimiter=' ')
                        !> norbitals number of orbitals
                        norbitals = size(word)-2 ! minus 2 to account for the ion, which is not an orbital, and the total!
                        if (present(norbitals)) norbitals = internal_norbitals
                        if (verbosity.ge.1) call am_print('number of orbitals',internal_norbitals,' ... ')
                        !> orbitals(norbitals) names of orbitals
                        if (present(orbitals)) then
                            allocate(character(len=15) :: orbitals(internal_norbitals))
                            do l = 1, internal_norbitals
                                orbitals(l) = word(l+1)
                            enddo
                            if (verbosity.ge.1) write(*,'(a5,a,100(x,a))' ) ' ... ', 'orbitals =', ( trim(orbitals(l)), l = 1, norbitals)
                        endif
                        !> lmproj(nspins,norbitals,nions,nbands,nkpts)
                        if (present(lmproj)) allocate(lmproj(internal_nspins,internal_norbitals,internal_nions,internal_nbands,internal_nkpts))
                        !
                    endif
                    ! (LINE 9,10) 
                    ! 1  0.224  0.003  0.001  0.009  0.000  0.000  0.000  0.000  0.000  0.237
                    ! 2  0.224  0.003  0.001  0.009  0.000  0.000  0.000  0.000  0.000  0.237
                    do l = 1, internal_nions
                        read(unit=fid,fmt='(a)') buffer
                        word = strsplit(buffer,delimiter=' ')
                        if (present(lmproj)) then
                            do m = 1, norbitals
                                !> lmproj(nspins,norbitals,nions,nbands,nkpts) (spins not implemented yet!)
                                read(word(m+1),*) lmproj(1,m,l,j,i)
                            enddo
                        endif
                    enddo
                    ! (LINE 11 = 60) tot  0.447  0.006  0.003  0.018  0.000  0.000  0.000  0.000  0.000  0.474
                    read(unit=fid,fmt=*)
                    !
                enddo
                ! (LINE 61) skip line
                read(unit=fid,fmt=*)
                !
            enddo
        close(fid)
    end subroutine read_procar

    subroutine     read_prjcar(nkpts,nbands,nspins,nkpts_prim,bas_prim,k_prim,w_prim,kpt,w,E,kproj,iopt_filename,iopt_verbosity)
        !>
        !> Reads PRJCAR into variables.
        !>
        !> nkpts_prim number of primitive kpoints (POSCAR.prim)
        !> nkpts number of kpoints (KPOINTS)
        !> nbands number of bands
        !> nspins number of spin components
        !> bas_prim(3,3) reciprocal basis (POSCAR.prim)
        !> k_prim(nkpts_prim) normalized kpoint weights (POSCAR.prim)
        !> kpt(3,nkpts) kvector (POSCAR)
        !> w(nkpts) normalized kpoint weights (POSCAR)
        !> E(nbands,nkpts) energies
        !> kproj(nspin,nkpts_prim,nbands,nkpts)
        !>
        implicit none
        !
        character(*), optional :: iopt_filename
        character(max_argument_length) :: fname
        integer, optional :: iopt_verbosity
        integer :: verbosity
        integer,intent(out) :: nkpts_prim      !> nkpts_prim number of primitive kpoints (POSCAR.prim)
        integer,intent(out) :: nkpts           !> nkpts number of kpoints (KPOINTS)
        integer,intent(out) :: nbands          !> nbands number of bands
        integer,intent(out) :: nspins          !> nspins number of spin components
        real(dp), intent(out) :: bas_prim(3,3) !> bas_prim(3,3) reciprocal basis (POSCAR.prim)
        real(dp), allocatable, intent(out) :: k_prim(:,:) !> k_prim(3,nkpts_prim) kvector (POSCAR.prim)
        real(dp), allocatable, intent(out) :: w_prim(:)   !> w_prim(nkpts_prim) weights (POSCAR.prim, normalized)
        real(dp), allocatable, intent(out) :: kpt(:,:)    !> kpt(3,nkpts) kvector (POSCAR)
        real(dp), allocatable, intent(out) :: w(:)        !> w(nkpts) weights (POSCAR, normalized)
        real(dp), allocatable, intent(out) :: E(:,:)      !> E(nbands,nkpts) energies
        real(dp), allocatable, intent(out) :: kproj(:,:,:,:) !> kproj(nspin,nkpts_prim,nbands,nkpts)
        integer :: n !> spin index
        ! i/o
        integer :: fid, iostat
        character(maximum_buffer_size) :: buffer
        character(len=:), allocatable :: word(:)
        ! loop variables
        integer :: i, j
        !        !
        fname = "PRJCAR"; if (present(iopt_filename)) fname = iopt_filename
        verbosity = 1; if (present(iopt_verbosity)) verbosity = iopt_verbosity
        !
        !
        !
        if (verbosity.ge.1) call print_title('Reading PRJCAR')
        !
        ! SKIM FILE TO DETEMRINE PARAMETERS NECESSARY TO ALLOCATE SPACE
        !
        fid = 1
        open(unit=fid,file=trim(fname),status="old",action='read')
            !
            if ( verbosity .ge. 1) call am_print('actual file read',trim(fname),' ... ')
            !
            nkpts_prim = 0
            nkpts  = 0
            nbands = 0
            nspins = 0
            !
            do
                read(unit=fid,fmt='(a)',iostat=iostat) buffer
                word = strsplit(buffer,delimiter=' ')
                !
                if ( size(word) .ge. 1 ) then
                select case (word(1))
                case('number')
                    read(word(8),*)  i
                    !> nkpts_prim number of primitive kpoints (POSCAR.prim)
                    if (i .gt. nkpts_prim) nkpts_prim = i
                case('spin')
                    read(word(3),*)  i
                    !> nspins number of spin components
                    if (i .gt. nspins) nspins = i
                case('kpt-point')
                    read(word(5),*)  i
                    !> nkpts number of kpoints (KPOINTS)
                    if (i .gt. nkpts) nkpts = i
                case('band:')
                    read(word(2),*) i
                    !> nbands number of bands
                    if (i .gt. nbands) nbands = i
                end select
            endif
            if (iostat .lt. 0) exit
            enddo
            !
        close(fid)
        !
        if (verbosity.ge.1) call am_print('number of kpoints to project onto',nkpts_prim,' ... ')
        if (verbosity.ge.1) call am_print('number of spin components',nspins,' ... ')
        if (verbosity.ge.1) call am_print('number of kpoints',nkpts,' ... ')
        if (verbosity.ge.1) call am_print('number of bands',nbands,' ... ')
        !> k_prim(3,nkpts_prim) kvector (POSCAR.prim)
        allocate(k_prim(3,nkpts_prim))
        !> w_prim(nkpts_prim) normalized kpoint weights (POSCAR.prim)
        allocate(w_prim(nkpts_prim))
        !> kpt(3,nkpts) kvector (POSCAR)
        allocate(kpt(3,nkpts))
        !> w(nkpts) normalized kpoint weights (POSCAR)
        allocate(w(nkpts))
        !> E(nbands,nkpts) energies
        allocate(E(nbands,nkpts))
        !> kproj(nspins,nkpts_prim,nbands,nkpts)
        allocate(kproj(nspins,nkpts_prim,nbands,nkpts))   
        !
        ! OPEN TO READ STUFF
        !
        fid = 1
        open(unit=fid,file=trim(fname),status="old",action='read')
            !
            ! (LINE 1) Basis vectors reciprocal space of POSCAR.prim (units of 2pi):
            read(fid,*)
            ! (LINE 2,3,4)     -0.2475248     0.2475248     0.2475248
            do i = 1,3
                read(unit=fid,fmt='(a)') buffer
                word = strsplit(buffer,delimiter=' ')
                !> bas_prim(3,3) reciprocal basis (POSCAR.prim)
                read(word(1),*) bas_prim(i,1)
                read(word(2),*) bas_prim(i,2)
                read(word(3),*) bas_prim(i,3)
            enddo
            ! (LINE 5)
            read(fid,*)
            ! (LINE 6) number of kpt-points in IBZ of POSCAR.prim:   120
                read(unit=fid,fmt='(a)') buffer
                word = strsplit(buffer,delimiter=' ')
                read(word(8),*) nkpts_prim
            ! (LINE 7)
            read(fid,*)
            ! (LINE 8)             b1            b2            b3      weight
            read(fid,*)
            ! (LINE 9)    1     0.0000000     0.0000000     0.0000000      1
            do i = 1,nkpts_prim
                read(unit=fid,fmt='(a)') buffer
                word = strsplit(buffer,delimiter=' ')
                !> k_prim(3,nkpts_prim) kvector (POSCAR.prim)
                read(word(2),*) k_prim(1,i)
                read(word(3),*) k_prim(2,i)
                read(word(4),*) k_prim(3,i)
                !> w_prim(nkpts_prim) normalized kpoint weights (POSCAR.prim)
                read(word(4),*) w_prim(i)
            enddo
            ! (LINE 129)
            read(fid,*)
            !
            do
                ! i = kpoint, j = bands, kpt = ions, m = orbital, n = spin
                read(unit=fid,fmt='(a)',iostat=iostat) buffer
                word = strsplit(buffer,delimiter=' ')
                !
                if ( size(word) .ge. 1 ) then
                select case (word(1))
                case('spin')
                    read(word(3),*)  n !> n is the index for spin
                case('kpt-point')
                    read(word(5),*)  i !> i is the index for kpoint
                    !> kpt(3,nkpts) kvector (POSCAR)
                    read(word(7),*)  kpt(1,i)
                    read(word(8),*)  kpt(2,i)
                    read(word(9),*)  kpt(3,i)
                    !> w(nkpts) normalized kpoint weights (POSCAR)
                    read(word(11),*) w(i)
                case('band')
                    read(word(2),*) j ! j is band index
                    read(word(4),*) E(j,i)
                    !> kproj(nspin,nkpts_prim,nbands,nkpts)
                    read(unit=fid,fmt=*) kproj(n,:,j,i)
                end select
            endif
            if (iostat .lt. 0) exit
            enddo
            !
        close(fid)
        !
    end subroutine read_prjcar

    subroutine     read_eigenval(nkpts,nbands,nspins,nelecs,kpt,w,E,lmproj,iopt_filename,iopt_verbosity)
        !>
        !> Read EIGENVAL into variables
        !>
        !> nkpts number of kpoints
        !> nbands number of bands
        !> nspins number of spins
        !> nelects number of electrons
        !> kpt(3,nkpts) kpoint
        !> w(nkpts) weights (normalized)
        !> E(nbands,nkpts) energies
        !> lmproj(nspins,norbitals,nions,nbands,nkpts) projections for spins 
        !>
        implicit none
        !
        integer :: nbands_without_spin
        integer :: ispin !> vasp ispin flag
        integer              , intent(out), optional :: nkpts      !> nkpts number of kpoints
        integer              , intent(out), optional :: nbands     !> nbands number of bands
        integer              , intent(out), optional :: nspins     !> nspins number of spins
        integer              , intent(out), optional :: nelecs     !> nelects number of electrons
        real(dp), allocatable, intent(out), optional :: kpt(:,:)   !> kpt(3,nkpts) kpoint
        real(dp), allocatable, intent(out), optional :: w(:)       !> w(nkpts) weights (normalized)
        real(dp), allocatable, intent(out), optional :: E(:,:)     !> E(nbands,nkpts) energies
        real(dp), allocatable, intent(out), optional :: lmproj(:,:,:,:,:) !> lmproj(nspins,norbitals,nions,nbands,nkpts) projections for spins
        ! <INTERNAL PARAMETERS>
        integer :: internal_nkpts     !> nkpts number of kpoints
        integer :: internal_nbands    !> nbands number of bands
        integer :: internal_nspins    !> nspins number of spins
        integer :: internal_norbitals !> norbitals number of ions - NO INFORMATION ABOUT THIS IN EIGEVAL, SET TO 1
        integer :: internal_nions     !> nions number of ions - NO INFORMATION ABOUT THIS IN EIGEVAL, SET TO 1
        ! i/o
        integer :: fid
        character(maximum_buffer_size) :: buffer
        character(len=:), allocatable :: word(:)
        ! loop variables
        integer :: i, j, m, j_and_m
        ! standard optionals
        character(*), intent(in), optional :: iopt_filename
        integer, intent(in), optional :: iopt_verbosity
        character(max_argument_length) :: fname
        integer :: verbosity
        !
        fname = "EIGENVAL"; if ( present(iopt_filename) ) fname = iopt_filename
        verbosity = 1; if ( present(iopt_verbosity) ) verbosity = iopt_verbosity
        !
        !
        !
        if (verbosity.ge.1) call print_title('Reading EIGENVAL')
        !
        fid = 1
        open(unit=fid,file=trim(fname),status="old",action='read')
            !
            if (verbosity.ge.1) call am_print('actual file read',trim(fname),' ... ')
            ! (LINE 1)     2    2    1    1
            read(unit=fid,fmt='(a)') buffer
            word = strsplit(buffer,delimiter=' ')
            !> vasp ispin flag
            read(word(4),*) ispin
            ! (LINE 2) skip
            read(fid,*)
            ! (LINE 3) skip
            read(fid,*)
            ! (LINE 4) skip
            read(fid,*)
            ! (LINE 5) skip
            read(fid,*)
            ! (LINE 6) skip
            read(unit=fid,fmt='(a)') buffer
            word = strsplit(buffer,delimiter=' ')
            !> nelects number of electrons
            write(*,'(a5,a,a)') ' ... ', 'electrons = ', trim(word(1))
            if (present(nelecs)) then
                read(word(1),*) nelecs
            endif
            !> nkpts number of kpoints
            read(word(2),*) internal_nkpts
            write(*,'(a5,a,a)') ' ... ', 'kpoints = ', trim(word(2))
            if (present(nkpts)) then
                nkpts = internal_nkpts
            endif
            !> nbands number of bands
            read(word(3),*) nbands_without_spin
            write(*,'(a5,a,a)') ' ... ', 'bands (neglecting spin) = ', trim(word(3))
            !
            if (present(kpt)) allocate(kpt(3,internal_nkpts))
            if (present(w))   allocate(w(internal_nkpts))
            !
            do i = 1, internal_nkpts
                ! (LINE 7) skip
                read(fid,*)
                ! (LINE 8)   0.0000000E+00  0.0000000E+00  0.0000000E+00  0.3356718E-04
                read(unit=fid,fmt='(a)') buffer
                word = strsplit(buffer,delimiter=' ')
                !> kpt(3,nkpts) kpoint
                if (present(kpt)) then
                    read(word(1),*) kpt(1,i)
                    read(word(2),*) kpt(2,i)
                    read(word(3),*) kpt(3,i)
                endif
                !> w(nkpts) weights
                if (present(w)) then
                    read(word(4),*) w(i)
                endif
                ! 
                j_and_m = 0
                do j = 1, nbands_without_spin
                    ! (LINE 9)    1      -52.220889   1.000000
                    read(unit=fid,fmt='(a)') buffer
                    word = strsplit(buffer,delimiter=' ')
                    if ( (i .eq. 1) .and. (j .eq. 1) ) then
                        !> nspins number of spins
                        internal_nspins = size(word)-1
                        ! <ORIGINAL>
                        ! internal_nspins = size(word)-2
                        ! </ORIGINAL>
                        write(*,'(a5,a,a)') ' ... ', 'spins = ', tostring(internal_nspins)
                        if (present(nspins)) then
                            nspins = internal_nspins
                        endif
                        internal_nbands = nbands_without_spin*internal_nspins
                        write(*,'(a5,a,a)') ' ... ', 'bands (with spin) = ', tostring(internal_nbands)
                        if (present(nbands)) then 
                            nbands = internal_nbands
                        endif
                        if (present(E)) allocate(E(internal_nbands,internal_nkpts))
                        !> lmproj(nspins,norbitals,nions,nbands,nkpts) projections for spins
                        !> this array contains the projection character onto everything from oribtals, to ions, to spins
                        !> for systems with 2 spins, the spin component of lmproj is binary with 0 or 1
                        !> for dirac spinors with 4 spin components (spin orbit coupling), the spin component of lmproj takes on non-integer values
                        internal_nions = 1
                        internal_norbitals = 1
                        if (present(lmproj)) then
                            allocate(lmproj(internal_nspins,internal_norbitals,internal_nions,internal_nbands,internal_nkpts))
                            lmproj = 0.0_dp
                        endif
                        !
                    endif
                    do m = 1, internal_nspins
                        j_and_m = j_and_m + 1
                        !> E(nbands,nkpts) energies
                        !> this array has one entry per band per spin
                        if (present(E)) then
                            read(word(m+1),*) E(j_and_m,i)
                        endif
                        !> lmproj(nspins,norbitals,nions,nbands,nkpts) projections for spins
                        if (present(lmproj)) lmproj(m,1,1,j_and_m,i) = 1.0_dp
                    enddo
                    !
                enddo
            enddo
        close(fid)
        !
    end subroutine read_eigenval

    subroutine     read_poscar(bas,natoms,nspecies,symb,tau_frac,atype,iopt_filename,iopt_verbosity)
        ! Reads poscar into variables performing basic checks on the way. Call it like this:
        !
        !    bas(3,3) column vectors a(1:3,i), a(1:3,j), a(1:3,kpt)
        !    natoms number of atoms
        !    nspecies number of unique atomic species
        !    symb(nspecies) list of atomic symbols
        !    natoms_per_species(natoms) number of atoms per species
        !    tau_frac(3,natoms) fractional atomic coordinates
        !    atype(natoms) type atom
        !
        !    call read_poscar( bas,natoms,nspecies,natoms_per_species,symb,tau_frac,atype )
        !
        implicit none
        ! subroutine i/o
        real(dp), intent(out) :: bas(3,3) ! column vectors a(1:3,i), a(1:3,j), a(1:3,kpt)
        integer,intent(out) :: natoms
        integer,intent(out) :: nspecies
        character(len=:), allocatable, intent(out) :: symb(:) ! 1 symb per species
        integer , allocatable :: natoms_per_species(:) ! number of atoms per species
        real(dp), allocatable, intent(out) :: tau_frac(:,:) ! atomic coordinates tau_frac_read(3,natoms) in fractional
        integer , allocatable, intent(out) :: atype(:) ! type of atom tau_frac_read(natoms)
        ! subroutine internal parameters
        real(dp) :: recbas(3,3)
        real(dp) :: vol
        real(dp), allocatable :: tau_frac_read(:,:) ! atomic coordinates tau_frac_read(3,natoms) read in
        logical :: direct
        logical :: cartesian
        logical :: selective_dynamics
        real(dp) :: scale
        ! file i/o
        integer :: fid
        character(maximum_buffer_size) :: buffer
        character(len=:), allocatable :: word(:)
        ! loop variables
        integer :: i, j, m, n
        ! optional i/o
        character(*), optional :: iopt_filename
        character(max_argument_length) :: fname
        integer, optional :: iopt_verbosity
        integer :: verbosity
        fname = "POSCAR"
        if ( present(iopt_filename) ) fname = iopt_filename
        verbosity = 1
        if ( present(iopt_verbosity) ) verbosity = iopt_verbosity
        !
        !
        !
        if (verbosity.ge.1) call print_title('Reading POSCAR')
        !
        fid = 1
        open(unit=fid,file=trim(fname),status="old",action='read')
            !
            if (verbosity.ge.1) call am_print('actual file read',trim(fname),' ... ')
            ! (LINE 1) header
             read(unit=fid,fmt='(a)')
            ! (LINE 2) lattice parameter scaling
            read(unit=fid,fmt=*) scale
            ! (LINES 3-5) basis
            do i = 1,3
                read(unit=fid,fmt='(a)') buffer
                word = strsplit(buffer,delimiter=' ')
                ! read primitive basis as column vectors
                ! note: vasp writes the basis vectors as row vectors
                do j = 1,3
                    read(word(j),*) bas(j,i)
                enddo
            enddo
            bas = bas*scale
            ! get cell volume
            vol=abs(bas(1,1)*(bas(2,2)*bas(3,3)-bas(3,2)*bas(2,3)) &
                   +bas(2,1)*(bas(3,2)*bas(1,3)-bas(1,2)*bas(3,3)) &
                   +bas(3,1)*(bas(1,2)*bas(2,3)-bas(1,3)*bas(2,2)))
            ! check that unit cell volume is non-zero
            if ( abs(vol) .lt. tiny ) then
                call am_print("ERROR","Basis vectors are not coplanar!"," >>> ")
                stop
            endif
            if (verbosity.ge.1) call am_print( "basis (column vectors)", bas, " ... " )
            ! get reciprocal basis as column vectors b
            recbas(1,1)=bas(2,2)*bas(3,3)-bas(3,2)*bas(2,3)
            recbas(1,2)=bas(3,2)*bas(1,3)-bas(1,2)*bas(3,3)
            recbas(1,3)=bas(1,2)*bas(2,3)-bas(1,3)*bas(2,2)
            recbas(2,1)=bas(2,3)*bas(3,1)-bas(2,1)*bas(3,3)
            recbas(2,2)=bas(1,1)*bas(3,3)-bas(3,1)*bas(1,3)
            recbas(2,3)=bas(2,1)*bas(1,3)-bas(1,1)*bas(2,3)
            recbas(3,1)=bas(2,1)*bas(3,2)-bas(2,2)*bas(3,1)
            recbas(3,2)=bas(3,1)*bas(1,2)-bas(1,1)*bas(3,2)
            recbas(3,3)=bas(1,1)*bas(2,2)-bas(1,2)*bas(2,1)
            recbas = recbas/vol
            ! (LINE 6) V  N
            read(unit=fid,fmt='(a)') buffer
            ! strslplit => allocate internally, character(500)::symb(j)
            symb = strsplit(buffer,delimiter=' ')
            nspecies = size(symb)
            if (verbosity.ge.1) call am_print( "number of unique atomic species", nspecies, " ... ")
            if (verbosity.ge.1) write(*,'(" ... ","atomic species =",100a3)') ( symb(i)(1:2), i = 1, nspecies)
            ! (LINE 6) 1  1
            allocate(natoms_per_species(nspecies))
            read(fid,*) natoms_per_species
            if (verbosity.ge.1) write(*,'(" ... ","number of atoms per species =",100i4)') ( natoms_per_species(i), i = 1, nspecies)
            ! get total number of atoms
            natoms=sum(natoms_per_species)
            ! write the total number of atoms
            if (verbosity.ge.1) call am_print("total number of atoms", natoms, " ... ")
            ! (LINE 8) Direct
            direct = .false.
            cartesian = .false.
            selective_dynamics = .false.
            !
            read(unit=fid,fmt='(a)') buffer
            select case (buffer(1:1))
            case('D','d')
                direct = .true.
            case('C','c','K','kpt')
                cartesian = .true.
            case('S','s')
                selective_dynamics = .true.
                !
                read(unit=fid,fmt='(a)') buffer
                select case (buffer(1:1))
                case('D','d')
                    direct = .true.
                case('C','c','K','kpt')
                    cartesian = .true.
                end select
                !
            end select
            !
            if ( selective_dynamics ) then
                if (verbosity.ge.1) write(*,'(a5,a)') ' ... ', 'selective dynamics'
            endif
            if ( direct ) then
                if (verbosity.ge.1) write(*,'(a5,a)') ' ... ', 'atomic coordinates read in with fractional units and converted to cartesian coordinates'
            endif
            if ( cartesian ) then
                if (verbosity.ge.1) write(*,'(a5,a)') ' ... ', 'atomic coordinates read in with cartesian coordinates'
            endif
            ! (LINE 9-10)   0.1250000000000000  0.1250000000000000  0.1250000000000000
            ! read atomic positions
            allocate(tau_frac_read(3,natoms))      ! coordinates of atoms read in
            allocate(tau_frac(3,natoms)) ! coordinates of atoms in cartesian  coordinates
            allocate(atype(natoms)) ! type of atom, integer 
            m = 0
            do i = 1, nspecies
                do j = 1, natoms_per_species(i)
                    !
                    m = m + 1
                    atype(m) = i
                    ! read atomic coordinates
                    read(unit=fid,fmt='(a)') buffer
                    word = strsplit(buffer,delimiter=' ')
                    do n = 1,3
                        read(word(n),*) tau_frac_read(n,m)
                    enddo
                    ! make sure everything is in fractional
                    if ( direct ) then
                        ! tau_frac(1:3,m) = matmul( bas,tau_frac_read(1:3,m) ) ! cartesian
                        tau_frac(1:3,m) = tau_frac_read(1:3,m)
                    endif
                    if ( cartesian ) then
                        ! tau_frac(1:3,m) = scale * tau_frac_read(1:3,m) ! cartesian
                        tau_frac(1:3,m) = matmul( recbas, scale * tau_frac_read(1:3,m) )
                    endif
                enddo
            enddo
            !
            if (verbosity.ge.1) then
                write(*,'(" ... ",a)') "atomic positions (fractional, cartesian)"
                m = 0
                do i = 1, nspecies
                    do j = 1, natoms_per_species(i)
                        !
                        m = m + 1
                        !
                        if (m.lt.11) then
                            write(*,'(5x,a5,a3,3f13.8,5x,3f13.8)') int2char(m), symb(atype(m)), tau_frac(1:3,m), matmul(bas,tau_frac(1:3,m))
                        elseif (m.eq.11) then
                            write(*,'(5x,a)') repeat(" ... ",15)
                        endif
                        !
                    enddo
                enddo
            endif
            !
        close(fid)
    end subroutine read_poscar

    subroutine     write_poscar(bas,natoms,nspecies,symb,tau_frac,atype,iopt_filename,iopt_header)
        !>
        !> Write poscar based on structure data.
        !>
        implicit none
        ! subroutine i/o
        real(dp), intent(in) :: bas(3,3) ! column vectors a(1:3,i), a(1:3,j), a(1:3,kpt)
        integer , intent(in) :: natoms
        integer , intent(in) :: nspecies
        character(len=*), intent(in) :: symb(:) ! 1 symb per species
        real(dp), intent(in) :: tau_frac(:,:) ! atomic coordinates tau_frac_read(3,natoms) in fractional
        integer , intent(in) :: atype(:) ! type of atom tau_frac_read(natoms)
        ! file i/o
        integer :: fid
        ! loop variables
        integer :: i, j, n
        ! optional i/o
        character(*), optional :: iopt_filename
        character(max_argument_length) :: fname
        character(*), optional :: iopt_header
        character(max_argument_length) :: header
        fname = 'outfile.POSCAR'
        if ( present(iopt_filename) ) fname = iopt_filename
        header = 'POSCAR'
        if ( present(iopt_header) )  header = iopt_header
        !
        !
        !
        call print_title('Writing POSCAR')
        !
        fid = 1
        open(unit=fid,file=trim(fname),status='replace',action='write')
            !
            call am_print('file',trim(fname),' ... ')
            ! (LINE 1) header
            write(unit=fid,fmt='(a)') trim(header)
            ! (LINE 2) lattice parameter scaling
            write(unit=fid,fmt='(f20.12)') 1.0_dp
            ! (LINES 3-5) basis
            do i = 1,3
            do j = 1,3
                write(unit=fid,fmt='(f20.12)',advance='no') bas(j,i)
            enddo
            write(unit=fid,fmt=*)
            enddo
            ! (LINE 6) V  N
            do i = 1, nspecies
            write(unit=fid,fmt='(x,a,x)',advance='no') trim(symb(i))
            enddo
            write(unit=fid,fmt=*)
            ! (LINE 7) 1  1
            do i = 1, nspecies
            write(unit=fid,fmt='(x,a,x)',advance='no') trim(int2char(count(atype.eq.i)))
            enddo
            write(unit=fid,fmt=*)
            ! (LINE 8) Direct
            write(unit=fid,fmt='(a)') 'Direct'
            ! (LINE 9-10)   0.1250000000000000  0.1250000000000000  0.1250000000000000
            do i = 1, nspecies
            do j = 1, natoms
            if (atype(j).eq.i) then
                !
                do n = 1,3
                    write(unit=fid,fmt='(f20.12)',advance='no') tau_frac(n,j)
                enddo
                write(unit=fid,fmt=*)
                !
            endif
            enddo
            enddo
            !
        close(fid)
    end subroutine write_poscar

    !
    ! wannier90
    !

    subroutine     read_amn(U,nbands,nkpts,nwanniers,iopt_filename,iopt_verbosity)
        !>
        !> U is equal to a_matrix if there is no disentanglment.
        !>
        !>
        implicit none
        !
        complex(dp), intent(out), allocatable :: U(:,:,:)
        integer, intent(out) :: nbands, nkpts, nwanniers
        integer  :: nwanprojs
        real(dp) :: a_real, a_imag
        integer  :: i, m, n, k
        !
        character(maximum_buffer_size) :: buffer ! read buffer
        character(len=:), allocatable :: word(:) ! read buffer
        !
        integer :: fid
        character(*), optional :: iopt_filename
        character(max_argument_length) :: fname
        integer, optional :: iopt_verbosity
        integer :: verbosity
        !
        fname = "wannier90.amn"
        if ( present(iopt_filename) ) fname = iopt_filename
        verbosity = 1
        if ( present(iopt_verbosity) ) verbosity = iopt_verbosity
        !
        !
        !
        if (verbosity.ge.1) call print_title('Reading wannier unitary matrix (U matrix)')
        !
        fid = 1
        open(unit=fid,file=trim(fname),status="old",action='read')
            !
            if (verbosity.ge.1) call am_print('actual file read is',trim(fname),' ... ')
            !
            ! (LINE 1) :File generated by VASP: unknown system
            read(fid,*)
            ! (LINE 2)            8          64           8
            read(unit=fid,fmt='(a)') buffer
            word = strsplit(buffer,delimiter=' ')
            read(word(1),*) nbands
            read(word(2),*) nkpts
            read(word(3),*) nwanniers
            nwanprojs=nbands*nwanniers*nkpts
            if (verbosity.ge.1) call am_print('number of bands',nbands,' ... ')
            if (verbosity.ge.1) call am_print('number of kpoints',nkpts,' ... ')
            if (verbosity.ge.1) call am_print('number of neighbors',nwanniers,' ... ')
            if (verbosity.ge.1) call am_print('number of projections',nwanprojs,' ... ')
            ! Read the projections
            allocate(U(nbands,nbands,nkpts))
            do i = 1, nwanprojs
                read(unit=fid,fmt='(a)') buffer
                word = strsplit(buffer,delimiter=' ')
                read(word(1),*) m ! loop over bands
                read(word(2),*) n ! loop over bands
                read(word(3),*) k ! loop over kpoints
                read(word(4),*) a_real
                read(word(4),*) a_imag
                U(m,n,k) = cmplx(a_real,a_imag,kind=dp)
            end do
            !
        close(fid)
        !
    end subroutine read_amn

    subroutine     read_eig(E,nbands,nkpts,iopt_filename,iopt_verbosity)
        !>
        !> it seems that the U is equal to the a_matrix if there is no disentanglment.
        !>
        !>
        implicit none
        !
        real(dp), intent(out), allocatable :: E(:,:) !> E(nbands,nkpts) energies
        integer, intent(out) :: nbands
        integer, intent(out) :: nkpts
        real(dp) :: Ein
        integer  :: i, j, k
        !
        character(maximum_buffer_size) :: buffer ! read buffer
        character(len=:), allocatable :: word(:) ! read buffer
        !
        integer :: fid
        integer :: iostat
        character(*), optional :: iopt_filename
        character(max_argument_length) :: fname
        integer, optional :: iopt_verbosity
        integer :: verbosity
        !
        fname = "wannier90.eig"
        if ( present(iopt_filename) ) fname = iopt_filename
        verbosity = 1
        if ( present(iopt_verbosity) ) verbosity = iopt_verbosity
        !
        if (verbosity.ge.1) call print_title('Reading eigenvalues for wannier')
        !
        ! SKIM FILE TO DETEMRINE PARAMETERS NECESSARY TO ALLOCATE SPACE
        !
        fid = 1
        open(unit=fid,file=trim(fname),status="old",action='read')
            !
            if (verbosity.ge.1) call am_print('actual file read is',trim(fname),' ... ')
            !
            nkpts  = 0
            nbands = 0
            !
            do
                ! (LINE 1)           1           1       -6.764130089716
                read(unit=fid,fmt='(a)',iostat=iostat) buffer
                word = strsplit(buffer,delimiter=' ')
                if ( iostat .eq. 0 ) then
                    !
                    if ( size(word) .ge. 1 ) then
                        !
                        read(word(1),*) i
                        if (i .gt. nbands) nbands = i
                        !
                        read(word(2),*) i
                        if (i .gt. nkpts) nkpts = i
                        !
                    endif
                    !
                else
                    exit
                endif
            enddo
        close(fid)
        !
        if (verbosity.ge.1) call am_print('number of kpoints',nkpts,' ... ')
        if (verbosity.ge.1) call am_print('number of bands',nbands,' ... ')
        !
        ! ACTUALLY READ IN THINGS NOW
        !
        allocate(E(nbands,nkpts)) !> E(nbands,nkpts)
        !
        fid = 1
        open(unit=fid,file=trim(fname),status="old",action='read')
            !
            do i = 1, (nbands*nkpts)
                ! (LINE 1)           1           1       -6.764130089716
                read(unit=fid,fmt='(a)') buffer
                word = strsplit(buffer,delimiter=' ')
                read(word(1),*) j
                read(word(2),*) k
                read(word(3),*) Ein
                E(j,k) = Ein
            enddo
        close(fid)
        !
    end subroutine read_eig

!     subroutine read_mmn(iopt_filename,iopt_verbosity)
!         !>
!         !> Definition of overlap m_matrix. Eq. (25) PhysRevB 56 12847
!         !>
!         !>    M_{m,n}^{(k,b)} = < u_{m,k} | u_{n,k+b} >
!         !> 
!         !> The overlap matrix is a rank 4 matrix which has indices running over bands m and n,
!         !> over kpoints k, and over neighboring kpoints to k. 
!         !>  
!         !> a_matrix      = projection of trial orbitals on bloch states
!         !> m_matrix_orig = overlap of bloch states
!         !> U      = the unitary rotations from the optimal subspace to the
!         !>                 optimally smooth states. 
!         !
!         implicit none
!         !
!         fname = "wannier90.mmn"
!         if ( present(iopt_filename) ) fname = iopt_filename
!         verbosity = 1
!         if ( present(iopt_verbosity) ) verbosity = iopt_verbosity
!         !
!         !
!         !
!         if (verbosity.ge.1) call print_title('Reading wannier overlaps (MMN)')
!         !
!         fid = 1
!         open(unit=fid,file=trim(fname),status="old",action='read')
!             !
!             if (verbosity.ge.1) call am_print('actual file read is',fname,' ... ')
!             !
!             ! (LINE 1) :File generated by VASP: unknown system
!             read(fid,*)
!             ! (LINE 2)            8          64           8
!             read(unit=fid,fmt='(a)') buffer
!             word = strsplit(buffer,delimiter=' ')
!             read(word(1),*) nbands
!             read(word(2),*) nkpts
!             read(word(3),*) nneighbors
!             noverlaps=nkpts*nneighbors
!             if (verbosity.ge.1) call am_print('number of bands',nbands,' ... ')
!             if (verbosity.ge.1) call am_print('number of kpoints',nkpts,' ... ')
!             if (verbosity.ge.1) call am_print('number of neighbors',nneighbors,' ... ')
!             if (verbosity.ge.1) call am_print('number of overlaps',noverlaps,' ... ')
!             !
!             allocate(mmn_tmp(nbands,nbands))
!             !
!             do i = 1, noverlaps
!                 ! (LINE 3)     1    2    0    0    0
!                 read(unit=fid,fmt='(a)') buffer
!                 word = strsplit(buffer,delimiter=' ')
!                 read(word(1),*) kpt1 ! kpoint 1
!                 read(word(2),*) kpt2 ! kpoint 2
!                 read(word(3),*) nnl  ! bz of nearest neighbor index l
!                 read(word(4),*) nnm  ! bz of nearest neighbor index m
!                 read(word(5),*) nnn  ! bz of nearest neighbor index n
!                 ! (LINE 3 - nbands^2+3) 0.208146341152   -0.972227642054
!                 do n=1,nbands
!                 do m=1,nbands
!                     read(unit=fid,fmt='(a)') buffer
!                     word = strsplit(buffer,delimiter=' ')
!                     read(word(1),*) m_real
!                     read(word(2),*) m_imag
!                     mmn_tmp(m,n) = cmplx(m_real,m_imag,kind=dp)
!                 enddo
!                 enddo

!                 ! nncell(3,nkpts,nneighbors)
!                 ! nncell(i,kpt,nn) tells us in which BZ is nn-th nearest-neighbour of kpoint kpt. 

!                 !
!                 nn=0
!                 nn_found=.false.
!                 do inn = 1, nneighbors
!                     if ((kpt2.eq.nnlist(kpt1,inn)).and. &
!                        (nnl.eq.nncell(1,kpt1,inn)).and. &
!                        (nnm.eq.nncell(2,kpt1,inn)).and. &
!                        (nnn.eq.nncell(3,kpt1,inn)) ) then
!                         if (.not.nn_found) then
!                             nn_found=.true.
!                             nn=inn
!                         else
!                             call io_error('Error reading '//trim(seedname)// &
!                                 '.mmn. More than one matching nearest neighbour found')
!                         endif
!                     endif
!                 end do
!                 if (nn.eq.0) then
!                     write(stdout,'(/a,i8,2i5,i4,2x,3i3)') &
!                     ' Error reading '//trim(seedname)//'.mmn:',i,kpt1,kpt2,nn,nnl,nnm,nnn
!                     call io_error('Neighbour not found')
!                 end if
!                 if (disentanglement) then
!                     m_matrix_orig(:,:,nn,kpt1) = mmn_tmp(:,:)
!                 else
!                     ! disentanglement=.false. means numbands=numwann, so no the dimensions are the same
!                     ! m_matrix(m,n,nn,kpt1) ! defined for each band 
!                     m_matrix(:,:,nn,kpt1) = mmn_tmp(:,:)
!                 end if
!             end do
!             !
!         close(fid)    

!     end subroutine read_mmn

end module
    