! #:include "fypp_macros.fpp"
module am_vasp_io

    use dispmodule
    use am_constants
    use am_matlab
    use am_mkl, only : det

    implicit none

    public

    contains

    ! vasp

    subroutine     read_ibzkpt(nkpts,kpt,w,ntets,vtet,tet,wtet,fname,verbosity) 
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
        integer ,              intent(out), optional :: nkpts   !> nkpts number of kpoints
        real(dp), allocatable, intent(out), optional :: kpt(:,:)!> kpt(3,nkpts) kpoint vectors
        real(dp), allocatable, intent(out), optional :: w(:)    !> weights normalized
        integer ,              intent(out), optional :: ntets   !> ntets number of tetrahedra
        real(dp), allocatable, intent(out), optional :: vtet(:) !> vtet(ntets) tetrahedron volume
        integer , allocatable, intent(out), optional :: tet(:,:)!> tet(4,ntets) tetrahedron connection table
        real(dp), allocatable, intent(out), optional :: wtet(:) !> wtet(ntets) weights
        character(len=*), intent(in) :: fname
        integer         , intent(in) :: verbosity
        integer  :: nkpts_internal, ntets_internal
        real(dp) :: vtet_tmp
        integer, allocatable :: w_int(:) ! w_int(nkpts) integer weights
        integer :: iostat
        integer :: fid ! file id
        character(maximum_buffer_size) :: buffer ! read buffer
        character(len=:), allocatable :: word(:) ! read buffer
        integer :: i ! loop variables
        !
        if (verbosity.ge.1) call print_title('IBZKPT')
        !
        fid = 1
        open(unit=fid,file=trim(fname),status='old',action='read',iostat=iostat)
            if (iostat.ne. 0) stop 'ERROR [read_ibzkpt]: problem reading ibzkpt. check if file exists'
            !
            if (verbosity.ge.1) write(*,'(a,a,a)') flare, 'file read = ', trim(fname)
            !
            ! (LINE 1) Automatically generated mesh
            read(fid,'(a)') buffer
            ! (LINE 2)      816
            read(unit=fid,fmt='(a)') buffer
            word = strsplit(buffer,delimiter=' ')
            !> nkpts number of kpoints
            read(word(1),*) nkpts_internal
            if (present(nkpts)) nkpts = nkpts_internal
            if (verbosity.ge.1) write(*,'(a,a,a)') flare, 'kpoints = ', tostring(nkpts_internal)
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
                    if (verbosity.ge.1) write(*,'(a,a,a)') flare, 'tetrahedra = ', tostring(ntets_internal)
                    !> vtet(ntets) tetrahedron volume
                    read(word(2),*) vtet_tmp
                    if (present(vtet)) then
                        allocate(vtet(ntets_internal))
                        vtet = vtet_tmp
                    endif
                    !
                    if (verbosity.ge.1) write(*,'(a,a,a)') flare, 'tetrahedra volume = ', tostring(vtet_tmp)
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
                            if (i.gt.10) then 
                                write(*,'(a,a)') flare, 'irreducible tetrahedra indicies (first ten)'
                                call disp_indent()
                                call disp(transpose(tet(:,1:10)))
                            endif
                            if (i.le.10) then 
                                write(*,'(a,a)') flare, 'irreducible tetrahedra indicies'
                                call disp_indent()
                                call disp(transpose(tet))
                            endif
                        endif
                    endif
                endif
                endif
                !
            endif
        close(fid)
    end subroutine read_ibzkpt 

    subroutine     read_procar(nkpts,nbands,nions,norbitals,nspins,E,occ,kpt,w,orbitals,lmproj,fname,verbosity)
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
        character(len=*), intent(in) :: fname
        integer         , intent(in) :: verbosity
        ! <INTERNAL>
        integer :: x_nkpts     !> nkpts number of kpoints
        integer :: x_nbands    !> nbands number of bands
        integer :: x_nions     !> nions number of ions
        integer :: x_norbitals !> norbitals number of orbitals
        integer :: x_nspins    !> nspins number of spins
        integer :: fid ! file id
        integer :: iostat
        character(maximum_buffer_size) :: buffer ! read buffer
        character(len=:), allocatable :: word(:) ! read buffer
        integer :: i, j, l, m ! loop variables
        !
        x_nspins = 1 ! not yet implemented!
        if (present(nspins)) nspins = x_nspins
        !
        if (verbosity.ge.1) call print_title('PROCAR')
        !
        fid = 1
        open(unit=fid,file=trim(fname),status='old',action='read',iostat=iostat)
            if (iostat.ne. 0) stop 'ERROR [read_procar]: problem reading procar. check if file exists'
            !
            if (verbosity.ge.1) write(*,'(a,a,a)') flare, 'tetrahedra = ', trim(fname)
            ! (LINE 1) PROCAR lm decomposed
            read(fid,'(a)')
            ! (LINE 2) # of kpt-points:  160         # of bands:   8         # of ions:   2
            read(unit=fid,fmt='(a)') buffer
            word = strsplit(buffer,delimiter=' ')
            !> nkpts number of kpoints
            read(word(4),*)  x_nkpts
            if (present(nkpts)) nkpts = x_nkpts
            !> nbands number of bands
            read(word(8),*)  x_nbands
            if (present(nbands)) nbands = x_nbands
            !> nions number of ions
            read(word(12),*) x_nions
            if (present(nions)) nions = x_nions
            !
            if (verbosity.ge.1) write(*,'(a,a,a)') flare, 'kpoints = ', tostring(x_nkpts)
            if (verbosity.ge.1) write(*,'(a,a,a)') flare, 'bands = ', tostring(x_nbands)
            if (verbosity.ge.1) write(*,'(a,a,a)') flare, 'ions = ', tostring(x_nions)
            !
            !> E(nbands,nkpts) energies
            !> occ(nbands,nkpts) occupancies
            !> kpt(3,nkpts) kpoint
            !> w(nkpts) normalized weights
            if (present(E))   allocate(E(x_nbands,x_nkpts))
            if (present(occ)) allocate(occ(x_nbands,x_nkpts))
            if (present(kpt)) allocate(kpt(3,x_nkpts))
            if (present(w))   allocate(w(x_nkpts))
            !
            do i = 1,x_nkpts
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
                do j = 1, x_nbands
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
                        if (present(norbitals)) norbitals = x_norbitals
                        if (verbosity.ge.1) write(*,'(a,a,a)') flare, 'orbitals = ', tostring(x_norbitals)
                        !> orbitals(norbitals) names of orbitals
                        if (present(orbitals)) then
                            allocate(character(len=15) :: orbitals(x_norbitals))
                            do l = 1, x_norbitals
                                orbitals(l) = word(l+1)
                            enddo
                            if (verbosity.ge.1) write(*,'(a,a,100(x,a))' ) flare, 'orbitals =', ( trim(orbitals(l)), l = 1, norbitals)
                        endif
                        !> lmproj(nspins,norbitals,nions,nbands,nkpts)
                        if (present(lmproj)) allocate(lmproj(x_nspins,x_norbitals,x_nions,x_nbands,x_nkpts))
                        !
                    endif
                    ! (LINE 9,10) 
                    ! 1  0.224  0.003  0.001  0.009  0.000  0.000  0.000  0.000  0.000  0.237
                    ! 2  0.224  0.003  0.001  0.009  0.000  0.000  0.000  0.000  0.000  0.237
                    do l = 1, x_nions
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

    subroutine     read_prjcar(nkpts,nbands,nspins,nkpts_prim,bas_prim,k_prim,w_prim,kpt,w,E,kproj,fname,verbosity)
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
        !> kproj(nspins,nkpts_prim,nbands,nkpts)
        !>
        implicit none
        !
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
        real(dp), allocatable, intent(out) :: kproj(:,:,:,:) !> kproj(nspins,nkpts_prim,nbands,nkpts)
        character(len=*), intent(in) :: fname
        integer         , intent(in) :: verbosity
        integer :: n !> spin index
        ! i/o
        integer :: fid, iostat
        character(maximum_buffer_size) :: buffer
        character(len=:), allocatable :: word(:)
        ! loop variables
        integer :: i, j
        !
        if (verbosity.ge.1) call print_title('PRJCAR')
        !
        ! SKIM FILE TO DETEMRINE PARAMETERS NECESSARY TO ALLOCATE SPACE
        !
        fid = 1
        open(unit=fid,file=trim(fname),status='old',action='read',iostat=iostat)
            if (iostat.ne. 0) stop 'ERROR [read_prjcar]: problem reading prjcar. check if file exists'
            !
            if ( verbosity .ge. 1) write(*,'(a,a,a)') flare, 'file = ', trim(fname)
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
        if (verbosity.ge.1) write(*,'(a,a,a)') flare, 'kpoint projections = ', tostring(nkpts_prim)
        if (verbosity.ge.1) write(*,'(a,a,a)') flare, 'kpoints = ', tostring(nkpts)
        if (verbosity.ge.1) write(*,'(a,a,a)') flare, 'spins = ', tostring(nspins)
        if (verbosity.ge.1) write(*,'(a,a,a)') flare, 'bands = ', tostring(nbands)
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
                    !> kproj(nspins,nkpts_prim,nbands,nkpts)
                    read(unit=fid,fmt=*) kproj(n,:,j,i)
                end select
            endif
            if (iostat .lt. 0) exit
            enddo
            !
        close(fid)
        !
    end subroutine read_prjcar

    subroutine     read_eigenval(nkpts,nbands,nspins,nelecs,kpt,w,E,lmproj,fname,verbosity)
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
        character(len=*)     , intent(in) :: fname
        integer              , intent(in) :: verbosity
        ! <INTERNAL PARAMETERS>
        integer :: x_nkpts     !> nkpts number of kpoints
        integer :: x_nbands    !> nbands number of bands
        integer :: x_nspins    !> nspins number of spins
        integer :: x_norbitals !> norbitals number of ions - NO INFORMATION ABOUT THIS IN EIGEVAL, SET TO 1
        integer :: x_nions     !> nions number of ions - NO INFORMATION ABOUT THIS IN EIGEVAL, SET TO 1
        ! i/o
        real(dp) :: try
        integer :: fid
        integer :: iostat
        character(maximum_buffer_size) :: buffer
        character(len=:), allocatable :: word(:)
        ! loop variables
        integer :: i, j, m, j_and_m
        !
        if (verbosity.ge.1) call print_title('EIGENVAL')
        !
        fid = 1
        open(unit=fid,file=trim(fname),status='old',action='read',iostat=iostat)
            if (iostat.ne. 0) stop 'ERROR [read_eigenval]: problem reading eigenval. check if file exists'
            !
            if (verbosity.ge.1) write(*,'(a,a,a)') flare, 'file = ', trim(fname)
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
            if (verbosity.ge.1) write(*,'(a,a,a)') flare, 'electrons = ', trim(word(1))
            if (present(nelecs)) then
                read(word(1),*) nelecs
            endif
            !> nkpts number of kpoints
            read(word(2),*) x_nkpts
            if (verbosity.ge.1) write(*,'(a,a,a)') flare, 'kpoints = ', trim(word(2))
            if (present(nkpts)) then
                nkpts = x_nkpts
            endif
            !> nbands number of bands
            read(word(3),*) nbands_without_spin
            if (verbosity.ge.1) write(*,'(a,a,a)') flare, 'bands (neglecting spin) = ', trim(word(3))
            !
            if (present(kpt)) allocate(kpt(3,x_nkpts))
            if (present(w))   allocate(w(x_nkpts))
            !
            do i = 1, x_nkpts
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
                    ! (LINE 9)    1      -46.290491  -- sometimes 1.00000 is missing.
                    read(unit=fid,fmt='(a)') buffer
                    word = strsplit(buffer,delimiter=' ')
                    if ( (i .eq. 1) .and. (j .eq. 1) ) then
                        ! determine which format
                        read(word(size(word)),*) try
                        if (abs(try-1.0_dp).lt.tiny) then
                            !> nspins number of spins
                            x_nspins = size(word)-2
                        else
                            !> nspins number of spins
                            x_nspins = size(word)-1
                        endif
                        if (verbosity.ge.1) write(*,'(a,a,a)') flare, 'spins = ', tostring(x_nspins)
                        if (present(nspins)) then
                            nspins = x_nspins
                        endif
                        x_nbands = nbands_without_spin*x_nspins
                        if (verbosity.ge.1) write(*,'(a,a,a)') flare, 'bands (with spin) = ', tostring(x_nbands)
                        if (present(nbands)) then 
                            nbands = x_nbands
                        endif
                        if (present(E)) allocate(E(x_nbands,x_nkpts))
                        !> lmproj(nspins,norbitals,nions,nbands,nkpts) projections for spins
                        !> this array contains the projection character onto everything from oribtals, to ions, to spins
                        !> for systems with 2 spins, the spin component of lmproj is binary with 0 or 1
                        !> for dirac spinors with 4 spin components (spin orbit coupling), the spin component of lmproj takes on non-integer values
                        x_nions = 1
                        x_norbitals = 1
                        if (present(lmproj)) then
                            allocate(lmproj(x_nspins,x_norbitals,x_nions,x_nbands,x_nkpts))
                            lmproj = 0.0_dp
                        endif
                        !
                    endif
                    do m = 1, x_nspins
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

    subroutine     query_wavecar(n_fft,nkpts,nbands,fname,verbosity)
        !
        implicit none
        !
        integer     , intent(out), optional :: n_fft(3)
        integer     , intent(out), optional :: nbands
        integer     , intent(out), optional :: nkpts
        character(*), intent(in) , optional :: fname
        integer     , intent(in) , optional :: verbosity
        integer :: kpt_id  ! kpt,band,spin index of wavefunction to read
        integer :: band_id ! kpt,band,spin index of wavefunction to read
        integer :: spin_id ! kpt,band,spin index of wavefunction to read
        real(dp) :: bas(3,3)
        real(dp) :: recbas(3,3) ! two pi convetion = keep it local to this function.
        real(dp) :: ecut
        real(dp) :: vol
        real(dp) :: xrecl,xnspins,xprec,xnkpts,xnbands
        integer  ::  recl, nspins, prec
        integer :: i, j ! loop variables
        integer :: fid
        integer :: iostat
        !
        ! initialize parameters
        fid = 1
        kpt_id = 1
        band_id = 1
        spin_id = 1
        ! 
        if (verbosity.ge.1) call print_title('WAVECAR')
        !
        if (verbosity.ge.1) write(*,'(a,a)') flare, 'file read = '//trim(fname)
        ! initialize record length
        recl=24
        open(unit=fid,file=trim(fname),access='direct',recl=recl,iostat=iostat,status='old')
            if (iostat.ne. 0) stop 'ERROR [read_wavecar]: problem reading WAVECAR. Check that file exists.'
            ! read record length, number of spins, and wavefunction precision
            read(unit=fid,rec=1) xrecl, xnspins, xprec
            ! recl
            recl = nint(xrecl)
            if (verbosity.ge.1) write(*,'(a,a)') flare, 'recl = '//tostring(recl)
            ! spins
            nspins = nint(xnspins)
            if (verbosity.ge.1) write(*,'(a,a)') flare, 'nspins = '//tostring(nspins)
            ! precision
            prec = nint(xprec)
            if (verbosity.ge.1) write(*,'(a,a)') flare, 'prec = '//tostring(prec)
            ! check if complex wavecar is sp
            if (prec.eq.45210) stop 'ERROR [read_wavecar]: WAVECAR must be in complex(sp), i.e. complex*16'
            if (spin_id.gt.nspins) stop 'ERROR [read_wavecar]: spin_id exceeds the number of spins'
        close(unit=fid)
        ! open again
        open(unit=fid,file=fname,access='direct',recl=recl,iostat=iostat,status='old')
            if (iostat.ne. 0) stop 'ERROR [read_wavecar]: problem reading WAVECAR. Check that file exists.'
            ! number of kpts, number of bands, cutoff energy, real-space basis
            read(unit=fid,rec=2) xnkpts, xnbands, ecut, ((bas(j,i),j=1,3),i=1,3)
            ! nkpts
            if (verbosity.ge.1) write(*,'(a,a)') flare, 'nkpts = '//tostring(nint(xnkpts))
            if (verbosity.ge.1) write(*,'(a,a)') flare, 'nbands = '//tostring(nint(xnbands))
            ! output
            if (present(nkpts)) nkpts = nint(xnkpts)
            if (present(nbands)) nbands = nint(xnbands)
            ! nband
            if (verbosity.ge.1) then
                write(*,'(a,a)') flare, 'bas = '
                call disp_indent()
                call disp(X=bas)
            endif
            ! get reciprocal cell volume
            vol = det(bas)
            ! stdout
            if (verbosity.ge.1) write(*,'(a,a)') flare, 'vol = '//tostring(vol)
            ! get reciprocal basis (two pi convetion used here)
            recbas = inv(bas) * 2*pi
            ! recbas
            if (verbosity.ge.1) then
                write(*,'(a,a)') flare, 'recbas (using 2pi convention) = '
                call disp_indent()
                call disp(X=recbas)
            endif
            ! estimate number of plane waves
            if (present(n_fft)) then
                n_fft = estimate_n_fft(ecut=ecut,recbas=recbas)
                ! stdout
                if (verbosity.ge.1) write(*,'(a,a)') flare, 'n_fft = '//tostring(n_fft)
            endif
        close(fid)
    end subroutine query_wavecar

    subroutine     read_wavecar(rpt_frac,kpt_id,band_id,spin_id,psi, occ,E,kpt, fname,verbosity)
        !
        implicit none
        !
        real(dp)    , intent(in), optional :: rpt_frac(:,:) ! only required for psi output
        integer     , intent(in) :: kpt_id  ! kpt,band,spin index of wavefunction to read
        integer     , intent(in) :: band_id ! kpt,band,spin index of wavefunction to read
        integer     , intent(in) :: spin_id ! kpt,band,spin index of wavefunction to read
        character(*), intent(in) :: fname
        integer     , intent(in) :: verbosity
        complex(sp) , intent(out), optional :: psi(:) !  wave function defined on rpt_frac
        real(dp)    , intent(out), optional :: kpt(3) ! kpoint coordinate
        real(dp)    , intent(out), optional :: occ(:) ! occ( nbands )
        real(dp)    , intent(out), optional :: E(:)   ! E( nbands )
        real(dp) :: a(3,3), recbas(3,3)
        real(dp) :: x_kpt(3)
        real(dp) :: ecut
        real(dp) :: vol
        real(dp)   , allocatable :: x_occ(:)
        complex(dp), allocatable :: x_E(:)
        complex(sp), allocatable :: coeff_nG(:)
        real(dp) :: xrecl,xnspins,xprec,x_nkpts,x_nbands,x_npws
        integer  :: recl, nspins, prec, nkpts, nbands, npws
        integer, allocatable :: G(:,:) ! G(1:3,nGs)
        integer :: n_fft(3)
        integer :: irec, iband, l
        integer :: i, j ! loop variables
        integer :: fid
        integer :: iostat
        !
        if (verbosity.ge.1) call print_title('WAVECAR')
        !
        if (verbosity.ge.1) write(*,'(a,a)') flare, 'file read = '//trim(fname)
        ! initialize FID
        fid = 1
        ! initialize record length
        recl=24
        open(unit=fid,file=trim(fname),access='direct',recl=recl,iostat=iostat,status='old')
            if (iostat.ne. 0) stop 'ERROR [read_wavecar]: problem reading WAVECAR. Check that file exists.'
            read(unit=fid,rec=1) xrecl, xnspins, xprec
            ! recl
            recl = nint(xrecl)
            if (verbosity.ge.1) write(*,'(a,a)') flare, 'recl = '//tostring(recl)
            ! spins
            nspins = nint(xnspins)
            if (verbosity.ge.1) write(*,'(a,a)') flare, 'nspins = '//tostring(nspins)
            ! precision
            prec = nint(xprec)
            if (verbosity.ge.1) write(*,'(a,a)') flare, 'prec = '//tostring(prec)
            ! check if complex wavecar is sp
            if (prec.eq.45210) stop 'ERROR [read_wavecar]: WAVECAR must be in complex(sp), i.e. complex*16'
            if (spin_id.gt.nspins) stop 'ERROR [read_wavecar]: spin_id exceeds the number of spins'
        close(unit=fid)
        ! open again
        open(unit=fid,file=fname,access='direct',recl=recl,iostat=iostat,status='old')
            if (iostat.ne. 0) stop 'ERROR [read_wavecar]: problem reading WAVECAR. Check that file exists.'
            ! number of kpts, number of bands, cutoff energy, real-space basis
            read(unit=fid,rec=2) x_nkpts, x_nbands, ecut, ((a(j,i),j=1,3),i=1,3)
            ! nkpts
            nkpts = nint(x_nkpts)
            if (verbosity.ge.1) write(*,'(a,a)') flare, 'nkpts = '//tostring(nkpts)
            ! nbands
            nbands = nint(x_nbands)
            if (verbosity.ge.1) write(*,'(a,a)') flare, 'nbands = '//tostring(nbands)
            ! checks
            if (kpt_id.gt.nkpts)   stop 'ERROR [read_wavecar]: kpt_id exceeds the number of kpoints'
            if (band_id.gt.nbands) stop 'ERROR [read_wavecar]: band_id exceeds the number of bands'
            ! stdout
            if (verbosity.ge.1) then
                write(*,'(a,a)') flare, 'bas = '
                call disp_indent()
                call disp(X=a)
            endif
            ! get reciprocal cell volume
            vol = det(a)
            ! stdout
            if (verbosity.ge.1) write(*,'(a,a)') flare, 'vol = '//tostring(vol)
            ! get reciprocal basis
            recbas = inv(a) * 2*pi
            ! recbas
            if (verbosity.ge.1) then
                write(*,'(a,a)') flare, 'recbas (using 2pi convention) = '
                call disp_indent()
                call disp(X=recbas)
            endif
            ! estimate number of plane waves
            n_fft = estimate_n_fft(ecut=ecut,recbas=recbas)
            ! stdout
            if (verbosity.ge.1) write(*,'(a,a)') flare, 'n_fft = '//tostring(n_fft)
            ! find the wave function
            irec=3+(kpt_id-1)*(nbands+1)+(spin_id-1)*nkpts*(nbands+1)
            ! allocate space for energies and occupancies
            allocate(x_occ(nbands))
            allocate(x_E(nbands))
            ! read the number of plane waves, the reciprocal lattice coordinate, and the band_id energy and occupancy
            read(unit=fid,rec=irec) x_npws, (x_kpt(i),i=1,3), (x_E(iband),x_occ(iband),iband=1,nbands)
            !
            if (present(occ)) occ = x_occ
            if (present(E))   E = real(x_E,dp)
            if (present(kpt)) kpt = x_kpt
            ! construct wave function
            if (present(psi)) then
                ! check for rpt_frac
                if (.not.present(rpt_frac)) stop 'ERROR [read_wavecar]: rpt_frac must be present when psi is called.'
                ! convert to integer
                npws = nint(x_npws)
                ! stdout
                if (verbosity.ge.1) write(*,'(a,a)') flare, 'npws = '//tostring(npws)
                ! get reciprocal lattice vectors offset by x_kpt which have kinetic energies below cutoff
                G = get_G(n_fft=n_fft, recbas=recbas, x_kpt=x_kpt, ecut=ecut, npws=npws)
                ! set record to planewave expansion coefficients corresponding to a particular band
                irec=irec+band_id
                ! allocate number of planewaves for coefficients
                allocate(coeff_nG(npws))
                ! plane-wave expansion coefficients depend on band index and kpoint, as expected wavefunctions corresponding to eigenvectors
                read(unit=fid,rec=irec) (coeff_nG(l), l=1 , npws)
                ! Fourier transform plane waves to obtain real-space wavefunctions
                psi = get_psi(vol=vol,x_kpt=x_kpt,coeff_nG=coeff_nG,G=G,rpt_frac=rpt_frac)
            endif
        close(fid)
        contains
        function       get_G(n_fft,recbas,x_kpt,ecut,npws) result(G)
            ! Count the number of plane waves that have energy below the plane-wave cutoff energy and checks 
            ! that it matches the value written by vasp; record reciprocal lattice points corresponding to these planewaves
            implicit none
            !
            integer , intent(in) :: n_fft(3)
            real(dp), intent(in) :: recbas(3,3)
            real(dp), intent(in) :: x_kpt(3)
            real(dp), intent(in) :: ecut ! planewave kinetic-energy cutoff
            integer , intent(in) :: npws ! number of planewaves
            integer, allocatable :: G(:,:) ! G(1:3,nGs)
            integer, allocatable :: gx(:)
            integer, allocatable :: gy(:)
            integer, allocatable :: gz(:)
            real(dp) :: kpt(3)
            integer :: i, j, k, n ! loop variables
            real(dp) :: c
            real(dp) :: etot
            ! allocate space
            allocate(G(3,npws))
            ! constant 'c' is 2m/hbar**2 in units of 1/eV Ang^2 (value is adjusted in final decimal places to agree with 
            ! VASP value; program checks for discrepancy of any results between this and VASP values)
            ! hbar/2m = 0.26246582250210965422d0 / eV Ang^2
            c = 0.262465831d0 ! to match vasp
            ! Generate the FFT grid used by VASP
            gx = fftshift([-n_fft(1):n_fft(1)])
            gy = fftshift([-n_fft(2):n_fft(2)])
            gz = fftshift([-n_fft(3):n_fft(3)])
            ! initialize planewave counter
            n=0
            ! loop over reciprocal lattice points
            do k = 1, 2*n_fft(3)+1
            do j = 1, 2*n_fft(2)+1
            do i = 1, 2*n_fft(1)+1
                kpt  = matmul(recbas, x_kpt + [gx(i),gy(j),gz(k)] )
                ! calculate the free-electron planewave energy
                ! etot = hbar^2 | k + G |^2 / 2m; note c = 1/(hbar^2/2m)
                etot = norm2(kpt)**2/c
                ! check whether it is below the planewave expansion cutoff energy
                if (etot.lt.ecut) then
                    n=n+1
                    G(1,n)=gx(i)
                    G(2,n)=gy(j)
                    G(3,n)=gz(k)
                endif
            enddo
            enddo
            enddo
            ! NOTE: THIS ROUTINE CAN TRIP AN OUT-OF-BOUNDS ERROR IF MORE PLANEWAVES ARE COUNTED THAN IDENTIFIED BY VASP
            ! ERROR: mismatched number of planewaves for which hbar^2 | k + G |^2 / 2m < E_cutoff
            ! the mismatch could be require a small adjustment in c value
            ! the mismatch could be due to using gamma version o vasp to generate wavecar
        end function   get_G
        function       get_psi(vol,x_kpt,coeff_nG,G,rpt_frac) result(psi)
            !
            implicit none
            !
            real(dp)   , intent(in) :: x_kpt(3)
            real(dp)   , intent(in) :: vol
            complex(sp), intent(in) :: coeff_nG(:)
            integer    , intent(in) :: G(:,:) ! G(1:3,nGs)
            real(dp)   , intent(in) :: rpt_frac(:,:)
            complex(sp), allocatable :: psi(:)
            integer  :: n
            integer  :: i, j
            ! get sizes
            n = size(rpt_frac,2)
            ! initialize wavefunction on realspace grid [0,1): rpt_frac
            allocate(psi(n))
            ! loop real space coordinates
            do i = 1, n
                psi(i) = 0.0_dp
                ! loop over planewaves
                do j = 1, npws
                    ! G is a list of reciprocal lattice points (integer kpt_id values in fractional coordinates; i.e. no forbidden points)
                    ! this is essentially performing the fourier transform to obtain the wave function in real space
                    psi(i) = psi(i) + coeff_nG(j) * exp(itwopi*dot_product( x_kpt + G(:,j), rpt_frac(1:3,i) ) )
                enddo
                psi(i) = psi(i)/sqrt(vol)
            enddo
        end function   get_psi
    end subroutine read_wavecar

    subroutine     read_doscar(efermi,nedos,dos,E,fname,verbosity)
        !
        implicit none
        !
        real(dp)             , intent(out), optional :: efermi
        integer              , intent(out), optional :: nedos
        real(dp), allocatable, intent(out), optional :: dos(:)
        real(dp), allocatable, intent(out), optional :: E(:)
        character(*)         , intent(in) , optional :: fname
        integer              , intent(in) , optional :: verbosity
        character(maximum_buffer_size) :: buffer
        character(len=:), allocatable  :: word(:)
        integer :: iostat
        integer :: x_nedos
        integer :: fid
        integer :: i
        !
        if (verbosity.ge.1) call print_title('DOSCAR')
        !
        fid = 1
        open(unit=fid,file=trim(fname),status='old',action='read',iostat=iostat)
            if (iostat.ne. 0) stop 'ERROR [read_doscar]: problem reading doscar. check if file exists'
            !
            if (verbosity.ge.1) write(*,'(a,a,a)') flare, 'file = ', trim(fname)
            ! (LINE 1) skip -    2    2    1    1
            read(fid,*)
            ! (LINE 2) skip - 0.8434557E+01  0.2878740E-09  0.2878740E-09  0.2878740E-09  0.5000000E-15
            read(fid,*)
            ! (LINE 3) skip - 1.000000000000000E-004
            read(fid,*)
            ! (LINE 4) skip - CAR
            read(fid,*)
            ! (LINE 5) skip - unknown system
            read(fid,*)
            ! (LINE 6) 53.80543798    -11.57073848  301      9.53971757      1.00000000
            read(unit=fid,fmt='(a)') buffer
            word = strsplit(buffer,delimiter=' ')
            ! save fermi energy
            if (present(efermi)) then
                read(word(4),*) efermi
                if (verbosity.ge.1) then
                    write(*,'(a,a,a)') flare, 'fermi energy = ', tostring(efermi)
                endif
            endif
            ! number of dos enegies
            read(word(3),*) x_nedos
            if (present(nedos)) nedos = x_nedos
            if (verbosity.ge.1) write(*,'(a,a,a)') flare, 'energies = ', tostring(x_nedos)
            ! allocate space for dos
            allocate(E(x_nedos))
            allocate(dos(x_nedos))
            ! (LINE 7) -11.571  0.0000E+00  0.0000E+00
            do i = 1, x_nedos
                read(unit=fid,fmt='(a)') buffer
                word = strsplit(buffer,delimiter=' ')
                read(word(1),*) E(i)
                read(word(1),*) dos(i)
            enddo
            !
        close(fid)
    end subroutine read_doscar

    subroutine     read_poscar(bas,natoms,nspecies,symb,tau_frac,atype,fname,verbosity)
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
        integer , intent(out) :: natoms
        integer , intent(out) :: nspecies
        character(len=:), allocatable, intent(out) :: symb(:) ! 1 symb per species
        integer , allocatable :: natoms_per_species(:) ! number of atoms per species
        real(dp), allocatable, intent(out) :: tau_frac(:,:) ! atomic coordinates tau_frac_read(3,natoms) in fractional
        integer , allocatable, intent(out) :: atype(:) ! type of atom tau_frac_read(natoms)
        character(*), intent(in) :: fname
        integer     , intent(in) :: verbosity
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
        ! output string
        character(:), allocatable :: disp_str(:)
        ! optional i/o
        integer :: iostat
        !
        call tostring_set(sep=' ')
        !
        if (verbosity.ge.1) call print_title('Unit cell')
        !
        fid = 1
        open(unit=fid,file=trim(fname),status='old',action='read',iostat=iostat)
            if (iostat.ne. 0) stop 'ERROR [read_poscar]: problem reading poscar. check if file exists'
            !
            if (verbosity.ge.1) write(*,'(a,a)') flare, 'file read = '//trim(fname)
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
            if ( abs(vol) .lt. tiny ) stop 'ERROR [read_poscar]: Basis vectors are not coplanar!'
            if (verbosity.ge.1) then 
                write(*,'(a,a)') flare, 'basis (column vectors) = '
                call disp_indent()
                call disp(X=bas)
            endif
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
            if (verbosity.ge.1) write(*,'(a,a)') flare, 'unique species = '//tostring(nspecies)
            if (verbosity.ge.1) write(*,'(a,a,100a)') flare, 'atomic species = ', (trim(symb(i)(1:2))//' ', i = 1, nspecies)
            ! (LINE 6) 1  1
            allocate(natoms_per_species(nspecies))
            read(fid,*) natoms_per_species
            if (verbosity.ge.1) write(*,'(a,a)') flare, 'atoms per species = '//tostring(natoms_per_species(1:nspecies))
            ! get total number of atoms
            natoms=sum(natoms_per_species)
            ! write the total number of atoms
            if (verbosity.ge.1) write(*,'(a,a)') flare, 'atoms = '//tostring(natoms)
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
                if (verbosity.ge.1) write(*,'(a,a)') flare, 'selective dynamics'
            endif
            if ( direct ) then
                if (verbosity.ge.1) write(*,'(a,a)') flare, 'atomic basis read as [frac]'
            endif
            if ( cartesian ) then
                if (verbosity.ge.1) write(*,'(a,a)') flare, 'atomic basis read as [cart]'
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
            ! write stdout
            if (verbosity.ge.1) then
                write(*,'(a,a)') flare, 'atomic basis = '
                call disp_indent()
                call disp(X=[1:natoms],style='underline',title='#',advance='no')
                allocate(character(2) :: disp_str(natoms))
                do i = 1, natoms
                    disp_str(i) = symb(atype(i))
                enddo
                call disp(X=disp_str                       ,style='underline',title='atom'  ,advance='no',fmt='a3')
                call disp(X=transpose(tau_frac)            ,style='underline',title='[frac]',advance='no')
                call disp(X=transpose(matmul(bas,tau_frac)),style='underline',title='[cart]',advance='yes')
            endif
        close(fid)
    end subroutine read_poscar

    subroutine     write_poscar(bas,natoms,nspecies,symb,tau_frac,atype,header,fname,verbosity)
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
        character(*), intent(in) :: header
        character(*), intent(in) :: fname
        integer     , intent(in) :: verbosity
        ! file i/o
        integer :: fid
        ! loop variables
        integer :: i, j, n
        ! optional i/o
        ! fname = 'outfile.POSCAR'
        ! header = 'POSCAR'
        !
        if (verbosity.ge.1) call print_title('Writing POSCAR')
        !
        fid = 1
        open(unit=fid,file=trim(fname),status='replace',action='write')
            !
            write(*,'(a,a,a)') flare, 'file = ', trim(fname)
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
            write(unit=fid,fmt='(x,a,x)',advance='no') tostring(count(atype.eq.i))
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

    ! wannier90

    subroutine     read_amn(U,nbands,nkpts,nwanniers,fname,verbosity)
        !>
        !> U is equal to a_matrix if there is no disentanglment.
        !>
        !>
        implicit none
        !
        complex(dp), intent(out), allocatable :: U(:,:,:)
        integer, intent(out) :: nbands, nkpts, nwanniers
        character(*), intent(in) :: fname
        integer     , intent(in) :: verbosity
        integer  :: nwanprojs
        real(dp) :: a_real, a_imag
        integer  :: i, m, n, k
        character(maximum_buffer_size) :: buffer ! read buffer
        character(len=:), allocatable :: word(:) ! read buffer
        integer :: iostat
        integer :: fid
        !
        ! fname = "wannier90.amn"
        !
        if (verbosity.ge.1) call print_title('Reading wannier unitary matrix (U matrix)')
        !
        fid = 1
        open(unit=fid,file=trim(fname),status='old',action='read',iostat=iostat)
            if (iostat.ne. 0) stop 'ERROR [read_amn]: problem reading wannier amn file. check if it exists'
            !
            if (verbosity.ge.1) write(*,'(a,a,a)') flare, 'file = ', trim(fname)
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
            if (verbosity.ge.1) then 
                write(*,'(a,a,a)') flare, 'bands = ', tostring(nbands)
                write(*,'(a,a,a)') flare, 'kpoints = ', tostring(nkpts)
                write(*,'(a,a,a)') flare, 'neighbors = ', tostring(nwanniers)
                write(*,'(a,a,a)') flare, 'projections = ', tostring(nwanprojs)
            endif
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

    subroutine     read_eig(E,nbands,nkpts,fname,verbosity)
        !>
        !> it seems that the U is equal to the a_matrix if there is no disentanglment.
        !>
        !>
        implicit none
        !
        real(dp), intent(out), allocatable :: E(:,:) !> E(nbands,nkpts) energies
        integer, intent(out) :: nbands
        integer, intent(out) :: nkpts
        character(*), intent(in) :: fname
        integer     , intent(in) :: verbosity
        real(dp) :: Ein
        integer  :: i, j, k
        character(maximum_buffer_size) :: buffer ! read buffer
        character(len=:), allocatable :: word(:) ! read buffer
        integer :: fid
        integer :: iostat
        !
        ! fname = "wannier90.eig"
        !
        if (verbosity.ge.1) call print_title('Reading eigenvalues for wannier')
        !
        ! SKIM FILE TO DETEMRINE PARAMETERS NECESSARY TO ALLOCATE SPACE
        !
        fid = 1
        open(unit=fid,file=trim(fname),status='old',action='read',iostat=iostat)
            if (iostat.ne. 0) stop 'ERROR [read_eig]: problem reading wannier eig file. check if it exists'
            !
            if (verbosity.ge.1) write(*,'(a,a,a)') flare, 'file = ', trim(fname)
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
        if (verbosity.ge.1) write(*,'(a,a,a)') flare, 'kpoints = ', tostring(nkpts)
        if (verbosity.ge.1) write(*,'(a,a,a)') flare, 'bands = ', tostring(nbands)
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


    ! auxiliary function

    function      estimate_n_fft(ecut,recbas) result(n_fft)
        !
        implicit none
        !
        real(dp), intent(in)  :: ecut
        real(dp), intent(in)  :: recbas(3,3)
        integer  :: n_fft(3)
        real(dp) :: vtmp(3)
        real(dp) :: c
        real(dp) :: phi12,phi13,phi23,phi123,sinphi123
        integer  :: x_nb1maxA, x_nb2maxA, x_nb3maxA
        integer  :: x_nb1maxB, x_nb2maxB, x_nb3maxB
        integer  :: x_nb1maxC, x_nb2maxC, x_nb3maxC
        integer  :: npmaxA , npmaxB , npmaxC
        real(dp) :: vmag
        ! constant 'c' is 2m/hbar**2 in units of 1/eV Ang^2 (value is adjusted in final decimal places to agree with 
        ! VASP value; program checks for discrepancy of any results between this and VASP values)
        ! hbar/2m = 0.26246582250210965422d0 / eV Ang^2
        c = 0.262465831d0 ! to match vasp
        ! estimate the number of planewaves along each primitive direction
        phi12=acos(dot_product(recbas(:,1),recbas(:,2))/(norm2(recbas(:,1))*norm2(recbas(:,2))))
        vtmp = cross_product(recbas(:,1),recbas(:,2))
        vmag=sqrt(vtmp(1)**2+vtmp(2)**2+vtmp(3)**2)
        sinphi123=(dot_product(recbas(:,3),Vtmp))/(vmag*norm2(recbas(:,3)))
        x_nb1maxA=(sqrt(ecut*c)/(norm2(recbas(:,1))*abs(sin(phi12))))+1
        x_nb2maxA=(sqrt(ecut*c)/(norm2(recbas(:,2))*abs(sin(phi12))))+1
        x_nb3maxA=(sqrt(ecut*c)/(norm2(recbas(:,3))*abs(sinphi123)))+1
        npmaxA=nint(4.*pi*x_nb1maxA*x_nb2maxA*x_nb3maxA/3.)
        ! repeat...
        phi13=acos(dot_product(recbas(:,1),recbas(:,3))/(norm2(recbas(:,1))*norm2(recbas(:,3))))
        vtmp = cross_product(recbas(:,1),recbas(:,3))
        vmag=sqrt(vtmp(1)**2+vtmp(2)**2+vtmp(3)**2)
        sinphi123=(dot_product(recbas(:,2),Vtmp))/(vmag*norm2(recbas(:,2)))
        phi123=abs(asin(sinphi123))
        x_nb1maxB=(sqrt(ecut*c)/(norm2(recbas(:,1))*abs(sin(phi13))))+1
        x_nb2maxB=(sqrt(ecut*c)/(norm2(recbas(:,2))*abs(sinphi123)))+1
        x_nb3maxB=(sqrt(ecut*c)/(norm2(recbas(:,3))*abs(sin(phi13))))+1
        npmaxB=nint(4.*pi*x_nb1maxB*x_nb2maxB*x_nb3maxB/3.)
        ! repeat...
        phi23=acos(dot_product(recbas(:,2),recbas(:,3))/(norm2(recbas(:,2))*norm2(recbas(:,3))))
        vtmp = cross_product(recbas(:,2),recbas(:,3))
        vmag=sqrt(vtmp(1)**2+vtmp(2)**2+vtmp(3)**2)
        sinphi123=(dot_product(recbas(:,1),Vtmp))/(vmag*norm2(recbas(:,1)))
        phi123=abs(asin(sinphi123))
        x_nb1maxC=(sqrt(ecut*c)/(norm2(recbas(:,1))*abs(sinphi123)))+1
        x_nb2maxC=(sqrt(ecut*c)/(norm2(recbas(:,2))*abs(sin(phi23))))+1
        x_nb3maxC=(sqrt(ecut*c)/(norm2(recbas(:,3))*abs(sin(phi23))))+1 
        npmaxC=nint(4.*pi*x_nb1maxC*x_nb2maxC*x_nb3maxC/3.)
        ! get the largest values in each of the three cases
        n_fft(1)=max(x_nb1maxA,x_nb1maxB,x_nb1maxC)
        n_fft(2)=max(x_nb2maxA,x_nb2maxB,x_nb2maxC)
        n_fft(3)=max(x_nb3maxA,x_nb3maxB,x_nb3maxC)
    end function  estimate_n_fft

!     subroutine read_mmn(fname,verbosity)
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
!         if ( present(fname) ) fname = fname
!         verbosity = 1
!         if ( present(verbosity) ) verbosity = verbosity
!         !
!         !
!         !
!         if (verbosity.ge.1) call print_title('Reading wannier overlaps (MMN)')
!         !
!         fid = 1
!         open(unit=fid,file=trim(fname),status="old",action='read')
!             !
!             if (verbosity.ge.1) call am_print('actual file read is',fname,flare)
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
!             if (verbosity.ge.1) call am_print('number of bands',nbands,flare)
!             if (verbosity.ge.1) call am_print('number of kpoints',nkpts,flare)
!             if (verbosity.ge.1) call am_print('number of neighbors',nneighbors,flare)
!             if (verbosity.ge.1) call am_print('number of overlaps',noverlaps,flare)
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
    