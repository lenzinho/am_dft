    !  src.f90
    !
    !  FUNCTIONS:
    !  src - Entry point of  application.
    !
    !
    !****************************************************************************
    !
    !  PROGRAM: src
    !
    !  PURPOSE:  Entry point for the console application.
    !
    !****************************************************************************
    program main
	!
	use am_helpers
    use am_options
	use am_unit_cell
    use am_brillouin_zone
    use am_symmetry
    use am_wannier
	!
	implicit none
	!
    type(am_class_symmetry) :: ss
    type(am_class_symmetry) :: pg
	type(am_class_unit_cell) :: uc
    type(am_class_options) :: opts
    type(am_class_tetrahedra) :: tet
    type(am_class_bz) :: bz
    type(am_class_bz) :: ibz
    type(am_class_wannier) :: wan
    ! testing
    type(am_class_bz) :: fbz
    !
    call opts%get
    !
    ! call uc%load_poscar(opts)
    !
    call bz%load_ibzkpt(opts=opts)
    !
    call wan%load_amn(opts=opts)
    !
    call wan%load_eig(opts=opts)

!     call determine_symmetry(uc=uc,ss=ss,pg=pg,opts=opts)
!     !
!     call fbz%full_monkhorst_pack_mesh(uc=uc,n=[31,31,31],s=real([0,0,0],dp),opts=opts)
!     !
!     call fbz%outfile_kpoints('outfile.full_monkhorst_pack')
!     !
!     call ibz%reduce_to_ibz(bz=fbz,uc=uc,pg=pg,opts=opts)
!     !
!     call ibz%outfile_kpoints('outfile.irre_monkhorst_pack')
!     !
    stop ! STOPPING HERE <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

    !>
    !> WORKS DECENTLY, the problems are:
    !> 1 ) does not reduce FBZ to IBZ properly
    !> 2 ) spurious peaks in DOS via tetrahedron method.
    !> 3 ) sum of partial density of states does not match total. It seems there is a band which is missing.
    !>

    call opts%get
    !
    call uc%load_poscar(opts)
    !
    call determine_symmetry(uc=uc,ss=ss,pg=pg,opts=opts)

    call bz%load_eigenval(opts)
    ! call bz%load_procar
    !
    call bz%outfile_bandcharacter(uc=uc,pg=pg,opts=opts)
    !
    ! <HAS IBZKPT>
    call tet%load_ibzkpt(iopt_bz=bz,opts=opts)
    ! <DOES NOT HAVE IBZKPT>
    ! call tet%get_tetrahedra(uc=uc,tet=tet,opts=opts)
    ! <END>
    !
    call bz%tetrahedra_pdos(tet=tet,opts=opts)
    !
    call bz%outfile_dosprojected(opts)
    !
    write(*,*) 'Done'
    

    end program