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
	!
	implicit none
	!
    type(am_class_symmetry) :: ss
    type(am_class_symmetry) :: pg
	type(am_class_unit_cell) :: uc
    type(am_class_options) :: opts
    type(am_class_bz) :: bz
    type(am_class_tetrahedra) :: tet
    !
    opts%verbosity = 1
    opts%sym_prec  = tiny
    !
    call uc%load_poscar('POSCAR')
    !
    call determine_symmetry(uc=uc,ss=ss,pg=pg,iopts=opts)
    !
    call bz%load_eigenval
    ! call bz%load_procar
    !
    call bz%outfile_bandcharacter(uc=uc,pg=pg,iopts=opts)
    !
    ! <HAS IBZKPT>
    call tet%load_ibzkpt(iopt_bz=bz,iopts=opts)
    ! <DOES NOT HAVE IBZKPT>
    ! call bz%get_tetrahedra(uc=uc,tet=tet,iopts=opts)
    ! <END>
    !
    call bz%tetrahedra_pdos(tet=tet,iopts=opts)
    !
    call bz%outfile_dosprojected
    !
    write(*,*) 'Done'
    

    end program