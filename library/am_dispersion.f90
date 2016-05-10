module am_dispersion

    use am_brillouin_zone
    use am_prim_cell
    use am_symmetry
    use am_constants
    use am_stdout
    use am_vasp_io
    use am_mkl
    use am_options

    implicit none

    type, public :: am_class_disp
        !
        integer :: nbands
        integer :: nkpts
        !
        real(dp), allocatable :: E(:,:)     ! E(nbands,nkpts) energies
        real(dp), allocatable :: w(:,:,:,:) ! w(nbands,nkpts,nspins,norbitals,nions) band character weights
        ! projection information:
        integer :: nspins
        integer :: norbitals
        integer :: nions
        character(:), allocatable :: orbitals(:) ! orbitals(norbitals) names of orbitals
    contains
        procedure :: load_procar
        procedure :: load_eigenval
        procedure :: write_bandcharacter ! requires load_eigenval or load_procar
        ! 
    end type am_class_disp

contains

    subroutine     load_procar(dr,bz,opts)
        ! 
        ! Reads the PROCAR file.
        !
        implicit none
        !
        class(am_class_disp), intent(out) :: dr
        class(am_class_bz)  , intent(out) :: bz
        type(am_class_options), intent(in) :: opts
        !
        call read_procar(nkpts    = dr%nkpts,&
                         nbands   = dr%nbands,&
                         nions    = dr%nions,&
                         nspins   = dr%nspins,&
                         norbitals= dr%norbitals,&
                         orbitals = dr%orbitals,&
                         lmproj   = dr%w,&
                         E        = bz%E,&
                         kpt      = bz%kpt,&
                         w        = bz%w,&
                         iopt_filename =opts%procar,&
                         iopt_verbosity=opts%verbosity)
        !
    end subroutine load_procar

    subroutine     load_eigenval(dr,bz,opts)
        ! 
        ! Reads the eigenval file.
        !
        implicit none
        !
        class(am_class_disp), intent(out) :: dr
        class(am_class_bz)  , intent(out) :: bz
        type(am_class_options), intent(in) :: opts
        !
        call read_eigenval(nkpts  = dr%nkpts,&
                           nbands = dr%nbands,&
                           nspins = dr%nspins,&
                           nelecs = dr%nelecs,&
                           lmproj = dr%w,&
                           kpt    = bz%kpt,&
                           w      = bz%w,&
                           E      = bz%E,&
                           iopt_filename=opts%eigenval,&
                           iopt_verbosity=opts%verbosity)
        !
        dr%nions=1
        dr%norbitals=1
        allocate(character(15) :: dr%orbitals(1))
        dr%orbitals='tot'
        !
    end subroutine load_eigenval

    subroutine     write_bandcharacter(bz,pc,pg,opts)
        ! 
        ! Creates write.bandcharacter, which contains the energy and character of the bands read either from EIGENVAL using bz%load_eigenval() or PROCAR using bz%load_procar()
        !
        use am_options
        !
        implicit none
        !
        class(am_class_disp), intent(in) :: dr ! brillouin zone class
        class(am_class_bz), intent(in) :: bz ! brillouin zone class
        type(am_class_prim_cell), intent(in) :: pc
        type(am_class_symmetry) , intent(in) :: pg
        type(am_class_options), intent(in) :: opts
        real(dp), allocatable :: grid_points(:,:) !> voronoi points (27=3^3)
        real(dp) :: k1(3) ! point for kdiff
        real(dp) :: k2(3) ! point for kdiff
        real(dp) :: kdiff ! for dispersion plot
        integer :: fid
        integer :: i, j, l, m, n ! loop variables
        !
        if (opts%verbosity.ge.1) call am_print_title('Writing bandcharacter file')
        !
        !
        ! generate voronoi points (cartesian), used to reduce points to wigner-seitz brillouin zone
        grid_points = matmul( inv(pc%bas) , meshgrid([-1:1],[-1:1],[-1:1]) )
        !
        fid = 1
        open(unit=fid,file="write.bandcharacter",status="replace",action='write')
            !
            ! STDOUT
            !
            if (opts%verbosity.ge.1) call am_print('number of kpoints', bz%nkpts,' ... ')
            if (opts%verbosity.ge.1) call am_print('number of bands', bz%nbands ,' ... ')
            if (opts%verbosity.ge.1) call am_print('number of ions', bz%nions ,' ... ')
            if (opts%verbosity.ge.1) call am_print('number of orbitals', bz%norbitals ,' ... ')
            if (opts%verbosity.ge.1) call am_print('number of spins', bz%nspins ,' ... ')
            if (opts%verbosity.ge.1) call am_print('number of columns', bz%nspins*bz%norbitals ,' ... ')
            !
            ! HEADER (first three lines)
            !
            write(fid,'(a)')  'Character-projected energy dispersion'
            write(fid,'(100a13)',advance='no') ' ', ' ',' ',' ','ions'
            do m = 1, bz%nions
                do l = 1, bz%norbitals
                do n = 1, bz%nspins
                    write(fid,'(i13)' ,advance='no') m
                enddo
                enddo
            enddo
            write(fid,*)
            !
            write(fid,'(100a13)',advance='no') ' ', ' ',' ',' ','spin'
            do m = 1, bz%nions
                do l = 1, bz%norbitals
                do n = 1, bz%nspins
                    write(fid,'(i13)' ,advance='no') n
                enddo
                enddo
            enddo
            write(fid,*)
            !
            write(fid,'(100a13)',advance='no') 'kpath','k_x','k_y','k_z','E [eV]'
            do m = 1, bz%nions
                do l = 1, bz%norbitals
                do n = 1, bz%nspins
                    write(fid,'(a13)',advance='no') trim(bz%orbitals(l))
                enddo
                enddo
            enddo 
            write(fid,*)
            !
            ! DATA
            !
            do j = 1, bz%nbands
                do i = 1, bz%nkpts
                    if (i .eq. 1) then
                        kdiff = 0
                    else
                        k1=bz%kpt(:,i-1)
                        k2=bz%kpt(:,i)
                        ! k1 = reduce_kpoint_to_ibz(kpoint=bz%kpt(:,i-1),grid_points=grid_points,bas=pc%bas,R=pg%R)
                        ! k2 = reduce_kpoint_to_ibz(kpoint=bz%kpt(:,i),grid_points=grid_points,bas=pc%bas,R=pg%R)
                        kdiff = kdiff + norm2( abs(k1-k2) )
                    endif
                    !
                    write(fid,'(100f13.5)',advance='no') kdiff, bz%kpt(1:3,i), bz%E(j,i)
                    do m = 1, bz%nions
                        do l = 1, bz%norbitals
                        do n = 1, bz%nspins
                            !> lmproj(nspins,norbitals,nions,nbands,nkpts)
                            write(fid,'(f13.5)',advance='no') bz%lmproj(n,l,m,j,i)
                        enddo
                        enddo
                    enddo
                    write(fid,*)
                    !
                enddo
                write(fid,'(a)') ' '
            enddo
        close(fid)
        !
    end subroutine write_bandcharacter


end module am_dispersion