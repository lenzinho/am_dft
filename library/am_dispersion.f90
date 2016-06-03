module am_dispersion

    use am_brillouin_zone
    use am_unit_cell
    use am_constants
    use am_stdout
    use am_vasp_io
    use am_mkl
    use am_options

    implicit none

    type, public :: am_class_dr
        !
        real(dp), allocatable :: E(:,:)            ! E(nbands,nkpts) energies
        real(dp), allocatable :: lmproj(:,:,:,:,:) ! lmproj(nspins,norbitals,nions,nbands,nkpts) band character weights
        !
        ! projection information:
        integer :: nspins
        integer :: norbitals
        integer :: nions
        integer :: nbands
        integer :: nkpts
        character(:), allocatable :: orbitals(:) ! orbitals(norbitals) names of orbitals
    contains
        procedure :: load_procar
        procedure :: load_eigenval
        procedure :: write_bandcharacter ! requires load_eigenval or load_procar
        ! 
    end type am_class_dr

contains

    subroutine     load_procar(dr,opts)
        ! 
        ! Reads the PROCAR file.
        !
        implicit none
        !
        class(am_class_dr), intent(out) :: dr
        type(am_class_options), intent(in) :: opts
        !
        call read_procar(nkpts    = dr%nkpts,&
                         nbands   = dr%nbands,&
                         nions    = dr%nions,&
                         nspins   = dr%nspins,&
                         norbitals= dr%norbitals,&
                         orbitals = dr%orbitals,&
                         lmproj   = dr%lmproj,&
                         E        = dr%E,&
                         iopt_filename =opts%procar,&
                         iopt_verbosity=opts%verbosity)
        !
    end subroutine load_procar

    subroutine     load_eigenval(dr,opts)
        ! 
        ! Reads the eigenval file.
        !
        implicit none
        !
        class(am_class_dr), intent(out) :: dr
        type(am_class_options), intent(in) :: opts
        !
        call read_eigenval(nkpts  = dr%nkpts,&
                           nbands = dr%nbands,&
                           nspins = dr%nspins,&
                           E      = dr%E,&
                           lmproj = dr%lmproj,&
                           iopt_filename=opts%eigenval,&
                           iopt_verbosity=opts%verbosity)
        !
        dr%nions=1
        dr%norbitals=1
        allocate(character(15) :: dr%orbitals(1))
        dr%orbitals='tot'
        !
    end subroutine load_eigenval

    subroutine     write_bandcharacter(dr,bz,uc,opts)
        ! 
        ! Creates outfile.bandcharacter, which contains the energy and character of the bands read either from EIGENVAL using bz%load_eigenval() or PROCAR using bz%load_procar()
        !
        use am_options
        !
        implicit none
        !
        class(am_class_dr), intent(in) :: dr ! dispersion relations
        class(am_class_bz), intent(in) :: bz ! brillouin zone class
        class(am_class_unit_cell), intent(in) :: uc
        type(am_class_options), intent(in) :: opts
        real(dp), allocatable :: grid_points(:,:) !> voronoi points (27=3^3)
        real(dp) :: k1(3) ! point for kdiff
        real(dp) :: k2(3) ! point for kdiff
        real(dp) :: kdiff ! for dispersion plot
        integer :: fid
        integer :: i, j, l, m, n ! loop variables
        !
        if (opts%verbosity.ge.1) call print_title('Writing bandcharacter file')
        !
        !
        ! generate voronoi points (cartesian), used to reduce points to wigner-seitz brillouin zone
        grid_points = matmul( inv(uc%bas) , meshgrid([-1:1],[-1:1],[-1:1]) )
        !
        fid = 1
        open(unit=fid,file="outfile.bandcharacter",status="replace",action='write')
            !
            ! STDOUT
            !
            if (opts%verbosity.ge.1) call am_print('kpoints',  bz%nkpts,' ... ')
            if (opts%verbosity.ge.1) call am_print('bands',    dr%nbands ,' ... ')
            if (opts%verbosity.ge.1) call am_print('ions',     dr%nions ,' ... ')
            if (opts%verbosity.ge.1) call am_print('orbitals', dr%norbitals ,' ... ')
            if (opts%verbosity.ge.1) call am_print('spins',    dr%nspins ,' ... ')
            if (opts%verbosity.ge.1) call am_print('columns',  dr%nspins*dr%norbitals ,' ... ')
            !
            ! HEADER (first three lines)
            !
            write(fid,'(a)')  'Character-projected energy dispersion'
            write(fid,'(100a13)',advance='no') ' ', ' ',' ',' ','ions'
            do m = 1, dr%nions
                do l = 1, dr%norbitals
                do n = 1, dr%nspins
                    write(fid,'(i13)' ,advance='no') m
                enddo
                enddo
            enddo
            write(fid,*)
            !
            write(fid,'(100a13)',advance='no') ' ', ' ',' ',' ','spin'
            do m = 1, dr%nions
                do l = 1, dr%norbitals
                do n = 1, dr%nspins
                    write(fid,'(i13)' ,advance='no') n
                enddo
                enddo
            enddo
            write(fid,*)
            !
            write(fid,'(100a13)',advance='no') 'kpath','k_x','k_y','k_z','E [eV]'
            do m = 1, dr%nions
                do l = 1, dr%norbitals
                do n = 1, dr%nspins
                    write(fid,'(a13)',advance='no') trim(dr%orbitals(l))
                enddo
                enddo
            enddo 
            write(fid,*)
            !
            ! DATA
            !
            do j = 1, dr%nbands
                do i = 1, bz%nkpts
                    if (i .eq. 1) then
                        kdiff = 0
                    else
                        k1=bz%kpt(:,i-1)
                        k2=bz%kpt(:,i)
                        kdiff = kdiff + norm2( abs(k1-k2) )
                    endif
                    !
                    write(fid,'(100f13.5)',advance='no') kdiff, bz%kpt(1:3,i), dr%E(j,i)
                    do m = 1, dr%nions
                    do l = 1, dr%norbitals
                    do n = 1, dr%nspins
                        !> lmproj(nbands,nkpts,nspins,norbitals,nions)
                        write(fid,'(f13.5)',advance='no') dr%lmproj(n,l,m,j,i)
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