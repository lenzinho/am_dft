module am_dispersion

    use dispmodule
    use am_brillouin_zone
    use am_unit_cell
    use am_constants
    use am_stdout
    use am_vasp_io
    use am_mkl
    use am_options
    use am_histogram

    implicit none

    private

    type, public :: am_class_dispersion
        integer :: nbands
        real(dp), allocatable :: E(:,:) ! E(nbands,nkpts) energies
    end type am_class_dispersion

   type, public, extends(am_class_dispersion) :: am_class_dispersion_vasp
        real(dp), allocatable :: lmproj(:,:,:,:,:) ! lmproj(nspins,norbitals,nions,nbands,nkpts) band character weights
        ! projection information:
        integer :: nspins
        integer :: norbitals
        integer :: nions
        character(:), allocatable :: orbitals(:) ! orbitals(norbitals) names of orbitals
    contains
        procedure :: load
        procedure :: write_bandcharacter ! requires load_eigenval or load_procar
    end type am_class_dispersion_vasp

contains

    subroutine     load(dr,opts,flags)
        !
        ! use text_viz
        !
        implicit none
        !
        class(am_class_dispersion_vasp), intent(out) :: dr
        type(am_class_options)         , intent(in)  :: opts
        character(*)                   , intent(in)  :: flags
        ! integer, allocatable :: canvas(:,:)
        ! integer :: i, height
        ! integer :: minmin
        ! integer :: maxmax
        ! integer :: minmax
        !
        !
        if (opts%verbosity.ge.1) call print_title('Input electronic band dispersion')
        !
        if     (index(flags,'eigenval')) then
            write(*,'(a,a,a)') flare, 'file = ', 'eigenval'
            call read_eigenval(E=dr%E, iopt_filename=opts%eigenval, iopt_verbosity=0, &
                nbands=dr%nbands, nspins=dr%nspins, lmproj=dr%lmproj)
        elseif (index(flags,'procar')) then
            write(*,'(a,a,a)') flare, 'file = ', 'procar'
            call read_procar(E=dr%E, iopt_filename=opts%procar, iopt_verbosity=0, &
                nbands=dr%nbands, nions=dr%nions, nspins=dr%nspins, norbitals=dr%norbitals, orbitals=dr%orbitals, lmproj=dr%lmproj)
        else
            stop 'ERROR: eigenval/procar flag required.'
        endif
        !
        if (opts%verbosity.ge.1) then
            write(*,'(a,a,a)') flare, 'bands = ', tostring(dr%nbands)
            write(*,'(a,a,a)') flare, 'spins = ', tostring(dr%nspins)
            write(*,'(a,a)  ') flare, 'band energy statistics = '
            write(*,'(5x,a,a)')          'lowest = ', tostring(minval(pack(dr%E,.true.)))
            write(*,'(5x,a,a)')          'highest = ', tostring(maxval(pack(dr%E,.true.)))
            write(*,'(5x,a,a)')          'mean = ', tostring( sum(pack(dr%E,.true.))/size(pack(dr%E,.true.)) )
            ! make nice dos histogram
            write(*,'(a,a)') flare, 'density of states'
            call plot_histogram(binned_data=histogram(x=pack(dr%E,.true.),m=98)) 
            ! make nice plot
            ! write(*,'(a,a)') flare, 'band ranges = '
            ! height = 98
            ! allocate(canvas(dr%nbands,height))
            ! canvas = 0
            ! minmin = minval(pack(dr%E,.true.))
            ! maxmax = maxval(pack(dr%E,.true.))
            ! minmax = maxmax-minmin
            ! do i = 1, dr%nbands
            ! canvas(i,min(max(nint((minval(dr%E(i,:))-minmin)/real(minmax,dp)*height),1),height)) = -1
            ! canvas(i,min(max(nint((maxval(dr%E(i,:))-minmin)/real(minmax,dp)*height),1),height)) = +1
            ! enddo
            ! ! make nice pot
            ! call viz_init(height,dr%nbands)
            ! call viz_clear
            ! call viz_plot(transpose(real(canvas)),height,dr%nbands)
            ! call viz_done
        endif
        !
        if (dr%nbands.eq.0) stop 'ERROR: nbands = 0'
        !
    end subroutine load

    subroutine     write_bandcharacter(dr,bz,uc,opts)
        ! 
        ! Creates outfile.bandcharacter, which contains the energy and character of the bands read either from EIGENVAL using bz%load_eigenval() or PROCAR using bz%load_procar()
        !
        use am_options
        !
        implicit none
        !
        class(am_class_dispersion_vasp), intent(in) :: dr ! dispersion relations
        class(am_class_bz)             , intent(in) :: bz ! brillouin zone class
        class(am_class_unit_cell)      , intent(in) :: uc
        type(am_class_options)         , intent(in) :: opts
        real(dp), allocatable :: grid_points(:,:) !> voronoi points (27=3^3)
        character(:), allocatable :: fname
        real(dp) :: k1(3) ! point for kdiff
        real(dp) :: k2(3) ! point for kdiff
        real(dp) :: kdiff ! for dispersion plot
        integer :: fid
        integer :: i, j, l, m, n ! loop variables
        !
        if (opts%verbosity.ge.1) call print_title('Bandcharacter output')
        !
        allocate(fname,source='outfile.bandcharacter')
        write(*,'(a,a)') flare, 'file = '//trim(fname)
        !
        ! generate voronoi points (cartesian), used to reduce points to wigner-seitz brillouin zone
        grid_points = matmul( inv(uc%bas) , meshgrid([-1:1],[-1:1],[-1:1]) )
        !
        fid = 1
        open(unit=fid,file=trim(fname),status="replace",action='write')
            ! STDOUT
            if (opts%verbosity.ge.1) write(*,'(a,a)') flare, 'kpoints' //tostring(bz%nkpts)
            if (opts%verbosity.ge.1) write(*,'(a,a)') flare, 'bands'   //tostring(dr%nbands)
            if (opts%verbosity.ge.1) write(*,'(a,a)') flare, 'ions'    //tostring(dr%nions)
            if (opts%verbosity.ge.1) write(*,'(a,a)') flare, 'orbitals'//tostring(dr%norbitals)
            if (opts%verbosity.ge.1) write(*,'(a,a)') flare, 'spins'   //tostring(dr%nspins)
            if (opts%verbosity.ge.1) write(*,'(a,a)') flare, 'columns' //tostring(dr%nspins*dr%norbitals)
            !
            ! HEADER (first three lines)
            write(fid,'(a)')  'Character-projected energy dispersion'
            write(fid,'(100a13)',advance='no') ' ', flare,' ','ions'
            do m = 1, dr%nions
                do l = 1, dr%norbitals
                do n = 1, dr%nspins
                    write(fid,'(i13)' ,advance='no') m
                enddo
                enddo
            enddo
            write(fid,*)
            !
            write(fid,'(100a13)',advance='no') ' ', flare,' ','spin'
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
                        k1=bz%kpt_frac(:,i-1)
                        k2=bz%kpt_frac(:,i)
                        kdiff = kdiff + norm2( abs(k1-k2) )
                    endif
                    !
                    write(fid,'(100f13.5)',advance='no') kdiff, bz%kpt_frac(1:3,i), dr%E(j,i)
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