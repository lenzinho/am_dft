module am_dos

    use am_dispersion
    use am_tetrahedra
    use am_brillouin_zone
    use am_prim_cell
    use am_symmetry
    use am_constants
    use am_stdout
    use am_vasp_io
    use am_mkl
    use am_options

    implicit none

    type, public :: am_class_dos
        integer :: nEs
        real(dp), allocatable :: E(:)          ! energies
        real(dp), allocatable :: dos(:)        ! dos(nEs)
        real(dp), allocatable :: lmproj(:,:,:,:) ! lmproj(nEs,nspins,norbitals,nions) pdos weights
    contains
        procedure :: tetrahedra_pdos
        procedure :: outfile_dosprojected
    end type am_class_dos

contains

    subroutine     tetrahedra_pdos(dos,dr,bz,tet,opts)
        !
        use am_histogram
        !
        implicit none
        !
        class(am_class_dos)      , intent(out) :: dos
        type(am_class_dr)        , intent(in)  :: dr
        type(am_class_bz)        , intent(in)  :: bz
        type(am_class_tetrahedra), intent(in)  :: tet
        type(am_class_options)   , intent(in)  :: opts
        real(dp) :: maxE ! maximum band energy
        real(dp) :: minE ! minimum band energy
        real(dp) :: dE   ! energy increment
        real(dp) :: w(4) ! weights on tetrahedron tet
        real(dp) :: Ec(4) ! tet energies
        integer :: i, j, k, m, n, o
        !
        if (opts%verbosity.ge.1) call print_title('Tetrahedron integration')
        !
        dE = 0.01_dp
        minE = minval(dr%E(:,:))
        maxE = maxval(dr%E(:,:))
        ! dos%E = linspace(minE,maxE,100)
        dos%E = regspace(minE,maxE,dE)
        dos%nEs = size(dos%E,1)
        !
        call plot_histogram( histogram(pack(dr%E(:,:),(dr%E(:,:).ge.minE).and.(dr%E(:,:).le.maxE)),column_width) )
        !
        if (opts%verbosity.ge.1) call am_print('number of bands',           dr%nbands     ,' ... ')
        if (opts%verbosity.ge.1) call am_print('number of tetrahedra',      tet%ntets     ,' ... ')
        if (opts%verbosity.ge.1) call am_print('lowest band energy',        minE          ,' ... ')
        if (opts%verbosity.ge.1) call am_print('highest band energy',       maxE          ,' ... ')
        if (opts%verbosity.ge.1) call am_print('minimum probing energy',    dos%E(1)      ,' ... ')
        if (opts%verbosity.ge.1) call am_print('maximum probing energy',    dos%E(dos%nEs),' ... ')
        if (opts%verbosity.ge.1) call am_print('energy increments dE',      dE            ,' ... ')
        if (opts%verbosity.ge.1) call am_print('number of probing energies',dos%nEs       ,' ... ')
        !
        allocate(dos%lmproj(dr%nspins,dr%norbitals,dr%nions,dos%nEs))
        allocate(dos%dos(dos%nEs))
        !
        dos%lmproj = 0.0_dp
        dos%dos  = 0.0_dp
        !
        ! two things
        ! 1) difference between total dos and sum of JDOS. ONE band is being left out in the JDOS summation. figure out which one.
        !  0.996425152311378        1.00000000031254
        ! 2) delta summation produces wild values. the divisions by 1/tiny are numerically unstable. Figure out how to get rid of them.!
        !
        do i = 1, dos%nEs
            do j = 1, dr%nbands
            do k = 1, tet%ntets
                !
                Ec(1:4) = dr%E(j,tet%tet(:,k))
                !
                w = tetrahedron_weights(x=dos%E(i),xc=Ec,integration_type='delta',apply_blochl_corrections=.true.)
                !
                dos%dos(i) = dos%dos(i) + tet%w(k)*tet%volume(k)*sum(w)
                !
                ! <projections>
                !
                do m = 1,dr%nspins
                do n = 1,dr%norbitals
                do o = 1,dr%nions
                    !> dos%lmproj(nEs,nspins,norbitals,nions)
                    !> dr%lmproj(nspins,norbitals,nions,nbands,nkpts)
                    dos%lmproj(i,m,n,o) = dos%lmproj(i,m,n,o) + tet%w(k)*tet%volume(k)*sum(w * dr%lmproj(m,n,o,j,tet%tet(:,k)) )/real(dr%nspins,dp) ! per spin
                enddo
                enddo
                enddo
                !
                ! </projections>
                !
            enddo
            enddo
            write(*,*) dos%E(i), sum(dos%lmproj(:,:,:,i)), dos%dos(i)
        enddo
        !
    end subroutine tetrahedra_pdos
    
    subroutine     outfile_dosprojected(dos,dr,opts)
        ! 
        ! Creates outfile.dosprojected
        !
        use am_options
        !
        implicit none
        !
        class(am_class_dos)   , intent(in) :: dos
        type(am_class_dr)     , intent(in)  :: dr
        type(am_class_options), intent(in) :: opts
        integer :: fid
        integer :: i, l, m, n ! loop variables
        !
        if (opts%verbosity.ge.1) call print_title('Writing dosprojected file')
        !
        fid = 1
        open(unit=fid,file="outfile.dosprojected",status="replace",action='write')
            !
            ! STDOUT
            !
            if (opts%verbosity.ge.1) write(*,*) ' ... Projection information'
            if (opts%verbosity.ge.1) call am_print('number of ions',     dr%nions ,' ... ')
            if (opts%verbosity.ge.1) call am_print('number of orbitals', dr%norbitals ,' ... ')
            if (opts%verbosity.ge.1) call am_print('number of spins',    dr%nspins ,' ... ')
            if (opts%verbosity.ge.1) call am_print('number of columns',  dr%nspins*dr%norbitals ,' ... ')
            !
            ! HEADER (first four lines)
            !
            write(fid,'(a)')  'Projected density of states [States/eV/spin/unit-cell]'
            write(fid,'(100a13)',advance='no')  'ions'
            do m = 1, dr%nions
            do l = 1, dr%norbitals
            do n = 1, dr%nspins
                write(fid,'(i13)' ,advance='no') m
            enddo
            enddo
            enddo
            write(fid,*)
            !
            write(fid,'(100a13)',advance='no') 'spin'
            do m = 1, dr%nions
            do l = 1, dr%norbitals
            do n = 1, dr%nspins
                write(fid,'(i13)' ,advance='no') n
            enddo
            enddo
            enddo
            write(fid,*)
            !
            write(fid,'(100a13)',advance='no') 'E [eV]'
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
            do i = 1, dos%nEs
                !
                write(fid,'(f13.5)',advance='no') dos%E(i)
                do m = 1, dr%nions
                do l = 1, dr%norbitals
                do n = 1, dr%nspins
                    !> lmproj(nspins,norbitals,nions,nbands,nkpts)
                    write(fid,'(f13.5)',advance='no') dos%lmproj(n,l,m,i)
                enddo
                enddo
                enddo
                write(fid,*)
                !
            enddo
        close(fid)
        !
    end subroutine outfile_dosprojected


end module am_dos