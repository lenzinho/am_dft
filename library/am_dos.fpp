#:include "fypp_macros.fpp"
module am_dos

    use am_dispersion
    use am_brillouin_zone
    use am_symmetry
    use am_constants
    use am_stdout
    use am_vasp_io
    use am_mkl
    use am_options

    implicit none

    type, public :: am_class_dos
        integer :: nEs
        real(dp), allocatable :: E(:) ! energies
        real(dp), allocatable :: D(:) ! dos(nEs)
        ! real(dp), allocatable :: lmproj(:,:,:,:) ! lmproj(nEs,nspins,norbitals,nions) pdos weights
        contains
        procedure :: get_dos
    end type am_class_dos

contains

    subroutine     get_dos(dos,dr,ibz,E,opts)
        !
        use am_histogram
        !
        implicit none
        !
        class(am_class_dos)       , intent(out) :: dos
        class(am_class_dispersion), intent(in)  :: dr
        type(am_class_ibz)        , intent(in)  :: ibz
        real(dp), optional        , intent(in)  :: E(:)
        type(am_class_options)    , intent(in)  :: opts
        real(dp) :: Emax ! maximum band energy
        real(dp) :: Emin ! minimum band energy
        real(dp) :: dE   ! energy increment
        real(dp) :: w(4) ! weights on tetrahedron tet
        real(dp) :: Ec(4) ! tet energies
        integer  :: nbands
        integer  :: i,j,k
        !
        if (opts%verbosity.ge.1) call print_title('Density of states')
        ! region over which there are bands
        Emin = minval(dr%E(:,:))
        Emax = maxval(dr%E(:,:))
        ! check if an energy range has been input
        if (present(E)) then
            allocate(dos%E,source=E)
            dE = E(2)-E(1)
            dos%nEs = size(E)
        else
            ! if not, create one an energy sampling with 1 meV spacing
            dE = 0.01_dp
            dos%E = regspace(Emin,Emax,dE)
            dos%nEs = size(dos%E,1)
        endif
        ! get number of bands
        nbands = size(dr%E,1)
        ! preliminar stuff to stdout
        if (opts%verbosity.ge.1) then 
            write(*,'(a,a)') flare, 'bands = '//tostring(nbands)
            write(*,'(a,a)') flare, 'tetrahedra = '//tostring(ibz%ntets)
            write(*,'(a,a)') flare, 'lowest band energy = '//tostring(Emin)
            write(*,'(a,a)') flare, 'highest band energy = '//tostring(Emax)
            write(*,'(a,a)') flare, 'probing energies = '//tostring(dos%nEs)
            write(*,'(a,a)') flare, 'lower energy range = '//tostring(dos%E(1))
            write(*,'(a,a)') flare, 'upper energy range = '//tostring(dos%E(dos%nEs))
            write(*,'(a,a)') flare, 'energy increment dE = '//tostring(dE)
        endif
        ! allocate space for dos
        allocate(dos%D(dos%nEs))
        !
        ! two things
        ! 1) difference between total dos and sum of JDOS. ONE band is being left out in the JDOS summation. figure out which one.
        !  0.996425152311378        1.00000000031254
        ! 2) delta summation produces wild values. the divisions by 1/tiny are numerically unstable. Figure out how to get rid of them.!
        !
        do i = 1, dos%nEs
            dos%D(i) = 0.0_dp
            do j = 1, nbands
            do k = 1, ibz%ntets
                ! get energies at the corner of the tetrahedra
                Ec(:) = dr%E(j,ibz%tet(:,k))
                ! get tetrahedron corner weights
                w = get_tetrahedron_weight(E=dos%E(i),Ec=Ec,flags='delta,blochl')
                ! increment DOS with contributions from band j in tetrahedron k
                dos%D(i) = dos%D(i) + ibz%tetw(k) * sum(w)
                ! !
                ! ! <projections>
                ! !
                ! do m = 1,dr%nspins
                ! do n = 1,dr%norbitals
                ! do o = 1,dr%nions
                !     !> dos%lmproj(nEs,nspins,norbitals,nions)
                !     !> dr%lmproj(nspins,norbitals,nions,nbands,nkpts)
                !     dos%lmproj(i,m,n,o) = dos%lmproj(i,m,n,o) + tet%w(k)*tet%volume(k)*sum(w * dr%lmproj(m,n,o,j,tet%tet(:,k)) )/real(dr%nspins,dp) ! per spin
                ! enddo
                ! enddo
                ! enddo
                ! !
                ! ! </projections>
                ! !
            enddo
            enddo
            write(*,*) dos%E(i), dos%D(i)
        enddo
        !
    end subroutine get_dos


end module am_dos













