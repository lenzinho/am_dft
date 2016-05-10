module am_dos

    implicit none

    type, public :: am_class_dos
        integer :: nEs
        real(dp), allocatable :: E(:)       ! energies
        real(dp), allocatable :: dos(:)     ! dos(nEs)
        real(dp), allocatable :: w(:,:,:)   ! w(nEs,norbitals,nions) pdos weights
    contains
        ! load vasp files into bz class
        procedure :: load_procar
        procedure :: load_eigenval
        procedure :: load_ibzkpt => load_ibzkpt_bz
        ! bz stuff
        procedure :: expand_to_fbz
        procedure :: reduce_to_ibz
        ! kpoints
        procedure :: full_monkhorst_pack_mesh
        procedure :: outfile_kpoints
        ! high level i/o
        procedure :: outfile_bandcharacter ! requires load_eigenval or load_procar
        procedure :: outfile_dosprojected
        !
    end type am_class_dos

contains

    subroutine     tetrahedra_pdos(dos,bz,dr,tet,opts)
        !
        use am_histogram
        use am_options
        !
        implicit none
        !
        class(am_class_bz), intent(inout) :: bz
        type(am_class_tetrahedra), intent(in) :: tet
        type(am_class_options), intent(in) :: opts
        real(dp) :: maxE ! maximum band energy
        real(dp) :: minE ! minimum band energy
        real(dp) :: dE   ! energy increment
        real(dp) :: w(4) ! weights on tetrahedron tet
        real(dp) :: Ec(4) ! tet energies
        integer :: i, j, k, m, n, o
        !
        if (opts%verbosity.ge.1) call am_print_title('Tetrahedron integration')
        !
        dE = 0.01_dp
        minE = minval(bz%E(:,:))
        maxE = maxval(bz%E(:,:))
        ! bz%Ep = linspace(minE,maxE,100)
        bz%Ep = regspace(minE,maxE,dE)
        bz%nEp = size(bz%Ep,1)
        !
        call plot_histogram( histogram(pack(bz%E(:,:),(bz%E(:,:).ge.minE).and.(bz%E(:,:).le.maxE)),column_width) )
        !
        if (opts%verbosity.ge.1) call am_print('number of bands',bz%nbands,' ... ')
        if (opts%verbosity.ge.1) call am_print('number of tetrahedra',tet%ntets,' ... ')
        if (opts%verbosity.ge.1) call am_print('lowest band energy',minE,' ... ')
        if (opts%verbosity.ge.1) call am_print('highest band energy',maxE,' ... ')
        if (opts%verbosity.ge.1) call am_print('minimum probing energy',bz%Ep(1),' ... ')
        if (opts%verbosity.ge.1) call am_print('maximum probing energy',bz%Ep(bz%nEp),' ... ')
        if (opts%verbosity.ge.1) call am_print('energy increments dE',dE,' ... ')
        if (opts%verbosity.ge.1) call am_print('number of probing energies',bz%nEp,' ... ')
        !
        allocate(bz%pdos(bz%nspins,bz%norbitals,bz%nions,bz%nEp))
        allocate(bz%dos(bz%nEp))
        !
        bz%pdos = 0.0_dp
        bz%dos = 0.0_dp
        !
        ! two things
        ! 1) difference between total dos and sum of JDOS. ONE band is being left out in the JDOS summation. figure out which one.
        !  0.996425152311378        1.00000000031254
        ! 2) delta summation produces wild values. the divisions by 1/tiny are numerically unstable. Figure out how to get rid of them.!
        !
        do i = 1, bz%nEp
            do j = 1, bz%nbands
            do k = 1, tet%ntets
                !
                Ec(1:4) = bz%E(j,tet%tet(:,k))
                !
                w = tetrahedron_weights(x=bz%Ep(i),xc=Ec,integration_type='delta',apply_blochl_corrections=.true.)
                !
                bz%dos(i) = bz%dos(i) + tet%w(k)*tet%volume(k)*sum(w)
                !
                ! <projections>
                !
                do m = 1,bz%nspins
                do n = 1,bz%norbitals
                do o = 1,bz%nions
                    !> pdos(nspins,norbitals,nions,nEP)
                    !> lmproj(nspins,norbitals,nions,nbands,nkpts)
                    bz%pdos(m,n,o,i) = bz%pdos(m,n,o,i) + tet%w(k)*tet%volume(k)*sum(w*bz%lmproj( m,n,o,j,tet%tet(:,k) ))/real(bz%nspins,dp) ! per spin
                enddo
                enddo
                enddo
                !
                ! </projections>
                !
            enddo
            enddo
            write(*,*) bz%Ep(i), sum(bz%pdos(:,:,:,i)), bz%dos(i)
        enddo
        !
    end subroutine tetrahedra_pdos
    
    subroutine     outfile_dosprojected(bz,opts)
        ! 
        ! Creates outfile.dosprojected
        !
        use am_options
        !
        implicit none
        !
        class(am_class_bz), intent(inout) :: bz !> brillouin zone class
        type(am_class_options), intent(in) :: opts
        integer :: fid
        integer :: i, l, m, n ! loop variables
        !
        if (opts%verbosity.ge.1) call am_print_title('Writing dosprojected file')
        !
        fid = 1
        open(unit=fid,file="outfile.dosprojected",status="replace",action='write')
            !
            ! STDOUT
            !
            if (opts%verbosity.ge.1) write(*,*) ' ... Projection information'
            if (opts%verbosity.ge.1) call am_print('number of ions', bz%nions ,' ... ')
            if (opts%verbosity.ge.1) call am_print('number of orbitals', bz%norbitals ,' ... ')
            if (opts%verbosity.ge.1) call am_print('number of spins', bz%nspins ,' ... ')
            if (opts%verbosity.ge.1) call am_print('number of columns', bz%nspins*bz%norbitals ,' ... ')
            !
            ! HEADER (first four lines)
            !
            write(fid,'(a)')  'Projected density of states [States/eV/spin/unit-cell]'
            write(fid,'(100a13)',advance='no')  'ions'
            do m = 1, bz%nions
                do l = 1, bz%norbitals
                do n = 1, bz%nspins
                    write(fid,'(i13)' ,advance='no') m
                enddo
                enddo
            enddo
            write(fid,*)
            !
            write(fid,'(100a13)',advance='no') 'spin'
            do m = 1, bz%nions
                do l = 1, bz%norbitals
                do n = 1, bz%nspins
                    write(fid,'(i13)' ,advance='no') n
                enddo
                enddo
            enddo
            write(fid,*)
            !
            write(fid,'(100a13)',advance='no') 'E [eV]'
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
            !real(dp), allocatable :: E(:)
            !real(dp), allocatable :: pdos(:,:,:,:) !> pdos(nspins,norbitals,nions,nEprobing)
            !real(dp), allocatable :: dos(:) !> dos(nEprobing)
            
            do i = 1, bz%nEs
                !
                write(fid,'(f13.5)',advance='no') bz%E(i)
                do m = 1, bz%nions
                    do l = 1, bz%norbitals
                        do n = 1, bz%nspins
                            !> lmproj(nspins,norbitals,nions,nbands,nkpts)
                            write(fid,'(f13.5)',advance='no') bz%pdos(n,l,m,i)
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