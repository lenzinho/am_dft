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