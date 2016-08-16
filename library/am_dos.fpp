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
        integer :: nprojections 
        real(dp), allocatable :: E(:)    ! energies
        real(dp), allocatable :: D(:)    ! dos(nEs)
        real(dp), allocatable :: pD(:,:) ! projected dos pD(nEs)
        contains
        procedure :: save => save_dos
        procedure :: load => load_dos
        procedure, private :: write_dos
        procedure, private :: read_dos
        generic :: write(formatted) => write_dos
        generic :: read(formatted)  => read_dos
    end type am_class_dos

contains


    ! i/o

    subroutine     write_dos(dtv, unit, iotype, v_list, iostat, iomsg)
        !
        implicit none
        !
        class(am_class_dos), intent(in) :: dtv
        integer     , intent(in)    :: unit
        character(*), intent(in)    :: iotype
        integer     , intent(in)    :: v_list(:)
        integer     , intent(out)   :: iostat
        character(*), intent(inout) :: iomsg
        !
        iostat = 0
        !
        if (iotype.eq.'LISTDIRECTED') then
            write(unit,'(a/)') '<dos>'
                ! non-allocatable
                #:for ATTRIBUTE in ['nEs','nprojections']
                    $:write_xml_attribute_nonallocatable(ATTRIBUTE)
                #:endfor
                ! allocatable        
                #:for ATTRIBUTE in ['E','D','pD']
                    $:write_xml_attribute_allocatable(ATTRIBUTE)
                #:endfor
            write(unit,'(a/)') '</dos>'
        else
            stop 'ERROR [write_tb]: iotype /= LISTDIRECTED'
        endif
    end subroutine write_dos

    subroutine     read_dos(dtv, unit, iotype, v_list, iostat, iomsg)
        !
        implicit none
        !
        class(am_class_dos), intent(inout) :: dtv
        integer     , intent(in)    :: unit
        character(*), intent(in)    :: iotype
        integer     , intent(in)    :: v_list(:)
        integer     , intent(out)   :: iostat
        character(*), intent(inout) :: iomsg
        logical :: isallocated
        integer :: dims_rank
        integer :: dims(10) ! read tensor up to rank 5
        !
        if (iotype.eq.'LISTDIRECTED') then
            read(unit=unit,fmt='(/)',iostat=iostat,iomsg=iomsg)
                ! non-allocatable
                #:for ATTRIBUTE in ['nEs','nprojections']
                    $:read_xml_attribute_nonallocatable(ATTRIBUTE)
                #:endfor
                ! allocatable        
                #:for ATTRIBUTE in ['E','D','pD']
                    $:read_xml_attribute_allocatable(ATTRIBUTE)
                #:endfor
            read(unit=unit,fmt='(/)',iostat=iostat,iomsg=iomsg)
            ! set iostat
            iostat=-1
        else
            stop 'ERROR [read_dos]: iotype /= LISTDIRECTED'
        endif
    end subroutine read_dos

    subroutine     save_dos(dos,fname)
        !
        implicit none
        !
        class(am_class_dos), intent(in) :: dos
        character(*), intent(in) :: fname
        integer :: fid
        integer :: iostat
        ! fid
        fid = 1
        ! save space group
        open(unit=fid, file=trim(fname), status='replace', action='write', iostat=iostat)
            if (iostat/=0) stop 'ERROR [dos:load]: opening file'
            write(fid,*) dos
        close(fid)
        !
    end subroutine save_dos

    subroutine     load_dos(dos,fname)
        !
        implicit none
        !
        class(am_class_dos), intent(inout) :: dos
        character(*), intent(in) :: fname
        integer :: fid
        integer :: iostat
        ! fid
        fid = 1
        ! save space group
        open(unit=fid, file=trim(fname), status='old', action='read', iostat=iostat)
            if (iostat/=0) stop 'ERROR [dos:load]: opening file'
            read(fid,*) dos
        close(fid)
        !
    end subroutine load_dos

    ! gaussian method

    function       get_dos_gauss(Ep,E,kptw,sigma) result(D)
        !
        implicit none
        !
        real(dp), intent(in) :: Ep(:)    ! probing energies
        real(dp), intent(in) :: E(:,:)   ! E(nbands,nkpts) band energies
        real(dp), intent(in) :: kptw(:)  ! kptw(nkpts) normalized kpoint weights
        real(dp), intent(in) :: sigma    ! smearing
        real(dp), allocatable :: D(:)
        integer  :: nEs
        integer  :: nbands
        integer  :: nkpts
        integer  :: i,j,k
        ! get number of probing energies
        nEs = size(Ep)
        ! get n of band
        nbands = size(E,1)
        ! get n of kpoints
        nkpts = size(kptw)
        ! allocate space for dos
        allocate(D(nEs))
        ! loop over energies/bands/kpoints
        !$OMP PARALLEL PRIVATE(i,j,k) SHARED(D,E,Ep,nEs,sigma,kptw)
        !$OMP DO
        do i = 1, nEs
            ! initialize D(i)
            D(i) = 0.0_dp
            ! print progress
            call show_progress(iteration=i,maximum=nEs)
            do j = 1, nbands
            do k = 1, nkpts
                ! get contribution
                D(i) = D(i) + kptw(k) * exp( -0.5_dp * ( (Ep(i)-E(j,k))/sigma )**2 )
            enddo
            enddo
        enddo
        !$OMP END DO
        !$OMP END PARALLEL
    end function   get_dos_gauss

    ! tetrahedra methods

    function       get_dos_tetra(Ep,E,tet,tetw) result(D)
        !
        implicit none
        !
        real(dp), intent(in) :: Ep(:) ! probing energies
        real(dp), intent(in) :: E(:,:) ! E(nbands,nkpts) band energies
        integer , intent(in) :: tet(:,:) ! tet(:,ntets) tetrahedra conenctivity
        real(dp), intent(in) :: tetw(:) ! tetw(ntets) weight of each tetrahedron
        real(dp), allocatable :: D(:)
        integer  :: nEs
        integer  :: nbands
        integer  :: ntets
        real(dp) :: w(4) ! weights on tetrahedron tet
        integer  :: i,j,k
        ! get number of probing energies
        nEs = size(Ep)
        ! get n of band
        nbands = size(E,1)
        ! number of tetrahedra
        ntets = size(tet,2)
        ! allocate space for dos
        allocate(D(nEs))
        !$OMP PARALLEL PRIVATE(i,j,k,w) SHARED(D,E,Ep,nEs,tetw)
        !$OMP DO
        do i = 1, nEs
            ! initialize D(i)
            D(i) = 0.0_dp
            ! print progress
            call show_progress(iteration=i,maximum=nEs)
            do j = 1, nbands
            do k = 1, ntets
                ! get tetrahedron corner weights
                w = get_tetrahedron_weight(E=Ep(i),Ec=E(j,tet(:,k)),flags='delta,blochl')
                ! increment DOS with contributions from band j in tetrahedron k
                D(i) = D(i) + tetw(k) * sum(w)
            enddo
            enddo
        enddo
        !$OMP END DO
        !$OMP END PARALLEL
    end function   get_dos_tetra

    function       get_pdos_tetra(Ep,E,tet,tetw,weight) result(pD)
        !
        implicit none
        !
        real(dp), intent(in) :: Ep(:) ! probing energies
        real(dp), intent(in) :: E(:,:) ! E(nbands,nkpts) band energies
        integer , intent(in) :: tet(:,:) ! tet(:,ntets) tetrahedra conenctivity
        real(dp), intent(in) :: tetw(:) ! tetw(ntets) weight of each tetrahedron
        real(dp), intent(in) :: weight(:,:,:) ! weights(nprojections,nbands,nkpts), weights (per band per kpoint) for projected dos: can be absolute square of TB or BvK eigenvectors
        real(dp), allocatable :: pD(:,:)
        integer  :: nEs
        integer  :: nbands
        integer  :: ntets
        integer  :: nprojections
        real(dp) :: w(4) ! weights on tetrahedron tet
        integer  :: i,j,k,m
        ! get number of probing energies
        nEs = size(Ep)
        ! get n of band
        nbands = size(E,1)
        ! get number of tetrahedra
        ntets = size(tet,2)
        ! get number of projections
        nprojections = size(weight,1)
        ! allocate space 
        allocate(pD(nprojections,nEs))
        ! loop over energies, bands, tetrahedra
        !$OMP PARALLEL PRIVATE(i,j,k,m,w) SHARED(pD,E,Ep,nEs,tetw,weight)
        !$OMP DO
        do i = 1, nEs
            ! initialize D(i)
            pD(:,i) = 0.0_dp
            ! print progress
            call show_progress(iteration=i,maximum=nEs)
            do j = 1, nbands
            do k = 1, ntets
                ! get tetrahedron corner weights
                w = get_tetrahedron_weight(E=Ep(i),Ec=E(j,tet(:,k)),flags='delta,blochl')
                ! loop over projections, increment pDOS with contributions from band j in tetrahedron k
                do m = 1, nprojections
                    pD(m,i) = pD(m,i) + tetw(k) * sum(w * weight(m,j,tet(:,k)) )
                enddo
            enddo
            enddo
        enddo
        !$OMP END DO
        !$OMP END PARALLEL
        !
    end function   get_pdos_tetra

end module am_dos













