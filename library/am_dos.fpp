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
        real(dp), allocatable :: iD(:)   ! integrated dos(nEs)
        real(dp), allocatable :: pD(:,:) ! projected dos pD(nEs)
        contains
        procedure :: save => save_dos
        procedure :: load => load_dos
        procedure, private :: write_dos
        procedure, private :: read_dos
        generic :: write(formatted) => write_dos
        generic :: read(formatted)  => read_dos
    end type am_class_dos

    public :: get_fermi_level, get_integrated_dos

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

    ! DOS properties

    function       get_integrated_dos(dos) result(idos)
        !
        implicit none
        !
        class(am_class_dos), intent(in) :: dos
        real(dp), allocatable :: idos(:) ! integrated dos
        ! space
        allocate(idos(dos%nEs))
        ! cumsum
        idos = cumsum(dos%D)*(dos%E(2)-dos%E(1))
        !
    end function   get_integrated_dos

    function       get_fermi_level(dos,nelecs) result(EF)
        !
        implicit none
        !
        class(am_class_dos), intent(in) :: dos
        real(dp), intent(in) :: nelecs ! number of electrons
        real(dp), allocatable :: idos(:) ! integrated dos
        real(dp) :: EF
        integer  :: i
        ! get integrated dos
        idos = get_integrated_dos(dos=dos)
        ! get fermi energy
        EF = dos%E( minpos(abs(idos-nelecs)) )
        !
    end function   get_fermi_level

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

    function       get_tetrahedron_weight(E,Ec,flags) result(w)
        !
        ! flags = heavi/delta, blochl
        !
        use am_rank_and_sort
        !
        implicit none
        !
        real(dp)    , intent(in) :: Ec(4)
        character(*), intent(in) :: flags
        real(dp) :: w(4)
        real(dp) :: E,e1,e2,e3,e4,f
        real(dp) :: Ei(4)
        integer :: bracket
        integer :: i
        !
        ! sort Ec into Ei
        ! Ei = Ec(sort_four_doubles(Ec))
        Ei = Ec
        if (Ei(1).gt.Ei(2)) then; Ei([1,2]) = Ei([2,1]); endif
        if (Ei(3).gt.Ei(4)) then; Ei([3,4]) = Ei([4,3]); endif
        if (Ei(1).gt.Ei(3)) then; Ei([1,3]) = Ei([3,1]); endif
        if (Ei(2).gt.Ei(4)) then; Ei([2,4]) = Ei([4,2]); endif
        if (Ei(2).gt.Ei(3)) then; Ei([2,3]) = Ei([3,2]); endif
        e1=Ei(1); e2=Ei(2); e3=Ei(3); e4=Ei(4)
        !
        if (e1.gt.e2) stop 'ERROR [tetrahedron_weights]: sorting failed'
        if (e2.gt.e3) stop 'ERROR [tetrahedron_weights]: sorting failed'
        if (e3.gt.e4) stop 'ERROR [tetrahedron_weights]: sorting failed'
        !
        bracket = 0
        if     (E.lt.e1) then; bracket = 1
        elseif (E.lt.e2) then; bracket = 2
        elseif (E.lt.e3) then; bracket = 3
        elseif (E.lt.e4) then; bracket = 4
        elseif (E.gt.e4) then; bracket = 5
        else
            stop 'ERROR [tetrahedron_weights]: bracketing failed'
        endif
        !
        if (index(flags,'heavi').ne.0) then
            select case(bracket)
            case (1)
                w(1) = 0.0_dp
                w(2) = 0.0_dp
                w(3) = 0.0_dp
                w(4) = 0.0_dp
                f    = 0.0_dp
            case (2)
                w(1) = -((E-e1)**3*((E-e1)*(1.0_dp/(e1-e2)+1.0_dp/(e1-e3)+1.0_dp/(e1-e4))+4.0_dp))/((e1-e2)*(e1-e3)*(e1-e4))
                w(2) =  ((E-e1)**4*1.0_dp/(e1-e2)**2)/((e1-e3)*(e1-e4))
                w(3) =  ((E-e1)**4*1.0_dp/(e1-e3)**2)/((e1-e2)*(e1-e4))
                w(4) =  ((E-e1)**4*1.0_dp/(e1-e4)**2)/((e1-e2)*(e1-e3))
                f    =  ((E-e1)**2*(-1.2D1))/((e1-e2)*(e1-e3)*(e1-e4))
            case (3)
                w(1) =  (E-e1)**2/((e1-e3)*(e1-e4))+(((E-e1)**2/((e1-e3)*(e1-e4))+((E-e1)*(E-e2)*(E-e3))/((e1-e3)*(e1-e4)*(e2-e3)))*(E-e3))/(e1-e3)+((E-e4)*((E-e1)**2/((e1-e3)*(e1-e4))+((E-e2)**2*(E-e4))/((e1-e4)*(e2-e3)*(e2-e4))+((E-e1)*(E-e2)*(E-e3))/((e1-e3)*(e1-e4)*(e2-e3))))/(e1-e4)
                w(2) =  ((((E-e2)**2*(E-e4))/((e1-e4)*(e2-e3)*(e2-e4))+((E-e1)*(E-e2)*(E-e3))/((e1-e3)*(e1-e4)*(e2-e3)))*(E-e3))/(e2-e3)+(E-e1)**2/((e1-e3)*(e1-e4))+((E-e2)**2*(E-e4))/((e1-e4)*(e2-e3)*(e2-e4))+((E-e2)**2*(E-e4)**2*1.0_dp/(e2-e4)**2)/((e1-e4)*(e2-e3))+((E-e1)*(E-e2)*(E-e3))/((e1-e3)*(e1-e4)*(e2-e3))
                w(3) = -((((E-e2)**2*(E-e4))/((e1-e4)*(e2-e3)*(e2-e4))+((E-e1)*(E-e2)*(E-e3))/((e1-e3)*(e1-e4)*(e2-e3)))*(E-e2))/(e2-e3)-(((E-e1)**2/((e1-e3)*(e1-e4))+((E-e1)*(E-e2)*(E-e3))/((e1-e3)*(e1-e4)*(e2-e3)))*(E-e1))/(e1-e3)
                w(4) = -((E-e1)*((E-e1)**2/((e1-e3)*(e1-e4))+((E-e2)**2*(E-e4))/((e1-e4)*(e2-e3)*(e2-e4))+((E-e1)*(E-e2)*(E-e3))/((e1-e3)*(e1-e4)*(e2-e3))))/(e1-e4)-((E-e2)**3*(E-e4)*1.0_dp/(e2-e4)**2)/((e1-e4)*(e2-e3))
                f    =  (E*2.4D1-e1*1.2D1-e2*1.2D1+((E-e2)**2*(e1*3.0_dp+e2*3.0_dp-e3*3.0_dp-e4*3.0_dp)*4.0_dp)/((e2-e3)*(e2-e4)))/((e1-e3)*(e1-e4))
            case (4)
                w(1) = -((E-e4)**4*1.0_dp/(e1-e4)**2)/((e2-e4)*(e3-e4))+1.0_dp
                w(2) = -((E-e4)**4*1.0_dp/(e2-e4)**2)/((e1-e4)*(e3-e4))+1.0_dp
                w(3) = -((E-e4)**4*1.0_dp/(e3-e4)**2)/((e1-e4)*(e2-e4))+1.0_dp
                w(4) =  ((E-e4)**3*((E-e4)*(1.0_dp/(e1-e4)+1.0_dp/(e2-e4)+1.0_dp/(e3-e4))-4.0_dp))/((e1-e4)*(e2-e4)*(e3-e4))+1.0_dp
                f    =  ((E-e4)**2*(-1.2D1))/((e1-e4)*(e2-e4)*(e3-e4))
            case (5)
                w(1) = 1.0_dp
                w(2) = 1.0_dp
                w(3) = 1.0_dp
                w(4) = 1.0_dp
                f    = 0.0_dp
            case default
                stop 'ERROR [tetrahedron_weights]: heaviside, invalid bracket'
            end select
        elseif (index(flags,'delta').ne.0) then
            select case(bracket)
            case (1)
                w(1) = 0.0_dp
                w(2) = 0.0_dp
                w(3) = 0.0_dp
                w(4) = 0.0_dp
                f    = 0.0_dp
            case (2)
                w(1) =   (E-e1)**2*1.0_dp/(e1-e2)**2*1.0_dp/(e1-e3)**2*1.0_dp/(e1-e4)**2*(E*(e4*(e2+e3)+e2*e3-e1*(e2+e3+e4)*2.0_dp+e1**2*3.0_dp)+e1*(e4*(e2+e3)+e2*e3)*2.0_dp-e1**2*(e2+e3+e4)-e2*e3*e4*3.0_dp)*(-4.0_dp)
                w(2) =  ((E-e1)**3*1.0_dp/(e1-e2)**2*4.0_dp)/((e1-e3)*(e1-e4))
                w(3) =  ((E-e1)**3*1.0_dp/(e1-e3)**2*4.0_dp)/((e1-e2)*(e1-e4))
                w(4) =  ((E-e1)**3*1.0_dp/(e1-e4)**2*4.0_dp)/((e1-e2)*(e1-e3))
                f    = -(E*2.4D1-e1*2.4D1)/((e1-e2)*(e1-e3)*(e1-e4))
            case (3)
                w(1) = (1.0_dp/(e1-e3)**2*1.0_dp/(e1-e4)**2*(E**2*(e1**2*e2*3.0_dp+e3*e4**2*3.0_dp+e3**2*e4*3.0_dp-e1*e3*e4*6.0_dp-e2*e3*e4*3.0_dp)-E**3*(e1*e2*2.0_dp-e1*e3*2.0_dp-e1*e4*2.0_dp-e2*e3-e2*e4+e3*e4+e1**2+e3**2+e4**2)-E*(e3**2*e4**2*3.0_dp+e1**2*e2*e3*3.0_dp+e1**2*e2*e4*3.0_dp-e1**2*e3*e4*3.0_dp-e1*e2*e3*e4*6.0_dp)+e1**2*e2*e3**2+e1**2*e2*e4**2+e1*e3**2*e4**2*2.0_dp-e1**2*e3*e4**2-e1**2*e3**2*e4+e2*e3**2*e4**2-e1*e2*e3*e4**2*2.0_dp-e1*e2*e3**2*e4*2.0_dp+e1**2*e2*e3*e4)*(-4.0_dp))/((e2-e3)*(e2-e4))
                w(2) = (1.0_dp/(e2-e3)**2*1.0_dp/(e2-e4)**2*(E**2*(e1*e2**2*3.0_dp+e3*e4**2*3.0_dp+e3**2*e4*3.0_dp-e1*e3*e4*3.0_dp-e2*e3*e4*6.0_dp)-E**3*(e1*e2*2.0_dp-e1*e3-e1*e4-e2*e3*2.0_dp-e2*e4*2.0_dp+e3*e4+e2**2+e3**2+e4**2)-E*(e3**2*e4**2*3.0_dp+e1*e2**2*e3*3.0_dp+e1*e2**2*e4*3.0_dp-e2**2*e3*e4*3.0_dp-e1*e2*e3*e4*6.0_dp)+e1*e2**2*e3**2+e1*e2**2*e4**2+e1*e3**2*e4**2+e2*e3**2*e4**2*2.0_dp-e2**2*e3*e4**2-e2**2*e3**2*e4-e1*e2*e3*e4**2*2.0_dp-e1*e2*e3**2*e4*2.0_dp+e1*e2**2*e3*e4)*(-4.0_dp))/((e1-e3)*(e1-e4))
                w(3) = (1.0_dp/(e1-e3)**2*1.0_dp/(e2-e3)**2*(E**2*(e1*e2**2*3.0_dp+e1**2*e2*3.0_dp+e3**2*e4*3.0_dp-e1*e2*e3*6.0_dp-e1*e2*e4*3.0_dp)-E**3*(e1*e2-e1*e3*2.0_dp-e1*e4-e2*e3*2.0_dp-e2*e4+e3*e4*2.0_dp+e1**2+e2**2+e3**2)-E*(e1**2*e2**2*3.0_dp-e1*e2*e3**2*3.0_dp+e1*e3**2*e4*3.0_dp+e2*e3**2*e4*3.0_dp-e1*e2*e3*e4*6.0_dp)-e1*e2**2*e3**2-e1**2*e2*e3**2+e1**2*e2**2*e3*2.0_dp+e1**2*e2**2*e4+e1**2*e3**2*e4+e2**2*e3**2*e4+e1*e2*e3**2*e4-e1*e2**2*e3*e4*2.0_dp-e1**2*e2*e3*e4*2.0_dp)*4.0_dp)/((e1-e4)*(e2-e4))
                w(4) = (1.0_dp/(e1-e4)**2*1.0_dp/(e2-e4)**2*(E**2*(e1*e2**2*3.0_dp+e1**2*e2*3.0_dp+e3*e4**2*3.0_dp-e1*e2*e3*3.0_dp-e1*e2*e4*6.0_dp)-E**3*(e1*e2-e1*e3-e1*e4*2.0_dp-e2*e3-e2*e4*2.0_dp+e3*e4*2.0_dp+e1**2+e2**2+e4**2)-E*(e1**2*e2**2*3.0_dp-e1*e2*e4**2*3.0_dp+e1*e3*e4**2*3.0_dp+e2*e3*e4**2*3.0_dp-e1*e2*e3*e4*6.0_dp)+e1**2*e2**2*e3-e1*e2**2*e4**2-e1**2*e2*e4**2+e1**2*e2**2*e4*2.0_dp+e1**2*e3*e4**2+e2**2*e3*e4**2+e1*e2*e3*e4**2-e1*e2**2*e3*e4*2.0_dp-e1**2*e2*e3*e4*2.0_dp)*4.0_dp)/((e1-e3)*(e2-e3))
                f    = (e1*e2*(-2.4D1)+e3*e4*2.4D1+E*(e1+e2-e3-e4)*2.4D1)/((e1-e3)*(e1-e4)*(e2-e3)*(e2-e4))
            case (4)
                w(1) =  ((E-e4)**3*1.0_dp/(e1-e4)**2*(-4.0_dp))/((e2-e4)*(e3-e4))
                w(2) =  ((E-e4)**3*1.0_dp/(e2-e4)**2*(-4.0_dp))/((e1-e4)*(e3-e4))
                w(3) =  ((E-e4)**3*1.0_dp/(e3-e4)**2*(-4.0_dp))/((e1-e4)*(e2-e4))
                w(4) =   (E-e4)**2*1.0_dp/(e1-e4)**2*1.0_dp/(e2-e4)**2*1.0_dp/(e3-e4)**2*(E*(e1*e2+e1*e3+e2*e3-e4*(e1*2.0_dp+e2*2.0_dp+e3*2.0_dp)+e4**2*3.0_dp)+e4*(e1*(e2+e3)*2.0_dp+e2*e3*2.0_dp)-e4**2*(e1+e2+e3)-e1*e2*e3*3.0_dp)*4.0_dp
                f    = -(E*2.4D1-e4*2.4D1)/((e1-e4)*(e2-e4)*(e3-e4))
            case (5)
                w(1) = 0.0_dp
                w(2) = 0.0_dp
                w(3) = 0.0_dp
                w(4) = 0.0_dp
                f    = 0.0_dp
            case default
                stop 'ERROR [tetrahedron_weights]: delta, invalid bracket'
            end select
        else
            stop 'ERROR [tetrahedron_weights]: integration type /= heavi or delta'
        endif
        ! check for NaNs
        if (any(isnan(w))) stop 'ERROR [tetrahedron_weights]: NaN obtained'
        ! apply Bloch correction
        if (index(flags,'blochl').ne.0) then
            f = f*0.25_dp
            do i = 1,4
                w(1)=w(1)+0.025_dp*f*(Ec(i)-e1)
                w(2)=w(2)+0.025_dp*f*(Ec(i)-e2)
                w(3)=w(3)+0.025_dp*f*(Ec(i)-e3)
                w(4)=w(4)+0.025_dp*f*(Ec(i)-e4)
            enddo
        endif
        ! divide by four because ... anyway seems right.
        w=w*0.25_dp
        !
    end function   get_tetrahedron_weight

end module am_dos













