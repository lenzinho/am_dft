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
        real(dp) :: nelecs ! number of electrons
        real(dp) :: Ef     ! fermi energy
        integer  :: nEs
        integer  :: nprojections
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

    public :: get_Ef

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
                #:for ATTRIBUTE in ['nelecs','Ef','nEs','nprojections']
                    $:write_xml_attribute_nonallocatable(ATTRIBUTE)
                #:endfor
                ! allocatable        
                #:for ATTRIBUTE in ['E','D','iD','pD']
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
                #:for ATTRIBUTE in ['nelecs','Ef','nEs','nprojections']
                    $:read_xml_attribute_nonallocatable(ATTRIBUTE)
                #:endfor
                ! allocatable        
                #:for ATTRIBUTE in ['E','D','iD','pD']
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

    function       get_dos_quick(Ep,E,kptw,degauss,flags) result(D)
        ! flags = fermi, mp, mv
        implicit none
        !
        real(dp), intent(in) :: Ep(:)    ! probing energies
        real(dp), intent(in) :: E(:,:)   ! E(nbands,nkpts) band energies
        real(dp), intent(in) :: kptw(:)  ! kptw(nkpts) normalized kpoint weights
        real(dp), intent(in) :: degauss
        character(*), intent(in) :: flags
        real(dp), allocatable :: D(:)
        integer  :: nEs
        integer  :: nbands
        integer  :: nkpts
        integer  :: i,j,k
        real(dp) :: xp
        ! get number of probing energies
        nEs = size(Ep)
        ! get n of band
        nbands = size(E,1)
        ! get n of kpoints
        nkpts = size(kptw)
        ! allocate space for dos
        allocate(D(nEs))
        !$OMP PARALLEL PRIVATE(i,j,k) SHARED(D,E,Ep,nEs,kptw)
        !$OMP DO
        do i = 1, nEs
            ! initialize D(i)
            D(i) = 0.0_dp
            ! print progress
            call show_progress(iteration=i,maximum=nEs)
            do j = 1, nbands
            do k = 1, nkpts
                xp = ( E(j,k)-Ep(i) )/degauss
                ! get contribution
                if     (index(flags,'fermi').ne.0) then
                    D(i) = D(i) + kptw(k) * fermi_dirac_dydx(x=xp)
                elseif (index(flags,'mp').ne.0) then
                    D(i) = D(i) + kptw(k) * methfessel_paxton_dydx(x=xp,n=1)
                elseif (index(flags,'mv').ne.0) then
                    D(i) = D(i) + kptw(k) * marzari_vanderbilt_dydx(x=xp)
                else
                    stop 'ERROR [get_dos_quick]: flag not recognized'
                endif
            enddo
            enddo
        enddo
        !$OMP END DO
        !$OMP END PARALLEL
        D = D / degauss
    end function   get_dos_quick

    ! DOS properties

    function       get_Ef(E,tet,tetw,nelecs) result(Ef)
        ! adapted from qe
        implicit none
        !
        real(dp), intent(in) :: E(:,:)
        integer , intent(in) :: tet(:,:)
        real(dp), intent(in) :: tetw(:)
        real(dp), intent(in) :: nelecs
        real(dp) :: Ef
        integer  :: i
        integer  :: ntets
        integer  :: nbands
        integer  :: nkpts
        real(dp) :: Emin    ! absolute lowest band energy
        real(dp) :: Emax    ! absolute highest band energy
        real(dp) :: elw     ! lower bound for EF
        real(dp) :: eup     ! upper bound for EF
        real(dp) :: idos_up ! number of states for Ef=eup
        real(dp) :: idos_lw ! number of states for Ef=elw
        real(dp) :: idos_ep ! number of states for Ef=(elw+eup)/2
        real(dp) :: try_Ef  ! new Ef to try
        real(dp) :: crit    ! stopping criterion
        ! get bands, kpts, tetrahedra
        nbands = size(E,1)
        nkpts  = size(E,2)
        ntets  = sizE(tet,2)
        ! get lower and upper bracketing bounds
        Emin = minval(E)
        Emax = maxval(E)
        elw = Emin
        eup = Emax
        ! initailize bracketing 
        idos_up = get_idos_at_ep(Ep=eup,E=E,tet=tet,tetw=tetw)
        idos_lw = get_idos_at_ep(Ep=elw,E=E,tet=tet,tetw=tetw)
        crit = 1.0d+10
        ! search for Ef by gradualling reducing bracketted region
        search : do i = 1, 1000 ! hard-coded 1000 maximum iterations
            try_Ef = (eup+elw) / 2.0_dp
            idos_ep = get_idos_at_ep(Ep=try_Ef,E=E,tet=tet,tetw=tetw)
            if (abs(idos_ep-nelecs).lt.crit) then
                crit = abs(idos_ep-nelecs)
                Ef = try_Ef
            endif
            ! converged
            if (abs(idos_ep-nelecs).lt.tiny) then
                exit search
            elseif((idos_ep-nelecs).lt.-tiny) then
                elw = try_Ef
            else
                eup = try_Ef
            endif
        enddo search
        ! check convergence
        if (abs(idos_ep-nelecs).gt.tiny) stop 'ERROR [XX] failed to converge on Ef.'
        ! check that Ef < highest band
        if (Ef.gt.Emax) stop 'ERROR [XX]: Ef is above the highest band energy'
        !
    end function   get_Ef

    ! tetrahedra methods

    function       get_dos_vs_Ep(Ep,E,tet,tetw) result(D)
        !
        implicit none
        !
        real(dp), intent(in) :: Ep(:) ! probing energies
        real(dp), intent(in) :: E(:,:) ! E(nbands,nkpts) band energies
        integer , intent(in) :: tet(:,:) ! tet(:,ntets) tetrahedra conenctivity
        real(dp), intent(in) :: tetw(:) ! tetw(ntets) weight of each tetrahedron
        real(dp), allocatable :: D(:)
        integer :: nEs
        integer :: i
        ! get number of probing energies
        nEs = size(Ep)
        ! allocate space for dos
        allocate(D(nEs))
        !$OMP PARALLEL PRIVATE(i) SHARED(D,nEs,E,Ep,tet,tetw)
        !$OMP DO
        do i = 1, nEs
            ! print progress
            call show_progress(iteration=i,maximum=nEs)
            ! get dos
            D(i) = get_dos_at_ep(Ep=Ep(i),E=E,tet=tet,tetw=tetw)
        enddo
        !$OMP END DO
        !$OMP END PARALLEL
    end function   get_dos_vs_Ep

    function       get_idos_vs_Ep(Ep,E,tet,tetw) result(iD)
        !
        implicit none
        !
        real(dp), intent(in) :: Ep(:) ! probing energies
        real(dp), intent(in) :: E(:,:) ! E(nbands,nkpts) band energies
        integer , intent(in) :: tet(:,:) ! tet(:,ntets) tetrahedra conenctivity
        real(dp), intent(in) :: tetw(:) ! tetw(ntets) weight of each tetrahedron
        real(dp), allocatable :: iD(:)
        integer :: nEs
        integer :: i
        ! get number of probing energies
        nEs = size(Ep)
        ! allocate space for dos
        allocate(iD(nEs))
        !$OMP PARALLEL PRIVATE(i) SHARED(iD,nEs,E,Ep,tet,tetw)
        !$OMP DO
        do i = 1, nEs
            ! print progress
            call show_progress(iteration=i,maximum=nEs)
            ! get dos
            iD(i) = get_idos_at_ep(Ep=Ep(i),E=E,tet=tet,tetw=tetw)
        enddo
        !$OMP END DO
        !$OMP END PARALLEL
    end function   get_idos_vs_Ep

    function       get_pdos_vs_Ep(Ep,E,tet,tetw,weight) result(pD)
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
        real(dp) :: wc(4) ! corner tetrahedron weights
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
        !$OMP PARALLEL PRIVATE(i,j,k,m,wc) SHARED(pD,E,Ep,nEs,tetw,weight)
        !$OMP DO
        do i = 1, nEs
            ! initialize D(i)
            pD(:,i) = 0.0_dp
            ! print progress
            call show_progress(iteration=i,maximum=nEs)
            do j = 1, nbands
            do k = 1, ntets
                ! get tetrahedron corner weights
                wc = get_delta_wc(Ep=Ep(i),Ec=E(j,tet(:,k)))
                ! loop over projections, increment pDOS with contributions from band j in tetrahedron k
                do m = 1, nprojections
                    pD(m,i) = pD(m,i) + tetw(k) * sum(wc * weight(m,j,tet(:,k)) )
                enddo
            enddo
            enddo
        enddo
        !$OMP END DO
        !$OMP END PARALLEL
    end function   get_pdos_vs_Ep

    pure function  get_theta_wc(Ep,Ec,ntets) result(wc)
        ! get linear tetrahedron corner weights for delta function with Blochl corrections
        implicit none
        !
        real(dp), intent(in) :: Ep    ! probing energy
        real(dp), intent(in) :: Ec(4) ! corner energies
        integer , intent(in) :: ntets
        real(dp) :: wc(4)
        real(dp) :: dosEp
        real(dp) :: etot
        real(dp) :: c1,c2,c3,c4
        real(dp) :: e1,e2,e3,e4
        integer  :: k1,k2,k3,k4
        integer  :: inds(4)
        ! rank corner weights in increasing order
        inds = fourrank(Ec)
        ! sort weights
        e1 = Ec(inds(1))
        e2 = Ec(inds(2))
        e3 = Ec(inds(3))
        e4 = Ec(inds(4))
        ! k1-k4 are the irreducible k-points corresponding to e1-e4
        k1 = inds(1)
        k2 = inds(2)
        k3 = inds(3)
        k4 = inds(4)
        ! calculate weights wc
        if     (Ep.le.e1) then
            ! Eq B1 from Blochl's PhysRevB.49.16223
            wc(k1)=0.0_dp
            wc(k2)=0.0_dp
            wc(k3)=0.0_dp
            wc(k4)=0.0_dp
        elseif (Ep.le.e2) then
            ! Eq B6 from Blochl's PhysRevB.49.16223
            c4=0.25_dp/ntets*(Ep-e1)**3/(e2-e1)/(e3-e1)/(e4-e1)
            ! Eq C2 from Blochl's PhysRevB.49.16223
            dosEp=3.0_dp/ntets*(Ep-e1)**2/(e2-e1)/(e3-e1)/(e4-e1)
            ! shortcut for Blochl corrections (Eq.22 PhysRevB.49.16223)
            etot=e1+e2+e3+e4
            ! Eq B2 from Blochl's PhysRevB.49.16223
            wc(k1)=c4*(4.0_dp-(Ep-e1)*(1.0_dp/(e2-e1)+1.0_dp/(e3-e1)+1.0_dp/(e4-e1)))+dosEp*(etot-4.0_dp*e1)*0.025_dp
            ! Eq B3 from Blochl's PhysRevB.49.16223
            wc(k2)=c4*(Ep-e1)/(e2-e1)+dosEp*(etot-4.0_dp*e2)*0.025_dp
            ! Eq B4 from Blochl's PhysRevB.49.16223
            wc(k3)=c4*(Ep-e1)/(e3-e1)+dosEp*(etot-4.0_dp*e3)*0.025_dp
            ! Eq B5 from Blochl's PhysRevB.49.16223
            wc(k4)=c4*(Ep-e1)/(e4-e1)+dosEp*(etot-4.0_dp*e4)*0.025_dp
        elseif (Ep.le.e3) then
            ! Eq B11 from Blochl's PhysRevB.49.16223
            c1=0.25_dp/ntets*(Ep-e1)**2/(e4-e1)/(e3-e1)
            ! Eq B12 from Blochl's PhysRevB.49.16223
            c2=0.25_dp/ntets*(Ep-e1)*(Ep-e2)*(e3-Ep)/(e4-e1)/(e3-e2)/(e3-e1)
            ! Eq B13 from Blochl's PhysRevB.49.16223
            c3=0.25_dp/ntets*(Ep-e2)**2*(e4-Ep)/(e4-e2)/(e3-e2)/(e4-e1)
            ! Eq C3 from Blochl's PhysRevB.49.16223
            dosEp=1.0_dp/ntets/(e3-e1)/(e4-e1)*(3.0_dp*(e2-e1)+6.0_dp*(Ep-e2)-3.0_dp*(e3-e1+e4-e2)*(Ep-e2)**2/(e3-e2)/(e4-e2))
            ! shortcut for Blochl corrections (Eq.22 PhysRevB.49.16223)
            etot=e1+e2+e3+e4
            ! Eq B7 from Blochl's PhysRevB.49.16223
            wc(k1)=c1+(c1+c2)*(e3-Ep)/(e3-e1)+(c1+c2+c3)*(e4-Ep)/(e4-e1)+dosEp*(etot-4.0_dp*e1)*0.025_dp
            ! Eq B8 from Blochl's PhysRevB.49.16223
            wc(k2)=c1+c2+c3+(c2+c3)*(e3-Ep)/(e3-e2)+c3*(e4-Ep)/(e4-e2)+dosEp*(etot-4.0_dp*e2)*0.025_dp
            ! Eq B9 from Blochl's PhysRevB.49.16223
            wc(k3)=(c1+c2)*(Ep-e1)/(e3-e1)+(c2+c3)*(Ep-e2)/(e3-e2)+dosEp*(etot-4.0_dp*e3)*0.025_dp
            ! Eq B10 from Blochl's PhysRevB.49.16223
            wc(k4)=(c1+c2+c3)*(Ep-e1)/(e4-e1)+c3*(Ep-e2)/(e4-e2)+dosEp*(etot-4.0_dp*e4)*0.025_dp
        elseif (Ep.le.e4) then
            ! Eq B18 from Blochl's PhysRevB.49.16223
            c4=0.25_dp/ntets*(e4-Ep)**3/(e4-e1)/(e4-e2)/(e4-e3)
            ! Eq C4 from Blochl's PhysRevB.49.16223
            dosEp=3.0_dp/ntets*(e4-Ep)**2/(e4-e1)/(e4-e2)/(e4-e3)
            ! shortcut for Blochl corrections (Eq.22 PhysRevB.49.16223)
            etot=e1+e2+e3+e4
            ! Eq B14 from Blochl's PhysRevB.49.16223
            wc(k1)=0.25_dp/ntets-c4*(e4-Ep)/(e4-e1)+dosEp*(etot-4.0_dp*e1)*0.025_dp
            ! Eq B15 from Blochl's PhysRevB.49.16223
            wc(k2)=0.25_dp/ntets-c4*(e4-Ep)/(e4-e2)+dosEp*(etot-4.0_dp*e2)*0.025_dp
            ! Eq B16 from Blochl's PhysRevB.49.16223
            wc(k3)=0.25_dp/ntets-c4*(e4-Ep)/(e4-e3)+dosEp*(etot-4.0_dp*e3)*0.025_dp
            ! Eq B17 from Blochl's PhysRevB.49.16223
            wc(k4)=0.25_dp/ntets-c4*(4.0_dp-(e4-Ep)*(1.0_dp/(e4-e1)+1.0_dp/(e4-e2)+1.0_dp/(e4-e3)))+dosEp*(etot-4.0_dp*e4)*0.025_dp
        elseif (Ep.ge.e4) then
            ! Eq B19 from Blochl's PhysRevB.49.16223
            wc(k1)=0.25_dp/ntets
            wc(k2)=0.25_dp/ntets
            wc(k3)=0.25_dp/ntets
            wc(k4)=0.25_dp/ntets
        endif
        ! 
        wc = wc * ntets
        !
    end function   get_theta_wc

    pure function  get_delta_wc(Ep,Ec) result(wc)
        ! get linear tetrahedron corner weights for delta function with Blochl corrections
        implicit none
        !
        real(dp), intent(in) :: Ep    ! probing energy
        real(dp), intent(in) :: Ec(4) ! corner energies
        real(dp) :: wc(4)
        real(dp) :: dosEp
        real(dp) :: e1,e2,e3,e4
        integer  :: k1,k2,k3,k4
        real(dp) :: o13,f12,f13,f14,f21,f23,f31,f32,f34,f24,f41,f42
        integer  :: inds(4)
        ! initialize shortcut
        o13  = 1.0_dp/3.0_dp
        ! rank corner weights in increasing order
        inds = fourrank(Ec)
        ! sort weights
        e1 = Ec(inds(1))
        e2 = Ec(inds(2))
        e3 = Ec(inds(3))
        e4 = Ec(inds(4))
        ! k1-k4 are the irreducible k-points corresponding to e1-e4
        k1 = inds(1)
        k2 = inds(2)
        k3 = inds(3)
        k4 = inds(4)
        ! calculate weights wc
        if     (Ep.le.e1) then
            ! Eq B1 from Blochl's PhysRevB.49.16223
            wc(k1)=0.0_dp
            wc(k2)=0.0_dp
            wc(k3)=0.0_dp
            wc(k4)=0.0_dp
        elseif (Ep.le.e2) then
            f12 = (Ep-e2)/(e1-e2)
            f21 = 1.0_dp - f12
            f13 = (Ep-e3)/(e1-e3)
            f31 = 1.0_dp - f13
            f14 = (Ep-e4)/(e1-e4)
            f41 = 1.0_dp - f14
            dosEp  = 3.0_dp * f21 * f31 * f41 / (Ep-e1)
            wc(k1) = o13 * (f12 + f13 + f14)
            wc(k2) = o13 * f21
            wc(k3) = o13 * f31
            wc(k4) = o13 * f41
            wc = wc * dosEp
        elseif (Ep.le.e3) then
            f13 = (Ep-e3)/(e1-e3)
            f31 = 1.0_dp - f13
            f14 = (Ep-e4)/(e1-e4)
            f41 = 1.0_dp-f14
            f23 = (Ep-e3)/(e2-e3)
            f32 = 1.0_dp - f23
            f24 = (Ep-e4)/(e2-e4)
            f42 = 1.0_dp - f24
            dosEp  = 3.0_dp * (f23*f31 + f32*f24)
            wc(k1) = f14 * o13 + f13*f31*f23 / dosEp
            wc(k2) = f23 * o13 + f24*f24*f32 / dosEp
            wc(k3) = f32 * o13 + f31*f31*f23 / dosEp
            wc(k4) = f41 * o13 + f42*f24*f32 / dosEp
            dosEp  = dosEp / (e4-e1)
            wc = wc * dosEp
        elseif (Ep.le.e4) then
            f14 = (Ep-e4)/(e1-e4)
            f24 = (Ep-e4)/(e2-e4)
            f34 = (Ep-e4)/(e3-e4)
            dosEp  = 3.0_dp * f14 * f24 * f34 / (e4-Ep)
            wc(k1) = f14 * o13
            wc(k2) = f24 * o13
            wc(k3) = f34 * o13
            wc(k4) = (3.0_dp - f14 - f24 - f34 ) * o13
            wc = wc * dosEp
        elseif (Ep.gt.e4) then
            wc(k1)=0.0_dp
            wc(k2)=0.0_dp
            wc(k3)=0.0_dp
            wc(k4)=0.0_dp
       endif
       !
    end function   get_delta_wc

    pure function  get_idos_at_Ep(Ep,E,tet,tetw) result(idos)
        ! return integrated DOS at Ep
        implicit none
        !
        real(dp), intent(in) :: Ep
        real(dp), intent(in) :: E(:,:)
        integer , intent(in) :: tet(:,:)
        real(dp), intent(in) :: tetw(:) ! weight of irreducible tetrahedron
        real(dp) :: idos
        real(dp) :: Ec(4)
        real(dp) :: e1,e2,e3,e4
        integer  :: i,j
        integer  :: nbands
        integer  :: ntets
        ! get bands, kpts, tetrahedra
        nbands = size(E,1)
        ntets  = sizE(tet,2)
        !
        idos = 0.0_dp
        do j = 1, ntets
            do i = 1, nbands
                ! get energies at the vertetxes of the j-th tethedron in ascending order
                Ec(1:4) = foursort(E(i,tet(:,j)))
                ! map to 
                e1 = Ec(1)
                e2 = Ec(2)
                e3 = Ec(3)
                e4 = Ec(4)
                ! calculate sum over k of the integrated charge
                if     (Ep.le.e1) then
                    ! Eq A1 from Blochl's PhysRevB.49.16223
                    ! do nothing
                elseif (Ep.le.e2) then
                    ! Eq A2 from Blochl's PhysRevB.49.16223
                    idos = idos + tetw(j)*(Ep-e1)**3/(e2-e1)/(e3-e1)/(e4-e1)
                elseif (Ep.le.e3) then
                    ! Eq A3 from Blochl's PhysRevB.49.16223
                    idos = idos + tetw(j)/(e3-e1)/(e4-e1)*((e2-e1)**2+3.0_dp*(e2-e1)*(Ep-e2)+3.0_dp*(Ep-e2)**2-(e3-e1+e4-e2)/(e3-e2)/(e4-e2)*(Ep-e2)**3)
                elseif (Ep.le.e4) then
                    ! Eq A4 from Blochl's PhysRevB.49.16223
                    idos = idos + tetw(j)*(1.0_dp-(e4-Ep)**3/(e4-e1)/(e4-e2)/(e4-e3))
                elseif (Ep.ge.e4) then
                    ! Eq A5 from Blochl's PhysRevB.49.16223
                    idos = idos + tetw(j) 
                endif
            enddo
        enddo
    end function   get_idos_at_Ep

    pure function  get_dos_at_Ep(Ep,E,tet,tetw) result(dosEp)
        ! returns DOS at Ep
        implicit none
        !
        real(dp), intent(in) :: Ep
        real(dp), intent(in) :: E(:,:)
        integer , intent(in) :: tet(:,:)
        real(dp), intent(in) :: tetw(:) ! weight of irreducible tetrahedron
        real(dp) :: dosEp
        real(dp) :: Ec(4)
        real(dp) :: e1,e2,e3,e4
        integer :: i,j
        integer :: nbands
        integer :: ntets
        ! get bands, kpts, tetrahedra
        nbands = size(E,1)
        ntets  = sizE(tet,2)
        !
        dosEp = 0.0_dp
        do j = 1, ntets
            do i = 1, nbands
                ! get energies at the vertetxes of the j-th tethedron in ascending order
                Ec(1:4) = foursort(E(i,tet(:,j)))
                ! map to 
                e1 = Ec(1)
                e2 = Ec(2)
                e3 = Ec(3)
                e4 = Ec(4)
                ! calculate sum over k of the integrated charge
                if     (Ep.le.e1) then
                    ! Eq C1 from Blochl's PhysRevB.49.16223
                    ! do nothing
                elseif (Ep.le.e2) then
                    ! Eq C2 from Blochl's PhysRevB.49.16223
                    dosEp = dosEp + tetw(j)*3.0_dp*(Ep-e1)**2/(e2-e1)/(e3-e1)/(e4-e1)
                elseif (Ep.le.e3) then
                    ! Eq C3 from Blochl's PhysRevB.49.16223
                    dosEp = dosEp + tetw(j)/(e3-e1)/(e4-e1)*(3.0_dp*(e2-e1)+6.0_dp*(Ep-e2)-3.0_dp*(e3-e1+e4-e2)/(e3-e2)/(e4-e2)*(Ep-e2)**2)
                elseif (Ep.le.e4) then
                    ! Eq C4 from Blochl's PhysRevB.49.16223
                    dosEp = dosEp + tetw(j)*(3.0_dp*(e4-Ep)**2/(e4-e1)/(e4-e2)/(e4-e3))
                elseif (Ep.ge.e4) then
                    ! Eq C1 from Blochl's PhysRevB.49.16223
                    ! do nothing
                endif
            enddo
        enddo
        !
    end function   get_dos_at_Ep

    ! heaviside (theta) integration functions

    pure function methfessel_paxton(x,n) result(y)
        ! theta function : PRB 40, 3616 (1989).
        implicit none
        !
        real(dp), intent(in) :: x
        integer , intent(in) :: n
        real(dp) :: y
        
        real(dp) :: a, hp, arg, hd
        integer :: i, ni
        real(dp), parameter :: maxarg = 200.0_dp
        
        ! Methfessel-Paxton
        y = gauss_freq(x * sqrt(2.0_dp) )
        if (n.eq.0) return
        hd = 0.0_dp
        arg = min(maxarg, x**2)
        hp = exp(- arg)
        ni = 0
        a = 1.0_dp / sqrt (pi)
        do i=1,n
            hd = 2.0_dp*x*hp-2.0_dp*real(ni,dp)*hd
            ni = ni+1
            a = -a/(real(i,dp)*4.0_dp)
            y = y-a*hd
            hp = 2.0_dp*x*hd-2.0_dp*real(ni,dp)*hp
            ni = ni+1
      enddo
    end function  methfessel_paxton

    pure function fermi_dirac(x) result(y)
        !
        implicit none
        !
        real(dp), intent(in) :: x
        real(dp) :: y
        real(dp) :: maxarg
        maxarg = 200.0_dp
        ! Fermi-Dirac smearing
        if (x.lt.-maxarg) then
            y = 0.0_dp
        elseif (x.gt.maxarg) then
            y = 1.0_dp
        else
            y = 1.0_dp/(1.0_dp + exp(-x))
        endif
    end function  fermi_dirac

    pure function marzari_vanderbilt(x) result(y)
        ! theta function: PRL 82, 3296 (1999)
        ! 1/2*erf(x-1/sqrt(2)) + 1/sqrt(2*pi)*exp(-(x-1/sqrt(2))**2) + 1/2
        implicit none
        !
        real(dp), intent(in) :: x
        real(dp) :: y
        real(dp) :: arg, xp
        real(dp) :: maxarg
        maxarg = 200.0_dp
        ! Cold smearing
         xp = x - 1.0_dp / sqrt (2.0_dp)
         arg = min (maxarg, xp**2)
         y = 0.5d0 * erf (xp) + 1.0_dp / sqrt (2.0_dp * pi) * exp (- arg) + 0.5d0
    end function  marzari_vanderbilt

    pure function erf(x)
      !---------------------------------------------------------------------
      !
      !     Error function - computed from the rational approximations of
      !     W. J. Cody, Math. Comp. 22 (1969), pages 631-637.
      !
      !     for abs(x) le 0.47 erf is calculated directly
      !     for abs(x) gt 0.47 erf is calculated via erf(x)=1-erfc(x)
      !
      implicit none  
      real(dp), intent(in) :: x
      real(dp) :: x2, p1(4), q1(4)
      real(dp) :: erf
      !
      p1 = [2.426679552305318E2_dp, 2.197926161829415E1_dp, 6.996383488619136_dp,  -3.560984370181538E-2_dp]
      q1 = [2.150588758698612E2_dp, 9.116490540451490E1_dp, 1.508279763040779E1_dp, 1.000000000000000_dp]
      !
      if (abs(x).gt.6.0_dp) then  
         !  erf(6)=1-10^(-17) cannot be distinguished from 1
         erf = sign (1.0_dp, x)  
      else  
         if (abs(x).le.0.47_dp) then  
            x2 = x**2  
            erf = x*(p1(1)+x2*(p1(2)+x2*(p1(3)+x2*p1(4))))/(q1(1)+x2*(q1(2)+x2*(q1(3)+x2*q1(4))))
         else  
            erf = 1.0_dp - erfc(x)  
         endif
      endif
    end function  erf

    pure function erfc(x)
      !     erfc(x) = 1-erf(x)  - See comments in erf
      implicit none  
      real(dp),intent(in) :: x
      real(dp)            :: erfc
      real(dp) :: ax, x2, xm2, p2(8), q2(8), p3(5), q3(5), pim1
      !
      p2 = [ 3.004592610201616E2_dp,  4.519189537118719E2_dp,  3.393208167343437E2_dp,  1.529892850469404E2_dp,  4.316222722205674E1_dp,  7.211758250883094_dp,    5.641955174789740E-1_dp,-1.368648573827167E-7_dp]
      q2 = [ 3.004592609569833E2_dp,  7.909509253278980E2_dp,  9.313540948506096E2_dp,  6.389802644656312E2_dp,  2.775854447439876E2_dp,  7.700015293522947E1_dp,  1.278272731962942E1_dp,  1.000000000000000_dp]
      p3 = [-2.996107077035422E-3_dp,-4.947309106232507E-2_dp,  -2.269565935396869E-1_dp,-2.786613086096478E-1_dp,  -2.231924597341847E-2_dp]
      q3 = [ 1.062092305284679E-2_dp, 1.913089261078298E-1_dp,  1.051675107067932_dp,    1.987332018171353_dp,     1.000000000000000_dp]
      pim1 = 0.56418958354775629_dp  ! sqrt(1/pi)
      ax = abs(x)  
      if (ax > 26.0_dp) then  
         !  erfc(26.0)=10^(-296); erfc(9.0)=10^(-37);
         erfc = 0.0_dp  
      elseif (ax.gt.4.0_dp) then
         xm2=(1.0_dp/ax)**2
         erfc=(1.0_dp/ax)*exp(-2.0_dp)*(pim1+xm2*(p3(1)+xm2*(p3(2)+xm2*(p3(3)+xm2*(p3(4)+xm2*p3(5)))))/(q3(1)+xm2*(q3(2)+xm2*(q3(3)+xm2*(q3(4)+xm2*q3(5))))))
      elseif(ax.gt.0.47_dp)then
         erfc=exp(-x**2)*(p2(1)+ax*(p2(2)+ax*(p2(3)+ax*(p2(4)+ax*(p2(5)+ax*(p2(6)+ax*(p2(7)+ax*p2(8))))))))/(q2(1)+ax*(q2(2)+ax*(q2(3)+ax*(q2(4)+ax*(q2(5)+ax*(q2(6)+ax*(q2(7)+ax*q2(8))))))))
      else
         erfc=1.0_dp-erf(ax)
      endif
      ! erf(-x)=-erf(x)  =>  erfc(-x) = 2-erfc(x)
      if (x < 0.0_dp) erfc = 2.0_dp - erfc 
    end function  erfc

    pure function gauss_freq(x)
        !     gauss_freq(x) = (1+erf(x/sqrt(2)))/2 = erfc(-x/sqrt(2))/2
        implicit none
        !
        real(dp),intent(in) :: x
        real(dp)            :: gauss_freq
        !
        gauss_freq = 0.5_dp * erfc(-x*0.7071067811865475_dp)
        !
    end function  gauss_freq

    ! derivatives of theta functions (= delta functions)

    pure function fermi_dirac_dydx(x) result(y)
        ! derivative of Fermi-Dirac function: 0.5/(1.0+cosh(x))
        implicit none
        real(dp), intent(in) :: x
        real(dp) :: y
        !
        if (abs(x).le.36.0_dp) then
            y = 1.0_dp/(2.0_dp+exp(-x)+exp(+x))
            ! in order to avoid problems for large values of x in the e
        else
            y = 0.0_dp
        endif
    end function  fermi_dirac_dydx

    pure function marzari_vanderbilt_dydx(x) result(y)
        ! 1/sqrt(pi)*exp(-(x-1/sqrt(2))**2)*(2-sqrt(2)*x)
        implicit none
        real(dp), intent(in) :: x
        real(dp) :: y
        real(dp) :: arg
        real(dp) :: sqrtpm1
        !
        sqrtpm1 = 0.564189583547756_dp ! 1/sqrt(pi)
        arg = min(200.0_dp,(x-1.0_dp/sqrt(2.0_dp))**2)
        y = sqrtpm1*exp(-arg)*(2.0_dp-sqrt(2.0_dp)*x)
    end function  marzari_vanderbilt_dydx

    pure function methfessel_paxton_dydx(x, n) result(y)
        ! derivative of the corresponding Methfessel-Paxton wgauss
        implicit none
        real(dp), intent(in) :: x
        integer , intent(in) :: n
        real(dp) :: y
        real(dp) :: a, arg, hp, hd
        integer :: i, ni
        real(dp) :: sqrtpm1
        !
        sqrtpm1 = 0.564189583547756_dp ! 1/sqrt(pi)
        arg = min(200.0_dp, x**2)
        y = exp(-arg)*sqrtpm1
        if (n.eq.0) return
        hd = 0.0_dp
        hp = exp(-arg)
        ni = 0
        a  = sqrtpm1
        do i = 1, n
            hd = 2.0_dp*x*hp-2.0_dp*real(ni,dp)*hd
            ni = ni+1
            a  = -a/(real(i,dp)*4.0_dp)
            hp = 2.0_dp*x*hd-2.0_dp*real(ni,dp)*hp
            ni = ni+1
            y  = y+a*hp
        enddo
    end function  methfessel_paxton_dydx


end module am_dos













