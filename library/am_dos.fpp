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
                #:for ATTRIBUTE in ['nelecs','Ef','nEs','nprojections']
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
        idos_up = get_idos_ep(Ep=eup,E=E,tet=tet,tetw=tetw)
        idos_lw = get_idos_ep(Ep=elw,E=E,tet=tet,tetw=tetw)
        crit = 1.0d+10
        ! search for Ef by gradualling reducing bracketted region
        search : do i = 1, 1000 ! hard-coded 1000 maximum iterations
            try_Ef = (eup+elw) / 2.0_dp
            idos_ep = get_idos_ep(Ep=try_Ef,E=E,tet=tet,tetw=tetw)
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
        real(dp) :: wc(4) ! weights on tetrahedron tet
        integer  :: i,j,k
        ! get number of probing energies
        nEs = size(Ep)
        ! get n of band
        nbands = size(E,1)
        ! number of tetrahedra
        ntets = size(tet,2)
        ! allocate space for dos
        allocate(D(nEs))
        !$OMP PARALLEL PRIVATE(i,j,k,wc) SHARED(D,E,Ep,nEs,tetw)
        !$OMP DO
        do i = 1, nEs
            ! initialize D(i)
            D(i) = 0.0_dp
            ! print progress
            call show_progress(iteration=i,maximum=nEs)
            do j = 1, nbands
            do k = 1, ntets
                ! get tetrahedron corner weights
                ! w = get_tetrahedron_weight(E=Ep(i),Ec=E(j,tet(:,k)),flags='delta,blochl')
                wc = get_wc(Ep=Ep(i),Ec=E(j,tet(:,k)),ntets=ntets)
                ! increment DOS with contributions from band j in tetrahedron k
                D(i) = D(i) + tetw(k) * sum(wc)
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
                wc = get_wc(Ep=Ep(i),Ec=E(j,tet(:,k)),ntets=ntets)
                ! loop over projections, increment pDOS with contributions from band j in tetrahedron k
                do m = 1, nprojections
                    pD(m,i) = pD(m,i) + tetw(k) * sum(wc * weight(m,j,tet(:,k)) )
                enddo
            enddo
            enddo
        enddo
        !$OMP END DO
        !$OMP END PARALLEL
    end function   get_pdos_tetra

    function       get_wc(Ep,Ec,ntets) result(wc)
        ! get linear tetrahedron corner weights for delta function with Blochl corrections
        implicit none
        !
        real(dp), intent(in) :: Ep    ! probing energy
        real(dp), intent(in) :: Ec(4) ! corner energies
        real(dp) :: wc(4)
        integer  :: ntets
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
    end function   get_wc

        ! function       get_tetrahedron_weight(E,Ec,flags) result(w)
        !     !
        !     ! flags = heavi/delta, blochl
        !     !
        !     use am_rank_and_sort
        !     !
        !     implicit none
        !     !
        !     real(dp)    , intent(in) :: Ec(4)
        !     character(*), intent(in) :: flags
        !     real(dp) :: w(4)
        !     real(dp) :: E,e1,e2,e3,e4,f
        !     real(dp) :: Ei(4)
        !     integer :: bracket
        !     integer :: i
        !     ! get sorted values
        !     Ei = foursort(Ec)
        !     e1 = Ei(1)
        !     e2 = Ei(2)
        !     e3 = Ei(3)
        !     e4 = Ei(4)
        !     !
        !     if (e1.gt.e2) stop 'ERROR [tetrahedron_weights]: sorting failed'
        !     if (e2.gt.e3) stop 'ERROR [tetrahedron_weights]: sorting failed'
        !     if (e3.gt.e4) stop 'ERROR [tetrahedron_weights]: sorting failed'
        !     !
        !     bracket = 0
        !     if     (E.lt.e1) then; bracket = 1
        !     elseif (E.lt.e2) then; bracket = 2
        !     elseif (E.lt.e3) then; bracket = 3
        !     elseif (E.lt.e4) then; bracket = 4
        !     elseif (E.gt.e4) then; bracket = 5
        !     else
        !         stop 'ERROR [tetrahedron_weights]: bracketing failed'
        !     endif
        !     !
        !     if (index(flags,'heavi').ne.0) then
        !         select case(bracket)
        !         case (1)
        !             w(1) = 0.0_dp
        !             w(2) = 0.0_dp
        !             w(3) = 0.0_dp
        !             w(4) = 0.0_dp
        !             f    = 0.0_dp
        !         case (2)
        !             w(1) = -((E-e1)**3*((E-e1)*(1.0_dp/(e1-e2)+1.0_dp/(e1-e3)+1.0_dp/(e1-e4))+4.0_dp))/((e1-e2)*(e1-e3)*(e1-e4))
        !             w(2) =  ((E-e1)**4*1.0_dp/(e1-e2)**2)/((e1-e3)*(e1-e4))
        !             w(3) =  ((E-e1)**4*1.0_dp/(e1-e3)**2)/((e1-e2)*(e1-e4))
        !             w(4) =  ((E-e1)**4*1.0_dp/(e1-e4)**2)/((e1-e2)*(e1-e3))
        !             f    =  ((E-e1)**2*(-1.2D1))/((e1-e2)*(e1-e3)*(e1-e4))
        !         case (3)
        !             w(1) =  (E-e1)**2/((e1-e3)*(e1-e4))+(((E-e1)**2/((e1-e3)*(e1-e4))+((E-e1)*(E-e2)*(E-e3))/((e1-e3)*(e1-e4)*(e2-e3)))*(E-e3))/(e1-e3)+((E-e4)*((E-e1)**2/((e1-e3)*(e1-e4))+((E-e2)**2*(E-e4))/((e1-e4)*(e2-e3)*(e2-e4))+((E-e1)*(E-e2)*(E-e3))/((e1-e3)*(e1-e4)*(e2-e3))))/(e1-e4)
        !             w(2) =  ((((E-e2)**2*(E-e4))/((e1-e4)*(e2-e3)*(e2-e4))+((E-e1)*(E-e2)*(E-e3))/((e1-e3)*(e1-e4)*(e2-e3)))*(E-e3))/(e2-e3)+(E-e1)**2/((e1-e3)*(e1-e4))+((E-e2)**2*(E-e4))/((e1-e4)*(e2-e3)*(e2-e4))+((E-e2)**2*(E-e4)**2*1.0_dp/(e2-e4)**2)/((e1-e4)*(e2-e3))+((E-e1)*(E-e2)*(E-e3))/((e1-e3)*(e1-e4)*(e2-e3))
        !             w(3) = -((((E-e2)**2*(E-e4))/((e1-e4)*(e2-e3)*(e2-e4))+((E-e1)*(E-e2)*(E-e3))/((e1-e3)*(e1-e4)*(e2-e3)))*(E-e2))/(e2-e3)-(((E-e1)**2/((e1-e3)*(e1-e4))+((E-e1)*(E-e2)*(E-e3))/((e1-e3)*(e1-e4)*(e2-e3)))*(E-e1))/(e1-e3)
        !             w(4) = -((E-e1)*((E-e1)**2/((e1-e3)*(e1-e4))+((E-e2)**2*(E-e4))/((e1-e4)*(e2-e3)*(e2-e4))+((E-e1)*(E-e2)*(E-e3))/((e1-e3)*(e1-e4)*(e2-e3))))/(e1-e4)-((E-e2)**3*(E-e4)*1.0_dp/(e2-e4)**2)/((e1-e4)*(e2-e3))
        !             f    =  (E*2.4D1-e1*1.2D1-e2*1.2D1+((E-e2)**2*(e1*3.0_dp+e2*3.0_dp-e3*3.0_dp-e4*3.0_dp)*4.0_dp)/((e2-e3)*(e2-e4)))/((e1-e3)*(e1-e4))
        !         case (4)
        !             w(1) = -((E-e4)**4*1.0_dp/(e1-e4)**2)/((e2-e4)*(e3-e4))+1.0_dp
        !             w(2) = -((E-e4)**4*1.0_dp/(e2-e4)**2)/((e1-e4)*(e3-e4))+1.0_dp
        !             w(3) = -((E-e4)**4*1.0_dp/(e3-e4)**2)/((e1-e4)*(e2-e4))+1.0_dp
        !             w(4) =  ((E-e4)**3*((E-e4)*(1.0_dp/(e1-e4)+1.0_dp/(e2-e4)+1.0_dp/(e3-e4))-4.0_dp))/((e1-e4)*(e2-e4)*(e3-e4))+1.0_dp
        !             f    =  ((E-e4)**2*(-1.2D1))/((e1-e4)*(e2-e4)*(e3-e4))
        !         case (5)
        !             w(1) = 1.0_dp
        !             w(2) = 1.0_dp
        !             w(3) = 1.0_dp
        !             w(4) = 1.0_dp
        !             f    = 0.0_dp
        !         case default
        !             stop 'ERROR [tetrahedron_weights]: heaviside, invalid bracket'
        !         end select
        !     elseif (index(flags,'delta').ne.0) then
        !         select case(bracket)
        !         case (1)
        !             w(1) = 0.0_dp
        !             w(2) = 0.0_dp
        !             w(3) = 0.0_dp
        !             w(4) = 0.0_dp
        !             f    = 0.0_dp
        !         case (2)
        !             w(1) =   (E-e1)**2*1.0_dp/(e1-e2)**2*1.0_dp/(e1-e3)**2*1.0_dp/(e1-e4)**2*(E*(e4*(e2+e3)+e2*e3-e1*(e2+e3+e4)*2.0_dp+e1**2*3.0_dp)+e1*(e4*(e2+e3)+e2*e3)*2.0_dp-e1**2*(e2+e3+e4)-e2*e3*e4*3.0_dp)*(-4.0_dp)
        !             w(2) =  ((E-e1)**3*1.0_dp/(e1-e2)**2*4.0_dp)/((e1-e3)*(e1-e4))
        !             w(3) =  ((E-e1)**3*1.0_dp/(e1-e3)**2*4.0_dp)/((e1-e2)*(e1-e4))
        !             w(4) =  ((E-e1)**3*1.0_dp/(e1-e4)**2*4.0_dp)/((e1-e2)*(e1-e3))
        !             f    = -(E*2.4D1-e1*2.4D1)/((e1-e2)*(e1-e3)*(e1-e4))
        !         case (3)
        !             w(1) = (1.0_dp/(e1-e3)**2*1.0_dp/(e1-e4)**2*(E**2*(e1**2*e2*3.0_dp+e3*e4**2*3.0_dp+e3**2*e4*3.0_dp-e1*e3*e4*6.0_dp-e2*e3*e4*3.0_dp)-E**3*(e1*e2*2.0_dp-e1*e3*2.0_dp-e1*e4*2.0_dp-e2*e3-e2*e4+e3*e4+e1**2+e3**2+e4**2)-E*(e3**2*e4**2*3.0_dp+e1**2*e2*e3*3.0_dp+e1**2*e2*e4*3.0_dp-e1**2*e3*e4*3.0_dp-e1*e2*e3*e4*6.0_dp)+e1**2*e2*e3**2+e1**2*e2*e4**2+e1*e3**2*e4**2*2.0_dp-e1**2*e3*e4**2-e1**2*e3**2*e4+e2*e3**2*e4**2-e1*e2*e3*e4**2*2.0_dp-e1*e2*e3**2*e4*2.0_dp+e1**2*e2*e3*e4)*(-4.0_dp))/((e2-e3)*(e2-e4))
        !             w(2) = (1.0_dp/(e2-e3)**2*1.0_dp/(e2-e4)**2*(E**2*(e1*e2**2*3.0_dp+e3*e4**2*3.0_dp+e3**2*e4*3.0_dp-e1*e3*e4*3.0_dp-e2*e3*e4*6.0_dp)-E**3*(e1*e2*2.0_dp-e1*e3-e1*e4-e2*e3*2.0_dp-e2*e4*2.0_dp+e3*e4+e2**2+e3**2+e4**2)-E*(e3**2*e4**2*3.0_dp+e1*e2**2*e3*3.0_dp+e1*e2**2*e4*3.0_dp-e2**2*e3*e4*3.0_dp-e1*e2*e3*e4*6.0_dp)+e1*e2**2*e3**2+e1*e2**2*e4**2+e1*e3**2*e4**2+e2*e3**2*e4**2*2.0_dp-e2**2*e3*e4**2-e2**2*e3**2*e4-e1*e2*e3*e4**2*2.0_dp-e1*e2*e3**2*e4*2.0_dp+e1*e2**2*e3*e4)*(-4.0_dp))/((e1-e3)*(e1-e4))
        !             w(3) = (1.0_dp/(e1-e3)**2*1.0_dp/(e2-e3)**2*(E**2*(e1*e2**2*3.0_dp+e1**2*e2*3.0_dp+e3**2*e4*3.0_dp-e1*e2*e3*6.0_dp-e1*e2*e4*3.0_dp)-E**3*(e1*e2-e1*e3*2.0_dp-e1*e4-e2*e3*2.0_dp-e2*e4+e3*e4*2.0_dp+e1**2+e2**2+e3**2)-E*(e1**2*e2**2*3.0_dp-e1*e2*e3**2*3.0_dp+e1*e3**2*e4*3.0_dp+e2*e3**2*e4*3.0_dp-e1*e2*e3*e4*6.0_dp)-e1*e2**2*e3**2-e1**2*e2*e3**2+e1**2*e2**2*e3*2.0_dp+e1**2*e2**2*e4+e1**2*e3**2*e4+e2**2*e3**2*e4+e1*e2*e3**2*e4-e1*e2**2*e3*e4*2.0_dp-e1**2*e2*e3*e4*2.0_dp)*4.0_dp)/((e1-e4)*(e2-e4))
        !             w(4) = (1.0_dp/(e1-e4)**2*1.0_dp/(e2-e4)**2*(E**2*(e1*e2**2*3.0_dp+e1**2*e2*3.0_dp+e3*e4**2*3.0_dp-e1*e2*e3*3.0_dp-e1*e2*e4*6.0_dp)-E**3*(e1*e2-e1*e3-e1*e4*2.0_dp-e2*e3-e2*e4*2.0_dp+e3*e4*2.0_dp+e1**2+e2**2+e4**2)-E*(e1**2*e2**2*3.0_dp-e1*e2*e4**2*3.0_dp+e1*e3*e4**2*3.0_dp+e2*e3*e4**2*3.0_dp-e1*e2*e3*e4*6.0_dp)+e1**2*e2**2*e3-e1*e2**2*e4**2-e1**2*e2*e4**2+e1**2*e2**2*e4*2.0_dp+e1**2*e3*e4**2+e2**2*e3*e4**2+e1*e2*e3*e4**2-e1*e2**2*e3*e4*2.0_dp-e1**2*e2*e3*e4*2.0_dp)*4.0_dp)/((e1-e3)*(e2-e3))
        !             f    = (e1*e2*(-2.4D1)+e3*e4*2.4D1+E*(e1+e2-e3-e4)*2.4D1)/((e1-e3)*(e1-e4)*(e2-e3)*(e2-e4))
        !         case (4)
        !             w(1) =  ((E-e4)**3*1.0_dp/(e1-e4)**2*(-4.0_dp))/((e2-e4)*(e3-e4))
        !             w(2) =  ((E-e4)**3*1.0_dp/(e2-e4)**2*(-4.0_dp))/((e1-e4)*(e3-e4))
        !             w(3) =  ((E-e4)**3*1.0_dp/(e3-e4)**2*(-4.0_dp))/((e1-e4)*(e2-e4))
        !             w(4) =   (E-e4)**2*1.0_dp/(e1-e4)**2*1.0_dp/(e2-e4)**2*1.0_dp/(e3-e4)**2*(E*(e1*e2+e1*e3+e2*e3-e4*(e1*2.0_dp+e2*2.0_dp+e3*2.0_dp)+e4**2*3.0_dp)+e4*(e1*(e2+e3)*2.0_dp+e2*e3*2.0_dp)-e4**2*(e1+e2+e3)-e1*e2*e3*3.0_dp)*4.0_dp
        !             f    = -(E*2.4D1-e4*2.4D1)/((e1-e4)*(e2-e4)*(e3-e4))
        !         case (5)
        !             w(1) = 0.0_dp
        !             w(2) = 0.0_dp
        !             w(3) = 0.0_dp
        !             w(4) = 0.0_dp
        !             f    = 0.0_dp
        !         case default
        !             stop 'ERROR [tetrahedron_weights]: delta, invalid bracket'
        !         end select
        !     else
        !         stop 'ERROR [tetrahedron_weights]: integration type /= heavi or delta'
        !     endif
        !     ! check for NaNs
        !     if (any(isnan(w))) stop 'ERROR [tetrahedron_weights]: NaN obtained'
        !     ! apply Bloch correction
        !     if (index(flags,'blochl').ne.0) then
        !         f = f*0.25_dp
        !         do i = 1,4
        !             w(1)=w(1)+0.025_dp*f*(Ec(i)-e1)
        !             w(2)=w(2)+0.025_dp*f*(Ec(i)-e2)
        !             w(3)=w(3)+0.025_dp*f*(Ec(i)-e3)
        !             w(4)=w(4)+0.025_dp*f*(Ec(i)-e4)
        !         enddo
        !     endif
        !     ! divide by four because ... anyway seems right.
        !     w=w*0.25_dp
        !     !
        ! end function   get_tetrahedron_weight

    pure function  get_idos_ep(Ep,E,tet,tetw) result(idos)
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
        !
    end function   get_idos_ep

    pure function  get_dos_ep(Ep,E,tet,tetw) result(dosEp)
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
    end function   get_dos_ep

end module am_dos













