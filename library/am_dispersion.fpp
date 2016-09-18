#:include "fypp_macros.fpp"
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
    use am_shells

    implicit none

    type, public :: am_class_dispersion
        integer :: nbands
        real(dp), allocatable :: E(:,:)   ! E(nbands,nkpts) energies
        contains
        procedure :: save => save_dr
        procedure :: load => load_dr
        procedure, private :: write_dr
        procedure, private :: read_dr
        generic :: write(formatted) => write_dr
        generic :: read(formatted)  => read_dr
    end type am_class_dispersion

    type, public, extends(am_class_dispersion) :: am_class_tightbinding_dispersion
        complex(dp), allocatable :: H(:,:,:)    ! H(:,:,nkpts) tight binding Hamiltonian
        complex(dp), allocatable :: C(:,:,:)    ! C(:,:,nkpts) tight binding coefficients (eigenvectors)
        real(dp)   , allocatable :: weight(:,:) ! weights(nbands,nkpts) whatever weights may be interesting to save
        complex(dp), allocatable :: Hk(:,:,:)   ! Hk(:,:,nkpts) wave-vector derivitive of Hamiltonian, i.e. grad_k H
    end type am_class_tightbinding_dispersion

contains

    ! i/o

    subroutine      write_dr(dtv, unit, iotype, v_list, iostat, iomsg)
        !
        implicit none
        !
        class(am_class_dispersion), intent(in) :: dtv
        integer     , intent(in)    :: unit
        character(*), intent(in)    :: iotype
        integer     , intent(in)    :: v_list(:)
        integer     , intent(out)   :: iostat
        character(*), intent(inout) :: iomsg
        !
        iostat = 0
        !
        if (iotype.eq.'LISTDIRECTED') then
            write(unit,'(a/)') '<dr>'
                ! non-allocatable
                #:for ATTRIBUTE in ['nbands']
                    $:write_xml_attribute_nonallocatable(ATTRIBUTE)
                #:endfor
                ! allocatable        
                #:for ATTRIBUTE in ['E']
                    $:write_xml_attribute_allocatable(ATTRIBUTE)
                #:endfor
                ! type-specific stuff
                select type (dtv)
                class is (am_class_tightbinding_dispersion)
                    ! allocatable
                    #:for ATTRIBUTE in ['C','weight']
                        $:write_xml_attribute_allocatable(ATTRIBUTE)
                    #:endfor
                class default
                    ! do nothing
                end select
            write(unit,'(a/)') '</dr>'
        else
            stop 'ERROR [write_dr]: iotype /= LISTDIRECTED'
        endif
    end subroutine  write_dr

    subroutine      read_dr(dtv, unit, iotype, v_list, iostat, iomsg)
        !
        implicit none
        !
        class(am_class_dispersion), intent(inout) :: dtv
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
                #:for ATTRIBUTE in ['nbands']
                    $:read_xml_attribute_nonallocatable(ATTRIBUTE)
                #:endfor
                ! allocatable        
                #:for ATTRIBUTE in ['E']
                    $:read_xml_attribute_allocatable(ATTRIBUTE)
                #:endfor
                ! type-specific stuff
                select type (dtv)
                class is (am_class_tightbinding_dispersion)
                    ! allocatable
                    #:for ATTRIBUTE in ['C','weight']
                        $:read_xml_attribute_allocatable(ATTRIBUTE)
                    #:endfor
                class default
                    ! do nothing
                end select
            read(unit=unit,fmt='(/)',iostat=iostat,iomsg=iomsg)
            ! set iostat
            iostat=-1
        else
            stop 'ERROR [read_dr]: iotype /= LISTDIRECTED'
        endif
    end subroutine  read_dr

    subroutine      save_dr(dr,fname)
        !
        implicit none
        !
        class(am_class_dispersion), intent(in) :: dr
        character(*), intent(in) :: fname
        integer :: fid
        integer :: iostat
        ! fid
        fid = 1
        ! save space group
        open(unit=fid, file=trim(fname), status='replace', action='write', iostat=iostat)
            if (iostat/=0) stop 'ERROR [dr:load]: opening file'
            write(fid,*) dr
        close(fid)
        !
    end subroutine  save_dr

    subroutine      load_dr(dr,fname)
        !
        implicit none
        !
        class(am_class_dispersion), intent(inout) :: dr
        character(*), intent(in) :: fname
        integer :: fid
        integer :: iostat
        ! fid
        fid = 1
        ! save space group
        open(unit=fid, file=trim(fname), status='old', action='read', iostat=iostat)
            if (iostat/=0) stop 'ERROR [dr:load]: opening file'
            read(fid,*) dr
        close(fid)
        !
    end subroutine  load_dr

    ! fourier interpolation

    pure function   interpolate_via_fft(V,kpt_recp,rpt_frac,kpt_recp_i) result(Vi)
        ! V(nkpts) - value to be interpolated, defined on full brillouin zone mesh 
        ! kpt_recp(3,nkpts) - full bz [0,1) reciprocal fractional coordinates 
        ! rpt_frac(r,nkpts) - real space integer fft mesh
        ! kpt_recp_i - points to interpolate V on
        implicit none
        !
        real(dp), intent(in) :: V(:) ! values on fbz
        real(dp), intent(in) :: kpt_recp(:,:)   ! kpt_recp(3,nkpts) fbz kpoint coordinates [0,1)
        real(dp), intent(in) :: rpt_frac(:,:)   ! rpt_frac(3,nkpts) fft mesh of fbz (primitive real-space lattice)
        real(dp), intent(in) :: kpt_recp_i(:,:) ! kpt_recp_i(3,nkpts_i) points to interpolate on [0,1)
        real(dp),allocatable :: Vi(:) ! interpolated values
        real(dp),allocatable :: Vr(:) ! fourier values
        complex(dp) :: wrk
        integer :: nkpts   ! number of kpoints on fbz
        integer :: nrpts   ! number of kpoints on real mesh (primitive lattice points)
        integer :: nkpts_i ! number of kpoints to interpolate
        integer :: i, j
        ! get number of real points
        nrpts = size(rpt_frac,2)
        ! get kpoints
        nkpts = size(kpt_recp)
        ! get interpolation points
        nkpts_i = size(kpt_recp_i)
        ! allocate real space
        allocate(Vr(nrpts))
        ! allocate interpolated space
        allocate(Vi(nkpts_i))
        ! Fourier transform Hamiltonian to real space
        !$OMP PARALLEL PRIVATE(i,j) SHARED(nrpts,nkpts,Vr,kpt_recp,rpt_frac)
        !$OMP DO
        do j = 1, nrpts
            ! initialize
            Vr(j) = 0.0_dp
            ! construct real-space hamiltonian (i.e. Fourier transform H)
            do i = 1, nkpts
                Vr(j) = Vr(j) + V(i) * exp(-itwopi*dot_product(kpt_recp(:,i),rpt_frac(:,j)))
            enddo
        enddo
        !$OMP END DO
        !$OMP END PARALLEL
        ! Fourier transform back to reciprocal space (only do point of interest this time)
        !$OMP PARALLEL PRIVATE(i,j,wrk) SHARED(nkpts_i,nrpts,Vr,kpt_recp_i,rpt_frac)
        !$OMP DO
        do i = 1, nkpts_i
            ! initialize complex variable
            wrk = 0.0_dp
            ! perform transform
            do j = 1, nrpts
                wrk = wrk + Vr(j) * exp(-itwopi*dot_product(kpt_recp_i(:,i),rpt_frac(:,j)))
            enddo
            ! get real valued function
            Vi(i) = real(wrk,dp)
        enddo ! ki
        !$OMP END DO
        !$OMP END PARALLEL
        ! normalize
        Vi = Vi/nkpts
    end function    interpolate_via_fft

end module am_dispersion







