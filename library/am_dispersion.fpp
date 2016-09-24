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
        real(dp), allocatable :: E(:,:)   ! E(nbands,nKs) energies
        contains
        procedure :: save => save_dr
        procedure :: load => load_dr
        procedure, private :: write_dr
        procedure, private :: read_dr
        generic :: write(formatted) => write_dr
        generic :: read(formatted)  => read_dr
    end type am_class_dispersion

    type, public, extends(am_class_dispersion) :: am_class_tightbinding_dispersion
        complex(dp), allocatable :: H(:,:,:)    ! H(:,:,nKs) tight binding Hamiltonian
        complex(dp), allocatable :: C(:,:,:)    ! C(:,:,nKs) tight binding coefficients (eigenvectors)
        real(dp)   , allocatable :: weight(:,:) ! weights(nbands,nKs) whatever weights may be interesting to save
        complex(dp), allocatable :: Hk(:,:,:)   ! Hk(:,:,nKs) wave-vector derivitive of Hamiltonian, i.e. grad_k H
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
        ! clock in
        call start_clock('load')
        ! save space group
        open(unit=fid, file=trim(fname), status='old', action='read', iostat=iostat)
            if (iostat/=0) stop 'ERROR [dr:load]: opening file'
            read(fid,*) dr
        close(fid)
        ! clock out
        call stop_clock('load')
        !
    end subroutine  load_dr

    ! fourier interpolation

    pure function   fourier_interpolation(V,K,R,Kq,A) result(Vq)
        !
        implicit none
        !
        complex(dp), intent(in) :: V(:,:,:)       ! values on fbz
        real(dp)   , intent(in) :: K(:,:)         ! K(3,nKs) fbz kpoint coordinates [0,1)
        real(dp)   , intent(in) :: R(:,:)         ! R(3,nKs) fft mesh of fbz (primitive real-space lattice)
        real(dp)   , intent(in) :: Kq(:,:)        ! Kq(3,nKqs) points to interpolate on [0,1)
        complex(dp), intent(in), optional :: A(:) ! Augment transform with kernel A: for differentiation use A = sum(R(:,j)) * cmplx_i
        complex(dp),allocatable :: Vq(:,:,:)      ! interpolated values
        complex(dp),allocatable :: Vr(:,:,:)      ! fourier values
        integer :: nKs                            ! number of kpoints on fbz
        integer :: nRs                            ! number of kpoints on real mesh (primitive lattice points)
        integer :: nKqs                           ! number of kpoints to interpolate
        integer :: i, j, n, m
        complex(dp) :: ffac
        ! nbands
        n = size(V,1)
        m = size(V,2)
        ! get number of real points
        nRs = size(R,2)
        ! get kpoints
        nKs = size(K,2)
        ! get interpolation points
        nKqs = size(Kq,2)
        ! allocate real space
        allocate(Vr(n,m,nRs))
        ! allocate interpolated space
        allocate(Vq(n,m,nKqs))
        ! Fourier transform Hamiltonian to real space
        do j = 1, nRs
            ! initialize
            Vr(:,:,j) = 0.0_dp
            ! construct real-space hamiltonian (i.e. Fourier transform H)
            do i = 1, nKs
                ffac = exp(-itwopi*dot_product(K(:,i),R(:,j)))
                Vr(:,:,j) = Vr(:,:,j) + V(:,:,i) * ffac
            enddo
        enddo
        ! apply kernel
        if (present(A)) then
            do i = 1, nRs
                Vr(:,:,i) = A(i) * Vr(:,:,i)
            enddo
        endif
        ! Fourier transform back to reciprocal space
        do i = 1, nKqs
            ! initialize complex variable
            Vq(:,:,i) = 0.0_dp
            ! perform transform
            do j = 1, nRs
                ffac = exp(-itwopi*dot_product(Kq(:,i),R(:,j)))
                Vq(:,:,i) = Vq(:,:,i) + Vr(:,:,j) * ffac
            enddo
        enddo ! ki
        ! normalize
        Vq = Vq/nKs
    end function    fourier_interpolation

    pure function   cosine_interpolation(V,K,R,Kq) result(Vq)
        !
        ! Requires inversion symmetry at Gamma. Okay in reciprocal space if system has time-reversal symmetry -- even if crystal is noncentrosymmetric.
        !
        implicit none
        !
        real(dp), intent(in) :: V(:,:,:)  ! values on fbz
        real(dp), intent(in) :: K(:,:)    ! K(3,nKs) fbz kpoint coordinates [0,1)
        real(dp), intent(in) :: R(:,:)    ! R(3,nKs) fft mesh of fbz (primitive real-space lattice)
        real(dp), intent(in) :: Kq(:,:)   ! Kq(3,nKqs) points to interpolate on [0,1)
        real(dp),allocatable :: Vq(:,:,:) ! interpolated values
        real(dp),allocatable :: Vr(:,:,:) ! fourier values
        integer :: nKs                    ! number of kpoints on fbz
        integer :: nRs                    ! number of kpoints on real mesh (primitive lattice points)
        integer :: nKqs                   ! number of kpoints to interpolate
        integer :: i, j, n, m
        complex(dp) :: ffac
        ! nbands
        n = size(V,1)
        m = size(V,2)
        ! get number of real points
        nRs = size(R,2)
        ! get kpoints
        nKs = size(K,2)
        ! get interpolation points
        nKqs = size(Kq,2)
        ! allocate real space
        allocate(Vr(n,m,nRs))
        ! allocate interpolated space
        allocate(Vq(n,m,nKqs))
        ! Fourier transform Hamiltonian to real space
        do j = 1, nRs
            ! initialize
            Vr(:,:,j) = 0.0_dp
            ! construct real-space hamiltonian (i.e. Fourier transform H)
            do i = 1, nKs
                ffac = cos(twopi*dot_product(K(:,i),R(:,j)))
                Vr(:,:,j) = Vr(:,:,j) + V(:,:,i) * ffac
            enddo
        enddo
        ! Fourier transform back to reciprocal space
        do i = 1, nKqs
            ! initialize complex variable
            Vq(:,:,i) = 0.0_dp
            ! perform transform
            do j = 1, nRs
                ffac = cos(twopi*dot_product(Kq(:,i),R(:,j)))
                Vq(:,:,i) = Vq(:,:,i) + Vr(:,:,j) * ffac
            enddo
        enddo ! ki
        ! normalize
        Vq = Vq/nKs
    end function    cosine_interpolation

end module am_dispersion







