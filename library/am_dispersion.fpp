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

    interface interpolate_via_fft
        module procedure :: c_interpolate_via_fft, r_interpolate_via_fft
    end interface ! interpolate_via_fft

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

    pure function   c_interpolate_via_fft(V,K,R,Kq,A) result(Vq)
        !
        implicit none
        !
        complex(dp), intent(in) :: V(:)    ! values on fbz
        real(dp)   , intent(in) :: K(:,:)  ! K(3,nKs) fbz kpoint coordinates [0,1)
        real(dp)   , intent(in) :: R(:,:)  ! R(3,nKs) fft mesh of fbz (primitive real-space lattice)
        real(dp)   , intent(in) :: Kq(:,:) ! Kq(3,nKqs) points to interpolate on [0,1)
        complex(dp), intent(in), optional :: A(:)    ! Augment transform with kernel A: for differentiation use A = sum(R(:,j)) * cmplx_i
        complex(dp),allocatable :: Vq(:)   ! interpolated values
        complex(dp),allocatable :: Vr(:)   ! fourier values
        integer :: nKs  ! number of kpoints on fbz
        integer :: nRs  ! number of kpoints on real mesh (primitive lattice points)
        integer :: nKqs ! number of kpoints to interpolate
        integer :: i, j
        complex(dp) :: wrk
        ! R (FFT MESH) is passed as input and can be generated like so:
        ! R = meshgrid( v1=real([1:fbz%n(1)]-1-floor(fbz%n(1)/2.0_dp),dp), &
        !             & v2=real([1:fbz%n(2)]-1-floor(fbz%n(2)/2.0_dp),dp), &
        !             & v3=real([1:fbz%n(3)]-1-floor(fbz%n(3)/2.0_dp),dp))
        ! get number of real points
        nRs = size(R,2)
        ! get kpoints
        nKs = size(K,2)
        ! get interpolation points
        nKqs = size(Kq,2)
        ! allocate real space
        allocate(Vr(nRs))
        ! allocate interpolated space
        allocate(Vq(nKqs))
        ! Fourier transform Hamiltonian to real space
        !$OMP PARALLEL PRIVATE(i,j) SHARED(nRs,nKs,Vr,K,R)
        !$OMP DO
        do j = 1, nRs
            ! initialize
            Vr(j) = 0.0_dp
            ! construct real-space hamiltonian (i.e. Fourier transform H)
            do i = 1, nKs
                Vr(j) = Vr(j) + V(i) * exp(-itwopi*dot_product(K(:,i),R(:,j)))
            enddo
        enddo
        !$OMP END DO
        !$OMP END PARALLEL
        ! apply kernel
        if (present(A)) then
            Vr = A * Vr
        endif
        ! Fourier transform back to reciprocal space
        !$OMP PARALLEL PRIVATE(i,j,wrk) SHARED(nKqs,nRs,Vr,Kq,R,Vq)
        !$OMP DO
        do i = 1, nKqs
            ! initialize complex variable
            wrk = 0.0_dp
            ! perform transform
            do j = 1, nRs
                wrk = wrk + Vr(j) * exp(-itwopi*dot_product(Kq(:,i),R(:,j)))
            enddo
            ! get real valued function
            Vq(i) = wrk
        enddo ! ki
        !$OMP END DO
        !$OMP END PARALLEL
        ! normalize
        Vq = Vq/nKs
    end function    c_interpolate_via_fft

    pure function   r_interpolate_via_fft(V,K,R,Kq,A) result(Vq)
        !
        implicit none
        !
        real(dp)   , intent(in) :: V(:)    ! values on fbz
        real(dp)   , intent(in) :: K(:,:)  ! K(3,nKs) fbz kpoint coordinates [0,1)
        real(dp)   , intent(in) :: R(:,:)  ! R(3,nKs) fft mesh of fbz (primitive real-space lattice)
        real(dp)   , intent(in) :: Kq(:,:) ! Kq(3,nKqs) points to interpolate on [0,1)
        complex(dp), intent(in), optional :: A(:)    ! Augment transform with kernel A: for differentiation use A = sum(R(:,j)) * cmplx_i
        real(dp)   ,allocatable :: Vq(:)   ! interpolated values
        complex(dp),allocatable :: Vr(:)   ! fourier values
        integer :: nKs  ! number of kpoints on fbz
        integer :: nRs  ! number of kpoints on real mesh (primitive lattice points)
        integer :: nKqs ! number of kpoints to interpolate
        integer :: i, j
        complex(dp) :: wrk
        ! R (FFT MESH) is passed as input and can be generated like so:
        ! R = meshgrid( v1=real([1:fbz%n(1)]-1-floor(fbz%n(1)/2.0_dp),dp), &
        !             & v2=real([1:fbz%n(2)]-1-floor(fbz%n(2)/2.0_dp),dp), &
        !             & v3=real([1:fbz%n(3)]-1-floor(fbz%n(3)/2.0_dp),dp))
        ! get number of real points
        nRs = size(R,2)
        ! get kpoints
        nKs = size(K,2)
        ! get interpolation points
        nKqs = size(Kq,2)
        ! allocate real space
        allocate(Vr(nRs))
        ! allocate interpolated space
        allocate(Vq(nKqs))
        ! Fourier transform Hamiltonian to real space
        !$OMP PARALLEL PRIVATE(i,j) SHARED(nRs,nKs,Vr,K,R)
        !$OMP DO
        do j = 1, nRs
            ! initialize
            Vr(j) = 0.0_dp
            ! construct real-space hamiltonian (i.e. Fourier transform H)
            do i = 1, nKs
                Vr(j) = Vr(j) + V(i) * exp(-itwopi*dot_product(K(:,i),R(:,j)))
            enddo
        enddo
        !$OMP END DO
        !$OMP END PARALLEL
        ! apply kernel
        if (present(A)) then
            Vr = A * Vr
        endif
        ! Fourier transform back to reciprocal space
        !$OMP PARALLEL PRIVATE(i,j,wrk) SHARED(nKqs,nRs,Vr,Kq,R,Vq)
        !$OMP DO
        do i = 1, nKqs
            ! initialize complex variable
            wrk = 0.0_dp
            ! perform transform
            do j = 1, nRs
                wrk = wrk + Vr(j) * exp(-itwopi*dot_product(Kq(:,i),R(:,j)))
            enddo
            ! get real valued function
            Vq(i) = real(wrk,dp)
        enddo ! ki
        !$OMP END DO
        !$OMP END PARALLEL
        ! normalize
        Vq = Vq/nKs
    end function    r_interpolate_via_fft

end module am_dispersion







