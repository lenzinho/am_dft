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

    #:setvar KINDS_interpolate_via_fft ['real(dp)', 'complex(dp)']
    #:setvar RANKS_interpolate_via_fft range(1,4)
    #:setvar PREFS_interpolate_via_fft ['{}{}_'.format(KIND[3], RANK) for KIND in KINDS_interpolate_via_fft for RANK in RANKS_interpolate_via_fft]

    interface interpolate_via_fft
        module procedure ${variants('interpolate_via_fft', prefixes=PREFS_interpolate_via_fft)}$
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

#:for KIND in KINDS_interpolate_via_fft
#:for RANK in RANKS_interpolate_via_fft
    pure function   ${KIND[3]}$${RANK}$_interpolate_via_fft(V,K,R,Kq,A) result(Vq)
        !
        implicit none
        !
        ${KIND}$, intent(in) :: V${ranksuffix(RANK)}$ ! values on fbz
        real(dp)   , intent(in) :: K(:,:)             ! K(3,nKs) fbz kpoint coordinates [0,1)
        real(dp)   , intent(in) :: R(:,:)             ! R(3,nKs) fft mesh of fbz (primitive real-space lattice)
        real(dp)   , intent(in) :: Kq(:,:)            ! Kq(3,nKqs) points to interpolate on [0,1)
        complex(dp), intent(in), optional :: A(:)     ! Augment transform with kernel A: for differentiation use A = sum(R(:,j)) * cmplx_i
        ${KIND}$   , allocatable :: Vq${ranksuffix(RANK)}$    ! interpolated values
        complex(dp), allocatable :: Vr${ranksuffix(RANK)}$    ! fourier values
        integer :: nKs  ! number of kpoints on fbz
        integer :: nRs  ! number of kpoints on real mesh (primitive lattice points)
        integer :: nKqs ! number of kpoints to interpolate
        integer :: i, j
        #:if RANK != 1
            integer :: n
            integer, allocatable :: vshape(:)
            complex(dp), allocatable :: wrk${ranksuffix(RANK-1)}$ ! workspace
        #:else
            complex(dp) :: wrk ! workspace
        #:endif
        ! get number of real points
        nRs = size(R,2)
        ! get kpoints
        nKs = size(K,2)
        ! get interpolation points
        nKqs = size(Kq,2)
        ! allocate space
        #:if RANK == 1
            ! for Vr
            allocate(Vr(nRs))
            ! for Vq
            allocate(Vq(nKqs))
        #:else
            ! get shape
            vshape = shape(V)
            ! rank
            n = size(vshape)
            ! set shape for Vr
            vshape(n) = nRs
            ! allocate real space
            call vector_allocate(A=Vr,vec=vshape)
            ! set shape for Vq
            vshape(n) = nKqs
            ! allocate interpolated space
            call vector_allocate(A=Vq,vec=vshape)
            ! wrkspace
            call vector_allocate(wrk,vshape(1:(n-1)))
        #:endif
        ! initialize
        Vr = 0.0_dp
        ! Fourier transform Hamiltonian to real space
        !$OMP PARALLEL PRIVATE(i,j) SHARED(nRs,nKs,Vr,K,R)
        !$OMP DO
        do j = 1, nRs
            ! construct real-space hamiltonian (i.e. Fourier transform H)
            do i = 1, nKs
                Vr${ranksuffix_last(RANK,last='j')}$ = Vr${ranksuffix_last(RANK,last='j')}$ &
                    & + V${ranksuffix_last(RANK,last='i')}$ * exp(-itwopi*dot_product(K(:,i),R(:,j)))
            enddo
        enddo
        !$OMP END DO
        !$OMP END PARALLEL
        ! apply kernel
        if (present(A)) then
            do j = 1, nRs
                Vr${ranksuffix_last(RANK,last='j')}$ = A(j) * Vr${ranksuffix_last(RANK,last='j')}$
            enddo
        endif
        ! Fourier transform back to reciprocal space
        !$OMP PARALLEL PRIVATE(i,j,wrk) SHARED(nKqs,nRs,Vr,Kq,R,Vq)
        !$OMP DO
        do i = 1, nKqs
            ! initialize complex variable
            wrk = 0.0_dp
            ! perform transform
            do j = 1, nRs
                wrk = wrk + Vr${ranksuffix_last(RANK,last='j')}$ * exp(-itwopi*dot_product(Kq(:,i),R(:,j)))
            enddo
            #:if KIND == 'real(dp)'
                ! get real valued function
                Vq${ranksuffix_last(RANK,last='i')}$ = real(wrk,dp)
            #:else
                Vq${ranksuffix_last(RANK,last='i')}$ = wrk
            #:endif
        enddo ! ki
        !$OMP END DO
        !$OMP END PARALLEL
        ! normalize
        Vq = Vq/nKs
    end function    ${KIND[3]}$${RANK}$_interpolate_via_fft
#:endfor
#:endfor

end module am_dispersion







