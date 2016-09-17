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
        complex(dp), allocatable :: Hk(:,:,:)   ! Hk(:,:,nkpts) momentum-derivitive of Hamiltonian, i.e. grad_k H
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

end module am_dispersion







