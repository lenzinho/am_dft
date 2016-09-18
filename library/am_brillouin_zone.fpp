#:include "fypp_macros.fpp"
module am_brillouin_zone

    use am_unit_cell
    use am_symmetry
    use am_matlab
    use am_constants
    use am_stdout
    use am_vasp_io
    use am_mkl
    use am_options
    use dispmodule
    use am_tet_mesh
    use am_rank_and_sort

    implicit none

    private

    ! BZ (primitive reciprocal-lattice cell defined with x,y,z (fractional) between [0,1).)

    type, public :: am_class_bz
        real(dp) :: bas(3,3)
        real(dp) :: recbas(3,3)
        integer  :: nkpts                       ! nkpts number of kpoints
        real(dp), allocatable :: kpt_cart(:,:)  ! kpt(3,nkpts) kpoint
        real(dp), allocatable :: kpt_recp(:,:)  ! kpt(3,nkpts) kpoint - notice that vasp IBZKPT is in reciprocal (fractional) units
        contains
        procedure :: create_bz
        procedure :: save => save_bz
        procedure :: load => load_bz
        procedure, private :: write_bz
        procedure, private :: read_bz
        generic :: write(formatted) => write_bz
        generic :: read(formatted) => read_bz
        ! add sort procedure
    end type am_class_bz

    ! FBZ (Monkhorst-pack grid; cart points defined in wigner-seitz cell; frac points defined between [0,1) )

    type, public, extends(am_class_bz) :: am_class_fbz
        integer  :: n(3) ! mesh dimensions
        real(dp) :: s(3) ! shift
        integer  :: ntets                     ! number of tetrahedra
        real(dp), allocatable :: w(:)         ! w(nkpts) normalized weights
        integer , allocatable :: kpt_int(:,:) ! kpt(3,nkpts) kpoint integer indices in regular grid
        integer , allocatable :: tet(:,:)     ! tetrahedra connectivity list (indices of corner kpoints)
        real(dp), allocatable :: tetv(:)      ! tetrahedra volume
        real(dp), allocatable :: tetw(:)      ! tetrahedron weight
        integer , allocatable :: fbz_id(:)
        integer , allocatable :: ibz_id(:)
        contains
        procedure :: get_fbz
    end type am_class_fbz

    ! IBZ (irreducible wedge)

    type, public, extends(am_class_fbz) :: am_class_ibz
        contains
        procedure :: get_irreducible ! requires fbz monkhorst-pack mesh
    end type am_class_ibz

    ! PATH 

    type, public, extends(am_class_bz) :: am_class_path
        integer :: ndivs ! number of divisions
        integer :: nticks ! umber of ticks
        character(:), allocatable :: kpt_symb(:)
        real(dp)    , allocatable :: tick(:)
        real(dp)    , allocatable :: x(:)
        contains
        procedure :: get_path
    end type am_class_path

    ! other stuff

    public :: get_pathsymb

contains

    ! brillouin-zone i/o

    subroutine     write_bz(dtv, unit, iotype, v_list, iostat, iomsg)
        !
        implicit none
        !
        class(am_class_bz), intent(in) :: dtv
        integer     , intent(in)    :: unit
        character(*), intent(in)    :: iotype
        integer     , intent(in)    :: v_list(:)
        integer     , intent(out)   :: iostat
        character(*), intent(inout) :: iomsg
        !
        iostat = 0
        !
        if (iotype.eq.'LISTDIRECTED') then
            write(unit,'(a/)') '<brillouin_zone>'
                ! non-allocatable
                #:for ATTRIBUTE in ['bas','recbas','nkpts']
                    $:write_xml_attribute_nonallocatable(ATTRIBUTE)
                #:endfor
                ! allocatable
                #:for ATTRIBUTE in ['kpt_cart','kpt_recp']
                    $:write_xml_attribute_allocatable(ATTRIBUTE)
                #:endfor
                ! type-specific stuff
                select type (dtv)
                class is (am_class_fbz)
                    ! non-allocatable
                    #:for ATTRIBUTE in ['n','s','ntets']
                        $:write_xml_attribute_nonallocatable(ATTRIBUTE)
                    #:endfor
                    ! allocatable
                    #:for ATTRIBUTE in ['kpt_int','tet','tetv','tetw','fbz_id','ibz_id','w']
                        $:write_xml_attribute_allocatable(ATTRIBUTE)
                    #:endfor
                class is (am_class_path)
                    ! non-allocatable
                    #:for ATTRIBUTE in ['ndivs','nticks']
                        $:write_xml_attribute_nonallocatable(ATTRIBUTE)
                    #:endfor
                    ! allocatable
                    #:for ATTRIBUTE in ['tick','x']
                        $:write_xml_attribute_allocatable(ATTRIBUTE)
                    #:endfor
                    ! alloctable string
                    #:for ATTRIBUTE in ['kpt_symb']
                        $:write_xml_attribute_allocatable_string(ATTRIBUTE)
                    #:endfor
                class default
                    ! do nothing
                end select
            write(unit,'(a/)') '</brillouin_zone>'
        else
            stop 'ERROR [write_bz]: iotype /= LISTDIRECTED'
        endif
    end subroutine write_bz

    subroutine     read_bz(dtv, unit, iotype, v_list, iostat, iomsg)
        !
        implicit none
        !
        class(am_class_bz), intent(inout) :: dtv
        integer     , intent(in)    :: unit
        character(*), intent(in)    :: iotype
        integer     , intent(in)    :: v_list(:)
        integer     , intent(out)   :: iostat
        character(*), intent(inout) :: iomsg
        logical :: isallocated
        integer :: dims_rank
        integer :: dims(10)
        integer :: str_length
        !
        if (iotype.eq.'LISTDIRECTED') then
            read(unit=unit,fmt='(/)',iostat=iostat,iomsg=iomsg)
                ! non-allocatable
                #:for ATTRIBUTE in ['bas','recbas','nkpts']
                    $:read_xml_attribute_nonallocatable(ATTRIBUTE)
                #:endfor
                ! allocatable
                #:for ATTRIBUTE in ['kpt_cart','kpt_recp']
                    $:read_xml_attribute_allocatable(ATTRIBUTE)
                #:endfor
                ! type-specific stuff
                select type (dtv)
                class is (am_class_fbz)
                    ! non-allocatable
                    #:for ATTRIBUTE in ['n','s','ntets']
                        $:read_xml_attribute_nonallocatable(ATTRIBUTE)
                    #:endfor
                    ! allocatable
                    #:for ATTRIBUTE in ['kpt_int','tet','tetv','tetw','fbz_id','ibz_id','w']
                        $:read_xml_attribute_allocatable(ATTRIBUTE)
                    #:endfor
                class is (am_class_path)
                    ! non-allocatable
                    #:for ATTRIBUTE in ['ndivs','nticks']
                        $:read_xml_attribute_nonallocatable(ATTRIBUTE)
                    #:endfor
                    ! allocatable
                    #:for ATTRIBUTE in ['tick','x']
                        $:read_xml_attribute_allocatable(ATTRIBUTE)
                    #:endfor
                    ! alloctable string
                    #:for ATTRIBUTE in ['kpt_symb']
                        $:read_xml_attribute_allocatable_string(ATTRIBUTE)
                    #:endfor
                class default
                    ! do nothing
                end select
            read(unit=unit,fmt='(/)',iostat=iostat,iomsg=iomsg)
            ! return iostat
            iostat=-1
        else
            stop 'ERROR [read_bz]: iotype /= LISTDIRECTED'
        endif
    end subroutine read_bz

    subroutine     save_bz(bz,fname)
        !
        implicit none
        !
        class(am_class_bz), intent(in) :: bz
        character(*), intent(in) :: fname
        integer :: fid
        integer :: iostat
        ! fid
        fid = 1
        ! save space group
        open(unit=fid, file=trim(fname), status='replace', action='write', iostat=iostat)
            if (iostat/=0) stop 'ERROR [bz:load]: opening file'
            write(fid,*) bz
        close(fid)
        !
    end subroutine save_bz

    subroutine     load_bz(bz,fname)
        !
        implicit none
        !
        class(am_class_bz), intent(inout) :: bz
        character(*), intent(in) :: fname
        integer :: fid
        integer :: iostat
        ! fid
        fid = 1
        ! save space group
        open(unit=fid, file=trim(fname), status='old', action='read', iostat=iostat)
            if (iostat/=0) stop 'ERROR [bz:load]: opening file'
            read(fid,*) bz
        close(fid)
        !
    end subroutine load_bz

    ! create

    subroutine     create_bz(bz,kpt_recp,kpt_cart,bas,prec)
        !
        ! kpt_recp points are shifted to be between [0,1)
        ! kpt_cart        are shifted to be in the wigner-seitz voronoi cell
        !
        implicit none
        !
        class(am_class_bz), intent(out) :: bz
        real(dp), intent(in), optional :: kpt_recp(:,:)
        real(dp), intent(in), optional :: kpt_cart(:,:)
        real(dp), intent(in) :: bas(3,3)
        real(dp), intent(in) :: prec
        real(dp) :: grid_points(3,27) ! voronoi points
        integer  :: i
        !
        bz%bas = bas
        bz%recbas= inv(bas)
        ! copy kpoint coordiantes [cart.] and [rec. frac.]
        if (present(kpt_recp)) then
            !
            bz%nkpts = size(kpt_recp,2)
            allocate(bz%kpt_cart(3,bz%nkpts))
            allocate(bz%kpt_recp(3,bz%nkpts))
            !
            bz%kpt_recp = kpt_recp
            bz%kpt_cart = matmul(bz%recbas,kpt_recp)
            !
        elseif (present(kpt_cart)) then
            !
            bz%nkpts = size(kpt_cart,2)
            allocate(bz%kpt_cart(3,bz%nkpts))
            allocate(bz%kpt_recp(3,bz%nkpts))
            !
            bz%kpt_recp = matmul(bz%bas,kpt_cart)
            bz%kpt_cart = kpt_cart
        else
            stop 'ERROR [create_bz]: either kpt_cart or kpt_recp must be present'
        endif
        ! generate voronoi points [cart]
        grid_points = meshgrid([-1:1],[-1:1],[-1:1])
        grid_points = matmul(bz%recbas,grid_points)
        ! make sure kpt_cart is in wigner-seitz cell 
        do i = 1, bz%nkpts
            bz%kpt_cart(:,i) = reduce_kpoint_to_fbz(kpoint_cart=bz%kpt_cart(:,i), grid_points=grid_points)
        enddo
        ! make sure kpt_recp is betwen [0,1)
        bz%kpt_recp = modulo(bz%kpt_recp+prec,1.0_dp)-prec
        !
        contains
        function       reduce_kpoint_to_fbz(kpoint_cart,grid_points) result(kpoint_fbz_cart)
            !> reduces kpoint (in fractional) to the first Brillouin zone (Wigner-Seitz cell, defined in cartesian coordinates)
            !> cartesian kpoint is returned! 
            implicit none
            !
            real(dp), intent(in) :: kpoint_cart(3) !> fractional
            real(dp), intent(in) :: grid_points(3,27) !> voronoi points (cartesian)
            real(dp) :: kpoint_fbz_cart(3) !> kpoint cartesian
            real(dp) :: G(3) !> reciprocal lattice vector
            real(dp) :: P !> bragg plane condition
            integer  :: i ! loop variable
            logical  :: is_not_done
            !
            ! copy kpoint
            kpoint_fbz_cart = kpoint_cart
            ! translating the k-point until the closest reciprocal lattice point is [0 0 0]
            is_not_done = .true.
            do while ( is_not_done )
                is_not_done = .false.
                do i = 1, 27
                    G = grid_points(:,i)
                    P = 2*dot_product(kpoint_fbz_cart,G) - dot_product(G,G)
                    if ( P.gt. tiny ) then
                        kpoint_fbz_cart = kpoint_fbz_cart - G
                        is_not_done = .true.
                    endif
                enddo
            end do
            !
        end function   reduce_kpoint_to_fbz
    end subroutine create_bz

    subroutine     get_fbz(fbz,pc,n,s,opts)
        !
        implicit none
        !
        class(am_class_fbz)    , intent(out):: fbz
        type(am_class_prim_cell), intent(in) :: pc
        integer                 , intent(in) :: n(3)
        real(dp)                , intent(in) :: s(3)
        type(am_class_options)  , intent(in) :: opts
        integer :: nkpts
        integer, allocatable :: kpt_int(:,:)
        real(dp),allocatable :: kpt_recp(:,:)
        real(dp) :: recvol
        integer :: i
        !
        if (opts%verbosity.ge.1) call print_title('Full Monkhorst-Pack mesh')
        ! save mesh dimensions and shift
        fbz%n = n
        fbz%s = s
        ! get regular integer monkhorst pack mesh
        kpt_int = meshgrid(v1=[1:fbz%n(1)]-1,v2=[1:fbz%n(2)]-1,v3=[1:fbz%n(3)]-1)
        ! nkpts
        nkpts = product(fbz%n)
        ! allocate space for frac
        allocate(kpt_recp(3,nkpts))
        ! convert integer mesh to fractional coordinates
        do i = 1, nkpts
            kpt_recp(:,i) = ( kpt_int(:,i) + s ) / real(n,dp)
        enddo 
        ! create instance
        call fbz%create_bz(kpt_recp=kpt_recp, bas=pc%bas, prec=opts%prec)
        ! save integer mesh
        allocate(fbz%kpt_int,source=kpt_int)
        ! allocate weights
        allocate(fbz%w(fbz%nkpts))
        ! get kpoint weights
        fbz%w = 1.0_dp/real(fbz%nkpts)
        ! stdout
        if (opts%verbosity.ge.1) then
            ! kpoint properties
            write(*,*) flare, 'kpoints = '//tostring(fbz%nkpts)
            call disp_indent()
            call disp(title='#'      ,style='underline',X=[1:fbz%nkpts]          ,fmt='i5'   ,advance='no')
            call disp(title='[cart.]',style='underline',X=transpose(fbz%kpt_cart),fmt='f10.5',advance='no')
            call disp(title='[frac.]',style='underline',X=transpose(fbz%kpt_recp),fmt='f10.5',advance='no')
            call disp(title='w'      ,style='underline',X=          fbz%w        ,fmt='f10.5',advance='yes')
            !
            write(*,'(a,a)') flare, 'Definitions:'
            write(*,'(5x,a)') 'cart points are reduced to voronoi wigner-seitz cell'
            write(*,'(5x,a)') 'frac points are reduced to between [0,1)'
        endif
        ! now get tetrahedra
        if (opts%verbosity.ge.1) call print_title('FBZ Tetrahedra')
        ! get tetrahedra
        fbz%tet = get_tetrahedra(recbas=pc%recbas,kpt_int=fbz%kpt_int,n=fbz%n)
        ! get number of tetrahedra
        fbz%ntets = size(fbz%tet,2)
        ! allocate tetra properties 
        allocate(fbz%tetv(fbz%ntets))
        allocate(fbz%tetw(fbz%ntets))
        ! primitive reciprocal cell volume
        recvol = det(fbz%recbas)
        ! get tetra volume = (volume/cell) / (nboxes/cell * ntets/boxes) = (volume/cell) / (tetra/cell) = volume/tetra
        fbz%tetv = recvol/(fbz%nkpts*6)
        ! get tetra weight
        fbz%tetw = fbz%tetv/recvol
        ! print stdout
        if (opts%verbosity.ge.1) then
            ! tetrahedra properties
            write(*,*) flare, 'tetrahedra = '//tostring(fbz%ntets)
            write(*,*) flare, 'tetrahedra volume = '//tostring(fbz%tetv(1))
            write(*,*) flare, 'tetrahedra weight = '//tostring(fbz%tetw(1))
            ! call disp_indent()
            ! call disp(title='#'                ,style='underline',X=[1:fbz%ntets]     ,fmt='i5' ,advance='no')
            ! call disp(title='connectivity list',style='underline',X=transpose(fbz%tet),fmt='i10',advance='yes')
        endif
    end subroutine get_fbz

    subroutine     get_irreducible(ibz,fbz,pg,opts)
        !
        implicit none
        !
        class(am_class_ibz)       , intent(out)   :: ibz
        type(am_class_fbz)        , intent(inout) :: fbz
        type(am_class_point_group), intent(in)    :: pg
        type(am_class_options)    , intent(in)    :: opts
        integer, allocatable :: PM(:,:)
        logical, allocatable :: mask(:)
        integer, allocatable :: ind(:)
        integer, allocatable :: ibz_tet(:,:)
        integer, allocatable :: ibz_tet_inds(:)
        integer, allocatable :: w(:)
        integer :: i,j,k
        !
        if (opts%verbosity.ge.1) call print_title('Irreducible Wedge (IBZ)')
        ! get permutation map which shows how k-points are permuted by each space symmetry operation, PM(fbz%nkpts,pg%nsyms) 
        PM = permutation_map(seitz=pg%seitz_recp, tau=fbz%kpt_recp, flags='', prec=opts%prec)
        ! determine irreducible kpoints
        allocate(mask(fbz%nkpts))
        mask = .true.
        allocate(ind(fbz%nkpts))
        ind = 0
        k=0
        do i = 1, fbz%nkpts
            if (mask(i)) then
                k=k+1
                ind(k)=i
                do j = 1, pg%nsyms
                    mask(PM(i,j))=.false.
                enddo
            endif
        enddo
        ! get number of k-points
        ibz%nkpts = k
        ! transfer irreducible k-points (frac)
        allocate(ibz%kpt_recp, source=fbz%kpt_recp(:,ind(1:k)))
        ! transfer irreducible k-points (cart)
        allocate(ibz%kpt_cart, source=fbz%kpt_cart(:,ind(1:k)))
        ! map irreducible k-point -> irreducible k-point
        allocate(ibz%ibz_id, source=[1:ibz%nkpts])
        ! map irreducible k-point onto -> primitive k-point
        allocate(ibz%fbz_id, source=ind(1:k))
        ! map primitive k-point onto -> irreducible k-point
        allocate(fbz%ibz_id(fbz%nkpts))
        fbz%ibz_id = 0
        do i = 1, ibz%nkpts
            ! PM(1,:) shows all k-points onto which k-point 1 is mapped by all point symmetry operations
            do j = 1, fbz%nkpts
                search : do k = 1, pg%nsyms
                    if (PM(j,k).eq.ibz%fbz_id(i)) then
                        fbz%ibz_id(j) = i
                        exit search
                    endif
                enddo search
            enddo
        enddo
        ! check for error
        if (any(fbz%ibz_id.eq.0)) stop 'ERROR [get_ibz]: not all fbz points were mapped onto an ibz point'
        ! integer weights
        allocate(w(ibz%nkpts))
        ! get kpoint weights (number of unique points onto which symmetries can map ibz points)
        do i = 1, ibz%nkpts
            w(i) = count(unique_inds(PM(ibz%fbz_id(i),:)).ne.0)
        enddo
        ! get normalized weights
        allocate(ibz%w,source=w/real(fbz%nkpts,dp))
        ! check weights
        if (abs(sum(ibz%w)-1.0_dp).gt.tiny) stop 'ERROR [get_ibz]: kpoint weights does not sum to 1'
        ! stdout
        if (opts%verbosity.ge.1) then
            !
            write(*,'(a,a)') flare, 'FBZ k-points (full monkhorst-pack mesh) = '//tostring(fbz%nkpts)
            write(*,'(a,a)') flare, 'IBZ k-points = '//tostring(ibz%nkpts)
            call disp_indent()
            call disp(title='#'       ,style='underline',X=[1:ibz%nkpts]          ,fmt='i5'   ,advance='no')
            call disp(title='[cart.]' ,style='underline',X=transpose(ibz%kpt_cart),fmt='f10.5',advance='no')
            call disp(title='[frac.]' ,style='underline',X=transpose(ibz%kpt_recp),fmt='f10.5',advance='no')
            call disp(title='w [int]' ,style='underline',X=              w        ,fmt='i5'   ,advance='no')
            call disp(title='w [norm]',style='underline',X=          ibz%w        ,fmt='f10.5',advance='yes')
            ! omit printing permutation map. Way too long for large kpoints
            write(*,'(a,a)',advance='no') flare, 'k-point map (to IBZ: full->irr)'
            call id_print_map(fbz%ibz_id)
            write(*,'(a,a)',advance='no') flare, 'k-point map (from IBZ: irr->full)'
            call id_print_map(ibz%fbz_id)
        endif
        ! now get tetrahedra
        if (opts%verbosity.ge.1) call print_title('Irreducible Tetrahedra')
        ! allocate space
        allocate(ibz_tet(4,fbz%ntets))
        ! sort to get irreducible tetrahedra
        do i = 1, fbz%ntets
            ibz_tet(:,i) = foursort( fbz%ibz_id(fbz%tet(:,i)) )
        enddo
        ! get unique tetrahedra
        ibz_tet_inds = trim_null(unique_inds(ibz_tet))
        ! get number of irreducible tetrahedra
        ibz%ntets = size(ibz_tet_inds)
        ! clear
        if (allocated(w)) deallocate(w)
        allocate(w(ibz%ntets))
        ! allocate stuff
        allocate(ibz%tet(4,ibz%ntets))
        allocate(ibz%tetv(ibz%ntets))
        allocate(ibz%tetw(ibz%ntets))
        ! get tetrahedra weight
        do i = 1, ibz%ntets
            ibz%tet(:,i) = ibz_tet(:,ibz_tet_inds(i))
            j = 0
            do k = 1, fbz%ntets
                if (isequal(ibz%tet(:,i),ibz_tet(:,k))) then
                    j=j+1
                endif
            enddo
            ibz%tetw(i) = j * fbz%tetw(ibz_tet_inds(i))
            ibz%tetv(i) =     fbz%tetv(ibz_tet_inds(i))
            ! save integer
            w(i) = j
        enddo
        ! print stdout
        if (opts%verbosity.ge.1) then
            ! tetrahedra properties
            write(*,*) flare, 'irreducible tetrahedra = '//tostring(ibz%ntets)
            write(*,*) flare, 'tetrahedron volume = '//tostring(ibz%tetv(1))
            call disp_indent()
            call disp(title='#'                ,style='underline',X=[1:ibz%ntets]      ,fmt='i5'   ,advance='no')
            call disp(title='connectivity list',style='underline',X=transpose(ibz%tet) ,fmt='i5'   ,advance='no')
            call disp(title='w [int]'          ,style='underline',X=              w    ,fmt='i5'   ,advance='no')
            call disp(title='w [norm]'         ,style='underline',X=          ibz%tetw ,fmt='f10.5',advance='yes')
        endif
        !
    end subroutine get_irreducible

    ! path

    subroutine     get_path(bzp,path_symb,ndivs,pc,opts)
        !
        implicit none
        !
        class(am_class_path)    , intent(out) :: bzp
        character(*)            , intent(in) :: path_symb
        type(am_class_prim_cell), intent(in) :: pc 
        integer                 , intent(in) :: ndivs ! divisions between kpoint pairs
        type(am_class_options)  , intent(in) :: opts
        character(:), allocatable :: kpt_symb(:)
        character(10000) :: kpt_symb_SE(2)
        real(dp) :: kpt_S(3)
        real(dp) :: kpt_E(3)
        real(dp) :: kpt_D(3)
        real(dp) :: d_SE ! distance from start to end
        integer :: j,k,z,x
        !
        if (opts%verbosity.ge.1) call print_title('Brillouin-zone Path')
        ! transfer bas
        bzp%bas = pc%bas
        ! transfer recbas
        bzp%recbas = pc%recbas
        ! split segment (G-W-K-G-L-U-W-L-K|U-X) => G, W, K, G, L, U, W, L, K|U, X
        kpt_symb = strsplit(path_symb,delimiter='-')
        ! get number of ticks 
        bzp%nticks = size(kpt_symb)
        ! allocate space for tick marks
        allocate(bzp%tick(bzp%nticks))
        ! transfer ndivs
        bzp%ndivs = ndivs
        ! total number of kpoints
        bzp%nkpts = (bzp%nticks-1)*bzp%ndivs
        ! allocate space for kpoints
        allocate(bzp%x(bzp%nkpts))
        allocate(bzp%kpt_recp(3,bzp%nkpts))
        ! loop over directions
        z = 0
        do j = 1, bzp%nticks-1
            ! check if it is a discontinuous jump
            x = index(kpt_symb(j),'|')
            if (x.ne.0) then; kpt_symb_SE(1) = trim(kpt_symb(j)((x+1):))
            else;             kpt_symb_SE(1) = trim(kpt_symb(j))
            endif
            ! check if it is a discontinuous jump
            x = index(kpt_symb(j+1),'|')
            if (x.ne.0) then; kpt_symb_SE(2) = trim(kpt_symb(j+1)(1:(x-1)))
            else;             kpt_symb_SE(2) = trim(kpt_symb(j+1))
            endif
            ! get path start and end coordinates
            kpt_S = symb2kpt(lattice_id=pc%lattice_id, kpt_symb=kpt_symb_SE(1))
            kpt_E = symb2kpt(lattice_id=pc%lattice_id, kpt_symb=kpt_symb_SE(2))
            kpt_D = kpt_E-kpt_S
            ! get normal distance along path
            d_SE = sqrt(sum(kpt_D**2))/bzp%ndivs
            ! increment along path to get x axis
            do k = 1, bzp%ndivs
                z = z + 1
                if (z.eq.1) then
                    bzp%x(z) = 0
                else
                    bzp%x(z) = bzp%x(z-1) + d_SE
                endif
                ! save tick location
                if (k.eq.1) bzp%tick(j) = bzp%x(z)
                ! save kpt_recp
                bzp%kpt_recp(:,z) = kpt_S + (k-1)/real(bzp%ndivs,dp) * kpt_D
            enddo
        enddo
        ! save last tick
        bzp%tick(j) = bzp%x(z)
        ! save kpt_symb, stupid fortran (can not use bzp%kpt_symb throughout beacuse of mysterious compiler problem)
        allocate(bzp%kpt_symb,source=kpt_symb)
        ! get kpt_cart
        allocate(bzp%kpt_cart,source=matmul(bzp%recbas,bzp%kpt_recp))
        ! stdout
        if (opts%verbosity.ge.1) then
            ! kpoint properties
            write(*,'(a,a)') flare, 'path = '//path_symb
            write(*,'(a,a)') flare, 'kpoints (per direction) = '//tostring(bzp%ndivs)
            write(*,'(a,a)') flare, 'kpoints (total) = '//tostring(bzp%nkpts)
            call disp_indent()
            call disp(title='#'      ,style='underline',X=[1:bzp%nkpts]          ,fmt='i5'   ,advance='no')
            call disp(title='x'      ,style='underline',X=          bzp%x        ,fmt='f10.5',advance='no')
            call disp(title='[cart.]',style='underline',X=transpose(bzp%kpt_cart),fmt='f10.5',advance='no')
            call disp(title='[frac.]',style='underline',X=transpose(bzp%kpt_recp),fmt='f10.5',advance='yes')
        endif
    end subroutine get_path

    function       symb2kpt(lattice_id,kpt_symb) result(kpt)
        !
        implicit none
        !
        integer, intent(in) :: lattice_id
        character(*), intent(in) :: kpt_symb
        real(dp) :: kpt(3)
        !
        select case(lattice_id)
        case(1) ! simple cubic
            select case(trim(kpt_symb))
            case('G'); kpt = [ 0.000_dp, 0.000_dp, 0.000_dp]
            case('M'); kpt = [ 0.500_dp, 0.500_dp, 0.000_dp]
            case('R'); kpt = [ 0.500_dp, 0.500_dp, 0.500_dp]
            case('X'); kpt = [ 0.000_dp, 0.500_dp, 0.000_dp]
            case default
            stop 'ERROR [symb2kpt]: unknown symbol (lattice_id = 1)'
            end select
        case(2) ! bcc
            select case(trim(kpt_symb))
            case('G'); kpt = [ 0.000_dp, 0.000_dp, 0.000_dp]
            case('H'); kpt = [ 0.500_dp,-0.500_dp, 0.500_dp]
            case('P'); kpt = [ 0.250_dp, 0.250_dp, 0.250_dp]
            case('N'); kpt = [ 0.000_dp, 0.000_dp, 0.500_dp]
            case default
            stop 'ERROR [symb2kpt]: unknown symbol (lattice_id = 2)'
            end select
        case(3) ! fcc
            select case(trim(kpt_symb))
            case('G'); kpt = [ 0.000_dp, 0.000_dp, 0.000_dp]
            case('K'); kpt = [ 0.375_dp, 0.375_dp, 0.750_dp]
            case('L'); kpt = [ 0.500_dp, 0.500_dp, 0.500_dp]
            case('U'); kpt = [ 0.625_dp, 0.250_dp, 0.625_dp]
            case('W'); kpt = [ 0.500_dp, 0.250_dp, 0.750_dp]
            case('X'); kpt = [ 0.500_dp, 0.000_dp, 0.500_dp]
            case default
            stop 'ERROR [symb2kpt]: unknown symbol (lattice_id = 3)'
            end select
        case(4) ! hexagonal
            select case(trim(kpt_symb))
            case('G'); kpt = [ 0.000_dp, 0.000_dp, 0.000_dp]
            case('A'); kpt = [ 0.000_dp, 0.000_dp, 0.500_dp]
            case('H'); kpt = [ 1/3.0_dp, 1/3.0_dp, 0.500_dp]
            case('K'); kpt = [ 1/3.0_dp, 1/3.0_dp, 0.000_dp]
            case('L'); kpt = [ 0.500_dp, 0.000_dp, 0.500_dp]
            case('M'); kpt = [ 0.500_dp, 0.000_dp, 0.000_dp]
            case default
            stop 'ERROR [symb2kpt]: unknown symbol (lattice_id = 4)'
            end select
        case(5) ! tetragonal (check if this is correct?)
            select case(trim(kpt_symb))
            case('G'); kpt = [ 0.000_dp, 0.000_dp, 0.000_dp]
            case('A'); kpt = [ 0.500_dp, 0.500_dp, 0.500_dp]
            case('M'); kpt = [ 0.500_dp, 0.500_dp, 0.000_dp]
            case('R'); kpt = [ 0.000_dp, 0.500_dp, 0.500_dp]
            case('X'); kpt = [ 0.000_dp, 0.500_dp, 0.000_dp]
            case('Z'); kpt = [ 0.000_dp, 0.000_dp, 0.500_dp]
            case default
            stop 'ERROR [symb2kpt]: unknown symbol (lattice_id = 5)'
            end select
        case(13) ! body-centered tetragonal
            select case(trim(kpt_symb))
            case('G'); kpt = [ 0.000_dp, 0.000_dp, 0.000_dp]
            case('M'); kpt = [ 0.500_dp, 0.500_dp,-0.500_dp]
            case('X'); kpt = [ 0.000_dp, 0.000_dp, 0.500_dp]
            case('P'); kpt = [ 0.250_dp, 0.250_dp, 0.250_dp]
            case('N'); kpt = [ 0.000_dp, 0.500_dp, 0.000_dp]
            case default
            stop 'ERROR [symb2kpt]: unknown symbol (lattice_id = 13)'
            end select
        case(7) ! trigonal (maybe this is triclinic?)
            select case(trim(kpt_symb))
            case('G'); kpt = [ 0.000_dp, 0.000_dp, 0.000_dp]
            case('L'); kpt = [ 0.500_dp,-0.500_dp, 0.000_dp]
            case('M'); kpt = [ 0.000_dp, 0.000_dp, 0.500_dp]
            case('N'); kpt = [-0.500_dp,-0.500_dp, 0.500_dp]
            case('R'); kpt = [ 0.000_dp,-0.500_dp,-0.500_dp]
            case('X'); kpt = [-0.500_dp, 0.000_dp, 0.000_dp]
            case('Y'); kpt = [ 0.500_dp, 0.000_dp, 0.000_dp]
            case('Z'); kpt = [-0.500_dp, 0.000_dp, 0.500_dp]
            case default
            stop 'ERROR [symb2kpt]: unknown symbol (lattice_id = 7)'
            end select
        case(8) ! orthorhombic (maybe this is not simple orthorhomnic but has another centering?)
            select case(trim(kpt_symb))
            case('G'); kpt = [ 0.000_dp, 0.000_dp, 0.000_dp]
            case('R'); kpt = [ 0.500_dp, 0.500_dp, 0.500_dp]
            case('S'); kpt = [ 0.500_dp, 0.500_dp, 0.000_dp]
            case('T'); kpt = [ 0.000_dp, 0.500_dp, 0.500_dp]
            case('U'); kpt = [ 0.500_dp, 0.000_dp, 0.500_dp]
            case('X'); kpt = [ 0.500_dp, 0.000_dp, 0.000_dp]
            case('Y'); kpt = [ 0.000_dp, 0.500_dp, 0.000_dp]
            case('Z'); kpt = [ 0.000_dp, 0.000_dp, 0.500_dp]
            case default
            stop 'ERROR [symb2kpt]: unknown symbol (lattice_id = 8)'
            end select
        case default
            stop 'ERROR [symb2kpt]: unknown lattice'
        end select
    end function   symb2kpt

    function       get_pathsymb(lattice_id) result(path_symb)
        !
        implicit none
        !
        integer, intent(in) :: lattice_id
        character(:), allocatable :: path_symb
        ! get path
        select case(lattice_id)
        case(1);  allocate(path_symb, source='G-M-G-R-X|M-R') ! cubic
        case(2);  allocate(path_symb, source='G-N-G-P-H|P-N') ! bcc
        case(3);  allocate(path_symb, source='G-W-K-G-L-U-W-L-K|U-X') ! fcc
        case(4);  allocate(path_symb, source='G-K-G-A-L-H|H-A|L-M|K-H') ! hexagonal
        case(5);  allocate(path_symb, source='G-M-G-Z-R-A-Z|X-R|M-A') ! tetragonal
        case(8);  allocate(path_symb, source='G-S-Y-G-Z-U-R-T|Z-T|Y-T|U-X|S-R') ! orthorhombic
        case(12); allocate(path_symb, source='G-H-C-E-M1-A-X-H1|M-D-Z|Y-D') ! monoclinic
        case default 
            stop 'ERROR [get_standard_path]: unknown lattice_id (probably not yet implemented)'
        end select
    end function   get_pathsymb

    ! tetrahedra functions

    function       get_tetrahedra(recbas,kpt_int,n) result(tet)
        !
        implicit none
        !
        real(dp),intent(in) :: recbas(3,3)
        integer, intent(in) :: kpt_int(:,:)
        integer, intent(in) :: n(3)
        integer, allocatable :: tet(:,:)
        integer, allocatable :: box(:,:) ! box(8,nboxes) vertices of boxes index kpoint in fbz 
        integer :: tetrahedron(4,6)      ! tetrahedron(nvertices,ntets) four verticies of the six tetrahedra making up one box
        integer :: nboxes
        integer :: t ! tetrahedron counter
        integer :: i,k
        ! divide mesh into boxes
        box = grid2box(kpt_int,n)
        ! divide a single box into six tetrahedron
        tetrahedron = box2tetrahedron(recbas=recbas)
        ! get number of boxes
        nboxes = size(box,2)
        ! allocate space for tetrahedra
        allocate(tet(4,6*nboxes))
        tet = 0
        ! initialize tetrahedron counter
        t = 0
        ! loop over boxes
        do i = 1, nboxes
            ! loop over tetrahedron/box
            do k = 1, 6
                ! augment tetrahedron counter
                t = t + 1
                ! define tetrahedra corners using indices of kpoints
                tet(:,t) = box(tetrahedron(:,k),i)
            enddo
        enddo
        !
    end function   get_tetrahedra

    function       grid2box(kpt_int,n) result(box)
        !
        implicit none
        !
        integer, intent(in) :: kpt_int(:,:)
        integer, intent(in) :: n(3)
        integer, allocatable :: box(:,:)
        integer :: nkpts
        integer :: sub(3)
        integer :: i,j
        integer :: ii,jj,kk
        ! get nkpts
        nkpts = size(kpt_int,2)
        ! there will be 1 box per kpoint and 8 vertices per box
        allocate(box(8,nkpts))
        ! loop over kpoints
        do i = 1, nkpts
            j=0
            ! m,n,o loops over 8 vertices
            do ii = 0,1
            do jj = 0,1
            do kk = 0,1
                j=j+1
                ! modulo applies periodic boundary conditions
                sub = [modulo(kpt_int(1,i)+ii,n(1)),modulo(kpt_int(2,i)+jj,n(2)),modulo(kpt_int(3,i)+kk,n(3))]
                ! Q: what is the index of the kpoint at this grid location?
                ! add [1,1,1] to account for the fact that kpoints start at [0,0,0] rather than [1,1,1]
                box(j,i) = sub2ind(dims=n,sub=sub+1)
            enddo
            enddo
            enddo
        enddo
    end function   grid2box

    function       box2tetrahedron(recbas) result(tetrahedron)
        !     7-------8
        !    /|      /|
        !   / |     / |
        !  5-------6  |
        !  |  3----|--4
        !  | /     | /
        !  |/      |/
        !  1-------2
        implicit none
        !
        real(dp), intent(in) :: recbas(3,3)
        integer  :: tetrahedron(4,6) ! tetrahedra indices
        integer  :: diags(2,4)  ! diagonal pairs
        real(dp) :: d(8)     ! diagonal pair distances 
        real(dp) :: sd       ! shortest diagonal pair
        integer  :: si       ! shortest diagonal pair index
        integer  :: i,j,k,l
        real(dp) :: box(3,8)
        ! get a box
        l=0
        do i=0,1
        do j=0,1
        do k=0,1
            l=l+1
            box(:,l)=matmul(recbas,[i,j,k])
        enddo
        enddo
        enddo
        ! get indices of diagonal pairs
        diags(:,1) = [1,8]
        diags(:,2) = [2,7]
        diags(:,3) = [3,6]
        diags(:,4) = [4,5]
        ! get distances across diagonals
        do i = 1, 4
            d(i) = norm2(box(:,diags(2,i))-box(:,diags(1,i)))
        enddo
        ! record smallest diagonal
        sd = 1.0D8
        ! initialize si
        si = 0
        ! loop over diagonals
        do i = 1, 4
        if (d(i).lt.sd) then
            ! update smallest diagonal distance
            sd = d(i)
            ! record index of smallest diagonal
            si = i
        endif
        enddo
        ! check
        if (si.eq.0) stop 'ERROR [divide_box_into_tetrahedra]: failed to get smallest diagonal'
        ! create connectivity list defining tetrahedra
        select case (si)
        case (1)
        tetrahedron(:,1) = [1,8,2,4]
        tetrahedron(:,2) = [1,8,2,6]
        tetrahedron(:,3) = [1,8,3,4]
        tetrahedron(:,4) = [1,8,3,7]
        tetrahedron(:,5) = [1,8,5,6]
        tetrahedron(:,6) = [1,8,5,7]
        case (2)
        tetrahedron(:,1) = [2,7,1,3]
        tetrahedron(:,2) = [2,7,1,5]
        tetrahedron(:,3) = [2,7,3,4]
        tetrahedron(:,4) = [2,7,4,8]
        tetrahedron(:,5) = [2,7,5,6]
        tetrahedron(:,6) = [2,7,6,8]
        case (3)
        tetrahedron(:,1) = [3,6,1,2]
        tetrahedron(:,2) = [3,6,1,5]
        tetrahedron(:,3) = [3,6,2,4]
        tetrahedron(:,4) = [3,6,4,8]
        tetrahedron(:,5) = [3,6,5,7]
        tetrahedron(:,6) = [3,6,7,8]
        case (4)
        tetrahedron(:,1) = [4,5,1,2]
        tetrahedron(:,2) = [4,5,1,3]
        tetrahedron(:,3) = [4,5,2,6]
        tetrahedron(:,4) = [4,5,3,7]
        tetrahedron(:,5) = [4,5,6,8]
        tetrahedron(:,6) = [4,5,7,8]
        end select
    end function   box2tetrahedron

end module am_brillouin_zone



