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
        real(dp), allocatable :: w(:)           ! w(nkpts) normalized weights
        ! get_irreducible
        integer , allocatable :: fbz_id(:)
        integer , allocatable :: ibz_id(:)
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
        integer  :: ntets                    ! number of tetrahedra
        integer , allocatable :: kpt_int(:,:)  ! kpt(3,nkpts) kpoint integer indices in regular grid
        integer , allocatable :: tet(:,:)   ! tetrahedra connectivity list (indices of corner kpoints)
        real(dp), allocatable :: tetv(:)    ! tetrahedra volume
        real(dp), allocatable :: tetw(:)    ! tetrahedron weight
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
    end type am_class_path

    ! other stuff

    public :: get_tetrahedron_weight

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
                #:for ATTRIBUTE in ['kpt_cart','kpt_recp','w','fbz_id','ibz_id']
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
                    #:for ATTRIBUTE in ['kpt_int','tet','tetv','tetw']
                        $:write_xml_attribute_allocatable(ATTRIBUTE)
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
        !
        if (iotype.eq.'LISTDIRECTED') then
            read(unit=unit,fmt='(/)',iostat=iostat,iomsg=iomsg)
                ! non-allocatable
                #:for ATTRIBUTE in ['bas','recbas','nkpts']
                    $:read_xml_attribute_nonallocatable(ATTRIBUTE)
                #:endfor
                ! allocatable
                #:for ATTRIBUTE in ['kpt_cart','kpt_recp','w','fbz_id','ibz_id']
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
                    #:for ATTRIBUTE in ['kpt_int','tet','tetv','tetw']
                        $:read_xml_attribute_allocatable(ATTRIBUTE)
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

    subroutine     create_bz(bz,kpt_recp,kpt_cart,bas,prec,w)
        !
        ! kpt_recp points are shifted to be between [0,1)
        ! kpt_cart        are shifted to be in the wigner-seitz voronoi cell
        !
        implicit none
        !
        class(am_class_bz), intent(out) :: bz
        real(dp), intent(in), optional :: kpt_recp(:,:)
        real(dp), intent(in), optional :: kpt_cart(:,:)
        real(dp), intent(in), optional :: w(:)
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
        ! weights
        if (present(w)) allocate(bz%w, source = w)
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
        kpt_int = generate_integer_monkhorst_pack_mesh(n=fbz%n,s=fbz%s)
        ! nkpts
        nkpts = size(kpt_int,2)
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
            call disp_indent()
            call disp(title='#'                ,style='underline',X=[1:fbz%ntets]     ,fmt='i5' ,advance='no')
            call disp(title='connectivity list',style='underline',X=transpose(fbz%tet),fmt='i10',advance='yes')
        endif
        !
        contains
        pure function  generate_integer_monkhorst_pack_mesh(n,s) result(kpt)
            !> returns kpoints in fractional coordinates for monkhorst pack mesh dimensions
            !> n(1:3)=[n1,n2,n3] and shift s(1:3)=[s1,s2,s3]
            !> according to http://cms.mpi.univie.ac.at/vasp/vasp/Automatic_k_mesh_generation.html
            implicit none
            !
            integer,  intent(in) :: n(3) !> monkhorst pack mesh dimensions
            real(dp), intent(in) :: s(3) !> monkhorst pack mesh shift
            integer, allocatable :: kpt(:,:) !> kpt in fractional
            integer :: i1,i2,i3,j
            ! create array in column-major order, essential because sub2ind is used later on
            allocate(kpt(3,product(n)))
            j=0
            do i3=1,n(3) ! z first
            do i2=1,n(2) ! y second
            do i1=1,n(1) ! x last
                j=j+1
                kpt(1:3,j) = [i1,i2,i3]-1
            enddo
            enddo
            enddo
            !
        end function   generate_integer_monkhorst_pack_mesh
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
        ! get irreducible tetrahedra
        do i = 1, fbz%ntets
            ibz_tet(:,i) = fbz%ibz_id(fbz%tet(:,i))
            ibz_tet(:,i) = ibz_tet(sort_four_integers(ibz_tet(:,i)),i)
        enddo
        ! get unique tetrahedra
        ibz_tet_inds = trim_null(unique_inds(ibz_tet))
        ! get number of irreducible tetrahedra
        ibz%ntets = size(ibz_tet_inds)
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
        enddo
        ! print stdout
        if (opts%verbosity.ge.1) then
            ! tetrahedra properties
            write(*,*) flare, 'tetrahedra = '//tostring(ibz%ntets)
            write(*,*) flare, 'tetrahedra volume = '//tostring(ibz%tetv(1))
            call disp_indent()
            call disp(title='#'                ,style='underline',X=[1:ibz%ntets]      ,fmt='i5'  ,advance='no')
            call disp(title='connectivity list',style='underline',X=transpose(ibz%tet) ,fmt='i5'  ,advance='no')
            call disp(title='w'                ,style='underline',X=          ibz%tetw ,fmt='f10.5',advance='yes')
        endif
        !
    end subroutine get_irreducible

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
        integer :: i,j,k
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
        ! divide by four because ... ???
        w=w*0.25_dp
        !
    end function   get_tetrahedron_weight







!     subroutine     get_path()
!         !
!         implicit none
!         !



! !     function     standardpath(bravais) result(points)
! !         !
! !         integer, intent(in) :: bravais
! !         character(len=3), allocatable :: points(:,:)
! !         !
! !         select case(trim(ucposcar%bravaislattice))
! !             case('CUB')
! !                 allocate(points(5,2))            
! !                 points(1,:) = ['GM','X ']
! !                 points(2,:) = ['X ','M ']
! !                 points(3,:) = ['M ','GM']
! !                 points(4,:) = ['GM','R ']
! !                 points(5,:) = ['R ','X ']
! !             case('FCC')
! !                 allocate(points(9,2))
! !                 points(1,:) = ['GM','X ']
! !                 points(2,:) = ['X ','W ']
! !                 points(3,:) = ['W ','K ']
! !                 points(4,:) = ['K ','GM']
! !                 points(5,:) = ['GM','L ']
! !                 points(6,:) = ['L ','U ']
! !                 points(7,:) = ['U ','W ']
! !                 points(8,:) = ['W ','L ']
! !                 points(9,:) = ['L ','K ']
! !             case('BCC')
! !                 allocate(points(5,2))
! !                 points(1,:) = ['GM','H ']
! !                 points(2,:) = ['H ','N ']
! !                 points(3,:) = ['N ','GM']
! !                 points(4,:) = ['GM','P ']
! !                 points(5,:) = ['P ','H ']
! !             case('HEX')
! !                allocate(points(7,2))
! !                 points(1,:) = ['GM','M ']
! !                 points(2,:) = ['M ','K ']
! !                 points(3,:) = ['K ','GM']
! !                 points(4,:) = ['GM','A ']
! !                 points(5,:) = ['A ','L ']
! !                 points(6,:) = ['L ','H ']
! !                 points(7,:) = ['H ','A ']
! !             case('RHL1')
! !                 allocate(points(9,2))
! !                 points(1,:) = ['GM','L ']
! !                 points(2,:) = ['L ','B1']
! !                 points(3,:) = ['B ','Z ']
! !                 points(4,:) = ['Z ','GM']
! !                 points(5,:) = ['GM','X ']
! !                 points(6,:) = ['Q ','F ']
! !                 points(7,:) = ['F ','P1']
! !                 points(8,:) = ['P1','Z ']
! !                 points(9,:) = ['L ','P ']            
! !             case('TRI')
! !                 allocate(points(7,2))
! !                 points(1,:) = ['X ','GM']
! !                 points(2,:) = ['GM','Y ']
! !                 points(3,:) = ['L ','GM']
! !                 points(4,:) = ['GM','Z ']
! !                 points(5,:) = ['N ','GM']
! !                 points(6,:) = ['GM','M ']
! !                 points(7,:) = ['R ','GM']
! !             case('ORC')
! !                 allocate(points(7,2))
! !                 points(1,:) = ['X ','GM']
! !                 points(2,:) = ['GM','Y ']
! !                 points(3,:) = ['L ','GM']
! !                 points(4,:) = ['GM','Z ']
! !                 points(5,:) = ['N ','GM']
! !                 points(6,:) = ['GM','M ']
! !                 points(7,:) = ['R ','GM']
! !             case('TET')
! !                allocate(points(9,2))
! !                 points(1,:) = ['GM','X ']
! !                 points(2,:) = ['X ','M ']
! !                 points(3,:) = ['M ','GM']
! !                 points(4,:) = ['GM','Z ']
! !                 points(5,:) = ['Z ','R ']
! !                 points(6,:) = ['R ','A ']
! !                 points(7,:) = ['A ','Z ']
! !                 points(8,:) = ['X ','R ']
! !                 points(9,:) = ['M ','A ']
! !             case default
! !                 stop 'could not find path for Bravais lattice'
! !         end select
! !         !
! !     end function standardpath

! !     subroutine lo_get_coordinates_from_symbolic_point(symbol,qv,ucposcar)
! !         !
! !         ! I have not bothered to add all Bravais lattices yet. I will fix that in due time.
! !         !
! !         character(len=*), intent(in) :: symbol
! !         real(dp), dimension(3), intent(out) :: qv
! !         type(lo_crystalstructure), intent(in) :: ucposcar
! !         !
! !         integer :: i
! !         real(dp), dimension(3) :: v
! !         real(dp) :: alpha,beta,gm,a,b,c
! !         real(dp) :: w,x,y,z
! !         !
! !         ! Get some parameters from the poscar:
! !         !
! !         a=lo_norm(ucposcar%basis(1,:))*ucposcar%latpar
! !         b=lo_norm(ucposcar%basis(2,:))*ucposcar%latpar
! !         c=lo_norm(ucposcar%basis(3,:))*ucposcar%latpar
! !         !
! !         alpha=lo_angle_between_vectors(ucposcar%basis(1,:),ucposcar%basis(2,:))*180/lo_pi
! !         beta=lo_angle_between_vectors(ucposcar%basis(2,:),ucposcar%basis(3,:))*180/lo_pi
! !         gm=lo_angle_between_vectors(ucposcar%basis(1,:),ucposcar%basis(3,:))*180/lo_pi
! !         !
! !         v=0.0_dp
! !         select case(trim(ucposcar%bravaislattice))
! !             case('CUB')
! !                 select case(trim(symbol))
! !                     case('GM'); v = [ 0.0_dp, 0.0_dp, 0.0_dp]
! !                     case('R');  v = [ 0.5_dp, 0.5_dp, 0.5_dp]
! !                     case('M');  v = [ 0.5_dp, 0.5_dp, 0.0_dp]
! !                     case('X');  v = [ 0.0_dp, 0.5_dp, 0.0_dp]
! !                 end select
! !             case('FCC')
! !                 select case(trim(symbol))
! !                 w = 0.25_dp
! !                 x = 3.0_dp/4.0_dp
! !                 y = 3.0_dp/8.0_dp
! !                 z = 5.0_dp/8.0_dp
! !                     case('GM'); v = [ 0.0_dp, 0.0_dp, 0.0_dp]
! !                     case('X');  v = [ 0.0_dp, 0.5_dp, 0.5_dp]
! !                     case('X2'); v = [ 0.5_dp, 0.5_dp, 1.0_dp]
! !                     case('U');  v = [ 0.0_dp, z     , y     ]
! !                     case('L');  v = [ 0.0_dp, 0.5_dp, 0.0_dp]
! !                     case('W');  v = [ w     , x     , 0.5_dp]
! !                     case('K');  v = [ y     , x     , y     ]
! !                 end select
! !             case('BCC')
! !                 select case(trim(symbol))
! !                     y = 0.25_dp
! !                     case('GM'); v = [ 0.0_dp, 0.0_dp, 0.0_dp]
! !                     case('H');  v = [ 0.5_dp,-0.5_dp, 0.5_dp]
! !                     case('P');  v = [ y     , y     , y     ]
! !                     case('N');  v = [ 0.0_dp, 0.0_dp, 0.5_dp]
! !                 end select
! !             case('HEX')
! !                 select case(trim(symbol))
! !                     y = 1.0_dp/3.0_dp
! !                     case('GM'); v = [ 0.0_dp, 0.0_dp, 0.0_dp]
! !                     case('A');  v = [ 0.0_dp, 0.0_dp, 0.5_dp]
! !                     case('K');  v = [ y     , y     , 0.0_dp]
! !                     case('H');  v = [ y     , y     , 0.5_dp]
! !                     case('M');  v = [ 0.5_dp, 0.0_dp, 0.0_dp]
! !                     case('L');  v = [ 0.5_dp, 0.0_dp, 0.5_dp]
! !                 end select
! !             case('BCT')
! !                 select case(trim(symbol))
! !                     y  = 0.25_dp
! !                     case('GM'); v = [ 0.0_dp, 0.0_dp, 0.0_dp]
! !                     case('M');  v = [ 0.5_dp, 0.5_dp,-0.5_dp]
! !                     case('X');  v = [ 0.0_dp, 0.0_dp, 0.5_dp]
! !                     case('P');  v = [ y     , y     , y     ]
! !                     case('N');  v = [ 0.0_dp, 0.5_dp, 0.0_dp]
! !                 end select
! !             case('TET') ! normal tetragonal
! !                 select case(trim(symbol))
! !                     case('GM'); v = [ 0.0_dp, 0.0_dp, 0.0_dp]
! !                     case('A');  v = [ 0.5_dp, 0.5_dp, 0.5_dp]
! !                     case('M');  v = [ 0.5_dp, 0.5_dp, 0.0_dp]
! !                     case('R');  v = [ 0.0_dp, 0.5_dp, 0.5_dp]
! !                     case('X');  v = [ 0.0_dp, 0.5_dp, 0.0_dp]
! !                     case('Z');  v = [ 0.0_dp, 0.0_dp, 0.5_dp]
! !                 end select
! !             case('RHL1')
! !                 x = (1.0_dp+4.0_dp*cos(alpha*pi/180.0_dp))/(2.0_dp+4.0_dp*cos(alpha*pi/180.0_dp))
! !                 y = 3.0_dp/4-x  /2
! !                 select case(trim(symbol))
! !                     case('GM'); v = [ 0.0_dp, 0.0_dp, 0.0_dp] ! Gamma
! !                     case('B');  v = [ x     , 0.5_dp, 1-x   ] ! B
! !                     case('B1'); v = [ 0.5_dp, 1-x   , x  -1 ] ! B1
! !                     case('F');  v = [ 0.5_dp, 0.5_dp, 0.0_dp] ! F
! !                     case('L');  v = [ 0.5_dp, 0.0_dp, 0.0_dp] ! L
! !                     case('L1'); v = [ 0.0_dp, 0.0_dp,-0.5_dp] ! L1
! !                     case('P');  v = [ x     , y     , y     ] ! P
! !                     case('P1'); v = [ 1-y   , 1-y   , 1-x   ] ! P1
! !                     case('P2'); v = [ y     , y     , x  -1 ] ! P2
! !                     case('X');  v = [ y     , 0.0_dp,-y     ] ! X
! !                     case('Z');  v = [ 0.5_dp, 0.5_dp, 0.5_dp] ! Z
! !                 end select
! !             case('TRI')
! !                 select case(trim(symbol))
! !                     case('GM'); v = [ 0.0_dp, 0.0_dp, 0.0_dp]
! !                     case('L');  v = [ 0.5_dp,-0.5_dp, 0.0_dp]
! !                     case('M');  v = [ 0.0_dp, 0.0_dp, 0.5_dp]
! !                     case('N');  v = [-0.5_dp,-0.5_dp, 0.5_dp]
! !                     case('R');  v = [ 0.0_dp,-0.5_dp,-0.5_dp]
! !                     case('X');  v = [-0.5_dp, 0.0_dp, 0.0_dp]
! !                     case('Y');  v = [ 0.5_dp, 0.0_dp, 0.0_dp]
! !                     case('Z');  v = [-0.5_dp, 0.0_dp, 0.5_dp]
! !                 end select
! !             case('ORC')
! !                 select case(trim(symbol))
! !                     case('GM'); v = [ 0.0_dp, 0.0_dp, 0.0_dp]
! !                     case('R');  v = [ 0.5_dp, 0.5_dp, 0.5_dp]
! !                     case('S');  v = [ 0.5_dp, 0.5_dp, 0.0_dp]
! !                     case('T');  v = [ 0.0_dp, 0.5_dp, 0.5_dp]
! !                     case('U');  v = [ 0.5_dp, 0.0_dp, 0.5_dp]
! !                     case('X');  v = [ 0.5_dp, 0.0_dp, 0.0_dp]
! !                     case('Y');  v = [ 0.0_dp, 0.5_dp, 0.0_dp]
! !                     case('Z');  v = [ 0.0_dp, 0.0_dp, 0.5_dp]
! !                 end select
! !         end select
! !         qv=0.0_dp
! !         do i=1,3
! !             qv=qv+v(i)*ucposcar%recbasis(i,:)
! !         enddo
! !     end subroutine
!     end subroutine get_path


!     subroutine     get_fbz(fbz,bz,pc,pg,opts)
!         !
!         implicit none
!         !
!         class(am_class_fbz), intent(out) :: fbz
!         class(am_class_bz) , intent(in)  :: bz
!         type(am_class_prim_cell), intent(in) :: pc
!         type(am_class_seitz_group) , intent(in) :: pg
!         type(am_class_options)  , intent(in) :: opts
!         real(dp), allocatable :: korbit(:,:)
!         real(dp), allocatable :: kpoints_fbz(:,:)
!         real(dp), allocatable :: grid_points(:,:)
!         integer :: i, j, m
!         integer :: w
!         !
!         if (opts%verbosity.ge.1) call print_title('Expand kpoints to FBZ')
!         if (opts%verbosity.ge.1) call am_print('number of point symmetries',pg%nsyms,flare)
!         if (opts%verbosity.ge.1) call am_print('number of kpoints in the original bz',bz%nkpts,flare)
!         if (opts%verbosity.ge.1) then
!             call am_print_two_matrices_side_by_side(name='original kpoints',&
!                 Atitle='fractional',A=transpose(bz%kpt),&
!                 Btitle='cartesian' ,B=transpose(matmul(inv(pc%bas),bz%kpt)),&
!                 iopt_emph=flare,iopt_teaser=.true.)
!         endif
!         !
!         allocate(kpoints_fbz(3,bz%nkpts*pg%nsyms)) ! wkrspace
!         m = 0
!         ! expand each kpoint onto the fbz
!         !$OMP PARALLEL PRIVATE(i,w,korbit) SHARED(kpoints_fbz,pg)
!         !$OMP DO
!         do i = 1, bz%nkpts
!             ! get its orbit
!             korbit = kpoint_orbit(pg%sym(1:3,1:3,:),bz%kpt(:,i))
!             ! reduce korbit to primitive reciprocal lattice
!             korbit = modulo(korbit+opts%prec,1.0_dp)-opts%prec
!             ! get unique values
!             korbit = unique(korbit,opts%prec)
!             ! get its weight
!             w = size(korbit,2)
!             ! check that weight is a factor of the number of symmetry operations
!             if ( modulo(pg%nsyms,w) .ne. 0 ) then
!                 call am_print('ERROR','Weight is not a factor of the number of symmetry operation',flags='E')
!                 call am_print('kpoint',bz%kpt(:,i))
!                 call am_print('weight',w)
!                 call am_print('number of symmetry operations',pg%nsyms)
!                 call am_print('kpoint orbit',transpose(korbit))
!                 stop
!             endif
!             ! save all unique points in orbit as part of the fbz
!             do j = 1,w
!                 m = m + 1
!                 kpoints_fbz(:,m) = korbit(:,j)
!             enddo
!         enddo
!         !$OMP END DO
!         !$OMP END PARALLEL
!         ! make sure there are no repetitions in the fbz
!         fbz%kpt = unique(kpoints_fbz(:,1:m),opts%prec)
!         ! get number of kpoints in the fbz
!         fbz%nkpts = size(fbz%kpt,2)
!         ! reduce kpoints to FBZ using voronoi points
!         grid_points = meshgrid([-1:1],[-1:1],[-1:1])
!         grid_points = matmul(inv(pc%bas),grid_points)
!         !$OMP PARALLEL PRIVATE(i) SHARED(fbz,grid_points,pc)
!         !$OMP DO
!         do i = 1,fbz%nkpts
!             fbz%kpt(:,i) = reduce_kpoint_to_fbz(kpoint=fbz%kpt(:,i),grid_points=grid_points,bas=pc%bas)
!         enddo
!         !$OMP END DO
!         !$OMP END PARALLEL
!         if (opts%verbosity.ge.1) call am_print('number of kpoints in the fbz',fbz%nkpts,flare)
!         if (opts%verbosity.ge.1) then
!             call am_print_two_matrices_side_by_side(name='full kpoints',&
!                 Atitle='fractional',A=transpose(fbz%kpt),&
!                 Btitle='cartesian' ,B=transpose(matmul(inv(pc%bas),fbz%kpt)),&
!                 iopt_emph=flare,iopt_teaser=.true.)
!         endif
!         ! set weights
!         allocate(fbz%w(fbz%nkpts))
!         fbz%w = 1/real(fbz%nkpts,dp)
!         ! set irreducible indices
!         if (allocated(fbz%ibz_id)) deallocate(fbz%ibz_id)
!         allocate(fbz%ibz_id(fbz%nkpts))
!         do i = 1, fbz%nkpts
!             get_indices_of_each_kpoint_in_bz : do j = 1, bz%nkpts
!                 if (all(abs(bz%kpt(:,j)-fbz%kpt(:,i)).lt.tiny)) then
!                     fbz%ibz_id(i) = j
!                     exit get_indices_of_each_kpoint_in_bz
!                 endif
!             enddo get_indices_of_each_kpoint_in_bz
!         enddo
!     end subroutine get_fbz


end module am_brillouin_zone



