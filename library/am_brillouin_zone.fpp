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
        real(dp), allocatable :: kpt_frac(:,:)  ! kpt(3,nkpts) kpoint - notice that vasp IBZKPT is in (reciprocal) fractional units.
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
        integer, allocatable :: kpt_int(:,:)  ! kpt(3,nkpts) kpoint integer indices in regular grid
        ! get_tetrahedra
        integer :: ntetra                   ! number of tetrahedra
        integer , allocatable :: tet(:,:)   ! tetrahedra connectivity list (indices of corner kpoints)
        real(dp), allocatable :: tetv(:)    ! tetrahedra volume
        real(dp), allocatable :: tetw(:)    ! tetrahedron weight
        contains
        procedure :: get_fbz
    end type am_class_fbz

    ! IBZ (irreducible wedge)

    type, public, extends(am_class_bz) :: am_class_ibz
        contains
        procedure :: get_irreducible ! requires fbz monkhorst-pack mesh
    end type am_class_ibz

    ! PATH 

    type, public, extends(am_class_bz) :: am_class_path
    end type am_class_path

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
                #:for ATTRIBUTE in ['kpt_cart','kpt_frac','w','fbz_id','ibz_id']
                    $:write_xml_attribute_allocatable(ATTRIBUTE)
                #:endfor
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
                #:for ATTRIBUTE in ['kpt_cart','kpt_frac','w','fbz_id','ibz_id']
                    $:read_xml_attribute_allocatable(ATTRIBUTE)
                #:endfor
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

    subroutine     create_bz(bz,kpt_frac,kpt_cart,bas,prec,w)
        !
        ! kpt_frac points are shifted to be between [0,1)
        ! kpt_cart        are shifted to be in the wigner-seitz voronoi cell
        !
        implicit none
        !
        class(am_class_bz), intent(out) :: bz
        real(dp), intent(in), optional :: kpt_frac(:,:)
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
        if (present(kpt_frac)) then
            !
            bz%nkpts = size(kpt_frac,2)
            allocate(bz%kpt_cart(3,bz%nkpts))
            allocate(bz%kpt_frac(3,bz%nkpts))
            !
            bz%kpt_frac = kpt_frac
            bz%kpt_cart = matmul(bz%recbas,kpt_frac)
            !
        elseif (present(kpt_cart)) then
            !
            bz%nkpts = size(kpt_cart,2)
            allocate(bz%kpt_cart(3,bz%nkpts))
            allocate(bz%kpt_frac(3,bz%nkpts))
            !
            bz%kpt_frac = matmul(bz%bas,kpt_cart)
            bz%kpt_cart = kpt_cart
        else
            stop 'ERROR [create_bz]: either kpt_cart or kpt_frac must be present'
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
        ! make sure kpt_frac is betwen [0,1)
        bz%kpt_frac = modulo(bz%kpt_frac+prec,1.0_dp)-prec
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
        real(dp),allocatable :: kpt_frac(:,:)
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
        allocate(kpt_frac(3,nkpts))
        ! convert integer mesh to fractional coordinates
        do i = 1, nkpts
            kpt_frac(:,i) = ( kpt_int(:,i) + s ) / real(n,dp)
        enddo 
        ! create instance
        call fbz%create_bz(kpt_frac=kpt_frac, bas=pc%bas, prec=opts%prec)
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
            call disp(title='[frac.]',style='underline',X=transpose(fbz%kpt_frac),fmt='f10.5',advance='no')
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
        fbz%ntetra = size(fbz%tet,2)
        ! allocate tetra properties 
        allocate(fbz%tetv(fbz%ntetra))
        allocate(fbz%tetw(fbz%ntetra))
        ! primitive reciprocal cell volume
        recvol = det(fbz%recbas)
        ! get tetra volume = (volume/cell) / (nboxes/cell * ntetra/boxes) = (volume/cell) / (tetra/cell) = volume/tetra
        fbz%tetv = recvol/(fbz%nkpts*6)
        ! get tetra weight
        fbz%tetw = fbz%tetv/recvol
        ! print stdout
        if (opts%verbosity.ge.1) then
            ! tetrahedra properties
            write(*,*) flare, 'tetrahedra = '//tostring(fbz%ntetra)
            write(*,*) flare, 'tetrahedra volume = '//tostring(fbz%tetv(1))
            write(*,*) flare, 'tetrahedra weight = '//tostring(fbz%tetw(1))
            call disp_indent()
            call disp(title='#'                ,style='underline',X=[1:fbz%ntetra]    ,fmt='i5' ,advance='no')
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

!     ! I would prefer this method, but it takes forever because it needs to do a double loop over all kpoints
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
        integer :: i,j,k
        !
        if (opts%verbosity.ge.1) call print_title('Irreducible Wedge (IBZ)')
        !
        ! get permutation map which shows how k-points are permuted by each space symmetry operation, PM(fbz%nkpts,pg%nsyms) 
        PM = permutation_map( permutation_rep(seitz=pg%seitz_frac, tau=fbz%kpt_frac, flags='', prec=opts%prec) )
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
        allocate(ibz%kpt_frac, source=fbz%kpt_frac(:,ind(1:k)))
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
        ! allocate space for weights
        allocate(ibz%w(ibz%nkpts))
        ! get kpoint weights (number of unique points onto which symmetries can map ibz points)
        do i = 1, ibz%nkpts
            ibz%w(i) = count(unique_inds(PM(ibz%fbz_id(i),:)).ne.0)/real(fbz%nkpts,dp)
        enddo
        ! check weights
        if (abs(sum(ibz%w)-1.0_dp).gt.tiny) stop 'ERROR [get_ibz]: kpoint weights does not sum to 1'
        ! print stdout
        if (opts%verbosity.ge.1) then
            !
            write(*,'(a,a)') flare, 'FBZ k-points (full monkhorst-pack mesh) = '//tostring(fbz%nkpts)
            write(*,'(a,a)') flare, 'IBZ k-points = '//tostring(ibz%nkpts)
            call disp_indent()
            call disp(title='#'      ,style='underline',X=[1:ibz%nkpts]          ,fmt='i5'   ,advance='no')
            call disp(title='[cart.]',style='underline',X=transpose(ibz%kpt_cart),fmt='f10.5',advance='no')
            call disp(title='[frac.]',style='underline',X=transpose(ibz%kpt_frac),fmt='f10.5',advance='no')
            call disp(title='w'      ,style='underline',X=          ibz%w        ,fmt='f10.5',advance='yes')
            !
            write(*,'(a,a)',advance='no') flare, 'k-point map (to IBZ: full->irr)'
            call id_print_map(fbz%ibz_id)
            !
            write(*,'(a,a)',advance='no') flare, 'k-point map (from IBZ: irr->full)'
            call id_print_map(ibz%fbz_id)
            !
        endif
        !
    end subroutine get_irreducible

    ! Using this version instead
!     subroutine     get_irreducible(ibz,fbz,pc,pg,opts)
!         !
!         use am_options
!         !
!         implicit none
!         !
!         class(am_class_ibz), intent(inout) :: ibz
!         type(am_class_fbz) , intent(inout) :: fbz
!         type(am_class_point_group) , intent(in) :: pg
!         type(am_class_prim_cell), intent(in) :: pc
!         type(am_class_options)  , intent(in) :: opts
!         real(dp), allocatable :: korbit(:,:)
!         real(dp), allocatable :: kpt_ibz(:,:)
!         integer , allocatable :: inds(:)
!         integer , allocatable :: w(:)
!         integer :: i, j, k
!         !
!         !
!         if (opts%verbosity.ge.1) call print_title('Irreducible wedge (IBZ)')
!         ! 
!         if (opts%verbosity.ge.1) write(*,'(a,a)') flare, 'number of input kpoints = '//tostring(fbz%nkpts)
!         ! allocate space for orbits (fractional)
!         allocate(korbit(3,pg%nsyms))
!         ! allocate space for inds used to sort
!         allocate(inds(pg%nsyms))
!         ! allocate space for irreducible weights
!         allocate(w(fbz%nkpts))
!         ! allocate space for irreducible kpoints
!         allocate(kpt_ibz(3,fbz%nkpts))
!         !$OMP PARALLEL PRIVATE(i,j,korbit,inds) SHARED(fbz,pg,kpt_ibz)
!         !$OMP DO
!         do i = 1, fbz%nkpts
!             ! get orbit of kpoint already in fbz (wigner-seitz cell)
!             korbit = get_orbit(R=pg%seitz_cart(1:3,1:3,:),kpoint=fbz%kpt_cart(:,i))
!             ! sort orbit
!             do j = 1,3
!                 ! rank, removing numerical noise
!                 call rank(nint(korbit(j,:)*1.0D6),inds)
!                 ! rearrange
!                 korbit = korbit(:,inds)
!             enddo
!             ! get representative point
!             kpt_ibz(:,i) = korbit(:,1)
!             ! get irreducible weights vector
!             w(i) = count(unique_inds(korbit).ne.0)
!             ! check weights
!             if (modulo(pg%nsyms,w(i)).ne.0) stop 'ERROR [get_irreducible]: weight is not a factor of the number of symmetry operation'
!         enddo
!         !$OMP END DO
!         !$OMP END PARALLEL
!         ! repurpose inds
!         deallocate(inds)
!         ! get inds of unique kpoints in ibz
!         inds = unique_inds(kpt_ibz)
!         inds = trim_null(inds)
!         ! save ibz kpoint coordinates
!         allocate(ibz%kpt_cart,source=kpt_ibz(:,inds))
!         ! get weights
!         allocate(ibz%w,source= w(inds)/real(fbz%nkpts,dp))
!         ! check weights
!         if (abs(sum(ibz%w)-1.0_dp).gt.tiny) then
!             write(*,*) sum(ibz%w)
!             stop 'ERROR [reduce_to_ibz]: kpoint weights does not sum to 1'
!         endif
!         ! The error above could mean that:
!                 ! 1) there was an internal error. this has been occuring when an even MP mesh is requested
!                 ! 2) the full fbz mesh is not a MP mesh
!                 ! 3) the MP mesh contains different dimensions along x,y,z, i.e. 3 3 5  
!         ! get number of points in ibz
!         ibz%nkpts = count(inds.ne.0)
!         ! allocate space for maps
!         allocate(fbz%ibz_id(fbz%nkpts))
!         fbz%ibz_id = 0
!         allocate(ibz%fbz_id(ibz%nkpts))
!         ibz%fbz_id = 0
!         ! get map from full to irreducible brillouin zone and vice versa
!         do j = 1, ibz%nkpts
!         do i = 1, fbz%nkpts
!             if (isequal(korbit(:,i),ibz%kpt_cart(:,j))) then
!                 fbz%ibz_id(i) = j
!                 if (ibz%fbz_id(j).eq.0) ibz%fbz_id(i) = i
!             endif
!         enddo
!         enddo
!         ! get fractional point
!         allocate(ibz%kpt_frac,source=matmul(pc%bas,ibz%kpt_cart)) 
!         ! make sure cart is between [0,1)
!         do i = 1, ibz%nkpts
!             ibz%kpt_frac(:,i) = modulo(ibz%kpt_frac(:,i)+opts%prec,1.0_dp) - opts%prec
!         enddo
!         ! stdout
!         if (opts%verbosity.ge.1) then
!             !
!             write(*,'(a,a)') flare, 'FBZ k-points (full monkhorst-pack mesh) = '//tostring(fbz%nkpts)
!             write(*,'(a,a)') flare, 'IBZ k-points = '//tostring(ibz%nkpts)
!             call disp_indent()
!             call disp(title='#'      ,style='underline',X=[1:ibz%nkpts]          ,fmt='i5'   ,advance='no')
!             call disp(title='[cart.]',style='underline',X=transpose(ibz%kpt_cart),fmt='f10.5',advance='no')
!             call disp(title='[frac.]',style='underline',X=transpose(ibz%kpt_frac),fmt='f10.5',advance='no')
!             call disp(title='w'      ,style='underline',X=          ibz%w        ,fmt='f10.5',advance='yes')
!             !
!             write(*,'(a,a)',advance='no') flare, 'k-point map (to IBZ: full->irr)'
!             call id_print_map(fbz%ibz_id)
!             !
!             write(*,'(a,a)',advance='no') flare, 'k-point map (from IBZ: irr->full)'
!             call id_print_map(ibz%fbz_id)
!             !
!         endif
!         !
!     end subroutine get_irreducible

    ! functions which act on kpoints

!     function       locate_kpoint_in_fbz(kpoint,grid_points,bas) result(position_in_fbz)
!         !> Given a kpoint in fractional coordinates, determine if ...
!         !> k is inside  ( 1) of FBZ, i.e. 2KG < G^2 for ALL Bragg planes
!         !> k is on edge ( 0) of FBZ, i.e. 2KG = G^2 for ONE Bragg plane
!         !>                        and NOT 2KG > G^2 for ALL others
!         !> k is outside (-1) of FBZ, i.e. 2KG > G^2 for ONE Bragg plane
!         implicit none
!         !
!         real(dp), intent(in) :: kpoint(3)
!         real(dp), intent(in) :: bas(3,3)
!         real(dp), intent(in) :: grid_points(3,27)
!         real(dp) :: k_cart(3)
!         real(dp) :: P(27)
!         integer  :: position_in_fbz
!         integer  :: i
!         ! convert kpoint to cartesian
!         k_cart = matmul(inv(bas),kpoint)
!         !
!         P = 0
!         do i=1,size(grid_points,2)
!             P(i) = 2*dot_product(k_cart,grid_points(1:3,i)) - dot_product(grid_points(1:3,i),grid_points(1:3,i))
!         enddo
!         !
!         ! NOTE THAT THE ORDER IS VERY IMPORTANT.
!         ! If the point is outside ANY Bragg plane than it is outside the BZ
!         ! The point can only be on the edge if it NOT outside. Apply the
!         ! outside test last converts any points on Bragg plane, but outside the
!         ! FBZ, to just pure simple outside.
!         !
!         ! assume inside FBZ
!         position_in_fbz = 1;
!         ! test if on edge of FBZ
!         do i = 1,27
!             if (abs(P(i)).lt.tiny) then
!                 position_in_fbz = 0
!             endif
!         enddo
!         ! test if outside FBZ
!         do i = 1,27
!             if (abs(P(i)).gt.tiny) then
!                 position_in_fbz = -1
!             endif
!         enddo
!     end function   locate_kpoint_in_fbz

!     function       reduce_kpoint_to_fbz(kpoint,grid_points,bas) result(kpoint_fbz)
!         !> reduces kpoint (in fractional) to the first Brillouin zone (Wigner-Seitz cell, defined in cartesian coordinates)
!         implicit none
!         !
!         real(dp), intent(in) :: kpoint(3) !> fractional
!         real(dp), intent(in) :: bas(3,3) !> real space basis (column vectors)
!         real(dp), intent(in) :: grid_points(3,27) !> voronoi points (cartesian)
!         real(dp) :: k_cart(3) !> kpoint cartesian
!         real(dp) :: G(3) !> reciprocal lattice vector
!         real(dp) :: P !> bragg plane condition
!         real(dp) :: kpoint_fbz(3)
!         integer :: i ! loop variable
!         logical :: is_not_done
!         !
!         ! take kpoint in fractional coordinates (will become cartesian later)
!         k_cart = kpoint
!         ! reduce to reciprocal unit cell
!         k_cart = modulo(k_cart+tiny,1.0_dp)-tiny
!         ! convert to cartesian
!         k_cart = matmul(inv(bas),k_cart)
!         ! reduce to Wigner-Seitz cell (aka first Brillouin zone) by translating the k-point
!         ! until the closest reciprocal lattice point is [0 0 0]
!         is_not_done = .true.
!         do while ( is_not_done )
!             is_not_done = .false.
!             do i = 1,27
!                 G = grid_points(:,i)
!                 P = 2*dot_product(k_cart,G) - dot_product(G,G)
!                 if ( P .gt. tiny ) then
!                     k_cart = k_cart - G
!                     is_not_done = .true.
!                 endif
!             enddo
!         end do
!         ! convert back to fractional
!         kpoint_fbz = matmul(bas,k_cart)
!         !
!     end function   reduce_kpoint_to_fbz

!     function       reduce_kpoint_to_ibz(kpoint,grid_points,bas,R,prec) result(kpoint_in_ibz)
!         !
!         use am_rank_and_sort
!         !
!         implicit none
!         !
!         real(dp), intent(in) :: kpoint(3) !> fractional
!         real(dp), intent(in) :: bas(3,3)  !> real space basis (column vectors)
!         real(dp), intent(in) :: grid_points(3,27) !> voronoi points (cartesian)
!         real(dp), intent(in) :: R(:,:,:)  !> point symmetries in cartesian
!         real(dp), intent(in) :: prec
!         real(dp) :: kpoint_in_fbz(3)
!         real(dp) :: kpoint_in_ibz(3)
!         real(dp), allocatable :: korbit(:,:) !> orbit
!         integer , allocatable :: indices(:)
!         integer :: i, nsyms
!         !
!         nsyms = size(R,3)
!         ! reduce kpoint to primitive reciprocal lattice
!         kpoint_in_fbz = modulo(kpoint+prec,1.0_dp) - prec
!         ! get kpoint orbit
!         korbit = kpoint_orbit(R,kpoint_in_fbz)
!         ! reduce korbit to primitive reciprocal lattice
!         korbit = modulo(korbit+prec,1.0_dp) - prec
!         ! get unique points in orbit
!         korbit = unique(korbit,prec)
!         ! sort
!         allocate(indices(size(korbit,2)))
!         do i = 1,3
!             ! rank, removing numerical noise
!             call rank(nint(korbit(i,:)*1.0D6),indices)
!             korbit=korbit(:,indices)
!         enddo
!         ! get representative point
!         kpoint_in_ibz = korbit(:,1)
!         ! reduce representative irreducible point to fbz
!         kpoint_in_ibz = reduce_kpoint_to_fbz(kpoint=kpoint_in_ibz,grid_points=grid_points,bas=bas)
!         !
!     end function   reduce_kpoint_to_ibz

    pure function  get_orbit(R,kpoint) result(korbit)
        !> Given a kpoint in fractional coordinates, return its orbit (star)
        !> Note: this routine does not check to see whether the orbits are unique. If two point
        !> symmetries produce the same orbital point, this routine returns the same point twice.
        implicit none
        !
        real(dp), intent(in) :: R(:,:,:) !> point symmetries
        real(dp), intent(in) :: kpoint(3) !> cartesian
        real(dp) :: korbit(3,size(R,3)) !> orbit
        integer  :: npntsym ! number of point symmetries
        integer  :: i ! loop variable
        !
        ! get number of point symmetries
        npntsym = size(R,3)
        !
        do i = 1,npntsym
            korbit(:,i) = matmul(R(:,:,i),kpoint)
        enddo
        !
    end function   get_orbit

!     function       kpoint_weight(korbit,prec) result(w)
!         !
!         implicit none
!         !
!         real(dp), intent(in) :: korbit(:,:)
!         real(dp), intent(in) :: prec
!         integer  :: w, nsyms
!         !
!         w = count(unique_inds(korbit,prec).ne.0)
!         ! get number of symmetries
!         nsyms = size(R,3)
!         ! check that weight divides nsyms
!         if ( modulo(nsyms,w).ne.0) stop 'ERROR [kpoint_weight]: weight is not a factor of the number of symmetry operation'
!         !
!     end function   kpoint_weight

    function       get_tetrahedra(recbas,kpt_int,n) result(tet)
        !
        implicit none
        !
        real(dp),intent(in) :: recbas(3,3)
        integer, intent(in) :: kpt_int(:,:)
        integer, intent(in) :: n(3)
        integer, allocatable :: tet(:,:)
        integer, allocatable :: box(:,:) ! box(8,nboxes) vertices of boxes index kpoint in fbz 
        integer :: tetrahedron(4,6)      ! tetrahedron(nvertices,ntetra) four verticies of the six tetrahedra making up one box
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



