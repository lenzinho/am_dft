module am_brillouin_zone

    use am_unit_cell
    use am_symmetry
    use am_constants
    use am_stdout
    use am_vasp_io
    use am_mkl
    use am_options
    use dispmodule

    implicit none

    private

    ! BZ (primitive reciprocal-lattice cell defined with x,y,z (fractional) between [0,1).)

    type, public :: am_class_bz
        real(dp) :: bas(3,3)
        real(dp) :: recbas(3,3)
        integer :: nkpts                        ! nkpts number of kpoints
        real(dp), allocatable :: kpt_cart(:,: ) ! kpt(3,nkpts) kpoint
        real(dp), allocatable :: kpt_frac(:,:)  ! kpt(3,nkpts) kpoint - notice that vasp IBZKPT is in (reciprocal) fractional units.
        real(dp), allocatable :: w(:)           ! w(nkpts) normalized weights (CURRENTLY NOT USED!!!!)
        integer , allocatable :: fbz_id(:)
        integer , allocatable :: ibz_id(:)
        contains
        procedure :: create_bz
        procedure :: load
        procedure :: write_kpoints
        procedure :: debug_dump => debug_dump_bz
        ! add sort procedure
    end type am_class_bz

    ! FBZ (wigner-seitz cell defined in cartesian coordinates, but kpoints are saved as fractional coordinates)

    type, public, extends(am_class_bz) :: am_class_fbz
        contains
!         procedure :: get_fbz
    end type am_class_fbz

    ! IBZ (irreducible wedge)

    type, public, extends(am_class_bz) :: am_class_ibz
        contains
!         procedure :: get_ibz
    end type am_class_ibz

    ! PATH 


contains

    ! create

    subroutine     debug_dump_bz(bz,fname)
        !
        implicit none
        !
        class(am_class_bz), intent(in) :: bz
        character(*), intent(in) :: fname
        !
                                                call dump(A=bz%nkpts             ,fname=trim(fname)//'.nkpts'            ) 
                                                call dump(A=bz%bas               ,fname=trim(fname)//'.bas'              )
                                                call dump(A=bz%recbas            ,fname=trim(fname)//'.recbas'           )
        if (allocated(bz%kpt_cart            )) call dump(A=bz%kpt_cart          ,fname=trim(fname)//'.kpt_cart'         )
        if (allocated(bz%kpt_frac            )) call dump(A=bz%kpt_frac          ,fname=trim(fname)//'.kpt_frac'         )
        if (allocated(bz%w                   )) call dump(A=bz%w                 ,fname=trim(fname)//'.w'                )
        if (allocated(bz%fbz_id              )) call dump(A=bz%fbz_id            ,fname=trim(fname)//'.fbz_id'           )
        if (allocated(bz%ibz_id              )) call dump(A=bz%ibz_id            ,fname=trim(fname)//'.ibz_id'           )
        !
    end subroutine debug_dump_bz

    subroutine     create_bz(bz,kpt_frac,kpt_cart,bas,prec)
        !
        ! kpt_frac points are shifted to be between [0,1)
        ! kpt_cart        are shifted to be in the wigner-seitz bz
        !
        implicit none
        !
        class(am_class_bz), intent(out) :: bz
        real(dp), intent(in), optional :: kpt_frac(:,:)
        real(dp), intent(in), optional :: kpt_cart(:,:)
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
        ! generate voronoi points [cart]
        grid_points = meshgrid([-1:1],[-1:1],[-1:1])
        grid_points = matmul(bz%recbas,grid_points)
        ! make sure kpt_cart is in wigner-seitz cell 
        do i = 1, bz%nkpts
        bz%kpt_cart(:,i) = reduce_kpoint_to_fbz(kpoint_cart=bz%kpt_cart(:,i), grid_points=grid_points)
        enddo
        ! make sure kpt_frac is betwen [0,1)
        bz%kpt_frac = modulo(bz%kpt_frac+prec,1.0_dp)-prec
        ! debug dump
        if (debug) then
            call execute_command_line ('mkdir -p '//trim(debug_dir)//'/bz')
            call bz%debug_dump(fname=               trim(debug_dir)//'/bz/outfile.bz')
        endif
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

    subroutine     load(bz,uc,pg,opts,flags)
        !
        implicit none
        !
        class(am_class_bz)        , intent(out):: bz
        type(am_class_unit_cell)  , intent(in) :: uc
        type(am_class_point_group), intent(in) :: pg
        type(am_class_options)    , intent(in) :: opts
        character(*)              , intent(in) :: flags
        real(dp), allocatable :: kpt(:,:) ! kpoint coordinates
        real(dp), allocatable :: w(:)     ! normalized weights
        !
        if (opts%verbosity.ge.1) call print_title('Input kpoints')
        !
        if     (index(flags,'ibzkpt')) then
            write(*,'(a,a,a)') flare, 'file = ', 'ibzkpt'
            call read_ibzkpt(  kpt=kpt, w=w, iopt_filename=opts%ibzkpt  , iopt_verbosity=0)
        elseif (index(flags,'procar')) then
            write(*,'(a,a,a)') flare, 'file = ', 'procar'
            call read_procar(  kpt=kpt, w=w, iopt_filename=opts%procar  , iopt_verbosity=0)
        elseif (index(flags,'eigenval')) then
            write(*,'(a,a,a)') flare, 'file = ', 'eigenval'
            call read_eigenval(kpt=kpt, w=w, iopt_filename=opts%eigenval, iopt_verbosity=0)
        else
            stop 'Unknown flag. Nothing read.'
        endif
        !
        call bz%create_bz(kpt_frac=kpt, bas=uc%bas, prec=opts%prec)
        !
        ! call bz%get_weights(pg=pg, opts=opts)
        ! if (.not.isequal(bz%w,w)) stop 'Computed weights do not match vasp input'
        !
        if (opts%verbosity.ge.1) then
            write(*,'(a,a,a)') flare, 'kpoints = ', tostring(bz%nkpts)
        endif
        !
    end subroutine load

    subroutine     write_kpoints(bz,fname)
        !
        implicit none
        !
        class(am_class_bz), intent(in) :: bz
        character(*)      , intent(in) :: fname
        integer :: fid
        !
        fid = 1
        open(unit=fid,file=trim(fname),status='replace',action='write')
            call disp(unit=fid,title='#'      ,style='underline',X=[1:bz%nkpts]          ,fmt='i5'   ,advance='no')
            call disp(unit=fid,title='[cart.]',style='underline',X=transpose(bz%kpt_cart),fmt='f10.5',advance='no')
            call disp(unit=fid,title='[frac.]',style='underline',X=transpose(bz%kpt_frac),fmt='f10.5',advance='yes')
            ! call disp(unit=fid,title='w'      ,style='underline',X=bz%w                  ,fmt='f10.5',advance='yes')
        close(fid)
        !
    end subroutine write_kpoints


!     subroutine     get_weights(bz,pg,opts)
!         !
!         implicit none
!         !
!         class(am_class_bz)        , intent(inout) :: bz
!         type(am_class_point_group), intent(in)    :: pg
!         type(am_options)          , intent(in)    :: opts
!         integer :: i, j
!         integer :: matches
!         !
!         allocate(bz%w(bz%nkpts))
!         !
!         do i = bz%nkpts
!             !
!             orbit = unique(kpoint_orbit(R=pg%seitz(1:3,1:3,:), kpoint=bz%kpt(:,i)))
!             !
!             matches = 0 
!             do j = 1, size(orbit,2)
!             do k = 1, bz%nkpts

!             enddo
!             enddo
!             !
!             bz%w(i) = kpoint_weight(R=pg%seitz(1:3,1:3,:), kpoint=, prec=opts%prec)
!         enddo
!         !
!         contains
!         function       kpoint_weight(R,kpoint,prec) result(w)
!             !
!             implicit none
!             !
!             real(dp), intent(in) :: R(:,:,:) !> point symmetries
!             real(dp), intent(in) :: kpoint(3) !> cartesian
!             real(dp), intent(in) :: prec
!             integer  :: w, nsyms
!             !
!             nsyms = size(R,3)
!             !
!             w = size(unique(kpoint_orbit(R,kpoint),prec),2)
!             !
!             if ( modulo(nsyms,w).ne.0) then
!                 call am_print('ERROR','Weight is not a factor of the number of symmetry operation',flags='E')
!                 call am_print('kpoint (fractional)',kpoint)
!                 call am_print('weight (integer)',w)
!                 call am_print('number of symmetry operations',nsyms)
!                 call am_print('kpoint orbit',transpose(unique(kpoint_orbit(R,kpoint),prec)))
!                 stop
!             endif
!             !
!         end function   kpoint_weight
!         pure function  kpoint_orbit(R,kpoint) result(korbit)
!             !> Given a kpoint in fractional coordinates, return its orbit (star)
!             !> Note: this routine does not check to see whether the orbits are unique. If two point
!             !> symmetries produce the same orbital point, this routine returns the same point twice.
!             implicit none
!             !
!             real(dp), intent(in) :: R(:,:,:) !> point symmetries
!             real(dp), intent(in) :: kpoint(3) !> cartesian
!             real(dp) :: korbit(3,size(R,3)) !> orbit
!             integer  :: npntsym ! number of point symmetries
!             integer  :: i ! loop variable
!             !
!             ! get number of point symmetries
!             npntsym = size(R,3)
!             !
!             do i = 1,npntsym
!                 korbit(:,i)=matmul(R(:,:,i),kpoint)
!             enddo
!             !
!         end function   kpoint_orbit
!     end subroutine get_weights

    ! goal: start with fbz (monkhorst-pack mesh) => generate ibz, done.
    !       alternatively, load vasp directly into fbz and reduce to ibz, or directly into ibz.

!     subroutine     get_fbz(fbz,uc,n,s,opts)
!         !
!         implicit none
!         !
!         class(am_class_fbz)     , intent(out):: fbz
!         type(am_class_unit_cell), intent(in) :: uc
!         integer                 , intent(in) :: n(3)
!         real(dp)                , intent(in) :: s(3)
!         type(am_class_options)  , intent(in) :: opts
!         !
!         if (opts%verbosity.ge.1) call print_title('Generating Monkhorst-Pack mesh on FBZ')
!         ! get k-points
!         call fbz%create_bz(kpt_frac=generate_monkhorst_pack_mesh(n=n,s=s), recbas=uc%recbas, bas=uc%bas, prec=opts%prec)
!         ! get weights
!         call fbz%get_weights()
!         ! print stdout
!         if (opts%verbosity.ge.1) then
!             ! am_print_two_matrices_side_by_side(name, Atitle, Btitle, A, B , in_emph, iopt_fid )
!             call am_print_two_matrices_side_by_side(name='kpoints',&
!                 Atitle='fractional',A=transpose(fbz%kpt_frac),&
!                 Btitle='cartesian' ,B=transpose(fbz%kpt_cart),&
!                 iopt_emph=flare,iopt_teaser=.true.)
!         endif
!         !
!         contains
!         pure function  generate_monkhorst_pack_mesh(n,s) result(kpt)
!             !> returns kpoints in fractional coordinates for monkhorst pack mesh dimensions
!             !> n(1:3)=[n1,n2,n3] and shift s(1:3)=[s1,s2,s3]
!             !> according to http://cms.mpi.univie.ac.at/vasp/vasp/Automatic_k_mesh_generation.html
!             implicit none
!             !
!             integer, intent(in) :: n(3) !> monkhorst pack mesh dimensions
!             real(dp), intent(in) :: s(3) !> monkhorst pack mesh shift
!             real(dp), allocatable :: kpt(:,:) !> kpt in fractional
!             integer :: i1,i2,i3,j
!             !
!             allocate(kpt(3,product(n)))
!             !
!             j=0
!             do i1=1,n(1)
!                 do i2=1,n(2)
!                     do i3=1,n(3)
!                         j=j+1
!                         kpt(1:3,j) = real([i1+s(1)-1,i2+s(2)-1,i3+s(3)-1],dp)
!                         kpt(1:3,j) = kpt(1:3,j)/real(n,dp)
!                     enddo
!                 enddo
!             enddo
!             !
!         end function   generate_monkhorst_pack_mesh
!     end subroutine get_fbz

!     subroutine     get_ibz(ibz,fbz,pg,opts)
!         !
!         implicit none
!         !
!         class(am_class_ibz)        , intent(out)   :: ibz
!         class(am_class_bz)         , intent(inout) :: fbz
!         class(am_class_point_group), intent(in)    :: pg
!         type(am_class_options)     , intent(in)    :: opts
!         integer, allocatable :: PM(:,:)
!         logical, allocatable :: mask(:)
!         integer, allocatable :: ind(:)
!         integer :: i,j,k
!         !
!         if (opts%verbosity.ge.1) call print_title('Reducing to irreducible brillouin zone')
!         !
!         ! get permutation map which shows how k-points are permuted by each space symmetry operation, PM(fbz%nkpts,pg%nsyms) 
!         PM = permutation_map( permutation_rep(seitz=pg%seitz_cart, tau=fbz%kpt_cart, flags='', prec=opts%prec) )
!         ! determine irreducible kpoints
!         allocate(mask(fbz%nkpts)); mask = .true.
!         allocate(ind(fbz%nkpts)) ; ind  = 0
!         k=0
!         do i = 1, fbz%nkpts
!         if (mask(i)) then
!             k=k+1
!             ind(k)=i
!             do j = 1, pg%nsyms
!                 mask(PM(i,j))=.false.
!             enddo
!         endif
!         enddo
!         ! get number of k-points
!         ibz%nkpts = k
!         ! transfer irreducible k-points (frac)
!         allocate(ibz%kpt_frac, source=fbz%kpt_frac(:,ind(1:k)))
!         ! transfer irreducible k-points (cart)
!         allocate(ibz%kpt_cart, source=fbz%kpt_cart(:,ind(1:k)))
!         ! map irreducible k-point -> irreducible k-point
!         allocate(ibz%ibz_id, source=[1:ibz%nkpts])
!         ! map irreducible k-point onto -> primitive k-point
!         allocate(ibz%fbz_id, source=ind(1:k))
!         ! map primitive k-point onto -> irreducible k-point
!         allocate(fbz%ibz_id(fbz%nkpts))
!         fbz%ibz_id = 0
!         do i = 1, ibz%nkpts
!             ! PM(1,:) shows all k-points onto which k-point 1 is mapped by all space symmetry operations
!             do j = 1, fbz%nkpts
!                 search : do k = 1, pg%nsyms
!                     if (PM(j,k).eq.ibz%fbz_id(i)) then
!                         fbz%ibz_id(j) = i
!                         exit search
!                     endif
!                 enddo search
!             enddo
!         enddo
!         if (any(fbz%ibz_id.eq.0)) stop 'ERROR: full->ibz mapping failed.'
!         !
!         ! print stdout
!         if (opts%verbosity.ge.1) then
!             !
!             call am_print('FBZ k-points',fbz%nkpts,flare)
!             !
!             call am_print('IBZ k-points',ibz%nkpts,flare)
!             !
!             call am_print_two_matrices_side_by_side(name='IBZ k-points',&
!                 Atitle='fractional',A=transpose(ibz%kpt_frac),&
!                 Btitle='cartesian' ,B=transpose(ibz%kpt_cart),&
!             iopt_emph=flare,iopt_teaser=.true.)
!             !
!             write(*,'(a,a)',advance='no') flare, 'atomic mapping (to IBZ: full->irr)'
!             call id_print_map(fbz%ibz_id)
!             !
!             write(*,'(a,a)',advance='no') flare, 'atomic mapping (from IBZ: irr->full)'
!             call id_print_map(ibz%fbz_id)
!             !
!         endif
!         !
!     end subroutine get_ibz

!     subroutine     get_ibz_v1(ibz,fbz,pc,pg,opts)
!         !
!         implicit none
!         !
!         class(am_class_ibz)    , intent(inout) :: ibz
!         class(am_class_bz)     , intent(inout) :: fbz
!         type(am_class_point_group), intent(in) :: pg
!         type(am_class_prim_cell)  , intent(in) :: pc
!         type(am_class_options)    , intent(in) :: opts
!         real(dp), allocatable :: grid_points(:,:)
!         integer :: i, j
!         !
!         !
!         if (opts%verbosity.ge.1) call print_title('Reducing to IBZ')
!         ! 
!         if (opts%verbosity.ge.1) call am_print('number of original kpoints',bz%nkpts,flare)
!         ! generate voronoi points (cartesian)
!         grid_points = matmul( inv(pc%bas), meshgrid([-1:1],[-1:1],[-1:1]) )
!         ! reduce kpoints to ibz
!         allocate(ibz%kpt(3,bz%nkpts))
!         !$OMP PARALLEL PRIVATE(i) SHARED(grid_points,bz,ibz,pc,pg)
!         !$OMP DO
!         do i = 1, bz%nkpts
!             ibz%kpt(:,i) = reduce_kpoint_to_ibz(kpoint=bz%kpt(:,i),grid_points=grid_points,bas=pc%bas,R=pg%sym(1:3,1:3,:),prec=opts%prec)
!         enddo
!         !$OMP END DO
!         !$OMP END PARALLEL
!         ! get unique points
!         ibz%kpt = unique(ibz%kpt,opts%prec)
!         ! get number of unique points
!         ibz%nkpts = size(ibz%kpt,2)
!         if (opts%verbosity.ge.1) call am_print('number of irreducible kpoints',ibz%nkpts,flare)
!         ! get weights
!         allocate(ibz%w(ibz%nkpts))
!         !$OMP PARALLEL PRIVATE(i) SHARED(grid_points,bz,ibz,pc,pg)
!         !$OMP DO
!         do i = 1, ibz%nkpts
!             ibz%w(i) = kpoint_weight(R=pg%sym(1:3,1:3,:),kpoint=ibz%kpt(:,i),prec=opts%prec)
!         enddo
!         !$OMP END DO
!         !$OMP END PARALLEL
!         if (opts%verbosity.ge.1) then
!             if (.not.nint(sum(ibz%w)).eq.bz%nkpts) then
!                 call am_print('WARNING','The sum of irreducible weights does not equal the number of k-points in the full zone.',flags='E')
!                 ! This could mean that:
!                 ! 1) there was an internal error. this has been occuring when an even MP mesh is requested
!                 ! 2) the full bz mesh is not a MP mesh
!                 ! 3) the MP mesh contains different dimensions along x,y,z, i.e. 3 3 5  
!             endif
!         endif
!         if (opts%verbosity.ge.1) call print_weights(ibz)
!         ! normalize weights
!         ibz%w = ibz%w/real(sum(ibz%w),dp)
!         ! set irreducible indices
!         if (.not.allocated(bz%ibz_id)) allocate(bz%ibz_id(bz%nkpts))
!         irrk : do i = 1,ibz%nkpts
!                do j = 1, bz%nkpts
!                    if ( all(abs(bz%kpt(:,j) - ibz%kpt(:,i)).lt.tiny) ) then
!                        bz%ibz_id(j) = i
!                        cycle irrk
!                    endif
!                enddo
!         enddo irrk
!         ! stdout
!         if (opts%verbosity.ge.1) then
!             call am_print_two_matrices_side_by_side(name='irreducible kpoints',&
!                 Atitle='fractional',A=transpose(ibz%kpt),&
!                 Btitle='cartesian' ,B=transpose(matmul(inv(pc%bas),ibz%kpt)),&
!                 iopt_emph=flare,iopt_teaser=.true.)
!         endif
!         !
!         contains
!         function       reduce_kpoint_to_ibz(kpoint,R,prec) result(kpoint_in_ibz)
!             !
!             use am_rank_and_sort
!             !
!             implicit none
!             !
!             real(dp), intent(in) :: kpoint(3) !> fractional
!             real(dp), intent(in) :: R(:,:,:)  !> point symmetries in cartesian
!             real(dp), intent(in) :: prec
!             real(dp) :: kpoint_in_fbz(3)
!             real(dp) :: kpoint_in_ibz(3)
!             real(dp), allocatable :: korbit(:,:) !> orbit
!             integer , allocatable :: indices(:)
!             integer :: i, nsyms
!             !
!             nsyms = size(R,3)
!             !
!             ! reduce kpoint to primitive reciprocal lattice
!             kpoint_in_fbz = modulo(kpoint+prec,1.0_dp) - prec
!             ! get kpoint orbit
!             korbit = kpoint_orbit(R,kpoint_in_fbz)
!             ! reduce korbit to primitive reciprocal lattice
!             korbit = modulo(korbit+prec,1.0_dp) - prec
!             ! get unique points in orbit
!             korbit = unique(korbit,prec)
!             ! sort
!             allocate(indices(size(korbit,2)))
!             do i = 1,3
!                 ! rank to remove numerical noise
!                 call rank(nint(korbit(i,:)*1.0D6),indices)
!                 korbit=korbit(:,indices)
!             enddo
!             ! get representative point
!             kpoint_in_ibz = korbit(:,1)
!             !
!         end function   reduce_kpoint_to_ibz
!     end subroutine get_ibz_v1

!     subroutine     print_weights(bz)
!         !
!         implicit none
!         !
!         class(am_class_bz) :: bz
!         integer :: i,j
!         !
!         call am_print('sum weight',nint(sum(bz%w)),flare)
!         call am_print('min weight',nint(minval(bz%w)),flare)
!         call am_print('max weight',nint(maxval(bz%w)),flare)
!         call am_print('avg weight',sum(bz%w)/real(bz%nkpts,dp),flare)
!         write(*,'(a,a)') flare, 'histogram of weights'
!         write(*,'(5x,a10,a10)' ) 'weights', '# kpts'
!         write(*,'(5x,a10,a10)' ) repeat('-',8), repeat('-',8)
!         do i = minval(bz%w), maxval(bz%w)
!             j = count(bz%w.eq.i)
!             if (j.ne.0) write(*,'(5x,i10,i10)' ) i, j
!         enddo
!         !
!     end subroutine print_weights

!     ! functions which operate on kpoint

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
!         !
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



