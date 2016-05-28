module am_brillouin_zone

    use am_prim_cell
    use am_symmetry
    use am_constants
    use am_stdout
    use am_vasp_io
    use am_mkl
    use am_options
    use am_unit_cell

    implicit none

    private

    ! BZ (primitive reciprocal-lattice cell defined with x,y,z (fractional) between [0,1).)

    type, public :: am_class_bz
        integer :: nkpts                       ! nkpts number of kpoints
        real(dp), allocatable :: kpt_cart(:,:) ! kpt(3,nkpts) kpt - notice that vasp IBZKPT is in fractional units.
        real(dp), allocatable :: kpt_frac(:,:) ! kpt(3,nkpts) kpt - notice that vasp IBZKPT is in fractional units.
        real(dp), allocatable :: w(:)          ! w(nkpts) normalized weights
        integer , allocatable :: fbz_id(:)
        integer , allocatable :: ibz_id(:)
    contains
        procedure :: create_reciprocal_lattice
        procedure :: load_ibzkpt
        procedure :: load_procar
        procedure :: load_eigenval
        procedure :: write_kpoints
    end type am_class_bz

    ! PRL (primitive reciprocal lattice), defined between [0,1)

    type, public, extends(am_class_bz) :: am_class_prl
    contains
        procedure :: get_monkhorst_pack_mesh
    end type am_class_prl

    ! FBZ (reciprocal wigner-seitz cell)

    type, public, extends(am_class_bz) :: am_class_fbz
    contains
        procedure :: get_fbz
    end type am_class_fbz

    ! IBZ (irreducible wedge)

    type, public, extends(am_class_bz) :: am_class_ibz
    contains
        procedure :: get_ibz
    end type am_class_ibz

    ! PATH 



contains

    ! load from vasp

    subroutine     load_ibzkpt(bz,uc,opts)
        ! Reads the IBZKPT file.
        implicit none
        !
        class(am_class_bz), intent(out) :: bz
        type(am_class_unit_cell) , intent(in) :: uc
        type(am_class_options), intent(in) :: opts
        !
        call read_ibzkpt(nkpts = bz%nkpts,&
                         kpt   = bz%kpt_frac,&
                         w     = bz%w,&
                 iopt_filename = opts%ibzkpt,&
                iopt_verbosity = opts%verbosity )
        !
        allocate(bz%kpt_cart,source=matmul(uc%recbas,bz%kpt_frac))
        !
    end subroutine load_ibzkpt

    subroutine     load_procar(bz,uc,opts)
        ! Reads the PROCAR file.
        implicit none
        !
        class(am_class_bz), intent(out) :: bz
        type(am_class_unit_cell), intent(in) :: uc
        type(am_class_options), intent(in) :: opts
        !
        call read_procar(nkpts    = bz%nkpts,&
                         kpt      = bz%kpt_frac,&
                         w        = bz%w,&
                         iopt_filename = opts%procar,&
                         iopt_verbosity= 0 )
        !
        allocate(bz%kpt_cart,source=matmul(uc%recbas,bz%kpt_frac))
        !
    end subroutine load_procar

    subroutine     load_eigenval(bz,uc,opts)
        ! Reads the eigenval file.
        implicit none
        !
        class(am_class_bz)  , intent(out) :: bz
        type(am_class_unit_cell), intent(in) :: uc
        type(am_class_options), intent(in) :: opts
        !
        call read_eigenval(nkpts  = bz%nkpts,&
                           kpt    = bz%kpt_frac,&
                           w      = bz%w,&
                           iopt_filename = opts%eigenval,&
                           iopt_verbosity = 0 )
        !
        allocate(bz%kpt_cart,source=matmul(uc%recbas,bz%kpt_frac))
        !
    end subroutine load_eigenval
   
    ! functions which operate on bz or similar derived types
    
    subroutine     create_reciprocal_lattice(bz,recbas,kpt_frac)
        !
        implicit none
        !
        class(am_class_bz), intent(out):: bz
        real(dp), intent(in) :: recbas(3,3)
        real(dp), intent(in) :: kpt_frac(:,:)
        !
        bz%nkpts = size(kpt_frac,2)
        allocate(bz%kpt_frac, source=kpt_frac)
        allocate(bz%kpt_cart, source=matmul(recbas,bz%kpt_frac))
        !
    end subroutine create_reciprocal_lattice

    subroutine     get_monkhorst_pack_mesh(prl,pc,n,s,opts)
        !
        implicit none
        !
        class(am_class_prl)     , intent(out):: prl
        type(am_class_prim_cell), intent(in) :: pc
        integer                 , intent(in) :: n(3)
        real(dp)                , intent(in) :: s(3)
        type(am_class_options)  , intent(in) :: opts
        !
        if (opts%verbosity.ge.1) call am_print_title('Generating Monkhorst-Pack mesh')
        !
        ! reduce korbit to primitive reciprocal lattice
        call prl%create_reciprocal_lattice(recbas=pc%recbas, kpt_frac=generate_monkhorst_pack_mesh(n=n,s=s))
        !
        ! print stdout
        if (opts%verbosity.ge.1) then
            call am_print('kpoints (prl)',prl%nkpts,' ... ')
            call am_print_two_matrices_side_by_side(name='kpoints (prl)',&
                Atitle='fractional',A=transpose(prl%kpt_frac),&
                Btitle='cartesian' ,B=transpose(prl%kpt_cart),&
                iopt_emph=' ... ',iopt_teaser=.true.)
        endif
        !
        contains
        pure function  generate_monkhorst_pack_mesh(n,s) result(kpt_frac)
            !> returns kpoints in fractional coordinates for monkhorst pack mesh dimensions
            !> n(1:3)=[n1,n2,n3] and shift s(1:3)=[s1,s2,s3]
            !> according to http://cms.mpi.univie.ac.at/vasp/vasp/Automatic_k_mesh_generation.html
            implicit none
            !
            integer , intent(in) :: n(3) !> monkhorst pack mesh dimensions
            real(dp), intent(in) :: s(3) !> monkhorst pack mesh shift
            real(dp), allocatable :: kpt_frac(:,:) !> kpt in fractional
            integer :: i1,i2,i3,j
            !
            allocate(kpt_frac(3,product(n)))
            !
            j=0
            do i1=1,n(1)
                do i2=1,n(2)
                    do i3=1,n(3)
                        j=j+1
                        kpt_frac(1:3,j) = real([i1+s(1)-1,i2+s(2)-1,i3+s(3)-1],dp)
                        kpt_frac(1:3,j) = kpt_frac(1:3,j)/real(n,dp)
                    enddo
                enddo
            enddo
            !
        end function   generate_monkhorst_pack_mesh
    end subroutine get_monkhorst_pack_mesh

    subroutine     get_prl(prl,bz,pc,opts)
        !
        implicit none
        !
        class(am_class_prl)     , intent(out) :: prl
        class(am_class_bz)      , intent(in)  :: bz
        type(am_class_prim_cell), intent(in)  :: pc
        type(am_class_options)  , intent(in)  :: opts
        !
        if (opts%verbosity.ge.1) then
            call am_print_title('Determining primitive reciprocal lattice')
        endif
        !
        ! reduce korbit to primitive reciprocal lattice
        call prl%create_reciprocal_lattice(recbas=pc%recbas, &
            kpt_frac=(modulo(bz%kpt_frac+opts%prec,1.0_dp)-opts%prec) )
        !
        ! print stdout
        if (opts%verbosity.ge.1) then
            call am_print('kpoints (prl)',prl%nkpts,' ... ')
            call am_print_two_matrices_side_by_side(name='kpoints (prl)',&
                Atitle='fractional',A=transpose(prl%kpt_frac),&
                Btitle='cartesian' ,B=transpose(prl%kpt_cart),&
                iopt_emph=' ... ',iopt_teaser=.true.)
        endif
    end subroutine get_prl

    subroutine     get_fbz(fbz,bz,pc,pg,opts)
        !
        implicit none
        !
        class(am_class_fbz)       , intent(out):: fbz
        class(am_class_bz)        , intent(in) :: bz
        type(am_class_prim_cell)  , intent(in) :: pc
        type(am_class_point_group), intent(in) :: pg
        type(am_class_options)    , intent(in) :: opts
        real(dp), allocatable :: kpt_orbit_cart(:,:)
        real(dp), allocatable :: fbz_kpt_cart(:,:)
        real(dp), allocatable :: vpt_cart(:,:)
        integer :: i, j, m
        integer :: w
        !
        if (opts%verbosity.ge.1) then
            !
            call am_print_title('Expand kpoints to FBZ')
            call am_print('point symmetries',pg%nsyms,' ... ')
            call am_print('kpoints (input)' ,bz%nkpts,' ... ')
            call am_print_two_matrices_side_by_side(name='kpoints (input)',&
                Atitle='fractional',A=transpose(bz%kpt_frac),&
                Btitle='cartesian' ,B=transpose(bz%kpt_cart),&
                iopt_emph=' ... ',iopt_teaser=.true.)
            !
        endif
        !
        ! get voronoi points
        vpt_cart = meshgrid([-1:1],[-1:1],[-1:1])
        vpt_cart = matmul(pc%recbas, vpt_cart)
        !
        allocate(fbz_kpt_cart(3,bz%nkpts*pg%nsyms)) ! wkrspace
        m = 0
        ! expand each kpt onto the fbz
        !$OMP PARALLEL PRIVATE(i,j,w,kpt_orbit_cart) SHARED(fbz_kpt_cart,pg,m)
        !$OMP DO
        do i = 1, bz%nkpts
            ! get its orbit
            kpt_orbit_cart = get_kpoint_orbit(pg%seitz_cart(1:3,1:3,:), bz%kpt_cart(:,i))
            ! reduce korbit to primitive reciprocal lattice
            kpt_orbit_cart = modulo(kpt_orbit_cart+opts%prec,1.0_dp)-opts%prec
            ! reduce to fbz
            do j = 1, pg%nsyms
            kpt_orbit_cart(:,j) = reduce_kpoint_to_fbz(kpt_cart=kpt_orbit_cart(:,j), vpt_cart=vpt_cart)
            enddo
            ! get unique values
            kpt_orbit_cart = unique(kpt_orbit_cart,opts%prec)
            ! get its weight
            w = size(kpt_orbit_cart,2)
            ! check that weight is a factor of the number of symmetry operations
            if ( modulo(pg%nsyms,w)/=0) stop 'Weight is not a factor of the number of symmetry operation'
            ! save all unique points in orbit as part of the fbz
            do j = 1, w
                m=m+1
                fbz_kpt_cart(:,m) = kpt_orbit_cart(:,j)
            enddo
        enddo
        !$OMP END DO
        !$OMP END PARALLEL
        !
        ! get kpt [cart]
        fbz%kpt_cart = unique(fbz_kpt_cart(:,1:m), opts%prec)
        ! get kpt [frac], defined between [0,1)
        fbz%kpt_frac = matmul(pc%recbas, fbz%kpt_frac)
        ! get number of kpoints in the fbz
        fbz%nkpts = size(fbz%kpt_cart,2)
        ! get weights
        allocate(fbz%w(fbz%nkpts))
        fbz%w = 1.0_dp/real(fbz%nkpts,dp)
        !
        ! print stdout
        if (opts%verbosity.ge.1) then
            !
            call am_print('kpoints (fbz)',fbz%nkpts,' ... ')
            call am_print_two_matrices_side_by_side(name='kpoints (fbz)',&
                Atitle='fractional',A=transpose(fbz%kpt_frac),&
                Btitle='cartesian' ,B=transpose(fbz%kpt_cart),&
                iopt_emph=' ... ',iopt_teaser=.true.)
            !
        endif
    end subroutine get_fbz

    subroutine     get_ibz(ibz,bz,pc,pg,opts)
        !
        implicit none
        !
        class(am_class_ibz)    , intent(inout) :: ibz
        class(am_class_bz)     , intent(inout) :: bz
        type(am_class_point_group), intent(in) :: pg
        type(am_class_prim_cell)  , intent(in) :: pc
        type(am_class_options)    , intent(in) :: opts
        real(dp), allocatable :: vpt_cart(:,:)
        integer :: i, j
        !
        !
        if (opts%verbosity.ge.1) then
            call am_print_title('Reducing to IBZ')
            call am_print('kpoints (input)',bz%nkpts,' ... ')
        endif
        ! generate voronoi points [cart]
        vpt_cart = meshgrid([-1:1],[-1:1],[-1:1])
        vpt_cart = matmul(pc%recbas, vpt_cart)
        ! reduce kpoints to ibz
        allocate(ibz%kpt(3,bz%nkpts))
        !$OMP PARALLEL PRIVATE(i) SHARED(vpt_cart,bz,ibz,pc,pg)
        !$OMP DO
        do i = 1, bz%nkpts
            ibz%kpt(:,i) = reduce_kpoint_to_ibz(kpt=bz%kpt(:,i),vpt_cart=vpt_cart,bas=pc%bas,R=pg%sym(1:3,1:3,:),prec=opts%prec)
        enddo
        !$OMP END DO
        !$OMP END PARALLEL
        ! get unique points
        ibz%kpt = unique(ibz%kpt,opts%prec)
        ! get number of unique points
        ibz%nkpts = size(ibz%kpt,2)
        if (opts%verbosity.ge.1) call am_print('number of irreducible kpoints',ibz%nkpts,' ... ')
        ! get weights
        allocate(ibz%w(ibz%nkpts))
        !$OMP PARALLEL PRIVATE(i) SHARED(vpt_cart,bz,ibz,pc,pg)
        !$OMP DO
        do i = 1, ibz%nkpts
            ibz%w(i) = kpoint_weight(R=pg%sym(1:3,1:3,:),kpt=ibz%kpt(:,i),prec=opts%prec)
        enddo
        !$OMP END DO
        !$OMP END PARALLEL
        if (opts%verbosity.ge.1) then
            if (.not.nint(sum(ibz%w)).eq.bz%nkpts) then
                call am_print('WARNING','The sum of irreducible weights does not equal the number of k-points in the full zone.',flags='E')
                ! This could mean that:
                ! 1) there was an internal error. this has been occuring when an even MP mesh is requested
                ! 2) the full bz mesh is not a MP mesh
                ! 3) the MP mesh contains different dimensions along x,y,z, i.e. 3 3 5  
            endif
        endif
        if (opts%verbosity.ge.1) call print_weights(ibz)
        ! normalize weights
        ibz%w = ibz%w/real(sum(ibz%w),dp)
        ! set irreducible indices
        if (.not.allocated(bz%ibz_id)) allocate(bz%ibz_id(bz%nkpts))
        irrk : do i = 1,ibz%nkpts
               do j = 1, bz%nkpts
                   if ( all(abs(bz%kpt(:,j) - ibz%kpt(:,i)).lt.tiny) ) then
                       bz%ibz_id(j) = i
                       cycle irrk
                   endif
               enddo
        enddo irrk
        ! stdout
        if (opts%verbosity.ge.1) then
            call am_print_two_matrices_side_by_side(name='irreducible kpoints',&
                Atitle='fractional',A=transpose(ibz%kpt),&
                Btitle='cartesian' ,B=transpose(matmul(inv(pc%bas),ibz%kpt)),&
                iopt_emph=' ... ',iopt_teaser=.true.)
        endif
        !
        contains
        function       reduce_kpoint_to_ibz(kpt,vpt_cart,R,prec) result(ibz_kpt_cart)
            !
            use am_rank_and_sort
            !
            implicit none
            !
            real(dp), intent(in) :: kpt(3)         ! fractional
            real(dp), intent(in) :: vpt_cart(3,27) ! voronoi points (cartesian)
            real(dp), intent(in) :: R(:,:,:)       ! point symmetries in cartesian
            real(dp), intent(in) :: prec
            real(dp) :: kpoint_in_fbz(3)
            real(dp) :: ibz_kpt_cart(3)
            real(dp), allocatable :: korbit(:,:) !> orbit
            integer , allocatable :: indices(:)
            integer :: i, nsyms
            !
            nsyms = size(R,3)
            ! reduce kpt to primitive reciprocal lattice
            kpoint_in_fbz = modulo(kpt_frac+prec,1.0_dp) - prec
            ! get kpt orbit
            korbit = kpoint_orbit(R,kpoint_in_fbz)
            ! reduce korbit to primitive reciprocal lattice
            korbit = modulo(korbit+prec,1.0_dp) - prec
            ! get unique points in orbit
            korbit = unique(korbit,prec)
            ! sort
            allocate(indices(size(korbit,2)))
            do i = 1,3
                ! rank to remove numerical noise
                call rank(nint(korbit(i,:)*1.0D6),indices)
                korbit=korbit(:,indices)
            enddo
            ! get representative point
            ibz_kpt_cart = korbit(:,1)
            ! make sure representative irreducible point is in fbz
            ibz_kpt_cart = reduce_kpoint_to_fbz(kpt_cart=ibz_kpt_cart, vpt_cart=vpt_cart)
            !
        end function   reduce_kpoint_to_ibz
    end subroutine get_ibz

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


    subroutine     print_weights(bz)
        !
        implicit none
        !
        class(am_class_bz) :: bz
        integer :: i,j
        !
        call am_print('sum weight',nint(sum(bz%w)),' ... ')
        call am_print('min weight',nint(minval(bz%w)),' ... ')
        call am_print('max weight',nint(maxval(bz%w)),' ... ')
        call am_print('avg weight',sum(bz%w)/real(bz%nkpts,dp),' ... ')
        write(*,'(a5,a)') ' ... ', 'histogram of weights'
        write(*,'(5x,a10,a10)' ) 'weights', '# kpts'
        write(*,'(5x,a10,a10)' ) repeat('-',8), repeat('-',8)
        do i = minval(bz%w), maxval(bz%w)
            j = count(bz%w.eq.i)
            if (j.ne.0) write(*,'(5x,i10,i10)' ) i, j
        enddo
        !
    end subroutine print_weights

    subroutine     write_kpoints(bz,fname)
        !
        ! fname = name of file to write kpoints to.
        !
        implicit none
        class(am_class_bz), intent(in) :: bz
        character(*), intent(in) :: fname
        integer :: fid
        integer :: i
        !
        fid = 1
        !
        open(unit=fid,file=trim(fname),status='replace',action='write')
            do i = 1, bz%nkpts
                write(fid,'(i7,3f20.6,f20.6)') i, bz%kpt(:,i), bz%w(i)
            enddo
        close(fid)
        !
    end subroutine write_kpoints

    ! functions which operate on kpt

    function       locate_kpoint_in_fbz(kpt,vpt_cart,bas) result(position_in_fbz)
        !> Given a kpt in fractional coordinates, determine if ...
        !> k is inside  ( 1) of FBZ, i.e. 2KG < G^2 for ALL Bragg planes
        !> k is on edge ( 0) of FBZ, i.e. 2KG = G^2 for ONE Bragg plane
        !>                        and NOT 2KG > G^2 for ALL others
        !> k is outside (-1) of FBZ, i.e. 2KG > G^2 for ONE Bragg plane
        implicit none
        !
        real(dp), intent(in) :: kpt(3)
        real(dp), intent(in) :: bas(3,3)
        real(dp), intent(in) :: vpt_cart(3,27)
        real(dp) :: k_cart(3)
        real(dp) :: P(27)
        integer  :: position_in_fbz
        integer  :: i
        !
        ! convert kpt to cartesian
        k_cart = matmul(inv(bas),kpt)
        !
        P = 0
        do i=1,size(vpt_cart,2)
            P(i) = 2*dot_product(k_cart,vpt_cart(1:3,i)) - dot_product(vpt_cart(1:3,i),vpt_cart(1:3,i))
        enddo
        !
        ! NOTE THAT THE ORDER IS VERY IMPORTANT.
        ! If the point is outside ANY Bragg plane than it is outside the BZ
        ! The point can only be on the edge if it NOT outside. Apply the
        ! outside test last converts any points on Bragg plane, but outside the
        ! FBZ, to just pure simple outside.
        !
        ! assume inside FBZ
        position_in_fbz = 1;
        ! test if on edge of FBZ
        do i = 1,27
            if (abs(P(i)).lt.tiny) then
                position_in_fbz = 0
            endif
        enddo
        ! test if outside FBZ
        do i = 1,27
            if (abs(P(i)).gt.tiny) then
                position_in_fbz = -1
            endif
        enddo
    end function   locate_kpoint_in_fbz

    function       reduce_kpoint_to_fbz(kpt_cart,vpt_cart) result(fbz_kpt_cart)
        !> reduces kpt (in fractional) to the first Brillouin zone (Wigner-Seitz cell, defined in cartesian coordinates)
        !> cartesian kpt is returned! 
        implicit none
        !
        real(dp), intent(in) :: kpt_cart(3)    ! kpt
        real(dp), intent(in) :: vpt_cart(3,27) ! voronoi points
        real(dp) :: fbz_kpt_cart(3) !> kpt cartesian
        real(dp) :: G(3) ! reciprocal lattice vector
        real(dp) :: P    ! bragg plane condition
        integer :: i     ! loop variable
        logical :: is_not_done
        !
        ! copy kpt_cart
        fbz_kpt_cart = kpt_cart
        ! translating k-point until the closest reciprocal lattice point is [0 0 0]
        is_not_done = .true.
        do while ( is_not_done )
            is_not_done = .false.
            do i = 1, 27
                G = vpt_cart(:,i)
                P = 2*dot_product(fbz_kpt_cart,G) - dot_product(G,G)
                if ( P .gt. tiny ) then
                    fbz_kpt_cart = fbz_kpt_cart - G
                    is_not_done = .true.
                endif
            enddo
        end do
        !
    end function   reduce_kpoint_to_fbz

    pure function  get_kpoint_orbit(R,kpt) result(korbit)
        !> Note: this routine does not check to see whether the orbits are unique. If two point
        !> symmetries produce the same orbital point, this routine returns the same point twice.
        implicit none
        !
        real(dp), intent(in)  :: R(:,:,:) !> point symmetries
        real(dp), intent(in)  :: kpt(3) !> cartesian
        real(dp), allocatable :: korbit(:,:) !> orbit
        integer  :: i ! loop variable
        !
        nsyms = size(R,3)
        !
        allocate(korbit(3,nsyms))
        do i = 1, nsyms
            korbit(:,i)=matmul(R(:,:,i),kpt)
        enddo
        !
    end function   get_kpoint_orbit

    function       kpoint_weight(R,kpt,prec) result(w)
        !
        implicit none
        !
        real(dp), intent(in) :: R(:,:,:) !> point symmetries
        real(dp), intent(in) :: kpt(3) !> cartesian
        real(dp), intent(in) :: prec
        integer  :: w, nsyms
        !
        nsyms = size(R,3)
        w = size(unique(get_kpoint_orbit(R,kpt),prec),2)
        if (modulo(nsyms,w)/=0) stop 'k-point weight is not a factor of the number of symmetry operation'
        !
    end function   kpoint_weight

end module am_brillouin_zone



