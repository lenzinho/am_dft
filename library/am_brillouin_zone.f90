module am_brillouin_zone

    use am_prim_cell
    use am_symmetry
    use am_constants
    use am_stdout
    use am_vasp_io
    use am_mkl
    use am_options

    implicit none

    private

    ! BZ (primitive reciprocal-lattice cell defined with x,y,z (fractional) between [0,1).)

    type, public :: am_class_bz
        integer :: nkpts                  ! nkpts number of kpoints
        real(dp), allocatable :: kpt(:,:) ! kpt(3,nkpts) kpoint - notice that vasp IBZKPT is in fractional units.
        real(dp), allocatable :: w(:)     ! w(nkpts) normalized weights
        integer , allocatable :: fbz_id(:)
        integer , allocatable :: ibz_id(:)
    contains
        procedure :: load_ibzkpt
        procedure :: load_procar
        procedure :: load_eigenval
        procedure :: write_kpoints
    end type am_class_bz

    ! FBZ (wigner-seitz cell defined in cartesian coordinates, but kpoints are saved as fractional coordinates)

    type, public, extends(am_class_bz) :: am_class_fbz
    contains
        procedure :: get_mp
        procedure :: get_fbz
    end type am_class_fbz

    ! IBZ (irreducible wedge)

    type, public, extends(am_class_bz) :: am_class_ibz
    contains
        procedure :: get_ibz
    end type am_class_ibz

contains

    ! load from vasp

    subroutine     load_ibzkpt(bz,opts)
        ! 
        ! Reads the IBZKPT file.
        !
        use am_options
        !
        implicit none
        !
        class(am_class_bz), intent(out) :: bz
        type(am_class_options), intent(in) :: opts
        !
        call read_ibzkpt(nkpts = bz%nkpts,&
                         kpt   = bz%kpt,&
                         w     = bz%w,&
                 iopt_filename = opts%ibzkpt,&
                iopt_verbosity = opts%verbosity )
        !
    end subroutine load_ibzkpt

    subroutine     load_procar(bz,opts)
        ! 
        ! Reads the PROCAR file.
        !
        implicit none
        !
        class(am_class_bz)  , intent(out) :: bz
        type(am_class_options), intent(in) :: opts
        !
        call read_procar(nkpts    = bz%nkpts,&
                         kpt      = bz%kpt,&
                         w        = bz%w,&
                         iopt_filename = opts%procar,&
                         iopt_verbosity= 0 )
        !
    end subroutine load_procar

    subroutine     load_eigenval(bz,opts)
        ! 
        ! Reads the eigenval file.
        !
        implicit none
        !
        class(am_class_bz)  , intent(out) :: bz
        type(am_class_options), intent(in) :: opts
        !
        call read_eigenval(nkpts  = bz%nkpts,&
                           kpt    = bz%kpt,&
                           w      = bz%w,&
                           iopt_filename = opts%eigenval,&
                           iopt_verbosity = 0 )
        !
    end subroutine load_eigenval
   
    ! functions which operate on bz or similar derived types
    
    subroutine     get_mp(fbz,pc,n,s,opts)
        !
        use am_rank_and_sort
        !
        implicit none
        !
        class(am_class_fbz), intent(out) :: fbz
        type(am_class_prim_cell), intent(in) :: pc
        integer , intent(in) :: n(3)
        real(dp), intent(in) :: s(3)
        type(am_class_options), intent(in) :: opts
        real(dp) :: grid_points(3,27) !> voronoi points (27=3^3)
        integer , allocatable :: sorted_indices(:)
        real(dp), allocatable :: sort_parameter(:)
        integer :: i
        !
        if (opts%verbosity.ge.1) call am_print_title('Generating Monkhorst-Pack mesh on FBZ')
        ! get kpoints in fractional
        fbz%kpt = generate_monkhorst_pack_mesh(n=n,s=s)
        ! determine how many kpoints there are
        fbz%nkpts = size(fbz%kpt,2)
        if (opts%verbosity.ge.1) call am_print('number of kpoints',fbz%nkpts,' ... ')
        ! determine weights
        allocate(fbz%w(fbz%nkpts))
        fbz%w = 1/real(fbz%nkpts,dp)
            ! write to file for debugging
            ! call fbz%write_kpoints(fname='outfile.primitive_monkhorst_pack_frac')
        ! generate voronoi points (cartesian)
        grid_points = meshgrid([-1:1],[-1:1],[-1:1])
        grid_points = matmul(inv(pc%bas),grid_points)
        ! reduce kpoints to FBZ using voronoi points
        !$OMP PARALLEL PRIVATE(i) SHARED(pc,fbz,grid_points)
        !$OMP DO
        do i = 1,fbz%nkpts
            fbz%kpt(:,i) = reduce_kpoint_to_fbz(kpoint=fbz%kpt(:,i),grid_points=grid_points,bas=pc%bas)
        enddo
        !$OMP END DO
        !$OMP END PARALLEL
            ! write to file for debugging
            ! call fbz%write_kpoints(fname='outfile.fbz_monkhorst_pack_frac')
        ! sort
        allocate(sorted_indices(fbz%nkpts))
        allocate(sort_parameter(fbz%nkpts))
        do i = 1, 3
            sort_parameter = fbz%kpt(i,:)
            call rank(sort_parameter,sorted_indices)
            fbz%w=fbz%w(sorted_indices)
            fbz%kpt=fbz%kpt(:,sorted_indices)
        enddo
        ! write to stdout
        if (opts%verbosity.ge.1) then
            ! am_print_two_matrices_side_by_side(name, Atitle, Btitle, A, B , in_emph, iopt_fid )
            call am_print_two_matrices_side_by_side(name='kpoints',&
                Atitle='fractional',A=transpose(fbz%kpt),&
                Btitle='cartesian' ,B=transpose(matmul(inv(pc%bas),fbz%kpt)),&
                iopt_emph=' ... ',iopt_teaser=.true.)
        endif
        !
        contains
        pure function  generate_monkhorst_pack_mesh(n,s) result(kpt)
            !> returns kpoints in fractional coordinates for monkhorst pack mesh dimensions
            !> n(1:3)=[n1,n2,n3] and shift s(1:3)=[s1,s2,s3]
            !> according to http://cms.mpi.univie.ac.at/vasp/vasp/Automatic_k_mesh_generation.html
            implicit none
            !
            integer, intent(in) :: n(3) !> monkhorst pack mesh dimensions
            real(dp), intent(in) :: s(3) !> monkhorst pack mesh shift
            real(dp), allocatable :: kpt(:,:) !> kpt in fractional
            integer :: i1,i2,i3,j
            !
            allocate(kpt(3,product(n)))
            !
            j=0
            do i1=1,n(1)
                do i2=1,n(2)
                    do i3=1,n(3)
                        j=j+1
                        kpt(1:3,j) = real([i1+s(1)-1,i2+s(2)-1,i3+s(3)-1],dp)
                        kpt(1:3,j) = kpt(1:3,j)/real(n,dp)
                    enddo
                enddo
            enddo
            !
        end function   generate_monkhorst_pack_mesh
    end subroutine get_mp

    subroutine     get_fbz(fbz,bz,pc,pg,opts)
        !
        implicit none
        !
        class(am_class_fbz), intent(out) :: fbz
        class(am_class_bz) , intent(in)  :: bz
        type(am_class_prim_cell), intent(in) :: pc
        type(am_class_symmetry) , intent(in) :: pg
        type(am_class_options)  , intent(in) :: opts
        real(dp), allocatable :: korbit(:,:)
        real(dp), allocatable :: kpoints_fbz(:,:)
        real(dp), allocatable :: grid_points(:,:)
        integer :: i, j, m
        integer :: w
        !
        if (opts%verbosity.ge.1) call am_print_title('Expand kpoints to FBZ')
        if (opts%verbosity.ge.1) call am_print('number of point symmetries',pg%nsyms,' ... ')
        if (opts%verbosity.ge.1) then
            call am_print('number of point symmetries compatible with original k-point mesh',&
                size(get_kpoint_compatible_symmetries(kpt=bz%kpt,R=pg%seitz(1:3,1:3,:),prec=opts%prec),3),' ... ')
        endif
        if (opts%verbosity.ge.1) call am_print('number of kpoints in the original bz',bz%nkpts,' ... ')
        if (opts%verbosity.ge.1) then
            call am_print_two_matrices_side_by_side(name='original kpoints',&
                Atitle='fractional',A=transpose(bz%kpt),&
                Btitle='cartesian' ,B=transpose(matmul(inv(pc%bas),bz%kpt)),&
                iopt_emph=' ... ',iopt_teaser=.true.)
        endif
        !
        allocate(kpoints_fbz(3,bz%nkpts*pg%nsyms)) ! wkrspace
        m = 0
        ! expand each kpoint onto the fbz
        !$OMP PARALLEL PRIVATE(i,w,korbit) SHARED(kpoints_fbz,pg)
        !$OMP DO
        do i = 1, bz%nkpts
            ! get its orbit
            korbit = kpoint_orbit(pg%R,bz%kpt(:,i))
            ! reduce korbit to primitive reciprocal lattice
            korbit = modulo(korbit+opts%prec,1.0_dp)-opts%prec
            ! get unique values
            korbit = unique(korbit,opts%prec)
            ! get its weight
            w = size(korbit,2)
            ! check that weight is a factor of the number of symmetry operations
            if ( modulo(pg%nsyms,w) .ne. 0 ) then
                call am_print('ERROR','Weight is not a factor of the number of symmetry operation',flags='E')
                call am_print('kpoint',bz%kpt(:,i))
                call am_print('weight',w)
                call am_print('number of symmetry operations',pg%nsyms)
                call am_print('kpoint orbit',transpose(korbit))
                stop
            endif
            ! save all unique points in orbit as part of the fbz
            do j = 1,w
                m = m + 1
                kpoints_fbz(:,m) = korbit(:,j)
            enddo
        enddo
        !$OMP END DO
        !$OMP END PARALLEL
        ! make sure there are no repetitions in the fbz
        fbz%kpt = unique(kpoints_fbz(:,1:m),opts%prec)
        ! get number of kpoints in the fbz
        fbz%nkpts = size(fbz%kpt,2)
        ! reduce kpoints to FBZ using voronoi points
        grid_points = meshgrid([-1:1],[-1:1],[-1:1])
        grid_points = matmul(inv(pc%bas),grid_points)
        !$OMP PARALLEL PRIVATE(i) SHARED(fbz,grid_points,pc)
        !$OMP DO
        do i = 1,fbz%nkpts
            fbz%kpt(:,i) = reduce_kpoint_to_fbz(kpoint=fbz%kpt(:,i),grid_points=grid_points,bas=pc%bas)
        enddo
        !$OMP END DO
        !$OMP END PARALLEL
        if (opts%verbosity.ge.1) call am_print('number of kpoints in the fbz',fbz%nkpts,' ... ')
        if (opts%verbosity.ge.1) then
            call am_print_two_matrices_side_by_side(name='full kpoints',&
                Atitle='fractional',A=transpose(fbz%kpt),&
                Btitle='cartesian' ,B=transpose(matmul(inv(pc%bas),fbz%kpt)),&
                iopt_emph=' ... ',iopt_teaser=.true.)
        endif
        ! set weights
        allocate(fbz%w(fbz%nkpts))
        fbz%w = 1/real(fbz%nkpts,dp)
        ! set irreducible indices
        if (allocated(fbz%ibz_id)) deallocate(fbz%ibz_id)
        allocate(fbz%ibz_id(fbz%nkpts))
        do i = 1, fbz%nkpts
            get_indices_of_each_kpoint_in_bz : do j = 1, bz%nkpts
                if (all(abs(bz%kpt(:,j)-fbz%kpt(:,i)).lt.tiny)) then
                    fbz%ibz_id(i) = j
                    exit get_indices_of_each_kpoint_in_bz
                endif
            enddo get_indices_of_each_kpoint_in_bz
        enddo
    end subroutine get_fbz

    subroutine     get_ibz(ibz,bz,pc,pg,opts)
        !
        implicit none
        !
        class(am_class_ibz), intent(inout) :: ibz
        class(am_class_bz) , intent(inout) :: bz
        type(am_class_symmetry) , intent(in) :: pg
        type(am_class_prim_cell), intent(in) :: pc
        type(am_class_options)  , intent(in) :: opts
        real(dp), allocatable :: grid_points(:,:)
        integer :: i, j
        !
        !
        if (opts%verbosity.ge.1) call am_print_title('Reducing to IBZ')
        ! 
        if (opts%verbosity.ge.1) call am_print('number of original kpoints',bz%nkpts,' ... ')
        ! generate voronoi points (cartesian)
        grid_points = matmul( inv(pc%bas), meshgrid([-1:1],[-1:1],[-1:1]) )
        ! reduce kpoints to ibz
        allocate(ibz%kpt(3,bz%nkpts))
        !$OMP PARALLEL PRIVATE(i) SHARED(grid_points,bz,ibz,pc,pg)
        !$OMP DO
        do i = 1, bz%nkpts
            ibz%kpt(:,i) = reduce_kpoint_to_ibz(kpoint=bz%kpt(:,i),grid_points=grid_points,bas=pc%bas,R=pg%R,prec=opts%prec)
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
        !$OMP PARALLEL PRIVATE(i) SHARED(grid_points,bz,ibz,pc,pg)
        !$OMP DO
        do i = 1, ibz%nkpts
            ibz%w(i) = kpoint_weight(R=pg%R,kpoint=ibz%kpt(:,i),prec=opts%prec)
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
    end subroutine get_ibz

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

    ! functions which operate on kpoint

    function       locate_kpoint_in_fbz(kpoint,grid_points,bas) result(position_in_fbz)
        !> Given a kpoint in fractional coordinates, determine if ...
        !> k is inside  ( 1) of FBZ, i.e. 2KG < G^2 for ALL Bragg planes
        !> k is on edge ( 0) of FBZ, i.e. 2KG = G^2 for ONE Bragg plane
        !>                        and NOT 2KG > G^2 for ALL others
        !> k is outside (-1) of FBZ, i.e. 2KG > G^2 for ONE Bragg plane
        implicit none
        !
        real(dp), intent(in) :: kpoint(3)
        real(dp), intent(in) :: bas(3,3)
        real(dp), intent(in) :: grid_points(3,27)
        real(dp) :: k_cart(3)
        real(dp) :: P(27)
        integer  :: position_in_fbz
        integer  :: i
        !
        ! convert kpoint to cartesian
        k_cart = matmul(inv(bas),kpoint)
        !
        P = 0
        do i=1,size(grid_points,2)
            P(i) = 2*dot_product(k_cart,grid_points(1:3,i)) - dot_product(grid_points(1:3,i),grid_points(1:3,i))
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
    end function  locate_kpoint_in_fbz

    function       reduce_kpoint_to_fbz(kpoint,grid_points,bas) result(kpoint_fbz)
        !> reduces kpoint (in fractional) to the first Brillouin zone (Wigner-Seitz cell, defined in cartesian coordinates)
        !> cartesian kpoint is returned! 
        implicit none
        !
        real(dp), intent(in) :: kpoint(3) !> fractional
        real(dp), intent(in) :: bas(3,3) !> real space basis (column vectors)
        real(dp), intent(in) :: grid_points(3,27) !> voronoi points (cartesian)
        real(dp) :: k_cart(3) !> kpoint cartesian
        real(dp) :: G(3) !> reciprocal lattice vector
        real(dp) :: P !> bragg plane condition
        real(dp) :: kpoint_fbz(3)
        integer :: i ! loop variable
        logical :: is_not_done
        !
        ! take kpoint in fractional coordinates (will become cartesian later)
        k_cart = kpoint
        ! reduce to reciprocal unit cell
        k_cart = modulo(k_cart+tiny,1.0_dp)-tiny
        ! convert to cartesian
        k_cart = matmul(inv(bas),k_cart)
        ! reduce to Wigner-Seitz cell (aka first Brillouin zone) by translating the k-point
        ! until the closest reciprocal lattice point is [0 0 0]
        is_not_done = .true.
        do while ( is_not_done )
            is_not_done = .false.
            do i = 1,27
                G = grid_points(:,i)
                P = 2*dot_product(k_cart,G) - dot_product(G,G)
                if ( P .gt. tiny ) then
                    k_cart = k_cart - G
                    is_not_done = .true.
                endif
            enddo
        end do
        ! convert back to fractional
        kpoint_fbz = matmul(bas,k_cart)
        !
    end function  reduce_kpoint_to_fbz

    function       reduce_kpoint_to_ibz(kpoint,grid_points,bas,R,prec) result(kpoint_in_ibz)
        !
        use am_rank_and_sort
        !
        implicit none
        !
        real(dp), intent(in) :: kpoint(3) !> fractional
        real(dp), intent(in) :: bas(3,3)  !> real space basis (column vectors)
        real(dp), intent(in) :: grid_points(3,27) !> voronoi points (cartesian)
        real(dp), intent(in) :: R(:,:,:)  !> point symmetries in cartesian
        real(dp), intent(in) :: prec
        real(dp) :: kpoint_in_fbz(3)
        real(dp) :: kpoint_in_ibz(3)
        real(dp), allocatable :: korbit(:,:) !> orbit
        integer , allocatable :: indices(:)
        integer :: i, nsyms
        !
        nsyms = size(R,3)
        ! reduce kpoint to primitive reciprocal lattice
        kpoint_in_fbz = modulo(kpoint+prec,1.0_dp) - prec
        ! get kpoint orbit
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
        kpoint_in_ibz = korbit(:,1)
        ! reduce representative irreducible point to fbz
        kpoint_in_ibz = reduce_kpoint_to_fbz(kpoint=kpoint_in_ibz,grid_points=grid_points,bas=bas)
        !
    end function  reduce_kpoint_to_ibz

    pure function  kpoint_orbit(R,kpoint) result(korbit)
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
            korbit(:,i)=matmul(R(:,:,i),kpoint)
        enddo
        !
    end function  kpoint_orbit

    function       kpoint_weight(R,kpoint,prec) result(w)
        !
        implicit none
        !
        real(dp), intent(in) :: R(:,:,:) !> point symmetries
        real(dp), intent(in) :: kpoint(3) !> cartesian
        real(dp), intent(in) :: prec
        integer  :: w, nsyms
        !
        nsyms = size(R,3)
        !
        w = size(unique(kpoint_orbit(R,kpoint),prec),2)
        !
        if ( modulo(nsyms,w).ne.0) then
            call am_print('ERROR','Weight is not a factor of the number of symmetry operation',flags='E')
            call am_print('kpoint (fractional)',kpoint)
            call am_print('weight (integer)',w)
            call am_print('number of symmetry operations',nsyms)
            call am_print('kpoint orbit',transpose(unique(kpoint_orbit(R,kpoint),prec)))
            stop
        endif
        !
    end function  kpoint_weight

end module am_brillouin_zone



