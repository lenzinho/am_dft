module am_brillouin_zone
    !
    use am_unit_cell
    use am_symmetry
    use am_constants
    use am_helpers
    use am_vasp_io
    !
    implicit none
    !
    private
    !
    public :: am_class_bz
    public :: am_class_tetrahedra
    !
    ! Tetrahedra
    !
    type am_class_tetrahedra
        integer :: ntets !> ntet number f tetrahedra
        real(dp), allocatable :: volume(:) !> volume(ntet) tetrahedra volume
        integer , allocatable :: corner(:,:) !> corner(4,ntet) indices of kpoints at the corner of the tetrahedra
        real(dp), allocatable :: w(:) !> w(ntet) tetrahedra weights
    contains
        procedure :: load_ibzkpt => load_ibzkpt_tet
!        procedure :: get_tetrahedra
    end type am_class_tetrahedra

    !
    ! Brillouin Zone
    !

    type am_class_bz
        integer  :: nkpts     !> nkpts number of kpoints
        integer  :: nbands    !> nbands number of bands
        integer  :: nions     !> nions number of ions
        integer  :: norbitals !> norbitals number of orbitals
        integer  :: nspins    !> nspins number of spins (1, 2 majority/minority, or 4 dirac spinors)
        integer  :: nelecs    !> number of electrons
        real(dp), allocatable :: kpt(:,:)  !> kpt(3,nkpts) kpoint - notice that vasp IBZKPT is in fractional units.
        real(dp), allocatable :: E(:,:)    !> E(nbands,nkpts) energies
        real(dp), allocatable :: occ(:,:)  !> occ(nbands,nkpts) occupancies
        real(dp), allocatable :: w(:)      !> w(nkpts) normalized weights
        ! integer , allocatable :: w_int(:)  !> w_integer(nkpts) integer weights
        character(len=:), allocatable :: orbitals(:) !> orbitals(norbitals) names of orbitals
        real(dp), allocatable :: lmproj(:,:,:,:,:) !> lmproj(nspins,norbitals,nions,nbands,nkpts)
        integer , allocatable :: ibz_indx(:) !> index of kpoint in ibz
        ! producd by integration techniques:
        integer :: nEp
        real(dp), allocatable :: Ep(:)
        real(dp), allocatable :: pdos(:,:,:,:) !> pdos(nspins,norbitals,nions,nEprobing)
        real(dp), allocatable :: dos(:) !> dos(nEprobing)
    contains
        ! load vasp files into bz class
        procedure :: load_procar
        procedure :: load_eigenval
        procedure :: load_ibzkpt => load_ibzkpt_bz
        ! bz stuff
        procedure :: expand_to_fbz
        procedure :: reduce_to_ibz
        ! kpoints
        procedure :: full_monkhorst_pack_mesh
        procedure :: outfile_kpoints
        ! high level i/o
        procedure :: outfile_bandcharacter ! requires load_eigenval or load_procar
        procedure :: outfile_dosprojected
        !
        procedure :: tetrahedra_pdos    
    end type am_class_bz
    

contains

    !
    ! functions which operate on kpoint
    !

    pure function  locate_kpoint_in_fbz(kpoint,grid_points,bas) result(position_in_fbz)
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
        ! convert kpoint to cartesian
        k_cart = matmul(reciprocal_basis(bas),kpoint)
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
    end function   locate_kpoint_in_fbz

    pure function  reduce_kpoint_to_fbz(kpoint,grid_points,bas) result(kpoint_fbz)
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
        k_cart = matmul(reciprocal_basis(bas),k_cart)
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
    end function   reduce_kpoint_to_fbz

    function       reduce_kpoint_to_ibz(kpoint,grid_points,bas,R,sym_prec) result(kpoint_in_ibz)
        !
        use am_rank_and_sort
        !
        implicit none
        !
        real(dp), intent(in) :: kpoint(3) !> fractional
        real(dp), intent(in) :: bas(3,3)  !> real space basis (column vectors)
        real(dp), intent(in) :: grid_points(3,27) !> voronoi points (cartesian)
        real(dp), intent(in) :: R(:,:,:)  !> point symmetries in cartesian
        real(dp), intent(in) :: sym_prec
        real(dp) :: kpoint_in_fbz(3)
        real(dp) :: kpoint_in_ibz(3)
        real(dp), allocatable :: korbit(:,:) !> orbit
        integer , allocatable :: indices(:)
        integer :: i, nsyms
        !
        nsyms = size(R,3)
        ! reduce kpoint to primitive reciprocal lattice
        kpoint_in_fbz = modulo(kpoint+sym_prec,1.0_dp) - sym_prec
        ! get kpoint orbit
        korbit = kpoint_orbit(R,kpoint_in_fbz)
        ! reduce korbit to primitive reciprocal lattice
        korbit = modulo(korbit+sym_prec,1.0_dp) - sym_prec
        ! get unique points in orbit
        korbit = unique(korbit,sym_prec)
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
    end function   reduce_kpoint_to_ibz

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
    end function   kpoint_orbit

    function       kpoint_weight(R,kpoint,sym_prec) result(w)
        !
        implicit none
        !
        real(dp), intent(in) :: R(:,:,:) !> point symmetries
        real(dp), intent(in) :: kpoint(3) !> cartesian
        real(dp), intent(in) :: sym_prec
        integer  :: w, nsyms
        !
        nsyms = size(R,3)
        !
        w = size(unique(kpoint_orbit(R,kpoint),sym_prec),2)
        !
        if ( modulo(nsyms,w).ne.0) then
            call am_print('ERROR','Weight is not a factor of the number of symmetry operation',' >>> ')
            call am_print('kpoint (fractional)',kpoint)
            call am_print('weight (integer)',w)
            call am_print('number of symmetry operations',nsyms)
            call am_print('kpoint orbit',transpose(unique(kpoint_orbit(R,kpoint),sym_prec)))
            stop
        endif
        !
    end function   kpoint_weight

    !
    ! functions which operate/generate on meshes/grids
    !

    pure function  mesh_monkhorstpack_grid(n,s) result(kpt)
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
    end function   mesh_monkhorstpack_grid

    !
    ! functions which operate on bz or similar derived types
    ! 

    subroutine     load_procar(bz,opts)
        ! 
        ! Reads the PROCAR file, look below: code is short and self explanatory.
        !
        use am_options
        use am_vasp_io
        use am_unit_cell
        !
        implicit none
        !
        class(am_class_bz), intent(inout) :: bz
        type(am_class_options), intent(in) :: opts
        !
        call read_procar(nkpts=bz%nkpts,&
            nbands=bz%nbands,&
            nions=bz%nions,&
            norbitals=bz%norbitals,&
            nspins=bz%nspins,&
            E=bz%E,&
            occ=bz%occ,&
            kpt=bz%kpt,&
            w=bz%w,&
            orbitals=bz%orbitals,&
            lmproj=bz%lmproj,&
            iopt_filename=opts%procar,&
            iopt_verbosity=opts%verbosity)
        !
    end subroutine load_procar

    subroutine     load_eigenval(bz,opts)
        ! 
        ! Reads the eigenval file, look below: code is short and self explanatory.
        !
        use am_vasp_io
        use am_options
        use am_unit_cell
        !
        implicit none
        !
        class(am_class_bz), intent(inout) :: bz
        type(am_class_options), intent(in) :: opts
        !
        call read_eigenval(nkpts=bz%nkpts,&
            nbands=bz%nbands,&
            nspins=bz%nspins,&
            nelecs=bz%nelecs,&
            kpt=bz%kpt,&
            w=bz%w,&
            E=bz%E,&
            lmproj=bz%lmproj,&
            iopt_filename=opts%eigenval,&
            iopt_verbosity=opts%verbosity)
        !
        bz%nions=1
        bz%norbitals=1
        allocate(character(15) :: bz%orbitals(1))
        bz%orbitals='tot'
        !
    end subroutine load_eigenval

    subroutine     load_ibzkpt_bz(bz,opts)
        ! 
        ! Reads the IBZKPT file, look below: code is short and self explanatory.
        !
        use am_options
        !
        implicit none
        !
        class(am_class_bz), intent(inout) :: bz
        type(am_class_options), intent(in) :: opts
        type(am_class_tetrahedra) :: tet
        !
        call read_ibzkpt(&
            nkpts = bz%nkpts,&
            kpt   = bz%kpt,&
            w     = bz%w,&
            ntets = tet%ntets,&
            vtet  = tet%volume,&
            tet   = tet%corner,&
            wtet  = tet%w,&
            iopt_filename = opts%ibzkpt,&
            iopt_verbosity = opts%verbosity )
        !
    end subroutine load_ibzkpt_bz

    subroutine     full_monkhorst_pack_mesh(fbz,uc,n,s,opts)
        !
        use am_rank_and_sort
        use am_options
        !
        implicit none
        !
        class(am_class_bz), intent(inout) :: fbz !> brillouin zone class
        type(am_class_unit_cell), intent(in) :: uc
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
        fbz%kpt = mesh_monkhorstpack_grid(n=n,s=s)
        ! determine how many kpoints there are
        fbz%nkpts = size(fbz%kpt,2)
        if (opts%verbosity.ge.1) call am_print('number of kpoints',fbz%nkpts,' ... ')
        ! determine weights
        allocate(fbz%w(fbz%nkpts))
        fbz%w = 1/real(fbz%nkpts,dp)
            ! write to file for debugging
            ! call fbz%outfile_kpoints(fname='outfile.primitive_monkhorst_pack_frac')
        ! generate voronoi points (cartesian)
        grid_points = mesh_grid([1,1,1])
        grid_points = matmul(reciprocal_basis(uc%bas),grid_points)
        ! reduce kpoints to FBZ using voronoi points
        !$OMP PARALLEL PRIVATE(i) SHARED(uc,fbz,grid_points)
        !$OMP DO
        do i = 1,fbz%nkpts
            fbz%kpt(:,i) = reduce_kpoint_to_fbz(kpoint=fbz%kpt(:,i),grid_points=grid_points,bas=uc%bas)
        enddo
        !$OMP END DO
        !$OMP END PARALLEL
            ! write to file for debugging
            ! call fbz%outfile_kpoints(fname='outfile.fbz_monkhorst_pack_frac')
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
            ! am_print_two_matrices_side_by_side(name, Atitle, Btitle, A, B , in_emph, in_fid )
            call am_print_two_matrices_side_by_side(name='kpoints',&
                Atitle='fractional',A=transpose(fbz%kpt),&
                Btitle='cartesian' ,B=transpose(matmul(reciprocal_basis(uc%bas),fbz%kpt)),&
                iopt_emph=' ... ',iopt_teaser=.true.)
        endif
        !
    end subroutine full_monkhorst_pack_mesh

    subroutine     reduce_to_ibz(ibz,bz,uc,pg,opts)
        !
        use am_options
        !
        implicit none
        !
        class(am_class_bz), intent(inout) :: ibz
        type(am_class_bz) , intent(inout) :: bz
        type(am_class_symmetry) , intent(in) :: pg
        type(am_class_unit_cell), intent(in) :: uc
        type(am_class_options)  , intent(in) :: opts
        real(dp), allocatable :: grid_points(:,:)
        integer :: i, j
        !
        !
        if (opts%verbosity.ge.1) call am_print_title('Reducing to IBZ')
        ! 
        if (opts%verbosity.ge.1) call am_print('number of original kpoints',bz%nkpts,' ... ')
        ! generate voronoi points (cartesian)
        grid_points = matmul( reciprocal_basis(uc%bas), mesh_grid([1,1,1]) )
        ! reduce kpoints to ibz
        allocate(ibz%kpt(3,bz%nkpts))
        !$OMP PARALLEL PRIVATE(i) SHARED(grid_points,bz,ibz,uc,pg)
        !$OMP DO
        do i = 1, bz%nkpts
            ibz%kpt(:,i) = reduce_kpoint_to_ibz(kpoint=bz%kpt(:,i),grid_points=grid_points,bas=uc%bas,R=pg%R,sym_prec=opts%sym_prec)
        enddo
        !$OMP END DO
        !$OMP END PARALLEL
        ! get unique points
        ibz%kpt = unique(ibz%kpt,opts%sym_prec)
        ! get number of unique points
        ibz%nkpts = size(ibz%kpt,2)
        if (opts%verbosity.ge.1) call am_print('number of irreducible kpoints',ibz%nkpts,' ... ')
        ! get weights
        allocate(ibz%w(ibz%nkpts))
        !$OMP PARALLEL PRIVATE(i) SHARED(grid_points,bz,ibz,uc,pg)
        !$OMP DO
        do i = 1, ibz%nkpts
            ibz%w(i) = kpoint_weight(R=pg%R,kpoint=ibz%kpt(:,i),sym_prec=opts%sym_prec)
        enddo
        !$OMP END DO
        !$OMP END PARALLEL
        if (opts%verbosity.ge.1) then
            if (.not.nint(sum(ibz%w)).eq.bz%nkpts) then
                call am_print('WARNING','The sum of irreducible weights does not equal the number of k-points in the full zone.',' >>> ')
                ! This could mean that:
                ! 1) there was an internal error. this has been occuring when an even MP mesh is requested
                ! 2) the full bz mesh is not a MP mesh
                ! 3) the MP mesh contains different dimensions along x,y,z, i.e. 3 3 5  
            endif
        endif
        if (opts%verbosity.ge.1) call stdout_weights(ibz)
        ! normalize weights
        ibz%w = ibz%w/real(sum(ibz%w),dp)
        ! set irreducible indices
        if (.not.allocated(bz%ibz_indx)) allocate(bz%ibz_indx(bz%nkpts))
        irrk : do i = 1,ibz%nkpts
               do j = 1, bz%nkpts
                   if ( all(abs(bz%kpt(:,j) - ibz%kpt(:,i)).lt.tiny) ) then
                       bz%ibz_indx(j) = i
                       cycle irrk
                   endif
               enddo
        enddo irrk
        ! stdout
        if (opts%verbosity.ge.1) then
            call am_print_two_matrices_side_by_side(name='irreducible kpoints',&
                Atitle='fractional',A=transpose(ibz%kpt),&
                Btitle='cartesian' ,B=transpose(matmul(reciprocal_basis(uc%bas),ibz%kpt)),&
                iopt_emph=' ... ',iopt_teaser=.true.)
        endif
        !
    end subroutine reduce_to_ibz

    subroutine     expand_to_fbz(fbz,bz,uc,pg,opts)
        !
        use am_options
        !
        implicit none
        !
        class(am_class_bz), intent(inout) :: fbz
        type(am_class_bz) , intent(in) :: bz
        type(am_class_unit_cell), intent(in) :: uc
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
                size(get_kpoint_compatible_symmetries(kpt=bz%kpt,R=pg%R,sym_prec=opts%sym_prec),3),' ... ')
        endif
        if (opts%verbosity.ge.1) call am_print('number of kpoints in the original bz',bz%nkpts,' ... ')
        if (opts%verbosity.ge.1) then
            call am_print_two_matrices_side_by_side(name='original kpoints',&
                Atitle='fractional',A=transpose(bz%kpt),&
                Btitle='cartesian' ,B=transpose(matmul(reciprocal_basis(uc%bas),bz%kpt)),&
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
            korbit = modulo(korbit+opts%sym_prec,1.0_dp)-opts%sym_prec
            ! get unique values
            korbit = unique(korbit,opts%sym_prec)
            ! get its weight
            w = size(korbit,2)
            ! check that weight is a factor of the number of symmetry operations
            if ( modulo(pg%nsyms,w) .ne. 0 ) then
                call am_print('ERROR','Weight is not a factor of the number of symmetry operation',' >>> ')
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
        fbz%kpt = unique(kpoints_fbz(:,1:m),opts%sym_prec)
        ! get number of kpoints in the fbz
        fbz%nkpts = size(fbz%kpt,2)
        ! reduce kpoints to FBZ using voronoi points
        grid_points = mesh_grid([1,1,1])
        grid_points = matmul(reciprocal_basis(uc%bas),grid_points)
        !$OMP PARALLEL PRIVATE(i) SHARED(fbz,grid_points,uc)
        !$OMP DO
        do i = 1,fbz%nkpts
            fbz%kpt(:,i) = reduce_kpoint_to_fbz(kpoint=fbz%kpt(:,i),grid_points=grid_points,bas=uc%bas)
        enddo
        !$OMP END DO
        !$OMP END PARALLEL
        if (opts%verbosity.ge.1) call am_print('number of kpoints in the fbz',fbz%nkpts,' ... ')
        if (opts%verbosity.ge.1) then
            call am_print_two_matrices_side_by_side(name='full kpoints',&
                Atitle='fractional',A=transpose(fbz%kpt),&
                Btitle='cartesian' ,B=transpose(matmul(reciprocal_basis(uc%bas),fbz%kpt)),&
                iopt_emph=' ... ',iopt_teaser=.true.)
        endif
        ! set weights
        allocate(fbz%w(fbz%nkpts))
        fbz%w = 1/real(fbz%nkpts,dp)
        ! set irreducible indices
        if (allocated(fbz%ibz_indx)) deallocate(fbz%ibz_indx)
        allocate(fbz%ibz_indx(fbz%nkpts))
        do i = 1, fbz%nkpts
            get_indices_of_each_kpoint_in_bz : do j = 1, bz%nkpts
                if (all(abs(bz%kpt(:,j)-fbz%kpt(:,i)).lt.tiny)) then
                    fbz%ibz_indx(i) = j
                    exit get_indices_of_each_kpoint_in_bz
                endif
            enddo get_indices_of_each_kpoint_in_bz
        enddo
    end subroutine expand_to_fbz
    
    subroutine     outfile_bandcharacter(bz,uc,pg,opts)
        ! 
        ! Creates outfile.bandcharacter, which contains the energy and character of the bands read either from EIGENVAL using bz%load_eigenval() or PROCAR using bz%load_procar()
        !
        use am_options
        !
        implicit none
        !
        class(am_class_bz), intent(in) :: bz !> brillouin zone class
        type(am_class_unit_cell), intent(in) :: uc
        type(am_class_symmetry) , intent(in) :: pg
        type(am_class_options), intent(in) :: opts
        real(dp), allocatable :: grid_points(:,:) !> voronoi points (27=3^3)
        real(dp) :: k1(3) ! point for kdiff
        real(dp) :: k2(3) ! point for kdiff
        real(dp) :: kdiff ! for dispersion plot
        integer :: fid
        integer :: i, j, l, m, n ! loop variables
        !
        if (opts%verbosity.ge.1) call am_print_title('Writing bandcharacter file')
        !
        !
        ! generate voronoi points (cartesian), used to reduce points to wigner-seitz brillouin zone
        grid_points = matmul( reciprocal_basis(uc%bas) , real(mesh_grid([1,1,1]),dp) )
        !
        fid = 1
        open(unit=fid,file="outfile.bandcharacter",status="replace",action='write')
            !
            ! STDOUT
            !
            if (opts%verbosity.ge.1) call am_print('number of kpoints', bz%nkpts,' ... ')
            if (opts%verbosity.ge.1) call am_print('number of bands', bz%nbands ,' ... ')
            if (opts%verbosity.ge.1) call am_print('number of ions', bz%nions ,' ... ')
            if (opts%verbosity.ge.1) call am_print('number of orbitals', bz%norbitals ,' ... ')
            if (opts%verbosity.ge.1) call am_print('number of spins', bz%nspins ,' ... ')
            if (opts%verbosity.ge.1) call am_print('number of columns', bz%nspins*bz%norbitals ,' ... ')
            !
            ! HEADER (first three lines)
            !
            write(fid,'(a)')  'Character-projected energy dispersion'
            write(fid,'(100a13)',advance='no') ' ', ' ',' ',' ','ions'
            do m = 1, bz%nions
                do l = 1, bz%norbitals
                do n = 1, bz%nspins
                    write(fid,'(i13)' ,advance='no') m
                enddo
                enddo
            enddo
            write(fid,*)
            !
            write(fid,'(100a13)',advance='no') ' ', ' ',' ',' ','spin'
            do m = 1, bz%nions
                do l = 1, bz%norbitals
                do n = 1, bz%nspins
                    write(fid,'(i13)' ,advance='no') n
                enddo
                enddo
            enddo
            write(fid,*)
            !
            write(fid,'(100a13)',advance='no') 'kpath','k_x','k_y','k_z','E [eV]'
            do m = 1, bz%nions
                do l = 1, bz%norbitals
                do n = 1, bz%nspins
                    write(fid,'(a13)',advance='no') trim(bz%orbitals(l))
                enddo
                enddo
            enddo 
            write(fid,*)
            !
            ! DATA
            !
            do j = 1, bz%nbands
                do i = 1, bz%nkpts
                    if (i .eq. 1) then
                        kdiff = 0
                    else
                        k1=bz%kpt(:,i-1)
                        k2=bz%kpt(:,i)
                        ! k1 = reduce_kpoint_to_ibz(kpoint=bz%kpt(:,i-1),grid_points=grid_points,bas=uc%bas,R=pg%R)
                        ! k2 = reduce_kpoint_to_ibz(kpoint=bz%kpt(:,i),grid_points=grid_points,bas=uc%bas,R=pg%R)
                        kdiff = kdiff + norm2( abs(k1-k2) )
                    endif
                    !
                    write(fid,'(100f13.5)',advance='no') kdiff, bz%kpt(1:3,i), bz%E(j,i)
                    do m = 1, bz%nions
                        do l = 1, bz%norbitals
                        do n = 1, bz%nspins
                            !> lmproj(nspins,norbitals,nions,nbands,nkpts)
                            write(fid,'(f13.5)',advance='no') bz%lmproj(n,l,m,j,i)
                        enddo
                        enddo
                    enddo
                    write(fid,*)
                    !
                enddo
                write(fid,'(a)') ' '
            enddo
        close(fid)
        !
    end subroutine outfile_bandcharacter

    subroutine     outfile_dosprojected(bz,opts)
        ! 
        ! Creates outfile.dosprojected
        !
        use am_options
        !
        implicit none
        !
        class(am_class_bz), intent(inout) :: bz !> brillouin zone class
        type(am_class_options), intent(in) :: opts
        integer :: fid
        integer :: i, l, m, n ! loop variables
        !
        if (opts%verbosity.ge.1) call am_print_title('Writing dosprojected file')
        !
        fid = 1
        open(unit=fid,file="outfile.dosprojected",status="replace",action='write')
            !
            ! STDOUT
            !
            if (opts%verbosity.ge.1) write(*,*) ' ... Projection information'
            if (opts%verbosity.ge.1) call am_print('number of ions', bz%nions ,' ... ')
            if (opts%verbosity.ge.1) call am_print('number of orbitals', bz%norbitals ,' ... ')
            if (opts%verbosity.ge.1) call am_print('number of spins', bz%nspins ,' ... ')
            if (opts%verbosity.ge.1) call am_print('number of columns', bz%nspins*bz%norbitals ,' ... ')
            !
            ! HEADER (first four lines)
            !
            write(fid,'(a)')  'Projected density of states [States/eV/spin/unit-cell]'
            write(fid,'(100a13)',advance='no')  'ions'
            do m = 1, bz%nions
                do l = 1, bz%norbitals
                do n = 1, bz%nspins
                    write(fid,'(i13)' ,advance='no') m
                enddo
                enddo
            enddo
            write(fid,*)
            !
            write(fid,'(100a13)',advance='no') 'spin'
            do m = 1, bz%nions
                do l = 1, bz%norbitals
                do n = 1, bz%nspins
                    write(fid,'(i13)' ,advance='no') n
                enddo
                enddo
            enddo
            write(fid,*)
            !
            write(fid,'(100a13)',advance='no') 'E [eV]'
            do m = 1, bz%nions
                do l = 1, bz%norbitals
                do n = 1, bz%nspins
                    write(fid,'(a13)',advance='no') trim(bz%orbitals(l))
                enddo
                enddo
            enddo 
            write(fid,*)
            !
            ! DATA
            !
            !real(dp), allocatable :: Ep(:)
            !real(dp), allocatable :: pdos(:,:,:,:) !> pdos(nspins,norbitals,nions,nEprobing)
            !real(dp), allocatable :: dos(:) !> dos(nEprobing)
            
            do i = 1, bz%nEp
                !
                write(fid,'(f13.5)',advance='no') bz%Ep(i)
                do m = 1, bz%nions
                    do l = 1, bz%norbitals
                        do n = 1, bz%nspins
                            !> lmproj(nspins,norbitals,nions,nbands,nkpts)
                            write(fid,'(f13.5)',advance='no') bz%pdos(n,l,m,i)
                        enddo
                    enddo
                enddo
                write(fid,*)
                !
            enddo
        close(fid)
        !
    end subroutine outfile_dosprojected

    subroutine     stdout_weights(bz)
        !
        implicit none
        !
        type(am_class_bz) :: bz
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
    end subroutine stdout_weights

    subroutine     outfile_kpoints(bz,fname)
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
    end subroutine outfile_kpoints

    !
    ! functions which operate on tetrahedra connectivity lists
    !

    function       tetrahedron_weights(x,xc,integration_type,apply_blochl_corrections) result(w)
        !
        ! integration_type = heavi/delta
        !
        use am_rank_and_sort
        !
        implicit none
        !
        real(dp), intent(in) :: xc(4)
        real(dp) :: w(4)
        character(len=5), intent(in) :: integration_type
        logical, intent(in) :: apply_blochl_corrections
        real(dp) :: x,x1,x2,x3,x4,f
        real(dp) :: xi(4)
        integer :: bracket
        integer :: i
        !
        ! sort xc into xi
        xi = xc([4:1:-1])
        if (xi(1).gt.xi(2)) then; xi([1,2]) = xi([2,1]); endif
        if (xi(3).gt.xi(4)) then; xi([3,4]) = xi([4,3]); endif
        if (xi(1).gt.xi(3)) then; xi([1,3]) = xi([3,1]); endif
        if (xi(2).gt.xi(4)) then; xi([2,4]) = xi([4,2]); endif
        if (xi(2).gt.xi(3)) then; xi([2,3]) = xi([3,2]); endif
        x1=xi(1); x2=xi(2); x3=xi(3); x4=xi(4)
        !
        ! if (x1.gt.x2) then; call am_print('ERROR','Sorting failed.',' >>> '); stop; endif
        ! if (x2.gt.x3) then; call am_print('ERROR','Sorting failed.',' >>> '); stop; endif
        ! if (x3.gt.x4) then; call am_print('ERROR','Sorting failed.',' >>> '); stop; endif
        !
        bracket = 0
        if     (x.lt.x1) then; bracket = 1 !; call am_print('bracket',xi-x)
        elseif (x.lt.x2) then; bracket = 2 !; call am_print('bracket',xi-x)
        elseif (x.lt.x3) then; bracket = 3 !; call am_print('bracket',xi-x)
        elseif (x.lt.x4) then; bracket = 4 !; call am_print('bracket',xi-x)
        elseif (x.gt.x4) then; bracket = 5 !; call am_print('bracket',xi-x)
        else
            call am_print('ERROR','Unable to bracket.',' >>> ')
            stop
        endif
        !
        select case(integration_type)
        case('heavi')
            select case(bracket)
                case (1)
                    w(1) = 0.0_dp
                    w(2) = 0.0_dp
                    w(3) = 0.0_dp
                    w(4) = 0.0_dp
                    f    = 0.0_dp
                case (2)
                    w(1) = -((x-x1)**3*((x-x1)*(1.0_dp/(x1-x2)+1.0_dp/(x1-x3)+1.0_dp/(x1-x4))+4.0_dp))/((x1-x2)*(x1-x3)*(x1-x4))
                    w(2) =  ((x-x1)**4*1.0_dp/(x1-x2)**2)/((x1-x3)*(x1-x4))
                    w(3) =  ((x-x1)**4*1.0_dp/(x1-x3)**2)/((x1-x2)*(x1-x4))
                    w(4) =  ((x-x1)**4*1.0_dp/(x1-x4)**2)/((x1-x2)*(x1-x3))
                    f    =  ((x-x1)**2*(-1.2D1))/((x1-x2)*(x1-x3)*(x1-x4))
                case (3)
                    w(1) =  (x-x1)**2/((x1-x3)*(x1-x4))+(((x-x1)**2/((x1-x3)*(x1-x4))+((x-x1)*(x-x2)*(x-x3))/((x1-x3)*(x1-x4)*(x2-x3)))*(x-x3))/(x1-x3)+((x-x4)*((x-x1)**2/((x1-x3)*(x1-x4))+((x-x2)**2*(x-x4))/((x1-x4)*(x2-x3)*(x2-x4))+((x-x1)*(x-x2)*(x-x3))/((x1-x3)*(x1-x4)*(x2-x3))))/(x1-x4)
                    w(2) =  ((((x-x2)**2*(x-x4))/((x1-x4)*(x2-x3)*(x2-x4))+((x-x1)*(x-x2)*(x-x3))/((x1-x3)*(x1-x4)*(x2-x3)))*(x-x3))/(x2-x3)+(x-x1)**2/((x1-x3)*(x1-x4))+((x-x2)**2*(x-x4))/((x1-x4)*(x2-x3)*(x2-x4))+((x-x2)**2*(x-x4)**2*1.0_dp/(x2-x4)**2)/((x1-x4)*(x2-x3))+((x-x1)*(x-x2)*(x-x3))/((x1-x3)*(x1-x4)*(x2-x3))
                    w(3) = -((((x-x2)**2*(x-x4))/((x1-x4)*(x2-x3)*(x2-x4))+((x-x1)*(x-x2)*(x-x3))/((x1-x3)*(x1-x4)*(x2-x3)))*(x-x2))/(x2-x3)-(((x-x1)**2/((x1-x3)*(x1-x4))+((x-x1)*(x-x2)*(x-x3))/((x1-x3)*(x1-x4)*(x2-x3)))*(x-x1))/(x1-x3)
                    w(4) = -((x-x1)*((x-x1)**2/((x1-x3)*(x1-x4))+((x-x2)**2*(x-x4))/((x1-x4)*(x2-x3)*(x2-x4))+((x-x1)*(x-x2)*(x-x3))/((x1-x3)*(x1-x4)*(x2-x3))))/(x1-x4)-((x-x2)**3*(x-x4)*1.0_dp/(x2-x4)**2)/((x1-x4)*(x2-x3))
                    f    =  (x*2.4D1-x1*1.2D1-x2*1.2D1+((x-x2)**2*(x1*3.0_dp+x2*3.0_dp-x3*3.0_dp-x4*3.0_dp)*4.0_dp)/((x2-x3)*(x2-x4)))/((x1-x3)*(x1-x4))
                case (4)
                    w(1) = -((x-x4)**4*1.0_dp/(x1-x4)**2)/((x2-x4)*(x3-x4))+1.0_dp
                    w(2) = -((x-x4)**4*1.0_dp/(x2-x4)**2)/((x1-x4)*(x3-x4))+1.0_dp
                    w(3) = -((x-x4)**4*1.0_dp/(x3-x4)**2)/((x1-x4)*(x2-x4))+1.0_dp
                    w(4) =  ((x-x4)**3*((x-x4)*(1.0_dp/(x1-x4)+1.0_dp/(x2-x4)+1.0_dp/(x3-x4))-4.0_dp))/((x1-x4)*(x2-x4)*(x3-x4))+1.0_dp
                    f    =  ((x-x4)**2*(-1.2D1))/((x1-x4)*(x2-x4)*(x3-x4))
                case (5)
                    w(1) = 1.0_dp
                    w(2) = 1.0_dp
                    w(3) = 1.0_dp
                    w(4) = 1.0_dp
                    f    = 0.0_dp
                case default
                    call am_print('ERROR','Error computing tetrahedron weight.',' >>> ')
                    call am_print('x',x)
                    call am_print('xc',[x1,x2,x3,x4])
                    call am_print('w',w)
                    call am_print('f',f)
                    stop                        
                end select
        case('delta')
            select case(bracket)
                case (1)
                    w(1) = 0.0_dp
                    w(2) = 0.0_dp
                    w(3) = 0.0_dp
                    w(4) = 0.0_dp
                    f    = 0.0_dp
                case (2)
                    w(1) =   (x-x1)**2*1.0_dp/(x1-x2)**2*1.0_dp/(x1-x3)**2*1.0_dp/(x1-x4)**2*(x*(x4*(x2+x3)+x2*x3-x1*(x2+x3+x4)*2.0_dp+x1**2*3.0_dp)+x1*(x4*(x2+x3)+x2*x3)*2.0_dp-x1**2*(x2+x3+x4)-x2*x3*x4*3.0_dp)*(-4.0_dp)
                    w(2) =  ((x-x1)**3*1.0_dp/(x1-x2)**2*4.0_dp)/((x1-x3)*(x1-x4))
                    w(3) =  ((x-x1)**3*1.0_dp/(x1-x3)**2*4.0_dp)/((x1-x2)*(x1-x4))
                    w(4) =  ((x-x1)**3*1.0_dp/(x1-x4)**2*4.0_dp)/((x1-x2)*(x1-x3))
                    f    = -(x*2.4D1-x1*2.4D1)/((x1-x2)*(x1-x3)*(x1-x4))
                case (3)
                    w(1) = (1.0_dp/(x1-x3)**2*1.0_dp/(x1-x4)**2*(x**2*(x1**2*x2*3.0_dp+x3*x4**2*3.0_dp+x3**2*x4*3.0_dp-x1*x3*x4*6.0_dp-x2*x3*x4*3.0_dp)-x**3*(x1*x2*2.0_dp-x1*x3*2.0_dp-x1*x4*2.0_dp-x2*x3-x2*x4+x3*x4+x1**2+x3**2+x4**2)-x*(x3**2*x4**2*3.0_dp+x1**2*x2*x3*3.0_dp+x1**2*x2*x4*3.0_dp-x1**2*x3*x4*3.0_dp-x1*x2*x3*x4*6.0_dp)+x1**2*x2*x3**2+x1**2*x2*x4**2+x1*x3**2*x4**2*2.0_dp-x1**2*x3*x4**2-x1**2*x3**2*x4+x2*x3**2*x4**2-x1*x2*x3*x4**2*2.0_dp-x1*x2*x3**2*x4*2.0_dp+x1**2*x2*x3*x4)*(-4.0_dp))/((x2-x3)*(x2-x4))
                    w(2) = (1.0_dp/(x2-x3)**2*1.0_dp/(x2-x4)**2*(x**2*(x1*x2**2*3.0_dp+x3*x4**2*3.0_dp+x3**2*x4*3.0_dp-x1*x3*x4*3.0_dp-x2*x3*x4*6.0_dp)-x**3*(x1*x2*2.0_dp-x1*x3-x1*x4-x2*x3*2.0_dp-x2*x4*2.0_dp+x3*x4+x2**2+x3**2+x4**2)-x*(x3**2*x4**2*3.0_dp+x1*x2**2*x3*3.0_dp+x1*x2**2*x4*3.0_dp-x2**2*x3*x4*3.0_dp-x1*x2*x3*x4*6.0_dp)+x1*x2**2*x3**2+x1*x2**2*x4**2+x1*x3**2*x4**2+x2*x3**2*x4**2*2.0_dp-x2**2*x3*x4**2-x2**2*x3**2*x4-x1*x2*x3*x4**2*2.0_dp-x1*x2*x3**2*x4*2.0_dp+x1*x2**2*x3*x4)*(-4.0_dp))/((x1-x3)*(x1-x4))
                    w(3) = (1.0_dp/(x1-x3)**2*1.0_dp/(x2-x3)**2*(x**2*(x1*x2**2*3.0_dp+x1**2*x2*3.0_dp+x3**2*x4*3.0_dp-x1*x2*x3*6.0_dp-x1*x2*x4*3.0_dp)-x**3*(x1*x2-x1*x3*2.0_dp-x1*x4-x2*x3*2.0_dp-x2*x4+x3*x4*2.0_dp+x1**2+x2**2+x3**2)-x*(x1**2*x2**2*3.0_dp-x1*x2*x3**2*3.0_dp+x1*x3**2*x4*3.0_dp+x2*x3**2*x4*3.0_dp-x1*x2*x3*x4*6.0_dp)-x1*x2**2*x3**2-x1**2*x2*x3**2+x1**2*x2**2*x3*2.0_dp+x1**2*x2**2*x4+x1**2*x3**2*x4+x2**2*x3**2*x4+x1*x2*x3**2*x4-x1*x2**2*x3*x4*2.0_dp-x1**2*x2*x3*x4*2.0_dp)*4.0_dp)/((x1-x4)*(x2-x4))
                    w(4) = (1.0_dp/(x1-x4)**2*1.0_dp/(x2-x4)**2*(x**2*(x1*x2**2*3.0_dp+x1**2*x2*3.0_dp+x3*x4**2*3.0_dp-x1*x2*x3*3.0_dp-x1*x2*x4*6.0_dp)-x**3*(x1*x2-x1*x3-x1*x4*2.0_dp-x2*x3-x2*x4*2.0_dp+x3*x4*2.0_dp+x1**2+x2**2+x4**2)-x*(x1**2*x2**2*3.0_dp-x1*x2*x4**2*3.0_dp+x1*x3*x4**2*3.0_dp+x2*x3*x4**2*3.0_dp-x1*x2*x3*x4*6.0_dp)+x1**2*x2**2*x3-x1*x2**2*x4**2-x1**2*x2*x4**2+x1**2*x2**2*x4*2.0_dp+x1**2*x3*x4**2+x2**2*x3*x4**2+x1*x2*x3*x4**2-x1*x2**2*x3*x4*2.0_dp-x1**2*x2*x3*x4*2.0_dp)*4.0_dp)/((x1-x3)*(x2-x3))
                    f    = (x1*x2*(-2.4D1)+x3*x4*2.4D1+x*(x1+x2-x3-x4)*2.4D1)/((x1-x3)*(x1-x4)*(x2-x3)*(x2-x4))
                case (4)
                    w(1) =  ((x-x4)**3*1.0_dp/(x1-x4)**2*(-4.0_dp))/((x2-x4)*(x3-x4))
                    w(2) =  ((x-x4)**3*1.0_dp/(x2-x4)**2*(-4.0_dp))/((x1-x4)*(x3-x4))
                    w(3) =  ((x-x4)**3*1.0_dp/(x3-x4)**2*(-4.0_dp))/((x1-x4)*(x2-x4))
                    w(4) =   (x-x4)**2*1.0_dp/(x1-x4)**2*1.0_dp/(x2-x4)**2*1.0_dp/(x3-x4)**2*(x*(x1*x2+x1*x3+x2*x3-x4*(x1*2.0_dp+x2*2.0_dp+x3*2.0_dp)+x4**2*3.0_dp)+x4*(x1*(x2+x3)*2.0_dp+x2*x3*2.0_dp)-x4**2*(x1+x2+x3)-x1*x2*x3*3.0_dp)*4.0_dp
                    f    = -(x*2.4D1-x4*2.4D1)/((x1-x4)*(x2-x4)*(x3-x4))
                case (5)
                    w(1) = 0.0_dp
                    w(2) = 0.0_dp
                    w(3) = 0.0_dp
                    w(4) = 0.0_dp
                    f    = 0.0_dp
                case default
                    call am_print('ERROR','Error computing tetrahedron weight.',' >>> ')
                    call am_print('x',x)
                    call am_print('xc',[x1,x2,x3,x4])
                    call am_print('w',w)
                    call am_print('f',f)
                    stop                        
                end select
            case default
                call am_print('ERROR','Integration type not valid',' >>> ')
                stop
        end select
        if (any(isnan(w))) then
            call am_print('ERROR','NaN weight returned.',' >>> ')
            call am_print('integration_type',integration_type)
            call am_print('bracket',bracket)
            call am_print('xc',[x1,x2,x3,x4])
            call am_print('x',x)
            call am_print('w',w)
            stop
        endif
        if ( apply_blochl_corrections ) then
            f = f*0.25_dp
            do i = 1,4
                w(1)=w(1)+0.025_dp*f*(xc(i)-x1)
                w(2)=w(2)+0.025_dp*f*(xc(i)-x2)
                w(3)=w(3)+0.025_dp*f*(xc(i)-x3)
                w(4)=w(4)+0.025_dp*f*(xc(i)-x4)
            enddo
        endif
        !
        w=w*0.25_dp
        !
    end function   tetrahedron_weights

    !
    ! functions which operate on tet or similar derived types
    !

    subroutine     load_ibzkpt_tet(tet,iopt_bz,opts)
        ! 
        ! Reads the IBZKPT file, look below: code is short and self explanatory.
        !
        use am_options
        !
        implicit none
        !
        class(am_class_tetrahedra), intent(inout) :: tet
        type(am_class_bz), intent(in), optional :: iopt_bz
        type(am_class_options), intent(in) :: opts
        type(am_class_bz) :: bz_tmp
        !
        ! tetrahedra parameters:
        ! ntets number f tetrahedra
        ! volume(ntet) tetrahedra volume
        ! corner(4,ntet) indices of kpoints at the corner of the 
        ! w(ntet) tetrahedra weights
        !
        call read_ibzkpt(&
            nkpts = bz_tmp%nkpts,&
            kpt   = bz_tmp%kpt,&
            w     = bz_tmp%w,&
            ntets = tet%ntets,&
            vtet  = tet%volume,&
            tet   = tet%corner,&
            wtet  = tet%w,&
            iopt_filename = opts%ibzkpt,&
            iopt_verbosity = opts%verbosity )
        !
        if (present(iopt_bz)) then
            !
            if (abs(bz_tmp%nkpts-iopt_bz%nkpts).gt.tiny) then
                call am_print('ERROR','mismatched number of kpoints with IBZKPT',' >>> ')
                stop
            endif
            if (any(abs(bz_tmp%kpt-iopt_bz%kpt).gt.tiny)) then
                call am_print('ERROR','mismatched kpoint coordinates with IBZKPT file',' >>> ')
                stop
            endif
            !
            if (any(abs(bz_tmp%w-iopt_bz%w).gt.tiny)) then
                call am_print('ERROR','mismatched kpoint weight with IBZKPT file',' >>> ')
                call am_print('bz_tmp%w and iopt_bz%w',reshape([bz_tmp%w,iopt_bz%w],[bz_tmp%nkpts,2]))
                stop
            endif
            !
        endif
        !
    end subroutine load_ibzkpt_tet

    subroutine     tetrahedra_pdos(bz,tet,opts)
        !
        use am_histogram
        use am_options
        !
        implicit none
        !
        class(am_class_bz), intent(inout) :: bz
        type(am_class_tetrahedra), intent(in) :: tet
        type(am_class_options), intent(in) :: opts
        real(dp) :: maxE ! maximum band energy
        real(dp) :: minE ! minimum band energy
        real(dp) :: dE   ! energy increment
        real(dp) :: w(4) ! weights on tetrahedron corner
        real(dp) :: Ec(4) ! corner energies
        integer :: i, j, k, m, n, o
        !
        if (opts%verbosity.ge.1) call am_print_title('Tetrahedron integration')
        !
        dE = 0.01_dp
        minE = minval(bz%E(:,:))
        maxE = maxval(bz%E(:,:))
        ! bz%Ep = linspace(minE,maxE,100)
        bz%Ep = regspace(minE,maxE,dE)
        bz%nEp = size(bz%Ep,1)
        !
        call plot_histogram( histogram(pack(bz%E(:,:),(bz%E(:,:).ge.minE).and.(bz%E(:,:).le.maxE)),column_width) )
        !
        if (opts%verbosity.ge.1) call am_print('number of bands',bz%nbands,' ... ')
        if (opts%verbosity.ge.1) call am_print('number of tetrahedra',tet%ntets,' ... ')
        if (opts%verbosity.ge.1) call am_print('lowest band energy',minE,' ... ')
        if (opts%verbosity.ge.1) call am_print('highest band energy',maxE,' ... ')
        if (opts%verbosity.ge.1) call am_print('minimum probing energy',bz%Ep(1),' ... ')
        if (opts%verbosity.ge.1) call am_print('maximum probing energy',bz%Ep(bz%nEp),' ... ')
        if (opts%verbosity.ge.1) call am_print('energy increments dE',dE,' ... ')
        if (opts%verbosity.ge.1) call am_print('number of probing energies',bz%nEp,' ... ')
        !
        allocate(bz%pdos(bz%nspins,bz%norbitals,bz%nions,bz%nEp))
        allocate(bz%dos(bz%nEp))
        !
        bz%pdos = 0.0_dp
        bz%dos = 0.0_dp
        !
        ! two things
        ! 1) difference between total dos and sum of JDOS. ONE band is being left out in the JDOS summation. figure out which one.
        !  0.996425152311378        1.00000000031254
        ! 2) delta summation produces wild values. the divisions by 1/tiny are numerically unstable. Figure out how to get rid of them.!
        !
        do i = 1, bz%nEp
            do j = 1, bz%nbands
            do k = 1, tet%ntets
                !
                Ec(1:4) = bz%E(j,tet%corner(:,k))
                !
                w = tetrahedron_weights(x=bz%Ep(i),xc=Ec,integration_type='delta',apply_blochl_corrections=.true.)
                !
                bz%dos(i) = bz%dos(i) + tet%w(k)*tet%volume(k)*sum(w)
                !
                ! <projections>
                !
                do m = 1,bz%nspins
                do n = 1,bz%norbitals
                do o = 1,bz%nions
                    !> pdos(nspins,norbitals,nions,nEP)
                    !> lmproj(nspins,norbitals,nions,nbands,nkpts)
                    bz%pdos(m,n,o,i) = bz%pdos(m,n,o,i) + tet%w(k)*tet%volume(k)*sum(w*bz%lmproj( m,n,o,j,tet%corner(:,k) ))/real(bz%nspins,dp) ! per spin
                enddo
                enddo
                enddo
                !
                ! </projections>
                !
            enddo
            enddo
            write(*,*) bz%Ep(i), sum(bz%pdos(:,:,:,i)), bz%dos(i)
        enddo
        !
    end subroutine tetrahedra_pdos


!     subroutine get_tetrahedra(tet,bz,uc,pg,tet,iopts)
!         !
!         use am_unit_cell
!         use am_tet_mesh
!         use am_options
!         !
!         implicit none
!         !
!         class(am_class_tetrahedra), intent(inout) :: tet
!         type(am_class_bz), intent(in) :: bz
!         type(am_class_unit_cell) , intent(in) :: uc
!         type(am_class_symmetry), intent(in) :: pg
!         type(am_class_bz) :: fbz
!         type(am_class_options), intent(in), optional :: iopts
!         type(am_class_options) :: opts
!         if (present(iopts)) then 
!            opts = iopts
!         else
!            call opts%defaults
!         endif
!         !
!         call bz%expand_to_fbz(uc=uc,fbz=fbz,pg=pg,iopts=opts)
!         !
!         ! call tessellate_mesh(kpt=fbz%kpt,corner=tet%corner,volume=tet%volume)
!         !
        
!     end subroutine get_tetrahedra

    
    !
    ! things up to here have been edited
    !
    ! INCOMPLETE THINGS: 



!         !
!         ! ibz part:
!         !
!         nR = size(R,3)
!         ! allocate space for ibz
!         allocate(kpoints_ibz(3,bz%nfbz))
!         ! reduce kpoints to ibz
!         do i = 1,bz%nfbz
!            kpoints_ibz(:,i) = reduce_kpoint_to_ibz(bas,recbas,grid_points,R,kpoints(:,i))
!         enddo
!         ! get unique points
!         kpoints_ibz_unique = unique(kpoints_ibz)
!         ! get number of unique points
!         bz%nibz = size(kpoints_ibz_unique,2)
!         ! allocate space for ibz points
!         allocate(bz%k_ibz(bz%nibz))
!         ! save ibz properties
!         do i = 1,bz%nibz
!            ! kpoint coordinate
!            bz%k_ibz(i)%k = kpoints_ibz_unique(:,i)
!            ! index connecting fbz to ibz
!            do j = 1, bz%nfbz
!                if ( all(abs(kpoints_ibz(:,j)-bz%k_ibz(i)%k).lt.tiny) ) bz%k_fbz(j)%irr_indx = i
!            enddo
!            ! orbit
!            korbit = kpoint_orbit(R,bz%k_ibz(i)%k)
!            ! get number of stabilizers (number of elements in little group) and allocate space
!            k = 0
!            do j = 1, nR
!                if ( all(abs(korbit(:,j)-bz%k_ibz(i)%k).lt.tiny) ) k=k+1
!            enddo
!            allocate(bz%k_ibz(i)%littlgrp(3,3,k))
!            ! get little group (stabilizers)
!            k = 0
!            do j = 1, nR
!                if ( all(abs(korbit(:,j)-bz%k_ibz(i)%k).lt.tiny) ) then
!                    k = k + 1
!                    bz%k_ibz(i)%littlgrp(:,:,k) = R(:,:,j)
!                endif
!            enddo
!            ! save orbit
!            bz%k_ibz(i)%korbitit = unique(korbit)
!            ! save weight
!            bz%k_ibz(i)%w = size(bz%k_ibz(i)%korbitit,2)
!         enddo

    !

    !
    !

    !
    !
    !
    !! make tetrahedra
    !subroutine tessellate_fbz(bz)
    !!
    !use am_rank_and_sort
    !use am_tet_mesh
    !!
    !implicit none
    !!
    !class(am_class_bz), intent(inout) :: bz !> brillouin zone class
    !integer , allocatable :: connectivity_list(:,:) ! connectivity_list(1:4,ntets)
    !integer , allocatable :: connectivity_list_ibz(:,:) ! unique connectivity_list for tetrahedra in the ibz.
    !integer , allocatable :: connectivity_list_ibz_unique(:,:) ! unique connectivity_list for tetrahedra in the ibz.
    !integer :: ntets
    !integer :: nitets
    !real(dp), allocatable :: vtet(:) ! tetrahedron volume
    !real(dp) :: kpts(3,bz%nfbz)
    !!
    !integer :: i
    !!
    !integer :: verbosity
    !!
    !verbosity = 1
    !!
    !! make a table of kpoints
    !do i = 1, bz%nfbz
    !    kpts(1:3,i) = bz%k_fbz(i)%k
    !enddo
    !! generate tetrahedra
    !call tessellate_mesh(kpts,connectivity_list,vtet)
    !! save number of tetrahedra
    !ntets = size(connectivity_list,2)
    !! get corner list for ibz points
    !allocate(connectivity_list_ibz(1:4,ntets))
    !do i = 1, ntets
    !    connectivity_list_ibz(1:4,i) = bz%k_fbz( connectivity_list(1:4,i) )%irr_indx
    !    call sort_using_insertion( connectivity_list_ibz(1:4,i) )
    !enddo
    !! get unique corner list for ibz points
    !connectivity_list_ibz_unique = unique(connectivity_list_ibz)
    !nitets = size(connectivity_list_ibz_unique,2)
    !! save irreducible tetrahedra
    !bz%nitets = nitets
    !allocate(bz%tet_ibz(nitets))
    !do i = 1,bz%nitets
    !    bz%tet_ibz(i)%ki = connectivity_list_ibz_unique(:,i)
    !    ! NEED TO MAKE SURE POINTS ON ALL BOUNDARIES ARE INCLUDED IN GENERATION OF FBZ
    !    ! need to figure out a way to determine the weight. compare to connectivity_list_ibz
    !    !
    !enddo
    !
    !
    !
    !end subroutine tessellate_fbz
    !
    !


! subroutine transfer_properties(bz_src,bz_dest)
        ! transfer properties to 
!         bz_dest%nbands    = bz_src%nbands
!         bz_dest%norbitals = bz_src%norbitals
!         bz_dest%nions     = bz_src%nions
!         bz_dest%nspins    = bz_src%nspins
!         bz_dest%nelecs    = bz_src%nelecs
!         allocate( bz_dest%E(bz_dest%nbands,bz_dest%nkpts) )
!         allocate( bz_dest%lmproj(bz_dest%nspins,bz_dest%norbitals,bz_dest%nions,bz_dest%nbands,bz_dest%nkpts) )
!         map_bz_point_properties_to_bz_dest : do i = 1, bz_dest%nkpts
!             do j = 1, bz%nkpts
!                 if ( all(abs(bz%kpt(:,j) - bz_dest%kpt(:,i)).lt.tiny) ) then
!                     bz_dest%E(:,i) = bz_dest%E(:,j)
!                     bz_dest%lmproj(:,:,:,:,i) = bz_dest%lmproj(:,:,:,:,j)
!                     cycle map_bz_point_properties_to_bz_dest
!                 endif
!             enddo
!         enddo map_bz_point_properties_to_ibz

    !
    !
    !

    end module



