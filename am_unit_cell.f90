    module am_unit_cell
    !
    use am_constants
    use am_helpers
    use am_options
    !
    implicit none
    !
    private
    public :: am_class_unit_cell
    
    type am_class_unit_cell
        ! produced by load_poscar:
        real(dp) :: bas(3,3) !> column vectors a(1:3,i), a(1:3,j), a(1:3,k)
        integer :: natoms !> number of atoms
        integer :: nspecies !> number of unique atomic species
        character(len=:), allocatable :: symbs(:) !> symbols of unique atomic species
        integer,allocatable :: natoms_per_species(:) !> number of atoms per species
        real(dp), allocatable :: tau(:,:) !> tau(3,natoms) cartesian atomic coordinates 
        integer,allocatable :: atype(:) !> index indentifying the atomic species for each atom
        ! produced by compute_basic_properties
        real(dp) :: recbas(3,3) !> reciprocal basis vectors,column vectors a(1:3,i), a(1:3,j), a(1:3,k)
        real(dp) :: metric(3,3) !> metric tensor of real space basis
        real(dp) :: vol !> cell volume
    contains
        procedure :: load_poscar
        procedure :: reduce_to_primitive
        procedure :: compute_basic_properties
    end type am_class_unit_cell

    contains

    ! functions which operate on bas(3,3)

    pure function reciprocal_basis(bas) result(recbas)
        !
        implicit none
        !
        real(dp), intent(in) :: bas(3,3)
        real(dp) :: recbas(3,3)
        !
        recbas = 0.0_dp
        recbas(1,1)=bas(2,2)*bas(3,3)-bas(3,2)*bas(2,3)
        recbas(1,2)=bas(3,2)*bas(1,3)-bas(1,2)*bas(3,3)
        recbas(1,3)=bas(1,2)*bas(2,3)-bas(1,3)*bas(2,2)
        recbas(2,1)=bas(2,3)*bas(3,1)-bas(2,1)*bas(3,3)
        recbas(2,2)=bas(1,1)*bas(3,3)-bas(3,1)*bas(1,3)
        recbas(2,3)=bas(2,1)*bas(1,3)-bas(1,1)*bas(2,3)
        recbas(3,1)=bas(2,1)*bas(3,2)-bas(2,2)*bas(3,1)
        recbas(3,2)=bas(3,1)*bas(1,2)-bas(1,1)*bas(3,2)
        recbas(3,3)=bas(1,1)*bas(2,2)-bas(1,2)*bas(2,1)
        recbas = recbas/(recbas(1,1)*bas(1,1)+recbas(1,2)*bas(2,1)+recbas(1,3)*bas(3,1))
    end function  reciprocal_basis

    pure function real_space_metric(bas) result(metric)
        !
        implicit none
        !
        real(dp), intent(in) :: bas(3,3)
        real(dp) :: metric(3,3)
        !
        metric = matmul(transpose(bas),bas)
        !
    end function  real_space_metric



    ! wrappers
    !

    subroutine load_poscar(uc,iopt_filename,iopts)
        ! 
        ! Reads the poscar file, look below: code is short and self explanatory.
        !
        use am_vasp_io
        !
        implicit none
        !
        class(am_class_unit_cell), intent(inout) :: uc
        ! optional i/o
        character(len=*), intent(in), optional :: iopt_filename
        character(1000) :: filename
        type(am_class_options), intent(in), optional :: iopts
        type(am_class_options) :: opts
        if (present(iopts)) then 
            opts = iopts
        else
            call opts%defaults
        endif
        filename = "POSCAR"; if ( present(iopt_filename) ) filename = trim(iopt_filename)
        !
        call read_poscar(&
            bas=uc%bas,&
            natoms=uc%natoms,&
            nspecies=uc%nspecies,&
            natoms_per_species=uc%natoms_per_species,&
            symbs=uc%symbs,&
            tau=uc%tau,&
            atype=uc%atype,&
            iopt_filename=filename,&
            iopt_verbosity=opts%verbosity)
        !
    end subroutine load_poscar

!     subroutine determine_symmetry(ss,uc,iopts)
!         !
!         ! 1) uses atomic positions to identify possible symmetry translations which could serve as primitive lattice vectors
!         ! 2) selects the shortest of these vectors as the basis for the primitive cell
!         ! 3) determines the point symmetries compatbile with the metric tensor corresponding to this primitive basis
!         !       Symmetry and Condensed Matter Physics: A Computational Approach. 1 edition. Cambridge,
!         !       UK; New York: Cambridge University Press, 2008. p 282, table 10.3; p 388, table 10.4
!         ! 4) determines the type of lattice from the number of point symmetries identified in the step above
!         !       Symmetry and Condensed Matter Physics: A Computational Approach. 1 edition.
!         !       Cambridge, UK; New York: Cambridge University Press, 2008. p 282, table 10.3
!         !       p 388, table 10.4
!         ! 5) of the point symmeties found, identify which ones are compatible with the atoms (position and type)
!         ! 6) determine basic point symmetry property properties 
!         ! 7) Sort point symmetries based on det and trace (first point symmetry is identity; last point symmetry is inversion, if present)s
!         !
!         implicit none
!         ! intput/outputs
!         class(am_class_symmetry), intent(inout) :: ss
!         type(am_class_unit_cell), intent(in) :: uc
!         type(am_class_unit_cell) :: prim
!         type(am_class_symmetry) :: ls ! symmetries compatbile with the lattice (not necessarily atomic basis)
!         ! generic point group properties
!         integer, allocatable :: indx(:) ! pairs with mask, see below.
!         logical, allocatable :: mask(:) ! mask which specifies which lattiec-compatible point symmetries are also compatbile with atomic positions
!         character(12) :: crystal_system
!         ! loop variables
!         integer :: i
!         ! options
!         type(am_class_options), intent(in), optional :: iopts
!         type(am_class_options) :: opts
!         if (present(iopts)) then 
!             opts = iopts
!         else
!             call opts%defaults
!         endif
!         !
!         ! write opening lines
!         if ( opts%verbosity .ge. 1 ) write(*,*)
!         if ( opts%verbosity .ge. 1 ) call am_print_title('Analyzing unit cell symmetry')
    
!         !
!         ! check that bas is not empty
!         if ( all( abs(uc%bas).lt.tiny) ) then
!             call am_print('ERROR','Primitive basis is empty! Most likely POSCAR was not loaded.',' >>> ')
!             call am_print('basis (column vectors)',uc%bas)
!             stop
!         endif
!         !
!         ! 1) uses atomic positions to identify possible symmetry translations which could serve as primitive lattice vectors
!         ! 2) selects the shortest of these vectors as the basis for the primitive cell
!         !
!         call uc%reduce_to_primitive(prim=prim,iopts=opts)
!         !
!         ! 3) determines the point symmetries compatbile with this primitive basis
!         !
!         call determine_lattice_point_symmetries(metric=prim%metric,R=ls%R)
!         ls%nptsym = size(ls%R,3)
!         if ( opts%verbosity .ge. 1 ) call am_print("point symmetries compatible with the metric tensor",ls%nptsym," ... ")
!         !
!         ! 4) determines the type of lattice from the number of point symmetries identified in the step above
!         ! 
!         lattice_type : select case (ls%nptsym)
!         case(2);  crystal_system='triclinic   '
!         case(4);  crystal_system='monoclinic  '
!         case(8);  crystal_system='orthorhombic'
!         case(12); crystal_system='trigonal    '
!         case(16); crystal_system='tetragonal  '
!         case(24); crystal_system='hexagonal   '
!         case(48); crystal_system='cubic       '
!         end select lattice_type
!         if ( opts%verbosity .ge. 1 ) call am_print('crystal system',trim(crystal_system)," ... ")
!         !
!         ! 5) Of the point symmeties found, identify which ones are compatible with the atoms (position and type)
!         !
!         allocate(mask(ls%nptsym))
!         allocate(indx(ls%nptsym))
!         mask = .false.
!         do i = 1,ls%nptsym
!             if (is_symmetry_valid(iopt_R=ls%R(:,:,i),tau_frac=uc%tau_cart,iopt_atype=uc%atype,iopt_sym_prec=opts%sym_prec)) then
!                 mask(i) = .true.
!             endif
!             indx(i) = i
!         enddo
!         ss%nptsym = count(mask)
!         if ( opts%verbosity .ge. 1 ) call am_print('point symmetries compatbile with atoms',ss%nptsym," ... ")
!         allocate(ss%R(3,3,ss%nptsym))
!         ss%R = ls%R(:,:,pack(indx,mask))
!         !
!         ! 6) Determine basic point symmetry property properties 
!         !
!         call determine_basic_point_symmetry_properties(R=ss%R,det=ss%det,trace=ss%trace,schoenflies=ss%schoenflies)
!         !
!         ! 7) Sort point symmetries based on det and trace (first point symmetry is identity; last point symmetry is inversion, if present)
!         !
!         call sort_point_symmetries_based_on_determinant_and_trace(R=ss%R,det=ss%det,trace=ss%trace,schoenflies=ss%schoenflies)
!         !
!         if ( .not. all(abs(ss%R(:,:,1)-eye(3)).lt.tiny) ) then
!             call am_print('ERROR','Something is wrong. Identity is not the first point symmetry.')
!             call am_print('R(:,:,1)',ss%R(:,:,i),'     ')
!             stop
!         endif
!         if ( all(abs(ss%R(:,:,ss%nptsym)+eye(3)).lt.tiny) ) then
!             ss%is_centrosymmetric = .true.
!         else
!             ss%is_centrosymmetric = .false.
!         endif
!         !
!         if ( opts%verbosity .ge. 1 ) call am_print('centrosymmetric?',ss%is_centrosymmetric,' ... ')
!         if ( opts%verbosity .ge. 1 ) call ss%write_stdout
!         !
!     end subroutine determine_symmetry
!     subroutine reduce_to_primitive(uc,prim,iopts)
!         ! 
!         ! 1) uses atomic positions to identify possible symmetry translations which could serve as primitive lattice vectors
!         ! 2) selects the shortest of these vectors as the basis for the primitive cell
!         ! 
!         ! essentially a wrapper for determine_primitive_basis, look below: code is short and self explanatory.
!         !
!         class(am_class_unit_cell), intent(in) :: uc
!         type(am_class_unit_cell), intent(out) :: prim
!         type(am_class_options), intent(in), optional :: iopts
!         type(am_class_options) :: opts
!         if (present(iopts)) then 
!             opts = iopts
!         else
!             call opts%defaults
!         endif
!         !
!         call determine_primitive_basis(&
!             bas=uc%bas,&
!             natoms=uc%natoms,&
!             atype=uc%atype,&
!             symbs=uc%symbs,&
!             tau_frac=uc%tau_frac,&
!             oopt_bas_prim=prim%bas,&
!             oopt_natoms_prim=prim%natoms,&
!             oopt_atype_prim=prim%atype,&
!             oopt_tau_prim_cart=prim%tau_cart,&
!             oopt_tau_prim_frac=prim%tau_frac,&
!             natoms_per_species_prim=prim%natoms_per_species,&
!             symbs_prim=prim%symbs,&
!             iopt_verbosity=opts%verbosity)
!         !
!         call prim%compute_basic_properties(iopts=opts)
!         !
!     end subroutine reduce_to_primitive
!     subroutine expand_to_supercell(uc,nTs,sc,iopts)
!         !
!         ! Creates a super cell with dimensions nTs(1),nTs(2),nTs(3) by adding all integer combinations of the vectors  [1...nTs(1),1...nTs(2),1...nTs(3)] to each atomic position. 
!         ! 
!         ! essentially a wrapper for determine_supercell, look below: code is short and self explanatory.
!         !
!         integer, intent(in) :: nTs(3) ! supercell dimensions
!         class(am_class_unit_cell), intent(in) :: uc
!         type(am_class_unit_cell), intent(out) :: sc
!         type(am_class_options), intent(in), optional :: iopts
!         type(am_class_options) :: opts
!         if (present(iopts)) then 
!             opts = iopts
!         else
!             call opts%defaults
!         endif
!         !
!         call determine_supercell(&
!             nTs=nTs,&
!             bas=uc%bas,&
!             natoms=uc%natoms,&
!             atype=uc%atype,&
!             tau_frac=uc%tau_frac,&
!             bas_sc=sc%bas,&
!             natoms_sc=sc%natoms,&
!             atype_sc=sc%atype,&
!             tau_frac_sc=sc%tau_frac,&
!             iopt_verbosity=opts%verbosity)
!         !
!         call sc%compute_basic_properties(iopts=opts)
!         !
!     end subroutine expand_to_supercell
!     subroutine compute_basic_properties(uc,iopts)
!         !
!         ! inputs:
!         ! bas : lattice basis
!         ! tau_frac : atomic coordinates (fractional)
!         ! 
!         ! calculates basic lattice properties
!         ! reciprocal basis, metric tensor, cell volume, atomic coordinates (cartesian)
!         !
!         !
!         class(am_class_unit_cell), intent(inout) :: uc
!         type(am_class_options), intent(in), optional :: iopts
!         type(am_class_options) :: opts
!         if (present(iopts)) then 
!             opts = iopts
!         else
!             call opts%defaults
!         endif
!         !
!         call determine_basic_lattice_properties(&
!             bas=uc%bas,&
!             tau_frac=uc%tau_frac,&
!             recbas=uc%recbas,&
!             metric=uc%metric,&
!             vol=uc%vol,&
!             tau_cart=uc%tau_cart,&
!             iopt_nspecies=uc%nspecies,&
!             iopt_natoms_per_species=uc%natoms_per_species,&
!             iopt_atype=uc%atype,&
!             iopt_symbs=uc%symbs,&
!             iopt_verbosity=opts%verbosity)
!         !
!     end subroutine compute_basic_properties
! 
!     subroutine write_stdout(ss)
!         !
!         implicit none
!         !
!         class(am_class_symmetry), intent(in) ::  ss
!         integer :: i
!         !
!         write(*,'(a5,a)') ' ... ','point symmetry properties'
!         write(*,'(5x,a3,a5,a6,a6)') '#', 'sch.', 'det', 'trace'
!         do i = 1, ss%nptsym
!             write(*,'(5x,i3,a5,f6.2,f6.2)') i, trim(ss%schoenflies(i)), ss%det(i), ss%trace(i)
!         enddo
!         !
!     end subroutine write_stdout

    !
    ! internal
    !


    ! GOAL : WRITE A FUCTION HERE WHICH OUTPUTS R AND T SPACE GROUPS. US IT IN AM_SYMMETRY MODULE.
    ! GOAL : WRITE A FUCTION HERE WHICH OUTPUTS R AND T SPACE GROUPS. US IT IN AM_SYMMETRY MODULE.
    ! GOAL : WRITE A FUCTION HERE WHICH OUTPUTS R AND T SPACE GROUPS. US IT IN AM_SYMMETRY MODULE.

    pure function determine_lattice_point_symmetries(bas) result(R)
        !>
        !> Given a metric tensor, this function returns the (arthimetic) a-holohodry:
        !> the group of point symmetries R (fractional) that are compatible with the
        !> metric tensor).
        !>
        !> Symmetry and Condensed Matter Physics: A Computational Approach. 1 edition.
        !> Cambridge, UK; New York: Cambridge University Press, 2008. pages 278-279,
        !> Eq. 10.60 and 10.61
        !>
        !> The Mathematical Theory of Symmetry in Solids: Representation Theory
        !> for Point Groups and Space Groups. 1 edition. Oxford?: New York:
        !> Oxford University Press, 2010. page 134, table 3.9
        !>
        !> Structure of Materials: An Introduction to Crystallography, Diffraction
        !> and Symmetry. 2 edition. New York: Cambridge University Press, 2012.
        !>
        !> page 169:
        !> "If [point symmetries] are described with respect to the 'primed' basis
        !>  vectors (fractional coordinates), then once again all [point symmetries]
        !>  can be represented by matrices which contain only -1, 0, +1"
        !>
        !> page 166:
        !> "The inverse of a unitary matrix (rotation) is equal to its transpose
        !>  (this is only true in cartesian coordinates.)"
        !>
        !> An example of the procedure:
        !>
        !> 1) Consider the matrix S which is composed of entries which are only -1, 0, and +1:
        !>    S =
        !>    [         1        -1         0 ]
        !>    [         1         0         0 ]
        !>    [         0         0         1 ]
        !>
        !> 2) Furthermore, consider the primitive basis for MoS2 (column vectors):
        !>    A =
        !>    [    3.3269   -1.6635         0 ]
        !>    [         0    2.8812         0 ]
        !>    [         0         0   15.4510 ]
        !>
        !> 3) Evidently, S^T*S is not the identity. That is to say: S is not an orthogonal matrix.
        !>
        !> 4) The basis transformation, C = A*S*B, for which B = inv(A), produces:
        !>    C =
        !>    [    0.5000   -0.8660         0 ]
        !>    [    0.8660    0.5000         0 ]
        !>    [         0         0    1.0000 ]
        !>    which, in fact, is an ORTHOGONAL MATRIX, corresponding to a unitary rotation!
        !>
        !> 5) Note that the order of the operation A*S*B (NOT B*S*A) is important!
        !>    The operation B*S*A produces, instead:
        !>    D =
        !>    [  1.5774  1.1547       0 ]
        !>    [ -1.6547 -0.5774       0 ]
        !>    [       0       0  1.0000 ]
        !>    A matrix which is neither orthogonal, NOR, composed of 0, +1, and -1.
        !>
        !> %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        !> % QUICK MATLAB SCRIPT TO CONFIRM THE ABOVE:
        !>     S=[
        !>          1    -1     0
        !>          1     0     0
        !>          0     0     1];
        !>     A=[
        !>         3.3269   -1.6635         0
        !>              0    2.8812         0
        !>              0         0   15.4510];
        !>     B=inv(A);
        !>     S'*S
        !>     C=A*S*B
        !>     C'*C-eye(3)
        !>
        !> %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        !>
        implicit none
        ! subroutine i/o
        real(dp), intent(in) :: bas(3,3)
        real(dp), allocatable :: R(:,:,:) !> point symmetries (cartesian coordinates)
        !
        real(dp) :: bas(3,3) 
        real(dp) :: recbas(3,3)
        real(dp) :: metric(3,3)
        real(dp) :: buffer(3,3,48) ! buffer for point symmetries
        integer  :: k ! point symmetry counter
        real(dp) :: o(3,3) ! point symmetry in fractional (converted to cartesian before output)
        integer  :: i, i11, i12, i13, i21, i22, i23, i31, i32, i33 ! used to generate unitary rotational matrices
        !
        recbas = reciprocal_basis(bas)
        !
        metric = real_space_metric(bas)
        !
        buffer = 0.0_dp
        !
        k = 0
        !
        do i11 = 1,3
        do i12 = 1,3
        do i13 = 1,3
            do i21 = 1,3
            do i22 = 1,3
            do i23 = 1,3
                do i31 = 1,3
                do i32 = 1,3
                do i33 = 1,3
                    ! can be either [-1,0,1] = [1:3]-2
                    o(1,1:3)=real([i11-2,i12-2,i13-2],dp)
                    o(2,1:3)=real([i21-2,i22-2,i23-2],dp)
                    o(3,1:3)=real([i31-2,i32-2,i33-2],dp)
                    ! check that metric is left unchanged (relaxed tiny bounds here; prone to rounding errors)
                    if ( all( abs(matmul(transpose(o),matmul(metric,o))-metric) .lt. 1E-4 ) ) then
                        k = k + 1
                        ! store point symmetry in cartesian coordinates
                        buffer(1:3,1:3,k) = matmul(bas,matmul(o,recbas))
                    endif
                    !
                enddo
                enddo
                enddo
            enddo
            enddo
            enddo
        enddo
        enddo
        enddo
        !
        allocate(R(3,3,k))
        do i = 1,k
            R(:,:,i)= buffer(:,:,i)
        enddo
        !
    end function  determine_lattice_point_symmetries


    subroutine determine_basic_point_symmetry_properties(R,det,trace,schoenflies)
        !
        implicit none
        !
        real(dp), intent(in) :: R(:,:,:)  !> point symmetries (fractional)
        real(dp), allocatable, intent(out) :: det(:)    !> point symmetry determinants
        real(dp), allocatable, intent(out) :: trace(:)  !> point symmetry traces
        character(4), allocatable, intent(out) :: schoenflies(:)
        integer :: nptsym
        integer :: ierror
        !
        integer :: i
        !
        nptsym = size(R,3)
        !
        allocate(det(nptsym))
        allocate(trace(nptsym))
        allocate(schoenflies(nptsym))
        !
        ! determine fold, trace, and determinant of point symmetries
        det=0.0_dp
        trace=0.0_dp
        do i = 1, nptsym
            ! determine point symmetry trace
            trace(i) = R(1,1,i)+R(2,2,i)+R(3,3,i)
            ! determine point symmetry determinant
            det(i) = -R(1,3,i)*R(2,2,i)*R(3,1,i)+R(1,2,i)*R(2,3,i)*R(3,1,i) &
                   & +R(1,3,i)*R(2,1,i)*R(3,2,i)-R(1,1,i)*R(2,3,i)*R(3,2,i) &
                   & -R(1,2,i)*R(2,1,i)*R(3,3,i)+R(1,1,i)*R(2,2,i)*R(3,3,i)
            ! determine schoenflies notation
            call identify_point_symmetry( trace(i),det(i),schoenflies(i),ierror )
            if ( ierror .ne. 0 ) then
                write(*,'(a5,a)') ' >>> ', 'ERROR: unable to indetify point symmetry'
                call am_print('index',i)
                call am_print('trace',trace(i))
                call am_print('determinant',det(i))
                stop
            endif
        enddo
        !
    end subroutine determine_basic_point_symmetry_properties

    subroutine determine_basic_lattice_properties(bas,tau_frac,recbas,metric,vol,tau_cart,iopt_nspecies,iopt_natoms_per_species,iopt_atype,iopt_symbs,iopt_verbosity)
        !
        implicit none
        ! subroutine i/o
        real(dp), intent(in) :: bas(3,3) ! column vectors a(1:3,i), a(1:3,j), a(1:3,k)
        real(dp), allocatable, intent(in) :: tau_frac(:,:) ! atomic coordinates tau(3,natoms) in fractional
        real(dp), intent(out) :: recbas(3,3)
        real(dp), intent(out) :: metric(3,3)
        real(dp), intent(out) :: vol
        real(dp), allocatable, intent(out) :: tau_cart(:,:) ! atomic coordinates tau(3,natoms) in cartesian
        ! subroutine internal parameters
        integer :: natoms
        ! loop variables
        integer :: i, j, m
        ! optional i/o, just to display on stdout nicely (doesn't make sense to set these if verbosity isn't .ge. 1):
        integer , intent(in), optional :: iopt_nspecies
        integer , intent(in), optional :: iopt_natoms_per_species(:)
        integer , intent(in), optional :: iopt_atype(:)
        character, intent(in), optional :: iopt_symbs(:)
        ! optional i/o
        integer, intent(in), optional :: iopt_verbosity
        integer :: verbosity
        if ( present(iopt_verbosity) ) then
            verbosity = iopt_verbosity
        else
            verbosity = 1
        endif
        !
        if (verbosity .ge. 1) call am_print_title('Determining basic lattice properties')
        ! get cell volume
        vol=abs(bas(1,1)*(bas(2,2)*bas(3,3)-bas(3,2)*bas(2,3)) &
               +bas(2,1)*(bas(3,2)*bas(1,3)-bas(1,2)*bas(3,3)) &
               +bas(3,1)*(bas(1,2)*bas(2,3)-bas(1,3)*bas(2,2)))
        if ( verbosity .ge. 1 ) call am_print( "volume", vol, " ... " )
        ! check that unit cell volume is non-zero
        if ( abs(vol) .lt. tiny ) then
            call am_print("ERROR","Basis vectors are not coplanar!"," >>> ")
            stop
        endif
        if ( verbosity .ge. 1 ) call am_print( "basis (column vectors)", bas, " ... " )
        ! get reciprocal basis as column vectors b
        recbas(1,1)=bas(2,2)*bas(3,3)-bas(3,2)*bas(2,3)
        recbas(1,2)=bas(3,2)*bas(1,3)-bas(1,2)*bas(3,3)
        recbas(1,3)=bas(1,2)*bas(2,3)-bas(1,3)*bas(2,2)
        recbas(2,1)=bas(2,3)*bas(3,1)-bas(2,1)*bas(3,3)
        recbas(2,2)=bas(1,1)*bas(3,3)-bas(3,1)*bas(1,3)
        recbas(2,3)=bas(2,1)*bas(1,3)-bas(1,1)*bas(2,3)
        recbas(3,1)=bas(2,1)*bas(3,2)-bas(2,2)*bas(3,1)
        recbas(3,2)=bas(3,1)*bas(1,2)-bas(1,1)*bas(3,2)
        recbas(3,3)=bas(1,1)*bas(2,2)-bas(1,2)*bas(2,1)
        recbas = recbas/vol
        ! write reciprocal basis
        if ( verbosity .ge. 1 ) call am_print( "rec. basis (column vectors)", recbas, " ... " )
        ! get metrix tensor
        metric = matmul(transpose(bas),bas)
        if ( verbosity .ge. 1 ) call am_print( "metric tensor", metric, " ... " )
        ! convert fractional to cartesian atomic positions
        natoms = size(tau_frac,2)
        allocate(tau_cart(3,natoms))
        do i = 1, natoms
            tau_cart(1:3,i) = matmul( bas ,tau_frac(1:3,i) )
        enddo
        !
        if ( verbosity .ge. 1 ) then
            if ( present(iopt_atype) .and. present(iopt_symbs) .and. present(iopt_natoms_per_species) .and. present(iopt_nspecies) ) then
                !
                write(*,'(" ... ",a)') "atomic positions (fractional)"
                m = 0
                do i = 1, iopt_nspecies
                    do j = 1, iopt_natoms_per_species(i)
                        m = m + 1
                        write(*,'(5x,a5,a3,3f13.8)') int2char(m), iopt_symbs(iopt_atype(m)), tau_frac(1:3,m)
                    enddo
                enddo
                !
                write(*,'(" ... ",a)') "atomic positions (cartesian)"
                m = 0
                do i = 1, iopt_nspecies
                    do j = 1, iopt_natoms_per_species(i)
                        m = m + 1
                        write(*,'(5x,a5,a3,3f13.8)') int2char(m), iopt_symbs(iopt_atype(m)), tau_cart(1:3,m)
                    enddo
                enddo
            else
                write(*,'(" ... ",a)') "atomic positions (fractional)"
                m = 0
                do i = 1, natoms
                    write(*,'(5x,a5,3f13.8)') int2char(i), tau_frac(1:3,i)
                enddo
                write(*,'(" ... ",a)') "atomic positions (cartesian)"
                m = 0
                do i = 1, natoms
                    write(*,'(5x,a5,3f13.8)') int2char(i), tau_cart(1:3,i)
                enddo
            endif
        endif
        !
    end subroutine determine_basic_lattice_properties

    subroutine determine_supercell(nTs,bas,natoms,atype,tau_frac,bas_sc,natoms_sc,atype_sc,tau_frac_sc,iopt_verbosity)
        !
        ! expands fractional atomic coordinates onto an nTs(1) x nTs(2) x nTs(3) supercell
        ! also maps the type of atom from atype (in the primitive cell) to atype_sc (in the supercell).
        !
        implicit none
        !
        integer,intent(in) :: nTs(3) ! supercell dimensions
        real(dp), intent(in) :: bas(3,3) !> primitive basis vectors
        integer,intent(in) :: atype(:) !> list identify type of atom
        real(dp), intent(in) :: tau_frac(:,:) !> atomic coordinates tau_frac(1:3,natoms)
        real(dp), intent(out) :: bas_sc(3,3) !> supercell basis vectors
        real(dp), allocatable, intent(out) :: tau_frac_sc(:,:) !> atomic coordinates tau_frac_sc(1:3,natoms) in the supercell
        integer,allocatable, intent(out) :: atype_sc(:) !> list identify type of atom in the supercell
        integer, intent(in) :: natoms
        integer, intent(out) :: natoms_sc
        integer :: i1,i2,i3
        integer :: i, j, m
        !
        integer, intent(in), optional :: iopt_verbosity
        integer :: verbosity
        if ( present(iopt_verbosity) ) then
            verbosity = iopt_verbosity
        else
            verbosity = 1
        endif
        !
        !
        if (verbosity .ge. 1) write(*,*)
        if (verbosity .ge. 1) write(*,'(a)') 'Constructing supercell'
        !
        if ( size(tau_frac,2) .ne. natoms ) then
            call am_print('ERROR','size(tau_frac,2) .ne. natoms',' >>> ')
            stop
        endif
        !
        if (verbosity .ge. 1) call am_print('dimensions',nTs,' ... ')
        !
        bas_sc = matmul(bas,1.0_dp*reshape([nTs(1),0,0,0,nTs(2),0,0,0,nTs(3)],[3,3]))
        !
        if (verbosity .ge. 1)call am_print('supercell basis',bas_sc,' ... ')
        !
        natoms_sc = nTs(1)*nTs(2)*nTs(3)*natoms
        allocate(tau_frac_sc(3,natoms_sc))
        allocate(atype_sc(natoms_sc))
        !
        if (verbosity .ge. 1)call am_print('number of atoms in supercell',natoms_sc,' ... ')
        !
        tau_frac_sc = 0.0_dp
        m=0
        do i1 = 1, nTs(1)
            do i2 = 1, nTs(2)
                do i3 = 1, nTs(3)
                    do j = 1,natoms
                        m=m+1
                        tau_frac_sc(1:3,m) = ( tau_frac(1:3,j) + ([i1,i2,i3]-1)*1.0_dp ) / ( nTs(1:3) *1.0_dp )
                        atype_sc(m) = atype(j)
                    enddo
                enddo
            enddo
        enddo
        !
        if (verbosity .ge. 1) then
            write(*,'(" ... ",a)') "supercell atomic positions (fractional)"
            do i = 1,natoms_sc
                write(*,'(5x,a5,3f13.8)') int2char(i), tau_frac_sc(1:3,i)
            enddo
            write(*,'(" ... ",a)') "supercell atomic positions (cartesian)"
            do i = 1,natoms_sc
                write(*,'(5x,a5,3f13.8)') int2char(i), matmul(bas_sc,tau_frac_sc(1:3,i))
            enddo
        endif
    end subroutine determine_supercell

    pure function determine_translations_from_atomic_basis(atype,tau_frac,iopt_sym_prec,iopt_verbosity) result(T_frac)
        !
        implicit none
        ! subroutine i/o
        integer , intent(in)  :: atype(:) !> list identify type of atom
        real(dp), intent(in)  :: tau_frac(:,:) ! atomic coordinates tau_frac(1:3,natoms)
        real(dp), allocatable :: T_frac(:,:)
        ! internal
        integer  :: natoms ! number of atoms
        integer  :: nTs ! number of primitive lattice vector candidates
        real(dp) :: tau_ref_frac(3) ! fractional coordinates of reference atom
        real(dp), allocatable :: wrk(:,:) ! wrk(1:3,nTs) list of possible lattice vectors that are symmetry-compatible volume
        real(dp) :: wrkr(3)
        ! loop variables
        integer :: i, j
        ! optionals
        real(dp), intent(in), optional :: iopt_sym_prec
        real(dp) :: sym_prec
        integer, intent(in), optional :: iopt_verbosity
        integer :: verbosity
        if ( present(iopt_verbosity) ) then
            verbosity = iopt_verbosity
        else
            verbosity = 1
        endif
        !
        if ( present(iopt_sym_prec) ) then
            sym_prec = iopt_sym_prec
        else
            sym_prec = tiny
        endif
        !
        ! get number of atoms
        natoms = size(tau_frac,2)
        ! setup workspace
        allocate(wrk(3,natoms-1+3))
        wrk = 0.0_dp
        ! choose one atom (ideally choose an atom corresponding to the species with the fewest number of atoms in unit cell)
        i = 1
        tau_ref_frac(1:3) =  tau_frac(1:3,i)
        ! search for lattice vectors using, as translational components, vectors connecting the choosen atom to other atoms of the same species
        nTs = 0
        do j = 1,natoms
            if ( i .ne. j) then
            if ( atype(i) .eq. atype(j) ) then
                ! shift to put reference atom at zero.
                wrkr(1:3) = tau_frac(1:3,j) - tau_ref_frac
                if ( is_symmetry_valid(iopt_T=wrkr,iopt_atype=atype, tau_frac=tau_frac, iopt_sym_prec=sym_prec) ) then
                    ! lattice vector found
                    nTs = nTs + 1
                    wrk(1:3,nTs) = wrkr
                endif
            endif
            endif
        enddo
        ! add the original lattice vectors
        nTs = nTs+1
        wrk(1:3,nTs) = real([1,0,0],dp)
        nTs = nTs+1
        wrk(1:3,nTs) = real([0,1,0],dp)
        nTs = nTs+1
        wrk(1:3,nTs) = real([0,0,1],dp)
        ! allocate output
        allocate(T_frac(3,nTs))
        T_frac(1:3,1:nTs) = wrk(1:3,1:nTs)
    end function determine_translations_from_atomic_basis






    subroutine determine_primitive_basis(bas,natoms,atype,tau_frac, &
        symbs,&
        oopt_bas_prim, &
        oopt_natoms_prim, &
        oopt_tau_prim_frac, &
        oopt_tau_prim_cart, &
        oopt_atype_prim, &
        symbs_prim,&
        natoms_per_species_prim, &
        iopt_sym_prec, &
        iopt_verbosity)
        !
        !The program first finds the primitive unit cell by looking for
        !additional lattice vectors within the unit cell given in the input.
        !
        ! It chooses one of the atoms and tests each vector from that atom to
        !every other atom of the same type in the unit cell. If that vector
        !will take us from each atom in the unit cell to another atom of the
        !same type (not necessarily in the same unit cell), then we have
        !found a new lattice vector.
        !
        ! This process is repeated until we have
        !found a complete list of lattice vectors. From that list, we find three
        !which can serve as primitive basis vectors ti of the lattice, i.e. every
        !lattice vector can be expressed as an integer linear combination of
        !these basis vectors.
        !
        use am_rank_and_sort
        !
        implicit none
        ! subroutine i/o
        real(dp), intent(in) :: bas(3,3) !> basis
        integer , intent(in) :: natoms
        integer , intent(in) :: atype(:) !> list identify type of atom
        real(dp), intent(in) :: tau_frac(:,:) ! atomic coordinates tau_frac(1:3,natoms)
        character(*), intent(in) :: symbs(*)
        character(:), intent(out), allocatable :: symbs_prim(:)
        real(dp), intent(out), optional :: oopt_bas_prim(3,3) !> basis of primitive cell
        integer , intent(out), optional :: oopt_natoms_prim
        real(dp), intent(out), optional, allocatable :: oopt_tau_prim_frac(:,:)
        real(dp), intent(out), optional, allocatable :: oopt_tau_prim_cart(:,:)
        integer , intent(out), optional, allocatable :: oopt_atype_prim(:)
        integer , allocatable, intent(out) :: natoms_per_species_prim(:)
        ! internal
        real(dp) :: bas_prim(3,3) !> basis of primitive cell
        integer  :: natoms_prim
        real(dp) :: tau_ref_frac(3) ! fractional coordinates of reference atom
        real(dp), allocatable :: tau_cart(:,:)
        real(dp), allocatable :: tau_prim_frac(:,:)
        real(dp), allocatable :: tau_prim_cart(:,:)
        integer , allocatable :: atype_prim(:)
        real(dp), allocatable :: T_frac(:,:) ! rlist(1:3,nTs) list of possible lattice vectors that are symmetry-compatible volume
        real(dp), allocatable :: T_norm(:) ! magnitude of possible lattice vectors
        integer , allocatable :: indices(:) ! used to sort based on magnitude
        real(dp) :: bpscf(3,3) ! primitive basis (fractional, in terms of the supercell basis)
        real(dp) :: bscpf(3,3) ! supercell basis (fractional, in terms of the primitive basis)
        real(dp) :: wrkbas(3,3)! workspace for determining lattice vectors corresponding to the smallest primitive cell volume
        real(dp) :: wrkvol     ! workspace for determining lattice vectors corresponding to the smallest primitive cell volume
        integer :: nTs ! number of primitive lattice vector candidates
        real(dp) :: vol ! volume of primitive cell
        integer :: i, j, k ! loop variables
        ! optionals
        real(dp), intent(in), optional :: iopt_sym_prec
        real(dp) :: sym_prec
        integer, intent(in), optional :: iopt_verbosity
        integer :: verbosity
        if ( present(iopt_verbosity) ) then
            verbosity = iopt_verbosity
        else
            verbosity = 1
        endif
        if ( present(iopt_sym_prec) ) then
            sym_prec = iopt_sym_prec
        else
            sym_prec = tiny
        endif
        !
        if ( natoms .ne. size(tau_frac,2) ) then
            call am_print('ERROR','natoms .ne. size(tau_frac,2)',' >>> ')
            stop
        endif
        !
        !
        !
        T_frac = determine_translations_from_atomic_basis(atype,tau_frac,iopt_sym_prec,iopt_verbosity)
        nTs = size(T_frac,2)
        ! get norm of each translation
        allocate(T_norm(nTs))
        allocate(indices(nTs))
        T_norm(1:nTs) = norm2(T_frac,1)
        ! sort based on magnitude (smallest last)
        call rank(T_norm,indices)
        indices = indices(nTs:1:-1)
        T_frac = T_frac(:,indices)
        ! writ to stdout
        if (verbosity .ge. 1) call am_print("possible primitive lattice translations found",nTs," ... ")
        if (verbosity .ge. 1) call am_print("possible primitive lattice vectors (fractional; column vectors, tranposed for clarity)",transpose(T_frac)," ... ")
        if (verbosity .ge. 1) call am_print("possible primitive lattice vectors (cartesian; column vectors, tranposed for clarity)",transpose(matmul(bas,T_frac))," ... ")
        !
        !
        ! select as primitive basis the lattice vectors with the smallest vector magnitude (i.e. the last three vectors from the sort above)
        !
        ! Notice that this procedure, for Si, rather than giving the standard fcc basis:
        !	0.50000        0.50000        0.00000
        !	0.50000        0.00000        0.50000
        !	0.00000        0.50000        0.50000
        ! produces this instead, which is the vectors with the smallest magnitude and the smallest cell volume.
        !	0.00000        0.00000        0.50000
        !	0.00000        0.50000        0.00000
        !   0.50000        0.00000        0.00000
        !
        if ( present(oopt_bas_prim) .or. present(oopt_tau_prim_frac) .or. present(oopt_natoms_prim) .or. present(oopt_atype_prim) ) then
            vol = 10000.0_dp
            ! set first lattice vector as the last vector candidate
            i = nTs
            ! loop until a second non-collinear lattice vector is found (non-zero cross product)
            do j = nTs,1,-1
                if ( any(cross_product(T_frac(1:3,i),T_frac(1:3,j)) .gt. tiny) ) exit
            enddo
            if (any([i,j].eq.0)) then
                call am_print('ERROR','No pair of non-collinear lattice vectors found',' >>> ')
                stop
            endif
            ! loop until third lattice vector is found which yields a non-zero cell volume
            do k = nTs,1,-1
                wrkbas = T_frac(1:3,[i,j,k])
                wrkvol = abs(wrkbas(1,1)*(wrkbas(2,2)*wrkbas(3,3)-wrkbas(3,2)*wrkbas(2,3)) &
                            +wrkbas(2,1)*(wrkbas(3,2)*wrkbas(1,3)-wrkbas(1,2)*wrkbas(3,3)) &
                            +wrkbas(3,1)*(wrkbas(1,2)*wrkbas(2,3)-wrkbas(1,3)*wrkbas(2,2)))
                if ( wrkvol .gt. tiny ) exit
            enddo
            if (any([i,j,k].eq.0)) then
                call am_print('ERROR','No primitive basis found',' >>> ')
                stop
            endif
            ! remember bpscf, primitive basis (fractional, in terms of the supercell basis)
            bpscf = T_frac(1:3,[i,j,k])
            !
            if (verbosity .ge. 1) call am_print("primitive basis (fractional, in terms of supercell basis)",bpscf," ... ")
            ! get the inverse of bpscf, corresponding to the supercell basis (fractional, in terms of the primitive basis)
            bscpf(1,1)=bpscf(2,2)*bpscf(3,3)-bpscf(3,2)*bpscf(2,3)
            bscpf(1,2)=bpscf(3,2)*bpscf(1,3)-bpscf(1,2)*bpscf(3,3)
            bscpf(1,3)=bpscf(1,2)*bpscf(2,3)-bpscf(1,3)*bpscf(2,2)
            bscpf(2,1)=bpscf(2,3)*bpscf(3,1)-bpscf(2,1)*bpscf(3,3)
            bscpf(2,2)=bpscf(1,1)*bpscf(3,3)-bpscf(3,1)*bpscf(1,3)
            bscpf(2,3)=bpscf(2,1)*bpscf(1,3)-bpscf(1,1)*bpscf(2,3)
            bscpf(3,1)=bpscf(2,1)*bpscf(3,2)-bpscf(2,2)*bpscf(3,1)
            bscpf(3,2)=bpscf(3,1)*bpscf(1,2)-bpscf(1,1)*bpscf(3,2)
            bscpf(3,3)=bpscf(1,1)*bpscf(2,2)-bpscf(1,2)*bpscf(2,1)
            bscpf = bscpf/(bpscf(1,1)*(bpscf(2,2)*bpscf(3,3)-bpscf(3,2)*bpscf(2,3)) &
                          +bpscf(2,1)*(bpscf(3,2)*bpscf(1,3)-bpscf(1,2)*bpscf(3,3)) &
                          +bpscf(3,1)*(bpscf(1,2)*bpscf(2,3)-bpscf(1,3)*bpscf(2,2)))
            if (verbosity .ge. 1) call am_print("supercell basis (fractional, in terms of primitive basis)",bscpf," ... ")
            ! compute the primitive basis for the primitive cell
            bas_prim = matmul(bas,bpscf)
            if (verbosity .ge. 1) call am_print("primitive basis (cartesian)",bas_prim," ... ")
            ! get primitive cell volume 
            vol = abs(bas_prim(1,1)*(bas_prim(2,2)*bas_prim(3,3)-bas_prim(3,2)*bas_prim(2,3)) &
                     +bas_prim(2,1)*(bas_prim(3,2)*bas_prim(1,3)-bas_prim(1,2)*bas_prim(3,3)) &
                     +bas_prim(3,1)*(bas_prim(1,2)*bas_prim(2,3)-bas_prim(1,3)*bas_prim(2,2)))
            if (verbosity .ge. 1) call am_print("primitive cell volume",vol," ... ")
            ! write to stdout
            if (present(oopt_bas_prim)) oopt_bas_prim = bas_prim
        endif
        ! 
        ! get atoms which are within the boundaries of bas_prim
        !
        if (present(oopt_tau_prim_frac) .or. present(oopt_natoms_prim) .or. present(oopt_tau_prim_cart) .or. present(oopt_atype_prim) ) then
            allocate(tau_prim_frac(3,natoms))
            do i = 1,natoms
                ! shift all atomic position based on reference atom tau_ref_frac
                tau_prim_frac(1:3,i) = tau_frac(1:3,i) - tau_ref_frac
                ! convert to primitive fractional coordinates
                tau_prim_frac(1:3,i) = matmul(bscpf,tau_prim_frac(1:3,i))
            enddo
            if (verbosity .ge. 1) call am_print('number of atoms in supercell',natoms,' ... ')
            if (verbosity .ge. 1) call am_print('atomic coordinates in supercell (fractional, in terms of primitive basis)',transpose(tau_prim_frac),' ... ')
            ! reduce to primitive cell
            do i = 1,natoms
                tau_prim_frac(1:3,i) = modulo(tau_prim_frac(1:3,i)+sym_prec,1.0_dp)-sym_prec
            enddo
            tau_prim_frac=unique(tau_prim_frac,sym_prec)
            natoms_prim = size(tau_prim_frac,2)
            !
            if (verbosity .ge. 1) call am_print('number of atoms in primitive cell',natoms_prim,' ... ')
            if (present(oopt_natoms_prim))   oopt_natoms_prim   = natoms_prim
            if (present(oopt_tau_prim_frac)) then
                oopt_tau_prim_frac = tau_prim_frac
                if (verbosity .ge. 1) call am_print('positions of atoms in primitive cell (fractional)',transpose(oopt_tau_prim_frac),' ... ')
            endif
            ! get atomic primitive atomic coordinates in cartesian
            tau_prim_cart = matmul(bas_prim,tau_prim_frac)
            if (present(oopt_tau_prim_cart)) then
                allocate(oopt_tau_prim_cart(3,natoms_prim))
                oopt_tau_prim_cart = tau_prim_cart
                if (verbosity .ge. 1) call am_print('positions of atoms in primitive cell (cartesian)',transpose(oopt_tau_prim_cart),' ... ')
            endif
            ! determine which atoms are of what type by comparing cartesian coordinates
            if (present(oopt_atype_prim)) then
                !
                allocate(tau_cart(3,natoms))
                tau_cart = matmul(bas,tau_frac)
                !
                allocate(atype_prim(natoms_prim))
                allocate(character(500) :: symbs_prim(natoms_prim))
                do i = 1,natoms_prim
                    atype_prim(i) = 0
                    search_for_similar_atom : do j = 1,natoms
                        if ( all(abs(tau_prim_cart(1:3,i)-tau_cart(1:3,j)).lt.sym_prec) ) then
                            atype_prim(i) = atype(j)
                            symbs_prim(i)  = symbs(j)
                            exit search_for_similar_atom
                        endif
                    enddo search_for_similar_atom
                    if (atype_prim(i) .eq. 0) then
                        call am_print('ERROR','Unable to identify type of atom in the primitive cell',' >>> ')
                        call am_print('atomic coordinates in original cell (cartesian)',tau_cart)
                        call am_print('atomic coordinates in reduced cell (cartesian)',tau_prim_cart)
                        stop
                    endif
                enddo
                allocate(oopt_atype_prim(natoms_prim))
                oopt_atype_prim = atype_prim
                !
                allocate(natoms_per_species_prim(size(unique(atype_prim))))
                do i = 1,maxval(atype_prim)
                    natoms_per_species_prim(i) = count(atype_prim.eq.i)
                enddo
            endif
        endif
    end subroutine determine_primitive_basis
    pure function is_symmetry_valid(iopt_R,iopt_T,iopt_atype,tau_frac,iopt_sym_prec)
        !
        ! check whether symmetry operation is valid
        !
        implicit none
        ! function i/o
        integer , intent(in), optional :: iopt_atype(:) !> list identify type of atom
        real(dp), intent(in) :: tau_frac(:,:) ! atomic coordinates tau_frac(1:3,natoms)
        real(dp), intent(in), optional :: iopt_R(3,3)
        real(dp), intent(in), optional :: iopt_T(3)
        real(dp), intent(in), optional :: iopt_sym_prec
        ! required for optional i/o
        integer , allocatable :: atype(:)
        real(dp) :: R(3,3)
        real(dp) :: T(3)
        real(dp) :: sym_prec
        ! internal
        real(dp) :: tau_rot_frac(3) ! rotated 
        logical :: point_symmetry
        logical :: translational_symmetry
        logical :: is_symmetry_valid
        integer :: natoms
        integer :: i, j, m
        logical :: overlap_found
        !
        if ( present(iopt_sym_prec) ) then
            sym_prec = iopt_sym_prec
        else
            sym_prec = tiny
        endif
        !
        if ( present(iopt_R) ) then
            point_symmetry = .true.
            R = iopt_R
        else
            point_symmetry = .false.
            R = reshape([1,0,0,0,1,0,0,0,1],[3,3])
        endif
        !
        !
        if ( present(iopt_T) ) then
            translational_symmetry = .true.
            T = iopt_T
        else
            translational_symmetry = .false.
            T = [0,0,0]*0.0_dp
        endif
        !
        if ( present(iopt_atype) ) then
            allocate(atype(size(iopt_atype)))
            atype = iopt_atype
        else
            allocate(atype(size(tau_frac)))
            atype = 1
        endif
        !
        natoms = size(tau_frac,2)
        !
        m = 0
        do i = 1, natoms
            ! apply symmetry operation
            tau_rot_frac(1:3) = matmul(R,tau_frac(1:3,i)) + T(1:3)
            ! reduce rotated+translated point to reciprocal unit cell
            tau_rot_frac(1:3) = modulo(tau_rot_frac(1:3)+sym_prec,1.0_dp)-sym_prec
            ! check that newly created point matches something already present
            overlap_found = .false.
            check_overlap : do j = 1,natoms
                if (atype(i) .eq. atype(j)) then
                    if ( all(abs(tau_rot_frac(1:3)-tau_frac(1:3,j)).lt.sym_prec) ) then
                        m = m + 1
                        overlap_found = .true.
                        exit check_overlap
                    endif
                endif
            enddo check_overlap
            if ( .not. overlap_found ) then
                is_symmetry_valid = .false.
                return
            endif
        enddo
        if (m .eq. natoms) then
            is_symmetry_valid = .true.
        endif
    end function is_symmetry_valid

    subroutine apply_elastic_deformation(bas,deformation_code,strain_max,nstrains,bas_def,iopt_verbosity)
        !
        ! creates at bas_def(3,3,nTs) matrix containing the nTs deformed primitive basises.
        !
        ! call apply_elastic_deformation(bas=uc%bas,deformation_code=6,strain_max=0.1_dp,nstrains=10,bas_def=bas_def,iopt_verbosity=verbosity)
        !
        implicit none
        !
        real(dp), intent(in) :: bas(3,3)
        integer , intent(in) :: deformation_code
        real(dp), intent(in) :: strain_max
        integer , intent(in) :: nstrains
        real(dp), intent(out), allocatable :: bas_def(:,:,:)
        !
        character(len=6) :: dc ! deformation code string
        integer :: i, j ! loop variable
        real(dp) :: eta_step ! the deformation step size
        real(dp) :: eta ! the magnitude of the deformation
        real(dp) :: e(6) ! strain voigt
        real(dp) :: eta_matrix(3,3) ! deformation matrix (lagrangian)
        real(dp) :: eps_matrix(3,3) ! deformation matrix (true)
        real(dp), allocatable :: def_mat(:,:,:) ! deformation matrix
        !
        integer, intent(in), optional :: iopt_verbosity
        integer :: verbosity
        if ( present(iopt_verbosity) ) then
            verbosity = iopt_verbosity
        else
            verbosity = 1
        endif
        !
        if (verbosity.ge.1) write(*,*)
        if (verbosity.ge.1) write(*,*) 'Elastic deformation'
        !
        if (nstrains .le. 1) then
            call am_print('ERROR','More than one strain point is required.',' >>> ')
            stop
        endif
        !
        if (strain_max .gt. 0.1) then
            call am_print('WARNI','Are you sure you want to use such large strain values? Continuing...',' >>> ')
        endif
        !
        dc = 'XXXXXX'
        select case(deformation_code)
            case(0);  dc='EEE000'
            case(1);  dc='E00000'
            case(2);  dc='0E0000'
            case(3);  dc='00E000'
            case(4);  dc='000E00'
            case(5);  dc='0000E0'
            case(6);  dc='00000E'
            case(7);  dc='000EEE'
            case(8);  dc='EE0000'
            case(9);  dc='Ee0000'
            case(10); dc='EEEEEE'
            case(11); dc='123456'
            case(12); dc='d14t6c'
            case(13); dc='3cu62q'
            case(14); dc='qs51t2'
            case(15); dc='546dut'
            case(16); dc='s3d541'
            case(20); dc='000EE0'
            case(40); dc='EEd000'
            case(41); dc='Ee0000'
        end  select
        !
        if ( dc .eq. 'XXXXXX' ) then
            call am_print('ERROR','Deformation code not found. Select one of the following:',' >>> ')
            write(*,'(5x,a)') "  0 =>  ( eta,  eta,  eta,    0,    0,    0)  | volume strain "
            write(*,'(5x,a)') "  1 =>  ( eta,    0,    0,    0,    0,    0)  | linear strain along x "
            write(*,'(5x,a)') "  2 =>  (   0,  eta,    0,    0,    0,    0)  | linear strain along y "
            write(*,'(5x,a)') "  3 =>  (   0,    0,  eta,    0,    0,    0)  | linear strain along z "
            write(*,'(5x,a)') "  4 =>  (   0,    0,    0,  eta,    0,    0)  | yz shear strain"
            write(*,'(5x,a)') "  5 =>  (   0,    0,    0,    0,  eta,    0)  | xz shear strain"
            write(*,'(5x,a)') "  6 =>  (   0,    0,    0,    0,    0,  eta)  | xy shear strain 2C44"
            write(*,'(5x,a)') "  7 =>  (   0,    0,    0,  eta,  eta,  eta)  | shear strain along (111)"
            write(*,'(5x,a)') "  8 =>  ( eta,  eta,    0,    0,    0,    0)  | xy in-plane strain "
            write(*,'(5x,a)') "  9 =>  ( eta, -eta,    0,    0,    0,    0)  | xy in-plane shear strain"
            write(*,'(5x,a)') " 10 =>  ( eta,  eta,  eta,  eta,  eta,  eta)  | global strain"
            write(*,'(5x,a)') " 11 =>  (  1e,   2e,   3e,   4e,   5e,   6e)  | universal strain 1"
            write(*,'(5x,a)') " 12 =>  ( -2e,   1e,   4e,  -3e,   6e,  -5e)  | universal strain 2"
            write(*,'(5x,a)') " 13 =>  (  3e,  -5e,  -1e,   6e,   2e,  -4e)  | universal strain 3"
            write(*,'(5x,a)') " 14 =>  ( -4e,  -6e,   5e,   1e,  -3e,   2e)  | universal strain 4"
            write(*,'(5x,a)') " 15 =>  (  5e,   4e,   6e,  -2e,  -1e,  -3e)  | universal strain 5"
            write(*,'(5x,a)') " 16 =>  ( -6e,   3e,  -2e,   5e,  -4e,   1e)  | universal strain 6"
            write(*,'(5x,a)') " 20 =>  (   0,    0,    0,  eta,  eta,    0)  | yz zx shear"
            write(*,'(5x,a)') " 40 =>  ( eta,  eta,  -2e,    0,    0,    0)  | tetragonal for 3/2(c11-c12)"
            write(*,'(5x,a)') " 41 =>  ( eta, -eta,    0,    0,    0,    0)  | orthorhombic for (c11-c12)"
            stop
        endif
        !
        eta_step = 2.0_dp*strain_max/(nstrains-1.0_dp)
        !
        allocate(bas_def(3,3,nstrains))
        allocate(def_mat(3,3,nstrains))
        !
        if (verbosity.ge.1) call am_print('Number of deformations',nstrains,' ... ')
        !
        do i = 1,nstrains
            ! set deformation magnitude
            eta = (i-1)*eta_step - strain_max
            if (abs(eta) < 1.0E-6_dp) eta=1.0E-6_dp ! why is this here?
            ! convert string to deformation
            e(1:6)=0.0_dp
            do j = 1,6
                select case (dc(j:j))
                    case('E'); e(j)= eta
                    case('e'); e(j)=-eta
                    case('0'); e(j)=0.0_dp
                    case('1'); e(j)= 1.0_dp*eta
                    case('2'); e(j)= 2.0_dp*eta
                    case('3'); e(j)= 3.0_dp*eta
                    case('4'); e(j)= 4.0_dp*eta
                    case('5'); e(j)= 5.0_dp*eta
                    case('6'); e(j)= 6.0_dp*eta
                    case('u'); e(j)=-1.0_dp*eta
                    case('d'); e(j)=-2.0_dp*eta
                    case('t'); e(j)=-3.0_dp*eta
                    case('q'); e(j)=-4.0_dp*eta
                    case('c'); e(j)=-5.0_dp*eta
                    case('s'); e(j)=-6.0_dp*eta
                    case default
                        call am_print('ERROR','Deformation not allowed',' >>> ')
                        stop
                end select
            enddo
            ! build deformation matrix
            eta_matrix(1:3,1) = [e(1)       ,e(6)/2.0_dp,e(5)/2.0_dp]
            eta_matrix(1:3,2) = [e(6)/2.0_dp,e(2)       ,e(4)/2.0_dp]
            eta_matrix(1:3,3) = [e(5)/2.0_dp,e(4)/2.0_dp,e(3)       ]
            !
            ! Solve for eps (true strain) from eta (lagrangian strain) by iteration
            !
            ! "ElaStic: A tool for calculating second-order elastic constants from first principles"
            ! Eq. 2 : tau = det(1+eps)(1+eps)^-1 . sigma . (1+eps)^-1 ; dot indicates tensor product
            ! Eq. 1 : eta = eps + 1/2 * eps*eps
            ! Also see Wallace "Thermodynamics of crystals"
            !
            ! % Matlab script
            ! % maximum difference is 1E-17 for e < 0.1, which is more than enough strain for elastic constants
            ! % if more strain is desired, then a differen method for solving for eps is required.
            ! for j = 1:200
            ! e = rand(6,1)*0.1;
            ! eta = [e(1),e(6)/2,e(5)/2;
            !        e(6)/2,e(2),e(4)/2;
            !        e(5)/2,e(4)/2,e(3)];
            ! eps = eta;
            ! for i = 1:300
            !     % eta = eps + 1/2eps*eps => rearrange:  eps = eta - 1/2eps*eps
            !     eps = eta - 1/2*eps*eps;
            ! end
            ! hold on; plot(j,max((eps+ 1/2*eps*eps-eta)),'s')
            ! end
            !
            eps_matrix = eta_matrix
            solve_for_eps : do j = 1,200
                eps_matrix = eta_matrix - 0.5_dp*matmul(eps_matrix,eps_matrix)
                if ( all(abs(eps_matrix).gt.10.0_dp) ) then
                    ! if this is the case then a differen method for solving for eps is required than iteration.
                    call am_print('ERROR','eps too large; possibly because it is diverging with each iteration: deformation may be too large.',' >>> ')
                    call am_print('eps (true strain matrix)',eps_matrix)
                    stop
                endif
                if ( all(abs(eta-eps_matrix).lt.1E-15_dp) ) exit solve_for_eps
            enddo solve_for_eps
            !
            def_mat(1:3,1:3,i) = eye(3) + eps_matrix
            bas_def(1:3,1:3,i)=matmul(def_mat(1:3,1:3,i),bas)
            !
        enddo
        !
        if (verbosity.ge.1) then
            write(*,'(a5,a)') ' ... ', 'Deformation matrices (true strain, Voigt notation)'
            write(*,'(a5,2x,6a16)') ' ... ', 'e1', 'e2', 'e3', 'e4', 'e5', 'e6'
            do i = 1,nstrains
                write(*,'(5x,i2,6f16.10)') i,def_mat(1,1,i),def_mat(2,2,i),def_mat(3,3,i),def_mat(2,3,i),def_mat(1,3,i),def_mat(1,2,i)
            enddo
        endif
    end subroutine apply_elastic_deformation

    subroutine sort_point_symmetries_based_on_determinant_and_trace(R,det,trace,schoenflies)
        !
        use am_rank_and_sort
        !
        implicit none
        !
        real(dp), intent(inout) :: R(:,:,:)
        real(dp), intent(inout) :: det(:)
        real(dp), intent(inout) :: trace(:)
        character(*), intent(inout), optional :: schoenflies(:)
        integer :: nptsym
        integer, allocatable :: indices(:)
        !
        nptsym = size(R,3)
        allocate(indices(nptsym))
        !
        call rank(det,indices);   indices = indices(nptsym:1:-1)
        R=R(:,:,indices); det=det(indices); trace=trace(indices); 
        if (present(schoenflies)) schoenflies=schoenflies(indices)
        !
        call rank(trace,indices); indices = indices(nptsym:1:-1)
        R=R(:,:,indices); det=det(indices); trace=trace(indices); 
        if (present(schoenflies)) schoenflies=schoenflies(indices)
        !
    end subroutine sort_point_symmetries_based_on_determinant_and_trace
    

    
    
    
    
    
    end module




