    module am_unit_cell
    !
    use am_constants
    use am_helpers
    use am_options
    !
    implicit none
    !
    private
    ! classes
    public :: am_class_unit_cell
    public :: am_class_symmetry
    !
    type am_class_symmetry 
        integer :: nptsym !> number of point symmetries
        real(dp), allocatable :: R(:,:,:) !> point symmetries
        real(dp), allocatable :: det(:) !> point symmetry determinants
        real(dp), allocatable :: trace(:) !> point symmetry traces
        character(4), allocatable :: schoenflies(:) !> point symmetry schoenflies notation
        logical :: is_centrosymmetric !> T/F if cell is centrosymmetric    
    contains
        procedure :: determine_symmetry
        procedure :: write_stdout
    end type
    !
    type am_class_unit_cell
        ! produced by load_poscar:
        real(dp) :: bas(3,3) !> column vectors a(1:3,i), a(1:3,j), a(1:3,k)
        integer :: natoms !> number of atoms
        integer :: nspecies !> number of unique atomic species
        character(len=:), allocatable :: symbs(:) !> symbols of unique atomic species
        integer,allocatable :: natoms_per_species(:) !> number of atoms per species
        real(dp), allocatable :: tau_frac(:,:)!> atomic coordinates tau(3,natoms) in fractional
        integer,allocatable :: atype(:) !> index indentifying the atomic species for each atom
        ! produced by compute_basic_properties
        real(dp) :: recbas(3,3) !> reciprocal basis vectors,column vectors a(1:3,i), a(1:3,j), a(1:3,k)
        real(dp) :: metric(3,3) !> metric tensor of real space basis
        real(dp) :: vol !> cell volume
        real(dp), allocatable :: tau_cart(:,:) !> atomic coordinates tau(3,natoms) in cartesian
    contains
        procedure :: load_poscar
        procedure :: reduce_to_primitive
        procedure :: compute_basic_properties
    end type am_class_unit_cell

    contains
  
    ! wrappers
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
        tau_frac=uc%tau_frac,&
        atype=uc%atype,&
        iopt_filename=filename,&
        iopt_verbosity=opts%verbosity)
    !
    end subroutine load_poscar
    subroutine determine_symmetry(ss,uc,iopts)
    !
    ! 1) uses atomic positions to identify possible symmetry translations which could serve as primitive lattice vectors
    ! 2) selects the shortest of these vectors as the basis for the primitive cell
    ! 3) determines the point symmetries compatbile with the metric tensor corresponding to this primitive basis
    !       Symmetry and Condensed Matter Physics: A Computational Approach. 1 edition. Cambridge,
    !       UK; New York: Cambridge University Press, 2008. p 282, table 10.3; p 388, table 10.4
    ! 4) determines the type of lattice from the number of point symmetries identified in the step above
    !       Symmetry and Condensed Matter Physics: A Computational Approach. 1 edition.
    !       Cambridge, UK; New York: Cambridge University Press, 2008. p 282, table 10.3
    !       p 388, table 10.4
    ! 5) of the point symmeties found, identify which ones are compatible with the atoms (position and type)
    ! 6) determine basic point symmetry property properties 
    ! 7) Sort point symmetries based on det and trace (first point symmetry is identity; last point symmetry is inversion, if present)s
    !
    implicit none
    ! intput/outputs
    class(am_class_symmetry), intent(inout) :: ss
    type(am_class_unit_cell), intent(in) :: uc
    type(am_class_unit_cell) :: prim
    type(am_class_symmetry) :: ls ! symmetries compatbile with the lattice (not necessarily atomic basis)
    ! generic point group properties
    integer, allocatable :: indx(:) ! pairs with mask, see below.
    logical, allocatable :: mask(:) ! mask which specifies which lattiec-compatible point symmetries are also compatbile with atomic positions
    character(12) :: crystal_system
    logical :: is_centrosymmetric
    ! loop variables
    integer :: i, j
    ! options
    type(am_class_options), intent(in), optional :: iopts
    type(am_class_options) :: opts
    if (present(iopts)) then 
        opts = iopts
    else
        call opts%defaults
    endif
    !
    ! write opening lines
    if ( opts%verbosity .ge. 1 ) write(*,*)
    if ( opts%verbosity .ge. 1 ) write(*,'(a)' ) 'Analyzing unit cell symmetry'
    !
    ! check that bas is not empty
    if ( all( abs(uc%bas).lt.tiny) ) then
        call am_print('ERROR','Primitive basis is empty! Most likely POSCAR was not loaded.',' >>> ')
        call am_print('basis (column vectors)',uc%bas)
        stop
    endif
    !
    ! 1) uses atomic positions to identify possible symmetry translations which could serve as primitive lattice vectors
    ! 2) selects the shortest of these vectors as the basis for the primitive cell
    !
    call uc%reduce_to_primitive(prim=prim,iopts=opts)
    !
    ! 3) determines the point symmetries compatbile with this primitive basis
    !
    call analyze_metric_tensor_symmetry(metric=prim%metric,R=ls%R)
    ls%nptsym = size(ls%R,3)
    if ( opts%verbosity .ge. 1 ) call am_print("point symmetries compatible with the metric tensor",ls%nptsym," ... ")
    !
    ! 4) determines the type of lattice from the number of point symmetries identified in the step above
    ! 
    lattice_type : select case (ls%nptsym)
    case(2);  crystal_system='triclinic   '
    case(4);  crystal_system='monoclinic  '
    case(8);  crystal_system='orthorhombic'
    case(12); crystal_system='trigonal    '
    case(16); crystal_system='tetragonal  '
    case(24); crystal_system='hexagonal   '
    case(48); crystal_system='cubic       '
    end select lattice_type
    if ( opts%verbosity .ge. 1 ) call am_print('crystal system',trim(crystal_system)," ... ")
    !
    ! 5) Of the point symmeties found, identify which ones are compatible with the atoms (position and type)
    !
    allocate(mask(ls%nptsym))
    allocate(indx(ls%nptsym))
    mask = .false.
    do i = 1,ls%nptsym
        if (is_symmetry_valid(iopt_R=ls%R(:,:,i),tau_frac=uc%tau_cart,iopt_atype=uc%atype,iopt_sym_prec=opts%sym_prec)) then
            mask(i) = .true.
        endif
        indx(i) = i
    enddo
    ss%nptsym = count(mask)
    if ( opts%verbosity .ge. 1 ) call am_print('point symmetries compatbile with atoms',ss%nptsym," ... ")
    allocate(ss%R(3,3,ss%nptsym))
    ss%R = ls%R(:,:,pack(indx,mask))
    !
    ! 6) Determine basic point symmetry property properties 
    !
    call determine_basic_point_symmetry_properties(R=ss%R,det=ss%det,trace=ss%trace,schoenflies=ss%schoenflies)
    !
    ! 7) Sort point symmetries based on det and trace (first point symmetry is identity; last point symmetry is inversion, if present)
    !
    call sort_point_symmetries_based_on_determinant_and_trace(R=ss%R,det=ss%det,trace=ss%trace,schoenflies=ss%schoenflies)
    !
    if ( .not. all(abs(ss%R(:,:,1)-eye(3)).lt.tiny) ) then
        call am_print('ERROR','Something is wrong. Identity is not the first point symmetry.')
        call am_print('R(:,:,1)',ss%R(:,:,i),'     ')
        stop
    endif
    if ( all(abs(ss%R(:,:,ss%nptsym)+eye(3)).lt.tiny) ) then
        ss%is_centrosymmetric = .true.
    else
        ss%is_centrosymmetric = .false.
    endif
    !
    if ( opts%verbosity .ge. 1 ) call am_print('centrosymmetric?',ss%is_centrosymmetric,' ... ')
    if ( opts%verbosity .ge. 1 ) call ss%write_stdout
    !
    end subroutine determine_symmetry
    subroutine reduce_to_primitive(uc,prim,iopts)
    ! 
    ! 1) uses atomic positions to identify possible symmetry translations which could serve as primitive lattice vectors
    ! 2) selects the shortest of these vectors as the basis for the primitive cell
    ! 
    ! essentially a wrapper for determine_primitive_basis, look below: code is short and self explanatory.
    !
    class(am_class_unit_cell), intent(in) :: uc
    type(am_class_unit_cell), intent(out) :: prim
    type(am_class_options), intent(in), optional :: iopts
    type(am_class_options) :: opts
    if (present(iopts)) then 
        opts = iopts
    else
        call opts%defaults
    endif
    !
    call determine_primitive_basis(&
        bas=uc%bas,&
        natoms=uc%natoms,&
        atype=uc%atype,&
        tau_frac=uc%tau_frac,&
        oopt_bas_prim=prim%bas,&
        oopt_natoms_prim=prim%natoms,&
        oopt_atype_prim=prim%atype,&
        oopt_tau_prim_cart=prim%tau_cart,&
        oopt_tau_prim_frac=prim%tau_frac,&
        iopt_verbosity=opts%verbosity)
    !
    call prim%compute_basic_properties(iopts=opts)
    !
    end subroutine reduce_to_primitive
    subroutine expand_to_supercell(uc,n,sc,iopts)
    !
    ! Creates a super cell with dimensions n(1),n(2),n(3) by adding all integer combinations of the vectors  [1...n(1),1...n(2),1...n(3)] to each atomic position. 
    ! 
    ! essentially a wrapper for determine_supercell, look below: code is short and self explanatory.
    !
    integer, intent(in) :: n(3) ! supercell dimensions
    class(am_class_unit_cell), intent(in) :: uc
    type(am_class_unit_cell), intent(out) :: sc
    type(am_class_options), intent(in), optional :: iopts
    type(am_class_options) :: opts
    if (present(iopts)) then 
        opts = iopts
    else
        call opts%defaults
    endif
    !
    call determine_supercell(&
        n=n,&
        bas=uc%bas,&
        natoms=uc%natoms,&
        atype=uc%atype,&
        tau_frac=uc%tau_frac,&
        bas_sc=sc%bas,&
        natoms_sc=sc%natoms,&
        atype_sc=sc%atype,&
        tau_frac_sc=sc%tau_frac,&
        iopt_verbosity=opts%verbosity)
    !
    call sc%compute_basic_properties(iopts=opts)
    !
    end subroutine expand_to_supercell
    subroutine compute_basic_properties(uc,iopts)
    !
    ! inputs:
    ! bas : lattice basis
    ! tau_frac : atomic coordinates (fractional)
    ! 
    ! calculates basic lattice properties
    ! reciprocal basis, metric tensor, cell volume, atomic coordinates (cartesian)
    !
    !
    class(am_class_unit_cell), intent(inout) :: uc
    type(am_class_options), intent(in), optional :: iopts
    type(am_class_options) :: opts
    if (present(iopts)) then 
        opts = iopts
    else
        call opts%defaults
    endif
    !
    call determine_basic_lattice_properties(&
        bas=uc%bas,&
        tau_frac=uc%tau_frac,&
        recbas=uc%recbas,&
        metric=uc%metric,&
        vol=uc%vol,&
        tau_cart=uc%tau_cart,&
        iopt_nspecies=uc%nspecies,&
        iopt_natoms_per_species=uc%natoms_per_species,&
        iopt_atype=uc%atype,&
        iopt_symbs=uc%symbs,&
        iopt_verbosity=opts%verbosity)
    !
    end subroutine compute_basic_properties
    
    subroutine write_stdout(ss)
    !
    implicit none
    !
    class(am_class_symmetry), intent(in) ::  ss
    integer :: i
    !
    write(*,'(a5,a)') ' ... ','point symmetry properties'
    write(*,'(5x,a3,a5,a6,a6)') '#', 'sch.', 'det', 'trace'
    do i = 1, ss%nptsym
        write(*,'(5x,i3,a5,f6.2,f6.2)') i, trim(ss%schoenflies(i)), ss%det(i), ss%trace(i)
    enddo
    !
    end subroutine write_stdout

    !
    ! internal
    !

    pure subroutine analyze_metric_tensor_symmetry(metric,R)
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
    real(dp), intent(in)  :: metric(3,3) !> lattice metric tensor
    real(dp), allocatable, intent(out) :: R(:,:,:) !> point symmetries (fractional coordinates)
    ! internal:
    real(dp) :: buffer(3,3,48) ! buffer for point symmetries
    integer  :: k       !> point symmetry counter
    real(dp) :: o(3,3)  !> point symmetry in fractional
    integer  :: i, i11, i12, i13, i21, i22, i23, i31, i32, i33 !> used to generate unitary rotational matrices
    !
    buffer=0
    k=0
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
                                        o(1,1:3)=[i11-2.0_dp,i12-2.0_dp,i13-2.0_dp]
                                        o(2,1:3)=[i21-2.0_dp,i22-2.0_dp,i23-2.0_dp]
                                        o(3,1:3)=[i31-2.0_dp,i32-2.0_dp,i33-2.0_dp]
                                        ! check that metric is left unchanged (relaxed tiny bounds here; prone to rounding errors)
                                        if ( all( abs(matmul(transpose(o),matmul(metric,o))-metric) .lt. 1E-4 ) ) then
                                            k = k + 1
                                            buffer(1:3,1:3,k) = o
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
    end subroutine analyze_metric_tensor_symmetry
    pure subroutine identify_point_symmetry(trace,det,schoenflies,ierror)
    !
    ! identify the point symmetry (schoenflies notation) given it's trace and determinant.
    ! if it is unable to identify the point symmetry ierror = 1 is returned.
    ! if it exists cleanly ierror = 0
    !
    implicit none
    !
    real(dp), intent(in) :: trace
    real(dp), intent(in) :: det
    character(4), intent(out) :: schoenflies
    integer, intent(out) :: ierror
    !
    ! The Mathematical Theory of Symmetry in Solids: Representation Theory for
    ! Point Groups and Space Groups. 1 edition. Oxford?: New York: Oxford University
    ! Press, 2010. page 138, table 3.8.
    !
    schoenflies = 'NL  '
    if     ( (abs(trace - 3).lt.tiny) .and. (abs(det - 1).lt.tiny) ) then; schoenflies = 'I   '
    elseif ( (abs(trace + 1).lt.tiny) .and. (abs(det - 1).lt.tiny) ) then; schoenflies = 'C2  '
    elseif ( (abs(trace - 0).lt.tiny) .and. (abs(det - 1).lt.tiny) ) then; schoenflies = 'C3  '
    elseif ( (abs(trace - 1).lt.tiny) .and. (abs(det - 1).lt.tiny) ) then; schoenflies = 'C4  '
    elseif ( (abs(trace - 2).lt.tiny) .and. (abs(det - 1).lt.tiny) ) then; schoenflies = 'C6  '
    elseif ( (abs(trace + 3).lt.tiny) .and. (abs(det + 1).lt.tiny) ) then; schoenflies = 'i   '
    elseif ( (abs(trace - 1).lt.tiny) .and. (abs(det + 1).lt.tiny) ) then; schoenflies = 's   '
    elseif ( (abs(trace - 0).lt.tiny) .and. (abs(det + 1).lt.tiny) ) then; schoenflies = 'S6  '
    elseif ( (abs(trace + 1).lt.tiny) .and. (abs(det + 1).lt.tiny) ) then; schoenflies = 'S4  '
    elseif ( (abs(trace + 2).lt.tiny) .and. (abs(det + 1).lt.tiny) ) then; schoenflies = 'S3  '
    endif
    !
    ierror = 0
    if (schoenflies .eq. 'NL  ' ) ierror = 1
    !
    end subroutine identify_point_symmetry
    pure subroutine identify_point_group(R,schoenflies)
    !>
    !> Refs:
    !>
    !> Applied Group Theory: For Physicists and Chemists. Reissue edition.
    !> Mineola, New York: Dover Publications, 2015. page 20.
    !>
    !> Casas, Ignasi, and Juan J. Pérez. “Modification to Flow Chart to
    !> Determine Point Groups.” Journal of Chemical Education 69, no. 1
    !> (January 1, 1992): 83. doi:10.1021/ed069p83.2.
    !>
    !> Breneman, G. L. “Crystallographic Symmetry Point Group Notation
    !> Flow Chart.” Journal of Chemical Education 64, no. 3 (March 1, 1987):
    !> 216. doi:10.1021/ed064p216.
    !>
    !> Adapted from VASP
    !>
    !>   pgind -->sfnames (Schoenflies)
    !>     1 --> c_1       9 --> c_3      17 --> d_4      25 --> c_6v
    !>     2 --> s_2      10 --> s_6      18 --> c_4v     26 --> d_3h
    !>     3 --> c_2      11 --> d_3      19 --> d_2d     27 --> d_6h
    !>     4 --> c_1h     12 --> c_3v     20 --> d_4h     28 --> t
    !>     5 --> c_2h     13 --> d_3d     21 --> c_6      29 --> t_h
    !>     6 --> d_2      14 --> c_4      22 --> c_3h     30 --> o
    !>     7 --> c_2v     15 --> s_4      23 --> c_6h     31 --> t_d
    !>     8 --> d_2h     16 --> c_4h     24 --> d_6      32 --> o_h
    !>
    implicit none
    !
    real(dp), intent(in) :: R(:,:,:)
    character(4), intent(out) :: schoenflies
    integer :: pgind
    integer :: trace, det, nrot
    integer :: invers, nc2, nc3, nc4, nc6, ns2, ns6, ns4, ns3, ir
    character(4), parameter :: schoenflies_database(32) = (/           &
        'c_1 ','s_2 ','c_2 ','c_1h','c_2h','d_2 ','c_2v',&
        'd_2h','c_3 ','s_6 ','d_3 ','c_3v','d_3d','c_4 ',&
        's_4 ','c_4h','d_4 ','c_4v','d_2d','d_4h','c_6 ',&
        'c_3h','c_6h','d_6 ','c_6v','d_3h','d_6h','t   ',&
        't_h ','o   ','t_d ','o_h ' /)
    !
    pgind = 0
    schoenflies='NL'
    !
    ! trivial cases which have unique numbers of rotation operations
    ! c_1 (pgind=1), o_h (pgind=32), d_4h (pgind=20), and c_3 (pgind=9)
    !
    nrot = size(R,3)
    if (nrot .eq. 1) then
        pgind=1
        schoenflies=schoenflies_database(pgind)
        return
    elseif (nrot .eq. 48) then
        pgind=32
        schoenflies=schoenflies_database(pgind)
        return
    elseif (nrot .eq. 16) then
        pgind=20
        schoenflies=schoenflies_database(pgind)
        return
    elseif (nrot .eq. 3) then
        pgind=9
        schoenflies=schoenflies_database(pgind)
        return
    endif
    !
    ! all other groups need further investigations and detailed analysis ...
    ! first determine the type of elements and count them. possible elements
    ! are e, i, c_2,3,4,6 and s_2,3,4,6 (s_2 = m), e is trivial and always
    ! present. the type of a symmetry operation can be identified simply by
    ! calculating the trace and the determinant of the rotation matrix. the
    ! combination of these two quantities is specific for specific elements:
    !
    ! element:         e    i  c_2  c_3  c_4  c_6  s_2  s_6  s_4  s_3
    ! trace:          +3   -3   -1    0   +1   +2   +1    0   -1   -2
    ! determinant:    +1   -1   +1   +1   +1   +1   -1   -1   -1   -1
    !
    invers=0; nc2=0; nc3=0; nc4=0; nc6=0; ns2=0; ns6=0; ns4=0; ns3=0
    !
    do ir = 1, nrot
        trace = R(1,1,ir)+R(2,2,ir)+R(3,3,ir)
        if (abs(trace-3).lt.tiny) then
            ! found unity operator (trivial).
        elseif (abs(trace+3).lt.tiny) then
            ! found inversion
            invers=1
        else
            det=R(1,1,ir)*(R(2,2,ir)*R(3,3,ir)-R(2,3,ir)*R(3,2,ir))+ &
                &   R(1,2,ir)*(R(2,3,ir)*R(3,1,ir)-R(2,1,ir)*R(3,3,ir))+ &
                &   R(1,3,ir)*(R(2,1,ir)*R(3,2,ir)-R(2,2,ir)*R(3,1,ir))
            if ( (abs(trace - -1) .lt. tiny) .and. (abs(det -  1) .lt. tiny) ) nc2=nc2+1 ! found c_2
            if ( (abs(trace -  1) .lt. tiny) .and. (abs(det - -1) .lt. tiny) ) ns2=ns2+1 ! found s_2
            if ( (abs(trace -  0) .lt. tiny) .and. (abs(det -  1) .lt. tiny) ) nc3=nc3+1 ! found c_3
            if ( (abs(trace -  0) .lt. tiny) .and. (abs(det - -1) .lt. tiny) ) ns6=ns6+1 ! found s_6
            if ( (abs(trace -  1) .lt. tiny) .and. (abs(det -  1) .lt. tiny) ) nc4=nc4+1 ! found c_4
            if ( (abs(trace - -1) .lt. tiny) .and. (abs(det - -1) .lt. tiny) ) ns4=ns4+1 ! found s_4
            if ( (abs(trace -  2) .lt. tiny) .and. (abs(det -  1) .lt. tiny) ) nc6=nc6+1 ! found c_6
            if ( (abs(trace - -2) .lt. tiny) .and. (abs(det - -1) .lt. tiny) ) ns3=ns3+1 ! found s_3
        endif
    enddo
    ! now we know which elements we have and so we know the group ... :
    if (nrot.eq.2) then
        ! groups with 2 elements:
        if (invers.eq.1) then
            ! contains inversion --> s_2 (pgind=2):
            pgind=2
            schoenflies=schoenflies_database(pgind)
            return
        end if
        if (nc2.eq.1) then
            ! contains twofold rotation --> c_2 (pgind=3):
            pgind=3
            schoenflies=schoenflies_database(pgind)
            return
        end if
        if (ns2.eq.1) then
            ! contains mirror plane --> c_1h (pgind=4):
            pgind=4
            schoenflies=schoenflies_database(pgind)
            return
        end if
    end if
    if (nrot.eq.4) then
        ! groups with 4 elements:
        if (invers.eq.1) then
            ! contains inversion --> c_2h (pgind=5):
            pgind=5
            schoenflies=schoenflies_database(pgind)
            return
        end if
        if (nc2.eq.3) then
            ! contains three twofold rotations --> d_2 (pgind=6):
            pgind=6
            schoenflies=schoenflies_database(pgind)
            return
        end if
        if (ns2.eq.2) then
            ! contains two mirror planes --> c_2v (pgind=7):
            pgind=7
            schoenflies=schoenflies_database(pgind)
            return
        end if
        if (nc4.eq.1) then
            ! contains fourfold rotation --> c_4 (pgind=14):
            pgind=14
            schoenflies=schoenflies_database(pgind)
            return
        end if
        if (ns4.eq.2) then
            ! contains fourfold improper rotation --> s_4 (pgind=15):
            pgind=15
            schoenflies=schoenflies_database(pgind)
            return
        end if
    end if
    if (nrot.eq.6) then
        ! groups with 6 elements:
        if (invers.eq.1) then
            ! contains inversion --> s_6 (pgind=10):
            pgind=10
            schoenflies=schoenflies_database(pgind)
            return
        end if
        if (nc2.eq.3) then
            ! contains three twofold rotations --> d_3 (pgind=11):
            pgind=11
            schoenflies=schoenflies_database(pgind)
            return
        end if
        if (ns2.eq.3) then
            ! contains three mirror planes --> c_3v (pgind=12):
            pgind=12
            schoenflies=schoenflies_database(pgind)
            return
        end if
        if (nc2.eq.1) then
            ! contains only (1._q,0._q) twofold rotations --> c_6 (pgind=21):
            pgind=21
            schoenflies=schoenflies_database(pgind)
            return
        end if
        if (ns2.eq.1) then
            ! contains only (1._q,0._q) mirror plane --> c_3h (pgind=22):
            pgind=22
            schoenflies=schoenflies_database(pgind)
            return
        end if
    end if
    if (nrot.eq.8) then
        ! groups with 8 elements:
        if (ns2.eq.3) then
            ! contains three mirror planes --> d_2h (pgind=8):
            pgind=8
            schoenflies=schoenflies_database(pgind)
            return
        end if
        if (ns2.eq.1) then
            ! contains (1._q,0._q) mirror planes --> c_4h (pgind=16):
            pgind=16
            schoenflies=schoenflies_database(pgind)
            return
        end if
        if (ns2.eq.0) then
            ! contains no mirror planes --> d_4 (pgind=17):
            pgind=17
            schoenflies=schoenflies_database(pgind)
            return
        end if
        if (ns2.eq.4) then
            ! contains four mirror planes --> c_4v (pgind=18):
            pgind=18
            schoenflies=schoenflies_database(pgind)
            return
        end if
        if (ns2.eq.2) then
            ! contains two mirror planes --> d_2d (pgind=19):
            pgind=19
            schoenflies=schoenflies_database(pgind)
            return
        end if
    end if
    if (nrot.eq.12) then
        ! groups with 12 elements:
        if (ns2.eq.3) then
            ! contains three mirror planes --> d_3d (pgind=13):
            pgind=13
            schoenflies=schoenflies_database(pgind)
            return
        end if
        if (ns2.eq.1) then
            ! contains (1._q,0._q) mirror planes --> c_6h (pgind=23):
            pgind=23
            schoenflies=schoenflies_database(pgind)
            return
        end if
        if (nc2.eq.7) then
            ! contains seven twofold rotations --> d_6 (pgind=24):
            pgind=24
            schoenflies=schoenflies_database(pgind)
            return
        end if
        if (ns2.eq.6) then
            ! contains six mirror planes --> c_6v (pgind=25):
            pgind=25
            schoenflies=schoenflies_database(pgind)
            return
        end if
        if (ns2.eq.4) then
            ! contains four mirror planes --> d_3h (pgind=26):
            pgind=26
            schoenflies=schoenflies_database(pgind)
            return
        end if
        if (nc3.eq.8) then
            ! contains eight threefold rotations --> t (pgind=28):
            pgind=28
            schoenflies=schoenflies_database(pgind)
            return
        end if
    end if
    if (nrot.eq.24) then
        ! groups with 24 elements:
        if (nc6.eq.2) then
            ! contains two sixfold rotations --> d_6h (pgind=27):
            pgind=27
            schoenflies=schoenflies_database(pgind)
            return
        end if
        if (invers.eq.1) then
            ! contains inversion --> t_h (pgind=29):
            pgind=29
            schoenflies=schoenflies_database(pgind)
            return
        end if
        if (nc4.eq.6) then
            ! contains six fourfold rotations --> o (pgind=30):
            pgind=30
            schoenflies=schoenflies_database(pgind)
            return
        end if
        if (ns4.eq.6) then
            ! contains six fourfold improper rotations --> t_d (pgind=31):
            pgind=31
            schoenflies=schoenflies_database(pgind)
            return
        end if
    end if
    end subroutine identify_point_group
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
    real(dp), allocatable :: tau(:,:) ! atomic coordinates tau(3,natoms) read in
    integer :: natoms
    ! loop variables
    integer :: i, j, m, n
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
    if (verbosity .ge. 1) write(*,*)
    if (verbosity .ge. 1) write(*,'(a)' ) 'Determining basic lattice properties'
    !
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
    subroutine determine_supercell(n,bas,natoms,atype,tau_frac,bas_sc,natoms_sc,atype_sc,tau_frac_sc,iopt_verbosity)
        !
        ! expands fractional atomic coordinates onto an n(1) x n(2) x n(3) supercell
        ! also maps the type of atom from atype (in the primitive cell) to atype_sc (in the supercell).
        !
        implicit none
        !
        integer,intent(in) :: n(3) ! supercell dimensions
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
        if (verbosity .ge. 1) call am_print('dimensions',n,' ... ')
        !
        bas_sc = matmul(bas,1.0_dp*reshape([n(1),0,0,0,n(2),0,0,0,n(3)],[3,3]))
        !
        if (verbosity .ge. 1)call am_print('supercell basis',bas_sc,' ... ')
        !
        natoms_sc = n(1)*n(2)*n(3)*natoms
        allocate(tau_frac_sc(3,natoms_sc))
        allocate(atype_sc(natoms_sc))
        !
        if (verbosity .ge. 1)call am_print('number of atoms in supercell',natoms_sc,' ... ')
        !
        tau_frac_sc = 0.0_dp
        m=0
        do i1 = 1, n(1)
            do i2 = 1, n(2)
                do i3 = 1, n(3)
                    do j = 1,natoms
                        m=m+1
                        tau_frac_sc(1:3,m) = ( tau_frac(1:3,j) + ([i1,i2,i3]-1)*1.0_dp ) / ( n(1:3) *1.0_dp )
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
    subroutine determine_primitive_basis(bas,natoms,atype,tau_frac, &
        oopt_bas_prim, &
        oopt_natoms_prim, &
        oopt_tau_prim_frac, &
        oopt_tau_prim_cart, &
        oopt_atype_prim, &
        oopt_rlist_frac, &
        oopt_rlist_cart, &
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
    real(dp), intent(out), optional :: oopt_bas_prim(3,3) !> basis of primitive cell
    integer , intent(out), optional :: oopt_natoms_prim
    real(dp), intent(out), optional, allocatable :: oopt_tau_prim_frac(:,:)
    real(dp), intent(out), optional, allocatable :: oopt_tau_prim_cart(:,:)
    real(dp), intent(out), optional, allocatable :: oopt_rlist_frac(:,:)
    real(dp), intent(out), optional, allocatable :: oopt_rlist_cart(:,:)
    integer , intent(out), optional, allocatable :: oopt_atype_prim(:)
    ! internal
    real(dp) :: bas_prim(3,3) !> basis of primitive cell
    integer  :: natoms_prim
    real(dp) :: tau_ref_frac(3) ! fractional coordinates of reference atom
    real(dp), allocatable :: tau_cart(:,:)
    real(dp), allocatable :: tau_prim_frac(:,:)
    real(dp), allocatable :: tau_prim_cart(:,:)
    integer , allocatable :: atype_prim(:)
    real(dp), allocatable :: rlist_frac(:,:) ! rlist(1:3,n) list of possible lattice vectors that are symmetry-compatible volume
    real(dp), allocatable :: rlist_cart_mag(:) ! magnitude of possible lattice vectors
    integer , allocatable :: indices(:) ! used to sort based on magnitude
    real(dp) :: bpscf(3,3) ! primitive basis (fractional, in terms of the supercell basis)
    real(dp) :: bscpf(3,3) ! supercell basis (fractional, in terms of the primitive basis)
    real(dp) :: wrkbas(3,3)! workspace for determining lattice vectors corresponding to the smallest primitive cell volume
    real(dp) :: wrkvol     ! workspace for determining lattice vectors corresponding to the smallest primitive cell volume
    real(dp) :: wrkr(3)    ! workspace for determining if the translation is a valid symmetry operation
    integer :: nlatvec ! number of primitive lattice vector candidates
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
    if (verbosity .ge. 1) write(*,*)
    if (verbosity .ge. 1) write(*,*) 'Determining primitive cell'
    !
    if ( natoms .ne. size(tau_frac,2) ) then
        call am_print('ERROR','natoms .ne. size(tau_frac,2)',' >>> ')
        stop
    endif
    !
    !
    !
    ! search for lattice vectors using as translational components the vectors connecting the choosen atom to all other atoms of the same species
    allocate(rlist_frac(3,natoms-1+3))
    allocate(rlist_cart_mag(natoms-1+3))
    rlist_frac = 0.0_dp
    rlist_cart_mag = 0.0_dp
    ! choose one atom (ideally choose an atom corresponding to the species with the fewest number of atoms in unit cell)
    i = 1
    tau_ref_frac(1:3) =  tau_frac(1:3,i)
    ! get the shifted cartesian coordinate (shifted relative to reference atom), used later when determining the type of atom
    allocate(tau_cart(3,natoms))
    do j = 1,natoms
        tau_cart(1:3,j) = matmul(bas,tau_frac(1:3,j) - tau_ref_frac)
    enddo
    ! perform the search
    nlatvec = 0
    do j = 1,natoms
        if ( i .ne. j) then
        if ( atype(i) .eq. atype(j) ) then
            ! shift to put reference atom at zero.
            wrkr(1:3) = tau_frac(1:3,j) - tau_ref_frac
            if ( is_symmetry_valid(iopt_T=wrkr,iopt_atype=atype, tau_frac=tau_frac, iopt_sym_prec=sym_prec) ) then
                ! lattice vector found
                nlatvec = nlatvec + 1
                rlist_frac(1:3,nlatvec) = wrkr
                ! rlist_cart_mag(nlatvec) = sqrt( dot_product( matmul(bas,rlist_frac(1:3,nlatvec)), matmul(bas,rlist_frac(1:3,nlatvec))) )
                rlist_cart_mag(nlatvec) = sqrt( dot_product( rlist_frac(1:3,nlatvec), rlist_frac(1:3,nlatvec)) )
            endif
        endif
        endif
    enddo
    ! add the original lattice vectors
    nlatvec = nlatvec+1
    rlist_frac(1:3,nlatvec) = [1,0,0]*1.0_dp
    rlist_cart_mag(nlatvec) = sqrt( dot_product( matmul(bas,rlist_frac(1:3,nlatvec)), matmul(bas,rlist_frac(1:3,nlatvec))) )
    nlatvec = nlatvec+1
    rlist_frac(1:3,nlatvec) = [0,1,0]*1.0_dp
    rlist_cart_mag(nlatvec) = sqrt( dot_product( matmul(bas,rlist_frac(1:3,nlatvec)), matmul(bas,rlist_frac(1:3,nlatvec))) )
    nlatvec = nlatvec+1
    rlist_frac(1:3,nlatvec) = [0,0,1]*1.0_dp
    rlist_cart_mag(nlatvec) = sqrt( dot_product( matmul(bas,rlist_frac(1:3,nlatvec)), matmul(bas,rlist_frac(1:3,nlatvec))) )
    ! sort lattice vectors found, smallest vector magntiude last.
    allocate(indices(nlatvec))
    call rank(rlist_cart_mag,indices)
    indices = indices(nlatvec:1:-1)
    rlist_frac = rlist_frac(1:3,indices)
    if (verbosity .ge. 1) call am_print("number of possible primitive lattice vectors found",nlatvec," ... ")
    if (verbosity .ge. 1) call am_print("possible primitive lattice vectors (fractional; column vectors, tranposed for clarity)",transpose(rlist_frac)," ... ")
    ! print some things to stdout
    if ( present(oopt_rlist_frac) ) then
        allocate(oopt_rlist_frac(3,nlatvec))
        oopt_rlist_frac = rlist_frac(1:3,1:nlatvec)
    endif
    if ( present(oopt_rlist_cart) ) then
        allocate(oopt_rlist_cart(3,nlatvec))
        do i = 1,nlatvec
            oopt_rlist_cart(1:3,i) = matmul(bas,rlist_frac(1:3,i))
        enddo
        if (verbosity .ge. 1) call am_print("possible primitive lattice vectors (cartesian; column vectors, tranposed for clarity)",transpose(oopt_rlist_cart)," ... ")
    endif
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
        i = nlatvec
        ! loop until a second non-collinear lattice vector is found (non-zero cross product)
        do j = nlatvec,1,-1
            if ( any(cross_product(rlist_frac(1:3,i),rlist_frac(1:3,j)) .gt. tiny) ) exit
        enddo
        if (any([i,j].eq.0)) then
            call am_print('ERROR','No pair of non-collinear lattice vectors found',' >>> ')
            stop
        endif
        ! loop until third lattice vector is found which yields a non-zero cell volume
        do k = nlatvec,1,-1
            wrkbas = rlist_frac(1:3,[i,j,k])
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
        bpscf = rlist_frac(1:3,[i,j,k])
        !
        if (verbosity .ge. 1) write(*,'(5x,a)') 'selecting three vectors with the smallest magnitudes'
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
            allocate(atype_prim(natoms_prim))
            do i = 1,natoms_prim
                atype_prim(i) = 0
                search_for_similar_atom : do j = 1,natoms
                    if ( all(abs(tau_prim_cart(1:3,i)-tau_cart(1:3,j)).lt.sym_prec) ) then
                        atype_prim(i) = atype(j)
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
        endif
    endif
    end subroutine determine_primitive_basis
    function is_symmetry_valid(iopt_R,iopt_T,iopt_atype,tau_frac,iopt_sym_prec)
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
    if ( (translational_symmetry .eq. .false.) .and. (point_symmetry .eq. .false.) ) then
        call am_print('ERROR','Translational and/or rotational symmetry must be supplied',' >>> ')
        stop
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
    else
        ! either atleast one atom was not found to overlap
        ! and the break statement in the check_overlap loop did not catch it.
        call am_print('ERROR','Something is wrong in function is_symmetry_valid',' >>> ')
        stop
    endif
    end function is_symmetry_valid
    subroutine apply_elastic_deformation(bas,deformation_code,strain_max,nstrains,bas_def,iopt_verbosity)
    !
    ! creates at bas_def(3,3,n) matrix containing the n deformed primitive basises.
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




