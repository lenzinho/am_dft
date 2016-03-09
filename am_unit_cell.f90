    module am_unit_cell
    !
    use am_constants
    use am_helpers
    use am_options
    !
    implicit none
    !
    private
    !
    public :: am_class_unit_cell
    public :: space_symmetries_from_basis ! used by am_symmetry
    public :: reciprocal_basis ! used by am_symmetry
    
    type am_class_unit_cell
        real(dp) :: bas(3,3) !> column vectors a(1:3,i), a(1:3,j), a(1:3,k)
        integer  :: natoms   !> number of atoms
        integer  :: nspecies !> number of unique atomic species
        character(len=:), allocatable :: symbs(:) !> symbs(nspecies) symbols of unique atomic species
        real(dp), allocatable :: tau(:,:) !> tau(3,natoms) cartesian atomic coordinates 
        integer , allocatable :: atype(:) !> index indentifying the atomic species for each atom
    contains
        procedure :: load_poscar
        procedure :: reduce_to_primitive
        procedure :: expand_to_supercell
        procedure :: primitive_to_conventional
    end type am_class_unit_cell

    contains

    ! functions which operate on bas(3,3)

    function      primitive_volume(bas,skip_volume_check) result(vol)
        !
        implicit none
        !
        real(dp), intent(in) :: bas(3,3)
        logical , intent(in), optional :: skip_volume_check
        real(dp) :: vol
        !
        vol=abs(bas(1,1)*(bas(2,2)*bas(3,3)-bas(3,2)*bas(2,3)) &
             & +bas(2,1)*(bas(3,2)*bas(1,3)-bas(1,2)*bas(3,3)) &
             & +bas(3,1)*(bas(1,2)*bas(2,3)-bas(1,3)*bas(2,2)))
        !
        if (.not.present(skip_volume_check)) then
        if ( abs(vol) .lt. tiny ) then
            call am_print("ERROR","Basis vectors are not coplanar!"," >>> ")
            stop
        endif
        endif
        !
        if (vol .lt. 0.0_dp) then
            call am_print("WARNING","Primitive cell volume is negative!"," >>> ")
        endif
        !
    end function  primitive_volume

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

    pure function metric_from_basis(bas) result(metric)
        !
        implicit none
        !
        real(dp), intent(in) :: bas(3,3)
        real(dp) :: metric(3,3)
        !
        metric = matmul(transpose(bas),bas)
        !
    end function  metric_from_basis

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
        real(dp) :: recbas(3,3)
        real(dp) :: metric(3,3)
        real(dp) :: buffer(3,3,48) ! buffer for point symmetries
        integer  :: k ! point symmetry counter
        real(dp) :: o(3,3) ! point symmetry in fractional (converted to cartesian before output)
        integer  :: i, i11, i12, i13, i21, i22, i23, i31, i32, i33 ! used to generate unitary rotational matrices
        !
        recbas = reciprocal_basis(bas)
        !
        metric = metric_from_basis(bas)
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

    function      apply_elastic_deformation(bas,deformation_code,strain_max,nstrains,iopt_verbosity) result(bas_def)
        !
        !
        ! needs editing to conform. make this a "pure function" , push all outputs except errors into its caller.
        !
        ! creates at bas_def(3,3,nTs) matrix containing the nTs deformed primitive basises.
        !
        ! bas_def = apply_elastic_deformation(bas=uc%bas,deformation_code=6,strain_max=0.1_dp,nstrains=10,iopt_verbosity=verbosity)
        !
        implicit none
        !
        real(dp), intent(in) :: bas(3,3)
        integer , intent(in) :: deformation_code
        real(dp), intent(in) :: strain_max
        integer , intent(in) :: nstrains
        real(dp), allocatable :: bas_def(:,:,:)
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
    end function  apply_elastic_deformation

    function      lattice_type(bas) result(lattice_code)
        !>
        !> get crystal system centering from the number of point symmetries and metric tensor. 
        !> For the formalism, see:
        !>
        !> Symmetry and Condensed Matter Physics: A Computational Approach. 1 edition. Cambridge,
        !> UK; New York: Cambridge University Press, 2008. p 282, table 10.3; p 288, table 10.4
        !>
        !>     code  lttice type
        !>        1  simple cubic cell                            
        !>        2  body-centered cubic cell                     
        !>        3  face-centered cubic cell                     
        !>        4  hexagonal cell                               
        !>        5  simple tetragonal cell                       
        !>        6  body-centered tetragonal cell                
        !>        7  rhombohedral (trigonal) cell                 
        !>        8  simple orthorhombic cell                     
        !>        9  body-centered orthorhombic cell              
        !>       10  face-centered orthorhombic cell              
        !>       11  base centered orthorhombic cell              
        !>       12  simple monoclinic cell                       
        !>       13  base-centered monoclinic cell                
        !>       14  triclinic cell                              
        !> 
        implicit none
        !
        real(dp), intent(in) :: bas(3,3) ! number of point group symmetries
        integer  :: lattice_code
        integer  :: nsyms ! number of symmetries
        real(dp) :: M(3,3) ! metric tensor
        !
        nsyms = size( determine_lattice_point_symmetries(bas), 3)
        !
        M = metric_from_basis(bas)
        !
        lattice_code = 0
        !
        select case (nsyms)
        case(2) ! triclinic
            lattice_code = 14
        case(4) ! monoclinic
            if (are_equal( [0.0_dp,M(1,3),M(2,3)] )) then
            if (are_different( [M(1,1),M(2,2),M(3,3)] )) then
                lattice_code = 12
                return
            endif
            endif
            if (are_equal( [M(1,1),M(2,2)] )) then
            if (are_different( [M(1,2),M(1,3),M(2,3),M(3,3)] )) then
                lattice_code = 13
                return
            endif
            endif
        case(8) 
            ! simple orthorhombic
            if (are_equal( [0.0_dp,M(1,2),M(1,3),M(2,3) ] )) then
            if (are_different( [M(1,1),M(2,2),M(3,3)] )) then
                lattice_code = 8
                return
            endif
            endif
            ! body-centered orthorhombic
            if (are_equal( [M(1,1),M(2,2),M(3,3) ] )) then
            if (are_different( [M(1,2),M(1,3),M(2,3)] )) then
                lattice_code = 9
                return
            endif
            endif
            ! based-centered orthorhombic
            if (are_equal( [M(1,1),M(2,2)] ).and.are_equal( [0.0_dp,M(1,3),M(2,3)] )) then
            if (are_different( [M(2,2),M(3,3),M(2,3),M(1,2)] )) then
                lattice_code = 11
                return
            endif
            endif
            ! face-centered orthorhombic
            if (are_different( [M(1,1),M(2,2),M(3,3),M(2,3),M(1,3),M(1,2)] )) then
                lattice_code = 10
                return
            endif
        case(12)! trigonal (rhombahedral)
            !  based-centered trigonal
            if (are_equal( [M(1,1),M(2,2),M(3,3)] ).and.are_equal( [M(1,2),M(1,3),M(2,3)] )) then
            if (are_different( [M(1,1),M(2,3)] )) then
                lattice_code = 7
                return
            endif
            endif
        case(16)! tetragonal
            ! simple tetragonal
            if (are_equal( [M(1,1),M(2,2)] ).and.are_equal( [0.0_dp,M(1,2),M(1,3),M(2,3)] )) then
            if (are_different( [M(1,1),M(3,3),M(1,2)] )) then
                lattice_code = 5
                return
            endif
            endif
            ! body-centered tetragonal
            if (are_equal( [M(1,1),M(2,2),M(3,3)] ).and.are_equal( [M(1,3),M(2,3)] )) then
            if (are_different( [M(1,1),M(1,2),M(1,3)] )) then
                lattice_code = 5
                return
            endif
            endif
        case(24)! hexagonal
            ! simple hexagonal
            if (are_equal( [M(1,1),M(2,2),2.0_dp*M(1,2)] ).and.are_equal( [0.0_dp,M(1,3),M(2,3)] )) then
            if (are_different( [M(2,2),M(3,3),M(2,3)] )) then
                lattice_code = 4
                return
            endif
            endif
        case(48)! cubic
            ! simple cubic
            if (are_equal( [M(1,1),M(2,2),M(3,3)] ).and.are_equal( [0.0_dp,M(1,2),M(1,3),M(2,3)] )) then
            if (are_different( [M(1,1),M(1,2)] )) then
                lattice_code = 1
                return
            endif
            endif
            ! face-centered cubic
            if (are_equal( [M(1,1),M(2,2),M(3,3),2.0_dp*M(1,2),2.0_dp*M(1,3),2.0_dp*M(2,3)] )) then
                lattice_code = 3
                return
            endif
            ! body-centered cubic
            if (are_equal( [M(1,1),M(2,2),M(3,3),-3.0_dp*M(1,2),-3.0_dp*M(1,3),-3.0_dp*M(2,3)] )) then
                lattice_code = 2
                return
            endif
        end select
        !
        if (lattice_code.eq.0) then
            call am_print('ERROR','Unable to identify lattice type.',' >>> ')
            call am_print('number of symmetry operations',nsyms)
            call am_print('primitive basis',bas)
            call am_print('metric tensor',M)
            stop
        endif
        !
    end function  lattice_type

    !
    ! functions which operate on lattice_code
    !

    function      lattice_centering(lattice_code) result(centering)
        ! 
        implicit none
        !
        integer  :: lattice_code
        real(dp), dimension(3,3) :: P,I,F,A,B,C,centering
        !
        ! simple/primitive
        P(1:3,1) = [1,0,0]
        P(1:3,2) = [0,1,0]
        P(1:3,3) = [0,0,1]
        ! body centered cell
        I(1:3,1) = [-1,1,1]/2.0_dp
        I(1:3,2) = [1,-1,1]/2.0_dp
        I(1:3,3) = [1,1,-1]/2.0_dp
        ! face centered
        F(1:3,1) = [0,1,1]/2.0_dp
        F(1:3,2) = [1,0,1]/2.0_dp
        F(1:3,3) = [1,1,0]/2.0_dp
        ! A
        A(1:3,1) = [ 0,0,2]/2.0_dp
        A(1:3,2) = [-1,1,0]/2.0_dp
        A(1:3,3) = [ 1,1,0]/2.0_dp
        ! B
        B(1:3,1) = [-1,1,0]/2.0_dp
        B(1:3,2) = [ 0,0,2]/2.0_dp
        B(1:3,3) = [ 1,1,0]/2.0_dp
        ! C
        C(1:3,1) = [-1,1,0]/2.0_dp
        C(1:3,2) = [ 1,1,0]/2.0_dp
        C(1:3,3) = [ 0,0,2]/2.0_dp

        select case(lattice_code)
        case(1);  centering = P ! simple cubic cell                            
        case(2);  centering = I ! body-centered cubic cell                     
        case(3);  centering = F ! face-centered cubic cell                     
        case(4);  centering = P ! hexagonal cell                               
        case(5);  centering = P ! simple tetragonal cell                       
        case(6);  centering = I ! body-centered tetragonal cell                
        case(7);  centering = A ! rhombohedral (trigonal) cell
        case(8);  centering = P ! simple orthorhombic cell                     
        case(9);  centering = I ! body-centered orthorhombic cell              
        case(10); centering = F ! face-centered orthorhombic cell              
        case(11); centering = A ! base centered orthorhombic cell              
        case(12); centering = P ! simple monoclinic cell                       
        case(13); centering = A ! base-centered monoclinic cell                
        case(14); centering = P ! triclinic cell      
        case default
            call am_print('ERROR','Unknown lattice code.')        
        end select
        !
    end function  lattice_centering

    !
    ! functions which operate on atomic basis
    !

    pure function translations_from_basis(bas,tau,atype,iopt_include,iopt_sym_prec) result(T)
        !
        implicit none
        ! subroutine i/o
        real(dp)        , intent(in)           :: bas(3,3)
        integer         , intent(in)           :: atype(:) !> atype(natoms) list identifing type of atom
        real(dp)        , intent(in)           :: tau(:,:) !> tau(3,natoms) cartesian atomic coordinates
        character(len=*), intent(in), optional :: iopt_include !> iopt_include can be 'primzero', 'prim', 'zero'; if present it adds [0,0,0] and the primitive basis to the T vectors returned. 
        real(dp)        , allocatable          :: T(:,:)   !> T(3,nTs) translation which leaves basis invariant
        ! internal
        integer  :: natoms ! natoms number of atoms
        integer  :: nTs ! nTs number of primitive lattice vector candidates
        real(dp) :: tau_ref(3) ! tau_ref(3) cartesian coordinates of reference atom
        real(dp), allocatable :: wrk(:,:) ! wrk(1:3,nTs) list of possible lattice vectors that are symmetry-compatible volume
        real(dp) :: wrkr(3)
        ! loop variables
        integer :: i, j
        ! optionals
        real(dp), intent(in), optional :: iopt_sym_prec
        real(dp) :: sym_prec
        if ( present(iopt_sym_prec) ) then
            sym_prec = iopt_sym_prec
        else
            sym_prec = tiny
        endif
        !
        !
        !
        natoms = size(tau,2)
        allocate(wrk(3,natoms-1+3)) 
        wrk = 0.0_dp
        ! choose one atom (ideally choose an atom corresponding to the species with the fewest number of atoms in unit cell)
        i = 1
        tau_ref(1:3) =  tau(1:3,i)
        ! search for lattice vectors using, as translational components, vectors connecting the choosen atom to other atoms of the same species
        nTs = 0
        do j = 1,natoms
            if ( i .ne. j) then
            if ( atype(i) .eq. atype(j) ) then
                ! shift to put reference atom at zero.
                wrkr(1:3) = tau(1:3,j) - tau_ref
                !
                if ( is_symmetry_valid(bas=bas, iopt_T=wrkr, atype=atype, tau=tau, iopt_sym_prec=sym_prec) ) then
                    ! lattice vector found
                    nTs = nTs + 1
                    wrk(1:3,nTs) = wrkr
                endif
            endif
            endif
        enddo
        if (present(iopt_include)) then
            if (index(iopt_include,"prim").ne.0) then
                ! add native basis to vectors
                nTs = nTs+1; wrk(1:3,nTs) = bas(:,1)
                nTs = nTs+1; wrk(1:3,nTs) = bas(:,2)
                nTs = nTs+1; wrk(1:3,nTs) = bas(:,3)
            endif
            if (index(iopt_include,"zero").ne.0) then
                ! add origin as a possibletranslation
                nTs = nTs+1; wrk(1:3,nTs) = real([0,0,0],dp)
            endif
        endif
        ! allocate output
        allocate(T(3,nTs))
        T = wrk(1:3,1:nTs)
    end function  translations_from_basis

    pure function is_symmetry_valid(bas,tau,atype,iopt_R,iopt_T,iopt_sym_prec)
        !
        ! check whether symmetry operation is valid
        !
        implicit none
        ! function i/o
        real(dp), intent(in)           :: bas(3,3)
        real(dp), intent(in)           :: tau(:,:) !> tau_frac(3,natoms) cartesian atomic basis
        integer , intent(in), optional :: atype(:) !> atype(natoms) list identify type of atom
        real(dp), intent(in), optional :: iopt_R(3,3)
        real(dp), intent(in), optional :: iopt_T(3)
        real(dp), intent(in), optional :: iopt_sym_prec
        !
        real(dp) :: recbas(3,3)
        real(dp), allocatable :: tau_frac(:,:) !> tau_frac_frac(3,natoms) fractional atomic basis
        real(dp) :: R(3,3)
        real(dp) :: T(3)
        real(dp) :: sym_prec
        real(dp) :: tau_frac_rot(3) ! rotated 
        logical :: is_symmetry_valid
        integer :: natoms
        integer :: i, j, m
        logical :: overlap_found
        !
        natoms = size(tau_frac,2)
        !
        if ( present(iopt_sym_prec) ) then
            sym_prec = iopt_sym_prec
        else
            sym_prec = tiny
        endif
        !
        if ( present(iopt_R) ) then
            R = iopt_R
        else
            R(1:3,1) = real([1,0,0],dp)
            R(1:3,2) = real([0,1,0],dp)
            R(1:3,3) = real([0,0,1],dp)
        endif
        !
        if ( present(iopt_T) ) then
            T = iopt_T
        else
            T = real([0,0,0],dp)
        endif
        !
        ! convert to fractional coordinates (easier to apply modulo)
        allocate(tau_frac(3,natoms))
        recbas = reciprocal_basis(bas)
        tau_frac = matmul(recbas,tau)
        !
        m = 0
        do i = 1, natoms
            ! apply symmetry operation
            tau_frac_rot(1:3) = matmul(R,tau_frac(1:3,i)) + T(1:3)
            ! reduce rotated+translated point to unit cell
            tau_frac_rot(1:3) = modulo(tau_frac_rot(1:3)+sym_prec,1.0_dp)-sym_prec
            ! check that newly created point matches something already present
            overlap_found = .false.
            check_overlap : do j = 1,natoms
                if (atype(i) .eq. atype(j)) then
                    if ( all(abs(tau_frac_rot(1:3)-tau_frac(1:3,j)).lt.sym_prec) ) then
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
    end function  is_symmetry_valid

    !
    ! functions which operate on uc or similar derived types
    !

    subroutine     load_poscar(uc,iopt_filename,iopts)
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
            symbs=uc%symbs,&
            tau=uc%tau,&
            atype=uc%atype,&
            iopt_filename=filename,&
            iopt_verbosity=opts%verbosity)
        !
    end subroutine load_poscar

    subroutine     reduce_to_primitive(prim,uc,iopts)
        !
        use am_rank_and_sort, only : rank
        !
        implicit none
        !
        class(am_class_unit_cell), intent(inout) :: prim
        type(am_class_unit_cell), intent(in) :: uc
        type(am_class_options), intent(in), optional :: iopts
        type(am_class_options) :: opts
        !
        integer , allocatable :: indices(:)
        real(dp), allocatable :: T(:,:)
        integer :: nTs
        real(dp), allocatable :: tau_prim_frac(:,:)
        real(dp) :: tau_ref(3)
        integer :: i, j, k 
        !
        if (present(iopts)) then 
            opts = iopts
        else
            call opts%defaults
        endif
        !
        if ( opts%verbosity .ge. 1) call am_print_title('Reducing to primitive cell')
        !
        !
        ! 1) get basis translations which could serve as primitive cell vectors
        !
        T = translations_from_basis(bas=uc%bas,tau=uc%tau,atype=uc%atype,iopt_include=trim('prim'),iopt_sym_prec=opts%sym_prec)
        nTs = size(T,2)
        if (opts%verbosity .ge. 1) call am_print("possible primitive lattice translations found",nTs," ... ")
        ! 2) sort them based on magnitude (smallest last)
        allocate(indices(nTs))
        call rank(norm2(T,1),indices)
        T=T(1:3, indices(nTs:1:-1) )
        if (opts%verbosity .ge. 1) call am_print("possible primitive lattice vectors (cart)",transpose(T)," ... ")
        if (opts%verbosity .ge. 1) call am_print("possible primitive lattice vectors (frac)",transpose(matmul(reciprocal_basis(uc%bas),T))," ... ")
        ! 3) select three primitive vectors which yield the smallest cell volume
        ! a) first vector: smallest vector
        i = nTs
        ! b) second vector: smallest pritimive vector which is noncollinear to the first (non-zero cross product)
        do j = nTs,1,-1
            if ( any(cross_product(T(1:3,i),T(1:3,j)) .gt. tiny) ) exit
        enddo
        if (any([i,j].eq.0)) then
            call am_print('ERROR','No pair of non-collinear lattice vectors found',' >>> ')
            stop
        endif
        ! c) third vector: smallest primitive vector which produces a non-zero cell volume
        do k = nTs,1,-1
            if ( abs(primitive_volume(bas=T(1:3,[i,j,k]),skip_volume_check=.true.)) .gt. tiny ) exit
        enddo
        if (any([i,j,k].eq.0)) then
            call am_print('ERROR','No primitive basis found',' >>> ')
            stop
        endif
        !
        prim%bas = T(1:3,[i,j,k])
        if (opts%verbosity .ge. 1) call am_print("volume",primitive_volume(prim%bas)," ... ")
        if (opts%verbosity .ge. 1) call am_print("primitive basis (cart)",prim%bas," ... ")
        if (opts%verbosity .ge. 1) call am_print("primitive basis (frac, in terms of original basis)",matmul(reciprocal_basis(uc%bas),prim%bas)," ... ")
        if (opts%verbosity .ge. 1) call am_print("supercell basis (frac, in terms of primitive basis)",matmul(reciprocal_basis(prim%bas),uc%bas)," ... ")
        !
        ! 4) reduce atoms to primitive cell
        !
        if (opts%verbosity .ge. 1) call am_print('number of atoms in original cell',uc%natoms,' ... ')
        ! a) set first atom at the origin and convert everything to fractional coordinates (easier for modulo)
        allocate(tau_prim_frac(3,uc%natoms))
        tau_ref = uc%tau(1:3,1)
        do i = 1, uc%natoms
            tau_prim_frac(1:3,i) = uc%tau(1:3,i) - tau_ref
        enddo
        ! b) convert to cartesian atomic basis into fractional primitive cell units
        tau_prim_frac = matmul(reciprocal_basis(prim%bas),tau_prim_frac)
        if (opts%verbosity .ge. 1) call am_print('original atomic basis (frac, in terms of primitive)',transpose(tau_prim_frac),' ... ')
        ! c) reduce to primitive cell
        do i = 1,uc%natoms
            tau_prim_frac(1:3,i) = modulo(tau_prim_frac(1:3,i)+opts%sym_prec,1.0_dp)-opts%sym_prec
        enddo
        ! d) get unique values
        tau_prim_frac = unique(tau_prim_frac,opts%sym_prec)
        prim%natoms = size(tau_prim_frac,2)
        if (opts%verbosity .ge. 1) call am_print('number of atoms in reduced cell',prim%natoms,' ... ')
        if (opts%verbosity .ge. 1) call am_print('reduced atomic basis (frac, in terms of primitive)',transpose(tau_prim_frac),' ... ')
        prim%tau = matmul(prim%bas,tau_prim_frac)
        if (opts%verbosity .ge. 1) call am_print('reduced atomic basis (cart)',transpose(prim%tau),' ... ')
        !
        ! 5) transfer atomic species by comparing atomic coordinates
        !
        allocate(prim%atype(prim%natoms))
        do i = 1, prim%natoms
            prim%atype(i) = 0
            search_for_atom : do j = 1, uc%natoms
                if ( all(abs(prim%tau(1:3,i)-uc%tau(1:3,j)).lt.opts%sym_prec) ) then
                    prim%atype(i) = uc%atype(j)
                    exit search_for_atom
                endif
            enddo search_for_atom
            if (prim%atype(i) .eq. 0) then
                call am_print('ERROR','Unable to identify atom in reduced cell',' >>> ')
                stop
            endif
        enddo
        !
        ! 6) transfer other stuff
        !
        allocate(character(500)::prim%symbs(uc%nspecies))
        prim%symbs = uc%symbs
        prim%nspecies = uc%nspecies
    end subroutine reduce_to_primitive

    subroutine     expand_to_supercell(sc,uc,bscfp,iopts)
        !
        ! bscfp = supercell basis in fractional primitive lattice coordinates
        ! for primitive -> conventional cell, use: bscfp = inv(centering_matrix)
        !
        implicit none
        !
        class(am_class_unit_cell), intent(inout) :: sc
        type(am_class_unit_cell) , intent(in) :: uc
        real(dp), intent(in) :: bscfp(3,3)
        real(dp) :: sc_recbas(3,3)
        real(dp) :: tau_wrk(3)
        integer :: i1, i2, i3, j, m
        type(am_class_options), intent(in), optional :: iopts
        type(am_class_options) :: opts
        if (present(iopts)) then 
            opts = iopts
        else
            call opts%defaults
        endif
        !
        !
        if (opts%verbosity .ge. 1) call am_print_title('Expanding to supercell')
        !
        if (opts%verbosity .ge. 1) call am_print('supercell basis (frac, in terms of primitive cell)',bscfp,' ... ')
        !
        if ((mod(det(bscfp)+opts%sym_prec,1.0_dp)-opts%sym_prec).gt.opts%sym_prec) then
            call am_print('ERROR','The determinant of the supercell basis (frac, primitive lattice coordinates) must be an integer.',' >>> ')
            stop
        endif
        !
        sc%natoms = det(bscfp)*uc%natoms
        if (opts%verbosity .ge. 1) call am_print('number of atoms in supercell',sc%natoms,' ... ')
        !
        sc%bas = matmul(uc%bas,bscfp)
        if (opts%verbosity .ge. 1) call am_print('supercell basis (cart)',sc%bas,' ... ')
        !
        allocate(sc%tau(3,sc%natoms))
        allocate(sc%atype(sc%natoms))
        !
        sc_recbas = reciprocal_basis(sc%bas)
        sc%tau = 0.0_dp
        m=0
        map : do i1 = 1, sc%natoms
              do i2 = 1, sc%natoms
              do i3 = 1, sc%natoms
                  do j = 1, uc%natoms
                      !
                      tau_wrk = uc%tau(1:3,j)
                      ! atomic basis in cartesian coordinates
                      tau_wrk = tau_wrk+real(i1,dp)*uc%bas(1:3,1)
                      tau_wrk = tau_wrk+real(i2,dp)*uc%bas(1:3,2)
                      tau_wrk = tau_wrk+real(i3,dp)*uc%bas(1:3,3)
                      ! convert to supercell fractional
                      tau_wrk = matmul(sc_recbas,tau_wrk)
                      ! reduce to supercell primitive
                      tau_wrk = modulo(tau_wrk+opts%sym_prec,1.0_dp) - opts%sym_prec
                      !
                      if (.not.issubset(sc%tau,tau_wrk,opts%sym_prec)) then
                          m = m+1 
                          sc%tau(1:3,m) = tau_wrk
                          sc%atype(m) = uc%atype(j)
                      endif
                      if (m.eq.sc%natoms) then
                          exit map
                      endif
                  enddo
              enddo
              enddo
        enddo map
        ! correct basic rounding error
        where(abs(sc%tau).lt.opts%sym_prec) sc%tau=0.0_dp 
        !
        if (opts%verbosity .ge. 1) call am_print("supercell atomic basis (frac)",transpose(sc%tau),' ... ')
        ! convert to cartesian
        sc%tau = matmul(sc%bas,sc%tau)
        if (opts%verbosity .ge. 1) call am_print("supercell atomic basis (cart)",transpose(sc%tau),' ... ')
        !
    end subroutine expand_to_supercell

    subroutine     primitive_to_conventional(conv,prim,iopts)
        !
        use am_rank_and_sort
        !
        implicit none
        !
        class(am_class_unit_cell), intent(inout) :: conv
        type(am_class_unit_cell), intent(in) :: prim
        real(dp) :: centering(3,3)
        integer :: lattice_code
        type(am_class_options), intent(in), optional :: iopts
        type(am_class_options) :: opts
        if (present(iopts)) then 
            opts = iopts
        else
            call opts%defaults
        endif
        !
        !
        if (opts%verbosity.ge.1) call am_print_title('Transforming to conventional cell') 
        !
        lattice_code = lattice_type(prim%bas)
        !
        if (opts%verbosity.ge.1) then
        select case(lattice_code)
        case(1);  call am_print('lattice','simple cubic cell')
        case(2);  call am_print('lattice','body-centered cubic cell')
        case(3);  call am_print('lattice','face-centered cubic cell')
        case(4);  call am_print('lattice','hexagonal cell')
        case(5);  call am_print('lattice','simple tetragonal cell')
        case(6);  call am_print('lattice','body-centered tetragonal cell')
        case(7);  call am_print('lattice','rhombohedral (trigonal) cell')
        case(8);  call am_print('lattice','simple orthorhombic cell')
        case(9);  call am_print('lattice','body-centered orthorhombic cell')
        case(10); call am_print('lattice','face-centered orthorhombic cell')
        case(11); call am_print('lattice','base centered orthorhombic cell')
        case(12); call am_print('lattice','simple monoclinic cell')
        case(13); call am_print('lattice','base-centered monoclinic cell')
        case(14); call am_print('lattice','triclinic cell')
        case default
            call am_print('ERROR','Unknown lattice.')
        end select
        endif
        !
        centering = lattice_centering(lattice_code)
        !
        if (opts%verbosity.ge.1) call am_print('centering matrix',centering)
        !
        call conv%expand_to_supercell(uc=prim,bscfp=inv(centering),iopts=opts)
        !
    end subroutine primitive_to_conventional

    pure function  space_symmetries_from_basis(uc,iopts) result(seitz)
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
        type(am_class_unit_cell), intent(in) :: uc
        real(dp), allocatable :: seitz(:,:,:)
        integer               :: nTs
        integer               :: nRs
        real(dp), allocatable :: T(:,:)
        real(dp), allocatable :: R(:,:,:)
        real(dp), allocatable :: wrkspace(:,:,:)
        integer :: i, j, m
        !
        type(am_class_options), intent(in), optional :: iopts
        type(am_class_options) :: opts
        if (present(iopts)) then 
            opts = iopts
        else
            call opts%defaults
        endif
        !
        R = determine_lattice_point_symmetries(bas=uc%bas)
        nRs = size(R,3)
        !
        T = translations_from_basis(bas=uc%bas,tau=uc%tau,atype=uc%atype,iopt_sym_prec=opts%sym_prec,iopt_include="zero")
        nTs = size(T,2)
        !
        allocate(wrkspace(4,4,nTs*nRs))
        wrkspace = 0.0_dp
        m=0
        do i = 1, nRs
        do j = 1, nTs
            if ( is_symmetry_valid(bas=uc%bas,tau=uc%tau,atype=uc%atype,&
                & iopt_R=R(1:3,1:3,i),iopt_T=T(1:3,j),iopt_sym_prec=opts%sym_prec)) then
                m = m + 1
                wrkspace(1:3,1:3,m)=R(1:3,1:3,i)
                wrkspace(1:3,4,m)=T(1:3,j)
                wrkspace(4,4,m)=1.0_dp
            endif
        enddo
        enddo
        !
        allocate(seitz(4,4,m))
        seitz = wrkspace(1:4,1:4,1:m)
        !
    end function   space_symmetries_from_basis
   


end module




