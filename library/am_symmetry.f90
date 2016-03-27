module am_symmetry
    !
    use am_constants
    use am_helpers
    use am_unit_cell
    use am_options
    use am_mkl
    !
    implicit none
    !
    private
    !
    public :: am_class_symmetry
    public :: determine_symmetry
    public :: bond_group
    public :: get_kpoint_compatible_symmetries
    public :: point_symmetry_convert_to_cartesian

    type am_class_symmetry 
        integer :: nsyms                  !> number of point symmetries
        real(dp), allocatable :: R(:,:,:) !> symmetry elements (operate on fractional atomic basis)
        real(dp), allocatable :: T(:,:)   !> symmetry elements (operate on fractional atomic basis)
        !
        integer , allocatable :: ps_identifier(:) ! look at decode_pointsymmetry to understand what each integer corresponds to
        integer               :: pg_identifier
        integer , allocatable :: conjugacy_class(:) !> indices assiging each element to a conjugacy class
    contains
        procedure :: point_group
        procedure :: space_group
        procedure :: stdout
        procedure :: name_symmetries
        procedure :: sort
        procedure :: create
        !
        procedure :: tensor_conductivity ! both electrical and thermal
        procedure :: tensor_thermoelectricity
        procedure :: tensor_elasticity
        procedure :: tensor_pizeoelectricity
        !
        procedure :: symmetry_adapted_2nd_order_force_constants
        procedure :: get_action_table
        procedure :: stabilizers
    end type

    interface transform_tensor
        module procedure transform_tensor_first_rank, transform_tensor_second_rank, transform_tensor_third_rank, transform_tensor_fourth_rank
    end interface ! transform_tensor

contains

    !
    ! functions which operate on R(3,3) in fractional coordinates
    !

    pure function  point_symmetry_convert_to_cartesian(R_frac,bas) result(R)
        !
        use am_unit_cell, only : reciprocal_basis
        !
        implicit none
        !
        real(dp), intent(in) :: R_frac(3,3)
        real(dp), intent(in) :: bas(3,3)
        real(dp) :: R(3,3)
        !
        R = matmul(bas,matmul(R_frac,reciprocal_basis(bas)))
        !
    end function   point_symmetry_convert_to_cartesian

    pure function  point_symmetry_convert_to_fractional(R,bas) result(R_frac)
        !
        use am_unit_cell, only : reciprocal_basis
        !
        implicit none
        !
        real(dp), intent(in) :: R(3,3)
        real(dp), intent(in) :: bas(3,3)
        real(dp) :: R_frac(3,3)
        !
        R_frac = matmul(reciprocal_basis(bas),matmul(R,bas))
        !
    end function   point_symmetry_convert_to_fractional

    pure function  point_symmetry_trace(R) result(trace)
        !
        implicit none
        !
        real(dp), intent(in) :: R(3,3)
        real(dp) :: trace
        !
        trace = R(1,1)+R(2,2)+R(3,3)
        !
    end function   point_symmetry_trace

    pure function  point_symmetry_determinant(R) result(determinant)
        !
        implicit none
        !
        real(dp), intent(in) :: R(3,3)
        real(dp) :: determinant
        !
        determinant=R(1,1)*(R(2,2)*R(3,3)-R(2,3)*R(3,2))+&
                  & R(1,2)*(R(2,3)*R(3,1)-R(2,1)*R(3,3))+&
                  & R(1,3)*(R(2,1)*R(3,2)-R(2,2)*R(3,1))
        !
    end function   point_symmetry_determinant

    function       point_symmetry_schoenflies(R) result(ps_identifier)
        !>
        !> Point symmetries in fractional coordinates so that they are nice integers which can be easily classified.
        !>
        !>    element:         e    i  c_2  c_3  c_4  c_6  s_2  s_6  s_4  s_3
        !>    trace:          +3   -3   -1    0   +1   +2   +1    0   -1   -2
        !>    determinant:    +1   -1   +1   +1   +1   +1   -1   -1   -1   -1
        !>
        implicit none
        !
        real(dp), intent(in) :: R(3,3)
        real(dp) :: tr
        real(dp) :: det
        integer :: ps_identifier
        !
        ! get trace and determinant (fractional)
        tr  = point_symmetry_trace(R)
        det = point_symmetry_determinant(R)
        !
        ! The Mathematical Theory of Symmetry in Solids: Representation Theory for
        ! Point Groups and Space Groups. 1 edition. Oxford?: New York: Oxford University
        ! Press, 2010. page 138, table 3.8.
        !
        ps_identifier = 0
        if     ( (abs(tr - 3).lt.tiny) .and. (abs(det - 1).lt.tiny) ) then; ps_identifier = 1  ! 'e'
        elseif ( (abs(tr + 1).lt.tiny) .and. (abs(det - 1).lt.tiny) ) then; ps_identifier = 2  ! 'c_2'
        elseif ( (abs(tr - 0).lt.tiny) .and. (abs(det - 1).lt.tiny) ) then; ps_identifier = 3  ! 'c_3'
        elseif ( (abs(tr - 1).lt.tiny) .and. (abs(det - 1).lt.tiny) ) then; ps_identifier = 4  ! 'c_4'
        elseif ( (abs(tr - 2).lt.tiny) .and. (abs(det - 1).lt.tiny) ) then; ps_identifier = 5  ! 'c_6'
        elseif ( (abs(tr + 3).lt.tiny) .and. (abs(det + 1).lt.tiny) ) then; ps_identifier = 6  ! 'i'
        elseif ( (abs(tr - 1).lt.tiny) .and. (abs(det + 1).lt.tiny) ) then; ps_identifier = 7  ! 's_2'
        elseif ( (abs(tr - 0).lt.tiny) .and. (abs(det + 1).lt.tiny) ) then; ps_identifier = 8  ! 's_6'
        elseif ( (abs(tr + 1).lt.tiny) .and. (abs(det + 1).lt.tiny) ) then; ps_identifier = 9  ! 's_4'
        elseif ( (abs(tr + 2).lt.tiny) .and. (abs(det + 1).lt.tiny) ) then; ps_identifier = 10 ! 's_3'
        endif
        !
        if (ps_identifier .eq. 0) then
            call am_print('ERROR','Unable to identify point symmetry.',' >>> ')
            call am_print('R (frac)',R)
            call am_print('tr (frac)',tr)
            call am_print('det (frac)',det)
            stop
        endif
        !
    end function   point_symmetry_schoenflies

    !
    ! functions which operate on S(4,4) seitz matrix in fractional coordinates
    !
    
    pure function  seitz_multiply(A,B) result(C)
        !
        ! in fractional coordinates rotational part of seitz operator is still unitary; it is not, however, an orthogonal matrix
        ! for which A^-1 = A^T. But that's okay over here. Reduces translational part to within primitive cell.
        !
        implicit none
        !
        real(dp), intent(in) :: A(4,4)
        real(dp), intent(in) :: B(4,4)
        real(dp) :: C(4,4)
        !
        C = matmul(A,B)
        !
        C(1:3,4) = modulo(C(1:3,4)+tiny,1.0_dp)-tiny
        !
    end function   seitz_multiply

    pure function  seitz_make(R,T) result(B)
        !
        implicit none
        !
        real(dp), intent(in), optional :: R(3,3)
        real(dp), intent(in), optional :: T(3)
        real(dp) :: B(4,4)
        !
        B      = 0.0_dp
        B(1,1) = 1.0_dp
        B(2,2) = 1.0_dp
        B(3,3) = 1.0_dp
        B(4,4) = 1.0_dp
        if (present(R)) B(1:3,1:3) = R
        if (present(T)) B(1:3,4) = T
        !
    end function   seitz_make


    !
    ! functions for determining symmetry-adapted second-order force consants
    !

    subroutine     symmetry_adapted_2nd_order_force_constants(sg,uc,opts)
        !
        use am_progress_bar
        !
        class(am_class_symmetry), intent(in) :: sg
        type(am_class_unit_cell), intent(in) :: uc
        type(am_class_options), intent(in) :: opts
        !
        real(dp) :: M(3,3,uc%natoms,uc%natoms)      !> IMPORTANT: M is the original tensor before the permutation (shape of tensor is very important!)
        real(dp) :: RM(3,3,uc%natoms,uc%natoms)     !> rotated/permuted matrix M
        integer  :: nterms                          !> nterms the number of terms
        integer  :: nsyms                           !> nsyms number of symmetry elements
        integer , allocatable :: PM(:,:)            !> PM(uc%natoms,sg%nsyms) permutation map; shows how atoms are permuted by each space symmetry operation
        real(dp), allocatable :: M1d(:)             !> M1d(nterms) an array which is used to build the matrix M
        real(dp), allocatable :: T(:,:)             !> T(nsyms,nterms) permutation matrx containing rows vectors describing how M is rotated by a point symmetry R
        real(dp), allocatable :: A(:,:)             !> A(nsyms*nterms,2*nterms) augmented matrix equation
        real(dp), allocatable :: LHS(:,:)           !> LHS(nterms,nterms) left hand side of augmented matrix equation (should be identity after reducing to row echlon form)
        real(dp), allocatable :: RHS(:,:)           !> RHS(nterms,nterms) right hand side of augmented matrix equation
        logical , allocatable :: is_dependent(:)    !> is_dependent(nterms) logical array which describes which terms depend on others
        logical , allocatable :: is_zero(:)         !> is_zero(nterms) logical array which shows terms equal to zero
        logical , allocatable :: is_independent(:)  !> is_independent(nterms) logical array which shows which terms are independent
        integer :: i, j, k, l, n
        integer :: nequations
        integer :: rank_of_T
        !
        if (opts%verbosity.ge.1) call am_print_title('Determining symmetry-adapted 2nd order force constants')
        !
        PM = permutation_map(sg=sg,uc=uc,opts=opts)
        !
        nsyms = size(sg%R,3)
        nterms = product(shape(M))
        !
        allocate(M1d(nterms))
        allocate(T(nsyms,nterms))
        allocate(A(nsyms*nterms,2*nterms))
        !
        A = 0
        nequations = 0
        !
        ! apply point symmetries
        !
        !$OMP PARALLEL PRIVATE(i,j,M,M1d,T) SHARED(A,uc,sg,PM,nequations)
        !$OMP DO
        do i = 1,nterms
            write(*,*) i/real(nterms)
            !
            M1d = 0
            M1d(i) = 1
            M = reshape(M1d,shape(M))
            do j = 1, nsyms
                !
                ! each point symmetry corresponds to nterms possible equations among terms
                !
                nequations = nequations + nterms
                !
                ! apply transformation
                !
                T(j,1:nterms) = pack( transform_2nd_order_force_constants(&
                    M=M,R=point_symmetry_convert_to_cartesian(R_frac=sg%R(:,:,j),bas=uc%bas),PM=PM(:,j)) ,.true.)
                !
            enddo
            ! save an augmented matrix A
            A([1:nsyms]+nsyms*(i-1),[1:nterms]+nterms*(1-1)) = T(1:nsyms,1:nterms)
            A([1:nsyms]+nsyms*(i-1),[1:nterms]+nterms*(2-1)) = 0.0_dp
            A([1:nsyms]+nsyms*(i-1),         i+nterms*(2-1)) = 1.0_dp
        enddo
        !$OMP END DO
        !$OMP END PARALLEL
        !
        ! reduce to row echlon form
        call rref(A)
        !
        !         ! now apply intrinsic symmetries, which are fundamental properties of the
        !         ! elastic/stress/strain tensors and, therefore, independent of the crystal
        !         ! Nye, J.F. "Physical properties of crystals: their representation by tensors and matrices". p 133 Eq 5 and 6
        !         n = nterms
        !         do l = 1,3
        !         do k = 1,3
        !         do j = 1,3
        !         do i = 1,3
        !             ! cijkl = cjikl ! permute first pair of indicies
        !             nequations = nequations+1
        !             n  = n+1
        !             M  = 0.0_dp
        !             RM = 0.0_dp
        !              M(i,j,l,k) = 1.0_dp
        !             RM(j,i,l,k) = 1.0_dp
        !             A(n,[1:nterms]+nterms*(1-1))= pack( M , .true. ) ! RHS
        !             A(n,[1:nterms]+nterms*(2-1))= pack( RM, .true. ) ! LHS
        !             ! cijkl = cijlk ! permute last  pair of indicies
        !             nequations = nequations+1
        !             n  = n+1
        !             M  = 0.0_dp; 
        !             RM = 0.0_dp
        !              M(i,j,k,l) = 1.0_dp
        !             RM(i,j,l,k) = 1.0_dp
        !             A(n,[1:nterms]+nterms*(1-1))= pack( M , .true. ) ! RHS
        !             A(n,[1:nterms]+nterms*(2-1))= pack( RM, .true. ) ! LHS
        !         enddo
        !         enddo
        !         enddo
        !         enddo
        !         ! reduce to echlon form once again to incorporate the intrinsic symmetries.
        !         call rref(A)
        !
        ! call am_print_sparse('A = [LHS|RHS]',A(1:nterms,:),' ... ')
        ! At this point A is an augmented matrix composed of [ LHS | RHS ]. The LHS
        ! should be the identity matrix, which, together with the RHS, completely 
        ! specifies all relationships between variables. 
        allocate(LHS(nterms,nterms))
        allocate(RHS(nterms,nterms))
        LHS = A(1:nterms,1:nterms)
        RHS = A(1:nterms,[1:nterms]+nterms)
        !
        if ( any(abs(LHS-eye(nterms)).gt.tiny) ) then
            call am_print('ERROR','Unable to reduce matrix to row echlon form.')
            call am_print_sparse('spy(LHS)',LHS)
            stop
        endif
        !
        call am_print('number of symmetry equations',nequations,' ... ')
        !
        call parse_symmetry_equations(LHS=LHS,RHS=RHS,is_zero=is_zero,&
            is_independent=is_independent,is_dependent=is_dependent,verbosity=opts%verbosity)
        !
    end subroutine symmetry_adapted_2nd_order_force_constants

    pure function  transform_2nd_order_force_constants(M,R,PM) result(RM)
        !
        implicit none
        !
        real(dp), intent(in) :: M(:,:,:,:) ! M(3,3,uc%natoms,uc%natoms) 
        real(dp), intent(in) :: R(3,3)  ! rotational part of space symmetry
        integer , intent(in) :: PM(:) ! ! PM(uc%natoms) permutation map contains a table of indices describing atoms space symmetry map other atoms to
        real(dp) :: RM(size(M,1),size(M,2),size(M,3),size(M,4))
        integer  :: i, j, alpha, beta, gamma, delta
        integer  :: natoms
        !
        natoms = size(PM,1)
        !
        ! second order force constants transform like 
        ! p 678 "Symmetry and condensed matter physics" Wooten
        !
        RM = 0
        !
        do i = 1,natoms
        do j = 1,natoms
        do beta  = 1,3
        do alpha = 1,3
            ! gamma -> alpha
            ! delta -> beta
            do gamma = 1,3
            do delta = 1,3
                RM(alpha,beta,PM(i),PM(j)) = RM(alpha,beta,PM(i),PM(j)) + R(alpha,gamma)*R(beta,delta)*M(gamma,delta,i,j)
            enddo
            enddo
        enddo
        enddo
        enddo
        enddo
        !
    end function   transform_2nd_order_force_constants

    !
    ! functions which operate on tensors or make symmetry-adapted tensors
    !

    subroutine     tensor_conductivity(pg,uc,opts)
        !
        class(am_class_symmetry), intent(in) :: pg
        type(am_class_unit_cell), intent(in) :: uc
        type(am_class_options), intent(in) :: opts
        !
        real(dp) :: M(3,3)                          !> IMPORTANT: M is the original tensor before the permutation (shape of tensor is very important!)
        real(dp) :: RM(3,3)                         !> rotated/permuted matrix M
        integer  :: nterms                          !> nterms the number of terms
        integer  :: nsyms                           !> nsyms number of symmetry elements
        integer  :: tensor_rank                     !> tensor_rank rank of tensor matrix M, determined automatically by code
        real(dp), allocatable :: M1d(:)             !> M1d(nterms) an array which is used to build the matrix M
        real(dp), allocatable :: T(:,:)             !> T(nsyms,nterms) permutation matrx containing rows vectors describing how M is rotated by a point symmetry R
        real(dp), allocatable :: A(:,:)             !> A(nsyms*nterms,2*nterms) augmented matrix equation
        real(dp), allocatable :: LHS(:,:)           !> LHS(nterms,nterms) left hand side of augmented matrix equation (should be identity after reducing to row echlon form)
        real(dp), allocatable :: RHS(:,:)           !> RHS(nterms,nterms) right hand side of augmented matrix equation
        logical , allocatable :: is_dependent(:)    !> is_dependent(nterms) logical array which describes which terms depend on others
        logical , allocatable :: is_zero(:)         !> is_zero(nterms) logical array which shows terms equal to zero
        logical , allocatable :: is_independent(:)  !> is_independent(nterms) logical array which shows which terms are independent
        integer :: i, j, n
        integer :: nequations
        !
        if (opts%verbosity.ge.1) call am_print_title('Determining conductivity tensor symmetry')
        !
        tensor_rank = size(shape(M))
        nsyms = size(pg%R,3)
        nterms = 3**tensor_rank
        !
        allocate(M1d(nterms))
        allocate(T(nsyms,nterms))
        allocate(A(nsyms*nterms,2*nterms)) ! augmented matrix
        !
        A = 0
        nequations = 0
        ! <CRYSTAL_SYMMETRIES>
        ! apply point symmetries (each point symmetry corresponds to nterms possible equations
        ! relating the terms to each other)
        do i = 1,nterms
            M1d = 0
            M1d(i) = 1
            M = reshape(M1d,shape(M))
            do j = 1, nsyms
                nequations = nequations + nterms
                ! RM shows how the M has been permuted by the symmetry operation j save the permutation obtained with RM
                ! as a row-vectors in T(nsyms,nterms)  Nye, J.F. "Physical properties of crystals: their representation
                ! by tensors and matrices"
                T(j,1:nterms) = pack( transform_tensor(M=M,R=point_symmetry_convert_to_cartesian(R_frac=pg%R(:,:,j),bas=uc%bas)) ,.true.)
            enddo
            ! save an augmented matrix A
            A([1:nsyms]+nsyms*(i-1),[1:nterms]+nterms*(1-1)) = T(1:nsyms,1:nterms)
            A([1:nsyms]+nsyms*(i-1),[1:nterms]+nterms*(2-1)) = 0.0_dp
            A([1:nsyms]+nsyms*(i-1),         i+nterms*(2-1)) = 1.0_dp
        enddo
        ! reduce to row echlon form to incorporate crystal symmetries
        call rref(A)
        ! </CRYSTAL_SYMMETRIES>
        !
        ! <INTRINSIC_SYMMETRIES>
        ! Now apply intrinsic symmetries, which are fundamental properties of the elastic/stress/strain tensors and,
        ! therefore, independent of the crystal Nye, J.F. "Physical properties of crystals: their representation by
        ! tensors and matrices". p 30. Note: Conductivity tensor is a symmetric second-rank tensor. (see p 22 and 227)
        ! Thermoelectric tensor is an anti-symmetric second-rank tensor. Also see p 194 Eq 4, which proves that the
        ! thermal and electrical conductivities are symmetric tensors explicitly.
        n = nterms
        do j = 1,3
        do i = 1,3
            ! sigma_ij = sigma_ji ! permute indices
            nequations = nequations+1
            n  = n+1
            M  = 0.0_dp
            RM = 0.0_dp
             M(i,j) = 1.0_dp
            RM(j,i) = 1.0_dp
            A(n,[1:nterms]+nterms*(1-1))= pack( M , .true. ) ! RHS
            A(n,[1:nterms]+nterms*(2-1))= pack( RM, .true. ) ! LHS
        enddo
        enddo
        ! reduce to echlon form once again to incorporate the intrinsic symmetries
        call rref(A)
        ! </INTRINSIC_SYMMETRIES>
        !
        ! At this point A is an augmented matrix composed of [ LHS | RHS ]. The LHS should be the identity matrix,
        ! which, together with the RHS, completely specifies all relationships between variables.
        allocate(LHS(nterms,nterms))
        allocate(RHS(nterms,nterms))
        LHS = A(1:nterms,1:nterms)
        RHS = A(1:nterms,[1:nterms]+nterms)
        !
        if ( any(abs(LHS-eye(nterms)).gt.tiny) ) then
            call am_print('ERROR','Unable to reduce matrix to row echlon form.')
            call am_print_sparse('spy(LHS)',LHS)
            stop
        endif
        !
        call am_print('number of symmetry equations',nequations,' ... ')
        !
        call parse_symmetry_equations(LHS=LHS,RHS=RHS,is_zero=is_zero,&
            is_independent=is_independent,is_dependent=is_dependent,verbosity=opts%verbosity)
        !
    end subroutine tensor_conductivity

    subroutine     tensor_thermoelectricity(pg,uc,opts)
        !
        ! Onsager’s Principle requires that the electric resistivity and thermal conductivity tensors be symmetric, but
        ! this does not hold for the Seebeck, Peltier (thermoelectric) coefficients which relate two different flows.
        ! Thus there are, at most, nine coefficients to be determined rather than six.
        !
        ! Newnham "Properties of Materials"
        !
        class(am_class_symmetry), intent(in) :: pg
        type(am_class_unit_cell), intent(in) :: uc
        type(am_class_options), intent(in) :: opts
        !
        real(dp) :: M(3,3)                          !> IMPORTANT: M is the original tensor before the permutation (shape of tensor is very important!)
        integer  :: nterms                          !> nterms the number of terms
        integer  :: nsyms                           !> nsyms number of symmetry elements
        integer  :: tensor_rank                     !> tensor_rank rank of tensor matrix M, determined automatically by code
        real(dp), allocatable :: M1d(:)             !> M1d(nterms) an array which is used to build the matrix M
        real(dp), allocatable :: T(:,:)             !> T(nsyms,nterms) permutation matrx containing rows vectors describing how M is rotated by a point symmetry R
        real(dp), allocatable :: A(:,:)             !> A(nsyms*nterms,2*nterms) augmented matrix equation
        real(dp), allocatable :: LHS(:,:)           !> LHS(nterms,nterms) left hand side of augmented matrix equation (should be identity after reducing to row echlon form)
        real(dp), allocatable :: RHS(:,:)           !> RHS(nterms,nterms) right hand side of augmented matrix equation
        logical , allocatable :: is_dependent(:)    !> is_dependent(nterms) logical array which describes which terms depend on others
        logical , allocatable :: is_zero(:)         !> is_zero(nterms) logical array which shows terms equal to zero
        logical , allocatable :: is_independent(:)  !> is_independent(nterms) logical array which shows which terms are independent
        integer :: i, j
        integer :: nequations
        !
        if (opts%verbosity.ge.1) call am_print_title('Determining thermoelectric tensor symmetry')
        !
        tensor_rank = size(shape(M))
        nsyms = size(pg%R,3)
        nterms = 3**tensor_rank
        !
        allocate(M1d(nterms))
        allocate(T(nsyms,nterms))
        allocate(A(nsyms*nterms,2*nterms)) ! augmented matrix
        !
        A = 0
        nequations = 0
        ! <CRYSTAL_SYMMETRIES>
        ! apply point symmetries (each point symmetry corresponds to nterms possible equations
        ! relating the terms to each other)
        do i = 1,nterms
            M1d = 0
            M1d(i) = 1
            M = reshape(M1d,shape(M))
            do j = 1, nsyms
                nequations = nequations + nterms
                ! RM shows how the M has been permuted by the symmetry operation j save the permutation obtained with RM
                ! as a row-vectors in T(nsyms,nterms)  Nye, J.F. "Physical properties of crystals: their representation
                ! by tensors and matrices"
                T(j,1:nterms) = pack( transform_tensor(M=M,R=point_symmetry_convert_to_cartesian(R_frac=pg%R(:,:,j),bas=uc%bas)) ,.true.)
            enddo
            ! save an augmented matrix A
            A([1:nsyms]+nsyms*(i-1),[1:nterms]+nterms*(1-1)) = T(1:nsyms,1:nterms)
            A([1:nsyms]+nsyms*(i-1),[1:nterms]+nterms*(2-1)) = 0.0_dp
            A([1:nsyms]+nsyms*(i-1),         i+nterms*(2-1)) = 1.0_dp
        enddo
        ! reduce to row echlon form to incorporate crystal symmetries
        call rref(A)
        ! </CRYSTAL_SYMMETRIES>
        !
        ! <INTRINSIC_SYMMETRIES>
        !
        ! THERMOELECTIRCITY TENSOR HAS NO INTRINSIC SYMMETRY. Nye p 227, Eq 47.
        !
        ! Now apply intrinsic symmetries, which are fundamental properties of the elastic/stress/strain tensors and,
        ! therefore, independent of the crystal Nye, J.F. "Physical properties of crystals: their representation by
        ! tensors and matrices". p 30. Note: Conductivity tensor is a symmetric second-rank tensor. (see p 22 and 227)
        ! Also see p 194 Eq 4, which proves that the thermal and electrical conductivities are symmetric tensors
        ! explicitly.
        ! </INTRINSIC_SYMMETRIES>
        !
        ! At this point A is an augmented matrix composed of [ LHS | RHS ]. The LHS should be the identity matrix,
        ! which, together with the RHS, completely specifies all relationships between variables.
        allocate(LHS(nterms,nterms))
        allocate(RHS(nterms,nterms))
        LHS = A(1:nterms,1:nterms)
        RHS = A(1:nterms,[1:nterms]+nterms)
        !
        if ( any(abs(LHS-eye(nterms)).gt.tiny) ) then
            call am_print('ERROR','Unable to reduce matrix to row echlon form.')
            call am_print_sparse('spy(LHS)',LHS)
            stop
        endif
        !
        call am_print('number of symmetry equations',nequations,' ... ')
        !
        call parse_symmetry_equations(LHS=LHS,RHS=RHS,is_zero=is_zero,&
            is_independent=is_independent,is_dependent=is_dependent,verbosity=opts%verbosity)
        !
    end subroutine tensor_thermoelectricity

    subroutine     tensor_pizeoelectricity(pg,uc,opts)
        !
        ! Thus piezoelectricity transforms as a polar third rank tensor. d'_imn = a_ij a_mk a_nl d_jkl. In general there
        ! are 33 = 27 tensor components, but because the stress tensor is symmetric (Xij = Xji), only 18 of the
        ! components are independent. 
        !
        ! Piezoelectric tensor is zero for crystals with inversion symmetry; i.e. the effect only exists in noncentrosymmetric crystals!
        !
        class(am_class_symmetry), intent(in) :: pg
        type(am_class_unit_cell), intent(in) :: uc
        type(am_class_options), intent(in) :: opts
        !
        real(dp) :: M(3,3,3)                        !> IMPORTANT: M is the original tensor before the permutation (shape of tensor is very important!)
        real(dp) :: RM(3,3,3)                       !> rotated/permuted matrix M
        integer  :: nterms                          !> nterms the number of terms
        integer  :: nsyms                           !> nsyms number of symmetry elements
        integer  :: tensor_rank                     !> tensor_rank rank of tensor matrix M, determined automatically by code
        real(dp), allocatable :: M1d(:)             !> M1d(nterms) an array which is used to build the matrix M
        real(dp), allocatable :: T(:,:)             !> T(nsyms,nterms) permutation matrx containing rows vectors describing how M is rotated by a point symmetry R
        real(dp), allocatable :: A(:,:)             !> A(nsyms*nterms,2*nterms) augmented matrix equation
        real(dp), allocatable :: LHS(:,:)           !> LHS(nterms,nterms) left hand side of augmented matrix equation (should be identity after reducing to row echlon form)
        real(dp), allocatable :: RHS(:,:)           !> RHS(nterms,nterms) right hand side of augmented matrix equation
        logical , allocatable :: is_dependent(:)    !> is_dependent(nterms) logical array which describes which terms depend on others
        logical , allocatable :: is_zero(:)         !> is_zero(nterms) logical array which shows terms equal to zero
        logical , allocatable :: is_independent(:)  !> is_independent(nterms) logical array which shows which terms are independent
        integer :: i, j, k, n
        integer :: nequations
        !
        if (opts%verbosity.ge.1) call am_print_title('Determining piezoelectric tensor symmetry')
        !
        tensor_rank = size(shape(M))
        nsyms = size(pg%R,3)
        nterms = 3**tensor_rank
        !
        allocate(M1d(nterms))
        allocate(T(nsyms,nterms))
        allocate(A(nsyms*nterms,2*nterms)) ! augmented matrix
        !
        A = 0
        nequations = 0
        ! <CRYSTAL_SYMMETRIES>
        ! apply point symmetries (each point symmetry corresponds to nterms possible equations
        ! relating the terms to each other)
        do i = 1,nterms
            M1d = 0
            M1d(i) = 1
            M = reshape(M1d,shape(M))
            do j = 1, nsyms
                nequations = nequations + nterms
                ! RM shows how the M has been permuted by the symmetry operation j save the permutation obtained with RM
                ! as a row-vectors in T(nsyms,nterms)  Nye, J.F. "Physical properties of crystals: their representation
                ! by tensors and matrices"
                T(j,1:nterms) = pack( transform_tensor(M=M,R=point_symmetry_convert_to_cartesian(R_frac=pg%R(:,:,j),bas=uc%bas)) ,.true.)
            enddo
            ! save an augmented matrix A
            A([1:nsyms]+nsyms*(i-1),[1:nterms]+nterms*(1-1)) = T(1:nsyms,1:nterms)
            A([1:nsyms]+nsyms*(i-1),[1:nterms]+nterms*(2-1)) = 0.0_dp
            A([1:nsyms]+nsyms*(i-1),         i+nterms*(2-1)) = 1.0_dp
        enddo
        ! reduce to row echlon form to incorporate crystal symmetries
        call rref(A)
        ! </CRYSTAL_SYMMETRIES>
        !
        ! <INTRINSIC_SYMMETRIES>
        ! Now apply intrinsic symmetries, which are fundamental properties, independent of crystal symmetries. Nye, J.F.
        ! "Physical properties of crystals: their representation by tensors and matrices". For the intrinsic
        ! thermoelectric tensor, see p 130. The two last indices, when permuted, leave the tensor invariant because
        ! these indices correspond to the electrical conductivity tensor which is also invariant under such
        ! permutations.
        n = nterms
        do k = 1,3
        do j = 1,3
        do i = 1,3
            ! d_ijk = d_ikj ! permute indices; Nye p 111, Eq 4
            nequations = nequations+1
            n  = n+1
            M  = 0.0_dp
            RM = 0.0_dp
             M(i,j,k) = 1.0_dp
            RM(i,k,j) = 1.0_dp
            A(n,[1:nterms]+nterms*(1-1))= pack( M , .true. ) ! RHS
            A(n,[1:nterms]+nterms*(2-1))= pack( RM, .true. ) ! LHS
        enddo
        enddo
        enddo
        call rref(A)
        ! </INTRINSIC_SYMMETRIES>
        !
        ! At this point A is an augmented matrix composed of [ LHS | RHS ]. The LHS should be the identity matrix,
        ! which, together with the RHS, completely specifies all relationships between variables.
        allocate(LHS(nterms,nterms))
        allocate(RHS(nterms,nterms))
        LHS = A(1:nterms,1:nterms)
        RHS = A(1:nterms,[1:nterms]+nterms)
        !
        if ( any(abs(LHS-eye(nterms)).gt.tiny) ) then
            call am_print('ERROR','Unable to reduce matrix to row echlon form.')
            call am_print_sparse('spy(LHS)',LHS)
            stop
        endif
        !
        call am_print('number of symmetry equations',nequations,' ... ')
        !
        call parse_symmetry_equations(LHS=LHS,RHS=RHS,is_zero=is_zero,&
            is_independent=is_independent,is_dependent=is_dependent,verbosity=opts%verbosity)
        !
    end subroutine tensor_pizeoelectricity

    subroutine     tensor_elasticity(pg,uc,opts)
        !
        class(am_class_symmetry), intent(in) :: pg
        type(am_class_unit_cell), intent(in) :: uc
        type(am_class_options), intent(in) :: opts
        !
        real(dp) :: M(3,3,3,3)                      !> IMPORTANT: M is the original tensor before the permutation (shape of tensor is very important!)
        real(dp) :: RM(3,3,3,3)                     !> rotated/permuted matrix M
        integer  :: nterms                          !> nterms the number of terms
        integer  :: nsyms                           !> nsyms number of symmetry elements
        integer  :: tensor_rank                     !> tensor_rank rank of tensor matrix M, determined automatically by code
        real(dp), allocatable :: M1d(:)             !> M1d(nterms) an array which is used to build the matrix M
        real(dp), allocatable :: T(:,:)             !> T(nsyms,nterms) permutation matrx containing rows vectors describing how M is rotated by a point symmetry R
        real(dp), allocatable :: A(:,:)             !> A(nsyms*nterms,2*nterms) augmented matrix equation
        real(dp), allocatable :: LHS(:,:)           !> LHS(nterms,nterms) left hand side of augmented matrix equation (should be identity after reducing to row echlon form)
        real(dp), allocatable :: RHS(:,:)           !> RHS(nterms,nterms) right hand side of augmented matrix equation
        logical , allocatable :: is_dependent(:)    !> is_dependent(nterms) logical array which describes which terms depend on others
        logical , allocatable :: is_zero(:)         !> is_zero(nterms) logical array which shows terms equal to zero
        logical , allocatable :: is_independent(:)  !> is_independent(nterms) logical array which shows which terms are independent
        integer :: i, j, k, l, n
        integer :: nequations
        !
        if (opts%verbosity.ge.1) call am_print_title('Determining elasticity tensor')
        !
        tensor_rank = size(shape(M))
        nsyms = size(pg%R,3)
        nterms = 3**tensor_rank
        !
        allocate(M1d(nterms))
        allocate(T(nsyms,nterms))
        allocate(A(nsyms*nterms,2*nterms)) ! augmented matrix
        !
        A = 0
        nequations = 0
        !
        ! apply point symmetries (each point symmetry corresponds to nterms possible equations
        ! relating the terms to each other)
        do i = 1,nterms
            M1d = 0
            M1d(i) = 1
            M = reshape(M1d,shape(M))
            do j = 1, nsyms
                nequations = nequations + nterms
                ! RM shows how the M has been permuted by the symmetry operation j
                ! save the permutation obtained with RM as a row-vectors in T(nsyms,nterms) 
                ! T(j,1:nterms) = pack( transform_tensor(M=M,R=pg%R(:,:,j)),.true.)
                ! Nye, J.F. "Physical properties of crystals: their representation by tensors and matrices". p 133 Eq 7
                T(j,1:nterms) = pack( transform_tensor(M=M,R=point_symmetry_convert_to_cartesian(R_frac=pg%R(:,:,j),bas=uc%bas)) ,.true.)
            enddo
            ! save an augmented matrix A
            A([1:nsyms]+nsyms*(i-1),[1:nterms]+nterms*(1-1)) = T(1:nsyms,1:nterms)
            A([1:nsyms]+nsyms*(i-1),[1:nterms]+nterms*(2-1)) = 0.0_dp
            A([1:nsyms]+nsyms*(i-1),         i+nterms*(2-1)) = 1.0_dp
        enddo
        ! reduce to row echlon form
        call rref(A)
        !
        ! now apply intrinsic symmetries, which are fundamental properties of the
        ! elastic/stress/strain tensors and, therefore, independent of the crystal
        ! Nye, J.F. "Physical properties of crystals: their representation by tensors and matrices". p 133 Eq 5 and 6
        n = nterms
        do l = 1,3
        do k = 1,3
        do j = 1,3
        do i = 1,3
            ! cijkl = cjikl ! permute first pair of indicies
            nequations = nequations+1
            n  = n+1
            M  = 0.0_dp
            RM = 0.0_dp
             M(i,j,l,k) = 1.0_dp
            RM(j,i,l,k) = 1.0_dp
            A(n,[1:nterms]+nterms*(1-1))= pack( M , .true. ) ! RHS
            A(n,[1:nterms]+nterms*(2-1))= pack( RM, .true. ) ! LHS
            ! cijkl = cijlk ! permute last  pair of indicies
            nequations = nequations+1
            n  = n+1
            M  = 0.0_dp; 
            RM = 0.0_dp
             M(i,j,k,l) = 1.0_dp
            RM(i,j,l,k) = 1.0_dp
            A(n,[1:nterms]+nterms*(1-1))= pack( M , .true. ) ! RHS
            A(n,[1:nterms]+nterms*(2-1))= pack( RM, .true. ) ! LHS
        enddo
        enddo
        enddo
        enddo
        ! reduce to echlon form once again to incorporate the intrinsic symmetries.
        call rref(A)
        !
        ! call am_print_sparse('A = [LHS|RHS]',A(1:nterms,:),' ... ')
        ! At this point A is an augmented matrix composed of [ LHS | RHS ]. The LHS
        ! should be the identity matrix, which, together with the RHS, completely 
        ! specifies all relationships between variables. 
        allocate(LHS(nterms,nterms))
        allocate(RHS(nterms,nterms))
        LHS = A(1:nterms,1:nterms)
        RHS = A(1:nterms,[1:nterms]+nterms)
        !
        if ( any(abs(LHS-eye(nterms)).gt.tiny) ) then
            call am_print('ERROR','Unable to reduce matrix to row echlon form.')
            call am_print_sparse('spy(LHS)',LHS)
            stop
        endif
        !
        call am_print('number of symmetry equations',nequations,' ... ')
        !
        call parse_symmetry_equations(LHS=LHS,RHS=RHS,is_zero=is_zero,&
            is_independent=is_independent,is_dependent=is_dependent,verbosity=opts%verbosity)
        !
    end subroutine tensor_elasticity

    pure function  transform_tensor_first_rank(M,R)  result(RM)
        !
        implicit none
        !
        real(dp), intent(in) :: M(:)
        real(dp), intent(in) :: R(:,:)
        real(dp) :: RM(size(M,1))
        integer :: i,ip
        !
        ! for a first-rank tensor M, this operation is equivalent to RM = R*M
        !
        RM = 0
        !
        do i = 1,3
        do ip = 1,3
            RM(i) = RM(i) + R(i,ip)*M(ip)
        enddo
        enddo
        !
    end function   transform_tensor_first_rank

    pure function  transform_tensor_second_rank(M,R) result(RM)
        !
        implicit none
        !
        real(dp), intent(in) :: M(:,:)
        real(dp), intent(in) :: R(:,:)
        real(dp) :: RM(size(M,1),size(M,2))
        integer :: i,j,ip,jp
        !
        ! for a second-rank tensor M, this operation is equivalent to RM = R*M*R^T
        ! Eq 13.9 Wootten, p 478 row major order accessed faster
        !
        RM = 0
        !
        do j = 1,3
        do i = 1,3
            do jp = 1,3
            do ip = 1,3
                RM(i,j) = RM(i,j) + R(i,ip)*R(j,jp)*M(ip,jp)
            enddo
            enddo
        enddo
        enddo
        !
    end function   transform_tensor_second_rank

    pure function  transform_tensor_third_rank(M,R)  result(RM)
        !
        implicit none
        !
        real(dp), intent(in) :: M(:,:,:)
        real(dp), intent(in) :: R(3,3)
        real(dp) :: RM(size(M,1),size(M,2),size(M,3))
        integer  :: i,j,k,ip,jp,kp
        !
        ! for a second-rank tensor M, this operation is equivalent to RM = R*M*R^T
        ! Eq 13.9 Wootten, p 478 row major order accessed faster
        !
        RM = 0
        !
        do k = 1,3
        do j = 1,3
        do i = 1,3
            do kp = 1,3
            do jp = 1,3
            do ip = 1,3
                RM(i,j,k) = RM(i,j,k) + R(i,ip)*R(j,jp)*R(k,kp)*M(ip,jp,kp)
            enddo
            enddo
            enddo
        enddo
        enddo
        enddo
        !
    end function   transform_tensor_third_rank

    pure function  transform_tensor_fourth_rank(M,R) result(RM)
        !
        implicit none
        !
        real(dp), intent(in) :: M(:,:,:,:)
        real(dp), intent(in) :: R(3,3)
        real(dp) :: RM(size(M,1),size(M,2),size(M,3),size(M,4))
        integer  :: i,j,k,l,ip,jp,kp,lp
        !
        ! for a second-rank tensor M, this operation is equivalent to RM = R*M*R^T
        ! Eq 13.9 Wootten, p 478 row major order accessed faster
        !
        RM = 0
        !
        do l = 1,3
        do k = 1,3
        do j = 1,3
        do i = 1,3
            do lp = 1,3
            do kp = 1,3
            do jp = 1,3
            do ip = 1,3
                RM(i,j,k,l) = RM(i,j,k,l) + R(i,ip)*R(j,jp)*R(k,kp)*R(l,lp)*M(ip,jp,kp,lp)
            enddo
            enddo
            enddo
            enddo
        enddo
        enddo
        enddo
        enddo
        !
    end function   transform_tensor_fourth_rank

    subroutine     parse_symmetry_equations(LHS,RHS,is_zero,is_independent,is_dependent,verbosity)
        !
        implicit none
        !
        real(dp), intent(in) :: LHS(:,:)
        real(dp), intent(in) :: RHS(:,:)
        logical , intent(out), allocatable :: is_dependent(:)
        logical , intent(out), allocatable :: is_zero(:)
        logical , intent(out), allocatable :: is_independent(:)
        integer :: nterms_zero
        integer :: nterms_independent
        integer :: nterms_dependent
        integer :: nterms
        integer :: verbosity
        integer :: i, j 
        !
        nterms = size(LHS,1)
        allocate(is_dependent(nterms))
        allocate(is_zero(nterms))
        allocate(is_independent(nterms))
        !
        if (verbosity.ge.1) call am_print('total number of terms',nterms,' ... ')
        !
        ! null terms (equal zero)
        !
        is_zero = (all(abs(RHS).lt.tiny,2))
        nterms_zero = count(all(abs(RHS).lt.tiny,2))
        if (verbosity.ge.1) call am_print('number of terms equal to zero',nterms_zero,' ... ')
        !
        ! independent terms (equal themselves and nothing else)
        !
        is_independent = (all(abs(RHS-LHS).lt.tiny,2))
        nterms_independent = count(is_independent)
        if (verbosity.ge.1) call am_print('number of independent terms',nterms_independent,' ... ')
        !
        ! dependent terms (can be written via independent terms)
        !
        is_dependent = (any(abs(RHS).gt.tiny,2))
        is_dependent = (is_dependent.and..not.is_independent)
        nterms_dependent = count(is_dependent)
        if (verbosity.ge.1) call am_print('number of dependent terms',nterms_dependent,' ... ')
        !
        if (nterms_zero+nterms_independent+nterms_dependent.ne.nterms) then
            call am_print('ERROR','The number of terms which are independent, null, and dependent do not add to the total number of terms.',' >>> ')
            stop
        endif
        !
        ! write symmetry equations to stdout
        !
        if (verbosity.ge.1) then
            if (nterms_zero.ne.nterms) then
                !
                write(*,'(a5,a)') ' ... ', 'irreducible symmetry equations'
                !
                ! write the independentterms (equal only to themselves)
                !
                do i = 1,nterms
                    if (is_independent(i)) then
                        write(*,'(5x,a1,a,a1,a1,a)',advance='no') 'a',trim(int2char(i)),'=', 'a', trim(int2char(i))
                        write(*,*)
                    endif
                enddo
                ! write the interrelated terms
                do i = 1,nterms
                    if (is_dependent(i)) then
                        write(*,'(5x,a1,a,a1)',advance='no') 'a',trim(int2char(i)),'='
                        do j = 1,nterms
                            if (abs(RHS(i,j)).gt.tiny) then
                                write(*,'(a,a,a)',advance='no') trim(dbl2charSP(RHS(i,j),5)), '*a', trim(int2char(j))
                            endif
                        enddo
                        write(*,*)
                    endif
                enddo
                ! !
                ! ! finally, write the terms that are equal to zero.
                ! !
                ! do i = 1,nterms
                !     if (is_zero(i)) then
                !         write(*,'(5x,a1,a,a)') 'a',trim(int2char(i)),'=+0.00'
                !     endif
                ! enddo
            endif
        endif
        !
    end subroutine parse_symmetry_equations

    !
    ! schoenflies decode
    !

    function       decode_pointsymmetry(ps_identifier) result(schoenflies)
        !
        implicit none
        !
        integer, intent(in) :: ps_identifier
        character(len=string_length_schoenflies) :: schoenflies
        !
        select case (ps_identifier)
        case(1);  schoenflies = 'e'
        case(2);  schoenflies = 'c_2'
        case(3);  schoenflies = 'c_3'
        case(4);  schoenflies = 'c_4'
        case(5);  schoenflies = 'c_6'
        case(6);  schoenflies = 'i'
        case(7);  schoenflies = 's_2'
        case(8);  schoenflies = 's_6'
        case(9);  schoenflies = 's_4'
        case(10); schoenflies = 's_3'
        case default
            call am_print('ERROR','Schoenflies code unknown.',' >>> ')
            stop
        end select
    end function   decode_pointsymmetry

    function       decode_pointgroup(pg_identifier) result(point_group_name)
        !
        implicit none
        !
        integer, intent(in) :: pg_identifier
        character(string_length_schoenflies) :: point_group_name
        !
        select case (pg_identifier)
                case(1); point_group_name = 'c_1'
                case(2); point_group_name = 's_2'
                case(3); point_group_name = 'c_2'
                case(4); point_group_name = 'c_1h'
                case(5); point_group_name = 'c_2h'
                case(6); point_group_name = 'd_2'
                case(7); point_group_name = 'c_2v'
                case(8); point_group_name = 'd_2h'
                case(9); point_group_name = 'c_3'
                case(10); point_group_name = 's_6'
                case(11); point_group_name = 'd_3'
                case(12); point_group_name = 'c_3v'
                case(13); point_group_name = 'd_3d'
                case(14); point_group_name = 'c_4'
                case(15); point_group_name = 's_4'
                case(16); point_group_name = 'c_4h'
                case(17); point_group_name = 'd_4'
                case(18); point_group_name = 'c_4v'
                case(19); point_group_name = 'd_2d'
                case(20); point_group_name = 'd_4h'
                case(21); point_group_name = 'c_6'
                case(22); point_group_name = 'c_3h'
                case(23); point_group_name = 'c_6h'
                case(24); point_group_name = 'd_6'
                case(25); point_group_name = 'c_6v'
                case(26); point_group_name = 'd_3h'
                case(27); point_group_name = 'd_6h'
                case(28); point_group_name = 't'
                case(29); point_group_name = 't_h'
                case(30); point_group_name = 'o'
                case(31); point_group_name = 't_d'
                case(32); point_group_name = 'o_h'
            case default
                call am_print('ERROR','Schoenflies code for point-group unknown.',' >>> ')
                call am_print('pg_identifier',pg_identifier)
                stop
            end select
    end function   decode_pointgroup

    function       point_group_schoenflies(ps_identifier) result(pg_identifier)
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
        !>   pg_identifier --> point_group_name
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
        integer, intent(in) :: ps_identifier(:)
        integer :: nsyms
        integer :: ni, nc2, nc3, nc4, nc6, ns2, ns6, ns4, ns3, ne
        integer :: pg_identifier
        !
        pg_identifier = 0
        !
        ! trivial cases first
        !
        nsyms = size(ps_identifier)
        if     (nsyms .eq. 1 ) then; pg_identifier=1;  return
        elseif (nsyms .eq. 48) then; pg_identifier=32; return
        elseif (nsyms .eq. 16) then; pg_identifier=20; return
        elseif (nsyms .eq. 3 ) then; pg_identifier=9;  return
        endif
        !
        ni=0; nc2=0; nc3=0; nc4=0; nc6=0; ns2=0; ns6=0; ns4=0; ns3=0
        ne  = count(ps_identifier.eq.1)  ! 'e' identity
        nc2 = count(ps_identifier.eq.2)  ! 'c_2'
        nc3 = count(ps_identifier.eq.3)  ! 'c_3'
        nc4 = count(ps_identifier.eq.4)  ! 'c_4'
        nc6 = count(ps_identifier.eq.5)  ! 'c_6'
        ni  = count(ps_identifier.eq.6)  ! 'i' inversion
        ns2 = count(ps_identifier.eq.7)  ! 's_2'
        ns6 = count(ps_identifier.eq.8)  ! 's_6'
        ns4 = count(ps_identifier.eq.9)  ! 's_4'
        ns3 = count(ps_identifier.eq.10) ! 's_3'
        !
        if (ni.gt.1) then
            call am_print('ERROR','More than one inversion operator found.',' >>> ')
            stop
        endif
        !
        if (ne.gt.1) then
            call am_print('ERROR','More than one identity operator found.',' >>> ')
        elseif (ne.eq.0) then
            call am_print('ERROR','No identity operator found.',' >>> ')
            stop
        endif
        !
        if (nsyms.eq.2)   then
            if (ni.eq.1)  then; pg_identifier=2;  return; endif
            if (nc2.eq.1) then; pg_identifier=3;  return; endif
            if (ns2.eq.1) then; pg_identifier=4;  return; endif
        endif
        if (nsyms.eq.4)   then
            if (ni.eq.1)  then; pg_identifier=5;  return; endif
            if (nc2.eq.3) then; pg_identifier=6;  return; endif
            if (ns2.eq.2) then; pg_identifier=7;  return; endif
            if (nc4.eq.1) then; pg_identifier=14; return; endif
            if (ns4.eq.2) then; pg_identifier=15; return; endif
        endif
        if (nsyms.eq.6)   then
            if (ni.eq.1)  then; pg_identifier=10; return; endif
            if (nc2.eq.3) then; pg_identifier=11; return; endif
            if (ns2.eq.3) then; pg_identifier=12; return; endif
            if (nc2.eq.1) then; pg_identifier=21; return; endif
            if (ns2.eq.1) then; pg_identifier=22; return; endif
        endif
        if (nsyms.eq.8)  then
            if (ns2.eq.3) then; pg_identifier=8;  return; endif
            if (ns2.eq.1) then; pg_identifier=16; return; endif
            if (ns2.eq.0) then; pg_identifier=17; return; endif
            if (ns2.eq.4) then; pg_identifier=18; return; endif
            if (ns2.eq.2) then; pg_identifier=19; return; endif
        endif
        if (nsyms.eq.12)  then
            if (ns2.eq.3) then; pg_identifier=13; return; endif
            if (ns2.eq.1) then; pg_identifier=23; return; endif
            if (nc2.eq.7) then; pg_identifier=24; return; endif
            if (ns2.eq.6) then; pg_identifier=25; return; endif
            if (ns2.eq.4) then; pg_identifier=26; return; endif
            if (nc3.eq.8) then; pg_identifier=28; return; endif
        endif
        if (nsyms.eq.24)  then
            if (nc6.eq.2) then; pg_identifier=27; return; endif
            if (ni.eq.1)  then; pg_identifier=29; return; endif
            if (nc4.eq.6) then; pg_identifier=30; return; endif
            if (ns4.eq.6) then; pg_identifier=31; return; endif
            endif
        ! if it makes it this far, it means nothing matches. return an error immediately. 
        call am_print('ERROR','Unable to identify point group', ' >>> ')
            write(*,'(" ... ",a)') 'number of symmetry operators'
            write(*,'(5x,a5,i5)') 'e  ', ne
            write(*,'(5x,a5,i5)') 'c_2', nc2
            write(*,'(5x,a5,i5)') 'c_3', nc3
            write(*,'(5x,a5,i5)') 'c_4', nc4
            write(*,'(5x,a5,i5)') 'c_6', nc6
            write(*,'(5x,a5,i5)') 'i  ', ni
            write(*,'(5x,a5,i5)') 's_2', ns2
            write(*,'(5x,a5,i5)') 's_6', ns6
            write(*,'(5x,a5,i5)') 's_4', ns4
            write(*,'(5x,a5,i5)') 's_3', ns3
            stop
    end function   point_group_schoenflies

    !
    ! functions which operate on sg
    !

    subroutine     determine_symmetry(uc,sg,pg,opts)
        !
        implicit none
        !
        type(am_class_symmetry), intent(out) :: sg
        type(am_class_symmetry), intent(out) :: pg
        type(am_class_unit_cell), intent(in) :: uc
        type(am_class_options)  , intent(in) :: opts
        !
        if (opts%verbosity.ge.1) call am_print_title('Analyzing symmetry')
        !
        call sg%space_group(uc=uc,opts=opts)
        !
        call sg%get_action_table(uc=uc,iopt_fname='outfile.action_space_group',opts=opts)
        !
        call pg%point_group(sg=sg,uc=uc,opts=opts)
        !
        call pg%get_action_table(uc=uc,iopt_fname='outfile.action_point_group',opts=opts)
        !
    end subroutine determine_symmetry

    subroutine     space_group(sg,uc,opts)
        !
        use am_rank_and_sort
        !
        implicit none
        !
        class(am_class_symmetry), intent(inout) :: sg
        type(am_class_unit_cell), intent(in) :: uc ! space groups are tabulated in the litearture for conventional cells, but primitive or arbitrary cell works just as wlel.
        type(am_class_options)  , intent(in) :: opts
        type(am_class_options) :: notalk
        real(dp), allocatable  :: rep(:,:,:)
        real(dp), allocatable  :: seitz(:,:,:)
        integer :: i
        !
        notalk=opts
        notalk%verbosity = 0
        !
        if (opts%verbosity.ge.1) call am_print_title('Determining space group symmetries')
        !
        ! determine space symmetries from atomic basis
        seitz = space_symmetries_from_basis(uc=uc,opts=opts)
        if (opts%verbosity.ge.1) call am_print('number of space group symmetries found',size(seitz,3),' ... ') 
        !
        ! create space group instance
        call sg%create(R=seitz(1:3,1:3,:),T=seitz(1:3,4,:))
        !
        ! determine conjugacy classes for representation
        sg%conjugacy_class = get_conjugacy_classes(rep=seitz,opts='seitz',verbosity=opts%verbosity)
        !
        ! sort space group symmetries based on parameters
        call sg%sort(sort_parameter=sg%T(1,:),iopt_direction='ascend')
        call sg%sort(sort_parameter=sg%T(2,:),iopt_direction='ascend')
        call sg%sort(sort_parameter=sg%T(3,:),iopt_direction='ascend')
        call sg%sort(sort_parameter=real(sg%conjugacy_class,dp),iopt_direction='ascend')
        !
        ! name the (im-)proper part of space symmetry
        call sg%name_symmetries(opts=notalk)
        !
        ! write to stdout and to file
        if (opts%verbosity.ge.1) call sg%stdout(iopt_uc=uc)
        call sg%stdout(iopt_uc=uc,iopt_filename=trim('outfile.spacegroup'))
        !
    end subroutine space_group

    subroutine     point_group(pg,uc,sg,opts)
        !
        ! requires only sg%R
        !
        use am_rank_and_sort
        !
        implicit none
        !
        class(am_class_symmetry), intent(inout) :: pg
        type(am_class_symmetry) , intent(in) :: sg
        type(am_class_unit_cell), intent(in) :: uc
        type(am_class_options), intent(in) :: opts
        integer , allocatable :: sorted_indices(:)
        real(dp), allocatable :: sort_parameter(:)
        type(am_class_options) :: notalk
        integer :: i
        !
        notalk=opts
        notalk%verbosity = 0
        !
        if (opts%verbosity.ge.1) call am_print_title('Determining point group symmetries')
        !
        ! create point group instance
        call pg%create(R=unique(sg%R))
        if (opts%verbosity.ge.1) call am_print('number of point symmetries',pg%nsyms,' ... ')
        !
        ! determine conjugacy classes for representation
        pg%conjugacy_class = get_conjugacy_classes(rep=pg%R,verbosity=0)
        !
        ! sort point group symmetries based on parameters
        allocate(sort_parameter(pg%nsyms))
        !
        do i = 1, pg%nsyms; sort_parameter(i) = point_symmetry_trace(pg%R(:,:,i)); enddo
        call pg%sort(sort_parameter=sort_parameter,iopt_direction='ascend')
        !
        do i = 1, pg%nsyms; sort_parameter(i) = point_symmetry_determinant(pg%R(:,:,i)); enddo
        call pg%sort(sort_parameter=sort_parameter,iopt_direction='ascend')
        !
        call pg%sort(sort_parameter=real(pg%conjugacy_class,dp),iopt_direction='ascend')
        !
        ! name point symmetries
        call pg%name_symmetries(opts=opts)
        !
        ! name point group
        pg%pg_identifier = point_group_schoenflies(pg%ps_identifier)
        if (opts%verbosity.ge.1) call am_print('point group',decode_pointgroup(pg%pg_identifier),' ... ')
        !
        ! write to stdout and to file
        if (opts%verbosity.ge.1) call pg%stdout(iopt_uc=uc)
        if (opts%verbosity.ge.1) call pg%stdout(iopt_uc=uc,iopt_filename=trim('outfile.pointgroup'))
        !
    end subroutine point_group

    subroutine     bond_group(sg,uc,opts)
        !
        use am_rank_and_sort, only : rank
        !
        implicit none
        !
        type(am_class_symmetry) , intent(in) :: sg ! space group
        type(am_class_unit_cell), intent(in) :: uc ! unit cell
        type(am_class_options)  , intent(in) :: opts
        type(am_class_symmetry) :: bg ! bond group
        type(am_class_symmetry), allocatable :: pg(:,:) ! point group of bond
        type(am_class_options) :: notalk
        integer  :: i, j, k
        real(dp) :: v(3)
        integer, allocatable :: indices(:)
        integer, allocatable :: pair(:,:)
        real(dp),allocatable :: paird(:)
        integer :: npairs
        !
        notalk = opts
        notalk%verbosity = 0
        !
        if (opts%verbosity.ge.1) call am_print_title('Determining bond symmetries')
        !
        npairs=floor(uc%natoms/2.0_dp)*(uc%natoms+1)+modulo(uc%natoms,2)*nint(uc%natoms/2.0_dp)
        if (opts%verbosity.ge.1) call am_print('number of atom pairs', npairs, ' ... ')
        !
        ! sort atom pairs based on distance
        !
        allocate(paird(npairs))
        allocate(pair(2,npairs))
        k=0
        do i = 1, uc%natoms
        do j = 1, i
            k=k+1
            ! get pair 
            pair(1,k) = i
            pair(2,k) = j
            ! bond vector
            v = uc%tau(:,i)-uc%tau(:,j)
            ! get distance between atoms (length of bond)
            paird(k) = norm2(matmul(uc%bas,v))
        enddo
        enddo
        !
        if (k.ne.npairs) then
            call am_print('ERROR','Inconsistent number of atom pairs. Something is wrong internally.',' >>> ')
            call am_print('number of atom pairs', k, ' ... ')
            stop
        endif
        !
        !
        !
        allocate(indices(npairs))
        call rank(paird,indices)
        pair(:,:) = pair(:,indices)
        paird(:)  = paird(indices)
        !
        allocate(pg(uc%natoms,uc%natoms))
        !
        if (opts%verbosity.ge.1) write(*,'(5x,2a5,3a10,a10,a10,a10)') 'i', 'j', 'v_x','v_y','v_z', '|v_frac|', '|v_cart|', 'group'
        do k = 1, npairs
            i = pair(1,k)
            j = pair(2,k)
            ! bond vector
            v = uc%tau(:,i)-uc%tau(:,j)
            ! get stabilizer of bond (group), symmetries which leave bond invariant
            call bg%stabilizers(sg=sg,v=v,opts=notalk)
            ! determine symmetry of bond
            call pg(i,j)%point_group(uc=uc,sg=bg,opts=notalk)
            !
            if (opts%verbosity.ge.1) then 
                write(*,'(5x)',advance='no') 
                write(*,'(2a5)',advance='no') trim(uc%symbs(uc%atype(i))), trim(uc%symbs(uc%atype(j)))
                write(*,'(3f10.2)',advance='no') v
                write(*,'(f10.2)',advance='no') norm2(v)
                write(*,'(f10.2)',advance='no') norm2(matmul(uc%bas,v))
                write(*,'(a10)',advance='no') trim(decode_pointgroup(pg(i,j)%pg_identifier))
                write(*,*)
            endif
            !
        enddo
        !
    end subroutine bond_group

    subroutine     stabilizers(bg,sg,v,opts)
        !
        implicit none
        !
        class(am_class_symmetry), intent(inout) :: bg ! bond group
        type(am_class_symmetry) , intent(in) :: sg ! space group
        type(am_class_options)  , intent(in) :: opts
        real(dp), intent(in) :: v(3) ! bond vector
        logical, allocatable :: mask(:)
        integer :: i, j
        !
        allocate(mask(sg%nsyms))
        mask = .false.
        !
        do i = 1, sg%nsyms
            ! take symmetry which leaves vector invariant (not even modulo a lattice vector). effectively, the
            ! symmetries which is able to map a point onto itself. ignoring translations here: is this correct?
            ! probablly not... need to confirm... one of the problems with the translations is that they are, for example, all
            ! positive, so while they move one atom up, the cannever move one atom back...
            if (is_symmetry_valid(tau=reshape(v,[3,1]),iopt_R=sg%R(:,:,i),&
                iopt_exact=.true.,iopt_sym_prec=opts%sym_prec)) then
                mask(i) = .true.
            endif
        enddo
        !
        bg%nsyms = count(mask)
        if (allocated(bg%R)) deallocate(bg%R)
        if (allocated(bg%T)) deallocate(bg%T)
        allocate(bg%R(3,3,bg%nsyms))
        allocate(bg%T(3,bg%nsyms))
        !
        j=0
        do i = 1, sg%nsyms
            if (mask(i)) then
                j=j+1
                bg%R(:,:,j)=sg%R(:,:,i)
                bg%T(:,j)=sg%T(:,i)
            endif
        enddo
        !
    end subroutine stabilizers

    subroutine     stdout(sg,iopt_uc,iopt_filename)
        !
        implicit none
        !
        class(am_class_symmetry), intent(in) ::  sg
        type(am_class_unit_cell), intent(in), optional :: iopt_uc
        character(*), intent(in), optional :: iopt_filename
        integer :: width
        real(dp) :: R(3,3), T(3)
        integer :: fid
        integer :: i,j,k
        character(15) :: fmt5, fmt4, fmt6, fmt1, fmt2, fmt3
        !
        if ( present(iopt_filename) ) then
            ! write to file
            fid = 1
            open(unit=fid,file=trim(iopt_filename),status="replace",action='write')
        else
            ! stdout
            fid = 6
        endif
        !
        ! write standard out to file
        !
        fmt5=   "(a4)"
        fmt4=   "(a4,3a7)"
        fmt6=       "(a7)"
        fmt1="(5x,i3,9a7)"
        fmt2="(5x,a3,9a7)"
        fmt3="(5x,a)"
        width = 9*7+3
        if (allocated(sg%T)) width=width + 3*7+4
        !
        !
        ! fractional
        !
        write(fid,fmt3,advance='no') centertitle('fractional',width)
        write(fid,*)
        write(fid,fmt3,advance='no') repeat('-',width)
        write(fid,*)
        !
        write(fid,fmt2,advance='no') '#', 'R11', 'R21', 'R31', 'R12', 'R22', 'R32', 'R13', 'R23', 'R33'
        if (allocated(sg%T))  write(fid,fmt4,advance='no') '  | ', 'T1', 'T2', 'T3'
        write(fid,*)
        !
        do i = 1, sg%nsyms
            R = sg%R(:,:,i)
            write(fid,fmt1,advance='no') i, (( trim(print_pretty(R(j,k))), j = 1,3) , k = 1,3)
            if (allocated(sg%T)) then
                T = sg%T(:,i)
                write(fid,fmt5,advance='no') '  | '
                do j = 1,3
                    write(fid,fmt6,advance='no') trim(print_pretty(T(j)))
                enddo
            endif
            write(fid,*)
            !
        enddo
        !
        ! cartesian
        !
        if (present(iopt_uc)) then
        !
        write(fid,fmt3,advance='no') centertitle('cartesian',width)
        write(fid,*)
        write(fid,fmt3,advance='no') repeat('-',width)
        write(fid,*)
        !
        write(fid,fmt2,advance='no') '#', 'R11', 'R21', 'R31', 'R12', 'R22', 'R32', 'R13', 'R23', 'R33'
        if (allocated(sg%T))  write(fid,fmt4,advance='no') '  | ', 'T1', 'T2', 'T3'
        write(fid,*)
        !
        do i = 1, sg%nsyms
            R = point_symmetry_convert_to_cartesian(R_frac=sg%R(:,:,i),bas=iopt_uc%bas)
            write(fid,fmt1,advance='no') i, (( trim(print_pretty(R(j,k))), j = 1,3) , k = 1,3)
            if (allocated(sg%T)) then
                T = matmul(iopt_uc%bas,sg%T(:,i))
                write(fid,fmt5,advance='no') '  | '
                do j = 1,3
                    write(fid,fmt6,advance='no') trim(print_pretty(T(j)))
                enddo
            endif
            write(fid,*)
        enddo
        !
        endif
        !
        ! close file
        !
        if ( present(iopt_filename) ) close(fid)
        !
    end subroutine stdout

    subroutine     name_symmetries(sg,opts)
        !
        implicit none
        !
        class(am_class_symmetry) , intent(inout) :: sg
        type(am_class_options), intent(in) :: opts
        integer :: i
        !
        if (allocated(sg%ps_identifier)) deallocate(sg%ps_identifier)
        allocate(sg%ps_identifier(sg%nsyms))
        !
        do i = 1, sg%nsyms
            sg%ps_identifier(i) = point_symmetry_schoenflies(sg%R(:,:,i))
        enddo
        !
        if (opts%verbosity.ge.1) then
            call am_print('number of e   point symmetries',count(sg%ps_identifier.eq.1),' ... ')
            call am_print('number of c_2 point symmetries',count(sg%ps_identifier.eq.2),' ... ')
            call am_print('number of c_3 point symmetries',count(sg%ps_identifier.eq.3),' ... ')
            call am_print('number of c_4 point symmetries',count(sg%ps_identifier.eq.4),' ... ')
            call am_print('number of c_6 point symmetries',count(sg%ps_identifier.eq.5),' ... ')
            call am_print('number of i   point symmetries',count(sg%ps_identifier.eq.6),' ... ')
            call am_print('number of s_2 point symmetries',count(sg%ps_identifier.eq.7),' ... ')
            call am_print('number of s_6 point symmetries',count(sg%ps_identifier.eq.8),' ... ')
            call am_print('number of s_4 point symmetries',count(sg%ps_identifier.eq.9),' ... ')
            call am_print('number of s_3 point symmetries',count(sg%ps_identifier.eq.10),' ... ')
        endif
        !
    end subroutine name_symmetries

    subroutine     create(sg,R,T)
        !
        implicit none
        !
        class(am_class_symmetry), intent(inout) :: sg
        real(dp), intent(in), optional :: R(:,:,:)
        real(dp), intent(in), optional :: T(:,:)
        integer  :: i
        !
        if (present(T).and.present(R)) then
        if (size(R,3).ne.size(T,2)) then
            call am_print('ERROR','R and T dimensions do not match.')
            call am_print('shape(R)',shape(R))
            call am_print('shape(T)',shape(T))
            stop
        endif
        endif
        !
        if (present(R)) then
            sg%nsyms = size(R,3)
            if (allocated(sg%R)) deallocate(sg%R)
            allocate(sg%R(3,3,sg%nsyms))
            do i = 1, sg%nsyms
                sg%R(:,:,i) = R(:,:,i)
            enddo
        endif
        !
        if (present(T)) then
            sg%nsyms = size(T,2)
            if (allocated(sg%T)) deallocate(sg%T)
            allocate(sg%T(3,sg%nsyms))
            do i = 1, sg%nsyms
                sg%T(:,i) = T(:,i)
            enddo
        endif
        !
    end subroutine create

    !
    ! group theoretical procedures which operate on sg
    !

    function       permutation_rep(sg,uc,opts) result(P)
       !
       ! find permutation representation; i.e. which atoms are connected by space symmetry oprations R, T.
       !
       implicit none
       !
       type(am_class_symmetry) , intent(in) :: sg
       type(am_class_unit_cell), intent(in) :: uc
       type(am_class_options)  , intent(in) :: opts
       integer, allocatable :: P(:,:,:)
       real(dp) :: tau_rot(3)
       integer :: i,j,k
       logical :: found
       !
       allocate(P(uc%natoms,uc%natoms,sg%nsyms))
       P = 0
       !
       do i = 1, sg%nsyms
           ! determine the permutations of atomic indicies which results from each space symmetry operation
           do j = 1,uc%natoms
                found = .false.
                ! apply rotational component
                tau_rot = matmul(sg%R(:,:,i),uc%tau(:,j))
                ! apply translational component
                if (allocated(sg%T)) tau_rot = tau_rot + sg%T(:,i)
                ! reduce rotated+translated point to unit cell
                tau_rot = modulo(tau_rot+opts%sym_prec,1.0_dp)-opts%sym_prec
                ! find matching atom
                do k = 1, uc%natoms
                    if (all(abs(tau_rot-uc%tau(:,k)).lt.opts%sym_prec)) then
                    P(j,k,i) = 1
                    found = .true.
                    exit ! break loop
                    endif
                enddo
                if (found.eq..false.) then
                    call am_print('ERROR','Unable to find matching atom.',' >>> ')
                    call am_print('tau (all atoms)',transpose(uc%tau))
                    call am_print('tau',uc%tau(:,j))
                    call am_print('R',sg%R(:,:,i))
                    call am_print('T',sg%T(:,i))
                    call am_print('tau_rot',tau_rot)
                    stop
                endif
           enddo
       enddo
       !
       ! check that each column and row of the P sums to 1; i.e. P(:,:,i) is orthonormal for all i.
       !
       do i = 1,sg%nsyms
       do j = 1,uc%natoms
           !
           if (sum(P(:,j,i)).ne.1) then
               call am_print('ERROR','Permutation matrix has a column which does not sum to 1.')
               call am_print('i',i)
               call am_print_sparse('spy(P_i)',P(:,:,i))
               call am_print('P',P(:,:,i))
               stop
           endif
           !
           if (sum(P(j,:,i)).ne.1) then
               call am_print('ERROR','Permutation matrix has a row which does not sum to 1.')
               call am_print('i',i)
               call am_print_sparse('spy(P_i)',P(:,:,i))
               call am_print('P',P(:,:,i))
               stop
           endif
           !
       enddo
       enddo
       !
    end function   permutation_rep

    function       permutation_map(sg,uc,opts) result(PM)
        !
        ! PM(uc%natoms,sg%nsyms) permutation map; shows how atoms are permuted by each space symmetry operation
        !
        implicit none
        !
        type(am_class_symmetry) , intent(in) :: sg
        type(am_class_unit_cell), intent(in) :: uc
        type(am_class_options)  , intent(in) :: opts
        integer, allocatable :: P(:,:,:)
        integer, allocatable :: PM(:,:)
        integer :: i
        !
        P = permutation_rep(sg=sg,uc=uc,opts=opts)
        !
        allocate(PM(uc%natoms,sg%nsyms))
        PM = 0
        do i = 1,sg%nsyms
            PM(:,i) = matmul(P(:,:,i),[1:uc%natoms])
        enddo
        !
    end function   permutation_map

    subroutine     get_action_table(sg,uc,iopt_fname,opts)
        !
        implicit none
        !
        class(am_class_symmetry), intent(inout) :: sg
        type(am_class_unit_cell), intent(in) :: uc
        type(am_class_options)  , intent(in) :: opts
        character(*), optional  , intent(in) :: iopt_fname
        integer, allocatable :: P(:,:,:)
        integer, allocatable :: PM(:,:)
        integer :: fid
        character(10) :: buffer
        character(3) :: xyz
        integer :: i,j,k
        real(dp) :: R(3,3)
        !
        if (opts%verbosity.ge.1) call am_print_title('Determining action table')
        !
        PM = permutation_map(sg=sg,uc=uc,opts=opts)
        !
        if (present(iopt_fname)) then 
            fid = 1
            open(unit=fid,file=trim(iopt_fname),status="replace",action='write')
        else
            fid = 6
        endif
            !
            write(fid,'(a)') 'Symmerties and coordinate reorientations are in fractional units.'
            ! SYMMETRY NUMBER
            write(fid,'(5x)',advance='no')
            write(fid,'(a5)',advance='no') '#'
            do i = 1, sg%nsyms
                write(fid,'(i10)',advance='no') i
            enddo
            write(fid,*)
            ! HEADER / SEPERATOR
            write(fid,'(5x)',advance='no')
            write(fid,'(5x)',advance='no')
            do i = 1, sg%nsyms
                write(fid,'(a10)',advance='no') ' '//repeat('-',9)
            enddo
            write(fid,*)
            ! CLASS IDENTIFIED
            write(fid,'(5x)',advance='no')
            write(fid,'(a5)',advance='no') 'class'
            do i = 1, sg%nsyms
                write(fid,'(i10)',advance='no') sg%conjugacy_class(i)
            enddo
            write(fid,*)
            ! POINT-SYMMMETRY IDENTIFIED
            write(fid,'(5x)',advance='no')
            write(fid,'(a5)',advance='no') 'SF'
            do i = 1, sg%nsyms
                write(fid,'(a10)',advance='no') trim(decode_pointsymmetry(sg%ps_identifier(i)))
            enddo
            write(fid,*)
            ! POINT SYMMETRY COMPONENTS (FRAC)
            do i = 1,3
            do j = 1,3
                write(fid,'(5x)',advance='no')
                write(fid,'(a5)',advance='no') adjustr('R'//trim(int2char(i))//trim(int2char(j)))
                do k = 1, sg%nsyms
                    write(fid,'(f10.2)',advance='no') sg%R(i,j,k)
                enddo
                write(fid,*)
            enddo
            enddo
            ! TRANSLATIONAL COMPONENTS (FRAC)
            if (allocated(sg%T)) then
            do j = 1,3
                write(fid,'(5x)',advance='no')
                write(fid,'(a5)',advance='no') adjustr('T'//trim(int2char(j)))
                do i = 1, sg%nsyms
                    write(fid,'(f10.2)',advance='no') sg%T(j,i)
                enddo
                write(fid,*)
            enddo
            endif
            ! ACTION TABLE (AXIS PERMUTATION)
            write(fid,'(5x)',advance='no')
            write(fid,'(a5)',advance='no') 'axes'
            do i = 1, sg%nsyms
                write(fid,'(a10)',advance='no') ' '//repeat('-',9)
            enddo
            write(fid,*)
            xyz = 'xyz'
            do j = 1,3
                write(fid,'(5x)',advance='no')
                write(fid,'(a5)',advance='no') trim(xyz(j:j))
                do i = 1, sg%nsyms
                    buffer=''
                    do k = 1,3
                        ! for fractional R; its elements will always be an integer...
                        ! inverse here because f(Rr) = R^-1 * f(r)
                        R=inv(sg%R(:,:,i))
                        if (R(j,k).ne.0) then
                            buffer = trim(buffer)//trim(int2char(nint(R(j,k)),'SP'))//xyz(j:j)
                        endif
                    enddo
                    write(fid,'(a10)',advance='no') trim(buffer)
                enddo
                write(fid,*)
            enddo
            ! ACTION TABLE (ATOM PERMUTATION)
            write(fid,'(5x)',advance='no')
            write(fid,'(a5)',advance='no') 'atoms'
            do i = 1, sg%nsyms
                write(fid,'(a10)',advance='no') ' '//repeat('-',9)
            enddo
            write(fid,*)
            !
            do j = 1, uc%natoms
                write(fid,'(5x)',advance='no')
                write(fid,'(i5)',advance='no') j
                do i = 1, sg%nsyms
                    write(fid,'(i10)',advance='no') PM(j,i)
                enddo
                write(fid,*)
            enddo
            !
        close(fid)
    end subroutine get_action_table

    subroutine     sort(sg,sort_parameter,iopt_direction)
        !
        ! iopt_direction = 'ascend'/'descend', only first character is important
        !
        use am_rank_and_sort
        !
        implicit none
        !
        class(am_class_symmetry) , intent(inout) :: sg
        character(*), intent(in), optional :: iopt_direction
        integer , allocatable :: sorted_indices(:)
        real(dp) :: sort_parameter(:)
        character(1) :: direction
        !
        if (present(iopt_direction)) then
            direction = iopt_direction(1:1)
        else 
            direction = 'a'
        endif
        !
        allocate(sorted_indices(sg%nsyms))
        !
        call rank(sort_parameter,sorted_indices)
        !
        if (direction.eq.'d') sorted_indices = sorted_indices(sg%nsyms:1:-1)
        !
        if (allocated(sg%R))                sg%R                = sg%R(:,:,sorted_indices)
        if (allocated(sg%T))                sg%T                = sg%T(:,sorted_indices)
        if (allocated(sg%ps_identifier)) sg%ps_identifier = sg%ps_identifier(sorted_indices)
        if (allocated(sg%conjugacy_class))  sg%conjugacy_class  = sg%conjugacy_class(sorted_indices)
        !
    end subroutine sort

    !
    ! procedures which operate on representations
    !

    function       get_conjugacy_classes(rep,opts,verbosity) result(conjugacy_class)
        !
        ! AX = XB conjugacy
        ! if rep and rep are conjugate pairs for some other element X in the group, then they are in the same class
        !
        implicit none
        !
        real(dp), intent(in) :: rep(:,:,:) ! list of 2D reps...
        character(*), intent(in), optional :: opts ! can be "seitz", in which case the translational part is reduced to zero
        integer , intent(in) :: verbosity
        real(dp), allocatable :: conjugate_elements(:,:,:)
        integer , allocatable :: conjugacy_class(:)
        integer , allocatable :: conjugacy_classes_relabled(:)
        integer  :: i, j, k, n
        !
        !
        if (verbosity.ge.1) call am_print_title('Determining conjugacy classes')
        !
        n = size(rep,3)
        !
        ! allocate space for conjugacy class
        !
        allocate(conjugacy_class(n))
        conjugacy_class = 0
        !
        ! allocate space for conjugate elements
        !
        allocate(conjugate_elements(size(rep,1),size(rep,2),n))
        !
        ! begin determination of conjugacy classes
        !
        k = 0
        do i = 1, n
        if (conjugacy_class(i).eq.0) then
            !
            k=k+1
            do j = 1, n
                conjugate_elements(:,:,j) = apply_conjugation( center_matrix=rep(:,:,i), side_matrix=rep(:,:,j) )
                if (index(opts,'seitz').ne.0) conjugate_elements(1:3,4,j) = modulo(conjugate_elements(1:3,4,j)+tiny,1.0_dp)-tiny
            enddo
            !
            do j = 1, n
            if (conjugacy_class(j).eq.0) then
            if (issubset( group=conjugate_elements, element=rep(:,:,j) )) then
                conjugacy_class(j) = k
            endif
            endif
            enddo
            !
        endif
        enddo
        !
        conjugacy_classes_relabled = relabel_based_on_occurances(conjugacy_class)
        conjugacy_class = conjugacy_classes_relabled
        !
        ! make sure identity is in the first class
        ! THIS PROCEDURE ONLY WORKS IF IDENTITY IS IN A CLASS BY ITSELF
        ! AND IF THE IDENTITY IS A DIAGONAL MATRIX
        do i = 1, n
            if (all(abs(rep(:,:,i)-eye(n)).lt.tiny)) then
                where (conjugacy_class.eq.1) conjugacy_class = conjugacy_class(i)
                conjugacy_class(i)=1
            endif
        enddo
        !
        if (verbosity.ge.1) then
            call am_print('conjugacy classes',k,' ... ')
            write(*,'(" ... ",a)') 'elements in each class'
            do i = 1, k
                write(*,'(5x,"class",i5," elements",i5)') i, count(i.eq.conjugacy_class)
            enddo
        endif
        !
        if (any(conjugacy_class.eq.0)) then
            call am_print('ERROR','Not every element in the group has been asigned a conjugacy class.',' >>> ')
            stop
        endif
        !
        contains
        function      relabel_based_on_occurances(list) result(A_relabeled)
            !
            use am_rank_and_sort
            !
            implicit none
            !
            integer , intent(in)  :: list(:)
            integer , allocatable :: A_sorted(:)
            integer , allocatable :: occurances(:) 
            integer , allocatable :: reverse_sort(:)
            integer , allocatable :: sorted_indices(:)
            integer , allocatable :: A_relabeled(:)
            integer :: n, i, j
            !
            n = size(list,1)
            !
            allocate(occurances(n))
            do i = 1, n
                occurances(i) = count(list(i).eq.list)
            enddo
            !
            ! this quick and dirty procedure lifts degeneracies
            !
            allocate(sorted_indices(n))
            allocate(A_sorted(n))
            call rank((1+maxval(list))*occurances+list,sorted_indices)
            A_sorted = list(sorted_indices)
            !
            allocate(A_relabeled(n))
            A_relabeled = 0
            !
            j=1
            A_relabeled(1) = 1
            do i = 2, n
                if (A_sorted(i).ne.A_sorted(i-1)) then
                    j=j+1
                endif
                A_relabeled(i) = j
            enddo
            !
            ! return everything to the original order at call
            !
            allocate(reverse_sort(n))
            reverse_sort(sorted_indices)=[1:n]
            A_relabeled=A_relabeled(reverse_sort)
            !
        end function  relabel_based_on_occurances
    end function   get_conjugacy_classes

    function       apply_conjugation(center_matrix,side_matrix) result(C)
            !
            implicit none
            !
            real(dp), intent(in) :: side_matrix(:,:)
            real(dp), intent(in) :: center_matrix(:,:)
            real(dp) :: C(size(side_matrix,1),size(side_matrix,2))
            !
            C = matmul(side_matrix,matmul(center_matrix,inv(side_matrix)))
            !
    end function   apply_conjugation
    !
    ! functions which operate on kpoints
    !

    function       get_kpoint_compatible_symmetries(kpt,R,sym_prec) result(R_syms)
        !
        implicit none
        !
        real(dp), intent(in) :: R(:,:,:)
        real(dp), intent(in) :: kpt(:,:)
        real(dp), intent(in) :: sym_prec
        real(dp), allocatable :: R_syms(:,:,:)
        real(dp) :: wrkspace(3,3,size(R,3))
        integer  :: i, j, nsyms
        !
        nsyms = size(R,3)
        !
        j=0
        do i = 1, nsyms
            if (is_symmetry_valid(tau=kpt,iopt_R=R,iopt_sym_prec=sym_prec)) then
                j=j+1
                wrkspace(:,:,j) = R(:,:,i)
            endif
        enddo
        allocate(R_syms(3,3,j))
        R_syms = wrkspace(:,:,1:j)
        !
    end function   get_kpoint_compatible_symmetries



!     subroutine     coset(C,H,G,index,coset_type)
!         !
!         ! In mathematics, if G is a group, and H is a subgroup of G, and g is an element of G, then:
!         ! gH = { gh : h an element of H } is the left coset of H in G with respect to g, and
!         ! Hg = { hg : h an element of H } is the right coset of H in G with respect to g.
!         ! Wikipedia.
!         !
!         implicit none
!         !
!         class(am_class_symmetry), intent(inout) :: C ! coset of H in G with respect to g
!         type(am_class_symmetry) , intent(in) :: H ! subgroup
!         type(am_class_symmetry) , intent(in) :: G
!         integer                 , intent(in) :: indx
!         real(dp)                , intent(in) :: gg_R(3,3)
!         real(dp)                , intent(in) :: gg_T(3)
!         character(len=1)        , intent(in) :: coset_type ! l (left) / r (right) / c (conjugate)
!         type(am_class_options)  , intent(in) :: opts
!         real(dp), allocatable :: seitz(:,:,:)
!         real(dp) :: hh(4,4)
!         real(dp) :: gg(4,4)
!         integer  :: i
!         !
!         allocate(seitz(4,4,G%nsyms*H%nsyms))
!         !
!         if (opts%verbosity.ge.1) call am_print_title('Determining cosets')
!         !
!         if (opts%verbosity.ge.1) call am_print('symmetries in group',G%nsyms)
!         if (opts%verbosity.ge.1) call am_print('symmetries in subgroup',H%nsyms)
!         !
!         if (modulo(G%nsyms,H%nsyms).ne.0) then
!             call am_print('ERROR','Symmetry subgroup order is not a factor of the group order.',' >>> ')
!             stop
!         endif
!         !
!         ! seitz symmetry operator for G subgroup element gg
!         !
!         gg = 0.0_dp
!         gg(1:3,1:3) = G%R(1:3,1:3,indx)
!         gg(1:3,4) = G%T(1:3,indx)
!         gg(4,4)   = 1.0_dp
!         !
!         do i = 1, H%nsyms
!             ! seitz symmetry operator for H subgroup element hh
!             hh = 0
!             hh(1:3,1:3) = H%R(:,:,i)
!             hh(1:3,4) = H%T(:,i)
!             hh(4,4) = 1.0_dp
!             ! left/right coset, coonjucate 
!             select case (coset_type)
!                 case ('r')
!                     seitz(1:4,1:4,i) = seitz_multiply(A=hh,B=gg)
!                 case ('l')
!                     seitz(1:4,1:4,i) = seitz_multiply(A=gg,B=hh)
!                 case ('c')
!                     seitz(1:4,1:4,i) = seitz_multiply(A=gg,B=hh)
!                 case default
!                     call am_print('ERROR','Coset type not valid')
!                     stop
!             end select
!         enddo
!         !
!         seitz=unique(seitz,iopt_tiny=opts%sym_prec)
!         !
!         C%nsyms = size(seitz,3)
!         if (opts%verbosity.ge.1) call am_print('symmetries in coset subgroup',C%nsyms)
!         !
!         if (modulo(G%nsyms,C%nsyms).ne.0) then
!             call am_print('ERROR','Coset symmetry subgroup order is not a factor of the group order.',' >>> ')
!             stop
!         endif
!         !
!         if (allocated(C%R)) deallocate(C%R)
!         if (allocated(C%T)) deallocate(C%T)
!         allocate(C%R(3,3,C%nsyms))
!         allocate(C%T(3,C%nsyms))
!         !
!         do i = 1, C%nsyms
!             C%R(:,:,i)=seitz(1:3,1:3,i)
!             C%T(:,i)=seitz(1:3,4,i)
!         enddo
!         !
!     end subroutine coset


end module am_symmetry


