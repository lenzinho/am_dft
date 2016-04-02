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
    public :: ps_frac2cart
    !

    type am_class_symmetry
        integer :: nsyms                  !> number of point symmetries
        real(dp), allocatable :: R(:,:,:) !> symmetry elements (operate on fractional atomic basis)
        real(dp), allocatable :: T(:,:)   !> symmetry elements (operate on fractional atomic basis)
        !
        integer , allocatable :: ps_identifier(:) !> indices labeling point symmetries according to their actions (see decode_pointsymmetry)
        integer               :: pg_identifier    !>
        integer , allocatable :: cc_identifier(:) !> indices assiging each element to a conjugacy class
        !
    contains
        procedure :: point_group
        procedure :: space_group
        procedure :: stdout
        procedure :: name_symmetries
        procedure :: sort
        procedure :: create
        !
        procedure :: symmetry_adapted_2nd_order_force_constants
        procedure :: symmetry_adapted_conducitvity
        procedure :: symmetry_adapted_thermoelectricity
        procedure :: symmetry_adapted_elasticity
        procedure :: symmetry_adapted_pizeoelectricity
        !
        procedure :: determine_character_table
        procedure :: action_table_get
        procedure :: stabilizers
    end type

contains

    !
    ! functions which operate on R(3,3) in fractional coordinates
    !

    pure function  ps_frac2cart(R_frac,bas) result(R)
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
    end function   ps_frac2cart

    pure function  ps_cart2frac(R,bas) result(R_frac)
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
    end function   ps_cart2frac

    pure function  ps_determinant(R) result(determinant)
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
    end function   ps_determinant

    function       ps_schoenflies(R) result(ps_identifier)
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
        tr  = trace(R)
        det = ps_determinant(R)
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
            call am_print('ERROR','Unable to identify point symmetry.',flags='E')
            call am_print('R (frac)',R)
            call am_print('tr (frac)',tr)
            call am_print('det (frac)',det)
            stop
        endif
        !
    end function   ps_schoenflies

    !
    ! functions for determining symmetry-adapted second-order force consants
    !

    subroutine     symmetry_adapted_2nd_order_force_constants(sg,uc,opts)
        !
        use am_rank_and_sort
        !
        class(am_class_symmetry), intent(in) :: sg
        type(am_class_unit_cell), intent(in) :: uc
        type(am_class_options), intent(in) :: opts
        !
        integer  :: npairs
        real(dp) :: metric(3,3)
        type(am_class_unit_cell) :: prim
        type(am_class_options) :: notalk
        real(dp), allocatable  :: d(:)
        integer , allocatable  :: p(:,:)
        integer , allocatable  :: uc2prim(:)
        integer , allocatable  :: indices(:)
        real(dp) :: v(3)
        !
        real(dp), allocatable :: M(:,:,:,:)         !> IMPORTANT: M is the original tensor before the permutation (shape of tensor is very important!)
        integer , allocatable :: PM(:,:)            !> PM(uc%natoms,sg%nsyms) permutation map; shows how atoms are permuted by each space symmetry operation
        integer  :: nterms                          !> nterms the number of terms
        integer  :: nsyms                           !> nsyms number of symmetry elements
        real(dp), allocatable :: M1d(:)             !> M1d(nterms) an array which is used to build the matrix M
        real(dp), allocatable :: T(:,:)             !> T(nsyms,nterms) permutation matrx containing rows vectors describing how M is rotated by a point symmetry R
        real(dp), allocatable :: A(:,:)             !> A(nsyms*nterms,2*nterms) augmented matrix equation
        real(dp), allocatable :: LHS(:,:)           !> LHS(nterms,nterms) left hand side of augmented matrix equation (should be identity after reducing to row echlon form)
        real(dp), allocatable :: RHS(:,:)           !> RHS(nterms,nterms) right hand side of augmented matrix equation
        logical , allocatable :: is_dependent(:)    !> is_dependent(nterms) logical array which describes which terms depend on others
        logical , allocatable :: is_zero(:)         !> is_zero(nterms) logical array which shows terms equal to zero
        logical , allocatable :: is_independent(:)  !> is_independent(nterms) logical array which shows which terms are independent
        logical , allocatable :: mask(:)            !> to speed up the symmetry determination process
        integer :: i, j, k, l, n
        integer :: nequations
        integer :: rank_of_T
        !
        if (opts%verbosity.ge.1) call am_print_title('Determining symmetry-adapted 2nd order force constants')
        !
        call am_print('total number of atoms',uc%natoms,' ... ')
        !
        ! get atoms in primitive cell
        notalk = opts 
        notalk%verbosity = 0
        call prim%reduce_to_primitive(uc=uc,oopts_uc2prim=uc2prim,opts=notalk)
        call am_print('number of primitive cell',prim%natoms,' ... ')
        !
        npairs = uc%natoms*prim%natoms
        call am_print('number of pairs formed with primitive cell atoms',npairs,' ... ')
        !
        ! get distances between atoms in primitive cell and in unit cell
        allocate(d(npairs)) ! distances
        allocate(p(2,npairs)) ! pairs
        k=0
        do n = 1, prim%natoms
            i = uc2prim(n)
            do j = 1, uc%natoms
                k=k+1
                ! get distance between atoms (length of bond)
                v = matmul(uc%bas,uc%tau(:,i)-uc%tau(:,j))
                d(k) = norm2(v)
                ! pair indices
                p(1,k) = i
                p(2,k) = j
            enddo
        enddo
        !
        ! sort based on distance
        allocate(indices(npairs))
        call rank(d,indices)
        d=d(indices)
        p(1,:)=p(1,indices)
        p(2,:)=p(2,indices)
        !
        ! write to stdout
        if (opts%verbosity.ge.1) then
            write(*,'(a5,a10,a10,a10)') ' ... ', 'pairs', 'distance'
            write(*,'(5x,a10,a10,a10)') repeat('-',9), repeat('-',9)
            do k = 1, npairs
                if (d(k).gt.tiny) then
                if (d(k).gt.d(k-1)) then
                    write(*,'(5x,a5,a5,f10.5,i10)') trim(uc%symb(uc%atype(p(1,k)))), trim(uc%symb(uc%atype(p(2,k)))), d(k) 
                endif
                endif
            enddo
        endif
        !
        ! now begin the symmetry determination 
        ! determine force constants for all bonds connecting primitive atoms to atoms within a user-specified real-space cutoff distance
        !
        ! permutation map shows how symmetry operations map what atoms to what... 
        PM = action_atom_permutation(sg=sg,uc=uc,opts=opts)
        call am_print('PM',PM)
        !
        stop
        !
        nsyms = size(sg%R,3)
        !
        allocate(M(3,3,prim%natoms,uc%natoms))
        nterms = product(shape(M))
        !
        allocate(M1d(nterms))
        allocate(T(nsyms,nterms))
        allocate(A(nsyms*nterms,2*nterms))
        !
        A = 0
        nequations = 0
        !
!         !$OMP PARALLEL PRIVATE(i,j,M,M1d,T) SHARED(A,uc,sg,PM,nequations)
!         !$OMP DO
!         do i = 1,nterms
!             write(*,*) i/real(nterms)
!             !
!             M1d = 0
!             M1d(i) = 1
!             M = reshape(M1d,shape(M))
!             do j = 1, nsyms
!                 !
!                 ! each point symmetry corresponds to nterms possible equations among terms
!                 !
!                 nequations = nequations + nterms
!                 !
!                 ! apply transformation
!                 !
!                 T(j,1:nterms) = pack( transform_2nd_order_force_constants(&
!                     M=M,R=ps_frac2cart(R_frac=sg%R(:,:,j),bas=uc%bas),PM=PM(:,j)) ,.true.)
!                 !
!             enddo
!             ! save an augmented matrix A
!             A([1:nsyms]+nsyms*(i-1),[1:nterms]+nterms*(1-1)) = T(1:nsyms,1:nterms)
!             A([1:nsyms]+nsyms*(i-1),[1:nterms]+nterms*(2-1)) = 0.0_dp
!             A([1:nsyms]+nsyms*(i-1),         i+nterms*(2-1)) = 1.0_dp
!         enddo
!         !$OMP END DO
!         !$OMP END PARALLEL
!         !
!         ! reduce to row echlon form
!         call rref(A)
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
        do i = 1, natoms
        do j = 1, natoms
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

    subroutine     symmetry_adapted_conducitvity(pg,uc,opts)
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
        logical , allocatable :: mask(:)            !> mask to speed up the symmetry determination process
        integer :: i, j, k, n
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
        !
        ! apply point symmetries (each point symmetry corresponds to nterms possible equations
        ! relating the terms to each other)
        !
        allocate(mask(nterms))
        mask = .true.
        ! <CRYSTAL_SYMMETRIES>
        do i = 1,nterms
            M1d = 0
            M1d(i) = 1
            M = reshape(M1d,shape(M))
            do j = 1, nsyms
                ! Track number of symmetry equations for fun.
                nequations = nequations + nterms
                ! Save the action of the symmetry operations as row-vectors in T(nsyms,nterms) 
                ! Nye, J.F. "Physical properties of crystals: their representation by tensors and matrices". p 133 Eq 7
                T(j,1:nterms) = matmul( kron_pow(ps_frac2cart(R_frac=pg%R(:,:,j),bas=uc%bas),tensor_rank) , M1d)
                ! Apply a mask to make this faster. Once the term is connected to an orbit, all other orbits follow
                ! immediately. There is no need to reapply the point symmetries to the term on which the mapping occurs
                do k = 1, nterms
                    if (abs(T(j,k)).gt.tiny) then
                        mask(k) = .false.
                    endif
                enddo
                !
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
        k=0
        do j = 1,3
        do i = 1,3
            ! sigma_ij = sigma_ji ! permute indices
            nequations = nequations+1
            n  = n+1
            M  = 0.0_dp
            RM = 0.0_dp
            M(i,j) = 1.0_dp
            RM(j,i) = 1.0_dp
            A(n,[1:nterms]+nterms*(1-1))=pack( M , .true. ) ! RHS
            A(n,[1:nterms]+nterms*(2-1))=pack( RM, .true. ) ! LHS
            k=k+1
            T(k,1:nterms) = pack(M,.true.)
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
    end subroutine symmetry_adapted_conducitvity

    pure function  T_op(n) result(M)
        ! c_ij -> c_ji in the flattened basis
        implicit none
        !
        integer, intent(in) :: n
        real(dp), allocatable :: M(:,:)
        integer  :: i, j
        !
        allocate(M(n**2,n**2))
        M=0
        !
        do i = 1, n
        do j = 1, n
           M(i+n*(j-1),j+n*(i-1)) = 1.0_dp
        enddo
        enddo
    end function   T_op

    subroutine     symmetry_adapted_thermoelectricity(pg,uc,opts)
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
        logical , allocatable :: mask(:)            !> mask to speed up the symmetry determination process
        integer :: i, j, k
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
        !
        ! apply point symmetries (each point symmetry corresponds to nterms possible equations
        ! relating the terms to each other)
        !
        allocate(mask(nterms))
        mask = .true.
        ! <CRYSTAL_SYMMETRIES>
        do i = 1,nterms
            M1d = 0
            M1d(i) = 1
            M = reshape(M1d,shape(M))
            do j = 1, nsyms
                !
                ! Track number of symmetry equations for fun.
                nequations = nequations + nterms
                !
                ! Save the action of the symmetry operations as row-vectors in T(nsyms,nterms) 
                ! Nye, J.F. "Physical properties of crystals: their representation by tensors and matrices". p 133 Eq 7
                T(j,1:nterms) = matmul( kron_pow(ps_frac2cart(R_frac=pg%R(:,:,j),bas=uc%bas),tensor_rank) , M1d)
                !
                ! Apply a mask to make this faster. Once the term is connected to an orbit, all other orbits follow
                ! immediately. There is no need to reapply the point symmetries to the term on which the mapping occurs
                do k = 1, nterms
                    if (abs(T(j,k)).gt.tiny) then
                        mask(k) = .false.
                    endif
                enddo
                !
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
    end subroutine symmetry_adapted_thermoelectricity

    subroutine     symmetry_adapted_pizeoelectricity(pg,uc,opts)
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
        logical , allocatable :: mask(:)            !> mask to speed up the symmetry determination process
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
        !
        ! apply point symmetries (each point symmetry corresponds to nterms possible equations
        ! relating the terms to each other)
        !
        allocate(mask(nterms))
        mask = .true.
        ! <CRYSTAL_SYMMETRIES>
        do i = 1,nterms
            M1d = 0
            M1d(i) = 1
            M = reshape(M1d,shape(M))
            do j = 1, nsyms
                !
                ! Track number of symmetry equations for fun.
                nequations = nequations + nterms
                !
                ! Save the action of the symmetry operations as row-vectors in T(nsyms,nterms) 
                ! Nye, J.F. "Physical properties of crystals: their representation by tensors and matrices". p 133 Eq 7
                T(j,1:nterms) = matmul( kron_pow(ps_frac2cart(R_frac=pg%R(:,:,j),bas=uc%bas),tensor_rank) , M1d)
                !
                ! Apply a mask to make this faster. Once the term is connected to an orbit, all other orbits follow
                ! immediately. There is no need to reapply the point symmetries to the term on which the mapping occurs
                do k = 1, nterms
                    if (abs(T(j,k)).gt.tiny) then
                        mask(k) = .false.
                    endif
                enddo
                !
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
    end subroutine symmetry_adapted_pizeoelectricity

    subroutine     symmetry_adapted_elasticity(pg,uc,opts)
        !
        class(am_class_symmetry), intent(in) :: pg
        type(am_class_unit_cell), intent(in) :: uc
        type(am_class_options), intent(in) :: opts
        !
        integer  :: nterms                          !> nterms the number of terms
        integer  :: nsyms                           !> nsyms number of symmetry elements
        integer  :: tensor_rank                     !> tensor_rank rank of tensor matrix M, determined automatically by code
        real(dp), allocatable :: M(:)
        real(dp), allocatable :: R(:,:)
        real(dp), allocatable :: T(:,:)             !> T(nsyms,nterms) permutation matrx containing rows vectors describing how M is rotated by a point symmetry R
        real(dp), allocatable :: A(:,:)             !> A(nsyms*nterms,2*nterms) augmented matrix equation
        real(dp), allocatable :: LHS(:,:)           !> LHS(nterms,nterms) left hand side of augmented matrix equation (should be identity after reducing to row echlon form)
        real(dp), allocatable :: RHS(:,:)           !> RHS(nterms,nterms) right hand side of augmented matrix equation
        logical , allocatable :: is_dependent(:)    !> is_dependent(nterms) logical array which describes which terms depend on others
        logical , allocatable :: is_zero(:)         !> is_zero(nterms) logical array which shows terms equal to zero
        logical , allocatable :: is_independent(:)  !> is_independent(nterms) logical array which shows which terms are independent
        logical , allocatable :: mask(:)            !> mask to speed up the symmetry determination process
        integer :: i, j, k, l, n
        integer :: nequations
        !
        if (opts%verbosity.ge.1) call am_print_title('Determining elasticity tensor')
        !
        tensor_rank = 4
        nsyms = size(pg%R,3)
        nterms = 3**tensor_rank
        !
        allocate(A(2*nterms,2*nterms)) ! augmented matrix
        A = 0
        !
        allocate(R(nterms,nterms))
        R = 0
        !
        nequations = 0
        !
        ! crystal symmetries
        A([1:nterms],[1:nterms]+nterms*0) = eye(nterms)
        A([1:nterms],[1:nterms]+nterms*1) = eye(nterms)
        call rref(A)
        do i = 1, nsyms
            ! Track number of symmetry equations for fun
            nequations = nequations + nterms
            !
            R = kron_pow(ps_frac2cart(R_frac=pg%R(:,:,i),bas=uc%bas),tensor_rank)
            ! Save the action of the symmetry operations as row-vectors in T(nsyms,nterms) 
            ! Nye, J.F. "Physical properties of crystals: their representation by tensors and matrices". p 133 Eq 7
            A([1:nterms]+nterms,[1:nterms]+nterms*(1-1)) = R
            A([1:nterms]+nterms,[1:nterms]+nterms*(2-1)) = eye(nterms)
            ! reduced to row echlon form
            call am_print('R',R,filename='debug_R'//trim(int2char(i))//'.txt',permission='w')
            call am_print('B',A,filename='debug_B'//trim(int2char(i))//'.txt',permission='w')
            call rref(A)
            call am_print('A',A,filename='debug_A'//trim(int2char(i))//'.txt',permission='w')
            ! call am_print('R',R,filename='debug_R'//trim(int2char(i))//'.txt',permission='w')
            ! call am_print('A',A,filename='debug_A'//trim(int2char(i))//'.txt',permission='w')
        enddo
        !
        ! now apply intrinsic symmetries, which are fundamental properties of the
        ! elastic/stress/strain tensors and, therefore, independent of the crystal
        ! Nye, J.F. "Physical properties of crystals: their representation by tensors and matrices". p 133 Eq 5 and 6
        ! cijkl = cjikl
        nequations = nequations+nterms
        call am_print('nterms+[1:nterms]',nterms+[1:nterms])
        call am_print('[1:nterms]+nterms*0)',[1:nterms]+nterms*0)
        A(nterms+[1:nterms],[1:nterms]+nterms*0) = kron(eye(9),T_op(3))
        A(nterms+[1:nterms],[1:nterms]+nterms*1) = eye(nterms)
        call rref(A)
        ! cijkl = cijlk
        nequations = nequations+nterms
        A(nterms+[1:nterms],[1:nterms]+nterms*0) = kron(eye(9),T_op(3))
        A(nterms+[1:nterms],[1:nterms]+nterms*1) = eye(nterms)
        call rref(A)
        ! cijkl = cjilk
        nequations = nequations+nterms
        A(nterms+[1:nterms],[1:nterms]+nterms*0) = kron(T_op(3),T_op(3))
        A(nterms+[1:nterms],[1:nterms]+nterms*1) = eye(nterms)
        call rref(A)
        ! cijkl = cijkl
        nequations = nequations+nterms
        A(nterms+[1:nterms],[1:nterms]+nterms*0) = eye(nterms)
        A(nterms+[1:nterms],[1:nterms]+nterms*1) = eye(nterms)
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
    end subroutine symmetry_adapted_elasticity

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
            call am_print('ERROR','The number of terms which are independent, null, and dependent do not add to the total number of terms.',flags='E')
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
            call am_print('ERROR','Schoenflies code unknown.',flags='E')
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
                call am_print('ERROR','Schoenflies code for point-group unknown.',flags='E')
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
            call am_print('ERROR','More than one inversion operator found.',flags='E')
            stop
        endif
        !
        if (ne.gt.1) then
            call am_print('ERROR','More than one identity operator found.',flags='E')
        elseif (ne.eq.0) then
            call am_print('ERROR','No identity operator found.',flags='E')
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
        call am_print('ERROR','Unable to identify point group',flags='E')
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
        call pg%point_group(sg=sg,uc=uc,opts=opts)
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
        real(dp), allocatable  :: seitz(:,:,:)
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
        sg%cc_identifier = cc_get(rep=seitz,flags='seitz')
        !
        ! sort space group symmetries based on parameters
        call sg%sort(sort_parameter=sg%T(1,:),iopt_direction='ascend')
        call sg%sort(sort_parameter=sg%T(2,:),iopt_direction='ascend')
        call sg%sort(sort_parameter=sg%T(3,:),iopt_direction='ascend')
        call sg%sort(sort_parameter=real(sg%cc_identifier,dp),iopt_direction='ascend')
        !
        ! name the (im-)proper part of space symmetry
        call sg%name_symmetries(opts=notalk)
        !
        ! write action table
        call sg%action_table_get(uc=uc,iopt_fname='outfile.action_space_group',opts=opts)
        !
        ! get character table
        call sg%determine_character_table(opts=opts)
        !
        ! write to stdout and to file
        ! if (opts%verbosity.ge.1) call sg%stdout(iopt_uc=uc)
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
        pg%cc_identifier = cc_get(rep=pg%R,flags='')
        !
        ! sort point group symmetries based on parameters
        allocate(sort_parameter(pg%nsyms))
        do i = 1, pg%nsyms; sort_parameter(i) = trace(pg%R(:,:,i)); enddo
        call pg%sort(sort_parameter=sort_parameter,iopt_direction='ascend')
        do i = 1, pg%nsyms; sort_parameter(i) = ps_determinant(pg%R(:,:,i)); enddo
        call pg%sort(sort_parameter=sort_parameter,iopt_direction='ascend')
        call pg%sort(sort_parameter=real(pg%cc_identifier,dp),iopt_direction='ascend')
        !
        ! name point symmetries
        call pg%name_symmetries(opts=opts)
        !
        ! name point group
        pg%pg_identifier = point_group_schoenflies(pg%ps_identifier)
        if (opts%verbosity.ge.1) call am_print('point group',decode_pointgroup(pg%pg_identifier),' ... ')
        !
        ! determine character table
        call pg%determine_character_table(opts=opts)
        !
        ! write action table 
        ! currently not working for nonsymorphic space groups. (i.e. those which necessairly have a translational component)
        ! actually it may just not be possible to do this for such a group because permutations linking atoms simply do not exist if only the rotational part is considered. the translational part is essential. an example is TiSe2.
        ! call pg%action_table_get(uc=uc,iopt_fname='outfile.action_point_group',opts=opts)
        !
        ! write to stdout and to file
        ! if (opts%verbosity.ge.1) call pg%stdout(iopt_uc=uc)
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
            call am_print('ERROR','Inconsistent number of atom pairs. Something is wrong internally.',flags='E')
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
                write(*,'(2a5)',advance='no') trim(uc%symb(uc%atype(i))), trim(uc%symb(uc%atype(j)))
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

    subroutine     determine_character_table(sg,opts)
        !
        use am_rank_and_sort
        !
        implicit none
        !
        class(am_class_symmetry), intent(inout) :: sg
        type(am_class_options)  , intent(in) :: opts
        real(dp), allocatable :: rep(:,:,:)
        real(dp), allocatable :: CT(:,:)
        character(100) :: flags
        integer :: i, j, k
        !
        integer :: nsyms
        integer :: nclasses
        integer :: nirreps  ! just to make things clearer; always equal to nclasses
        integer, allocatable :: nelements(:) ! number of elements in each class
        integer, allocatable :: member(:,:) ! members(nclass,maxval(nelements))
        integer, allocatable :: H(:,:,:) ! class multiplication coefficients
        integer, allocatable :: CCT(:,:) ! class constant table which becomes character table (can be complex!)
        integer, allocatable :: irrep_dim(:) ! rep dimensions
        integer, allocatable :: indices(:)
        real(dp) :: wrk

        !
        if (opts%verbosity.ge.1) call am_print_title('Determining character table')
        !
        ! seitz rep
        if (allocated(sg%T)) then
            rep = rep_seitz(sg%R,sg%T)
            flags='seitz'
        else
            rep = sg%R
            flags=''
        endif
        !
        nsyms = size(rep,3)
        !
        ! identify which class symmetry belongs to
        ! cc_identifier(nsyms)
        sg%cc_identifier = cc_get(rep=rep,flags=flags)
        !
        ! get number of classes
        nclasses = maxval(sg%cc_identifier)
        ! if (opts%verbosity.ge.1) call am_print('conjugacy classes',nclasses,' ... ')
        !
        ! get number of elements in each class
        ! nelements(nclasses)
        nelements = cc_nelements(sg%cc_identifier)
        !
        ! record indicies of each member element for each class (use the first member of class as class representative)
        ! members(nclass,maxval(nelements)) 
        member = cc_member(sg%cc_identifier)
        ! call am_print('class memebers',member,' ... ')
        !
        ! construct cayley table for quick access to symmetry multiplication
        CT = cayley_table(rep=rep,flags=flags)
        !
        ! get class mutliplication coefficients
        allocate(H(nclasses,nclasses,nclasses))
        do i = 1,nclasses
        do j = 1,nclasses
        do k = 1,nclasses
            H(j,k,i) = count( pack( CT(member(i,1:nelements(i)),member(j,1:nelements(j))) , .true. ) .eq. member(k,1) )
        enddo
        enddo
        enddo
        !
        ! determine eigenvector which simultaneously diagonalize i Hjk(:,:) matrices
        ! CCT(nirreps,nclasses)
        CCT = cc_constant_table(real(H,dp))
        !
        ! get dimensions of representations using Eq. 4.53, p 91 of Wooten
        ! Note: The absolute square is the result of the conjugate multiplication
        nirreps = nclasses
        allocate(irrep_dim(nirreps))
        !
        do i = 1, nirreps
            ! sum on classes
            wrk = 0
            do j = 1, nclasses
                wrk = wrk + abs(CCT(i,j))**2/real(nelements(j),dp)
            enddo
            ! compute irrep dimension
            irrep_dim(i) = nint((nsyms/wrk)**0.5)
        enddo
        !
        ! convert class constant table into character table
        do i = 1, nirreps
        do j = 1, nclasses
            CCT(i,j) = nint(irrep_dim(i)/real(nelements(j),dp)*CCT(i,j))
        enddo
        enddo
        !
        ! sort irreps based on dimension
        allocate(indices(nirreps))
        indices = [1:nirreps]
        call rank(irrep_dim,indices)
        irrep_dim = irrep_dim(indices)
        do j = 1,nclasses
            CCT(:,j) = CCT(indices,j)
        enddo
        !
        ! write class properties to stdout
        !
        if (opts%verbosity.ge.1) then
            !
            write(*,'(5x,a10)',advance='no') 'class'
            do i = 1, nclasses
                write(*,'(i5)',advance='no') i
            enddo
            write(*,*)
            !
            write(*,'(5x,a10)',advance='no') 'elements'
            do i = 1, nclasses
                write(*,'(i5)',advance='no') nelements(i)
            enddo
            write(*,*)
            !
            write(*,'(5x,a10)',advance='no') 'class rep'
            do i = 1, nclasses
                write(*,'(i5)',advance='no') member(i,1)
            enddo
            write(*,*)
            !
            write(*,'(5x,a10)',advance='no') repeat('-',10)
            do i = 1, nclasses
                write(*,'(a5)',advance='no') repeat('-',5)
            enddo
            write(*,*)
            !
            do i = 1, nirreps
                write(*,'(5x,a6,i4)',advance='no') 'irrep', i
                do j = 1, nclasses
                    write(*,'(i5)',advance='no') CCT(i,j)
                enddo
                write(*,*)
            enddo
            !
        endif
        ! check that the number of elements ri in class Ci is a divisor of theorder of the group
        ! Symmetry and Condensed Matter Physics: A Computational Approach. 1 edition. Cambridge, UK?; New York: Cambridge University Press, 2008. p 35.
        do i = 1, nclasses
            if (modulo(nsyms,nelements(i)).ne.0) then
                call am_print('ERROR','Number of elements in class subgroup is not a divisor of the order of the group.')
                call am_print('class',i)
                call am_print('elements in class',nelements(i))
                call am_print('elements in group',nsyms)
                stop
            endif
        enddo
        !
    contains
    function       cc_nelements(cc_identifier) result(nelements)
        !
        implicit none
        !
        integer, intent(in)  :: cc_identifier(:)
        integer, allocatable :: nelements(:)
        integer :: nclasses
        integer :: i
        !
        nclasses = maxval(cc_identifier)
        !
        allocate(nelements(nclasses))
        !
        do i = 1, nclasses
            nelements(i) = count(i.eq.cc_identifier)
        enddo
    end function   cc_nelements
    function       cc_member(cc_identifier) result(member)
        !
        implicit none
        !
        integer, allocatable, intent(in) :: cc_identifier(:)
        integer, allocatable :: nelements(:) ! number of elements in each class
        integer, allocatable :: member(:,:) ! members(nclass,maxval(nelements))
        integer :: i, j, k
        integer :: nsyms
        integer :: nclasses
        !
        nsyms = size(cc_identifier,1)
        !
        nclasses = maxval(cc_identifier)
        !
        nelements = cc_nelements(cc_identifier)
        !
        allocate(member(nclasses,maxval(nelements)))
        member = 0
        do i = 1, nclasses
            k=0
            do j = 1, nsyms
                if (cc_identifier(j).eq.i) then
                    k=k+1
                    member(i,k) = j
                endif
            enddo
        enddo
    end function   cc_member
    function       cc_constant_table(A) result(D)
        !
        ! determines characters which simultaneously diagonalizes i square matrices A(:,:,i)
        ! only possible if A(:,:,i) commute with each other. Assumes all dimensions of A are the same.
        !
        implicit none
        !
        real(dp), intent(in) :: A(:,:,:)
        real(dp), allocatable :: V(:,:) ! these two should match in this particular situation (character table determination from class multiplication coefficients)
        integer , allocatable :: D(:,:) ! these two should match in this particular situation (character table determination from class multiplication coefficients)
        real(dp), allocatable :: M(:,:)
        integer , allocatable :: p(:)
        integer , allocatable :: small(:)
        integer :: n, i
        !
        call am_print('WARNING','The procedure used to determine the character table does not allow for complex character',flags='W')
        !
        n = size(A,3)
        p = primes(n)
        !
        ! matrix pencil based on sqrt(primes)
        ! all this does is lift the degeneracy of eigenvalues if necessary
        allocate(M(n,n))
        M = 0
        do i = 1,n
            M = M + p(i)**0.5*A(:,:,i)
        enddo
        !
        ! get left and right eigenvectors: transpose(VL)*A*VR = D
        call am_dgeev(A=M,VR=V)
        !
        ! check that matrices have been diagonalized properly and save eigenvalues
        allocate(d(n,n))
        do i = 1, n
            !
            M = matmul(inv(V),matmul(A(:,:,i),V))
            ! tiny may be too small a criteria here... should probably use a value which scales with the size of the matrix.
            if (any( abs(diag(diag(M))-M).gt. (tiny*n**2) )) then
                call am_print('ERROR','Unable to perform simultaneous matrix diagonalization. Check that whether they commute.',flags='E')
                call am_print('M',M)
                call am_print('diag(M)',diag(M))
                call am_print('diag(diag(M))',diag(diag(M)))
                call am_print('abs(diag(diag(M))-M)',abs(diag(diag(M))-M))
                call am_print('sum of error',sum(abs(diag(diag(M))-M)))
                stop
            endif
            !
            ! save diagonal elements as row vectors
            D(:,i) = nint(diag( M ))
            !
        enddo
        !
        ! normalize each row of V to the smallest value in that row
        allocate(small(n))
        small = minloc( abs(V) , dim=1, mask = abs(V).gt.tiny )
        do i = 1, n
            V(:,i)=V(:,i)/V(small(i),i)
        enddo
        V = nint(transpose(V))
        !
        ! at this point V and D should be identical. The elements are to within +/- sign.
        ! Return D which is more robust than V. 
        !
    end function   cc_constant_table
    end subroutine determine_character_table

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
            R = ps_frac2cart(R_frac=sg%R(:,:,i),bas=iopt_uc%bas)
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
            sg%ps_identifier(i) = ps_schoenflies(sg%R(:,:,i))
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

    function       action_atom_permutation(sg,uc,opts) result(PM)
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
        P = rep_permutation(sg=sg,uc=uc,opts=opts)
        !
        allocate(PM(uc%natoms,sg%nsyms))
        PM = 0
        do i = 1,sg%nsyms
            PM(:,i) = matmul(P(:,:,i),[1:uc%natoms])
        enddo
        !
    end function   action_atom_permutation

    subroutine     action_table_get(sg,uc,iopt_fname,opts)
        !
        implicit none
        !
        class(am_class_symmetry), intent(inout) :: sg
        type(am_class_unit_cell), intent(in) :: uc
        type(am_class_options)  , intent(in) :: opts
        character(*), optional  , intent(in) :: iopt_fname
        integer, allocatable :: PM(:,:)
        integer :: fid
        character(10) :: buffer
        character(3) :: xyz
        integer :: i,j,k
        real(dp) :: R(3,3)
        !
        if (opts%verbosity.ge.1) call am_print_title('Determining action table')
        !
        PM = action_atom_permutation(sg=sg,uc=uc,opts=opts)
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
            ! CLASS IDENTIFIED
            write(fid,'(5x)',advance='no')
            write(fid,'(a5)',advance='no') 'class'
            do i = 1, sg%nsyms
                write(fid,'(i10)',advance='no') sg%cc_identifier(i)
            enddo
            write(fid,*)
            ! HEADER / SEPERATOR
            write(fid,'(5x)',advance='no')
            write(fid,'(5x)',advance='no')
            do i = 1, sg%nsyms
                write(fid,'(a10)',advance='no') ' '//repeat('-',9)
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
    end subroutine action_table_get

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
        if (allocated(sg%R)) sg%R = sg%R(:,:,sorted_indices)
        if (allocated(sg%T)) sg%T = sg%T(:,sorted_indices)
        if (allocated(sg%ps_identifier)) sg%ps_identifier = sg%ps_identifier(sorted_indices)
        if (allocated(sg%cc_identifier)) sg%cc_identifier = sg%cc_identifier(sorted_indices)
        !
    end subroutine sort

    !
    ! procedures which create representations
    !

    function       rep_seitz(R,T) result(seitz)
        !
        implicit none
        !
        real(dp), intent(in) :: R(:,:,:)
        real(dp), intent(in) :: T(:,:)
        real(dp), allocatable :: seitz(:,:,:)
        integer :: i, n
        !
        n = size(R,3)
        !
        allocate(seitz(4,4,n))
        seitz = 0
        seitz(4,4,:)=1
        do i = 1, n
            seitz(1:3,1:3,i)=R(:,:,i)
            seitz(1:3,4,i)  =T(:,i)
        enddo
        !
    end function   rep_seitz

    function       rep_permutation(sg,uc,opts) result(rep)
       !
       ! find permutation representation; i.e. which atoms are connected by space symmetry oprations R, T.
       !
       implicit none
       !
       type(am_class_symmetry) , intent(in) :: sg
       type(am_class_unit_cell), intent(in) :: uc
       type(am_class_options)  , intent(in) :: opts
       integer, allocatable :: rep(:,:,:)
       real(dp) :: tau_rot(3)
       integer :: i,j,k
       logical :: found
       !
       allocate(rep(uc%natoms,uc%natoms,sg%nsyms))
       rep = 0
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
                    rep(j,k,i) = 1
                    found = .true.
                    exit ! break loop
                    endif
                enddo
                if (found.eq..false.) then
                    call am_print('ERROR','Unable to find matching atom.',flags='E')
                    call am_print('tau (all atoms)',transpose(uc%tau))
                    call am_print('tau',uc%tau(:,j))
                    call am_print('R',sg%R(:,:,i))
                    if (allocated(sg%T)) call am_print('T',sg%T(:,i))
                    call am_print('tau_rot',tau_rot)
                    stop
                endif
           enddo
       enddo
       !
       ! check that each column and row of the rep sums to 1; i.e. rep(:,:,i) is orthonormal for all i.
       !
       do i = 1,sg%nsyms
       do j = 1,uc%natoms
           !
           if (sum(rep(:,j,i)).ne.1) then
               call am_print('ERROR','Permutation matrix has a column which does not sum to 1.')
               call am_print('i',i)
               call am_print_sparse('spy(P_i)',rep(:,:,i))
               call am_print('rep',rep(:,:,i))
               stop
           endif
           !
           if (sum(rep(j,:,i)).ne.1) then
               call am_print('ERROR','Permutation matrix has a row which does not sum to 1.')
               call am_print('i',i)
               call am_print_sparse('spy(P_i)',rep(:,:,i))
               call am_print('rep',rep(:,:,i))
               stop
           endif
           !
       enddo
       enddo
       !
    end function   rep_permutation

    function       rep_regular(rep,flags) result(reg_rep)
        !>
        !> Generates regular representation for each operator
        !>
        !> For details, see page 73-74, especially Eqs. 4.18 Symmetry and Condensed Matter Physics: A Computational
        !> Approach. 1 edition. Cambridge, UK; New York: Cambridge  University Press, 2008.
        !>
        !> Requires identity to be the first element of group. regular representation is built from multiplication table
        !> by making sure that the identity appears along the diagonal. this sometimes causes the row to be permuted with
        !> respect to the columns. so, first find the permutation necessary to bring the identities to the center and
        !> then build regular representations from this new multiplication table.
        !>
        implicit none
        !
        real(dp), intent(in) :: rep(:,:,:)
        character(*), intent(in) :: flags ! can be nothing flags='' at call
        integer, allocatable :: reg_rep(:,:,:)
        integer, allocatable :: CT(:,:)
        integer  :: n
        integer  :: i
        !
        ! obtain cayley table
        CT = cayley_table(rep=rep,flags=flags//'reg')
        !
        ! construct regular representation from cayley table
        n = size(rep,3)
        allocate(reg_rep(n,n,n))
        reg_rep = 0
        do i = 1, n
            where (CT.eq.i) reg_rep(:,:,i) = 1
        enddo
        !
    end function   rep_regular

    !
    ! procedures which operate on representations
    !

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

    function       cayley_table(rep,flags) result(CT)
        !
        implicit none
        !
        real(dp), intent(in) :: rep(:,:,:) ! list of 2D reps...
        character(*), intent(in) :: flags ! 'reg' sorts rows to put identity along diagonals (useful for constructing regular representation), 'seitz' reduces translational part to primitive cell
        integer, allocatable :: CT(:,:)
        real(dp) :: W(size(rep,1),size(rep,2))
        integer, allocatable :: sortmat(:,:)
        integer, allocatable :: indices(:)
        integer  :: n
        integer  :: i, j
        !
        ! before doing anything, confirm that first element is the identity
        if (any(abs(rep(:,:,1)-eye(size(rep,1))).gt.tiny)) then
            call am_print('ERROR','First element of rep is not the identity.')
            call am_print('first element of rep',rep(:,:,1))
            stop
        endif
        !
        n=size(rep,3)
        !
        allocate(CT(n,n))
        !
        do i = 1, n
        do j = 1, n
            !
            W = matmul(rep(:,:,i),rep(:,:,j))
            !
            ! if 'seitz' flagged, reduce translational part to primitive cell 
            if (index(flags,'seitz').ne.0) then
                W(1:3,4) = modulo(W(1:3,4)+tiny,1.0_dp)-tiny
            endif
            !
            CT(i,j) = get_matching_element_index_in_list(list=rep,elem=W)
        enddo
        enddo
        !
        ! if 'reg' flagged, re-order rows so that identities are along diagonal
        if (index(flags,'reg').ne.0) then
            !
            allocate(sortmat(n,n))
            allocate(indices(n))
            indices = [1:n]
            sortmat = 0
            where (CT.eq.1) sortmat(:,:) = 1
            indices = matmul(sortmat(:,:),indices)
            CT = CT(indices,:)
        endif
        !
    contains
    function       get_matching_element_index_in_list(list,elem) result(i)
        !
        implicit none
        real(dp), intent(in) :: list(:,:,:)
        real(dp), intent(in) :: elem(:,:)
        integer :: i
        !
        do i = 1,size(list,3)
            if (all(abs(list(:,:,i)-elem).lt.tiny)) return
        enddo
    end function   get_matching_element_index_in_list
    end function   cayley_table

    !
    ! functions which operate on conjugacy classes
    !

    function       cc_get(rep,flags) result(cc_identifier)
        !
        ! AX = XB conjugacy
        ! if rep and rep are conjugate pairs for some other element X in the group, then they are in the same class
        !
        implicit none
        !
        real(dp), intent(in) :: rep(:,:,:)    ! list of matrices representing group reps
        character(*), intent(in) :: flags     ! can be 'seitz', in which case the translational part is reduced to zero
        real(dp), allocatable :: celem(:,:,:) ! conjugate elements
        integer , allocatable :: cc_identifier(:)
        integer  :: i, j, k, n, nclasses
        !
        n = size(rep,3)
        !
        ! allocate space for conjugacy class
        allocate(cc_identifier(n))
        cc_identifier = 0
        !
        ! allocate space for conjugate elements
        allocate(celem(size(rep,1),size(rep,2),n))
        !
        ! determine conjugacy classes
        k = 0
        do i = 1, n
        if (cc_identifier(i).eq.0) then
            !
            ! conjugate each element with all other group elements
            k=k+1
            do j = 1, n
                celem(:,:,j) = apply_conjugation( center_matrix=rep(:,:,i), side_matrix=rep(:,:,j) )
                if (index(flags,'seitz').ne.0) then
                    celem(1:3,4,j) = modulo(celem(1:3,4,j)+tiny,1.0_dp)-tiny
                endif
            enddo
            !
            ! for each subgroup element created by conjugation find the corresponding index of the element in the group
            ! in order to save the class identifier number
            do j = 1, n
            if (cc_identifier(j).eq.0) then
            if (issubset( group=celem, element=rep(:,:,j) )) then
                cc_identifier(j) = k
            endif
            endif
            enddo
            !
        endif
        enddo
        nclasses = k
        !
        ! relabel classes based on number of elements in each class
        call relabel_based_on_occurances(cc_identifier)
        !
        ! make sure identity is in the first class
        do i = 1, n
            if (all(abs(rep(:,:,i)-eye(n)).lt.tiny)) then
                where (cc_identifier.eq.1) cc_identifier = cc_identifier(i)
                cc_identifier(i)=1
            endif
        enddo
        !
        ! consistency checks
        ! 1)
        if (any(cc_identifier.eq.0)) then
            call am_print('ERROR','Not every element in the group has been asigned a conjugacy class.',flags='E')
            stop
        endif
        !
        contains
        subroutine     relabel_based_on_occurances(list)
            !
            use am_rank_and_sort
            !
            implicit none
            !
            integer , intent(inout) :: list(:)
            integer , allocatable :: A_sorted(:)
            integer , allocatable :: occurances(:) 
            integer , allocatable :: reverse_sort(:)
            integer , allocatable :: sorted_indices(:)
            integer , allocatable :: list_relabled(:)
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
            allocate(list_relabled(n))
            list_relabled = 0
            !
            j=1
            list_relabled(1) = 1
            do i = 2, n
                if (A_sorted(i).ne.A_sorted(i-1)) then
                    j=j+1
                endif
                list_relabled(i) = j
            enddo
            !
            ! return everything to the original order at call
            !
            allocate(reverse_sort(n))
            reverse_sort(sorted_indices)=[1:n]
            list=list_relabled(reverse_sort)
        end subroutine relabel_based_on_occurances
    end function   cc_get

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
!             call am_print('ERROR','Symmetry subgroup order is not a factor of the group order.',flags='E')
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
!             call am_print('ERROR','Coset symmetry subgroup order is not a factor of the group order.',flags='E')
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


