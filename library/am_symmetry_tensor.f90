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
        !
        if (opts%verbosity.ge.1) call print_title('Determining second-rank tensor')
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
        !
        do i = 1,nterms
            M1d = 0
            M1d(i) = 1
            M = reshape(M1d,shape(M))
            do j = 1, nsyms
                ! RM shows how the M has been permuted by the symmetry operation j
                ! save the permutation obtained with RM as a row-vectors in T(nsyms,nterms) 
                ! T(j,1:nterms) = pack( transform_tensor(M=M,R=pg%R(:,:,j)),.true.)
                T(j,1:nterms) = pack( transform_tensor(M=M,R=point_symmetry_convert_to_cartesian(R_frac=pg%R(:,:,j),bas=uc%bas)) ,.true.)
            enddo
            ! save an augmented matrix A
            A([1:nsyms]+nsyms*(i-1),[1:nterms]+nterms*(1-1)) = T(1:nsyms,1:nterms)
            A([1:nsyms]+nsyms*(i-1),[1:nterms]+nterms*(2-1)) = 0.0_dp
            A([1:nsyms]+nsyms*(i-1),         i+nterms*(2-1)) = 1.0_dp
        enddo
        !
        call rref(A)
        ! At this point A is an augmented matrix composed of [ LHS | RHS ]. The LHS
        ! should be the identity matrix, which, together with the RHS, completely 
        ! specifies all relationships between variables. 
        allocate(LHS(nterms,nterms))
        allocate(RHS(nterms,nterms))
        LHS = A(1:nterms,1:nterms)
        RHS = A(1:nterms,[1:nterms]+nterms)
        !
        call print_sparse('A = [LHS|RHS]',A(1:nterms,:),' ... ')
        !
        if ( any(abs(LHS-eye(nterms)).gt.tiny) ) then
            call am_print('ERROR','Unable to reduce matrix to row echlon form.')
            call print_sparse('spy(LHS)',LHS)
            stop
        endif
        !
        call parse_symmetry_equations(LHS=LHS,RHS=RHS,is_zero=is_zero,&
            is_independent=is_independent,is_dependent=is_dependent,verbosity=opts%verbosity)
        !