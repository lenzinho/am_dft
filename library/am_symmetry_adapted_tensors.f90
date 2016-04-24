module am_symmetry_adapted_tensors

    use am_constants
    use am_stdout
    use am_unit_cell
    use am_options
    use am_mkl
    use am_symmetry

	implicit none
	
	private

	public :: get_symmetry_adapted_tensor

contains

    subroutine     get_symmetry_adapted_tensor(pg,uc,opts,relations,property)
        !
        type(am_class_symmetry), intent(in) :: pg
        type(am_class_unit_cell), intent(in) :: uc
        type(am_class_options), intent(in) :: opts
        character(*), intent(in) :: property
        real(dp), allocatable, optional, intent(out) :: relations(:,:)
        !
        integer  :: ndim
        integer  :: nterms                          !> nterms the number of terms
        integer  :: tensor_rank                     !> tensor_rank rank of tensor matrix M, determined automatically by code
        integer , allocatable :: class_member(:,:)
        real(dp), allocatable :: R(:,:)
        real(dp), allocatable :: S(:,:,:)           !> intrinsic symmetries
        real(dp), allocatable :: A(:,:)             !> A(2*nterms,2*nterms) augmented matrix equation
        real(dp), allocatable :: LHS(:,:)           !> LHS(nterms,nterms) left hand side of augmented matrix equation (should be identity after reducing to row echlon form)
        real(dp), allocatable :: RHS(:,:)           !> RHS(nterms,nterms) right hand side of augmented matrix equation
        logical , allocatable :: is_dependent(:)    !> is_dependent(nterms) logical array which describes which terms depend on others
        logical , allocatable :: is_zero(:)         !> is_zero(nterms) logical array which shows terms equal to zero
        logical , allocatable :: is_independent(:)  !> is_independent(nterms) logical array which shows which terms are independent
        integer , allocatable :: indices(:)
        character(10) :: flags
        integer :: i, j, k
        integer :: nequations
        !
        !
        if (opts%verbosity.ge.1) call am_print_title('Determining symmetry-adapted tensor: '//trim(property))
        !------------------------------------------------- FIRST-RANK TENSORS --------------------------------------------------
        if     (index(property,'pyroelectricity')          .ne.0) then; tensor_rank = 1 !    ! P_{i}     = p_{i} \Delta T
        !------------------------------------------------- SECOND-RANK TENSORS -------------------------------------------------
        elseif (index(property,'pair force-constants')     .ne.0) then; tensor_rank = 2; flags = 'polar' ! second-order force-constants
        elseif (index(property,'pair reversal')            .ne.0) then; tensor_rank = 2; flags = 'polar' ! second-order force-constants corresponding to reversal group (will apply flip operator!)
        elseif (index(property,'electrical susceptibility').ne.0) then; tensor_rank = 2; flags = 'polar' ! S  ! P_{i}     = \alpha_{ij}  E_{j}
        elseif (index(property,'magnetic susceptibility')  .ne.0) then; tensor_rank = 2; flags = 'axial' ! S  ! M_{i}     = \mu_{ij}     H_{j}
        elseif (index(property,'magneto-electric')         .ne.0) then; tensor_rank = 2; flags = 'axial' 
        elseif (index(property,'thermal expansion')        .ne.0) then; tensor_rank = 2; flags = 'polar' ! S  ! \eps_{ij} = \alpha_{ij}  \Delta T
        elseif (index(property,'electrical conductivity')  .ne.0) then; tensor_rank = 2; flags = 'polar' ! S  ! J_{i}     = \sigma_{ij}  E_{i}
        elseif (index(property,'electrical resistivity')   .ne.0) then; tensor_rank = 2; flags = 'polar' ! S  ! E_{i}     = \rho_{ij}    J_{j}
        elseif (index(property,'thermal conductivity')     .ne.0) then; tensor_rank = 2; flags = 'polar' ! S  ! q_{i}     = \kappa_{ij}  \frac{\partial T}/{\partial r_{j}}
        elseif (index(property,'thermoelectricity')        .ne.0) then; tensor_rank = 2; flags = 'polar' ! N  ! 
        elseif (index(property,'seebeck')                  .ne.0) then; tensor_rank = 2; flags = 'polar' ! N  ! E_{i}     = \beta_{ij}   \frac{\partial T}/{\partial r_{j}}
        elseif (index(property,'peltier')                  .ne.0) then; tensor_rank = 2; flags = 'polar' ! N  ! q_{i}     = \pi_{ij}     J_{j}
        !------------------------------------------------- THIRD-RANK TENSORS -------------------------------------------------
        elseif (index(property,'hall')                     .ne.0) then; tensor_rank = 3;                 !    ! E_{i}     = h_{ijk}      J_{j} H_{k} 
        elseif (index(property,'piezoelectricity')         .ne.0) then; tensor_rank = 3; flags = 'polar' !    ! P_{i}     = d_{ijk}      \sigma_{jk}
        elseif (index(property,'piezomagnetic')            .ne.0) then; tensor_rank = 3; flags = 'axial' !    ! M_{i}     = Q_{ijk}      \sigma_{jk}
        !------------------------------------------------- FOURTH-RANK TENSORS ------------------------------------------------
        elseif (index(property,'elasticity')               .ne.0) then; tensor_rank = 4; flags = 'polar' !    ! 
        elseif (index(property,'piezo-optic')              .ne.0) then; tensor_rank = 4 !    ! 
        elseif (index(property,'kerr')                     .ne.0) then; tensor_rank = 4 !    ! 
        elseif (index(property,'electrostriction')         .ne.0) then; tensor_rank = 4 !    ! 
        !------------------------------------------------- SXITH-RANK TENSORS -------------------------------------------------
        elseif (index(property,'third-order elasticity')   .ne.0) then; tensor_rank = 6; flags = 'polar' !    ! 
        !----------------------------------------------------------------------------------------------------------------------
        else
            call am_print('ERROR','Unknown property.',flags='E')
            stop
        endif
        !----------------------------------------------------------------------------------------------------------------------
        !
        ! NOTES:
        !
        ! ... Wooten page 449; sflag specifies whether the neighboring indices are symmetric with respect to each other. NEED TO IMPLEMENT STILL. 
        ! N = no intrinsic symmetry on this pair of indices
        ! S = symmetric pair of indices
        ! A = antisymmetric pair of indices
        !
        ! Tensors transferm differently whether they are axial or polar:
        ! POLAR :  T_{i,j,k,l,...} =        sum_{ip,jp,kp,lp,...} R_{i,ip} R_{j,jp} R_{k,kp} R_{l,lp} T_{ip,jp,kp,lp,...}
        ! AXIAL :  T_{i,j,k,l,...} = det(R) sum_{ip,jp,kp,lp,...} R_{i,ip} R_{j,jp} R_{k,kp} R_{l,lp} T_{ip,jp,kp,lp,...}
        ! Thus, all axial tensors of even rank and polar tensors of odd rank are null are null.  Wooten p 485. Eq. 13.21. 
        !
        ! Onsager’s Principle requires that the electric resistivity and thermal conductivity tensors be symmetric.
        ! This does not hold for the Seebeck and Peltier (thermoelectric) tensors which relate two different flows. Thus
        ! there are, at most, nine independent parameters rather than six. [Newnham "Properties of Materials"]
        !
        !----------------------------------------------------------------------------------------------------------------------
        !
        ndim = 3
        nterms = ndim**tensor_rank
        !
        ! get conjugacy class members
        ! members(nclass,maxval(class_nelements))
        class_member = member(pg%cc_identifier)
        !
        ! initialize A. 
        ! using LU factorization instead of applying rref in order to incorporate effect of symmetry on tensor at each step.
        ! LU needs two additional work spaces.
        allocate(A(3*nterms,2*nterms)) ! augmented matrix
        A = 0
        !
        allocate(indices(nterms))
        indices = [1:nterms]
        !
        ! start counter for number of equations
        nequations = 0
        !
        ! crystal symmetries
        !
        do j = 1, size(class_member,1) ! loop over classes
            ! track number of symmetry equations for fun
            nequations = nequations + nterms*count(class_member(j,:).ne.0)
            ! get index of class representative 
            i=class_member(j,1)
            ! construct symmetry operator in the flattend basis
            ! Nye, J.F. "Physical properties of crystals: their representation by tensors and matrices". p 133 Eq 7
            R = kron_pow(ps_frac2cart(R_frac=pg%R(:,:,i),bas=uc%bas),tensor_rank)
            ! apply flip operator. S_{alpha,beta,n,m} = S_{beta,alpha,m,n}
            ! reversal group have operations S which flip n->m and m->n. To return n and m to their correct indicies, apply "flip" operator which flips n<->m AND transposes alpha and beta cartesian indices.
            if (index(property,'pair reversal').ne.0) then
                R = matmul(T(ndim),R)
            endif
            ! if the quantity corresponds to an axial tensor
            if (index(flags,'axial').ne.0) R = det(pg%R(:,:,i)) * R
            ! Save the action of the symmetry operations
            A(2*nterms+indices,0*nterms+indices) = R
            A(2*nterms+indices,1*nterms+indices) = eye(nterms)
            ! At this point A is an augmented matrix of the form
            !
            !        [ A1 , B1 ]
            !   A =  [ A2 , B2 ]
            !        [ A3 , I  ]
            !
            ! in which the square matrix [A1,B1;A2,B2] is an augmented upper triangular matrix corresponding to 2*nterms
            ! coupled equations which describe how the terms are interrelated. The last augmented matrix [A3,I]
            ! describes how the tensor rotates under the most recently considered symmetry operation R. In order
            ! incorporate the effect of symmetry operation R on the set of linearly coupled 2*nterms equations, LU
            ! factorization is performed. By doing so, the augmented matrix A is converted into its factored upper
            ! triangular matrix U, which, by definition, only occupies rows corresponding to the smallest dimension of
            ! the matrix [either row or column, i.e. min(m,n)]. For this case, it occupies the top 2x2 augmented blocks.
            ! This is also why the symmetry equations are added as components in the [A3,I] block -- so they don't
            ! overlap with U.
            !
            ! perform LU factorization (saving only U) to incorporate effect
            call lu(A)
            ! debug flags
            if (opts%verbosity.ge.2) then
                call am_print('ps',ps_frac2cart(R_frac=pg%R(:,:,i),bas=uc%bas),filename='debug_ps'//trim(int2char(i))//'.txt',permission='w')
                call am_print('R',R,filename='debug_R'//trim(int2char(i))//'.txt',permission='w')
                call am_print('A',A,filename='debug_A'//trim(int2char(i))//'.txt',permission='w')
            endif
        enddo
        !
        ! intrinsic symmetries
        !
        allocate(S(nterms,nterms,100))
        k=0
        if     (index(property,'conductivity'    ).ne.0 &
         & .or. index(property,'resistivity'     ).ne.0 &
         & .or. index(property,'voigt'           ).ne.0) then
            k=k+1; S(:,:,k) = T(ndim)                     ! s_ij  = s_ji
        elseif (index(property,'piezoelectricity').ne.0) then
            k=k+1; S(:,:,k) = kron(T(ndim),eye(ndim))     ! d_ijk = d_ikj
        elseif (index(property,'elasticity'      ).ne.0) then
            k=k+1; S(:,:,k) = eye(nterms)                 ! cijkl = cijkl
            k=k+1; S(:,:,k) = kron(eye(ndim**2),T(ndim))  ! cijkl = cjikl
            k=k+1; S(:,:,k) = kron(T(ndim),eye(ndim**2))  ! cijkl = cjilk
            k=k+1; S(:,:,k) = kron(T(ndim),T(ndim))       ! cijkl = cjilk
        endif
        !
        do i = 1, k
            nequations = nequations+nterms
            A(2*nterms+indices,0*nterms+indices) = S(:,:,i)
            A(2*nterms+indices,1*nterms+indices) = eye(nterms)
            call lu(A)
            ! if (opts%verbosity.ge.2) call am_print('I',S(:,:,i),filename='debug_I'//trim(int2char(i))//'.txt',permission='w')
        enddo
        !
        ! Apply Gram-Schmidt orthogonalization by converting A into reduced row echelon form
        !
        call rref(A)
        !
        ! At this point A is an augmented matrix composed of [ LHS | RHS ]. The LHS should be the identity matrix,
        ! which, together with the RHS, completely  specifies all relationships between variables.
        !
        allocate(LHS(nterms,nterms))
        allocate(RHS(nterms,nterms))
        LHS = A(0*nterms+indices,0*nterms+indices)
        RHS = A(0*nterms+indices,1*nterms+indices)
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
        if (present(relations)) allocate(relations, source=RHS)
        !
    end subroutine get_symmetry_adapted_tensor

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
                                write(*,'(a,a,a)',advance='no') trim(dbl2charSP(RHS(i,j),7)), '*a', trim(int2char(j))
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

    ! transpose operator in the flattened basis

    pure function  T(n) result(M)
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
    end function   T

end module am_symmetry_adapted_tensors