module am_symmetry
    !
    use am_constants
    use am_helpers
    use am_unit_cell
    use am_options
    !
    implicit none
    !
    private
    !
    public :: am_class_symmetry
    public :: determine_symmetry
    public :: get_kpoint_compatible_symmetries
    public :: point_symmetry_convert_to_cartesian

    type am_class_symmetry 
        integer :: nsyms                  !> number of point symmetries
        real(dp), allocatable :: R(:,:,:) !> symmetry elements (operate on fractional atomic basis)
        real(dp), allocatable :: T(:,:)   !> symmetry elements (operate on fractional atomic basis)
        !
        integer , allocatable :: schoenflies_code(:) ! look at schoenflies_decode to understand what each integer corresponds to
        character(string_length_schoenflies) , allocatable :: schoenflies(:) ! decoded string
        !
    contains
        ! procedure :: determine_symmetry
        procedure :: point_group
        procedure :: space_group
        procedure :: stdout
        procedure :: tensor_2nd_rank
        procedure :: tensor_3rd_rank
        procedure :: elasticity_tensor
    end type

    interface transform_tensor
        module procedure transform_tensor_fourth_rank, transform_tensor_second_rank, transform_tensor_third_rank
    end interface ! transform_tensor

contains

    !
    ! functions which operate on R(3,3) in cartesian coordinates
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

    function       point_symmetry_schoenflies(R) result(schoenflies_code)
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
        integer :: schoenflies_code
        !
        ! get trace and determinant (fractional)
        tr  = point_symmetry_trace(R)
        det = point_symmetry_determinant(R)
        !
        ! The Mathematical Theory of Symmetry in Solids: Representation Theory for
        ! Point Groups and Space Groups. 1 edition. Oxford?: New York: Oxford University
        ! Press, 2010. page 138, table 3.8.
        !
        schoenflies_code = 0
        if     ( (abs(tr - 3).lt.tiny) .and. (abs(det - 1).lt.tiny) ) then; schoenflies_code = 1  ! 'e'
        elseif ( (abs(tr + 1).lt.tiny) .and. (abs(det - 1).lt.tiny) ) then; schoenflies_code = 2  ! 'c_2'
        elseif ( (abs(tr - 0).lt.tiny) .and. (abs(det - 1).lt.tiny) ) then; schoenflies_code = 3  ! 'c_3'
        elseif ( (abs(tr - 1).lt.tiny) .and. (abs(det - 1).lt.tiny) ) then; schoenflies_code = 4  ! 'c_4'
        elseif ( (abs(tr - 2).lt.tiny) .and. (abs(det - 1).lt.tiny) ) then; schoenflies_code = 5  ! 'c_6'
        elseif ( (abs(tr + 3).lt.tiny) .and. (abs(det + 1).lt.tiny) ) then; schoenflies_code = 6  ! 'i'
        elseif ( (abs(tr - 1).lt.tiny) .and. (abs(det + 1).lt.tiny) ) then; schoenflies_code = 7  ! 's_2'
        elseif ( (abs(tr - 0).lt.tiny) .and. (abs(det + 1).lt.tiny) ) then; schoenflies_code = 8  ! 's_6'
        elseif ( (abs(tr + 1).lt.tiny) .and. (abs(det + 1).lt.tiny) ) then; schoenflies_code = 9  ! 's_4'
        elseif ( (abs(tr + 2).lt.tiny) .and. (abs(det + 1).lt.tiny) ) then; schoenflies_code = 10 ! 's_3'
        endif
        !
        if (schoenflies_code .eq. 0) then
            call am_print('ERROR','Unable to identify point symmetry.',' >>> ')
            call am_print('R (frac)',R)
            call am_print('tr (frac)',tr)
            call am_print('det (frac)',det)
            stop
        endif
        !
    end function   point_symmetry_schoenflies

    !
    ! functions which operate on tensors or make symmetry-adapted tensors
    !

    subroutine     tensor_2nd_rank(pg,uc,opts)
        !
        class(am_class_symmetry), intent(in) :: pg
        type(am_class_unit_cell), intent(in) :: uc
        type(am_class_options), intent(in) :: opts
        !
        real(dp) :: M(3,3)                          !> IMPORTANT: M is the original tensor before the permutation 
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
        if (opts%verbosity.ge.1) call am_print_title('Determining second-rank tensor')
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
        call am_print_double_matrix(A=A,iopt_fname='outfile.sym.A_prerref')
        call rref(A)
        !
        call print_sparse('A = [LHS|RHS]',A(1:nterms,:),' ... ')
        call am_print_double_matrix(A=A,iopt_fname='outfile.sym.A_postrref')
        ! At this point A is an augmented matrix composed of [ LHS | RHS ]. The LHS
        ! should be the identity matrix, which, together with the RHS, completely 
        ! specifies all relationships between variables. 
        allocate(LHS(nterms,nterms))
        allocate(RHS(nterms,nterms))
        LHS = A(1:nterms,1:nterms)
        RHS = A(1:nterms,[1:nterms]+nterms)
        !
        call am_print_double_matrix(A=LHS,iopt_fname='outfile.sym.LHS')
        call am_print_double_matrix(A=RHS,iopt_fname='outfile.sym.RHS')
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
    end subroutine tensor_2nd_rank

    subroutine     tensor_3rd_rank(pg,uc,opts)
        !
        class(am_class_symmetry), intent(in) :: pg
        type(am_class_unit_cell), intent(in) :: uc
        type(am_class_options), intent(in) :: opts
        !
        real(dp) :: M(3,3,3) !> IMPORTANT: M is the original tensor before the permutation 
#include "am_symmetry_tensor.f90"
    end subroutine tensor_3rd_rank

    subroutine     elasticity_tensor(pg,uc,opts)
        !
        class(am_class_symmetry), intent(in) :: pg
        type(am_class_unit_cell), intent(in) :: uc
        type(am_class_options), intent(in) :: opts
        !
        real(dp) :: M(3,3,3,3)                      !> IMPORTANT: M is the original tensor before the permutation 
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
    end subroutine elasticity_tensor

    pure function  transform_tensor_second_rank(M,R) resulT(RM)
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

    pure function  transform_tensor_third_rank(M,R) resulT(RM)
        !
        implicit none
        !
        real(dp), intent(in) :: M(:,:,:)
        real(dp), intent(in) :: R(3,3)
        real(dp) :: RM(size(M,1),size(M,2),size(M,3))
        integer :: i,j,k,ip,jp,kp
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

    pure function  transform_tensor_fourth_rank(M,R) resulT(RM)
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
            write(*,'(a5,a)') ' ... ', 'symmetry equations'
            !
            ! write the terms that are independent (equal only to themselves)
            !
            do i = 1,nterms
                if (is_independent(i)) then
                    write(*,'(5x,a1,a,a1,a1,a)',advance='no') 'a',trim(int2char(i)),'=', 'a', trim(int2char(i))
                    write(*,*)
                endif
            enddo
            !
            ! write the terms that are interrelated.
            !
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
            !
            ! finally, write the terms that are equal to zero.
            !
            do i = 1,nterms
                if (is_zero(i)) then
                    write(*,'(5x,a1,a,a)') 'a',trim(int2char(i)),'=+0.00'
                endif
            enddo
        endif
        !
    end subroutine parse_symmetry_equations

    !
    ! schoenflies decode
    !

    function       schoenflies_decode(schoenflies_code) result(schoenflies)
        !
        implicit none
        !
        integer, intent(in) :: schoenflies_code
        character(len=string_length_schoenflies) :: schoenflies
        !
        select case (schoenflies_code)
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
    end function   schoenflies_decode

    function       point_group_schoenflies_decode(pg_schoenflies_code) result(pg_schoenflies)
        !
        implicit none
        !
        integer, intent(in) :: pg_schoenflies_code
        character(string_length_schoenflies) :: pg_schoenflies
        !
        select case (pg_schoenflies_code)
                case(1); pg_schoenflies = 'c_1'
                case(2); pg_schoenflies = 's_2'
                case(3); pg_schoenflies = 'c_2'
                case(4); pg_schoenflies = 'c_1h'
                case(5); pg_schoenflies = 'c_2h'
                case(6); pg_schoenflies = 'd_2'
                case(7); pg_schoenflies = 'c_2v'
                case(8); pg_schoenflies = 'd_2h'
                case(9); pg_schoenflies = 'c_3'
                case(10); pg_schoenflies = 's_6'
                case(11); pg_schoenflies = 'd_3'
                case(12); pg_schoenflies = 'c_3v'
                case(13); pg_schoenflies = 'd_3d'
                case(14); pg_schoenflies = 'c_4'
                case(15); pg_schoenflies = 's_4'
                case(16); pg_schoenflies = 'c_4h'
                case(17); pg_schoenflies = 'd_4'
                case(18); pg_schoenflies = 'c_4v'
                case(19); pg_schoenflies = 'd_2d'
                case(20); pg_schoenflies = 'd_4h'
                case(21); pg_schoenflies = 'c_6'
                case(22); pg_schoenflies = 'c_3h'
                case(23); pg_schoenflies = 'c_6h'
                case(24); pg_schoenflies = 'd_6'
                case(25); pg_schoenflies = 'c_6v'
                case(26); pg_schoenflies = 'd_3h'
                case(27); pg_schoenflies = 'd_6h'
                case(28); pg_schoenflies = 't'
                case(29); pg_schoenflies = 't_h'
                case(30); pg_schoenflies = 'o'
                case(31); pg_schoenflies = 't_d'
                case(32); pg_schoenflies = 'o_h'
            case default
                call am_print('ERROR','Schoenflies code for point-group unknown.',' >>> ')
                call am_print('pg_schoenflies_code',pg_schoenflies_code)
                stop
            end select
    end function   point_group_schoenflies_decode

    function       point_group_schoenflies(schoenflies_code) result(pg_schoenflies_code)
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
        !>   pg_schoenflies_code --> pg_schoenflies
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
        integer, intent(in) :: schoenflies_code(:)
        integer :: nsyms
        integer :: ni, nc2, nc3, nc4, nc6, ns2, ns6, ns4, ns3, ne
        integer :: pg_schoenflies_code
        !
        pg_schoenflies_code = 0
        !
        ! trivial cases first
        !
        nsyms = size(schoenflies_code)
        if     (nsyms .eq. 1 ) then; pg_schoenflies_code=1;  return
        elseif (nsyms .eq. 48) then; pg_schoenflies_code=32; return
        elseif (nsyms .eq. 16) then; pg_schoenflies_code=20; return
        elseif (nsyms .eq. 3 ) then; pg_schoenflies_code=9;  return
        endif
        !
        ni=0; nc2=0; nc3=0; nc4=0; nc6=0; ns2=0; ns6=0; ns4=0; ns3=0
        ne  = count(schoenflies_code.eq.1)  ! 'e' identity
        nc2 = count(schoenflies_code.eq.2)  ! 'c_2'
        nc3 = count(schoenflies_code.eq.3)  ! 'c_3'
        nc4 = count(schoenflies_code.eq.4)  ! 'c_4'
        nc6 = count(schoenflies_code.eq.5)  ! 'c_6'
        ni  = count(schoenflies_code.eq.6)  ! 'i' inversion
        ns2 = count(schoenflies_code.eq.7)  ! 's_2'
        ns6 = count(schoenflies_code.eq.8)  ! 's_6'
        ns4 = count(schoenflies_code.eq.9)  ! 's_4'
        ns3 = count(schoenflies_code.eq.10) ! 's_3'
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
            if (ni.eq.1)  then; pg_schoenflies_code=2;  return; endif
            if (nc2.eq.1) then; pg_schoenflies_code=3;  return; endif
            if (ns2.eq.1) then; pg_schoenflies_code=4;  return; endif
        endif
        if (nsyms.eq.4)   then
            if (ni.eq.1)  then; pg_schoenflies_code=5;  return; endif
            if (nc2.eq.3) then; pg_schoenflies_code=6;  return; endif
            if (ns2.eq.2) then; pg_schoenflies_code=7;  return; endif
            if (nc4.eq.1) then; pg_schoenflies_code=14; return; endif
            if (ns4.eq.2) then; pg_schoenflies_code=15; return; endif
        endif
        if (nsyms.eq.6)   then
            if (ni.eq.1)  then; pg_schoenflies_code=10; return; endif
            if (nc2.eq.3) then; pg_schoenflies_code=11; return; endif
            if (ns2.eq.3) then; pg_schoenflies_code=12; return; endif
            if (nc2.eq.1) then; pg_schoenflies_code=21; return; endif
            if (ns2.eq.1) then; pg_schoenflies_code=22; return; endif
        endif
        if (nsyms.eq.8)  then
            if (ns2.eq.3) then; pg_schoenflies_code=8;  return; endif
            if (ns2.eq.1) then; pg_schoenflies_code=16; return; endif
            if (ns2.eq.0) then; pg_schoenflies_code=17; return; endif
            if (ns2.eq.4) then; pg_schoenflies_code=18; return; endif
            if (ns2.eq.2) then; pg_schoenflies_code=19; return; endif
        endif
        if (nsyms.eq.12)  then
            if (ns2.eq.3) then; pg_schoenflies_code=13; return; endif
            if (ns2.eq.1) then; pg_schoenflies_code=23; return; endif
            if (nc2.eq.7) then; pg_schoenflies_code=24; return; endif
            if (ns2.eq.6) then; pg_schoenflies_code=25; return; endif
            if (ns2.eq.4) then; pg_schoenflies_code=26; return; endif
            if (nc3.eq.8) then; pg_schoenflies_code=28; return; endif
        endif
        if (nsyms.eq.24)  then
            if (nc6.eq.2) then; pg_schoenflies_code=27; return; endif
            if (ni.eq.1)  then; pg_schoenflies_code=29; return; endif
            if (nc4.eq.6) then; pg_schoenflies_code=30; return; endif
            if (ns4.eq.6) then; pg_schoenflies_code=31; return; endif
        endif
    end function   point_group_schoenflies

    !
    ! functions which operate on ss
    !

    subroutine     point_group(pg,uc,ss,opts)
        !
        ! requires only ss%R
        !
        use am_rank_and_sort
        !
        implicit none
        !
        class(am_class_symmetry), intent(inout) :: pg
        type(am_class_symmetry) , intent(in) :: ss
        type(am_class_unit_cell), intent(in) :: uc
        type(am_class_options), intent(in) :: opts
        integer :: pg_schoenflies_code
        character(string_length_schoenflies) :: pg_schoenflies
        integer , allocatable :: sorted_indices(:)
        real(dp), allocatable :: sort_parameter(:)
        integer :: i
        !
        if (opts%verbosity.ge.1) call am_print_title('Determining point group symmetries')
        !
        pg%R = unique(ss%R)
        !
        pg%nsyms = size(pg%R,3)
        !
        ! <sort>
        !
        allocate(sorted_indices(pg%nsyms))
        allocate(sort_parameter(pg%nsyms))
        !
        do i = 1, pg%nsyms
            sort_parameter(i) = point_symmetry_determinant(pg%R(:,:,i))
        enddo
        call rank(sort_parameter,sorted_indices)
        sorted_indices = sorted_indices(pg%nsyms:1:-1)
        pg%R=pg%R(:,:,sorted_indices)
        !
        do i = 1, pg%nsyms
            sort_parameter(i) = point_symmetry_trace(pg%R(:,:,i))
        enddo
        call rank(sort_parameter,sorted_indices)
        sorted_indices = sorted_indices(pg%nsyms:1:-1)
        pg%R=pg%R(:,:,sorted_indices)
        !
        ! </sort>
        !
        allocate(pg%schoenflies_code(pg%nsyms))
        allocate(pg%schoenflies(pg%nsyms))
        do i = 1, pg%nsyms
            pg%schoenflies_code(i) = point_symmetry_schoenflies(pg%R(:,:,i))
            pg%schoenflies(i)      = schoenflies_decode(pg%schoenflies_code(i))
        enddo
        !
        pg_schoenflies_code = point_group_schoenflies(pg%schoenflies_code)
        pg_schoenflies      = point_group_schoenflies_decode(pg_schoenflies_code)
        !
        if (opts%verbosity.ge.1) call am_print('point group',pg_schoenflies,' ... ')
        !
        if (opts%verbosity.ge.1) call am_print('number of point symmetries',pg%nsyms,' ... ')
        !
        if (opts%verbosity.ge.1) then
            call am_print('number of e   point symmetries',count(pg%schoenflies_code.eq.1),' ... ')
            call am_print('number of c_2 point symmetries',count(pg%schoenflies_code.eq.2),' ... ')
            call am_print('number of c_3 point symmetries',count(pg%schoenflies_code.eq.3),' ... ')
            call am_print('number of c_4 point symmetries',count(pg%schoenflies_code.eq.4),' ... ')
            call am_print('number of c_6 point symmetries',count(pg%schoenflies_code.eq.5),' ... ')
            call am_print('number of i   point symmetries',count(pg%schoenflies_code.eq.6),' ... ')
            call am_print('number of s_2 point symmetries',count(pg%schoenflies_code.eq.7),' ... ')
            call am_print('number of s_6 point symmetries',count(pg%schoenflies_code.eq.8),' ... ')
            call am_print('number of s_4 point symmetries',count(pg%schoenflies_code.eq.9),' ... ')
            call am_print('number of s_3 point symmetries',count(pg%schoenflies_code.eq.10),' ... ')
        endif
        !
        if (opts%verbosity.ge.1) call pg%stdout(iopt_uc=uc)
        !
        call pg%stdout(iopt_uc=uc,iopt_filename=trim('outfile.pointgroup'))
        !
    end subroutine point_group

    subroutine     space_group(ss,uc,opts)
        !
        use am_rank_and_sort
        !
        implicit none
        !
        class(am_class_symmetry), intent(inout) :: ss
        type(am_class_unit_cell), intent(in) :: uc ! space groups are tabulated in the litearture for conventional cells, but primitive or arbitrary cell works just as wlel.
        type(am_class_options), intent(in) :: opts
        real(dp), allocatable :: seitz(:,:,:)
        integer , allocatable :: sorted_indices(:)
        real(dp), allocatable :: sort_parameter(:)
        integer :: i, j
        !
        if (opts%verbosity.ge.1) call am_print_title('Determining space group symmetries')
        !
        seitz = space_symmetries_from_basis(uc=uc,opts=opts)
        !
        ss%nsyms=size(seitz,3)
        !
        if (opts%verbosity.ge.1) call am_print('number of space group symmetries found',ss%nsyms,' ... ') 
        !
        allocate(ss%R(3,3,ss%nsyms))
        allocate(ss%T(3,ss%nsyms))
        ss%R(1:3,1:3,1:ss%nsyms) = seitz(1:3,1:3,1:ss%nsyms)
        ss%T(1:3,1:ss%nsyms)     = seitz(1:3,4,1:ss%nsyms)
        !
        ! <sort>
        !
        allocate(sorted_indices(ss%nsyms))
        allocate(sort_parameter(ss%nsyms))
        !
        do i = 1, ss%nsyms; sort_parameter(i) = point_symmetry_determinant(ss%R(:,:,i)); enddo
        call rank(sort_parameter,sorted_indices)
        sorted_indices = sorted_indices(ss%nsyms:1:-1)
        ss%R=ss%R(:,:,sorted_indices)
        ss%T=ss%T(:,sorted_indices)
        !
        do i = 1, ss%nsyms; sort_parameter(i) = point_symmetry_trace(ss%R(:,:,i)); enddo
        call rank(sort_parameter,sorted_indices)
        sorted_indices = sorted_indices(ss%nsyms:1:-1)
        ss%R=ss%R(:,:,sorted_indices)
        ss%T=ss%T(:,sorted_indices)
        !
        do j = 1,3
            sort_parameter = ss%T(j,:)
            call rank(sort_parameter,sorted_indices)
            ss%R=ss%R(:,:,sorted_indices)
            ss%T=ss%T(:,sorted_indices)
        enddo
        !
        ! </sort>
        !
        if (opts%verbosity.ge.1) call ss%stdout(iopt_uc=uc)
        !
        call ss%stdout(iopt_uc=uc,iopt_filename=trim('outfile.spacegroup'))
        !
    end subroutine space_group

    subroutine     stdout(ss,iopt_uc,iopt_filename)
        !
        implicit none
        !
        class(am_class_symmetry), intent(in) ::  ss
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
        if (allocated(ss%T)) width=width + 3*7+4
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
        if (allocated(ss%T))  write(fid,fmt4,advance='no') '  | ', 'T1', 'T2', 'T3'
        write(fid,*)
        !
        do i = 1, ss%nsyms
            R = ss%R(:,:,i)
            write(fid,fmt1,advance='no') i, (( trim(print_pretty(R(j,k))), j = 1,3) , k = 1,3)
            if (allocated(ss%T)) then
                T = ss%T(:,i)
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
        !
        if (present(iopt_uc)) then
        !
        write(fid,fmt3,advance='no') centertitle('cartesian',width)
        write(fid,*)
        write(fid,fmt3,advance='no') repeat('-',width)
        write(fid,*)
        !
        write(fid,fmt2,advance='no') '#', 'R11', 'R21', 'R31', 'R12', 'R22', 'R32', 'R13', 'R23', 'R33'
        if (allocated(ss%T))  write(fid,fmt4,advance='no') '  | ', 'T1', 'T2', 'T3'
        write(fid,*)
        !
        do i = 1, ss%nsyms
            R = point_symmetry_convert_to_cartesian(R_frac=ss%R(:,:,i),bas=iopt_uc%bas)
            write(fid,fmt1,advance='no') i, (( trim(print_pretty(R(j,k))), j = 1,3) , k = 1,3)
            if (allocated(ss%T)) then
                T = matmul(iopt_uc%bas,ss%T(:,i))
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

    subroutine     determine_symmetry(uc,ss,pg,opts)
        !
        implicit none
        !
        type(am_class_symmetry), intent(out) :: ss
        type(am_class_symmetry), intent(out) :: pg
        type(am_class_unit_cell), intent(in) :: uc
        type(am_class_options)  , intent(in) :: opts
        !
        if (opts%verbosity.ge.1) call am_print_title('Analyzing symmetry')
        !
        call ss%space_group(uc=uc,opts=opts)
        !
        call pg%point_group(ss=ss,uc=uc,opts=opts)
        !
        !
    end subroutine determine_symmetry

    !
    ! functions which operate on kpoints
    !

    function        get_kpoint_compatible_symmetries(kpt,R,sym_prec) result(R_syms)
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
    end function    get_kpoint_compatible_symmetries

end module am_symmetry


