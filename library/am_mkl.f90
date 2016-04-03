module am_mkl
    !
    use am_constants
    !
    integer, parameter :: lwmax = 1000
    !
    public

contains

    ! inverse of square matrix
    
    function       inv(A) result(Ainv)
        !
        use am_helpers
        !
        implicit none
        !
        real(dp), intent(in)  :: A(:,:)
        integer,  allocatable :: ipiv(:)
        real(dp), allocatable :: Ainv(:,:)
        real(dp), allocatable :: WORK(:)
        integer  :: m, n, lda, info, lwork
        !
        m = size(A,1)
        n = size(A,2)
        allocate(ipiv(n))
        allocate(Ainv(n,n))
        allocate(WORK(lwmax))
        Ainv = A
        lda  = m
        !
        ! factorize 
        call dgetrf( m, n, ainv, lda, ipiv, info )
        if (info.ne.0) then
            call am_print('ERROR',info)
            stop
        endif
        !
        ! querery workspace
        lwork=-1
        call dgetri( n, ainv, lda, ipiv, WORK, lwork, info )
        lwork = min( lwmax, int( WORK( 1 ) ) )
        !
        ! get inverse
        !
        call dgetri( n, ainv, lda, ipiv, WORK, lwork, info )
        if (info.ne.0) then
            call am_print('ERROR',info)
            stop
        endif
        !
    end function   inv

    function       pinv(A) result(Ainv)
        !
        implicit none
        !
        real(dp), intent(in)  :: A(:,:)
        real(dp), allocatable :: B(:,:)
        real(dp), allocatable :: Ainv(:,:)
        real(dp), allocatable :: S(:) ! singular values of A
        real(dp), allocatable :: C(:,:)
        real(dp), allocatable :: WORK(:)
        real(dp) :: rcond
        integer :: m, n, sz, i, j
        integer :: lda, ldb, nrhs
        integer :: info
        integer :: lwork
        integer :: rank ! effective rank of A: the number of singular values which are greater than rcond*s(1).
        !
        rcond = 1.0D-8
        !
        m = size(A,1)
        n = size(A,2)
        lda  = m
        !
        allocate(C,source=A)
        allocate(WORK(lwmax))
        allocate(S(max(1,min(m,n))))
        !
        ! copy variables
        !
        sz = max(m,n)
        allocate(B(sz,sz))
        do i = 1, sz
        do j = 1, sz
            if (i.eq.j) then
                B(i,j) = 1.0_dp
            else
                B(i,j) = 0.0_dp
            endif
        enddo
        enddo
        ldb  = sz
        nrhs = sz
        !
        ! query the optimal workspace
        lwork = -1
        call dgelss(m, n, nrhs, C, lda, B, ldb, S, rcond, rank, WORK, lwork, info)
        lwork = min(lwmax,int(WORK(1)))
        !
        ! solve the least squares problem min( norm2(B - Ax) ) for the x of minimum norm.
        call dgelss(m, n, nrhs, C, lda, B, ldb, S, rcond, rank, WORK, lwork, info)
        !
        ! check for convergence
        if (info.gt.0) then
            write(*,*) 'algorithm failed', info
            stop
        end if
        !
        allocate(Ainv,source=B(1:n,1:m))
        !
    end function   pinv

    ! diagonalize real symmetric square matrix

    subroutine     am_dsyev(A,V,D)
        !
        implicit none
        !
        real(dp), intent(in)  :: A(:,:)
        complex(dp), intent(out), allocatable :: D(:)   ! eigenvalues of the matrix A in ascending order
        real(dp), intent(out), allocatable :: V(:,:) ! orthonormal eigenvectors of the matrix A
        real(dp), allocatable :: WORK(:)
        integer :: n
        integer :: lda
        integer :: info
        integer :: lwork
        !
        ! initialize variables
        !
        lda = size(A,1)
        n   = size(A,2)
        allocate(D(n))
        allocate(V(lda,n))
        V = A
        !
        ! query the optimal workspace
        !
        allocate(WORK(lwmax))
        lwork = -1
        call dsyev( 'vectors', 'upper', n, V, lda, D, WORK, lwork, info )
        lwork = min( lwmax, int( WORK( 1 ) ) )
        !
        ! solve eigenproblem
        !
        call dsyev( 'vectors', 'upper', n, V, lda, D, WORK, lwork, info )
        !
        ! check for convergence.
        !
        if( info.gt.0 ) then
            write(*,*) 'the algorithm failed to compute eigenvalues.'
            stop
        end if
    end subroutine am_dsyev
    
    ! diagonalize real general square matrix

    subroutine     am_dgeev(A,VL,D,VR)
        !
        ! transpose(VL)*A*VR = D
        ! this may be so, but VL bares no obvious relationship to VR. For example, VL is not the inverse of VR! 
        ! 
        implicit none
        !
        real(dp)   , intent(in) :: A(:,:)
        real(dp)   , intent(out), allocatable, optional :: VR(:,:), VL(:,:) ! right and left eigenvectors
        complex(dp), intent(out), allocatable, optional :: D(:)
        real(dp), allocatable :: CA(:,:)
        real(dp), allocatable :: WR(:), WI(:), WORK(:)
        real(dp), allocatable :: VL_internal(:,:), VR_internal(:,:)
        integer :: i
        integer :: n
        integer :: lda, ldvl, ldvr
        integer :: info
        integer :: lwork
        !
        lda = size(A,1)
        n   = size(A,2)
        ldvl = n
        ldvr = n
        !
        allocate(VL_internal(ldvl,n))
        allocate(VR_internal(ldvr,n))
        allocate(WR(n))
        allocate(WI(n))
        allocate(WORK(lwmax))
        !
        ! copy variables
        !
        allocate(CA, source=A)
        !
        ! query the optimal workspace
        !
        lwork = -1
        call dgeev( 'N', 'V', n, CA, lda, WR, WI, VL_internal, ldvl, VR_internal, ldvr, WORK, lwork, info )
        lwork = min( lwmax, int( WORK( 1 ) ) )
        !
        ! solve eigenproblem
        !
        if (present(VL).and.present(VR)) then
            call dgeev( 'V', 'V', n, CA, lda, WR, WI, VL_internal, ldvl, VR_internal, ldvr, WORK, lwork, info )
        elseif (present(VR)) then
            call dgeev( 'N', 'V', n, CA, lda, WR, WI, VL_internal, ldvl, VR_internal, ldvr, WORK, lwork, info )
        elseif (present(VL)) then
            call dgeev( 'V', 'N', n, CA, lda, WR, WI, VL_internal, ldvl, VR_internal, ldvr, WORK, lwork, info )
        endif
        !
        ! check for convergence
        !
        if( info.gt.0 ) then
            write(*,*)'the algorithm failed to compute eigenvalues.'
            stop
        end if
        !
        if (present(D)) then
            allocate(D(n))
            do i = 1, n
                D(i) = cmplx(WR(i),WI(i))
            enddo
        endif
        !
        if (present(VR)) then
            allocate(VR,source=VR_internal)
        endif
        !
        if (present(VL)) then
            allocate(VL,source=VL_internal)
        endif
    end subroutine am_dgeev
    
    ! diagonalize complex hermitian matrix

    subroutine     am_zheev(A,V,D)
        !
        implicit none
        !
        complex(dp), intent(in)  :: A(:,:)
        complex(dp), intent(out), allocatable :: D(:)   ! eigenvalues of the matrix A in ascending order
        complex(dp), intent(out), allocatable :: V(:,:) ! orthonormal eigenvectors of the matrix A
        complex(dp), allocatable :: WORK(:)
        real(dp)   , allocatable :: RWORK(:)
        integer :: n
        integer :: lda
        integer :: info
        integer :: lwork
        !
        ! initialize variables
        !
        lda = size(A,1)
        n   = size(A,2)
        allocate(D(n))
        allocate(RWORK(3*n-2))
        !
        ! copy variables
        !

        allocate(V, source=A)
        !
        ! query the optimal workspace
        !
        allocate(WORK(lwmax))
        lwork = -1
        call zheev( 'vectors', 'upper', n, V, lda, D, WORK, lwork, RWORK, info )
        lwork = min( lwmax, int( WORK( 1 ) ) )
        !
        ! solve eigenproblem
        !
        call zheev( 'vectors', 'upper', n, V, lda, D, WORK, lwork, RWORK, info )
        !
        ! check for convergence.
        !
        if( info.gt.0 ) then
            write(*,*) 'the algorithm failed to compute eigenvalues.'
            stop
        end if
    end subroutine am_zheev

    ! perform svd decomposition on real general rectangular matrix

    subroutine     am_dgesvd(A,U,S,VT)
        !
        !  The routine computes the singular value decomposition (SVD) of a real
        !  m-by-n matrix A, optionally computing the left and/or right singular
        !  vectors. The SVD is written as
        !
        !  A = U*S*VT
        !
        !  where S is an m-by-n matrix which is zero except for its min(m,n)
        !  diagonal elements, U is an m-by-m orthogonal matrix and VT (V transposed)
        !  is an n-by-n orthogonal matrix. The diagonal elements of S
        !  are the singular values of A; they are real and non-negative, and are
        !  returned in descending order. The first min(m, n) columns of U and V are
        !  the left and right singular vectors of A.
        !
        !  Note that the routine returns VT, not V.
        real(dp), intent(in)  :: A(:,:)
        real(dp), intent(out), allocatable :: U(:,:)
        real(dp), intent(out), allocatable :: S(:)
        real(dp), intent(out), allocatable :: VT(:,:)
        real(dp), allocatable :: CA(:,:)
        real(dp), allocatable :: WORK(:)
        integer :: m, n
        integer :: lda, ldu, ldvt
        integer :: info
        integer :: lwork
        !
        m = size(A,1)
        n = size(A,2)
        lda  = m
        ldu  = m
        ldvt = n
        !
        allocate(U(ldu,m))
        allocate(VT(ldvt,n))
        allocate(S(n))
        allocate(WORK(lwmax))
        !
        ! copy variables
        !
        allocate(CA, source=A)
        !
        ! query the optimal workspace
        !
        lwork = -1
        call dgesvd( 'all', 'all', m, n, CA, lda, S, U, ldu, VT, ldvt, WORK, lwork, info )
        lwork = min( lwmax, int( WORK( 1 ) ) )
        !
        ! solve eigenproblem
        !
        call dgesvd( 'all', 'all', m, n, CA, lda, S, U, ldu, VT, ldvt, WORK, lwork, info )
        !
        ! check for convergence
        !
        if( info.gt.0 ) then
            write(*,*) 'the algorithm computing svd failed to converge.'
            stop
        end if
    end subroutine am_dgesvd

    ! solve Ax = B by QR factorization

    subroutine     am_dgels(A,X,B)
        ! Program computes the least squares solution to the overdetermined linear
        ! system A*X = B with full rank matrix A using QR factorization,
        !
        ! The routine solves overdetermined or underdetermined real linear systems
        ! involving an m-by-n matrix A, or its transpose, using a QR or LQ
        ! factorization of A. It is assumed that A has full rank.
        !
        ! Several right hand side vectors b and solution vectors x can be handled
        ! in a single call; they are stored as the columns of the m-by-nrhs right
        ! hand side matrix B and the n-by-nrhs solution matrix X.
        !
        implicit none
        !
        real(dp), intent(in)  :: A(:,:)
        real(dp), intent(in)  :: B(:,:)
        real(dp), intent(out), allocatable :: X(:,:)
        real(dp), allocatable :: CA(:,:)
        real(dp), allocatable :: WORK(:)
        integer :: m, n, nrhs
        integer :: lda, ldb
        integer :: info
        integer :: lwork
        !
        ! initialize variables
        !
        lda = size(A,1)
        n   = size(A,2)
        ldb = size(B,1)
        nrhs= size(B,2)
        allocate(WORK(lwmax))
        !
        ! copy variables
        !
        allocate(CA, source = A)
        allocate(X , source = B)
        !
        ! query the optimal workspace
        !
        allocate(WORK(lwmax))
        lwork = -1
        call dgels( 'No transpose', m, n, nrhs, CA, lda, X, ldb, WORK, lwork, info )
        lwork = min( lwmax, int( WORK( 1 ) ) )
        !
        ! call solver
        !
        call dgels( 'No transpose', m, n, nrhs, CA, lda, X, ldb, WORK, lwork, info )
        !
        ! check for convergence.
        !
        if( info.gt.0 ) then
            write(*,*) 'The diagonal element ',info,' of the triangular '
            write(*,*) 'factor of a is zero, so that a does not have full '
            write(*,*) 'rank; the least squares solution could not be '
            write(*,*) 'computed.'
            stop
        end if
    end subroutine am_dgels

    ! solve Ax = B by SVD decomposition

    subroutine     am_dgelss(A,X,B)
        !
        real(dp), intent(in) :: A(:,:)
        real(dp), intent(in) :: B(:,:)
        real(dp), intent(out), allocatable :: X(:,:)
        real(dp), allocatable :: S(:) ! singular values of A
        real(dp), allocatable :: CA(:,:)
        real(dp), allocatable :: WORK(:)
        real(dp) :: rcond
        integer :: m, n
        integer :: lda, ldb, nrhs
        integer :: info
        integer :: lwork
        integer :: rank ! effective rank of A: the number of singular values which are greater than rcond*s(1).
        !
        rcond = 1.0D-8
        !
        m = size(A,1)
        n = size(A,2)
        ldb  = size(B,1)
        nrhs = size(B,2)
        lda  = m
        !
        allocate(WORK(lwmax))
        allocate(S(max(1, min(m, n))))
        !
        ! copy variables
        !
        allocate(CA, source = A)
        allocate(X , source = B)
        !
        ! query the optimal workspace
        !
        lwork = -1
        call dgelss(m, n, nrhs, CA, lda, X, ldb, S, rcond, rank, WORK, lwork, info)
        lwork = min( lwmax, int( WORK( 1 ) ) )
        !
        ! solve the least squares problem min( norm2(B - Ax) ) for the x of minimum norm.
        !
        call dgelss(m, n, nrhs, CA, lda, X, ldb, S, rcond, rank, WORK, lwork, info)
        !
        ! check for convergence
        !
        if( info.gt.0 ) then
            write(*,*) 'the algorithm computing svd failed to converge.'
            stop
        end if
    end subroutine am_dgelss


end module