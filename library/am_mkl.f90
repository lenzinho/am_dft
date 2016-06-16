module am_mkl
    
    use am_constants
    
    integer , parameter :: lwmax = 1000  ! maximum workspace size
    real(dp), parameter :: eps = 1.0D-14 ! used for regularization
    
    public
    
    private :: lwmax
    private :: eps
    
    interface norm
        module procedure am_dznrm2, am_dnrm2    
    end interface ! norm

    interface dot
        module procedure am_zdotc, am_ddot
    end interface ! dot

    interface det
        module procedure zdet, ddet
    end interface ! det

    interface inv
        module procedure am_zinv, am_dinv
    end interface ! det

    interface rand
        module procedure :: d_rand, dv_rand, dm_rand, dt_rand
    end interface ! rand

    interface null_space
        module procedure null_space_zgesvd, null_space_dgesvd
    end interface ! null_space

contains

    ! random number

    function       d_rand() result(a)
        !
        use ifport, only : rand_interal => rand
        !
        implicit none
        !
        real(dp) :: a
        !
        a = rand_interal()
        !
    end function   d_rand

    function       dv_rand(n) result(a)
        !
        use ifport, only : rand_interal => rand
        !
        implicit none
        !
        integer, intent(in) :: n
        real(dp) :: a(n)
        integer :: i
        !
        do i = 1,n
            a(i) = rand_interal()
        enddo
        !
    end function   dv_rand

    function       dm_rand(n,m) result(a)
        !
        use ifport, only : rand_interal => rand
        !
        implicit none
        !
        integer, intent(in) :: n,m
        real(dp) :: a(n,m)
        integer :: i,j
        !
        do i = 1,n
        do j = 1,m
            a(i,j) = rand_interal()
        enddo
        enddo
        !
    end function   dm_rand

    function       dt_rand(n,m,o) result(a)
        !
        use ifport, only : rand_interal => rand
        !
        implicit none
        !
        integer, intent(in) :: n,m,o
        real(dp) :: a(n,m,o)
        integer :: i,j,k
        !
        do i = 1,n
        do j = 1,m
        do k = 1,o
            a(i,j,k) = rand_interal()
        enddo
        enddo
        enddo
        !
    end function   dt_rand

    ! eucledian norm of vector

    function       am_dnrm2(V) result(n)
        !
        implicit none
        !
        real(dp), intent(in) :: V(:)
        real(dp) :: n
        real(dp), external :: dnrm2
        !
        n = dnrm2(size(V), V, 1)
        !
    end function   am_dnrm2 

    function       am_dznrm2(V) result(n)
        !
        implicit none
        !
        complex(dp), intent(in) :: V(:)
        real(dp) :: n
        real(dp), external :: dznrm2
        !
        n = dznrm2(size(V), V, 1)
        !
    end function   am_dznrm2

    ! dot products

    function       am_zdotc(X,Y) result(res)
        !
        ! res = dot(conjg(x), y)
        ! 
        implicit none
        !
        complex(dp), intent(in) :: X(:),Y(:)
        real(dp) :: res
        complex(dp), external :: zdotc
        !
        res = zdotc(size(X), X, 1, Y, 1)
        !
    end function   am_zdotc

    function       am_ddot(X,Y) result(res)
        !
        ! res = dot(x, y)
        ! 
        implicit none
        !
        real(dp), intent(in) :: X(:),Y(:)
        real(dp) :: res
        real(dp), external :: ddot
        !
        res = ddot(size(X), X, 1, Y, 1)
        !
    end function   am_ddot

    ! sort array in increasing order
    
    function       sort(A,flags) result(d)
        !
        implicit none
        !
        real(dp), intent(in) :: A(:)
        character(*), intent(in), optional :: flags
        real(dp), allocatable :: d(:)
        character(1) :: id
        integer :: info
        !
        ! copy input 
        allocate(d,source=A)
        ID = 'I'
        if (present(flags)) then
            if     (index(flags,'ascend')) then 
                id = 'I'
            elseif (index(flags,'descend')) then
                id = 'D'
            else
                stop 'ERROR [sort]: flags /= I or D'
            endif
        endif
        !
        call dlasrt( id, size(d), d, info )
        !
    end function   sort 
    
    ! get determinant by LU facorization

    function       zdet(A) result(d)
        ! 
        implicit none
        !
        complex(dp), intent(in)  :: A(:,:)
        complex(dp) :: d
        complex(dp), allocatable :: C(:,:)
        integer , allocatable :: ipiv(:)
        integer :: pivots
        integer :: m, n, info
        integer :: i 
        !
        !
        m = size(A,1)
        n = size(A,2)
        !
        if (m/=n) stop 'Determinant are only defined for square matrices'

        allocate(ipiv(m))
        !
        ! remove numerical singularites by regularization
        allocate(C,source=A)
        do i = 1, m
            C(i,i) = C(i,i) + eps
        enddo
        !
        ! perform LU factorization
        call zgetrf( m, m, C, m, ipiv, info )
        !
        if (info/=0) stop 'Error in zgetrf'
        !
        ! get determinant
        d=1.0_dp
        do i = 1, m
            d = d*C(i,i)
        enddo
        !
        ! correct for number of pivots
        pivots = count( abs([1:m]-ipiv) .ne. 0 )
        d = d*(-1)**pivots
        !
        ! LAPACK does not contain routines for computing determinants, although determinants can be obtained following
        ! the use of a factorization routine. For example, for the PLU factorization the determinant is the product of
        ! the diagonal elements of U and the sign of P, that is - for an odd number of permutations and + for an even
        ! number. But aware that determinants can readliy overflow, or underflow.
        !
    end function   zdet

    function       ddet(A) result(d)
        ! 
        implicit none
        !
        real(dp), intent(in)  :: A(:,:)
        real(dp) :: d
        real(dp), allocatable :: C(:,:)
        integer , allocatable :: ipiv(:)
        integer :: pivots
        integer :: m, n, info
        integer :: i 
        !
        !
        m = size(A,1)
        n = size(A,2)
        !
        if (m/=n) stop 'Determinant are only defined for square matrices'
        !
        allocate(ipiv(m))
        !
        ! remove numerical singularites by regularization
        allocate(C,source=A)
        do i = 1, m
            C(i,i) = C(i,i) + eps
        enddo
        !
        ! perform LU factorization
        call dgetrf( m, m, C, m, ipiv, info )
        !
        if (info/=0) stop 'Error in zgetrf'
        !
        ! get determinant
        d=1.0_dp
        do i = 1, m
            d = d*C(i,i)
        enddo
        !
        ! correct for number of pivots
        pivots = count( abs([1:m]-ipiv) .ne. 0 )
        d = d*(-1)**pivots
        !
        ! LAPACK does not contain routines for computing determinants, although determinants can be obtained following
        ! the use of a factorization routine. For example, for the PLU factorization the determinant is the product of
        ! the diagonal elements of U and the sign of P, that is - for an odd number of permutations and + for an even
        ! number. But aware that determinants can readliy overflow, or underflow.
        !
    end function   ddet
    
    ! LU decomposition

    subroutine     lu(A)
        ! returns the upper triangular matrix U
        implicit none
        !
        real(dp), intent(inout)  :: A(:,:)
        integer , allocatable :: ipiv(:)
        integer :: m, n, info
        integer :: i, j 
        !
        m = size(A,1)
        n = size(A,2)
        !
        allocate(ipiv(max(1,min(m,n))))
        !
        ! remove numerical singularites by regularization
        do i = 1, min(n,m)
            A(i,i) = A(i,i) + eps
        enddo
        !
        ! perform LU factorization
        call dgetrf( m, n, A, m, ipiv, info )
        !
        if (info/=0) stop 'Error in dgetrf'
        !
        ! return upper triangular matrix U
        do i = 1, m
        do j = 1, min(i-1,n)
            A(i,j) = 0.0_dp
        enddo
        enddo
    end subroutine lu

    ! inverse of real square matrix
    
    function       am_dinv(A) result(Ainv)
        !
        implicit none
        !
        real(dp), intent(in)  :: A(:,:)
        integer,  allocatable :: ipiv(:)
        real(dp), allocatable :: Ainv(:,:)
        real(dp), allocatable :: WORK(:)
        integer  :: m, n, info, lwork
        integer :: i
        !
        m = size(A,1)
        n = size(A,2)
        !
        if (m/=n) stop 'Determinant are only defined for square matrices'
        !
        allocate(ipiv(n))
        allocate(WORK(lwmax))
        !
        ! copy A and remove numerical singularites by regularization
        allocate(Ainv,source=A)
        do i = 1, min(n,m)
            Ainv(i,i) = Ainv(i,i) + eps
        enddo
        !
        ! factorize 
        call dgetrf( m, n, Ainv, m, ipiv, info )
        !
        if (info/=0) stop 'Error in dgetrf'
        !
        ! query workspace
        lwork=-1
        call dgetri( n, Ainv, m, ipiv, WORK, lwork, info )
        lwork = min( lwmax, int( WORK( 1 ) ) )
        !
        ! get inverse
        call dgetri( n, Ainv, m, ipiv, WORK, lwork, info )
        !
        if (info.ne.0) stop 'Error in dgetri'
        !
    end function   am_dinv

    ! inverse of complex square matrixg

    function       am_zinv(A) result(Ainv)
        !
        implicit none
        !
        complex(dp), intent(in)  :: A(:,:)
        integer,  allocatable :: ipiv(:)
        complex(dp), allocatable :: Ainv(:,:)
        complex(dp), allocatable :: WORK(:)
        integer  :: m, n, info, lwork
        integer :: i
        !
        m = size(A,1)
        n = size(A,2)
        !
        if (m/=n) stop 'Determinant are only defined for square matrices'
        !
        allocate(ipiv(n))
        allocate(WORK(lwmax))
        !
        ! copy A and remove numerical singularites by regularization
        allocate(Ainv,source=A)
        do i = 1, min(n,m)
            Ainv(i,i) = Ainv(i,i) + eps
        enddo
        !
        ! factorize 
        call zgetrf( m, n, Ainv, m, ipiv, info )
        !
        if (info/=0) stop 'Error in dgetrf'
        !
        ! query workspace
        lwork=-1
        call zgetri( n, Ainv, m, ipiv, WORK, lwork, info )
        lwork = min( lwmax, int( WORK( 1 ) ) )
        !
        ! get inverse
        call zgetri( n, Ainv, m, ipiv, WORK, lwork, info )
        !
        if (info.ne.0) stop 'Error in dgetri'
        !
    end function   am_zinv

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
        integer :: ldb, nrhs
        integer :: info
        integer :: lwork
        integer :: rank ! effective rank of A: the number of singular values which are greater than rcond*s(1).
        !
        rcond = 1.0D-8
        !
        m = size(A,1)
        n = size(A,2)
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
        call dgelss(m, n, nrhs, C, m, B, ldb, S, rcond, rank, WORK, lwork, info)
        lwork = min(lwmax,int(WORK(1)))
        !
        ! solve the least squares problem min( norm2(B - Ax) ) for the x of minimum norm.
        call dgelss(m, n, nrhs, C, m, B, ldb, S, rcond, rank, WORK, lwork, info)
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
        ! DGEEV computes for an N-by-N real nonsymmetric matrix A, the
        !  eigenvalues and, optionally, the left and/or right eigenvectors.
        !
        !  The right eigenvector v(j) of A satisfies
        !                   A * v(j) = lambda(j) * v(j)
        !  where lambda(j) is its eigenvalue.
        !
        !  The left eigenvector u(j) of A satisfies
        !                u(j)**H * A = lambda(j) * u(j)**H
        !  where u(j)**H denotes the conjugate-transpose of u(j).
        !
        !  The computed eigenvectors are normalized to have Euclidean norm
        !  equal to 1 and largest component real.
        !
        !  If the j-th eigenvalue is real, then uj = vl(:,j), the j-th column of vl. If the j-th and
        !  (j+1)-st eigenvalues form a complex conjugate pair, then for i = sqrt(-1), uj = vl(:,j) +
        !  i*vl(:,j+1) and uj + 1 = vl(:,j)- i*vl(:,j+1).
        !
        !  If the j-th eigenvalue is real, then vj = vr(:,j), the j-th column of vr. If the j-th and 
        !  (j+1)-st eigenvalues form a complex conjugate pair, then for i = sqrt(-1), vj = vr(:,j) + 
        !  i*vr(:,j+1) and vj + 1 = vr(:,j) - i*vr(:,j+1).
        !
        implicit none
        !
        real(dp)   , intent(in) :: A(:,:)
        complex(dp), intent(out), optional :: VR(:,:), VL(:,:) ! right and left eigenvectors
        complex(dp), intent(out), optional :: D(:)
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
            if (size(D).ne.n) stop 'ERROR [am_dgeev]: D dimension mismatch'
            do i = 1, n
                D(i) = cmplx(WR(i),WI(i))
            enddo
        endif
        !
        if (present(VR)) then
            if (size(VR,1).ne.n) stop 'ERROR [am_dgeev]: VR dimension 1 mismatch'
            if (size(VR,2).ne.n) stop 'ERROR [am_dgeev]: VR dimension 2 mismatch'
            i = 0
            do
                i=i+1
                if (abs(WI(i)).gt.tiny) then
                    VR(:,i  ) = cmplx(VR_internal(:,i),+VR_internal(:,i+1),dp)
                    VR(:,i+1) = cmplx(VR_internal(:,i),-VR_internal(:,i+1),dp)
                    i=i+1
                else
                    VR(:,i  ) = VR_internal(:,i)
                endif
                if (i.eq.n) exit
            enddo
        endif
        !
        if (present(VL)) then
            if (size(VL,1).ne.n) stop 'ERROR [am_dgeev]: VL dimension 1 mismatch'
            if (size(VL,2).ne.n) stop 'ERROR [am_dgeev]: VL dimension 2 mismatch'
            i = 0
            do
                i=i+1
                if (abs(WI(i)).gt.tiny) then
                    VL(:,i  ) = cmplx(VL_internal(:,i),+VL_internal(:,i+1),dp)
                    VL(:,i+1) = cmplx(VL_internal(:,i),-VL_internal(:,i+1),dp)
                    i=i+1
                else
                    VL(:,i  ) = VL_internal(:,i)
                endif
                if (i.eq.n) exit
            enddo
        endif
        !
    end subroutine am_dgeev

    ! diagonalize complex general square matrix

    subroutine     am_zgeev(A,VL,D,VR)
        !
        ! ZGEEV computes for an N-by-N complex nonsymmetric matrix A, the
        !  eigenvalues and, optionally, the left and/or right eigenvectors.
        !
        !  The right eigenvector v(j) of A satisfies
        !                   A * v(j) = lambda(j) * v(j)
        !  where lambda(j) is its eigenvalue.
        !
        !  The left eigenvector u(j) of A satisfies
        !                u(j)**H * A = lambda(j) * u(j)**H
        !  where u(j)**H denotes the conjugate transpose of u(j).
        !
        !  The computed eigenvectors are normalized to have Euclidean norm
        !  equal to 1 and largest component real.
        ! 
        implicit none
        !
        complex(dp), intent(in) :: A(:,:)
        complex(dp), intent(out), allocatable, optional :: VR(:,:) ! right and left eigenvectors
        complex(dp), intent(out), allocatable, optional :: VL(:,:) ! right and left eigenvectors
        complex(dp), intent(out), allocatable, optional :: D(:)
        complex(dp), allocatable :: CA(:,:)
        complex(dp), allocatable :: W(:)
        complex(dp), allocatable :: WORK(:)
        real(dp)   , allocatable :: RWORK(:)
        complex(dp), allocatable :: VL_internal(:,:)
        complex(dp), allocatable :: VR_internal(:,:)
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
        allocate(W(n))
        allocate(WORK(lwmax))
        allocate(RWORK(max(lwmax,n)))
        !
        ! copy variables
        !
        allocate(CA, source=A)
        !
        ! query the optimal workspace
        !
        lwork = -1
        call zgeev( 'N', 'V', n, CA, lda, W, VL_internal, ldvl, VR_internal, ldvr, WORK, lwork, RWORK, info )
        lwork = min( lwmax, int( WORK( 1 ) ) )
        !
        ! solve eigenproblem
        !
        if (present(VL).and.present(VR)) then
            call zgeev( 'V', 'V', n, CA, lda, W, VL_internal, ldvl, VR_internal, ldvr, WORK, lwork, RWORK, info )
        elseif (present(VR)) then
            call zgeev( 'N', 'V', n, CA, lda, W, VL_internal, ldvl, VR_internal, ldvr, WORK, lwork, RWORK, info )
        elseif (present(VL)) then
            call zgeev( 'V', 'N', n, CA, lda, W, VL_internal, ldvl, VR_internal, ldvr, WORK, lwork, RWORK, info )
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
            allocate(D(n),source=W)
        endif
        !
        if (present(VR)) then
            allocate(VR,source=VR_internal)
        endif
        !
        if (present(VL)) then
            allocate(VL,source=VL_internal)
        endif
    end subroutine am_zgeev
    
    ! diagonalize complex hermitian matrix

    subroutine     am_zheev(A,V,D)
        !
        implicit none
        !
        complex(dp), intent(in)  :: A(:,:)
        real(dp)   , intent(out), allocatable :: D(:)   ! eigenvalues of the matrix A in ascending order
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

    subroutine     am_zheev_eig(A,D)
        !
        implicit none
        !
        complex(dp), intent(in)  :: A(:,:)
        real(dp)   , allocatable, intent(out) :: D(:)   ! eigenvalues of the matrix A in ascending order
        complex(dp), allocatable :: V(:,:) ! orthonormal eigenvectors of the matrix A
        complex(dp), allocatable :: WORK(:)
        real(dp)   , allocatable :: RWORK(:)
        integer :: n
        integer :: lda
        integer :: info
        integer :: lwork
        !
        !
        ! initialize variables
        lda = size(A,1)
        n   = size(A,2)
        allocate(D(n))
        allocate(RWORK(3*n-2))
        ! copy variables
        allocate(V, source=A)
        ! query the optimal workspace
        allocate(WORK(lwmax))
        lwork = -1
        call zheev( 'N', 'U', n, V, lda, D, WORK, lwork, RWORK, info )
        lwork = min( lwmax, int( WORK( 1 ) ) )
        ! solve eigenproblem
        call zheev( 'N', 'U', n, V, lda, D, WORK, lwork, RWORK, info )
        ! check for convergence.
        if( info.gt.0 ) stop 'ERROR: ZHEEV failed to compute eigenvalues'
        !
    end subroutine am_zheev_eig

    ! diagonalize banded complex hermitian matrix

    subroutine     am_zhbev(A_b,V,D)
        !
        ! noted A is a banded matrix containing either upper or lower triangular part of the Hermitian matrix A.
        ! If A contains the lower triangular part, as it hard coded to do so here, it should be stored like so:
        !
        ! [ a_{1,1} a_{2,2} a_{3,3} ... a_{n-2,n-2} a_{n-1,n-1} a_{n,n}   ] diagonal
        ! [ a_{1,2} a_{2,3} a_{2,4} ... a_{n-1,n  } a_{n-1,n  } 0         ] subdiagonal # 1
        ! [ a_{1,3} a_{2,4} a_{2,5} ... a_{n-1,n+1} 0           0         ] subdiagonal # 2
        ! 
        implicit none
        !
        complex(dp), intent(in)  :: A_b(:,:)
        real(dp)   , intent(out), allocatable :: D(:)   ! eigenvalues of the matrix A in ascending order
        complex(dp), intent(out), allocatable :: V(:,:) ! orthonormal eigenvectors of the matrix A
        complex(dp), allocatable :: A(:,:)
        complex(dp), allocatable :: work(:)
        real(dp)   , allocatable :: rwork(:)
        integer :: kd
        integer :: n
        integer :: lda
        integer :: info
        integer :: ldz
        !
        ! copy Ab so that the original is not overwritten by zhbev
        allocate(A,source=A_b)
        ! initialize variables
        lda = size(A,1) ! leading dimension of ab
        kd  = lda - 1   ! The number of super- or sub-diagonals in A (subtract because of diagonal)
        n   = size(A,2) ! order of A (number of diagonal = number of rows = number of columns = order of matrix)
        ldz = n         ! leading dimension of output array (square matrix on output of order equal to the order of A)
        ! allocate work space
        allocate(work(max(1,n)))
        allocate(rwork(max(1,3*n-2)))
        ! allocate output arrays
        allocate(V(n,n))
        allocate(D(n))
        ! perform checks
        if (lda.le.0) stop 'zhbev : lda .le. 0'
        if (kd .le.0) stop 'zhbev : kd  .le. 0'
        if (n  .le.0) stop 'zhbev : n   .le. 0'
        ! solve eigenproblem
        call zhbev('V', 'L', n, kd, A, lda, D, V, ldz, work, rwork, info)
        ! check for success
        if(info/=0) stop 'zhbev : failed to compute eigenvalues'
        !
    end subroutine am_zhbev

    ! perform svd decomposition on real general rectangular matrix

    subroutine     am_dgesvd(A,U,S,V)
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
        real(dp), intent(in)  :: A(:,:)
        real(dp), intent(out), allocatable, optional :: S(:)
        real(dp), intent(out), allocatable, optional :: U(:,:)
        real(dp), intent(out), allocatable, optional :: V(:,:)
        real(dp), allocatable ::  S_internal(:)
        real(dp), allocatable ::  U_internal(:,:)
        real(dp), allocatable :: VT_internal(:,:)
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
        ! allocate space
        allocate( S_internal(n))
        allocate( U_internal(ldu,m))
        allocate(VT_internal(ldvt,n))
        allocate(WORK(lwmax))
        ! copy variables
        allocate(CA, source=A)
        ! query the optimal workspace
        lwork = -1
        call dgesvd( 'all', 'all', m, n, CA, lda, S_internal, U_internal, ldu, VT_internal, ldvt, WORK, lwork, info )
        lwork = min( lwmax, int( WORK( 1 ) ) )
        ! solve eigenproblem
        call dgesvd( 'all', 'all', m, n, CA, lda, S_internal, U_internal, ldu, VT_internal, ldvt, WORK, lwork, info )
        ! check for convergence
        if( info.gt.0 ) stop 'ERROR [am_dgesvd]: SVD failed to converge'
        ! output
        if (present(S)) allocate(S,source=S_internal)
        if (present(U)) allocate(U,source=U_internal)
        if (present(V)) allocate(V,source=transpose(VT_internal))
    end subroutine am_dgesvd

    subroutine     am_zgesvd(A,U,S,V)
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
        complex(dp), intent(in)  :: A(:,:)
        complex(dp), intent(out), allocatable, optional :: S(:)
        complex(dp), intent(out), allocatable, optional :: U(:,:)
        complex(dp), intent(out), allocatable, optional :: V(:,:)
        complex(dp), allocatable ::  S_internal(:)
        complex(dp), allocatable ::  U_internal(:,:)
        complex(dp), allocatable :: VT_internal(:,:)
        complex(dp), allocatable :: CA(:,:)
        complex(dp), allocatable :: WORK(:)
        real(dp)   , allocatable :: RWORK(:)
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
        ! allocate space
        allocate( S_internal(n))
        allocate( U_internal(ldu,m))
        allocate(VT_internal(ldvt,n))
        allocate(  WORK(lwmax) )
        allocate( RWORK(max(1,5*min(m,n))) )
        ! copy variables
        allocate(CA, source=A)
        ! query the optimal workspace
        lwork = -1
        call zgesvd( 'all', 'all', m, n, CA, lda, S_internal, U_internal, ldu, VT_internal, ldvt, WORK, lwork, RWORK, info )
        lwork = min(lwmax,int(WORK(1)))
        ! solve eigenproblem
        call zgesvd( 'all', 'all', m, n, CA, lda, S_internal, U_internal, ldu, VT_internal, ldvt, WORK, lwork, RWORK, info )
        ! check for convergence
        if( info.gt.0 ) stop 'ERROR [am_dgesvd]: SVD failed to converge'
        ! output
        if (present(S)) allocate(S,source=S_internal)
        if (present(U)) allocate(U,source=U_internal)
        if (present(V)) allocate(V,source=conjg(transpose((VT_internal))))
    end subroutine am_zgesvd

    ! get null space using svd decomposition

    function       null_space_dgesvd(A) result(null_space)
        !
        real(dp), intent(in)  :: A(:,:)
        real(dp), allocatable :: null_space(:,:)
        real(dp), allocatable :: S_internal(:)
        real(dp), allocatable :: V_internal(:,:)
        logical , allocatable :: T(:)
        integer , allocatable :: inds(:)
        ! perform svd decomposition
        call am_dgesvd(A=A,S=S_internal,V=V_internal)
        ! allocate mask
        allocate(T( size(S_internal) ))
        ! create mask
        where (abs(S_internal).lt.tiny)
            T = .true.
        elsewhere
            T = .false.
        endwhere
        ! create indices
        allocate(inds, source = [1:size(V_internal,2)])
        ! return null_space
        allocate(null_space, source = V_internal(:, pack(inds,T) ))
    end function   null_space_dgesvd

    function       null_space_zgesvd(A) result(null_space)
        !
        complex(dp), intent(in)  :: A(:,:)
        complex(dp), allocatable :: null_space(:,:)
        complex(dp), allocatable :: S_internal(:)
        complex(dp), allocatable :: V_internal(:,:)
        logical , allocatable :: T(:)
        integer , allocatable :: inds(:)
        ! perform svd decomposition
        call am_zgesvd(A=A,S=S_internal,V=V_internal)
        ! allocate mask
        allocate(T( size(S_internal) ))
        ! create mask
        where (abs(S_internal).lt.tiny)
            T = .true.
        elsewhere
            T = .false.
        endwhere
        ! create indices
        allocate(inds, source = [1:size(V_internal,2)])
        ! return null_space
        allocate(null_space, source = V_internal(:, pack(inds,T) ))
    end function   null_space_zgesvd

    ! solve Ax = B by QR factorization

    function       am_dgels(A,B) result(X)
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
        real(dp), allocatable :: X(:,:)
        real(dp), allocatable :: CA(:,:)
        real(dp), allocatable :: CB(:,:)
        real(dp), allocatable :: WORK(:)
        integer :: m, n, nrhs
        integer :: lda, ldb
        integer :: info
        integer :: lwork
        !
        ! initialize variables
        !
        m   = size(A,1)
        n   = size(A,2)
        lda = m          ! at least max(1, m).
        ldb = max(1,m,n) ! must be at least max(1, m, n).
        nrhs= size(B,2)  ! The second dimension of b must be at least max(1, nrhs).
        if (lda.ne.m) stop 'ERROR [am_dgels]: lda /= ldb'
        ! copy variables
        allocate(CA, source = A)
        allocate(CB(ldb,nrhs))
        CB = 0
        CB(1:size(B,1),1:nrhs) = B
        ! query the optimal workspace
        allocate(WORK(1))
        lwork = -1
        call dgels( 'No transpose', m, n, nrhs, CA, lda, CB, ldb, WORK, lwork, info )
        lwork = max( min(m,n)+max(1,m,n,nrhs), min( lwmax, int(WORK(1))) )
        ! reallocate work space
        deallocate(WORK)
        allocate(WORK(lwmax))
        ! call solver
        call dgels( 'No transpose', m, n, nrhs, CA, lda, CB, ldb, WORK, lwork, info )
        ! check for convergence.
        if( info.gt.0 ) then
            write(*,*) 'The diagonal element ',info,' of the triangular '
            write(*,*) 'factor of a is zero, so that a does not have full '
            write(*,*) 'rank; the least squares solution could not be '
            write(*,*) 'computed.'
            stop
        end if
        ! return X
        allocate(X,source=CB)
    end function   am_dgels

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

    ! LU decomposition

    subroutine     am_dgetrf(A,L,U)
        ! NOT YET TESTED.
        implicit none
        !
        real(dp), intent(in)  :: A(:,:)
        real(dp), allocatable :: L(:,:)
        real(dp), allocatable :: U(:,:)
        integer , allocatable :: ipiv(:)
        integer :: m, n, info
        integer :: i, j 
        !
        m = size(A,1)
        n = size(A,2)
        allocate(ipiv(max(1,min(m,n))))
        !
        ! copy A and remove numerical singularites by regularization
        allocate(U,source=A)
        do i = 1, min(n,m)
            U(i,i) = U(i,i) + eps
        enddo
        !
        ! perform LU factorization
        call dgetrf( m, n, U, m, ipiv, info )
        !
        if (info/=0) stop 'Error in dgetrf'
        !
        ! return lower triangular matrix L
        allocate(L,source=U)
        do j = 1, n
        do i = 1, min(i,m)
            if (i.eq.j) then
                U(i,j) = 1.0_dp
            else
                U(i,j) = 0.0_dp
            endif
        enddo
        enddo
        !
        ! return upper triangular matrix U
        do i = 1, m
        do j = 1, min(i-1,n)
            U(i,j) = 0.0_dp
        enddo
        enddo
    end subroutine am_dgetrf

    ! multiply A and B matrices

    function       am_dgemm(A,B,flags) result(C)
        !
        ! multiply matrices A and B
        ! 
        implicit none
        !
        real(dp), intent(in) :: A(:,:)
        real(dp), intent(in) :: B(:,:)
        character(*), intent(in), optional :: flags
        real(dp), allocatable :: C(:,:)
        character(len=1) :: transa
        character(len=1) :: transb
        integer :: Am,An,Bm,Bn
        integer :: m,n,k
        !
        Am = size(A,1)
        An = size(A,2)
        Bm = size(B,1)
        Bn = size(B,2)
        !
        transa='N'
        transb='N'
        if (present(flags)) then
            if (index(flags,'AT').ne.0) then
                transa = 'T'
            endif
            if (index(flags,'AC').ne.0) then
                transa = 'C'
            endif
            if (index(flags,'BT').ne.0) then
                transb = 'T'
            endif
            if (index(flags,'BC').ne.0) then
                transb = 'C'
            endif
        endif
        !
        ! m Specifies the number of rows of the matrix op(A)
        m = Am
        if ((transa(1:1).eq.'T').or.(transa(1:1).eq.'C')) m = An
        !
        ! n Specifies the number of columns of the matrix op(B)
        n = Bn
        if ((transb(1:1).eq.'T').or.(transb(1:1).eq.'C')) n = Bm
        !
        ! k Specifies the number of columns of the matrix op(A)
        k = An
        if ((transa(1:1).eq.'T').or.(transa(1:1).eq.'C')) k = Am
        !
        allocate(C(m,n))
        C = 0
        !
        ! C := alpha*op(A)*op(B) + beta*C,
        !            [m*k] [k*n]
        !
        !    dgemm(transa, transb, m, n, k, 1.0_dp, A, lda, B, ldb, 0.0_dp , C, ldc)
        call dgemm(transa, transb, m, n, k, 1.0_dp, A,  Am, B,  Bm, 0.0_dp , C,   m)
        !
    end function   am_dgemm

    ! multiply A and B and get a symmetric matrix C (only supported in ifort v 16)

!     function       am_dgemmt(A,B,flags) result(C)
!         !
!         implicit none
!         !
!         real(dp), intent(in) :: A(:,:)
!         real(dp), intent(in) :: B(:,:)
!         character(*), intent(in), optional :: flags
!         real(dp), allocatable :: C(:,:)
!         character(len=1) :: uplo
!         character(len=1) :: transa
!         character(len=1) :: transb
!         integer :: Am,An,Bm,Bn
!         integer :: m,n,k
!         integer :: i, j
!         !
!         Am = size(A,1)
!         An = size(A,2)
!         Bm = size(B,1)
!         Bn = size(B,2)
!         !
!         uplo  ='U'
!         transa='N'
!         transb='N'
!         if (present(flags)) then
!             if (index(flags,'U').ne.0) then
!                 uplo    = 'T'
!             endif
!             if (index(flags,'L').ne.0) then
!                 uplo    = 'L'
!             endif
!             if (index(flags,'AT').ne.0) then
!                 transa = 'T'
!             endif
!             if (index(flags,'AC').ne.0) then
!                 transa = 'C'
!             endif
!             if (index(flags,'BT').ne.0) then
!                 transb = 'T'
!             endif
!             if (index(flags,'BC').ne.0) then
!                 transb = 'C'
!             endif
!         endif
!         !
!         ! m Specifies the number of rows of the matrix op(A)
!         m = Am
!         if ((transa(1:1).eq.'T').or.(transa(1:1).eq.'C')) m = An
!         !
!         ! n Specifies the number of columns of the matrix op(B)
!         n = Bn
!         if ((transb(1:1).eq.'T').or.(transb(1:1).eq.'C')) n = Bm
!         !
!         if (n.ne.m) then
!             write(*,*) 'ERROR: The final matrix must be square. n /= m.'
!             stop
!         endif
!         !
!         ! k Specifies the number of columns of the matrix op(A)
!         k = An
!         if ((transa(1:1).eq.'T').or.(transa(1:1).eq.'C')) k = Am
!         !
!         allocate(C(m,n))
!         C = 0
!         !
!         ! C := alpha*op(A)*op(B) + beta*C,
!         !            [m*k] [k*n]
!         !
!         !    dgemmt(uplo, transa, transb, n, k,  alpha, a, lda, b, ldb,   beta, c, ldc)
!         call dgemmt(uplo, transa, transb, n, k, 1.0_dp, A,  Am, B,  Bm, 0.0_dp, C,   m)
!         !
!         do i = 1, m
!             do j = 1, (i-1)
!                 C(i,j) = C(j,i)
!             enddo
!         enddo
!         !
!     end function   am_dgemmt

end module am_mkl
 















