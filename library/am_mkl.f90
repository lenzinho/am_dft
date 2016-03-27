module am_mkl
    !
    use am_constants
    use am_helpers
    !
    public

contains

    
    function inv(A) result(Ainv)
        !
        ! 1) factorize: dgetrf( m, n, a, lda, ipiv, info )
        ! The routine computes the LU factorization of a general m-by-n matrix A as A = P*L*U, where P is a permutation
        ! matrix, L is lower triangular with unit diagonal elements (lower trapezoidal if m > n) and U is upper
        ! triangular (upper trapezoidal if m < n). The routine uses partial pivoting, with row interchanges.
        !
        !  m    INTEGER. The number of rows in the matrix A (m ≥ 0).
        !  n    INTEGER. The number of columns in A; n ≥ 0.
        !  a    REAL for sgetrf
        !       DOUBLE PRECISION for dgetrf
        !       COMPLEX for cgetrf
        !       DOUBLE COMPLEX for zgetrf.
        !       Array, DIMENSION (lda,*). Contains the matrix A. The second dimension of a must be at least max(1, n).
        !  lda  INTEGER. Leading dimension of array A.
        !
        ! 2) invert
        !
        !
        !
        !
        !
        !
        ! NOTE: Input arguments such as array dimensions are not required in Fortran 95 and are skipped from the calling
        ! sequence. Array dimensions are reconstructed from the user data that must exactly follow the required array
        ! shape.
        !
        ! NOTE: Argument info is declared as optional in the Fortran 95 interface. If it is present in the calling
        ! sequence, the value assigned to info is interpreted as follows:
        !     - If this value is more than -1000, its meaning is the same as in the FORTRAN 77 routine.
        !     - If this value is equal to -1000, it means that there is not enough work memory.
        !     - If this value is equal to -1001, incompatible arguments are present in the calling sequence.
        ! 
        ! NOTE: Another type of generic arguments that are skipped in the Fortran 95 interface are arguments that
        ! represent workspace arrays (such as work, rwork, and so on). The only exception are cases when workspace
        ! arrays return significant information on output.
        !
        implicit none
        !
        real(dp), intent(in)  :: A(:,:)
        integer,  allocatable :: ipiv(:)
        real(dp), allocatable :: Ainv(:,:)
        real(dp), allocatable :: work(:)
        integer  :: m, n, lda, info, lwork
        !
        m = size(A,1)
        n = size(A,2)
        allocate(ipiv(n))
        allocate(Ainv(n,n))
        Ainv = A
        lda  = m
        !
        if (m.ne.n) then
            call am_print('ERROR','Matrix is not square.',' >>> ')
            call am_print('shape(A)',shape(A))
            stop
        endif
        !
        ! 1) factorize 
        !
        ! call dgetrf( m, n, a, lda, ipiv, info )
        call dgetrf( m, n, ainv, lda, ipiv, info )
        if (info.ne.0) then
            call am_print('ERROR','Internal problem with factorization.')
            call am_print('A',A)
            call am_print('Ainv',Ainv)
            call am_print('info',info)
            stop
        endif
        !
        ! 2) get the inverse
        !
        ! call dgetri( n, a, lda, ipiv, work, lwork, info )
        !
        ! 2a) query optimal workspace
        !
        ! lwork INTEGER. The size of the work array; lwork ≥ n.
        !       If lwork = -1, then a workspace query is assumed; the routine only calculates the optimal size of the work
        !       array, returns this value as the first entry of the work array, and no error message related to lwork is
        !       issued by xerbla.
        allocate(work(1))
        lwork=-1
        call dgetri( n, ainv, lda, ipiv, work, lwork, info )
        !
        ! 2b) perform inversion
        !
        lwork = nint(work(1))
        deallocate(work)
        allocate(work(lwork))
        call dgetri( n, ainv, lda, ipiv, work, lwork, info )
        !
        if (info.ne.0) then
            call am_print('ERROR','Internal problem with inversion.')
            call am_print('A',A)
            call am_print('Ainv',Ainv)
            call am_print('info',info)
            stop
        endif
        !
    end function inv




end module