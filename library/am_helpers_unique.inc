
    pure function unique_matrices_double(A,iopt_tiny) result(B)
        !> returns unique matrices of A(:,:,i) within numerical precision
        implicit none
        !
        real(dp), intent(in) :: A(:,:,:)
        real(dp) :: wrkspace(size(A,1),size(A,2),size(A,3))
        real(dp), allocatable :: B(:,:,:)
        integer :: i,j,k
        real(dp), intent(in), optional :: iopt_tiny
        real(dp) :: tiny_prec
        !
        if ( present(iopt_tiny) ) then
            tiny_prec = iopt_tiny
        else
            tiny_prec = tiny
        endif
        !
        k=1
        wrkspace(:,:,1) = A(:,:,1)
        !
        try_loop : do i = 1, size(A,3)
            do j = 1, k
                if ( all(abs(A(:,:,i)-wrkspace(:,:,j)).lt.tiny_prec) ) cycle try_loop
            enddo
            k = k + 1
            wrkspace(:,:,k) = A(:,:,i)
        enddo try_loop
        !
        allocate(B(size(A,1),size(A,2),k))
        B(:,:,1:k) = wrkspace(:,:,1:k)
        !
    end function  unique_matrices_double

    pure function unique_columns_double(A,iopt_tiny) result(B)
        !> returns unique columns of double matrix A(:,i) within numerical precision
        implicit none
        !
        real(dp), intent(in) :: A(:,:)
        real(dp) :: wrkspace(size(A,1),size(A,2))
        real(dp), allocatable :: B(:,:)
        integer :: i,j,k
        real(dp), intent(in), optional :: iopt_tiny
        real(dp) :: tiny_prec
        !
        if ( present(iopt_tiny) ) then
            tiny_prec = iopt_tiny
        else
            tiny_prec = tiny
        endif
        !
        k=1
        wrkspace(:,1) = A(:,1)
        !
        try_loop : do i = 1,size(A,2)
            do j = 1,k
                if ( all(abs(A(:,i)-wrkspace(:,j)).lt.tiny_prec) ) cycle try_loop
            enddo
            k = k + 1
            wrkspace(:,k) = A(:,i)
        enddo try_loop
        !
        allocate(B(size(A,1),k))
        !
        do i = 1,k
            B(:,i) = wrkspace(:,i)
        enddo
    end function  unique_columns_double

    pure function unique_columns_integer(A) result(B)
        !> returns unique columns of double matrix A(:,i) within numerical precision
        implicit none
        !
        integer, intent(in) :: A(:,:)
        integer :: wrkspace(size(A,1),size(A,2))
        integer, allocatable :: B(:,:)
        integer :: i,j,k
        !
        k=1
        wrkspace(:,1) = A(:,1)
        !
        try_loop : do i = 1,size(A,2)
            do j = 1,k
                if ( all(abs(A(:,i)-wrkspace(:,j)).eq.0) ) cycle try_loop
            enddo
            k = k + 1
            wrkspace(:,k) = A(:,i)
        enddo try_loop
        !
        allocate(B(size(A,1),k))
        !
        do i = 1,k
            B(:,i) = wrkspace(:,i)
        enddo
    end function  unique_columns_integer
    
    pure function unique_integer(A) result(B)
        !> returns unique values of integer array A(:)
        implicit none
        !
        integer, intent(in) :: A(:)
        integer :: wrkspace(size(A))
        integer, allocatable :: B(:)
        integer :: i, k
        !
        k=1
        wrkspace(1) = A(1)
        !
        try_loop : do i = 1, size(A)
            if (.not.any(A(i).eq.wrkspace(:))) then
                k = k + 1
                wrkspace(k) = A(i)
            endif
        enddo try_loop
        !
        allocate(B(k))
        !
        do i = 1, k
            B(i) = wrkspace(i)
        enddo
    end function  unique_integer