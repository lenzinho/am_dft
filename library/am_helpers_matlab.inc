    !
    ! matlab-inspired functions
    !
    
    function      exists(fname)
        !
        implicit none
        !
        character(max_argument_length) :: fname
        logical :: exists
        !
        inquire(file=trim(fname),exist=exists)
        !
    end function  exists

    pure function strsplit(str,delimiter) result(word)
        !
        character(maximum_buffer_size), intent(in) :: str
        character(1), intent(in) :: delimiter
        character(len=:), allocatable :: word(:) ! character(500) :: word(500)
        integer :: i, j, k
        !
        j = 0; k = 0
        do i = 1, len(str)
            if (i .gt. k) then
            if ( str(i:i) .ne. delimiter ) then
                ! determine the next position of delimiter
                k = i+index(str(i:),delimiter)-1
                ! increase word count
                j=j+1
            endif
            endif
        enddo
        !
        allocate(character(maximum_buffer_size) :: word(j))
        !
        j = 0; k = 0
        do i = 1, len(str)
            if (i .gt. k) then
            if ( str(i:i) .ne. delimiter ) then
                ! determine the next position of delimiter
                k = i+index(str(i:),delimiter)-1
                ! increase word count
                j=j+1
                word(j) = str(i:k)
                ! write(*,"(I,A)") j, word(j)
            endif
            endif
        enddo
        !
    end function  strsplit

    pure function trace(R) result(tr)
        !
        implicit none
        !
        real(dp), intent(in) :: R(:,:)
        real(dp) :: tr
        integer :: i
        !
        tr = 0
        do i = 1,size(R,2)
            tr = tr + R(i,i)
        enddo
        !
    end function  trace
    
    pure function cross_product(a,b) result(c)
        !
        implicit none
        !
        real(dp) :: c(3)
        real(dp), INTENT(IN) :: a(3), b(3)
        !
        c(1) = a(2) * b(3) - a(3) * b(2)
        c(2) = a(3) * b(1) - a(1) * b(3)
        c(3) = a(1) * b(2) - a(2) * b(1)
        !
    end function  cross_product

    pure function lcm(a,b)
        !
        implicit none
        !
        integer, intent(in) :: a,b
        integer :: lcm
        !
        lcm = a*b / gcd(a,b)
        !
    end function  lcm
 
    pure function gcd(a,b)
        !
        implicit none
        !
        integer, intent(in) :: a, b
        integer :: t, bc, ac
        integer :: gcd
        !
        bc = b
        ac = a
        !
        do while (bc/=0)
            t = bc
            bc = mod(ac,bc)
            ac = t
        end do
        gcd = abs(ac)
        !
    end function  gcd

    pure function linspace_double(d1,d2,n) result(y)
        !
        implicit none
        !
        real(dp), intent(in) :: d1
        real(dp), intent(in) :: d2
        integer, intent(in)  :: n
        real(dp) :: d
        real(dp) :: y(n)
        integer :: i
        !
        d = (d2-d1)/(n-1.0_dp)
        do i = 1,n
            y(i) = real(d1,dp) + (i-1.0_dp)*d;
        enddo
        !
    end function  linspace_double

    pure function linspace_integer(d1,d2,n) result(y)
        !
        implicit none
        !
        integer, intent(in) :: d1
        integer, intent(in) :: d2
        integer, intent(in)  :: n
        real(dp) :: d
        real(dp) :: y(n)
        integer :: i
        !
        d = (d2-d1)/(n-1.0_dp)
        do i = 1,n
            y(i) = real(d1,dp) + (i-1.0_dp)*d;
        enddo
        !
    end function  linspace_integer

    pure function primes(nprimes)
	    !
    	! naive approach to generating the first n prime numbers
    	!
		implicit none
		!
		integer, intent(in)  :: nprimes
		integer, allocatable :: primes(:) ! array that will hold the primes
		integer :: at, found, i
		logical :: is_prime
		!
		allocate (primes(nprimes))
		!
		primes(1) = 2
		at = 2
		found = 1
		do
		    is_prime = .true. ! assume prime
		    do i = 1, found
		        if (modulo(at,primes(i)).eq.0) then ! if divisible by any other element
		            is_prime = .false.               ! in the array, then not prime.
		            at = at + 1
		            continue
		        end if
		    end do
		    found = found + 1
		    primes(found) = at
		    at = at + 1
		    if (found == nprimes) then ! stop when all primes are found
		        exit
		    endif
		end do
		!
	end function  primes

    pure function diag1(M) result(d)
        !
        ! get diagonal elements of matrix M
        implicit none
        !
        real(dp), intent(in)  :: M(:,:)
        real(dp), allocatable :: d(:)
        integer :: n
        integer :: i
        !
        n = min(size(M,1),size(M,2))
        !
        allocate(d(n))
        !
        do i = 1, n
            d(i) = M(i,i)
        enddo
        !
    end function  diag1

    pure function diag2(d) result(M)
        !
        ! get diagonal elements of matrix M
        implicit none
        !
        real(dp), intent(in)  :: d(:)
        real(dp), allocatable :: M(:,:)
        integer :: n
        integer :: i
        !
        n = size(d,1)
        !
        allocate(M(n,n))
        M=0
        !
        do i = 1, n
            M(i,i) = d(i)
        enddo
        !
    end function  diag2






