module am_matlab

    implicit none
    
    use am_constants

    interface linspace
        module procedure linspace_double, linspace_integer
    end interface ! linspace

    interface diag
        module procedure diag1, diag2
    end interface ! diag
    
contains

    !
    ! matlab-inspired functions
    !

    pure function Ylm(l,m,theta,phi)
        !
        ! computes spherical harmonics. Theta and phi in radians. Theta and phi should be the same size.
        !
        ! W. H. Press, B. P. Flannery, S. A. Teukolsky, and W. T. Vetterling, Numerical Recipes in Fortran 77: The Art
        ! of Scientific Computing, 2 edition (Cambridge University Press, Cambridge England ; New York, 1992), p 246.
        ! 
        implicit none
        !
        integer , intent(in)  :: l, m
        real(dp), intent(in) :: theta(:), phi(:)
        complex(dp) :: Ylm(size(theta))
        !
        Ylm = sqrt( real(2*l+1,dp)/fourpi * factorial(l-m)/real(factorial(l+m),dp) ) * legendre(l,m,cos(theta)) * exp(cmplx_i*m*phi)
        !
    end function  Ylm

    pure function factorial(n) result(y)
        !
        implicit none
        !
        integer, intent(in) :: n
        integer :: y
        !
        if (n.ge.1) then
            y = product([1:n])
        else
            y = 0
        endif
    end function  factorial

    pure function legendre(l,m,x) result(y)
        !
        ! computes the associated legrende polynomial, x = cos(theta)
        !
        ! W. H. Press, B. P. Flannery, S. A. Teukolsky, and W. T. Vetterling, Numerical Recipes in Fortran 77: The Art
        ! of Scientific Computing, 2 edition (Cambridge University Press, Cambridge England ; New York, 1992), p 246.
        !
        implicit none
        !
        integer , intent(in) :: l,m
        real(dp), intent(in) :: x(:)
        real(dp) :: y(size(x))
        integer( :: ll
        real(dp) :: pll(size(x))
        real(dp) :: pmm(size(x))
        real(dp) :: pmmp1(size(x))
        real(dp) :: somx2(size(x))
        !
        ! call assert(m >= 0, m <= l, all(abs(x) <= 1.0), 'y args')
        !
        pmm=1.0
        if (m > 0) then
            somx2=sqrt((1.0_dp-x)*(1.0_dp+x))
            pmm=product(arth(1.0_dp,2.0_dp,m))*somx2**m
            if (mod(m,2) == 1) pmm=-pmm
        end if
        if (l == m) then
            y=pmm
        else
            pmmp1=x*(2*m+1)*pmm
            if (l == m+1) then
                y=pmmp1
            else
                do ll=m+2,l
                    pll=(x*(2*ll-1)*pmmp1-(ll+m-1)*pmm)/(ll-m)
                    pmm=pmmp1
                    pmmp1=pll
                end do
                y=pll
            end if
        end if
        contains
        pure function arth(first,increment,n)
            !
            implicit none
            !
            real(dp), intent(in) :: first
            real(dp), intent(in) :: increment
            integer , intent(in) :: n
            real(dp) :: arth_d(n)
            integer  :: k,k2
            real(dp) :: temp
            !
            if (n > 0) arth_d(1)=first
            if (n <= npar_arth) then
                do k=2,n
                    arth_d(k)=arth_d(k-1)+increment
                end do
            else
                do k=2,npar2_arth
                    arth_d(k)=arth_d(k-1)+increment
                end do
                temp=increment*npar2_arth
                k=npar2_arth
                do
                    if (k >= n) exit
                    k2=k+k
                    arth_d(k+1:min(k2,n))=temp+arth_d(1:min(k,n-k))
                    temp=temp+temp
                    k=k2
                end do
            end if
        end function arth
    end function  legendre

    function      perms(n) result(PT)
        !
        implicit none
        !
        integer, intent(in) :: n
        integer, allocatable :: P(:,:)
        integer, allocatable :: PT(:,:)
        integer :: m
        !
        m = product([1:n])
        !
        allocate(P(m,n))
        !
        call permutate([1:n],P)
        !
        allocate(PT,source=transpose(P))
        !
        contains
        recursive subroutine permutate(E, P)
            !
            implicit none
            !
            integer, intent(in)  :: E(:) ! array of objects 
            integer, intent(out) :: P(:,:) ! permutations of E 
            integer :: N, Nfac, i, k, S(size(P,1)/size(E), size(E)-1) 
            N = size(E); Nfac = size(P,1); 
            do i = 1, N
              if( N>1 ) call permutate((/E(:i-1), E(i+1:)/), S) 
              forall(k=1:Nfac/N) P((i-1)*Nfac/N+k,:) = (/E(i), S(k,:)/) 
            enddo 
        end subroutine permutate 
    end function  perms

    function      fopen(filename,permission) result(fid)
        !
        implicit none
        !
        character(*),intent(in) :: filename
        character(*),intent(in) :: permission
        integer :: fid
        integer :: status
        !
        status=0
        !
        fid=unit_number(filename)
        !
        if (index(permission,'r')) then
            ! read
            open(unit=fid, file=filename, status='old'    , action='read',iostat=status)
        elseif (index(permission,'w')) then
            ! write
            open(unit=fid, file=filename, status='replace', action='write',iostat=status)
        elseif (index(permission,'a')) then 
            ! append
            open(unit=fid, file=filename, status='old'    , action='write',iostat=status)
        endif
        if (status.ne.0) then
            write(*,*) 'Could not open ', filename
            stop
        endif
        contains
        integer function unit_number(filename)
            !
            implicit none
            !
            character(len=*),intent(in) :: filename
            !
            integer :: unit, i=100
            logical :: isopen
            !
            inquire(file=filename,number=unit,opened=isopen)
            if (isopen.eqv..true.) then
                unit_number=unit
            else
                do
                    inquire(unit=i,opened=isopen)
                    if ( isopen .eqv. .false. ) then
                        unit_number=i
                        exit
                    endif
                    i=i+1
                    if ( i == 1000 ) then
                        write(*,*) 'No available units between 100 and 1000. Are you sure you want that'
                        write(*,*) 'many open files?'
                    endif
                enddo
            endif
        end function   unit_number
    end function  fopen

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

    pure function kron(A,B) result(C)
        ! kronecker product
        implicit none
        !
        real(dp), intent(in) :: A(:,:)
        real(dp), intent(in) :: B(:,:)
        real(dp), allocatable :: C(:,:)
        integer :: Am, An, Bm, Bn, i, j
        !
        Am = size(A,1)
        An = size(A,2)
        Bm = size(B,1)
        Bn = size(B,2)
        !
        allocate(C(Am*Bm,An*Bn))
        !
        do i = 1,Am
        do j = 1,An
        C( [1:Bm]+Bm*(i-1), [1:Bn]+Bn*(j-1) ) = A(i,j)*B
        enddo
        enddo
    end function  kron

    pure function kron_pow(A,n) result(C)
        ! kronecker power
        implicit none
        !
        real(dp), intent(in) :: A(:,:)
        integer , intent(in) :: n
        real(dp), allocatable :: Cin(:,:)
        real(dp), allocatable :: C(:,:)
        integer :: i
        !
        allocate(C,source=A)
        !
        ! if n = 1, C = A^1 
        ! if n = 2, C = A^2 
        ! etc... 
        do i = 1, (n-1)
            Cin = kron(C,A)
            deallocate(C)
            allocate(C,source=Cin)
        enddo
        !
    end function  kron_pow

    pure function eye(n)
        !> nxn identity matrix
        implicit none
        !
        integer, intent(in) :: n
        real(dp), dimension(:,:), allocatable :: eye
        integer :: i
        !
        allocate(eye(n,n))
        eye=0.0_dp
        do i = 1, n
            eye(i,i) = 1.0_dp
        enddo
    end function  eye

    pure function ones(n)
        !> nxn identity matrix
        implicit none
        !
        integer, intent(in) :: n
        real(dp), dimension(:,:), allocatable :: ones
        !
        allocate(ones(n,n))
        ones=1.0_dp
        !
    end function  ones

    pure function eye_nxm(n,m) result(id)
        !> nxm identity matrix
        implicit none
        !
        integer, intent(in) :: n, m
        real(dp), dimension(:,:), allocatable :: id
        integer :: i
        !
        allocate(id(n,m))
        id=0.0_dp
        do i = 1, min(n,m)
            id(i,i) = 1.0_dp
        enddo
    end function  eye_nxm

    pure subroutine rref(matrix)
        !
        ! note: algorithm on roseta code is broken.
        ! this has been fixed by:
        ! 1) changed: "(n.lt.pivot)" and "(pivot.gt.n)"
        ! 2) changed: "(i.gt.m)"
        !
        implicit none
        !
        real(dp), intent(inout) :: matrix(:,:)
        real(dp), allocatable :: temp(:)
        integer :: pivot, m, n
        integer :: r, i
        !
        pivot = 1
        m=size(matrix,1)
        n=size(matrix,2)
        !
        allocate(temp(n))
        !
        do r = 1, m
        if (n.lt.pivot) exit
        i = r
        do while (abs(matrix(i,pivot)).lt.tiny)
            i=i+1
            if (i.gt.m) then
                i=r
                pivot=pivot+1
                if (pivot.gt.n) return
            end if
        end do
        ! swap rows i and r
        temp = matrix(i,:)
        matrix(i,:) = matrix(r,:)
        matrix(r,:) = temp
        !
        matrix(r,:) = matrix(r,:)/matrix(r,pivot)
        do i = 1, m
            if (i.ne.r) matrix(i,:)=matrix(i,:)-matrix(r,:)*matrix(i,pivot)
        end do
        pivot = pivot + 1
        end do
        deallocate(temp)
    end subroutine  rref

end module am_matlab
