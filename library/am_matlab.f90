module am_matlab

    use am_constants

    implicit none

    interface linspace
        module procedure linspace_double, linspace_integer
    end interface ! linspace

    interface diag
        module procedure diag1, diag2
    end interface ! diag

    interface eye
        module procedure eye, eye_nxm
    end interface ! eye

    interface ones
        module procedure ones, ones_nxm
    end interface ! ones

    interface unique
        module procedure unique_columns_double, unique_matrices_double, unique_columns_integer, unique_integer
    end interface ! unique

    interface meshgrid
        module procedure dmeshgrid, imeshgrid
    end interface ! meshgrid
    
contains

    !
    ! matlab-inspired functions
    !

    pure function rotmat(a,b) result(R)
        !
        ! Rotation matrix R which aligns unit vector (Ra) to unit vector b (that is, B = R*A)
        !
        ! based on Rodrigues' Rotation Formula
        ! "A Mathematical Introduction to Robotic Manipulation", Richard M. Murray, Zexiang Li, S. Shankar Sastry, pp. 26-28
        !
        ! % quick matlab implementation:
        ! a = rand(3,1); a=a./norm(a);
        ! b = rand(3,1); b=b./norm(a);
        ! v = cross(a,b);
        ! s = norm(v);
        ! c = dot(a,b);
        ! vx(1:3,1) = [0.0,v(3),-v(2)];
        ! vx(1:3,2) = [-v(3),0.0,v(1)];
        ! vx(1:3,3) = [v(2),-v(1),0.0];
        ! R = eye(3) + vx + (vx*vx)*(1.0-c)/(s.^2)
        ! R*a-b
        !
        implicit none
        !
        real(dp), intent(in) :: a(3)
        real(dp), intent(in) :: b(3)
        real(dp) :: R(3,3)
        real(dp) :: s ! sine of angle
        real(dp) :: c ! cosine of angle
        real(dp) :: v(3) ! cross product of a and b
        real(dp) :: vx(3,3) ! skew symmetric cross product of v
        ! 
        v = cross_product(a,b)
        s = norm2(v)
        c = dot_product(a,b)
        vx(1:3,1) = [0.0_dp,v(3),-v(2)]
        vx(1:3,2) = [-v(3),0.0_dp,v(1)]
        vx(1:3,3) = [v(2),-v(1),0.0_dp]
        !
        R = eye(3) + vx + matmul(vx,vx)*(1.0_dp - c)/(s**2)
        !
    end function  rotmat

    pure function dmeshgrid(n1,n2,n3) result(grid_points)
        !> 
        !> Generates a mesh of lattice vectors (fractional coordinates) around the origin. 
        !> n(1), n(2), n(3) specifies the number of mesh points away from 0, i.e. [-n:1:n]
        !> To obtain a mesh of reciprocal lattice vectors: matmul(recbas,grid_points)
        !> To obtain a mesh of real-space lattice vectors: matmul(bas,grid_points)
        !>
        implicit none
        !
        real(dp), intent(in) :: n1(:)
        real(dp), intent(in) :: n2(:)
        real(dp), intent(in) :: n3(:)
        real(dp), allocatable :: grid_points(:,:) !> voronoi points (27=3^3)
        integer :: i1,i2,i3,j
        !
        allocate(grid_points(3,size(n1)*size(n2)*size(n3)))
        !
        j=0
        do i1 = 1, size(n1)
        do i2 = 1, size(n2)
        do i3 = 1, size(n3)
            j=j+1
            grid_points(1,j)=n1(i1)
            grid_points(2,j)=n2(i2)
            grid_points(3,j)=n3(i3)
        enddo
        enddo
        enddo
    end function  dmeshgrid

    pure function imeshgrid(n1,n2,n3) result(grid_points)
        !
        implicit none
        !
        integer, intent(in) :: n1(:)
        integer, intent(in) :: n2(:)
        integer, intent(in) :: n3(:)
        real(dp), allocatable :: grid_points(:,:)
        !
        grid_points=dmeshgrid(real(n1,dp),real(n2,dp),real(n3,dp))
        !
    end function  imeshgrid

    pure function factorial(n) result(y)
        !
        implicit none
        !
        integer, intent(in) :: n
        real(dp) :: y
        !
        if (n.ge.1) then
            y = product([1:n])
        else
            y = 1
        endif
    end function  factorial

    pure function legendre(l,m,x) result(y)
        !
        ! Associated legrende polynomial
        !
        ! W. H. Press, B. P. Flannery, S. A. Teukolsky, and W. T. Vetterling, Numerical Recipes in Fortran 77: The Art
        ! of Scientific Computing, 2 edition (Cambridge University Press, Cambridge England ; New York, 1992), p 246.
        !
        implicit none
        !
        integer , intent(in) :: l,m
        real(dp), intent(in) :: x(:)
        real(dp) :: y(size(x))
        integer  :: ll
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
            real(dp) :: arth(n)
            integer  :: k,k2
            real(dp) :: temp
            !
            if (n > 0) arth(1)=first
            if (n <= 16) then
                do k=2,n
                    arth(k)=arth(k-1)+increment
                enddo
            else
                do k=2,8
                    arth(k)=arth(k-1)+increment
                end do
                temp=increment*8
                k=8
                do
                    if (k >= n) exit
                    k2=k+k
                    arth(k+1:min(k2,n))=temp+arth(1:min(k,n-k))
                    temp=temp+temp
                    k=k2
                enddo
            end if
        end function arth
    end function  legendre

    pure function laguerre(k,p,x) result(y)
        ! associated Laguerre polynomial L_k^p(x)
        ! Using the expression from Samuel Shaw Ming Wong "Computational Methods in Physics and Engineering", Eq 4-86 p 139
        ! Note, there is a typographical error on http://mathworld.wolfram.com/AssociatedLaguerrePolynomial.html Eq 10
        ! Also see Linus Pauling "Introuction to Quantum Mechancs"
        implicit none 
        !
        integer , intent(in) :: k, p
        real(dp), intent(in) :: x(:)
        real(dp), allocatable :: y(:)
        integer :: j
        !
        allocate(y(size(x)))
        y=0
        !
        do j = 0, k
            y = y + nchoosek(k+p,k-j) * (-x)**j / factorial(j)
        enddo
        !
    end function  laguerre

    pure function nchoosek(n,k) result (res)
        !
        implicit none
        !
        integer, intent(in) :: n
        integer, intent(in) :: k
        integer :: res
        !
        res=factorial(n)/(factorial(k)*factorial(n-k))
        !
    end function  nchoosek

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

    pure function eye_nxm(n) result(id)
        !> nxm identity matrix
        implicit none
        !
        integer, intent(in) :: n(2)
        real(dp), dimension(:,:), allocatable :: id
        integer :: i
        !
        allocate(id(n(1),n(2)))
        id=0.0_dp
        do i = 1, minval(n)
            id(i,i) = 1.0_dp
        enddo
    end function  eye_nxm

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

    pure function ones_nxm(n) result(M)
        !> nxn identity matrix
        implicit none
        !
        integer, intent(in) :: n(2)
        real(dp), dimension(:,:), allocatable :: M
        !
        allocate(M(n(1),n(2)))
        M=1.0_dp
        !
    end function  ones_nxm

    pure function heavi(m)
        !
        ! A. V. Podolskiy and P. Vogl, Phys. Rev. B 69, 233101 (2004). Eq 15
        !
        implicit none
        !
        integer, intent(in) :: m
        real(dp) :: heavi
        !
        if (m.ge.0) then
            heavi = 1.0_dp
        else
            heavi = 0.0_dp
        endif
        !
    end function  heavi

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
        logical :: found
        !
        k=1
        wrkspace(:,1) = A(:,1)
        !
        do i = 1, size(A,2)
            !
            found = .false.
            !
            do j = 1, k
            if ( all(A(:,i).eq.wrkspace(:,j)) ) then
                found = .true.
                exit
            endif
            enddo
            !
            if (.not.found) then
                k = k + 1
                wrkspace(:,k) = A(:,i)
            endif
        enddo 
        !
        allocate(B,source=wrkspace(:,1:k))
        !
    end function  unique_columns_integer

    pure function unique_integer(A) result(B)
        !> returns unique values of integer array A(:)
        implicit none
        !
        integer, intent(in) :: A(:)
        integer, allocatable :: wrkspace(:)
        integer, allocatable :: B(:)
        integer :: i, k, n
        !
        n = size(A)
        !
        allocate(wrkspace(n))
        wrkspace=0
        !
        !
        k=0
        do i = 1, n
        if (.not. any(A(i).eq.wrkspace)) then
            k=k+1
            wrkspace(k) = A(i)
        endif
        enddo 
        !
        allocate(B,source=wrkspace(1:k))
        !
    end function  unique_integer

end module am_matlab
