#:include "fypp_macros.fpp"
module am_matlab

    use am_constants
    use am_mkl, only : am_zheev, inv, rand, am_zgeev, am_dgeev, null_svd, orth_svd

    implicit none
    
    public

    #:if COLOR > 0
        logical, private :: use_escape_codes = .true.
    #:else
        logical, private :: use_escape_codes = .false.
    #:endif

    interface adjoint
        module procedure zadjoint, dadjoint
    end interface ! adjoint
    
    interface linspace
        module procedure linspace_double, linspace_integer
    end interface ! linspace

    interface diag
        module procedure ddiag1, ddiag2, zdiag1, zdiag2
    end interface ! diag

    interface eye
        module procedure eye, eye_nxm
    end interface ! eye

    interface ones
        module procedure ones, ones_nxm
    end interface ! ones

    interface zeros
        module procedure zeros_nxm, zeros
    end interface ! zeros

    interface cross_product
        module procedure :: d_cross_product, z_cross_product
    end interface ! cross_product

    interface meshgrid
        module procedure dmeshgrid, imeshgrid
    end interface ! meshgrid

    interface fftshift
        module procedure i_fftshift, d_fftshift
    end interface ! fftshift

    interface rot2irrep
        module procedure rot2irrep_l, rot2irrep_j
    end interface ! rot2irrep
    
    interface trace
        module procedure dtrace, ztrace
    end interface ! trace

    interface isequal
        module procedure r0_isequal, r1_isequal, r2_isequal, &
                         c0_isequal, c1_isequal, c2_isequal, &
                         i0_isequal, i1_isequal, i2_isequal, &
                         l0_isequal, l1_isequal, l2_isequal
    end interface ! isequal
    
    interface trim_null
        module procedure i_trim_null, vi_trim_null, vz_trim_null, vd_trim_null, vs_trim_null, s_trim_null
    end interface ! trim_null

    #:setvar KINDS_issubset ['real(dp)', 'complex(dp)', 'integer']
    #:setvar RANKS_issubset range(1,3)
    #:setvar PREFS_issubset ['{}{}_'.format(KIND[3], RANK) for KIND in KINDS_issubset for RANK in RANKS_issubset]

    interface issubset
        module procedure ${variants('issubset', prefixes=PREFS_issubset)}$
    end interface ! issubset

    #:setvar KINDS_reallocate ['real(dp)', 'complex(dp)', 'integer', 'logical']
    #:setvar RANKS_reallocate range(1,4)
    #:setvar PREFS_reallocate ['{}{}_'.format(KIND[3], RANK) for KIND in KINDS_reallocate for RANK in RANKS_reallocate]

    interface reallocate
        module procedure ${variants('reallocate', prefixes=PREFS_reallocate)}$
    end interface ! reallocate

    #:setvar KINDS_vector_allocate ['real(dp)', 'complex(dp)', 'integer', 'logical', 'character(:)']
    #:setvar RANKS_vector_allocate range(1,4)
    #:setvar PREFS_vector_allocate ['{}{}_'.format(KIND[3], RANK) for KIND in KINDS_vector_allocate for RANK in RANKS_vector_allocate]

    interface vector_allocate
        module procedure ${variants('vector_allocate', prefixes=PREFS_vector_allocate)}$
    end interface ! vector_allocate

    #:setvar KINDS_write_xml_attribute ['real(dp)', 'complex(dp)', 'integer', 'logical', 'character(*)']
    #:setvar RANKS_write_xml_attribute range(0,4)
    #:setvar PREFS_write_xml_attribute ['{}{}_'.format(KIND[3], RANK) for KIND in KINDS_write_xml_attribute for RANK in RANKS_write_xml_attribute]

    interface write_xml_attribute
        module procedure ${variants('write_xml_attribute', prefixes=PREFS_write_xml_attribute)}$
    end interface ! write_xml_attribute

    #:setvar KINDS_read_xml_attribute ['real(dp)', 'complex(dp)', 'integer', 'logical', 'character(*)']
    #:setvar RANKS_read_xml_attribute range(0,4)
    #:setvar PREFS_read_xml_attribute ['{}{}_'.format(KIND[3], RANK) for KIND in KINDS_read_xml_attribute for RANK in RANKS_read_xml_attribute]

    interface read_xml_attribute
        module procedure ${variants('read_xml_attribute', prefixes=PREFS_read_xml_attribute)}$
    end interface ! read_xml_attribute

    #:setvar KINDS_unique ['real(dp)', 'complex(dp)', 'integer', 'logical']
    #:setvar RANKS_unique range(1,4)
    #:setvar PREFS_unique ['{}{}_'.format(KIND[3], RANK) for KIND in KINDS_unique for RANK in RANKS_unique]

    interface unique_inds
        module procedure ${variants('unique_inds', prefixes=PREFS_unique)}$
    end interface ! unique_inds

    interface unique
        module procedure ${variants('unique', prefixes=PREFS_unique)}$
    end interface ! unique

    interface cumsum
        module procedure dcumsum, icumsum
    end interface ! cumsum

    interface cumprod
        module procedure dcumprod, icumprod
    end interface ! cumprod

    interface foursort
        module procedure d_foursort, i_foursort
    end interface ! foursort

    interface fourrank
        module procedure d_fourrank, i_fourrank
    end interface ! fourrank

    interface rref
        module procedure :: z_rref, d_rref
    end interface ! rref

    interface eigenspace
        module procedure :: z_eigenspace, d_eigenspace
    end interface ! eigenspace

    interface subspace_intersection
        module procedure z_subspace_intersection, d_subspace_intersection
    end interface ! subspace_intersection

    contains

    ! basic print routine for title

    pure function  centertitle(title,length) result(title_centered)
        !
        implicit none
        !
        character(len=*), intent(in) :: title
        integer, intent(in) :: length
        character(len=length) :: title_centered
        integer :: blanks
        integer :: i 
        !
        title_centered = title
        blanks = 0 
        do  i = length,1,-1                      ! starting on the right side 
            if (title_centered(i:i).ne.' ') then ! check for trailing blanks 
                exit
            endif 
            blanks = blanks + 1                     ! count the blanks 
        enddo 
        !
        if ( blanks .gt. 1 ) then 
            blanks = blanks/2                       ! cut half of the blanks 
            do i = length , 1 , -1 
              if ( i-blanks .gt. 0 ) title_centered(i:i) = title_centered(i-blanks:i-blanks) 
              if ( i-blanks .le. 0 ) title_centered(i:i) = ' '     ! blank front half 
            enddo
        endif
    end function   centertitle

    subroutine     print_title(title)
        !
        implicit none
        !
        character(len=*), intent(in) :: title
        integer :: fid
        !
        fid = 6
        !
        write(unit=fid,fmt='(a,a,a)') style('reserve,bright,txt:blue')// "+", repeat("-",column_width-2), "+"
        write(unit=fid,fmt='(a,a,a)') "|", style('bright,txt:yellow')//centertitle(title,column_width-2)//style('bright,txt:blue'), "|"
        write(unit=fid,fmt='(a,a,a)') "+", repeat("-",column_width-2), "+"// style('clear')
        !
    end subroutine print_title

    function       style(fmt) result(seq_out)
        !
        implicit none
        !
        character(*), intent(in) :: fmt
        character(:), allocatable :: seq_out
        character(256) :: seq
        ! start sequence with command to clear any previous attributes
        if (use_escape_codes) then
            seq = char(27)//'[0'
            if (index(fmt,'clear')      .ne.0) seq = trim(seq)//';0'
            if (index(fmt,'bright')     .ne.0) seq = trim(seq)//';1'
            if (index(fmt,'faint')      .ne.0) seq = trim(seq)//';2'
            if (index(fmt,'italic')     .ne.0) seq = trim(seq)//';3'
            if (index(fmt,'underline')  .ne.0) seq = trim(seq)//';4'
            if (index(fmt,'blink_slow') .ne.0) seq = trim(seq)//';5'
            if (index(fmt,'blink_fast') .ne.0) seq = trim(seq)//';6'
            if (index(fmt,'txt:black')  .ne.0) seq = trim(seq)//';30'
            if (index(fmt,'txt:red')    .ne.0) seq = trim(seq)//';31'
            if (index(fmt,'txt:green')  .ne.0) seq = trim(seq)//';32'
            if (index(fmt,'txt:yellow') .ne.0) seq = trim(seq)//';33'
            if (index(fmt,'txt:blue')   .ne.0) seq = trim(seq)//';34'
            if (index(fmt,'txt:magenta').ne.0) seq = trim(seq)//';35'
            if (index(fmt,'txt:cyan')   .ne.0) seq = trim(seq)//';36'
            if (index(fmt,'txt:white')  .ne.0) seq = trim(seq)//';37'
            if (index(fmt,'bg:black')   .ne.0) seq = trim(seq)//';40'
            if (index(fmt,'bg:red')     .ne.0) seq = trim(seq)//';41'
            if (index(fmt,'bg:green')   .ne.0) seq = trim(seq)//';42'
            if (index(fmt,'bg:yellow')  .ne.0) seq = trim(seq)//';43'
            if (index(fmt,'bg:blue')    .ne.0) seq = trim(seq)//';44'
            if (index(fmt,'bg:magenta') .ne.0) seq = trim(seq)//';45'
            if (index(fmt,'bg:cyan')    .ne.0) seq = trim(seq)//';46'
            if (index(fmt,'bg:white')   .ne.0) seq = trim(seq)//';47'
            ! append end of sequence marker
            seq = trim(seq)//'m'
            allocate(seq_out, source=trim(seq))
        else 
            allocate(seq_out, source='')
        endif
    end function   style

    pure function  lowercase(str) result(str_lwr)
        !
        implicit none
        !
        character(len=*), intent(in) :: str
        character(len=maximum_buffer_size) :: str_lwr
        !
        integer :: iA,iZ,idiff,ipos,ilett
        !
        iA = ichar('A')
        iZ = ichar('Z')
        idiff = iZ-ichar('z')
        !
        str_lwr = str
        !
        do ipos=1,len(str)
           ilett = ichar(str(ipos:ipos))
           if ((ilett.ge.iA).and.(ilett.le.iZ)) &
                str_lwr(ipos:ipos)=char(ilett-idiff)
        enddo
        !
        str_lwr = trim(adjustl(str_lwr))
        !
    end function   lowercase

    ! progress bar

    subroutine     show_progress(iteration,maximum)
        !
        implicit none
        !
        integer, intent(in) :: iteration
        integer, intent(in) :: maximum
        !

        WRITE(*,"(5x,f6.2,A,$)") iteration/real(maximum,dp)*100.0_dp, CHAR(13)
        !
    end subroutine show_progress

    ! physics

    ! angular momenta operators

    pure function  Lz(l)
        !
        ! Construct  Ly matrix in the Lz basis from raising and lower operators
        ! Laloe, p 666; Martin, p 573, Eq N4; Sakurai, p 207. hbar is set to 1.
        !
        ! R. Shankar, Principles of Quantum Mechanics, Softcover reprint of the original 1st ed.
        ! 1980 edition (Springer, 2013), p 327, Eq 12.5.21b
        !
        implicit none
        !
        integer, intent(in) :: l
        complex(dp), allocatable :: Lz(:,:)
        !
        allocate(Lz(-l:l,-l:l))
        !
        Lz = (Lm(l)-Lp(l))*0.5_dp*cmplx_i
        !
    end function   Lz

    pure function  Lx(l)
        !
        ! Construct  Ly matrix in the Lz basis from raising and lower operators
        ! Laloe, p 666; Martin, p 573, Eq N4; Sakurai, p 207. hbar is set to 1.
        !
        ! R. Shankar, Principles of Quantum Mechanics, Softcover reprint of the original 1st ed.
        ! 1980 edition (Springer, 2013), p 327, Eq 12.5.21a
        !
        implicit none
        !
        integer, intent(in) :: l
        complex(dp), allocatable :: Lx(:,:)
        !
        allocate(Lx(-l:l,-l:l))
        !
        Lx = (Lp(l)+Lm(l))*0.5_dp
        !
    end function   Lx

    pure function  Lp(l)
        !
        ! Raising angular momentum operator.
        !
        ! R. Shankar, Principles of Quantum Mechanics, Softcover reprint of the original 1st ed.
        ! 1980 edition (Springer, 2013), p 327, Eq 12.5.20
        !
        implicit none
        !
        integer, intent(in) :: l
        complex(dp), allocatable :: Lp(:,:)
        integer :: m, mp
        !
        allocate(Lp(-l:l,-l:l))
        Lp = 0.0_dp
        !
        do m = -l, l
        do mp= -l, l
            if (mp==m+1) Lp(m,mp) = sqrt(real( (l-m)*(l+m+1) ,dp))
        enddo
        enddo
        !
    end function   Lp

    pure function  Lm(l)
        !
        ! Lowering angular momentum operator.
        !
        ! R. Shankar, Principles of Quantum Mechanics, Softcover reprint of the original 1st ed.
        ! 1980 edition (Springer, 2013), p 327, Eq 12.5.20
        !
        implicit none
        !
        integer, intent(in) :: l
        complex(dp), allocatable :: Lm(:,:)
        integer :: m, mp
        !
        allocate(Lm(-l:l,-l:l))
        Lm = 0.0_dp
        !
        do m = -l, l
        do mp= -l, l
            if (mp==m-1) Lm(m,mp) = sqrt(real( (l+m)*(l-m+1) ,dp))
        enddo
        enddo
        !
    end function   Lm

    function       Jz(j)
        !
        ! Construct  Jy matrix in the Jz basis from raising and lower operators
        ! Jaloe, p 666; Martin, p 573, Eq N4; Sakurai, p 207. hbar is set to 1.
        !
        ! R. Shankar, Principles of Quantum Mechanics, Softcover reprint of the original 1st ed.
        ! 1980 edition (Springer, 2013), p 327, Eq 12.5.21b
        !
        implicit none
        !
        real(dp), intent(in) :: j
        complex(dp), allocatable :: Jz(:,:)
        integer :: n
        !
        if (abs(modulo(j+tiny,0.5_dp)-tiny).gt.tiny) stop 'Quantum number j invalid.'
        !
        ! dimension of irreducible matrix representation
        n = nint(2.0_dp*j+1.0_dp)
        !
        allocate(Jz(1:n,1:n))
        Jz = 0.0_dp
        !
        Jz = (Jm(j)-Jp(j))*0.5_dp*cmplx_i
        !
    end function   Jz
    
    function       Jx(j)
        !
        ! Construct  Jy matrix in the Jz basis from raising and lower operators
        ! Jaloe, p 666; Martin, p 573, Eq N4; Sakurai, p 207. hbar is set to 1.
        !
        ! R. Shankar, Principles of Quantum Mechanics, Softcover reprint of the original 1st ed.
        ! 1980 edition (Springer, 2013), p 327, Eq 12.5.21a
        !
        implicit none
        !
        real(dp), intent(in) :: j
        complex(dp), allocatable :: Jx(:,:)
        integer :: n
        !
        if (abs(modulo(j+tiny,0.5_dp)-tiny).gt.tiny) stop 'Quantum number j invalid.'
        !
        ! dimension of irreducible matrix representation
        n = nint(2.0_dp*j+1.0_dp)
        !
        allocate(Jx(1:n,1:n))
        Jx = 0.0_dp
        !
        Jx = (Jp(j)+Jm(j))*0.5_dp
        !
    end function   Jx

    function       Jp(j)
        !
        ! Raising angular momentum operator.
        !
        ! R. Shankar, Principles of Quantum Mechanics, Softcover reprint of the original 1st ed.
        ! 1980 edition (Springer, 2013), p 327, Eq 12.5.20
        !
        implicit none
        !
        real(dp), intent(in) :: j
        complex(dp), allocatable :: Jp(:,:)
        real(dp), allocatable :: mlist(:)
        integer :: m, mp, n, i
        !
        if (abs(modulo(j+tiny,0.5_dp)-tiny).gt.tiny) stop 'Quantum number j invalid.'
        !
        ! dimension of irreducible matrix representation
        n = nint(2.0_dp*j+1.0_dp)
        !
        ! create mlist
        allocate(mlist(n))
        do i = 1, n
            mlist(i) =  -j + real(i-1,dp)
        enddo
        !
        allocate(Jp(1:n,1:n))
        Jp = 0.0_dp
        !
        do m = 1, n
        do mp= 1, n
            if (mp==m+1) Jp(m,mp) = sqrt(real( (j-mlist(m))*(j+mlist(m)+1) ,dp))
        enddo
        enddo
        !
    end function   Jp

    function       Jm(j)
        !
        ! Lowering angular momentum operator.
        !
        ! R. Shankar, Principles of Quantum Mechanics, Softcover reprint of the original 1st ed.
        ! 1980 edition (Springer, 2013), p 327, Eq 12.5.20
        !
        implicit none
        !
        real(dp), intent(in) :: j
        complex(dp), allocatable :: Jm(:,:)
        real(dp), allocatable :: mlist(:)
        integer :: m, mp, n, i
        !
        !
        if (abs(modulo(j+tiny,0.5_dp)-tiny).gt.tiny) stop 'Quantum number j invalid.'
        !
        ! dimension of irreducible matrix representation
        n = nint(2.0_dp*j+1.0_dp)
        !
        ! create mlist
        allocate(mlist(n))
        do i = 1, n
            mlist(i) =  -j + real(i-1,dp)
        enddo
        !
        allocate(Jm(1:n,1:n))
        Jm = 0.0_dp
        !
        do m = 1, n
        do mp= 1, n
            if (mp==m-1) Jm(m,mp) = sqrt(real( (j+mlist(m))*(j-mlist(m)+1) ,dp))
        enddo
        enddo
        !
    end function   Jm

    ! special functions

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
        real(dp), allocatable :: y(:)
        integer  :: ll
        real(dp) :: pll(size(x))
        real(dp) :: pmm(size(x))
        real(dp) :: pmmp1(size(x))
        real(dp) :: somx2(size(x))
        !
        allocate(y(size(x)))
        !
        if ((m.lt.0).or.(m.gt.l).or.(all(abs(x).gt.1.0_dp))) then
            y = 0
        endif
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
        ! associated Laguerre polynomial L_k^p(x)(k,p,x)
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
        allocate(y,mold=x)
        y=0.0_dp
        !
        do j = 0, k
            y = y + nchoosek(k+p,k-j) * (-x)**j / factorial(j)
        enddo
        !
    end function  laguerre

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

    ! rotation functions

    elemental subroutine correct_rounding_error(sym)
        !
        implicit none
        !
        real(dp), intent(inout) :: sym
        !
        ! correct rounding error in cart
        ! Possible values are: cos(pi/n), sin(pi/n) for n = 1, 2, 3, 6
        ! Symmetry and Condensed Matter Physics: A Computational Approach. 1 edition. 
        ! Cambridge, UK ; New York: Cambridge University Press, 2008. page 275.
        !           n  =      1         2         3         6
        !    sin(pi/n) =   0.0000    1.0000    0.8660    0.5000
        !    cos(pi/n) =  -1.0000    0.0000    0.5000    0.8660 - sqrt(3)/2
        !
        if    (abs(sym - nint(sym)           ).lt.tiny) then; sym = nint(sym)            ! +1,-1,0
        elseif(abs(sym + 0.500000000000000_dp).lt.tiny) then; sym =-0.500000000000000_dp ! -0.5
        elseif(abs(sym - 0.500000000000000_dp).lt.tiny) then; sym = 0.500000000000000_dp ! +0.5
        elseif(abs(sym + 1.154700538379252_dp).lt.tiny) then; sym =-1.154700538379252_dp ! -sqrt(3)*(2/3)
        elseif(abs(sym - 1.154700538379252_dp).lt.tiny) then; sym = 1.154700538379252_dp ! +sqrt(3)*(2/3)
        elseif(abs(sym + 0.866025403784439_dp).lt.tiny) then; sym =-0.866025403784439_dp ! -sqrt(3)*(1/2)
        elseif(abs(sym - 0.866025403784439_dp).lt.tiny) then; sym = 0.866025403784439_dp ! +sqrt(3)*(1/2)
        elseif(abs(sym + 1.732050807568877_dp).lt.tiny) then; sym =-1.732050807568877_dp ! -sqrt(3)
        elseif(abs(sym - 1.732050807568877_dp).lt.tiny) then; sym = 1.732050807568877_dp ! +sqrt(3)
        endif
        !
    end subroutine       correct_rounding_error

    pure function vec2dcosines(vec) result(dcosines)
        !
        implicit none
        !
        real(dp), intent(in) :: vec(3)
        real(dp) :: dcosines(3)
        !
        if (norm2(vec).gt.tiny) then
            dcosines = vec/norm2(vec)
        else
            dcosines = real([0,0,1],dp)
        endif
    end function  vec2dcosines

    pure function dcosines2euler(dcosines) result(euler)
        !
        ! returns the euler angles which correspond to the rotation which alignes the dcosines vector with z 
        !
        ! R = X(alpha) * Y(beta) * Z(gamma)
        ! euler = [alpha, beta, gamma]
        !
        implicit none
        !
        real(dp), intent(in) :: dcosines(3)
        real(dp) :: euler(3)
        ! 
        !
        ! angle is not used in spherical coordinates
        euler(1) = 0.0_dp
        ! phi  (azimuthal angle) in spherical coordinates
        euler(2) = atan(dcosines(2)/dcosines(1)+1.0D-14)
        ! theta (polar angle) in spherical coordinates
        euler(3) = acos(dcosines(3))
        !
    end function  dcosines2euler

    pure function rot2axis_angle(R) result(aa)
        !
        implicit none
        !
        real(dp), intent(in) :: R(3,3)
        real(dp) :: Rcp(3,3)
        real(dp) :: aa(4) ! axis & angle [x,y,z,th]
        real(dp) :: tr, phi, axis(3), d
        integer  :: i
        !
        ! if R is a rotoinversion, get the angle and axis of the rotational part only (without the inversion)
        d = R(1,1)*R(2,2)*R(3,3)-R(1,1)*R(2,3)*R(3,2)-R(1,2)*R(2,1)*R(3,3)+R(1,2)*R(2,3)*R(3,1)+R(1,3)*R(2,1)*R(3,2)-R(1,3)*R(2,2)*R(3,1)
        Rcp = R*sign(1.0_dp,d)
        !
        tr = trace(Rcp)
        !
        if (abs(tr-3.0_dp).lt.tiny) then
            aa = real([0,1,0,0],dp)
        elseif (abs(tr+1.0_dp).lt.tiny) then
            do i = 1,3
                axis(i) = sqrt(max(0.5_dp*(Rcp(i,i)+1),0.0_dp))
            enddo
            aa(1:3) = axis/norm2(axis)
            aa(4)   = pi
        else
            phi = acos( (tr-1.0_dp)/2.0_dp )
            axis = [Rcp(3,2)-Rcp(2,3),Rcp(1,3)-Rcp(3,1),Rcp(2,1)-Rcp(1,2)]/(2.0_dp*sin(phi))
            aa(1:3) = axis/norm2(axis)
            aa(4)   = phi
        endif
        !
    end function  rot2axis_angle

    pure function axis_angle2rot(aa) result(R)
        !
        ! Note the different sign convention from matlab. aa(4) here = -aa(4) in matlab.
        ! 
        ! When aa = [0,0,1,45/180*pi] ~ 45 deg rotation around z axis in positive direction.
        !   x -> (2^(1/2)*x)/2 + (2^(1/2)*y)/2
        !   y -> (2^(1/2)*y)/2 - (2^(1/2)*x)/2
        !   z ->                             z
        !
        implicit none
        !
        real(dp), intent(in) :: aa(4) ! axis & angle [x,y,z,th]
        real(dp) :: R(3,3)
        real(dp) :: s, c, t, n(3), x, y, z
        !
        s = sin(aa(4))
        c = cos(aa(4))
        t = 1.0_dp - c
        n = aa(1:3)/(norm2(aa(1:3))+1.0D-14)
        x = n(1)
        y = n(2)
        z = n(3)
        R(:,1) = [ t*x*x + c,   t*x*y - s*z, t*x*z + s*y ]
        R(:,2) = [ t*x*y + s*z, t*y*y + c,   t*y*z - s*x ]
        R(:,3) = [ t*x*z - s*y, t*y*z + s*x, t*z*z + c   ]
        !
    end function  axis_angle2rot

    pure function rotmat(A,B) result(R)
        !
        ! Rotation matrix R which aligns unit vector A to unit vector B: B = rotmat(A,B) * A
        !
        ! based on Rodrigues' Rotation Formula
        ! "A Mathematical Introduction to Robotic Manipulation", Richard M. Murray, Zexiang Li, S. Shankar Sastry, pp. 26-28
        !
        ! % quick matlab implementation:
        ! A = rand(3,1); A=A./norm(A);
        ! B = rand(3,1); B=B./norm(A);
        ! v = cross(A,B);
        ! s = norm(v);
        ! c = dot(A,B);
        ! vx(1:3,1) = [0.0,v(3),-v(2)];
        ! vx(1:3,2) = [-v(3),0.0,v(1)];
        ! vx(1:3,3) = [v(2),-v(1),0.0];
        ! R = eye(3) + vx + (vx*vx)*(1.0-c)/(s.^2)
        ! R*A-B
        !
        implicit none
        !
        real(dp), intent(in) :: A(3)
        real(dp), intent(in) :: B(3)
        real(dp) :: R(3,3)
        real(dp) :: s ! sine of angle
        real(dp) :: c ! cosine of angle
        real(dp) :: v(3) ! cross product of A and B
        real(dp) :: vx(3,3) ! skew symmetric cross product of v
        ! 
        v = cross_product(A,B)
        s = norm2(v)
        c = dot_product(A,B)
        vx(1:3,1) = [0.0_dp,v(3),-v(2)]
        vx(1:3,2) = [-v(3),0.0_dp,v(1)]
        vx(1:3,3) = [v(2),-v(1),0.0_dp]
        !
        R = eye(3) + vx + matmul(vx,vx)*(1.0_dp - c)/(s**2)
        !
    end function  rotmat

    pure function euler2rot(euler) result(R)
        !
        ! Tait–Bryan "pitch-roll-yaw" convetion
        ! R = X(alpha) * Y(beta) * Z(gamma)
        ! euler = [alpha, beta, gamma]
        !
        implicit none
        !
        real(dp), intent(in) :: euler(3)
        real(dp) :: X(3,3), Z(3,3), Y(3,3)
        real(dp) :: R(3,3)
        !
        X = axis_angle2rot([1.0_dp, 0.0_dp, 0.0_dp, euler(1)]) ! rot around X
        Y = axis_angle2rot([0.0_dp, 1.0_dp, 0.0_dp, euler(2)]) ! rot around Y
        Z = axis_angle2rot([0.0_dp, 0.0_dp, 1.0_dp, euler(3)]) ! rot around Z
        !
        R = matmul(X,matmul(Y,Z))
        !
    end function  euler2rot

    pure function rot2euler(R) result(euler)
        !
        ! Tait–Bryan "pitch-roll-yaw" convetion
        ! R = X(alpha) * Y(beta) * Z(gamma)
        ! euler = [alpha, beta, gamma]
        ! From : Extracting Euler Angles from a Rotation Matrix
        !        Mike Day, Insomniac Games
        implicit none
        !
        real(dp), intent(in) :: R(3,3)
        real(dp) :: euler(3), s, c
        integer :: i,j,k
        !
        i = 1
        j = 2
        k = 3
        !
        euler(i) = atan2(R(j,k),R(k,k))
        euler(j) = atan2(-R(i,k),sqrt(R(i,i)**2+R(i,j)**2))
        s = sin(euler(1))
        c = cos(euler(1))
        euler(k) = atan2(s*R(k,i)-c*R(j,i),c*R(j,j)-s*R(k,j))
        !
    end function  rot2euler

    function      euler2SO3(l,euler) result(SO3)
        !
        ! The group SO(3) is the set of all three dimensional, real orthogonal matrices with unit
        ! determinant. [Requiring that the determinant equal 1 and not −1, excludes inversions.]
        !
        ! The group SO(3) represents the set of all possible rotations of a three dimensional real
        ! vector; it is essentially the collection of all proper rotations. The orthogonal group O(3) 
        ! includes improper rotations and it is given by the direct product of SO(3) and the inversion
        ! group i.
        !
        ! setup rotation matrices, "pitch-roll-yaw" convetion
        ! alpha (           around X)
        ! beta  (polar,     around Y)
        ! gamma (azimuthal, around Z)
        !
        ! R = X(alpha) * Y(beta) * Z(gamma)
        !
        ! Determines 3D irrep rotation by calculating the eigenvalues D and eigevectors V of Ly/Lx
        ! operator in the Lz basis using raising/lowering operators.
        !
        !   - R. M. Martin, Electronic Structure: Basic Theory and Practical Methods, 1 edition
        !        (Cambridge University Press, Cambridge, UK ; New York, 2008), p 573.
        !   - Romero Nichols "Density Functional Study of Fullerene-Based Solids: Crystal Structure,
        !        Doping, and Electron- Phonon Interaction", Ph.D. Thesis UIUC, p 96.
        !   - C. Cohen-Tannoudji, B. Diu, and F. Laloe, Quantum Mechanics, 1 edition (Wiley-VCH, New
        !        York; Paris, 1992), p 666.
        !   - J. J. Sakurai, Modern Quantum Mechanics, Revised edition (Addison Wesley, Reading, Mass,
        !        1993). p 207
        !
        implicit none
        !
        integer , intent(in) :: l
        real(dp), intent(in) :: euler(3) ! around X, Y, Z
        integer :: n
        real(dp)   , allocatable :: SO3(:,:) ! tesseral harmonics rotation matrix, (2*l+1) x (2*l+1) irrep of rotation in 3 dimensional space
        complex(dp), allocatable :: V(:,:) ! eigenvector of Lz/Lx in Ly basis
        real(dp)   , allocatable :: D(:)   ! eigenvalues of Lz/Lx in Ly basis = -l, -l+1, ... -1, 0, 1, ... l-1, l  -- Note: Lz/Lx are Hermitian.
        complex(dp), allocatable :: A(:,:) ! matrix corresponding to exp(-i*alpha*Ly)
        complex(dp), allocatable :: B(:,:) ! matrix corresponding to exp(-i*beta*Ly)
        complex(dp), allocatable :: G(:,:) ! matrix corresponding to exp(-i*gamma*Ly)
        complex(dp), allocatable :: C(:,:) ! similarity transform which maps complex spherical harmonics onto real tesseral harmonics
        real(dp)   , allocatable :: X(:,:) ! counter-clockwise rotation (right-hand rule) around X axis
        real(dp)   , allocatable :: Y(:,:) ! counter-clockwise rotation (right-hand rule) around Y axis
        real(dp)   , allocatable :: Z(:,:) ! counter-clockwise rotation (right-hand rule) around Z axis
        complex(dp), allocatable :: H(:,:) ! helper matrix
        complex(dp), allocatable :: Q(:,:) ! helper matrix
        ! get dimensions of matrices
        n = 2*l+1
        ! determine similarity transform to convert spherical into tesseral harmonics (complex to real)
        C = tesseral(l=l)
        ! if this Q is not set deliberately, sometimes the code will crash
        Q = Lx(l=l)
        ! determine rotation around X = real( (B*V)*A*(B*V)' )
        call am_zheev(A=Q,V=V,D=D)
        H = matmul(C,V)
        A = diag(exp(-cmplx_i*euler(1)*D))
        X = matmul(H,matmul(A,adjoint(H)))
        ! determine rotation around Y = real( B*P*B' )
        allocate(Y(n,n))
        B = diag(exp( cmplx_i*euler(2)*D)) ! there should be a minus sign here...
        Y = matmul(C,matmul(B,adjoint(C)))
        ! determine rotation around Z = real( (B*V)*T*(B*V)' )
        call am_zheev(A=Lz(l),V=V,D=D)
        H = matmul(C,V)
        G = diag(exp(-cmplx_i*euler(3)*D))
        Z = matmul(H,matmul(G,adjoint(H)))
        ! generate rotation
        allocate(SO3(n,n))
        SO3 = matmul(X,matmul(Y,Z))
        !        
        contains
        pure function  tesseral(l) result(B)
            !
            ! latex equations
            !
            ! $$
            ! Y_{lm} = 
            ! \begin{cases}
            ! Y_l^m                                               & \text{if } m = 0 \\
            ! \frac{i}{\sqrt{2}} ( (-1)^{m} Y_l^{-m} - Y_l^{ m} ) & \text{if } m > 0 \\
            ! \frac{1}{\sqrt{2}} ( (-1)^{m} Y_l^{ m} + Y_l^{-m} ) & \text{if } m < 0
            ! \end{cases}
            ! $$
            !
            ! R. R. Sharma, Phys. Rev. B. 19, 2813 (1979).
            ! Also, see Wikipedia.
            implicit none
            !
            integer, intent(in) :: l
            complex(dp), allocatable :: B(:,:)
            integer :: m,mp
            !
            allocate(B(-l:l,-l:l))
            B = 0.0_dp
            !
            do m = -l, l
            do mp= -l, l
                if (m.eq.0) then
                   if (m.eq. mp) B(m,mp) = 1.0_dp
                elseif (m.gt.0) then
                   if (m.eq. mp) B(m,mp) = -cmplx_i/sqrt(2.0_dp)
                   if (m.eq.-mp) B(m,mp) = +cmplx_i/sqrt(2.0_dp) * (-1.0_dp)**m
                elseif (m.lt.0) then
                   if (m.eq. mp) B(m,mp) =   1.0_dp/sqrt(2.0_dp) * (-1.0_dp)**m
                   if (m.eq.-mp) B(m,mp) =   1.0_dp/sqrt(2.0_dp)
                endif
            enddo
            enddo
            !
        end function   tesseral
    end function  euler2SO3

    function      euler2SU2(j,euler) result(SU2)
        !
        implicit none
        !
        real(dp), intent(in) :: j
        real(dp), intent(in) :: euler(3) ! around X, Y, Z
        integer :: n
        complex(dp), allocatable :: SU2(:,:) ! tesseral harmonics rotation matrix, (2*j+1) x (2*j+1) irrep of rotation in 3 dimensional space
        complex(dp), allocatable :: V(:,:)   ! eigenvector of Jz/Jx in Jy basis
        real(dp)   , allocatable :: D(:)     ! eigenvalues of Jz/Jx in Jy basis = -j, -j+1, ... -1, 0, 1, ... j-1, j  -- Note: Jz/Jx are Hermitian.
        complex(dp), allocatable :: A(:,:)   ! matrix corresponding to exp(-i*alpha*Jy)
        complex(dp), allocatable :: B(:,:)   ! matrix corresponding to exp(-i*beta*Jy)
        complex(dp), allocatable :: G(:,:)   ! matrix corresponding to exp(-i*gamma*Jy)
        complex(dp), allocatable :: X(:,:)   ! counter-clockwise rotation (right-hand rule) around X axis
        complex(dp), allocatable :: Y(:,:)   ! counter-clockwise rotation (right-hand rule) around Y axis
        complex(dp), allocatable :: Z(:,:)   ! counter-clockwise rotation (right-hand rule) around Z axis
        ! get dimensions of matrices
        n = nint(2.0_dp*j+1.0_dp)
        ! determine rotation around X = real( (B*V)*A*(B*V)' )
        call am_zheev(A=Jx(j),V=V,D=D)
        A = diag(exp(-cmplx_i*euler(1)*D ))
        X = matmul(V,matmul(A,adjoint(V)))
        ! determine rotation around Y = real( B*P*B' )
        allocate(Y(n,n))
        B = diag(exp( cmplx_i*euler(2)*D ))
        Y = B
        ! determine rotation around Z = real( (B*V)*T*(B*V)' )
        call am_zheev(A=Jz(j),V=V,D=D)
        G = diag(exp( cmplx_i*euler(3)*D ))
        Z = matmul(V,matmul(G,adjoint(V)))
        ! generate rotation
        allocate(SU2(n,n))
        SU2 = matmul(X,matmul(Y,Z))
        !
    end function  euler2SU2

    function      rot2O3(l,R) result(O3)
        !
        implicit none
        !
        integer, intent(in) :: l
        real(dp), intent(in) :: R(3,3)
        real(dp), allocatable :: O3(:,:)
        real(dp) :: d ! det
        !
        ! computed determinant
        d = R(1,1)*R(2,2)*R(3,3)-R(1,1)*R(2,3)*R(3,2)-R(1,2)*R(2,1)*R(3,3)+R(1,2)*R(2,3)*R(3,1)+R(1,3)*R(2,1)*R(3,2)-R(1,3)*R(2,2)*R(3,1)
        !
        ! check that the rotation is unitary, with determinant equalt to plus or minus one
        if (abs(abs(d)-1.0_dp).gt.tiny) stop 'Rotation is not unitary. det /= +1 or -1'
        !
        ! remove rouding errors
        d = sign(1.0_dp,d)
        !
        ! convert rotoinversion to pure rotation then back to rotoinversion
        O3 = euler2SO3(l=l,euler=rot2euler(R*d)) * ( d ** l )
        !
    end function  rot2O3

    function      rot2irrep_l(l,R) result(irrep)
        ! generates n-dimensional irreducible representation matrix corresponding to (im)proper rotation R 
        implicit none
        !
        integer , intent(in) :: l
        real(dp), intent(in) :: R(3,3)
        complex(dp), allocatable :: irrep(:,:)
        integer  :: improper_fac
        real(dp) :: d ! det
        ! convert rotoinversion to pure rotation then back to rotoinversion
        d = R(1,1)*R(2,2)*R(3,3)-R(1,1)*R(2,3)*R(3,2)-R(1,2)*R(2,1)*R(3,3)+R(1,2)*R(2,3)*R(3,1)+R(1,3)*R(2,1)*R(3,2)-R(1,3)*R(2,2)*R(3,1)
        ! check that the rotation is unitary, with determinant equalt to plus or minus one
        if (abs(abs(d)-1.0_dp).gt.tiny) stop 'Rotation is not unitary. det /= +1 or -1'
        ! remove rouding errors
        d = nint(sign(1.0_dp,d))
        ! rotoinversions factor : (-1)**l <= this corresponds to inversion matrix, use d instead of conditional statement
        improper_fac = (d)**l
        ! converts (im)rotation -> proper rotation -> builds S03 rpresentation -> tack on (im)proper factor to recover inversion
        irrep = euler2SO3(l=l, euler=rot2euler(R*d)) * improper_fac
        !
    end function  rot2irrep_l

    function      rot2irrep_j(j,R) result(irrep)
        ! generates n-dimensional irreducible representation matrix corresponding to (im)proper rotation R 
        implicit none
        !
        real(dp), intent(in) :: j
        real(dp), intent(in) :: R(3,3)
        complex(dp), allocatable :: irrep(:,:)
        real(dp) :: d ! det
        ! check that j is a half or full integer
        if (mod(j+tiny,0.5_dp)-tiny.gt.tiny) then
            stop 'ERROR [rot2irrep_j]: j is not a half/full integer'
        endif
        ! convert rotoinversion to pure rotation then back to rotoinversion
        d = R(1,1)*R(2,2)*R(3,3)-R(1,1)*R(2,3)*R(3,2)-R(1,2)*R(2,1)*R(3,3)+R(1,2)*R(2,3)*R(3,1)+R(1,3)*R(2,1)*R(3,2)-R(1,3)*R(2,2)*R(3,1)
        ! check that the rotation is unitary, with determinant equalt to plus or minus one
        if (abs(abs(d)-1.0_dp).gt.tiny) stop 'Rotation is not unitary. det /= +1 or -1'
        ! remove rouding errors
        d = nint(sign(1.0_dp,d))
        ! Improper rotations in SU(2): Multiply the matrix elements by +1.
        ! pg 58 of Simon L. Altmann, Peter Herzig "Point-Group Theory Tables".
        ! Also see: Altmann, S. L. (1986). Rotations, quaternions, and double groups. Clarendon Press, Oxford.
        irrep = euler2SU2(j=j, euler=rot2euler(R*d))
        !
    end function  rot2irrep_j

    ! statistics

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

    ! file io functions

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

    function      fexists(fname)
        !
        implicit none
        !
        character(*) :: fname
        logical :: fexists
        !
        inquire(file=trim(fname),exist=fexists)
        !
    end function  fexists

    pure function strsplit(str,delimiter) result(word)
        !
        character(*), intent(in) :: str
        character(1), intent(in) :: delimiter
        character(len=:), allocatable :: word(:)
        integer :: i, j, k, m
        !
        j = 0; k = 0
        do i = 1, len(str)
            if (i .gt. k) then
            if ( str(i:i) .ne. delimiter ) then
                ! determine the next position of delimiter
                k=i+index(str(i:),delimiter)-2
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
                m = index(str(i:),delimiter)
                ! determine the next position of delimiter
                k=i+index(str(i:),delimiter)-2
                ! increase word count
                j=j+1
                if (m.ne.0) then
                    word(j) = str(i:k)
                else
                    word(j) = str(i:)
                endif
            endif
            endif
        enddo
        !
    end function  strsplit

    elemental function strrep(string,place,ins) result(new)
        !
        implicit none
        !
        character(len=*), intent(in)                          :: string,place,ins
        character(len=len(string)+max(0,len(ins)-len(place))) :: new
        integer                                               :: idx
        idx = index(string,place)
        if ( idx == 0 ) then
          new = string
        else
          new = string(1:idx-1)//ins//string(idx+len(place):len(string))
        end if
    end function       strrep 

    ! type conversion

    function     str2int(str) result(num)
        !
        implicit none
        !
        character(*), intent(in) :: str
        integer :: num
        !
        num = nint(str2dbl(str))
        !
    end function str2int

    function     str2dbl(str) result(num)
        !
        implicit none
        !
        character(*), intent(in) :: str
        real(dp) :: num
        !
        read(str,*) num
        !
    end function str2dbl
    
    ! vector

    pure function d_cross_product(a,b) result(c)
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
    end function  d_cross_product

    pure function z_cross_product(a,b) result(c)
        !
        implicit none
        !
        complex(dp), INTENT(IN) :: a(3), b(3)
        complex(dp) :: c(3)
        !
        c(1) = a(2) * b(3) - a(3) * b(2)
        c(2) = a(3) * b(1) - a(1) * b(3)
        c(3) = a(1) * b(2) - a(2) * b(1)
        !
    end function  z_cross_product

    ! math

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

    ! grid functions

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
        real(dp), intent(in), optional :: n3(:)
        real(dp), allocatable :: grid_points(:,:) !> voronoi points (27=3^3)
        integer :: i1,i2,i3,j
        !
        if (present(n3)) then
            !
            allocate(grid_points(3,size(n1)*size(n2)*size(n3)))
            grid_points = 0
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
            !
        else
            !
            allocate(grid_points(2,size(n1)*size(n2)))
            grid_points = 0
            !
            j=0
            do i1 = 1, size(n1)
            do i2 = 1, size(n2)
                j=j+1
                grid_points(1,j)=n1(i1)
                grid_points(2,j)=n2(i2)
            enddo
            enddo
            !
        endif
        !
    end function  dmeshgrid

    pure function imeshgrid(n1,n2,n3) result(grid_points)
        !
        implicit none
        !
        integer, intent(in) :: n1(:)
        integer, intent(in) :: n2(:)
        integer, intent(in), optional :: n3(:)
        real(dp), allocatable :: grid_points(:,:)
        !
        if (present(n3)) then
            grid_points=dmeshgrid(real(n1,dp),real(n2,dp),real(n3,dp))
        else
            grid_points=dmeshgrid(real(n1,dp),real(n2,dp))
        endif
        !
    end function  imeshgrid

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

    pure function regspace(d1,d2,d) result(y)
        !
        implicit none
        !
        real(dp), intent(in) :: d1
        real(dp), intent(in) :: d2
        real(dp), intent(in) :: d
        real(dp), allocatable :: y(:)
        integer :: n
        !
        n=nint((d2-d1)/d)
        !
        allocate(y(n))
        y = d1+d*[0:n:1]
        !
    end function  regspace

    ! matrix properties

    pure function ishermitian(X) result(bool)
        !
        implicit none
        !
        complex(dp), intent(in) :: x(:,:)
        logical :: bool
        !
        if (isequal(X,adjoint(X))) then
            bool = .true.
        else
            bool = .false.
        endif
        !        
    end function  ishermitian

    pure function zadjoint(A) result(adj)
        !
        implicit none
        !
        complex(dp), intent(in) :: A(:,:)
        complex(dp), allocatable :: adj(:,:)
        !
        allocate(adj(size(A,2),size(A,1)))
        adj = transpose(conjg(A))
        !
    end function  zadjoint

    pure function dadjoint(A) result(adj)
        !
        implicit none
        !
        real(dp), intent(in) :: A(:,:)
        real(dp), allocatable :: adj(:,:)
        !
        allocate(adj(size(A,2),size(A,1)))
        adj = transpose(A)
        !
    end function  dadjoint

    pure function dtrace(R) result(tr)
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
    end function  dtrace

    pure function ztrace(R) result(tr)
        !
        implicit none
        !
        complex(dp), intent(in) :: R(:,:)
        complex(dp) :: tr
        integer :: i
        !
        tr = 0
        do i = 1,size(R,2)
            tr = tr + R(i,i)
        enddo
        !
    end function  ztrace

    function      z_eigenspace(A) result(V)
        ! Determine eigenvectors which simultaneously diagonalize n square matrices A(:,:,i)
        ! Note: only possible if A(:,:,i) commute with each other.
        !       only possible if the matrix pencil produces exclusively unique (no nondegenerate) eigenvalues!
        implicit none
        !
        complex(dp), intent(in)  :: A(:,:,:)
        complex(dp), allocatable :: V(:,:)
        complex(dp), allocatable :: Q(:,:)
        complex(dp), allocatable :: D(:)
        integer :: m, n, i
        ! get dimensions
        n = size(A,1)
        m = size(A,3)
        if (m.eq.0) stop 'ERROR [eigenanalysis]: m = 0'
        ! allocate space
        allocate(V(n,n))
        ! use matrix pencil approach to lift the degeneracy of eigenvalues
        allocate(Q(n,n))
        Q = 0
        do i = 1,m
            Q = Q + rand() * A(:,:,i)
        enddo
        ! get eigenvectors
        call am_zgeev(A=Q,VR=V,D=D)
        ! check that matrices have been diagonalized properly and save eigenvalues
        do i = 1, m
            ! get Q
            Q = matmul(inv(V),matmul(A(:,:,i),V))
            ! test Q
            if (.not.isequal(diag(diag(Q)),Q,tiny*m**2)) stop 'ERROR [eigenanalysis]: simultaneous matrix diagonalization failed.'
        enddo
        !
    end function  z_eigenspace

    function      d_eigenspace(A) result(V)
        ! Determine eigenvectors which simultaneously diagonalize n square matrices A(:,:,i)
        ! Note: only possible if A(:,:,i) commute with each other.
        !       only possible if the matrix pencil produces exclusively unique (no nondegenerate) eigenvalues!
        implicit none
        !
        real(dp), intent(in)  :: A(:,:,:)
        complex(dp), allocatable :: V(:,:)
        real(dp), allocatable :: Q(:,:)
        integer :: m, n, i
        ! get dimensions
        n = size(A,1)
        m = size(A,3)
        if (m.eq.0) stop 'ERROR [eigenanalysis]: m = 0'
        ! allocate space
        allocate(V(n,n))
        ! use matrix pencil approach to lift the degeneracy of eigenvalues
        allocate(Q(n,n))
        Q = 0
        do i = 1,m
            Q = Q + rand() * A(:,:,i)
        enddo
        ! get eigenvectors
        call am_dgeev(A=Q,VR=V)
        ! check that matrices have been diagonalized properly and save eigenvalues
        do i = 1, m
            ! get Q
            Q = matmul(inv(V),matmul(A(:,:,i),V))
            ! test Q
            if (.not.isequal(diag(diag(Q)),Q,tiny*m**2)) stop 'ERROR [eigenanalysis]: simultaneous matrix diagonalization failed.'
        enddo
        !
    end function  d_eigenspace

    pure function d_rref(A) result(matrix)
        !
        ! note: algorithm on roseta code is broken.
        ! this has been fixed by:
        ! 1) changed: "(n.lt.pivot)" and "(pivot.gt.n)"
        ! 2) changed: "(i.gt.m)"
        !
        implicit none
        !
        real(dp), intent(in) :: A(:,:)
        real(dp),allocatable :: matrix(:,:)
        real(dp),allocatable :: temp(:)
        integer :: pivot, m, n
        integer :: r, i
        !
        allocate(matrix,source=A)
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
                endif
            enddo
            ! swap rows i and r
            temp = matrix(i,:)
            matrix(i,:) = matrix(r,:)
            matrix(r,:) = temp
            !
            matrix(r,:) = matrix(r,:)/matrix(r,pivot)
            do i = 1, m
                if (i.ne.r) matrix(i,:)=matrix(i,:)-matrix(r,:)*matrix(i,pivot)
            enddo
            pivot = pivot + 1
        enddo
        deallocate(temp)
    end function  d_rref

    pure function z_rref(A) result(matrix)
        !
        ! note: algorithm on roseta code is broken.
        ! this has been fixed by:
        ! 1) changed: "(n.lt.pivot)" and "(pivot.gt.n)"
        ! 2) changed: "(i.gt.m)"
        !
        implicit none
        !
        complex(dp), intent(in) :: A(:,:)
        complex(dp),allocatable :: matrix(:,:)
        complex(dp),allocatable :: temp(:)
        integer :: pivot, m, n
        integer :: r, i
        !
        allocate(matrix,source=A)
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
                endif
            enddo
            ! swap rows i and r
            temp = matrix(i,:)
            matrix(i,:) = matrix(r,:)
            matrix(r,:) = temp
            !
            matrix(r,:) = matrix(r,:)/matrix(r,pivot)
            do i = 1, m
                if (i.ne.r) matrix(i,:)=matrix(i,:)-matrix(r,:)*matrix(i,pivot)
            enddo
            pivot = pivot + 1
        enddo
        deallocate(temp)
    end function  z_rref

    ! trim

    pure function i_trim_null(A) result(B)
        !
        implicit none
        !
        integer, intent(in) :: A(:)
        integer, allocatable :: B(:)
        !
        allocate(B, source=A(selector(A.ne.0)) )
        !
    end function  i_trim_null

    pure function vi_trim_null(A) result(B)
        !
        implicit none
        !
        integer, intent(in) :: A(:,:)
        integer,allocatable :: B(:,:)
        integer,allocatable :: zeros(:)
        logical,allocatable :: mask(:)
        integer,allocatable :: inds(:)
        integer :: n, i
        !
        n = size(A,2)
        allocate(inds, source = [1:n] )
        allocate(zeros(size(A,1)))
        zeros = 0
        !
        allocate(mask(n))
        do i = 1, n
            if (isequal(A(:,i),zeros)) then
                mask(i) = .false.
            else
                mask(i) = .true.
            endif
        enddo
        !
        allocate(B, source = A(:,pack(inds,mask)))
        !
    end function  vi_trim_null

    pure function vd_trim_null(A) result(B)
        !
        implicit none
        !
        real(dp), intent(in) :: A(:,:)
        real(dp), allocatable :: B(:,:)
        real(dp), allocatable :: zeros(:)
        logical , allocatable :: mask(:)
        integer , allocatable :: inds(:)
        integer :: n, i
        !
        n = size(A,2)
        allocate(inds, source = [1:n] )
        allocate(zeros(size(A,1)))
        zeros = 0
        !
        allocate(mask(n))
        do i = 1, n
            if (isequal( A(:,i) ,zeros)) then
                mask(i) = .false.
            else
                mask(i) = .true.
            endif
        enddo
        !
        allocate(B, source = A(:,pack(inds,mask)))
        !
    end function  vd_trim_null

    pure function vz_trim_null(A) result(B)
        !
        implicit none
        !
        complex(dp), intent(in) :: A(:,:)
        complex(dp), allocatable :: B(:,:)
        complex(dp), allocatable :: zeros(:)
        logical , allocatable :: mask(:)
        integer , allocatable :: inds(:)
        integer :: n, i
        !
        n = size(A,2)
        allocate(inds, source = [1:n] )
        allocate(zeros(size(A,1)))
        zeros = 0
        !
        allocate(mask(n))
        do i = 1, n
            if (isequal( A(:,i), zeros)) then
                mask(i) = .false.
            else
                mask(i) = .true.
            endif
        enddo
        !
        allocate(B, source = A(:,pack(inds,mask)))
        !
    end function  vz_trim_null

    pure function vs_trim_null(A) result(B)
        !
        implicit none
        !
        character(*), intent(in) :: A(:)
        character(:),allocatable :: B(:)
        integer :: m,n,lt,i,j
        integer :: maxlen
        integer, allocatable :: selector_array(:)
        !
        m = size(A)
        ! initialize selector
        allocate(selector_array(m))
        selector_array = 0
        ! initialize maxlen to get largest string length
        maxlen = 0
        ! loop over string array
        n = 0
        do j = 1, m
            lt = len_trim(A(j))
            ! check if array is not empty
            if (lt.ne.0) then
                n = n + 1
                selector_array(n) = j
            endif
            ! check whether it is the largest array
            if (lt.gt.maxlen) maxlen = len_trim(A(j))
        enddo
        ! allocate new array
        allocate(character(maxlen+1)::B(n))
        ! create new array
        do i = 1, n
            j = selector_array(i)
            B(i) = trim(A(j))
            B(i) = adjustl(B(i))
        enddo
        !
    end function  vs_trim_null

    pure function  s_trim_null(A) result(B)
        !
        implicit none
        !
        character(*), intent(in) :: A
        character(:),allocatable :: B
        !
        allocate(B, source=trim(A))
        !
    end function   s_trim_null

    ! get vector subspace intersection 

    function      d_subspace_intersection(A,B) result(C)
        !
        implicit none
        !
        real(dp), intent(in) :: A(:,:)
        real(dp), intent(in) :: B(:,:)
        real(dp),allocatable :: C(:,:)
        real(dp),allocatable :: ns(:,:)
        integer :: na,nb,m,n
        !
        m = size(A,1)
        na= size(A,2)
        nb= size(B,2)
        n = na+nb
        if (m.ne.size(B,1)) stop 'ERROR [subspace_intersection]: m /= size(B,1)'
        !
        allocate(ns(m,n))
        ns(1:m,[1:na]   ) = A
        ns(1:m,[1:nb]+na) = B
        ! get null space (coefficients)
        ns = null_svd(ns)
        ! check if there is an intersection
        if (size(ns,2).ne.0) then
            ! apply coefficients to construct interception subspace
            C = orth_svd( matmul(A,ns(1:na,:)) )
        else
            allocate(C(0,0))
        endif
        !
    end function  d_subspace_intersection

    function      z_subspace_intersection(A,B) result(C)
        !
        implicit none
        !
        complex(dp), intent(in) :: A(:,:)
        complex(dp), intent(in) :: B(:,:)
        complex(dp),allocatable :: C(:,:)
        complex(dp),allocatable :: ns(:,:)
        integer :: na,nb,m,n
        !
        m = size(A,1)
        na= size(A,2)
        nb= size(B,2)
        n = na+nb
        if (m.ne.size(B,1)) stop 'ERROR [subspace_intersection]: m /= size(B,1)'
        !
        allocate(ns(m,n))
        ns(1:m,[1:na]   ) = A
        ns(1:m,[1:nb]+na) = B
        ! get null space (coefficients)
        ns = null_svd(ns)
        ! check if there is an intersection
        if (size(ns,2).ne.0) then
            ! apply coefficients to construct interception subspace
            C = orth_svd( matmul(A,ns(1:na,:)) ) 
        else
            allocate(C(0,0))
        endif
        !
    end function  z_subspace_intersection

    ! matrix indexing

    pure function selector(mask) result(inds)
        !
        implicit none
        !
        logical, intent(in) :: mask(:)
        integer,allocatable :: inds(:)
        integer :: n
        !
        n = size(mask)
        !
        allocate(inds, source=pack([1:n],mask))
        !
    end function  selector

    pure function ind2sub(dims,ind) result(sub)
        ! test with: ind2sub(dims=[5,3],ind=sub2ind(dims=[5,3],sub=[1,3]))
        ! Map a scalar index of a flat 1D array to the equivalent
        ! d-dimensional index
        ! Example:
        !     ind      -->        sub
        ! | 1  4  7 |      | 1,1  1,2  1,3 |
        ! | 2  5  8 |  --> | 2,1  2,2  2,3 |
        ! | 3  6  9 |      | 3,1  3,2  3,3 |
        integer, intent(in)  :: ind
        integer, intent(in)  :: dims(:)
        integer, allocatable :: sub(:)
        integer :: i, ndims, ind1, ind2
        !
        ind1=ind
        ndims = size(dims)
        allocate(sub(ndims))
        !
        do i=1,ndims-1
            ind2 = (ind1-1)/dims(i) + 1
            sub(i) = ind1 - dims(i)*(ind2-1)
            ind1 = ind2
        enddo
        sub(ndims) = ind1
    end function  ind2sub

    function      sub2ind(dims,sub) result(ind)
        ! test with: ind2sub(dims=[5,3],ind=sub2ind(dims=[5,3],sub=[1,3]))
        ! Map a d-dimensional index to the scalar index of the equivalent flat array
        ! Example:
        !       sub         -->     ind
        ! | 1,1  1,2  1,3 |     | 1  4  7 |
        ! | 2,1  2,2  2,3 | --> | 2  5  8 |
        ! | 3,1  3,2  3,3 |     | 3  6  9 |
        integer, intent(in) :: dims(:)
        integer, intent(in) :: sub(:)
        integer :: ind
        integer, allocatable :: k(:)
        integer :: i, ndims
        !
        ndims = size(dims)
        allocate(k,source=[1,cumprod(dims(1:(ndims-1)))])
        ind = 1
        do i = 1, ndims
            ind = ind + (sub(i)-1)*k(i);
        enddo
    end function  sub2ind

    ! basic sorting

    pure function d_foursort(x) result(y)
        !
        implicit none
        !
        real(dp), intent(in) :: x(4)
        real(dp) :: y(4)
        !
        y = x
        if (y(1).gt.y(2)) then; y([1,2]) = y([2,1]); endif
        if (y(3).gt.y(4)) then; y([3,4]) = y([4,3]); endif
        if (y(1).gt.y(3)) then; y([1,3]) = y([3,1]); endif
        if (y(2).gt.y(4)) then; y([2,4]) = y([4,2]); endif
        if (y(2).gt.y(3)) then; y([2,3]) = y([3,2]); endif
        !
    end function  d_foursort

    pure function i_foursort(x) result(y)
        !
        implicit none
        !
        integer, intent(in) :: x(4)
        integer :: y(4)
        !
        y = x
        if (y(1).gt.y(2)) then; y([1,2]) = y([2,1]); endif
        if (y(3).gt.y(4)) then; y([3,4]) = y([4,3]); endif
        if (y(1).gt.y(3)) then; y([1,3]) = y([3,1]); endif
        if (y(2).gt.y(4)) then; y([2,4]) = y([4,2]); endif
        if (y(2).gt.y(3)) then; y([2,3]) = y([3,2]); endif
        !
    end function  i_foursort

    pure function d_fourrank(n) result(ind)
        ! borrowed from olle, who took it from stack overflow
        implicit none
        !
        real(dp),intent(in) :: n(4)
        integer :: ind(4)
        integer :: low1,high1,low2,high2,highest,lowest,middle1,middle2
        !    
        if ( n(1) <= n(2) ) then
            low1 = 1
            high1 = 2
        else 
            low1 = 2
            high1 = 1
        endif

        if ( n(3) <= n(4) ) then
            low2 = 3
            high2 = 4
        else
            low2 = 4
            high2 = 3
        endif

        if ( n(low1) <= n(low2) ) then
            lowest = low1
            middle1 = low2
        else
            lowest = low2
            middle1 = low1
        endif

        if ( n(high1) >= n(high2) ) then
            highest = high1
            middle2 = high2
        else
            highest = high2
            middle2 = high1
        endif

        if ( n(middle1) < n(middle2) ) then
            ind=(/lowest,middle1,middle2,highest/)
        else
            ind=(/lowest,middle2,middle1,highest/)
        endif
        !
    end function  d_fourrank

    pure function i_fourrank(n) result(ind)
        ! borrowed from olle, who took it from stack overflow
        implicit none
        !
        integer,intent(in) :: n(4)
        integer :: ind(4)
        integer :: low1,high1,low2,high2,highest,lowest,middle1,middle2
        !    
        if ( n(1) <= n(2) ) then
            low1 = 1
            high1 = 2
        else 
            low1 = 2
            high1 = 1
        endif

        if ( n(3) <= n(4) ) then
            low2 = 3
            high2 = 4
        else
            low2 = 4
            high2 = 3
        endif

        if ( n(low1) <= n(low2) ) then
            lowest = low1
            middle1 = low2
        else
            lowest = low2
            middle1 = low1
        endif

        if ( n(high1) >= n(high2) ) then
            highest = high1
            middle2 = high2
        else
            highest = high2
            middle2 = high1
        endif

        if ( n(middle1) < n(middle2) ) then
            ind=(/lowest,middle1,middle2,highest/)
        else
            ind=(/lowest,middle2,middle1,highest/)
        endif
        !
    end function  i_fourrank

    ! cumulative sum / product

    pure function dcumprod(x) result(y)
        !
        implicit none
        !
        real(dp), intent(in) :: x(:)
        real(dp), allocatable :: y(:)
        integer :: i,n
        !
        n = size(x)
        ! 
        allocate(y(n))
        ! initialize
        y(1) = x(1)
        ! loop
        do i = 2, size(x)
            y(i) = y(i-1)*x(i)
        enddo
        !
    end function  dcumprod

    pure function icumprod(x) result(y)
        !
        implicit none
        !
        integer, intent(in) :: x(:)
        integer, allocatable :: y(:)
        integer :: i,n
        !
        n = size(x)
        ! 
        allocate(y(n))
        ! initialize
        y(1) = x(1)
        ! loop
        do i = 2, size(x)
            y(i) = y(i-1)*x(i)
        enddo
        !
    end function  icumprod

    pure function dcumsum(x) result(y)
        !
        implicit none
        !
        real(dp), intent(in) :: x(:)
        real(dp), allocatable :: y(:)
        integer :: i,n
        !
        n = size(x)
        ! 
        allocate(y(n))
        ! initialize
        y(1) = x(1)
        ! loop
        do i = 2, size(x)
            y(i) = y(i-1) + x(i)
        enddo
        !
    end function  dcumsum

    pure function icumsum(x) result(y)
        !
        implicit none
        !
        integer, intent(in) :: x(:)
        integer, allocatable :: y(:)
        integer :: i,n
        !
        n = size(x)
        ! 
        allocate(y(n))
        ! initialize
        y(1) = x(1)
        ! loop
        do i = 2, size(x)
            y(i) = y(i-1) + x(i)
        enddo
        !
    end function  icumsum

    ! minimum 

    pure function minpos(x) result(id)
        !
        implicit none
        !
        real(dp), intent(in) :: x(:)
        integer  :: id
        integer  :: i
        real(dp) :: temp
        ! loop
        temp = 1.0D30
        do i = 1, size(x)
            if (x(i).lt.temp) then
                temp = x(i)
                id = i
            endif
        enddo
        !
    end function  minpos

    ! matrix generation

    pure function ddiag1(M,j) result(d)
        !
        ! gets jth (sub-,super-)diagonal elements of matrix M
        ! if j is not passed, retunrs diagonal elements (j = 0)
        ! j > 0 super
        ! j < 0 sub
        !
        ! for example, 
        ! [ a_{1,1} a_{2,2} a_{3,3} ... a_{n-2,n-2} a_{n-1,n-1} a_{n,n}   ]  ==>  j =  0   diagonal
        ! [ a_{1,2} a_{2,3} a_{3,4} ... a_{n-2,n-1} a_{n-1,n  } 0         ]  ==>  j =  1   superdiagonal # 1
        ! [ a_{1,3} a_{2,4} a_{3,5} ... a_{n-2,n}   0           0         ]  ==>  j =  2   superdiagonal # 2
        ! 
        ! [ a_{1,1} a_{2,2} a_{3,3} ... a_{n-2,n-2} a_{n-1,n-1} a_{n,n}   ]  ==>  j =  0   diagonal
        ! [ a_{2,1} a_{3,2} a_{4,3} ... a_{n-1,n-2} a_{n  ,n-1} 0         ]  ==>  j = -1   subdiagonal # 1
        ! [ a_{3,1} a_{4,2} a_{5,3} ... a_{n  ,n-2} 0           0         ]  ==>  j = -2   subdiagonal # 2
        !
        implicit none
        !
        real(dp), intent(in)  :: M(:,:)
        integer , optional, intent(in) :: j ! selects a sub- or super- diagonal
        real(dp), allocatable :: d(:)
        integer :: n
        integer :: i, k
        !
        if (present(j)) then
            k = j
        else 
            k = 0
        endif
        !
        ! number of elements in off-diagonal
        n = min(size(M,1),size(M,2))-abs(k)
        allocate(d(n))
        !
        ! construct matrix
        !
        if     (k.eq.0) then
            ! diagonal
            do i = 1, n
                d(i) = M(i,i)
            enddo
        elseif (k.gt.0) then
            ! superdiagonal (upper triangular part)
            do i = 1, n
                d(i) = M(i,i+k)
            enddo
        elseif (k.lt.0) then
            ! subdiagonal (lower triangular part)
            do i = 1, n
                d(i) = M(i-k,i)
            enddo
        endif
        !
    end function  ddiag1

    pure function ddiag2(d,j) result(M)
        !
        ! get diagonal elements of matrix M
        implicit none
        !
        real(dp), intent(in)  :: d(:)
        integer , optional, intent(in) :: j ! selects a sub- or super- diagonal
        real(dp), allocatable :: M(:,:)
        integer :: n
        integer :: i, k
        !
        if (present(j)) then
            k = j
        else 
            k = 0
        endif
        !
        n = size(d) + abs(k)
        !
        if (.not.allocated(M)) then
            allocate(M(n,n))
            M=0
        endif
        !
        ! construct matrix
        !
        if     (k.eq.0) then
            ! diagonal
            do i = 1, n
                M(i,i) = d(i)
            enddo
        elseif (k.gt.0) then
            ! superdiagonal (upper triangular part)
            do i = 1, n-abs(k)
                M(i,i+k) = d(i)
            enddo
        elseif (k.lt.0) then
            ! subdiagonal (lower triangular part)
            do i = 1, n-abs(k)
                M(i-k,i) = d(i)
            enddo
        endif
        !
    end function  ddiag2 

    pure function zdiag1(M,j) result(d)
        !
        ! gets jth (sub-,super-)diagonal elements of matrix M
        ! if j is not passed, retunrs diagonal elements (j = 0)
        ! j > 0 super
        ! j < 0 sub
        !
        ! for example, 
        ! [ a_{1,1} a_{2,2} a_{3,3} ... a_{n-2,n-2} a_{n-1,n-1} a_{n,n}   ]  ==>  j =  0   diagonal
        ! [ a_{1,2} a_{2,3} a_{3,4} ... a_{n-2,n-1} a_{n-1,n  } 0         ]  ==>  j =  1   superdiagonal # 1
        ! [ a_{1,3} a_{2,4} a_{3,5} ... a_{n-2,n}   0           0         ]  ==>  j =  2   superdiagonal # 2
        ! 
        ! [ a_{1,1} a_{2,2} a_{3,3} ... a_{n-2,n-2} a_{n-1,n-1} a_{n,n}   ]  ==>  j =  0   diagonal
        ! [ a_{2,1} a_{3,2} a_{4,3} ... a_{n-1,n-2} a_{n  ,n-1} 0         ]  ==>  j = -1   subdiagonal # 1
        ! [ a_{3,1} a_{4,2} a_{5,3} ... a_{n  ,n-2} 0           0         ]  ==>  j = -2   subdiagonal # 2
        !
        implicit none
        !
        complex(dp), intent(in)  :: M(:,:)
        integer , optional, intent(in) :: j ! selects a sub- or super- diagonal
        complex(dp), allocatable :: d(:)
        integer :: n
        integer :: i, k
        !
        if (present(j)) then
            k = j
        else 
            k = 0
        endif
        !
        ! number of elements in off-diagonal
        n = min(size(M,1),size(M,2))-abs(k)
        allocate(d(n))
        !
        ! construct matrix
        !
        if     (k.eq.0) then
            ! diagonal
            do i = 1, n
                d(i) = M(i,i)
            enddo
        elseif (k.gt.0) then
            ! superdiagonal (upper triangular part)
            do i = 1, n
                d(i) = M(i,i+k)
            enddo
        elseif (k.lt.0) then
            ! subdiagonal (lower triangular part)
            do i = 1, n
                d(i) = M(i-k,i)
            enddo
        endif
        !
    end function  zdiag1

    pure function zdiag2(d,j) result(M)
        !
        ! get diagonal elements of matrix M
        implicit none
        !
        complex(dp), intent(in)  :: d(:)
        integer    , optional, intent(in) :: j ! selects a sub- or super- diagonal
        complex(dp), allocatable :: M(:,:)
        integer :: n
        integer :: i, k
        !
        if (present(j)) then
            k = j
        else 
            k = 0
        endif
        !
        n = size(d) + abs(k)
        !
        if (.not.allocated(M)) then
            allocate(M(n,n))
            M=0
        endif
        !
        ! construct matrix
        !
        if     (k.eq.0) then
            ! diagonal
            do i = 1, n
                M(i,i) = d(i)
            enddo
        elseif (k.gt.0) then
            ! superdiagonal (upper triangular part)
            do i = 1, n-abs(k)
                M(i,i+k) = d(i)
            enddo
        elseif (k.lt.0) then
            ! subdiagonal (lower triangular part)
            do i = 1, n-abs(k)
                M(i-k,i) = d(i)
            enddo
        endif
        !
    end function  zdiag2

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

    pure function zeros(n) result(M)
        !> nxn identity matrix
        implicit none
        !
        integer, intent(in) :: n
        integer, dimension(:,:), allocatable :: M
        !
        allocate(M(n,n))
        M = 0.0_dp
        !
    end function  zeros

    pure function zeros_nxm(n) result(M)
        !> nxn identity matrix
        implicit none
        !
        integer, intent(in) :: n(2)
        integer, allocatable :: M(:,:)
        !
        allocate(M(n(1),n(2)))
        M=0.0_dp
        !
    end function  zeros_nxm

    ! matrix-matrix operations

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

    pure function direct_sum(A,B) result(C)
        ! direct sum of matrices A and B
        ! C = [ A , 0 ]
        !     [ 0 , B ]
        implicit none
        !
        real(dp), intent(in) :: A(:,:)
        real(dp), intent(in) :: B(:,:)
        real(dp), allocatable :: C(:,:)
        integer :: Am, An, Bm, Bn
        !
        Am = size(A,1)
        An = size(A,2)
        Bm = size(B,1)
        Bn = size(B,2)
        !
        allocate(C(Am+Bm,An+Bn))
        !
        C = 0.0_dp
        !
        C( 0+1:An, 0+1:Am) = A
        C(An+1:Bn,Am+1:Bm) = B
        !
    end function  direct_sum

    pure function d_fftshift(x) result(y)
        ! converts array:
        !    -5    -4    -3    -2    -1     0     1     2     3     4     5
        ! to:
        !     1     2     3     4     5    -5    -4    -3    -2    -1     0
        ! by permuting positions. The latter is used by fft routines.
        implicit none
        !
        real(dp), intent(in) :: x(:)
        real(dp),allocatable :: y(:)
        integer :: n
        !
        n = size(x)
        !
        allocate(y(n))
        ! 
        y = cshift(x,floor(n/2.0_dp))
        !
    end function  d_fftshift

    pure function i_fftshift(x) result(y)
        ! converts array:
        !    -5    -4    -3    -2    -1     0     1     2     3     4     5
        ! to:
        !     1     2     3     4     5    -5    -4    -3    -2    -1     0
        ! by permuting positions. The latter is used by fft routines.
        implicit none
        !
        integer, intent(in) :: x(:)
        integer,allocatable :: y(:)
        integer :: n
        !
        n = size(x)
        !
        allocate(y(n))
        ! 
        y = cshift(x,floor(n/2.0_dp))
        !
    end function  i_fftshift

    ! logicals

    pure function isint(x) result(bool)
        !
        implicit none
        !
        real(dp), intent(in) :: x
        logical :: bool
        !
        if (abs(nint(x)-x).lt.tiny) then
            bool = .true.
        else
            bool = .false.
        endif
        !        
    end function  isint

    pure function iszero(x) result(bool)
        !
        implicit none
        !
        real(dp), intent(in) :: x
        logical :: bool
        !
        if (abs(x).lt.tiny) then
            bool = .true.
        else
            bool = .false.
        endif
        !        
    end function  iszero

    subroutine     spy(mask)
        !
        implicit none
        !
        logical, intent(in) :: mask(:,:)
        integer :: i,j,n,m
        !
        m = size(mask,1)
        n = size(mask,2)
        !
        write(unit=*,fmt='(5x,a,a,a)') ,'+', repeat('--',n), '-+'
        !
        do i = 1, m
            write(unit=*,fmt='(5x,a)',advance='no') '|'
            do j = 1, n
                if (mask(i,j)) then
                    write(unit=*,fmt='(a)',advance='no') ' X'
                else
                    write(unit=*,fmt='(a)',advance='no') ' .'
                endif
            enddo
            write(unit=*,fmt='(a)',advance='yes') ' |'
        enddo
        write(unit=*,fmt='(5x,a,a,a)') ,'+', repeat('--',n), '-+'
        !
    end subroutine spy

    ! min/max functions

    pure function smallest_nonzero(a) result(y)
        !
        implicit none
        !
        real(dp), intent(in) :: a(:)
        real(dp) :: y
        integer :: i
        !
        y = maxval(abs(a))
        !
        do i = 1, size(a)
            if (.not.iszero(a(i))) then
                if ( abs(a(i)).lt.y ) then
                    y = a(i)
                endif
            endif
        enddo
        !
    end function  smallest_nonzero

    ! isequal

    pure function r0_isequal(x,y,iopt_prec) result(bool)
        !
        implicit none
        !
        real(dp), intent(in) :: x, y
        logical :: bool
        real(dp), intent(in), optional :: iopt_prec
        real(dp) :: prec
        !
        if ( present(iopt_prec) ) then
            prec = iopt_prec
        else
            prec = tiny
        endif
        !
        !
        if (abs(x-y).lt.prec) then
            bool = .true.
        else
            bool = .false.
        endif
        !        
    end function  r0_isequal

    pure function r1_isequal(x,y,iopt_prec) result(bool)
        !
        implicit none
        !
        real(dp), intent(in) :: x(:), y(:)
        logical :: bool
        real(dp), intent(in), optional :: iopt_prec
        real(dp) :: prec
        integer :: i, m
        !
        if ( present(iopt_prec) ) then
            prec = iopt_prec
        else
            prec = tiny
        endif
        !
        m = size(x,1)
        !
        bool = .true.
        do i = 1, m
            if (abs(x(i)-y(i)).gt.prec) then
                bool = .false.
                return
            endif
        enddo
        !          
    end function  r1_isequal

    pure function r2_isequal(x,y,iopt_prec) result(bool)
        !
        implicit none
        !
        real(dp), intent(in) :: x(:,:), y(:,:)
        logical :: bool
        real(dp), intent(in), optional :: iopt_prec
        real(dp) :: prec
        integer :: i, j, m, n
        !
        if ( present(iopt_prec) ) then
            prec = iopt_prec
        else
            prec = tiny
        endif
        !
        m = size(x,1)
        n = size(x,2)
        !
        bool = .true.
        do i = 1, m
        do j = 1, n
            if (abs(x(i,j)-y(i,j)).gt.prec) then
                bool = .false.
                return
            endif
        enddo
        enddo
        !        
    end function  r2_isequal

    pure function c0_isequal(x,y,iopt_prec) result(bool)
        !
        implicit none
        !
        complex(dp), intent(in) :: x, y
        logical :: bool
        real(dp), intent(in), optional :: iopt_prec
        real(dp) :: prec
        !
        if ( present(iopt_prec) ) then
            prec = iopt_prec
        else
            prec = tiny
        endif
        !
        if     (isequal( real(x), real(y),prec)) then
            if (isequal(aimag(x),aimag(y),prec)) then
                bool = .true.
            else
                bool = .false.
            endif
        else
            bool = .false.
        endif
        !
    end function  c0_isequal

    pure function c1_isequal(x,y,iopt_prec) result(bool)
        !
        implicit none
        !
        complex(dp), intent(in) :: x(:), y(:)
        logical :: bool
        real(dp), intent(in), optional :: iopt_prec
        real(dp) :: prec
        !
        if ( present(iopt_prec) ) then
            prec = iopt_prec
        else
            prec = tiny
        endif
        !
        if     (isequal( real(x), real(y),prec)) then
            if (isequal(aimag(x),aimag(y),prec)) then
                bool = .true.
            else
                bool = .false.
            endif
        else
            bool = .false.
        endif
    end function  c1_isequal

    pure function c2_isequal(x,y,iopt_prec) result(bool)
        !
        implicit none
        !
        complex(dp), intent(in) :: x(:,:), y(:,:)
        logical :: bool
        real(dp), intent(in), optional :: iopt_prec
        real(dp) :: prec
        !
        if ( present(iopt_prec) ) then
            prec = iopt_prec
        else
            prec = tiny
        endif
        !
        if     (isequal( real(x), real(y),prec)) then
            if (isequal(aimag(x),aimag(y),prec)) then
                bool = .true.
            else
                bool = .false.
            endif
        else
            bool = .false.
        endif
        !
    end function  c2_isequal

    pure function i0_isequal(x,y) result(bool)
        !
        implicit none
        !
        integer, intent(in) :: x, y
        logical :: bool
        !
        if (x.eq.y) then
            bool = .true.
        else
            bool = .false.
        endif
        !        
    end function  i0_isequal

    pure function i1_isequal(x,y) result(bool)
        !
        implicit none
        !
        integer, intent(in) :: x(:), y(:)
        logical :: bool
        integer :: i, m
        !
        m = size(x,1)
        !
        bool = .true.
        do i = 1, m
            if (x(i).ne.y(i)) then
                bool = .false.
                return
            endif
        enddo
        !          
    end function  i1_isequal

    pure function i2_isequal(x,y) result(bool)
        !
        implicit none
        !
        integer, intent(in) :: x(:,:), y(:,:)
        logical :: bool
        integer :: i, j, m, n
        !
        m = size(x,1)
        n = size(x,2)
        !
        bool = .true.
        do i = 1, m
        do j = 1, n
            if (x(i,j).ne.y(i,j)) then
                bool = .false.
                return
            endif
        enddo
        enddo
        !        
    end function  i2_isequal

    pure function l0_isequal(x,y) result(bool)
        !
        implicit none
        !
        logical, intent(in) :: x, y
        logical :: bool
        !
        if (x.eq.y) then
            bool = .true.
        else
            bool = .false.
        endif
        !        
    end function  l0_isequal

    pure function l1_isequal(x,y) result(bool)
        !
        implicit none
        !
        logical, intent(in) :: x(:), y(:)
        logical :: bool
        integer :: i, m
        !
        m = size(x,1)
        !
        bool = .true.
        do i = 1, m
            if (x(i).ne.y(i)) then
                bool = .false.
                return
            endif
        enddo
        !          
    end function  l1_isequal

    pure function l2_isequal(x,y) result(bool)
        !
        implicit none
        !
        logical, intent(in) :: x(:,:), y(:,:)
        logical :: bool
        integer :: i, j, m, n
        !
        m = size(x,1)
        n = size(x,2)
        !
        bool = .true.
        do i = 1, m
        do j = 1, n
            if (x(i,j).ne.y(i,j)) then
                bool = .false.
                return
            endif
        enddo
        enddo
        !        
    end function  l2_isequal

    ! issubset

#:for KIND in KINDS_issubset
#:for RANK in RANKS_issubset
    #:if KIND == 'real(dp)' or KIND == 'complex(dp)'
    pure function ${KIND[3]}$${RANK}$_issubset(A,B,iopt_prec) result(bool)
    #:else
    pure function ${KIND[3]}$${RANK}$_issubset(A,B) result(bool)
    #:endif
        !
        ${KIND}$, intent(in) :: A${ranksuffix(RANK)}$
        ${KIND}$, intent(in) :: B${ranksuffix(RANK-1)}$
        logical :: bool
        integer :: i,n
    #:if KIND == 'real(dp)' or KIND == 'complex(dp)'
        real(dp), intent(in), optional :: iopt_prec
        real(dp) :: prec
        !
        if ( present(iopt_prec) ) then
            prec = iopt_prec
        else
            prec = tiny
        endif
    #:endif
        !
        n = size(A,${RANK}$)
        bool = .false.
        do i = 1, n
            #:if KIND == 'real(dp)' or KIND == 'complex(dp)'
            if (isequal(A${ranksuffix_last(RANK,last='i')}$,B,iopt_prec=prec)) then
            #:else
            if (isequal(A${ranksuffix_last(RANK,last='i')}$,B)) then
            #:endif
                bool = .true.
                return 
            endif
        enddo
        !
    end function  ${KIND[3]}$${RANK}$_issubset
#:endfor
#:endfor

    ! reallocate

#:for KIND in KINDS_reallocate
#:for RANK in RANKS_reallocate
    pure function ${KIND[3]}$${RANK}$_reallocate(A) result(B)
        !
        ${KIND}$, intent(in) :: A${ranksuffix(RANK)}$
        ${KIND}$,allocatable :: B${ranksuffix(RANK)}$
        !
        allocate(B,source=A)
        !
    end function  ${KIND[3]}$${RANK}$_reallocate
#:endfor
#:endfor

    ! vector_allocate

#:for KIND in KINDS_vector_allocate
#:for RANK in RANKS_vector_allocate
    #:if KIND == 'character(:)' 
    pure subroutine ${KIND[3]}$${RANK}$_vector_allocate(A,vec,str_length)
    #:else
    pure subroutine ${KIND[3]}$${RANK}$_vector_allocate(A,vec)
    #:endif
        !
        ${KIND}$, allocatable, intent(inout) :: A${ranksuffix(RANK)}$
        integer , intent(in) :: vec(${RANK}$)
        #:if KIND == 'character(:)' 
        integer , intent(in) :: str_length
        #:endif
        !
        #:if KIND == 'character(:)' 
        allocate(character(str_length)::A(${ 'vec(' + '),vec('.join([str(x+1) for x in range(0,RANK)]) + ')' }$))
        #:else
        allocate(A(${ 'vec(' + '),vec('.join([str(x+1) for x in range(0,RANK)]) + ')' }$))
        #:endif
        !
    end subroutine  ${KIND[3]}$${RANK}$_vector_allocate
#:endfor
#:endfor

    ! write_xml_attribute

#:for KIND in KINDS_write_xml_attribute
#:for RANK in RANKS_write_xml_attribute
    subroutine     ${KIND[3]}$${RANK}$_write_xml_attribute(unit,value,attribute)
        !
        integer     , intent(in) :: unit
        ${KIND}$    , intent(in) :: value${ranksuffix(RANK)}$
        character(*), intent(in) :: attribute
        !
        write(unit,'(a,/)') '<'//trim(attribute)//'>'
        #:if   KIND == 'character(*)' 
            write(unit,'(*(x,a))') value
        #:elif KIND == 'real(dp)' 
            write(unit,'(*(x,d))') value
        #:elif KIND == 'complex(dp)' 
            write(unit,'(*(x,SP,d,d,"i"))') value
        #:elif KIND == 'integer'
            write(unit,'(*(x,i))') value
        #:elif KIND == 'logical'
            write(unit,'(*(x,l))') value
        #:endif
        write(unit,'(/)')
        write(unit,'(a,/)') '</'//trim(attribute)//'>'
        !
    end subroutine  ${KIND[3]}$${RANK}$_write_xml_attribute
#:endfor
#:endfor

    ! read_xml_attribute

#:for KIND in KINDS_read_xml_attribute
#:for RANK in RANKS_read_xml_attribute
    subroutine     ${KIND[3]}$${RANK}$_read_xml_attribute(unit,value,iostat,iomsg)
        !
        integer     , intent(in) :: unit
        ${KIND}$    , intent(inout) :: value${ranksuffix(RANK)}$
        integer     , intent(out)   :: iostat
        character(*), intent(inout) :: iomsg
        !
        read(unit=unit,fmt='(/)',iostat=iostat,iomsg=iomsg)
        #:if   KIND == 'character(*)'
            read(unit,'(*(x,a))') value
        #:elif KIND == 'real(dp)' 
            read(unit,'(*(x,d))') value
        #:elif KIND == 'complex(dp)' 
            read(unit,'(*(x,SP,d,d,"i"))') value
        #:elif KIND == 'integer'
            read(unit,'(*(x,i))') value
        #:elif KIND == 'logical'
            read(unit,'(*(x,l))') value
        #:endif
        read(unit=unit,fmt='(/)',iostat=iostat,iomsg=iomsg)
        read(unit=unit,fmt='(/)',iostat=iostat,iomsg=iomsg)
        !
    end subroutine  ${KIND[3]}$${RANK}$_read_xml_attribute
#:endfor
#:endfor

    ! unique_inds

#:for KIND in KINDS_unique
#:for RANK in RANKS_unique
    #:if KIND == 'real(dp)' or KIND == 'complex(dp)'
    pure function ${KIND[3]}$${RANK}$_unique_inds(A,iopt_prec) result(inds)
    #:else
    pure function ${KIND[3]}$${RANK}$_unique_inds(A) result(inds)
    #:endif
        !
        ${KIND}$, intent(in) :: A${ranksuffix(RANK)}$
        integer , allocatable :: inds(:)
        integer :: i,j,k,n
    #:if KIND == 'real(dp)' or KIND == 'complex(dp)'
        real(dp), intent(in), optional :: iopt_prec
        real(dp) :: prec
        !
        if ( present(iopt_prec) ) then
            prec = iopt_prec
        else
            prec = tiny
        endif
    #:endif
        !
        n = size(A,${RANK}$)
        allocate(inds(n))
        inds = 0
        !
        k=0
        try_loop : do i = 1, n
            do j = 1, k
                #:if KIND == 'real(dp)' or KIND == 'complex(dp)'
                    if ( isequal(A${ranksuffix_last(RANK,last='i')}$,A${ranksuffix_last(RANK,last='inds(j)')}$,prec) ) cycle try_loop
                #:else
                    if ( isequal(A${ranksuffix_last(RANK,last='i')}$,A${ranksuffix_last(RANK,last='inds(j)')}$) ) cycle try_loop
                #:endif
            enddo
            k = k + 1
            inds(k) = i
        enddo try_loop
        !
        inds = trim_null(inds)
        !
    end function  ${KIND[3]}$${RANK}$_unique_inds
#:endfor
#:endfor

    ! unique

#:for KIND in KINDS_unique
#:for RANK in RANKS_unique
    #:if KIND == 'real(dp)' or KIND == 'complex(dp)'
    pure function ${KIND[3]}$${RANK}$_unique(A,prec) result(B)
    #:else
    pure function ${KIND[3]}$${RANK}$_unique(A) result(B)
    #:endif
        !
        implicit none
        !
        ${KIND}$, intent(in) :: A${ranksuffix(RANK)}$
        ${KIND}$, allocatable :: B${ranksuffix(RANK)}$
    #:if KIND == 'real(dp)' or KIND == 'complex(dp)'
        real(dp), intent(in), optional :: prec
        if ( present(prec) ) then
            allocate(B, source=A${ranksuffix_last(RANK,last='unique_inds(A,prec)')}$)
        else
            allocate(B, source=A${ranksuffix_last(RANK,last='unique_inds(A,tiny)')}$)
        endif
    #:else
        allocate(B, source=A${ranksuffix_last(RANK,last='unique_inds(A)')}$)
    #:endif
        !
    end function  ${KIND[3]}$${RANK}$_unique
#:endfor
#:endfor

end module am_matlab
