module am_matlab

    use am_constants
    use am_mkl, only : am_zheev

    implicit none
    
    public

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

    interface unique
        module procedure dvunique, dmunique, ivunique, iunique
    end interface ! unique

    interface issubset
        module procedure dm_issubset, im_issubset, dv_issubset, iv_issubset, d_issubset, i_issubset
    end interface ! issubset

    interface meshgrid
        module procedure dmeshgrid, imeshgrid
    end interface ! meshgrid

    interface isequal
        module procedure :: d_isequal, dv_isequal, dm_isequal, z_isequal, zv_isequal, zm_isequal
    end interface ! isequal

    interface trace
        module procedure :: dtrace, ztrace
    end interface ! trace
    
    contains

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

    ! rotation functions

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
            aa(1:3) = axis
            aa(4)   = pi
        else
            phi = acos( (tr-1.0_dp)/2.0_dp )
            axis = [Rcp(3,2)-Rcp(2,3),Rcp(1,3)-Rcp(3,1),Rcp(2,1)-Rcp(1,2)]/(2.0_dp*sin(phi))
            aa(1:3) = axis
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
        !
        ! get dimensions of matrices
        n = 2*l+1
        !
        ! determine similarity transform to convert spherical into tesseral harmonics (complex to real)
        C = tesseral(l=l)
        !
        ! determine rotation around X = real( (B*V)*A*(B*V)' )
        call am_zheev(A=Lx(l),V=V,D=D)
        H = matmul(C,V)
        A = diag(exp(-cmplx_i*euler(1)*D))
        X = matmul(H,matmul(A,adjoint(H)))
        !
        ! determine rotation around Y = real( B*P*B' )
        allocate(Y(n,n))
        B = diag(exp( cmplx_i*euler(2)*D)) ! there should be a minus sign here...
        Y = matmul(C,matmul(B,adjoint(C)))
        !
        ! determine rotation around Z = real( (B*V)*T*(B*V)' )
        call am_zheev(A=Lz(l),V=V,D=D)
        H = matmul(C,V)
        G = diag(exp(-cmplx_i*euler(3)*D))
        Z = matmul(H,matmul(G,adjoint(H)))
        !
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
        !
        ! get dimensions of matrices
        n = nint(2.0_dp*j+1.0_dp)
        !
        ! determine rotation around X = real( (B*V)*A*(B*V)' )
        call am_zheev(A=Jx(j),V=V,D=D)
        A = diag(exp(-cmplx_i*euler(1)*D))
        X = matmul(V,matmul(A,adjoint(V)))
        !
        ! determine rotation around Y = real( B*P*B' )
        allocate(Y(n,n))
        B = diag(exp( cmplx_i*euler(2)*D))
        Y = B
        !
        ! determine rotation around Z = real( (B*V)*T*(B*V)' )
        call am_zheev(A=Jz(j),V=V,D=D)
        G = diag(exp( cmplx_i*euler(3)*D))
        Z = matmul(V,matmul(G,adjoint(V)))
        !
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

    function      rot2irrep(l,R) result(irrep)
        !
        ! generates n-dimensional irreducible representation matrix corresponding to (im)proper rotation R 
        !
        implicit none
        !
        integer , intent(in) :: l
        real(dp), intent(in) :: R(3,3)
        complex(dp), allocatable :: irrep(:,:)
        integer  :: improper_fac
        real(dp) :: d ! det
        real(dp) :: n
        !
        ! convert rotoinversion to pure rotation then back to rotoinversion
        d = R(1,1)*R(2,2)*R(3,3)-R(1,1)*R(2,3)*R(3,2)-R(1,2)*R(2,1)*R(3,3)+R(1,2)*R(2,3)*R(3,1)+R(1,3)*R(2,1)*R(3,2)-R(1,3)*R(2,2)*R(3,1)
        !
        ! check that the rotation is unitary, with determinant equalt to plus or minus one
        if (abs(abs(d)-1.0_dp).gt.tiny) stop 'Rotation is not unitary. det /= +1 or -1'
        !
        ! remove rouding errors
        d = nint(sign(1.0_dp,d))
        !
        ! determine dimension n
        n = (2*l+1)
        !
        ! rotoinversions factor : (-1)**l <= this corresponds to inversion matrix, use d instead of conditional statement
        improper_fac = (d)**l
        !
        ! converts (im)rotation -> proper rotation -> builds S03 rpresentation -> tack on (im)proper factor to recover inversion
        irrep = euler2SO3(l=l, euler=rot2euler(R*d)) * improper_fac
        !
    end function  rot2irrep

    function      proper_th2chi(n,th) result(chi)
        !
        ! Character of proper rotation irrep with dimension n and angle th
        !
        ! T. Wolfram and Ş. Ellialtıoğlu, Applications of Group Theory to Atoms, Molecules,
        ! and Solids, 1 edition (Cambridge University Press, Cambridge, 2014), p 74, Eq. 3.20.
        !
        ! M. Tinkham, Group Theory and Quantum Mechanics (McGraw-Hill, 1964), p 66 Eq 4.6 
        ! and p 100 bottom.
        !
        implicit none
        !
        integer,  intent(in) :: n  ! dimension of irrep
        real(dp), intent(in) :: th ! rotation angle
        real(dp)    :: j           ! azimuthal + spin quantum number
        complex(dp) :: chi
        !
        j = (n-1.0_dp)/2.0_dp
        !
        ! check that j is an integer or half integer number
        if (abs(modulo(j+tiny,0.5_dp)-tiny).gt.tiny) stop 'Quantum number j invalid.'
        !
        ! evaluate the character of the irrep
        if (abs(th).lt.tiny) then
            ! limiting case is obtained by expanding sin( (j+1/2)*x ) / sin( x/2 ) at zero
            chi = (2*j+1)
        else
            ! general case
            chi = sin( (j+0.5_dp)*th ) / sin( th*0.5_dp )
        endif
        !
    end function  proper_th2chi

    function      irrep_character(n,R) result(chi)
        !
        ! Character of n-dimensional orthogonal group O(3) irrep parameterized by 3x3 rotoinversion matrix R
        !       
        ! T. Wolfram and Ş. Ellialtıoğlu, Applications of Group Theory to Atoms, Molecules, 
        ! and Solids, 1 edition (Cambridge University Press, Cambridge, 2014), p 74, Eq. 3.20.
        !
        ! M. Tinkham, Group Theory and Quantum Mechanics (McGraw-Hill, 1964), p 66 Eq 4.6 
        ! and p 100 bottom.
        !
        implicit none
        !
        integer,  intent(in) :: n ! dimension of irrep
        real(dp), intent(in) :: R(3,3) ! rotation angle
        real(dp) :: aa(4)
        real(dp) :: d
        integer :: l
        complex(dp) :: chi
        !
        ! if R is a rotoinversion, get the angle and axis of the rotational part only (without the inversion)
        d = R(1,1)*R(2,2)*R(3,3)-R(1,1)*R(2,3)*R(3,2)-R(1,2)*R(2,1)*R(3,3)+R(1,2)*R(2,3)*R(3,1)+R(1,3)*R(2,1)*R(3,2)-R(1,3)*R(2,2)*R(3,1)
        !
        aa = rot2axis_angle(R)
        !
        l = nint((n-1.0_dp)/2.0_dp)
        !
        chi = proper_th2chi(n=n,th=aa(4)) * sign(1.0_dp,d)**(l)
        !
    end function  irrep_character

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

    ! statistics functions

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
        character(max_argument_length) :: fname
        logical :: fexists
        !
        inquire(file=trim(fname),exist=fexists)
        !
    end function  fexists

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

    ! matrix properties

    pure function adjoint(A)
        !
        implicit none
        !
        complex(dp), intent(in) :: A(:,:)
        complex(dp), allocatable :: adjoint(:,:)
        !
        allocate(adjoint(size(A,2),size(A,1)))
        adjoint = transpose(conjg(A))
        !
    end function  adjoint

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

    ! unique

    pure function dmunique(A,iopt_prec) result(B)
        !> returns unique matrices of A(:,:,i) within numerical precision
        implicit none
        !
        real(dp), intent(in) :: A(:,:,:)
        real(dp) :: wrkspace(size(A,1),size(A,2),size(A,3))
        real(dp), allocatable :: B(:,:,:)
        integer :: i,j,k
        real(dp), intent(in), optional :: iopt_prec
        real(dp) :: prec
        !
        if ( present(iopt_prec) ) then
            prec = iopt_prec
        else
            prec = tiny
        endif
        !
        k=1
        wrkspace(:,:,1) = A(:,:,1)
        !
        try_loop : do i = 1, size(A,3)
            do j = 1, k
                if ( all(abs(A(:,:,i)-wrkspace(:,:,j)).lt.prec) ) cycle try_loop
            enddo
            k = k + 1
            wrkspace(:,:,k) = A(:,:,i)
        enddo try_loop
        !
        allocate(B(size(A,1),size(A,2),k))
        B(:,:,1:k) = wrkspace(:,:,1:k)
        !
    end function  dmunique

    pure function dvunique(A,iopt_prec) result(B)
        !> returns unique columns of double matrix A(:,i) within numerical precision
        implicit none
        !
        real(dp), intent(in) :: A(:,:)
        real(dp) :: wrkspace(size(A,1),size(A,2))
        real(dp), allocatable :: B(:,:)
        integer :: i,j,k
        real(dp), intent(in), optional :: iopt_prec
        real(dp) :: prec
        !
        if ( present(iopt_prec) ) then
            prec = iopt_prec
        else
            prec = tiny
        endif
        !
        k=1
        wrkspace(:,1) = A(:,1)
        !
        try_loop : do i = 1,size(A,2)
            do j = 1,k
                if ( all(abs(A(:,i)-wrkspace(:,j)).lt.prec) ) cycle try_loop
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
    end function  dvunique

    pure function ivunique(A) result(B)
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
    end function  ivunique

    pure function iunique(A) result(B)
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
    end function  iunique

    ! issubset

    pure function dm_issubset(A,B,iopt_prec) result(bool)
        ! returns true if B(:,:) is a subset of A(:,:,i)
        implicit none
        !
        real(dp), intent(in) :: A(:,:,:)
        real(dp), intent(in) :: B(:,:)
        logical :: bool
        integer :: i,n
        real(dp), intent(in), optional :: iopt_prec
        real(dp) :: prec
        !
        if ( present(iopt_prec) ) then
            prec = iopt_prec
        else
            prec = tiny
        endif
        !
        n = size(A,3)
        !
        bool = .false.
        do i = 1, n
            if (all(abs(A(:,:,i)-B(:,:)).lt.prec)) then
                bool = .true.
                return 
            endif
        enddo
        !
    end function  dm_issubset

    pure function im_issubset(A,B) result(bool)
        ! returns true if B(:,:) is a subset of A(:,:,i)
        implicit none
        !
        integer, intent(in) :: A(:,:,:)
        integer, intent(in) :: B(:,:)
        logical :: bool
        integer :: i,n
        !
        n = size(A,3)
        !
        bool = .false.
        do i = 1, n
            if (all(A(:,:,i).eq.B(:,:))) then
                bool = .true.
                return 
            endif
        enddo
        !
    end function  im_issubset

    pure function dv_issubset(A,B,iopt_prec) result(bool)
        ! returns true if B(:,:) is a subset of A(:,:,i)
        implicit none
        !
        real(dp), intent(in) :: A(:,:)
        real(dp), intent(in) :: B(:)
        logical :: bool
        integer :: i,n
        real(dp), intent(in), optional :: iopt_prec
        real(dp) :: prec
        !
        if ( present(iopt_prec) ) then
            prec = iopt_prec
        else
            prec = tiny
        endif
        !
        n = size(A,2)
        !
        bool = .false.
        do i = 1, n
            if (all(abs(A(:,i)-B(:)).lt.prec)) then
                bool = .true.
                return 
            endif
        enddo
        !
    end function  dv_issubset

    pure function iv_issubset(A,B) result(bool)
        ! returns true if B(:,:) is a subset of A(:,:,i)
        implicit none
        !
        integer, intent(in) :: A(:,:)
        integer, intent(in) :: B(:)
        logical :: bool
        integer :: i,n
        !
        n = size(A,2)
        !
        bool = .false.
        do i = 1, n
            if (all(A(:,i).eq.B(:))) then
                bool = .true.
                return 
            endif
        enddo
        !
    end function  iv_issubset

    pure function d_issubset(A,B,iopt_prec) result(bool)
        ! returns true if B is a subset of Ai)
        implicit none
        !
        real(dp), intent(in) :: A(:)
        real(dp), intent(in) :: B
        logical :: bool
        integer :: i,n
        real(dp), intent(in), optional :: iopt_prec
        real(dp) :: prec
        !
        if ( present(iopt_prec) ) then
            prec = iopt_prec
        else
            prec = tiny
        endif
        !
        n = size(A,1)
        !
        bool = .false.
        do i = 1, n
            if (abs(A(i)-B).lt.prec) then
                bool = .true.
                return 
            endif
        enddo
        !
    end function  d_issubset

    pure function i_issubset(A,B) result(bool)
        ! returns true if B is a subset of A(i)
        implicit none
        !
        integer, intent(in) :: A(:)
        integer, intent(in) :: B
        logical :: bool
        integer :: i,n
        !
        n = size(A,1)
        !
        bool = .false.
        do i = 1, n
            if (A(i).eq.B) then
                bool = .true.
                return 
            endif
        enddo
        !
    end function  i_issubset

    ! _isequal

    pure function d_isequal(x,y,iopt_prec) result(bool)
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
    end function  d_isequal

    pure function dv_isequal(x,y,iopt_prec) result(bool)
        !
        implicit none
        !
        real(dp), intent(in) :: x(:), y(:)
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
        if (all(abs(x-y).lt.prec)) then
            bool = .true.
        else
            bool = .false.
        endif
        !        
    end function  dv_isequal

    pure function dm_isequal(x,y,iopt_prec) result(bool)
        !
        implicit none
        !
        real(dp), intent(in) :: x(:,:), y(:,:)
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
        if (all(abs(x-y).lt.prec)) then
            bool = .true.
        else
            bool = .false.
        endif
        !        
    end function  dm_isequal

    pure function z_isequal(x,y,iopt_prec) result(bool)
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
        bool = d_isequal(real(x),real(y),iopt_prec) * d_isequal(aimag(x),aimag(y),iopt_prec)
        !
    end function  z_isequal

    pure function zv_isequal(x,y,iopt_prec) result(bool)
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
        bool = dv_isequal(real(x),real(y),iopt_prec) * dv_isequal(aimag(x),aimag(y),iopt_prec)
        !
    end function  zv_isequal

    pure function zm_isequal(x,y,iopt_prec) result(bool)
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
        bool = dm_isequal(real(x),real(y),iopt_prec) * dm_isequal(aimag(x),aimag(y),iopt_prec)
        !
    end function  zm_isequal

    ! pure function dv_isequal(a) result(bool)
    ! !
    ! ! are all elements of a equal to each other
    ! !
    ! implicit none
    ! !
    ! real(dp), intent(in) :: a(:)
    ! logical :: bool
    ! !
    ! if ( all(abs(a(1)-a).lt.tiny) ) then
    !     bool = .true.
    ! else
    !     bool = .false.
    ! endif
    ! !
    ! end function  dv_isequal
 
    ! pure function dv_distinct(a) result(bool)
    !     !
    !     ! are all elements of a different from each other?
    !     !
    !     implicit none
    !     !
    !     real(dp), intent(in) :: a(:)
    !     logical :: bool
    !     integer :: i, j, n
    !     !
    !     n = size(a)
    !     !
    !     if (n.eq.1) then
    !         bool = .false.
    !         return
    !     endif    
    !     !
    !     bool = .true.
    !     do i = 1, n
    !     do j = 1, i
    !         if (i.ne.j) then
    !         if (abs(a(i)-a(j)).lt.tiny) then
    !             bool = .false.
    !             return
    !         endif
    !         endif
    !     enddo
    !     enddo        
    ! end function  dv_distinct

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

end module am_matlab
