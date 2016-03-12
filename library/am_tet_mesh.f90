    module am_tet_mesh
    !
    use am_constants
    !
    implicit none
    !
    public :: tessellate_mesh
    !
    private
    !
    contains
    !
    subroutine tessellate_mesh(kpt,corner,volume)
    ! wrapper for tesselation. 
    ! 
    implicit none
    !
    real(dp), intent(in) :: kpt(:,:) ! kpt(1:3,nkpts)
    integer , intent(out), allocatable :: corner(:,:) ! corner(1:4,ntets)
    real(dp), intent(out), allocatable :: volume(:) ! tetrahedron volume
    ! dimensions
    integer :: nkpts
    integer :: ndim
    ! work space dimensions
    integer, parameter :: bf_max = 8000
    integer, parameter :: fc_max = 40000
    integer,allocatable,dimension(:) :: ht
    ! other stuff
    integer :: bf_num
    integer,dimension(1:3,1:bf_max) :: bf
    integer :: face_num
    integer,dimension(1:7,1:fc_max) :: fc
    integer :: fc_num
    integer :: ht_num
    integer :: ierror
    integer :: ntets
    integer :: tetra_num2
    integer, allocatable :: vm(:)
    ! loop variable
    integer :: i
    !
    integer :: verbosity
    !
    verbosity = 1
    !
    if (verbosity .ge. 1) write(*,*)
    if (verbosity .ge. 1) write(*,'(a)' ) 'Tessellating tetrahedra'
    !
    ierror = 0
    !
    ndim   = size(kpt,1)
    if (verbosity .ge. 1) write(*,'(a5,a20,i6)' ) ' ... ', 'spatial dimensions', ndim
    nkpts  = size(kpt,2)
    if (verbosity .ge. 1) write(*,'(a5,a20,i6)' ) ' ... ', 'kpoints', nkpts
    ht_num = (3 * nkpts) / 2  ! size of workspace
    if (verbosity .ge. 1) write(*,'(a5,a20,i6)' ) ' ... ', 'hash table size', ht_num
    !
    allocate(ht(ht_num))
    allocate(vm(1:nkpts))
    !
    do i = 1,nkpts
        vm(i) = i
    enddo
    !
    call dtris3(nkpts,ht_num,bf_max,fc_max,kpt,vm,bf_num,fc_num,face_num,ntets,bf,fc,ht,ierror)
    !
    if(ierror /= 0) write(*,'(a5,a)' ) ' ERR ', 'Fatal error encounted!'
    !
    ! if (verbosity .ge. 1) write(*,'(a5,a20,i6)' ) ' ... ', 'max number of bf', bf_max
    if (verbosity .ge. 1) write(*,'(a5,a20,i6)' ) ' ... ', 'boundary faces', bf_num
    ! if (verbosity .ge. 1) write(*,'(a5,a20,i6)' ) ' ... ', 'max number of fc', fc_max
    if (verbosity .ge. 1) write(*,'(a5,a20,i6)' ) ' ... ', 'face records', fc_num
    !
    if (verbosity .ge. 1) write(*,'(a5,a20,i6)' ) ' ... ', 'tetrahedra', ntets
    !
    allocate(corner(1:4,1:ntets))
    !
    call tetlst(fc_max,fc_num,vm,fc,ntets,tetra_num2,corner)
    !
    ! compute volume (IN PROGRESS)
    !
    allocate(volume(ntets))

    do i = 1, ntets
         volume(i) = tetrahedron_volume( kpt(1:3,corner(1:4,i)) )
    end do
    !
    if (verbosity .ge. 1) write(*,'(a5,a20,f6.2)' ) ' ... ', 'maximum volume', maxval(volume)
    if (verbosity .ge. 1) write(*,'(a5,a20,f6.2)' ) ' ... ', 'minimum volume', minval(volume)
    if (verbosity .ge. 1) write(*,'(a5,a20,f6.2)' ) ' ... ', 'average volume', sum(volume)/(1.0_dp*ntets)
    if (verbosity .ge. 1) write(*,'(a5,a20,f6.2)' ) ' ... ', 'total volume', sum(volume)
    !
    !

    end subroutine tessellate_mesh
    !
    ! Stuff below is borrowed from John Burkardt's TABLE_TET_MESH.
    !
    subroutine availf(hdavfc,fc_num,fc_max,fc,ind,ierr)

    !*****************************************************************************80
    !
    !! AVAILF returns the index of the next available record in the FC array.
    !
    !  Discussion: 
    !
    !    This routine returns the index of the next available record in the
    !    FC array,either HDAVFC or FC_NUM+1.
    !
    !  Licensing:
    !
    !    This code is distributed under the GNU LGPL license. 
    !
    !  Modified:
    !
    !    07 September 2005
    !
    !  Author:
    !
    !    Original FORTRAN77 version by Barry Joe.
    !    FORTRAN90 version by John Burkardt.
    !
    !  Reference:
    !
    !    Barry Joe,
    !    GEOMPACK - a software package for the generation of meshes
    !      using geometric algorithms,
    !    Advances in Engineering Software,
    !    Volume 13,pages 325-331,1991.
    !
    !  Parameters:
    !
    !    Input/output,integer :: HDAVFC,head pointer of available
    !    records in FC.
    !
    !    Input/output,integer :: FC_NUM,current number of records 
    !    used in FC.
    !
    !    Input,integer :: FC_MAX,the maximum number of records 
    !    available in FC.
    !
    !    Input,integer :: FC(1:7,1:*),array of face records; see 
    !    routine DTRIS3.
    !
    !    Output,integer :: IND,the index of available record(if FC 
    !    not full).
    !
    !    Output,integer :: IERR,error flag,which is zero unless an 
    !    error occurred.
    !
      implicit none

      integer :: fc(7,*)
      integer :: fc_num
      integer :: hdavfc
      integer :: ierr
      integer :: ind
      integer :: fc_max

      ierr = 0

      if(hdavfc /= 0) then
        ind = hdavfc
        hdavfc = -fc(1,hdavfc)
      else if(fc_max <= fc_num) then
        ierr = 11
        write(*,'(a)') ' '
        write(*,'(a)') 'AVAILF - Fatal error!'
        write(*,'(a)') '  Memory requirements for array FC exceed the'
        write(*,'(a,i12)') '  current limit of FC_MAX = ',fc_max
      else
        fc_num = fc_num + 1
        ind = fc_num
      end if

      return
    end subroutine
    subroutine baryth(a,b,c,d,e,alpha,degen)

    !*****************************************************************************80
    !
    !! BARYTH computes barycentric coordinates of a point in 3D.
    !
    !  Discussion: 
    !
    !    This routine computes the barycentric coordinates of a 3D point with
    !    respect to the four vertices of a tetrahedron.
    !
    !  Licensing:
    !
    !    This code is distributed under the GNU LGPL license. 
    !
    !  Modified:
    !
    !    24 January 2009
    !
    !  Author:
    !
    !    Original FORTRAN77 version by Barry Joe.
    !    FORTRAN90 version by John Burkardt.
    !
    !  Reference:
    !
    !    Barry Joe,
    !    GEOMPACK - a software package for the generation of meshes
    !    using geometric algorithms,
    !    Advances in Engineering Software,
    !    Volume 13,pages 325-331,1991.
    !
    !  Parameters:
    !
    !    Input,real(dp) :: A(1:3),B(1:3),C(1:3),D(1:3),4 vertices
    !    of tetrahedron.
    !
    !    Input,real(dp) :: E(1:3),fifth point for which 
    !    barycentric coordinates found
    !
    !    Output,real(dp) :: ALPHA(1:4),the scaled barycentric coordinates
    !   (if DEGEN = .FALSE.) such that 
    !      E =(ALPHA(1)*A + ALPHA(2)*B + ALPHA(3)*C +ALPHA(4)*D)/DET 
    !    where DET = 6 *(volume of tetra ABCD);  an ALPHA(I) may be set to 0 
    !    after tolerance test to indicate that E is coplanar with a face,so 
    !    sum of ALPHA(I)/DET may not be 1; if the actual barycentric
    !    coordinates rather than just their signs are needed,
    !    modify this routine to divide ALPHA(I) by DET.
    !
    !    Output,logical :: DEGEN,TRUE iff A,B,C,D are coplanar.
    !
      implicit none

      real(dp) :: a(3)
      real(dp) :: alpha(4)
      real(dp) :: amax
      real(dp) :: b(3)
      real(dp) :: bmax
      real(dp) :: c(3)
      real(dp) :: cmax
      real(dp) :: cp1
      real(dp) :: cp2
      real(dp) :: cp3
      real(dp) :: d(3)
      real(dp) :: da(3)
      real(dp) :: db(3)
      real(dp) :: dc(3)
      real(dp) :: de(3)
      logical :: degen
      real(dp) :: det
      real(dp) :: dmax
      real(dp) :: e(3)
      real(dp) :: ea(3)
      real(dp) :: eb(3)
      real(dp) :: ec(3)
      real(dp) :: emax
      real(dp) :: tol

      tol = 100.0D+00 * epsilon(tol)
      degen = .false.

      da(1:3) = a(1:3) - d(1:3)
      db(1:3) = b(1:3) - d(1:3)
      dc(1:3) = c(1:3) - d(1:3)

      amax = max(abs(a(1)),abs(a(2)),abs(a(3)))
      bmax = max(abs(b(1)),abs(b(2)),abs(b(3)))
      cmax = max(abs(c(1)),abs(c(2)),abs(c(3)))
      dmax = max(abs(d(1)),abs(d(2)),abs(d(3)))

      cp1 = db(2) * dc(3) - db(3) * dc(2)
      cp2 = db(3) * dc(1) - db(1) * dc(3)
      cp3 = db(1) * dc(2) - db(2) * dc(1)
      det = da(1) * cp1 + da(2) * cp2 + da(3) * cp3

      if(abs(det) <= 0.01D+00 * tol * max(amax,bmax,cmax,dmax)) then
        degen = .true. 
        return
      end if

      de(1:3) = e(1:3) - d(1:3)
      ea(1:3) = a(1:3) - e(1:3)
      eb(1:3) = b(1:3) - e(1:3)
      ec(1:3) = c(1:3) - e(1:3)

      alpha(1) = de(1) * cp1 + de(2) * cp2 + de(3) * cp3

      cp1 = da(2) * de(3) - da(3) * de(2)
      cp2 = da(3) * de(1) - da(1) * de(3)
      cp3 = da(1) * de(2) - da(2) * de(1)

      alpha(2) = dc(1) * cp1 + dc(2) * cp2 + dc(3) * cp3
      alpha(3) = db(1) * cp1 + db(2) * cp2 + db(3) * cp3

      alpha(4) = ea(1) *(eb(2) * ec(3) - eb(3) * ec(2)) &
               + ea(2) *(eb(3) * ec(1) - eb(1) * ec(3)) &
               + ea(3) *(eb(1) * ec(2) - eb(2) * ec(1))

      if(det < 0.0D+00) then
        alpha(1) = -alpha(1)
        alpha(2) = -alpha(2)
        alpha(4) = -alpha(4)
      else
        alpha(3) = -alpha(3)
      end if

      emax = max(abs(e(1)),abs(e(2)),abs(e(3)))

      if(abs(alpha(1)) <= tol * max(bmax,cmax,dmax,emax)) then
        alpha(1) = 0.0D+00
      end if

      if(abs(alpha(2)) <= tol * max(amax,cmax,dmax,emax)) then
        alpha(2) = 0.0D+00
      end if

      if(abs(alpha(3)) <= tol * max(amax,bmax,dmax,emax)) then
        alpha(3) = 0.0D+00
      end if

      if(abs(alpha(4)) <= tol * max(amax,bmax,cmax,emax)) then
        alpha(4) = 0.0D+00
      end if

      return
    end subroutine
    subroutine ccsph(intest,a,b,c,d,e,center,radsq,inv)

    !*****************************************************************************80
    !
    !! CCSPH finds the circumsphere through the vertices of a tetrahedron.
    !
    !  Discussion: 
    !
    !    This routine finds the center and the square of the radius of 
    !    the circumsphere through four vertices of a tetrahedron,and 
    !    possibly determines whether a fifth 3D point is inside the sphere.
    !
    !  Licensing:
    !
    !    This code is distributed under the GNU LGPL license. 
    !
    !  Modified:
    !
    !    24 January 2009
    !
    !  Author:
    !
    !    Original FORTRAN77 version by Barry Joe.
    !    FORTRAN90 version by John Burkardt.
    !
    !  Reference:
    !
    !    Barry Joe,
    !    GEOMPACK - a software package for the generation of meshes
    !    using geometric algorithms,
    !    Advances inv Engineering Software,
    !    Volume 13,pages 325-331,1991.
    !
    !  Parameters:
    !
    !    Input,logical :: INTEST,is TRUE,if and only if the test for fifth point 
    !    inv sphere is to be made.
    !
    !    Input,real(dp) :: A(1:3),B(1:3),C(1:3),D(1:3),vertices 
    !    of tetrahedron.
    !
    !    Input,real(dp) :: E(1:3),a fifth point; referenced if 
    !    and only if INTEST is TRUE.
    !
    !    Output,real(dp) :: CENTER(1:3),center of sphere; undefined 
    !    if A,B,C,D coplanar.
    !
    !    Output,real(dp) :: RADSQ,the square of radius of sphere; 
    !    -1 if A,B,C,D coplanar.
    !
    !    Output,integer :: IN,contains following value if INTEST 
    !    is .TRUE.:
    !     2 if A,B,C,D coplanar; 
    !     1 if E inside sphere;
    !     0 if E on sphere; 
    !    -1 if E outside sphere
    !
      implicit none

      real(dp) :: a(3)
      real(dp) :: b(3)
      real(dp) :: c(3)
      real(dp) :: center(3)
      real(dp) :: cmax
      real(dp) :: cp1
      real(dp) :: cp2
      real(dp) :: cp3
      real(dp) :: d(3)
      real(dp) :: da(3)
      real(dp) :: db(3)
      real(dp) :: dc(3)
      real(dp) :: det
      real(dp) :: dsq
      real(dp) :: e(3)
      integer :: inv
      logical :: intest
      real(dp) :: radsq
      real(dp) :: rhs(3)
      real(dp) :: tol

      tol = 100.0D+00 * epsilon(tol)

      da(1:3) = a(1:3) - d(1:3)
      db(1:3) = b(1:3) - d(1:3)
      dc(1:3) = c(1:3) - d(1:3)

      rhs(1) = 0.5D+00 * sum(da(1:3)**2)
      rhs(2) = 0.5D+00 * sum(db(1:3)**2)
      rhs(3) = 0.5D+00 * sum(dc(1:3)**2)

      cmax = max(&
        abs(a(1)),abs(a(2)),abs(a(3)),&
        abs(b(1)),abs(b(2)),abs(b(3)),&
        abs(c(1)),abs(c(2)),abs(c(3)),&
        abs(d(1)),abs(d(2)),abs(d(3)))

      cp1 = db(2) * dc(3) - dc(2) * db(3)
      cp2 = dc(2) * da(3) - da(2) * dc(3)
      cp3 = da(2) * db(3) - db(2) * da(3)

      det = da(1) * cp1 + db(1) * cp2 + dc(1) * cp3

      if(abs(det) <= 0.01D+00 * tol * cmax) then
        radsq = -1.0D+00
        inv = 2
        return
      end if

      center(1) =(rhs(1) * cp1 + rhs(2) * cp2 + rhs(3) * cp3) / det

      cp1 = db(1) * rhs(3) - dc(1) * rhs(2)
      cp2 = dc(1) * rhs(1) - da(1) * rhs(3)
      cp3 = da(1) * rhs(2) - db(1) * rhs(1)

      center(2) = (da(3) * cp1 + db(3) * cp2 + dc(3) * cp3) / det
      center(3) = -(da(2) * cp1 + db(2) * cp2 + dc(2) * cp3) / det

      radsq = sum(center(1:3)**2)

      center(1:3) = center(1:3) + d(1:3)

      if(intest) then

        dsq = sum((e(1:3) - center(1:3))**2)

        if((1.0D+00 + tol) * radsq < dsq) then
          inv = -1
        else if(dsq <(1.0D+00 - tol) * radsq) then
          inv = 1
        else
          inv = 0
        end if

      end if

      return
    end subroutine
    subroutine dhpsrt(k,n,lda,a,mapv)

    !*****************************************************************************80
    !
    !! DHPSRT sorts a list of double precision points in KD.
    !
    !  Discussion: 
    !
    !    This routine uses heapsort to obtain the permutation of N K-dimensional
    !    double precision points so that the points are in lexicographic
    !    increasing order.
    !
    !  Licensing:
    !
    !    This code is distributed under the GNU LGPL license. 
    !
    !  Modified:
    !
    !    31 August 2005
    !
    !  Author:
    !
    !    Original FORTRAN77 version by Barry Joe.
    !    FORTRAN90 version by John Burkardt.
    !
    !  Reference:
    !
    !    Barry Joe,
    !    GEOMPACK - a software package for the generation of meshes
    !      using geometric algorithms,
    !    Advances in Engineering Software,
    !    Volume 13,pages 325-331,1991.
    !
    !  Parameters:
    !
    !    Input,integer :: K,the dimension of points.
    !
    !    Input,integer :: N,the number of points.
    !
    !    Input,integer :: LDA,the leading dimension of array A in 
    !    calling routine; K <= LDA.
    !
    !    Input,real(dp) :: A(1:K,1:*),array of points.
    !
    !    Input/output,integer :: MAP(1:N).  On input,the points of A 
    !    with indices MAP(1),MAP(2),...,MAP(N) are to be sorted.  On output,
    !    the elements are permuted so that 
    !      A(*,MAP(1)) <= A(*,MAP(2)) <= ... <= A(*,MAP(N)).
    !
      implicit none

      integer :: lda
      integer :: n

      real(dp) :: a(lda,*)
      integer :: i
      integer :: k
      integer :: mapv(n)
      integer :: t

      do i = n/2,1,-1
        call dsftdw(i,n,k,lda,a,mapv)
      end do

      do i = n,2,-1
        t = mapv(1)
        mapv(1) = mapv(i)
        mapv(i) = t
        call dsftdw(1,i-1,k,lda,a,mapv)
      end do

      return
    end subroutine
    function dless(k,p,q)

    !*****************************************************************************80
    !
    !! DLESS determines the lexicographically lesser of two double precision values.
    !
    !  Discussion: 
    !
    !    This routine determines whether P is lexicographically less than Q in
    !    floating point arithmetic.
    !
    !  Licensing:
    !
    !    This code is distributed under the GNU LGPL license. 
    !
    !  Modified:
    !
    !    24 January 2009
    !
    !  Author:
    !
    !    Original FORTRAN77 version by Barry Joe.
    !    FORTRAN90 version by John Burkardt.
    !
    !  Reference:
    !
    !    Barry Joe,
    !    GEOMPACK - a software package for the generation of meshes
    !    using geometric algorithms,
    !    Advances in Engineering Software,
    !    Volume 13,pages 325-331,1991.
    !
    !  Parameters:
    !
    !    Input,integer :: K,dimension of points.
    !
    !    Input,real(dp) :: P(1:K),Q(1:K),two points.
    !
    !    Output,logical :: DLESS,TRUE if P < Q,FALSE otherwise.
    !
      implicit none

      real(dp) :: cmax
      logical :: dless
      integer :: i
      integer :: k
      real(dp) :: p(k)
      real(dp) :: q(k)
      real(dp) :: tol

      tol = 100.0D+00 * epsilon(tol)

      do i = 1,k

        cmax = max(abs(p(i)),abs(q(i)))

        if(abs(p(i) - q(i)) <= tol * cmax .or. cmax <= tol) then
          cycle
        end if

         if(p(i) < q(i)) then
           dless = .true.
         else
           dless = .false.
         end if

         return

      end do

      dless = .false.

      return
    end function
    subroutine dsftdw(l,u,k,lda,a,mapv)

    !*****************************************************************************80
    !
    !! DSFTDW does one step of the heap sort algorithm for double precision data.
    !
    !  Discussion: 
    !
    !    This routine sifts A(*,MAP(L)) down a heap of size U.
    !
    !  Licensing:
    !
    !    This code is distributed under the GNU LGPL license. 
    !
    !  Modified:
    !
    !    24 January 2009
    !
    !  Author:
    !
    !    Original FORTRAN77 version by Barry Joe.
    !    FORTRAN90 version by John Burkardt.
    !
    !  Reference:
    !
    !    Barry Joe,
    !    GEOMPACK - a software package for the generation of meshes
    !    using geometric algorithms,
    !    Advances in Engineering Software,
    !    Volume 13,pages 325-331,1991.
    !
    !  Parameters:
    !
    !    Input,integer :: L,U,the lower and upper index of
    !    part of heap.
    !
    !    Input,integer :: K,the dimension of points.
    !
    !    Input,integer :: LDA,the leading dimension of array A in 
    !    calling routine.
    !
    !    Input,real(dp) :: A(1:K,1:*),see routine DHPSRT.
    !
    !    Input/output,integer :: MAP(1:*),see routine DHPSRT.
    !
      implicit none

      integer :: lda

      integer :: k
      integer :: l
      integer :: mapv(*)
      integer :: u

      real(dp) :: a(lda,*)
      integer :: i
      integer :: j
      integer :: t

      i = l
      j = 2 * i
      t = mapv(i)

      do

        if(u < j) then
          exit
        end if

        if(j < u) then
          if(dless(k,a(1,mapv(j)),a(1,mapv(j+1)))) then
            j = j + 1
          end if
        end if

        if(dless(k,a(1,mapv(j)),a(1,t))) then
          exit
        end if

        mapv(i) = mapv(j)
        i = j
        j = 2 * i

      end do

      mapv(i) = t

      return
    end subroutine
    subroutine dtris3(npt,sizht,bf_max,fc_max,vcl,vm,bf_num,fc_num,&
      nface,ntetra,bf,fc,ht,ierr)

    !*****************************************************************************80
    !
    !! DTRIS3 constructs a Delaunay triangulation of vertices in 3D.
    !
    !  Discussion:
    !
    !    This routine constructs a Delaunay triangulation of 3D vertices using
    !    an incremental approach and local transformations.  Vertices are
    !    first sorted in lexicographically increasing(x,y,z) order,and
    !    then are inserted one at a time from outside the convex hull.
    !
    !  Licensing:
    !
    !    This code is distributed under the GNU LGPL license. 
    !
    !  Modified:
    !
    !    02 September 2005
    !
    !  Author:
    !
    !    Original FORTRAN77 version by Barry Joe.
    !    FORTRAN90 version by John Burkardt.
    !
    !  Reference:
    !
    !    Barry Joe,
    !    GEOMPACK - a software package for the generation of meshes
    !    using geometric algorithms,
    !    Advances in Engineering Software,
    !    Volume 13,pages 325-331,1991.
    !
    !  Parameters:
    !
    !    Input,integer :: NPT,the number of 3D vertices.
    !
    !    Input,integer :: SIZHT,the size of the hash table HT; a good 
    !    choice is a prime number which is about 1/8 * NFACE(or 3/2 * NPT for 
    !    random points from the uniform distribution).
    !
    !    Input,integer :: BF_MAX,the maximum size available for BF 
    !    array.  This needs to be at least as big as the number of boundary faces.
    !
    !    Input,integer :: FC_MAX,the maximum size available for 
    !    FC array.  This needs to be at least as big as the number of faces.
    !
    !    Input,real(dp) :: VCL(1:3,1:NPT),the vertex coordinates.
    !    In the general case,VCL may contain the coordinates for more
    !    than NPT vertices,and the VM array is used to select them.
    !
    !    Input/output,integer :: VM(1:NPT),the vertices of VCL to be 
    !    triangulated.  On output,these indices are permuted,so that 
    !    VCL(*,VM(1)),... ,VCL(*,VM(NPT)) are in lexicographic increasing order,
    !    with possible slight reordering so first 4 vertices are
    !    non-coplanar.  Typically,the input value of VM might be 1 through
    !    NPT.
    !
    !    Output,integer :: BF_NUM,the number of positions used in BF.
    !    BF_NUM <= BF_MAX.
    !
    !    Output,integer :: FC_NUM,the number of positions used in FC.
    !    FC_NUM <= FC_MAX.
    !
    !    Output,integer :: NFACE,the number of faces in triangulation; 
    !    NFACE <= FC_NUM.
    !
    !    Output,integer :: NTETRA,the number of tetrahedra in 
    !    the triangulation.
    !
    !    Output,integer :: BF(1:3,1:BF_NUM),boundary face records 
    !    containing pointers(indices) to FC; if FC(5,I) = -J < 0 and 
    !    FC(1:3,I) = ABC,then BF(1,J) points to other boundary face with edge BC,
    !    BF(2,J) points to other boundary face with edge AC,and
    !    BF(3,J) points to other boundary face with edge AB;
    !    if BF(1,J) <= 0,record is not used and is in avail list.
    !
    !    Output,integer :: FC(1:7,1:FC_NUM),face records which are in 
    !    linked lists in hash table with direct chaining. Fields are:
    !    FC(1:3,*) - A,B,C with 1<=A<B<C<=NPT; indices in VM of 3
    !    vertices of face; if A <= 0,record is not used(it is
    !    in linked list of available records with indices <= FC_NUM);
    !    internal use: if B <= 0,face in queue,not in triangulation.
    !    FC(4:5,*) - D,E; indices in VM of 4th vertex of 1 or 2
    !    tetrahedra with face ABC; if ABC is boundary face
    !    then E < 0 and |E| is an index of BF array
    !    FC(6,*) - HTLINK; pointer(index in FC) of next element
    !    in linked list(or NULL = 0)
    !    FC(7,*) - used internally for QLINK(link for queues or
    !    stacks); pointer(index in FC) of next face in queue/
    !    stack(or NULL = 0); QLINK = -1 indicates face is not
    !    in any queue/stack,and is output value(for records
    !    not in avail list),except:
    !    FC(7,1:2) - HDAVBF,HDAVFC : head pointers of avail list in BF,FC.
    !
    !    Output,integer :: HT(0:SIZHT-1),a hash table using direct
    !    chaining; entries are head pointers of linked lists(indices of FC array)
    !    containing the faces and tetrahedra of the triangulation.
    !
    !    Output,integer :: IERR,error flag,which is zero unless an 
    !    error occurred.
    !
      implicit none

      integer :: bf_max
      integer :: fc_max
      integer :: npt
      integer :: sizht

      integer :: a
      integer :: b
      integer :: back
      integer :: bf(3,bf_max)
      integer :: bfi
      real(dp) :: ctr(3)
      integer :: e
      integer :: fc(7,fc_max)
      integer :: fc_num
      integer :: front
      integer :: hdavbf
      integer :: hdavfc
      integer :: ht(0:sizht-1)
      integer :: i
      integer :: i3
      integer :: i4
      integer :: ierr
      integer :: ip
      integer :: j
      integer :: k
      integer,parameter :: msglvl = 0
      integer :: bf_num
      integer :: nface
      integer :: ntetra
      integer :: op
      ! integer :: opside
      integer :: ptr
      integer :: topnv
      integer :: top
      integer :: va
      integer :: vb
      integer :: vc
      real(dp) :: vcl(3,*)
      integer :: vi
      integer :: vm(npt)

      ierr = 0
    !
    !  Permute elements of VM so that vertices are in lexicographic order.
    !
      call dhpsrt(3,npt,3,vcl,vm)
    !
    !  Reorder points so that first four points are in general position.
    !
      call frstet(.true.,npt,vcl,vm,i3,i4,ierr)

      if(ierr /= 0) then
        write(*,'(a)') ' '
        write(*,'(a)') 'DTRIS3 - Error!'
        write(*,'(a,i8)') '  FRSTET returned IERR = ',ierr
        return
      end if
    !
    !  Initialize data structures.
    !
      do i = 1,3
        ctr(i) = sum(vcl(i,vm(1:4))) / 4.0D+00
      end do

      ht(0:sizht-1) = 0
      hdavbf = 0
      hdavfc = 0
      bf_num = 4
      fc_num = 4
      ntetra = 1

      call htins(1,1,2,3,4,-1,npt,sizht,fc,ht)
      call htins(2,1,2,4,3,-2,npt,sizht,fc,ht)
      call htins(3,1,3,4,2,-3,npt,sizht,fc,ht)
      call htins(4,2,3,4,1,-4,npt,sizht,fc,ht)

      bf(1:3,1) =(/ 4,3,2 /)
      bf(1:3,2) =(/ 4,3,1 /)
      bf(1:3,3) =(/ 4,2,1 /)
      bf(1:3,4) =(/ 3,2,1 /)

      if(msglvl == 4) then
        write(*,'(a)') ' '
        write(*,'(a)') 'DTRIS3:'
        write(*,'(a)') '  First tetrahedron:'
        write(*,'(2x,i8,2x,i8,2x,i8,2x,i8)') vm(1:4)
        write(*,'(a,i8,a,i8)') '  I3 = ',i3,'  I4 = ',i4
      end if
    !
    !  Insert the I-th vertex into Delaunay triangle of first I-1 vertices.
    !
      do i = 5,npt

        vi = vm(i)

        if(msglvl == 4) then
          write(*,'(a,i8,a,i8)') '  Step: ',i,'  Vertex: ',vi
        end if

        if(i == 5) then
          ip = 2
        else
          ip = i - 1
        end if

        if(i == i3 + 2) then
          ip = 3
        end if

        if(i == i4 + 1) then
          ip = 4
        end if
    !
    !  Form stacks of boundary faces involving vertex IP.
    !  TOP is for stack of boundary faces to be tested for visibility.
    !  FRONT is for stack of boundary faces visible from vertex I.
    !  TOPNV is for stack of boundary faces not visible from I.
    !
        front = 0
        topnv = 0

        if(i == 5) then

          top = 4

          if(ip == 2) then
            a = 2
          else
            a = 3
          end if

          if(ip <= 3) then
            b = 1
          else
            b = 2
          end if

          fc(7,top) = a
          fc(7,a) = b
          fc(7,b) = 0

        else if(ip == i - 1) then

          top = bfi
          fc(7,bfi) = 0
          b = fc(2,bfi)
          ptr = bf(1,-fc(5,bfi))

          do

            if(fc(1,ptr) == b) then
              b = fc(2,ptr)
              j = 1
            else
              b = fc(1,ptr)
              j = 2
            end if

            fc(7,ptr) = top
            top = ptr
            ptr = bf(j,-fc(5,ptr))

            if(ptr == bfi) then
              exit
            end if

          end do

        else

          j = 0

          do k = 1,bf_num

            if(bf(1,k) <= 0) then
              cycle
            end if

            do e = 1,3

              ptr = bf(e,k)

              if(fc(1,ptr) == ip) then
                b = fc(2,ptr)
                j = 3
                exit
              else if(fc(2,ptr) == ip) then
                b = fc(1,ptr)
                j = 3
                exit
              else if(fc(3,ptr) == ip) then
                b = fc(1,ptr)
                j = 2
                exit
              end if

            end do

            if(j /= 0) then
              exit
            end if

          end do

          bfi = ptr
          top = bfi
          fc(7,bfi) = 0
          ptr = bf(j,-fc(5,bfi))

          do

            if(fc(1,ptr) == b) then
              j = 1
              if(fc(2,ptr) == ip) then
                b = fc(3,ptr)
              else
                b = fc(2,ptr)
              end if
            else if(fc(2,ptr) == b) then
              j = 2
              if(fc(1,ptr) == ip) then
                b = fc(3,ptr)
              else
                b = fc(1,ptr)
              end if
            else
              j = 3
              if(fc(1,ptr) == ip) then
                b = fc(2,ptr)
              else
                b = fc(1,ptr)
              end if
            end if

            fc(7,ptr) = top
            top = ptr
            ptr = bf(j,-fc(5,ptr))

            if(ptr == bfi) then
              exit
            end if

          end do

        end if
    !
    !  Find a boundary face visible from vertex I.
    !
        do while(top /= 0)

          ptr = top
          top = fc(7,ptr)
          va = vm(fc(1,ptr))
          vb = vm(fc(2,ptr))
          vc = vm(fc(3,ptr))
          op = opside(vcl(1,va),vcl(1,vb),vcl(1,vc),ctr,vcl(1,vi))

          if(op == 2) then
            write(*,'(a)') ' '
            write(*,'(a)') 'DTRIS3 - Error!'
            write(*,'(a)') '  Unexpected return value from OPSIDE.'
            ierr = 301
            return
          end if

          if(op == 1) then

            front = ptr

            do while(top /= 0)

              ptr = top
              top = fc(7,ptr)
              fc(7,ptr) = -1

            end do

          else

            fc(7,ptr) = topnv
            topnv = ptr

          end if

        end do

        if(front == 0) then
          write(*,'(a)') ' '
          write(*,'(a)') 'DTRIS3 - Error!'
          write(*,'(a)') '  FRONT = 0.'
          ierr = 306
          return
        end if
    !
    !  Find remaining visible boundary faces,add new tetrahedra with
    !  vertex I,apply local transformation based on empty sphere criterion.
    !
        call vbfac(vcl(1,vi),ctr,vcl,vm,bf,fc,front,topnv)

        call nwthou(i,npt,sizht,bf_num,fc_num,bf_max,fc_max,bf,fc,ht,&
          ntetra,hdavbf,hdavfc,front,back,bfi,ierr)

        if(ierr /= 0) then
          write(*,'(a)') ' '
          write(*,'(a)') 'DTRIS3 - Error!'
          write(*,'(a,i8)') '  NWTHOU returned IERR = ',ierr
          return
        end if

        call swapes(.false.,i,npt,sizht,fc_num,fc_max,vcl,vm,bf,fc,ht,&
          ntetra,hdavfc,front,back,j,ierr)

        if(ierr /= 0) then
          write(*,'(a)') ' '
          write(*,'(a)') 'DTRIS3 - Error!'
          write(*,'(a,i8)') '  SWAPES returned IERR = ',ierr
          return
        end if

      end do

      nface = fc_num
      ptr = hdavfc

      do while(ptr /= 0)
        nface = nface - 1
        ptr = -fc(1,ptr)
      end do

      fc(7,1) = hdavbf
      fc(7,2) = hdavfc

      return
    end subroutine
    subroutine frstet(shift,nv,vcl,mapv,i3,i4,ierr)

    !*****************************************************************************80
    !
    !! FRSTET shifts vertices so the first 4 vertices are in general position in 3D.
    !
    !  Discussion: 
    !
    !    This routine shifts or swaps vertices if necessary so first 3 vertices
    !   (according to MAP) are not collinear and first 4 vertices are
    !    not coplanar(so that first tetrahedron is valid).
    !
    !  Licensing:
    !
    !    This code is distributed under the GNU LGPL license. 
    !
    !  Modified:
    !
    !    24 January 2009
    !
    !  Author:
    !
    !    Original FORTRAN77 version by Barry Joe.
    !    FORTRAN90 version by John Burkardt.
    !
    !  Reference:
    !
    !    Barry Joe,
    !    GEOMPACK - a software package for the generation of meshes
    !    using geometric algorithms,
    !    Advances in Engineering Software,
    !    Volume 13,pages 325-331,1991.
    !
    !  Parameters:
    !
    !    Input,logical :: SHIFT,if TRUE,MAP(3),MAP(4) may be updated due to shift,
    !    else they may be updated due to swaps; in former case,
    !    it is assumed MAP gives vertices in lexicographic order.
    !
    !    Input,integer :: NV,the number of vertices.
    !
    !    Input,real(dp) :: VCL(1:3,1:*),the vertex coordinate list.
    !
    !    Input/output,integer :: MAP(1:NV),on input,contains vertex 
    !    indices of VCL.  On output,shifted or 2 swaps applied if necessary so 
    !    that vertices indexed by MAP(1),MAP(2),MAP(3),MAP(4) not coplanar.
    !
    !    Output,integer :: I3,I4,the indices such that
    !    MAP_in(I3) = MAP_out(3) and MAP_in(I4) = MAP_out(4).
    !
    !    Output,integer :: IERR,error flag,which is zero unless 
    !    an error occurred.
    !
      implicit none

      integer :: nv

      real(dp) :: cmax
      real(dp) :: cp1
      real(dp) :: cp2
      real(dp) :: cp3
      real(dp) :: dmax
      real(dp) :: dotp
      real(dp) :: dv2(3)
      real(dp) :: dvk(3)
      real(dp) :: dvl(3)
      integer :: i
      integer :: i3
      integer :: i4
      integer :: ierr
      integer :: k
      integer :: l
      integer :: m
      integer :: m1
      integer :: m2
      integer :: mapv(nv)
      logical :: shift
      real(dp) :: tol
      real(dp) :: vcl(3,*)
    !
    !  First check that consecutive vertices are not identical.
    !
      ierr = 0
      tol = 100.0D+00 * epsilon(tol)

      if(shift) then
        l = nv - 1
      else
        l = 1
      end if

      m1 = mapv(1)

      do i = 1,l

        m = m1
        m1 = mapv(i+1)

        do k = 1,3
          cmax = max(abs(vcl(k,m)),abs(vcl(k,m1)))
          if(tol * cmax < abs(vcl(k,m) - vcl(k,m1)) .and. tol < cmax) then
            go to 20
          end if
        end do

        ierr = 302
        return

    20  continue

      end do
    !
    !  Find index K = I3 and L = I4.
    !
      m1 = mapv(1)
      m2 = mapv(2)

      dv2(1:3) = vcl(1:3,m2) - vcl(1:3,m1)

      cmax = max(abs(vcl(1,m1)),abs(vcl(2,m1)),abs(vcl(3,m1)),&
        abs(vcl(1,m2)),abs(vcl(2,m2)),abs(vcl(3,m2)))
      k = 2

      do

        k = k + 1

        if(nv < k) then
          ierr = 303
          return
        end if

        m = mapv(k)

        dvk(1:3) = vcl(1:3,m) - vcl(1:3,m1)

        dmax = max(cmax,abs(vcl(1,m)),abs(vcl(2,m)),abs(vcl(3,m)))

        cp1 = dv2(2) * dvk(3) - dv2(3) * dvk(2)
        cp2 = dv2(3) * dvk(1) - dv2(1) * dvk(3)
        cp3 = dv2(1) * dvk(2) - dv2(2) * dvk(1)

        if(tol * dmax < max(abs(cp1),abs(cp2),abs(cp3))) then
          exit
        end if

      end do

      cmax = dmax
      l = k

      do

        l = l + 1

        if(nv < l) then
          ierr = 304
          return
        end if

        m = mapv(l)

        dvl(1:3) = vcl(1:3,m) - vcl(1:3,m1)

        dmax = max(cmax,abs(vcl(1,m)),abs(vcl(2,m)),abs(vcl(3,m)))

        dotp = dvl(1) * cp1 + dvl(2) * cp2 + dvl(3) * cp3

        if(tol * dmax < abs(dotp)) then
          exit
        end if

      end do
    !
    !  Shift or swap elements of MAP if necessary.
    !
      if(shift) then

        if(3 < k) then
          m1 = mapv(k)
        end if

        if(4 < l) then
          m2 = mapv(l)
          do i = l,k+2,-1
             mapv(i) = mapv(i-1)
          end do
          do i = k+1,5,-1
            mapv(i) = mapv(i-2)
          end do
          mapv(4) = m2
        end if

        if(3 < k) then
          mapv(3) = m1
        end if

      else

        if(3 < k) then
          m = mapv(3)
          mapv(3) = mapv(k)
          mapv(k) = m
        end if

        if(4 < l) then
          m = mapv(4)
          mapv(4) = mapv(l)
          mapv(l) = m
        end if

      end if

      i3 = k
      i4 = l

      return
    end subroutine
    subroutine htdel(ind,n,p,fc,ht)

    !*****************************************************************************80
    !
    !! HTDEL deletes a record from the hash table.
    !
    !  Discussion: 
    !
    !    This routine deletes record FC(1:7,IND) from the hash table HT.
    !
    !  Licensing:
    !
    !    This code is distributed under the GNU LGPL license. 
    !
    !  Modified:
    !
    !    24 January 2009
    !
    !  Author:
    !
    !    Original FORTRAN77 version by Barry Joe.
    !    FORTRAN90 version by John Burkardt.
    !
    !  Parameters:
    !
    !    Input,integer :: IND,the index of FC array.
    !
    !    Input,integer :: N,the upper bound on vertex indices.
    !
    !    Input,integer :: P,the size of hash table.
    !
    !    Input/output,integer :: FC(1:7,1:*),the array of face records; 
    !    see routine DTRIS3.  On output,one link in FC is updated.
    !
    !    Input/output,integer :: HT(0:P-1),the hash table using 
    !    direct chaining.  On output,one link in HT is updated.
    !
      implicit none

      integer :: p

      integer :: fc(7,*)
      integer :: ht(0:p-1)
      integer :: ind
      integer :: k
      integer :: n
      integer :: ptr

      k = mod(fc(1,ind) * n + fc(2,ind),p)
      k = mod(k * n + fc(3,ind),p)
      ptr = ht(k)

      if(ptr == ind) then

        ht(k) = fc(6,ind)

      else

        do while(fc(6,ptr) /= ind)
          ptr = fc(6,ptr)
        end do
        fc(6,ptr) = fc(6,ind)

      end if

      return
    end subroutine
    subroutine htins(ind,a,b,c,d,e,n,p,fc,ht)

    !*****************************************************************************80
    !
    !! HTINS inserts a record into the hash table.
    !
    !  Discussion: 
    !
    !    This routine inserts record FC(1:7,IND) containing A,B,C,D,E,HTLINK,-1
    !    into hash table HT.
    !
    !  Licensing:
    !
    !    This code is distributed under the GNU LGPL license. 
    !
    !  Modified:
    !
    !    24 January 2009
    !
    !  Author:
    !
    !    Original FORTRAN77 version by Barry Joe.
    !    FORTRAN90 version by John Burkardt.
    !
    !  Parameters:
    !
    !    Input,integer :: IND,the index of FC array.
    !
    !    Input,integer :: A,B,C,D,E,the first 5 fields of FC 
    !    record(or column).
    !
    !    Input,integer :: N,the upper bound on vertex indices.
    !
    !    Input,integer :: P,the size of hash table.
    !
    !    Input/output,integer :: FC(1:7,1:*),the array of face records; 
    !    see routine DTRIS3.
    !
    !    Input/output,integer :: HT(0:P-1),the hash table using 
    !    direct chaining.
    !
      implicit none

      integer :: p

      integer :: a
      integer :: aa
      integer :: b
      integer :: bb
      integer :: c
      integer :: cc
      integer :: d
      integer :: e
      integer :: fc(7,*)
      integer :: ht(0:p-1)
      integer :: ind
      integer :: k
      integer :: n

      aa = a
      bb = b
      cc = c
      call order3(aa,bb,cc)
      k = mod(aa * n + bb,p)
      k = mod(k * n + cc,p)
      fc(1,ind) = aa
      fc(2,ind) = bb
      fc(3,ind) = cc
      fc(4,ind) = d
      fc(5,ind) = e
      fc(6,ind) = ht(k)
      fc(7,ind) = -1
      ht(k) = ind

      return
    end subroutine
    function htsrc(a,b,c,n,p,fc,ht)

    !*****************************************************************************80
    !
    !! HTSRC searches for a record in the hash table.
    !
    !  Discussion: 
    !
    !    This routine searches for record FC(1:7,IND) containing key A,B,C
    !    in hash table HT.
    !
    !  Licensing:
    !
    !    This code is distributed under the GNU LGPL license. 
    !
    !  Modified:
    !
    !    24 January 2009
    !
    !  Author:
    !
    !    Original FORTRAN77 version by Barry Joe.
    !    FORTRAN90 version by John Burkardt.
    !
    !  Parameters:
    !
    !    Input,integer :: A,B,C,first 3 fields of FC record 
    !   (in any order).
    !
    !    Input,integer :: N,upper bound on vertex indices.
    !
    !    Input,integer :: P,size of hash table.
    !
    !    Input,integer :: FC(1:7,1:*),array of face records; 
    !    see routine DTRIS3.
    !
    !    Input,integer :: HT(0:P-1),hash table using direct chaining.
    !
    !    Output,integer :: HTSRC,index of FC record with key A,B,C 
    !    if found,or 0 if not found.
    !
      implicit none

      integer :: p

      integer :: a
      integer :: aa
      integer :: b
      integer :: bb
      integer :: c
      integer :: cc
      integer :: fc(7,*)
      integer :: ht(0:p-1)
      integer :: htsrc
      integer :: ind
      integer :: k
      integer :: n

      aa = a
      bb = b
      cc = c
      call order3(aa,bb,cc)
      k = mod(aa * n + bb,p)
      k = mod(k * n + cc,p)
      ind = ht(k)

      do

        if(ind == 0) then
          exit
        end if

        if(fc(1,ind) == aa .and. fc(2,ind) == bb .and. fc(3,ind) == cc) then
          exit
        end if

        ind = fc(6,ind)

      end do

      htsrc = ind

      return
    end function
    function i4_modp(i,j)

    !*****************************************************************************80
    !
    !! I4_MODP returns the nonnegative remainder of integer division.
    !
    !  Discussion:
    !
    !    If
    !      NREM = I4_MODP(I,J)
    !      NMULT =(I - NREM) / J
    !    then
    !      I = J * NMULT + NREM
    !    where NREM is always nonnegative.
    !
    !    The MOD function computes a result with the same sign as the
    !    quantity being divided.  Thus,suppose you had an angle A,
    !    and you wanted to ensure that it was between 0 and 360.
    !    Then mod(A,360) would do,if A was positive,but if A
    !    was negative,your result would be between -360 and 0.
    !
    !    On the other hand,I4_MODP(A,360) is between 0 and 360,always.
    !
    !  Example:
    !
    !        I     J     MOD  I4_MODP    Factorization
    !
    !      107    50       7       7    107 =  2 *  50 + 7
    !      107   -50       7       7    107 = -2 * -50 + 7
    !     -107    50      -7      43   -107 = -3 *  50 + 43
    !     -107   -50      -7      43   -107 =  3 * -50 + 43
    !
    !  Licensing:
    !
    !    This code is distributed under the GNU LGPL license. 
    !
    !  Modified:
    !
    !    02 March 1999
    !
    !  Author:
    !
    !    John Burkardt
    !
    !  Parameters:
    !
    !    Input,integer :: I,the number to be divided.
    !
    !    Input,integer :: J,the number that divides I.
    !
    !    Output,integer :: I4_MODP,the nonnegative remainder when I is
    !    divided by J.
    !
      implicit none

      integer :: i
      integer :: i4_modp
      integer :: j

      if(j == 0) then
        write(*,'(a)') ' '
        write(*,'(a)') 'I4_MODP - Fatal error!'
        write(*,'(a,i8)') '  I4_MODP(I,J) called with J = ',j
        stop
      end if

      i4_modp = mod(i,j)

      if(i4_modp < 0) then
        i4_modp = i4_modp + abs(j)
      end if

      return
    end function
    subroutine i4_swap(i,j)

    !*****************************************************************************80
    !
    !! I4_SWAP swaps two I4's.
    !
    !  Licensing:
    !
    !    This code is distributed under the GNU LGPL license. 
    !
    !  Modified:
    !
    !    30 November 1998
    !
    !  Author:
    !
    !    John Burkardt
    !
    !  Parameters:
    !
    !    Input/output,integer :: I,J.  On output,the values of I and
    !    J have been interchanged.
    !
      implicit none

      integer :: i
      integer :: j
      integer :: k

      k = i
      i = j
      j = k

      return
    end subroutine
    function i4_wrap(ival,ilo,ihi)

    !*****************************************************************************80
    !
    !! I4_WRAP forces an integer to lie between given limits by wrapping.
    !
    !  Example:
    !
    !    ILO = 4,IHI = 8
    !
    !    I  I4_WRAP
    !
    !    -2     8
    !    -1     4
    !     0     5
    !     1     6
    !     2     7
    !     3     8
    !     4     4
    !     5     5
    !     6     6
    !     7     7
    !     8     8
    !     9     4
    !    10     5
    !    11     6
    !    12     7
    !    13     8
    !    14     4
    !
    !  Licensing:
    !
    !    This code is distributed under the GNU LGPL license. 
    !
    !  Modified:
    !
    !    19 August 2003
    !
    !  Author:
    !
    !    John Burkardt
    !
    !  Parameters:
    !
    !    Input,integer :: IVAL,an integer value.
    !
    !    Input,integer :: ILO,IHI,the desired bounds.
    !
    !    Output,integer :: I4_WRAP,a "wrapped" version of IVAL.
    !
      implicit none

   !   integer :: i4_modp
      integer :: i4_wrap
      integer :: ihi
      integer :: ilo
      integer :: ival
      integer :: jhi
      integer :: jlo
      integer :: wide

      jlo = min(ilo,ihi)
      jhi = max(ilo,ihi)

      wide = jhi - jlo + 1

      if(wide == 1) then
        i4_wrap = jlo
      else
        i4_wrap = jlo + i4_modp(ival - jlo,wide)
      end if

      return
    end function i4_wrap
    subroutine nwthou(i,npt,sizht,bf_num,fc_num,bf_max,fc_max,bf,fc,ht,&
      ntetra,hdavbf,hdavfc,front,back,bfi,ierr)

    !*****************************************************************************80
    !
    !! NWTHOU creates new tetrahedra outside the current convex hull.
    !
    !  Discussion: 
    !
    !    This routine creates new tetrahedra in a 3D triangulation outside the
    !    convex hull by joining vertex I to visible boundary faces.
    !
    !  Licensing:
    !
    !    This code is distributed under the GNU LGPL license. 
    !
    !  Modified:
    !
    !    24 January 2009
    !
    !  Author:
    !
    !    Original FORTRAN77 version by Barry Joe.
    !    FORTRAN90 version by John Burkardt.
    !
    !  Parameters:
    !
    !    Input,integer :: I,the(local) index of next vertex inserted in
    !    triangulation.
    !
    !    Input,integer :: NPT,the number of 3D points to be 
    !    triangulated.
    !
    !    Input,integer :: SIZHT,the size of hash table HT.
    !
    !    Input/output,integer :: BF_NUM,the number of positions in BF.
    !
    !    Input/output,integer :: FC_NUM,the number of positions in FC.
    !
    !    Input,integer :: BF_MAX,the maximum size available for BF.
    !
    !    Input,integer :: FC_MAX,the maximum size available for FC.
    !
    !    Input/output,integer :: BF(1:3,1:BF_MAX),the array of 
    !    boundary face records; see DTRIS3.
    !
    !    Input/output,integer :: FC(1:7,1:FC_MAX),the array of face 
    !    records; see routine DTRIS3.
    !
    !    Input/output,integer :: HT(0:SIZHT-1),the hash table using 
    !    direct chaining.
    !
    !    Input/output,integer :: NTETRA,the number of tetrahedra in 
    !    triangulation.
    !
    !    Input/output,integer :: HDAVBF,the head pointer to available 
    !    BF records.
    !
    !    Input/output,integer :: HDAVFC,the head pointer to available 
    !    FC records.
    !
    !    Input,integer :: FRONT,the index of front of queue(or top 
    !    of stack) of visible boundary faces.
    !
    !    Output,integer :: BACK,the index of back of queue(or bottom 
    !    of stack) of visible boundary faces(which become interior faces).
    !
    !    Output,integer :: BFI,the index of FC of a boundary face 
    !    containing vertex I.
    !
    !    Output,integer :: IERR,error flag,which is zero unless an 
    !    error occurred.
    !
      implicit none

      integer :: bf_max
      integer :: fc_max
      integer :: sizht

      integer :: a
      integer :: b
      integer :: back
      integer :: bf(3,bf_max)
      integer :: bfi
      integer :: bfnew
      integer :: c
      integer :: d
      integer :: e
      integer :: fc(7,fc_max)
      integer :: fc_num
      integer :: front
      integer :: hdavbf
      integer :: hdavfc
      integer :: ht(0:sizht-1)
   !    integer :: htsrc
      integer :: i
      integer :: ierr
      integer :: ind
      integer :: j
      integer :: k
      integer :: l
      integer,parameter :: msglvl = 0
      integer :: bf_num
      integer :: nbr
      integer :: npt
      integer :: ntetra
      integer :: ptr
    !
    !  For ABC in queue,form tetrahedron ABCI + add faces ABI,ACI,BCI.
    !  PTR,NBR,IND are indices of FC; K,L,BFNEW indices of BF.
    !
      ierr = 0
      bfi = 0
      ptr = front

      do

        back = ptr
        a = fc(1,ptr)
        b = fc(2,ptr)
        c = fc(3,ptr)
        k = -fc(5,ptr)
        fc(5,ptr) = i
        ntetra = ntetra + 1

        if(msglvl == 4) then
          write(*,600) a,b,c,i
        end if

        do e = 1,3

          if(e == 2) then
            call i4_swap(a,b)
          else if(e == 3) then
            call i4_swap(a,c)
          end if

          nbr = bf(e,k)

          if(fc(7,nbr) /= -1) then
            if(fc(5,nbr) == i) then
              cycle
            end if
          end if

          call availf(hdavfc,fc_num,fc_max,fc,ind,ierr)

          if(ierr /= 0) then
            write(*,'(a)') ' '
            write(*,'(a)') 'NWTHOU - Error!'
            write(*,'(a,i8)') '  AVAILF returned IERR = ',ierr
            return
          end if

          l = -fc(5,nbr)

          if(bf(1,l) == ptr) then
            j = 1
          else if(bf(2,l) == ptr) then
            j = 2
          else
            j = 3
          end if

          if(fc(7,nbr) /= -1) then

            call htins(ind,b,c,i,a,fc(j,nbr),npt,sizht,fc,ht)

          else

            if(hdavbf /= 0) then

              bfnew = hdavbf
              hdavbf = -bf(1,hdavbf)

            else

              if(bf_max <= bf_num) then
                write(*,'(a)') ' '
                write(*,'(a)') 'NWTHOU - Error!'
                write(*,'(a)') '  BF_MAX <= BF_NUM.'
                write(*,'(a)') '  BF_MAX must be increased to proceed.'
                ierr = 12
                return
              end if

              bf_num = bf_num + 1
              bfnew = bf_num

            end if

            if(bfi == 0) then
              bfi = ind
            end if

            call htins(ind,b,c,i,a,-bfnew,npt,sizht,fc,ht)
            bf(j,l) = ind
            bf(3,bfnew) = nbr

          end if

        end do

        if(k == bf_num) then
          bf_num = bf_num - 1
        else
          bf(1,k) = -hdavbf
          hdavbf = k
        end if

        ptr = fc(7,ptr)

        if(ptr == 0) then
          exit
        end if

      end do
    !
    !  Set BF(1:2,BFNEW) fields for new boundary faces.
    !
      ptr = bfi
      a = fc(1,ptr)
      j = 2

      do

        b = fc(j,ptr)
        c = fc(4,ptr)

        do

          nbr = htsrc(a,c,i,npt,sizht,fc,ht)
 
          if(nbr <= 0) then
            write(*,'(a)') ' '
            write(*,'(a)') 'NWTHOU - Error!'
            write(*,'(a,i8)') '  HTSRC returned NBR <= 0.'
            ierr = 300
            return
          end if

          if(fc(5,nbr) <= 0) then
            exit
          end if

          if(fc(4,nbr) == b) then
            d = fc(5,nbr)
          else
            d = fc(4,nbr)
          end if

          b = c
          c = d

        end do

        k = -fc(5,ptr)
        l = -fc(5,nbr)

        if(fc(1,ptr) == a) then
          bf(2,k) = nbr
        else
          bf(1,k) = nbr
        end if

        if(fc(1,nbr) == a) then
          j = 1
        else
          j = 2
        end if

        bf(3-j,l) = ptr
        a = fc(3-j,nbr)
        ptr = nbr

        if(ptr == bfi) then
          exit
        end if

      end do

      600 format('  New tetra: ',4i7)

      return
    end subroutine
    function opside(a,b,c,d,e)

    !*****************************************************************************80
    !
    !! OPSIDE tests if points are on opposite sides of a triangular face.
    !
    !  Discussion: 
    !
    !    This routine tests if points D,E are on opposite sides of triangular
    !    face with vertices A,B,C.
    !
    !  Licensing:
    !
    !    This code is distributed under the GNU LGPL license. 
    !
    !  Modified:
    !
    !    31 August 2005
    !
    !  Author:
    !
    !    Original FORTRAN77 version by Barry Joe.
    !    FORTRAN90 version by John Burkardt.
    !
    !  Parameters:
    !
    !    Input,real(dp) :: A(1:3),B(1:3),C(1:3),D(1:3),E(1:3),
    !    five 3D points.
    !
    !    Output,integer :: OPSIDE,the result of the test:
    !    +1 if D,E on opposite sides; 
    !    -1 if on same side;
    !     2 if D is coplanar with face ABC(ABCD is a degenerate tetrahedron); 
    !     0 if E is coplanar with face ABC
    !
      implicit none

      real(dp) :: a(3)
      real(dp) :: ab(3)
      real(dp) :: ac(3)
      real(dp) :: b(3)
      real(dp) :: c(3)
      real(dp) :: d(3)
      real(dp) :: ddp
      real(dp) :: dmax
      real(dp) :: e(3)
      real(dp) :: edp
      real(dp) :: emax
      real(dp) :: nrml1
      real(dp) :: nrml2
      real(dp) :: nrml3
      integer :: opside
      real(dp) :: tol

      tol = 100.0D+00 * epsilon(tol)

      ab(1:3) = b(1:3) - a(1:3)
      ac(1:3) = c(1:3) - a(1:3)

      emax = max(&
        abs(a(1)),abs(a(2)),abs(a(3)),&
        abs(b(1)),abs(b(2)),abs(b(3)),&
        abs(c(1)),abs(c(2)),abs(c(3)))

      dmax = max(emax,abs(d(1)),abs(d(2)),abs(d(3)))

      nrml1 = ab(2) * ac(3) - ab(3) * ac(2)
      nrml2 = ab(3) * ac(1) - ab(1) * ac(3)
      nrml3 = ab(1) * ac(2) - ab(2) * ac(1)

      ddp =(d(1) - a(1)) * nrml1 &
          +(d(2) - a(2)) * nrml2 &
          +(d(3) - a(3)) * nrml3

      if(abs(ddp) <= tol * dmax) then
        opside = 2
        return
      end if

      emax = max(emax,abs(e(1)),abs(e(2)),abs(e(3)))

      edp =(e(1) - a(1)) * nrml1 &
          +(e(2) - a(2)) * nrml2 &
          +(e(3) - a(3)) * nrml3

      if(abs(edp) <= tol * emax) then
        opside = 0
      else if(ddp * edp < 0.0D+00) then
        opside = 1
      else
        opside = -1
      end if

      return
    end function
    subroutine order3(i,j,k)

    !*****************************************************************************80
    !
    !! ORDER3 reorders 3 integers into ascending order.
    !
    !  Discussion: 
    !
    !    This routine reorders I,J,K so that I <= J <= K.
    !
    !  Licensing:
    !
    !    This code is distributed under the GNU LGPL license. 
    !
    !  Modified:
    !
    !    24 January 2009
    !
    !  Author:
    !
    !    Original FORTRAN77 version by Barry Joe.
    !    FORTRAN90 version by John Burkardt.
    !
    !  Parameters:
    !
    !    Input/output,integer :: I,J,K,on output are sorted into
    !    nondecreasing order.
    !
      implicit none

      integer :: i
      integer :: j
      integer :: k
      integer :: t

      if(j < i) then
        if(k < j) then
          call i4_swap(i,k)
        else if(k < i) then
          t = i
          i = j
          j = k
          k = t
        else
          call i4_swap(i,j)
        end if
      else
        if(k < i) then
          t = i
          i = k
          k = j
          j = t
        else if(k < j) then
          call i4_swap(j,k)
        end if
      end if

      return
    end subroutine
    subroutine swapes(bndcon,i,npt,sizht,fc_num,fc_max,vcl,vm,bf,fc,&
      ht,ntetra,hdavfc,front,back,ifac,ierr)

    !*****************************************************************************80
    !
    !! SWAPES swaps faces in a 3D triangulation.
    !
    !  Discussion: 
    !
    !    This routine swaps faces,applying local transformations,in a 3D 
    !    triangulation based on the empty circumsphere criterion until(nearly)
    !    all faces are locally optimal.  I is the index of the new vertex
    !    added to the triangulation,or 0 if an initial triangulation is given.
    !
    !  Licensing:
    !
    !    This code is distributed under the GNU LGPL license. 
    !
    !  Modified:
    !
    !    07 September 2005
    !
    !  Author:
    !
    !    Original FORTRAN77 version by Barry Joe.
    !    FORTRAN90 version by John Burkardt.
    !
    !  Parameters:
    ! 
    !    Input,logical :: BNDCON,TRUE iff boundary faces are constrained(i.e. not
    !    swapped by local transformations).
    !
    !    Input,integer :: I,the local index of next vertex inserted in
    !    triangulation,or 0; if positive,it is assumed I is largest index so far.
    !
    !    Input,integer :: NPT,the number of 3D points to 
    !    be triangulated.
    !
    !    Input,integer :: SIZHT,the size of hash table HT.
    !
    !    Input/output,integer :: FC_NUM,the number of positions used 
    !    in FC array.
    !
    !    Input,integer :: FC_MAX,the maximum size available for FC.
    !
    !    Input,real(dp) :: VCL(1:3,1:*),the vertex coordinate list.
    !
    !    Input,integer :: VM(1:NPT),the indices of vertices of VCL being
    !    triangulated.
    !
    !    Input/output,integer :: BF(1:3,1:*),the  array of boundary 
    !    face records; see DTRIS3.
    !
    !    Input/output,integer :: FC(1:7,1:FC_MAX),the array of 
    !    face records; see routine DTRIS3.
    !
    !    Input/output,integer :: HT(0:SIZHT-1),the hash table using 
    !    direct chaining.
    !
    !    Input/output,integer :: NTETRA,the number of tetrahedra 
    !    in triangulation.
    !
    !    Input/output,integer :: HDAVFC,the head pointer to available 
    !    FC records.
    !
    !    Input/output,integer :: FRONT,BACK,the indices of front and 
    !    back of queue of interior faces for which sphere test is applied.
    !
    !    Output,integer :: IFAC,the index of last face for which sphere
    !    test applied,or 0.
    !
    !    Output,integer :: IERR,error flag,which is zero unless an 
    !    error occurred.
    !
      implicit none

      integer :: fc_max
      integer :: npt
      integer :: sizht

      integer :: a
      integer :: aa
      real(dp) :: alpha(4)
      integer :: b
      integer :: back
      integer :: bf(3,*)
      integer :: bfx(2)
      logical :: bndcon
      integer :: c
      real(dp) :: center(3)
      integer :: d
      integer :: dd
      logical :: degen
      integer :: e
      integer :: f
      integer :: fc(7,fc_max)
      integer :: fc_num
      integer :: front
      integer :: g
      integer :: hdavfc
      integer :: ht(0:sizht-1)
 !     integer :: htsrc
      integer :: i
      integer :: ierr
      integer :: ifac
      integer :: in
      integer :: ind
      integer :: ind1
      integer :: ind2
      integer :: indx(2)
      integer :: j
      integer :: k
      integer :: kneg
      integer :: kzero
      integer,parameter :: msglvl = 0
      integer :: nbr(2,2)
      integer :: ntetra
      real(dp) :: radsq
      real(dp) :: vcl(3,*)
      integer :: va
      integer :: vb
      integer :: vc
      integer :: vd
      integer :: ve
      integer :: vm(npt)

      ierr = 0
      ifac = 0

      do

        do

          if(front == 0) then
            return
          end if

          ind = front
          front = fc(7,ind)

          if(fc(2,ind) /= 0) then
            exit
          end if

          if(ind == fc_num) then
            fc_num = fc_num - 1
          else
            fc(1,ind) = -hdavfc
            hdavfc = ind
          end if

        end do

        ifac = ind
        fc(7,ind) = -1
        a = fc(1,ind)
        b = fc(2,ind)
        c = fc(3,ind)
        d = fc(4,ind)
        e = fc(5,ind)
        va = vm(a)
        vb = vm(b)
        vc = vm(c)
        vd = vm(d)
        ve = vm(e)

        if(msglvl == 4) then
          write(*,600) ind,a,b,c,d,e
        end if

        call ccsph(.true.,vcl(1,va),vcl(1,vb),vcl(1,vc),vcl(1,vd),&
          vcl(1,ve),center,radsq,in)

        if(in == 2) then
          write(*,'(a)') ' '
          write(*,'(a)') 'SWAPES - Fatal error!'
          write(*,'(a)') '  CCSPH returned IN = 2.'
          ierr = 301
          return
        end if

        if(1 <= in) then

          call baryth(vcl(1,va),vcl(1,vb),vcl(1,vc),vcl(1,vd),&
            vcl(1,ve),alpha,degen)

          if(degen) then
            write(*,'(a)') ' '
            write(*,'(a)') 'SWAPES - Fatal error!'
            write(*,'(a)') '  BARYTH detected a degenerate tetrahedron.'
            ierr = 301
            return
          else if(0.0D+00 < alpha(4)) then
            write(*,'(a)') ' '
            write(*,'(a)') 'SWAPES - Fatal error!'
            write(*,'(a)') '  BARYTH detected 0 < ALPHA(4).'
            ierr = 309
            return
          end if

          kneg = 1
          kzero = 0

          do j = 1,3
            if(alpha(j) < 0.0D+00) then
              kneg = kneg + 1
            else if(alpha(j) == 0.0D+00) then
              kzero = kzero + 1
            end if
          end do
    !
    !  Swap 2 tetrahedra for 3.
    !
          if(kneg == 1 .and. kzero == 0) then

            call updatf(a,b,d,c,e,i,npt,sizht,front,back,fc,ht,ierr)
            call updatf(a,c,d,b,e,i,npt,sizht,front,back,fc,ht,ierr)
            call updatf(b,c,d,a,e,i,npt,sizht,front,back,fc,ht,ierr)
            call updatf(a,b,e,c,d,i,npt,sizht,front,back,fc,ht,ierr)
            call updatf(a,c,e,b,d,i,npt,sizht,front,back,fc,ht,ierr)
            call updatf(b,c,e,a,d,i,npt,sizht,front,back,fc,ht,ierr)

            if(ierr /= 0) then
              write(*,'(a)') ' '
              write(*,'(a)') 'SWAPES - Fatal error!'
              write(*,'(a,i8)') '  UPDATF returned IERR = ',ierr
              return
            end if

            call htdel(ind,npt,sizht,fc,ht)
            call htins(ind,a,d,e,b,c,npt,sizht,fc,ht)
            call availf(hdavfc,fc_num,fc_max,fc,ind,ierr)

            if(ierr /= 0) then
              write(*,'(a)') ' '
              write(*,'(a)') 'SWAPES - Fatal error!'
              write(*,'(a,i8)') '  AVAILF returned IERR = ',ierr
              return
            end if

            call htins(ind,b,d,e,a,c,npt,sizht,fc,ht)
            call availf(hdavfc,fc_num,fc_max,fc,ind,ierr)

            if(ierr /= 0) then
              write(*,'(a)') ' '
              write(*,'(a)') 'SWAPES - Fatal error!'
              write(*,'(a,i8)') '  AVAILF returned IERR = ',ierr
              return
            end if

            call htins(ind,c,d,e,a,b,npt,sizht,fc,ht)
            ntetra = ntetra + 1

            if(msglvl == 4) then
              write(*,610)
            end if
    !
    !  Swap 3 tetrahedra for 2 if possible. Relabel so edge
    !  AB would be deleted. Swap if ABDE is in current triangulation.
    !
          else if(kneg == 2 .and. kzero == 0) then

            if(alpha(1) < 0.0D+00) then
              call i4_swap(a,c)
            else if(alpha(2) < 0.0D+00) then
              call i4_swap(b,c)
            end if

            ind1 = htsrc(a,b,d,npt,sizht,fc,ht)

            if(ind1 <= 0) then
              write(*,'(a)') ' '
              write(*,'(a)') 'SWAPES - Fatal error!'
              write(*,'(a)') '  HTSRC returned IND1 <= 0.'
              ierr = 300
              return
            end if

            if(fc(4,ind1) == e .or. fc(5,ind1) == e) then

              call updatf(a,c,d,b,e,i,npt,sizht,front,back,fc,ht,&
                ierr)
              call updatf(a,c,e,b,d,i,npt,sizht,front,back,fc,ht,&
                ierr)
              call updatf(a,d,e,b,c,i,npt,sizht,front,back,fc,ht,&
                ierr)
              call updatf(b,c,d,a,e,i,npt,sizht,front,back,fc,ht,&
                ierr)
              call updatf(b,c,e,a,d,i,npt,sizht,front,back,fc,ht,&
                ierr)
              call updatf(b,d,e,a,c,i,npt,sizht,front,back,fc,ht,&
                ierr)

              if(ierr /= 0) then
                write(*,'(a)') ' '
                write(*,'(a)') 'SWAPES - Fatal error!'
                write(*,'(a,i8)') '  UPDATF returned IERR = ',ierr
                return
              end if

              call htdel(ind,npt,sizht,fc,ht)
              call htins(ind,c,d,e,a,b,npt,sizht,fc,ht)
              call htdel(ind1,npt,sizht,fc,ht)

              if(0 <= fc(7,ind1)) then
                fc(2,ind1) = 0
              else
                if(ind1 == fc_num) then
                  fc_num = fc_num - 1
                else
                  fc(1,ind1) = -hdavfc
                  hdavfc = ind1
                end if
              end if

              ind1 = htsrc(a,b,e,npt,sizht,fc,ht)

              if(ind1 <= 0) then
                write(*,'(a)') ' '
                write(*,'(a)') 'SWAPES - Fatal error!'
                write(*,'(a)') '  HTSRC returned IND1 <= 0.'
                ierr = 300
                return
              end if

              call htdel(ind1,npt,sizht,fc,ht)

              if(0 <= fc(7,ind1)) then

                fc(2,ind1) = 0

              else

                if(ind1 == fc_num) then
                  fc_num = fc_num - 1
                else
                  fc(1,ind1) = -hdavfc
                  hdavfc = ind1
                end if

              end if

              ntetra = ntetra - 1

              if(msglvl == 4) then
                write(*,620) c,d,e
              end if

            else

              if(msglvl == 4) then
                write(*,630) a,b,d,e
              end if

            end if
    !
    !  Coplanar faces: swap 2 tetrahedra for 2 if boundary faces
    ! (and BNDCON is .FALSE.),else do pair of 2 for 2 swaps if
    !  possible.  Relabel vertices so that DE intersects AB.
    !  Also swap if necessary to make A < B and D < E.
    !
          else if(kneg == 1 .and. kzero == 1) then

            if(alpha(1) == 0.0D+00) then

              call i4_swap(a,c)

            else if(alpha(2) == 0.0D+00) then

              call i4_swap(b,c)

            end if

            if(b < a) then
              call i4_swap(a,b)
            end if

            if(e < d) then
              call i4_swap(d,e)
            end if

            ind1 = htsrc(a,b,d,npt,sizht,fc,ht)

            if(ind1 <= 0) then
              write(*,'(a)') ' '
              write(*,'(a)') 'SWAPES - Fatal error!'
              write(*,'(a)') '  HTSRC returned IND1 <= 0.'
              ierr = 300
              return
            end if

            ind2 = htsrc(a,b,e,npt,sizht,fc,ht)

            if(ind2 <= 0) then
              write(*,'(a)') ' '
              write(*,'(a)') 'SWAPES - Fatal error!'
              write(*,'(a)') '  HTSRC returned IND2 <= 0.'
              ierr = 300
              return
            end if

            if(fc(4,ind1) == c) then
              f = fc(5,ind1)
            else
              f = fc(4,ind1)
            end if

            if(fc(4,ind2) == c) then
              g = fc(5,ind2)
            else
              g = fc(4,ind2)
            end if

            if(f <= 0 .and. g <= 0) then

              if(.not. bndcon) then

                call updatf(a,c,d,b,e,i,npt,sizht,front,back,fc,ht,ierr)
                call updatf(a,c,e,b,d,i,npt,sizht,front,back,fc,ht,ierr)
                call updatf(b,c,d,a,e,i,npt,sizht,front,back,fc,ht,ierr)
                call updatf(b,c,e,a,d,i,npt,sizht,front,back,fc,ht,ierr)

                if(ierr /= 0) then
                  write(*,'(a)') ' '
                  write(*,'(a)') 'SWAPES - Fatal error!'
                  write(*,'(a,i8)') '  UPDATF returned IERR = ',ierr
                  return
                end if

                call htdel(ind,npt,sizht,fc,ht)
                call htins(ind,c,d,e,a,b,npt,sizht,fc,ht)
                call htdel(ind1,npt,sizht,fc,ht)
                call htins(ind1,a,d,e,c,fc(5,ind1),npt,sizht,fc,ht)
                call htdel(ind2,npt,sizht,fc,ht)
                call htins(ind2,b,d,e,c,fc(5,ind2),npt,sizht,fc,ht)
                indx(1) = ind1
                indx(2) = ind2
                bfx(1) = -fc(5,ind1)
                bfx(2) = -fc(5,ind2)
                dd = d

                do j = 1,2

                  if(j == 2) then
                    dd = e
                  end if

                  if(dd < a) then
                    nbr(j,1) = bf(3,bfx(j))
                    nbr(j,2) = bf(2,bfx(j))
                  else if(dd < b) then
                    nbr(j,1) = bf(3,bfx(j))
                    nbr(j,2) = bf(1,bfx(j))
                  else
                    nbr(j,1) = bf(2,bfx(j))
                    nbr(j,2) = bf(1,bfx(j))
                  end if

                end do

                aa = a
                k = -fc(5,nbr(1,2))
  
                do j = 1,2

                  if(j == 2) then
                    aa = b
                    k = -fc(5,nbr(2,1))
                  end if

                  if(aa < d) then
                    bf(1,bfx(j)) = indx(3-j)
                    bf(2,bfx(j)) = nbr(2,j)
                    bf(3,bfx(j)) = nbr(1,j)
                  else if(aa < e) then
                    bf(1,bfx(j)) = nbr(2,j)
                    bf(2,bfx(j)) = indx(3-j)
                    bf(3,bfx(j)) = nbr(1,j)
                  else
                    bf(1,bfx(j)) = nbr(2,j)
                    bf(2,bfx(j)) = nbr(1,j)
                    bf(3,bfx(j)) = indx(3-j)
                  end if

                  if(bf(1,k) == indx(j)) then
                    bf(1,k) = indx(3-j)
                  else if(bf(2,k) == indx(j)) then
                    bf(2,k) = indx(3-j)
                  else
                    bf(3,k) = indx(3-j)
                  end if

                end do

                if(msglvl == 4) then
                  write(*,640) a,b,d,e
                end if
  
              end if

            else if(f == g) then

              call updatf(a,c,d,b,e,i,npt,sizht,front,back,fc,ht,ierr)
              call updatf(a,c,e,b,d,i,npt,sizht,front,back,fc,ht,ierr)
              call updatf(b,c,d,a,e,i,npt,sizht,front,back,fc,ht,ierr)
              call updatf(b,c,e,a,d,i,npt,sizht,front,back,fc,ht,ierr)
              call updatf(a,d,f,b,e,i,npt,sizht,front,back,fc,ht,ierr)
              call updatf(a,e,f,b,d,i,npt,sizht,front,back,fc,ht,ierr)
              call updatf(b,d,f,a,e,i,npt,sizht,front,back,fc,ht,ierr)
              call updatf(b,e,f,a,d,i,npt,sizht,front,back,fc,ht,ierr)

              if(ierr /= 0) then
                write(*,'(a)') ' '
                write(*,'(a)') 'SWAPES - Fatal error!'
                write(*,'(a,i8)') '  UPDATF returned IERR = ',ierr
                return
              end if

              call htdel(ind,npt,sizht,fc,ht)
              call htins(ind,c,d,e,a,b,npt,sizht,fc,ht)

              ind = htsrc(a,b,f,npt,sizht,fc,ht)

              if(ind <= 0) then
                write(*,'(a)') ' '
                write(*,'(a)') 'SWAPES - Fatal error!'
                write(*,'(a)') '  HTSRC returned IND <= 0.'
                ierr = 300
                return
              end if

              call htdel(ind,npt,sizht,fc,ht)

              if(0 <= fc(7,ind)) then
                fc(2,ind) = 0
                call availf(hdavfc,fc_num,fc_max,fc,ind,ierr)
                if(ierr /= 0) then
                  write(*,'(a)') ' '
                  write(*,'(a)') 'SWAPES - Fatal error!'
                  write(*,'(a,i8)') '  AVAILF returned IERR = ',ierr
                  return
                end if
              end if

              call htins(ind,d,e,f,a,b,npt,sizht,fc,ht)
              call htdel(ind1,npt,sizht,fc,ht)
              j = fc(7,ind1)
              call htins(ind1,a,d,e,c,f,npt,sizht,fc,ht)
              fc(7,ind1) = j
              call htdel(ind2,npt,sizht,fc,ht)
              j = fc(7,ind2)
              call htins(ind2,b,d,e,c,f,npt,sizht,fc,ht)
              fc(7,ind2) = j

              if(i <= 0 .and. fc(7,ind1) == -1) then
                fc(7,ind1) = 0
                if(front == 0) then
                  front = ind1
                else
                  fc(7,back) = ind1
                end if
                back = ind1
              end if

              if(i <= 0 .and. fc(7,ind2) == -1) then
                fc(7,ind2) = 0
                if(front == 0) then
                  front = ind2
                else
                  fc(7,back) = ind2
                end if
                back = ind2
              end if

              if(msglvl == 4) then
                write(*,650) a,b,d,e,f
              end if

            else

              if(msglvl == 4) then
                write(*,660) a,b,d,e,f,g
              end if

            end if

          end if

        end if

      end do

      600 format(1x,'index =',i7,' : ',5i7)
      610 format(4x,'swap 2-3')
      620 format(4x,'swap 3-2 with new common face:',3i7)
      630 format(4x,'swap 3-2 not poss,tetra missing:',4i7)
      640 format(4x,'swap 2-2: edge ',2i7,' repl by ',2i7)
      650 format(4x,'swap 4-4: edge ',2i7,' repl by ',2i7,'   f =',i7)
      660 format(4x,'swap 4-4 not poss: a,b,d,e,f,g =',6i7)

      return
    end subroutine
    subroutine tetlst(fc_max,fc_num,vm,fc,tetra_num,tetra_num2,tetra)

    !*****************************************************************************80
    !
    !! TETLST constructs a list of tetrahedra from the FC array.
    !
    !  Discussion: 
    !
    !    This routine constructs a list of tetrahedra from the FC array. 
    !
    !    Global vertex indices from VM are produced.  The vertex indices for each
    !    tetrahedron are sorted in increasing order.
    !
    !  Licensing:
    !
    !    This code is distributed under the GNU LGPL license. 
    !
    !  Modified:
    !
    !    05 September 2005
    !
    !  Author:
    !
    !    Original FORTRAN77 version by Barry Joe.
    !    FORTRAN90 version by John Burkardt.
    !
    !  Parameters:
    !
    !    Input,integer :: FC_MAX,the maximum number of positions 
    !    in the FC array.
    !
    !    Input,integer :: FC_NUM,the number of positions used in FC.
    !
    !    Input,integer :: VM(1:*),the indices of vertices of VCL that
    !    are triangulated.
    !
    !    Input,integer :: FC(7,FC_MAX),array of face records; 
    !    see routine DTRIS3.
    !
    !    Input,integer :: TETRA_NUM,the number of tetrahedrons expected.
    !
    !    Output,integer :: TETRA_NUM2,the number of tetrahedrons found.
    !
    !    Output,integer :: TETRA(4,TETRA_NUM),contains global 
    !    tetrahedron indices; it is assumed there is enough space.
    !
      implicit none

      integer :: fc_max
      integer :: tetra_num

      integer :: a
      integer :: fc(7,fc_max)
      integer :: fc_num
      integer :: i
      integer :: j
      integer :: k
      integer :: l
      integer :: t(4)
      integer :: tetra(4,tetra_num)
      integer :: tetra_num2
      integer :: vm(*)

      tetra_num2 = 0

      do i = 1,fc_num

        if(fc(1,i) <= 0)  then
          cycle
        end if

        do k = 4,5

          if(fc(3,i) < fc(k,i)) then

            tetra_num2 = tetra_num2 + 1

            if(tetra_num2 <= tetra_num) then
              tetra(1,tetra_num2) = fc(1,i)
              tetra(2,tetra_num2) = fc(2,i)
              tetra(3,tetra_num2) = fc(3,i)
              tetra(4,tetra_num2) = fc(k,i)
            end if

          end if

        end do

      end do

      if(tetra_num2 /= tetra_num) then
        write(*,'(a)') ' '
        write(*,'(a)') 'TETLST - Warning!'
        write(*,'(a)') '  Inconsistent tetrahedron information.'
        write(*,'(a,i8)') '  Was expecting TETRA_NUM = ',tetra_num
        write(*,'(a,i8)') '  Found TETRA_NUM2 = ',tetra_num2
        return
      end if

      do k = 1,tetra_num2

        t(1:4) = vm(tetra(1:4,k))

        do i = 1,3
          l = i
          do j = i+1,4
            if(t(j) < t(l)) then
              l = j
            end if
          end do
          a = t(i)
          t(i) = t(l)
          t(l) = a
        end do

        tetra(1:4,k) = t(1:4)

      end do

      return
    end subroutine
    subroutine updatf(a,b,c,d,e,i,n,p,front,back,fc,ht,ierr)

    !*****************************************************************************80
    !
    !! UPDATF updates a record in FC after a local transformation.
    !
    !  Discussion: 
    !
    !    This routine updates a record in FC due to a local transformation.
    !
    !    Tetrahedron ABCD becomes ABCE. Add face ABC to queue if it is
    !    interior face,not yet in queue,and its largest index isn't I.
    !
    !  Licensing:
    !
    !    This code is distributed under the GNU LGPL license. 
    !
    !  Modified:
    !
    !    24 January 2009
    !
    !  Author:
    !
    !    Original FORTRAN77 version by Barry Joe.
    !    FORTRAN90 version by John Burkardt.
    !
    !  Parameters:
    !
    !    Input,integer :: A,B,C,the first 3 fields of FC record 
    !   (in any order).
    !
    !    Input,integer :: D,E,the fourth vertex indices of old and 
    !    new tetrahedrons.
    !
    !    Input,integer :: I,the vertex index determining whether 
    !    face put on queue.
    !
    !    Input,integer :: N,the upper bound on vertex indices.
    !
    !    Input,integer :: P,the size of hash table.
    !
    !    Input/output,integer :: FRONT,BACK,the front and back 
    !    pointers of queue.
    !
    !    Input,integer :: FC(1:7,1:*),the array of face records; 
    !    see routine DTRIS3.
    !
    !    Input,integer :: HT(0:P-1),the hash table using direct 
    !    chaining.
    !
    !    Output,integer :: IERR,error flag,which is zero unless 
    !    an error occurred.
    !
      implicit none

      integer :: p

      integer :: a
      integer :: b
      integer :: back
      integer :: c
      integer :: d
      integer :: e
      integer :: fc(7,*)
      integer :: front
      integer :: ht(0:p-1)
   !   integer :: htsrc
      integer :: i
      integer :: ierr
      integer :: ind
      integer :: n

      ierr = 0
      ind = htsrc(a,b,c,n,p,fc,ht)

      if(ind <= 0) then
        ierr = 300
        return
      end if

      if(fc(4,ind) == d) then
        fc(4,ind) = e
      else
        fc(5,ind) = e
      end if

      if(fc(7,ind) == -1 .and. fc(3,ind) /= i .and. 0 < fc(5,ind)) then

        fc(7,ind) = 0

        if(front == 0) then
          front = ind
        else
          fc(7,back) = ind
        end if

        back = ind
      end if

      return
    end subroutine
    subroutine vbfac(pt,ctr,vcl,vm,bf,fc,topv,topnv)

    !*****************************************************************************80
    !
    !! VBFAC determines the boundary faces of a 3D triangulation.
    !
    !  Discussion: 
    !
    !    This routine determines boundary faces of a 3D triangulation visible
    !    from point PT,given a starting visible boundary face.
    !
    !  Licensing:
    !
    !    This code is distributed under the GNU LGPL license. 
    !
    !  Modified:
    !
    !    24 January 2009
    !
    !  Author:
    !
    !    Original FORTRAN77 version by Barry Joe.
    !    FORTRAN90 version by John Burkardt.
    !
    !  Parameters:
    !
    !    Input,real(dp) :: PT(1:3),the 3D point.
    !
    !    Input,real(dp) :: CTR(1:3),the 3D point in interior of
    !    triangulation.
    !
    !    Input,real(dp) :: VCL(1:3,1:*),the vertex coordinate list.
    !
    !    Input,integer :: VM(1:*),the indices of vertices of VCL 
    !    being triangulated.
    !
    !    Input,integer :: BF(1:3,1:*),the array of boundary face 
    !    records; see DTRIS3.
    !
    !    Input/output,integer :: FC(1:7,1:*),the array of face records;
    !    see routine DTRIS3; row 7 is used for links of 3 stacks in this routine.  
    !    On output,FC(7,*) has been updated,so that only stack of visible 
    !    boundary faces remains.
    !
    !    Input/output,integer :: TOPV.  On input,index of FC of visible
    !    boundary face.  On output,index of top of stack of visible boundary faces.
    !
    !    Input,integer :: TOPNV,the index of top of stack of boundary 
    !    faces already found to be not visible from PT,or 0 for empty stack.
    !
    !    Output,integer :: IERR,error flag,which is zero unless an 
    !    error occurred.
    !
      implicit none

      integer :: bf(3,*)
      real(dp) :: ctr(3)
      integer :: fc(7,*)
      integer :: ierr
      integer :: j
      integer :: k
      integer :: nbr
      integer :: op
 !     integer :: opside
      real(dp) :: pt(3)
      integer :: ptr
      integer :: topn
      integer :: topnv
      integer :: topt
      integer :: topv
      integer :: va
      integer :: vb
      integer :: vc
      integer :: vm(*)
      real(dp) :: vcl(3,*)
    !
    !  TOPN is index of top of stack of non-visible boundary faces.
    !  TOPT is index of top of stack of boundary faces to be tested.
    !
      topn = topnv
      topt = 0
      fc(7,topv) = 0
      k = -fc(5,topv)

      do j = 1,3
        nbr = bf(j,k)
        if(fc(7,nbr) == -1) then
          fc(7,nbr) = topt
          topt = nbr
        end if
      end do

      do while(topt /= 0)

        ptr = topt
        topt = fc(7,ptr)
        va = vm(fc(1,ptr))
        vb = vm(fc(2,ptr))
        vc = vm(fc(3,ptr))
        op = opside(vcl(1,va),vcl(1,vb),vcl(1,vc),ctr,pt)

        if(op == 2) then
          ierr = 301
          return
        end if

        if(op == 1) then

          fc(7,ptr) = topv
          topv = ptr
          k = -fc(5,ptr)

          do j = 1,3
            nbr = bf(j,k)
            if(fc(7,nbr) == -1) then
              fc(7,nbr) = topt
              topt = nbr
            end if
          end do

        else

          fc(7,ptr) = topn
          topn = ptr

        end if

      end do
    !
    !  For boundary faces not visible from PT,set FC(7,*) = -1.
    !
      do while(topn /= 0) 
        ptr = topn
        topn = fc(7,ptr)
        fc(7,ptr) = -1
      end do

      return
    end subroutine
    !
    ! for volume
    !
    pure function tetrahedron_volume(tetra) result(volume)
    !
    ! TETRAHEDRON_VOLUME_3D computes the volume of a tetrahedron.
    !
    !  Licensing:
    !
    !    This code is distributed under the GNU LGPL license.
    !
    !  Modified:
    !
    !    30 December 2004
    !
    !  Author:
    !
    !    John Burkardt
    !
    !  Parameters:
    !
    !    Input, real(dp) :: TETRA(3,4), the vertices of the tetrahedron.
    !
    !    Output, real(dp) :: VOLUME, the volume of the tetrahedron.
    !
    !
    implicit none
    !
    integer, parameter :: dim_num = 3
    real(dp) :: a(4,4)
    real(dp) :: det
    real(dp), intent(in) :: tetra(dim_num,4)
    real(dp) :: volume
    !
    a(1:dim_num,1:4) = tetra(1:dim_num,1:4)
    a(4,1:4) = 1.0_dp
    !
    ! get determinant
    !
    det = &
        a(1,1) * ( &
        a(2,2) * ( a(3,3) * a(4,4) - a(3,4) * a(4,3) ) &
        - a(2,3) * ( a(3,2) * a(4,4) - a(3,4) * a(4,2) ) &
        + a(2,4) * ( a(3,2) * a(4,3) - a(3,3) * a(4,2) ) ) &
        - a(1,2) * ( &
        a(2,1) * ( a(3,3) * a(4,4) - a(3,4) * a(4,3) ) &
        - a(2,3) * ( a(3,1) * a(4,4) - a(3,4) * a(4,1) ) &
        + a(2,4) * ( a(3,1) * a(4,3) - a(3,3) * a(4,1) ) ) &
        + a(1,3) * ( &
        a(2,1) * ( a(3,2) * a(4,4) - a(3,4) * a(4,2) ) &
        - a(2,2) * ( a(3,1) * a(4,4) - a(3,4) * a(4,1) ) &
        + a(2,4) * ( a(3,1) * a(4,2) - a(3,2) * a(4,1) ) ) &
        - a(1,4) * ( &
        a(2,1) * ( a(3,2) * a(4,3) - a(3,3) * a(4,2) ) &
        - a(2,2) * ( a(3,1) * a(4,3) - a(3,3) * a(4,1) ) &
        + a(2,3) * ( a(3,1) * a(4,2) - a(3,2) * a(4,1) ) )
    !
    ! volume
    !
    volume = abs(det)/6.0_dp
    !
    end function tetrahedron_volume
   
end module