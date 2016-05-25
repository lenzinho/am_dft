module am_character_table

	implicit none

	private


	type :: am_class_chartab_rep
		integer :: nreps
		character(10) :: flags ! rep/irrep/orb
		complex(dp) , allocatable :: chi(:,:)
		character(:), allocatable :: label(:) 
	end type am_class_chartab_rep

	type, public :: am_class_chartab
		integer :: nreps
		type(am_class_chartab_rep), allocatable :: rep(:)
		contains
           ! procedure :: add_subgroup
           ! procedure :: add_orbitals
           ! procedure :: add_kpt_path
	end type am_class_chartab

	contains

	subroutine 	   initialize_character_table(ct,pg)
		!
		implicit none
		!
		class(am_class_chartab), intent(out) :: ct
        class(am_class_group)  , intent(in) :: grp
        !

        ! initialize rep counter
        ct%nreps = 0
        ! add orbital rep
		ct%nreps = ct%nreps + 1
        allocate(ct%rep(1))
        call get_orbital_rep(ct%rep(1),grp)
        !
    end subroutine initialize_character_table

    subroutine     add_irrep(ct,pg,prepend_label)
		!
		implicit none
		!
		class(am_class_chartab) , intent(inout) :: ct
        class(am_class_group)   , intent(in) :: grp
        character(*) 		    , intent(in) :: prepend_label
        class(am_class_chartab_rep) :: rep
		!
		! copy rep
		rep = ct%rep
		! increase rep size
		ct%nreps = ct%nreps + 1 
		! deallocate
		deallocate(ct%rep)
		! reallocate
		allocate(ct%rep(ct%nreps))
		! copy back
		do i = 1, ct%nreps - 1
			ct%rep(i) = rep(i)
		enddo
		!
		! add irrep
		ct%rep(ct%nreps)%flags = 'irrep'
		! add irrep character table
        ct%rep(ct%nreps)%chi = get_chartab(multab=grp%mt%multab, nclasses=grp%cc%nclasses, &
            class_nelements=grp%cc%nelements, class_member=grp%cc%member, ps_id=grp%ps_id)
        ! get number of irreps
        ct%rep(ct%nreps)%nreps = grp%cc%nclasses
        ! get muliken labels for irreps
        ct%rep(ct%nreps)%label = get_muliken(chartab=grp%ct%chartab, nclasses=grp%cc%nclasses, &
            class_nelements=grp%cc%nelements, class_member=grp%cc%member, ps_id=grp%ps_id)
        !
		call get_orbital_rep(ct%rep(ct%nreps), grp, prepend_label)
		!
        contains
		subroutine 	   add_orbital_rep(rep,grp)
			!
			implicit none
			!
			class(am_class_chartab_rep), intent(out) :: rep
	        class(am_class_group)      , intent(in)  :: grp
			integer :: norbitals
			!
	        rep%nreps = 4
	        rep%flags = 'orb'
	        allocate(rep%chi(rep%nreps,grp%cc%nclasses))
	        allocate(rep%label(rep%nreps))
	        ! representation based on s orbitals
	        rep%label(1) = 's D^(0)'
	        rep%chi(1,:) = get_orb_characters(l=0, R=grp%sym(1:3,1:3,:), nclasses=grp%cc%nclasses, class_representative=grp%cc%representative)
	        ! representation based on p orbitals
	        rep%label(2) = 'p D^(1)'
	        rep%chi(2,:) = get_orb_characters(l=1, R=grp%sym(1:3,1:3,:), nclasses=grp%cc%nclasses, class_representative=grp%cc%representative)
	        ! representation based on d orbitals
	        rep%label(3) = 'd D^(2)'
	        rep%chi(3,:) = get_orb_characters(l=2, R=grp%sym(1:3,1:3,:), nclasses=grp%cc%nclasses, class_representative=grp%cc%representative)
	        ! representation based on f orbitals
	        rep%label(4) = 'f D^(3)'
	        rep%chi(4,:) = get_orb_characters(l=3, R=grp%sym(1:3,1:3,:), nclasses=grp%cc%nclasses, class_representative=grp%cc%representative)
	        !
	    end subroutine add_orbital_rep
        function       get_orb_characters(l,R,nclasses,class_representative) result(orbchi)
            !
            implicit none
            !
            integer,  intent(in) :: l
            real(dp), intent(in) :: R(:,:,:)
            integer,  intent(in) :: nclasses
            integer,  intent(in) :: class_representative(:)
            complex(dp), allocatable :: orbchi(:)
            integer :: i
            !
            allocate(orbchi(nclasses))
            !
            do i = 1, nclasses
                orbchi(i) = get_orbchi(l=l,R=R(:,:,class_representative(i)))
            enddo
            !
        end function   get_orb_characters
        function       get_orbchi(l,R) result(chi)
            !
            ! Character of (2l+1)-dimensional orthogonal group O(3) irrep parameterized by 3x3 rotoinversion matrix R
            !       
            ! T. Wolfram and Ş. Ellialtıoğlu, Applications of Group Theory to Atoms, Molecules, 
            ! and Solids, 1 edition (Cambridge University Press, Cambridge, 2014), p 74, Eq. 3.20.
            !
            ! M. Tinkham, Group Theory and Quantum Mechanics (McGraw-Hill, 1964), p 66 Eq 4.6 
            ! and p 100 bottom.
            !
            implicit none
            !
            integer,  intent(in) :: l ! dimension of irrep
            real(dp), intent(in) :: R(3,3) ! rotation angle
            real(dp) :: aa(4)
            real(dp) :: d
            complex(dp) :: chi
            !
            ! if R is a rotoinversion, get the angle and axis of the rotational part only (without the inversion)
            d = R(1,1)*R(2,2)*R(3,3)-R(1,1)*R(2,3)*R(3,2)-R(1,2)*R(2,1)*R(3,3)+R(1,2)*R(2,3)*R(3,1)+R(1,3)*R(2,1)*R(3,2)-R(1,3)*R(2,2)*R(3,1)
            ! convert to proper rotation
            aa = rot2axis_angle(d*R)
            if (abs(aa(4)).lt.tiny) then
                ! limiting case is obtained by expanding sin( (j+1/2)*x ) / sin( x/2 ) at zero
                chi = (2*l+1)
            else
                ! general case
                chi = sin( (l+0.5_dp)*aa(4) ) / sin( aa(4)*0.5_dp )
            endif
            ! convert back to improper rotation
            chi = chi * sign(1.0_dp,d)**(l)
            !
        end function   get_orbchi
	end subroutine add_irrep

	subroutine 	   add_rep(ct,grp,prepend_label)
		!
		implicit none
		!
		class(am_class_chartab) , intent(inout) :: ct
        class(am_class_group)   , intent(in) :: grp
        character(*) 		    , intent(in) :: prepend_label
        class(am_class_chartab_rep) :: rep
		!
		! copy rep
		rep = ct%rep
		! increase rep size
		ct%nreps = ct%nreps + 1 
		! deallocate
		deallocate(ct%rep)
		! reallocate
		allocate(ct%rep(ct%nreps))
		! copy back
		do i = 1, ct%nreps - 1
			ct%rep(i) = rep(i)
		enddo
		! add new rep
		call get_orbital_rep(ct%rep(ct%nreps), grp, prepend_label)
		!
		contains
		subroutine 	   get_orbital_rep(rep, grp, prepend_label)
			!
			implicit none
			!
			class(am_class_chartab_rep), intent(out) :: rep
	        class(am_class_group)      , intent(in)  :: grp
	        character(*) 		       , intent(in)  :: prepend_label
			integer :: norbitals
			!
	        rep%nreps = grp%nsyms
	        rep%flags = 'rep'
	        allocate(rep%chi(rep%nreps,grp%cc%nclasses))
	        allocate(rep%label(rep%nreps))
	        !
			do i = 1, grp%nsyms
				rep%label(i) = prepend_label//trim(int2char(i))
				rep%chi(i,:) = get_rep_characters(sym=grp%sym, nclasses=grp%cc%nclasses, class_representative=grp%cc%representative)
			enddo
	    end subroutine get_orbital_rep
    end subroutine add_rep

    subroutine     print_master_character_table(master_chartab)
        !
        implicit none
        !
        class(am_class_chartab), intent(in) :: master_chartab
        integer :: i
        integer :: nreps
        complex(dp) , allocatable :: rep_chi(:,:)
        character(7), allocatable :: rep_label(:)
        !
        ! start rep counter
        i=0
        !
        select type (grp)
        class is (am_class_point_group) 
            nreps = 5
            allocate(rep_chi(nreps,grp%cc%nclasses))
            allocate(rep_label(nreps))
            ! representation based on basis functions
            i=i+1; rep_label(i) = 'rep'
            rep_chi(i,:) = get_rep_characters(sym=grp%sym(1:3,1:3,:), nclasses=grp%cc%nclasses, class_representative=grp%cc%representative)
            ! representation based on s orbitals
            i=i+1; rep_label(i) = 's D^(0)'
            rep_chi(i,:) = get_orb_characters(l=0, R=grp%sym(1:3,1:3,:), nclasses=grp%cc%nclasses, class_representative=grp%cc%representative)
            ! representation based on p orbitals
            i=i+1; rep_label(i) = 'p D^(1)'
            rep_chi(i,:) = get_orb_characters(l=1, R=grp%sym(1:3,1:3,:), nclasses=grp%cc%nclasses, class_representative=grp%cc%representative)
            ! representation based on d orbitals
            i=i+1; rep_label(i) = 'd D^(2)'
            rep_chi(i,:) = get_orb_characters(l=2, R=grp%sym(1:3,1:3,:), nclasses=grp%cc%nclasses, class_representative=grp%cc%representative)
            ! representation based on f orbitals
            i=i+1; rep_label(i) = 'f D^(3)'
            rep_chi(i,:) = get_orb_characters(l=3, R=grp%sym(1:3,1:3,:), nclasses=grp%cc%nclasses, class_representative=grp%cc%representative)
        class is (am_class_space_group) 
            ! do nothing.
        class default
            nreps = 1
            allocate(rep_chi(nreps,grp%cc%nclasses))
            allocate(rep_label(nreps))
            ! representation based on basis functions
            i=i+1; rep_label(i) = 'rep'
            rep_chi(i,:) = get_rep_characters(sym=grp%sym, nclasses=grp%cc%nclasses, class_representative=grp%cc%representative)
        end select
        !
        !
        if (allocated(rep_chi).and.allocated(rep_label)) then
            call print_chartab(chartab=grp%ct%chartab, nclasses=grp%cc%nclasses, class_nelements=grp%cc%nelements, &
                class_member=grp%cc%member, irrep_label=grp%ct%irrep_label, ps_id=grp%ps_id, rep_chi=rep_chi, rep_label=rep_label)
        else
            call print_chartab(chartab=grp%ct%chartab, nclasses=grp%cc%nclasses, class_nelements=grp%cc%nelements, &
                class_member=grp%cc%member, irrep_label=grp%ct%irrep_label, ps_id=grp%ps_id)
        endif
        !
        contains
        function       get_rep_characters(sym,nclasses,class_representative) result(repchi)
            !
            implicit none
            !
            real(dp),    intent(in) :: sym(:,:,:)
            integer,     intent(in) :: nclasses
            integer,     intent(in) :: class_representative(:)
            integer, allocatable :: repchi(:)
            integer :: i
            !
            allocate(repchi(nclasses))
            !
            do i = 1, nclasses
                repchi(i) = trace( sym(:,:,class_representative(i)) )
            enddo
            !
        end function   get_rep_characters
    end subroutine print_master_character_table


    subroutine     get_character_table(grp)
        !
        implicit none
        !
        class(am_class_group), intent(inout) :: grp
        !
        ! get character table (ps_id is required for sorting the irreps: proper before improper)
        grp%ct%chartab = get_chartab(multab=grp%mt%multab, nclasses=grp%cc%nclasses, &
            class_nelements=grp%cc%nelements, class_member=grp%cc%member, ps_id=grp%ps_id)
        ! get muliken labels for irreps
        grp%ct%irrep_label = get_muliken(chartab=grp%ct%chartab, nclasses=grp%cc%nclasses, &
            class_nelements=grp%cc%nelements, class_member=grp%cc%member, ps_id=grp%ps_id)
        !
    end subroutine get_character_table

    subroutine     print_character_table(grp)
        !
        implicit none
        !
        class(am_class_group), intent(in) :: grp
        integer :: i
        integer :: nreps
        complex(dp) , allocatable :: rep_chi(:,:)
        character(7), allocatable :: rep_label(:)
        !
        ! start rep counter
        i=0
        !
        select type (grp)
        class is (am_class_point_group) 
            nreps = 5
            allocate(rep_chi(nreps,grp%cc%nclasses))
            allocate(rep_label(nreps))
            ! representation based on basis functions
            i=i+1; rep_label(i) = 'rep'
            rep_chi(i,:) = get_rep_characters(sym=grp%sym(1:3,1:3,:), nclasses=grp%cc%nclasses, class_representative=grp%cc%representative)
            ! representation based on s orbitals
            i=i+1; rep_label(i) = 's D^(0)'
            rep_chi(i,:) = get_orb_characters(l=0, R=grp%sym(1:3,1:3,:), nclasses=grp%cc%nclasses, class_representative=grp%cc%representative)
            ! representation based on p orbitals
            i=i+1; rep_label(i) = 'p D^(1)'
            rep_chi(i,:) = get_orb_characters(l=1, R=grp%sym(1:3,1:3,:), nclasses=grp%cc%nclasses, class_representative=grp%cc%representative)
            ! representation based on d orbitals
            i=i+1; rep_label(i) = 'd D^(2)'
            rep_chi(i,:) = get_orb_characters(l=2, R=grp%sym(1:3,1:3,:), nclasses=grp%cc%nclasses, class_representative=grp%cc%representative)
            ! representation based on f orbitals
            i=i+1; rep_label(i) = 'f D^(3)'
            rep_chi(i,:) = get_orb_characters(l=3, R=grp%sym(1:3,1:3,:), nclasses=grp%cc%nclasses, class_representative=grp%cc%representative)
        class is (am_class_space_group) 
            ! do nothing.
        class default
            nreps = 1
            allocate(rep_chi(nreps,grp%cc%nclasses))
            allocate(rep_label(nreps))
            ! representation based on basis functions
            i=i+1; rep_label(i) = 'rep'
            rep_chi(i,:) = get_rep_characters(sym=grp%sym, nclasses=grp%cc%nclasses, class_representative=grp%cc%representative)
        end select
        !
        !
        if (allocated(rep_chi).and.allocated(rep_label)) then
            call print_chartab(chartab=grp%ct%chartab, nclasses=grp%cc%nclasses, class_nelements=grp%cc%nelements, &
                class_member=grp%cc%member, irrep_label=grp%ct%irrep_label, ps_id=grp%ps_id, rep_chi=rep_chi, rep_label=rep_label)
        else
            call print_chartab(chartab=grp%ct%chartab, nclasses=grp%cc%nclasses, class_nelements=grp%cc%nelements, &
                class_member=grp%cc%member, irrep_label=grp%ct%irrep_label, ps_id=grp%ps_id)
        endif
        !
        contains
        function       get_rep_characters(sym,nclasses,class_representative) result(repchi)
            !
            implicit none
            !
            real(dp),    intent(in) :: sym(:,:,:)
            integer,     intent(in) :: nclasses
            integer,     intent(in) :: class_representative(:)
            integer, allocatable :: repchi(:)
            integer :: i
            !
            allocate(repchi(nclasses))
            !
            do i = 1, nclasses
                repchi(i) = trace( sym(:,:,class_representative(i)) )
            enddo
            !
        end function   get_rep_characters
        function       get_orb_characters(l,R,nclasses,class_representative) result(orbchi)
            !
            implicit none
            !
            integer,  intent(in) :: l
            real(dp), intent(in) :: R(:,:,:)
            integer,  intent(in) :: nclasses
            integer,  intent(in) :: class_representative(:)
            complex(dp), allocatable :: orbchi(:)
            integer :: i
            !
            allocate(orbchi(nclasses))
            !
            do i = 1, nclasses
                orbchi(i) = get_orbchi(l=l,R=R(:,:,class_representative(i)))
            enddo
            !
        end function   get_orb_characters
        function       get_orbchi(l,R) result(chi)
            !
            ! Character of (2l+1)-dimensional orthogonal group O(3) irrep parameterized by 3x3 rotoinversion matrix R
            !       
            ! T. Wolfram and Ş. Ellialtıoğlu, Applications of Group Theory to Atoms, Molecules, 
            ! and Solids, 1 edition (Cambridge University Press, Cambridge, 2014), p 74, Eq. 3.20.
            !
            ! M. Tinkham, Group Theory and Quantum Mechanics (McGraw-Hill, 1964), p 66 Eq 4.6 
            ! and p 100 bottom.
            !
            implicit none
            !
            integer,  intent(in) :: l ! dimension of irrep
            real(dp), intent(in) :: R(3,3) ! rotation angle
            real(dp) :: aa(4)
            real(dp) :: d
            complex(dp) :: chi
            !
            ! if R is a rotoinversion, get the angle and axis of the rotational part only (without the inversion)
            d = R(1,1)*R(2,2)*R(3,3)-R(1,1)*R(2,3)*R(3,2)-R(1,2)*R(2,1)*R(3,3)+R(1,2)*R(2,3)*R(3,1)+R(1,3)*R(2,1)*R(3,2)-R(1,3)*R(2,2)*R(3,1)
            ! convert to proper rotation
            aa = rot2axis_angle(d*R)
            if (abs(aa(4)).lt.tiny) then
                ! limiting case is obtained by expanding sin( (j+1/2)*x ) / sin( x/2 ) at zero
                chi = (2*l+1)
            else
                ! general case
                chi = sin( (l+0.5_dp)*aa(4) ) / sin( aa(4)*0.5_dp )
            endif
            ! convert back to improper rotation
            chi = chi * sign(1.0_dp,d)**(l)
            !
        end function   get_orbchi
    end subroutine print_character_table

    ! low level stuff

    function       get_chartab(multab,nclasses,class_nelements,class_member,ps_id) result(chartab)
        !
	    use am_rank_and_sort
        !
        implicit none
        !
        integer, intent(in) :: multab(:,:)
        integer, intent(in) :: ps_id(:)
        integer, intent(in) :: nclasses
		integer, intent(in) :: class_nelements(:) ! number of elements in each class
		integer, intent(in) :: class_member(:,:) ! members(nclass,maxval(class_nelements))
        complex(dp), allocatable :: chartab(:,:) ! class constant chartab which becomes character chartab (can be complex!)
        integer :: i, j, k
        integer :: nsyms
        integer :: nirreps  ! just to make things clearer; always equal to nclasses
        integer, allocatable :: H(:,:,:) ! class multiplication 
        integer, allocatable :: irrep_dim(:) ! sym dimensions
        integer, allocatable :: indices(:)
        real(dp) :: wrk
        !
        nsyms = size(multab,2)
        !
        ! get class coefficients (thenumber of times class k appears in the pdocut o class j and k )
        ! Wooten p 35, Eq 2.10; p 40, Example 2.8; p 89, Example 4.6
        allocate(H(nclasses,nclasses,nclasses))
        do i = 1,nclasses
        do j = 1,nclasses
        do k = 1,nclasses
            H(j,k,i) = count( pack( multab(class_member(i,1:class_nelements(i)),class_member(j,1:class_nelements(j))) , .true. ) .eq. class_member(k,1) )
        enddo
        enddo
        enddo
        !
        ! chartab(nirreps,nclasses) determine eigenvector which simultaneously diagonalizes i Hjk(:,:) matrices
        chartab = eigen_analysis(cmplx(H,0,dp))
        !
        ! get dimensions of representations using Eq. 4.53, p 91 of Wooten
        ! Note: The absolute square is the result of the conjugate multiplication
        nirreps = nclasses
        allocate(irrep_dim(nirreps))
        !
        do i = 1, nirreps
            ! sum on classes
            wrk = 0
            do j = 1, nclasses
                wrk = wrk + abs(chartab(i,j))**2/real(class_nelements(j),dp)
            enddo
            ! compute irrep dimension
            irrep_dim(i) = nint((nsyms/wrk)**0.5_dp)
        enddo
        !
        ! convert class constant chartab into character chartab
        do i = 1, nirreps
        do j = 1, nclasses
            chartab(i,j) = irrep_dim(i)/real(class_nelements(j),dp)*chartab(i,j)
        enddo
        enddo
        !
        ! initialize indices sort sorting
        allocate(indices(nirreps))
        !
        ! sort irreps based on character of class containing identity
        ! find class containing inversion (1 is the ps_id for identity)
        id_search : do i = 1, nclasses
        do j = 1, class_nelements(i)
            if (ps_id(class_member(i,j)).eq.1) then
                ! sort based on inversion-containing class character
                call rank(nint(real(chartab(:,i))),indices)
                ! indices = indices(nirreps:1:-1)
                chartab = chartab(indices,:)
                ! exit loop
                exit id_search
            endif
        enddo
        enddo id_search
        !
        ! sort irreps based on character of class containing inversion
        ! find class containing inversion (6 is the ps_id for inversion)
        inv_search : do i = 1, nclasses
        do j = 1, class_nelements(i)
            if (ps_id(class_member(i,j)).eq.6) then
                ! sort based on inversion-containing class character
                call rank(-sign(1,nint(real(chartab(:,i)))),indices)
                ! indices = indices(nirreps:1:-1)
                chartab = chartab(indices,:)
                ! exit loop
                exit inv_search
            endif
        enddo
        enddo inv_search
        !
        !
        ! check that the number of elements ri in class Ci is a divisor of theorder of the group
        ! Symmetry and Condensed Matter Physics: A Computational Approach. 1 edition. Cambridge, UK?; New York: Cambridge University Press, 2008. p 35.
        do i = 1, nclasses
            if (modulo(nsyms,class_nelements(i)).ne.0) then
                stop 'Number of elements in class is not a divisor of the group order.'
            endif
        enddo
        !
        contains
        function       eigen_analysis(A) result(D)
            !
            ! determines characters which simultaneously diagonalizes i square matrices A(:,:,i)
            ! only possible if A(:,:,i) commute with each other. Assumes all dimensions of A are the same.
            ! 
            !
            implicit none
            !
            complex(dp), intent(in)  :: A(:,:,:)
            complex(dp), allocatable :: V(:,:) ! these two should match in this particular situation (character chartab determination from class multiplication )
            complex(dp), allocatable :: D(:,:) ! these two should match in this particular situation (character chartab determination from class multiplication )
            complex(dp), allocatable :: M(:,:)
            integer :: n, i
            !
            n = size(A,3)
            !
            ! matrix pencil approach to lift the degeneracy of eigenvalues
            allocate(M(n,n))
            M = 0
            do i = 1,n
                M = M + rand() * A(:,:,i)
            enddo
            !
            ! get left and right eigenvectors: transpose(VL)*A*VR = D
            call am_zgeev(A=M,VR=V)
            !
            ! check that matrices have been diagonalized properly and save eigenvalues
            allocate(d(n,n))
            do i = 1, n
                !
                M = matmul(inv(V),matmul(A(:,:,i),V))
                !
                if (any( abs(diag(diag(M))-M).gt. (tiny*n**2) )) then
                    ! tiny may be too small a criterion here... should probably use a value which scales with the size of the matrix.
                    stop 'Unable to perform simultaneous matrix diagonalization. Check that whether they commute.'
                endif
                !
                ! save diagonal elements as row vectors
                D(:,i) = diag( M )
                !
            enddo
            !
            ! should enforce check here...
            !
            ! ! normalize each eigencolumn to the element of the column (that is, the first element, i.e. row, corresponds to the class containing the identity element)
            ! do i = 1, n
            !     V(:,i)=V(:,i)
            ! enddo
            ! V = transpose(V)
            !
            ! at this point V and D should be identical. The elements are to within +/- sign.
            ! D is more robust 
            !
        end function   eigen_analysis
    end function   get_chartab

    function       get_muliken(chartab,nclasses,class_nelements,class_member,ps_id) result(irrep_label)
        !
        implicit none
        !
        complex(dp), intent(in) :: chartab(:,:) ! class constant chartab which becomes character chartab (can be complex!)
        integer, intent(in) :: nclasses
		integer, intent(in) :: class_nelements(:)
		integer, intent(in) :: class_member(:,:)
        integer, intent(in) :: ps_id(:)
        character(:), allocatable :: irrep_label(:)
        integer :: i, j, k
        integer :: nirreps
        integer :: class_containing_identity
        integer :: class_containing_inversion
        integer :: class_containing_highsym
        integer, allocatable :: highsymlist(:)
        !
        ! they are always identical...
        nirreps = nclasses
        !
        ! create irrep label
        allocate(character(7) :: irrep_label(nirreps))
        !
        class_containing_identity  = find_class_containing_ps(nclasses=nclasses, id=1, ps_id=ps_id, class_nelements=class_nelements, class_member=class_member)
        class_containing_inversion = find_class_containing_ps(nclasses=nclasses, id=6, ps_id=ps_id, class_nelements=class_nelements, class_member=class_member)
        allocate(highsymlist, source=[5,10,4,9,3,8,2,7,1,6])
        do i = 1,10
	        class_containing_highsym = find_class_containing_ps(nclasses=nclasses, id=highsymlist(i), ps_id=ps_id, class_nelements=class_nelements, class_member=class_member)
	        if (class_containing_highsym.ne.0) exit
        enddo
        !
        do j = 1, nirreps
            ! 0) nullify irrep label
            irrep_label(j) = ''
            ! 1) find class containing identity (ps_id=1)
            i = class_containing_identity
            select case (nint(real(chartab(j,i))))
                case (1)
                	k = class_containing_highsym
                	if (real(chartab(j,k)).ge.0) then 
                	 	irrep_label(j) = trim(irrep_label(j))//'A'
                	else
                		irrep_label(j) = trim(irrep_label(j))//'B'
                	endif
                case (2); irrep_label(j) = trim(irrep_label(j))//'E'
                case (3); irrep_label(j) = trim(irrep_label(j))//'T'
            end select
            ! 2) find class containing inversion (ps_id=6)
            i = class_containing_inversion
            if (i.ne.0) then
            select case (sign(1,nint(real(chartab(j,i)))))
                case (+1); irrep_label(j) = trim(irrep_label(j))//'_g'
                case (-1); irrep_label(j) = trim(irrep_label(j))//'_u'
            end select
            endif
            !
        enddo
        !
        contains
            pure function find_class_containing_ps(nclasses, id, ps_id, class_nelements, class_member) result(i)
                !
                implicit none
                !
                integer, intent(in) :: id
                integer, intent(in) :: nclasses
                integer, intent(in) :: ps_id(:)
                integer, intent(in) :: class_nelements(:)
                integer, intent(in) :: class_member(:,:)
                integer :: i, j
                !
                do i = 1, nclasses
                do j = 1, class_nelements(i)
                    if (ps_id(class_member(i,j)).eq.id) then
                        return
                    endif
                enddo
                enddo
                ! if not found
                i = 0
            end function  find_class_containing_ps
    end function   get_muliken

    subroutine     print_chartab(chartab,nclasses,class_nelements,class_member,irrep_label,ps_id,rep_chi,rep_label)
        !
        implicit none
        !
        complex(dp) , intent(in) :: chartab(:,:) ! class constant chartab which becomes character chartab (can be complex!)
        integer		, intent(in) :: nclasses
		integer		, intent(in) :: class_nelements(:)
		integer		, intent(in) :: class_member(:,:)
        character(*), intent(in) :: irrep_label(:)
        integer		, intent(in) :: ps_id(:)
        complex(dp) , intent(in), optional :: rep_chi(:,:)
        character(*), intent(in), optional :: rep_label(:)
        integer :: i, j, k
        integer :: nirreps
         ! rep to irrep decompostion (used only if rep_chi nd rep_label are present)
        integer :: nreps
        integer :: nsyms
        integer, allocatable :: beta(:,:)
        complex(dp) :: try
        ! complex to symbol
        complex(dp), allocatable :: s(:)
        character(:), allocatable :: str
        integer :: char_start
        ! print formats
        character(50) :: fmts(5)
        !
        ! they are always identical...
        nirreps = nclasses
        !
        ! left-side headers (sum to 10)
        fmts(2) = '(5x,a10)'
        fmts(4) = '(5x,a6,i4)'
        ! chartab (sum to 6)
        fmts(1) = '(i6)'
        fmts(3) = '(a6)'
        fmts(5) = '(f6.2)'
        !
        call print_chartab_header(nclasses,class_nelements,class_member,ps_id,fmts)
        !
        allocate(s(100)) ! value of symbolic output, 30 possible options
        allocate(character(4)::str)
        char_start = 96 ! 97 = a
        k=0
        !
        do i = 1, nirreps
            write(*,'(5x,i2,a8)',advance='no') i, irrep_label(i)
            do j = 1, nclasses
                call print_chartab_symb(Z=chartab(i,j),char_start=char_start,s=s,k=k,fmts=fmts,str=str)
                write(*,fmts(3),advance='no') str
            enddo
            write(*,*)
        enddo
        ! 
        if (present(rep_chi).and.present(rep_label)) then
        	if (size(rep_chi,2).ne.nclasses) stop 'dimensions of rep character does not match number of classes'
        	!
	        ! rep characters
	        call print_bar(fmts=fmts, nclasses=nclasses)
    		nreps = size(rep_chi,1)
	        do i = 1, nreps
	            !
	            write(*,'(5x,i2,a8)',advance='no') i, rep_label(i)
	            do j = 1, nclasses
	                call print_chartab_symb(Z=rep_chi(i,j),char_start=char_start,s=s,k=k,fmts=fmts,str=str)
	                write(*,fmts(3),advance='no') str
	            enddo
	            write(*,*)
	        enddo
	    endif
	    !
        ! rep decompositions
	    if (present(rep_chi).and.present(rep_label)) then
	        call print_bar(fmts=fmts, nclasses=nclasses)
	        !
	        write(*,'(5x,a)') 'Decompositions:'
	        ! initialize
	        allocate(beta(nreps,nirreps))
	        beta = 0
	        nsyms = sum(class_nelements)
	        ! calculate
	        do i = 1, nreps
        	do j = 1, nirreps
        		try = dot(rep_chi(i,:),chartab(j,:)*class_nelements)/real(nsyms,dp)
        		if (abs(aimag(try)).gt.tiny) then
        			stop 'irrep decomposition coefficient is not real'
        		endif
        		if (abs(nint(real(try))-real(try)).gt.tiny) then
        			stop 'irrep decomposition coefficient is not an integer'
        		endif
	        	beta(i,j) = nint(real(try))
        	enddo
        	enddo
        	! print
	        do i = 1, nreps
	        if (any(beta(i,:).ne.0)) then
	        	write(*,'(5x,a9,a)',advance='no') trim(rep_label(i)), ' = '
	        	do j = 1, nirreps
        		if (beta(i,j).ne.0) then
	        		write(*,'(a)',advance='no') trim(int2char(beta(i,j),'SP'))//'*'//trim(irrep_label(j))//'('//trim(int2char(j))//') '
	        	endif
	        	enddo
	        	write(*,*)
	        endif
        	enddo
	    endif
	    !
	    !
	    !
	    !
        !
        !
        if (k.ne.0) then
            write(*,'(5x,a)') 'Definitions:'
            do i = 1, k
                write(*,'(5x)',advance='no')
                write(*,'(a)' ,advance='no') char(char_start+k)//' = '
                write(*,'(a)' ,advance='no') trim(dbl2char(real(s(i)),8))//trim(dbl2charSP(aimag(s(i)),9))//'i'
                write(*,*)
            enddo
            !
        endif
        !
	    contains
	    subroutine     print_chartab_header(nclasses,class_nelements,class_member,ps_id,fmts)
	    	!
	    	implicit none
	    	!
	        integer		, intent(in) :: nclasses
			integer		, intent(in) :: class_nelements(:)
			integer		, intent(in) :: class_member(:,:)
	        integer		, intent(in) :: ps_id(:)
	        character(*), intent(in) :: fmts(:)
	        integer :: i
	    	!
	        ! start printing
	        write(*,fmts(2),advance='no') 'class'
	        do i = 1, nclasses
	            write(*,fmts(1),advance='no') i
	        enddo
	        write(*,*)
	        !
	        write(*,fmts(2),advance='no') 'elements'
	        do i = 1, nclasses
	            write(*,fmts(1),advance='no') class_nelements(i)
	        enddo
	        write(*,*)
	        !
	        write(*,fmts(2),advance='no') 'repr.'
	        do i = 1, nclasses
	            write(*,fmts(3),advance='no') trim(decode_pointsymmetry(ps_id(class_member(i,1))))
	        enddo
	        write(*,*)
	        !
	        call print_bar(fmts=fmts, nclasses=nclasses)
	        !
	    end subroutine print_chartab_header
	    subroutine     print_chartab_symb(Z,char_start,s,k,fmts,str)
	    	!
	    	implicit none
	    	!
	        complex(dp), intent(in) :: Z
	        integer    , intent(in) :: char_start
	        complex(dp), intent(inout) :: s(:)
	        integer    , intent(inout) :: k
	        character(*), intent(in) :: fmts(:)
	        character(:), allocatable, intent(out) :: str
	        !
	        integer :: kk, k_exp
	        complex(dp) :: s_exp
	        real(dp) :: Zr,Zi
	        logical :: strmatch
	        !
		    strmatch = .false.
		    !
		    Zr= real(Z)
		    Zi= aimag(Z)
		    !
		    if ( isint(Zr) .and. iszero(Zi) ) then
		        ! no imaginary, integer real
		        strmatch = .true.
		        str = trim(int2char(nint(Zr)))
		        !
		    elseif ( iszero(Zr) .and. isint(Zi) ) then
		        ! no real, imaginary integer
		        strmatch = .true.
		        str = trim(int2char(nint(Zi)))//'i'
		        !
		    elseif ( (.not. isint(Zr)) .and. (.not. isint(Zi)) ) then
		        ! complex number
		        !
		        if (k.ge.1) then
		        search : do kk = 1, k
		            !
		            s_exp = cmplx(0.0_dp,0.0_dp)
		            k_exp = 0 
		            do while ( .not. isequal(s_exp,cmplx(1,0,dp)) )
		                ! do a full loop. complex numbers form a cyclic abelian group.
		                ! exponentiate it until it loops back to one, the identity
		                k_exp = k_exp+1
		                s_exp = s(k)**k_exp
		                ! check positive
		                if ( isequal(Z,s_exp) ) then
		                    !
		                    strmatch = .true.
		                    if (k_exp.eq.1) then
		                        str = char(char_start+k)
		                    else
		                        str = char(char_start+k)//trim(int2char(k_exp))
		                    endif
		                    exit search
		                    !
		                endif
		                ! check negative
		                if ( isequal(Z,-s_exp) ) then
		                    !
		                    strmatch = .true.
		                    if (k_exp.eq.1) then
		                        str = '-'//char(char_start+k)
		                    else
		                        str = '-'//char(char_start+k)//trim(int2char(k_exp))
		                    endif
		                    exit search
		                    !
		                endif
		            enddo
		        enddo search
		        endif
		        !
		        ! if match is not found assign a new character to variable
		        if (.not.strmatch) then
		            !
		            strmatch = .true.
		            k = k + 1
		            s(k) = Z
		            str = char(char_start+k)
		            !
		        endif
		    endif
		    !
		    if (.not.strmatch) then
		        write(str,fmts(5)) real(Z)
		    endif
		    !
	    end subroutine print_chartab_symb
	    subroutine     print_bar(fmts,nclasses)
	    	!
	    	implicit none
	    	!
	        integer		, intent(in) :: nclasses
	        character(*), intent(in) :: fmts(:)
	        integer :: i
	    	!
	        write(*,fmts(2),advance='no') repeat('-',10)
	        do i = 1, nclasses
	            write(*,fmts(3),advance='no') repeat('-',6)
	        enddo
	        write(*,*)
        end subroutine print_bar
    end subroutine print_chartab


end module am_character_table