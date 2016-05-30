module am_symmetry

    use am_constants
    use am_stdout
    use am_unit_cell
    use am_options
    use am_mkl
    use am_symmetry_tables
    use am_matlab
    use dispmodule

    implicit none

    public

    type, public :: am_class_group
        integer :: nsyms                              ! number of symmetries
        integer :: nbases                             ! number of bases functions (dimensions of rep)
        integer , allocatable :: ps_id(:)             ! integer which identifies point symmetries (see decode_pointsymmetry)
        type(am_class_conjugacy_class)      :: cc     ! conjugacy classes
        type(am_class_character_table)      :: ct     ! character table
        type(am_class_multiplication_table) :: mt     ! multiplication table
    contains
        procedure :: sort_symmetries
        procedure :: dump
        procedure :: get_multiplication_table
        procedure :: get_conjugacy_classes
        procedure :: get_character_table
        procedure :: print_character_table
    end type am_class_group

    type, public, extends(am_class_group)       :: am_class_symrep_group
        ! generic symmetry group
        real(dp), allocatable :: sym(:,:,:)           ! symmetry elements [unitless rep]
    end type am_class_symrep_group

    type, public, extends(am_class_group)       :: am_class_seitz_group
        real(dp) :: bas(3,3)                          ! basis in which [frac.] are defined
        real(dp) :: recbas(3,3)                       ! basis in which [recfrac] are defined
        real(dp), allocatable :: seitz_cart(:,:,:)    ! symmetry elements [cart.]
        real(dp), allocatable :: seitz_frac(:,:,:)    ! symmetry elements [frac, direct]; acts on atoms
        real(dp), allocatable :: seitz_recfrac(:,:,:) ! symmetry elements [frac, reciprocal]; acts on kpoints
        ! sym(:,:,:) are seitz operators
        contains
        procedure :: create_seitz_group
        procedure :: write_outfile
        procedure :: write_action_table
    end type am_class_seitz_group

    type, public, extends(am_class_seitz_group) :: am_class_point_group
        ! sym(:,:,:) are seitz operators
        contains
        procedure :: get_point_group
        procedure :: get_rotational_group
        procedure :: get_stabilizer_group
        procedure :: get_reversal_group
    end type am_class_point_group

    type, public, extends(am_class_seitz_group) :: am_class_space_group
        ! sym(:,:,:) are seitz operators
        contains
        procedure :: get_space_group
    end type am_class_space_group

contains

    ! operate on group representation

    subroutine     sort_symmetries(grp,criterion,flags)
        ! flags = 'ascend'/'descend'
        use am_rank_and_sort
        !
        implicit none
        !
        class(am_class_group), intent(inout) :: grp
        real(dp)    , intent(in) :: criterion(:)
        character(*), intent(in) :: flags
        integer , allocatable ::  inds(:)
        integer , allocatable :: rinds(:) ! reverse inds
        integer :: i, j
        !
        allocate(inds(grp%nsyms))
        call rank(criterion,inds)
        !
        if     (index(flags,'descend').ne.0) then
            inds = inds(grp%nsyms:1:-1)
        elseif (index(flags,'ascend' ).ne.0) then
            inds = inds(1:grp%nsyms:+1)
        else
            stop 'sort_symmetries: ascend/descend?'
        endif
        !
        ! get reverse indices
        allocate(rinds(grp%nsyms))
        rinds(inds) = [1:grp%nsyms]
        !
        ! sort symmetries
        select type (grp)
        class is (am_class_symrep_group)
            ! [undefined units]  
            grp%sym = grp%sym(:,:,inds)
        class is (am_class_seitz_group)
            ! symmerties [frac.]
            grp%seitz_frac = grp%seitz_frac(:,:,inds)
            ! symmerties [cart.]
            grp%seitz_cart = grp%seitz_cart(:,:,inds)
            ! symmetries [recfrac]
            grp%seitz_recfrac = grp%seitz_recfrac(:,:,inds)
        class default
            stop 'Class unknown!'
        end select
        !
        ! symmetry id
        if (allocated(grp%ps_id)) then
            grp%ps_id = grp%ps_id(inds)
        endif
        !
        ! conjugacy class
        if (allocated(grp%cc%id)) then
            !
            grp%cc%id = grp%cc%id(inds)
            !
            do i = 1, grp%cc%nclasses
            do j = 1, grp%cc%nelements(i)
            grp%cc%member(i,j) = rinds(grp%cc%member(i,j))
            enddo
            enddo
            !
            do i = 1, grp%cc%nclasses
            grp%cc%representative(i) = rinds(grp%cc%representative(i))
            enddo
            !
        endif
        !
        ! multiplication table
        if (allocated(grp%mt%multab)) then
            do i = 1, grp%nsyms
            do j = 1, grp%nsyms
                grp%mt%multab(i,j) = rinds(grp%mt%multab(i,j))
            enddo
            enddo
            grp%mt%multab = grp%mt%multab(:,inds)
            grp%mt%multab = grp%mt%multab(inds,:)
        endif
        !
    end subroutine sort_symmetries

    subroutine     dump(grp,iopt_filename)
        !
        implicit none
        !
        class(am_class_group), intent(in) :: grp
        character(*), intent(in), optional :: iopt_filename
        character(100) :: fname
        integer :: fid
        !
        ! set default
        fname = 'dump.symmetries'
        if (present(iopt_filename)) fname = iopt_filename
        !
        fid = 1
        open(unit=fid,file=trim(iopt_filename),status="replace",action='write')
            select type (grp)
            class is (am_class_symrep_group)
                ! [unitless]  
                call dump_sym_to_file(fid=fid, sym=grp%sym, header='')
            class is (am_class_seitz_group)
                ! symmerties [frac.]
                call dump_sym_to_file(fid=fid, sym=grp%seitz_frac,    header='[frac.]')
                call dump_sym_to_file(fid=fid, sym=grp%seitz_cart,    header='[cart.]')
                call dump_sym_to_file(fid=fid, sym=grp%seitz_recfrac, header='[rec. frac.]')
            class default
                stop 'Class unknown!'
            end select
        close(fid)
        !
        contains
        subroutine     dump_sym_to_file(fid,sym,header)
            !
            implicit none
            !
            integer, intent(in) :: fid
            real(dp), intent(in) :: sym(:,:,:)
            character(*), intent(in) :: header
            integer :: nbases, nsyms
            integer :: i,j,k
            !
            nsyms  = size(sym,3)
            nbases = size(sym,1)
            !
            do i = 1, nsyms
                write(fid,*) 'symmetry', i, header
                do j = 1, nbases
                    do k = 1, nbases
                        write(fid,"(f)",advance='no') sym(j,k,i)
                    enddo
                    write(fid,*)
                enddo
            enddo
        end subroutine dump_sym_to_file
    end subroutine dump

    subroutine     get_multiplication_table(grp)
        !
        implicit none
        !
        class(am_class_group), intent(inout) :: grp
        character(100) :: flags
        !
        select type(grp)
        class is(am_class_seitz_group)
            flags='seitz'
        class default
            flags=''
        end select
        !
        ! get multiplication table
        select type (grp)
        class is (am_class_symrep_group)
            ! [unitless]  
            grp%mt%multab = get_multab(sym=grp%sym, flags=flags) !//'prog')
        class is (am_class_seitz_group)
            ! symmerties [frac.]
            grp%mt%multab = get_multab(sym=grp%seitz_frac, flags=flags) !//'prog')
        class default
            stop 'unkown class'
        end select
        ! get inverse indices
        grp%mt%inv_id = get_inverse_indices(multab=grp%mt%multab)
        ! get cyclic order
        grp%mt%order = get_cyclic_order(multab=grp%mt%multab)
        ! 
    end subroutine get_multiplication_table

    subroutine     get_conjugacy_classes(grp)
        !
        implicit none
        !
        class(am_class_group), intent(inout) :: grp
        !
        ! check
        if (.not.allocated(grp%mt%multab)) stop 'get_conjugacy_classes: multiplication table is required.'
        ! get conjugacy classes (ps_id is required for sorting the classes: proper before improper)
        grp%cc%id = get_class_id(multab=grp%mt%multab, inv_id=grp%mt%inv_id, ps_id=grp%ps_id)
        ! number of classes
        grp%cc%nclasses = maxval(grp%cc%id)
        ! get number of elements in each class
        grp%cc%nelements = id_nelements(id=grp%cc%id)
        ! get members
        grp%cc%member = id_member(id=grp%cc%id)
        ! get representatives
        grp%cc%representative = id_representative(id=grp%cc%id)
        !
    end subroutine get_conjugacy_classes

    subroutine     get_character_table(grp)
        !
        implicit none
        !
        class(am_class_group), intent(inout) :: grp
        !
        ! check
        if (.not.allocated(grp%mt%multab)) stop 'get_character_table: multiplication table is required.'
        if (.not.allocated(grp%cc%member)) stop 'get_character_table: conjugacy class analyses is required.'
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
        integer :: nreps
        complex(dp) , allocatable :: rep_chi(:,:)
        character(7), allocatable :: rep_label(:)
        !
        select type (grp)
        class is (am_class_point_group) 
            nreps = 15
            allocate(rep_chi(nreps,grp%cc%nclasses))
            allocate(rep_label(nreps))
            ! representation based on basis functions
            rep_label(1) = 'rep'
            rep_chi(1,:) = get_rep_characters(   sym=grp%seitz_frac(1:3,1:3,:), nclasses=grp%cc%nclasses, class_representative=grp%cc%representative)
            ! representation based on s orbitals
            rep_label(2) = 's D^(0)'
            rep_chi(2,:) = get_orb_characters(l=0, R=grp%seitz_cart(1:3,1:3,:), nclasses=grp%cc%nclasses, class_representative=grp%cc%representative)
            ! representation based on p orbitals
            rep_label(3) = 'p D^(1)'
            rep_chi(3,:) = get_orb_characters(l=1, R=grp%seitz_cart(1:3,1:3,:), nclasses=grp%cc%nclasses, class_representative=grp%cc%representative)
            ! representation based on d orbitals
            rep_label(4) = 'd D^(2)'
            rep_chi(4,:) = get_orb_characters(l=2, R=grp%seitz_cart(1:3,1:3,:), nclasses=grp%cc%nclasses, class_representative=grp%cc%representative)
            ! representation based on f orbitals
            rep_label(5) = 'f D^(3)'
            rep_chi(5,:) = get_orb_characters(l=3, R=grp%seitz_cart(1:3,1:3,:), nclasses=grp%cc%nclasses, class_representative=grp%cc%representative)
            ! orbital products involving s orbitals
            rep_label(6) = 's * s'
            rep_chi(6,:) = rep_chi(2,:) * rep_chi(2,:)
            rep_label(7) = 's * p'
            rep_chi(7,:) = rep_chi(2,:) * rep_chi(3,:)
            rep_label(8) = 's * d'
            rep_chi(8,:) = rep_chi(2,:) * rep_chi(4,:)
            rep_label(9) = 's * f'
            rep_chi(9,:) = rep_chi(2,:) * rep_chi(5,:)
            ! orbital products involving p orbitals
            rep_label(10)= 'p * p'
            rep_chi(10,:)= rep_chi(3,:) * rep_chi(3,:)
            rep_label(11)= 'p * d'
            rep_chi(11,:)= rep_chi(3,:) * rep_chi(4,:)
            rep_label(12)= 'p * f'
            rep_chi(12,:)= rep_chi(3,:) * rep_chi(5,:)
            ! orbital products involving d orbitals
            rep_label(13)= 'd * d'
            rep_chi(13,:)= rep_chi(4,:) * rep_chi(4,:)
            rep_label(14)= 'd * f'
            rep_chi(14,:)= rep_chi(4,:) * rep_chi(5,:)
            ! orbital products involving f orbitals
            rep_label(15)= 'f * f'
            rep_chi(15,:)= rep_chi(5,:) * rep_chi(5,:)
        class is (am_class_space_group) 
            ! do nothing.
        class is (am_class_symrep_group)
            nreps = 1
            allocate(rep_chi(nreps,grp%cc%nclasses))
            allocate(rep_label(nreps))
            ! representation based on basis functions
            rep_label(1) = 'rep'
            rep_chi(1,:) = get_rep_characters(sym=grp%sym, nclasses=grp%cc%nclasses, class_representative=grp%cc%representative)
        class default
            stop 'unkown class'
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

    ! high level routines which operate on sg/pg

    subroutine     create_seitz_group(sg,bas,seitz_frac)
        !
        implicit none
        !
        class(am_class_seitz_group), intent(out):: sg
        real(dp), intent(in) :: bas(3,3)
        real(dp), intent(in) :: seitz_frac(:,:,:)
        integer  :: i
        !
        ! basis in which [frac.] are defined
        sg%bas = bas
        ! get reciprocal basis in which [recfrac] are defined
        sg%recbas = inv(bas)
        ! number of "bases functions", doesn't have any real meaning here. Using 3 rather than 4.
        sg%nbases = 3
        ! get number of symmetries
        sg%nsyms = size(seitz_frac,3)
        ! get symmetries [frac.]
        allocate(sg%seitz_frac,source=seitz_frac)
        ! put identity first
        call put_identity_first(sg%seitz_frac)
        ! get symmetries [cart.]
        allocate(sg%seitz_cart(4,4,sg%nsyms))
        do i = 1, sg%nsyms
            sg%seitz_cart(:,:,i) = seitz_frac2cart(seitz_frac=sg%seitz_frac(:,:,i), bas=sg%bas, recbas=sg%recbas)
        enddo
        ! get symmetries [rec. frac.]
        allocate(sg%seitz_recfrac(4,4,sg%nsyms))
        do i = 1, sg%nsyms
            sg%seitz_recfrac(:,:,i) = seitz_cart2recfrac(seitz_cart=sg%seitz_cart(:,:,i), bas=sg%bas, recbas=sg%recbas)
        enddo
        ! correct rounding error
        where(abs(sg%seitz_frac).lt.tiny) sg%seitz_frac = 0
        where(abs(sg%seitz_cart).lt.tiny) sg%seitz_cart = 0
        ! identify point symmetries
        allocate(sg%ps_id(sg%nsyms))
        do i = 1, sg%nsyms
            sg%ps_id(i) = ps_schoenflies(R=sg%seitz_frac(1:3,1:3,i))
        enddo
        ! get multiplication table
        call sg%get_multiplication_table()
        ! get conjugacy classes
        call sg%get_conjugacy_classes()
        ! sort symmetries based on parameters
        call sg%sort_symmetries(criterion=sg%seitz_frac(1,4,:), flags='ascend')
        call sg%sort_symmetries(criterion=sg%seitz_frac(2,4,:), flags='ascend')
        call sg%sort_symmetries(criterion=sg%seitz_frac(3,4,:), flags='ascend')
        call sg%sort_symmetries(criterion=real(sg%ps_id,dp)   , flags='ascend')
        call sg%sort_symmetries(criterion=real(sg%cc%id,dp)   , flags='ascend')
        !get character table
        call sg%get_character_table()
        !
        contains
        function       seitz_frac2cart(seitz_frac,bas,recbas) result(seitz_cart)
            !
            implicit none
            !
            real(dp), intent(in) :: seitz_frac(4,4)
            real(dp), intent(in) :: bas(3,3)
            real(dp), intent(in) :: recbas(3,3)
            real(dp) :: seitz_cart(4,4)
            !
            seitz_cart(1:3,1:3) = matmul(bas,matmul(seitz_frac(1:3,1:3),recbas))
            seitz_cart(1:3,4)   = matmul(bas,seitz_frac(1:3,4))
            seitz_cart(4,4)     = 1.0_dp
            seitz_cart(4,1:3)   = 0.0_dp
            !
        end function   seitz_frac2cart
        function       seitz_cart2recfrac(seitz_cart,bas,recbas) result(seitz_recfrac)
            !
            implicit none
            !
            real(dp), intent(in) :: seitz_cart(4,4)
            real(dp), intent(in) :: bas(3,3)
            real(dp), intent(in) :: recbas(3,3)
            real(dp) :: seitz_recfrac(4,4)
            !
            seitz_recfrac(1:3,1:3) = matmul(bas,matmul(seitz_cart(1:3,1:3),recbas))
            seitz_recfrac(1:3,4)   = matmul(bas,seitz_cart(1:3,4))
            seitz_recfrac(4,4)     = 1.0_dp
            seitz_recfrac(4,1:3)   = 0.0_dp
            !
        end function   seitz_cart2recfrac
    end subroutine create_seitz_group

    subroutine     get_space_group(sg,pc,opts)
        !
        implicit none
        !
        class(am_class_space_group), intent(inout) :: sg
        class(am_class_unit_cell)  , intent(in) :: pc ! space groups are tabulated in the litearture for conventional cells, but primitive or arbitrary cell works just as well.
        type(am_class_options)     , intent(in) :: opts
        real(dp), allocatable  :: seitz_frac(:,:,:)
        character(20) :: str
        integer :: mpl
        integer :: i
        !
        if (opts%verbosity.ge.1) call print_title('Determining space group symmetries')
        !
        ! determine space symmetries from atomic basis [frac.]
        seitz_frac = space_symmetries_from_basis(bas=pc%bas, tau=pc%tau_frac,Z=pc%Z, prec=opts%prec)
        ! put identity first
        call put_identity_first(seitz=seitz_frac)
        ! create space group instance
        call sg%create_seitz_group(seitz_frac=seitz_frac, bas=pc%bas)
        ! print stdout
        if (opts%verbosity.ge.1) then
            write(*,'(a5,a,a)') ' ... ', 'space symmetries = ' , tostring(sg%nsyms)
            write(*,'(a5,a,a)') ' ... ', 'conjugacy classes = ', tostring(sg%cc%nclasses)
            write(*,'(a5,a)')   ' ... ', 'seitz symmetries [frac/cart] = '
            mpl=2 ! matrices per line
            do i = 1,sg%nsyms
                str = tostring(i)//': '//decode_pointsymmetry(sg%ps_id(i))
                if ((mod(i,mpl).eq.1)) then
                    call disp(title=' '      , X=zeros([1,1]),              style='underline',fmt='i2'  ,zeroas=' ',advance='no' , trim='no')
                    call disp(title=trim(str), X=nint(sg%seitz_frac(:,:,i)),style='underline',fmt='i2'  ,zeroas='0',advance='no' , trim='no')
                    call disp(title=trim(str), X=(sg%seitz_cart(:,:,i))    ,style='underline',fmt='f5.2',zeroas='0',advance='no' , trim='no')
                elseif ((mod(i,mpl).eq.0).or.(i.eq.sg%nsyms)) then
                    call disp(title=' '      , X=zeros([1,1]),              style='underline',fmt='i2'  ,zeroas=' ',advance='no' , trim='no')
                    call disp(title=trim(str), X=nint(sg%seitz_frac(:,:,i)),style='underline',fmt='i2'  ,zeroas='0',advance='no' , trim='no')
                    call disp(title=trim(str), X=(sg%seitz_cart(:,:,i))    ,style='underline',fmt='f5.2',zeroas='0',advance='yes', trim='no')
                else
                    call disp(title=' '      , X=zeros([1,1]),              style='underline',fmt='i2'  ,zeroas=' ',advance='no' , trim='no')
                    call disp(title=trim(str), X=nint(sg%seitz_frac(:,:,i)),style='underline',fmt='i2'  ,zeroas='0',advance='no' , trim='no')
                    call disp(title=trim(str), X=(sg%seitz_cart(:,:,i))    ,style='underline',fmt='f5.2',zeroas='0',advance='no' , trim='no')
                endif
            enddo
        endif
        ! write to write_outfile and to file
        call sg%write_outfile(iopt_filename=trim('outfile.spacegroup'))
        call sg%write_action_table(uc=pc,fname='outfile.space_group_action',opts=opts)
        !
        contains
        function        space_symmetries_from_basis(bas,tau,Z,prec) result(seitz)
            !
            !The program first finds the primitive unit cell by looking for
            !additional lattice vectors within the unit cell given in the input.
            !
            ! It chooses one of the atoms and tests each vector from that atom to
            !every other atom of the same type in the unit cell. If that vector
            !will take us from each atom in the unit cell to another atom of the
            !same type (not necessarily in the same unit cell), then we have
            !found a new lattice vector.
            !
            ! This process is repeated until we have
            !found a complete list of lattice vectors. From that list, we find three
            !which can serve as primitive basis vectors ti of the lattice, i.e. every
            !lattice vector can be expressed as an integer linear combination of
            !these basis vectors.
            !
            use am_rank_and_sort
            !
            implicit none
            ! subroutine i/o
            real(dp), intent(in) :: bas(3,3)
            real(dp), intent(in) :: tau(:,:)
            integer , intent(in) :: Z(:)
            real(dp), intent(in) :: prec
            real(dp), allocatable :: seitz(:,:,:)
            real(dp), allocatable :: try(:,:)
            real(dp), allocatable :: T(:,:)
            real(dp), allocatable :: R(:,:,:)
            real(dp), allocatable :: wrkspace(:,:,:)
            real(dp) :: shift(3)
            real(dp) :: T_shifted(3)
            integer :: nTs
            integer :: nRs
            integer :: i, j, k
            !
            R = lattice_symmetries(bas=bas, prec=prec)
            nRs = size(R,3)
            !
            T = translations_from_basis(tau=tau, Z=Z, prec=prec, flags='zero,relax')
            nTs = size(T,2)
            !
            allocate(wrkspace(4,4,nTs*nRs))
            wrkspace = 0.0_dp
            k=0
            do i = 1, nRs
            do j = 1, nTs
                ! Shift tau to put the first atom at the origin. This is important because the
                ! rotational part of the symmetry operation must leave at least one atom in it's
                ! original location. If the atom already has a displacement, then it will be rotated
                ! to another position.
                if (.not.any(all(abs(tau).lt.prec,2))) then
                    shift = tau(:,1)
                    ! tau' = R*tau + T
                    ! tau' - s = R*(tau-s) + T
                    ! tau' = R*tau + (R*S+T-s) = R*tau + T_shifted ! sign is wrong here for some reason.
                    ! thus, T_shifted = T - R*s + s
                    T_shifted = T(:,j) - matmul(R(:,:,i),shift) + shift
                    ! reduce to primitive
                    T_shifted = modulo(T_shifted+opts%prec,1.0_dp)-opts%prec
                else
                    T_shifted = T(:,j)
                endif
                !
                try = eye(4)
                try(1:3,1:3) = R(1:3,1:3,i)
                try(1:3,4) = T_shifted
                !
                if ( is_symmetry_valid(tau=tau, Z=Z, seitz=try, prec=opts%prec)) then
                    k = k + 1
                    wrkspace(:,:,k)=try
                endif
            enddo
            enddo
            !
            allocate(seitz(4,4,k))
            seitz = wrkspace(1:4,1:4,1:k)
            !
            ! print stdout
            if (opts%verbosity.ge.1) then
                call am_print('(im-)proper generators',nRs,' ... ')
                call am_print('translation generators',nTs,' ... ')
                ! call am_print_two_matrices_side_by_side(name='translations',&
                !     Atitle='fractional',A=transpose(T),&
                !     Btitle='cartesian' ,B=transpose(matmul(bas,T)),&
                !     iopt_emph=' ... ',iopt_teaser=.true.)
            endif
        end function    space_symmetries_from_basis
        function        lattice_symmetries(bas,prec) result(R_frac)
            !>
            !> Given a metric tensor, this function returns the (arthimetic) a-holohodry:
            !> the group of point symmetries R (fractional) that are compatible with the
            !> metric tensor).
            !>
            !> Symmetry and Condensed Matter Physics: A Computational Approach. 1 edition.
            !> Cambridge, UK; New York: Cambridge University Press, 2008. pages 278-279,
            !> Eq. 10.60 and 10.61
            !>
            !> The Mathematical Theory of Symmetry in Solids: Representation Theory
            !> for Point Groups and Space Groups. 1 edition. Oxford?: New York:
            !> Oxford University Press, 2010. page 134, table 3.9
            !>
            !> Structure of Materials: An Introduction to Crystallography, Diffraction
            !> and Symmetry. 2 edition. New York: Cambridge University Press, 2012.
            !>
            !> page 169:
            !> "If [point symmetries] are described with respect to the 'primed' basis
            !>  vectors (fractional coordinates), then once again all [point symmetries]
            !>  can be represented by matrices which contain only -1, 0, +1"
            !>
            !> page 166:
            !> "The inverse of a unitary matrix (rotation) is equal to its transpose
            !>  (this is only true in cartesian coordinates.)"
            !>
            implicit none
            ! subroutine i/o
            real(dp), intent(in) :: bas(3,3)
            real(dp), intent(in) :: prec
            real(dp), allocatable :: R_frac(:,:,:) !> point symmetries (fractional)
            !
            real(dp) :: id(3,3)
            real(dp) :: recbas(3,3)
            real(dp) :: metric(3,3)
            real(dp) :: buffer(3,3,48) ! buffer for point symmetries
            integer  :: k ! point symmetry counter
            real(dp) :: o(3,3)  ! point symmetry in fractional
            integer  :: i11, i12, i13, i21, i22, i23, i31, i32, i33 ! used to generate unitary rotational matrices
            !
            id = eye(3)
            !
            recbas = inv(bas)
            !
            metric = matmul(transpose(bas),bas)
            !
            buffer = 0.0_dp
            !
            k = 0
            !
            do i11 = -1,1
            do i12 = -1,1
            do i13 = -1,1
                do i21 = -1,1
                do i22 = -1,1
                do i23 = -1,1
                    do i31 = -1,1
                    do i32 = -1,1
                    do i33 = -1,1
                        !
                        o(1,1:3)=real([i11,i12,i13],dp)
                        o(2,1:3)=real([i21,i22,i23],dp)
                        o(3,1:3)=real([i31,i32,i33],dp)
                        !
                        ! Check that metric is left unchanged by symmetry opreation in fractional coordinates.
                        ! i.e. that metric tensor commutes with point symmetry.
                        if ( all( abs(matmul(transpose(o),matmul(metric,o))-metric) .lt. prec ) ) then
                            k = k + 1
                            ! store point symmetry in fractional coordinates
                            buffer(:,:,k) = o
                        endif
                        !
                    enddo
                    enddo
                    enddo
                enddo
                enddo
                enddo
            enddo
            enddo
            enddo
            !
            allocate(R_frac,source=buffer(:,:,1:k))
            !
        end function    lattice_symmetries
    end subroutine get_space_group

    subroutine     get_point_group(pg,pc,sg,opts)
        !
        implicit none
        !
        class(am_class_point_group), intent(inout) :: pg
        type(am_class_space_group) , intent(in) :: sg
        class(am_class_unit_cell)  , intent(in) :: pc  ! accepts unit cell, primitive cell, conventional cell.
        type(am_class_options)     , intent(in) :: opts
        real(dp), allocatable  :: wrk_frac(:,:,:)
        character(20) :: str
        integer :: mpl
        integer :: i
        !
        if (opts%verbosity.ge.1) call print_title('Point group symmetries')
        !
        ! create point group instance
        allocate(wrk_frac,source=sg%seitz_frac)
        wrk_frac(1:3,4,:) = 0
        call pg%create_seitz_group(seitz_frac=unique(wrk_frac), bas=pc%bas)
        !
        if (opts%verbosity.ge.1) then
            write(*,'(a5,a,a)') ' ... ', 'point group = '      , decode_pointgroup(point_group_schoenflies(pg%ps_id))
            write(*,'(a5,a,a)') ' ... ', 'point symmetries = ' , tostring(pg%nsyms)
            write(*,'(a5,a,a)') ' ... ', 'conjugacy classes = ', tostring(pg%cc%nclasses)
            write(*,'(a5,a)')   ' ... ', 'point symmetries [frac/cart] = '
            mpl=2 ! matrices per line
            do i = 1,pg%nsyms
                str = tostring(i)//': '//decode_pointsymmetry(pg%ps_id(i))
                if ((mod(i,mpl).eq.1)) then
                    call disp(title=' '      , X=zeros([1,1]),              style='underline',fmt='i2'  ,zeroas=' ',advance='no' , trim='no')
                    call disp(title=trim(str), X=nint(pg%seitz_frac(:,:,i)),style='underline',fmt='i2'  ,zeroas='0',advance='no' , trim='no')
                    call disp(title=trim(str), X=(pg%seitz_cart(:,:,i))    ,style='underline',fmt='f5.2',zeroas='0',advance='no' , trim='no')
                elseif ((mod(i,mpl).eq.0).or.(i.eq.pg%nsyms)) then
                    call disp(title=' '      , X=zeros([1,1]),              style='underline',fmt='i2'  ,zeroas=' ',advance='no' , trim='no')
                    call disp(title=trim(str), X=nint(pg%seitz_frac(:,:,i)),style='underline',fmt='i2'  ,zeroas='0',advance='no' , trim='no')
                    call disp(title=trim(str), X=(pg%seitz_cart(:,:,i))    ,style='underline',fmt='f5.2',zeroas='0',advance='yes', trim='no')
                else
                    call disp(title=' '      , X=zeros([1,1]),              style='underline',fmt='i2'  ,zeroas=' ',advance='no' , trim='no')
                    call disp(title=trim(str), X=nint(pg%seitz_frac(:,:,i)),style='underline',fmt='i2'  ,zeroas='0',advance='no' , trim='no')
                    call disp(title=trim(str), X=(pg%seitz_cart(:,:,i))    ,style='underline',fmt='f5.2',zeroas='0',advance='no' , trim='no')
                endif
            enddo
            !
            call print_title('Point group character table')
            call pg%print_character_table()
        endif
        !
        call pg%write_outfile(iopt_filename='outfile.pointgroup')
        !
        call pg%write_action_table(uc=pc,fname='outfile.point_group_action',opts=opts)
        !
    end subroutine get_point_group

    subroutine     get_rotational_group(rotg,uc,pg,opts,flags)
        !
        ! Get point symmetries which are compatible with the atomic basis (essentially keeps
        ! rotational part of space symmetries which are able to map all atoms onto each other)
        !
        ! This is also useful for determining the point symmetries compatible with kpoint mesh
        !
        ! Name "point group" - NOTE: not really a point group because translational components of
        ! space symmetries were ignored in the genereation process this would be the name of the
        ! group, if the group were symmorphic (in which case it will match the real point group
        ! name)...
        !
        implicit none
        ! subroutine i/o
        class(am_class_point_group),intent(out):: rotg
        class(am_class_unit_cell),  intent(in) :: uc
        type(am_class_point_group), intent(in) :: pg
        type(am_class_options),     intent(in) :: opts
        character(*)              , intent(in) :: flags
        integer, allocatable   :: indicies(:)
        logical, allocatable   :: mask(:)
        integer :: i
        !
        allocate(mask(pg%nsyms))
        mask = .false. 
        !
        if      (index(flags,'cart').ne.0) then
            do i = 1, pg%nsyms
                if ( is_symmetry_valid(tau=uc%tau_cart, Z=uc%Z, seitz=pg%seitz_cart, prec=opts%prec, flags='cart')) then
                    mask(i) = .true.
                endif
            enddo
        else
            stop 'Note: rotational group only makes sense for [cart], for which symmetry operations are strictly unitary.'
        endif
        ! create indices
        allocate(indicies(pg%nsyms))
        indicies=[1:pg%nsyms]
        !
        ! create group
        call rotg%create_seitz_group(seitz_frac=pg%seitz_frac(:,:,pack(indicies,mask)), bas=pg%bas)
        !
    end subroutine get_rotational_group

    subroutine     get_stabilizer_group(stab,pg,v,opts,flags)
        !
        implicit none
        !
        class(am_class_point_group), intent(out):: stab ! stabilizer group associated with vector v
        type(am_class_point_group) , intent(in) :: pg
        real(dp)                   , intent(in) :: v(3) ! vector which is stabilized (should be same units as R, which is fractional)
        type(am_class_options)     , intent(in) :: opts
        character(*)               , intent(in) :: flags
        integer, allocatable :: indicies(:)
        logical, allocatable :: mask(:)
        integer :: i
        !
        allocate(mask(pg%nsyms))
        mask = .false.
        !
        if      (index(flags,'frac').ne.0) then
            do i = 1, pg%nsyms
                if (is_symmetry_valid(tau=reshape(v,[3,1]), Z=[1], seitz=pg%seitz_frac(:,:,i), prec=opts%prec, flags=flags//'exact')) then
                    mask(i) = .true.
                endif
            enddo
        elseif  (index(flags,'cart').ne.0) then
            do i = 1, pg%nsyms
                if (is_symmetry_valid(tau=reshape(v,[3,1]), Z=[1], seitz=pg%seitz_cart(:,:,i), prec=opts%prec, flags=flags//'exact')) then
                    mask(i) = .true.
                endif
            enddo
        else
            stop 'flags must contain frac/cart'
        endif
        !
        allocate(indicies(pg%nsyms))
        indicies=[1:pg%nsyms]
        !
        ! create group
        call stab%create_seitz_group(seitz_frac=pg%seitz_frac(:,:,pack(indicies,mask)), bas=pg%bas)
        !
    end subroutine get_stabilizer_group
    
    subroutine     get_reversal_group(revg,pg,v,opts,flags)
        !
        ! there is a bug in the reversal group. Si-Si first naerest neighbor shells on different primitive cell atoms does not match.
        !
        !  shell Zi-Zj   i-j   m-n    m   stab.    rot.    rev. |v(cart)|           v(cart)             |v(frac)|           v(frac)
        !  ----- ----- ----- ----- ---- ------- ------- ------- --------- ----------------------------- --------- -----------------------------
        !      1 Si-Si   1-1   1-1    1     o_h     o_h     t_d     0.000    -0.000    -0.000    -0.000     0.000    -0.000    -0.000    -0.000
        !      2 Si-Si   1-1   1-2    4    c_3v     o_h    d_3d     2.351     1.357     1.357     1.357     0.217     0.125     0.125     0.125
        !  primitive atom 2 at 0.12,0.12,0.12 (frac) has 2 nearest-neighbor shells
        !  shell Zi-Zj   i-j   m-n    m   stab.    rot.    rev. |v(cart)|           v(cart)             |v(frac)|           v(frac)
        !  ----- ----- ----- ----- ---- ------- ------- ------- --------- ----------------------------- --------- -----------------------------
        !      3 Si-Si   1-1   2-2    1     o_h     o_h     t_d     0.000     0.000     0.000     0.000     0.000     0.000     0.000     0.000
        !      4 Si-Si   1-1   2-1    4    c_3v     o_h    c_3v     2.351    -1.357    -1.357    -1.357     0.217    -0.125    -0.125    -0.125
        !
        implicit none
        !
        class(am_class_point_group), intent(out):: revg ! reversal group associated with vector v, i.e. the space symmetries with bring v to [0,0,0], essentially flipping the atoms at the corners of the vector
        type(am_class_point_group) , intent(in) :: pg ! space group
        real(dp)                   , intent(in) :: v(3)
        type(am_class_options)     , intent(in) :: opts
        character(*)               , intent(in) :: flags
        integer, allocatable :: indicies(:)
        logical, allocatable :: mask(:)
        real(dp),allocatable :: seitz_frac(:,:,:)
        integer :: i
        !
        allocate(mask(pg%nsyms))
        mask = .false.
        !
        if      (index(flags,'frac').ne.0) then
            do i = 1, pg%nsyms
                if (is_symmetry_valid(tau=reshape(v,[3,1]), Z=[1], seitz=pg%seitz_frac(:,:,i), prec=opts%prec, flags=flags//'reversal')) then
                    mask(i) = .true.
                endif
            enddo
        elseif  (index(flags,'cart').ne.0) then
            do i = 1, pg%nsyms
                if (is_symmetry_valid(tau=reshape(v,[3,1]), Z=[1], seitz=pg%seitz_cart(:,:,i), prec=opts%prec, flags=flags//'reversal')) then
                    mask(i) = .true.
                endif
            enddo
        else
            stop 'flags must contain frac/cart'
        endif
        !
        allocate(indicies(pg%nsyms))
        indicies=[1:pg%nsyms]
        !
        allocate(seitz_frac, source=pg%seitz_frac(:,:,pack(indicies,mask)))
        seitz_frac(1:3,4,:) = 0 ! set all translations to zero
        !
        ! create group
        call revg%create_seitz_group(seitz_frac=unique(seitz_frac), bas=pg%bas)
        !
    end subroutine get_reversal_group

    subroutine     write_outfile(sg,iopt_filename)
        !
        implicit none
        !
        class(am_class_seitz_group), intent(in) ::  sg
        character(*), intent(in), optional :: iopt_filename
        integer :: width
        real(dp) :: R(3,3), T(3)
        integer :: fid
        integer :: i,j,k
        character(15) :: fmt5, fmt4, fmt6, fmt1, fmt2, fmt3
        !
        if ( present(iopt_filename) ) then
            ! write to file
            fid = 1
            open(unit=fid,file=trim(iopt_filename),status="replace",action='write')
        else
            ! write_outfile
            fid = 6
        endif
            !
            ! write standard out to file
            fmt5=   "(a4)"
            fmt4=   "(a4,3a7)"
            fmt6=       "(a7)"
            fmt1="(5x,i3,9a7)"
            fmt2="(5x,a3,9a7)"
            fmt3="(5x,a)"
            width = 9*7+3
            width = width + 3*7+4
            !
            !
            ! fractional
            write(fid,fmt3) centertitle('seitz [frac.]',width)
            write(fid,fmt3) repeat('-',width)
            write(fid,fmt2,advance='no') '#', 'R11', 'R21', 'R31', 'R12', 'R22', 'R32', 'R13', 'R23', 'R33'
            write(fid,fmt4,advance='no') '  | ', 'T1', 'T2', 'T3'
            write(fid,*)
            do i = 1, sg%nsyms
                !
                R = sg%seitz_frac(1:3,1:3,i)
                write(fid,fmt1,advance='no') i, (( trim(print_pretty(R(j,k))), j = 1,3) , k = 1,3)
                !
                T = sg%seitz_frac(1:3,4,i)
                write(fid,fmt5,advance='no') '  | '
                do j = 1,3
                    write(fid,fmt6,advance='no') trim(print_pretty(T(j)))
                enddo
                write(fid,*)
                !
            enddo
            !
            ! cartesian
            write(fid,fmt3) centertitle('seitz [cart.]',width)
            write(fid,fmt3) repeat('-',width)
            write(fid,fmt2,advance='no') '#', 'R11', 'R21', 'R31', 'R12', 'R22', 'R32', 'R13', 'R23', 'R33'
            write(fid,fmt4,advance='no') '  | ', 'T1', 'T2', 'T3'
            write(fid,*)
            !
            do i = 1, sg%nsyms
                !
                R = sg%seitz_cart(1:3,1:3,i)
                write(fid,fmt1,advance='no') i, (( trim(print_pretty(R(j,k))), j = 1,3) , k = 1,3)
                !
                T = sg%seitz_cart(1:3,4,i)
                write(fid,fmt5,advance='no') '  | '
                do j = 1,3
                    write(fid,fmt6,advance='no') trim(print_pretty(T(j)))
                enddo
                write(fid,*)
            enddo
            !
            ! recfrac
            write(fid,fmt3) centertitle('seitz [rec. frac.]',width)
            write(fid,fmt3) repeat('-',width)
            write(fid,fmt2,advance='no') '#', 'R11', 'R21', 'R31', 'R12', 'R22', 'R32', 'R13', 'R23', 'R33'
            write(fid,fmt4,advance='no') '  | ', 'T1', 'T2', 'T3'
            write(fid,*)
            !
            do i = 1, sg%nsyms
                !
                R = sg%seitz_cart(1:3,1:3,i)
                write(fid,fmt1,advance='no') i, (( trim(print_pretty(R(j,k))), j = 1,3) , k = 1,3)
                !
                T = sg%seitz_cart(1:3,4,i)
                write(fid,fmt5,advance='no') '  | '
                do j = 1,3
                    write(fid,fmt6,advance='no') trim(print_pretty(T(j)))
                enddo
                write(fid,*)
            enddo
        ! close file
        if ( present(iopt_filename) ) close(fid)
        !
    end subroutine write_outfile

    subroutine     write_action_table(sg,uc,fname,opts)
        !
        implicit none
        !
        class(am_class_seitz_group), intent(in) :: sg
        class(am_class_unit_cell),   intent(in) :: uc
        character(*),                intent(in) :: fname
        type(am_class_options),      intent(in) :: opts
        real(dp), allocatable :: R(:,:)   ! used for basis functions output
        character(:), allocatable :: xyz(:) ! basis functions output
        integer, allocatable :: PM(:,:)
        integer, allocatable :: P(:,:,:)
        real(dp) :: aa(4)  ! axis and angle 
        integer  :: fid
        integer  :: i,j,k ! loop variables
        character(10) :: buffer
        character(1)  :: s ! sign 
        character(10) :: fmt1
        character(10) :: fmt2
        character(10) :: fmt3
        character(10) :: fmt4
        character(10) :: fmt5
        character(10) :: fmt6
        !
        fmt1 = '(a5)'
        fmt2 = '(i5)'
        fmt3 = '(a10)'
        fmt4 = '(f10.2)'
        fmt5 = '(i10)'
        fmt6 = '(l10)'
        !
        fid = 1
        open(unit=fid,file=trim(fname),status="replace",action='write')
            !
            ! SYMMETRY NUMBER
            write(fid,'(5x)',advance='no')
            write(fid,fmt1,advance='no') '#'
            do i = 1, sg%nsyms
                write(fid,fmt5,advance='no') i
            enddo
            write(fid,*)
            ! CLASS IDENTIFIED
            write(fid,'(5x)',advance='no')
            write(fid,fmt1,advance='no') 'class'
            do i = 1, sg%nsyms
                write(fid,fmt5,advance='no') sg%cc%id(i)
            enddo
            write(fid,*)
            ! HEADER / SEPERATOR
            write(fid,'(5x)',advance='no')
            write(fid,'(5x)',advance='no')
            do i = 1, sg%nsyms
                write(fid,fmt3,advance='no') ' '//repeat('-',9)
            enddo
            write(fid,*)
            ! POINT-SYMMETRY IDENTIFIED
            write(fid,'(5x)',advance='no')
            write(fid,fmt1,advance='no') 'schf.'
            do i = 1, sg%nsyms
                write(fid,fmt3,advance='no') trim(decode_pointsymmetry(sg%ps_id(i)))
            enddo
            write(fid,*)
            ! POINT SYMMETRY COMPONENTS
            do i = 1,3
            do j = 1,3
                write(fid,'(5x)',advance='no')
                write(fid,fmt1,advance='no') adjustr('R'//trim(int2char(i))//trim(int2char(j)))
                do k = 1, sg%nsyms
                    write(fid,fmt4,advance='no') sg%seitz_cart(i,j,k)
                enddo
                write(fid,*)
            enddo
            enddo
            ! HEADER / SEPERATOR
            write(fid,'(5x)',advance='no')
            write(fid,fmt1,advance='no') 'ax/th'
            do i = 1, sg%nsyms
                write(fid,fmt3,advance='no') ' '//repeat('-',9)
            enddo
            write(fid,*)
            ! POINT-SYMMETRY ROTATION AXIS
            allocate(character(1)::xyz(3))
            i=0
            i=i+1; xyz(i)='x'
            i=i+1; xyz(i)='y'
            i=i+1; xyz(i)='z'
            do j = 1,3
            write(fid,'(5x)',advance='no')
            write(fid,fmt1,advance='no') 'ax_'//trim(xyz(j))
            do i = 1, sg%nsyms
                aa = rot2axis_angle(R=sg%seitz_cart(1:3,1:3,i))
                write(fid,fmt4,advance='no') aa(j)
            enddo
            write(fid,*)
            enddo
            deallocate(xyz)
            ! POINT-SYMMETRY ROTATION ANGLE
            write(fid,'(5x)',advance='no')
            write(fid,fmt1,advance='no') 'R_th'
            do i = 1, sg%nsyms
                aa = rot2axis_angle(R=sg%seitz_cart(1:3,1:3,i))
                write(fid,fmt4,advance='no') aa(4)*180.0_dp/pi
            enddo
            write(fid,*)
            ! POINT-SYMMETRY ROTOINVERSION VS INVERSION
            write(fid,'(5x)',advance='no')
            write(fid,fmt1,advance='no') 'inv'
            do i = 1, sg%nsyms
                write(fid,fmt6,advance='no') (abs(det(sg%seitz_cart(1:3,1:3,i))+1).lt.tiny)
            enddo
            write(fid,*)
            ! HEADER / SEPERATOR
            write(fid,'(5x)',advance='no')
            write(fid,fmt1,advance='no') 'trans'
            do i = 1, sg%nsyms
                write(fid,fmt3,advance='no') ' '//repeat('-',9)
            enddo
            write(fid,*)
            ! TRANSLATIONAL COMPONENTS (FRAC)
            do j = 1,3
                write(fid,'(5x)',advance='no')
                write(fid,fmt1,advance='no') adjustr('T'//trim(int2char(j)))
                do i = 1, sg%nsyms
                    write(fid,fmt4,advance='no') sg%seitz_cart(j,4,i)
                enddo
                write(fid,*)
            enddo
            ! ACTION chartab (AXES)
            write(fid,'(5x)',advance='no')
            write(fid,fmt1,advance='no') 'axes'
            do i = 1, sg%nsyms
                write(fid,fmt3,advance='no') ' '//repeat('-',9)
            enddo
            write(fid,*)
            allocate(character(1)::xyz(3))
            i=0
            i=i+1; xyz(i)='x'
            i=i+1; xyz(i)='y'
            i=i+1; xyz(i)='z'
            do j = 1,3
                write(fid,'(5x)',advance='no')
                write(fid,fmt1,advance='no') trim(xyz(j))
                do i = 1, sg%nsyms
                    ! inverse here because f(Rr) = R^-1 * f(r)
                    ! for fractional R; its elements will always be an integer...
                    R=inv(sg%seitz_cart(1:3,1:3,i))
                    ! if statement are for the cases in which output needs to be trimmed
                    if (count(abs(R(j,:)).gt.tiny).gt.3) then
                        buffer=' ... ' 
                    else
                        buffer=''
                        do k = 1,3
                        if (abs(R(j,k)).gt.tiny) then
                            ! s is only the sign
                            s = trim(int2char(nint(R(j,k)),'SP'))
                            buffer = trim(buffer)//trim(s)//trim(xyz(k))
                            !
                        endif
                        enddo
                    endif
                    write(fid,fmt3,advance='no') trim(buffer)
                enddo
                write(fid,*)
            enddo
            deallocate(xyz)
            ! ACTION chartab (BASIS FUNCTIONS FOR P ORBITALS)
            write(fid,'(5x)',advance='no')
            write(fid,fmt1,advance='no') 'p orbs.'
            do i = 1, sg%nsyms
                write(fid,fmt3,advance='no') ' '//repeat('-',9)
            enddo
            write(fid,*)
            allocate(character(1)::xyz(3))
            i=0
            i=i+1; xyz(i)='x'
            i=i+1; xyz(i)='y'
            i=i+1; xyz(i)='z'
            do j = 1, size(xyz)
                write(fid,'(5x)',advance='no')
                write(fid,fmt1,advance='no') trim(xyz(j))
                do i = 1, sg%nsyms
                    !
                    R = rot2O3(l=1,R=sg%seitz_cart(1:3,1:3,i))
                    ! if statement are for the cases in which output needs to be trimmed
                    if (count(abs(R(j,:)).gt.tiny).gt.3) then
                        buffer=' ... ' 
                    else
                        buffer=''
                        do k = 1, size(xyz)
                            if (abs(R(j,k)).gt.tiny) then
                                ! s is only the sign
                                s = trim(int2char(nint(R(j,k)),'SP'))
                                buffer = trim(buffer)//trim(s)//trim(xyz(k))
                                !
                            endif
                        enddo
                    endif
                    write(fid,fmt3,advance='no') trim(buffer)
                enddo
                write(fid,*)
            enddo
            deallocate(xyz)
            ! ACTION chartab (BASIS FUNCTIONS FOR D ORBITALS)
            write(fid,'(5x)',advance='no')
            write(fid,fmt1,advance='no') 'd orbs.'
            do i = 1, sg%nsyms
                write(fid,fmt3,advance='no') ' '//repeat('-',9)
            enddo
            write(fid,*)
            allocate(character(2)::xyz(5))
            i=0
            i=i+1; xyz(i)='xx'
            i=i+1; xyz(i)='xy'
            i=i+1; xyz(i)='yy'
            i=i+1; xyz(i)='yz'
            i=i+1; xyz(i)='zz'
            do j = 1,5
                write(fid,'(5x)',advance='no')
                write(fid,fmt1,advance='no') trim(xyz(j))
                do i = 1, sg%nsyms
                    !
                    R = rot2O3(l=2,R=sg%seitz_cart(1:3,1:3,i))
                    !
                    if (count(abs(R(j,:)).gt.tiny).gt.3) then
                        buffer=' ... ' 
                    else
                        buffer=''
                        do k = 1,5
                            if (abs(R(j,k)).gt.tiny) then
                                ! s is only the sign
                                s = trim(int2char(nint(R(j,k)),'SP'))
                                buffer = trim(buffer)//trim(s)//trim(xyz(k))
                                !
                            endif
                        enddo
                    endif
                    write(fid,fmt3,advance='no') trim(buffer)
                enddo
                write(fid,*)
            enddo
            deallocate(xyz)
            ! ACTION chartab (BASIS FUNCTIONS FOR D ORBITALS)
            write(fid,'(5x)',advance='no')
            write(fid,fmt1,advance='no') 'f orbs.'
            do i = 1, sg%nsyms
                write(fid,fmt3,advance='no') ' '//repeat('-',9)
            enddo
            write(fid,*)
            allocate(character(1)::xyz(7))
            i=0
            i=i+1; xyz(i)='1'
            i=i+1; xyz(i)='2'
            i=i+1; xyz(i)='3'
            i=i+1; xyz(i)='4'
            i=i+1; xyz(i)='5'
            i=i+1; xyz(i)='6'
            i=i+1; xyz(i)='7'
            do j = 1,7
                write(fid,'(5x)',advance='no')
                write(fid,fmt1,advance='no') trim(xyz(j))
                do i = 1, sg%nsyms
                    !
                    R = rot2O3(l=3,R=sg%seitz_cart(1:3,1:3,i))
                    !
                    if (count(abs(R(j,:)).gt.tiny).gt.3) then
                        buffer=' ... ' 
                    else
                        buffer=''
                        do k = 1,7
                            if (abs(R(j,k)).gt.tiny) then
                                ! s is only the sign
                                s = trim(int2char(nint(R(j,k)),'SP'))
                                buffer = trim(buffer)//trim(s)//trim(xyz(k))
                                !
                            endif
                        enddo
                    endif
                    write(fid,fmt3,advance='no') trim(buffer)
                enddo
                write(fid,*)
            enddo
            deallocate(xyz)
            ! ACTION chartab (ATOM PERMUTATION)
            select type (sg)
            type is (am_class_space_group)
                !
                P  = permutation_rep(seitz=sg%seitz_frac,tau=uc%tau_frac,flags='',prec=opts%prec)
                PM = permutation_map(P)
                !
                write(fid,'(5x)',advance='no')
                write(fid,fmt1,advance='no') 'atoms'
                do i = 1, sg%nsyms
                    write(fid,fmt3,advance='no') ' '//repeat('-',9)
                enddo
                write(fid,*)
                !
                do j = 1, uc%natoms
                    write(fid,'(5x)',advance='no')
                    write(fid,fmt2,advance='no') j
                    do i = 1, sg%nsyms
                        write(fid,fmt5,advance='no') PM(j,i)
                    enddo
                    write(fid,*)
                enddo
            end select
            !
            write(fid,*)
            write(fid,'(a)') 'p,d,f orbital expansion coefficients are omited for clarity. Also omitted are orbitals which decompose into more than three terms.'
            write(fid,'(a)') 'Cartesian coordinates are used throughout this output.'
            !
        ! close file
        close(fid)
    end subroutine write_action_table

!     subroutine     print_tensor(seitz)
!         !
!         implicit none
!         !
!         real(dp), target, intent(in) :: seitz(:,:,:)
!         integer :: i,j,k
!         integer :: m,n,o
!         integer :: t
!         integer :: matrix_per_line
!         real(dp), pointer :: p(:,:,:)
!         !
!         m=size(seitz,1)
!         n=size(seitz,2)
!         o=size(seitz,3)
!         !
!         t = 0
!         do i = 1, m
!             do k = 1, o ! get three tensors per line
!                 do j = 1, n
!                     write(*,'(f8.2)') seitz(i,j,k)
!                 enddo
!                     write(*,'(5x)')
!                 enddo

!             enddo

!         enddo
!         !
!     end subroutine print_tensor

    ! decoding/naming functions 

    function       ps_schoenflies(R) result(ps_id)
        !
        ! Point symmetries in fractional coordinates so that they are nice integers which can be easily classified.
        !
        !    element:         e    i  c_2  c_3  c_4  c_6  s_2  s_6  s_4  s_3
        !    trace:          +3   -3   -1    0   +1   +2   +1    0   -1   -2
        !    determinant:    +1   -1   +1   +1   +1   +1   -1   -1   -1   -1
        !
        implicit none
        !
        real(dp), intent(in) :: R(3,3)
        real(dp) :: tr
        real(dp) :: d
        integer :: ps_id
        !
        ! get trace and determinant (fractional)
        tr  = trace(R)
        d = det(R)
        !
        ! The Mathematical Theory of Symmetry in Solids: Representation Theory for
        ! Point Groups and Space Groups. 1 edition. Oxford?: New York: Oxford University
        ! Press, 2010. page 138, chartab 3.8.
        !
        ps_id = 0
        if     ( (abs(tr - 3).lt.tiny) .and. (abs(d - 1).lt.tiny) ) then; ps_id = 1  ! 'e'
        elseif ( (abs(tr + 1).lt.tiny) .and. (abs(d - 1).lt.tiny) ) then; ps_id = 2  ! 'c_2'
        elseif ( (abs(tr - 0).lt.tiny) .and. (abs(d - 1).lt.tiny) ) then; ps_id = 3  ! 'c_3'
        elseif ( (abs(tr - 1).lt.tiny) .and. (abs(d - 1).lt.tiny) ) then; ps_id = 4  ! 'c_4'
        elseif ( (abs(tr - 2).lt.tiny) .and. (abs(d - 1).lt.tiny) ) then; ps_id = 5  ! 'c_6'
        elseif ( (abs(tr + 3).lt.tiny) .and. (abs(d + 1).lt.tiny) ) then; ps_id = 6  ! 'i'
        elseif ( (abs(tr - 1).lt.tiny) .and. (abs(d + 1).lt.tiny) ) then; ps_id = 7  ! 's_2'
        elseif ( (abs(tr - 0).lt.tiny) .and. (abs(d + 1).lt.tiny) ) then; ps_id = 8  ! 's_6'
        elseif ( (abs(tr + 1).lt.tiny) .and. (abs(d + 1).lt.tiny) ) then; ps_id = 9  ! 's_4'
        elseif ( (abs(tr + 2).lt.tiny) .and. (abs(d + 1).lt.tiny) ) then; ps_id = 10 ! 's_3'
        endif
        !
        if (ps_id.eq.0) stop 'Unable to identify point symmetry.'
        !
    end function   ps_schoenflies

    function       decode_pointgroup(pg_code) result(point_group_name)
        ! Hermann-Mauguin     Shubnikov[1]        Schoenflies     Orbifold    Coxeter 
        ! 1                   1\                  C1              11          [ ]+    
        ! -1                  \tilde{2}           S2 (=Ci)        ×           [2+,2+] 
        ! 2                   2\                  C2              22          [2]+    
        ! m                   m\                  C1h             *           [ ]     
        ! 2/m                 2:m\                C2h             2*          [2,2+]  
        ! 222                 2:2\                D2              222         [2,2]+  
        ! mm2                 2 \cdot m\          C2v             *22         [2]     
        ! mmm                 m \cdot 2:m\        D2h             *222        [2,2]   
        ! 4                   4\                  C4              44          [4]+    
        ! -4                  \tilde{4}           S4              2×          [2+,4+] 
        ! 4/m                 4:m\                C4h             4*          [2,4+]  
        ! 422                 4:2\                D4              422         [4,2]+  
        ! 4mm                 4 \cdot m\          C4v             *44         [4]     
        ! 42m                 \tilde{4}\cdot m    D2d             2*2         [2+,4]  
        ! 4/mmm               m \cdot 4:m\        D4h             *422        [4,2]   
        ! 3                   3\                  C3              33          [3]+    
        ! -3                  \tilde{6}           S6 (=C3i)       3×          [2+,6+] 
        ! 32                  3:2\                D3              322         [3,2]+  
        ! 3m                  3 \cdot m\          C3v             *33         [3]     
        ! -3m                 \tilde{6}\cdot m    D3d             2*3         [2+,6]  
        ! 6                   6\                  C6              66          [6]+    
        ! -6                  3:m\                C3h             3*          [2,3+]  
        ! 6/m                 6:m\                C6h             6*          [2,6+]  
        ! 622                 6:2\                D6              622         [6,2]+  
        ! 6mm                 6 \cdot m\          C6v             *66         [6]     
        ! 6m2                 m \cdot 3:m\        D3h             *322        [3,2]   
        ! 6/mmm               m \cdot 6:m\        D6h             *622        [6,2]   
        ! 23                  3/2\                T               332         [3,3]+  
        ! m3                  \tilde{6}/2         Th              3*2         [3+,4]  
        ! 432                 3/4\                O               432         [4,3]+  
        ! 43m                 3/\tilde{4}         Td              *332        [3,3]   
        ! m3m                 \tilde{6}/4         Oh              *432        [4,3]   
        !
        implicit none
        !
        integer, intent(in) :: pg_code
        character(string_length_schoenflies) :: point_group_name
        !
        select case (pg_code)
                case(1);  point_group_name = 'c_1'
                case(2);  point_group_name = 's_2'
                case(3);  point_group_name = 'c_2'
                case(4);  point_group_name = 'c_1h'
                case(5);  point_group_name = 'c_2h'
                case(6);  point_group_name = 'd_2'
                case(7);  point_group_name = 'c_2v'
                case(8);  point_group_name = 'd_2h'
                case(9);  point_group_name = 'c_3'
                case(10); point_group_name = 's_6'
                case(11); point_group_name = 'd_3'
                case(12); point_group_name = 'c_3v'
                case(13); point_group_name = 'd_3d'
                case(14); point_group_name = 'c_4'
                case(15); point_group_name = 's_4'
                case(16); point_group_name = 'c_4h'
                case(17); point_group_name = 'd_4'
                case(18); point_group_name = 'c_4v'
                case(19); point_group_name = 'd_2d'
                case(20); point_group_name = 'd_4h'
                case(21); point_group_name = 'c_6'
                case(22); point_group_name = 'c_3h'
                case(23); point_group_name = 'c_6h'
                case(24); point_group_name = 'd_6'
                case(25); point_group_name = 'c_6v'
                case(26); point_group_name = 'd_3h'
                case(27); point_group_name = 'd_6h'
                case(28); point_group_name = 't'
                case(29); point_group_name = 't_h'
                case(30); point_group_name = 'o'
                case(31); point_group_name = 't_d'
                case(32); point_group_name = 'o_h'
            case default
                call am_print('ERROR','Schoenflies code for point-group unknown.',flags='E')
                call am_print('pg_code',pg_code)
                stop
            end select
    end function   decode_pointgroup

    function       point_group_schoenflies(ps_id) result(pg_code)
        !>
        !> Refs:
        !>
        !> Applied Group Theory: For Physicists and Chemists. Reissue edition.
        !> Mineola, New York: Dover Publications, 2015. page 20.
        !>
        !> Casas, Ignasi, and Juan J. Pérez. “Modification to Flow Chart to
        !> Determine Point Groups.” Journal of Chemical Education 69, no. 1
        !> (January 1, 1992): 83. doi:10.1021/ed069p83.2.
        !>
        !> Breneman, G. L. “Crystallographic Symmetry Point Group Notation
        !> Flow Chart.” Journal of Chemical Education 64, no. 3 (March 1, 1987):
        !> 216. doi:10.1021/ed064p216.
        !>
        !> Adapted from VASP
        !>
        !>   pg_code --> point_group_name
        !>     1 --> c_1       9 --> c_3      17 --> d_4      25 --> c_6v
        !>     2 --> s_2      10 --> s_6      18 --> c_4v     26 --> d_3h
        !>     3 --> c_2      11 --> d_3      19 --> d_2d     27 --> d_6h
        !>     4 --> c_1h     12 --> c_3v     20 --> d_4h     28 --> t
        !>     5 --> c_2h     13 --> d_3d     21 --> c_6      29 --> t_h
        !>     6 --> d_2      14 --> c_4      22 --> c_3h     30 --> o
        !>     7 --> c_2v     15 --> s_4      23 --> c_6h     31 --> t_d
        !>     8 --> d_2h     16 --> c_4h     24 --> d_6      32 --> o_h
        !>
        implicit none
        !
        integer, intent(in) :: ps_id(:)
        integer :: nsyms
        integer :: ni, nc2, nc3, nc4, nc6, ns2, ns6, ns4, ns3, ne
        integer :: pg_code
        !
        pg_code = 0
        !
        ! trivial cases first
        !
        nsyms = size(ps_id)
        if     (nsyms .eq. 1 ) then; pg_code=1;  return
        elseif (nsyms .eq. 48) then; pg_code=32; return
        elseif (nsyms .eq. 16) then; pg_code=20; return
        elseif (nsyms .eq. 3 ) then; pg_code=9;  return
        endif
        !
        ni=0; nc2=0; nc3=0; nc4=0; nc6=0; ns2=0; ns6=0; ns4=0; ns3=0
        ne  = count(ps_id.eq.1)  ! 'e' identity
        nc2 = count(ps_id.eq.2)  ! 'c_2'
        nc3 = count(ps_id.eq.3)  ! 'c_3'
        nc4 = count(ps_id.eq.4)  ! 'c_4'
        nc6 = count(ps_id.eq.5)  ! 'c_6'
        ni  = count(ps_id.eq.6)  ! 'i' inversion
        ns2 = count(ps_id.eq.7)  ! 's_2'
        ns6 = count(ps_id.eq.8)  ! 's_6'
        ns4 = count(ps_id.eq.9)  ! 's_4'
        ns3 = count(ps_id.eq.10) ! 's_3'
        !
        if (ni.gt.1) then
            call am_print('ERROR','More than one inversion operator found.',flags='E')
            stop
        endif
        !
        if (ne.gt.1) then
            call am_print('ERROR','More than one identity operator found.',flags='E')
        elseif (ne.eq.0) then
            call am_print('ERROR','No identity operator found.',flags='E')
            stop
        endif
        !
        if (nsyms.eq.2)   then
            if (ni.eq.1)  then; pg_code=2;  return; endif
            if (nc2.eq.1) then; pg_code=3;  return; endif
            if (ns2.eq.1) then; pg_code=4;  return; endif
        endif
        if (nsyms.eq.4)   then
            if (ni.eq.1)  then; pg_code=5;  return; endif
            if (nc2.eq.3) then; pg_code=6;  return; endif
            if (ns2.eq.2) then; pg_code=7;  return; endif
            if (nc4.eq.1) then; pg_code=14; return; endif
            if (ns4.eq.2) then; pg_code=15; return; endif
        endif
        if (nsyms.eq.6)   then
            if (ni.eq.1)  then; pg_code=10; return; endif
            if (nc2.eq.3) then; pg_code=11; return; endif
            if (ns2.eq.3) then; pg_code=12; return; endif
            if (nc2.eq.1) then; pg_code=21; return; endif
            if (ns2.eq.1) then; pg_code=22; return; endif
        endif
        if (nsyms.eq.8)  then
            if (ns2.eq.3) then; pg_code=8;  return; endif
            if (ns2.eq.1) then; pg_code=16; return; endif
            if (ns2.eq.0) then; pg_code=17; return; endif
            if (ns2.eq.4) then; pg_code=18; return; endif
            if (ns2.eq.2) then; pg_code=19; return; endif
        endif
        if (nsyms.eq.12)  then
            if (ns2.eq.3) then; pg_code=13; return; endif
            if (ns2.eq.1) then; pg_code=23; return; endif
            if (nc2.eq.7) then; pg_code=24; return; endif
            if (ns2.eq.6) then; pg_code=25; return; endif
            if (ns2.eq.4) then; pg_code=26; return; endif
            if (nc3.eq.8) then; pg_code=28; return; endif
        endif
        if (nsyms.eq.24)  then
            if (nc6.eq.2) then; pg_code=27; return; endif
            if (ni.eq.1)  then; pg_code=29; return; endif
            if (nc4.eq.6) then; pg_code=30; return; endif
            if (ns4.eq.6) then; pg_code=31; return; endif
            endif
        ! if it makes it this far, it means nothing matches. return an error immediately. 
        call am_print('ERROR','Unable to identify point group',flags='E')
            write(*,'(" ... ",a)') 'number of symmetry operators'
            write(*,'(5x,a5,i5)') 'e  ', ne
            write(*,'(5x,a5,i5)') 'c_2', nc2
            write(*,'(5x,a5,i5)') 'c_3', nc3
            write(*,'(5x,a5,i5)') 'c_4', nc4
            write(*,'(5x,a5,i5)') 'c_6', nc6
            write(*,'(5x,a5,i5)') 'i  ', ni
            write(*,'(5x,a5,i5)') 's_2', ns2
            write(*,'(5x,a5,i5)') 's_6', ns6
            write(*,'(5x,a5,i5)') 's_4', ns4
            write(*,'(5x,a5,i5)') 's_3', ns3
            stop
    end function   point_group_schoenflies

    ! needed

    subroutine     put_identity_first(seitz)
        !
        implicit none
        !
        real(dp), intent(inout) :: seitz(:,:,:)
        real(dp) :: id(4,4)
        integer  :: nsyms
        integer  :: i
        !
        id = eye(4)
        nsyms = size(seitz,3)
        !
        search : do i = 1, nsyms
            if (isequal(seitz(:,:,i),id)) then
                seitz(:,:,i) = seitz(:,:,1)
                seitz(:,:,1) = id
                exit search
            endif
        enddo search
        ! 
    end subroutine put_identity_first

end module am_symmetry


