module am_symmetry

    use am_constants
    use am_stdout
    use am_unit_cell
    use am_options
    use am_mkl
    use am_rank_and_sort

    implicit none

    private

    public :: get_multiplication_table
    
    public :: permutation_rep   ! am_irre_cell
    public :: permutation_map   ! am_irre_cell

    public :: member            ! am_symmetry_adapted_tensors
    public :: decode_pointgroup ! am_shell
    public :: nelements         ! am_shell

    public :: ps_frac2cart
    
    type, public :: am_class_symmetry 
        !
        integer :: nsyms
        real(dp), allocatable :: seitz(:,:,:) ! symmetry elements (operate on fractional atomic basis)
        integer, allocatable :: class_id(:)   ! integer assigning each element to a conjugacy class
        integer, allocatable :: ps_id(:)      ! integer which identifies point symmetries (see decode_pointsymmetry)
        integer, allocatable :: multab(:,:)   ! multiplication table
        complex(dp), allocatable :: chartab(:,:) ! character table
        !
        contains
        !
        procedure :: copy
        procedure :: create
        procedure :: write_outfile
        procedure :: write_action_table
        !
    end type am_class_symmetry

    type, public, extends(am_class_symmetry) :: am_class_point_group
        integer   :: pg_id        !> integer which identifies the point group
        contains
        procedure :: get_point_group
        procedure :: get_rotational_group
        procedure :: get_stabilizer_group
        procedure :: get_reversal_group
    end type am_class_point_group

    type, public, extends(am_class_symmetry) :: am_class_space_group
        contains
        procedure :: get_space_group
    end type am_class_space_group

contains

    ! high level routines which operate on sg

    subroutine     get_space_group(sg,uc,opts)
        !
        implicit none
        !
        class(am_class_space_group), intent(inout) :: sg
        class(am_class_unit_cell)  , intent(in) :: uc ! space groups are tabulated in the litearture for conventional cells, but primitive or arbitrary cell works just as well.
        type(am_class_options)     , intent(in) :: opts
        type(am_class_options) :: notalk
        real(dp), allocatable  :: seitz(:,:,:)
        !
        notalk=opts
        notalk%verbosity = 0
        !
        if (opts%verbosity.ge.1) call am_print_title('Determining space group symmetries')
        !
        ! determine space symmetries from atomic basis
        seitz = space_symmetries_from_basis(uc=uc,opts=opts)
        if (opts%verbosity.ge.1) call am_print('number of space group symmetries found',size(seitz,3),' ... ') 
        !
        ! put identity first
        call put_identity_first(seitz=seitz)
        !
        ! create space group instance (puts identity first)
        call sg%create(seitz=seitz)
        !
        ! stdout
        if (opts%verbosity.ge.1) then
            !
            call am_print('space symmetries',sg%nsyms,' ... ')
            !
            call am_print('classes',maxval(sg%class_id))
            !
            write(*,*) maxval(sg%class_id)
            !
            call print_character_table(chartab=sg%chartab, class_nelements=nelements(sg%class_id), class_member=member(sg%class_id), ps_id=sg%ps_id)
            !
        endif
        !
        ! write to write_outfile and to file
        call sg%write_outfile(iopt_uc=uc,iopt_filename=trim('outfile.spacegroup'))
        !
        call sg%write_action_table(uc=uc,fname='outfile.space_group_action',opts=opts)
        !
        contains
            subroutine     put_identity_first(seitz)
                !
                implicit none
                !
                real(dp), intent(inout) :: seitz(:,:,:)
                real(dp), allocatable :: id(:,:)
                integer :: i, nsyms
                !
                nsyms = size(seitz,3)
                id = eye(size(seitz,1))
                do i = 1, nsyms
                    if (all((abs(seitz(:,:,i)-id)).lt.tiny)) then
                        seitz(:,:,i) = seitz(:,:,1)
                        seitz(:,:,1) = id
                    endif
                enddo
                ! 
            end subroutine put_identity_first
    end subroutine get_space_group

    subroutine     get_point_group(pg,uc,sg,opts)
        !
        implicit none
        !
        class(am_class_point_group), intent(inout) :: pg
        type(am_class_space_group) , intent(in) :: sg
        class(am_class_unit_cell)  , intent(in) :: uc ! accepts unit cell, primitive cell, conventional cell.
        type(am_class_options)     , intent(in) :: opts
        real(dp), allocatable  :: seitz(:,:,:)
        type(am_class_options) :: notalk
        !
        notalk=opts
        notalk%verbosity = 0
        !
        if (opts%verbosity.ge.1) call am_print_title('Determining point group symmetries')
        !
        ! create point group instance (puts identity first, inversion second)
        !
        allocate(seitz,source=sg%seitz)
        seitz(1:3,4,:) = 0
        !
        call pg%create(seitz=unique(seitz))
        !
        if (opts%verbosity.ge.1) then
            !
            call am_print('point group',decode_pointgroup(pg%pg_id),' ... ')
            !
            call am_print('point symmetries',pg%nsyms,' ... ')
            !
            call print_character_table(chartab=pg%chartab, class_nelements=nelements(pg%class_id), class_member=member(pg%class_id), ps_id=pg%ps_id)
            !
        endif
        !
        call pg%write_outfile(iopt_uc=uc,iopt_filename='outfile.pointgroup')
        !
        call pg%write_action_table(uc=uc,fname='outfile.point_group_action',opts=opts)
        !
    end subroutine get_point_group

    subroutine     get_rotational_group(rg,uc,pg,opts)
        !
        ! Get point symmetries which are compatible with the atomic basis (essentially keeps
        ! rotational part of space symmetries which are able to map all atoms onto each other)
        !
        ! This is also useful for determining the point symmetries compatible with kpoint mesh.
        !
        ! Name "point group" - NOTE: not really a point group because translational components of
        ! space symmetries were ignored in the genereation process this would be the name of the
        ! group, if the group were symmorphic (in which case it will match the real point group
        ! name)...
        !
        implicit none
        ! subroutine i/o
        class(am_class_point_group),intent(out) :: rg
        class(am_class_unit_cell),  intent(in) :: uc
        type(am_class_point_group), intent(in) :: pg
        type(am_class_options),     intent(in) :: opts
        type(am_class_options) :: notalk
        integer, allocatable   :: indicies(:)
        logical, allocatable   :: mask(:)
        integer :: i
        !
        notalk=opts
        notalk%verbosity = 0
        !
        allocate(mask(pg%nsyms))
        !
        do i = 1, pg%nsyms
            if ( is_symmetry_valid(tau=uc%tau, Z=uc%Z, seitz=pg%seitz, prec=opts%prec)) then
                mask(i) = .true.
            else 
                mask(i) = .false.
            endif
        enddo
        ! create indices
        allocate(indicies(pg%nsyms))
        indicies=[1:pg%nsyms]
        !
        ! create group
        call rg%create(seitz=pg%seitz(:,:,pack(indicies,mask)))
        !
    end subroutine get_rotational_group

    subroutine     get_stabilizer_group(vg,pg,v,opts)
        !
        implicit none
        !
        class(am_class_point_group), intent(out) :: vg ! stabilizer group associated with vector v
        type(am_class_point_group) , intent(in) :: pg
        real(dp)                   , intent(in) :: v(3) ! vector which is stabilized (should be same units as R, which is fractional)
        type(am_class_options)     , intent(in) :: opts
        integer, allocatable :: indicies(:)
        logical, allocatable :: mask(:)
        integer :: i
        !
        allocate(mask(pg%nsyms))
        mask = .false.
        !
        do i = 1, pg%nsyms
            ! mark symmetries which leaves vector strictly invariant (not even modulo a lattice vector)
            if (is_symmetry_valid(tau=reshape(v,[3,1]), Z=[1], seitz=pg%seitz(:,:,i), prec=opts%prec, flags='exact')) then
                mask(i) = .true.
            endif
        enddo
        !
        allocate(indicies(pg%nsyms))
        indicies=[1:pg%nsyms]
        !
        ! create group
        call vg%create(seitz=pg%seitz(:,:,pack(indicies,mask)))
        !
    end subroutine get_stabilizer_group
    
    subroutine     get_reversal_group(revg,sg,v,opts)
        !
        implicit none
        !
        class(am_class_point_group), intent(inout) :: revg ! reversal group associated with vector v, i.e. the space symmetries with bring v to [0,0,0], essentially flipping the atoms at the corners of the vector
        type(am_class_space_group) , intent(in) :: sg ! space group
        real(dp)                   , intent(in) :: v(3)
        type(am_class_options)     , intent(in) :: opts
        integer, allocatable :: indicies(:)
        logical, allocatable :: mask(:)
        integer :: i
        !
        allocate(mask(sg%nsyms))
        mask = .false.
        !
        do i = 1, sg%nsyms
            if (is_symmetry_valid(tau=reshape(v,[3,1]), Z=[1], seitz=sg%seitz(:,:,i), prec=opts%prec, flags='zero')) then
                mask(i) = .true.
            endif
        enddo
        !
        allocate(indicies(sg%nsyms))
        indicies=[1:sg%nsyms]
        !
        ! create group
        call revg%create(seitz=sg%seitz(:,:,pack(indicies,mask)))
        !
    end subroutine get_reversal_group

    subroutine     write_outfile(sg,iopt_uc,iopt_filename)
        !
        ! if iopt_uc is present -> also print symmetries in cartesian coordinates.
        !
        implicit none
        !
        class(am_class_symmetry), intent(in) ::  sg
        type(am_class_unit_cell), intent(in), optional :: iopt_uc
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
        !
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
        !
        write(fid,fmt3,advance='no') centertitle('fractional',width)
        write(fid,*)
        write(fid,fmt3,advance='no') repeat('-',width)
        write(fid,*)
        !
        write(fid,fmt2,advance='no') '#', 'R11', 'R21', 'R31', 'R12', 'R22', 'R32', 'R13', 'R23', 'R33'
        write(fid,fmt4,advance='no') '  | ', 'T1', 'T2', 'T3'
        write(fid,*)
        !
        do i = 1, sg%nsyms
            !
            R = sg%seitz(1:3,1:3,i)
            write(fid,fmt1,advance='no') i, (( trim(print_pretty(R(j,k))), j = 1,3) , k = 1,3)
            !
            T = sg%seitz(1:3,4,i)
            write(fid,fmt5,advance='no') '  | '
            do j = 1,3
                write(fid,fmt6,advance='no') trim(print_pretty(T(j)))
            enddo
            !
            write(fid,*)
            !
        enddo
        !
        ! cartesian
        !
        if (present(iopt_uc)) then
        !
        write(fid,fmt3,advance='no') centertitle('cartesian',width)
        write(fid,*)
        write(fid,fmt3,advance='no') repeat('-',width)
        write(fid,*)
        !
        write(fid,fmt2,advance='no') '#', 'R11', 'R21', 'R31', 'R12', 'R22', 'R32', 'R13', 'R23', 'R33'
        write(fid,fmt4,advance='no') '  | ', 'T1', 'T2', 'T3'
        write(fid,*)
        !
        do i = 1, sg%nsyms
            !
            R = ps_frac2cart(R_frac=sg%seitz(1:3,1:3,i),bas=iopt_uc%bas)
            write(fid,fmt1,advance='no') i, (( trim(print_pretty(R(j,k))), j = 1,3) , k = 1,3)
            !
            T = matmul(iopt_uc%bas,sg%seitz(1:3,4,i))
            write(fid,fmt5,advance='no') '  | '
            do j = 1,3
                write(fid,fmt6,advance='no') trim(print_pretty(T(j)))
            enddo
            !
            write(fid,*)
        enddo
        !
        endif
        !
        ! close file
        !
        if ( present(iopt_filename) ) close(fid)
        !
    end subroutine write_outfile

    subroutine     create(sg,seitz)
        !
        implicit none
        !
        class(am_class_symmetry), intent(out) :: sg
        real(dp)                , intent(in) :: seitz(:,:,:)
        real(dp) :: id(4,4)
        integer  :: i
        !
        ! get number of symmetries
        sg%nsyms = size(seitz,3)
        ! transfer seitz
        allocate(sg%seitz,source=seitz)
        !
        ! put identity first
        id = eye(4)
        search : do i = 1, sg%nsyms
            if (all(abs(sg%seitz(:,:,i)-id).lt.tiny)) then
                sg%seitz(:,:,i) = sg%seitz(:,:,1)
                sg%seitz(:,:,1) = id
                exit search
            endif
        enddo search
        !
        ! identify point symmetries
        allocate(sg%ps_id(sg%nsyms))
        do i = 1, sg%nsyms
            sg%ps_id(i) = ps_schoenflies(R=sg%seitz(1:3,1:3,i))
        enddo
        !
        ! name point group
        select type (sg)
        type is (am_class_point_group)
            sg%pg_id = point_group_schoenflies(sg%ps_id)
        end select
        !
        ! get multiplication table
        sg%multab = get_multiplication_table(seitz=sg%seitz)
        ! determine conjugacy classes (needs identity first, inversion second)
        sg%class_id = get_conjugacy_classes(multab=sg%multab,ps_id=sg%ps_id)
        !
        ! sort symmetries based on parameters
        ! SORT HERE IS BROKEN. NEED TO RECALCULATE MATMUL AND CLASS ID AFTER SORT.
        call sort_symmetries(sg=sg, criterion=sg%seitz(1,4,:), flags='ascend')
        call sort_symmetries(sg=sg, criterion=sg%seitz(2,4,:), flags='ascend')
        call sort_symmetries(sg=sg, criterion=sg%seitz(3,4,:), flags='ascend')
        call sort_symmetries(sg=sg, criterion=real(sg%ps_id,dp), flags='acsend')
        call sort_symmetries(sg=sg, criterion=real(sg%class_id,dp), flags='ascend')
        !
        ! get character table
        sg%chartab = get_character_table(multab=sg%multab,ps_id=sg%ps_id)
        !
        contains
        subroutine     sort_symmetries(sg,criterion,flags)
            ! flags = 'ascend'/'descend'
            implicit none
            !
            class(am_class_symmetry), intent(inout) :: sg
            character(*), intent(in) :: flags
            integer , allocatable :: inds(:)
            integer , allocatable ::rinds(:) ! reverse inds
            real(dp) :: criterion(:)
            integer :: total
            integer :: i, j
            !
            allocate(inds(sg%nsyms))
            !
            call rank(criterion,inds)
            !
            if (index(flags,'descend').ne.0) then
                inds = inds(sg%nsyms:1:-1)
            endif
            !
            if (allocated(sg%seitz   )) sg%seitz    = sg%seitz(:,:,inds)
            if (allocated(sg%ps_id   )) sg%ps_id    = sg%ps_id(inds)
            if (allocated(sg%class_id)) sg%class_id = sg%class_id(inds)
            !
            if (allocated(sg%multab)) then
                !
                allocate(rinds(sg%nsyms))
                rinds(inds) = [1:sg%nsyms]
                !
                do i = 1, sg%nsyms
                do j = 1, sg%nsyms
                    sg%multab(i,j) = rinds(sg%multab(i,j))
                enddo
                enddo
                !
                sg%multab = sg%multab(:,inds)
                sg%multab = sg%multab(inds,:)
                !
            endif
            !
        end subroutine sort_symmetries
    end subroutine create
    
    subroutine     copy(cp,sg)
        !
        implicit none
        !
        class(am_class_symmetry), intent(out) :: cp
        class(am_class_symmetry), intent(in) :: sg
        !
        cp%nsyms = sg%nsyms
        !
        if (allocated(sg%seitz)) allocate(cp%seitz, source=sg%seitz)
        if (allocated(sg%class_id)) allocate(cp%class_id, source=sg%class_id)
        if (allocated(sg%ps_id)) allocate(cp%ps_id, source=sg%ps_id)
        !
    end subroutine copy

    subroutine     write_action_table(sg,uc,fname,opts)
        !
        implicit none
        !
        class(am_class_symmetry), intent(in) :: sg
        class(am_class_unit_cell),intent(in) :: uc
        character(*),             intent(in) :: fname
        type(am_class_options),   intent(in) :: opts
        real(dp), allocatable :: seitz_cart(:,:,:)
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
        allocate(seitz_cart(4,4,sg%nsyms))
        do i = 1, sg%nsyms
            seitz_cart(:,:,i) = seitz_frac2cart(seitz_frac=sg%seitz(:,:,i),bas=uc%bas)
        enddo
        !
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
                write(fid,fmt5,advance='no') sg%class_id(i)
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
                    write(fid,fmt4,advance='no') seitz_cart(i,j,k)
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
                aa = rot2axis_angle(R=seitz_cart(1:3,1:3,i))
                write(fid,fmt4,advance='no') aa(j)
            enddo
            write(fid,*)
            enddo
            deallocate(xyz)
            ! POINT-SYMMETRY ROTATION ANGLE
            write(fid,'(5x)',advance='no')
            write(fid,fmt1,advance='no') 'R_th'
            do i = 1, sg%nsyms
                aa = rot2axis_angle(R=seitz_cart(1:3,1:3,i))
                write(fid,fmt4,advance='no') aa(4)*180.0_dp/pi
            enddo
            write(fid,*)
            ! POINT-SYMMETRY ROTOINVERSION VS INVERSION
            write(fid,'(5x)',advance='no')
            write(fid,fmt1,advance='no') 'inv'
            do i = 1, sg%nsyms
                write(fid,fmt6,advance='no') (abs(det(seitz_cart(1:3,1:3,i))+1).lt.tiny)
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
                    write(fid,fmt4,advance='no') seitz_cart(j,4,i)
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
                    R=inv(seitz_cart(1:3,1:3,i))
                    ! if statement are for the cases in which output needs to be trimmed
                    if (count(abs(R(j,:)).gt.tiny).gt.3) then
                        buffer=' ... ' 
                    else
                        buffer=''
                        do k = 1,3
                        if (abs(R(j,k)).gt.tiny) then
                            ! s is only the sign
                            s = trim(int2char(nint(R(j,k)),'SP'))
                            !
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
                    R = rot2O3(l=1,R=seitz_cart(1:3,1:3,i))
                    !
                    if (any(abs(R-seitz_cart(1:3,1:3,i)).gt.tiny)) then
                        call am_print('seitz_cart(1:3,1:3,i)',seitz_cart(1:3,1:3,i))
                        call am_print('R',R)
                        call am_print('R',R-seitz_cart(1:3,1:3,i))
                    endif
                    ! if statement are for the cases in which output needs to be trimmed
                    if (count(abs(R(j,:)).gt.tiny).gt.3) then
                        buffer=' ... ' 
                    else
                        buffer=''
                        do k = 1, size(xyz)
                            if (abs(R(j,k)).gt.tiny) then
                                ! s is only the sign
                                s = trim(int2char(nint(R(j,k)),'SP'))
                                !
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
                    R = rot2O3(l=2,R=seitz_cart(1:3,1:3,i))
                    !
                    if (count(abs(R(j,:)).gt.tiny).gt.3) then
                        buffer=' ... ' 
                    else
                        buffer=''
                        do k = 1,5
                            if (abs(R(j,k)).gt.tiny) then
                                ! s is only the sign
                                s = trim(int2char(nint(R(j,k)),'SP'))
                                !
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
                    R = rot2O3(l=3,R=seitz_cart(1:3,1:3,i))
                    !
                    if (count(abs(R(j,:)).gt.tiny).gt.3) then
                        buffer=' ... ' 
                    else
                        buffer=''
                        do k = 1,7
                            if (abs(R(j,k)).gt.tiny) then
                                ! s is only the sign
                                s = trim(int2char(nint(R(j,k)),'SP'))
                                !
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
                P  = permutation_rep(seitz=sg%seitz,tau=uc%tau,flags='',prec=opts%prec)
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
        close(fid)
    end subroutine write_action_table

    ! functions which operate on R(3,3) in fractional coordinates

    function       seitz_frac2cart(seitz_frac,bas) result(seitz_cart)
        !
        implicit none
        !
        real(dp), intent(in) :: seitz_frac(4,4)
        real(dp), intent(in) :: bas(3,3)
        real(dp) :: seitz_cart(4,4)
        !
        seitz_cart(1:3,1:3) = matmul(bas,matmul(seitz_frac(1:3,1:3),inv(bas)))
        seitz_cart(1:3,4)   = matmul(bas,seitz_frac(1:3,4))
        seitz_cart(4,4)     = 1.0_dp
        seitz_cart(4,1:3)   = 0.0_dp
        !
    end function   seitz_frac2cart

    function       seitz_cart2frac(seitz_cart,bas) result(seitz_frac)
        !
        implicit none
        !
        real(dp), intent(in) :: seitz_cart(4,4)
        real(dp), intent(in) :: bas(3,3)
        real(dp) :: seitz_frac(4,4)
        real(dp) :: recbas(3,3)
        !
        recbas = inv(bas)
        !
        seitz_frac(1:3,1:3) = matmul(recbas,matmul(seitz_cart(1:3,1:3),bas))
        seitz_frac(1:3,4)   = matmul(recbas,seitz_cart(1:3,4))
        seitz_frac(4,4)     = 1.0_dp
        seitz_frac(4,1:3)   = 0.0_dp
        !
    end function   seitz_cart2frac

    function       ps_frac2cart(R_frac,bas) result(R)
        !
        ! The values which are possible are: cos(pi/n), sin(pi/n) for n = 1, 2, 3, 6 Symmetry and
        ! Condensed Matter Physics: A Computational Approach. 1 edition. Cambridge, UK ; New York:
        ! Cambridge University Press, 2008. page 275.
        !
        !           n  =      1         2         3         6
        !    sin(pi/n) =   0.0000    1.0000    0.8660    0.5000
        !    cos(pi/n) =  -1.0000    0.0000    0.5000    0.8660
        !
        implicit none
        !
        real(dp), intent(in) :: R_frac(3,3)
        real(dp), intent(in) :: bas(3,3)
        real(dp) :: R(3,3)
        ! real(dp) :: wdv(7) ! used to correct rouding errors
        ! integer :: i, j, k
        ! !
        ! wdv(1)=+0.0_dp
        ! wdv(2)=+1.0_dp
        ! wdv(3)=-1.0_dp
        ! wdv(4)=+0.86602540378443862 !  sqrt(3.0_dp)/2.0_dp
        ! wdv(5)=-0.86602540378443862 ! -sqrt(3.0_dp)/2.0_dp
        ! wdv(6)=+0.5_dp
        ! wdv(7)=-0.5_dp
        !
        R = matmul(bas,matmul(R_frac,inv(bas)))
        !
        ! do i = 1,3
        ! do j = 1,3
        !     do k = 1,7
        !     if (abs(wdv(k)-R(i,j)).lt.tiny) then
        !         R(i,j) = wdv(k)
        !         exit
        !     endif
        !     enddo
        ! enddo
        ! enddo
    end function   ps_frac2cart

    function       ps_cart2frac(R,bas) result(R_frac)
        !
        implicit none
        !
        real(dp), intent(in) :: R(3,3)
        real(dp), intent(in) :: bas(3,3)
        real(dp) :: R_frac(3,3)
        !
        R_frac = matmul(inv(bas),matmul(R,bas))
        !
        ! correct rounding errors
        ! R_frac = real(nint(R_frac),dp)
        !
    end function   ps_cart2frac

    ! decoding/naming functions 

    function       ps_schoenflies(R) result(ps_id)
        !>
        !> Point symmetries in fractional coordinates so that they are nice integers which can be easily classified.
        !>
        !>    element:         e    i  c_2  c_3  c_4  c_6  s_2  s_6  s_4  s_3
        !>    trace:          +3   -3   -1    0   +1   +2   +1    0   -1   -2
        !>    determinant:    +1   -1   +1   +1   +1   +1   -1   -1   -1   -1
        !>
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
        if (ps_id .eq. 0) then
            call am_print('ERROR','Unable to identify point symmetry.',flags='E')
            call am_print('R (frac)',R)
            call am_print('tr (frac)',tr)
            call am_print('det (frac)',d)
            stop
        endif
        !
    end function   ps_schoenflies

    function       decode_pointsymmetry(ps_id) result(schoenflies)
        !
        implicit none
        !
        integer, intent(in) :: ps_id
        character(len=string_length_schoenflies) :: schoenflies
        !
        select case (ps_id)
        case(1);  schoenflies = 'e'
        case(2);  schoenflies = 'c_2'
        case(3);  schoenflies = 'c_3'
        case(4);  schoenflies = 'c_4'
        case(5);  schoenflies = 'c_6'
        case(6);  schoenflies = 'i'
        case(7);  schoenflies = 's_2'
        case(8);  schoenflies = 's_6'
        case(9);  schoenflies = 's_4'
        case(10); schoenflies = 's_3'
        case default
            call am_print('ERROR','Schoenflies code unknown.',flags='E')
            stop
        end select
    end function   decode_pointsymmetry

    function       decode_pointgroup(pg_id) result(point_group_name)
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
        integer, intent(in) :: pg_id
        character(string_length_schoenflies) :: point_group_name
        !
        select case (pg_id)
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
                call am_print('pg_id',pg_id)
                stop
            end select
    end function   decode_pointgroup

    function       point_group_schoenflies(ps_id) result(pg_id)
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
        !>   pg_id --> point_group_name
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
        integer :: pg_id
        !
        pg_id = 0
        !
        ! trivial cases first
        !
        nsyms = size(ps_id)
        if     (nsyms .eq. 1 ) then; pg_id=1;  return
        elseif (nsyms .eq. 48) then; pg_id=32; return
        elseif (nsyms .eq. 16) then; pg_id=20; return
        elseif (nsyms .eq. 3 ) then; pg_id=9;  return
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
            if (ni.eq.1)  then; pg_id=2;  return; endif
            if (nc2.eq.1) then; pg_id=3;  return; endif
            if (ns2.eq.1) then; pg_id=4;  return; endif
        endif
        if (nsyms.eq.4)   then
            if (ni.eq.1)  then; pg_id=5;  return; endif
            if (nc2.eq.3) then; pg_id=6;  return; endif
            if (ns2.eq.2) then; pg_id=7;  return; endif
            if (nc4.eq.1) then; pg_id=14; return; endif
            if (ns4.eq.2) then; pg_id=15; return; endif
        endif
        if (nsyms.eq.6)   then
            if (ni.eq.1)  then; pg_id=10; return; endif
            if (nc2.eq.3) then; pg_id=11; return; endif
            if (ns2.eq.3) then; pg_id=12; return; endif
            if (nc2.eq.1) then; pg_id=21; return; endif
            if (ns2.eq.1) then; pg_id=22; return; endif
        endif
        if (nsyms.eq.8)  then
            if (ns2.eq.3) then; pg_id=8;  return; endif
            if (ns2.eq.1) then; pg_id=16; return; endif
            if (ns2.eq.0) then; pg_id=17; return; endif
            if (ns2.eq.4) then; pg_id=18; return; endif
            if (ns2.eq.2) then; pg_id=19; return; endif
        endif
        if (nsyms.eq.12)  then
            if (ns2.eq.3) then; pg_id=13; return; endif
            if (ns2.eq.1) then; pg_id=23; return; endif
            if (nc2.eq.7) then; pg_id=24; return; endif
            if (ns2.eq.6) then; pg_id=25; return; endif
            if (ns2.eq.4) then; pg_id=26; return; endif
            if (nc3.eq.8) then; pg_id=28; return; endif
        endif
        if (nsyms.eq.24)  then
            if (nc6.eq.2) then; pg_id=27; return; endif
            if (ni.eq.1)  then; pg_id=29; return; endif
            if (nc4.eq.6) then; pg_id=30; return; endif
            if (ns4.eq.6) then; pg_id=31; return; endif
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

    ! permutation rep

    function       permutation_rep(seitz,tau,prec,flags) result(rep)
        !
        ! find permutation representation; i.e. which atoms are connected by space symmetry oprations R, T.
        ! also works to find which kpoint are connected by point group operations
        !
        implicit none
        !
        real(dp), intent(in) :: seitz(:,:,:)
        real(dp), intent(in) :: tau(:,:)
        real(dp), intent(in) :: prec
        character(*), intent(in) :: flags
        integer, allocatable :: rep(:,:,:)
        real(dp) :: tau_rot(3)
        integer :: i,j,k
        logical :: found
        integer :: ntaus ! number of tau points tau(1:3,ntaus)
        integer :: nsyms
        !
        nsyms = size(seitz,3)
        !
        ntaus = size(tau,2)
        !
        allocate(rep(ntaus,ntaus,nsyms))
        rep = 0
        !
        do i = 1, nsyms
            ! determine the permutations of atomic indicies which results from each space symmetry operation
            do j = 1,ntaus
                found = .false.
                ! apply rotational component
                tau_rot = matmul(seitz(1:3,1:3,i),tau(:,j))
                ! apply translational component
                tau_rot = tau_rot + seitz(1:3,4,i)
                ! reduce rotated+translated point to unit cell
                if (index(flags,'relax_pbc').eq.0) then
                    tau_rot = modulo(tau_rot+prec,1.0_dp)-prec
                endif
                ! find matching atom
                search : do k = 1, ntaus
                    if (all(abs(tau_rot-tau(:,k)).lt.prec)) then
                    rep(j,k,i) = 1
                    found = .true.
                    exit search ! break loop
                    endif
                enddo search
                ! if "relax_pbc" is present (i.e. periodic boundary conditions are relax), perform check
                ! to ensure atoms must permute onto each other
                if (index(flags,'relax_pbc').eq.0) then
                if (found.eq..false.) then
                    call am_print('ERROR','Unable to find matching atom.',flags='E')
                    call am_print('tau (all atoms)',transpose(tau))
                    call am_print('tau',tau(:,j))
                    call am_print('R',seitz(1:3,1:3,i))
                    call am_print('T',seitz(1:3,4,i))
                    call am_print('tau_rot',tau_rot)
                    stop
                endif
                endif
            enddo
        enddo
        !
        ! check that each column and row of the rep sums to 1; i.e. rep(:,:,i) is orthonormal for all i.
        !
        ! if "relax_pbc" is present (i.e. periodic boundary conditions are relax), perform check
        ! to ensure atoms must permute onto each other
        if (index(flags,'relax_pbc').eq.0) then
        do i = 1, nsyms
        do j = 1, ntaus
            !
            if (sum(rep(:,j,i)).ne.1) then
               call am_print('ERROR','Permutation matrix has a column which does not sum to 1.')
               call am_print('i',i)
               call am_print_sparse('spy(P_i)',rep(:,:,i))
               call am_print('rep',rep(:,:,i))
               stop
            endif
            !
            if (sum(rep(j,:,i)).ne.1) then
               call am_print('ERROR','Permutation matrix has a row which does not sum to 1.')
               call am_print('i',i)
               call am_print_sparse('spy(P_i)',rep(:,:,i))
               call am_print('rep',rep(:,:,i))
               stop
            endif
            !
        enddo
        enddo
        endif
       !
    end function   permutation_rep
    
    function       permutation_map(P) result(PM)
        !
        ! PM(uc%natoms,sg%nsyms) permutation map; shows how atoms are permuted by each space symmetry operation
        !
        implicit none
        !
        integer, allocatable, intent(in) :: P(:,:,:)
        integer, allocatable :: PM(:,:)
        integer, allocatable :: ind(:)
        integer :: natoms,nsyms
        integer :: i
        !
        natoms = size(P,1)
        nsyms  = size(P,3)
        !
        ind = [1:natoms]
        !
        allocate(PM(natoms,nsyms))
        PM=0
        !
        do i = 1, nsyms
            PM(:,i) = matmul(P(:,:,i),ind)
        enddo
        !
    end function   permutation_map

    ! identifier functions which operate on identifiers

    function       member(id) result(class_member)
        !
        implicit none
        !
        integer, intent(in) :: id(:)
        integer, allocatable :: class_nelements(:) ! number of elements in each class
        integer, allocatable :: class_member(:,:) ! members(nclass,maxval(class_nelements))
        integer :: i, j, k
        integer :: nsyms
        integer :: nclasses
        !
        nsyms = size(id,1)
        !
        nclasses = maxval(id)
        !
        class_nelements = nelements(id)
        !
        allocate(class_member(nclasses,maxval(class_nelements)))
        class_member = 0
        do i = 1, nclasses
            k=0
            do j = 1, nsyms
                if (id(j).eq.i) then
                    k=k+1
                    class_member(i,k) = j
                endif
            enddo
        enddo
    end function   member

    function       nelements(id) result(class_nelements)
        !
        implicit none
        !
        integer, intent(in)  :: id(:)
        integer, allocatable :: class_nelements(:)
        integer :: nclasses
        integer :: i
        !
        nclasses = maxval(id)
        !
        allocate(class_nelements(nclasses))
        !
        do i = 1, nclasses
            class_nelements(i) = count(i.eq.id)
        enddo
    end function   nelements

    ! multiplication chartab

    function       get_multiplication_table(seitz,flags) result(multab)
        !
        ! flag options:
        !  reg    -  sorts rows to put identity along diagonals (useful for constructing regular representation)
        !
        implicit none
        !
        real(dp), intent(in) :: seitz(:,:,:) ! list of 2D reps...
        character(*), intent(in), optional :: flags
        integer, allocatable :: multab(:,:)
        integer, allocatable :: sortmat(:,:)
        integer, allocatable :: indices(:)
        real(dp) :: W(size(seitz,1),size(seitz,2)) ! workspace
        integer  :: n
        integer  :: i, j
        !
        ! before doing anything, confirm that first element is the identity
        if (any(abs(seitz(:,:,1)-eye(size(seitz,1))).gt.tiny)) then
            call am_print('ERROR','First element of seitz is not the identity.')
            call am_print('first element of seitz',seitz(:,:,1))
            stop
        endif
        !
        n=size(seitz,3)
        !
        allocate(multab(n,n))
        !
        do i = 1, n
        do j = 1, n
            ! multiply the two seitz operators
            W = matmul(seitz(:,:,i),seitz(:,:,j))
            ! if 'seitz' flagged, reduce translational part to primitive cell 
            W(1:3,4) = modulo(W(1:3,4)+tiny,1.0_dp)-tiny
            ! get matching lment
            multab(i,j) = get_matching_element_index_in_list(list=seitz,elem=W)
        enddo
        enddo
        !
        ! if 'reg' flagged, re-order rows so that identities are along diagonal
        if (present(flags)) then
        if (index(flags,'reg').ne.0) then
            !
            allocate(sortmat(n,n))
            allocate(indices(n))
            indices = [1:n]
            sortmat = 0
            where (multab.eq.1) sortmat(:,:) = 1
            indices = matmul(sortmat(:,:),indices)
            multab = multab(indices,:)
        endif
        endif
        !
        ! quick consitency check. not comprehensive
        do i = 1, n
            if (sum(multab(:,i)).ne.sum(multab(:,1))) then
                call am_print('ERROR','multiplication chartab sums along rows/columns are not consistent.')
                stop
            endif
        enddo
        !
        contains
        function       get_matching_element_index_in_list(list,elem) result(i)
            !
            implicit none
            real(dp), intent(in) :: list(:,:,:)
            real(dp), intent(in) :: elem(:,:)
            integer :: i
            !
            do i = 1,size(list,3)
                if (all(abs(list(:,:,i)-elem).lt.tiny)) return
            enddo
        end function   get_matching_element_index_in_list
    end function   get_multiplication_table

    pure function  get_inverse_indices(multab) result(inv_ind)
        ! returns the inv_ind of the inverse elements given the multiplication chartab.
        ! requires identity as first element
        implicit none
        !
        integer, intent(in) :: multab(:,:)
        integer, allocatable :: inv_ind(:)
        integer, allocatable :: P(:,:)
        !
        inv_ind = [1:size(multab,1)]
        !
        allocate(P,source=multab)
        where (P.ne.1) P = 0
        !
        inv_ind = matmul(P,inv_ind)
        !
    end function   get_inverse_indices

    function       get_cyclic_order(multab) result(order)
        !
        implicit none
        !
        integer, intent(in) :: multab(:,:)
        integer, allocatable :: order(:)
        integer :: power
        integer :: i, j
        integer :: nsyms
        !
        nsyms = size(multab,2)
        !
        allocate(order(nsyms))
        !
        do i = 1, nsyms
            ! initialize
            j = 0
            power = i
            ! apply operation until identity (1) is returned
            do while (power.ne.1)
                !
                j=j+1
                power = multab(i,power)
                !
                if (j.eq.50) stop 'Error: unable to determine symmetry order.'
            enddo
            !
            order(i) = j
        enddo
        !
    end function   get_cyclic_order

    function       get_conjugacy_classes(multab,ps_id) result(class_id)
        !
        ! for AX = XB, if elements A and B are conjugate pairs for some other element X in the group, then they are in the same class
        !
        implicit none
        !
        integer, intent(in)  :: multab(:,:) ! multiplication chartab
        integer, intent(in)  :: ps_id(:)    ! sort classes basd on point symmetry id: puts identity first and orders classes containing inversion nicely in the second half
        integer, allocatable :: inv_ind(:)  ! inv_ind(nsyms) index of the inverse of symmetry element i
        integer, allocatable :: conjugates(:)
        integer, allocatable :: class_id(:)
        integer :: i, j, k
        integer :: nsyms
        integer :: nclasses
        !
        nsyms = size(multab,2)
        !
        ! allocate space for conjugacy class
        allocate(class_id(nsyms))
        class_id = 0
        !
        ! get element inverse
        inv_ind = get_inverse_indices(multab)
        !
        ! allocate space for conjugate elements
        allocate(conjugates(nsyms))
        !
        ! determine conjugacy classes
        k = 0
        do i = 1, nsyms
        if (class_id(i).eq.0) then
            k=k+1
            ! conjugate each element with all other group elements
            ! A = X(j) * B * X(j)^-1
            do j = 1, nsyms
                conjugates(j) = multab(j,multab(i,inv_ind(j)))
            enddo
            ! for each subgroup element created by conjugation find the corresponding index of the element in the group
            ! in order to save the class class_id number
            do j = 1, nsyms
                class_id( conjugates(j) ) = k
            enddo
            !
        endif
        enddo
        nclasses = k
        !
        ! relabel based on how many elements each class has
        call relabel_based_on_occurances(class_id)
        ! relabel classes based on point symmetry id of representative class element
        call relabel_based_on_ps_id(class_id,ps_id)
        !
        ! make sure the identity is in the first class
        ! class_id(i=1) is the class of the first element (the identity); swap it's location with whaterver elements are in the first class
        where (class_id.eq.1) class_id = class_id(1)
        class_id(1)=1
        !
        ! check that classes are disjoint and complete
        if (any(class_id.eq.0)) then
            call am_print('ERROR','Not every element in the group has been asigned a conjugacy class.',flags='E')
            stop
        endif
        !
        contains
        subroutine     relabel_based_on_occurances(class_id)
            !
            implicit none
            !
            integer , intent(inout) :: class_id(:)
            integer , allocatable :: A_sorted(:)
            integer , allocatable :: occurances(:) 
            integer , allocatable :: reverse_sort(:)
            integer , allocatable :: sorted_indices(:)
            integer , allocatable :: list_relabled(:)
            integer :: nsyms, i, j
            !
            nsyms = size(class_id,1)
            !
            allocate(occurances(nsyms))
            do i = 1, nsyms
                occurances(i) = count(class_id(i).eq.class_id)
            enddo
            !
            ! this quick and dirty procedure lifts degeneracies
            !
            allocate(sorted_indices(nsyms))
            allocate(A_sorted(nsyms))
            call rank((1+maxval(class_id))*occurances+class_id,sorted_indices)
            A_sorted = class_id(sorted_indices)
            !
            allocate(list_relabled(nsyms))
            list_relabled = 0
            !
            j=1
            list_relabled(1) = 1
            do i = 2, nsyms
                if (A_sorted(i).ne.A_sorted(i-1)) then
                    j=j+1
                endif
                list_relabled(i) = j
            enddo
            !
            ! return everything to the original order at call
            !
            allocate(reverse_sort(nsyms))
            reverse_sort(sorted_indices)=[1:nsyms]
            class_id=list_relabled(reverse_sort)
        end subroutine relabel_based_on_occurances
        subroutine     relabel_based_on_ps_id(class_id,ps_id)
            !
            implicit none
            !
            integer, intent(inout) :: class_id(:)
            integer, intent(in)  :: ps_id(:)
            integer :: nclasses
            integer, allocatable :: s(:)
            integer, allocatable :: inds(:)
            integer, allocatable :: rinds(:)
            integer, allocatable :: class_member(:,:)
            integer :: nsyms, i, j
            !
            nsyms = size(class_id)
            !
            nclasses = maxval(class_id)
            !
            allocate(inds,source=[1:nclasses])
            !
            ! identify class members 
            ! members(nclass,maxval(class_nelements))
            class_member = member(class_id)
            !
            ! for each class, identify the representative element
            allocate(s(nclasses))
            do i = 1, nclasses
                s(i) = ps_id( class_member(i,1) )
            enddo
            ! 
            ! sort classes based on representative element id
            call rank(s,inds)
            !
            !
            allocate(rinds(nclasses))
            rinds(inds) = [1:nclasses] 
            !
            ! at this point inds contains the new labeling for class_id
            ! [1:nclasses] => inds; relabel according to new labeling
            do i = 1, nsyms
                class_id(i) = rinds( class_id(i) )
            enddo
            !
        end subroutine relabel_based_on_ps_id
    end function   get_conjugacy_classes

    function       get_character_table(multab,ps_id) result(chartab)
        !
        implicit none
        !
        integer, intent(in) :: multab(:,:)
        integer, intent(in) :: ps_id(:)
        integer :: i, j, k
        integer :: nsyms
        integer :: nclasses
        integer :: nirreps  ! just to make things clearer; always equal to nclasses
        integer, allocatable :: class_nelements(:) ! number of elements in each class
        integer, allocatable :: class_member(:,:) ! members(nclass,maxval(class_nelements))
        integer, allocatable :: H(:,:,:) ! class multiplication 
        complex(dp), allocatable :: chartab(:,:) ! class constant chartab which becomes character chartab (can be complex!)
        integer, allocatable :: irrep_dim(:) ! sym dimensions
        integer, allocatable :: indices(:)
        integer, allocatable :: class_id(:)
        real(dp) :: wrk
        !
        nsyms = size(multab,2)
        !
        ! id(nsyms) identify conjugacy classes
        class_id = get_conjugacy_classes(multab=multab,ps_id=ps_id)
        !
        ! get number of classes
        nclasses = maxval(class_id)
        !
        ! class_nelements(nclasses) get number of elements in each class
        class_nelements = nelements(class_id)
        !
        ! members(nclass,maxval(class_nelements)) record indicies of each class_member element for each class 
        ! (use the first class_member of class as class representative)
        class_member = member(class_id)
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
    end function   get_character_table

    subroutine     print_character_table(chartab,class_nelements,class_member,ps_id)
        !
        implicit none
        !
        complex(dp), intent(in) :: chartab(:,:) ! class constant chartab which becomes character chartab (can be complex!)
        integer, intent(in) :: class_nelements(:) ! number of elements in each class
        integer, intent(in) :: class_member(:,:) ! members(nclass,maxval(class_nelements))
        integer, intent(in) :: ps_id(:)
        integer :: i, j, k
        integer :: nirreps, nclasses
        ! complex to symbol
        integer :: kk, k_exp
        complex(dp), allocatable :: s(:)
        complex(dp) :: s_exp
        complex(dp) :: Z
        real(dp) :: Zr,Zi
        character(:), allocatable :: str
        integer :: char_start
        logical :: strmatch
        ! labels
        character(:), allocatable :: irrep_label(:)
        integer :: class_containing_identity
        integer :: class_containing_inversion
        integer :: class_containing_Cn
        ! print formats
        character(50) :: fmt1
        character(50) :: fmt2
        character(50) :: fmt3
        character(50) :: fmt4
        character(50) :: fmt5
        !
        ! they are always identical...
        nirreps = size(chartab,1)
        nclasses = size(chartab,2)
        !
        ! left-side headers
        fmt2 = '(5x,a10)'
        fmt4 = '(5x,a6,i4)'
        ! chartab
        fmt1 = '(i6)'
        fmt3 = '(a6)'
        fmt5 = '(f6.2)'
        !
        !
        ! create irrep label
        allocate(character(7) :: irrep_label(nirreps))
        !
        class_containing_identity = get_class_with_ps_id(nclasses=nclasses, id=1, ps_id=ps_id, class_nelements=class_nelements, class_member=class_member)
        class_containing_inversion= get_class_with_ps_id(nclasses=nclasses, id=6, ps_id=ps_id, class_nelements=class_nelements, class_member=class_member)
        class_containing_Cn       = get_class_with_ps_id(nclasses=nclasses, id=1, ps_id=ps_id, class_nelements=class_nelements, class_member=class_member)
        do k = 1, nirreps
            ! 0) nullify irrep label
            irrep_label(k) = ''
            ! 1) find class containing identity (ps_id=1)
            i = class_containing_identity
            select case (nint(real(chartab(k,i))))
                case (1)
                    irrep_label(k) = trim(irrep_label(k))//' A'
                case (2)
                    irrep_label(k) = trim(irrep_label(k))//' E'
                case (3)
                    irrep_label(k) = trim(irrep_label(k))//' T'
            end select
            ! 2) find class containing inversion (ps_id=6)
            i = class_containing_inversion
            if (i.ne.0) then
            select case (sign(1,nint(real(chartab(k,i)))))
                case (+1); irrep_label(k) = trim(irrep_label(k))//'_g'
                case (-1); irrep_label(k) = trim(irrep_label(k))//'_u'
            end select
            endif
            !
        enddo
        !
        !
        !
        write(*,fmt2,advance='no') 'class'
        do i = 1, nclasses
            write(*,fmt1,advance='no') i
        enddo
        write(*,*)
        !
        write(*,fmt2,advance='no') 'elements'
        do i = 1, nclasses
            write(*,fmt1,advance='no') class_nelements(i)
        enddo
        write(*,*)
        !
        write(*,fmt2,advance='no') 'class rep'
        do i = 1, nclasses
            write(*,fmt3,advance='no') trim(decode_pointsymmetry(ps_id(class_member(i,1))))
        enddo
        write(*,*)
        !
        write(*,fmt2,advance='no') repeat('-',10)
        do i = 1, nclasses
            write(*,fmt3,advance='no') repeat('-',6)
        enddo
        write(*,*)
        !
        allocate(s(30)) ! value of symbolic output, 30 possible options
        allocate(character(4)::str)
        char_start = 96 ! 97 = a
        k=0
        !
        do i = 1, nirreps
            !
            ! write(*,fmt4,advance='no') 'irrep', i
            write(*,'(5x,i2,a8)',advance='no') i, irrep_label(i)
            !
            do j = 1, nclasses
                !
                strmatch = .false.
                !
                Z = chartab(i,j)
                Zr= real(Z)
                Zi= aimag(Z)
                !
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
                        do while ( .not. equals(s_exp,cmplx(1,0,dp)) )
                            ! do a full loop. complex numbers form a cyclic abelian group.
                            ! exponentiate it until it loops back to one, the identity
                            k_exp = k_exp+1
                            s_exp = s(k)**k_exp
                            ! check positive
                            if ( equals(Z,s_exp) ) then
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
                            if ( equals(Z,-s_exp) ) then
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
                if (strmatch) then
                    write(*,fmt3,advance='no') str
                else
                    write(*,fmt5,advance='no') real(chartab(i,j))
                endif
                !
            enddo
            write(*,*)
        enddo
        !
        if (k.ne.0) then
            !
            write(*,*)
            !
            do i = 1, k
                write(*,'(5x,a,a)') char(char_start+k)//' = ', trim(dbl2char(real(s(i)),8))//trim(dbl2charSP(aimag(s(i)),9))//'i'
            enddo
            !
        endif
        !
        contains
            pure function get_class_with_ps_id(nclasses, id, ps_id, class_nelements, class_member) result(i)
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
            end function  get_class_with_ps_id
    end subroutine print_character_table

end module am_symmetry


