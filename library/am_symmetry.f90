module am_symmetry
    !
    use am_constants
    use am_helpers
    use am_unit_cell
    use am_options
    !
    implicit none
    !
    private
    public :: am_class_symmetry
    public :: determine_symmetry

    type am_class_symmetry 
        integer :: nsyms                  !> number of point symmetries
        real(dp), allocatable :: R(:,:,:) !> symmetry elements (operate on cartesian atomic basis)
        real(dp), allocatable :: T(:,:)   !> symmetry elements (operate on cartesian atomic basis)
        !
        integer , allocatable :: schoenflies_code(:) ! look at schoenflies_decode to understand what each integer corresponds to
        character(string_length_schoenflies) , allocatable :: schoenflies(:) ! decoded string
        !
    contains
        ! procedure :: determine_symmetry
        procedure :: point_group
        procedure :: space_group
        procedure :: stdout
    end type

contains

    !
    ! functions which operate on R(3,3) in cartesian coordinates
    !

    pure function  point_symmetry_convert_to_cartesian(R_frac,bas) result(R)
        !
        use am_unit_cell, only : reciprocal_basis
        !
        implicit none
        !
        real(dp), intent(in) :: R_frac(3,3)
        real(dp), intent(in) :: bas(3,3)
        real(dp) :: R(3,3)
        !
        R = matmul(bas,matmul(R_frac,reciprocal_basis(bas)))
        !
    end function   point_symmetry_convert_to_cartesian

    pure function  point_symmetry_convert_to_fractional(R,bas) result(R_frac)
        !
        use am_unit_cell, only : reciprocal_basis
        !
        implicit none
        !
        real(dp), intent(in) :: R(3,3)
        real(dp), intent(in) :: bas(3,3)
        real(dp) :: R_frac(3,3)
        !
        R_frac = matmul(reciprocal_basis(bas),matmul(R,bas))
        !
    end function   point_symmetry_convert_to_fractional

    pure function  point_symmetry_trace(R) result(trace)
        !
        implicit none
        !
        real(dp), intent(in) :: R(3,3)
        real(dp) :: trace
        !
        trace = R(1,1)+R(2,2)+R(3,3)
        !
    end function   point_symmetry_trace

    pure function  point_symmetry_determinant(R) result(determinant)
        !
        implicit none
        !
        real(dp), intent(in) :: R(3,3)
        real(dp) :: determinant
        !
        determinant=R(1,1)*(R(2,2)*R(3,3)-R(2,3)*R(3,2))+&
                  & R(1,2)*(R(2,3)*R(3,1)-R(2,1)*R(3,3))+&
                  & R(1,3)*(R(2,1)*R(3,2)-R(2,2)*R(3,1))
        !
    end function   point_symmetry_determinant

    function       point_symmetry_schoenflies(R) result(schoenflies_code)
        !>
        !> Point symmetries are input into this function in cartesian coordinates
        !> and converted inside to fractional coordinates so that they are nice
        !> integers which can be easily classified.
        !>
        !>    element:         e    i  c_2  c_3  c_4  c_6  s_2  s_6  s_4  s_3
        !>    trace:          +3   -3   -1    0   +1   +2   +1    0   -1   -2
        !>    determinant:    +1   -1   +1   +1   +1   +1   -1   -1   -1   -1
        !>
        implicit none
        !
        real(dp), intent(in) :: R(3,3)
        real(dp) :: tr
        real(dp) :: det
        integer :: schoenflies_code
        !
        ! get trace and determinant (fractional)
        tr  = point_symmetry_trace(R)
        det = point_symmetry_determinant(R)
        !
        ! The Mathematical Theory of Symmetry in Solids: Representation Theory for
        ! Point Groups and Space Groups. 1 edition. Oxford?: New York: Oxford University
        ! Press, 2010. page 138, table 3.8.
        !
        schoenflies_code = 0
        if     ( (abs(tr - 3).lt.tiny) .and. (abs(det - 1).lt.tiny) ) then; schoenflies_code = 1  ! 'e'
        elseif ( (abs(tr + 1).lt.tiny) .and. (abs(det - 1).lt.tiny) ) then; schoenflies_code = 2  ! 'c_2'
        elseif ( (abs(tr - 0).lt.tiny) .and. (abs(det - 1).lt.tiny) ) then; schoenflies_code = 3  ! 'c_3'
        elseif ( (abs(tr - 1).lt.tiny) .and. (abs(det - 1).lt.tiny) ) then; schoenflies_code = 4  ! 'c_4'
        elseif ( (abs(tr - 2).lt.tiny) .and. (abs(det - 1).lt.tiny) ) then; schoenflies_code = 5  ! 'c_6'
        elseif ( (abs(tr + 3).lt.tiny) .and. (abs(det + 1).lt.tiny) ) then; schoenflies_code = 6  ! 'i'
        elseif ( (abs(tr - 1).lt.tiny) .and. (abs(det + 1).lt.tiny) ) then; schoenflies_code = 7  ! 's_2'
        elseif ( (abs(tr - 0).lt.tiny) .and. (abs(det + 1).lt.tiny) ) then; schoenflies_code = 8  ! 's_6'
        elseif ( (abs(tr + 1).lt.tiny) .and. (abs(det + 1).lt.tiny) ) then; schoenflies_code = 9  ! 's_4'
        elseif ( (abs(tr + 2).lt.tiny) .and. (abs(det + 1).lt.tiny) ) then; schoenflies_code = 10 ! 's_3'
        endif
        !
        if (schoenflies_code .eq. 0) then
            call am_print('ERROR','Unable to identify point symmetry.',' >>> ')
            call am_print('R (frac)',R)
            call am_print('tr (frac)',tr)
            call am_print('det (frac)',det)
            stop
        endif
        !
    end function   point_symmetry_schoenflies

    !
    ! schoenflies decode
    !

    function       schoenflies_decode(schoenflies_code) result(schoenflies)
        !
        implicit none
        !
        integer, intent(in) :: schoenflies_code
        character(len=string_length_schoenflies) :: schoenflies
        !
        select case (schoenflies_code)
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
            call am_print('ERROR','Schoenflies code unknown.',' >>> ')
            stop
        end select
    end function   schoenflies_decode

    function       point_group_schoenflies_decode(pg_schoenflies_code) result(pg_schoenflies)
        !
        implicit none
        !
        integer, intent(in) :: pg_schoenflies_code
        character(string_length_schoenflies) :: pg_schoenflies
        !
        select case (pg_schoenflies_code)
                case(1); pg_schoenflies = 'c_1'
                case(2); pg_schoenflies = 's_2'
                case(3); pg_schoenflies = 'c_2'
                case(4); pg_schoenflies = 'c_1h'
                case(5); pg_schoenflies = 'c_2h'
                case(6); pg_schoenflies = 'd_2'
                case(7); pg_schoenflies = 'c_2v'
                case(8); pg_schoenflies = 'd_2h'
                case(9); pg_schoenflies = 'c_3'
                case(10); pg_schoenflies = 's_6'
                case(11); pg_schoenflies = 'd_3'
                case(12); pg_schoenflies = 'c_3v'
                case(13); pg_schoenflies = 'd_3d'
                case(14); pg_schoenflies = 'c_4'
                case(15); pg_schoenflies = 's_4'
                case(16); pg_schoenflies = 'c_4h'
                case(17); pg_schoenflies = 'd_4'
                case(18); pg_schoenflies = 'c_4v'
                case(19); pg_schoenflies = 'd_2d'
                case(20); pg_schoenflies = 'd_4h'
                case(21); pg_schoenflies = 'c_6'
                case(22); pg_schoenflies = 'c_3h'
                case(23); pg_schoenflies = 'c_6h'
                case(24); pg_schoenflies = 'd_6'
                case(25); pg_schoenflies = 'c_6v'
                case(26); pg_schoenflies = 'd_3h'
                case(27); pg_schoenflies = 'd_6h'
                case(28); pg_schoenflies = 't'
                case(29); pg_schoenflies = 't_h'
                case(30); pg_schoenflies = 'o'
                case(31); pg_schoenflies = 't_d'
                case(32); pg_schoenflies = 'o_h'
            case default
                call am_print('ERROR','Schoenflies code for point-group unknown.',' >>> ')
                call am_print('pg_schoenflies_code',pg_schoenflies_code)
                stop
            end select
    end function   point_group_schoenflies_decode

    function       point_group_schoenflies(schoenflies_code) result(pg_schoenflies_code)
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
        !>   pg_schoenflies_code --> pg_schoenflies
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
        integer, intent(in) :: schoenflies_code(:)
        integer :: nsyms
        integer :: ni, nc2, nc3, nc4, nc6, ns2, ns6, ns4, ns3, ne
        integer :: pg_schoenflies_code
        !
        pg_schoenflies_code = 0
        !
        ! trivial cases first
        !
        nsyms = size(schoenflies_code)
        if     (nsyms .eq. 1 ) then; pg_schoenflies_code=1;  return
        elseif (nsyms .eq. 48) then; pg_schoenflies_code=32; return
        elseif (nsyms .eq. 16) then; pg_schoenflies_code=20; return
        elseif (nsyms .eq. 3 ) then; pg_schoenflies_code=9;  return
        endif
        !
        ni=0; nc2=0; nc3=0; nc4=0; nc6=0; ns2=0; ns6=0; ns4=0; ns3=0
        ne  = count(schoenflies_code.eq.1)  ! 'e' identity
        nc2 = count(schoenflies_code.eq.2)  ! 'c_2'
        nc3 = count(schoenflies_code.eq.3)  ! 'c_3'
        nc4 = count(schoenflies_code.eq.4)  ! 'c_4'
        nc6 = count(schoenflies_code.eq.5)  ! 'c_6'
        ni  = count(schoenflies_code.eq.6)  ! 'i' inversion
        ns2 = count(schoenflies_code.eq.7)  ! 's_2'
        ns6 = count(schoenflies_code.eq.8)  ! 's_6'
        ns4 = count(schoenflies_code.eq.9)  ! 's_4'
        ns3 = count(schoenflies_code.eq.10) ! 's_3'
        !
        if (ni.gt.1) then
            call am_print('ERROR','More than one inversion operator found.',' >>> ')
            stop
        endif
        !
        if (ne.gt.1) then
            call am_print('ERROR','More than one identity operator found.',' >>> ')
        elseif (ne.eq.0) then
            call am_print('ERROR','No identity operator found.',' >>> ')
            stop
        endif
        !
        if (nsyms.eq.2)   then
            if (ni.eq.1)  then; pg_schoenflies_code=2;  return; endif
            if (nc2.eq.1) then; pg_schoenflies_code=3;  return; endif
            if (ns2.eq.1) then; pg_schoenflies_code=4;  return; endif
        endif
        if (nsyms.eq.4)   then
            if (ni.eq.1)  then; pg_schoenflies_code=5;  return; endif
            if (nc2.eq.3) then; pg_schoenflies_code=6;  return; endif
            if (ns2.eq.2) then; pg_schoenflies_code=7;  return; endif
            if (nc4.eq.1) then; pg_schoenflies_code=14; return; endif
            if (ns4.eq.2) then; pg_schoenflies_code=15; return; endif
        endif
        if (nsyms.eq.6)   then
            if (ni.eq.1)  then; pg_schoenflies_code=10; return; endif
            if (nc2.eq.3) then; pg_schoenflies_code=11; return; endif
            if (ns2.eq.3) then; pg_schoenflies_code=12; return; endif
            if (nc2.eq.1) then; pg_schoenflies_code=21; return; endif
            if (ns2.eq.1) then; pg_schoenflies_code=22; return; endif
        endif
        if (nsyms.eq.8)  then
            if (ns2.eq.3) then; pg_schoenflies_code=8;  return; endif
            if (ns2.eq.1) then; pg_schoenflies_code=16; return; endif
            if (ns2.eq.0) then; pg_schoenflies_code=17; return; endif
            if (ns2.eq.4) then; pg_schoenflies_code=18; return; endif
            if (ns2.eq.2) then; pg_schoenflies_code=19; return; endif
        endif
        if (nsyms.eq.12)  then
            if (ns2.eq.3) then; pg_schoenflies_code=13; return; endif
            if (ns2.eq.1) then; pg_schoenflies_code=23; return; endif
            if (nc2.eq.7) then; pg_schoenflies_code=24; return; endif
            if (ns2.eq.6) then; pg_schoenflies_code=25; return; endif
            if (ns2.eq.4) then; pg_schoenflies_code=26; return; endif
            if (nc3.eq.8) then; pg_schoenflies_code=28; return; endif
        endif
        if (nsyms.eq.24)  then
            if (nc6.eq.2) then; pg_schoenflies_code=27; return; endif
            if (ni.eq.1)  then; pg_schoenflies_code=29; return; endif
            if (nc4.eq.6) then; pg_schoenflies_code=30; return; endif
            if (ns4.eq.6) then; pg_schoenflies_code=31; return; endif
        endif
    end function   point_group_schoenflies

    !
    ! functions which operate on ss
    !

    subroutine     point_group(pg,ss,opts)
        !
        ! requires only ss%R
        !
        use am_rank_and_sort
        !
        implicit none
        !
        class(am_class_symmetry), intent(inout) :: pg
        type(am_class_symmetry) , intent(in) :: ss
        type(am_class_options), intent(in) :: opts
        integer :: pg_schoenflies_code
        character(string_length_schoenflies) :: pg_schoenflies
        integer , allocatable :: sorted_indices(:)
        real(dp), allocatable :: sort_parameter(:)
        integer :: i
        !
        if (opts%verbosity.ge.1) call am_print_title('Determining point group symmetries')
        !
        pg%R = unique(ss%R)
        !
        pg%nsyms = size(pg%R,3)
        !
        ! <sort>
        !
        allocate(sorted_indices(pg%nsyms))
        allocate(sort_parameter(pg%nsyms))
        !
        do i = 1, pg%nsyms; sort_parameter(i) = point_symmetry_determinant(pg%R(:,:,i)); enddo
        call rank(sort_parameter,sorted_indices)
        sorted_indices = sorted_indices(pg%nsyms:1:-1)
        pg%R=pg%R(:,:,sorted_indices)
        !
        do i = 1, pg%nsyms; sort_parameter(i) = point_symmetry_trace(pg%R(:,:,i)); enddo
        call rank(sort_parameter,sorted_indices)
        sorted_indices = sorted_indices(pg%nsyms:1:-1)
        pg%R=pg%R(:,:,sorted_indices)
        !
        ! </sort>
        !
        allocate(pg%schoenflies_code(pg%nsyms))
        allocate(pg%schoenflies(pg%nsyms))
        do i = 1, pg%nsyms
            pg%schoenflies_code(i) = point_symmetry_schoenflies(pg%R(:,:,i))
            pg%schoenflies(i)      = schoenflies_decode(pg%schoenflies_code(i))
        enddo
        !
        pg_schoenflies_code = point_group_schoenflies(pg%schoenflies_code)
        pg_schoenflies      = point_group_schoenflies_decode(pg_schoenflies_code)
        !
        if (opts%verbosity.ge.1) call am_print('point group',pg_schoenflies,' ... ')
        !
        if (opts%verbosity.ge.1) call am_print('number of point symmetries',pg%nsyms,' ... ')
        !
        if (opts%verbosity.ge.1) then
            call am_print('number of e   point symmetries',count(pg%schoenflies_code.eq.1),' ... ')
            call am_print('number of c_2 point symmetries',count(pg%schoenflies_code.eq.2),' ... ')
            call am_print('number of c_3 point symmetries',count(pg%schoenflies_code.eq.3),' ... ')
            call am_print('number of c_4 point symmetries',count(pg%schoenflies_code.eq.4),' ... ')
            call am_print('number of c_6 point symmetries',count(pg%schoenflies_code.eq.5),' ... ')
            call am_print('number of i   point symmetries',count(pg%schoenflies_code.eq.6),' ... ')
            call am_print('number of s_2 point symmetries',count(pg%schoenflies_code.eq.7),' ... ')
            call am_print('number of s_6 point symmetries',count(pg%schoenflies_code.eq.8),' ... ')
            call am_print('number of s_4 point symmetries',count(pg%schoenflies_code.eq.9),' ... ')
            call am_print('number of s_3 point symmetries',count(pg%schoenflies_code.eq.10),' ... ')
        endif
        !
        if (opts%verbosity.ge.1) call pg%stdout
        !
        call pg%stdout(trim('outfile.pointgroup'))
        !
    end subroutine point_group

    subroutine     space_group(ss,conv,opts)
        !
        use am_rank_and_sort
        !
        implicit none
        !
        class(am_class_symmetry), intent(inout) :: ss
        type(am_class_unit_cell), intent(in) :: conv
        type(am_class_options), intent(in) :: opts
        real(dp), allocatable :: seitz(:,:,:)
        integer , allocatable :: sorted_indices(:)
        real(dp), allocatable :: sort_parameter(:)
        integer :: i
        !
        if (opts%verbosity.ge.1) call am_print_title('Determining space group symmetries')
        !
        seitz = space_symmetries_from_basis(uc=conv,opts=opts)
        !
        ss%nsyms=size(seitz,3)
        !
        if (opts%verbosity.ge.1) call am_print('number of space group symmetries found',ss%nsyms,' ... ') 
        !
        allocate(ss%R(3,3,ss%nsyms))
        allocate(ss%T(3,ss%nsyms))
        ss%R(1:3,1:3,1:ss%nsyms) = seitz(1:3,1:3,1:ss%nsyms)
        ss%T(1:3,1:ss%nsyms)     = seitz(1:3,4,1:ss%nsyms)
        !
        ! <sort>
        !
        allocate(sorted_indices(ss%nsyms))
        allocate(sort_parameter(ss%nsyms))
        !
        do i = 1, ss%nsyms; sort_parameter(i) = point_symmetry_determinant(ss%R(:,:,i)); enddo
        call rank(sort_parameter,sorted_indices)
        sorted_indices = sorted_indices(ss%nsyms:1:-1)
        ss%R=ss%R(:,:,sorted_indices)
        ss%T=ss%T(:,sorted_indices)
        !
        do i = 1, ss%nsyms; sort_parameter(i) = point_symmetry_trace(ss%R(:,:,i)); enddo
        call rank(sort_parameter,sorted_indices)
        sorted_indices = sorted_indices(ss%nsyms:1:-1)
        ss%R=ss%R(:,:,sorted_indices)
        ss%T=ss%T(:,sorted_indices)
        !
        sort_parameter = ss%T(1,:)
        call rank(sort_parameter,sorted_indices)
        ss%R=ss%R(:,:,sorted_indices)
        ss%T=ss%T(:,sorted_indices)
        !
        sort_parameter = ss%T(2,:)
        call rank(sort_parameter,sorted_indices)
        ss%R=ss%R(:,:,sorted_indices)
        ss%T=ss%T(:,sorted_indices)
        !
        sort_parameter = ss%T(3,:)
        call rank(sort_parameter,sorted_indices)
        ss%R=ss%R(:,:,sorted_indices)
        ss%T=ss%T(:,sorted_indices)
        !
        ! </sort>
        !
        if (opts%verbosity.ge.1) call ss%stdout
        !
        call ss%stdout(trim('outfile.spacegroup'))
        !
    end subroutine space_group

    subroutine     stdout(ss,iopt_filename)
        !
        implicit none
        !
        class(am_class_symmetry), intent(in) ::  ss
        character(*), intent(in), optional :: iopt_filename
        integer :: fid
        integer :: i,j,k
        !
        if ( present(iopt_filename) ) then
            ! write to file
            fid = 1
            open(unit=fid,file=trim(iopt_filename),status="replace",action='write')
        else
            ! stdout
            fid = 6
        endif
        !
        ! write standard out to file
        !
        write(fid,'(5x,a3,9a5)',advance='no') '#', 'R11', 'R21', 'R31', 'R12', 'R22', 'R32', 'R13', 'R23', 'R33'
        if (allocated(ss%T))  write(fid,'(a4,3a5)',advance='no') '  | ', 'T1', 'T2', 'T3'
        write(fid,*)
        do i = 1, ss%nsyms
            write(fid,'(5x,i3,9i5)',advance='no') i, (( nint(ss%R(j,k,i)), j = 1,3) , k = 1,3)
            if (allocated(ss%T)) then
                write(fid,'(a4)',advance='no') '  | '
                do j = 1,3
                    write(fid,'(a5)',advance='no') trim(print_pretty(ss%T(j,i)))
                enddo
            endif
            write(fid,*)
        enddo
        !
        ! close file
        !
        if ( present(iopt_filename) ) close(fid)
        !
    end subroutine stdout

    subroutine     determine_symmetry(uc,ss,pg,opts)
        !
        implicit none
        !
        type(am_class_symmetry), intent(out) :: ss
        type(am_class_symmetry), intent(out) :: pg
        type(am_class_unit_cell), intent(in) :: uc
        type(am_class_unit_cell) :: conv, prim
        type(am_class_options), intent(in) :: opts
        !
        if (opts%verbosity.ge.1) call am_print_title('Analyzing symmetry')
        !
        call prim%reduce_to_primitive(uc=uc,opts=opts)
        !
        call conv%primitive_to_conventional(prim=prim,opts=opts)
        !
        call ss%space_group(conv=conv,opts=opts)
        !
        call pg%point_group(ss=ss,opts=opts)
        !
    end subroutine determine_symmetry

end module am_symmetry


