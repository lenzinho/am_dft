module am_symmetry
    !
    use am_constants
    use am_helpers
    use am_unit_cell
    use am_options
    use am_mkl
    !
    implicit none
    !
    private
    !
    public :: am_class_symmetry
    public :: determine_symmetry
    public :: get_kpoint_compatible_symmetries
    public :: ps_frac2cart
    public :: irreducible_cell

    public :: member
    public :: nelements
    public :: decode_pointgroup
    
    type am_class_symmetry 
        integer :: nsyms
        real(dp), allocatable :: R(:,:,:)         !> symmetry elements (operate on fractional atomic basis)
        real(dp), allocatable :: T(:,:)           !> symmetry elements (operate on fractional atomic basis)
        integer , allocatable :: cc_identifier(:) !> indices assiging each element to a conjugacy class
        integer , allocatable :: ps_identifier(:) !> integers which identiy point symmetries (see decode_pointsymmetry)
        integer               :: pg_identifier    !> integer which identifies the point group
    contains
        !
        procedure :: copy
        procedure :: create
        procedure :: space_group
        procedure :: point_group
        procedure :: rotational_group
        procedure :: stabilizer_group
        !
        procedure :: name_symmetries
        procedure :: stdout
        procedure :: sort
        !
        procedure :: determine_character_table
        procedure :: symmetry_action
    end type am_class_symmetry

contains

    !
    ! functions which operate on R(3,3) in fractional coordinates
    !

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

    function       ps_schoenflies(R) result(ps_identifier)
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
        integer :: ps_identifier
        !
        ! get trace and determinant (fractional)
        tr  = trace(R)
        d = det(R)
        !
        ! The Mathematical Theory of Symmetry in Solids: Representation Theory for
        ! Point Groups and Space Groups. 1 edition. Oxford?: New York: Oxford University
        ! Press, 2010. page 138, table 3.8.
        !
        ps_identifier = 0
        if     ( (abs(tr - 3).lt.tiny) .and. (abs(d - 1).lt.tiny) ) then; ps_identifier = 1  ! 'e'
        elseif ( (abs(tr + 1).lt.tiny) .and. (abs(d - 1).lt.tiny) ) then; ps_identifier = 2  ! 'c_2'
        elseif ( (abs(tr - 0).lt.tiny) .and. (abs(d - 1).lt.tiny) ) then; ps_identifier = 3  ! 'c_3'
        elseif ( (abs(tr - 1).lt.tiny) .and. (abs(d - 1).lt.tiny) ) then; ps_identifier = 4  ! 'c_4'
        elseif ( (abs(tr - 2).lt.tiny) .and. (abs(d - 1).lt.tiny) ) then; ps_identifier = 5  ! 'c_6'
        elseif ( (abs(tr + 3).lt.tiny) .and. (abs(d + 1).lt.tiny) ) then; ps_identifier = 6  ! 'i'
        elseif ( (abs(tr - 1).lt.tiny) .and. (abs(d + 1).lt.tiny) ) then; ps_identifier = 7  ! 's_2'
        elseif ( (abs(tr - 0).lt.tiny) .and. (abs(d + 1).lt.tiny) ) then; ps_identifier = 8  ! 's_6'
        elseif ( (abs(tr + 1).lt.tiny) .and. (abs(d + 1).lt.tiny) ) then; ps_identifier = 9  ! 's_4'
        elseif ( (abs(tr + 2).lt.tiny) .and. (abs(d + 1).lt.tiny) ) then; ps_identifier = 10 ! 's_3'
        endif
        !
        if (ps_identifier .eq. 0) then
            call am_print('ERROR','Unable to identify point symmetry.',flags='E')
            call am_print('R (frac)',R)
            call am_print('tr (frac)',tr)
            call am_print('det (frac)',d)
            stop
        endif
        !
    end function   ps_schoenflies

    !
    ! schoenflies decode
    !

    function       decode_pointsymmetry(ps_identifier) result(schoenflies)
        !
        implicit none
        !
        integer, intent(in) :: ps_identifier
        character(len=string_length_schoenflies) :: schoenflies
        !
        select case (ps_identifier)
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

    function       decode_pointgroup(pg_identifier) result(point_group_name)
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
        integer, intent(in) :: pg_identifier
        character(string_length_schoenflies) :: point_group_name
        !
        select case (pg_identifier)
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
                call am_print('pg_identifier',pg_identifier)
                stop
            end select
    end function   decode_pointgroup

    function       point_group_schoenflies(ps_identifier) result(pg_identifier)
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
        !>   pg_identifier --> point_group_name
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
        integer, intent(in) :: ps_identifier(:)
        integer :: nsyms
        integer :: ni, nc2, nc3, nc4, nc6, ns2, ns6, ns4, ns3, ne
        integer :: pg_identifier
        !
        pg_identifier = 0
        !
        ! trivial cases first
        !
        nsyms = size(ps_identifier)
        if     (nsyms .eq. 1 ) then; pg_identifier=1;  return
        elseif (nsyms .eq. 48) then; pg_identifier=32; return
        elseif (nsyms .eq. 16) then; pg_identifier=20; return
        elseif (nsyms .eq. 3 ) then; pg_identifier=9;  return
        endif
        !
        ni=0; nc2=0; nc3=0; nc4=0; nc6=0; ns2=0; ns6=0; ns4=0; ns3=0
        ne  = count(ps_identifier.eq.1)  ! 'e' identity
        nc2 = count(ps_identifier.eq.2)  ! 'c_2'
        nc3 = count(ps_identifier.eq.3)  ! 'c_3'
        nc4 = count(ps_identifier.eq.4)  ! 'c_4'
        nc6 = count(ps_identifier.eq.5)  ! 'c_6'
        ni  = count(ps_identifier.eq.6)  ! 'i' inversion
        ns2 = count(ps_identifier.eq.7)  ! 's_2'
        ns6 = count(ps_identifier.eq.8)  ! 's_6'
        ns4 = count(ps_identifier.eq.9)  ! 's_4'
        ns3 = count(ps_identifier.eq.10) ! 's_3'
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
            if (ni.eq.1)  then; pg_identifier=2;  return; endif
            if (nc2.eq.1) then; pg_identifier=3;  return; endif
            if (ns2.eq.1) then; pg_identifier=4;  return; endif
        endif
        if (nsyms.eq.4)   then
            if (ni.eq.1)  then; pg_identifier=5;  return; endif
            if (nc2.eq.3) then; pg_identifier=6;  return; endif
            if (ns2.eq.2) then; pg_identifier=7;  return; endif
            if (nc4.eq.1) then; pg_identifier=14; return; endif
            if (ns4.eq.2) then; pg_identifier=15; return; endif
        endif
        if (nsyms.eq.6)   then
            if (ni.eq.1)  then; pg_identifier=10; return; endif
            if (nc2.eq.3) then; pg_identifier=11; return; endif
            if (ns2.eq.3) then; pg_identifier=12; return; endif
            if (nc2.eq.1) then; pg_identifier=21; return; endif
            if (ns2.eq.1) then; pg_identifier=22; return; endif
        endif
        if (nsyms.eq.8)  then
            if (ns2.eq.3) then; pg_identifier=8;  return; endif
            if (ns2.eq.1) then; pg_identifier=16; return; endif
            if (ns2.eq.0) then; pg_identifier=17; return; endif
            if (ns2.eq.4) then; pg_identifier=18; return; endif
            if (ns2.eq.2) then; pg_identifier=19; return; endif
        endif
        if (nsyms.eq.12)  then
            if (ns2.eq.3) then; pg_identifier=13; return; endif
            if (ns2.eq.1) then; pg_identifier=23; return; endif
            if (nc2.eq.7) then; pg_identifier=24; return; endif
            if (ns2.eq.6) then; pg_identifier=25; return; endif
            if (ns2.eq.4) then; pg_identifier=26; return; endif
            if (nc3.eq.8) then; pg_identifier=28; return; endif
        endif
        if (nsyms.eq.24)  then
            if (nc6.eq.2) then; pg_identifier=27; return; endif
            if (ni.eq.1)  then; pg_identifier=29; return; endif
            if (nc4.eq.6) then; pg_identifier=30; return; endif
            if (ns4.eq.6) then; pg_identifier=31; return; endif
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

    !
    ! very high level routines which operate on sg
    !

    subroutine     determine_symmetry(uc,sg,pg,opts)
        !
        implicit none
        !
        class(am_class_unit_cell), intent(in) :: uc
        type(am_class_symmetry), intent(out) :: sg
        type(am_class_symmetry), intent(out) :: pg
        type(am_class_options)  , intent(in) :: opts
        !
        if (opts%verbosity.ge.1) call am_print_title('Analyzing symmetry')
        !
        call sg%space_group(uc=uc,opts=opts)
        !
        call pg%point_group(sg=sg,uc=uc,opts=opts)
        !
    end subroutine determine_symmetry

    subroutine     space_group(sg,uc,opts)
        !
        use am_rank_and_sort
        !
        implicit none
        !
        class(am_class_symmetry), intent(inout) :: sg
        type(am_class_unit_cell), intent(in) :: uc ! space groups are tabulated in the litearture for conventional cells, but primitive or arbitrary cell works just as wlel.
        type(am_class_options)  , intent(in) :: opts
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
        call put_identity_first(rep=seitz)
        !
        ! create space group instance (puts identity first)
        call sg%create(R=seitz(1:3,1:3,:),T=seitz(1:3,4,:))
        !
        ! determine conjugacy classes (needs identity first)
        sg%cc_identifier = get_conjugacy_classes(rep=seitz,flags='seitz')
        !
        ! sort space group symmetries based on parameters
        call sg%sort(sort_parameter=sg%T(1,:),iopt_direction='ascend')
        call sg%sort(sort_parameter=sg%T(2,:),iopt_direction='ascend')
        call sg%sort(sort_parameter=sg%T(3,:),iopt_direction='ascend')
        call sg%sort(sort_parameter=real(sg%cc_identifier,dp),iopt_direction='ascend')
        !
        ! name the (im-)proper part of space symmetry
        call sg%name_symmetries(opts=notalk)
        !
        ! write action table
        call sg%symmetry_action(uc=uc,flags='',iopt_fname='outfile.action_space_group',opts=opts)
        !
        ! get character table
        call sg%determine_character_table(opts=opts)
        !
        ! write to stdout and to file
        ! if (opts%verbosity.ge.1) call sg%stdout(iopt_uc=uc)
        call sg%stdout(iopt_uc=uc,iopt_filename=trim('outfile.spacegroup'))
        !
    end subroutine space_group

    subroutine     point_group(pg,uc,sg,opts)
        !
        ! requires only sg%R
        !
        use am_rank_and_sort
        !
        implicit none
        !
        class(am_class_symmetry), intent(inout) :: pg
        type(am_class_symmetry) , intent(in) :: sg
        type(am_class_unit_cell), intent(in) :: uc
        type(am_class_options), intent(in) :: opts
        real(dp), allocatable :: sort_parameter(:)
        type(am_class_options) :: notalk
        integer :: i
        !
        notalk=opts
        notalk%verbosity = 0
        !
        if (opts%verbosity.ge.1) call am_print_title('Determining point group symmetries')
        !
        ! create point group instance (puts identity first, inversion second)
        call pg%create(R=unique(sg%R))
        if (opts%verbosity.ge.1) call am_print('number of point symmetries',pg%nsyms,' ... ')
        !
        ! determine conjugacy classes (needs identity first, inversion second)
        pg%cc_identifier = get_conjugacy_classes(rep=pg%R)
        !
        ! sort point group symmetries based on parameters
        allocate(sort_parameter(pg%nsyms))
        do i = 1, pg%nsyms; sort_parameter(i) = trace(pg%R(:,:,i)); enddo
        call pg%sort(sort_parameter=sort_parameter,iopt_direction='ascend')
        do i = 1, pg%nsyms; sort_parameter(i) = det(pg%R(:,:,i)); enddo
        call pg%sort(sort_parameter=sort_parameter,iopt_direction='ascend')
        call pg%sort(sort_parameter=real(pg%cc_identifier,dp),iopt_direction='ascend')
        !
        ! name point symmetries
        call pg%name_symmetries(opts=opts)
        !
        ! name point group
        pg%pg_identifier = point_group_schoenflies(pg%ps_identifier)
        if (opts%verbosity.ge.1) call am_print('point group',decode_pointgroup(pg%pg_identifier),' ... ')
        !
        ! determine character table
        call pg%determine_character_table(opts=opts)
        !
        ! write action table 
        ! currently not working for nonsymorphic space groups. (i.e. those which necessairly have a translational component)
        ! actually it may just not be possible to do this for such a group because permutations linking atoms simply do not exist if only the rotational part is considered. the translational part is essential. an example is TiSe2.
        ! call pg%action_table(uc=uc,iopt_fname='outfile.action_point_group',opts=opts)
        !
        ! write to stdout and to file
        ! if (opts%verbosity.ge.1) call pg%stdout(iopt_uc=uc)
        if (opts%verbosity.ge.1) call pg%stdout(iopt_uc=uc,iopt_filename=trim('outfile.pointgroup'))
        !
    end subroutine point_group

    subroutine     rotational_group(rg,uc,pg,opts)
        ! get point symmetries which are compatible with the atomic basis (essentially keeps rotational part of space symmetries which are able to map all atoms onto each other)
        implicit none
        ! subroutine i/o
        class(am_class_symmetry), intent(inout) :: rg
        type(am_class_symmetry) , intent(in) :: pg
        class(am_class_unit_cell), intent(in) :: uc
        type(am_class_options)  , intent(in) :: opts
        type(am_class_options) :: notalk
        real(dp), allocatable  :: wrkspace(:,:,:)
        integer :: i, m
        !
        notalk=opts
        notalk%verbosity = 0
        !
        allocate(wrkspace(3,3,pg%nsyms))
        m=0
        do i = 1, pg%nsyms
            if ( is_symmetry_valid(tau=uc%tau,iopt_atype=uc%atype, &
               & iopt_R=pg%R(1:3,1:3,i),iopt_sym_prec=opts%sym_prec)) then
                m = m + 1
                wrkspace(1:3,1:3,m)=pg%R(1:3,1:3,i)
            endif
        enddo
        !
        call rg%create(R=unique(wrkspace(1:3,1:3,1:m)))
        !
        ! determine conjugacy classes (needs identity first)
        rg%cc_identifier = get_conjugacy_classes(rep=rg%R)
        !
        ! name the (im-)proper part of space symmetry
        call rg%name_symmetries(opts=notalk)
        !
        ! name "point group" - NOTE: not really a point group because translational components of space symmetries were ignored in the genereation process
        ! this would be the name of the group, if the group were symmorphic (in which case it will match the real point group name)... 
        rg%pg_identifier = point_group_schoenflies(rg%ps_identifier)
        !
    end subroutine rotational_group

    subroutine     stabilizer_group(vg,pg,v,opts)
        !
        implicit none
        !
        class(am_class_symmetry), intent(inout) :: vg ! stabilizer group associated with vector v
        type(am_class_symmetry) , intent(in) :: pg   ! space group
        real(dp)                , intent(in) :: v(3) ! vector which is stabilized (should be same units as R, which is fractional)
        type(am_class_options)  , intent(in) :: opts
        integer, allocatable :: indicies(:)
        logical, allocatable :: mask(:)
        integer :: i
        !
        allocate(mask(pg%nsyms))
        mask = .false.
        !
        do i = 1, pg%nsyms
            ! mark symmetries which leaves vector strictly invariant (not even modulo a lattice vector)
            if (is_symmetry_valid(tau=reshape(v,[3,1]),iopt_R=pg%R(:,:,i),&
                iopt_exact=.true.,iopt_sym_prec=opts%sym_prec)) then
                mask(i) = .true.
            endif
        enddo
        !
        allocate(indicies(pg%nsyms))
        indicies=[1:pg%nsyms]
        ! create group
        call vg%create(R=pg%R(:,:,pack(indicies,mask)))
        ! determine conjugacy classes (needs identity first, inversion second)
        vg%cc_identifier = get_conjugacy_classes(rep=vg%R)
        ! name point symmetries
        call vg%name_symmetries(opts=opts)
        ! name point group
        vg%pg_identifier = point_group_schoenflies(vg%ps_identifier)
        !
    end subroutine stabilizer_group

    subroutine     irreducible_cell(ic,pc,sg,opts)
        !
        implicit none
        !
        class(am_class_unit_cell), intent(inout) :: ic ! irreducible cell
        type(am_class_unit_cell), intent(in) :: pc
        type(am_class_symmetry), intent(in) :: sg
        type(am_class_options), intent(in) :: opts
        type(am_class_options) :: notalk
        integer, allocatable :: PM(:,:)
        logical, allocatable :: mask(:)
        integer, allocatable :: ind(:)
        integer :: i,j,k
        !
        if (opts%verbosity.ge.1) call am_print_title('Reducing to irreducible cell')
        !
        notalk = opts
        notalk%verbosity=0
        !
        if (opts%verbosity.ge.1) call am_print('number of atoms in primitive cell',pc%natoms,' ... ')
        !
        call sg%symmetry_action(uc=pc,flags='',oopt_PM=PM,opts=notalk)
        !
        allocate(mask(pc%natoms))
        allocate(ind(pc%natoms))
        mask = .true.
        !
        k=0
        do i = 1, pc%natoms
        if (mask(i)) then
            k=k+1
            ind(k)=i
            ! PM(uc%natoms,sg%nsyms) shows how atoms are permuted by each space symmetry operation
            do j = 1,sg%nsyms
                mask(PM(i,j))=.false.
            enddo
        endif
        enddo
        !
        ic%bas = pc%bas
        ic%natoms = k
        ic%nspecies = pc%nspecies
        !
        allocate(ic%symb ,source=pc%symb)
        allocate(ic%tau  ,source=pc%tau(:,ind(1:k)))
        allocate(ic%atype,source=pc%atype(ind(1:k)))
        !
        if (opts%verbosity.ge.1) call am_print('number of atoms in irreducible cell',ic%natoms,' ... ')
        if (opts%verbosity.ge.1) then
            call am_print_two_matrices_side_by_side(name='irreducible atomic basis',&
                Atitle='fractional',A=transpose(ic%tau),&
                Btitle='cartesian' ,B=transpose(matmul(ic%bas,ic%tau)),&
                iopt_emph=' ... ',iopt_teaser=.true.)
        endif
        !
    end subroutine irreducible_cell

    !
    ! medium level routines which operate on sg
    !

    subroutine     determine_character_table(sg,opts)
        !
        use am_rank_and_sort
        !
        implicit none
        !
        class(am_class_symmetry), intent(inout) :: sg
        type(am_class_options)  , intent(in) :: opts
        real(dp), allocatable :: rep(:,:,:)
        integer , allocatable :: cayley_table(:,:)
        character(100) :: flags
        integer :: i, j, k
        !
        integer :: nsyms
        integer :: nclasses
        integer :: nirreps  ! just to make things clearer; always equal to nclasses
        integer, allocatable :: class_nelements(:) ! number of elements in each class
        integer, allocatable :: class_member(:,:) ! members(nclass,maxval(class_nelements))
        integer, allocatable :: H(:,:,:) ! class multiplication 
        integer, allocatable :: CCT(:,:) ! class constant table which becomes character table (can be complex!)
        integer, allocatable :: irrep_dim(:) ! rep dimensions
        integer, allocatable :: indices(:)
        integer  :: index_of_class_containing_inversion
        real(dp) :: wrk
        !
        if (opts%verbosity.ge.1) call am_print_title('Determining character table')
        !
        ! seitz rep
        if (allocated(sg%T)) then
            rep = rep_seitz(sg%R,sg%T)
            flags='seitz'
        else
            rep = sg%R
            flags=''
        endif
        !
        nsyms = size(rep,3)
        !
        ! construct cayley table for quick access to symmetry multiplication
        cayley_table = get_cayley_table(rep=rep,flags=flags)
        !
        ! identify which class symmetry belongs to
        ! identifier(nsyms)
        sg%cc_identifier = get_conjugacy_classes(rep=rep,flags=flags)
        !
        ! get number of classes
        nclasses = maxval(sg%cc_identifier)
        ! if (opts%verbosity.ge.1) call am_print('conjugacy classes',nclasses,' ... ')
        !
        ! get number of elements in each class
        ! class_nelements(nclasses)
        class_nelements = nelements(sg%cc_identifier)
        !
        ! record indicies of each class_member element for each class (use the first class_member of class as class representative)
        ! members(nclass,maxval(class_nelements)) 
        class_member = member(sg%cc_identifier)
        ! call am_print('class memebers',class_member,' ... ')
        !
        ! get class coefficients (thenumber of times class k appears in the pdocut o class j and k )
        ! Wooten p 35, Eq 2.10; p 40, Example 2.8; p 89, Example 4.6
        allocate(H(nclasses,nclasses,nclasses))
        do i = 1,nclasses
        do j = 1,nclasses
        do k = 1,nclasses
            H(j,k,i) = count( pack( cayley_table(class_member(i,1:class_nelements(i)),class_member(j,1:class_nelements(j))) , .true. ) .eq. class_member(k,1) )
        enddo
        enddo
        enddo
        !
        ! determine eigenvector which simultaneously diagonalizes i Hjk(:,:) matrices
        ! CCT(nirreps,nclasses)
        CCT = class_constant_table(real(H,dp))
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
                wrk = wrk + abs(CCT(i,j))**2/real(class_nelements(j),dp)
            enddo
            ! compute irrep dimension
            irrep_dim(i) = nint((nsyms/wrk)**0.5)
        enddo
        !
        ! convert class constant table into character table
        do i = 1, nirreps
        do j = 1, nclasses
            CCT(i,j) = nint(irrep_dim(i)/real(class_nelements(j),dp)*CCT(i,j))
        enddo
        enddo
        !
        ! sort irreps
        allocate(indices(nirreps))
        ! sort irreps based on dimension
        call rank(irrep_dim,indices)
        irrep_dim = irrep_dim(indices)
        do j = 1,nclasses
            CCT(:,j) = CCT(indices,j)
        enddo
        ! sort irreps based on inversion (requires inversion class to be second)
        index_of_class_containing_inversion = 2
        call rank(-sign(1,CCT(:,index_of_class_containing_inversion)),indices)
        irrep_dim = irrep_dim(indices)
        do j = 1,nclasses
            CCT(:,j) = CCT(indices,j)
        enddo
        !
        ! write class properties to stdout
        !
        if (opts%verbosity.ge.1) then
            !
            write(*,'(5x,a10)',advance='no') 'class'
            do i = 1, nclasses
                write(*,'(i5)',advance='no') i
            enddo
            write(*,*)
            !
            write(*,'(5x,a10)',advance='no') 'elements'
            do i = 1, nclasses
                write(*,'(i5)',advance='no') class_nelements(i)
            enddo
            write(*,*)
            !
            write(*,'(5x,a10)',advance='no') 'class rep'
            do i = 1, nclasses
                write(*,'(i5)',advance='no') class_member(i,1)
            enddo
            write(*,*)
            !
            write(*,'(5x,a10)',advance='no') ' '
            do i = 1, nclasses
                write(*,'(a5)',advance='no') trim(decode_pointsymmetry(sg%ps_identifier(class_member(i,1))))
            enddo
            write(*,*)
            !
            write(*,'(5x,a10)',advance='no') repeat('-',10)
            do i = 1, nclasses
                write(*,'(a5)',advance='no') repeat('-',5)
            enddo
            write(*,*)
            !
            do i = 1, nirreps
                write(*,'(5x,a6,i4)',advance='no') 'irrep', i
                do j = 1, nclasses
                    write(*,'(i5)',advance='no') CCT(i,j)
                enddo
                write(*,*)
            enddo
            !
        endif
        ! check that the number of elements ri in class Ci is a divisor of theorder of the group
        ! Symmetry and Condensed Matter Physics: A Computational Approach. 1 edition. Cambridge, UK?; New York: Cambridge University Press, 2008. p 35.
        do i = 1, nclasses
            if (modulo(nsyms,class_nelements(i)).ne.0) then
                call am_print('ERROR','Number of elements in class subgroup is not a divisor of the order of the group.')
                call am_print('class',i)
                call am_print('elements in class',class_nelements(i))
                call am_print('elements in group',nsyms)
                stop
            endif
        enddo
        !
        contains
        function       class_constant_table(A) result(D)
            !
            ! determines characters which simultaneously diagonalizes i square matrices A(:,:,i)
            ! only possible if A(:,:,i) commute with each other. Assumes all dimensions of A are the same.
            ! 
            !
            implicit none
            !
            real(dp), intent(in)  :: A(:,:,:)
            real(dp), allocatable :: V(:,:) ! these two should match in this particular situation (character table determination from class multiplication )
            integer , allocatable :: D(:,:) ! these two should match in this particular situation (character table determination from class multiplication )
            real(dp), allocatable :: M(:,:)
            integer , allocatable :: p(:)
            integer :: n, i
            !
            call am_print('WARNING','The procedure used to determine the character table does not allow for complex character',flags='W')
            !
            n = size(A,3)
            p = primes(n)
            !
            ! matrix pencil based on sqrt(primes)
            ! all this does is lift the degeneracy of eigenvalues if necessary
            allocate(M(n,n))
            M = 0
            do i = 1,n
                M = M + p(i)**0.5*A(:,:,i)
            enddo
            !
            ! get left and right eigenvectors: transpose(VL)*A*VR = D
            call am_dgeev(A=M,VR=V)
            !
            ! check that matrices have been diagonalized properly and save eigenvalues
            allocate(d(n,n))
            do i = 1, n
                !
                M = matmul(inv(V),matmul(A(:,:,i),V))
                ! tiny may be too small a criteria here... should probably use a value which scales with the size of the matrix.
                if (any( abs(diag(diag(M))-M).gt. (tiny*n**2) )) then
                    call am_print('ERROR','Unable to perform simultaneous matrix diagonalization. Check that whether they commute.',flags='E')
                    call am_print('M',M)
                    call am_print('diag(M)',diag(M))
                    call am_print('diag(diag(M))',diag(diag(M)))
                    call am_print('abs(diag(diag(M))-M)',abs(diag(diag(M))-M))
                    call am_print('sum of error',sum(abs(diag(diag(M))-M)))
                    stop
                endif
                !
                ! save diagonal elements as row vectors
                D(:,i) = nint(diag( M ))
                !
            enddo
            !
            ! normalize each eigencolumn to the element of the column (that is, the first element, i.e. row, corresponds to the class containing the identity element)
            do i = 1, n
                V(:,i)=V(:,i)/V(1,i)
            enddo
            V = nint(transpose(V))
            !
            ! at this point V and D should be identical. The elements are to within +/- sign.
            ! Return D which is more robust than V. 
            !
        end function   class_constant_table
    end subroutine determine_character_table

    subroutine     name_symmetries(sg,opts)
        !
        implicit none
        !
        class(am_class_symmetry) , intent(inout) :: sg
        type(am_class_options), intent(in) :: opts
        integer :: i
        !
        if (allocated(sg%ps_identifier)) deallocate(sg%ps_identifier)
        allocate(sg%ps_identifier(sg%nsyms))
        !
        do i = 1, sg%nsyms
            sg%ps_identifier(i) = ps_schoenflies(sg%R(:,:,i))
        enddo
        !
        if (opts%verbosity.ge.1) then
            call am_print('number of e   point symmetries',count(sg%ps_identifier.eq.1),' ... ')
            call am_print('number of c_2 point symmetries',count(sg%ps_identifier.eq.2),' ... ')
            call am_print('number of c_3 point symmetries',count(sg%ps_identifier.eq.3),' ... ')
            call am_print('number of c_4 point symmetries',count(sg%ps_identifier.eq.4),' ... ')
            call am_print('number of c_6 point symmetries',count(sg%ps_identifier.eq.5),' ... ')
            call am_print('number of i   point symmetries',count(sg%ps_identifier.eq.6),' ... ')
            call am_print('number of s_2 point symmetries',count(sg%ps_identifier.eq.7),' ... ')
            call am_print('number of s_6 point symmetries',count(sg%ps_identifier.eq.8),' ... ')
            call am_print('number of s_4 point symmetries',count(sg%ps_identifier.eq.9),' ... ')
            call am_print('number of s_3 point symmetries',count(sg%ps_identifier.eq.10),' ... ')
        endif
        !
    end subroutine name_symmetries

    subroutine     symmetry_action(sg,uc,flags,iopt_fname,oopt_PM,oopt_P,opts)
        !
        ! PM(uc%natoms,sg%nsyms) permutation map; shows how atoms are permuted by each space symmetry operation
        !
        implicit none
        !
        class(am_class_symmetry) , intent(in) :: sg
        type(am_class_unit_cell), intent(in) :: uc
        type(am_class_options)  , intent(in) :: opts
        character(*)            , intent(in) :: flags ! if flags='relax_pbc', do not return atoms to primitive cell. if no matching atom is found, return index 0 for that orbit.
        character(*), optional  , intent(in) :: iopt_fname
        integer, optional, allocatable, intent(out) :: oopt_PM(:,:)
        integer, optional, allocatable, intent(out) :: oopt_P(:,:,:)
        integer, allocatable :: P(:,:,:)
        integer, allocatable :: PM(:,:)
        integer :: fid
        character(10) :: buffer
        character(3) :: xyz
        integer :: i,j,k
        real(dp) :: R(3,3)
        !
        if (opts%verbosity.ge.1) call am_print_title('Determining symmetry action')
        !
        P = rep_permutation(sg=sg,uc=uc,flags=flags,opts=opts)
        if (present(oopt_P)) allocate(oopt_P,source=P)
        !
        allocate(PM(uc%natoms,sg%nsyms))
        !
        PM = 0
        do i = 1,sg%nsyms
            PM(:,i) = matmul(P(:,:,i),[1:uc%natoms])
        enddo
        if (present(oopt_PM)) allocate(oopt_PM,source=PM)
        !
        if (present(iopt_fname)) then 
            !
            fid = 1
            open(unit=fid,file=trim(iopt_fname),status="replace",action='write')
                !
                ! SYMMETRY NUMBER
                write(fid,'(5x)',advance='no')
                write(fid,'(a5)',advance='no') '#'
                do i = 1, sg%nsyms
                    write(fid,'(i10)',advance='no') i
                enddo
                write(fid,*)
                ! CLASS IDENTIFIED
                write(fid,'(5x)',advance='no')
                write(fid,'(a5)',advance='no') 'class'
                do i = 1, sg%nsyms
                    write(fid,'(i10)',advance='no') sg%cc_identifier(i)
                enddo
                write(fid,*)
                ! HEADER / SEPERATOR
                write(fid,'(5x)',advance='no')
                write(fid,'(5x)',advance='no')
                do i = 1, sg%nsyms
                    write(fid,'(a10)',advance='no') ' '//repeat('-',9)
                enddo
                write(fid,*)
                ! POINT-SYMMMETRY IDENTIFIED
                write(fid,'(5x)',advance='no')
                write(fid,'(a5)',advance='no') 'SF'
                do i = 1, sg%nsyms
                    write(fid,'(a10)',advance='no') trim(decode_pointsymmetry(sg%ps_identifier(i)))
                enddo
                write(fid,*)
                ! POINT SYMMETRY COMPONENTS (FRAC)
                do i = 1,3
                do j = 1,3
                    write(fid,'(5x)',advance='no')
                    write(fid,'(a5)',advance='no') adjustr('R'//trim(int2char(i))//trim(int2char(j)))
                    do k = 1, sg%nsyms
                        write(fid,'(f10.2)',advance='no') sg%R(i,j,k)
                    enddo
                    write(fid,*)
                enddo
                enddo
                ! TRANSLATIONAL COMPONENTS (FRAC)
                if (allocated(sg%T)) then
                do j = 1,3
                    write(fid,'(5x)',advance='no')
                    write(fid,'(a5)',advance='no') adjustr('T'//trim(int2char(j)))
                    do i = 1, sg%nsyms
                        write(fid,'(f10.2)',advance='no') sg%T(j,i)
                    enddo
                    write(fid,*)
                enddo
                endif
                ! ACTION TABLE (AXIS PERMUTATION)
                write(fid,'(5x)',advance='no')
                write(fid,'(a5)',advance='no') 'axes'
                do i = 1, sg%nsyms
                    write(fid,'(a10)',advance='no') ' '//repeat('-',9)
                enddo
                write(fid,*)
                xyz = 'xyz'
                do j = 1,3
                    write(fid,'(5x)',advance='no')
                    write(fid,'(a5)',advance='no') trim(xyz(j:j))
                    do i = 1, sg%nsyms
                        buffer=''
                        do k = 1,3
                            ! for fractional R; its elements will always be an integer...
                            ! inverse here because f(Rr) = R^-1 * f(r)
                            R=inv(sg%R(:,:,i))
                            if (abs(R(j,k)).gt.tiny) then
                                buffer = trim(buffer)//trim(int2char(nint(R(j,k)),'SP'))//xyz(j:j)
                            endif
                        enddo
                        write(fid,'(a10)',advance='no') trim(buffer)
                    enddo
                    write(fid,*)
                enddo
                ! ACTION TABLE (ATOM PERMUTATION)
                write(fid,'(5x)',advance='no')
                write(fid,'(a5)',advance='no') 'atoms'
                do i = 1, sg%nsyms
                    write(fid,'(a10)',advance='no') ' '//repeat('-',9)
                enddo
                write(fid,*)
                !
                do j = 1, uc%natoms
                    write(fid,'(5x)',advance='no')
                    write(fid,'(i5)',advance='no') j
                    do i = 1, sg%nsyms
                        write(fid,'(i10)',advance='no') PM(j,i)
                    enddo
                    write(fid,*)
                enddo
                !
            close(fid)
        endif
    end subroutine symmetry_action

    subroutine     create(sg,R,T)
        !
        implicit none
        !
        class(am_class_symmetry), intent(inout) :: sg
        real(dp), intent(in), optional :: R(:,:,:)
        real(dp), intent(in), optional :: T(:,:)
        real(dp) :: id(3,3)
        integer  :: i
        !
        if (present(T).and.present(R)) then
        if (size(R,3).ne.size(T,2)) then
            call am_print('ERROR','R and T dimensions do not match.')
            call am_print('shape(R)',shape(R))
            call am_print('shape(T)',shape(T))
            stop
        endif
        endif
        !
        if (present(R)) then
            sg%nsyms = size(R,3)
            if (allocated(sg%R)) deallocate(sg%R)
            allocate(sg%R(3,3,sg%nsyms))
            do i = 1, sg%nsyms
                sg%R(:,:,i) = R(:,:,i)
            enddo
        endif
        !
        if (present(T)) then
            sg%nsyms = size(T,2)
            if (allocated(sg%T)) deallocate(sg%T)
            allocate(sg%T(3,sg%nsyms))
            do i = 1, sg%nsyms
                sg%T(:,i) = T(:,i)
            enddo
        endif
        !
        ! put identity first
        !
        id = eye(3)
        do i = 1, sg%nsyms
            if (present(R)) then
                if (present(T)) then
                    ! R and T
                    if (all(abs(sg%R(:,:,i)-id).lt.tiny)) then
                    if (all(abs(sg%T(:,i)).lt.tiny)) then
                        sg%R(:,:,i) = sg%R(:,:,1)
                        sg%T(:,i)   = sg%T(:,1)
                        sg%R(:,:,1) = id
                        sg%T(:,1)   = real([0,0,0],dp)
                    endif
                    endif
                else
                    ! just R
                    if (all(abs(sg%R(:,:,i)-id).lt.tiny)) then
                        sg%R(:,:,i) = sg%R(:,:,1)
                        sg%R(:,:,1) = id
                    endif
                endif
            else
                ! just T
                if (all(abs(sg%T(:,i)).lt.tiny)) then
                    sg%T(:,i) = sg%T(:,1)
                    sg%T(:,1) = real([0,0,0],dp)
                endif
            endif
        enddo
    end subroutine create

    subroutine     copy(cp,sg)
        !
        implicit none
        !
        class(am_class_symmetry), intent(inout) :: cp
        type(am_class_symmetry) , intent(in) :: sg
        !
        cp%nsyms         = sg%nsyms
        cp%pg_identifier = sg%pg_identifier
        if (allocated(sg%cc_identifier)) allocate(cp%cc_identifier  , source=sg%cc_identifier)
        if (allocated(sg%R))             allocate(cp%R              , source=sg%R)
        if (allocated(sg%T))             allocate(cp%T              , source=sg%T)
        if (allocated(sg%ps_identifier)) allocate(cp%ps_identifier  , source=sg%ps_identifier)
        !
    end subroutine copy
   
    subroutine     sort(sg,sort_parameter,iopt_direction)
        !
        ! iopt_direction = 'ascend'/'descend', only first character is important
        !
        use am_rank_and_sort
        !
        implicit none
        !
        class(am_class_symmetry) , intent(inout) :: sg
        character(*), intent(in), optional :: iopt_direction
        integer , allocatable :: sorted_indices(:)
        real(dp) :: sort_parameter(:)
        character(1) :: direction
        !
        if (present(iopt_direction)) then
            direction = iopt_direction(1:1)
        else 
            direction = 'a'
        endif
        !
        allocate(sorted_indices(sg%nsyms))
        !
        call rank(sort_parameter,sorted_indices)
        !
        if (direction.eq.'d') sorted_indices = sorted_indices(sg%nsyms:1:-1)
        !
        if (allocated(sg%R)) sg%R = sg%R(:,:,sorted_indices)
        if (allocated(sg%T)) sg%T = sg%T(:,sorted_indices)
        if (allocated(sg%ps_identifier)) sg%ps_identifier = sg%ps_identifier(sorted_indices)
        if (allocated(sg%cc_identifier)) sg%cc_identifier = sg%cc_identifier(sorted_indices)
        !
    end subroutine sort

    subroutine     stdout(sg,iopt_uc,iopt_filename)
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
            ! stdout
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
        if (allocated(sg%T)) width=width + 3*7+4
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
        if (allocated(sg%T))  write(fid,fmt4,advance='no') '  | ', 'T1', 'T2', 'T3'
        write(fid,*)
        !
        do i = 1, sg%nsyms
            R = sg%R(:,:,i)
            write(fid,fmt1,advance='no') i, (( trim(print_pretty(R(j,k))), j = 1,3) , k = 1,3)
            if (allocated(sg%T)) then
                T = sg%T(:,i)
                write(fid,fmt5,advance='no') '  | '
                do j = 1,3
                    write(fid,fmt6,advance='no') trim(print_pretty(T(j)))
                enddo
            endif
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
        if (allocated(sg%T))  write(fid,fmt4,advance='no') '  | ', 'T1', 'T2', 'T3'
        write(fid,*)
        !
        do i = 1, sg%nsyms
            R = ps_frac2cart(R_frac=sg%R(:,:,i),bas=iopt_uc%bas)
            write(fid,fmt1,advance='no') i, (( trim(print_pretty(R(j,k))), j = 1,3) , k = 1,3)
            if (allocated(sg%T)) then
                T = matmul(iopt_uc%bas,sg%T(:,i))
                write(fid,fmt5,advance='no') '  | '
                do j = 1,3
                    write(fid,fmt6,advance='no') trim(print_pretty(T(j)))
                enddo
            endif
            write(fid,*)
        enddo
        !
        endif
        !
        ! close file
        !
        if ( present(iopt_filename) ) close(fid)
        !
    end subroutine stdout

    subroutine     convert(sg,uc,flags)
        !
        ! one way to remember which to multiply on the left and right is to make form the identity; to convert V and R
        ! in fractional coordinates to cartesian: V -> BV; R -> B*R*inv(B); because: V_rot_cart = B*R*inv(B) * B*V =
        ! B*R*I*V = B*R_frac*V_frac
        !
        implicit none
        !
        type(am_class_symmetry) , intent(inout) :: sg
        type(am_class_unit_cell), intent(inout) :: uc
        character(*), intent(in) :: flags
        real(dp) :: recbas(3,3), bas(3,3)
        integer :: i
        !
        if     (index(flags,'frac2cart').ne.0) then
            bas    = uc%bas
            recbas = inv(uc%bas)
        elseif (index(flags,'cart2frac').ne.0) then
            bas    = inv(uc%bas)
            recbas = uc%bas
        endif
        !
        ! atomic basis
        do i = 1,uc%natoms
            uc%tau(:,i) = matmul(bas,uc%tau(:,i))
        enddo
        !
        ! rotational component
        if (allocated(sg%R)) then
        do i = 1,sg%nsyms
            sg%R(:,:,i) = matmul(bas,matmul(sg%R(:,:,i),recbas))
        enddo
        endif
        !
        ! translational component
        if (allocated(sg%T)) then
        do i = 1,sg%nsyms
            sg%T(:,i)   = matmul(bas,sg%T(:,i))
        enddo
        endif
        !
    end subroutine convert

    !
    ! procedures which create representations
    !

    function       rep_seitz(R,T) result(seitz)
        !
        implicit none
        !
        real(dp), intent(in) :: R(:,:,:)
        real(dp), intent(in) :: T(:,:)
        real(dp), allocatable :: seitz(:,:,:)
        integer :: i, n
        !
        n = size(R,3)
        !
        allocate(seitz(4,4,n))
        seitz = 0
        seitz(4,4,:)=1
        do i = 1, n
            seitz(1:3,1:3,i)=R(:,:,i)
            seitz(1:3,4,i)  =T(:,i)
        enddo
        !
    end function   rep_seitz

    function       rep_permutation(sg,uc,flags,opts) result(rep)
        !
        ! find permutation representation; i.e. which atoms are connected by space symmetry oprations R, T.
        !
        implicit none
        !
        type(am_class_symmetry) , intent(in) :: sg
        type(am_class_unit_cell), intent(in) :: uc
        character(*)            , intent(in) :: flags
        type(am_class_options)  , intent(in) :: opts
        integer, allocatable :: rep(:,:,:)
        real(dp) :: tau_rot(3)
        integer :: i,j,k
        logical :: found
        !
        allocate(rep(uc%natoms,uc%natoms,sg%nsyms))
        rep = 0
        !
        do i = 1, sg%nsyms
            ! determine the permutations of atomic indicies which results from each space symmetry operation
            do j = 1,uc%natoms
                found = .false.
                ! apply rotational component
                tau_rot = matmul(sg%R(:,:,i),uc%tau(:,j))
                ! apply translational component
                if (allocated(sg%T)) tau_rot = tau_rot + sg%T(:,i)
                ! reduce rotated+translated point to unit cell
                if (index(flags,'relax_pbc').eq.0) then
                    tau_rot = modulo(tau_rot+opts%sym_prec,1.0_dp)-opts%sym_prec
                endif
                ! find matching atom
                do k = 1, uc%natoms
                    if (all(abs(tau_rot-uc%tau(:,k)).lt.opts%sym_prec)) then
                    rep(j,k,i) = 1
                    found = .true.
                    exit ! break loop
                    endif
                enddo
                ! if "relax_pbc" is present (i.e. periodic boundary conditions are relax), perform check
                ! to ensure atoms must permute onto each other
                if (index(flags,'relax_pbc').eq.0) then
                if (found.eq..false.) then
                    call am_print('ERROR','Unable to find matching atom.',flags='E')
                    call am_print('tau (all atoms)',transpose(uc%tau))
                    call am_print('tau',uc%tau(:,j))
                    call am_print('R',sg%R(:,:,i))
                    if (allocated(sg%T)) call am_print('T',sg%T(:,i))
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
        do i = 1,sg%nsyms
        do j = 1,uc%natoms
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
    end function   rep_permutation

    function       rep_regular(rep,flags) result(reg_rep)
        !>
        !> Generates regular representation for each operator
        !>
        !> For details, see page 73-74, especially Eqs. 4.18 Symmetry and Condensed Matter Physics: A Computational
        !> Approach. 1 edition. Cambridge, UK; New York: Cambridge  University Press, 2008.
        !>
        !> Requires identity to be the first element of group. regular representation is built from multiplication table
        !> by making sure that the identity appears along the diagonal. this sometimes causes the row to be permuted with
        !> respect to the columns. so, first find the permutation necessary to bring the identities to the center and
        !> then build regular representations from this new multiplication table.
        !>
        implicit none
        !
        real(dp), intent(in) :: rep(:,:,:)
        character(*), intent(in) :: flags ! can be nothing flags='' at call
        integer, allocatable :: reg_rep(:,:,:)
        integer, allocatable :: cayley_table(:,:)
        integer :: n
        integer :: i
        !
        ! obtain cayley table
        cayley_table = get_cayley_table(rep=rep,flags=flags//'reg')
        !
        ! construct regular representation from cayley table
        n = size(rep,3)
        allocate(reg_rep(n,n,n))
        reg_rep = 0
        do i = 1, n
            where (cayley_table.eq.i) reg_rep(:,:,i) = 1
        enddo
        !
    end function   rep_regular

    !
    ! procedures which operate on matrix representations
    !

    function       get_cayley_table(rep,flags) result(cayley_table)
        !
        ! flag options:
        !  reg    -  sorts rows to put identity along diagonals (useful for constructing regular representation)
        !  seitz  -  reduces translational part to primitive cell
        !
        implicit none
        !
        real(dp), intent(in) :: rep(:,:,:) ! list of 2D reps...
        character(*), intent(in), optional :: flags
        integer, allocatable :: cayley_table(:,:)
        integer, allocatable :: sortmat(:,:)
        integer, allocatable :: indices(:)
        real(dp) :: W(size(rep,1),size(rep,2)) ! workspace
        integer  :: n
        integer  :: i, j
        !
        ! before doing anything, confirm that first element is the identity
        if (any(abs(rep(:,:,1)-eye(size(rep,1))).gt.tiny)) then
            call am_print('ERROR','First element of rep is not the identity.')
            call am_print('first element of rep',rep(:,:,1))
            stop
        endif
        !
        n=size(rep,3)
        !
        allocate(cayley_table(n,n))
        !
        do i = 1, n
        do j = 1, n
            !
            W = matmul(rep(:,:,i),rep(:,:,j))
            !
            ! if 'seitz' flagged, reduce translational part to primitive cell 
            if (present(flags)) then
            if (index(flags,'seitz').ne.0) then
                W(1:3,4) = modulo(W(1:3,4)+tiny,1.0_dp)-tiny
            endif
            endif
            !
            cayley_table(i,j) = get_matching_element_index_in_list(list=rep,elem=W)
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
            where (cayley_table.eq.1) sortmat(:,:) = 1
            indices = matmul(sortmat(:,:),indices)
            cayley_table = cayley_table(indices,:)
        endif
        endif
        !
        ! quick consitency check. not comprehensive
        do i = 1, n
            if (sum(cayley_table(:,i)).ne.sum(cayley_table(:,1))) then
                call am_print('ERROR','Cayley table sums along rows/columns are not consistent.')
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
    end function   get_cayley_table

    subroutine     put_identity_first(rep)
        !
        implicit none
        !
        real(dp), intent(inout) :: rep(:,:,:)
        real(dp), allocatable :: id(:,:)
        integer :: i, nsyms
        !
        nsyms = size(rep,3)
        id = eye(size(rep,1))
        do i = 1, nsyms
            if (all((abs(rep(:,:,i)-id)).lt.tiny)) then
                rep(:,:,i) = rep(:,:,1)
                rep(:,:,1) = id
            endif
        enddo
        ! 
    end subroutine put_identity_first

    function       get_conjugacy_classes(rep,flags) result(cc_identifier)
        !
        ! for AX = XB, if elements A and B are conjugate pairs for some other element X in the group, then they are in the same class
        !
        use am_rank_and_sort
        !
        implicit none
        !
        real(dp), intent(in) :: rep(:,:,:)
        character(*), intent(in), optional :: flags
        integer, allocatable :: cayley_table(:,:) ! cayley table
        integer, allocatable :: sinv(:) ! sinv(nsyms) index of the inverse of symmetry element i
        integer, allocatable :: cc_identifier(:)
        integer, allocatable :: celem(:)
        ! integer, allocatable :: class_member(:,:) ! used to sort based on determinant
        ! integer, allocatable :: indices(:)        ! used to sort based on determinant
        ! integer, allocatable :: relabel_indices(:)! used to sort based on determinant
        ! real(dp),allocatable :: ccdet(:)          ! used to sort based on determinant
        integer :: i, j, k
        integer :: nsyms
        integer :: nclasses
        !
        nsyms = size(rep,3)
        !
        ! get cayley table
        if (present(flags)) then
            cayley_table = get_cayley_table(rep=rep,flags=flags)
        else
            cayley_table = get_cayley_table(rep=rep)
        endif
        !
        ! allocate space for conjugacy class
        allocate(cc_identifier(nsyms))
        cc_identifier = 0
        !
        ! get element inverse
        sinv = get_inverse_indices(cayley_table)
        !
        ! allocate space for conjugate elements
        allocate(celem(nsyms))
        !
        ! determine conjugacy classes
        k = 0
        do i = 1, nsyms
        if (cc_identifier(i).eq.0) then
            k=k+1
            ! conjugate each element with all other group elements
            ! A = X(j) * B * X(j)^-1
            do j = 1, nsyms
                celem(j) = cayley_table(j,cayley_table(i,sinv(j)))
            enddo
            ! for each subgroup element created by conjugation find the corresponding index of the element in the group
            ! in order to save the class cc_identifier number
            do j = 1, nsyms
                cc_identifier( celem(j) ) = k
            enddo
            !
        endif
        enddo
        nclasses = k
        !
        ! relabel classes based on number of elements in each class
        call relabel_based_on_occurances(cc_identifier)
        !
        ! relabel classes based on det of class representative
        ! class_member = member(cc_identifier)
        ! allocate(ccdet(nclasses)) 
        ! do i = 1,nclasses
        !     ccdet(i) = det( rep(:,:,class_member(i,1)) )
        ! enddo
        ! !
        ! allocate(indices(nclasses))
        ! call rank(ccdet,indices)
        ! !
        ! allocate(relabel_indices(nclasses))
        ! do i = 1, nclasses
        !     relabel_indices(indices(i)) = i
        ! enddo
        ! do i = 1, nsyms
        !     cc_identifier(i) = relabel_indices( cc_identifier(i) )
        ! enddo
        !
        ! make sure the identity is in the first class
        ! cc_identifier(i=1) is the class of the first element (the identity); swap it's location with whaterver elements are in the first class
        where (cc_identifier.eq.1) cc_identifier = cc_identifier(1)
        cc_identifier(1)=1
        !
        ! check that classes are disjoint and complete
        if (any(cc_identifier.eq.0)) then
            call am_print('ERROR','Not every element in the group has been asigned a conjugacy class.',flags='E')
            stop
        endif
        !
        contains
        subroutine     relabel_based_on_occurances(list)
            !
            use am_rank_and_sort
            !
            implicit none
            !
            integer , intent(inout) :: list(:)
            integer , allocatable :: A_sorted(:)
            integer , allocatable :: occurances(:) 
            integer , allocatable :: reverse_sort(:)
            integer , allocatable :: sorted_indices(:)
            integer , allocatable :: list_relabled(:)
            integer :: nsyms, i, j
            !
            nsyms = size(list,1)
            !
            allocate(occurances(nsyms))
            do i = 1, nsyms
                occurances(i) = count(list(i).eq.list)
            enddo
            !
            ! this quick and dirty procedure lifts degeneracies
            !
            allocate(sorted_indices(nsyms))
            allocate(A_sorted(nsyms))
            call rank((1+maxval(list))*occurances+list,sorted_indices)
            A_sorted = list(sorted_indices)
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
            list=list_relabled(reverse_sort)
        end subroutine relabel_based_on_occurances
    end function   get_conjugacy_classes

    !
    ! procedures which operate on cayley table 
    ! (the idea should be to get the rep, transform it to the cayley table and perform as many operations on the cayley table as necessary)

    pure function  get_inverse_indices(cayley_table) result(ind)
        ! returns the ind of the inverse elements given the cayley table.
        ! requires identity as first element
        implicit none
        !
        integer, intent(in) :: cayley_table(:,:)
        integer, allocatable :: ind(:)
        integer, allocatable :: P(:,:)
        !
        ind = [1:size(cayley_table,1)]
        !
        allocate(P,source=cayley_table)
        where (P.ne.1) P = 0
        !
        ind = matmul(P,ind)
        !
    end function   get_inverse_indices

    function       get_coset(g_identifier,H_identifier,cayley_table,flags) result(ind)
        !
        ! In mathematics, if G is a group, and H is a subgroup of G, and g is an element of G, then:
        ! gH = { gh : h an element of H } is the left coset of H in G with respect to g, and
        ! Hg = { hg : h an element of H } is the right coset of H in G with respect to g.
        ! Wikipedia.
        !
        implicit none
        !
        integer, intent(in) :: g_identifier
        integer, intent(in) :: H_identifier(:)
        integer, intent(in) :: cayley_table(:,:)
        character(*), intent(in) :: flags ! left/right
        integer, allocatable :: ind(:)
        integer :: i, m
        !
        m = size(H_identifier)
        !
        allocate(ind(m))
        !
        do i = 1, m
            ! left/right coset, coonjucate 
            if     (index(flags,'right')) then; ind(i) = cayley_table(H_identifier(i),g_identifier)
            elseif (index(flags,'left' )) then; ind(i) = cayley_table(g_identifier,H_identifier(i))
            else
                call am_print('ERROR','Flag '//trim(flags)//' not valid.',flags='E')
                stop
            endif
        enddo
        !
    end function   get_coset

    ! 
    ! functions which operate on conjugate classes identifiers
    ! 

    function       member(identifier) result(class_member)
        !
        implicit none
        !
        integer, allocatable, intent(in) :: identifier(:)
        integer, allocatable :: class_nelements(:) ! number of elements in each class
        integer, allocatable :: class_member(:,:) ! members(nclass,maxval(class_nelements))
        integer :: i, j, k
        integer :: nsyms
        integer :: nclasses
        !
        nsyms = size(identifier,1)
        !
        nclasses = maxval(identifier)
        !
        class_nelements = nelements(identifier)
        !
        allocate(class_member(nclasses,maxval(class_nelements)))
        class_member = 0
        do i = 1, nclasses
            k=0
            do j = 1, nsyms
                if (identifier(j).eq.i) then
                    k=k+1
                    class_member(i,k) = j
                endif
            enddo
        enddo
    end function   member

    function       nelements(identifier) result(class_nelements)
        !
        implicit none
        !
        integer, intent(in)  :: identifier(:)
        integer, allocatable :: class_nelements(:)
        integer :: nclasses
        integer :: i
        !
        nclasses = maxval(identifier)
        !
        allocate(class_nelements(nclasses))
        !
        do i = 1, nclasses
            class_nelements(i) = count(i.eq.identifier)
        enddo
    end function   nelements

    !
    ! functions which operate on kpoints
    !

    function       get_kpoint_compatible_symmetries(kpt,R,sym_prec) result(R_syms)
        !
        implicit none
        !
        real(dp), intent(in) :: R(:,:,:)
        real(dp), intent(in) :: kpt(:,:)
        real(dp), intent(in) :: sym_prec
        real(dp), allocatable :: R_syms(:,:,:)
        real(dp) :: wrkspace(3,3,size(R,3))
        integer  :: i, j, nsyms
        !
        nsyms = size(R,3)
        !
        j=0
        do i = 1, nsyms
            if (is_symmetry_valid(tau=kpt,iopt_R=R,iopt_sym_prec=sym_prec)) then
                j=j+1
                wrkspace(:,:,j) = R(:,:,i)
            endif
        enddo
        allocate(R_syms(3,3,j))
        R_syms = wrkspace(:,:,1:j)
        !
    end function   get_kpoint_compatible_symmetries


        ! pure function  reduce_to_wigner_seitz(pnt,grid_points,bas,sym_prec) result(pnt_reduced)
        !     !> reduces pnt (in fractional) to the first Brillouin zone (Wigner-Seitz cell, defined in cartesian coordinates)
        !     !> cartesian pnt is returned! 
        !     implicit none
        !     !
        !     real(dp), intent(in) :: pnt(3) !> fractional
        !     real(dp), intent(in) :: bas(3,3) !> real space basis (column vectors)
        !     real(dp), intent(in) :: grid_points(3,27) !> voronoi points (cartesian)
        !     real(dp), intent(in) :: sym_prec
        !     real(dp) :: pnt_cart(3) !> pnt cartesian
        !     real(dp) :: G(3) !> reciprocal lattice vector
        !     real(dp) :: pnt_reduced(3)
        !     integer :: i ! loop variable
        !     logical :: is_not_done
        !     ! real(dp) :: small_offset(3)
        !     !
        !     ! take pnt in fractional coordinates (will become cartesian later)
        !     pnt_cart = pnt
        !     ! reduce to reciprocal unit cell
        !     pnt_cart = modulo(pnt_cart+sym_prec,1.0_dp)-sym_prec
        !     ! convert to cartesian
        !     pnt_cart = matmul(bas,pnt_cart)
        !     ! reduce to Wigner-Seitz cell (aka first Brillouin zone) by translating the k-point
        !     ! until the closest reciprocal lattice point is [0 0 0]
        !     is_not_done = .true.
        !     do while ( is_not_done )
        !         is_not_done = .false.
        !         do i = 1,27
        !             G = grid_points(:,i)
        !             ! bragg plane
        !             if ( 2*dot_product(pnt_cart,G) .gt. dot_product(G,G)+sym_prec ) then
        !                 pnt_cart = pnt_cart - G
        !                 is_not_done = .true.
        !             endif
        !         enddo
        !     end do
        !     ! convert back to fractional
        !     pnt_reduced = matmul(inv(bas),pnt_cart)
        !     !
        ! end function   reduce_to_wigner_seitz


end module am_symmetry


