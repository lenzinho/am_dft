#:include "fypp_macros.fpp"
module am_unit_cell

    use am_constants
    use am_stdout
    use am_options
    use am_mkl
    use dispmodule
    use am_symmetry

    implicit none

    public

    ! unit cell (generic)

    type, public, extends(am_class_cell) :: am_class_unit_cell
        integer , allocatable :: pc_id(:) ! identifies corresponding atom in primitive cell
        integer , allocatable :: ic_id(:) ! identifies corresponding atom in irreducible cell
        integer :: lattice_id
        contains
        procedure :: write_poscar
        procedure :: get_supercell
        procedure :: copy
        procedure :: filter
        procedure :: sort_atoms
        procedure :: initialize
        procedure :: save => save_uc
        procedure :: load => load_uc
        procedure, private :: write_uc
        procedure, private :: read_uc
        generic :: write(formatted) => write_uc
        generic :: read(formatted) => read_uc
    end type am_class_unit_cell

    ! primitive

    type, public, extends(am_class_unit_cell) :: am_class_prim_cell
        contains
        procedure :: get_primitive => get_primitive_cell
    end type am_class_prim_cell

    type, public, extends(am_class_unit_cell) :: am_class_conv_cell
        real(dp) :: centering(3,3) ! centering matrix
        contains
        ! procedure :: get_conventional => get_conventional_cell
    end type am_class_conv_cell

    ! irreducible

    type, public, extends(am_class_prim_cell) :: am_class_irre_cell
        ! atom is used by tight binding
        class(am_class_atom), allocatable :: atom(:)
        !
        contains
        procedure :: get_irreducible => get_irreducible_cell
        procedure :: initialize_orbitals
    end type am_class_irre_cell

    ! shell

    type, public, extends(am_class_unit_cell) :: am_shell_cell
        logical :: isonsite ! onsite
        real(dp):: center(3) ! center of shell [cart]
        integer :: i ! identifies irreducible atoms (center)
        integer :: j ! identifies irreducible atoms (shell)
        integer :: m ! identifies primitive atoms (center)
        integer :: n ! identifies primitive atoms (shell)
        integer, allocatable :: pg_id(:)   ! identifies the symmetry in the point group which takes atom tau_frac(:,1) to atom(:,i)
        type(am_class_point_group) :: rotg ! local point group as seen by rotating the shell
        type(am_class_point_group) :: revg ! reversal group of a typical bond in the shell v
        type(am_class_point_group) :: stab ! stabilizer of a typical bond in the shell v
    end type am_shell_cell

contains

    ! i/o

    subroutine      write_uc(dtv, unit, iotype, v_list, iostat, iomsg)
        !
        implicit none
        !
        class(am_class_unit_cell), intent(in) :: dtv
        integer     , intent(in)    :: unit
        character(*), intent(in)    :: iotype
        integer     , intent(in)    :: v_list(:)
        integer     , intent(out)   :: iostat
        character(*), intent(inout) :: iomsg
        !
        iostat = 0
        !
        if (iotype.eq.'LISTDIRECTED') then
            write(unit,'(a/)') '<unit_cell>'
                ! non-allocatable
                #:for ATTRIBUTE in ['bas','recbas','natoms','lattice_id']
                    $:write_xml_attribute_nonallocatable(ATTRIBUTE)
                #:endfor
                ! allocatable        
                #:for ATTRIBUTE in ['tau_frac','tau_cart','Z','uc_id','pc_id','ic_id']
                    $:write_xml_attribute_allocatable(ATTRIBUTE)
                #:endfor
                ! type-specific stuff
                select type (dtv)
                class is (am_class_irre_cell)
                    ! nested-allocatable object
                    #:for ATTRIBUTE in ['atom']
                        $:write_xml_attribute_allocatable_derivedtype(ATTRIBUTE)
                    #:endfor
                class is (am_shell_cell)
                    ! non-allocatable
                    #:for ATTRIBUTE in ['center','i','j','m','n']
                        $:write_xml_attribute_nonallocatable(ATTRIBUTE)
                    #:endfor
                    ! allocatable
                    #:for ATTRIBUTE in ['pg_id']
                        $:write_xml_attribute_allocatable(ATTRIBUTE)
                    #:endfor
                    ! nested objects
                    #:for ATTRIBUTE in ['rotg','revg','stab']
                        write(unit,*) dtv%${ATTRIBUTE}$
                    #:endfor
                class default
                    ! do nothing
                end select
            write(unit,'(a/)') '</unit_cell>'
        else
            stop 'ERROR [write_uc]: iotype /= LISTDIRECTED'
        endif
    end subroutine  write_uc

    subroutine      read_uc(dtv, unit, iotype, v_list, iostat, iomsg)
        !
        implicit none
        !
        class(am_class_unit_cell), intent(inout) :: dtv
        integer     , intent(in)    :: unit
        character(*), intent(in)    :: iotype
        integer     , intent(in)    :: v_list(:)
        integer     , intent(out)   :: iostat
        character(*), intent(inout) :: iomsg
        logical :: isallocated
        integer :: dims_rank
        integer :: dims(10) ! read tensor up to rank 5
        !
        if (iotype.eq.'LISTDIRECTED') then
            read(unit=unit,fmt='(/)',iostat=iostat,iomsg=iomsg)
                ! non-allocatable
                #:for ATTRIBUTE in ['bas','recbas','natoms','lattice_id']
                    $:read_xml_attribute_nonallocatable(ATTRIBUTE)
                #:endfor
                ! allocatable        
                #:for ATTRIBUTE in ['tau_frac','tau_cart','Z','uc_id','pc_id','ic_id']
                    $:read_xml_attribute_allocatable(ATTRIBUTE)
                #:endfor
                ! type-specific stuff
                select type (dtv)
                class is (am_class_irre_cell)
                    ! nested-allocatable object
                    #:for ATTRIBUTE in ['atom']
                       $:read_xml_attribute_allocatable_derivedtype(ATTRIBUTE)
                    #:endfor
                class is (am_shell_cell)
                    ! non-allocatable
                    #:for ATTRIBUTE in ['center','i','j','m','n']
                        $:read_xml_attribute_nonallocatable(ATTRIBUTE)
                    #:endfor
                    ! allocatable
                    #:for ATTRIBUTE in ['pg_id']
                        $:read_xml_attribute_allocatable(ATTRIBUTE)
                    #:endfor
                    ! nested objects
                    #:for ATTRIBUTE in ['rotg','revg','stab']
                        read(unit,*) dtv%${ATTRIBUTE}$
                    #:endfor
                class default
                    ! do nothing
                end select
            read(unit=unit,fmt='(/)',iostat=iostat,iomsg=iomsg)
            ! without the iostat=-1 here the following error is produced at compile time:
            ! tb(67203,0x7fff7e4dd300) malloc: *** error for object 0x10c898cec: pointer being freed was not allocated
            ! *** set a breakpoint in malloc_error_break to debug
            iostat=-1
        else
            stop 'ERROR [read_uc]: iotype /= LISTDIRECTED'
        endif
    end subroutine  read_uc

    subroutine      save_uc(uc,fname)
        !
        implicit none
        !
        class(am_class_unit_cell), intent(in) :: uc
        character(*), intent(in) :: fname
        integer :: fid
        integer :: iostat
        ! fid
        fid = 1
        ! save space group
        open(unit=fid, file=trim(fname), status='replace', action='write', iostat=iostat)
            if (iostat/=0) stop 'ERROR [uc:load]: opening file'
            write(fid,*) uc
        close(fid)
        !
    end subroutine  save_uc

    subroutine      load_uc(uc,fname)
        !
        implicit none
        !
        class(am_class_unit_cell), intent(inout) :: uc
        character(*), intent(in) :: fname
        integer :: fid
        integer :: iostat
        ! fid
        fid = 1
        ! clock in
        call start_clock('load')
        ! save space group
        open(unit=fid, file=trim(fname), status='old', action='read', iostat=iostat)
            if (iostat/=0) stop 'ERROR [uc:load]: opening file'
            read(fid,*) uc
        close(fid)
        ! clock out
        call stop_clock('load')
        !
    end subroutine  load_uc

    ! unit cell

    subroutine      write_poscar(uc,file_output_poscar)
        ! 
        ! Reads the poscar file, look below: code is short and self explanatory.
        !
        use am_vasp_io, write_poscar_internal => write_poscar
        !
        implicit none
        !
        class(am_class_unit_cell), intent(in) :: uc
        character(*), intent(in) :: file_output_poscar
        !
        integer, allocatable :: unique_Z(:)
        integer :: nspecies
        character(:), allocatable :: symb(:)
        integer, allocatable :: atype(:)
        integer :: i
        !
        if (uc%natoms.lt.1) then
            call am_print('ERROR','Number of atoms < 1.',flags='E')
            stop
        endif
        !
        ! nspecies
        unique_Z = unique(uc%Z)
        nspecies = size(unique_Z)
        if (nspecies.lt.1) then
            call am_print('ERROR','Number of atomic species < 1.',flags='E')
            stop
        endif
        !
        ! symbs
        allocate(character(500)::symb(nspecies))
        do i = 1, nspecies
            symb(i) = atm_symb(unique_Z(i))
        enddo
        !
        ! atype
        allocate(atype(uc%natoms))
        do i = 1, nspecies
            where(uc%Z(:).eq.unique_Z(i)) atype = i
        enddo
        !
        call write_poscar_internal(bas=uc%bas, natoms=uc%natoms, nspecies=nspecies, symb=symb, &
            & tau_frac=uc%tau_frac, atype=atype, header='', fname=file_output_poscar, verbosity=1)
        !
    end subroutine  write_poscar

    subroutine      get_supercell(sc,uc,bscfp,opts)
        !
        ! bscfp = supercell basis in fractional primitive lattice coordinates
        ! for primitive -> conventional cell, use: bscfp = inv(centering_matrix)
        !
        implicit none
        !
        class(am_class_unit_cell), intent(out) :: sc
        class(am_class_unit_cell), intent(in)  :: uc
        real(dp), intent(in) :: bscfp(3,3)
        real(dp) :: inv_bscfp(3,3)
        real(dp) :: tau_frac_wrk(3)
        integer :: i1, i2, i3, j, m, i
        type(am_class_options), intent(in) :: opts
        !
        if (opts%verbosity.ge.1) call print_title('Expanding to supercell')
        !
        inv_bscfp = inv(bscfp)
        sc%bas    = matmul(uc%bas,bscfp)
        sc%recbas = inv(sc%bas)
        !
        if ((mod(det(bscfp)+opts%prec,1.0_dp)-opts%prec).gt.opts%prec) then
            call am_print('ERROR','The determinant of the supercell basis in primitive fractional coordinates must be an integer.',flags='E')
            call am_print('bscfp',bscfp)
            call am_print('det(bscfp)',det(bscfp))
            call am_print('check',modulo(det(bscfp),1.0_dp).lt.tiny)
            write(*,*) opts%prec
            write(*,*) mod(det(bscfp),1.0_dp)
            write(*,*) mod(det(bscfp)+opts%prec,1.0_dp)-opts%prec
            stop
        endif
        !
        sc%natoms = nint(det(bscfp)*uc%natoms)
        !
        allocate(sc%tau_frac(3,sc%natoms))
        allocate(sc%Z(sc%natoms))
        allocate(sc%uc_id(sc%natoms))
        !
        sc%tau_frac = 1.0D5
        m=0
        map : do i1 = 0, nint(sum(abs(bscfp(:,1)))) ! sc%natoms ! shouldn't this be column vectors here?
              do i2 = 0, nint(sum(abs(bscfp(:,2)))) ! sc%natoms ! shouldn't this be column vectors here?
              do i3 = 0, nint(sum(abs(bscfp(:,3)))) ! sc%natoms ! shouldn't this be column vectors here?
                  do j = 1, uc%natoms
                      !
                      tau_frac_wrk = uc%tau_frac(1:3,j)
                      ! atomic basis in primitive fractional
                      tau_frac_wrk = tau_frac_wrk+real([i1,i2,i3],dp)
                      ! convert to supercell fractional
                      tau_frac_wrk = matmul(inv_bscfp,tau_frac_wrk)
                      ! reduce to primitive supercell
                      tau_frac_wrk = modulo(tau_frac_wrk+opts%prec,1.0_dp) - opts%prec
                      !
                      if (.not.issubset(sc%tau_frac,tau_frac_wrk,opts%prec)) then
                          m = m+1 
                          sc%tau_frac(1:3,m) = tau_frac_wrk
                          sc%Z(m) = uc%Z(j)
                          sc%uc_id(m) = j
                      endif
                      if (m.eq.sc%natoms) then
                          exit map
                      endif
                  enddo
              enddo
              enddo
        enddo map
        ! correct basic rounding error
        where(abs(sc%tau_frac).lt.opts%prec) sc%tau_frac=0.0_dp
        ! get cartesaian atomic basis
        allocate(sc%tau_cart(3,uc%natoms))
        sc%tau_cart = matmul(sc%bas,sc%tau_frac)
        ! correct basic rounding error
        where(abs(sc%tau_cart).lt.opts%prec) sc%tau_cart=0.0_dp
        ! mapping of each atom in the sphere to the corresponding unit cell atom was already performed in the creation of the supercell
        ! now map each supercell atom the corresponding primitive atom if the mapping from uc->pc was already established
        if (allocated(uc%pc_id)) then
            allocate(sc%pc_id(sc%natoms))
            do i = 1, sc%natoms
                sc%pc_id(i) = uc%pc_id( sc%uc_id(i) )
            enddo
        endif
        ! now map each supercell atom the corresponding irreducible cell atom if the mapping from uc->ic was already established
        if (allocated(uc%ic_id)) then
            allocate(sc%ic_id(sc%natoms))
            do i = 1, sc%natoms
                sc%ic_id(i) = uc%ic_id( sc%uc_id(i) )
            enddo
        endif
        ! print stdout
        if (opts%verbosity.ge.1) then
            !
            call am_print_two_matrices_side_by_side(name='primitive basis',&
                Atitle='fractional (supercell)',A=inv_bscfp,&
                Btitle='cartesian'             ,B=uc%bas,&
                iopt_emph=flare,iopt_teaser=.true.)
            !
            call am_print_two_matrices_side_by_side(name='supercell basis',&
                Atitle='fractional (primitive)',A=bscfp,&
                Btitle='cartesian'             ,B=sc%bas,&
                iopt_emph=flare,iopt_teaser=.true.)
            !
            call am_print('number of atoms (supercell)',sc%natoms,flare)
            !
            call am_print('number of atoms (primitive)',uc%natoms,flare)
            !
            call am_print_two_matrices_side_by_side(name='primitive atomic basis',&
                Atitle='fractional (primitive)',A=transpose(uc%tau_frac),&
                Btitle='cartesian'             ,B=transpose(uc%tau_cart),&
                iopt_emph=flare,iopt_teaser=.true.)
            !
            call am_print_two_matrices_side_by_side(name='supercell atomic basis',&
                Atitle='fractional (supercell)',A=transpose(sc%tau_frac),&
                Btitle='cartesian'             ,B=transpose(sc%tau_cart),&
                iopt_emph=flare,iopt_teaser=.true.)
            !
            write(*,'(a,a)',advance='no') flare, 'atomic mapping (to intput: sc->uc)'
            call id_print_map(sc%uc_id)
            !
            if (allocated(uc%pc_id)) then
            write(*,'(a,a)',advance='no') flare, 'atomic mapping (to primitive: sc->pc)'
            call id_print_map(sc%ic_id)
            endif
            !
            if (allocated(uc%ic_id)) then
            write(*,'(a,a)',advance='no') flare, 'atomic mapping (to irreducible: sc->ic)'
            call id_print_map(sc%ic_id)
            endif
        endif
        !
    end subroutine  get_supercell

    pure subroutine copy(cp,uc)
        !
        implicit none
        !
        class(am_class_unit_cell), intent(inout) :: cp
        type(am_class_unit_cell) , intent(in) :: uc
        !
        cp%natoms = uc%natoms
        cp%bas    = uc%bas
        cp%recbas = uc%recbas
        allocate(cp%tau_frac, source=uc%tau_frac)
        allocate(cp%tau_cart, source=uc%tau_cart)
        allocate(cp%Z       , source=uc%Z)
        if (allocated(uc%uc_id)) allocate(cp%uc_id, source=uc%uc_id)
        if (allocated(uc%pc_id)) allocate(cp%pc_id, source=uc%pc_id)
        if (allocated(uc%ic_id)) allocate(cp%ic_id, source=uc%ic_id)
        !
    end subroutine  copy

    subroutine      filter(uc,ind)
        !
        implicit none
        !
        class(am_class_unit_cell), intent(inout) :: uc
        integer     , intent(in) :: ind(:)
        !
        ! count how many atoms should be passed
        uc%natoms = size(ind)
        ! filter stuff
        if (allocated(uc%tau_frac)) call dvfilter(values=uc%tau_frac, ind=ind)
        if (allocated(uc%tau_cart)) call dvfilter(values=uc%tau_cart, ind=ind)
        if (allocated(uc%Z))        call  ifilter(values=uc%Z       , ind=ind)
        if (allocated(uc%uc_id))    call  ifilter(values=uc%uc_id   , ind=ind)
        if (allocated(uc%pc_id))    call  ifilter(values=uc%pc_id   , ind=ind)
        if (allocated(uc%ic_id))    call  ifilter(values=uc%ic_id   , ind=ind)
        !
        contains
        subroutine     ifilter(values,ind)
            !
            implicit none
            !
            integer , allocatable, intent(inout) :: values(:)
            integer , intent(in)    :: ind(:)
            integer , allocatable   :: Z(:) ! temporary variable
            integer :: i, ninds
            !
            ninds = size(ind)
            ! setup temporary variable
            allocate(Z,source=values)
            ! clear original variable
            deallocate(values)
            ! reallocate space
            allocate(values(ninds))
            ! transfer content through filter
            do i = 1, ninds
                values(i) = Z(ind(i))
            enddo
            ! deallocate
            deallocate(Z)
        end subroutine ifilter
        subroutine     dvfilter(values,ind)
            !
            implicit none
            !
            real(dp), allocatable, intent(inout) :: values(:,:)
            integer , intent(in)    :: ind(:)
            real(dp), allocatable   :: Z(:,:) ! temporary variable
            integer :: i, ninds, d1
            !
            ninds = size(ind)
            d1    = size(values,1)
            ! setup temporary variable
            allocate(Z,source=values)
            ! clear original variable
            deallocate(values)
            ! reallocate space
            allocate(values(d1,ninds))
            ! transfer content through filter
            do i = 1, ninds
                values(:,i) = Z(:,ind(i))
            enddo
            ! deallocate
            deallocate(Z)
        end subroutine dvfilter
    end subroutine  filter

    subroutine      sort_atoms(uc,criterion,flags)
        !
        use am_rank_and_sort
        !
        implicit none
        !
        class(am_class_unit_cell), intent(inout) :: uc
        real(dp)    , intent(in) :: criterion(:)
        character(*), intent(in) :: flags
        integer , allocatable ::  inds(:)
        !
        if (size(criterion).ne.uc%natoms) stop 'ERROR [sort_atoms]: size(criterion) /= number of atoms?'
        !
        allocate(inds(uc%natoms))
        if     (index(flags,'descend').ne.0) then
            call rank(-criterion,inds)
        elseif (index(flags,'ascend' ).ne.0) then
            call rank(criterion,inds)
        else
            stop 'ERROR [sort_atoms]: ascend/descend?'
        endif
        !
        ! sort stuff
        uc%tau_frac = uc%tau_frac(:,inds)
        uc%tau_cart = uc%tau_cart(:,inds)
        uc%Z        = uc%Z(inds)
        if (allocated(uc%uc_id)) uc%uc_id = uc%uc_id(inds)
        if (allocated(uc%pc_id)) uc%pc_id = uc%pc_id(inds)
        if (allocated(uc%ic_id)) uc%ic_id = uc%ic_id(inds)
        !
    end subroutine  sort_atoms

    subroutine      initialize(uc,bas,tau_frac,Z,uc_id,pc_id,ic_id)
        !
        implicit none
        !
        class(am_class_unit_cell) , intent(out) :: uc
        real(dp), intent(in) :: bas(3,3)
        real(dp), intent(in) :: tau_frac(:,:)
        integer , intent(in) :: Z(:)
        integer , intent(in), optional :: uc_id(:)
        integer , intent(in), optional :: pc_id(:)
        integer , intent(in), optional :: ic_id(:)
        !
        uc%bas    = bas
        uc%recbas = inv(bas)
        !
        uc%natoms = size(tau_frac,2)
        !
        allocate(uc%tau_frac,source=tau_frac)
        allocate(uc%tau_cart,source=matmul(bas,tau_frac))
        !
        allocate(uc%Z,source=Z)
        !
        if (present(uc_id)) allocate(uc%uc_id,source=uc_id)
        if (present(pc_id)) allocate(uc%pc_id,source=pc_id)
        if (present(ic_id)) allocate(uc%ic_id,source=ic_id)
        !
    end subroutine  initialize

    subroutine      id_print_map(id)
        !
        implicit none
        !
        integer, intent(in) :: id(:)
        integer :: i
        !
        do i = 1, size(id)
            if (modulo(i,10).eq.1) then
                write(*,*)
                write(*,'(5x)',advance='no')
            endif
            write(*,'(a9)',advance='no') trim(int2char(i))//'->'//trim(int2char(id(i)))
        enddo
        write(*,*)
    end subroutine  id_print_map

    ! primitive cell

    subroutine      get_primitive_cell(pc,uc,opts)
        !
        use am_rank_and_sort, only : rank
        !
        implicit none
        !
        class(am_class_prim_cell), intent(inout) :: pc
        class(am_class_unit_cell), intent(inout) :: uc
        type(am_class_options)   , intent(in) :: opts
        character(:), allocatable :: str(:)
        real(dp) :: uc2prim(3,3)
        integer  :: i
        !
        if (opts%verbosity.ge.1) call print_title('Primitive cell')
        !
        ! get primitive cell basis from translations which leave unit cell invariant
        pc%bas = get_primitive_basis(bas=uc%bas,tau_frac=uc%tau_frac,Z=uc%Z,prec=opts%prec)
        ! check if singular
        if (abs(det(pc%bas)).lt.tiny) stop 'ERROR [get_primitive]: singular basis'
        ! get reciprocal basis
        pc%recbas = inv(pc%bas)
        ! get lattice id
        pc%lattice_id = get_lattice_id(bas=pc%bas,prec=opts%prec)
        ! get basis transformation (uc fractional to pc fractional)
        uc2prim = matmul(pc%recbas,uc%bas)
        ! reduce atoms to primitive cell
        pc%tau_frac = reduce_atoms_to_primitive_cell(uc2prim=uc2prim, uc_tau_frac=uc%tau_frac, prec=opts%prec)
        ! get cartesian coordiantes
        allocate(pc%tau_cart, source=matmul(pc%bas,pc%tau_frac))
        ! get numer of atoms
        pc%natoms = size(pc%tau_frac,2)
        ! get map: pc -> pc
        allocate(pc%pc_id(pc%natoms),source=[1:pc%natoms])
        ! get map: uc -> uc
        if (.not.allocated(uc%uc_id)) stop 'ERROR [get_primitive]: uc_id is not allocated'
        ! allocate(uc%uc_id(uc%natoms),source=[1:uc%natoms])
        ! get map: uc -> pc
        uc%pc_id = get_uc2pc_map(uc2prim=uc2prim, pc_natoms=pc%natoms, pc_tau_frac=pc%tau_frac, uc_natoms=uc%natoms, uc_tau_frac=uc%tau_frac, prec=opts%prec)
        ! get map: pc -> uc (same algorithm as used in irreducible pair shell)
        allocate(pc%uc_id(pc%natoms))
        pc%uc_id = 0
        do i = 1, uc%natoms
        if (pc%uc_id(uc%pc_id(i)).eq.0) then
            pc%uc_id(uc%pc_id(i)) = i
        endif
        enddo
        ! transfer Z
        allocate(pc%Z(pc%natoms))
        pc%Z = uc%Z(pc%uc_id)
        !
        ! print stdout
        if (opts%verbosity.ge.1) then
            allocate(character(2) :: str(pc%natoms))
            !
            write(*,'(a,a)') flare, 'lattice = '//trim(get_lattice_type(lattice_id=pc%lattice_id))
            !
            write(*,'(a,a)') flare, 'original basis ='
            call disp_indent()
            call disp(title='fractional (primitive)', X=matmul(pc%recbas,uc%bas), style='underline', fmt='f13.8',advance='no' , trim = 'no')
            call disp(title='cartesian'             , X=uc%bas                  , style='underline', fmt='f13.8',advance='yes', trim = 'no')
            !
            write(*,'(a,a)') flare, 'primitive basis ='
            call disp_indent()
            call disp(title='fractional (input)'    , X=matmul(uc%recbas,pc%bas), style='underline', fmt='f13.8',advance='no' , trim = 'no')
            call disp(title='cartesian'             , X=pc%bas                  , style='underline', fmt='f13.8',advance='yes', trim = 'no')
            !
            write(*,'(a,a,a)') flare, 'primitive cell atoms = ', tostring(pc%natoms)
            do i = 1, pc%natoms
                str(i) = atm_symb(uc%Z(pc%uc_id(i)))
            enddo
            call disp_indent()
            call disp(title='#'                     , X=[1:pc%natoms]           , style='underline', fmt='i7'   ,advance='no' , trim = 'yes')
            call disp(title='Z'                     , X=str                     , style='underline', fmt='a3'   ,advance='no' , trim = 'no')
            call disp(title='fractional (primitive)', X=transpose(pc%tau_frac)  , style='underline', fmt='f11.8',advance='no' , trim = 'no')
            call disp(title='cartesian'             , X=transpose(pc%tau_cart)  , style='underline', fmt='f11.8',advance='yes', trim = 'no')
            !
            write(*,'(a,a)',advance='no') flare, 'atomic mapping (to primitive: uc->pc)'
            call id_print_map(uc%pc_id)
            !
            write(*,'(a,a)',advance='no') flare, 'atomic mapping (to input: pc->uc)'
            call id_print_map(pc%uc_id)
            !
        endif
        !
        contains
        function       get_primitive_basis(bas,tau_frac,Z,prec) result(pc_bas)
            !
            implicit none
            !
            real(dp), intent(in) :: bas(3,3)
            real(dp), intent(in) :: tau_frac(:,:)
            integer , intent(in) :: Z(:)
            real(dp), intent(in) :: prec
            real(dp) :: pc_bas(3,3)
            integer , allocatable :: indices(:)
            real(dp), allocatable :: T(:,:)
            integer :: nTs
            integer :: i,j,k
            !
            ! get basis translations which could serve as primitive cell vectors
            T = translations_from_basis(tau=tau_frac, Z=Z, prec=prec, flags='prim')
            T = matmul(bas,T)
            nTs = size(T,2)
            ! sort primitive vectors based on magnitude (smallest last)
            allocate(indices(nTs))
            call rank(norm2(T,1),indices)
            T=T(1:3, indices(nTs:1:-1))
            ! select three primitive vectors which yield the smallest cell volume by ...
            ! ... selecting smallest vector as first vector
            i = nTs
            ! ... selecting smallest pritimive vector which is noncollinear to the first (non-zero cross product) as second vector
            do j = nTs,1,-1
                if ( any(cross_product(T(1:3,i),T(1:3,j)) .gt. tiny) ) exit
                enddo
            if (any([i,j].eq.0)) stop 'ERROR [get_primitive_basis]: No pair of non-collinear lattice vectors found'
            ! ... selecting smallest primitive vector which produces a non-zero cell volume as third vector
            do k = nTs,1,-1
                ! regularization added to prevent division by zero
                if ( abs( det( T(1:3,[i,j,k])+eye(3)*1.0D-14 ) ) .gt. tiny ) exit
            enddo
            if (any([i,j,k].eq.0)) stop 'ERROR [get_primitive_basis]: No primitive basis found'
            !
            pc_bas = T(1:3,[i,j,k])
            ! 
        end function   get_primitive_basis
        function       reduce_atoms_to_primitive_cell(uc2prim,uc_tau_frac,prec) result(pc_tau_frac)
            !
            implicit none
            !
            real(dp), intent(in) :: uc2prim(3,3)
            real(dp), intent(in) :: uc_tau_frac(:,:)
            real(dp), intent(in) :: prec
            real(dp), allocatable :: pc_tau_frac(:,:)
            !
            ! allocating space
            allocate(pc_tau_frac,source=uc_tau_frac)
            ! converting atoms to pc fractional coordinates
            pc_tau_frac = matmul(uc2prim, pc_tau_frac)
            ! reduce atoms to pc
            pc_tau_frac = modulo(pc_tau_frac+prec,1.0_dp)-prec
            ! get unique values
            pc_tau_frac = unique(pc_tau_frac,prec)
            !
        end function   reduce_atoms_to_primitive_cell
        function       get_uc2pc_map(uc2prim,pc_natoms,pc_tau_frac,uc_natoms,uc_tau_frac,prec) result(uc_pc_id)
            !
            implicit none
            !
            real(dp), intent(in) :: uc2prim(3,3)
            integer , intent(in) :: pc_natoms
            integer , intent(in) :: uc_natoms
            real(dp), intent(in) :: uc_tau_frac(:,:)
            real(dp), intent(in) :: pc_tau_frac(:,:)
            real(dp), intent(in) :: prec
            integer :: uc_pc_id(uc_natoms)
            real(dp):: uc_tau_prim_frac(3,uc_natoms)
            integer :: i,j
            !
            ! intiailize
            uc_pc_id=0
            ! get uc atoms in primitive fractional coordinates
            uc_tau_prim_frac = matmul(uc2prim,uc_tau_frac)
            ! reduce to primitive cell between [0,1)
            uc_tau_prim_frac = modulo(uc_tau_prim_frac+prec,1.0_dp)-prec
            ! make map
            do i = 1, pc_natoms
            do j = 1, uc_natoms
                if ( isequal(pc_tau_frac(:,i),uc_tau_prim_frac(:,j)) ) uc_pc_id(j) = i
            enddo
            enddo
            !
            if (any(uc_pc_id.eq.0)) stop 'Unit cell -> primitive cell map failed.'
            !
        end function   get_uc2pc_map
    end subroutine  get_primitive_cell

    ! conventional cell

    function       get_lattice_id(bas,prec) result(lattice_id)
        ! Symmetry and Condensed Matter Physics: A Computational Approach. 1 edition. Cambridge,
        ! UK; New York: Cambridge University Press, 2008. p 282, table 10.3; p 288, table 10.4
        !
        !  code  lattice type
        !     1  simple cubic cell                            
        !     2  body-centered cubic cell                     
        !     3  face-centered cubic cell                     
        !     4  hexagonal cell                               
        !     5  simple tetragonal cell                       
        !     6  body-centered tetragonal cell                
        !     7  rhombohedral (trigonal) cell                 
        !     8  simple orthorhombic cell                     
        !     9  body-centered orthorhombic cell              
        !    10  face-centered orthorhombic cell              
        !    11  base centered orthorhombic cell              
        !    12  simple monoclinic cell                       
        !    13  body-centered monoclinic cell                
        !    14  triclinic cell
        implicit none
        !
        real(dp), intent(in) :: bas(3,3) ! number of point group symmetries
        real(dp), intent(in) :: prec
        integer  :: lattice_id
        real(dp) :: M_tensor(3,3) ! metric tensor
        real(dp) :: M(6) ! metric tensor values (voigt notation)
        integer  :: n_unique_diags
        integer  :: n_unique_offdiags
        integer  :: n_zeros_off_diags
        ! get metric tensor
        M_tensor = matmul(transpose(bas),bas)
        ! voigt notation
        M(1) = M_tensor(1,1) ! diags
        M(2) = M_tensor(2,2) ! diags
        M(3) = M_tensor(3,3) ! diags
        M(4) = M_tensor(2,3) ! angle
        M(5) = M_tensor(1,3) ! angle
        M(6) = M_tensor(1,2) ! angle
        ! initialize lattice_id
        lattice_id = 0
        ! start with high symmetry stuff and work way down
        n_unique_diags    = count(unique_inds(M(1:3),prec).ne.0)
        n_unique_offdiags = count(unique_inds(M(4:6),prec).ne.0)
        n_zeros_off_diags = count(abs(M(4:6)).lt.prec)
        ! figure lattice type
        select case (n_unique_diags)
        case(1)
            select case(n_unique_offdiags)
            case(1)
                !    BRAV    CEN DI  OFF additional conditions
                ! 2  CUBIC   I   1   1   (multiply offdiagonals by -3 and compare: 1 unique total)
                ! 3  CUBIC   F   1   1   (multiply offdiagonals by 2 and compare: 1 unique total)
                ! 1  CUBIC   P   1   1   *three zeros
                ! 7  TRI     R   1   1   
                if     (count(unique_inds([M(1:3),-3.0_dp*M(4:6)],prec).ne.0).eq.1) then
                    lattice_id = 2
                elseif (count(unique_inds([M(1:3),+2.0_dp*M(4:6)],prec).ne.0).eq.1) then
                    lattice_id = 3
                elseif (n_zeros_off_diags.eq.3) then
                    lattice_id = 1
                else
                    lattice_id = 7
                endif
            case(2)
                !    BRAV    CEN DI  OFF additional conditions
                ! 6  TETRA   I   1   2   
                lattice_id = 6
            case(3)
                !    BRAV    CEN DI  OFF additional conditions
                ! 9  ORTH    I   1   3   
                lattice_id = 9
            end select
        case(2)
            select case(n_unique_offdiags)
            case(1)
                !    BRAV    CEN DI  OFF additional conditions
                ! 5  TETRA   P   2   1   *three zeros
                lattice_id = 5
            case(2)
                !    BRAV    CEN DI  OFF additional conditions
                ! 11 ORTH    ABC 2   2   *two zeros
                ! 4  HEX     P   2   2   *two zeros (multiply offdiagonals by 2 and compare: 3 unique total)
                ! 13 MONO    I   2   3   
                if     (n_zeros_off_diags.eq.2) then
                    if (count(unique_inds([M(1:3),+2.0_dp*M(4:6)],prec).ne.0).eq.3) then
                        lattice_id = 3
                    else
                        lattice_id = 2
                    endif
                else
                    lattice_id = 13
                endif
            end select
        case(3)
            select case(n_unique_offdiags)
            case(2)
                !    BRAV    CEN DI  OFF additional conditions
                ! 8  ORTH    P   3   1   *three zeros
                ! 12 MONO    P   3   2   *two zeros
                if     (n_zeros_off_diags.eq.3) then
                    lattice_id = 8
                else
                    lattice_id = 12
                endif
            case(3)
                !    BRAV    CEN DI  OFF additional conditions
                ! 14 TRI     P   3   3   
                ! 10 ORTH    F   3   3   
                stop 'not implemented'
            end select
        end select
        if (lattice_id.eq.0) stop 'ERROR [lattice_type]: unkown lattice type'
        !
    end function   get_lattice_id

    function       get_lattice_type(lattice_id) result(lattice_type)
        !
        implicit none
        !
        integer, intent(in) :: lattice_id
        character(:), allocatable :: lattice_type
        !
        select case(lattice_id)
            case(1);  allocate(lattice_type,source='simple cubic')
            case(2);  allocate(lattice_type,source='body-centered cubic')
            case(3);  allocate(lattice_type,source='face-centered cubic')
            case(4);  allocate(lattice_type,source='hexagonal')
            case(5);  allocate(lattice_type,source='simple tetragonal')
            case(6);  allocate(lattice_type,source='body-centered tetragonal')
            case(7);  allocate(lattice_type,source='rhombohedral (trigonal)')
            case(8);  allocate(lattice_type,source='simple orthorhombic')
            case(9);  allocate(lattice_type,source='body-centered orthorhombic')
            case(10); allocate(lattice_type,source='face-centered orthorhombic')
            case(11); allocate(lattice_type,source='base-centered orthorhombic')
            case(12); allocate(lattice_type,source='simple monoclinic')
            case(13); allocate(lattice_type,source='base-centered monoclinic')
            case(14); allocate(lattice_type,source='triclinic')
        case default
            stop 'ERROR [get_lattice_type]: Unknown lattice_id'
        end select
        !
    end function   get_lattice_type

    function       get_bravais_id(lattice_id) result(bravais_id)
        !
        implicit none
        !
        integer, intent(in) :: lattice_id
        integer :: bravais_id
        !
        select case(lattice_id)
        case(1:3);   bravais_id = 1 ! cubic 
        case(4);     bravais_id = 2 ! hexagonal
        case(5:6);   bravais_id = 3 ! tetragonal
        case(7);     bravais_id = 4 ! trigonal
        case(8:11);  bravais_id = 5 ! orthorhombic
        case(12:13); bravais_id = 6 ! monoclinic
        case(14);    bravais_id = 7 ! triclinic
        case default
            stop 'ERROR [get_bravais_id]: Unknown lattice_id'
        end select
    end function   get_bravais_id

    pure function  get_celldm(bas) result(celldm)
        !
        implicit none
        !
        real(dp), intent(in) :: bas(3,3)
        real(dp) :: M(3,3) ! metric tensor
        real(dp) :: celldm(6)
        real(dp) :: rad2deg
        !
        rad2deg = 180.0_dp/pi
        ! 
        M = matmul(transpose(bas),bas)
        !
        celldm(1) = sqrt(M(1,1))
        celldm(2) = sqrt(M(2,2))
        celldm(3) = sqrt(M(3,3))
        celldm(4) = acos(M(2,3)/celldm(2)*celldm(3))*rad2deg
        celldm(5) = acos(M(1,3)/celldm(1)*celldm(3))*rad2deg
        celldm(6) = acos(M(1,2)/celldm(1)*celldm(2))*rad2deg
        !
    end function   get_celldm

    ! irreducible cell

    subroutine     get_irreducible_cell(ic,pc,uc,sg,opts)
        !
        implicit none
        !
        class(am_class_irre_cell), intent(out)   :: ic ! irreducible cell
        class(am_class_prim_cell), intent(inout) :: pc
        class(am_class_unit_cell), intent(inout), optional :: uc ! if present, mapping from ic <-> uc is obtained
        type(am_class_space_group),intent(in) :: sg
        type(am_class_options)   , intent(in) :: opts
        character(:), allocatable :: str(:)
        integer, allocatable :: PM(:,:)
        logical, allocatable :: mask(:)
        integer, allocatable :: ind(:)
        integer :: i,j,k
        !
        if (opts%verbosity.ge.1) call print_title('Irreducible cell')
        !
        ! get primitive basis
        ic%bas = pc%bas
        ! get reciprocal basis
        ic%recbas = pc%recbas
        ! get permutation map which shows how atoms are permuted by each space symmetry operation, PM(pc%natoms,sg%nsyms) 
        PM = permutation_map(seitz=sg%seitz_frac, tau=pc%tau_frac, flags='', prec=opts%prec)
        ! determine irreducible atoms, i.e. get all atoms, in increments, which have not been already mapped onto by a space symmetry operation
        ! thus, primitive cell atom i corresponds to irreducible atom k ...
        allocate(mask(pc%natoms))
        mask = .true.
        ! ind(1:k) are indices on primitive atoms corresponding to irreducible atoms
        allocate(ind(pc%natoms)) 
        ind  = 0
        !
        k=0
        do i = 1, pc%natoms
        if (mask(i)) then
            k=k+1
            ind(k)=i
            do j = 1, sg%nsyms
                ! mask keeps track of pc atoms which were already mapped onto irreducible atoms
                mask(PM(i,j))=.false.
            enddo
        endif
        enddo
        ! get number of atoms
        ic%natoms = k
        ! transfer irreducible atoms (frac)
        allocate(ic%tau_frac,source=pc%tau_frac(:,ind(1:k)))
        ! transfer irreducible atoms (cart)
        allocate(ic%tau_cart,source=pc%tau_cart(:,ind(1:k)))
        ! transfer Z
        allocate(ic%Z,source=pc%Z(ind(1:k)))
        ! map irreducible atom -> irreducible atom
        allocate(ic%ic_id,source=[1:ic%natoms])
        ! map irreducible atom onto -> primitive atom
        allocate(ic%pc_id,source=ind(1:k))
        ! map irreducible atom onto -> unit atom
        allocate(ic%uc_id,source=pc%uc_id(ind(1:k)))
        ! map primitive atom onto -> irreducible atom
        allocate(pc%ic_id(pc%natoms))
        pc%ic_id = 0
        do i = 1, ic%natoms
            ! PM(1,:) shows all atoms onto which atom 1 is mapped by all space symmetry operations
            do j = 1, pc%natoms
                search : do k = 1, sg%nsyms
                    if (PM(j,k).eq.ic%pc_id(i)) then
                        pc%ic_id(j) = i
                        exit search
                    endif
                enddo search
            enddo
        enddo
        if (any(pc%ic_id.eq.0)) stop 'ERROR: pc->ic mapping failed.'
        ! maps (input) unit cell atom onto -> irreducible atom
        if (present(uc)) then
            allocate(uc%ic_id(uc%natoms))
            uc%ic_id = 0
            do i = 1, uc%natoms
                uc%ic_id(i) = pc%ic_id(uc%pc_id(i))
            enddo
        endif
        ! print stdout
        if (opts%verbosity.ge.1) then
            allocate(character(4) :: str(ic%natoms))
            !
            write(*,'(a,a,a)') flare, 'input atoms = ', tostring(uc%natoms)
            write(*,'(a,a,a)') flare, 'primitive atoms = ', tostring(pc%natoms)
            write(*,'(a,a,a)') flare, 'irreducible atoms = ', tostring(ic%natoms)
            do i = 1, ic%natoms
                str(i) = atm_symb(ic%Z(i))
            enddo
            call disp(title=' '                     , X=zeros([1,1])            , style='underline', fmt='i2'   ,advance='no' , zeroas=' ' )
            call disp(title=' id '                  , X=[1:ic%natoms]           , style='underline', fmt='i7'   ,advance='no' , trim = 'yes')
            call disp(title='atom'                  , X=str                     , style='underline', fmt='a3'   ,advance='no' , trim = 'no')
            call disp(title='fractional (primitive)', X=transpose(ic%tau_frac)  , style='underline', fmt='f11.8',advance='no' , trim = 'no')
            call disp(title='cartesian'             , X=transpose(ic%tau_cart)  , style='underline', fmt='f11.8',advance='yes', trim = 'no')
            !
            write(*,'(a,a)',advance='no') flare, 'atomic mapping (to input: ic->uc)'
            call id_print_map(ic%uc_id)
            !
            write(*,'(a,a)',advance='no') flare, 'atomic mapping (to primitive: ic->pc)'
            call id_print_map(ic%pc_id)
            !
            write(*,'(a,a)',advance='no') flare, 'atomic mapping (from primitive: pc->ic)'
            call id_print_map(pc%ic_id)
            !
            if (present(uc)) then
            write(*,'(a,a)',advance='no') flare, 'atomic mapping (from input: uc->ic)'
            call id_print_map(uc%ic_id)
            endif
            !
        endif
        !
    end subroutine get_irreducible_cell

    subroutine     initialize_orbitals(ic,opts)
        !
        ! Useful if no states on atoms are defined. Gives each atom all the possible states
        !
        ! Correspoding to the L shell (n=2) [ s p   states (1+3=4   states) ]
        ! Correspoding to the M shell (n=3) [ s p d states (1+3+5=9 states) ]
        !
        ! There will be 1 matrix element for every irreducible pair of orbitals
        !
        implicit none
        !
        class(am_class_irre_cell), intent(inout) :: ic
        type(am_class_options)   , intent(in) :: opts
        integer :: i
        character(100) :: fname
        !
        if (opts%verbosity.ge.1) call print_title('Atomic orbitals')
        !
        fname = 'infile.orbitals'
        if (fexists(fname)) then
            ! read it
            call read_orbital_basis(ic,fname,opts%verbosity)
        else
            if (opts%verbosity.ge.1) write(*,'(a,a)') flare, trim(fname)//' not found.'
            ! create basis on each irreducible atoms
            allocate(ic%atom(ic%natoms))
            do i = 1, ic%natoms
                call ic%atom(i)%gen_orbitals(orbital_flags='1s,2p,3d,4f')
            enddo
            call template_orbital_basis(ic)
            !
            if (opts%verbosity.ge.1)  then
                write(*,'(a,a)') flare, 'Template produced: template.orbitals'
                call print_title('Done!')
            endif
            !
            stop 
            !************************************************************************************
            ! END PROGRAM HERE IF INFILE.TB_ORBITALS DOES NOT EXIST!
            !************************************************************************************
        endif
        ! check for error
        do i = 1, ic%natoms
            if (ic%atom(i)%norbitals.le.0) then
                stop 'ERROR [initialize_orbitals]: no orbitals set on an atom. check input.'
            endif
        enddo
        !
        if (opts%verbosity.ge.1) call print_orbital_basis(ic)
        !
        contains
        subroutine     template_orbital_basis(ic)
            !
            implicit none
            !
            class(am_class_irre_cell), intent(inout) :: ic
            integer :: i, j
            integer :: fid
            ! export file
            fid = 1
            open(unit=fid,file='template.orbitals',status='replace',action='write')
                ! spin polarized?
                write(fid,'(a)') 'spin: off'
                ! irreducible atoms
                write(fid,'(a,a)') 'irreducible atoms: ', tostring(ic%natoms)
                ! loop over irreducible atoms
                do i = 1, ic%natoms
                    write(fid,'(i10)',advance='no') i
                    ! loop over azimuthal quantum numbers
                    do j = 1, ic%atom(i)%nazimuthals
                        write(fid,'(a)',advance='no') ' '//trim(ic%atom(i)%orbname(j))
                    enddo
                    write(fid,*)
                enddo
            close(fid)
            !
        end subroutine template_orbital_basis
        subroutine     read_orbital_basis(ic,fname,verbosity)
            !
            implicit none
            !
            class(am_class_irre_cell), intent(inout) :: ic
            character(*), intent(in) :: fname
            integer, intent(in) :: verbosity
            integer :: i,j,k
            integer :: fid
            integer :: iostat
            character(maximum_buffer_size) :: buffer ! read buffer
            character(len=:), allocatable :: word(:) ! read buffer
            integer :: ic_natoms, nazimuthals
            character(len=100) :: orbital_flags
            character(len=100) :: spin_polarized_flag
            !
            fid = 1
            open(unit=fid,file=trim(fname),status="old",action='read')
                !
                if (verbosity.ge.1) write(*,'(a,a,a)') flare, 'atomic orbitals read = ', trim(fname)
                ! spin polarized?
                read(unit=fid,fmt='(a)') buffer
                word = strsplit(buffer,delimiter=' ')
                read(word(2),*) spin_polarized_flag
                ! irreducible atoms
                read(unit=fid,fmt='(a)') buffer
                word = strsplit(buffer,delimiter=' ')
                read(word(3),*) ic_natoms
                if (verbosity.ge.1) write(*,'(a,a,a)') flare, 'irreducible atoms = ', tostring(ic_natoms)
                ! set irreducible atoms
                if (ic%natoms.ne.ic_natoms) stop 'number of irreducible atoms input does not match internally calculated.'
                if (allocated(ic%atom)) deallocate(ic%atom)
                allocate(ic%atom(ic_natoms))
                ! loop over irreducible atoms
                do j = 1, ic_natoms
                    read(unit=fid,fmt='(a)') buffer
                    word = strsplit(buffer,delimiter=' ')
                    read(word(1),*) i
                    !
                    nazimuthals = size(word) - 1
                    if (verbosity.ge.1) write(*,'(a,a,a,a)', advance='no') flare, 'atom ', tostring(j), ' azimuthals ('//tostring(nazimuthals)//'):'
                    if (nazimuthals.le.0) stop 'number of azimuthals < 0'
                    !
                    orbital_flags=trim(spin_polarized_flag)
                    do k = 1, nazimuthals
                        orbital_flags=trim(orbital_flags)//','//trim(word(k+1))
                        if (verbosity.ge.1) write(*,'(a)',advance='no') ' '//trim(word(k+1))
                    enddo
                    if (verbosity.ge.1) write(*,*)
                    !
                    call ic%atom(i)%gen_orbitals(orbital_flags=orbital_flags)
                    !
                enddo
            close(fid)
            !
        end subroutine read_orbital_basis
        subroutine     print_orbital_basis(ic)
            !
            implicit none
            !
            class(am_class_irre_cell) , intent(inout) :: ic
            integer :: n,l,m,s
            integer :: i,j,a
            !
            do i = 1, ic%natoms
                write(*,'(a,a)') flare, 'irreducible atom '//tostring(i)//' contributes '//tostring((ic%atom(i)%norbitals))//' orbitals |n,l,m,s> :'
                do a = 1, ic%atom(i)%nazimuthals
                    write(*,'(5x,a,a)') trim(ic%atom(i)%orbname(a)), ':'
                do j = 1, ic%atom(i)%norbitals
                    n = ic%atom(i)%orbital(1,j)
                    l = ic%atom(i)%orbital(2,j)
                    m = ic%atom(i)%orbital(3,j)
                    s = ic%atom(i)%orbital(4,j)
                    if (ic%atom(i)%azimuthal(a).eq.l) write(*,'(5x,a2,i2,a1,i2,a1,i2,a1,i2,a2)') ' |', n, ',', l, ',', m, ',', s, '> '
                enddo
                enddo
            enddo
        end subroutine print_orbital_basis
    end subroutine initialize_orbitals

    ! deform (strain related stuff)

    subroutine     deform(def,uc,opts)
        !
        implicit none
        !
        type(am_class_unit_cell), allocatable, intent(out) :: def(:)
        type(am_class_unit_cell), intent(in) :: uc
        type(am_class_options)  , intent(in) :: opts
        real(dp), allocatable :: bas_def(:,:,:)
        integer :: i
        !
        if (opts%verbosity.ge.1) call print_title('Elastically deforming structure')
        !
        bas_def = apply_elastic_deformation(bas=uc%bas,deformation_code=opts%deformation_code,strain_max=opts%maxstrain,nstrains=opts%nstrains,iopt_verbosity=opts%verbosity)
        !
        allocate(def(opts%nstrains))
        !
        do i = 1,opts%nstrains
            !
            call def(i)%copy(uc)
            !
            def(i)%bas = bas_def(:,:,i)
            !
        enddo
        !
    end subroutine deform

    function       apply_elastic_deformation(bas,deformation_code,strain_max,nstrains,iopt_verbosity) result(bas_def)
        !
        !
        ! needs editing to conform. make this a "pure function" , push all outputs except errors into its caller.
        !
        ! creates at bas_def(3,3,nTs) matrix containing the nTs deformed primitive basises.
        !
        ! bas_def = apply_elastic_deformation(bas=uc%bas,deformation_code=6,strain_max=0.1_dp,nstrains=10,iopt_verbosity=verbosity)
        !
        implicit none
        !
        real(dp), intent(in) :: bas(3,3)
        integer , intent(in) :: deformation_code
        real(dp), intent(in) :: strain_max
        integer , intent(in) :: nstrains
        real(dp), allocatable :: bas_def(:,:,:)
        !
        character(len=6) :: dc ! deformation code string
        integer :: i, j ! loop variable
        real(dp) :: eta_step ! the deformation step size
        real(dp) :: eta ! the magnitude of the deformation
        real(dp) :: e(6) ! strain voigt
        real(dp) :: eta_matrix(3,3) ! deformation matrix (lagrangian)
        real(dp) :: eps_matrix(3,3) ! deformation matrix (true)
        real(dp), allocatable :: def_mat(:,:,:) ! deformation matrix
        !
        integer, intent(in), optional :: iopt_verbosity
        integer :: verbosity
        if ( present(iopt_verbosity) ) then
            verbosity = iopt_verbosity
        else
            verbosity = 1
        endif
        !
        if (nstrains .le. 1) then
            call am_print('ERROR','More than one strain point is required.',flags='E')
            stop
        endif
        !
        if (strain_max .gt. 0.1) then
            call am_print('ERROR','Strain value too large.',flags='E')
            stop
        endif
        !
        dc = 'XXXXXX'
        select case(deformation_code)
            case(0);  dc='EEE000'
            case(1);  dc='E00000'
            case(2);  dc='0E0000'
            case(3);  dc='00E000'
            case(4);  dc='000E00'
            case(5);  dc='0000E0'
            case(6);  dc='00000E'
            case(7);  dc='000EEE'
            case(8);  dc='EE0000'
            case(9);  dc='Ee0000'
            case(10); dc='EEEEEE'
            case(11); dc='123456'
            case(12); dc='d14t6c'
            case(13); dc='3cu62q'
            case(14); dc='qs51t2'
            case(15); dc='546dut'
            case(16); dc='s3d541'
            case(20); dc='000EE0'
            case(40); dc='EEd000'
            case(41); dc='Ee0000'
            case default
            call am_print('ERROR','Deformation code not found. Select one of the following:',flags='E')
            write(*,'(5x,a)') "  0 =>  ( eta,  eta,  eta,    0,    0,    0)  | volume strain "
            write(*,'(5x,a)') "  1 =>  ( eta,    0,    0,    0,    0,    0)  | linear strain along x "
            write(*,'(5x,a)') "  2 =>  (   0,  eta,    0,    0,    0,    0)  | linear strain along y "
            write(*,'(5x,a)') "  3 =>  (   0,    0,  eta,    0,    0,    0)  | linear strain along z "
            write(*,'(5x,a)') "  4 =>  (   0,    0,    0,  eta,    0,    0)  | yz shear strain"
            write(*,'(5x,a)') "  5 =>  (   0,    0,    0,    0,  eta,    0)  | xz shear strain"
            write(*,'(5x,a)') "  6 =>  (   0,    0,    0,    0,    0,  eta)  | xy shear strain 2C44"
            write(*,'(5x,a)') "  7 =>  (   0,    0,    0,  eta,  eta,  eta)  | shear strain along (111)"
            write(*,'(5x,a)') "  8 =>  ( eta,  eta,    0,    0,    0,    0)  | xy in-plane strain "
            write(*,'(5x,a)') "  9 =>  ( eta, -eta,    0,    0,    0,    0)  | xy in-plane shear strain"
            write(*,'(5x,a)') " 10 =>  ( eta,  eta,  eta,  eta,  eta,  eta)  | global strain"
            write(*,'(5x,a)') " 11 =>  (  1e,   2e,   3e,   4e,   5e,   6e)  | universal strain 1"
            write(*,'(5x,a)') " 12 =>  ( -2e,   1e,   4e,  -3e,   6e,  -5e)  | universal strain 2"
            write(*,'(5x,a)') " 13 =>  (  3e,  -5e,  -1e,   6e,   2e,  -4e)  | universal strain 3"
            write(*,'(5x,a)') " 14 =>  ( -4e,  -6e,   5e,   1e,  -3e,   2e)  | universal strain 4"
            write(*,'(5x,a)') " 15 =>  (  5e,   4e,   6e,  -2e,  -1e,  -3e)  | universal strain 5"
            write(*,'(5x,a)') " 16 =>  ( -6e,   3e,  -2e,   5e,  -4e,   1e)  | universal strain 6"
            write(*,'(5x,a)') " 20 =>  (   0,    0,    0,  eta,  eta,    0)  | yz zx shear"
            write(*,'(5x,a)') " 40 =>  ( eta,  eta,  -2e,    0,    0,    0)  | tetragonal for 3/2(c11-c12)"
            write(*,'(5x,a)') " 41 =>  ( eta, -eta,    0,    0,    0,    0)  | orthorhombic for (c11-c12)"
            stop
        end select
        !
        eta_step = 2.0_dp*strain_max/(nstrains-1.0_dp)
        !
        allocate(bas_def(3,3,nstrains))
        allocate(def_mat(3,3,nstrains))
        !
        if (verbosity.ge.1) call am_print('Number of deformations',nstrains,flare)
        !
        do i = 1, nstrains
            ! set deformation magnitude
            eta = (i-1)*eta_step - strain_max
            if (abs(eta) < 1.0E-6_dp) eta=1.0E-6_dp ! why is this here?
            ! convert string to deformation
            e(1:6)=0.0_dp
            do j = 1,6
                select case (dc(j:j))
                    case('E'); e(j)= eta
                    case('e'); e(j)=-eta
                    case('0'); e(j)= 0.0_dp
                    case('1'); e(j)= 1.0_dp*eta
                    case('2'); e(j)= 2.0_dp*eta
                    case('3'); e(j)= 3.0_dp*eta
                    case('4'); e(j)= 4.0_dp*eta
                    case('5'); e(j)= 5.0_dp*eta
                    case('6'); e(j)= 6.0_dp*eta
                    case('u'); e(j)=-1.0_dp*eta
                    case('d'); e(j)=-2.0_dp*eta
                    case('t'); e(j)=-3.0_dp*eta
                    case('q'); e(j)=-4.0_dp*eta
                    case('c'); e(j)=-5.0_dp*eta
                    case('s'); e(j)=-6.0_dp*eta
                    case default
                        call am_print('ERROR','Deformation not allowed',flags='E')
                        stop
                end select
            enddo
            ! build deformation matrix
            eta_matrix(1:3,1) = [e(1)       ,e(6)/2.0_dp,e(5)/2.0_dp]
            eta_matrix(1:3,2) = [e(6)/2.0_dp,e(2)       ,e(4)/2.0_dp]
            eta_matrix(1:3,3) = [e(5)/2.0_dp,e(4)/2.0_dp,e(3)       ]
            !
            ! Solve for eps (true strain) from eta (lagrangian strain) by iteration
            !
            ! "ElaStic: A tool for calculating second-order elastic constants from first principles"
            ! Eq. 2 : tau = det(1+eps)(1+eps)^-1 . sigma . (1+eps)^-1 ; dot indicates tensor product
            ! Eq. 1 : eta = eps + 1/2 * eps*eps
            ! Also see Wallace "Thermodynamics of crystals"
            !
            ! % Matlab script
            ! % maximum difference is 1E-17 for e < 0.1, which is more than enough strain for elastic constants
            ! % if more strain is desired, then a differen method for solving for eps is required.
            ! for j = 1:200
            ! e = rand(6,1)*0.1;
            ! eta = [e(1),e(6)/2,e(5)/2;
            !        e(6)/2,e(2),e(4)/2;
            !        e(5)/2,e(4)/2,e(3)];
            ! eps = eta;
            ! for i = 1:300
            !     % eta = eps + 1/2eps*eps => rearrange:  eps = eta - 1/2eps*eps
            !     eps = eta - 1/2*eps*eps;
            ! end
            ! hold on; plot(j,max((eps+ 1/2*eps*eps-eta)),'s')
            ! end
            !
            eps_matrix = eta_matrix
            solve_for_eps : do j = 1,200
                eps_matrix = eta_matrix - 0.5_dp*matmul(eps_matrix,eps_matrix)
                if ( all(abs(eps_matrix).gt.10.0_dp) ) stop 'ERROR [apply_elastic_deformation]: eps is diverging, try a smaller value.'
                if ( all(abs(eta-eps_matrix).lt.1E-15_dp) ) exit solve_for_eps
            enddo solve_for_eps
            !
            def_mat(1:3,1:3,i) = eye(3) + eps_matrix
            bas_def(1:3,1:3,i)=matmul(def_mat(1:3,1:3,i),bas)
            !
        enddo
        !
        if (verbosity.ge.1) then
            write(*,'(a,a)') flare, 'Deformation matrices (true strain, Voigt notation)'
            write(*,'(a,2x,6a16)') flare, 'e1', 'e2', 'e3', 'e4', 'e5', 'e6'
            do i = 1,nstrains
                write(*,'(5x,i2,6f16.10)') i,def_mat(1,1,i),def_mat(2,2,i),def_mat(3,3,i),def_mat(2,3,i),def_mat(1,3,i),def_mat(1,2,i)
            enddo
        endif
    end function   apply_elastic_deformation

end module am_unit_cell




