#:include "fypp_macros.fpp"
module am_tight_binding

    use dispmodule
    use am_constants
    use am_mkl
    use am_matlab
    use am_options
    use am_shells
    use am_symmetry
    use am_unit_cell
    use am_symmetry_tensor
    use am_symmetry_relations
    use am_dispersion
    use am_brillouin_zone
    use am_dft

    implicit none

    private

    type :: am_class_tightbinding_fitter
        integer  :: maxiter  ! maximum number of iterations  
        integer  :: nbands   ! number of bands
        integer  :: nkpts    ! numnbr of kpoints
        integer  :: nxs      ! number of parameters to fit
        integer  :: nrs      ! size of residual vector
        real(dp) :: eps      ! stopping criterion
        integer               :: skip_band
        integer , allocatable :: selector_shell(:)
        integer , allocatable :: selector_kpoint(:)
        integer , allocatable :: selector_x(:)
        real(dp), allocatable :: x(:)
        real(dp), allocatable :: r(:)
        real(dp)              :: rms
    end type am_class_tightbinding_fitter

    type, extends(am_class_dispersion) :: am_class_tightbinding_dispersion
        complex(dp), allocatable :: C(:,:,:) ! tight binding coefficients
    end type am_class_tightbinding_dispersion

    type, public :: am_class_tightbinding
        integer :: nshells                             ! how many shells irreducible atoms
        integer :: nVs                                 ! number of irreducible matrix element values
        real(dp), allocatable :: V(:)                  ! their values
        integer , allocatable :: V_ind(:,:)            ! their indices: [ip_id, alpha, beta]
        type(am_class_tensor), allocatable :: tens(:)  ! tens(nshells) describes symmetry-adapted matrix elementes (corresponds to irreducible pair shells)
        type(am_class_tightbinding_pointgroup) :: pg   ! point group in tight binding representation
        type(am_class_tightbinding_dispersion) :: dr   ! band dispersion computed using the tight binding model
        type(am_class_tightbinding_fitter)     :: ft   ! fitter for optimizing tight binding matrix elements
        contains
        procedure :: initialize_tb
        procedure :: get_hamiltonian
        procedure :: get_dispersion
        procedure :: initialize_ft
        procedure :: optimize_matrix_element
        procedure :: read_irreducible_matrix_element
        procedure :: write_irreducible_matrix_element
        procedure :: write_all_matrix_element
        procedure :: set_matrix_element
        procedure :: export_to_matlab
        procedure :: save => save_tb
        procedure :: load => load_tb
        procedure, private :: write_tb
        procedure, private :: read_tb
        generic :: write(formatted) => write_tb
        generic :: read(formatted)  => read_tb
    end type am_class_tightbinding

contains

    ! i/o

    subroutine      write_tb(dtv, unit, iotype, v_list, iostat, iomsg)
        !
        implicit none
        !
        class(am_class_tightbinding), intent(in) :: dtv
        integer     , intent(in)    :: unit
        character(*), intent(in)    :: iotype
        integer     , intent(in)    :: v_list(:)
        integer     , intent(out)   :: iostat
        character(*), intent(inout) :: iomsg
        !
        iostat = 0
        !
        if (iotype.eq.'LISTDIRECTED') then
            write(unit,'(a/)') '<tb>'
                ! non-allocatable
                #:for ATTRIBUTE in ['nshells','nVs']
                    $:write_xml_attribute_nonallocatable(ATTRIBUTE)
                #:endfor
                ! allocatable        
                #:for ATTRIBUTE in ['V','V_ind']
                    $:write_xml_attribute_allocatable(ATTRIBUTE)
                #:endfor
                ! nested objects
                write(unit,*) dtv%pg
                ! write(unit,*) dtv%dr
                ! write(unit,*) dtv%ft
                ! nested-allocatable object
                #:for ATTRIBUTE in ['tens']
                    $:write_xml_attribute_allocatable_derivedtype(ATTRIBUTE)
                #:endfor
            write(unit,'(a/)') '</tb>'
        else
            stop 'ERROR [write_tb]: iotype /= LISTDIRECTED'
        endif
    end subroutine  write_tb

    subroutine      read_tb(dtv, unit, iotype, v_list, iostat, iomsg)
        !
        implicit none
        !
        class(am_class_tightbinding), intent(inout) :: dtv
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
                #:for ATTRIBUTE in ['nshells','nVs']
                    $:read_xml_attribute_nonallocatable(ATTRIBUTE)
                #:endfor
                ! allocatable        
                #:for ATTRIBUTE in ['V','V_ind']
                    $:read_xml_attribute_allocatable(ATTRIBUTE)
                #:endfor
                ! nested objects
                read(unit,*) dtv%pg
                ! read(unit,*) dtv%dr
                ! read(unit,*) dtv%ft
                ! nested-allocatable object
                #:for ATTRIBUTE in ['tens']
                   $:read_xml_attribute_allocatable_derivedtype(ATTRIBUTE)
                #:endfor
            read(unit=unit,fmt='(/)',iostat=iostat,iomsg=iomsg)
            ! without the iostat=-1 here the following error is produced at compile time:
            ! tb(67203,0x7fff7e4dd300) malloc: *** error for object 0x10c898cec: pointer being freed was not allocated
            ! *** set a breakpoint in malloc_error_break to debug
            iostat=-1
        else
            stop 'ERROR [read_tb]: iotype /= LISTDIRECTED'
        endif
    end subroutine  read_tb

    subroutine      save_tb(tb,fname)
        !
        implicit none
        !
        class(am_class_tightbinding), intent(in) :: tb
        character(*), intent(in) :: fname
        integer :: fid
        integer :: iostat
        ! fid
        fid = 1
        ! save space group
        open(unit=fid, file=trim(fname), status='replace', action='write', iostat=iostat)
            if (iostat/=0) stop 'ERROR [tb:load]: opening file'
            write(fid,*) tb
        close(fid)
        !
    end subroutine  save_tb

    subroutine      load_tb(tb,fname)
        !
        implicit none
        !
        class(am_class_tightbinding), intent(inout) :: tb
        character(*), intent(in) :: fname
        integer :: fid
        integer :: iostat
        ! fid
        fid = 1
        ! save space group
        open(unit=fid, file=trim(fname), status='old', action='read', iostat=iostat)
            if (iostat/=0) stop 'ERROR [tb:load]: opening file'
            read(fid,*) tb
        close(fid)
        !
    end subroutine  load_tb

    ! tight binding

    subroutine     initialize_tb(tb,pg,ip,ic,pc,opts)
        !
        implicit none
        !
        class(am_class_tightbinding), intent(out) :: tb   ! tight binding parameters
        type(am_class_point_group)  , intent(in)  :: pg   ! seitz point group
        type(am_class_prim_cell)    , intent(in)  :: pc   ! primitive cell
        type(am_class_irre_cell)    , intent(in)  :: ic   ! irreducible cell
        type(am_class_irre_pair)    , intent(in)  :: ip   ! irreducible pairs
        type(am_class_options)      , intent(in)  :: opts
        logical, allocatable :: is_independent(:)
        integer :: sub(2)
        integer :: i,j,k
        !
        if (opts%verbosity.ge.1) call print_title('Symmetry-adapted tight-binding parameters')
        ! get point group in tight binding basis
        tb%pg = get_tightbinding_pointgroup(pg=pg, pc=pc, ic=ic)
        ! set number of shells
        tb%nshells = ip%nshells
        ! get symmeterized tight binding matrix elements
        allocate(tb%tens(tb%nshells))
        do k = 1, tb%nshells
            call tb%tens(k)%symmetrize(pg=pg, pc=pc, ic=ic, shell=ip%shell(k), opts=opts, property='irreducible tight binding shell '//tostring(k))
        enddo
        ! get number of independent (irreducible) matrix elements
        tb%nVs = 0
        do i = 1, tb%nshells
            tb%nVs = tb%nVs + count(get_independent(tb%tens(i)%relations)) 
        enddo
        ! allocate space for independent (irreducible) matrix elements V
        allocate(tb%V(tb%nVs))
        tb%V = 0
        ! set indices : V_ind( [i, subd2ind(dims=tb%tens(i)%dims,sub=[alpha,beta])], nVs)
        allocate(tb%V_ind(3,tb%nVs))
        tb%V_ind = 0
        k=0
        do i = 1, tb%nshells
            is_independent = get_independent(tb%tens(i)%relations)
            do j = 1, product(tb%tens(i)%dims)
            if (is_independent(j)) then
                k=k+1
                sub = ind2sub(dims=tb%tens(i)%dims, ind=j)
                tb%V_ind(1,k) = i      ! shell
                tb%V_ind(2,k) = sub(1) ! alpha
                tb%V_ind(3,k) = sub(2) ! beta
            endif
            enddo
        enddo
        ! write template for irreducible matrix elements
        if (fexists('infile.tb_matrix_elements_irreducible').ne.0) then 
            ! print title
            call print_title('Input tight-binding matrix elements')
            ! read matrix elements
            call tb%read_irreducible_matrix_element()
            ! display input file on stdout
            call disp(x=[0],zeroas=' ',advance='no')
            call disp(x=tb%V,title='V',style='underline',advance='no')
            call disp(x=tb%V_ind(1,:),title='shell',style='underline',advance='no')
            call disp(x=tb%V_ind(2,:),title='alpha',style='underline',advance='no')
            call disp(x=tb%V_ind(3,:),title='beta',style='underline',advance='yes')
        else
            ! set matrix elements at zero
            call tb%set_matrix_element(flags='zero')
            ! write template for irreducible matrix elements
            call tb%write_irreducible_matrix_element()
        endif
        !
        contains
        function     get_tightbinding_pointgroup(pg,pc,ic) result(tb_pg)
            !
            implicit none
            !
            type(am_class_point_group), intent(in)  :: pg   ! seitz point group (rev stab rot groups as well)
            type(am_class_prim_cell)  , intent(in)  :: pc   ! primitive cell
            type(am_class_irre_cell)  , intent(in)  :: ic   ! irreducible cell
            type(am_class_tightbinding_pointgroup)  :: tb_pg ! point symmetries in tight binding basis
            integer  :: i
            !
            ! get number of bases functions in representation (see below) ...
            tb_pg%nbases = 0
            ! ... and also determine subsections of rotations in Hamiltonian basis corresponding to each atom
            allocate(tb_pg%S(pc%natoms))
            allocate(tb_pg%E(pc%natoms))
            do i = 1, pc%natoms
                tb_pg%S(i)   = tb_pg%nbases + 1
                tb_pg%nbases = tb_pg%nbases + ic%atom(pc%ic_id(i))%norbitals
                tb_pg%E(i)   = tb_pg%nbases
            enddo
            ! number of symmetries
            tb_pg%nsyms = pg%nsyms
            ! generate intrinsic symmetries
            allocate(tb_pg%sym(tb_pg%nbases,tb_pg%nbases,tb_pg%nsyms))
            tb_pg%sym = 0
            ! determine rotation in the hamiltonian basis
            ! Nye, J.F. "Physical properties of crystals: their representation by tensors and matrices". p 133 Eq 7
            do i = 1, pg%nsyms
                ! convert rotation R to symmetry representation in tight-binding basis (direct sum of wigner matrices)
                tb_pg%sym(:,:,i) = ps2D(S=tb_pg%S, E=tb_pg%E, R_cart=pg%seitz_cart(1:3,1:3,i), pc=pc, ic=ic)
            enddo
            ! check that identity is first
            if (.not.isequal(tb_pg%sym(:,:,1),eye(tb_pg%nbases))) stop 'ERROR [get_tight_binding_point_group]: Identity is not first.'
            ! correct basic rounding errors
            call correct_rounding_error(tb_pg%sym)
            ! copy symmetry ids
            allocate(tb_pg%ps_id, source=pg%ps_id)
            ! copy classes
            allocate(tb_pg%cc%id, source=pg%cc%id)
            ! check that multiplication table is identical to pg%mt
            call tb_pg%get_multiplication_table()
            if (.not.isequal(tb_pg%mt%multab,pg%mt%multab)) stop 'ERROR [get_tight_binding_point_group]: Multiplication table mismatch.'
            ! get conjugacy clases
            call tb_pg%get_conjugacy_classes()
            ! get character table
            call tb_pg%get_character_table()
            !
        end function get_tightbinding_pointgroup
        function     ps2D(S,E,R_cart,pc,ic) result(Dsum)
            ! produces rotation which commutes with the entire Hamiltonian (useful building Hamiltonian and probably later for kpoints stuff too)
            implicit none
            !
            integer , intent(in) :: S(:)
            integer , intent(in) :: E(:)
            real(dp), intent(in) :: R_cart(3,3)
            type(am_class_prim_cell), intent(in) :: pc ! primitive cell
            type(am_class_irre_cell), intent(in) :: ic ! irreducible cell
            real(dp), allocatable :: Dsum(:,:)
            integer :: nbases
            integer :: i
            ! get hamil
            nbases = maxval(E)
            ! allocate and initialize space for rotations in Hamiltonin bais
            allocate(Dsum(nbases,nbases))
            Dsum = 0.0_dp
            ! construct rotation in the Hamiltonian basis
            do i = 1, pc%natoms
                Dsum(S(i):E(i), S(i):E(i)) = ps2tb(R_cart=R_cart, atom=ic%atom(pc%ic_id(i)) )
            enddo
            !
        end function ps2D
    end subroutine initialize_tb

    function       get_hamiltonian(tb,pp,kpt) result(H)
        ! 
        ! Get tight binding Hamiltonian at kpt.
        ! 
        implicit none
        !
        class(am_class_tightbinding), intent(in) :: tb     ! tight binding matrix elements
        type(am_class_prim_pair)    , intent(in) :: pp     ! primitive pairs
        real(dp)                    , intent(in) :: kpt(3) ! cart
        complex(dp), allocatable, target :: Hsub_target(:,:)
        real(dp)   , allocatable, target :: pg_target(:,:,:)
        integer    , allocatable :: selector_shell(:)
        real(dp)   , pointer     :: Dm(:,:), Dn(:,:)
        integer    , allocatable :: S(:),E(:)
        complex(dp), allocatable :: H(:,:)
        complex(dp), pointer     :: Hsub(:,:)
        integer :: m ! primitive atom 1 index 
        integer :: n ! primitive atom 2 index
        integer :: k ! shell (primitive)
        integer :: l ! shell (irreducible)
        integer :: p ! atoms
        integer :: ip_nshells
        ! get number of irreducile shell
        ip_nshells = maxval(abs(pp%ip_id(:)))
        ! selector
        if (allocated(tb%ft%selector_shell)) then
            allocate(selector_shell, source=tb%ft%selector_shell)
        else
            allocate(selector_shell, source=[1:ip_nshells])
        endif
        ! allocate space for vectors demarking start and end of Hamiltonian subsection
        allocate(S, source=tb%pg%S)
        allocate(E, source=tb%pg%E)
        ! allocate workspace for H subsection (initialized later)
        allocate(Hsub_target(tb%pg%nbases,tb%pg%nbases))
        ! allocate workspace for H subsection (initialized later)
        allocate(pg_target, source=tb%pg%sym)
        ! allocate and initialize
        allocate(H(tb%pg%nbases,tb%pg%nbases))
        H = cmplx(0,0,dp)
        ! construct Hamiltonian
        do l = 1, ip_nshells
        if ( any(selector_shell.eq.l) ) then
            do k = 1, pp%nshells
            if (abs(pp%ip_id(k)).eq.l) then
                ! primitive atom indicies
                m = pp%shell(k)%m
                n = pp%shell(k)%n
                ! compute bloch sum by loop over atoms in shell
                do p = 1, pp%shell(k)%natoms
                    ! set pointers
                    Hsub => Hsub_target(S(m):E(m), S(n):E(n))
                    Dm   =>   pg_target(S(m):E(m), S(m):E(m), pp%shell(k)%pg_id(p))
                    Dn   =>   pg_target(S(n):E(n), S(n):E(n), pp%shell(k)%pg_id(p))
                    ! get matrix elements (initialize Hsub)
                    Hsub = get_matrix_element(tb=tb, ip_id=pp%ip_id(k))
                    ! rotate matrix elements as needed to get from the tau_frac(:,1) => tau_frac(:,x)
                    Hsub = matmul(matmul(transpose(Dm), Hsub), Dn)
                    ! multiply exponential factor from Bloch sum
                    Hsub = Hsub * exp(-itwopi*dot_product(pp%shell(k)%tau_cart(1:3,p), kpt)) ! kpt [rec. cart.]
                    ! this pair's contribution to the Hamiltonian
                    H(S(m):E(m), S(n):E(n)) = H(S(m):E(m), S(n):E(n)) + Hsub
                enddo
            endif
            enddo
            !
            if (debug) then
                if (.not.ishermitian(H)) stop 'ERROR [get_hamiltonian]: H is not Hermitian!'
            endif
        endif
        enddo
    end function   get_hamiltonian

    subroutine     get_dispersion(tb,pp,bz)
        !
        implicit none
        !
        class(am_class_tightbinding), intent(inout) :: tb
        type(am_class_prim_pair)    , intent(in) :: pp
        type(am_class_bz)           , intent(in) :: bz
        integer    , allocatable :: selector_kpoint(:)
        complex(dp), allocatable :: H(:,:), V(:,:)
        real(dp)   , allocatable :: D(:)
        integer :: i
        ! selector
        if (allocated(tb%ft%selector_kpoint)) then
            allocate(selector_kpoint, source=tb%ft%selector_kpoint)
        else
            allocate(selector_kpoint, source=[1:bz%nkpts])
        endif
        ! tight binding coefficients
        if (.not.allocated(tb%dr%C)) allocate(tb%dr%C(tb%pg%nbases,tb%pg%nbases,bz%nkpts))
        ! eigenvalues
        if (.not.allocated(tb%dr%E)) allocate(tb%dr%E(tb%pg%nbases,bz%nkpts))
        ! loop over kpoints
        do i = 1, bz%nkpts
        if ( any(selector_kpoint.eq.i) ) then
            ! construct hamiltonian
            H = tb%get_hamiltonian(pp=pp, kpt=bz%kpt_cart(:,i))
            ! diagonalize hamiltonian
            call am_zheev(A=H,V=V,D=D)
            tb%dr%C(:,:,i) = V
            tb%dr%E(:,i) = D
        endif
        enddo
    end subroutine get_dispersion

    ! export to matlab

    subroutine     export_to_matlab(tb,ip,pp,flags)
        ! flags = symbolic/numerical cart/frac
        implicit none
        !
        class(am_class_tightbinding), intent(in) :: tb     ! tight binding matrix elements
        type(am_class_prim_pair)    , intent(in) :: pp     ! primitive pairs
        type(am_class_irre_pair)    , intent(in) :: ip     ! irreducible pairs
        character(*)                , intent(in) :: flags  ! symbolic/numerical cart/frac
        integer, allocatable :: Sv(:) ! start and end of irreducible vector
        integer, allocatable :: Ev(:) ! start and end of irreducible vector
        integer, allocatable :: S(:)  ! start and end of hamiltonian subsection
        integer, allocatable :: E(:)  ! start and end of hamiltonian subsection
        integer :: i,j,l,k,p,m,n,end,fid
        ! parameters
        character(100) :: H_fnc_name
        character(100) :: V_fnc_name
        character(100) :: str_H
        character(100) :: str_Dm
        character(100) :: str_Dn
        character(100) :: str_tau
        character(100) :: str_Vab
        !
        ! set function names
        if     (index(flags,'symbolic').ne.0) then
            H_fnc_name = 'get_H_symbolic'
        elseif (index(flags, 'numeric').ne.0) then
            H_fnc_name = 'get_H_numeric'
        else
            stop 'ERROR [export_to_matlab]: flags != symbolic/numeric'
        endif
        !
        if     (index(flags,'cart').ne.0) then
            H_fnc_name=trim(H_fnc_name)//'_cart'
        elseif (index(flags,'frac').ne.0) then
            H_fnc_name=trim(H_fnc_name)//'_frac'
        else
            stop 'ERROR [export_to_matlab]: flags != cart/frac'
        endif
        !
        V_fnc_name = 'getV'
        ! create abbreviations
        allocate(S, source=tb%pg%S)
        allocate(E, source=tb%pg%E)
        end = size(tb%pg%E)
        ! allocate start and end vectors for irreducible matrix elements
        allocate(Sv(ip%nshells))
        allocate(Ev(ip%nshells))
        j = 0
        do i = 1, ip%nshells
            j = j + 1
            Sv(i) = j
            j = j + count(get_independent(tb%tens(i)%relations)) - 1
            Ev(i) = j
        enddo
        ! export hamiltonian
        fid = 1
        ! create directory
        call execute_command_line ('mkdir -p '//trim(outfile_dir_tb))
        call execute_command_line( 'mkdir -p '//trim(outfile_dir_tb)//'/matlab')
        ! export file
        open(unit=fid,file=trim(outfile_dir_tb)//'/matlab/'//trim(H_fnc_name)//'.m',status='replace',action='write')
            write(fid,'(a,a,a)') 'function [H] = ', trim(H_fnc_name), '(pg,v,kpt,shell)'
            write(fid,'(a)') 'i2pi = 2*sqrt(-1)*pi;'
            write(fid,'(a)') 'H(1:'//tostring(E(end))//',1:'//tostring(E(end))//') = 0;'
            if (index(flags,'symb').ne.0) write(fid,'(a)') 'H = sym(H);'
                ! construct Hamiltonian
                do l = 1, ip%nshells
                write(fid,'(a)') 'if any(shell=='//tostring(l)//')'
                    do k = 1, pp%nshells
                    if (abs(pp%ip_id(k)).eq.l) then
                        ! compute bloch sum by loop over atoms in shell
                        do p = 1, pp%shell(k)%natoms
                            ! determine how to write stuff
                            m = pp%shell(k)%m
                            n = pp%shell(k)%n
                            str_H  =      'H('//tostring(S(m))//':'//tostring(E(m))//','//tostring(S(n))//':'//tostring(E(n))//')'
                            str_Dm = 'pg.sym('//tostring(S(m))//':'//tostring(E(m))//','//tostring(S(m))//':'//tostring(E(m))//','//tostring(pp%shell(k)%pg_id(p))//')'
                            str_Dn = 'pg.sym('//tostring(S(n))//':'//tostring(E(n))//','//tostring(S(n))//':'//tostring(E(n))//','//tostring(pp%shell(k)%pg_id(p))//')'
                            if     (index(flags,'cart').ne.0) then
                                str_tau= '['//tostring(pp%shell(k)%tau_cart(1:3,p),fmt='SP,f10.5')//']'
                            elseif (index(flags,'frac').ne.0) then
                                str_tau= '['//tostring(pp%shell(k)%tau_frac(1:3,p),fmt='SP,f10.5')//']'
                            endif
                            str_Vab= trim(V_fnc_name)//tostring(l)//'(v('//tostring(Sv(l))//':'//tostring(Ev(l))//'))'
                            if (pp%ip_id(k).lt.0) str_Vab=trim(str_Vab)//"'" ! adjoint
                            ! write stuff
                            write(fid,'(a)',advance='no') trim(str_H)       // ' = ' // trim(str_H)   // ' + '
                            write(fid,'(a)',advance='no') trim(str_Dm)//"'" // ' * ' // trim(str_Vab) // ' * ' // trim(str_Dn)
                            write(fid,'(a)',advance='no')     ' * exp(i2pi*dot(kpt,' // trim(str_tau) // '));'
                            write(fid,*)
                        enddo
                    endif
                    enddo
                write(fid,'(a)') 'end'
                enddo
            write(fid,'(a)') 'end'
            ! export symmetry relations for each irreducible pair (append to the same file)
            do i = 1, ip%nshells
                call export_relations2matlab(relations=tb%tens(i)%relations, dims=tb%tens(i)%dims, fnc_name=trim(V_fnc_name)//tostring(i), fid_append=fid, flags='append')
            enddo
        close(fid)
    end subroutine export_to_matlab

    ! matrix element stuff

    subroutine     set_matrix_element(tb,V,flags)
        ! flags = seq/rand/zero/selector/decompress
        implicit none
        !
        class(am_class_tightbinding)      , intent(inout)        :: tb
        real(dp)    , optional            , intent(in)           :: V(:)
        character(*), optional            , intent(in)           :: flags
        integer :: i, j, k, alpha, beta
        !
        ! set irreducible matrix elements
        if      (index(flags,   'seq').ne.0) then
            tb%V = [1:tb%nVs]
        elseif  (index(flags,  'rand').ne.0) then
            stop 'ERROR [set_matrix_element]: rand does not seem to be working here'
            do i = 1, tb%nVs
                tb%V(i) = rand()
            enddo
        elseif  (index(flags,  'zero').ne.0) then
            tb%V = 0
        elseif  (present(V)                ) then
            if (size(V).ne.tb%nVs) stop 'ERROR [set_matrix_element]: V /= nVs dimension mismatch'
            if (index(flags,'selector').ne.0) then
                if (.not.allocated(tb%ft%selector_x)) stop 'ERROR [set_matrix_element]: selector must be allocated'
                ! initialize as zero
                tb%V = 0
                ! copy matrix elements belonging to select shells only
                tb%V( tb%ft%selector_x ) = V
            else
                ! copy all matrix elements
                tb%V = V
            endif
        else
            stop 'ERROR [set_matrix_element]: unknown flag set_matrix_element_dummies'
        endif
        !
        ! clear irreducible matrix elements in each shell
        do k = 1, tb%nshells
            tb%tens(k)%V = 0
        enddo
        ! transfer irreducible matrix elements to each shell
        do i = 1, tb%nVs
            k     = tb%V_ind(1,i)
            alpha = tb%V_ind(2,i)
            beta  = tb%V_ind(3,i)
            j     = sub2ind(dims=tb%tens(k)%dims, sub=[alpha,beta])
            !
            tb%tens(k)%V(j) = tb%V(i)
        enddo
        ! once the irreducible matrix elements have been copied, symmetrize
        do k = 1, tb%nshells
            tb%tens(k)%V(:) = matmul(tb%tens(k)%relations,tb%tens(k)%V(:))
        enddo
        ! things are looking good up to this point.
        !  ... tb%tens(k)%V(:) =
        !              1.00000        0.00000        0.00000        0.00000
        !              0.00000        2.00000        0.00000        0.00000
        !              0.00000        0.00000        2.00000        0.00000
        !              0.00000        0.00000        0.00000        2.00000
        !  ... tb%tens(k)%V(:) =
        !              3.00000        4.00000        4.00000        4.00000
        !             -4.00000        6.00000        5.00000        5.00000
        !             -4.00000        5.00000        6.00000        5.00000
        !             -4.00000        5.00000        5.00000        6.00000
    end subroutine set_matrix_element

    function       get_matrix_element(tb,ip_id) result(V)
        ! irreducible pair (ip_id), can be negative (corresponds to pair n-m rathet than m-n, on the
        ! opposite [upper/lower] side of the Hamiltonian), in which case adjoint of V is returned
        implicit none
        !
        class(am_class_tightbinding), intent(in) :: tb
        integer , intent(in) :: ip_id
        real(dp), allocatable :: V(:,:)
        integer :: id
        integer :: m,n
        !
        ! note the absolute value
        id = abs(ip_id)
        ! get primitive atom indices
        m = tb%tens(id)%dims(1)
        n = tb%tens(id)%dims(2)
        ! if irreducible pair id is negative, it means the pair was flipped
        if (ip_id.lt.0) then
            allocate(V, source=adjoint(reshape(tb%tens(id)%V,[m,n])) )
        else
            allocate(V, source=        reshape(tb%tens(id)%V,[m,n])  )
        endif
        !
    end function   get_matrix_element

    subroutine     write_irreducible_matrix_element(tb)
        !
        implicit none
        !
        class(am_class_tightbinding), intent(in) :: tb
        integer :: k ! shell index
        integer :: alpha
        integer :: beta
        integer :: i
        integer :: fid
        ! create tb dir
        call execute_command_line ('mkdir -p '//trim(outfile_dir_tb))
        ! write point group
        fid = 1
        open(unit=fid,file=trim(outfile_dir_tb)//'/'//'outfile.tb_matrix_elements_irreducible',status='replace',action='write')
            !
            write(fid,'(a,a16,3a6)') '#', 'V', 'shell', 'alpha', 'beta'
            do i = 1, tb%nVs
                ! get shell index
                k     = tb%V_ind(1,i)
                ! get orbital indices
                alpha = tb%V_ind(2,i)
                beta  = tb%V_ind(3,i)
                ! write 
                write(fid,'(i5,f16.8,3i6)') i, tb%V(i), k, alpha, beta
            enddo
            !
        close(fid)
    end subroutine write_irreducible_matrix_element

    subroutine     read_irreducible_matrix_element(tb)
        !
        implicit none
        !
        class(am_class_tightbinding), intent(inout) :: tb
        character(maximum_buffer_size) :: buffer ! read buffer
        character(len=:), allocatable :: word(:) ! read buffer
        integer :: fid
        real(dp), allocatable :: V(:)
        integer :: k
        integer :: alpha
        integer :: beta
        integer :: i
        integer :: j
        !
        allocate(V(tb%nVs))
        !
        fid = 1
        open(unit=fid,file='infile.tb_matrix_elements_irreducible',status="old",action='read')
            ! skip header
            read(fid,*)
            ! read matrix elements
            do j = 1, tb%nVs
                read(unit=fid,fmt='(a)') buffer
                word = strsplit(buffer,delimiter=' ')
                !
                read(word(1),*) i
                read(word(2),*) V(i)
                read(word(3),*) k
                read(word(4),*) alpha
                read(word(5),*) beta
                !
                if (.not.isequal(    k,tb%V_ind(1,i))) stop 'ERROR [read_irreducible_matrix_element]: k /= V_ind(1,i)'
                if (.not.isequal(alpha,tb%V_ind(2,i))) stop 'ERROR [read_irreducible_matrix_element]: alpha /= V_ind(2,i)'
                if (.not.isequal( beta,tb%V_ind(3,i))) stop 'ERROR [read_irreducible_matrix_element]: beta /= V_ind(3,i)'
            enddo
            call tb%set_matrix_element(V)
        close(fid)
    end subroutine read_irreducible_matrix_element

    ! DEPRECATED: 
    ! TODO, convert write_all_matrix_element into an extenral subroutine for getting the hamiltonian which takes as input a vector of irreducible matrix elements ... 
    ! make it such that the subroutine can be compiled afterwards...
    ! if this is done, the model is set. do not need to go through everything again.

    subroutine     write_all_matrix_element(tb,pp)
        !
        implicit none
        !
        class(am_class_tightbinding), intent(in) :: tb    ! irreducible tight binding matrix elements
        type(am_class_prim_pair)    , intent(in) :: pp    ! primitive pairs
        real(dp), allocatable, target :: Hsub_target(:,:)
        real(dp), allocatable, target :: pg_target(:,:,:)
        integer , allocatable :: S(:) ,E(:)
        real(dp), pointer :: Dm(:,:), Dn(:,:)
        real(dp), pointer :: Hsub(:,:)
        integer :: fid
        integer :: k,p, m,n
        !
        ! allocate space for vectors demarking start and end of Hamiltonian subsection
        allocate(Hsub_target(tb%pg%nbases,tb%pg%nbases))
        allocate(pg_target, source=tb%pg%sym)
        allocate(S, source=tb%pg%S)
        allocate(E, source=tb%pg%E)
        !
        call tostring_set(sep=' ')
            ! create directory
            call execute_command_line ('mkdir -p '//trim(outfile_dir_tb))
            ! export symmetry
            fid = 1
            open(unit=fid,file=trim(outfile_dir_tb)//'/'//'outfile.tb_matrix_elements',status='replace',action='write')
                ! print stuff
                write(fid,'(a,a)')                 tostring(pp%nshells), ' primitive shells'
                write(fid,'(a,a)')                          tostring(S), ' subregion start'
                write(fid,'(a,a)')                          tostring(E), ' subregion end'
                do k = 1, pp%nshells
                ! abbreviations
                m = pp%shell(k)%m
                n = pp%shell(k)%n
                ! print stuff
                write(fid,'(a,a)')                       repeat('=',79), repeat('=',79)
                write(fid,'(a,a)')              centertitle(tostring(k)//' shell',158)
                write(fid,'(a,a)')                          tostring(m), ' primitive m'
                write(fid,'(a,a)')                          tostring(n), ' primitive n'
                write(fid,'(a,a)')         tostring(pp%shell(k)%natoms), ' atoms in shell'
                do p = 1, pp%shell(k)%natoms
                ! abbreviations
                Dm   => pg_target(S(m):E(m), S(m):E(m), pp%shell(k)%pg_id(p) )
                Dn   => pg_target(S(n):E(n), S(n):E(n), pp%shell(k)%pg_id(p) )
                Hsub => Hsub_target(S(m):E(m), S(n):E(n))
                ! print stuff
                Hsub = get_matrix_element(tb=tb, ip_id=pp%ip_id(k)) 
                ! print
                call disp(unit=fid, fmt='f18.10',x=matmul(matmul(transpose(Dm), Hsub), Dn),title=tostring(pp%shell(k)%tau_cart(1:3,p),fmt='f10.5')//' [cart.]',style='above')
                enddo
                enddo
            close(fid)
        ! return delimiter to ","
        call tostring_set(sep=',')
        !
    end subroutine write_all_matrix_element

    ! optimizer

    subroutine     optimize_matrix_element(tb,dft,pp,opts,flags)
        !
        implicit none
        !   
        class(am_class_tightbinding), intent(inout) :: tb
        type(am_class_dft)          , intent(in) :: dft
        type(am_class_prim_pair)    , intent(in) :: pp
        type(am_class_options)      , intent(in) :: opts
        character(*)                , intent(in) :: flags
        type(am_class_options) :: notalk
        real(dp) :: best_rms
        real(dp), allocatable :: best_V(:)
        integer, allocatable :: selector_kpoint(:)
        integer :: i,j,k
        !
        if (opts%verbosity.ge.1) call print_title('Optimized matrix elements')
        ! supress output
        notalk = opts
        notalk%verbosity = 0
        !
        if      (index(flags,'input').ne.0) then
            ! optimize matrix elements starting from input matrix elements
            ! [1:size(tb%tens)] ! all shells
            ! [1:dft%bz%nkpts]  ! all kpoints
            call tb%initialize_ft(selector_shell=[1:tb%nshells],selector_kpoint=[1:dft%bz%nkpts],opts=opts)
            ! 
            call perform_optimization(tb=tb, dft=dft, pp=pp, opts=opts)
            !
        elseif  (index(flags, 'shell_progressive').ne.0) then
            ! course optimization: use the same number of kpoints as there are matrix elements in the total model 
            ! progressively add on one shell at a time to the model 
            ! each time, search optimal parameters starting at 10 stochastically choosen staring conditions
            ! keep the best one. the goal here is to deal with cases where bands are crossing when they shouldn't and vice versa
            ! fine optimization: perform 1 minimization using all kpoints. use as starting condition the best value of the course optimization
            ! ----------------------------------------------------
            ! intialize rms at a large number
            best_rms = 1.0D15
            ! intialize best matrix elements
            allocate(best_V(tb%nVs))
            best_V = 0
            ! check for enough kpoints
            write(*,*)
            if (dft%bz%nkpts.lt.tb%nVs) stop 'ERROR [optimize_matrix_element]: not enough kpoints'
            ! set kpoint selector
            allocate(selector_kpoint(tb%nVs))
            selector_kpoint(1:tb%nVs) = floor(linspace(1, dft%bz%nkpts, tb%nVs ))
            ! loop to find best (semi-stochastic + NNLS)
            do k = 1, tb%nshells
                ! initialize fitter
                call tb%initialize_ft(selector_shell=[1:k],selector_kpoint=selector_kpoint,opts=notalk)
                ! print stdout
                if (k.eq.1) then
                    write(*,'(5x,a8,a10,a10,a)') 'iter', 'rms ', 'log(d rms)',centertitle('parameters',tb%nVs*10)
                endif
                ! stochastic part
                do i = 1, 10
                    if (i.eq.1) then
                        ! write stdout header
                        if (opts%verbosity.ge.1) then
                            write(*,'(5x,a8,a10,a10,a,a)') 'shell '//tostring(k), ' '//repeat('-',10-1), ' '//repeat('-',10-1), ' '//repeat('-',tb%ft%nxs*10-1)
                        endif
                    else
                        ! stochastically determine a new starting condition based on best value
                        do j = 1, tb%nVs
                            if ( any(tb%ft%selector_x.eq.j) ) then
                                tb%V(j) = tb%V(j) + 5.0_dp*(rand() - 0.50_dp)*best_rms
                                ! tb%V(j) = tb%V(j) + (rand() - 0.50_dp)*abs(tb%V(j))
                            endif
                        enddo
                    endif
                    ! perform optimization utilizing jacobian
                    call perform_optimization(tb=tb, dft=dft, pp=pp, opts=notalk)
                    ! check rms
                    if (tb%ft%rms.lt.best_rms) then
                        ! if the current rms is better, update best value
                        best_rms = tb%ft%rms
                        best_V   = tb%V
                    else
                        ! otherwise reset to best value
                        tb%V = best_V
                    endif
                    ! write stdout
                    if (opts%verbosity.ge.1) write(*,'(5x,i8,f10.2, 10x ,SP,1000f10.2)') i, tb%ft%rms, tb%ft%x
                enddo
            enddo
            ! initialize final optimization using best matrix elements found thus far
            tb%V = best_V
            ! initialize final optimization using all kpoints
            call tb%initialize_ft(selector_shell=[1:tb%nshells],selector_kpoint=[1:dft%bz%nkpts],opts=opts)
            ! perform final optimization
            call perform_optimization(tb=tb, dft=dft, pp=pp, opts=opts)
            ! output results
            write(*,'(a,a)') flare, 'final rms = '//tostring(tb%ft%rms)
            !
        endif
        !
        if (debug) then
            write(*,'(a,a)')   flare, 'optimization parameters and conditions:'
            write(*,'(5x,a,a)') 'irreducible matrix elements = '//tostring(tb%ft%nxs)
            write(*,'(5x,a,a)') 'selected k-points = '          //tostring(tb%ft%nkpts)
            write(*,'(5x,a,a)') 'bands = '                      //tostring(tb%ft%nbands)
            write(*,'(5x,a,a)') 'bands skipped = '              //tostring(tb%ft%skip_band)
            write(*,'(5x,a,a)') 'residual vector length = '     //tostring(tb%ft%nrs)
            write(*,'(5x,a,a)') 'max iterations = '             //tostring(tb%ft%maxiter)
        endif
    end subroutine optimize_matrix_element

    subroutine     initialize_ft(tb,selector_shell,selector_kpoint,opts)
        !
        implicit none
        !
        class(am_class_tightbinding), intent(inout) :: tb
        integer                     , intent(in) :: selector_shell(:)
        integer                     , intent(in) :: selector_kpoint(:)
        type(am_class_options)      , intent(in) :: opts
        integer :: i
        ! defaults
        if (allocated(tb%ft%selector_shell)) deallocate(tb%ft%selector_shell)
        allocate(tb%ft%selector_shell, source=selector_shell)
        ! defaults
        if (allocated(tb%ft%selector_kpoint)) deallocate(tb%ft%selector_kpoint)
        allocate(tb%ft%selector_kpoint,source=selector_kpoint)
        ! set number of bands to fit
        tb%ft%nbands = maxval(tb%pg%E)
        ! set parameter selector based on irreducible shells
        if (allocated(tb%ft%selector_x)) deallocate(tb%ft%selector_x)
        allocate(tb%ft%selector_x(tb%nVs))
        do i = 1, tb%nVs
            if ( any(tb%ft%selector_shell.eq.tb%V_ind(1,i)) ) then
                tb%ft%selector_x(i) = i
            else
                tb%ft%selector_x(i) = 0
            endif
        enddo
        tb%ft%selector_x = trim_null(tb%ft%selector_x)
        ! get number of fitting parameters
        tb%ft%nxs = count(tb%ft%selector_x.ne.0)
        ! set number of iterations
        tb%ft%maxiter = 10
        ! initialize parameter vector
        if (allocated(tb%ft%x)) deallocate(tb%ft%x)
        allocate(tb%ft%x(tb%ft%nxs))
        tb%ft%x     = 0             
        ! get number of kpoints to fit
        tb%ft%nkpts = size(tb%ft%selector_kpoint)
        ! get size of residual vector
        tb%ft%nrs   = tb%ft%nkpts * tb%ft%nbands
        ! initialize residual vector
        if (allocated(tb%ft%r)) deallocate(tb%ft%r)
        allocate(tb%ft%r(tb%ft%nrs))
        tb%ft%r     = 0             
        ! initialize rms error
        tb%ft%rms   = 0              ! rms error
        ! set number of bands to skip
        tb%ft%skip_band = opts%skip_band
    end subroutine initialize_ft

    subroutine     perform_optimization(tb,dft,pp,opts)
        !
        use mkl_rci
        use mkl_rci_type
        !
        implicit none
        !
        type(am_class_tightbinding)  , intent(inout) :: tb
        type(am_class_dft)           , intent(in) :: dft
        type(am_class_prim_pair)     , intent(in) :: pp
        type(am_class_options)       , intent(in) :: opts
        real(dp), allocatable :: x(:)
        real(dp), allocatable :: FVEC(:)
        real(dp), allocatable :: FJAC(:,:)
        type(handle_tr) :: handle
        real(dp) :: eps(6)
        integer  :: info(6)
        integer  :: rci_request
        integer  :: successful
        integer  :: i
        ! precisions for stop-criteria
        eps(1:6) = 1.0D-6
        ! initialize fitting parameters (irreducible matrix elements) 
        allocate(x(tb%ft%nxs))
        x = tb%V 
        ! initialize residual vector
        allocate(FVEC(tb%ft%nrs))
        FVEC = 0.0_dp
        ! initialize jacobian matrix
        allocate(FJAC(tb%ft%nrs,tb%ft%nxs))
        FJAC = 0
        ! initializes the solver of a nonlinear least squares problem
        if (dtrnlsp_init(handle=handle, n=tb%ft%nxs, m=tb%ft%nrs, x=x, eps=eps, iter1=tb%ft%maxiter, iter2=100, rs=100.0D0) /= tr_success) then
            call mkl_free_buffers
            stop 'ERROR [perform_optimization]: dtrnlsp_init'
        end if
        ! check the correctness of handle and arrays containing Jacobian matrix, objective function, and stopping criteria
        if (dtrnlsp_check(handle=handle, n=tb%ft%nxs, m=tb%ft%nrs, FJAC=FJAC, FVEC=FVEC, eps=eps, info=info) /= tr_success) then
            call mkl_free_buffers
            stop 'ERROR [perform_optimization]: dtrnlspbc_init'
        else
            if ( info(1) /= 0 .or. info(2) /= 0 .or. info(3) /= 0 .or. info(4) /= 0 ) then
                call mkl_free_buffers
                stop 'ERROR [perform_optimization]: dtrnlspbc_init, invalid input parameters'
            endif
        endif
        ! set initial rci variables
        rci_request = 0
        successful = 0
        i = 0
        ! write header
        if (opts%verbosity.ge.1) write(*,'(5x,a8,a10,a10,a)') 'iter', 'rms ', 'log(d rms)',centertitle('parameters',tb%ft%nxs*10)
        if (opts%verbosity.ge.1) write(*,'(5x,a8,a10,a10,a)') ' '//repeat('-',8-1), ' '//repeat('-',10-1), ' '//repeat('-',10-1), ' '//repeat('-',tb%ft%nxs*10-1)
        ! enter optimization loop
        do while (successful == 0)
            ! solve nonlinear least squares problem using the TR algorithm
            if (dtrnlsp_solve(handle=handle, FVEC=FVEC, FJAC=FJAC, rci_request=rci_request) /= tr_success) then
                call mkl_free_buffers
                stop 'ERROR [perform_optimization]: dtrnlsp_solve'
            endif
            select case (rci_request)
            case (-1, -2, -3, -4, -5, -6)
                successful = 1
            case (1)
                ! update x in tb model
                call tb%set_matrix_element(V=x,flags='selector,decompress')
                ! recalculate function
                FVEC = compute_residual(tb=tb, dft=dft, pp=pp)
                ! save result
                tb%ft%r   = FVEC 
                tb%ft%rms = sqrt(sum(FVEC**2)/tb%ft%nrs)
                ! increase counter
                i = i + 1
                ! print results
                if (opts%verbosity.ge.1) then
                    if (i.eq.1) then
                        write(*,'(5x,i8,f10.2, 10x ,SP,1000f10.2)') i, tb%ft%rms, x
                    else
                        write(*,'(5x,i8,f10.2,g10.2,SP,1000f10.2)') i, tb%ft%rms, abs(norm(1.0D-14 + x-tb%ft%x )), x
                    endif
                endif
                ! save results
                tb%ft%x = x
            case (2)
                ! update x in tb model
                call tb%set_matrix_element(V=x,flags='selector,decompress')
                ! compute jacobian matrix (uses central difference)
                FJAC = compute_jacobian(tb=tb, dft=dft, pp=pp)
            end select
        end do
        ! clean up
        if (dtrnlsp_delete(handle) /= tr_success) then
            call mkl_free_buffers
            stop 'ERROR [perform_optimization]: dtrnlsp_delete'
        else
            call mkl_free_buffers
        endif
    end subroutine perform_optimization

    function       compute_residual(tb,dft,pp) result(R)
        !
        implicit none
        !
        type(am_class_tightbinding), intent(inout) :: tb
        type(am_class_prim_pair)   , intent(in) :: pp
        type(am_class_dft)         , intent(in) :: dft
        real(dp), allocatable :: R(:)
        ! get tb dispersion 
        call tb%get_dispersion(pp=pp, bz=dft%bz)
        ! calculate residual vector indices
        R = pack(dft%dr%E([1:tb%ft%nbands]+tb%ft%skip_band, tb%ft%selector_kpoint) - tb%dr%E(:,tb%ft%selector_kpoint) , .true.)
        !
    end function   compute_residual

    function       compute_jacobian(tb,dft,pp) result(FJAC)
        !
        use mkl_rci
        !
        implicit none 
        !
        type(am_class_tightbinding), intent(inout) :: tb
        type(am_class_prim_pair)   , intent(in) :: pp
        type(am_class_dft)         , intent(in) :: dft
        real(dp), allocatable :: FJAC(:,:) ! fjac(m,n) jacobian matrix
        real(dp), allocatable :: f1(:)     ! f1(m)     residual vector
        real(dp), allocatable :: f2(:)     ! f1(m)     residual vector
        real(dp), allocatable :: x(:)      ! x(n)      parameter vector
        real(dp) :: eps_jac
        integer*8:: handle 
        integer  :: successful
        integer  :: rci_request
        ! set epsilon for jacobian calc
        eps_jac = 1.d-5
        ! allocate space
        allocate(  f1(tb%ft%nrs))
        allocate(  f2(tb%ft%nrs))
        allocate(FJAC(tb%ft%nrs,tb%ft%nxs))
        ! initialize x
        allocate(x, source=tb%V)
        ! begin jacobian computation
        if (djacobi_init(handle=handle, n=tb%ft%nxs, m=tb%ft%nrs, x=x, FJAC=FJAC, eps=eps_jac) .ne. tr_success) then
            print *, '#fail: error in djacobi_init' 
            call mkl_free_buffers
            stop 1
        end if
        rci_request = 0
        successful  = 0
        do while (successful.eq.0)
            if (djacobi_solve(handle=handle, f1=f1, f2=f2, rci_request=rci_request) .ne. tr_success) then
                print *, '#fail: error in djacobi_solve'
                call mkl_free_buffers
                stop 1
            end if
            if (rci_request .eq. 1) then
                ! update x in TB model
                call tb%set_matrix_element(V=x,flags='selector,decompress')
                ! compute jacobian
                f1 = compute_residual(tb=tb,dft=dft,pp=pp)
            else if (rci_request .eq. 2) then
                ! update x in TB model
                call tb%set_matrix_element(V=x,flags='selector,decompress')
                ! compute jacobian
                f2 = compute_residual(tb=tb,dft=dft,pp=pp)
            else if (rci_request .eq. 0) then
                successful = 1 
            end if
        end do
        if (djacobi_delete(handle) .ne. tr_success) then
            print *, '#fail: error in djacobi_delete' 
            call mkl_free_buffers
            stop 1
        else
            call mkl_free_buffers
        end if
    end function   compute_jacobian

end module am_tight_binding











