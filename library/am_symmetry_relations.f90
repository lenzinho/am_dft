module am_symmetry_relations

    use am_mkl
    use am_matlab
    use am_stdout
    use dispmodule

	implicit none

	public

	contains

    function       combine_relations(relationsA,relationsB,relationsC) result(relations)
        !
        ! combins up to three relations
        !
        implicit none
        !
        real(dp), intent(in) :: relationsA(:,:)
        real(dp), intent(in) :: relationsB(:,:)
        real(dp), intent(in), optional :: relationsC(:,:)
        real(dp), allocatable :: relations(:,:)     !
        real(dp), allocatable :: A(:,:)             ! A(2*flat%nbases,2*flat%nbases) - augmented matrix equation
        real(dp), allocatable :: LHS(:,:)           ! LHS(flat%nbases,flat%nbases)   - left hand side of augmented matrix equation (should be identity after reducing to row echlon form)
        real(dp), allocatable :: RHS(:,:)           ! RHS(flat%nbases,flat%nbases)   - right hand side of augmented matrix equation
        integer , allocatable :: indices(:)         ! used for clarity
        integer :: nbases
        !
        !
        if (size(relationsA,1).ne.size(relationsA,2)) stop 'Dimension mismatch: A1 vs A2'
        if (size(relationsA,1).ne.size(relationsB,1)) stop 'Dimension mismatch: A1 vs B1'
        if (size(relationsA,1).ne.size(relationsB,2)) stop 'Dimension mismatch: A1 vs B2'
        if (present(relationsC)) then
        if (size(relationsA,1).ne.size(relationsC,1)) stop 'Dimension mismatch: A1 vs C1'
        if (size(relationsA,1).ne.size(relationsC,2)) stop 'Dimension mismatch: A1 vs C2'
        endif
        !
        nbases = size(relationsA,1)
        !
        ! get indices
        allocate(indices,source=[1:nbases])
        ! initialize augmented workspace matrix A
        allocate(A(3*nbases,2*nbases))
        A = 0
        ! construct slice of A
        A(0*nbases+indices,1*nbases+indices) = relationsA
        A(0*nbases+indices,0*nbases+indices) = eye(nbases)
        A(1*nbases+indices,1*nbases+indices) = relationsB
        A(1*nbases+indices,0*nbases+indices) = eye(nbases)
        if (present(relationsC)) then
        A(2*nbases+indices,1*nbases+indices) = relationsC
        A(2*nbases+indices,0*nbases+indices) = eye(nbases)
        endif
        ! incorporate symmetry via lu factorization (equivalent to applying rref)
        call lu(A)
        ! Apply Gram-Schmidt orthogonalization to obtain A in reduced row echelon form
        call rref(A)
        ! correct basic rounding error
        where (abs(nint(A)-A).lt.tiny) A = nint(A)
        ! At this point, A = [ LHS | RHS ], in which LHS = E, identity matrix; A completely specifies all relationships between variables: LHS = RHS.
        allocate(LHS(nbases,nbases))
        allocate(RHS(nbases,nbases))
        LHS = A(0*nbases+indices,0*nbases+indices)
        RHS = A(0*nbases+indices,1*nbases+indices)
        !
        ! checks
        if (.not.isequal(LHS,eye(nbases))) then
            stop 'Failed to reduce matrix to row echlon form.'
        endif
        if (count(get_null(RHS))+count(get_independent(RHS))+count(get_depenent(RHS)).ne.nbases) then
            stop 'Number of null, independent, and dependent terms do not sum to the number of terms.'
        endif
        !
        allocate(relations, source=RHS)
        !
    end function   combine_relations

    pure function  get_null(relations) result(is_null)
        !
        implicit none
        !
        real(dp), intent(in) :: relations(:,:)
        logical , allocatable :: is_null(:)
        !
        ! null terms (equal zero)
        allocate(is_null, source=all(abs(relations).lt.tiny,2))
        !
    end function   get_null

    pure function  get_independent(relations) result(is_independent)
        !
        implicit none
        !
        real(dp), intent(in) :: relations(:,:)
        logical , allocatable :: is_independent(:)
        !
        ! independent terms (equal themselves and nothing else)
        allocate(is_independent, source=all(abs(relations-eye(size(relations,1))).lt.tiny,2))
        !
    end function   get_independent

    pure function  get_depenent(relations) result(is_dependent)
        !
        implicit none
        !
        real(dp), intent(in) :: relations(:,:)
        logical , allocatable :: is_dependent(:)
        !
        ! dependent terms (can be written via independent terms)
        allocate(is_dependent, source=any(abs(relations).gt.tiny,2))
        is_dependent = (is_dependent.and..not.get_independent(relations))
        !
    end function   get_depenent

    function       get_basis_label(dims,ind,flags) result(str)
        !
        implicit none
        !
        integer     , intent(in) :: dims(:)
        integer     , intent(in) :: ind
        character(*), intent(in) :: flags
        character(:), allocatable :: str
        integer     , allocatable :: sub(:)
        integer :: i, n
        !
        allocate(character(100) :: str)
        !
        n = size(dims)
        !
        sub = ind2sub(dims=dims,ind=ind)
        !
        i = 1
        if     (index(flags,'bra').ne.0) then
            str = '<'//trim(int2char(sub(i)))
        elseif (index(flags,'ket').ne.0) then
            str = '|'//trim(int2char(sub(i)))
        else
            str = 'a('//trim(int2char(sub(i)))
        endif
        !
        if (n.ge.2) then
        do i = 2, n
            str = trim(str)//','//trim(int2char(sub(i)))
        enddo
        endif
        !
        if     (index(flags,'bra').ne.0) then
            str = trim(str)//'|'
        elseif (index(flags,'ket').ne.0) then
            str = trim(str)//'>'
        else
            str = trim(str)//')'
        endif
        !
    end function   get_basis_label

    subroutine     print_relations(relations,dims,flags)
        !
        ! flags - null, dependent, independent, header
        !       - bra, ket
        !
        implicit none
        !
        real(dp), intent(in) :: relations(:,:)
        integer , optional, intent(in) :: dims(:)
        character(*), intent(in) :: flags
        integer :: i, j 
        logical, allocatable :: is_null(:)
        logical, allocatable :: is_independent(:)
        logical, allocatable :: is_dependent(:)
        character(:), allocatable :: str
        integer :: nterms
        !
        ! get terms
        nterms = size(relations,1)
        is_null = get_null(relations)
        is_independent = get_independent(relations)
        is_dependent = get_depenent(relations)
        !
        ! print header
        if (index(flags,'header').ne.0) then
            write(*,'(a5,a,a)',advance='no') ' ... ', trim(int2char(nterms)), ' terms = '
            write(*,'(i4,a,f5.1,a)',advance='no') count(is_null)       , ' null ('       , count(is_null)       /real(nterms,dp)*100.0_dp , '%) '
            write(*,'(i4,a,f5.1,a)',advance='no') count(is_dependent)  , ' dependent ('  , count(is_dependent)  /real(nterms,dp)*100.0_dp , '%) '
            write(*,'(i4,a,f5.1,a)',advance='no') count(is_independent), ' independent (', count(is_independent)/real(nterms,dp)*100.0_dp , '%) '
            writE(*,*)
        endif
        !
        ! print irreducible symmetry relations
        if ((index(flags,'independent').ne.0).or.(index(flags,'dependent').ne.0).or.(index(flags,'null').ne.0)) then
            !
            if (.not.present(dims)) stop 'dims is required for printing relations'
            !
            ! write the independent terms (equal only to themselves)
            if (index(flags,'independent').ne.0) then
            do i = 1, nterms
                if (is_independent(i)) then
                    str = get_basis_label(dims=dims, ind=i ,flags=flags)
                    write(*,'(5x,a,a,a)',advance='no') trim(str),' = ', trim(str)
                    write(*,*)
                endif
            enddo
            endif
            !
            ! write the dependent terms
            if (index(flags,'dependent').ne.0) then
            do i = 1,nterms
                if (is_dependent(i)) then
                    str = get_basis_label(dims=dims, ind=i ,flags=flags)
                    write(*,'(5x,a,a)',advance='no') trim(str), ' = '
                    do j = 1,nterms
                        if (abs(relations(i,j)).gt.tiny) then
                            str = get_basis_label(dims=dims, ind=j ,flags=flags)
                            write(*,'(a,a,a)',advance='no') trim(dbl2charSP(relations(i,j),7)), '*', trim(str)
                        endif
                    enddo
                    write(*,*)
                endif
            enddo
            endif
            !
            ! write null terms
            if (index(flags,'null').ne.0) then
            do i = 1, nterms
                if (is_null(i)) then
                    str = get_basis_label(dims=dims, ind=i ,flags=flags)
                    write(*,'(5x,a,a)',advance='no') trim(str),' = 0'
                    write(*,*)
                endif
            enddo
            endif
        endif
        !
    end subroutine print_relations

    subroutine     dump_relations(relations,iopt_filename)
        !
        implicit none
        !
        real(dp), intent(in) :: relations(:,:)
        character(*), intent(in), optional :: iopt_filename
        character(100) :: fname
        integer :: fid
        integer :: m,n
        integer :: j,k
        !
        ! set default
        fname = 'dump.relations'
        if (present(iopt_filename)) fname = iopt_filename
        !
        m = size(relations,1)
        n = size(relations,2)
        !
        fid = 1
        open(unit=fid,file=trim(fname),status="replace",action='write')
            !
            do j = 1, m
            do k = 1, n
                write(fid,"(f)",advance='no') relations(j,k)
            enddo
            write(fid,*)
            enddo
            !
        close(fid)
        !
    end subroutine dump_relations

    subroutine     export_relations2matlab(relations,dims,fnc_name,fid_append,flags)
        !
        ! creates .m file based on relations
        ! flags = append/standalone
        implicit none
        !
        real(dp),           intent(in) :: relations(:,:)
        integer , optional, intent(in) :: dims(:)
        character(*),       intent(in) :: fnc_name
        character(*),       intent(in) :: flags
        integer, optional , intent(in) :: fid_append
        integer :: fid
        integer :: i, j 
        logical, allocatable :: is_null(:)
        logical, allocatable :: is_independent(:)
        logical, allocatable :: is_dependent(:)
        character(:), allocatable :: str
        integer :: nterms
        !
        ! get terms
        nterms         = size(relations,1)
        is_null        = get_null(relations)
        is_independent = get_independent(relations)
        is_dependent   = get_depenent(relations)
        !
        if     (index(flags,'append').ne.0) then
            if (.not.present(fid_append)) stop 'ERROR [export_relations2matlab]: fid_append required.'
            fid = (fid_append)
        elseif (index(flags,'standalone').ne.0) then
            fid = 1
            open(unit=fid,file=trim(fnc_name)//'.m',status='replace',action='write')
        else
            stop 'ERROR [export_relations2matlab]: flags != standalone/append.'
        endif
            !
            write(fid,'(a,a,a)') 'function [a] = ', trim(fnc_name), '(v)'
                ! write the independent terms
                j = 0
                do i = 1, nterms
                    if (is_independent(i)) then
                        j=j+1
                        str = get_basis_label(dims=dims, ind=i ,flags='')
                        write(fid,'(a)',advance='no') trim(str)//' = v('//tostring(j)//');'
                        write(fid,*)
                    endif
                enddo
                ! write the dependent terms
                do i = 1,nterms
                    if (is_dependent(i)) then
                        str = get_basis_label(dims=dims, ind=i ,flags='')
                        write(fid,'(a)',advance='no') trim(str)//' = '
                        do j = 1,nterms
                            if (abs(relations(i,j)).gt.tiny) then
                                str = get_basis_label(dims=dims, ind=j ,flags='')
                                if (present(fid_append)) then
                                    ! HIGH PRECISION
                                    write(fid,'(a)',advance='no') tostring(relations(i,j),fmt='SP,f24.16')//'*'//trim(str)//';'
                                else
                                    write(fid,'(a)',advance='no') tostring(relations(i,j),fmt='SP,f10.5')//'*'//trim(str)//';'
                                endif
                            endif
                        enddo
                        write(fid,*)
                    endif
                enddo
                ! write null terms
                do i = 1, nterms
                    if (is_null(i)) then
                        str = get_basis_label(dims=dims, ind=i ,flags='')
                        write(fid,'(a)',advance='no') trim(str)//' = 0;'
                        write(fid,*)
                    endif
                enddo
                !
            write(fid,'(a)') 'end'
            !
        ! only close file if standalone
        if (index(flags,'standalone').ne.0) then
        close(fid)
        endif
    end subroutine export_relations2matlab

end module am_symmetry_relations











