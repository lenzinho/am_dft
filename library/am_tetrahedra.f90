module am_tetrahedra

	use am_brillouin_zone
    use am_prim_cell
    use am_symmetry
    use am_constants
    use am_stdout
    use am_vasp_io
    use am_mkl
    use am_options

	implicit none

	private

    type, public :: am_class_tetrahedra
        integer :: ntets !> ntet number f tetrahedra
        real(dp), allocatable :: volume(:) !> volume(ntet) tetrahedra volume
        integer , allocatable :: tet(:,:) !> tet(4,ntet) indices of kpoints at the tet of the tetrahedra
        real(dp), allocatable :: w(:) !> w(ntet) tetrahedra weights
    contains
        procedure :: load_ibzkpt => load_ibzkpt
    end type am_class_tetrahedra

contains

    ! functions which operate on tet or similar derived types
    
    subroutine     load_ibzkpt(tet,opts)
        ! 
        ! Reads the IBZKPT file.
        !
        implicit none
        !
        class(am_class_tetrahedra), intent(inout) :: tet
        type(am_class_options), intent(in) :: opts
        !
        call read_ibzkpt(ntets = tet%ntets,&
                         vtet  = tet%volume,&
                         tet   = tet%tet,&
                         wtet  = tet%w,&
            	 iopt_filename = opts%ibzkpt,&
            	iopt_verbosity = opts%verbosity )
        !
    end subroutine load_ibzkpt

    subroutine     tetrahedra_pdos(bz,tet,opts)
        !
        use am_histogram
        use am_options
        !
        implicit none
        !
        class(am_class_bz), intent(inout) :: bz
        type(am_class_tetrahedra), intent(in) :: tet
        type(am_class_options), intent(in) :: opts
        real(dp) :: maxE ! maximum band energy
        real(dp) :: minE ! minimum band energy
        real(dp) :: dE   ! energy increment
        real(dp) :: w(4) ! weights on tetrahedron tet
        real(dp) :: Ec(4) ! tet energies
        integer :: i, j, k, m, n, o
        !
        if (opts%verbosity.ge.1) call am_print_title('Tetrahedron integration')
        !
        dE = 0.01_dp
        minE = minval(bz%E(:,:))
        maxE = maxval(bz%E(:,:))
        ! bz%Ep = linspace(minE,maxE,100)
        bz%Ep = regspace(minE,maxE,dE)
        bz%nEp = size(bz%Ep,1)
        !
        call plot_histogram( histogram(pack(bz%E(:,:),(bz%E(:,:).ge.minE).and.(bz%E(:,:).le.maxE)),column_width) )
        !
        if (opts%verbosity.ge.1) call am_print('number of bands',bz%nbands,' ... ')
        if (opts%verbosity.ge.1) call am_print('number of tetrahedra',tet%ntets,' ... ')
        if (opts%verbosity.ge.1) call am_print('lowest band energy',minE,' ... ')
        if (opts%verbosity.ge.1) call am_print('highest band energy',maxE,' ... ')
        if (opts%verbosity.ge.1) call am_print('minimum probing energy',bz%Ep(1),' ... ')
        if (opts%verbosity.ge.1) call am_print('maximum probing energy',bz%Ep(bz%nEp),' ... ')
        if (opts%verbosity.ge.1) call am_print('energy increments dE',dE,' ... ')
        if (opts%verbosity.ge.1) call am_print('number of probing energies',bz%nEp,' ... ')
        !
        allocate(bz%pdos(bz%nspins,bz%norbitals,bz%nions,bz%nEp))
        allocate(bz%dos(bz%nEp))
        !
        bz%pdos = 0.0_dp
        bz%dos = 0.0_dp
        !
        ! two things
        ! 1) difference between total dos and sum of JDOS. ONE band is being left out in the JDOS summation. figure out which one.
        !  0.996425152311378        1.00000000031254
        ! 2) delta summation produces wild values. the divisions by 1/tiny are numerically unstable. Figure out how to get rid of them.!
        !
        do i = 1, bz%nEp
            do j = 1, bz%nbands
            do k = 1, tet%ntets
                !
                Ec(1:4) = bz%E(j,tet%tet(:,k))
                !
                w = tetrahedron_weights(x=bz%Ep(i),xc=Ec,integration_type='delta',apply_blochl_corrections=.true.)
                !
                bz%dos(i) = bz%dos(i) + tet%w(k)*tet%volume(k)*sum(w)
                !
                ! <projections>
                !
                do m = 1,bz%nspins
                do n = 1,bz%norbitals
                do o = 1,bz%nions
                    !> pdos(nspins,norbitals,nions,nEP)
                    !> lmproj(nspins,norbitals,nions,nbands,nkpts)
                    bz%pdos(m,n,o,i) = bz%pdos(m,n,o,i) + tet%w(k)*tet%volume(k)*sum(w*bz%lmproj( m,n,o,j,tet%tet(:,k) ))/real(bz%nspins,dp) ! per spin
                enddo
                enddo
                enddo
                !
                ! </projections>
                !
            enddo
            enddo
            write(*,*) bz%Ep(i), sum(bz%pdos(:,:,:,i)), bz%dos(i)
        enddo
        !
    end subroutine tetrahedra_pdos
    
    function       tetrahedron_weights(x,xc,integration_type,apply_blochl_corrections) result(w)
        !
        ! integration_type = heavi/delta
        !
        use am_rank_and_sort
        !
        implicit none
        !
        real(dp), intent(in) :: xc(4)
        real(dp) :: w(4)
        character(len=5), intent(in) :: integration_type
        logical, intent(in) :: apply_blochl_corrections
        real(dp) :: x,x1,x2,x3,x4,f
        real(dp) :: xi(4)
        integer :: bracket
        integer :: i
        !
        ! sort xc into xi
        xi = xc([4:1:-1])
        if (xi(1).gt.xi(2)) then; xi([1,2]) = xi([2,1]); endif
        if (xi(3).gt.xi(4)) then; xi([3,4]) = xi([4,3]); endif
        if (xi(1).gt.xi(3)) then; xi([1,3]) = xi([3,1]); endif
        if (xi(2).gt.xi(4)) then; xi([2,4]) = xi([4,2]); endif
        if (xi(2).gt.xi(3)) then; xi([2,3]) = xi([3,2]); endif
        x1=xi(1); x2=xi(2); x3=xi(3); x4=xi(4)
        !
        ! if (x1.gt.x2) then; call am_print('ERROR','Sorting failed.',flags='E'); stop; endif
        ! if (x2.gt.x3) then; call am_print('ERROR','Sorting failed.',flags='E'); stop; endif
        ! if (x3.gt.x4) then; call am_print('ERROR','Sorting failed.',flags='E'); stop; endif
        !
        bracket = 0
        if     (x.lt.x1) then; bracket = 1 !; call am_print('bracket',xi-x)
        elseif (x.lt.x2) then; bracket = 2 !; call am_print('bracket',xi-x)
        elseif (x.lt.x3) then; bracket = 3 !; call am_print('bracket',xi-x)
        elseif (x.lt.x4) then; bracket = 4 !; call am_print('bracket',xi-x)
        elseif (x.gt.x4) then; bracket = 5 !; call am_print('bracket',xi-x)
        else
            call am_print('ERROR','Unable to bracket.',flags='E')
            stop
        endif
        !
        select case(integration_type)
        case('heavi')
            select case(bracket)
                case (1)
                    w(1) = 0.0_dp
                    w(2) = 0.0_dp
                    w(3) = 0.0_dp
                    w(4) = 0.0_dp
                    f    = 0.0_dp
                case (2)
                    w(1) = -((x-x1)**3*((x-x1)*(1.0_dp/(x1-x2)+1.0_dp/(x1-x3)+1.0_dp/(x1-x4))+4.0_dp))/((x1-x2)*(x1-x3)*(x1-x4))
                    w(2) =  ((x-x1)**4*1.0_dp/(x1-x2)**2)/((x1-x3)*(x1-x4))
                    w(3) =  ((x-x1)**4*1.0_dp/(x1-x3)**2)/((x1-x2)*(x1-x4))
                    w(4) =  ((x-x1)**4*1.0_dp/(x1-x4)**2)/((x1-x2)*(x1-x3))
                    f    =  ((x-x1)**2*(-1.2D1))/((x1-x2)*(x1-x3)*(x1-x4))
                case (3)
                    w(1) =  (x-x1)**2/((x1-x3)*(x1-x4))+(((x-x1)**2/((x1-x3)*(x1-x4))+((x-x1)*(x-x2)*(x-x3))/((x1-x3)*(x1-x4)*(x2-x3)))*(x-x3))/(x1-x3)+((x-x4)*((x-x1)**2/((x1-x3)*(x1-x4))+((x-x2)**2*(x-x4))/((x1-x4)*(x2-x3)*(x2-x4))+((x-x1)*(x-x2)*(x-x3))/((x1-x3)*(x1-x4)*(x2-x3))))/(x1-x4)
                    w(2) =  ((((x-x2)**2*(x-x4))/((x1-x4)*(x2-x3)*(x2-x4))+((x-x1)*(x-x2)*(x-x3))/((x1-x3)*(x1-x4)*(x2-x3)))*(x-x3))/(x2-x3)+(x-x1)**2/((x1-x3)*(x1-x4))+((x-x2)**2*(x-x4))/((x1-x4)*(x2-x3)*(x2-x4))+((x-x2)**2*(x-x4)**2*1.0_dp/(x2-x4)**2)/((x1-x4)*(x2-x3))+((x-x1)*(x-x2)*(x-x3))/((x1-x3)*(x1-x4)*(x2-x3))
                    w(3) = -((((x-x2)**2*(x-x4))/((x1-x4)*(x2-x3)*(x2-x4))+((x-x1)*(x-x2)*(x-x3))/((x1-x3)*(x1-x4)*(x2-x3)))*(x-x2))/(x2-x3)-(((x-x1)**2/((x1-x3)*(x1-x4))+((x-x1)*(x-x2)*(x-x3))/((x1-x3)*(x1-x4)*(x2-x3)))*(x-x1))/(x1-x3)
                    w(4) = -((x-x1)*((x-x1)**2/((x1-x3)*(x1-x4))+((x-x2)**2*(x-x4))/((x1-x4)*(x2-x3)*(x2-x4))+((x-x1)*(x-x2)*(x-x3))/((x1-x3)*(x1-x4)*(x2-x3))))/(x1-x4)-((x-x2)**3*(x-x4)*1.0_dp/(x2-x4)**2)/((x1-x4)*(x2-x3))
                    f    =  (x*2.4D1-x1*1.2D1-x2*1.2D1+((x-x2)**2*(x1*3.0_dp+x2*3.0_dp-x3*3.0_dp-x4*3.0_dp)*4.0_dp)/((x2-x3)*(x2-x4)))/((x1-x3)*(x1-x4))
                case (4)
                    w(1) = -((x-x4)**4*1.0_dp/(x1-x4)**2)/((x2-x4)*(x3-x4))+1.0_dp
                    w(2) = -((x-x4)**4*1.0_dp/(x2-x4)**2)/((x1-x4)*(x3-x4))+1.0_dp
                    w(3) = -((x-x4)**4*1.0_dp/(x3-x4)**2)/((x1-x4)*(x2-x4))+1.0_dp
                    w(4) =  ((x-x4)**3*((x-x4)*(1.0_dp/(x1-x4)+1.0_dp/(x2-x4)+1.0_dp/(x3-x4))-4.0_dp))/((x1-x4)*(x2-x4)*(x3-x4))+1.0_dp
                    f    =  ((x-x4)**2*(-1.2D1))/((x1-x4)*(x2-x4)*(x3-x4))
                case (5)
                    w(1) = 1.0_dp
                    w(2) = 1.0_dp
                    w(3) = 1.0_dp
                    w(4) = 1.0_dp
                    f    = 0.0_dp
                case default
                    call am_print('ERROR','Error computing tetrahedron weight.',flags='E')
                    call am_print('x',x)
                    call am_print('xc',[x1,x2,x3,x4])
                    call am_print('w',w)
                    call am_print('f',f)
                    stop                        
                end select
        case('delta')
            select case(bracket)
                case (1)
                    w(1) = 0.0_dp
                    w(2) = 0.0_dp
                    w(3) = 0.0_dp
                    w(4) = 0.0_dp
                    f    = 0.0_dp
                case (2)
                    w(1) =   (x-x1)**2*1.0_dp/(x1-x2)**2*1.0_dp/(x1-x3)**2*1.0_dp/(x1-x4)**2*(x*(x4*(x2+x3)+x2*x3-x1*(x2+x3+x4)*2.0_dp+x1**2*3.0_dp)+x1*(x4*(x2+x3)+x2*x3)*2.0_dp-x1**2*(x2+x3+x4)-x2*x3*x4*3.0_dp)*(-4.0_dp)
                    w(2) =  ((x-x1)**3*1.0_dp/(x1-x2)**2*4.0_dp)/((x1-x3)*(x1-x4))
                    w(3) =  ((x-x1)**3*1.0_dp/(x1-x3)**2*4.0_dp)/((x1-x2)*(x1-x4))
                    w(4) =  ((x-x1)**3*1.0_dp/(x1-x4)**2*4.0_dp)/((x1-x2)*(x1-x3))
                    f    = -(x*2.4D1-x1*2.4D1)/((x1-x2)*(x1-x3)*(x1-x4))
                case (3)
                    w(1) = (1.0_dp/(x1-x3)**2*1.0_dp/(x1-x4)**2*(x**2*(x1**2*x2*3.0_dp+x3*x4**2*3.0_dp+x3**2*x4*3.0_dp-x1*x3*x4*6.0_dp-x2*x3*x4*3.0_dp)-x**3*(x1*x2*2.0_dp-x1*x3*2.0_dp-x1*x4*2.0_dp-x2*x3-x2*x4+x3*x4+x1**2+x3**2+x4**2)-x*(x3**2*x4**2*3.0_dp+x1**2*x2*x3*3.0_dp+x1**2*x2*x4*3.0_dp-x1**2*x3*x4*3.0_dp-x1*x2*x3*x4*6.0_dp)+x1**2*x2*x3**2+x1**2*x2*x4**2+x1*x3**2*x4**2*2.0_dp-x1**2*x3*x4**2-x1**2*x3**2*x4+x2*x3**2*x4**2-x1*x2*x3*x4**2*2.0_dp-x1*x2*x3**2*x4*2.0_dp+x1**2*x2*x3*x4)*(-4.0_dp))/((x2-x3)*(x2-x4))
                    w(2) = (1.0_dp/(x2-x3)**2*1.0_dp/(x2-x4)**2*(x**2*(x1*x2**2*3.0_dp+x3*x4**2*3.0_dp+x3**2*x4*3.0_dp-x1*x3*x4*3.0_dp-x2*x3*x4*6.0_dp)-x**3*(x1*x2*2.0_dp-x1*x3-x1*x4-x2*x3*2.0_dp-x2*x4*2.0_dp+x3*x4+x2**2+x3**2+x4**2)-x*(x3**2*x4**2*3.0_dp+x1*x2**2*x3*3.0_dp+x1*x2**2*x4*3.0_dp-x2**2*x3*x4*3.0_dp-x1*x2*x3*x4*6.0_dp)+x1*x2**2*x3**2+x1*x2**2*x4**2+x1*x3**2*x4**2+x2*x3**2*x4**2*2.0_dp-x2**2*x3*x4**2-x2**2*x3**2*x4-x1*x2*x3*x4**2*2.0_dp-x1*x2*x3**2*x4*2.0_dp+x1*x2**2*x3*x4)*(-4.0_dp))/((x1-x3)*(x1-x4))
                    w(3) = (1.0_dp/(x1-x3)**2*1.0_dp/(x2-x3)**2*(x**2*(x1*x2**2*3.0_dp+x1**2*x2*3.0_dp+x3**2*x4*3.0_dp-x1*x2*x3*6.0_dp-x1*x2*x4*3.0_dp)-x**3*(x1*x2-x1*x3*2.0_dp-x1*x4-x2*x3*2.0_dp-x2*x4+x3*x4*2.0_dp+x1**2+x2**2+x3**2)-x*(x1**2*x2**2*3.0_dp-x1*x2*x3**2*3.0_dp+x1*x3**2*x4*3.0_dp+x2*x3**2*x4*3.0_dp-x1*x2*x3*x4*6.0_dp)-x1*x2**2*x3**2-x1**2*x2*x3**2+x1**2*x2**2*x3*2.0_dp+x1**2*x2**2*x4+x1**2*x3**2*x4+x2**2*x3**2*x4+x1*x2*x3**2*x4-x1*x2**2*x3*x4*2.0_dp-x1**2*x2*x3*x4*2.0_dp)*4.0_dp)/((x1-x4)*(x2-x4))
                    w(4) = (1.0_dp/(x1-x4)**2*1.0_dp/(x2-x4)**2*(x**2*(x1*x2**2*3.0_dp+x1**2*x2*3.0_dp+x3*x4**2*3.0_dp-x1*x2*x3*3.0_dp-x1*x2*x4*6.0_dp)-x**3*(x1*x2-x1*x3-x1*x4*2.0_dp-x2*x3-x2*x4*2.0_dp+x3*x4*2.0_dp+x1**2+x2**2+x4**2)-x*(x1**2*x2**2*3.0_dp-x1*x2*x4**2*3.0_dp+x1*x3*x4**2*3.0_dp+x2*x3*x4**2*3.0_dp-x1*x2*x3*x4*6.0_dp)+x1**2*x2**2*x3-x1*x2**2*x4**2-x1**2*x2*x4**2+x1**2*x2**2*x4*2.0_dp+x1**2*x3*x4**2+x2**2*x3*x4**2+x1*x2*x3*x4**2-x1*x2**2*x3*x4*2.0_dp-x1**2*x2*x3*x4*2.0_dp)*4.0_dp)/((x1-x3)*(x2-x3))
                    f    = (x1*x2*(-2.4D1)+x3*x4*2.4D1+x*(x1+x2-x3-x4)*2.4D1)/((x1-x3)*(x1-x4)*(x2-x3)*(x2-x4))
                case (4)
                    w(1) =  ((x-x4)**3*1.0_dp/(x1-x4)**2*(-4.0_dp))/((x2-x4)*(x3-x4))
                    w(2) =  ((x-x4)**3*1.0_dp/(x2-x4)**2*(-4.0_dp))/((x1-x4)*(x3-x4))
                    w(3) =  ((x-x4)**3*1.0_dp/(x3-x4)**2*(-4.0_dp))/((x1-x4)*(x2-x4))
                    w(4) =   (x-x4)**2*1.0_dp/(x1-x4)**2*1.0_dp/(x2-x4)**2*1.0_dp/(x3-x4)**2*(x*(x1*x2+x1*x3+x2*x3-x4*(x1*2.0_dp+x2*2.0_dp+x3*2.0_dp)+x4**2*3.0_dp)+x4*(x1*(x2+x3)*2.0_dp+x2*x3*2.0_dp)-x4**2*(x1+x2+x3)-x1*x2*x3*3.0_dp)*4.0_dp
                    f    = -(x*2.4D1-x4*2.4D1)/((x1-x4)*(x2-x4)*(x3-x4))
                case (5)
                    w(1) = 0.0_dp
                    w(2) = 0.0_dp
                    w(3) = 0.0_dp
                    w(4) = 0.0_dp
                    f    = 0.0_dp
                case default
                    call am_print('ERROR','Error computing tetrahedron weight.',flags='E')
                    call am_print('x',x)
                    call am_print('xc',[x1,x2,x3,x4])
                    call am_print('w',w)
                    call am_print('f',f)
                    stop                        
                end select
            case default
                call am_print('ERROR','Integration type not valid',flags='E')
                stop
        end select
        if (any(isnan(w))) then
            call am_print('ERROR','NaN weight returned.',flags='E')
            call am_print('integration_type',integration_type)
            call am_print('bracket',bracket)
            call am_print('xc',[x1,x2,x3,x4])
            call am_print('x',x)
            call am_print('w',w)
            stop
        endif
        if ( apply_blochl_corrections ) then
            f = f*0.25_dp
            do i = 1,4
                w(1)=w(1)+0.025_dp*f*(xc(i)-x1)
                w(2)=w(2)+0.025_dp*f*(xc(i)-x2)
                w(3)=w(3)+0.025_dp*f*(xc(i)-x3)
                w(4)=w(4)+0.025_dp*f*(xc(i)-x4)
            enddo
        endif
        !
        w=w*0.25_dp
        !
    end function   tetrahedron_weights


!     subroutine get_tetrahedra(tet,bz,pc,pg,tet,iopts)
!         !
!         use am_unit_cell
!         use am_tet_mesh
!         use am_options
!         !
!         implicit none
!         !
!         class(am_class_tetrahedra), intent(inout) :: tet
!         type(am_class_bz), intent(in) :: bz
!         type(am_class_prim_cell) , intent(in) :: pc
!         type(am_class_symmetry), intent(in) :: pg
!         type(am_class_bz) :: fbz
!         type(am_class_options), intent(in), optional :: iopts
!         type(am_class_options) :: opts
!         if (present(iopts)) then 
!            opts = iopts
!         else
!            call opts%defaults
!         endif
!         !
!         call bz%expand_to_fbz(pc=pc,fbz=fbz,pg=pg,iopts=opts)
!         !
!         ! call tessellate_mesh(kpt=fbz%kpt,tet=tet%tet,volume=tet%volume)
!         !
        
!     end subroutine get_tetrahedra

end module am_tetrahedra
