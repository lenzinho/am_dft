#:include "fypp_macros.fpp"
module am_tetrahedra

	use am_brillouin_zone
    use am_symmetry
    use am_constants
    use am_stdout
    use am_vasp_io
    use am_mkl
    use am_options

	implicit none

	private

	public :: tetrahedron_weights

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
        if (x1.gt.x2) stop 'ERROR [tetrahedron_weights]: sorting failed'
        if (x2.gt.x3) stop 'ERROR [tetrahedron_weights]: sorting failed'
        if (x3.gt.x4) stop 'ERROR [tetrahedron_weights]: sorting failed'
        !
        bracket = 0
        if     (x.lt.x1) then; bracket = 1
        elseif (x.lt.x2) then; bracket = 2
        elseif (x.lt.x3) then; bracket = 3
        elseif (x.lt.x4) then; bracket = 4
        elseif (x.gt.x4) then; bracket = 5
        else
            stop 'ERROR [tetrahedron_weights]: bracketing failed'
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
                    stop 'ERROR [tetrahedron_weights]: heaviside, invalid bracket'
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
                    stop 'ERROR [tetrahedron_weights]: delta, invalid bracket'
                end select
            case default
                stop 'ERROR [tetrahedron_weights]: integration type /= heavi or delta'
        end select
        ! check for NaNs
        if (any(isnan(w))) stop 'ERROR [tetrahedron_weights]: NaN obtained'
        ! apply Bloch correction
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
!         type(am_class_seitz_group), intent(in) :: pg
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
