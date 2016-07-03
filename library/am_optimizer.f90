module am_optimizer

    use dispmodule
    use am_constants
    use am_stdout
    use am_matlab
    use am_options
    use am_shells
    use am_symmetry_rep
    use am_symmetry_relations
    use am_mkl
    use am_brillouin_zone
    use am_dispersion
    use am_tight_binding

    type, private :: am_class_tb_optimizer
        integer  :: maxiter  ! maximum number of iterations  
        integer  :: nbands   ! number of bands
        integer  :: nkpts    ! numnbr of kpoints
        integer  :: nxs      ! number of parameters to fit
        integer  :: nrs  	 ! size of residual vector
        real(dp) :: eps      ! stopping criterion
        integer               :: skip_band
        integer , allocatable :: selector_shell(:)
        integer , allocatable :: selector_kpoint(:)
        real(dp), allocatable :: x(:) 
        real(dp), allocatable :: r(:)
        real(dp)              :: rms
    end type am_class_tb_optimizer

    public :: optimize_matrix_elements

contains

    subroutine     optimize_matrix_elements(tb,tbpg,bz,dr_dft,pp,opts)
        !
        implicit none
        !   
        class(am_class_tight_binding) , intent(inout) :: tb
        type(am_class_tb_group)       , intent(in) :: tbpg
        type(am_class_bz)             , intent(in) :: bz
        type(am_class_dispersion_dft) , intent(in) :: dr_dft
        type(am_class_prim_pair)      , intent(in) :: pp
        type(am_class_options)        , intent(in) :: opts
        type(am_class_tb_optimizer) :: ft
        !
        if (opts%verbosity.ge.1) call print_title('Optimized matrix elements')
        !
        ! initialize optimizer
        ft%maxiter    = 10
        ft%nbands     = maxval(tbpg%E)
        ft%nkpts      = bz%nkpts
        ft%nxs        = tb%nVs
        ft%skip_band  = opts%skip_band
        ft%nrs        = bz%nkpts * ft%nbands
        allocate(ft%selector_shell,  source=[1:50]) ! first fifty shells
        allocate(ft%selector_kpoint, source=[1:bz%nkpts])
        allocate(ft%r(ft%nrs))
        ft%r = 0
        allocate(ft%x(ft%nxs))
        ft%x = 0
        !
        if (debug) then
        write(*,'(a,a,a)') flare, 'irreducible matrix elements = '//tostring(ft%nxs)
        write(*,'(a,a,a)') flare, 'k-points = '                   //tostring(ft%nkpts)
        write(*,'(a,a,a)') flare, 'bands = '                      //tostring(ft%nbands)
        write(*,'(a,a,a)') flare, 'number of bands to skip = '    //tostring(ft%skip_band)
        write(*,'(a,a,a)') flare, 'residual vector length = '     //tostring(ft%nrs)
        write(*,'(a,a,a)') flare, 'max iterations = '             //tostring(ft%maxiter)
        endif
        !
        call perform_optimization(ft=ft, tb=tb, tbpg=tbpg, bz=bz, dr_dft=dr_dft, pp=pp)
        !
        ! call dump(ft%x,'outfile.debug.V')
        !
        contains
        subroutine     perform_optimization(ft,tb,tbpg,bz,dr_dft,pp)
            !
            use mkl_rci
            use mkl_rci_type
            !
            implicit none
            !
            class(am_class_tb_optimizer)  , intent(inout) :: ft
            type(am_class_tight_binding)  , intent(inout) :: tb
            type(am_class_tb_group)       , intent(in) :: tbpg
            type(am_class_bz)             , intent(in) :: bz
            type(am_class_dispersion_dft) , intent(in) :: dr_dft
            type(am_class_prim_pair)      , intent(in) :: pp
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
            allocate(x(ft%nxs))
            x = tb%V 
            ! initialize residual vector
            allocate(FVEC(ft%nrs))
            FVEC = 0.0_dp
            ! initialize jacobian matrix
            allocate(FJAC(ft%nrs,ft%nxs))
            FJAC = 0
            ! initializes the solver of a nonlinear least squares problem
            if (dtrnlsp_init(handle=handle, n=ft%nxs, m=ft%nrs, x=x, eps=eps, iter1=ft%maxiter, iter2=100, rs=100.0D0) /= tr_success) then
                call mkl_free_buffers
                stop 'ERROR [perform_optimization]: dtrnlsp_init'
            end if
            ! check the correctness of handle and arrays containing Jacobian matrix, objective function, and stopping criteria
            if (dtrnlsp_check(handle=handle, n=ft%nxs, m=ft%nrs, FJAC=FJAC, FVEC=FVEC, eps=eps, info=info) /= tr_success) then
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
            write(*,'(5x,a8,a10,a10,a)') 'iter', 'rms ', 'log(d rms)',centertitle('parameters',ft%nxs*10)
            write(*,'(5x,a8,a10,a10,a)') ' '//repeat('-',8-1), ' '//repeat('-',10-1), ' '//repeat('-',10-1), ' '//repeat('-',ft%nxs*10-1)
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
                    call tb%set_Vsk(V=x)
                    ! recalculate function
                    call compute_residual(ft=ft, tb=tb, tbpg=tbpg, bz=bz, dr_dft=dr_dft, pp=pp, R=FVEC)
                    ! save result
                    ft%r   = FVEC 
                    ft%rms = sqrt(sum(FVEC**2)/ft%nrs)
                    ! increase counter
                    i = i + 1
                    ! print results
                    if (i.eq.1) then
                        write(*,'(5x,i8,f10.2,10x,SP,1000f10.2)')   i, ft%rms, x
                    else
                        write(*,'(5x,i8,f10.2,f10.2,SP,1000f10.2)') i, ft%rms, log10(abs(norm(1.0D-14 + x-ft%x ))), x
                    endif
                    ! save results
                    ft%x = x
                case (2)
                    ! update x in tb model
                    call tb%set_Vsk(V=x)
                    ! compute jacobian matrix (uses central difference)
                    call compute_jacobian(ft=ft, tb=tb, tbpg=tbpg, bz=bz, dr_dft=dr_dft,  pp=pp, FJAC=FJAC)
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
        subroutine     compute_residual(ft,tb,tbpg,bz,dr_dft,pp,R)
            !
            implicit none
            !
            class(am_class_tb_optimizer)  , intent(in) :: ft
            type(am_class_tight_binding)  , intent(in) :: tb
            type(am_class_tb_group)       , intent(in) :: tbpg
            type(am_class_bz)             , intent(in) :: bz
            type(am_class_dispersion_dft) , intent(in) :: dr_dft 
            type(am_class_prim_pair)      , intent(in) :: pp
            real(dp), allocatable         , intent(out):: R(:)
            type(am_class_dispersion_tb) :: dr_tb 
            integer, allocatable :: inds(:)
            integer :: i
            ! index vector
            allocate(inds(ft%nbands))
            ! create residual vector
            allocate(R(ft%nrs))
            R = 0
            ! get tb dispersion 
            call dr_tb%get_dispersion(tb=tb,tbpg=tbpg,bz=bz,pp=pp)
            ! loop over kpoints
            do i = 1, bz%nkpts
                ! get indices
                inds = [1:ft%nbands] + (i-1)*ft%nbands
                ! calculate residual vector indices
                R(inds) = dr_dft%E([1:ft%nbands]+ft%skip_band, i ) - dr_tb%E(:,i)
            enddo
            !
        end subroutine compute_residual
        subroutine     compute_jacobian(ft,tb,tbpg,bz,dr_dft,pp,FJAC)
            !
            use mkl_rci
            !
            implicit none 
            !
            class(am_class_tb_optimizer)  , intent(in) :: ft
            type(am_class_tight_binding)  , intent(inout) :: tb
            type(am_class_tb_group)       , intent(in) :: tbpg
            type(am_class_bz)             , intent(in) :: bz
            type(am_class_dispersion_dft) , intent(in) :: dr_dft
            type(am_class_prim_pair)      , intent(in) :: pp
            real(dp), allocatable         , intent(out):: FJAC(:,:) ! fjac(m,n) jacobian matrix
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
            allocate(  f1(ft%nrs))
            allocate(  f2(ft%nrs))
            allocate(FJAC(ft%nrs,ft%nxs))
            ! initialize x
            allocate(x, source=tb%V)
            ! begin jacobian computation
            if (djacobi_init(handle=handle, n=ft%nxs, m=ft%nrs, x=x, FJAC=FJAC, eps=eps_jac) .ne. tr_success) then
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
                    call tb%set_Vsk(V=x)
                    ! compute jacobian
                    call compute_residual(ft=ft,tb=tb,tbpg=tbpg,bz=bz,dr_dft=dr_dft,pp=pp,R=f1)
                else if (rci_request .eq. 2) then
                    ! update x in TB model
                    call tb%set_Vsk(V=x)
                    ! compute jacobian
                    call compute_residual(ft=ft,tb=tb,tbpg=tbpg,bz=bz,dr_dft=dr_dft,pp=pp,R=f2)
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
        end subroutine compute_jacobian
    end subroutine optimize_matrix_elements


end module am_optimizer