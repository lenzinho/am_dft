#:include "fypp_macros.fpp"
module am_dispersion

    use dispmodule
    use am_brillouin_zone
    use am_unit_cell
    use am_constants
    use am_stdout
    use am_vasp_io
    use am_mkl
    use am_options
    use am_histogram
    use am_shells

    implicit none

    private

    type, public :: am_class_dispersion
        integer :: nbands
        integer :: nspins
        real(dp), allocatable :: E(:,:) ! E(nbands,nkpts) energies
    contains
        procedure :: debug_dump => debug_dump_dr
    end type am_class_dispersion

contains

    subroutine     debug_dump_dr(dr,fname)
        !
        implicit none
        !
        class(am_class_dispersion), intent(in) :: dr
        character(*), intent(in) :: fname
        !
                                                call dump(A=dr%nbands            ,fname=trim(fname)//'.nbands'           )
        if (allocated(dr%E                   )) call dump(A=dr%E                 ,fname=trim(fname)//'.E'                )
        ! ! dump dr specific data                                                                                        
        ! select type (dr)                                                                                         
        ! class is (am_class_dispersion)                                                                                    
        ! class is (am_class_seitz_group)                                                                           
        ! class default
        !     stop 'ERROR [debug_dump]: invalid group class'
        ! end select

        !
    end subroutine debug_dump_dr







end module am_dispersion