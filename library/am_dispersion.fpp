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
        real(dp), allocatable :: E(:,:) ! E(nbands,nkpts) energies
    end type am_class_dispersion

contains

end module am_dispersion