#:include "fypp_macros.fpp"
module am_constants
    !
    implicit none
	!
	! define fundamental parameters
	!
    integer    , parameter :: sp      = kind(0.0E0) !> single precision
    integer    , parameter :: dp      = kind(0.0D0) !> double precision
    real(dp)   , parameter :: tiny    =  1.0D-5     ! ultimately this is determined by to how many decimals atomic basis in poscar is specified to 
    real(dp)   , parameter :: halfpi  =  1.570796326794897
    real(dp)   , parameter :: sqrtpi  =  1.772453850905516
    real(dp)   , parameter :: pi      =  3.141592653589793
    real(dp)   , parameter :: twopi   =  6.283185307179586
    real(dp)   , parameter :: fourpi  = 12.566370614359172
    complex(dp), parameter :: cmplx_i = cmplx(0,1,dp)
    complex(dp), parameter :: itwopi  = cmplx_i*twopi
    ! magic numbers
    integer    , parameter :: default_ps_id_value = 100 ! must be equal to or larger than 11
    integer    , parameter :: column_width = 150
    integer    , parameter :: string_length_schoenflies = 6
    integer    , parameter :: max_argument_length = 500
    integer    , parameter :: maximum_buffer_size = 500

    ! color

    #:if DEBUG > 0
    character(18) , parameter :: flare = char(27)//'[1;31m'//' ... '//char(27)//'[0;0m'
    #:else
    character(5)  , parameter :: flare = ' ... '
    #:endif

    ! DEBUG

    #:if DEBUG > 0
    logical, parameter :: debug = .true.
    #:else
    logical, parameter :: debug = .false.
    #:endif

end module am_constants