module am_constants
    !
    implicit none
	!
	! define fundamental parameters
	!
    integer    , parameter :: sp      = kind(0.0E0) !> single precision
    integer    , parameter :: dp      = kind(0.0D0) !> double precision
    real(dp)   , parameter :: tiny    = 1.0D-5 ! ultimately this is determined by to how many decimals atomic basis in poscar is specified to 
    real(dp)   , parameter :: halfpi  = 1.570796326794897
    real(dp)   , parameter :: pi      = 3.141592653589793
    real(dp)   , parameter :: twopi   = 6.283185307179586
    real(dp)   , parameter :: fourpi  = 12.566370614359172
    complex(dp), parameter :: cmplx_i = (0.0_dp,1.0_dp)
    complex(dp), parameter :: itwopi  = cmplx_i*twopi
    !
    ! magic numbers
    integer    , parameter :: column_width = 100
    integer    , parameter :: string_length_schoenflies = 6
    integer    , parameter :: max_argument_length = 500
    integer    , parameter :: maximum_buffer_size = 500
    
    
end module am_constants