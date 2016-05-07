module am_constants
    !
    implicit none
	!
	! define fundamental parameters
	!
    integer    , parameter :: sp=kind(0.0E0) !> single precision
    integer    , parameter :: dp=kind(0.0D0) !> double precision
    integer    , parameter :: qp=selected_real_kind(32) !> quadrupole precision (not in use)
    real(dp)   , parameter :: tiny=1E-6_dp
    real(dp)   , parameter :: halfpi=1.570796326794897
    real(dp)   , parameter :: pi=3.141592653589793
    real(dp)   , parameter :: twopi=6.283185307179586
    real(dp)   , parameter :: fourpi=12.566370614359172
    complex(dp), parameter :: cmplx_i=(0.0_dp,1.0_dp)
    complex(dp), parameter :: itwopi=cmplx_i*twopi
    !
    ! magic numbers
    integer    , parameter :: column_width = 100
    integer    , parameter :: string_length_schoenflies = 6
    integer    , parameter :: max_argument_length = 500
    integer    , parameter :: maximum_buffer_size = 500
    
    
end module am_constants