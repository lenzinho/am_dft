module am_constants
    !
    implicit none
	!
	! define fundamental parameters
	!
    integer , parameter :: sp=kind(0.0E0) !> single precision
    integer , parameter :: dp=kind(0.0D0) !> double precision
    real(dp), parameter :: tiny=1E-6_dp
    real(dp), parameter :: twopi=6.283185307179586
    complex , parameter :: cmplx_i = (0.0_dp,1.0_dp)
    !
    ! magic numbers
    integer , parameter :: column_width = 90
    integer , parameter :: string_length_schoenflies = 6
    integer , parameter :: max_argument_length = 500
    integer , parameter :: maximum_buffer_size = 500
    
    
end module am_constants