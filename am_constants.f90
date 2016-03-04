module am_constants
    !
    implicit none
	!
	! define fundamental parameters
	!
    integer , parameter :: sp=kind(0.0E0) !> single precision
    integer , parameter :: dp=kind(0.0D0) !> double precision
    real(dp), parameter :: tiny=1E-6_dp
    !
end module am_constants