module am_stopwatch
	!
	! borrowed from wannier90
	!
	use am_constants, only : dp
  use am_helpers
	!
	implicit none
	!
	private
	!
	public :: io_stopwatch
  public :: io_print_timings
	!
  integer, parameter :: max_timers = 100
  integer, parameter :: max_calls_to_stopwatch = 100
  !
  type am_class_timers
     integer :: ncalls
     real(kind=dp) :: ctime
     real(kind=dp) :: ptime
     character(len=100) :: label
  end type am_class_timers
  
	contains

    subroutine initialize(timers)
      !
      implicit none
      !
      class(am_class_timers) :: timers(max_timers)
      !


    subroutine stopwatch(tag,mode)
        !
        implicit none
        class(am_class_timers) ;: timer
        character(len=*), intent(in) :: tag
        integer, intent(in)          :: mode
        integer :: i
        real(kind=dp) :: t
        !
        call cpu_time(t)
        !
        select case (mode)
        case (1)
            do i=1, nnames
                if (timer(i)%label .eq. tag) then
                    timer(i)%ptime  = t
                    timer(i)%ncalls = timer(i)%ncalls + 1
                    return
                endif
            enddo
            nnames = nnames + 1
            if (nnames .gt. max_calls_to_stopwatch) call io_error('Maximum number of calls to io_stopwatch exceeded')
            timer(nnames)%label = tag
            timer(nnames)%ctime = 0.0_dp
            timer(nnames)%ptime = t
            timer(nnames)%ncalls = 1
        case (2)
            do i=1,nnames
                if (timer(i)%label .eq. tag) then
                    timer(i)%ctime = timer(i)%ctime + t - timer(i)%ptime
                    return
                endif
            end do
            write(stdout,'(1x,3a)') 'WARNING: name = ',trim(tag),' not found in io_stopwatch' 
        case default
            write(stdout,*) ' Name = ',trim(tag),' mode = ',mode
            call am_print('ERROR','Value of mode not recognised in io_stopwatch',' >>> ')
            stop
        end select
    end subroutine stopwatch

    subroutine print_timings(timers)
        implicit none
        integer :: i
        !
        call am_print_title('Timing information')
        write(stdout,'(a5,a60,a20,a20)') ' ... ', 'Tag','Ncalls','Time [s]'
        do i=1,nnames
            write(stdout,'(a5,a60,i20,f20.3)') &
            timer(i)%label,timer(i)%ncalls,timer(i)%ctime 
        enddo
    end subroutine print_timings

end module am_stopwatch

