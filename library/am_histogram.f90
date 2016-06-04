module am_histogram
    !
    use am_constants
    use am_stdout
    !
    implicit none
    !
    private
    !
    public :: histogram
    public :: plot_histogram
    !
    interface histogram
        module procedure histogram_double, histogram_integer
    end interface ! histogram
    !
    contains

    pure function   histogram_double(x,m) result(binned_data)
        !
        ! m number of bins
        ! x vector to be binned
        !
        implicit none
        real(dp), intent(in) :: x(:) ! data to convert to histogram
        integer , intent(in) :: m ! number of divisors
        integer, allocatable :: binned_data(:)
        integer :: n ! size(x)
        integer, allocatable :: ranges(:)
        integer :: i, j
        !
        n = size(x)
        !
        ranges = linspace(minval(x),maxval(x),m)
        !
        allocate(binned_data(m))
        !
        binned_data = 0
        !
        do i = 1, n ! for each input score
        do j = 1, m ! determine the binned_data
           if (x(i) .lt. ranges(j)) then
              binned_data(j) = binned_data(j) + 1
              exit
           end if
        end do
        end do
        !
        binned_data = nint(10.0_dp*binned_data/maxval(binned_data))
        !
    end function    histogram_double

    function        histogram_integer(x,m) result(binned_data)
        ! --------------------------------------------------------------------
        ! subroutine  distribute() :
        !    this subroutine receives a score array and a range array, binned_datas
        ! the number of each scores in each range, and calls plot() to print
        ! a histogram.
        ! --------------------------------------------------------------------
        implicit none
        integer, intent(in) :: x(:) ! data to convert to histogram
        integer, intent(in) :: m ! number of divisors
        integer, allocatable :: binned_data(:)
        integer :: n ! size(x)
        integer, allocatable :: ranges(:)
        integer :: i, j
        !
        n = size(x)
        !
        ranges = linspace(minval(x),maxval(x),m)
        !
        allocate(binned_data(m))
        !
        binned_data = 0
        !
        do i = 1, n ! for each input score
        do j = 1, m ! determine the binned_data
           if (x(i) .lt. ranges(j)) then
              binned_data(j) = binned_data(j) + 1
              exit
           end if
        end do
        end do
        !
        binned_data = nint(10.0_dp*binned_data/maxval(binned_data))
        !
    end function    histogram_integer

    subroutine      plot_histogram(binned_data)
        ! --------------------------------------------------------------------
        ! subroutine  plot() :
        !    this subroutine receives a binned_data array and prints a vertical
        ! bar histogram.
        ! --------------------------------------------------------------------
        !
        implicit none
        !
        integer, intent(in) :: binned_data(:)
        integer :: m ! size(binned_data)
        character(len=1), allocatable :: line(:)
        character(len=1) :: division
        character(len=1) :: empty
        character(len=1) :: datahas
        character(len=1) :: last
        integer          :: i, j, maximum
        character(len=99):: fmt
        !
        m = size(binned_data)+1
        allocate(line(m))
        !
        division  = "-" ! = "-+"
        empty     = " " ! = "  "
        datahas   = "*" ! = "* "
        last      = "|" ! = "*|"
        !
        maximum = binned_data(1) ! find the maximum of the binned_data
        line(1) = division ! clear the print line
        do i = 2, (m-1)
           line(i) = division
           if (maximum < binned_data(i)) maximum = binned_data(i)
        end do
        line(m:m) = "+"
        !
        fmt = '(1000a)'
        !
        ! print the top border
        write(*,fmt) "+", line
        ! print from the top
        do i = maximum,1,-1
         ! for each binned_data value
           do j = 1,m
                ! if >= current value, show ***
                if (j.lt.m) then
                  if (binned_data(j) >= i) then
                    line(j) = datahas
                else   
                    line(j) = empty
                endif
              else
                  line(j) = last
              endif
           end do
           write(*,fmt) "|", (line(j), j=1,m)
        end do
        !
        do j = 1, m
           line(j) = division
        end do
        write(*,fmt) "+", line
    end subroutine  plot_histogram



end module am_histogram