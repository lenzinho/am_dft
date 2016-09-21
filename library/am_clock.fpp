module am_clock
    !
    use am_constants,only : dp
    !
    implicit none
    !
    save
    !
    integer, parameter :: maxclock = 101
    real(dp),parameter :: notrunning = - 1.0_dp
    real(dp) :: cputime(maxclock),t0cpu(maxclock)
    real(dp) :: walltime(maxclock),t0wall(maxclock)
    character(len=12) :: clock_label(maxclock)
    integer :: called(maxclock)
    integer :: nclock = 0
    logical :: no
    integer :: trace_depth = 0
    !
end module am_clock

    function       t_wall()
        !
        use am_constants,only : dp
        !
        implicit none
        !
        integer :: t_wall_
        real(dp) :: t_wall
        !
        call system_clock(t_wall_)
        t_wall = real(t_wall_,dp)*1.0e-4_dp
        !
    end function   t_wall

    function       t_cpu()
        !
        use am_constants,only : dp
        !
        implicit none
        !
        real(dp) :: t_cpu
        !
        call cpu_time(t_cpu)
        !
    end function   t_cpu

    subroutine     init_clocks(go)
        ! ... go = .true.   : clocks will run
        ! ... go = .false. : only clock #1 will run
        use am_constants,   only : dp
        use am_clock,only : called,t0cpu,cputime,no,notrunning,maxclock,clock_label,walltime,t0wall,nclock
        implicit none
        logical :: go
        integer :: n
        no = .not. go
        nclock = 0
        do n = 1,maxclock
            called(n)   = 0
            cputime(n)  = 0.0_dp
            t0cpu(n)    = notrunning
            walltime(n) = 0.0_dp
            t0wall(n)   = notrunning
            clock_label(n) = ' '
        enddo
    end subroutine init_clocks

    subroutine     start_clock(label)
        !
        use am_constants,   only : dp
        use am_clock,       only : trace_depth
        use am_clock,       only : nclock,clock_label,notrunning,no,maxclock,t0cpu,t0wall
        !
        implicit none
        !
        character(len=*)  :: label
        character(len=12) :: label_
        integer           :: n
        real(dp),external :: t_cpu,t_wall
        !
        ! write(*,'("trace (depth=",i2,") start: ",a12)') trace_depth,label
        !
        trace_depth = trace_depth + 1
        !
        if (no .and. (nclock == 1)) return
        ! ... prevent trouble if label is longer than 12 characters
        label_ = trim (label)
        do n = 1,nclock
            if (clock_label(n) == label_) then
                if (t0cpu(n) /= notrunning) then
                    write(*,'("start_clock: clock # ",i2," for ",a12," already started")') n,label_
                else
                    t0cpu(n) = t_cpu()
                    t0wall(n) = t_wall()
                endif
                return
            endif
        enddo
        !
        if (nclock == maxclock) then
            write(*,'("start_clock(",a,"): too many clocks! call ignored")') label
        else
            nclock              = nclock + 1
            clock_label(nclock) = label_
            t0cpu(nclock)       = t_cpu()
            t0wall(nclock)      = t_wall()
        endif
    end subroutine start_clock

    subroutine     stop_clock(label)
        !
        use am_constants,  only : dp
        use am_clock,      only : trace_depth
        use am_clock,      only : no,nclock,clock_label,cputime,walltime,notrunning,called,t0cpu,t0wall
        !
        implicit none
        !
        character(len=*)  :: label
        character(len=12) :: label_
        integer           :: n
        real(dp),external :: t_cpu,t_wall
        !
        trace_depth = trace_depth - 1
        ! write(*,'("trace (depth=",i2,") end: ",a12)') trace_depth,label
        if (no) return
        ! ... prevent trouble if label is longer than 12 characters
        label_ = trim (label)
        do n = 1,nclock
            if (clock_label(n) == label_) then
                if (t0cpu(n) == notrunning) then
                    write(*,'("stop_clock: clock # ",i2," for ",a12," not running")') n,label
                else
                    cputime(n)  = cputime(n) + t_cpu() - t0cpu(n)
                    walltime(n) = walltime(n) + t_wall() - t0wall(n)
                    t0cpu(n)    = notrunning
                    t0wall(n)   = notrunning
                    called(n)   = called(n) + 1
                endif
                return
            endif
        enddo
        ! ... clock not found
        write(*,'("stop_clock: no clock for ",a12," found !")') label
    end subroutine stop_clock

    subroutine     print_clock(label)
        !
        use am_clock,       only : nclock,clock_label
        !
        implicit none
        !
        character(len=*)   :: label
        character(len=12)  :: label_
        integer            :: n
        !
        if (trim(label).eq.'') then
            do n = 1,nclock
                call print_this_clock(n)
            enddo
        else
            ! ... prevent trouble if label is longer than 12 characters
            label_ = trim (label)
            do n = 1,nclock
                if (clock_label(n) == label_) then
                    call print_this_clock(n)
                    exit
                endif
            enddo
        endif
    end subroutine print_clock

    subroutine     print_this_clock(n)
        !
        use am_constants,   only : dp
        use am_clock,       only : clock_label,cputime,walltime,notrunning,called,t0cpu,t0wall
        !
        implicit none
        !
        integer :: n
        real(dp) :: elapsed_cpu_time,elapsed_wall_time,nsec,msec
        integer :: nday,nhour,nmin,nmax,mday,mhour,mmin
        real(dp),external :: t_cpu,t_wall
        !
        if (t0cpu(n) == notrunning) then
            ! ... clock stopped,print the stored value for the cpu time
            elapsed_cpu_time = cputime(n)
            elapsed_wall_time= walltime(n)
        else
            ! ... clock not stopped,print the current value of the cpu time
            elapsed_cpu_time   = cputime(n) + t_cpu() - t0cpu(n)
            elapsed_wall_time  = walltime(n) + t_wall() - t0wall(n)
            called(n)          = called(n) + 1
        endif
        nmax = called(n)
        if (n == 1) then
            nday    = elapsed_cpu_time / 86400
            nsec    = elapsed_cpu_time - 86400 * nday
            nhour   = nsec / 3600
            nsec    = nsec - 3600 * nhour
            nmin    = nsec / 60
            nsec    = nsec - 60 * nmin
            ! ... the first clock writes elapsed (wall) time as well
            mday    = elapsed_wall_time / 86400
            msec    = elapsed_wall_time - 86400 * mday
            mhour   = msec / 3600
            msec    = msec - 3600 * mhour
            mmin    = msec / 60
            msec    = msec - 60 * mmin
            if (nday > 0 .or. mday > 0) then
                write(*,'(5x,a12," : ",3x,i2,"d",3x,i2,"h",i2,"m cpu ","     ",3x,i2,"d",3x,i2,"h",i2,"m wall")') clock_label(n),nday,nhour,nmin,mday,mhour,mmin
            elseif (nhour > 0 .or. mhour > 0) then
                write(*,'(5x,a12," : ",3x,i2,"h",i2,"m cpu ","   ",3x,i2,"h",i2,"m wall")') clock_label(n),nhour,nmin,mhour,mmin
            elseif (nmin > 0 .or. mmin > 0) then
                write(*,'(5x,a12," : ",i2,"m",f5.2,"s cpu ","    ",i2,"m",f5.2,"s wall")') clock_label(n),nmin,nsec,mmin,msec
            else
                write(*,'(5x,a12," : ",f9.2,"s cpu ",f9.2,"s wall")') clock_label(n),nsec,msec
            endif
        elseif (nmax == 1 .or. t0cpu(n) /= notrunning) then
            ! ... for clocks that have been called only once
            write(*,'(5x,a12," : ",f9.2,"s cpu ",f9.2,"s wall",i10," calls")') clock_label(n),elapsed_cpu_time,elapsed_wall_time,nmax
        elseif (nmax == 0) then
            ! ... for clocks that have never been called
            write(*,'("clock # ",i2," for ",a12," never called !"/)') n,clock_label(n)
        else
            ! ... for all other clocks
            write(*,'(5x,a12," : ",f9.2,"s cpu ",f9.2,"s wall",i10," calls")') clock_label(n),elapsed_cpu_time,elapsed_wall_time,nmax
        endif
    end subroutine print_this_clock



