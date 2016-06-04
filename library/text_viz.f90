! small module to do text mode graphics using ANSI terminal escape sequences
! Copyright (c) 2013 Axel Kohlmeyer <akohlmey@gmail.com>
! v0.1 2013-06-22, initial public version
! v0.2 2013-06-24, better autoscale and overflow protection

MODULE text_viz
  IMPLICIT NONE

  PRIVATE
  CHARACTER, PARAMETER :: esc = ACHAR(27)
  INTEGER,   PARAMETER :: maxplot = 10
  CHARACTER, PARAMETER :: plot(0:maxplot) = (/ &
      ' ', '.', ',', ':', '=', '+', 'o', 'x', 'X', '#', '@' /)

  INTEGER :: rows, cols

  PUBLIC :: viz_init, viz_clear, viz_pos, viz_plot, viz_done

CONTAINS
  ! initialize package and set size of plot area
  SUBROUTINE viz_init(x,y)
    INTEGER, INTENT(in):: x,y
    cols = x
    rows = y
    ! hide the cursor to reduce flicker
    WRITE(*,FMT='(A1,A)',ADVANCE='NO') esc,'[?25l'
  END SUBROUTINE viz_init

  ! restore settings to something sane.
  SUBROUTINE viz_done
    ! place cursor in last line
    CALL viz_pos(1,rows)
    ! set forground color to principal color
    WRITE(*,FMT='(A1,A,A1)',ADVANCE='NO') esc,'[30m'
    ! re-enable the cursor
    WRITE(*,FMT='(A1,A)', ADVANCE='NO') esc,'[?25h'
    ! and call for a reset
    WRITE(*,FMT='(A1,A)') esc,'[0m'
  END SUBROUTINE viz_done


  ! clear the screen
  SUBROUTINE viz_clear
    WRITE(*,FMT='(A1,A)',ADVANCE='NO') esc,'[2J'
  END SUBROUTINE viz_clear

  SUBROUTINE add_code(code,c,n)
    CHARACTER(len=1),INTENT(inout) :: code(*)
    CHARACTER(len=1),INTENT(in)    :: c
    INTEGER,INTENT(inout)          :: n

    n = n+1
    code(n) = c
  END SUBROUTINE add_code


  ! position cursor at a given location within screen.
  ! top left corner is (1,1)
  SUBROUTINE viz_pos(x,y)
    INTEGER, INTENT(in) :: x,y
    INTEGER :: i,n
    CHARACTER(len=1) :: code(7)

    n = 0

    i = y
    if (i < 1) i = 1
    if (i > rows) i = rows
    CALL add_code(code,'[',n)
    IF (i > 9) CALL add_code(code,ACHAR(48+i/10),n)
    CALL add_code(code,ACHAR(48+MOD(i,10)),n)
    CALL add_code(code,';',n)

    i = x
    if (i < 1) i = 1
    if (i > cols) i = cols
    IF (i > 9) CALL add_code(code,ACHAR(48+i/10),n)
    CALL add_code(code,ACHAR(48+MOD(i,10)),n)
    CALL add_code(code,'H',n)
    WRITE(*,FMT='(A1,A)',ADVANCE='NO') esc,code(1:n)
  END SUBROUTINE viz_pos

  SUBROUTINE viz_plot(val,nx,ny,max)
    REAL, INTENT(in)          :: val(:,:)
    INTEGER, INTENT(in)       :: nx,ny
    REAL, INTENT(in),OPTIONAL :: max
    INTEGER :: i, j, k, l, m, n, dx, dy
    REAL    :: vmax, scalef, tmp

    ! set blocksize for averaging
    dx = nx / cols
    if (dx < 1) dx = 1
    dy = ny / rows
    if (dy < 1) dy = 1

    ! set or determine scaling factor for data points
    vmax=1.0e-30
    IF (PRESENT(max)) THEN
        vmax = ABS(max)
    ELSE
        ! find absolute maximum value for scaling
        DO j=1,rows
            DO i=1,cols
                ! average over cells
                tmp = 0.0
                n = 0
                DO k=(j-1)*dy+1,j*dy
                    DO l=(i-1)*dx+1,i*dx
                        tmp = tmp + val(l,k)
                        n = n + 1
                    END DO
                END DO

                tmp = ABS(tmp)/REAL(n)
                IF (vmax < tmp) vmax = tmp
            END DO
        END DO
    END IF
    scalef = REAL(maxplot)/vmax

    ! now plot
    DO j=1,rows
        CALL viz_pos(1,j)
        DO i=1,cols
            ! average over cells
            tmp = 0.0
            n = 0
            DO k=(j-1)*dy+1,j*dy
                DO l=(i-1)*dx+1,i*dx
                    tmp = tmp + val(l,k)
                    n = n + 1
                END DO
            END DO
            ! convert absolute value into character
            m = INT(scalef*ABS(tmp)/real(n)+0.5)
            if (m > maxplot) m = 10
            IF (tmp < 0.0) THEN
                IF (m > 5) THEN
                    WRITE(*,FMT='(A1,A,A1)',ADVANCE='NO') esc,'[36m',plot(m)
                ELSE
                    WRITE(*,FMT='(A1,A,A1)',ADVANCE='NO') esc,'[34m',plot(m)
                END IF
            ELSE
                IF (m > 5) THEN
                    WRITE(*,FMT='(A1,A,A1)',ADVANCE='NO') esc,'[31m',plot(m)
                ELSE
                    WRITE(*,FMT='(A1,A,A1)',ADVANCE='NO') esc,'[30m',plot(m)
                END IF
            END IF
        END DO
        WRITE(*,FMT='(A1,A)',ADVANCE='NO') esc,'[30m|'
    END DO
    CALL viz_pos(1,rows)
    WRITE(*,FMT='(A1,A)',ADVANCE='NO') esc,'[30m'
    DO i=1,cols+1
        WRITE(*,FMT='(A1)',ADVANCE='no') '-'
    END DO

  END SUBROUTINE viz_plot
END MODULE text_viz

