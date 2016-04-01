module am_rank_and_sort
    !
    ! SORT/RANK routines(BASED ON ORDERPACK 2.0 ROUTINES)
    !
    !
    ! Ranking: In some instances, one is not actually interested in modifying the order of the
    ! elements in a set, but only in knowing how to access them in increasing -- or decreasing --
    ! order. Ranking, as it is called, provides the index array I(:) such as the set S(I(:)) is
    ! ordered. One of the advantages of carrying out ranking rather than sorting is that the index
    ! array can be computed without the performance penalty of moving the elements around when
    ! they are of large sizes. A similar point is that the index array can be used to index other
    ! data.
    !
    use am_constants
    !
    implicit none
    !
    private
    !
    public :: rank
    public :: unique_rank
    public :: sort_using_quicksort
    public :: sort_using_insertion
    
    !

    interface rank
        module procedure d_rank, r_rank, i_rank
    end interface rank

    interface unique_rank
    module procedure d_unique_rank, r_unique_rank, i_unique_rank
    end interface unique_rank

    interface nearless
    module procedure d_nearless, r_nearless, i_nearless
    end interface nearless

    interface sort_using_quicksort
    ! This subroutine uses quick sort.
    module procedure d_sort_using_quicksort, r_sort_using_quicksort, i_sort_using_quicksort
    end interface sort_using_quicksort

    interface sort_using_insertion
    ! This subroutine uses insertion sort. Faster for small arrays (<20).
    module procedure d_sort_using_insertion, r_sort_using_insertion, i_sort_using_insertion
    end interface sort_using_insertion

    contains
    subroutine d_rank(xdont, irngt)
        real(kind=dp), dimension(:), intent(in) :: xdont
        integer, dimension(:), intent(out) :: irngt
        real(kind=dp) :: xvala, xvalb
        integer, dimension(size(irngt)) :: jwrkt
        integer :: lmtna, lmtnc, irng1, irng2
        integer :: nval, iind, iwrkd, iwrk, iwrkf, jinda, iinda, iindb
        nval = min(size(xdont), size(irngt))
        select case(nval)
        case(:0)
            return
        case(1)
            irngt(1) = 1
            return
        case default
            continue
        end select
        do iind = 2, nval, 2
            if(xdont(iind-1) <= xdont(iind)) then
                irngt(iind-1) = iind - 1
                irngt(iind) = iind
            else
                irngt(iind-1) = iind
                irngt(iind) = iind - 1
            end if
        end do
        if(modulo(nval, 2) /= 0) then
            irngt(nval) = nval
        end if
        lmtna = 2
        lmtnc = 4
        do
            if(nval <= 2) exit
            do iwrkd = 0, nval - 1, 4
                if((iwrkd+4) > nval) then
                    if((iwrkd+2) >= nval) exit
                    if(xdont(irngt(iwrkd+2)) <= xdont(irngt(iwrkd+3))) exit
                    if(xdont(irngt(iwrkd+1)) <= xdont(irngt(iwrkd+3))) then
                        irng2 = irngt(iwrkd+2)
                        irngt(iwrkd+2) = irngt(iwrkd+3)
                        irngt(iwrkd+3) = irng2
                    else
                        irng1 = irngt(iwrkd+1)
                        irngt(iwrkd+1) = irngt(iwrkd+3)
                        irngt(iwrkd+3) = irngt(iwrkd+2)
                        irngt(iwrkd+2) = irng1
                    end if
                    exit
                end if
                if(xdont(irngt(iwrkd+2)) <= xdont(irngt(iwrkd+3))) cycle
                if(xdont(irngt(iwrkd+1)) <= xdont(irngt(iwrkd+3))) then
                    irng2 = irngt(iwrkd+2)
                    irngt(iwrkd+2) = irngt(iwrkd+3)
                    if(xdont(irng2) <= xdont(irngt(iwrkd+4))) then
                        irngt(iwrkd+3) = irng2
                    else
                        irngt(iwrkd+3) = irngt(iwrkd+4)
                        irngt(iwrkd+4) = irng2
                    end if
                else
                    irng1 = irngt(iwrkd+1)
                    irng2 = irngt(iwrkd+2)
                    irngt(iwrkd+1) = irngt(iwrkd+3)
                    if(xdont(irng1) <= xdont(irngt(iwrkd+4))) then
                        irngt(iwrkd+2) = irng1
                        if(xdont(irng2) <= xdont(irngt(iwrkd+4))) then
                            irngt(iwrkd+3) = irng2
                        else
                            irngt(iwrkd+3) = irngt(iwrkd+4)
                            irngt(iwrkd+4) = irng2
                        end if
                    else
                        irngt(iwrkd+2) = irngt(iwrkd+4)
                        irngt(iwrkd+3) = irng1
                        irngt(iwrkd+4) = irng2
                    end if
                end if
            end do
            lmtna = 4
            exit
        end do
        do
            if(lmtna >= nval) exit
            iwrkf = 0
            lmtnc = 2 * lmtnc
            do
                iwrk = iwrkf
                iwrkd = iwrkf + 1
                jinda = iwrkf + lmtna
                iwrkf = iwrkf + lmtnc
                if(iwrkf >= nval) then
                    if(jinda >= nval) exit
                    iwrkf = nval
                end if
                iinda = 1
                iindb = jinda + 1
                jwrkt(1:lmtna) = irngt(iwrkd:jinda)
                xvala = xdont(jwrkt(iinda))
                xvalb = xdont(irngt(iindb))
                do
                    iwrk = iwrk + 1
                    if(xvala > xvalb) then
                        irngt(iwrk) = irngt(iindb)
                        iindb = iindb + 1
                        if(iindb > iwrkf) then
                            irngt(iwrk+1:iwrkf) = jwrkt(iinda:lmtna)
                            exit
                        end if
                        xvalb = xdont(irngt(iindb))
                    else
                        irngt(iwrk) = jwrkt(iinda)
                        iinda = iinda + 1
                        if(iinda > lmtna) exit! only b still with unprocessed values
                        xvala = xdont(jwrkt(iinda))
                    end if
                end do
            end do
            lmtna = 2 * lmtna
        end do
        return
        end subroutine d_rank
        subroutine r_rank(xdont, irngt)
        real, dimension(:), intent(in) :: xdont
        integer, dimension(:), intent(out) :: irngt
        real :: xvala, xvalb
        integer, dimension(size(irngt)) :: jwrkt
        integer :: lmtna, lmtnc, irng1, irng2
        integer :: nval, iind, iwrkd, iwrk, iwrkf, jinda, iinda, iindb
        nval = min(size(xdont), size(irngt))
        select case(nval)
        case(:0)
            return
        case(1)
            irngt(1) = 1
            return
            case default
            continue
        end select
        do iind = 2, nval, 2
            if(xdont(iind-1) <= xdont(iind)) then
                irngt(iind-1) = iind - 1
                irngt(iind) = iind
            else
                irngt(iind-1) = iind
                irngt(iind) = iind - 1
            end if
        end do
        if(modulo(nval, 2) /= 0) then
            irngt(nval) = nval
        end if
        lmtna = 2
        lmtnc = 4
        do
            if(nval <= 2) exit
            do iwrkd = 0, nval - 1, 4
                if((iwrkd+4) > nval) then
                    if((iwrkd+2) >= nval) exit
                    if(xdont(irngt(iwrkd+2)) <= xdont(irngt(iwrkd+3))) exit
                    if(xdont(irngt(iwrkd+1)) <= xdont(irngt(iwrkd+3))) then
                        irng2 = irngt(iwrkd+2)
                        irngt(iwrkd+2) = irngt(iwrkd+3)
                        irngt(iwrkd+3) = irng2
                    else
                        irng1 = irngt(iwrkd+1)
                        irngt(iwrkd+1) = irngt(iwrkd+3)
                        irngt(iwrkd+3) = irngt(iwrkd+2)
                        irngt(iwrkd+2) = irng1
                    end if
                    exit
                end if
                if(xdont(irngt(iwrkd+2)) <= xdont(irngt(iwrkd+3))) cycle
                if(xdont(irngt(iwrkd+1)) <= xdont(irngt(iwrkd+3))) then
                    irng2 = irngt(iwrkd+2)
                    irngt(iwrkd+2) = irngt(iwrkd+3)
                    if(xdont(irng2) <= xdont(irngt(iwrkd+4))) then
                        irngt(iwrkd+3) = irng2
                    else
                        irngt(iwrkd+3) = irngt(iwrkd+4)
                        irngt(iwrkd+4) = irng2
                    end if
                else
                    irng1 = irngt(iwrkd+1)
                    irng2 = irngt(iwrkd+2)
                    irngt(iwrkd+1) = irngt(iwrkd+3)
                    if(xdont(irng1) <= xdont(irngt(iwrkd+4))) then
                        irngt(iwrkd+2) = irng1
                        if(xdont(irng2) <= xdont(irngt(iwrkd+4))) then
                            irngt(iwrkd+3) = irng2
                        else
                            irngt(iwrkd+3) = irngt(iwrkd+4)
                            irngt(iwrkd+4) = irng2
                        end if
                    else
                        irngt(iwrkd+2) = irngt(iwrkd+4)
                        irngt(iwrkd+3) = irng1
                        irngt(iwrkd+4) = irng2
                    end if
                end if
            end do
            lmtna = 4
            exit
        end do
        do
            if(lmtna >= nval) exit
            iwrkf = 0
            lmtnc = 2 * lmtnc
            do
                iwrk = iwrkf
                iwrkd = iwrkf + 1
                jinda = iwrkf + lmtna
                iwrkf = iwrkf + lmtnc
                if(iwrkf >= nval) then
                    if(jinda >= nval) exit
                    iwrkf = nval
                end if
                iinda = 1
                iindb = jinda + 1
                jwrkt(1:lmtna) = irngt(iwrkd:jinda)
                xvala = xdont(jwrkt(iinda))
                xvalb = xdont(irngt(iindb))
                do
                    iwrk = iwrk + 1
                    if(xvala > xvalb) then
                        irngt(iwrk) = irngt(iindb)
                        iindb = iindb + 1
                        if(iindb > iwrkf) then
                            irngt(iwrk+1:iwrkf) = jwrkt(iinda:lmtna)
                            exit
                        end if
                        xvalb = xdont(irngt(iindb))
                    else
                        irngt(iwrk) = jwrkt(iinda)
                        iinda = iinda + 1
                        if(iinda > lmtna) exit! only b still with unprocessed values
                        xvala = xdont(jwrkt(iinda))
                    end if
                end do
            end do
            lmtna = 2 * lmtna
        end do
        return
        end subroutine r_rank
        subroutine i_rank(xdont, irngt)
        integer, dimension(:), intent(in)  :: xdont
        integer, dimension(:), intent(out) :: irngt
        integer :: xvala, xvalb
        integer, dimension(size(irngt)) :: jwrkt
        integer :: lmtna, lmtnc, irng1, irng2
        integer :: nval, iind, iwrkd, iwrk, iwrkf, jinda, iinda, iindb
        nval = min(size(xdont), size(irngt))
        select case(nval)
        case(:0)
            return
        case(1)
            irngt(1) = 1
            return
            case default
            continue
        end select
        do iind = 2, nval, 2
            if(xdont(iind-1) <= xdont(iind)) then
                irngt(iind-1) = iind - 1
                irngt(iind) = iind
            else
                irngt(iind-1) = iind
                irngt(iind) = iind - 1
            end if
        end do
        if(modulo(nval, 2) /= 0) then
            irngt(nval) = nval
        end if
        lmtna = 2
        lmtnc = 4
        do
            if(nval <= 2) exit
            do iwrkd = 0, nval - 1, 4
                if((iwrkd+4) > nval) then
                    if((iwrkd+2) >= nval) exit
                    if(xdont(irngt(iwrkd+2)) <= xdont(irngt(iwrkd+3))) exit
                    if(xdont(irngt(iwrkd+1)) <= xdont(irngt(iwrkd+3))) then
                        irng2 = irngt(iwrkd+2)
                        irngt(iwrkd+2) = irngt(iwrkd+3)
                        irngt(iwrkd+3) = irng2
                    else
                        irng1 = irngt(iwrkd+1)
                        irngt(iwrkd+1) = irngt(iwrkd+3)
                        irngt(iwrkd+3) = irngt(iwrkd+2)
                        irngt(iwrkd+2) = irng1
                    end if
                    exit
                end if
                if(xdont(irngt(iwrkd+2)) <= xdont(irngt(iwrkd+3))) cycle
                if(xdont(irngt(iwrkd+1)) <= xdont(irngt(iwrkd+3))) then
                    irng2 = irngt(iwrkd+2)
                    irngt(iwrkd+2) = irngt(iwrkd+3)
                    if(xdont(irng2) <= xdont(irngt(iwrkd+4))) then
                        irngt(iwrkd+3) = irng2
                    else
                        irngt(iwrkd+3) = irngt(iwrkd+4)
                        irngt(iwrkd+4) = irng2
                    end if
                else
                    irng1 = irngt(iwrkd+1)
                    irng2 = irngt(iwrkd+2)
                    irngt(iwrkd+1) = irngt(iwrkd+3)
                    if(xdont(irng1) <= xdont(irngt(iwrkd+4))) then
                        irngt(iwrkd+2) = irng1
                        if(xdont(irng2) <= xdont(irngt(iwrkd+4))) then
                            irngt(iwrkd+3) = irng2
                        else
                            irngt(iwrkd+3) = irngt(iwrkd+4)
                            irngt(iwrkd+4) = irng2
                        end if
                    else
                        irngt(iwrkd+2) = irngt(iwrkd+4)
                        irngt(iwrkd+3) = irng1
                        irngt(iwrkd+4) = irng2
                    end if
                end if
            end do
            lmtna = 4
            exit
        end do
        do
            if(lmtna >= nval) exit
            iwrkf = 0
            lmtnc = 2 * lmtnc
            do
                iwrk = iwrkf
                iwrkd = iwrkf + 1
                jinda = iwrkf + lmtna
                iwrkf = iwrkf + lmtnc
                if(iwrkf >= nval) then
                    if(jinda >= nval) exit
                    iwrkf = nval
                end if
                iinda = 1
                iindb = jinda + 1
                jwrkt(1:lmtna) = irngt(iwrkd:jinda)
                xvala = xdont(jwrkt(iinda))
                xvalb = xdont(irngt(iindb))
                do
                    iwrk = iwrk + 1
                    if(xvala > xvalb) then
                        irngt(iwrk) = irngt(iindb)
                        iindb = iindb + 1
                        if(iindb > iwrkf) then
                            irngt(iwrk+1:iwrkf) = jwrkt(iinda:lmtna)
                            exit
                        end if
                        xvalb = xdont(irngt(iindb))
                    else
                        irngt(iwrk) = jwrkt(iinda)
                        iinda = iinda + 1
                        if(iinda > lmtna) exit! only b still with unprocessed values
                        xvala = xdont(jwrkt(iinda))
                    end if
                end do
            end do
            lmtna = 2 * lmtna
        end do
        return
    end subroutine i_rank
    subroutine d_unique_rank(xvalt, irngt, nuni)
        real(kind=dp), dimension(:), intent(in) :: xvalt
        integer, dimension(:), intent(out) :: irngt
        integer, intent(out) :: nuni
        integer, dimension(size(irngt)) :: jwrkt
        integer :: lmtna, lmtnc, irng, irng1, irng2
        integer :: nval, iind, iwrkd, iwrk, iwrkf, jinda, iinda, iindb
        real(kind=dp) :: xtst, xvala, xvalb
        nval = min(size(xvalt), size(irngt))
        nuni = nval
        select case(nval)
        case(:0)
            return
        case(1)
            irngt(1) = 1
            return
            case default
            continue
        end select
        do iind = 2, nval, 2
            if(xvalt(iind-1) < xvalt(iind)) then
                irngt(iind-1) = iind - 1
                irngt(iind) = iind
            else
                irngt(iind-1) = iind
                irngt(iind) = iind - 1
            end if
        end do
        if(modulo(nval, 2) /= 0) then
            irngt(nval) = nval
        end if
        lmtna = 2
        lmtnc = 4
        do
            if(nval <= 4) exit
            do iwrkd = 0, nval - 1, 4
                if((iwrkd+4) > nval) then
                    if((iwrkd+2) >= nval) exit
                    if(xvalt(irngt(iwrkd+2)) <= xvalt(irngt(iwrkd+3))) exit
                    if(xvalt(irngt(iwrkd+1)) <= xvalt(irngt(iwrkd+3))) then
                        irng2 = irngt(iwrkd+2)
                        irngt(iwrkd+2) = irngt(iwrkd+3)
                        irngt(iwrkd+3) = irng2
                    else
                        irng1 = irngt(iwrkd+1)
                        irngt(iwrkd+1) = irngt(iwrkd+3)
                        irngt(iwrkd+3) = irngt(iwrkd+2)
                        irngt(iwrkd+2) = irng1
                    end if
                    exit
                end if
                if(xvalt(irngt(iwrkd+2)) <= xvalt(irngt(iwrkd+3))) cycle
                if(xvalt(irngt(iwrkd+1)) <= xvalt(irngt(iwrkd+3))) then
                    irng2 = irngt(iwrkd+2)
                    irngt(iwrkd+2) = irngt(iwrkd+3)
                    if(xvalt(irng2) <= xvalt(irngt(iwrkd+4))) then
                        irngt(iwrkd+3) = irng2
                    else
                        irngt(iwrkd+3) = irngt(iwrkd+4)
                        irngt(iwrkd+4) = irng2
                    end if
                else
                    irng1 = irngt(iwrkd+1)
                    irng2 = irngt(iwrkd+2)
                    irngt(iwrkd+1) = irngt(iwrkd+3)
                    if(xvalt(irng1) <= xvalt(irngt(iwrkd+4))) then
                        irngt(iwrkd+2) = irng1
                        if(xvalt(irng2) <= xvalt(irngt(iwrkd+4))) then
                            irngt(iwrkd+3) = irng2
                        else
                            irngt(iwrkd+3) = irngt(iwrkd+4)
                            irngt(iwrkd+4) = irng2
                        end if
                    else
                        irngt(iwrkd+2) = irngt(iwrkd+4)
                        irngt(iwrkd+3) = irng1
                        irngt(iwrkd+4) = irng2
                    end if
                end if
            end do
            lmtna = 4
            exit
        end do
        do
            if(2*lmtna >= nval) exit
            iwrkf = 0
            lmtnc = 2 * lmtnc
            do
                iwrk = iwrkf
                iwrkd = iwrkf + 1
                jinda = iwrkf + lmtna
                iwrkf = iwrkf + lmtnc
                if(iwrkf >= nval) then
                    if(jinda >= nval) exit
                    iwrkf = nval
                end if
                iinda = 1
                iindb = jinda + 1
                jwrkt(1:lmtna) = irngt(iwrkd:jinda)
                xvala = xvalt(jwrkt(iinda))
                xvalb = xvalt(irngt(iindb))
                do
                    iwrk = iwrk + 1
                    if(xvala > xvalb) then
                        irngt(iwrk) = irngt(iindb)
                        iindb = iindb + 1
                        if(iindb > iwrkf) then
                            irngt(iwrk+1:iwrkf) = jwrkt(iinda:lmtna)
                            exit
                        end if
                        xvalb = xvalt(irngt(iindb))
                    else
                        irngt(iwrk) = jwrkt(iinda)
                        iinda = iinda + 1
                        if(iinda > lmtna) exit! only b still with unprocessed values
                        xvala = xvalt(jwrkt(iinda))
                    end if
                end do
            end do
            lmtna = 2 * lmtna
        end do
        iinda = 1
        iindb = lmtna + 1
        nuni = 0
        jwrkt(1:lmtna) = irngt(1:lmtna)
        if(iindb <= nval) then
            xtst = nearless(min(xvalt(jwrkt(1)), xvalt(irngt(iindb))))
        else
            xtst = nearless(xvalt(jwrkt(1)))
        endif
        do iwrk = 1, nval
            if(iinda <= lmtna) then
                if(iindb <= nval) then
                    if(xvalt(jwrkt(iinda)) > xvalt(irngt(iindb))) then
                        irng = irngt(iindb)
                        iindb = iindb + 1
                    else
                        irng = jwrkt(iinda)
                        iinda = iinda + 1
                    end if
                else
                    irng = jwrkt(iinda)
                    iinda = iinda + 1
                end if
            else
                irng = irngt(iwrk)
            end if
            if(xvalt(irng) > xtst) then
                xtst = xvalt(irng)
                nuni = nuni + 1
                irngt(nuni) = irng
            end if
        end do
        return
        end subroutine d_unique_rank
        subroutine r_unique_rank(xvalt, irngt, nuni)
        real, dimension(:), intent(in) :: xvalt
        integer, dimension(:), intent(out) :: irngt
        integer, intent(out) :: nuni
        integer, dimension(size(irngt)) :: jwrkt
        integer :: lmtna, lmtnc, irng, irng1, irng2
        integer :: nval, iind, iwrkd, iwrk, iwrkf, jinda, iinda, iindb
        real :: xtst, xvala, xvalb
        nval = min(size(xvalt), size(irngt))
        nuni = nval
        select case(nval)
        case(:0)
            return
        case(1)
            irngt(1) = 1
            return
            case default
            continue
        end select
        do iind = 2, nval, 2
            if(xvalt(iind-1) < xvalt(iind)) then
                irngt(iind-1) = iind - 1
                irngt(iind) = iind
            else
                irngt(iind-1) = iind
                irngt(iind) = iind - 1
            end if
        end do
        if(modulo(nval, 2) /= 0) then
            irngt(nval) = nval
        end if
        lmtna = 2
        lmtnc = 4
        do
            if(nval <= 4) exit
            do iwrkd = 0, nval - 1, 4
                if((iwrkd+4) > nval) then
                    if((iwrkd+2) >= nval) exit
                    if(xvalt(irngt(iwrkd+2)) <= xvalt(irngt(iwrkd+3))) exit
                    if(xvalt(irngt(iwrkd+1)) <= xvalt(irngt(iwrkd+3))) then
                        irng2 = irngt(iwrkd+2)
                        irngt(iwrkd+2) = irngt(iwrkd+3)
                        irngt(iwrkd+3) = irng2
                    else
                        irng1 = irngt(iwrkd+1)
                        irngt(iwrkd+1) = irngt(iwrkd+3)
                        irngt(iwrkd+3) = irngt(iwrkd+2)
                        irngt(iwrkd+2) = irng1
                    end if
                    exit
                end if
                if(xvalt(irngt(iwrkd+2)) <= xvalt(irngt(iwrkd+3))) cycle
                if(xvalt(irngt(iwrkd+1)) <= xvalt(irngt(iwrkd+3))) then
                    irng2 = irngt(iwrkd+2)
                    irngt(iwrkd+2) = irngt(iwrkd+3)
                    if(xvalt(irng2) <= xvalt(irngt(iwrkd+4))) then
                        irngt(iwrkd+3) = irng2
                    else
                        irngt(iwrkd+3) = irngt(iwrkd+4)
                        irngt(iwrkd+4) = irng2
                    end if
                else
                    irng1 = irngt(iwrkd+1)
                    irng2 = irngt(iwrkd+2)
                    irngt(iwrkd+1) = irngt(iwrkd+3)
                    if(xvalt(irng1) <= xvalt(irngt(iwrkd+4))) then
                        irngt(iwrkd+2) = irng1
                        if(xvalt(irng2) <= xvalt(irngt(iwrkd+4))) then
                            irngt(iwrkd+3) = irng2
                        else
                            irngt(iwrkd+3) = irngt(iwrkd+4)
                            irngt(iwrkd+4) = irng2
                        end if
                    else
                        irngt(iwrkd+2) = irngt(iwrkd+4)
                        irngt(iwrkd+3) = irng1
                        irngt(iwrkd+4) = irng2
                    end if
                end if
            end do
            lmtna = 4
            exit
        end do
        do
            if(2*lmtna >= nval) exit
            iwrkf = 0
            lmtnc = 2 * lmtnc
            do
                iwrk = iwrkf
                iwrkd = iwrkf + 1
                jinda = iwrkf + lmtna
                iwrkf = iwrkf + lmtnc
                if(iwrkf >= nval) then
                    if(jinda >= nval) exit
                    iwrkf = nval
                end if
                iinda = 1
                iindb = jinda + 1
                jwrkt(1:lmtna) = irngt(iwrkd:jinda)
                xvala = xvalt(jwrkt(iinda))
                xvalb = xvalt(irngt(iindb))
                do
                    iwrk = iwrk + 1
                    if(xvala > xvalb) then
                        irngt(iwrk) = irngt(iindb)
                        iindb = iindb + 1
                        if(iindb > iwrkf) then
                            irngt(iwrk+1:iwrkf) = jwrkt(iinda:lmtna)
                            exit
                        end if
                        xvalb = xvalt(irngt(iindb))
                    else
                        irngt(iwrk) = jwrkt(iinda)
                        iinda = iinda + 1
                        if(iinda > lmtna) exit! only b still with unprocessed values
                        xvala = xvalt(jwrkt(iinda))
                    end if
                end do
            end do
            lmtna = 2 * lmtna
        end do
        iinda = 1
        iindb = lmtna + 1
        nuni = 0
        jwrkt(1:lmtna) = irngt(1:lmtna)
        if(iindb <= nval) then
            xtst = nearless(min(xvalt(jwrkt(1)), xvalt(irngt(iindb))))
        else
            xtst = nearless(xvalt(jwrkt(1)))
        endif
        do iwrk = 1, nval
            if(iinda <= lmtna) then
                if(iindb <= nval) then
                    if(xvalt(jwrkt(iinda)) > xvalt(irngt(iindb))) then
                        irng = irngt(iindb)
                        iindb = iindb + 1
                    else
                        irng = jwrkt(iinda)
                        iinda = iinda + 1
                    end if
                else
                    irng = jwrkt(iinda)
                    iinda = iinda + 1
                end if
            else
                irng = irngt(iwrk)
            end if
            if(xvalt(irng) > xtst) then
                xtst = xvalt(irng)
                nuni = nuni + 1
                irngt(nuni) = irng
            end if
        end do
        return
    end subroutine r_unique_rank
    subroutine i_unique_rank(xvalt, irngt, nuni)
        integer, dimension(:), intent(in) :: xvalt
        integer, dimension(:), intent(out) :: irngt
        integer, intent(out) :: nuni
        integer, dimension(size(irngt)) :: jwrkt
        integer :: lmtna, lmtnc, irng, irng1, irng2
        integer :: nval, iind, iwrkd, iwrk, iwrkf, jinda, iinda, iindb
        integer :: xtst, xvala, xvalb
        nval = min(size(xvalt), size(irngt))
        nuni = nval
        select case(nval)
        case(:0)
            return
        case(1)
            irngt(1) = 1
            return
            case default
            continue
        end select
        do iind = 2, nval, 2
            if(xvalt(iind-1) < xvalt(iind)) then
                irngt(iind-1) = iind - 1
                irngt(iind) = iind
            else
                irngt(iind-1) = iind
                irngt(iind) = iind - 1
            end if
        end do
        if(modulo(nval, 2) /= 0) then
            irngt(nval) = nval
        end if
        lmtna = 2
        lmtnc = 4
        do
            if(nval <= 4) exit
            do iwrkd = 0, nval - 1, 4
                if((iwrkd+4) > nval) then
                    if((iwrkd+2) >= nval) exit
                    if(xvalt(irngt(iwrkd+2)) <= xvalt(irngt(iwrkd+3))) exit
                    if(xvalt(irngt(iwrkd+1)) <= xvalt(irngt(iwrkd+3))) then
                        irng2 = irngt(iwrkd+2)
                        irngt(iwrkd+2) = irngt(iwrkd+3)
                        irngt(iwrkd+3) = irng2
                    else
                        irng1 = irngt(iwrkd+1)
                        irngt(iwrkd+1) = irngt(iwrkd+3)
                        irngt(iwrkd+3) = irngt(iwrkd+2)
                        irngt(iwrkd+2) = irng1
                    end if
                    exit
                end if
                if(xvalt(irngt(iwrkd+2)) <= xvalt(irngt(iwrkd+3))) cycle
                if(xvalt(irngt(iwrkd+1)) <= xvalt(irngt(iwrkd+3))) then
                    irng2 = irngt(iwrkd+2)
                    irngt(iwrkd+2) = irngt(iwrkd+3)
                    if(xvalt(irng2) <= xvalt(irngt(iwrkd+4))) then
                        irngt(iwrkd+3) = irng2
                    else
                        irngt(iwrkd+3) = irngt(iwrkd+4)
                        irngt(iwrkd+4) = irng2
                    end if
                else
                    irng1 = irngt(iwrkd+1)
                    irng2 = irngt(iwrkd+2)
                    irngt(iwrkd+1) = irngt(iwrkd+3)
                    if(xvalt(irng1) <= xvalt(irngt(iwrkd+4))) then
                        irngt(iwrkd+2) = irng1
                        if(xvalt(irng2) <= xvalt(irngt(iwrkd+4))) then
                            irngt(iwrkd+3) = irng2
                        else
                            irngt(iwrkd+3) = irngt(iwrkd+4)
                            irngt(iwrkd+4) = irng2
                        end if
                    else
                        irngt(iwrkd+2) = irngt(iwrkd+4)
                        irngt(iwrkd+3) = irng1
                        irngt(iwrkd+4) = irng2
                    end if
                end if
            end do
            lmtna = 4
            exit
        end do
        do
            if(2*lmtna >= nval) exit
            iwrkf = 0
            lmtnc = 2 * lmtnc
            do
                iwrk = iwrkf
                iwrkd = iwrkf + 1
                jinda = iwrkf + lmtna
                iwrkf = iwrkf + lmtnc
                if(iwrkf >= nval) then
                    if(jinda >= nval) exit
                    iwrkf = nval
                end if
                iinda = 1
                iindb = jinda + 1
                jwrkt(1:lmtna) = irngt(iwrkd:jinda)
                xvala = xvalt(jwrkt(iinda))
                xvalb = xvalt(irngt(iindb))
                do
                    iwrk = iwrk + 1
                    if(xvala > xvalb) then
                        irngt(iwrk) = irngt(iindb)
                        iindb = iindb + 1
                        if(iindb > iwrkf) then
                            irngt(iwrk+1:iwrkf) = jwrkt(iinda:lmtna)
                            exit
                        end if
                        xvalb = xvalt(irngt(iindb))
                    else
                        irngt(iwrk) = jwrkt(iinda)
                        iinda = iinda + 1
                        if(iinda > lmtna) exit! only b still with unprocessed values
                        xvala = xvalt(jwrkt(iinda))
                    end if
                end do
            end do
            lmtna = 2 * lmtna
        end do
        iinda = 1
        iindb = lmtna + 1
        nuni = 0
        jwrkt(1:lmtna) = irngt(1:lmtna)
        if(iindb <= nval) then
            xtst = nearless(min(xvalt(jwrkt(1)), xvalt(irngt(iindb))))
        else
            xtst = nearless(xvalt(jwrkt(1)))
        endif
        do iwrk = 1, nval
            if(iinda <= lmtna) then
                if(iindb <= nval) then
                    if(xvalt(jwrkt(iinda)) > xvalt(irngt(iindb))) then
                        irng = irngt(iindb)
                        iindb = iindb + 1
                    else
                        irng = jwrkt(iinda)
                        iinda = iinda + 1
                    end if
                else
                    irng = jwrkt(iinda)
                    iinda = iinda + 1
                end if
            else
                irng = irngt(iwrk)
            end if
            if(xvalt(irng) > xtst) then
                xtst = xvalt(irng)
                nuni = nuni + 1
                irngt(nuni) = irng
            end if
        end do
        return
        end subroutine i_unique_rank

        function d_nearless(xval) result(d_nl)
        real(kind=dp), intent(in) :: xval
        real(kind=dp) :: d_nl
        d_nl = nearest(xval, -1.0_dp)
        return
        end function d_nearless
        function r_nearless(xval) result(r_nl)
        real, intent(in) :: xval
        real :: r_nl
        r_nl = nearest(xval, -1.0)
        return
        end function r_nearless
        function i_nearless(xval) result(i_nl)
        integer, intent(in) :: xval
        integer :: i_nl
        i_nl = xval - 1
        return
    end function i_nearless

    subroutine d_sort_using_quicksort(xdont)
        real(kind=dp), dimension(:), intent(inout) :: xdont
        call d_subsor(xdont, 1, size(xdont))
        call d_sort_using_insertion_for_sort_using_quicksort(xdont)
        return
    end subroutine d_sort_using_quicksort
    recursive subroutine d_subsor(xdont, ideb1, ifin1)
        real(kind=dp), dimension(:), intent(inout) :: xdont
        integer, intent(in) :: ideb1, ifin1
        integer, parameter :: nins = 16 ! max for insertion sort
        integer :: icrs, ideb, idcr, ifin, imil
        real(kind=dp) :: xpiv, xwrk
        ideb = ideb1
        ifin = ifin1
        if((ifin - ideb) > nins) then
            imil =(ideb+ifin) / 2
            if(xdont(imil) < xdont(ideb)) then
                xwrk = xdont(ideb)
                xdont(ideb) = xdont(imil)
                xdont(imil) = xwrk
            end if
            if(xdont(imil) > xdont(ifin)) then
                xwrk = xdont(ifin)
                xdont(ifin) = xdont(imil)
                xdont(imil) = xwrk
                if(xdont(imil) < xdont(ideb)) then
                    xwrk = xdont(ideb)
                    xdont(ideb) = xdont(imil)
                    xdont(imil) = xwrk
                end if
            end if
            xpiv = xdont(imil)
            icrs = ideb
            idcr = ifin
            ech2: do
                do
                    icrs = icrs + 1
                    if(icrs >= idcr) then
                        exit ech2
                    end if
                    if(xdont(icrs) > xpiv) exit
                end do
                do
                    if(xdont(idcr) <= xpiv) exit
                    idcr = idcr - 1
                    if(icrs >= idcr) then
                        exit ech2
                    end if
                end do
                xwrk = xdont(idcr)
                xdont(idcr) = xdont(icrs)
                xdont(icrs) = xwrk
            end do ech2
            call d_subsor(xdont, ideb1, icrs-1)
            call d_subsor(xdont, idcr, ifin1)
        end if
        return
    end subroutine d_subsor
    subroutine d_sort_using_insertion_for_sort_using_quicksort(xdont)
        real(kind=dp), dimension(:), intent(inout) :: xdont
        integer :: icrs, idcr
        real(kind=dp) :: xwrk
        do icrs = 2, size(xdont)
            xwrk = xdont(icrs)
            if(xwrk >= xdont(icrs-1)) cycle
            xdont(icrs) = xdont(icrs-1)
            do idcr = icrs - 2, 1, - 1
                if(xwrk >= xdont(idcr)) exit
                xdont(idcr+1) = xdont(idcr)
            end do
            xdont(idcr+1) = xwrk
        end do
        return
        end subroutine d_sort_using_insertion_for_sort_using_quicksort
        subroutine r_sort_using_quicksort(xdont)
        real, dimension(:), intent(inout) :: xdont
        call r_subsor(xdont, 1, size(xdont))
        call r_sort_using_insertion_for_sort_using_quicksort(xdont)
        return
    end subroutine r_sort_using_quicksort
    recursive subroutine r_subsor(xdont, ideb1, ifin1)
        real, dimension(:), intent(inout) :: xdont
        integer, intent(in) :: ideb1, ifin1
        integer, parameter :: nins = 16 ! max for insertion sort
        integer :: icrs, ideb, idcr, ifin, imil
        real :: xpiv, xwrk
        ideb = ideb1
        ifin = ifin1
        if((ifin - ideb) > nins) then
            imil =(ideb+ifin) / 2
            if(xdont(imil) < xdont(ideb)) then
                xwrk = xdont(ideb)
                xdont(ideb) = xdont(imil)
                xdont(imil) = xwrk
            end if
            if(xdont(imil) > xdont(ifin)) then
                xwrk = xdont(ifin)
                xdont(ifin) = xdont(imil)
                xdont(imil) = xwrk
                if(xdont(imil) < xdont(ideb)) then
                    xwrk = xdont(ideb)
                    xdont(ideb) = xdont(imil)
                    xdont(imil) = xwrk
                end if
            end if
            xpiv = xdont(imil)
            icrs = ideb
            idcr = ifin
            ech2: do
                do
                    icrs = icrs + 1
                    if(icrs >= idcr) then
                        exit ech2
                    end if
                    if(xdont(icrs) > xpiv) exit
                end do
                do
                    if(xdont(idcr) <= xpiv) exit
                    idcr = idcr - 1
                    if(icrs >= idcr) then
                        exit ech2
                    end if
                end do
                xwrk = xdont(idcr)
                xdont(idcr) = xdont(icrs)
                xdont(icrs) = xwrk
            end do ech2
            call r_subsor(xdont, ideb1, icrs-1)
            call r_subsor(xdont, idcr, ifin1)
        end if
        return
    end subroutine r_subsor
    subroutine r_sort_using_insertion_for_sort_using_quicksort(xdont)
        real, dimension(:), intent(inout) :: xdont
        integer :: icrs, idcr
        real :: xwrk
        do icrs = 2, size(xdont)
            xwrk = xdont(icrs)
            if(xwrk >= xdont(icrs-1)) cycle
            xdont(icrs) = xdont(icrs-1)
            do idcr = icrs - 2, 1, - 1
                if(xwrk >= xdont(idcr)) exit
                xdont(idcr+1) = xdont(idcr)
            end do
            xdont(idcr+1) = xwrk
        end do
        return
    end subroutine r_sort_using_insertion_for_sort_using_quicksort
    subroutine i_sort_using_quicksort(xdont)
        integer, dimension(:), intent(inout)  :: xdont
        call i_subsor(xdont, 1, size(xdont))
        call i_sort_using_insertion_for_sort_using_quicksort(xdont)
        return
    end subroutine i_sort_using_quicksort
    recursive subroutine i_subsor(xdont, ideb1, ifin1)
        integer, dimension(:), intent(inout) :: xdont
        integer, intent(in) :: ideb1, ifin1
        integer, parameter :: nins = 16 ! max for insertion sort
        integer :: icrs, ideb, idcr, ifin, imil
        integer :: xpiv, xwrk
        ideb = ideb1
        ifin = ifin1
        if((ifin - ideb) > nins) then
            imil =(ideb+ifin) / 2
            if(xdont(imil) < xdont(ideb)) then
                xwrk = xdont(ideb)
                xdont(ideb) = xdont(imil)
                xdont(imil) = xwrk
            end if
            if(xdont(imil) > xdont(ifin)) then
                xwrk = xdont(ifin)
                xdont(ifin) = xdont(imil)
                xdont(imil) = xwrk
                if(xdont(imil) < xdont(ideb)) then
                    xwrk = xdont(ideb)
                    xdont(ideb) = xdont(imil)
                    xdont(imil) = xwrk
                end if
            end if
            xpiv = xdont(imil)
            icrs = ideb
            idcr = ifin
            ech2: do
                do
                    icrs = icrs + 1
                    if(icrs >= idcr) then
                        exit ech2
                    end if
                    if(xdont(icrs) > xpiv) exit
                end do
                do
                    if(xdont(idcr) <= xpiv) exit
                    idcr = idcr - 1
                    if(icrs >= idcr) then
                        exit ech2
                    end if
                end do
                xwrk = xdont(idcr)
                xdont(idcr) = xdont(icrs)
                xdont(icrs) = xwrk
            end do ech2
            call i_subsor(xdont, ideb1, icrs-1)
            call i_subsor(xdont, idcr, ifin1)
        end if
        return
    end subroutine i_subsor
    subroutine i_sort_using_insertion_for_sort_using_quicksort(xdont)
        integer, dimension(:), intent(inout) :: xdont
        integer :: icrs, idcr
        integer :: xwrk
        do icrs = 2, size(xdont)
            xwrk = xdont(icrs)
            if(xwrk >= xdont(icrs-1)) cycle
            xdont(icrs) = xdont(icrs-1)
            do idcr = icrs - 2, 1, - 1
                if(xwrk >= xdont(idcr)) exit
                xdont(idcr+1) = xdont(idcr)
            end do
            xdont(idcr+1) = xwrk
        end do
        return
    end subroutine i_sort_using_insertion_for_sort_using_quicksort

    subroutine d_sort_using_insertion(xdont)
        real(kind=dp), dimension(:), intent(inout) :: xdont
        real(kind=dp) :: xwrk, xmin
        integer :: icrs, idcr, ndon
        ndon = size(xdont)
        if(xdont(1) < xdont(ndon)) then
            xmin = xdont(1)
        else
            xmin = xdont(ndon)
            xdont(ndon) = xdont(1)
        endif
        do idcr = ndon-1, 2, -1
            xwrk = xdont(idcr)
            if(xwrk < xmin) then
                xdont(idcr) = xmin
                xmin = xwrk
            end if
        end do
        xdont(1) = xmin
        do icrs = 3, ndon
            xwrk = xdont(icrs)
            idcr = icrs - 1
            if(xwrk < xdont(idcr)) then
                xdont(icrs) = xdont(idcr)
                idcr = idcr - 1
                do
                    if(xwrk >= xdont(idcr)) exit
                    xdont(idcr+1) = xdont(idcr)
                    idcr = idcr - 1
                end do
                xdont(idcr+1) = xwrk
            end if
        end do
        return
    end subroutine d_sort_using_insertion
    subroutine r_sort_using_insertion(xdont)
        real, dimension(:), intent(inout) :: xdont
        real :: xwrk, xmin
        integer :: icrs, idcr, ndon
        ndon = size(xdont)
        if(xdont(1) < xdont(ndon)) then
            xmin = xdont(1)
        else
            xmin = xdont(ndon)
            xdont(ndon) = xdont(1)
        endif
        do idcr = ndon-1, 2, -1
            xwrk = xdont(idcr)
            if(xwrk < xmin) then
                xdont(idcr) = xmin
                xmin = xwrk
            end if
        end do
        xdont(1) = xmin
        do icrs = 3, ndon
            xwrk = xdont(icrs)
            idcr = icrs - 1
            if(xwrk < xdont(idcr)) then
                xdont(icrs) = xdont(idcr)
                idcr = idcr - 1
                do
                    if(xwrk >= xdont(idcr)) exit
                    xdont(idcr+1) = xdont(idcr)
                    idcr = idcr - 1
                end do
                xdont(idcr+1) = xwrk
            end if
        end do
        return
    end subroutine r_sort_using_insertion
    subroutine i_sort_using_insertion(xdont)
        integer, dimension(:), intent(inout)  :: xdont
        integer :: xwrk, xmin
        integer :: icrs, idcr, ndon
        ndon = size(xdont)
        if(xdont(1) < xdont(ndon)) then
            xmin = xdont(1)
        else
            xmin = xdont(ndon)
            xdont(ndon) = xdont(1)
        endif
        do idcr = ndon-1, 2, -1
            xwrk = xdont(idcr)
            if(xwrk < xmin) then
                xdont(idcr) = xmin
                xmin = xwrk
            end if
        end do
        xdont(1) = xmin
        do icrs = 3, ndon
            xwrk = xdont(icrs)
            idcr = icrs - 1
            if(xwrk < xdont(idcr)) then
                xdont(icrs) = xdont(idcr)
                idcr = idcr - 1
                do
                    if(xwrk >= xdont(idcr)) exit
                    xdont(idcr+1) = xdont(idcr)
                    idcr = idcr - 1
                end do
                xdont(idcr+1) = xwrk
            end if
        end do
        return
    end subroutine i_sort_using_insertion

end module