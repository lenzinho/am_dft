module am_sort

   use am_constants

   contains

   function     dlasrt2(v) result(key)
      !
      implicit none
      !
      real(dp), intent(in) :: v(:)
      real(dp), allocatable :: d(:)
      integer, allocatable :: key(:)
      integer :: select, endd, i, j, start, stkpnt, tmpkey, stack(2,32), n
      real(dp) :: d1, d2, d3, dmnmx, tmp
      !
      select = 20
      !
      n = size(v)
      !
      allocate(d,source=v)
      !
      allocate(key,source=[1:n])
      !
      if (n.le.1) return
      !
      stkpnt = 1
      stack(1,1) = 1
      stack(2,1) = n
      !
      do
         start = stack(1, stkpnt)
         endd = stack(2, stkpnt)
         stkpnt = stkpnt - 1
         if(endd-start.gt.0) then
            do i = start + 1, endd
               do j = i, start + 1, -1
                  if(d(j).lt.d(j-1)) then
                     dmnmx = d(j)
                     d(j) = d(j-1)
                     d(j-1) = dmnmx
                     tmpkey = key(j)
                     key(j) = key(j-1)
                     key(j-1) = tmpkey
                  else
                     exit
                  end if
               enddo
            enddo
         else if(endd-start.gt.select) then
            d1 = d(start)
            d2 = d(endd)
            i = (start+endd) / 2
            d3 = d(i)
            if(d1.lt.d2) then
               if(d3.lt.d1) then
                  dmnmx = d1
               else if(d3.lt.d2) then
                  dmnmx = d3
               else
                  dmnmx = d2
               end if
            else
               if(d3.lt.d2) then
                  dmnmx = d2
               else if(d3.lt.d1) then
                  dmnmx = d3
               else
                  dmnmx = d1
               end if
            end if
            i = start - 1
            j = endd + 1
            do
               do
                  j = j - 1
                  if (.not.(d(j).gt.dmnmx)) exit
               enddo
               do
                  i = i + 1
                  if (.not.(d(i).lt.dmnmx)) exit
               enddo
               if(i.gt.j) then
                  exit   
               else 
                  tmp = d(i)
                  d(i) = d(j)
                  d(j) = tmp
                  tmpkey = key(j)
                  key(j) = key(i)
                  key(i) = tmpkey
               end if
            enddo
            if(j-start.gt.endd-j-1) then
               stkpnt = stkpnt + 1
               stack(1, stkpnt) = start
               stack(2, stkpnt) = j
               stkpnt = stkpnt + 1
               stack(1, stkpnt) = j + 1
               stack(2, stkpnt) = endd
            else
               stkpnt = stkpnt + 1
               stack(1, stkpnt) = j + 1
               stack(2, stkpnt) = endd
               stkpnt = stkpnt + 1
               stack(1, stkpnt) = start
               stack(2, stkpnt) = j
            end if
         end if
         if (.not.(stkpnt.gt.0)) exit
      enddo
      return
   end function dlasrt2

end module am_sort

program main
   !
   use am_constants
   use am_sort
   !
   implicit none
   !
   integer :: i
   real(dp) :: d(55)

   d(1:5)=[5,4,2,3,1]
   d(6:55)=[6:55]
   d([16,10]) =d([10,16])

   write(*,*) d
   write(*,*) dlasrt2(d)
   write(*,*) d(dlasrt2(d))

end program main

