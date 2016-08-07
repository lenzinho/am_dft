#:def assertTrue(cond)
    if (.not. ${cond}$) then
      print *, "Assert failed in file ${_FILE_}$, line ${_LINE_}$"
      error stop
    end if
#:enddef



#:def write_xml_attribute_allocatable(ATTRIBUTE)
    write(unit,'(a/)') '<${ATTRIBUTE}$>'
        if (allocated(dtv%${ATTRIBUTE}$)) then
            call write_xml_attribute(unit=unit,value=.true.                        ,attribute='allocated')
            call write_xml_attribute(unit=unit,value=size(shape(dtv%${ATTRIBUTE}$)),attribute='rank')
            call write_xml_attribute(unit=unit,value=shape(dtv%${ATTRIBUTE}$)      ,attribute='shape')
            call write_xml_attribute(unit=unit,value=      dtv%${ATTRIBUTE}$       ,attribute='value')
        else
            call write_xml_attribute(unit=unit,value=.false.                       ,attribute='allocated')
        endif
    write(unit,'(a/)') '</${ATTRIBUTE}$>'
#:enddef

#:def write_xml_attribute_nonallocatable(ATTRIBUTE)
	call write_xml_attribute(unit=unit,value=dtv%${ATTRIBUTE}$,attribute='${ATTRIBUTE}$')
#:enddef

#:def write_xml_attribute_allocatable_string(ATTRIBUTE)
    write(unit,'(a/)') '<${ATTRIBUTE}$>'
        if (allocated(dtv%${ATTRIBUTE}$)) then
            call write_xml_attribute(unit=unit,value=.true.                        ,attribute='allocated')
            call write_xml_attribute(unit=unit,value= len(dtv%${ATTRIBUTE}$)       ,attribute='length')
            call write_xml_attribute(unit=unit,value=size(dtv%${ATTRIBUTE}$)       ,attribute='size')
            call write_xml_attribute(unit=unit,value=     dtv%${ATTRIBUTE}$        ,attribute='value')
        else
            call write_xml_attribute(unit=unit,value=.false.                       ,attribute='allocated')
        endif
    write(unit,'(a/)') '</${ATTRIBUTE}$>'
#:enddef



#:def read_xml_attribute_allocatable(ATTRIBUTE)
    #! write(unit,'(a/)') '<${ATTRIBUTE}$>'
    read(unit=unit,fmt='(/)',iostat=iostat,iomsg=iomsg)
        ! call write_xml_attribute(unit=unit,value=.true.,attribute='allocated')
        call read_xml_attribute(unit=unit,value=isallocated,iostat=iostat,iomsg=iomsg)
        if (isallocated) then
            #! call write_xml_attribute(unit=unit,value=size(shape(dtv%${ATTRIBUTE}$)),attribute='rank')
            call read_xml_attribute(unit=unit,value=dims_rank,iostat=iostat,iomsg=iomsg)
            #! write_xml_attribute(unit=unit,value=shape(dtv%${ATTRIBUTE}$),attribute='shape')
            call read_xml_attribute(unit=unit,value=dims(1:dims_rank),iostat=iostat,iomsg=iomsg)
            ! allocate attribute
            if (allocated(dtv%${ATTRIBUTE}$)) deallocate(dtv%${ATTRIBUTE}$)
            call vector_allocate(A=dtv%${ATTRIBUTE}$, vec=dims(1:dims_rank))
            #! call write_xml_attribute(unit=unit,value=dtv%${ATTRIBUTE}$,attribute='value')
            call read_xml_attribute(unit=unit,value=dtv%${ATTRIBUTE}$,iostat=iostat,iomsg=iomsg)
            ! write
            write(*,'(5x,a)') '${ATTRIBUTE}$('//tostring(dims(1:dims_rank))//')'
        endif
    #! write(unit,'(a)') '</${ATTRIBUTE}$>'
    read(unit=unit,fmt='(/)',iostat=iostat,iomsg=iomsg)
#:enddef

#:def read_xml_attribute_nonallocatable(ATTRIBUTE)
	call read_xml_attribute(unit=unit,value=dtv%${ATTRIBUTE}$,iostat=iostat,iomsg=iomsg)
#:enddef

#:def read_xml_attribute_allocatable_string(ATTRIBUTE)
    #! write(unit,'(a/)') '<${ATTRIBUTE}$>'
    read(unit=unit,fmt='(/)',iostat=iostat,iomsg=iomsg)
        ! call write_xml_attribute(unit=unit,value=.true.,attribute='allocated')
        call read_xml_attribute(unit=unit,value=isallocated,iostat=iostat,iomsg=iomsg)
        if (isallocated) then
            #! call write_xml_attribute(unit=unit,value=size(shape(dtv%${ATTRIBUTE}$)),attribute='rank')
            call read_xml_attribute(unit=unit,value=dims_rank,iostat=iostat,iomsg=iomsg)
            #! write_xml_attribute(unit=unit,value=shape(dtv%${ATTRIBUTE}$),attribute='shape')
            call read_xml_attribute(unit=unit,value=dims(1:1),iostat=iostat,iomsg=iomsg)
            ! allocate attribute
            if (allocated(dtv%${ATTRIBUTE}$)) deallocate(dtv%${ATTRIBUTE}$)
            allocate(character(dims_rank) :: dtv%${ATTRIBUTE}$(dims(1)))
            #! call write_xml_attribute(unit=unit,value=dtv%${ATTRIBUTE}$,attribute='value')
            call read_xml_attribute(unit=unit,value=dtv%${ATTRIBUTE}$,iostat=iostat,iomsg=iomsg)
            ! write
            write(*,'(5x,a)') '${ATTRIBUTE}$('//tostring(dims(1:dims_rank))//')'
        endif
    #! write(unit,'(a)') '</${ATTRIBUTE}$>'
    read(unit=unit,fmt='(/)',iostat=iostat,iomsg=iomsg)
#:enddef



