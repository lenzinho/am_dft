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
            call write_xml_attribute(unit=unit,value=size(shape(dtv%${ATTRIBUTE}$)),attribute='rank')
            call write_xml_attribute(unit=unit,value=     shape(dtv%${ATTRIBUTE}$) ,attribute='shape')
            call write_xml_attribute(unit=unit,value=       len(dtv%${ATTRIBUTE}$) ,attribute='length')
            call write_xml_attribute(unit=unit,value=          dtv%${ATTRIBUTE}$   ,attribute='value')
        else
            call write_xml_attribute(unit=unit,value=.false.                       ,attribute='allocated')
        endif
    write(unit,'(a/)') '</${ATTRIBUTE}$>'
#:enddef

#:def write_xml_attribute_allocatable_derivedtype(ATTRIBUTE)
    write(unit,'(a/)') '<${ATTRIBUTE}$>'
        if (allocated(dtv%${ATTRIBUTE}$)) then
            call write_xml_attribute(unit=unit,value=.true.                        ,attribute='allocated')
            call write_xml_attribute(unit=unit,value=size(shape(dtv%${ATTRIBUTE}$)),attribute='rank')
            call write_xml_attribute(unit=unit,value=shape(dtv%${ATTRIBUTE}$)      ,attribute='shape')
            write(unit,'(a,/)') '<value>'
                write(unit,*) dtv%${ATTRIBUTE}$
                write(unit,'(/)')
            write(unit,'(a,/)') '</value>'
        else
            call write_xml_attribute(unit=unit,value=.false.                       ,attribute='allocated')
        endif
    write(unit,'(a/)') '</${ATTRIBUTE}$>'
#:enddef



#:def read_xml_attribute_allocatable(ATTRIBUTE)
    read(unit=unit,fmt='(/)',iostat=iostat,iomsg=iomsg)
        call read_xml_attribute(unit=unit,value=isallocated,iostat=iostat,iomsg=iomsg)			! allocated
        if (isallocated) then
            call read_xml_attribute(unit=unit,value=dims_rank,iostat=iostat,iomsg=iomsg) 		! rank
            call read_xml_attribute(unit=unit,value=dims(1:dims_rank),iostat=iostat,iomsg=iomsg)! shape
            if (allocated(dtv%${ATTRIBUTE}$)) deallocate(dtv%${ATTRIBUTE}$)
            call vector_allocate(A=dtv%${ATTRIBUTE}$, vec=dims(1:dims_rank))
            call read_xml_attribute(unit=unit,value=dtv%${ATTRIBUTE}$,iostat=iostat,iomsg=iomsg)! value
            #! write(*,'(5x,a)') '${ATTRIBUTE}$('//tostring(dims(1:dims_rank))//')'
        endif
    read(unit=unit,fmt='(/)',iostat=iostat,iomsg=iomsg)
#:enddef

#:def read_xml_attribute_nonallocatable(ATTRIBUTE)
	call read_xml_attribute(unit=unit,value=dtv%${ATTRIBUTE}$,iostat=iostat,iomsg=iomsg)
#:enddef

#:def read_xml_attribute_allocatable_string(ATTRIBUTE)
    read(unit=unit,fmt='(/)',iostat=iostat,iomsg=iomsg)
        call read_xml_attribute(unit=unit,value=isallocated,iostat=iostat,iomsg=iomsg)
        if (isallocated) then
            call read_xml_attribute(unit=unit,value=dims_rank,iostat=iostat,iomsg=iomsg) 	! rank
            call read_xml_attribute(unit=unit,value=dims(1:1),iostat=iostat,iomsg=iomsg) 	! shape
        	call read_xml_attribute(unit=unit,value=str_length,iostat=iostat,iomsg=iomsg) 	! length
            if (allocated(dtv%${ATTRIBUTE}$)) deallocate(dtv%${ATTRIBUTE}$)
            call vector_allocate(A=dtv%${ATTRIBUTE}$, vec=dims(1:dims_rank), str_length=str_length)
            call read_xml_attribute(unit=unit,value=dtv%${ATTRIBUTE}$,iostat=iostat,iomsg=iomsg) ! values
            #! write(*,'(5x,a)') '${ATTRIBUTE}$('//tostring(dims(1:dims_rank))//')'
        endif
    read(unit=unit,fmt='(/)',iostat=iostat,iomsg=iomsg)
#:enddef

#:def read_xml_attribute_allocatable_derivedtype(ATTRIBUTE)
    read(unit=unit,fmt='(/)',iostat=iostat,iomsg=iomsg)
        call read_xml_attribute(unit=unit,value=isallocated,iostat=iostat,iomsg=iomsg)          ! allocated
        if (isallocated) then
            !
            call read_xml_attribute(unit=unit,value=dims_rank,iostat=iostat,iomsg=iomsg)        ! rank
            if (dims_rank.gt.1) stop 'ERROR [read_uc]: allocation of user-defined derived-type tensor not yet implemented'
            ! allocate
            call read_xml_attribute(unit=unit,value=dims(1:dims_rank),iostat=iostat,iomsg=iomsg)! shape
            if (allocated(dtv%${ATTRIBUTE}$)) deallocate(dtv%${ATTRIBUTE}$)
            allocate(dtv%${ATTRIBUTE}$(dims(1)))
            ! read nested object
            read(unit=unit,fmt='(/)',iostat=iostat,iomsg=iomsg)
                read(unit,*) dtv%${ATTRIBUTE}$
                read(unit=unit,fmt='(/)',iostat=iostat,iomsg=iomsg)
            read(unit=unit,fmt='(/)',iostat=iostat,iomsg=iomsg)
        endif
    read(unit=unit,fmt='(/)',iostat=iostat,iomsg=iomsg)
#:enddef


