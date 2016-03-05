module am_options
    !
    use am_constants 
    !
    type am_class_options
        integer :: verbosity
        real(dp) :: sym_prec
    contains
        procedure :: defaults
    end type am_class_options    
    !
    contains
    
    pure subroutine defaults(opts)
        !
        implicit none
        !
        class(am_class_options), intent(inout) :: opts
        !
        opts%verbosity = 1
        opts%sym_prec = tiny
        !
    end subroutine defaults
end module
    