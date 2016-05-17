module am_symmetry_rep

    use am_constants
    use am_stdout
    use am_mkl
    use am_options
    use am_symmetry

    implicit none

    public

    type am_class_rep
        integer :: nsyms ! dimensions of reps
        integer :: ndims ! dimensions of reps
        real(dp), allocatable :: rep(:,:,:)
    end type am_class_rep

    type, extends(am_class_rep) :: am_class_rep_reg
	    contains
        procedure :: create => create_regular
    end type am_class_rep_reg

	contains

	! procedures which create representations

	subroutine     create_regular(rr,seitz)
		!
		implicit none
		!
		class(am_class_rep_reg), intent(out) :: rr
        real(dp), intent(in) :: seitz(:,:,:)
		!
		rr%rep   = rep_regular(seitz=seitz)
		!
		rr%ndims = size(rr%rep,1)
        !
		rr%nsyms = size(rr%rep,3)
		!
        contains
        function       rep_regular(seitz) result(reg_rep)
            !>
            !> Generates regular representation for each operator
            !>
            !> For details, see page 73-74, especially Eqs. 4.18 Symmetry and Condensed Matter Physics: A Computational
            !> Approach. 1 edition. Cambridge, UK; New York: Cambridge  University Press, 2008.
            !>
            !> Requires identity to be the first element of group. regular representation is built from multiplication table
            !> by making sure that the identity appears along the diagonal. this sometimes causes the row to be permuted with
            !> respect to the columns. so, first find the permutation necessary to bring the identities to the center and
            !> then build regular representations from this new multiplication table.
            !>
            implicit none
            !
            real(dp), intent(in) :: seitz(:,:,:)
            integer, allocatable :: reg_rep(:,:,:)
            integer, allocatable :: multab(:,:)
            integer :: n
            integer :: i
            !
            ! obtain multiplication table
            multab = get_multiplication_table(seitz=seitz,flags='reg')
            !
            ! construct regular representation from multiplication table
            n = size(seitz,3)
            allocate(reg_rep(n,n,n))
            reg_rep = 0
            do i = 1, n
                where (multab.eq.i) reg_rep(:,:,i) = 1
            enddo
            !
        end function   rep_regular
	end subroutine create_regular

    ! procedures which operate on matrix-rep symmetry elements

    subroutine     get_irreps(rr)
    	!
    	implicit none
    	!
    	class(am_class_rep) :: rr ! any type of representation... regular, permutation, except for seitz
        !
    end subroutine get_irreps

end module am_symmetry_rep









