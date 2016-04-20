module am_conv_cell

    use am_constants
    use am_helpers
    use am_options
    use am_unit_cell
    use am_prim_cell
    use am_symmetry
    use am_mkl

    implicit none

    private

    type, public, extends(am_class_unit_cell) :: am_class_conv_cell
        real(dp) :: centering(3,3) ! centering matrix
    contains
!         procedure :: get_conventional
    end type am_class_conv_cell

contains

! !     subroutine     get_conventional(conv,prim,opts)
! !         !
! !         use am_rank_and_sort
! !         !
! !         implicit none
! !         !
! !         class(am_class_conv_cell), intent(inout) :: conv
! !         type(am_class_prim_cell), intent(in) :: prim
! !         type(am_class_options), intent(in) :: opts
! !         type(am_class_options)  :: notalk
! !         real(dp) :: S(3,3,4) ! four types of centering matrices
! !         integer , allocatable :: P(:,:) ! permutation matrices
! !         real(dp), allocatable :: C(:,:,:)
! !         real(dp) :: bas(3,3)
! !         integer :: ID(3,3)
! !         integer :: m, j, i, k, n
! !         logical :: found
! !         !
! !         if (opts%verbosity.ge.1) call am_print_title('Transforming to conventional cell') 
! !         !
! !         notalk = opts
! !         notalk%verbosity = 0
! !         !
! !         ! construct all possible centering matrices
! !         n = 0
! !         n=n+1 ! BCC
! !         S(1:3,1:3,n) = (ones(3)-eye(3)*2.0_dp)/real(2,dp)
! !         n=n+1 ! ABC
! !         S(1,1:3,n) = [ 0, 0, 2]
! !         S(2,1:3,n) = [-1, 1, 0]
! !         S(3,1:3,n) = [ 1, 1, 0]
! !         n=n+1 ! FCC
! !         S(1:3,1:3,n) = (ones(3)-eye(3))/real(2,dp)
! !         !
! !         ID = eye(3)
! !         P = perms(3)
! !         m = sine(P,2)
! !         allocate(C(3,3,6*6*n+1))
! !         n = 0
! !         do k = 1, 3
! !         do i = 1, m
! !         do j = 1, m
! !             !
! !             n=n+1
! !             C(:,:,n) = matmul(ID(1:3,P(1:3,i)),S(P(1:3,j),1:3,k))
! !             !
! !         enddo
! !         enddo
! !         enddo
! !         !
! !         found = .false.
! !         do i = 1, n
! !             bas = matmul(inv(C(:,:,i)),prim%bas)
! !             if ( all(abs(diag(diag(bas))-bas).lt.tiny)) then
! !                 found = .true.
! !                 conv%centering = C(:,:,i)
! !                 exit
! !             endif
! !         enddo
! !         !
! !         if (.not.found) then
! !             conv%centering = eye(3)
! !         endif
! !         !
! !         if (opts%verbosity.ge.1) call am_print('centering matrix',conv%centering,' ... ')
! !         !
! !         call conv%get_supercell(uc=prim,bscfp=inv(conv%centering),opts=notalk)
! !         !
! !     end subroutine get_conventional

!     subroutine     get_conventional(conv,prim,opts)
!         !
!         use am_rank_and_sort
!         !
!         implicit none
!         !
!         class(am_class_conv_cell), intent(inout) :: conv
!         type(am_class_prim_cell), intent(in) :: prim
!         type(am_class_options), intent(in) :: opts
!         real(dp) :: centering(3,3)
!         integer :: lattice_code
!         integer, allocatable :: P(:,:) ! permutation matrices
!         integer :: i, j
! !         !
! !         if (opts%verbosity.ge.1) call am_print_title('Transforming to conventional cell') 
! !         !
! !         !
! !         P = perms(3)
! !         s : do i = 1, sine(P,2)
! !         do j = 1, sine(P,2)
! !             lattice_code = lattice_type(bas=prim%bas(P(1:3,j),P(1:3,i)),sym_prec=opts%sym_prec)
! !             if (lattice_code.ne.0) exit s
! !         enddo
! !         enddo s
! !         if (lattice_code.eq.0) then
! !             call am_print('ERROR','Unable to determine lattice type.',flags='E')
! !             stop
! !         endif
! !         !
! !         !
! !         if (opts%verbosity.ge.1) then
! !         select case(lattice_code)
! !         case(1);  call am_print('lattice','simple cubic cell',' ... ')
! !         case(2);  call am_print('lattice','body-centered cubic cell',' ... ')
! !         case(3);  call am_print('lattice','face-centered cubic cell',' ... ')
! !         case(4);  call am_print('lattice','hexagonal cell',' ... ')
! !         case(5);  call am_print('lattice','simple tetragonal cell',' ... ')
! !         case(6);  call am_print('lattice','body-centered tetragonal cell',' ... ')
! !         case(7);  call am_print('lattice','rhombohedral (trigonal) cell',' ... ')
! !         case(8);  call am_print('lattice','simple orthorhombic cell',' ... ')
! !         case(9);  call am_print('lattice','body-centered orthorhombic cell',' ... ')
! !         case(10); call am_print('lattice','face-centered orthorhombic cell',' ... ')
! !         case(11); call am_print('lattice','base centered orthorhombic cell',' ... ')
! !         case(12); call am_print('lattice','simple monoclinic cell',' ... ')
! !         case(13); call am_print('lattice','base-centered monoclinic cell',' ... ')
! !         case(14); call am_print('lattice','triclinic cell',' ... ')
! !         case default
! !             call am_print('ERROR','Unknown lattice code.')
! !             stop
! !         end select
! !         endif
! !         !
! !         centering = lattice_centering(lattice_code)
! !         !
! !         if (opts%verbosity.ge.1) call am_print('centering matrix',centering,' ... ')
! !         !
! !         call conv%get_supercell(uc=prim,bscfp=inv(centering),opts=opts)
! !         !
!     end subroutine get_conventional

!     function       lattice_type(bas,sym_prec) result(lattice_code)
!         !>
!         !> Returns 0 if unable to identify lattice type.
!         !>
!         !> get crystal system centering from the number of point symmetries and metric tensor. 
!         !> For the formalism, see:
!         !>
!         !> Symmetry and Condensed Matter Physics: A Computational Approach. 1 edition. Cambridge,
!         !> UK; New York: Cambridge University Press, 2008. p 282, table 10.3; p 288, table 10.4
!         !>
!         !>     code  lttice type
!         !>        1  simple cubic cell                            
!         !>        2  body-centered cubic cell                     
!         !>        3  face-centered cubic cell                     
!         !>        4  hexagonal cell                               
!         !>        5  simple tetragonal cell                       
!         !>        6  body-centered tetragonal cell                
!         !>        7  rhombohedral (trigonal) cell                 
!         !>        8  simple orthorhombic cell                     
!         !>        9  body-centered orthorhombic cell              
!         !>       10  face-centered orthorhombic cell              
!         !>       11  base centered orthorhombic cell              
!         !>       12  simple monoclinic cell                       
!         !>       13  base-centered monoclinic cell                
!         !>       14  triclinic cell                              
!         !> 
!         implicit none
!         !
!         real(dp), intent(in) :: bas(3,3) ! number of point group symmetries
!         real(dp), intent(in) :: sym_prec
!         integer  :: lattice_code
!         integer  :: nsyms ! number of symmetries
!         real(dp) :: M(3,3) ! metric tensor
!         !
!         nsyms = sine( lattice_symmetries(bas,sym_prec), 3)
!         !
!         M = matmul(transpose(bas),bas)
!         !
!         lattice_code = 0
!         !
!         select case (nsyms)
!         case(2) ! triclinic
!                 lattice_code = 14
!                 return
!         case(4) ! monoclinic
!             if (are_equal( [0.0_dp,M(1,3),M(2,3)] )) then
!             if (are_different( [M(1,1),M(2,2),M(3,3)] )) then
!                 lattice_code = 12
!                 return
!             endif
!             endif
!             if (are_equal( [M(1,1),M(2,2)] )) then
!             if (are_different( [M(1,2),M(1,3),M(2,3),M(3,3)] )) then
!                 lattice_code = 13
!                 return
!             endif
!             endif
!         case(8) 
!             ! simple orthorhombic
!             if (are_equal( [0.0_dp,M(1,2),M(1,3),M(2,3) ] )) then
!             if (are_different( [M(1,1),M(2,2),M(3,3)] )) then
!                 lattice_code = 8
!                 return
!             endif
!             endif
!             ! body-centered orthorhombic
!             if (are_equal( [M(1,1),M(2,2),M(3,3) ] )) then
!             if (are_different( [M(1,2),M(1,3),M(2,3)] )) then
!                 lattice_code = 9
!                 return
!             endif
!             endif
!             ! based-centered orthorhombic
!             if (are_equal( [M(1,1),M(2,2)] ).and.are_equal( [0.0_dp,M(1,3),M(2,3)] )) then
!             if (are_different( [M(2,2),M(3,3),M(2,3),M(1,2)] )) then
!                 lattice_code = 11
!                 return
!             endif
!             endif
!             ! face-centered orthorhombic
!             if (are_different( [M(1,1),M(2,2),M(3,3),M(2,3),M(1,3),M(1,2)] )) then
!                 lattice_code = 10
!                 return
!             endif
!         case(12)! trigonal (rhombahedral)
!             !  based-centered trigonal
!             if (are_equal( [M(1,1),M(2,2),M(3,3)] ).and.are_equal( [M(1,2),M(1,3),M(2,3)] )) then
!             if (are_different( [M(1,1),M(2,3)] )) then
!                 lattice_code = 7
!                 return
!             endif
!             endif
!         case(16)! tetragonal
!             ! simple tetragonal
!             if (are_equal( [M(1,1),M(2,2)] ).and.are_equal( [0.0_dp,M(1,2),M(1,3),M(2,3)] )) then
!             if (are_different( [M(1,1),M(3,3),M(1,2)] )) then
!                 lattice_code = 5
!                 return
!             endif
!             endif
!             ! body-centered tetragonal
!             if (are_equal( [M(1,1),M(2,2),M(3,3)] ).and.are_equal( [M(1,3),M(2,3)] )) then
!             if (are_different( [M(1,1),M(1,2),M(1,3)] )) then
!                 lattice_code = 5
!                 return
!             endif
!             endif
!         case(24)! hexagonal
!             ! simple hexagonal
!                 lattice_code = 4
!                 return
!         case(48)! cubic
!             ! simple cubic
!             if (are_equal( [M(1,1),M(2,2),M(3,3)] ).and.are_equal( [0.0_dp,M(1,2),M(1,3),M(2,3)] )) then
!             if (are_different( [M(1,1),M(1,2)] )) then
!                 lattice_code = 1
!                 return
!             endif
!             endif
!             ! face-centered cubic
!             if (are_equal( [M(1,1),M(2,2),M(3,3),2.0_dp*M(1,2),2.0_dp*M(1,3),2.0_dp*M(2,3)] )) then
!                 lattice_code = 3
!                 return
!             endif
!             ! body-centered cubic
!             if (are_equal( [M(1,1),M(2,2),M(3,3),-3.0_dp*M(1,2),-3.0_dp*M(1,3),-3.0_dp*M(2,3)] )) then
!                 lattice_code = 2
!                 return
!             endif
!         end select
!         !
!     end function   lattice_type

!     function       lattice_centering(lattice_code) result(centering)
!         ! 
!         implicit none
!         !
!         integer  :: lattice_code
!         real(dp), dimension(3,3) :: P,I,F,A,B,C,centering
!         !
!         call am_print('WARNING','The method used to get the centering matrix is not guarenteed to work.')
!         call am_print('WARNING','The method used to get the centering matrix is not guarenteed to work.')
!         call am_print('WARNING','The method used to get the centering matrix is not guarenteed to work.')
!         ! Not guarneteed to work because not all the matrices possible are listed down here. There can also be their permutations.
!         !
!         !
!         ! simple/primitive
!         P(1:3,1) = [1,0,0]
!         P(1:3,2) = [0,1,0]
!         P(1:3,3) = [0,0,1]
!         ! body centered cell
!         I(1:3,1) = [-1,1,1]/2.0_dp
!         I(1:3,2) = [1,-1,1]/2.0_dp
!         I(1:3,3) = [1,1,-1]/2.0_dp
!         ! face centered
!         F(1:3,1) = [0,1,1]/2.0_dp
!         F(1:3,2) = [1,0,1]/2.0_dp
!         F(1:3,3) = [1,1,0]/2.0_dp
!         ! A
!         A(1:3,1) = [ 0,0,2]/2.0_dp
!         A(1:3,2) = [-1,1,0]/2.0_dp
!         A(1:3,3) = [ 1,1,0]/2.0_dp
!         ! B
!         B(1:3,1) = [-1,1,0]/2.0_dp
!         B(1:3,2) = [ 0,0,2]/2.0_dp
!         B(1:3,3) = [ 1,1,0]/2.0_dp
!         ! C
!         C(1:3,1) = [-1,1,0]/2.0_dp
!         C(1:3,2) = [ 1,1,0]/2.0_dp
!         C(1:3,3) = [ 0,0,2]/2.0_dp

!         select case(lattice_code)
!         case(1);  centering = P ! simple cubic cell                            
!         case(2);  centering = I ! body-centered cubic cell                     
!         case(3);  centering = F ! face-centered cubic cell                     
!         case(4);  centering = P ! hexagonal cell                               
!         case(5);  centering = P ! simple tetragonal cell                       
!         case(6);  centering = I ! body-centered tetragonal cell                
!         case(7);  centering = A ! rhombohedral (trigonal) cell
!         case(8);  centering = P ! simple orthorhombic cell                     
!         case(9);  centering = I ! body-centered orthorhombic cell              
!         case(10); centering = F ! face-centered orthorhombic cell              
!         case(11); centering = A ! base centered orthorhombic cell              
!         case(12); centering = P ! simple monoclinic cell                       
!         case(13); centering = A ! base-centered monoclinic cell                
!         case(14); centering = P ! triclinic cell      
!         case default
!             call am_print('ERROR','Unknown lattice code.')        
!         end select
!         !
!     end function   lattice_centering


end module am_conv_cell