module am_atom
	!
	use am_constants
	use am_stdout
	use am_options
	!
	implicit none

	public

	type am_class_atom
		!
        integer :: norbitals
        integer, allocatable :: orbital(:,:) ! quantum numbers [n,l,m,s,#], last number in nlms array is the state index starting with H 1s
        !
    contains
        procedure :: generate_orbitals
	end type am_class_atom

    contains
    
! 	subroutine     generate_orbitals(atom,Z,n_valence_states,n_excited_states)
! 		!
! 		! nvalence = number of valence orbitals
! 		! nexcited = number of excited orbitals
! 		!
! 		use am_rank_and_sort
! 		!
! 		implicit none
! 	  	!
!         class(am_class_atom), intent(inout) :: atom
! 		integer, intent(in) :: Z, n_valence_states, n_excited_states
! 		integer, allocatable :: ind(:)
! 		integer, allocatable :: nlms(:,:) ! quantum numbers, [n,l,m,s,#] state index starting from the first electron on H 1s] last index is electron count
! 		integer, allocatable :: state(:,:)
! 		integer :: nmax, nstates
! 		integer :: n,l,m,s,k
! 		!
! 		!
! 		! generate all possible states for atoms in the periodic table
! 		nmax = 6
! 		nstates = 182
! 		allocate(state(5,nstates))
! 		k=0
! 		do n = 1, nmax 
! 		do l = 0, (n-1)
! 		do m = -l,l
! 		do s = 1,2
! 			k = k+1
! 			state(1:4,k) = [n,l,m,s]
! 		enddo
! 		enddo
! 		enddo
! 		enddo
! 		!
! 		! Orbitals are filled in the order of increasing n+l;
! 		! Where two orbitals have the same value of n+l, they are filled in order of increasing n.
! 		allocate(ind(k))
! 		call rank(state(1,:)+state(2,:),ind)
! 		state = state(:,ind)
! 		state(5,:)=[1:k]
! 		!
! 		allocate(atom%state,source=state(:,max(1,(Z-n_valence_states)):(Z+n_excited_states)))
! 		!
! 	end subroutine generate_orbitals

	pure function  atm_symb(Z) result(symb)
		!
		implicit none
		!
		integer, intent(in) :: Z
		character(3), dimension(1:116) :: symbs
		character(3) :: symb
		!
		symbs = [ &
			 "H  ","He ","Li ","Be ","B  ","C  ","N  ","O  ","F  ","Ne ","Na ","Mg ","Al ","Si ","P  ","S  ", &
			 "Cl ","Ar ","K  ","Ca ","Sc ","Ti ","V  ","Cr ","Mn ","Fe ","Co ","Ni ","Cu ","Zn ","Ga ","Ge ", &
			 "As ","Se ","Br ","Kr ","Rb ","Sr ","Y  ","Zr ","Nb ","Mo ","Tc ","Ru ","Rh ","Pd ","Ag ","Cd ", &
			 "In ","Sn ","Sb ","Te ","I  ","Xe ","Cs ","Ba ","La ","Ce ","Pr ","Nd ","Pm ","Sm ","Eu ","Gd ", &
			 "Tb ","Dy ","Ho ","Er ","Tm ","Yb ","Lu ","Hf ","Ta ","W  ","Re ","Os ","Ir ","Pt ","Au ","Hg ", &
			 "Tl ","Pb ","Bi ","Po ","At ","Rn ","Fr ","Ra ","Ac ","Th ","Pa ","U  ","Np ","Pu ","Am ","Cm ", &
			 "Bk ","Cf ","Es ","Fm ","Md ","No ","Lr ","Rf ","Db ","Sg ","Bh ","Hs ","Mt ","Ds ","Rg ","Uub", &
			 "Uut","Uuq","Uup","Uuh"]
		!
		symb=symbs(Z)
		!
	end function   atm_symb

	pure function  atm_Z(symb) result(Z)
		!
		implicit none
		!
		character(*), intent(in) :: symb
		character(3), dimension(1:116) :: symbs
		integer :: Z
		!
		symbs = [ &
			 "h  ","he ","li ","be ","b  ","c  ","n  ","o  ","f  ","ne ","na ","mg ","al ","si ","p  ","s  ", &
			 "cl ","ar ","k  ","ca ","sc ","ti ","v  ","cr ","mn ","fe ","co ","ni ","cu ","zn ","ga ","ge ", &
			 "as ","se ","br ","kr ","rb ","sr ","y  ","zr ","nb ","mo ","tc ","ru ","rh ","pd ","ag ","cd ", &
			 "in ","sn ","sb ","te ","i  ","xe ","cs ","ba ","la ","ce ","pr ","nd ","pm ","sm ","eu ","gd ", &
			 "tb ","dy ","ho ","er ","tm ","yb ","lu ","hf ","ta ","w  ","re ","os ","ir ","pt ","au ","hg ", &
			 "tl ","pb ","bi ","po ","at ","rn ","fr ","ra ","ac ","th ","pa ","u  ","np ","pu ","am ","cm ", &
			 "bk ","cf ","es ","fm ","md ","no ","lr ","rf ","db ","sg ","bh ","hs ","mt ","ds ","rg ","uub", &
			 "uut","uuq","uup","uuh"]
		!
		do Z = 1, 116
		if (index(symbs(Z),trim(lowercase(symb))).ne.0) then
			return
		 endif
		enddo
		!
	end function   atm_Z

end module