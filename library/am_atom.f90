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
        integer, allocatable :: orbital(:,:) ! quantum numbers [n,l,m,s]
        integer, allocatable :: azimuthal(:) ! list of l oribtal angular momenta, used for the construction of O3 irrep rotation
        !
    contains
    	procedure :: gen_orbitals
	end type am_class_atom

    contains

	pure function   atm_symb(Z) result(symb)
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
	end function    atm_symb

	pure function   atm_Z(symb) result(Z)
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
	end function    atm_Z

	pure function   spdf2l(spdf) result(l)
		!
		character(1), intent(in) :: spdf
		integer :: l
		!
		if     (spdf(1:1).eq.'s') then; l = 0
		elseif (spdf(1:1).eq.'p') then; l = 1
		elseif (spdf(1:1).eq.'d') then; l = 2
		elseif (spdf(1:1).eq.'f') then; l = 3
		endif
	end function    spdf2l

	pure function   l2spdf(l) result(spdf)
		!
		integer, intent(in) :: l
		character(1) :: spdf
		!
		if     (l.eq.0) then; spdf = 's'
		elseif (l.eq.1) then; spdf = 'p'
		elseif (l.eq.2) then; spdf = 'd'
		elseif (l.eq.3) then; spdf = 'f'
		endif
	end function    l2spdf

	pure function   m2spdf(m) result(sigma_pi_delta_phi)
		!
		integer, intent(in) :: m
		character(5) :: sigma_pi_delta_phi
		!
		if     (m.eq.0) then; sigma_pi_delta_phi = 'sigma'
		elseif (m.eq.1) then; sigma_pi_delta_phi = 'pi'
		elseif (m.eq.2) then; sigma_pi_delta_phi = 'delta'
		elseif (m.eq.3) then; sigma_pi_delta_phi = 'phi'
		endif
	end function    m2spdf

	subroutine      gen_orbitals(atom,orbital_flags)
		!
		implicit none
		!
		class(am_class_atom) :: atom
		character(*), intent(in) :: orbital_flags
		character(2) :: orbital_name ! 1s, 2s, 2p, 3s, 3p, etc
		integer, allocatable :: slist(:)
		integer :: n, m, l, s, smax, k
		integer :: nazimuthals
		!
		!
		if (index(orbital_flags,'spin_polarized')) then
            allocate(slist,source=[-1,1])
            smax = size(slist)
		else
            allocate(slist,source=[0])
            smax = size(slist)
		endif
		!
		! just determine the sizes of atom%azimuthal and atom%orbital here 
		nazimuthals = 0
		atom%norbitals = 0
		do n = 1,4
		do l = 0, (n-1)
			!
			orbital_name(1:2) = trim(int2char(n))//trim(l2spdf(l))
			!
			if (index(orbital_flags,orbital_name).ne.0) then
				!
				nazimuthals = nazimuthals + 1
				!
				do m = -l, l
		        do s = 1, smax
			        atom%norbitals = atom%norbitals + 1
	            enddo
				enddo
			endif
		enddo
		enddo
		!
		! 
		allocate(atom%orbital(4,atom%norbitals))
		allocate(atom%azimuthal(nazimuthals))
		!
		nazimuthals = 0
		k = 0
		do n = 1, 4
		do l = 0, (n-1)
			!
			orbital_name(1:2) = trim(int2char(n))//trim(l2spdf(l))
			!
			if (index(orbital_flags,orbital_name).ne.0) then
				!
				nazimuthals = nazimuthals + 1
				atom%azimuthal(nazimuthals) = l
				!
				do m = -l, l
		        do s = 1, smax
			        k = k+1
			        atom%orbital(1,k) = n
			        atom%orbital(2,k) = l
			        atom%orbital(3,k) = m
			        atom%orbital(4,k) = slist(s)
	            enddo
				enddo
			endif
		enddo
		enddo
		!
	end subroutine  gen_orbitals

end module








