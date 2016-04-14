module am_atom
	!
	use am_constants
	!
	implicit none

	public

	contains

    pure function   symb(Z)
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
    end function    symb

    pure function   Z(symb)
        !
        implicit none
        !
        character(3), intent(in) :: symb
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
    end function    Z

end module