#!/bin/bash

for i in "scf" "nscf"; do
	# prepare run
	ln -sf INCAR.$i INCAR
	ln -sf KPOINTS.$i KPOINTS
	# perform run
	vasp_std | tee OUTPUT.$i
	# save results
	mv EIGENVAL EIGENVAL.$i
	mv OUTCAR OUTCAR.$i
done

# ln -sf INCAR.scf INCAR
# ln -sf KPOINTS.scf KPOINTS
# vasp_std | tee OUTPUT.scf
# mv EIGENVAL EIGENVAL.scf
# mv OUTCAR OUTCAR.scf

# ln -sf INCAR.nscf INCAR
# ln -sf KPOINTS.nscf KPOINTS
# vasp_std
# vasp_std | tee OUTPUT.nscf
# mv EIGENVAL EIGENVAL.nscf
# mv OUTCAR OUTCAR.nscf

