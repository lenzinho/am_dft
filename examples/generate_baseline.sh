#!/bin/bash


odir=`pwd`

./quick_make.sh sym
./quick_make.sh tb

dlist=`find . -maxdepth 1 -mindepth 1 -type d `
for d in ${dlist}
do
	echo ${odir}/${d}
	cd ${d}
	sym -poscar outfile.supercell | tee baseline.sym
	tb  -poscar outfile.supercell | tee baseline.tb
	cd ${odir}
done
