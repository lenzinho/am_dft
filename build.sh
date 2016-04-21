#!/bin/bash

odir=`pwd`

#
# make library
#
cd ${odir}/library
	if [ ! -f makefile.inc ]; then
		ln -s ../makefile.inc .
	fi
	make clean
	make
cd ${odir}

#
# make programs
#
for program in `find . -name main.f90` ; do
	#
	d=`dirname ${program}`
	cd ${odir}/${d}
		make clean
		#
		echo ""
		echo "+----------------------------------------------------------------------------------------+"
		#
		if [ ! -f makefile.inc ]; then
			ln -s ../makefile.inc .
		fi
		#
		program_name=${PWD##*/}
		program_name=`echo ${program_name} | sed 's/prog_//'`
		#
		echo ${program_name}
		pwd
		#
		make
		#
		if [ -f ${program_name} ]; then
			if [ ! -f ../bin/${program_name} ]; then
			    ln -s ../${d}/${program_name} ../bin/
			fi
		else
			echo "ERROR! Program did not compile correctly."
			echo ${program_name}
		fi
		echo "+----------------------------------------------------------------------------------------+"
		echo ""
		#
	cd ${odir}
	#
done
