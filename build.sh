#!/bin/bash

odir=`pwd`

#
# make library
#
cd ${odir}/library
	if [ ! -f makefile.inc ]; then
		ln -s ../makefile.inc .
	fi
	make
cd ${odir}

#
# make programs
#
for program in `find . -name main.f90` ; do
	#
	d=`dirname ${program}`
	cd ${odir}/${d}
		#
		echo ""
		echo "+----------------------------------------------------------------------------------------+"
		echo ""
		pwd
		echo ""
		echo "+----------------------------------------------------------------------------------------+"
		echo ""
		#
		if [ ! -f makefile.inc ]; then
			ln -s ../makefile.inc .
		fi
		#
		make
		#
		program_name=${PWD##*/}
		program_name=`echo ${program_name} | sed 's/prog_//'`
		if [ -f ${program_name} ]; then
		if [ ! -f ../${d}/${program_name} ]; then
		    ln -s ../${d}/${program_name} ../bin/
		fi
		fi
		#
	cd ${odir}
	#
done
