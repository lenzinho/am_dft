#!/bin/bash
#
# Excude command $2 in all folders containing $1:
#
#	rexec <$1:FILE> <$2:COMMAND>
#
# See which runs did not finish: 
# grep -LR "Voluntary" --include="OUTCAR" . 
# runs 
odir=$PWD
dlist=$(find . -name "$1" -exec sh -c '(cd `dirname {}` && pwd  )' ';' )
for d in ${dlist}; do
	$2 $d
# cd ${d}
#  if ! grep -sq "Total CPU time" OUTCAR; then
#   if [[ -e INCAR   ]]; then
#    if [[ -e POSCAR  ]]; then
#     if [[ -e KPOINTS ]]; then
#      if [[ -e POTCAR  ]]; then
#       echo $PWD
#       # if hash qsub   2>/dev/null; then qsub   QSUB; fi
#       # if hash sbatch 2>/dev/null; then sbatch QSUB; fi
#      fi
#     fi
#    fi
#   fi
#  fi
# cd ${odir}
done
