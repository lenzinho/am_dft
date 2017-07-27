#!/bin/bash
# See which runs did not finish: 
# grep -LR "Voluntary" --include="OUTCAR" . 
# 
odir=$PWD
dlist=$(find . -name "QSUB" -exec sh -c '(cd `dirname {}` && pwd  )' ';' )
for d in ${dlist}; do
cd ${d}
 if ! grep -sq "Total CPU time" OUTCAR; then
  if [[ -e INCAR   ]]; then
   if [[ -e POSCAR  ]]; then
    if [[ -e KPOINTS ]]; then
     if [[ -e POTCAR  ]]; then
      echo $PWD
      if hash qsub   2>/dev/null; then qsub   QSUB; fi
      if hash sbatch 2>/dev/null; then sbatch QSUB; fi
     fi
    fi
   fi
  fi
 fi
cd ${odir}
done
