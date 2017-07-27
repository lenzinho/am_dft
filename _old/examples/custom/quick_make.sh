#!/bin/bash

ODIR=$PWD
cd /Volumes/Lenzinho/MyLinux/calc.vasp.5.3.3/development/am_lib/
make $1
cd ${ODIR}

