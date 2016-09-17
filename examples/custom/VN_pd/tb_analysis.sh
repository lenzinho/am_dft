#!/bin/bash


time sym | tee stdout.sym
time tbbuild -pair_cutoff 3.2 | tee stdout.tbbuild
time tbfit -file eigenval -shift fermi -skip 5 | tee stdout.tbfit
time ibz -n 10 10 10 | tee stdout.ibz
time tbdos -nelecs $((18-5*2)) | tee stdout.tbdos
# time tbforce -calc dos  -def "./vasp/def_0/10/outfile.tb_vsk" -ref "./outfile.tb_vsk" -delta 0.05 | tee stdout.tbforce.dos
# time tbforce -calc path -def "./vasp/def_0/10/outfile.tb_vsk" -ref "./outfile.tb_vsk" -delta 0.05 | tee stdout.tbforce.path
time tbdr | tee stdout.tbdr
time tbdf -nelecs $((18-5*2)) -degauss 0.25 | tee stdout.tbdf
