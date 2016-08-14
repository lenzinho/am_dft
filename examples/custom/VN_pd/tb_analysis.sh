#!/bin/bash


time sym | tee stdout.sym
time tbvsk -pair_cutoff 3.2 | tee stdout.tbvsk
time tbfit -file eigenval -shift fermi -skip 5 | tee stdout.tbfit
time ibz -n 10 10 10 | tee stdout.ibz
time tbdos | tee stdout.tbdos
time tbforce -def "./vasp/def_0/10/outfile.tb_vsk" -ref "./outfile.tb_vsk" -delta 0.05 | tee stdout.tbforce
