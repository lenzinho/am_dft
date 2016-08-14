#!/bin/bash


sym
tbvsk -pair_cutoff 3.2
tbfit -file eigenval -shift fermi -skip 5
ibz -n 10 10 10


