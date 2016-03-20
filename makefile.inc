# assume realloc_lhs Tells the compiler that when the left-hand side of an assignment 
# is an allocatable object, it should be reallocated to the shape of the right-hand side 
# of the assignment before the assignment occurs. This is the Fortran 2003 definition. 
# This feature may cause extra overhead at run time. Required for the times when an 
# allocatable array is returned by a function. This seems to be the default behavior 
# and will only cause problems when "-check all" or "-check pointer" is included in 
# the flag.

FC     = ifort
FFLAGS = -fpp -implicitnone -assume realloc_lhs -mkl -O5 -openmp -Difortprogressbar  
DEBUG  = -g -traceback  -warn all -traceback -check all -gen-interfaces -warn interfaces -check uninit -ftrapuv -fp-stack-check -fpe0 -check bounds -check noarg_temp_created
LPATH  = -L../library
IPATH  = -I../library
