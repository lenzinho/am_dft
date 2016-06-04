MAKE    = make -s
FC      = ifort
FFLAGS  = -fpp -O5 -mkl -openmp -implicitnone -assume realloc_lhs
DEBUG   = -g -traceback -warn all -traceback -check all -gen-interfaces -warn interfaces -check uninit -ftrapuv -fp-stack-check -fpe0 -check bounds -check noarg_temp_created -diag-disable 10182 -error_limit 10
AMLIB   = libam.a
DIRLIB  = library
LPATH   = -L../$(DIRLIB)
IPATH   = -I../$(DIRLIB)