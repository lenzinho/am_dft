MAKE      = make
FC        = ifort
OPTS      = -DCOLOR 
MKL       = -mkl
FFLAGS    = $(MKL) -fpp -O5 -openmp -implicitnone -assume realloc_lhs
FYPPOPTS  = -DCOLOR=1
AMLIB     = libam.a
DIRLIB    = library
LPATH     = -L../$(DIRLIB)
IPATH     = -I../$(DIRLIB)

FYPPFLAGS = $(FYPPOPTS) -DDEBUG=1
DEBUG     = -g -traceback -warn all -traceback -check all -gen-interfaces \
		    -warn interfaces -check uninit -ftrapuv -fp-stack-check -fpe0 -check bounds \
		    -check noarg_temp_created -debug all -diag-disable 10182 -error_limit 10
