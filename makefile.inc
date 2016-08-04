MAKE      = make
FC        = ifort
OPTS      = -DCOLOR
MKL       = -mkl
FFLAGS    = -fpp -O5 -openmp -implicitnone -assume realloc_lhs $(MKL) 
FYPPFLAGS = -DDEBUG=1
DEBUG     = -DDEBUG -g -traceback -warn all -traceback -check all -gen-interfaces \
		    -warn interfaces -check uninit -ftrapuv -fp-stack-check -fpe0 -check bounds \
		    -check noarg_temp_created -debug all -diag-disable 10182 -error_limit 10
#DEBUG     = 
AMLIB     = libam.a
DIRLIB    = library
LPATH     = -L../$(DIRLIB)
IPATH     = -I../$(DIRLIB)