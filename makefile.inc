MAKE      = make
FC        = mpifort
OPTS      = -DCOLOR 
MKL       = -mkl -lpthread -ldl -lmkl_lapack95_lp64 -lmkl_scalapack_lp64 -lmkl_blacs_mpich_lp64
FFLAGS    = $(MKL) -fpp -O3 -implicitnone -assume realloc_lhs -assume byterecl
FYPPOPTS  = -DCOLOR=1
AMLIB     = libam.a
DIRLIB    = library
LPATH     = -L../$(DIRLIB)
IPATH     = -I../$(DIRLIB)

FYPPFLAGS = $(FYPPOPTS) -DDEBUG=1
DEBUG     = -g -traceback -warn all -traceback -check all -gen-interfaces \
		    -warn interfaces -check uninit -ftrapuv -fp-stack-check -fpe0 -check bounds \
		    -check noarg_temp_created -debug all -diag-disable 10182 -error_limit 10