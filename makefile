

FC = ifort 
DEBUG = -g -traceback  -warn all -traceback -check all -gen-interfaces -warn interfaces -check uninit -ftrapuv -fp-stack-check -fpe0 -check bounds

# assume realloc_lhs Tells the compiler that when the left-hand side of an assignment 
# is an allocatable object, it should be reallocated to the shape of the right-hand side 
# of the assignment before the assignment occurs. This is the Fortran 2003 definition. 
# This feature may cause extra overhead at run time. Required for the times when an 
# allocatable array is returned by a function. This seems to be the default behavior 
# and will only cause problems when "-check all" or "-check pointer" is included in 
# the flag.

FFLAGS = -implicitnone -assume realloc_lhs -mkl
CPP = icc -E 
CPFLAGS = -Difortprogressbar -Ddebug 

PROG = am_sym
SRCS = am_constants.f90 am_helpers.f90 am_rank_and_sort.f90 am_options.f90 am_vasp_io.f90 am_tet_mesh.f90 am_unit_cell.f90 am_brillouin_zone.f90 am_toy_models.f90 main.f90
OBJS = am_constants.o   am_helpers.o   am_rank_and_sort.o   am_options.o   am_vasp_io.o   am_tet_mesh.o   am_unit_cell.o   am_brillouin_zone.o   am_toy_models.o   main.o  



all: $(PROG)

$(PROG): $(OBJS)
	$(FC) $(FFLAGS) $(DEBUG) -o $@ $(OBJS) $(LIBS)

clean:
	rm -f $(PROG) $(OBJS) *.mod *.o

%: %.o
	$(FC) $(FFLAGS) $(DEBUG) -o $@ $^ $(LDFLAGS)

%.o: %.f90
	$(FC) $(FFLAGS) $(DEBUG) -c $<

%.o: %.F90
	$(FC) $(FFLAGS) $(DEBUG) -c $<



