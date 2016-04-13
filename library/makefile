include makefile.inc

LIBS = libam.a
SRCS = am_constants.f90 am_progress_bar.f90 am_helpers.f90 am_mkl.f90 am_histogram.f90 am_rank_and_sort.f90 am_options.f90 am_vasp_io.f90 am_tet_mesh.f90 am_unit_cell.f90 am_symmetry.f90 am_prim_cell.f90 am_irre_cell.f90 am_conv_cell.f90 am_symmetry_adapted_tensors.f90 am_shells.f90 am_brillouin_zone.f90 am_wannier.f90
OBJS = am_constants.o   am_progress_bar.o   am_helpers.o   am_mkl.o   am_histogram.o   am_rank_and_sort.o   am_options.o   am_vasp_io.o   am_tet_mesh.o   am_unit_cell.o   am_symmetry.o   am_prim_cell.o   am_irre_cell.o   am_conv_cell.o   am_symmetry_adapted_tensors.o   am_shells.o   am_brillouin_zone.o   am_wannier.o  

all: $(LIBS)

$(LIBS): $(OBJS)
	-rm libam.a
	ar rcv libam.a $(OBJS)
	ranlib libam.a

clean:
	rm -f $(OBJS) *.mod *.o *.i90

%: %.o
	$(FC) $(FFLAGS) $(LFLAGS) $(DEBUG) -o $@ $^ $(LDFLAGS)

%.o: %.f90
	$(FC) $(FFLAGS) $(LFLAGS) $(DEBUG) -c $<

%.o: %.F90
	$(FC) $(FFLAGS) $(LFLAGS) $(DEBUG) -c $<



