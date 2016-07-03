include makefile.inc

LIBS  = libam.a

OBJS  = mkl_rci.o 						\
		text_viz.o 						\
		dispmodule.o 					\
		am_constants.o              	\
		am_mkl.o                    	\
		am_matlab.o                 	\
		am_stdout.o                 	\
		am_histogram.o              	\
		am_rank_and_sort.o          	\
		am_options.o                	\
		am_vasp_io.o                	\
		am_atom.o                   	\
		am_unit_cell.o              	\
		am_prim_cell.o              	\
		am_conv_cell.o              	\
		am_symmetry_tables.o  			\
		am_symmetry.o               	\
		am_irre_cell.o              	\
		am_shells.o                 	\
		am_symmetry_relations.o 		\
		am_symmetry_rep.o 				\
		am_brillouin_zone.o             \
        am_tight_binding.o              \
		am_dispersion.o                 \
		am_optimizer.o

		# am_tet_mesh.o                   \
		# am_tetrahedra.o                 \
		# am_dos.o                        \
		# am_wannier.o                    \

all: $(LIBS)
	@echo > /dev/null

$(LIBS): $(OBJS)
	touch  libam.a
	rm     libam.a
	ar rcv libam.a $(OBJS)
	ranlib libam.a
	touch  a.i90
	rm     *.i90
	

clean:
	rm -f $(OBJS) *.mod *.o *.i90

%: %.o
	$(FC) $(FFLAGS) $(LFLAGS) $(DEBUG) $(OPTS) -o $@ $^ $(LDFLAGS)

%.o: %.f90
	$(FC) $(FFLAGS) $(LFLAGS) $(DEBUG) $(OPTS) -c $<

%.o: %.F90
	$(FC) $(FFLAGS) $(LFLAGS) $(DEBUG) $(OPTS) -c $<

.PHONY: all



















