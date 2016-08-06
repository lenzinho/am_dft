include makefile.inc

# first suffix call is mandatory; clears suffix list
.SUFFIXES: 
.SUFFIXES: .fpp .i90 .f90

LIBS  = libam.a

# goal get rid of am_stdout

FYPPINC = fypp_inc.py

OBJS  = mkl_rci.o 						\
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
		am_symmetry_tensor.o 			\
		am_brillouin_zone.o             \
		am_dispersion.o                 \
		am_dos.o                        \
		am_dft.o 						\
        am_tight_binding.o              

		# am_tet_mesh.o                   \
		# am_tetrahedra.o                 \
		# am_wannier.o                    \

$(LIBS): $(OBJS)
	touch  libam.a
	rm     libam.a
	ar rcv libam.a $(OBJS)
	ranlib libam.a
	# touch  a.i90
	# rm     *.i90
	
clean:
	rm -f $(OBJS) *.mod *.o *.i90

%.o: %.f90
	$(FC) $(FFLAGS) $(LFLAGS) $(DEBUG) $(OPTS) -c $<

%.f90: %.fpp
	fypp $(FYPPFLAGS) -i $(FYPPINC) $< $@

.PHONY: all



















