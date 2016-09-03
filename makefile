PROGS   = uc sym tensor tbvsk tbfit tbforce tbdr ibz tbdos wave

include makefile.inc

TIME_STAMP=`date +%b_%d_%H:%M`

all: $(PROGS)

uc: $(AMLIB) 
	(cd prog_uc; $(MAKE) )

sym: $(AMLIB) 
	(cd prog_sym; $(MAKE) )

tensor: $(AMLIB) 
	(cd prog_tensor; $(MAKE) )

tbvsk: $(AMLIB) 
	(cd prog_tbvsk; $(MAKE) )

tbfit: $(AMLIB) 
	(cd prog_tbfit; $(MAKE) )

tbdos: $(AMLIB) 
	(cd prog_tbdos; $(MAKE) )

tbforce: $(AMLIB) 
	(cd prog_tbforce; $(MAKE) )

tbdr: $(AMLIB) 
	(cd prog_tbdr; $(MAKE) )

ibz: $(AMLIB) 
	(cd prog_ibz; $(MAKE) )

wave: $(AMLIB) 
	(cd prog_wave; $(MAKE) )

$(AMLIB): 
	(cd $(DIRLIB); $(MAKE) )

clean:
	(cd prog_uc;   	 $(MAKE) clean)
	(cd prog_sym;  	 $(MAKE) clean)
	(cd prog_tensor; $(MAKE) clean)
	(cd prog_tbvsk;  $(MAKE) clean)
	(cd prog_tbfit;  $(MAKE) clean)
	(cd prog_tbdos;  $(MAKE) clean)
	(cd prog_tbforce;$(MAKE) clean)
	(cd prog_tbdr;   $(MAKE) clean)
	(cd prog_ibz; 	 $(MAKE) clean)
	(cd prog_wave; 	 $(MAKE) clean)
	(cd $(DIRLIB); 	 $(MAKE) clean)





















