PROGS   = tbvsk tbfit uc sym

include makefile.inc

TIME_STAMP=`date +%b_%d_%H:%M`

all: $(PROGS)

uc: $(AMLIB) 
	(cd prog_uc; $(MAKE) )

sym: $(AMLIB) 
	(cd prog_sym; $(MAKE) )

tbvsk: $(AMLIB) 
	(cd prog_tbvsk; $(MAKE) )

tbfit: $(AMLIB) 
	(cd prog_tbfit; $(MAKE) )

# dos: $(AMLIB) 
# 	(cd prog_dos; $(MAKE) )

# ibz: $(AMLIB) 
# 	(cd prog_ibz; $(MAKE) )

# bc: $(AMLIB) 
# 	(cd prog_bc; $(MAKE) )

# wan: $(AMLIB) 
# 	(cd prog_wan; $(MAKE) )

$(AMLIB): 
	(cd $(DIRLIB); $(MAKE) )

clean:
	(cd prog_bc; 	$(MAKE) clean)
	(cd prog_wan;   $(MAKE) clean)
	(cd prog_dos;  	$(MAKE) clean)
	(cd prog_ibz;  	$(MAKE) clean)
	(cd prog_sym;  	$(MAKE) clean)
	(cd prog_tb;   	$(MAKE) clean)
	(cd prog_uc;   	$(MAKE) clean)
	(cd prog_2fc;  	$(MAKE) clean)
	(cd prog_toy;  	$(MAKE) clean)
	(cd $(DIRLIB); 	$(MAKE) clean)






