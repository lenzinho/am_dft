PROGS   = tb 2fc uc wan sym ibz dos bc test

include makefile.inc

TIME_STAMP=`date +%b_%d_%H:%M`


all: $(PROGS)

tb: $(AMLIB) 
	(cd prog_tb; $(MAKE) )

2fc: $(AMLIB) 
	(cd prog_2fc; $(MAKE) )

uc: $(AMLIB) 
	(cd prog_uc; $(MAKE) )

wan: $(AMLIB) 
	(cd prog_wan; $(MAKE) )

sym: $(AMLIB) 
	(cd prog_sym; $(MAKE) )

dos: $(AMLIB) 
	(cd prog_dos; $(MAKE) )

ibz: $(AMLIB) 
	(cd prog_ibz; $(MAKE) )

bc: $(AMLIB) 
	(cd prog_bc; $(MAKE) )

test: $(AMLIB) 
	(cd prog_test; $(MAKE) )

$(AMLIB): 
	(cd prog_test; $(MAKE) )

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

# save:
# 	(cd $(COMDIR); $(MAKE) "SUF=$(SUF)" "TAR=$(TAR)" save)
# 	(cd $(PWDIR); $(MAKE) "SUF=$(SUF)" "TAR=$(TAR)" save)
# 	(cd $(TBDIR); $(MAKE) "SUF=$(SUF)" "TAR=$(TAR)" save)
# 	$(TAR) $(PROJ)_$(TIME_STAMP).tar $(COMDIR)/$(COMPROJ).tar \
# 		$(PWDIR)/$(PWPROJ).tar $(TBDIR)/$(TBPROJ).tar makefile
# 	rm $(COMDIR)/$(COMPROJ).tar $(PWDIR)/$(PWPROJ).tar \
# 		$(TBDIR)/$(TBPROJ).tar
# 	gzip $(PROJ)_*.tar
# #	mv $(PROJ)_*.tar.gz $(HOME)/save/.






