all:

%:
	make -C src           $@
	make -C example       $@
	@echo " > > > Make done < < <"