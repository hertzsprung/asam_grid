################################################################################
###                                                                          ###
###                          M A K E F I L E                                 ###
###                                                                          ###
###                    for the ASAM Grid Generator                           ###
###                                                                          ###
###      Institute for Tropospheric Research (TROPOS), Leipzig, Germany      ###
###                                                                          ###
################################################################################


IN = M_DEF_$(MACH)
include ${IN}

GRIDLIB = -L$(LIBGRID) -lgrid

OBJGRID = $(GRIDLIB)       

GridGen: GridD 
	$(LINK) $(OPT) $(INCL) $(KGRIDINCL) -o GridGen\
        GRID/MainProg.f90 $(OP) $(OBJGRID) $(CURFIT) $(SURFIT) $(GK) ;\
        mv GridGen ../LIBGRID_$(MACH)/Grid.$(MACH);\
        rm -f ~/bin/Grid.$(MACH);\
        cp $(DEFPATH)/LIBGRID_$(MACH)/Grid.$(MACH) ~/bin/Grid.$(MACH);

GridCompose:
	$(LINK) -o GridCompose GRID/GridCompose.f90;\
        mv GridCompose ../LIBGRID_$(MACH)/GridCompose.$(MACH);\
        rm -f ~/bin/GridCompose.$(MACH);\
        cp $(DEFPATH)/LIBGRID_$(MACH)/GridCompose.$(MACH) ~/bin/GridCompose.$(MACH);

GridD: 
	@cd GRID; make -f Make_src "IN=$(IN)"
