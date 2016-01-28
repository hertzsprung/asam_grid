################################################################################
###                                                                          ###
###                          M A K E F I L E                                 ###
###                                                                          ###
###                    for the ASAM Grid Generator                           ###
###                                                                          ###
###      Institute for Tropospheric Research (TROPOS), Leipzig, Germany      ###
###                                                                          ###
################################################################################


IN = M_DEF
include ${IN}

GRIDLIB = -L./$(LIBGRID) -lgrid

OBJGRID = $(GRIDLIB)       

GridGen: GridD GRID/LIB/lib_cur.a GRID/LIB/lib_surf.a GRID/LIB/lib_gk.a
	$(LINK) $(OPT) $(INCL) -I./LIBGRID -o GridGen\
	  GRID/MainProg.f90 $(OP) $(OBJGRID) $(CURFIT) $(SURFIT) $(GK)

GridD: 
	@cd GRID; $(MAKE) -f Make_src "IN=$(IN)"

GRID/LIB/lib_cur.a:
	$(MAKE) -C GRID/LIB/CURFIT_R8 -f Make_gf

GRID/LIB/lib_surf.a:
	$(MAKE) -C GRID/LIB/SURFIT_R8 -f Make_gf

GRID/LIB/lib_gk.a:
	$(MAKE) -C GRID/LIB/GaussKrueGeo -f Make_gf
