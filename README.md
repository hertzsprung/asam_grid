ASAM cut cell grid generator
============================

[![DOI](https://zenodo.org/badge/48172254.svg)](https://zenodo.org/badge/latestdoi/48172254)

Edit `DEFPATH` in `SRC/M_DEF_Mac` to specify the ASAM_GRID directory.  Then compile with
    
    make MACH=Mac GridGen

which installs `Grid.Mac` into `$HOME/bin`.  To generate a grid, ensure the cwd is the same as the grid file (such as `Agnesi2D.grid`), then execute

    Grid.Mac Agnesi2D.grid
    
More documentation is available at http://asamwiki.tropos.de/index.php/Grid_Generator
