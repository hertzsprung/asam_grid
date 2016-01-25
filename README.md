Edit `DEFPATH` in `SRC/M_DEF_Mac` to specify the ASAM_GRID directory.  Then compile with
    
    make MACH=Mac GridGen

which installs `Grid.Mac` into `$HOME/bin`.  To generate a grid, ensure the cwd is the same as the grid file (such as `Agnesi2D.grid`), then execute

    Grid.Mac Agnesi2D.grid
    
