The grid-metric binary files correspond to MITgcm standard cs32 grid
(`grid_cs32.face00[1-6].bin`, from verification experiment `tutorial_held_suarez_cs/input`)
and halo-regions were added by running exp. advect_cs with minor changes
to `model/src/ini_curvilinear_grid.F` (un-commenting lines 422-426 and adding
missing `WRITE_FULLARRAY_RS` calls).

(by J.-M. C.)
