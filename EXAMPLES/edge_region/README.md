# Edge-region impurity tracing (i_option = 17)

Markers in the edge region of the extended EFIT mesh (GORILLA branch `impurity_region`:
`boole_extend_to_wall`, `s_edge_region`) until they hit the wall, cross the inner region boundary
or reach `t_max_s`. Background of constant density and temperature, collisions with `stost`.

- `edge_region.inp`: fixed species (W5+), point source at the AUG outer midplane.
- `edge_region_be_chain.inp`: charge-state changes (`marker_species_mod`), Be1+ -> 2+ -> 3+ -> 4+ with
  constant rates; final charge fractions follow the Bateman solution. Velocity is kept at a transition,
  the guiding centre jumps by (b x v_perp)(1/Omega_old - 1/Omega_new) with random gyrophase.

Plots: `python region_plots.py RUNDIR [--prefix edge_region]` (density, boundary hits coloured on the
contours, hits vs arc length; one density map per charge state in charge-state mode).
`python region_check.py MESHDIR` checks the region boundary of a mesh.
