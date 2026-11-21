# Impurity tracing (i_option = 17)

Impurity guiding centres in the edge region of the extended EFIT mesh (GORILLA: `boole_extend_to_wall`,
`s_edge_region` in `tetra_grid.inp`) are traced until they hit the wall, cross the inner region boundary
or reach `t_max_s`. The background has constant density and temperature; collisions use `stost`.
Source: `SRC/impurity_tracing/impurity_tracing_mod.f90`, charge states in `marker_species_mod.f90`.

Input: `impurity_tracing.inp`, namelist `impurity_tracing_nml`.
- Default: Be0 from a point at the AUG outer midplane, 25 eV in -R (Harting et al. 2025, Sec. 5). The
  neutral flies on a straight line (`neutral_flight_mod`) until it ionises with `ionisation_rate_s(0)`,
  leaves the mesh or reaches `t_max_s`. At ionisation the ion keeps the velocity and its guiding centre
  is x - (b x v)/Omega (gyrophase known). Face test `neutral_face_test`: `hybrid` (exact, any toroidal
  resolution), `cartesian` (planar Cartesian faces, needs fine toroidal slices), `auto` (cartesian if
  the face sag / inradius <= `face_sag_tol`).
  Then the charge-state chain Be1+ -> 2+ -> 3+ -> 4+ with constant rates (`species_element`, `ionisation_rate_s`, `recombination_rate_s`). At a transition
  the velocity is kept and the guiding centre jumps by (b x v_perp)(1/Omega_old - 1/Omega_new) with a
  random gyrophase (`boole_gc_jump`).
- Fixed species: `species_element = ''` with `test_mass_amu`, `test_charge_number` (commented W5+ line).

Output (prefix `output_prefix`): `<prefix>_exits.dat` (status, time and position per marker),
`<prefix>_residence.dat` (time per 2D mesh triangle, per charge state 0 ... Z in charge-state mode) and,
with neutrals, `<prefix>_neutrals.dat` (end of each neutral flight: ionised / wall / time, position).

Plots: `python region_plots.py RUNDIR [--prefix impurity_tracing]` (density, boundary hits coloured on the
contours, hits vs arc length; one density map per charge state in charge-state mode).
`python region_check.py MESHDIR` checks the region boundary of a mesh.
