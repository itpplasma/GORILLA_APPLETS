# Impurity tracing (i_option = 17)

Impurity guiding centres in the edge region of the extended EFIT mesh (GORILLA: `boole_extend_to_wall`,
`s_edge_region` in `tetra_grid.inp`) are traced until they hit the wall, cross the inner region boundary
or reach `t_max_s`. The background has constant density and temperature; collisions use `stost`.
Source: `SRC/impurity_tracing/impurity_tracing_mod.f90`, charge states in `marker_species_mod.f90`.

Input: `impurity_tracing.inp`, namelist `impurity_tracing_nml`.
- Default: Be1+ from a point source at the AUG outer midplane, charge-state chain Be1+ -> 2+ -> 3+ -> 4+
  with constant rates (`species_element`, `ionisation_rate_s`, `recombination_rate_s`). At a transition
  the velocity is kept and the guiding centre jumps by (b x v_perp)(1/Omega_old - 1/Omega_new) with a
  random gyrophase (`boole_gc_jump`).
- Fixed species: `species_element = ''` with `test_mass_amu`, `test_charge_number` (commented W5+ line).

Output (prefix `output_prefix`): `<prefix>_exits.dat` (status, time and position per marker) and
`<prefix>_residence.dat` (time per 2D mesh triangle, per charge state in charge-state mode).

Plots: `python region_plots.py RUNDIR [--prefix impurity_tracing]` (density, boundary hits coloured on the
contours, hits vs arc length; one density map per charge state in charge-state mode).
`python region_check.py MESHDIR` checks the region boundary of a mesh.

## Atomic data (ADAS)

Ionisation and recombination rates for Be come from OPEN-ADAS (https://open.adas.ac.uk), ADF11 files
of year 96 (the ADAS 96 data used in Harting et al., NME 42 (2025) 101887). With `adas_scd_file` and/or `adas_acd_file`
in the namelist, the rates become n_e <sigma v> at the background n_e and T_e (bilinear interpolation in
log10 n_e, log10 T_e, edge values outside the table; `SRC/impurity_tracing/adf11_mod.f90`).
Recombination to neutrals is set to zero (neutrals are not traced yet).

Rates per cell: with `background_cells_file` (one line `n_e[cm^-3] T_e[eV]` per 2D mesh triangle, in
the triangle order of `<prefix>_residence.dat`; `#` lines are comments) the ADAS rates are precomputed
per triangle and charge state. Along the orbit the rate is then piecewise constant: an exponentially
distributed optical depth is drawn and nu_cell * dt is summed per tetrahedron until it is reached, so
the transition point is sampled exactly for the given cell values. Collisions still use the constant
`electron_density_cm3`, `temperature_ev`. A background that varies inside a cell (e.g. linear per
tetrahedron) would need the rates evaluated along the orbit at the local n_e, T_e instead (not
implemented; time-step error from summing nu dt per step).
- `scd96_be.dat`: effective ionisation coefficients, Be0 -> Be1+ ... Be3+ -> Be4+
- `acd96_be.dat`: effective recombination coefficients (radiative, dielectronic, three-body), Be1+ -> Be0 ... Be4+ -> Be3+

Each file (about 27 kB) contains log10 of the rate coefficient in cm^3/s per charge state on a grid of
24 log10 n_e [cm^-3] x 25 log10 T_e [eV] values. Download (no registration needed) from the ADF11 page on
the OPEN-ADAS website or directly:

    curl -O https://open.adas.ac.uk/download/adf11/scd96/scd96_be.dat
    curl -O https://open.adas.ac.uk/download/adf11/acd96/acd96_be.dat

The files are NOT part of this repository: the OPEN-ADAS terms (at the end of each file) do not allow
redistributing the data together with a modelling code. Every user has to download them.
