# Impurity tracing (i_option = 17)

Impurity guiding centres in the edge region of the extended EFIT mesh (GORILLA: `boole_extend_to_wall`,
`s_edge_region` in `tetra_grid.inp`) are traced until they hit the wall, cross the inner region boundary
or reach `t_max_s`. The background has constant density and temperature; collisions use `stost`;
optional anomalous diffusion `d_perp_m2s`.
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

Cross-field transport: besides the guiding-centre drifts, `d_perp_m2s > 0` adds an anomalous random walk
perpendicular to B after every orbit step of an ion (the same displacement as i_option 13,
`anomalous_transport_displacement`: Delta x = sqrt(2 dt) alpha xi + V_c dt with alpha alpha^T = D_perp (g - h h)).
The steps are limited to dt <= `d_perp_max_step_cm`^2 / (2 D_perp) (default 0.1 cm). A step that leaves the
mesh ends the marker at the wall (status 1), a step into the core at the inner boundary (status 2).
Neutrals do not diffuse.

Toroidally resolved output: `residence_output = '3d'` writes the residence time per prism (2D triangle x
toroidal slice), rows slice by slice; `'2d'` (default) sums over the slices of one field period. Use
`start_phi_rad` (>= 0) for a toroidally localised source (e.g. a puff at phi = 0); negative = random phi.
A periodic sector, e.g. the 40 degree sector of Harting et al. (9 puffs in the torus), is set in
`tetra_grid.inp` with `boole_n_field_periods = .false.`, `n_field_periods_manual = 9`, and `n2` slices per sector.
`<prefix>_cells.dat` lists the triangle corners, area, centroid R and the cell volumes (prism of one slice
and ring of one field period).

Absolute density for a source of S particles/s: `python impurity_density.py RUNDIR --source S [--phi DEG]
[--cut R1 Z1 R2 Z2]`, n = (S/N) T_res / V [m^-3], map of the plane phi (3d; on a slice boundary the two
adjacent slices) or the toroidal average (2d), and with `--cut` the profile per charge state along a line.
Markers that reach `t_max_s` are counted only up to t_max_s; their fraction is printed (should be small).

Plots: `python region_plots.py RUNDIR [--prefix impurity_tracing]` (density, boundary hits coloured on the
contours, hits vs arc length; one density map per charge state in charge-state mode).
`python region_check.py MESHDIR` checks the region boundary of a mesh.

## Atomic data (ADAS)

Ionisation and recombination rates for Be come from OPEN-ADAS (https://open.adas.ac.uk), ADF11 files
of year 96 (the ADAS 96 data used in Harting et al., NME 42 (2025) 101887). With `adas_scd_file` and/or `adas_acd_file`
in the namelist, the rates become n_e <sigma v> at the background n_e and T_e (bilinear interpolation in
log10 n_e, log10 T_e, edge values outside the table; `SRC/impurity_tracing/adf11_mod.f90`).
Recombination Be1+ -> Be0 (acd block Z1 = 1, or `recombination_rate_s(1)`) turns the ion back into a neutral: particle position and velocity from the guiding centre with a random gyrophase, then the straight-line flight as for the launch. `<prefix>_neutrals.dat` lists the last neutral flight of each marker.

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
