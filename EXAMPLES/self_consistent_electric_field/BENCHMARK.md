# SCEF benchmark: self-consistent E_r in W7-X

Work in progress on the branch `scef_benchmark`.

## Goal

Two other codes have computed the converged, self-consistent radial electric field
E_r(r/a) for prescribed density and temperature profiles in W7-X. Their E_r profiles
differ. The aim is to compute E_r with the self-consistent electric field (SCEF)
applet of GORILLA_APPLETS (`i_option = 12`), using the same profiles. The result
either confirms one of the two profiles or gives a third one.

Both reference profiles are positive in the core (electron root) and negative
further out (ion root), with values of up to about ±15 kV/m.

## Files that are not in the repository

For privacy, these files are not on GitHub. Ask the branch owner for them.

| File | What it is |
|---|---|
| `Er_kV_m.eps` | Plot of the two reference E_r profiles, E_r [kV/m] over r/a |
| `er_michal.dat`, `er_fort401.dat` | The two curves recovered from the `.eps` (purple: 32 points, green: 128 points). Columns: r/a, E_r [kV/m] |
| `wout_EIM_1000_1000_1000_1000_+0000_+0000.01.00jh_l_ref3_w7x-sc1.nc` | W7-X EIM VMEC equilibrium (a = 0.5127 m), referenced by `tetra_grid.inp`. Place it in this folder |

## Plasma profiles

All profiles are functions of rad = r/a = sqrt(s), where s is the normalised toroidal
flux. They are defined in
`SRC/self_consistent_electric_field/scef_target_profiles_mod.f90`:

    n(rad)   = 0.37e14 * (0.74 + 0.26*(1 - rad^2.5)^1.5 - 0.06*(1 - exp(-4*rad^2)))   [cm^-3]
    T_i(rad) = 1770 * (0.06 + 0.94*(1 - rad^2.9)^2)                                 [eV]
    T_e(rad) = T_i(rad) + 3190 * (1 - rad^2.1)^5.4                                  [eV]

| r/a | T_e [eV] | T_i [eV] | n [cm^-3] |
|---|---|---|---|
| 0   | 4960 | 1770 | 3.70e13 |
| 0.5 | 2114 | 1354 | 3.32e13 |
| 0.7 |  898 |  797 | 2.98e13 |
| 0.9 |  222 |  222 | 2.63e13 |
| 1   |  106 |  106 | 2.52e13 |

The ions are deuterium (`ispecies = 2` in `gorilla.inp`).

## Model

The applet alternates between pushing markers and updating the potential and the
particle source.

- **Collision background.** Test-particle collisions with a fixed Maxwellian background
  that has the target n, T_i and T_e. Energy and momentum conservation is off, so the
  background acts as a heat bath for the markers.
- **Ions.** Traced honestly with GORILLA until they leave the domain at s = 1
  (`boole_honest_tracing(1) = .true.`).
- **Electrons.** Full tracing is too expensive, so electrons do a 1D random walk in s
  (`boole_honest_tracing(2) = .false.`). Its coefficients come from a calibration:
  - On each flux surface, electrons are started from a Maxwellian at the local T_e
    and traced with the full collision operator.
  - Each surface is traced for `diffusion_tau_collisions` local collision times
    tau_e (NRL formulary, at the target n and T_e).
  - Fitting <Δs> and <Δs²> over time gives the drift A and the diffusion coefficient B.
  - The random walk uses a quadratic fit of ln B(s) and the drift B·qE_r/T + dB/ds.
  - The calibration is done without electric field.
- **Markers.** Energies are drawn from a Maxwellian at `energy_eV` and then rescaled to
  the local target temperature of each species. The marker weights use the same local
  temperature.
- **Potential update.** A local linearised Boltzmann response is used (no Poisson
  solve). At each vertex:

      ΔΦ = ρ · T_h / (e² n) / 3,   with T_h = 2 T_i T_e / (T_i + T_e)

  For T_i = T_e this reduces to the former T/(e² n₀). The factor 1/3 is a relaxation
  factor.
- **Source tuning.** The particle source is iterated until the marker density matches
  the target n. It is not tuned in the outer buffer.
- **Outer buffer** (`edge_buffer_width`, as a width in r/a). The boundary at s = 1
  absorbs particles, so the marker density always drops to zero there. In the outer
  buffer, each flux-shell density is therefore blended with the target density:
  - the weight goes from 0 at the inner edge of the buffer to 1 at the boundary,
    following a cosine;
  - the weight is evaluated at the outer edge of each shell, so the outermost shell
    equals the target exactly and carries no charge.

  As a consequence, Φ hardly evolves in the buffer. E_r there is not a result and
  should be left out of the comparison.

## Settings of the first benchmark run (low statistics)

Differences from the example inputs in `self_consistent_ef.inp`:

| Key | Value | Note |
|---|---|---|
| `n_particles` | 500 | ion markers |
| `boole_honest_tracing(1:2)` | `.true., .false.` | ions honest, electrons random walk |
| `boole_recompute_D(1:2)` | `.false., .true.` | calibrate the electrons |
| `n_particles_diffusion` | 500 | electrons per surface in the calibration |
| `diffusion_tau_collisions` | 25 | tracing time per surface in local collision times |
| `edge_buffer_width` | 0.1 | in r/a, i.e. s > 0.81 |
| `n_source_updates`, `n_electric_potential_updates` | 2, 4 | |
| `time_step` | 1000 s | maximum tracing time per marker |

The grid is 30×30×30, uniform in s, with `sfc_s_min = 1e-4`. On one VSC-5 zen3 node
(128 threads), the calibration takes about 5 minutes and one ion pass about 2 minutes.

## How to run

1. Build the branch against GORILLA at the commit in `GORILLA_REF` (see `docs/building.md`).
2. Put the W7-X wout file into this folder.
3. Set the values from the table above in `self_consistent_ef.inp`.
4. Run `./gorilla_applets_main.x`.

The outputs `phi_elec_1d_<source>_<update>.dat`, `one_d_densities*_*.dat`,
`particle_source_<k>.dat` and `A_and_B_species2.dat` are written per iteration.
`A_and_B_species2.dat` holds the calibrated A and B per flux surface. A later run can
reuse it with `boole_recompute_D(2) = .false.`, but only with the same profiles and
the same grid.

## Status and known limitations

- The first low-statistics run is in progress. In the first update, 9974 of 10000
  electron walkers left the domain (mean exit time 15.7 s), while the ions left after
  1.07 s on average. The resulting electron excess (quasineutrality mismatch 0.72)
  has to be removed by the potential iterations.
- Calibrated B falls by about a factor of 60 from the core to the edge. The quadratic
  fit of ln B is within ±20–30% of the values at 500 electrons per surface, apart from
  a few single-surface outliers. It is about 20–30% too high for r/a < 0.3.
- D is calibrated without E_r. Once the field has converged, it should be calibrated
  once more with the field: a finite E_r can reduce electron transport considerably.
- The marker temperature is not adjusted toward the target. Collisions with the
  background are expected to relax the markers to the target T, but this has not been
  checked yet.
- The edge buffer is implemented only for `update_dimension = 1`.
- Statistics are low on purpose. Raise `n_particles` and `n_particles_diffusion`
  (2000–4000) near convergence.

## Fixes made on this branch along the way

- **B fit.** The random walk used a quadratic fit of B itself, which turned negative
  towards the edge. Walkers there got NaN positions and never left the domain. The fit
  is now of ln B.
- **Debug output.** When the pusher loses a particle it should not lose, it writes a
  debug dump. After a bad potential update, these dumps filled 33 GB in one run. Only
  the first 1000 are now written.
- **Profile coordinate.** The target density and the old analytic temperature model
  were evaluated in s instead of r/a = sqrt(s).
