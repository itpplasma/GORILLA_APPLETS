# Central marker-weight routine — plan and status

Branch: `central_marker_weights` (created from `scef_3d_updates`). Last updated 2026-10-17.

## Goal

One central `calc_particle_weights_and_jperp` used by all applets, instead of near-identical copies
(base version used by Boltzmann / divertor / anomalous transport, helical-core version, SCEF version).
The central routine uses the clean form

    w = f · J_x · J_y / (pdf_position · pdf_velocity_space)

with `f` the physical distribution, `J_x`, `J_y` the Jacobians of position and velocity space, and the
pdfs the sampling densities of the markers (`start%dist_position`, `dist_lambda`, `dist_energy` from
`marker_distribution_mod`).

## Decisions taken

- **Location:** the central routine replaces the base version in `SRC/utils_orbit_timestep_mod.f90`
  (same name, base callers unchanged; helical core and SCEF will `use` it).
- **Linear profiles:** one shape for all applets, `linear_profile_factor(x) = 1 − 0.9·x` (0.1 at the
  boundary) and `linear_profile_slope = −0.9`, in `marker_distribution_mod` next to `pdf_boltzmann`.
  Done in commit `aacc18a`. The module is not renamed (used on ~25 other branches).
- **Base normalisation:** shown analytically to equal the clean form for uniform position (box
  R × 2π × Z), uniform pitch on [−1, 1] and uniform energy on [0, 5T]: Boltzmann
  `w = n·V_box·(10/√π)·r·√(E/T)·e^(−(E+qΦ)/T)`, monoenergetic `w = n·V_box·r`. So the old prefactor
  (`set_weights`) and the `energy_eV·ev2erg` trick become unnecessary.
- **`J_x`** chosen by coordinate system: flux coordinates `|sqg1 + gsqg·z|` (√g), otherwise refined √g
  or `r`. All base/helical examples and tests use cylindrical coordinates, so they are unaffected.
- **`c%weight_factor`** (only used by the experimental `boole_preserve_energy_and_momentum_during_collisions`,
  off in all inputs): compute it directly from the inputs, `1/(n·V_box·(10/√π if Boltzmann)·amax)`, in the
  collision setup (`calc_collision_coefficients_for_all_tetrahedra`), instead of from the `set_weights` prefactor.
- **Point source:** degenerate position distribution evaluates as 1; the central routine multiplies `f` by
  the box volume (as helical core did), which reproduces the base prefactor.
- **SCEF specifics** (to do): `Φ = 0` in `f` via an explicit switch (optional argument, default `.true.`);
  Maxwell energy pdf (√E·e^(−E/T)) added to `marker_distribution_mod` so SCEF samples through `dist_energy`;
  the diffusion calibration sets its weights to `1/n_particles` itself and only calls a separate jperp
  routine (no calibration flag in the weight routine).
- **Out of scope:** field-line tracing (own weight arrays), delta-f in SCEF, the experimental
  energy/momentum-conserving collisions.

## Status

Done (commit "WIP: central marker-weight routine replacing the base version"):

- Central routine in `utils_orbit_timestep_mod.f90`, plus helper `calc_radial_coordinate`. Also fixes the
  implicit-SAVE bug of `integer :: species = 1` and drops the `pdf_energy_ev == 1` workaround.
- `init_velocity_space_distributions` split out of the shared `set_rest_of_start_type`
  (`utils_data_pre_and_post_processing_mod.f90`).
- Divertor (deterministic marker placement) initialises the uniform distributions for evaluation only
  (`init_marker_distributions` in `divertor_heat_loads_mod.f90`).
- Builds. **Not yet regression-tested.** `set_weights` and `c%weight_factor` still unchanged.

Remaining (estimate ~5–7 h):

1. Regression for Boltzmann, divertor, anomalous transport: run the pre-refactoring build (commit
   `aacc18a`) and the current build on small inputs with `OMP_NUM_THREADS=1` and a fixed seed, compare all
   output `.dat` files numerically (expect round-off differences only). Divertor may not run locally:
   `DATA/ASDEX/convexwall.dat` is missing in the example.
2. Move `c%weight_factor` into the collision setup; remove `set_weights` and divertor's copied prefactor.
3. Helical core: delete `calc_particle_weights_and_jperp_helical_core`, call the central routine and
   `adapt_weights_delta_f` at the call site; check with the CI test `helical_core`.
4. SCEF: Maxwell energy pdf, initialise the distributions in `calc_starting_conditions`, `Φ = 0` switch,
   split off the calibration, `J_x` already handled; regression with honest tracing for both species.
5. Cleanup (remove this document once done), final build, push, merge after `scef_3d_updates`.
