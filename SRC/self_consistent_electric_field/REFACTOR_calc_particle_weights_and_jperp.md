# Central marker-weight routine — status and remaining work

Branch: `central_marker_weights` (from `scef_3d_updates`). Pre-refactor baseline: commit `aacc18a`
("Use one linear profile shape (0.1 at the boundary) in all applets"). The last commit that
reproduces the old numbers bit-exactly on every applet is `aacc18a`.

## Goal

One central `calc_particle_weights_and_jperp` used by all applets, instead of near-identical copies
(Boltzmann / divertor / anomalous base, helical core, SCEF). Form:

    w = f · J_x · J_y / (pdf_position · pdf_velocity_space)

with `f` the physical distribution, `J_x`/`J_y` the position/velocity-space Jacobians, and the pdfs
the sampling densities of the markers, evaluated from `start%dist_position`, `start%dist_lambda` and
`start%dist_energy` (`marker_distribution_mod`).

## What is already done on this branch

- **Central routine** in [`SRC/utils_orbit_timestep_mod.f90`](../utils_orbit_timestep_mod.f90)
  — the only `calc_particle_weights_and_jperp` left in the tree. Reads `density` from
  `start%density(n, species)` so per-marker physical density (e.g. SCEF particle source) is a
  one-line change at the call site. Picks `J_x` by coordinate system (flux: `|sqg1 + gsqg·z|`;
  cylindrical: refined √g or `r`). Has an optional `boole_include_potential_in` (default `.true.`)
  that SCEF passes as `.false.` to drop `qΦ` from the Boltzmann `f`.
- **`calc_particle_jperp`** split out of the central routine so callers that assign weights
  themselves (SCEF's diffusion calibration) can still set `start%jperp`.
- **`calc_radial_coordinate`** helper: flux `s` in flux coordinates, normalised poloidal flux in
  cylindrical/axisymmetric cases, hard stop in cylindrical/stellarator.
- **`start%density(n, species)`** added to the `start_t` type and allocated by every
  `allocate_start_type` (base, divertor, SCEF, anomalous local). Default-filled with `in%density`;
  SCEF's `apply_particle_source_weights` overwrites it with `particle_source(s_marker)`.
- **`init_velocity_space_distributions`** split out of `set_rest_of_start_type` so applets that
  place markers deterministically (divertor, SCEF's diffusion calibration) can initialise the
  distributions for evaluation only.
- **Linear profiles**: single shape across all applets,
  `linear_profile_factor(x) = 1 − 0.9·x` (0.1 at the boundary), slope `−0.9`, in
  `marker_distribution_mod` next to `pdf_boltzmann`.
- **Maxwell energy pdf** `pdf_maxwell(E) ∝ √E · e^(−E/T)` with module-level
  `maxwell_temperature_eV` in `marker_distribution_mod`. SCEF samples energies through
  `start%dist_energy` initialised from it.
- **Base normalisation**: `set_weights` deleted. `c%weight_factor` (only consumed by the
  experimental energy/momentum-conserving collisions) now computed directly from the inputs
  inside `calc_collision_coefficients_for_all_tetrahedra`, gated on
  `in%boole_preserve_energy_and_momentum_during_collisions` so helical_core and divertor can call
  the collision setup before `g%amin/amax/cmin/cmax` are populated.
- **Divertor** (deterministic markers): the copied weight prefactor in `set_start_type` is gone;
  `init_marker_distributions` seeds the uniform position/pitch/energy distributions used only for
  evaluation.
- **Helical core**: local `calc_particle_weights_and_jperp_helical_core` deleted; the single call
  site now does `call calc_particle_weights_and_jperp(...)` followed by
  `if (boole_delta_f) call adapt_weights_delta_f(...)` + `weights%original = weights%w`.
- **SCEF**: local `calc_particle_weights_and_jperp` and `set_weight` deleted. The main orbit
  timestep calls the central routine with `boole_include_potential_in = .false.`; the diffusion
  calibration path sets `weights%w = 1.0_dp/s%n_particles` and calls `calc_particle_jperp`. The
  `boole_diffusion_coefficient` flag on the weight routine is gone. `calc_density_via_random_walk`
  calls the central routine (same `Φ = 0` switch). `init_marker_distributions(temperature_eV)` is
  called once from `calc_starting_conditions` (at `in%energy_eV`) and once from
  `calc_diffusion_coefficients` (at `s%temperature`); it cleans up and re-initialises
  `dist_energy` when the temperature changes between the two.
- **Also fixed as part of the refactor**: implicit-SAVE bug of `integer :: species = 1` in
  `update_local_tetr_moments`; dropped the `pdf_energy_ev == 1` workaround (unused by the central
  routine).

Builds. All six CI tests pass. The helical_core CI run is **bit-exact** against the pre-refactor
baseline on all eight output `.dat` files (`exit_data`, `prism_volumes`, …):

    cd BUILD/TESTS/helical_core_work && sha256sum *.dat   # matches baseline from aacc18a

## What is NOT yet done

1. **Numerical regression on Boltzmann / divertor / anomalous_transport / SCEF.** The helical_core
   CI test confirms the central routine faithfully reproduces that applet, but the other applets
   have no in-tree test. The weight derivation shows round-off differences only (the central
   normalisation equals the base one analytically — see "Decisions taken" below), so running
   small deterministic examples at `OMP_NUM_THREADS=1` with a fixed seed and comparing all output
   `.dat` files numerically should show round-off only.

   Suggested small inputs: the examples in `EXAMPLES/<applet>/` are the obvious starting point.
   `EXAMPLES/divertor/` is likely missing `DATA/ASDEX/convexwall.dat` (noted in the pre-refactor
   version of this doc); if so, run divertor on a different example or skip.

   Procedure for each applet:

       # baseline
       git checkout aacc18a && make -j
       cd <work_dir> && OMP_NUM_THREADS=1 ./gorilla_applets_main.x
       cp -r <work_dir> /tmp/baseline_<applet>

       # current branch
       git checkout central_marker_weights && make -j
       cd <work_dir> && OMP_NUM_THREADS=1 ./gorilla_applets_main.x

       # numerical diff of every .dat file
       for f in /tmp/baseline_<applet>/*.dat; do
           numdiff -r 1e-12 "$f" "<work_dir>/$(basename $f)" | head -20
       done

   SCEF needs honest tracing for both species to exercise the Φ=0 and source-density paths.

2. **Remove this document and merge.** Once the regression comparisons look clean, delete this
   file, confirm `make test` still passes, push `central_marker_weights`, and merge after
   `scef_3d_updates` (the parent branch).

## Design decisions worth preserving

- **Location:** central routine lives in `utils_orbit_timestep_mod.f90` next to the shared push
  utilities (`identify_particles_entering_annulus`, `update_local_tetr_moments`,
  `initialize_constants_of_motion`, `compute_radial_fluxes`, `calc_particle_jperp`).
- **Base normalisation = clean form, analytically.** Shown that for uniform position
  (box `R × 2π × Z`), uniform pitch on `[−1, 1]` and uniform energy on `[0, 5T]`:
  Boltzmann `w = n·V_box·(10/√π)·r·√(E/T)·e^{−(E+qΦ)/T}`, monoenergetic `w = n·V_box·r`. The old
  `set_weights` prefactor and the `energy_eV·ev2erg` trick are therefore unnecessary. For SCEF
  with Φ = 0 and non-uniform source: `w = source(s)·V_box_flux·√g` for both monoenergetic and
  Boltzmann (the velocity-space factors cancel when the sampling pdf matches the physical
  Maxwellian).
- **`c%weight_factor`**: only consumer is `boole_preserve_energy_and_momentum_during_collisions`
  (off in all shipped inputs), so compute it directly from the inputs in the collision setup and
  gate on that switch. Decouples it from `weights%w` and from `g%amax` being populated.
- **Point source** (`in%boole_point_source`): degenerate position distribution evaluates as 1;
  the central routine multiplies `f` by the box volume, reproducing the old prefactor behaviour.
- **Delta-f** (`in%boole_delta_f`, only helical_core so far): adapt is called at the helical_core
  call site after the central routine, and `weights%original = weights%w` is set from the adapted
  weight, matching the old behaviour where the old helical-core routine finished with
  `adapt_weights_delta_f` before `weights%original` was set outside.
- **Out of scope for this branch:** field-line tracing (own weight arrays), delta-f in SCEF, the
  experimental energy/momentum-conserving collisions.
