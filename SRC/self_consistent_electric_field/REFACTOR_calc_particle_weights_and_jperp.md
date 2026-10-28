# Plan — refactor SCEF `calc_particle_weights_and_jperp` in the helical-core style

Drafted 2026-08-27. Original scratch copy: `~/.claude/plans/giggly-beaming-unicorn.md`.

## Context

Three variants of `calc_particle_weights_and_jperp` currently live in the tree:

- base — `SRC/utils_orbit_timestep_mod.f90:86` (used by boltzmann / divertor / anomalous)
- helical-core — `SRC/helical_core/utils_helical_core_mod.f90:524`
- SCEF — `SRC/self_consistent_electric_field/utils_scef_particle_pushing_mod.f90:438`

The SCEF variant accumulates the weight in place (`weights%w = weights%w * …` scattered across ~40 lines), mixes physical distribution with sampling Jacobians, hard-overrides `phi_elec = 0` mid-routine, and carries commented-out dead code. The helical-core rewrite already assembles the weight in the mathematically clean form

    w = f · J_x · J_y / (pdf_position · pdf_velocity_space)

and delegates evaluation of the sampling PDFs to `marker_distribution_mod`. The goal is to port the SCEF routine to that same skeleton, preserve the SCEF-specific semantics that are load-bearing (fixed-source `Phi_elec = 0`, `boole_diffusion_coefficient` mode), extract the pieces the two routines duplicate into one shared module, and verify with an honest-tracing run that the physics has not been perturbed.

## Non-goals

- Not unifying the base variant into the same skeleton. That would touch boltzmann / divertor / anomalous and is out of scope. Only the three helpers are shared.
- Not enabling delta-f in SCEF. `boole_delta_f` and `adapt_weights_delta_f` (`utils_helical_core_mod.f90:471`) stay a helical-core-only path.
- Not changing the `1 − 0.9·x1` vs. `(1.1 − s)/1.1` mismatch between SCEF and the other two applets. If we want to reconcile them, do it separately with its own regression run.

## Step 1 — new shared helper module `SRC/utils_marker_weights_mod.f90`

Rationale: `utils_orbit_timestep_mod.f90` is named for timestep/orbit driving, and pulling three per-particle primitives into it would force `utils_scef_particle_pushing_mod` to newly `use` it and drag in unrelated symbols. A dedicated small module keeps the coupling narrow.

Contents (three pure helpers, no I/O, no OMP state):

- `pure function s_value_from_coordinates(ind_tetr, z_save) result(s_value)` —
  the `coord_system == 2` / `grid_kind /= 3` / else branch currently duplicated at `utils_orbit_timestep_mod.f90:122-137` and `utils_helical_core_mod.f90:558-572`. Returns `real(dp)`. Stops with the same error message if a stellarator is combined with cylindrical coordinates.

- `pure function spatial_jacobian_x(ind_tetr, z_save) result(J_x)` —
  dispatches on `in%boole_refined_sqrt_g` and on `coord_system`.
  - `coord_system == 2` (flux): `|sqg1 + gsqg·z_save|` — the SCEF form at `utils_scef_particle_pushing_mod.f90:455`.
  - refined sqrt_g in cylindrical: the ratio form at `utils_orbit_timestep_mod.f90:116-118` / `utils_helical_core_mod.f90:592-593`.
  - plain cylindrical: `r + tetra_physics(ind_tetr)%x1(1)` (`utils_helical_core_mod.f90:595`, `utils_orbit_timestep_mod.f90:120`).

- `pure function jperp_from_vperp(vperp, ind_tetr, z_save, species) result(jperp)` —
  the one-liner `m * vperp**2 * cm_over_e / (2*B) * (-1)` repeated verbatim at all three sites.

Only dependencies: `tetra_physics_mod`, `volume_integrals_and_sqrt_g_mod`, `supporting_functions_mod`, `gorilla_settings_mod`, `tetra_grid_settings_mod`, `gorilla_applets_types_mod`, `constants`. All three already imported by every consumer.

## Step 2 — rewrite `calc_particle_weights_and_jperp` in `utils_scef_particle_pushing_mod.f90`

Target shape (mirrors `utils_helical_core_mod.f90:524-630`, minus delta-f, plus SCEF hooks):

    subroutine calc_particle_weights_and_jperp(n, z_save, vpar, vperp, ind_tetr, &
                                               species, boole_diffusion_coefficient)
        ! diffusion-coefficient short-circuit: reuse the current semantics literally
        if (boole_diffusion_coefficient) then
            weights%w(n, species) = 1.0_dp / s%n_particles
            start%jperp(n, species) = jperp_from_vperp(vperp, ind_tetr, z_save, species)
            return
        endif

        ! sampling PDFs (from marker_distribution_mod)
        x_global = z_save + tetra_physics(ind_tetr)%x1
        pdf_position = evaluate_distribution_3d(start%dist_position, x_global)
        pdf_energy   = evaluate_distribution_1d(start%dist_energy,   energy_ev)
        pdf_lambda   = evaluate_distribution_1d(start%dist_lambda,   start%pitch(n,species))
        pdf_gyro     = 1.0_dp / (2.0_dp * pi)

        ! Jacobians
        J_x = spatial_jacobian_x(ind_tetr, z_save)
        ! J_y: identical to helical_core lines 617-620 (Boltzmann case) / 606-608 (monoenergetic)

        ! physical distribution f
        if (in%boole_monoenergetic) then
            f = density * pdf_gyro / (2.0_dp * start%v0(species))
        else
            ! SCEF-SPECIFIC: fixed-source semantics -> Phi_elec = 0
            Phi_elec = 0.0_dp
            T_local = in%energy_eV * ev2erg
            if (in%boole_linear_temperature_simulation) T_local = T_local * (1.0_dp - 0.9_dp*x_global(1))
            density_local = in%density
            if (in%boole_linear_density_simulation) density_local = density_local * (1.0_dp - 0.9_dp*x_global(1))
            f = pdf_boltzmann(density_local, T_local, m, energy_erg, q, Phi_elec)
        endif

        weights%w(n, species) = f * J_x * J_y / (pdf_position * pdf_lambda * pdf_energy * pdf_gyro)
        start%jperp(n, species) = jperp_from_vperp(vperp, ind_tetr, z_save, species)
    end subroutine

Preserved SCEF-specific semantics, all consciously kept:

1. **`Phi_elec = 0`** — matches the intent of the current line 476 comment ("fixed sources"). Documented with a one-line comment pointing at the fact that if the electric potential feedback is desired in the weight, this line becomes `Phi_elec = tetra_physics(ind_tetr)%Phi1 + sum(...)` (as helical-core does at line 611).

2. **`boole_diffusion_coefficient` short-circuit** — the current line 463 override is preserved verbatim. It is placed early to avoid touching any of the Boltzmann/PDF math; the diffusion-coefficient path never needed it.

3. **`1 − 0.9·x1` linear profile** — kept exactly as at lines 459 and 489-492. Deliberately not migrating to `(1.1 − s)/1.1` in this refactor — see Non-goals.

4. **Mandatory `species` argument** — retained (matches current signature at line 438). Not making it optional, because every SCEF call site (`utils_scef_particle_pushing_mod.f90:216`, `utils_scef_diffusion_mod.f90:258`) already passes it.

Dead code to drop:
- lines 470-473 (commented old formula)
- lines 480-482 (commented `epsilon_max` branch)
- lines 501-503 (commented print)
- the local `temperature` reassignments guarded by `if (boole_diffusion_coefficient)` that live inside an `else` branch keyed on `.not. boole_diffusion_coefficient` (unreachable — lines 467, 485, 490).

## Step 3 — SCEF start-distribution initialisation

Missing prerequisite: SCEF currently allocates `start%pitch` and `start%energy` (`utils_scef_particle_init_mod.f90:218-219`) and samples them directly (`:338-347`) without touching `start%dist_position`, `start%dist_energy`, `start%dist_lambda`. If the rewritten routine calls `evaluate_distribution_*` on those, they must be initialised or it will trip the "uninitialized distribution" error at `marker_distribution_mod.f90:131,276`.

Add three init calls in `utils_scef_particle_init_mod.f90`, ideally in `calc_starting_conditions` before the first sample (so re-initialising a new source doesn't leak state — check `cleanup_distribution_*` is called first if the dist is already `initialized`):

- `dist_position`: mirror the existing helical-core init pattern at `utils_data_pre_and_post_processing_mod.f90:211-233`. `pdf_flat` over the flux-coordinate cube `[sfc_s_min, 1] × [0, 2π] × [0, 2π/n_field_periods]`, matching the sampling box actually used at `utils_scef_particle_init_mod.f90:338-339` and its neighbours. (Re-read `calc_starting_conditions` in that file to lift the exact bounds — do not re-derive.)
- `dist_lambda`: `pdf_flat` on `[-1, 1]`, matching line 338 (`start%pitch = 2*rand - 1`).
- `dist_energy`: this is the delicate one. Line 339 sets `start%energy = in%energy_eV` (monoenergetic default), and the commented line 341 / active line 347 (`start%energy = temperature * radial_transport_energies`) is a *weighted* energy sampling — not flat. Two acceptable resolutions:
  - **Preferred:** `pdf_flat` on `[energy_eV, energy_eV]` for the monoenergetic case, `pdf_flat` on `[0, 5*energy_eV]` for the Boltzmann case (same choice as `utils_data_pre_and_post_processing_mod.f90:274-280`). This is only correct if SCEF also switches its actual sampling to flat. If SCEF keeps `radial_transport_energies`, the sampling PDF is not flat and using `pdf_flat` here will bake a wrong Jacobian into the weight.
  - **Alternative:** in the monoenergetic path only, skip `evaluate_distribution_1d(start%dist_energy, ...)` in the rewritten routine (helical-core's line 605 monoenergetic branch already does not use it). For the non-monoenergetic path, resolve the sampling-PDF question explicitly — either flatten the SCEF energy sampling or add the appropriate `pdf_*` to `marker_distribution_mod`.

  **Decision needed at implementation time**, not now: read `radial_transport_energies` in `utils_scef_particle_init_mod.f90` around lines 340-350 and decide whether we can flatten the SCEF energy sampling without changing physics. If unsure, take the alternative path (skip `pdf_energy` for the monoenergetic run — which is what the test config uses anyway).

## Step 4 — switch on honest tracing for verification

The current input `EXAMPLES/self_consistent_electric_field/self_consistent_ef.inp:86-87` has

    boole_honest_tracing(1) = .false. ,
    boole_honest_tracing(2) = .false. ,

so `calc_particle_weights_and_jperp` is not exercised on the pushing path — only via the random-walk substitute in `utils_scef_diffusion_mod.f90:258`. **Any test of this refactor must set both to `.true.`** so `parallelised_particle_pushing` actually calls the rewritten routine at `utils_scef_particle_pushing_mod.f90:216`.

## Verification

Run from `EXAMPLES/self_consistent_electric_field/` with a scaled-down config:

1. **Regression config** — copy `self_consistent_ef.inp`, set `boole_honest_tracing = .true., .true.`, `n_particles = 1000`, `n_electric_potential_updates = 1`, `n_source_updates = 1`, `seed_option = 2`. Build with `make CONFIG=Release`. Run before the refactor, save `1_prism_moments.dat`, `2_prism_moments.dat`, `phi_elec_1d_1_1.dat`, `A_and_B*.dat`, and the two `sum(weights)` printouts at `SRC/self_consistent_electric_field/self_consistent_electric_field_mod.f90:89-90`. Run after the refactor. Diff.

   Expected: differences are floating-point noise from arithmetic reordering (relative ~1e-12) if the refactor preserves math, or a clear systematic sign/scale change if a factor was dropped or added. `sum(weights)/(N·density·V)` and `sum(volumes·moments)/(V·density)` must both stay ≈ 1.

2. **Diffusion-coefficient path** — a second small run with `boole_recompute_D = .true., .true.`. Compare `A_and_B_species1.dat` and `A_and_B_species2.dat` byte-for-byte (the short-circuit path is preserved verbatim, so identical output is expected).

3. **Cross-check against helical-core skeleton** — pick one tetrahedron and one particle, evaluate the rewritten SCEF routine with `Phi_elec = 0` overridden to the helical-core evaluation of `Phi_elec`, no linear profiles, monoenergetic mode; separately evaluate `calc_particle_weights_and_jperp_helical_core` on the same inputs. They should agree to FP tolerance. Doing this in a small standalone Fortran driver (or a temporary `program` in a scratch file) is cleaner than parsing outputs.

4. **Compile the whole tree** — `make` from repo root. All existing callers of the base `calc_particle_weights_and_jperp` (`boltzmann_mod`, `divertor_heat_loads_mod`, anomalous transport) and of the helical-core variant must still build; the shared helpers change no signatures on their side.

5. **Full end-to-end example** — restore the user's production config (`n_particles = 10 000`, `n_electric_potential_updates = 10`, `n_source_updates = 100`, honest tracing on for at least species 1) and run to completion. Compare the final `phi_elec_1d_100_10.dat` against a saved pre-refactor baseline. Systematic offset in the boundary potential (visible in the plotting script's second figure) is a red flag.

## Files touched

- **new**: `SRC/utils_marker_weights_mod.f90`
- **rewrite**: `SRC/self_consistent_electric_field/utils_scef_particle_pushing_mod.f90` (subroutine at line 438 only)
- **augment**: `SRC/self_consistent_electric_field/utils_scef_particle_init_mod.f90` (add `init_distribution_*` calls in `calc_starting_conditions`)
- **build**: `CMakeLists.txt` — add the new module file
- **test-only**: `EXAMPLES/self_consistent_electric_field/self_consistent_ef.inp` (temporary honest-tracing / small-n edits for regression runs — not committed)

## Deferred / follow-ups

- Migrating the base `calc_particle_weights_and_jperp` in `utils_orbit_timestep_mod.f90` to use the same helpers (safe: same math, and it already contains the `s_value` block verbatim).
- Reconciling the `1 − 0.9·x1` vs. `(1.1 − s)/1.1` linear-profile inconsistency. Needs its own physics decision.
- Optional delta-f in SCEF via `adapt_weights_delta_f`.
