# TODO: GORILLA Kramers-Moyal D11 Benchmark

Branch: `transport-km-benchmark` (stacked on `multi-species-collision-background` PR #25)
Benchmark target: GORILLA vs NEO-2-QL on AUG 39461 neon discharge

## Current status

GORILLA D11 is 3-5x below NEO-2-QL across the full radius (see comparison plot).
The radial trend matches (D11 increases with s) but the absolute magnitude is off.
Per-surface collisionality matching and 5-species background are implemented.

Comparison plot: https://litter.catbox.moe/e8vyqc.png

## Done

- [x] KM unit tests (1D random walk, fit_transport_coefficients)
- [x] KM module (`kramers_moyal_transport_mod.f90`)
- [x] i_option=14 driver (`km_benchmark_mod.f90`, `km_benchmark_settings_mod.f90`)
- [x] N-species collision background (PR #25, `set_custom_background`)
- [x] Per-surface collision profile (reads NEO-2 densities + collisionality-matched energy)
- [x] Comparison script (`convert_and_compare_d11.py`)
- [x] Run script (`run_km_benchmark_39461.sh`)
- [x] Benchmark run: 5000 particles, 20 surfaces, 5 species, per-surface energy
- [x] Bugs fixed: start%t allocation order, s%temperature init, mono-energetic energy flag,
      aliasing in set_custom_background call

## Remaining: investigate 3-5x D11 offset

The remaining offset is systematic and must be understood before publishing. Candidates:

### 1. Coulomb logarithm mismatch

GORILLA computes ln(Lambda) via `lambda_alpha_beta` in `collis_ions_mod.f90`.
NEO-2 uses its own Coulomb logarithm from `libneo` (`collision_freqs.f90`).
The formulas differ. At T=500 eV, n=7e13: GORILLA gives ln(Lambda)~14.
Need to verify NEO-2 gives the same value. A factor of 2 in ln(Lambda)
gives a factor of 2 in D11.

Action: compare GORILLA and NEO-2 Coulomb logarithm at same parameters.

### 2. Collision frequency definition

NEO-2's `collpar_spec` may include geometric factors (bounce averaging,
trapped fraction) that the simple Spitzer formula in GORILLA does not.
The implied mono-energetic energy from NEO-2 collpar varies 87-999 eV
across the radius, which I match empirically. But the exact relationship
between NEO-2 collpar and GORILLA efcolf is not verified analytically.

Action: read NEO-2 source code for collpar definition. Or run NEO-2 at
a single surface with known parameters and compare collision frequency.

### 3. Energy convolution

NEO-2 D11_AX is the FULL energy-integrated transport coefficient
(velocity-space integral of mono-energetic D11 weighted by the
distribution function). GORILLA computes mono-energetic D11 at a single
energy. The energy convolution typically gives a factor of 1-3x enhancement.

Action: run GORILLA at multiple energies per surface to scan nu_star,
compute the energy convolution integral, compare with NEO-2.

### 4. Tracing time estimate

The hardcoded tau_c_ei = 1.7e-4 s may be too short at some surfaces
(especially at high energy where collision time is longer). Although
convergence was verified at s=0.5, other surfaces may need longer tracing.

Action: run convergence test at multiple surfaces with varying total_time.

### 5. Orbit dynamics validation

GORILLA's polynomial pusher (ipusher=2, poly_order=2) may not capture
electron banana orbits accurately enough. The D11/Dp0 ratio from GORILLA
(~0.1 at banana-plateau transition) is lower than expected (~0.5).

Action: compare GORILLA banana orbit width with analytical estimate.
Run single trapped particle orbit and measure radial excursion.

## Files

| File | Purpose |
|------|---------|
| `SRC/kramers_moyal_transport_mod.f90` | KM D11 profile, per-surface collision support |
| `SRC/km_benchmark_mod.f90` | i_option=14 driver |
| `SRC/km_benchmark_settings_mod.f90` | settings + collision profile reader |
| `SRC/transport_benchmark_utils_mod.f90` | fit_transport_coefficients |
| `SRC/utils_data_pre_and_post_processing_mod.f90` | set_custom_background, set_c, collision init |
| `SRC/collis_ions_mod.f90` | Coulomb collision operator (N-species) |
| `TESTS/test_kramers_moyal.f90` | KM unit tests |
| `TESTS/test_multi_species_collision.f90` | N-species collision test |
| `INPUT/km_benchmark.inp` | benchmark input (5 species, profile file) |
| `INPUT/neo2_collision_profile_39461.dat` | per-surface density+energy from NEO-2 |
| `scripts/run_km_benchmark_39461.sh` | end-to-end benchmark runner |
| `scripts/convert_and_compare_d11.py` | D11 conversion and comparison plot |

## Reference data

| Data | Location |
|------|----------|
| NEO-2 D11 | `$DATA/AUG/NEO-2/39461/neoart_benchmark_neon_discharge/neon_discharge_out.h5` |
| NEOART D11 | `$DATA/AUG/NEO-ART/39461/neoart_benchmark_neon_discharge/neon_discharge_neoart.nc` |
| EQDSK | `$DATA/AUG/EQDSK/39461/eqdsk_39461_5.38s` |
| Convex wall | `$DATA/AUG/BOOZER/39461/convexwall.dat` |
| Benchmark tools | `$DATA/AUG/NEO-2/neo2_neoart_benchmark_tools/` |
| GORILLA output | `$DATA/AUG/GORILLA/39461/neoart_benchmark_neon_discharge/` |

## GitHub issues

- itpplasma/GORILLA_APPLETS#23 -- N-species collision support (PR #25)
- itpplasma/GORILLA_APPLETS#24 -- read collision profiles from HDF5
- itpplasma/libneo#271 -- expose collision modules for GORILLA
- itpplasma/libneo#272 -- Zeff and collision frequency utilities
