# Homma–Hatayama thermal force in GORILLA_APPLETS

**Date:** 2026-11-05
**Branch:** `feat/homma-thermal-force` (GORILLA_APPLETS)
**Benchmark:** `thermal-force-benchmark` (shared with kin6d), cases B01–B07 first

**Scope:** Add the Fokker–Planck test-ion collision operator of Homma &
Hatayama (2014, Contrib. Plasma Phys. 54, 394) to GORILLA_APPLETS. The
background ions follow a Maxwellian distorted by their classical heat flux.
The operator therefore contains the thermal force from parallel and
perpendicular temperature gradients. The thermal force enters the
guiding-centre motion through a gyrophase-resolved collision kick. A new
applet runs the benchmark cases. Existing applets and `stost` are unchanged.

## 1. Independence rule

The benchmark compares independent implementations. The GORILLA side
therefore uses no kin6d results, kernels or systematic findings, and none
from the benchmark repository's full-orbit Monte Carlo. Specifically:

- The coefficient algebra is derived here with Wolfram Language (`math`).
  kin6d and the benchmark Monte Carlo both use FortSym, so we do not.
- The oracle values in the GORILLA unit tests come from our own derivation,
  not from `reference/` or `results/` of the benchmark repository.
- Numerical choices follow from our own convergence scans. One example is
  whether an implicit friction step is needed.
- The benchmark's case files, `docs/conventions.md` (SI units, Itô form,
  `2D = BB^T`) and its physical parameters are used, since they define the
  common problem.

We compare with the other codes only at the end, through the benchmark's
`scripts/compare.py`.

## 2. Model

Test ion `a` (mass `m_a`, charge `q_a`) in background ions `b`
(`m_b`, `q_b`, `n_b`, `T_b(x)`, flow `u_b`). Let `w = v − u_b` and
`x̃ = w / sqrt(2 T_b/m_b)`. The background distribution is Homma's Eq. (1):

    f_b = f_M [1 + (2 m_b/(5 n_b T_b²)) (x̃² − 5/2) q·w]

The heat flux is the classical closure printed below Eq. (1). We implement
all three terms; the perpendicular term is not dropped:

    q = −κ_∥ ∇_∥T_b + κ_∧ b×∇_⊥T_b − κ_⊥ ∇_⊥T_b
    κ_∥ = 3.9 n_b T_b τ_b/m_b,  κ_∧ = 5 n_b T_b/(2 m_b Ω_b),  κ_⊥ = n_b T_b/(m_b Ω_b² τ_b)
    τ_b = 12 π^{3/2} ε0² m_b^{1/2} T_b^{3/2}/(n_b q_b⁴ ln Λ)   (as printed)

At the B07 parameters `Ω_b τ_b = 471`, so `κ_⊥/κ_∧ = 8.5e-4`.

The Fokker–Planck operator is in Itô form with `dv = A dt + B dW` and
`2D = B Bᵀ`. In terms of the Rosenbluth potentials,

    A = Γ (1 + m_a/m_b) ∇_v H,   D = (Γ/2) ∇_v∇_v G,   Γ = q_a² q_b² ln Λ/(4π ε0² m_a²)
    H = ∫ f_b(v')/|v−v'| d³v',  G = ∫ f_b(v')|v−v'| d³v'

Splitting into the Maxwellian part and the part linear in `q` gives
`A = A^(M) + A^(T)` and `D = D^(M) + D^(T)`. These are Homma's Eqs. (2), (5), (6).

## 3. Algebra (step 1)

`DERIVATIONS/homma_hatayama_proofs.wl` (run by hand, about 3 min) does items 1-3 below; `DERIVATIONS/homma_hatayama.wl` (about 6 s, rerun by CTest) does item 4 from the proved kernel functions in `homma_hatayama_common.wl`. Together they:

1. Computes `H` and `G` of Eq. (1) directly. The Maxwellian part is the
   ℓ = 0 Legendre component and the `q` part is ℓ = 1, each a radial integral.
2. Forms `A` and `D` from them and proves that they equal the printed
   Eqs. (2), (5), (6). This includes the prefactors and the factor
   `sqrt(m_b/(2T_b))` in `F^∇T`, which is ambiguous in the scanned PDF.
3. Derives the small-`x̃` Taylor series. The closed forms cancel
   catastrophically there (`Φ/x̃⁵`, `Φ/x̃⁷`), so the kernel switches to the
   series below a threshold that the script chooses. The script writes the
   coefficients as a Fortran include file. A test checks that the include
   file is current.
4. Writes high-precision oracle values (30 digits) for the unit tests:
   - coefficients at `x̃ = 0, 0.1, 0.5, 1, 2, 5, 10` for the B01/B02 gradient
     variants, with the formulas proved in item 2;
   - force moments over test distributions (B03);
   - the gyrophase and Maxwellian averages needed for the B07 analytic
     targets.

## 4. Collision kernel: `SRC/homma_collision_mod.f90`

This is a standalone module. It does not depend on the GORILLA state, the
grid or `constants_mod`, so it can be unit-tested directly.

- **Unit system:** one setup routine takes the physical constants in the
  caller's unit system and precomputes the prefactors. The kernels work in
  `x̃` and the dimensionless heat flux `ĝ = m_b v_tb q/(n_b T_b²)`, where
  `v_tb = sqrt(2T_b/m_b)`. Tests use CODATA SI. GORILLA runs use its CGS
  `constants_mod`; the deviation (~2e-5 relative, from 5-digit constants)
  is declared in result headers.
- `homma_heat_flux(grad_T, b, background) → q`
- `homma_coefficients(v, background, A, D)` returns the Cartesian `A` and
  `D`, Maxwellian plus thermal parts, with an optional split output for
  diagnostics. The sum over several background species is built in;
  electrons are deferred.
- `homma_kick(v, background, dt, xi, dv, ierr)` performs the Euler–Maruyama
  step `dv = A dt + L xi sqrt(dt)`, where `L Lᵀ = 2D` is the Cholesky factor
  of the total `D` and `xi` is a standard normal 3-vector. Splitting off
  `D^(T)` is not possible because it is not positive definite on its own.
  `ierr ≠ 0` reports a non-positive `D` or a negative `f_b` region. Such
  failures are counted and reported, never silently repaired.

Only Euler–Maruyama, as in Homma, is implemented at first. A
higher-order friction step is added only if our own time-step scan (B08)
shows a bias at practical `Δt`.

## 5. Guiding-centre coupling

GORILLA has no gyrophase. A collision step at guiding centre `X` with
`(v_∥, v_⊥)` proceeds as follows:

1. Draw a gyrophase `φ ~ U[0, 2π)` and build
   `v = v_∥ b + v_⊥ (cos φ e₁ + sin φ e₂)`, where `(e₁, e₂, b)` is a local
   orthonormal frame.
2. Evaluate the background at:
   - `collision_background_at = 'gyro_orbit'` (default): the particle
     position `x = X + ρ`, with `ρ = b×v/Ω_a`. To first order this is
     `T_b(X) + ρ·∇T_b`, which is exact for the linear benchmark profiles.
     This is the consistent guiding-centre collision operator: the collision
     operator is evaluated at the particle and gyro-averaged by the random
     phase.
   - `'guiding_centre'` (diagnostic): `X`. This isolates the gyro-averaged
     Eq. (2) thermal force. It drops a drift that is first order in `∇T`
     and as large as the thermal drift itself for low charge states.
3. Apply `homma_kick` to obtain `v' = v + Δv`. Then set
   `v_∥' = v'·b` and `v_⊥' = |v' − v_∥' b|`.
4. Shift the guiding centre by `ΔX = −b×Δv/Ω_a`, keeping the particle
   position fixed. In GORILLA this goes through `displace_by_straight_line`,
   as in the anomalous-transport applet. That routine assumes cylindrical
   coordinates, so version 1 supports `coord_system = 1` only.

The expected value of the step contains the following. Each item is derived
in the Wolfram script and checked by the statistical kick test:

- `⟨Δv_∥⟩ = (⟨A^(M)⟩_φ + ⟨A^(T)⟩_φ)·b Δt`
- `⟨ΔX⟩_th = ⟨F^∇T⟩_φ,⊥ × b/(q_a B) Δt`, the F×B drift of the
  gyro-averaged thermal force
- in `gyro_orbit` mode only, `⟨ΔX⟩_flr = (∂ν/∂T_b) v_⊥²/(2Ω_a²) ∇_⊥T_b Δt`,
  where `A^(M) = −ν(|v|, T_b) v`. Averaging over φ in
  `−b×⟨A^(M)(T_b(X)+ρ·∇T_b, v)⟩/Ω_a` gives this in two lines.
- zero-mean fluctuations `−b×(L xi)/Ω_a sqrt(Δt)`, which give classical
  collisional diffusion without extra modelling.

The deterministic drift formulas are used only as test expectations. The
code always executes the random-phase kick.

## 6. Applet `i_option = 16`: `SRC/thermal_force/`

It reads its own input file, `thermal_force.inp`, with group
`&thermal_force_nml`. Names carry units, in GORILLA's CGS + eV convention.

| entry | meaning |
|---|---|
| `driver` | `'velocity'` (no spatial motion: B04, B05/B06 short time), `'slab'` (uniform `B`, exact streaming `X += v_∥ b Δt`), `'tokamak'` (GORILLA orbit, `grid_kind = 5`) |
| `test_mass_amu`, `test_charge_number` | test ion; sets `particle_mass`, `particle_charge` and `cm_over_e` together |
| `background_mass_amu`, `background_charge_number` | one background ion species, no electrons |
| `density_cm3`, `temperature_ev`, `flow_cm_s(3)` | background at the reference point |
| `grad_temperature_ev_cm(3)` | constant gradient; slab: `(x, y, z)`; tokamak: physical `(R, φ, Z)` components |
| `coulomb_log` | fixed `ln Λ` |
| `heat_flux_model` | `'homma2014'` or `'none'` (Maxwellian operator only) |
| `collision_background_at` | `'gyro_orbit'` (default) or `'guiding_centre'` |
| `test_temperature_ev`, `test_drift_cm_s(3)`, `test_temperature_par_ev` | initial test distribution (B04 uses a drifting bi-Maxwellian) |
| `b_gauss` | slab field strength, `b = e_z` |
| `start_x_cm(3)` | start position |
| `dt_coll_over_tau_s` | fixed collision step, `τ_s` at `v = sqrt(3T_a/m_a)` (conventions.md) |
| `t_end_gyro`, `n_output` | run length in `Ω_a t`, number of output times |
| `n_particles`, `seed` | Monte Carlo size, base seed |
| `boole_control_run` | also run `∇T = 0` with identical random numbers |

**Random numbers.** Each particle has its own generator (xoshiro256** with
splitmix64 seeding from `seed` and the particle index). Results are
therefore independent of OpenMP scheduling, and gradient and control runs
use identical streams. Normals come from Box–Muller.

**Output.** The applet writes time series of `⟨X⟩, ⟨Y⟩, ⟨v_∥⟩, ⟨v⟩` and the
velocity moments of B04, with standard errors and the counts of kernel
failures. A small script converts this to `results/gorilla/<case>.csv` in
the benchmark format, with header metadata (commit, case-file sha256,
deviations, step and particle count).

**Tokamak driver.** It uses `grid_kind = 5` with `R0 = 1e5 cm`, `a = 100 cm`,
`B0 = 1e4 G` and `q0 = 10`. The start is on the outboard midplane, with
`e_z ↔ e_φ`, `e_x ↔ e_R`, `e_y ↔ −e_Z`. The background is
`T_b(R) = T0 + g (R − R_s)`, evaluated analytically, not from the
per-tetrahedron table. The curvature and ∇B drift is
`2T/(Z e B R0) = 0.1/Z m/s` along Z, perpendicular to the screening drift.
A control run subtracts it.

## 7. Tests (CTest, `TESTS/`)

| test | content |
|---|---|
| `homma_derivation_current` | regenerates the Wolfram include file and compares it (only when `math` is available) |
| `test_homma_coefficients.f90` | B01–B03 against our Wolfram oracle (rel. 1e-10 where far from zeros). Also: `∇T → 0` limit, symmetry of `D`, finite-difference identities (`A` is a gradient, `D` is a Hessian), continuity at the series/closed-form switch, positivity of `D` over the benchmark range |
| `test_stost_cartesian.f90` | B01 for the existing `coleff`/`stost` through the Itô conversion `D_∥ = v0³ dpp`, `D_⊥ = v0³ p² dhh`, `A = (v0² fpeff − 2D_⊥/|v|) v̂` |
| `test_homma_kick.f90` | sample mean and covariance of `Δv` equal `A Δt` and `2D Δt` (fixed seed, tolerance from the standard error). Random-phase averages of `ΔX` equal the drifts in §5 |
| `test_thermal_force_applet.py` | small slab run: drift sign, ratio between 300 and 100 eV/m, control subtraction |

## 8. Order of work

1. Wolfram derivation and oracle values (§3)
2. Kernel module and `test_homma_coefficients`, `test_stost_cartesian` (B01–B03)
3. `homma_kick` and `test_homma_kick`
4. Applet skeleton, `velocity` driver (B04–B06)
5. `slab` driver (B07), then the time-step and particle-number scan (B08)
6. `tokamak` driver (B05–B07 in GORILLA geometry; start of B12)
7. Result CSVs, then `compare.py` in the benchmark repository

## 9. Review of the benchmark team's suggestions (2026-11-05)

| item | verdict |
|---|---|
| Eq. (1) rewritten as `f_M[1 + (2m_b/(5n_bT_b²))(x̃²−5/2) q·w]` | equal to Homma's bracket, since `w²/(5v_th²) = 2x̃²/5`. Correct |
| `q ≃ κ_∧ b×∇T` | correct as the dominant ⊥ term. We keep all three closure terms |
| `F^∇T` with `sqrt(m_b/(2T_b))` | dimensionally consistent (`C n/T` and `C q sqrt(m/T)/T²` are forces). Re-proved independently in step 1 |
| `stost` → Cartesian `D_∥, D_⊥, A` | checked by Itô calculus for `|v|` and `λ`. Correct, including `−2D_⊥/|v|`. It is used only in `test_stost_cartesian`; the new operator does not use `stost` |
| B07 step as three deterministic increments | equal to the expectation of our random-phase kick (re-derived for the FLR term, §5). We execute the kick, which also gives classical diffusion, and use the increments as test expectations |
| target table (kin6d / benchmark MC) | not used during development. Our own analytic targets come from step 1; this table is compared only at the end |
| tokamak geometry and namelist | adopted as driver 3. Per our decision, the slab driver comes first |
| explicit-friction bias at 0.05 τ_s | not adopted as an input; our own B08 scan decides |
| `boltzmann_mod.f90:94` sets charge 15 e but leaves `cm_over_e` at Z = 74 | confirmed inconsistent: the pusher uses `cm_over_e`, so orbits see Z = 74. Not touched here; the new applet sets all three. To be raised separately |
