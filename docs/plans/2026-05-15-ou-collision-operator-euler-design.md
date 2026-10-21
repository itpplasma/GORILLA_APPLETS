# Ornstein–Uhlenbeck collision operator on v_∥ — Euler–Maruyama form

**Date:** 2026-05-15
**Branch:** `screening-current` (GORILLA_APPLETS)
**Supersedes:** `docs/plans/2026-05-08-ou-collision-operator-design.md`
(exact-OU formulation; replaced by the Euler–Maruyama form below at user
request)

**Scope:** Add an OU collision operator that updates the parallel velocity
v_∥ of a guiding-centre marker via an Euler–Maruyama step, as a cheaper
drop-in alternative to the existing pitch-angle + energy Coulomb operator
in `stost`. Used by the `rmp_response_currents` applet (delta-f screening
current runs); other applets are unaffected.

## 1. Motivation

The full Coulomb collision operator (`stost`, `iswmode = 1`) is an
Euler–Maruyama integrator of the Spitzer–Härm pitch-angle + energy SDE.
For electrons in AUG 33353 it dominates the per-step cost and forces a
small `dtau` to keep `dλ ≪ 1`. For δ-f screening-current studies we
only need a finite, tunable v_∥ decorrelation time τ_c = 1/ν — not a
faithful pitch-angle Coulomb integral. A 1-D OU operator on v_∥ alone
achieves this at much lower per-step cost.

## 2. Operator (Euler–Maruyama, uniform ξ)

Per step of length `dtau`:

```
ξ_s ∼ U[0,1]                          ! uniform draw, one per call
ξ   = sqrt(12) · (ξ_s − 1/2)           ! mean 0, variance 1
Δv_∥ = sqrt(2 · ν · v_T² · dtau) · ξ  −  v_∥ · ν · dtau
v_∥(t + dtau) = v_∥(t) + Δv_∥
```

Choices:

- **Drift toward zero (u = 0):** relaxation to a Maxwellian centred at
  zero parallel velocity.
- **ν = 2 · sum(dhh_vec):** inherited from the existing Coulomb sum so
  the OU damping rate is calibrated to the same v_∥ relaxation rate as
  the full operator. No new free parameter; full local n(s), T(s), |v|
  dependence.
- **v_T² = T_bg / m_bg of the background species.** In stost's
  normalised state (z(4) = |v|/(√2 v_T_test), v_∥_norm = z(4)·z(5))
  this becomes the equilibrium variance
  `sigma_eq² = 0.5 / enrat(i_bg)`, where `enrat(i) = T_test / T_i`.
- **v_⊥ untouched.** Energy is not conserved by this operator.

**Stability constraint.** Euler–Maruyama on OU is stable for
`|1 − ν dtau| < 1`, i.e. `ν dtau < 2`. For acceptable accuracy of the
single-step variance (with a bounded uniform ξ) we want `ν dtau ≲ 0.1`.
The existing `stost` dtau logic (`dhh dtau < q²` with the `z(4) < q`
reduced-step branch) is already tighter than this for the Coulomb
mode, so we add no new step-size cap and keep that logic untouched.

## 3. Placement

Inside `SRC/collis_ions_mod.f90 → stost`, add a new branch after the
`coleff` call (so `dhh_vec` and `enrat` are populated) and **before**
the existing `iswmode == 1 .or. iswmode == 4` block, with an early
`return`:

```fortran
! iswmode = 5: Euler-Maruyama OU step on v_par, leaves v_perp untouched
if (iswmode == 5) then
    ...
    return
end if
```

State-vector handling (z(1:5) = (x, p, λ); 4 = |v|_norm, 5 = pitch):

```fortran
vpar_norm  = z(4) * z(5)
vperp_norm = z(4) * sqrt(max(0.0_dp, 1.0_dp - z(5)**2))
nu_step    = 2.0_dp * sum(dhh_vec)
if (enrat(i_bg) > 0.0_dp) then
    sigma_eq2 = 0.5_dp / enrat(i_bg)
else
    sigma_eq2 = 0.5_dp
end if
if (present(randnum)) then
    ur = randnum(1)
else
    call getran(0, ur)
end if
xi = sqrt(12.0_dp) * (ur - 0.5_dp)
dvpar = sqrt(2.0_dp * nu_step * sigma_eq2 * dtau) * xi &
      - vpar_norm * nu_step * dtau
vpar_norm = vpar_norm + dvpar
z(4) = sqrt(vpar_norm**2 + vperp_norm**2)
z(5) = vpar_norm / z(4)
! standard pmin reflection on z(4) (existing convention)
```

**i_bg choice.** Default `i_bg = num_background_species`, i.e. the
electron entry under the standard `collis_init` ordering for an e/D
plasma. Captured as a module-level constant; no namelist exposure for
now.

**Random number / antithetic.** One uniform draw per step. Antithetic
flipping (ur → 1 − ur) automatically flips ξ → −ξ, so the existing
antithetic-variate plumbing in `carry_out_collisions` works unchanged.

**Edge cases.**

- `nu_step <= 0` → skip update, return ierr = 0.
- `enrat(i_bg) <= 0` → fallback `sigma_eq2 = 0.5`.
- new `z(4) < pmin` → `ierr = ierr + 10`, reflect (existing convention).

## 4. Caller wiring

Add `i_collision_mode` to `&rmp_response_currents_nml` (default `1`):

- `1` → existing Coulomb (`iswmode = 1` in `stost`).
- `5` → Euler–Maruyama OU on v_∥ (this design).

Plumb through `utils_rmp_response_currents_mod.f90` and replace the
hard-coded `iswmode_in = 1` at the `carry_out_collisions` call site with
`in%i_collision_mode`. No other applet reads this flag.

## 5. Files touched

- `SRC/collis_ions_mod.f90` — new `iswmode == 5` branch in `stost`
  (~25 lines).
- `SRC/rmp_response_currents/utils_rmp_response_currents_mod.f90` —
  add `i_collision_mode = 1` at module level, add to the namelist,
  plumb into `in%i_collision_mode`, replace the hard-coded
  `iswmode_in = 1` at the `carry_out_collisions` call site.
- `SRC/gorilla_applets_types_mod.f90` — add
  `integer :: i_collision_mode = 1` to `input_t`.
- `runs/.../rmp_response_currents.inp` — example: add
  `i_collision_mode = 5` for OU runs. Default of 1 keeps legacy runs
  byte-identical.

## 6. Validation plan

1. **Relaxation test (scratch driver).** Single tetra, 10⁵ markers,
   fixed (n, T), `i_collision_mode = 5`. After `t = 5 / ν_step`, expect
     `< vpar_norm² > = sigma_eq² = 0.5 / enrat(i_bg)`  to ~1 %
     `< vpar_norm > = 0`.
   Visual: histogram of vpar_norm against a Gaussian of that variance.
   Deviation in the tails is fine — the single-step ξ is uniform-derived
   so the bounded distribution shows up only on one-step statistics; CLT
   restores Gaussian after a few decorrelation times.
2. **dtau sensitivity.** Same setup, sweep dtau over decades within
   the Euler stability range. Equilibrium variance should hold to within
   MC error until `ν dtau ≳ 0.3`, then drift. Documents the practical
   dtau ceiling and confirms `stost`'s existing dtau logic is conservative
   enough.
3. **End-to-end on AUG #33353 collisional.** Re-run
   `runs/33353_deltaf_colls/` with `i_collision_mode = 5` and compare
   to the Coulomb baseline (`= 1`):
     - Wall-clock: expect close to collisionless.
     - (m,n) = (−6,+2) j_∥ peak amplitude and location: small change
       if `ν_step · t_tot ≲ 1`.
     - Sideband content (−5,+2)/(−7,+2): should not get worse; possibly
       slightly cleaner since the v_∥ relaxation has no aliased
       pitch-angle component.

## 7. Out of scope (YAGNI)

- OU on v_⊥.
- Shifted Maxwellian (u ≠ 0).
- Drag-only / diffusion-only sub-modes (analogues of iswmode 2,3).
- Energy / momentum conservation across antithetic pairs.
- Per-background-species separate ν, T (one combined ν from
  `sum(dhh_vec)`, one `enrat(i_bg)`).

## 8. Rollback

Set `i_collision_mode = 1` (default) → byte-identical to today's
behaviour. The new branch is reached only with explicit opt-in.
