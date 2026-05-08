# Ornstein–Uhlenbeck collision operator on v_∥

**Date:** 2026-05-08
**Branch:** `screening-current` (GORILLA_APPLETS)
**Scope:** Add an OU collision operator that updates the parallel velocity v_∥
of a guiding-centre marker via an exact Ornstein–Uhlenbeck step, as a cheaper
drop-in alternative to the existing pitch-angle + energy Coulomb operator in
`stost`. Used by the `rmp_response_currents` applet (delta-f screening-current
runs); other applets are unaffected.

## 1. Motivation

The Coulomb collision operator (`stost`, `iswmode=1`) is an Euler–Maruyama
integrator of the linearised Spitzer–Härm pitch-angle + energy SDE. For
electrons in AUG 33353 it dominates the per-step cost and forces small
`dtau` to keep `dλ << 1`. For delta-f screening-current studies we mostly
need a finite, tunable v_∥ decorrelation time τ_c = 1/ν — not a faithful
Coulomb collision integral. An exact OU step on v_∥ achieves this at
essentially zero stability constraint on `ν dtau`.

## 2. Operator

For each step of length `dtau`:

```
v_∥(t + dtau) = u + (v_∥(t) - u) · exp(-ν · dtau)
              + sqrt( (T/m) · (1 - exp(-2 ν · dtau)) ) · ξ
```

with ξ ~ N(0,1). Choices for this design:

- **u = 0** (relaxation to Maxwellian centred at zero).
- **ν = 2 · sum(dhh_vec)**: inherited from the existing Coulomb collision
  rate so the OU damping rate is calibrated to match the v_∥ relaxation
  rate of the existing operator. No new free parameter, full local
  n(s), T(s), v dependence.
- **Equilibrium variance** in the units `stost` uses
  (z(4) = |v|/(√2 v_T_test), v_∥_norm = z(4)·z(5)):
    `sigma_eq² = 0.5 / enrat(i_bg)`
  where `enrat(i) = e0/T_i`. The test species relaxes to the Maxwellian of
  background `i_bg` (default: electron entry, the last in `collis_init`'s
  ordering for typical setups).
- **v_⊥ untouched.** Energy is **not** conserved by this operator.

The exact OU step is unconditionally stable in `ν dtau`. The only
remaining `dtau` constraint is keeping ν up-to-date with the local
plasma parameters; we use `dtau = min(τ, upper_limit, c_ou / ν_step)`
with `c_ou = 1.0`. The `z(4) < q` reduced-step branch from the Coulomb
operator is dropped — it existed to control `dhh ∝ 1/v` blow-up in the
Euler scheme, which the exact OU step does not have.

## 3. Placement

Inside `SRC/collis_ions_mod.f90:stost`, add a new branch:

```fortran
! iswmode = 5: exact OU step on v_∥, leaves v_⊥ untouched
if (iswmode == 5) then
    ! coleff has already been called above; use dhh_vec, enrat
    ...
    return
end if
```

placed after the `coleff` call and the `dpp/dhh/fpeff` sums, before the
`iswmode == 1 .or. iswmode == 4` branch. Early `return` so it does not
fall through to the existing pitch / energy branches.

State-vector handling (z(1:5) = (x, p, λ)):

1. `vpar_norm = z(4) * z(5)`
2. `vperp_norm = z(4) * sqrt(max(0, 1 - z(5)**2))`
3. Apply the OU step above to `vpar_norm`.
4. `z(4) = sqrt(vpar_norm**2 + vperp_norm**2)`
5. `z(5) = vpar_norm / z(4)`
6. Apply the standard `pmin` reflection on `z(4)`.

Random number: one ξ per step. Use `randnum(1)` if `present(randnum)`,
otherwise `call getran(0, ur)` (same convention as the energy branch).
This keeps the existing antithetic-variate plumbing in
`carry_out_collisions` working unchanged.

Edge cases:

- `enrat(i_bg) <= 0` → fallback `sigma_eq² = 0.5d0`.
- `ν_step <= 0` → skip the update, return with `ierr = 0`.
- New `z(4) < pmin` → `ierr = ierr + 10`, reflect (existing convention).

## 4. Caller wiring

Add `i_collision_mode` to the `rmp_response_currents.inp` namelist
(default `1`):

- `1` → existing Coulomb (`iswmode = 1` in `stost`).
- `5` → OU on v_∥.

Read it in
`SRC/rmp_response_currents/utils_rmp_response_currents_mod.f90` along
with the other namelist entries, store in `in%i_collision_mode`, and
pass it through at the existing `carry_out_collisions` call:

```fortran
call carry_out_collisions(i, n, t, x, vpar, vperp, ind_tetr, iface, &
                          species, iswmode_in = in%i_collision_mode)
```

(currently hard-coded to `iswmode_in = 1`).

No other applet reads this flag, so no other `i_option` is affected.

## 5. Files touched

- `SRC/collis_ions_mod.f90`
  - `stost`: add `iswmode == 5` branch (~30 lines), guarded by an early
    `return`. May factor out a private `ou_step_vpar(...)` helper if the
    branch grows past ~40 lines.
- `SRC/rmp_response_currents/utils_rmp_response_currents_mod.f90`
  - Add `i_collision_mode` to the module-level integer block, default 1.
  - Add to namelist `&rmp_response_currents_nml`.
  - Plumb into `in%i_collision_mode`.
  - Replace the hard-coded `iswmode_in = 1` at the `carry_out_collisions`
    call site with `in%i_collision_mode`.
- `SRC/rmp_response_currents/rmp_response_currents_mod.f90`
  - No code change expected; the namelist is already read via
    `utils_rmp_response_currents_mod`.
- `runs/.../rmp_response_currents.inp` — example: add
  `i_collision_mode = 5` for OU runs (no change needed for legacy runs;
  the default of 1 reproduces current behaviour).

## 6. Validation plan

1. **Unit-style relaxation test (manual).** A scratch driver with a
   single tetra (or just z(4), z(5) state, no orbit), 10⁵ particles, fixed
   (n, T), `i_collision_mode = 5`. After `t = 5 / ν_step`, confirm
   `< (z(4) · z(5))² > = 0.5 / enrat(i_bg)` to ~1%. Visual: histogram of
   v_∥_norm against `N(0, sqrt(sigma_eq²))`.
2. **dtau-independence.** Same setup, sweep `dtau` over three decades.
   The equilibrium variance must be flat — this is the point of the
   exact step. If it drifts, the implementation is wrong.
3. **End-to-end on AUG 33353.** Re-run the electron delta-f screening
   case `runs/33353_deltaf` with `i_collision_mode = 5`; compare to the
   existing collisionless run (`boole_collisions = .false.`) and to the
   Coulomb run that was too slow. Expectations:
     - Throughput close to collisionless.
     - n=2 j_∥ peak still present if `ν_step · t_tot ≲ 1`.
     - Over-damped (peak vanishes) if `ν_step · t_tot ≫ 1` —
       a useful sensitivity diagnostic in itself.

## 7. Out of scope (YAGNI)

- OU on v_⊥.
- Shifted Maxwellian (u ≠ 0). Easy to add later: one extra namelist
  entry for `u_par` (or a profile lookup). Not needed for the
  screening-current physics goal.
- Drag-only / diffusion-only OU sub-modes (analogous to iswmode=2,3 of
  the Coulomb operator).
- Energy or momentum conservation across antithetic pairs.
- Generalising to multiple background species with separate ν, T (the
  design uses one combined ν from `sum(dhh_vec)` and one `enrat(i_bg)`).

## 8. Rollback

If the OU mode misbehaves, set `i_collision_mode = 1` (default) in the
input file. The Coulomb path is byte-identical to today's behaviour.
