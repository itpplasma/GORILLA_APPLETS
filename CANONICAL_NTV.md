# Canonical delta-f and p_phi torque estimators

The RMP applet (`i_option = 15`) has an optional estimator for the toroidal torque
(NTV). It is off by default. Turn it on with a second namelist group in
`rmp_response_currents.inp`:

```
&canonical_ntv_nml
  boole_canonical_ntv = .true. ,
  ntv_m_pol = 3 , ntv_n_tor = 2 ,   ! harmonic cos(m theta - n phi)
  ntv_amplitude = 1.0 ,             ! sigma (CRN runs: +sigma, -sigma, 0)
  ntv_bmn = 10.0 ,                  ! b_mn [G] when ntv_profile_file is empty
  ntv_phimn = 0.0 ,                 ! Phi_mn [statV]
  ntv_profile_file = '' ,           ! optional table: s, b_mn [G], Phi_mn [statV]
  ntv_theta_kind = 1 ,              ! 1 straight-field-line angle, 0 geometric
  ntv_output_file = 'canonical_ntv_markers.dat' ,
  ntv_chartmap_file = '' ,          ! chartmap whose |B| harmonics give dB1 (grid_kind 6)
  ntv_m_max = -1 ,                  ! chartmap modes m = -ntv_m_max..ntv_m_max; <0: m = ntv_m_pol
  loss_wall_file = '' ,             ! optional (R, Z) [cm] wall polygon as loss boundary
  ntv_marker_rng = .true. ,         ! one random stream per marker (CRN)
  ntv_bmn_relative = .false. ,      ! .true.: b_mn is dB/B, dB1 = sigma b_mn B(x) cos
  ntv_lorentz_nu = 0.0 ,            ! > 0: Lorentz kicks at this rate [1/s] (see below)
  ntv_lorentz_dt = 0.0 ,            ! kick step [s]; <= 0: 1e-4/nu
  ntv_lorentz_mh = .true. ,         ! MH-correct the Lorentz kicks to f0(p_phi, H0)
  ntv_n_check = 0 ,                 ! per-marker checkpoints at k*ntv_check_dt
  ntv_check_dt = 0.0 ,
  ntv_w_reset_dt = 0.0 ,            ! weight control: reset w every ntv_w_reset_dt [s]
  ntv_w_reset_burn = 0.0 ,          ! and tally the torque only after this burn-in [s]
/
```

The toroidal mode number must be a multiple of the mesh field-period count
(`n_field_periods_manual` in `tetra_grid.inp` for `grid_kind = 5`): the mesh
covers one period and markers are relocated into it, so a harmonic that is not
periodic on that interval is wrong (for `n = 3` on the default two-period mesh,
`int sin(3 phi)` over a period does not vanish and every marker picks up a secular
`p_phi` drift). `canonical_ntv_init` stops in that case.

## Model

- **Background.** The background is the canonical Maxwellian `f0(p_phi, H0)` of
  the axisymmetric field (see `ln_f0_and_H`). It uses
  `psi* = psi0 + (c m/e) v_par h_phi` and `H0 = K + e Phi0(s0)`.
- **Loading.** Markers are loaded by the nonlinear loader: uniform physical
  volume and a local Maxwellian. Each marker carries the birth factor
  `P = f0/g0`. Collisions use mode-5 OU with the Metropolis-Hastings
  correction, which leaves `f0` invariant. The temperature is `T_i` for ions
  and `T_e` otherwise.
- **Perturbation.** The perturbation is `dB1 = sigma b_mn(s) cos(m theta - n phi)`
  and the same for `dPhi1`. `B_k` and `A` stay axisymmetric. The perturbation is
  evaluated analytically along the orbits of the mesh field, using the
  generated kernels in `SRC/generated`.
- **Weights.** The weights obey `dw/dt = (d ln f0/dp_phi)|_H0 dH1/dphi` with
  `dH1/dphi = (mu + m v_par^2/B) dB1/dphi + e dPhi1/dphi`. Collision kicks leave
  `w` unchanged. Lorentz scattering is therefore a momentum sink of the test
  species with no restoring term. Its `p_phi` change is tallied separately.
- **Time.** Push times are physical when `boole_time_hamiltonian` is set
  (polynomial pusher). Otherwise the tracer time is used.
- **Charts.** `grid_kind = 5` with `coord_system = 1` (analytic circular field,
  cylindrical coordinates, harmonic in the angle of `ntv_theta_kind`) and
  `grid_kind = 6` with `coord_system = 2` (Boozer chartmap, `x = (s, theta_B,
  phi_B)`, harmonic in the Boozer angles); see "Chartmap route".
- **Restrictions.** The estimator cannot be combined with nonlinear weights,
  respawn, the profile deposit or anomalous kicks.
- **Common random numbers.** With `ntv_marker_rng` (default) every marker
  reseeds its random stream from `seed.inp` and its index, and the MH
  collision step draws its acceptance number on every proposal. The
  `+sigma/-sigma/0` runs then share starts and collision streams marker by
  marker, also when orbits differ (full-f) or earlier markers are lost. The
  number of kicks still depends on the local collision time step, so the
  streams of perturbed full-f orbits slowly lose alignment.

## Monoenergetic Lorentz mode

For monoenergetic transport coefficients (comparison with NEO-RT/NEO-2 at fixed
`E/T`) set `energy_dist_kind = 'mono'` in `rmp_response_currents_nml`: the
canonical loader then gives every marker the energy `energy_eV` and an isotropic
pitch. The temperature and `Phi0` profiles must be flat (checked); then the birth
factor `P = exp(ln F0 - ln M_T(s0))` is `f0` on the energy shell up to a
constant. With `ntv_lorentz_nu > 0` the stost kicks are replaced by energy-
and position-conserving pitch-angle scattering at the constant rate `nu`,
`D_lambda_lambda = nu (1 - lambda^2)/2`. One kick over `dt` proposes
`lambda' = lambda (1 - nu dt) + sqrt(nu (1 - lambda^2) dt) xi` (`xi` normal,
`ntv_lorentz_proposal`). With `ntv_lorentz_mh` it is accepted with
`min(1, pi(lambda') q(lambda|lambda') / (pi(lambda) q(lambda'|lambda)))`, where
`pi = F0` (times the GC phase-space measure with `boole_gc_phase_measure`) per
unit `lambda`. Uniform `lambda` is the invariant measure of the Lorentz
operator at fixed speed, so the kick leaves `F0(p_phi, H0)` exactly invariant.
Proposals with `|lambda'| >= 1` are rejected. The proposal and rejection counts
are printed (`Canonical NTV MH collisions`).

### Collisions and the canonical background

The weights obey the linear equation along the stochastic marker paths,
`dw/dt = S = (d ln f0/dp_phi) dH1/dphi`, and the kicks leave `w` and `P`
unchanged. This is exact if the kick generator `C` leaves `f0` invariant:
- `P` only reweights the initial law, `E_g0[P F(path)] = E_f0[F(path)]` for
  every path functional. The unweighted markers relax from `g0` towards the
  stationary law, but the `P`-weighted ensemble stays `f0` at all times.
- With `C* f0 = 0`, `f0` is stationary for orbits plus kicks. The linear
  response `d(df)/dt = L* df + f0 S` is then
  `int df A = E_f0[w(t) A(z(t))]`, with `w = int_0^t S dt'` along the path
  (time reversal of the stationary process).

The physical Lorentz operator `C_L` conserves a local isotropic distribution,
not `f0`. `f0` depends on `lambda` through
`p_phi = (e/c) psi + m v lambda h_phi`, and
`d ln f0/d lambda = m v h_phi d ln f0/dp_phi ~ rho_pol/L_n`, so `C_L f0 != 0`.
There are two consistent treatments:

1. **MH kicks** (`ntv_lorentz_mh = .true.`, default; also `boole_mh_collisions`
   for the OU kicks). Here `C_MH = C_L - d_lambda[nu (1 - lambda^2)/2 (d_lambda
   ln f0) .]` + O(dt). This is `C_L` plus an `O(rho_pol/L_n)` pitch friction
   that makes `f0` exactly stationary. It matches the thin-orbit theory
   (NEO-RT, NEO-2), where the background is stationary and its axisymmetric
   neoclassical distortion is not part of the NTV drive.
2. **Plain kicks** (`ntv_lorentz_mh = .false.`). The `P`-weighted ensemble then
   evolves as `f0 + df_A(t)`, where `df_A` is the axisymmetric neoclassical
   response to the source `C_L f0`. Its direct torque `int df_A (-dH1/dphi)`
   vanishes because `df_A` does not depend on `phi`. The weight equation still
   drives with `d ln f0/dp_phi` and misses `-{df_A, H1}`, an `O(rho_pol/L_n)`
   relative term.

Both treatments agree to `O(rho_pol/L_n)`. An explicit source weight
`dw_A/dt = C_L f0/f0` together with plain kicks would count `df_A` twice: the
plain kicks already carry it in the marker density. The effect is measured as
a common-random-number difference of the two settings (bench_aug30835
`ntv_circular/outputs/deltaf_rung1_*`). The MH rejection fraction is the size
of the operator change.

## Checkpoints, conversion to D_s, weight control

`ntv_n_check > 0` writes `canonical_ntv_checkpoints.dat`. The header holds the
times `k ntv_check_dt` (physical, Hamiltonian time). Each row holds the marker
index, its physical trace time, then `dh_int`, `tq_deltaf` and `w` at the K
checkpoints. Values are interpolated linearly within the push that crosses
the checkpoint. After the marker's end or loss the final values are kept.

With a constant `kappa = d ln f0/dp_phi` (flat `T`, density exponential in the
poloidal flux), `w = -kappa X` with `X = dh_int`, and the trapezoidal tally gives
exactly `tq_deltaf = -kappa X^2/2` (CTest `canonical_ntv_lorentz`). The delta-f
`p_phi` flux `d<tq>/dt = -kappa D_pp` thus defines the monoenergetic diffusion
coefficient `D_pp(t) = d<X^2>/dt / 2`, the running Green-Kubo integral of
`dp_phi/dt = -dH1/dphi`. With `p_phi = (e/c) psi_pol` (psi per radian) the
flux-surface coefficient is `D_s = D_pp (c/(e dpsi_pol/ds))^2`. This is the
NEO-RT `D11 Dp (ds/dr)^2`, i.e. `Gamma^s = -D_s dn/ds`.

`<w^2> = kappa^2 <X^2>` grows linearly in time once the resonant phase has
decorrelated, and so does the variance of the flux over a fixed window. With
`ntv_w_reset_dt > 0`, `w` is reset to zero at the first push after each
multiple of the reset period. Each cycle is a fresh linear response of the
`f0`-stationary ensemble. `canonical_ntv_reset.dat` holds per marker the torque
integral and the time after the burn-in `ntv_w_reset_burn` of each cycle. The
estimate is unbiased once the burn-in exceeds the time needed for the running
Green-Kubo integral to converge. The variance per cycle is then bounded.

## Output

`canonical_ntv_markers.dat` has one row per marker. The header gives the
spawn volume and the trace time. The columns are:

| column | meaning |
|---|---|
| `n`, `lost` | marker index, 1 if the marker left the mesh before the trace time |
| `P` | birth factor f0/g0 |
| `pphi_birth`, `s_birth`, `sstar_birth` | birth p_phi (mesh A_phi, h_phi), s(psi0), s(psi*) |
| `w_end` | linear delta-f weight at the end of the trace |
| `tq_deltaf` | `int w (-dH1/dphi) dt` |
| `dh_int` | full-f `int (-dH1/dphi) dt` |
| `dpphi_orbit` | p_phi change along the orbit pushes (mesh A_phi, h_phi) |
| `dpphi_coll` | p_phi change by the collision kicks, `m h_phi Delta v_par` |
| `pphi_end`, `s_end`, `sstar_end`, `t_end` | state at the end of the trace or at the loss |
| `n_push`, `dlnf0_dpphi_birth` | number of pushes, `d ln f0/dp_phi` at birth |
| `dh_mesh` | full-f `int (-dH/dphi) dt` with the toroidal derivatives of the mesh `|B|` and `Phi` |
| `psi_birth`, `psi_end` | `psi0` (the `A_phi` of f0) at birth and at the end or loss |
| `loss_kind` | 0 confined, 1 grid exit (last mesh surface or inner boundary), 2 wall, 3 `s_outer_cut` |

The delta-f torque in the spawn window is `V/N sum P tq_deltaf / T_trace`.

**Full-f decomposition.** Paired runs at `+sigma`, `-sigma` and `0` share the
seed (CRN). For each marker,
`Delta p_phi = dpphi_orbit + dpphi_coll`, and for a lost marker `pphi_end` is
the p_phi it carries out. The NTV, delta-B-loss and axisymmetric-loss parts
are computed outside the code (bench_aug30835 `ntv_circular/ntv_decomposition.py`).
For this, the mesh must carry the perturbation: use a perturbed chartmap
(`grid_kind = 6`, "Chartmap route"). In `grid_kind = 5` the mesh field is
unperturbed, so full-f runs give only the zero-sigma baseline and the drift
gate.

**Mesh drift (unperturbed field, no collisions).** On the piecewise-linear
mesh, p_phi is not an exact invariant. The canonical label drifts by
`|Delta s*|`. For 1 keV electrons at `n2 = 24` the rms is 5.7e-5 after 2e-4 s
and 3.8e-4 after 2e-3 s, growing about as t^0.8. The drift falls as
`n2^-1.6` to `n2^-2`: after 2e-4 s it is 1.3e-5 at 48 and 4.2e-6 at 96, after
2e-3 s it is 1.2e-4 and 4.0e-5. It does not depend on the toroidal resolution
`n3`.

Tests: `ntv_lorentz_proposal` checks the MH Lorentz chain against
`pi ~ exp(a lambda)` and the moments of one proposal. `canonical_ntv_lorentz`
checks the mono/Lorentz mode: the identity `tq = -kappa X^2/2`, the checkpoint
end values, the weight-reset bookkeeping and bound, and the field-period guard.
`canonical_ntv_kernel` checks the generated torque sample against a
finite-difference dH/dphi at fixed p_phi along an orbit. `canonical_ntv` checks
the zero-sigma weights and the drift, and that the torque is quadratic and even
in sigma with CRN.

## Chartmap route (grid_kind 6)

The Boozer chartmap (GORILLA `boozer_chartmap_mod`, libneo `neo2bc2chartmap`,
or the circular chartmaps of bench_aug30835 `ntv_circular`) has axisymmetric
`A_theta = torflux s`, `A_phi(s)`, `B_theta(rho)`, `B_phi(rho)` and a 3D `|B|`
on `(rho, theta_B, zeta_B)` nodes (`SRC/rmp_response_currents/chartmap_ntv_mod.f90`).

- **Flux label.** At start-up the applet writes `equil_mapping_file` from the
  chartmap nodes (`psi_pol = A_phi(s)`, `psi_tor = torflux s`,
  `q = -torflux/A_phi'`, mean minor radius `r_eff`), as `grid_kind = 5` writes
  its analytic table. `profile_data_mod` inverts it, `psi* -> s*`. Profiles
  `n(r)`, `T(r)` refer to this `r_eff`. `psi0` is the mesh `A_phi`, which is
  axisymmetric in every run of a set.
- **Loading.** Uniform physical volume in the spawn window by rejection with
  the axisymmetric Jacobian `J0 = (B_phi torflux - B_theta A_phi')/B0^2`,
  `B0` the zeta average of the mesh `|B|` on the chartmap zeta nodes. The
  starts depend only on `J0` and the seed, so the `+sigma/-sigma/0` meshes
  share them.
- **Delta-f.** Mesh: the axisymmetric chartmap. `ntv_chartmap_file`: a
  perturbed chartmap; its `|B|` is projected onto
  `bc cos(m theta - n zeta) + bs sin(m theta - n zeta)` for `n = ntv_n_tor` on
  the rho nodes (exact discrete Fourier projection), splined in `s` and scaled
  by `ntv_amplitude`. `dB1/dphi` comes from the generated harmonic kernel (the
  sine part at a phase shifted by `pi/2`). Without `ntv_chartmap_file` the
  single harmonic `ntv_bmn`/`ntv_profile_file` is used in the Boozer angles.
- **Full-f.** Mesh: the perturbed chartmap (`netcdf_filename`), `ntv_amplitude = 0`.
  The orbits move in the perturbed `|B|`; `dh_mesh` accumulates
  `int -(mu + m v_par^2/B) dB/dphi|mesh dt` (generated `deltaf_weight_kernel`),
  `dpphi_orbit` the p_phi change from the mesh `A_phi`, `h_phi`. All runs of a
  set start from the same physical distribution, the f0 of the unperturbed
  invariants (the delta-f initial condition): the birth factor of the mesh f0
  is corrected by `ln f0(B0) - ln f0(B) + 2 ln(B0/B)` at fixed `v_par` and
  `mu` (zero on an axisymmetric mesh). The torque follows from the CRN second
  difference of `P dh_mesh` (NTV part) and of `P (pphi_end - pphi_birth)`.
- **E_r.** With `ntv_mesh_phi0 = .true.` (default off) the mesh potential on
  the chartmap vertices is set to the profile potential `Phi0(s)` of f0
  (`Er.dat`: `Phi0 = -int E_r dr_eff`, see profile_data_mod), and the mesh
  physics is rebuilt (`flux_mesh_potential_mod`). Orbits (ExB) and
  `f0(p_phi, H0)` then use the same `Phi0`; the toroidal ExB frequency is
  `Omega_E = c dPhi0/dA_phi`. The potential is linear in each tetrahedron, so
  the mesh `dPhi0/ds` is the secant of the radial cell (error
  `|Phi0''| ds_cell/2`); `n1` must resolve the E_r profile. CTest
  `flux_mesh_potential_chartmap`: field-line-label rate equals the cell secant
  to 2e-3 at s = 0.2, 0.4, 0.8 and `Omega_E(s)` the analytic value within the
  element bound; `Phi0 = 0` leaves orbits bitwise unchanged. Off: the mesh
  keeps `eps_Phi A_theta` (`eps_Phi = 0`, empty `er_profile_file`;
  GORILLA's prescribed-E_r loader serves EFIT grids only), while f0 still
  sees `Phi0` from `Er.dat`, so a nonzero `Er.dat` needs the option.
- **Mesh drift.** Sigma = 0, no collisions, 1 keV-scale electrons, 2e-4 s,
  16384 markers: rms `|Delta s*|` 6.1e-5, 1.9e-5, 5.9e-6 for `n2` = 24, 48,
  96 (about `n2^-1.7`), as on `grid_kind = 5` (5.7e-5, 1.3e-5, 4.2e-6). On
  the chartmap mesh `A_phi` depends on `s` only, so it has no spurious phi
  gradient; the drift comes from `h_phi` and `|B|` in the skewed cells.
- **Full-f statistics.** The full-f quadratic part is a small mean of a
  phase-mixed second-order orbit response. At 65536 markers and 2e-5 s its
  standard error is 44 times that of the delta-f torque (-0.19 +- 4.0 against
  2.02 +- 0.09 in the units of the test): consistent, but not a measurement.
  The linear order agrees per marker (corr 0.993 for `dh_mesh`).
- **Not supported on chartmaps.** The refined sqrt(g) prism volumes
  (`boole_refined_sqrt_g`, cylindrical only; moments only).

Tests: `canonical_ntv_chartmap_unit` writes the Boozer chartmap of the analytic
circular field (`write_circular_boozer_chartmap.x`; for this field
`(r, theta*, phi)` are Boozer coordinates) and checks the field, the Jacobian
against the geometry, the Fourier projection, the flux table against GORILLA's
analytic table, the `J0` window volume and the wall polygon.
`canonical_ntv_chartmap` runs the applet: (a) zero weights and the mesh drift
at sigma = 0; (b) chartmap-file harmonic equals the analytic harmonic (same
markers, rel. 4e-14); (c) analytic route and chartmap route agree
statistically; (d) full-f on the perturbed meshes against delta-f: per marker
the odd part of `dh_mesh` equals the linear delta-f response (corr 0.9997 at
2e-6 s, 0.993 at 2e-5 s with collisions; slope 0.94-0.95, because the mesh
resolves the m = 3 harmonic with 24 poloidal points); the quadratic part
agrees within its large statistical error; (e) the wall loss boundary.

## Loss boundary

A marker is lost (`lost = 1`) when its guiding centre leaves the mesh before
the trace time: through the last mesh surface `sfc_s_max` (or the inner
boundary `sfc_s_min`), `loss_kind = 1`; or when it enters a cell beyond
`s_outer_cut`, `loss_kind = 3`. With `loss_wall_file` (default off) a closed
polygon of `R Z` vertices [cm] (a g-file limiter `RLIM, ZLIM` times 100)
is a further loss boundary: a push ending outside it (even-odd rule) is the
loss, `loss_kind = 2`. Coordinates: `(x1, x3)` on cylindrical grids, the
chartmap geometry on `grid_kind = 6`. The wall must lie inside the mesh to
act; outside it the grid exit comes first. `pphi_end`, `psi_end`, `s_end`
and `t_end` record the loss.
