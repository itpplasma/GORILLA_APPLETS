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
`pi = F0` per unit `lambda`. Uniform `lambda` is the invariant measure of the Lorentz
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
Two kick models keep the weight equation consistent:

1. **MH kicks** (`ntv_lorentz_mh = .true.`, default; also `boole_mh_collisions`
   for the OU kicks). Here `C_MH = C_L - d_lambda[nu (1 - lambda^2)/2 (d_lambda
   ln f0) .]` + O(dt). This is `C_L` plus an `O(rho_pol/L_n)` pitch friction
   that makes `f0` exactly stationary.
2. **Plain kicks** (`ntv_lorentz_mh = .false.`). The `P`-weighted ensemble then
   evolves as `f0 + df_A(t)`, where `df_A` is the axisymmetric neoclassical
   response to the source `C_L f0`. Its direct torque `int df_A (-dH1/dphi)`
   vanishes because `df_A` does not depend on `phi`, and the weight equation
   still drives with `d ln f0/dp_phi` only.

In both the NTV drive is `-{f0, H1}` alone. This is not the complete
collisional drive. The stationary axisymmetric solution of the physical
operator is `f_M + f_AX` with `f_AX = G + h`: `G = f0 - f_M` is the canonical
shift, `h` the collisional remainder, which is non-zero in the passing region
and in the trapped-passing boundary layer. The term `-{h, H1}` is of the same
order in `rho*` as the drive `-{G, H1}`. NEO-2-QL includes it in production (`nobounceaver`,
Kernbichler et al. PPCF 58 (2016) 104001 Eq. (28), Kasilov et al. PoP 21
(2014) 092506 Eq. (65)). The term vanishes as `nu -> 0`: there `h` is a function
of the invariants in the passing region and zero for trapped particles. Plain
kicks do not supply it: the weights are still driven by `d ln f0/dp_phi`. An
explicit source weight `dw_A/dt = C_L f0/f0` together with plain kicks would
count `df_A` twice, because the plain kicks already carry it in the marker
density. The MH rejection fraction is the size of the operator change.

### Collisional background drive (`ntv_fax_file`, opt-in)

`ntv_fax_file` adds the mirror-force part of `-{h, H1}` to the weight equation
(default off; with it unset the code path and results are those of the
estimator without it). In `(x, v, lambda)` the kinetic energy is conserved
and `h` is axisymmetric, so only the perturbed mirror force acts on it:

    d lambda/dt = -(v (1 - lambda^2)/2) b.grad ln B,
    dw/dt += -delta(d lambda/dt) d(h/f0)/d lambda
           = (d ln f0/dp_phi) mu (h_phi/R)^2 (d chi/d lambda) dB1/dphi,

with `h/f0 = (d ln f0/dp_phi) m v h_phi chi`, so `chi = lambda h/G` and
`chi = lambda` for `h = G`. For an m = 0 harmonic `b.grad = b^phi d/dphi` and
`b^phi = h_phi/R^2` on grid_kind 5; other `m` remain rejected there. The perturbed radial drift
acting on `h` is `O(rho*)` smaller and is dropped, as in NEO-2. The kernels
`deltaf_fax_coordinates_kernel` and `deltaf_fax_weight_kernel` are generated
(`gen/deltaf_kernels`, identities F1-F3). F3 checks the normalisation and sign
against the existing source: for `chi = lambda` (`d chi/d lambda = 1`) the
kernel returns `(d ln f0/dp_phi) mu (h_phi/R)^2 dB1/dphi`. This is the `mu`
part of the canonical source, i.e. the mirror force on `G`. The remaining
part, `m v_par^2/B dB1/dphi`, is the perturbed radial drift on `f_M`.

The table holds `d chi/d lambda` on a grid of
- `vartheta = branch arccos(1 - 2 (B - Bmin)/(Bmax - Bmin))`, where `branch`
  is the sign of `dB` along `b`, i.e. along the motion of `lambda > 0`, and
- the signed pitch `lambda = v_par/v` in `[-1, 1]`.

At fixed position `chi` is one function of the signed `lambda`, continuous
through the turning point `lambda = 0`, so `d chi/d lambda` stays finite there
(`d chi/d eta` grows as `1/lambda`). These labels do not depend on the
poloidal angle or orientation convention of the code that solved for `h`.
On grid_kind 5, `Bmin` and `Bmax` are those of the circular surface through
the marker (outboard and inboard midplane of the analytic field). On
`grid_kind=6, coord_system=2`, extrema are sampled from the mesh chartmap
through GORILLA's field evaluator at all unique angular nodes on 2001
surfaces spanning its rho range, indexed by `s=rho^2` and interpolated in
rho as in the core. Angular extrema of its multilinear interpolant occur
at nodes; changes of the extremizing node between tabulated surfaces leave
a radial tabulation error. `branch` uses the local orbit mesh's `b.grad B`,
including both angular derivatives, rather than a cylindrical height.

The chartmap source supports every retained poloidal mode (cosine and sine):

    dw/dt += (d ln f0/dp_phi) mu h_phi (d chi/d lambda)
             * (b^theta dB1/dtheta + b^phi dB1/dphi),
    b^i = curl(A)^i / (h_j curl(A)^j).

Here the components and gradients are those of the orbit tetrahedron.
For fixed covariant Boozer `B_k` and `A_k`, `b^i` scales with `B`, so the
variation of `b.grad ln B` is `b^i d_i B1/B`: the variation of `b^i` cancels
the denominator derivative. A single helical mode contains `m*iota-n`;
retaining just `d/dphi` loses the `m*iota` term. The generated
`deltaf_fax_chartmap_kernel` evaluates the contraction and parallel gradient.
With relative amplitudes, the derivative includes the local `B` factor.
The cylindrical m=0 path and its generated kernel are unchanged.

This remains local thin-orbit use of a single-surface table. The labels
cannot distinguish different magnetic wells with the same B and branch;
a multiwell background needs a richer table. Device/species tables and
transport agreement with NEO-2 require separate qualification.
Format: `#` header lines,
`n_theta n_lambda`, the `vartheta` grid (increasing, spanning `[-pi, pi]`),
the `lambda` grid (increasing, spanning `[-1, 1]`), then `n_theta` rows of
`n_lambda` values. Lookup is bilinear and clamped to the table.

CTest `canonical_ntv_fax_chartmap_unit` compares the circular chartmap's
extrema, position labels and m=0 source with cylindrical analytic values
(source tolerance 0.5%), and m=3 with the analytic helical mirror derivative
(tolerance 1%). The fine radial fixture resolves A_phi' near m*iota=n;
resonant cancellation, field reversal and relative amplitudes are also checked.
`canonical_ntv_fax_chartmap` checks a zero table, unchanged orbits and the
m=3 spectral versus explicit source with eight markers. Its nonzero table
is artificial and supplies no transport reference.

The table is built outside the code from the axisymmetric NEO-2-QL solution
`g_k0` (drive k = 1) as `h = g_k0 - G`. In NEO-2 the `sigma = +1` solution
moves towards `+phi_mfl`, so `branch = sign(d B/d phi_mfl)`.

NEO-2 also has the term `-(B_n/B_0) L_C g_k0`, which is not included. Our
reading, not yet checked separately, is that it comes from the Boozer `|B|`
perturbation changing `b.grad phi`, which GORILLA's perturbation model does
not have. Removing it changes the NEO-2 result by about 1 % (`coll0` patch).

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
are computed outside the code.
For this, the mesh must carry the perturbation: use a perturbed chartmap
(`grid_kind = 6`, "Chartmap route"). In `grid_kind = 5` the mesh field is
unperturbed, so full-f runs give only the zero-sigma baseline and the drift
gate.

**Mesh drift (unperturbed field, no collisions).** On the piecewise-linear
mesh, p_phi is not an exact invariant. The canonical label `s*` drifts slowly
(about as `t^0.8`); the drift falls as `n2^-1.6` to `n2^-2` with the poloidal
resolution and does not depend on the toroidal resolution `n3`.

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

### Energy and species resolved background tables

The original numeric `n_theta n_lambda` format remains the monoenergetic block-0
`d chi/d lambda` route. Its reader, interpolation and weight arithmetic are
unchanged. For a physical multispecies NEO-2 background, use `FAX_ENERGY_V1`:

```text
# optional leading comments
FAX_ENERGY_V1
n_theta n_lambda n_basis n_species degree n_segments extend_tail
mass_g(n_species)
charge_number(n_species)
temperature_erg(n_species)
source_amplitude(3,n_species)
basis_edges(n_segments+1)
basis_centers(n_segments)
basis_coeff(0:degree,n_basis,n_segments)
vartheta(n_theta)
lambda(n_lambda)
g(n_lambda,n_theta,n_basis,3,n_species,n_species)
dg_dlambda(same dimensions as g)
```

All arrays use list-directed input and Fortran order (leftmost index fastest).
`extend_tail` is a Fortran logical. Position spans `[-pi,pi]` and pitch
`[-1,1]`. The last two axes are **driving species, field species**. Species are
matched to markers by mass and charge; absent or ambiguous species fail.
Energies are kinetic energies in erg; each marker uses
`x=sqrt(E/T_field_species)`, not the driving species' temperature. Polynomial
basis segments are supplied by NEO-2, retaining its normalized Laguerre basis
or actual B-spline knots, normalization and extrapolation. A zero-tail basis
returns zero outside its final edge; continuation uses the final polynomial.

The values are NEO-2 response functions normalized to the field-species
Maxwellian. `source_amplitude` holds `beta_b A_k^b`, in internal NEO-2 order
(density, parallel electric field, temperature). Each source and driving
species is summed independently. `ntv_fax_response` returns `f_AX/f_M` or its
pitch derivative at fixed energy and position. The marker weight uses
`d(h/f_M)/d lambda = d(f_AX/f_M)/d lambda - kappa m v h_phi`, with
`kappa=d ln f0/dp_phi`, so a vanishing profile gradient never divides the
multispecies response. Horner evaluation, speed conversion and canonical-shift
subtraction come from the FortSym generator.

Generate these tables with NEO-2's `lsw_export_fax` feature and
`python/src/neo2_ql/fax.py`. Its raw HDF5 set retains the distribution band
integrals and source-force normalization. The converter preserves native `avbhat` factors in the source amplitudes and
requires every field species and an explicit electric field (or the completed native output).
A single table covers one radial surface and one single-well position chart;
using it across different profiles/surfaces requires separate qualification.
The table does not assert physical transport convergence.

`canonical_ntv_fax_energy` tests two field species, both driving species,
all sources, energy dependence and the polynomial tail against independent
analytic values. For the native NEO-2 distribution oracle, configure
`FAX_NEO2_READER` with NEO-2's `fax.py` path and `FAX_NEO2_EXPORTS` with its
semicolon-separated export files. CTest then builds a real table and compares
the APPLETS reader with native NEO-2 band moments (2% matrix-scale tolerance,
including position/pitch interpolation).

For admission of an existing device table, set `FAX_NATIVE_TABLE` and
`FAX_NATIVE_TABLE_EXPORTS` in addition to `FAX_NEO2_READER`. Native CTest
checks its checksum and evaluates every field species against the native band
moments without running marker orbits or a physical pilot.

`GORILLA_REF` records the compatible core commit used by CI and the native
checks. The energy/species feature uses the same frozen core as #80.

### Local collision mixture profiles

With `boole_local_background`, the legacy input without a species file describes
one singly charged ion and electrons, both at `n.dat`. For a mixture, provide
`profile_dir/collision_species.dat`: the species count on the first line, then
mass [g], charge number, and a density-profile filename relative to `profile_dir`
for each background, ions first and electrons last. For example:

```
3
3.345228e-24  1  nD.dat
5.017842e-24  1  nT.dat
9.109100e-28 -1  n.dat
```

Each density file has the same two-column radius [cm], density [cm^-3] format
and coordinate mapping as `n.dat`. Set `n_background_species=3` for this example.
The source supplies each density independently: D/T fractions are never inferred
from electron density. All ions use `Ti.dat`, electrons use `Te.dat`. The loader
and every centroid evaluation check `sum(Z*n)=0` to relative tolerance `1e-10`
(normalized by `sum(abs(Z*n))`), nonnegative finite densities, and agreement of
the supplied electron profile with `n.dat` to relative tolerance `1e-10`. Invalid
inputs stop the run rather than being renormalized. The file is optional for the
legacy pure plasma, whose collision coefficients retain their previous bits.

This changes the backgrounds passed to `collis_init`, retaining its native CGS
constants and pair Coulomb logarithms. The canonical Metropolis-Hastings path
still requires mode-5 OU collisions; these profile checks do not establish
full-Fokker-Planck transport parity.
