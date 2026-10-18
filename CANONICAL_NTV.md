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
/
```

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
- **Restrictions.** Only `grid_kind = 5` is supported. The estimator cannot be
  combined with nonlinear weights, respawn, the profile deposit or anomalous
  kicks.

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

The delta-f torque in the spawn window is `V/N sum P tq_deltaf / T_trace`.

**Full-f decomposition.** Paired runs at `+sigma`, `-sigma` and `0` share the
seed (CRN). For each marker,
`Delta p_phi = dpphi_orbit + dpphi_coll`, and for a lost marker `pphi_end` is
the p_phi it carries out. The NTV, delta-B-loss and axisymmetric-loss parts
are computed outside the code (bench_aug30835 `ntv_circular/ntv_decomposition.py`).
For this, the mesh must carry the same perturbation as the analytic harmonic.
In `grid_kind = 5` the mesh field is unperturbed, so full-f runs give only the
zero-sigma baseline and the drift gate.

**Mesh drift (unperturbed field, no collisions).** On the piecewise-linear
mesh, p_phi is not an exact invariant. The canonical label drifts by
`|Delta s*|` with rms 5.7e-5 at `n2 = 24` over 2e-4 s for 1 keV electrons. The
drift falls roughly as `n2^-2` (1.3e-5 at 48, 4.2e-6 at 96) and does not
depend on the toroidal resolution `n3`.

Tests: `canonical_ntv_kernel` checks the generated torque sample against a
finite-difference dH/dphi at fixed p_phi along an orbit. `canonical_ntv` checks
the zero-sigma weights and the drift, and that the torque is quadratic and even
in sigma with CRN.
