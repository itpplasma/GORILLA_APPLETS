# Energy- and momentum-conserving collisions in helical_core

Notes from testing `boole_preserve_energy_and_momentum_during_collisions` in helical_core
(branch `helical_core_collision_test`, November 2026). The operator is described in
section 4.3.2 of J. Schatzlmayr's dissertation (Eqs. 268, 269). It had never been run in
helical_core before.

## How the background update works

Each tetrahedron holds a background reservoir with density `n` (`c%dens_mat`),
parallel flow `u` (`c%vpar_mat`) and temperature `T` (`c%temp_mat`). A marker collides in
the frame moving with `u`. Its changes in parallel velocity `Δv_∥` and energy `Δε`
(lab frame) are then given to the reservoir (`update_background_reservoir`):

    u_new = u_old - r (m_t/m_b) Δv_∥
    T_new = T_old - (2 r Δε + m_b (u_new² - u_old²)) / 3

`r` is the marker-to-background ratio, `w/n_V` in Eq. (269).

## What r has to be

Let a marker stand for `p` real particles. In a collision it gains momentum
`p m_t Δv_∥`. The reservoir of its tetrahedron holds `n V_tet` particles. Momentum and
energy are conserved only if the reservoir takes exactly what the marker loses, that is

    r = p / (n V_tet).

`p` itself is a free choice; choosing it is a relabelling of how many real particles a
marker represents. That only works under two conditions:

1. **`p` must be constant along the marker's trajectory.** A marker cannot represent
   more particles in one tetrahedron than in another. A fixed `r` for every collision
   means `p = r n V_tet`, which changes with the volume of the cell the marker is in.
   Each collision is then consistent with itself, but no fixed `p` exists for which total
   momentum and energy are conserved. On the shipped grid, tetrahedron volumes range from
   2e-5 to 870 cm³ (a factor of about 1e7).
2. **`p` must be the same in the feedback and in the moments.** The moments count a
   marker as `w/n_particles` real particles, because they are divided by `n_particles`.
   If the feedback uses a different `p`, the output describes one plasma and the momentum
   exchange another. Also, if `r` does not depend on `n_particles`, doubling the number of
   markers doubles the momentum fed to the background, so results change with the marker
   number instead of converging.

The code therefore now uses

    r = (w / n_particles) / (n V_tet)

with the exact physical tetrahedron volume `c%tetr_volume`. Their sum matches
`prism_volumes` to 1e-11. The volume is the coordinate volume times R at the centroid,
which is exact because the Jacobian R is linear.

To let the markers represent a different number of real particles, scale `w`, not a
separate coupling factor. Then moments and feedback change together. For the beam test,
`beam_density_fraction` does this.

Changing `p` is not only a noise setting. It sets the ratio of marker particles to
background particles, which is physics:

- **Small `p`:** the background hardly changes, so this is the pure test-particle limit.
- **Large `p`:** strong coupling, but a noisy background in each cell.

In helical_core the markers are full-f samples of the same ions as the background (both
at density `n`). The consistent choice there is `p = w/n_particles`.

### The previous formula

    r = c%weight_factor * w * 1e-4,    c%weight_factor = 1/(n V_box R_max)

With the full-f weight `w = n R V_box` this is `r = 1e-4 R/R_max`, about 1e-4 regardless
of `n_particles`, tetrahedron volume and density. It violates both conditions above. In
addition, helical_core computed `c%weight_factor` before `calc_starting_conditions` had
set the box sizes (`g%amax` etc.). That divided by zero, so the run crashed at once
(SIGFPE).

### Orders of magnitude for the shipped example

Shipped example: n = 3e13 cm⁻³, 1 keV deuterium, 20x30x50 grid with `i_radial_spacing=2`
and `sfc_s_min=1e-7`, 1000 markers.

| quantity | value |
|---|---|
| marker weight `w = n R V_box` | 4e20 to 7e20 |
| real particles per marker `w/N` | about 5e17 |
| check: `Σ w/N` against `n V_plasma` | 3.91e20 = 3.91e20 |
| tetrahedron volume: median / smallest | 54 cm³ / 2e-5 cm³ |
| background particles per tetrahedron: median / smallest | 1.6e15 / 6e8 |
| consistent `r`: median cell / smallest cell | about 340 / about 1e9 |
| previous `r` | about 1e-4 (about 3e6 times too small in a typical cell) |

## Noise: the dilemma of section 4.3.2.2

With a consistent `r`, a single marker can carry many times the particles of a cell. Each
collision then kicks that cell's flow and temperature hard. With the shipped grid (axis
cells 1e-6 times the median) this produced cell flows of about 6 v0 and temperatures of
about 1 eV. The guard that rejects `T <= 0` then fired, which breaks conservation.

The noise can no longer be reduced with a free coupling factor. Only these help:

- **More markers.** `r` per marker scales as 1/`n_particles`.
- **Reservoir cells** (`reservoir_min_volume_fraction`, default 0 = off). Small tetrahedra
  near the axis are merged into cells of at least this fraction of the median tetrahedron
  volume (`build_reservoir_cells`):
  1. **Radial:** starting at the axis, the innermost flux-surface layers are combined
     radially until the complete poloidal ring of every toroidal slice reaches the
     threshold. Their thinness comes from the grid geometry near the axis, not from a
     wish for resolution.
  2. **Poloidal:** the tetrahedra of each layer and toroidal slice are merged in poloidal
     order into cells of at least the threshold volume. This starts at the axis and stops
     at the first layer without small tetrahedra.

  There is no toroidal merging. All tetrahedra of a cell share flow and temperature, and
  the deferred update sums the exchange over the cell. A layer is identified by the lowest
  vertex flux of a tetrahedron, a slice by its lowest vertex φ.

  On the small-hole grid (`sfc_s_min = 1e-7`, `i_radial_spacing = 2`, 10x10x20) with
  fraction 0.3, the 8 innermost of 48 layers form one radial layer, 8170 of 19200
  tetrahedra (all at s < 0.012) are merged, and every cell reaches the threshold.
- **Large axis hole: not recommended.** An earlier workaround used `sfc_s_min = 0.05` with
  `i_radial_spacing = 1`. The axis-hole reflection (x → 3·axis − 2x) then moves markers by
  several cm, which breaks canonical toroidal angular momentum: a 9.5 % leak of the beam's
  angular momentum over 500 µs. With a small hole (`sfc_s_min = 1e-7`) the reflection
  hardly ever happens (0 or 1 events per run), and the angular momentum is conserved to
  below 1e-3.
- **Not implemented:** a reservoir per complete flux surface. Parallel flows even out along
  field lines much faster than collisions act, so this would be physically reasonable and
  would cut the noise further.

## Parallelisation and common time (`n_time_slices`)

Markers are traced in parallel with OpenMP over markers, with a static schedule and
per-thread seeds (dissertation section 4.3.2.3). Each thread traces its markers one after
another, each over the full `time_step`. At any moment only one marker per thread is
being traced, and these markers are at unrelated physical times. A thread's 1000th marker
starts at t = 0 against a background that has already absorbed the full histories of all
earlier markers.

For a stationary end state this acts like an iteration, which is the argument of section
4.3.2.3. For a relaxation followed in time it is wrong. Single-slice test: at t = 500 µs
the beam markers end at 0.099 v0 and the background at 0.154 v0. Markers slower than the
background are impossible, because markers traced late relaxed toward a background that
had already been spun up.

With `n_time_slices > 1`, all markers are pushed by `time_step/n_time_slices`. Then the
next slice starts from where each marker stopped (`update_start_type`). The
parallelisation over markers is unchanged, and markers and background stay in step to
within one slice. Within a run, the weights are computed once and kept fixed. A
zero-length step sets each marker's weight before its first collision.

**Random numbers.** The per-thread seeding runs at the start of every parallel region. With
several slices, every marker therefore got the same random numbers in every slice. That
noise adds up coherently and heated a beam relaxing in a fixed 1000 eV background to about
1130 eV. Each slice now uses its own seed factor (the slice index). A single slice keeps
factor 1, so results are unchanged.

**Deferred background updates.** In helical_core the background is frozen during a slice.
Each thread adds the momentum and energy its markers gain in each tetrahedron into its own
slot of `background_exchange` (`start_deferred_background_updates`). After the slice,
`apply_deferred_background_updates` sums the slots in thread order and updates every
reservoir once with the totals.

- Runs are reproducible for a fixed number of threads; two runs with 16 threads gave
  bit-identical output.
- There is no critical section per collision.
- There is no race between reading and updating the background.
- The same structure works for MPI: one `MPI_Allreduce` of the exchange per slice, then
  the same update on every rank.

Updates that would drive a temperature to T <= 0 are rejected and reported per slice as a
warning, because conservation is broken there. They signal too much exchange in one cell
during one slice, so the remedy is shorter slices or larger reservoir cells. In a test with
1e4 markers, deferred and immediate updates agreed within the noise. Boltzmann and SCEF
keep the immediate update after every collision.

**Performance note.** Every marker adds an ntetr-sized moment array to the output inside
a critical section. With slices this happens once per marker and slice and dominates the
run time: 1 slice took 40 s, 20 slices took 7 min.

## Other observations

- Only ion background species exchange momentum and energy (`j = 1 .. c%n-1`). Electrons
  are left out in conserving mode, while the non-conserving operator includes them.
- With immediate updates (Boltzmann, SCEF), the background flow and temperature are read
  outside the critical section and updated inside it. Other threads can change them in
  between, as noted in section 4.3.2.3. This does not affect conservation, because the
  update uses the current values. With the deferred updates of helical_core the background
  is frozen during a slice, so this race does not occur.
- The time step does not set the noise level. A smaller collision step makes each kick
  smaller but more frequent. The drag adds up to the same total and the random parts add up
  to the same variance, so the noise in a cell's flow scales like `p/(n V_cell)` regardless
  of the step. Limiting the step in tiny cells would only remove single large jumps.
- In a torus, Σ m v_∥ is not conserved by the orbits (the direction of b changes, and the
  mirror force acts once μ > 0). The conserved quantity is the canonical toroidal angular
  momentum, m v_∥ h_φ + (q/c) A_φ for the markers plus m_b u h_φ for the background.
  Orbits alone keep it to 1e-4 over 500 µs.

## Beam relaxation test

Input keys (helical_core.inp):

| key | default | meaning |
|---|---|---|
| `n_time_slices` | 1 | push all markers in common slices |
| `boole_unidirectional_beam` | false | all markers start with pitch +1; needs `boole_antithetic_variate = .false.` |
| `beam_density_fraction` | 1 | marker density / background density |
| `reservoir_min_volume_fraction` | 0 | merge small tetrahedra near the axis into reservoir cells of at least this fraction of the median volume |

Output with conserving collisions:

- `collision_relaxation.dat`, one row per slice: mean marker v_∥, mean background flow,
  momentum, energy and canonical angular momentum of markers, background and lost
  markers, and temperatures.
- `background_reservoir.dat`: the reservoir state per tetrahedron at the end.
- `markers_final.dat`: the marker states at the end.

Expected result for a cold beam (density α n, energy E0) in a background at rest at
temperature T0, same mass:

- common flow u_f = α/(1+α) v0
- common temperature T_f = (α E0 + 1.5 T0 − (1+α) E0 (u_f/v0)²) / (1.5 (1+α))

For α = 0.2 and E0 = T0 = 1 keV: u_f = 0.167 v0 and T_f = 926 eV.

This estimate assumes conservation of Σ m v_∥, which does not hold in a torus. In the
runs about 11 % of Σ m v_∥ goes to the field, and markers and background meet at about
0.14–0.15 v0. The quantity that is conserved is the canonical toroidal angular momentum:
to below 1e-3 with a small axis hole. Total energy is conserved to 1e-4 or better.
Results (beam inside s < 0.7, 20 slices):

| run | markers | meet at [v0] (markers / bg) | T_m / T_b [eV] | ΔE/E | ΔL/L_beam | min T_b [eV] |
|---|---|---|---|---|---|---|
| hole 0.05, equidistant s | 6e4 | 0.139 / 0.143 | 897 / 926 | -1e-5 | 0.095 | 330 |
| hole 1e-7, no cells | 3e4 | 0.150 / 0.145 | 880 / 912 | 3e-3 | -0.001 | 0.1 |
| hole 1e-7, poloidal cells | 3e4 | 0.139 / 0.147 | 889 / 909 | 2e-4 | 0.000 | 8 |

The two small-hole runs still had 20 rejected updates (T <= 0) in the innermost rings.
The radial merging removes those cells (run r26).

Checks of the parallel-velocity handling:

- The final v_∥ distributions are Maxwellians shifted by the mean flow.
- In a fixed background (α = 1e-4) the beam relaxes towards rest and 1000 eV.
- A standalone `stost` ensemble without orbits gives the same relaxation.

Run setup and plots: `helical/runs/make_run.py` and `helical/runs/plot_relaxation.py` in
the working area (not in the repository).
