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

## Coupling: physical content or free capacities

A collision of a marker that stands for `p = w/n_particles` real particles changes its cell by

    u_new = u_old - r_V (m_t/m_b) Δv_∥,     r_V = p / N_V
    T_new = T_old - (2 r_T Δε + m_b (u_new² - u_old²)) / 3,     r_T = p / N_T

`N_V` and `N_T` are the capacities of the cell for flow and temperature, in particles.
The cell then conserves `N_V m_b u` and `N_T (3/2 T + m_b u²/2)` together with the markers.
This is Eq. (28) of Runov, Kasilov & Helander, J. Comput. Phys. 300 (2015) 605, with
`N_V = W_V`, `N_T = W_T` and the marker weight `w`.

What the capacities must be depends on the problem.

- **Relaxation in time** (e.g. the beam test): the background is a second, physical
  population, and how fast and where markers and background meet depends on how many
  particles it holds. Then `N_V = N_T = n V_cell`, the physical content of the cell. This
  is the default (`reservoir_capacity_flow = reservoir_capacity_temperature = 0`).
- **Stationary state** (e.g. neoclassical flows and currents): the background only stands
  for the field particles in the collision operator. In a stationary state its flow and
  temperature no longer change, so the summed kicks in every cell vanish:
  Σ w Δv_∥ = 0 and Σ w Δε = 0 per cell. These are the conservation laws of the collision
  operator (Runov et al., Eqs. 25–27), and the capacities drop out of them. `N_V` and
  `N_T` are then free numerical parameters, and they may even differ between cells. They
  set only how noisy the cells are and how fast they relax:
  - variance of a cell's flow ≈ v_T² p/N_V, of its temperature ≈ T² p/N_T (Runov et al.,
    App. C);
  - the markers must spend a total time of about τ_c N/p in a cell before it has relaxed.

  Runov et al. use `N_T/p` = 50–100 and, for the parallel flow in a tokamak core (small
  compared with v_T), `N_V/p` = 1e5.

In the stationary case the capacities are set in helical_core.inp in units of the mean
marker weight, `reservoir_capacity_flow` = W_V/w and `reservoir_capacity_temperature` = W_T/w
(both > 0, or both 0 for the physical content). The mean marker weight is the physical
particle number of the first background species in the mesh divided by `n_particles`
(`set_reservoir_capacities`), i.e. the mean weight if the markers sample that species. Every
cell gets the same capacities, independent of its volume, for every background species.
A collision of a marker of weight w_n then changes its cell by r = (w_n/w_mean)/(W/w). Small
cells are visited less often, so they relax more slowly, but they are not noisier.

**Correction (2026-11-18).** Earlier versions of this file stated two conditions that
supposedly forced the physical content in all cases: the number of particles a marker
represents had to be the same in every cell, and the same in the feedback as in the moments,
otherwise results would change with the marker number. Both treat the background as a
physical population whose momentum has to balance the markers' at every moment. That is
true for a relaxation in time, but not for a stationary state, where only the vanishing
summed exchange per cell matters. The statement "the noise can no longer be reduced with a
free coupling factor" was wrong for stationary problems.

### The previous formula

    r = c%weight_factor * w * 1e-4,    c%weight_factor = 1/(n V_box R_max)

With the full-f weight `w = n R V_box` this is `r = 1e-4 R/R_max`, about 1e-4 regardless
of `n_particles`, tetrahedron volume and density. For a stationary state this was a free
capacity of about 1e4 marker weights (the same for flow and temperature), which is allowed; for
the beam relaxation it was wrong. In
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
| `r` with physical content: median cell / smallest cell | about 340 / about 1e9 |
| previous `r` | about 1e-4 (capacity about 3e6 times the physical content of a typical cell) |

## Noise: the dilemma of section 4.3.2.2

With the physical content as capacity, a single marker can carry many times the particles of
a cell. Each collision then kicks that cell's flow and temperature hard. With the shipped grid
(axis cells 1e-6 times the median) this produced cell flows of about 6 v0 and temperatures of
about 1 eV. The guard that rejects `T <= 0` then fired, which breaks conservation.

For a stationary state the remedy is larger capacities (`reservoir_capacity_flow`,
`reservoir_capacity_temperature`, see above). For a relaxation in time only these help:

- **More markers.** `r` per marker scales as 1/`n_particles`.
- **Reservoir cells** (`reservoir_min_volume_fraction`, default 0.5, 0 = off). Small tetrahedra
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

Markers are traced in parallel with OpenMP over markers, with a static schedule
(dissertation section 4.3.2.3). Since November 2026 each marker has its own random stream
(`initialise_seed_for_marker`, from seed.inp, the marker index and the slice).
Before, each thread had its own seed. Each thread traces its markers one after
another, each over the full `time_step`. At any moment only one marker per thread is
being traced, and these markers are at unrelated physical times. A thread's 1000th marker
starts at t = 0 against a background that has already absorbed the full histories of all
earlier markers.

For a stationary end state this acts like an iteration, which is the argument of section
4.3.2.3 and of the single-particle algorithm of Runov et al. For a relaxation followed in time it
is wrong. Single-slice test: at t = 500 µs
the beam markers end at 0.099 v0 and the background at 0.154 v0. Markers slower than the
background are impossible, because markers traced late relaxed toward a background that
had already been spun up.

With `n_time_slices > 1`, all markers are pushed by `time_step/n_time_slices`. Then the
next slice starts from where each marker stopped (`update_start_type`). The
parallelisation over markers is unchanged, and markers and background stay in step to
within one slice. Within a run, the weights are computed once and kept fixed. A
zero-length step sets each marker's weight before its first collision.

**Stationary states.** Common time is not needed there, but with the deferred updates below
the background only changes between slices. With `n_time_slices = 1` it is updated once,
after all markers have been traced, so the markers never feel the exchange (helical_core
prints a warning). For a stationary state, use many slices: each slice is then one step of
the iteration between markers and background. Every slice prints the largest change of a
cell (`|du_par|/v0`, `|dT|/T`); the per-collision updates of Runov et al. are reproduced
while these stay small. Large capacities make them small.

**Random numbers.** The seed must differ between slices. The per-thread seeding first ran at
the start of every parallel region with the same factor, so with several slices every
marker got the same random numbers in every slice. That noise adds up coherently and heated
a beam relaxing in a fixed 1000 eV background to about 1130 eV. Each slice therefore uses
its own seed factor (the slice index). The per-marker streams keep this: the stream of a
marker depends on the marker index and the slice.

**Deferred background updates.** In helical_core the background is frozen during a slice.
Each thread adds the momentum and energy its markers gain in each tetrahedron into its own
slot of `background_exchange` (`start_deferred_background_updates`). After the slice,
`apply_deferred_background_updates` sums the slots in thread order and updates every
reservoir once with the totals.

- Runs are reproducible. The marker histories do not depend on the number of threads.
  The background exchange is summed in thread order, so with several slices its rounding
  can differ between thread counts. Two runs with 16 threads gave bit-identical output.
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

- Electrons are included since November 2026; see "Electrons" below. Before that, the
  conserving loop ran over the ion species only (`j = 1 .. c%n-1`), so markers felt no
  electron drag or scattering at all.
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

## Electrons (instantaneous e-i friction)

The marker collides with every background species, electrons included (species `c%n`).
The electrons get a reservoir temperature `T_e`. They do not get an independent flow:
e-i friction (microseconds) is much faster than the slowing-down, so the electron flow is
locked to the mass-weighted mean ion flow. `update_electron_reservoir` handles the
exchange with the electrons:

- **Momentum:** all background flows (ions and electrons) shift by the same amount, so that
  the electron flow ends at the mean ion flow and the total momentum, including the
  electron flow, is conserved.
- **Energy:** the exchanged energy minus the change of the flow kinetic energy goes into `T_e`.
- **Rejection:** an update that would make `T_e <= 0` is rejected, as for the ions.

With deferred updates, this runs once per reservoir cell and slice, after the ion
updates. The momentum the markers lose to electrons therefore ends up in the ion
rotation.

**Not modelled (possible later switch, "option A"):** a finite e-i friction time and
T_e-T_i equipartition. Option A would give the electrons their own flow `u_e` and, in
`apply_deferred_background_updates`, exchange momentum between `u_e` and `u_i` (time τ_ei)
and energy between `T_e` and `T_i` (equipartition time). Both exchanges can be made exactly
conservative. Option A matters when T_e ≠ T_i must relax, or when the slice length is
not long compared with τ_ei. An independent `u_e` would be very noisy: Δu = Δp/(n m_e) is
m_i/m_e larger than for ions. Option A therefore needs the friction step, and a free `u_e`
alone is not enough.

**Common collision step.** `stost` chooses the step `dtau` itself when it is given `tau`,
and overwrites the caller's `t%step`. Calling it once per species therefore let the last
species (electrons, weak pitch-angle scattering, long step) set the time step for all of
them. Ion collisions were then applied for a shorter time than the marker advanced; in a
test, T_m reached only 703 eV instead of 900 eV. The step is now chosen once from the
pitch-angle diffusion summed over all species (`collision_step`, the same rule as in
`stost`), and each species is applied with that step.

(The precalculated random numbers, `boole_precalc_collisions`, which would have given all
species the same three numbers, have been removed in favour of per-marker streams.)

**Checks** (binary of branch `conserving_electrons`, small hole, cells threshold 0.5, 24 threads):

| run | result |
|---|---|
| 1 keV beam, α = 0.2, 6e4 markers, 500 µs (as the threshold-0.5 run above) | meet 0.143 / 0.147 v0, T_m 900 eV, T_b 926 eV, T_e 1000 → 996 eV; same as without electrons within noise (expected: E_c ≈ 19 keV, so electrons take about 1 % of the drag); E_tot conserved to 1e-10, 0 rejections |
| 30 keV beam on a 1 keV background (`background_temperature_eV = 1000`), 3e4 markers, 4 slices of 50 µs | electrons take 64 % of the beam's energy loss in the first slice (estimate 1/(1+(E_c/E)^1.5) = 0.665 with E_c = 18.7 keV); loss rate 324/s against 2(1+(E_c/E)^1.5)/τ_s ≈ 376/s (lnΛ approximated); E_tot conserved to 1e-10, 0 rejections |

`background_temperature_eV` (helical_core.inp, default 0 = use `energy_eV`) sets the
background temperature of all species independently of the marker energy. It needs the
conserving collisions, because the other operator precomputes its coefficients from
`energy_eV`.

## Beam relaxation test

Input keys (helical_core.inp):

| key | default | meaning |
|---|---|---|
| `n_time_slices` | 1 | push all markers in common slices |
| `boole_unidirectional_beam` | false | all markers start with pitch +1; needs `boole_antithetic_variate = .false.` |
| `beam_density_fraction` | 1 | marker density / background density |
| `reservoir_min_volume_fraction` | 0.5 | merge small tetrahedra near the axis into reservoir cells of at least this fraction of the median volume |
| `background_temperature_eV` | 0 | background temperature of all species if > 0, else `energy_eV` (conserving collisions only) |
| `reservoir_capacity_flow` | 0 | capacity W_V/w of every reservoir cell for the flow, in mean marker weights; 0 = physical content (stationary states only, see "Coupling") |
| `reservoir_capacity_temperature` | 0 | same for the temperature (W_T/w); set both or neither |

Output with conserving collisions:

- `collision_relaxation.dat`, one row per slice: mean marker v_∥, mean background flow,
  momentum, energy and canonical angular momentum of markers, background and lost
  markers, and temperatures. P_b, E_b and L_b include all background species. The last two
  columns are the mean and minimum electron temperature.
- `background_reservoir.dat`: the reservoir state per tetrahedron at the end (last column T_e).
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
| hole 1e-7, radial + poloidal cells | 3e4 | 0.139 / 0.146 | 883 / 909 | 1.4e-4 | 0.0006 | 1.4 |
| hole 1e-7, cells, threshold 1.0 | 3e4 | 0.136 / 0.147 | 885 / 913 | 1e-10 | 0.0002 | 129 |
| hole 1e-7, cells, threshold 0.3 | 6e4 | 0.146 / 0.146 | 892 / 923 | 0 | 0.00004 | 9 |
| hole 1e-7, cells, threshold 0.5 | 6e4 | 0.140 / 0.148 | 895 / 922 | 1e-10 | 0.00008 | 184 |
| hole 1e-7, cells, threshold 0.5 | 1.2e5 | 0.148 / 0.147 | 901 / 929 | -2e-10 | 0.0001 | 322 |

The first three small-hole runs had rejected updates (T <= 0): 20, 20 and 16. With radial plus
poloidal cells, every cell reaches the threshold of 0.3 times the median volume. The
coldest cells are then merged cells near s ≈ 0.001–0.002 with volumes just above the
threshold (T_b = 17 eV, u_b = −0.5 v0 at the end). A threshold of 0.3 times the median
volume is therefore not enough for 3e4 markers near the axis. Larger cells or more markers
remove the rejections:

- **Threshold 1.0:** safe even with 3e4 markers, but merges about half of the grid (the
  five innermost rings), which costs resolution of the rotation near the axis.
- **Threshold 0.3:** no rejections with 6e4 markers, but the coldest cell still reaches 9 eV.
- **Threshold 0.5 (the default):** merges one ring more than 0.3, needs about 6e4 markers
  (2.3e4 inside the beam region), and keeps every cell above 180 eV.

With 1.2e5 markers the cell-to-cell scatter of u_b drops by √2 (0.148 to 0.105 v0). It is
still about 70 % of the mean flow, so rotation profiles should be read as flux-surface
averages.

Checks of the parallel-velocity handling:

- The final v_∥ distributions are Maxwellians shifted by the mean flow.
- In a fixed background (α = 1e-4) the beam relaxes towards rest and 1000 eV.
- A standalone `stost` ensemble without orbits gives the same relaxation.

Run setup and plots: `helical/runs/make_run.py` and `helical/runs/plot_relaxation.py` in
the working area (not in the repository).
