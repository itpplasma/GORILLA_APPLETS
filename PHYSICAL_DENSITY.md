# Guiding-center occupation diagnostics

The RMP applet provides two optional settings, both disabled by default:

- `boole_physical_profile_time` uses the Hamiltonian-time occupation moment
  for the harmonic density deposit. Parallel current uses its corresponding
  integrated velocity moment.
- `boole_gc_phase_measure` includes the guiding-center phase-space Jacobian
  relative to a local physical-volume Maxwellian proposal in birth weights
  and in the Metropolis collision acceptance ratio.

The physical-clock lane requires the polynomial pusher, Hamiltonian tracing
and both Hamiltonian-time and parallel-current moments. Forward time and
zero profile burn-in are currently required. The phase-measure lane also
requires nonlinear delta-f weights. With collisions, use zero-flow mode-5
OU proposals, MH correction and a fixed collision background. Unsupported
combinations stop with an explicit error.

The strong-electric-field reference includes the electric-drift kinetic
energy and its contribution to canonical toroidal momentum. Random kinetic
energy remains the Maxwellian proposal variable. A configured equilibrium
Er profile requires the companion GORILLA potential-loading repair in the
weak-field ordering as well as the Hamiltonian-clock repair in the strong
ordering.

The tests independently integrate the phase-space measure and compare the
Metropolis stationary mean and variance with analytic Gaussian moments.
These settings supply guiding-center diagnostics. Finite-Larmor-radius
particle-density pullback, polarization and self-consistent field closure
require separate physical validation.
