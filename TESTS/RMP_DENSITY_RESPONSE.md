# RMP harmonic density diagnostic

The RMP profile file appends `Re(dens_mn) Im(dens_mn)` after its original eight
columns. These are complex marker-density sums demodulated by
`exp(-i(m theta + n phi))`, at the same midpoint and with the same dwell-time,
burn-in, loading and batch normalization as current. The old density columns
remain the axisymmetric marker moment. Current readers selecting the original
columns retain their meaning.

Multiply the density sums by the same volume ratio and inverse batch marker
count used for current, and by the background number density when it is not
already included in the nonlinear marker weight. Number density omits charge
and velocity factors. Preserve both complex components.

This diagnostic is the moment of the evolved marker weight. It does not define
whether that weight represents the full distribution perturbation or its
nonadiabatic part, nor add an adiabatic response. That definition and the
normalization must be fixed before constructing response matrices or solving
quasineutrality. This change supplies no self-consistent potential, physical
transport closure, or validated nonlinear estimator.

The analytic test integrates a sine/cosine density over a full phase period.
It checks the known complex Fourier coefficient, the vanishing mean, constant
density, untouched bins and the zero harmonic. It exercises the helper used
by production deposition, without validating loading or orbit discretization.
