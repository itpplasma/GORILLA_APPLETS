# OU collision invariant

Collision mode 5 advances parallel velocity and retains perpendicular velocity.
Its old post-step speed reflection at p_min=0.05 changes both components and
therefore is incompatible with that contract and with the Gaussian OU transition.

Remove the reflection only from mode 5. Retain the other collision modes' speed
floor. At exactly zero speed, set pitch to zero rather than dividing by zero.

The test calls the production stost routine with an externally specified rate.
It checks a deterministic low-speed perpendicular invariant, a zero-speed and
zero-step limit, and the independently known Gaussian second moment of 100000
exact one-step OU draws at fixed perpendicular velocity. The statistical bound
is seven standard errors of the Gaussian second-moment estimator.

Without an external rate, mode 5 used nu = 4 dhh(|v|) at the marker's own
speed, and the step dtau = ou_nu_dtau/nu (or epsilon^2/(2 dhh) for the
Euler-Maruyama branch) therefore grew with speed. Each step keeps the Gaussian
v_par law, but the time spent per step does not: a time-sampled Maxwellian is
weighted by 1/nu(v) and heats up (<E/T> = 2.17 instead of 1.5 after twelve
collision times in the applet stationarity test). Mode 5 now evaluates the
rate and the step at the background thermal speed sqrt(2 sigma_eq2), so both
depend on position only; an external rate (nu_override) was already speed
independent and now also sets the Euler-Maruyama step. A deterministic check
requires equal steps for markers at p = 0.2 and p = 3 in both branches; the
previous source fails it.

The legacy getran stub errors if called: the tested paths explicitly supply
increments or use intrinsic Gaussian sampling. Constants are the core's actual
constants module. No orbit tracing or field closure is simulated here.

On aCluster, job 21804866 (2 CPUs, 2 GiB, 5 minutes) first built both original
and repaired kernels. The original production source failed the deterministic
invariant oracle; the repaired kernel passed all checks. A local fo pipeline
also passed. The first isolated job could not compile without the core constants
module; it produced no physics evidence and was repaired before comparison.

This fixes a kernel defect. It does not establish an unbiased nonlinear marker
measure, validate the Metropolis proposal or prove layer-response convergence.
