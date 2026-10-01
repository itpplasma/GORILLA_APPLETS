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
