program test_deltaf_3d_kernel
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use deltaf_dke_characteristics_kernel, only: evaluate_dke_characteristics
    use deltaf_dke_source_kernel, only: evaluate_dke_source
    use deltaf_dke_step_kernel, only: evaluate_dke_step
    implicit none
    real(dp) :: s, theta, lam, h, u(3), want(3), got, fd, ea, eb, em, w, w_next
    real(dp) :: eps, delta, alpha, minimum, maximum, b, bt, bz, phase, q, pi
    integer :: k, j, mismatches, cases
    ! Dimensionless manufactured chart. Independent cross product of h and
    ! grad B gives radial drift; parallel conservation of mu fixes mirror rate.
    s = 0.4_dp
    theta = 0.7_dp
    lam = -0.3_dp
    h = 1.0e-6_dp
    b = 2.0_dp + s*cos(theta)
    bt = -s*sin(theta)
    call evaluate_dke_characteristics(2.0_dp, 3.0_dp, 5.0_dp, 7.0_dp, lam, &
   b, cos(theta), bt, 0.2_dp, 0.4_dp, 4.0_dp, 6.0_dp, -3.0_dp, 0.0_dp, u(1), u(2), u(3))
    want(1) = (5.0_dp/3.0_dp)*(49.0_dp*(1 + lam**2)/b) &
              *(0.4_dp*0.2_dp - 4.0_dp*bt)/(4.0_dp*6.0_dp + 0.4_dp*3.0_dp)
    want(3) = -7.0_dp*(1 - lam**2)*(0.5_dp*bt + 0.2_dp)/(2*(4 + 0.4_dp*0.5_dp))
    if (abs(u(1) - want(1)) > 1.0e-13_dp .or. &
        abs(u(3) - want(3)) > 1.0e-13_dp) error stop 'Boozer drift/mirror oracle'
    ! A reference F=exp(2s)*(1+0.01 sin(theta)*lambda); the thin-orbit
    ! angular source is independently checked by a directional difference.
    u = [0.0_dp, 0.4_dp, -0.7_dp]
    call evaluate_dke_source(2.0_dp, u(1), u(2), u(3), &
                             0.01_dp*cos(theta)*lam, 0.01_dp*sin(theta), got)
    fd = -(reference(s, theta + h*u(2), lam + h*u(3)) &
           - reference(s, theta - h*u(2), lam - h*u(3)))/(2*h*exp(2*s))
    if (abs(got - fd) > 3.0e-11_dp) error stop 'DKE directional derivative'
    call evaluate_dke_source(2.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 1.0_dp, 1.0_dp, got)
    if (abs(got) > tiny(got)) error stop 'AX source reduction'
  call evaluate_dke_step(0.0_dp, log(4.0_dp), 1.0_dp, 0.0_dp, 3.0_dp, 0.2_dp, em, eb, w)
    if (abs(eb - 4.0_dp) + abs(em - 2.0_dp) + abs(w - 1.2_dp) > 1.0e-14_dp) &
        error stop 'two-weight evolution'
    call evaluate_dke_step(log(4.0_dp), 0.0_dp, eb, w, 0.0_dp, 0.2_dp, em, ea, w_next)
    if (abs(ea - 1.0_dp) > 1.0e-14_dp) error stop 'two-weight reversibility'
    ! Independent large-aspect-ratio circular well oracle: extrema of the
    ! generated parallel mirror force exist iff alpha*=eps|sin theta|/(Nq delta)<1.
    eps = 0.1_dp
    q = 2.0_dp
    pi = acos(-1.0_dp)
    mismatches = 0
    cases = 0
    do j = 1, 4
        delta = 0.0005_dp*real(2**j, dp)
        do k = 1, 37
            theta = pi*real(k - 1, dp)/36
            alpha = eps*abs(sin(theta))/(18*q*delta)
            minimum = huge(1.0_dp)
            maximum = -huge(1.0_dp)
            block
                integer :: z
                do z = 0, 127
                    phase = 2*pi*real(z, dp)/128
                    b = 1 - eps*cos(theta) + delta*cos(phase)
                    bt = eps*sin(theta)
                    bz = -18*delta*sin(phase)
                    call evaluate_dke_characteristics(1.0_dp, 1.0_dp, 0.0_dp, 1.0_dp, &
                      0.0_dp, b, 0.0_dp, bt, bz, 0.0_dp, 1.0_dp, 1.0_dp, -1/q, 0.0_dp, &
                                                      u(1), u(2), u(3))
                    minimum = min(minimum, u(3))
                    maximum = max(maximum, u(3))
                end do
            end block
            cases = cases + 1
     if ((minimum < 0 .and. maximum > 0) .neqv. (alpha < 1)) mismatches = mismatches + 1
        end do
    end do
   print '(a,i0,a,i0)', 'alpha* well oracle: ', cases, ' cases, mismatches ', mismatches
    if (mismatches /= 0) error stop 'ripple well boundary'
    print '(a)', 'PASS: independent drift, source, two-weight and ripple-well oracles'
contains
    function reference(s, t, l) result(f)
        real(dp), intent(in) :: s, t, l
        real(dp) :: f
        f = exp(2*s)*(1 + 0.01_dp*sin(t)*l)
    end function reference
end program test_deltaf_3d_kernel
