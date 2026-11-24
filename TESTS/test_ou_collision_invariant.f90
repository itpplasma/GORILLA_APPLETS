program test_ou_collision_invariant
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use collis_ions, only: onseff, stost, ou_nu_dtau
    implicit none
    real(dp) :: z(5), dt, ef(1), vr(1), er(1), draws(3)
    real(dp) :: perpendicular, vpar, moment, expected_variance, dt_slow, dt_fast
    integer :: ierr, j, seed_size
    integer, allocatable :: seed(:)
    integer, parameter :: samples = 100000

    call test_collision_coefficients()

    ef = 0.0_dp
    vr = 1.0_dp
    er = 1.0_dp
    draws = 0.0_dp
    ou_nu_dtau = 0.0_dp
    ! With v_parallel=0 and zero supplied increment, perpendicular speed is invariant.
    z = [0.0_dp, 0.0_dp, 0.0_dp, 0.01_dp, 0.0_dp]
    dt = 0.01_dp
    call stost(ef, vr, er, z, dt, 5, ierr, randnum=draws, nu_override=1.0_dp)
    perpendicular = z(4) * sqrt(max(0.0_dp, 1.0_dp-z(5)**2))
    if (abs(perpendicular-0.01_dp) > 1.0e-14_dp) error stop 1
    if (ierr /= 0) error stop 2

    ! A valid zero-speed limit must not divide by zero or reflect to a finite speed.
    z = 0.0_dp
    dt = 0.0_dp
    call stost(ef, vr, er, z, dt, 5, ierr, randnum=draws, nu_override=1.0_dp)
    if (z(4) /= 0.0_dp .or. z(5) /= 0.0_dp) error stop 3

    ! Exact OU one-step second moment: sigma^2 (1-exp(-2 nu dt)).
    ! Low variance makes the old p_min reflection conspicuous, not a rare-tail test.
    ou_nu_dtau = 1.0_dp
    er = 5000.0_dp
    call random_seed(size=seed_size)
    allocate(seed(seed_size))
    seed = 271828
    call random_seed(put=seed)
    moment = 0.0_dp
    do j = 1, samples
        z = [0.0_dp, 0.0_dp, 0.0_dp, 0.01_dp, 0.0_dp]
        dt = 1.0_dp
        call stost(ef, vr, er, z, dt, 5, ierr, nu_override=1.0_dp)
        perpendicular = z(4) * sqrt(max(0.0_dp, 1.0_dp-z(5)**2))
        if (abs(perpendicular-0.01_dp) > 1.0e-12_dp) error stop 4
        vpar = z(4) * z(5)
        moment = moment + vpar**2
    end do
    moment = moment / samples
    expected_variance = 0.0001_dp * (1.0_dp-exp(-2.0_dp))
    ! Seven standard errors of an independent Gaussian second-moment estimator.
    if (abs(moment/expected_variance-1.0_dp) > &
        7.0_dp*sqrt(2.0_dp/real(samples, dp))) error stop 5

    ! Without an external rate the OU rate and step must not depend on the
    ! marker speed: a v_par-dependent step length biases the time-sampled
    ! distribution towards fast markers (a Maxwellian heats up).
    ef = 1.0e-3_dp
    vr = 1.0_dp
    er = 1.0_dp
    do j = 1, 2
        if (j == 1) ou_nu_dtau = 0.1_dp
        if (j == 2) ou_nu_dtau = 0.0_dp
        z = [0.0_dp, 0.0_dp, 0.0_dp, 0.2_dp, 0.5_dp]
        call stost(ef, vr, er, z, dt_slow, 5, ierr, tau=1.0e30_dp, randnum=draws)
        z = [0.0_dp, 0.0_dp, 0.0_dp, 3.0_dp, 0.5_dp]
        call stost(ef, vr, er, z, dt_fast, 5, ierr, tau=1.0e30_dp, randnum=draws)
        if (abs(dt_fast/dt_slow - 1.0_dp) > 1.0e-12_dp) error stop 6
    end do

    print *, 'PASS: OU perpendicular invariant, zero speed and Gaussian moment', &
        moment, expected_variance

contains

    subroutine test_collision_coefficients()
        real(dp), parameter :: speeds(*) = [0.0_dp, 1.0e-5_dp, 0.00999_dp, &
                                           0.01_dp, 0.01001_dp, 0.1_dp, &
                                           1.0_dp, 6.0_dp, 10.0_dp]
        integer, parameter :: intervals = 4096
        real(dp) :: actual(3), expected(3), v, x, t2, weight
        integer :: i, k

        do i = 1, size(speeds)
            v = speeds(i)
            x = v**2
            expected = 0.0_dp
            ! Integrate the Maxwellian and its speed derivative without subtracting
            ! nearly equal error-function and Gaussian terms.
            do k = 0, intervals
                t2 = (real(k, dp) / intervals)**2
                weight = real(2 + 2 * mod(k, 2), dp)
                if (k == 0 .or. k == intervals) weight = 1.0_dp
                expected = expected + weight * exp(-x*t2) * &
                           [2*t2, 1-t2, 4*t2*(1-x*t2)]
            end do
            expected = expected * (2 / sqrt(acos(-1.0_dp))) / (3 * intervals)
            call onseff(v, actual(1), actual(2), actual(3))
            if (.not. all(abs(actual/expected - 1) <= 1.0e-11_dp)) then
                print *, 'Collision coefficients at normalized speed', v
                print *, 'Relative errors (parallel, perpendicular, derivative)', &
                    abs(actual/expected - 1)
                error stop 'Collision coefficients disagree with Maxwellian integrals'
            end if
        end do
    end subroutine test_collision_coefficients
end program test_ou_collision_invariant

subroutine getran(mode, value)
    use, intrinsic :: iso_fortran_env, only: dp => real64
    implicit none
    integer, intent(in) :: mode
    real(dp), intent(out) :: value
    ! Both tested OU branches supply their own random numbers.
    value = real(mode, dp)
    error stop 'unexpected legacy random source'
end subroutine getran
