program test_ou_collision_invariant
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use collis_ions, only: stost, ou_nu_dtau
    implicit none
    real(dp) :: z(5), dt, ef(1), vr(1), er(1), draws(3)
    real(dp) :: perpendicular, vpar, moment, expected_variance
    integer :: ierr, j, seed_size
    integer, allocatable :: seed(:)
    integer, parameter :: samples = 100000

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
    print *, 'PASS: OU perpendicular invariant, zero speed and Gaussian moment', &
        moment, expected_variance
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
