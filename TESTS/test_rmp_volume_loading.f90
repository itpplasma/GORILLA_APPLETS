program test_rmp_volume_loading
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use rmp_volume_loading_mod, only: volume_candidate, analytic_window_volume
    implicit none
    real(dp), parameter :: r0 = 10.0_dp, a = 3.0_dp, lo = 0.2_dp, hi = 0.8_dp
    real(dp) :: x(3), u(3), probability, norm, mean_r, mean_z, rho2_mean
    real(dp) :: edge, rho2_lo, rho2_hi, expected, pi
    integer :: i, j

    pi = acos(-1.0_dp)
    ! Independent geometry: rho^2 = R0^2 - (R0-s*edge)^2.
    edge = r0 - sqrt(r0*r0-a*a)
    rho2_lo = r0*r0 - (r0-lo*edge)**2
    rho2_hi = r0*r0 - (r0-hi*edge)**2
    norm = 0.0_dp
    mean_r = 0.0_dp
    mean_z = 0.0_dp
    rho2_mean = 0.0_dp
    do i = 1, 32
        do j = 1, 64
            u = [(real(i, dp)-0.5_dp)/32.0_dp, &
                (real(j, dp)-0.5_dp)/64.0_dp, 0.37_dp]
            call volume_candidate(r0, a, lo, hi, u, x, probability)
            norm = norm + probability
            mean_r = mean_r + probability*x(1)
            mean_z = mean_z + probability*x(3)
            rho2_mean = rho2_mean + probability*((x(1)-r0)**2 + x(3)**2)
        end do
    end do
    expected = r0 + (rho2_lo+rho2_hi)/(4.0_dp*r0)
    if (abs(mean_r/norm-expected) > 1.0e-12_dp) error stop 1
    if (abs(mean_z/norm) > 1.0e-12_dp) error stop 2
    if (abs(rho2_mean/norm-(rho2_lo+rho2_hi)/2.0_dp) > 1.0e-12_dp) error stop 3
    expected = 2.0_dp*pi*pi*r0*(rho2_hi-rho2_lo)
    if (abs(analytic_window_volume(r0, a, lo, hi)/expected-1.0_dp) &
        > 1.0e-13_dp) error stop 4
    if (abs(norm/2048.0_dp-r0/(r0+sqrt(rho2_hi))) > 1.0e-13_dp) error stop 5
    print *, 'PASS: physical-volume moments, acceptance and analytic normalization'
end program test_rmp_volume_loading
