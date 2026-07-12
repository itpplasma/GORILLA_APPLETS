program test_helical_core_profile
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use helical_core_profile_mod, only: linear_density_factor, &
        delta_f_density_weight

    implicit none

    real(dp), parameter :: tolerance = 2.0e-13_dp
    real(dp), parameter :: step = 1.0e-3_dp
    real(dp) :: derivative, source, s, s_zero, velocity

    s = 0.4_dp
    s_zero = 1.25_dp
    velocity = -0.03_dp
    derivative = (linear_density_factor(s + step, s_zero) &
        - linear_density_factor(s - step, s_zero))/(2.0_dp*step)
    source = delta_f_density_weight(1.0_dp, velocity, s_zero)

    call check('configured profile', linear_density_factor(s, s_zero), 0.68_dp)
    call check('profile derivative', derivative, -1.0_dp/s_zero)
    call check('delta-f source', source, velocity/s_zero)
    call check('legacy default', linear_density_factor(0.4_dp, 1.1_dp), &
        0.7_dp/1.1_dp)
    if (delta_f_density_weight(0.7_dp, velocity, 1.1_dp) &
        /= 0.7_dp*velocity/1.1_dp) error stop 'legacy operation order'

    print '(A)', 'PASS: helical-core density profile and delta-f source agree'

contains

    subroutine check(label, actual, expected)
        character(len=*), intent(in) :: label
        real(dp), intent(in) :: actual, expected

        if (abs(actual - expected) > tolerance*max(1.0_dp, abs(expected))) &
            error stop label
    end subroutine check

end program test_helical_core_profile
