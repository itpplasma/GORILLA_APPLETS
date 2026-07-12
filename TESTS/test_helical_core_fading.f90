program test_helical_core_fading
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use helical_core_fading_mod, only: delta_f_fade_factor

    implicit none

    real(dp), parameter :: tolerance = 2.0e-14_dp

    call check('before fade', delta_f_fade_factor(5.0_dp, 8.0_dp, 0.25_dp), &
        1.0_dp)
    call check('fade start', delta_f_fade_factor(6.0_dp, 8.0_dp, 0.25_dp), &
        1.0_dp)
    call check('fade midpoint', delta_f_fade_factor(7.0_dp, 8.0_dp, 0.25_dp), &
        0.5_dp)
    call check('fade end', delta_f_fade_factor(8.0_dp, 8.0_dp, 0.25_dp), &
        0.0_dp)
    call check('after end', delta_f_fade_factor(8.1_dp, 8.0_dp, 0.25_dp), &
        0.0_dp)

    print '(A)', 'PASS: helical-core delta-f fade window'

contains

    subroutine check(label, actual, expected)
        character(len=*), intent(in) :: label
        real(dp), intent(in) :: actual, expected

        if (abs(actual - expected) > tolerance) error stop label
    end subroutine check

end program test_helical_core_fading
