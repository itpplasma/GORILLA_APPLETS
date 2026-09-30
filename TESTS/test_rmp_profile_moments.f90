program test_rmp_profile_moments
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use rmp_profile_moments_mod, only: deposit_profile_moments
    implicit none
    complex(dp) :: moments(6), density, current, expected
    real(dp) :: phase, pi
    integer :: j

    pi = acos(-1.0_dp)
    moments = cmplx(0.0_dp, 0.0_dp, kind=dp)
    ! The Fourier coefficient of 2 cos(alpha) + 3 sin(alpha) is 1-1.5i.
    ! An undemodulated density deposit would instead give zero.
    do j = 0, 31
        phase = 2.0_dp * pi * real(j, dp) / 32.0_dp
        density = cmplx(2.0_dp * cos(phase) + 3.0_dp * sin(phase), 0.0_dp, kind=dp)
        current = 7.0_dp * density
        call deposit_profile_moments(moments, 2, 1, phase, current, density)
    end do
    expected = cmplx(1.0_dp, -1.5_dp, kind=dp)
    if (abs(moments(5) / 32.0_dp - expected) > 1.0e-12_dp) error stop 1
    if (abs(moments(1) / 32.0_dp - 7.0_dp * expected) > 1.0e-12_dp) error stop 2
    if (abs(moments(3)) > 1.0e-12_dp) error stop 3
    if (any(abs(moments([2, 4, 6])) > 1.0e-12_dp)) error stop 4

    ! Constant density has zero helical component and a finite mean.
    moments = cmplx(0.0_dp, 0.0_dp, kind=dp)
    do j = 0, 31
        phase = 2.0_dp * pi * real(j, dp) / 32.0_dp
        call deposit_profile_moments(moments, 2, 2, phase, &
            cmplx(0.0_dp, 0.0_dp, kind=dp), &
            cmplx(4.0_dp, 0.0_dp, kind=dp))
    end do
    if (abs(moments(6)) > 1.0e-12_dp) error stop 5
    if (abs(moments(4) - 128.0_dp) > 1.0e-12_dp) error stop 6
    ! Zero phase is the axisymmetric limit, including complex weights.
    call deposit_profile_moments(moments, 2, 1, 0.0_dp, &
        cmplx(2.0_dp, -3.0_dp, kind=dp), expected)
    if (abs(moments(5) - expected) > 1.0e-12_dp) error stop 7
    if (abs(moments(3) - expected) > 1.0e-12_dp) error stop 8
    print *, 'PASS: density, current, bin isolation and axisymmetric limit'
end program test_rmp_profile_moments
