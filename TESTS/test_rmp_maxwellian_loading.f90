program test_rmp_maxwellian_loading
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use rmp_maxwellian_loading_mod, only: maxwellian_from_uniforms
    implicit none
    integer, parameter :: ns = 200000
    integer, allocatable :: seed(:)
    integer :: i, nseed, ntail
    real(dp) :: u(4), energy, pitch, sum_e, sum_pitch, sum_pitch2, expected, pi

    call random_seed(size=nseed)
    allocate(seed(nseed))
    seed = 61003
    call random_seed(put=seed)
    sum_e = 0.0_dp
    sum_pitch = 0.0_dp
    sum_pitch2 = 0.0_dp
    ntail = 0
    do i = 1, ns
        call random_number(u)
        call maxwellian_from_uniforms(u, energy, pitch)
        sum_e = sum_e + energy
        sum_pitch = sum_pitch + pitch
        sum_pitch2 = sum_pitch2 + pitch*pitch
        if (energy > 5.0_dp) ntail = ntail + 1
    end do
    ! Independent Maxwell-Boltzmann energy Gamma(3/2,1) and isotropic pitch.
    pi = acos(-1.0_dp)
    expected = erfc(sqrt(5.0_dp))+2.0_dp*sqrt(5.0_dp/pi)*exp(-5.0_dp)
    if (abs(sum_e/ns-1.5_dp) > 7.0_dp*sqrt(1.5_dp/ns)) error stop 1
    if (abs(sum_pitch/ns) > 7.0_dp/sqrt(3.0_dp*ns)) error stop 2
    if (abs(sum_pitch2/ns-1.0_dp/3.0_dp) &
        > 7.0_dp*sqrt((1.0_dp/5.0_dp-1.0_dp/9.0_dp)/ns)) error stop 3
    if (abs(real(ntail, dp)/ns-expected) &
        > 7.0_dp*sqrt(expected*(1.0_dp-expected)/ns)) error stop 4
    call maxwellian_from_uniforms([0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp], energy, pitch)
    if (energy /= 0.0_dp .or. pitch /= 0.0_dp) error stop 5
    print *, 'PASS: Maxwellian energy, angular moments and untruncated tail'
end program test_rmp_maxwellian_loading
