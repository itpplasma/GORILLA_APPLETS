program test_self_consistent_potential
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use constants, only: echarge, ev2erg
    use utils_self_consistent_ef_mod, only: charge_density_to_potential

    implicit none

    real(dp) :: expected

    expected = 2.0_dp*3.5e3_dp*ev2erg/(echarge**2*4.0_dp)
    if (abs(charge_density_to_potential(2.0_dp, 3.5e3_dp, 4.0_dp, &
        .true., 7.0_dp) - expected) > epsilon(expected)) &
        error stop 'static-density potential used the dynamic factor'
    if (abs(charge_density_to_potential(2.0_dp, 3.5e3_dp, 4.0_dp, &
        .false., 7.0_dp) - 7.0_dp*expected) > 8.0_dp*epsilon(expected)) &
        error stop 'dynamic-density potential ignored the configured factor'

    print '(A)', 'PASS: self-consistent potential normalization'
end program test_self_consistent_potential
