program test_marker_species_reject

! init_species_set must stop with a message for ionisation of the bare nucleus (rate_ion(z_nuc) > 0).
! ctest passes on the message (PASS_REGULAR_EXPRESSION), not on the exit code.

    use, intrinsic :: iso_fortran_env, only: dp => real64
    use marker_species_mod, only: species_set_t, init_species_set

    implicit none

    type(species_set_t) :: ss
    real(dp) :: rate_ion(0:4), rate_rec(0:4)

    rate_ion = 0.0_dp; rate_rec = 0.0_dp
    rate_ion(1) = 1.0e5_dp
    rate_ion(4) = 1.0e5_dp
    call init_species_set(ss, 'Be', [2.01410177812_dp*1.66053906660e-24_dp], [1.0_dp], [1.0e13_dp], [20.0_dp], &
                          20.0_dp, rate_ion, rate_rec, .true.)
    print '(a)', 'test_marker_species_reject: init_species_set accepted ionisation of the bare nucleus'

end program test_marker_species_reject
