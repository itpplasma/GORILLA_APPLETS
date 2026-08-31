program test_collision_conservation

    use, intrinsic :: iso_fortran_env, only: dp => real64
    use collision_conservation_mod, only: update_background_reservoir

    implicit none

    call test_unequal_masses_and_flow
    call test_zero_ratio

    print *, 'PASS: collision reservoir updates conserve momentum and energy'

contains

    subroutine test_unequal_masses_and_flow

        real(dp), parameter :: test_mass = 2.0_dp
        real(dp), parameter :: background_mass = 5.0_dp
        real(dp), parameter :: ratio = 0.125_dp
        real(dp), parameter :: delta_vpar = -0.7_dp
        real(dp), parameter :: delta_energy = 1.3_dp
        real(dp), parameter :: tolerance = 64.0_dp*epsilon(1.0_dp)
        real(dp) :: background_vpar
        real(dp) :: background_temperature
        real(dp) :: momentum_before, momentum_after
        real(dp) :: energy_before, energy_after

        background_vpar = 0.4_dp
        background_temperature = 3.2_dp
        momentum_before = background_mass*background_vpar
        energy_before = 0.5_dp*background_mass*background_vpar**2 + &
            1.5_dp*background_temperature

        call update_background_reservoir(test_mass, background_mass, ratio, &
            delta_vpar, delta_energy, background_vpar, background_temperature)

        momentum_after = ratio*test_mass*delta_vpar + &
            background_mass*background_vpar
        energy_after = ratio*delta_energy + &
            0.5_dp*background_mass*background_vpar**2 + &
            1.5_dp*background_temperature

        if (abs(momentum_after - momentum_before) > tolerance) &
            error stop 'momentum balance failed'
        if (abs(energy_after - energy_before) > tolerance) &
            error stop 'energy balance failed'

    end subroutine test_unequal_masses_and_flow

    subroutine test_zero_ratio

        real(dp) :: background_vpar
        real(dp) :: background_temperature

        background_vpar = -0.2_dp
        background_temperature = 4.0_dp

        call update_background_reservoir(2.0_dp, 3.0_dp, 0.0_dp, 1.0_dp, &
            2.0_dp, background_vpar, background_temperature)

        if (background_vpar /= -0.2_dp) error stop 'zero ratio changed flow'
        if (background_temperature /= 4.0_dp) &
            error stop 'zero ratio changed temperature'

    end subroutine test_zero_ratio

end program test_collision_conservation
