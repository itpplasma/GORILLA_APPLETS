module rmp_response_currents_mod

    use, intrinsic :: iso_fortran_env, only: dp => real64

    implicit none

contains

subroutine calc_rmp_response_currents

    use orbit_timestep_gorilla_mod, only: initialize_gorilla
    use constants, only: ev2erg
    use tetra_grid_settings_mod, only: grid_kind
    use tetra_grid_mod, only: verts_rphiz
    use gorilla_applets_settings_mod, only: i_option
    use field_mod, only: ipert
    use volume_integrals_and_sqrt_g_mod, only: calc_square_root_g, calc_volume_integrals
    use gorilla_applets_types_mod, only: moment_specs, counter, c, in, start, output, weights
    use utils_write_data_to_files_mod, only: write_data_to_files, give_file_names, unlink_files
    use utils_data_pre_and_post_processing_mod, only: set_seed_for_random_numbers, &
        get_ipert, set_moment_specifications, initialise_output, initialize_exit_data, calc_poloidal_flux, &
        calc_collision_coefficients_for_all_tetrahedra, normalise_prism_moments_and_prism_moments_squared, &
        fourier_transform_moments, calc_starting_conditions
    use utils_helical_core_mod, only: eliminate_particles_outside_flux_threshold
    use utils_rmp_response_currents_mod, only: read_rmp_response_currents_inp_into_type, &
        parallelised_particle_pushing_rmp_response_currents, &
        profile_dir, equil_mapping_file, boole_constant_delta_B_psi, delta_B_psi_const, &
        pert_m_fourier, pert_n_fourier, delta_B_psi_file, boole_skip_phase_for_test, &
        bias_starting_positions_to_s_window
    use profile_data_mod, only: load_profiles
    use perturbation_field_mod, only: init_constant_perturbation, load_perturbation_field, boole_skip_phase

    call set_seed_for_random_numbers
    call read_rmp_response_currents_inp_into_type
    call get_ipert()
    call initialize_gorilla(i_option,ipert)

    call set_moment_specifications
    call initialise_output
    call calc_square_root_g
    call calc_volume_integrals

    call initialize_exit_data
    call calc_poloidal_flux(verts_rphiz)
    call give_file_names
    call unlink_files

    if (in%boole_collisions) call calc_collision_coefficients_for_all_tetrahedra

    ! Load plasma profiles and initialise the perturbation field for the
    ! delta-f weight evolution (Albert 2016, Eq. 4).
    if (in%boole_delta_f) then
        call load_profiles(trim(profile_dir), trim(equil_mapping_file))
        if (boole_constant_delta_B_psi) then
            call init_constant_perturbation(delta_B_psi_const, pert_m_fourier, pert_n_fourier)
            boole_skip_phase = boole_skip_phase_for_test
        else
            call load_perturbation_field(trim(delta_B_psi_file), trim(equil_mapping_file), &
                                         pert_m_fourier, pert_n_fourier)
        end if
    end if

    call calc_starting_conditions
    call eliminate_particles_outside_flux_threshold
    if (in%boole_delta_f) call bias_starting_positions_to_s_window

    call parallelised_particle_pushing_rmp_response_currents(species=1)

    call normalise_prism_moments_and_prism_moments_squared(boole_skip_time_normalisation_in=in%boole_delta_f)

    if (moment_specs%n_moments.gt.0) call fourier_transform_moments
    call write_data_to_files

    if (in%boole_precalc_collisions) print*, "maxcol = ", c%maxcol
    print*, 'Number of lost ions', counter%lost_particles
    print*, 'average number of pushings = ', counter%tetr_pushings/in%n_particles
    print*, 'average number of toroidal revolutions = ', counter%phi_0_mappings/in%n_particles
    print*, 'average number of integration steps = ', counter%integration_steps/in%n_particles
    PRINT*, 'ion mass = ', start%particle_mass(1)
    PRINT*, 'absolute value of ion velocity = ', start%v0(1)
    PRINT*, 'ion charge = ', start%particle_charge(1)
    PRINT*, 'temperature = ', ev2erg*in%energy_eV
    print*, 'energy in eV = ', in%energy_eV
    print*, 'tracing time in seconds = ', in%time_step
    if((grid_kind.eq.2).or.(grid_kind.eq.3)) then
         print*, 'number of times that particles were pushed across the inside hole = ', counter%lost_inside
    endif

end subroutine calc_rmp_response_currents

end module rmp_response_currents_mod
