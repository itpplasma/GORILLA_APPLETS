module anomalous_transport_mod

    use, intrinsic :: iso_fortran_env, only: dp => real64

    implicit none

contains

subroutine calc_anomalous_transport

    use orbit_timestep_gorilla_mod, only: initialize_gorilla
    use constants, only: ev2erg
    use tetra_grid_settings_mod, only: grid_kind
    use tetra_grid_mod, only: verts_sthetaphi
    use gorilla_applets_settings_mod, only: i_option
    use field_mod, only: ipert
    use volume_integrals_and_sqrt_g_mod, only: calc_volume_integrals
    use gorilla_applets_types_mod, only: moment_specs, counter, c, in, start, s, output
    use utils_write_data_to_files_mod, only: write_data_to_files, give_file_names, unlink_files
    use utils_data_pre_and_post_processing_mod, only: set_seed_for_random_numbers, &
        get_ipert, set_moment_specifications, initialise_output, initialize_exit_data, calc_poloidal_flux, &
        calc_collision_coefficients_for_all_tetrahedra, normalise_prism_moments_and_prism_moments_squared, &
        fourier_transform_moments
    use utils_self_consistent_ef_mod, only: print_errors_for_bad_inputs, calc_starting_conditions
    use utils_anomalous_transport_mod, only: read_anomalous_transport_inp_into_type, parallelised_particle_pushing_anomalous_transport

    call set_seed_for_random_numbers
    call read_anomalous_transport_inp_into_type
    call get_ipert()
    call initialize_gorilla(i_option,ipert)

    call set_moment_specifications
    call initialise_output
    call calc_volume_integrals(in%boole_boltzmann_energies,in%boole_refined_sqrt_g, in%density, in%energy_eV)
    s%temperature = in%energy_eV

    call initialize_exit_data
    call calc_poloidal_flux(verts_sthetaphi)
    call give_file_names
    call unlink_files
    call print_errors_for_bad_inputs
    call calc_starting_conditions

    if (in%boole_collisions) call calc_collision_coefficients_for_all_tetrahedra
    call parallelised_particle_pushing_anomalous_transport(species=1)
    call normalise_prism_moments_and_prism_moments_squared


    if (moment_specs%n_moments.gt.0) call fourier_transform_moments
    call write_data_to_files

    if (in%boole_precalc_collisions) print*, "maxcol = ", c%maxcol
    print*, 'Number of lost ions',counter%lost_particles
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

end subroutine calc_anomalous_transport

end module anomalous_transport_mod
