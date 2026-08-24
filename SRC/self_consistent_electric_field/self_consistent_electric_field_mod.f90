module self_consistent_electric_field_mod

    use, intrinsic :: iso_fortran_env, only: dp => real64

    implicit none
   
contains

subroutine calc_self_consistent_electric_field

    use orbit_timestep_gorilla_mod, only: initialize_gorilla
    use constants, only: ev2erg,pi,echarge
    use omp_lib, only: omp_get_num_threads, omp_get_thread_num
    use tetra_grid_settings_mod, only: grid_kind
    use tetra_grid_mod, only: ntetr, verts_sthetaphi
    use gorilla_settings_mod, only: ispecies
    use gorilla_applets_settings_mod, only: i_option
    use field_mod, only: ipert
    use volume_integrals_and_sqrt_g_mod, only: calc_volume_integrals_in_flux_coordinates
    use gorilla_applets_types_mod, only: moment_specs, counter, c, in, start, s, output, weights, ep
    use utils_write_data_to_files_mod, only: write_data_to_files, give_file_names, unlink_files
    use utils_data_pre_and_post_processing_mod, only: set_seed_for_random_numbers, &
    get_ipert, set_moment_specifications, initialise_output, initialize_exit_data, calc_poloidal_flux, &
    fourier_transform_moments
    use utils_scef_input_mod, only: read_self_consistent_electric_field_inp_into_type, print_errors_for_bad_inputs
    use utils_scef_particle_init_mod, only: calc_starting_conditions, calc_desired_density_profile, &
    initialize_particle_source, update_particle_source, write_particle_source
    use utils_scef_electric_potential_mod, only: allocate_electric_potential_type, initialize_flux_shell_grid
    use utils_scef_diffusion_mod, only: calc_diffusion_coefficients

    integer :: species, source_step

    call set_seed_for_random_numbers
    call read_self_consistent_electric_field_inp_into_type
    call get_ipert()
    call initialize_gorilla(i_option,ipert)
    
    call set_moment_specifications
    call initialise_output
    call calc_volume_integrals_in_flux_coordinates
    s%temperature = in%energy_eV

    !if (.not.in%boole_static_ne) call calc_electron_diffusion_coefficients

    call initialize_exit_data
    call initialize_flux_shell_grid
    call calc_desired_density_profile
    call initialize_particle_source
    call calc_poloidal_flux(verts_sthetaphi)
    call allocate_electric_potential_type
    call give_file_names
    call unlink_files
    call print_errors_for_bad_inputs

    !Prep diffusion coefficients once per species that will use the random walk.
    !Skipping this step reuses the cached A_and_B_species<species>.dat if present.
    do species = 1,2
        if ((.not.in%boole_honest_tracing(species)).and.in%boole_recompute_D(species)) then
            call calc_diffusion_coefficients(species)
        endif
    enddo

    do source_step = 1, max(in%n_source_updates, 1)
        if (in%boole_static_ne) call calc_starting_conditions   !fresh weights (source may have changed since last iter)
        call write_particle_source(source_step)
        call run_scef_iterations(source_step)
        call update_particle_source
    enddo
    call write_particle_source(max(in%n_source_updates, 1) + 1)   !state after the final update

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
    print*, 'Average abs Delta Phi at all the electric potential updates = ', ep%average_abs_phi_elec_from_rho
    print*, sum(weights%w(:,1))/(in%num_particles*in%density*sum(output%prism_volumes(:)))
    print*, sum(output%prism_volumes(:)*output%prism_moments(1,:,1))/(sum(output%prism_volumes(:))*in%density)

end subroutine calc_self_consistent_electric_field

subroutine run_scef_iterations(source_step)

    !One round of self-consistent electric-field iterations: for the current
    !source and particle initialisation, push both species n_electric_potential_updates
    !times and update the electric potential each time. source_step is the source-update
    !iteration index, passed to perform_electric_potential_update for file naming.

    use gorilla_applets_types_mod, only: in, ep, output, start, exit_data
    use utils_data_pre_and_post_processing_mod, only: prepare_next_round_of_parallelised_particle_pushing, &
    calc_collision_coefficients_for_all_tetrahedra, normalise_prism_moments_and_prism_moments_squared
    use utils_scef_particle_init_mod, only: calc_starting_conditions
    use utils_scef_particle_pushing_mod, only: parallelised_particle_pushing
    use utils_scef_electric_potential_mod, only: perform_electric_potential_update
    use utils_scef_diffusion_mod, only: calc_density_via_random_walk

    integer, intent(in) :: source_step
    integer :: i, species

    do i = 1, max(in%n_electric_potential_updates,1)
        if (.not.in%boole_static_ne) call calc_starting_conditions
        do species = 1,2 !trace electrons and ions
            if ((species.eq.2).and.(in%boole_static_ne)) cycle
            call prepare_next_round_of_parallelised_particle_pushing(species)
            if (.not.in%boole_honest_tracing(species)) then
                call calc_density_via_random_walk(species, i)
            else
                if (in%boole_collisions) call calc_collision_coefficients_for_all_tetrahedra(species)
                call parallelised_particle_pushing(species,i,boole_diffusion_coefficient=.false.)
                call normalise_prism_moments_and_prism_moments_squared(species)
                ep%mean_exit_time(species) = sum(exit_data%t_confined(:, species))/in%num_particles
            endif
        enddo

        !prism_moments are normalised by in%time_step, but each particle only lives its own exit time.
        !Rescale by in%time_step / max(mean_exit_time) uniformly across both species so the resulting
        !densities represent what a continuous source would produce at steady state. Uniform scaling
        !(not per-species) is essential: a per-species correction would force n_i = n_e and short-circuit
        !the ambipolar feedback that drives the self-consistent electric field.
        ep%exit_time_correction = 1.0_dp
        if (max(ep%mean_exit_time(1), ep%mean_exit_time(2)) > 0.0_dp) then
            ep%exit_time_correction = in%time_step / max(ep%mean_exit_time(1), ep%mean_exit_time(2))
        endif
        ep%rho_prism = 0
        do species = 1,2
            if ((species.eq.2).and.(in%boole_static_ne)) cycle
            output%prism_moments(1,:,species) = output%prism_moments(1,:,species) * ep%exit_time_correction
            ep%rho_prism = ep%rho_prism + real(output%prism_moments(1,:,species))*start%particle_charge(species)
        enddo

        call perform_electric_potential_update(i, source_step)
    enddo

end subroutine run_scef_iterations

end module self_consistent_electric_field_mod