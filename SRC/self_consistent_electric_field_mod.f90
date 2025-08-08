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
    use gorilla_applets_types_mod, only: moment_specs, counter, c, in, start, s, output
    use utils_write_data_to_files_mod, only: write_data_to_files, give_file_names, unlink_files
    use utils_data_pre_and_post_processing_mod, only: set_seed_for_random_numbers, &
    get_ipert, set_moment_specifications, initialise_output, initialize_exit_data, calc_poloidal_flux, &
    calc_collision_coefficients_for_all_tetrahedra, normalise_prism_moments_and_prism_moments_squared, fourier_transform_moments, &
    find_minimal_angle_between_curlA_and_tetrahedron_faces, analyse_particle_weight_distribution, &
    set_weights, prepare_next_round_of_parallelised_particle_pushing
    use utils_self_consistent_ef_mod, only: allocate_electric_potential_type, perform_electric_potential_update, &
    associate_flux_labels_with_tetrahedra_and_vertices, print_errors_for_bad_inputs, &
    read_self_consistent_electric_field_inp_into_type, calc_starting_conditions, calc_electron_diffusion_coefficients, &
    parallelised_particle_pushing, calc_s_shell_volumes
    use gorilla_applets_types_mod, only: output, ep
    use tetra_physics_mod, only: particle_mass

    integer :: i, species

    call set_seed_for_random_numbers
    call read_self_consistent_electric_field_inp_into_type
    call get_ipert()
    call initialize_gorilla(i_option,ipert)
    
    call set_moment_specifications
    call initialise_output
    call calc_volume_integrals_in_flux_coordinates

    if (.not.in%boole_static_ne) call calc_electron_diffusion_coefficients

    call initialize_exit_data
    call associate_flux_labels_with_tetrahedra_and_vertices
    call calc_poloidal_flux(verts_sthetaphi)
    call allocate_electric_potential_type
    call calc_s_shell_volumes
    call give_file_names
    call unlink_files
    call print_errors_for_bad_inputs    

    !call perform_electric_potential_update(0)
    do i = 1, max(in%n_electric_potential_updates,1)
        call calc_starting_conditions
        ep%rho_prism = 0
        do species = 1,2 !trace electrons and ions
            if ((species.eq.2).and.(in%boole_static_ne)) cycle
            call prepare_next_round_of_parallelised_particle_pushing(species)
            if (in%boole_collisions) call calc_collision_coefficients_for_all_tetrahedra(species)
            call parallelised_particle_pushing(species,i,boole_diffusion_coefficient=.false.)
            call normalise_prism_moments_and_prism_moments_squared(species)
            ep%rho_prism = ep%rho_prism + real(output%prism_moments(1,:,species))*start%particle_charge(species)
        enddo
        call perform_electric_potential_update(i)
        print*, 'average ion weight is ', sum(start%weight(:,1))/in%n_particles
        print*, 'average electron weight is ', sum(start%weight(:,2))/in%n_particles
    enddo

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
    print*, sum(start%weight(:,1))/(in%num_particles*in%density*sum(output%prism_volumes(:)))
    print*, sum(output%prism_volumes(:)*output%prism_moments(1,:,1))/(sum(output%prism_volumes(:))*in%density)

end subroutine calc_self_consistent_electric_field

end module self_consistent_electric_field_mod