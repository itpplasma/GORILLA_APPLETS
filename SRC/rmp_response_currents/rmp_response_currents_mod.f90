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
        fourier_transform_moments
    use utils_helical_core_mod, only: eliminate_particles_outside_flux_threshold
    use utils_rmp_response_currents_mod, only: read_rmp_response_currents_inp_into_type, &
        parallelised_particle_pushing_rmp_response_currents, &
        calc_starting_conditions_rmp_response_currents, &
        profile_dir, equil_mapping_file, boole_constant_delta_B_r, delta_B_r_const, &
        pert_m_fourier, pert_n_fourier, delta_B_r_file, boole_skip_phase_for_test, &
        boole_step_delta_B_r, step_center_reff, step_halfwidth_reff, &
        bias_starting_positions_to_s_window, dump_start_positions, &
        spawn_equidistant_in_s, boole_equidistant_s_sampling, &
        boole_dump_orbit_n1, traj_dump_unit, traj_step_count, &
        boole_dump_collisions_n1, coll_dump_unit, coll_event_count, &
        coll_dt_sum, coll_dist_sum, &
        boole_use_kim_nu, kim_nu_file, &
        compute_spawn_volume, filter_markers_by_trapping, trapping_filter_mode
    use profile_data_mod, only: load_profiles, load_kim_nu
    use perturbation_field_mod, only: init_constant_perturbation, load_perturbation_field, &
                                      init_step_perturbation, boole_skip_phase

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
        if (boole_use_kim_nu) then
            if (len_trim(kim_nu_file) == 0) then
                print *, 'ERROR: boole_use_kim_nu=.true. but kim_nu_file is empty'
                stop
            end if
            call load_kim_nu(trim(kim_nu_file))
        end if
        if (boole_step_delta_B_r) then
            call init_step_perturbation(step_center_reff, step_halfwidth_reff, &
                                        delta_B_r_const, pert_m_fourier, pert_n_fourier, &
                                        trim(equil_mapping_file))
            boole_skip_phase = boole_skip_phase_for_test
        else if (boole_constant_delta_B_r) then
            call init_constant_perturbation(delta_B_r_const, pert_m_fourier, pert_n_fourier)
            boole_skip_phase = boole_skip_phase_for_test
        else
            call load_perturbation_field(trim(delta_B_r_file), trim(equil_mapping_file), &
                                         pert_m_fourier, pert_n_fourier)
        end if
    end if

    call calc_starting_conditions_rmp_response_currents
    call eliminate_particles_outside_flux_threshold
    if (in%boole_delta_f) then
        if (boole_equidistant_s_sampling) then
            block
                integer :: n_spawned
                call spawn_equidistant_in_s(species=1, n_spawned=n_spawned)
            end block
        else
            call bias_starting_positions_to_s_window
        end if
    end if

    ! Optional trapped/passing population filter (post-spawn).
    if (trim(trapping_filter_mode) /= 'none') then
        call filter_markers_by_trapping(species=1)
    end if

    ! Diagnostic: dump per-particle initial positions to start_positions.dat.
    call dump_start_positions('start_positions.dat', species=1)

    ! Diagnostic: open trajectory dump for marker n=1 (q-vs-s comparison).
    if (boole_dump_orbit_n1) then
        open(newunit=traj_dump_unit, file='traj_n1.dat', status='unknown', action='write')
        write(traj_dump_unit, '(a)') '# t  R  phi  Z  vpar  Re(w)  Im(w)  Re(w*vpar)  Im(w*vpar)  ind_tetr  iper_phi'
        traj_step_count = 0
    end if

    ! Diagnostic: open collision-event dump for marker n=1.
    if (boole_dump_collisions_n1) then
        open(newunit=coll_dump_unit, file='coll_n1.dat', status='unknown', action='write')
        write(coll_dump_unit, '(a)') '# event  t  ind_tetr  R  phi  Z  vpar  vperp  v_tot  dt_to_next_collision'
        coll_event_count = 0
    end if

    call parallelised_particle_pushing_rmp_response_currents(species=1)

    if (boole_dump_collisions_n1 .and. coll_dump_unit /= 0) then
        close(coll_dump_unit)
        coll_dump_unit = 0
        print*, ''
        print*, '=== Collision-event diagnostics for marker n=1 ==='
        print*, 'Total collision events: ', coll_event_count
        if (coll_event_count > 0) then
            print*, 'Mean dt between collisions [s]   : ', coll_dt_sum   / real(coll_event_count, 8)
            print*, 'Mean distance between events [cm]: ', coll_dist_sum / real(coll_event_count, 8)
        end if
        print*, ''
    end if

    if (boole_dump_orbit_n1 .and. traj_dump_unit /= 0) then
        close(traj_dump_unit)
        traj_dump_unit = 0
    end if

    ! For delta-f the per-particle time average is already applied at
    ! deposit time inside orbit_timestep_rmp_response_currents (each
    ! contribution is divided by that marker's own T_n = N_collisions *
    ! tau_c(s_n), since the trace time varies across the population). So
    ! we skip the global time divisor here -- otherwise we would divide by
    ! in%time_step a second time. The volume / N_particles part of the
    ! normalisation still runs.
    call normalise_prism_moments_and_prism_moments_squared( &
        boole_skip_time_normalisation_in=in%boole_delta_f)

    ! Apply the spawn-window volume V_W. With markers placed uniformly in
    ! (R, phi, Z) inside the spawn region, each marker represents
    ! (n(s) * V_W / N_p) physical particles, so the per-cell density
    ! moment carries an extra factor V_W relative to what the deposit
    ! sums to. See docs/2026-05-06-spawn-volume-normalisation.md.
    if (in%boole_delta_f) then
        block
            real(dp) :: V_W
            V_W = compute_spawn_volume()
            output%prism_moments(:,:,:) = output%prism_moments(:,:,:) * V_W
            if (moment_specs%boole_squared_moments) then
                output%prism_moments_squared(:,:,:) = output%prism_moments_squared(:,:,:) * V_W**2
            end if
            print '(a, es12.5, a)', ' Spawn-window volume V_W = ', V_W, &
                                     ' cm^3 (applied as moment multiplier)'
        end block
    end if

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
