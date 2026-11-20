module helical_core_mod
!
! Module for helical core particle tracing.
!
! This is a simplified version of anomalous_transport_mod that:
!   - Traces particles in a helical core geometry
!   - Does NOT include anomalous transport displacement (straight-line displacements)
!   - Does NOT include electric potential perturbation scans
!
! The module uses the same infrastructure as anomalous transport but removes
! all anomalous-specific features.
!

    use, intrinsic :: iso_fortran_env, only: dp => real64

    implicit none

contains

subroutine calc_helical_core

    use utils_data_pre_and_post_processing_mod, only: setup_applet, print_applet_summary
    use tetra_grid_mod, only: verts_rphiz
    use volume_integrals_and_sqrt_g_mod, only: calc_square_root_g, calc_volume_integrals
    use gorilla_applets_types_mod, only: moment_specs, in, start, output, weights, c
    use utils_write_data_to_files_mod, only: write_data_to_files, give_file_names, unlink_files
    use utils_data_pre_and_post_processing_mod, only: set_seed_for_random_numbers, initialize_exit_data, &
    calc_poloidal_flux, calc_collision_coefficients_for_all_tetrahedra, build_reservoir_cells, set_reservoir_capacities, &
    normalise_prism_moments_and_prism_moments_squared, fourier_transform_moments, calc_starting_conditions
    use utils_parallelised_particle_pushing_mod, only: start_deferred_background_updates, &
        apply_deferred_background_updates
    use utils_helical_core_mod, only: read_helical_core_inp_into_type, &
        parallelised_particle_pushing_helical_core, eliminate_particles_outside_flux_threshold, &
        set_unidirectional_beam, write_collision_relaxation_diagnostics, write_collision_relaxation_profiles, &
        write_marker_weights, boole_df_time_dependent, df_burn_in_time, df_time_offset, w_td, w1_td
    use delta_f_restoring_mod, only: boole_restoring, setup_restoring, compute_restoring_coefficients, rest_loss, &
        write_restoring_diagnostics, n_rest_exact, n_rest_fallback, n_rest_empty, rest_given

    integer :: i_slice, relaxation_unit, n_rejected, restoring_unit
    real(dp) :: dp_tot, de_tot, dp_abs
    real(dp), dimension(:), allocatable :: base_w, w1_use
    real(dp) :: max_delta_u, max_delta_temp_rel
    logical :: boole_relaxation_diagnostics
    real(dp), dimension(:), allocatable :: pitch_0, energy_0
    real(dp), dimension(:,:), allocatable :: x_0
    logical, dimension(:), allocatable :: lost_0
    real(dp), dimension(:,:), allocatable :: u_b_0, temp_b_0

    call set_seed_for_random_numbers
    call read_helical_core_inp_into_type
    call setup_applet

    call calc_square_root_g
    call calc_volume_integrals

    call initialize_exit_data
    call calc_poloidal_flux(verts_rphiz)
    call give_file_names
    call unlink_files

    if (in%boole_collisions) call calc_collision_coefficients_for_all_tetrahedra
    if (in%background_temperature_eV.gt.0.0_dp) then
        !the conserving collisions read temp_mat at every collision; the other operator uses precomputed coefficients
        if (.not.(in%boole_collisions.and.in%boole_preserve_energy_and_momentum_during_collisions)) then
            print*, 'Error: background_temperature_eV requires boole_preserve_energy_and_momentum_during_collisions'
            stop
        endif
        if (in%boole_profiles_from_vertex_data) then
            print*, 'Error: background_temperature_eV would overwrite the temperatures of boole_profiles_from_vertex_data'
            stop
        endif
        c%temp_mat = in%background_temperature_eV
    endif

    call calc_starting_conditions
    call eliminate_particles_outside_flux_threshold
    call set_unidirectional_beam

    ! Standard particle tracing (no scan, no diffusion coefficient calculation). With n_time_slices > 1,
    ! all markers are pushed slice by slice. The background reservoir of the conserving collisions is
    ! updated after every slice, so markers and background evolve in common time (needed for a
    ! relaxation in time), and for a stationary state each slice is one step of the iteration between
    ! markers and background.
    boole_relaxation_diagnostics = in%boole_collisions.and.in%boole_preserve_energy_and_momentum_during_collisions
    if (boole_relaxation_diagnostics) then
        x_0 = start%x(:,:,1)
        pitch_0 = start%pitch(:,1)
        energy_0 = start%energy(:,1)
        call build_reservoir_cells(in%reservoir_min_volume_fraction)
        call set_reservoir_capacities(in%reservoir_capacity_flow, in%reservoir_capacity_temperature, in%n_particles)
        if (in%n_time_slices.eq.1) print*, 'WARNING: n_time_slices = 1: the background reservoir is updated only ', &
            'after all markers are traced, so the markers do not feel the momentum and energy exchange'
        call start_deferred_background_updates
        lost_0 = start%lost(:,1)
        u_b_0 = c%vpar_mat
        temp_b_0 = c%temp_mat
        open(newunit=relaxation_unit, file='collision_relaxation.dat', status='replace')
        write(relaxation_unit,'(A)') '#time <v_par>_m/v0 <u_par>_b/v0 P_m P_b P_lost E_m E_b E_lost E_tot ' // &
                                     'L_m L_b L_lost L_tot <T>_b T_m min_T_b max|u_b|/v0 n_confined <T>_e min_T_e'
    endif
    if (boole_df_time_dependent) then
        allocate(w_td(in%num_particles), w1_td(in%num_particles))
        w_td = 0.0_dp
        w1_td = 0.0_dp
        print*, 'time-dependent delta-f weights: ', in%n_time_slices, ' slices, burn-in ', df_burn_in_time, &
                ' s of ', in%time_step, ' s'
        if (boole_restoring) then
            call setup_restoring(in%num_particles)
            allocate(base_w(in%num_particles), w1_use(in%num_particles))
            open(newunit=restoring_unit, file='delta_f_restoring.dat', status='replace')
            write(restoring_unit,'(A)') '#slice time sum_dp sum_dE sum|dp| given_dp given_dE n_bins_exact n_bins_fallback n_bins_empty'
        endif
    endif
    start%t = in%time_step/in%n_time_slices
    do i_slice = 1, in%n_time_slices
        df_time_offset = (i_slice - 1)*start%t(1)
        call parallelised_particle_pushing_helical_core(species=1, boole_reuse_weights_in=(i_slice.gt.1), i_slice_in=i_slice)
        if (boole_restoring) then
            !momentum and energy lost by delta-f in this slice are returned during the next slice
            dp_tot = sum(rest_loss(1, :))
            de_tot = sum(rest_loss(2, :))
            dp_abs = sum(abs(rest_loss(1, :)))
            base_w = real(weights%original(:, 1), dp)
            w1_use = 0.0_dp
            if (allocated(w1_td)) w1_use = w1_td
            call compute_restoring_coefficients(start%t(1), 1, base_w, w1_use)
            call write_restoring_diagnostics(restoring_unit, i_slice, i_slice*start%t(1), dp_tot, de_tot, dp_abs)
            rest_given = 0.0_dp
            print*, 'slice ', i_slice, ': restoring bins exact/fallback/empty = ', n_rest_exact, n_rest_fallback, &
                    n_rest_empty
        endif
        if (boole_relaxation_diagnostics) then
            call apply_deferred_background_updates(n_rejected, max_delta_u, max_delta_temp_rel)
            if (n_rejected.gt.0) print*, 'WARNING: background updates rejected (T <= 0) in ', n_rejected, &
                                         ' tetrahedra in slice ', i_slice, ', energy and momentum not conserved there'
            print*, 'slice ', i_slice, ': largest background change of a cell: |du_par|/v0 = ', &
                    max_delta_u/start%v0(1), ', |dT|/T = ', max_delta_temp_rel
            !weights are known after the first slice and stay fixed, so the initial row can be written now
            if (i_slice.eq.1) call write_collision_relaxation_diagnostics(relaxation_unit, 0.0_dp, x_0, pitch_0, energy_0, &
                                                                          lost_0, u_b_0, temp_b_0)
            call write_collision_relaxation_diagnostics(relaxation_unit, i_slice*start%t(1), start%x(:,:,1), start%pitch(:,1), &
                                                        start%energy(:,1), start%lost(:,1), c%vpar_mat, c%temp_mat)
        endif
    enddo
    start%t = in%time_step
    if (boole_restoring) close(restoring_unit)
    if (boole_relaxation_diagnostics) then
        close(relaxation_unit)
        call write_collision_relaxation_profiles
    endif
    call write_marker_weights

    if (boole_df_time_dependent) then
        !time average over the deposit window [burn-in, time_step]: the moments are then densities (as in full-f)
        start%t = in%time_step - df_burn_in_time
        call normalise_prism_moments_and_prism_moments_squared(boole_skip_time_normalisation_in=.false.)
        start%t = in%time_step
    else
        call normalise_prism_moments_and_prism_moments_squared(boole_skip_time_normalisation_in=in%boole_delta_f)
    endif

    if (moment_specs%n_moments.gt.0) call fourier_transform_moments
    call write_data_to_files

    call print_applet_summary('ion', start%particle_mass(1), start%v0(1), start%particle_charge(1))

end subroutine calc_helical_core

end module helical_core_mod
