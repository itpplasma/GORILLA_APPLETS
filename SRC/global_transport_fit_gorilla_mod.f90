module global_transport_fit_gorilla_mod

    use, intrinsic :: iso_fortran_env, only: dp => real64
    use transport_statistics_mod, only: sanitize_transport_output_triplet

    implicit none

    private

    public :: run_gorilla_mc_global_transport_fit
    public :: select_local_reference_boundaries

contains

subroutine run_gorilla_mc_global_transport_fit(control, result)

    use field_mod, only: ipert
    use global_transport_fit_core_mod, only: fit_global_transport
    use global_transport_fit_io_mod, only: write_convergence_history, write_global_transport_fit_outputs
    use global_transport_fit_math_mod, only: compute_geometry_from_boundaries
    use global_transport_fit_settings_mod, only: filename_comparison_profiles, filename_convergence_history, &
        filename_experiment_1_summary, filename_experiment_2_summary, filename_fit_profiles, filename_fit_summary, &
        filename_local_profiles, local_density_gradient_1, local_density_gradient_2, local_n_particles, &
        local_total_time, local_v_E, max_source_trials, min_density_source_relstd, min_supported_flux_boundaries, &
        n_global_batches, n_local_batches, n_local_gradients, n_local_reference_surfaces, source_gradient_1, &
        source_gradient_2, source_reference_s_1, source_reference_s_2, source_time_growth, source_width_1, source_width_2
    use global_transport_fit_types_mod, only: global_transport_experiment_t, global_transport_fit_control_t, &
        global_transport_fit_result_t, global_transport_quality_t
    use gorilla_applets_settings_mod, only: i_option
    use gorilla_applets_types_mod, only: in, s
    use orbit_timestep_gorilla_mod, only: initialize_gorilla
    use transport_benchmark_utils_mod, only: prepare_minimal_transport_output, set_active_species_parameters
    use utils_data_pre_and_post_processing_mod, only: get_ipert, set_seed_for_random_numbers
    use utils_self_consistent_ef_mod, only: allocate_electric_potential_type, associate_flux_labels_with_tetrahedra_and_vertices, &
        calc_s_shell_volumes, print_errors_for_bad_inputs, read_self_consistent_electric_field_inp_into_type
    use volume_integrals_and_sqrt_g_mod, only: calc_volume_integrals_in_flux_coordinates

    type(global_transport_fit_control_t), intent(in) :: control
    type(global_transport_experiment_t), allocatable :: experiments(:)
    type(global_transport_fit_result_t), intent(out) :: result
    type(global_transport_quality_t), allocatable :: qualities(:)
    real(dp), allocatable :: boundary_areas(:)
    real(dp), allocatable :: boundary_s(:)
    real(dp), allocatable :: cell_centers(:)
    real(dp), allocatable :: local_a(:)
    real(dp), allocatable :: local_a_std(:)
    real(dp), allocatable :: local_b(:)
    real(dp), allocatable :: local_b_std(:)
    real(dp), allocatable :: local_boundary_s(:)
    real(dp), allocatable :: shell_volumes(:)
    logical, allocatable :: local_valid(:)
    integer, allocatable :: local_boundary_indices(:)
    real(dp) :: base_time_step
    real(dp) :: global_tracing_time
    real(dp) :: local_time_step
    integer :: local_particle_count

    call set_seed_for_random_numbers()
    call read_self_consistent_electric_field_inp_into_type()
    in%tracer_species = 1
    in%n_species = 1
    in%boole_custom_background = .true.
    in%background_density_cm3 = (/in%density, in%density/)
    in%background_temperature_eV = (/in%energy_eV, in%energy_eV/)
    in%collision_operator = 4
    s%temperature = in%energy_eV
    call get_ipert()
    call set_active_species_parameters(in%tracer_species)
    call initialize_gorilla(i_option, ipert)
    call set_active_species_parameters(in%tracer_species)
    call prepare_minimal_transport_output()
    call calc_volume_integrals_in_flux_coordinates()
    call associate_flux_labels_with_tetrahedra_and_vertices()
    call allocate_electric_potential_type()
    call calc_s_shell_volumes()
    call print_errors_for_bad_inputs()

    base_time_step = in%time_step
    call extract_boundary_geometry(boundary_s, shell_volumes)
    call compute_geometry_from_boundaries(boundary_s, shell_volumes, cell_centers, boundary_areas)

    allocate(experiments(2))
    allocate(qualities(2))
    call run_electron_source_experiment(source_gradient_1, source_reference_s_1, boundary_s, shell_volumes, boundary_areas, &
        source_width_1, base_time_step, n_global_batches, max_source_trials, source_time_growth, &
        min_supported_flux_boundaries, min_density_source_relstd, control%support_sigma_multiplier, &
        filename_experiment_1_summary, experiments(1), qualities(1))
    call run_electron_source_experiment(source_gradient_2, source_reference_s_2, boundary_s, shell_volumes, boundary_areas, &
        source_width_2, base_time_step, n_global_batches, max_source_trials, source_time_growth, &
        min_supported_flux_boundaries, min_density_source_relstd, control%support_sigma_multiplier, &
        filename_experiment_2_summary, experiments(2), qualities(2))

    call fit_global_transport(experiments, control, result)
    call write_global_transport_fit_outputs(boundary_s, result, filename_fit_summary, filename_fit_profiles)
    call write_convergence_history(result, filename_convergence_history)

    global_tracing_time = max(qualities(1)%tracing_time, qualities(2)%tracing_time)
    local_time_step = global_tracing_time
    if (local_total_time > 0.0_dp) local_time_step = local_total_time
    local_particle_count = in%num_particles
    if (local_n_particles > 0) local_particle_count = local_n_particles
    call build_local_transport_reference(boundary_s, n_local_reference_surfaces, local_density_gradient_1, &
        local_density_gradient_2, n_local_gradients, local_v_E, local_time_step, &
        local_particle_count, n_local_batches, local_boundary_indices, local_boundary_s, local_a, local_a_std, local_b, &
        local_b_std, local_valid)
    call write_local_transport_profiles(local_boundary_indices, local_boundary_s, local_valid, local_a, local_a_std, local_b, &
        local_b_std, filename_local_profiles)
    call write_transport_comparison(local_boundary_indices, local_boundary_s, local_valid, local_a, local_a_std, local_b, &
        local_b_std, result%a_profile, result%a_std, result%b_profile, result%b_std, filename_comparison_profiles)

    in%time_step = base_time_step

end subroutine run_gorilla_mc_global_transport_fit

subroutine extract_boundary_geometry(boundary_s, shell_volumes)

    use gorilla_applets_types_mod, only: ep
    use tetra_grid_mod, only: verts_sthetaphi
    use tetra_grid_settings_mod, only: grid_size

    real(dp), allocatable, intent(out) :: boundary_s(:)
    real(dp), allocatable, intent(out) :: shell_volumes(:)
    integer :: i

    allocate(boundary_s(grid_size(1) + 1))
    allocate(shell_volumes(grid_size(1)))
    do i = 1, grid_size(1) + 1
        boundary_s(i) = verts_sthetaphi(1, grid_size(3) * (i - 1) + 1)
    end do
    shell_volumes = ep%s_shell_volumes

end subroutine extract_boundary_geometry

subroutine run_electron_source_experiment(source_gradient, source_reference_s, boundary_s, shell_volumes, boundary_areas, &
    source_width, initial_total_time, n_batches, max_trials, time_growth, min_supported_flux_boundaries, &
    min_density_source_relstd, &
    sigma_multiplier, summary_filename, experiment, quality)

    use global_transport_fit_types_mod, only: global_transport_experiment_t, global_transport_quality_t

    real(dp), intent(in) :: boundary_areas(:)
    real(dp), intent(in) :: boundary_s(:)
    real(dp), intent(in) :: initial_total_time
    real(dp), intent(in) :: min_density_source_relstd
    real(dp), intent(in) :: shell_volumes(:)
    real(dp), intent(in) :: sigma_multiplier
    real(dp), intent(in) :: source_gradient
    real(dp), intent(in) :: source_reference_s
    real(dp), intent(in) :: source_width
    real(dp), intent(in) :: time_growth
    integer, intent(in) :: max_trials
    integer, intent(in) :: min_supported_flux_boundaries
    integer, intent(in) :: n_batches
    character(len=*), intent(in) :: summary_filename
    type(global_transport_experiment_t), intent(out) :: experiment
    type(global_transport_quality_t), intent(out) :: quality

    type(global_transport_experiment_t), allocatable :: batch_experiments(:)
    real(dp), allocatable :: batch_lost_fraction(:)
    real(dp), allocatable :: batch_tracing_time(:)
    real(dp) :: current_total_time
    integer :: i_batch
    integer :: i_trial

    current_total_time = initial_total_time
    do i_trial = 1, max_trials
        allocate(batch_experiments(n_batches))
        allocate(batch_lost_fraction(n_batches))
        allocate(batch_tracing_time(n_batches))
        do i_batch = 1, n_batches
            call run_electron_source_experiment_batch(source_gradient, source_reference_s, boundary_s, shell_volumes, &
                boundary_areas, source_width, current_total_time, batch_experiments(i_batch), batch_lost_fraction(i_batch), &
                batch_tracing_time(i_batch))
        end do

        call aggregate_batched_experiments(batch_experiments, experiment)
        call evaluate_transport_signal(experiment, batch_lost_fraction, batch_tracing_time, sigma_multiplier, &
            min_supported_flux_boundaries, min_density_source_relstd, i_trial, quality)
        print '(A,1X,I0,1X,A,1X,ES12.4,1X,A,1X,ES12.4,1X,A,1X,I0,1X,A,1X,L1)', &
            'global_transport_trial', i_trial, 'time', quality%tracing_time, 'density_source_relstd', &
            quality%density_source_relstd, 'supported_flux_boundaries', quality%n_supported_flux_boundaries, 'accepted', &
            quality%transport_supported
        if (quality%transport_supported) exit

        current_total_time = current_total_time * time_growth
        deallocate(batch_experiments, batch_lost_fraction, batch_tracing_time)
    end do

    call write_experiment_summary(boundary_s, shell_volumes, boundary_areas, experiment, quality, summary_filename, &
        source_gradient, source_reference_s, sigma_multiplier)
    if (.not.quality%transport_supported) then
        print *, 'Error: global transport source experiment did not develop enough radial transport signal.'
        print *, 'source_gradient = ', source_gradient
        print *, 'source_reference_s = ', source_reference_s
        stop
    end if
    deallocate(batch_experiments, batch_lost_fraction, batch_tracing_time)

end subroutine run_electron_source_experiment

subroutine run_electron_source_experiment_batch(source_gradient, source_reference_s, boundary_s, shell_volumes, boundary_areas, &
    source_width, total_time, experiment, lost_fraction, mean_tracing_time)

    use global_transport_fit_types_mod, only: global_transport_experiment_t, global_transport_samples_t
    use gorilla_applets_types_mod, only: in, output, start
    use start_position_sampling_mod, only: set_random_start_positions, validate_current_start_positions
    use tetra_grid_settings_mod, only: grid_size, sfc_s_min
    use utils_data_pre_and_post_processing_mod, only: calc_collision_coefficients_for_all_tetrahedra, initialize_exit_data, &
        prepare_next_round_of_parallelised_particle_pushing
    use utils_self_consistent_ef_mod, only: allocate_start_type, parallelised_particle_pushing, &
        set_particle_type_specifications, set_rest_of_individual_particle_specifications

    real(dp), intent(in) :: boundary_areas(:)
    real(dp), intent(in) :: boundary_s(:)
    real(dp), intent(in) :: shell_volumes(:)
    real(dp), intent(in) :: source_gradient
    real(dp), intent(in) :: source_reference_s
    real(dp), intent(in) :: source_width
    real(dp), intent(in) :: total_time
    type(global_transport_experiment_t), intent(out) :: experiment
    real(dp), intent(out) :: lost_fraction
    real(dp), intent(out) :: mean_tracing_time

    type(global_transport_samples_t) :: samples
    real(dp), allocatable :: rand_matrix(:, :, :)
    integer :: failed_positions
    real(dp) :: initial_weight_sum
    real(dp) :: source_s_max
    real(dp) :: source_s_min

    in%time_step = total_time
    allocate(samples%source_weight_sum(grid_size(1)))
    allocate(samples%shell_time_sum(grid_size(1)))
    allocate(samples%boundary_weighted_flux_sum(grid_size(1) + 1))
    samples%source_weight_sum = 0.0_dp
    samples%shell_time_sum = 0.0_dp
    samples%boundary_weighted_flux_sum = 0.0_dp

    allocate(rand_matrix(5, in%num_particles, 1))
    call random_number(rand_matrix)

    output%radial_flux = 0.0_dp
    output%weighted_radial_flux = 0.0_dp

    in%boole_custom_source_profile = source_width > 0.0_dp
    in%source_profile_reference_s = source_reference_s
    in%source_profile_width = source_width
    if (in%boole_custom_source_profile) then
        in%density_log_gradient_per_s = 0.0_dp
    else
        in%density_log_gradient_per_s = source_gradient
    end if
    in%density_profile_reference_s = source_reference_s

    call allocate_start_type(in%num_particles)
    call set_particle_type_specifications()
    call initialize_exit_data(in%num_particles)
    source_s_min = sfc_s_min
    source_s_max = 1.0_dp
    if (source_width > 0.0_dp) then
        source_s_min = max(sfc_s_min, source_reference_s - 4.0_dp * source_width)
        source_s_max = min(1.0_dp, source_reference_s + 4.0_dp * source_width)
    end if
    call set_random_start_positions(rand_matrix, species_in=(/1/), s_min_in=source_s_min, s_max_in=source_s_max)
    call set_rest_of_individual_particle_specifications(rand_matrix, boole_diffusion_coefficient_in=.true., &
        species_in=(/1/), n_particles_in=in%num_particles)
    call validate_current_start_positions(species_in=(/1/), s_min_in=source_s_min, s_max_in=source_s_max, &
        max_validation_attempts=256, failed_positions=failed_positions)
    if (failed_positions /= 0) then
        print *, 'Error: unable to initialize validated source starts on the GORILLA mesh.'
        print *, 'failed_positions = ', failed_positions
        print *, 'source_reference_s = ', source_reference_s
        print *, 'source_width = ', source_width
        stop
    end if
    initial_weight_sum = sum(start%weight(:, 1))
    call prepare_next_round_of_parallelised_particle_pushing(1)
    if (in%boole_collisions) call calc_collision_coefficients_for_all_tetrahedra(1)
    call parallelised_particle_pushing(species=1, j=1, boole_diffusion_coefficient=.false., n_particles_in=in%num_particles, &
        transport_samples=samples)
    if (sum(samples%source_weight_sum) <= 0.0_dp) then
        print *, 'Error: global transport sampling produced zero source support.'
        print *, 'source_reference_s = ', source_reference_s
        print *, 'source_width = ', source_width
        print *, 'initial_weight_sum = ', initial_weight_sum
        print *, 'sum_start_weight = ', sum(start%weight(:, 1))
        print *, 'min_start_weight = ', minval(start%weight(:, 1))
        print *, 'max_start_weight = ', maxval(start%weight(:, 1))
        print *, 'lost_particles = ', samples%lost_particles
        print *, 'mean_tracing_time = ', samples%mean_tracing_time
        stop
    end if
    if (sum(samples%shell_time_sum) <= 0.0_dp) then
        print *, 'Error: global transport sampling produced zero residence time.'
        print *, 'source_reference_s = ', source_reference_s
        print *, 'source_width = ', source_width
        print *, 'sum_source_weight = ', sum(samples%source_weight_sum)
        print *, 'lost_particles = ', samples%lost_particles
        print *, 'mean_tracing_time = ', samples%mean_tracing_time
        stop
    end if

    call build_experiment_from_samples(boundary_s, shell_volumes, boundary_areas, samples, experiment)
    mean_tracing_time = samples%mean_tracing_time
    lost_fraction = real(samples%lost_particles, dp) / real(max(samples%n_particles, 1), dp)

end subroutine run_electron_source_experiment_batch

subroutine build_experiment_from_samples(boundary_s, shell_volumes, boundary_areas, samples, experiment)

    use global_transport_fit_types_mod, only: global_transport_experiment_t, global_transport_samples_t

    real(dp), intent(in) :: boundary_areas(:)
    real(dp), intent(in) :: boundary_s(:)
    real(dp), intent(in) :: shell_volumes(:)
    type(global_transport_samples_t), intent(in) :: samples
    type(global_transport_experiment_t), intent(out) :: experiment

    real(dp), allocatable :: safe_boundary_areas(:)
    real(dp), allocatable :: safe_shell_volumes(:)
    real(dp) :: normalization

    normalization = max(sum(samples%source_weight_sum), 1.0d-20)

    experiment%boundary_s = boundary_s
    experiment%shell_volumes = shell_volumes
    allocate(safe_shell_volumes(size(shell_volumes)))
    safe_shell_volumes = max(shell_volumes, 1.0d-20)
    experiment%source = samples%source_weight_sum / (normalization * safe_shell_volumes)
    experiment%density = samples%shell_time_sum / (normalization * safe_shell_volumes)
    allocate(experiment%density_variance(size(shell_volumes)))
    experiment%density_variance = 1.0d-30
    allocate(experiment%integrated_flux(size(boundary_s)))
    experiment%integrated_flux = samples%boundary_weighted_flux_sum / normalization
    allocate(experiment%integrated_flux_variance(size(boundary_s)))
    experiment%integrated_flux_variance = 1.0d-30
    allocate(experiment%flux(size(boundary_s)))
    experiment%flux = 0.0_dp
    allocate(safe_boundary_areas(size(boundary_areas)))
    safe_boundary_areas = max(boundary_areas, 1.0d-20)
    experiment%flux(2:) = experiment%integrated_flux(2:) / safe_boundary_areas(2:)
    allocate(experiment%flux_variance(size(boundary_s)))
    experiment%flux_variance = 1.0d-30
    experiment%flux_variance(1) = 0.0_dp

end subroutine build_experiment_from_samples

subroutine aggregate_batched_experiments(batch_experiments, experiment)

    use global_transport_fit_types_mod, only: global_transport_experiment_t
    use transport_statistics_mod, only: compute_mean_and_variance

    type(global_transport_experiment_t), intent(in) :: batch_experiments(:)
    type(global_transport_experiment_t), intent(out) :: experiment

    real(dp), allocatable :: density_samples(:, :)
    real(dp), allocatable :: flux_samples(:, :)
    real(dp), allocatable :: integrated_flux_samples(:, :)
    integer :: i_batch

    experiment%boundary_s = batch_experiments(1)%boundary_s
    experiment%shell_volumes = batch_experiments(1)%shell_volumes

    allocate(density_samples(size(batch_experiments(1)%density), size(batch_experiments)))
    allocate(flux_samples(size(batch_experiments(1)%flux), size(batch_experiments)))
    allocate(integrated_flux_samples(size(batch_experiments(1)%integrated_flux), size(batch_experiments)))
    do i_batch = 1, size(batch_experiments)
        if (i_batch == 1) then
            allocate(experiment%source(size(batch_experiments(i_batch)%source)))
            experiment%source = 0.0_dp
        end if
        density_samples(:, i_batch) = batch_experiments(i_batch)%density
        flux_samples(:, i_batch) = batch_experiments(i_batch)%flux
        integrated_flux_samples(:, i_batch) = batch_experiments(i_batch)%integrated_flux
        experiment%source = experiment%source + batch_experiments(i_batch)%source
    end do

    experiment%source = experiment%source / real(size(batch_experiments), dp)
    call compute_mean_and_variance(density_samples, experiment%density, experiment%density_variance)
    call compute_mean_and_variance(flux_samples, experiment%flux, experiment%flux_variance)
    call compute_mean_and_variance(integrated_flux_samples, experiment%integrated_flux, experiment%integrated_flux_variance)
    experiment%flux(1) = 0.0_dp
    experiment%flux_variance(1) = 0.0_dp
    experiment%integrated_flux(1) = 0.0_dp
    experiment%integrated_flux_variance(1) = 0.0_dp

end subroutine aggregate_batched_experiments

subroutine evaluate_transport_signal(experiment, batch_lost_fraction, batch_tracing_time, sigma_multiplier, &
    min_supported_flux_boundaries, min_density_source_relstd, n_trials, quality)

    use global_transport_fit_types_mod, only: global_transport_experiment_t, global_transport_quality_t
    use transport_statistics_mod, only: compute_density_source_relstd, count_supported_boundaries, transport_signal_supported

    type(global_transport_experiment_t), intent(in) :: experiment
    real(dp), intent(in) :: batch_lost_fraction(:)
    real(dp), intent(in) :: batch_tracing_time(:)
    real(dp), intent(in) :: min_density_source_relstd
    real(dp), intent(in) :: sigma_multiplier
    integer, intent(in) :: min_supported_flux_boundaries
    integer, intent(in) :: n_trials
    type(global_transport_quality_t), intent(out) :: quality

    call compute_density_source_relstd(experiment%source, experiment%density, quality%density_source_relstd)
    quality%n_supported_flux_boundaries = count_supported_boundaries(experiment%integrated_flux, &
        experiment%integrated_flux_variance, sigma_multiplier)
    quality%lost_fraction = sum(batch_lost_fraction) / real(size(batch_lost_fraction), dp)
    quality%tracing_time = sum(batch_tracing_time) / real(size(batch_tracing_time), dp)
    quality%n_batches = size(batch_lost_fraction)
    quality%n_trials = n_trials
    quality%transport_supported = transport_signal_supported(quality%n_supported_flux_boundaries, min_supported_flux_boundaries, &
        quality%density_source_relstd, min_density_source_relstd)

end subroutine evaluate_transport_signal

subroutine write_experiment_summary(boundary_s, shell_volumes, boundary_areas, experiment, quality, filename, source_gradient, &
    source_reference_s, sigma_multiplier)

    use global_transport_fit_types_mod, only: global_transport_experiment_t, global_transport_quality_t
    use transport_statistics_mod, only: build_supported_boundary_mask

    real(dp), intent(in) :: boundary_areas(:)
    real(dp), intent(in) :: boundary_s(:)
    real(dp), intent(in) :: shell_volumes(:)
    character(len=*), intent(in) :: filename
    type(global_transport_experiment_t), intent(in) :: experiment
    type(global_transport_quality_t), intent(in) :: quality
    real(dp), intent(in) :: source_gradient
    real(dp), intent(in) :: source_reference_s
    real(dp), intent(in) :: sigma_multiplier

    integer :: i
    integer :: io_unit
    logical, allocatable :: supported_mask(:)

    call build_supported_boundary_mask(experiment%integrated_flux, experiment%integrated_flux_variance, sigma_multiplier, &
        supported_mask)

    open(newunit=io_unit, file=filename, status='replace', action='write')
    write(io_unit, '(A)') 'shell_index shell_center_s source_density density density_2sigma right_boundary_flux_density ' // &
        'right_boundary_flux_2sigma right_boundary_integrated_flux right_boundary_integrated_flux_2sigma right_boundary_supported'
    do i = 1, size(shell_volumes)
        write(io_unit, '(I0,8(1X,ES24.16),1X,L1)') i, 0.5_dp * (boundary_s(i) + boundary_s(i + 1)), experiment%source(i), &
            experiment%density(i), 2.0_dp * sqrt(experiment%density_variance(i)), experiment%flux(i + 1), &
            2.0_dp * sqrt(experiment%flux_variance(i + 1)), experiment%integrated_flux(i + 1), &
            2.0_dp * sqrt(experiment%integrated_flux_variance(i + 1)), supported_mask(i + 1)
    end do
    write(io_unit, '(A)') &
        'source_gradient source_reference_s tracing_time lost_fraction n_batches n_trials density_source_relstd ' // &
        'supported_flux_boundaries transport_supported'
    write(io_unit, '(4(1X,ES24.16),2(1X,I0),1X,ES24.16,1X,I0,1X,L1)') source_gradient, source_reference_s, &
        quality%tracing_time, quality%lost_fraction, quality%n_batches, quality%n_trials, quality%density_source_relstd, &
        quality%n_supported_flux_boundaries, quality%transport_supported
    close(io_unit)

end subroutine write_experiment_summary

subroutine build_local_transport_reference(boundary_s, n_reference_surfaces, density_gradient_1, density_gradient_2, &
    n_local_gradients, v_E, total_time, n_particles, n_batches, boundary_indices, &
    local_boundary_s, local_a, local_a_std, local_b, local_b_std, local_valid)

    use gorilla_applets_types_mod, only: in
    use transport_statistics_mod, only: build_gradient_grid, compute_mean_and_variance, count_supported_values, &
        fit_linear_transport_response
    use usual_transport_benchmark_mod, only: run_usual_transport_benchmark_case, usual_transport_result_t

    real(dp), intent(in) :: boundary_s(:)
    real(dp), intent(in) :: density_gradient_1
    real(dp), intent(in) :: density_gradient_2
    real(dp), intent(in) :: total_time
    real(dp), intent(in) :: v_E
    integer, intent(in) :: n_particles
    integer, intent(in) :: n_batches
    integer, intent(in) :: n_local_gradients
    integer, intent(in) :: n_reference_surfaces
    integer, allocatable, intent(out) :: boundary_indices(:)
    real(dp), allocatable, intent(out) :: local_a(:)
    real(dp), allocatable, intent(out) :: local_a_std(:)
    real(dp), allocatable, intent(out) :: local_b(:)
    real(dp), allocatable, intent(out) :: local_b_std(:)
    real(dp), allocatable, intent(out) :: local_boundary_s(:)
    logical, allocatable, intent(out) :: local_valid(:)

    type(usual_transport_result_t) :: result
    real(dp), allocatable :: gradients(:)
    real(dp), allocatable :: normalized_flux_mean(:)
    real(dp), allocatable :: normalized_flux_samples(:, :)
    real(dp), allocatable :: normalized_flux_var(:)
    real(dp) :: local_reference_s
    integer :: i
    integer :: i_batch
    integer :: i_gradient
    integer :: n_supported_gradients

    call select_local_reference_boundaries(boundary_s, n_reference_surfaces, boundary_indices)
    call build_gradient_grid(density_gradient_1, density_gradient_2, n_local_gradients, gradients)
    allocate(local_boundary_s(size(boundary_indices)))
    allocate(local_a(size(boundary_indices)))
    allocate(local_a_std(size(boundary_indices)))
    allocate(local_b(size(boundary_indices)))
    allocate(local_b_std(size(boundary_indices)))
    allocate(local_valid(size(boundary_indices)))
    allocate(normalized_flux_samples(n_local_gradients, n_batches))

    do i = 1, size(boundary_indices)
        local_reference_s = boundary_s(boundary_indices(i))
        do i_gradient = 1, n_local_gradients
            do i_batch = 1, n_batches
                call run_usual_transport_benchmark_case(local_reference_s, 1, &
                    in%density, in%density, in%energy_eV, &
                    in%density, in%energy_eV, in%energy_eV, &
                    total_time, n_particles, gradients(i_gradient), &
                    local_reference_s, v_E, in%collision_operator, &
                    in%boole_precalc_collisions, in%i_integrator_type, &
                    result, reset_random_seed=.false., &
                    boole_boltzmann_energies_in=.true., &
                    boole_static_ne_in=in%boole_static_ne)
                normalized_flux_samples(i_gradient, i_batch) = &
                    result%weighted_flux_density / in%density
            end do
        end do
        call compute_mean_and_variance(normalized_flux_samples, normalized_flux_mean, normalized_flux_var)
        call fit_linear_transport_response(gradients, normalized_flux_mean, normalized_flux_var, local_a(i), local_b(i), &
            local_a_std(i), local_b_std(i), local_valid(i))
        n_supported_gradients = count_supported_values(normalized_flux_mean, normalized_flux_var, 2.0_dp)
        if (n_supported_gradients < 2) then
            local_valid(i) = .false.
            local_a(i) = 0.0_dp
            local_b(i) = 0.0_dp
            local_a_std(i) = huge(1.0_dp)
            local_b_std(i) = huge(1.0_dp)
        end if
        local_boundary_s(i) = boundary_s(boundary_indices(i))
    end do

end subroutine build_local_transport_reference

subroutine select_local_reference_boundaries(boundary_s, n_reference_surfaces, boundary_indices)

    real(dp), intent(in) :: boundary_s(:)
    integer, intent(in) :: n_reference_surfaces
    integer, allocatable, intent(out) :: boundary_indices(:)

    integer :: end_index
    integer :: i
    integer :: margin
    integer :: n_interior
    integer :: n_selected
    integer :: start_index

    n_interior = max(size(boundary_s) - 2, 1)
    n_selected = min(max(n_reference_surfaces, 1), n_interior)
    allocate(boundary_indices(n_selected))

    margin = max(1, n_interior / 10)
    start_index = min(1 + margin, size(boundary_s) - 1)
    end_index = max(start_index, size(boundary_s) - 1 - margin)

    if (n_selected == 1) then
        boundary_indices(1) = (start_index + end_index) / 2
        return
    end if

    do i = 1, n_selected
        boundary_indices(i) = start_index + int(real(i - 1, dp) * real(end_index - start_index, dp) / &
            real(n_selected - 1, dp))
    end do

end subroutine select_local_reference_boundaries

subroutine write_local_transport_profiles(boundary_indices, boundary_s, local_valid, local_a, local_a_std, local_b, local_b_std, &
    filename)

    integer, intent(in) :: boundary_indices(:)
    real(dp), intent(in) :: boundary_s(:)
    logical, intent(in) :: local_valid(:)
    real(dp), intent(in) :: local_a(:)
    real(dp), intent(in) :: local_a_std(:)
    real(dp), intent(in) :: local_b(:)
    real(dp), intent(in) :: local_b_std(:)
    character(len=*), intent(in) :: filename

    integer :: i
    integer :: io_unit
    real(dp) :: local_a_2sigma
    real(dp) :: local_b_2sigma
    real(dp) :: sanitized_a
    real(dp) :: sanitized_a_std
    real(dp) :: sanitized_b
    real(dp) :: sanitized_b_std

    open(newunit=io_unit, file=filename, status='replace', action='write')
    write(io_unit, '(A)') 'boundary_index boundary_s local_valid local_A local_A_std local_A_2sigma local_B local_B_std local_B_2sigma'
    do i = 1, size(boundary_indices)
        call sanitize_transport_output_triplet(local_valid(i), local_a(i), local_a_std(i), sanitized_a, sanitized_a_std, local_a_2sigma)
        call sanitize_transport_output_triplet(local_valid(i), local_b(i), local_b_std(i), sanitized_b, sanitized_b_std, local_b_2sigma)
        write(io_unit, '(I0,1X,ES24.16,1X,L1,6(1X,ES24.16))') boundary_indices(i), boundary_s(i), local_valid(i), &
            sanitized_a, sanitized_a_std, local_a_2sigma, sanitized_b, sanitized_b_std, local_b_2sigma
    end do
    close(io_unit)

end subroutine write_local_transport_profiles

subroutine write_transport_comparison(boundary_indices, boundary_s, local_valid, local_a, local_a_std, local_b, local_b_std, &
    fit_a, fit_a_std, fit_b, fit_b_std, filename)

    integer, intent(in) :: boundary_indices(:)
    real(dp), intent(in) :: boundary_s(:)
    logical, intent(in) :: local_valid(:)
    real(dp), intent(in) :: local_a(:)
    real(dp), intent(in) :: local_a_std(:)
    real(dp), intent(in) :: local_b(:)
    real(dp), intent(in) :: local_b_std(:)
    real(dp), intent(in) :: fit_a(:)
    real(dp), intent(in) :: fit_a_std(:)
    real(dp), intent(in) :: fit_b(:)
    real(dp), intent(in) :: fit_b_std(:)
    character(len=*), intent(in) :: filename

    integer :: i
    integer :: io_unit
    real(dp) :: fit_a_2sigma
    real(dp) :: fit_b_2sigma
    real(dp) :: local_a_2sigma
    real(dp) :: local_b_2sigma
    real(dp) :: sanitized_fit_a
    real(dp) :: sanitized_fit_a_std
    real(dp) :: sanitized_fit_b
    real(dp) :: sanitized_fit_b_std
    real(dp) :: sanitized_local_a
    real(dp) :: sanitized_local_a_std
    real(dp) :: sanitized_local_b
    real(dp) :: sanitized_local_b_std

    open(newunit=io_unit, file=filename, status='replace', action='write')
    write(io_unit, '(A)') &
        'boundary_index boundary_s local_valid local_A local_A_std local_A_2sigma fit_A fit_A_std fit_A_2sigma ' // &
        'local_B local_B_std local_B_2sigma fit_B fit_B_std fit_B_2sigma'
    do i = 1, size(boundary_indices)
        call sanitize_transport_output_triplet(local_valid(i), local_a(i), local_a_std(i), sanitized_local_a, sanitized_local_a_std, &
            local_a_2sigma)
        call sanitize_transport_output_triplet(local_valid(i), local_b(i), local_b_std(i), sanitized_local_b, sanitized_local_b_std, &
            local_b_2sigma)
        call sanitize_transport_output_triplet(.true., fit_a(boundary_indices(i)), fit_a_std(boundary_indices(i)), sanitized_fit_a, &
            sanitized_fit_a_std, fit_a_2sigma)
        call sanitize_transport_output_triplet(.true., fit_b(boundary_indices(i)), fit_b_std(boundary_indices(i)), sanitized_fit_b, &
            sanitized_fit_b_std, fit_b_2sigma)
        write(io_unit, '(I0,1X,ES24.16,1X,L1,12(1X,ES24.16))') boundary_indices(i), boundary_s(i), local_valid(i), &
            sanitized_local_a, sanitized_local_a_std, local_a_2sigma, sanitized_fit_a, sanitized_fit_a_std, fit_a_2sigma, &
            sanitized_local_b, sanitized_local_b_std, local_b_2sigma, sanitized_fit_b, sanitized_fit_b_std, fit_b_2sigma
    end do
    close(io_unit)

end subroutine write_transport_comparison

end module global_transport_fit_gorilla_mod
